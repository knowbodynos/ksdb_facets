#!/usr/bin/env python3
"""
run_cone_nf.py
--------------
Run cone_nf on every row of a parquet file and write results to a new parquet.

For each polytope the output contains:
  idx            – original row index in the source parquet
  vertices       – original List(List(Int32)) vertices (pass-through)
  cone_nfs       – List(List(List(Int32))): one NF vertex matrix per unique cone type
  multiplicities – List(Int32): how many facets share each cone NF (parallel to cone_nfs)

Usage:
    uv run python run_cone_nf.py INPUT OUTPUT [options]

    INPUT   parquet path, glob, or hf:// URL
    OUTPUT  output parquet path

Examples:
    # Full dataset
    uv run python run_cone_nf.py \
        "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
        cone_nf_results.parquet

    # Quick test on first 1000 rows
    uv run python run_cone_nf.py --limit 1000 \
        "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
        cone_nf_test.parquet
"""

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_BIN        = Path(__file__).resolve().parent / "cone_nf" / "cone_nf"
DEFAULT_BATCH_SIZE = 10_000   # rows read into memory per iteration
DEFAULT_CHUNK_SIZE = 500      # rows sent to one cone_nf subprocess
DEFAULT_WORKERS    = os.cpu_count() or 4

# ---------------------------------------------------------------------------
# Worker (serialisable: runs inside ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _worker(args):
    """
    Pipe one sub-chunk through cone_nf.
    args = ([(global_idx, verts), ...], bin_path_str)
    Returns list of parsed JSON result dicts.
    """
    rows, bin_path = args
    lines = [
        json.dumps({"POLYID": idx, "verts": verts}, separators=(',', ':'))
        for idx, verts in rows
    ]
    proc = subprocess.Popen(
        [bin_path],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, _ = proc.communicate('\n'.join(lines))
    results = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            results.append(json.loads(line))
        except json.JSONDecodeError:
            pass
    return results


# ---------------------------------------------------------------------------
# PyArrow output schema
# ---------------------------------------------------------------------------

OUT_SCHEMA = pa.schema([
    pa.field("idx",            pa.int64()),
    pa.field("vertices",       pa.list_(pa.list_(pa.int32()))),
    pa.field("cone_nfs",       pa.list_(pa.list_(pa.list_(pa.int32())))),
    pa.field("multiplicities", pa.list_(pa.int32())),
])


# ---------------------------------------------------------------------------
# Process one in-memory batch
# ---------------------------------------------------------------------------

def process_batch(global_indices, verts_list, bin_path, n_workers, chunk_size):
    """
    global_indices : list of int  (source parquet row indices)
    verts_list     : list of list-of-list-of-int
    Returns a pyarrow RecordBatch (or None if empty).
    """
    rows = list(zip(global_indices, verts_list))
    chunks = [
        (rows[i:i + chunk_size], str(bin_path))
        for i in range(0, len(rows), chunk_size)
    ]

    all_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futs = [pool.submit(_worker, c) for c in chunks]
        for fut in as_completed(futs):
            try:
                all_results.extend(fut.result())
            except Exception as exc:
                print(f"\n  [warn] worker error: {exc}", file=sys.stderr)

    if not all_results:
        return None

    all_results.sort(key=lambda r: r["POLYID"])

    # Build a lookup from global_idx → original verts
    verts_by_idx = dict(rows)

    out_idx   = []
    out_verts = []
    out_nfs   = []   # List[List[List[int]]]  — one NF matrix per unique cone
    out_mults = []   # List[int]              — multiplicity per unique cone

    for r in all_results:
        gidx = r["POLYID"]
        out_idx.append(gidx)
        out_verts.append(verts_by_idx[gidx])
        out_nfs.append([u["cone_nf"]      for u in r["unique_cones"]])
        out_mults.append([u["multiplicity"] for u in r["unique_cones"]])

    return pa.record_batch(
        [
            pa.array(out_idx,   type=pa.int64()),
            pa.array(out_verts, type=pa.list_(pa.list_(pa.int32()))),
            pa.array(out_nfs,   type=pa.list_(pa.list_(pa.list_(pa.int32())))),
            pa.array(out_mults, type=pa.list_(pa.int32())),
        ],
        schema=OUT_SCHEMA,
    )


def upload_to_s3(local_path: str, s3_uri: str) -> None:
    if not s3_uri.startswith("s3://"):
        raise ValueError("s3_uri must start with s3://")
    bucket, key = s3_uri[5:].split("/", 1)
    import boto3
    boto3.client("s3").upload_file(local_path, bucket, key)


def terminate_instance() -> None:
    # Requires instance role permission: ec2:TerminateInstances
    subprocess.run(["sudo", "shutdown", "-h", "now"], check=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Run cone_nf on every row of a parquet file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("input",  help="Source parquet path, glob, or hf:// URL")
    ap.add_argument("output", help="Output parquet path")
    ap.add_argument("--verts-col",  default="vertices",
                    help="Name of the vertices column in the source parquet")
    ap.add_argument("--bin",        default=str(DEFAULT_BIN),
                    help="Path to the cone_nf binary")
    ap.add_argument("--workers",    type=int, default=DEFAULT_WORKERS,
                    help="Number of parallel worker processes")
    ap.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE,
                    help="Rows loaded into memory per iteration")
    ap.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
                    help="Rows per worker subprocess call")
    ap.add_argument("--limit",      type=int, default=None,
                    help="Only process the first N rows (for testing)")
    ap.add_argument("--s3-uri",      type=str, default=None,
                    help="Upload output to this S3 URI then shut down the instance")
    args = ap.parse_args()

    bin_path = Path(args.bin)
    if not bin_path.exists():
        sys.exit(f"ERROR: cone_nf binary not found: {bin_path}\n"
                 f"  Build it with:  cd cone_nf && make")

    # --- Count rows ---
    print(f"Scanning {args.input} …", flush=True)
    lf = pl.scan_parquet(args.input)
    total = lf.select(pl.len()).collect().item()
    if args.limit:
        total = min(total, args.limit)
    print(f"  {total:,} rows to process  |  {args.workers} workers  |  "
          f"batch={args.batch_size}  chunk={args.chunk_size}", flush=True)

    # --- Streaming write ---
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = pq.ParquetWriter(str(out_path), OUT_SCHEMA, compression="snappy")

    n_written = 0
    with tqdm(total=total, unit="poly", smoothing=0.05) as bar:
        offset = 0
        while offset < total:
            batch_n = min(args.batch_size, total - offset)
            df = lf.slice(offset, batch_n).collect()
            global_indices = list(range(offset, offset + len(df)))
            verts_list = df[args.verts_col].to_list()

            rb = process_batch(
                global_indices, verts_list,
                bin_path, args.workers, args.chunk_size,
            )
            if rb is not None:
                writer.write_batch(rb)
                n_written += rb.num_rows

            bar.update(len(df))
            offset += len(df)

    writer.close()
    print(f"\nWrote {n_written:,} rows → {out_path}", flush=True)

    if args.s3_uri is not None:
        try:
            upload_to_s3(str(out_path), args.s3_uri)
            print(f"Uploaded {out_path} to {args.s3_uri}")
        finally:
            terminate_instance()


if __name__ == "__main__":
    main()
