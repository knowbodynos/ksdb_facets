#!/usr/bin/env python3
"""
run_compute_facets.py
--------------
Run compute_facets on every row of a parquet file and write results to a new parquet.

For each polytope the output contains:
  index            – original row index in the source parquet (pass-through)
  verts            – List(List(Int32)): original 4D polytope vertices
  dual_verts       – List(List(Int32)): vertices of the dual polytope (facet normals);
                     dual_verts[i] is the inner normal of facets[i] (inner product == -1)
  facets           – List(List(List(Int32))): 3D facets in dual-vertex order;
                     facets[i] corresponds 1:1 to dual_verts[i]
  facet_nfs        – List(List(List(Int32))): GL(3,Z) normal form of each 3D facet;
                     canonical invariant of the facet as an abstract 3D lattice polytope
  maximal_cone_nfs – List(List(List(Int32))): GL(4,Z) normal form of conv(facets[i] ∪ {0});
                     finer invariant capturing the embedding of the facet in Z^4

Usage:
    uv run python run_compute_facets.py INPUT OUTPUT [options]

    INPUT   parquet path, glob, or hf:// URL
    OUTPUT  output parquet path

Examples:
    # Full dataset
    uv run python run_compute_facets.py \
        "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
        facet_results.parquet

    # Quick test on first 1000 rows
    uv run python run_compute_facets.py --limit 1000 \
        "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
        facet_results.parquet
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
DEFAULT_BIN        = Path(__file__).resolve().parent / "compute_facets" / "compute_facets"
DEFAULT_BATCH_SIZE = 10_000   # rows read into memory per iteration
DEFAULT_CHUNK_SIZE = 500      # rows sent to one compute_facets subprocess
DEFAULT_WORKERS    = os.cpu_count() or 4

# ---------------------------------------------------------------------------
# Worker (serialisable: runs inside ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _worker(args):
    """
    Pipe one sub-chunk through compute_facets.
    args = ([(global_idx, verts), ...], bin_path_str)
    Returns list of parsed JSON result dicts.
    """
    rows, bin_path = args
    lines = [
        json.dumps({"index": idx, "verts": verts}, separators=(',', ':'))
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
    pa.field("index",            pa.int64()),
    pa.field("verts",            pa.list_(pa.list_(pa.int32()))),
    pa.field("dual_verts",       pa.list_(pa.list_(pa.int32()))),
    pa.field("facets",           pa.list_(pa.list_(pa.list_(pa.int32())))),
    pa.field("facet_nfs",        pa.list_(pa.list_(pa.list_(pa.int32())))),
    pa.field("maximal_cone_nfs", pa.list_(pa.list_(pa.list_(pa.int32())))),
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

    all_results.sort(key=lambda r: r["index"])

    out_index        = []
    out_verts        = []
    out_dual_verts   = []
    out_facets       = []
    out_facet_nfs    = []
    out_maxcone_nfs  = []

    for r in all_results:
        out_index.append(r["index"])
        out_verts.append(r["verts"])
        out_dual_verts.append(r["dual_verts"])
        out_facets.append(r["facets"])
        out_facet_nfs.append(r["facet_nfs"])
        out_maxcone_nfs.append(r["maximal_cone_nfs"])

    return pa.record_batch(
        [
            pa.array(out_index,       type=pa.int64()),
            pa.array(out_verts,       type=pa.list_(pa.list_(pa.int32()))),
            pa.array(out_dual_verts,  type=pa.list_(pa.list_(pa.int32()))),
            pa.array(out_facets,      type=pa.list_(pa.list_(pa.list_(pa.int32())))),
            pa.array(out_facet_nfs,   type=pa.list_(pa.list_(pa.list_(pa.int32())))),
            pa.array(out_maxcone_nfs, type=pa.list_(pa.list_(pa.list_(pa.int32())))),
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
        description="Run compute_facets on every row of a parquet file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("input",  help="Source parquet path, glob, or hf:// URL")
    ap.add_argument("output", help="Output parquet path")
    ap.add_argument("--verts-col",  default="vertices",
                    help="Name of the vertices column in the source parquet")
    ap.add_argument("--bin",        default=str(DEFAULT_BIN),
                    help="Path to the compute_facets binary")
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
        sys.exit(f"ERROR: compute_facets binary not found: {bin_path}\n"
                 f"  Build it with:  cd compute_facets && make")

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
