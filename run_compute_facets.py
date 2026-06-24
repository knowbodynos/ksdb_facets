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
        facet_results_test.parquet
"""

import argparse
import datetime
import json
import os
import queue
import signal
import subprocess
import sys
import threading
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
DEFAULT_WORKERS    = os.cpu_count() or 4

# ---------------------------------------------------------------------------
# Persistent worker pool
# ---------------------------------------------------------------------------

class WorkerPool:
    """
    N persistent compute_facets subprocesses, each owned by one thread.
    The binary speaks line-by-line JSON (one row in → one row out), so each
    thread sends rows one at a time and reads responses without restarting the
    process between chunks.
    """
    def __init__(self, bin_path: str, n_workers: int, log_path: str = None):
        self._task_q   = queue.Queue(maxsize=n_workers * 4)
        self._result_q = queue.Queue()
        self._log_lock = threading.Lock()
        self._log_file = open(log_path, 'a') if log_path else None
        self._threads  = []
        for _ in range(n_workers):
            t = threading.Thread(target=self._run, args=(bin_path,), daemon=True)
            t.start()
            self._threads.append(t)

    def _warn(self, msg: str):
        print(f"\n  [warn] {msg}", file=sys.stderr)
        if self._log_file is not None:
            with self._log_lock:
                self._log_file.write(f"{datetime.datetime.now().isoformat()} [warn] {msg}\n")
                self._log_file.flush()

    @staticmethod
    def _exit_str(returncode: int) -> str:
        if returncode >= 0:
            return f"exit {returncode}"
        try:
            name = signal.Signals(-returncode).name
        except ValueError:
            name = f"signal {-returncode}"
        s = f"killed by {name}"
        if name == "SIGKILL":
            s += " (possible OOM — check: dmesg | grep -i 'killed process')"
        return s

    @staticmethod
    def _start_proc(bin_path: str):
        return subprocess.Popen(
            [bin_path],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

    def _run(self, bin_path: str):
        proc = self._start_proc(bin_path)
        while True:
            item = self._task_q.get()
            if item is None:
                self._task_q.task_done()
                break
            idx, verts = item
            try:
                proc.stdin.write(
                    json.dumps({"index": idx, "verts": verts}, separators=(',', ':')) + '\n'
                )
                proc.stdin.flush()
                line = proc.stdout.readline().strip()
                if line:
                    self._result_q.put(json.loads(line))
                elif proc.poll() is not None:
                    stderr_out = proc.stderr.read().strip()
                    self._warn(f"worker crashed on index {idx} ({self._exit_str(proc.returncode)})"
                               + (f": {stderr_out}" if stderr_out else ""))
                    proc = self._start_proc(bin_path)
            except BrokenPipeError:
                proc.wait()
                stderr_out = proc.stderr.read().strip()
                self._warn(f"broken pipe on index {idx} ({self._exit_str(proc.returncode)})"
                           + (f": {stderr_out}" if stderr_out else ""))
                proc = self._start_proc(bin_path)
            except Exception as exc:
                self._warn(f"worker error on index {idx}: {exc}")
            self._task_q.task_done()
        proc.stdin.close()
        proc.wait()

    def submit(self, idx: int, verts) -> None:
        self._task_q.put((idx, verts))

    def join(self) -> list:
        """Block until all submitted tasks finish; return collected results."""
        self._task_q.join()
        results = []
        while not self._result_q.empty():
            results.append(self._result_q.get_nowait())
        return results

    def close(self):
        for _ in self._threads:
            self._task_q.put(None)
        for t in self._threads:
            t.join()
        if self._log_file is not None:
            self._log_file.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


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

def process_batch(pool, global_indices, verts_list):
    """
    global_indices : list of int  (source parquet row indices)
    verts_list     : list of list-of-list-of-int
    Returns a pyarrow RecordBatch (or None if empty).
    """
    for idx, verts in zip(global_indices, verts_list):
        pool.submit(idx, verts)

    all_results = pool.join()
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
    ap.add_argument("--limit",      type=int, default=None,
                    help="Only process the first N rows (for testing)")
    ap.add_argument("--s3-uri",      type=str, default=None,
                    help="Upload output to this S3 URI then shut down the instance")
    ap.add_argument("--log",         type=str, default=None,
                    help="Append worker warnings to this file")
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
          f"batch={args.batch_size}", flush=True)

    # --- Streaming write ---
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    writer = pq.ParquetWriter(str(out_path), OUT_SCHEMA, compression="snappy")

    n_written = 0
    with WorkerPool(str(bin_path), args.workers, log_path=args.log) as pool:
        with tqdm(total=total, unit="poly", smoothing=0.05) as bar:
            offset = 0
            while offset < total:
                batch_n = min(args.batch_size, total - offset)
                df = lf.slice(offset, batch_n).collect()
                global_indices = list(range(offset, offset + len(df)))
                verts_list = df[args.verts_col].to_list()

                rb = process_batch(pool, global_indices, verts_list)
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
