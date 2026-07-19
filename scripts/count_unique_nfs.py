#!/usr/bin/env python3
"""
count_unique_nfs.py
-------------------
Count unique facet normal forms (GL(3,Z)) and maximal cone normal forms
(GL(4,Z)) across all parquet files in a local directory or S3 prefix.

Usage:
    uv run python count_unique_nfs.py PATH

Examples:
    uv run python count_unique_nfs.py facet_results/
    uv run python count_unique_nfs.py s3://toriccy/shared/4d-polytope-facets/
"""

import argparse
import datetime
import hashlib
import os
import resource
import struct
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pyarrow.fs as pafs
import pyarrow.parquet as pq
from tqdm import tqdm

DEFAULT_REGION  = "us-east-2"
DEFAULT_WORKERS = os.cpu_count() or 1
FACET_NF_COL    = "facet_nfs"
CONE_NF_COL     = "maximal_cone_nfs"


def nf_to_key(nf):
    rows = len(nf)
    cols = len(nf[0]) if rows else 0
    flat = [v for row in nf for v in row]
    data = struct.pack(f'HH{len(flat)}i', rows, cols, *flat)
    return int.from_bytes(hashlib.blake2b(data, digest_size=8).digest(), 'little')


def get_fs_and_path(path, region):
    if path.startswith("s3://"):
        return pafs.S3FileSystem(region=region), path.removeprefix("s3://").rstrip("/"), True
    return pafs.LocalFileSystem(), path.rstrip("/"), False


def list_parquets(fs, path):
    infos = fs.get_file_info(pafs.FileSelector(path, recursive=True))
    return [
        info.path for info in infos
        if info.type == pafs.FileType.File and info.path.endswith(".parquet")
    ]


def _read_rg_count(args):
    fpath, fs = args
    with fs.open_input_file(fpath) as f:
        return fpath, pq.ParquetFile(f).metadata.num_row_groups


def _process_file(fpath, fs, progress_callback=None):
    """Read one parquet file; return deduplicated uint64 hash arrays + counts."""
    facet_hashes = []
    cone_hashes  = []
    n_polys      = 0
    n_facets     = 0

    with fs.open_input_file(fpath) as f:
        pf = pq.ParquetFile(f)
        cols = [c for c in [FACET_NF_COL, CONE_NF_COL] if c in pf.schema_arrow.names]
        if not cols:
            if progress_callback:
                for _ in range(pf.metadata.num_row_groups):
                    progress_callback()
            return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.uint64), 0, 0

        for i in range(pf.metadata.num_row_groups):
            batch = pf.read_row_group(i, columns=cols)
            n_polys += batch.num_rows

            if FACET_NF_COL in cols:
                for row_nfs in batch.column(FACET_NF_COL).to_pylist():
                    n_facets += len(row_nfs)
                    for nf in row_nfs:
                        facet_hashes.append(nf_to_key(nf))

            if CONE_NF_COL in cols:
                for row_nfs in batch.column(CONE_NF_COL).to_pylist():
                    if FACET_NF_COL not in cols:
                        n_facets += len(row_nfs)
                    for nf in row_nfs:
                        cone_hashes.append(nf_to_key(nf))

            if progress_callback:
                progress_callback()

    return (
        np.unique(np.array(facet_hashes, dtype=np.uint64)),
        np.unique(np.array(cone_hashes,  dtype=np.uint64)),
        n_polys,
        n_facets,
    )


def rss_gb():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != "darwin":  # Linux reports KB; macOS reports bytes
        rss *= 1024
    return rss / 1e9


def main():
    ap = argparse.ArgumentParser(
        description="Count unique GL(3,Z) and GL(4,Z) normal forms across parquet files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("path", help="Local directory or s3://bucket/prefix/ of output parquet files")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                    help="Number of parallel threads")
    ap.add_argument("--region", default=DEFAULT_REGION, help="AWS region (S3 only)")
    args = ap.parse_args()

    fs, path, _ = get_fs_and_path(args.path, args.region)

    print("Listing files…", flush=True)
    files = list_parquets(fs, path)
    if not files:
        sys.exit(f"No parquet files found under {args.path}")
    print(f"  {len(files)} file(s) | {args.workers} worker(s)")

    print("Reading file metadata…", flush=True)
    with ThreadPoolExecutor(max_workers=min(32, len(files))) as tex:
        rg_counts = dict(tqdm(
            tex.map(_read_rg_count, [(fpath, fs) for fpath in files]),
            total=len(files), unit="file", desc="  metadata",
        ))
    total_rgs = sum(rg_counts.values())
    print(f"  {total_rgs:,} row groups total | RSS {rss_gb():.2f} GB", flush=True)

    facet_nf_arr = np.empty(0, dtype=np.uint64)
    cone_nf_arr  = np.empty(0, dtype=np.uint64)
    total_facets = 0
    total_polys  = 0

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_process_file, fpath, fs): fpath for fpath in files}
        with tqdm(total=total_rgs, unit="rg") as bar:
            for future in as_completed(futures):
                fpath = futures[future]
                bar.set_postfix_str(fpath.rsplit("/", 1)[-1], refresh=False)
                try:
                    facet_arr, cone_arr, n_polys, n_facets = future.result()
                    facet_nf_arr = np.unique(np.concatenate([facet_nf_arr, facet_arr]))
                    cone_nf_arr  = np.unique(np.concatenate([cone_nf_arr,  cone_arr]))
                    del facet_arr, cone_arr
                    total_polys  += n_polys
                    total_facets += n_facets
                    bar.update(rg_counts.get(fpath, 1))
                    bar.set_postfix(rss=f"{rss_gb():.2f}GB", refresh=False)
                except Exception as exc:
                    tqdm.write(f"[warn] {fpath.rsplit('/', 1)[-1]}: {exc}")
                    bar.update(rg_counts.get(fpath, 1))

    n_facet_nfs = len(facet_nf_arr)
    n_cone_nfs  = len(cone_nf_arr)

    print(f"\n  polytopes processed        : {total_polys:,}")
    print(f"  total facets seen          : {total_facets:,}")
    print()
    print(f"  unique facet NFs (GL(3,Z)) : {n_facet_nfs:,}")
    print(f"  unique cone NFs  (GL(4,Z)) : {n_cone_nfs:,}")
    print(f"\n  peak RSS: {rss_gb():.2f} GB")

    info_path = path + "/INFO.md"
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    info_md = f"""\
# 4D Polytope Facets — Normal Form Statistics

Generated: {timestamp}

| Metric | Count |
|--------|------:|
| Files | {len(files):,} |
| Polytopes | {total_polys:,} |
| Total facets | {total_facets:,} |
| Unique facet NFs (GL(3,Z)) | {n_facet_nfs:,} |
| Unique maximal cone NFs (GL(4,Z)) | {n_cone_nfs:,} |
"""
    with fs.open_output_stream(info_path) as out:
        out.write(info_md.encode())

    prefix = "s3://" if args.path.startswith("s3://") else ""
    print(f"\nWrote {prefix}{info_path}")


if __name__ == "__main__":
    main()
