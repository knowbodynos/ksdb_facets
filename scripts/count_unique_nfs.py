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
import os
import struct
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import pyarrow.fs as pafs
import pyarrow.parquet as pq
from tqdm import tqdm

DEFAULT_REGION  = "us-east-2"
DEFAULT_WORKERS = os.cpu_count() or 1
FACET_NF_COL    = "facet_nfs"
CONE_NF_COL     = "maximal_cone_nfs"


def nf_to_key(nf):
    """Pack a normal form (list of lists of int) into a bytes key for hashing."""
    rows = len(nf)
    cols = len(nf[0]) if rows else 0
    flat = [v for row in nf for v in row]
    return struct.pack(f'HH{len(flat)}i', rows, cols, *flat)


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


def _process_file(args):
    """Worker: read one file, return (facet_key_set, cone_key_set, n_polys, n_facets)."""
    fpath, is_s3, region = args
    fs = pafs.S3FileSystem(region=region) if is_s3 else pafs.LocalFileSystem()

    facet_keys = set()
    cone_keys  = set()
    n_polys    = 0
    n_facets   = 0

    with fs.open_input_file(fpath) as f:
        pf = pq.ParquetFile(f)
        cols = [c for c in [FACET_NF_COL, CONE_NF_COL] if c in pf.schema_arrow.names]
        if not cols:
            return facet_keys, cone_keys, 0, 0

        for i in range(pf.metadata.num_row_groups):
            batch = pf.read_row_group(i, columns=cols)
            n_polys += batch.num_rows

            if FACET_NF_COL in cols:
                for row_nfs in batch.column(FACET_NF_COL).to_pylist():
                    n_facets += len(row_nfs)
                    for nf in row_nfs:
                        facet_keys.add(nf_to_key(nf))

            if CONE_NF_COL in cols:
                for row_nfs in batch.column(CONE_NF_COL).to_pylist():
                    if FACET_NF_COL not in cols:
                        n_facets += len(row_nfs)
                    for nf in row_nfs:
                        cone_keys.add(nf_to_key(nf))

    return facet_keys, cone_keys, n_polys, n_facets


def main():
    ap = argparse.ArgumentParser(
        description="Count unique GL(3,Z) and GL(4,Z) normal forms across parquet files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("path", help="Local directory or s3://bucket/prefix/ of output parquet files")
    ap.add_argument("--workers", type=int, default=DEFAULT_WORKERS,
                    help="Number of parallel worker processes")
    ap.add_argument("--region", default=DEFAULT_REGION, help="AWS region (S3 only)")
    args = ap.parse_args()

    fs, path, is_s3 = get_fs_and_path(args.path, args.region)

    print("Listing files…", flush=True)
    files = list_parquets(fs, path)
    if not files:
        sys.exit(f"No parquet files found under {args.path}")
    print(f"  {len(files)} file(s)  |  {args.workers} worker(s)", flush=True)

    facet_nf_set = set()
    cone_nf_set  = set()
    total_facets = 0
    total_polys  = 0

    worker_args = [(fpath, is_s3, args.region) for fpath in files]

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_process_file, arg): arg[0] for arg in worker_args}
        with tqdm(total=len(files), unit="file") as bar:
            for future in as_completed(futures):
                fpath = futures[future]
                bar.set_postfix_str(fpath.rsplit("/", 1)[-1], refresh=False)
                facet_keys, cone_keys, n_polys, n_facets = future.result()
                facet_nf_set.update(facet_keys)
                cone_nf_set.update(cone_keys)
                total_polys  += n_polys
                total_facets += n_facets
                bar.update(1)

    n_facet_nfs = len(facet_nf_set)
    n_cone_nfs  = len(cone_nf_set)

    print(f"\n  polytopes processed : {total_polys:,}")
    print(f"  total facets seen   : {total_facets:,}")
    print()
    print(f"  unique facet NFs (GL(3,Z)) : {n_facet_nfs:,}")
    print(f"  unique cone NFs  (GL(4,Z)) : {n_cone_nfs:,}")

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
