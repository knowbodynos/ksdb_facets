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
import sys

import pyarrow.fs as pafs
import pyarrow.parquet as pq
from tqdm import tqdm

DEFAULT_REGION = "us-east-2"
FACET_NF_COL   = "facet_nfs"
CONE_NF_COL    = "maximal_cone_nfs"


def to_hashable(v):
    if isinstance(v, list):
        return tuple(to_hashable(x) for x in v)
    return v


def get_fs_and_path(path, region):
    if path.startswith("s3://"):
        return pafs.S3FileSystem(region=region), path.removeprefix("s3://").rstrip("/")
    return pafs.LocalFileSystem(), path.rstrip("/")


def list_parquets(fs, path):
    infos = fs.get_file_info(pafs.FileSelector(path, recursive=True))
    return [
        info.path for info in infos
        if info.type == pafs.FileType.File and info.path.endswith(".parquet")
    ]


def main():
    ap = argparse.ArgumentParser(
        description="Count unique GL(3,Z) and GL(4,Z) normal forms across parquet files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("path", help="Local directory or s3://bucket/prefix/ of output parquet files")
    ap.add_argument("--region", default=DEFAULT_REGION, help="AWS region (S3 only)")
    args = ap.parse_args()

    fs, path = get_fs_and_path(args.path, args.region)

    print("Listing files…", flush=True)
    files = list_parquets(fs, path)
    if not files:
        sys.exit(f"No parquet files found under {args.path}")
    print(f"  {len(files)} file(s)", flush=True)

    print("Reading file metadata…", flush=True)
    file_meta = []
    for fpath in files:
        with fs.open_input_file(fpath) as f:
            file_meta.append((fpath, pq.ParquetFile(f).metadata))
    total_rgs = sum(m.num_row_groups for _, m in file_meta)

    facet_nf_set  = set()
    cone_nf_set   = set()
    total_facets  = 0
    total_polys   = 0

    with tqdm(total=total_rgs, unit="rg") as bar:
        for fpath, _ in file_meta:
            bar.set_postfix_str(fpath.rsplit("/", 1)[-1], refresh=False)
            with fs.open_input_file(fpath) as f:
                pf = pq.ParquetFile(f)
                schema_names = pf.schema_arrow.names
                cols = [c for c in [FACET_NF_COL, CONE_NF_COL] if c in schema_names]
                if not cols:
                    bar.update(pf.metadata.num_row_groups)
                    continue
                for i in range(pf.metadata.num_row_groups):
                    batch = pf.read_row_group(i, columns=cols)
                    total_polys += batch.num_rows

                    if FACET_NF_COL in cols:
                        for row_nfs in batch.column(FACET_NF_COL).to_pylist():
                            total_facets += len(row_nfs)
                            for nf in row_nfs:
                                facet_nf_set.add(to_hashable(nf))

                    if CONE_NF_COL in cols:
                        for row_nfs in batch.column(CONE_NF_COL).to_pylist():
                            if FACET_NF_COL not in cols:
                                total_facets += len(row_nfs)
                            for nf in row_nfs:
                                cone_nf_set.add(to_hashable(nf))

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
