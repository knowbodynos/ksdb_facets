#!/usr/bin/env python3
"""
find_missing.py
---------------
Compare 1:1 mapped parquet files across two S3 prefixes and report which
files and row indexes from S3_PATH_1 are missing in S3_PATH_2.

Rows are matched by the value of a key column (default: "vertices").
Only the key column is loaded per chunk, so memory stays flat.

Usage:
    uv run python find_missing.py S3_PATH_1 S3_PATH_2

Example:
    uv run python find_missing.py \\
        s3://toriccy/datasets--calabi-yau-data--polytopes-4d/ \\
        s3://toriccy/4d-polytope-facets/
"""

import argparse
import re

import pyarrow.fs as pafs
import pyarrow.parquet as pq
from tqdm import tqdm

DEFAULT_REGION = "us-east-2"
DEFAULT_KEY_COL = "vertices"


def to_hashable(v):
    if isinstance(v, list):
        return tuple(to_hashable(x) for x in v)
    return v


def list_parquets(s3, path, pattern):
    """Return {matched_key: full_s3_path} for parquet files whose name matches pattern."""
    infos = s3.get_file_info(pafs.FileSelector(path, recursive=True))
    result = {}
    for info in infos:
        if info.type != pafs.FileType.File or not info.path.endswith(".parquet"):
            continue
        m = pattern.search(info.base_name)
        if m:
            result[m.group(1)] = info.path
    return result


def read_key_column(s3, fpath, col, bar):
    """Stream key column row-group by row-group; return {hashable_key: row_index}."""
    keys = {}
    row_offset = 0
    with s3.open_input_file(fpath) as f:
        pf = pq.ParquetFile(f)
        for i in range(pf.metadata.num_row_groups):
            batch = pf.read_row_group(i, columns=[col])
            for local_idx, val in enumerate(batch.column(col).to_pylist()):
                keys[to_hashable(val)] = row_offset + local_idx
            row_offset += batch.num_rows
            bar.update(1)
    return keys


def main():
    ap = argparse.ArgumentParser(
        description="Find rows present in S3_PATH_1 but missing in S3_PATH_2.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("s3_path_1", help="Source S3 prefix")
    ap.add_argument("s3_path_2", help="Derived S3 prefix")
    ap.add_argument(
        "--key-col",
        default=DEFAULT_KEY_COL,
        help="Column used to match rows across files",
    )
    ap.add_argument(
        "--re-1",
        default=r"polytopes-4d-0*(\d+)-vertices\.parquet$",
        help="Regex with one capture group applied to path-1 filenames",
    )
    ap.add_argument(
        "--re-2",
        default=r"4d-polytope-facets-(\d+)-vertices\.parquet$",
        help="Regex with one capture group applied to path-2 filenames (defaults to --re-1)",
    )
    ap.add_argument("--region", default=DEFAULT_REGION, help="AWS region")
    args = ap.parse_args()

    re_1 = re.compile(args.re_1)
    re_2 = re.compile(args.re_2 if args.re_2 else args.re_1)

    s3 = pafs.S3FileSystem(region=args.region)
    path_1 = args.s3_path_1.removeprefix("s3://").rstrip("/")
    path_2 = args.s3_path_2.removeprefix("s3://").rstrip("/")

    print("Listing files…", flush=True)
    files_1 = list_parquets(s3, path_1, re_1)
    files_2 = list_parquets(s3, path_2, re_2)

    missing_files = sorted(name for name in files_1 if name not in files_2)
    matched_files = sorted(name for name in files_1 if name in files_2)

    print(
        f"  {len(files_1)} file(s) in path 1 | "
        f"{len(files_2)} file(s) in path 2 | "
        f"{len(missing_files)} file(s) missing entirely | "
        f"{len(matched_files)} to compare\n"
    )

    # Pre-read footers to get total row group count (reads only metadata, no row data)
    print("Reading file metadata…", flush=True)
    matched_meta = {}
    for name in matched_files:
        with s3.open_input_file(files_1[name]) as f:
            matched_meta[name] = pq.ParquetFile(f).metadata
    # Each matched file is read twice (path 1 + path 2)
    total_rgs = 2 * sum(m.num_row_groups for m in matched_meta.values())

    results = []

    # Files entirely absent from path 2
    for name in missing_files:
        with s3.open_input_file(files_1[name]) as f:
            n = pq.ParquetFile(f).metadata.num_rows
        results.append(
            {"file": name, "missing_rows": list(range(n)), "entirely_absent": True}
        )

    # Files present in both — compare row by row via key column
    with tqdm(total=total_rgs, unit="rg", desc="Comparing") as bar:
        for name in matched_files:
            bar.set_postfix_str(name, refresh=False)
            a_keys = read_key_column(s3, files_1[name], args.key_col, bar)
            b_keys = read_key_column(s3, files_2[name], args.key_col, bar)

            missing_rows = sorted(idx for key, idx in a_keys.items() if key not in b_keys)
            if missing_rows:
                results.append(
                    {"file": name, "missing_rows": missing_rows, "entirely_absent": False}
                )

    # Report
    if not results:
        print("No missing files or rows.")
        return

    total_rows = sum(len(r["missing_rows"]) for r in results)
    print(f"\n{total_rows} missing row(s) across {len(results)} file(s):\n")
    for r in results:
        tag = "ABSENT" if r["entirely_absent"] else "partial"
        rows = r["missing_rows"]
        preview = str(rows[:10]) + ("…" if len(rows) > 10 else "")
        print(f"  [{tag}] {r['file']}: {len(rows)} missing row(s)  {preview}")


if __name__ == "__main__":
    main()
