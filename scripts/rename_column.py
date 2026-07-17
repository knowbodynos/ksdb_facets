#!/usr/bin/env python3
"""
rename_column.py
----------------
Rename a column in every parquet file under an S3 prefix, in-place.
Streams one row group at a time to keep memory flat.

Usage:
    uv run python rename_column.py s3://bucket/prefix/ OLD_COL NEW_COL

Example:
    uv run python rename_column.py s3://toriccy/shared/4d-polytope-facets/ verts vertices
"""

import argparse
import sys

import boto3
import pyarrow.fs as pafs
import pyarrow.parquet as pq
from tqdm import tqdm

DEFAULT_REGION = "us-east-2"


def renamed_schema(schema, old_name, new_name):
    idx = schema.get_field_index(old_name)
    return schema.set(idx, schema.field(idx).with_name(new_name))


def process_file(s3, s3_client, fpath, old_name, new_name, bar):
    """Stream-rename via a .tmp key, then server-side copy+delete. Returns False if column not found."""
    tmp_path = fpath + ".tmp"

    with s3.open_input_file(fpath) as f:
        pf = pq.ParquetFile(f)
        if old_name not in pf.schema_arrow.names:
            bar.update(pf.metadata.num_row_groups)
            return False
        out_schema = renamed_schema(pf.schema_arrow, old_name, new_name)
        n_row_groups = pf.metadata.num_row_groups

        with s3.open_output_stream(tmp_path) as out:
            writer = pq.ParquetWriter(out, out_schema, compression="zstd", compression_level=3)
            for i in range(n_row_groups):
                writer.write_table(pf.read_row_group(i).rename_columns(out_schema.names))
                bar.update(1)
            writer.close()

    bucket, src_key = tmp_path.split("/", 1)
    _, dst_key = fpath.split("/", 1)
    s3_client.copy({"Bucket": bucket, "Key": src_key}, bucket, dst_key)
    s3_client.delete_object(Bucket=bucket, Key=src_key)
    return True


def main():
    ap = argparse.ArgumentParser(description="Rename a parquet column in-place on S3.")
    ap.add_argument("s3_prefix", help="S3 prefix to scan (e.g. s3://bucket/prefix/)")
    ap.add_argument("old_name",  help="Column name to rename from")
    ap.add_argument("new_name",  help="Column name to rename to")
    ap.add_argument("--region",  default=DEFAULT_REGION, help="AWS region")
    ap.add_argument("--dry-run", action="store_true",
                    help="List files that would be changed without writing")
    args = ap.parse_args()

    s3        = pafs.S3FileSystem(region=args.region)
    s3_client = boto3.client("s3", region_name=args.region)
    path      = args.s3_prefix.removeprefix("s3://").rstrip("/")

    infos = s3.get_file_info(pafs.FileSelector(path, recursive=True))
    files = [f.path for f in infos if f.path.endswith(".parquet")]

    if not files:
        sys.exit(f"No parquet files found under s3://{path}")

    print(f"Found {len(files)} parquet file(s) under s3://{path}")
    if args.dry_run:
        print("(dry run — no files will be written)\n")

    # Read footers only (no row data) to get total row group count for the progress bar
    print("Reading file metadata…", flush=True)
    file_meta = []
    for fpath in files:
        with s3.open_input_file(fpath) as f:
            meta = pq.ParquetFile(f).metadata
        file_meta.append((fpath, meta))
    total_rgs = sum(m.num_row_groups for _, m in file_meta)

    skipped = 0
    with tqdm(total=total_rgs, unit="rg") as bar:
        for fpath, meta in file_meta:
            if args.dry_run:
                has_col = args.old_name in meta.schema.to_arrow_schema().names
                if has_col:
                    tqdm.write(f"  would rename '{args.old_name}' → '{args.new_name}' in {fpath}")
                else:
                    skipped += 1
                bar.update(meta.num_row_groups)
                continue

            if not process_file(s3, s3_client, fpath, args.old_name, args.new_name, bar):
                skipped += 1

    if skipped:
        print(f"\nSkipped {skipped} file(s) that did not have column '{args.old_name}'.")
    if not args.dry_run:
        print(f"\nDone. Renamed '{args.old_name}' → '{args.new_name}' in {len(files) - skipped} file(s).")


if __name__ == "__main__":
    main()
