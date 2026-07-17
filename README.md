# ksdb-facets

Compute PALP GL(3,Z) and GL(4,Z) normal forms for every facet of every polytope in the [Kreuzer-Skarke database](https://hf.co/datasets/calabi-yau-data/polytopes-4d) of 4D reflexive lattice polytopes.

For each polytope, the output records:
- The 3D facet vertices (one per facet, in dual-vertex order)
- The **facet normal form** (`facet_nf`): GL(3,Z) normal form of the facet as an abstract 3D lattice polytope
- The **maximal cone normal form** (`maximal_cone_nf`): GL(4,Z) normal form of the cone `conv(facet ∪ {0})`, a finer invariant that captures the facet's embedding in Z⁴

## Components

| Path | Language | Role |
|------|----------|------|
| `compute_facets/compute_facets.cpp` | C++ | Core binary — reads one JSON object per stdin line, emits one per stdout line |
| `scripts/run_compute_facets.py` | Python | Drives the binary over a set of Parquet files with a persistent parallel worker pool |
| `scripts/find_missing.py` | Python | Compares two S3 prefixes and reports files/rows present in one but absent in the other |
| `scripts/rename_column.py` | Python | Renames a Parquet column in-place across all files under an S3 prefix |

## Requirements

- **C++ build**: `g++` (C++11); `GlobalP.h` (PALP header) is bundled in `compute_facets/`
- **Python**: Python ≥ 3.12 and [uv](https://github.com/astral-sh/uv)

## Build

```bash
cd compute_facets
make
```

This produces `compute_facets/compute_facets`.

## Usage

### Python runner (recommended)

```bash
# Full Kreuzer-Skarke dataset from HuggingFace
uv run python scripts/run_compute_facets.py \
    "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
    facet_results/

# Quick test on the first 1000 polytopes
uv run python scripts/run_compute_facets.py --limit 1000 \
    "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
    facet_results_test/
```

Each input file produces one output file in the output directory with `-facets` appended to the stem (e.g. `train-00000-of-00042.parquet` → `train-00000-of-00042-facets.parquet`). Output files are compressed with zstd at level 3.

**Options**

| Flag | Default | Description |
|------|---------|-------------|
| `--verts-col NAME` | `vertices` | Column name for vertex data in the source Parquet |
| `--bin PATH` | `compute_facets/compute_facets` | Path to the compiled binary |
| `--workers N` | `cpu_count()` | Number of parallel worker processes |
| `--batch-size N` | `10000` | Rows loaded into memory per iteration |
| `--limit N` | — | Stop after N rows total (for testing) |
| `--s3-uri URI` | — | Upload each output file (and the log) to this S3 directory, then shut down the instance |
| `--log PATH` | — | Append worker warnings to this file |

### EC2 spot instance use

```bash
uv run python scripts/run_compute_facets.py \
    "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
    /tmp/facet_results/ \
    --s3-uri s3://my-bucket/facet_results \
    --log /tmp/facets.log
```

When `--s3-uri` is set, each output file is uploaded to the S3 directory as it completes, the log file is uploaded before shutdown, and the instance terminates itself via `sudo shutdown -h now`. The EC2 instance role needs `s3:PutObject` on the output bucket.

### Binary directly

The binary reads one JSON object per line from stdin and writes one per line to stdout:

```bash
echo '{"index":0,"verts":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,-1]]}' \
    | ./compute_facets/compute_facets
```

Both `"verts"` (JSON array) and `"NVERTS"` (Mathematica-style `{{...}}` string) are accepted. The `"index"` key is optional and passed through to the output unchanged.

**Output** (one JSON per line):

```json
{
  "index": 0,
  "verts": [[1,0,0,0], ...],
  "dual_verts": [[1,1,1,1], ...],
  "facets": [[[0,0,0], ...], ...],
  "facet_nfs": [[[0,0,0],[1,0,0],[0,1,0],[0,0,1]], ...],
  "maximal_cone_nfs": [[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]], ...]
}
```

## Output Parquet schema

| Column | Type | Description |
|--------|------|-------------|
| `index` | `int64` | Row index within the source Parquet file (0-based) |
| `verts` | `list<list<int32>>` | Original 4D polytope vertices |
| `dual_verts` | `list<list<int32>>` | Facet inner normals; `dual_verts[i]` corresponds to `facets[i]` |
| `facets` | `list<list<list<int32>>>` | Vertex coordinates of each 3D facet |
| `facet_nfs` | `list<list<list<int32>>>` | GL(3,Z) normal form of each facet as an abstract 3D lattice polytope |
| `maximal_cone_nfs` | `list<list<list<int32>>>` | GL(4,Z) normal form of `conv(facets[i] ∪ {0})` |

## Install Python dependencies

```bash
uv sync
```
