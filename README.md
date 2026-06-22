# ksdb-facets

Compute the PALP normal forms of the maximal cones (3D facets joined to the origin) for every polytope in the [Kreuzer-Skarke database](https://hf.co/datasets/calabi-yau-data/polytopes-4d) of 4D reflexive lattice polytopes, and record how many facets share each unique cone type.

## Components

| Path | Language | Role |
|------|----------|------|
| `cone_nf/cone_nf.cpp` | C++ | Core binary — reads one JSON object per stdin line, emits one JSON object per stdout line |
| `run_cone_nf.py` | Python | Drives the binary at scale over a Parquet file with parallel workers |

## Requirements

- **C++ build**: `g++` (C++11) and the PALP header `GlobalP.h` (expected at `modules/facets/code/GlobalP.h` relative to the `cone_nf/` directory)
- **Python**: Python ≥ 3.12 and [uv](https://github.com/astral-sh/uv)

## Build

```bash
cd cone_nf
make
```

This produces `cone_nf/cone_nf`.

## Usage

### Python runner (recommended)

```bash
# Full Kreuzer-Skarke dataset from HuggingFace
uv run python run_cone_nf.py \
    "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
    cone_nf_results.parquet

# Quick test on the first 1000 polytopes
uv run python run_cone_nf.py --limit 1000 \
    "hf://datasets/calabi-yau-data/polytopes-4d/*.parquet" \
    cone_nf_test.parquet
```

**Options**

| Flag | Default | Description |
|------|---------|-------------|
| `--verts-col NAME` | `vertices` | Column name for vertex data in the source Parquet |
| `--bin PATH` | `cone_nf/cone_nf` | Path to the compiled binary |
| `--workers N` | `cpu_count()` | Parallel worker processes |
| `--batch-size N` | `10000` | Rows loaded into memory per iteration |
| `--chunk-size N` | `500` | Rows per worker subprocess call |
| `--limit N` | — | Stop after N rows (for testing) |
| `--s3-uri URI` | — | Upload output to S3 then shut down the instance (EC2 spot use) |

### Binary directly

The binary reads one JSON object per line from stdin:

```bash
# JSON array of 4D integer coordinates
echo '{"POLYID":1,"verts":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1],[-1,-1,-1,-1]]}' \
    | ./cone_nf/cone_nf

# Mathematica-style NVERTS string
echo '{"POLYID":1,"NVERTS":"{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{-1,-1,-1,-1}}"}' \
    | ./cone_nf/cone_nf
```

**Output** (one JSON per line):

```json
{
  "POLYID": 1,
  "cones": [
    { "facet_verts": [[...], ...], "cone_nf": [[...], ...] },
    ...
  ],
  "unique_cones": [
    { "cone_nf": [[...], ...], "multiplicity": 3, "example_facet": [[...], ...] },
    ...
  ]
}
```

## Output Parquet schema

| Column | Type | Description |
|--------|------|-------------|
| `idx` | `int64` | Row index in the source Parquet |
| `vertices` | `list<list<int32>>` | Original polytope vertices (pass-through) |
| `cone_nfs` | `list<list<list<int32>>>` | PALP normal-form vertex matrix for each unique cone type |
| `multiplicities` | `list<int32>` | Number of facets that map to each cone NF (parallel to `cone_nfs`) |

## Install Python dependencies

```bash
uv sync
```
