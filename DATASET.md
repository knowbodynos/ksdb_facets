# 4D Reflexive Polytope Facets

One row per polytope. A **4D reflexive lattice polytope** has vertices in Z⁴, a unique interior lattice point at the origin, and a dual polytope that is also a lattice polytope. Each row contains the polytope's vertices, its facet decomposition, and two normal-form invariants per facet.

**Format:** Parquet, zstd-compressed. Files are named by vertex count (e.g. `4d-polytope-facets-5-vertices.parquet`).

## Schema

| Column | Type | Shape | Description |
|--------|------|-------|-------------|
| `index` | `int64` | scalar | 0-based row index within the source file |
| `vertices` | `list<list<int32>>` | `(V, 4)` | Polytope vertices in Z⁴ |
| `dual_vertices` | `list<list<int32>>` | `(F, 4)` | Dual polytope vertices; `dual_vertices[i]` is the inner normal of `facets[i]` |
| `facets` | `list<list<list<int32>>>` | `(F, *, 4)` | Vertex coordinates of each 3D facet |
| `facet_nfs` | `list<list<list<int32>>>` | `(F, *, 3)` | GL(3,Z) normal form of each facet |
| `maximal_cone_nfs` | `list<list<list<int32>>>` | `(F, *, 4)` | GL(4,Z) normal form of `conv(facets[i] ∪ {0})` |

`F` = number of facets; `V`, `*` = variable length.

## Column relationships

All list columns are **aligned**: `dual_vertices[i]`, `facets[i]`, `facet_nfs[i]`, and `maximal_cone_nfs[i]` all refer to the same facet.

**`dual_vertices[i]`** is the inner normal of `facets[i]` — the corresponding vertex of the dual polytope — satisfying:

```
⟨dual_vertices[i], v⟩ = −1   for all v ∈ facets[i]
⟨dual_vertices[i], v⟩ ≥ −1   for all v ∈ vertices
```

**`facet_nfs[i]`** is the GL(3,Z) normal form of `facets[i]` projected into its own 3D affine span — an isomorphism invariant of the facet as an abstract lattice polytope, independent of its embedding in Z⁴.

**`maximal_cone_nfs[i]`** is the GL(4,Z) normal form of the cone `conv(facets[i] ∪ {0})`, capturing how the facet embeds in Z⁴. It is a strictly finer invariant: equal `maximal_cone_nf` implies equal `facet_nf`, but not vice versa — two facets can be abstractly isomorphic yet have inequivalent embeddings.

## Loading

### PyArrow

```python
import pyarrow.parquet as pq

table = pq.read_table("4d-polytope-facets-5-vertices.parquet")

# Access one row
row = table.slice(0, 1)
vertices         = row["vertices"][0].as_py()         # list of [x,y,z,w]
dual_vertices    = row["dual_vertices"][0].as_py()    # list of [x,y,z,w], one per facet
facets           = row["facets"][0].as_py()           # list of facets, each a list of [x,y,z,w]
facet_nfs        = row["facet_nfs"][0].as_py()        # list of NFs, each a list of [x,y,z]
maximal_cone_nfs = row["maximal_cone_nfs"][0].as_py() # list of NFs, each a list of [x,y,z,w]

# Verify dual_vertices / facets correspondence for facet 0
n = dual_vertices[0]
for v in facets[0]:
    assert sum(n[j] * v[j] for j in range(4)) == -1
```

Read specific columns only:

```python
table = pq.read_table("4d-polytope-facets-5-vertices.parquet",
                      columns=["index", "facet_nfs", "maximal_cone_nfs"])
```

Stream row-group by row-group:

```python
pf = pq.ParquetFile("4d-polytope-facets-5-vertices.parquet")
for i in range(pf.metadata.num_row_groups):
    batch = pf.read_row_group(i)
```

### Pandas

```python
import pandas as pd

df = pd.read_parquet("4d-polytope-facets-5-vertices.parquet")
facets_row0 = df.loc[0, "facets"]  # nested list columns are Python objects
```

### From S3

```python
import pyarrow.fs as pafs, pyarrow.parquet as pq

s3 = pafs.S3FileSystem(region="us-east-2")
with s3.open_input_file("your-bucket/prefix/4d-polytope-facets-5-vertices.parquet") as f:
    table = pq.read_table(f)
```
