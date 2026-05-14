# PRIDE

PRIDE (Pivot-Ring Indexable Diversified Exploration) is a diversified similarity retrieval approach designed for Relational Database Management Systems (RDBMSs), particularly PostgreSQL.

This repository contains PostgreSQL implementations of diversified similarity queries:

- Range queries (`PRIDE_r`)
- k-Nearest Neighbors queries (`PRIDE_k`)

using different execution strategies:

- Sequential scan
- Cube + GiST indexing
- Omni-based pruning using arrays
- Omni-based pruning using one column per pivot

The implementations are based on pivot-ring radicals, where each object is mapped to a discrete metric-space region according to its distances to a set of pivots.

---

# Requirements

The implementations based on:
- `cube` feature vectors;
- GiST indexing;
- cube distance operators;

require the PostgreSQL `cube` extension.

To enable it:

```sql
CREATE EXTENSION IF NOT EXISTS cube;
```

---

# Dataset Structure

PRIDE assumes a dataset containing:
- an object identifier;
- `float8[]` feature vectors for sequential and Omni-based implementations;
- `cube` feature vectors for GiST-indexed implementations.

Example table:

```sql
create table data_table  (
    id serial primary key,
    features_array float8[] not null,
    features_cube cube not null
);
```

Recommended GiST index:

```sql
create index idx_data_table_gist
on data_table
using gist (features_cube);
```

---

# Intrinsic Dimensionality Estimation

Before selecting pivots, the intrinsic dimensionality (`D2`) of the dataset must be estimated.

The intrinsic dimensionality measures how the data is distributed in the metric space and is used to define the recommended number of pivots (`h`).

After estimating `D2`, the repository uses the following heuristic:

```math
h_{recommended} = \lceil D2 \rceil + 1
```

where:

- `D2` = estimated intrinsic dimensionality of the dataset;
- `h_recommended` = recommended number of pivots.

The computed `h_recommended` value is then used during the pivot selection phase.

# Pivot Selection

PRIDE requires a precomputed set of pivots.

The pivots are selected following the same pivot-selection strategy used by the OMNI-family approaches.

After the intrinsic dimensionality (`D2`) estimation, the recommended number of pivots (`h_recommended`) is used during the pivot selection phase.

The selected pivots are stored in a dedicated pivot table:

```sql
create table pivots_data_table (
    pivot_pos      int primary key,
    obj_id         bigint not null unique,
    features_array float8[] not null,
    features_cube  cube not null
);
```

## Column Description

| Column | Description |
|---|---|
| `pivot_pos` | Fixed pivot ordering (1..h) |
| `obj_id` | Identifier of the original object selected as pivot |
| `features_array` | Array representation of the pivot |
| `features_cube` | Cube representation of the pivot |

---

# Omni Coordinates

Some implementations use Omni-based pruning.

Two Omni representations are provided:

---

## Array-based Omni Coordinates

```sql
create table omni_data_table (
    obj_id       bigint primary key,
    omni_coords  float8[] not null
);
```

Each position in `omni_coords` corresponds to the distance between the object and one pivot.

---

## Column-based Omni Coordinates

This representation creates one column per pivot. For example, for 3 pivots:

```sql
create table omni_data_table_cols (
    obj_id    bigint primary key,
    pivot_1   float8 not null,
    pivot_2   float8 not null,
    pivot_3   float8 not null
);
```

B+-tree indexes can be created over the pivot columns to support pruning during query processing.

# Dataset Metadata

PRIDE requires a precomputed metadata table storing the maximum distance (`max_dist`) between any dataset object and any pivot.

```sql
create table dataset_pivot_metadata (
    data_table    regclass not null,
    pivots_table  regclass not null,
    max_dist      float8 not null,
    primary key (data_table, pivots_table)
);
```
---

# Implementations

## Range Queries

- `PRIDE_r_seq.sql`: Sequential diversified range query;
- `PRIDE_r_cube.sql`: GiST-indexed diversified range query using `cube`;
- `PRIDE_r_omni_array.sql`: Diversified range query using array-based Omni pruning;
- `PRIDE_r_omni_cols.sql`: Diversified range query using column-based Omni pruning.

---

## k-NN Queries

- `PRIDE_k_seq.sql`: Sequential diversified k-NN query;
- `PRIDE_k_cube.sql`: GiST-indexed diversified k-NN query using `cube`;
- `PRIDE_k_omni_array.sql`: Diversified k-NN query using array-based Omni pruning;
- `PRIDE_k_omni_cols.sql`: Diversified k-NN query using column-based Omni pruning.

