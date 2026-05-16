# Build MCML from Raw Transition Data

Builds a Multi-Cluster Multi-Level (MCML) model from raw transition data
(edge lists or sequences) by recoding node labels to cluster labels and
counting actual transitions. Unlike
[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
which aggregates a pre-computed weight matrix, this function works from
the original transition data to produce the TRUE Markov chain over
cluster states.

## Usage

``` r
build_mcml(
  x,
  clusters = NULL,
  method = c("sum", "mean", "median", "max", "min", "density", "geomean"),
  type = c("tna", "frequency", "cooccurrence", "semi_markov", "raw"),
  directed = TRUE,
  compute_within = TRUE,
  actor = NULL,
  action = NULL,
  time = NULL,
  order = NULL,
  session = NULL,
  time_threshold = 900,
  labels = NULL
)
```

## Arguments

- x:

  Input data. Accepts multiple formats:

  data.frame with from/to columns

  :   Edge list. Columns named from/source/src/v1/node1/i and
      to/target/tgt/v2/node2/j are auto-detected. Optional weight column
      (weight/w/value/strength).

  data.frame without from/to columns

  :   Sequence data. Each row is a sequence, columns are time steps.
      Consecutive pairs (t, t+1) become transitions.

  tna object

  :   If `x$data` is non-NULL, uses sequence path on the raw data.
      Otherwise falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  netobject

  :   If `x$data` is non-NULL, detects edge list vs sequence data.
      Otherwise falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  cluster_summary

  :   Returns as-is.

  square numeric matrix

  :   Falls back to
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md).

  non-square or character matrix

  :   Treated as sequence data.

- clusters:

  Cluster/group assignments. Accepts:

  named list

  :   Direct mapping. List names = cluster names, values = character
      vectors of node labels. Example:
      `list(A = c("N1","N2"), B = c("N3","N4"))`

  data.frame

  :   A data frame where the first column contains node names and the
      second column contains group/cluster names. Example:
      `data.frame(node = c("N1","N2","N3"), group = c("A","A","B"))`

  membership vector

  :   Character or numeric vector. Node names are extracted from the
      data. Example: `c("A","A","B","B")`

  column name string

  :   For edge list data.frames, the name of a column containing cluster
      labels. The mapping is built from unique (node, group) pairs in
      both from and to columns. **Limitation:** this mode assigns the
      row's group label to *both* endpoints, so it only makes sense for
      edge lists where source and target nodes always share the same
      group (within-group edges only). For general edge lists where a
      single node may be source in some rows and target in others, or
      where source and target belong to different groups, pass an
      explicit named list (`list(G1 = c("N1","N2"), ...)`) or a
      two-column data frame `data.frame(node, group)` instead.

  NULL

  :   Auto-detect from `netobject$nodes` or `$node_groups` (same logic
      as
      [`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)).

- method:

  Aggregation method for combining edge weights: "sum", "mean",
  "median", "max", "min", "density", "geomean". Default "sum". For raw
  sequence/event-log inputs the function is counting observed
  transitions, so `"sum"` is the only interpretation that preserves the
  count semantics – the other methods are useful when aggregating
  weighted edge lists or pre-existing weight matrices, where each row
  already represents a measurement rather than a single observation.

- type:

  Post-processing: "tna" (row-normalize), "cooccurrence" (symmetrize),
  "semi_markov", or "raw". Default "tna".

- directed:

  Logical. If `TRUE` (default), treat transitions as directed. If
  `FALSE`, symmetrize sequence- and edge-derived weights before
  returning raw/frequency weights or before row-normalizing transition
  probabilities.

- compute_within:

  Logical. Compute within-cluster matrices? Default TRUE.

- actor, action, time, order, session, time_threshold:

  Long-format event-log shortcut. When `action` is supplied on a
  data.frame input, the data is passed through
  [`prepare()`](https://saqr.me/Nestimate/reference/prepare.md) to
  derive a wide sequence, which is then routed to the existing sequence
  path. Behaves identically to
  `prepare(...) |> build_network() |> build_mcml()`.

- labels:

  Optional name -\> label remap applied to within-cluster nodes (the
  macro layer is left untouched because its labels are cluster names).
  Accepts a 2-column data.frame `(name, label)`, a named character
  vector `c(name = "label")`, or a named list. Unmapped names pass
  through unchanged.

## Value

A `cluster_summary` object with `meta$source = "transitions"`, fully
compatible with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`as_tna()`](https://saqr.me/Nestimate/reference/as_tna.md), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## See also

[`cluster_summary`](https://saqr.me/Nestimate/reference/cluster_summary.md)
for matrix-based aggregation, `net_as_tna()` to convert to tna objects,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
visualization

## Examples

``` r
# Edge list with clusters
edges <- data.frame(
  from = c("A", "A", "B", "C", "C", "D"),
  to   = c("B", "C", "A", "D", "D", "A"),
  weight = c(1, 2, 1, 3, 1, 2)
)
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
cs <- build_mcml(edges, clusters)
cs$macro$weights
#>           G1        G2
#> G1 0.5000000 0.5000000
#> G2 0.3333333 0.6666667

# Sequence data with clusters
seqs <- data.frame(
  T1 = c("A", "C", "B"),
  T2 = c("B", "D", "A"),
  T3 = c("C", "C", "D"),
  T4 = c("D", "A", "C")
)
cs <- build_mcml(seqs, clusters, type = "raw")
cs$macro$weights
#>    G1 G2
#> G1  2  2
#> G2  1  4
```
