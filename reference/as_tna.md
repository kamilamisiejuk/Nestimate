# Convert cluster_summary to tna Objects

Converts a `cluster_summary` object to proper tna objects that can be
used with all functions from the tna package. Creates both a
between-cluster tna model (cluster-level transitions) and within-cluster
tna models (internal transitions within each cluster).

## Usage

``` r
as_tna(x)

# S3 method for class 'mcml'
as_tna(x)

# Default S3 method
as_tna(x)
```

## Arguments

- x:

  A `cluster_summary` object created by
  [`cluster_summary`](https://mohsaqr.github.io/Nestimate/reference/cluster_summary.md).
  The aggregated weights are passed to
  [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md), which
  row-normalises them as needed.

## Value

A `cluster_tna` object (S3 class) containing:

- between:

  A tna object representing cluster-level transitions. Contains
  `$weights` (k x k transition matrix), `$inits` (initial distribution),
  and `$labels` (cluster names). Use this for analyzing how
  learners/entities move between high-level groups or phases.

- within:

  Named list of tna objects, one per cluster. Each tna object represents
  internal transitions within that cluster. Contains `$weights` (n_i x
  n_i matrix), `$inits` (initial distribution), and `$labels` (node
  labels). Clusters with single nodes or zero-row nodes are excluded
  (tna requires positive row sums).

A `netobject_group` with data preserved from each sub-network.

A `tna` object constructed from the input.

## Details

This is the final step in the MCML workflow, enabling full integration
with the tna package for centrality analysis, bootstrap validation,
permutation tests, and visualization.

### Requirements

The tna package must be installed. If not available, the function throws
an error with installation instructions.

### Workflow


    # Full MCML workflow
    net <- build_network(data, method = "relative")
    net$nodes$clusters <- group_assignments
    cs <- cluster_summary(net)
    tna_models <- as_tna(cs)

    # Now use tna package functions
    plot(tna_models$macro)
    tna::centralities(tna_models$macro)
    tna::bootstrap(tna_models$macro, iter = 1000)

    # Analyze within-cluster patterns
    plot(tna_models$clusters$ClusterA)
    tna::centralities(tna_models$clusters$ClusterA)

### Excluded Clusters

A within-cluster tna cannot be created when:

- The cluster has only 1 node (no internal transitions possible)

- Some nodes in the cluster have no outgoing edges (row sums to 0)

These clusters are silently excluded from `$clusters`. The
between-cluster model still includes all clusters.

## See also

[`cluster_summary`](https://mohsaqr.github.io/Nestimate/reference/cluster_summary.md)
to create the input object,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
visualization without conversion,
[`tna::tna`](http://sonsoles.me/tna/reference/build_model.md) for the
underlying tna constructor

## Examples

``` r
mat <- matrix(runif(36), 6, 6)
rownames(mat) <- colnames(mat) <- LETTERS[1:6]
clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"), G3 = c("E", "F"))
cs <- cluster_summary(mat, clusters)
tna_models <- as_tna(cs)
tna_models
#> Group Networks (4 groups)
#> 
#>   Group  Nodes  Edges  Weights
#>   macro  3      9      [1.259, 3.086]
#>   G1     2      4      [0.196, 0.678]
#>   G2     2      4      [0.432, 0.696]
#>   G3     2      4      [0.389, 0.971]
tna_models$macro$weights
#>          G1       G2       G3
#> G1 1.566271 1.518955 1.258863
#> G2 1.398476 2.295576 2.939396
#> G3 3.086310 1.911803 2.661667
```
