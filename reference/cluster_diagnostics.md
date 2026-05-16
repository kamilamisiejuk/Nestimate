# Cluster Diagnostics

Unified entry point for clustering quality information. Returns a
`net_cluster_diagnostics` object that normalises the diagnostic surface
across distance-based and model-based clusterings – you no longer have
to know which fields live on `net_clustering` vs. `net_mmm` vs. the slim
`net_mmm_clustering` attribute of a `netobject_group`.

## Usage

``` r
cluster_diagnostics(x, ...)

# S3 method for class 'net_cluster_diagnostics'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `net_clustering`, `net_mmm`, `netobject_group` (with
  `attr(, "clustering")` attached by
  [`cluster_network()`](https://saqr.me/Nestimate/reference/cluster_network.md)
  or
  [`cluster_mmm()`](https://saqr.me/Nestimate/reference/cluster_mmm.md)),
  or `net_mmm_clustering`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

- row.names, optional:

  Standard `as.data.frame` arguments (ignored).

## Value

A `net_cluster_diagnostics` object.

## Details

The returned object carries:

- family:

  Either `"distance"` or `"mmm"`.

- k, n, sizes:

  Number of clusters, number of sequences, sizes vector.

- per_cluster:

  A `data.frame` – one row per cluster, columns differ by family.
  Distance: `cluster`, `size`, `pct`, `mean_within_dist`, `sil_mean`.
  MMM: `cluster`, `size`, `pct`, `mix_pct`, `avepp`, `class_err_pct`.

- overall:

  A named list of family-specific summary metrics (`silhouette` for
  distance; `avepp_overall`, `entropy`, `classification_error` for MMM).

- ics:

  For MMM: a list with `BIC`, `AIC`, `ICL`. `NULL` for distance.

- metadata:

  Method / dissimilarity / weighted / lambda etc.

- source:

  The original clustering object, kept by reference so
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) can delegate
  without recomputing anything.

## See also

[`print.net_cluster_diagnostics`](https://saqr.me/Nestimate/reference/print.net_cluster_diagnostics.md),
[`plot.net_cluster_diagnostics`](https://saqr.me/Nestimate/reference/plot.net_cluster_diagnostics.md),
[`compare_mmm`](https://saqr.me/Nestimate/reference/compare_mmm.md) for
k-sweep model selection (MMM only).

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
cl <- build_clusters(seqs, k = 2, method = "ward.D2")
cluster_diagnostics(cl)
#> Cluster Diagnostics (distance) [ward.D2 / hamming]
#>   Sequences: 30  |  Clusters: 2
#>   Quality: silhouette = 0.403
#> 
#>   Cluster  N           Mean within-dist  Silhouette
#>   1        20 (66.7%)  1.137             0.326
#>   2        10 (33.3%)  0.733             0.558
# \donttest{
grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
cluster_diagnostics(grp)
#> Cluster Diagnostics (mmm) [k = 2]
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   Quality: AvePP = 0.989  |  Entropy = 0.090  |  Class.Err = 0.0%
#>   ICs: LL = -61.147  |  BIC = 180.114  |  AIC = 156.293  |  ICL = 180.806
#> 
#>   Cluster  N           Mix%   AvePP  Class.Err%
#>   1        24 (80.0%)  79.4%  0.989   0.0%
#>   2        6 (20.0%)   20.6%  0.987   0.0%
as.data.frame(cluster_diagnostics(grp))
#>   cluster size pct  mix_pct     avepp class_err_pct
#> 1       1   24  80 79.44647 0.9889875             0
#> 2       2    6  20 20.55353 0.9867095             0
# }
```
