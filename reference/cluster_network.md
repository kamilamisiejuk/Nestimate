# Cluster data and build per-cluster networks in one step

Combines sequence clustering and network estimation into a single call.
Clusters the data using the specified algorithm, then calls
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
on each cluster subset.

## Usage

``` r
cluster_network(data, k, cluster_by = "pam", dissimilarity = "hamming", ...)
```

## Arguments

- data:

  Sequence data. Accepts a data frame, matrix, or `netobject`. See
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)
  for supported formats.

- k:

  Integer. Number of clusters.

- cluster_by:

  Character. Clustering algorithm passed to
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)'s
  `method` parameter (`"pam"`, `"ward.D2"`, `"ward.D"`, `"complete"`,
  `"average"`, `"single"`, `"mcquitty"`, `"median"`, `"centroid"`), or
  `"mmm"` for Mixed Markov Model clustering. Default: `"pam"`.

- dissimilarity:

  Character. Distance metric for sequence clustering. Only valid when
  `cluster_by != "mmm"`. Default: `"hamming"`.

- ...:

  Routed to two stages. For distance clustering (`cluster_by != "mmm"`),
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)
  arguments `na_syms`, `weighted`, `lambda`, `seed`, `q`, `p`, and
  `covariates` are intercepted and forwarded to the clusterer;
  recognised network-estimation arguments flow to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  Unknown arguments error before either stage is run. When
  `cluster_by = "mmm"`, the recognised
  [`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md)
  arguments (`n_starts`, `max_iter`, `tol`, `smooth`, `seed`,
  `covariates`) are intercepted instead, and the rest flows to
  `build_network`. In both modes, when `data` is a `netobject`, its
  `build_args` are merged into the `build_network` side (caller's
  explicit values take precedence).

## Value

A `netobject_group`.

## Details

If `data` is a `netobject` and `method` is not provided in `...`, the
original network method is inherited automatically so the per-cluster
networks match the type of the input network.

## See also

[`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md),
[`cluster_mmm`](https://saqr.me/Nestimate/reference/cluster_mmm.md),
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
grp <- cluster_network(seqs, k = 2)
grp
#> Group Networks (2 clusters via pam / hamming)
#> 
#>   Group      Nodes  Edges  Weights         N
#>   Cluster 1  3      2      [1.000, 1.000]  2 (40.0%)
#>   Cluster 2  3      4      [0.500, 1.000]  3 (60.0%)
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 50, TRUE), V2 = sample(LETTERS[1:4], 50, TRUE),
  V3 = sample(LETTERS[1:4], 50, TRUE), V4 = sample(LETTERS[1:4], 50, TRUE)
)
# Default: PAM clustering, relative (transition) networks
grp <- cluster_network(seqs, k = 3)

# Specify network method (cor requires numeric panel data)
if (FALSE) { # \dontrun{
panel <- as.data.frame(matrix(rnorm(1500), nrow = 300, ncol = 5))
grp <- cluster_network(panel, k = 2, method = "cor")
} # }

# MMM-based clustering
grp <- cluster_network(seqs, k = 2, cluster_by = "mmm")
# }
```
