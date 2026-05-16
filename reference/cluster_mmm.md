# Cluster sequences using Mixed Markov Models

Fits a mixture of Markov chains to sequence data and returns a
`netobject_group` containing per-cluster transition networks. This is
the MMM equivalent of
[`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md)
(which uses distance-based clustering); both functions share the
`cluster_by = ...` surface argument so the call shape stays uniform
across clustering families.

## Usage

``` r
cluster_mmm(
  data,
  k = 2L,
  n_starts = 50L,
  max_iter = 200L,
  tol = 1e-06,
  smooth = 0.01,
  seed = NULL,
  covariates = NULL,
  estimator = c("firth", "multinom", "chisq"),
  cluster_by = "mmm",
  ...
)
```

## Arguments

- data:

  A data.frame (wide format), `netobject`, or `tna` model. For tna
  objects, extracts the stored data.

- k:

  Integer. Whole finite number of mixture components, \>= 2. Default: 2.

- n_starts:

  Integer. Positive whole finite number of random restarts. Default: 50.

- max_iter:

  Integer. Positive whole finite maximum EM iterations per start.
  Default: 200.

- tol:

  Numeric. Finite positive convergence tolerance. Default: 1e-6.

- smooth:

  Numeric. Finite non-negative Laplace smoothing constant. Default:
  0.01.

- seed:

  Integer or NULL. Random seed.

- covariates:

  Optional. Covariates integrated into the EM algorithm to model
  covariate-dependent mixing proportions. Accepts a string, character
  vector, formula, or data.frame (same forms as
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)).
  For `netobject` or `cograph_network` input, names are resolved against
  `$metadata` first, so a typical call is
  `build_mmm(net, k = 3, covariates = "session_label")`. Unlike the
  post-hoc analysis in
  [`build_clusters()`](https://saqr.me/Nestimate/reference/build_clusters.md),
  these covariates directly influence cluster membership during EM
  estimation.

- estimator:

  Multinomial fitter for the post-hoc covariate analysis (does not
  affect EM): `"firth"` (default, via
  [`brglm2::brmultinom`](https://rdrr.io/pkg/brglm2/man/brmultinom.html);
  finite under separation), `"multinom"`
  ([`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html); warns
  about separation risk), or `"chisq"` (descriptive tests, no logit).
  See
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)
  for full details.

- cluster_by:

  Character. Accepted only as `"mmm"` (the default). Present so
  `cluster_mmm()` and
  [`cluster_network()`](https://saqr.me/Nestimate/reference/cluster_network.md)
  share the same call shape; any other value raises an error pointing at
  [`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md).

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `netobject_group` (list of `netobject`s, one per cluster).
MMM-specific information is stored in `attr(, "clustering")` (class
`"net_mmm_clustering"`):

- assignments:

  Integer vector of cluster assignments.

- k:

  Number of clusters.

- posterior:

  N x k matrix of posterior probabilities.

- mixing:

  Mixing proportions.

- quality:

  List with AvePP, entropy, classification error.

- BIC, AIC, ICL:

  Model fit statistics.

- data:

  The full N-row sequence frame, matching `$assignments` – so
  [`sequence_plot`](https://saqr.me/Nestimate/reference/sequence_plot.md)
  and
  [`distribution_plot`](https://saqr.me/Nestimate/reference/distribution_plot.md)
  can recover both.

## Details

For the full `net_mmm` object with posterior probabilities, model fit
statistics, and S3 methods, use
[`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md) instead.

## See also

[`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md) for the
full MMM object,
[`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md)
for distance-based clustering

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
grp[[1]]$weights
#>           A         B          C
#> A 0.1089275 0.8624067 0.02866584
#> B 0.4564744 0.4564744 0.08705128
#> C 0.1384553 0.5770726 0.28447209
attr(grp, "clustering")$assignments
#>  [1] 1 2 1 2 1 1 1 1 1 1 1 2 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1
# \donttest{
# Visualise with sequence_plot
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 40, TRUE),
  V2 = sample(LETTERS[1:3], 40, TRUE),
  V3 = sample(LETTERS[1:3], 40, TRUE)
)
grp <- cluster_mmm(seqs, k = 2)
sequence_plot(grp, type = "index")

# }
```
