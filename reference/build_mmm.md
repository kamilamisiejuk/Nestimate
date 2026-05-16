# Fit a Mixed Markov Model

Discovers latent subgroups with different transition dynamics using
Expectation-Maximization. Each mixture component has its own transition
matrix. Sequences are probabilistically assigned to components.

## Usage

``` r
build_mmm(
  data,
  k = 2L,
  n_starts = 50L,
  max_iter = 200L,
  tol = 1e-06,
  smooth = 0.01,
  seed = NULL,
  covariates = NULL,
  estimator = c("firth", "multinom", "chisq")
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

## Value

An object of class `net_mmm` with components:

- data:

  The full N-row sequence frame used for estimation.

- models:

  List of `netobject`s, one per component. Each component carries the
  rows assigned to that component in its `$data` slot, while its
  transition matrix is the EM-estimated component transition matrix.

- k:

  Number of components.

- mixing:

  Numeric vector of mixing proportions.

- posterior:

  N x k matrix of posterior probabilities.

- assignments:

  Integer vector of hard assignments (1..k).

- quality:

  List: `avepp` (per-class), `avepp_overall`, `entropy`,
  `relative_entropy`, `classification_error`.

- log_likelihood, BIC, AIC, ICL:

  Model fit statistics.

- states:

  Character vector of state names.

## Initial states

The first sequence column has special status: it is read directly as the
per-sequence initial state
(`init_state[i] <- match(raw_data[i, state_cols[1L]], states)`). The
function does **not** scan forward to the first non-missing position,
and it does not apply any `na_syms`-style symbol conversion (unlike
[`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)).
The state vocabulary is built from the unique non-`NA` values across all
columns, so if your data uses a sentinel character such as `"*"` or
`"%"` for missing cells, that sentinel becomes a real state and the
first column reads it as a valid initial state. If you want padded
leading missings to be treated as missing, recode them to `NA` before
calling `build_mmm()` (then
[`match()`](https://rdrr.io/r/base/match.html) returns `NA`, which the
EM treats as an uninformative initial distribution), or left-trim the
leading missings so each sequence's first column carries an observed
state.

## See also

[`compare_mmm`](https://saqr.me/Nestimate/reference/compare_mmm.md),
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
mmm
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -57.596  |  BIC = 173.013  |  AIC = 149.192  |  ICL = 174.931
#>   Quality: AvePP = 0.969  |  Entropy = 0.192  |  Class.Err = 0.0%
#>   Status: did not converge in 10 iterations
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  78.5%  0.971
#>   2        6 (20.0%)   21.5%  0.961
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, seed = 42)
print(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -124.400  |  BIC = 306.621  |  AIC = 282.800  |  ICL = 313.645
#>   Quality: AvePP = 0.904  |  Entropy = 0.304  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        20 (66.7%)  62.4%  0.896
#>   2        10 (33.3%)  37.6%  0.919
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -124.400  |  BIC = 306.621  |  AIC = 282.800  |  ICL = 313.645
#>   Quality: AvePP = 0.904  |  Entropy = 0.304  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        20 (66.7%)  62.4%  0.896
#>   2        10 (33.3%)  37.6%  0.919
#> 
#> --- Cluster 1 (62.4%, n=20) ---
#>       A     B     C
#> A 0.001 0.568 0.431
#> B 0.485 0.203 0.313
#> C 0.189 0.629 0.181
#> 
#> --- Cluster 2 (37.6%, n=10) ---
#>       A     B     C
#> A 0.593 0.261 0.145
#> B 0.149 0.397 0.454
#> C 0.394 0.002 0.604
#> 
#>   component     prior n_assigned mean_posterior     avepp
#> 1         1 0.6241614         20      0.8959195 0.8959195
#> 2         2 0.3758386         10      0.9193350 0.9193350
# }
```
