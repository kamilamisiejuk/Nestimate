# Bootstrap for Regularized Partial Correlation Networks

Fast, single-call bootstrap for EBICglasso partial correlation networks.
Combines nonparametric edge/centrality bootstrap, case-dropping
stability analysis, edge/centrality difference tests, predictability
CIs, and thresholded network into one function. Designed as a faster
alternative to bootnet with richer output.

## Usage

``` r
boot_glasso(
  x,
  iter = 1000L,
  cs_iter = 500L,
  cs_drop = seq(0.1, 0.9, by = 0.1),
  alpha = 0.05,
  gamma = 0.5,
  nlambda = 100L,
  centrality = c("strength", "expected_influence", "betweenness", "closeness"),
  centrality_fn = NULL,
  cor_method = "pearson",
  ncores = 1L,
  seed = NULL
)
```

## Arguments

- x:

  A data frame, numeric matrix (observations x variables), or a
  `netobject` with `method = "glasso"`.

- iter:

  Integer. Number of nonparametric bootstrap iterations (default: 1000).

- cs_iter:

  Integer. Number of case-dropping iterations per drop proportion
  (default: 500).

- cs_drop:

  Numeric vector. Drop proportions for CS-coefficient computation
  (default: `seq(0.1, 0.9, by = 0.1)`).

- alpha:

  Numeric. Significance level for CIs (default: 0.05).

- gamma:

  Numeric. EBIC hyperparameter (default: 0.5).

- nlambda:

  Integer. Number of lambda values in the regularization path (default:
  100).

- centrality:

  Character vector. Centrality measures to compute. Built-in:
  `"strength"`, `"expected_influence"`, `"betweenness"`, `"closeness"`.
  Custom measures beyond these require `centrality_fn`. Default:
  `c("strength", "expected_influence", "betweenness", "closeness")`.

- centrality_fn:

  Optional function. A custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors. When `NULL`
  (default), only `"strength"` and `"expected_influence"` are computed
  via `rowSums`/ `colSums`. When provided, the function is called as
  `centrality_fn(mat)` and should return a named list (e.g.,
  `list(closeness = ..., betweenness = ...)`).

- cor_method:

  Character. Correlation method: `"pearson"` (default), `"spearman"`, or
  `"kendall"`.

- ncores:

  Integer. Number of parallel cores for mclapply (default: 1,
  sequential).

- seed:

  Integer or NULL. RNG seed for reproducibility.

## Value

An object of class `"boot_glasso"` containing:

- original_pcor:

  Original partial correlation matrix.

- original_precision:

  Original precision matrix.

- original_centrality:

  Named list of original centrality vectors.

- original_predictability:

  Named numeric vector of node R-squared.

- edge_ci:

  Data frame of edge CIs (edge, weight, ci_lower, ci_upper, inclusion).

- edge_inclusion:

  Named numeric vector of edge inclusion probabilities.

- thresholded_pcor:

  Partial correlation matrix with non-significant edges zeroed.

- centrality_ci:

  Named list of data frames (node, value, ci_lower, ci_upper) per
  centrality measure.

- cs_coefficient:

  Named numeric vector of CS-coefficients per centrality measure.

- cs_data:

  Data frame of case-dropping correlations (drop_prop, measure,
  correlation).

- edge_diff_p:

  Symmetric matrix of pairwise edge difference p-values.

- centrality_diff_p:

  Named list of symmetric p-value matrices per centrality measure.

- predictability_ci:

  Data frame of node predictability CIs (node, r2, ci_lower, ci_upper).

- boot_edges:

  iter x n_edges matrix of bootstrap edge weights.

- boot_centrality:

  Named list of iter x p bootstrap centrality matrices.

- boot_predictability:

  iter x p matrix of bootstrap R-squared.

- nodes:

  Character vector of node names.

- n:

  Sample size.

- p:

  Number of variables.

- iter:

  Number of nonparametric iterations.

- cs_iter:

  Number of case-dropping iterations.

- cs_drop:

  Drop proportions used.

- alpha:

  Significance level.

- gamma:

  EBIC hyperparameter.

- nlambda:

  Lambda path length.

- centrality_measures:

  Character vector of centrality measures.

- cor_method:

  Correlation method.

- lambda_path:

  Lambda sequence used.

- lambda_selected:

  Selected lambda for original data.

- timing:

  Named numeric vector with timing in seconds.

## See also

[`build_network`](https://mohsaqr.github.io/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://mohsaqr.github.io/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
set.seed(1)
dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
net <- build_network(dat, method = "glasso")
bg <- boot_glasso(net, iter = 10, cs_iter = 5, centrality = "strength")
# \donttest{
set.seed(42)
mat <- matrix(rnorm(60), ncol = 4)
colnames(mat) <- LETTERS[1:4]
net <- build_network(as.data.frame(mat), method = "glasso")
boot <- boot_glasso(net, iter = 100, cs_iter = 50, seed = 42,
  centrality = c("strength", "expected_influence"))
print(boot)
#> GLASSO Bootstrap (100 iterations, 50 case-drop per proportion)
#>   Data: 15 x 4  |  Alpha: 0.05  |  Gamma: 0.50
#>   Edges: 0/6 significant (CI excludes zero)
#>   Mean inclusion probability: 0.32
#> 
#>   Centrality Stability (CS-coefficient):
#>     strength:              0.00 [Unstable]
#>     expected_influence:    0.00 [Unstable]
#> 
#>   Edge differences: 1/15 pairs significantly different
#>   Timing: 0.4s (bootstrap: 0.3s, case-drop: 0.1s)
summary(boot, type = "edges")
#>     edge weight    ci_lower  ci_upper inclusion
#> 1 A -- B      0 -0.15265878 0.3016145      0.24
#> 2 A -- C      0 -0.41694566 0.0000000      0.45
#> 3 B -- C      0 -0.37041174 0.2982661      0.43
#> 4 A -- D      0  0.00000000 0.3870390      0.41
#> 5 B -- D      0 -0.19093168 0.2384451      0.21
#> 6 C -- D      0 -0.01732718 0.2188424      0.18
# }
```
