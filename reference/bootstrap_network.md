# Bootstrap a Network Estimate

Non-parametric bootstrap for any network estimated by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
Works with all built-in methods (transition and association) as well as
custom registered estimators.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses a fast pre-computation strategy: per-sequence count matrices are
computed once, and each bootstrap iteration only resamples sequences via
`colSums` (C-level) plus lightweight post-processing. Data must be in
wide format for transition bootstrap; use
[`convert_sequence_format`](https://saqr.me/Nestimate/reference/convert_sequence_format.md)
to convert long-format data first.

For association methods (`"cor"`, `"pcor"`, `"glasso"`, and custom
estimators), the full estimator is called on resampled rows each
iteration.

If a transition network contains only one sequence, the function warns
that such a network is not recommended for bootstrap or other
confirmatory testing.

## Usage

``` r
bootstrap_network(
  x,
  iter = 1000L,
  ci_level = 0.05,
  inference = "stability",
  consistency_range = c(0.75, 1.25),
  edge_threshold = NULL,
  seed = NULL,
  boundary = c("inclusive", "strict")
)
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).
  The data, method, params, scaling, threshold, and level are all
  extracted from this object.

- iter:

  Integer. Number of bootstrap iterations (default: 1000).

- ci_level:

  Numeric. Significance level for CIs and p-values (default: 0.05).

- inference:

  Character. `"stability"` (default) tests whether bootstrap replicates
  fall within a multiplicative consistency range around the original
  weight. `"threshold"` tests whether replicates exceed a fixed edge
  threshold.

- consistency_range:

  Numeric vector of length 2. Multiplicative bounds for stability
  inference (default: `c(0.75, 1.25)`).

- edge_threshold:

  Numeric or NULL. Fixed threshold for `inference = "threshold"`. If
  NULL, defaults to the 10th percentile of absolute original edge
  weights.

- seed:

  Integer or NULL. RNG seed for reproducibility.

- boundary:

  Character. Comparison rule when computing the consistency-range
  p-value. `"inclusive"` (default, tna-compatible) counts iterations
  that meet the bound (\\\le\\ / \\\ge\\); `"strict"` counts only
  iterations strictly outside (\\\<\\ / \\\>\\).

## Value

An object of class `"net_bootstrap"` containing:

- original:

  The original `netobject`.

- mean:

  Bootstrap mean weight matrix.

- sd:

  Bootstrap SD matrix.

- p_values:

  P-value matrix.

- significant:

  Original weights where p \< ci_level, else 0.

- ci_lower:

  Lower CI bound matrix.

- ci_upper:

  Upper CI bound matrix.

- cr_lower:

  Consistency range lower bound (stability only).

- cr_upper:

  Consistency range upper bound (stability only).

- summary:

  Long-format data frame of edge-level statistics.

- model:

  Pruned `netobject` (non-significant edges zeroed).

- method, params, iter, ci_level, inference:

  Bootstrap config.

- consistency_range, edge_threshold:

  Inference parameters.

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`print.net_bootstrap`](https://saqr.me/Nestimate/reference/print.net_bootstrap.md),
[`summary.net_bootstrap`](https://saqr.me/Nestimate/reference/summary.net_bootstrap.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
  method = "relative")
boot <- bootstrap_network(net, iter = 10)
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
boot <- bootstrap_network(net, iter = 100)
print(boot)
#> Bootstrap Network  [Transition Network (relative) | directed]
#>   Iterations : 100  |  Nodes : 4
#>   Edges      : 0 significant / 16 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
summary(boot)
#>    from to    weight       mean         sd   p_value   sig   ci_lower  ci_upper
#> 1     A  A 0.2380952 0.22549744 0.10350003 0.5742574 FALSE 0.00000000 0.4091452
#> 2     A  B 0.2380952 0.24114731 0.09590287 0.5049505 FALSE 0.10000000 0.4523864
#> 3     A  C 0.1904762 0.19542089 0.08640450 0.5940594 FALSE 0.04449405 0.3481731
#> 4     A  D 0.3333333 0.33793436 0.11154625 0.5247525 FALSE 0.15230179 0.5402206
#> 5     B  A 0.1304348 0.13866145 0.07187177 0.6336634 FALSE 0.04252717 0.3107143
#> 6     B  B 0.1739130 0.16900511 0.06392779 0.5544554 FALSE 0.02065217 0.2663462
#> 7     B  C 0.3043478 0.30361963 0.09915202 0.4752475 FALSE 0.11764368 0.5000000
#> 8     B  D 0.3913043 0.38871381 0.11498719 0.4158416 FALSE 0.17010870 0.5724256
#> 9     C  A 0.2105263 0.22757228 0.11419320 0.7227723 FALSE 0.02159091 0.4735294
#> 10    C  B 0.1052632 0.09744051 0.05974522 0.6633663 FALSE 0.00000000 0.2000000
#> 11    C  C 0.2631579 0.25532349 0.09327931 0.4455446 FALSE 0.06153846 0.4362092
#> 12    C  D 0.4210526 0.41966372 0.10250544 0.3366337 FALSE 0.23424908 0.6133547
#> 13    D  A 0.2962963 0.30959223 0.08589135 0.4059406 FALSE 0.15967433 0.4642308
#> 14    D  B 0.2592593 0.25870020 0.07826515 0.3960396 FALSE 0.12691532 0.4136798
#> 15    D  C 0.2592593 0.26226022 0.08068165 0.3663366 FALSE 0.12895833 0.4268519
#> 16    D  D 0.1851852 0.16944736 0.07886852 0.5544554 FALSE 0.03771368 0.3337131
#>      cr_lower  cr_upper
#> 1  0.17857143 0.2976190
#> 2  0.17857143 0.2976190
#> 3  0.14285714 0.2380952
#> 4  0.25000000 0.4166667
#> 5  0.09782609 0.1630435
#> 6  0.13043478 0.2173913
#> 7  0.22826087 0.3804348
#> 8  0.29347826 0.4891304
#> 9  0.15789474 0.2631579
#> 10 0.07894737 0.1315789
#> 11 0.19736842 0.3289474
#> 12 0.31578947 0.5263158
#> 13 0.22222222 0.3703704
#> 14 0.19444444 0.3240741
#> 15 0.19444444 0.3240741
#> 16 0.13888889 0.2314815
# }
```
