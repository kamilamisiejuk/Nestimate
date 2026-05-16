# Permutation Test for Network Comparison

Compares two networks estimated by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
using a permutation test. Works with all built-in methods (transition
and association) as well as custom registered estimators. The test
shuffles which observations belong to which group, re-estimates
networks, and tests whether observed edge-wise differences exceed
chance.

For transition methods (`"relative"`, `"frequency"`, `"co_occurrence"`),
uses a fast pre-computation strategy: per-sequence count matrices are
computed once, and each permutation iteration only shuffles group labels
and computes group-wise `colSums`.

For association methods (`"cor"`, `"pcor"`, `"glasso"`, and custom
estimators), the full estimator is called on each permuted group split.

If either transition network contains only one sequence, the function
warns that such a network is not recommended for permutation or other
confirmatory testing.

## Usage

``` r
permutation(
  x,
  y = NULL,
  iter = 1000L,
  alpha = 0.05,
  paired = FALSE,
  adjust = "none",
  nlambda = 50L,
  seed = NULL
)
```

## Arguments

- x:

  A `netobject` (from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).

- y:

  A `netobject` (from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).
  Must use the same method and have the same nodes as `x`.

- iter:

  Integer. Number of permutation iterations (default: 1000).

- alpha:

  Numeric. Significance level (default: 0.05).

- paired:

  Logical. If `TRUE`, permute within pairs (requires equal number of
  observations in `x` and `y`). Default: FALSE.

- adjust:

  Character. p-value adjustment method passed to
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) (default:
  `"none"`). Common choices: `"holm"`, `"BH"`, `"bonferroni"`.

- nlambda:

  Integer. Number of lambda values for the `glassopath` regularisation
  path (only used when `method = "glasso"`). Higher values give finer
  lambda resolution at the cost of speed. Default: 50.

- seed:

  Integer or NULL. RNG seed for reproducibility.

## Value

An object of class `"net_permutation"` containing:

- x:

  The first `netobject`.

- y:

  The second `netobject`.

- diff:

  Observed difference matrix (`x - y`).

- diff_sig:

  Observed difference where `p < alpha`, else 0.

- p_values:

  P-value matrix (adjusted if `adjust != "none"`).

- effect_size:

  Effect size matrix (observed diff / SD of permutation diffs).

- summary:

  Long-format data frame of edge-level results.

- method:

  The network estimation method.

- iter:

  Number of permutation iterations.

- alpha:

  Significance level used.

- paired:

  Whether paired permutation was used.

- adjust:

  p-value adjustment method used.

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md),
[`print.net_permutation`](https://saqr.me/Nestimate/reference/print.net_permutation.md),
[`summary.net_permutation`](https://saqr.me/Nestimate/reference/summary.net_permutation.md)

## Examples

``` r
s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
n1 <- build_network(s1, method = "relative")
n2 <- build_network(s2, method = "relative")
perm <- permutation(n1, n2, iter = 10)
# \donttest{
set.seed(1)
d1 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
                 V2 = sample(LETTERS[1:4], 20, TRUE),
                 V3 = sample(LETTERS[1:4], 20, TRUE))
d2 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
                 V2 = sample(LETTERS[1:4], 20, TRUE),
                 V3 = sample(LETTERS[1:4], 20, TRUE))
net1 <- build_network(d1, method = "relative")
net2 <- build_network(d2, method = "relative")
perm <- permutation(net1, net2, iter = 100, seed = 42)
print(perm)
#> Permutation Test:Transition Network (relative probabilities) [directed]
#>   Iterations: 100  |  Alpha: 0.05
#>   Nodes: 4  |  Edges tested: 16  |  Significant: 3
summary(perm)
#>    from to   weight_x  weight_y        diff effect_size    p_value   sig
#> 1     A  A 0.33333333 0.2222222  0.11111111   0.6504300 0.49504950 FALSE
#> 2     A  B 0.06666667 0.3333333 -0.26666667  -1.9871473 0.04950495  TRUE
#> 3     A  C 0.20000000 0.2222222 -0.02222222  -0.1329339 0.93069307 FALSE
#> 4     A  D 0.40000000 0.2222222  0.17777778   1.0449465 0.25742574 FALSE
#> 5     B  A 0.30769231 0.3333333 -0.02564103  -0.1106148 0.98019802 FALSE
#> 6     B  B 0.38461538 0.0000000  0.38461538   1.9215279 0.03960396  TRUE
#> 7     B  C 0.23076923 0.2500000 -0.01923077  -0.1103038 0.99009901 FALSE
#> 8     B  D 0.07692308 0.4166667 -0.33974359  -1.9910325 0.07920792 FALSE
#> 9     C  A 0.33333333 0.1111111  0.22222222   1.0147102 0.38613861 FALSE
#> 10    C  B 0.66666667 0.3333333  0.33333333   1.3787271 0.24752475 FALSE
#> 11    C  C 0.00000000 0.2222222 -0.22222222  -1.3964532 0.26732673 FALSE
#> 12    C  D 0.00000000 0.3333333 -0.33333333  -1.8265261 0.12871287 FALSE
#> 13    D  A 0.00000000 0.5000000 -0.50000000  -1.8205956 0.10891089 FALSE
#> 14    D  B 0.00000000 0.2000000 -0.20000000  -0.9828935 0.45544554 FALSE
#> 15    D  C 0.33333333 0.2000000  0.13333333   0.5436513 0.56435644 FALSE
#> 16    D  D 0.66666667 0.1000000  0.56666667   2.5117967 0.01980198  TRUE
# }
```
