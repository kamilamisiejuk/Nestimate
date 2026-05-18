# Centrality Stability Coefficient (CS-coefficient)

Estimates the stability of centrality indices under case-dropping. For
each drop proportion, sequences are randomly removed and the network is
re-estimated. The correlation between the original and subset centrality
values is computed. The CS-coefficient is the maximum proportion of
cases that can be dropped while maintaining a correlation above
`threshold` in at least `certainty` of bootstrap samples.

For transition methods, uses pre-computed per-sequence count matrices
for fast resampling. Strength centralities (InStrength, OutStrength) are
computed directly from the matrix without igraph.

## Usage

``` r
centrality_stability(
  x,
  measures = c("InStrength", "OutStrength", "Betweenness"),
  iter = 1000L,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  threshold = 0.7,
  certainty = 0.95,
  method = "pearson",
  centrality_fn = NULL,
  loops = FALSE,
  seed = NULL
)
```

## Arguments

- x:

  A `netobject` from
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).

- measures:

  Character vector. Centrality measures to assess. Built-in:
  `"InStrength"`, `"OutStrength"`, `"Betweenness"`, `"InCloseness"`,
  `"OutCloseness"`, `"Closeness"`. `"Closeness"` is defined only for
  undirected networks; `"InCloseness"`/`"OutCloseness"` only for
  directed networks (requesting the wrong one for the network's
  directedness is an error). Custom measures beyond these are valid only
  when a `centrality_fn` is supplied to resolve them. Default:
  `c("InStrength", "OutStrength", "Betweenness")`.

- iter:

  Integer. Number of bootstrap iterations per drop proportion (default:
  1000).

- drop_prop:

  Numeric vector. Proportions of cases to drop (default:
  `seq(0.1, 0.9, by = 0.1)`).

- threshold:

  Numeric. Minimum correlation to consider stable (default: 0.7).

- certainty:

  Numeric. Required proportion of iterations above threshold (default:
  0.95).

- method:

  Character. Correlation method: `"pearson"`, `"spearman"`, or
  `"kendall"` (default: `"pearson"`).

- centrality_fn:

  Optional function. A custom centrality function that takes a weight
  matrix and returns a named list of centrality vectors. When `NULL`
  (default), all built-in measures are computed internally:
  `"InStrength"`/`"OutStrength"` via `colSums`/`rowSums`, and
  `"Betweenness"`/ `"InCloseness"`/`"OutCloseness"`/`"Closeness"` via an
  internal Floyd-Warshall shortest-path routine. When provided, the
  function is called as `centrality_fn(mat)` and is used only for
  requested measures that are not one of the six built-ins; it should
  return a named list (e.g., `list(my_metric = ...)`).

- loops:

  Logical. If `FALSE` (default), self-loops (diagonal) are excluded from
  centrality computation. This does not modify the stored matrix.

- seed:

  Integer or NULL. RNG seed for reproducibility.

## Value

An object of class `"net_stability"` containing:

- cs:

  Named numeric vector of CS-coefficients per measure.

- correlations:

  Named list of matrices (iter x n_prop) of correlation values per
  measure.

- measures:

  Character vector of measures assessed.

- drop_prop:

  Drop proportions used.

- threshold:

  Stability threshold.

- certainty:

  Required certainty level.

- iter:

  Number of iterations.

- method:

  Correlation method.

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`network_reliability`](https://saqr.me/Nestimate/reference/network_reliability.md)

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#> Warning: All centrality measures have zero variance. No stability can be assessed.
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
cs <- centrality_stability(net, iter = 100, seed = 42,
  measures = c("InStrength", "OutStrength"))
print(cs)
#> Centrality Stability (100 iterations, threshold = 0.7)
#>   Drop proportions: 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9
#> 
#>   CS-coefficients:
#>     InStrength       0.30
#>     OutStrength      0.00
# }
```
