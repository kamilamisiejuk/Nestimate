# Aggregate Edge Weights

Aggregates a vector of edge weights using various methods. Compatible
with igraph's edge.attr.comb parameter.

## Usage

``` r
net_aggregate_weights(w, method = "sum", n_possible = NULL)
```

## Arguments

- w:

  Numeric vector of finite edge weights. `NA` and zero weights are
  excluded before aggregation, so every method (including `"density"`,
  `"min"`, `"max"`, `"prod"`, `"geomean"`) operates on the non-zero,
  non-`NA` subset.

- method:

  Single aggregation method: "sum", "mean", "median", "max", "min",
  "prod", "density", or "geomean". Because zeros are stripped first,
  `"density"` (`sum(w) / n_possible`) and `"mean"`
  (`sum(w) / number of non-zero edges`) return the *same* value whenever
  the block is fully dense – i.e. when the count of non-zero edges
  equals `n_possible`. They diverge only when zero/`NA` edges are
  present (then `"density"` divides by the larger `n_possible`, `"mean"`
  by the smaller non-zero count).

- n_possible:

  Optional single finite numeric number of possible edges for density
  calculation. When omitted, `"density"` falls back to
  `sum(w) / length(w)` on the non-zero subset (equivalent to `"mean"`);
  supply `n_possible` (e.g. the block size `n_i * n_j`) for a true edge
  density.

## Value

Single aggregated value

## Examples

``` r
w <- c(0.5, 0.8, 0.3, 0.9)
net_aggregate_weights(w, "sum")   # 2.5
#> [1] 2.5
net_aggregate_weights(w, "mean")  # 0.625
#> [1] 0.625
net_aggregate_weights(w, "max")   # 0.9
#> [1] 0.9
net_aggregate_weights(w, "density", n_possible = 9)  # 2.5 / 9
#> [1] 0.2777778
```
