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
  excluded before aggregation.

- method:

  Single aggregation method: "sum", "mean", "median", "max", "min",
  "prod", "density", or "geomean".

- n_possible:

  Optional single finite numeric number of possible edges for density
  calculation.

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
