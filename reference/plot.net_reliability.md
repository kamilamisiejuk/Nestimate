# Plot Method for net_reliability

Density plots of split-half metrics faceted by metric type. Multi-model
comparisons show overlaid densities colored by model.

## Usage

``` r
# S3 method for class 'net_reliability'
plot(x, bins = 60L, combined = TRUE, ...)
```

## Arguments

- x:

  A `net_reliability` object.

- bins:

  Integer. Number of histogram bins per panel (default 60).

- combined:

  When `TRUE` (default), all four metrics are shown in one ggplot via
  `facet_wrap(~ metric)`. When `FALSE`, returns a named list of four
  single-panel ggplots, one per metric.

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object (invisibly), or a named list of four ggplots when
`combined = FALSE`.

## Examples

``` r
net <- build_network(data.frame(V1 = c("A","B","C","A"),
  V2 = c("B","C","A","B")), method = "relative")
rel <- network_reliability(net, iter = 10)
plot(rel)

# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
net <- build_network(seqs, method = "relative")
rel <- network_reliability(net, iter = 20, seed = 1)
plot(rel)

# }
```
