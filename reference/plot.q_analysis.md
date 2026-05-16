# Plot Q-Analysis

Two panels: Q-vector (components at each connectivity level) and
structure vector (max simplex dimension per node).

## Usage

``` r
# S3 method for class 'q_analysis'
plot(x, combined = TRUE, ...)
```

## Arguments

- x:

  A `q_analysis` object.

- combined:

  When `TRUE` (default), the two panels are stitched side-by-side via
  [`gridExtra::arrangeGrob`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html).
  When `FALSE`, returns a named list (`q_vector`, `structure_vector`) of
  ggplots.

- ...:

  Ignored.

## Value

A grid grob (invisibly) when `combined = TRUE`; a named list of two
ggplots when `combined = FALSE`.

## Examples

``` r
# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
net <- build_network(seqs, method = "relative")
sc  <- build_simplicial(net, type = "clique")
qa  <- q_analysis(sc)
plot(qa)

# }
```
