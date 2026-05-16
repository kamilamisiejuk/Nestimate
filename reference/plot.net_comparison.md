# Plot a network comparison

Visualises a `net_comparison` object. Currently supports the edge-weight
scatterplot (default), with the diagonal reference (perfect agreement)
and the OLS regression line annotated by Pearson, Spearman, and Kendall
correlations.

## Usage

``` r
# S3 method for class 'net_comparison'
plot(
  x,
  type = c("scatter", "heatmap", "diff_hist", "weight_dist", "all"),
  combined = TRUE,
  ...
)
```

## Arguments

- x:

  A `net_comparison` object from
  [`compare_model()`](https://saqr.me/Nestimate/reference/compare_model.md).

- type:

  Character. One of `"scatter"` (default — edge-weight scatter with OLS
  fit and correlation overlay), `"heatmap"` (n by n grid of x - y
  differences using the diverging palette), `"diff_hist"` (histogram of
  \|x - y\| absolute differences with rug + density), `"weight_dist"`
  (overlaid distributions of \|x\| and \|y\| edge weights), or `"all"`
  (2 by 2 grid of all four panels; requires the gridExtra package).

- combined:

  When `type = "all"` and `combined = TRUE` (default), the four panels
  are stitched into a 2x2 gtable. When `FALSE`, returns a named list of
  the four ggplots so each can be printed, saved, or re-laid-out
  independently. Ignored for other `type` values.

- ...:

  Ignored.

## Value

A `ggplot` object; for `type = "all"` with `combined = TRUE` a `gtable`
arranged 2 by 2; for `type = "all"` with `combined = FALSE` a named list
of four ggplots.
