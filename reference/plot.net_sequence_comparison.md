# Plot Method for net_sequence_comparison

Visualizes pattern-level standardized residuals across groups. Two
styles are available:

- `"pyramid"`:

  Back-to-back bars of pattern proportions, shaded by each side's
  standardized residual. Requires exactly 2 groups.

- `"heatmap"`:

  One tile per (pattern, group) cell, colored by standardized residual.
  Works for any number of groups.

Residuals are read directly from the `resid_<group>` columns in
`$patterns`, which are always populated regardless of the inference
method chosen in `sequence_compare`.

## Usage

``` r
# S3 method for class 'net_sequence_comparison'
plot(
  x,
  top_n = 10L,
  style = c("pyramid", "heatmap"),
  sort = c("statistic", "frequency"),
  alpha = 0.05,
  show_residuals = FALSE,
  ...
)
```

## Arguments

- x:

  A `net_sequence_comparison` object.

- top_n:

  Integer. Show top N patterns. Default: 10.

- style:

  Character. `"pyramid"` (default) or `"heatmap"`.

- sort:

  Character. `"statistic"` (default) ranks patterns by test statistic or
  residual magnitude. `"frequency"` ranks by total occurrence count
  across all groups.

- alpha:

  Numeric. Significance threshold for p-value display in the pyramid:
  patterns with `p_value < alpha` are starred and drawn in bold dark
  text, the rest stay plain grey. Default: 0.05.

- show_residuals:

  Logical. If `TRUE`, print the standardized residual value inside each
  pyramid bar. Default: `FALSE`. Ignored for the heatmap (which always
  shows residuals).

- ...:

  Additional arguments (ignored).

## Value

A `ggplot` object, invisibly.

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 60, TRUE),
  V2 = sample(LETTERS[1:4], 60, TRUE),
  V3 = sample(LETTERS[1:4], 60, TRUE),
  V4 = sample(LETTERS[1:4], 60, TRUE)
)
grp <- rep(c("X", "Y"), 30)
net <- build_network(seqs, method = "relative")
res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
```
