# Mosaic Plot of a Network's Transition or Co-occurrence Counts

Draws a Hartigan-Friendly mosaic (marimekko geometry, chi-square
standardized-residual fill) for an integer-weighted network. Equivalent
in algorithm and appearance to
[`tna::plot_mosaic()`](http://sonsoles.me/tna/reference/plot_mosaic.md);
named differently to avoid an export clash when both packages are
attached.

## Usage

``` r
mosaic_plot(x, ...)

# Default S3 method
mosaic_plot(x, ...)

# S3 method for class 'netobject'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = 90,
  left_angle = 0,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ...
)

# S3 method for class 'htna'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = 90,
  left_angle = 0,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ...
)

# S3 method for class 'mcml'
mosaic_plot(
  x,
  level = c("macro", "clusters"),
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = 90,
  left_angle = 0,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ncol = 2L,
  ...
)

# S3 method for class 'netobject_group'
mosaic_plot(
  x,
  xlab = NULL,
  ylab = NULL,
  range = NULL,
  top_angle = 90,
  left_angle = 0,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ncol = 2L,
  ...
)

# S3 method for class 'table'
mosaic_plot(
  x,
  xlab = "Row",
  ylab = "Column",
  range = NULL,
  top_angle = NULL,
  left_angle = NULL,
  residuals = c("permutation", "asymptotic"),
  n_perm = 500L,
  seed = NULL,
  ...
)

# S3 method for class 'matrix'
mosaic_plot(x, ...)
```

## Arguments

- x:

  One of the four data-bearing Nestimate classes: `netobject` (single
  mosaic of `$weights`), `netobject_group` (one panel per group), `mcml`
  (between-cluster mosaic by default; per-cluster panels with
  `level = "within"`), or `htna` (single mosaic of `$weights`; htna
  inherits netobject so the geometry matches). Also accepts a
  contingency `table` or plain numeric `matrix` for ad-hoc plotting.

- ...:

  Ignored.

- xlab, ylab:

  Axis labels. `NULL` (default) draws no axis title. Pass any string to
  add one.

- range:

  Numeric of length 2 giving the lower and upper colour-scale limits for
  the standardized residual. `NULL` (default) auto-fits the limits to
  the symmetric range `c(-M, M)` where `M = max(|stdres|)`, so no signal
  is squished. Pass an explicit range (e.g. `c(-4, 4)` for tna-style
  display, `c(-6, 6)` for moderate clipping) to clamp the colour scale.

- top_angle, left_angle:

  Rotation in degrees for the top (x) and left (y) tick labels. `NULL`
  (default) uses the auto rule `90 if n_levels > 3 else 0` on each axis.
  Pass any numeric to override (e.g. `top_angle = 45, left_angle = 0`).

- residuals:

  One of `"permutation"` (default) or `"asymptotic"`. `"permutation"`
  computes empirical-null z-scores by shuffling one variable's labels
  against the other for `n_perm` draws and reporting
  `(O - mean_perm) / sd_perm` per cell. Robust on sparse tables.
  `"asymptotic"` returns `stats::chisq.test()$stdres` (the closed-form
  `(O - E) / sqrt(E*(1 - p_row)*(1 - p_col))` that vcd and tna use).

- n_perm:

  Number of permutations when `residuals = "permutation"`. Default 500;
  use `>= 1000` for stable tail estimates.

- seed:

  Optional integer seed for the permutation RNG. Use for reproducible
  plots; ignored when `residuals = "asymptotic"`.

- level:

  For `mcml` only. `"macro"` (default) draws a single mosaic of
  `x$macro$weights` (the cluster-by-cluster aggregate); `"clusters"`
  draws one mosaic per cluster from `x$clusters[[k]]$weights`, faceted
  into one combined ggplot.

- ncol:

  For `netobject_group`: number of columns in the small- multiples
  layout. Default 2.

- values:

  Logical. When `TRUE`, overlay each cell's standardized residual as a
  numeric label (one decimal). Text colour switches to white on
  saturated cells (\|stdres\| \> 1.5) and dark grey otherwise. Default
  `FALSE` – the colour bar legend already conveys the sign and
  magnitude.

## Value

A `ggplot` object (or a `gtable` from
[`gridExtra::arrangeGrob`](https://rdrr.io/pkg/gridExtra/man/arrangeGrob.html)
for `netobject_group` when gridExtra is available).

## Details

Column widths are proportional to row marginals of the weight matrix
(incoming totals when the matrix is transposed, as for transitions).
Within each column, segment heights are proportional to that row's
conditional distribution. Cell fill is the standardized residual from
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html), with a
diverging palette clipped to \\\pm 4\\. Mosaics only make sense for
integer-valued matrices, so this function rejects relative / glasso /
correlation networks.

## See also

[`plot_mosaic`](https://saqr.me/Nestimate/reference/plot_mosaic.md) for
the lower-level data.frame primitive.

## Examples

``` r
if (FALSE) { # \dontrun{
  net <- build_network(group_regulation, method = "frequency")
  mosaic_plot(net)
} # }
```
