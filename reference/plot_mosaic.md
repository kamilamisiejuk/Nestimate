# Draw a Marimekko / Mosaic Plot from a Tidy Data Frame

Low-level rectangle-coordinate builder for marimekko (mosaic) plots.
Column widths are proportional to the per-column total of `weight`;
within each column, segments stack to height 1 with sub-heights
proportional to each row's share of that column's total.

## Usage

``` r
plot_mosaic(
  data,
  x,
  y,
  weight,
  fill = "y",
  colors = NULL,
  show_labels = TRUE,
  label_size = 3.5,
  x_label = NULL,
  y_label = NULL
)
```

## Arguments

- data:

  A data.frame in long form. Must contain the columns named in `x`, `y`,
  and `weight`.

- x:

  Column name for the X (column) variable.

- y:

  Column name for the Y (segment) variable.

- weight:

  Column name for the cell weight (e.g. count).

- fill:

  Either `"y"` (color by Y category, e.g. state – default) or the name
  of another column to map to fill (e.g. a residual column for diverging
  color).

- colors:

  Optional character vector of fill colors. When `fill = "y"`, length
  must be at least the number of distinct y levels. Defaults to recycled
  Okabe-Ito.

- show_labels:

  If `TRUE`, draw within-segment percentage labels.

- label_size:

  Numeric size for segment labels.

- x_label, y_label:

  Optional axis labels.

## Value

A `ggplot` object.

## Details

Used internally by
[`plot_state_frequencies`](https://mohsaqr.github.io/Nestimate/reference/plot_state_frequencies.md);
exposed so that other plot methods (e.g. permutation-residual
visualisations) can reuse the same geometry by supplying a different
fill column.

## Examples

``` r
df <- data.frame(
  group = rep(c("A", "B", "C"), each = 3),
  state = rep(c("s1", "s2", "s3"), 3),
  count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
)
plot_mosaic(df, x = "group", y = "state", weight = "count")
```
