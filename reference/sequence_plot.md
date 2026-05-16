# Sequence Plot (heatmap, index, or distribution)

Single entry point for three categorical-sequence visualisations.

- `type = "heatmap"` (default): dense carpet, rows reordered by `sort` /
  dendrogram (single panel).

- `type = "index"`: same data layout, but rows separated by thin gaps
  (no dendrogram). Supports grouping via `group` or a `net_clustering`,
  plus a `ncol` x `nrow` facet grid.

- `type = "distribution"`: dispatches to
  [`distribution_plot`](https://saqr.me/Nestimate/reference/distribution_plot.md).

## Usage

``` r
sequence_plot(
  x,
  type = c("heatmap", "index", "distribution"),
  sort = c("lcs", "frequency", "start", "end", "hamming", "osa", "lv", "dl", "qgram",
    "cosine", "jaccard", "jw"),
  tree = NULL,
  group = NULL,
  scale = c("proportion", "count"),
  geom = c("area", "bar"),
  na = TRUE,
  row_gap = 0,
  dendrogram_width = 1.2,
  k = NULL,
  k_color = "white",
  k_line_width = 2.5,
  state_colors = NULL,
  na_color = "grey90",
  cell_border = NA,
  frame = FALSE,
  width = NULL,
  height = NULL,
  main = NULL,
  show_n = TRUE,
  time_label = "Time",
  xlab = NULL,
  y_label = NULL,
  ylab = NULL,
  tick = NULL,
  ncol = NULL,
  nrow = NULL,
  combined = TRUE,
  legend = NULL,
  legend_size = NULL,
  legend_title = NULL,
  legend_ncol = NULL,
  legend_border = NA,
  legend_bty = "n"
)
```

## Arguments

- x:

  Wide-format sequence data. Accepts:

  data.frame / matrix

  :   Rows = sequences, columns = time points.

  netobject

  :   Extracts `$data`.

  net_clustering

  :   From
      [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md).
      Uses `$data`, `$assignments` for grouping, and `$distance` for
      dendrogram.

  netobject_group

  :   From
      [`cluster_network`](https://saqr.me/Nestimate/reference/cluster_network.md)
      or
      [`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
      on a clustering. Extracts data and assignments from
      `attr(, "clustering")`.

  net_mmm

  :   From
      [`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md).
      Uses `$models[[1]]$data` and `$assignments`.

  tna

  :   From the tna package. Decodes integer-encoded sequences.

- type:

  One of `"heatmap"` (default), `"index"`, or `"distribution"`.

- sort:

  Row-ordering strategy for heatmap / within-panel for index. One of
  `"lcs"` (default), `"frequency"`, `"start"`, `"end"`, or any
  [`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)
  distance (`"hamming"`, `"osa"`, `"lv"`, `"dl"`, `"qgram"`, `"cosine"`,
  `"jaccard"`, `"jw"`).

- tree:

  Optional `hclust`/`dendrogram`/`agnes` object to supply row ordering
  (heatmap only; overrides `sort`).

- group:

  Optional grouping vector (length `nrow(x)`) producing one facet per
  group. Index/distribution only. Ignored for heatmap.

- scale, geom, na:

  Passed to
  [`distribution_plot`](https://saqr.me/Nestimate/reference/distribution_plot.md)
  when `type = "distribution"`.

- row_gap:

  Fraction of row height used as vertical gap between sequences in index
  plots. `0` (default) = dense like heatmap. Try `0.15` for visible
  separators at low row counts.

- dendrogram_width:

  Width ratio of the dendrogram panel (heatmap).

- k:

  Optional integer. When supplied in `type = "heatmap"`, cuts the
  dendrogram into `k` clusters and draws thin horizontal separators
  between them in the carpet. Ignored when there is no dendrogram (e.g.
  `sort = "start"`) or for other types.

- k_color:

  Colour for the cluster separator lines. Default `"white"`.

- k_line_width:

  Line width for the cluster separators. Default `2.5`.

- state_colors:

  Vector of colours, one per state.

- na_color:

  Colour for `NA` cells.

- cell_border:

  Cell border colour. `NA` = off.

- frame:

  If `TRUE` (default), draw a box around each panel. If `FALSE`, no
  box - axis ticks and labels still appear.

- width, height:

  Optional device dimensions in inches. When supplied, opens a new
  graphics device via
  [`grDevices::dev.new()`](https://rdrr.io/r/grDevices/dev.html). In
  knitr chunks use the `fig.width` / `fig.height` chunk options instead.

- main:

  Plot title.

- show_n:

  Append `"(n = N)"` to the title.

- time_label, xlab:

  X-axis label. `xlab` is an alias.

- y_label, ylab:

  Y-axis label (distribution only). `ylab` alias.

- tick:

  Show every Nth x-axis label. `NULL` = auto.

- ncol, nrow:

  Facet grid dimensions (index + distribution). Ignored when
  `combined = FALSE`.

- combined:

  Index and distribution types only. When `TRUE` (default), groups are
  arranged on one figure via
  [`graphics::layout()`](https://rdrr.io/r/graphics/layout.html). When
  `FALSE`, each group is drawn on its own page (one full-size figure per
  group, with its own legend). Single-group calls (`G == 1`) ignore this
  argument. Heatmap is always single-figure.

- legend:

  Legend position: `"bottom"`, `"right"`, or `"none"`. Default varies by
  type.

- legend_size:

  Legend text size. `NULL` (default) auto-scales from the device width
  so the legend looks proportional at 5 in vs 12 in figures (clamped to
  `[0.65, 1.2]`).

- legend_title:

  Optional legend title.

- legend_ncol:

  Number of legend columns.

- legend_border:

  Swatch border colour.

- legend_bty:

  `"n"` or `"o"`.

## Value

Invisibly, a list describing the plot (shape depends on `type`).

## See also

[`distribution_plot`](https://saqr.me/Nestimate/reference/distribution_plot.md),
[`build_clusters`](https://saqr.me/Nestimate/reference/build_clusters.md)

## Examples

``` r
# \donttest{
sequence_plot(trajectories)

sequence_plot(trajectories, type = "index")

sequence_plot(trajectories, type = "distribution")

# }
```
