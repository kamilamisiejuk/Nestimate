# Plot Method for cluster_choice

Six explicit chart types plus a smart `"auto"` default. The user picks
the shape; the function does not editorialise (no "best" annotation, no
interpretive subtitles, no inferred recommendation).

## Usage

``` r
# S3 method for class 'cluster_choice'
plot(
  x,
  type = c("auto", "lines", "bars", "heatmap", "tradeoff", "facet"),
  abbrev = FALSE,
  combined = TRUE,
  ...
)
```

## Arguments

- x:

  A `cluster_choice` object.

- type:

  Character. One of `"auto"` (default), `"lines"`, `"bars"`,
  `"heatmap"`, `"tradeoff"`, `"facet"`.

- abbrev:

  Logical. If `TRUE`, dissimilarity and method names shown on tick
  labels and point labels are shortened (e.g. `"hamming"` -\> `"ham"`,
  `"ward.D2"` -\> `"wD2"`). The legend shows the full canonical name.
  Default `FALSE`.

- combined:

  Only meaningful for `type = "facet"`. When `TRUE` (default), all
  methods are shown in one ggplot via `facet_wrap(~ method)`. When
  `FALSE`, returns a named list of single-panel ggplots, one per method.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `ggplot` object, invisibly; for `type = "facet"` with
`combined = FALSE`, a named list of ggplots.

## Details

Type cheat-sheet:

- `"auto"`:

  Default. Picks one of the others based on which axes were swept.
  k-only -\> `"lines"`; one categorical axis swept -\> `"bars"`; k plus
  one categorical -\> `"lines"`; k plus two categoricals -\> `"facet"`;
  both categoricals without k -\> `"heatmap"`.

- `"lines"`:

  Silhouette across k (and `mean_within_dist` when `k` is the only swept
  axis), one line per non-k axis when present.

- `"bars"`:

  Horizontal bar chart of silhouette per axis level. Bars sorted by
  silhouette.

- `"heatmap"`:

  Tiled silhouette across two categorical axes. Requires both
  `dissimilarity` and `method` swept.

- `"tradeoff"`:

  Scatter: silhouette (y) vs `size_ratio` (x). Works for any sweep;
  labels each point.

- `"facet"`:

  Lines vs k, colour by one categorical axis, facet by another. Requires
  `k` plus two categoricals.

Asking for a type the data can't support raises an error pointing at the
alternatives.
