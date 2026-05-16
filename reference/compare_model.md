# Compare two networks descriptively

Computes a battery of descriptive comparison metrics between two
networks or two weight matrices: weight deviations (mean / median / RMS
/ max absolute difference, relative mean absolute difference,
coefficient-of- variation ratio), four correlation measures (Pearson,
Spearman, Kendall, distance correlation), five dissimilarity measures
(Euclidean, Manhattan, Canberra, Bray-Curtis, Frobenius), five
similarity measures (Cosine, Jaccard, Dice, Overlap, RV), pattern
agreements, and side-by-side network metrics. Optionally adds centrality
differences and centrality correlations.

## Usage

``` r
compare_model(x, ...)

# S3 method for class 'netobject'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'cograph_network'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)

# S3 method for class 'matrix'
compare_model(
  x,
  y,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)
```

## Arguments

- x:

  A `netobject`, `cograph_network`, or numeric square matrix.

- ...:

  Ignored.

- y:

  A `netobject`, `cograph_network`, or numeric square matrix.

- scaling:

  Scaling applied to both weight matrices before comparison. One of:

  `"none"`

  :   Identity (default).

  `"minmax"`

  :   \\(w - \min) / (\max - \min)\\; maps to \\\[0, 1\]\\.

  `"max"`

  :   \\w / \max(\|w\|)\\; preserves sign.

  `"rank"`

  :   Min-max of average ranks; ordinal scaling.

  `"zscore"`

  :   \\(w - \bar w) / s_w\\; standard score.

  `"robust"`

  :   \\(w - \mathrm{med}(w)) / \mathrm{mad}(w)\\; Huber-style robust
      z-score, resists outliers.

  `"log"`

  :   \\\log(w)\\; requires \\w \> 0\\.

  `"log1p"`

  :   \\\log(1 + w)\\; admits \\w \ge 0\\.

  `"softmax"`

  :   Numerically stable softmax over the flattened vector.

  `"quantile"`

  :   Empirical CDF of the flattened vector.

  `"frobenius"`

  :   Divide the matrix by its Frobenius norm \\\\W\\\_F = \sqrt{\sum
      w\_{ij}^2}\\; matrix-level normalisation.

  `"row"`

  :   Row-stochastic normalisation (each row's absolute values sum to
      1). Only meaningful for non-negative matrices; rows summing to
      zero are left unchanged.

  Scalings that produce negative weights (`zscore`, `robust`) are
  compatible with `network = TRUE` because the side-by-side metrics use
  Nestimate's base-R Floyd-Warshall, which handles negative weights.

- measures:

  Character vector of centrality measures to compare. Empty by default
  (no centrality block).

- network:

  Logical. Include side-by-side network metrics from
  [`summary()`](https://rdrr.io/r/base/summary.html)? Default `TRUE`.

## Value

A `net_comparison` object: a named list with `matrices`,
`difference_matrix`, `edge_metrics`, `summary_metrics`, optionally
`network_metrics`, `centrality_differences`, `centrality_correlations`.

## Details

Mirrors [`tna::compare()`](http://sonsoles.me/tna/reference/compare.md)
numerically. Inputs are converted to weight matrices and scaled before
comparison; the choice of scaling determines how weights from different
estimators are placed on a common footing.
