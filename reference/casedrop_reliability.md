# Edge-weight Case-dropping Stability

Computes a **CS-coefficient for the edge-weight vector** of a network:
the maximum proportion of cases (rows of `x$data`) that can be dropped
while the flattened edge-weight vector of the re-estimated network still
correlates with the original above `threshold` in at least `certainty`
of iterations.

Plots the four model-level reliability metrics across drop proportions:
`correlation`, `mean_abs_dev`, `median_abs_dev`, `max_abs_dev`. Each
panel shows the per-iteration mean with a ribbon at mean +/- sd. The
`correlation` panel includes a dashed horizontal line at the user's
`threshold` (default 0.7).

Overlay of per-cluster correlation curves across drop proportions. One
colour per sub-network; ribbons show mean +/- sd across iterations.
Dashed horizontal line marks the stability threshold (default 0.7).

## Usage

``` r
casedrop_reliability(
  x,
  iter = 1000L,
  drop_prop = seq(0.1, 0.9, by = 0.1),
  threshold = 0.7,
  certainty = 0.95,
  method = c("spearman", "pearson", "kendall"),
  include_diag = FALSE,
  seed = NULL
)

# S3 method for class 'net_casedrop_reliability'
print(x, digits = 3, ...)

# S3 method for class 'net_casedrop_reliability'
summary(object, ...)

# S3 method for class 'net_casedrop_reliability_group'
print(x, ...)

# S3 method for class 'net_casedrop_reliability_group'
summary(object, drop_prop = 0.7, ...)

# S3 method for class 'summary.net_casedrop_reliability_group'
print(x, ...)

# S3 method for class 'net_casedrop_reliability'
plot(x, combined = TRUE, ...)

# S3 method for class 'net_casedrop_reliability_group'
plot(
  x,
  metric = c("correlation", "mean_abs_dev", "median_abs_dev", "max_abs_dev"),
  ...
)
```

## Arguments

- x:

  A `net_casedrop_reliability_group` object.

- iter:

  Integer. Iterations per drop proportion. Default `1000`.

- drop_prop:

  Drop proportion at which to report the four metrics (mean +/- sd per
  network). Default `0.7`.

- threshold:

  Numeric in `[0, 1]`. Minimum edge-vector correlation for an iteration
  to count as stable. Default `0.7`.

- certainty:

  Numeric in `[0, 1]`. Required fraction of iterations whose correlation
  must exceed `threshold` for a drop proportion to qualify. Default
  `0.95`.

- method:

  Correlation method: `"pearson"` (weight magnitudes), `"spearman"`
  (ranks, robust to scale), or `"kendall"`. Default `"spearman"` because
  edge weights often span several orders of magnitude and rank stability
  is the typical target.

- include_diag:

  Logical. Include diagonal (self-loop) edges in the edge vector.
  Default `FALSE`.

- seed:

  Optional integer for reproducibility.

- digits:

  Digits to display. Default `3`.

- ...:

  Additional arguments (ignored).

- object:

  A `net_casedrop_reliability_group`.

- combined:

  When `TRUE` (default), all four metrics are shown in one ggplot via
  `facet_wrap(~ metric)`. When `FALSE`, returns a named list of four
  single-panel ggplots, one per metric.

- metric:

  Which metric to plot. One of `"correlation"` (default),
  `"mean_abs_dev"`, `"median_abs_dev"`, `"max_abs_dev"`.

## Value

An object of class `net_casedrop_reliability` with:

- `cs`:

  Scalar CS-coefficient — the maximum drop proportion for which the
  edge-vector correlation remains \>= `threshold` in at least
  `certainty` of iterations. Zero if no proportion qualifies.

- `correlations`:

  `iter` x `length(drop_prop)` matrix of per- iteration correlations.

- `drop_prop`, `threshold`, `certainty`, `iter`, `method`:

  Inputs.

The input `x` invisibly.

A tidy data frame with columns `metric`, `drop_prop`, `mean`, `sd`
summarising edge-weight stability across case-dropping iterations.

A data frame with one row per network containing `cor`, `mean_abs_dev`,
`median_abs_dev`, `max_abs_dev` formatted as "mean +/- sd".

A `ggplot` object, or a named list of four ggplots when
`combined = FALSE`.

A `ggplot` object.

## Details

Complements
[`centrality_stability()`](https://saqr.me/Nestimate/reference/centrality_stability.md):
that function asks whether centrality *rankings* are stable; this one
asks whether the *edge-weight structure itself* is stable. For
MCML-derived networks where each row of `$data` is one transition, this
is case-dropping of **edges**.

For each `drop_prop` p and each iteration, a size `n_cases * (1 - p)`
subset of `$data` rows is selected **without replacement**, the network
is re-estimated using the same method/scaling/threshold as the input,
and the upper/lower-triangle (directed: all off-diagonal entries) of the
new weight matrix is flattened and correlated with the corresponding
vector of the original matrix. The correlation method defaults to
Spearman for robustness to the wide dynamic range of transition
probabilities.

Unlike bootstrap CIs, case-dropping does not estimate sampling variance
and so does not rely on the i.i.d. assumption. This makes it the
appropriate robustness check for **edgelist-derived** networks (where
rows of `$data` lack actor grouping), since dropping rows at random is a
well-posed operation regardless of within-actor correlation.

## References

Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
psychological networks and their accuracy: A tutorial paper. *Behavior
Research Methods* 50(1), 195-212.
[doi:10.3758/s13428-017-0862-1](https://doi.org/10.3758/s13428-017-0862-1)

## See also

[`centrality_stability()`](https://saqr.me/Nestimate/reference/centrality_stability.md),
[`bootstrap_network()`](https://saqr.me/Nestimate/reference/bootstrap_network.md).

## Examples

``` r
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE),
  V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
es  <- casedrop_reliability(net, iter = 50, drop_prop = c(0.1, 0.3, 0.5),
                      seed = 1)
print(es)
#> Edge-weight Case-dropping Stability
#>   Cases (rows of $data) : 30
#>   Edges assessed        : 12 (diagonal excluded)
#>   Iterations / prop     : 50
#>   Correlation method    : spearman
#>   CS-coefficient (r)    : 0.10  (threshold=0.70, certainty=0.95)
#> 
#> Model-level reliability across iterations (mean +/- sd per drop):
#>   drop_prop      p=0.1        p=0.3        p=0.5      
#>   mean|diff|      0.027+- 0.006   0.060+- 0.013   0.094+- 0.022
#>   MAD             0.021+- 0.008   0.051+- 0.016   0.078+- 0.018
#>   cor             0.921+- 0.035   0.787+- 0.109   0.674+- 0.190
#>   max|diff|       0.076+- 0.023   0.159+- 0.046   0.250+- 0.070
```
