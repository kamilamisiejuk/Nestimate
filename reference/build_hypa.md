# Detect Path Anomalies via HYPA

Constructs a k-th order De Bruijn graph from sequential trajectory data
and uses a hypergeometric null model to detect paths with anomalous
frequencies. Paths occurring more or less often than expected under the
null model are flagged as over- or under-represented.

## Usage

``` r
build_hypa(data, k = 3L, alpha = 0.05, min_count = 5L, p_adjust = "BH")
```

## Arguments

- data:

  A data.frame (rows = trajectories), list of character vectors, `tna`
  object, or `netobject` with sequence data. For `tna`/`netobject`,
  numeric state IDs are automatically converted to label names.

- k:

  Integer. Order of the De Bruijn graph (default 2). Detects anomalies
  in paths of length k.

- alpha:

  Numeric. Significance threshold for anomaly classification (default
  0.05). Paths with HYPA score \< alpha are under-represented; paths
  with score \> 1-alpha are over-represented.

- min_count:

  Integer. Minimum observed count for a path to be classified as
  anomalous (default 2). Paths with fewer observations are always
  classified as `"normal"` regardless of their HYPA score, since single
  occurrences are unreliable.

- p_adjust:

  Character. Method for multiple testing correction of p-values. Default
  `"BH"` (Benjamini-Hochberg FDR control). Accepts any method from
  [`p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html) or
  `"none"` to skip correction. Under- and over-representation p-values
  are adjusted separately (two-sided testing).

## Value

An object of class `net_hypa` with components:

- scores:

  Data frame with path, from, to, observed, expected, ratio, p_value,
  p_under, p_over, p_adjusted_under, p_adjusted_over, anomaly columns.
  The `path` column shows the full state sequence (e.g., "A -\> B -\>
  C"); `from` is the context (conditioning states); `to` is the next
  state; `ratio` is observed / expected; `p_value` is retained as an
  alias for `p_under`, the raw lower-tail hypergeometric CDF value;
  `p_over` is the inclusive upper-tail probability `P(X >= observed)`;
  `p_adjusted_under` and `p_adjusted_over` are the corrected p-values
  for under- and over-representation tests respectively.

- adjacency:

  Weighted adjacency matrix of the De Bruijn graph.

- xi:

  Fitted propensity matrix.

- k:

  Order of the De Bruijn graph.

- alpha:

  Significance threshold used.

- p_adjust:

  Multiple testing correction method used.

- n_anomalous:

  Number of anomalous paths detected.

- n_over:

  Number of over-represented paths.

- n_under:

  Number of under-represented paths.

- n_edges:

  Total number of edges.

- nodes:

  Node names in the De Bruijn graph.

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
& Schweitzer, F. (2020). HYPA: Efficient Detection of Path Anomalies in
Time Series Data on Networks. *SDM 2020*, 460–468.

## Examples

``` r
seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
hyp <- build_hypa(seqs, k = 2)
#> Warning: 'k' is deprecated; use 'order' instead.

# \donttest{
trajs <- list(c("A","B","C"), c("A","B","C"), c("A","B","C"),
              c("A","B","D"), c("C","B","D"), c("C","B","A"))
h <- build_hypa(trajs, k = 2)
#> Warning: 'k' is deprecated; use 'order' instead.
print(h)
#> HYPA: Path Anomaly Detection
#>   Order(s):     2
#>   Edges:        4
#>   Anomalous:    0 (alpha=0.05, p_adjust=BH)
#>     Over-repr:  0
#>     Under-repr: 0
# }
```
