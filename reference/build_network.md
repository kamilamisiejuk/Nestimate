# Build a Network

Universal network estimation function that supports both transition
networks (relative, frequency, co-occurrence) and association networks
(correlation, partial correlation, graphical lasso). Uses the global
estimator registry, so custom estimators can also be used.

## Usage

``` r
build_network(
  data,
  method,
  actor = NULL,
  action = NULL,
  time = NULL,
  session = NULL,
  order = NULL,
  codes = NULL,
  group = NULL,
  format = "auto",
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
  scaling = NULL,
  threshold = 0,
  level = NULL,
  time_threshold = 900,
  predictability = TRUE,
  state_cols = NULL,
  metadata_cols = NULL,
  params = list(),
  labels = NULL,
  ...
)
```

## Arguments

- data:

  Data frame (sequences or per-observation frequencies) or a square
  symmetric matrix (correlation or covariance).

- method:

  Character. Required. Name of a registered estimator. Built-in methods:
  `"relative"`, `"frequency"`, `"co_occurrence"`, `"cor"`, `"pcor"`,
  `"glasso"`, `"ising"`, `"mgm"`, `"attention"`, `"wtna"`,
  `"wtna_cooccurrence"`. Aliases: `"tna"` and `"transition"` map to
  `"relative"`; `"ftna"` and `"counts"` map to `"frequency"`; `"cna"`
  and `"wcna"` map to `"co_occurrence"`; `"corr"` and `"correlation"`
  map to `"cor"`; `"partial"` maps to `"pcor"`; `"ebicglasso"` and
  `"regularized"` map to `"glasso"`; `"isingfit"` maps to `"ising"`;
  `"atna"` maps to `"attention"`; `"mixed"` and `"mixed_graphical"` map
  to `"mgm"`; `"wtna_transition"` maps to `"wtna"`.

- actor:

  Character. Name of the actor/person ID column for sequence grouping.
  Default: `NULL`.

- action:

  Character. Name of the action/state column (long format). Default:
  `NULL`.

- time:

  Character. Name of the time column (long format). Default: `NULL`.

- session:

  Character. Name of the session column. Default: `NULL`.

- order:

  Character. Name of the ordering column. Default: `NULL`.

- codes:

  Character vector. Column names of one-hot encoded states (for onehot
  format). Default: `NULL`.

- group:

  Character. Name of a grouping column for per-group networks. Returns a
  `netobject_group` (named list of netobjects). Default: `NULL`.

- format:

  Character. Input format: `"auto"`, `"wide"`, `"long"`, or `"onehot"`.
  Default: `"auto"`.

- window_size:

  Integer. Window size for one-hot windowing. Default: `3L`.

- mode:

  Character. Windowing mode for one-hot input only: `"non-overlapping"`
  or `"overlapping"`. Has no effect on wide or long sequence data (only
  the one-hot/wtna path reads it). Default: `"non-overlapping"`.

- scaling:

  Character vector or NULL. Post-estimation scaling to apply (in order).
  Options: `"minmax"`, `"max"`, `"rank"`, `"normalize"`. Can combine:
  `c("rank", "minmax")`. Default: `NULL` (no scaling).

- threshold:

  Numeric. Absolute values below this are set to zero in the result
  matrix. Default: 0 (no thresholding).

- level:

  Character or NULL. Multilevel decomposition for association methods.
  One of `NULL`, `"between"`, `"within"`, `"both"`. Requires `id_col`.
  Default: `NULL`.

- time_threshold:

  Numeric. Maximum time gap (seconds) for long format session splitting.
  Default: `900`.

- predictability:

  Logical. If `TRUE` (default), compute and store node predictability
  (R-squared) for undirected association methods (glasso, pcor, cor).
  Stored in `$predictability` and auto-displayed as donuts by
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

- state_cols:

  Character vector or `NULL`. Explicit names of columns to classify as
  state columns in the returned netobject's `$data` slot. When provided,
  all other columns of the cleaned input go to `$metadata`.
  Auto-detection (values-in-nodes heuristic) is bypassed. Use this when
  a metadata column happens to contain values that overlap with node
  names (e.g. condition labels `"A","B","C"` and nodes `"A","B","C"`)
  and auto-detection would misclassify it. Default: `NULL`
  (auto-detect).

- metadata_cols:

  Character vector or `NULL`. Explicit names of columns to force into
  the `$metadata` slot. The remaining columns are auto-detected as state
  via the values-in-nodes rule. Cannot overlap with `state_cols`.
  Default: `NULL`.

- params:

  Named list. Method-specific parameters passed to the estimator
  function (e.g. `list(gamma = 0.5)` for glasso, or
  `list(format = "wide")` for transition methods). This is the key
  composability feature: downstream functions like bootstrap or grid
  search can store and replay the full params list without knowing
  method internals. Transition estimators accept tna-style sequence
  options such as `weighted`, `begin_state`, `end_state`, and `concat`.
  Column-like entries in `params` (`action`, `id`, `id_col`, `time`,
  `session`, `order`, `codes`, and `group`) are resolved before format
  detection and must name existing columns. If the same column role is
  supplied both directly and through `params`, the names must agree.

- labels:

  Optional name -\> label remap applied after construction. Accepts a
  2-column data.frame `(name, label)`, a named character vector
  `c(name = "label")`, or a named list. Rewrites `$nodes$label` and
  `dimnames(weights)`. Unmapped names pass through unchanged.

- ...:

  Additional arguments passed to the estimator function.

## Value

An object of class `c("netobject", "cograph_network")` containing:

- data:

  The input data used for estimation, as a data frame.

- weights:

  The estimated network weight matrix.

- nodes:

  Data frame with columns `id`, `label`, `name`, `x`, `y`. Node labels
  are in `$nodes$label`.

- edges:

  Data frame of non-zero edges with integer `from`/`to` (node IDs) and
  numeric `weight`.

- directed:

  Logical. Whether the network is directed.

- method:

  The resolved method name.

- params:

  The params list used (for reproducibility).

- scaling:

  The scaling applied (or NULL).

- threshold:

  The threshold applied.

- n_nodes:

  Number of nodes.

- n_edges:

  Number of non-zero edges.

- level:

  Decomposition level used (or NULL).

- meta:

  List with `source`, `layout`, and `tna` metadata (cograph-compatible).

- node_groups:

  Node groupings data frame, or NULL.

- predictability:

  Named numeric vector of R-squared predictability values per node (for
  undirected association methods when `predictability = TRUE`). NULL for
  directed methods.

Method-specific extras (e.g. `precision_matrix`, `cor_matrix`,
`frequency_matrix`, `lambda_selected`, etc.) are preserved from the
estimator output.

When `level = "both"`, returns an object of class `"netobject_ml"` with
`$between` and `$within` sub-networks and a `$method` field.

## Details

The function works as follows:

1.  Resolves method aliases to canonical names.

2.  Validates explicit column arguments before any format guessing.

3.  Retrieves the estimator function from the global registry.

4.  For association methods with `level` specified, decomposes the data
    (between-person means or within-person centering).

5.  Calls the estimator: `do.call(fn, c(list(data = data), params))`.

6.  Applies scaling and thresholding to the result matrix.

7.  Extracts edges and constructs the `netobject`.

For long-format transition data, supplying `action` without `actor` is
allowed and treats all rows as one sequence in row/time order. The
function warns because a one-sequence transition network is not
recommended and cannot be validated by bootstrap or other confirmatory
tests.

## See also

[`register_estimator`](https://saqr.me/Nestimate/reference/register_estimator.md),
[`list_estimators`](https://saqr.me/Nestimate/reference/list_estimators.md),
[`bootstrap_network`](https://saqr.me/Nestimate/reference/bootstrap_network.md)

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
net <- build_network(seqs, method = "relative")
net
#> Transition Network (relative probabilities) [directed]
#>   Weights: [1.000, 1.000]  |  mean: 1.000
#> 
#>   Weight matrix:
#>     A B C
#>   A 0 1 0
#>   B 0 0 1
#>   C 1 0 0 
#> 
#>   Initial probabilities:
#>   A             0.500  ████████████████████████████████████████
#>   B             0.250  ████████████████████
#>   C             0.250  ████████████████████
# \donttest{
# Transition network (relative probabilities)
seqs <- data.frame(
  V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
  V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
print(net)
#> Transition Network (relative probabilities) [directed]
#>   Weights: [0.111, 0.556]  |  mean: 0.250
#> 
#>   Weight matrix:
#>         A     B     C     D
#>   A 0.423 0.192 0.192 0.192
#>   B 0.304 0.304 0.217 0.174
#>   C 0.217 0.304 0.304 0.174
#>   D 0.222 0.111 0.556 0.111 
#> 
#>   Initial probabilities:
#>   A             0.333  ████████████████████████████████████████
#>   B             0.300  ████████████████████████████████████
#>   D             0.200  ████████████████████████
#>   C             0.167  ████████████████████

# Association network (glasso)
freq_data <- convert_sequence_format(seqs, format = "frequency")
net_glasso <- build_network(freq_data, method = "glasso",
                             params = list(gamma = 0.5, nlambda = 50))

# With scaling
net_scaled <- build_network(seqs, method = "relative",
                             scaling = c("rank", "minmax"))
# }
```
