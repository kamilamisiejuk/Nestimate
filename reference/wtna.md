# Window-based Transition Network Analysis

Computes networks from one-hot (binary indicator) data using temporal
windowing. Supports transition (directed), co-occurrence (undirected),
or both network types.

## Usage

``` r
wtna(
  data,
  method = c("transition", "cooccurrence", "both"),
  type = c("frequency", "relative"),
  codes = NULL,
  window_size = 3L,
  mode = c("non-overlapping", "overlapping"),
  actor = NULL
)
```

## Arguments

- data:

  Data frame with one-hot encoded columns (0/1 binary).

- method:

  Character. Network type: `"transition"` (directed), `"cooccurrence"`
  (undirected), or `"both"` (returns list of two networks). Default:
  `"transition"`.

- type:

  Character. Output type: `"frequency"` (raw counts) or `"relative"`
  (row-normalized probabilities). Default: `"frequency"`. Note that
  `type = "relative"` applied to `method = "cooccurrence"` produces an
  asymmetric matrix (conditional co-occurrence given row state), not a
  symmetric undirected weight matrix — use `type = "frequency"` if
  symmetric co-occurrence counts are required.

- codes:

  Character vector or NULL. Names of the one-hot columns to use. If
  NULL, auto-detects binary columns. Default: NULL.

- window_size:

  Integer (\>= 1). Number of consecutive rows to aggregate per window.
  Default: 3 (windowed pairwise between-window counting). Set
  `window_size = 1` for ordinary consecutive (t -\> t+1) transitions
  with no windowing.

- mode:

  Character. Window mode: `"non-overlapping"` (fixed, separate windows)
  or `"overlapping"` (rolling, step = 1). Default: `"non-overlapping"`.

- actor:

  Character or NULL. Name of the actor/ID column for per-group
  computation. If NULL, treats all rows as one group. Default: NULL.

## Value

For `method = "transition"` or `"cooccurrence"`: a `netobject` (see
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).

For `method = "both"`: a `wtna_mixed` object with elements `$transition`
and `$cooccurrence`, each a `netobject`.

## Details

**Transitions**: Uses `crossprod(X[-n,], X[-1,])` to count how often
state i is active at time t AND state j at time t+1.

**Co-occurrence**: Uses `crossprod(X)` to count states that are
simultaneously active in the same row.

**Windowing**: For `window_size > 1`, rows are aggregated into windows
before computing networks. Non-overlapping windows are fixed, separate
blocks; overlapping windows roll forward one row at a time. Within each
window, any active indicator (1) in any row makes that state active for
the window.

**Per-actor**: When `actor` is specified, networks are computed per
group and summed.

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md),
[`prepare_onehot`](https://saqr.me/Nestimate/reference/prepare_onehot.md)

## Examples

``` r
oh <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,0,0), nrow = 4, byrow = TRUE,
             dimnames = list(NULL, c("A","B","C")))
w <- wtna(oh)

# \donttest{
# Simple one-hot data
df <- data.frame(
  A = c(1, 0, 1, 0, 1),
  B = c(0, 1, 0, 1, 0),
  C = c(0, 0, 1, 0, 0)
)

# Transition network
net <- wtna(df)
print(net)
#> Network (method: wtna_transition) [directed]
#>   Weights: [1.000, 2.000]  |  mean: 1.333
#> 
#>   Weight matrix:
#>     A B C
#>   A 2 2 0
#>   B 1 1 0
#>   C 1 1 0 
#> 
#>   Initial probabilities:
#>   A             1.000  ████████████████████████████████████████
#>   B             0.000  
#>   C             0.000  

# Both networks
nets <- wtna(df, method = "both")
print(nets$transition)
#> Network (method: wtna_transition) [directed]
#>   Weights: [1.000, 2.000]  |  mean: 1.333
#> 
#>   Weight matrix:
#>     A B C
#>   A 2 2 0
#>   B 1 1 0
#>   C 1 1 0 
#> 
#>   Initial probabilities:
#>   A             1.000  ████████████████████████████████████████
#>   B             0.000  
#>   C             0.000  
print(nets$cooccurrence)
#> Network (method: wtna_cooccurrence) [undirected]
#>   Weights: [1.000, 5.000]  |  mean: 2.333
#> 
#>   Weight matrix:
#>     A B C
#>   A 5 3 2
#>   B 3 2 1
#>   C 2 1 1 

# With windowing
net <- wtna(df, window_size = 2, mode = "non-overlapping")
# }
```
