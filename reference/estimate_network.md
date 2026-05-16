# Estimate a Network (Deprecated)

This function is deprecated. Use
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
instead.

## Usage

``` r
estimate_network(
  data,
  method = "relative",
  params = list(),
  scaling = NULL,
  threshold = 0,
  level = NULL,
  ...
)
```

## Arguments

- data:

  Data frame (sequences or per-observation frequencies) or a square
  symmetric matrix (correlation or covariance).

- method:

  Character. Defaults to `"relative"` for backward compatibility.

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

- ...:

  Additional arguments passed to
  [`build_network`](https://saqr.me/Nestimate/reference/build_network.md).

## Value

A `netobject` (see
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)).

## See also

[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)

## Examples

``` r
data <- data.frame(A = c("x","y","z","x"), B = c("y","x","z","y"))
net <- estimate_network(data, method = "relative")
#> Warning: 'estimate_network' is deprecated.
#> Use 'build_network' instead.
#> See help("Deprecated")
```
