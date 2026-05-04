# Per-Class State Distribution as a Tidy Data Frame

Returns a tidy `data.frame(group, state, count, proportion)` with one
row per (group, state) cell. Companion to
[`state_frequencies`](https://mohsaqr.github.io/Nestimate/reference/state_frequencies.md)
(which counts unique states in raw sequence input);
`state_distribution()` pulls the same shape of frame from a fitted
Nestimate object so analyses don't have to reach for the underlying
`$data` slot directly.

## Usage

``` r
state_distribution(x, ...)

# S3 method for class 'netobject'
state_distribution(x, ...)

# S3 method for class 'htna'
state_distribution(x, ...)

# S3 method for class 'mcml'
state_distribution(x, include_macro = FALSE, ...)

# S3 method for class 'netobject_group'
state_distribution(x, ...)

# Default S3 method
state_distribution(x, ...)
```

## Arguments

- x:

  A `netobject`, `netobject_group`, `mcml`, or `htna` object.

- ...:

  Currently unused.

- include_macro:

  For `mcml`: when `TRUE`, prepend a `group = "macro"` block aggregating
  across clusters. Ignored for the other classes.

## Value

A `data.frame` with columns `group` (character), `state` (character),
`count` (integer), and `proportion` (numeric, within-group share).

## Details

Used internally by
[`plot_state_frequencies`](https://mohsaqr.github.io/Nestimate/reference/plot_state_frequencies.md)
as the data layer behind every chart, and surfaced as the `$table` slot
of the returned `state_freq` object.

## Examples

``` r
if (FALSE) { # \dontrun{
  data(ai_long)
  net <- build_network(ai_long, method = "frequency",
                       id_col = "session_id",
                       time_col = "order_in_session", action = "code")
  state_distribution(net)
} # }
```
