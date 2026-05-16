# Import One-Hot Encoded Data into Sequence Format

Converts binary indicator (one-hot) data into the wide sequence format
expected by
[`build_network`](https://saqr.me/Nestimate/reference/build_network.md)
and [`tna::tna()`](http://sonsoles.me/tna/reference/build_model.md).
Each binary column represents a state; rows where the value is 1 are
marked with the column name. Supports optional windowed aggregation.

Simultaneous active states are preserved using the same window-span
representation as
[`tna::import_onehot()`](http://sonsoles.me/tna/reference/import_onehot.md):
each input row/window is expanded to one sequence slot per code and
transition counting occurs between windows, not between simultaneous
states inside the same row.

## Usage

``` r
prepare_onehot(
  data,
  cols,
  actor = NULL,
  session = NULL,
  interval = NULL,
  window_size = 3L,
  window_type = c("non-overlapping", "overlapping"),
  aggregate = FALSE
)
```

## Arguments

- data:

  Data frame with binary (0/1) indicator columns.

- cols:

  Character vector. Names of the one-hot columns to use.

- actor:

  Character or NULL. Name of the actor/ID column. If NULL, all rows are
  treated as a single sequence. Default: NULL.

- session:

  Character or NULL. Name of the session column for sub-grouping within
  actors. Default: NULL.

- interval:

  Integer or NULL. Number of rows per time point in the output. If NULL,
  all rows become a single time point group. Default: NULL.

- window_size:

  Integer. Number of consecutive rows to aggregate into each window.
  Default: 1 (no windowing).

- window_type:

  Character. `"non-overlapping"` (fixed, separate windows) or
  `"overlapping"` (rolling, step = 1). Default: `"non-overlapping"`.

- aggregate:

  Logical. If TRUE, aggregate within each window by taking the first
  non-NA indicator per column. Default: FALSE.

## Value

A data frame in wide format with columns named `W{window}_T{time}` where
each cell contains a state name or NA. Attributes `windowed`,
`window_size`, `window_span` are set on the result.

## See also

[`action_to_onehot`](https://saqr.me/Nestimate/reference/action_to_onehot.md)
for the reverse conversion.

## Examples

``` r
# Simple binary data
df <- data.frame(
  A = c(1, 0, 1, 0, 1),
  B = c(0, 1, 0, 1, 0),
  C = c(0, 0, 0, 0, 0)
)
seq_data <- prepare_onehot(df, cols = c("A", "B", "C"))

# With actor grouping
df$actor <- c(1, 1, 1, 2, 2)
seq_data <- prepare_onehot(df, cols = c("A", "B", "C"), actor = "actor")

# With windowing
seq_data <- prepare_onehot(df, cols = c("A", "B", "C"),
                          window_size = 2, window_type = "non-overlapping")
```
