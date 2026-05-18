# Convert Sequence Data to Different Formats

Convert wide or long sequence data into frequency counts, one-hot
encoding, edge lists, or follows format.

## Usage

``` r
convert_sequence_format(
  data,
  seq_cols = NULL,
  id_col = NULL,
  action = NULL,
  time = NULL,
  format = c("frequency", "onehot", "edgelist", "follows")
)
```

## Arguments

- data:

  Data frame containing sequence data.

- seq_cols:

  Character vector. Names of columns containing sequential states (for
  wide format input). If NULL, all columns except `id_col` are used.
  Default: NULL.

- id_col:

  Character vector. Name(s) of the ID column(s). For long format,
  required. For wide format, optional: if NULL, the first column is used
  as the id only when it is a genuine identifier (its values are
  disjoint from the states in the remaining columns); for canonical wide
  sequence data with no id column (e.g. `V1..Vn` / `T1..Tn`), a
  row-index id is synthesized and every column is treated as a sequence
  column. Default: NULL.

- action:

  Character or NULL. Name of the column containing actions/states (for
  long format input). If provided, data is treated as long format.
  Default: NULL.

- time:

  Character or NULL. Name of the time column for ordering actions within
  sequences (for long format). Default: NULL.

- format:

  Character. Output format:

  "frequency"

  :   Count of each action per sequence (wide, one column per state).

  "onehot"

  :   Binary presence/absence of each action per sequence.

  "edgelist"

  :   Consecutive transition pairs (from, to) per sequence.

  "follows"

  :   Each action paired with the action that preceded it.

## Value

A data frame in the requested format:

- frequency:

  ID columns + one integer column per state with counts.

- onehot:

  ID columns + one binary column per state (0/1).

- edgelist:

  ID columns + `from` and `to` columns.

- follows:

  ID columns + `act` and `follows` columns.

## See also

[`frequencies`](https://saqr.me/Nestimate/reference/frequencies.md) for
building transition frequency matrices.

## Examples

``` r
# Wide format input
seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
convert_sequence_format(seqs, format = "frequency")
#>   rid .rid A B C
#> 1   1    1 2 1 0
#> 2   2    2 1 1 1
#> 3   3    3 1 1 1
convert_sequence_format(seqs, format = "edgelist")
#>   .rid from to
#> 1    1    A  B
#> 2    1    B  A
#> 3    2    B  A
#> 4    2    A  C
#> 5    3    A  C
#> 6    3    C  B
```
