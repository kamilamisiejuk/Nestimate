# Plot method for `chain_structure`.

Renders the hitting-probability matrix as a heatmap, with rows and
columns ordered by communicating class so the block structure is visible
at a glance. State labels along both axes are coloured by classification
(absorbing / recurrent / transient). The subtitle summarises the
chain-level properties (regular, reversible).

## Usage

``` r
# S3 method for class 'chain_structure'
plot(x, show_values = TRUE, digits = 2L, ...)
```

## Arguments

- x:

  A `chain_structure` object.

- show_values:

  Logical. If `TRUE` (default), prints the numeric probability inside
  each cell. Set `FALSE` for large state spaces (n \> 10) where labels
  overlap.

- digits:

  Integer. Decimal places for in-cell labels.

- ...:

  Ignored.

## Value

A `ggplot` object.

## Details

Cell colour encodes `P(ever reach j | start at i)`. The diagonal uses
the return-time convention (`P(return to j in >= 1 steps)`), matching
[`markovchain::hittingProbabilities`](https://rdrr.io/pkg/markovchain/man/hittingProbabilities.html).
A non-irreducible chain shows zero off-block entries – visual evidence
of one-way doors between behavioural phases. An absorbing chain shows a
column of 1's for the absorbing state.
