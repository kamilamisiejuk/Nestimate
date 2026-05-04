# Print, Plot, and Convert a state_freq Object

[`plot_state_frequencies()`](https://mohsaqr.github.io/Nestimate/reference/plot_state_frequencies.md)
returns a `state_freq` object holding both the rendered chart and the
tidy frequency table. [`print()`](https://rdrr.io/r/base/print.html)
shows the table in the console,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) renders the
chart, and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the tidy table for downstream piping.

## Usage

``` r
# S3 method for class 'state_freq'
print(x, digits = 1, max_states = 20L, ...)

# S3 method for class 'state_freq'
plot(x, ...)

# S3 method for class 'state_freq'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `state_freq` object.

- digits:

  Number of decimal places for proportion / share columns.

- max_states:

  Cap on rows shown per group in the per-state table. The full table
  remains available via `x$table`.

- ...:

  Unused.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `invisible(x)`;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns
`invisible(NULL)` after drawing;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
`x$table`.
