# Print Method for net_cluster_diagnostics

Prints a uniform header, family-specific quality / IC line, and a
per-cluster table. Layout matches
[`print.net_clustering`](https://saqr.me/Nestimate/reference/print.net_clustering.md)
and
[`print.net_mmm`](https://saqr.me/Nestimate/reference/print.net_mmm.md).

## Usage

``` r
# S3 method for class 'net_cluster_diagnostics'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `net_cluster_diagnostics` object.

- digits:

  Integer. Decimal places for floating-point statistics. Default `3L`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.
