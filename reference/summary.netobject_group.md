# Network metrics for a netobject_group

Returns one summary per constituent network. With `combined = TRUE`
(default) the per-group tables are joined into a single wide
`data.frame` with one column per group; with `combined = FALSE` returns
a named list.

## Usage

``` r
# S3 method for class 'netobject_group'
summary(object, combined = TRUE, ...)
```

## Arguments

- object:

  A `netobject_group`.

- combined:

  Logical. Combine into one wide data.frame? Default `TRUE`.

- ...:

  Ignored.

## Value

Either a `data.frame` (one column per group) or a named list of
`summary.netobject` objects, of class
`c("summary.netobject_group", ...)`.
