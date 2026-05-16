# Compare two networks within a netobject_group

Selects two members of a `netobject_group` (by index or name) and
dispatches to
[`compare_model.netobject()`](https://saqr.me/Nestimate/reference/compare_model.md).

## Usage

``` r
# S3 method for class 'netobject_group'
compare_model(
  x,
  i = 1L,
  j = 2L,
  scaling = "none",
  measures = character(0),
  network = TRUE,
  ...
)
```

## Arguments

- x:

  A `netobject_group`.

- i:

  Integer or character. Index or name of the first network. Default
  `1L`.

- j:

  Integer or character. Index or name of the second network. Default
  `2L`.

- scaling:

  See
  [`compare_model()`](https://saqr.me/Nestimate/reference/compare_model.md).

- measures:

  See
  [`compare_model()`](https://saqr.me/Nestimate/reference/compare_model.md).

- network:

  See
  [`compare_model()`](https://saqr.me/Nestimate/reference/compare_model.md).

- ...:

  Passed to
  [`compare_model.netobject()`](https://saqr.me/Nestimate/reference/compare_model.md).

## Value

A `net_comparison` object.
