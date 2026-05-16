# Rename the models of a `netobject_group`

Replaces the names of the constituent networks in a `netobject_group`
(or any object inheriting from it). Useful when
[`build_network()`](https://saqr.me/Nestimate/reference/build_network.md)
produced generic labels (e.g. `"Cluster 1"`, `"Cluster 2"`) and you want
to substitute meaningful ones (e.g. `"High engagement"`,
`"Low engagement"`).

## Usage

``` r
rename_models(x, new_names)

# S3 method for class 'netobject_group'
rename_models(x, new_names)

# Default S3 method
rename_models(x, new_names)
```

## Arguments

- x:

  A `netobject_group` (or any object inheriting from it, such as
  `net_mlvar`).

- new_names:

  A character vector of new names. Must have the same length as `x`,
  contain no `NA` or empty strings, and be unique.

## Value

A `netobject_group` of the same class with renamed members.

## Examples

``` r
if (FALSE) { # \dontrun{
  d   <- tna::group_regulation
  grp <- build_network(d, method = "tna",
                       group = sample(c("a", "b"), nrow(d), TRUE))
  grp <- rename_models(grp, c("High", "Low"))
  names(grp)
} # }
```
