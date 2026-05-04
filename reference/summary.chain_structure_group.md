# Cross-group comparison of `chain_structure_group`.

Produces a single tidy data.frame with one row per (group, state)
combination, combining classification, persistence, sojourn, and – when
applicable – stationary or mean-absorption-time columns. Useful for
side-by-side reporting of
[`chain_structure()`](https://mohsaqr.github.io/Nestimate/reference/chain_structure.md)
across the members of a `netobject_group`.

## Usage

``` r
# S3 method for class 'chain_structure_group'
summary(object, ...)
```

## Arguments

- object:

  A `chain_structure_group`.

- ...:

  Ignored.

## Value

A `data.frame` with columns `group`, `state`, `classification`,
`period`, `persistence`, `return_probability`, `sojourn_steps`, plus
`stationary_probability` if all groups are irreducible and
`mean_absorption_time` if any group has absorbing states.
