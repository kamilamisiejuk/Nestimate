# Network metrics for a netobject

Computes node count, edge count, density, mean shortest-path distance,
mean and SD of in/out strength, mean and SD of in/out degree, in/out
degree centralization (Freeman), and reciprocity. Mirrors the metric set
returned by
[`tna::summary.tna()`](http://sonsoles.me/tna/reference/summary.tna.md)
so a Nestimate netobject and the equivalent tna model report numerically
identical descriptive metrics.

## Usage

``` r
# S3 method for class 'netobject'
summary(object, ...)
```

## Arguments

- object:

  A `netobject` (or `cograph_network`) object.

- ...:

  Ignored.

## Value

A `data.frame` with columns `metric` and `value`, of class
`c("summary.netobject", "data.frame")`.
