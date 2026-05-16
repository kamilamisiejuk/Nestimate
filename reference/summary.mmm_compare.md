# Summary Method for mmm_compare

Summary Method for mmm_compare

## Usage

``` r
# S3 method for class 'mmm_compare'
summary(object, ...)
```

## Arguments

- object:

  An `mmm_compare` object (a data.frame subclass).

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A tidy data frame with one row per `k`, plus a `best` character column
flagging the minimum-BIC and minimum-ICL solutions.
