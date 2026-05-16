# Plot Method for net_mmm

Plot Method for net_mmm

## Usage

``` r
# S3 method for class 'net_mmm'
plot(x, type = c("posterior", "covariates"), combined = TRUE, ...)
```

## Arguments

- x:

  A `net_mmm` object.

- type:

  Character. Plot type: `"posterior"` (default) or `"covariates"`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `ggplot` object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
plot(mmm, type = "posterior")

# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
plot(mmm, type = "posterior")

# }
```
