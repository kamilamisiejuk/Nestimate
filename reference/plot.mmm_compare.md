# Plot Method for mmm_compare

Plot Method for mmm_compare

## Usage

``` r
# S3 method for class 'mmm_compare'
plot(x, ...)
```

## Arguments

- x:

  An `mmm_compare` object.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `ggplot` object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
cmp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
plot(cmp)

# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
cmp <- compare_mmm(seqs, k = 2:3, n_starts = 5, seed = 1)
plot(cmp)

# }
```
