# Print Method for net_mmm

Compact summary of a Mixed Markov Model fit. Header carries dimensions
and information criteria; cluster table carries N, mixing share, and
per-cluster average posterior probability (AvePP). Layout matches
[`print.net_clustering`](https://saqr.me/Nestimate/reference/print.net_clustering.md)
so distance- and model-based clusterings can be compared at a glance.

## Usage

``` r
# S3 method for class 'net_mmm'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `net_mmm` object.

- digits:

  Integer. Decimal places for floating-point statistics. Default `3`.
  Non-breaking: `print(x)` keeps the same alignment as before.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
print(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -65.033  |  BIC = 187.886  |  AIC = 164.066  |  ICL = 190.065
#>   Quality: AvePP = 0.964  |  Entropy = 0.217  |  Class.Err = 0.0%
#>   Status: did not converge in 10 iterations
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        25 (83.3%)  81.4%  0.966
#>   2        5 (16.7%)   18.6%  0.954
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
print(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -89.292  |  BIC = 236.405  |  AIC = 212.584  |  ICL = 241.989
#>   Quality: AvePP = 0.928  |  Entropy = 0.173  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  72.9%  0.911
#>   2        6 (20.0%)   27.1%  0.999
# }
```
