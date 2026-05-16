# Summary Method for net_mmm

Summary Method for net_mmm

## Usage

``` r
# S3 method for class 'net_mmm'
summary(object, ...)
```

## Arguments

- object:

  A `net_mmm` object.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -65.033  |  BIC = 187.886  |  AIC = 164.066  |  ICL = 190.065
#>   Quality: AvePP = 0.964  |  Entropy = 0.217  |  Class.Err = 0.0%
#>   Status: did not converge in 10 iterations
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        25 (83.3%)  81.4%  0.966
#>   2        5 (16.7%)   18.6%  0.954
#> 
#> --- Cluster 1 (81.4%, n=25) ---
#>       A     B     C
#> A 0.333 0.333 0.333
#> B 0.556 0.411 0.033
#> C 0.447 0.220 0.333
#> 
#> --- Cluster 2 (18.6%, n=5) ---
#>       A     B     C
#> A 0.333 0.333 0.333
#> B 0.038 0.037 0.925
#> C 0.337 0.329 0.334
#> 
#>   component     prior n_assigned mean_posterior     avepp
#> 1         1 0.8137522         25      0.9664648 0.9664648
#> 2         2 0.1862478          5      0.9543665 0.9543665
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
summary(mmm)
#> Mixed Markov Model
#>   Sequences: 30  |  Clusters: 2  |  States: 3
#>   ICs: LL = -89.292  |  BIC = 236.405  |  AIC = 212.584  |  ICL = 241.989
#>   Quality: AvePP = 0.928  |  Entropy = 0.173  |  Class.Err = 0.0%
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  72.9%  0.911
#>   2        6 (20.0%)   27.1%  0.999
#> 
#> --- Cluster 1 (72.9%, n=24) ---
#>       A     B     C
#> A 0.001 0.624 0.375
#> B 0.430 0.219 0.351
#> C 0.599 0.001 0.400
#> 
#> --- Cluster 2 (27.1%, n=6) ---
#>       A     B     C
#> A 0.988 0.008 0.004
#> B 0.256 0.445 0.299
#> C 0.003 0.994 0.003
#> 
#>   component    prior n_assigned mean_posterior     avepp
#> 1         1 0.728925         24      0.9108722 0.9108722
#> 2         2 0.271075          6      0.9988113 0.9988113
# }
```
