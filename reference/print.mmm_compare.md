# Print Method for mmm_compare

Print Method for mmm_compare

## Usage

``` r
# S3 method for class 'mmm_compare'
print(x, ...)
```

## Arguments

- x:

  An `mmm_compare` object.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
cmp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
print(cmp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -64.44053      162.8811 186.7014 188.8147 0.9655187 0.2118353 FALSE    
#>  3 -64.44054      180.8811 217.3122 221.4919 0.9331484 0.2601578 FALSE    
#>  best   
#>  <-- BIC
#>         
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 30, TRUE),
  V2 = sample(c("A","B","C"), 30, TRUE),
  V3 = sample(c("A","B","C"), 30, TRUE)
)
cmp <- compare_mmm(seqs, k = 2:3, n_starts = 5, seed = 1)
print(cmp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -89.29212      212.5842 236.4046 241.9893 0.9284601 0.1726947 TRUE     
#>  3 -87.61855      227.2371 263.6682 270.9108 0.9094695 0.1893851 TRUE     
#>  best   
#>  <-- BIC
#>         
# }
```
