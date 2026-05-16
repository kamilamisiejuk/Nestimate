# Compare MMM fits across different k

Compare MMM fits across different k

## Usage

``` r
compare_mmm(data, k = 2:5, return_fits = FALSE, ...)
```

## Arguments

- data:

  Data frame, netobject, or tna model.

- k:

  Integer vector of component counts. Values must be whole finite
  numbers \>= 2. Default: 2:5.

- return_fits:

  Logical. When `TRUE` the fitted models are retained on the result via
  `attr(result, "fits")` (a list of `net_mmm` objects, named by `k`), so
  the user can pick the chosen model without re-running the EM. Default
  `FALSE` keeps the historical lightweight return shape – only the
  comparison table is allocated.

- ...:

  Arguments passed to
  [`build_mmm`](https://saqr.me/Nestimate/reference/build_mmm.md).

## Value

A `mmm_compare` data frame with BIC, AIC, ICL, AvePP, entropy per k.
When `return_fits = TRUE`, the fitted models are attached as
`attr(result, "fits")`.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
comp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
comp
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP    Entropy   converged
#>  2 -63.27791      160.5558 184.3762 186.5265 0.964899 0.2153937 FALSE    
#>  3 -63.27792      178.5558 214.9870 218.8336 0.938229 0.2455097  TRUE    
#>  best   
#>  <-- BIC
#>         
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
  V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
)
comp <- compare_mmm(seqs, k = 2:3, seed = 42)
print(comp)
#> MMM Model Comparison
#> 
#>  k log_likelihood AIC      BIC      ICL      AvePP     Entropy   converged
#>  2 -122.4385      278.8771 302.6974 310.2070 0.8974786 0.3160632 TRUE     
#>  3 -118.1764      288.3527 324.7838 332.3711 0.9006835 0.1882551 TRUE     
#>  best   
#>  <-- BIC
#>         

# Retain fits to avoid a re-fit after picking the BIC-min model.
comp_with_fits <- compare_mmm(seqs, k = 2:3, seed = 42, return_fits = TRUE)
best_k <- comp_with_fits$k[which.min(comp_with_fits$BIC)]
best_fit <- attr(comp_with_fits, "fits")[[as.character(best_k)]]
# }
```
