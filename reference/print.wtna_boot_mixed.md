# Print Method for wtna_boot_mixed

Print Method for wtna_boot_mixed

## Usage

``` r
# S3 method for class 'wtna_boot_mixed'
print(x, ...)
```

## Arguments

- x:

  A `wtna_boot_mixed` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 10)
print(boot)
#> Mixed Window TNA Bootstrap
#> -- Transition --
#> Bootstrap Network  [Window TNA | directed]
#>   Iterations : 10  |  Nodes : 3
#>   Edges      : 0 significant / 3 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
#> -- Co-occurrence --
#> Bootstrap Network  [Network (wtna_cooccurrence) | undirected]
#>   Iterations : 10  |  Nodes : 3
#>   Edges      : 0 significant / 6 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
# \donttest{
set.seed(1)
oh <- data.frame(
  A = c(1,0,1,0,1,0,1,0),
  B = c(0,1,0,1,0,1,0,1),
  C = c(1,1,0,0,1,1,0,0)
)
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 20)
print(boot)
#> Mixed Window TNA Bootstrap
#> -- Transition --
#> Bootstrap Network  [Window TNA | directed]
#>   Iterations : 20  |  Nodes : 3
#>   Edges      : 0 significant / 9 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
#> -- Co-occurrence --
#> Bootstrap Network  [Network (wtna_cooccurrence) | undirected]
#>   Iterations : 20  |  Nodes : 3
#>   Edges      : 0 significant / 6 total
#>   CI         : 95%  |  Inference: stability  |  CR [0.75, 1.25]
# }
```
