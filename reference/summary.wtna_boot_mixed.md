# Summary Method for wtna_boot_mixed

Summary Method for wtna_boot_mixed

## Usage

``` r
# S3 method for class 'wtna_boot_mixed'
summary(object, ...)
```

## Arguments

- object:

  A `wtna_boot_mixed` object.

- ...:

  Additional arguments (ignored).

## Value

A list with `$transition` and `$cooccurrence` summary data frames.

## Examples

``` r
oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 10)
summary(boot)
#> $transition
#>   from to weight mean        sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  B      2  1.3 1.1595018 0.5454545 FALSE        0    2.775     1.50
#> 2    B  B      1  0.8 0.9189366 0.5454545 FALSE        0    2.550     0.75
#> 3    C  B      2  1.2 0.9189366 0.5454545 FALSE        0    2.000     1.50
#>   cr_upper
#> 1     2.50
#> 2     1.25
#> 3     2.50
#> 
#> $cooccurrence
#>   from to weight mean       sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  A      4  3.8 2.780887 1.0000000 FALSE    1.000    8.875      3.0
#> 2    A  B      2  1.6 0.843274 0.2727273 FALSE    0.000    2.000      1.5
#> 3    A  C      4  1.6 1.264911 0.9090909 FALSE    0.000    3.775      3.0
#> 4    B  B      2  3.0 2.748737 1.0000000 FALSE    0.225    8.100      1.5
#> 5    B  C      2  2.1 2.378141 0.9090909 FALSE    0.000    6.000      1.5
#> 6    C  C      4  2.7 2.790858 0.8181818 FALSE    0.225    8.100      3.0
#>   cr_upper
#> 1      5.0
#> 2      2.5
#> 3      5.0
#> 4      2.5
#> 5      2.5
#> 6      5.0
#> 
# \donttest{
set.seed(1)
oh <- data.frame(
  A = c(1,0,1,0,1,0,1,0),
  B = c(0,1,0,1,0,1,0,1),
  C = c(1,1,0,0,1,1,0,0)
)
mixed <- wtna(oh, method = "both")
boot  <- bootstrap_network(mixed, iter = 20)
summary(boot)
#> $transition
#>   from to weight mean       sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  A      3 3.65 2.739093 1.0000000 FALSE    0.000    9.050     2.25
#> 2    A  B      5 3.75 1.996708 0.6666667 FALSE    0.950    8.050     3.75
#> 3    A  C      4 4.35 2.978431 0.9523810 FALSE    0.000    9.525     3.00
#> 4    B  A      3 3.85 2.300458 0.9047619 FALSE    0.000    8.050     2.25
#> 5    B  B      4 3.75 3.258592 0.9523810 FALSE    0.000   11.050     3.00
#> 6    B  C      2 3.30 2.657660 0.9523810 FALSE    0.000    9.200     1.50
#> 7    C  A      4 4.20 3.053902 0.8571429 FALSE    0.000   10.525     3.00
#> 8    C  B      6 4.10 3.143916 0.8571429 FALSE    0.475   11.050     4.50
#> 9    C  C      4 4.20 3.778053 0.9523810 FALSE    0.000   13.575     3.00
#>   cr_upper
#> 1     3.75
#> 2     6.25
#> 3     5.00
#> 4     3.75
#> 5     5.00
#> 6     2.50
#> 7     5.00
#> 8     7.50
#> 9     5.00
#> 
#> $cooccurrence
#>   from to weight mean       sd   p_value   sig ci_lower ci_upper cr_lower
#> 1    A  A      6 7.85 4.356181 0.8095238 FALSE    1.475   16.625     4.50
#> 2    A  B      5 3.40 1.429022 0.3809524 FALSE    0.475    5.000     3.75
#> 3    A  C      6 6.00 3.308681 0.6190476 FALSE    2.000   13.525     4.50
#> 4    B  B      6 7.35 4.498830 0.8571429 FALSE    1.475   15.100     4.50
#> 5    B  C      6 5.95 3.677456 0.8571429 FALSE    1.475   13.525     4.50
#> 6    C  C      8 8.00 4.142209 0.7619048 FALSE    2.475   14.000     6.00
#>   cr_upper
#> 1     7.50
#> 2     6.25
#> 3     7.50
#> 4     7.50
#> 5     7.50
#> 6    10.00
#> 
# }
```
