# Summary Method for mcml

Summary Method for mcml

## Usage

``` r
# S3 method for class 'mcml'
summary(object, ...)
```

## Arguments

- object:

  An `mcml` object.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A tidy data frame with one row per cluster and columns `cluster`,
`size`, `within_total`, `between_out`, `between_in`. For undirected
macro networks the in/out split is not meaningful, so `between_out`
reports total incident weight and `between_in` is `NA`. The data frame
is returned silently *without* printing the full object – call
`print(object)` explicitly if you want the verbose dump.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
clusters <- list(G1 = c("A","B"), G2 = c("C"))
cs <- build_mcml(seqs, clusters)
summary(cs)
#>   cluster size within_total between_out between_in
#> 1      G1    2            1   0.3333333  1.0000000
#> 2      G2    1            0   1.0000000  0.3333333
# \donttest{
seqs <- data.frame(
  T1 = c("A","B","A"), T2 = c("B","C","B"),
  T3 = c("C","A","C"), T4 = c("A","B","A")
)
clusters <- c("Alpha", "Beta", "Alpha")
cs <- build_mcml(seqs, clusters, type = "raw")
summary(cs)
#>   cluster size within_total between_out between_in
#> 1   Alpha    2            3           3          3
#> 2    Beta    1            0           3          3
# }
```
