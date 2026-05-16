# Print Method for mcml

Print Method for mcml

## Usage

``` r
# S3 method for class 'mcml'
print(x, ...)
```

## Arguments

- x:

  An `mcml` object.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
clusters <- list(G1 = c("A","B"), G2 = c("C"))
cs <- build_mcml(seqs, clusters)
print(cs)
#> MCML Network
#> ============
#> Type: tna  | Method: sum 
#> Nodes: 3  | Clusters: 2 
#> Transitions: 4 
#>   Macro: 2  | Per-cluster: 2 
#> 
#> Clusters:
#>   G1 (2): A, B
#>   G2 (1): C
#> 
#> Macro (cluster-level) weights:
#>        G1     G2
#> G1 0.6667 0.3333
#> G2 1.0000 0.0000
# \donttest{
seqs <- data.frame(
  T1 = c("A","B","A"), T2 = c("B","C","B"),
  T3 = c("C","A","C"), T4 = c("A","B","A")
)
clusters <- c("Alpha", "Beta", "Alpha")
cs <- build_mcml(seqs, clusters, type = "raw")
print(cs)
#> MCML Network
#> ============
#> Type: raw  | Method: sum 
#> Nodes: 3  | Clusters: 2 
#> Transitions: 9 
#>   Macro: 6  | Per-cluster: 3 
#> 
#> Clusters:
#>   Alpha (2): A, C
#>   Beta (1): B
#> 
#> Macro (cluster-level) weights:
#>       Alpha Beta
#> Alpha     3    3
#> Beta      3    0
# }
```
