# Print Method for net_clustering

Compact, fixed-width summary of a sequence-clustering result. The header
carries the clustering method and dissimilarity; the per-cluster table
carries cluster size (count and percentage) and mean within-cluster
distance when available. Optional medoid and covariate lines surface
only when those fields are populated.

## Usage

``` r
# S3 method for class 'net_clustering'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `net_clustering` object.

- digits:

  Integer. Decimal places used for floating-point statistics in the
  printout. Default `3`. Non-breaking: existing `print(x)` calls keep
  their previous formatting.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
cl <- build_clusters(seqs, k = 2)
print(cl)
#> Sequence Clustering [pam]
#>   Sequences: 5  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.600
#> 
#>   Cluster  N          Mean within-dist  Medoid
#>   1        2 (40.0%)  0.000             4
#>   2        3 (60.0%)  2.000             5
# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 20, TRUE),
  V2 = sample(c("A","B","C"), 20, TRUE),
  V3 = sample(c("A","B","C"), 20, TRUE)
)
cl <- build_clusters(seqs, k = 2)
print(cl)
#> Sequence Clustering [pam]
#>   Sequences: 20  |  Clusters: 2
#>   Dissimilarity: hamming
#>   Quality: silhouette = 0.341
#> 
#>   Cluster  N           Mean within-dist  Medoid
#>   1        11 (55.0%)  1.600             14
#>   2        9 (45.0%)   1.444             18
# }
```
