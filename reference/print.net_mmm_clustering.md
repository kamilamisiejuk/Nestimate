# Print Method for MMM Clustering Attribute

Prints the clustering metadata that
[`cluster_mmm`](https://saqr.me/Nestimate/reference/cluster_mmm.md)
attaches to its `netobject_group` return value
(`attr(grp, "clustering")`). Layout mirrors
[`print.net_clustering`](https://saqr.me/Nestimate/reference/print.net_clustering.md):
a one-line dimension header, a quality line with AvePP / entropy /
classification error, information criteria, and a per-cluster table.

## Usage

``` r
# S3 method for class 'net_mmm_clustering'
print(x, digits = 3L, ...)
```

## Arguments

- x:

  A `net_mmm_clustering` object.

- digits:

  Integer. Decimal places for floating-point statistics. Default `3`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

The input object, invisibly.

## Examples

``` r
seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
                   V2 = sample(c("A","B","C"), 30, TRUE))
grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
print(attr(grp, "clustering"))
#> MMM Clustering [k = 2]
#>   Sequences: 30  |  Clusters: 2
#>   Quality: AvePP = 0.960  |  Entropy = 0.234  |  Class.Err = 0.0%
#>   ICs: BIC = 178.831  |  AIC = 155.011  |  ICL = 181.336
#> 
#>   Cluster  N           Mix%   AvePP
#>   1        24 (80.0%)  79.5%  0.971
#>   2        6 (20.0%)   20.5%  0.912
```
