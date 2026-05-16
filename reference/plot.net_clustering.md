# Plot Sequence Clustering Results

Plot Sequence Clustering Results

## Usage

``` r
# S3 method for class 'net_clustering'
plot(
  x,
  type = c("silhouette", "mds", "heatmap", "predictors"),
  combined = TRUE,
  ...
)
```

## Arguments

- x:

  A `net_clustering` object.

- type:

  Character. Plot type: `"silhouette"` (per-observation silhouette
  bars), `"mds"` (2D MDS projection), or `"heatmap"` (distance matrix
  heatmap ordered by cluster). Default: `"silhouette"`.

- ...:

  Unsupported. Supplying unused arguments raises an error.

## Value

A `ggplot` object (invisibly).

## Examples

``` r
seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
                   V3 = c("C","A","B","C","B"))
cl <- build_clusters(seqs, k = 2)
plot(cl, type = "silhouette")

# \donttest{
set.seed(1)
seqs <- data.frame(
  V1 = sample(c("A","B","C"), 20, TRUE),
  V2 = sample(c("A","B","C"), 20, TRUE),
  V3 = sample(c("A","B","C"), 20, TRUE)
)
cl <- build_clusters(seqs, k = 2)
plot(cl, type = "silhouette")

# }
```
