# Evaluate Link Predictions Against Known Edges

Computes AUC-ROC, precision\\@\\k, and average precision for link
predictions against a set of known true edges.

## Usage

``` r
evaluate_links(pred, true_edges, k = c(5L, 10L, 20L))
```

## Arguments

- pred:

  A `net_link_prediction` object.

- true_edges:

  A data frame with columns `from` and `to`, or a binary matrix where 1
  indicates a true edge.

- k:

  Integer vector. Values of k for precision\\@\\k. Default:
  `c(5, 10, 20)`.

## Value

A data frame with columns: method, auc, average_precision, and one
precision_at_k column per k value.

## Examples

``` r
set.seed(42)
seqs <- data.frame(
  V1 = sample(LETTERS[1:5], 50, TRUE),
  V2 = sample(LETTERS[1:5], 50, TRUE),
  V3 = sample(LETTERS[1:5], 50, TRUE)
)
net <- build_network(seqs, method = "relative")
pred <- predict_links(net, exclude_existing = FALSE)

# Evaluate: predict the network's own edges
true <- data.frame(from = pred$predictions$from[1:5],
                   to = pred$predictions$to[1:5])
evaluate_links(pred, true)
#>                    method       auc average_precision precision_at_5
#> 1        common_neighbors 0.9933333         1.0000000            1.0
#> 2     resource_allocation 0.9933333         1.0000000            1.0
#> 3             adamic_adar 0.9933333         1.0000000            1.0
#> 4                 jaccard 0.9400000         0.8100000            0.8
#> 5 preferential_attachment 0.6066667         0.5110731            0.4
#> 6                    katz 0.6400000         0.3671032            0.4
#>   precision_at_10 precision_at_20
#> 1             0.5            0.25
#> 2             0.5            0.25
#> 3             0.5            0.25
#> 4             0.5            0.25
#> 5             0.3            0.25
#> 6             0.4            0.25
```
