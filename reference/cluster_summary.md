# Cluster Summary Statistics

Aggregates node-level network weights to cluster-level summaries.
Computes both between-cluster transitions (how clusters connect to each
other) and within-cluster transitions (how nodes connect within each
cluster).

## Usage

``` r
cluster_summary(
  x,
  clusters = NULL,
  method = c("sum", "mean", "median", "max", "min", "density", "geomean"),
  directed = TRUE,
  compute_within = TRUE
)
```

## Arguments

- x:

  Network input. Accepts multiple formats:

  matrix

  :   Numeric adjacency/weight matrix. Row and column names are used as
      node labels. Values represent edge weights (e.g., transition
      counts, co-occurrence frequencies, or probabilities).

  netobject

  :   A cograph network object. The function extracts the weight matrix
      from `x$weights` or converts via `to_matrix()`. Clusters can be
      auto-detected from node attributes.

  tna

  :   A tna object from the tna package. Extracts `x$weights`.

  cluster_summary

  :   If already a cluster_summary, returns unchanged.

- clusters:

  Cluster/group assignments for nodes. Accepts multiple formats:

  NULL

  :   (default) Auto-detect from netobject. Looks for columns named
      'clusters', 'cluster', 'groups', or 'group' in `x$nodes`. Throws
      an error if no cluster column is found. This option only works
      when `x` is a netobject.

  vector

  :   Cluster membership for each node, in the same order as the matrix
      rows/columns. Can be numeric (1, 2, 3) or character ("A", "B").
      Cluster names will be derived from unique values. Example:
      `c(1, 1, 2, 2, 3, 3)` assigns first two nodes to cluster 1.

  data.frame

  :   A data frame where the first column contains node names and the
      second column contains group/cluster names. Example:
      `data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))`

  named list

  :   Explicit mapping of cluster names to node labels. List names
      become cluster names, values are character vectors of node labels
      that must match matrix row/column names. Example:
      `list(Alpha = c("A", "B"), Beta = c("C", "D"))`

- method:

  Aggregation method for combining edge weights within/between clusters.
  Controls how multiple node-to-node edges are summarized:

  "sum"

  :   (default) Sum of all edge weights. Best for count data (e.g.,
      transition frequencies). Preserves total flow.

  "mean"

  :   Average edge weight. Best when cluster sizes differ and you want
      to control for size. Note: when input is already a transition
      matrix (rows sum to 1), "mean" avoids size bias. Example: cluster
      with 5 nodes won't have 5x the weight of cluster with 1 node.

  "median"

  :   Median edge weight. Robust to outliers.

  "max"

  :   Maximum edge weight. Captures strongest connection.

  "min"

  :   Minimum edge weight. Captures weakest connection.

  "density"

  :   Sum divided by number of possible edges. Normalizes by cluster
      size combinations.

  "geomean"

  :   Geometric mean of positive weights. Useful for multiplicative
      processes.

- directed:

  Logical. If `TRUE` (default), treat network as directed. A-\>B and
  B-\>A are separate edges. If `FALSE`, edges are undirected and the
  matrix is symmetrized before processing.

- compute_within:

  Logical. If `TRUE` (default), compute within-cluster transition
  matrices for each cluster. Each cluster gets its own n_i x n_i matrix
  showing internal node-to-node transitions. Set to `FALSE` to skip this
  computation for better performance when only between-cluster summary
  is needed.

## Value

A `cluster_summary` object (S3 class) containing:

- between:

  List with two elements:

  weights

  :   k x k matrix of cluster-to-cluster weights, where k is the number
      of clusters. Row i, column j contains the elementwise aggregation
      (per `method`) of all edges from nodes in cluster i to nodes in
      cluster j. Diagonal contains within-cluster totals. Pure
      arithmetic – no row normalization.

  inits

  :   Numeric vector of length k. Initial state distribution across
      clusters, computed from column sums of the original matrix.
      Represents the proportion of incoming edges to each cluster.

- within:

  Named list with one element per cluster. Each element contains:

  weights

  :   n_i x n_i matrix for nodes within that cluster. Shows internal
      transitions between nodes in the same cluster.

  inits

  :   Initial distribution within the cluster.

  NULL if `compute_within = FALSE`.

- clusters:

  Named list mapping cluster names to their member node labels. Example:
  `list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))`

- meta:

  List of metadata:

  method

  :   The `method` argument used ("sum", "mean", etc.)

  directed

  :   Logical, whether network was treated as directed

  n_nodes

  :   Total number of nodes in original network

  n_clusters

  :   Number of clusters

  cluster_sizes

  :   Named vector of cluster sizes

## Details

This is the core function for Multi-Cluster Multi-Level (MCML) analysis.
Use
[`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md) to
convert results to tna objects for further analysis with the tna
package.

### Workflow

Typical MCML analysis workflow:


    # 1. Create network
    net <- build_network(data, method = "relative")
    net$nodes$clusters <- group_assignments

    # 2. Compute cluster summary (arithmetic aggregation over edges)
    cs <- cluster_summary(net, method = "sum")

    # 3. Convert to tna models (normalization happens in as_tna)
    tna_models <- as_tna(cs)

    # 4. Analyze/visualize
    plot(tna_models$macro)
    tna::centralities(tna_models$macro)

### Between-Cluster Matrix Structure

The `macro$weights` matrix has clusters as both rows and columns:

- Off-diagonal (row i, col j): Aggregated weight from cluster i to
  cluster j

- Diagonal (row i, col i): Within-cluster total (aggregation of internal
  edges)

Rows are NOT normalized. Entries are elementwise aggregates produced by
`method`. If the caller wants probabilities, they should normalize
downstream (e.g. via
[`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md)).
Mixing an arithmetic aggregation with row-normalization here (the old
`type = "tna"` combined with `method = "min"` / `"mean"` etc.) produces
numbers that sum to 1 per row but are not a probability distribution
over any process; that silently-wrong combination is why `type` was
removed from the matrix path. The sequence and edgelist paths of
[`build_mcml()`](https://mohsaqr.github.io/Nestimate/reference/build_mcml.md)
keep `type`, where the aggregation is always counts and the
post-processing chooses between well-defined network constructions.

### Choosing method

|  |  |  |
|----|----|----|
| **Input data** | **Recommended method** | **Reason** |
| Edge counts | `"sum"` | Preserves total flow between clusters |
| Transition matrix | `"mean"` | Avoids cluster size bias |
| Correlation matrix | `"mean"` | Average correlations |
| Dense weighted | `"max"` / `"median"` | Robust summary |

## See also

[`as_tna()`](https://mohsaqr.github.io/Nestimate/reference/as_tna.md) to
convert results to tna objects,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for two-layer
visualization, [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
for flat cluster visualization

## Examples

``` r
# -----------------------------------------------------
# Basic usage with matrix and cluster vector
# -----------------------------------------------------
mat <- matrix(runif(100), 10, 10)
rownames(mat) <- colnames(mat) <- LETTERS[1:10]

clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
cs <- cluster_summary(mat, clusters)

# Access results
cs$macro$weights    # 3x3 cluster transition matrix
#>          1        2        3
#> 1 4.784827 5.331427 6.186107
#> 2 3.368078 4.033927 6.611524
#> 3 5.279989 6.773042 9.423943
cs$macro$inits      # Initial distribution
#>        1        2        3 
#> 0.259358 0.311595 0.429047 
cs$clusters$`1`$weights # Within-cluster 1 transitions
#>           A         B         C
#> A 0.7720692 0.6796373 0.5709870
#> B 0.6576550 0.3574632 0.7623484
#> C 0.7057731 0.1342718 0.1446221
cs$meta               # Metadata
#> $type
#> [1] "aggregate"
#> 
#> $method
#> [1] "sum"
#> 
#> $directed
#> [1] TRUE
#> 
#> $n_nodes
#> [1] 10
#> 
#> $n_clusters
#> [1] 3
#> 
#> $cluster_sizes
#> 1 2 3 
#> 3 3 4 
#> 
#> $source
#> [1] "matrix"
#> 

# -----------------------------------------------------
# Named list clusters (more readable)
# -----------------------------------------------------
clusters <- list(
  Alpha = c("A", "B", "C"),
  Beta = c("D", "E", "F"),
  Gamma = c("G", "H", "I", "J")
)
cs <- cluster_summary(mat, clusters)
cs$macro$weights    # Rows/cols named Alpha, Beta, Gamma
#>          Alpha     Beta    Gamma
#> Alpha 4.784827 5.331427 6.186107
#> Beta  3.368078 4.033927 6.611524
#> Gamma 5.279989 6.773042 9.423943
cs$clusters$Alpha       # Within Alpha cluster
#> MCML layer (transition probabilities)  [directed]
#>   Nodes: 3  |  Non-zero edges: 9
#>   Weights: [0.134, 0.772]  |  mean: 0.532
#> 
#>   Weight matrix:
#>         A     B     C
#>   A 0.772 0.680 0.571
#>   B 0.658 0.357 0.762
#>   C 0.706 0.134 0.145 
#> 
#>   Initial probabilities:
#>   A               0.446  ########################################
#>   C               0.309  ############################
#>   B               0.245  ######################

# -----------------------------------------------------
# Auto-detect clusters from netobject
# -----------------------------------------------------
# \donttest{
seqs <- data.frame(
  V1 = sample(LETTERS[1:10], 30, TRUE), V2 = sample(LETTERS[1:10], 30, TRUE),
  V3 = sample(LETTERS[1:10], 30, TRUE)
)
net <- build_network(seqs, method = "relative")
cs2 <- cluster_summary(net, c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
# }

# -----------------------------------------------------
# Different aggregation methods
# -----------------------------------------------------
cs_sum <- cluster_summary(mat, clusters, method = "sum")   # Total flow
cs_mean <- cluster_summary(mat, clusters, method = "mean") # Average
cs_max <- cluster_summary(mat, clusters, method = "max")   # Strongest

# -----------------------------------------------------
# Skip within-cluster computation for speed
# -----------------------------------------------------
cs_fast <- cluster_summary(mat, clusters, compute_within = FALSE)
cs_fast$clusters  # NULL
#> NULL

# -----------------------------------------------------
# Convert to tna objects for tna package
# (as_tna() applies its own row normalisation)
# -----------------------------------------------------
cs <- cluster_summary(mat, clusters, method = "sum")
tna_models <- as_tna(cs)
# tna_models$macro      # tna object
# tna_models$clusters$Alpha # tna object
```
