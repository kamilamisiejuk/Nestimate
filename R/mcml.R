# Cluster Metrics for Network Analysis
# Summary measures for between/within clusters and multilayer networks


# Tagged-list constructor for an mcml layer (macro or one cluster). The
# class buys us a `print.mcml_layer` method so `print(mc$macro)` shows a
# tidy summary instead of dumping every named element of the list.
# Field access (`$weights`, `$inits`, `$labels`, `$data`) is unchanged.
.mcml_layer <- function(weights, inits, labels, data = NULL) {
  if (!is.matrix(weights) || !is.numeric(weights)) {
    stop("'weights' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(weights) != ncol(weights)) {
    stop("'weights' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(weights)
  if (!is.character(labels) || length(labels) != nrow(weights) ||
      any(is.na(labels)) || any(!nzchar(labels))) {
    stop("'labels' must be a non-missing character vector with one value per node.",
         call. = FALSE)
  }
  if (anyDuplicated(labels)) {
    stop("'labels' must be unique.", call. = FALSE)
  }
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(labels) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per node.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(labels, names(inits))
      extra_inits <- setdiff(names(inits), labels)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match layer labels.", call. = FALSE)
      }
      inits <- inits[labels]
    } else {
      names(inits) <- labels
    }
  }
  obj <- list(weights = weights, inits = inits, labels = labels, data = data)
  class(obj) <- "mcml_layer"
  obj
}


#' Print Method for an mcml Layer
#'
#' Compact summary of one mcml layer (macro or a single within-cluster
#' network) -- nodes, edges, weight matrix dimensions -- without spilling
#' the full \code{$weights}, \code{$inits}, or \code{$data} contents.
#'
#' @param x An \code{mcml_layer}.
#' @param ... Unused.
#' @return The input, invisibly.
#' @export
print.mcml_layer <- function(x, ...) {
  .mcml_check_unused_dots("print.mcml_layer", ...)
  w  <- x$weights
  if (is.null(w) || nrow(w) == 0L) {
    cat("MCML layer  [empty]\n"); return(invisible(x))
  }

  cat("MCML layer (transition probabilities)  [directed]\n")
  cat(sprintf("  Nodes: %d  |  Non-zero edges: %d\n",
              nrow(w), sum(w != 0)))

  # ---- Weight summary line ----
  nz <- w[w != 0]
  if (length(nz) > 0) {
    cat(sprintf("  Weights: [%.3f, %.3f]  |  mean: %.3f\n",
                min(nz), max(nz), mean(nz)))
  }

  # ---- Weight matrix ----
  cat("\n  Weight matrix:\n")
  digits <- if (all(nz == floor(nz))) 0L else 3L
  mat_r <- round(w, digits)
  if (!is.null(x$labels)) dimnames(mat_r) <- list(x$labels, x$labels)
  formatted <- utils::capture.output(print(mat_r))
  cat(paste0("  ", formatted, collapse = "\n"), "\n")

  # ---- Initial probabilities (bar plot, same style as netobject) ----
  if (!is.null(x$inits) && length(x$inits) > 0L) {
    cat("\n  Initial probabilities:\n")
    init  <- x$inits
    nm    <- names(init)
    if (is.null(nm) && !is.null(x$labels)) nm <- x$labels
    if (is.null(nm)) nm <- as.character(seq_along(init))
    ord   <- order(init, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(init, na.rm = TRUE)
    for (i in ord) {
      bars <- if (is.finite(max_v) && max_v > 0)
        strrep("#", round(init[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-14s  %.3f  %s\n", nm[i], init[i], bars))
    }
  }

  # ---- Data dimensions ----
  d <- x$data
  if (!is.null(d)) {
    data_dim <- if (is.data.frame(d) || is.matrix(d))
      sprintf("%d x %d", nrow(d), ncol(d))
    else sprintf("length %d", length(d))
    cat(sprintf("\n  Data: %s\n", data_dim))
  }

  invisible(x)
}

# ==============================================================================
# 1. Edge Weight Aggregation
# ==============================================================================

#' Aggregate Edge Weights
#'
#' Aggregates a vector of edge weights using various methods.
#' Compatible with igraph's edge.attr.comb parameter.
#'
#' @param w Numeric vector of finite edge weights. \code{NA} and zero weights
#'   are excluded before aggregation.
#' @param method Single aggregation method: "sum", "mean", "median", "max",
#'   "min", "prod", "density", or "geomean".
#' @param n_possible Optional single finite numeric number of possible edges
#'   for density calculation.
#' @return Single aggregated value
#' @export
#' @examples
#' w <- c(0.5, 0.8, 0.3, 0.9)
#' net_aggregate_weights(w, "sum")   # 2.5
#' net_aggregate_weights(w, "mean")  # 0.625
#' net_aggregate_weights(w, "max")   # 0.9
#' net_aggregate_weights(w, "density", n_possible = 9)  # 2.5 / 9
net_aggregate_weights <- function(w, method = "sum", n_possible = NULL) {
  if (!is.numeric(w)) {
    stop("'w' must be a numeric vector.", call. = FALSE)
  }
  if (any(!is.na(w) & !is.finite(w))) {
    stop("'w' must contain only finite values or NA.", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("'method' must be a single non-missing character value.", call. = FALSE)
  }
  valid_methods <- c("sum", "mean", "median", "max", "min",
                     "prod", "density", "geomean")
  if (!method %in% valid_methods) {
    stop("Unknown method: ", method, call. = FALSE)
  }
  if (!is.null(n_possible) &&
      (!is.numeric(n_possible) || length(n_possible) != 1L ||
       is.na(n_possible) || !is.finite(n_possible))) {
    stop("'n_possible' must be a single finite numeric value or NULL.",
         call. = FALSE)
  }

  # Remove NA and zero weights
  w <- w[!is.na(w) & w != 0]
  if (length(w) == 0) return(0)

  switch(method,
    "sum"     = sum(w),
    "mean"    = mean(w),
    "median"  = stats::median(w),
    "max"     = max(w),
    "min"     = min(w),
    "prod"    = prod(w),
    "density" = if (!is.null(n_possible) && n_possible > 0) {
      sum(w) / n_possible
    } else {
      sum(w) / length(w)
    },
    "geomean" = {
      pos_w <- w[w > 0]
      if (length(pos_w) == 0) 0 else exp(mean(log(pos_w)))
    },
    stop("Unknown method: ", method, call. = FALSE)
  )
}


# ==============================================================================
# 2. Cluster Summary (Between/Within Aggregates)
# ==============================================================================

#' Cluster Summary Statistics
#'
#' Aggregates node-level network weights to cluster-level summaries. Computes
#' both between-cluster transitions (how clusters connect to each other) and
#' within-cluster transitions (how nodes connect within each cluster).
#'
#' This is the core function for Multi-Cluster Multi-Level (MCML) analysis.
#' Use \code{as_tna()} to convert results to tna objects for further
#' analysis with the tna package.
#'
#' @param x Network input. Accepts multiple formats:
#'   \describe{
#'     \item{matrix}{Numeric adjacency/weight matrix. Row and column names are
#'       used as node labels. Values represent edge weights (e.g., transition
#'       counts, co-occurrence frequencies, or probabilities).}
#'     \item{netobject}{A cograph network object. The function extracts
#'       the weight matrix from \code{x$weights} or converts via
#'       \code{to_matrix()}. Clusters can be auto-detected from node attributes.}
#'     \item{tna}{A tna object from the tna package. Extracts \code{x$weights}.}
#'     \item{cluster_summary}{If already a cluster_summary, returns unchanged.}
#'   }
#'
#' @param clusters Cluster/group assignments for nodes. Accepts multiple formats:
#'   \describe{
#'     \item{NULL}{(default) Auto-detect from netobject. Looks for columns
#'       named 'clusters', 'cluster', 'groups', or 'group' in \code{x$nodes}.
#'       Throws an error if no cluster column is found.
#'       This option only works when \code{x} is a netobject.}
#'     \item{vector}{Cluster membership for each node, in the same order as the

#'       matrix rows/columns. Can be numeric (1, 2, 3) or character ("A", "B").
#'       Cluster names will be derived from unique values.
#'       Example: \code{c(1, 1, 2, 2, 3, 3)} assigns first two nodes to cluster 1.}
#'     \item{data.frame}{A data frame where the first column contains node names
#'       and the second column contains group/cluster names.
#'       Example: \code{data.frame(node = c("A", "B", "C"), group = c("G1", "G1", "G2"))}}
#'     \item{named list}{Explicit mapping of cluster names to node labels.
#'       List names become cluster names, values are character vectors of node
#'       labels that must match matrix row/column names.
#'       Example: \code{list(Alpha = c("A", "B"), Beta = c("C", "D"))}}
#'   }
#'
#' @param method Aggregation method for combining edge weights within/between
#'   clusters. Controls how multiple node-to-node edges are summarized:
#'   \describe{
#'     \item{"sum"}{(default) Sum of all edge weights. Best for count data
#'       (e.g., transition frequencies). Preserves total flow.}
#'     \item{"mean"}{Average edge weight. Best when cluster sizes differ and
#'       you want to control for size. Note: when input is already a transition
#'       matrix (rows sum to 1), "mean" avoids size bias.
#'       Example: cluster with 5 nodes won't have 5x the weight of cluster with 1 node.}
#'     \item{"median"}{Median edge weight. Robust to outliers.}
#'     \item{"max"}{Maximum edge weight. Captures strongest connection.}
#'     \item{"min"}{Minimum edge weight. Captures weakest connection.}
#'     \item{"density"}{Sum divided by number of possible edges. Normalizes
#'       by cluster size combinations.}
#'     \item{"geomean"}{Geometric mean of positive weights. Useful for
#'       multiplicative processes.}
#'   }
#'
#' @param directed Logical. If \code{TRUE} (default), treat network as directed.
#'   A->B and B->A are separate edges. If \code{FALSE}, edges are undirected
#'   and the matrix is symmetrized before processing.
#'
#' @param compute_within Logical. If \code{TRUE} (default), compute within-cluster
#'   transition matrices for each cluster. Each cluster gets its own n_i x n_i
#'   matrix showing internal node-to-node transitions.
#'   Set to \code{FALSE} to skip this computation for better performance when
#'   only between-cluster summary is needed.
#'
#' @return A \code{cluster_summary} object (S3 class) containing:
#'   \describe{
#'     \item{between}{List with two elements:
#'       \describe{
#'         \item{weights}{k x k matrix of cluster-to-cluster weights, where k is
#'           the number of clusters. Row i, column j contains the elementwise
#'           aggregation (per \code{method}) of all edges from nodes in cluster
#'           i to nodes in cluster j. Diagonal contains within-cluster totals.
#'           Pure arithmetic -- no row normalization.}
#'         \item{inits}{Numeric vector of length k. Initial state distribution
#'           across clusters, computed from column sums of the original matrix.
#'           Represents the proportion of incoming edges to each cluster.}
#'       }
#'     }
#'     \item{within}{Named list with one element per cluster. Each element contains:
#'       \describe{
#'         \item{weights}{n_i x n_i matrix for nodes within that cluster.
#'           Shows internal transitions between nodes in the same cluster.}
#'         \item{inits}{Initial distribution within the cluster.}
#'       }
#'       NULL if \code{compute_within = FALSE}.}
#'     \item{clusters}{Named list mapping cluster names to their member node labels.
#'       Example: \code{list(A = c("n1", "n2"), B = c("n3", "n4", "n5"))}}
#'     \item{meta}{List of metadata:
#'       \describe{
#'         \item{method}{The \code{method} argument used ("sum", "mean", etc.)}
#'         \item{directed}{Logical, whether network was treated as directed}
#'         \item{n_nodes}{Total number of nodes in original network}
#'         \item{n_clusters}{Number of clusters}
#'         \item{cluster_sizes}{Named vector of cluster sizes}
#'       }
#'     }
#'   }
#'
#' @details
#' ## Workflow
#'
#' Typical MCML analysis workflow:
#' \preformatted{
#' # 1. Create network
#' net <- build_network(data, method = "relative")
#' net$nodes$clusters <- group_assignments
#'
#' # 2. Compute cluster summary (arithmetic aggregation over edges)
#' cs <- cluster_summary(net, method = "sum")
#'
#' # 3. Convert to tna models (normalization happens in as_tna)
#' tna_models <- as_tna(cs)
#'
#' # 4. Analyze/visualize
#' plot(tna_models$macro)
#' tna::centralities(tna_models$macro)
#' }
#'
#' ## Between-Cluster Matrix Structure
#'
#' The \code{macro$weights} matrix has clusters as both rows and columns:
#' \itemize{
#'   \item Off-diagonal (row i, col j): Aggregated weight from cluster i to cluster j
#'   \item Diagonal (row i, col i): Within-cluster total (aggregation of internal edges)
#' }
#'
#' Rows are NOT normalized. Entries are elementwise aggregates produced by
#' \code{method}. If the caller wants probabilities, they should normalize
#' downstream (e.g. via \code{as_tna()}). Mixing an arithmetic aggregation
#' with row-normalization here (the old \code{type = "tna"} combined with
#' \code{method = "min"} / \code{"mean"} etc.) produces numbers that sum to 1
#' per row but are not a probability distribution over any process; that
#' silently-wrong combination is why \code{type} was removed from the matrix
#' path. The sequence and edgelist paths of \code{build_mcml()} keep
#' \code{type}, where the aggregation is always counts and the post-processing
#' chooses between well-defined network constructions.
#'
#' ## Choosing method
#'
#' \tabular{lll}{
#'   \strong{Input data} \tab \strong{Recommended method} \tab \strong{Reason} \cr
#'   Edge counts \tab \code{"sum"} \tab Preserves total flow between clusters \cr
#'   Transition matrix \tab \code{"mean"} \tab Avoids cluster size bias \cr
#'   Correlation matrix \tab \code{"mean"} \tab Average correlations \cr
#'   Dense weighted \tab \code{"max"} / \code{"median"} \tab Robust summary \cr
#' }
#'
#' @export
#' @seealso
#'   \code{as_tna()} to convert results to tna objects,
#'   \code{plot()} for two-layer visualization,
#'   \code{plot()} for flat cluster visualization
#'
#' @examples
#' # -----------------------------------------------------
#' # Basic usage with matrix and cluster vector
#' # -----------------------------------------------------
#' mat <- matrix(runif(100), 10, 10)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:10]
#'
#' clusters <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3)
#' cs <- cluster_summary(mat, clusters)
#'
#' # Access results
#' cs$macro$weights    # 3x3 cluster transition matrix
#' cs$macro$inits      # Initial distribution
#' cs$clusters$`1`$weights # Within-cluster 1 transitions
#' cs$meta               # Metadata
#'
#' # -----------------------------------------------------
#' # Named list clusters (more readable)
#' # -----------------------------------------------------
#' clusters <- list(
#'   Alpha = c("A", "B", "C"),
#'   Beta = c("D", "E", "F"),
#'   Gamma = c("G", "H", "I", "J")
#' )
#' cs <- cluster_summary(mat, clusters)
#' cs$macro$weights    # Rows/cols named Alpha, Beta, Gamma
#' cs$clusters$Alpha       # Within Alpha cluster
#'
#' # -----------------------------------------------------
#' # Auto-detect clusters from netobject
#' # -----------------------------------------------------
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:10], 30, TRUE), V2 = sample(LETTERS[1:10], 30, TRUE),
#'   V3 = sample(LETTERS[1:10], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' cs2 <- cluster_summary(net, c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3))
#' }
#'
#' # -----------------------------------------------------
#' # Different aggregation methods
#' # -----------------------------------------------------
#' cs_sum <- cluster_summary(mat, clusters, method = "sum")   # Total flow
#' cs_mean <- cluster_summary(mat, clusters, method = "mean") # Average
#' cs_max <- cluster_summary(mat, clusters, method = "max")   # Strongest
#'
#' # -----------------------------------------------------
#' # Skip within-cluster computation for speed
#' # -----------------------------------------------------
#' cs_fast <- cluster_summary(mat, clusters, compute_within = FALSE)
#' cs_fast$clusters  # NULL
#'
#' # -----------------------------------------------------
#' # Convert to tna objects for tna package
#' # (as_tna() applies its own row normalisation)
#' # -----------------------------------------------------
#' cs <- cluster_summary(mat, clusters, method = "sum")
#' tna_models <- as_tna(cs)
#' # tna_models$macro      # tna object
#' # tna_models$clusters$Alpha # tna object
cluster_summary <- function(x,
                            clusters = NULL,
                            method = c("sum", "mean", "median", "max",
                                       "min", "density", "geomean"),
                            directed = TRUE,
                            compute_within = TRUE) {
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop("'compute_within' must be TRUE or FALSE.", call. = FALSE)
  }

  # If already an mcml object, return as-is
  if (inherits(x, "mcml")) {
    return(x) # nocov
  }

  method <- match.arg(method)

  # Store original for cluster extraction
  x_orig <- x

  # Extract matrix from various input types
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
    mat <- x$weights
  } else if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    mat <- x$weights
  } else if (inherits(x, "tna")) {
    mat <- x$weights
  } else {
    mat <- x
  }

  if (is.null(clusters)) {
    stop("clusters argument is required", call. = FALSE)
  }

  # Validate input matrix
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("x must be a netobject, tna object, or numeric matrix", call. = FALSE) # nocov
  }
  if (nrow(mat) != ncol(mat)) {
    stop("x must be a square matrix", call. = FALSE)
  }
  .validate_mcml_matrix(mat)

  # Symmetrize matrix when treating it as undirected, matching the docstring
  # contract. Without this, callers passing `directed = FALSE` got the same
  # output as `directed = TRUE`, which silently violated the contract.
  if (!isTRUE(directed)) mat <- (mat + t(mat)) / 2

  n <- nrow(mat)
  node_names <- rownames(mat)
  if (is.null(node_names)) node_names <- as.character(seq_len(n))

  # Convert clusters to list format
  cluster_list <- .normalize_clusters(clusters, node_names)
  n_clusters <- length(cluster_list)
  cluster_names <- names(cluster_list)
  names(cluster_list) <- cluster_names

  # Get node indices for each cluster
  cluster_indices <- lapply(cluster_list, function(nodes_vec) {
    match(nodes_vec, node_names)
  })

  # ============================================================================
  # Between-cluster computation (always computed)
  # ============================================================================

  # Aggregate between-cluster weights
  between_raw <- matrix(0, n_clusters, n_clusters,
                        dimnames = list(cluster_names, cluster_names))

  for (i in seq_len(n_clusters)) {
    idx_i <- cluster_indices[[i]]
    n_i <- length(idx_i)

    for (j in seq_len(n_clusters)) {
      idx_j <- cluster_indices[[j]]
      n_j <- length(idx_j)

      # Both diagonal (within-cluster) and off-diagonal (between-cluster)
      w_ij <- mat[idx_i, idx_j]
      n_possible <- n_i * n_j
      between_raw[i, j] <- net_aggregate_weights(as.vector(w_ij), method, n_possible)
    }
  }

  # Matrix path is aggregation only -- no post-processing. Caller normalises
  # downstream if they need a Markov chain (e.g. via as_tna()).
  between_weights <- between_raw

  # Compute inits from column sums
  col_sums <- colSums(between_raw)
  total <- sum(col_sums)
  if (total > 0) {
    between_inits <- col_sums / total
  } else {
    between_inits <- rep(1 / n_clusters, n_clusters)
  }
  names(between_inits) <- cluster_names

  # Build $macro
  between <- .mcml_layer(
    weights = between_weights,
    inits   = between_inits,
    labels  = cluster_names,
    data    = NULL
  )

  # ============================================================================
  # Within-cluster computation (optional)
  # ============================================================================

  within_data <- NULL
  if (isTRUE(compute_within)) {
    within_data <- lapply(seq_len(n_clusters), function(i) {
      idx_i <- cluster_indices[[i]]
      n_i <- length(idx_i)
      cl_nodes <- cluster_list[[i]]

      if (n_i <= 1) {
        # Single node: preserve self-loop value
        within_raw <- mat[idx_i, idx_i, drop = FALSE]
        dimnames(within_raw) <- list(cl_nodes, cl_nodes)
        within_weights_i <- within_raw
        within_inits_i <- setNames(1, cl_nodes)
      } else {
        # Extract within-cluster raw weights (preserves self-loops)
        within_raw <- mat[idx_i, idx_i]
        dimnames(within_raw) <- list(cl_nodes, cl_nodes)

        within_weights_i <- within_raw

        # Within-cluster inits (handle NAs)
        col_sums_w <- colSums(within_raw, na.rm = TRUE)
        total_w <- sum(col_sums_w, na.rm = TRUE)
        within_inits_i <- if (!is.na(total_w) && total_w > 0) {
          col_sums_w / total_w
        } else {
          rep(1 / n_i, n_i)
        }
        names(within_inits_i) <- cl_nodes
      }

      .mcml_layer(
        weights = within_weights_i,
        inits   = within_inits_i,
        labels  = cl_nodes,
        data    = NULL
      )
    })
    names(within_data) <- cluster_names
  }

  # ============================================================================
  # Build result
  # ============================================================================

  result <- structure(
    list(
      macro = between,
      clusters = within_data,
      cluster_members = cluster_list,
      edges = NULL,    # NULL placeholder; transitions paths fill this in
      meta = list(
        type = "aggregate",
        method = method,
        directed = directed,
        n_nodes = n,
        n_clusters = n_clusters,
        cluster_sizes = vapply(cluster_list, length, integer(1)),
        source = "matrix"
      )
    ),
    class = "mcml"
  )

  result
}

.validate_mcml_matrix <- function(mat) {
  if (any(is.na(mat) | !is.finite(mat))) {
    stop("x matrix must contain finite non-missing weights.", call. = FALSE)
  }
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  has_row_names <- !is.null(row_names)
  has_col_names <- !is.null(col_names)
  if (has_row_names && (any(is.na(row_names)) || any(!nzchar(row_names)))) {
    stop("x matrix row names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_col_names && (any(is.na(col_names)) || any(!nzchar(col_names)))) {
    stop("x matrix column names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_row_names && anyDuplicated(row_names)) {
    stop("x matrix row names must be unique.", call. = FALSE)
  }
  if (has_col_names && anyDuplicated(col_names)) {
    stop("x matrix column names must be unique.", call. = FALSE)
  }
  if (has_row_names && has_col_names && !identical(row_names, col_names)) {
    stop("x matrix row and column names must be identical and in the same order.",
         call. = FALSE)
  }
  invisible(TRUE)
}


# ==============================================================================
# 2b. Build MCML from Raw Transition Data
# ==============================================================================

#' Build MCML from Raw Transition Data
#'
#' Builds a Multi-Cluster Multi-Level (MCML) model from raw transition data
#' (edge lists or sequences) by recoding node labels to cluster labels and
#' counting actual transitions. Unlike \code{\link{cluster_summary}} which
#' aggregates a pre-computed weight matrix, this function works from the
#' original transition data to produce the TRUE Markov chain over cluster states.
#'
#' @param x Input data. Accepts multiple formats:
#'   \describe{
#'     \item{data.frame with from/to columns}{Edge list. Columns named
#'       from/source/src/v1/node1/i and to/target/tgt/v2/node2/j are
#'       auto-detected. Optional weight column (weight/w/value/strength).}
#'     \item{data.frame without from/to columns}{Sequence data. Each row is a
#'       sequence, columns are time steps. Consecutive pairs (t, t+1) become
#'       transitions.}
#'     \item{tna object}{If \code{x$data} is non-NULL, uses sequence path on
#'       the raw data. Otherwise falls back to \code{\link{cluster_summary}}.}
#'     \item{netobject}{If \code{x$data} is non-NULL, detects edge list
#'       vs sequence data. Otherwise falls back to \code{\link{cluster_summary}}.}
#'     \item{cluster_summary}{Returns as-is.}
#'     \item{square numeric matrix}{Falls back to \code{\link{cluster_summary}}.}
#'     \item{non-square or character matrix}{Treated as sequence data.}
#'   }
#'
#' @param clusters Cluster/group assignments. Accepts:
#'   \describe{
#'     \item{named list}{Direct mapping. List names = cluster names, values =
#'       character vectors of node labels.
#'       Example: \code{list(A = c("N1","N2"), B = c("N3","N4"))}}
#'     \item{data.frame}{A data frame where the first column contains node names
#'       and the second column contains group/cluster names.
#'       Example: \code{data.frame(node = c("N1","N2","N3"), group = c("A","A","B"))}}
#'     \item{membership vector}{Character or numeric vector. Node names are
#'       extracted from the data.
#'       Example: \code{c("A","A","B","B")}}
#'     \item{column name string}{For edge list data.frames, the name of a
#'       column containing cluster labels. The mapping is built from unique
#'       (node, group) pairs in both from and to columns. \strong{Limitation:}
#'       this mode assigns the row's group label to \emph{both} endpoints, so
#'       it only makes sense for edge lists where source and target nodes
#'       always share the same group (within-group edges only). For general
#'       edge lists where a single node may be source in some rows and target
#'       in others, or where source and target belong to different groups,
#'       pass an explicit named list (\code{list(G1 = c("N1","N2"), ...)})
#'       or a two-column data frame \code{data.frame(node, group)} instead.}
#'     \item{NULL}{Auto-detect from \code{netobject$nodes} or
#'       \code{$node_groups} (same logic as \code{\link{cluster_summary}}).}
#'   }
#'
#' @param method Aggregation method for combining edge weights: "sum", "mean",
#'   "median", "max", "min", "density", "geomean". Default "sum". For raw
#'   sequence/event-log inputs the function is counting observed transitions,
#'   so \code{"sum"} is the only interpretation that preserves the count
#'   semantics -- the other methods are useful when aggregating
#'   weighted edge lists or pre-existing weight matrices, where each row
#'   already represents a measurement rather than a single observation.
#' @param type Post-processing: "tna" (row-normalize), "cooccurrence"
#'   (symmetrize), "semi_markov", or "raw". Default "tna".
#' @param directed Logical. If \code{TRUE} (default), treat transitions as
#'   directed. If \code{FALSE}, symmetrize sequence- and edge-derived weights
#'   before returning raw/frequency weights or before row-normalizing
#'   transition probabilities.
#' @param compute_within Logical. Compute within-cluster matrices? Default TRUE.
#' @param actor,action,time,order,session,time_threshold Long-format event-log
#'   shortcut. When \code{action} is supplied on a data.frame input, the data
#'   is passed through \code{prepare()} to derive a wide sequence, which is
#'   then routed to the existing sequence path. Behaves identically to
#'   \code{prepare(...) |> build_network() |> build_mcml()}.
#' @param labels Optional name -> label remap applied to within-cluster nodes
#'   (the macro layer is left untouched because its labels are cluster
#'   names). Accepts a 2-column data.frame \code{(name, label)}, a named
#'   character vector \code{c(name = "label")}, or a named list. Unmapped
#'   names pass through unchanged.
#'
#' @return A \code{cluster_summary} object with \code{meta$source = "transitions"},
#'   fully compatible with \code{plot()}, \code{as_tna()}, and
#'   \code{plot()}.
#'
#' @export
#' @seealso \code{\link{cluster_summary}} for matrix-based aggregation,
#'   \code{net_as_tna()} to convert to tna objects,
#'   \code{plot()} for visualization
#'
#' @examples
#' # Edge list with clusters
#' edges <- data.frame(
#'   from = c("A", "A", "B", "C", "C", "D"),
#'   to   = c("B", "C", "A", "D", "D", "A"),
#'   weight = c(1, 2, 1, 3, 1, 2)
#' )
#' clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
#' cs <- build_mcml(edges, clusters)
#' cs$macro$weights
#'
#' # Sequence data with clusters
#' seqs <- data.frame(
#'   T1 = c("A", "C", "B"),
#'   T2 = c("B", "D", "A"),
#'   T3 = c("C", "C", "D"),
#'   T4 = c("D", "A", "C")
#' )
#' cs <- build_mcml(seqs, clusters, type = "raw")
#' cs$macro$weights
build_mcml <- function(x,
                       clusters = NULL,
                       method = c("sum", "mean", "median", "max",
                                  "min", "density", "geomean"),
                       type = c("tna", "frequency", "cooccurrence",
                                "semi_markov", "raw"),
                       directed = TRUE,
                       compute_within = TRUE,
                       actor = NULL,
                       action = NULL,
                       time = NULL,
                       order = NULL,
                       session = NULL,
                       time_threshold = 900,
                       labels = NULL) {
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop("'compute_within' must be TRUE or FALSE.", call. = FALSE)
  }

  # If already an mcml object, return as-is
  if (inherits(x, "mcml")) {
    return(x)
  }

  # Remember whether the caller passed `type` explicitly -- used below to
  # warn when type is applied to matrix input, where it has no effect.
  type_explicit <- !missing(type)

  type <- match.arg(type)
  method <- match.arg(method)

  # Long-format event-log shortcut: if `action` is supplied on a data.frame,
  # delegate to prepare() and feed the resulting wide sequence into the
  # existing sequence path. Behaves identically to
  # prepare(...) |> build_network() |> build_mcml().
  if (is.data.frame(x) && !is.null(action)) {
    prepared <- prepare(x, actor = actor, action = action, time = time,
                        order = order, session = session,
                        time_threshold = time_threshold)
    x <- prepared$sequence_data
  }

  # Coerce cograph_network so downstream branches see a netobject
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)

  input_type <- .detect_mcml_input(x)

  is_matrix_path <- input_type %in% c("matrix", "tna_matrix",
                                       "netobject_matrix")
  if (type_explicit && is_matrix_path) {
    warning(
      "`type` is ignored for matrix input: the matrix path of build_mcml() ",
      "is aggregation-only. Pass `method` to choose the aggregation and ",
      "normalise downstream if a Markov chain is needed (e.g. via as_tna()).",
      call. = FALSE
    )
  }

  result <- switch(input_type,
    "edgelist" = .build_mcml_edgelist(x, clusters, method, type,
                                       directed, compute_within),
    "sequence" = .build_mcml_sequence(x, clusters, method, type,
                                       directed, compute_within),
    "tna_data" = .build_mcml_sequence(x$data, clusters, method, type,
                                       directed, compute_within),
    "tna_matrix" = cluster_summary(x, clusters, method = method,
                                    directed = directed,
                                    compute_within = compute_within),
    "netobject_data" = {
      # Auto-detect clusters from network if not provided
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x) # nocov
      }
      # Respect the netobject's method. Count-based methods (relative,
      # frequency, co_occurrence) are well-defined to re-derive from raw
      # data, so we route through the data path for an exact conditional
      # Markov chain. Model-derived methods (attention, ising, glasso,
      # mgm, pcor, cor, mlvar_*) carry weights that re-counting would
      # silently discard, so we route through the matrix path on $weights.
      net_method <- x$method %||% "relative"
      count_methods <- c("relative", "frequency", "co_occurrence",
                         "tna", "ftna", "cna", "wcna", "wtna",
                         "wtna_cooccurrence", "wtna_transition")
      if (net_method %in% count_methods) {
        data <- x$data
        sub_type <- .detect_mcml_input(data)
        if (sub_type == "edgelist") {
          .build_mcml_edgelist(data, clusters, method, type,
                                directed, compute_within)
        } else {
          .build_mcml_sequence(data, clusters, method, type,
                                directed, compute_within)
        }
      } else {
        cluster_summary(x, clusters, method = method,
                        directed = directed,
                        compute_within = compute_within)
      }
    },
    "netobject_matrix" = {
      if (is.null(clusters)) {
        clusters <- .auto_detect_clusters(x) # nocov
      }
      cluster_summary(x, clusters, method = method,
                       directed = directed, compute_within = compute_within)
    },
    "matrix" = cluster_summary(x, clusters, method = method,
                                directed = directed,
                                compute_within = compute_within),
    stop("Cannot build MCML from input of class '", class(x)[1], "'",
         call. = FALSE)
  )

  if (!is.null(labels)) {
    result <- .apply_node_labels(result, labels)
  }
  result
}

#' Detect input type for build_mcml
#' @keywords internal
.detect_mcml_input <- function(x) {
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) return("tna_data")
    return("tna_matrix")
  }

  if (inherits(x, "netobject") || inherits(x, "netobject_ml")) {
    if (!is.null(x$data)) return("netobject_data")
    return("netobject_matrix")
  }

  if (is.data.frame(x)) {
    col_names <- tolower(names(x))
    from_cols <- c("from", "source", "src", "v1", "node1", "i")
    to_cols <- c("to", "target", "tgt", "v2", "node2", "j")
    weight_cols <- c("weight", "w", "value", "strength")
    has_from <- any(from_cols %in% col_names)
    has_to <- any(to_cols %in% col_names)
    has_weight <- any(weight_cols %in% col_names)
    if ((has_from || has_to || has_weight) && ncol(x) >= 2L) return("edgelist")
    if (ncol(x) == 2L && !.mcml_has_time_step_names(names(x))) {
      return("edgelist")
    }
    return("sequence")
  }

  if (is.matrix(x)) {
    if (is.numeric(x) && nrow(x) == ncol(x)) return("matrix")
    return("sequence")
  }

  "unknown"
}

.mcml_has_time_step_names <- function(nms) {
  if (is.null(nms) || any(is.na(nms)) || any(!nzchar(nms))) {
    return(FALSE)
  }
  all(grepl("^(t|time|step)[0-9]+$", tolower(nms)))
}

#' Auto-detect clusters from netobject
#' @keywords internal
.auto_detect_clusters <- function(x) {
  # Prefer x$nodes-embedded cluster column. x$nodes is the canonical node
  # table — if it carries a cluster column, the row order is by definition
  # aligned with the weight matrix, so positional use is safe.
  clusters <- NULL
  if (!is.null(x$nodes)) {
    cluster_cols <- c("clusters", "cluster", "groups", "group")
    for (col in cluster_cols) {
      if (col %in% names(x$nodes)) {
        clusters <- x$nodes[[col]]
        .validate_auto_clusters(clusters, "x$nodes")
        break
      }
    }
  }

  # Fall back to x$node_groups, but require explicit label alignment.
  # audit_mcml #1: node_groups was previously read positionally, which
  # silently mis-assigned nodes whenever the rows were in a different
  # order than x$nodes (a common case for externally constructed
  # netobjects). The accepted shapes are now:
  #   * data.frame with a node identifier column AND a cluster column
  #   * named character/factor vector (names = node labels, values = cluster)
  if (is.null(clusters) && !is.null(x$node_groups)) {
    ng <- x$node_groups

    target_labels <- if (!is.null(x$nodes) && "label" %in% names(x$nodes)) {
      as.character(x$nodes$label)
    } else if (!is.null(rownames(x$weights))) {
      rownames(x$weights)
    } else {
      stop("Cannot align node_groups: x$nodes$label and ",
           "rownames(x$weights) are both unavailable.", call. = FALSE)
    }

    if (is.data.frame(ng)) {
      cluster_col <- intersect(c("cluster", "group", "layer"), names(ng))
      if (length(cluster_col) == 0L) {
        stop("'node_groups' data.frame is missing a recognised cluster ",
             "column (one of: cluster, group, layer).", call. = FALSE)
      }
      node_col <- intersect(c("node", "name", "label", "id"), names(ng))
      if (length(node_col) == 0L) {
        stop("'node_groups' must include a node identifier column ",
             "(one of: node, name, label, id) so cluster assignments can ",
             "be aligned with x$nodes by label rather than by row order. ",
             "Add the column or pass clusters explicitly as a named list ",
             "or vector.", call. = FALSE)
      }
      node_values <- as.character(ng[[node_col[1L]]])
      cluster_values <- as.character(ng[[cluster_col[1L]]])
      if (any(is.na(node_values) | !nzchar(node_values))) {
        stop("'node_groups' node column must not contain missing or empty values.",
             call. = FALSE)
      }
      .validate_auto_clusters(cluster_values, "node_groups")
      duplicated_nodes <- unique(node_values[duplicated(node_values)])
      if (length(duplicated_nodes) > 0L) {
        stop("'node_groups' assigns duplicate rows to node(s): ",
             paste(utils::head(duplicated_nodes, 5L), collapse = ", "),
             call. = FALSE)
      }
      lookup <- setNames(cluster_values, node_values)
    } else if (is.atomic(ng) && !is.null(names(ng))) {
      node_values <- names(ng)
      cluster_values <- as.character(ng)
      if (any(is.na(node_values) | !nzchar(node_values))) {
        stop("'node_groups' names must not contain missing or empty values.",
             call. = FALSE)
      }
      .validate_auto_clusters(cluster_values, "node_groups")
      duplicated_nodes <- unique(node_values[duplicated(node_values)])
      if (length(duplicated_nodes) > 0L) {
        stop("'node_groups' assigns duplicate rows to node(s): ",
             paste(utils::head(duplicated_nodes, 5L), collapse = ", "),
             call. = FALSE)
      }
      lookup <- setNames(cluster_values, node_values)
    } else {
      stop("'node_groups' must be a data.frame with node + cluster ",
           "columns, or a named atomic vector keyed by node label. ",
           "Unnamed vectors cannot be safely aligned to x$nodes.",
           call. = FALSE)
    }

    missing_nodes <- setdiff(target_labels, names(lookup))
    if (length(missing_nodes) > 0L) {
      stop(sprintf(
        paste0("node_groups is missing cluster assignments for %d node(s) ",
               "present in x$nodes/weights: %s"),
        length(missing_nodes),
        paste(utils::head(missing_nodes, 5L), collapse = ", ")
      ), call. = FALSE)
    }
    clusters <- unname(lookup[target_labels])
  }

  if (is.null(clusters)) {
    stop("No clusters found in netobject. ",
         "Add a 'clusters' column to nodes or provide clusters argument.",
         call. = FALSE)
  }
  clusters
}

.validate_auto_clusters <- function(clusters, source) {
  clusters <- as.character(clusters)
  if (any(is.na(clusters) | !nzchar(clusters))) {
    stop(source, " cluster assignments must not contain missing or empty values.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Build node-to-cluster lookup from cluster specification
#' @keywords internal
.build_cluster_lookup <- function(clusters, all_nodes) {
  if (is.data.frame(clusters)) { # nocov start
    # Defensive: .normalize_clusters converts df to list before this is called
    stopifnot(ncol(clusters) >= 2)
    nodes <- as.character(clusters[[1]])
    groups <- as.character(clusters[[2]])
    clusters <- split(nodes, groups)
  } # nocov end

  if (is.list(clusters) && !is.data.frame(clusters)) {
    # Named list: cluster_name -> node vector
    .validate_cluster_partition(clusters, all_nodes)
    lookup <- character(0)
    for (cl_name in names(clusters)) {
      nodes <- clusters[[cl_name]]
      lookup[nodes] <- cl_name
    }
    # Verify all nodes are mapped
    unmapped <- setdiff(all_nodes, names(lookup))
    if (length(unmapped) > 0) {
      stop("Unmapped nodes: ",
           paste(utils::head(unmapped, 5), collapse = ", "),
           call. = FALSE)
    }
    return(lookup)
  }

  if (is.character(clusters) || is.factor(clusters)) {
    clusters <- as.character(clusters)
    if (length(clusters) != length(all_nodes)) {
      stop("Membership vector length (", length(clusters),
           ") must equal number of unique nodes (", length(all_nodes), ")",
           call. = FALSE)
    }
    lookup <- setNames(clusters, all_nodes)
    return(lookup)
  }

  if (is.numeric(clusters) || is.integer(clusters)) {
    if (length(clusters) != length(all_nodes)) {
      stop("Membership vector length (", length(clusters),
           ") must equal number of unique nodes (", length(all_nodes), ")",
           call. = FALSE)
    }
    lookup <- setNames(as.character(clusters), all_nodes)
    return(lookup)
  }

  stop("clusters must be a named list, character/numeric vector, or column name",
       call. = FALSE)
}

#' Empirical first-state distribution from a wide sequence data.frame.
#' Matches `tna::tna()`'s inits: tabulate the values of column 1 (the first
#' time-step), dropping NAs, optionally restricted to a `restrict_to`
#' whitelist (anything outside the whitelist is treated as NA so it is
#' dropped). Normalised to sum to 1. Returns a zero-vector with
#' `length(levels)` entries when no rows contribute.
#' @keywords internal
#' @noRd
.first_state_distribution <- function(wide_df, levels, restrict_to = NULL) {
  empty <- setNames(rep(0, length(levels)), levels)
  if (!is.data.frame(wide_df) || nrow(wide_df) == 0L || ncol(wide_df) == 0L) {
    return(empty)
  }
  first <- as.character(wide_df[[1]])
  if (!is.null(restrict_to)) first[!first %in% restrict_to] <- NA_character_
  first <- first[!is.na(first)]
  if (length(first) == 0L) return(empty)
  fr <- table(factor(first, levels = levels))
  out <- as.numeric(fr) / sum(fr)
  names(out) <- levels
  out
}

#' Build cluster_summary from transition vectors
#' @keywords internal
.build_from_transitions <- function(from_nodes, to_nodes, weights,
                                     cluster_lookup, cluster_list,
                                     method, type, directed,
                                     compute_within, data = NULL) {
  if (!is.character(from_nodes) || !is.character(to_nodes)) {
    stop("'from_nodes' and 'to_nodes' must be character vectors.",
         call. = FALSE)
  }
  if (length(from_nodes) != length(to_nodes) ||
      length(from_nodes) != length(weights)) {
    stop("'from_nodes', 'to_nodes', and 'weights' must have the same length.",
         call. = FALSE)
  }
  if (any(is.na(from_nodes) | !nzchar(from_nodes) |
          is.na(to_nodes) | !nzchar(to_nodes))) {
    stop("'from_nodes' and 'to_nodes' must not contain missing or empty values.",
         call. = FALSE)
  }
  if (!is.numeric(weights) || any(is.na(weights) | !is.finite(weights))) {
    stop("'weights' must be a finite non-missing numeric vector.",
         call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("'weights' must not contain negative values.", call. = FALSE)
  }
  if (!is.list(cluster_list) || is.data.frame(cluster_list)) {
    stop("'cluster_list' must be a named list.", call. = FALSE)
  }
  cluster_nodes <- sort(unique(unlist(cluster_list, use.names = FALSE)))
  .validate_cluster_partition(cluster_list, cluster_nodes)
  if (!is.character(cluster_lookup) || is.null(names(cluster_lookup))) {
    stop("'cluster_lookup' must be a named character vector.", call. = FALSE)
  }
  if (any(is.na(names(cluster_lookup)) | !nzchar(names(cluster_lookup))) ||
      any(is.na(cluster_lookup) | !nzchar(cluster_lookup))) {
    stop("'cluster_lookup' names and values must not be missing or empty.",
         call. = FALSE)
  }
  if (anyDuplicated(names(cluster_lookup))) {
    stop("'cluster_lookup' names must be unique.", call. = FALSE)
  }
  missing_lookup <- setdiff(unique(c(from_nodes, to_nodes)),
                            names(cluster_lookup))
  if (length(missing_lookup) > 0L) {
    stop("'cluster_lookup' is missing node(s): ",
         paste(utils::head(missing_lookup, 5L), collapse = ", "),
         call. = FALSE)
  }
  unknown_clusters <- setdiff(unique(unname(cluster_lookup)),
                              names(cluster_list))
  if (length(unknown_clusters) > 0L) {
    stop("'cluster_lookup' contains unknown cluster(s): ",
         paste(utils::head(unknown_clusters, 5L), collapse = ", "),
         call. = FALSE)
  }
  if (!is.logical(compute_within) || length(compute_within) != 1L ||
      is.na(compute_within)) {
    stop("'compute_within' must be TRUE or FALSE.", call. = FALSE)
  }

  # Sort clusters alphabetically (TNA convention)
  cluster_list <- cluster_list[order(names(cluster_list))]
  cluster_names <- names(cluster_list)
  n_clusters <- length(cluster_names)

  # Recode to cluster labels
  from_clusters <- cluster_lookup[from_nodes]
  to_clusters <- cluster_lookup[to_nodes]

  # ---- Between-cluster matrix (includes diagonal = within-cluster loops) ----
  between_raw <- matrix(0, n_clusters, n_clusters,
                        dimnames = list(cluster_names, cluster_names))

  # Include ALL transitions -- node-level self-loops (A->A) are valid
  # cluster-level self-loops (e.g. discuss->discuss = Social->Social)
  b_from <- from_clusters
  b_to <- to_clusters
  b_w <- weights

  if (length(b_from) > 0) {
    cluster_sizes <- vapply(cluster_list, length, integer(1L))
    n_possible_mat <- if (method == "density") {
      outer(cluster_sizes, cluster_sizes, "*")
    } else NULL

    fc <- factor(b_from, levels = cluster_names)
    tc <- factor(b_to, levels = cluster_names)
    agg <- tapply(b_w, list(fc, tc), function(w) {
      net_aggregate_weights(w, method)
    }, default = NA_real_)
    agg[is.na(agg)] <- 0
    between_raw[] <- as.numeric(agg)
    if (!is.null(n_possible_mat)) {
      sums <- tapply(b_w, list(fc, tc), sum, default = 0)
      sums[is.na(sums)] <- 0
      between_raw[] <- as.numeric(sums) / as.numeric(n_possible_mat)
      between_raw[!is.finite(between_raw)] <- 0
    }
  }

  # Process between weights
  between_weights <- .process_weights(between_raw, type, directed)

  # Compute initial-state probabilities. When sequence data is available we
  # use the empirical first-state distribution (matches tna::tna()'s inits);
  # otherwise (edgelist or no data) we fall back to the column-sum proxy
  # documented earlier -- there is no first-state to read off an edgelist.
  is_edgelist_data <- identical(attr(data, "source"), "edgelist")
  is_seq_for_inits <- is.data.frame(data) && !is_edgelist_data &&
    !any(tolower(names(data)) %in%
    c("from", "source", "src", "v1", "node1", "i",
      "to", "target", "tgt", "v2", "node2", "j"))
  between_inits <- if (is_seq_for_inits) {
    # tna-equivalent macro inits: tabulate column 1 (first time-step) of the
    # wide sequence, drop NAs, map each first state to its cluster, normalise.
    first_state <- as.character(data[[1]])
    first_state <- first_state[!is.na(first_state)]
    if (length(first_state) == 0L) {
      setNames(rep(0, n_clusters), cluster_names)
    } else {
      first_clu <- cluster_lookup[first_state]
      fr <- table(factor(first_clu, levels = cluster_names))
      out <- as.numeric(fr) / sum(fr)
      names(out) <- cluster_names
      out
    }
  } else {
    col_sums <- colSums(between_raw)
    total <- sum(col_sums)
    out <- if (total > 0) col_sums / total else rep(1 / n_clusters, n_clusters)
    names(out) <- cluster_names
    out
  }

  # ---- Build recoded sequence data for bootstrap compatibility ----
  # The bootstrap fast path (.bootstrap_transition) needs a data.frame
  # that it can resample row-wise. Two cases:
  #
  #   * Wide-format sequence input (rows = actors, cols = time-steps):
  #     recode every cell to the cluster label (macro) or keep only
  #     nodes of that cluster (per-cluster).
  #
  #   * Edgelist input (each row is a single from -> to transition):
  #     build 2-column pseudo-sequences, each row representing one
  #     transition. Resampling rows = resampling transitions, which
  #     is the correct bootstrap at the edge level.
  between_seq_data <- NULL
  within_seq_data_list <- NULL
  is_seq <- is.data.frame(data) && !is_edgelist_data &&
    !any(tolower(names(data)) %in%
    c("from", "source", "src", "v1", "node1", "i",
      "to", "target", "tgt", "v2", "node2", "j"))

  if (is_seq) {
    # Recode sequences to cluster labels
    between_seq_data <- as.data.frame(
      lapply(data, function(col) {
        recoded <- unname(cluster_lookup[as.character(col)])
        recoded
      }),
      stringsAsFactors = FALSE
    )

    # Filter sequences per cluster (keep only that cluster's nodes, NA others)
    within_seq_data_list <- lapply(cluster_list, function(cl_nodes) {
      as.data.frame(
        lapply(data, function(col) {
          vals <- as.character(col)
          vals[!vals %in% cl_nodes] <- NA_character_
          vals
        }),
        stringsAsFactors = FALSE
      )
    })
  } else if (length(from_nodes) > 0L) {
    # Edgelist branch: every transition -> one 2-row pseudo-sequence.
    # Columns named `from`/`to` so the data is self-describing on inspection.
    # Tagged with source = "edgelist" so bootstrap_network() can warn --
    # edgelist bootstrap treats transitions as independent, ignoring
    # within-actor correlation, so CIs are anti-conservative.
    between_seq_data <- data.frame(
      from = unname(from_clusters),
      to = unname(to_clusters),
      stringsAsFactors = FALSE
    )
    attr(between_seq_data, "source") <- "edgelist"

    within_seq_data_list <- lapply(cluster_list, function(cl_nodes) {
      keep <- from_nodes %in% cl_nodes & to_nodes %in% cl_nodes
      df <- data.frame(
        from = from_nodes[keep],
        to = to_nodes[keep],
        stringsAsFactors = FALSE
      )
      attr(df, "source") <- "edgelist"
      df
    })
  }

  between <- .mcml_layer(
    weights = between_weights,
    inits   = between_inits,
    labels  = cluster_names,
    data    = between_seq_data
  )

  # ---- Within-cluster matrices ----
  within_data <- NULL
  if (isTRUE(compute_within)) {
    # Filter within-cluster transitions (same cluster, including self-loops)
    is_within <- from_clusters == to_clusters
    w_from <- from_nodes[is_within]
    w_to <- to_nodes[is_within]
    w_w <- weights[is_within]

    within_data <- lapply(seq_len(n_clusters), function(i) {
      cl_name <- cluster_names[i]
      cl_nodes <- cluster_list[[cl_name]]
      n_i <- length(cl_nodes)

      if (n_i <= 1) {
        # Singleton cluster: aggregate self-loop weight, then route through
        # .process_weights so type = "tna"/"semi_markov" row-normalises just
        # like the multi-node branch (1x1 [N] -> [1.0]; 1x1 [0] stays [0]).
        keep_self <- w_from %in% cl_nodes & w_to %in% cl_nodes
        self_weight <- if (any(keep_self)) {
          net_aggregate_weights(w_w[keep_self], method)
        } else { 0 }
        within_raw <- matrix(self_weight, 1, 1,
                             dimnames = list(cl_nodes, cl_nodes))
        within_weights_i <- .process_weights(within_raw, type, directed)
        within_inits_i <- if (is_seq_for_inits) {
          .first_state_distribution(data, levels = cl_nodes, restrict_to = cl_nodes)
        } else {
          setNames(1, cl_nodes)
        }
      } else {
        # Filter transitions for this cluster (keep self-loops)
        keep <- w_from %in% cl_nodes & w_to %in% cl_nodes
        cf <- w_from[keep]
        ct <- w_to[keep]
        cw <- w_w[keep]

        within_raw <- matrix(0, n_i, n_i,
                              dimnames = list(cl_nodes, cl_nodes))

        if (length(cf) > 0) {
          fc <- factor(cf, levels = cl_nodes)
          tc <- factor(ct, levels = cl_nodes)
          agg <- tapply(cw, list(fc, tc), function(w) {
            net_aggregate_weights(w, method)
          }, default = NA_real_)
          agg[is.na(agg)] <- 0
          within_raw[] <- as.numeric(agg)
        }

        within_weights_i <- .process_weights(within_raw, type, directed)

        # Within-cluster inits: first cluster-k state per row, NA-masked
        # against other clusters (matches tna::tna() on the masked sequence).
        # For edgelist or matrix-only input, fall back to the column-sum proxy.
        within_inits_i <- if (is_seq_for_inits) {
          .first_state_distribution(data, levels = cl_nodes,
                                 restrict_to = cl_nodes)
        } else {
          col_sums_w <- colSums(within_raw, na.rm = TRUE)
          total_w <- sum(col_sums_w, na.rm = TRUE)
          out <- if (!is.na(total_w) && total_w > 0) col_sums_w / total_w
                 else rep(1 / n_i, n_i)
          names(out) <- cl_nodes
          out
        }
      }

      # Attach filtered sequence data for this cluster
      cl_seq_data <- if (!is.null(within_seq_data_list)) {
        within_seq_data_list[[cl_name]]
      } else {
        NULL
      }

      .mcml_layer(
        weights = within_weights_i,
        inits   = within_inits_i,
        labels  = cl_nodes,
        data    = cl_seq_data
      )
    })
    names(within_data) <- cluster_names
  }

  # ---- Edges data.frame ----
  edge_type <- ifelse(from_clusters == to_clusters, "within", "between")
  edges <- data.frame(
    from = from_nodes,
    to = to_nodes,
    weight = weights,
    cluster_from = unname(from_clusters),
    cluster_to = unname(to_clusters),
    type = edge_type,
    stringsAsFactors = FALSE
  )

  # ---- Assemble result ----
  all_nodes <- sort(unique(unlist(cluster_list, use.names = FALSE)))
  n_nodes <- length(all_nodes)

  structure(
    list(
      macro = between,
      clusters = within_data,
      cluster_members = cluster_list,
      edges = edges,
      meta = list(
        type = type,
        method = method,
        directed = directed,
        n_nodes = n_nodes,
        n_clusters = n_clusters,
        cluster_sizes = vapply(cluster_list, length, integer(1)),
        source = "transitions"
      )
    ),
    class = "mcml"
  )
}

#' Build MCML from edge list data.frame
#' @keywords internal
.build_mcml_edgelist <- function(df, clusters, method, type,
                                  directed, compute_within) {

  col_names <- tolower(names(df))

  # Detect from/to columns
  from_col <- which(col_names %in% c("from", "source", "src",
                                       "v1", "node1", "i"))[1]
  to_col <- which(col_names %in% c("to", "target", "tgt",
                                     "v2", "node2", "j"))[1]
  if (is.na(from_col)) from_col <- 1L
  if (is.na(to_col)) {
    to_col <- setdiff(seq_len(ncol(df)), from_col)[1L]
  }
  if (is.na(to_col) || from_col == to_col) {
    stop("Edge-list input must include at least two endpoint columns.",
         call. = FALSE)
  }

  # Detect weight column
  weight_col <- which(col_names %in% c("weight", "w", "value", "strength"))[1]
  has_weight <- !is.na(weight_col)

  from_vals <- as.character(df[[from_col]])
  to_vals <- as.character(df[[to_col]])
  if (any(is.na(from_vals) | !nzchar(from_vals) |
          is.na(to_vals) | !nzchar(to_vals))) {
    stop("Edge-list source and target columns must not contain missing or empty values.",
         call. = FALSE)
  }

  if (has_weight) {
    weights <- df[[weight_col]]
    if (!is.numeric(weights)) {
      stop("Edge-list weight column must be numeric.", call. = FALSE)
    }
    if (any(is.na(weights) | !is.finite(weights))) {
      stop("Edge-list weight column must contain finite non-missing values.",
           call. = FALSE)
    }
    if (any(weights < 0)) {
      stop("Edge-list weight column must not contain negative values.",
           call. = FALSE)
    }
  } else {
    weights <- rep(1, nrow(df))
  }

  all_nodes <- sort(unique(c(from_vals, to_vals)))

  # Handle clusters parameter
  if (is.character(clusters) && length(clusters) == 1 &&
      clusters %in% names(df)) {
    # Column name: build lookup from both from+group and to+group
    group_col <- as.character(df[[clusters]])
    if (any(is.na(group_col) | !nzchar(group_col))) {
      stop("Edge-list cluster column must not contain missing or empty values.",
           call. = FALSE)
    }

    node_group <- data.frame(
      node = c(from_vals, to_vals),
      group = c(group_col, group_col),
      stringsAsFactors = FALSE
    )
    node_group <- unique(node_group)
    conflicts <- unique(node_group$node[duplicated(node_group$node)])
    if (length(conflicts) > 0L) {
      stop("Edge-list cluster column assigns nodes to multiple groups: ",
           paste(utils::head(conflicts, 5L), collapse = ", "),
           ". Pass clusters as a named list or node-group data.frame instead.",
           call. = FALSE)
    }

    full_map <- setNames(node_group$group, node_group$node)

    # Build cluster_list
    cluster_list <- split(names(full_map), unname(full_map))
    cluster_list <- lapply(cluster_list, sort)
    cluster_lookup <- full_map

    # Re-derive all_nodes from the lookup
    all_nodes <- sort(names(cluster_lookup))
  } else {
    # List or vector clusters
    if (is.null(clusters)) {
      stop("clusters argument is required for data.frame input", call. = FALSE)
    }

    if (is.list(clusters) && !is.data.frame(clusters)) {
      cluster_list <- clusters
    } else {
      # Membership vector
      cluster_list <- .normalize_clusters(clusters, all_nodes)
    }

    cluster_lookup <- .build_cluster_lookup(cluster_list, all_nodes)
  }

  attr(df, "source") <- "edgelist"

  .build_from_transitions(from_vals, to_vals, weights,
                            cluster_lookup, cluster_list,
                            method, type, directed, compute_within,
                            data = df)
}

#' Build MCML from sequence data.frame
#' @keywords internal
.build_mcml_sequence <- function(df, clusters, method, type,
                                  directed, compute_within) {

  if (is.matrix(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)

  stopifnot(is.data.frame(df))

  nc <- ncol(df)
  if (nc < 2) {
    stop("Sequence data must have at least 2 columns (time steps)",
         call. = FALSE)
  }

  # Extract consecutive pairs: (col[t], col[t+1]) for all rows
  pairs <- lapply(seq_len(nc - 1), function(t) {
    from_t <- as.character(df[[t]])
    to_t <- as.character(df[[t + 1]])
    data.frame(from = from_t, to = to_t, stringsAsFactors = FALSE)
  })
  pairs <- as.data.frame(data.table::rbindlist(pairs))

  # Remove NA pairs
  valid <- !is.na(pairs$from) & !is.na(pairs$to)
  from_vals <- pairs$from[valid]
  to_vals <- pairs$to[valid]
  weights <- rep(1, length(from_vals))

  observed_nodes <- as.character(unlist(df, use.names = FALSE))
  observed_nodes <- observed_nodes[!is.na(observed_nodes)]
  all_nodes <- sort(unique(observed_nodes))

  if (is.null(clusters)) {
    stop("clusters argument is required for sequence data", call. = FALSE)
  }

  if (is.list(clusters) && !is.data.frame(clusters)) {
    cluster_list <- clusters
  } else {
    cluster_list <- .normalize_clusters(clusters, all_nodes)
  }

  cluster_lookup <- .build_cluster_lookup(cluster_list, all_nodes)

  .build_from_transitions(from_vals, to_vals, weights,
                            cluster_lookup, cluster_list,
                            method, type, directed, compute_within,
                            data = df)
}

#' Process weights based on type
#' @keywords internal
.process_weights <- function(raw_weights, type, directed = TRUE) {
  if (!is.matrix(raw_weights) || !is.numeric(raw_weights)) {
    stop("'raw_weights' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(raw_weights) != ncol(raw_weights)) {
    stop("'raw_weights' must be a square matrix.", call. = FALSE)
  }
  if (!is.character(type) || length(type) != 1L || is.na(type)) {
    stop("'type' must be a single non-missing character value.", call. = FALSE)
  }
  valid_types <- c("raw", "frequency", "cooccurrence", "tna", "semi_markov")
  if (!type %in% valid_types) {
    stop("'type' must be one of: ",
         paste(valid_types, collapse = ", "), call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }

  if (!isTRUE(directed) || type == "cooccurrence") {
    raw_weights <- (raw_weights + t(raw_weights)) / 2
  }

  if (type == "raw" || type == "frequency") {
    return(raw_weights)
  }

  if (type == "cooccurrence") {
    return(raw_weights)
  }

  if (type == "tna" || type == "semi_markov") {
    # Row-normalize so rows sum to 1
    rs <- rowSums(raw_weights, na.rm = TRUE)
    processed <- raw_weights / ifelse(rs == 0 | is.na(rs), 1, rs)
    processed[is.na(processed)] <- 0
    return(processed)
  }

  raw_weights # nocov
}

#' Convert cluster_summary to tna Objects
#'
#' Converts a \code{cluster_summary} object to proper tna objects that can be
#' used with all functions from the tna package. Creates both a between-cluster
#' tna model (cluster-level transitions) and within-cluster tna models (internal
#' transitions within each cluster).
#'
#' This is the final step in the MCML workflow, enabling full integration with
#' the tna package for centrality analysis, bootstrap validation, permutation
#' tests, and visualization.
#'
#' @param x A \code{cluster_summary} object created by
#'   \code{\link{cluster_summary}}. The aggregated weights are passed to
#'   \code{tna::tna()}, which row-normalises them as needed.
#'
#' @return A \code{cluster_tna} object (S3 class) containing:
#'   \describe{
#'     \item{between}{A tna object representing cluster-level transitions.
#'       Contains \code{$weights} (k x k transition matrix), \code{$inits}
#'       (initial distribution), and \code{$labels} (cluster names).
#'       Use this for analyzing how learners/entities move between high-level
#'       groups or phases.}
#'     \item{within}{Named list of tna objects, one per cluster. Each tna object
#'       represents internal transitions within that cluster. Contains
#'       \code{$weights} (n_i x n_i matrix), \code{$inits} (initial distribution),
#'       and \code{$labels} (node labels). Clusters with single nodes or zero-row
#'       nodes are excluded (tna requires positive row sums).}
#'   }
#'
#' @details
#' ## Requirements
#'
#' The tna package must be installed. If not available, the function throws
#' an error with installation instructions.
#'
#' ## Workflow
#'
#' \preformatted{
#' # Full MCML workflow
#' net <- build_network(data, method = "relative")
#' net$nodes$clusters <- group_assignments
#' cs <- cluster_summary(net)
#' tna_models <- as_tna(cs)
#'
#' # Now use tna package functions
#' plot(tna_models$macro)
#' tna::centralities(tna_models$macro)
#' tna::bootstrap(tna_models$macro, iter = 1000)
#'
#' # Analyze within-cluster patterns
#' plot(tna_models$clusters$ClusterA)
#' tna::centralities(tna_models$clusters$ClusterA)
#' }
#'
#' ## Excluded Clusters
#'
#' A within-cluster network is dropped only when row-normalisation would
#' fail. Specifically, when the recorded \code{net_method} is
#' \code{"relative"} (row-stochastic transitions) and any node in the
#' cluster has zero outgoing weight, that cluster is excluded from
#' \code{$clusters} and a \code{warning()} is emitted listing the dropped
#' cluster names. For \code{net_method = "frequency"} (raw counts), a
#' zero-row node is a legitimate sink and the cluster is retained.
#' The macro / between-cluster network always includes every cluster
#' regardless of the per-cluster drop decisions.
#'
#' If a cluster you expect to see is missing from the returned
#' \code{$clusters}, check the warning output and consider building with
#' \code{type = "raw"} (which carries through to a frequency-method
#' netobject and skips the drop) or inspect \code{rowSums(x$clusters[[cl]]$weights)}.
#'
#' @export
#' @seealso
#'   \code{\link{cluster_summary}} to create the input object,
#'   \code{plot()} for visualization without conversion,
#'   \code{tna::tna} for the underlying tna constructor
#'
#' @examples
#' mat <- matrix(runif(36), 6, 6)
#' rownames(mat) <- colnames(mat) <- LETTERS[1:6]
#' clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"), G3 = c("E", "F"))
#' cs <- cluster_summary(mat, clusters)
#' tna_models <- as_tna(cs)
#' tna_models
#' tna_models$macro$weights
as_tna <- function(x) {
  UseMethod("as_tna")
}

#' @rdname as_tna
#' @return A \code{netobject_group} with data preserved from each sub-network.
#' @export
as_tna.mcml <- function(x) {
  # Determine the method to record on the wrapped netobjects:
  #   * "raw"/"frequency"/"aggregate" -> "frequency" (matrix-path or counts)
  #   * "tna"/"semi_markov" etc.      -> "relative" (already row-normalised)
  meta_type <- x$meta$type
  net_method <- if (is.null(meta_type) ||
                    meta_type %in% c("raw", "frequency", "aggregate")) {
    "frequency"
  } else {
    "relative"
  }
  directed <- if (!is.null(x$meta$directed)) x$meta$directed else TRUE

  # Macro
  macro_net <- .wrap_netobject(x$macro$weights, data = x$macro$data,
                               method = net_method, directed = directed,
                               inits = x$macro$inits)

  # Per-cluster. Drop a cluster only when the wrapped netobject would be
  # un-normalisable (relative-method requires positive row sums). For
  # frequency-method counts a zero row is a legitimate sink and is kept.
  cluster_nets <- if (!is.null(x$clusters)) {
    drop_zero_rows <- net_method == "relative"
    dropped <- character()
    nets <- lapply(names(x$clusters), function(cl) {
      w <- x$clusters[[cl]]$weights
      d <- x$clusters[[cl]]$data
      i <- x$clusters[[cl]]$inits
      if (drop_zero_rows && any(rowSums(w) == 0)) {
        dropped[[length(dropped) + 1L]] <<- cl
        return(NULL)
      }
      .wrap_netobject(w, data = d, method = net_method, directed = directed,
                      inits = i)
    })
    names(nets) <- names(x$clusters)
    if (length(dropped) > 0L) {
      warning("Dropped clusters with zero row sums (cannot row-normalise): ",
              paste(dropped, collapse = ", "), call. = FALSE)
    }
    nets[!vapply(nets, is.null, logical(1))]
  } else {
    list()
  }

  result <- c(list(macro = macro_net), cluster_nets)
  class(result) <- "netobject_group"
  result
}


#' Wrap a weight matrix + optional data into a minimal netobject
#' @noRd
.wrap_netobject <- function(mat, data = NULL, method = "relative",
                            directed = TRUE, inits = NULL) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("'mat' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(mat)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("'method' must be a single non-missing character value.",
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  states <- rownames(mat)
  if (is.null(states)) {
    states <- as.character(seq_len(nrow(mat)))
    dimnames(mat) <- list(states, states)
  }
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  nodes_df <- data.frame(
    id = seq_along(states), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(states) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per state.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(states, names(inits))
      extra_inits <- setdiff(names(inits), states)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match matrix state names.", call. = FALSE)
      }
      inits <- inits[states]
    } else {
      names(inits) <- states
    }
  }

  structure(
    list(
      data       = data,
      weights    = mat,
      inits      = inits,
      nodes      = nodes_df,
      edges      = edges,
      directed   = directed,
      method     = method,
      params     = list(),
      scaling    = NULL,
      threshold  = 0,
      n_nodes    = length(states),
      n_edges    = nrow(edges),
      level      = NULL,
      meta       = list(source = "nestimate", layout = NULL,
                        tna = list(method = method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}

#' @rdname as_tna
#' @return A \code{tna} object constructed from the input.
#' @export
as_tna.default <- function(x) {
 if (inherits(x, "tna")) {
    return(x)
  }
  stop("Cannot convert object of class '", class(x)[1], "' to tna", call. = FALSE)
}

#' Normalize cluster specification to list format
#' @keywords internal
#' @noRd
.normalize_clusters <- function(clusters, node_names) {
  if (is.data.frame(clusters)) {
    # Data frame with node and group columns
    if (ncol(clusters) < 2L) {
      stop("clusters data.frame must have at least two columns.",
           call. = FALSE)
    }
    nodes <- as.character(clusters[[1]])
    groups <- as.character(clusters[[2]])
    if (any(is.na(nodes) | !nzchar(nodes))) {
      stop("clusters data.frame node column must not contain missing or empty values.",
           call. = FALSE)
    }
    if (any(is.na(groups) | !nzchar(groups))) {
      stop("clusters data.frame group column must not contain missing or empty values.",
           call. = FALSE)
    }
    clusters <- split(nodes, groups)
  }

  if (is.list(clusters)) {
    # Already a list - validate node names
    .validate_cluster_partition(clusters, node_names)
    return(clusters)
  }

  if (is.vector(clusters) && (is.numeric(clusters) || is.integer(clusters))) {
    # Membership vector
    clusters <- .align_cluster_membership(clusters, node_names)
    if (any(is.na(clusters) | !is.finite(clusters))) {
      stop("Membership vector must not contain missing or non-finite values.",
           call. = FALSE)
    }
    # Convert to list
    unique_clusters <- sort(unique(clusters))
    cluster_list <- lapply(unique_clusters, function(k) {
      node_names[clusters == k]
    })
    names(cluster_list) <- as.character(unique_clusters)
    return(cluster_list)
  }

  if (is.factor(clusters) || is.character(clusters)) {
    # Named membership
    clusters <- .align_cluster_membership(clusters, node_names)
    clusters <- as.character(clusters)
    if (any(is.na(clusters) | !nzchar(clusters))) {
      stop("Membership vector must not contain missing or empty values.",
           call. = FALSE)
    }
    unique_clusters <- unique(clusters)
    cluster_list <- lapply(unique_clusters, function(k) {
      node_names[clusters == k]
    })
    names(cluster_list) <- unique_clusters
    return(cluster_list)
  }

  stop("clusters must be a list, numeric vector, or factor", call. = FALSE)
}

.align_cluster_membership <- function(clusters, node_names) {
  if (length(clusters) != length(node_names)) {
    stop("Membership vector length (", length(clusters),
         ") must equal number of nodes (", length(node_names), ")",
         call. = FALSE)
  }
  nm <- names(clusters)
  if (is.null(nm)) {
    return(clusters)
  }
  if (any(is.na(nm) | !nzchar(nm))) {
    stop("Named membership vector names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop("Named membership vector names must be unique.", call. = FALSE)
  }
  missing_nodes <- setdiff(node_names, nm)
  extra_nodes <- setdiff(nm, node_names)
  if (length(missing_nodes) > 0L || length(extra_nodes) > 0L) {
    stop("Named membership vector names must match node names exactly.",
         call. = FALSE)
  }
  clusters[node_names]
}

.validate_cluster_partition <- function(clusters, node_names) {
  if (!is.list(clusters) || is.data.frame(clusters)) {
    stop("clusters must be a named list.", call. = FALSE)
  }
  cluster_names <- names(clusters)
  if (is.null(cluster_names) || any(!nzchar(cluster_names)) ||
      any(is.na(cluster_names))) {
    stop("clusters list must have non-empty cluster names.", call. = FALSE)
  }
  cluster_sizes <- vapply(clusters, length, integer(1L))
  if (any(cluster_sizes == 0L)) {
    stop("clusters list contains empty clusters: ",
         paste(names(clusters)[cluster_sizes == 0L], collapse = ", "),
         call. = FALSE)
  }

  nodes <- unlist(clusters, use.names = FALSE)
  nodes <- as.character(nodes)
  if (any(is.na(nodes)) || any(!nzchar(nodes))) {
    stop("clusters list contains missing or empty node names.",
         call. = FALSE)
  }

  unknown <- setdiff(nodes, node_names)
  if (length(unknown) > 0L) {
    stop("Unknown nodes in clusters: ",
         paste(utils::head(unknown, 5L), collapse = ", "), call. = FALSE)
  }

  duplicated_nodes <- unique(nodes[duplicated(nodes)])
  if (length(duplicated_nodes) > 0L) {
    stop("Nodes assigned to multiple clusters: ",
         paste(utils::head(duplicated_nodes, 5L), collapse = ", "),
         call. = FALSE)
  }

  unmapped <- setdiff(node_names, nodes)
  if (length(unmapped) > 0L) {
    stop("Unmapped nodes: ",
         paste(utils::head(unmapped, 5L), collapse = ", "),
         call. = FALSE)
  }

  invisible(TRUE)
}

# ==============================================================================
# S3 Methods
# ==============================================================================

#' Print Method for mcml
#'
#' @param x An \code{mcml} object.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
#' clusters <- list(G1 = c("A","B"), G2 = c("C"))
#' cs <- build_mcml(seqs, clusters)
#' print(cs)
#' \donttest{
#' seqs <- data.frame(
#'   T1 = c("A","B","A"), T2 = c("B","C","B"),
#'   T3 = c("C","A","C"), T4 = c("A","B","A")
#' )
#' clusters <- c("Alpha", "Beta", "Alpha")
#' cs <- build_mcml(seqs, clusters, type = "raw")
#' print(cs)
#' }
#'
#' @export
print.mcml <- function(x, ...) {
  .mcml_check_unused_dots("print.mcml", ...)
  n_clusters <- x$meta$n_clusters
  n_nodes <- x$meta$n_nodes
  cluster_sizes <- x$meta$cluster_sizes

  cat("MCML Network\n")
  cat("============\n")
  if (!is.null(x$meta$type)) {
    cat("Type:", x$meta$type, " | Method:", x$meta$method, "\n")
  } else {
    cat("Method:", x$meta$method, " (matrix path, aggregation-only)\n")
  }
  cat("Nodes:", n_nodes, " | Clusters:", n_clusters, "\n")
  if (!is.null(x$edges)) {
    cat("Transitions:", nrow(x$edges), "\n")
    n_between <- sum(x$edges$type == "between")
    n_within <- sum(x$edges$type == "within")
    cat("  Macro:", n_between, " | Per-cluster:", n_within, "\n\n")
  } else {
    cat("\n")
  }

  cat("Clusters:\n")
  for (cl in names(x$cluster_members)) {
    cat("  ", cl, " (", cluster_sizes[cl], "): ",
        paste(x$cluster_members[[cl]], collapse = ", "), "\n", sep = "")
  }

  cat("\nMacro (cluster-level) weights:\n")
  print(round(x$macro$weights, 4))

  invisible(x)
}

#' Summary Method for mcml
#'
#' @param object An \code{mcml} object.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#'
#' @return A tidy data frame with one row per cluster and columns
#'   \code{cluster}, \code{size}, \code{within_total}, \code{between_out},
#'   \code{between_in}. For undirected macro networks the in/out split is
#'   not meaningful, so \code{between_out} reports total incident weight
#'   and \code{between_in} is \code{NA}. The data frame is returned
#'   silently \emph{without} printing the full object -- call
#'   \code{print(object)} explicitly if you want the verbose dump.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
#' clusters <- list(G1 = c("A","B"), G2 = c("C"))
#' cs <- build_mcml(seqs, clusters)
#' summary(cs)
#' \donttest{
#' seqs <- data.frame(
#'   T1 = c("A","B","A"), T2 = c("B","C","B"),
#'   T3 = c("C","A","C"), T4 = c("A","B","A")
#' )
#' clusters <- c("Alpha", "Beta", "Alpha")
#' cs <- build_mcml(seqs, clusters, type = "raw")
#' summary(cs)
#' }
#'
#' @export
summary.mcml <- function(object, ...) {
  .mcml_check_unused_dots("summary.mcml", ...)
  as_mat <- function(x) {
    if (is.null(x)) return(NULL)
    if (is.matrix(x) && is.numeric(x)) return(x)
    if (is.list(x) && is.numeric(x$weights)) return(x$weights)
    if (is.list(x) && is.numeric(x$matrix))  return(x$matrix)
    NULL
  }

  macro    <- as_mat(object$macro)
  cluster_mats <- lapply(object$clusters, as_mat)
  cluster_names <- names(cluster_mats)
  if (is.null(cluster_names))
    cluster_names <- paste0("Cluster_", seq_along(cluster_mats))
  sizes <- object$meta$cluster_sizes
  if (is.null(sizes)) sizes <- vapply(object$cluster_members, length, integer(1L))

  within_total <- vapply(cluster_mats, function(mat) {
    if (is.null(mat)) return(0)
    sum(abs(mat[row(mat) != col(mat)]))
  }, numeric(1L))

  directed <- isTRUE(object$meta$directed)
  if (!is.null(macro)) {
    diag0 <- macro
    diag(diag0) <- 0
    if (directed) {
      between_out <- rowSums(abs(diag0))
      between_in  <- colSums(abs(diag0))
    } else {
      # For undirected macro networks, in/out split is not meaningful;
      # report total incident weight under between_out and leave between_in NA.
      between_out <- rowSums(abs(diag0))
      between_in  <- rep(NA_real_, length(cluster_names))
    }
  } else {
    between_out <- rep(NA_real_, length(cluster_names))
    between_in  <- rep(NA_real_, length(cluster_names))
  }

  data.frame(
    cluster      = cluster_names,
    size         = as.integer(sizes),
    within_total = as.numeric(within_total),
    between_out  = as.numeric(between_out),
    between_in   = as.numeric(between_in),
    stringsAsFactors = FALSE,
    row.names    = NULL
  )
}

.mcml_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(
    method, "() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", "),
    call. = FALSE
  )
}
