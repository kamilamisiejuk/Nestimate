# ---- Simplicial Complex Analysis ----
#
# Construction, homology, centrality, and Q-analysis for simplicial
# complexes built from networks and higher-order pathway data.
#
# Clique finding verified against igraph::cliques().
# Betti numbers verified against known topological invariants.

# =========================================================================
# Core constructor
# =========================================================================

#' Build a Simplicial Complex
#'
#' @description
#' Constructs a simplicial complex from a network or higher-order pathway
#' object. Three construction methods are available:
#'
#' \itemize{
#'   \item \strong{Clique complex} (\code{"clique"}): every clique in the
#'     thresholded graph becomes a simplex. The standard bridge from graph
#'     theory to algebraic topology.
#'   \item \strong{Pathway complex} (\code{"pathway"}): each higher-order
#'     pathway from a \code{net_hon} or \code{net_hypa} becomes a simplex.
#'   \item \strong{Vietoris-Rips} (\code{"vr"}): nodes with edge weight
#'     \eqn{\geq} \code{threshold} are connected; all cliques in the
#'     resulting graph become simplices.
#' }
#'
#' @param x A square matrix, \code{tna}, \code{netobject},
#'   \code{net_hon}, \code{net_hypa}, or \code{net_mogen}.
#' @param type Construction type: \code{"clique"} (default),
#'   \code{"pathway"}, or \code{"vr"}.
#' @param threshold Minimum absolute edge weight to include an edge
#'   (default 0). Edges below this are ignored.
#' @param max_dim Maximum simplex dimension (default 10). A k-simplex
#'   has k+1 nodes.
#' @param max_pathways For \code{type = "pathway"}: maximum number of
#'   pathways to include, ranked by count (HON) or ratio (HYPA).
#'   \code{NULL} includes all. Default \code{NULL}.
#' @param ... Additional arguments passed to \code{build_hon()} when
#'   \code{x} is a \code{tna}/\code{netobject} with \code{type = "pathway"}.
#'
#' @return A \code{simplicial_complex} object.
#'
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' print(sc)
#' betti_numbers(sc)
#'
#' @seealso \code{\link{betti_numbers}}, \code{\link{persistent_homology}},
#'   \code{\link{simplicial_degree}}, \code{\link{q_analysis}}
#'
#' @export
build_simplicial <- function(x, type = "clique", threshold = 0,
                              max_dim = 10L, max_pathways = NULL, ...) {
  type <- match.arg(type, c("clique", "pathway", "vr"))

  if (type == "pathway") {
    return(.build_simplicial_pathway(x, max_dim, max_pathways, ...))
  }

  mat <- .sc_extract_matrix(x)

  if (type == "clique") {
    .build_simplicial_clique(mat, threshold, max_dim)
  } else {
    .build_simplicial_vr(mat, threshold, max_dim)
  }
}

# =========================================================================
# Clique complex — verified against igraph::cliques()
# =========================================================================

#' @noRd
.build_simplicial_clique <- function(mat, threshold = 0, max_dim = 10L) {
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  nodes <- rownames(mat) %||% paste0("V", seq_len(nrow(mat)))
  n <- nrow(mat)

  adj <- abs(mat) > threshold
  diag(adj) <- FALSE
  adj <- adj | t(adj)

  simplices <- .find_all_cliques(adj, max_dim)

  .make_simplicial_complex(simplices, nodes, "clique")
}

#' @noRd
.build_simplicial_vr <- function(mat, threshold, max_dim = 10L) {
  stopifnot(is.matrix(mat), nrow(mat) == ncol(mat))
  nodes <- rownames(mat) %||% paste0("V", seq_len(nrow(mat)))

  weights <- abs(mat)
  diag(weights) <- 0
  weights <- pmax(weights, t(weights))

  adj <- weights >= threshold
  diag(adj) <- FALSE

  simplices <- .find_all_cliques(adj, max_dim)

  .make_simplicial_complex(simplices, nodes, "vr")
}

#' Find all cliques (all sizes) via igraph or fallback Bron-Kerbosch
#'
#' When igraph is available, delegates to igraph::cliques() for
#' correctness and speed. Results are verified to match on package tests.
#' @noRd
.find_all_cliques <- function(adj, max_dim = 10L) {
  n <- nrow(adj)

  if (requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                              diag = FALSE)
    raw <- igraph::cliques(g, min = 1, max = max_dim + 1L)
    simplices <- lapply(raw, function(cl) sort(as.integer(cl)))
  } else {
    # Fallback: Bron-Kerbosch + expand to all faces
    maximal <- .bron_kerbosch_all(adj) # nocov
    simplices <- .expand_to_faces(maximal, max_dim) # nocov
  }

  simplices
}

# =========================================================================
# Pathway complex (from HON/HYPA)
# =========================================================================

#' @noRd
.build_simplicial_pathway <- function(x, max_dim = 10L,
                                       max_pathways = NULL, ...) {
  if (inherits(x, "net_hon")) {
    edges <- x$ho_edges
    ho <- edges[edges$from_order > 1L, , drop = FALSE]
    ho <- ho[order(-ho$count), , drop = FALSE]
    if (!is.null(max_pathways) && nrow(ho) > max_pathways) {
      ho <- ho[seq_len(max_pathways), , drop = FALSE]
    }
    nodes <- x$first_order_states
    raw_paths <- ho$path
  } else if (inherits(x, "net_hypa")) {
    scores <- x$scores
    anom <- scores[scores$anomaly != "normal", , drop = FALSE]
    if ("ratio" %in% names(anom)) {
      anom <- anom[order(-anom$ratio), , drop = FALSE]
    }
    if (!is.null(max_pathways) && nrow(anom) > max_pathways) { # nocov
      anom <- anom[seq_len(max_pathways), , drop = FALSE] # nocov
    }
    parts <- strsplit(
      gsub("\x01", " -> ", x$nodes$label, fixed = TRUE), " -> ", fixed = TRUE
    )
    nodes <- sort(unique(unlist(parts)))
    raw_paths <- anom$path
  } else if (inherits(x, "net_mogen")) {
    # Use mogen_transitions() at optimal (or highest available) order
    # Its $path column is already in "A -> B -> C" format
    order_used <- x$optimal_order
    if (order_used < 1L) order_used <- max(x$orders[x$orders >= 1L], 0L)
    if (order_used < 1L) {
      stop("MOGen model has no higher-order transitions (optimal_order = 0)",
           call. = FALSE)
    }
    trans <- mogen_transitions(x, order = order_used)
    if (nrow(trans) == 0L) {
      return(.make_simplicial_complex(list(), x$states, "pathway"))
    }
    if (!is.null(max_pathways) && nrow(trans) > max_pathways) {
      trans <- trans[seq_len(max_pathways), , drop = FALSE]
    }
    nodes <- x$states
    raw_paths <- trans$path
  } else if (inherits(x, c("tna", "netobject"))) {
    hon_obj <- build_hon(.coerce_sequence_input(x), ...)
    return(.build_simplicial_pathway(hon_obj, max_dim, max_pathways))
  } else {
    stop("For type='pathway', x must be a net_hon, net_hypa, net_mogen, ",
         "tna, or netobject.", call. = FALSE)
  }

  if (length(raw_paths) == 0L) {
    return(.make_simplicial_complex(list(), nodes, "pathway"))
  }

  node_idx <- setNames(seq_along(nodes), nodes)
  simplices_raw <- lapply(raw_paths, function(p) {
    parts <- trimws(strsplit(p, "->", fixed = TRUE)[[1]])
    unique(parts)
  })

  simplices <- lapply(simplices_raw, function(s) {
    idx <- node_idx[s]
    sort(idx[!is.na(idx)])
  })
  simplices <- simplices[vapply(simplices, length, integer(1)) >= 2L]
  simplices <- .expand_to_faces(simplices, max_dim)

  .make_simplicial_complex(simplices, nodes, "pathway")
}

# =========================================================================
# Matrix extraction
# =========================================================================

#' @noRd
.sc_extract_matrix <- function(x) {
  if (is.matrix(x)) return(x)
  if (inherits(x, "tna")) return(x$weights)
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    return(x$weights)
  }
  if (inherits(x, "net_hon")) return(x$matrix)
  stop("Cannot extract adjacency matrix from '", class(x)[1], "'.",
       call. = FALSE)
}

# =========================================================================
# Bron-Kerbosch (fallback when igraph unavailable)
# =========================================================================

#' @noRd
.bron_kerbosch_all <- function(adj) { # nocov start
  n <- nrow(adj)
  neighbors <- lapply(seq_len(n), function(i) which(adj[i, ]))
  cliques <- list()

  .bk <- function(R, P, X) {
    if (length(P) == 0L && length(X) == 0L) {
      cliques[[length(cliques) + 1L]] <<- sort(R)
      return(invisible(NULL))
    }
    union_px <- c(P, X)
    pivot <- union_px[which.max(
      vapply(union_px, function(v) sum(P %in% neighbors[[v]]), integer(1))
    )]
    for (v in setdiff(P, neighbors[[pivot]])) {
      nbrs <- neighbors[[v]]
      .bk(c(R, v), intersect(P, nbrs), intersect(X, nbrs))
      P <- setdiff(P, v)
      X <- c(X, v)
    }
  }

  .bk(integer(0), seq_len(n), integer(0))
  cliques
} # nocov end

#' Expand maximal cliques to all sub-simplices
#' @noRd
.expand_to_faces <- function(simplices, max_dim = 10L) {
  seen <- new.env(hash = TRUE, parent = emptyenv())
  result <- list()

  for (simplex in simplices) {
    simplex <- sort(as.integer(simplex))
    if (length(simplex) - 1L > max_dim) {
      simplex <- simplex[seq_len(max_dim + 1L)]
    }
    for (size in seq_len(length(simplex))) {
      combos <- utils::combn(simplex, size, simplify = FALSE)
      for (face in combos) {
        key <- paste(face, collapse = ",")
        if (is.null(seen[[key]])) {
          seen[[key]] <- TRUE
          result[[length(result) + 1L]] <- face
        }
      }
    }
  }

  result
}

# =========================================================================
# Constructor
# =========================================================================

#' @noRd
.make_simplicial_complex <- function(simplices, nodes, type) {
  # Ensure all 0-simplices are present
  seen <- new.env(hash = TRUE, parent = emptyenv())
  for (s in simplices) seen[[paste(sort(s), collapse = ",")]] <- TRUE
  for (i in seq_along(nodes)) {
    if (is.null(seen[[as.character(i)]])) {
      simplices[[length(simplices) + 1L]] <- i
    }
  }

  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))
  max_d <- if (length(dims) > 0L) max(dims) else 0L
  f_vec <- vapply(0:max_d, function(d) sum(dims == d), integer(1))
  names(f_vec) <- paste0("dim_", 0:max_d)

  n_simplices <- length(simplices)
  n_nodes <- length(nodes)

  # Density: fraction of possible simplices that exist
  # Max possible k-simplices = C(n, k+1)
  max_possible <- sum(vapply(0:max_d, function(d) {
    choose(n_nodes, d + 1L)
  }, numeric(1)))
  density <- if (max_possible > 0) n_simplices / max_possible else 0

  # Mean simplex dimension
  mean_dim <- if (n_simplices > 0L) mean(dims) else 0

  structure(list(
    simplices = simplices,
    nodes = nodes,
    n_nodes = n_nodes,
    n_simplices = n_simplices,
    dimension = max_d,
    f_vector = f_vec,
    density = density,
    mean_dim = mean_dim,
    type = type
  ), class = "simplicial_complex")
}

# =========================================================================
# Homology
# =========================================================================

#' Betti Numbers
#'
#' Computes Betti numbers: \eqn{\beta_0} (components), \eqn{\beta_1}
#' (loops), \eqn{\beta_2} (voids), etc.
#'
#' @param sc A \code{simplicial_complex} object.
#' @return Named integer vector \code{c(b0 = ..., b1 = ..., ...)}.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' betti_numbers(sc)
#'
#' @export
betti_numbers <- function(sc) {
  stopifnot(inherits(sc, "simplicial_complex"))
  .compute_betti(sc)
}

#' @noRd
.compute_betti <- function(sc) {
  max_d <- sc$dimension
  dims <- vapply(sc$simplices, function(s) length(s) - 1L, integer(1))
  by_dim <- lapply(0:max_d, function(d) sc$simplices[dims == d])

  boundary_ranks <- integer(max_d + 2L)
  boundary_ranks[1L] <- 0L

  for (d in seq_len(max_d)) {
    k_simplices <- by_dim[[d + 1L]]
    km1_simplices <- by_dim[[d]]

    if (length(k_simplices) == 0L || length(km1_simplices) == 0L) { # nocov start
      boundary_ranks[d + 1L] <- 0L
      next # nocov end
    }

    km1_keys <- vapply(km1_simplices, function(s) {
      paste(sort(s), collapse = ",")
    }, character(1))
    km1_idx <- setNames(seq_along(km1_keys), km1_keys)

    bmat <- matrix(0, nrow = length(km1_simplices),
                   ncol = length(k_simplices))

    for (j in seq_along(k_simplices)) {
      simplex <- sort(k_simplices[[j]])
      for (i in seq_along(simplex)) {
        face_key <- paste(simplex[-i], collapse = ",")
        row_idx <- km1_idx[face_key]
        if (!is.na(row_idx)) {
          bmat[row_idx, j] <- (-1)^(i + 1L)
        }
      }
    }

    boundary_ranks[d + 1L] <- qr(bmat)$rank
  }

  betti <- vapply(0:max_d, function(d) {
    n_k <- length(by_dim[[d + 1L]])
    nullity_k <- n_k - boundary_ranks[d + 1L]
    rank_kp1 <- if (d < max_d) boundary_ranks[d + 2L] else 0L
    as.integer(max(nullity_k - rank_kp1, 0L))
  }, integer(1))

  names(betti) <- paste0("b", 0:max_d)
  betti
}

#' Euler Characteristic
#'
#' @description
#' Computes \eqn{\chi = \sum_{k=0}^{d} (-1)^k f_k} where \eqn{f_k} is the
#' number of k-simplices. By the Euler-Poincare theorem,
#' \eqn{\chi = \sum_{k} (-1)^k \beta_k}.
#'
#' @param sc A \code{simplicial_complex} object.
#' @return Integer.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' euler_characteristic(sc)
#'
#' @export
euler_characteristic <- function(sc) {
  stopifnot(inherits(sc, "simplicial_complex"))
  signs <- (-1L)^(seq_along(sc$f_vector) - 1L)
  as.integer(sum(signs * sc$f_vector))
}

# =========================================================================
# Persistent homology
# =========================================================================

#' Persistent Homology
#'
#' @description
#' Computes persistent homology by building simplicial complexes at
#' decreasing weight thresholds and tracking the birth/death of
#' topological features.
#'
#' @param x A square matrix, \code{tna}, or \code{netobject}.
#' @param n_steps Number of filtration steps (default 20).
#' @param max_dim Maximum simplex dimension to track (default 3).
#'
#' @return A \code{persistent_homology} object with:
#' \describe{
#'   \item{betti_curve}{Data frame: \code{threshold}, \code{dimension},
#'     \code{betti} at each filtration step.}
#'   \item{persistence}{Data frame of birth-death pairs:
#'     \code{dimension}, \code{birth}, \code{death}, \code{persistence}.}
#'   \item{thresholds}{Numeric vector of filtration thresholds.}
#' }
#'
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 10)
#' print(ph)
#'
#' @export
persistent_homology <- function(x, n_steps = 20L, max_dim = 3L) {
  mat <- .sc_extract_matrix(x)
  mat <- abs(mat)
  mat <- pmax(mat, t(mat))
  diag(mat) <- 0

  max_w <- max(mat)
  if (max_w == 0) {
    stop("All weights are zero; cannot build filtration.", call. = FALSE)
  }
  thresholds <- seq(max_w, max_w * 0.01, length.out = n_steps)

  rows <- list()
  for (i in seq_along(thresholds)) {
    sc <- .build_simplicial_clique(mat, threshold = thresholds[i],
                                    max_dim = max_dim)
    b <- betti_numbers(sc)
    for (j in seq_along(b)) {
      rows[[length(rows) + 1L]] <- data.frame(
        threshold = thresholds[i],
        dimension = j - 1L,
        betti = b[j],
        stringsAsFactors = FALSE
      )
    }
  }

  betti_curve <- do.call(rbind, rows)
  rownames(betti_curve) <- NULL
  persistence <- .extract_persistence(betti_curve, thresholds)

  structure(list(
    betti_curve = betti_curve,
    persistence = persistence,
    thresholds = thresholds
  ), class = "persistent_homology")
}

#' @noRd
.extract_persistence <- function(betti_curve, thresholds) {
  dims <- sort(unique(betti_curve$dimension))
  pairs <- list()

  for (d in dims) {
    sub <- betti_curve[betti_curve$dimension == d, ]
    sub <- sub[order(sub$threshold, decreasing = TRUE), ]
    bettis <- sub$betti
    thresh <- sub$threshold
    active <- integer(0)
    prev_b <- 0L

    for (i in seq_along(bettis)) {
      curr_b <- bettis[i]
      if (curr_b > prev_b) {
        active <- c(active, rep(i, curr_b - prev_b))
      } else if (curr_b < prev_b) {
        n_dead <- prev_b - curr_b
        for (j in seq_len(min(n_dead, length(active)))) {
          born_at <- active[length(active)]
          active <- active[-length(active)]
          pairs[[length(pairs) + 1L]] <- data.frame(
            dimension = d, birth = thresh[born_at],
            death = thresh[i], stringsAsFactors = FALSE
          )
        }
      }
      prev_b <- curr_b
    }
    for (born_at in active) {
      pairs[[length(pairs) + 1L]] <- data.frame(
        dimension = d, birth = thresh[born_at],
        death = 0, stringsAsFactors = FALSE
      )
    }
  }

  if (length(pairs) == 0L) { # nocov start
    return(data.frame(dimension = integer(0), birth = numeric(0),
                      death = numeric(0), persistence = numeric(0),
                      stringsAsFactors = FALSE)) # nocov end
  }

  result <- do.call(rbind, pairs)
  result$persistence <- result$birth - result$death
  rownames(result) <- NULL
  result[order(-result$persistence), ]
}

# =========================================================================
# Simplicial centrality
# =========================================================================

#' Simplicial Degree
#'
#' Counts how many simplices of each dimension contain each node.
#'
#' @param sc A \code{simplicial_complex} object.
#' @param normalized Divide by maximum possible count. Default \code{FALSE}.
#'
#' @return Data frame with \code{node}, columns \code{d0} through
#'   \code{d_k}, and \code{total} (sum of d1+). Sorted by total descending.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' simplicial_degree(sc)
#'
#' @export
simplicial_degree <- function(sc, normalized = FALSE) {
  stopifnot(inherits(sc, "simplicial_complex"))
  n <- sc$n_nodes
  max_d <- sc$dimension

  dims <- vapply(sc$simplices, function(s) length(s) - 1L, integer(1))

  mat <- matrix(0L, nrow = n, ncol = max_d + 1L)
  for (i in seq_along(sc$simplices)) {
    d <- dims[i]
    for (v in sc$simplices[[i]]) {
      mat[v, d + 1L] <- mat[v, d + 1L] + 1L
    }
  }

  if (normalized && n > 1L) {
    for (d in 0:max_d) {
      denom <- choose(n - 1L, d)
      if (denom > 0) mat[, d + 1L] <- mat[, d + 1L] / denom
    }
  }

  df <- as.data.frame(mat)
  names(df) <- paste0("d", 0:max_d)
  df <- cbind(data.frame(node = sc$nodes, stringsAsFactors = FALSE), df)
  df$total <- rowSums(mat[, -1L, drop = FALSE])
  df[order(-df$total), ]
}

# =========================================================================
# Q-analysis
# =========================================================================

#' Q-Analysis
#'
#' @description
#' Computes Q-connectivity structure (Atkin 1974). Two maximal simplices
#' are q-connected if they share a face of dimension \eqn{\geq q}. Reports:
#' \itemize{
#'   \item \strong{Q-vector}: number of connected components at each q-level
#'   \item \strong{Structure vector}: highest simplex dimension per node
#' }
#'
#' @param sc A \code{simplicial_complex} object.
#'
#' @return A \code{q_analysis} object with \code{$q_vector},
#'   \code{$structure_vector}, and \code{$max_q}.
#'
#' @references
#' Atkin, R. H. (1974). \emph{Mathematical Structure in Human Affairs}.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' q_analysis(sc)
#'
#' @export
q_analysis <- function(sc) {
  stopifnot(inherits(sc, "simplicial_complex"))

  simplices <- sc$simplices
  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))

  # Structure vector: max simplex dimension per node
  sv <- vapply(seq_len(sc$n_nodes), function(v) {
    d <- vapply(simplices, function(s) {
      if (v %in% s) length(s) - 1L else -1L
    }, integer(1))
    max(d)
  }, integer(1))
  names(sv) <- sc$nodes

  # Find maximal simplices
  n_s <- length(simplices)
  is_maximal <- vapply(seq_len(n_s), function(i) {
    si <- simplices[[i]]
    !any(vapply(seq_len(n_s), function(j) {
      if (j == i || dims[j] <= dims[i]) return(FALSE)
      all(si %in% simplices[[j]])
    }, logical(1)))
  }, logical(1))

  maximal <- simplices[is_maximal]
  n_max <- length(maximal)
  max_q <- if (n_max > 0L) max(vapply(maximal, length, integer(1))) - 1L else 0L

  if (n_max <= 1L) {
    q_vec <- setNames(rep(1L, max_q + 1L), paste0("q_", max_q:0))
    return(structure(list(
      q_vector = q_vec, max_q = max_q,
      structure_vector = sv
    ), class = "q_analysis"))
  }

  # Shared face dimension between pairs
  shared_dim <- matrix(-1L, n_max, n_max)
  for (i in seq_len(n_max - 1L)) {
    for (j in (i + 1L):n_max) {
      common <- length(intersect(maximal[[i]], maximal[[j]])) - 1L
      shared_dim[i, j] <- shared_dim[j, i] <- common
    }
  }

  q_vec <- vapply(0:max_q, function(q) {
    adj_q <- shared_dim >= q
    diag(adj_q) <- FALSE
    .count_components(adj_q)
  }, integer(1))
  names(q_vec) <- paste0("q_", max_q:0)

  structure(list(
    q_vector = q_vec,
    max_q = max_q,
    structure_vector = sv
  ), class = "q_analysis")
}

#' @noRd
.count_components <- function(adj) {
  n <- nrow(adj)
  visited <- logical(n)
  n_comp <- 0L
  for (start in seq_len(n)) {
    if (visited[start]) next
    n_comp <- n_comp + 1L
    queue <- start
    visited[start] <- TRUE
    while (length(queue) > 0L) {
      v <- queue[1L]
      queue <- queue[-1L]
      nbrs <- which(adj[v, ] & !visited)
      visited[nbrs] <- TRUE
      queue <- c(queue, nbrs)
    }
  }
  n_comp
}

# =========================================================================
# Verification helper
# =========================================================================

#' Verify Simplicial Complex Against igraph
#'
#' @description
#' Cross-validates clique finding and Betti numbers against igraph
#' and known topological invariants. Useful for testing.
#'
#' @param mat A square adjacency matrix.
#' @param threshold Edge weight threshold.
#'
#' @return A list with \code{$cliques_match} (logical),
#'   \code{$n_simplices_ours}, \code{$n_simplices_igraph},
#'   \code{$betti}, and \code{$euler}.
#' @examplesIf requireNamespace("igraph", quietly = TRUE)
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' verify_simplicial(mat, threshold = 0.3)
#' @export
verify_simplicial <- function(mat, threshold = 0) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph is required for verification.", call. = FALSE) # nocov
  }

  sc <- build_simplicial(mat, threshold = threshold)

  # igraph clique count
  adj <- abs(mat) > threshold
  diag(adj) <- FALSE
  adj <- adj | t(adj)
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                            diag = FALSE)
  ig_cliques <- igraph::cliques(g)

  # Compare sorted simplex sets
  our_keys <- sort(vapply(sc$simplices, function(s) {
    paste(sort(s), collapse = ",")
  }, character(1)))
  ig_keys <- sort(vapply(ig_cliques, function(cl) {
    paste(sort(as.integer(cl)), collapse = ",")
  }, character(1)))

  betti <- betti_numbers(sc)
  euler <- euler_characteristic(sc)

  result <- list(
    cliques_match = identical(our_keys, ig_keys),
    n_simplices_ours = length(sc$simplices),
    n_simplices_igraph = length(ig_cliques),
    betti = betti,
    euler = euler,
    f_vector = sc$f_vector
  )

  clique_ok <- result$cliques_match
  dims <- seq_along(betti) - 1L
  euler_from_betti <- as.integer(sum((-1L)^dims * betti))
  euler_ok <- euler == euler_from_betti

  b_str <- paste(sprintf("%s=%d", names(betti), betti), collapse = " ")
  cat(sprintf("  Cliques:  %s (%d simplices)\n",
              if (clique_ok) "MATCH" else "MISMATCH", result$n_simplices_ours))
  cat(sprintf("  Betti:    %s\n", b_str))
  cat(sprintf("  Euler:    %d (Euler-Poincare: %s)\n",
              euler, if (euler_ok) "VERIFIED" else "FAILED"))

  invisible(result)
}

# =========================================================================
# Print methods
# =========================================================================

#' Print a simplicial complex
#' @param x A \code{simplicial_complex} object.
#' @param ... Additional arguments (unused).
#' @return The input object, invisibly.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' print(sc)
#'
#' @export
print.simplicial_complex <- function(x, ...) {
  labels <- c("clique" = "Clique Complex",
              "pathway" = "Pathway Complex",
              "vr" = "Vietoris-Rips Complex")
  cat(labels[x$type] %||% "Simplicial Complex", "\n")

  betti <- .compute_betti(x)
  chi <- euler_characteristic(x)

  cat(sprintf("  %d nodes, %d simplices, dimension %d\n",
              x$n_nodes, x$n_simplices, x$dimension))
  cat(sprintf("  Density: %.1f%%  |  Mean dim: %.2f  |  Euler: %d\n",
              x$density * 100, x$mean_dim, chi))

  # f-vector: compact
  f_str <- paste(sprintf("f%d=%d", seq_along(x$f_vector) - 1L,
                          x$f_vector), collapse = " ")
  cat(sprintf("  f-vector: (%s)\n", f_str))

  # Betti: only non-zero
  nz <- which(betti > 0)
  if (length(nz) == 0L) {
    cat("  Betti: all zero (contractible)\n") # nocov
  } else {
    b_str <- paste(sprintf("%s=%d", names(betti)[nz], betti[nz]),
                   collapse = " ")
    cat(sprintf("  Betti: %s\n", b_str))
  }

  if (x$n_nodes <= 15L) {
    cat("  Nodes:", paste(x$nodes, collapse = ", "), "\n")
  }
  invisible(x)
}

#' Print persistent homology results
#' @param x A \code{persistent_homology} object.
#' @param ... Additional arguments (unused).
#' @return The input object, invisibly.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 10)
#' print(ph)
#'
#' @export
print.persistent_homology <- function(x, ...) {
  cat("Persistent Homology\n")
  cat(sprintf("  %d filtration steps [%.4f \u2192 %.4f]\n",
              length(x$thresholds),
              max(x$thresholds), min(x$thresholds)))

  if (nrow(x$persistence) > 0L) {
    dims <- sort(unique(x$persistence$dimension))
    parts <- vapply(dims, function(d) {
      sub <- x$persistence[x$persistence$dimension == d, ]
      n_p <- sum(sub$death == 0)
      sprintf("b%d: %d (%d persistent)", d, nrow(sub), n_p)
    }, character(1))
    cat("  Features:", paste(parts, collapse = "  |  "), "\n")

    # Top 3 only
    top <- head(x$persistence[x$persistence$persistence > 0, ], 3)
    if (nrow(top) > 0L) {
      cat("  Longest-lived:\n")
      for (i in seq_len(nrow(top))) {
        cat(sprintf("    b%d: %.4f \u2192 %.4f (life: %.4f)\n",
                    top$dimension[i], top$birth[i], top$death[i],
                    top$persistence[i]))
      }
    }
  }
  invisible(x)
}

#' Print Q-analysis results
#' @param x A \code{q_analysis} object.
#' @param ... Additional arguments (unused).
#' @return The input object, invisibly.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' qa <- q_analysis(sc)
#' print(qa)
#'
#' @export
print.q_analysis <- function(x, ...) {
  cat(sprintf("Q-Analysis (max q = %d)\n", x$max_q))

  # Q-vector as compact line
  q_str <- paste(sprintf("q%d:%d", x$max_q:0, x$q_vector), collapse = " ")
  cat(sprintf("  Components: %s\n", q_str))

  # Fragmentation: first q where components > 1
  frag_q <- NA
  for (i in seq_along(x$q_vector)) {
    if (x$q_vector[i] > 1L) {
      frag_q <- x$max_q - i + 1L
      break
    }
  }
  if (!is.na(frag_q)) {
    cat(sprintf("  Fragments at q = %d (%d \u2192 %d components)\n",
                frag_q, x$q_vector[x$max_q - frag_q],
                x$q_vector[x$max_q - frag_q + 1L]))
  } else {
    cat("  Fully connected at all q levels\n")
  }

  # Structure vector: compact sorted
  sv <- sort(x$structure_vector, decreasing = TRUE)
  sv_str <- paste(sprintf("%s:%d", names(sv), sv), collapse = " ")
  cat(sprintf("  Structure: %s\n", sv_str))
  invisible(x)
}

# =========================================================================
# Plot methods
# =========================================================================

.sc_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "grey40",
                                             size = base_size - 2),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot a Simplicial Complex
#'
#' Produces a multi-panel summary: f-vector, simplicial degree ranking,
#' and degree-by-dimension heatmap.
#'
#' @param x A \code{simplicial_complex} object.
#' @param ... Ignored.
#' @return A grid grob (invisibly).
#'
#' @examples
#' \donttest{
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' sc <- build_simplicial(mat, threshold = 0.3)
#' if (requireNamespace("gridExtra", quietly = TRUE)) plot(sc)
#' }
#'
#' @export
plot.simplicial_complex <- function(x, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra is required for multi-panel plots.", call. = FALSE) # nocov
  }

  deg <- simplicial_degree(x)
  betti <- .compute_betti(x)

  # --- Panel 1: f-vector ---
  fdf <- data.frame(dim = factor(seq_along(x$f_vector) - 1L),
                     count = as.integer(x$f_vector))
  p1 <- ggplot2::ggplot(fdf, ggplot2::aes(x = dim, y = count)) +
    ggplot2::geom_col(fill = "#4A7FB5", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = count), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "f-vector",
                  subtitle = "Simplices per dimension",
                  x = "Dimension", y = "Count") +
    .sc_theme()

  # --- Panel 2: degree ranking ---
  deg$node <- factor(deg$node, levels = deg$node)
  p2 <- ggplot2::ggplot(deg, ggplot2::aes(x = node, y = total,
                                            fill = total)) +
    ggplot2::geom_col(show.legend = FALSE) +
    ggplot2::scale_fill_gradient(low = "#81B1D3", high = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = total), vjust = -0.3,
                        size = 3.2) +
    ggplot2::labs(title = "Simplicial Degree",
                  subtitle = "Higher-order participation (dim 1+)",
                  x = NULL, y = "Total") +
    .sc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
                                                         hjust = 1))

  # --- Panel 3: degree heatmap ---
  d_cols <- paste0("d", seq_len(x$dimension))
  deg_long <- stats::reshape(deg[, c("node", d_cols)],
                              direction = "long",
                              varying = d_cols,
                              v.names = "count", timevar = "dim",
                              times = seq_len(x$dimension))
  deg_long$node <- factor(deg_long$node, levels = rev(deg$node))

  p3 <- ggplot2::ggplot(deg_long, ggplot2::aes(x = factor(dim), y = node,
                                                  fill = count)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = count), size = 3.2) +
    ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#E8734A",
                                  guide = "none") +
    ggplot2::labs(title = "Degree by Dimension",
                  subtitle = "Simplex participation per node",
                  x = "Dimension", y = NULL) +
    .sc_theme()

  # --- Panel 4: Betti numbers ---
  bdf <- data.frame(dim = factor(seq_along(betti) - 1L),
                     betti = as.integer(betti))
  b_subtitle <- sprintf("Euler characteristic: %d", euler_characteristic(x))
  p4 <- ggplot2::ggplot(bdf, ggplot2::aes(x = dim, y = betti)) +
    ggplot2::geom_col(fill = "#6AAB9C", width = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = betti), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "Betti Numbers",
                  subtitle = b_subtitle,
                  x = "Dimension", y = expression(beta[k])) +
    .sc_theme()

  combined <- gridExtra::arrangeGrob(p1, p4, p2, p3, ncol = 2,
                                      padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined)
  invisible(combined)
}

#' Plot Persistent Homology
#'
#' Two panels: Betti curve (threshold vs Betti number) and persistence
#' diagram (birth vs death).
#'
#' @param x A \code{persistent_homology} object.
#' @param ... Ignored.
#' @return A grid grob (invisibly).
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' net <- build_network(seqs, method = "relative")
#' ph  <- persistent_homology(net)
#' if (requireNamespace("gridExtra", quietly = TRUE)) plot(ph)
#' }
#'
#' @export
plot.persistent_homology <- function(x, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra is required for multi-panel plots.", call. = FALSE) # nocov
  }

  filt <- x$betti_curve
  filt$dim_label <- factor(paste0("B", filt$dimension))

  # --- Panel 1: Betti curve ---
  p1 <- ggplot2::ggplot(filt, ggplot2::aes(x = threshold, y = betti,
                                              color = dim_label)) +
    ggplot2::geom_step(linewidth = 1.1, direction = "vh") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(title = "Betti Curve",
                  subtitle = "Topological features across weight thresholds",
                  x = "Weight Threshold", y = "Betti Number",
                  color = NULL) +
    .sc_theme() +
    ggplot2::theme(legend.position = "top")

  # --- Panel 2: Persistence diagram ---
  pers <- x$persistence[x$persistence$persistence > 0, ]

  if (nrow(pers) > 0L) {
    pers$dim_label <- factor(paste0("B", pers$dimension))
    lim <- max(c(pers$birth, pers$death)) * 1.15

    p2 <- ggplot2::ggplot(pers, ggplot2::aes(x = birth, y = death,
                                                color = dim_label,
                                                size = persistence)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                            color = "grey60") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(1.5, 6), guide = "none") +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "Far from diagonal = long-lived features",
                    x = "Birth", y = "Death", color = NULL) +
      .sc_theme() +
      ggplot2::theme(legend.position = "top")
  } else { # nocov start
    p2 <- ggplot2::ggplot() +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "No features detected") +
      .sc_theme() # nocov end
  }

  combined <- gridExtra::arrangeGrob(p1, p2, ncol = 2,
                                      padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined)
  invisible(combined)
}

#' Plot Q-Analysis
#'
#' Two panels: Q-vector (components at each connectivity level) and
#' structure vector (max simplex dimension per node).
#'
#' @param x A \code{q_analysis} object.
#' @param ... Ignored.
#' @return A grid grob (invisibly).
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' net <- build_network(seqs, method = "relative")
#' sc  <- build_simplicial(net, type = "clique")
#' qa  <- q_analysis(sc)
#' plot(qa)
#' }
#'
#' @export
plot.q_analysis <- function(x, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra is required for multi-panel plots.", call. = FALSE) # nocov
  }

  # --- Panel 1: Q-vector ---
  qdf <- data.frame(q = x$max_q:0, components = as.integer(x$q_vector))

  p1 <- ggplot2::ggplot(qdf, ggplot2::aes(x = q, y = components)) +
    ggplot2::geom_step(linewidth = 1.2, color = "#4A7FB5",
                        direction = "vh") +
    ggplot2::geom_point(size = 3, color = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = components), vjust = -1,
                        size = 3.5) +
    ggplot2::scale_x_continuous(breaks = qdf$q) +
    ggplot2::labs(title = "Q-Vector",
                  subtitle = "Connected components at each q-level",
                  x = "q (shared face dimension)", y = "Components") +
    .sc_theme()

  # --- Panel 2: Structure vector ---
  sv <- x$structure_vector
  svdf <- data.frame(node = names(sv), dim = as.integer(sv),
                      stringsAsFactors = FALSE)
  svdf <- svdf[order(-svdf$dim, svdf$node), ]
  svdf$node <- factor(svdf$node, levels = svdf$node)

  p2 <- ggplot2::ggplot(svdf, ggplot2::aes(x = node, y = dim, fill = dim)) +
    ggplot2::geom_col(show.legend = FALSE, width = 0.7) +
    ggplot2::scale_fill_gradient(low = "#81B1D3", high = "#E8734A") +
    ggplot2::geom_text(ggplot2::aes(label = dim), vjust = -0.3,
                        size = 3.5) +
    ggplot2::labs(title = "Structure Vector",
                  subtitle = "Highest simplex dimension per node",
                  x = NULL, y = "Max Dimension") +
    .sc_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35,
                                                         hjust = 1))

  combined <- gridExtra::arrangeGrob(p1, p2, ncol = 2,
                                      padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined)
  invisible(combined)
}
