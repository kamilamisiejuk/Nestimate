# ---- Centrality Stability ----

#' Centrality Stability Coefficient (CS-coefficient)
#'
#' @description
#' Estimates the stability of centrality indices under case-dropping.
#' For each drop proportion, sequences are randomly removed and the
#' network is re-estimated. The correlation between the original and
#' subset centrality values is computed. The CS-coefficient is the
#' maximum proportion of cases that can be dropped while maintaining
#' a correlation above \code{threshold} in at least \code{certainty}
#' of bootstrap samples.
#'
#' For transition methods, uses pre-computed per-sequence count matrices
#' for fast resampling. Strength centralities (InStrength, OutStrength)
#' are computed directly from the matrix without igraph.
#'
#' @param x A \code{netobject} from \code{\link{build_network}}.
#' @param measures Character vector. Centrality measures to assess.
#'   Built-in: \code{"InStrength"}, \code{"OutStrength"}, \code{"Betweenness"},
#'   \code{"InCloseness"}, \code{"OutCloseness"}, \code{"Closeness"}.
#'   \code{"Closeness"} is defined only for undirected networks;
#'   \code{"InCloseness"}/\code{"OutCloseness"} only for directed
#'   networks (requesting the wrong one for the network's directedness
#'   is an error). Custom measures beyond these are valid only when a
#'   \code{centrality_fn} is supplied to resolve them.
#'   Default: \code{c("InStrength", "OutStrength", "Betweenness")}.
#' @param iter Integer. Number of bootstrap iterations per drop
#'   proportion (default: 1000).
#' @param drop_prop Numeric vector. Proportions of cases to drop
#'   (default: \code{seq(0.1, 0.9, by = 0.1)}).
#' @param threshold Numeric. Minimum correlation to consider stable
#'   (default: 0.7).
#' @param certainty Numeric. Required proportion of iterations above
#'   threshold (default: 0.95).
#' @param method Character. Correlation method: \code{"pearson"},
#'   \code{"spearman"}, or \code{"kendall"} (default: \code{"pearson"}).
#' @param centrality_fn Optional function. A custom centrality function
#'   that takes a weight matrix and returns a named list of centrality
#'   vectors. When \code{NULL} (default), all built-in measures are
#'   computed internally: \code{"InStrength"}/\code{"OutStrength"} via
#'   \code{colSums}/\code{rowSums}, and \code{"Betweenness"}/
#'   \code{"InCloseness"}/\code{"OutCloseness"}/\code{"Closeness"} via an
#'   internal Floyd-Warshall shortest-path routine. When provided, the
#'   function is called as \code{centrality_fn(mat)} and is used only for
#'   requested measures that are not one of the six built-ins; it should
#'   return a named list (e.g., \code{list(my_metric = ...)}).
#' @param loops Logical. If \code{FALSE} (default), self-loops (diagonal)
#'   are excluded from centrality computation. This does not modify the
#'   stored matrix.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_stability"} containing:
#' \describe{
#'   \item{cs}{Named numeric vector of CS-coefficients per measure.}
#'   \item{correlations}{Named list of matrices (iter x n_prop) of
#'     correlation values per measure.}
#'   \item{measures}{Character vector of measures assessed.}
#'   \item{drop_prop}{Drop proportions used.}
#'   \item{threshold}{Stability threshold.}
#'   \item{certainty}{Required certainty level.}
#'   \item{iter}{Number of iterations.}
#'   \item{method}{Correlation method.}
#' }
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' cs <- centrality_stability(net, iter = 100, seed = 42,
#'   measures = c("InStrength", "OutStrength"))
#' print(cs)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{network_reliability}}
#'
#' @importFrom stats cor sd
#' @export
centrality_stability <- function(x,
                                 measures = c("InStrength", "OutStrength",
                                              "Betweenness"),
                                 iter = 1000L,
                                 drop_prop = seq(0.1, 0.9, by = 0.1),
                                 threshold = 0.7,
                                 certainty = 0.95,
                                 method = "pearson",
                                 centrality_fn = NULL,
                                 loops = FALSE,
                                 seed = NULL) {

  # ---- Input validation ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      centrality_stability(net, measures = measures, iter = iter,
                           drop_prop = drop_prop, threshold = threshold,
                           certainty = certainty, method = method,
                           centrality_fn = centrality_fn, loops = loops,
                           seed = seed)
    })
    class(out) <- c("net_stability_group", "list")
    return(out)
  }
  if (!inherits(x, "netobject")) {
    stop("'x' must be a netobject from build_network().", call. = FALSE)
  }
  if (is.null(x$data)) {
    stop("netobject does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(drop_prop), all(drop_prop > 0), all(drop_prop < 1),
    is.numeric(threshold), length(threshold) == 1,
    threshold >= 0, threshold <= 1,
    is.numeric(certainty), length(certainty) == 1,
    certainty >= 0, certainty <= 1
  )
  iter <- as.integer(iter)
  method <- match.arg(method, c("pearson", "spearman", "kendall"))

  if (!is.null(centrality_fn)) {
    stopifnot("centrality_fn must be a function" = is.function(centrality_fn))
  }

  valid_measures <- c("InStrength", "OutStrength", "Betweenness",
                       "Closeness", "InCloseness", "OutCloseness")
  bad <- setdiff(measures, valid_measures)
  # A custom centrality_fn resolves non-builtin measure names (mirrors
  # net_centrality, which has no whitelist). Only reject unknown names
  # when no centrality_fn is supplied; otherwise .compute_centralities
  # routes them through centrality_fn.
  if (length(bad) > 0L && is.null(centrality_fn)) {
    stop("Unknown measures: ", paste(bad, collapse = ", "),
         ". Options: ", paste(valid_measures, collapse = ", "),
         ". Custom measures require centrality_fn.",
         call. = FALSE)
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Setup ----
  net_method <- .resolve_method_alias(x$method)
  states <- x$nodes$label
  n_states <- length(states)
  is_transition <- net_method %in% c("relative", "frequency", "co_occurrence")
  is_relative <- net_method == "relative"
  scaling <- x$scaling
  thresh <- x$threshold

  # ---- Original centralities ----
  orig_cents <- .compute_centralities(x$weights, states, x$directed, measures,
                                       loops, centrality_fn)

  # Drop measures with zero variance (e.g. OutStrength for relative networks)
  keep <- vapply(measures, function(m) sd(orig_cents[[m]]) > 0, logical(1))
  if (!any(keep)) {
    warning("All centrality measures have zero variance. ",
            "No stability can be assessed.", call. = FALSE)
    result <- list(
      cs = stats::setNames(rep(0, length(measures)), measures),
      correlations = stats::setNames(
        lapply(measures, function(m) {
          matrix(NA_real_, nrow = iter, ncol = length(drop_prop))
        }), measures),
      measures = measures,
      drop_prop = drop_prop,
      threshold = threshold,
      certainty = certainty,
      iter = iter,
      method = method
    )
    class(result) <- "net_stability"
    return(result)
  }
  measures <- measures[keep]
  orig_cents <- orig_cents[measures]

  # ---- Pre-compute for transitions ----
  if (is_transition) {
    trans_2d <- .precompute_per_sequence(x$data, net_method, x$params, states)
    n_seq <- nrow(trans_2d)
  } else {
    data <- x$data
    n_seq <- nrow(data)
    estimator <- get_estimator(net_method)
    params <- x$params
    level <- x$level
    id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")
  }

  # ---- Build matrix from subset (transition fast path) ----
  build_matrix_transition <- function(idx) {
    counts <- colSums(trans_2d[idx, , drop = FALSE])
    mat <- matrix(counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov
    if (thresh > 0) mat[abs(mat) < thresh] <- 0 # nocov
    dimnames(mat) <- list(states, states)
    mat
  }

  # ---- Build matrix from subset (association path) ----
  build_matrix_association <- function(idx) {
    sub_data <- data[idx, , drop = FALSE]
    if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
      sub_data <- tryCatch( # nocov start
        .decompose_multilevel(sub_data, id_col = id_col, level = level),
        error = function(e) NULL
      )
      if (is.null(sub_data)) return(NULL) # nocov end
    }
    est <- tryCatch(
      do.call(estimator$fn, c(list(data = sub_data), params)),
      error = function(e) NULL
    )
    if (is.null(est)) return(NULL)
    mat <- est$matrix[states, states] # nocov start
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
    if (thresh > 0) mat[abs(mat) < thresh] <- 0
    mat # nocov end
  }

  build_matrix <- if (is_transition) build_matrix_transition
                  else build_matrix_association

  # ---- Run stability assessment ----
  n_prop <- length(drop_prop)
  n_measures <- length(measures)

  # Pre-allocate correlation storage
  corr_storage <- lapply(measures, function(m) {
    matrix(NA_real_, nrow = iter, ncol = n_prop)
  })
  names(corr_storage) <- measures

  for (p_idx in seq_len(n_prop)) {
    prop <- drop_prop[p_idx]
    n_keep <- n_seq - floor(n_seq * prop)
    if (n_keep < 2L) next # nocov

    # vapply: each iteration returns correlation vector (one per measure)
    corr_mat <- vapply(seq_len(iter), function(i) {
      idx <- sample.int(n_seq, n_keep, replace = FALSE)
      mat <- build_matrix(idx)
      if (is.null(mat)) return(rep(NA_real_, n_measures))

      sub_cents <- .compute_centralities(mat, states, x$directed, measures,
                                          loops, centrality_fn)

      vapply(measures, function(m) {
        sc <- sub_cents[[m]]
        oc <- orig_cents[[m]]
        if (sd(sc) == 0) return(NA_real_)
        cor(sc, oc, method = method, use = "complete.obs")
      }, numeric(1))
    }, numeric(n_measures))

    # corr_mat is n_measures x iter, store each measure's column
    if (n_measures == 1L) {
      corr_storage[[1]][, p_idx] <- corr_mat
    } else {
      for (k in seq_len(n_measures)) {
        corr_storage[[measures[k]]][, p_idx] <- corr_mat[k, ]
      }
    }
  }

  # ---- Compute CS coefficients ----
  cs <- vapply(measures, function(m) {
    .calculate_cs(corr_storage[[m]], threshold, certainty, drop_prop)
  }, numeric(1))

  result <- list(
    cs = cs,
    correlations = corr_storage,
    measures = measures,
    drop_prop = drop_prop,
    threshold = threshold,
    certainty = certainty,
    iter = iter,
    method = method
  )
  class(result) <- "net_stability"
  result
}


# ---- Helpers ----

#' Compute centralities from a weight matrix
#'
#' InStrength and OutStrength use colSums/rowSums (no dependencies).
#' Other measures (Betweenness, Closeness, InCloseness, OutCloseness)
#' require a user-supplied \code{centrality_fn}.
#'
#' @param mat Weight matrix.
#' @param states Character vector of state names.
#' @param directed Logical.
#' @param measures Character vector of requested measures.
#' @param loops Logical. If FALSE, zero out diagonal before computing.
#' @param centrality_fn Optional function taking a weight matrix and
#'   returning a named list of centrality vectors.
#' @noRd
.compute_centralities <- function(mat, states, directed, measures,
                                  loops = FALSE, centrality_fn = NULL) {
  n <- length(states)
  if (!loops) diag(mat) <- 0
  result <- list()

  # .closeness() is documented to return a single "Closeness" only for
  # UNDIRECTED networks and "InCloseness"/"OutCloseness" only for DIRECTED
  # networks. Requesting the wrong one for the network's directedness
  # would otherwise silently produce no column (or crash a downstream
  # sd() check). Error cleanly with the correct measure name instead.
  if (directed && "Closeness" %in% measures) {
    stop("'Closeness' is defined only for undirected networks. ",
         "For directed networks use 'InCloseness' and/or 'OutCloseness'.",
         call. = FALSE)
  }
  if (!directed && any(c("InCloseness", "OutCloseness") %in% measures)) {
    stop("'InCloseness'/'OutCloseness' are defined only for directed ",
         "networks. For undirected networks use 'Closeness'.",
         call. = FALSE)
  }

  # Built-in matrix-based centralities (no dependencies)
  builtin <- c("InStrength", "OutStrength",
               "Betweenness", "InCloseness", "OutCloseness", "Closeness")

  if ("InStrength"  %in% measures) result[["InStrength"]]  <- colSums(mat)
  if ("OutStrength" %in% measures) result[["OutStrength"]] <- rowSums(mat)

  path_measures <- intersect(measures,
                             c("Betweenness", "InCloseness", "OutCloseness",
                               "Closeness"))
  if (length(path_measures) > 0L) {
    path_mat <- abs(mat)
    if ("Betweenness" %in% path_measures) {
      result[["Betweenness"]] <- .betweenness(path_mat, directed = directed,
                                               invert = TRUE)
    }
    cl_measures <- intersect(path_measures, c("InCloseness", "OutCloseness",
                                               "Closeness"))
    if (length(cl_measures) > 0L) {
      cl <- .closeness(path_mat, directed = directed, invert = TRUE)
      for (m in cl_measures) result[[m]] <- cl[[m]]
    }
  }

  # External centralities via centrality_fn (user-supplied measures only)
  external <- setdiff(measures, builtin)
  if (length(external) > 0L) {
    if (is.null(centrality_fn)) {
      stop("centrality_fn is required for measures: ",
           paste(external, collapse = ", "),
           ". Provide a function that takes a weight matrix and returns ",
           "a named list of centrality vectors.", call. = FALSE)
    }
    custom <- centrality_fn(mat)
    if (!is.list(custom)) { # nocov start
      stop("centrality_fn must return a named list.", call. = FALSE)
    } # nocov end
    for (m in external) {
      if (is.null(custom[[m]])) { # nocov start
        stop("centrality_fn did not return measure '", m, "'.",
             call. = FALSE)
      } # nocov end
      result[[m]] <- custom[[m]]
    }
  }

  result
}


#' Calculate CS-coefficient from correlation matrix
#' @noRd
.calculate_cs <- function(corr_mat, threshold, certainty, drop_prop) {
  prop_above <- colMeans(corr_mat >= threshold, na.rm = TRUE)
  valid <- which(prop_above >= certainty)
  if (length(valid) == 0L) 0 else drop_prop[max(valid)]
}


# ---- S3 Methods ----

#' Print Method for net_stability
#'
#' @param x A \code{net_stability} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#' print(cs)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' stab <- centrality_stability(net, measures = c("InStrength","OutStrength"),
#'                               iter = 10)
#' print(stab)
#' }
#'
#' @export
print.net_stability <- function(x, ...) {
  cat(sprintf("Centrality Stability (%d iterations, threshold = %.1f)\n",
              x$iter, x$threshold))
  cat(sprintf("  Drop proportions: %s\n",
              paste(x$drop_prop, collapse = ", ")))
  cat("\n  CS-coefficients:\n")
  for (m in names(x$cs)) {
    cat(sprintf("    %-15s  %.2f\n", m, x$cs[m]))
  }
  invisible(x)
}

#' Print Method for net_stability_group
#'
#' @param x A `net_stability_group` (returned by `centrality_stability()`
#'   when called on a `netobject_group` or an `mcml`).
#' @param ... Additional arguments (ignored).
#' @return The input `x` invisibly.
#' @export
print.net_stability_group <- function(x, ...) {
  measures <- unique(unlist(lapply(x, function(e) names(e$cs))))
  cs_mat <- sapply(x, function(e) e$cs[measures])
  rownames(cs_mat) <- measures
  cat(sprintf("Centrality Stability (%d networks, threshold = %.2f)\n",
              length(x), x[[1]]$threshold))
  print(round(t(cs_mat), 2))
  invisible(x)
}


#' Summary Method for net_stability_group
#'
#' @description
#' Per-network stability as a tidy data frame. Stacks \code{summary()}
#' results for each network with a \code{group} column.
#'
#' @param object A \code{net_stability_group} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with columns \code{group}, \code{measure},
#'   \code{drop_prop}, \code{mean_cor}, \code{sd_cor}, \code{prop_above}.
#' @export
summary.net_stability_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df <- summary(object[[nm]])
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' Summary Method for net_stability
#'
#' @description
#' Returns the mean correlation at each drop proportion for each measure.
#'
#' @param object A \code{net_stability} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with columns \code{measure}, \code{drop_prop},
#'   \code{mean_cor}, \code{sd_cor}, \code{prop_above}.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#' summary(cs)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' stab <- centrality_stability(net, measures = c("InStrength","OutStrength"),
#'                               iter = 10)
#' summary(stab)
#' }
#'
#' @export
summary.net_stability <- function(object, ...) {
  rows <- do.call(rbind, lapply(object$measures, function(m) {
    corr_mat <- object$correlations[[m]]
    do.call(rbind, lapply(seq_along(object$drop_prop), function(j) {
      vals <- corr_mat[, j]
      data.frame(
        measure = m,
        drop_prop = object$drop_prop[j],
        mean_cor = mean(vals, na.rm = TRUE),
        sd_cor = sd(vals, na.rm = TRUE),
        prop_above = mean(vals >= object$threshold, na.rm = TRUE),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }))
  }))
  rownames(rows) <- NULL
  rows
}


#' Plot Method for net_stability
#'
#' @description
#' Plots mean correlation vs drop proportion for each centrality measure.
#' The CS-coefficient is marked where the curve crosses the threshold.
#'
#' @param x A \code{net_stability} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' cs <- centrality_stability(net, iter = 10, drop_prop = 0.3)
#' plot(cs)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' stab <- centrality_stability(net, measures = c("InStrength","OutStrength"),
#'                               iter = 10)
#' plot(stab)
#' }
#'
#' @export
plot.net_stability <- function(x, ...) {
  summ <- summary(x)

  p <- ggplot2::ggplot(summ, ggplot2::aes(
    x = .data$drop_prop, y = .data$mean_cor,
    color = .data$measure)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = .data$mean_cor - .data$sd_cor,
        ymax = pmin(.data$mean_cor + .data$sd_cor, 1),
        fill = .data$measure),
      alpha = 0.15, color = NA
    ) +
    ggplot2::geom_hline(yintercept = x$threshold,
                        linetype = "dashed", color = "grey40") +
    ggplot2::annotate("text", x = max(x$drop_prop), y = x$threshold,
                      label = sprintf("threshold = %.1f", x$threshold),
                      hjust = 1, vjust = -0.5, size = 3, color = "grey40") +
    ggplot2::scale_x_continuous(
      breaks = x$drop_prop,
      labels = x$drop_prop
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Proportion dropped",
      y = sprintf("Mean correlation (%s)", x$method),
      title = "Centrality Stability",
      color = "Measure", fill = "Measure"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}
