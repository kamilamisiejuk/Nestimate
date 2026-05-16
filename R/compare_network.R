# ---- Network comparison (descriptive, no inference) ----
#
# Mirrors `tna::compare()` exactly: same scaling menu, edge-level metrics,
# summary metrics (weight deviations, correlations, dissimilarities,
# similarities, pattern agreements), network metrics (via summary.netobject
# from R/summary_network.R), and optional centrality differences.
#
# All metrics are descriptive. Inferential comparisons (NCT, permutation
# tests, sequence comparisons) live in their own files.

#' Compare two networks descriptively
#'
#' Computes a battery of descriptive comparison metrics between two networks
#' or two weight matrices: weight deviations (mean / median / RMS / max
#' absolute difference, relative mean absolute difference, coefficient-of-
#' variation ratio), four correlation measures (Pearson, Spearman, Kendall,
#' distance correlation), five dissimilarity measures (Euclidean, Manhattan,
#' Canberra, Bray-Curtis, Frobenius), five similarity measures (Cosine,
#' Jaccard, Dice, Overlap, RV), pattern agreements, and side-by-side network
#' metrics. Optionally adds centrality differences and centrality
#' correlations.
#'
#' Mirrors `tna::compare()` numerically. Inputs are converted to weight
#' matrices and scaled before comparison; the choice of scaling determines
#' how weights from different estimators are placed on a common footing.
#'
#' @param x A `netobject`, `cograph_network`, or numeric square matrix.
#' @param y A `netobject`, `cograph_network`, or numeric square matrix.
#' @param scaling Scaling applied to both weight matrices before comparison.
#'   One of:
#'   \describe{
#'     \item{`"none"`}{Identity (default).}
#'     \item{`"minmax"`}{\eqn{(w - \min) / (\max - \min)}; maps to \eqn{[0, 1]}.}
#'     \item{`"max"`}{\eqn{w / \max(|w|)}; preserves sign.}
#'     \item{`"rank"`}{Min-max of average ranks; ordinal scaling.}
#'     \item{`"zscore"`}{\eqn{(w - \bar w) / s_w}; standard score.}
#'     \item{`"robust"`}{\eqn{(w - \mathrm{med}(w)) / \mathrm{mad}(w)};
#'       Huber-style robust z-score, resists outliers.}
#'     \item{`"log"`}{\eqn{\log(w)}; requires \eqn{w > 0}.}
#'     \item{`"log1p"`}{\eqn{\log(1 + w)}; admits \eqn{w \ge 0}.}
#'     \item{`"softmax"`}{Numerically stable softmax over the flattened
#'       vector.}
#'     \item{`"quantile"`}{Empirical CDF of the flattened vector.}
#'     \item{`"frobenius"`}{Divide the matrix by its Frobenius norm
#'       \eqn{\|W\|_F = \sqrt{\sum w_{ij}^2}}; matrix-level normalisation.}
#'     \item{`"row"`}{Row-stochastic normalisation (each row's absolute
#'       values sum to 1). Only meaningful for non-negative matrices; rows
#'       summing to zero are left unchanged.}
#'   }
#'   Scalings that produce negative weights (`zscore`, `robust`) are
#'   compatible with `network = TRUE` because the side-by-side metrics use
#'   Nestimate's base-R Floyd-Warshall, which handles negative weights.
#' @param measures Character vector of centrality measures to compare. Empty
#'   by default (no centrality block).
#' @param network Logical. Include side-by-side network metrics from
#'   `summary()`? Default `TRUE`.
#' @param ... Ignored.
#' @return A `net_comparison` object: a named list with `matrices`,
#'   `difference_matrix`, `edge_metrics`, `summary_metrics`, optionally
#'   `network_metrics`, `centrality_differences`, `centrality_correlations`.
#' @export
compare_model <- function(x, ...) UseMethod("compare_model")

#' @rdname compare_model
#' @export
compare_model.netobject <- function(x, y, scaling = "none", measures = character(0),
                              network = TRUE, ...) {
  .compare_impl(.weights_of(x), .weights_of(y), scaling, measures, network)
}

#' @rdname compare_model
#' @export
compare_model.cograph_network <- function(x, y, scaling = "none",
                                    measures = character(0), network = TRUE,
                                    ...) {
  .compare_impl(.weights_of(x), .weights_of(y), scaling, measures, network)
}

#' @rdname compare_model
#' @export
compare_model.matrix <- function(x, y, scaling = "none", measures = character(0),
                           network = TRUE, ...) {
  .compare_impl(.weights_of(x), .weights_of(y), scaling, measures, network)
}

#' Compare two networks within a netobject_group
#'
#' Selects two members of a `netobject_group` (by index or name) and
#' dispatches to `compare_model.netobject()`.
#'
#' @param x A `netobject_group`.
#' @param i Integer or character. Index or name of the first network.
#'   Default `1L`.
#' @param j Integer or character. Index or name of the second network.
#'   Default `2L`.
#' @param scaling See `compare_model()`.
#' @param measures See `compare_model()`.
#' @param network See `compare_model()`.
#' @param ... Passed to `compare_model.netobject()`.
#' @return A `net_comparison` object.
#' @export
compare_model.netobject_group <- function(x, i = 1L, j = 2L, scaling = "none",
                                    measures = character(0), network = TRUE,
                                    ...) {
  stopifnot(inherits(x, "netobject_group"))
  if (is.character(i)) i <- match(i, names(x))
  if (is.character(j)) j <- match(j, names(x))
  stopifnot(!is.na(i), !is.na(j), i >= 1L, i <= length(x),
            j >= 1L, j <= length(x), i != j)
  compare_model(x[[i]], x[[j]], scaling = scaling, measures = measures,
          network = network, ...)
}


# ---- Internal: extract weight matrix ----

.weights_of <- function(z) {
  if (is.matrix(z)) {
    stopifnot(nrow(z) == ncol(z), is.numeric(z))
    if (is.null(rownames(z)) || is.null(colnames(z))) {
      labs <- paste0("V", seq_len(nrow(z)))
      rownames(z) <- colnames(z) <- labs
    }
    return(z)
  }
  if (inherits(z, "netobject") || inherits(z, "cograph_network")) {
    return(z$weights)
  }
  stop("`compare_model()` expects a netobject, cograph_network, or numeric matrix.",
       call. = FALSE)
}


# ---- Internal: the comparison machinery ----

.compare_impl <- function(x, y, scaling, measures, network) {
  stopifnot(identical(dim(x), dim(y)))
  stopifnot(!anyNA(x), !anyNA(y))
  n <- nrow(x)

  scaling <- match.arg(scaling, .compare_scaling_names())
  x_scaled <- .apply_compare_scaling(x, scaling)
  y_scaled <- .apply_compare_scaling(y, scaling)
  d <- x_scaled - y_scaled
  x_vec <- as.vector(x_scaled)
  y_vec <- as.vector(y_scaled)
  abs_diff <- abs(x_vec - y_vec)
  abs_x <- abs(x_vec)
  abs_y <- abs(y_vec)
  pos <- abs_x > 0 & abs_y > 0

  rn <- rownames(x_scaled)
  edges <- expand.grid(source = rn, target = rn,
                       KEEP.OUT.ATTRS = FALSE,
                       stringsAsFactors = FALSE)
  edges$weight_x <- x_vec
  edges$weight_y <- y_vec
  edges$raw_difference <- as.vector(d)
  edges$absolute_difference <- abs_diff
  edges$squared_difference <- abs_diff^2
  edges$relative_difference <- abs_diff / (x_vec + y_vec)
  edges$similarity_strength_index <- x_vec / y_vec
  edges$difference_index <- (x_vec - y_vec) / y_vec
  edges$rank_difference <- abs(rank(x_vec, na.last = "keep") -
                               rank(y_vec, na.last = "keep"))
  edges$percentile_difference <- abs(stats::ecdf(x_vec)(x_vec) -
                                     stats::ecdf(y_vec)(y_vec))
  edges$logarithmic_ratio <- log1p(x_vec) - log1p(y_vec)
  edges$standardized_weight_x <- (x_vec - mean(x_scaled)) / stats::sd(x_scaled)
  edges$standardized_weight_y <- (y_vec - mean(y_scaled)) / stats::sd(y_scaled)
  edges$standardized_score_inflation <-
    edges$standardized_weight_x / edges$standardized_weight_y

  # Single source of truth: the worker computes all 22 metrics; this block
  # just adds category labels and pretty names for the data.frame view.
  sim <- .network_similarity(x_scaled, y_scaled)
  summary_metrics <- data.frame(
    category = c(
      rep("Weight Deviations", 6L),
      rep("Correlations", 4L),
      rep("Dissimilarities", 5L),
      rep("Similarities", 5L),
      rep("Pattern Similarities", 2L)
    ),
    metric = c(
      "Mean Abs. Diff.", "Median Abs. Diff.", "RMS Diff.",
      "Max Abs. Diff.", "Rel. Mean Abs. Diff.", "CV Ratio",
      "Pearson", "Spearman", "Kendall", "Distance",
      "Euclidean", "Manhattan", "Canberra", "Bray-Curtis", "Frobenius",
      "Cosine", "Jaccard", "Dice", "Overlap", "RV",
      "Rank Agreement", "Sign Agreement"
    ),
    value = unname(sim[c(
      "mean_abs_diff", "median_abs_diff", "rms_diff",
      "max_abs_diff", "rel_mean_abs", "cv_ratio",
      "pearson", "spearman", "kendall", "distance_cor",
      "euclidean", "manhattan", "canberra", "bray_curtis", "frobenius",
      "cosine", "jaccard", "dice", "overlap", "rv",
      "rank_agreement", "sign_agreement"
    )]),
    stringsAsFactors = FALSE
  )

  out <- list(
    matrices = list(x = x_scaled, y = y_scaled),
    difference_matrix = d,
    edge_metrics = edges,
    summary_metrics = summary_metrics
  )

  if (network) {
    sx <- .summary_metrics_from_weights(x_scaled,
                                        directed = !isSymmetric(x_scaled))
    sy <- .summary_metrics_from_weights(y_scaled,
                                        directed = !isSymmetric(y_scaled))
    out$network_metrics <- data.frame(metric = sx$metric,
                                      x = sx$value,
                                      y = sy$value,
                                      stringsAsFactors = FALSE)
  }

  if (length(measures) > 0L) {
    cents <- .compare_centralities(x_scaled, y_scaled, measures)
    out$centrality_differences <- cents$differences
    out$centrality_correlations <- cents$correlations
  }

  structure(out, class = "net_comparison")
}


# ---- Scaling catalogue ----

# Names of every scaling option understood by `compare_model()`.
.compare_scaling_names <- function() {
  c("none", "minmax", "max", "rank", "zscore", "robust",
    "log", "log1p", "softmax", "quantile",
    "frobenius", "row")
}

# Apply a single scaling option to a weight matrix. Most options operate on
# the flattened vector and preserve dimnames; `frobenius` and `row` are
# matrix-level normalisations.
.apply_compare_scaling <- function(W, scaling) {
  if (scaling == "row") {
    rs <- rowSums(abs(W))
    nz <- rs > 0
    out <- W
    out[nz, ] <- W[nz, ] / rs[nz]
    return(out)
  }
  if (scaling == "frobenius") {
    fro <- sqrt(sum(W^2))
    return(if (fro == 0) W else W / fro)
  }
  vec_fn <- switch(scaling,
    none     = identity,
    minmax   = function(w) {
      rng <- range(w)
      if (diff(rng) == 0) w else (w - rng[1L]) / diff(rng)
    },
    max      = function(w) {
      m <- max(w, na.rm = TRUE)
      if (m == 0) w else w / m
    },
    rank     = function(w) {
      r <- rank(w, ties.method = "average")
      rng <- range(r)
      if (diff(rng) == 0) r else (r - rng[1L]) / diff(rng)
    },
    zscore   = function(w) {
      s <- stats::sd(w)
      if (s == 0) w - mean(w) else (w - mean(w)) / s
    },
    robust   = function(w) {
      m <- stats::mad(w)
      if (m == 0) w - stats::median(w) else (w - stats::median(w)) / m
    },
    log      = log,
    log1p    = log1p,
    softmax  = function(w) {
      mx <- max(w); e <- exp(w - mx); e / sum(e)
    },
    quantile = function(w) stats::ecdf(w)(w),
    stop("Unknown scaling: ", scaling, call. = FALSE)
  )
  out <- W
  out[] <- vec_fn(as.vector(W))
  out
}


# ---- Network similarity worker ----

# Single source of truth for descriptive similarity between two networks.
# Returns a named numeric vector. Inputs can be matrices (full 22 metrics)
# or vectors (subset that does not need matrix structure — RV and
# rank_agreement are skipped).
#
# Used by `.compare_impl()`, `.split_half_metrics()` (reliability), and
# the casedrop drop-fraction loop. When you add a metric, add it once
# here and every caller picks it up.
.network_similarity <- function(x, y, metrics = NULL) {
  all_metrics <- c(
    "mean_abs_diff", "median_abs_diff", "rms_diff", "max_abs_diff",
    "rel_mean_abs", "cv_ratio",
    "pearson", "spearman", "kendall", "distance_cor",
    "euclidean", "manhattan", "canberra", "bray_curtis", "frobenius",
    "cosine", "jaccard", "dice", "overlap", "rv",
    "rank_agreement", "sign_agreement"
  )
  if (is.null(metrics)) metrics <- all_metrics
  metrics <- match.arg(metrics, all_metrics, several.ok = TRUE)

  has_matrix <- is.matrix(x) && is.matrix(y)
  if (has_matrix) {
    stopifnot(identical(dim(x), dim(y)))
    x_mat <- x; y_mat <- y
    x_vec <- as.vector(x); y_vec <- as.vector(y)
    n <- nrow(x_mat)
  } else {
    x_vec <- as.numeric(x); y_vec <- as.numeric(y)
    stopifnot(length(x_vec) == length(y_vec))
    n <- length(x_vec)
  }

  # Cheap shared quantities
  abs_diff <- abs(x_vec - y_vec)
  needs_abs_xy <- any(metrics %in% c("canberra", "bray_curtis",
                                     "jaccard", "dice", "overlap"))
  if (needs_abs_xy) {
    abs_x <- abs(x_vec); abs_y <- abs(y_vec); pos <- abs_x > 0 & abs_y > 0
  }
  needs_sd <- any(metrics %in% c("pearson", "spearman", "kendall", "cv_ratio"))
  if (needs_sd) { sd_x <- stats::sd(x_vec); sd_y <- stats::sd(y_vec) }

  # Compute requested metrics on demand. Each branch only fires if asked.
  one <- function(m) switch(m,
    mean_abs_diff   = mean(abs_diff),
    median_abs_diff = stats::median(abs_diff),
    rms_diff        = sqrt(mean(abs_diff^2)),
    max_abs_diff    = max(abs_diff),
    rel_mean_abs    = {
      my <- mean(abs(y_vec))
      if (my == 0) NA_real_ else mean(abs_diff) / my
    },
    cv_ratio        = {
      mx <- mean(x_vec); my <- mean(y_vec)
      if (mx == 0 || sd_y == 0) NA_real_ else sd_x * my / (mx * sd_y)
    },
    pearson         = if (sd_x == 0 || sd_y == 0) NA_real_
                      else stats::cor(x_vec, y_vec, method = "pearson",
                                      use = "complete.obs"),
    spearman        = if (sd_x == 0 || sd_y == 0) NA_real_
                      else stats::cor(x_vec, y_vec, method = "spearman",
                                      use = "complete.obs"),
    kendall         = if (sd_x == 0 || sd_y == 0) NA_real_
                      else stats::cor(x_vec, y_vec, method = "kendall",
                                      use = "complete.obs"),
    distance_cor    = .distance_correlation(x_vec, y_vec),
    euclidean       = sqrt(sum(abs_diff^2)),
    manhattan       = sum(abs_diff),
    canberra        = sum(abs_diff[pos] / (abs_x[pos] + abs_y[pos])),
    bray_curtis     = {
      s <- sum(abs_x + abs_y)
      if (s == 0) NA_real_ else sum(abs_diff) / s
    },
    frobenius       = if (has_matrix) sqrt(sum(abs_diff^2)) / sqrt(n / 2)
                      else NA_real_,
    cosine          = {
      d <- sqrt(sum(x_vec^2)) * sqrt(sum(y_vec^2))
      if (d == 0) NA_real_ else sum(x_vec * y_vec) / d
    },
    jaccard         = {
      s <- sum(pmax(abs_x, abs_y))
      if (s == 0) NA_real_ else sum(pmin(abs_x, abs_y)) / s
    },
    dice            = {
      s <- sum(abs_x) + sum(abs_y)
      if (s == 0) NA_real_ else 2 * sum(pmin(abs_x, abs_y)) / s
    },
    overlap         = {
      m <- min(sum(abs_x), sum(abs_y))
      if (m == 0) NA_real_ else sum(pmin(abs_x, abs_y)) / m
    },
    rv              = if (has_matrix) .rv_coefficient(x_mat, y_mat)
                      else NA_real_,
    rank_agreement  = if (has_matrix)
                        mean(sign(diff(x_mat)) == sign(diff(y_mat)))
                      else NA_real_,
    sign_agreement  = mean(sign(x_vec) == sign(y_vec))
  )

  out <- vapply(metrics, one, numeric(1L))
  names(out) <- metrics
  out
}


# ---- Helpers ----

# Distance correlation per Szekely et al. (2007).
.distance_correlation <- function(x, y) {
  dist_x <- as.matrix(stats::dist(x, diag = TRUE, upper = TRUE))
  dist_y <- as.matrix(stats::dist(y, diag = TRUE, upper = TRUE))
  n <- ncol(dist_x)
  rm_x <- matrix(.rowMeans(dist_x, n, n), n, n)
  rm_y <- matrix(.rowMeans(dist_y, n, n), n, n)
  cm_x <- matrix(.colMeans(dist_x, n, n), n, n, byrow = TRUE)
  cm_y <- matrix(.colMeans(dist_y, n, n), n, n, byrow = TRUE)
  m_x <- mean(dist_x); m_y <- mean(dist_y)
  xx <- dist_x - rm_x - cm_x + m_x
  yy <- dist_y - rm_y - cm_y + m_y
  v_xy <- n^-2 * sum(xx * yy)
  v_x  <- n^-2 * sum(xx^2)
  v_y  <- n^-2 * sum(yy^2)
  v_xy / sqrt(v_x * v_y)
}

# RV coefficient per Robert & Escoufier (1976).
.rv_coefficient <- function(x, y) {
  x <- scale(x, scale = FALSE)
  y <- scale(y, scale = FALSE)
  xx <- tcrossprod(x)
  yy <- tcrossprod(y)
  sum(diag(xx %*% yy)) /
    sqrt(sum(diag(xx %*% xx)) * sum(diag(yy %*% yy)))
}

.compare_centralities <- function(x, y, measures) {
  cx <- centrality(.wrap_netobject(x, method = "compare_model", directed = TRUE,
                                   data = NULL))
  cy <- centrality(.wrap_netobject(y, method = "compare_model", directed = TRUE,
                                   data = NULL))
  cx_df <- as.data.frame(cx)
  cy_df <- as.data.frame(cy)
  keep <- intersect(measures, names(cx_df))
  if (length(keep) == 0L) {
    return(list(differences = data.frame(),
                correlations = data.frame()))
  }
  state <- if ("state" %in% names(cx_df)) cx_df$state
           else if ("name" %in% names(cx_df)) cx_df$name
           else seq_len(nrow(cx_df))
  diffs <- do.call(rbind, lapply(keep, function(m) {
    data.frame(state = state, centrality = m,
               x = cx_df[[m]], y = cy_df[[m]],
               difference = cx_df[[m]] - cy_df[[m]],
               stringsAsFactors = FALSE)
  }))
  cors <- do.call(rbind, lapply(keep, function(m) {
    r <- tryCatch(stats::cor(cx_df[[m]], cy_df[[m]], use = "complete.obs"),
                  error = function(e) NA_real_)
    data.frame(centrality = m, correlation = r, stringsAsFactors = FALSE)
  }))
  list(differences = diffs, correlations = cors)
}


# ---- print method ----

#' @export
print.net_comparison <- function(x, ...) {
  cat("Network comparison\n")
  cat("==================\n")
  cat("Summary metrics:\n")
  sm <- x$summary_metrics
  sm$value <- vapply(sm$value,
                     function(v) formatC(v, digits = 4L, format = "g"),
                     character(1L))
  print.data.frame(sm, row.names = FALSE)
  if (!is.null(x$network_metrics)) {
    cat("\nNetwork metrics (x vs y):\n")
    nm <- x$network_metrics
    nm$x <- .format_metric_values(nm$metric, nm$x)
    nm$y <- .format_metric_values(nm$metric, nm$y)
    print.data.frame(nm, row.names = FALSE)
  }
  if (!is.null(x$centrality_correlations)) {
    cat("\nCentrality correlations:\n")
    cc <- x$centrality_correlations
    cc$correlation <- vapply(cc$correlation,
                             function(v) formatC(v, digits = 4L, format = "g"),
                             character(1L))
    print.data.frame(cc, row.names = FALSE)
  }
  invisible(x)
}


# ---- plot method ----

#' Plot a network comparison
#'
#' Visualises a `net_comparison` object. Currently supports the edge-weight
#' scatterplot (default), with the diagonal reference (perfect agreement)
#' and the OLS regression line annotated by Pearson, Spearman, and Kendall
#' correlations.
#'
#' @param x A `net_comparison` object from [compare_model()].
#' @param type Character. One of `"scatter"` (default — edge-weight scatter
#'   with OLS fit and correlation overlay), `"heatmap"` (n by n grid of
#'   x - y differences using the diverging palette), `"diff_hist"`
#'   (histogram of |x - y| absolute differences with rug + density),
#'   `"weight_dist"` (overlaid distributions of |x| and |y| edge weights),
#'   or `"all"` (2 by 2 grid of all four panels; requires the gridExtra
#'   package).
#' @param combined When `type = "all"` and `combined = TRUE` (default),
#'   the four panels are stitched into a 2x2 gtable. When `FALSE`, returns
#'   a named list of the four ggplots so each can be printed, saved, or
#'   re-laid-out independently. Ignored for other `type` values.
#' @param ... Ignored.
#' @return A `ggplot` object; for `type = "all"` with `combined = TRUE`
#'   a `gtable` arranged 2 by 2; for `type = "all"` with `combined = FALSE`
#'   a named list of four ggplots.
#' @export
plot.net_comparison <- function(x,
                                type = c("scatter", "heatmap",
                                         "diff_hist", "weight_dist",
                                         "all"),
                                combined = TRUE,
                                ...) {
  type <- match.arg(type)
  stopifnot(is.logical(combined), length(combined) == 1L)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot.net_comparison().", call. = FALSE)
  }
  switch(type,
    scatter     = .plot_comparison_scatter(x),
    heatmap     = .plot_comparison_heatmap(x),
    diff_hist   = .plot_comparison_diff_hist(x),
    weight_dist = .plot_comparison_weight_dist(x),
    all         = .plot_comparison_all(x, combined = combined)
  )
}

# 2x2 grid of all four panels, stripped of titles so the layout reads as
# a single figure. Uses gridExtra::arrangeGrob (already in Suggests).
# combined = FALSE returns the four ggplots as a named list instead of
# stitching them — useful when each panel needs its own device or the
# user wants to save them separately.
.plot_comparison_all <- function(x, combined = TRUE) {
  strip <- function(p, label) {
    p +
      ggplot2::labs(title = label, subtitle = NULL) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 11,
                                                        face = "bold"))
  }
  panels <- list(
    scatter     = strip(.plot_comparison_scatter(x),     "Scatter"),
    heatmap     = strip(.plot_comparison_heatmap(x),     "Heatmap (x - y)"),
    diff_hist   = strip(.plot_comparison_diff_hist(x),   "|x - y|"),
    weight_dist = strip(.plot_comparison_weight_dist(x), "Weight distributions")
  )
  if (!combined) return(invisible(panels))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("type = 'all' with combined = TRUE requires the gridExtra package.",
         call. = FALSE)
  }
  combined_grob <- gridExtra::arrangeGrob(panels$scatter, panels$heatmap,
                                          panels$diff_hist, panels$weight_dist,
                                          nrow = 2L, ncol = 2L)
  grid::grid.draw(combined_grob)
  invisible(combined_grob)
}

# Edge-weight scatter with diagonal, OLS line, and correlation annotations.
# House defaults: theme_minimal(base_size = 12), linewidth = 0.4. The
# regression band uses the project diverging palette's high-end blue
# (#4A6FE3) so it is consistent with the heatmap palette to follow.
.plot_comparison_scatter <- function(x) {
  em <- x$edge_metrics
  sm <- x$summary_metrics
  pick <- function(label) {
    v <- sm$value[sm$metric == label]
    if (length(v) == 0L) NA_real_ else v[[1L]]
  }
  r_p <- pick("Pearson")
  r_s <- pick("Spearman")
  r_k <- pick("Kendall")
  ann <- sprintf("Pearson = %.3f\nSpearman = %.3f\nKendall = %.3f",
                 r_p, r_s, r_k)

  rng <- range(c(em$weight_x, em$weight_y), na.rm = TRUE)
  pad <- 0.04 * diff(rng)
  lo  <- rng[1L] - pad
  hi  <- rng[2L] + pad

  ggplot2::ggplot(em, ggplot2::aes(x = .data$weight_x, y = .data$weight_y)) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey60",
                         linewidth = 0.4) +
    ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x,
                         colour = "#4A6FE3", fill = "#4A6FE3", alpha = 0.15,
                         linewidth = 0.4) +
    ggplot2::geom_point(alpha = 0.7, size = 1.6, colour = "#333333") +
    ggplot2::annotate("text", x = lo, y = hi,
                      label = ann, hjust = 0, vjust = 1, size = 3.5,
                      family = "mono") +
    ggplot2::coord_cartesian(xlim = c(lo, hi), ylim = c(lo, hi)) +
    ggplot2::labs(x = "Network x: edge weight",
                  y = "Network y: edge weight",
                  title = "Edge-weight comparison",
                  subtitle = "Dashed: y = x. Blue: OLS fit + 95% CI.") +
    ggplot2::theme_minimal(base_size = 12)
}

# Heatmap of (x - y), the scaled difference matrix. Diverging palette
# matches the project house defaults (red -> white -> blue, midpoint = 0).
# Per CLAUDE.md, never touch the diagonal — it stays as estimated.
.plot_comparison_heatmap <- function(x) {
  d <- x$difference_matrix
  rn <- rownames(d); cn <- colnames(d)
  if (is.null(rn)) rn <- as.character(seq_len(nrow(d)))
  if (is.null(cn)) cn <- as.character(seq_len(ncol(d)))

  long <- data.frame(
    source = factor(rep(rn, times = ncol(d)), levels = rn),
    target = factor(rep(cn, each  = nrow(d)), levels = rev(cn)),
    diff   = as.vector(d),
    stringsAsFactors = FALSE
  )
  lim <- max(abs(long$diff), na.rm = TRUE)

  ggplot2::ggplot(long,
                  ggplot2::aes(x = .data$target, y = .data$source,
                               fill = .data$diff)) +
    ggplot2::geom_tile(colour = "grey90", linewidth = 0.4) +
    ggplot2::scale_fill_gradient2(low = "#D33F6A", mid = "white",
                                  high = "#4A6FE3", midpoint = 0,
                                  limits = c(-lim, lim),
                                  name = "x - y") +
    ggplot2::labs(x = "to", y = "from",
                  title = "Difference matrix (x - y)",
                  subtitle = "Red: x > y. Blue: y > x. White: agreement.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1))
}

# Histogram of absolute edge-weight differences |x - y| with rug and
# density overlay. Magnitude-of-disagreement view: 0 means perfect
# agreement at that cell; larger means bigger discrepancy. Rug guarantees
# each individual edge is visible even when most cells cluster near zero.
.plot_comparison_diff_hist <- function(x) {
  abs_diffs <- abs(as.vector(x$difference_matrix))
  .magnitude_density_plot(
    abs_diffs,
    title    = "Distribution of |edge-weight differences|",
    subtitle = "Mean (red), median (grey). Rug: individual edges.",
    x_lab    = "|x - y|"
  )
}

# Overlaid distributions of |edge weight| for the two scaled networks.
# Lets the reader see whether the two networks have similar weight
# magnitudes overall (same dynamic range, similar concentration of mass
# near zero) or whether one is systematically heavier-tailed than the
# other.
.plot_comparison_weight_dist <- function(x) {
  vx <- abs(as.vector(x$matrices$x))
  vy <- abs(as.vector(x$matrices$y))
  df <- data.frame(
    weight  = c(vx, vy),
    network = factor(rep(c("x", "y"), c(length(vx), length(vy))),
                     levels = c("x", "y"))
  )

  ggplot2::ggplot(df,
                  ggplot2::aes(x = .data$weight,
                               fill   = .data$network,
                               colour = .data$network)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 30, alpha = 0.4,
                            position = "identity",
                            colour = "white", linewidth = 0.4) +
    ggplot2::geom_density(linewidth = 0.6, fill = NA) +
    ggplot2::geom_rug(sides = "b", alpha = 0.5,
                      length = ggplot2::unit(0.025, "npc")) +
    ggplot2::scale_fill_manual(values = c(x = "#4A6FE3", y = "#D33F6A")) +
    ggplot2::scale_colour_manual(values = c(x = "#1d3a8a", y = "#8a1d3a")) +
    ggplot2::labs(x = "|edge weight|", y = "density",
                  title = "Distribution of |edge weights| per network",
                  subtitle = "Blue: x. Red: y. Rug: individual edges.") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

# Shared layout for the magnitude-style histogram with rug + density.
.magnitude_density_plot <- function(values, title, subtitle, x_lab) {
  values <- values[is.finite(values)]
  med <- stats::median(values)
  mu  <- mean(values)
  iqr <- stats::IQR(values)
  rng <- range(values)
  bw  <- if (iqr > 0) 2 * iqr / (length(values)^(1 / 3))
         else diff(rng) / 12
  nbins <- max(8L, min(30L, ceiling(diff(rng) / max(bw, .Machine$double.eps))))

  ggplot2::ggplot(data.frame(value = values),
                  ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = nbins, fill = "#4A6FE3", alpha = 0.55,
                            colour = "white", linewidth = 0.4) +
    ggplot2::geom_density(colour = "#1d3a8a", linewidth = 0.6,
                          fill = NA) +
    ggplot2::geom_rug(sides = "b", alpha = 0.6,
                      length = ggplot2::unit(0.03, "npc")) +
    ggplot2::geom_vline(xintercept = mu, colour = "#D33F6A",
                        linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = med, colour = "grey50",
                        linetype = "dashed", linewidth = 0.4) +
    ggplot2::annotate("text", x = mu, y = Inf,
                      label = sprintf(" mean = %.3f", mu),
                      hjust = 0, vjust = 1.5, size = 3.5,
                      colour = "#D33F6A", family = "mono") +
    ggplot2::annotate("text", x = med, y = Inf,
                      label = sprintf(" median = %.3f", med),
                      hjust = 0, vjust = 3, size = 3.5,
                      colour = "grey40", family = "mono") +
    ggplot2::labs(x = x_lab, y = "density",
                  title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12)
}
