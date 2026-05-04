# ==============================================================================
# cluster_choice() -- k + dissimilarity + method sweep for distance clustering.
#
# Parallel to compare_mmm() in R/mmm.R: one call, table back, print/summary/plot
# methods on the table. Where compare_mmm sweeps k for one MMM family,
# cluster_choice sweeps the cartesian product of (k, dissimilarity, method)
# for distance clustering. "all" sentinels expand to .clustering_metrics /
# .clustering_methods so users don't have to type out every option.
#
# Reuses build_clusters() per row -- no new metric formulas. Silhouette is
# already computed inside build_clusters; mean within-cluster distance is the
# same lower.tri block used in summary.net_clustering.
# ==============================================================================

#' Cluster Choice -- sweep k, dissimilarity and method
#'
#' One-call sweep across any combination of k, dissimilarity metric, and
#' clustering algorithm for distance-based sequence clustering. Mirrors
#' \code{\link{compare_mmm}} for model-based clustering: returns a data
#' frame with one row per swept configuration, a \code{best} marker on
#' the silhouette-max row in the print method, and a \code{plot()} that
#' adapts to the swept axes.
#'
#' @param data Sequence data (data frame or matrix) -- forwarded to
#'   \code{\link{build_clusters}}.
#' @param k Integer vector of cluster counts to sweep. Default
#'   \code{2:5}. Each value must be >= 2 and <= n - 1.
#' @param dissimilarity Character vector of dissimilarity metrics. Use
#'   \code{"all"} to expand to every supported metric:
#'   \code{c("hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine",}
#'   \code{"jaccard", "jw")}. Default \code{"hamming"}.
#' @param method Character vector of clustering algorithms. Use
#'   \code{"all"} to expand to every supported method:
#'   \code{c("pam", "ward.D2", "ward.D", "complete", "average", "single",}
#'   \code{"mcquitty", "median", "centroid")}. Default \code{"ward.D2"}.
#' @param ... Other arguments forwarded to
#'   \code{\link{build_clusters}} (\code{weighted}, \code{lambda},
#'   \code{q}, \code{p}, \code{seed}, \code{na_syms}, \code{covariates}).
#'   Note: \code{weighted = TRUE} only works with
#'   \code{dissimilarity = "hamming"} and is rejected up-front when
#'   sweeping mixed dissimilarities.
#' @return A \code{cluster_choice} object (a data.frame subclass) with
#'   one row per (k, dissimilarity, method) combination and columns:
#'   \describe{
#'     \item{k, dissimilarity, method}{The configuration for that row.}
#'     \item{silhouette}{Overall average silhouette width (from
#'       \code{cluster::silhouette}, computed inside
#'       \code{\link{build_clusters}}).}
#'     \item{mean_within_dist}{Size-weighted mean of within-cluster
#'       distances, in the units of the row's dissimilarity.}
#'     \item{min_size, max_size, size_ratio}{Cluster-size balance bounds
#'       and their ratio (\code{max / min}).}
#'   }
#' @seealso \code{\link{build_clusters}}, \code{\link{compare_mmm}} for
#'   the model-based equivalent, \code{\link{cluster_diagnostics}} for
#'   the post-fit diagnostic surface on a single clustering.
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 40, TRUE),
#'                    V2 = sample(c("A","B","C"), 40, TRUE))
#' cluster_choice(seqs, k = 2:4)
#' \donttest{
#' # Sweep dissimilarities at fixed k
#' cluster_choice(seqs, k = 3, dissimilarity = c("hamming", "lcs", "jaccard"))
#'
#' # Full grid of k x dissimilarity
#' cluster_choice(seqs, k = 2:4, dissimilarity = c("hamming", "lcs"))
#'
#' # "all" sentinel
#' cluster_choice(seqs, k = 3, dissimilarity = "all")
#' }
#' @export
cluster_choice <- function(data,
                            k             = 2:5,
                            dissimilarity = "hamming",
                            method        = "ward.D2",
                            ...) {
  # ---- Expand "all" sentinels ----------------------------------------
  dissimilarity <- .expand_all(dissimilarity, .clustering_metrics,
                                "dissimilarity")
  method        <- .expand_all(method,        .clustering_methods,
                                "method")
  k             <- as.integer(k)

  # ---- Pre-validate weighted/dissimilarity combination ---------------
  dots <- list(...)
  if (isTRUE(dots$weighted) &&
      any(dissimilarity != "hamming")) {
    stop("weighted = TRUE requires dissimilarity = \"hamming\". ",
         "Got: ", paste(dissimilarity, collapse = ", "), ".",
         call. = FALSE)
  }

  # ---- Build the cartesian-product grid ------------------------------
  # Order: k inside dissimilarity inside method (so the printed table
  # reads like rows of a method block with dissim sub-blocks).
  grid <- expand.grid(
    k             = k,
    dissimilarity = dissimilarity,
    method        = method,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS  = FALSE
  )

  # ---- Sweep ---------------------------------------------------------
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fit <- do.call(build_clusters, c(
      list(data          = data,
           k             = grid$k[[i]],
           dissimilarity = grid$dissimilarity[[i]],
           method        = grid$method[[i]]),
      dots))
    rows[[i]] <- .cluster_choice_row(fit)
  }
  out <- cbind(grid, do.call(rbind, rows))
  rownames(out) <- NULL
  class(out) <- c("cluster_choice", "data.frame")
  swept_axes <- character()
  if (length(k)             > 1L) swept_axes <- c(swept_axes, "k")
  if (length(dissimilarity) > 1L) swept_axes <- c(swept_axes, "dissimilarity")
  if (length(method)        > 1L) swept_axes <- c(swept_axes, "method")
  attr(out, "swept") <- swept_axes
  out
}

# Expand a user-supplied character vector that may contain "all" into the
# canonical full set. Deduplicates while preserving order. Empty / NA
# inputs fall through to build_clusters() which will raise the existing
# input-validation error.
#' @noRd
.expand_all <- function(values, full, name) {
  if (!is.character(values)) {
    values <- as.character(values)
  }
  if (any(values == "all")) {
    values <- unique(c(values[values != "all"], full))
  }
  if (!length(values)) {
    stop("'", name, "' is empty.", call. = FALSE)
  }
  values
}

#' @noRd
.cluster_choice_row <- function(fit) {
  sizes <- as.integer(fit$sizes %||% tabulate(fit$assignments))
  n     <- sum(sizes)
  k_    <- length(sizes)

  mean_within <- NA_real_
  if (!is.null(fit$distance) && n > 0L) {
    per_cluster <- .per_cluster_within_dist(fit$distance, fit$assignments,
                                              k_)
    mean_within <- stats::weighted.mean(per_cluster, sizes)
  }

  data.frame(
    silhouette       = as.numeric(fit$silhouette),
    mean_within_dist = mean_within,
    min_size         = as.integer(min(sizes)),
    max_size         = as.integer(max(sizes)),
    size_ratio       = if (min(sizes) > 0L) max(sizes) / min(sizes) else
                         NA_real_,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for cluster_choice
#'
#' @param x A \code{cluster_choice} object.
#' @param digits Integer. Decimal places for floating-point columns.
#'   Default \code{3L}.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.cluster_choice <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  swept  <- attr(x, "swept") %||%
            c("k", "dissimilarity", "method")

  if (length(swept) == 0L) {
    only_method <- unique(x$method)[1L]
    only_dissim <- unique(x$dissimilarity)[1L]
    cat(sprintf("Cluster Choice [%s / %s]\n\n", only_method, only_dissim))
  } else {
    cat(sprintf("Cluster Choice (sweep: %s)\n\n",
                paste(swept, collapse = " x ")))
  }

  disp <- data.frame(
    k             = x$k,
    dissimilarity = x$dissimilarity,
    method        = x$method,
    silhouette    = round(x$silhouette,       digits),
    within_dist   = round(x$mean_within_dist, digits),
    sizes         = sprintf("[%d, %d]", x$min_size, x$max_size),
    ratio         = round(x$size_ratio,       digits),
    stringsAsFactors = FALSE
  )

  best_idx <- which.max(x$silhouette)
  best <- rep("", nrow(disp))
  best[best_idx] <- "<-- best"
  disp$best <- best

  for (axis in c("k", "dissimilarity", "method")) {
    if (!(axis %in% swept)) disp[[axis]] <- NULL
  }

  print.data.frame(disp, row.names = FALSE, right = FALSE)

  invisible(x)
}

#' Summary Method for cluster_choice
#'
#' @param object A \code{cluster_choice} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with the swept configurations, all metrics, and
#'   a \code{best} character column flagging the silhouette-max row.
#' @export
summary.cluster_choice <- function(object, ...) {
  best_idx <- which.max(object$silhouette)
  best <- rep("", nrow(object))
  best[best_idx] <- "silhouette"
  out <- as.data.frame(object)
  out$best <- best
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Abbreviation helpers (display-only). The underlying $dissimilarity and
# $method columns always carry the canonical names; abbreviation is a
# rendering concern triggered by `abbrev = TRUE` on plot().
# ---------------------------------------------------------------------------

.dissimilarity_abbrev <- c(
  hamming = "ham", osa = "osa", lv  = "lv",  dl  = "dl",  lcs = "lcs",
  qgram   = "qgr", cosine = "cos", jaccard = "jac", jw  = "jw"
)

.method_abbrev <- c(
  pam = "pam", ward.D2 = "wD2", ward.D = "wD", complete = "cmp",
  average = "avg", single = "sng", mcquitty = "mcq", median = "mdn",
  centroid = "cen"
)

# Map canonical names -> abbreviations; pass through any unknown value.
#' @noRd
.abbrev_dissimilarity <- function(x) {
  out <- unname(.dissimilarity_abbrev[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

#' @noRd
.abbrev_method <- function(x) {
  out <- unname(.method_abbrev[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

#' Plot Method for cluster_choice
#'
#' Six explicit chart types plus a smart \code{"auto"} default. The user
#' picks the shape; the function does not editorialise (no "best"
#' annotation, no interpretive subtitles, no inferred recommendation).
#'
#' Type cheat-sheet:
#' \describe{
#'   \item{\code{"auto"}}{Default. Picks one of the others based on which
#'     axes were swept. k-only -> \code{"lines"}; one categorical axis
#'     swept -> \code{"bars"}; k plus one categorical -> \code{"lines"};
#'     k plus two categoricals -> \code{"facet"}; both categoricals
#'     without k -> \code{"heatmap"}.}
#'   \item{\code{"lines"}}{Silhouette across k (and \code{mean_within_dist}
#'     when \code{k} is the only swept axis), one line per non-k axis
#'     when present.}
#'   \item{\code{"bars"}}{Horizontal bar chart of silhouette per axis
#'     level. Bars sorted by silhouette.}
#'   \item{\code{"heatmap"}}{Tiled silhouette across two categorical
#'     axes. Requires both \code{dissimilarity} and \code{method} swept.}
#'   \item{\code{"tradeoff"}}{Scatter: silhouette (y) vs \code{size_ratio}
#'     (x). Works for any sweep; labels each point.}
#'   \item{\code{"facet"}}{Lines vs k, colour by one categorical axis,
#'     facet by another. Requires \code{k} plus two categoricals.}
#' }
#'
#' Asking for a type the data can't support raises an error pointing at
#' the alternatives.
#'
#' @param x A \code{cluster_choice} object.
#' @param type Character. One of \code{"auto"} (default), \code{"lines"},
#'   \code{"bars"}, \code{"heatmap"}, \code{"tradeoff"}, \code{"facet"}.
#' @param abbrev Logical. If \code{TRUE}, dissimilarity and method names
#'   shown on tick labels and point labels are shortened (e.g.
#'   \code{"hamming"} -> \code{"ham"}, \code{"ward.D2"} -> \code{"wD2"}).
#'   The legend shows the full canonical name. Default \code{FALSE}.
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object, invisibly.
#' @export
plot.cluster_choice <- function(x,
                                 type   = c("auto", "lines", "bars",
                                             "heatmap", "tradeoff",
                                             "facet"),
                                 abbrev = FALSE,
                                 ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) { # nocov start
    stop("Package 'ggplot2' required.", call. = FALSE)
  } # nocov end
  type  <- match.arg(type)
  swept <- attr(x, "swept") %||% c("k", "dissimilarity", "method")

  if (type == "auto") type <- .auto_choice_type(swept)
  .require_type_supported(type, swept)

  # Apply abbreviations to a working copy if requested. The original
  # object is never mutated -- callers keep canonical names in their
  # data frame.
  df <- as.data.frame(x)
  if (abbrev) {
    df$dissimilarity <- .abbrev_dissimilarity(df$dissimilarity)
    df$method        <- .abbrev_method(df$method)
  }

  p <- switch(type,
              lines    = .plot_choice_lines(df, swept),
              bars     = .plot_choice_bars(df, swept),
              heatmap  = .plot_choice_heatmap(df),
              tradeoff = .plot_choice_tradeoff(df, swept),
              facet    = .plot_choice_facet(df, swept))

  print(p)
  invisible(p)
}

#' @noRd
.auto_choice_type <- function(swept) {
  n_axes <- length(swept)
  if (n_axes <= 1L) {
    if ("k" %in% swept) "lines" else "bars"
  } else if (all(c("k", "dissimilarity", "method") %in% swept)) {
    "facet"
  } else if ("k" %in% swept) {
    "lines"
  } else {
    "heatmap"
  }
}

#' @noRd
.require_type_supported <- function(type, swept) {
  msg <- function(...) stop(..., call. = FALSE)
  switch(type,
         lines = if (!("k" %in% swept))
           msg("type = \"lines\" requires k to be swept (length > 1). ",
               "Use type = \"bars\" or type = \"heatmap\" for a fixed-k ",
               "sweep."),
         heatmap = if (!all(c("dissimilarity", "method") %in% swept))
           msg("type = \"heatmap\" requires both dissimilarity and method ",
               "to be swept. Use type = \"bars\" for one categorical axis ",
               "or type = \"lines\" for a k sweep."),
         facet = if (!("k" %in% swept) || length(swept) < 3L)
           msg("type = \"facet\" requires k plus two categorical axes ",
               "(dissimilarity and method) to be swept. Use type = ",
               "\"lines\" for k + one categorical, or type = \"heatmap\" ",
               "for two categoricals at fixed k."),
         bars     = invisible(),
         tradeoff = invisible())
}

# ---------------------------------------------------------------------------
# Per-type plot builders (each returns a ggplot, no editorialising)
# ---------------------------------------------------------------------------

#' @noRd
.plot_choice_lines <- function(df, swept) {
  # k-only sweep: silhouette and within_dist as two y-series.
  if (identical(swept, "k")) {
    long <- data.frame(
      k = rep(df$k, 2L),
      value = c(df$silhouette, df$mean_within_dist),
      metric = rep(c("silhouette", "mean_within_dist"),
                   each = nrow(df)),
      stringsAsFactors = FALSE
    )
    return(
      ggplot2::ggplot(long, ggplot2::aes(x = .data$k, y = .data$value,
                                           colour = .data$metric,
                                           group  = .data$metric)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_point(size = 2.8) +
        ggplot2::scale_x_continuous(breaks = unique(df$k)) +
        ggplot2::labs(title = "Cluster Choice", x = "k", y = NULL,
                      colour = NULL) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "bottom")
    )
  }

  # k + one categorical axis: silhouette vs k, one line per axis level.
  colour_axis <- if ("dissimilarity" %in% swept) "dissimilarity" else
                 if ("method"        %in% swept) "method"        else NA

  aes_obj <- if (!is.na(colour_axis))
    ggplot2::aes(x = .data$k, y = .data$silhouette,
                 colour = .data[[colour_axis]],
                 group  = .data[[colour_axis]])
  else
    ggplot2::aes(x = .data$k, y = .data$silhouette)

  ggplot2::ggplot(df, aes_obj) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.8) +
    ggplot2::scale_x_continuous(breaks = unique(df$k)) +
    ggplot2::scale_colour_viridis_d(option = "D", end = 0.85) +
    ggplot2::labs(title = "Cluster Choice", x = "k", y = "silhouette",
                  colour = colour_axis) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

#' @noRd
.plot_choice_bars <- function(df, swept) {
  # Pick the categorical axis to bar against. If multiple axes vary,
  # flatten them into one composite "k / dissim / method" label.
  axis <- if (identical(swept, "dissimilarity")) "dissimilarity"
          else if (identical(swept, "method"))   "method"
          else if (identical(swept, "k"))        "k"
          else NA

  if (is.na(axis)) {
    label_parts <- list()
    if ("k"             %in% swept) label_parts$k    <- as.character(df$k)
    if ("dissimilarity" %in% swept) label_parts$diss <- df$dissimilarity
    if ("method"        %in% swept) label_parts$meth <- df$method
    df$.row <- do.call(paste, c(label_parts, sep = " / "))
    axis <- ".row"
  }

  ggplot2::ggplot(df,
                   ggplot2::aes(x = stats::reorder(.data[[axis]],
                                                     .data$silhouette),
                                  y = .data$silhouette)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f",
                                                      .data$silhouette)),
                        hjust = -0.15, size = 3.2, colour = "grey25") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult =
                                                              c(0, 0.18))) +
    ggplot2::labs(title = "Cluster Choice",
                  x = if (axis == ".row") NULL else axis,
                  y = "silhouette") +
    ggplot2::theme_minimal(base_size = 12)
}

#' @noRd
.plot_choice_heatmap <- function(df) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$dissimilarity,
                                     y = .data$method,
                                     fill = .data$silhouette)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f",
                                                      .data$silhouette)),
                        size = 3.4, colour = "grey15") +
    ggplot2::scale_fill_viridis_c(option = "D", end = 0.95) +
    ggplot2::labs(title = "Cluster Choice", x = "dissimilarity",
                  y = "method", fill = "silhouette") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' @noRd
.plot_choice_tradeoff <- function(df, swept) {
  # Build a label that names whichever axes vary, so the scatter is
  # readable without a legend.
  parts <- list()
  if ("k"             %in% swept) parts$k    <- sprintf("k=%d", df$k)
  if ("dissimilarity" %in% swept) parts$diss <- df$dissimilarity
  if ("method"        %in% swept) parts$meth <- df$method
  df$.label <- do.call(paste, c(parts, sep = "/"))

  colour_axis <- if ("dissimilarity" %in% swept) "dissimilarity" else
                 if ("method"        %in% swept) "method"        else
                 if ("k"             %in% swept) "k"             else NA

  aes_obj <- if (!is.na(colour_axis))
    ggplot2::aes(x = .data$size_ratio, y = .data$silhouette,
                 colour = .data[[colour_axis]])
  else
    ggplot2::aes(x = .data$size_ratio, y = .data$silhouette)

  p <- ggplot2::ggplot(df, aes_obj) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = .data$.label),
                        hjust = -0.15, size = 3, colour = "grey25") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult =
                                                              c(0.05, 0.18))) +
    ggplot2::labs(title = "Cluster Choice",
                  x = "size_ratio", y = "silhouette",
                  colour = if (!is.na(colour_axis)) colour_axis else NULL) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.na(colour_axis) &&
      colour_axis %in% c("dissimilarity", "method")) {
    p <- p + ggplot2::scale_colour_viridis_d(option = "D", end = 0.85)
  }
  p
}

#' @noRd
.plot_choice_facet <- function(df, swept) {
  # k + dissim + method: lines vs k, colour = dissim, facet = method.
  ggplot2::ggplot(df, ggplot2::aes(x = .data$k, y = .data$silhouette,
                                     colour = .data$dissimilarity,
                                     group  = .data$dissimilarity)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_continuous(breaks = unique(df$k)) +
    ggplot2::scale_colour_viridis_d(option = "D", end = 0.85) +
    ggplot2::facet_wrap(~ method) +
    ggplot2::labs(title = "Cluster Choice", x = "k", y = "silhouette",
                  colour = "dissimilarity") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"))
}
