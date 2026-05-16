# ==============================================================================
# cluster_diagnostics() -- unified clustering quality surface.
#
# Single entry point that returns a structured `net_cluster_diagnostics`
# object regardless of whether the clustering came from build_clusters()
# (distance), build_mmm() (mmm), or one of the cluster_*() wrappers that
# return a netobject_group with an attr(, "clustering") of class
# `net_clustering` or `net_mmm_clustering`.
#
# This module only PACKAGES information that's already computed:
#  - .mmm_quality()    in R/mmm.R         -- avepp, entropy, classification_error
#  - cluster::silhouette() (already used inside plot.net_clustering)
#  - cl$silhouette / mmm$BIC / etc.       -- read straight off the source
#
# Helpers it reuses (do not reimplement):
#  - .cluster_table_lines() / .fmt_size_pct() in R/cluster_data.R
#  - .plot_mmm_posterior() / .plot_covariate_forest() if a method needs them
# ==============================================================================

#' Cluster Diagnostics
#'
#' Unified entry point for clustering quality information. Returns a
#' \code{net_cluster_diagnostics} object that normalises the diagnostic
#' surface across distance-based and model-based clusterings -- you no
#' longer have to know which fields live on \code{net_clustering} vs.
#' \code{net_mmm} vs. the slim \code{net_mmm_clustering} attribute of a
#' \code{netobject_group}.
#'
#' The returned object carries:
#' \describe{
#'   \item{family}{Either \code{"distance"} or \code{"mmm"}.}
#'   \item{k, n, sizes}{Number of clusters, number of sequences, sizes
#'     vector.}
#'   \item{per_cluster}{A \code{data.frame} -- one row per cluster, columns
#'     differ by family. Distance: \code{cluster}, \code{size},
#'     \code{pct}, \code{mean_within_dist}, \code{sil_mean}. MMM:
#'     \code{cluster}, \code{size}, \code{pct}, \code{mix_pct},
#'     \code{avepp}, \code{class_err_pct}.}
#'   \item{overall}{A named list of family-specific summary metrics
#'     (\code{silhouette} for distance; \code{avepp_overall},
#'     \code{entropy}, \code{classification_error} for MMM).}
#'   \item{ics}{For MMM: a list with \code{BIC}, \code{AIC}, \code{ICL}.
#'     \code{NULL} for distance.}
#'   \item{metadata}{Method / dissimilarity / weighted / lambda etc.}
#'   \item{source}{The original clustering object, kept by reference so
#'     \code{plot()} can delegate without recomputing anything.}
#' }
#'
#' @param x A \code{net_clustering}, \code{net_mmm}, \code{netobject_group}
#'   (with \code{attr(, "clustering")} attached by \code{cluster_network()}
#'   or \code{cluster_mmm()}), or \code{net_mmm_clustering}.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#' @return A \code{net_cluster_diagnostics} object.
#' @seealso \code{\link{print.net_cluster_diagnostics}},
#'   \code{\link{plot.net_cluster_diagnostics}},
#'   \code{\link{compare_mmm}} for k-sweep model selection (MMM only).
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' cl <- build_clusters(seqs, k = 2, method = "ward.D2")
#' cluster_diagnostics(cl)
#' \donttest{
#' grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
#' cluster_diagnostics(grp)
#' as.data.frame(cluster_diagnostics(grp))
#' }
#' @export
cluster_diagnostics <- function(x, ...) {
  UseMethod("cluster_diagnostics")
}

#' @export
cluster_diagnostics.default <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  stop("cluster_diagnostics() has no method for class '",
       paste(class(x), collapse = "/"), "'. Supported inputs: ",
       "net_clustering (from build_clusters), net_mmm (from build_mmm), ",
       "netobject_group (from cluster_network / cluster_mmm), ",
       "or net_mmm_clustering (the clustering attribute itself).",
       call. = FALSE)
}

# ---------------------------------------------------------------------------
# Distance-based: net_clustering
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_clustering <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  k <- as.integer(x$k)
  assignments <- as.integer(x$assignments)
  sizes <- as.integer(x$sizes %||% tabulate(assignments, nbins = k))
  n_total <- sum(sizes)

  mean_within <- if (!is.null(x$distance)) {
    .per_cluster_within_dist(x$distance, assignments, k)
  } else {
    rep(NA_real_, k)
  }

  # Per-cluster silhouette mean -- compute once, aggregate. Skip when the
  # cluster package isn't available (it's a Suggests-style dep here, used
  # only inside plot.net_clustering already).
  sil_mean <- rep(NA_real_, k)
  overall_sil <- as.numeric(x$silhouette %||% NA_real_)
  if (!is.null(x$distance) && requireNamespace("cluster", quietly = TRUE)) {
    sil <- tryCatch(
      cluster::silhouette(assignments, dist = x$distance),
      error = function(e) NULL
    )
    if (!is.null(sil) && nrow(sil) > 0L) {
      widths <- sil[, 3L]
      cls    <- sil[, 1L]
      for (cl in seq_len(k)) {
        sel <- cls == cl
        if (any(sel)) sil_mean[cl] <- mean(widths[sel])
      }
      if (is.na(overall_sil)) overall_sil <- mean(widths)
    }
  }

  per_cluster <- data.frame(
    cluster          = seq_len(k),
    size             = sizes,
    pct              = if (n_total > 0L) sizes / n_total * 100 else
                        rep(0, k),
    mean_within_dist = mean_within,
    sil_mean         = sil_mean,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      family      = "distance",
      k           = k,
      n           = as.integer(n_total),
      sizes       = sizes,
      per_cluster = per_cluster,
      overall     = list(silhouette = overall_sil),
      ics         = NULL,
      metadata    = list(
        method        = x$method,
        dissimilarity = x$dissimilarity,
        weighted      = isTRUE(x$weighted),
        lambda        = x$lambda,
        medoids       = x$medoids,
        covariates    = if (is.null(x$covariates)) NULL else
                          setdiff(unique(x$covariates$coefficients$variable),
                                  "(Intercept)")
      ),
      source      = x
    ),
    class = "net_cluster_diagnostics"
  )
}

# ---------------------------------------------------------------------------
# Model-based: net_mmm
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_mmm <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  k <- as.integer(x$k)
  assignments <- as.integer(x$assignments)
  sizes <- as.integer(tabulate(assignments, nbins = k))
  n_total <- as.integer(x$n_sequences %||% sum(sizes))

  # Per-cluster classification-error decomposition: of obs assigned to
  # cluster m, what fraction have max(posterior) < 0.5? This is the
  # missing breakdown of quality$classification_error.
  class_err_pct <- rep(NA_real_, k)
  if (!is.null(x$posterior)) {
    max_post <- apply(x$posterior, 1L, max)
    for (cl in seq_len(k)) {
      members <- which(assignments == cl)
      if (length(members) > 0L) {
        class_err_pct[cl] <- mean(max_post[members] < 0.5) * 100
      }
    }
  }

  avepp <- as.numeric(x$quality$avepp %||% rep(NA_real_, k))
  mix   <- as.numeric(x$mixing %||% rep(NA_real_, k))

  per_cluster <- data.frame(
    cluster       = seq_len(k),
    size          = sizes,
    pct           = if (n_total > 0L) sizes / n_total * 100 else rep(0, k),
    mix_pct       = mix * 100,
    avepp         = avepp,
    class_err_pct = class_err_pct,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      family      = "mmm",
      k           = k,
      n           = n_total,
      sizes       = sizes,
      per_cluster = per_cluster,
      overall     = list(
        avepp_overall        = as.numeric(x$quality$avepp_overall),
        entropy              = as.numeric(x$quality$entropy),
        classification_error = as.numeric(x$quality$classification_error)
      ),
      ics         = list(BIC = as.numeric(x$BIC),
                         AIC = as.numeric(x$AIC),
                         ICL = as.numeric(x$ICL),
                         log_likelihood = as.numeric(x$log_likelihood)),
      metadata    = list(
        states     = x$states,
        n_states   = length(x$states),
        iterations = x$iterations,
        converged  = x$converged,
        covariates = if (is.null(x$covariates)) NULL else
                       setdiff(unique(x$covariates$coefficients$variable),
                               "(Intercept)")
      ),
      source      = x
    ),
    class = "net_cluster_diagnostics"
  )
}

# ---------------------------------------------------------------------------
# net_mmm_clustering -- the slim attribute object on a netobject_group.
# Same field shape as net_mmm minus $models, plus $data. We reuse the
# net_mmm method by re-classing in place (no copy of the heavy slots).
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.net_mmm_clustering <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  shadow <- structure(unclass(x), class = "net_mmm")
  out <- cluster_diagnostics.net_mmm(shadow)
  # Keep the source pointing at the actual input class so plot delegation
  # routes through plot.net_mmm_clustering rather than plot.net_mmm.
  out$source <- x
  out
}

# ---------------------------------------------------------------------------
# netobject_group -- delegate to attr(, "clustering")
# ---------------------------------------------------------------------------

#' @export
cluster_diagnostics.netobject_group <- function(x, ...) {
  .cluster_diag_check_unused_dots(...)
  cl <- attr(x, "clustering")
  if (is.null(cl)) {
    stop("cluster_diagnostics() requires a clustering attribute on the ",
         "netobject_group. Build with cluster_network() or cluster_mmm() ",
         "(or attach an existing net_clustering / net_mmm_clustering as ",
         "attr(grp, \"clustering\")).", call. = FALSE)
  }
  cluster_diagnostics(cl)
}

.cluster_diag_check_unused_dots <- function(...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(
    "cluster_diagnostics() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", "),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# print method
# ---------------------------------------------------------------------------

#' Print Method for net_cluster_diagnostics
#'
#' Prints a uniform header, family-specific quality / IC line, and a
#' per-cluster table. Layout matches \code{\link{print.net_clustering}}
#' and \code{\link{print.net_mmm}}.
#'
#' @param x A \code{net_cluster_diagnostics} object.
#' @param digits Integer. Decimal places for floating-point statistics.
#'   Default \code{3L}.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#' @return The input object, invisibly.
#' @export
print.net_cluster_diagnostics <- function(x, digits = 3L, ...) {
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop("print.net_cluster_diagnostics() got unsupported argument",
         if (length(dots) == 1L) ": " else "s: ",
         paste(dot_names, collapse = ", "), call. = FALSE)
  }
  digits <- as.integer(digits)
  k <- x$k
  n <- x$n

  if (x$family == "distance") {
    cat(sprintf("Cluster Diagnostics (distance) [%s / %s]\n",
                x$metadata$method %||% "?", x$metadata$dissimilarity %||% "?"))
    cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n, k))
    if (!is.null(x$overall$silhouette) && !is.na(x$overall$silhouette)) {
      cat(sprintf("  Quality: silhouette = %.*f\n",
                  digits, x$overall$silhouette))
    }
  } else {
    cat(sprintf("Cluster Diagnostics (mmm) [k = %d]\n", k))
    cat(sprintf("  Sequences: %d  |  Clusters: %d  |  States: %d\n",
                n, k, x$metadata$n_states %||% NA_integer_))
    if (!is.null(x$overall$avepp_overall)) {
      cat(sprintf(
        "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
        digits, x$overall$avepp_overall,
        digits, x$overall$entropy,
        x$overall$classification_error * 100))
    }
    if (!is.null(x$ics)) {
      cat(sprintf(
        "  ICs: LL = %.*f  |  BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
        digits, x$ics$log_likelihood, digits, x$ics$BIC,
        digits, x$ics$AIC,            digits, x$ics$ICL))
    }
  }

  cat("\n")
  pc <- x$per_cluster
  if (x$family == "distance") {
    cols <- list(
      Cluster            = sprintf("%d", pc$cluster),
      N                  = .fmt_size_pct(pc$size, n),
      `Mean within-dist` = ifelse(is.na(pc$mean_within_dist), "--",
                                   sprintf(paste0("%.", digits, "f"),
                                           pc$mean_within_dist)),
      Silhouette         = ifelse(is.na(pc$sil_mean), "--",
                                   sprintf(paste0("%.", digits, "f"),
                                           pc$sil_mean))
    )
  } else {
    cols <- list(
      Cluster      = sprintf("%d", pc$cluster),
      N            = .fmt_size_pct(pc$size, n),
      `Mix%`       = ifelse(is.na(pc$mix_pct), "--",
                             sprintf("%4.1f%%", pc$mix_pct)),
      AvePP        = ifelse(is.na(pc$avepp), "--",
                             sprintf(paste0("%.", digits, "f"), pc$avepp)),
      `Class.Err%` = ifelse(is.na(pc$class_err_pct), "--",
                             sprintf("%4.1f%%", pc$class_err_pct))
    )
  }
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  cov_names <- x$metadata$covariates
  if (!is.null(cov_names) && length(cov_names) > 0L) {
    label <- if (x$family == "distance") "post-hoc" else "integrated"
    cat(sprintf("\n  Covariates: %s (%s, %d predictors)\n",
                paste(cov_names, collapse = ", "), label,
                length(cov_names)))
  }

  invisible(x)
}

# ---------------------------------------------------------------------------
# plot method -- delegates to the source's plot method
# ---------------------------------------------------------------------------

#' Plot Method for net_cluster_diagnostics
#'
#' Delegates to the original clustering object's plot method
#' (\code{\link{plot.net_clustering}} for distance-based diagnostics,
#' \code{\link{plot.net_mmm_clustering}} or \code{\link{plot.net_mmm}}
#' for model-based). The diagnostics object itself stores no plot
#' geometry -- it just keeps a reference to the source so the existing
#' visual layer is reused.
#'
#' @param x A \code{net_cluster_diagnostics} object.
#' @param type Character. Forwarded to the underlying plot method. Valid
#'   values for distance: \code{"silhouette"} (default), \code{"mds"},
#'   \code{"heatmap"}, \code{"predictors"}. Valid values for mmm:
#'   \code{"posterior"} (default), \code{"covariates"} /
#'   \code{"predictors"}.
#' @param ... Forwarded to the underlying plot method.
#' @return A \code{ggplot} object, invisibly.
#' @export
plot.net_cluster_diagnostics <- function(x, type = NULL, ...) {
  if (is.null(type)) {
    type <- if (x$family == "distance") "silhouette" else "posterior"
  }
  plot(x$source, type = type, ...)
}

# ---------------------------------------------------------------------------
# as.data.frame method
# ---------------------------------------------------------------------------

#' @rdname cluster_diagnostics
#' @method as.data.frame net_cluster_diagnostics
#' @param row.names,optional Standard \code{as.data.frame} arguments
#'   (ignored).
#' @export
as.data.frame.net_cluster_diagnostics <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
  dots <- list(...)
  if (length(dots)) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop("as.data.frame.net_cluster_diagnostics() got unsupported argument",
         if (length(dots) == 1L) ": " else "s: ",
         paste(dot_names, collapse = ", "), call. = FALSE)
  }
  x$per_cluster
}
