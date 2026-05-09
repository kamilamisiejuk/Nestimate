# ---- Split-Half Reliability ----

#' Split-Half Reliability for Network Estimates
#'
#' @description
#' Assesses the stability of network estimates by repeatedly splitting
#' sequences into two halves, building networks from each half, and
#' comparing them. Supports single-model reliability assessment and
#' multi-model comparison with optional scaling for cross-method
#' comparability.
#'
#' For transition methods (\code{"relative"}, \code{"frequency"},
#' \code{"co_occurrence"}), uses pre-computed per-sequence count matrices
#' for fast resampling (same infrastructure as
#' \code{\link{bootstrap_network}}).
#'
#' @param ... One or more \code{netobject}s (from \code{\link{build_network}}).
#'   If unnamed, each model is auto-named from its \code{$method}.
#'   A \code{netobject_group} is flattened into its constituent models.
#' @param iter Integer. Number of split-half iterations (default: 1000).
#' @param split Numeric. Fraction of sequences assigned to the first half
#'   (default: 0.5).
#' @param scale Character. Scaling applied to both split-half matrices
#'   before computing metrics. One of \code{"none"} (default),
#'   \code{"minmax"}, \code{"standardize"}, or \code{"proportion"}.
#'   Use scaling when comparing models on different scales (e.g. frequency
#'   vs relative).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_reliability"} containing:
#' \describe{
#'   \item{iterations}{Data frame with columns \code{model}, \code{mean_dev},
#'     \code{median_dev}, \code{cor}, \code{max_dev} (one row per iteration
#'     per model).}
#'   \item{summary}{Data frame with columns \code{model}, \code{metric},
#'     \code{mean}, \code{sd}.}
#'   \item{models}{Named list of the original \code{netobject}s.}
#'   \item{iter}{Number of iterations.}
#'   \item{split}{Split fraction.}
#'   \item{scale}{Scaling method used.}
#' }
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' rel <- network_reliability(net, iter = 10)
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rel <- network_reliability(net, iter = 100, seed = 42)
#' print(rel)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}}
#'
#' @importFrom stats cor median sd
#' @export
network_reliability <- function(..., iter = 1000L, split = 0.5,
                        scale = "none", seed = NULL) {

  dots <- list(...)

  # ---- Input validation ----
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(split), length(split) == 1, split > 0, split < 1,
    is.character(scale), length(scale) == 1
  )
  iter <- as.integer(iter)
  scale <- match.arg(scale, c("none", "minmax", "standardize", "proportion"))

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Flatten netobject_group elements ----
  objs <- list()
  labels <- character(0)
  for (i in seq_along(dots)) {
    obj <- dots[[i]]
    nm <- names(dots)[i]
    if (inherits(obj, "mcml")) obj <- as_tna(obj)
    if (inherits(obj, "cograph_network")) obj <- .as_netobject(obj)
    if (inherits(obj, "netobject_group")) {
      for (g in names(obj)) {
        objs <- c(objs, list(obj[[g]]))
        labels <- c(labels, g)
      }
    } else if (inherits(obj, "netobject")) {
      label <- if (!is.null(nm) && nzchar(nm)) nm else obj$method
      objs <- c(objs, list(obj))
      labels <- c(labels, label)
    } else {
      stop("All arguments must be netobject or netobject_group objects.",
           call. = FALSE)
    }
  }

  if (length(objs) == 0L) {
    stop("At least one netobject is required.", call. = FALSE)
  }

  # Deduplicate names
  labels <- make.unique(labels, sep = "_")
  model_list <- stats::setNames(objs, labels)

  # ---- Warn if different methods without scaling ----
  methods <- vapply(model_list, function(m) m$method, character(1))
  if (length(unique(methods)) > 1L && scale == "none") {
    warning(
      "Models use different methods (",
      paste(unique(methods), collapse = ", "),
      "). Consider setting scale = 'minmax' for comparable results.",
      call. = FALSE
    )
  }

  # ---- Run split-half per model ----
  all_iters <- lapply(names(model_list), function(model_name) {
    net <- model_list[[model_name]]

    if (is.null(net$data)) {
      stop("netobject '", model_name,
           "' does not contain $data. Rebuild with build_network().",
           call. = FALSE)
    }

    method <- .resolve_method_alias(net$method)
    is_transition <- method %in% c("relative", "frequency", "co_occurrence")

    if (is_transition) {
      .reliability_transition(net, method, iter, split, scale)
    } else {
      .reliability_association(net, method, iter, split, scale)
    }
  })
  names(all_iters) <- names(model_list)

  # ---- Assemble iterations data frame ----
  iterations <- do.call(rbind, lapply(names(all_iters), function(nm) {
    df <- all_iters[[nm]]
    df$model <- nm
    df[, c("model", "mean_dev", "median_dev", "cor", "max_dev")]
  }))
  rownames(iterations) <- NULL

  # ---- Compute summary ----
  metric_names <- c("mean_dev", "median_dev", "cor", "max_dev")
  summary_rows <- do.call(rbind, lapply(names(model_list), function(nm) {
    sub <- iterations[iterations$model == nm, , drop = FALSE]
    do.call(rbind, lapply(metric_names, function(met) {
      vals <- sub[[met]]
      data.frame(
        model = nm,
        metric = met,
        mean = mean(vals, na.rm = TRUE),
        sd = sd(vals, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
  rownames(summary_rows) <- NULL

  result <- list(
    iterations = iterations,
    summary = summary_rows,
    models = model_list,
    iter = iter,
    split = split,
    scale = scale
  )
  class(result) <- "net_reliability"
  result
}


# ---- Transition fast path ----

#' Split-half reliability for transition methods
#' @noRd
.reliability_transition <- function(net, method, iter, split, scale) {
  states <- net$nodes$label
  n_states <- length(states)
  nbins <- n_states * n_states
  is_relative <- method == "relative"
  scaling <- net$scaling
  threshold <- net$threshold

  # Pre-compute per-sequence count matrix (reuse bootstrap infrastructure)
  trans_2d <- .precompute_per_sequence(net$data, method, net$params, states)
  n_seq <- nrow(trans_2d)
  n_half <- max(1L, round(n_seq * split))

  # Post-process raw counts into the appropriate matrix form
  postprocess <- function(counts) {
    mat <- matrix(counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov
    if (threshold > 0) mat[abs(mat) < threshold] <- 0 # nocov
    mat
  }

  # Run iterations
  results <- vapply(seq_len(iter), function(i) {
    idx_a <- sample.int(n_seq, n_half, replace = FALSE)
    idx_b <- setdiff(seq_len(n_seq), idx_a)

    counts_a <- colSums(trans_2d[idx_a, , drop = FALSE])
    counts_b <- colSums(trans_2d[idx_b, , drop = FALSE])

    mat_a <- postprocess(counts_a)
    mat_b <- postprocess(counts_b)

    # Apply reliability scaling
    mat_a <- .scale_matrix(mat_a, scale)
    mat_b <- .scale_matrix(mat_b, scale)

    .split_half_metrics(mat_a, mat_b)
  }, numeric(4))

  data.frame(
    mean_dev = results[1, ],
    median_dev = results[2, ],
    cor = results[3, ],
    max_dev = results[4, ],
    stringsAsFactors = FALSE
  )
}


# ---- Association path ----

#' Split-half reliability for association methods
#' @noRd
.reliability_association <- function(net, method, iter, split, scale) {
  data <- net$data
  states <- net$nodes$label
  n_states <- length(states)
  nbins <- n_states * n_states
  scaling <- net$scaling
  threshold <- net$threshold
  params <- net$params
  level <- net$level
  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")

  estimator <- get_estimator(method)
  n <- nrow(data)
  n_half <- max(1L, round(n * split))

  postprocess <- function(mat) {
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov start
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    mat # nocov end
  }

  results <- vapply(seq_len(iter), function(i) {
    idx_a <- sample.int(n, n_half, replace = FALSE)
    idx_b <- setdiff(seq_len(n), idx_a)

    build_half <- function(sub_data) {
      if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
        sub_data <- tryCatch( # nocov start
          .decompose_multilevel(sub_data, id_col = id_col, level = level),
          error = function(e) NULL
        )
        if (is.null(sub_data)) return(NULL) # nocov end
      }
      tryCatch(
        do.call(estimator$fn, c(list(data = sub_data), params)),
        error = function(e) NULL
      )
    }

    est_a <- build_half(data[idx_a, , drop = FALSE])
    est_b <- build_half(data[idx_b, , drop = FALSE])

    if (is.null(est_a) || is.null(est_b)) return(rep(NA_real_, 4))

    mat_a <- est_a$matrix[states, states] # nocov start
    mat_b <- est_b$matrix[states, states]

    mat_a <- postprocess(mat_a)
    mat_b <- postprocess(mat_b)

    mat_a <- .scale_matrix(mat_a, scale)
    mat_b <- .scale_matrix(mat_b, scale)

    .split_half_metrics(mat_a, mat_b) # nocov end
  }, numeric(4))

  data.frame(
    mean_dev = results[1, ],
    median_dev = results[2, ],
    cor = results[3, ],
    max_dev = results[4, ],
    stringsAsFactors = FALSE
  )
}


# ---- Helpers ----

#' Compute split-half metrics between two matrices
#'
#' Thin wrapper around `.network_similarity()` (R/compare_network.R) that
#' returns the four metrics the reliability summary table reports, in the
#' fixed order `c(mean_dev, median_dev, cor, max_dev)`.
#' @noRd
.split_half_metrics <- function(mat_a, mat_b) {
  unname(.network_similarity(mat_a, mat_b,
                             metrics = c("mean_abs_diff", "median_abs_diff",
                                         "pearson", "max_abs_diff")))
}


#' Scale a matrix for cross-method comparison
#' @noRd
.scale_matrix <- function(mat, method) {
  switch(method,
    none = mat,
    minmax = {
      rng <- range(mat)
      if (rng[1] == rng[2]) mat
      else (mat - rng[1]) / (rng[2] - rng[1])
    },
    standardize = {
      s <- sd(as.vector(mat))
      if (s == 0) mat
      else (mat - mean(mat)) / s
    },
    proportion = {
      total <- sum(mat)
      if (total == 0) mat
      else mat / total
    }
  )
}


# ---- S3 Methods ----

#' Print Method for net_reliability
#'
#' @param x A \code{net_reliability} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' rel <- network_reliability(net, iter = 10)
#' print(rel)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rel <- network_reliability(net, iter = 20, seed = 1)
#' print(rel)
#' }
#'
#' @export
print.net_reliability <- function(x, ...) {
  cat(sprintf("Split-Half Reliability (%d iterations, split = %.0f%%",
              x$iter, x$split * 100))
  if (x$scale != "none") {
    cat(sprintf(", scale = %s", x$scale))
  }
  cat(")\n")

  models <- unique(x$summary$model)
  metric_labels <- c(
    mean_dev = "Mean Abs. Diff.",
    median_dev = "Median Abs. Diff.",
    cor = "Pearson",
    max_dev = "Max Abs. Diff."
  )

  for (m in models) {
    if (length(models) > 1L) cat(sprintf("\n  %s:\n", m))
    sub <- x$summary[x$summary$model == m, , drop = FALSE]
    for (i in seq_len(nrow(sub))) {
      label <- metric_labels[sub$metric[i]]
      prefix <- if (length(models) > 1L) "    " else "  "
      cat(sprintf("%s%-18s  mean = %.4f  sd = %.4f\n",
                  prefix, label, sub$mean[i], sub$sd[i]))
    }
  }

  invisible(x)
}


#' Summary Method for net_reliability
#'
#' @param object A \code{net_reliability} object.
#' @param ... Ignored.
#' @return A tidy data frame with columns \code{model}, \code{metric},
#'   \code{mean}, \code{sd} summarising the split-half iterations.
#' @inherit network_reliability examples
#' @export
summary.net_reliability <- function(object, ...) {
  object$summary
}


#' Plot Method for net_reliability
#'
#' @description
#' Density plots of split-half metrics faceted by metric type.
#' Multi-model comparisons show overlaid densities colored by model.
#'
#' @param x A \code{net_reliability} object.
#' @param bins Integer. Number of histogram bins per panel (default 60).
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C","A"),
#'   V2 = c("B","C","A","B")), method = "relative")
#' rel <- network_reliability(net, iter = 10)
#' plot(rel)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rel <- network_reliability(net, iter = 20, seed = 1)
#' plot(rel)
#' }
#'
#' @export
plot.net_reliability <- function(x, bins = 60L, ...) {
  iters <- x$iterations
  models <- unique(iters$model)
  multi <- length(models) > 1L

  metric_labels <- c(
    max_dev = "Max Abs. Diff.",
    mean_dev = "Mean Abs. Diff.",
    median_dev = "Median Abs. Diff.",
    cor = "Pearson"
  )

  metric_cols <- names(metric_labels)
  long <- do.call(rbind, lapply(metric_cols, function(met) {
    data.frame(
      model = iters$model,
      metric = metric_labels[met],
      value = iters[[met]],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  long$metric <- factor(long$metric, levels = metric_labels)

  means <- aggregate(value ~ model + metric, data = long, FUN = mean)
  means$label <- sprintf("%.2f", means$value)

  palette <- c("#4DA167", "#E89AB4", "#2B6CB0", "#DD8452",
               "#8172B2", "#937860", "#DA8BC3", "#8C8C8C")
  fill_vals <- palette[seq_along(models)]
  names(fill_vals) <- models

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_histogram(
      ggplot2::aes(fill = .data$model),
      bins = bins, alpha = 0.65, position = "identity",
      color = NA
    ) +
    ggplot2::geom_vline(
      data = means,
      ggplot2::aes(xintercept = .data$value, color = .data$model),
      linetype = "dashed", linewidth = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = means,
      ggplot2::aes(x = .data$value, y = Inf, label = .data$label,
                   color = .data$model),
      vjust = 1.2, size = 3.3, label.size = 0.4,
      fill = "white", show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, name = "Model") +
    ggplot2::scale_color_manual(values = fill_vals, guide = "none") +
    ggplot2::facet_wrap(~ metric, scales = "free", ncol = 2L) +
    ggplot2::labs(x = "Metric Value", y = "Frequency",
                  title = if (multi) "Split-Half Reliability" else NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = if (multi) "bottom" else "none",
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      panel.grid.minor = ggplot2::element_blank()
    )

  print(p)
  invisible(p)
}
