# ---- Edge-weight case-dropping stability --------------------------------
# Mirrors centrality_stability() but the stability target is the vector of
# edge weights (flattened adjacency) rather than per-node centralities.
# Answers: "How many cases can be dropped before the full set of edge
# weights loses rank correlation with the original?"

#' Edge-weight Case-dropping Stability
#'
#' Computes a **CS-coefficient for the edge-weight vector** of a network:
#' the maximum proportion of cases (rows of `x$data`) that can be dropped
#' while the flattened edge-weight vector of the re-estimated network
#' still correlates with the original above `threshold` in at least
#' `certainty` of iterations.
#'
#' Complements [centrality_stability()]: that function asks whether
#' centrality *rankings* are stable; this one asks whether the *edge-weight
#' structure itself* is stable. For MCML-derived networks where each row
#' of `$data` is one transition, this is case-dropping of **edges**.
#'
#' @param x A `netobject`, `cograph_network`, `netobject_group`, or `mcml`.
#'   For the group types this function iterates over each element and
#'   returns a named list.
#' @param iter Integer. Iterations per drop proportion. Default `1000`.
#' @param drop_prop Numeric vector of proportions to evaluate. Each entry
#'   must lie strictly between 0 and 1. Default `seq(0.1, 0.9, by = 0.1)`.
#' @param threshold Numeric in `[0, 1]`. Minimum edge-vector correlation
#'   for an iteration to count as stable. Default `0.7`.
#' @param certainty Numeric in `[0, 1]`. Required fraction of iterations
#'   whose correlation must exceed `threshold` for a drop proportion to
#'   qualify. Default `0.95`.
#' @param method Correlation method: `"pearson"` (weight magnitudes),
#'   `"spearman"` (ranks, robust to scale), or `"kendall"`. Default
#'   `"spearman"` because edge weights often span several orders of
#'   magnitude and rank stability is the typical target.
#' @param include_diag Logical. Include diagonal (self-loop) edges in the
#'   edge vector. Default `FALSE`.
#' @param seed Optional integer for reproducibility.
#'
#' @return An object of class `net_casedrop_reliability` with:
#' \describe{
#'   \item{`cs`}{Scalar CS-coefficient — the maximum drop proportion for
#'     which the edge-vector correlation remains >= `threshold` in at
#'     least `certainty` of iterations. Zero if no proportion qualifies.}
#'   \item{`correlations`}{`iter` x `length(drop_prop)` matrix of per-
#'     iteration correlations.}
#'   \item{`drop_prop`, `threshold`, `certainty`, `iter`, `method`}{Inputs.}
#' }
#'
#' @details
#' For each `drop_prop` p and each iteration, a size `n_cases * (1 - p)`
#' subset of `$data` rows is selected **without replacement**, the network
#' is re-estimated using the same method/scaling/threshold as the input,
#' and the upper/lower-triangle (directed: all off-diagonal entries) of
#' the new weight matrix is flattened and correlated with the
#' corresponding vector of the original matrix. The correlation method
#' defaults to Spearman for robustness to the wide dynamic range of
#' transition probabilities.
#'
#' Unlike bootstrap CIs, case-dropping does not estimate sampling variance
#' and so does not rely on the i.i.d. assumption. This makes it the
#' appropriate robustness check for **edgelist-derived** networks (where
#' rows of `$data` lack actor grouping), since dropping rows at random is
#' a well-posed operation regardless of within-actor correlation.
#'
#' @seealso [centrality_stability()], [bootstrap_network()].
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE),
#'   V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' es  <- casedrop_reliability(net, iter = 50, drop_prop = c(0.1, 0.3, 0.5),
#'                       seed = 1)
#' print(es)
#'
#' @references
#' Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#' psychological networks and their accuracy: A tutorial paper.
#' \emph{Behavior Research Methods} 50(1), 195-212.
#' \doi{10.3758/s13428-017-0862-1}
#'
#' @export
casedrop_reliability <- function(x,
                           iter        = 1000L,
                           drop_prop   = seq(0.1, 0.9, by = 0.1),
                           threshold   = 0.7,
                           certainty   = 0.95,
                           method      = c("spearman", "pearson", "kendall"),
                           include_diag = FALSE,
                           seed        = NULL) {

  # ---- Dispatch for mcml / group inputs ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      casedrop_reliability(net, iter = iter, drop_prop = drop_prop,
                     threshold = threshold, certainty = certainty,
                     method = method, include_diag = include_diag,
                     seed = seed)
    })
    class(out) <- c("net_casedrop_reliability_group", "list")
    return(out)
  }
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }
  if (!inherits(x, "netobject")) {
    stop("'x' must be a netobject from build_network() (or an mcml / ",
         "netobject_group that contains one).", call. = FALSE)
  }
  if (is.null(x$data)) {
    stop("netobject does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  # ---- Input validation ----
  stopifnot(
    is.numeric(iter), length(iter) == 1L, iter >= 2,
    is.numeric(drop_prop), all(drop_prop > 0), all(drop_prop < 1),
    is.numeric(threshold), length(threshold) == 1L,
    threshold >= 0, threshold <= 1,
    is.numeric(certainty), length(certainty) == 1L,
    certainty >= 0, certainty <= 1,
    is.logical(include_diag), length(include_diag) == 1L
  )
  iter <- as.integer(iter)
  method <- match.arg(method)
  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1L)
    set.seed(seed)
  }

  # ---- Resolve estimator + reuse bootstrap precompute for transition methods ----
  net_method <- .resolve_method_alias(x$method)
  states     <- x$nodes$label
  n_states   <- length(states)
  is_transition <- net_method %in% c("relative", "frequency", "co_occurrence")
  is_relative   <- net_method == "relative"
  scaling <- x$scaling
  thresh  <- x$threshold

  # Original edge vector (flattened)
  idx_mask <- if (include_diag) {
    matrix(TRUE, n_states, n_states)
  } else {
    diag_mask <- matrix(TRUE, n_states, n_states)
    diag(diag_mask) <- FALSE
    diag_mask
  }
  orig_edges <- as.vector(x$weights[idx_mask])
  if (stats::sd(orig_edges, na.rm = TRUE) == 0) {
    warning("Original edge vector has zero variance; CS undefined.",
            call. = FALSE)
    empty <- matrix(NA_real_, iter, length(drop_prop))
    result <- list(
      cs            = 0,
      summary       = data.frame(),
      metrics       = list(mean_abs_dev = empty, median_abs_dev = empty,
                           correlation = empty,  max_abs_dev = empty),
      correlations  = empty,
      drop_prop     = drop_prop,
      threshold     = threshold,
      certainty     = certainty,
      iter          = iter,
      method        = method,
      include_diag  = include_diag,
      n_cases       = nrow(x$data),
      n_edges       = length(orig_edges)
    )
    class(result) <- "net_casedrop_reliability"
    return(result)
  }

  # ---- Build-matrix closure: transition fast path or association slow path ----
  if (is_transition) {
    trans_2d <- .precompute_per_sequence(x$data, net_method, x$params, states)
    n_cases  <- nrow(trans_2d)

    build_matrix <- function(idx) {
      counts <- colSums(trans_2d[idx, , drop = FALSE])
      mat <- matrix(counts, n_states, n_states, byrow = TRUE)
      if (is_relative) {
        rs <- rowSums(mat)
        nz <- rs > 0
        mat[nz, ] <- mat[nz, ] / rs[nz]
      }
      if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
      if (thresh > 0) mat[abs(mat) < thresh] <- 0
      mat
    }
  } else {
    data      <- x$data
    n_cases   <- nrow(data)
    estimator <- get_estimator(net_method)
    params    <- x$params
    level     <- x$level
    id_col    <- .param_get(x$params, "id") %||%
      .param_get(x$params, "id_col")

    build_matrix <- function(idx) {
      sub <- data[idx, , drop = FALSE]
      if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
        sub <- tryCatch(
          .decompose_multilevel(sub, id_col = id_col, level = level),
          error = function(e) NULL
        )
        if (is.null(sub)) return(NULL)
      }
      est <- tryCatch(
        do.call(estimator$fn, c(list(data = sub), params)),
        error = function(e) NULL
      )
      if (is.null(est)) return(NULL)
      mat <- est$matrix[states, states]
      if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
      if (thresh > 0) mat[abs(mat) < thresh] <- 0
      mat
    }
  }

  # ---- Main case-dropping loop ----
  n_prop <- length(drop_prop)
  # Per-iteration, per-drop_prop storage for the four model-level metrics
  metric_names <- c("mean_abs_dev", "median_abs_dev", "correlation",
                    "max_abs_dev")
  metrics <- lapply(metric_names, function(nm) {
    matrix(NA_real_, nrow = iter, ncol = n_prop,
           dimnames = list(NULL, paste0("p", drop_prop)))
  })
  names(metrics) <- metric_names

  case_seq <- seq_len(n_cases)
  for (p_idx in seq_len(n_prop)) {
    n_drop <- floor(n_cases * drop_prop[p_idx])
    if (n_drop == 0L || n_drop >= n_cases) next
    keep_n <- n_cases - n_drop

    cor_metric <- match.arg(method, c("pearson", "spearman", "kendall"))
    for (it in seq_len(iter)) {
      idx <- sample(case_seq, keep_n, replace = FALSE)
      mat <- build_matrix(idx)
      if (is.null(mat)) next
      sub_edges <- as.vector(mat[idx_mask])
      sim <- .network_similarity(orig_edges, sub_edges,
                                 metrics = c("mean_abs_diff",
                                             "median_abs_diff",
                                             cor_metric,
                                             "max_abs_diff"))
      metrics$mean_abs_dev  [it, p_idx] <- sim[["mean_abs_diff"]]
      metrics$median_abs_dev[it, p_idx] <- sim[["median_abs_diff"]]
      metrics$correlation   [it, p_idx] <- sim[[cor_metric]]
      metrics$max_abs_dev   [it, p_idx] <- sim[["max_abs_diff"]]
    }
  }

  # ---- Per-drop-proportion summary across iterations ----
  summary_rows <- lapply(metric_names, function(nm) {
    m <- metrics[[nm]]
    data.frame(
      metric    = nm,
      drop_prop = drop_prop,
      mean      = colMeans(m, na.rm = TRUE),
      sd        = apply(m, 2L, stats::sd,       na.rm = TRUE),
      median    = apply(m, 2L, stats::median,   na.rm = TRUE),
      mad       = apply(m, 2L, stats::mad,      na.rm = TRUE),
      q025      = apply(m, 2L, stats::quantile,
                        probs = 0.025, na.rm = TRUE),
      q975      = apply(m, 2L, stats::quantile,
                        probs = 0.975, na.rm = TRUE),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  summary_df <- do.call(rbind, summary_rows)

  # ---- CS-coefficient: max drop_prop where corr certainty-fraction >= threshold ----
  prop_above <- apply(metrics$correlation, 2L, function(c) {
    mean(c >= threshold, na.rm = TRUE)
  })
  qualifying <- which(prop_above >= certainty)
  cs <- if (length(qualifying) > 0L) max(drop_prop[qualifying]) else 0

  structure(
    list(
      cs            = cs,
      summary       = summary_df,
      metrics       = metrics,
      correlations  = metrics$correlation,
      drop_prop     = drop_prop,
      threshold     = threshold,
      certainty     = certainty,
      iter          = iter,
      method        = method,
      include_diag  = include_diag,
      n_cases       = n_cases,
      n_edges       = length(orig_edges)
    ),
    class = "net_casedrop_reliability"
  )
}

#' @param x An edge-stability object.
#' @param digits Digits to display. Default `3`.
#' @param ... Additional arguments (ignored).
#' @return The input `x` invisibly.
#' @rdname casedrop_reliability
#' @export
print.net_casedrop_reliability <- function(x, digits = 3, ...) {
  cat(sprintf("Edge-weight Case-dropping Stability\n"))
  cat(sprintf("  Cases (rows of $data) : %d\n", x$n_cases))
  cat(sprintf("  Edges assessed        : %d%s\n", x$n_edges,
              if (!x$include_diag) " (diagonal excluded)" else ""))
  cat(sprintf("  Iterations / prop     : %d\n", x$iter))
  cat(sprintf("  Correlation method    : %s\n", x$method))
  cat(sprintf("  CS-coefficient (r)    : %.2f  (threshold=%.2f, certainty=%.2f)\n",
              x$cs, x$threshold, x$certainty))
  cat("\nModel-level reliability across iterations (mean +/- sd per drop):\n")
  # Pivot summary into a compact matrix for display
  props  <- x$drop_prop
  fmt <- function(v) formatC(v, digits = digits, format = "f", flag = " ")
  show <- function(metric_name, display) {
    rows <- x$summary[x$summary$metric == metric_name, ]
    cat(sprintf("  %-14s ", display))
    cat(paste(sprintf("%s+-%s",
                      fmt(rows$mean), fmt(rows$sd)),
              collapse = "  "))
    cat("\n")
  }
  cat(sprintf("  %-14s %s\n", "drop_prop",
              paste(sprintf("%-11s", sprintf("p=%.1f", props)),
                    collapse = "  ")))
  show("mean_abs_dev",   "mean|diff|")
  show("median_abs_dev", "MAD")
  show("correlation",    "cor")
  show("max_abs_dev",    "max|diff|")
  invisible(x)
}

#' Summary method for net_casedrop_reliability
#'
#' @param object A `net_casedrop_reliability`.
#' @param ... Additional arguments (ignored).
#' @return A tidy data frame with columns \code{metric}, \code{drop_prop},
#'   \code{mean}, \code{sd} summarising edge-weight stability across
#'   case-dropping iterations.
#' @rdname casedrop_reliability
#' @export
summary.net_casedrop_reliability <- function(object, ...) {
  object$summary
}

#' @rdname casedrop_reliability
#' @export
print.net_casedrop_reliability_group <- function(x, ...) {
  cs        <- vapply(x, function(e) e$cs,       numeric(1))
  n_edges   <- vapply(x, function(e) e$n_edges,  integer(1))
  n_cases   <- vapply(x, function(e) e$n_cases,  integer(1))
  out <- data.frame(
    n_cases = n_cases,
    n_edges = n_edges,
    CS      = round(cs, 2),
    row.names = names(x)
  )
  cat(sprintf("Edge-weight Case-dropping Stability (%d networks, threshold = %.2f)\n",
              length(x), x[[1]]$threshold))
  print(out)
  invisible(x)
}

#' Summary method for net_casedrop_reliability_group
#'
#' @param object A `net_casedrop_reliability_group`.
#' @param drop_prop Drop proportion at which to report the four metrics
#'   (mean +/- sd per network). Default `0.7`.
#' @param ... Additional arguments (ignored).
#' @return A data frame with one row per network containing
#'   `cor`, `mean_abs_dev`, `median_abs_dev`, `max_abs_dev` formatted as
#'   "mean +/- sd".
#' @rdname casedrop_reliability
#' @export
summary.net_casedrop_reliability_group <- function(object, drop_prop = 0.7, ...) {
  fmt <- function(m, s) sprintf("%.3f +/- %.3f", m, s)
  one_row <- function(e) {
    s <- e$summary
    r <- s[abs(s$drop_prop - drop_prop) < 1e-9, , drop = FALSE]
    stats::setNames(
      list(fmt(r$mean[r$metric == "correlation"],    r$sd[r$metric == "correlation"]),
           fmt(r$mean[r$metric == "mean_abs_dev"],   r$sd[r$metric == "mean_abs_dev"]),
           fmt(r$mean[r$metric == "median_abs_dev"], r$sd[r$metric == "median_abs_dev"]),
           fmt(r$mean[r$metric == "max_abs_dev"],    r$sd[r$metric == "max_abs_dev"])),
      c("cor", "mean_abs_dev", "median_abs_dev", "max_abs_dev")
    )
  }
  rows <- lapply(object, one_row)
  out  <- do.call(rbind.data.frame, rows)
  out  <- cbind(
    n_edges = vapply(object, function(e) e$n_edges, integer(1)),
    out
  )
  rownames(out) <- names(object)
  structure(out, class = c("summary.net_casedrop_reliability_group", "data.frame"),
            drop_prop = drop_prop)
}

#' @param x A `summary.net_casedrop_reliability_group` object.
#' @param ... Additional arguments (ignored).
#' @rdname casedrop_reliability
#' @export
print.summary.net_casedrop_reliability_group <- function(x, ...) {
  cat(sprintf("Edge-weight reliability at drop = %.2f (mean +/- sd over iterations)\n",
              attr(x, "drop_prop")))
  print.data.frame(x)
  invisible(x)
}

#' Plot method for edge-stability result
#'
#' Plots the four model-level reliability metrics across drop
#' proportions: `correlation`, `mean_abs_dev`, `median_abs_dev`,
#' `max_abs_dev`. Each panel shows the per-iteration mean with a ribbon
#' at mean +/- sd. The `correlation` panel includes a dashed horizontal
#' line at the user's `threshold` (default 0.7).
#'
#' @param x A `net_casedrop_reliability` object from [casedrop_reliability()].
#' @param ... Additional arguments (ignored).
#'
#' @return A `ggplot` object.
#' @rdname casedrop_reliability
#' @export
plot.net_casedrop_reliability <- function(x, ...) {
  df <- x$summary
  df$metric <- factor(df$metric,
                      levels = c("correlation", "mean_abs_dev",
                                 "median_abs_dev", "max_abs_dev"),
                      labels = c("Correlation",
                                 "Mean |diff|",
                                 "Median |diff| (MAD)",
                                 "Max |diff|"))

  ggplot2::ggplot(df,
                  ggplot2::aes(x = .data$drop_prop,
                               y = .data$mean)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$mean - .data$sd,
                   ymax = .data$mean + .data$sd),
      fill = "#2196F3", alpha = 0.2
    ) +
    ggplot2::geom_line(color = "#0D47A1", linewidth = 0.8) +
    ggplot2::geom_point(color = "#0D47A1", size = 1.5) +
    ggplot2::geom_hline(
      data = data.frame(metric = factor("Correlation",
                                         levels = levels(df$metric)),
                        yint = x$threshold),
      ggplot2::aes(yintercept = .data$yint),
      linetype = "dashed", color = "grey40"
    ) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y") +
    ggplot2::labs(
      x = "Proportion of cases dropped",
      y = "Mean \u00b1 SD across iterations",
      title = sprintf("Edge-weight reliability (CS = %.2f, %d iter, %d edges)",
                      x$cs, x$iter, x$n_edges),
      subtitle = sprintf("Case-dropping on %d cases, %s correlation",
                         x$n_cases, x$method)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      strip.text    = ggplot2::element_text(face = "bold")
    )
}

#' Plot method for grouped edge-stability result
#'
#' Overlay of per-cluster correlation curves across drop proportions.
#' One colour per sub-network; ribbons show mean +/- sd across
#' iterations. Dashed horizontal line marks the stability threshold
#' (default 0.7).
#'
#' @param x A `net_casedrop_reliability_group` object.
#' @param metric Which metric to plot. One of `"correlation"`
#'   (default), `"mean_abs_dev"`, `"median_abs_dev"`, `"max_abs_dev"`.
#' @param ... Additional arguments (ignored).
#'
#' @return A `ggplot` object.
#' @rdname casedrop_reliability
#' @export
plot.net_casedrop_reliability_group <- function(x,
                                          metric = c("correlation",
                                                     "mean_abs_dev",
                                                     "median_abs_dev",
                                                     "max_abs_dev"),
                                          ...) {
  metric <- match.arg(metric)

  rows <- lapply(names(x), function(nm) {
    s <- x[[nm]]$summary
    s <- s[s$metric == metric, , drop = FALSE]
    s$network <- nm
    s$cs <- x[[nm]]$cs
    s
  })
  df <- do.call(rbind, rows)
  df$network <- factor(df$network, levels = names(x))

  # Use a discrete palette, recycling if there are many networks
  palette <- grDevices::hcl.colors(length(levels(df$network)), "Dynamic")

  p <- ggplot2::ggplot(df,
                       ggplot2::aes(x = .data$drop_prop, y = .data$mean,
                                    colour = .data$network,
                                    fill   = .data$network)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$mean - .data$sd,
                   ymax = .data$mean + .data$sd),
      alpha = 0.15, colour = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::scale_fill_manual(values = palette, drop = FALSE) +
    ggplot2::labs(
      x = "Proportion of cases dropped",
      y = switch(metric,
                 correlation    = "Edge-vector correlation",
                 mean_abs_dev   = "Mean |diff|",
                 median_abs_dev = "Median |diff| (MAD)",
                 max_abs_dev    = "Max |diff|"),
      colour = "Network", fill = "Network",
      title = "Edge-weight case-dropping reliability",
      subtitle = sprintf("Metric: %s; mean \u00b1 sd across iterations",
                         metric)
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      legend.position = "right"
    )

  if (metric == "correlation") {
    threshold <- x[[1]]$threshold
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                                 linetype = "dashed", color = "grey40")
  }

  p
}
