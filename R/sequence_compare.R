# ---- Sequence Pattern Comparison ----

utils::globalVariables(c(
  "freq", "p_label", "pattern", "resid", "resid_label",
  "side", "x_center", "xmax", "xmin", "yidx", "ymax", "ymin"
))

#' Compare Subsequence Patterns Between Groups
#'
#' @description
#' Extracts all k-gram patterns (subsequences of length k) from sequences
#' in each group, computes standardized residuals against the independence
#' model, and optionally runs a permutation or chi-square test of
#' group differences.
#'
#' @details
#' Standardized residuals are always computed from a 2xG contingency table
#' of (this pattern vs. everything else) using the textbook formula
#' \code{(o - e) / sqrt(e * (1 - r/N) * (1 - c/N))}. They describe how much
#' each group's count for a given pattern deviates from expectation under
#' independence, scaled to be approximately N(0,1) under the null.
#'
#' The optional \code{test} argument chooses an inference method:
#' \describe{
#'   \item{\code{"permutation"}}{Shuffles group labels across sequences and
#'     recomputes a per-pattern statistic (row-wise Euclidean residual norm).
#'     Answers: "is this pattern's distribution associated with group
#'     membership at the \emph{actor} level?" Respects the sequence as the
#'     unit of analysis; can be underpowered when the number of sequences
#'     is small.}
#'   \item{\code{"chisq"}}{Runs \code{chisq.test} on the 2xG table per
#'     pattern. Answers: "do the group \emph{streams} generate this pattern
#'     at different rates?" Treats each k-gram occurrence as an event; fast
#'     and powerful even with few sequences, but the iid assumption it makes
#'     is optimistic when sequences are strongly autocorrelated.}
#'   \item{\code{"none"}}{Skip inference. Only residuals, frequencies, and
#'     proportions are returned.}
#' }
#'
#' P-values are adjusted once across all patterns (not per-pattern) using
#' any method supported by \code{\link[stats]{p.adjust}}. The default is
#' \code{"fdr"} (Benjamini-Hochberg).
#'
#' @param x A \code{netobject_group} (from grouped \code{build_network}),
#'   a \code{netobject} (requires \code{group}), or a wide-format
#'   \code{data.frame} (requires \code{group}).
#' @param group Character or vector. Column name or vector of group labels.
#'   Not needed for \code{netobject_group}.
#' @param sub Integer vector. Pattern lengths to analyze. Default: \code{3:5}.
#' @param min_freq Integer. Minimum frequency in each group for a pattern
#'   to be included. Default: 5.
#' @param test Character. Inference method: one of \code{"permutation"}
#'   (default), \code{"chisq"}, or \code{"none"}. See Details.
#' @param iter Integer. Permutation iterations. Only used when
#'   \code{test = "permutation"}. Default: 1000.
#' @param adjust Character. P-value correction method (see
#'   \code{\link[stats]{p.adjust}}). Default: \code{"fdr"}.
#'
#' @return An object of class \code{"net_sequence_comparison"} containing:
#' \describe{
#'   \item{patterns}{Tidy data.frame. Always present:
#'     \code{pattern}, \code{length}, \code{freq_<group>},
#'     \code{prop_<group>}, \code{resid_<group>}. If
#'     \code{test = "permutation"}: \code{effect_size}, \code{p_value}.
#'     If \code{test = "chisq"}: \code{statistic}, \code{p_value}.}
#'   \item{groups}{Character vector of group names.}
#'   \item{n_patterns}{Integer. Number of patterns passing min_freq.}
#'   \item{params}{List of sub, min_freq, test, iter, adjust.}
#' }
#'
#' @examples
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 60, TRUE),
#'   V2 = sample(LETTERS[1:4], 60, TRUE),
#'   V3 = sample(LETTERS[1:4], 60, TRUE),
#'   V4 = sample(LETTERS[1:4], 60, TRUE)
#' )
#' grp <- rep(c("X", "Y"), 30)
#' net <- build_network(seqs, method = "relative")
#' res <- sequence_compare(net, group = grp, sub = 2:3, test = "chisq")
#'
#' @importFrom data.table data.table rbindlist setDT setorderv
#' @importFrom stats p.adjust sd chisq.test
#' @export
sequence_compare <- function(x, group = NULL, sub = 3:5,
                             min_freq = 5L,
                             test = c("permutation", "chisq", "none"),
                             iter = 1000L, adjust = "fdr") {

  test <- match.arg(test)
  min_freq <- as.integer(min_freq)
  iter <- as.integer(iter)
  sub <- as.integer(sub)
  stopifnot(
    length(sub) >= 1L, all(sub >= 1L),
    min_freq >= 1L,
    iter >= 1L
  )
  adjust <- match.arg(adjust, stats::p.adjust.methods)

  # ---- Extract sequences and group labels ----
  parsed <- .cs_parse_input(x, group)
  seq_mat <- parsed$seq_mat
  group_vec <- parsed$group
  groups <- sort(unique(group_vec))
  n_groups <- length(groups)
  stopifnot("Need at least 2 groups" = n_groups >= 2L)

  n <- nrow(seq_mat)
  p <- ncol(seq_mat)
  sub <- sub[sub <= p]
  stopifnot("No valid pattern lengths (sub > sequence length)" = length(sub) >= 1L)

  pat_list <- .cs_extract_patterns(seq_mat, sub)

  group_idx <- lapply(groups, function(g) which(group_vec == g))
  names(group_idx) <- groups

  row_list <- list()

  for (si in seq_along(sub)) {
    k <- sub[si]
    pat_mat <- pat_list[[si]]
    unique_pats <- sort(unique(pat_mat[nzchar(pat_mat)]))
    if (!length(unique_pats)) next

    counts_all <- vapply(groups, function(g) {
      vals <- pat_mat[, group_idx[[g]], drop = FALSE]
      idx <- match(vals, unique_pats, nomatch = 0L)
      tabulate(idx, nbins = length(unique_pats))
    }, integer(length(unique_pats)))
    if (!is.matrix(counts_all)) counts_all <- matrix(counts_all, ncol = n_groups)
    rownames(counts_all) <- unique_pats
    colnames(counts_all) <- groups

    total_by_group <- colSums(counts_all)
    N_grand <- sum(total_by_group)
    if (N_grand == 0) next

    # min_freq filter: every group must meet threshold
    keep <- apply(counts_all, 1, min) >= min_freq
    if (!any(keep)) next
    counts <- counts_all[keep, , drop = FALSE]
    kept_pats <- rownames(counts)
    n_pat <- nrow(counts)

    props <- sweep(counts, 2, total_by_group, `/`)
    props[!is.finite(props)] <- NA_real_

    # Standardized residuals: (o - e) / sqrt(e * (1 - r/N) * (1 - c/N))
    # 2xG table is (pattern, total - pattern). Group totals come from the
    # unfiltered count matrix so "rest" includes rare patterns.
    resids <- matrix(0, nrow = n_pat, ncol = n_groups,
                     dimnames = list(kept_pats, groups))
    col_props <- total_by_group / N_grand
    row_sums <- rowSums(counts)
    row_props <- row_sums / N_grand
    for (i in seq_len(n_pat)) {
      row1 <- counts[i, ]
      exp1 <- row_sums[i] * col_props
      denom <- sqrt(exp1 * (1 - row_props[i]) * (1 - col_props))
      resids[i, ] <- (row1 - exp1) / denom
    }
    resids[!is.finite(resids)] <- 0

    freq_df  <- as.data.frame(counts, stringsAsFactors = FALSE)
    names(freq_df)  <- paste0("freq_",  groups)
    prop_df  <- as.data.frame(props,  stringsAsFactors = FALSE)
    names(prop_df)  <- paste0("prop_",  groups)
    resid_df <- as.data.frame(resids, stringsAsFactors = FALSE)
    names(resid_df) <- paste0("resid_", groups)

    rows <- data.table(pattern = kept_pats, length = k)
    rows <- cbind(rows, freq_df, prop_df, resid_df)

    # ---- Inference ----
    if (test == "chisq") {
      stat_vec <- numeric(n_pat)
      p_vec <- numeric(n_pat)
      for (i in seq_len(n_pat)) {
        row1 <- counts[i, ]
        row2 <- total_by_group - row1
        tab <- rbind(row1, row2)
        ch <- suppressWarnings(stats::chisq.test(tab))
        stat_vec[i] <- unname(ch$statistic)
        p_vec[i] <- ch$p.value
      }
      rows$statistic <- stat_vec
      rows$p_value <- p_vec
    } else if (test == "permutation") {
      stat_obs <- .cs_pattern_statistic(counts)
      stat_perm <- matrix(0, nrow = n_pat, ncol = iter)
      stat_ge <- integer(n_pat)

      for (it in seq_len(iter)) {
        perm_group <- group_vec[sample(n)]
        perm_idx <- lapply(groups, function(g) which(perm_group == g))
        perm_counts_all <- vapply(seq_along(groups), function(gi) {
          vals <- pat_mat[, perm_idx[[gi]], drop = FALSE]
          idx <- match(vals, unique_pats, nomatch = 0L)
          tabulate(idx, nbins = length(unique_pats))
        }, integer(length(unique_pats)))
        if (!is.matrix(perm_counts_all)) {
          perm_counts_all <- matrix(perm_counts_all, ncol = n_groups)
        }
        rownames(perm_counts_all) <- unique_pats
        perm_counts <- perm_counts_all[kept_pats, , drop = FALSE]
        stat_perm[, it] <- .cs_pattern_statistic(perm_counts)
        stat_ge <- stat_ge + (stat_perm[, it] >= stat_obs)
      }

      perm_mean <- rowMeans(stat_perm)
      perm_sd <- apply(stat_perm, 1, sd)
      perm_sd[perm_sd == 0] <- 1
      rows$effect_size <- (stat_obs - perm_mean) / perm_sd
      rows$p_value <- (stat_ge + 1) / (iter + 1)
    }
    # test == "none" adds nothing

    row_list[[length(row_list) + 1L]] <- rows
  }

  if (!length(row_list)) {
    stop(sprintf("No patterns with min_freq >= %d in all groups.", min_freq),
         call. = FALSE)
  }

  results <- rbindlist(row_list)

  # Adjust p-values ONCE across the full pattern set
  if ("p_value" %in% names(results) && adjust != "none") {
    results$p_value <- stats::p.adjust(results$p_value, method = adjust)
  }

  # Sort: by p_value if tested, else by max |residual|
  if ("p_value" %in% names(results)) {
    if ("effect_size" %in% names(results)) {
      results <- results[order(results$p_value, -abs(results$effect_size)), ]
    } else {
      results <- results[order(results$p_value, -results$statistic), ]
    }
  } else {
    resid_cols <- paste0("resid_", groups)
    max_r <- do.call(pmax, lapply(resid_cols, function(cn) abs(results[[cn]])))
    results <- results[order(-max_r), ]
  }
  rownames(results) <- NULL

  structure(list(
    patterns   = as.data.frame(results, stringsAsFactors = FALSE),
    groups     = groups,
    n_patterns = nrow(results),
    params     = list(sub = sub, min_freq = min_freq,
                      test = test, iter = iter, adjust = adjust)
  ), class = "net_sequence_comparison")
}


# ---- Input parsing ----

#' @noRd
.cs_parse_input <- function(x, group) {
  # netobject_group: each element has $data
  if (inherits(x, "netobject_group")) {
    group_names <- names(x)
    mats <- lapply(seq_along(x), function(i) {
      d <- x[[i]]$data
      if (is.null(d)) stop("netobject_group element has no $data.", call. = FALSE)
      m <- as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
      # Decode integer-encoded data
      lbl <- rownames(x[[i]]$weights)
      if (!is.null(lbl) && length(lbl) > 0L &&
          (is.integer(m[, 1]) || is.numeric(m[, 1]))) {
        storage.mode(m) <- "character"
        m[] <- lbl[as.integer(m)]
      }
      m
    })
    # Align column count
    max_cols <- max(vapply(mats, ncol, integer(1)))
    mats <- lapply(mats, function(m) {
      if (ncol(m) < max_cols) {
        pad <- matrix(NA_character_, nrow(m), max_cols - ncol(m))
        cbind(m, pad)
      } else m
    })
    seq_mat <- do.call(rbind, mats)
    group_vec <- rep(group_names, vapply(mats, nrow, integer(1)))
    return(list(seq_mat = seq_mat, group = group_vec))
  }

  # netobject or tna: extract $data
  if (inherits(x, "tna") || inherits(x, "netobject") ||
      inherits(x, "cograph_network")) {
    if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
      x <- .as_netobject(x)
    }
    if (inherits(x, "tna")) {
      d <- x$data
      lbl <- x$labels
    } else {
      d <- x$data
      lbl <- rownames(x$weights)
    }
    if (is.null(d)) stop("Object has no $data.", call. = FALSE)
    seq_mat <- as.matrix(as.data.frame(d, stringsAsFactors = FALSE))
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(seq_mat[, 1]) || is.numeric(seq_mat[, 1]))) {
      storage.mode(seq_mat) <- "character"
      seq_mat[] <- lbl[as.integer(seq_mat)]
    }
  } else if (is.data.frame(x) || is.matrix(x)) {
    seq_mat <- as.matrix(as.data.frame(x, stringsAsFactors = FALSE))
  } else {
    stop("x must be a netobject_group, netobject, tna, data.frame, or matrix.",
         call. = FALSE)
  }

  # Resolve group
  if (is.null(group)) {
    stop("'group' is required when x is not a netobject_group.", call. = FALSE)
  }
  if (length(group) == 1L && is.character(group) && group %in% colnames(seq_mat)) {
    group_vec <- seq_mat[, group]
    seq_mat <- seq_mat[, colnames(seq_mat) != group, drop = FALSE]
  } else {
    stopifnot("'group' length must match nrow(x)" = length(group) == nrow(seq_mat))
    group_vec <- as.character(group)
  }

  list(seq_mat = seq_mat, group = group_vec)
}


# ---- Pattern extraction (vectorized) ----

#' @noRd
.cs_extract_patterns <- function(m, sub) {
  n <- nrow(m)
  p <- ncol(m)
  mis <- is.na(m)

  lapply(sub, function(k) {
    n_pos <- p - k + 1L
    if (n_pos < 1L) {
      return(matrix("", nrow = 0, ncol = n))
    }
    # Build pattern strings for each position × sequence
    out <- matrix("", nrow = n_pos, ncol = n)
    for (i in seq_len(n_pos)) {
      idx <- i:(i + k - 1L)
      has_na <- rowSums(mis[, idx, drop = FALSE]) > 0L
      cols <- m[, idx, drop = FALSE]
      pat <- do.call(paste, c(as.data.frame(cols, stringsAsFactors = FALSE),
                              list(sep = "->")))
      pat[has_na] <- ""
      out[i, ] <- pat
    }
    out
  })
}


# ---- Test statistic ----

#' @noRd
.cs_pattern_statistic <- function(x) {
  n <- sum(x)
  nr <- nrow(x)
  nc <- ncol(x)
  rs <- rowSums(x)
  cs <- colSums(x)
  expected <- outer(rs, cs) / n
  sqrt(rowSums((x - expected)^2))
}


# ---- S3 Methods ----

#' Print Method for net_sequence_comparison
#'
#' @param x A \code{net_sequence_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @inherit sequence_compare examples
#' @export
print.net_sequence_comparison <- function(x, ...) {
  cat(sprintf("Sequence Comparison  [%d patterns | %d groups: %s]\n",
              x$n_patterns, length(x$groups), paste(x$groups, collapse = ", ")))
  cat(sprintf("  Lengths: %s  |  min_freq: %d",
              paste(x$params$sub, collapse = ", "), x$params$min_freq))
  test <- x$params$test
  if (identical(test, "permutation")) {
    cat(sprintf("  |  permutation: %d iter (%s)",
                x$params$iter, x$params$adjust))
  } else if (identical(test, "chisq")) {
    cat(sprintf("  |  chi-square  (%s)", x$params$adjust))
  }
  cat("\n")

  if (x$n_patterns > 0) {
    top <- utils::head(x$patterns, 10)
    cat("\n")
    print(top, row.names = FALSE)
    if (x$n_patterns > 10) {
      cat(sprintf("  ... and %d more patterns\n", x$n_patterns - 10))
    }
  }
  invisible(x)
}


#' Summary Method for net_sequence_comparison
#'
#' @param object A \code{net_sequence_comparison} object.
#' @param ... Additional arguments (ignored).
#' @return The patterns data.frame (tidy: one row per k-gram pattern, per
#'   group, with frequency and proportion columns; includes p-values when a
#'   permutation test was run).
#' @inherit sequence_compare examples
#' @export
summary.net_sequence_comparison <- function(object, ...) {
  object$patterns
}


#' Plot Method for net_sequence_comparison
#'
#' @description
#' Visualizes pattern-level standardized residuals across groups. Two styles
#' are available:
#' \describe{
#'   \item{\code{"pyramid"}}{Back-to-back bars of pattern proportions, shaded
#'     by each side's standardized residual. Requires exactly 2 groups.}
#'   \item{\code{"heatmap"}}{One tile per (pattern, group) cell, colored by
#'     standardized residual. Works for any number of groups.}
#' }
#' Residuals are read directly from the \code{resid_<group>} columns in
#' \code{$patterns}, which are always populated regardless of the inference
#' method chosen in \code{sequence_compare}.
#'
#' @param x A \code{net_sequence_comparison} object.
#' @param top_n Integer. Show top N patterns. Default: 10.
#' @param style Character. \code{"pyramid"} (default) or \code{"heatmap"}.
#' @param sort Character. \code{"statistic"} (default) ranks patterns by test
#'   statistic or residual magnitude. \code{"frequency"} ranks by total
#'   occurrence count across all groups.
#' @param alpha Numeric. Significance threshold for p-value display in the
#'   pyramid: patterns with \code{p_value < alpha} are starred and drawn in
#'   bold dark text, the rest stay plain grey. Default: 0.05.
#' @param show_residuals Logical. If \code{TRUE}, print the standardized
#'   residual value inside each pyramid bar. Default: \code{FALSE}. Ignored
#'   for the heatmap (which always shows residuals).
#' @param ... Additional arguments (ignored).
#' @return A \code{ggplot} object, invisibly.
#' @inherit sequence_compare examples
#' @import ggplot2
#' @export
plot.net_sequence_comparison <- function(x, top_n = 10L,
                                          style = c("pyramid", "heatmap"),
                                          sort = c("statistic", "frequency"),
                                          alpha = 0.05,
                                          show_residuals = FALSE, ...) {
  style <- match.arg(style)
  sort  <- match.arg(sort)
  pat <- x$patterns
  groups <- x$groups
  if (nrow(pat) == 0) {
    message("No patterns to plot.")
    return(invisible(NULL))
  }

  resid_cols <- paste0("resid_", groups)
  freq_cols  <- paste0("freq_",  groups)

  if (sort == "frequency") {
    total_freq <- Reduce("+", lapply(freq_cols, function(cn) pat[[cn]]))
    pat <- pat[order(-total_freq), ]
  } else {
    if ("effect_size" %in% names(pat)) {
      pat <- pat[order(-abs(pat$effect_size)), ]
    } else if ("statistic" %in% names(pat)) {
      pat <- pat[order(-pat$statistic), ]
    } else {
      max_r <- do.call(pmax, lapply(resid_cols, function(cn) abs(pat[[cn]])))
      pat <- pat[order(-max_r), ]
    }
  }
  pat <- utils::head(pat, top_n)

  if (style == "heatmap") {
    return(.sc_plot_heatmap(pat, groups, resid_cols))
  }
  .sc_plot_pyramid(pat, groups, alpha, show_residuals)
}

#' @noRd
.sc_plot_pyramid <- function(pat, groups, alpha, show_residuals = FALSE) {
  if (length(groups) != 2L) {
    stop("style = 'pyramid' requires exactly 2 groups; got ",
         length(groups), ". Use style = 'heatmap' instead.",
         call. = FALSE)
  }
  g1 <- groups[1]; g2 <- groups[2]
  prop1 <- paste0("prop_", g1); prop2 <- paste0("prop_", g2)
  freq1 <- paste0("freq_", g1); freq2 <- paste0("freq_", g2)
  rcol1 <- paste0("resid_", g1); rcol2 <- paste0("resid_", g2)

  # Preserve the rank order passed in (top-ranked first) so the pyramid
  # agrees with the heatmap. Reverse factor levels so row 1 appears at
  # the top of the y-axis rather than the bottom.
  resid1 <- pat[[rcol1]]
  resid2 <- pat[[rcol2]]

  pyramid <- data.frame(
    pattern = factor(pat$pattern, levels = rev(pat$pattern)),
    left    = -pat[[prop1]],
    right   =  pat[[prop2]],
    freq_left   = pat[[freq1]],
    freq_right  = pat[[freq2]],
    resid_left  = resid1,
    resid_right = resid2,
    stringsAsFactors = FALSE
  )

  r_max <- max(abs(c(resid1, resid2)), na.rm = TRUE)
  x_max <- max(abs(c(pyramid$left, pyramid$right)), na.rm = TRUE)
  gutter <- x_max * 0.06
  has_pval <- "p_value" %in% names(pat)
  bar_hw <- 0.42
  y_num <- as.numeric(pyramid$pattern)

  rect_df <- rbind(
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = pyramid$left - gutter, xmax = -gutter,
               freq = pyramid$freq_left, resid = pyramid$resid_left,
               side = "left", stringsAsFactors = FALSE),
    data.frame(yidx = y_num, ymin = y_num - bar_hw, ymax = y_num + bar_hw,
               xmin = gutter, xmax = pyramid$right + gutter,
               freq = pyramid$freq_right, resid = pyramid$resid_right,
               side = "right", stringsAsFactors = FALSE)
  )

  if (has_pval) {
    p_num <- ifelse(pat$p_value < 0.001, "<.001",
              ifelse(pat$p_value < 0.01,
                     sprintf("%.3f", pat$p_value),
                     sprintf("%.2f", pat$p_value)))
    # `alpha` is the significance threshold for the pyramid's p-value
    # display: patterns with p < alpha are starred and drawn bold/dark,
    # the rest stay plain grey so the reader sees which patterns clear
    # the chosen level.
    p_sig <- !is.na(pat$p_value) & pat$p_value < alpha
    p_df <- data.frame(
      yidx       = y_num,
      p_label    = ifelse(p_sig, paste0(p_num, "*"), p_num),
      sig        = p_sig,
      p_face     = ifelse(p_sig, "bold.italic", "italic"),
      p_color    = ifelse(p_sig, "grey15", "grey55"),
      stringsAsFactors = FALSE
    )
  }

  # Residual labels inside each bar (optional)
  if (show_residuals) {
    rect_df$x_center <- (rect_df$xmin + rect_df$xmax) / 2
    rect_df$resid_label <- sprintf("%.2f", rect_df$resid)
    rect_df$text_color <- ifelse(abs(rect_df$resid) > r_max * 0.35,
                                  "white", "grey30")
  }

  p <- ggplot() +
    geom_rect(data = rect_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                  fill = resid)) +
    geom_text(data = rect_df,
              aes(x = ifelse(side == "left", xmin, xmax),
                  y = yidx, label = freq),
              hjust = ifelse(rect_df$side == "left", 1.15, -0.15),
              size = 3, fontface = "bold",
              color = ifelse(abs(rect_df$resid) > r_max * 0.5, "grey25", "grey45")) +
    {if (show_residuals)
      geom_text(data = rect_df,
                aes(x = x_center, y = yidx, label = resid_label),
                size = 2.8, fontface = "bold",
                color = rect_df$text_color)
    } +
    {if (has_pval)
      geom_text(data = p_df,
                aes(x = 0, y = yidx, label = p_label),
                size = 2.5, color = p_df$p_color, fontface = p_df$p_face)
    } +
    scale_fill_gradient2(low = "#CB181D", mid = "white", high = "#2171B5",
                         midpoint = 0, limits = c(-3, 3),
                         oob = scales::squish,
                         name = "Standardized\nresidual") +
    scale_y_continuous(breaks = y_num,
                       labels = as.character(pat$pattern)) +
    scale_x_continuous(
      labels = function(x) format(abs(x), scientific = FALSE)
    ) +
    labs(x = "Proportion", y = NULL,
         title = "Sequence Pattern Comparison",
         subtitle = paste0(g1, "  \u2190\u2190    \u2192\u2192  ", g2)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y          = element_text(family = "mono", size = 8),
      axis.ticks.y         = element_blank(),
      panel.grid.major.y   = element_blank(),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle        = element_text(hjust = 0.5, color = "grey40",
                                           size = 11),
      legend.position      = "right",
      plot.margin          = margin(5, 10, 5, 5)
    )

  print(p)
  invisible(p)
}

#' @noRd
.sc_plot_heatmap <- function(pat, groups, resid_cols) {
  long <- as.data.frame(data.table::rbindlist(lapply(seq_along(groups), function(gi) {
    data.frame(
      pattern = pat$pattern,
      group   = groups[gi],
      resid   = pat[[resid_cols[gi]]],
      stringsAsFactors = FALSE
    )
  })))
  long$pattern <- factor(long$pattern, levels = rev(pat$pattern))
  long$group   <- factor(long$group, levels = groups)

  r_max <- max(abs(long$resid), na.rm = TRUE)
  # Saturate at ±3 (|z|>3 is the conventional "very strong" threshold).
  # Anything beyond clips to full color via scales::squish so mid-range
  # residuals (±1 to ±2) pick up visible color instead of fading to white.
  r_lim <- 3

  p <- ggplot(long, aes(x = group, y = pattern, fill = resid)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", resid)),
              color = "white", size = 5, fontface = "bold") +
    scale_fill_gradient2(low = "#CB181D", mid = "white", high = "#2171B5",
                         midpoint = 0, limits = c(-r_lim, r_lim),
                         oob = scales::squish,
                         name = "Standardized\nresidual") +
    labs(x = "Groups", y = NULL, title = "Sequence Pattern Comparison") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y        = element_text(family = "mono", size = 10),
      axis.text.x        = element_text(angle = 0, size = 11,
                                         face = "bold"),
      axis.title.x       = element_text(size = 14),
      panel.grid         = element_blank(),
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 15),
      legend.position    = "right",
      legend.title       = element_text(size = 12),
      legend.text        = element_text(size = 11),
      legend.key.height  = unit(1.2, "cm"),
      plot.margin        = margin(10, 15, 10, 10)
    )

  print(p)
  invisible(p)
}
