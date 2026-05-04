# ==============================================================================
# plot_state_frequencies: marimekko + colored-bar state-frequency plotter
#
# Native Nestimate plotter for state (node) frequency distributions across
# the group-bearing classes (netobject, netobject_group, mcml, htna).
#
# Distinct from tna::plot_frequencies (different name, different geometry).
# Defaults to a marimekko (mosaic) layout where column widths reflect
# group totals and segment heights reflect within-group state proportions.
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

# Drop-in safe state_frequencies() that ignores NA / void markers. Returns
# a data.frame(state, count, proportion) sorted by count descending.
.freq_count_states <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0L) {
    return(data.frame(state = character(0), count = integer(0),
                      proportion = numeric(0), stringsAsFactors = FALSE))
  }
  tbl <- sort(table(values), decreasing = TRUE)
  data.frame(state = names(tbl),
             count = as.integer(tbl),
             proportion = as.numeric(tbl) / sum(tbl),
             stringsAsFactors = FALSE)
}


# Tidy a single (group, data) pair into a freq data.frame. Returns NULL when
# the group has no countable states.
.freq_one_group <- function(group_label, data) {
  if (is.null(data)) return(NULL)
  if (is.data.frame(data) || is.matrix(data)) {
    values <- as.vector(as.matrix(data))
  } else if (is.list(data)) {
    values <- unlist(data, use.names = FALSE)
  } else {
    values <- data
  }
  fr <- .freq_count_states(values)
  if (nrow(fr) == 0L) return(NULL)
  fr$group <- as.character(group_label)
  fr[, c("group", "state", "count", "proportion")]
}


# ------------------------------------------------------------------------------
# Per-class extractors. Each returns data.frame(group, state, count, proportion).
# Public generic state_distribution() at the bottom of the file forwards here.
# ------------------------------------------------------------------------------

.freq_df_netobject <- function(x) {
  if (is.null(x$data)) {
    stop("netobject has no `$data` to compute frequencies from.",
         call. = FALSE)
  }

  # If $node_groups present, group states by node_groups$group; otherwise
  # treat the whole data as a single group ("all").
  ng <- x$node_groups
  if (is.null(ng) || nrow(ng) == 0L) {
    return(.freq_one_group("all", x$data))
  }

  # Compute frequency over $data, then attach group via node->group lookup.
  fr <- .freq_count_states(as.vector(as.matrix(x$data)))
  if (nrow(fr) == 0L) return(fr[, integer(0)])
  lookup <- setNames(as.character(ng$group), as.character(ng$node))
  fr$group <- unname(lookup[fr$state])
  fr$group[is.na(fr$group)] <- "(ungrouped)"
  # Re-aggregate per (group, state). State is unique per row already because
  # .freq_count_states returns one row per unique state, but a state could
  # belong to one group only -- keep as-is. Recompute within-group proportion.
  totals <- ave(fr$count, fr$group, FUN = sum)
  fr$proportion <- fr$count / totals
  fr[, c("group", "state", "count", "proportion")]
}


.freq_df_htna <- function(x) {
  # htna inherits from netobject; node_groups is always populated by build_htna.
  .freq_df_netobject(x)
}


.freq_df_mcml <- function(x, include_macro = FALSE) {
  if (is.null(x$clusters)) {
    stop("mcml object has no `$clusters` field.", call. = FALSE)
  }
  cluster_names <- names(x$clusters)
  parts <- lapply(cluster_names, function(cn) {
    cl <- x$clusters[[cn]]
    .freq_one_group(cn, cl$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop("mcml clusters have no `$data` to compute frequencies from.",
         call. = FALSE)
  }
  out <- do.call(rbind, parts)

  # Apply the labels remap stored by `.apply_node_labels` if present, so
  # raw codes in $data (e.g. "c_evaluation") are shown as their human-
  # readable labels (e.g. "Evaluation").
  map <- attr(x, "labels_map")
  if (!is.null(map)) {
    hits <- out$state %in% names(map)
    out$state[hits] <- unname(map[out$state[hits]])
  }

  if (isTRUE(include_macro)) {
    macro <- aggregate(out$count, by = list(state = out$state), FUN = sum)
    macro$group <- "macro"
    macro$proportion <- macro$x / sum(macro$x)
    macro <- data.frame(group = macro$group, state = macro$state,
                        count = macro$x, proportion = macro$proportion,
                        stringsAsFactors = FALSE)
    out <- rbind(macro, out)
  }
  out
}


.freq_df_netobject_group <- function(x) {
  group_names <- names(x)
  if (is.null(group_names)) group_names <- paste0("Group ", seq_along(x))
  parts <- lapply(seq_along(x), function(i) {
    .freq_one_group(group_names[i], x[[i]]$data)
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0L) {
    stop("netobject_group has no `$data` in any constituent.",
         call. = FALSE)
  }
  do.call(rbind, parts)
}


# ------------------------------------------------------------------------------
# State ordering and color assignment
# ------------------------------------------------------------------------------

.order_states <- function(states, counts, sort_states) {
  if (sort_states == "none") return(unique(states))
  if (sort_states == "alpha") return(sort(unique(states)))
  agg <- aggregate(counts, by = list(state = states), FUN = sum)
  agg <- agg[order(-agg$x), , drop = FALSE]
  agg$state
}


# ------------------------------------------------------------------------------
# Mosaic primitive (exported)
# ------------------------------------------------------------------------------

#' Draw a Marimekko / Mosaic Plot from a Tidy Data Frame
#'
#' Low-level rectangle-coordinate builder for marimekko (mosaic) plots.
#' Column widths are proportional to the per-column total of \code{weight};
#' within each column, segments stack to height 1 with sub-heights
#' proportional to each row's share of that column's total.
#'
#' Used internally by \code{\link{plot_state_frequencies}}; exposed so that
#' other plot methods (e.g. permutation-residual visualisations) can reuse
#' the same geometry by supplying a different fill column.
#'
#' @param data A data.frame in long form. Must contain the columns named
#'   in \code{x}, \code{y}, and \code{weight}.
#' @param x Column name for the X (column) variable.
#' @param y Column name for the Y (segment) variable.
#' @param weight Column name for the cell weight (e.g. count).
#' @param fill Either \code{"y"} (color by Y category, e.g. state -- default)
#'   or the name of another column to map to fill (e.g. a residual column
#'   for diverging color).
#' @param colors Optional character vector of fill colors. When
#'   \code{fill = "y"}, length must be at least the number of distinct y
#'   levels. Defaults to recycled Okabe-Ito.
#' @param show_labels If \code{TRUE}, draw within-segment percentage labels.
#' @param label_size Numeric size for segment labels.
#' @param x_label,y_label Optional axis labels.
#'
#' @return A \code{ggplot} object.
#' @import ggplot2
#' @export
#' @examples
#' df <- data.frame(
#'   group = rep(c("A", "B", "C"), each = 3),
#'   state = rep(c("s1", "s2", "s3"), 3),
#'   count = c(10, 5, 2,  7, 8, 3,  4, 6, 12)
#' )
#' plot_mosaic(df, x = "group", y = "state", weight = "count")
plot_mosaic <- function(data,
                        x,
                        y,
                        weight,
                        fill        = "y",
                        colors      = NULL,
                        show_labels = TRUE,
                        label_size  = 3.5,
                        x_label     = NULL,
                        y_label     = NULL) {
  stopifnot(is.data.frame(data),
            is.character(x), length(x) == 1L,
            is.character(y), length(y) == 1L,
            is.character(weight), length(weight) == 1L)
  needed <- c(x, y, weight)
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop("plot_mosaic: missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Build cross-tab. Column = x, row = y. Order respects existing factor
  # levels if present; otherwise alphabetical.
  xv <- if (is.factor(data[[x]])) data[[x]] else factor(data[[x]])
  yv <- if (is.factor(data[[y]])) data[[y]] else factor(data[[y]])
  tab <- tapply(data[[weight]], list(xv, yv), FUN = sum, default = 0)
  tab[is.na(tab)] <- 0
  x_levels <- rownames(tab)
  y_levels <- colnames(tab)

  # Column widths from row sums, segment heights from within-row proportions.
  col_totals <- rowSums(tab)
  total      <- sum(col_totals)
  if (total <= 0) {
    stop("plot_mosaic: total weight is zero -- nothing to draw.", call. = FALSE)
  }
  widths_cum <- c(0, cumsum(col_totals)) / total

  # Build rect coordinates: one row per (x_level, y_level) cell.
  rects <- do.call(rbind, lapply(seq_along(x_levels), function(i) {
    row    <- tab[i, , drop = TRUE]
    rs     <- sum(row)
    if (rs <= 0) return(NULL)
    heights_cum <- c(0, cumsum(row / rs))
    data.frame(
      x_level = rep(x_levels[i], length(y_levels)),
      y_level = y_levels,
      xmin    = widths_cum[i],
      xmax    = widths_cum[i + 1L],
      ymin    = heights_cum[seq_len(length(y_levels))],
      ymax    = heights_cum[seq_len(length(y_levels)) + 1L],
      count   = as.numeric(row),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  rects <- rects[rects$ymax > rects$ymin, , drop = FALSE]

  # Attach extra fill column if requested
  if (!identical(fill, "y") && fill %in% names(data)) {
    extra <- aggregate(data[[fill]],
                       by = list(x_level = as.character(xv),
                                 y_level = as.character(yv)),
                       FUN = function(z) z[1L])
    names(extra)[3L] <- fill
    rects <- merge(rects, extra,
                   by = c("x_level", "y_level"), all.x = TRUE)
  }

  # Map fill aesthetic
  fill_var <- if (identical(fill, "y")) "y_level" else fill
  rects[[fill_var]] <- if (identical(fill, "y")) {
    factor(rects$y_level, levels = y_levels)
  } else {
    rects[[fill_var]]
  }

  # Build palette
  pal <- if (identical(fill, "y")) {
    .state_palette(colors, length(y_levels))
  } else {
    NULL
  }

  # Column midpoint x positions for axis ticks
  mid_x <- (widths_cum[-1L] + widths_cum[-length(widths_cum)]) / 2

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data[[fill_var]])) +
    ggplot2::geom_rect(color = "white", linewidth = 0.3) +
    ggplot2::scale_x_continuous(breaks = mid_x, labels = x_levels,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                labels = function(v) sprintf("%d%%", round(v * 100))) +
    ggplot2::labs(x = x_label %||% x, y = y_label %||% "Proportion",
                  fill = if (identical(fill, "y")) y else fill) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (!is.null(pal)) {
    p <- p + ggplot2::scale_fill_manual(values = pal, drop = FALSE)
  }

  if (isTRUE(show_labels)) {
    rects$pct <- rects$ymax - rects$ymin
    rects$x_mid <- (rects$xmin + rects$xmax) / 2
    rects$y_mid <- (rects$ymin + rects$ymax) / 2
    rects$lab   <- sprintf("%.1f%%", 100 * rects$pct)
    rects$lab[rects$pct < 0.04] <- ""
    p <- p + .geom_fit_label(rects, label_size)
  }
  p
}


# Tiny null-coalesce helper (avoid rlang dep)
`%||%` <- function(a, b) if (is.null(a)) b else a


# Fit-aware tile labels. Uses ggfittext::geom_fit_text() when available so
# each tile's text auto-shrinks (and reflows onto multiple lines) to fit
# its rectangle, with `min.size` falling back to dropping the label when
# nothing legible fits. When ggfittext is not installed, returns a regular
# ggplot2::geom_text() at midpoint with the (unfitted) label_size; the
# upstream caller is expected to have already nulled out labels for tiles
# below `tile_area < 0.05`.
.geom_fit_label <- function(rects, label_size, color = "grey15") {
  if (is.null(rects$angle)) {
    rects$angle <- ifelse((rects$ymax - rects$ymin) >
                            (rects$xmax - rects$xmin), 90, 0)
  }
  if (requireNamespace("ggfittext", quietly = TRUE)) {
    return(ggfittext::geom_fit_text(
      data = rects,
      ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                   ymin = .data$ymin, ymax = .data$ymax,
                   label = .data$lab,
                   angle = .data$angle),
      inherit.aes = FALSE,
      reflow      = TRUE,
      min.size    = 1,
      grow        = FALSE,
      padding.x   = grid::unit(0.6, "mm"),
      padding.y   = grid::unit(0.6, "mm"),
      colour      = color,
      size        = label_size * 3.2
    ))
  }
  ggplot2::geom_text(
    data = rects,
    ggplot2::aes(x = .data$x_mid, y = .data$y_mid, label = .data$lab,
                 angle = .data$angle),
    inherit.aes = FALSE, size = label_size, color = color
  )
}


# ------------------------------------------------------------------------------
# mosaic_plot: tna-equivalent chi-square mosaic for netobjects
# ------------------------------------------------------------------------------

#' Mosaic Plot of a Network's Transition or Co-occurrence Counts
#'
#' Draws a Hartigan-Friendly mosaic (marimekko geometry, chi-square
#' standardized-residual fill) for an integer-weighted network. Equivalent in
#' algorithm and appearance to \code{tna::plot_mosaic()}; named differently to
#' avoid an export clash when both packages are attached.
#'
#' Column widths are proportional to row marginals of the weight matrix
#' (incoming totals when the matrix is transposed, as for transitions). Within
#' each column, segment heights are proportional to that row's conditional
#' distribution. Cell fill is the standardized residual from
#' \code{stats::chisq.test()}, with a diverging palette clipped to \eqn{\pm 4}.
#' Mosaics only make sense for integer-valued matrices, so this function rejects
#' relative / glasso / correlation networks.
#'
#' @param x A \code{netobject} (single network) or \code{netobject_group} (one
#'   mosaic per group, arranged with \code{gridExtra::grid.arrange}).
#' @param xlab,ylab Axis labels. Defaults match tna's wording for transition
#'   matrices: \code{"Incoming edges"} on x, \code{"Outgoing edges"} on y.
#' @param range Numeric of length 2 giving the lower and upper colour-scale
#'   limits for the standardized residual. \code{NULL} (default) auto-fits
#'   the limits to the symmetric range \code{c(-M, M)} where
#'   \code{M = max(|stdres|)}, so no signal is squished. Pass an explicit
#'   range (e.g. \code{c(-4, 4)} for tna-style display, \code{c(-6, 6)} for
#'   moderate clipping) to clamp the colour scale.
#' @param top_angle,left_angle Rotation in degrees for the top (x) and left
#'   (y) tick labels. \code{NULL} (default) uses the auto rule
#'   \code{90 if n_levels > 3 else 0} on each axis. Pass any numeric to
#'   override (e.g. \code{top_angle = 45, left_angle = 0}).
#' @param residuals One of \code{"permutation"} (default) or
#'   \code{"asymptotic"}. \code{"permutation"} computes empirical-null
#'   z-scores by shuffling one variable's labels against the other for
#'   \code{n_perm} draws and reporting
#'   \code{(O - mean_perm) / sd_perm} per cell. Robust on sparse tables.
#'   \code{"asymptotic"} returns \code{stats::chisq.test()$stdres} (the
#'   closed-form \code{(O - E) / sqrt(E*(1 - p_row)*(1 - p_col))} that vcd
#'   and tna use).
#' @param n_perm Number of permutations when \code{residuals = "permutation"}.
#'   Default 500; use \code{>= 1000} for stable tail estimates.
#' @param seed Optional integer seed for the permutation RNG. Use for
#'   reproducible plots; ignored when \code{residuals = "asymptotic"}.
#' @param ncol For \code{netobject_group}: number of columns in the small-
#'   multiples layout. Default 2.
#' @param ... Ignored.
#'
#' @return A \code{ggplot} object (or a \code{gtable} from
#'   \code{gridExtra::arrangeGrob} for \code{netobject_group} when
#'   \pkg{gridExtra} is available).
#' @seealso \code{\link{plot_mosaic}} for the lower-level data.frame primitive.
#' @export
#' @examples
#' \dontrun{
#'   net <- build_network(group_regulation, method = "frequency")
#'   mosaic_plot(net)
#' }
mosaic_plot <- function(x, ...) UseMethod("mosaic_plot")

#' @export
#' @rdname mosaic_plot
mosaic_plot.default <- function(x, ...) {
  stop("mosaic_plot: no method for class ",
       paste(class(x), collapse = "/"),
       ". Use plot_mosaic() for a tidy data.frame.", call. = FALSE)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.netobject <- function(x,
                                  xlab = "Incoming edges",
                                  ylab = "Outgoing edges",
                                  range = NULL,
                                  top_angle = 90,
                                  left_angle = 0,
                                  residuals = c("permutation", "asymptotic"),
                                  n_perm = 500L,
                                  seed = NULL,
                                  ...) {
  residuals <- match.arg(residuals)
  w <- x$weights
  if (!is.matrix(w) || !is.numeric(w) || !all(is.finite(w))) {
    stop("mosaic_plot: $weights must be a finite numeric matrix.",
         call. = FALSE)
  }
  if (any(w < 0) || any(abs(w - round(w)) > 1e-8)) {
    stop("mosaic_plot is only defined for integer-valued weight matrices ",
         "(method = 'frequency' or 'co_occurrence'). Got method = ",
         x$method %||% "<unknown>", ".", call. = FALSE)
  }
  if (sum(w) <= 0) {
    stop("mosaic_plot: total weight is zero -- nothing to draw.",
         call. = FALSE)
  }
  .mosaic_plot_tab(as.table(t(w)), xlab = xlab, ylab = ylab, range = range,
                   top_angle = top_angle, left_angle = left_angle,
                   residuals = residuals, n_perm = n_perm, seed = seed)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.table <- function(x,
                              xlab = "Row",
                              ylab = "Column",
                              range = NULL,
                              top_angle = NULL,
                              left_angle = NULL,
                              residuals = c("permutation", "asymptotic"),
                              n_perm = 500L,
                              seed = NULL,
                              ...) {
  residuals <- match.arg(residuals)
  .mosaic_plot_tab(x, xlab = xlab, ylab = ylab, range = range,
                   top_angle = top_angle, left_angle = left_angle,
                   residuals = residuals, n_perm = n_perm, seed = seed)
}

#' @export
#' @rdname mosaic_plot
mosaic_plot.matrix <- function(x, ...) mosaic_plot.table(as.table(x), ...)

#' @export
#' @rdname mosaic_plot
mosaic_plot.netobject_group <- function(x,
                                        xlab = "Incoming edges",
                                        ylab = "Outgoing edges",
                                        range = NULL,
                                        top_angle = 90,
                                        left_angle = 0,
                                        residuals = c("permutation",
                                                      "asymptotic"),
                                        n_perm = 500L,
                                        seed = NULL,
                                        ncol = 2L,
                                        ...) {
  residuals <- match.arg(residuals)
  group_names <- names(x) %||% paste0("Group ", seq_along(x))
  plots <- lapply(seq_along(x), function(i) {
    p <- mosaic_plot.netobject(x[[i]], xlab = xlab, ylab = ylab,
                               range = range, top_angle = top_angle,
                               left_angle = left_angle,
                               residuals = residuals, n_perm = n_perm,
                               seed = seed, ...)
    p + ggplot2::ggtitle(group_names[i])
  })
  if (length(plots) == 1L) return(plots[[1L]])
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::arrangeGrob(grobs = plots, ncol = ncol))
  }
  warning("Install 'gridExtra' to combine grouped mosaics; returning a list.",
          call. = FALSE)
  plots
}

# Choose symmetric colour-scale limits: explicit `range` if supplied, else
# auto-fit to the actual stdres range so no signal is squished. Floors at
# +/-1 to keep the legend readable on near-independent tables.
.mosaic_residual_limits <- function(stdres, range) {
  if (!is.null(range)) {
    stopifnot(is.numeric(range), length(range) == 2L, range[1L] < range[2L])
    return(range)
  }
  m <- max(abs(stdres), na.rm = TRUE)
  m <- max(m, 1)
  c(-m, m)
}

# Pick at most 5 breaks symmetric around 0, rounded to a tidy step (1, 2, 4,
# 5, or 10) so the colour bar is scannable at any data scale.
.mosaic_residual_breaks <- function(stdres, range) {
  lim <- .mosaic_residual_limits(stdres, range)
  m <- lim[2L]
  step <- if      (m <= 4)  1
          else if (m <= 8)  2
          else if (m <= 16) 4
          else if (m <= 25) 5
          else              10
  s <- seq(0, m, by = step)
  unique(c(-rev(s), s))
}

# Build a y-axis scale for the mosaic. All labels shown at the cell
# midpoint of the leftmost column.
.mosaic_y_scale <- function(d, col_labels, ...) {
  breaks_all <- d$ycent[d$xmin == 0]
  o <- order(breaks_all)
  ggplot2::scale_y_continuous(
    breaks = breaks_all[o],
    labels = col_labels[o],
    expand = c(0.01, 0)
  )
}

# Internal: vectorized port of tna::plot_mosaic_(). Builds a marimekko data
# frame from a contingency table and renders it with chi-square stdres fill.
# `tab` rows are the x-axis categories; columns are the within-column stack.
# Permutation-based standardized residuals. Shuffles one variable's labels
# against the other (preserves both marginals under the independence null),
# tabulates n_perm times, and returns a per-cell empirical z-score
# (O - mean_perm) / sd_perm. Vectorized: each iteration is a single
# tabulate() over a column-major linear index, stacked into a (n*m) x B
# matrix; row means and row sds in closed form.
.mosaic_perm_stdres <- function(tab, n_perm = 500L, seed = NULL) {
  stopifnot(n_perm >= 2L)
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(tab); m <- ncol(tab)
  cells <- expand.grid(row_i = seq_len(n), col_j = seq_len(m))
  freqs <- as.vector(tab)
  long_row <- rep(cells$row_i, times = freqs)
  long_col <- rep(cells$col_j, times = freqs)

  perm_counts <- vapply(seq_len(n_perm), function(b) {
    pc <- sample.int(length(long_col))
    tabulate(long_row + (long_col[pc] - 1L) * n, nbins = n * m)
  }, numeric(n * m))

  mean_perm <- rowMeans(perm_counts)
  sd_perm   <- sqrt(rowMeans((perm_counts - mean_perm)^2) *
                      n_perm / (n_perm - 1L))
  obs <- as.vector(tab)
  z <- (obs - mean_perm) / sd_perm
  z[!is.finite(z)] <- 0
  matrix(z, n, m, dimnames = dimnames(tab))
}

.mosaic_plot_tab <- function(tab, xlab, ylab, range = NULL,
                             top_angle = NULL, left_angle = NULL,
                             residuals = "permutation", n_perm = 500L,
                             seed = NULL) {
  n <- nrow(tab)
  m <- ncol(tab)
  if (n < 1L || m < 1L) {
    stop("mosaic_plot: contingency table must have >= 1 row and column.",
         call. = FALSE)
  }
  rs <- rowSums(tab)
  total <- sum(rs)
  if (total <= 0) {
    stop("mosaic_plot: total weight is zero -- nothing to draw.",
         call. = FALSE)
  }
  widths <- c(0, cumsum(rs)) / total
  # heights[, i] = c(0, cumsum(tab[i, ] / rs[i])); zero-row safe.
  heights <- vapply(seq_len(n), function(i) {
    if (rs[i] <= 0) return(c(0, seq_len(m) / m))
    c(0, cumsum(as.numeric(tab[i, ]) / rs[i]))
  }, numeric(m + 1L))

  i_idx <- rep(seq_len(n), each = m)
  j_idx <- rep(seq_len(m), times = n)
  row_offset <- (i_idx - 1L) * n * 0.0025
  col_offset <- (j_idx - 1L) * m * 0.0025

  row_labels <- rownames(tab) %||% as.character(seq_len(n))
  col_labels <- colnames(tab) %||% as.character(seq_len(m))

  stdres <- if (identical(residuals, "permutation")) {
    .mosaic_perm_stdres(tab, n_perm = n_perm, seed = seed)
  } else {
    suppressWarnings(stats::chisq.test(tab))$stdres
  }
  if (is.null(stdres) || !all(is.finite(stdres))) {
    stdres <- matrix(0, n, m)
  }

  d <- data.frame(
    xmin   = widths[i_idx] + row_offset,
    xmax   = widths[i_idx + 1L] + row_offset,
    ymin   = heights[cbind(j_idx, i_idx)] + col_offset,
    ymax   = heights[cbind(j_idx + 1L, i_idx)] + col_offset,
    freq   = as.numeric(tab[cbind(i_idx, j_idx)]),
    row    = row_labels[i_idx],
    col    = col_labels[j_idx],
    stdres = as.numeric(stdres[cbind(i_idx, j_idx)]),
    stringsAsFactors = FALSE
  )
  d$xcent <- (d$xmin + d$xmax) / 2
  d$ycent <- (d$ymin + d$ymax) / 2

  ggplot2::ggplot(d,
    ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                 ymin = .data$ymin, ymax = .data$ymax,
                 fill = .data$stdres)) +
    ggplot2::geom_rect(color = "black", linewidth = 0.4,
                       show.legend = TRUE) +
    ggplot2::scale_fill_gradient2(
      name   = "Standardized\nresidual",
      oob    = scales::oob_squish,
      low    = "#D33F6A",
      high   = "#4A6FE3",
      limits = .mosaic_residual_limits(d$stdres, range),
      breaks = .mosaic_residual_breaks(d$stdres, range)
    ) +
    ggplot2::scale_x_continuous(
      breaks   = unique(d$xcent),
      labels   = row_labels,
      position = "top",
      expand   = c(0.01, 0)
    ) +
    .mosaic_y_scale(d, col_labels,
                    left_angle = left_angle %||% (if (m > 3) 90 else 0)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      panel.grid    = ggplot2::element_blank(),
      axis.ticks    = ggplot2::element_blank(),
      axis.line     = ggplot2::element_blank(),
      axis.text.x   = ggplot2::element_text(
        angle = top_angle %||% (if (n > 3) 90 else 0),
        hjust = if ((top_angle %||% (if (n > 3) 90 else 0)) == 0) 0.5 else 0,
        vjust = if ((top_angle %||% (if (n > 3) 90 else 0)) == 0) 0   else 0.5
      ),
      axis.text.y   = ggplot2::element_text(
        angle = left_angle %||% (if (m > 3) 90 else 0),
        hjust = if ((left_angle %||% (if (m > 3) 90 else 0)) == 0) 1   else 0.5,
        vjust = if ((left_angle %||% (if (m > 3) 90 else 0)) == 0) 0.4 else 0.5
      )
    ) +
    ggplot2::labs(x = xlab, y = ylab)
}


# ------------------------------------------------------------------------------
# Internal renderers
# ------------------------------------------------------------------------------

# Hierarchical 2D marimekko: column widths = group total, segments stacked
# vertically with heights = within-group state proportions. Used for mcml,
# where states partition cleanly into clusters so each rectangle's area
# carries unique information.
.plot_marimekko_hierarchical <- function(freq_df, sort_states, colors,
                                          label, label_size,
                                          legend = "bottom",
                                          legend_dir = "auto",
                                          legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  freq_df$state <- factor(freq_df$state, levels = state_levels)
  group_levels <- unique(freq_df$group)
  freq_df$group <- factor(freq_df$group, levels = group_levels)

  pal <- .state_palette(colors, length(state_levels))
  names(pal) <- state_levels

  # Build cumulative-width / cumulative-height rectangle coordinates.
  tab <- tapply(freq_df$count,
                list(freq_df$group, freq_df$state),
                FUN = sum, default = 0)
  tab[is.na(tab)] <- 0
  group_levels <- rownames(tab)
  state_levels <- colnames(tab)
  col_totals <- rowSums(tab)
  total <- sum(col_totals)
  widths_cum <- c(0, cumsum(col_totals)) / total

  rects <- do.call(rbind, lapply(seq_along(group_levels), function(i) {
    row <- tab[i, , drop = TRUE]
    rs  <- sum(row)
    if (rs <= 0) return(NULL)
    heights_cum <- c(0, cumsum(row / rs))
    data.frame(
      group = rep(group_levels[i], length(state_levels)),
      state = state_levels,
      xmin  = widths_cum[i],
      xmax  = widths_cum[i + 1L],
      ymin  = heights_cum[seq_len(length(state_levels))],
      ymax  = heights_cum[seq_len(length(state_levels)) + 1L],
      count = as.numeric(row),
      proportion = as.numeric(row) / rs,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  rects <- rects[rects$ymax > rects$ymin, , drop = FALSE]
  rects$state <- factor(rects$state, levels = state_levels)
  rects$x_mid <- (rects$xmin + rects$xmax) / 2
  rects$y_mid <- (rects$ymin + rects$ymax) / 2
  rects$tile_height <- rects$ymax - rects$ymin

  mid_x <- (widths_cum[-1L] + widths_cum[-length(widths_cum)]) / 2

  legend_ncol <- min(length(state_levels), 5L)

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = mid_x, labels = group_levels,
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0),
                                labels = function(v) sprintf("%d%%", round(v * 100))) +
    ggplot2::labs(x = NULL, y = "Within-cluster proportion", fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count, rects$proportion, label)
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size)
  }
  p
}


# Squarified slice-and-dice treemap: split along the longer side at each
# step, placing the largest remaining rectangle as a full slice. Produces
# rectangles whose AREAS are exactly proportional to `values`, with a
# reasonable aspect ratio for small N (3-9 states).
.simple_treemap <- function(values, x = 0, y = 0, w = 1, h = 1) {
  n <- length(values)
  if (n == 0L) return(data.frame())
  if (n == 1L) {
    return(data.frame(xmin = x, xmax = x + w,
                      ymin = y, ymax = y + h, idx = 1L))
  }
  total <- sum(values)
  share <- values[1L] / total

  if (w >= h) {
    first_w <- w * share
    rest <- .simple_treemap(values[-1L], x + first_w, y, w - first_w, h)
    rest$idx <- rest$idx + 1L
    rbind(
      data.frame(xmin = x, xmax = x + first_w,
                 ymin = y, ymax = y + h, idx = 1L),
      rest
    )
  } else {
    first_h <- h * share
    rest <- .simple_treemap(values[-1L], x, y + first_h, w, h - first_h)
    rest$idx <- rest$idx + 1L
    rbind(
      data.frame(xmin = x, xmax = x + w,
                 ymin = y, ymax = y + first_h, idx = 1L),
      rest
    )
  }
}


# Format a tile / segment label according to user choice. All formats
# render inline on a single line (never two-line).
# "none"  -> ""
# "prop"  -> "66%"
# "freq"  -> "1,234"
# "both"  -> "1,234 (66%)"
# "state" -> "Average"
# "all"   -> "Average (66%)"
# Build the guide_legend layer for a given (position, direction) pair.
# legend_dir = "auto" derives from position (bottom/top -> horizontal,
# left/right -> vertical). "horizontal" / "vertical" force the layout.
.legend_layer <- function(legend, n_states, legend_dir = "auto") {
  if (identical(legend, "none")) {
    return(ggplot2::guides(fill = "none"))
  }
  effective_dir <- if (legend_dir == "auto") {
    if (legend %in% c("bottom", "top")) "horizontal" else "vertical"
  } else legend_dir

  if (effective_dir == "horizontal") {
    ncol_legend <- min(n_states, 5L)
    ggplot2::guides(fill = ggplot2::guide_legend(
      direction = "horizontal", ncol = ncol_legend, byrow = TRUE
    ))
  } else {
    ggplot2::guides(fill = ggplot2::guide_legend(
      direction = "vertical", ncol = 1L
    ))
  }
}


# Theme additions controlling legend position, internal layout, and an
# optional border ("frame") around the legend box.
.legend_theme <- function(legend, legend_dir = "auto", legend_frame = "none") {
  if (identical(legend, "none")) {
    return(ggplot2::theme(legend.position = "none"))
  }
  effective_dir <- if (legend_dir == "auto") {
    if (legend %in% c("bottom", "top")) "horizontal" else "vertical"
  } else legend_dir

  bg_rect <- if (identical(legend_frame, "border")) {
    ggplot2::element_rect(color = "grey40", fill = "white", linewidth = 0.4)
  } else {
    ggplot2::element_blank()
  }

  ggplot2::theme(
    legend.position   = legend,
    legend.box        = effective_dir,
    legend.title      = ggplot2::element_text(face = "bold"),
    legend.background = bg_rect,
    legend.box.background = bg_rect,
    legend.margin     = if (identical(legend_frame, "border"))
                          ggplot2::margin(4, 6, 4, 6) else ggplot2::margin()
  )
}


.format_value_label <- function(state, count, proportion, label) {
  state <- as.character(state)
  freq_str <- format(count, big.mark = ",", trim = TRUE)
  switch(label,
    none  = rep("", length(count)),
    prop  = sprintf("%.1f%%", 100 * proportion),
    freq  = freq_str,
    both  = sprintf("%s (%.1f%%)", freq_str, 100 * proportion),
    state = state,
    all   = sprintf("%s (%.1f%%)", state, 100 * proportion))
}


# Build a single-panel treemap plot with its OWN legend, restricted to
# the states actually present in `sub`. Color is drawn from the global
# named palette so the same state always uses the same color across the
# isolated plots that share a `gridExtra::arrangeGrob` layout.
.single_treemap_plot <- function(sub, pal_named, label, label_size,
                                  legend, legend_dir, legend_frame,
                                  title = NULL, panel_frame = FALSE) {
  states_here <- as.character(sub$state)
  rects <- .simple_treemap(sub$count)
  rects$state      <- states_here[rects$idx]
  rects$count      <- sub$count[rects$idx]
  rects$proportion <- sub$proportion[rects$idx]
  rects$x_mid <- (rects$xmin + rects$xmax) / 2
  rects$y_mid <- (rects$ymin + rects$ymax) / 2
  rects$tile_w <- rects$xmax - rects$xmin
  rects$tile_h <- rects$ymax - rects$ymin
  rects$tile_area <- rects$tile_w * rects$tile_h
  rects$state <- factor(rects$state, levels = states_here)

  pal_local <- pal_named[states_here]

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = pal_local, drop = FALSE,
                               breaks = states_here) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "State") +
    .legend_layer(legend, length(states_here), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )

  if (isTRUE(panel_frame)) {
    p <- p + ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = NA, color = "grey60", linewidth = 0.4),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )
  }

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count,
                                     rects$proportion, label)
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size)
  }
  p
}


# Per-facet renderer: one separate treemap plot per group, each with its
# OWN legend, arranged via gridExtra. Each panel's legend shows only the
# states present in that panel; colors stay consistent across panels via
# a single shared named palette. Returns a gtable.
.plot_per_facet_grid <- function(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend = "right",
                                  legend_dir = "auto",
                                  legend_frame = "border",
                                  combine = TRUE,
                                  ncol = NULL) {
  groups <- unique(as.character(freq_df$group))
  if (length(groups) == 0L) groups <- "all"

  # Per-panel legend placement: right when there are <= 2 panels (each
  # panel still has plenty of horizontal room, e.g. htna AI + Human);
  # bottom otherwise (3+ panels in a combined gtable squeeze each panel
  # too narrow for a right-side legend, which overflows the column).
  if (legend == "per_facet") {
    legend <- if (length(groups) <= 2L) "right" else "bottom"
  }

  # Resolve combine = "auto": when there are 4+ panels, return a list so
  # each panel renders as its own full-size figure under knitr (the chunk
  # `fig.width`/`fig.height` per panel is much more readable than cramming
  # 6+ panels into a 3x2 gtable). 1-3 panels still combine cleanly.
  if (identical(combine, "auto")) {
    combine <- length(groups) <= 3L
  }
  if (isTRUE(combine) && !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("legend = 'per_facet' with combine = TRUE requires the ",
         "'gridExtra' package. Install it, set combine = FALSE, or pick ",
         "a different legend value.", call. = FALSE)
  }

  # Per-panel palettes: each group is independent, so colors restart from
  # the first palette entry inside every panel. This is the right choice
  # when state vocabularies are disjoint (htna actors, mcml clusters) --
  # there is no cross-panel matching to preserve, and starting fresh
  # avoids palette recycling once a panel has many states.
  plots <- lapply(groups, function(g) {
    sub <- freq_df[as.character(freq_df$group) == g, , drop = FALSE]
    local_levels <- .order_states(sub$state, sub$count, sort_states)
    sub <- sub[match(local_levels, as.character(sub$state)), , drop = FALSE]
    sub <- sub[!is.na(sub$state) & sub$count > 0, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)

    local_pal <- .state_palette(colors, length(local_levels))
    names(local_pal) <- local_levels

    .single_treemap_plot(sub, local_pal, label, label_size,
                          legend, legend_dir, legend_frame,
                          title = if (length(groups) > 1L) g else NULL,
                          panel_frame = FALSE)
  })
  plots <- Filter(Negate(is.null), plots)

  if (!isTRUE(combine)) {
    class(plots) <- c("nestimate_facet_list", "list")
    return(plots)
  }

  ncol_arrange <- if (!is.null(ncol)) {
    as.integer(ncol)
  } else if (length(plots) <= 2L) {
    length(plots)
  } else {
    # Per-panel legends are wide, so default to 2 columns for any
    # 3+ panel layout. User can override with ncol = 3 explicitly.
    2L
  }

  gt <- gridExtra::arrangeGrob(grobs = plots, ncol = ncol_arrange)
  class(gt) <- c("nestimate_facet_plot", class(gt))
  gt
}


#' @export
#' @keywords internal
print.nestimate_facet_plot <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x)
  invisible(x)
}


#' @export
#' @keywords internal
print.nestimate_facet_list <- function(x, ...) {
  for (p in x) print(p)
  invisible(x)
}


# knitr's auto-print captures only the last `grid.newpage` from a single
# evaluated expression. To get one rendered figure per panel we must hook
# into knitr::knit_print and emit each plot through its own knit_print
# call (which knitr tracks individually). Falls back to print() outside
# knit context.
#' @rawNamespace if (getRversion() >= "3.6.0") S3method(knitr::knit_print, nestimate_facet_list)
#' @keywords internal
knit_print.nestimate_facet_list <- function(x, ...) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    print(x); return(invisible(NULL))
  }
  parts <- lapply(x, function(p) knitr::knit_print(p, ...))
  knitr::asis_output(paste(unlist(parts), collapse = "\n\n"))
}


# Treemap renderer: one panel per group, each panel a squarified treemap
# whose tile areas are proportional to within-group state counts.
# Single-panel when no groups present.
.plot_treemap_panels <- function(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend = "bottom",
                                  legend_dir = "auto",
                                  legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  pal <- .state_palette(colors, length(state_levels))
  names(pal) <- state_levels

  groups <- unique(as.character(freq_df$group))
  has_groups <- length(groups) > 1L || !identical(groups, "all")

  # Use the SAME global state order in every facet so the leftmost (and
  # second, third, ...) tile always represents the same state across
  # panels. Tile sizes still reflect within-panel proportions, but the
  # topological position of each color stays put -- making cross-facet
  # comparison by color and tile location possible.
  rects_list <- lapply(groups, function(g) {
    sub <- freq_df[as.character(freq_df$group) == g, , drop = FALSE]
    sub <- sub[match(state_levels, as.character(sub$state)), , drop = FALSE]
    sub <- sub[!is.na(sub$state) & sub$count > 0, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    rects <- .simple_treemap(sub$count)
    rects$state      <- as.character(sub$state[rects$idx])
    rects$count      <- sub$count[rects$idx]
    rects$proportion <- sub$proportion[rects$idx]
    rects$group      <- g
    rects$x_mid      <- (rects$xmin + rects$xmax) / 2
    rects$y_mid      <- (rects$ymin + rects$ymax) / 2
    rects$tile_area  <- (rects$xmax - rects$xmin) * (rects$ymax - rects$ymin)
    rects
  })
  rects <- do.call(rbind, Filter(Negate(is.null), rects_list))
  rects$state <- factor(rects$state, levels = state_levels)

  p <- ggplot2::ggplot(rects,
                       ggplot2::aes(xmin = .data$xmin, xmax = .data$xmax,
                                    ymin = .data$ymin, ymax = .data$ymax,
                                    fill = .data$state)) +
    ggplot2::geom_rect(color = "white", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(x = NULL, y = NULL, fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (has_groups) {
    p <- p + ggplot2::facet_wrap(~ group)
  }

  if (label != "none") {
    rects$lab <- .format_value_label(rects$state, rects$count, rects$proportion, label)
    rects$tile_w <- rects$xmax - rects$xmin
    rects$tile_h <- rects$ymax - rects$ymin
    rects$angle <- ifelse(rects$tile_h > rects$tile_w, 90, 0)
    p <- p + .geom_fit_label(rects, label_size)
  }
  p
}




# Horizontal bars: state on Y, count/proportion on X, faceted by group.
# Bar length encodes whichever metric you pass via `metric` ("freq" or
# "prop"). The `label` arg controls the inline numeric annotation.
.plot_state_bars <- function(freq_df, sort_states, colors,
                              label, label_size, metric,
                              legend = "bottom",
                              legend_dir = "auto",
                              legend_frame = "none") {
  state_levels <- .order_states(freq_df$state, freq_df$count, sort_states)
  # Reverse so the largest count appears at the top of the y-axis
  freq_df$state <- factor(freq_df$state, levels = rev(state_levels))
  pal <- .state_palette(colors, length(state_levels))
  names(pal) <- state_levels

  x_var <- if (metric == "freq") "count" else "proportion"
  x_lab <- if (metric == "freq") "Count" else "Proportion"

  freq_df$lab <- .format_value_label(freq_df$state, freq_df$count, freq_df$proportion, label)

  groups <- unique(freq_df$group)
  has_groups <- length(groups) > 1L || !identical(groups, "all")

  p <- ggplot2::ggplot(freq_df,
                       ggplot2::aes(y = .data$state, x = .data[[x_var]],
                                    fill = .data$state)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(values = pal, drop = FALSE,
                               breaks = state_levels) +
    ggplot2::labs(x = x_lab, y = NULL, fill = "State") +
    .legend_layer(legend, length(state_levels), legend_dir) +
    ggplot2::theme_minimal(base_size = 12) +
    .legend_theme(legend, legend_dir, legend_frame) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())

  if (has_groups) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x")
  }

  if (label != "none") {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$lab), hjust = -0.15, size = label_size,
      color = "grey20"
    ) + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.28)))
  }
  p
}


# ------------------------------------------------------------------------------
# Generic + methods (exported)
# ------------------------------------------------------------------------------

#' Plot State Frequency Distributions
#'
#' Visualise state (node) frequency distributions across groups for any
#' Nestimate object that carries sequence data: a single \code{netobject},
#' a \code{netobject_group}, an \code{mcml} model, or an \code{htna} network.
#'
#' The marimekko layout is dispatched per class:
#' \itemize{
#'   \item For \code{mcml}, where states partition cleanly into clusters,
#'     the chart is a hierarchical 2D marimekko: cluster columns of width
#'     proportional to cluster total, segments stacked vertically with
#'     heights proportional to within-cluster state proportions.
#'   \item For all other classes (\code{netobject}, \code{netobject_group},
#'     \code{htna}), each group is rendered as its own panel containing a
#'     squarified treemap: each state becomes a rectangular tile whose
#'     AREA is exactly proportional to the state's share within that group.
#'     Single-panel when no groups exist; faceted when groups are present.
#' }
#'
#' The bar style produces horizontal bars (state on the y-axis), faceted by
#' group when groups exist. All variants use the Okabe-Ito palette.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml}, or
#'   \code{htna} object.
#' @param style One of:
#'   \itemize{
#'     \item \code{"marimekko"} (default) -- per-group treemap panels with
#'       cumulative-width geometry; tile area = within-group state share.
#'     \item \code{"bars"} -- horizontal bars sorted by frequency, faceted
#'       per group.
#'   }
#'   For chi-square mosaics of a (group x state) contingency table, use
#'   \code{\link{mosaic_plot}} directly -- it is kept as a separate
#'   function with its own dispatch surface.
#' @param metric For \code{style = "bars"}: which value the bar length
#'   encodes -- \code{"prop"} (default) or \code{"freq"}. Treemap and
#'   hierarchical-marimekko areas always encode proportion within group.
#' @param label Inline tile / bar annotation. All formats render on a
#'   single line.
#'   \itemize{
#'     \item \code{"prop"} (default) -- proportion only, e.g. \code{"66\%"}
#'     \item \code{"freq"} -- count only, e.g. \code{"1,234"}
#'     \item \code{"both"} -- count + proportion, e.g. \code{"1,234 (66\%)"}
#'     \item \code{"state"} -- state name only, e.g. \code{"Average"}
#'     \item \code{"all"} -- state + proportion, e.g. \code{"Average (66\%)"}
#'     \item \code{"none"} -- no inline labels
#'   }
#' @param legend Legend position. \code{"auto"} (default) resolves per
#'   style: \code{"none"} for \code{style = "bars"} (the y-axis already
#'   names every state, so a colour legend is redundant);
#'   \code{"per_facet"} for \code{htna}/\code{mcml} treemaps (state
#'   vocabularies differ per panel, so each gets its own legend);
#'   \code{"bottom"} for single-network and \code{netobject_group}
#'   treemaps (shared state vocabulary, one shared legend).
#'   Override with any of \code{"bottom"}, \code{"top"}, \code{"right"},
#'   \code{"left"}, \code{"none"}, or \code{"per_facet"}. The
#'   \code{"per_facet"} option requires the \pkg{gridExtra} package and
#'   returns a \code{gtable}.
#' @param legend_dir Legend internal layout: \code{"auto"} (default --
#'   horizontal for top/bottom, vertical for left/right), or force
#'   \code{"horizontal"} or \code{"vertical"} regardless of position.
#' @param legend_frame \code{"none"} (default) for an unframed legend, or
#'   \code{"border"} to draw a thin grey rectangle around the legend
#'   ("legend enclosed in a square").
#' @param sort_states One of \code{"frequency"} (default -- most frequent
#'   first), \code{"alpha"}, or \code{"none"}.
#' @param colors Optional character vector overriding the default
#'   Okabe-Ito state palette. Length must be at least the number of unique
#'   states.
#' @param label_size Numeric size of inline labels (max size when
#'   \pkg{ggfittext} is installed -- text auto-shrinks per tile).
#' @param abbreviate Abbreviate state names. \code{FALSE} (default) shows
#'   full names; \code{TRUE} truncates to the first 3 characters via
#'   \code{base::abbreviate()} (which extends the truncation as needed to
#'   keep names unique after collision); a positive integer sets the
#'   target minimum length explicitly (e.g. \code{abbreviate = 4}).
#'   Affects tile labels, legend, and the returned \code{$table}.
#' @param include_macro For \code{mcml} only: prepend a \code{"macro"}
#'   reference column showing aggregate state frequencies across all
#'   clusters. Default \code{FALSE}.
#' @param combine For \code{legend = "per_facet"} only. \code{"auto"}
#'   (default) returns a single combined gtable for 1-3 panels and a
#'   list of ggplots (one per panel) for 4+ panels -- many-cluster
#'   \code{mcml} layouts read better as separate figures than as a tile
#'   grid. \code{TRUE} forces a combined gtable via \pkg{gridExtra};
#'   \code{FALSE} forces a list (knitr renders each at the chunk's full
#'   \code{fig.width} / \code{fig.height}).
#' @param node_groups Optional named character vector mapping node labels to
#'   semantic groups. When supplied, panels (or bars) are coloured / annotated
#'   by group rather than by individual state, so state-level palettes can
#'   collapse onto a smaller categorical legend.
#' @param ncol For \code{legend = "per_facet"} with \code{combine = TRUE}:
#'   number of columns in the grid arrangement. \code{NULL} (default)
#'   picks 1, 2, or 3 columns based on the number of panels.
#' @param ... Reserved for future use.
#'
#' @return A \code{state_freq} object: a list with the rendered \code{$plot}
#'   (a \code{ggplot} or \code{gtable}), the tidy \code{$table} (a
#'   \code{data.frame} with columns \code{group}, \code{state}, \code{count},
#'   \code{proportion}), and the call's \code{$style}, \code{$metric},
#'   \code{$source_class}. The class supports \code{print()} (shows the
#'   tidy table in the console), \code{plot()} (renders the chart), and
#'   \code{as.data.frame()} (returns the table).
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   data(group_regulation_long, package = "Nestimate")
#'   nw <- build_network(group_regulation_long,
#'                       method = "relative", format = "long",
#'                       actor = "Actor", action = "Action",
#'                       order = "Time", group = "Course")
#'   res <- plot_state_frequencies(nw)
#'   print(res)            # tidy frequency table in the console
#'   plot(res)             # ggplot chart
#'   head(as.data.frame(res))
#' }
#' }
#' @export
plot_state_frequencies <- function(x, ...) {
  UseMethod("plot_state_frequencies")
}


#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject <- function(x, legend = "auto", ...) {
  .plot_state_frequencies_impl(x, source_class = "netobject",
                               hierarchical = FALSE, legend = legend, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.htna <- function(x, legend = "auto", ...) {
  .plot_state_frequencies_impl(x, source_class = "htna",
                               hierarchical = FALSE, legend = legend, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.mcml <- function(x, legend = "auto", ...) {
  .plot_state_frequencies_impl(x, source_class = "mcml",
                               hierarchical = TRUE, legend = legend, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.netobject_group <- function(x, legend = "auto", ...) {
  .plot_state_frequencies_impl(x, source_class = "netobject_group",
                               hierarchical = FALSE, legend = legend, ...)
}

#' @export
#' @rdname plot_state_frequencies
plot_state_frequencies.default <- function(x, ...) {
  cls <- paste(class(x), collapse = "/")
  stop("plot_state_frequencies(): no method for class '", cls, "'.\n",
       "Supported: netobject, netobject_group, mcml, htna.\n",
       "For tna objects, use tna::plot_frequencies(x).",
       call. = FALSE)
}


# Single worker for all four dispatch methods. Validates args once, pulls
# the tidy frame via state_distribution(), runs the existing renderer
# pipeline, wraps both into a state_freq object.
.plot_state_frequencies_impl <- function(x, source_class, hierarchical,
                                          style       = c("marimekko",
                                                          "bars"),
                                          metric      = c("prop", "freq"),
                                          label       = c("prop", "freq",
                                                          "both", "state",
                                                          "all", "none"),
                                          legend      = c("auto", "bottom",
                                                          "right", "top",
                                                          "left", "none",
                                                          "per_facet"),
                                          legend_dir  = c("auto",
                                                          "horizontal",
                                                          "vertical"),
                                          legend_frame = c("none", "border"),
                                          sort_states = c("frequency",
                                                          "alpha", "none"),
                                          colors      = NULL,
                                          label_size  = 3.5,
                                          abbreviate  = FALSE,
                                          include_macro = FALSE,
                                          combine     = "auto",
                                          ncol        = NULL,
                                          node_groups = NULL,
                                          ...) {
  style        <- match.arg(style)
  metric       <- match.arg(metric)
  label        <- match.arg(label)
  legend       <- match.arg(legend)
  legend_dir   <- match.arg(legend_dir)
  legend_frame <- match.arg(legend_frame)
  sort_states  <- match.arg(sort_states)

  if (identical(legend, "auto")) {
    legend <- if (style == "bars") {
      "none"
    } else if (source_class %in% c("htna", "mcml")) {
      "per_facet"
    } else {
      "bottom"
    }
  }

  freq_df <- if (identical(source_class, "mcml")) {
    state_distribution(x, include_macro = include_macro)
  } else {
    state_distribution(x)
  }

  freq_df <- .abbreviate_states(freq_df, abbreviate)

  # Demote per_facet to a shared bottom legend when every group has the
  # same state vocabulary -- repeating the same legend in every panel is
  # pure redundancy.
  if (identical(legend, "per_facet") &&
      .vocab_is_shared(freq_df)) {
    legend <- "bottom"
  }

  p <- .render_freq(freq_df, style, hierarchical = hierarchical,
                    metric, label, sort_states, colors, label_size,
                    legend, legend_dir, legend_frame,
                    combine = combine, ncol = ncol,
                    node_groups = node_groups)

  .new_state_freq(plot = p, table = freq_df,
                  style = style, metric = metric,
                  source_class = source_class)
}


# Returns TRUE when every group in freq_df contains the same set of
# states (any order). Used to detect per_facet calls that would render
# the same legend k times.
.vocab_is_shared <- function(freq_df) {
  groups <- as.character(freq_df$group)
  if (length(unique(groups)) < 2L) return(FALSE)
  vocabs <- split(as.character(freq_df$state), groups)
  vocabs <- lapply(vocabs, function(v) sort(unique(v)))
  all(vapply(vocabs[-1L], identical, logical(1L), vocabs[[1L]]))
}


# Abbreviate state names with base R `abbreviate()` so duplicates after
# truncation get extra letters automatically. `abbreviate = FALSE` (default)
# is a no-op; `TRUE` uses minlength = 3; a numeric value sets minlength.
# Aggregates the freq_df after truncation in case different long names
# collapse to the same short name despite the uniqueness guarantee
# (rare; happens with single-letter prefixes or after manual relabel).
.abbreviate_states <- function(freq_df, abbreviate) {
  if (isFALSE(abbreviate) || is.null(abbreviate)) return(freq_df)
  minlen <- if (isTRUE(abbreviate)) 3L else as.integer(abbreviate[[1L]])
  stopifnot(minlen >= 1L)
  short <- abbreviate(as.character(freq_df$state), minlength = minlen)
  freq_df$state <- unname(short)
  agg <- aggregate(count ~ group + state, data = freq_df, FUN = sum)
  totals <- ave(agg$count, agg$group, FUN = sum)
  agg$proportion <- agg$count / totals
  agg[, c("group", "state", "count", "proportion")]
}


# ------------------------------------------------------------------------------
# state_distribution(): public tidy frequency extractor
# ------------------------------------------------------------------------------

#' Per-Class State Distribution as a Tidy Data Frame
#'
#' Returns a tidy \code{data.frame(group, state, count, proportion)} with one
#' row per (group, state) cell. Companion to \code{\link{state_frequencies}}
#' (which counts unique states in raw sequence input);
#' \code{state_distribution()} pulls the same shape of frame from a fitted
#' Nestimate object so analyses don't have to reach for the underlying
#' \code{$data} slot directly.
#'
#' Used internally by \code{\link{plot_state_frequencies}} as the data layer
#' behind every chart, and surfaced as the \code{$table} slot of the
#' returned \code{state_freq} object.
#'
#' @param x A \code{netobject}, \code{netobject_group}, \code{mcml}, or
#'   \code{htna} object.
#' @param include_macro For \code{mcml}: when \code{TRUE}, prepend a
#'   \code{group = "macro"} block aggregating across clusters. Ignored for
#'   the other classes.
#' @param ... Currently unused.
#'
#' @return A \code{data.frame} with columns \code{group} (character),
#'   \code{state} (character), \code{count} (integer), and
#'   \code{proportion} (numeric, within-group share).
#' @export
#' @examples
#' \dontrun{
#'   data(ai_long)
#'   net <- build_network(ai_long, method = "frequency",
#'                        id_col = "session_id",
#'                        time_col = "order_in_session", action = "code")
#'   state_distribution(net)
#' }
state_distribution <- function(x, ...) UseMethod("state_distribution")

#' @export
#' @rdname state_distribution
state_distribution.netobject <- function(x, ...) .freq_df_netobject(x)

#' @export
#' @rdname state_distribution
state_distribution.htna <- function(x, ...) .freq_df_htna(x)

#' @export
#' @rdname state_distribution
state_distribution.mcml <- function(x, include_macro = FALSE, ...) {
  .freq_df_mcml(x, include_macro = include_macro)
}

#' @export
#' @rdname state_distribution
state_distribution.netobject_group <- function(x, ...) {
  .freq_df_netobject_group(x)
}

#' @export
#' @rdname state_distribution
state_distribution.default <- function(x, ...) {
  cls <- paste(class(x), collapse = "/")
  stop("state_distribution(): no method for class '", cls, "'.\n",
       "Supported: netobject, netobject_group, mcml, htna.",
       call. = FALSE)
}


# ------------------------------------------------------------------------------
# state_freq S3 class: holds plot + tidy table together
# ------------------------------------------------------------------------------

.new_state_freq <- function(plot, table, style, metric, source_class) {
  out <- list(
    plot         = plot,
    table        = table,
    style        = style,
    metric       = metric,
    source_class = source_class
  )
  class(out) <- "state_freq"
  out
}

#' Print, Plot, and Convert a state_freq Object
#'
#' \code{plot_state_frequencies()} returns a \code{state_freq} object holding
#' both the rendered chart and the tidy frequency table. \code{print()} shows
#' the table in the console, \code{plot()} renders the chart, and
#' \code{as.data.frame()} returns the tidy table for downstream piping.
#'
#' @param x A \code{state_freq} object.
#' @param digits Number of decimal places for proportion / share columns.
#' @param max_states Cap on rows shown per group in the per-state table.
#'   The full table remains available via \code{x$table}.
#' @param ... Unused.
#' @return \code{print()} returns \code{invisible(x)}; \code{plot()} returns
#'   \code{invisible(NULL)} after drawing; \code{as.data.frame()} returns
#'   \code{x$table}.
#' @name state_freq
NULL

#' @export
#' @rdname state_freq
print.state_freq <- function(x, digits = 1, max_states = 20L, ...) {
  tbl <- x$table
  groups <- unique(as.character(tbl$group))
  total_events <- sum(tbl$count)
  n_states <- length(unique(tbl$state))
  pct <- function(p) sprintf(paste0("%.", digits, "f%%"), 100 * p)

  cat(sprintf(
    "State frequencies (style = %s, source = %s)\n",
    x$style, x$source_class
  ))
  cat(sprintf(
    "  Total events: %s  |  Groups: %d  |  States: %d\n\n",
    format(total_events, big.mark = ","), length(groups), n_states
  ))

  group_totals <- tapply(tbl$count, tbl$group, sum)
  group_totals <- group_totals[groups]
  totals_lines <- .cluster_table_lines(list(
    group  = groups,
    events = format(as.integer(group_totals), big.mark = ","),
    share  = pct(as.numeric(group_totals) / total_events)
  ))
  cat("Per-group totals\n")
  cat(paste(paste0("  ", totals_lines), collapse = "\n"), "\n\n")

  cat("Per-state proportions (within group)\n")
  parts <- lapply(groups, function(g) {
    sub <- tbl[as.character(tbl$group) == g, , drop = FALSE]
    sub <- sub[order(-sub$count), , drop = FALSE]
    if (nrow(sub) > max_states) {
      kept <- sub[seq_len(max_states), , drop = FALSE]
      kept$state <- as.character(kept$state)
      kept <- rbind(kept, data.frame(
        group = g, state = sprintf("(+%d more)", nrow(sub) - max_states),
        count = sum(sub$count[-seq_len(max_states)]),
        proportion = sum(sub$proportion[-seq_len(max_states)]),
        stringsAsFactors = FALSE
      ))
      sub <- kept
    }
    sub
  })
  combined <- do.call(rbind, parts)
  body_lines <- .cluster_table_lines(list(
    group = as.character(combined$group),
    state = as.character(combined$state),
    count = format(as.integer(combined$count), big.mark = ","),
    share = pct(combined$proportion)
  ))
  cat(paste(paste0("  ", body_lines), collapse = "\n"), "\n")

  # Render the chart to the active graphics device so the table + plot
  # appear together (console: opens a window; knitr: embeds inline).
  if (inherits(x$plot, "ggplot")) {
    print(x$plot)
  } else if (inherits(x$plot, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(x$plot)
  }
  invisible(x)
}

#' @export
#' @rdname state_freq
plot.state_freq <- function(x, ...) {
  if (inherits(x$plot, "ggplot")) {
    print(x$plot)
  } else if (inherits(x$plot, "gtable")) {
    grid::grid.newpage()
    grid::grid.draw(x$plot)
  } else {
    print(x$plot)
  }
  invisible(NULL)
}

#' @export
#' @rdname state_freq
as.data.frame.state_freq <- function(x, ...) x$table


# Shared dispatcher.
# `hierarchical = TRUE` (mcml) -> 2D cumulative marimekko.
# `hierarchical = FALSE` (everything else) -> per-panel squarified treemap
# (single panel when only one group exists).
.render_freq <- function(freq_df, style, hierarchical,
                          metric, label, sort_states, colors, label_size,
                          legend, legend_dir, legend_frame,
                          combine = TRUE, ncol = NULL,
                          node_groups = NULL) {
  if (nrow(freq_df) == 0L) {
    stop("No state observations available to plot.", call. = FALSE)
  }

  # Per-facet legend mode applies to the marimekko/treemap path only.
  # For style = "bars", per_facet doesn't make sense (bars already
  # facet-wrap with shared legend), so fall back to a right-side legend.
  if (identical(legend, "per_facet")) {
    if (style == "bars") {
      return(.plot_state_bars(freq_df, sort_states, colors,
                              label, label_size, metric,
                              "right", legend_dir, legend_frame))
    }
    return(.plot_per_facet_grid(freq_df, sort_states, colors,
                                 label, label_size,
                                 legend = "per_facet",
                                 legend_dir = legend_dir,
                                 legend_frame = legend_frame,
                                 combine = combine,
                                 ncol = ncol))
  }

  if (style == "bars") {
    return(.plot_state_bars(freq_df, sort_states, colors,
                            label, label_size, metric,
                            legend, legend_dir, legend_frame))
  }

  if (isTRUE(hierarchical)) {
    .plot_marimekko_hierarchical(freq_df, sort_states, colors,
                                  label, label_size,
                                  legend, legend_dir, legend_frame)
  } else {
    .plot_treemap_panels(freq_df, sort_states, colors,
                          label, label_size,
                          legend, legend_dir, legend_frame)
  }
}
