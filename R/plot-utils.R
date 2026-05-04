# ==============================================================================
# Internal plotting helpers shared by sequence_plot() and distribution_plot().
# Not exported. Keep this file thin: only logic used in 2+ places belongs here.
# ==============================================================================

.okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")


# Encode a wide sequence matrix (states as character/factor/integer) into
# integer codes and return the sorted level set. Returns a list with:
#   codes      - integer matrix, same shape as input, NA preserved
#   levels     - sort.int(unique(non-NA))
.encode_states <- function(x) {
  stopifnot(is.data.frame(x) || is.matrix(x))
  if (is.data.frame(x)) x <- as.matrix(x)
  levels <- sort.int(unique(as.vector(x)), na.last = NA)
  if (length(levels) == 0L) {
    stop("`x` contains no non-NA values.", call. = FALSE)
  }
  codes <- matrix(match(x, levels),
                  nrow = nrow(x), ncol = ncol(x),
                  dimnames = dimnames(x))
  list(codes = codes, levels = levels)
}


# Compute only the dissimilarity matrix for a sequence of states using the
# same dispatch as build_clusters(), but without running PAM / hclust /
# silhouette. Used by sequence_plot() when sort = any distance metric.
.sequence_distance <- function(x, dissimilarity = "lcs") {
  enc <- .encode_sequences(x)
  use_stringdist <- requireNamespace("stringdist", quietly = TRUE) &&
    enc$n_states < 62L &&
    !(dissimilarity == "hamming")
  if (use_stringdist) {
    .dissimilarity_matrix_stringdist(enc, dissimilarity, lambda = 0)
  } else {
    .dissimilarity_matrix_r(enc, dissimilarity, lambda = 0)
  }
}


# Translate the legend knobs into outer-margin sizes (in text-lines). Uses
# the device's actual char-width / line-height ratio so long labels or
# titles don't clip. Returns a 2-vector c(oma_bottom, oma_right); all
# other edges get 0.3 elsewhere.
.legend_oma_size <- function(labels, position, size, ncol, title) {
  label_chars <- max(nchar(as.character(labels)))
  title_chars <- if (!is.null(title)) nchar(title) else 0
  wide_chars  <- max(label_chars, title_chars)
  cw_per_line <- graphics::par("cin")[1] / graphics::par("csi")
  rows <- if (!is.null(ncol) && ncol > 0) {
    ceiling(length(labels) / ncol)
  } else if (position == "bottom") 1L else length(labels)
  cols <- if (!is.null(ncol) && ncol > 0) {
    ncol
  } else if (position == "bottom") length(labels) else 1L
  # text line factor 1.6 (was 1.25): reserves enough height that cex = 1.2
  # descenders + padding fit inside the oma strip without clipping at the
  # device bottom edge.
  h <- rows * size * 1.6 + (if (!is.null(title)) 1.2 else 0.6)
  w <- cols * ((wide_chars + 3) * cw_per_line * size * 1.15) + 2.5
  c(oma_b = if (position == "bottom") h else 0.3,
    oma_r = if (position == "right")  w else 0.3)
}


# Draw a legend anchored to the inner-figure edge, growing outward into
# the reserved oma strip. Called after all plot panels are drawn.
.draw_legend_in_oma <- function(labels, colors, position, size, ncol,
                                title, border, bty) {
  args <- list(legend = as.character(labels),
               fill   = colors,
               border = border,
               bty    = bty,
               cex    = size,
               pt.cex = size * 1.3,
               title  = title,
               xpd    = NA)
  if (is.null(ncol)) {
    args$horiz <- position == "bottom"
  } else {
    args$ncol <- ncol
  }
  omd <- graphics::par("omd")
  gap <- 0.005
  anchor <- switch(position,
                   bottom = list(x = 0.5,            y = omd[3] - gap,
                                 xj = 0.5,           yj = 1),
                   right  = list(x = omd[2] + gap,   y = 0.5,
                                 xj = 0,             yj = 0.5))
  x <- graphics::grconvertX(anchor$x, from = "ndc", to = "user")
  y <- graphics::grconvertY(anchor$y, from = "ndc", to = "user")
  do.call(graphics::legend,
          c(list(x = x, y = y, xjust = anchor$xj, yjust = anchor$yj), args))
}


# Every-Nth tick positions. When `tick` is NULL, thin automatically so at
# most ~15 labels are drawn.
.tick_positions <- function(n_cols, tick, colnames_ = NULL) {
  by <- if (!is.null(tick)) {
    stopifnot(is.numeric(tick), length(tick) == 1L, tick >= 1)
    as.integer(tick)
  } else {
    max(1L, as.integer(ceiling(n_cols / 15)))
  }
  at  <- seq(1L, n_cols, by = by)
  lab <- if (!is.null(colnames_)) colnames_[at] else at
  list(at = at, labels = lab)
}


# Adaptive legend text size. NULL = auto-scale from device width so the
# legend visually matches the plot: ~0.7 on a 5-inch figure, ~1.0 on a
# 10-inch figure, clamped to [0.65, 1.2] to stay readable.
.auto_legend_size <- function(user_size) {
  if (!is.null(user_size)) return(user_size)
  din_w <- graphics::par("din")[1]
  max(0.65, min(1.2, din_w / 10))
}


# Build the state palette from an optional user-supplied vector, falling
# back to recycled Okabe-Ito. Returns a character vector length = n_states.
.state_palette <- function(user_colors, n_states) {
  if (is.null(user_colors)) return(rep_len(.okabe_ito, n_states))
  stopifnot(is.character(user_colors), length(user_colors) >= n_states)
  user_colors[seq_len(n_states)]
}


# ==============================================================================
# Plot-family classification.
#
# The "data-bearing" classes are the four kinds of object the Nestimate plot
# family operates on: netobject, netobject_group, mcml, htna. Every plot-family
# verb (plot_state_frequencies, mosaic_plot, state_distribution, ...) must
# either accept one of these and dispatch through the shared *_impl worker, or
# error with a default method that names this set.
#
# .is_data_bearing(x) -> single logical
# .data_bearing_class(x) -> one of "netobject" / "netobject_group" / "mcml" /
#   "htna", or NA_character_ when not data-bearing. The order matters: htna
#   inherits from netobject, so check htna first.
# ==============================================================================

.data_bearing_classes <- c("htna", "mcml", "netobject_group", "netobject")

.is_data_bearing <- function(x) {
  any(.data_bearing_classes %in% class(x))
}

.data_bearing_class <- function(x) {
  cls <- class(x)
  hit <- .data_bearing_classes[.data_bearing_classes %in% cls]
  if (length(hit) == 0L) return(NA_character_)
  hit[[1L]]
}


# ==============================================================================
# Shared layout helpers used across plot-family verbs (plot_state_frequencies,
# mosaic_plot, ...). All operate on tidy ggplot inputs; none assume a specific
# data-bearing source class.
# ==============================================================================

# Fit-aware tile labels. Uses ggfittext::geom_fit_text() when available so each
# tile's text auto-shrinks (and reflows onto multiple lines) to fit its
# rectangle, with min.size dropping the label when nothing legible fits. When
# ggfittext is not installed, falls back to ggplot2::geom_text() at midpoint.
# Caller is expected to have already nulled out labels for tiles too small for
# legible rendering.
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


# TRUE when every group in freq_df shares the same state set (any order).
# Used by per-facet legend logic to demote per_facet -> bottom when repeating
# the same legend in every panel would be pure redundancy.
.vocab_is_shared <- function(freq_df) {
  groups <- as.character(freq_df$group)
  if (length(unique(groups)) < 2L) return(FALSE)
  vocabs <- split(as.character(freq_df$state), groups)
  vocabs <- lapply(vocabs, function(v) sort(unique(v)))
  all(vapply(vocabs[-1L], identical, logical(1L), vocabs[[1L]]))
}


# Build the guide_legend layer for a given (position, direction) pair.
# legend_dir = "auto" derives from position (bottom/top -> horizontal, else
# vertical). "horizontal" / "vertical" force the layout.
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
