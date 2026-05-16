# ==============================================================================
# sequence_plot() - unified entry point for categorical sequence plots.
# type = "heatmap"      - dense carpet with optional dendrogram (single panel).
# type = "index"        - row-gapped carpet without dendrogram; supports group
#                         + ncol/nrow facet grid (like distribution_plot).
# type = "distribution" - dispatches to distribution_plot().
# ==============================================================================

# ---- Input extraction helper ------------------------------------------------
#' Extract sequence data and clustering info from various input types
#'
#' Supports: data.frame, matrix, net_clustering, netobject, netobject_group,
#' net_mmm. Returns a list with data, group (assignments), and clustering info.
#' @noRd
.extract_seqplot_input <- function(x, group = NULL) {
  clustering <- NULL
  data <- NULL


  # --- net_clustering (from build_clusters) ---
  if (inherits(x, "net_clustering")) {
    data <- x$data
    if (is.null(group)) group <- x$assignments
    clustering <- x
  }
  # --- net_mmm (from build_mmm) ---
  else if (inherits(x, "net_mmm")) {
    data <- x$data %||% x$models[[1L]]$data
    if (is.null(group)) group <- x$assignments
    # No distance matrix in MMM; clustering info for reference only
    clustering <- list(assignments = x$assignments, k = x$k)
  }
  # --- netobject_group (from cluster_network / build_network on clustering) ---
  else if (inherits(x, "netobject_group")) {
    cl <- attr(x, "clustering")
    # Prefer the full N-row data carried by the clustering attribute --
    # cluster_network() splits per-cluster $data subsets into each member,
    # so x[[1L]]$data has only n1 rows whereas $assignments has all N.
    # cluster_mmm() now stashes the same full data to make this invariant.
    if (!is.null(cl) && !is.null(cl$data)) {
      data <- cl$data
      if (is.null(group)) group <- cl$assignments
      clustering <- cl
    } else {
      # Fallback: no clustering attribute (plain group_col split). Rebuild
      # the full sequence frame from per-group subsets and label rows by
      # their group membership.
      parts <- lapply(seq_along(x), function(i) x[[i]]$data)
      if (any(vapply(parts, is.null, logical(1L)))) {
        stop("netobject_group members have no $data; rebuild with ",
             "build_network() carrying data, or pass the original ",
             "clustering object instead.", call. = FALSE)
      }
      part_nrow <- vapply(parts, NROW, integer(1L))
      data <- do.call(rbind, parts)
      if (is.null(group)) {
        labs  <- if (!is.null(names(x))) names(x) else as.character(seq_along(x))
        group <- factor(rep(labs, part_nrow), levels = labs)
      }
    }
  }
  # --- netobject (single network) ---
  else if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    if (is.null(x$data)) {
      stop("netobject has no $data. Pass the original sequence data.",
           call. = FALSE)
    }
    data <- x$data
  }
  # --- tna model ---
  else if (inherits(x, "tna")) {
    if (is.null(x$data) || is.null(x$labels)) {
      stop("tna object missing $data or $labels.", call. = FALSE)
    }
    # Decode integer matrix back to state names
    decoded <- matrix(x$labels[x$data], nrow = nrow(x$data),
                      ncol = ncol(x$data))
    colnames(decoded) <- colnames(x$data)
    data <- as.data.frame(decoded, stringsAsFactors = FALSE)
  }
  # --- data.frame or matrix (pass through) ---
  else if (is.data.frame(x) || is.matrix(x)) {
    data <- x
  }
  else {
    stop("Unsupported input type: ", paste(class(x), collapse = ", "),
         ". Expected data.frame, matrix, netobject, netobject_group, ",
         "net_clustering, net_mmm, or tna.", call. = FALSE)
  }

  list(data = data, group = group, clustering = clustering)
}

#' Sequence Plot (heatmap, index, or distribution)
#'
#' Single entry point for three categorical-sequence visualisations.
#' \itemize{
#'   \item \code{type = "heatmap"} (default): dense carpet, rows reordered
#'     by \code{sort} / dendrogram (single panel).
#'   \item \code{type = "index"}: same data layout, but rows separated by
#'     thin gaps (no dendrogram). Supports grouping via \code{group} or a
#'     \code{net_clustering}, plus a \code{ncol} x \code{nrow} facet grid.
#'   \item \code{type = "distribution"}: dispatches to
#'     \code{\link{distribution_plot}}.
#' }
#'
#' @param x Wide-format sequence data. Accepts:
#'   \describe{
#'     \item{data.frame / matrix}{Rows = sequences, columns = time points.}
#'     \item{netobject}{Extracts \code{$data}.}
#'     \item{net_clustering}{From \code{\link{build_clusters}}. Uses
#'       \code{$data}, \code{$assignments} for grouping, and \code{$distance}
#'       for dendrogram.}
#'     \item{netobject_group}{From \code{\link{cluster_network}} or
#'       \code{\link{build_network}} on a clustering. Extracts data and
#'       assignments from \code{attr(, "clustering")}.
#'     }
#'     \item{net_mmm}{From \code{\link{build_mmm}}. Uses \code{$models[[1]]$data}
#'       and \code{$assignments}.}
#'     \item{tna}{From the tna package. Decodes integer-encoded sequences.}
#'   }
#' @param type One of \code{"heatmap"} (default), \code{"index"}, or
#'   \code{"distribution"}.
#' @param sort Row-ordering strategy for heatmap / within-panel for index.
#'   One of \code{"lcs"} (default), \code{"frequency"}, \code{"start"},
#'   \code{"end"}, or any \code{\link{build_clusters}} distance
#'   (\code{"hamming"}, \code{"osa"}, \code{"lv"}, \code{"dl"},
#'   \code{"qgram"}, \code{"cosine"}, \code{"jaccard"}, \code{"jw"}).
#' @param tree Optional \code{hclust}/\code{dendrogram}/\code{agnes}
#'   object to supply row ordering (heatmap only; overrides \code{sort}).
#' @param group Optional grouping vector (length \code{nrow(x)}) producing
#'   one facet per group. Index/distribution only. Ignored for heatmap.
#' @param scale,geom,na Passed to \code{\link{distribution_plot}} when
#'   \code{type = "distribution"}.
#' @param row_gap Fraction of row height used as vertical gap between
#'   sequences in index plots. \code{0} (default) = dense like heatmap.
#'   Try \code{0.15} for visible separators at low row counts.
#' @param dendrogram_width Width ratio of the dendrogram panel (heatmap).
#' @param k Optional integer. When supplied in \code{type = "heatmap"},
#'   cuts the dendrogram into \code{k} clusters and draws thin horizontal
#'   separators between them in the carpet. Ignored when there is no
#'   dendrogram (e.g. \code{sort = "start"}) or for other types.
#' @param k_color Colour for the cluster separator lines. Default
#'   \code{"white"}.
#' @param k_line_width Line width for the cluster separators. Default
#'   \code{2.5}.
#' @param state_colors Vector of colours, one per state.
#' @param na_color Colour for \code{NA} cells.
#' @param cell_border Cell border colour. \code{NA} = off.
#' @param frame If \code{TRUE} (default), draw a box around each panel.
#'   If \code{FALSE}, no box - axis ticks and labels still appear.
#' @param width,height Optional device dimensions in inches. When supplied,
#'   opens a new graphics device via \code{grDevices::dev.new()}. In knitr
#'   chunks use the \code{fig.width} / \code{fig.height} chunk options
#'   instead.
#' @param main Plot title.
#' @param show_n Append \code{"(n = N)"} to the title.
#' @param time_label,xlab X-axis label. \code{xlab} is an alias.
#' @param y_label,ylab Y-axis label (distribution only). \code{ylab} alias.
#' @param tick Show every Nth x-axis label. \code{NULL} = auto.
#' @param ncol,nrow Facet grid dimensions (index + distribution).
#'   Ignored when \code{combined = FALSE}.
#' @param combined Index and distribution types only. When \code{TRUE}
#'   (default), groups are arranged on one figure via
#'   \code{graphics::layout()}. When \code{FALSE}, each group is drawn
#'   on its own page (one full-size figure per group, with its own
#'   legend). Single-group calls (\code{G == 1}) ignore this argument.
#'   Heatmap is always single-figure.
#' @param legend Legend position: \code{"bottom"}, \code{"right"}, or
#'   \code{"none"}. Default varies by type.
#' @param legend_size Legend text size. \code{NULL} (default) auto-scales
#'   from the device width so the legend looks proportional at 5 in vs
#'   12 in figures (clamped to \code{[0.65, 1.2]}).
#' @param legend_title Optional legend title.
#' @param legend_ncol Number of legend columns.
#' @param legend_border Swatch border colour.
#' @param legend_bty \code{"n"} or \code{"o"}.
#'
#' @return Invisibly, a list describing the plot (shape depends on
#'   \code{type}).
#' @seealso \code{\link{distribution_plot}}, \code{\link{build_clusters}}
#' @examples
#' \donttest{
#' sequence_plot(trajectories)
#' sequence_plot(trajectories, type = "index")
#' sequence_plot(trajectories, type = "distribution")
#' }
#' @export
sequence_plot <- function(x,
                          type             = c("heatmap", "index",
                                               "distribution"),
                          sort             = c("lcs", "frequency",
                                               "start", "end",
                                               "hamming", "osa", "lv", "dl",
                                               "qgram", "cosine",
                                               "jaccard", "jw"),
                          tree             = NULL,
                          group            = NULL,
                          scale            = c("proportion", "count"),
                          geom             = c("area", "bar"),
                          na               = TRUE,
                          row_gap          = 0,
                          dendrogram_width = 1.2,
                          k                = NULL,
                          k_color          = "white",
                          k_line_width     = 2.5,
                          state_colors     = NULL,
                          na_color         = "grey90",
                          cell_border      = NA,
                          frame            = FALSE,
                          width            = NULL,
                          height           = NULL,
                          main             = NULL,
                          show_n           = TRUE,
                          time_label       = "Time",
                          xlab             = NULL,
                          y_label          = NULL,
                          ylab             = NULL,
                          tick             = NULL,
                          ncol             = NULL,
                          nrow             = NULL,
                          combined         = TRUE,
                          legend           = NULL,
                          legend_size      = NULL,
                          legend_title     = NULL,
                          legend_ncol      = NULL,
                          legend_border    = NA,
                          legend_bty       = "n") {

  type <- match.arg(type)
  sort <- match.arg(sort)
  stopifnot(is.logical(combined), length(combined) == 1L)
  if (is.null(legend)) legend <- "right"
  legend <- match.arg(legend, c("bottom", "right", "none"))
  if (!is.null(xlab)) time_label <- xlab

  # Open a new device when width/height supplied (interactive use). In
  # knitr, set fig.width / fig.height in the chunk header instead - this
  # call would open a detached device the chunk won't capture.
  if (interactive() && (!is.null(width) || !is.null(height))) {
    grDevices::dev.new(width  = if (!is.null(width))  width  else 7,
                       height = if (!is.null(height)) height else 7,
                       noRStudioGD = TRUE)
  }

  # Adaptive legend_size: measured against the *actual* device dimensions,
  # so the same call produces a proportional legend on a 7-inch knitr figure
  # and on a 12-inch interactive window.
  legend_size <- .auto_legend_size(legend_size)

  if (type == "distribution") {
    return(distribution_plot(
      x, group = group, scale = match.arg(scale), geom = match.arg(geom),
      na = na, state_colors = state_colors, na_color = na_color,
      frame = frame,
      main = main, show_n = show_n,
      time_label = time_label, y_label = y_label, ylab = ylab,
      tick = tick, ncol = ncol, nrow = nrow, combined = combined,
      legend = legend, legend_size = legend_size,
      legend_title = legend_title, legend_ncol = legend_ncol,
      legend_border = legend_border, legend_bty = legend_bty))
  }

  if (type == "heatmap") {
    return(.sequence_plot_heatmap(
      x, sort, tree, dendrogram_width, k, k_color, k_line_width,
      state_colors, na_color, cell_border, frame,
      main, show_n, time_label, tick,
      legend, legend_size, legend_title, legend_ncol,
      legend_border, legend_bty))
  }

  .sequence_plot_index(
    x, sort, group, row_gap,
    state_colors, na_color, cell_border, frame,
    main, show_n, time_label, tick, ncol, nrow, combined,
    legend, legend_size, legend_title, legend_ncol,
    legend_border, legend_bty)
}


# ---- Heatmap (single-panel dense carpet with optional dendrogram) ---------
.sequence_plot_heatmap <- function(x, sort, tree, dendrogram_width,
                                   k, k_color, k_line_width,
                                   state_colors, na_color, cell_border, frame,
                                   main, show_n, time_label, tick,
                                   legend, legend_size, legend_title,
                                   legend_ncol, legend_border, legend_bty) {

  # Extract data from various input types
  extracted <- .extract_seqplot_input(x)
  clustering <- extracted$clustering
  sort_used <- sort

  # Build dendrogram from clustering distance matrix if available
  if (!is.null(clustering) && is.null(tree) && !is.null(clustering$distance)) {
    hc_method <- if (!is.null(clustering$method) &&
                     clustering$method %in% c("ward.D", "ward.D2", "single",
                                              "complete", "average", "mcquitty",
                                              "median", "centroid")) {
      clustering$method
    } else "ward.D2"
    tree <- stats::hclust(stats::as.dist(clustering$distance), method = hc_method)
    sort_used <- "net_clustering"
  }

  x <- extracted$data

  enc        <- .encode_states(x)
  codes      <- enc$codes
  levels_all <- enc$levels
  n_rows     <- nrow(codes); n_cols <- ncol(codes)

  ord  <- .row_order(codes, sort, tree, levels_all)
  tree <- attr(ord, "tree")

  if (length(ord) != n_rows) {
    stop(sprintf("`tree` has %d leaves but `x` has %d rows.",
                 length(ord), n_rows), call. = FALSE)
  }

  palette <- .state_palette(state_colors, length(levels_all))
  z       <- t(codes[ord, , drop = FALSE])

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  oma <- .legend_oma_size(levels_all, legend, legend_size,
                          legend_ncol, legend_title)
  graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                        0.3, 0.3,
                        if (legend == "right")  oma[["oma_r"]] else 0.3))

  if (!is.null(tree)) {
    graphics::layout(matrix(c(1L, 2L), nrow = 1L),
                     widths = c(dendrogram_width, 4))
  }
  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n)) 2.5 else 1.5

  if (!is.null(tree)) {
    # mar[4] = 0.3: leave a sliver of room for the root branch, otherwise
    # the joining line where all clades merge gets clipped at the panel
    # boundary. `ylim = c(0.5, n_rows + 0.5)` aligns leaves with image
    # cells (which span y = 0.5 .. n_rows + 0.5).
    graphics::par(mar = c(mar_bottom, 0.5, mar_top, 0.3))
    graphics::plot(tree, horiz = TRUE, axes = FALSE, yaxs = "i",
                   ylim = c(0.5, n_rows + 0.5),
                   leaflab = "none", xlab = "", ylab = "")
  }
  graphics::par(mar = c(mar_bottom, if (!is.null(tree)) 0 else 2,
                        mar_top, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0.5, n_cols + 0.5),
                        ylim = c(0.5, n_rows + 0.5), yaxs = "i", xaxs = "i")
  graphics::rect(0.5, 0.5, n_cols + 0.5, n_rows + 0.5,
                 col = na_color, border = NA)
  # useRaster = TRUE: render as bitmap instead of one polygon per cell.
  # Massively smaller output, crisper on hi-DPI displays, fine for our
  # uniform grid. Safe: image() itself decides not to raster if the
  # device can't handle it.
  graphics::image(x = seq_len(n_cols), y = seq_len(n_rows),
                  z = z, col = palette, add = TRUE, useRaster = TRUE)
  if (isTRUE(frame)) graphics::box()
  ticks <- .tick_positions(n_cols, tick, colnames(codes))
  graphics::axis(1, at = ticks$at, labels = ticks$labels,
                 mgp = c(1, 0.5, 0), lwd = if (isTRUE(frame)) 1 else 0,
                 lwd.ticks = 1)
  if (!is.null(time_label) && nzchar(time_label)) {
    graphics::mtext(time_label, side = 1, line = 1.7)
  }
  if (!is.na(cell_border)) {
    graphics::abline(v = seq(0.5, n_cols + 0.5, 1),
                     col = cell_border, lwd = 0.3)
    graphics::abline(h = seq(0.5, n_rows + 0.5, 1),
                     col = cell_border, lwd = 0.3)
    if (isTRUE(frame)) graphics::box()
  }

  # Cluster separators: cut the dendrogram (if any) into k groups, find
  # where group membership changes in the plotted row order, draw thin
  # lines at those row boundaries.
  if (!is.null(k) && !is.null(tree)) {
    stopifnot(is.numeric(k), length(k) == 1L, k >= 2L, k < n_rows)
    hc  <- stats::as.hclust(tree)
    mem <- stats::cutree(hc, k = as.integer(k))[ord]
    boundaries <- which(diff(mem) != 0L)
    if (length(boundaries) > 0L) {
      graphics::abline(h = boundaries + 0.5, col = k_color, lwd = k_line_width)
    }
  }
  t_line <- .title_line(main, show_n, n_rows)
  if (!is.null(t_line)) graphics::mtext(t_line, side = 3, line = 1, font = 2)

  if (legend != "none") {
    .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  invisible(list(ord = ord, codes = codes, palette = palette,
                 levels = levels_all, sort_used = sort_used))
}


# ---- Index (gapped per-row bars, optional facet grid) ---------------------
.sequence_plot_index <- function(x, sort, group, row_gap,
                                 state_colors, na_color, cell_border, frame,
                                 main, show_n, time_label, tick, ncol, nrow,
                                 combined,
                                 legend, legend_size, legend_title,
                                 legend_ncol, legend_border, legend_bty) {

  stopifnot(is.numeric(row_gap), length(row_gap) == 1L,
            row_gap >= 0, row_gap < 1)

  # Extract data and group from various input types
  extracted <- .extract_seqplot_input(x, group)
  x <- extracted$data
  group <- extracted$group

  if (!is.null(group)) {
    stopifnot(length(group) == nrow(x))
    group <- as.factor(group)
  }

  enc        <- .encode_states(x)
  codes      <- enc$codes
  levels_all <- enc$levels
  n_cols     <- ncol(codes)
  palette    <- .state_palette(state_colors, length(levels_all))

  groups       <- if (is.null(group)) factor(rep("all", nrow(codes))) else group
  group_levels <- levels(groups)
  G            <- length(group_levels)
  group_sizes  <- vapply(group_levels,
                         function(g) sum(groups == g), integer(1))

  ords <- lapply(group_levels, function(g) {
    rows <- which(groups == g)
    sub  <- codes[rows, , drop = FALSE]
    rows[.row_order(sub, sort, tree = NULL, levels_all)]
  })

  if (is.null(ncol) && is.null(nrow)) {
    ncol <- as.integer(ceiling(sqrt(G)))
    nrow <- as.integer(ceiling(G / ncol))
  } else if (is.null(ncol)) {
    ncol <- as.integer(ceiling(G / nrow))
  } else if (is.null(nrow)) {
    nrow <- as.integer(ceiling(G / ncol))
  }
  stopifnot(ncol * nrow >= G)

  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  oma <- .legend_oma_size(levels_all, legend, legend_size,
                          legend_ncol, legend_title)
  graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                        0.3, 0.3,
                        if (legend == "right")  oma[["oma_r"]] else 0.3))
  use_layout <- combined && nrow * ncol > 1L
  if (use_layout) {
    lm <- matrix(c(seq_len(G), integer(nrow * ncol - G)),
                 nrow = nrow, ncol = ncol, byrow = TRUE)
    graphics::layout(lm)
  }

  mar_bottom <- if (!is.null(time_label) && nzchar(time_label)) 3 else 1
  mar_top    <- if (!is.null(main) || isTRUE(show_n) || G > 1L) 2.5 else 1.5
  ticks      <- .tick_positions(n_cols, tick, colnames(codes))

  invisible(lapply(seq_len(G), function(g_idx) {
    g    <- group_levels[g_idx]
    ord  <- ords[[g_idx]]
    sub  <- codes[ord, , drop = FALSE]
    nN   <- nrow(sub)

    if (!combined && g_idx > 1L) {
      graphics::par(oma = c(if (legend == "bottom") oma[["oma_b"]] else 0.3,
                            0.3, 0.3,
                            if (legend == "right")  oma[["oma_r"]] else 0.3))
    }
    graphics::par(mar = c(mar_bottom, 2, mar_top, 1))
    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, n_cols + 0.5),
                          ylim = c(0.5, nN + 0.5), yaxs = "i", xaxs = "i")

    # Fill whole panel with na_color so NA cells show through.
    graphics::rect(0.5, 0.5, n_cols + 0.5, nN + 0.5,
                   col = na_color, border = NA)

    # Vectorised per-cell rects with a vertical gap.
    ts <- rep(seq_len(n_cols), each  = nN)
    ss <- rep(seq_len(nN),     times = n_cols)
    cell_col <- palette[as.vector(sub)]
    keep <- !is.na(cell_col)
    graphics::rect(xleft   = ts[keep] - 0.5,
                   ybottom = ss[keep] - 0.5 + row_gap / 2,
                   xright  = ts[keep] + 0.5,
                   ytop    = ss[keep] + 0.5 - row_gap / 2,
                   col     = cell_col[keep],
                   border  = if (is.na(cell_border)) NA else cell_border,
                   lwd     = 0.3)
    if (isTRUE(frame)) graphics::box()
    graphics::axis(1, at = ticks$at, labels = ticks$labels,
                   mgp = c(1, 0.5, 0),
                   lwd = if (isTRUE(frame)) 1 else 0, lwd.ticks = 1)
    if (!is.null(time_label) && nzchar(time_label)) {
      graphics::mtext(time_label, side = 1, line = 1.7)
    }

    panel_title <- if (G > 1L) {
      base <- sprintf("Cluster %s", g)
      if (isTRUE(show_n)) sprintf("%s  (n = %d)", base,
                                  group_sizes[g_idx]) else base
    } else .title_line(main, show_n, nN)
    if (!is.null(panel_title)) {
      graphics::mtext(panel_title, side = 3, line = 0.5, font = 2,
                      cex = if (G > 1L) 0.9 else 1)
    }
    if (!combined && legend != "none") {
      .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                          legend_ncol, legend_title, legend_border, legend_bty)
    }
  }))
  if (combined && G > 1L && !is.null(main)) {
    graphics::mtext(main, side = 3, line = 1.5, font = 2,
                    outer = TRUE, cex = 1)
  }

  if (combined && legend != "none") {
    .draw_legend_in_oma(levels_all, palette, legend, legend_size,
                        legend_ncol, legend_title, legend_border, legend_bty)
  }

  invisible(list(codes = codes, palette = palette, levels = levels_all,
                 orders = ords, groups = group_levels))
}


# ---- internal: row ordering for a codes matrix ----------------------------
.row_order <- function(codes, sort, tree, levels_all) {
  n_rows <- nrow(codes); n_cols <- ncol(codes)
  if (!is.null(tree)) {
    if (!inherits(tree, "dendrogram")) tree <- stats::as.dendrogram(tree)
    return(structure(stats::order.dendrogram(tree), tree = tree))
  }
  if (sort == "start") {
    ord <- do.call(order, c(lapply(seq_len(n_cols), function(t) codes[, t]),
                            list(na.last = TRUE)))
    return(structure(ord, tree = NULL))
  }
  if (sort == "end") {
    ord <- do.call(order, c(lapply(rev(seq_len(n_cols)),
                                   function(t) codes[, t]),
                            list(na.last = TRUE)))
    return(structure(ord, tree = NULL))
  }
  if (sort == "frequency") {
    counts <- vapply(seq_along(levels_all),
                     function(k) rowSums(codes == k, na.rm = TRUE),
                     numeric(n_rows))
    tr <- stats::as.dendrogram(stats::hclust(stats::dist(counts),
                                             method = "ward.D2"))
    return(structure(stats::order.dendrogram(tr), tree = tr))
  }
  # distance-based
  d  <- .sequence_distance(codes_to_char(codes, levels_all),
                           dissimilarity = sort)
  tr <- stats::as.dendrogram(stats::hclust(stats::as.dist(d),
                                           method = "ward.D2"))
  structure(stats::order.dendrogram(tr), tree = tr)
}


# ---- internal: decode integer codes back to character (for .sequence_distance)
codes_to_char <- function(codes, levels_all) {
  out <- matrix(NA_character_, nrow = nrow(codes), ncol = ncol(codes))
  nz  <- !is.na(codes)
  out[nz] <- levels_all[codes[nz]]
  out
}


# ---- internal: standard title with "(n = N)" suffix ----------------------
.title_line <- function(main, show_n, n) {
  if (!is.null(main)) {
    if (isTRUE(show_n)) sprintf("%s  (n = %d)", main, n) else main
  } else if (isTRUE(show_n)) sprintf("n = %d", n) else NULL
}
