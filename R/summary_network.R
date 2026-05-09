# ---- Network-level descriptive summary ----
#
# Same metric set as `tna::summary.tna()` but computed from `$weights`
# directly so we do not depend on igraph at runtime. Distances reuse
# `.floyd_warshall_sp()` from R/centrality_measures.R.

#' Network metrics for a netobject
#'
#' Computes node count, edge count, density, mean shortest-path distance,
#' mean and SD of in/out strength, mean and SD of in/out degree, in/out
#' degree centralization (Freeman), and reciprocity. Mirrors the metric set
#' returned by `tna::summary.tna()` so a Nestimate netobject and the
#' equivalent tna model report numerically identical descriptive metrics.
#'
#' @param object A `netobject` (or `cograph_network`) object.
#' @param ... Ignored.
#' @return A `data.frame` with columns `metric` and `value`, of class
#'   `c("summary.netobject", "data.frame")`.
#' @export
summary.netobject <- function(object, ...) {
  stopifnot(inherits(object, "netobject") || inherits(object, "cograph_network"))
  W <- object$weights
  stopifnot(is.matrix(W), nrow(W) == ncol(W))
  directed <- isTRUE(object$directed)

  .summary_metrics_from_weights(W, directed)
}

#' Network metrics for a netobject_group
#'
#' Returns one summary per constituent network. With `combined = TRUE`
#' (default) the per-group tables are joined into a single wide
#' `data.frame` with one column per group; with `combined = FALSE`
#' returns a named list.
#'
#' @param object A `netobject_group`.
#' @param combined Logical. Combine into one wide data.frame? Default `TRUE`.
#' @param ... Ignored.
#' @return Either a `data.frame` (one column per group) or a named list of
#'   `summary.netobject` objects, of class `c("summary.netobject_group", ...)`.
#' @export
summary.netobject_group <- function(object, combined = TRUE, ...) {
  stopifnot(inherits(object, "netobject_group"))
  stopifnot(is.logical(combined), length(combined) == 1L)

  per_group <- lapply(object, summary)

  if (!combined) {
    return(structure(per_group, class = c("summary.netobject_group", "list"),
                     combined = FALSE))
  }

  metric <- per_group[[1L]]$metric
  values <- vapply(per_group, function(s) s$value,
                   numeric(length(metric)))
  out <- data.frame(metric = metric, values, check.names = FALSE,
                    stringsAsFactors = FALSE)
  names(out)[-1L] <- names(per_group)
  structure(out, class = c("summary.netobject_group", "data.frame"),
            combined = TRUE)
}


# ---- Internal: compute the 13 metrics from a weight matrix ----

.summary_metrics_from_weights <- function(W, directed) {
  # Diagonal is preserved as estimated. tna's `summary.tna()` counts edges
  # via `sum(weights > 0)` without diagonal exclusion, computes density via
  # igraph's `edge_density` (which uses the same edge total), and labels
  # in/out strength following the assignment used in tna (`in_strength <-
  # igraph::strength(g, mode = "out")` and the symmetric swap), so the
  # values reported under "In-Strength" come from row sums and the values
  # reported under "Out-Strength" come from column sums.
  n <- nrow(W)
  pos <- W > 0
  edge_count <- sum(pos)

  density <- if (n > 1L) {
    if (directed) edge_count / (n * (n - 1))
    else edge_count / (n * (n - 1) / 2)
  } else 0
  density <- min(density, 1)

  D <- .floyd_warshall_sp(W, invert = FALSE)$D
  finite_off <- is.finite(D) & D > 0
  mean_distance <- if (any(finite_off)) mean(D[finite_off]) else NaN

  abs_W <- abs(W)
  row_strength <- rowSums(abs_W)
  col_strength <- colSums(abs_W)
  out_degree   <- rowSums(W != 0)
  in_degree    <- colSums(W != 0)

  # Centralization uses degrees without self-loops (matches
  # `igraph::centr_degree(..., loops = FALSE)`).
  W_no_loops <- W
  diag(W_no_loops) <- 0
  out_degree_nl <- rowSums(W_no_loops != 0)
  in_degree_nl  <- colSums(W_no_loops != 0)

  # Mirror tna's labelling: row sums are reported as In-Strength,
  # column sums as Out-Strength.
  out_strength <- col_strength
  in_strength  <- row_strength

  cent_out <- .freeman_degree_centralization(out_degree_nl, n, directed)
  cent_in  <- .freeman_degree_centralization(in_degree_nl,  n, directed)

  # igraph's default `reciprocity()` excludes self-loops, so we measure
  # mutual-edge fraction over off-diagonal positions only.
  pos_off <- pos
  diag(pos_off) <- FALSE
  reciprocity <- if (directed && any(pos_off)) {
    mutual <- pos_off & t(pos_off)
    sum(mutual) / sum(pos_off)
  } else if (!directed) {
    1
  } else {
    NaN
  }

  metric <- c(
    "Node Count",
    "Edge Count",
    "Network Density",
    "Mean Distance",
    "Mean Out-Strength",
    "SD Out-Strength",
    "Mean In-Strength",
    "SD In-Strength",
    "Mean Out-Degree",
    "SD Out-Degree",
    "Centralization (Out-Degree)",
    "Centralization (In-Degree)",
    "Reciprocity"
  )
  value <- c(
    n,
    edge_count,
    density,
    mean_distance,
    mean(out_strength, na.rm = TRUE),
    stats::sd(out_strength, na.rm = TRUE),
    mean(in_strength, na.rm = TRUE),
    stats::sd(in_strength, na.rm = TRUE),
    mean(out_degree),
    stats::sd(out_degree),
    cent_out,
    cent_in,
    reciprocity
  )
  out <- data.frame(metric = metric, value = unname(value),
                    stringsAsFactors = FALSE)
  structure(out, class = c("summary.netobject", "data.frame"))
}

# Freeman degree centralization. For directed graphs the maximum possible
# sum of (max - d) across n nodes is (n - 1)^2 — the star graph attains it.
# For undirected graphs the corresponding bound is (n - 1)(n - 2).
.freeman_degree_centralization <- function(d, n, directed) {
  if (n < 2L) return(0)
  num <- sum(max(d) - d)
  denom <- if (directed) (n - 1)^2 else (n - 1) * (n - 2)
  if (denom <= 0) return(0)
  num / denom
}


# ---- print methods ----

#' @export
print.summary.netobject <- function(x, ...) {
  cat("Network metrics:\n")
  out <- x
  out$value <- .format_metric_values(out$metric, out$value)
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}

# Counts print as integers; other metrics as 4 significant digits, formatted
# value-by-value so `format()` cannot promote the whole column to scientific
# notation when one cell is near machine epsilon.
.format_metric_values <- function(metric, value) {
  is_count <- metric %in% c("Node Count", "Edge Count")
  vapply(seq_along(value), function(i) {
    if (is_count[i]) formatC(value[i], format = "d", big.mark = "")
    else formatC(value[i], digits = 4L, format = "g")
  }, character(1L))
}

#' @export
print.summary.netobject_group <- function(x, ...) {
  if (isFALSE(attr(x, "combined"))) {
    for (nm in names(x)) {
      cat("Group:", nm, "\n")
      print(x[[nm]])
      cat("\n")
    }
    return(invisible(x))
  }
  cat("Network metrics by group:\n")
  out <- x
  numeric_cols <- vapply(out, is.numeric, logical(1L))
  out[numeric_cols] <- lapply(out[numeric_cols],
                              function(v) .format_metric_values(out$metric, v))
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}
