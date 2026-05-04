# chain_structure() -- qualitative analysis of a discrete-time Markov chain.
#
# Computes properties that are independent of any starting distribution:
# state classification (absorbing / recurrent / transient), communicating
# classes (SCCs of the support graph), per-state period, hitting probability
# matrix H[i,j] = P(ever reach j | start at i), and absorption analysis when
# absorbing states exist. Used as a diagnostic before trusting passage_time()
# / markov_stability(): if a chain is not regular, the stationary distribution
# those report can mix multiple behavioural phases.
#
# All linear algebra; no new dependencies.

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Reachability matrix of a directed adjacency matrix.
#'
#' `R[i, j] = 1` iff there is a path of any positive length from i to j (note:
#' a self-loop or a cycle through i marks `R[i, i] = 1`; a state with no
#' outgoing edges has `R[i, ]` all zero except itself which is also zero).
#' Computed via repeated boolean multiplication; O(n^4) but n is small for
#' Markov chains in this package.
#' @noRd
.cs_reachability <- function(A) {
  n <- nrow(A)
  R <- A > 0
  for (k in seq_len(n - 1L)) {
    R_new <- (R | (R %*% A) > 0)
    if (identical(R_new, R)) break
    R <- R_new
  }
  R * 1L
}

#' Strongly connected components (communicating classes) of A.
#'
#' Two states are in the same SCC iff each is reachable from the other.
#' Returns a list of integer index vectors. Self-loop is counted as
#' "reachable from self".
#' @noRd
.cs_scc <- function(A) {
  n <- nrow(A)
  if (n == 0L) return(list())
  R <- .cs_reachability(A)
  diag(R) <- 1L
  M <- R & t(R)
  keys <- apply(M, 1L, function(row) paste(which(row), collapse = ","))
  unique_keys <- unique(keys)
  lapply(unique_keys, function(k) which(keys == k))
}

#' GCD of a non-empty integer vector via Euclid (no external deps).
#' @noRd
.cs_gcd <- function(x) {
  x <- x[x > 0L]
  if (length(x) == 0L) return(NA_integer_)
  out <- x[1L]
  for (v in x[-1L]) {
    while (v != 0L) { tmp <- out %% v; out <- v; v <- tmp }
  }
  as.integer(out)
}

#' Period of a recurrent class given its support submatrix.
#'
#' Period of state i = gcd{n >= 1 : P^n[i, i] > 0}. For an irreducible
#' aperiodic class on a self-looped state, returns 1. Returns NA if no
#' return is possible (which shouldn't happen for a recurrent class).
#' @noRd
.cs_class_period <- function(A_sub) {
  n <- nrow(A_sub)
  if (n == 0L) return(NA_integer_)
  if (n == 1L) return(if (A_sub[1L, 1L] > 0) 1L else NA_integer_)
  An <- A_sub > 0
  An_running <- An
  cycle_lengths <- integer()
  if (any(diag(An_running))) cycle_lengths <- c(cycle_lengths, 1L)
  for (k in seq_len(2L * n)) {
    An_running <- ((An_running %*% An) > 0)
    if (any(diag(An_running))) {
      cycle_lengths <- c(cycle_lengths, k + 1L)
      if (length(cycle_lengths) >= 3L && .cs_gcd(cycle_lengths) == 1L) break
    }
  }
  .cs_gcd(cycle_lengths)
}

#' Hitting probability matrix `H[i, j] = P(T_j < infty | X_0 = i)`.
#'
#' For i != j: T_j = inf{n >= 0 : X_n = j}, so the off-diagonal hitting
#' probability is the standard "eventually reach j" probability. For i == j:
#' uses T_j = inf{n >= 1 : X_n = j} (return-time convention, matching
#' `markovchain::hittingProbabilities`), so `H[j, j] = sum_k P[j, k] H[k, j]`.
#'
#' For each column j, the off-diagonal system is solved over the states from
#' which j is reachable; states that cannot reach j get `H[i, j] = 0` (the
#' standard *minimal non-negative solution*; cf. Norris, *Markov Chains*,
#' Theorem 1.3.2). Without the reachability restriction the system is
#' rank-deficient whenever a closed class disjoint from j exists.
#' @noRd
.cs_hitting <- function(P, reach = NULL) {
  n <- nrow(P)
  state_names <- rownames(P)
  H <- matrix(0, n, n, dimnames = list(state_names, state_names))
  if (is.null(reach)) {
    A <- (P > 0) * 1L
    reach <- .cs_reachability(A)
    diag(reach) <- 1L
  }
  for (j in seq_len(n)) {
    if (n == 1L) {
      H[j, j] <- if (P[1L, 1L] > 0) 1 else 0
      next
    }
    can_reach <- setdiff(which(reach[, j] > 0), j)
    if (length(can_reach) > 0L) {
      P_sub <- P[can_reach, can_reach, drop = FALSE]
      rhs   <- P[can_reach, j]
      I_sub <- diag(length(can_reach))
      h <- tryCatch(
        drop(solve(I_sub - P_sub, rhs)),
        error = function(e) rep(NA_real_, length(can_reach))
      )
      H[can_reach, j] <- h
    }
  }
  # Diagonal via return-time recursion. H[j, j] is the probability of
  # returning to j in >= 1 steps: the self-loop contributes P[j, j],
  # and any path leaving j must eventually hit j again with probability
  # H[k, j] (already computed off-diagonal).
  for (j in seq_len(n)) {
    off_idx <- setdiff(seq_len(n), j)
    H[j, j] <- P[j, j] + sum(P[j, off_idx] * H[off_idx, j])
  }
  H[!is.na(H)] <- pmin(pmax(H[!is.na(H)], 0), 1)
  H
}

#' Absorption analysis (transient -> absorbing) via fundamental matrix.
#'
#' Returns NULL when the canonical transient/absorbing partition is empty
#' (no absorbing states or no transient states feeding into them).
#' @noRd
.cs_absorption <- function(P, transient_idx, absorbing_idx, state_names) {
  if (length(absorbing_idx) == 0L || length(transient_idx) == 0L) {
    return(list(probabilities = NULL, mean_time = NULL))
  }
  Q <- P[transient_idx, transient_idx, drop = FALSE]
  R <- P[transient_idx, absorbing_idx,  drop = FALSE]
  I <- diag(nrow(Q))
  N <- tryCatch(solve(I - Q),
                error = function(e) NULL)
  if (is.null(N)) return(list(probabilities = NULL, mean_time = NULL))
  abs_probs <- N %*% R
  rownames(abs_probs) <- state_names[transient_idx]
  colnames(abs_probs) <- state_names[absorbing_idx]
  abs_time <- as.numeric(N %*% rep(1, nrow(N)))
  names(abs_time) <- state_names[transient_idx]
  list(probabilities = abs_probs, mean_time = abs_time)
}

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Qualitative structure of a discrete-time Markov chain.
#'
#' Computes properties that depend only on the transition matrix support,
#' not on any starting distribution: state classification, communicating
#' classes, periods, irreducibility / aperiodicity / regularity /
#' reversibility, hitting probabilities, and absorption analysis when
#' absorbing states exist.
#'
#' @param x A `netobject`, `cograph_network`, `tna` model, transition
#'   matrix, or sequence data.frame (passed through `build_network()` with
#'   `method = "relative"`).
#' @param normalize Logical. If `TRUE` (default), rows of the transition
#'   matrix are renormalized to sum to 1 before analysis (see
#'   [passage_time()] for the same convention).
#' @param tol Numerical tolerance for the reversibility check (detailed
#'   balance) and for treating near-zero entries as zero when building the
#'   support graph. Default `1e-10`.
#' @return A `chain_structure` object: a list with elements
#'   \describe{
#'     \item{`states`}{Character vector of state names.}
#'     \item{`classification`}{Named character vector. One of `"absorbing"`,
#'       `"recurrent"`, `"transient"` per state.}
#'     \item{`communicating_classes`}{List of state-name vectors. Each
#'       sublist is a strongly connected component of the support graph.}
#'     \item{`recurrent_classes`}{Subset of `communicating_classes` that are
#'       closed (no transitions leaving the class).}
#'     \item{`transient_classes`}{Subset that are not closed.}
#'     \item{`absorbing_states`}{Character vector of states with `P[i, i] = 1`.}
#'     \item{`period`}{Named integer vector. Period of each recurrent state;
#'       `NA` for transient states.}
#'     \item{`is_irreducible`}{Logical. `TRUE` iff there is exactly one
#'       communicating class.}
#'     \item{`is_aperiodic`}{Logical. `TRUE` iff every recurrent state has
#'       period 1.}
#'     \item{`is_regular`}{Logical. `is_irreducible && is_aperiodic`.}
#'     \item{`is_reversible`}{Logical or `NA`. `TRUE` iff the chain
#'       satisfies detailed balance against its stationary distribution.
#'       `NA` for non-irreducible chains (no unique stationary).}
#'     \item{`hitting_probabilities`}{`n x n` matrix. `[i, j] = P(ever reach j
#'       starting from i)`.}
#'     \item{`absorption_probabilities`}{`n_transient x n_absorbing` matrix
#'       or `NULL` if no transient -> absorbing pathway exists. `[i, j] =
#'       P(eventual absorption in j | start in i)`.}
#'     \item{`mean_absorption_time`}{Named numeric vector or `NULL`. Expected
#'       number of steps until absorption from each transient state.}
#'     \item{`P`}{The (possibly normalized) transition matrix used.}
#'   }
#' @details
#'   Built specifically as a diagnostic to run *before* trusting the output
#'   of [passage_time()] or [markov_stability()]. Both implicitly assume a
#'   regular chain (irreducible + aperiodic) so that the stationary
#'   distribution is unique and meaningful. Use `is_regular` to check.
#'
#'   The fundamental-matrix absorption math follows Kemeny & Snell (1976);
#'   the hitting-probability linear system follows Norris (1997).
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' cs  <- chain_structure(net)
#' print(cs)
#' \donttest{
#' summary(cs)
#' }
#'
#' @seealso [passage_time()], [markov_stability()], [build_network()]
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' Norris, J. R. (1997). \emph{Markov Chains}. Cambridge University Press.
#'
#' @export
chain_structure <- function(x, normalize = TRUE, tol = 1e-10) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      chain_structure(net, normalize = normalize, tol = tol)
    })
    class(out) <- c("chain_structure_group", "list")
    return(out)
  }
  P <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }
  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)
  n <- nrow(P)

  # ---- Support graph and SCCs ----
  A <- (P > tol) * 1L
  classes <- .cs_scc(A)

  # A class is closed (recurrent) iff no outgoing edges leave it.
  is_closed <- vapply(classes, function(cl) {
    if (n == length(cl)) return(TRUE)
    out_rows <- A[cl, , drop = FALSE]
    !any(out_rows[, setdiff(seq_len(n), cl), drop = FALSE] > 0)
  }, logical(1))
  recurrent_classes <- classes[is_closed]
  transient_classes <- classes[!is_closed]

  recurrent_idx <- if (length(recurrent_classes) > 0L)
    sort(unlist(recurrent_classes)) else integer(0)
  transient_idx <- if (length(transient_classes) > 0L)
    sort(unlist(transient_classes)) else integer(0)
  absorbing_mask <- vapply(seq_len(n), function(i)
    abs(P[i, i] - 1) < tol, logical(1))
  absorbing_idx <- which(absorbing_mask)

  classification <- character(n)
  classification[recurrent_idx] <- "recurrent"
  classification[transient_idx] <- "transient"
  classification[absorbing_idx] <- "absorbing"
  names(classification) <- state_names

  # ---- Period per state (only meaningful for recurrent states) ----
  period <- rep(NA_integer_, n)
  for (cl in recurrent_classes) {
    cl_period <- .cs_class_period(A[cl, cl, drop = FALSE])
    period[cl] <- cl_period
  }
  names(period) <- state_names

  # ---- Chain-level properties ----
  is_irreducible <- length(classes) == 1L
  is_aperiodic <- length(recurrent_classes) > 0L &&
                  all(period[recurrent_idx] == 1L, na.rm = TRUE)
  is_regular <- is_irreducible && is_aperiodic

  is_reversible <- if (is_irreducible) {
    pi <- tryCatch(.mpt_stationary(P), error = function(e) NULL)
    if (is.null(pi)) NA else {
      flux <- pi * P
      max_dev <- max(abs(flux - t(flux)))
      max_dev < tol * max(abs(flux), 1)
    }
  } else NA

  # ---- Hitting probabilities and absorption analysis ----
  H <- .cs_hitting(P)
  abs_info <- .cs_absorption(P, transient_idx, absorbing_idx, state_names)

  structure(list(
    states                   = state_names,
    classification           = classification,
    communicating_classes    = lapply(classes, function(cl) state_names[cl]),
    recurrent_classes        = lapply(recurrent_classes, function(cl) state_names[cl]),
    transient_classes        = lapply(transient_classes, function(cl) state_names[cl]),
    absorbing_states         = state_names[absorbing_idx],
    period                   = period,
    is_irreducible           = is_irreducible,
    is_aperiodic             = is_aperiodic,
    is_regular               = is_regular,
    is_reversible            = is_reversible,
    hitting_probabilities    = H,
    absorption_probabilities = abs_info$probabilities,
    mean_absorption_time     = abs_info$mean_time,
    P                        = P
  ), class = "chain_structure")
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print method for `chain_structure`.
#'
#' Prints a compact chain-level header. For the full per-state table,
#' call `summary()` on the same object.
#' @param x A `chain_structure` object.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.chain_structure <- function(x, ...) {
  rev <- x$is_reversible
  cat(sprintf("Chain structure  [%d states, %d communicating classes]\n",
              length(x$states), length(x$communicating_classes)))
  cat(sprintf(
    "  irreducible: %s   aperiodic: %s   regular: %s   reversible: %s\n",
    x$is_irreducible, x$is_aperiodic, x$is_regular,
    if (is.na(rev)) "NA" else rev))
  cat(sprintf("  recurrent classes: %d   transient classes: %d\n",
              length(x$recurrent_classes), length(x$transient_classes)))
  if (length(x$absorbing_states) > 0L) {
    cat(sprintf("  absorbing states: %s\n",
                paste(x$absorbing_states, collapse = ", ")))
  }
  cat("\nUse summary(x) for the per-state table, plot(x) for the heatmap.\n")
  invisible(x)
}

#' Plot method for `chain_structure`.
#'
#' Renders the hitting-probability matrix as a heatmap, with rows and
#' columns ordered by communicating class so the block structure is
#' visible at a glance. State labels along both axes are coloured by
#' classification (absorbing / recurrent / transient). The subtitle
#' summarises the chain-level properties (regular, reversible).
#'
#' Cell colour encodes `P(ever reach j | start at i)`. The diagonal
#' uses the return-time convention (`P(return to j in >= 1 steps)`),
#' matching `markovchain::hittingProbabilities`. A non-irreducible chain
#' shows zero off-block entries -- visual evidence of one-way doors
#' between behavioural phases. An absorbing chain shows a column of 1's
#' for the absorbing state.
#'
#' @param x A `chain_structure` object.
#' @param show_values Logical. If `TRUE` (default), prints the numeric
#'   probability inside each cell. Set `FALSE` for large state spaces
#'   (n > 10) where labels overlap.
#' @param digits Integer. Decimal places for in-cell labels.
#' @param ... Ignored.
#' @return A `ggplot` object.
#' @export
plot.chain_structure <- function(x, show_values = TRUE, digits = 2L, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting.", call. = FALSE)
  }

  # Order states: members of each communicating class kept contiguous,
  # classes ordered by [recurrent first, then transient], and within each
  # group classes ordered by size (descending) for stable layout.
  classes <- x$communicating_classes
  is_rec  <- vapply(classes, function(cl) {
    any(x$classification[cl] %in% c("recurrent", "absorbing"))
  }, logical(1))
  rec_classes <- classes[is_rec]
  tra_classes <- classes[!is_rec]
  rec_classes <- rec_classes[order(-vapply(rec_classes, length, integer(1)))]
  tra_classes <- tra_classes[order(-vapply(tra_classes, length, integer(1)))]
  ordered_states <- unlist(c(rec_classes, tra_classes), use.names = FALSE)
  if (length(ordered_states) != length(x$states)) {
    ordered_states <- x$states  # fallback
  }

  H <- x$hitting_probabilities[ordered_states, ordered_states, drop = FALSE]
  cls <- x$classification[ordered_states]

  # Long form for ggplot
  df <- expand.grid(from = ordered_states, to = ordered_states,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$prob <- as.numeric(H)
  df$from <- factor(df$from, levels = ordered_states)
  df$to   <- factor(df$to,   levels = ordered_states)

  # Okabe-Ito-flavoured classification palette.
  pal <- c(absorbing  = "#D55E00",   # vermillion
           recurrent  = "#0072B2",   # blue
           transient  = "#E69F00")   # orange
  axis_cols <- pal[cls]

  subtitle <- sprintf(
    "irreducible: %s   aperiodic: %s   regular: %s   reversible: %s   classes: %d",
    x$is_irreducible, x$is_aperiodic, x$is_regular,
    if (is.na(x$is_reversible)) "NA" else x$is_reversible,
    length(x$communicating_classes)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$to, y = .data$from,
                                         fill = .data$prob)) +
    ggplot2::geom_tile(colour = "grey85", linewidth = 0.4) +
    ggplot2::scale_fill_gradient(low = "white", high = "#08306B",
                                  limits = c(0, 1),
                                  name = "P(reach)",
                                  na.value = "grey90") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::scale_y_discrete(limits = rev(ordered_states)) +
    ggplot2::labs(
      title = "Hitting probabilities",
      subtitle = subtitle,
      x = "to (target state)",
      y = "from (source state)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x.top  = ggplot2::element_text(angle = 45, hjust = 0,
                                                vjust = 0,
                                                colour = axis_cols),
      axis.text.y      = ggplot2::element_text(colour = rev(axis_cols)),
      plot.subtitle    = ggplot2::element_text(size = 10, colour = "grey30")
    )

  if (show_values && nrow(df) <= 400L) {
    df$label <- ifelse(df$prob > 0,
                       formatC(df$prob, digits = digits, format = "f"),
                       "")
    p <- p + ggplot2::geom_text(
      data = df,
      ggplot2::aes(label = .data$label,
                    colour = .data$prob > 0.5),
      size = 3, show.legend = FALSE
    ) +
      ggplot2::scale_colour_manual(
        values = c(`TRUE` = "white", `FALSE` = "black")
      )
  }

  p
}

#' Tidy per-state summary of a `chain_structure`.
#'
#' Returns a single data.frame with one row per state, combining every
#' per-state metric `chain_structure()` computes. Always includes
#' `state`, `classification`, `period`, `return_probability` (the
#' diagonal of the hitting matrix) and `persistence` (the diagonal of
#' the transition matrix); adds `sojourn` whenever it is finite, the
#' chain's `stationary_probability` when irreducible, and absorption
#' columns when the chain has any absorbing states.
#'
#' Columns are ordered for readability: identifiers first, classification
#' second, dynamic per-state metrics last.
#'
#' @param object A `chain_structure` object.
#' @param ... Ignored.
#' @return A `data.frame` with one row per state. Columns described above.
#' @export
summary.chain_structure <- function(object, ...) {
  state_names <- object$states
  P <- object$P
  H <- object$hitting_probabilities

  persistence <- diag(P)
  sojourn <- ifelse(persistence < 1 - .Machine$double.eps,
                    1 / (1 - persistence), Inf)

  out <- data.frame(
    state              = state_names,
    classification     = unname(object$classification),
    period             = unname(object$period),
    persistence        = round(unname(persistence), 4),
    return_probability = round(unname(diag(H)), 4),
    sojourn_steps      = round(unname(sojourn), 2),
    stringsAsFactors   = FALSE
  )

  if (object$is_irreducible) {
    pi <- tryCatch(.mpt_stationary(P), error = function(e) NULL)
    if (!is.null(pi)) {
      out$stationary_probability <- round(unname(pi), 4)
    }
  }

  if (!is.null(object$absorption_probabilities)) {
    abs_states <- colnames(object$absorption_probabilities)
    abs_probs <- object$absorption_probabilities
    for (a in abs_states) {
      colname <- if (length(abs_states) == 1L) "absorption_probability"
                 else sprintf("absorbed_in_%s", a)
      vals <- rep(NA_real_, length(state_names))
      vals[match(rownames(abs_probs), state_names)] <- round(abs_probs[, a], 4)
      out[[colname]] <- vals
    }
    abs_time <- rep(NA_real_, length(state_names))
    abs_time[match(names(object$mean_absorption_time), state_names)] <-
      round(object$mean_absorption_time, 2)
    out$mean_absorption_time <- abs_time
  }

  rownames(out) <- NULL
  class(out) <- c("summary_chain_structure", "data.frame")
  attr(out, "is_regular") <- object$is_regular
  attr(out, "is_irreducible") <- object$is_irreducible
  attr(out, "is_aperiodic") <- object$is_aperiodic
  attr(out, "is_reversible") <- object$is_reversible
  attr(out, "n_classes") <- length(object$communicating_classes)
  attr(out, "absorbing_states") <- object$absorbing_states
  out
}

#' Print method for `chain_structure_group`.
#'
#' One header line per group, followed by each group's per-state
#' table (via `summary.chain_structure`).
#' @param x A `chain_structure_group` (named list of `chain_structure`).
#' @param ... Forwarded to `print.chain_structure`.
#' @return `x` invisibly.
#' @export
print.chain_structure_group <- function(x, ...) {
  cat(sprintf("Chain structure -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}

#' Cross-group comparison of `chain_structure_group`.
#'
#' Produces a single tidy data.frame with one row per (group, state)
#' combination, combining classification, persistence, sojourn, and --
#' when applicable -- stationary or mean-absorption-time columns. Useful
#' for side-by-side reporting of `chain_structure()` across the
#' members of a `netobject_group`.
#'
#' @param object A `chain_structure_group`.
#' @param ... Ignored.
#' @return A `data.frame` with columns `group`, `state`, `classification`,
#'   `period`, `persistence`, `return_probability`, `sojourn_steps`, plus
#'   `stationary_probability` if all groups are irreducible and
#'   `mean_absorption_time` if any group has absorbing states.
#' @export
summary.chain_structure_group <- function(object, ...) {
  parts <- lapply(names(object), function(nm) {
    s <- summary(object[[nm]])
    s$group <- nm
    s
  })
  # Union of columns across groups, preserving order
  all_cols <- unique(unlist(lapply(parts, colnames)))
  parts <- lapply(parts, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  out <- do.call(rbind, parts)
  out <- out[, c("group", setdiff(all_cols, "group")), drop = FALSE]
  rownames(out) <- NULL
  # The merged frame is a plain comparison table -- drop the
  # per-chain `summary_chain_structure` class (and its attributes) so
  # `print.data.frame` is dispatched, not the per-chain pretty-printer.
  class(out) <- "data.frame"
  out
}

#' Print method for `summary.chain_structure`.
#'
#' Prints a one-line chain header followed by the tidy per-state table.
#' @param x A `summary_chain_structure` object.
#' @param ... Forwarded to `print.data.frame`.
#' @return `x` invisibly.
#' @export
print.summary_chain_structure <- function(x, ...) {
  rev <- attr(x, "is_reversible")
  cat(sprintf(
    "Chain structure summary  [%d states, %d classes]\n",
    nrow(x), attr(x, "n_classes")))
  cat(sprintf(
    "  irreducible: %s   aperiodic: %s   regular: %s   reversible: %s\n",
    attr(x, "is_irreducible"), attr(x, "is_aperiodic"),
    attr(x, "is_regular"), if (is.na(rev)) "NA" else rev))
  abs_states <- attr(x, "absorbing_states")
  if (length(abs_states) > 0L) {
    cat(sprintf("  absorbing states: %s\n",
                paste(abs_states, collapse = ", ")))
  }
  cat("\n")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}
