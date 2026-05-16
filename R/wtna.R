#' Window-based Transition Network Analysis
#'
#' @description
#' Computes networks from one-hot (binary indicator) data using temporal
#' windowing. Supports transition (directed), co-occurrence (undirected),
#' or both network types.
#'
#' @param data Data frame with one-hot encoded columns (0/1 binary).
#' @param method Character. Network type: \code{"transition"} (directed),
#'   \code{"cooccurrence"} (undirected), or \code{"both"} (returns list of
#'   two networks). Default: \code{"transition"}.
#' @param type Character. Output type: \code{"frequency"} (raw counts) or
#'   \code{"relative"} (row-normalized probabilities). Default: \code{"frequency"}.
#'   Note that \code{type = "relative"} applied to \code{method = "cooccurrence"}
#'   produces an asymmetric matrix (conditional co-occurrence given row state),
#'   not a symmetric undirected weight matrix — use \code{type = "frequency"}
#'   if symmetric co-occurrence counts are required.
#' @param codes Character vector or NULL. Names of the one-hot columns to use.
#'   If NULL, auto-detects binary columns. Default: NULL.
#' @param window_size Integer. Number of consecutive rows to aggregate per
#'   window. Default: 1 (no windowing).
#' @param mode Character. Window mode: \code{"non-overlapping"} (fixed, separate
#'   windows) or \code{"overlapping"} (rolling, step = 1).
#'   Default: \code{"non-overlapping"}.
#' @param actor Character or NULL. Name of the actor/ID column for per-group
#'   computation. If NULL, treats all rows as one group. Default: NULL.
#'
#' @return For \code{method = "transition"} or \code{"cooccurrence"}: a
#'   \code{netobject} (see \code{\link{build_network}}).
#'
#'   For \code{method = "both"}: a \code{wtna_mixed} object with elements
#'   \code{$transition} and \code{$cooccurrence}, each a \code{netobject}.
#'
#' @details
#' \strong{Transitions}: Uses \code{crossprod(X[-n,], X[-1,])} to count
#' how often state i is active at time t AND state j at time t+1.
#'
#' \strong{Co-occurrence}: Uses \code{crossprod(X)} to count states that are
#' simultaneously active in the same row.
#'
#' \strong{Windowing}: For \code{window_size > 1}, rows are aggregated into
#' windows before computing networks. Non-overlapping windows are fixed,
#' separate blocks; overlapping windows roll forward one row at a time.
#' Within each window, any active indicator (1) in any row makes that state
#' active for the window.
#'
#' \strong{Per-actor}: When \code{actor} is specified, networks are computed
#' per group and summed.
#'
#' @examples
#' oh <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,0,0), nrow = 4, byrow = TRUE,
#'              dimnames = list(NULL, c("A","B","C")))
#' w <- wtna(oh)
#'
#' \donttest{
#' # Simple one-hot data
#' df <- data.frame(
#'   A = c(1, 0, 1, 0, 1),
#'   B = c(0, 1, 0, 1, 0),
#'   C = c(0, 0, 1, 0, 0)
#' )
#'
#' # Transition network
#' net <- wtna(df)
#' print(net)
#'
#' # Both networks
#' nets <- wtna(df, method = "both")
#' print(nets$transition)
#' print(nets$cooccurrence)
#'
#' # With windowing
#' net <- wtna(df, window_size = 2, mode = "non-overlapping")
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{prepare_onehot}}
#'
#' @export
wtna <- function(data,
                 method = c("transition", "cooccurrence", "both"),
                 type = c("frequency", "relative"),
                 codes = NULL,
                 window_size = 3L,
                 mode = c("non-overlapping", "overlapping"),
                 actor = NULL) {
  method <- match.arg(method)
  type <- match.arg(type)
  mode <- match.arg(mode)

  df <- as.data.frame(data)
  codes <- .resolve_codes(df, codes, exclude = actor)

  stopifnot(length(codes) >= 2L)

  if (is.null(actor)) {
    X_raw <- .wtna_code_matrix(df, codes)
    weights <- .wtna_compute_weights(X_raw, method, window_size, mode)
  } else {
    stopifnot(all(actor %in% names(df)))
    weights <- .wtna_compute_by_actor(df, codes, window_size, mode, actor, method)
  }

  # Initial state probs for transition networks (directed only)
  initial <- if (method %in% c("transition", "both"))
    .wtna_initial_probs(df, codes, actor) else NULL

  if (method == "both") {
    result <- list(
      transition   = .wtna_finalize(weights$transition,   type, codes, data, "transition",
                                    initial = initial, window_size = window_size,
                                    mode = mode, actor = actor),
      cooccurrence = .wtna_finalize(weights$cooccurrence, type, codes, data, "cooccurrence",
                                    window_size = window_size, mode = mode, actor = actor),
      method = "wtna_both"
    )
    class(result) <- "wtna_mixed"
    return(result)
  }

  .wtna_finalize(weights, type, codes, data, method, initial = initial,
                 window_size = window_size, mode = mode, actor = actor)
}


# ---- Private helpers ----

#' Prepare one-hot code matrix
#'
#' Missing one-hot indicators are treated as inactive cells. This keeps sparse
#' indicator exports from propagating NA through cross-product counts.
#' @noRd
.wtna_code_matrix <- function(df, codes) {
  X <- as.matrix(df[, codes, drop = FALSE])
  X[is.na(X)] <- 0L
  storage.mode(X) <- "integer"
  X
}


#' Compute transition counts
#'
#' For window_size <= 1: consecutive crossprod (t->t+1).
#' For window_size > 1: pairwise between-window transitions matching tna's
#' windowed algorithm — every position in window_i paired with every position
#' in window_{i+1}.
#' @noRd
.wtna_transitions <- function(X, window_size = 3L, mode = "non-overlapping") {
  n <- nrow(X)
  k <- ncol(X)
  if (n < 2L) return(matrix(0, k, k))

  if (window_size <= 1L) {
    return(crossprod(X[-n, , drop = FALSE], X[-1L, , drop = FALSE]))
  }

  weights <- matrix(0, k, k)

  if (mode == "non-overlapping") {
    # Match tna's compute_transitions_windowed: pair every position in
    # window_i with every position in window_{i+1}
    divides <- n %% window_size == 0L
    q <- n %/% window_size - 1L * divides
    for (i in seq_len(q)) {
      j_idx <- seq((i - 1L) * window_size + 1L, i * window_size)
      k_idx <- seq(i * window_size + 1L, min(n, (i + 1L) * window_size))
      # sum_j sum_ki tcrossprod(X[j,], X[ki,]) = tcrossprod(colSums, colSums)
      s_curr <- colSums(X[j_idx, , drop = FALSE])
      s_next <- colSums(X[k_idx, , drop = FALSE])
      weights <- weights + tcrossprod(s_curr, s_next)
    }
  } else {
    # Overlapping: consecutive windows shifted by 1
    n_windows <- n - window_size + 1L
    if (n_windows < 2L) return(weights)
    for (i in seq_len(n_windows - 1L)) {
      j_idx <- seq(i, i + window_size - 1L)
      k_idx <- seq(i + 1L, i + window_size)
      s_curr <- colSums(X[j_idx, , drop = FALSE])
      s_next <- colSums(X[k_idx, , drop = FALSE])
      weights <- weights + tcrossprod(s_curr, s_next)
    }
  }
  weights
}


#' Compute co-occurrence counts
#'
#' For window_size <= 1: standard crossprod across all rows.
#' For window_size > 1: within-window pairwise co-occurrence matching tna's
#' windowed algorithm — every position in a window paired with every other
#' position in the same window.
#' @noRd
.wtna_cooccurrence <- function(X, window_size = 3L, mode = "non-overlapping") {
  n <- nrow(X)
  k <- ncol(X)

  if (window_size <= 1L) return(crossprod(X))

  weights <- matrix(0, k, k)

  if (mode == "non-overlapping") {
    # Match tna: pair every position with every other position within
    # the same window (including self-pairs at same position)
    n_windows <- ceiling(n / window_size)
    for (i in seq_len(n_windows)) {
      idx <- seq((i - 1L) * window_size + 1L, min(n, i * window_size))
      s <- colSums(X[idx, , drop = FALSE])
      weights <- weights + tcrossprod(s)
    }
  } else {
    n_windows <- n - window_size + 1L # nocov start
    if (n_windows < 1L) return(weights)
    for (i in seq_len(n_windows)) {
      idx <- seq(i, i + window_size - 1L)
      s <- colSums(X[idx, , drop = FALSE])
      weights <- weights + tcrossprod(s) # nocov end
    }
  }
  weights
}


#' Dispatch weight computation by method
#'
#' @param X_raw Raw one-hot matrix (not collapsed). Transitions use this
#'   directly with pairwise between-window counting. Co-occurrence collapses
#'   windows first via \code{.wtna_to_matrix}.
#' @noRd
.wtna_compute_weights <- function(X_raw, method, window_size = 3L,
                                   mode = "non-overlapping") {
  switch(method,
    transition = .wtna_transitions(X_raw, window_size, mode),
    cooccurrence = .wtna_cooccurrence(X_raw, window_size, mode),
    both = list(
      transition = .wtna_transitions(X_raw, window_size, mode),
      cooccurrence = .wtna_cooccurrence(X_raw, window_size, mode)
    )
  )
}


#' Compute weights per actor and sum
#' @noRd
.wtna_compute_by_actor <- function(df, codes, window_size, mode, actor, method) {
  n_codes <- length(codes)
  init_mat <- matrix(0, n_codes, n_codes)
  if (length(actor) == 1L) {
    groups <- split(df, df[[actor]])
  } else {
    grp_key <- interaction(df[, actor, drop = FALSE], drop = TRUE)
    groups <- split(df, grp_key)
  }

  matrices <- lapply(groups, function(g) {
    X_raw <- .wtna_code_matrix(g, codes)
    .wtna_compute_weights(X_raw, method, window_size, mode)
  })

  if (method == "both") {
    list(
      transition = Reduce(`+`, lapply(matrices, `[[`, "transition"), init_mat),
      cooccurrence = Reduce(`+`, lapply(matrices, `[[`, "cooccurrence"), init_mat)
    )
  } else {
    Reduce(`+`, matrices, init_mat)
  }
}


#' Resolve codes specification to column names
#'
#' Supports multiple selection styles:
#' \itemize{
#'   \item \code{NULL} — auto-detect binary 0/1 columns (excluding \code{exclude})
#'   \item Character vector — column names (e.g. \code{c("A", "B", "C")})
#'   \item Numeric vector — column indices (e.g. \code{2:9})
#'   \item Single string with \code{:} — column name range
#'     (e.g. \code{"Planning:Evaluating"})
#' }
#' @param df Data frame.
#' @param codes Column specification (NULL, character, numeric, or range string).
#' @param exclude Character vector of column names to exclude from auto-detection.
#' @return Character vector of resolved column names.
#' @noRd
.resolve_codes <- function(df, codes = NULL, exclude = NULL) {
  if (is.null(codes)) {
    # Auto-detect binary columns, excluding specified columns
    check_df <- df[, setdiff(names(df), exclude), drop = FALSE]
    return(.wtna_auto_detect_codes(check_df))
  }

  if (is.numeric(codes)) {
    # Column indices
    stopifnot(all(codes >= 1L & codes <= ncol(df)))
    return(names(df)[codes])
  }

  if (is.character(codes) && length(codes) == 1L && grepl(":", codes)) {
    # Name range "start:end"
    parts <- strsplit(codes, ":")[[1]]
    stopifnot(length(parts) == 2L)
    start <- match(parts[1], names(df))
    end <- match(parts[2], names(df))
    if (is.na(start)) stop("Column '", parts[1], "' not found.")
    if (is.na(end)) stop("Column '", parts[2], "' not found.")
    return(names(df)[start:end])
  }

  # Character vector of names
  stopifnot(is.character(codes))
  missing <- setdiff(codes, names(df))
  if (length(missing) > 0) {
    stop("Columns not found: ", paste(missing, collapse = ", "))
  }
  codes
}


#' Auto-detect one-hot binary columns
#' @noRd
.wtna_auto_detect_codes <- function(df) {
  is_onehot <- vapply(df, function(x) {
    if (!is.numeric(x) && !is.logical(x)) return(FALSE)
    vals <- unique(x[!is.na(x)])
    length(vals) > 0 && all(vals %in% c(0, 1))
  }, logical(1L))

  codes <- names(df)[is_onehot]
  if (length(codes) == 0L) stop("No one-hot columns found.") # nocov
  codes
}


#' Compute initial state probabilities for wtna transition networks
#'
#' For each actor (or the whole dataset if no actor), finds the first active
#' row and distributes probability uniformly across all simultaneously active
#' states. Results are averaged across actors to produce a distribution summing
#' to 1. This handles binary data where multiple states can be active at once.
#' @noRd
.wtna_initial_probs <- function(df, codes, actor) {
  X <- .wtna_code_matrix(df, codes)

  .first_row_probs <- function(mat) {
    active_rows <- which(rowSums(mat) > 0L)
    if (length(active_rows) == 0L) return(NULL)
    active_cols <- which(mat[active_rows[1L], ] > 0L)
    if (length(active_cols) == 0L) return(NULL)
    init <- setNames(numeric(length(codes)), codes)
    init[active_cols] <- 1.0 / length(active_cols)
    init
  }

  if (is.null(actor)) {
    return(.first_row_probs(X))
  }

  grp       <- df[[actor[1L]]]
  actor_ids <- unique(grp)

  prob_accum <- setNames(numeric(length(codes)), codes)
  n_valid    <- 0L

  for (id in actor_ids) {
    p <- .first_row_probs(X[grp == id, , drop = FALSE])
    if (is.null(p)) next
    prob_accum <- prob_accum + p
    n_valid    <- n_valid + 1L
  }

  if (n_valid == 0L) return(NULL)
  prob <- prob_accum / n_valid
  prob / sum(prob)   # ensure exact sum = 1 despite floating point
}

#' Finalize: row-normalize and build netobject
#' @noRd
.wtna_finalize <- function(weights, type, codes, data, method, initial = NULL,
                            window_size = 3L, mode = "non-overlapping",
                            actor = NULL) {
  if (type == "relative") {
    rs <- rowSums(weights)
    rs[rs == 0] <- 1
    weights <- weights / rs
  }

  dimnames(weights) <- list(codes, codes)
  directed <- method == "transition"

  # Extract edges
  edges <- .extract_edges_from_matrix(weights, directed = directed)

  nodes_df <- data.frame(
    id = seq_along(codes),
    label = codes,
    name = codes,
    x = NA_real_,
    y = NA_real_,
    stringsAsFactors = FALSE
  )
  wtna_method <- paste0("wtna_", method)

  structure(
    list(
      data = data,
      weights = weights,
      nodes = nodes_df,
      edges = edges,
      directed = directed,
      method = wtna_method,
      params = list(type = type, window_size = window_size, mode = mode,
                    codes = codes, actor = actor),
      scaling = NULL,
      threshold = 0,
      n_nodes = length(codes),
      n_edges = nrow(edges),
      level = NULL,
      initial = initial,
      meta = list(source = "nestimate", layout = NULL,
                  tna = list(method = wtna_method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}


# ---- Registry estimator wrappers ----
# These allow wtna to be used via build_network() and gain
# bootstrap_network() / permutation() support automatically.

#' Core wtna estimator for the registry
#'
#' @param data Data frame with one-hot columns.
#' @param codes Character vector or NULL. One-hot column names.
#' @param window_size Integer. Window size. Default: 1.
#' @param mode Character. "non-overlapping" or "overlapping".
#' @param actor Character or NULL. Actor grouping column.
#' @param wtna_method Character. "transition" or "cooccurrence".
#' @param ... Ignored.
#' @return Standard estimator list (matrix, nodes, directed, cleaned_data).
#' @noRd
.estimator_wtna_core <- function(data, codes = NULL, window_size = 3L,
                                  mode = "non-overlapping", actor = NULL,
                                  wtna_method = "transition",
                                  type = "frequency", ...) {
  df <- as.data.frame(data)
  codes <- .resolve_codes(df, codes, exclude = actor)

  stopifnot(length(codes) >= 2L)
  mode <- match.arg(mode, c("non-overlapping", "overlapping"))
  window_size <- as.integer(window_size)

  if (is.null(actor)) {
    X_raw <- .wtna_code_matrix(df, codes)
    weights <- .wtna_compute_weights(X_raw, wtna_method, window_size, mode)
  } else {
    stopifnot(all(actor %in% names(df)))
    weights <- .wtna_compute_by_actor(df, codes, window_size, mode,
                                       actor, wtna_method)
  }

  if (type == "relative") {
    rs <- rowSums(weights)
    rs[rs == 0] <- 1
    weights <- weights / rs
  }

  dimnames(weights) <- list(codes, codes)
  directed <- wtna_method == "transition"

  list(
    matrix = weights,
    nodes = codes,
    directed = directed,
    cleaned_data = data
  )
}


#' WTNA transition estimator (directed)
#'
#' Accepts both one-hot binary data and sequence data (wide/long).
#' One-hot: windowing + crossprod. Sequence: existing transition counting.
#' @noRd
.estimator_wtna <- function(data,
                             format = "auto",
                             action = "Action",
                             id = NULL,
                             time = "Time",
                             cols = NULL,
                             codes = NULL,
                             window_size = 3L,
                             mode = "non-overlapping",
                             actor = NULL,
                             ...) {
  df <- as.data.frame(data)

  # Detect format: one-hot vs sequence
  if (format == "auto") {
    format <- if (action %in% names(df)) "long" else "wide"
  }

  if (format == "wide") {
    state_cols <- .select_state_cols(df, c(id, actor), cols)
    is_binary <- length(state_cols) >= 2L && all(vapply(
      df[, state_cols, drop = FALSE],
      function(x) is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)),
      logical(1)
    ))

    if (is_binary || !is.null(codes)) {
      # One-hot path: windowing + crossprod
      return(.estimator_wtna_core(
        data, codes = codes, window_size = window_size,
        mode = mode, actor = actor, wtna_method = "transition", ...
      ))
    }
  }

  # Sequence path: use existing transition counting directly
  .estimator_frequency(
    data, format = format, action = action, id = id,
    time = time, cols = cols, ...
  )
}




#' Print Method for wtna_mixed
#'
#' @param x A \code{wtna_mixed} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @examples
#' oh <- matrix(c(1,0,0, 0,1,0, 0,0,1, 1,0,0), nrow = 4, byrow = TRUE,
#'              dimnames = list(NULL, c("A","B","C")))
#' mixed <- wtna(oh, method = "both")
#' print(mixed)
#'
#' \donttest{
#' oh <- data.frame(
#'   A = c(1,0,1,0,1,0,1,0),
#'   B = c(0,1,0,1,0,1,0,1),
#'   C = c(1,1,0,0,1,1,0,0)
#' )
#' mixed <- wtna(oh, method = "both")
#' print(mixed)
#' }
#' @export
print.wtna_mixed <- function(x, ...) {
  cat("Mixed Window TNA (transition + co-occurrence)\n")
  cat("-- Transition (directed) --\n")
  t <- x$transition
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", t$n_nodes, t$n_edges))
  cat("-- Co-occurrence (undirected) --\n")
  co <- x$cooccurrence
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", co$n_nodes, co$n_edges))
  invisible(x)
}
