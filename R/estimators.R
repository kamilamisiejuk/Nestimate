# ---- Built-in Estimator Implementations ----

# ---- Shared helpers ----

# TraMineR/tna void and missing markers to treat as NA
.void_markers <- c("%", "*", "", "NA", "NaN")

#' Remove void/missing markers from a character vector
#'
#' Strips TraMineR void (\code{\%}), missing (\code{*}), empty strings,
#' \code{"NA"}, \code{"NaN"}, and actual \code{NA}s from state values.
#'
#' @param x Character vector of state values.
#' @return Character vector with void/missing markers removed.
#' @noRd
.clean_states <- function(x) {
  x[x %in% .void_markers] <- NA
  x
}

#' Select state columns from a data frame
#'
#' Resolves which columns contain state data based on explicit \code{cols},
#' exclusion of \code{id} columns, or all columns as fallback.
#'
#' @param data Data frame.
#' @param id Character vector or NULL. ID column(s) to exclude.
#' @param cols Character vector or NULL. Explicit state columns.
#' @return Character vector of state column names.
#' @noRd
.select_state_cols <- function(data, id = NULL, cols = NULL) {
  if (!is.null(cols)) {
    cols
  } else if (!is.null(id)) {
    setdiff(names(data), id)
  } else {
    names(data)
  }
}

#' Prepare wide sequence data with tna-compatible sequence options
#' @noRd
.prepare_sequence_wide <- function(data, id = NULL, cols = NULL,
                                   alphabet = NULL, concat = 1L,
                                   begin_state = NULL, end_state = NULL) {
  state_cols <- .select_state_cols(data, id, cols)
  missing_cols <- setdiff(state_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found: ", paste(missing_cols, collapse = ", "))
  }
  if (length(state_cols) < 1L) {
    stop("At least 1 state column is required for wide format.")
  }

  mat <- as.matrix(data[, state_cols, drop = FALSE])
  mat[] <- .clean_states(mat)

  if (!is.null(alphabet)) {
    alphabet <- as.character(alphabet)
  } else {
    vals <- as.vector(mat)
    alphabet <- sort(unique(vals[!is.na(vals)]))
  }

  concat <- as.integer(concat %||% 1L)
  if (length(concat) != 1L || is.na(concat) || concat < 1L) {
    stop("'concat' must be a positive integer.", call. = FALSE)
  }
  if (concat > 1L && nrow(mat) > 0L) {
    n <- nrow(mat)
    k <- ncol(mat)
    m <- n %/% concat
    extra <- n %% concat
    out <- matrix(NA_character_, nrow = m + as.integer(extra > 0L),
                  ncol = concat * k)
    if (m > 0L) {
      idx <- rep(seq_len(m), each = concat)
      block <- split(seq_len(m * concat), idx)
      out[seq_len(m), ] <- do.call(rbind, lapply(block, function(ii) {
        as.vector(t(mat[ii, , drop = FALSE]))
      }))
    }
    if (extra > 0L) {
      rows <- seq.int(m * concat + 1L, n)
      out[m + 1L, seq_len(length(rows) * k)] <- as.vector(t(mat[rows, , drop = FALSE]))
    }
    mat <- out
  }

  if (!is.null(begin_state)) {
    begin_state <- as.character(begin_state)[1L]
    mat <- cbind(begin_state, mat)
    alphabet <- c(begin_state, alphabet)
  }
  if (!is.null(end_state) && nrow(mat) > 0L) {
    end_state <- as.character(end_state)[1L]
    last_obs <- max.col(!is.na(mat), ties.method = "last")
    mat <- cbind(mat, NA_character_)
    has_obs <- rowSums(!is.na(mat[, -ncol(mat), drop = FALSE])) > 0L
    idx <- cbind(which(has_obs), last_obs[has_obs] + 1L)
    mat[idx] <- end_state
    alphabet <- c(alphabet, end_state)
  }

  alphabet <- unique(alphabet)
  list(matrix = mat, states = alphabet)
}

#' Validate or transform matrix input for transition estimators
#' @noRd
.transition_matrix_input <- function(data, type = c("relative", "frequency",
                                                    "co_occurrence")) {
  type <- match.arg(type)
  x <- try(data.matrix(data), silent = TRUE)
  if (inherits(x, "try-error")) {
    stop("'data' must be coercible to a numeric matrix.", call. = FALSE)
  }
  if (nrow(x) != ncol(x) || nrow(x) < 1L) {
    stop("'data' must be a square matrix.", call. = FALSE)
  }
  if (is.null(colnames(x))) {
    dimnames(x) <- list(as.character(seq_len(ncol(x))),
                        as.character(seq_len(ncol(x))))
  } else if (is.null(rownames(x))) {
    rownames(x) <- colnames(x)
  }
  if (type == "relative") {
    rs <- rowSums(x)
    if (!all(rs > 0)) {
      stop("At least one element of each row of 'data' must be positive.",
           call. = FALSE)
    }
    x <- x / rs
  } else {
    xi <- as.integer(x)
    if (!all(x == xi) || any(x < 0)) {
      stop("All elements of 'data' must be non-negative integers.",
           call. = FALSE)
    }
    x[] <- xi
  }
  x
}


# ---- Core transition counting engine ----

#' Count transitions from sequence data
#'
#' Dispatches to wide or long format counting. Returns a square integer matrix
#' of transition frequencies.
#'
#' @param data Data frame of sequence data.
#' @param format Character: \code{"auto"}, \code{"wide"}, or \code{"long"}.
#' @param action Character. Action column name (long format).
#' @param id Character vector. ID column(s).
#' @param time Character. Time column name (long format).
#' @param cols Character vector. State columns (wide format).
#'
#' @return Square integer matrix with row/column names = sorted unique states.
#' @noRd
.count_transitions <- function(data,
                               format = "auto",
                               action = "Action",
                               id = NULL,
                               time = "Time",
                               cols = NULL,
                               alphabet = NULL,
                               weighted = FALSE,
                               concat = 1L,
                               begin_state = NULL,
                               end_state = NULL) {
  stopifnot(is.data.frame(data))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "wide") {
    .count_transitions_wide(
      data, id = id, cols = cols, alphabet = alphabet,
      weighted = weighted, concat = concat,
      begin_state = begin_state, end_state = end_state
    )
  } else {
    .count_transitions_long(
      data, action = action, id = id, time = time, alphabet = alphabet,
      weighted = weighted, begin_state = begin_state, end_state = end_state
    )
  }
}


#' Count transitions from wide format (vectorized base R)
#'
#' Each row is a sequence, columns are consecutive time points.
#' Uses matrix slicing + tabulate for speed.
#'
#' @noRd
.count_transitions_wide <- function(data, id = NULL, cols = NULL,
                                    alphabet = NULL, weighted = FALSE,
                                    concat = 1L, begin_state = NULL,
                                    end_state = NULL) {
  prepared <- .prepare_sequence_wide(
    data, id = id, cols = cols, alphabet = alphabet, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  mat <- prepared$matrix
  all_states <- prepared$states
  state_cols <- seq_len(ncol(mat))
  windowed <- isTRUE(attr(data, "windowed"))
  window_size <- attr(data, "window_size")
  window_span <- attr(data, "window_span")
  if (length(state_cols) < 2L) {
    n_states <- length(all_states)
    return(matrix(0L, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }
  nc <- ncol(mat)

  # from = all columns except last, to = all columns except first
  n_states <- length(all_states)
  if (n_states == 0L) {
    return(matrix(0L, 0, 0))
  }

  if (windowed) {
    window_size <- if (is.null(window_size)) 2L else as.integer(window_size)
    window_span <- if (is.null(window_span)) 1L else as.integer(window_span)
    effective_window <- window_size * window_span
    divides <- nc %% effective_window == 0L
    q <- nc %/% effective_window - as.integer(divides)
    out <- if (isTRUE(weighted)) {
      matrix(0, nrow = n_states, ncol = n_states,
             dimnames = list(all_states, all_states))
    } else {
      matrix(0L, nrow = n_states, ncol = n_states,
             dimnames = list(all_states, all_states))
    }
    if (q <= 0L) return(out)
    # Per-row weights mirror the non-windowed branch: each transition
    # observed in a sequence of length L gets weight 1/L. seq_lengths is
    # computed once for the whole matrix so the windowed and non-windowed
    # paths agree on row weight when the same row is fed through either.
    if (isTRUE(weighted)) {
      seq_lengths <- rowSums(!is.na(mat))
    }
    nbins <- n_states * n_states

    for (i in seq_len(q)) {
      j_idx <- seq((i - 1L) * effective_window + 1L,
                   i * effective_window)
      k_idx <- seq(i * effective_window + 1L,
                   min(nc, (i + 1L) * effective_window))
      for (j in j_idx) {
        for (k in k_idx) {
          from <- mat[, j]
          to <- mat[, k]
          valid <- !is.na(from) & !is.na(to)
          if (!any(valid)) next
          from_int <- match(from[valid], all_states)
          to_int <- match(to[valid], all_states)
          pair_idx <- (from_int - 1L) * n_states + to_int
          if (isTRUE(weighted)) {
            w <- (1 / seq_lengths)[valid]
            add <- rowsum(w, pair_idx, reorder = FALSE)
            tmp <- numeric(nbins)
            tmp[as.integer(rownames(add))] <- add[, 1L]
            out <- out + matrix(tmp, nrow = n_states, ncol = n_states,
                                byrow = TRUE)
          } else {
            counts <- tabulate(pair_idx, nbins = nbins)
            out <- out + matrix(counts, nrow = n_states, ncol = n_states,
                                byrow = TRUE)
          }
        }
      }
    }
    return(out)
  }

  from_mat <- mat[, -nc, drop = FALSE]
  to_mat <- mat[, -1L, drop = FALSE]
  valid <- !is.na(from_mat) & !is.na(to_mat)
  if (!any(valid)) {
    return(matrix(0L, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }
  from_vec <- as.vector(from_mat[valid])
  to_vec <- as.vector(to_mat[valid])
  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int
  if (isTRUE(weighted)) {
    seq_lengths <- rowSums(!is.na(mat))
    w_mat <- matrix(1 / seq_lengths, nrow = nrow(mat), ncol = nc - 1L)
    weights <- as.vector(w_mat[valid])
    counts <- rowsum(weights, pair_idx, reorder = FALSE)
    out <- numeric(n_states * n_states)
    out[as.integer(rownames(counts))] <- counts[, 1L]
    counts <- out
  } else {
    counts <- tabulate(pair_idx, nbins = n_states * n_states)
  }

  matrix(
    counts,
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


#' Count transitions from long format (data.table)
#'
#' Uses data.table for fast grouping and lag computation.
#'
#' @importFrom data.table setDT setorderv
#' @noRd
.count_transitions_long <- function(data, action = "Action", id = NULL,
                                    time = "Time", alphabet = NULL,
                                    weighted = FALSE, begin_state = NULL,
                                    end_state = NULL) {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }
  if (!is.null(id)) {
    missing_ids <- setdiff(id, names(data))
    if (length(missing_ids) > 0) {
      stop("ID column(s) not found: ", paste(missing_ids, collapse = ", "))
    }
  }

  dt <- data.table::as.data.table(data)

  # Clean void/missing markers in action column
  data.table::set(dt, j = action, value = .clean_states(as.character(dt[[action]])))

  # Preserve original row order as tiebreaker for duplicate timestamps
  data.table::set(dt, j = ".orig_row", value = seq_len(nrow(dt)))

  # Order by id + time + original row order
  order_cols <- c(id, if (time %in% names(dt)) time, ".orig_row")
  data.table::setorderv(dt, order_cols)

  action_col <- action

  # Build group key for sequences
  if (is.null(id)) {
    # Single sequence: consecutive pairs from all rows
    actions <- dt[[action_col]]
    if (!is.null(begin_state)) actions <- c(as.character(begin_state)[1L], actions)
    if (!is.null(end_state)) actions <- c(actions, as.character(end_state)[1L])
    n <- length(actions)
    if (n < 2L) {
      all_states <- alphabet %||% sort(unique(c(actions[!is.na(actions)],
                                                begin_state, end_state)))
      n_states <- length(all_states)
      return(matrix(0L, nrow = n_states, ncol = n_states,
                    dimnames = list(all_states, all_states)))
    }
    from_vec <- actions[-n]
    to_vec <- actions[-1L]
    # Filter out pairs where either side is NA
    valid <- !is.na(from_vec) & !is.na(to_vec)
    from_vec <- from_vec[valid]
    to_vec <- to_vec[valid]
  } else {
    # Group by ID columns, extract consecutive pairs
    # Use data.table's fast grouping
    if (length(id) == 1L) {
      grp_col <- id
    } else {
      # Create composite group key
      dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")),
         .SDcols = id]
      grp_col <- ".grp_key"
    }

    # Extract from/to pairs per group using data.table
    # NAs are kept in position; pairs with NA on either side are filtered out
    pairs <- dt[, {
      a <- get(action_col)
      if (!is.null(begin_state)) a <- c(as.character(begin_state)[1L], a)
      if (!is.null(end_state)) a <- c(a, as.character(end_state)[1L])
      n <- length(a)
      if (n < 2L) {
        list(from = character(0), to = character(0), weight = numeric(0))
      } else {
        f <- a[-n]
        t <- a[-1L]
        ok <- !is.na(f) & !is.na(t)
        w <- if (isTRUE(weighted)) rep(1 / sum(!is.na(a)), sum(ok)) else rep(1, sum(ok))
        list(from = f[ok], to = t[ok], weight = w)
      }
    }, by = grp_col]

    from_vec <- pairs$from
    to_vec <- pairs$to
    pair_weights <- pairs$weight
  }

  if (length(from_vec) == 0L) {
    # Collect all states to set matrix dimensions
    all_states <- alphabet %||% sort(unique(c(dt[[action_col]][!is.na(dt[[action_col]])],
                                              begin_state, end_state)))
    n_states <- length(all_states)
    return(matrix(0L, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }

  # Integer encode + tabulate
  all_states <- alphabet %||% sort(unique(c(from_vec, to_vec,
                                            begin_state, end_state)))
  n_states <- length(all_states)

  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int

  if (isTRUE(weighted)) {
    if (!exists("pair_weights")) {
      pair_weights <- rep(1 / sum(!is.na(actions)), length(pair_idx))
    }
    rs <- rowsum(pair_weights, pair_idx, reorder = FALSE)
    counts <- numeric(n_states * n_states)
    counts[as.integer(rownames(rs))] <- rs[, 1L]
  } else {
    counts <- tabulate(pair_idx, nbins = n_states * n_states)
  }

  matrix(
    counts,
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


# ---- Transition estimators ----

#' Compute initial state probability distribution
#'
#' Returns a named numeric vector (summing to 1) of the proportion of sequences
#' starting in each state. States present in the transition matrix but never
#' observed as a first state receive probability 0.
#'
#' @param data Data frame of sequence data.
#' @param states Character vector of all state names (from the transition matrix).
#' @param format Character: \code{"wide"} or \code{"long"}.
#' @param action Character. Action column name (long format).
#' @param id Character vector or NULL. ID column(s).
#' @param time Character. Time column name (long format).
#' @param cols Character vector or NULL. State columns (wide format).
#'
#' @return Named numeric vector of initial probabilities aligned to \code{states}.
#' @noRd
.compute_initial_probs <- function(data, states,
                                    format = "wide",
                                    action = "Action",
                                    id = NULL,
                                    time = "Time",
                                    cols = NULL,
                                    concat = 1L,
                                    begin_state = NULL,
                                    end_state = NULL) {
  if (format == "wide") {
    prepared <- .prepare_sequence_wide(
      data, id = id, cols = cols, alphabet = states, concat = concat,
      begin_state = begin_state, end_state = end_state
    )
    mat <- prepared$matrix
    first_states <- apply(mat, 1L, function(r) {
      r <- r[!is.na(r)]
      if (length(r)) r[[1L]] else NA_character_
    })
  } else {
    dt <- data.table::as.data.table(data)
    data.table::set(dt, j = action, value = .clean_states(as.character(dt[[action]])))
    data.table::set(dt, j = ".orig_row", value = seq_len(nrow(dt)))
    order_cols <- c(id, if (time %in% names(dt)) time, ".orig_row")
    data.table::setorderv(dt, order_cols)
    if (is.null(id)) {
      vals <- dt[[action]]
      if (!is.null(begin_state)) vals <- c(as.character(begin_state)[1L], vals)
      if (!is.null(end_state)) vals <- c(vals, as.character(end_state)[1L])
      vals <- vals[!is.na(vals)]
      first_states <- if (length(vals)) vals[[1L]] else character(0)
    } else {
      grp_col <- if (length(id) == 1L) {
        id
      } else {
        dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")), .SDcols = id]
        ".grp_key"
      }
      firsts <- dt[, {
        a <- get(action)
        if (!is.null(begin_state)) a <- c(as.character(begin_state)[1L], a)
        if (!is.null(end_state)) a <- c(a, as.character(end_state)[1L])
        a <- a[!is.na(a)]
        list(first = if (length(a)) a[[1L]] else NA_character_)
      }, by = grp_col]
      first_states <- firsts$first
    }
  }

  first_states <- first_states[!is.na(first_states)]
  initial <- setNames(numeric(length(states)), states)
  if (length(first_states) == 0L) return(initial)
  counts <- tabulate(match(first_states, states), nbins = length(states))
  total <- sum(counts)
  if (total > 0L) initial <- counts / total
  names(initial) <- states
  initial
}


#' Frequency estimator: raw transition counts
#' @noRd
.estimator_frequency <- function(data,
                                 format = "auto",
                                 action = "Action",
                                 id = NULL,
                                 time = "Time",
                                 cols = NULL,
                                 codes = NULL,
                                 window_size = 3L,
                                 mode = "non-overlapping",
                                 actor = NULL,
                                 alphabet = NULL,
                                 weighted = FALSE,
                                 concat = 1L,
                                 begin_state = NULL,
                                 end_state = NULL,
                                 ...) {
  if (is.matrix(data)) {
    freq_mat <- .transition_matrix_input(data, type = "frequency")
    states <- rownames(freq_mat)
    return(list(matrix = freq_mat, nodes = states, directed = TRUE,
                cleaned_data = NULL, frequency_matrix = freq_mat))
  }
  if (identical(format, "onehot") || !is.null(codes)) {
    out <- .estimator_wtna_core(
      data, codes = codes, window_size = window_size, mode = mode,
      actor = actor, wtna_method = "transition", type = "frequency"
    )
    out$frequency_matrix <- out$matrix
    return(out)
  }
  freq_mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols,
    alphabet = alphabet, weighted = weighted, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  states <- rownames(freq_mat)
  resolved_format <- if (format == "auto") {
    if (action %in% names(data)) "long" else "wide"
  } else format
  initial <- .compute_initial_probs(
    data, states, format = resolved_format,
    action = action, id = id, time = time, cols = cols, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  list(
    matrix = freq_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = freq_mat,
    initial = initial
  )
}


#' Relative estimator: row-normalized transition probabilities
#' @noRd
.estimator_relative <- function(data,
                                format = "auto",
                                action = "Action",
                                id = NULL,
                                time = "Time",
                                cols = NULL,
                                codes = NULL,
                                window_size = 1L,
                                mode = "non-overlapping",
                                actor = NULL,
                                alphabet = NULL,
                                weighted = FALSE,
                                concat = 1L,
                                begin_state = NULL,
                                end_state = NULL,
                                ...) {
  if (is.matrix(data)) {
    rel_mat <- .transition_matrix_input(data, type = "relative")
    states <- rownames(rel_mat)
    return(list(matrix = rel_mat, nodes = states, directed = TRUE,
                cleaned_data = NULL))
  }
  if (identical(format, "onehot") || !is.null(codes)) {
    out <- .estimator_wtna_core(
      data, codes = codes, window_size = window_size, mode = mode,
      actor = actor, wtna_method = "transition", type = "relative"
    )
    out$frequency_matrix <- .estimator_wtna_core(
      data, codes = codes, window_size = window_size, mode = mode,
      actor = actor, wtna_method = "transition", type = "frequency"
    )$matrix
    return(out)
  }
  freq_mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols,
    alphabet = alphabet, weighted = weighted, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  states <- rownames(freq_mat)

  # Row-normalize
  row_sums <- rowSums(freq_mat)
  rel_mat <- freq_mat
  storage.mode(rel_mat) <- "double"
  nonzero <- row_sums > 0
  rel_mat[nonzero, ] <- rel_mat[nonzero, ] / row_sums[nonzero]

  resolved_format <- if (format == "auto") {
    if (action %in% names(data)) "long" else "wide"
  } else format
  initial <- .compute_initial_probs(
    data, states, format = resolved_format,
    action = action, id = id, time = time, cols = cols, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  list(
    matrix = rel_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    frequency_matrix = freq_mat,
    initial = initial
  )
}


#' Co-occurrence estimator: positional co-occurrence within sequences
#'
#' Counts all positional column pairs (i, j) where i < j. For each pair,
#' if both positions have non-NA states, the co-occurrence count is incremented
#' for both (from, to) and (to, from). This matches tna::ctna() semantics.
#'
#' When the input is one-hot binary data (all values 0/1), automatically
#' uses crossprod-based co-occurrence counting instead. This means
#' \code{method = "cna"} works for both sequence and one-hot data.
#'
#' @noRd
.estimator_co_occurrence <- function(data,
                                     format = "auto",
                                     action = "Action",
                                     id = NULL,
                                     time = "Time",
                                     cols = NULL,
                                     codes = NULL,
                                     window_size = 3L,
                                     mode = "non-overlapping",
                                     actor = NULL,
                                     alphabet = NULL,
                                     weighted = FALSE,
                                     concat = 1L,
                                     begin_state = NULL,
                                     end_state = NULL,
                                     ...) {
  if (is.matrix(data)) {
    cooc_mat <- .transition_matrix_input(data, type = "co_occurrence")
    return(list(matrix = cooc_mat, nodes = rownames(cooc_mat),
                directed = FALSE, cleaned_data = NULL))
  }
  stopifnot(is.data.frame(data))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  # One-hot binary input detection for wide format
  if (format == "wide") {
    # If codes explicitly provided, treat as one-hot
    if (!is.null(codes)) {
      return(.estimator_wtna_core(
        data, codes = codes, window_size = window_size,
        mode = mode, actor = actor,
        wtna_method = "cooccurrence", ...
      ))
    }

    # Auto-detect: check if all state columns are binary 0/1
    state_cols <- .select_state_cols(data, c(id, actor), cols)
    is_binary <- length(state_cols) >= 2L && all(vapply(
      data[, state_cols, drop = FALSE],
      function(x) is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)),
      logical(1)
    ))

    if (is_binary) { # nocov start
      return(.estimator_wtna_core(
        data, codes = state_cols, window_size = window_size,
        mode = mode, actor = actor,
        wtna_method = "cooccurrence", ...
      )) # nocov end
    }
  }

  # Sequence co-occurrence (column-pair counting)
  if (format == "wide") {
    cooc_mat <- .count_cooccurrence_wide(
      data, id = id, cols = cols, alphabet = alphabet, weighted = weighted,
      concat = concat, begin_state = begin_state, end_state = end_state
    )
  } else {
    cooc_mat <- .count_cooccurrence_long( # nocov start
      data, action = action, id = id, time = time, alphabet = alphabet,
      weighted = weighted, begin_state = begin_state, end_state = end_state
    ) # nocov end
  }

  list(
    matrix = cooc_mat,
    nodes = rownames(cooc_mat),
    directed = FALSE,
    cleaned_data = data
  )
}


#' Count positional co-occurrences from wide format
#'
#' For each pair of column positions (i, j) where i < j, counts how many
#' sequences have non-NA values at both positions. Symmetric: both (from, to)
#' and (to, from) are incremented.
#'
#' Strategy: integer-encode the matrix once, then iterate over column pairs
#' with lightweight tabulate accumulation. Avoids materializing huge
#' n_rows x n_pairs matrices.
#' @noRd
.count_cooccurrence_wide <- function(data, id = NULL, cols = NULL,
                                     alphabet = NULL, weighted = FALSE,
                                     concat = 1L, begin_state = NULL,
                                     end_state = NULL) {
  prepared <- .prepare_sequence_wide(
    data, id = id, cols = cols, alphabet = alphabet, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  mat <- prepared$matrix
  nc <- ncol(mat)
  all_states <- prepared$states
  n_states <- length(all_states)
  windowed <- isTRUE(attr(data, "windowed"))
  window_size <- attr(data, "window_size")
  window_span <- attr(data, "window_span")

  if (n_states == 0L || nc < 2L) {
    cooc <- matrix(0, n_states, n_states,
                   dimnames = list(all_states, all_states))
    return(cooc)
  }

  nbins <- n_states * n_states

  if (windowed) {
    window_size <- if (is.null(window_size)) 2L else as.integer(window_size)
    window_span <- if (is.null(window_span)) 1L else as.integer(window_span)
    effective_window <- window_size * window_span
    divides <- nc %% effective_window == 0L
    q <- nc %/% effective_window - as.integer(divides)
    counts <- if (isTRUE(weighted)) numeric(nbins) else integer(nbins)
    if (isTRUE(weighted)) {
      seq_lengths <- rowSums(!is.na(mat))
    }

    for (i in seq_len(q + 1L)) {
      j_idx <- seq((i - 1L) * effective_window + 1L,
                   min(nc, i * effective_window))
      for (j in j_idx) {
        from <- mat[, j]
        for (k in j_idx) {
          to <- mat[, k]
          valid <- !is.na(from) & !is.na(to)
          if (!any(valid)) next
          from_int <- match(from[valid], all_states)
          to_int <- match(to[valid], all_states)
          pair_idx <- (from_int - 1L) * n_states + to_int
          if (isTRUE(weighted)) {
            w <- (1 / seq_lengths)[valid]
            add <- rowsum(w, pair_idx, reorder = FALSE)
            tmp <- numeric(nbins)
            tmp[as.integer(rownames(add))] <- add[, 1L]
            counts <- counts + tmp
          } else {
            counts <- counts + tabulate(pair_idx, nbins)
          }
        }
      }
    }

    return(matrix(
      as.numeric(counts),
      nrow = n_states,
      ncol = n_states,
      byrow = TRUE,
      dimnames = list(all_states, all_states)
    ))
  }

  # Integer-encode the entire matrix once (NA stays NA)
  int_mat <- matrix(match(mat, all_states), nrow = nrow(mat), ncol = nc)

  # Accumulate counts across column pairs
  counts <- integer(nbins)

  for (i in seq_len(nc - 1L)) {
    col_i <- int_mat[, i]
    for (j in seq(i + 1L, nc)) {
      col_j <- int_mat[, j]
      valid <- !is.na(col_i) & !is.na(col_j)
      fi <- col_i[valid]
      tj <- col_j[valid]
      # Forward direction: i -> j
      idx_fwd <- (fi - 1L) * n_states + tj
      # Reverse direction: j -> i
      idx_rev <- (tj - 1L) * n_states + fi
      if (isTRUE(weighted)) {
        w <- (1 / rowSums(!is.na(mat)))[valid]
        add_f <- rowsum(w, idx_fwd, reorder = FALSE)
        add_r <- rowsum(w, idx_rev, reorder = FALSE)
        tmp <- numeric(nbins)
        tmp[as.integer(rownames(add_f))] <- tmp[as.integer(rownames(add_f))] + add_f[, 1L]
        tmp[as.integer(rownames(add_r))] <- tmp[as.integer(rownames(add_r))] + add_r[, 1L]
        counts <- counts + tmp
      } else {
        counts <- counts + tabulate(idx_fwd, nbins) + tabulate(idx_rev, nbins)
      }
    }
  }

  cooc <- matrix(
    as.numeric(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )

  # Self-pairs (A,A) are double-counted by the bidirectional approach
  diag(cooc) <- diag(cooc) / 2

  cooc
}


#' Count positional co-occurrences from long format
#'
#' Converts to wide-like structure per group, then counts column-pair
#' co-occurrences.
#' @noRd
.count_cooccurrence_long <- function(data, action = "Action", id = NULL,
                                     time = "Time", alphabet = NULL,
                                     weighted = FALSE, begin_state = NULL,
                                     end_state = NULL) {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }

  dt <- data.table::as.data.table(data)

  # Clean void/missing markers in action column
  data.table::set(dt, j = action, value = .clean_states(as.character(dt[[action]])))

  # Preserve original row order as tiebreaker
  data.table::set(dt, j = ".orig_row", value = seq_len(nrow(dt)))

  # Order by id + time + original row order
  order_cols <- c(id, if (time %in% names(dt)) time, ".orig_row")
  data.table::setorderv(dt, order_cols)

  # Build group key
  if (is.null(id)) {
    data.table::set(dt, j = ".seq_grp", value = rep(1L, nrow(dt)))
    grp_col <- ".seq_grp"
  } else if (length(id) == 1L) {
    grp_col <- id
  } else {
    dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")),
       .SDcols = id]
    grp_col <- ".grp_key"
  }

  action_col <- action

  # For each group, create all position pairs and collect (from, to)
  pairs <- dt[!is.na(get(action_col)), {
    a <- get(action_col)
    if (!is.null(begin_state)) a <- c(as.character(begin_state)[1L], a)
    if (!is.null(end_state)) a <- c(a, as.character(end_state)[1L])
    n <- length(a)
    if (n < 2L) {
      list(from = character(0), to = character(0), weight = numeric(0))
    } else {
      cp <- utils::combn(n, 2)
      f <- a[cp[1, ]]
      t <- a[cp[2, ]]
      # Both directions
      w <- if (isTRUE(weighted)) rep(1 / sum(!is.na(a)), 2L * length(f)) else rep(1, 2L * length(f))
      list(from = c(f, t), to = c(t, f), weight = w)
    }
  }, by = grp_col]

  from_vec <- pairs$from
  to_vec <- pairs$to
  pair_weights <- pairs$weight

  # Collect all states
  all_states <- alphabet %||% sort(unique(c(dt[[action_col]][!is.na(dt[[action_col]])],
                                            begin_state, end_state)))
  n_states <- length(all_states)

  if (length(from_vec) == 0L || n_states == 0L) {
    return(matrix(0, nrow = n_states, ncol = n_states,
                  dimnames = list(all_states, all_states)))
  }

  from_int <- match(from_vec, all_states)
  to_int <- match(to_vec, all_states)
  pair_idx <- (from_int - 1L) * n_states + to_int

  if (isTRUE(weighted)) {
    rs <- rowsum(pair_weights, pair_idx, reorder = FALSE)
    counts <- numeric(n_states * n_states)
    counts[as.integer(rownames(rs))] <- rs[, 1L]
  } else {
    counts <- tabulate(pair_idx, nbins = n_states * n_states)
  }

  cooc <- matrix(
    as.numeric(counts),
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )

  # Self-pairs (A,A) are double-counted by the bidirectional approach
  diag(cooc) <- diag(cooc) / 2

  cooc
}


# ---- Attention (decay-weighted) estimator ----

#' Count attention-weighted transitions from wide format
#'
#' For each pair of column positions (i, j), computes a decay weight based on
#' the temporal distance, then accumulates weighted counts per state pair.
#' Direction controls which pairs are considered.
#'
#' @param data Data frame with sequences in rows.
#' @param id Character vector or NULL. ID columns to exclude.
#' @param cols Character vector or NULL. State columns.
#' @param lambda Numeric. Decay rate parameter. Default: 1.
#' @param direction Character. "forward" (i < j), "backward" (i > j),
#'   or "both" (all pairs). Default: "forward".
#' @param decay Function or NULL. Custom decay function(ti, tj, lambda).
#'   Default: exp(-abs(ti - tj) / lambda).
#' @param time_matrix Matrix or NULL. Custom time values per cell.
#' @param duration Numeric vector or NULL. Per-column durations to convert
#'   to cumulative time.
#' @noRd
.count_attention_wide <- function(data, id = NULL, cols = NULL,
                                   lambda = 1, direction = "forward",
                                   decay = NULL, time_matrix = NULL,
                                   duration = NULL, alphabet = NULL,
                                   concat = 1L, begin_state = NULL,
                                   end_state = NULL) {
  prepared <- .prepare_sequence_wide(
    data, id = id, cols = cols, alphabet = alphabet, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  mat <- prepared$matrix
  nr <- nrow(mat)
  nc <- ncol(mat)

  all_states <- prepared$states
  n_states <- length(all_states)

  if (n_states == 0L || nc < 2L) {
    return(matrix(0, n_states, n_states,
                  dimnames = list(all_states, all_states)))
  }

  # Build time matrix
  if (!is.null(time_matrix)) {
    stopifnot(nrow(time_matrix) == nr, ncol(time_matrix) == nc)
    tmat <- time_matrix
  } else if (!is.null(duration)) {
    stopifnot(length(duration) == nc)
    cum_time <- cumsum(duration)
    tmat <- matrix(rep(cum_time, each = nr), nr, nc)
  } else {
    tmat <- matrix(rep(seq_len(nc), each = nr), nr, nc)
  }

  # Default decay function
  if (is.null(decay)) {
    decay <- function(ti, tj, lam) exp(-abs(ti - tj) / lam)
  }

  # Integer-encode states
  int_mat <- matrix(match(mat, all_states), nrow = nr, ncol = nc)

  # Accumulate weighted counts
  counts <- numeric(n_states * n_states)

  for (i in seq_len(nc)) {
    for (j in seq_len(nc)) {
      if (i == j) next
      if (direction == "forward" && i >= j) next
      if (direction == "backward" && i <= j) next

      col_i <- int_mat[, i]
      col_j <- int_mat[, j]
      valid <- !is.na(col_i) & !is.na(col_j)

      if (!any(valid)) next # nocov

      fi <- col_i[valid]
      tj <- col_j[valid]
      d <- decay(tmat[valid, i], tmat[valid, j], lambda)

      pair_idx <- (fi - 1L) * n_states + tj
      agg <- tapply(d, pair_idx, sum)
      idx <- as.integer(names(agg))
      counts[idx] <- counts[idx] + as.numeric(agg)
    }
  }

  matrix(
    counts,
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


#' Count attention-weighted transitions from long format
#'
#' Converts to per-group wide format, then applies attention counting.
#' @noRd
.count_attention_long <- function(data, action = "Action", id = NULL,
                                   time = "Time", lambda = 1,
                                   direction = "forward", decay = NULL,
                                   time_matrix = NULL, duration = NULL,
                                   alphabet = NULL, begin_state = NULL,
                                   end_state = NULL) {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }

  dt <- data.table::as.data.table(data)

  # Clean void/missing markers in action column
  data.table::set(dt, j = action, value = .clean_states(as.character(dt[[action]])))

  # Order by id + time
  order_cols <- c(id, if (time %in% names(dt)) time)
  if (length(order_cols) > 0) {
    data.table::setorderv(dt, order_cols)
  }

  # Build group key
  if (is.null(id)) {
    data.table::set(dt, j = ".seq_grp", value = rep(1L, nrow(dt)))
    grp_col <- ".seq_grp"
  } else if (length(id) == 1L) {
    grp_col <- id
  } else {
    dt[, .grp_key := do.call(paste, c(.SD, sep = "\x1f")),
       .SDcols = id]
    grp_col <- ".grp_key"
  }

  action_col <- action
  all_states <- alphabet %||% sort(unique(c(dt[[action_col]][!is.na(dt[[action_col]])],
                                            begin_state, end_state)))
  n_states <- length(all_states)

  if (n_states == 0L) {
    return(matrix(0, 0, 0))
  }

  # Default decay function
  if (is.null(decay)) {
    decay <- function(ti, tj, lam) exp(-abs(ti - tj) / lam)
  }

  # Process per group
  counts <- numeric(n_states * n_states)

  groups <- split(dt, dt[[grp_col]])
  for (g in groups) {
    a <- g[[action_col]]
    if (!is.null(begin_state)) a <- c(as.character(begin_state)[1L], a)
    if (!is.null(end_state)) a <- c(a, as.character(end_state)[1L])
    n <- length(a)
    if (n < 2L) next

    # Time positions
    if (time %in% names(g)) {
      t_pos <- as.numeric(g[[time]])
      if (!is.null(begin_state)) t_pos <- c(min(t_pos, na.rm = TRUE) - 1, t_pos)
      if (!is.null(end_state)) t_pos <- c(t_pos, max(t_pos, na.rm = TRUE) + 1)
    } else {
      t_pos <- seq_len(n)
    }

    a_int <- match(a, all_states)

    for (i in seq_len(n)) {
      if (is.na(a_int[i])) next
      for (j in seq_len(n)) {
        if (i == j) next
        if (is.na(a_int[j])) next
        if (direction == "forward" && i >= j) next
        if (direction == "backward" && i <= j) next # nocov

        d <- decay(t_pos[i], t_pos[j], lambda)
        idx <- (a_int[i] - 1L) * n_states + a_int[j]
        counts[idx] <- counts[idx] + d
      }
    }
  }

  matrix(
    counts,
    nrow = n_states,
    ncol = n_states,
    byrow = TRUE,
    dimnames = list(all_states, all_states)
  )
}


#' Attention (decay-weighted) estimator
#'
#' Computes decay-weighted attention transitions across all position pairs.
#' For each pair of positions (i, j) in a sequence, computes a weight based
#' on temporal distance and accumulates weighted counts per state pair.
#'
#' @param data Data frame of sequence data.
#' @param format Character. "auto", "wide", or "long".
#' @param action Character. Action column name (long format).
#' @param id Character vector. ID column(s).
#' @param time Character. Time column name (long format).
#' @param cols Character vector. State columns (wide format).
#' @param lambda Numeric. Decay rate parameter. Higher = slower decay.
#'   Default: 1.
#' @param direction Character. Which position pairs to consider:
#'   "forward" (i < j), "backward" (i > j), "both". Default: "forward".
#' @param decay Function or NULL. Custom decay function(ti, tj, lambda).
#'   Default: exp(-abs(ti - tj) / lambda).
#' @param time_matrix Matrix or NULL. Custom time values per cell.
#' @param duration Numeric vector or NULL. Per-column durations.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with matrix, nodes, directed, cleaned_data.
#' @noRd
.estimator_attention <- function(data,
                                  format = "auto",
                                  action = "Action",
                                  id = NULL,
                                  time = "Time",
                                  cols = NULL,
                                  alphabet = NULL,
                                  concat = 1L,
                                  begin_state = NULL,
                                  end_state = NULL,
                                  lambda = 1,
                                  direction = "forward",
                                  decay = NULL,
                                  time_matrix = NULL,
                                  duration = NULL,
                                  ...) {
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(lambda), length(lambda) == 1, lambda > 0)
  direction <- match.arg(direction, c("forward", "backward", "both"))

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide" # nocov
  }

  if (format == "wide") {
    attn_mat <- .count_attention_wide(
      data, id = id, cols = cols, lambda = lambda,
      direction = direction, decay = decay,
      time_matrix = time_matrix, duration = duration, alphabet = alphabet,
      concat = concat, begin_state = begin_state, end_state = end_state
    )
  } else {
    attn_mat <- .count_attention_long(
      data, action = action, id = id, time = time,
      lambda = lambda, direction = direction, decay = decay,
      time_matrix = time_matrix, duration = duration, alphabet = alphabet,
      begin_state = begin_state, end_state = end_state
    )
  }

  states <- rownames(attn_mat)
  initial <- .compute_initial_probs(
    data, states, format = format,
    action = action, id = id, time = time, cols = cols, concat = concat,
    begin_state = begin_state, end_state = end_state
  )
  list(
    matrix = attn_mat,
    nodes = states,
    directed = TRUE,
    cleaned_data = data,
    initial = initial
  )
}


# ---- Association estimators ----

#' Prepare association input: clean data frame or validate matrix
#'
#' Handles data frame cleaning (drop NA, zero-variance, non-syntactic cols)
#' and matrix input (symmetric check, cor/cov detection).
#'
#' @return List with \code{S} (correlation matrix), \code{n} (sample size).
#' @noRd
.prepare_association_input <- function(data, id_col = NULL, n = NULL,
                                       cor_method = "pearson",
                                       input_type = "auto") {
  if (is.data.frame(data)) {
    # Exclude id columns, "rid", and non-numeric columns
    exclude <- c(id_col, "rid")
    numeric_cols <- vapply(data, is.numeric, logical(1))
    keep <- setdiff(names(data)[numeric_cols], exclude)

    # Drop columns with non-syntactic names
    syntactic <- make.names(keep) == keep
    if (any(!syntactic)) {
      dropped <- keep[!syntactic]
      message("Dropping non-syntactic columns: ",
              paste(dropped, collapse = ", "))
      keep <- keep[syntactic]
    }

    if (length(keep) < 2) {
      stop("At least 2 numeric columns are required after cleaning.")
    }

    mat <- as.matrix(data[, keep, drop = FALSE])

    # Drop all-NA columns
    all_na <- apply(mat, 2, function(x) all(is.na(x)))
    if (any(all_na)) {
      message("Dropping all-NA columns: ",
              paste(colnames(mat)[all_na], collapse = ", "))
      mat <- mat[, !all_na, drop = FALSE]
    }

    # Drop rows with any NA
    complete <- complete.cases(mat)
    if (!all(complete)) {
      n_dropped <- sum(!complete)
      message("Dropping ", n_dropped, " rows with NA values.")
      mat <- mat[complete, , drop = FALSE]
    }

    if (nrow(mat) < 3) {
      stop("Fewer than 3 complete rows remain after removing NAs.")
    }

    # Drop zero-variance columns
    col_vars <- apply(mat, 2, stats::var)
    zero_var <- colnames(mat)[col_vars == 0]
    if (length(zero_var) > 0) {
      message("Dropping zero-variance columns: ",
              paste(zero_var, collapse = ", "))
      mat <- mat[, col_vars > 0, drop = FALSE]
    }

    if (ncol(mat) < 2) {
      stop("At least 2 variable columns are required after cleaning.")
    }

    n <- nrow(mat)
    S <- cor(mat, method = cor_method)

  } else if (is.matrix(data) && is.numeric(data)) {
    # Non-square numeric matrix is raw data (rows = observations,
    # cols = variables) — treat it like a data.frame. This path is hit by
    # resampling callers (centrality_stability, bootstrap_network) that
    # re-invoke estimators on row-subsetted copies of netobject$data.
    if (nrow(data) != ncol(data)) {
      return(Recall(as.data.frame(data), id_col = id_col, n = n,
                    cor_method = cor_method, input_type = input_type))
    }
    if (!isSymmetric(unname(data), tol = 1e-8)) {
      stop("Matrix input must be symmetric.")
    }
    if (is.null(n)) {
      stop("Sample size 'n' is required when data is a matrix.")
    }
    stopifnot(is.numeric(n), length(n) == 1, n > 0)

    if (input_type == "auto") {
      diag_vals <- diag(data)
      input_type <- if (all(abs(diag_vals - 1) < 1e-8)) "cor" else "cov"
    }

    if (input_type == "cov") {
      d <- sqrt(diag(data))
      S <- data / outer(d, d)
    } else {
      S <- data
    }

    if (is.null(colnames(S))) {
      colnames(S) <- rownames(S) <- paste0("V", seq_len(ncol(S)))
    }
    mat <- NULL
  } else {
    stop("data must be a data frame or a square symmetric matrix.")
  }

  list(S = S, n = n, mat = mat)
}


#' Correlation estimator
#' @noRd
.estimator_cor <- function(data,
                           id_col = NULL,
                           n = NULL,
                           cor_method = "pearson",
                           input_type = "auto",
                           threshold = 0,
                           ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n

  net <- S
  diag(net) <- 0
  net[abs(net) < threshold] <- 0

  nodes <- colnames(net)
  list(
    matrix = net,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    cor_matrix = S,
    n = n_obs,
    p = ncol(S)
  )
}


#' Partial correlation estimator (unregularized)
#' @noRd
.estimator_pcor <- function(data,
                            id_col = NULL,
                            n = NULL,
                            cor_method = "pearson",
                            input_type = "auto",
                            threshold = 0,
                            ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n

  rc <- rcond(S)
  if (rc < .Machine$double.eps) {
    stop(
      "Correlation matrix is singular (p >= n or collinear variables). ",
      "Use method = 'glasso' for regularised estimation.",
      call. = FALSE
    )
  }
  if (rc < 1e-12) {
    warning(sprintf(
      "Correlation matrix is near-singular (rcond = %.2e). ",
      rc
    ), "Results may be numerically unstable. Consider method = 'glasso'.",
    call. = FALSE)
  }

  Wi <- solve(S)
  colnames(Wi) <- rownames(Wi) <- colnames(S)

  pcor <- .precision_to_pcor(Wi, threshold)
  colnames(pcor) <- rownames(pcor) <- colnames(S)

  nodes <- colnames(pcor)
  list(
    matrix = pcor,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    precision_matrix = Wi,
    cor_matrix = S,
    n = n_obs,
    p = ncol(S)
  )
}


# ---- Shared association helpers ----

#' Compute log-spaced lambda path
#' @noRd
.compute_lambda_path <- function(S, nlambda, lambda.min.ratio) {
  lambda_max <- max(abs(S[upper.tri(S)]))
  if (lambda_max <= 0) {
    stop("All off-diagonal correlations are zero; nothing to regularise.")
  }
  lambda_min <- lambda_max * lambda.min.ratio
  exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))
}


#' Select best lambda via EBIC using glasso fits with warm starts
#' @noRd
.select_ebic <- function(S, lambda_path, n, gamma, penalize_diagonal) {
  p <- ncol(S)
  n_lambda <- length(lambda_path)
  ebic_vals <- numeric(n_lambda)

  w_prev <- NULL
  wi_prev <- NULL
  best_idx <- 1L
  best_ebic <- Inf
  best_wi <- NULL

  for (i in seq_along(lambda_path)) {
    lam <- lambda_path[i]

    fit <- tryCatch(
      glasso::glasso(
        s = S,
        rho = lam,
        penalize.diagonal = penalize_diagonal,
        start = if (is.null(w_prev)) "cold" else "warm",
        w.init = w_prev,
        wi.init = wi_prev,
        trace = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      ebic_vals[i] <- Inf # nocov start
      next # nocov end
    }

    w_prev <- fit$w
    wi_prev <- fit$wi

    log_det <- determinant(fit$wi, logarithm = TRUE)
    if (log_det$sign <= 0) {
      ebic_vals[i] <- Inf # nocov start
      next # nocov end
    }
    log_det_val <- as.numeric(log_det$modulus)

    loglik <- (n / 2) * (log_det_val - sum(diag(S %*% fit$wi)))
    npar <- sum(abs(fit$wi[upper.tri(fit$wi)]) > 1e-10)
    ebic_vals[i] <- -2 * loglik + npar * log(n) +
      4 * npar * gamma * log(p)

    if (ebic_vals[i] < best_ebic) {
      best_ebic <- ebic_vals[i]
      best_idx <- i
      best_wi <- fit$wi
    }
  }

  if (is.null(best_wi)) {
    stop("All glasso fits failed. Check your input data.") # nocov
  }

  colnames(best_wi) <- rownames(best_wi) <- colnames(S)

  list(
    wi        = best_wi,
    lambda    = lambda_path[best_idx],
    ebic      = best_ebic,
    ebic_path = ebic_vals
  )
}


#' Convert precision matrix to partial correlations (qgraph-compatible)
#' Uses cov2cor for numerical stability, matching qgraph::wi2net.
#' @noRd
.wi2net <- function(x) {
  x <- -stats::cov2cor(x)
  diag(x) <- 0
  # forceSymmetric: copy upper triangle to lower (matches qgraph)
  x[lower.tri(x)] <- t(x)[lower.tri(x)]
  x
}

.precision_to_pcor <- function(Wi, threshold) {
  d <- sqrt(diag(Wi))
  pcor <- -Wi / outer(d, d)
  diag(pcor) <- 0
  pcor <- (pcor + t(pcor)) / 2  # symmetrize (glasso convergence can leave tiny asymmetry)
  pcor[abs(pcor) < threshold] <- 0
  pcor
}


#' EBICglasso estimator
#' @noRd
.estimator_glasso <- function(data,
                              id_col = NULL,
                              n = NULL,
                              gamma = 0.5,
                              nlambda = 100L,
                              lambda.min.ratio = 0.01,
                              penalize.diagonal = FALSE,
                              refit = FALSE,
                              cor_method = "pearson",
                              input_type = "auto",
                              threshold = 0,
                              ...) {
  prepared <- .prepare_association_input(
    data, id_col = id_col, n = n,
    cor_method = cor_method, input_type = input_type
  )
  S <- prepared$S
  n_obs <- prepared$n
  p <- ncol(S)

  stopifnot(is.numeric(gamma), length(gamma) == 1, gamma >= 0)
  stopifnot(is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2)
  stopifnot(is.numeric(lambda.min.ratio), lambda.min.ratio > 0,
            lambda.min.ratio < 1)
  stopifnot(is.logical(penalize.diagonal), length(penalize.diagonal) == 1)

  lambda_path <- .compute_lambda_path(S, nlambda, lambda.min.ratio)
  selected <- .select_ebic(S, lambda_path, n_obs, gamma, penalize.diagonal)

  wi <- selected$wi

  if (isTRUE(refit)) {
    # Refit with zero-constrained unregularized glasso for unbiased estimates
    net_pattern <- -.wi2net(wi)
    zero_idx <- which(net_pattern == 0 & upper.tri(net_pattern), arr.ind = TRUE)
    if (nrow(zero_idx) > 0L) {
      refit_result <- suppressWarnings(glasso::glasso(
        S, rho = 0, zero = zero_idx, trace = 0,
        penalize.diagonal = penalize.diagonal))
    } else {
      refit_result <- suppressWarnings(glasso::glasso( # nocov start
        S, rho = 0, trace = 0,
        penalize.diagonal = penalize.diagonal)) # nocov end
    }
    wi <- refit_result$wi
  }

  pcor <- .wi2net(wi)
  pcor[abs(pcor) < threshold] <- 0
  colnames(pcor) <- rownames(pcor) <- colnames(S)

  nodes <- colnames(pcor)
  list(
    matrix = pcor,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = prepared$mat,
    precision_matrix = wi,
    cor_matrix = S,
    lambda_selected = selected$lambda,
    ebic_selected = selected$ebic,
    lambda_path = lambda_path,
    ebic_path = selected$ebic_path,
    gamma = gamma,
    n = n_obs,
    p = p
  )
}


# ---- Ising model estimator ----

#' Prepare input for Ising model estimation
#'
#' Validates and cleans a data frame for Ising model estimation.
#' Drops non-numeric, ID, non-syntactic, zero-variance, and all-NA columns.
#' Validates that remaining values are all binary (0 or 1).
#'
#' @param data Data frame of binary (0/1) variables.
#' @param id_col Character or NULL. ID column(s) to exclude.
#'
#' @return A list with:
#'   \describe{
#'     \item{mat}{Numeric matrix of binary values (complete cases).}
#'     \item{n}{Number of observations (rows).}
#'     \item{p}{Number of variables (columns).}
#'     \item{nodes}{Character vector of variable names.}
#'   }
#' @noRd
.prepare_ising_input <- function(data, id_col = NULL) {
  stopifnot(is.data.frame(data))

  # Exclude id columns and "rid"

  exclude <- c(id_col, "rid")
  numeric_cols <- vapply(data, is.numeric, logical(1))
  keep <- setdiff(names(data)[numeric_cols], exclude)

  # Drop columns with non-syntactic names
  syntactic <- make.names(keep) == keep
  if (any(!syntactic)) {
    dropped <- keep[!syntactic]
    message("Dropping non-syntactic columns: ",
            paste(dropped, collapse = ", "))
    keep <- keep[syntactic]
  }

  if (length(keep) < 2L) {
    stop("At least 2 numeric columns are required after cleaning.")
  }

  mat <- as.matrix(data[, keep, drop = FALSE])

  # Drop all-NA columns
  all_na <- apply(mat, 2, function(x) all(is.na(x)))
  if (any(all_na)) {
    message("Dropping all-NA columns: ",
            paste(colnames(mat)[all_na], collapse = ", "))
    mat <- mat[, !all_na, drop = FALSE]
  }

  # Drop rows with any NA
  complete <- complete.cases(mat)
  if (!all(complete)) {
    n_dropped <- sum(!complete)
    message("Dropping ", n_dropped, " rows with NA values.")
    mat <- mat[complete, , drop = FALSE]
  }

  if (nrow(mat) < 3L) {
    stop("Fewer than 3 complete rows remain after removing NAs.")
  }

  # Drop zero-variance columns
  col_vars <- apply(mat, 2, stats::var)
  zero_var <- colnames(mat)[col_vars == 0]
  if (length(zero_var) > 0L) {
    message("Dropping zero-variance columns: ",
            paste(zero_var, collapse = ", "))
    mat <- mat[, col_vars > 0, drop = FALSE]
  }

  if (ncol(mat) < 2L) {
    stop("At least 2 variable columns are required after cleaning.")
  }

  # Validate binary: all values must be 0 or 1
  unique_vals <- unique(as.vector(mat))
  if (!all(unique_vals %in% c(0, 1))) {
    bad <- setdiff(unique_vals, c(0, 1))
    stop("Ising model requires binary (0/1) data. Found non-binary values: ",
         paste(head(bad, 5), collapse = ", "))
  }

  list(
    mat = mat,
    n = nrow(mat),
    p = ncol(mat),
    nodes = colnames(mat)
  )
}


#' Numerically stable log(1 + exp(x))
#'
#' Avoids overflow for large x and precision loss for small x.
#'
#' @param x Numeric vector.
#' @return Numeric vector of log(1 + exp(x)).
#' @noRd
.log1pexp <- function(x) {
  out <- numeric(length(x))
  big <- x > 20
  small <- x < -20
  mid <- !big & !small
  out[big] <- x[big]
  out[small] <- exp(x[small])
  out[mid] <- log1p(exp(x[mid]))
  out
}


#' Nodewise L1-penalized logistic regression with EBIC selection
#'
#' Core algorithm for Ising model estimation (IsingFit approach).
#' For each node j, fits L1-penalized logistic regression of \code{X[,j]} on \code{X[,-j]}
#' using glmnet, then selects lambda via EBIC.
#'
#' @param mat Numeric matrix of binary (0/1) values (n x p).
#' @param gamma Numeric. EBIC hyperparameter (0 = BIC, higher = sparser).
#' @param nlambda Integer. Number of lambda values in the regularization path.
#'
#' @return A list with:
#'   \describe{
#'     \item{coef_matrix}{p x p asymmetric coefficient matrix (row j = regression
#'       of node j on others).}
#'     \item{thresholds}{Numeric vector of intercepts (length p).}
#'     \item{lambda_selected}{Numeric vector of selected lambda per node.}
#'   }
#' @noRd
.ising_nodewise_ebic <- function(mat, gamma = 0.25, nlambda = 100L) {
  n <- nrow(mat)
  p <- ncol(mat)
  n_predictors <- p - 1L
  node_names <- colnames(mat)

  coef_matrix <- matrix(0, nrow = p, ncol = p,
                         dimnames = list(node_names, node_names))
  thresholds <- numeric(p)
  names(thresholds) <- node_names
  lambda_selected <- numeric(p)
  names(lambda_selected) <- node_names

  for (j in seq_len(p)) {
    y <- mat[, j]
    X <- mat[, -j, drop = FALSE]

    # Fit L1-penalized logistic regression
    fit <- glmnet::glmnet(X, y, family = "binomial", nlambda = nlambda)

    # Compute EBIC for each lambda in the path
    # Nodewise EBIC (IsingFit): -2*loglik + k*log(n) + 2*gamma*k*log(p-1)
    n_lam <- length(fit$lambda)
    ebic_vals <- numeric(n_lam)

    for (k in seq_len(n_lam)) {
      beta_k <- fit$beta[, k]
      intercept_k <- fit$a0[k]

      # Linear predictor
      eta <- as.vector(X %*% beta_k) + intercept_k

      # Log-likelihood: sum(y * eta - log(1 + exp(eta)))
      loglik <- sum(y * eta - .log1pexp(eta))

      # Number of nonzero coefficients (excluding intercept)
      n_edges <- sum(abs(beta_k) > 0)

      # Nodewise EBIC = -2*loglik + k*log(n) + 2*gamma*k*log(p-1)
      ebic_vals[k] <- -2 * loglik + n_edges * log(n) +
        2 * n_edges * gamma * log(n_predictors)
    }

    # Select lambda minimizing EBIC
    best_idx <- which.min(ebic_vals)
    best_beta <- fit$beta[, best_idx]
    best_intercept <- fit$a0[best_idx]

    # Place coefficients in the row for node j
    other_idx <- seq_len(p)[-j]
    coef_matrix[j, other_idx] <- as.vector(best_beta)
    thresholds[j] <- best_intercept
    lambda_selected[j] <- fit$lambda[best_idx]
  }

  list(
    coef_matrix = coef_matrix,
    thresholds = thresholds,
    lambda_selected = lambda_selected
  )
}


#' Symmetrize asymmetric Ising coefficient matrix
#'
#' @param coef_matrix p x p asymmetric coefficient matrix from nodewise
#'   regression.
#' @param rule Character. Symmetrization rule: \code{"AND"} (default) or
#'   \code{"OR"}.
#'
#' @return Symmetric p x p weight matrix with zero diagonal.
#' @noRd
.symmetrize_ising <- function(coef_matrix, rule = "AND") {
  p <- nrow(coef_matrix)
  sym <- matrix(0, nrow = p, ncol = p,
                dimnames = dimnames(coef_matrix))

  if (rule == "AND") {
    # Edge only if BOTH directions nonzero; weight = average
    both_nonzero <- (coef_matrix != 0) & (t(coef_matrix) != 0)
    sym[both_nonzero] <- (coef_matrix[both_nonzero] +
                            t(coef_matrix)[both_nonzero]) / 2
  } else if (rule == "OR") {
    # Simple average of both directions (matches IsingFit)
    sym <- (coef_matrix + t(coef_matrix)) / 2
  }

  diag(sym) <- 0
  sym
}


#' Ising Model Network Estimator
#'
#' Estimates an Ising model network using nodewise L1-penalized logistic
#' regression with EBIC model selection (van Borkulo et al., 2014). Produces
#' an undirected weighted network of conditional dependencies between binary
#' (0/1) variables.
#'
#' @param data Data frame of binary (0/1) variables.
#' @param id_col Character or NULL. ID column(s) to exclude.
#' @param gamma Numeric. EBIC hyperparameter. Higher values produce sparser
#'   networks. Default: 0.25 (IsingFit convention).
#' @param nlambda Integer. Number of lambda values for the regularization path.
#'   Default: 100.
#' @param rule Character. Symmetrization rule: \code{"AND"} (conservative,
#'   edge only if both directions nonzero) or \code{"OR"} (liberal, edge if
#'   either direction nonzero). Default: \code{"AND"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with:
#'   \describe{
#'     \item{matrix}{Symmetric weighted adjacency matrix.}
#'     \item{nodes}{Character vector of variable names.}
#'     \item{directed}{Logical: always FALSE.}
#'     \item{cleaned_data}{Cleaned binary data matrix.}
#'     \item{thresholds}{Numeric vector of node thresholds (intercepts).}
#'     \item{asymm_weights}{Asymmetric coefficient matrix before symmetrization.}
#'     \item{rule}{Symmetrization rule used.}
#'     \item{gamma}{EBIC hyperparameter used.}
#'     \item{n}{Sample size.}
#'     \item{p}{Number of variables.}
#'     \item{lambda_selected}{Per-node selected lambda values.}
#'   }
#'
#' @references
#' van Borkulo, C. D., Borsboom, D., Epskamp, S., Blanken, T. F.,
#' Boschloo, L., Schoevers, R. A., & Waldorp, L. J. (2014). A new method
#' for constructing networks from binary data. \emph{Scientific Reports},
#' 4, 5918. \doi{10.1038/srep05918}
#'
#' @noRd
.estimator_ising <- function(data,
                              id_col = NULL,
                              gamma = 0.25,
                              nlambda = 100L,
                              rule = "AND",
                              ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop( # nocov start
      "Package 'glmnet' is required for Ising model estimation. ",
      "Install it with: install.packages('glmnet')",
      call. = FALSE
    ) # nocov end
  }

  # Validate parameters
  stopifnot(is.numeric(gamma), length(gamma) == 1, gamma >= 0)
  nlambda <- as.integer(nlambda)
  stopifnot(is.integer(nlambda), length(nlambda) == 1, nlambda >= 2L)
  rule <- match.arg(rule, c("AND", "OR"))

  # Prepare input
  prepared <- .prepare_ising_input(data, id_col = id_col)
  mat <- prepared$mat
  n <- prepared$n
  p <- prepared$p
  nodes <- prepared$nodes

  # Run nodewise logistic regression with EBIC
  nodewise <- .ising_nodewise_ebic(mat, gamma = gamma, nlambda = nlambda)

  # Symmetrize
  sym_matrix <- .symmetrize_ising(nodewise$coef_matrix, rule = rule)

  list(
    matrix = sym_matrix,
    nodes = nodes,
    directed = FALSE,
    cleaned_data = mat,
    thresholds = nodewise$thresholds,
    asymm_weights = nodewise$coef_matrix,
    rule = rule,
    gamma = gamma,
    n = n,
    p = p,
    lambda_selected = nodewise$lambda_selected
  )
}
