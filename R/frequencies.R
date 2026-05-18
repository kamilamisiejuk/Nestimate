#' Sequence Data Conversion Functions
#'
#' @description
#' Functions for converting sequence data (long or wide format) into
#' transition frequency matrices and other useful representations.
#'
#' @name frequencies
#' @keywords internal
NULL


#' Build a Transition Frequency Matrix
#'
#' @description
#' Convert long or wide format sequence data into a transition frequency matrix.
#' Counts how many times each transition from state_i to state_j occurs across
#' all sequences.
#'
#' @param data Data frame containing sequence data in long or wide format.
#' @param action Character. Name of the column containing actions/states
#'   (for long format). Default: "Action".
#' @param id Character vector. Name(s) of the column(s) identifying sequences.
#'   For long format, each unique combination of ID values defines a sequence.
#'   For wide format, used to exclude non-state columns. Default: NULL.
#' @param time Character. Name of the time column used to order actions within
#'   sequences (for long format). Default: "Time".
#' @param cols Character vector. Names of columns containing states (for wide
#'   format). If NULL, all non-ID columns are used. Default: NULL.
#' @param format Character. Format of input data: "auto" (detect automatically),
#'   "long", or "wide". Default: "auto".
#'
#' @return A square integer matrix of transition frequencies where
#'   \code{mat[i, j]} is the number of times state i was followed by state j.
#'   Row and column names are the sorted unique states. Can be passed directly
#'   to \code{tna::tna()}.
#'
#' @details
#' For \strong{long format} data, each row is a single action/event. Sequences
#' are defined by the \code{id} column(s), and actions are ordered by the
#' \code{time} column within each sequence. Consecutive actions within a
#' sequence form transition pairs.
#'
#' For \strong{wide format} data, each row is a sequence and columns represent
#' consecutive time points. Transitions are counted across consecutive columns,
#' skipping any \code{NA} values.
#'
#' @examples
#' # Wide format
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' freq <- frequencies(seqs, format = "wide")
#'
#' # Long format
#' long <- data.frame(
#'   Actor = rep(1:2, each = 3), Time = rep(1:3, 2),
#'   Action = c("A","B","C","B","A","C")
#' )
#' freq <- frequencies(long, action = "Action", id = "Actor")
#'
#' @seealso \code{\link{convert_sequence_format}} for converting to other
#'   representations (frequency counts, one-hot, edge lists).
#'
#' @export
frequencies <- function(data,
                        action = "Action",
                        id = NULL,
                        time = "Time",
                        cols = NULL,
                        format = c("auto", "long", "wide")) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(action), length(action) == 1)
  stopifnot(is.character(time), length(time) == 1)
  stopifnot(is.null(id) || is.character(id))
  stopifnot(is.null(cols) || is.character(cols))
  format <- match.arg(format)

  mat <- .count_transitions(
    data, format = format, action = action, id = id, time = time, cols = cols
  )
  if (is.matrix(mat)) {
    class(mat) <- c("nest_transition_counts", class(mat))
  }
  mat
}


#' Summary Method for Transition Count Matrices
#'
#' @param object A \code{nest_transition_counts} matrix returned by
#'   \code{\link{frequencies}()}.
#' @param ... Additional arguments (ignored).
#' @return A tidy data frame with columns \code{from}, \code{to},
#'   \code{count}, with one row per non-zero transition.
#' @export
summary.nest_transition_counts <- function(object, ...) {
  .matrix_to_long_df(unclass(object), value_col = "count",
                     include = "positive")
}


#' Convert Sequence Data to Different Formats
#'
#' @description
#' Convert wide or long sequence data into frequency counts, one-hot encoding,
#' edge lists, or follows format.
#'
#' @param data Data frame containing sequence data.
#' @param seq_cols Character vector. Names of columns containing sequential
#'   states (for wide format input). If NULL, all columns except \code{id_col}
#'   are used. Default: NULL.
#' @param id_col Character vector. Name(s) of the ID column(s). For long
#'   format, required. For wide format, optional: if NULL, the first column
#'   is used as the id only when it is a genuine identifier (its values are
#'   disjoint from the states in the remaining columns); for canonical wide
#'   sequence data with no id column (e.g. \code{V1..Vn} / \code{T1..Tn}), a
#'   row-index id is synthesized and every column is treated as a sequence
#'   column. Default: NULL.
#' @param action Character or NULL. Name of the column containing actions/states
#'   (for long format input). If provided, data is treated as long format.
#'   Default: NULL.
#' @param time Character or NULL. Name of the time column for ordering actions
#'   within sequences (for long format). Default: NULL.
#' @param format Character. Output format:
#' \describe{
#'   \item{"frequency"}{Count of each action per sequence (wide, one column per
#'     state).}
#'   \item{"onehot"}{Binary presence/absence of each action per sequence.}
#'   \item{"edgelist"}{Consecutive transition pairs (from, to) per sequence.}
#'   \item{"follows"}{Each action paired with the action that preceded it.}
#' }
#'
#' @return A data frame in the requested format:
#' \describe{
#'   \item{frequency}{ID columns + one integer column per state with counts.}
#'   \item{onehot}{ID columns + one binary column per state (0/1).}
#'   \item{edgelist}{ID columns + \code{from} and \code{to} columns.}
#'   \item{follows}{ID columns + \code{act} and \code{follows} columns.}
#' }
#'
#' @examples
#' # Wide format input
#' seqs <- data.frame(V1 = c("A","B","A"), V2 = c("B","A","C"), V3 = c("A","C","B"))
#' convert_sequence_format(seqs, format = "frequency")
#' convert_sequence_format(seqs, format = "edgelist")
#'
#' @seealso \code{\link{frequencies}} for building transition frequency matrices.
#'
#' @export
convert_sequence_format <- function(data,
                                    seq_cols = NULL,
                                    id_col = NULL,
                                    action = NULL,
                                    time = NULL,
                                    format = c("frequency", "onehot",
                                               "edgelist", "follows")) {
  stopifnot(is.data.frame(data))
  format <- match.arg(format)

  is_long <- !is.null(action) && action %in% names(data)

  if (is_long) {
    long <- .standardize_long(data, id_col, action, time)
    id_col <- long$id_col
    long <- long$data
  } else {
    if (is.null(id_col)) {
      # No id column declared. Never infer one from the data: value-based id
      # detection silently drops the first time point whenever a start-only
      # state does not recur later (audit A11-F01; Codex review). Synthesize a
      # row-index id (the same safe pattern as wide_to_long()) and treat every
      # column as a state, so no event is ever dropped. Callers with a genuine
      # leading id column must pass `id_col` explicitly.
      data[[".rid"]] <- seq_len(nrow(data))
      id_col <- ".rid"
    }
    if (is.null(seq_cols)) seq_cols <- setdiff(names(data), id_col)
    .validate_cols(data, c(id_col, seq_cols))
    long <- .pivot_wide_to_long(data, id_col, seq_cols)
  }

  switch(format,
    frequency = .fmt_frequency(long, id_col),
    onehot    = .fmt_onehot(long, id_col),
    edgelist  = .fmt_edgelist(long, id_col),
    follows   = .fmt_follows(long, id_col)
  )
}


# ---- Internal helpers for frequencies() ----

#' @noRd
.build_pairs_long <- function(data, action, id, time) {
  if (!action %in% names(data)) {
    stop("Action column '", action, "' not found in data.")
  }
  if (!is.null(id)) {
    missing_ids <- setdiff(id, names(data))
    if (length(missing_ids) > 0) {
      stop("ID column(s) not found: ", paste(missing_ids, collapse = ", "))
    }
  }

  if (is.null(id)) {
    groups <- list(data)
  } else if (length(id) == 1L) {
    groups <- split(data, data[[id]], drop = TRUE)
  } else {
    group_key <- interaction(data[, id, drop = FALSE], drop = TRUE)
    groups <- split(data, group_key, drop = TRUE)
  }

  pairs <- lapply(groups, function(g) {
    if (time %in% names(g)) {
      g <- g[order(g[[time]]), ]
    }
    actions <- g[[action]]
    n <- length(actions)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    from <- actions[-n]
    to <- actions[-1L]
    valid <- !is.na(from) & !is.na(to)
    data.frame(from = from[valid], to = to[valid], stringsAsFactors = FALSE)
  })

  do.call(rbind, pairs)
}

#' @noRd
.build_pairs_wide <- function(data, id, cols) {
  if (!is.null(cols)) {
    state_cols <- cols
  } else if (!is.null(id)) {
    state_cols <- setdiff(names(data), id)
  } else {
    state_cols <- names(data)
  }

  missing_cols <- setdiff(state_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columns not found: ", paste(missing_cols, collapse = ", "))
  }
  if (length(state_cols) < 2L) {
    stop("At least 2 state columns are required for wide format.")
  }

  state_data <- data[, state_cols, drop = FALSE]

  pairs <- lapply(seq_len(nrow(state_data)), function(i) {
    row_vals <- unlist(state_data[i, ], use.names = FALSE)
    row_vals <- row_vals[!is.na(row_vals)]
    n <- length(row_vals)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    data.frame(from = row_vals[-n], to = row_vals[-1L],
               stringsAsFactors = FALSE)
  })

  do.call(rbind, pairs)
}


# ---- Internal helpers for convert_sequence_format() ----

#' @noRd
.validate_cols <- function(data, cols) {
  missing <- setdiff(cols, names(data))
  if (length(missing) > 0) {
    stop("Columns not found: ", paste(missing, collapse = ", "))
  }
}

#' Standardize long format input for convert_sequence_format
#' @noRd
.standardize_long <- function(data, id_col, action, time) {
  if (is.null(id_col)) { # nocov start
    non_meta <- setdiff(names(data), c(action, time))
    if (length(non_meta) == 0) {
      stop("id_col is required for long format data.")
    }
    id_col <- non_meta[1]
    message("Using '", id_col, "' as id_col.")
  } # nocov end
  .validate_cols(data, c(id_col, action))

  # Order by id columns + time
  order_cols <- id_col
  if (!is.null(time) && time %in% names(data)) {
    order_cols <- c(order_cols, time)
  }
  data <- data[do.call(order, data[, order_cols, drop = FALSE]), ]

  # Create sequence identifier
  if (length(id_col) == 1L) {
    rid_vals <- data[[id_col]]
  } else {
    rid_vals <- interaction(data[, id_col, drop = FALSE], drop = TRUE)
  }

  long <- data[, id_col, drop = FALSE]
  long$rid <- rid_vals
  long$act <- data[[action]]
  long <- long[!is.na(long$act) & long$act != "", ]

  list(data = long, id_col = id_col)
}

#' Pivot wide data to standardized long format
#' @noRd
.pivot_wide_to_long <- function(data, id_col, seq_cols) {
  rid_vals <- seq_len(nrow(data))
  parts <- lapply(seq_cols, function(col) {
    df <- data[, id_col, drop = FALSE]
    df$rid <- rid_vals
    df$act <- data[[col]]
    df
  })
  long <- do.call(rbind, parts)
  long <- long[!is.na(long$act) & long$act != "", ]
  as.data.frame(long)
}

#' @noRd
.fmt_frequency <- function(long, id_col) {
  tab <- table(long$rid, long$act)
  freq_df <- as.data.frame.matrix(tab)
  freq_df$rid <- as.integer(rownames(freq_df))
  id_map <- unique(long[, c(id_col, "rid"), drop = FALSE])
  result <- merge(id_map, freq_df, by = "rid")
  as.data.frame(result)
}

#' @noRd
.fmt_onehot <- function(long, id_col) {
  uniq <- unique(long[, c(id_col, "rid", "act"), drop = FALSE])
  all_acts <- sort(unique(uniq$act))
  id_map <- unique(uniq[, c(id_col, "rid"), drop = FALSE])
  mat <- matrix(0L, nrow = nrow(id_map), ncol = length(all_acts),
                dimnames = list(NULL, all_acts))
  for (i in seq_len(nrow(id_map))) {
    acts_i <- uniq$act[uniq$rid == id_map$rid[i]]
    mat[i, acts_i] <- 1L
  }
  result <- cbind(id_map, as.data.frame(mat))
  as.data.frame(result)
}

#' @noRd
.fmt_edgelist <- function(long, id_col) {
  groups <- split(long, long$rid)
  parts <- lapply(groups, function(g) {
    n <- nrow(g)
    if (n < 2L) {
      return(data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE))
    }
    edge <- g[seq_len(n - 1L), id_col, drop = FALSE]
    edge$from <- g$act[-n]
    edge$to <- g$act[-1L]
    edge
  })
  result <- do.call(rbind, parts)
  result <- result[, c(id_col, "from", "to"), drop = FALSE]
  rownames(result) <- NULL
  as.data.frame(result)
}

#' @noRd
.fmt_follows <- function(long, id_col) {
  groups <- split(long, long$rid)
  parts <- lapply(groups, function(g) {
    n <- nrow(g)
    if (n < 2L) {
      return(data.frame(act = character(0), follows = character(0),
                        stringsAsFactors = FALSE))
    }
    fol <- g[2:n, c(id_col), drop = FALSE]
    fol$act <- g$act[-1L]
    fol$follows <- g$act[-n]
    fol
  })
  result <- do.call(rbind, parts)
  result <- result[, c(id_col, "act", "follows"), drop = FALSE]
  rownames(result) <- NULL
  as.data.frame(result)
}
