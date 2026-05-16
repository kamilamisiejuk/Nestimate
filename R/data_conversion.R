#' Data Format Conversion Functions
#'
#' @description
#' Functions for converting between wide and long data formats commonly used
#' in Temporal Network Analysis.
#'
#' @name data_conversion
#' @keywords internal
NULL

#' Convert Wide Sequences to Long Format
#'
#' @description
#' Convert sequence data from wide format (one row per sequence, columns as
#' time points) to long format (one row per action).
#'
#' @param data Data frame in wide format with sequences in rows.
#' @param id_col Character. Name of the ID column, or NULL to auto-generate IDs.
#'   Default: NULL.
#' @param time_prefix Character. Prefix for time point columns (e.g., "V" for V1,
#'   V2, ...). Default: "V".
#' @param action_col Character. Name of the action column in output.
#'   Default: "Action".
#' @param time_col Character. Name of the time column in output.
#'   Default: "Time".
#' @param drop_na Logical. Whether to drop NA values. Default: TRUE.
#'
#' @return A data frame in long format with columns:
#' \describe{
#'   \item{id}{Sequence identifier (integer).}
#'   \item{Time}{Time point within the sequence (integer).}
#'   \item{Action}{The action/state at that time point (character).}
#' }
#' Any additional columns from the original data are preserved.
#'
#' @details
#' This function converts data from the format produced by `simulate_sequences()`
#' to the long format used by many TNA functions and analyses.
#'
#' @examples
#' wide_data <- data.frame(
#'   V1 = c("A", "B", "C"), V2 = c("B", "C", "A"), V3 = c("C", "A", "B")
#' )
#' long_data <- wide_to_long(wide_data)
#' head(long_data)
#'
#' @seealso \code{\link{long_to_wide}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @export
wide_to_long <- function(data,
                         id_col = NULL,
                         time_prefix = "V",
                         action_col = "Action",
                         time_col = "Time",
                         drop_na = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  stopifnot(is.character(time_prefix), length(time_prefix) == 1)
  stopifnot(is.character(action_col), length(action_col) == 1)
  stopifnot(is.character(time_col), length(time_col) == 1)
  stopifnot(is.logical(drop_na), length(drop_na) == 1)

  # Identify time columns
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(data), value = TRUE)
  if (length(time_cols) == 0) {
    stop("No columns matching time prefix '", time_prefix, "' found.")
  }

  # Add ID column if needed
  if (is.null(id_col)) {
    data$id <- seq_len(nrow(data))
    id_col <- "id"
  } else if (!id_col %in% names(data)) {
    stop("ID column '", id_col, "' not found in data.")
  }

  # Identify other columns to preserve
  other_cols <- setdiff(names(data), c(time_cols, id_col))

  # Extract time indices from column names
  time_nums <- as.integer(gsub(paste0("^", time_prefix), "", time_cols))

  # Manual stack: one data frame per time column, then rbind
  parts <- lapply(seq_along(time_cols), function(i) {
    row <- data[, c(id_col, other_cols), drop = FALSE]
    row[[time_col]] <- time_nums[i]
    row[[action_col]] <- data[[time_cols[i]]]
    row
  })
  result <- do.call(rbind, parts)

  # Drop NA values if requested
  if (drop_na) {
    result <- result[!is.na(result[[action_col]]), ]
  }

  # Sort by ID and time
  result <- result[order(result[[id_col]], result[[time_col]]), ]

  # Reset row names
  rownames(result) <- NULL

  return(result)
}


#' Convert Long Format to Wide Sequences
#'
#' @description
#' Convert sequence data from long format (one row per action) to wide format
#' (one row per sequence, columns as time points).
#'
#' @param data Data frame in long format.
#' @param id_col Character. Name of the column identifying sequences.
#'   Default: "Actor".
#' @param time_col Character. Name of the column identifying time points.
#'   Default: "Time".
#' @param action_col Character. Name of the column containing actions/states.
#'   Default: "Action".
#' @param time_prefix Character. Prefix for time point columns in output.
#'   Default: "V".
#' @param fill_na Logical. Whether to fill missing time points with NA.
#'   Default: TRUE.
#'
#' @return A data frame in wide format where each row is a sequence and
#'   columns V1, V2, ... contain the actions at each time point.
#'
#' @details
#' This function converts long format data (like that from `simulate_long_data()`)
#' to the wide format expected by `tna::tna()` and related functions.
#'
#' If `time_col` contains non-integer values (e.g., timestamps), the function
#' will use the ordering within each sequence to create time indices.
#'
#' @examples
#' long_data <- data.frame(
#'   Actor = rep(1:3, each = 4),
#'   Time = rep(1:4, 3),
#'   Action = sample(c("A", "B", "C"), 12, replace = TRUE)
#' )
#' wide_data <- long_to_wide(long_data, id_col = "Actor")
#' head(wide_data)
#'
#' @seealso \code{\link{wide_to_long}} for the reverse conversion,
#'   \code{\link{prepare_for_tna}} for preparing data for TNA analysis.
#'
#' @export
long_to_wide <- function(data,
                         id_col = "Actor",
                         time_col = "Time",
                         action_col = "Action",
                         time_prefix = "V",
                         fill_na = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  stopifnot(is.character(id_col), length(id_col) == 1)
  stopifnot(is.character(time_col), length(time_col) == 1)
  stopifnot(is.character(action_col), length(action_col) == 1)
  stopifnot(is.character(time_prefix), length(time_prefix) == 1)

  # Check required columns exist
  required_cols <- c(id_col, action_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Sort by id and time (stable sort preserves original order for ties)
  data$.orig_order <- seq_len(nrow(data))
  if (time_col %in% names(data)) {
    data <- data[order(data[[id_col]], data[[time_col]], data$.orig_order), ]
  } else {
    data <- data[order(data[[id_col]], data$.orig_order), ]
  }
  data$.orig_order <- NULL

  # Compute within-group sequential index
  data[[".time_idx"]] <- ave(
    seq_len(nrow(data)), data[[id_col]], FUN = seq_along
  )

  # Create column names
  data[[".time_name"]] <- paste0(time_prefix, data[[".time_idx"]])

  # Spread manually: split by id, build one row per id
  groups <- split(data, data[[id_col]])
  wide_list <- lapply(groups, function(g) {
    row <- data.frame(g[[id_col]][1], stringsAsFactors = FALSE)
    names(row) <- id_col
    for (i in seq_len(nrow(g))) {
      row[[g[[".time_name"]][i]]] <- g[[action_col]][i]
    }
    row
  })
  # Align columns (sequences may differ in length)
  all_cols <- unique(unlist(lapply(wide_list, names)))
  wide_list <- lapply(wide_list, function(row) {
    missing <- setdiff(all_cols, names(row))
    for (m in missing) row[[m]] <- NA_character_
    row[, all_cols, drop = FALSE]
  })
  result <- do.call(rbind, wide_list)

  # Reorder columns to ensure V1, V2, V3... order
  time_cols <- grep(paste0("^", time_prefix, "[0-9]+$"), names(result), value = TRUE)
  time_nums <- as.integer(gsub(time_prefix, "", time_cols))
  time_cols <- time_cols[order(time_nums)]
  other_cols <- setdiff(names(result), time_cols)
  result <- result[, c(other_cols, time_cols), drop = FALSE]

  # Convert to data frame
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Prepare Data for TNA Analysis
#'
#' @description
#' Prepare simulated or real data for use with `tna::tna()` and related
#' functions. Handles various input formats and ensures the output is
#' compatible with TNA models.
#'
#' @param data Data frame containing sequence data.
#' @param type Character. Type of input data:
#' \describe{
#'   \item{"sequences"}{Wide format with one row per sequence (default).}
#'   \item{"long"}{Long format with one row per action.}
#'   \item{"auto"}{Automatically detect format based on column names.}
#' }
#' @param state_names Character vector. Expected state names, or NULL to
#'   extract from data. Default: NULL.
#' @param id_col Character. Name of ID column for long format data.
#'   Default: "Actor".
#' @param time_col Character. Name of time column for long format data.
#'   Default: "Time".
#' @param action_col Character. Name of action column for long format data.
#'   Default: "Action".
#' @param validate Logical. Whether to validate that all actions are in
#'   state_names. Default: TRUE.
#'
#' @return A data frame ready for use with TNA functions. For "sequences" type,
#'   returns a data frame where each row is a sequence and columns are time
#'   points (V1, V2, ...). For "long" type, converts to wide format first.
#'
#' @details
#' This function performs several preparations:
#' \enumerate{
#'   \item Converts long format to wide format if needed.
#'   \item Validates that all actions/states are recognized.
#'   \item Removes any non-sequence columns (e.g., id, metadata).
#'   \item Converts factors to characters.
#'   \item Ensures consistent column naming (V1, V2, ...).
#' }
#'
#' @examples
#' # From wide format sequences
#' sequences <- data.frame(
#'   V1 = c("A","B","C","A"), V2 = c("B","C","A","B"),
#'   V3 = c("C","A","B","C"), V4 = c("A","B","A","B")
#' )
#' tna_data <- prepare_for_tna(sequences, type = "sequences")
#'
#' @seealso \code{\link{wide_to_long}}, \code{\link{long_to_wide}} for
#'   format conversions.
#'
#' @export
prepare_for_tna <- function(data,
                            type = c("sequences", "long", "auto"),
                            state_names = NULL,
                            id_col = "Actor",
                            time_col = "Time",
                            action_col = "Action",
                            validate = TRUE) {
  # Input validation
  stopifnot(is.data.frame(data))
  type <- match.arg(type)

  # Auto-detect format
  if (type == "auto") {
    has_time_cols <- any(grepl("^V[0-9]+$", names(data)))
    has_action_col <- action_col %in% names(data)

    if (has_time_cols && !has_action_col) {
      type <- "sequences"
    } else if (has_action_col && !has_time_cols) {
      type <- "long"
    } else if (has_time_cols && has_action_col) {
      # If both, use long if nrow >> number of unique IDs
      if (id_col %in% names(data)) {
        n_ids <- length(unique(data[[id_col]]))
        if (nrow(data) > n_ids * 2) {
          type <- "long"
        } else {
          type <- "sequences"
        }
      } else {
        type <- "sequences"
      }
    } else {
      stop("Cannot auto-detect data format. Please specify 'type' explicitly.")
    }
  }

  # Convert long to wide if needed
  if (type == "long") {
    if (!action_col %in% names(data)) {
      stop("Action column '", action_col, "' not found in data.")
    }
    if (!id_col %in% names(data)) {
      stop("ID column '", id_col, "' not found in data.")
    }

    data <- long_to_wide(
      data,
      id_col = id_col,
      time_col = time_col,
      action_col = action_col
    )
  }

  # Identify sequence columns (V1, V2, ...)
  seq_cols <- grep("^V[0-9]+$", names(data), value = TRUE)
  if (length(seq_cols) == 0) {
    stop("No sequence columns (V1, V2, ...) found in data.")
  }

  # Extract just the sequence columns
  result <- data[, seq_cols, drop = FALSE]

  # Convert factors to characters
  for (col in names(result)) {
    if (is.factor(result[[col]])) {
      result[[col]] <- as.character(result[[col]])
    }
  }

  # Validate state names if requested
  if (validate && !is.null(state_names)) {
    all_states <- unique(unlist(result, use.names = FALSE))
    all_states <- all_states[!is.na(all_states)]
    unknown_states <- setdiff(all_states, state_names)
    if (length(unknown_states) > 0) {
      warning("Unknown states found: ", paste(unknown_states, collapse = ", "))
    }
  }

  # Ensure proper column ordering
  col_nums <- as.integer(gsub("V", "", seq_cols))
  result <- result[, order(col_nums), drop = FALSE]

  # Convert to data frame
  result <- as.data.frame(result)
  rownames(result) <- NULL

  return(result)
}


#' Convert Action Column to One-Hot Encoding
#'
#' @description
#' Convert a categorical Action column to one-hot (binary indicator) columns.
#'
#' @param data Data frame containing an action column.
#' @param action_col Character. Name of the action column. Default: "Action".
#' @param states Character vector or NULL. States to include as columns.
#'   If NULL, uses all unique values. Default: NULL.
#' @param drop_action Logical. Remove the original action column. Default: TRUE.
#' @param sort_states Logical. Sort state columns alphabetically. Default: FALSE.
#' @param prefix Character. Prefix for state column names. Default: "".
#'
#' @return Data frame with one-hot encoded columns (0/1 integers).
#'
#' @examples
#' long_data <- data.frame(
#'   Actor = rep(1:3, each = 4),
#'   Time = rep(1:4, 3),
#'   Action = sample(c("A", "B", "C"), 12, replace = TRUE)
#' )
#' onehot_data <- action_to_onehot(long_data)
#' head(onehot_data)
#'
#' @export
action_to_onehot <- function(data, action_col = "Action", states = NULL,
                             drop_action = TRUE, sort_states = FALSE,
                             prefix = "") {
  stopifnot(is.data.frame(data))
  stopifnot(action_col %in% names(data))

  if (is.null(states)) {
    states <- unique(data[[action_col]])
    states <- states[!is.na(states)]
  }

  if (sort_states) {
    states <- sort(states)
  }

  for (state in states) {
    col_name <- paste0(prefix, state)
    data[[col_name]] <- as.integer(data[[action_col]] == state)
  }

  if (drop_action) {
    data[[action_col]] <- NULL
  }

  return(data)
}


#' Import One-Hot Encoded Data into Sequence Format
#'
#' @description
#' Converts binary indicator (one-hot) data into the wide sequence format
#' expected by \code{\link{build_network}} and \code{tna::tna()}. Each binary
#' column represents a state; rows where the value is 1 are marked with the
#' column name. Supports optional windowed aggregation.
#'
#' Simultaneous active states are preserved using the same window-span
#' representation as \code{tna::import_onehot()}: each input row/window is
#' expanded to one sequence slot per code and transition counting occurs between
#' windows, not between simultaneous states inside the same row.
#'
#' @param data Data frame with binary (0/1) indicator columns.
#' @param cols Character vector. Names of the one-hot columns to use.
#' @param actor Character or NULL. Name of the actor/ID column. If NULL,
#'   all rows are treated as a single sequence. Default: NULL.
#' @param session Character or NULL. Name of the session column for
#'   sub-grouping within actors. Default: NULL.
#' @param interval Integer or NULL. Number of rows per time point in the
#'   output. If NULL, all rows become a single time point group.
#'   Default: NULL.
#' @param window_size Integer. Number of consecutive rows to aggregate
#'   into each window. Default: 1 (no windowing).
#' @param window_type Character. \code{"non-overlapping"} (fixed, separate
#'   windows) or \code{"overlapping"} (rolling, step = 1).
#'   Default: \code{"non-overlapping"}.
#' @param aggregate Logical. If TRUE, aggregate within each window by
#'   taking the first non-NA indicator per column. Default: FALSE.
#'
#' @return A data frame in wide format with columns named
#'   \code{W{window}_T{time}} where each cell contains a state name or NA.
#'   Attributes \code{windowed}, \code{window_size}, \code{window_span}
#'   are set on the result.
#'
#' @examples
#' # Simple binary data
#' df <- data.frame(
#'   A = c(1, 0, 1, 0, 1),
#'   B = c(0, 1, 0, 1, 0),
#'   C = c(0, 0, 0, 0, 0)
#' )
#' seq_data <- prepare_onehot(df, cols = c("A", "B", "C"))
#'
#' # With actor grouping
#' df$actor <- c(1, 1, 1, 2, 2)
#' seq_data <- prepare_onehot(df, cols = c("A", "B", "C"), actor = "actor")
#'
#' # With windowing
#' seq_data <- prepare_onehot(df, cols = c("A", "B", "C"),
#'                           window_size = 2, window_type = "non-overlapping")
#'
#' @seealso \code{\link{action_to_onehot}} for the reverse conversion.
#'
#' @export
prepare_onehot <- function(data,
                           cols,
                           actor = NULL,
                           session = NULL,
                           interval = NULL,
                           window_size = 3L,
                           window_type = c("non-overlapping", "overlapping"),
                           aggregate = FALSE) {
  stopifnot(is.data.frame(data))
  stopifnot(is.character(cols), length(cols) >= 1)
  stopifnot(all(cols %in% names(data)))
  window_size <- as.integer(window_size)
  stopifnot(window_size >= 1L)
  window_type <- match.arg(window_type)

  interval <- if (is.null(interval)) ncol(data) * nrow(data) else interval
  tna_window_type <- switch(window_type,
                            "non-overlapping" = "tumbling",
                            "overlapping" = "sliding")

  work <- data
  code_values <- lapply(cols, function(col) {
    ifelse(work[[col]] == 1, col, NA_character_)
  })
  names(code_values) <- cols
  work[cols] <- code_values

  actor_col <- if (is.null(actor)) ".actor" else actor
  session_col <- if (is.null(session)) ".session" else session
  if (!is.null(actor)) stopifnot(actor %in% names(data))
  if (!is.null(session)) stopifnot(session %in% names(data))
  if (is.null(actor)) work[[actor_col]] <- 1L
  if (is.null(session)) work[[session_col]] <- 1L

  group_key <- interaction(work[, c(actor_col, session_col), drop = FALSE],
                           drop = TRUE)
  group_key <- as.character(group_key)
  groups <- lapply(unique(group_key), function(key) work[group_key == key, , drop = FALSE])

  lag_chr <- function(x, n) {
    if (n == 0L) return(x)
    if (n >= length(x)) return(rep(NA_character_, length(x)))
    c(rep(NA_character_, n), head(x, -n))
  }

  first_non_missing <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) NA_character_ else x[1L]
  }

  work <- do.call(rbind, lapply(groups, function(g) {
    rownames(g) <- NULL
    if (tna_window_type == "sliding") {
      # seq_len(window_size - 1L) is cleanly empty when window_size == 1L;
      # seq(1, window_size - 1) on R's "decreasing implied" rule would be
      # c(1, 0), which iterates twice with an off-by-one shift.
      for (w in seq_len(window_size - 1L)) {
        g[cols] <- lapply(cols, function(col) {
          x <- g[[col]]
          ifelse(!is.na(x) | !is.na(lag_chr(x, w)), col, NA_character_)
        })
        names(g[cols]) <- cols
      }
      # Drop only the leading row that the lag propagation invalidates.
      # window_size == 1L has no overlap, so no row should be dropped.
      if (window_size > 1L && nrow(g) > 0L) {
        g <- g[-1L, , drop = FALSE]
      }
      g$.window <- floor(seq_len(nrow(g)) - 1L)
    } else {
      g$.window <- floor((seq_len(nrow(g)) - 1L) / window_size)
    }
    g
  }))

  if (aggregate && nrow(work) > 0L) {
    ord <- do.call(order, work[, c(actor_col, session_col, ".window"),
                                drop = FALSE])
    work <- work[ord, , drop = FALSE]
    agg_key <- interaction(work[, c(actor_col, session_col, ".window"),
                                drop = FALSE], drop = TRUE)
    agg_key <- as.character(agg_key)
    work <- do.call(rbind, lapply(unique(agg_key), function(key) {
      g <- work[agg_key == key, , drop = FALSE]
      out <- g[1L, , drop = FALSE]
      out[cols] <- lapply(cols, function(col) first_non_missing(g[[col]]))
      names(out[cols]) <- cols
      out
    }))
    rownames(work) <- NULL
  }

  if (is.null(work) || nrow(work) == 0L) {
    result <- as.data.frame(setNames(replicate(0, character(0), simplify = FALSE),
                                     character(0)))
  } else {
    work$.window_grp <- work$.window %/% interval
    work$.window_idx <- work$.window %% interval

    long_rows <- lapply(seq_len(nrow(work)), function(i) {
      data.frame(
        .actor = work[[actor_col]][i],
        .session = work[[session_col]][i],
        .window_grp = work$.window_grp[i],
        .window_idx = work$.window_idx[i],
        value = unname(unlist(work[i, cols, drop = FALSE], use.names = FALSE)),
        stringsAsFactors = FALSE
      )
    })
    long <- do.call(rbind, long_rows)
    obs_key <- interaction(long[, c(".actor", ".session", ".window_grp",
                                    ".window_idx"), drop = FALSE],
                           drop = TRUE)
    long$.obs <- ave(seq_len(nrow(long)), obs_key, FUN = seq_along)
    id_key <- interaction(long[, c(".actor", ".session", ".window_grp"),
                               drop = FALSE], drop = TRUE)
    id_key <- as.character(id_key)
    wide_list <- lapply(unique(id_key), function(key) {
      g <- long[id_key == key, , drop = FALSE]
      row <- data.frame(.row = 1L)
      for (i in seq_len(nrow(g))) {
        col_name <- sprintf("W%d_T%d", g$.window_idx[i], g$.obs[i])
        row[[col_name]] <- g$value[i]
      }
      row$.row <- NULL
      row
    })

    all_cols <- unique(unlist(lapply(wide_list, names), use.names = FALSE))
    wide_list <- lapply(wide_list, function(row) {
      missing <- setdiff(all_cols, names(row))
      for (m in missing) {
        row[[m]] <- NA_character_
      }
      row[, all_cols, drop = FALSE]
    })
    result <- do.call(rbind, wide_list)
    rownames(result) <- NULL
  }

  if (ncol(result) > 0L) {
    result[] <- lapply(result, as.character)
  }

  attr(result, "windowed") <- TRUE
  attr(result, "window_size") <- window_size^(!aggregate)
  attr(result, "window_span") <- length(cols)
  attr(result, "codes") <- cols

  result
}
