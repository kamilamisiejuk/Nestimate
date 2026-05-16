#' Build a Network
#'
#' @description
#' Universal network estimation function that supports both transition
#' networks (relative, frequency, co-occurrence) and association networks
#' (correlation, partial correlation, graphical lasso). Uses the global
#' estimator registry, so custom estimators can also be used.
#'
#' @param data Data frame (sequences or per-observation frequencies) or a
#'   square symmetric matrix (correlation or covariance).
#' @param method Character. Required. Name of a registered estimator.
#'   Built-in methods: \code{"relative"}, \code{"frequency"},
#'   \code{"co_occurrence"}, \code{"cor"}, \code{"pcor"}, \code{"glasso"}.
#'   Aliases: \code{"tna"} and \code{"transition"} map to \code{"relative"};
#'   \code{"ftna"} and \code{"counts"} map to \code{"frequency"};
#'   \code{"cna"} maps to \code{"co_occurrence"};
#'   \code{"corr"} and \code{"correlation"} map to \code{"cor"};
#'   \code{"partial"} maps to \code{"pcor"};
#'   \code{"ebicglasso"} and \code{"regularized"} map to \code{"glasso"}.
#' @param params Named list. Method-specific parameters passed to the estimator
#'   function (e.g. \code{list(gamma = 0.5)} for glasso, or
#'   \code{list(format = "wide")} for transition methods). This is the key
#'   composability feature: downstream functions like bootstrap or grid search
#'   can store and replay the full params list without knowing method internals.
#'   Transition estimators accept tna-style sequence options such as
#'   \code{weighted}, \code{begin_state}, \code{end_state}, and \code{concat}.
#'   Column-like entries in \code{params} (\code{action}, \code{id},
#'   \code{id_col}, \code{time}, \code{session}, \code{order}, \code{codes},
#'   and \code{group}) are resolved before format detection and must name
#'   existing columns. If the same column role is supplied both directly and
#'   through \code{params}, the names must agree.
#' @param labels Optional name -> label remap applied after construction.
#'   Accepts a 2-column data.frame \code{(name, label)}, a named character
#'   vector \code{c(name = "label")}, or a named list. Rewrites
#'   \code{$nodes$label} and \code{dimnames(weights)}. Unmapped names pass
#'   through unchanged.
#' @param scaling Character vector or NULL. Post-estimation scaling to apply
#'   (in order). Options: \code{"minmax"}, \code{"max"}, \code{"rank"},
#'   \code{"normalize"}. Can combine: \code{c("rank", "minmax")}.
#'   Default: \code{NULL} (no scaling).
#' @param threshold Numeric. Absolute values below this are set to zero in the
#'   result matrix. Default: 0 (no thresholding).
#' @param level Character or NULL. Multilevel decomposition for association
#'   methods. One of \code{NULL}, \code{"between"}, \code{"within"},
#'   \code{"both"}. Requires \code{id_col}. Default: \code{NULL}.
#' @param actor Character. Name of the actor/person ID column for sequence
#'   grouping. Default: \code{NULL}.
#' @param action Character. Name of the action/state column (long format).
#'   Default: \code{NULL}.
#' @param time Character. Name of the time column (long format).
#'   Default: \code{NULL}.
#' @param session Character. Name of the session column. Default: \code{NULL}.
#' @param order Character. Name of the ordering column. Default: \code{NULL}.
#' @param codes Character vector. Column names of one-hot encoded states
#'   (for onehot format). Default: \code{NULL}.
#' @param group Character. Name of a grouping column for per-group networks.
#'   Returns a \code{netobject_group} (named list of netobjects).
#'   Default: \code{NULL}.
#' @param format Character. Input format: \code{"auto"}, \code{"wide"},
#'   \code{"long"}, or \code{"onehot"}. Default: \code{"auto"}.
#' @param window_size Integer. Window size for one-hot windowing.
#'   Default: \code{3L}.
#' @param mode Character. Windowing mode: \code{"non-overlapping"} or
#'   \code{"overlapping"}. Default: \code{"non-overlapping"}.
#' @param time_threshold Numeric. Maximum time gap (seconds) for long format
#'   session splitting. Default: \code{900}.
#' @param predictability Logical. If \code{TRUE} (default), compute and store
#'   node predictability (R-squared) for undirected association methods
#'   (glasso, pcor, cor). Stored in \code{$predictability} and auto-displayed
#'   as donuts by \code{cograph::splot()}.
#' @param state_cols Character vector or \code{NULL}. Explicit names of columns
#'   to classify as state columns in the returned netobject's \code{$data}
#'   slot. When provided, all other columns of the cleaned input go to
#'   \code{$metadata}. Auto-detection (values-in-nodes heuristic) is bypassed.
#'   Use this when a metadata column happens to contain values that overlap
#'   with node names (e.g. condition labels \code{"A","B","C"} and nodes
#'   \code{"A","B","C"}) and auto-detection would misclassify it. Default:
#'   \code{NULL} (auto-detect).
#' @param metadata_cols Character vector or \code{NULL}. Explicit names of
#'   columns to force into the \code{$metadata} slot. The remaining columns
#'   are auto-detected as state via the values-in-nodes rule. Cannot overlap
#'   with \code{state_cols}. Default: \code{NULL}.
#' @param ... Additional arguments passed to the estimator function.
#'
#' @return An object of class \code{c("netobject", "cograph_network")} containing:
#' \describe{
#'   \item{data}{The input data used for estimation, as a data frame.}
#'   \item{weights}{The estimated network weight matrix.}
#'   \item{nodes}{Data frame with columns \code{id}, \code{label}, \code{name},
#'     \code{x}, \code{y}. Node labels are in \code{$nodes$label}.}
#'   \item{edges}{Data frame of non-zero edges with integer \code{from}/\code{to}
#'     (node IDs) and numeric \code{weight}.}
#'   \item{directed}{Logical. Whether the network is directed.}
#'   \item{method}{The resolved method name.}
#'   \item{params}{The params list used (for reproducibility).}
#'   \item{scaling}{The scaling applied (or NULL).}
#'   \item{threshold}{The threshold applied.}
#'   \item{n_nodes}{Number of nodes.}
#'   \item{n_edges}{Number of non-zero edges.}
#'   \item{level}{Decomposition level used (or NULL).}
#'   \item{meta}{List with \code{source}, \code{layout}, and \code{tna} metadata
#'     (cograph-compatible).}
#'   \item{node_groups}{Node groupings data frame, or NULL.}
#'   \item{predictability}{Named numeric vector of R-squared predictability
#'     values per node (for undirected association methods when
#'     \code{predictability = TRUE}). NULL for directed methods.}
#' }
#' Method-specific extras (e.g. \code{precision_matrix}, \code{cor_matrix},
#' \code{frequency_matrix}, \code{lambda_selected}, etc.) are preserved
#' from the estimator output.
#'
#' When \code{level = "both"}, returns an object of class
#' \code{"netobject_ml"} with \code{$between} and \code{$within}
#' sub-networks and a \code{$method} field.
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Resolves method aliases to canonical names.
#'   \item Validates explicit column arguments before any format guessing.
#'   \item Retrieves the estimator function from the global registry.
#'   \item For association methods with \code{level} specified, decomposes
#'     the data (between-person means or within-person centering).
#'   \item Calls the estimator: \code{do.call(fn, c(list(data = data), params))}.
#'   \item Applies scaling and thresholding to the result matrix.
#'   \item Extracts edges and constructs the \code{netobject}.
#' }
#'
#' For long-format transition data, supplying \code{action} without
#' \code{actor} is allowed and treats all rows as one sequence in row/time
#' order. The function warns because a one-sequence transition network is not
#' recommended and cannot be validated by bootstrap or other confirmatory
#' tests.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
#' net <- build_network(seqs, method = "relative")
#' net
#' \donttest{
#' # Transition network (relative probabilities)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' print(net)
#'
#' # Association network (glasso)
#' freq_data <- convert_sequence_format(seqs, format = "frequency")
#' net_glasso <- build_network(freq_data, method = "glasso",
#'                              params = list(gamma = 0.5, nlambda = 50))
#'
#' # With scaling
#' net_scaled <- build_network(seqs, method = "relative",
#'                              scaling = c("rank", "minmax"))
#' }
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{bootstrap_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @importFrom utils capture.output
#' @export
build_network <- function(data,
                          method,
                          actor = NULL,
                          action = NULL,
                          time = NULL,
                          session = NULL,
                          order = NULL,
                          codes = NULL,
                          group = NULL,
                          format = "auto",
                          window_size = 3L,
                          mode = c("non-overlapping", "overlapping"),
                          scaling = NULL,
                          threshold = 0,
                          level = NULL,
                          time_threshold = 900,
                          predictability = TRUE,
                          state_cols = NULL,
                          metadata_cols = NULL,
                          params = list(),
                          labels = NULL,
                          ...) {
  # --- Coerce sequence matrices (character / logical / integer) to data.frame ---
  # Numeric square matrices are left alone — they may already be transition
  # matrices for downstream estimators that accept them.
  if (is.matrix(data) && !is.numeric(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }

  # --- Early dispatch for net_clustering objects ---
  if (inherits(data, "net_clustering")) {
    if (missing(method)) method <- data$network_method %||% "relative"
    return(.build_network_clustering(data, method = method, ...))
  }

  # --- Early dispatch for net_mmm objects ---
  if (inherits(data, "net_mmm")) {
    if (missing(method)) method <- data$network_method %||% "relative"
    resolved <- .resolve_method_alias(method)
    dots <- list(...)
    if (resolved != "relative") {
      # Re-build per-component networks from hard assignments using requested method
      raw_data    <- data$data %||% data$models[[1L]]$data
      assignments <- data$assignments
      k_comp      <- data$k
      # Merge stored build_args with caller's ...; caller takes precedence
      call_args <- if (!is.null(data$build_args)) modifyList(data$build_args, dots) else dots
      nets <- lapply(seq_len(k_comp), function(m) {
        sub <- raw_data[assignments == m, , drop = FALSE]
        net <- do.call(build_network, c(list(data = sub, method = method), call_args))
        # Inject EM-fitted initials only for directed sequence methods
        if (resolved %in% c("relative", "frequency", "attention")) {
          net$initial <- data$models[[m]]$initial
        }
        net
      })
      names(nets) <- paste0("Cluster ", seq_len(k_comp))
      class(nets) <- "netobject_group"
      return(.attach_mmm_clustering(nets, data, full_data = raw_data))
    }
    if (length(dots)) {
      stop(
        "Unsupported argument(s) for build_network(net_mmm) with the ",
        "pre-built relative MMM networks: ",
        paste(names(dots), collapse = ", "), ". ",
        "Use a non-relative `method` to rebuild per-cluster networks, or ",
        "set these options before fitting the MMM.",
        call. = FALSE
      )
    }
    # Default: wrap pre-built "relative" models (retain $initial from EM)
    nets <- data$models
    if (is.null(names(nets))) names(nets) <- paste0("Cluster ", seq_along(nets))
    class(nets) <- "netobject_group"
    return(.attach_mmm_clustering(nets, data,
                                  full_data = data$data %||% data$models[[1L]]$data))
  }

  stopifnot(is.character(method), length(method) == 1)
  stopifnot(is.list(params))
  stopifnot(is.numeric(threshold), length(threshold) == 1, threshold >= 0)
  mode <- match.arg(mode)

  if (!is.null(state_cols)) {
    stopifnot(is.character(state_cols))
  }
  if (!is.null(metadata_cols)) {
    stopifnot(is.character(metadata_cols))
  }
  if (!is.null(state_cols) && !is.null(metadata_cols)) {
    overlap <- intersect(state_cols, metadata_cols)
    if (length(overlap)) {
      stop(
        "state_cols and metadata_cols overlap: ",
        paste(overlap, collapse = ", "), ". A column must be either state ",
        "or metadata, not both.",
        call. = FALSE
      )
    }
  }
  if (is.data.frame(data)) {
    missing_state <- setdiff(state_cols,    names(data))
    missing_meta  <- setdiff(metadata_cols, names(data))
    if (length(missing_state)) {
      stop("state_cols not found in data: ",
           paste(missing_state, collapse = ", "), call. = FALSE)
    }
    if (length(missing_meta)) {
      stop("metadata_cols not found in data: ",
           paste(missing_meta, collapse = ", "), call. = FALSE)
    }
  }

  # Merge ... into params (... takes precedence)
  dots <- list(...)
  if (length(dots) > 0) {
    params <- modifyList(params, dots)
  }

  # Resolve method aliases early (needed for format detection)
  method <- .resolve_method_alias(method)

  if (is.data.frame(data)) {
    canonical <- .canonicalize_build_network_params(
      params,
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, group = group, format = format
    )
    actor <- canonical$actor
    action <- canonical$action
    time <- canonical$time
    session <- canonical$session
    order <- canonical$order
    codes <- canonical$codes
    group <- canonical$group
    format <- canonical$format

    .validate_build_network_columns(
      data,
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, group = group,
      state_cols = state_cols, metadata_cols = metadata_cols
    )
    .validate_build_network_params(data, params, actor = actor,
                                   action = action, time = time,
                                   codes = codes)
  }

  # ---- Auto-match standard column names (case-insensitive) ----
  if (is.data.frame(data)) {
    col_lower <- tolower(names(data))
    .match1 <- function(name) {
      hit <- which(col_lower == name)
      if (length(hit) == 1L) names(data)[hit] else NULL
    }
    if (is.null(action))  action  <- .match1("action")
    if (is.null(time))    time    <- .match1("time")
    if (is.null(session)) session <- .match1("session") %||% .match1("session_id")
  }

  # ---- Group dispatch: per-group networks ----
  if (!is.null(group)) {
    stopifnot(is.character(group))
    transition_methods <- c("relative", "frequency", "co_occurrence", "attention")
    if (method %in% transition_methods && is.data.frame(data) &&
        is.null(.param_get(params, "alphabet"))) {
      if (!is.null(action) && action %in% names(data)) {
        vals <- .clean_states(as.character(data[[action]]))
      } else {
        exclude <- c(group, actor, session)
        state_names <- setdiff(names(data), exclude)
        vals <- .clean_states(as.character(unlist(data[, state_names, drop = FALSE])))
      }
      params$alphabet <- sort(unique(c(
        vals[!is.na(vals)],
        .param_get(params, "begin_state"),
        .param_get(params, "end_state")
      )))
    }
    if (length(group) == 1L) {
      grp_key <- data[[group]]
    } else {
      grp_key <- interaction(data[, group, drop = FALSE], sep = "-",
                             drop = TRUE)
    }
    grp_levels <- unique(grp_key)
    # Drop group column(s) from sub-data so they don't become state cols
    drop_cols <- intersect(group, names(data))
    # Strip group from params so the recursive per-group call does not
    # try to re-canonicalize a `group` column that has been dropped from
    # `sub`. The parent call already consumed it.
    child_params <- if (is.list(params)) params[setdiff(names(params), "group")] else params
    nets <- lapply(grp_levels, function(g) {
      sub <- data[grp_key == g, , drop = FALSE]
      if (length(drop_cols) > 0L) {
        sub <- sub[, setdiff(names(sub), drop_cols), drop = FALSE]
      }
      build_network(
        sub, method = method, actor = actor, action = action,
        time = time, session = session, order = order, codes = codes,
        group = NULL, format = format, window_size = window_size,
        mode = mode, scaling = scaling, threshold = threshold,
        level = level, time_threshold = time_threshold,
        predictability = predictability,
        state_cols = state_cols, metadata_cols = metadata_cols,
        params = child_params, ...
      )
    })
    names(nets) <- as.character(grp_levels)
    attr(nets, "group_col") <- group
    class(nets) <- "netobject_group"
    return(nets)
  }

  # ---- Auto-detect input format ----
  is_onehot <- FALSE
  onehot_codes <- NULL
  if (format == "auto" && is.data.frame(data)) {
    if (!is.null(codes)) {
      # Explicit codes = one-hot
      is_onehot <- TRUE
      format <- "onehot"
    } else if (!is.null(action) && action %in% names(data)) {
      # Has action column = long event log
      format <- "long"
    } else {
      # Check if all columns (excluding actor/session) are binary
      exclude <- c(actor, session)
      check_cols <- setdiff(names(data), exclude)
      if (length(check_cols) >= 2L && all(vapply(
        data[, check_cols, drop = FALSE],
        function(x) is.numeric(x) && all(x[!is.na(x)] %in% c(0, 1)),
        logical(1)
      ))) {
        is_onehot <- TRUE
        format <- "onehot"
      } else {
        format <- "wide"
      }
    }
  }

  # ---- Long format: prepare event log data ----
  if (format == "long" && !is.null(action) && is.data.frame(data) &&
      action %in% names(data)) {
    if (is.null(actor)) {
      warning(
        "A network with one long sequence is not recommended and can't be ",
        "validated using bootstrap and other confirmatory testings.",
        call. = FALSE
      )
    }
    prep_args <- list(data = data, action = action)
    if (!is.null(actor)) prep_args$actor <- actor
    if (!is.null(time)) prep_args$time <- time
    if (!is.null(session)) prep_args$session <- session # nocov start
    if (!is.null(order)) prep_args$order <- order # nocov end
    prep_args$time_threshold <- time_threshold

    prepared <- do.call(prepare, prep_args)
    data <- prepared$sequence_data
    format <- "wide"
    params$format <- "wide"
  } else {
    prepared <- NULL
  }

  # ---- One-hot format: dispatch to wtna/cna path ----
  if (format == "onehot" || is_onehot) {
    resolved_codes <- .resolve_codes(
      as.data.frame(data), codes,
      exclude = c(actor, session)
    )

    # Warning for one-hot without windowing or session
    if (window_size <= 1L && is.null(session)) {
      warning(
        "One-hot data without windowing or sessions produces sparse networks. ",
        "Consider setting window_size > 1 or providing session.",
        call. = FALSE
      )
    }

    # Build grouping from actor + session
    grp_col <- NULL
    if (!is.null(actor) || !is.null(session)) {
      grp_parts <- c(actor, session)
      grp_col <- grp_parts
    }

    params$codes <- resolved_codes
    onehot_codes <- resolved_codes
    params$window_size <- window_size
    params$mode <- mode
    if (!is.null(grp_col)) params$actor <- grp_col
    # Let the estimator handle it
  } else {
    # Wide sequence data: pass format through
    if (!"format" %in% names(params)) params$format <- format
  }

  # ---- Auto-convert sequences to frequencies for association methods ----
  # Association methods (cor, pcor, glasso, ising) require numeric data.
  # If the data is character/factor sequence data, convert to per-row state
  # frequency counts automatically so users can pass sequences directly.
  assoc_methods <- c("cor", "pcor", "glasso")
  if (method %in% assoc_methods && is.data.frame(data)) {
    exclude <- intersect(c(actor, session), names(data))
    check_cols <- setdiff(names(data), exclude)
    has_char <- length(check_cols) >= 2L && any(vapply(
      data[, check_cols, drop = FALSE],
      function(x) is.character(x) || is.factor(x), logical(1)
    ))
    if (has_char) {
      freq_data <- convert_sequence_format(
        data, id_col = if (length(exclude) == 0L) character(0) else exclude,
        format = "frequency"
      )
      # Drop rid and ID columns, keep only numeric frequency columns
      drop <- c("rid", exclude)
      data <- freq_data[, setdiff(names(freq_data), drop), drop = FALSE]
      params$format <- "wide"
    }
  }

  # ---- Multilevel decomposition ----
  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col") %||%
    actor

  # Validate level parameter
  if (!is.null(level)) {
    level <- match.arg(level, c("between", "within", "both"))
    if (is.null(id_col)) {
      stop("'id' is required when 'level' is specified.", call. = FALSE)
    }
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame when 'level' is specified.",
           call. = FALSE)
    }
  }

  # Validate scaling
  if (!is.null(scaling)) {
    valid_scaling <- c("minmax", "max", "rank", "normalize")
    bad <- setdiff(scaling, valid_scaling)
    if (length(bad) > 0) {
      stop("Unknown scaling method(s): ", paste(bad, collapse = ", "),
           ". Options: ", paste(valid_scaling, collapse = ", "),
           call. = FALSE)
    }
  }

  # Get estimator from registry
  estimator <- get_estimator(method)

  # level = "both": recursive dispatch
  if (identical(level, "both")) {
    between <- build_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "between",
      predictability = predictability
    )
    within_net <- build_network(
      data, method = method, params = params, scaling = scaling,
      threshold = threshold, level = "within",
      predictability = predictability
    )
    result <- list(between = between, within = within_net, method = method)
    class(result) <- "netobject_ml"
    return(result)
  }

  # Multilevel decomposition for association methods
  if (!is.null(level) && !estimator$directed) {
    data <- .decompose_multilevel(data, id_col = id_col, level = level)
  }

  # Call estimator
  est_result <- do.call(estimator$fn, c(list(data = data), params))

  # Validate estimator output
  if (!is.list(est_result) ||
      is.null(est_result$matrix) ||
      is.null(est_result$nodes) ||
      is.null(est_result$directed)) {
    stop("Estimator '", method,
         "' must return a list with 'matrix', 'nodes', and 'directed'.",
         call. = FALSE)
  }

  net_matrix <- est_result$matrix
  nodes <- est_result$nodes
  directed <- est_result$directed

  # Apply scaling
  if (!is.null(scaling)) {
    net_matrix <- .apply_scaling(
      net_matrix, scaling,
      include_zeros = method %in% c("relative", "frequency", "co_occurrence",
                                    "attention")
    )
  }

  # Apply threshold
  if (threshold > 0) {
    net_matrix[abs(net_matrix) < threshold] <- 0
  }

  # Extract edges
  edges <- .extract_edges_from_matrix(net_matrix, directed = directed)

  # Split cleaned_data into state-only $data and $metadata
  # A column is a state column if all its non-void/non-NA values are in nodes
  raw_data <- est_result$cleaned_data
  metadata <- NULL
  # Use prepared$meta_data when available (long format path)
  if (!is.null(prepared) && !is.null(prepared$meta_data)) {
    md <- prepared$meta_data
    # Drop internal .session_id column
    md <- md[, setdiff(names(md), ".session_id"), drop = FALSE]
    if (ncol(md) > 0L) metadata <- md
  }
  if (is.data.frame(raw_data)) {
    # If the user passed explicit overrides that name columns present in
    # raw_data, those take priority over the values-in-nodes heuristic.
    auto_onehot_state <- if (is.null(state_cols)) {
      intersect(onehot_codes, names(raw_data))
    } else {
      character(0)
    }
    user_state <- union(intersect(state_cols, names(raw_data)),
                        auto_onehot_state)
    user_meta  <- intersect(metadata_cols, names(raw_data))

    if (length(user_state) > 0L) {
      # Explicit state_cols: those ARE state; rest of raw_data goes to metadata.
      resolved_state <- user_state
      resolved_extra <- setdiff(names(raw_data), user_state)
    } else {
      is_state_col <- vapply(raw_data, function(col) {
        vals <- .clean_states(as.character(col))
        vals <- vals[!is.na(vals)]
        length(vals) > 0L && all(vals %in% nodes)
      }, logical(1))
      # Honor metadata_cols by forcing those FALSE in the auto-detection.
      if (length(user_meta) > 0L) is_state_col[user_meta] <- FALSE
      resolved_state <- names(raw_data)[is_state_col]
      resolved_extra <- names(raw_data)[!is_state_col]
    }

    if (length(resolved_extra) > 0L) {
      extra_df <- raw_data[, resolved_extra, drop = FALSE]
      metadata <- if (is.null(metadata)) extra_df else cbind(metadata, extra_df)
      raw_data <- raw_data[, resolved_state, drop = FALSE]
    }
    # Clean void/missing markers in character/factor state columns
    if (is.data.frame(raw_data)) {
      raw_data[] <- lapply(raw_data, function(col) {
        if (is.character(col) || is.factor(col)) {
          .clean_states(as.character(col))
        } else {
          col
        }
      })
    }
  }

  # Build unified netobject / cograph_network
  nodes_df <- data.frame(
    id = seq_along(nodes),
    label = nodes,
    name = nodes,
    stringsAsFactors = FALSE
  )

  result <- list(
    data = raw_data,
    metadata = metadata,
    weights = net_matrix,
    nodes = nodes_df,
    edges = edges,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = length(nodes),
    n_edges = nrow(edges),
    level = level,
    build_args = list(
      actor = actor, action = action, time = time, session = session,
      order = order, codes = codes, format = format,
      window_size = window_size, mode = mode
    ),
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna = list(method = method)
    ),
    node_groups = NULL
  )

  # Carry over method-specific extras (exclude cleaned_data — identical to $data)
  extras <- setdiff(names(est_result),
                    c("matrix", "nodes", "directed", "cleaned_data"))
  for (key in extras) {
    result[[key]] <- est_result[[key]]
  }

  # Auto-compute predictability (R²) for undirected association methods.
  # Stored as a named numeric vector for backward compatibility with downstream
  # consumers (e.g. cograph::splot pie ring). Users calling predictability()
  # directly get the full tidy data.frame with R2 + RMSE.
  if (isTRUE(predictability) && !directed) {
    class(result) <- c("netobject", "cograph_network")
    pred_df <- tryCatch(predictability(result), error = function(e) NULL)
    if (!is.null(pred_df)) {
      vec <- pred_df$R2
      names(vec) <- pred_df$node
      result$predictability <- vec
    }
  }

  result <- structure(result, class = c("netobject", "cograph_network"))
  if (!is.null(labels)) result <- .apply_node_labels(result, labels)
  result
}


#' @noRd
.canonicalize_build_network_params <- function(params,
                                               actor = NULL,
                                               action = NULL,
                                               time = NULL,
                                               session = NULL,
                                               order = NULL,
                                               codes = NULL,
                                               group = NULL,
                                               format = "auto") {
  .take_param <- function(current, param_name, formal_name = param_name) {
    value <- .param_get(params, param_name)
    if (is.null(value)) {
      return(current)
    }
    if (is.null(current)) {
      return(value)
    }
    if (!identical(as.character(current), as.character(value))) {
      stop(
        "'", formal_name, "' and 'params$", param_name,
        "' specify different columns: ",
        paste(as.character(current), collapse = ", "), " vs ",
        paste(as.character(value), collapse = ", "),
        call. = FALSE
      )
    }
    current
  }

  params_format <- .param_get(params, "format")
  if (!is.null(params_format)) {
    if (!is.character(params_format) || length(params_format) != 1L) {
      stop("'params$format' must be a single character value.", call. = FALSE)
    }
    params_format <- match.arg(params_format, c("auto", "wide", "long", "onehot"))
    if (!identical(format, "auto") && !identical(format, params_format)) {
      stop(
        "'format' and 'params$format' specify different values: ",
        format, " vs ", params_format,
        call. = FALSE
      )
    }
    if (identical(format, "auto")) {
      format <- params_format
    }
  }

  list(
    actor = .take_param(.take_param(actor, "id", formal_name = "actor"),
                        "actor"),
    action = .take_param(action, "action"),
    time = .take_param(time, "time"),
    session = .take_param(session, "session"),
    order = .take_param(order, "order"),
    codes = .take_param(codes, "codes"),
    group = .take_param(group, "group"),
    format = format
  )
}


#' @noRd
.param_get <- function(params, name, default = NULL) {
  if (!is.list(params) || !name %in% names(params)) {
    return(default)
  }
  params[[name]]
}


#' @noRd
.validate_build_network_columns <- function(data, ...) {
  args <- list(...)
  cols <- names(data)

  missing_by_arg <- lapply(names(args), function(arg) {
    value <- args[[arg]]
    if (is.null(value)) {
      return(character(0))
    }
    if (!is.character(value)) {
      stop("'", arg, "' must be a character vector of column name(s).",
           call. = FALSE)
    }
    setdiff(value, cols)
  })
  names(missing_by_arg) <- names(args)

  missing_by_arg <- missing_by_arg[lengths(missing_by_arg) > 0L]
  if (length(missing_by_arg) == 0L) {
    return(invisible(TRUE))
  }

  first_arg <- names(missing_by_arg)[1L]
  missing_cols <- missing_by_arg[[1L]]
  available <- paste(cols, collapse = ", ")
  stop(
    "'", first_arg, "' column",
    if (length(missing_cols) == 1L) " " else "s ",
    "not found in data: ", paste(missing_cols, collapse = ", "),
    ". Available columns: ", available,
    call. = FALSE
  )
}


#' @noRd
.validate_build_network_params <- function(data, params,
                                           actor = NULL,
                                           action = NULL,
                                           time = NULL,
                                           codes = NULL) {
  column_params <- intersect(
    names(params),
    c("action", "id", "time", "cols", "codes", "actor", "id_col",
      "session", "order", "group")
  )
  missing_by_arg <- lapply(column_params, function(arg) {
    value <- params[[arg]]
    if (is.null(value)) {
      return(character(0))
    }
    if (!is.character(value)) {
      stop("'params$", arg, "' must be a character vector of column name(s).",
           call. = FALSE)
    }
    setdiff(value, names(data))
  })
  names(missing_by_arg) <- column_params
  missing_by_arg <- missing_by_arg[lengths(missing_by_arg) > 0L]
  if (length(missing_by_arg) > 0L) {
    first_arg <- names(missing_by_arg)[1L]
    missing_cols <- missing_by_arg[[1L]]
    stop(
      "'params$", first_arg, "' column",
      if (length(missing_cols) == 1L) " " else "s ",
      "not found in data: ", paste(missing_cols, collapse = ", "),
      ". Available columns: ", paste(names(data), collapse = ", "),
      call. = FALSE
    )
  }

  .check_param_conflict <- function(formal_value, param_name,
                                    formal_name = param_name) {
    if (is.null(formal_value) || is.null(params[[param_name]])) {
      return(invisible(TRUE))
    }
    if (!identical(as.character(formal_value), as.character(params[[param_name]]))) {
      stop(
        "'", formal_name, "' and 'params$", param_name,
        "' specify different columns: ",
        paste(as.character(formal_value), collapse = ", "), " vs ",
        paste(as.character(params[[param_name]]), collapse = ", "),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  .check_param_conflict(action, "action")
  .check_param_conflict(time, "time")
  .check_param_conflict(codes, "codes")
  .check_param_conflict(actor, "actor")
  .check_param_conflict(actor, "id", formal_name = "actor")
  .check_param_conflict(actor, "id_col", formal_name = "actor")
}


# ---- S3 methods ----

#' Print Method for Network Object
#'
#' @param x A \code{netobject}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
#' net <- build_network(seqs, method = "relative")
#' print(net)
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B")
#' )
#' net <- build_network(seqs, method = "relative")
#' print(net)
#' }
#'
#' @export
print.netobject <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    ising         = "Ising Model Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)"
  )

  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label   <- if (x$directed) " [directed]" else " [undirected]"
  level_label <- if (!is.null(x$level)) sprintf(" [%s-person]", x$level) else ""
  cat(label, dir_label, level_label, "\n", sep = "")

  if (!is.null(x$n)) cat(sprintf("  Sample size: %d\n", x$n))

  # ---- Weight summary (one line) ----
  mat <- x$weights
  if (x$directed) {
    nz <- mat[mat != 0]
  } else {
    nz <- mat[mat != 0 & row(mat) <= col(mat)]
  }
  if (length(nz) > 0) {
    is_assoc <- x$method %in% c("cor", "pcor", "glasso", "ising")
    if (is_assoc) {
      cat(sprintf("  Weights: [%.3f, %.3f]  |  +%d / -%d edges\n",
                  min(nz), max(nz), sum(nz > 0), sum(nz < 0)))
    } else {
      cat(sprintf("  Weights: [%.3f, %.3f]  |  mean: %.3f\n",
                  min(nz), max(nz), mean(nz)))
    }
  }

  # ---- Weight matrix ----
  cat("\n  Weight matrix:\n")
  digits <- if (all(nz == floor(nz))) 0L else 3L
  mat_r  <- round(mat, digits)
  labels <- x$nodes$label
  dimnames(mat_r) <- list(labels, labels)
  formatted <- capture.output(print(mat_r))
  cat(paste0("  ", formatted, collapse = "\n"), "\n")

  # ---- Initial probabilities ----
  if (!is.null(x$initial) && length(x$initial) > 0) {
    cat("\n  Initial probabilities:\n")
    init  <- x$initial
    ord   <- order(init, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(init)
    vapply(ord, function(i) {
      bars <- if (max_v > 0) strrep("\u2588", round(init[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-12s  %.3f  %s\n", names(init)[i], init[i], bars))
      invisible("")
    }, character(1L))
  }

  # ---- Predictability (R\u00b2) ----
  if (!is.null(x$predictability) && length(x$predictability) > 0) {
    cat("\n  Predictability (R\u00b2):\n")
    pred <- x$predictability
    ord  <- order(pred, decreasing = TRUE)
    bar_w <- 40L
    max_v <- max(pred)
    vapply(ord, function(i) {
      bars <- if (max_v > 0) strrep("\u2588", round(pred[i] / max_v * bar_w)) else ""
      cat(sprintf("  %-12s  %.3f  %s\n", names(pred)[i], pred[i], bars))
      invisible("")
    }, character(1L))
  }

  # ---- Method-specific params ----
  if (x$method == "glasso" && !is.null(x$gamma)) {
    cat(sprintf("\n  Gamma: %.2f  |  Lambda: %.4f\n", x$gamma, x$lambda_selected))
  }
  if (x$method == "ising") {
    cat(sprintf("\n  Gamma: %.2f  |  Rule: %s\n", x$gamma, x$rule))
    if (!is.null(x$thresholds)) {
      thr_rng <- range(x$thresholds)
      cat(sprintf("  Thresholds: [%.3f, %.3f]\n", thr_rng[1], thr_rng[2]))
    }
  }
  if (!is.null(x$scaling))  cat(sprintf("\n  Scaling: %s\n",   paste(x$scaling, collapse = " -> ")))
  if (x$threshold > 0)      cat(sprintf("  Threshold: %g\n",  x$threshold))

  invisible(x)
}


#' Print Method for Group Network Object
#'
#' Compact summary of a \code{netobject_group}. Header surfaces the source
#' (a clustering attached by \code{\link{cluster_network}} or
#' \code{\link{cluster_mmm}}, or a plain split by \code{group_col}). The
#' per-group table carries node and edge counts, weight range, and -- when
#' a clustering attribute is present -- N and percentage of sequences per
#' cluster (matching the layout used by \code{\link{print.net_clustering}}
#' and \code{\link{print.net_mmm}}).
#'
#' @param x A \code{netobject_group}.
#' @param digits Integer. Decimal places for the weight summary. Default
#'   \code{3}. Non-breaking: \code{print(x)} keeps the same shape as
#'   before, with the addition of a weight-range column.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A","B"), V2 = c("B","A","B","A"),
#'                    grp = c("X","X","Y","Y"))
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' print(nets)
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B","A"),
#'   V2 = c("B","C","B","A","C","B"),
#'   V3 = c("C","A","C","B","A","C"),
#'   grp = c("X","X","X","Y","Y","Y")
#' )
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' print(nets)
#' }
#'
#' @export
print.netobject_group <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  grps   <- names(x)
  if (is.null(grps) || any(grps == "")) {
    grps <- paste0("Group ", seq_along(x))
  }
  k <- length(grps)
  cl <- attr(x, "clustering")

  # ---- Header --------------------------------------------------------
  if (inherits(cl, "net_mmm_clustering")) {
    cat(sprintf("Group Networks (%d clusters from MMM)\n", k))
  } else if (inherits(cl, "net_clustering")) {
    cat(sprintf("Group Networks (%d clusters via %s / %s)\n",
                k, cl$method, cl$dissimilarity))
  } else {
    group_col <- attr(x, "group_col")
    if (!is.null(group_col)) {
      cat(sprintf("Group Networks (%d groups, group_col: %s)\n",
                  k, group_col))
    } else {
      cat(sprintf("Group Networks (%d groups)\n", k))
    }
  }

  # ---- Per-group table ----------------------------------------------
  nodes_v <- integer(k)
  edges_v <- integer(k)
  wmin_v  <- numeric(k)
  wmax_v  <- numeric(k)
  for (i in seq_len(k)) {
    net <- x[[i]]
    nodes_v[i] <- as.integer(net$n_nodes %||% 0L)
    edges_v[i] <- as.integer(net$n_edges %||% 0L)
    mat <- net$weights
    if (!is.null(mat)) {
      nz <- if (isTRUE(net$directed)) mat[mat != 0] else
        mat[mat != 0 & row(mat) <= col(mat)]
      if (length(nz) > 0L) {
        wmin_v[i] <- min(nz); wmax_v[i] <- max(nz)
      } else {
        wmin_v[i] <- NA_real_; wmax_v[i] <- NA_real_
      }
    } else {
      wmin_v[i] <- NA_real_; wmax_v[i] <- NA_real_
    }
  }

  weights_str <- ifelse(
    is.na(wmin_v),
    "--",
    sprintf(paste0("[%.", digits, "f, %.", digits, "f]"), wmin_v, wmax_v)
  )

  cat("\n")
  cols <- list(
    Group   = grps,
    Nodes   = sprintf("%d", nodes_v),
    Edges   = sprintf("%d", edges_v),
    Weights = weights_str
  )

  # If a clustering attribute provides assignments, append N + %.
  if (!is.null(cl) && !is.null(cl$assignments)) {
    sizes <- as.integer(tabulate(cl$assignments, nbins = k))
    n_total <- sum(sizes)
    if (n_total > 0L) {
      cols[["N"]] <- .fmt_size_pct(sizes, n_total)
    }
  }

  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  invisible(x)
}


#' Print Method for Multilevel Network Object
#'
#' @param x A \code{netobject_ml}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' set.seed(1)
#' obs <- data.frame(id = rep(1:3, each = 5),
#'                   A = rnorm(15), B = rnorm(15), C = rnorm(15))
#' net_ml <- build_network(obs, method = "cor",
#'                          params = list(id = "id"), level = "both")
#' print(net_ml)
#' \donttest{
#' set.seed(1)
#' obs <- data.frame(
#'   id  = rep(1:5, each = 8),
#'   A   = rnorm(40), B = rnorm(40),
#'   C   = rnorm(40), D = rnorm(40)
#' )
#' net_ml <- build_network(obs, method = "cor",
#'                          params = list(id = "id"), level = "both")
#' print(net_ml)
#' }
#'
#' @export
print.netobject_ml <- function(x, ...) {
  cat(sprintf("Multilevel Network (method: %s)\n", x$method))
  cat("-- Between-person --\n")
  b <- x$between
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", b$n_nodes, b$n_edges))
  if (!is.null(b$n)) {
    cat(sprintf("  Sample size: %d (unique persons)\n", b$n))
  }
  cat("-- Within-person --\n")
  w <- x$within
  cat(sprintf("  Nodes: %d  |  Edges: %d\n", w$n_nodes, w$n_edges))
  if (!is.null(w$n)) {
    cat(sprintf("  Sample size: %d (observations)\n", w$n))
  }
  invisible(x)
}


# ---- Predictability ----

#' Compute Node Predictability
#'
#' @description
#' Computes the proportion of variance explained (R\eqn{^2}) for each node in
#' the network, following Haslbeck & Waldorp (2018).
#'
#' For \code{method = "glasso"} or \code{"pcor"}, predictability is computed
#' analytically from the precision matrix:
#' \deqn{R^2_j = 1 - 1 / \Omega_{jj}}
#' where \eqn{\Omega} is the precision (inverse correlation) matrix.
#'
#' For \code{method = "cor"}, predictability is the multiple R\eqn{^2} from
#' regressing each node on its network neighbors (nodes with non-zero edges).
#'
#' @param object A \code{netobject} or \code{netobject_ml} object.
#' @param data Optional data frame of the original variables used to estimate
#'   the network. Required for \code{method = "cor"} (multiple-R\eqn{^2} regression
#'   of each node on its neighbours); ignored for the precision-matrix path used
#'   by \code{glasso}/\code{pcor}, which has no need of the raw data.
#' @param ... Additional arguments (ignored).
#'
#' @return For \code{netobject}: a named numeric vector of R\eqn{^2} values
#'   (one per node, between 0 and 1).
#'
#'   For \code{netobject_ml}: a list with elements \code{$between} and
#'   \code{$within}, each a named numeric vector.
#'
#' @references
#' Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models
#' predict observations? On the importance of predictability in network models.
#' \emph{Behavior Research Methods}, 50(2), 853--861.
#' \doi{10.3758/s13428-017-0910-x}
#'
#' @examples
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' net <- build_network(as.data.frame(mat), method = "glasso")
#' predictability(net)
#'
#' @export
predictability <- function(object, ...) {
  UseMethod("predictability")
}


#' @rdname predictability
#' @return A named numeric vector of predictability values per node.
#' @export
predictability.netobject <- function(object, data = NULL, ...) {
  labels <- object$nodes$label
  p <- length(labels)

  # ---- R² ----
  if (!is.null(object$precision_matrix)) {
    omega_diag <- diag(object$precision_matrix)
    r2 <- 1 - 1 / omega_diag
  } else if (!is.null(object$cor_matrix)) {
    S   <- object$cor_matrix
    net <- object$weights
    r2 <- vapply(seq_len(p), function(j) {
      neighbors <- which(net[j, ] != 0)
      if (length(neighbors) == 0L) return(0)
      if (length(neighbors) == 1L) return(S[neighbors, j]^2)
      r_vec <- S[neighbors, j]
      R_nn  <- S[neighbors, neighbors]
      tryCatch(as.numeric(crossprod(r_vec, solve(R_nn, r_vec))),
               error = function(e) 0)
    }, numeric(1))
  } else {
    stop("predictability() requires a precision or correlation matrix ",
         "(methods: glasso, pcor, cor). Method '", object$method,
         "' does not support predictability.", call. = FALSE)
  }
  r2 <- pmin(pmax(r2, 0), 1)

  # ---- RMSE ----
  if (is.null(data)) data <- object$data
  rmse <- rep(NA_real_, p)

  if (is.data.frame(data) || is.matrix(data)) {
    X <- tryCatch(as.matrix(data[, labels, drop = FALSE]),
                  error = function(e) NULL)
    if (!is.null(X) && is.numeric(X)) {
      X  <- X[complete.cases(X), , drop = FALSE]
      mu <- colMeans(X)
      sd <- apply(X, 2, stats::sd)
      Z  <- scale(X, center = mu, scale = sd)

      if (!is.null(object$precision_matrix)) {
        Omega <- object$precision_matrix
        rmse <- vapply(seq_len(p), function(j) {
          beta_z <- -Omega[-j, j] / Omega[j, j]
          z_hat  <- Z[, -j, drop = FALSE] %*% beta_z
          y_hat  <- mu[j] + sd[j] * z_hat
          sqrt(mean((X[, j] - y_hat)^2))
        }, numeric(1))
      } else {
        S   <- object$cor_matrix
        net <- object$weights
        rmse <- vapply(seq_len(p), function(j) {
          nbrs <- which(net[j, ] != 0)
          if (length(nbrs) == 0L) return(sd[j])
          beta_z <- tryCatch(solve(S[nbrs, nbrs, drop = FALSE], S[nbrs, j]),
                             error = function(e) NULL)
          if (is.null(beta_z)) return(NA_real_)
          z_hat <- Z[, nbrs, drop = FALSE] %*% beta_z
          y_hat <- mu[j] + sd[j] * z_hat
          sqrt(mean((X[, j] - y_hat)^2))
        }, numeric(1))
      }
    }
  }

  data.frame(node = labels, R2 = r2, RMSE = rmse,
             stringsAsFactors = FALSE, row.names = NULL)
}


#' @rdname predictability
#' @return A list with \code{within} and \code{between} predictability vectors.
#' @export
predictability.netobject_ml <- function(object, ...) {
  list(
    between = predictability(object$between),
    within  = predictability(object$within)
  )
}


#' @rdname predictability
#' @return A named list of per-group predictability vectors.
#' @export
predictability.netobject_group <- function(object, ...) {
  lapply(object, predictability)
}
