# ---- Bootstrap Network Estimation ----

#' Bootstrap a Network Estimate
#'
#' @description
#' Non-parametric bootstrap for any network estimated by
#' \code{\link{build_network}}. Works with all built-in methods
#' (transition and association) as well as custom registered estimators.
#'
#' For transition methods (\code{"relative"}, \code{"frequency"},
#' \code{"co_occurrence"}), uses a fast pre-computation strategy:
#' per-sequence count matrices are computed once, and each bootstrap
#' iteration only resamples sequences via \code{colSums} (C-level)
#' plus lightweight post-processing. Data must be in wide format for
#' transition bootstrap; use \code{\link{convert_sequence_format}} to
#' convert long-format data first.
#'
#' For association methods (\code{"cor"}, \code{"pcor"}, \code{"glasso"},
#' and custom estimators), the full estimator is called on resampled rows
#' each iteration.
#'
#' If a transition network contains only one sequence, the function warns that
#' such a network is not recommended for bootstrap or other confirmatory
#' testing.
#'
#' @param x A \code{netobject} from \code{\link{build_network}}.
#'   The data, method, params, scaling, threshold, and level are all
#'   extracted from this object.
#' @param iter Integer. Number of bootstrap iterations (default: 1000).
#' @param ci_level Numeric. Significance level for CIs and p-values
#'   (default: 0.05).
#' @param inference Character. \code{"stability"} (default) tests whether
#'   bootstrap replicates fall within a multiplicative consistency range
#'   around the original weight. \code{"threshold"} tests whether
#'   replicates exceed a fixed edge threshold.
#' @param consistency_range Numeric vector of length 2. Multiplicative
#'   bounds for stability inference (default: \code{c(0.75, 1.25)}).
#' @param edge_threshold Numeric or NULL. Fixed threshold for
#'   \code{inference = "threshold"}. If NULL, defaults to the 10th
#'   percentile of absolute original edge weights.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#' @param boundary Character. Comparison rule when computing the consistency-range
#'   p-value. \code{"inclusive"} (default, tna-compatible) counts iterations
#'   that meet the bound (\eqn{\le} / \eqn{\ge}); \code{"strict"} counts only
#'   iterations strictly outside (\eqn{<} / \eqn{>}).
#'
#' @return An object of class \code{"net_bootstrap"} containing:
#' \describe{
#'   \item{original}{The original \code{netobject}.}
#'   \item{mean}{Bootstrap mean weight matrix.}
#'   \item{sd}{Bootstrap SD matrix.}
#'   \item{p_values}{P-value matrix.}
#'   \item{significant}{Original weights where p < ci_level, else 0.}
#'   \item{ci_lower}{Lower CI bound matrix.}
#'   \item{ci_upper}{Upper CI bound matrix.}
#'   \item{cr_lower}{Consistency range lower bound (stability only).}
#'   \item{cr_upper}{Consistency range upper bound (stability only).}
#'   \item{summary}{Long-format data frame of edge-level statistics.}
#'   \item{model}{Pruned \code{netobject} (non-significant edges zeroed).}
#'   \item{method, params, iter, ci_level, inference}{Bootstrap config.}
#'   \item{consistency_range, edge_threshold}{Inference parameters.}
#' }
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
#'   method = "relative")
#' boot <- bootstrap_network(net, iter = 10)
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 30, TRUE), V2 = sample(LETTERS[1:4], 30, TRUE),
#'   V3 = sample(LETTERS[1:4], 30, TRUE), V4 = sample(LETTERS[1:4], 30, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' boot <- bootstrap_network(net, iter = 100)
#' print(boot)
#' summary(boot)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{print.net_bootstrap}},
#'   \code{\link{summary.net_bootstrap}}
#'
#' @importFrom stats quantile sd
#' @export
bootstrap_network <- function(x,
                              iter = 1000L,
                              ci_level = 0.05,
                              inference = "stability",
                              consistency_range = c(0.75, 1.25),
                              edge_threshold = NULL,
                              seed = NULL,
                              boundary = c("inclusive", "strict")) {
  boundary <- match.arg(boundary)

  # ---- wtna_mixed dispatch: bootstrap both components ----
  if (inherits(x, "wtna_mixed")) {
    if (!is.null(seed)) set.seed(seed)
    result <- list(
      transition   = bootstrap_network(x$transition,   iter = iter,
                                       ci_level = ci_level, inference = inference,
                                       consistency_range = consistency_range,
                                       edge_threshold = edge_threshold,
                                       boundary = boundary),
      cooccurrence = bootstrap_network(x$cooccurrence, iter = iter,
                                       ci_level = ci_level, inference = inference,
                                       consistency_range = consistency_range,
                                       edge_threshold = edge_threshold,
                                       boundary = boundary)
    )
    class(result) <- "wtna_boot_mixed"
    return(result)
  }

  # ---- mcml dispatch: convert to netobject_group via as_tna ----
  if (inherits(x, "mcml")) {
    x <- as_tna(x)
  }

  # ---- netobject_group dispatch: bootstrap each element ----
  if (inherits(x, "netobject_group")) {
    results <- lapply(x, function(net) {
      bootstrap_network(net, iter = iter, ci_level = ci_level,
                        inference = inference,
                        consistency_range = consistency_range,
                        edge_threshold = edge_threshold, seed = seed,
                        boundary = boundary)
    })
    class(results) <- c("net_bootstrap_group", "list")
    return(results)
  }

  # Accept netobject or cograph_network
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (!inherits(x, "netobject")) {
    stop("'x' must be a netobject from build_network().", call. = FALSE) # nocov
  }
  if (is.null(x$data)) {
    stop("netobject does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  # Note: edgelist-derived data has no actor / session grouping, so
  # row-resampling uses individual transitions as the unit -- which may
  # under-represent within-actor correlation. Case-dropping and actor-
  # level sequence reconstruction are both more robust alternatives.
  if (identical(attr(x$data, "source"), "edgelist")) {
    warning(
      "Bootstrapping a network built from edgelist input (no per-actor ",
      "sequence data). Each transition is treated as an independent case, ",
      "which has limitations when transitions are correlated within ",
      "actors or sessions.\n",
      "Better alternatives:\n",
      "  * Rebuild as a sequence-based network: build_network(data, ",
      "method = \"tna\", actor = <...>, action = <...>, time = <...>) ",
      "then feed the result to build_mcml(). This reconstructs per-actor ",
      "sequences (via prepare()) so bootstrap and case-dropping operate at ",
      "the actor level natively.\n",
      "  * Use case-dropping rather than bootstrap:\n",
      "      - casedrop_reliability(x)  - model-level edge-weight reliability\n",
      "      - centrality_stability(x)  - per-centrality CS coefficient\n",
      "  * permutation_test(x, y) for comparing two networks.",
      call. = FALSE
    )
  }

  data <- x$data
  method <- x$method
  params <- x$params
  scaling <- x$scaling
  threshold <- x$threshold
  level <- x$level

  # ---- Input validation ----
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(ci_level), length(ci_level) == 1, ci_level > 0, ci_level < 1,
    is.character(inference), length(inference) == 1,
    is.numeric(consistency_range), length(consistency_range) == 2
  )
  iter <- as.integer(iter)
  inference <- match.arg(inference, c("stability", "threshold"))

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Resolve method and get estimator ----
  method <- .resolve_method_alias(method)
  estimator <- get_estimator(method)
  directed <- estimator$directed

  # ---- Original network is x itself ----
  original <- x
  states <- original$nodes$label
  n_states <- length(states)

  # ---- Dispatch bootstrap ----
  # Fast precomputed path only works when data has usable columns
  has_data <- is.data.frame(data) && ncol(data) > 0L && nrow(data) > 0L
  if (method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data) {
      stop("Bootstrap requires the original data stored in the netobject. ",
           "For wtna/cna networks, use wtna() directly instead of ",
           "build_network(method='cna').", call. = FALSE)
    }
    transition_data <- .resampling_transition_data(data, x$metadata, params)
    .warn_single_sequence_confirmatory_network(transition_data, params)
    boot_matrices <- .bootstrap_transition(
      data = transition_data, method = method, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter
    )
  } else {
    boot_matrices <- .bootstrap_association(
      data = data, estimator = estimator, params = params, states = states,
      scaling = scaling, threshold = threshold, iter = iter,
      level = level, id_col = params$id
    )
  }

  # ---- Auto edge_threshold for threshold inference ----
  if (inference == "threshold" && is.null(edge_threshold)) {
    abs_weights <- abs(as.vector(original$weights))
    nz_weights <- abs_weights[abs_weights > 0]
    edge_threshold <- if (length(nz_weights) > 0) {
      quantile(nz_weights, probs = 0.10)
    } else {
      0
    }
  }

  # ---- Compute statistics ----
  stats <- .compute_bootstrap_stats(
    boot_matrices = boot_matrices,
    original_matrix = original$weights,
    states = states,
    directed = directed,
    iter = iter,
    ci_level = ci_level,
    inference = inference,
    consistency_range = consistency_range,
    edge_threshold = edge_threshold,
    boundary = boundary
  )

  # ---- Build summary data frame ----
  summary_df <- .build_bootstrap_summary(
    stats = stats,
    original_matrix = original$weights,
    states = states,
    directed = directed,
    ci_level = ci_level,
    inference = inference
  )

  # ---- Build pruned model ----
  pruned_matrix <- stats$significant
  pruned_edges <- .extract_edges_from_matrix(pruned_matrix, directed = directed)

  model <- list(
    weights = pruned_matrix,
    nodes = original$nodes,
    edges = pruned_edges,
    directed = directed,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    n_nodes = n_states,
    n_edges = nrow(pruned_edges),
    level = level,
    meta = list(source = "nestimate", layout = NULL,
                tna = list(method = method)),
    node_groups = NULL
  )
  class(model) <- c("netobject", "cograph_network")

  # ---- Assemble result ----
  result <- list(
    original          = original,
    mean              = stats$mean,
    sd                = stats$sd,
    p_values          = stats$p_values,
    significant       = stats$significant,
    ci_lower          = stats$ci_lower,
    ci_upper          = stats$ci_upper,
    cr_lower          = stats$cr_lower,
    cr_upper          = stats$cr_upper,
    summary           = summary_df,
    model             = model,
    method            = method,
    params            = params,
    iter              = iter,
    ci_level          = ci_level,
    inference         = inference,
    consistency_range = consistency_range,
    edge_threshold    = edge_threshold
  )
  class(result) <- "net_bootstrap"
  result
}


#' Reattach grouping columns needed by transition resampling
#' @noRd
.resampling_transition_data <- function(data, metadata, params) {
  if (!is.data.frame(data) || !is.data.frame(metadata)) return(data)
  needed <- unique(c(params$actor, params$id, params$id_col))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  add <- setdiff(intersect(needed, names(metadata)), names(data))
  if (length(add) == 0L) return(data)
  cbind(metadata[, add, drop = FALSE], data)
}


#' Warn when confirmatory testing is requested for one sequence
#' @noRd
.warn_single_sequence_confirmatory_network <- function(data, params,
                                                       label = NULL) {
  n_seq <- .transition_resampling_n_sequences(data, params)
  if (!is.na(n_seq) && n_seq <= 1L) {
    prefix <- if (is.null(label)) "" else paste0(label, ": ")
    warning(
      prefix,
      "A network with one long sequence is not recommended and can't be ",
      "validated using bootstrap and other confirmatory testings.",
      call. = FALSE
    )
  }
  invisible(n_seq)
}


#' @noRd
.transition_resampling_n_sequences <- function(data, params) {
  if (!is.data.frame(data)) {
    return(NA_integer_)
  }
  actor <- .param_get(params, "actor") %||% .param_get(params, "id") %||%
    .param_get(params, "id_col")
  if (!is.null(actor) && all(actor %in% names(data))) {
    if (length(actor) == 1L) {
      return(length(unique(data[[actor]])))
    }
    return(length(unique(interaction(data[, actor, drop = FALSE],
                                     drop = TRUE))))
  }
  nrow(data)
}


# ---- Transition fast path ----

#' Bootstrap transition networks via pre-computed per-sequence counts
#' @noRd
.bootstrap_transition <- function(data, method, params, states,
                                  scaling, threshold, iter) {
  n_states <- length(states)
  nbins <- n_states * n_states
  is_relative <- method == "relative"

  # Pre-compute per-sequence count matrix ONCE
  trans_2d <- .precompute_per_sequence(data, method, params, states)
  n_seq <- nrow(trans_2d)

  # vapply: each iteration resamples sequences + sums + post-processes
  boot_flat <- vapply(seq_len(iter), function(i) {
    idx <- sample.int(n_seq, n_seq, replace = TRUE)
    boot_counts <- colSums(trans_2d[idx, , drop = FALSE])
    # Inline post-processing (no dimnames -- we flatten immediately)
    mat <- matrix(boot_counts, n_states, n_states, byrow = TRUE)
    if (is_relative) {
      rs <- rowSums(mat)
      nz <- rs > 0
      mat[nz, ] <- mat[nz, ] / rs[nz]
    }
    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling)
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat)
  }, numeric(nbins))

  t(boot_flat)
}


#' Pre-compute per-sequence count matrices
#'
#' Returns an n_seq x (n_states^2) matrix where each row holds the
#' flattened transition (or co-occurrence) counts for one sequence.
#'
#' Only wide-format data is supported for the bootstrap fast path.
#' For long-format data, convert first using
#' \code{\link{convert_sequence_format}}.
#' @noRd
.precompute_per_sequence <- function(data, method, params, states) {
  # Determine format (mirrors estimator dispatch)
  format <- .param_get(params, "format", "auto")
  action <- .param_get(params, "action", "Action")

  if (format == "auto") {
    format <- if (action %in% names(data)) "long" else "wide"
  }

  if (format == "long") {
    stop(
      "Bootstrap fast path requires wide-format data. ",
      "Convert with convert_sequence_format() first.",
      call. = FALSE
    )
  }

  id_col <- .param_get(params, "id") %||% .param_get(params, "id_col")
  codes <- .param_get(params, "codes")
  if (identical(format, "onehot") || !is.null(codes)) {
    return(.precompute_per_sequence_onehot(data, method, params, states))
  }
  cols <- params$cols
  .precompute_per_sequence_wide(data, method, cols, id_col, states)
}


#' Pre-compute per-sequence counts from one-hot format
#' @noRd
.precompute_per_sequence_onehot <- function(data, method, params, states) {
  codes <- .param_get(params, "codes", states)
  actor <- .param_get(params, "actor")
  window_size <- .param_get(params, "window_size", 1L)
  mode <- .param_get(params, "mode", "non-overlapping")
  wtna_method <- if (method == "co_occurrence") "cooccurrence" else "transition"
  n_states <- length(states)
  nbins <- n_states * n_states

  groups <- if (is.null(actor)) {
    list(data)
  } else if (length(actor) == 1L) {
    split(data, data[[actor]], drop = TRUE)
  } else {
    split(data, interaction(data[, actor, drop = FALSE], drop = TRUE),
          drop = TRUE)
  }

  rows <- lapply(groups, function(g) {
    mat <- .estimator_wtna_core(
      g, codes = codes, window_size = window_size, mode = mode,
      actor = NULL, wtna_method = wtna_method, type = "frequency"
    )$matrix
    mat <- mat[states, states, drop = FALSE]
    as.vector(t(mat))
  })

  if (length(rows) == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = nbins))
  }
  do.call(rbind, rows)
}


#' Pre-compute per-sequence counts from wide format (transition + co-occurrence)
#'
#' Single function handling both consecutive-pair transitions
#' (relative/frequency) and all-column-pair co-occurrences.
#' @noRd
.precompute_per_sequence_wide <- function(data, method, cols, id_col, states) {
  state_cols <- .select_state_cols(data, id_col, cols)
  mat <- as.matrix(data[, state_cols, drop = FALSE])
  nc <- ncol(mat)
  nr <- nrow(mat)
  n_states <- length(states)
  nbins <- n_states * n_states

  # Integer-encode the entire matrix once
  int_mat <- matrix(match(mat, states), nrow = nr, ncol = nc)

  if (method %in% c("relative", "frequency")) {
    # Consecutive column pairs: (col[i], col[i+1])
    from_mat <- int_mat[, -nc, drop = FALSE]
    to_mat <- int_mat[, -1L, drop = FALSE]

    row_ids <- rep(seq_len(nr), times = nc - 1L)
    from_vec <- as.vector(from_mat)
    to_vec <- as.vector(to_mat)

    valid <- !is.na(from_vec) & !is.na(to_vec)
    row_ids <- row_ids[valid]
    from_vec <- from_vec[valid]
    to_vec <- to_vec[valid]

    pair_idx <- (from_vec - 1L) * n_states + to_vec
    combined_idx <- (row_ids - 1L) * nbins + pair_idx

    counts <- tabulate(combined_idx, nbins = nr * nbins)
    matrix(as.numeric(counts), nrow = nr, ncol = nbins, byrow = TRUE)

  } else {
    # Co-occurrence: all column pairs (i, j) where i < j
    result <- matrix(0, nrow = nr, ncol = nbins)

    for (i in seq_len(nc - 1L)) {
      col_i <- int_mat[, i]
      for (j in seq(i + 1L, nc)) {
        col_j <- int_mat[, j]
        valid <- !is.na(col_i) & !is.na(col_j)
        fi <- col_i[valid]
        tj <- col_j[valid]
        rows_valid <- which(valid)

        # Bidirectional: forward + reverse
        idx_fwd <- (fi - 1L) * n_states + tj
        idx_rev <- (tj - 1L) * n_states + fi
        combined <- c(
          (rows_valid - 1L) * nbins + idx_fwd,
          (rows_valid - 1L) * nbins + idx_rev
        )

        tab <- tabulate(combined, nbins = nr * nbins)
        result <- result + matrix(tab, nrow = nr, ncol = nbins, byrow = TRUE)
      }
    }

    # Self-pairs are double-counted by bidirectional approach
    diag_indices <- (seq_len(n_states) - 1L) * n_states + seq_len(n_states)
    result[, diag_indices] <- result[, diag_indices] / 2

    result
  }
}


# ---- Association path ----

#' Bootstrap association networks via full estimator calls
#' @noRd
.bootstrap_association <- function(data, estimator, params, states,
                                   scaling, threshold, iter,
                                   level, id_col) {
  n <- nrow(data)
  n_states <- length(states)
  nbins <- n_states * n_states

  boot_flat <- vapply(seq_len(iter), function(i) {
    boot_data <- data[sample.int(n, n, replace = TRUE), , drop = FALSE]

    # For multilevel: apply decomposition to bootstrap sample
    if (!is.null(level) && !is.null(id_col) && !estimator$directed) {
      boot_data <- tryCatch( # nocov start
        .decompose_multilevel(boot_data, id_col = id_col, level = level),
        error = function(e) NULL
      )
      if (is.null(boot_data)) return(rep(NA_real_, nbins)) # nocov end
    }

    est <- tryCatch(
      do.call(estimator$fn, c(list(data = boot_data), params)),
      error = function(e) NULL
    )
    if (is.null(est)) return(rep(NA_real_, nbins))

    mat <- est$matrix
    # Align to expected states order
    if (!identical(rownames(mat), states)) {
      common <- intersect(states, rownames(mat))
      if (length(common) < n_states) return(rep(NA_real_, nbins))
      mat <- mat[states, states] # nocov
    }

    if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov start
    if (threshold > 0) mat[abs(mat) < threshold] <- 0
    as.vector(mat) # nocov end
  }, numeric(nbins))

  t(boot_flat)
}


# ---- Vectorized statistics ----

#' Compute bootstrap statistics from the iter x nbins matrix
#' @noRd
.compute_bootstrap_stats <- function(boot_matrices, original_matrix, states,
                                     directed, iter, ci_level, inference,
                                     consistency_range, edge_threshold,
                                     boundary = "inclusive") {
  n_states <- length(states)
  orig_flat <- as.vector(original_matrix)

  # Mean and SD (vectorized)
  weights_mean <- colMeans(boot_matrices, na.rm = TRUE)
  weights_sd <- apply(boot_matrices, 2, sd, na.rm = TRUE)

  # Percentile CIs (vectorized via apply)
  ci_lower <- apply(boot_matrices, 2, quantile,
                    probs = ci_level / 2, na.rm = TRUE)
  ci_upper <- apply(boot_matrices, 2, quantile,
                    probs = 1 - ci_level / 2, na.rm = TRUE)

  # p-values via vectorized sweep
  valid_rows <- rowSums(is.na(boot_matrices)) == 0
  bm <- boot_matrices[valid_rows, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 1) {
    p_values <- rep(1, length(orig_flat))
  } else if (inference == "stability") {
    cr_low <- pmin(orig_flat * consistency_range[1],
                   orig_flat * consistency_range[2])
    cr_high <- pmax(orig_flat * consistency_range[1],
                    orig_flat * consistency_range[2])
    # boundary = "inclusive" (default, tna-compatible): <=, >= count an
    # iteration that lands EXACTLY on the band edge as outside the band.
    # boundary = "strict": <, > count only iterations strictly outside.
    cmp_lo <- if (identical(boundary, "strict")) "<"  else "<="
    cmp_hi <- if (identical(boundary, "strict")) ">"  else ">="
    below <- sweep(bm, 2, cr_low,  cmp_lo)
    above <- sweep(bm, 2, cr_high, cmp_hi)
    p_counts <- colSums(below | above)
    p_values <- (p_counts + 1) / (n_valid + 1)
  } else {
    below_thresh <- sweep(abs(bm), 2, edge_threshold, "<")
    p_counts <- colSums(below_thresh)
    p_values <- (p_counts + 1) / (n_valid + 1)
  }

  # Reshape all to n_states x n_states matrices
  to_mat <- function(x) {
    matrix(x, n_states, n_states, dimnames = list(states, states))
  }

  p_mat <- to_mat(p_values)
  sig_mask <- (p_mat < ci_level) * 1
  sig_mat <- original_matrix * sig_mask

  list(
    mean        = to_mat(weights_mean),
    sd          = to_mat(weights_sd),
    ci_lower    = to_mat(ci_lower),
    ci_upper    = to_mat(ci_upper),
    p_values    = p_mat,
    significant = sig_mat,
    cr_lower    = to_mat(pmin(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2])),
    cr_upper    = to_mat(pmax(orig_flat * consistency_range[1],
                              orig_flat * consistency_range[2]))
  )
}


#' Build long-format summary data frame from bootstrap stats
#' @noRd
.build_bootstrap_summary <- function(stats, original_matrix, states, directed,
                                     ci_level, inference) {
  n <- length(states)
  dt <- data.table::data.table(
    from     = rep(states, each = n),
    to       = rep(states, times = n),
    weight   = as.vector(t(original_matrix)),
    mean     = as.vector(t(stats$mean)),
    sd       = as.vector(t(stats$sd)),
    p_value  = as.vector(t(stats$p_values)),
    sig      = as.vector(t(stats$p_values)) < ci_level,
    ci_lower = as.vector(t(stats$ci_lower)),
    ci_upper = as.vector(t(stats$ci_upper))
  )
  if (inference == "stability") {
    dt[, cr_lower := as.vector(t(stats$cr_lower))]
    dt[, cr_upper := as.vector(t(stats$cr_upper))]
  }

  # Filter: keep only non-zero original edges (including self-loops)
  if (directed) {
    dt <- dt[weight != 0]
  } else {
    dt <- dt[weight != 0 & from <= to]
  }

  as.data.frame(dt)
}


# ---- S3 Methods ----

#' Print Method for net_bootstrap
#'
#' @param x A \code{net_bootstrap} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
#'   method = "relative")
#' boot <- bootstrap_network(net, iter = 10)
#' print(boot)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
#'   V3 = c("C","A","C","B","A")
#' )
#' net  <- build_network(seqs, method = "relative")
#' boot <- bootstrap_network(net, iter = 20)
#' print(boot)
#' }
#'
#' @export
print.net_bootstrap <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative)",
    frequency     = "Transition Network (frequency)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation (EBICglasso)",
    pcor          = "Partial Correlation",
    cor           = "Correlation Network",
    attention     = "Attention Network",
    wtna          = "Window TNA"
  )
  lbl_raw   <- method_labels[x$method]
  label     <- if (is.na(lbl_raw)) sprintf("Network (%s)", x$method) else unname(lbl_raw)
  dir_label <- if (x$original$directed) "directed" else "undirected"
  ci_pct    <- sprintf("%d%%", round((1 - x$ci_level) * 100))
  n_orig    <- x$original$n_edges
  n_sig     <- x$model$n_edges
  n_nonsig  <- n_orig - n_sig

  # Top significant edges first
  s <- x$summary
  if (!is.null(s) && nrow(s) > 0L) {
    sig_s <- s[s$sig, , drop = FALSE]
    if (nrow(sig_s) > 0L) {
      sig_s <- sig_s[order(abs(sig_s$mean), decreasing = TRUE), , drop = FALSE]
      top   <- head(sig_s, 5L)
      cat("  Edge                   Mean     95% CI          p\n")
      cat("  -----------------------------------------------\n")
      for (i in seq_len(nrow(top))) {
        r <- top[i, ]
        lbl <- if (!is.null(r$from) && !is.null(r$to))
          sprintf("%-20s", paste0(r$from, " \u2192 ", r$to))
        else
          sprintf("%-20s", r$edge)
        stars <- if (r$p_value < 0.001) "***" else if (r$p_value < 0.01) "** " else "*  "
        cat(sprintf("  %s  %6.3f  [%5.3f, %5.3f]  %s\n",
                    lbl, r$mean, r$ci_lower, r$ci_upper, stars))
      }
      if (nrow(sig_s) > 5L)
        cat(sprintf("  ... and %d more significant edges\n", nrow(sig_s) - 5L))
      cat("\n")
    }
  }

  cat(sprintf("Bootstrap Network  [%s | %s]\n", label, dir_label))
  cat(sprintf("  Iterations : %d  |  Nodes : %d\n", x$iter, x$original$n_nodes))
  cat(sprintf("  Edges      : %d significant / %d total\n", n_sig, n_orig))
  cat(sprintf("  CI         : %s  |  Inference: %s", ci_pct, x$inference))
  if (x$inference == "stability") {
    cat(sprintf("  |  CR [%.2f, %.2f]", x$consistency_range[1L], x$consistency_range[2L]))
  } else {
    cat(sprintf("  |  Threshold: %g", x$edge_threshold))
  }
  cat("\n")

  invisible(x)
}


#' Summary Method for net_bootstrap
#'
#' @param object A \code{net_bootstrap} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with edge-level bootstrap statistics.
#'
#' @examples
#' net <- build_network(data.frame(V1 = c("A","B","C"), V2 = c("B","C","A")),
#'   method = "relative")
#' boot <- bootstrap_network(net, iter = 10)
#' summary(boot)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B"), V2 = c("B","C","B","A","C"),
#'   V3 = c("C","A","C","B","A")
#' )
#' net  <- build_network(seqs, method = "relative")
#' boot <- bootstrap_network(net, iter = 20)
#' summary(boot)
#' }
#'
#' @export
summary.net_bootstrap <- function(object, ...) {
  object$summary
}

#' Print Method for net_bootstrap_group
#' @param x A \code{net_bootstrap_group} object.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","C","A"),
#'   V3 = c("C","A","B","B"), grp = c("X","X","Y","Y"))
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' boot <- bootstrap_network(nets, iter = 10)
#' print(boot)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B","A"),
#'   V2 = c("B","C","B","A","C","B"),
#'   V3 = c("C","A","C","B","A","C"),
#'   grp = c("X","X","X","Y","Y","Y")
#' )
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' boot <- bootstrap_network(nets, iter = 20)
#' print(boot)
#' }
#' @export
print.net_bootstrap_group <- function(x, ...) {
  grp_names <- names(x)
  iter      <- x[[1L]]$iter
  ci_pct    <- sprintf("%d%%", round((1 - x[[1L]]$ci_level) * 100))

  # Per-group edge counts
  grp_stats <- vapply(grp_names, function(nm) {
    b <- x[[nm]]
    c(total = b$original$n_edges, sig = b$model$n_edges)
  }, numeric(2L))

  # Helper: make edge keys "from->to" from a summary df
  .edge_keys <- function(s) {
    if (is.null(s) || nrow(s) == 0L) return(character(0L))
    idx <- which(s$sig)
    if (length(idx) == 0L) return(character(0L))
    paste0(s$from[idx], "\u2192", s$to[idx])
  }

  # Shared significant edges (present in all groups)
  sig_keys  <- lapply(grp_names, function(nm) .edge_keys(x[[nm]]$summary))
  shared    <- Reduce(intersect, sig_keys)

  # Detect disjoint state spaces (e.g. mcml: macro carries cluster
  # names, each within-cluster carries different code names). When the
  # union of every group's FULL edge set has empty pairwise overlap,
  # the "shared edges" metric is meaningless -- suppress its line.
  full_keys <- lapply(grp_names, function(nm) {
    s <- x[[nm]]$summary
    if (is.null(s) || nrow(s) == 0L) character(0L)
    else paste0(s$from, "->", s$to)
  })
  has_shared_space <- length(Reduce(intersect, full_keys)) > 0L

  # Top shared edges table first
  if (length(shared) > 0L) {
    .mean_for_key <- function(nm, key) {
      s <- x[[nm]]$summary
      parts <- strsplit(key, "\u2192", fixed = TRUE)[[1L]]
      idx <- which(s$from == parts[1L] & s$to == parts[2L])
      if (length(idx) == 0L) NA_real_ else s$mean[idx[1L]]
    }
    shared_max <- vapply(shared, function(k) {
      max(abs(vapply(grp_names, .mean_for_key, numeric(1L), key = k)), na.rm = TRUE)
    }, numeric(1L))
    top_shared <- shared[order(shared_max, decreasing = TRUE)][seq_len(min(5L, length(shared)))]
    hdr <- paste(sprintf("%-8s", grp_names), collapse = "  ")
    cat(sprintf("  Edge                   %s\n", hdr))
    cat(sprintf("  %s\n", strrep("-", 24L + 10L * length(grp_names))))
    for (k in top_shared) {
      grp_vals <- vapply(grp_names, .mean_for_key, numeric(1L), key = k)
      vals_str <- paste(sprintf("%-8.3f", grp_vals), collapse = "  ")
      cat(sprintf("  %-22s  %s\n", k, vals_str))
    }
    if (length(shared) > 5L)
      cat(sprintf("  ... and %d more shared significant edges\n", length(shared) - 5L))
    cat("\n")
  }

  cat(sprintf("Grouped Bootstrap  [%d groups | %d iterations | %s CI]\n",
              length(grp_names), iter, ci_pct))
  for (nm in grp_names) {
    cat(sprintf("  %-20s  %d sig / %d total\n",
                nm, grp_stats["sig", nm], grp_stats["total", nm]))
  }
  if (has_shared_space) {
    cat(sprintf("  Shared (all groups)   %d edges\n", length(shared)))
  }

  invisible(x)
}

#' Summary Method for net_bootstrap_group
#' @param object A \code{net_bootstrap_group} object.
#' @param ... Ignored.
#' @return A data frame with group, edge, and bootstrap statistics columns.
#' @examples
#' seqs <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","C","A"),
#'   V3 = c("C","A","B","B"), grp = c("X","X","Y","Y"))
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' boot <- bootstrap_network(nets, iter = 10)
#' summary(boot)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = c("A","B","A","C","B","A"),
#'   V2 = c("B","C","B","A","C","B"),
#'   V3 = c("C","A","C","B","A","C"),
#'   grp = c("X","X","X","Y","Y","Y")
#' )
#' nets <- build_network(seqs, method = "relative", group = "grp")
#' boot <- bootstrap_network(nets, iter = 20)
#' summary(boot)
#' }
#' @export
summary.net_bootstrap_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df       <- object[[nm]]$summary
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' Print Method for wtna_boot_mixed
#'
#' @param x A \code{wtna_boot_mixed} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
#' mixed <- wtna(oh, method = "both")
#' boot  <- bootstrap_network(mixed, iter = 10)
#' print(boot)
#' \donttest{
#' set.seed(1)
#' oh <- data.frame(
#'   A = c(1,0,1,0,1,0,1,0),
#'   B = c(0,1,0,1,0,1,0,1),
#'   C = c(1,1,0,0,1,1,0,0)
#' )
#' mixed <- wtna(oh, method = "both")
#' boot  <- bootstrap_network(mixed, iter = 20)
#' print(boot)
#' }
#'
#' @export
print.wtna_boot_mixed <- function(x, ...) {
  cat("Mixed Window TNA Bootstrap\n")
  cat("-- Transition --\n")
  print(x$transition)
  cat("-- Co-occurrence --\n")
  print(x$cooccurrence)
  invisible(x)
}


#' Summary Method for wtna_boot_mixed
#'
#' @param object A \code{wtna_boot_mixed} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A list with \code{$transition} and \code{$cooccurrence} summary data frames.
#'
#' @examples
#' oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
#' mixed <- wtna(oh, method = "both")
#' boot  <- bootstrap_network(mixed, iter = 10)
#' summary(boot)
#' \donttest{
#' set.seed(1)
#' oh <- data.frame(
#'   A = c(1,0,1,0,1,0,1,0),
#'   B = c(0,1,0,1,0,1,0,1),
#'   C = c(1,1,0,0,1,1,0,0)
#' )
#' mixed <- wtna(oh, method = "both")
#' boot  <- bootstrap_network(mixed, iter = 20)
#' summary(boot)
#' }
#'
#' @export
summary.wtna_boot_mixed <- function(object, ...) {
  list(
    transition   = summary(object$transition),
    cooccurrence = summary(object$cooccurrence)
  )
}
