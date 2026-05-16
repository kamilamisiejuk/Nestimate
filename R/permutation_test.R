# ---- Permutation Test for Network Comparison ----

#' Permutation Test for Network Comparison
#'
#' @description
#' Compares two networks estimated by \code{\link{build_network}} using a
#' permutation test. Works with all built-in methods (transition and
#' association) as well as custom registered estimators. The test shuffles
#' which observations belong to which group, re-estimates networks, and tests
#' whether observed edge-wise differences exceed chance.
#'
#' For transition methods (\code{"relative"}, \code{"frequency"},
#' \code{"co_occurrence"}), uses a fast pre-computation strategy: per-sequence
#' count matrices are computed once, and each permutation iteration only
#' shuffles group labels and computes group-wise \code{colSums}.
#'
#' For association methods (\code{"cor"}, \code{"pcor"}, \code{"glasso"},
#' and custom estimators), the full estimator is called on each permuted
#' group split.
#'
#' If either transition network contains only one sequence, the function warns
#' that such a network is not recommended for permutation or other
#' confirmatory testing.
#'
#' @param x A \code{netobject} (from \code{\link{build_network}}).
#' @param y A \code{netobject} (from \code{\link{build_network}}).
#'   Must use the same method and have the same nodes as \code{x}.
#' @param iter Integer. Number of permutation iterations (default: 1000).
#' @param alpha Numeric. Significance level (default: 0.05).
#' @param paired Logical. If \code{TRUE}, permute within pairs (requires
#'   equal number of observations in \code{x} and \code{y}). Default: FALSE.
#' @param adjust Character. p-value adjustment method passed to
#'   \code{\link[stats]{p.adjust}} (default: \code{"none"}). Common choices:
#'   \code{"holm"}, \code{"BH"}, \code{"bonferroni"}.
#' @param nlambda Integer. Number of lambda values for the \code{glassopath}
#'   regularisation path (only used when \code{method = "glasso"}).
#'   Higher values give finer lambda resolution at the cost of speed.
#'   Default: 50.
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"net_permutation"} containing:
#' \describe{
#'   \item{x}{The first \code{netobject}.}
#'   \item{y}{The second \code{netobject}.}
#'   \item{diff}{Observed difference matrix (\code{x - y}).}
#'   \item{diff_sig}{Observed difference where \code{p < alpha}, else 0.}
#'   \item{p_values}{P-value matrix (adjusted if \code{adjust != "none"}).}
#'   \item{effect_size}{Effect size matrix (observed diff / SD of permutation diffs).}
#'   \item{summary}{Long-format data frame of edge-level results.}
#'   \item{method}{The network estimation method.}
#'   \item{iter}{Number of permutation iterations.}
#'   \item{alpha}{Significance level used.}
#'   \item{paired}{Whether paired permutation was used.}
#'   \item{adjust}{p-value adjustment method used.}
#' }
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
#'                  V2 = sample(LETTERS[1:4], 20, TRUE),
#'                  V3 = sample(LETTERS[1:4], 20, TRUE))
#' d2 <- data.frame(V1 = sample(LETTERS[1:4], 20, TRUE),
#'                  V2 = sample(LETTERS[1:4], 20, TRUE),
#'                  V3 = sample(LETTERS[1:4], 20, TRUE))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 100, seed = 42)
#' print(perm)
#' summary(perm)
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}},
#'   \code{\link{print.net_permutation}},
#'   \code{\link{summary.net_permutation}}
#'
#' @importFrom stats p.adjust sd
#' @export
permutation <- function(x, y = NULL,
                             iter = 1000L,
                             alpha = 0.05,
                             paired = FALSE,
                             adjust = "none",
                             nlambda = 50L,
                             seed = NULL) {

  # ---- wtna_mixed dispatch: permute both components ----
  if (inherits(x, "wtna_mixed") || inherits(y, "wtna_mixed")) {
    if (!inherits(x, "wtna_mixed") || !inherits(y, "wtna_mixed")) {
      stop("Both x and y must be wtna_mixed objects.", call. = FALSE)
    }
    result <- list(
      transition = permutation(
        x$transition, y$transition, iter = iter, alpha = alpha,
        paired = paired, adjust = adjust, nlambda = nlambda, seed = seed
      ),
      cooccurrence = permutation(
        x$cooccurrence, y$cooccurrence, iter = iter, alpha = alpha,
        paired = paired, adjust = adjust, nlambda = nlambda, seed = seed
      )
    )
    class(result) <- "wtna_perm_mixed"
    return(result)
  }

  # ---- mcml dispatch: convert to netobject_group via as_tna ----
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(y, "mcml")) y <- as_tna(y)

  # ---- Single netobject_group: all-pairs permutation tests ----
  if (inherits(x, "netobject_group") && is.null(y)) {
    grp_names <- names(x)
    n_grps <- length(grp_names)
    if (n_grps < 2L) {
      stop("Need at least 2 groups for pairwise permutation tests.",
           call. = FALSE)
    }
    pairs <- combn(n_grps, 2L)
    results <- lapply(seq_len(ncol(pairs)), function(k) {
      i <- pairs[1L, k]
      j <- pairs[2L, k]
      permutation(x[[i]], x[[j]], iter = iter, alpha = alpha,
                  paired = paired, adjust = adjust,
                  nlambda = nlambda, seed = seed)
    })
    pair_labels <- vapply(seq_len(ncol(pairs)), function(k) {
      paste(grp_names[pairs[1L, k]], "vs", grp_names[pairs[2L, k]])
    }, character(1))
    names(results) <- pair_labels
    class(results) <- c("net_permutation_group", "list")
    return(results)
  }

  # ---- netobject_group dispatch: permute each matching element ----
  if (inherits(x, "netobject_group") && inherits(y, "netobject_group")) {
    common <- intersect(names(x), names(y))
    if (length(common) == 0L) {
      stop("No matching group names between x and y.", call. = FALSE)
    }
    results <- lapply(common, function(nm) {
      permutation(x[[nm]], y[[nm]], iter = iter, alpha = alpha,
                  paired = paired, adjust = adjust, nlambda = nlambda,
                  seed = seed)
    })
    names(results) <- common
    class(results) <- c("net_permutation_group", "list")
    return(results)
  }

  # ---- Coerce cograph_network inputs ----
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(y, "cograph_network")) y <- .as_netobject(y)

  # ---- Input validation ----
  stopifnot(
    inherits(x, "netobject"),
    inherits(y, "netobject"),
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.logical(paired), length(paired) == 1,
    is.character(adjust), length(adjust) == 1
  )
  iter <- as.integer(iter)

  if (is.null(x$data)) {
    stop("'x' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }
  if (is.null(y$data)) {
    stop("'y' does not contain $data. Rebuild with build_network().",
         call. = FALSE)
  }

  if (x$method != y$method) {
    stop("Methods must match: x uses '", x$method,
         "', y uses '", y$method, "'.", call. = FALSE)
  }

  if (!setequal(x$nodes$label, y$nodes$label)) {
    stop("Nodes must be the same in both networks.", call. = FALSE)
  }

  # Ensure same node order
  nodes <- x$nodes$label
  if (!identical(x$nodes$label, y$nodes$label)) {
    y$weights <- y$weights[nodes, nodes]
  }

  method <- .resolve_method_alias(x$method)
  directed <- x$directed
  n_nodes <- length(nodes)

  if (paired) {
    if (nrow(x$data) != nrow(y$data)) {
      stop("Paired test requires equal number of observations in x and y.",
           call. = FALSE)
    }
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Observed difference ----
  obs_diff <- x$weights - y$weights

  # ---- Dispatch permutation ----
  has_data_x <- is.data.frame(x$data) && ncol(x$data) > 0L
  has_data_y <- is.data.frame(y$data) && ncol(y$data) > 0L
  if (method %in% c("relative", "frequency", "co_occurrence")) {
    if (!has_data_x || !has_data_y) {
      stop("Permutation test requires the original data stored in the netobject. ",
           "For wtna/cna networks, use wtna() directly instead of ",
           "build_network(method='cna').", call. = FALSE)
    }
    x$data <- .resampling_transition_data(x$data, x$metadata, x$params)
    y$data <- .resampling_transition_data(y$data, y$metadata, y$params)
    n_seq_x <- .transition_resampling_n_sequences(x$data, x$params)
    n_seq_y <- .transition_resampling_n_sequences(y$data, y$params)
    if ((!is.na(n_seq_x) && n_seq_x <= 1L) ||
        (!is.na(n_seq_y) && n_seq_y <= 1L)) {
      warning(
        "A network with one long sequence is not recommended and can't be ",
        "validated using bootstrap and other confirmatory testings.",
        call. = FALSE
      )
    }
    perm_result <- .permutation_transition(
      x = x, y = y, nodes = nodes, method = method,
      iter = iter, paired = paired
    )
  } else {
    perm_result <- .permutation_association(
      x = x, y = y, nodes = nodes, method = method,
      iter = iter, paired = paired, nlambda = nlambda
    )
  }

  # ---- P-values ----
  # (sum(|perm_diff| >= |obs_diff|) + 1) / (iter + 1)
  obs_flat <- as.vector(obs_diff)
  p_values_flat <- (perm_result$exceed_counts + 1L) / (iter + 1L)

  # Apply multiple comparison correction
  p_values_flat <- p.adjust(p_values_flat, method = adjust)

  p_mat <- matrix(p_values_flat, n_nodes, n_nodes,
                  dimnames = list(nodes, nodes))

  # ---- Effect size ----
  # Cohen's d style: observed_diff / sd(perm_diffs)
  perm_sd <- perm_result$perm_sd
  perm_sd[perm_sd == 0] <- NA_real_
  es_flat <- obs_flat / perm_sd
  es_flat[is.na(es_flat)] <- 0
  es_mat <- matrix(es_flat, n_nodes, n_nodes,
                   dimnames = list(nodes, nodes))

  # ---- Significant diff ----
  sig_mask <- (p_mat < alpha) * 1
  diff_sig <- obs_diff * sig_mask

  # ---- Summary ----
  summary_df <- .build_permutation_summary(
    obs_diff = obs_diff,
    p_mat = p_mat,
    es_mat = es_mat,
    x_matrix = x$weights,
    y_matrix = y$weights,
    nodes = nodes,
    directed = directed,
    alpha = alpha
  )

  # ---- Assemble result ----
  result <- list(
    x           = x,
    y           = y,
    diff        = obs_diff,
    diff_sig    = diff_sig,
    p_values    = p_mat,
    effect_size = es_mat,
    summary     = summary_df,
    method      = method,
    iter        = iter,
    alpha       = alpha,
    paired      = paired,
    adjust      = adjust
  )
  class(result) <- "net_permutation"
  result
}


# ---- Transition fast path ----

#' Permutation test for transition networks via pre-computed counts
#' @noRd
.permutation_transition <- function(x, y, nodes, method, iter, paired) {
  n_nodes <- length(nodes)
  nbins <- n_nodes * n_nodes
  is_relative <- method == "relative"

  # Pre-compute per-sequence counts for both groups
  trans_x <- .precompute_per_sequence(x$data, method, x$params, nodes)
  trans_y <- .precompute_per_sequence(y$data, method, y$params, nodes)

  n_x <- nrow(trans_x)
  n_y <- nrow(trans_y)

  # Pool sequences
  pooled <- rbind(trans_x, trans_y)
  n_total <- n_x + n_y

  # Observed diff (recomputed from counts for consistency)
  obs_flat <- as.vector(x$weights - y$weights)

  # Running counters
  exceed_counts <- integer(nbins)
  sum_diffs <- numeric(nbins)
  sum_diffs_sq <- numeric(nbins)

  for (i in seq_len(iter)) {
    if (paired) {
      # Paired: randomly swap x/y within each pair
      swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
      idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
      idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
      counts_x <- colSums(pooled[idx_x, , drop = FALSE])
      counts_y <- colSums(pooled[idx_y, , drop = FALSE])
    } else {
      # Unpaired: shuffle group labels
      idx_x <- sample.int(n_total, n_x)
      counts_x <- colSums(pooled[idx_x, , drop = FALSE])
      counts_y <- colSums(pooled[-idx_x, , drop = FALSE])
    }

    # Post-process each group to get network matrix
    mat_x <- .postprocess_counts(counts_x, n_nodes, is_relative,
                                 x$scaling, x$threshold)
    mat_y <- .postprocess_counts(counts_y, n_nodes, is_relative,
                                 y$scaling, y$threshold)

    perm_diff <- as.vector(mat_x) - as.vector(mat_y)

    # Accumulate
    exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
    sum_diffs <- sum_diffs + perm_diff
    sum_diffs_sq <- sum_diffs_sq + perm_diff^2
  }

  # SD of permutation diffs
  perm_mean <- sum_diffs / iter
  perm_sd <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))

  list(
    exceed_counts = exceed_counts,
    perm_sd = perm_sd
  )
}


#' Convert flat count vector to network matrix with post-processing
#' @noRd
.postprocess_counts <- function(counts, n_nodes, is_relative,
                                scaling, threshold) {
  mat <- matrix(counts, n_nodes, n_nodes, byrow = TRUE)
  if (is_relative) {
    rs <- rowSums(mat)
    nz <- rs > 0
    mat[nz, ] <- mat[nz, ] / rs[nz]
  }
  if (!is.null(scaling)) mat <- .apply_scaling(mat, scaling) # nocov start
  if (threshold > 0) mat[abs(mat) < threshold] <- 0 # nocov end
  mat
}


# ---- Association path (optimized) ----

#' Permutation test for association networks
#'
#' Pre-cleans pooled data once, then uses lightweight per-iteration
#' estimation (direct cor/solve/glasso calls) to avoid repeated
#' input validation overhead. For custom/unknown estimators, falls
#' back to full estimator calls.
#' @noRd
.permutation_association <- function(x, y, nodes, method, iter, paired,
                                    nlambda = 50L) {
  n_nodes <- length(nodes)
  nbins <- n_nodes * n_nodes

  # $data is already cleaned by the estimator (numeric matrix, no NAs,
  # no zero-variance columns) — just pool directly
  n_x <- nrow(x$data)
  n_y <- nrow(y$data)
  pooled_mat <- rbind(x$data, y$data)
  n_total <- n_x + n_y

  # Extract params
  params_x <- x$params
  cor_method <- .param_get(params_x, "cor_method", "pearson")
  threshold_x <- x$threshold
  threshold_y <- y$threshold
  scaling_x <- x$scaling
  scaling_y <- y$scaling
  obs_flat <- as.vector(x$weights - y$weights)

  # Select fast path based on method
  use_fast <- method %in% c("cor", "pcor", "glasso")

  # Pre-compute glasso lambda path: narrowed around original lambdas,
  # solved via glassopath (single Fortran call for entire path)
  if (use_fast && method == "glasso") {
    gamma <- .param_get(params_x, "gamma", 0.5)
    penalize_diag <- .param_get(params_x, "penalize.diagonal", FALSE)

    # Compute lambda path from pooled correlation, using glassopath
    # for the per-iteration solve (single Fortran call for full path)
    S_pooled <- cor(pooled_mat, method = cor_method)
    perm_rholist <- .compute_lambda_path(S_pooled, nlambda, 0.01)
    p_glasso <- ncol(pooled_mat)
  }

  # Build the per-iteration estimator function
  if (use_fast) {
    estimate_from_rows <- switch(method,
      cor = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        diag(S) <- 0
        S
      },
      pcor = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        Wi <- tryCatch(solve(S), error = function(e) NULL)
        if (is.null(Wi)) return(NULL) # nocov
        .precision_to_pcor(Wi, threshold = 0)
      },
      glasso = function(mat_subset) {
        S <- cor(mat_subset, method = cor_method)
        n_obs <- nrow(mat_subset)
        gp <- tryCatch(
          glasso::glassopath(s = S, rholist = perm_rholist, trace = 0,
                             penalize.diagonal = penalize_diag),
          error = function(e) NULL
        )
        if (is.null(gp)) return(NULL) # nocov
        # EBIC selection across the path
        best_wi <- .select_ebic_from_path(
          gp, S, n_obs, gamma, p_glasso, perm_rholist
        )
        if (is.null(best_wi)) return(NULL) # nocov
        .precision_to_pcor(best_wi, threshold = 0)
      }
    )
  } else {
    # Fallback: full estimator for custom methods
    estimator <- get_estimator(method)
    estimate_from_rows <- function(mat_subset) {
      df <- as.data.frame(mat_subset)
      est <- tryCatch(
        do.call(estimator$fn, c(list(data = df), params_x)),
        error = function(e) NULL
      )
      if (is.null(est)) return(NULL)
      mat <- est$matrix
      if (!identical(rownames(mat), nodes)) {
        common <- intersect(nodes, rownames(mat))
        if (length(common) < n_nodes) return(NULL)
        mat <- mat[nodes, nodes] # nocov
      }
      mat
    }
  }

  # Running counters
  exceed_counts <- integer(nbins)
  sum_diffs <- numeric(nbins)
  sum_diffs_sq <- numeric(nbins)

  for (i in seq_len(iter)) {
    if (paired) {
      swaps <- sample(c(TRUE, FALSE), n_x, replace = TRUE)
      idx_x <- ifelse(swaps, seq(n_x + 1L, n_total), seq_len(n_x))
      idx_y <- ifelse(swaps, seq_len(n_x), seq(n_x + 1L, n_total))
    } else {
      idx_x <- sample.int(n_total, n_x)
      idx_y <- seq_len(n_total)[-idx_x]
    }

    mat_x <- estimate_from_rows(pooled_mat[idx_x, , drop = FALSE])
    mat_y <- estimate_from_rows(pooled_mat[idx_y, , drop = FALSE])

    if (is.null(mat_x) || is.null(mat_y)) next

    # Apply scaling and threshold
    if (!is.null(scaling_x)) mat_x <- .apply_scaling(mat_x, scaling_x) # nocov
    if (threshold_x > 0) mat_x[abs(mat_x) < threshold_x] <- 0
    if (!is.null(scaling_y)) mat_y <- .apply_scaling(mat_y, scaling_y) # nocov
    if (threshold_y > 0) mat_y[abs(mat_y) < threshold_y] <- 0

    perm_diff <- as.vector(mat_x) - as.vector(mat_y)

    exceed_counts <- exceed_counts + (abs(perm_diff) >= abs(obs_flat))
    sum_diffs <- sum_diffs + perm_diff
    sum_diffs_sq <- sum_diffs_sq + perm_diff^2
  }

  perm_mean <- sum_diffs / iter
  perm_sd <- sqrt(pmax(sum_diffs_sq / iter - perm_mean^2, 0))

  list(
    exceed_counts = exceed_counts,
    perm_sd = perm_sd
  )
}


#' Select best precision matrix from glassopath output via EBIC
#'
#' Vectorized EBIC computation over the 3D wi array from glassopath.
#' @noRd
.select_ebic_from_path <- function(gp, S, n, gamma, p, rholist) {
  n_lambda <- length(rholist)
  best_ebic <- Inf
  best_wi <- NULL

  for (k in seq_len(n_lambda)) {
    wi_k <- gp$wi[, , k]

    log_det <- determinant(wi_k, logarithm = TRUE)
    if (log_det$sign <= 0) next # nocov
    log_det_val <- as.numeric(log_det$modulus)

    loglik <- (n / 2) * (log_det_val - sum(diag(S %*% wi_k)))
    npar <- sum(abs(wi_k[upper.tri(wi_k)]) > 1e-10)
    ebic <- -2 * loglik + npar * log(n) + 4 * npar * gamma * log(p)

    if (ebic < best_ebic) {
      best_ebic <- ebic
      best_wi <- wi_k
    }
  }

  best_wi
}


# ---- Summary builder ----

#' Build long-format summary data frame from permutation test results
#' @noRd
.build_permutation_summary <- function(obs_diff, p_mat, es_mat,
                                       x_matrix, y_matrix,
                                       nodes, directed, alpha) {
  n <- length(nodes)
  dt <- data.table::data.table(
    from        = rep(nodes, each = n),
    to          = rep(nodes, times = n),
    weight_x    = as.vector(t(x_matrix)),
    weight_y    = as.vector(t(y_matrix)),
    diff        = as.vector(t(obs_diff)),
    effect_size = as.vector(t(es_mat)),
    p_value     = as.vector(t(p_mat)),
    sig         = as.vector(t(p_mat)) < alpha
  )

  # Filter: keep edges present in either network
  if (directed) {
    dt <- dt[weight_x != 0 | weight_y != 0]
  } else {
    dt <- dt[(weight_x != 0 | weight_y != 0) & from <= to]
  }

  as.data.frame(dt)
}


# ---- S3 Methods ----

#' Print Method for net_permutation
#'
#' @param x A \code{net_permutation} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' print(perm)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
#'                  V3 = c("C","A","C"))
#' d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
#'                  V3 = c("B","C","B"))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 20, seed = 1)
#' print(perm)
#' }
#'
#' @export
print.net_permutation <- function(x, ...) {
  method_labels <- c(
    relative      = "Transition Network (relative probabilities)",
    frequency     = "Transition Network (frequency counts)",
    co_occurrence = "Co-occurrence Network",
    glasso        = "Partial Correlation Network (EBICglasso)",
    pcor          = "Partial Correlation Network (unregularised)",
    cor           = "Correlation Network",
    attention     = "Attention Network (decay-weighted transitions)",
    wtna          = "Window TNA (transitions)"
  )
  label <- if (x$method %in% names(method_labels)) {
    method_labels[[x$method]]
  } else {
    sprintf("Network (method: %s)", x$method)
  }

  dir_label <- if (x$x$directed) " [directed]" else " [undirected]"

  cat("Permutation Test:", label, dir_label, "\n", sep = "")
  cat(sprintf("  Iterations: %d  |  Alpha: %.2f",
              x$iter, x$alpha))
  if (x$paired) cat("  |  Paired")
  if (x$adjust != "none") cat(sprintf("  |  Adjust: %s", x$adjust))
  cat("\n")

  n_sig <- sum(x$summary$sig)
  n_total <- nrow(x$summary)
  cat(sprintf("  Nodes: %d  |  Edges tested: %d  |  Significant: %d\n",
              x$x$n_nodes, n_total, n_sig))

  invisible(x)
}


#' Summary Method for net_permutation
#'
#' @param object A \code{net_permutation} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with edge-level permutation test results.
#'
#' @examples
#' s1 <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' s2 <- data.frame(V1 = c("A","C","B"), V2 = c("C","B","A"))
#' n1 <- build_network(s1, method = "relative")
#' n2 <- build_network(s2, method = "relative")
#' perm <- permutation(n1, n2, iter = 10)
#' summary(perm)
#' \donttest{
#' set.seed(1)
#' d1 <- data.frame(V1 = c("A","B","A"), V2 = c("B","C","B"),
#'                  V3 = c("C","A","C"))
#' d2 <- data.frame(V1 = c("C","A","C"), V2 = c("A","B","A"),
#'                  V3 = c("B","C","B"))
#' net1 <- build_network(d1, method = "relative")
#' net2 <- build_network(d2, method = "relative")
#' perm <- permutation(net1, net2, iter = 20, seed = 1)
#' summary(perm)
#' }
#'
#' @export
summary.net_permutation <- function(object, ...) {
  object$summary
}


#' Print Method for net_permutation_group
#'
#' @param x A \code{net_permutation_group} object.
#' @param ... Additional arguments (ignored).
#' @return \code{x} invisibly.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'   V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 10)
#' print(perm)
#' \donttest{
#' set.seed(1)
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'                  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
#' print(perm)
#' }
#' @export
print.net_permutation_group <- function(x, ...) {
  cat("Grouped Permutation Test\n")
  cat("Groups:", paste(names(x), collapse = ", "), "\n")
  cat("Use summary() on each element for edge-level results.\n")
  invisible(x)
}

#' Summary Method for net_permutation_group
#'
#' Returns a combined summary data frame across all groups.
#'
#' @param object A \code{net_permutation_group} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame with group, edge, p_value, and sig columns.
#' @examples
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'   V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'   V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 10)
#' summary(perm)
#' \donttest{
#' set.seed(1)
#' s1 <- data.frame(V1 = c("A","B","A","C"), V2 = c("B","C","B","A"),
#'                  V3 = c("C","A","C","B"), grp = c("X","X","Y","Y"))
#' s2 <- data.frame(V1 = c("C","A","C","B"), V2 = c("A","B","A","C"),
#'                  V3 = c("B","C","B","A"), grp = c("X","X","Y","Y"))
#' nets1 <- build_network(s1, method = "relative", group = "grp")
#' nets2 <- build_network(s2, method = "relative", group = "grp")
#' perm  <- permutation(nets1, nets2, iter = 20, seed = 1)
#' summary(perm)
#' }
#' @export
summary.net_permutation_group <- function(object, ...) {
  do.call(rbind, lapply(names(object), function(nm) {
    df      <- object[[nm]]$summary
    df$group <- nm
    df[c("group", setdiff(names(df), "group"))]
  }))
}


#' Print Method for wtna_perm_mixed
#'
#' @param x A \code{wtna_perm_mixed} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.wtna_perm_mixed <- function(x, ...) {
  cat("Mixed WTNA Permutation Test (transition + co-occurrence)\n")
  cat("-- Transition (directed) --\n")
  print(x$transition)
  cat("-- Co-occurrence (undirected) --\n")
  print(x$cooccurrence)
  invisible(x)
}


#' Summary Method for wtna_perm_mixed
#'
#' @param object A \code{wtna_perm_mixed} object.
#' @param ... Additional arguments (ignored).
#' @return A list with transition and co-occurrence permutation summaries.
#' @export
summary.wtna_perm_mixed <- function(object, ...) {
  list(
    transition = summary(object$transition),
    cooccurrence = summary(object$cooccurrence)
  )
}
