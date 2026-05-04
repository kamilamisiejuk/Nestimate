# Mixed Markov Model — latent-class mixture of Markov chains
#
# Discovers subgroups with different transition dynamics via EM.
# Uses precomputed per-sequence transition counts for speed.

# ---------------------------------------------------------------------------
# Core EM
# ---------------------------------------------------------------------------

#' Run one EM from a given initialization
#' @param counts N x K^2 per-sequence transition counts.
#'   Column j = (from-1)*K + to encodes the (from, to) transition pair.
#' @param init_state N-length integer vector of initial states (1..K)
#' @param n_comp Number of mixture components
#' @param max_iter Max iterations
#' @param tol Log-likelihood convergence tolerance
#' @param smooth Laplace smoothing
#' @param K Number of states
#' @param from_ind Pre-computed K^2 x K indicator grouping by "from" state
#' @return List with P_all, init_all, pi_mix, posterior, ll, iterations, converged
#' @noRd
.mmm_em <- function(counts, init_state, n_comp, max_iter, tol, smooth, K,
                    init_posterior = NULL, from_ind = NULL, cov_df = NULL) {
  N <- nrow(counts)
  K2 <- K * K

  # Initialize posterior (N x M)
  if (is.null(init_posterior)) {
    post <- matrix(runif(N * n_comp), nrow = N, ncol = n_comp)
    post <- post / .rowSums(post, N, n_comp)
  } else {
    post <- init_posterior
  }

  # Index of "from" state for each K^2 pair (K^2-length integer vector).
  # pair j = (from-1)*K + to, so from = (j-1) %/% K + 1.
  from_idx <- rep(seq_len(K), each = K)

  # Initial state indicator (N x K, sparse: one 1 per row).
  # Used in the M-step only — the E-step indexes log(init_all) directly.
  init_ind <- matrix(0, N, K)
  valid_init <- !is.na(init_state)
  init_ind[cbind(which(valid_init), init_state[valid_init])] <- 1
  init_state_safe <- init_state
  init_state_safe[!valid_init] <- 1L

  ll_prev <- -Inf
  converged <- FALSE
  cov_beta <- NULL

  # Pre-compute design matrix for covariates (once)
  X_cov <- NULL
  if (!is.null(cov_df)) {
    X_cov <- stats::model.matrix(~ ., data = cov_df)
  }

  log_lik <- matrix(0, N, n_comp)

  for (iter in seq_len(max_iter)) {
    # ---- M-step ----
    # Compute transition table directly in K^2 x M shape so the E-step
    # can use it without a transpose. Pure-R BLAS one-shot:
    #   counts is N x K^2, post is N x M, so crossprod(counts, post)
    #   = t(counts) %*% post is K^2 x M.
    P_unnorm <- crossprod(counts, post) + smooth        # K^2 x M
    # Sum over from-state groups via 3D array view. P_unnorm has columns
    # indexed by j = (f-1)*K + t (column-major), so reshaping to
    # (K_to, K_from, M) and colSums-ing over the to-axis gives K_from x M.
    dim(P_unnorm) <- c(K, K, n_comp)
    from_sums_t <- colSums(P_unnorm)                    # K x M
    dim(P_unnorm) <- c(K2, n_comp)
    from_sums_t[from_sums_t == 0] <- 1
    # Expand K x M to K^2 x M by repeating each row K times in place.
    P_all <- P_unnorm / from_sums_t[from_idx, , drop = FALSE]   # K^2 x M

    # Initial state probabilities: K x M directly (skip the t() at the end).
    init_unnorm <- crossprod(init_ind, post) + smooth     # K x M
    init_sums <- .colSums(init_unnorm, K, n_comp)         # M
    init_sums[init_sums == 0] <- 1
    init_all <- init_unnorm / rep(init_sums, each = K)    # K x M

    # Mixing proportions
    if (is.null(cov_df)) {
      pi_mix <- .colMeans(post, N, n_comp)
      log_pi <- log(pi_mix + 1e-300)
    } else {
      sm <- .mmm_softmax_mstep_fast(post, X_cov, beta_prev = cov_beta,
                                      n_steps = 3L)
      log_pi_mat <- sm$log_pi_mat
      cov_beta <- sm$beta
      pi_mix <- .colMeans(exp(log_pi_mat), N, n_comp)
    }

    # ---- E-step ----
    # log_P: K^2 x M. log_init: K x M. log_pi: M (or N x M with covariates).
    log_P <- log(P_all + 1e-300)
    log_init <- log(init_all + 1e-300)

    # Transition log-likelihood: N x M = (N x K^2) %*% (K^2 x M)
    log_lik <- counts %*% log_P
    # Initial-state log-likelihood: each row i needs log_init[init_state[i], ].
    # init_ind %*% log_init reduces to indexed extraction.
    log_lik <- log_lik + log_init[init_state_safe, , drop = FALSE]
    # Mixing
    if (is.null(cov_df)) {
      log_lik <- log_lik + rep(log_pi, each = N)
    } else {
      log_lik <- log_lik + log_pi_mat
    }

    # Log-sum-exp row max: pmax.int is variadic, so a single C call covers
    # any n_comp. Special-cased only because column extraction is unrolled.
    if (n_comp == 2L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L])
    } else if (n_comp == 3L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L], log_lik[, 3L])
    } else if (n_comp == 4L) {
      log_max <- pmax.int(log_lik[, 1L], log_lik[, 2L], log_lik[, 3L], log_lik[, 4L])
    } else {
      log_max <- do.call(pmax.int, lapply(seq_len(n_comp), function(m) log_lik[, m]))
    }
    exp_lik <- exp(log_lik - log_max)
    row_sums <- .rowSums(exp_lik, N, n_comp)
    row_sums[row_sums == 0] <- 1e-300
    post <- exp_lik / row_sums

    ll <- sum(log_max + log(row_sums))

    if (abs(ll - ll_prev) < tol) {
      converged <- TRUE
      break
    }
    # For covariate models: also check relative LL change (handles
    # perfect separation where LL grows toward 0 but never converges)
    if (!is.null(cov_df) && abs(ll) > 0 &&
        abs((ll - ll_prev) / ll) < tol) {
      converged <- TRUE
      break
    }
    ll_prev <- ll
  }

  list(
    P_all = P_all,
    init_all = init_all,
    pi_mix = pi_mix,
    posterior = post,
    ll = ll,
    iterations = iter,
    converged = converged,
    cov_beta = cov_beta
  )
}

# ---------------------------------------------------------------------------
# Initialization via k-means on count vectors
# ---------------------------------------------------------------------------

#' @noRd
.mmm_init_kmeans <- function(counts, n_comp) {
  N <- nrow(counts)
  if (N <= n_comp) {
    post <- matrix(0, nrow = N, ncol = n_comp)
    for (i in seq_len(N)) post[i, ((i - 1L) %% n_comp) + 1L] <- 1
    return(post)
  }

  km <- tryCatch(
    stats::kmeans(counts, centers = n_comp, nstart = 1L, iter.max = 20L),
    error = function(e) NULL
  )

  if (is.null(km)) {
    post <- matrix(runif(N * n_comp), nrow = N, ncol = n_comp)
    return(post / rowSums(post))
  }

  # Convert cluster assignments to soft posterior (vectorized)
  post <- matrix(0.01 / n_comp, nrow = N, ncol = n_comp)
  post[cbind(seq_len(N), km$cluster)] <- 0.99
  post / rowSums(post)
}

# ---------------------------------------------------------------------------
# Cluster quality metrics
# ---------------------------------------------------------------------------

#' Compute cluster quality metrics from posterior matrix
#' @noRd
.mmm_quality <- function(posterior, assignments, M) {
  N <- nrow(posterior)

  # Average Posterior Probability per class
  avepp <- vapply(seq_len(M), function(m) {
    idx <- which(assignments == m)
    if (length(idx) == 0L) return(NA_real_)
    mean(posterior[idx, m])
  }, numeric(1))

  # Max posterior per sequence — direct pmax.int beats `do.call(pmax,
  # as.data.frame(posterior))` which copies the matrix into a list.
  max_post <- if (M == 2L) {
    pmax.int(posterior[, 1L], posterior[, 2L])
  } else if (M == 3L) {
    pmax.int(posterior[, 1L], posterior[, 2L], posterior[, 3L])
  } else if (M == 4L) {
    pmax.int(posterior[, 1L], posterior[, 2L], posterior[, 3L], posterior[, 4L])
  } else {
    do.call(pmax.int, lapply(seq_len(M), function(m) posterior[, m]))
  }

  # Overall AvePP (mean of max posteriors)
  avepp_overall <- mean(max_post)

  # Entropy (normalized: 0 = perfect, 1 = random)
  ent_raw <- -sum(posterior * log(posterior + 1e-300))
  ent_max <- N * log(M)
  entropy <- if (ent_max > 0) ent_raw / ent_max else 0
  relative_entropy <- 1 - entropy

  # Classification error (proportion with max posterior < 0.5)
  classification_error <- mean(max_post < 0.5)

  # ICL classification entropy (vectorized via cbind indexing)
  class_ent <- -sum(log(posterior[cbind(seq_len(N), assignments)] + 1e-300))

  list(
    avepp = avepp,
    avepp_overall = avepp_overall,
    entropy = entropy,
    relative_entropy = relative_entropy,
    classification_error = classification_error,
    class_entropy = class_ent
  )
}

# ---------------------------------------------------------------------------
# Covariate M-step helpers
# ---------------------------------------------------------------------------

#' Fast softmax regression M-step via Newton-Raphson
#'
#' Replaces nnet::multinom() inside the EM loop. Does a few Newton steps
#' on the pre-computed design matrix. Much faster than calling multinom()
#' at every EM iteration.
#'
#' @param post N x k posterior matrix
#' @param X N x p design matrix (with intercept column)
#' @param beta_prev (k-1) x p coefficient matrix from previous iteration (or NULL)
#' @param n_steps Number of Newton-Raphson steps (default: 5)
#' @return List with log_pi_mat (N x k) and beta ((k-1) x p)
#' @noRd
.mmm_softmax_mstep_fast <- function(post, X, beta_prev = NULL, n_steps = 5L) {
  N <- nrow(X)
  p <- ncol(X)
  k <- ncol(post)
  km1 <- k - 1L

  # Initialize beta
  if (is.null(beta_prev)) {
    beta <- matrix(0, km1, p)
  } else {
    beta <- beta_prev
  }

  for (step in seq_len(n_steps)) {
    # Compute softmax probabilities: N x k
    eta <- X %*% t(beta)  # N x (k-1)
    eta_full <- cbind(0, eta)  # N x k (reference class = 0)
    eta_max <- do.call(pmax, as.data.frame(eta_full))
    exp_eta <- exp(eta_full - eta_max)
    pi_mat <- exp_eta / rowSums(exp_eta)

    # Gradient for each non-reference class
    residual <- post[, -1L, drop = FALSE] - pi_mat[, -1L, drop = FALSE]  # N x (k-1)

    if (km1 == 1L) {
      # k=2: logistic regression Newton step
      grad <- crossprod(X, residual)  # p x 1
      W <- pi_mat[, 1L] * pi_mat[, 2L]  # N-length
      H <- crossprod(X * W, X) + diag(1e-6, p)
      H_inv <- tryCatch(solve(H), error = function(e) NULL)
      if (is.null(H_inv)) break # nocov
      beta <- beta + t(H_inv %*% grad)
    } else {
      # General k: block Newton step
      beta_vec <- as.numeric(t(beta))
      grad_vec <- as.numeric(t(crossprod(X, residual)))  # km1*p

      delta <- numeric(km1 * p)
      for (m in seq_len(km1)) {
        W <- pi_mat[, m + 1L] * (1 - pi_mat[, m + 1L])
        H <- crossprod(X * W, X) + diag(1e-6, p)
        H_inv <- tryCatch(solve(H), error = function(e) NULL)
        if (is.null(H_inv)) next # nocov
        idx <- ((m - 1L) * p + 1L):(m * p)
        delta[idx] <- H_inv %*% grad_vec[idx]
      }
      beta <- matrix(beta_vec + delta, km1, p, byrow = TRUE)
    }
  }

  log_pi_mat <- log(pi_mat + 1e-300)
  list(log_pi_mat = log_pi_mat, beta = beta)
}


# ---------------------------------------------------------------------------
# Data decoding helper
# ---------------------------------------------------------------------------

#' Convert integer-encoded matrix (from tna/cograph) to character labels
#' @noRd
.mmm_decode_int_data <- function(raw_data, states) {
  raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
  if (is.integer(raw_data[[1]]) || is.numeric(raw_data[[1]])) {
    raw_data[] <- lapply(raw_data, function(col) {
      idx <- as.integer(col)
      ifelse(is.na(idx) | idx < 1L | idx > length(states), NA_character_,
             states[idx])
    })
  }
  raw_data
}

# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' Fit a Mixed Markov Model
#'
#' Discovers latent subgroups with different transition dynamics using
#' Expectation-Maximization. Each mixture component has its own transition
#' matrix. Sequences are probabilistically assigned to components.
#'
#' @param data A data.frame (wide format), \code{netobject}, or
#'   \code{tna} model. For tna objects, extracts the stored data.
#' @param k Integer. Number of mixture components. Default: 2.
#' @param n_starts Integer. Number of random restarts. Default: 50.
#' @param max_iter Integer. Maximum EM iterations per start. Default: 200.
#' @param tol Numeric. Convergence tolerance. Default: 1e-6.
#' @param smooth Numeric. Laplace smoothing constant. Default: 0.01.
#' @param seed Integer or NULL. Random seed.
#' @param covariates Optional. Covariates integrated into the EM algorithm
#'   to model covariate-dependent mixing proportions. Accepts formula,
#'   character vector, string, or data.frame (same forms as
#'   \code{\link{build_clusters}}). Unlike the post-hoc analysis in
#'   \code{build_clusters()}, these covariates directly influence cluster
#'   membership during estimation. Requires the \pkg{nnet} package.
#'
#' @return An object of class \code{net_mmm} with components:
#'   \describe{
#'     \item{models}{List of \code{netobject}s, one per component.}
#'     \item{k}{Number of components.}
#'     \item{mixing}{Numeric vector of mixing proportions.}
#'     \item{posterior}{N x k matrix of posterior probabilities.}
#'     \item{assignments}{Integer vector of hard assignments (1..k).}
#'     \item{quality}{List: \code{avepp} (per-class), \code{avepp_overall},
#'       \code{entropy}, \code{relative_entropy},
#'       \code{classification_error}.}
#'     \item{log_likelihood, BIC, AIC, ICL}{Model fit statistics.}
#'     \item{states}{Character vector of state names.}
#'   }
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' mmm
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
#'   V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
#' )
#' mmm <- build_mmm(seqs, k = 2, seed = 42)
#' print(mmm)
#' summary(mmm)
#' }
#'
#' @seealso \code{\link{compare_mmm}}, \code{\link{build_network}}
#'
#' @export
build_mmm <- function(data,
                      k = 2L,
                      n_starts = 50L,
                      max_iter = 200L,
                      tol = 1e-6,
                      smooth = 0.01,
                      seed = NULL,
                      covariates = NULL) {

  k <- as.integer(k)
  n_starts <- as.integer(n_starts)
  max_iter <- as.integer(max_iter)
  stopifnot(
    "'k' must be >= 2" = k >= 2L,
    "'n_starts' must be >= 1" = n_starts >= 1L,
    "'max_iter' must be >= 1" = max_iter >= 1L,
    "'tol' must be > 0" = tol > 0,
    "'smooth' must be >= 0" = smooth >= 0
  )

  # ---- Extract data and states ----
  network_method <- NULL
  build_args     <- NULL
  if (inherits(data, "tna") || inherits(data, "ftna")) {
    raw_data <- data$data
    states <- data$labels
    if (is.null(raw_data)) {
      stop("tna model does not contain $data. ", # nocov start
           "Rebuild with tna::tna(data) where data is kept.", call. = FALSE) # nocov end
    }
    raw_data <- .mmm_decode_int_data(raw_data, states)
  } else if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    # Pure cograph_network (from cograph package): coerce to decode integer data
    data <- .as_netobject(data)
    raw_data <- data$data
    states <- data$nodes$label
  } else if (inherits(data, "netobject")) {
    network_method <- data$method
    build_args     <- data$build_args
    raw_data <- data$data
    states <- data$nodes$label
  } else {
    raw_data <- as.data.frame(data)
    state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
  }

  n_states <- length(states)
  stopifnot("Need at least 2 states" = n_states >= 2L)

  # ---- Resolve covariates ----
  cov_df <- NULL
  if (!is.null(covariates)) {
    cov_resolved <- .resolve_covariates(covariates, data, nrow(raw_data))
    cov_df <- cov_resolved$cov_df
  }

  # Strip covariate columns from raw_data before state extraction
  if (!is.null(cov_df) && is.data.frame(raw_data)) {
    cov_col_names <- names(cov_df)
    keep <- setdiff(names(raw_data), cov_col_names)
    if (length(keep) > 0L) raw_data <- raw_data[, keep, drop = FALSE]
    # Re-extract states after removing covariate columns
    state_cols_clean <- .select_state_cols(raw_data, id = NULL, cols = NULL)
    mat <- as.matrix(raw_data[, state_cols_clean, drop = FALSE])
    states <- sort(unique(as.vector(mat)))
    states <- states[!is.na(states)]
    n_states <- length(states)
    stopifnot("Need at least 2 states" = n_states >= 2L)
  }

  # ---- Precompute per-sequence transition counts ----
  counts <- .precompute_per_sequence_wide(
    raw_data, method = "relative", cols = NULL, id_col = NULL, states = states
  )
  N <- nrow(counts)
  K2 <- n_states * n_states

  # Extract initial states (first non-NA state per sequence)
  state_cols <- .select_state_cols(raw_data, id = NULL, cols = NULL)
  first_col <- as.character(raw_data[[state_cols[1L]]])
  init_state <- match(first_col, states)

  # Handle covariate NAs: subset counts, init_state, cov_df
  if (!is.null(cov_df)) {
    complete <- stats::complete.cases(cov_df)
    n_dropped <- sum(!complete)
    if (n_dropped > 0L) {
      warning(sprintf(
        "Dropped %d rows with NA covariates (%d remaining).",
        n_dropped, sum(complete)
      ), call. = FALSE)
      counts <- counts[complete, , drop = FALSE]
      init_state <- init_state[complete]
      cov_df <- cov_df[complete, , drop = FALSE]
      raw_data <- raw_data[complete, , drop = FALSE]
      N <- nrow(counts)
    }
  }

  stopifnot("Need more sequences than components" = N > k)

  # ---- Run EM with screen-then-refine strategy ----
  if (!is.null(seed)) set.seed(seed)

  # Pre-compute from-state indicator (shared across all EM runs)
  from_idx <- rep(seq_len(n_states), each = n_states)
  from_ind <- matrix(0, K2, n_states)
  from_ind[cbind(seq_len(K2), from_idx)] <- 1

  # Screen phase
  screen_iter <- min(25L, max_iter)
  n_refine <- max(3L, n_starts %/% 10L)

  # Collect all initializations (before forking for seed control)
  inits <- vector("list", n_starts)
  inits[[1L]] <- .mmm_init_kmeans(counts, k)
  for (s in seq_len(n_starts - 1L)) {
    inits[[s + 1L]] <- matrix(runif(N * k), nrow = N, ncol = k)
    inits[[s + 1L]] <- inits[[s + 1L]] / rowSums(inits[[s + 1L]])
  }

  # Screen: short EM runs (parallel on Unix)
  n_cores <- 1L
  if (.Platform$OS.type == "unix" && n_starts > 1L) {
    detected <- parallel::detectCores(logical = FALSE)
    n_cores <- min(if (is.finite(detected)) detected else 1L, n_starts)
    if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_", "FALSE")))) { # nocov start
      n_cores <- 1L
    } # nocov end
  }

  .run_screen <- function(init) {
    # Screen phase: no covariates (fast), just find good transition inits
    .mmm_em(counts, init_state, k, screen_iter, tol, smooth, n_states,
            init_posterior = init, from_ind = from_ind, cov_df = NULL)
  }

  if (n_cores > 1L) {
    screen_results <- parallel::mclapply(inits, .run_screen,
                                          mc.cores = n_cores)
  } else {
    screen_results <- lapply(inits, .run_screen)
  }

  # Handle failed parallel runs
  screen_lls <- vapply(screen_results, function(r) {
    if (is.null(r)) -Inf else r$ll # nocov
  }, numeric(1))
  screen_conv <- vapply(screen_results, function(r) {
    if (is.null(r)) FALSE else r$converged # nocov
  }, logical(1))

  # Refine: top candidates + any that converged early
  top_idx <- order(screen_lls, decreasing = TRUE)[seq_len(n_refine)]
  conv_idx <- which(screen_conv)
  refine_idx <- unique(c(top_idx, conv_idx))

  best <- list(ll = -Inf)
  for (idx in refine_idx) {
    if (is.null(screen_results[[idx]])) next
    # Always refine when covariates present (screen ran without them)
    if (screen_results[[idx]]$converged && is.null(cov_df)) {
      run <- screen_results[[idx]]
    } else {
      run <- .mmm_em(counts, init_state, k, max_iter, tol, smooth, n_states,
                     init_posterior = screen_results[[idx]]$posterior,
                     from_ind = from_ind, cov_df = cov_df)
    }
    if (run$ll > best$ll) best <- run
  }

  # Fallback: if all parallel screen runs returned NULL (e.g. macOS arm64 fork
  # failures), run one guaranteed sequential EM before touching best$P_all.
  if (is.null(best$P_all)) {
    run <- .mmm_em(counts, init_state, k, max_iter, tol, smooth, n_states,
                   init_posterior = inits[[1L]], from_ind = from_ind,
                   cov_df = cov_df)
    if (!is.null(run)) best <- run
  }
  if (is.null(best$P_all)) {
    stop("EM algorithm failed to produce a valid result in all starting configurations.",
         call. = FALSE)
  }

  # ---- Build netobjects for each component ----
  models <- lapply(seq_len(k), function(m) {
    P_vec <- best$P_all[, m]
    P_mat <- t(matrix(P_vec, nrow = n_states, ncol = n_states))
    dimnames(P_mat) <- list(states, states)

    edges <- .extract_edges_from_matrix(P_mat, directed = TRUE)

    nodes_df <- data.frame(
      id = seq_along(states), label = states, name = states,
      x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
    )

    # Initial state probabilities from EM M-step (states x components matrix)
    init <- setNames(best$init_all[, m], states)

    structure(list(
      data = raw_data,
      weights = P_mat,
      nodes = nodes_df,
      edges = edges,
      directed = TRUE,
      method = "relative",
      params = list(),
      scaling = NULL,
      threshold = 0,
      n_nodes = n_states,
      n_edges = nrow(edges),
      level = NULL,
      initial = init,
      meta = list(source = "nestimate", layout = NULL,
                  tna = list(method = "relative")),
      node_groups = NULL
    ), class = c("netobject", "cograph_network"))
  })

  names(models) <- paste0("Cluster ", seq_len(k))

  # ---- Assignments & quality ----
  # max.col is the vectorized row-wise argmax — avoids `apply` overhead.
  assignments <- max.col(best$posterior, ties.method = "first")
  quality <- .mmm_quality(best$posterior, assignments, k)

  # ---- Information criteria ----
  n_params <- k * n_states * (n_states - 1L) +
              k * (n_states - 1L) +
              if (is.null(cov_df)) (k - 1L) else (k - 1L) * (ncol(cov_df) + 1L)
  BIC_val <- -2 * best$ll + n_params * log(N)
  AIC_val <- -2 * best$ll + 2 * n_params
  ICL_val <- BIC_val + 2 * quality$class_entropy

  # ---- Covariate output ----
  cov_result <- NULL
  if (!is.null(cov_df) && !is.null(best$cov_beta)) {
    # Run final nnet::multinom once for SEs and proper inference
    cov_result <- .run_covariate_analysis(
      assignments, cov_df, paste(names(cov_df), collapse = " + "), k
    )
    # Store EM-estimated beta
    cov_result$beta <- best$cov_beta
  }

  structure(list(
    models = models,
    k = k,
    mixing = best$pi_mix,
    posterior = best$posterior,
    assignments = assignments,
    quality = quality,
    log_likelihood = best$ll,
    BIC = BIC_val,
    AIC = AIC_val,
    ICL = ICL_val,
    n_params = n_params,
    iterations = best$iterations,
    converged = best$converged,
    states = states,
    n_sequences = N,
    covariates = cov_result,
    network_method = network_method,
    build_args     = build_args
  ), class = "net_mmm")
}

# ---------------------------------------------------------------------------
# Compare multiple k values
# ---------------------------------------------------------------------------

#' Compare MMM fits across different k
#'
#' @param data Data frame, netobject, or tna model.
#' @param k Integer vector of component counts. Default: 2:5.
#' @param ... Arguments passed to \code{\link{build_mmm}}.
#'
#' @return A \code{mmm_compare} data frame with BIC, AIC, ICL, AvePP,
#'   entropy per k.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' comp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
#' comp
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 30, TRUE), V2 = sample(LETTERS[1:3], 30, TRUE),
#'   V3 = sample(LETTERS[1:3], 30, TRUE), V4 = sample(LETTERS[1:3], 30, TRUE)
#' )
#' comp <- compare_mmm(seqs, k = 2:3, seed = 42)
#' print(comp)
#' }
#'
#' @export
compare_mmm <- function(data, k = 2:5, ...) {
  results <- lapply(k, function(m) {
    fit <- build_mmm(data, k = m, ...)
    data.frame(
      k = m,
      log_likelihood = fit$log_likelihood,
      AIC = fit$AIC,
      BIC = fit$BIC,
      ICL = fit$ICL,
      AvePP = fit$quality$avepp_overall,
      Entropy = fit$quality$entropy,
      converged = fit$converged,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, results)
  class(result) <- c("mmm_compare", "data.frame")
  result
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_mmm
#'
#' Compact summary of a Mixed Markov Model fit. Header carries dimensions
#' and information criteria; cluster table carries N, mixing share, and
#' per-cluster average posterior probability (AvePP). Layout matches
#' \code{\link{print.net_clustering}} so distance- and model-based
#' clusterings can be compared at a glance.
#'
#' @param x A \code{net_mmm} object.
#' @param digits Integer. Decimal places for floating-point statistics.
#'   Default \code{3}. Non-breaking: \code{print(x)} keeps the same
#'   alignment as before.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' print(mmm)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
#' print(mmm)
#' }
#'
#' @export
print.net_mmm <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  n_total <- as.integer(x$n_sequences)

  cat("Mixed Markov Model\n")
  cat(sprintf("  Sequences: %d  |  Clusters: %d  |  States: %d\n",
              n_total, k, length(x$states)))
  cat(sprintf("  ICs: LL = %.*f  |  BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
              digits, x$log_likelihood, digits, x$BIC,
              digits, x$AIC, digits, x$ICL))
  if (!is.null(x$quality)) {
    cat(sprintf(
      "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
      digits, x$quality$avepp_overall, digits, x$quality$entropy,
      x$quality$classification_error * 100))
  }
  if (isFALSE(x$converged)) {
    cat(sprintf("  Status: did not converge in %d iterations\n",
                as.integer(x$iterations)))
  }

  cat("\n")
  sizes <- as.integer(tabulate(x$assignments, nbins = k))
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total),
    `Mix%`  = sprintf("%4.1f%%", as.numeric(x$mixing) * 100),
    AvePP   = sprintf(paste0("%.", digits, "f"),
                      as.numeric(x$quality$avepp))
  )
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("\n  Covariates: %s (integrated, %d predictors)\n",
                paste(cov_names, collapse = ", "), length(cov_names)))
  }

  invisible(x)
}

#' Summary Method for net_mmm
#'
#' @param object A \code{net_mmm} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' summary(mmm)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
#' summary(mmm)
#' }
#'
#' @export
summary.net_mmm <- function(object, ...) {
  print(object)
  cat("\n")
  for (m in seq_len(object$k)) {
    n_in <- sum(object$assignments == m)
    cat(sprintf("--- Cluster %d (%.1f%%, n=%d) ---\n",
                m, object$mixing[m] * 100, n_in))
    print(round(object$models[[m]]$weights, 3))
    cat("\n")
  }

  # Covariate analysis
  if (!is.null(object$covariates)) {
    cov <- object$covariates
    .print_covariate_profiles(
      cov, "Covariate Analysis (integrated into EM -- influences cluster membership)"
    )
  }

  k <- object$k
  mean_posterior <- vapply(seq_len(k), function(m) {
    idx <- which(object$assignments == m)
    if (length(idx) == 0L) return(NA_real_)
    mean(object$posterior[idx, m])
  }, numeric(1L))

  data.frame(
    component      = seq_len(k),
    prior          = as.numeric(object$mixing),
    n_assigned     = as.integer(tabulate(object$assignments, nbins = k)),
    mean_posterior = mean_posterior,
    avepp          = as.numeric(object$quality$avepp),
    stringsAsFactors = FALSE,
    row.names      = NULL
  )
}

#' Plot Method for net_mmm
#'
#' @param x A \code{net_mmm} object.
#' @param type Character. Plot type: \code{"posterior"} (default) or \code{"covariates"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' mmm <- build_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' plot(mmm, type = "posterior")
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' mmm <- build_mmm(seqs, k = 2, n_starts = 5, seed = 1)
#' plot(mmm, type = "posterior")
#' }
#'
#' @export
plot.net_mmm <- function(x, type = c("posterior", "covariates"), ...) {
  type <- match.arg(type)

  if (type == "covariates") {
    if (is.null(x$covariates)) {
      stop("No covariate analysis found. Run build_mmm() with covariates.",
           call. = FALSE)
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects (ref: Cluster %s)",
              x$covariates$fit$reference_cluster)
    ))
  }

  if (type == "posterior") {
    return(.plot_mmm_posterior(x))
  }
}

#' Plot posterior probability distribution
#' @noRd
.plot_mmm_posterior <- function(x) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required.", call. = FALSE) # nocov
  }

  max_post <- apply(x$posterior, 1L, max)
  df <- data.frame(
    cluster = factor(x$assignments),
    max_posterior = max_post
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = max_posterior, fill = cluster)) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
    ggplot2::labs(x = "Max Posterior Probability", y = "Count",
                  title = "Classification Certainty", fill = "Cluster") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}

#' Print Method for mmm_compare
#'
#' @param x An \code{mmm_compare} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' cmp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
#' print(cmp)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' cmp <- compare_mmm(seqs, k = 2:3, n_starts = 5, seed = 1)
#' print(cmp)
#' }
#'
#' @export
print.mmm_compare <- function(x, ...) {
  cat("MMM Model Comparison\n\n")
  best_bic <- which.min(x$BIC)
  best_icl <- which.min(x$ICL)
  x$best <- ""
  x$best[best_bic] <- "<-- BIC"
  if (best_icl != best_bic) x$best[best_icl] <- paste(x$best[best_icl], "<-- ICL")
  print.data.frame(x, row.names = FALSE, right = FALSE)
  invisible(x)
}

#' Summary Method for mmm_compare
#'
#' @param object An \code{mmm_compare} object (a data.frame subclass).
#' @param ... Additional arguments (ignored).
#' @return A tidy data frame with one row per \code{k}, plus a \code{best}
#'   character column flagging the minimum-BIC and minimum-ICL solutions.
#' @export
summary.mmm_compare <- function(object, ...) {
  best_bic <- which.min(object$BIC)
  best_icl <- which.min(object$ICL)
  best <- rep("", nrow(object))
  best[best_bic] <- "BIC"
  if (length(best_icl) && best_icl != best_bic) {
    best[best_icl] <- if (nzchar(best[best_icl]))
      paste(best[best_icl], "ICL", sep = "+") else "ICL"
  }
  out <- as.data.frame(object)
  out$best <- best
  row.names(out) <- NULL
  out
}

#' Plot Method for mmm_compare
#'
#' @param x An \code{mmm_compare} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' cmp <- compare_mmm(seqs, k = 2:3, n_starts = 1, max_iter = 10, seed = 1)
#' plot(cmp)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 30, TRUE),
#'   V2 = sample(c("A","B","C"), 30, TRUE),
#'   V3 = sample(c("A","B","C"), 30, TRUE)
#' )
#' cmp <- compare_mmm(seqs, k = 2:3, n_starts = 5, seed = 1)
#' plot(cmp)
#' }
#'
#' @export
plot.mmm_compare <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) { # nocov start
    stop("Package 'ggplot2' required.", call. = FALSE)
  } # nocov end

  df <- data.frame(
    k = rep(x$k, 3),
    value = c(x$BIC, x$AIC, x$ICL),
    criterion = rep(c("BIC", "AIC", "ICL"), each = nrow(x)),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = k, y = value, color = criterion)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_continuous(breaks = x$k) +
    ggplot2::labs(x = "k (components)", y = "Information Criterion",
                  title = "MMM Model Selection", color = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")

  print(p)
  invisible(p)
}

# ---------------------------------------------------------------------------
# cluster_mmm — wrapper returning netobject_group (parallel to cluster_network)
# ---------------------------------------------------------------------------

# Attach MMM clustering metadata (assignments, posterior, mixing, ICs, full
# sequence data) to a netobject_group built from a net_mmm. Used by both
# cluster_mmm() and build_network.net_mmm so the netobject_group invariant
# is the same regardless of how the user got there: every member has its
# own $weights/$nodes/$edges, and `attr(, "clustering")` carries N-row
# $data matching $assignments.
#' @noRd
.attach_mmm_clustering <- function(nets, mmm, full_data = NULL) {
  clustering_info <- mmm[setdiff(names(mmm), "models")]
  if (is.null(full_data) && length(nets) > 0L) {
    full_data <- nets[[1L]]$data
  }
  if (!is.null(full_data)) clustering_info$data <- full_data
  class(clustering_info) <- "net_mmm_clustering"
  attr(nets, "clustering") <- clustering_info
  attr(nets, "group_col")  <- "cluster"
  nets
}

#' Cluster sequences using Mixed Markov Models
#'
#' Fits a mixture of Markov chains to sequence data and returns a
#' \code{netobject_group} containing per-cluster transition networks.
#' This is the MMM equivalent of \code{\link{cluster_network}} (which uses
#' distance-based clustering); both functions share the
#' \code{cluster_by = ...} surface argument so the call shape stays
#' uniform across clustering families.
#'
#' For the full \code{net_mmm} object with posterior probabilities, model
#' fit statistics, and S3 methods, use \code{\link{build_mmm}} instead.
#'
#' @inheritParams build_mmm
#' @param cluster_by Character. Accepted only as \code{"mmm"} (the
#'   default). Present so \code{cluster_mmm()} and \code{cluster_network()}
#'   share the same call shape; any other value raises an error pointing
#'   at \code{\link{cluster_network}}.
#' @param ... Reserved for forward compatibility with the unified
#'   \code{cluster_*} surface. Currently unused.
#' @return A \code{netobject_group} (list of \code{netobject}s, one per
#'   cluster). MMM-specific information is stored in
#'   \code{attr(, "clustering")} (class \code{"net_mmm_clustering"}):
#'   \describe{
#'     \item{assignments}{Integer vector of cluster assignments.}
#'     \item{k}{Number of clusters.}
#'     \item{posterior}{N x k matrix of posterior probabilities.}
#'     \item{mixing}{Mixing proportions.}
#'     \item{quality}{List with AvePP, entropy, classification error.}
#'     \item{BIC, AIC, ICL}{Model fit statistics.}
#'     \item{data}{The full N-row sequence frame, matching
#'       \code{$assignments} -- so \code{\link{sequence_plot}} and
#'       \code{\link{distribution_plot}} can recover both.}
#'   }
#' @seealso \code{\link{build_mmm}} for the full MMM object,
#'   \code{\link{cluster_network}} for distance-based clustering
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' grp[[1]]$weights
#' attr(grp, "clustering")$assignments
#' \donttest{
#' # Visualise with sequence_plot
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 40, TRUE),
#'   V2 = sample(LETTERS[1:3], 40, TRUE),
#'   V3 = sample(LETTERS[1:3], 40, TRUE)
#' )
#' grp <- cluster_mmm(seqs, k = 2)
#' sequence_plot(grp, type = "index")
#' }
#' @export
cluster_mmm <- function(data, k = 2L, n_starts = 50L, max_iter = 200L,
                        tol = 1e-6, smooth = 0.01, seed = NULL,
                        covariates = NULL,
                        cluster_by = "mmm", ...) {
  # cluster_by exists for API parity with cluster_network() so a single
  # surface argument toggles the clustering family. Only "mmm" is valid
  # here; anything else is a programming error worth catching loudly.
  if (!identical(as.character(cluster_by), "mmm")) {
    stop("cluster_mmm() only supports cluster_by = \"mmm\". For other ",
         "clustering algorithms, use cluster_network(..., cluster_by = ...).",
         call. = FALSE)
  }
  # Quietly ignore further `...` so the cluster_network()-style call
  # `cluster_mmm(x, k, dissimilarity = "hamming")` doesn't fail noisily.
  # (build_mmm has no `...` itself; we don't want to leak strange args
  # into the EM core either.)

  mmm <- build_mmm(data = data, k = k, n_starts = n_starts,
                   max_iter = max_iter, tol = tol, smooth = smooth,
                   seed = seed, covariates = covariates)

  grp <- mmm$models
  if (is.null(names(grp))) names(grp) <- paste0("Cluster ", seq_along(grp))
  class(grp) <- "netobject_group"
  .attach_mmm_clustering(grp, mmm, full_data = grp[[1L]]$data)
}

#' Print Method for MMM Clustering Attribute
#'
#' Prints the clustering metadata that \code{\link{cluster_mmm}} attaches
#' to its \code{netobject_group} return value (\code{attr(grp, "clustering")}).
#' Layout mirrors \code{\link{print.net_clustering}}: a one-line dimension
#' header, a quality line with AvePP / entropy / classification error,
#' information criteria, and a per-cluster table.
#'
#' @param x A \code{net_mmm_clustering} object.
#' @param digits Integer. Decimal places for floating-point statistics.
#'   Default \code{3}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 30, TRUE),
#'                    V2 = sample(c("A","B","C"), 30, TRUE))
#' grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 10, seed = 1)
#' print(attr(grp, "clustering"))
#'
#' @export
print.net_mmm_clustering <- function(x, digits = 3L, ...) {
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  n_total <- as.integer(x$n_sequences)

  cat("MMM Clustering [k = ", k, "]\n", sep = "")
  cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n_total, k))
  if (!is.null(x$quality)) {
    cat(sprintf(
      "  Quality: AvePP = %.*f  |  Entropy = %.*f  |  Class.Err = %.1f%%\n",
      digits, x$quality$avepp_overall, digits, x$quality$entropy,
      x$quality$classification_error * 100))
  }
  if (!is.null(x$BIC)) {
    cat(sprintf("  ICs: BIC = %.*f  |  AIC = %.*f  |  ICL = %.*f\n",
                digits, x$BIC, digits, x$AIC, digits, x$ICL))
  }

  cat("\n")
  sizes <- as.integer(tabulate(x$assignments, nbins = k))
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total),
    `Mix%`  = sprintf("%4.1f%%", as.numeric(x$mixing) * 100),
    AvePP   = sprintf(paste0("%.", digits, "f"),
                      as.numeric(x$quality$avepp))
  )
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("\n  Covariates: %s (integrated, %d predictors)\n",
                paste(cov_names, collapse = ", "), length(cov_names)))
  }

  invisible(x)
}

#' Plot Method for MMM Clustering Attribute
#'
#' Plot routines for the MMM clustering metadata attached to a
#' \code{netobject_group} by \code{\link{cluster_mmm}} (or
#' \code{\link{cluster_network}} with \code{cluster_by = "mmm"}).
#' Mirrors the type-driven surface of
#' \code{\link{plot.net_clustering}} but covers only the metrics the EM
#' fit produces -- there is no distance matrix on an MMM clustering, so
#' \code{"silhouette"} / \code{"mds"} / \code{"heatmap"} aren't defined
#' here. The \code{netobject_group} dispatcher
#' (\code{\link{plot.netobject_group}}) raises a clear error if you ask
#' for one of those on an MMM result.
#'
#' @param x A \code{net_mmm_clustering} object.
#' @param type Character. One of \code{"posterior"} (default; histogram of
#'   max posterior probability per sequence, coloured by cluster),
#'   \code{"covariates"} or its alias \code{"predictors"} (covariate
#'   forest plot when \code{cluster_mmm()} was run with \code{covariates}).
#' @param ... Currently unused.
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(V1 = sample(c("A","B","C"), 40, TRUE),
#'                    V2 = sample(c("A","B","C"), 40, TRUE))
#' grp <- cluster_mmm(seqs, k = 2, n_starts = 1, max_iter = 20, seed = 1)
#' plot(attr(grp, "clustering"), type = "posterior")
#' }
#' @export
plot.net_mmm_clustering <- function(x, type = c("posterior", "covariates",
                                                 "predictors"), ...) {
  type <- match.arg(type)
  if (type == "predictors") type <- "covariates"

  if (type == "covariates") {
    if (is.null(x$covariates)) {
      stop("No covariate analysis on this clustering. Re-run ",
           "cluster_mmm() (or build_mmm()) with covariates = ...",
           call. = FALSE)
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects (ref: Cluster %s)",
              x$covariates$fit$reference_cluster)
    ))
  }

  # type == "posterior": .plot_mmm_posterior() reads only $posterior +
  # $assignments, both carried by net_mmm_clustering.
  .plot_mmm_posterior(x)
}
