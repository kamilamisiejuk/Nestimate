# ---- Bootstrap for Regularized Partial Correlation Networks ----

#' Bootstrap for Regularized Partial Correlation Networks
#'
#' @description
#' Fast, single-call bootstrap for EBICglasso partial correlation networks.
#' Combines nonparametric edge/centrality bootstrap, case-dropping stability
#' analysis, edge/centrality difference tests, predictability CIs, and
#' thresholded network into one function. Designed as a faster alternative
#' to bootnet with richer output.
#'
#' @param x A data frame, numeric matrix (observations x variables), or
#'   a \code{netobject} with \code{method = "glasso"}.
#' @param iter Integer. Number of nonparametric bootstrap iterations
#'   (default: 1000).
#' @param cs_iter Integer. Number of case-dropping iterations per drop
#'   proportion (default: 500).
#' @param cs_drop Numeric vector. Drop proportions for CS-coefficient
#'   computation (default: \code{seq(0.1, 0.9, by = 0.1)}).
#' @param alpha Numeric. Significance level for CIs (default: 0.05).
#' @param gamma Numeric. EBIC hyperparameter (default: 0.5).
#' @param nlambda Integer. Number of lambda values in the regularization
#'   path (default: 100).
#' @param centrality Character vector. Centrality measures to compute.
#'   Built-in: \code{"strength"}, \code{"expected_influence"},
#'   \code{"betweenness"}, \code{"closeness"}.
#'   Custom measures beyond these require \code{centrality_fn}.
#'   Default: \code{c("strength", "expected_influence", "betweenness", "closeness")}.
#' @param centrality_fn Optional function. A custom centrality function
#'   that takes a weight matrix and returns a named list of centrality
#'   vectors. When \code{NULL} (default), only \code{"strength"} and
#'   \code{"expected_influence"} are computed via \code{rowSums}/
#'   \code{colSums}. When provided, the function is called as
#'   \code{centrality_fn(mat)} and should return a named list (e.g.,
#'   \code{list(closeness = ..., betweenness = ...)}).
#' @param cor_method Character. Correlation method: \code{"pearson"}
#'   (default), \code{"spearman"}, or \code{"kendall"}.
#' @param ncores Integer. Number of parallel cores for mclapply
#'   (default: 1, sequential).
#' @param seed Integer or NULL. RNG seed for reproducibility.
#'
#' @return An object of class \code{"boot_glasso"} containing:
#' \describe{
#'   \item{original_pcor}{Original partial correlation matrix.}
#'   \item{original_precision}{Original precision matrix.}
#'   \item{original_centrality}{Named list of original centrality vectors.}
#'   \item{original_predictability}{Named numeric vector of node R-squared.}
#'   \item{edge_ci}{Data frame of edge CIs (edge, weight, ci_lower, ci_upper,
#'     inclusion).}
#'   \item{edge_inclusion}{Named numeric vector of edge inclusion probabilities.}
#'   \item{thresholded_pcor}{Partial correlation matrix with non-significant
#'     edges zeroed.}
#'   \item{centrality_ci}{Named list of data frames (node, value, ci_lower,
#'     ci_upper) per centrality measure.}
#'   \item{cs_coefficient}{Named numeric vector of CS-coefficients per
#'     centrality measure.}
#'   \item{cs_data}{Data frame of case-dropping correlations (drop_prop,
#'     measure, correlation).}
#'   \item{edge_diff_p}{Symmetric matrix of pairwise edge difference p-values.}
#'   \item{centrality_diff_p}{Named list of symmetric p-value matrices per
#'     centrality measure.}
#'   \item{predictability_ci}{Data frame of node predictability CIs (node,
#'     r2, ci_lower, ci_upper).}
#'   \item{boot_edges}{iter x n_edges matrix of bootstrap edge weights.}
#'   \item{boot_centrality}{Named list of iter x p bootstrap centrality
#'     matrices.}
#'   \item{boot_predictability}{iter x p matrix of bootstrap R-squared.}
#'   \item{nodes}{Character vector of node names.}
#'   \item{n}{Sample size.}
#'   \item{p}{Number of variables.}
#'   \item{iter}{Number of nonparametric iterations.}
#'   \item{cs_iter}{Number of case-dropping iterations.}
#'   \item{cs_drop}{Drop proportions used.}
#'   \item{alpha}{Significance level.}
#'   \item{gamma}{EBIC hyperparameter.}
#'   \item{nlambda}{Lambda path length.}
#'   \item{centrality_measures}{Character vector of centrality measures.}
#'   \item{cor_method}{Correlation method.}
#'   \item{lambda_path}{Lambda sequence used.}
#'   \item{lambda_selected}{Selected lambda for original data.}
#'   \item{timing}{Named numeric vector with timing in seconds.}
#' }
#'
#' @examples
#' set.seed(1)
#' dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
#' net <- build_network(dat, method = "glasso")
#' bg <- boot_glasso(net, iter = 10, cs_iter = 5, centrality = "strength")
#' \donttest{
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' net <- build_network(as.data.frame(mat), method = "glasso")
#' boot <- boot_glasso(net, iter = 100, cs_iter = 50, seed = 42,
#'   centrality = c("strength", "expected_influence"))
#' print(boot)
#' summary(boot, type = "edges")
#' }
#'
#' @seealso \code{\link{build_network}}, \code{\link{bootstrap_network}}
#'
#' @importFrom stats cor quantile
#' @importFrom parallel mclapply
#' @export
boot_glasso <- function(x,
                        iter = 1000L,
                        cs_iter = 500L,
                        cs_drop = seq(0.1, 0.9, by = 0.1),
                        alpha = 0.05,
                        gamma = 0.5,
                        nlambda = 100L,
                        centrality = c("strength", "expected_influence",
                                       "betweenness", "closeness"),
                        centrality_fn = NULL,
                        cor_method = "pearson",
                        ncores = 1L,
                        seed = NULL) {

  # ---- Input validation ----
  stopifnot(
    is.numeric(iter), length(iter) == 1, iter >= 2,
    is.numeric(cs_iter), length(cs_iter) == 1, cs_iter >= 1,
    is.numeric(cs_drop), all(cs_drop > 0), all(cs_drop < 1),
    is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
    is.numeric(gamma), length(gamma) == 1, gamma >= 0,
    is.numeric(nlambda), length(nlambda) == 1, nlambda >= 2,
    is.character(centrality), length(centrality) >= 1,
    is.character(cor_method), length(cor_method) == 1,
    is.numeric(ncores), length(ncores) == 1, ncores >= 1
  )
  iter <- as.integer(iter)
  cs_iter <- as.integer(cs_iter)
  nlambda <- as.integer(nlambda)
  ncores <- as.integer(ncores)
  cor_method <- match.arg(cor_method, c("pearson", "spearman", "kendall"))
  centrality <- match.arg(centrality, c("strength", "expected_influence",
                                         "closeness", "betweenness"),
                          several.ok = TRUE)

  if (!is.null(centrality_fn)) {
    stopifnot("centrality_fn must be a function" = is.function(centrality_fn))
  }

  if (!is.null(seed)) {
    stopifnot(is.numeric(seed), length(seed) == 1)
    set.seed(seed)
  }

  # ---- Input dispatch ----
  if (inherits(x, "cograph_network")) x <- .as_netobject(x)
  if (inherits(x, "netobject")) {
    if (x$method != "glasso") {
      stop("netobject must use method = 'glasso'.", call. = FALSE)
    }
    if (is.null(x$data)) {
      stop("netobject does not contain $data. Rebuild with build_network().",
           call. = FALSE)
    }
    data_mat <- x$data
    if (!is.matrix(data_mat)) data_mat <- as.matrix(data_mat)
    # Extract params from netobject if available
    if (!is.null(x$params$gamma)) gamma <- x$params$gamma
    if (!is.null(x$params$nlambda)) nlambda <- as.integer(x$params$nlambda)
    if (!is.null(x$params$cor_method)) cor_method <- x$params$cor_method
  } else if (is.data.frame(x) || is.matrix(x)) {
    # Use .prepare_association_input for cleaning
    # Convert matrix to data frame so it goes through the raw-data path
    input <- if (is.matrix(x)) as.data.frame(x) else x
    prepared <- .prepare_association_input(input, cor_method = cor_method)
    data_mat <- prepared$mat
  } else {
    stop("x must be a data frame, numeric matrix, or netobject.",
         call. = FALSE)
  }

  n <- nrow(data_mat)
  p <- ncol(data_mat)
  nodes <- colnames(data_mat)

  if (p < 2) stop("At least 2 variables are required.", call. = FALSE) # nocov start
  if (n < 3) stop("At least 3 observations are required.", call. = FALSE) # nocov end

  t_start <- proc.time()["elapsed"]

  # ---- Phase 1: Original estimation ----
  S <- cor(data_mat, method = cor_method)
  lambda_path <- .compute_lambda_path(S, nlambda, 0.01)

  gp_orig <- tryCatch(
    glasso::glassopath(s = S, rholist = lambda_path, trace = 0,
                       penalize.diagonal = FALSE),
    error = function(e) stop("glassopath failed on original data: ",
                              e$message, call. = FALSE)
  )

  Wi_orig <- .select_ebic_from_path(gp_orig, S, n, gamma, p, lambda_path)
  if (is.null(Wi_orig)) {
    stop("EBIC selection failed on original data.", call. = FALSE) # nocov
  }

  lambda_selected <- .bg_find_selected_lambda(gp_orig, Wi_orig, lambda_path)

  pcor_orig <- .precision_to_pcor(Wi_orig, threshold = 0)
  pcor_orig <- (pcor_orig + t(pcor_orig)) / 2
  colnames(pcor_orig) <- rownames(pcor_orig) <- nodes

  # Original centrality
  cent_orig <- .bg_compute_centrality(pcor_orig, p, nodes, centrality,
                                       centrality_fn)


  # Original predictability
  pred_orig <- pmin(pmax(1 - 1 / diag(Wi_orig), 0), 1)
  names(pred_orig) <- nodes

  # Pre-allocate bootstrap storage
  n_upper <- p * (p - 1) / 2
  ut <- .bg_upper_tri_indices(p)

  t_phase1 <- proc.time()["elapsed"]

  # ---- Phase 2: Nonparametric bootstrap ----
  boot_fn <- function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    .bg_estimate_once(data_mat[idx, , drop = FALSE], p, gamma, nlambda,
                      penalize_diag = FALSE, cor_method, lambda_path,
                      centrality, centrality_fn)
  }

  if (ncores > 1L) {
    boot_results <- parallel::mclapply(seq_len(iter), boot_fn,
                                        mc.cores = ncores)
  } else {
    boot_results <- lapply(seq_len(iter), boot_fn)
  }

  # Unpack results into matrices
  boot_edges <- matrix(NA_real_, nrow = iter, ncol = n_upper)
  boot_centrality <- lapply(centrality, function(m) {
    matrix(NA_real_, nrow = iter, ncol = p,
           dimnames = list(NULL, nodes))
  })
  names(boot_centrality) <- centrality
  boot_pred <- matrix(NA_real_, nrow = iter, ncol = p,
                      dimnames = list(NULL, nodes))

  for (i in seq_len(iter)) {
    res <- boot_results[[i]]
    if (!is.null(res)) {
      boot_edges[i, ] <- res$edges
      for (m in centrality) {
        boot_centrality[[m]][i, ] <- res$centrality[[m]]
      }
      boot_pred[i, ] <- res$predictability
    }
  }

  edge_names <- .bg_build_edge_names(nodes, ut$row_idx, ut$col_idx)
  colnames(boot_edges) <- edge_names

  t_phase2 <- proc.time()["elapsed"]

  # ---- Phase 3: Case-dropping ----
  cs_result <- .bg_case_drop(
    data_mat, p, gamma, nlambda, cor_method, lambda_path, centrality,
    cs_iter, cs_drop, cent_orig, ncores, centrality_fn
  )

  cs_coef <- vapply(centrality, function(m) {
    .bg_cs_coefficient(cs_result$cors[[m]], cs_drop, 0.7, 0.95)
  }, numeric(1))
  names(cs_coef) <- centrality

  t_phase3 <- proc.time()["elapsed"]

  # ---- Phase 4: Statistics ----

  # Edge CIs
  orig_edge_vec <- pcor_orig[upper.tri(pcor_orig)]
  edge_ci_lower <- apply(boot_edges, 2, quantile, probs = alpha / 2,
                          na.rm = TRUE)
  edge_ci_upper <- apply(boot_edges, 2, quantile, probs = 1 - alpha / 2,
                          na.rm = TRUE)

  # Edge inclusion probability
  edge_inclusion <- colMeans(boot_edges != 0, na.rm = TRUE)
  names(edge_inclusion) <- edge_names

  # Thresholded network: zero out edges where CI includes zero
  thresholded <- .bg_threshold_network(pcor_orig, edge_ci_lower,
                                        edge_ci_upper, ut)

  # Edge CI data frame
  edge_ci_df <- data.frame(
    edge = edge_names,
    weight = orig_edge_vec,
    ci_lower = edge_ci_lower,
    ci_upper = edge_ci_upper,
    inclusion = edge_inclusion,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Centrality CIs
  centrality_ci <- lapply(centrality, function(m) {
    ci_lo <- apply(boot_centrality[[m]], 2, quantile, probs = alpha / 2,
                    na.rm = TRUE)
    ci_hi <- apply(boot_centrality[[m]], 2, quantile, probs = 1 - alpha / 2,
                    na.rm = TRUE)
    data.frame(
      node = nodes,
      value = cent_orig[[m]],
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  names(centrality_ci) <- centrality

  # Edge difference test
  edge_diff_p <- .bg_edge_diff_test(boot_edges)

  # Centrality difference test
  centrality_diff_p <- lapply(centrality, function(m) {
    .bg_centrality_diff_test(boot_centrality[[m]])
  })
  names(centrality_diff_p) <- centrality

  # Predictability CIs
  pred_ci_lower <- apply(boot_pred, 2, quantile, probs = alpha / 2,
                          na.rm = TRUE)
  pred_ci_upper <- apply(boot_pred, 2, quantile, probs = 1 - alpha / 2,
                          na.rm = TRUE)
  pred_ci_df <- data.frame(
    node = nodes,
    r2 = pred_orig,
    ci_lower = pred_ci_lower,
    ci_upper = pred_ci_upper,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # CS data frame
  cs_data <- cs_result$cs_data

  t_end <- proc.time()["elapsed"]

  timing <- c(
    total = unname(t_end - t_start),
    phase1 = unname(t_phase1 - t_start),
    bootstrap = unname(t_phase2 - t_phase1),
    case_drop = unname(t_phase3 - t_phase2),
    statistics = unname(t_end - t_phase3)
  )

  # ---- Assemble result ----
  result <- list(
    original_pcor          = pcor_orig,
    original_precision     = Wi_orig,
    original_centrality    = cent_orig,
    original_predictability = pred_orig,
    edge_ci                = edge_ci_df,
    edge_inclusion         = edge_inclusion,
    thresholded_pcor       = thresholded,
    centrality_ci          = centrality_ci,
    cs_coefficient         = cs_coef,
    cs_data                = cs_data,
    edge_diff_p            = edge_diff_p,
    centrality_diff_p      = centrality_diff_p,
    predictability_ci      = pred_ci_df,
    boot_edges             = boot_edges,
    boot_centrality        = boot_centrality,
    boot_predictability    = boot_pred,
    nodes                  = nodes,
    n                      = n,
    p                      = p,
    iter                   = iter,
    cs_iter                = cs_iter,
    cs_drop                = cs_drop,
    alpha                  = alpha,
    gamma                  = gamma,
    nlambda                = nlambda,
    centrality_measures    = centrality,
    cor_method             = cor_method,
    lambda_path            = lambda_path,
    lambda_selected        = lambda_selected,
    timing                 = timing
  )

  class(result) <- "boot_glasso"
  result
}


# ---- Internal helpers ----

#' Single bootstrap iteration: cor -> glassopath -> EBIC -> pcor -> centrality
#' @noRd
.bg_estimate_once <- function(data_mat, p, gamma, nlambda, penalize_diag,
                               cor_method, lambda_path, centrality,
                               centrality_fn = NULL) {
  n_boot <- nrow(data_mat)
  nodes <- colnames(data_mat)

  S_boot <- tryCatch(
    cor(data_mat, method = cor_method),
    error = function(e) NULL
  )
  if (is.null(S_boot)) return(NULL) # nocov

  # Check for NAs in correlation matrix (can happen with zero-variance
  # columns in bootstrap sample)
  if (any(is.na(S_boot))) return(NULL) # nocov

  gp <- tryCatch(
    glasso::glassopath(s = S_boot, rholist = lambda_path, trace = 0,
                       penalize.diagonal = penalize_diag),
    error = function(e) NULL
  )
  if (is.null(gp)) return(NULL) # nocov

  Wi <- .select_ebic_from_path(gp, S_boot, n_boot, gamma, p, lambda_path)
  if (is.null(Wi)) return(NULL) # nocov

  pcor <- .precision_to_pcor(Wi, threshold = 0)
  pcor <- (pcor + t(pcor)) / 2

  # Extract upper triangle
  edges <- pcor[upper.tri(pcor)]

  # Centrality
  cent <- .bg_compute_centrality(pcor, p, nodes, centrality, centrality_fn)

  # Predictability
  pred <- pmin(pmax(1 - 1 / diag(Wi), 0), 1)

  list(edges = edges, centrality = cent, predictability = pred)
}


#' Compute requested centrality measures
#'
#' strength and expected_influence use rowSums (fast, no dependencies).
#' closeness and betweenness require a user-supplied \code{centrality_fn}.
#'
#' @param pcor Partial correlation matrix.
#' @param p Number of nodes.
#' @param nodes Character vector of node names.
#' @param measures Character vector of requested measures.
#' @param centrality_fn Optional function taking a weight matrix and
#'   returning a named list of centrality vectors. Required for
#'   measures other than strength/expected_influence.
#' @noRd
.bg_compute_centrality <- function(pcor, p, nodes, measures,
                                    centrality_fn = NULL) {
  # Built-in measures (no dependencies)
  builtin  <- c("strength", "expected_influence",
                "betweenness", "closeness")
  external <- setdiff(measures, builtin)

  result <- vector("list", length(measures))
  names(result) <- measures

  # Strength measures (rowSums based)
  if ("strength" %in% measures) {
    v <- rowSums(abs(pcor)); names(v) <- nodes
    result[["strength"]] <- v
  }
  if ("expected_influence" %in% measures) {
    v <- rowSums(pcor); names(v) <- nodes
    result[["expected_influence"]] <- v
  }

  # Path-based measures (Floyd-Warshall, no igraph needed)
  path_mat <- pcor
  rownames(path_mat) <- colnames(path_mat) <- nodes
  if ("betweenness" %in% measures) {
    result[["betweenness"]] <- .betweenness(abs(path_mat),
                                             directed = FALSE, invert = TRUE)
  }
  if ("closeness" %in% measures) {
    cl <- .closeness(abs(path_mat), directed = FALSE, invert = TRUE)
    result[["closeness"]] <- cl[["Closeness"]]
  }

  if (length(external) > 0L) {
    if (is.null(centrality_fn)) {
      stop("centrality_fn is required for measures: ",
           paste(external, collapse = ", "),
           ". Provide a function that takes a weight matrix and returns ",
           "a named list of centrality vectors.", call. = FALSE)
    }
    custom <- centrality_fn(pcor)
    if (!is.list(custom)) {
      stop("centrality_fn must return a named list.", call. = FALSE) # nocov
    }
    for (m in external) {
      if (is.null(custom[[m]])) {
        stop("centrality_fn did not return measure '", m, "'.", # nocov start
             call. = FALSE) # nocov end
      }
      v <- custom[[m]]
      names(v) <- nodes
      result[[m]] <- v
    }
  }

  result
}


#' Run case-dropping iterations and compute CS data
#'
#' Matches bootnet's approach: each of cs_iter iterations randomly picks one
#' drop proportion from cs_drop, draws a subsample of that size, estimates
#' the network, and computes centrality. The CS-coefficient is then computed
#' as the maximum drop proportion where >95% of samples have correlation
#' with the original centrality exceeding 0.7.
#'
#' @noRd
.bg_case_drop <- function(data_mat, p, gamma, nlambda, cor_method,
                           lambda_path, centrality, cs_iter, cs_drop,
                           original_centralities, ncores,
                           centrality_fn = NULL) {
  n <- nrow(data_mat)
  nodes <- colnames(data_mat)
  n_drop <- length(cs_drop)

  # Pre-compute subsample sizes for each drop proportion
  sub_sizes <- vapply(cs_drop, function(dp) {
    as.integer(max(floor(n * (1 - dp)), p + 1L))
  }, integer(1))

  # Each iteration randomly picks a drop proportion (bootnet approach)
  cs_fn <- function(j) {
    d_idx <- sample.int(n_drop, 1L)
    n_sub <- sub_sizes[d_idx]
    idx <- sample.int(n, n_sub, replace = FALSE)
    res <- .bg_estimate_once(data_mat[idx, , drop = FALSE], p, gamma,
                              nlambda, FALSE, cor_method, lambda_path,
                              centrality, centrality_fn)
    if (is.null(res)) return(list(d_idx = d_idx, cors = NULL)) # nocov

    # Compute Pearson correlation with original for each measure
    # (matches bootnet::cor0 which uses Pearson by default)
    cors_m <- vapply(centrality, function(m) {
      sub_vals <- res$centrality[[m]]
      orig_vals <- original_centralities[[m]]
      if (any(!is.finite(sub_vals)) || any(!is.finite(orig_vals)) ||
          sd(sub_vals, na.rm = TRUE) == 0 || sd(orig_vals, na.rm = TRUE) == 0) {
        return(0)
      }
      tryCatch(
        cor(sub_vals, orig_vals),
        warning = function(w) NA_real_,
        error = function(e) NA_real_
      )
    }, numeric(1))
    names(cors_m) <- centrality

    list(d_idx = d_idx, cors = cors_m)
  }

  if (ncores > 1L) {
    cs_results <- parallel::mclapply(seq_len(cs_iter), cs_fn,
                                      mc.cores = ncores)
  } else {
    cs_results <- lapply(seq_len(cs_iter), cs_fn)
  }

  # Organize results by drop proportion
  # cors: list by measure, each a list of vectors by drop proportion
  cors <- lapply(centrality, function(m) {
    mat <- matrix(NA_real_, nrow = 0, ncol = n_drop)
    # Collect per-proportion vectors
    per_prop <- vector("list", n_drop)
    for (d_idx in seq_len(n_drop)) per_prop[[d_idx]] <- numeric(0)

    for (res in cs_results) {
      if (!is.null(res$cors)) {
        per_prop[[res$d_idx]] <- c(per_prop[[res$d_idx]], res$cors[[m]])
      }
    }
    per_prop
  })
  names(cors) <- centrality

  # Build cs_data: mean correlation and P(cor > 0.7) for each proportion
  cs_rows <- vector("list", n_drop)
  for (d_idx in seq_len(n_drop)) {
    for (m in centrality) {
      vals <- cors[[m]][[d_idx]]
      mean_cor <- if (length(vals) > 0) mean(vals, na.rm = TRUE) else NA_real_
      prop_above <- if (length(vals) > 0) {
        mean(vals > 0.7, na.rm = TRUE)
      } else {
        NA_real_
      }
      cs_rows[[d_idx]] <- rbind(cs_rows[[d_idx]], data.frame(
        drop_prop = cs_drop[d_idx],
        measure = m,
        mean_cor = mean_cor,
        prop_above = prop_above,
        n_samples = length(vals),
        stringsAsFactors = FALSE
      ))
    }
  }

  cs_data <- do.call(rbind, cs_rows)
  rownames(cs_data) <- NULL

  list(cors = cors, cs_data = cs_data)
}


#' Compute CS-coefficient from case-dropping correlations
#'
#' Matches bootnet::corStability: CS = max drop proportion where the
#' proportion of bootstrap samples with cor > cor_threshold exceeds
#' prop_threshold (default 0.95).
#'
#' @param cors_per_prop List of numeric vectors, one per drop proportion.
#'   Each vector contains per-iteration Spearman correlations.
#' @param cs_drop Numeric vector of drop proportions.
#' @param cor_threshold Numeric. Correlation threshold (default 0.7).
#' @param prop_threshold Numeric. Required proportion above threshold
#'   (default 0.95).
#' @return Numeric scalar: CS-coefficient.
#' @noRd
.bg_cs_coefficient <- function(cors_per_prop, cs_drop,
                                cor_threshold = 0.7,
                                prop_threshold = 0.95) {
  # For each drop proportion, compute proportion of samples exceeding threshold
  prop_above <- vapply(cors_per_prop, function(vals) {
    if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
    mean(vals > cor_threshold, na.rm = TRUE)
  }, numeric(1))

  # CS = max drop proportion where prop_above > prop_threshold
  above <- !is.na(prop_above) & prop_above > prop_threshold
  if (!any(above)) return(0)
  max(cs_drop[above])
}


#' Find the selected lambda value
#' @noRd
.bg_find_selected_lambda <- function(gp, Wi_orig, lambda_path) {
  # Find which lambda in the path produced the selected Wi
  for (k in seq_along(lambda_path)) {
    wi_k <- gp$wi[, , k]
    if (max(abs(wi_k - Wi_orig)) < 1e-8) {
      return(lambda_path[k])
    }
  }
  # Fallback: return NA
  NA_real_
}


#' Pairwise edge difference p-values
#'
#' For each pair of edges (i, j), compute the proportion of bootstrap
#' samples where edge_i > edge_j (and vice versa), return 2-sided p-value.
#'
#' @noRd
.bg_edge_diff_test <- function(boot_edges) {
  n_edges <- ncol(boot_edges)
  # For efficiency, only compute for manageable number of edges
  if (n_edges > 500) {
    # Return NULL for very large networks
    return(NULL)
  }

  p_mat <- matrix(1, nrow = n_edges, ncol = n_edges,
                   dimnames = list(colnames(boot_edges),
                                    colnames(boot_edges)))

  valid <- complete.cases(boot_edges)
  bm <- boot_edges[valid, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 2) return(p_mat)

  for (i in seq_len(n_edges - 1L)) {
    js <- seq(i + 1L, n_edges)
    diffs <- bm[, i] - bm[, js, drop = FALSE]
    p_greater <- colMeans(diffs > 0)
    p_less    <- colMeans(diffs < 0)
    p_vals    <- 2 * pmin(p_greater, p_less)
    p_mat[i, js] <- p_vals
    p_mat[js, i] <- p_vals
  }

  diag(p_mat) <- 0
  p_mat
}


#' Pairwise centrality difference p-values
#' @noRd
.bg_centrality_diff_test <- function(boot_cent) {
  p_nodes <- ncol(boot_cent)
  node_names <- colnames(boot_cent)

  p_mat <- matrix(1, nrow = p_nodes, ncol = p_nodes,
                   dimnames = list(node_names, node_names))

  valid <- complete.cases(boot_cent)
  bm <- boot_cent[valid, , drop = FALSE]
  n_valid <- nrow(bm)

  if (n_valid < 2) return(p_mat)

  for (i in seq_len(p_nodes - 1L)) {
    js <- seq(i + 1L, p_nodes)
    diffs <- bm[, i] - bm[, js, drop = FALSE]
    p_greater <- colMeans(diffs > 0)
    p_less    <- colMeans(diffs < 0)
    p_vals    <- 2 * pmin(p_greater, p_less)
    p_mat[i, js] <- p_vals
    p_mat[js, i] <- p_vals
  }

  diag(p_mat) <- 0
  p_mat
}


#' Threshold network: zero out edges whose CI includes zero
#' @noRd
.bg_threshold_network <- function(pcor, ci_lower_vec, ci_upper_vec, ut) {
  result <- pcor
  # An edge is significant if both CI bounds have the same sign and neither is 0
  not_sig <- (sign(ci_lower_vec) != sign(ci_upper_vec)) |
    (ci_lower_vec == 0) | (ci_upper_vec == 0)

  for (k in which(not_sig)) {
    i <- ut$row_idx[k]
    j <- ut$col_idx[k]
    result[i, j] <- 0
    result[j, i] <- 0
  }

  result
}


#' Pre-compute upper triangle row/col indices
#' @noRd
.bg_upper_tri_indices <- function(p) {
  ut <- which(upper.tri(matrix(0, p, p)), arr.ind = TRUE)
  list(row_idx = ut[, 1], col_idx = ut[, 2])
}


#' CS-coefficient label
#' @noRd
.bg_cs_label <- function(cs_value) {
  if (is.na(cs_value)) return("Unknown")
  if (cs_value >= 0.5) return("Stable")
  if (cs_value >= 0.25) return("Acceptable")
  "Unstable"
}


#' Build edge names from upper triangle indices
#' @noRd
.bg_build_edge_names <- function(nodes, row_idx, col_idx) {
  paste(nodes[row_idx], "--", nodes[col_idx])
}


# ---- S3 Methods ----

#' Print Method for boot_glasso
#'
#' @param x A \code{boot_glasso} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' set.seed(1)
#' dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
#' bg <- boot_glasso(dat, iter = 10, cs_iter = 5, centrality = "strength")
#' print(bg)
#' \donttest{
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' boot <- boot_glasso(as.data.frame(mat), iter = 20, cs_iter = 10,
#'   centrality = "strength", seed = 42)
#' print(boot)
#' }
#'
#' @export
print.boot_glasso <- function(x, ...) {
  cat(sprintf(
    "GLASSO Bootstrap (%d iterations, %d case-drop per proportion)\n",
    x$iter, x$cs_iter
  ))
  cat(sprintf("  Data: %d x %d  |  Alpha: %.2f  |  Gamma: %.2f\n",
              x$n, x$p, x$alpha, x$gamma))

  # Count significant edges (CI excludes zero)
  n_upper <- x$p * (x$p - 1) / 2
  n_sig <- sum(x$thresholded_pcor[upper.tri(x$thresholded_pcor)] != 0)
  cat(sprintf("  Edges: %d/%d significant (CI excludes zero)\n",
              n_sig, n_upper))
  cat(sprintf("  Mean inclusion probability: %.2f\n",
              mean(x$edge_inclusion)))

  cat("\n  Centrality Stability (CS-coefficient):\n")
  for (m in x$centrality_measures) {
    cs_val <- x$cs_coefficient[[m]]
    label <- .bg_cs_label(cs_val)
    cat(sprintf("    %-22s %.2f [%s]\n", paste0(m, ":"), cs_val, label))
  }

  # Edge differences
  if (!is.null(x$edge_diff_p)) {
    n_edge_pairs <- sum(upper.tri(x$edge_diff_p))
    n_diff_sig <- sum(x$edge_diff_p[upper.tri(x$edge_diff_p)] < x$alpha)
    cat(sprintf("\n  Edge differences: %d/%d pairs significantly different\n",
                n_diff_sig, n_edge_pairs))
  }

  cat(sprintf("  Timing: %.1fs (bootstrap: %.1fs, case-drop: %.1fs)\n",
              x$timing["total"], x$timing["bootstrap"],
              x$timing["case_drop"]))

  invisible(x)
}


#' Summary Method for boot_glasso
#'
#' @param object A \code{boot_glasso} object.
#' @param type Character. Summary type: \code{"edges"} (default),
#'   \code{"centrality"}, \code{"cs"}, \code{"predictability"}, or
#'   \code{"all"}.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame or list of data frames depending on \code{type}.
#'
#' @examples
#' set.seed(1)
#' dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
#' bg <- boot_glasso(dat, iter = 10, cs_iter = 5, centrality = "strength")
#' summary(bg, type = "edges")
#' \donttest{
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' boot <- boot_glasso(as.data.frame(mat), iter = 20, cs_iter = 10,
#'   centrality = "strength", seed = 42)
#' summary(boot, type = "edges")
#' }
#'
#' @export
summary.boot_glasso <- function(object, type = "edges", ...) {
  type <- match.arg(type, c("edges", "centrality", "cs", "predictability",
                             "all"))

  switch(type,
    edges = {
      df <- object$edge_ci
      df[order(-abs(df$weight)), ]
    },
    centrality = {
      object$centrality_ci
    },
    cs = {
      object$cs_data
    },
    predictability = {
      object$predictability_ci
    },
    all = {
      list(
        edges = {
          df <- object$edge_ci
          df[order(-abs(df$weight)), ]
        },
        centrality = object$centrality_ci,
        cs = object$cs_data,
        predictability = object$predictability_ci
      )
    }
  )
}


#' Plot Method for boot_glasso
#'
#' @description
#' Plots bootstrap results for GLASSO networks.
#'
#' @param x A \code{boot_glasso} object.
#' @param type Character. Plot type: \code{"edges"} (default),
#'   \code{"stability"}, \code{"edge_diff"}, \code{"centrality_diff"},
#'   or \code{"inclusion"}.
#' @param measure Character. Centrality measure for
#'   \code{type = "centrality_diff"} (default: first available measure).
#' @param ... Additional arguments passed to plotting functions. For
#'   \code{type = "edge_diff"} and \code{type = "centrality_diff"},
#'   accepts \code{order}: \code{"sample"} (default, sorted by value)
#'   or \code{"id"} (alphabetical).
#'
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' set.seed(1)
#' dat <- as.data.frame(matrix(rnorm(60), ncol = 3))
#' bg <- boot_glasso(dat, iter = 10, cs_iter = 5, centrality = "strength")
#' plot(bg, type = "edges")
#' \donttest{
#' set.seed(42)
#' mat <- matrix(rnorm(60), ncol = 4)
#' colnames(mat) <- LETTERS[1:4]
#' boot <- boot_glasso(as.data.frame(mat), iter = 20, cs_iter = 10,
#'   centrality = "strength", seed = 42)
#' plot(boot, type = "edges")
#' }
#'
#' @export
plot.boot_glasso <- function(x, type = "edges", measure = NULL, ...) {
  type <- match.arg(type, c("edges", "stability", "edge_diff",
                              "centrality_diff", "inclusion"))

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for this plot type.", call. = FALSE) # nocov
  }
  switch(type,
    edges           = .bg_plot_edges(x, ...),
    stability       = .bg_plot_stability(x, ...),
    edge_diff       = .bg_plot_edge_diff(x, ...),
    centrality_diff = .bg_plot_centrality_diff(x, measure, ...),
    inclusion       = .bg_plot_inclusion(x, ...)
  )
}


#' Edge CI plot: sorted pointrange
#' @noRd
.bg_plot_edges <- function(x, ...) {
  df <- x$edge_ci
  df <- df[order(abs(df$weight)), ]
  df$edge <- factor(df$edge, levels = df$edge)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$weight, y = .data$edge)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      size = 0.3
    ) +
    ggplot2::labs(
      title = "Edge Weight Bootstrap CIs",
      x = "Edge Weight",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' CS stability plot: drop curves with threshold line
#' @noRd
.bg_plot_stability <- function(x, ...) {
  cs_data <- x$cs_data

  ggplot2::ggplot(cs_data,
    ggplot2::aes(x = .data$drop_prop, y = .data$mean_cor,
                 color = .data$measure, group = .data$measure)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed",
                         color = "red") +
    ggplot2::scale_x_continuous(breaks = x$cs_drop) +
    ggplot2::labs(
      title = "Centrality Stability (Case-Dropping Bootstrap)",
      x = "Proportion Dropped",
      y = "Mean Correlation with Original",
      color = "Measure"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Edge difference heatmap
#'
#' Upper-triangle heatmap coloured by the mean bootstrap difference between
#' each edge pair. Red = row edge stronger; blue = column edge stronger.
#' Significant cells (p < alpha) are shown at full intensity with p-value
#' stars; non-significant cells are nearly invisible. Diagonal shows the
#' sample edge weight.
#'
#' @param x boot_glasso object.
#' @param order Character. \code{"sample"} (default) orders edges by
#'   absolute sample weight; \code{"id"} orders alphabetically.
#' @param pos_color Colour for positive differences (row edge > col edge).
#'   Default \code{"#C0392B"} (crimson).
#' @param neg_color Colour for negative differences (col edge > row edge).
#'   Default \code{"#2C6E8A"} (teal).
#' @param ... Ignored.
#' @noRd
.bg_plot_edge_diff <- function(x,
                                order     = c("sample", "id"),
                                pos_color = "#C0392B",
                                neg_color = "#2C6E8A",
                                ...) {
  if (is.null(x$edge_diff_p))
    stop("Edge difference test was not computed (too many edges).",
         call. = FALSE)

  order <- match.arg(order)

  p_mat      <- x$edge_diff_p
  alpha      <- x$alpha
  boot_mat   <- x$boot_edges
  edge_names <- colnames(p_mat)
  n_e        <- length(edge_names)

  # Mean bootstrap differences: E[edge_i - edge_j] = mean(edge_i) - mean(edge_j)
  edge_means    <- colMeans(boot_mat, na.rm = TRUE)
  mean_diff_mat <- outer(edge_means, edge_means, `-`)
  rownames(mean_diff_mat) <- colnames(mean_diff_mat) <- edge_names

  # Edge weight lookup (diagonal labels)
  edge_weights <- x$edge_ci$weight
  names(edge_weights) <- x$edge_ci$edge

  # Ordering
  if (order == "sample") {
    ordered_names <- edge_names[order(abs(edge_weights[edge_names]))]
  } else {
    ordered_names <- sort(edge_names)
  }

  # Build full grid
  grid <- expand.grid(edge1 = ordered_names, edge2 = ordered_names,
                      stringsAsFactors = FALSE)
  idx1 <- match(grid$edge1, ordered_names)
  idx2 <- match(grid$edge2, ordered_names)

  is_upper <- idx1 < idx2
  is_diag  <- idx1 == idx2

  grid$p_value   <- p_mat[cbind(grid$edge1, grid$edge2)]
  grid$mean_diff <- mean_diff_mat[cbind(grid$edge1, grid$edge2)]
  grid$sig       <- grid$p_value < alpha

  # Fill value: mean_diff for upper triangle, 0 (midpoint) otherwise
  grid$fill_val  <- ifelse(is_upper | is_diag, grid$mean_diff, NA_real_)
  grid$fill_val[is_diag] <- 0

  # Alpha: sig upper = 1, non-sig upper = 0.08, diagonal = 0.06, lower = 0
  grid$alpha_val <- 0
  grid$alpha_val[is_upper &  grid$sig] <- 1
  grid$alpha_val[is_upper & !grid$sig] <- 0.08
  grid$alpha_val[is_diag]              <- 0.06

  # Significance stars for upper-triangle significant cells
  grid$stars <- ""
  grid$stars[is_upper & grid$sig] <- vapply(
    grid$p_value[is_upper & grid$sig],
    function(p) if (p < 0.001) "***" else if (p < 0.01) "**" else "*",
    character(1)
  )

  # Diagonal weight labels
  grid$label <- ifelse(is_diag,
                        sprintf("%.2f", edge_weights[grid$edge1]), "")

  grid$edge1 <- factor(grid$edge1, levels = ordered_names)
  grid$edge2 <- factor(grid$edge2, levels = rev(ordered_names))

  label_size <- if (n_e <= 10) 3.2 else if (n_e <= 20) 2.5 else 1.8
  star_size  <- label_size * 0.80

  # Subtitle
  ut     <- which(upper.tri(p_mat), arr.ind = TRUE)
  n_sig  <- sum(p_mat[ut] < alpha)
  n_pairs <- nrow(ut)
  subtitle <- sprintf(
    "%d of %d pairs significantly different (p < %s)  |  Red: row > col  |  Blue: col > row",
    n_sig, n_pairs, alpha
  )

  # Abs max for symmetric scale
  diff_max <- max(abs(grid$fill_val), na.rm = TRUE)
  if (diff_max == 0) diff_max <- 0.1

  ggplot2::ggplot(grid, ggplot2::aes(x = .data$edge1, y = .data$edge2)) +
    # White background for all cells
    ggplot2::geom_tile(fill = "white", colour = "#F0F0F0", linewidth = 0.25) +
    # Coloured tiles (significant: full; non-sig: faint; diagonal: near-white)
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$fill_val, alpha = .data$alpha_val),
      colour = "white", linewidth = 0.25
    ) +
    # Stars on significant cells
    ggplot2::geom_text(
      data = grid[is_upper & grid$sig, ],
      ggplot2::aes(label = .data$stars),
      size = star_size, colour = "white", fontface = "bold"
    ) +
    # Edge weight on diagonal
    ggplot2::geom_text(
      data = grid[is_diag, ],
      ggplot2::aes(label = .data$label),
      size = label_size * 0.88, colour = "#444444"
    ) +
    ggplot2::scale_fill_gradient2(
      low      = neg_color,
      high     = pos_color,
      mid      = "white",
      midpoint = 0,
      limits   = c(-diff_max, diff_max),
      na.value = "transparent",
      name     = "Mean\ndifference"
    ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(
      title    = "Bootstrap Edge Difference Test",
      subtitle = subtitle,
      x = NULL, y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(
        angle  = 45, hjust = 1, vjust = 1,
        size   = label_size * 2.8, colour = "#333333"
      ),
      axis.text.y      = ggplot2::element_text(
        size   = label_size * 2.8, colour = "#333333"
      ),
      plot.title       = ggplot2::element_text(
        size = 13, face = "bold", colour = "#1A1A1A",
        margin = ggplot2::margin(b = 4)
      ),
      plot.subtitle    = ggplot2::element_text(
        size = 9, colour = "#666666",
        margin = ggplot2::margin(b = 10)
      ),
      legend.title     = ggplot2::element_text(size = 8.5, colour = "#444444"),
      legend.text      = ggplot2::element_text(size = 8, colour = "#555555"),
      legend.key.width = ggplot2::unit(0.5, "cm"),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin      = ggplot2::margin(12, 12, 12, 12)
    )
}


#' Centrality difference heatmap
#'
#' Upper triangle heatmap with red (significant) / green (non-significant)
#' tiles and sample centrality values on the diagonal.
#'
#' @param x boot_glasso object.
#' @param measure Character. Centrality measure to plot.
#' @param order Character. \code{"sample"} (default) orders nodes by
#'   centrality value (ascending, highest at top-right); \code{"id"}
#'   orders alphabetically.
#' @param ... Ignored.
#' @noRd
.bg_plot_centrality_diff <- function(x, measure = NULL,
                                      order = c("sample", "id"), ...) {
  if (is.null(measure)) measure <- x$centrality_measures[1]
  measure <- match.arg(measure, x$centrality_measures)
  order <- match.arg(order)

  p_mat <- x$centrality_diff_p[[measure]]
  alpha <- x$alpha
  node_names <- colnames(p_mat)
  cent_vals <- x$original_centrality[[measure]]

  # Determine ordering
  if (order == "sample") {
    ordered_names <- node_names[order(cent_vals[node_names])]
  } else {
    ordered_names <- sort(node_names)
  }

  # Build grid: upper triangle + diagonal only
  grid <- expand.grid(
    node1 = ordered_names,
    node2 = ordered_names,
    stringsAsFactors = FALSE
  )
  idx1 <- match(grid$node1, ordered_names)
  idx2 <- match(grid$node2, ordered_names)

  grid$p_value <- vapply(seq_len(nrow(grid)), function(i) {
    p_mat[grid$node1[i], grid$node2[i]]
  }, numeric(1))

  is_diag <- idx1 == idx2
  is_upper <- idx1 < idx2
  is_lower <- idx1 > idx2

  grid$fill <- ifelse(
    is_diag, "diagonal",
    ifelse(is_lower, "blank",
      ifelse(grid$p_value < alpha, "significant", "non-significant")
    )
  )

  # Diagonal labels: sample centrality value
  grid$label <- ifelse(is_diag,
    sprintf("%.2f", cent_vals[grid$node1]), "")

  grid$node1 <- factor(grid$node1, levels = ordered_names)
  grid$node2 <- factor(grid$node2, levels = rev(ordered_names))

  n_p <- length(ordered_names)
  label_size <- if (n_p <= 8) 3.5 else if (n_p <= 15) 2.5 else 1.8

  # Count significant pairs
  ut <- which(upper.tri(p_mat), arr.ind = TRUE)
  n_sig <- sum(p_mat[ut] < alpha)
  n_pairs <- nrow(ut)
  subtitle <- sprintf("%d of %d pairs significantly different (p < %s)",
                       n_sig, n_pairs, alpha)

  ggplot2::ggplot(grid,
    ggplot2::aes(x = .data$node1, y = .data$node2, fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$label),
      size = label_size, colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "significant"     = "#E74C3C",
        "non-significant" = "#2ECC71",
        "diagonal"        = "white",
        "blank"           = "white"
      ),
      labels = c(
        "significant"     = sprintf("p < %s", alpha),
        "non-significant" = sprintf("p >= %s", alpha)
      ),
      breaks = c("significant", "non-significant"),
      name = NULL
    ) +
    ggplot2::labs(
      title = sprintf("Centrality Difference Test (%s)", measure),
      subtitle = subtitle, x = NULL, y = NULL
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text()
    )
}


#' Inclusion probability bar plot
#' @noRd
.bg_plot_inclusion <- function(x, ...) {
  df <- data.frame(
    edge = names(x$edge_inclusion),
    inclusion = as.numeric(x$edge_inclusion),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$inclusion), ]
  df$edge <- factor(df$edge, levels = df$edge)

  ggplot2::ggplot(df,
    ggplot2::aes(x = .data$inclusion, y = .data$edge)) +
    ggplot2::geom_col(fill = "#377EB8") +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed",
                         color = "grey50") +
    ggplot2::labs(
      title = "Edge Inclusion Probability",
      x = "Inclusion Probability",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


