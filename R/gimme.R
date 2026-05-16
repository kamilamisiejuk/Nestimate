# ---- Group Iterative Multiple Model Estimation (GIMME) ----

#' GIMME: Group Iterative Multiple Model Estimation
#'
#' @description
#' Estimates person-specific directed networks from intensive longitudinal data
#' using the unified Structural Equation Modeling (uSEM) framework. Implements
#' a data-driven search that identifies:
#' \enumerate{
#'   \item \strong{Group-level paths}: Directed edges present for a majority
#'     (default 75\%) of individuals.
#'   \item \strong{Individual-level paths}: Additional edges specific to each
#'     person, found after group paths are established.
#' }
#'
#' Uses \code{lavaan} for SEM estimation and modification indices.
#' Accepts a single data frame with an ID column (not CSV directories).
#'
#' @param data A \code{data.frame} in long format with columns for person ID,
#'   time-varying variables, and optionally a time/beep column.
#' @param vars Character vector of variable names to model.
#' @param id Character string naming the person-ID column.
#' @param time Character string naming the time/order column, or \code{NULL}.
#'   When provided, data is sorted by \code{id} then \code{time} before lagging.
#' @param ar Logical. If \code{TRUE} (default), autoregressive paths
#'   (each variable predicting itself at lag 1) are included as fixed paths.
#' @param standardize Logical. If \code{TRUE} (default \code{FALSE}),
#'   variables are standardized per person before estimation.
#' @param groupcutoff Numeric between 0 and 1. Proportion of individuals for
#'   whom a path must be significant to be added at group level.
#'   Default \code{0.75}.
#' @param subcutoff Numeric. Not used (reserved for future subgrouping).
#' @param paths Character vector of lavaan-syntax paths to force into the model
#'   (e.g., \code{"V2~V1lag"}). Default \code{NULL}.
#' @param exogenous Character vector of variable names to treat as exogenous.
#'   Default \code{NULL}.
#' @param hybrid Logical. If \code{TRUE}, also searches residual covariances.
#'   Default \code{FALSE}.
#' @param rmsea_cutoff Numeric. RMSEA threshold for excellent fit (default 0.05).
#' @param srmr_cutoff Numeric. SRMR threshold for excellent fit (default 0.05).
#' @param nnfi_cutoff Numeric. NNFI/TLI threshold for excellent fit (default 0.95).
#' @param cfi_cutoff Numeric. CFI threshold for excellent fit (default 0.95).
#' @param n_excellent Integer. Number of fit indices that must be excellent to
#'   stop individual search. Default \code{2}.
#' @param seed Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return An S3 object of class \code{"net_gimme"} containing:
#' \describe{
#'   \item{\code{temporal}}{p x p matrix of group-level temporal (lagged)
#'     path counts -- entry \code{[i,j]} = number of individuals with path j(t-1)->i(t).}
#'   \item{\code{contemporaneous}}{p x p matrix of group-level contemporaneous
#'     path counts -- entry \code{[i,j]} = number of individuals with path j(t)->i(t).}
#'   \item{\code{coefs}}{List of per-person p x 2p coefficient matrices
#'     (rows = endogenous, cols = \code{[lagged, contemporaneous]}).}
#'   \item{\code{psi}}{List of per-person residual covariance matrices.}
#'   \item{\code{fit}}{Data frame of per-person fit indices (chisq, df, pvalue,
#'     rmsea, srmr, nnfi, cfi, bic, aic, logl, status).}
#'   \item{\code{path_counts}}{p x 2p matrix: how many individuals have each path.}
#'   \item{\code{paths}}{List of per-person character vectors of lavaan path syntax.}
#'   \item{\code{group_paths}}{Character vector of group-level paths found.}
#'   \item{\code{individual_paths}}{List of per-person character vectors of
#'     individual-level paths (beyond group).}
#'   \item{\code{syntax}}{List of per-person full lavaan syntax strings.}
#'   \item{\code{labels}}{Character vector of variable names.}
#'   \item{\code{n_subjects}}{Integer. Number of individuals.}
#'   \item{\code{n_obs}}{Integer vector. Time points per individual.}
#'   \item{\code{config}}{List of configuration parameters.}
#' }
#'
#' @examples
#' \donttest{
#' # Create simple panel data (3 subjects, 4 variables, 50 time points).
#' set.seed(42)
#' n_sub <- 3; n_t <- 50; vars <- paste0("V", 1:4)
#' rows <- lapply(seq_len(n_sub), function(i) {
#'   d <- as.data.frame(matrix(rnorm(n_t * 4), ncol = 4))
#'   names(d) <- vars; d$id <- i; d
#' })
#' panel <- do.call(rbind, rows)
#' res <- build_gimme(panel, vars = vars, id = "id")
#' print(res)
#' }
#'
#' @seealso \code{\link{build_network}}
#'
#' @export
build_gimme <- function(data,
                        vars,
                        id,
                        time = NULL,
                        ar = TRUE,
                        standardize = FALSE,
                        groupcutoff = 0.75,
                        subcutoff = 0.50,
                        paths = NULL,
                        exogenous = NULL,
                        hybrid = FALSE,
                        rmsea_cutoff = 0.05,
                        srmr_cutoff = 0.05,
                        nnfi_cutoff = 0.95,
                        cfi_cutoff = 0.95,
                        n_excellent = 2L,
                        seed = NULL) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Package 'lavaan' is required for build_gimme(). ", # nocov
         "Install it with install.packages('lavaan').", call. = FALSE) # nocov
  }

  if (!is.null(seed)) set.seed(seed)

  # --- Input validation ---
  stopifnot(is.data.frame(data))
  stopifnot(is.character(vars), length(vars) >= 2L)
  stopifnot(is.character(id), length(id) == 1L, id %in% names(data))
  if (!all(vars %in% names(data))) {
    missing_v <- setdiff(vars, names(data))
    stop("Variables not found in data: ", paste(missing_v, collapse = ", "),
         call. = FALSE)
  }
  stopifnot(is.numeric(groupcutoff), groupcutoff > 0, groupcutoff <= 1)

  # --- Prepare per-person data ---
  ts_list <- .gimme_prepare_data(data, vars, id, time, standardize, exogenous)
  n_subj <- length(ts_list)
  varnames <- vars
  p <- length(varnames)

  if (n_subj < 2L) {
    stop("build_gimme() requires at least 2 individuals.", call. = FALSE)
  }

  # --- Build lavaan syntax components ---
  lag_names <- paste0(varnames, "lag")
  endo_names <- varnames
  exog_names <- lag_names

  syntax_info <- .gimme_build_syntax(varnames, lag_names, endo_names,
                                      exog_names, ar, paths, exogenous,
                                      hybrid)

  base_syntax <- syntax_info$base_syntax
  candidate_paths <- syntax_info$candidate_paths
  candidate_corr <- syntax_info$candidate_corr
  fixed_paths <- syntax_info$fixed_paths

  if (hybrid) {
    elig_paths <- c(candidate_paths, candidate_corr)
  } else {
    elig_paths <- candidate_paths
  }

  # --- Group-level search ---
  # Bonferroni correction for group search
  grp_cutoff <- stats::qchisq(1 - 0.05 / n_subj, 1)
  grp_z_cutoff <- abs(stats::qnorm(0.025 / n_subj))

  group_paths <- character(0)
  group_search_converged <- FALSE

  while (!group_search_converged) {
    current_syntax <- c(base_syntax, group_paths)

    # Fit model per person, collect modification indices
    mi_list <- lapply(ts_list, function(d) {
      .gimme_fit_and_mi(current_syntax, d, elig_paths)
    })

    # Select best path
    add_path <- .gimme_select_path(mi_list, elig_paths, groupcutoff,
                                    n_subj, grp_cutoff, hybrid)

    if (is.na(add_path)) {
      group_search_converged <- TRUE
    } else {
      group_paths <- c(group_paths, add_path)
    }
  }

  # --- Group-level pruning ---
  if (length(group_paths) > 0) {
    group_paths <- .gimme_prune_paths(base_syntax, group_paths, ts_list,
                                       n_subj, groupcutoff, grp_z_cutoff)
  }

  # --- Individual-level search ---
  # Bonferroni across eligible paths for individual level
  ind_cutoff <- stats::qchisq(1 - 0.05 / length(elig_paths), 1)
  ind_z_cutoff <- abs(stats::qnorm(0.05 / length(elig_paths)))

  fit_indices <- c("rmsea_cutoff" = rmsea_cutoff, "srmr_cutoff" = srmr_cutoff,
                   "nnfi_cutoff" = nnfi_cutoff, "cfi_cutoff" = cfi_cutoff)

  ind_results <- lapply(seq_len(n_subj), function(k) {
    .gimme_individual_search(
      base_syntax = base_syntax,
      group_paths = group_paths,
      data_k = ts_list[[k]],
      elig_paths = elig_paths,
      ind_cutoff = ind_cutoff,
      ind_z_cutoff = ind_z_cutoff,
      fit_indices = fit_indices,
      n_excellent = n_excellent,
      hybrid = hybrid,
      endo_names = endo_names,
      lag_names = lag_names
    )
  })
  names(ind_results) <- names(ts_list)

  # --- Extract results ---
  result <- .gimme_extract_results(ind_results, varnames, lag_names,
                                    group_paths, base_syntax, fixed_paths,
                                    n_subj, ts_list, hybrid)

  result$labels <- varnames
  result$n_subjects <- n_subj
  result$n_obs <- vapply(ts_list, nrow, integer(1))
  result$config <- list(
    ar = ar, standardize = standardize, groupcutoff = groupcutoff,
    hybrid = hybrid, rmsea_cutoff = rmsea_cutoff, srmr_cutoff = srmr_cutoff,
    nnfi_cutoff = nnfi_cutoff, cfi_cutoff = cfi_cutoff,
    n_excellent = n_excellent, exogenous = exogenous,
    fixed_paths = fixed_paths, seed = seed
  )

  class(result) <- "net_gimme"
  result
}


# ============================================================================
# Data preparation
# ============================================================================

#' Prepare per-person time series with lagged columns
#' @noRd
.gimme_prepare_data <- function(data, vars, id, time, standardize, exogenous) {
  ids <- unique(data[[id]])

  ts_list <- lapply(ids, function(pid) {
    d <- data[data[[id]] == pid, , drop = FALSE]

    # Sort by time if provided
    if (!is.null(time)) {
      d <- d[order(d[[time]]), , drop = FALSE]
    }

    # Extract variables
    mat <- d[, vars, drop = FALSE]

    # Standardize per person if requested
    if (standardize) {
      mat <- as.data.frame(lapply(mat, function(x) {
        s <- stats::sd(x, na.rm = TRUE)
        if (is.na(s) || s == 0) return(x) # nocov
        (x - mean(x, na.rm = TRUE)) / s
      }))
    }

    n <- nrow(mat)
    if (n < 3L) return(NULL)

    # Create lagged version
    lag_mat <- mat[-n, , drop = FALSE]
    colnames(lag_mat) <- paste0(colnames(lag_mat), "lag")

    # Current (drop first row)
    cur_mat <- mat[-1L, , drop = FALSE]

    # Combine: current + lagged
    result <- cbind(cur_mat, lag_mat)
    rownames(result) <- NULL
    result
  })

  names(ts_list) <- as.character(ids)
  # Remove NULLs (subjects with too few observations)
  ts_list <- ts_list[!vapply(ts_list, is.null, logical(1))]

  if (length(ts_list) == 0) {
    stop("No subjects have enough time points (minimum 3).", call. = FALSE)
  }

  ts_list
}


# ============================================================================
# Syntax building
# ============================================================================

#' Build base lavaan syntax and candidate path lists
#' @noRd
.gimme_build_syntax <- function(varnames, lag_names, endo_names, exog_names,
                                 ar, paths, exogenous, hybrid) {
  # Endogenous variances and intercepts
  var_endo <- paste0(endo_names, "~~", endo_names)
  int_endo <- paste0(endo_names, "~1")

  # Exogenous (lagged) covariances and intercepts
  exog_pairs <- outer(exog_names, exog_names, function(x, y) paste0(x, "~~", y))
  cov_exog <- exog_pairs[lower.tri(exog_pairs, diag = TRUE)]
  int_exog <- paste0(exog_names, "~1")

  # Nonsense paths: lagged cannot be predicted by current
  nons_reg <- c(t(outer(exog_names, endo_names, function(x, y) {
    paste0(x, "~0*", y)
  })))

  # Fixed paths (AR if requested + user-specified)
  fixed_paths <- paths
  if (ar) {
    ar_paths <- paste0(varnames, "~", varnames, "lag")
    fixed_paths <- c(fixed_paths, ar_paths)
  }

  # All possible directed paths: endo ~ (endo + exog)
  # Contemporaneous: endo_i ~ endo_j (i != j)
  # Lagged (cross): endo_i ~ exog_j (where j is not i's own lag, if ar handles that)
  all_poss <- c(t(outer(endo_names, c(endo_names, exog_names), function(x, y) {
    paste0(x, "~", y)
  })))
  # Remove self-regression (endo_i ~ endo_i) -- these are variances
  self_reg <- paste0(endo_names, "~", endo_names)
  all_poss <- setdiff(all_poss, self_reg)

  # All possible residual covariances
  corr_pairs <- outer(endo_names, endo_names, function(x, y) paste0(x, "~~", y))
  all_corr <- c(corr_pairs[lower.tri(corr_pairs)],
                corr_pairs[upper.tri(corr_pairs)])

  # Candidate paths = all possible minus fixed and nonsense
  candidate_paths <- setdiff(all_poss, c(fixed_paths, nons_reg))
  candidate_corr <- all_corr

  base_syntax <- c(var_endo, int_endo, cov_exog, int_exog, nons_reg,
                   fixed_paths)

  list(
    base_syntax = base_syntax,
    candidate_paths = candidate_paths,
    candidate_corr = candidate_corr,
    fixed_paths = fixed_paths
  )
}


# ============================================================================
# Lavaan fitting helpers
# ============================================================================

#' Fit lavaan model and extract modification indices for eligible paths
#' @noRd
.gimme_fit_and_mi <- function(syntax, data_k, elig_paths) {
  fit <- .gimme_fit_lavaan(syntax, data_k)

  if (is.null(fit)) return(NA) # nocov start
  if (!lavaan::lavInspect(fit, "converged")) return(NA)
  if (any(is.na(lavaan::lavInspect(fit, what = "list")$se))) return(NA) # nocov end

  mi <- tryCatch({
    mis <- lavaan::modindices(fit, standardized = FALSE, sort. = FALSE)
    mis$param <- paste0(mis$lhs, mis$op, mis$rhs)
    mis[mis$param %in% elig_paths, , drop = FALSE]
  }, error = function(e) NA)

  mi
}


#' Fit lavaan model and extract z-values for specified paths
#' @noRd
.gimme_fit_and_z <- function(syntax, data_k, elig_paths) {
  fit <- .gimme_fit_lavaan(syntax, data_k)

  if (is.null(fit)) return(NA) # nocov start
  if (!lavaan::lavInspect(fit, "converged")) return(NA) # nocov end

  # Use standardizedSolution for z-values (matches gimme's return.zs)
  ss <- tryCatch(lavaan::standardizedSolution(fit), error = function(e) NULL)
  if (is.null(ss)) return(NA) # nocov

  ss$param <- paste0(ss$lhs, ss$op, ss$rhs)
  ss[ss$param %in% elig_paths, , drop = FALSE]
}


#' Fit final lavaan model and extract all results
#' @noRd
.gimme_fit_final <- function(syntax, data_k, varnames, lag_names) {
  fit <- .gimme_fit_lavaan(syntax, data_k)

  if (is.null(fit)) { # nocov start
    p <- length(varnames)
    return(list(
      coefs = matrix(0, p, 2 * p,
                     dimnames = list(varnames, c(lag_names, varnames))),
      psi = matrix(0, p, p, dimnames = list(varnames, varnames)),
      fit_indices = data.frame(
        chisq = NA, df = NA, pvalue = NA, rmsea = NA, srmr = NA,
        nnfi = NA, cfi = NA, bic = NA, aic = NA, logl = NA
      ),
      status = "failed to converge"
    )) # nocov end
  }

  converged <- lavaan::lavInspect(fit, "converged")

  p <- length(varnames)
  all_names <- c(lag_names, varnames)

  # Extract standardized coefficient matrix (betas) -- matches gimme's output
  std_est <- tryCatch(lavaan::lavInspect(fit, "std"), error = function(e) NULL)

  coef_mat <- matrix(0, p, 2 * p, dimnames = list(varnames, all_names))
  psi_mat <- matrix(0, 2 * p, 2 * p, dimnames = list(all_names, all_names))

  if (!is.null(std_est)) {
    beta <- std_est$beta
    avail_rows <- intersect(varnames, rownames(beta))
    avail_cols <- intersect(all_names, colnames(beta))
    coef_mat[avail_rows, avail_cols] <- round(
      beta[avail_rows, avail_cols, drop = FALSE], digits = 4
    )

    psi <- std_est$psi
    avail_psi_r <- intersect(all_names, rownames(psi))
    avail_psi_c <- intersect(all_names, colnames(psi))
    psi_mat[avail_psi_r, avail_psi_c] <- round(
      psi[avail_psi_r, avail_psi_c, drop = FALSE], digits = 4
    )
  }

  # Extract fit indices
  fi <- tryCatch(
    lavaan::fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "srmr",
                                "nnfi", "cfi", "bic", "aic", "logl")),
    error = function(e) rep(NA, 10)
  )
  fi_df <- as.data.frame(as.list(fi))

  # Status
  status <- if (converged) "converged normally" else "failed to converge"

  list(
    coefs = coef_mat,
    psi = psi_mat,
    fit_indices = fi_df,
    status = status
  )
}


# ============================================================================
# Group search
# ============================================================================

#' Select best candidate path from modification indices across subjects
#' @noRd
.gimme_select_path <- function(mi_list, elig_paths, prop_cutoff, n_subj,
                                chisq_cutoff, hybrid) {
  # Remove NAs (non-converged subjects)
  mi_valid <- mi_list[!vapply(mi_list, function(x) is.atomic(x) && length(x) == 1L && is.na(x), logical(1))]
  n_converge <- length(mi_valid)

  if (n_converge <= (n_subj / 2)) return(NA_character_)
  if (n_converge == 0) return(NA_character_) # nocov

  # Combine all modification indices
  mi_all <- do.call(rbind, mi_valid)
  if (is.null(mi_all) || nrow(mi_all) == 0) return(NA_character_) # nocov

  mi_all$param <- paste0(mi_all$lhs, mi_all$op, mi_all$rhs)
  mi_all <- mi_all[mi_all$param %in% elig_paths, , drop = FALSE]
  if (nrow(mi_all) == 0) return(NA_character_) # nocov

  # Count significant MIs per path
  mi_all$sig <- ifelse(mi_all$mi >= chisq_cutoff, 1L, 0L)

  # Aggregate per path
  agg <- stats::aggregate(
    cbind(mi, sig) ~ param,
    data = mi_all,
    FUN = function(x) c(sum = sum(x), count_or_mean = NA)
  )
  # Redo properly
  param_stats <- data.frame(
    param = unique(mi_all$param),
    stringsAsFactors = FALSE
  )
  param_stats$sum_mi <- vapply(param_stats$param, function(p) {
    sum(mi_all$mi[mi_all$param == p])
  }, numeric(1))
  param_stats$count_sig <- vapply(param_stats$param, function(p) {
    sum(mi_all$sig[mi_all$param == p])
  }, numeric(1))
  param_stats$mean_mi <- vapply(param_stats$param, function(p) {
    mean(mi_all$mi[mi_all$param == p])
  }, numeric(1))

  # Sort: most people significant first, then highest mean MI (matches gimme)
  param_stats <- param_stats[order(-param_stats$count_sig,
                                    -param_stats$mean_mi), , drop = FALSE]

  # Check if top path meets group cutoff
  if (param_stats$count_sig[1] > (prop_cutoff * n_converge)) {
    return(param_stats$param[1])
  }

  NA_character_
}


# ============================================================================
# Group pruning
# ============================================================================

#' Prune group paths by checking z-values across subjects
#' @noRd
.gimme_prune_paths <- function(base_syntax, group_paths, ts_list,
                                n_subj, prop_cutoff, z_cutoff) {
  pruning <- TRUE

  while (pruning) {
    current_syntax <- c(base_syntax, group_paths)

    # Get z-values per person
    z_list <- lapply(ts_list, function(d) {
      .gimme_fit_and_z(current_syntax, d, group_paths)
    })

    # Find weakest path
    drop_path <- .gimme_find_weakest(z_list, group_paths, prop_cutoff,
                                      n_subj, z_cutoff)

    if (is.na(drop_path)) {
      pruning <- FALSE
    } else {
      group_paths <- setdiff(group_paths, drop_path) # nocov
      if (length(group_paths) == 0) { # nocov start
        pruning <- FALSE # nocov end
      }
    }
  }

  group_paths
}


#' Find the weakest path that should be pruned
#' @noRd
.gimme_find_weakest <- function(z_list, elig_paths, prop_cutoff, n_subj,
                                 z_cutoff) {
  z_valid <- z_list[!vapply(z_list, function(x) is.atomic(x) && length(x) == 1L && is.na(x), logical(1))]
  n_converge <- length(z_valid)
  if (n_converge == 0) return(NA_character_)

  z_all <- do.call(rbind, z_valid)
  if (is.null(z_all) || nrow(z_all) == 0) return(NA_character_)

  z_all$param <- paste0(z_all$lhs, z_all$op, z_all$rhs)
  z_all <- z_all[z_all$param %in% elig_paths, , drop = FALSE]
  if (nrow(z_all) == 0) return(NA_character_)

  # Count non-significant z-values per path
  z_all$nonsig <- ifelse(abs(z_all$z) < z_cutoff, 1L, 0L)

  param_stats <- data.frame(
    param = unique(z_all$param),
    stringsAsFactors = FALSE
  )
  param_stats$count_nonsig <- vapply(param_stats$param, function(p) {
    sum(z_all$nonsig[z_all$param == p])
  }, numeric(1))
  param_stats$mean_abs_z <- vapply(param_stats$param, function(p) {
    mean(abs(z_all$z[z_all$param == p]))
  }, numeric(1))

  # Sort: most non-significant first, then lowest mean |z|
  param_stats <- param_stats[order(-param_stats$count_nonsig,
                                    param_stats$mean_abs_z), , drop = FALSE]

  # Prune if the weakest path is non-significant for > (1 - prop_cutoff) of subjects
  if (param_stats$count_nonsig[1] > ((1 - prop_cutoff) * n_converge)) {
    return(param_stats$param[1]) # nocov
  }

  NA_character_
}


# ============================================================================
# Individual search
# ============================================================================

#' Run individual-level path search for one person
#' @noRd
.gimme_individual_search <- function(base_syntax, group_paths, data_k,
                                      elig_paths, ind_cutoff, ind_z_cutoff,
                                      fit_indices, n_excellent, hybrid,
                                      endo_names, lag_names) {
  ind_paths <- character(0)
  nonconv_path <- character(0)  # paths that caused instability
  dropped_param <- character(0) # paths dropped by z-pruning

  # --- Phase 1: Forward search ---
  ind_paths <- .gimme_ind_forward_search(
    base_syntax, group_paths, ind_paths, data_k, elig_paths,
    ind_cutoff, fit_indices, n_excellent,
    exclude = character(0)
  )

  # --- Phase 2: Stability + prune + resume cycle ---
  # Mirrors gimme's search.paths.ind post-search loop
  outer_done <- FALSE

  while (!outer_done) {
    # 2a: Check stability -- pop unstable paths
    stable_result <- .gimme_stabilize(
      base_syntax, group_paths, ind_paths, data_k,
      endo_names, lag_names
    )
    ind_paths <- stable_result$ind_paths
    nonconv_path <- c(nonconv_path, stable_result$removed)
    converged <- stable_result$converged
    stable <- stable_result$stable

    if (!converged || !stable) break # nocov

    # 2b: Z-prune individual paths
    pruned_any <- FALSE
    if (length(ind_paths) > 0) {
      prune_done <- FALSE
      while (!prune_done) {
        current_syntax <- c(base_syntax, group_paths, ind_paths)
        z_info <- .gimme_fit_and_z(current_syntax, data_k, ind_paths)

        if (!identical(z_info, NA) && is.data.frame(z_info) && nrow(z_info) > 0) {
          z_info$param <- paste0(z_info$lhs, z_info$op, z_info$rhs)
          # Find weakest non-significant path
          nonsig <- z_info[abs(z_info$z) < ind_z_cutoff, , drop = FALSE]
          if (nrow(nonsig) > 0) {
            weakest <- nonsig$param[which.min(abs(nonsig$z))] # nocov start
            dropped_param <- c(dropped_param, weakest)
            ind_paths <- setdiff(ind_paths, weakest)
            pruned_any <- TRUE # nocov end
          } else {
            prune_done <- TRUE
          }
        } else {
          prune_done <- TRUE # nocov
        }
      }
    }

    # 2c: If stability or z-pruning removed paths, try to resume forward search
    stability_removed <- length(stable_result$removed) > 0
    if (pruned_any || stability_removed) {
      exclude <- c(dropped_param, nonconv_path) # nocov start
      new_paths <- .gimme_ind_forward_search(
        base_syntax, group_paths, ind_paths, data_k,
        elig_paths, ind_cutoff, fit_indices, n_excellent,
        exclude = exclude
      ) # nocov end
      if (length(new_paths) > length(ind_paths)) { # nocov start
        ind_paths <- new_paths
        # Loop back to stability check
        next # nocov end
      }
    }

    outer_done <- TRUE
  }

  list(
    group_paths = group_paths,
    ind_paths = ind_paths,
    full_syntax = c(base_syntax, group_paths, ind_paths)
  )
}


#' Forward search: add individual paths one at a time via MI
#' @noRd
.gimme_ind_forward_search <- function(base_syntax, group_paths, ind_paths,
                                       data_k, elig_paths, ind_cutoff,
                                       fit_indices, n_excellent, exclude) {
  search <- TRUE

  while (search) {
    current_syntax <- c(base_syntax, group_paths, ind_paths)

    # Check if current fit is excellent enough
    fit <- .gimme_fit_lavaan(current_syntax, data_k)

    if (!is.null(fit) && lavaan::lavInspect(fit, "converged")) {
      fi <- tryCatch(
        lavaan::fitMeasures(fit, c("rmsea", "srmr", "nnfi", "cfi")),
        error = function(e) NULL
      )

      if (!is.null(fi)) {
        n_exc <- sum(c(
          !is.na(fi["rmsea"]) && fi["rmsea"] <= fit_indices["rmsea_cutoff"],
          !is.na(fi["srmr"]) && fi["srmr"] <= fit_indices["srmr_cutoff"],
          !is.na(fi["nnfi"]) && fi["nnfi"] >= fit_indices["nnfi_cutoff"],
          !is.na(fi["cfi"]) && fi["cfi"] >= fit_indices["cfi_cutoff"]
        ))

        if (n_exc >= n_excellent) {
          search <- FALSE
          next
        }
      }
    }

    # Get modification indices
    mi <- .gimme_fit_and_mi(current_syntax, data_k, elig_paths)
    if (is.null(mi) || (is.atomic(mi) && length(mi) == 1L && is.na(mi)) || nrow(mi) == 0) { # nocov start
      search <- FALSE
      next # nocov end
    }

    mi$param <- paste0(mi$lhs, mi$op, mi$rhs)
    # Filter to paths not already in model and not excluded
    already_in <- c(group_paths, ind_paths, exclude)
    mi <- mi[!mi$param %in% already_in, , drop = FALSE]
    if (nrow(mi) == 0) { # nocov start
      search <- FALSE
      next # nocov end
    }

    # Find significant path with highest MI
    mi_sig <- mi[mi$mi >= ind_cutoff, , drop = FALSE]
    if (nrow(mi_sig) == 0) {
      search <- FALSE
      next
    }

    mi_sig <- mi_sig[order(-mi_sig$mi), , drop = FALSE]
    ind_paths <- c(ind_paths, mi_sig$param[1])
  }

  ind_paths
}


#' Test whether standardized beta eigenvalues indicate instability
#' @details Matches gimme's testWeights: checks if any eigenvalue of
#'   the contemporaneous or lagged standardized beta block has Re >= 1.
#' @noRd
.gimme_test_weights <- function(fit, endo_names, lag_names) {
  std_beta <- tryCatch(
    lavaan::lavInspect(fit, "std")$beta,
    error = function(e) NULL
  )
  if (is.null(std_beta)) return(TRUE)  # treat errors as unstable

  # Extract and reorder: rows = endogenous, cols = c(lag_names, endo_names)
  coln <- c(lag_names, endo_names)
  avail_rows <- intersect(endo_names, rownames(std_beta))
  avail_cols <- intersect(coln, colnames(std_beta))
  if (length(avail_rows) == 0 || length(avail_cols) == 0) return(TRUE) # nocov

  betas <- round(std_beta[avail_rows, avail_cols, drop = FALSE], digits = 4)
  p <- length(endo_names)

  # Lagged block = first p columns, contemporaneous block = next p columns
  lag_block <- betas[, seq_len(min(p, ncol(betas))), drop = FALSE]
  cont_start <- p + 1L
  cont_end <- min(2L * p, ncol(betas))

  unstable <- any(Re(eigen(lag_block, only.values = TRUE)$values) >= 1)
  if (cont_start <= cont_end) {
    cont_block <- betas[, cont_start:cont_end, drop = FALSE]
    unstable <- unstable || any(Re(eigen(cont_block, only.values = TRUE)$values) >= 1)
  }

  unstable
}


#' Pop unstable individual paths until model is stable
#' @noRd
.gimme_stabilize <- function(base_syntax, group_paths, ind_paths, data_k,
                              endo_names, lag_names) {
  removed <- character(0)

  repeat {
    current_syntax <- c(base_syntax, group_paths, ind_paths)
    fit <- .gimme_fit_lavaan(current_syntax, data_k)

    if (is.null(fit)) { # nocov start
      return(list(ind_paths = ind_paths, removed = removed,
                  converged = FALSE, stable = FALSE))
    } # nocov end

    converged <- lavaan::lavInspect(fit, "converged")
    zero_se <- tryCatch(
      sum(lavaan::lavInspect(fit, "se")$beta, na.rm = TRUE) == 0,
      error = function(e) TRUE
    )
    unstable <- .gimme_test_weights(fit, endo_names, lag_names)

    if (converged && !zero_se && !unstable) {
      return(list(ind_paths = ind_paths, removed = removed,
                  converged = TRUE, stable = TRUE))
    }

    # Model is unstable/non-converged -- pop last individual path
    if (length(ind_paths) == 0) {
      return(list(ind_paths = ind_paths, removed = removed,
                  converged = converged, stable = !unstable))
    }

    removed <- c(removed, ind_paths[length(ind_paths)]) # nocov start
    ind_paths <- ind_paths[-length(ind_paths)] # nocov end
  }
}


#' Fit a lavaan model with standard gimme settings
#' @noRd
.gimme_fit_lavaan <- function(syntax, data_k) {
  tryCatch(
    lavaan::lavaan(
      model = paste(syntax, collapse = "\n"),
      data = data_k,
      model.type = "sem",
      missing = "fiml",
      estimator = "ml",
      int.ov.free = FALSE,
      int.lv.free = TRUE,
      auto.fix.first = TRUE,
      auto.var = TRUE,
      auto.cov.lv.x = TRUE,
      auto.th = TRUE,
      auto.delta = TRUE,
      auto.cov.y = FALSE,
      auto.fix.single = TRUE,
      warn = FALSE
    ),
    error = function(e) NULL
  )
}


# ============================================================================
# Result extraction
# ============================================================================

#' Extract structured results from individual search outputs
#' @noRd
.gimme_extract_results <- function(ind_results, varnames, lag_names,
                                    group_paths, base_syntax, fixed_paths,
                                    n_subj, ts_list, hybrid) {
  p <- length(varnames)
  all_names <- c(lag_names, varnames)
  subj_names <- names(ind_results)

  # Fit final models and extract per-person results
  coefs_list <- vector("list", n_subj)
  psi_list <- vector("list", n_subj)
  fit_list <- vector("list", n_subj)
  syntax_list <- vector("list", n_subj)
  ind_paths_list <- vector("list", n_subj)

  for (k in seq_len(n_subj)) {
    syntax_k <- ind_results[[k]]$full_syntax
    syntax_list[[k]] <- syntax_k
    ind_paths_list[[k]] <- ind_results[[k]]$ind_paths

    final <- .gimme_fit_final(syntax_k, ts_list[[k]], varnames, lag_names)

    coefs_list[[k]] <- final$coefs
    psi_list[[k]] <- final$psi
    fit_list[[k]] <- cbind(
      data.frame(file = subj_names[k], stringsAsFactors = FALSE),
      final$fit_indices,
      data.frame(status = final$status, stringsAsFactors = FALSE)
    )
  }
  names(coefs_list) <- subj_names
  names(psi_list) <- subj_names
  names(syntax_list) <- subj_names
  names(ind_paths_list) <- subj_names

  # Build fit data.frame
  fit_df <- do.call(rbind, fit_list)
  rownames(fit_df) <- NULL

  # Build path count matrices
  path_counts <- matrix(0L, p, 2 * p,
                        dimnames = list(varnames, all_names))
  for (k in seq_len(n_subj)) {
    mat <- coefs_list[[k]]
    path_counts <- path_counts + (mat != 0) * 1L
  }

  fixed_reg <- fixed_paths[grepl("~", fixed_paths, fixed = TRUE) &
                             !grepl("~~", fixed_paths, fixed = TRUE)]
  if (length(fixed_reg) > 0L) {
    fixed_lhs <- sub("\\s*~.*$", "", fixed_reg)
    fixed_rhs <- sub("^.*?~\\s*", "", fixed_reg)
    fixed_rhs <- sub("^[+-]?[0-9.]+\\*", "", fixed_rhs)
    fixed_ok <- fixed_lhs %in% rownames(path_counts) &
      fixed_rhs %in% colnames(path_counts)
    path_counts[cbind(fixed_lhs[fixed_ok], fixed_rhs[fixed_ok])] <- n_subj
  }

  # Separate temporal and contemporaneous count matrices
  temporal_counts <- path_counts[, lag_names, drop = FALSE]
  colnames(temporal_counts) <- varnames  # Remove "lag" suffix for display
  contemp_counts <- path_counts[, varnames, drop = FALSE]

  # Build group-level average coefficient matrices
  temporal_avg <- Reduce("+", lapply(coefs_list, function(m) m[, lag_names, drop = FALSE])) / n_subj
  colnames(temporal_avg) <- varnames
  contemp_avg <- Reduce("+", lapply(coefs_list, function(m) m[, varnames, drop = FALSE])) / n_subj

  list(
    temporal = temporal_counts,
    temporal_avg = temporal_avg,
    contemporaneous = contemp_counts,
    contemporaneous_avg = contemp_avg,
    coefs = coefs_list,
    psi = psi_list,
    fit = fit_df,
    path_counts = path_counts,
    paths = syntax_list,
    group_paths = group_paths,
    individual_paths = ind_paths_list,
    syntax = syntax_list
  )
}


# ============================================================================
# S3 Methods
# ============================================================================

#' Print Method for net_gimme
#'
#' @param x A \code{net_gimme} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' panel <- data.frame(
#'   id = rep(1:5, each = 20),
#'   t  = rep(seq_len(20), 5),
#'   A  = rnorm(100), B = rnorm(100), C = rnorm(100)
#' )
#' gm <- build_gimme(panel, vars = c("A","B","C"), id = "id", time = "t")
#' print(gm)
#' }
#'
#' @export
print.net_gimme <- function(x, ...) {
  cat("GIMME Network Analysis\n")
  cat(strrep("-", 30), "\n")
  cat("Subjects:  ", x$n_subjects, "\n")
  cat("Variables: ", length(x$labels), " (",
      paste(x$labels, collapse = ", "), ")\n")
  cat("AR paths:  ", ifelse(x$config$ar, "yes", "no"), "\n")
  cat("Hybrid:    ", ifelse(x$config$hybrid, "yes", "no"), "\n\n")

  cat("Group-level paths found:", length(x$group_paths), "\n")
  if (length(x$group_paths) > 0) {
    for (gp in x$group_paths) cat("  ", gp, "\n")
  }

  n_ind <- vapply(x$individual_paths, length, integer(1))
  cat("\nIndividual-level paths: ",
      sprintf("mean %.1f, range %d-%d\n", mean(n_ind), min(n_ind), max(n_ind)))

  cat("\nTemporal path counts (lagged):\n")
  print(x$temporal)

  cat("\nContemporaneous path counts:\n")
  print(x$contemporaneous)

  invisible(x)
}


#' Summary Method for net_gimme
#'
#' @param object A \code{net_gimme} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' panel <- data.frame(
#'   id = rep(1:5, each = 20),
#'   t  = rep(seq_len(20), 5),
#'   A  = rnorm(100), B = rnorm(100), C = rnorm(100)
#' )
#' gm <- build_gimme(panel, vars = c("A","B","C"), id = "id", time = "t")
#' summary(gm)
#' }
#'
#' @export
summary.net_gimme <- function(object, ...) {
  cat("GIMME Network Analysis -- Summary\n")
  cat(strrep("=", 40), "\n\n")

  # Fit summary
  cat("FIT INDICES (per subject)\n")
  cat(strrep("-", 30), "\n")
  fit_nums <- object$fit[, c("rmsea", "srmr", "nnfi", "cfi")]
  fit_summary <- data.frame(
    metric = c("RMSEA", "SRMR", "NNFI", "CFI"),
    mean = vapply(fit_nums, mean, numeric(1), na.rm = TRUE),
    sd = vapply(fit_nums, stats::sd, numeric(1), na.rm = TRUE),
    min = vapply(fit_nums, min, numeric(1), na.rm = TRUE),
    max = vapply(fit_nums, max, numeric(1), na.rm = TRUE)
  )
  print(fit_summary, row.names = FALSE, digits = 3)

  # Group paths
  cat("\n\nGROUP-LEVEL PATHS (", length(object$group_paths), ")\n")
  cat(strrep("-", 30), "\n")
  if (length(object$group_paths) > 0) {
    for (gp in object$group_paths) { # nocov start
      cat("  ", gp, "\n") # nocov end
    }
  } else {
    cat("  (none)\n")
  }

  # Average temporal coefficients
  cat("\nAVERAGE TEMPORAL COEFFICIENTS\n")
  cat(strrep("-", 30), "\n")
  print(round(object$temporal_avg, 3))

  # Average contemporaneous coefficients
  cat("\nAVERAGE CONTEMPORANEOUS COEFFICIENTS\n")
  cat(strrep("-", 30), "\n")
  print(round(object$contemporaneous_avg, 3))

  fit_summary
}


#' Plot Method for net_gimme
#'
#' @param x A \code{net_gimme} object.
#' @param type Character. Plot type: \code{"temporal"}, \code{"contemporaneous"},
#'   \code{"individual"}, \code{"counts"}, or \code{"fit"}.
#' @param subject Integer or character. Subject to plot for \code{type = "individual"}.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' panel <- data.frame(
#'   id = rep(1:5, each = 20),
#'   t  = rep(seq_len(20), 5),
#'   A  = rnorm(100), B = rnorm(100), C = rnorm(100)
#' )
#' gm <- build_gimme(panel, vars = c("A","B","C"), id = "id", time = "t")
#' plot(gm, type = "temporal")
#' }
#'
#' @export
plot.net_gimme <- function(x, type = c("temporal", "contemporaneous",
                                          "individual", "counts", "fit"),
                             subject = NULL, ...) {
  type <- match.arg(type)

  if (type == "temporal") {
    mat <- x$temporal_avg
    # Threshold: only show paths present in > 50% of subjects
    threshold_mat <- x$temporal
    mat[threshold_mat < (x$n_subjects / 2)] <- 0
    .gimme_plot_matrix(mat, "Group Temporal Network")

  } else if (type == "contemporaneous") {
    mat <- x$contemporaneous_avg # nocov start
    threshold_mat <- x$contemporaneous
    mat[threshold_mat < (x$n_subjects / 2)] <- 0
    .gimme_plot_matrix(mat, "Group Contemporaneous Network") # nocov end

  } else if (type == "individual") {
    if (is.null(subject)) {
      subject <- 1L
      message("No subject specified, plotting subject 1.")
    }
    if (is.character(subject)) {
      subj_idx <- which(names(x$coefs) == subject)
    } else {
      subj_idx <- subject
    }
    if (length(subj_idx) == 0 || subj_idx > length(x$coefs)) {
      stop("Subject not found.", call. = FALSE)
    }

    subj_name <- names(x$coefs)[subj_idx]
    mat <- x$coefs[[subj_idx]]
    p <- length(x$labels)

    # Split into temporal and contemporaneous
    temp_mat <- mat[, seq_len(p), drop = FALSE]
    colnames(temp_mat) <- x$labels
    cont_mat <- mat[, (p + 1):(2 * p), drop = FALSE]

    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    .gimme_plot_matrix(temp_mat, paste0(subj_name, " -- Temporal"))
    .gimme_plot_matrix(cont_mat, paste0(subj_name, " -- Contemporaneous"))

  } else if (type == "counts") {
    old_par <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    .gimme_plot_matrix(x$temporal, "Temporal Path Counts")
    .gimme_plot_matrix(x$contemporaneous, "Contemporaneous Path Counts")

  } else if (type == "fit") {
    fit_data <- x$fit
    metrics <- c("rmsea", "srmr", "cfi", "nnfi")
    old_par <- graphics::par(mfrow = c(2, 2))
    on.exit(graphics::par(old_par), add = TRUE)
    for (m in metrics) {
      vals <- fit_data[[m]]
      if (all(is.na(vals))) next
      graphics::hist(vals, main = toupper(m), xlab = m, col = "steelblue",
                     border = "white", breaks = 10)
    }
  }

  invisible(x)
}


#' Simple matrix heatmap fallback when cograph is not available
#' @noRd
.gimme_plot_matrix <- function(mat, title = "") {
  graphics::image(t(mat[nrow(mat):1, , drop = FALSE]),
                  axes = FALSE, main = title,
                  col = grDevices::heat.colors(50))
  graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)),
                 labels = colnames(mat), las = 2, cex.axis = 0.7)
  graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)),
                 labels = rev(rownames(mat)), las = 2, cex.axis = 0.7)
}
