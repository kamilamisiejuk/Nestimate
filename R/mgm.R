# ---- Mixed Graphical Models (MGM) ----

#' Mixed Graphical Model estimator (internal)
#'
#' Estimates a Mixed Graphical Model via nodewise L1-regularized regression
#' with EBIC lambda selection and LW thresholding. Implementation matches
#' \code{mgm::mgm()} with defaults
#' \code{lambdaSel = "EBIC"}, \code{ruleReg = "AND"}, \code{threshold = "LW"},
#' \code{overparameterize = FALSE}, \code{scale = TRUE} at machine precision.
#'
#' @param data A numeric matrix or data.frame with one column per variable.
#'   Continuous columns are coerced to numeric; categorical columns are
#'   coerced to factor.
#' @param type Character vector of length \code{ncol(data)}: \code{"g"} for
#'   Gaussian (continuous), \code{"c"} for categorical.
#' @param level Integer vector of length \code{ncol(data)}: number of levels
#'   for categorical columns, \code{1} for continuous. Retained for
#'   \code{mgm::mgm()} API parity and length-validated only; it does
#'   \emph{not} enter the model fit, the EBIC penalty, or thresholding.
#'   The EBIC neighbour count is derived from the non-zero coefficient
#'   columns (\code{colSums(beta != 0)}), so categorical dummies are
#'   penalised individually without needing an explicit level count.
#'   Changing \code{level} has no effect on the returned network.
#' @param lambdaGam EBIC gamma parameter. Default 0.25.
#' @param ruleReg Symmetrization rule. \code{"AND"} (default) or \code{"OR"}.
#' @param threshold Coefficient thresholding rule. \code{"LW"} (default) or
#'   \code{"none"}.
#' @param scale Logical. Standardize continuous columns before fitting.
#'   Default \code{TRUE}.
#' @return A list with the symmetric weighted adjacency matrix \code{wadj}
#'   and per-node fit metadata.
#' @noRd
.mgm_estimate <- function(data, type, level,
                          lambdaGam = 0.25,
                          ruleReg   = c("AND", "OR"),
                          threshold = c("LW", "none"),
                          scale     = TRUE) {
  ruleReg   <- match.arg(ruleReg)
  threshold <- match.arg(threshold)
  stopifnot(
    requireNamespace("glmnet", quietly = TRUE),
    is.character(type), length(type) == ncol(data),
    all(type %in% c("g", "c")),
    is.numeric(level), length(level) == ncol(data)
  )

  data <- as.data.frame(data)
  n <- nrow(data)
  p <- ncol(data)

  # Scale ALL continuous columns of the data once (matches mgm main line 128).
  # Categorical columns become factors.
  if (scale) {
    for (i in which(type == "g")) {
      data[[i]] <- as.numeric(scale(data[[i]]))
    }
  }
  for (i in which(type == "c")) data[[i]] <- as.factor(data[[i]])

  # ---- per-node fits ----
  fits <- lapply(seq_len(p), function(v) {
    form <- stats::as.formula(paste(colnames(data)[v], "~ ."))
    mm   <- stats::model.matrix(form, data = data)
    X    <- mm[, -1, drop = FALSE]
    col_to_pred <- attr(mm, "assign")[-1]   # 1-based over rhs predictors
    rhs_idx <- setdiff(seq_len(p), v)
    pred_var <- rhs_idx[col_to_pred]
    npar_v <- ncol(X)
    y <- as.numeric(data[[v]])

    if (type[v] == "c") {
      fit <- glmnet::glmnet(X, y, family = "multinomial", alpha = 1,
                            intercept = TRUE)
      beta_list <- lapply(fit$beta, as.matrix)
      n_lambda  <- length(fit$lambda)

      # n_neighbors: count predictor columns with any non-zero across classes
      nz_per_col <- Reduce("+", lapply(beta_list, function(B) (B != 0) * 1)) > 0
      n_neighbors <- colSums(nz_per_col)

      # Null gaussian-style LL constant (drops out of argmin); replicate exactly
      tab  <- tabulate(y, nbins = max(y))
      pj   <- tab / n
      LL_null <- n * sum(pj[pj > 0] * log(pj[pj > 0]))
      LL_sat  <- 0.5 * fit$nulldev + LL_null

      deviance_seq <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda    <- -0.5 * deviance_seq + LL_sat
      EBIC <- -2 * LL_lambda + n_neighbors * log(n) +
              2 * lambdaGam * n_neighbors * log(npar_v)
      idx <- which.min(EBIC)

      beta_sel <- lapply(beta_list, function(B) B[, idx])
      if (threshold == "LW") {
        for (cl in seq_along(beta_sel)) {
          b <- beta_sel[[cl]]
          tau <- sqrt(sum(b ^ 2)) * sqrt(log(npar_v) / n)
          b[abs(b) < tau] <- 0
          beta_sel[[cl]] <- b
        }
      }

      list(beta_sel = beta_sel, npar = npar_v, pred_var = pred_var,
           lambda = fit$lambda[idx], multinomial = TRUE)

    } else {
      fit <- glmnet::glmnet(X, y, family = "gaussian", alpha = 1,
                            intercept = TRUE)
      beta_path <- as.matrix(fit$beta)
      n_lambda  <- length(fit$lambda)

      n_neighbors <- colSums(beta_path != 0)

      LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y)) ^ 2)) + 1)
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      deviance_seq <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda    <- -0.5 * deviance_seq + LL_sat
      EBIC <- -2 * LL_lambda + n_neighbors * log(n) +
              2 * lambdaGam * n_neighbors * log(npar_v)
      idx <- which.min(EBIC)

      b <- beta_path[, idx]
      if (threshold == "LW") {
        tau <- sqrt(sum(b ^ 2)) * sqrt(log(npar_v) / n)
        b[abs(b) < tau] <- 0
      }

      list(beta_sel = b, npar = npar_v, pred_var = pred_var,
           lambda = fit$lambda[idx], multinomial = FALSE)
    }
  })

  # ---- pair extraction + symmetrization ----
  side_mag <- function(node_fit, j) {
    cols <- which(node_fit$pred_var == j)
    if (!length(cols)) return(0)
    if (node_fit$multinomial) {
      vals <- unlist(lapply(node_fit$beta_sel, function(b) b[cols]))
    } else {
      vals <- node_fit$beta_sel[cols]
    }
    mean(abs(vals))
  }

  wadj <- matrix(0, p, p, dimnames = list(colnames(data), colnames(data)))
  for (i in seq_len(p - 1)) {
    for (j in (i + 1):p) {
      mag_ij <- side_mag(fits[[i]], j)
      mag_ji <- side_mag(fits[[j]], i)
      m_par_seq <- c(mag_ij, mag_ji)
      if (ruleReg == "AND") {
        edge <- if (any(m_par_seq == 0)) 0 else mean(m_par_seq)
      } else {
        edge <- mean(m_par_seq)
      }
      wadj[i, j] <- wadj[j, i] <- edge
    }
  }

  list(wadj = wadj, fits = fits)
}


#' MGM estimator hook for the registry
#'
#' @noRd
.estimator_mgm <- function(data, type = NULL, level = NULL,
                            lambdaGam = 0.25, ruleReg = "AND",
                            threshold = "LW", scale = TRUE, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Method 'mgm' requires the 'glmnet' package.", call. = FALSE)
  }
  data <- as.data.frame(data)
  p <- ncol(data)

  # Auto-detect types if not provided: factor / character / integer with
  # few unique values -> categorical; otherwise gaussian
  if (is.null(type)) {
    type <- vapply(data, function(col) {
      if (is.factor(col) || is.character(col)) "c"
      else if (is.numeric(col) &&
               length(unique(col[!is.na(col)])) <= 10 &&
               all(col == round(col), na.rm = TRUE)) "c"
      else "g"
    }, character(1))
    # Detection mirrors mgm::mgm()'s default (kept exactly so the
    # machine-precision equivalence holds), but the surprising case -- a
    # numeric column auto-classified CATEGORICAL via the <=10-integer rule
    # (e.g. a Likert/count item) -- must never be silent: it changes the
    # node's model from Gaussian to multinomial. Announce it; the user
    # overrides with explicit `type=`/`level=`.
    num_cat <- vapply(data, is.numeric, logical(1)) & type == "c"
    if (any(num_cat)) {
      warning(
        "mgm: numeric column(s) auto-detected as CATEGORICAL ",
        "(integer with <=10 distinct values): ",
        paste(names(data)[num_cat], collapse = ", "),
        ". If these are continuous, pass explicit `type=`/`level=`.",
        call. = FALSE
      )
    }
  }
  if (is.null(level)) {
    level <- vapply(seq_along(type), function(i) {
      if (type[i] == "c") length(unique(data[[i]])) else 1L
    }, integer(1))
  }

  res <- .mgm_estimate(
    data       = data,
    type       = type,
    level      = level,
    lambdaGam  = lambdaGam,
    ruleReg    = ruleReg,
    threshold  = threshold,
    scale      = scale
  )

  list(
    matrix       = res$wadj,
    nodes        = colnames(data),
    directed     = FALSE,
    cleaned_data = data,
    type         = type,
    level        = level
  )
}


# ---- Moderated MGM ----

# Helpers: magnitude extraction from per-node fits
# Used for computing AND indicators after fitting.

#' Per-side pairwise magnitude for variable `target` in a node's fit
#' @noRd
.mgm_main_mag <- function(node_fit, target) {
  cn <- node_fit$col_names
  idx <- which(!grepl(":", cn, fixed = TRUE) &
                 grepl(paste0("V", target, "."), cn, fixed = TRUE))
  if (!length(idx)) return(0)
  if (node_fit$multinomial) {
    vals <- unlist(lapply(node_fit$beta_list, function(b) b[idx]))
  } else {
    vals <- node_fit$beta[idx]
  }
  mean(abs(vals))
}

#' Per-side interaction magnitude for (target, other) pair in a node's fit
#' @noRd
.mgm_int_mag <- function(node_fit, target, other) {
  cn <- node_fit$col_names
  has_int <- grepl(":", cn, fixed = TRUE)
  has_t <- grepl(paste0("V", target, "."), cn, fixed = TRUE) & has_int
  has_o <- grepl(paste0("V", other, "."), cn, fixed = TRUE) & has_int
  idx <- which(has_t & has_o)
  if (!length(idx)) return(0)
  if (node_fit$multinomial) {
    vals <- unlist(lapply(node_fit$beta_list, function(b) b[idx]))
  } else {
    vals <- node_fit$beta[idx]
  }
  mean(abs(vals))
}

#' Conditioning indicator for a single interaction coefficient name
#'
#' For categorical moderator: 1 if the dummy level matches mod_value, else 0.
#' For continuous moderator: mod_value itself.
#' @noRd
.mgm_cond_indicator <- function(coef_name, moderator, mod_value, type_mod) {
  parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
  mod_pat <- paste0("V", moderator, ".")
  mod_part <- parts[grepl(mod_pat, parts, fixed = TRUE)]
  if (!length(mod_part)) return(0)
  stripped <- sub("^V", "", mod_part[1])
  if (type_mod == "c") {
    dot_parts <- strsplit(stripped, ".", fixed = TRUE)[[1]]
    if (length(dot_parts) < 2 || dot_parts[2] == "") return(mod_value)
    lev <- as.numeric(dot_parts[2])
    if (is.na(lev)) return(0)
    return(if (lev == mod_value) 1 else 0)
  }
  mod_value
}


#' Moderated Mixed Graphical Model
#'
#' Fits a single Mixed Graphical Model in which a chosen variable acts as a
#' moderator of every pairwise edge.  Group differences in network structure
#' show up as non-zero interaction terms between the moderator and the edge
#' endpoints.  Inspired by Haslbeck's "Group differences via moderation"
#' approach (\url{https://jonashaslbeck.com/Groupdifferences-via-Moderation/}).
#'
#' @section Equivalence:
#' Matches \code{mgm::mgm(..., moderators = k)} + \code{mgm::condition()} at
#' machine precision.  The pipeline mirrors mgm's internal flow:
#' \enumerate{
#'   \item Per-node nodewise lasso with interaction terms (same design matrix
#'         as \code{mgm::ModelMatrix_standard} moderator branch)
#'   \item \code{Reg2Graph}-equivalent 3-side symmetrization: for each
#'         (i, j, moderator) triple, evidence from i-as-outcome, j-as-outcome,
#'         AND moderator-as-outcome regressions is aggregated
#'   \item \code{applyTauAND}-equivalent pre-filter: pairwise main effects and
#'         3-way interactions are zeroed unless ALL sides survived
#'   \item \code{condition_core}-equivalent absorption: moderator-level
#'         indicators multiply interaction coefficients, which are then
#'         absorbed into the corresponding pairwise main-effect slots
#' }
#'
#' @param data Data frame or matrix.
#' @param type Character vector of "g"/"c" per column.
#' @param level Integer vector of levels per column (1 for continuous).
#'   Retained for \code{mgm::mgm()} API parity and echoed back in the
#'   returned structure only; it does \emph{not} enter the fit, the EBIC
#'   penalty, or the conditioning. Changing it has no effect on the result.
#' @param moderator Integer column index of the moderator variable.
#' @param lambdaGam EBIC tuning parameter. Default 0.25.
#' @param ruleReg Symmetrization rule. Default \code{"AND"}.
#' @param threshold Coefficient thresholding rule. Default \code{"LW"}.
#' @param scale Logical. Standardize continuous columns. Default \code{TRUE}.
#' @return A list with the per-node fits and AND indicators. Use
#'   \code{condition_moderated()} to extract the effective network at a
#'   specific moderator value.
#' @noRd
.mgm_estimate_moderated <- function(data, type, level, moderator,
                                     lambdaGam = 0.25,
                                     ruleReg   = c("AND", "OR"),
                                     threshold = c("LW", "none"),
                                     scale     = TRUE) {
  ruleReg   <- match.arg(ruleReg)
  threshold <- match.arg(threshold)
  stopifnot(
    requireNamespace("glmnet", quietly = TRUE),
    is.numeric(moderator), length(moderator) == 1L,
    moderator >= 1, moderator <= ncol(data)
  )

  data <- as.data.frame(data)
  p <- ncol(data)
  if (scale) {
    for (i in which(type == "g")) data[[i]] <- as.numeric(scale(data[[i]]))
  }
  for (i in which(type == "c")) data[[i]] <- as.factor(data[[i]])
  colnames(data) <- paste0("V", seq_len(p), ".")
  mod_name <- colnames(data)[moderator]
  n <- nrow(data)

  # ---- per-node nodewise lasso with interaction terms ----
  fits <- lapply(seq_len(p), function(v) {
    is_mod <- (v == moderator)
    main_rhs <- paste(colnames(data)[-v], collapse = " + ")
    if (is_mod) {
      other <- colnames(data)[-v]
      pairs_rhs <- utils::combn(other, 2, simplify = FALSE)
      int_terms <- vapply(pairs_rhs,
                          function(pr) paste(pr, collapse = "*"),
                          character(1))
    } else {
      nonmod_preds <- setdiff(colnames(data)[-v], mod_name)
      int_terms <- paste0(nonmod_preds, "*", mod_name)
    }
    form <- stats::as.formula(
      paste(colnames(data)[v], "~", main_rhs, "+",
            paste(int_terms, collapse = " + "))
    )
    mm <- stats::model.matrix(form, data = data)
    X <- mm[, -1, drop = FALSE]

    # Scale interaction columns where ALL involved variables are Gaussian
    # (matches mgm main lines 216-226)
    if (scale && any(type == "g")) {
      cn_x <- colnames(X)
      l_split <- strsplit(cn_x, ":", fixed = TRUE)
      all_gauss <- vapply(l_split, function(parts) {
        nums <- as.numeric(sub("\\.[0-9]*$", "",  sub("^V", "", parts)))
        all(!is.na(nums) & type[nums] == "g")
      }, logical(1))
      if (any(all_gauss))
        X[, all_gauss] <- apply(X[, all_gauss, drop = FALSE], 2, scale)
    }

    npar <- ncol(X)
    y <- as.numeric(data[[v]])

    if (type[v] == "c") {
      fit <- glmnet::glmnet(X, y, family = "multinomial", alpha = 1,
                             intercept = TRUE)
      beta_list <- lapply(fit$beta, as.matrix)
      n_lambda  <- length(fit$lambda)
      nz_per_col <- Reduce("+", lapply(beta_list,
                                        function(B) (B != 0) * 1)) > 0
      n_neighbors <- colSums(nz_per_col)

      tab <- tabulate(y, nbins = max(y))
      pj  <- tab / n
      LL_null <- n * sum(pj[pj > 0] * log(pj[pj > 0]))
      LL_sat  <- 0.5 * fit$nulldev + LL_null
      deviance_seq <- (1 - fit$dev.ratio) * fit$nulldev
      LL_lambda    <- -0.5 * deviance_seq + LL_sat
      EBIC <- -2 * LL_lambda + n_neighbors * log(n) +
              2 * lambdaGam * n_neighbors * log(npar)
      idx <- which.min(EBIC)
      beta_sel <- lapply(beta_list, function(B) B[, idx])
      if (threshold == "LW") {
        for (cl in seq_along(beta_sel)) {
          bb <- beta_sel[[cl]]
          tau <- sqrt(2L) * sqrt(sum(bb ^ 2)) * sqrt(log(npar) / n)
          bb[abs(bb) < tau] <- 0
          beta_sel[[cl]] <- bb
        }
      }
      return(list(beta = NULL, beta_list = beta_sel,
                  col_names = colnames(X), multinomial = TRUE,
                  lambda = fit$lambda[idx]))
    }

    fit <- glmnet::glmnet(X, y, family = "gaussian", alpha = 1,
                          intercept = TRUE)
    beta_path <- as.matrix(fit$beta)
    n_neighbors <- colSums(beta_path != 0)
    LL_null <- -n / 2 * (log(2 * pi * mean((y - mean(y)) ^ 2)) + 1)
    LL_sat  <- 0.5 * fit$nulldev + LL_null
    deviance_seq <- (1 - fit$dev.ratio) * fit$nulldev
    LL_lambda    <- -0.5 * deviance_seq + LL_sat
    EBIC <- -2 * LL_lambda + n_neighbors * log(n) +
            2 * lambdaGam * n_neighbors * log(npar)
    idx <- which.min(EBIC)
    b <- beta_path[, idx]
    if (threshold == "LW") {
      tau <- sqrt(2L) * sqrt(sum(b ^ 2)) * sqrt(log(npar) / n)
      b[abs(b) < tau] <- 0
    }
    list(beta = b, col_names = colnames(X), multinomial = FALSE,
         lambda = fit$lambda[idx])
  })

  # ---- Compute AND indicators (mirrors mgm::Reg2Graph aggregation) ----
  # Pairwise alive: 2-side AND on main-effect magnitudes
  pw_alive <- matrix(FALSE, p, p)
  # 3-way interaction alive: 3-side AND on interaction magnitudes
  int_alive <- matrix(FALSE, p, p)

  if (ruleReg == "AND") {
    all_pairs <- utils::combn(p, 2)
    for (k in seq_len(ncol(all_pairs))) {
      i <- all_pairs[1, k]; j <- all_pairs[2, k]
      pw_alive[i, j] <- pw_alive[j, i] <-
        (.mgm_main_mag(fits[[i]], j) > 0 && .mgm_main_mag(fits[[j]], i) > 0)
    }
    non_mod <- setdiff(seq_len(p), moderator)
    if (length(non_mod) >= 2L) {
      nm_pairs <- utils::combn(non_mod, 2)
      for (k in seq_len(ncol(nm_pairs))) {
        i <- nm_pairs[1, k]; j <- nm_pairs[2, k]
        mag_i <- .mgm_int_mag(fits[[i]], j, moderator)
        mag_j <- .mgm_int_mag(fits[[j]], i, moderator)
        mag_m <- .mgm_int_mag(fits[[moderator]], i, j)
        int_alive[i, j] <- int_alive[j, i] <-
          (mag_i > 0 && mag_j > 0 && mag_m > 0)
      }
    }
  } else {
    pw_alive[] <- TRUE
    non_mod <- setdiff(seq_len(p), moderator)
    int_alive[non_mod, non_mod] <- TRUE
  }

  structure(
    list(
      fits      = fits,
      moderator = moderator,
      p         = p,
      type      = type,
      level     = level,
      pw_alive  = pw_alive,
      int_alive = int_alive,
      params    = list(lambdaGam = lambdaGam, ruleReg = ruleReg,
                       threshold = threshold, scale = scale)
    ),
    class = "mgm_moderated"
  )
}


#' Condition a moderated MGM at a specific moderator value
#'
#' Mirrors \code{mgm::condition()}: applies the AND-rule pre-filter
#' (\code{applyTauAND}), absorbs the moderator value into main-effect
#' coefficients (\code{condition_core}), then re-aggregates pairwise edges
#' (\code{Reg2Graph} with \code{thresholding = FALSE}).
#'
#' @param fit An \code{mgm_moderated} object from \code{.mgm_estimate_moderated}.
#' @param mod_value The moderator value to condition on (numeric, e.g. 0 or 1
#'   for a binary moderator).
#' @param ruleReg Symmetrization rule. Defaults to the rule used during fit.
#' @return A symmetric \code{p x p} numeric matrix of effective edge weights.
#' @noRd
condition_moderated <- function(fit, mod_value, ruleReg = NULL) {
  stopifnot(inherits(fit, "mgm_moderated"))
  if (is.null(ruleReg)) ruleReg <- fit$params$ruleReg
  fits      <- fit$fits
  p         <- fit$p
  moderator <- fit$moderator
  pw_alive  <- fit$pw_alive
  int_alive <- fit$int_alive
  type_mod  <- fit$type[moderator]

  # Per-(node, target) conditioned magnitude: mirrors applyTauAND +
  # condition_core + Reg2Graph coefficient extraction per side.
  side_mag <- matrix(0, p, p)

  for (v in seq_len(p)) {
    if (v == moderator) next
    nf <- fits[[v]]
    cn <- nf$col_names

    for (j in seq_len(p)) {
      if (j == v || j == moderator) next
      # Main-effect indices for V_j (no ":" in name)
      m_idx <- which(!grepl(":", cn, fixed = TRUE) &
                       grepl(paste0("V", j, "."), cn, fixed = TRUE))
      # Interaction indices for V_j x V_moderator
      has_int <- grepl(":", cn, fixed = TRUE)
      has_j   <- grepl(paste0("V", j, "."), cn, fixed = TRUE) & has_int
      has_m   <- grepl(paste0("V", moderator, "."), cn, fixed = TRUE) & has_int
      i_idx   <- which(has_j & has_m)

      # Conditioned per-dummy magnitude (handles categorical targets correctly)
      .cond_dummies <- function(beta_vec) {
        vapply(m_idx, function(d) {
          val <- beta_vec[d]
          if (!pw_alive[v, j]) val <- 0
          # Find interactions absorbing into this exact dummy
          # model.matrix may order terms either way (V5.:V4.1 or V4.1:V5.)
          if (length(i_idx) && int_alive[v, j]) {
            dn <- gsub(".", "\\.", cn[d], fixed = TRUE)
            assoc <- i_idx[grepl(paste0("^", dn, ":"), cn[i_idx]) |
                             grepl(paste0(":", dn, "$"), cn[i_idx])]
            for (ii in assoc) {
              ind <- .mgm_cond_indicator(cn[ii], moderator, mod_value, type_mod)
              val <- val + beta_vec[ii] * ind
            }
          }
          val
        }, numeric(1))
      }

      if (!nf$multinomial) {
        cvals <- .cond_dummies(nf$beta)
        side_mag[v, j] <- if (length(cvals)) mean(abs(cvals)) else 0
      } else {
        per_class <- vapply(nf$beta_list, function(b) {
          cvals <- .cond_dummies(b)
          if (length(cvals)) mean(abs(cvals)) else 0
        }, numeric(1))
        side_mag[v, j] <- mean(per_class)
      }
    }
  }

  # Final 2-side pairwise aggregation (Reg2Graph on conditioned coefficients)
  wadj <- matrix(0, p, p)
  pairs <- utils::combn(p, 2)
  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k]; j <- pairs[2, k]
    if (i == moderator || j == moderator) next
    m_par <- c(side_mag[i, j], side_mag[j, i])
    edge <- if (ruleReg == "AND") {
      if (any(m_par == 0)) 0 else mean(m_par)
    } else {
      mean(m_par)
    }
    wadj[i, j] <- wadj[j, i] <- edge
  }
  wadj
}
