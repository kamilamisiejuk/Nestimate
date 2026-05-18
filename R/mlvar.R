# ---- Multilevel Vector Autoregression (mlVAR) ----

#' Build a Multilevel Vector Autoregression (mlVAR) network
#'
#' @description Estimates three networks from ESM/EMA panel data, matching
#'   `mlVAR::mlVAR()` with `estimator = "lmer"`, `temporal = "fixed"`,
#'   `contemporaneous = "fixed"` at machine precision: (1) a directed
#'   temporal network of fixed-effect lagged regression coefficients, (2)
#'   an undirected contemporaneous network of partial correlations among
#'   residuals, and (3) an undirected between-subjects network of partial
#'   correlations derived from the person-mean fixed effects.
#'
#' #' @details The algorithm follows mlVAR's lmer pipeline exactly:
#' \enumerate{
#'   \item Drop rows with NA in id/day/beep and optionally grand-mean
#'         standardize each variable.
#'   \item Expand the per-(id, day) beep grid and right-join original
#'         values, producing the augmented panel (`augData`).
#'   \item Add within-person lagged predictors (`L1_*`) and person-mean
#'         predictors (`PM_*`).
#'   \item For each outcome variable fit
#'         `lmer(y ~ within + between-except-own-PM + (1 | id))` with
#'         `REML = FALSE`. Collect the fixed-effect temporal matrix `B`,
#'         between-effect matrix `Gamma`, random-intercept SDs (`mu_SD`),
#'         and lmer residual SDs.
#'   \item Contemporaneous network:
#'         `cor2pcor(D %*% cov2cor(cor(resid)) %*% D)`.
#'   \item Between-subjects network:
#'         `cor2pcor(pseudoinverse(forcePositive(D (I - Gamma))))`.
#' }
#'
#' Validated to machine precision (max_diff < 1e-10) against
#' `mlVAR::mlVAR()` on 25 real ESM datasets from `openesm` and 20 simulated
#' configurations (seeds 201-220). See `tmp/mlvar_equivalence_real20.R` and
#' `tmp/mlvar_equivalence_20seeds.R`.
#'
#' @param data A `data.frame` containing the panel data.
#' @param vars Character vector of variable column names to model.
#' @param id Character string naming the person-ID column.
#' @param day Character string naming the day/session column, or `NULL`.
#'   When provided, lag pairs are only formed within the same day.
#' @param beep Character string naming the measurement-occasion column, or
#'   `NULL`. When `NULL`, row position within each (id, day) is used.
#' @param lag Integer. The lag order (default 1).
#' @param standardize Logical. If `TRUE`, each variable is grand-mean
#'   centered and divided by its pooled SD *before* augmentation. Default
#'   `FALSE`, matching `mlVAR::mlVAR(scale = FALSE)` — the only setting for
#'   which numerical equivalence has been validated.
#'
#' @return A dual-class `c("net_mlvar", "netobject_group")` object — a
#'   named list of three full netobjects, one per network, plus
#'   model-level metadata stored as attributes. Each element is a
#'   standard `c("netobject", "cograph_network")` weight-matrix wrapper
#'   (no raw `$data`), so `print()`, `summary()`, [coefs()], and
#'   `cograph::splot(fit$temporal)` work directly. The three constituents
#'   are matrix-wrapped and carry no underlying panel data, so
#'   data-resampling verbs such as [bootstrap_network()] (and
#'   reliability/stability) cannot iterate over them — extract a single
#'   constituent and rebuild via [build_network()] if you need those.
#'   Structure:
#'   \describe{
#'     \item{`fit$temporal`}{Directed netobject for the `d x d` matrix of
#'       fixed-effect lagged coefficients. `$weights[i, j]` is the effect
#'       of variable j at t-lag on variable i at t. `method =
#'       "mlvar_temporal"`, `directed = TRUE`.}
#'     \item{`fit$contemporaneous`}{Undirected netobject for the `d x d`
#'       partial-correlation network of within-person lmer residuals.
#'       `method = "mlvar_contemporaneous"`, `directed = FALSE`.}
#'     \item{`fit$between`}{Undirected netobject for the `d x d`
#'       partial-correlation network of person means, derived from
#'       `D (I - Gamma)`. `method = "mlvar_between"`, `directed = FALSE`.}
#'     \item{`attr(fit, "coefs")` / [coefs()]}{Tidy `data.frame` with one
#'       row per `(outcome, predictor)` pair and columns `outcome`,
#'       `predictor`, `beta`, `se`, `t`, `p`, `ci_lower`, `ci_upper`,
#'       `significant`. Filter, sort, or plot with base R or the tidyverse.
#'       Retrieve with `coefs(fit)`.}
#'     \item{`attr(fit, "n_obs")`}{Number of rows in the augmented panel
#'       after na.omit.}
#'     \item{`attr(fit, "n_subjects")`}{Number of unique subjects remaining.}
#'     \item{`attr(fit, "lag")`}{Lag order used.}
#'     \item{`attr(fit, "standardize")`}{Logical; whether pre-augmentation
#'       standardization was applied.}
#'   }
#'
#' @examples
#' \dontrun{
#' d <- simulate_data("mlvar", seed = 1)
#' fit <- build_mlvar(d, vars = attr(d, "vars"),
#'                    id = "id", day = "day", beep = "beep")
#' print(fit)
#' summary(fit)
#' }
#'
#' @seealso [build_network()]
#' @export
build_mlvar <- function(data, vars, id,
                        day = NULL, beep = NULL,
                        lag = 1L,
                        standardize = FALSE) {
  # ---- Input validation ----
  stopifnot(
    is.data.frame(data),
    is.character(vars), length(vars) >= 2L,
    is.character(id), length(id) == 1L,
    is.numeric(lag), length(lag) == 1L, lag >= 1L,
    is.logical(standardize), length(standardize) == 1L
  )

  required <- c(vars, id)
  if (!is.null(day))  required <- c(required, day)
  if (!is.null(beep)) required <- c(required, beep)
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (!requireNamespace("lme4", quietly = TRUE)) {
    stop("Package 'lme4' is required for build_mlvar().", call. = FALSE)
  }
  if (!requireNamespace("corpcor", quietly = TRUE)) {
    stop("Package 'corpcor' is required for build_mlvar().", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for build_mlvar().", call. = FALSE)
  }

  lag <- as.integer(lag)

  prepared <- .mlvar_prepare_data(data, vars, id, day, beep, standardize)
  aug      <- .mlvar_augment_data(prepared, vars, id, day, beep, lag)
  Res      <- .mlvar_estimate_lmer(aug$data, aug$predModel, vars, id)

  # Wrap each of the three matrices as a full cograph_network netobject via
  # the package-wide `.wrap_netobject()` constructor. Nestimate never calls
  # cograph — plotting is handled by cograph's existing splot.netobject /
  # splot.cograph_network dispatch, which fires automatically because each
  # constituent here is a standard netobject.
  temporal_net        <- .wrap_netobject(Res$temporal$B,
                                         method   = "mlvar_temporal",
                                         directed = TRUE)
  contemporaneous_net <- .wrap_netobject(Res$contemporaneous,
                                         method   = "mlvar_contemporaneous",
                                         directed = FALSE)
  between_net         <- .wrap_netobject(Res$between,
                                         method   = "mlvar_between",
                                         directed = FALSE)

  nets <- list(
    temporal        = temporal_net,
    contemporaneous = contemporaneous_net,
    between         = between_net
  )

  # Model-level metadata lives in attributes so the list stays a pure
  # netobject_group (each element is a netobject). Use coefs(fit) to
  # retrieve the tidy coefs data.frame.
  attr(nets, "coefs")       <- Res$temporal$coefs
  attr(nets, "n_obs")       <- nrow(aug$data)
  attr(nets, "n_subjects")  <- length(unique(aug$data[[id]]))
  attr(nets, "lag")         <- lag
  attr(nets, "standardize") <- standardize
  attr(nets, "group_col")   <- "network_type"

  class(nets) <- c("net_mlvar", "netobject_group")
  nets
}

#' Tidy coefficients from a fitted mlvar model
#'
#' Generic accessor for the tidy coefficient table stored on a
#' [build_mlvar()] result. Returns a `data.frame` with one row per
#' `(outcome, predictor)` pair and columns `outcome`, `predictor`,
#' `beta`, `se`, `t`, `p`, `ci_lower`, `ci_upper`, `significant`.
#'
#' Only the within-person (temporal) coefficients are tabulated —
#' these are the lagged fixed effects that populate `fit$temporal`.
#' The between-subjects effects that go into `fit$between` are handled
#' via the `D (I - Gamma)` transformation and are not exposed as a
#' separate tidy table.
#'
#' @param x A fitted model object — currently only `net_mlvar` is supported.
#' @param ... Unused.
#' @return A tidy `data.frame` of coefficient estimates.
#' @inherit build_mlvar examples
#' @export
coefs <- function(x, ...) {
  UseMethod("coefs")
}

#' @rdname coefs
#' @export
coefs.net_mlvar <- function(x, ...) {
  attr(x, "coefs")
}

#' @rdname coefs
#' @export
coefs.default <- function(x, ...) {
  stop("No coefs() method for object of class '",
       class(x)[1], "'", call. = FALSE)
}

# ---- Internal helpers --------------------------------------------------

#' Drop rows with NA metadata and optionally grand-mean standardize
#' @noRd
.mlvar_prepare_data <- function(data, vars, id, day, beep, scale) {
  df <- as.data.frame(data)

  md_cols <- c(id,
               if (!is.null(day))  day,
               if (!is.null(beep)) beep)
  df <- df[stats::complete.cases(df[, md_cols, drop = FALSE]), , drop = FALSE]

  if (isTRUE(scale)) {
    for (v in vars) {
      x <- as.numeric(df[[v]])
      sd_val <- stats::sd(x, na.rm = TRUE)
      if (is.na(sd_val) || sd_val == 0) {
        df[[v]] <- 0
      } else {
        df[[v]] <- (x - mean(x, na.rm = TRUE)) / sd_val
      }
    }
  }
  df
}

#' Beep-grid augmentation + within/between predictor construction
#'
#' Hybrid implementation: `data.table::CJ` + keyed join for the expensive
#' grid construction, base R `ave()` for the within-group lag/center/mean
#' arithmetic. The split matters because `data.table`'s optimized group
#' aggregators (`gmean` et al.) accumulate sums in a different order than
#' base R, which drifts at 1e-16 and amplifies through lmer into 1e-10
#' coefficient diffs against mlVAR.
#' @noRd
.mlvar_augment_data <- function(data, vars, id, day, beep, lag) {
  # Silence R CMD check NOTE for data.table NSE symbols
  .first <- .last <- NULL

  id_col   <- id
  day_col  <- if (is.null(day))  ".day"  else day
  beep_col <- if (is.null(beep)) ".beep" else beep

  dt <- data.table::as.data.table(data)
  if (is.null(day)) {
    dt[, (day_col) := 1L]
  }
  if (is.null(beep)) {
    dt[, (beep_col) := seq_len(.N), by = c(id_col, day_col)]
  }

  # Per-(id, day) beep range
  first_last <- dt[, .(.first = min(get(beep_col), na.rm = TRUE),
                       .last  = max(get(beep_col), na.rm = TRUE)),
                   by = c(id_col, day_col)]

  # Global beep grid, filtered down to each (id, day)'s actual range
  gb_min <- min(dt[[beep_col]], na.rm = TRUE)
  gb_max <- max(dt[[beep_col]], na.rm = TRUE)
  allBeeps <- data.table::CJ(
    V1 = unique(dt[[id_col]]),
    V2 = unique(dt[[day_col]]),
    V3 = seq(gb_min, gb_max),
    sorted = FALSE
  )
  data.table::setnames(allBeeps, c("V1", "V2", "V3"),
                       c(id_col, day_col, beep_col))
  allBeeps <- allBeeps[first_last, on = c(id_col, day_col), nomatch = NULL]
  allBeeps <- allBeeps[get(beep_col) >= .first & get(beep_col) <= .last]
  allBeeps[, c(".first", ".last") := NULL]

  # Right-join original data onto the grid (keeps all grid rows)
  data.table::setkeyv(dt, c(id_col, day_col, beep_col))
  data.table::setkeyv(allBeeps, c(id_col, day_col, beep_col))
  augData <- dt[allBeeps, on = c(id_col, day_col, beep_col),
                allow.cartesian = TRUE]
  data.table::setkeyv(augData, c(id_col, day_col, beep_col))

  # Drop back to a plain data.frame so subsequent `ave()` calls match
  # mlVAR's accumulation order bit-for-bit.
  augData <- as.data.frame(augData)
  rownames(augData) <- NULL

  predModel <- list()

  # Within (lagged, person-centered) predictors
  for (v in vars) {
    p_id <- paste0("L", lag, "_", v)
    augData[[p_id]] <- stats::ave(
      augData[[v]], augData[[id_col]], augData[[day_col]],
      FUN = function(x) .mlvar_aveLag(x, lag)
    )
    augData[[p_id]] <- stats::ave(
      augData[[p_id]], augData[[id_col]],
      FUN = function(x) x - mean(x, na.rm = TRUE)
    )
    predModel[[length(predModel) + 1L]] <- list(
      dep = vars, pred = v, id = p_id, type = "within"
    )
  }

  # Between (person-mean) predictors
  for (v in vars) {
    p_id <- paste0("PM_", v)
    augData[[p_id]] <- stats::ave(
      augData[[v]], augData[[id_col]],
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    predModel[[length(predModel) + 1L]] <- list(
      dep = vars, pred = v, id = p_id, type = "between"
    )
  }

  involved <- unique(c(vars, vapply(predModel, `[[`, character(1), "id")))
  augData <- stats::na.omit(
    augData[, c(involved, id_col, day_col, beep_col), drop = FALSE]
  )
  rownames(augData) <- NULL

  list(data = augData, predModel = predModel)
}

#' Fit d outcome-specific lmer models and assemble the three networks
#'
#' Matches `mlVAR:::lmer_mlVAR` with `temporal = "fixed"`,
#' `contemporaneous = "fixed"`. For each outcome k fits
#' `outcome_k ~ L1_v1 + ... + L1_vd + PM_v_{-k} + (1 | id)` with
#' `REML = FALSE`, then assembles Beta, Gamma, mu_SD, residuals.
#' @noRd
.mlvar_estimate_lmer <- function(augData, predModel, vars, id) {
  d <- length(vars)
  n_obs <- nrow(augData)

  B             <- matrix(0, d, d, dimnames = list(vars, vars))
  Gamma         <- matrix(0, d, d, dimnames = list(vars, vars))
  mu_SD         <- stats::setNames(numeric(d), vars)
  sigma_vec     <- stats::setNames(numeric(d), vars)
  residuals_mat <- matrix(NA_real_, n_obs, d, dimnames = list(NULL, vars))

  within_ids <- vapply(Filter(function(m) m$type == "within",  predModel),
                       `[[`, character(1), "id")
  between_ids <- vapply(Filter(function(m) m$type == "between", predModel),
                        `[[`, character(1), "id")
  var_to_within  <- stats::setNames(within_ids,  vars)
  var_to_between <- stats::setNames(between_ids, vars)

  z975 <- stats::qnorm(0.975)

  # Tidy coefs: one row per (outcome, predictor) pair — fills d * d rows.
  # Faster and cleaner than growing a list of per-outcome data.frames and
  # `do.call(rbind, ...)` at the end.
  n_coef_rows <- d * d
  coefs_tidy <- data.frame(
    outcome     = rep(vars, each = d),
    predictor   = rep(vars, times = d),
    beta        = numeric(n_coef_rows),
    se          = numeric(n_coef_rows),
    t           = numeric(n_coef_rows),
    p           = numeric(n_coef_rows),
    ci_lower    = numeric(n_coef_rows),
    ci_upper    = numeric(n_coef_rows),
    significant = logical(n_coef_rows),
    stringsAsFactors = FALSE
  )

  for (k in seq_len(d)) {
    outcome <- vars[k]
    # Own PM excluded — matches mlVAR's `getModel` filter on `dep == outcome`
    fixed_preds <- c(within_ids, var_to_between[-k])

    # Random intercept only — matches mlVAR temporal="fixed".
    fm_str <- paste0(outcome, " ~ ",
                     paste(fixed_preds, collapse = " + "),
                     " + (1 | ", id, ")")
    fit <- suppressMessages(suppressWarnings(
      lme4::lmer(stats::as.formula(fm_str), data = augData, REML = FALSE)
    ))

    # Convergence diagnostics — warn but don't stop (matches mlVAR behaviour)
    if (lme4::isSingular(fit)) {
      warning(sprintf(
        "Model for '%s': singular fit (random-effects variance near zero).",
        outcome
      ), call. = FALSE)
    }
    conv_msgs <- fit@optinfo$conv$lme4$messages
    if (length(conv_msgs) > 0L) {
      warning(sprintf(
        "Model for '%s': %s", outcome, paste(conv_msgs, collapse = "; ")
      ), call. = FALSE)
    }

    fe <- lme4::fixef(fit)
    B[k, ]        <- fe[var_to_within[vars]]
    Gamma[k, -k]  <- fe[var_to_between[vars[-k]]]

    vc <- lme4::VarCorr(fit)
    ri_var <- as.numeric(vc[[id]][1, 1])
    mu_SD[k] <- if (!is.na(ri_var) && ri_var > 0) sqrt(ri_var) else 0

    sigma_vec[k] <- stats::sigma(fit)

    # Align residuals to augData row order (lmer drops any NA rows)
    res <- stats::residuals(fit)
    row_names <- rownames(augData)
    if (!is.null(row_names) && !is.null(names(res))) {
      residuals_mat[, k] <- res[match(row_names, names(res))]
    } else {
      residuals_mat[, k] <- res
    }

    sfe    <- summary(fit)$coefficients
    beta_k <- as.numeric(fe[var_to_within[vars]])
    se_k   <- as.numeric(sfe[var_to_within[vars], "Std. Error"])
    t_k    <- as.numeric(sfe[var_to_within[vars], "t value"])
    p_k    <- 2 * (1 - stats::pnorm(abs(t_k)))

    rows <- ((k - 1L) * d + 1L):(k * d)
    coefs_tidy$beta[rows]        <- beta_k
    coefs_tidy$se[rows]          <- se_k
    coefs_tidy$t[rows]           <- t_k
    coefs_tidy$p[rows]           <- p_k
    coefs_tidy$ci_lower[rows]    <- beta_k - z975 * se_k
    coefs_tidy$ci_upper[rows]    <- beta_k + z975 * se_k
    coefs_tidy$significant[rows] <- p_k < 0.05
  }

  contemporaneous <- .mlvar_contemporaneous_fixed(residuals_mat, sigma_vec, vars)
  between         <- .mlvar_compute_between_from_gamma(Gamma, mu_SD, vars)

  list(temporal = list(B = B, coefs = coefs_tidy, residuals = residuals_mat),
       contemporaneous = contemporaneous,
       between = between)
}

#' Within-group lag (matches mlVAR:::aveLag)
#'
#' Uses a logical `NA` (not `NA_real_`) for the prepended entries so that
#' integer input columns retain integer type. Preserving the integer type
#' is critical because base R's `mean()` uses a two-pass summation
#' correction for numeric input but a simple sum/n for integer input — the
#' two paths drift by ~1.4e-14, which then amplifies through lmer into
#' ~1e-10 coefficient diffs against mlVAR's integer-typed pipeline.
#' @noRd
.mlvar_aveLag <- function(x, lag = 1L) {
  n <- length(x)
  if (lag >= n) return(rep(NA, n))
  c(rep(NA, lag), x[seq_len(n - lag)])
}

#' Force a symmetric matrix to be positive-definite — byte-for-byte replica
#' of `mlVAR:::forcePositive`.
#'
#' Note the scalar-recycling quirk in the upstream implementation. In
#' `x - (diag(n) * min_ev - 0.001)`, the `0.001` scalar is subtracted from
#' every element of the diagonal matrix — so the final operation adds
#' `|min_ev|` to the diagonal *and* `+0.001` to every off-diagonal element.
#' This looks unintentional upstream but has to be replicated exactly for
#' equivalence with `mlVAR::mlVAR()`.
#' @noRd
.mlvar_force_positive <- function(x) {
  x <- (x + t(x)) / 2
  ev <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev < 0)) {
    x - (diag(nrow(x)) * min(ev) - 0.001)
  } else {
    x
  }
}

#' Between-subjects partial correlation from Gamma + mu_SD
#'
#' Matches mlVAR's Omega_mu branch:
#'   `D = diag(1 / mu_SD^2)`
#'   `inv = forcePositive(D (I - Gamma))`
#'   `cov = corpcor::pseudoinverse(inv)`
#'   `pcor = corpcor::cor2pcor(cov)`
#' @noRd
.mlvar_compute_between_from_gamma <- function(Gamma, mu_SD, vars) {
  d <- length(vars)
  if (any(mu_SD == 0)) {
    return(matrix(0, d, d, dimnames = list(vars, vars)))
  }

  D <- diag(1 / mu_SD^2)
  inv <- D %*% (diag(d) - Gamma)
  inv <- (inv + t(inv)) / 2
  inv <- .mlvar_force_positive(inv)

  mu_cov <- corpcor::pseudoinverse(inv)
  pcor <- corpcor::cor2pcor(mu_cov)
  diag(pcor) <- 0
  rownames(pcor) <- colnames(pcor) <- vars
  pcor
}

#' Contemporaneous partial correlation via mlVAR's "fixed" path
#'
#' Replicates `mlVAR:::lmer_mlVAR` Theta assembly for
#' `contemporaneous = "fixed"`: rescale the residual correlation by the
#' per-outcome lmer residual SDs and take `cor2pcor` directly. No
#' EBIC-GLASSO regularization. Note `cor2pcor` is scale-invariant, so the
#' `D %*% . %*% D` rescaling does not affect the pcor output — it is kept
#' only for parity with mlVAR's `cov`/`prec` slots.
#' @noRd
.mlvar_contemporaneous_fixed <- function(residuals_mat, sigma_vec, vars) {
  d <- length(vars)
  R <- stats::cor(residuals_mat, use = "pairwise.complete.obs")
  if (any(is.na(R))) {
    return(matrix(0, d, d, dimnames = list(vars, vars)))
  }
  D <- diag(sigma_vec)
  Theta_cov <- D %*% stats::cov2cor(R) %*% D
  pcor <- corpcor::cor2pcor(Theta_cov)
  diag(pcor) <- 0
  rownames(pcor) <- colnames(pcor) <- vars
  pcor
}


# ---- S3 methods --------------------------------------------------------

#' Print method for net_mlvar
#'
#' @param x A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `x`.
#' @inherit build_mlvar examples
#' @export
print.net_mlvar <- function(x, ...) {
  coef_df <- attr(x, "coefs")
  d <- nrow(x$temporal$weights)
  n_sig <- sum(coef_df$significant, na.rm = TRUE)
  n_tot <- nrow(coef_df)

  cat(sprintf("mlVAR result: %d subjects, %d observations, %d variables (lag %d)\n",
              attr(x, "n_subjects"),
              attr(x, "n_obs"),
              d,
              attr(x, "lag")))
  cat(sprintf("  $temporal        %d x %d directed    (%d/%d edges significant at p<0.05)\n",
              d, d, n_sig, n_tot))
  cat(sprintf("  $contemporaneous %d x %d undirected\n", d, d))
  cat(sprintf("  $between         %d x %d undirected\n", d, d))
  cat("Access: fit$temporal, fit$contemporaneous, fit$between (each a netobject);",
      "coefs(fit) for tidy coefs.\n")
  invisible(x)
}

#' Summary method for net_mlvar
#'
#' @param object A `net_mlvar` object returned by [build_mlvar()].
#' @param ... Unused; present for S3 consistency.
#' @return Invisibly returns `object`.
#' @inherit build_mlvar examples
#' @export
summary.net_mlvar <- function(object, ...) {
  coef_df <- attr(object, "coefs")
  vars    <- rownames(object$temporal$weights)
  d       <- length(vars)

  cat("=== mlVAR Summary ===\n")
  cat("Subjects:",      attr(object, "n_subjects"),
      " | Observations:", attr(object, "n_obs"),
      " | Variables:",  d,
      " | Lag:",        attr(object, "lag"), "\n")
  cat("Standardized:",  attr(object, "standardize"), "\n\n")

  cat("--- Temporal Network (B matrix) ---\n")
  print(round(object$temporal$weights, 4))
  cat("\n")

  sig_rows <- coef_df[coef_df$significant, , drop = FALSE]
  if (nrow(sig_rows) > 0L) {
    cat("Significant temporal edges (p < 0.05):\n")
    sig_print <- sig_rows[, c("outcome", "predictor",
                              "beta", "se", "t", "p"), drop = FALSE]
    sig_print$beta <- round(sig_print$beta, 4)
    sig_print$se   <- round(sig_print$se, 4)
    sig_print$t    <- round(sig_print$t, 3)
    sig_print$p    <- round(sig_print$p, 4)
    rownames(sig_print) <- NULL
    print(sig_print)
  } else {
    cat("No significant temporal edges at p < 0.05.\n")
  }
  cat("\n")

  cat("--- Contemporaneous Network ---\n")
  print(round(object$contemporaneous$weights, 4))
  cat("\n")

  cat("--- Between-Subjects Network ---\n")
  print(round(object$between$weights, 4))
  cat("\n")

  coef_df
}

