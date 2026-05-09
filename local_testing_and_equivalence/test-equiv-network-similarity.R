# Equivalence: .network_similarity() worker
#
#   (A) refactor invariant — every metric matches its pre-refactor inline
#       formula on 100 random matrix pairs.
#   (B) cross-tool — every metric matches tna::compare()'s corresponding
#       summary-metric value on 30 random matrix pairs (using
#       scaling = "none" so we are comparing raw weight matrices, not
#       a scaling pipeline).
#
# Gated by NESTIMATE_EQUIV_TESTS=true. Tolerance 1e-12.

skip_unless_equiv <- function() {
  if (!identical(tolower(Sys.getenv("NESTIMATE_EQUIV_TESTS")), "true")) {
    testthat::skip("NESTIMATE_EQUIV_TESTS != true")
  }
  testthat::skip_if_not_installed("tna")
}

# Frozen pre-refactor formulas — these are the inline implementations as
# they existed in `.compare_impl()` before extraction. Treat as a contract:
# the worker must reproduce these values bit-for-bit on shared inputs.
.reference_inline_metrics <- function(W_x, W_y) {
  x <- W_x; y <- W_y
  n <- nrow(x)
  x_vec <- as.vector(x); y_vec <- as.vector(y)
  abs_diff <- abs(x_vec - y_vec)
  abs_x <- abs(x_vec); abs_y <- abs(y_vec)
  pos <- abs_x > 0 & abs_y > 0
  c(
    mean_abs_diff   = mean(abs_diff),
    median_abs_diff = stats::median(abs_diff),
    rms_diff        = sqrt(mean(abs_diff^2)),
    max_abs_diff    = max(abs_diff),
    rel_mean_abs    = mean(abs_diff) / mean(abs_y),
    cv_ratio        = stats::sd(x_vec) * mean(y_vec) /
                      (mean(x_vec) * stats::sd(y_vec)),
    pearson         = stats::cor(x_vec, y_vec, method = "pearson",
                                 use = "complete.obs"),
    spearman        = stats::cor(x_vec, y_vec, method = "spearman",
                                 use = "complete.obs"),
    kendall         = stats::cor(x_vec, y_vec, method = "kendall",
                                 use = "complete.obs"),
    euclidean       = sqrt(sum(abs_diff^2)),
    manhattan       = sum(abs_diff),
    canberra        = sum(abs_diff[pos] / (abs_x[pos] + abs_y[pos])),
    bray_curtis     = sum(abs_diff) / sum(abs_x + abs_y),
    frobenius       = sqrt(sum(abs_diff^2)) / sqrt(n / 2),
    cosine          = sum(x * y) /
                      (sqrt(sum(x^2)) * sqrt(sum(y^2))),
    jaccard         = sum(pmin(abs_x, abs_y)) / sum(pmax(abs_x, abs_y)),
    dice            = 2 * sum(pmin(abs_x, abs_y)) /
                      (sum(abs_x) + sum(abs_y)),
    overlap         = sum(pmin(abs_x, abs_y)) /
                      min(sum(abs_x), sum(abs_y)),
    rank_agreement  = mean(sign(diff(x)) == sign(diff(y))),
    sign_agreement  = mean(sign(x_vec) == sign(y_vec))
  )
}

testthat::test_that(".network_similarity matches pre-refactor inline formulas (100 reps)", {
  skip_unless_equiv()
  set.seed(2025051)

  metric_keys <- c("mean_abs_diff", "median_abs_diff", "rms_diff",
                   "max_abs_diff", "rel_mean_abs", "cv_ratio",
                   "pearson", "spearman", "kendall",
                   "euclidean", "manhattan", "canberra", "bray_curtis",
                   "frobenius", "cosine", "jaccard", "dice", "overlap",
                   "rank_agreement", "sign_agreement")

  for (rep in seq_len(100)) {
    n <- sample(4:12, 1)
    W_x <- matrix(stats::runif(n * n, 0, 1), n, n)
    W_y <- matrix(stats::runif(n * n, 0, 1), n, n)
    rownames(W_x) <- colnames(W_x) <- paste0("S", seq_len(n))
    rownames(W_y) <- colnames(W_y) <- paste0("S", seq_len(n))

    inline <- .reference_inline_metrics(W_x, W_y)
    worker <- Nestimate:::.network_similarity(W_x, W_y, metrics = metric_keys)
    delta  <- max(abs(inline - worker), na.rm = TRUE)
    testthat::expect_lt(delta, 1e-12,
                        label = paste("rep", rep, "inline vs worker"))
  }
})

testthat::test_that(".network_similarity matches tna::compare() metrics (30 reps)", {
  skip_unless_equiv()

  # Map worker keys -> tna's (category, metric) labels used in
  # tna::compare()$summary_metrics. RV requires matrix structure (covered).
  tna_label <- list(
    mean_abs_diff   = c("Weight Deviations",   "Mean Abs. Diff."),
    median_abs_diff = c("Weight Deviations",   "Median Abs. Diff."),
    rms_diff        = c("Weight Deviations",   "RMS Diff."),
    max_abs_diff    = c("Weight Deviations",   "Max Abs. Diff."),
    rel_mean_abs    = c("Weight Deviations",   "Rel. Mean Abs. Diff."),
    cv_ratio        = c("Weight Deviations",   "CV Ratio"),
    pearson         = c("Correlations",        "Pearson"),
    spearman        = c("Correlations",        "Spearman"),
    kendall         = c("Correlations",        "Kendall"),
    distance_cor    = c("Correlations",        "Distance"),
    euclidean       = c("Dissimilarities",     "Euclidean"),
    manhattan       = c("Dissimilarities",     "Manhattan"),
    canberra        = c("Dissimilarities",     "Canberra"),
    bray_curtis     = c("Dissimilarities",     "Bray-Curtis"),
    frobenius       = c("Dissimilarities",     "Frobenius"),
    cosine          = c("Similarities",        "Cosine"),
    jaccard         = c("Similarities",        "Jaccard"),
    dice            = c("Similarities",        "Dice"),
    overlap         = c("Similarities",        "Overlap"),
    rv              = c("Similarities",        "RV"),
    rank_agreement  = c("Pattern Similarities","Rank Agreement"),
    sign_agreement  = c("Pattern Similarities","Sign Agreement")
  )

  set.seed(2025052)
  for (rep in seq_len(30)) {
    n <- sample(4:12, 1)
    W_x <- matrix(stats::runif(n * n, 0, 1), n, n)
    W_y <- matrix(stats::runif(n * n, 0, 1), n, n)
    rownames(W_x) <- colnames(W_x) <- paste0("S", seq_len(n))
    rownames(W_y) <- colnames(W_y) <- paste0("S", seq_len(n))

    worker <- Nestimate:::.network_similarity(W_x, W_y)
    tna_out <- tna::compare(W_x, W_y, scaling = "none")$summary_metrics
    tna_out <- as.data.frame(tna_out)

    for (key in names(tna_label)) {
      lbl <- tna_label[[key]]
      tna_val <- tna_out$value[tna_out$category == lbl[1] &
                               tna_out$metric  == lbl[2]]
      our_val <- worker[[key]]
      delta <- abs(tna_val - our_val)
      testthat::expect_lt(
        delta, 1e-12,
        label = sprintf("rep %d %s (%s/%s) ours=%g tna=%g",
                        rep, key, lbl[1], lbl[2], our_val, tna_val)
      )
    }
  }
})
