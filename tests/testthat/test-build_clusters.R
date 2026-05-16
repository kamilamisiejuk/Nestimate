testthat::skip_on_cran()

# ==============================================================================
# Tests for build_clusters()
# ==============================================================================

# ---- Test data ----

make_test_data <- function(n = 50, k = 26, n_states = 4, seed = 42) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  mat <- matrix(sample(states, n * k, replace = TRUE), nrow = n, ncol = k)
  colnames(mat) <- paste0("T", seq_len(k))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

make_data_with_na <- function(n = 30, k = 20, n_states = 3, seed = 42) {
  set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  mat <- matrix(sample(states, n * k, replace = TRUE), nrow = n, ncol = k)
  # Add void markers to last few columns
  for (i in seq_len(n)) {
    trail_start <- sample(k - 5, 1) + 5
    if (trail_start <= k) mat[i, trail_start:k] <- "%"
  }
  # Add some mid-sequence NAs
  mat[sample(n * k, 10)] <- "*"
  colnames(mat) <- paste0("T", seq_len(k))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# ==============================================================================
# 1. Input validation
# ==============================================================================

test_that("build_clusters validates inputs", {
  df <- make_test_data(n = 20, k = 10)
  # audit_clustering #2: explicit messages instead of bare predicate dump.
  expect_error(build_clusters(df, k = 1), "must be at least 2")
  expect_error(build_clusters(df, k = 20), "must be <= n - 1")
  expect_error(build_clusters(df, k = 2, dissimilarity = "invalid"))
  expect_error(build_clusters(df, k = 2, method = "invalid"))
  expect_error(build_clusters(df, k = 2, dissimilarity = "lv", weighted = TRUE),
               "Weighting is only supported")
  expect_error(build_clusters("not a df", k = 2))
})

test_that("build_clusters rejects all-missing input early (audit #4)", {
  # Every cell is a configured NA symbol, so .encode_sequences() finds 0
  # observed states. Without the early check the failure would be indirect
  # (downstream pam/hclust crash on a degenerate distance matrix).
  df <- data.frame(
    V1 = c("*", "*", "*"),
    V2 = c("%", "%", "%"),
    V3 = c("*", "%", "*"),
    stringsAsFactors = FALSE
  )
  expect_error(build_clusters(df, k = 2),
               "No observed sequence states were found")
})

test_that("build_clusters error messages name the offending argument", {
  # audit_clustering #2: failures should tell the user what to fix.
  df <- make_test_data(n = 20, k = 10)
  expect_error(build_clusters(df, k = "two"),
               "'k' must be a single numeric value")
  expect_error(build_clusters(df, k = 2, na_syms = 42),
               "'na_syms' must be a character vector")
  expect_error(build_clusters(df, k = 2, weighted = "yes"),
               "'weighted' must be a single logical")
})

test_that("build_clusters rejects non-finite or silently truncated numeric arguments", {
  df <- make_test_data(n = 20, k = 10)

  expect_error(build_clusters(df, k = 2.9), "'k' must be a whole")
  expect_error(build_clusters(df, k = NA_real_), "'k' must be a whole")
  expect_error(build_clusters(df, k = 2, weighted = NA),
               "'weighted' must be TRUE or FALSE")
  expect_error(build_clusters(df, k = 2, lambda = -1),
               "'lambda' must be a finite non-negative")
  expect_error(build_clusters(df, k = 2, q = 2.5),
               "'q' must be a positive whole")
  expect_error(build_clusters(df, k = 2, q = 0),
               "'q' must be a positive whole")
  expect_error(build_clusters(df, k = 2, p = 0.5),
               "'p' must be a finite number between 0 and 0.25")
})

test_that("build_clusters rejects unsupported extra arguments", {
  df <- make_test_data(n = 20, k = 10)
  expect_error(build_clusters(df, k = 2, typo_arg = TRUE),
               "unsupported argument: typo_arg")
})

# ==============================================================================
# 2. Basic clustering works for all metrics
# ==============================================================================

test_that("build_clusters runs for all 9 metrics", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  for (metric in c("hamming", "osa", "lv", "dl", "lcs",
                    "qgram", "cosine", "jaccard", "jw")) {
    cl <- build_clusters(df, k = 3, dissimilarity = metric)
    expect_true(inherits(cl, "net_clustering"), info = metric)
    expect_equal(cl$k, 3L, info = metric)
    expect_equal(sum(cl$sizes), 30L, info = metric)
    expect_length(cl$assignments, 30L)
    expect_true(inherits(cl$distance, "dist"), info = metric)
    expect_true(is.numeric(cl$silhouette), info = metric)
  }
})

# ==============================================================================
# 3. Clustering methods
# ==============================================================================

test_that("PAM returns medoids", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 3, method = "pam")
  expect_false(is.null(cl$medoids))
  expect_length(cl$medoids, 3L)
  expect_true(all(cl$medoids %in% seq_len(20)))
})

test_that("hierarchical methods work", {
  df <- make_test_data(n = 20, k = 10)
  for (m in c("ward.D2", "ward.D", "complete", "average", "single")) {
    cl <- build_clusters(df, k = 3, method = m)
    expect_true(inherits(cl, "net_clustering"), info = m)
    expect_true(is.null(cl$medoids), info = m)
    expect_equal(sum(cl$sizes), 20L, info = m)
  }
})

# ==============================================================================
# 4. Weighted Hamming
# ==============================================================================

test_that("weighted Hamming produces different distances than unweighted", {
  df <- make_test_data(n = 20, k = 10)
  cl_uw <- build_clusters(df, k = 2, dissimilarity = "hamming", weighted = FALSE)
  cl_w <- build_clusters(df, k = 2, dissimilarity = "hamming",
                       weighted = TRUE, lambda = 1)
  # Distance matrices should differ
  d_uw <- as.matrix(cl_uw$distance)
  d_w <- as.matrix(cl_w$distance)
  expect_false(isTRUE(all.equal(d_uw, d_w)))
  # Weighted distances should be <= unweighted (weights <= 1)
  expect_true(all(d_w <= d_uw + 1e-10))
})

test_that("lambda = 0 weighted matches unweighted", {
  df <- make_test_data(n = 20, k = 10)
  cl_uw <- build_clusters(df, k = 2, dissimilarity = "hamming", weighted = FALSE)
  # weighted = TRUE with lambda = 0 should give same distances
  # (lambda forced to 0 when weighted = FALSE, but lambda = 0 means uniform)
  # Actually testing: explicit lambda = 0 with weighted = TRUE
  # In our impl: lambda = 0 → weights = exp(0 * ...) / max = rep(1)
  # But the function sets lambda <- if (weighted) lambda else 0
  # With weighted = TRUE, lambda = 0 → weights = 1
  cl_w0 <- build_clusters(df, k = 2, dissimilarity = "hamming",
                        weighted = TRUE, lambda = 0)
  d_uw <- as.matrix(cl_uw$distance)
  d_w0 <- as.matrix(cl_w0$distance)
  expect_equal(d_uw, d_w0, tolerance = 1e-10)
})

# ==============================================================================
# 5. NA / void marker handling
# ==============================================================================

test_that("na_syms are treated as missing", {
  df <- make_data_with_na(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  expect_s3_class(cl, "net_clustering")
  expect_equal(sum(cl$sizes), 20L)
})

test_that("custom na_syms work", {
  df <- make_test_data(n = 20, k = 10, n_states = 3)
  # Replace some values with custom NA symbol
  df[1:5, 8:10] <- "MISSING"
  cl <- build_clusters(df, k = 2, na_syms = c("*", "%", "MISSING"))
  expect_s3_class(cl, "net_clustering")
})

# ==============================================================================
# 6. Seed reproducibility
# ==============================================================================

test_that("seed produces reproducible results", {
  df <- make_test_data(n = 30, k = 10)
  cl1 <- build_clusters(df, k = 3, seed = 123)
  cl2 <- build_clusters(df, k = 3, seed = 123)
  expect_equal(cl1$assignments, cl2$assignments)
  expect_equal(as.matrix(cl1$distance), as.matrix(cl2$distance))
})

# ==============================================================================
# 7. S3 methods
# ==============================================================================

test_that("print.net_clustering works", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  out <- capture.output(print(cl))
  expect_true(any(grepl("Sequence Clustering", out)))
  expect_true(any(grepl("pam", out)))
  expect_true(any(grepl("hamming", out)))
  expect_true(any(grepl("silhouette", out, ignore.case = TRUE)))
})

test_that("summary.net_clustering works", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  stats <- capture.output(res <- summary(cl))
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 2L)
  expect_true("size" %in% names(res))
  expect_true("mean_within_dist" %in% names(res))
})

test_that("summary.net_clustering matches direct within-distance computation", {
  df <- make_test_data(n = 24, k = 8, n_states = 3, seed = 31)
  cl <- build_clusters(df, k = 3, method = "ward.D2")
  capture.output(res <- summary(cl))

  expect_equal(res$size, as.integer(tabulate(cl$assignments, nbins = cl$k)))
  expect_equal(
    res$mean_within_dist,
    .per_cluster_within_dist(cl$distance, cl$assignments, cl$k),
    tolerance = 1e-12
  )
})

test_that("net_clustering methods reject unsupported dots and bad digits", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)

  expect_error(print(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(summary(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(plot(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(print(cl, digits = NA_real_),
               "'digits' must be a single non-negative whole")
  expect_error(print(cl, digits = 2.5),
               "'digits' must be a single non-negative whole")
})

test_that("plot.net_clustering silhouette works", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  p <- plot(cl, type = "silhouette")
  expect_s3_class(p, "ggplot")
})

test_that("plot.net_clustering mds works", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  p <- plot(cl, type = "mds")
  expect_s3_class(p, "ggplot")
})

test_that("plot.net_clustering heatmap works", {
  df <- make_test_data(n = 20, k = 10)
  cl <- build_clusters(df, k = 2)
  p <- plot(cl, type = "heatmap")
  expect_s3_class(p, "ggplot")
})

# ==============================================================================
# 9. Cross-validation against tna
# ==============================================================================

test_that("R fallback matches stringdist for all applicable metrics", {
  skip_if_not_installed("stringdist")
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  enc <- Nestimate:::.encode_sequences(df, c("*", "%"))

  for (metric in c("hamming", "osa", "lv", "dl", "lcs",
                    "qgram", "cosine", "jaccard", "jw")) {
    d_r <- as.matrix(
      Nestimate:::.dissimilarity_matrix_r(enc, metric, lambda = 0, q = 2L, p = 0.1)
    )
    d_sd <- as.matrix(
      Nestimate:::.dissimilarity_matrix_stringdist(enc, metric, lambda = 0, q = 2L, p = 0.1)
    )
    expect_equal(d_r, d_sd, tolerance = 1e-10, info = metric)
  }
})

# ==============================================================================
# 11. Edge cases
# ==============================================================================

test_that("single state data works", {
  df <- data.frame(T1 = rep("A", 10), T2 = rep("A", 10), T3 = rep("A", 10))
  cl <- build_clusters(df, k = 2, dissimilarity = "hamming")
  # All distances should be 0
  expect_true(all(as.matrix(cl$distance) == 0))
})

test_that("two-state binary data works", {
  set.seed(1)
  df <- data.frame(
    T1 = sample(c("X", "Y"), 20, replace = TRUE),
    T2 = sample(c("X", "Y"), 20, replace = TRUE),
    T3 = sample(c("X", "Y"), 20, replace = TRUE)
  )
  cl <- build_clusters(df, k = 2, dissimilarity = "hamming")
  expect_s3_class(cl, "net_clustering")
})

test_that("data with all trailing NAs handled", {
  df <- data.frame(
    T1 = c("A", "B", "A", "B", "C"),
    T2 = c("B", "A", "C", "A", "A"),
    T3 = c("%", "%", "%", "%", "%"),
    T4 = c("%", "%", "%", "%", "%")
  )
  cl <- build_clusters(df, k = 2, dissimilarity = "hamming")
  expect_s3_class(cl, "net_clustering")
  expect_equal(sum(cl$sizes), 5L)
})

test_that("matrix input works", {
  mat <- matrix(sample(LETTERS[1:3], 30, replace = TRUE), nrow = 10)
  colnames(mat) <- paste0("T", 1:3)
  cl <- build_clusters(mat, k = 2)
  expect_s3_class(cl, "net_clustering")
  expect_equal(sum(cl$sizes), 10L)
})

test_that("q-gram parameter affects result", {
  df <- make_test_data(n = 20, k = 10, n_states = 4)
  cl_q2 <- build_clusters(df, k = 2, dissimilarity = "qgram", q = 2L)
  cl_q3 <- build_clusters(df, k = 2, dissimilarity = "qgram", q = 3L)
  # Different q should give different distance matrices
  expect_false(isTRUE(all.equal(
    as.matrix(cl_q2$distance), as.matrix(cl_q3$distance)
  )))
})

test_that("p parameter affects jw result", {
  df <- make_test_data(n = 20, k = 10, n_states = 4)
  cl_p1 <- build_clusters(df, k = 2, dissimilarity = "jw", p = 0.1)
  cl_p2 <- build_clusters(df, k = 2, dissimilarity = "jw", p = 0.2)
  # Different p should give different distance matrices
  expect_false(isTRUE(all.equal(
    as.matrix(cl_p1$distance), as.matrix(cl_p2$distance)
  )))
})

# ==============================================================================
# 12. Encoding tests
# ==============================================================================

test_that(".encode_sequences handles na_syms correctly", {
  df <- data.frame(T1 = c("A", "*", "B"), T2 = c("B", "A", "%"))
  enc <- Nestimate:::.encode_sequences(df, na_syms = c("*", "%"))
  # df row 1 = c("A", "B"), row 2 = c("*", "A"), row 3 = c("B", "%")
  expect_true(is.na(enc$int_mat[2, 1]))  # "*" → NA
  expect_true(is.na(enc$int_mat[3, 2]))  # "%" → NA
  expect_equal(enc$states, c("A", "B"))
  expect_equal(enc$n_states, 2L)
  expect_equal(enc$len, c(2L, 2L, 1L))  # row 3: last obs at col 1
})

test_that(".encode_sequences integer encoding is correct", {
  df <- data.frame(T1 = c("C", "A", "B"), T2 = c("A", "B", "C"))
  enc <- Nestimate:::.encode_sequences(df, na_syms = c("*", "%"))
  # States sorted: A=1, B=2, C=3
  expect_equal(enc$states, c("A", "B", "C"))
  expect_equal(enc$int_mat[1, ], c(3L, 1L))  # C=3, A=1
  expect_equal(enc$int_mat[2, ], c(1L, 2L))  # A=1, B=2
  expect_equal(enc$int_mat[3, ], c(2L, 3L))  # B=2, C=3
})

# ==============================================================================
# 13. Individual distance function tests
# ==============================================================================

test_that("hamming distance is correct", {
  # Same sequence → 0
  expect_equal(Nestimate:::.hamming_dist_r(c(1,2,3), c(1,2,3), 3, 3, 1), 0)
  # All different → length
  expect_equal(Nestimate:::.hamming_dist_r(c(1,2,3), c(4,5,6), 3, 3, 1), 3)
  # One difference
  expect_equal(Nestimate:::.hamming_dist_r(c(1,2,3), c(1,2,4), 3, 3, 1), 1)
})

test_that("hamming distance with weights is correct", {
  w <- exp(-1 * c(0, 1, 2))
  w <- w / max(w)
  # All different: sum of weights
  expect_equal(
    Nestimate:::.hamming_dist_r(c(1,2,3), c(4,5,6), 3, 3, w),
    sum(w)
  )
  # First position different only
  expect_equal(
    Nestimate:::.hamming_dist_r(c(1,2,3), c(4,2,3), 3, 3, w),
    w[1]
  )
})

test_that("levenshtein distance is correct", {
  # kitten → sitting = 3
  x <- c(1, 2, 3, 3, 4, 5)  # k i t t e n
  y <- c(6, 2, 3, 3, 2, 5, 7)  # s i t t i n g
  expect_equal(Nestimate:::.levenshtein_dist_r(x, y, 6, 7), 3)
  # Same → 0
  expect_equal(Nestimate:::.levenshtein_dist_r(c(1,2), c(1,2), 2, 2), 0)
  # Empty vs non-empty
  expect_equal(Nestimate:::.levenshtein_dist_r(c(1,2,3), c(1,2,3), 0, 3), 3)
})

test_that("lcs distance is correct", {
  # lcs of (1,2,3,4) and (2,4,3) → lcs length = 2 (2,3 or 2,4)
  # distance = 4 + 3 - 2*2 = 3
  expect_equal(Nestimate:::.lcs_dist_r(c(1,2,3,4), c(2,4,3), 4, 3), 3)
  # Same → 0
  expect_equal(Nestimate:::.lcs_dist_r(c(1,2,3), c(1,2,3), 3, 3), 0)
})

# ==============================================================================
# 14. Input extraction tests
# ==============================================================================

test_that("build_clusters works on netobject from sequence methods", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  net <- build_network(df, method = "relative")
  cl_net <- build_clusters(net, k = 3)
  cl_df <- build_clusters(df, k = 3)
  expect_s3_class(cl_net, "net_clustering")
  expect_equal(as.matrix(cl_net$distance), as.matrix(cl_df$distance),
               tolerance = 1e-10)
})

test_that("build_clusters works on frequency netobject", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  net <- build_network(df, method = "frequency")
  cl <- build_clusters(net, k = 3)
  expect_s3_class(cl, "net_clustering")
  expect_equal(sum(cl$sizes), 30L)
})

test_that("build_clusters rejects association-method netobjects", {
  set.seed(42)
  ndf <- data.frame(matrix(rnorm(100), 20, 5))
  colnames(ndf) <- paste0("V", 1:5)
  net <- build_network(ndf, method = "cor")
  expect_error(build_clusters(net, k = 2), "sequence data")
})

test_that("build_network dispatches on net_clustering", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  cl <- build_clusters(df, k = 3)
  grp <- build_network(cl)
  expect_s3_class(grp, "netobject_group")
  expect_length(grp, 3L)
  expect_equal(names(grp), paste("Cluster", 1:3))
  # Each sub-network should be a netobject
  for (nm in names(grp)) {
    expect_s3_class(grp[[nm]], "netobject")
  }
})

test_that("build_network(clustering) default method is relative", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  cl <- build_clusters(df, k = 3)
  grp <- build_network(cl)
  expect_equal(grp[[1]]$method, "relative")
})

test_that("build_network(clustering) respects method arg", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  cl <- build_clusters(df, k = 3)
  grp <- build_network(cl, method = "frequency")
  expect_equal(grp[[1]]$method, "frequency")
})

test_that("build_network(clustering) sub-networks have correct sizes", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  cl <- build_clusters(df, k = 3)
  grp <- build_network(cl)
  total_seqs <- sum(vapply(grp, function(net) nrow(net$data), integer(1L)))
  expect_equal(total_seqs, 30L)
  # Sizes match clustering sizes
  for (i in seq_along(grp)) {
    expect_equal(nrow(grp[[i]]$data), unname(cl$sizes[i]))
  }
})

test_that("build_network(clustering) stores clustering metadata", {
  df <- make_test_data(n = 30, k = 10, n_states = 4)
  cl <- build_clusters(df, k = 3)
  grp <- build_network(cl)
  expect_s3_class(attr(grp, "clustering"), "net_clustering")
  expect_equal(attr(grp, "clustering")$k, 3L)
})

# ==============================================================================
# 16. Covariate analysis — input forms
# ==============================================================================

make_cov_data <- function(n = 40, seed = 42) {
  set.seed(seed)
  data.frame(
    T1 = sample(LETTERS[1:3], n, replace = TRUE),
    T2 = sample(LETTERS[1:3], n, replace = TRUE),
    T3 = sample(LETTERS[1:3], n, replace = TRUE),
    Age = rnorm(n, 25, 5),
    Gender = sample(c("M", "F"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

test_that("all 5 covariate input forms produce identical results", {
  df <- make_cov_data()
  cl_formula <- build_clusters(df, k = 3, covariates = ~ Age + Gender)
  cl_char <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  cl_string <- build_clusters(df, k = 3, covariates = "Age + Gender")
  cl_df <- build_clusters(df, k = 3, covariates = df[, c("Age", "Gender")])

  ref <- cl_char$covariates$coefficients
  expect_equal(cl_formula$covariates$coefficients, ref)
  expect_equal(cl_string$covariates$coefficients, ref)
  expect_equal(cl_df$covariates$coefficients, ref)
})

test_that("covariates = NULL produces no covariate analysis", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = NULL)
  expect_null(cl$covariates)
})

# ==============================================================================
# 17. Covariate analysis — profiles
# ==============================================================================

test_that("numeric profiles are correct", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  np <- cl$covariates$profiles$numeric
  expect_true(is.data.frame(np))
  expect_true(all(c("cluster", "n", "pct", "variable", "mean", "sd", "median")
                   %in% names(np)))
  # Verify against manual computation
  for (clust in unique(np$cluster)) {
    mask <- cl$assignments == clust
    manual_mean <- mean(df$Age[mask])
    manual_sd <- sd(df$Age[mask])
    manual_med <- median(df$Age[mask])
    row <- np[np$cluster == clust & np$variable == "Age", ]
    expect_equal(row$mean, manual_mean, tolerance = 1e-10)
    expect_equal(row$sd, manual_sd, tolerance = 1e-10)
    expect_equal(row$median, manual_med, tolerance = 1e-10)
  }
})

test_that("categorical profiles are correct", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  cp <- cl$covariates$profiles$categorical
  expect_true(is.data.frame(cp))
  expect_true(all(c("cluster", "n", "variable", "level", "count", "pct")
                   %in% names(cp)))
  # Verify against manual computation
  for (clust in unique(cp$cluster)) {
    mask <- cl$assignments == clust
    tab <- table(df$Gender[mask])
    for (lev in names(tab)) {
      row <- cp[cp$cluster == clust & cp$variable == "Gender" &
                cp$level == lev, ]
      expect_equal(row$count, as.integer(tab[[lev]]))
    }
  }
})

test_that("all-numeric covariates produce no categorical profile", {
  df <- make_cov_data()
  df$Score <- rnorm(nrow(df))
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Score"))
  expect_null(cl$covariates$profiles$categorical)
  expect_true(is.data.frame(cl$covariates$profiles$numeric))
})

test_that("all-categorical covariates produce no numeric profile", {
  df <- make_cov_data()
  df$Group <- sample(c("X", "Y"), nrow(df), replace = TRUE)
  cl <- build_clusters(df, k = 3, covariates = c("Gender", "Group"))
  expect_null(cl$covariates$profiles$numeric)
  expect_true(is.data.frame(cl$covariates$profiles$categorical))
})

# ==============================================================================
# 18. Covariate analysis — k=2 cross-validation against glm(binomial)
# ==============================================================================

test_that("k=2 multinomial matches glm(binomial)", {
  set.seed(99)
  df <- data.frame(
    T1 = sample(LETTERS[1:4], 60, replace = TRUE),
    T2 = sample(LETTERS[1:4], 60, replace = TRUE),
    T3 = sample(LETTERS[1:4], 60, replace = TRUE),
    T4 = sample(LETTERS[1:4], 60, replace = TRUE),
    Age = rnorm(60, 30, 8),
    Score = runif(60, 0, 100),
    stringsAsFactors = FALSE
  )
  cl <- build_clusters(df, k = 2, covariates = c("Age", "Score"))
  our <- cl$covariates$coefficients

  # glm reference
  cov_df <- df[, c("Age", "Score")]
  cov_df$cluster <- factor(cl$assignments)
  ref_fit <- glm(cluster ~ Age + Score, data = cov_df, family = binomial)
  ref <- summary(ref_fit)$coefficients

  # Compare coefficients (excluding intercept from our table for variable match)
  for (vname in c("Age", "Score")) {
    our_row <- our[our$variable == vname, ]
    expect_equal(our_row$estimate[[1L]], unname(ref[vname, "Estimate"]),
                 tolerance = 1e-3, info = vname)
    expect_equal(our_row$std_error[[1L]], unname(ref[vname, "Std. Error"]),
                 tolerance = 1e-3, info = vname)
    expect_equal(our_row$z[[1L]], unname(ref[vname, "z value"]),
                 tolerance = 1e-3, info = vname)
    expect_equal(our_row$p[[1L]], unname(ref[vname, "Pr(>|z|)"]),
                 tolerance = 1e-2, info = vname)
  }
})

# ==============================================================================
# 19. Covariate analysis — k>2 manual verification
# ==============================================================================

test_that("k>2 coefficients match manual extraction from multinom", {
  set.seed(77)
  df <- data.frame(
    T1 = sample(LETTERS[1:3], 50, replace = TRUE),
    T2 = sample(LETTERS[1:3], 50, replace = TRUE),
    T3 = sample(LETTERS[1:3], 50, replace = TRUE),
    X1 = rnorm(50),
    X2 = rnorm(50),
    stringsAsFactors = FALSE
  )
  cl <- build_clusters(df, k = 3, covariates = c("X1", "X2"))
  our <- cl$covariates$coefficients

  # Manual reference from raw multinom
  fit_df <- data.frame(cluster = factor(cl$assignments), X1 = df$X1, X2 = df$X2)
  fit <- nnet::multinom(cluster ~ X1 + X2, data = fit_df, trace = FALSE)
  s <- summary(fit)
  coefs <- s$coefficients
  ses <- s$standard.errors

  for (i in seq_len(nrow(coefs))) {
    for (j in seq_len(ncol(coefs))) {
      est <- coefs[i, j]
      se <- ses[i, j]
      z_manual <- est / se
      p_manual <- 2 * (1 - pnorm(abs(z_manual)))
      or_manual <- exp(est)
      ci_lo_manual <- exp(est - 1.96 * se)
      ci_hi_manual <- exp(est + 1.96 * se)

      row <- our[our$cluster == rownames(coefs)[i] &
                 our$variable == colnames(coefs)[j], ]
      expect_equal(row$estimate, est, tolerance = 1e-10)
      expect_equal(row$std_error, se, tolerance = 1e-10)
      expect_equal(row$odds_ratio, or_manual, tolerance = 1e-10)
      expect_equal(row$ci_lower, ci_lo_manual, tolerance = 1e-10)
      expect_equal(row$ci_upper, ci_hi_manual, tolerance = 1e-10)
      expect_equal(row$z, z_manual, tolerance = 1e-10)
      expect_equal(row$p, p_manual, tolerance = 1e-10)
    }
  }
})

# ==============================================================================
# 20. Covariate analysis — model fit
# ==============================================================================

test_that("model fit stats are correct", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  fit <- cl$covariates$fit
  expect_true(is.numeric(fit$aic))
  expect_true(is.numeric(fit$bic))
  expect_true(is.numeric(fit$deviance))
  expect_true(fit$mcfadden_r2 >= 0 && fit$mcfadden_r2 <= 1)
  expect_equal(fit$reference_cluster, "1")
})

# ==============================================================================
# 21. Covariate analysis — S3 methods
# ==============================================================================

test_that("print shows covariates line", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  out <- capture.output(print(cl))
  expect_true(any(grepl("Covariates:", out)))
  expect_true(any(grepl("post-hoc", out)))
})

test_that("summary shows covariate analysis", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  out <- capture.output(res <- summary(cl))
  expect_true(any(grepl("Post-hoc Covariate Analysis", out)))
  expect_true(any(grepl("Predictors of Membership", out)))
  expect_true(any(grepl("McFadden", out)))
  expect_true(any(grepl("does not influence", out)))
  # Return value is the cluster-stats data.frame; covariate block is an attribute.
  expect_s3_class(res, "data.frame")
  expect_true(!is.null(attr(res, "covariates")))
})

test_that("summary without covariates returns data.frame (backwards compat)", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3)
  out <- capture.output(res <- summary(cl))
  expect_true(is.data.frame(res))
})

test_that("plot type='predictors' returns ggplot", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  p <- plot(cl, type = "predictors")
  expect_s3_class(p, "ggplot")
})

test_that("plot type='predictors' errors without covariates", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3)
  expect_error(plot(cl, type = "predictors"), "No covariate analysis")
})

# ==============================================================================
# 22. Covariate analysis — edge cases
# ==============================================================================

test_that("single covariate works", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = "Age")
  expect_true(!is.null(cl$covariates))
  expect_true(all(cl$covariates$coefficients$variable %in%
                  c("(Intercept)", "Age")))
})

test_that("covariate with NAs warns and works", {
  df <- make_cov_data()
  df$Age[c(1, 5, 10)] <- NA
  expect_warning(
    cl <- build_clusters(df, k = 3, covariates = "Age"),
    "Dropped 3 rows"
  )
  expect_true(!is.null(cl$covariates))
})

test_that("constant covariate errors", {
  df <- make_cov_data()
  df$Const <- 5
  expect_error(
    build_clusters(df, k = 3, covariates = "Const"),
    "constant"
  )
})

test_that("covariate data.frame row count mismatch errors", {
  df <- make_cov_data()
  bad_cov <- data.frame(Age = rnorm(10))
  expect_error(
    build_clusters(df, k = 3, covariates = bad_cov),
    "rows"
  )
})

test_that("netobject metadata covariates work", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(LETTERS[1:3], 40, replace = TRUE),
    T2 = sample(LETTERS[1:3], 40, replace = TRUE),
    T3 = sample(LETTERS[1:3], 40, replace = TRUE),
    Age = rnorm(40, 25, 5),
    Score = runif(40, 0, 100),
    stringsAsFactors = FALSE
  )
  net <- build_network(df, method = "relative")
  # Age and Score are numeric → stored in $metadata
  expect_true(!is.null(net$metadata))
  cl <- build_clusters(net, k = 3, covariates = c("Age", "Score"))
  expect_true(!is.null(cl$covariates))
  # Should match direct data.frame input
  cl_df <- build_clusters(df, k = 3, covariates = c("Age", "Score"))
  expect_equal(cl$covariates$coefficients, cl_df$covariates$coefficients)
})

test_that("build_network round-trip preserves covariates", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  grp <- build_network(cl)
  stored <- attr(grp, "clustering")
  expect_true(!is.null(stored$covariates))
  expect_equal(stored$covariates$fit$aic, cl$covariates$fit$aic)
})

# ==============================================================================
# 23. dl distance function — zero-length edge cases (L87-88)
# ==============================================================================

test_that(".dl_dist_r zero-length edge cases work", {
  # n == 0 → return m
  expect_equal(Nestimate:::.dl_dist_r(integer(), c(1L, 2L), 0L, 2L), 2)
  # m == 0 → return n
  expect_equal(Nestimate:::.dl_dist_r(c(1L, 2L), integer(), 2L, 0L), 2)
  # Same sequence → 0
  expect_equal(Nestimate:::.dl_dist_r(c(1L, 2L, 3L), c(1L, 2L, 3L), 3L, 3L), 0)
  # Transposition: one transposition = 1 in DL
  expect_equal(Nestimate:::.dl_dist_r(c(1L, 2L), c(2L, 1L), 2L, 2L), 1)
})

# ==============================================================================
# 24. levenshtein zero-length edge cases (L103)
# ==============================================================================

test_that(".levenshtein_dist_r zero-length handles both empty", {
  # both empty → 0
  expect_equal(Nestimate:::.levenshtein_dist_r(integer(), integer(), 0L, 0L), 0)
})

# ==============================================================================
# 25. qgram distance functions direct tests (L141-150, L155-160, L165-173, L178-181, L186-187, L195)
# ==============================================================================

test_that(".get_qgram_r returns empty for n < q", {
  x <- c(1L, 2L)
  result <- Nestimate:::.get_qgram_r(x, 2L, 3L)
  expect_length(result, 0L)
})

test_that(".qgram_dist_r computes correct L1 distance", {
  qx <- c("A\x01B" = 2L, "B\x01C" = 1L)
  qy <- c("A\x01B" = 1L, "C\x01D" = 2L)
  # |2-1| + |1-0| + |0-2| = 1 + 1 + 2 = 4
  expect_equal(Nestimate:::.qgram_dist_r(NULL, NULL, qx, qy), 4L)
})

test_that(".cosine_dist_r is zero for identical profiles", {
  qx <- c("AB" = 3L, "BC" = 2L)
  expect_equal(Nestimate:::.cosine_dist_r(NULL, NULL, qx, qx), 0)
})

test_that(".cosine_dist_r returns 1 when denominator is zero", {
  # Both zero-norm profiles: den == 0
  qx <- integer()
  qy <- integer()
  expect_equal(Nestimate:::.cosine_dist_r(NULL, NULL, qx, qy), 1)
})

test_that(".jaccard_dist_r returns 0 for empty vs empty", {
  expect_equal(Nestimate:::.jaccard_dist_r(integer(), integer()), 0)
})

test_that(".jaccard_dist_r is correct for known inputs", {
  qx <- c("AB" = 1L, "BC" = 1L)
  qy <- c("BC" = 1L, "CD" = 1L)
  # |intersection| / |union| = 1 / 3 → distance = 1 - 1/3 = 2/3
  expect_equal(Nestimate:::.jaccard_dist_r(qx, qy), 2/3)
})

test_that(".jaro_dist_r returns 0 for identical", {
  x <- c(1L, 2L, 3L)
  expect_equal(Nestimate:::.jaro_dist_r(x, x, 3L, 3L), 0)
})

# ==============================================================================
# 26. Cosine distance matrix with zero-norm rows (L319-322)
# ==============================================================================

test_that("cosine distance R path handles all-NA sequences (zero-norm rows)", {
  # Tests the zero_rows path in .cosine_matrix_tdm via the R fallback
  df <- data.frame(
    T1 = c("A", "%", "B"),
    T2 = c("B", "%", "C"),
    T3 = c("C", "%", "A"),
    stringsAsFactors = FALSE
  )
  enc <- Nestimate:::.encode_sequences(df, c("*", "%"))
  dm <- as.matrix(
    Nestimate:::.dissimilarity_matrix_r(enc, "cosine", lambda = 0, q = 2L, p = 0.1)
  )
  # Row 2 is all NA → distance to rows 1 and 3 should be 1
  expect_equal(dm[2, 1], 1)
  expect_equal(dm[2, 3], 1)
  expect_equal(dm[2, 2], 0)
})

# ==============================================================================
# 27. tna / cograph_network covariates rejection (L430, L440)
# ==============================================================================

test_that("cograph_network input with column-name covariates errors", {
  cg <- structure(
    list(data = make_test_data(n = 20, k = 5),
         weights = matrix(0, 4, 4), directed = TRUE),
    class = c("cograph_network", "list")
  )
  expect_error(
    build_clusters(cg, k = 2, covariates = "some_col"),
    "tna/cograph_network"
  )
})

# ==============================================================================
# 28. netobject covariate lookup: missing column error (L477, L485)
# ==============================================================================

test_that("netobject covariate with missing column errors", {
  set.seed(1)
  df <- data.frame(
    T1 = sample(LETTERS[1:3], 30, replace = TRUE),
    T2 = sample(LETTERS[1:3], 30, replace = TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(df, method = "relative")
  expect_error(
    build_clusters(net, k = 2, covariates = "NonExistent"),
    "not found"
  )
})

# ==============================================================================
# 29. Raw data.frame covariate with missing column error (L497-498)
# ==============================================================================

test_that("data.frame covariate missing column errors with helpful message", {
  df <- make_cov_data()
  expect_error(
    build_clusters(df, k = 3, covariates = c("Age", "NotAColumn")),
    "not found|NotAColumn"
  )
})

# ==============================================================================
# 30. Unsupported input type for build_clusters (L755 via .extract_sequence_data)
# ==============================================================================

test_that(".extract_sequence_data errors on unsupported input", {
  expect_error(
    build_clusters(list(a = 1, b = 2), k = 2),
    "Unsupported input"
  )
})

# ==============================================================================
# 31. summary.net_clustering single-member cluster (L755 singleton path)
# ==============================================================================

test_that("summary.net_clustering handles singleton clusters", {
  # Force a single-member cluster by using PAM with extreme data
  df <- data.frame(
    T1 = c("A", "A", "A", "Z"),
    T2 = c("A", "A", "A", "Z"),
    T3 = c("A", "A", "A", "Z"),
    stringsAsFactors = FALSE
  )
  cl <- build_clusters(df, k = 2, dissimilarity = "hamming")
  out <- capture.output(res <- summary(cl))
  expect_true(is.data.frame(res))
})

# ==============================================================================
# 32. .run_covariate_analysis: small-cluster warning (L895-896)
# ==============================================================================

test_that("covariate analysis warns when cluster too small for params", {
  # Create highly imbalanced clusters: 1 observation in one cluster
  # Use tiny n and 2 clusters to force min_cl < n_params
  set.seed(2)
  df <- data.frame(
    T1 = sample(LETTERS[1:3], 8, replace = TRUE),
    T2 = sample(LETTERS[1:3], 8, replace = TRUE),
    T3 = sample(LETTERS[1:3], 8, replace = TRUE),
    X1 = rnorm(8),
    X2 = rnorm(8),
    X3 = rnorm(8),
    X4 = rnorm(8),
    stringsAsFactors = FALSE
  )
  # Use many covariates and small k to trigger the warning
  expect_warning(
    build_clusters(df, k = 2, covariates = c("X1", "X2", "X3", "X4")),
    "parameters|Estimates may be unreliable"
  )
})

# ==============================================================================
# 33. ordered factor covariate gets unordered (L925-926, L930)
# ==============================================================================

test_that("ordered factor covariates are coerced to unordered", {
  df <- make_cov_data()
  df$Level <- factor(sample(c("Low", "Med", "High"), nrow(df), replace = TRUE),
                     levels = c("Low", "Med", "High"), ordered = TRUE)
  cl <- build_clusters(df, k = 3, covariates = "Level")
  expect_true(!is.null(cl$covariates))
})

# ==============================================================================
# 34. .print_covariate_profiles called directly via summary (L950, L959-964)
# ==============================================================================

test_that(".print_covariate_profiles prints both numeric and categorical", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 3, covariates = c("Age", "Gender"))
  out <- capture.output(summary(cl))
  expect_true(any(grepl("Cluster Profiles.*numeric", out)))
  expect_true(any(grepl("Cluster Profiles.*categorical", out)))
})

# ==============================================================================
# 35. .print_covariate_profiles: categorical level missing in a cluster (L1242)
# ==============================================================================

test_that(".print_covariate_profiles handles missing level in a cluster", {
  # Use a rare level that may not appear in every cluster
  set.seed(9)
  df <- data.frame(
    T1 = sample(LETTERS[1:4], 30, replace = TRUE),
    T2 = sample(LETTERS[1:4], 30, replace = TRUE),
    Group = c(rep("X", 25), rep("Y", 5)),  # Y is rare
    stringsAsFactors = FALSE
  )
  cl <- build_clusters(df, k = 3, covariates = "Group")
  out <- capture.output(summary(cl))
  expect_true(any(grepl("0 \\(0%\\)", out)))
})

# ==============================================================================
# 36. .resolve_covariates: no variables specified (L971-972, L982-987, L992-993)
# ==============================================================================

test_that(".resolve_covariates formula with no vars errors", {
  df <- make_cov_data()
  # An empty formula-like string would parse to no vars; use empty char vector
  expect_error(
    build_clusters(df, k = 3, covariates = character(0)),
    "No covariate|specified"
  )
})

# ==============================================================================
# 37. .resolve_covariates: bad type (L997-1000)
# ==============================================================================

test_that(".resolve_covariates rejects non-supported type", {
  df <- make_cov_data()
  expect_error(
    build_clusters(df, k = 3, covariates = 42),
    "formula|character|data.frame"
  )
})

# ==============================================================================
# 38. .resolve_covariates: nrow mismatch via resolved cov_df (L997-1000)
# ==============================================================================

test_that(".resolve_covariates errors on row mismatch from raw data extraction", {
  # This tests the path where raw_data is df but cov_names not present → error
  # with "Available:" message
  df <- make_cov_data()
  expect_error(
    build_clusters(df, k = 3, covariates = "MissingColumn"),
    "MissingColumn|not found"
  )
})

# ==============================================================================
# 39. .run_covariate_analysis: k=2 returns vector not matrix (L1085)
# ==============================================================================

test_that("k=2 covariate analysis still produces a proper coef data.frame", {
  df <- make_cov_data()
  cl <- build_clusters(df, k = 2, covariates = c("Age", "Gender"))
  coefs <- cl$covariates$coefficients
  expect_true(is.data.frame(coefs))
  expect_true("cluster" %in% names(coefs))
  expect_true("variable" %in% names(coefs))
  # For k=2, there is only one row of clusters
  expect_equal(length(unique(coefs$cluster)), 1L)
})

# ==============================================================================
# 40. .run_covariate_analysis: n_dropped > 0 → profiles on complete-case (L1121-1124)
# ==============================================================================

test_that("covariate NA rows use complete-case for profiles", {
  df <- make_cov_data()
  df$Age[1:5] <- NA
  expect_warning(
    cl <- build_clusters(df, k = 3, covariates = "Age"),
    "Dropped 5"
  )
  # Profiles should reflect 35 complete cases
  np <- cl$covariates$profiles$numeric
  total_n <- sum(unique(np[, c("cluster", "n")])$n)
  expect_equal(total_n, 35L)
})


# ==============================================================================
# cluster_network() convenience wrapper (L1435-1452)
# ==============================================================================

test_that("cluster_network with default pam returns netobject_group (L1435-1452)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 50, TRUE), V2 = sample(states, 50, TRUE),
    V3 = sample(states, 50, TRUE), V4 = sample(states, 50, TRUE),
    stringsAsFactors = FALSE
  )
  grp <- cluster_network(data, k = 2)
  expect_true(inherits(grp, "netobject_group"))
  expect_equal(length(grp), 2)
})

test_that("cluster_network with mmm returns netobject_group (L1445-1448)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 50, TRUE), V2 = sample(states, 50, TRUE),
    V3 = sample(states, 50, TRUE), V4 = sample(states, 50, TRUE),
    stringsAsFactors = FALSE
  )
  grp <- cluster_network(data, k = 2, cluster_by = "mmm")
  expect_true(inherits(grp, "netobject_group"))
})

test_that("cluster_network from netobject inherits build_args (L1438-1443)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 50, TRUE), V2 = sample(states, 50, TRUE),
    V3 = sample(states, 50, TRUE), V4 = sample(states, 50, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(data, method = "relative")
  grp <- cluster_network(net, k = 2)
  expect_true(inherits(grp, "netobject_group"))
})

# ------------------------------------------------------------------
# audit_clustering finding #1: cluster_network() must forward
# build_clusters() args (weighted, lambda, seed, q, p, na_syms,
# covariates) instead of silently passing them to build_network().
# ------------------------------------------------------------------

test_that("cluster_network forwards weighted+lambda to build_clusters", {
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 11)
  grp <- cluster_network(data, k = 2, weighted = TRUE, lambda = 2.5)
  cls <- attr(grp, "clustering")
  expect_true(inherits(cls, "net_clustering"))
  expect_true(cls$weighted)
  expect_equal(cls$lambda, 2.5)
})

test_that("cluster_network forwards seed and yields reproducible assignments", {
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 12)
  grp1 <- cluster_network(data, k = 2, seed = 99L)
  grp2 <- cluster_network(data, k = 2, seed = 99L)
  cls1 <- attr(grp1, "clustering")
  cls2 <- attr(grp2, "clustering")
  expect_equal(cls1$seed, 99L)
  expect_equal(cls1$assignments, cls2$assignments)
})

test_that("cluster_network forwards covariates to build_clusters", {
  skip_if_not_installed("nnet")
  data <- make_test_data(n = 40, k = 6, n_states = 3, seed = 13)
  cov_df <- data.frame(grp = factor(rep(c("x", "y"), each = 20)))
  grp <- cluster_network(data, k = 2, covariates = cov_df)
  cls <- attr(grp, "clustering")
  expect_false(is.null(cls$covariates))
})

test_that("cluster_network still routes build_network args to the network step", {
  # Regression guard: the new arg-split must NOT pull build_network params
  # (method, threshold, scaling) into build_clusters(). If a build_network
  # arg leaked into the cluster step, build_clusters() would either ignore
  # it (silently dropping the user intent) or error on an unknown arg.
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 14)
  grp <- cluster_network(data, k = 2, method = "frequency")
  expect_true(inherits(grp, "netobject_group"))
  expect_equal(grp[[1]]$method, "frequency")
  expect_equal(grp[[2]]$method, "frequency")
})

test_that("cluster_network preserves build_args when input is a netobject", {
  # Caller dots should still flow to build_clusters; netobject build_args
  # should still flow to build_network (caller wins on conflicts). This
  # locks in the contract that the split operates on caller dots only,
  # not on build_args.
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 15)
  net <- build_network(data, method = "frequency")
  grp <- cluster_network(net, k = 2, weighted = TRUE, lambda = 1.5)
  cls <- attr(grp, "clustering")
  expect_true(cls$weighted)
  expect_equal(cls$lambda, 1.5)
  # Inherited method survives
  expect_equal(grp[[1]]$method, "frequency")
})

test_that("cluster_network distance branch is equivalent to explicit two-step call", {
  data <- make_test_data(n = 36, k = 6, n_states = 3, seed = 18)

  one_step <- cluster_network(
    data, k = 3, seed = 180L, weighted = TRUE, lambda = 1.25,
    method = "frequency"
  )
  clustered <- build_clusters(
    data, k = 3, seed = 180L, weighted = TRUE, lambda = 1.25
  )
  two_step <- build_network(clustered, method = "frequency")

  expect_equal(attr(one_step, "clustering")$assignments, clustered$assignments)
  expect_equal(lapply(one_step, `[[`, "weights"),
               lapply(two_step, `[[`, "weights"))
})

test_that("cluster_network MMM branch is equivalent to explicit two-step call", {
  data <- make_test_data(n = 36, k = 6, n_states = 3, seed = 19)

  one_step <- cluster_network(
    data, k = 2, cluster_by = "mmm", n_starts = 3, max_iter = 50,
    seed = 190L, method = "frequency"
  )
  mmm <- build_mmm(data, k = 2, n_starts = 3, max_iter = 50, seed = 190L)
  two_step <- build_network(mmm, method = "frequency")

  expect_equal(attr(one_step, "clustering")$assignments, mmm$assignments)
  expect_equal(lapply(one_step, `[[`, "weights"),
               lapply(two_step, `[[`, "weights"))
})

test_that("cluster_network rejects unknown routed arguments", {
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 20)

  expect_error(
    cluster_network(data, k = 2, typo_arg = TRUE),
    "Unknown argument"
  )
  expect_error(
    cluster_network(
      data, k = 2, cluster_by = "mmm", n_starts = 1, max_iter = 5,
      typo_arg = TRUE
    ),
    "Unknown argument"
  )
})

test_that("cluster_network rejects distance options on MMM branch", {
  data <- make_test_data(n = 30, k = 6, n_states = 3, seed = 21)

  expect_error(
    cluster_network(
      data, k = 2, cluster_by = "mmm", dissimilarity = "lv",
      n_starts = 1, max_iter = 5
    ),
    "dissimilarity"
  )
})

# ---- Branch-matrix coverage (task #17) ----
# Crosses dissimilarity x method. All 9 dissimilarities x 5 representative
# methods = 45 cells. Each must produce a valid net_clustering with k
# non-empty clusters. Detects regressions where one dissimilarity silently
# returns NA-filled distance or a clustering backend drops support for a
# particular metric. weighted=FALSE here; weighted-path error is tested
# separately because it only applies to hamming.

test_that("build_clusters branch matrix: dissimilarity x method all succeed", {
  set.seed(17)
  data <- make_test_data(n = 30, k = 8, n_states = 3)
  k_target <- 3L

  # Keep the method list small and representative: one partitional (pam),
  # one agglomerative per linkage family (ward.D2 minvar, complete, average,
  # single). Skipping centroid/median/mcquitty — they exercise the same
  # hclust branch as the ones already covered.
  methods_subset <- c("pam", "ward.D2", "complete", "average", "single")

  grid <- expand.grid(
    dissimilarity = Nestimate:::.clustering_metrics,
    method        = methods_subset,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    cfg <- grid[i, ]
    info <- sprintf("dissim=%s method=%s", cfg$dissimilarity, cfg$method)
    fit <- tryCatch(
      build_clusters(
        data, k = k_target,
        dissimilarity = cfg$dissimilarity,
        method        = cfg$method,
        seed          = 17L
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) {
      fail(sprintf("%s -> %s", info, conditionMessage(fit)))
      next
    }
    expect_true(inherits(fit, "net_clustering"), info = info)
    # The fit must produce exactly k_target clusters and every cluster
    # must have at least one assigned sequence.
    clust <- fit$clusters %||% fit$assignments %||% fit$cluster
    if (is.null(clust)) clust <- fit[[which(vapply(fit, function(el)
      is.integer(el) && length(el) == nrow(data), logical(1)))[1]]]
    expect_equal(length(unique(clust)), k_target, info = info)
    expect_true(all(tabulate(clust, nbins = k_target) > 0), info = info)
  }
})

test_that("build_clusters weighted=TRUE errors on non-hamming dissimilarity", {
  data <- make_test_data(n = 20, k = 6, n_states = 3)
  # Only hamming supports weighted=TRUE — every other metric must error.
  for (d in setdiff(Nestimate:::.clustering_metrics, "hamming")) {
    expect_error(
      build_clusters(data, k = 2, dissimilarity = d, weighted = TRUE),
      "Weighting is only supported for Hamming",
      info = sprintf("dissim=%s", d)
    )
  }
})
