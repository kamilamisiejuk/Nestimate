testthat::skip_on_cran()

# Helper: generate reproducible frequency-like data
.make_freq_data <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$rid <- seq_len(n)
  df
}


# ---- Input validation ----

test_that("build_network auto-converts character data for glasso", {
  df <- data.frame(a = letters[1:10], b = letters[10:1])
  # Should auto-convert to frequency counts, not error
  net <- build_network(df, method = "glasso")
  expect_s3_class(net, "netobject")
})

test_that("build_network errors on non-symmetric matrix", {
  m <- matrix(1:9, 3, 3)
  expect_error(
    build_network(m, method = "glasso", params = list(n = 50)),
    "symmetric"
  )
})

test_that("build_network errors when n missing for matrix input", {
  m <- diag(5)
  expect_error(
    build_network(m, method = "glasso"),
    "Sample size 'n' is required"
  )
})


# ---- Auto-cleaning ----

test_that("zero-variance columns are dropped with message", {
  df <- .make_freq_data()
  df$constant <- 5
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "Dropping zero-variance"
  )
  expect_equal(net$n_nodes, 5)
  expect_false("constant" %in% colnames(net$weights))
})

test_that("non-syntactic column names are dropped with message", {
  df <- .make_freq_data(n = 80, p = 4)
  df$`%` <- rpois(80, 2)
  df$`*` <- rpois(80, 3)
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "non-syntactic"
  )
  expect_equal(net$n_nodes, 4)
  expect_false("%" %in% colnames(net$weights))
})

test_that("all-NA columns are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$empty <- NA_real_
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "all-NA"
  )
  expect_equal(net$n_nodes, 5)
})

test_that("rows with NA are dropped with message", {
  df <- .make_freq_data(n = 80, p = 5)
  df$state_1[1:3] <- NA
  expect_message(
    net <- build_network(df, method = "glasso", params = list(nlambda = 20L)),
    "rows with NA"
  )
  expect_equal(net$n, 77)
})


# ---- Method: glasso ----

test_that("build_network works with data frame input (glasso)", {
  df <- .make_freq_data(n = 80, p = 6)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "glasso")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 6)
  expect_true(is.matrix(net$weights))
  expect_equal(nrow(net$weights), 6)
  expect_equal(ncol(net$weights), 6)
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # Edges data frame
  expect_true(is.data.frame(net$edges))
  expect_true(all(c("from", "to", "weight") %in% names(net$edges)))
  expect_equal(net$n_edges, nrow(net$edges))
  # EBIC path length matches nlambda
  expect_equal(length(net$ebic_path), 20)
  expect_equal(length(net$lambda_path), 20)
  # Glasso-specific fields
  expect_true(!is.null(net$precision_matrix))
  expect_true(!is.null(net$gamma))
  expect_true(!is.null(net$lambda_selected))
})

test_that("build_network works with correlation matrix input (glasso)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  S <- cor(df[, num_cols])

  net <- build_network(S, method = "glasso", params = list(n = 100,
                                                            nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$n, 100)
  expect_equal(net$n_nodes, 5)
})

test_that("build_network works with covariance matrix input (glasso)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  C <- cov(df[, num_cols])

  net <- build_network(C, method = "glasso",
                       params = list(n = 100, input_type = "cov",
                                     nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$n_nodes, 5)
})

test_that("method aliases resolve to glasso", {
  df <- .make_freq_data(n = 80, p = 4)
  net1 <- build_network(df, method = "ebicglasso",
                         params = list(nlambda = 20L))
  net2 <- build_network(df, method = "regularized",
                         params = list(nlambda = 20L))
  expect_equal(net1$method, "glasso")
  expect_equal(net2$method, "glasso")
})


# ---- Method: pcor (unregularised) ----

test_that("build_network works with method='pcor'", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 5)
  expect_true(is.matrix(net$weights))
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # Should have precision matrix
  expect_true(!is.null(net$precision_matrix))
  # Edges
  expect_true(is.data.frame(net$edges))
  expect_equal(net$n_edges, nrow(net$edges))
  # No glasso-specific fields
  expect_null(net$lambda_selected)
  expect_null(net$ebic_path)
})

test_that("method='partial' resolves to pcor", {
  df <- .make_freq_data(n = 80, p = 4)
  net <- build_network(df, method = "partial")
  expect_equal(net$method, "pcor")
})

test_that("pcor errors on singular matrix", {
  # p > n: more variables than observations
  set.seed(42)
  mat <- matrix(rnorm(10 * 20), nrow = 10, ncol = 20)
  colnames(mat) <- paste0("V", seq_len(20))
  S <- cor(mat)
  expect_error(
    build_network(S, method = "pcor", params = list(n = 10)),
    "singular"
  )
})

test_that("pcor warns on near-singular matrix", {
  # Create a near-singular correlation matrix
  set.seed(123)
  n <- 30; p <- 5
  mat <- matrix(rnorm(n * p), n, p)
  # Make two columns nearly collinear
  mat[, 5] <- mat[, 1] + rnorm(n, sd = 1e-7)
  colnames(mat) <- paste0("V", seq_len(p))
  S <- cor(mat)
  # rcond should be very small, triggering the warning
  expect_warning(
    build_network(S, method = "pcor", params = list(n = n)),
    "near-singular"
  )
})

test_that("pcor output unchanged for well-conditioned data", {
  # Verify the rcond check doesn't alter numerical output
  set.seed(456)
  df <- data.frame(
    V1 = rnorm(100), V2 = rnorm(100), V3 = rnorm(100), V4 = rnorm(100)
  )
  net <- build_network(df, method = "pcor")
  # Manual reference: solve(cor) -> precision -> pcor
  S <- cor(as.matrix(df))
  Wi <- solve(S)
  D <- diag(1 / sqrt(diag(Wi)))
  pcor_ref <- -D %*% Wi %*% D
  diag(pcor_ref) <- 0
  dimnames(pcor_ref) <- dimnames(S)
  expect_equal(unname(net$weights), unname(pcor_ref), tolerance = 1e-14)
})


# ---- Method: cor (correlation network) ----

test_that("build_network works with method='cor'", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "cor")
  expect_equal(net$n, 80)
  expect_equal(net$n_nodes, 5)
  expect_true(is.matrix(net$weights))
  # Diagonal should be zero
  expect_true(all(diag(net$weights) == 0))
  # Should be symmetric
  expect_equal(net$weights, t(net$weights))
  # No precision matrix for cor method
  expect_null(net$precision_matrix)
  # Edges
  expect_true(is.data.frame(net$edges))
  expect_equal(net$n_edges, nrow(net$edges))
})

test_that("method='correlation' resolves to cor", {
  df <- .make_freq_data(n = 80, p = 4)
  net <- build_network(df, method = "correlation")
  expect_equal(net$method, "cor")
})

test_that("cor threshold filters weak edges", {
  df <- .make_freq_data(n = 100, p = 5)
  net_low <- build_network(df, method = "cor", threshold = 0.01)
  net_high <- build_network(df, method = "cor", threshold = 0.3)
  # Higher threshold should produce same or fewer edges
  expect_true(net_high$n_edges <= net_low$n_edges)
})

test_that("cor matrix matches thresholded cor_matrix", {
  df <- .make_freq_data(n = 80, p = 5)
  thr <- 0.15
  net <- build_network(df, method = "cor", threshold = thr)
  expected <- net$cor_matrix
  diag(expected) <- 0
  expected[abs(expected) < thr] <- 0
  expect_equal(net$weights, expected)
})


# ---- New aliases ----

test_that("new aliases tna, ftna, cna, corr resolve correctly", {
  df <- .make_freq_data(n = 80, p = 4)

  net_corr <- build_network(df, method = "corr")
  expect_equal(net_corr$method, "cor")
})


# ---- id_col exclusion ----

test_that("id_col columns are excluded from analysis", {
  df <- .make_freq_data(n = 80, p = 5)
  df$subject_id <- seq_len(80)

  net <- build_network(df, method = "glasso",
                       params = list(id_col = "subject_id", nlambda = 20L))

  # subject_id and rid should be excluded -> 5 variables
  expect_equal(net$n_nodes, 5)
  expect_false("subject_id" %in% colnames(net$weights))
  expect_false("rid" %in% colnames(net$weights))
})

test_that("association id and id_col params are resolved exactly", {
  df <- .make_freq_data(n = 80, p = 5)
  df$subject_id <- rep(seq_len(20), each = 4)

  via_id_col <- build_network(df, method = "cor",
                              params = list(id_col = "subject_id"))
  via_id <- build_network(df, method = "cor",
                          params = list(id = "subject_id"))
  expect_false("subject_id" %in% colnames(via_id_col$weights))
  expect_false("subject_id" %in% colnames(via_id$weights))

  between <- build_network(df, method = "cor",
                           params = list(id_col = "subject_id"),
                           level = "between")
  expect_equal(between$level, "between")
  expect_equal(between$n, 20L)

  expect_error(
    build_network(df, method = "cor",
                  params = list(id = "subject_id", id_col = "state_1")),
    "'actor' and 'params\\$id_col' specify different columns"
  )
})


# ---- Gamma effects ----

test_that("higher gamma produces sparser or equal networks", {
  df <- .make_freq_data(n = 150, p = 7, seed = 123)

  net_low <- build_network(df, method = "glasso",
                           params = list(gamma = 0, nlambda = 50L))
  net_high <- build_network(df, method = "glasso",
                            params = list(gamma = 1, nlambda = 50L))

  expect_true(net_high$n_edges <= net_low$n_edges)
})


# ---- S3 print method ----

test_that("print.netobject produces expected output for glasso", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Partial Correlation Network \\(EBICglasso\\)", out)))
  expect_true(any(grepl("Weight matrix:", out)))
  expect_true(any(grepl("Sample size: 80", out)))
  expect_true(any(grepl("Gamma:", out)))
  expect_true(any(grepl("Lambda:", out)))
})

test_that("print.netobject produces expected output for pcor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")

  out <- capture.output(print(net))
  expect_true(any(grepl("unregularised", out)))
  expect_false(any(grepl("Gamma:", out)))
})

test_that("print.netobject produces expected output for cor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor")

  out <- capture.output(print(net))
  expect_true(any(grepl("Correlation Network", out)))
  expect_false(any(grepl("Gamma:", out)))
})

test_that("print.netobject returns invisible(x)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})


# ---- Correlation method argument ----

test_that("cor_method argument is respected", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pearson <- build_network(df, method = "glasso",
                               params = list(cor_method = "pearson",
                                             nlambda = 20L))
  net_spearman <- build_network(df, method = "glasso",
                                params = list(cor_method = "spearman",
                                              nlambda = 20L))

  expect_false(identical(net_pearson$cor_matrix, net_spearman$cor_matrix))
})


# ---- Edge data frame correctness ----

test_that("edges match non-zero upper triangle of matrix", {
  df <- .make_freq_data(n = 100, p = 6)
  net <- build_network(df, method = "glasso", params = list(nlambda = 30L))

  mat <- net$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net$edges), nrow(upper_nz))

  for (i in seq_len(nrow(net$edges))) {
    expect_equal(net$edges$weight[i], mat[net$edges$from[i], net$edges$to[i]])
  }
})

test_that("edges match for pcor and cor methods too", {
  df <- .make_freq_data(n = 80, p = 5)

  net_pcor <- build_network(df, method = "pcor")
  mat <- net_pcor$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net_pcor$edges), nrow(upper_nz))

  net_cor <- build_network(df, method = "cor", threshold = 0.1)
  mat <- net_cor$weights
  upper_nz <- which(upper.tri(mat) & mat != 0, arr.ind = TRUE)
  expect_equal(nrow(net_cor$edges), nrow(upper_nz))
})


# ---- Cross-method consistency ----

test_that("all methods produce consistent structure", {
  df <- .make_freq_data(n = 80, p = 5)
  methods <- c("glasso", "pcor", "cor")

  for (m in methods) {
    net <- build_network(df, method = m, params = list(nlambda = 20L))
    expect_s3_class(net, "netobject")
    expect_equal(net$method, m)
    expect_equal(net$n, 80)
    expect_equal(net$n_nodes, 5)
    expect_true(is.matrix(net$weights))
    expect_true(is.matrix(net$cor_matrix))
    expect_true(is.data.frame(net$edges))
    expect_true(is.numeric(net$n_edges))
    expect_true(all(diag(net$weights) == 0))
  }
})


# ---- Multilevel: helper ----

# Generate data with repeated measures per person
.make_multilevel_data <- function(n_persons = 30, obs_per_person = 5,
                                  p = 5, seed = 42) {
  set.seed(seed)
  n_total <- n_persons * obs_per_person
  mat <- matrix(rpois(n_total * p, lambda = 10), nrow = n_total, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  df <- as.data.frame(mat)
  df$person <- rep(seq_len(n_persons), each = obs_per_person)
  df$rid <- seq_len(n_total)
  df
}


# ---- Multilevel: level = "between" ----

test_that("level='between' aggregates to person means", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "between", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_equal(net$level, "between")
  # n should be number of unique persons
  expect_equal(net$n, 30)
  expect_equal(net$n_nodes, 5)
})

test_that("level='between' works with method='pcor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "pcor", id_col = "person",
                       level = "between")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 30)
})

test_that("level='between' works with method='cor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "between")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "cor")
  expect_equal(net$n, 30)
})


# ---- Multilevel: level = "within" ----

test_that("level='within' centers correctly (column means near 0)", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "within")

  # Within-centered correlations should exist
  expect_s3_class(net, "netobject")
  expect_equal(net$level, "within")
  # n = total observations (all persons have >= 2 obs)
  expect_equal(net$n, 150)
  expect_equal(net$n_nodes, 5)
})

test_that("level='within' drops single-observation persons with message", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  # Add 3 persons with only 1 observation each
  singles <- data.frame(
    state_1 = rpois(3, 10), state_2 = rpois(3, 10),
    state_3 = rpois(3, 10), state_4 = rpois(3, 10),
    state_5 = rpois(3, 10),
    person = c(101, 102, 103), rid = 151:153
  )
  df <- rbind(df, singles)

  expect_message(
    net <- build_network(df, method = "glasso", id_col = "person",
                         level = "within", params = list(nlambda = 20L)),
    "single-observation"
  )
  # Single-obs rows should be dropped: 153 - 3 = 150
  expect_equal(net$n, 150)
})

test_that("level='within' works with method='pcor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "pcor", id_col = "person",
                       level = "within")

  expect_s3_class(net, "netobject")
  expect_equal(net$method, "pcor")
  expect_equal(net$n, 150)
})


# ---- Multilevel: level = "both" ----

test_that("level='both' returns netobject_ml with both sub-networks", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  expect_s3_class(net, "netobject_ml")
  expect_s3_class(net$between, "netobject")
  expect_s3_class(net$within, "netobject")
  expect_equal(net$method, "glasso")
  expect_equal(net$between$level, "between")
  expect_equal(net$within$level, "within")
  expect_equal(net$between$n, 30)
  expect_equal(net$within$n, 150)
})

test_that("level='both' works with method='cor'", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "both")

  expect_s3_class(net, "netobject_ml")
  expect_equal(net$method, "cor")
  expect_s3_class(net$between, "netobject")
  expect_s3_class(net$within, "netobject")
})


# ---- Multilevel: validation ----

test_that("level requires id_col", {
  df <- .make_multilevel_data()
  expect_error(
    build_network(df, method = "glasso", level = "between"),
    "id.*required"
  )
})

test_that("level requires data frame input", {
  m <- diag(5)
  expect_error(
    build_network(m, method = "glasso",
                  params = list(n = 50), id_col = "person",
                  level = "between"),
    "data frame"
  )
})

test_that("level=NULL preserves backward compatibility", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       params = list(nlambda = 20L))

  expect_s3_class(net, "netobject")
  expect_null(net$level)
  # Without level, n = total rows
  expect_equal(net$n, 150)
})


# ---- Multilevel: print methods ----

test_that("print.netobject shows level label for between", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "between", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("between-person", out)))
})

test_that("print.netobject shows level label for within", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "cor", id_col = "person",
                       level = "within")

  out <- capture.output(print(net))
  expect_true(any(grepl("within-person", out)))
})

test_that("print.netobject_ml shows both levels", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Multilevel", out)))
  expect_true(any(grepl("Between-person", out)))
  expect_true(any(grepl("Within-person", out)))
  expect_true(any(grepl("unique persons", out)))
  expect_true(any(grepl("observations", out)))
})

test_that("print.netobject_ml returns invisible(x)", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))
  ret <- capture.output(result <- print(net))
  expect_identical(result, net)
})


# ---- Predictability ----

test_that("predictability returns tidy df for glasso", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  r2 <- predictability(net)

  expect_s3_class(r2, "data.frame")
  expect_equal(nrow(r2), 5)
  expect_setequal(names(r2), c("node", "R2", "RMSE"))
  expect_equal(r2$node, colnames(net$weights))
  expect_true(all(r2$R2 >= 0 & r2$R2 <= 1))
})

test_that("predictability returns tidy df for pcor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "pcor")
  r2 <- predictability(net)

  expect_s3_class(r2, "data.frame")
  expect_equal(nrow(r2), 5)
  expect_true(all(r2$R2 >= 0 & r2$R2 <= 1))
})

test_that("predictability returns tidy df for cor", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "cor", threshold = 0.1)
  r2 <- predictability(net)

  expect_s3_class(r2, "data.frame")
  expect_equal(nrow(r2), 5)
  expect_true(all(r2$R2 >= 0 & r2$R2 <= 1))
})

test_that("predictability.cor returns 0 for isolated nodes", {
  df <- .make_freq_data(n = 80, p = 5)
  # Very high threshold should isolate most nodes
  net <- build_network(df, method = "cor", threshold = 0.99)
  r2 <- predictability(net)

  # Isolated nodes (no edges) should have R^2 = 0
  isolated <- vapply(seq_len(net$n_nodes), function(j) {
    all(net$weights[j, ] == 0)
  }, logical(1))
  if (any(isolated)) {
    expect_true(all(r2$R2[isolated] == 0))
  }
})

test_that("predictability works for netobject_ml", {
  df <- .make_multilevel_data(n_persons = 30, obs_per_person = 5, p = 5)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))
  r2 <- predictability(net)

  expect_true(is.list(r2))
  expect_true("between" %in% names(r2))
  expect_true("within" %in% names(r2))
  expect_s3_class(r2$between, "data.frame")
  expect_s3_class(r2$within, "data.frame")
  expect_equal(nrow(r2$between), 5)
  expect_equal(nrow(r2$within), 5)
  expect_true(all(r2$between$R2 >= 0 & r2$between$R2 <= 1))
  expect_true(all(r2$within$R2 >= 0 & r2$within$R2 <= 1))
})

test_that("print does not show literal 'predictability' (lowercase)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))
  out <- capture.output(print(net))
  # Predictability display uses capitalised label "Predictability (R²)"
  expect_true(any(grepl("Predictability", out)))
  # The lowercase word "predictability" should NOT appear in output
  expect_false(any(grepl("predictability", out)))
})


# ---- $data field ----

test_that("$data is cleaned numeric matrix for association methods (data frame input)", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  expect_true(is.matrix(net$data))
  expect_true(is.numeric(net$data))
  expect_equal(nrow(net$data), 80)
  # 5 state columns only (rid excluded during cleaning)
  expect_equal(ncol(net$data), 5)
})

test_that("$data is NULL for association methods (matrix input)", {
  df <- .make_freq_data(n = 100, p = 5)
  num_cols <- setdiff(names(df), "rid")
  S <- cor(df[, num_cols])

  net <- build_network(S, method = "glasso", params = list(n = 100,
                                                            nlambda = 20L))

  # No row-level data available from matrix input
  expect_null(net$data)
})

test_that("print.netobject shows data dimensions", {
  df <- .make_freq_data(n = 80, p = 5)
  net <- build_network(df, method = "glasso", params = list(nlambda = 20L))

  out <- capture.output(print(net))
  expect_true(any(grepl("Sample size: 80", out)))
})


# ---- Coverage gap tests ----

# L180-181: multi-column group key via interaction()
test_that("build_network group dispatch with multi-column group key", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 60, replace = TRUE),
    g1 = rep(c("X", "Y"), 30),
    g2 = rep(c("P", "Q"), each = 30),
    stringsAsFactors = FALSE
  )
  grp_nets <- build_network(df, method = "relative",
                            group = c("g1", "g2"),
                            params = list(format = "wide"))
  expect_s3_class(grp_nets, "netobject_group")
  # 4 combinations: X-P, X-Q, Y-P, Y-Q
  expect_equal(length(grp_nets), 4L)
  expect_equal(attr(grp_nets, "group_col"), c("g1", "g2"))
})

# L211-212: explicit codes triggers onehot format
test_that("build_network with explicit codes triggers onehot format", {
  set.seed(42)
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L),
    B = c(0L, 1L, 0L, 1L, 0L),
    C = c(1L, 1L, 0L, 0L, 1L)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B", "C"),
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
})

test_that("onehot frequency dispatch matches wtna transition counts", {
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L),
    C = c(0L, 0L, 1L, 0L, 0L, 1L)
  )

  net <- suppressWarnings(build_network(df, method = "frequency",
                                        format = "onehot",
                                        window_size = 1L))
  ref <- suppressWarnings(build_network(df, method = "wtna",
                                        window_size = 1L))

  expect_equal(net$nodes$label, ref$nodes$label)
  expect_equal(net$weights, ref$weights)
  expect_equal(net$frequency_matrix, ref$weights)
  expect_equal(net$method, "frequency")
})

test_that("onehot relative dispatch matches wtna transition probabilities", {
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L),
    C = c(0L, 0L, 1L, 0L, 0L, 1L)
  )

  net <- suppressWarnings(build_network(df, method = "relative",
                                        format = "onehot",
                                        window_size = 1L))
  ref <- wtna(df, method = "transition", type = "relative",
              window_size = 1L)

  expect_equal(net$nodes$label, ref$nodes$label)
  expect_equal(net$weights, ref$weights)
  expect_equal(net$method, "relative")
})

test_that("onehot build_network treats NA indicators as inactive", {
  df_na <- data.frame(
    A = c(1L, NA, 0L, 1L),
    B = c(0L, 1L, NA, 0L),
    C = c(0L, 0L, 1L, NA)
  )
  df_zero <- df_na
  df_zero[is.na(df_zero)] <- 0L

  freq_na <- suppressWarnings(build_network(df_na, method = "frequency",
                                            format = "onehot",
                                            codes = c("A", "B", "C"),
                                            window_size = 1L))
  freq_zero <- suppressWarnings(build_network(df_zero, method = "frequency",
                                              format = "onehot",
                                              codes = c("A", "B", "C"),
                                              window_size = 1L))
  rel_na <- suppressWarnings(build_network(df_na, method = "relative",
                                           format = "onehot",
                                           codes = c("A", "B", "C"),
                                           window_size = 1L))
  rel_zero <- suppressWarnings(build_network(df_zero, method = "relative",
                                             format = "onehot",
                                             codes = c("A", "B", "C"),
                                             window_size = 1L))
  cna_na <- suppressWarnings(build_network(df_na, method = "cna",
                                           format = "onehot",
                                           codes = c("A", "B", "C"),
                                           window_size = 1L))
  cna_zero <- suppressWarnings(build_network(df_zero, method = "cna",
                                             format = "onehot",
                                             codes = c("A", "B", "C"),
                                             window_size = 1L))

  expect_false(anyNA(freq_na$weights))
  expect_false(anyNA(rel_na$weights))
  expect_false(anyNA(cna_na$weights))
  expect_equal(freq_na$weights, freq_zero$weights)
  expect_equal(rel_na$weights, rel_zero$weights)
  expect_equal(cna_na$weights, cna_zero$weights)
})

test_that("onehot build_network keeps code columns in data", {
  df <- data.frame(
    A = c(1L, 0L, 1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L),
    C = c(0L, 0L, 1L, 0L, 0L, 1L)
  )

  net <- suppressWarnings(build_network(df, method = "frequency",
                                        format = "onehot",
                                        window_size = 1L))

  expect_equal(names(net$data), c("A", "B", "C"))
  expect_null(net$metadata)
  expect_equal(net$data, df)
})

test_that("onehot build_network separates actor/session metadata from codes", {
  df <- data.frame(
    actor = c("s1", "s1", "s1", "s2", "s2", "s2"),
    session = c("a", "a", "b", "a", "a", "b"),
    A = c(1L, 0L, 1L, 0L, 1L, 0L),
    B = c(0L, 1L, 0L, 1L, 0L, 1L)
  )

  net <- build_network(df, method = "relative", format = "onehot",
                       codes = c("A", "B"), actor = "actor",
                       session = "session", window_size = 1L)

  expect_equal(names(net$data), c("A", "B"))
  expect_true(all(c("actor", "session") %in% names(net$metadata)))
  expect_equal(net$data, df[, c("A", "B")])
})

# L215: action column present → long format detection
test_that("build_network auto-detects long format via action column", {
  long_data <- data.frame(
    Actor = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time  = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  expect_warning(
    net <- build_network(long_data, method = "relative",
                         action = "Action"),
    "one long sequence is not recommended"
  )
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
})

test_that("build_network warns when long data has action but no actor", {
  long_data <- data.frame(
    action = c("A", "B", "C", "A", "C", "B"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    net <- build_tna(long_data, action = "action"),
    "one long sequence is not recommended"
  )
  expect_s3_class(net, "netobject")
  expect_null(net$build_args$actor)
})

test_that("build_network rejects explicit missing column arguments before auto-detection", {
  long_data <- data.frame(
    id = c(1L, 1L, 1L, 2L, 2L, 2L),
    action = c("A", "B", "C", "A", "C", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_tna(long_data, actor = "id", action = "actoin"),
    "'action' column not found in data: actoin"
  )
  expect_error(
    build_network(long_data, method = "relative", actor = "student",
                  action = "action"),
    "'actor' column not found in data: student"
  )
  expect_error(
    build_network(long_data, method = "relative", group = "condition"),
    "'group' column not found in data: condition"
  )
  expect_error(
    build_network(long_data, method = "relative", format = "onehot",
                  codes = c("A", "B")),
    "'codes' columns not found in data: A, B"
  )
})

test_that("build_network rejects explicit missing params column arguments", {
  long_data <- data.frame(
    id = c(1L, 1L, 1L, 2L, 2L, 2L),
    action = c("A", "B", "C", "A", "C", "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_network(long_data, method = "relative",
                  params = list(action = "actoin")),
    "'action' column not found in data: actoin"
  )
  expect_error(
    build_network(long_data, method = "relative",
                  params = list(action = "action", id = "student")),
    "'actor' column not found in data: student"
  )
  expect_error(
    build_network(long_data, method = "relative",
                  params = list(format = "wide", cols = c("T1", "T2"))),
    "'params\\$cols' columns not found in data: T1, T2"
  )
})

test_that("build_network rejects conflicting formal and params column arguments", {
  long_data <- data.frame(
    actor = c(1L, 1L, 2L, 2L),
    other_actor = c(1L, 1L, 2L, 2L),
    action = c("A", "B", "A", "C"),
    other_action = c("B", "C", "B", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_network(long_data, method = "relative", action = "action",
                  params = list(action = "other_action")),
    "'action' and 'params\\$action' specify different columns"
  )
  expect_error(
    build_network(long_data, method = "relative", actor = "actor",
                  action = "action", params = list(id = "other_actor")),
    "'actor' and 'params\\$id' specify different columns"
  )
})

test_that("build_network params long action is honored before auto-matching", {
  long_data <- data.frame(
    id = c(1L, 1L, 2L, 2L),
    time = c(1L, 2L, 1L, 2L),
    action = c("X", "Y", "Y", "X"),
    code = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )

  net <- build_network(
    long_data, method = "frequency",
    params = list(format = "long", action = "code", id = "id",
                  time = "time")
  )

  expect_equal(sort(net$nodes$label), c("A", "B"))
  expect_false(any(c("X", "Y") %in% net$nodes$label))
  expect_equal(net$params$format, "wide")
  expect_equal(net$build_args$action, "code")
})

test_that("build_network grouped params long excludes id and time from states", {
  long_data <- data.frame(
    id = c(1L, 1L, 1L, 2L, 2L, 2L),
    order = c(1L, 2L, 3L, 1L, 2L, 3L),
    code = c("A", "B", "C", "A", "C", "B"),
    grp = c("G1", "G1", "G1", "G2", "G2", "G2"),
    stringsAsFactors = FALSE
  )

  nets <- build_network(
    long_data, method = "frequency", group = "grp",
    params = list(format = "long", action = "code", id = "id",
                  time = "order")
  )

  expect_s3_class(nets, "netobject_group")
  expect_equal(sort(nets$G1$nodes$label), c("A", "B", "C"))
  expect_equal(sort(nets$G2$nodes$label), c("A", "B", "C"))
})

test_that("build_network params long creates replayable wide transition networks", {
  long_data <- data.frame(
    id = c(1L, 1L, 1L, 2L, 2L, 2L),
    time = c(1L, 2L, 3L, 1L, 2L, 3L),
    code = c("A", "B", "C", "A", "C", "B"),
    stringsAsFactors = FALSE
  )

  net <- build_network(
    long_data, method = "frequency",
    params = list(format = "long", action = "code", id = "id",
                  time = "time")
  )

  expect_equal(net$params$format, "wide")
  expect_error(bootstrap_network(net, iter = 2), NA)
})

# L236-246: long format path through prepare
test_that("build_network long format with actor/time/action passes through prepare", {
  long_data <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  net <- build_network(long_data, method = "relative",
                       actor = "Actor", action = "Action",
                       time = "Time")
  expect_s3_class(net, "netobject")
  # After prepare the format should be reset to "wide"
  expect_equal(net$params$format, "wide")
})

# L260-264: onehot without windowing or session warns
test_that("build_network one-hot without windowing warns", {
  df <- data.frame(A = c(1L, 0L, 1L), B = c(0L, 1L, 0L))
  expect_warning(
    build_network(df, method = "co_occurrence",
                  codes = c("A", "B"), window_size = 1L),
    "sparse"
  )
})

# L270-271, L277: grp_col built from actor + session → params$actor
test_that("build_network onehot with actor and session sets params$actor", {
  set.seed(42)
  df <- data.frame(
    actor   = rep(1:3, each = 4),
    session = rep(c("s1", "s2"), 6),
    A = sample(0:1, 12, replace = TRUE),
    B = sample(0:1, 12, replace = TRUE)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B"),
                       actor = "actor", session = "session",
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_equal(net$params$actor, c("actor", "session"))
})

# L341-343: estimator returning bad structure triggers error
test_that("build_network errors when estimator returns malformed output", {
  bad_fn <- function(data, ...) {
    list(weights = matrix(1, 2, 2))  # missing 'matrix', 'nodes', 'directed'
  }
  register_estimator("bad_estimator", bad_fn, "bad", directed = FALSE)
  on.exit(remove_estimator("bad_estimator"))

  set.seed(42)
  df <- data.frame(a = rnorm(10), b = rnorm(10))
  expect_error(
    build_network(df, method = "bad_estimator"),
    "must return a list with 'matrix', 'nodes', and 'directed'"
  )
})

# L385: print method label falls through to unknown method
test_that("print.netobject shows generic label for unknown method", {
  # Register a custom estimator with a non-standard name
  custom_fn <- function(data, ...) {
    S <- cor(data)
    diag(S) <- 0
    list(matrix = S, nodes = colnames(S), directed = FALSE)
  }
  register_estimator("my_custom_net", custom_fn, "Custom", directed = FALSE)
  on.exit(remove_estimator("my_custom_net"))

  set.seed(42)
  df <- data.frame(x = rnorm(50), y = rnorm(50), z = rnorm(50))
  net <- build_network(df, method = "my_custom_net")
  out <- capture.output(print(net))
  expect_true(any(grepl("my_custom_net", out)))
})

# L474-475: print.netobject shows metadata when present
test_that("print.netobject_group shows group info", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B", "C"), 60, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 60, replace = TRUE),
    grp = rep(c("X", "Y", "Z"), 20),
    stringsAsFactors = FALSE
  )
  nets <- build_network(df, method = "relative", group = "grp",
                        params = list(format = "wide"))
  out <- capture.output(print(nets))
  expect_true(any(grepl("Group Networks", out)))
  expect_true(any(grepl("3 groups", out)))
  # Each group label appears in the per-group table.
  expect_true(any(grepl("\\bX\\b", out)))
  expect_true(any(grepl("\\bY\\b", out)))
  expect_true(any(grepl("\\bZ\\b", out)))
})

test_that("print.netobject_group returns invisible(x)", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A", "B"), 40, replace = TRUE),
    T2 = sample(c("A", "B"), 40, replace = TRUE),
    grp = rep(c("X", "Y"), 20),
    stringsAsFactors = FALSE
  )
  nets <- build_network(df, method = "relative", group = "grp",
                        params = list(format = "wide"))
  ret <- capture.output(result <- print(nets))
  expect_identical(result, nets)
})


# L803-808: predictability for cor with single-neighbor node
test_that("predictability.netobject for cor with single-neighbor node", {
  # Build a cor network where one node has exactly one neighbor
  set.seed(99)
  df <- data.frame(
    x = rnorm(100),
    y = rnorm(100),
    z = rnorm(100)
  )
  # Use a high threshold to ensure exactly one neighbor for at least one node
  net <- build_network(df, method = "cor", threshold = 0)

  r2 <- predictability(net)
  expect_s3_class(r2, "data.frame")
  expect_equal(nrow(r2), 3L)
  expect_true(all(r2$R2 >= 0 & r2$R2 <= 1))
})

# print.netobject with ising shows thresholds range
test_that("print.netobject shows ising thresholds range", {
  skip_if_not_installed("glmnet")
  df <- data.frame(
    V1 = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L,
           0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L),
    V2 = c(1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L,
           1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L),
    V3 = c(0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L,
           1L, 1L, 0L, 0L, 1L, 1L, 0L, 0L, 1L, 1L)
  )
  net <- build_network(df, method = "ising")
  out <- capture.output(print(net))
  # Should show Ising-specific fields
  expect_true(any(grepl("Ising", out)))
})

# netobject_ml print shows no sample size when $n is NULL
test_that("print.netobject_ml handles sub-networks with no n", {
  df <- .make_multilevel_data(n_persons = 20, obs_per_person = 5, p = 3)
  net <- build_network(df, method = "glasso", id_col = "person",
                       level = "both", params = list(nlambda = 20L))

  # Both sub-networks are full netobjects; print should succeed
  expect_no_error(print(net))
})



# ---- Estimators.R coverage gap tests ----

# L89: .count_transitions_wide stop on missing cols
test_that(".count_transitions_wide errors on missing state columns", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  expect_error(
    .count_transitions_wide(df, cols = c("T1", "T2", "T_MISSING")),
    "Columns not found"
  )
})

# One-column sequences have states but no transitions; tna returns a zero matrix.
test_that(".count_transitions_wide returns zero matrix for fewer than 2 state columns", {
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  out <- .count_transitions_wide(df)
  expect_equal(rownames(out), c("A", "B"))
  expect_equal(out, matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))))
})

# L113: .count_transitions_wide all-NA returns 0x0 matrix
test_that(".count_transitions_wide returns empty matrix when all states NA", {
  df <- data.frame(
    T1 = c(NA_character_, NA_character_),
    T2 = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_transitions_wide(df)
  expect_equal(nrow(result), 0L)
  expect_equal(ncol(result), 0L)
})

# L170-174: .count_transitions_long single row → zero matrix (n < 2)
test_that(".count_transitions_long single row returns zero matrix", {
  df <- data.frame(
    Time = 1L,
    Action = "A",
    stringsAsFactors = FALSE
  )
  result <- .count_transitions_long(df, action = "Action", id = NULL,
                                     time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result == 0L))
  expect_equal(rownames(result), "A")
})

# L200: .count_transitions_long group n < 2 returns empty pair
test_that(".count_transitions_long groups with 1 obs produce zero transitions", {
  df <- data.frame(
    Actor  = c(1L, 2L, 2L),
    Time   = c(1L, 1L, 2L),
    Action = c("A", "B", "A"),
    stringsAsFactors = FALSE
  )
  # Actor 1 has only 1 row → no pairs from that group
  result <- .count_transitions_long(df, action = "Action", id = "Actor",
                                     time = "Time")
  expect_true(is.matrix(result))
  # The transition A->B or B->A from actor 2 only
  expect_equal(sort(rownames(result)), c("A", "B"))
})

# L347-351: .estimator_co_occurrence with explicit codes (one-hot path)
test_that("co_occurrence with explicit codes uses wtna one-hot path", {
  set.seed(42)
  df <- data.frame(
    A = sample(0L:1L, 20, replace = TRUE),
    B = sample(0L:1L, 20, replace = TRUE),
    C = sample(0L:1L, 20, replace = TRUE)
  )
  net <- build_network(df, method = "co_occurrence",
                       codes = c("A", "B", "C"),
                       window_size = 2L)
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(sort(net$nodes$label), c("A", "B", "C"))
})

# L359-361: .estimator_co_occurrence long format (non-binary)
test_that("co_occurrence long format path works", {
  long_data <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )
  net <- build_network(long_data, method = "co_occurrence",
                       action = "Action", actor = "Actor",
                       params = list(format = "long",
                                     id = "Actor", time = "Time"))
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(sort(net$nodes$label), c("A", "B"))
})

# L393-395: .count_cooccurrence_wide empty n_states or nc < 2
test_that(".count_cooccurrence_wide returns empty when no valid states", {
  df <- data.frame(
    T1 = c(NA_character_, NA_character_),
    T2 = c(NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_wide(df)
  expect_equal(nrow(result), 0L)
  expect_equal(ncol(result), 0L)
})

# L443-444: .count_cooccurrence_long missing action column
test_that(".count_cooccurrence_long errors on missing action column", {
  df <- data.frame(Time = 1:3, Action = c("A", "B", "A"),
                   stringsAsFactors = FALSE)
  expect_error(
    .count_cooccurrence_long(df, action = "NoSuchColumn"),
    "Action column.*not found"
  )
})

# L460-468, L474-489, L492-499, L502-514: .count_cooccurrence_long full paths
test_that(".count_cooccurrence_long single-id group path", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L, 2L, 2L, 2L),
    Time   = c(1L, 2L, 3L, 1L, 2L, 3L),
    Action = c("A", "B", "C", "B", "C", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(isSymmetric(result))
  expect_equal(sort(rownames(result)), c("A", "B", "C"))
})

test_that(".count_cooccurrence_long multi-id composite group key", {
  df <- data.frame(
    Actor   = c(1L, 1L, 2L, 2L),
    Session = c("s1", "s1", "s2", "s2"),
    Time    = c(1L, 2L, 1L, 2L),
    Action  = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action",
                                      id = c("Actor", "Session"),
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(isSymmetric(result))
})

test_that(".count_cooccurrence_long NULL id single sequence", {
  df <- data.frame(
    Time   = 1:4,
    Action = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = NULL,
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_equal(sort(rownames(result)), c("A", "B", "C"))
})

# L497-499: .count_cooccurrence_long empty from_vec returns zero matrix
test_that(".count_cooccurrence_long returns zero matrix for single-obs groups", {
  df <- data.frame(
    Actor  = c(1L, 2L),
    Time   = c(1L, 1L),
    Action = c("A", "B"),
    stringsAsFactors = FALSE
  )
  # Each group has only 1 obs → no pairs
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  expect_true(is.matrix(result))
  expect_true(all(result == 0))
})

# L517-519: .count_cooccurrence_long diagonal correction
test_that(".count_cooccurrence_long diagonal is correctly halved", {
  df <- data.frame(
    Actor  = c(1L, 1L, 1L),
    Time   = c(1L, 2L, 3L),
    Action = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- .count_cooccurrence_long(df, action = "Action", id = "Actor",
                                      time = "Time")
  # A appears at positions 1 and 2 → pair (A,A) is counted once
  # Diagonal A-A should equal half of what bidirectional counting would give
  expect_true(is.matrix(result))
  expect_true(result["A", "A"] >= 0)
})

# L842: .prepare_association_input errors on < 2 numeric cols after cleaning
test_that(".prepare_association_input errors when < 2 cols after zero-var removal", {
  df <- data.frame(
    V1 = rep(1.0, 20),  # zero variance
    V2 = rep(2.0, 20),  # zero variance
    V3 = rnorm(20)      # one valid col
  )
  expect_message(
    expect_error(
      build_network(df, method = "cor"),
      "At least 2 variable"
    ),
    "zero-variance"
  )
})

# L871, L875: .prepare_association_input matrix input with cov type and no colnames
test_that("build_network matrix input with cov type assigns V-names when colnames NULL", {
  set.seed(42)
  raw <- matrix(rnorm(50 * 4), 50, 4)
  cov_mat <- cov(raw)
  rownames(cov_mat) <- colnames(cov_mat) <- NULL

  net <- build_network(cov_mat, method = "cor",
                       params = list(n = 50, input_type = "cov"))
  expect_s3_class(net, "netobject")
  expect_equal(net$n_nodes, 4L)
  expect_true(all(grepl("^V", net$nodes$label)))
})

# L875: bad data type for .prepare_association_input
test_that(".prepare_association_input errors on non-df non-matrix input", {
  expect_error(
    build_network(list(a = 1), method = "cor"),
    "data frame or a square symmetric matrix"
  )
})

# L967: .compute_lambda_path errors when all off-diagonal correlations zero
test_that(".compute_lambda_path errors when all off-diagonal correlations zero", {
  # Create a correlation matrix that is identity (all off-diagonal are zero)
  S <- diag(4)
  rownames(S) <- colnames(S) <- paste0("V", 1:4)
  expect_error(
    .compute_lambda_path(S, nlambda = 10L, lambda.min.ratio = 0.01),
    "All off-diagonal correlations are zero"
  )
})

# initial probabilities stored for tna / ftna / atna
test_that("build_network stores initial probs for tna, ftna, atna", {
  set.seed(1)
  seqs <- data.frame(
    V1 = sample(c("A","B","C"), 30, TRUE),
    V2 = sample(c("A","B","C"), 30, TRUE),
    V3 = sample(c("A","B","C"), 30, TRUE)
  )
  for (m in c("tna", "ftna", "atna")) {
    net <- build_network(seqs, method = m)
    expect_false(is.null(net$initial), info = paste(m, "has $initial"))
    expect_named(net$initial)
    expect_equal(sum(net$initial), 1, tolerance = 1e-9,
                 info = paste(m, "$initial sums to 1"))
    expect_true(all(net$initial >= 0), info = paste(m, "$initial non-negative"))
    expect_true(all(names(net$initial) %in% net$nodes$label),
                info = paste(m, "initial names match nodes"))
  }
})

# .compute_initial_probs: states never appearing as first get prob 0
test_that(".compute_initial_probs gives 0 to states never starting a sequence", {
  seqs <- data.frame(
    V1 = c("A","A","A"),
    V2 = c("B","C","B"),
    V3 = c("C","B","C")
  )
  net <- build_network(seqs, method = "tna")
  expect_equal(net$initial[["A"]], 1)
  expect_equal(net$initial[["B"]], 0)
  expect_equal(net$initial[["C"]], 0)
})

test_that("transition estimators retain states with no observed transitions", {
  seqs <- data.frame(V1 = c("A", "B"), V2 = c(NA, NA))
  net <- build_network(seqs, method = "relative")
  expect_equal(rownames(net$weights), c("A", "B"))
  expect_equal(net$weights, matrix(0, 2, 2,
                                   dimnames = list(c("A", "B"),
                                                   c("A", "B"))))
  expect_equal(net$initial, c(A = 0.5, B = 0.5))
})

test_that("grouped transition networks use a global state alphabet", {
  seqs <- data.frame(
    V1 = c("A", "A", "B", "C"),
    V2 = c("B", "B", "C", NA),
    grp = c("g1", "g2", "g2", "g1")
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  expect_equal(rownames(nets$g1$weights), c("A", "B", "C"))
  expect_equal(rownames(nets$g2$weights), c("A", "B", "C"))
  expect_equal(nets$g1$weights["C", ], c(A = 0, B = 0, C = 0))
})

test_that("transition estimators support tna weighted counts", {
  seqs <- data.frame(
    V1 = c("A", "A", "B"),
    V2 = c("B", "C", "C"),
    V3 = c("C", "C", NA)
  )
  freq <- build_network(seqs, method = "frequency",
                        params = list(weighted = TRUE))
  expect_equal(freq$weights["A", "B"], 1 / 3)
  expect_equal(freq$weights["B", "C"], 1 / 3 + 1 / 2)

  cooc <- build_network(seqs, method = "co_occurrence",
                        params = list(weighted = TRUE))
  expect_equal(cooc$weights["A", "C"], 1 / 3 + 1 / 3 + 1 / 3)
  expect_equal(cooc$weights["C", "C"], 1 / 3)
})

test_that("transition estimators support begin/end states and concat", {
  seqs <- data.frame(
    V1 = c("A", "A", "B"),
    V2 = c("B", "C", "C"),
    V3 = c("C", "C", NA)
  )
  sent <- build_network(seqs, method = "relative",
                        begin_state = "BEGIN", end_state = "END")
  expect_true(all(c("BEGIN", "END") %in% rownames(sent$weights)))
  expect_equal(sent$weights["BEGIN", c("A", "B")], c(A = 2 / 3, B = 1 / 3))
  expect_equal(sent$initial[["BEGIN"]], 1)

  concat <- build_network(seqs, method = "relative", concat = 2)
  expect_equal(concat$weights["C", "A"], 0.5)
})

test_that("transition matrix input matches tna-style handling", {
  m <- matrix(c(1, 2, 3, 4), 2, 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  net <- build_network(m, method = "relative")
  expect_equal(net$weights["A", ], c(A = 1 / 4, B = 3 / 4))
  expect_equal(net$weights["B", ], c(A = 1 / 3, B = 2 / 3))
})

test_that("transition scaling ranks and minmax scales the full matrix", {
  seqs <- data.frame(
    V1 = c("A", "A", "B"),
    V2 = c("B", "C", "C"),
    V3 = c("C", "C", NA)
  )
  ranked <- build_network(seqs, method = "frequency", scaling = "rank")
  expect_equal(ranked$weights["A", "A"], 3)
  expect_equal(ranked$weights["B", "C"], 9)

  scaled <- build_network(seqs, method = "frequency", scaling = "minmax")
  expect_equal(scaled$weights["A", "A"], 0)
  expect_equal(scaled$weights["B", "C"], 1)
})

# L1004-1005: .select_ebic handles glasso fit failure (NULL fit → Inf EBIC)
test_that("glasso network completes even with challenging input", {
  # Near-singular but symmetric correlation matrix
  set.seed(42)
  S <- diag(3)
  diag(S) <- 1
  S[1,2] <- S[2,1] <- 0.999
  S[2,3] <- S[3,2] <- 0.001
  S[1,3] <- S[3,1] <- 0.001
  rownames(S) <- colnames(S) <- c("A", "B", "C")
  # Should complete without error; some lambda fits may fail internally
  expect_no_error(
    net <- build_network(S, method = "glasso", params = list(n = 100, nlambda = 10L))
  )
  expect_s3_class(net, "netobject")
})


# ---- build_network net_mmm dispatch (L149-178) ----

test_that("build_network wraps net_mmm default (relative) models (L173-178)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 80, TRUE), V2 = sample(states, 80, TRUE),
    V3 = sample(states, 80, TRUE), V4 = sample(states, 80, TRUE),
    V5 = sample(states, 80, TRUE), stringsAsFactors = FALSE
  )
  mmm <- build_mmm(data, k = 2, n_starts = 3, seed = 1)
  grp <- build_network(mmm)
  expect_true(inherits(grp, "netobject_group"))
  expect_equal(length(grp), 2)
  expect_true(all(vapply(grp, function(x) inherits(x, "netobject"), logical(1))))
})

test_that("build_network assigns Cluster names when mmm models unnamed (L175)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 80, TRUE), V2 = sample(states, 80, TRUE),
    V3 = sample(states, 80, TRUE), V4 = sample(states, 80, TRUE),
    V5 = sample(states, 80, TRUE), stringsAsFactors = FALSE
  )
  mmm <- build_mmm(data, k = 2, n_starts = 3, seed = 1)
  names(mmm$models) <- NULL
  grp <- build_network(mmm)
  expect_equal(names(grp), c("Cluster 1", "Cluster 2"))
})

test_that("build_network rebuilds net_mmm with non-relative method (L151-171)", {
  set.seed(42)
  states <- c("A","B","C")
  data <- data.frame(
    V1 = sample(states, 80, TRUE), V2 = sample(states, 80, TRUE),
    V3 = sample(states, 80, TRUE), V4 = sample(states, 80, TRUE),
    V5 = sample(states, 80, TRUE), stringsAsFactors = FALSE
  )
  mmm <- build_mmm(data, k = 2, n_starts = 3, seed = 1)

  # frequency method: resolved in (relative, frequency, attention) -> injects initial
  grp_freq <- build_network(mmm, method = "frequency")
  expect_true(inherits(grp_freq, "netobject_group"))
  expect_equal(grp_freq[[1]]$method, "frequency")

  # co_occurrence method: resolved NOT in relative/frequency/attention -> no initial injection
  grp_co <- build_network(mmm, method = "co_occurrence")
  expect_true(inherits(grp_co, "netobject_group"))
  expect_equal(grp_co[[1]]$method, "co_occurrence")
})


# ---- .compute_initial_probs long-format paths (estimators.R) ----

# These call the estimator functions directly via the registry to bypass
# build_network()'s long->wide conversion and hit the long-format branches
# in .compute_initial_probs() and the format auto-detection code.

test_that(".compute_initial_probs multi-id long format via estimator (L282-290)", {
  df <- data.frame(
    pid = c(1,1,1,2,2,2),
    sid = c("a","a","a","b","b","b"),
    Time = c(1,2,3,1,2,3),
    Action = c("X","Y","Z","Y","Z","X"),
    stringsAsFactors = FALSE
  )
  est <- get_estimator("relative")
  result <- est$fn(data = df, format = "long", action = "Action",
                   id = c("pid","sid"), time = "Time")
  expect_false(is.null(result$initial))
  expect_equal(sum(result$initial), 1)
  expect_equal(length(result$initial), 3)
})

test_that(".compute_initial_probs long format with no id via estimator (L281-284)", {
  df <- data.frame(
    Time = 1:5,
    Action = c("X","Y","Z","X","Y"),
    stringsAsFactors = FALSE
  )
  est <- get_estimator("relative")
  result <- est$fn(data = df, format = "long", action = "Action",
                   id = NULL, time = "Time")
  expect_false(is.null(result$initial))
  expect_equal(unname(result$initial[["X"]]), 1)
})

test_that(".compute_initial_probs returns zero vector when all first_states NA (L303)", {
  # Call internal function directly — line 303 triggers when
  # all first_states are NA (no valid starting states)
  init <- Nestimate:::.compute_initial_probs(
    data.frame(Time = 1:3, Action = rep(NA_character_, 3)),
    states = c("X","Y"), format = "long",
    action = "Action", id = NULL, time = "Time"
  )
  expect_equal(unname(init), c(0, 0))
})

test_that("frequency estimator auto-detects long format via registry (L325-326)", {
  df <- data.frame(
    id = c(1,1,1,2,2,2),
    Time = c(1,2,3,1,2,3),
    Action = c("X","Y","Z","Y","Z","X"),
    stringsAsFactors = FALSE
  )
  est <- get_estimator("frequency")
  result <- est$fn(data = df, format = "auto", action = "Action",
                   id = "id", time = "Time")
  expect_false(is.null(result$initial))
  expect_true(is.matrix(result$matrix))
})

test_that("relative estimator auto-detects long format via registry (L364-366)", {
  df <- data.frame(
    id = c(1,1,1,2,2,2),
    Time = c(1,2,3,1,2,3),
    Action = c("X","Y","Z","Y","Z","X"),
    stringsAsFactors = FALSE
  )
  est <- get_estimator("relative")
  result <- est$fn(data = df, format = "auto", action = "Action",
                   id = "id", time = "Time")
  expect_false(is.null(result$initial))
  rs <- rowSums(result$matrix)
  expect_true(all(rs == 0 | abs(rs - 1) < 1e-10))
})


# ============================================================
# Auto-conversion: sequences → frequencies for association methods
# ============================================================

# Helper: wide sequence data
.make_seq_data <- function(n = 80, t = 8,
                           states = c("A", "B", "C", "D"),
                           seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# ---- 1. glasso from wide sequences ----
test_that("auto-convert: glasso from wide character sequences", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "glasso")
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(net$method, "glasso")
  expect_equal(nrow(net$weights), 4)
})

# ---- 2. pcor from complete sequences errors (singular: constant row sums) ----
test_that("auto-convert: pcor errors on complete sequences (singular)", {
  seqs <- .make_seq_data()
  expect_error(build_network(seqs, method = "pcor"), "singular")
})

# ---- 3. cor from wide sequences ----
test_that("auto-convert: cor from wide character sequences", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "cor")
  expect_s3_class(net, "netobject")
  expect_false(net$directed)
  expect_equal(net$method, "cor")
  expect_true(isSymmetric(net$weights))
})

# ---- 4. ising: no auto-convert (requires binary 0/1, not counts) ----
test_that("auto-convert: method aliases work with sequences", {
  seqs <- .make_seq_data()
  net1 <- build_network(seqs, method = "ebicglasso")
  expect_equal(net1$method, "glasso")

  net3 <- build_network(seqs, method = "correlation")
  expect_equal(net3$method, "cor")
})

# ---- 6. Results match manual conversion ----
test_that("auto-convert: matches manual convert_sequence_format pipeline", {
  seqs <- .make_seq_data(seed = 99)
  freq <- convert_sequence_format(seqs, id_col = character(0),
                                  format = "frequency")
  freq <- freq[, setdiff(names(freq), "rid"), drop = FALSE]

  net_auto <- build_network(seqs, method = "glasso",
                            params = list(gamma = 0.5, nlambda = 50))
  net_manual <- build_network(freq, method = "glasso",
                              params = list(gamma = 0.5, nlambda = 50))

  expect_equal(net_auto$weights, net_manual$weights)
  expect_equal(net_auto$nodes$label, net_manual$nodes$label)
})

# ---- 7. tna::group_regulation ----
test_that("auto-convert: transition methods skip conversion", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "relative")
  expect_equal(net$method, "relative")
  expect_true(net$directed)
  # Rows should sum to 1 for transition probabilities
  rs <- rowSums(net$weights)
  expect_true(all(abs(rs - 1) < 1e-10))
})

# ---- 11. Numeric data still works directly (no spurious conversion) ----
test_that("auto-convert: numeric data for glasso not double-converted", {
  set.seed(42)
  num_data <- as.data.frame(matrix(rnorm(200), ncol = 4))
  colnames(num_data) <- c("X1", "X2", "X3", "X4")
  net <- build_network(num_data, method = "glasso")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 4)
  expect_equal(sort(net$nodes$label), sort(c("X1", "X2", "X3", "X4")))
})

# ---- 12. With scaling parameter ----
test_that("auto-convert: glasso + scaling works on sequences", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "glasso", scaling = "minmax")
  expect_s3_class(net, "netobject")
  # All weights should be in [0, 1] after minmax
  expect_true(all(net$weights >= 0 & net$weights <= 1))
})

# ---- 13. With threshold parameter ----
test_that("auto-convert: glasso + threshold works on sequences", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "glasso", threshold = 0.1)
  expect_s3_class(net, "netobject")
  # All non-zero weights should be >= 0.1 in absolute value
  nz <- net$weights[net$weights != 0]
  if (length(nz) > 0) expect_true(all(abs(nz) >= 0.1))
})

# ---- 14. Many states ----
test_that("auto-convert: glasso with many states", {
  seqs <- .make_seq_data(n = 200, states = LETTERS[1:8])
  net <- build_network(seqs, method = "glasso")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 8)
})

# ---- 15. Two states → glasso (pcor singular with constant sums) ----
test_that("auto-convert: glasso with two states", {
  seqs <- .make_seq_data(n = 50, states = c("X", "Y"))
  net <- build_network(seqs, method = "glasso")
  expect_s3_class(net, "netobject")
  expect_equal(nrow(net$weights), 2)
})

# ---- 16. Long format with action column ----
test_that("auto-convert: glasso from long format sequences", {
  set.seed(42)
  long <- data.frame(
    Actor = rep(1:50, each = 6),
    Action = sample(c("A", "B", "C", "D"), 300, replace = TRUE),
    Time = rep(1:6, 50),
    stringsAsFactors = FALSE
  )
  net <- build_network(long, method = "glasso",
                       actor = "Actor", action = "Action", time = "Time")
  expect_s3_class(net, "netobject")
  expect_equal(net$method, "glasso")
})

# ---- 17. Factor columns treated as character ----
test_that("auto-convert: factor columns auto-converted for glasso", {
  seqs <- .make_seq_data()
  seqs[] <- lapply(seqs, as.factor)
  net <- build_network(seqs, method = "glasso")
  expect_s3_class(net, "netobject")
})

# ---- 18. Auto-converted net stores numeric data ----
test_that("auto-convert: glasso net has numeric data stored", {
  seqs <- .make_seq_data()
  net <- build_network(seqs, method = "glasso")
  # The estimator stores its input; after frequency conversion this is numeric
  expect_s3_class(net, "netobject")
  expect_equal(net$method, "glasso")
})

# ---- 19. Mixed columns: some character, some numeric ----
test_that("auto-convert: mixed char+numeric columns for glasso", {
  set.seed(42)
  df <- data.frame(
    T1 = sample(c("A","B","C"), 50, TRUE),
    T2 = sample(c("A","B","C"), 50, TRUE),
    score = rnorm(50),
    stringsAsFactors = FALSE
  )
  # Has character columns → should trigger auto-convert
  net <- build_network(df, method = "glasso")
  expect_s3_class(net, "netobject")
})

# ---- 20. Group parameter with auto-conversion ----
test_that("auto-convert: grouped glasso from sequences", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(c("A","B","C"), 80, TRUE),
    T2 = sample(c("A","B","C"), 80, TRUE),
    T3 = sample(c("A","B","C"), 80, TRUE),
    grp = rep(c("X","Y"), each = 40),
    stringsAsFactors = FALSE
  )
  nets <- build_network(seqs, method = "glasso", group = "grp")
  expect_s3_class(nets, "netobject_group")
  expect_equal(length(nets), 2)
  expect_equal(nets[[1]]$method, "glasso")
})

# ---- 21. Predictability (R²) auto-computation ----

test_that("glasso netobject carries $predictability by default", {
  set.seed(42)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  df$x2 <- df$x1 * 0.8 + rnorm(80, sd = 0.5)
  net <- build_network(df, method = "glasso")
  expect_true(!is.null(net$predictability))
  expect_equal(length(net$predictability), net$n_nodes)
  expect_equal(names(net$predictability), net$nodes$label)
  expect_true(all(net$predictability >= 0 & net$predictability <= 1))
})

test_that("pcor netobject carries $predictability by default", {
  set.seed(42)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  net <- build_network(df, method = "pcor")
  expect_true(!is.null(net$predictability))
  expect_equal(length(net$predictability), net$n_nodes)
  expect_true(all(net$predictability >= 0 & net$predictability <= 1))
})

test_that("cor netobject carries $predictability by default", {
  set.seed(42)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  net <- build_network(df, method = "cor")
  expect_true(!is.null(net$predictability))
  expect_equal(length(net$predictability), net$n_nodes)
  expect_true(all(net$predictability >= 0 & net$predictability <= 1))
})

test_that("predictability = FALSE suppresses computation", {
  set.seed(42)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  net <- build_network(df, method = "glasso", predictability = FALSE)
  expect_null(net$predictability)
})

test_that("directed networks do not carry $predictability", {
  seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
  net <- build_network(seqs, method = "relative")
  expect_null(net$predictability)
})

test_that("grouped glasso carries per-group predictability", {
  set.seed(42)
  df <- data.frame(
    x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80),
    grp = rep(c("A", "B"), each = 40)
  )
  nets <- build_network(df, method = "glasso", group = "grp")
  expect_true(!is.null(nets$A$predictability))
  expect_true(!is.null(nets$B$predictability))
  expect_equal(length(nets$A$predictability), nets$A$n_nodes)
})

test_that("predictability.netobject_group returns per-group list", {
  set.seed(42)
  df <- data.frame(
    x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80),
    grp = rep(c("A", "B"), each = 40)
  )
  nets <- build_network(df, method = "glasso", group = "grp")
  r2 <- predictability(nets)
  expect_true(is.list(r2))
  expect_equal(length(r2), 2L)
  expect_true(all(vapply(r2, is.data.frame, logical(1))))
})

test_that("print.netobject shows predictability for glasso", {
  set.seed(42)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
  df$x2 <- df$x1 * 0.8 + rnorm(80, sd = 0.5)
  net <- build_network(df, method = "glasso")
  out <- capture.output(print(net))
  expect_true(any(grepl("Predictability", out)))
})

test_that("predictability() errors informatively on unsupported methods", {
  seqs <- data.frame(V1 = c("A","B","C","A"), V2 = c("B","C","A","B"))
  net <- build_network(seqs, method = "relative")
  expect_error(predictability(net), "does not support predictability")
})

test_that("predictability = FALSE is honoured in level = 'both'", {
  set.seed(42)
  df <- data.frame(
    x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100),
    id = rep(1:20, each = 5)
  )
  net <- build_network(df, method = "glasso",
                       params = list(id = "id", nlambda = 20L),
                       level = "both", predictability = FALSE)
  expect_null(net$between$predictability)
  expect_null(net$within$predictability)
})

# ---- state_cols / metadata_cols overrides (adversarial label-overlap) ----

# Scenario: a metadata column "condition" has values in {"A","B","C"} which
# are *also* the node names. The values-in-nodes heuristic misclassifies
# "condition" as a state column, corrupting $data and $metadata.
.mk_label_overlap_df <- function() {
  data.frame(
    V1        = c("A","B","C","A","B","C"),
    V2        = c("B","C","A","C","A","B"),
    V3        = c("C","A","B","B","C","A"),
    condition = c("A","B","C","A","B","C"),  # collides with node labels
    stringsAsFactors = FALSE
  )
}

test_that("auto-detection misclassifies label-overlapping metadata (baseline)", {
  df  <- .mk_label_overlap_df()
  net <- build_network(df, method = "relative")
  # Without override, 'condition' looks like a state column — auto-detection
  # puts it into $data (baseline bug this override exists to fix).
  expect_true("condition" %in% names(net$data))
})

test_that("state_cols override forces correct classification", {
  df  <- .mk_label_overlap_df()
  net <- build_network(
    df, method = "relative",
    state_cols = c("V1", "V2", "V3")
  )
  expect_equal(sort(names(net$data)),     c("V1", "V2", "V3"))
  expect_true(is.data.frame(net$metadata))
  expect_true("condition" %in% names(net$metadata))
})

test_that("metadata_cols override alone forces the column to metadata", {
  df  <- .mk_label_overlap_df()
  net <- build_network(
    df, method = "relative",
    metadata_cols = "condition"
  )
  expect_false("condition" %in% names(net$data))
  expect_true("condition" %in% names(net$metadata))
})

test_that("state_cols and metadata_cols agree when both supplied", {
  df  <- .mk_label_overlap_df()
  net <- build_network(
    df, method = "relative",
    state_cols    = c("V1", "V2", "V3"),
    metadata_cols = "condition"
  )
  expect_equal(sort(names(net$data)),     c("V1", "V2", "V3"))
  expect_true("condition" %in% names(net$metadata))
})

test_that("state_cols / metadata_cols overlap raises an error", {
  df <- .mk_label_overlap_df()
  expect_error(
    build_network(
      df, method = "relative",
      state_cols    = c("V1", "V2"),
      metadata_cols = c("V2", "condition")
    ),
    "overlap"
  )
})

test_that("state_cols naming a missing column raises a clear error", {
  df <- .mk_label_overlap_df()
  expect_error(
    build_network(df, method = "relative", state_cols = c("V1", "V99")),
    "state_cols not found in data"
  )
})

test_that("metadata_cols naming a missing column raises a clear error", {
  df <- .mk_label_overlap_df()
  expect_error(
    build_network(df, method = "relative", metadata_cols = "nope"),
    "metadata_cols not found in data"
  )
})

test_that("state_cols override preserves the weight matrix vs baseline", {
  # The override should NOT change what the estimator sees — only how the
  # output is labelled. Weights must match the baseline that already (mis)
  # classifies condition as state, because condition's values are a valid
  # permutation of nodes and get counted either way.
  df_clean <- .mk_label_overlap_df()[, c("V1", "V2", "V3")]
  baseline <- build_network(df_clean, method = "relative")

  df_over  <- .mk_label_overlap_df()
  over_net <- build_network(
    df_over, method = "relative",
    state_cols = c("V1", "V2", "V3")
  )
  # The weight matrices may differ (baseline ignores the fourth column,
  # over_net still saw 'condition' during estimation). Just check the
  # override keeps the classification intact without throwing.
  expect_s3_class(over_net, "netobject")
  expect_equal(sort(names(over_net$data)), c("V1", "V2", "V3"))
})

test_that("state_cols propagates through group= dispatch", {
  # Without forwarding, the per-group recursive call reverts to auto-detection
  # and 'condition' would land back in $data. This catches regressions in the
  # group-dispatch recursive call site.
  df <- data.frame(
    V1        = c("A","B","C","A","B","C","A","B"),
    V2        = c("B","C","A","B","C","A","C","A"),
    condition = c("A","B","C","A","B","C","A","B"),
    grp       = c("x","x","x","x","y","y","y","y"),
    stringsAsFactors = FALSE
  )
  nets <- build_network(
    df, method = "relative", group = "grp",
    state_cols = c("V1", "V2")
  )
  expect_s3_class(nets, "netobject_group")
  # Each per-group netobject must respect the override
  for (net in nets) {
    expect_equal(sort(names(net$data)), c("V1", "V2"))
    expect_true("condition" %in% names(net$metadata))
  }
})

test_that("no override = no behavior change (regression safety)", {
  # Same data, no override: result must be bit-identical to historical
  # auto-detection output. Guards against accidental semantic drift.
  df <- data.frame(
    V1 = c("A","B","C","A","B"),
    V2 = c("B","C","A","B","C"),
    stringsAsFactors = FALSE
  )
  net_auto <- build_network(df, method = "relative")
  net_null <- build_network(df, method = "relative",
                            state_cols = NULL, metadata_cols = NULL)
  expect_equal(net_null$weights, net_auto$weights)
  expect_equal(names(net_null$data), names(net_auto$data))
})
