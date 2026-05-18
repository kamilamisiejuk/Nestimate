# Regression: mgm numeric-vs-categorical type detection must not be silent.
# The auto-detect logic is unchanged (mgm::mgm() machine-precision
# equivalence preserved); the *surprising* numeric->categorical
# reclassification (the <=10-distinct-integer rule, e.g. a Likert/count
# item) is now announced, and the explicit type=/level= override suppresses
# it. (Item 2 from the column-guess inventory.)

test_that("mgm warns, naming the column, when a numeric col is auto-typed categorical", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 120
  df <- data.frame(likert = sample(1:5, n, TRUE),
                   g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  expect_warning(build_network(df, method = "mgm"),
                 "auto-detected as CATEGORICAL.*likert")
})

test_that("mgm type warning does NOT fire for genuinely continuous data", {
  skip_if_not_installed("glmnet")
  set.seed(2); n <- 120
  dfg <- data.frame(a = rnorm(n), b = rnorm(n), c = rnorm(n))
  expect_no_warning(build_network(dfg, method = "mgm"))
})

test_that("explicit type=/level= suppresses the auto-detect warning (override reachable)", {
  skip_if_not_installed("glmnet")
  set.seed(3); n <- 120
  df <- data.frame(likert = sample(1:5, n, TRUE),
                   g1 = rnorm(n), g2 = rnorm(n), g3 = rnorm(n))
  expect_no_warning(
    suppressMessages(.estimator_mgm(df, type = c("c", "g", "g", "g"),
                                    level = c(5L, 1L, 1L, 1L)))
  )
})

test_that("warning is side-effect only: detection deterministic, weights unchanged", {
  skip_if_not_installed("glmnet")
  set.seed(4); n <- 110
  df <- data.frame(k = sample(1:4, n, TRUE), x = rnorm(n), y = rnorm(n))
  n1 <- suppressWarnings(build_network(df, method = "mgm"))
  n2 <- suppressWarnings(build_network(df, method = "mgm"))
  expect_equal(n1$weights, n2$weights)
  expect_s3_class(n1, "netobject")
})
