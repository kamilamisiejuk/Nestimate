# Numerical equivalence: summary.netobject and compare_model() vs tna
#
# Gated by NESTIMATE_EQUIV_TESTS=true (see other test-equiv-* files).
# Runs only when the heavy reference packages are installed.
#
# Equivalence covers the 7 scaling methods that exist in BOTH packages
# without producing negative weights (which tna's igraph-backed network
# metrics reject). Nestimate-only scalings (`robust`, `frobenius`, `row`)
# are smoke-tested separately for stability, not equivalence.

skip_unless_equiv <- function() {
  if (!identical(tolower(Sys.getenv("NESTIMATE_EQUIV_TESTS")), "true")) {
    testthat::skip("NESTIMATE_EQUIV_TESTS != true")
  }
  testthat::skip_if_not_installed("tna")
  testthat::skip_if_not_installed("igraph")
}

testthat::test_that("summary.netobject matches tna::summary.tna across methods", {
  skip_unless_equiv()

  d <- tna::group_regulation
  cases <- list(
    list(method = "tna",  ref = function(x) tna::tna(x)),
    list(method = "ftna", ref = function(x) tna::ftna(x)),
    list(method = "atna", ref = function(x) tna::atna(x))
  )
  for (cs in cases) {
    nest <- Nestimate::build_network(d, method = cs$method)
    tmod <- cs$ref(d)
    s_n <- summary(nest)
    s_t <- as.data.frame(summary(tmod))
    m <- merge(s_n, s_t, by = "metric", suffixes = c("_n", "_t"))
    delta <- abs(m$value_n - m$value_t)
    testthat::expect_lt(max(delta, na.rm = TRUE), 1e-10,
                        label = paste("summary delta for", cs$method))
  }
})

testthat::test_that("compare_model() matches tna::compare() across shared scalings", {
  skip_unless_equiv()

  set.seed(20260509)
  # zscore produces negative weights and tna's network-metrics block hands
  # them to igraph::mean_distance which errors out — skip it for the
  # equivalence sweep. Tested separately as a Nestimate-only sanity check.
  shared <- c("none", "minmax", "max", "rank", "log1p", "softmax", "quantile")

  d <- tna::group_regulation
  n <- nrow(d)
  for (rep in seq_len(30)) {
    i1 <- sample.int(n, 200)
    i2 <- sample.int(n, 200)
    m1 <- Nestimate::build_network(d[i1, ], method = "tna")
    m2 <- Nestimate::build_network(d[i2, ], method = "tna")
    t1 <- tna::tna(d[i1, ]); t2 <- tna::tna(d[i2, ])
    sc <- sample(shared, 1)

    cn <- Nestimate::compare_model(m1, m2, scaling = sc)
    ct <- tna::compare(t1, t2, scaling = sc)

    sn <- cn$summary_metrics
    st <- as.data.frame(ct$summary_metrics)
    msum <- merge(sn, st, by = c("category", "metric"),
                  suffixes = c("_n", "_t"))
    testthat::expect_lt(
      max(abs(msum$value_n - msum$value_t), na.rm = TRUE),
      1e-10,
      label = paste("summary_metrics rep", rep, "scaling", sc)
    )

    nn <- cn$network_metrics
    tt <- as.data.frame(ct$network_metrics)
    mnet <- merge(nn, tt, by = "metric", suffixes = c("_n", "_t"))
    testthat::expect_lt(
      max(c(abs(mnet$x_n - mnet$x_t), abs(mnet$y_n - mnet$y_t)),
          na.rm = TRUE),
      1e-10,
      label = paste("network_metrics rep", rep, "scaling", sc)
    )
  }
})

testthat::test_that("Nestimate-only scalings produce finite, dimensioned output", {
  skip_unless_equiv()

  d <- tna::group_regulation
  m1 <- Nestimate::build_network(d[1:200, ], method = "tna")
  m2 <- Nestimate::build_network(d[1001:1200, ], method = "tna")

  for (sc in c("zscore", "robust", "frobenius", "row")) {
    cn <- Nestimate::compare_model(m1, m2, scaling = sc)
    testthat::expect_identical(dim(cn$matrices$x), dim(m1$weights),
                               label = paste("dims preserved for", sc))
    fin <- vapply(cn$summary_metrics$value, is.finite, logical(1L))
    testthat::expect_true(all(fin),
                          label = paste("finite metrics for", sc))
  }
})
