testthat::skip_on_cran()

# Regression tests for FIX2 (audit findings A02-F01, A02-F02, A09-F01,
# A09-F02, A09-F03, A09-F04). Each test reproduces the documented bug
# and asserts the corrected behavior.

# Helper: strong-correlation continuous data for a saturated glasso net.
.fix2_cont <- function(n = 250, p = 6, rho = 0.45, seed = 55) {
  set.seed(seed)
  z <- rnorm(n)
  mat <- vapply(seq_len(p), function(j) {
    sqrt(rho) * z + sqrt(1 - rho) * rnorm(n)
  }, numeric(n))
  colnames(mat) <- paste0("V", seq_len(p))
  as.data.frame(mat)
}

# Helper: directed transition net from bundled group_regulation_long.
.fix2_directed_net <- function(method = "relative") {
  e <- new.env()
  data("group_regulation_long", package = "Nestimate", envir = e)
  suppressWarnings(build_network(e$group_regulation_long, method = method))
}


# ---- A02-F01: boot_glasso(centrality_fn=) is now reachable ----

test_that("A02-F01: boot_glasso invokes a user-supplied centrality_fn", {
  cd <- .fix2_cont()
  hit <- new.env(); hit$n <- 0L
  my_fn <- function(mat) {
    hit$n <- hit$n + 1L
    list(my_metric = setNames(seq_len(ncol(mat)) * 100, colnames(mat)))
  }

  b <- boot_glasso(cd, iter = 8, cs_iter = 3, seed = 2,
                   centrality = c("strength", "my_metric"),
                   centrality_fn = my_fn)

  # The custom measure is now produced and the fn was actually called.
  expect_true("my_metric" %in% names(b$original_centrality))
  expect_gt(hit$n, 0L)
  expect_equal(unname(b$original_centrality$my_metric),
               c(100, 200, 300, 400, 500, 600))
  expect_true(all(c("strength", "my_metric") %in% b$centrality_measures))
  expect_true("my_metric" %in% names(b$cs_coefficient))
})

test_that("A02-F01: builtins still use builtin path even with centrality_fn", {
  cd <- .fix2_cont()
  hit <- new.env(); hit$n <- 0L
  fn <- function(mat) {
    hit$n <- hit$n + 1L
    list(closeness = setNames(rep(7, ncol(mat)), colnames(mat)))
  }
  # Pure builtin request + a centrality_fn: builtin wins, fn NOT called.
  b <- boot_glasso(cd, iter = 8, cs_iter = 3, seed = 2,
                   centrality = "closeness", centrality_fn = fn)
  expect_false(all(b$original_centrality$closeness == 7))
  expect_equal(hit$n, 0L)
})

test_that("A02-F01: garbage centrality name still rejected when no fn", {
  cd <- .fix2_cont()
  expect_error(
    boot_glasso(cd, iter = 5, cs_iter = 2, centrality = "garbage"),
    "should be one of"
  )
})


# ---- A02-F02: roxygen now matches behavior (betweenness/closeness
#       compute with centrality_fn = NULL) ----

test_that("A02-F02: betweenness & closeness compute without centrality_fn", {
  cd <- .fix2_cont()
  b <- boot_glasso(cd, iter = 8, cs_iter = 3, seed = 2,
                   centrality = c("betweenness", "closeness"))
  expect_false(is.null(b$original_centrality$betweenness))
  expect_false(is.null(b$original_centrality$closeness))
  expect_true("betweenness" %in% names(b$centrality_ci))
  expect_true("closeness" %in% names(b$centrality_ci))
})


# ---- A09-F01: centrality_stability(centrality_fn=) is now reachable ----

test_that("A09-F01: centrality_stability invokes a user-supplied fn", {
  netd <- .fix2_directed_net()
  flag <- new.env(); flag$hit <- FALSE
  fn <- function(m) {
    flag$hit <- TRUE
    list(Eigen = setNames(abs(rnorm(nrow(m))) + 1, rownames(m)))
  }
  cs <- centrality_stability(netd, measures = "Eigen", centrality_fn = fn,
                             iter = 8, drop_prop = c(0.2, 0.4), seed = 1)
  expect_true(flag$hit)
  expect_s3_class(cs, "net_stability")
  expect_true("Eigen" %in% names(cs$cs))
})

test_that("A09-F01: builtin measure does not call centrality_fn", {
  netd <- .fix2_directed_net()
  flag <- new.env(); flag$hit <- FALSE
  fn <- function(m) {
    flag$hit <- TRUE
    list(Eigen = setNames(rnorm(nrow(m)), rownames(m)))
  }
  invisible(centrality_stability(netd, measures = "Betweenness",
                                 centrality_fn = fn, iter = 8,
                                 drop_prop = c(0.2, 0.4), seed = 1))
  expect_false(flag$hit)
})

test_that("A09-F01: unknown measure still rejected when no centrality_fn", {
  netd <- .fix2_directed_net()
  expect_error(
    centrality_stability(netd, measures = "Eigen", iter = 5,
                         drop_prop = 0.3),
    "Unknown measures"
  )
})


# ---- A09-F02: centrality_stability("Closeness") on directed errors
#       cleanly (no cryptic crash) ----

test_that("A09-F02: directed Closeness gives an actionable error", {
  netd <- .fix2_directed_net()
  expect_error(
    centrality_stability(netd, measures = "Closeness", iter = 8,
                         drop_prop = c(0.2, 0.5), seed = 1),
    "only for undirected networks"
  )
  # The default directed measures still work.
  cs <- centrality_stability(netd, iter = 6, drop_prop = c(0.2, 0.4),
                             seed = 1)
  expect_s3_class(cs, "net_stability")
  expect_true(all(c("InStrength", "OutStrength", "Betweenness") %in%
                    names(cs$cs)))
})


# ---- A09-F03: net_centrality("Closeness") on directed errors cleanly
#       (no silent 0-column data.frame) ----

test_that("A09-F03: directed net_centrality Closeness errors clearly", {
  netd <- .fix2_directed_net()
  expect_error(
    suppressMessages(net_centrality(netd, measures = "Closeness")),
    "only for undirected networks"
  )
  expect_error(
    suppressMessages(net_centrality(netd,
                                    measures = c("InStrength", "Closeness"))),
    "only for undirected networks"
  )
  # Undirected Closeness still works; defaults still work both ways.
  set.seed(1)
  panel <- as.data.frame(matrix(rnorm(600), 150, 4,
                                dimnames = list(NULL, LETTERS[1:4])))
  und <- build_network(panel, method = "cor")
  uc <- suppressMessages(net_centrality(und, measures = "Closeness"))
  expect_true("Closeness" %in% colnames(uc))
  expect_equal(nrow(uc), 4L)
  dd <- suppressMessages(net_centrality(netd))
  expect_true(all(c("InStrength", "OutStrength", "Betweenness") %in%
                    colnames(dd)))
})

test_that("A09-F03: undirected In/OutCloseness errors clearly", {
  set.seed(1)
  panel <- as.data.frame(matrix(rnorm(600), 150, 4,
                                dimnames = list(NULL, LETTERS[1:4])))
  und <- build_network(panel, method = "cor")
  expect_error(
    suppressMessages(net_centrality(und, measures = "InCloseness")),
    "only for directed networks"
  )
})


# ---- A09-F04: summary.net_casedrop_reliability_group() no longer
#       crashes when the grid omits 0.7 ----

test_that("A09-F04: group summary uses the object's grid, not 0.7", {
  netd <- .fix2_directed_net("relative")
  netf <- .fix2_directed_net("frequency")
  grp <- structure(list("Cluster 1" = netd, "Cluster 2" = netf),
                    class = c("netobject_group", "list"))

  rg <- suppressWarnings(
    casedrop_reliability(grp, iter = 8, drop_prop = c(0.2, 0.5), seed = 1)
  )

  # Default summary no longer crashes; it reports at a grid point.
  s <- summary(rg)
  expect_s3_class(s, "summary.net_casedrop_reliability_group")
  expect_equal(nrow(s), 2L)
  expect_true(any(abs(rg[[1]]$drop_prop - attr(s, "drop_prop")) < 1e-9))

  # Explicit in-grid value works.
  s2 <- summary(rg, drop_prop = 0.5)
  expect_equal(nrow(s2), 2L)

  # Explicit out-of-grid value errors cleanly (not a cryptic dimnames crash).
  expect_error(summary(rg, drop_prop = 0.7),
               "not in the object's grid")
})

test_that("A09-F04: documented example grid c(0.1,0.3,0.5) summarizes", {
  netd <- .fix2_directed_net("relative")
  netf <- .fix2_directed_net("frequency")
  grp <- structure(list("Cluster 1" = netd, "Cluster 2" = netf),
                    class = c("netobject_group", "list"))
  rg <- suppressWarnings(
    casedrop_reliability(grp, iter = 6, drop_prop = c(0.1, 0.3, 0.5),
                         seed = 1)
  )
  s <- summary(rg)
  expect_s3_class(s, "summary.net_casedrop_reliability_group")
  expect_true(any(abs(c(0.1, 0.3, 0.5) - attr(s, "drop_prop")) < 1e-9))
})
