# Regression tests for FIX4 — A06 (mcml / mmm) confirmed findings.
#
# A06-F01: build_mcml(type="semi_markov") was byte-identical to type="tna"
#          (no holding-time model exists). Now rejected, not silently aliased.
# A06-F02: type="frequency" is a legitimate raw-count alias of type="raw" —
#          behavior kept, documented as an explicit alias.
# A06-F03: summary.net_mmm @return now documents the data.frame it really
#          returns (regression locks the actual return shape/class/visibility).
# A06-F04: net_aggregate_weights "density" must equal sum / n_possible (true
#          edge density), and differ from "mean" when zero edges are present.

test_that("A06-F01: build_mcml(type='semi_markov') errors instead of aliasing tna", {
  sq <- simulate_sequences(30, 5, 10, seed = 7)
  st <- sort(unique(unlist(sq)))
  cl <- list(G1 = st[1:2], G2 = st[3:5])

  # The honest behaviour: a documented-but-unimplemented construction errors
  # clearly rather than returning the tna numbers under a false label.
  expect_error(
    build_mcml(sq, cl, type = "semi_markov"),
    "semi_markov.*not implemented"
  )

  # Internal root cause also rejects it (was a shared if-|| branch with tna).
  raw <- matrix(c(0, 3, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(
    Nestimate:::.process_weights(raw, "semi_markov"),
    "'type' must be one of: ", fixed = TRUE
  )

  # Sanity: the real constructions still work and tna is still row-stochastic.
  a <- build_mcml(sq, cl, type = "tna")
  expect_equal(unname(rowSums(a$macro$weights)), rep(1, nrow(a$macro$weights)))
  expect_false("semi_markov" %in% eval(formals(build_mcml)$type))
})

test_that("A06-F02: type='frequency' is an explicit, behaviour-preserving alias of 'raw'", {
  sq <- simulate_sequences(30, 5, 10, seed = 3)
  st <- sort(unique(unlist(sq)))
  cl <- list(A = st[1:2], B = st[3:5])

  f <- build_mcml(sq, cl, type = "frequency")
  r <- build_mcml(sq, cl, type = "raw")

  # Documented as a synonym -> behaviour intentionally identical (raw counts).
  expect_equal(f$macro$weights, r$macro$weights)
  expect_equal(f$clusters$A$weights, r$clusters$A$weights)

  # Still a reachable enum value (alias kept, not dropped).
  expect_true("frequency" %in% eval(formals(build_mcml)$type))
})

test_that("A06-F03: summary.net_mmm returns the documented data.frame shapes", {
  sq <- simulate_sequences(30, 5, 10, seed = 2)
  m <- build_mmm(sq, k = 2, n_starts = 2, max_iter = 20, seed = 1)

  s <- summary(m)

  # No-covariate path: a plain data.frame with the documented columns,
  # returned VISIBLY, and NOT the input object (the old @return was wrong).
  expect_s3_class(s, "data.frame")
  expect_false("net_mmm" %in% class(s))
  expect_false(identical(s, m))
  expect_setequal(
    names(s),
    c("component", "prior", "n_assigned", "mean_posterior", "avepp")
  )
  expect_identical(nrow(s), 2L)
  expect_true(withVisible(summary(m))$visible)
})

test_that("A06-F04: density divides by n_possible and differs from mean with zeros", {
  # Documented contract: density = sum / number of possible edges.
  w <- c(0.5, 0.8, 0.3, 0.9)
  expect_equal(net_aggregate_weights(w, "density", n_possible = 9), 2.5 / 9)

  # With a zero edge present, density (÷ n_possible) must be STRICTLY smaller
  # than mean (÷ non-zero count). They were byte-identical before because the
  # bare call fell back to sum/length; the documented density needs n_possible.
  w0 <- c(0.5, 0.8, 0, 0.9)
  dens <- net_aggregate_weights(w0, "density", n_possible = 4)
  mn   <- net_aggregate_weights(w0, "mean")
  expect_equal(dens, sum(w0) / 4)            # true edge density
  expect_true(dens < mn)                     # diverges when zeros present
  expect_false(isTRUE(all.equal(dens, mn)))

  # cluster_summary matrix path: density == sum(block) / (n_i * n_j),
  # and differs from mean once a zero edge is injected into the block.
  set.seed(7)
  mat <- matrix(runif(64, 0.2, 1), 8, 8)
  rownames(mat) <- colnames(mat) <- LETTERS[1:8]
  cl <- c(1, 1, 1, 1, 2, 2, 2, 2)

  mat0 <- mat
  mat0[1, 5] <- 0                            # one zero edge in block (1 -> 2)
  cd0 <- cluster_summary(mat0, cl, method = "density")
  cm0 <- cluster_summary(mat0, cl, method = "mean")
  cs0 <- cluster_summary(mat0, cl, method = "sum")

  # density block-(1,2) is exactly the documented sum / number-of-possible.
  expect_equal(cd0$macro$weights[1, 2], cs0$macro$weights[1, 2] / 16)
  # ... and now differs from mean because a zero edge is present.
  expect_false(isTRUE(all.equal(cd0$macro$weights[1, 2],
                                cm0$macro$weights[1, 2])))
})
