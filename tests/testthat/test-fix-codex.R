# Regression tests for the 3 Codex P2 findings raised against the audit-fix
# diff: frequencies.R wide-id inference, plot_state_frequencies.R mosaic
# recount, hypa.R post-Filter order metadata.

test_that("Codex P2 #1: id_col=NULL never drops a start-only state", {
  # "START" appears only in the first column. The old value-overlap heuristic
  # treated column 1 as an id (its values were disjoint from later columns)
  # and silently dropped every first event. No id must be inferred from data.
  wide <- data.frame(V1 = c("START", "START"), V2 = c("A", "B"),
                     V3 = c("B", "A"), stringsAsFactors = FALSE)
  el <- convert_sequence_format(wide, format = "edgelist")
  expect_equal(nrow(el), 4L)                       # START->A,A->B,START->B,B->A
  expect_true("START" %in% c(el$from, el$to))
  fr <- convert_sequence_format(wide, format = "frequency")
  expect_true("START" %in% names(fr))
})

test_that("Codex P2 #2: mosaic uses the estimator's stored frequency_matrix", {
  set.seed(1)
  d <- data.frame(id = rep(1:12, each = 6), t = rep(1:6, 12),
                  s = sample(c("A", "B", "C"), 72, TRUE),
                  stringsAsFactors = FALSE)
  rn <- build_network(d, method = "relative", id = "id", time = "t")
  mc <- .mosaic_count_or_stop(rn)
  # The stored count construction is preferred over a naive raw recount, so
  # the mosaic agrees with the network instead of disagreeing / erroring.
  expect_identical(mc, rn$frequency_matrix)
  expect_true(all(mc == round(mc)))
  expect_no_error(mosaic_plot(rn))
})

test_that("Codex P2 #3: build_hypa $order/$k report only built layers", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "C"), c("C", "B", "D"),
                c("A", "B", "D"), c("C", "B", "A"), c("A", "B", "C"))
  h <- build_hypa(trajs, order = c(2L, 10L))   # order 10 builds no layer
  expect_identical(as.integer(h$order), as.integer(names(h$by_order)))
  expect_identical(as.integer(h$k), as.integer(names(h$by_order)))
  expect_false(10L %in% h$order)
})
