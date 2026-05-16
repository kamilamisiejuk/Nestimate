testthat::skip_on_cran()

# ==============================================================================
# Tests for the unified clustering print methods:
#   - print.net_clustering           (R/cluster_data.R)
#   - print.net_mmm                  (R/mmm.R)
#   - print.net_mmm_clustering       (R/mmm.R)
#   - print.netobject_group          (R/build_network.R)
#
# Goals:
#   1. Headers / column shapes match the agreed layout.
#   2. New `digits` arg is honoured but optional (non-breaking).
#   3. Old behaviours preserved (returns invisibly, no error on minimal input).
#   4. Regression: cluster_network/cluster_mmm output flows through
#      sequence_plot()/distribution_plot() without length mismatch (the
#      .extract_seqplot_input bug fixed in this change).
# ==============================================================================

# ---- Synthetic data ----------------------------------------------------------

.make_clust_data <- function(n_per = 20, n_cols = 8, seed = 7) {
  set.seed(seed)
  g1 <- matrix(sample(c("A", "B", "C"), n_per * n_cols, replace = TRUE,
                      prob = c(0.7, 0.2, 0.1)), ncol = n_cols)
  g2 <- matrix(sample(c("A", "B", "C"), n_per * n_cols, replace = TRUE,
                      prob = c(0.1, 0.2, 0.7)), ncol = n_cols)
  out <- as.data.frame(rbind(g1, g2), stringsAsFactors = FALSE)
  colnames(out) <- paste0("T", seq_len(n_cols))
  out
}

# ==============================================================================
# print.net_clustering
# ==============================================================================

test_that("print.net_clustering: header shape and dimension line", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2", dissimilarity = "hamming")
  out <- capture.output(print(cl))

  expect_match(out[1L], "^Sequence Clustering \\[ward\\.D2\\]")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("Clusters: 2", out)))
  expect_true(any(grepl("Dissimilarity: hamming", out)))
  expect_true(any(grepl("silhouette = ", out)))
})

test_that("print.net_clustering: cluster table has Cluster, N, sizes sum to total", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  out <- capture.output(print(cl))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b", out[hdr]))

  # extract numeric sizes from the "N (pct%)" column on the next 2 rows
  body <- out[seq.int(hdr + 1L, hdr + 2L)]
  sizes <- as.integer(sub("^\\s*\\d+\\s+(\\d+)\\s+.*", "\\1", body))
  expect_length(sizes, 2L)
  expect_equal(sum(sizes), 40L)
})

test_that("print.net_clustering: PAM emits Medoid column", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "pam")
  out <- capture.output(print(cl))
  hdr <- out[grep("^\\s*Cluster\\b", out)]
  expect_true(grepl("Medoid", hdr))
})

test_that("print.net_clustering: digits arg controls float decimals", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")

  o3 <- capture.output(print(cl, digits = 3))
  o6 <- capture.output(print(cl, digits = 6))

  # The silhouette line should differ between digits = 3 and digits = 6
  s3 <- o3[grep("silhouette", o3)]
  s6 <- o6[grep("silhouette", o6)]
  expect_false(identical(s3, s6))
  expect_match(s3, "silhouette = -?\\d+\\.\\d{3}")
  expect_match(s6, "silhouette = -?\\d+\\.\\d{6}")
})

test_that("print.net_clustering: returns invisibly with no extra args (non-breaking)", {
  d <- .make_clust_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  res <- withVisible(print(cl))
  expect_identical(res$value, cl)
  expect_false(res$visible)
})

# ==============================================================================
# print.net_mmm
# ==============================================================================

test_that("print.net_mmm: header / IC line / no leftover '%%'", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  out <- capture.output(print(mmm))

  expect_match(out[1L], "^Mixed Markov Model$")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("Clusters: 2", out)))
  expect_true(any(grepl("States: 3", out)))
  expect_true(any(grepl("BIC = ", out)))
  expect_true(any(grepl("AIC = ", out)))
  expect_true(any(grepl("ICL = ", out)))

  # Regression for the literal "%%" bug in the previous print method.
  expect_false(any(grepl("%%", out)))
})

test_that("print.net_mmm: cluster table shape (Cluster, N, Mix%, AvePP)", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  out <- capture.output(print(mmm))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b",     out[hdr]))
  expect_true(grepl("Mix%",        out[hdr]))
  expect_true(grepl("AvePP",       out[hdr]))

  body <- out[seq.int(hdr + 1L, hdr + 2L)]
  sizes <- as.integer(sub("^\\s*\\d+\\s+(\\d+)\\s+.*", "\\1", body))
  expect_equal(sum(sizes), 40L)
})

test_that("print.net_mmm: digits arg honoured + invisible return", {
  d <- .make_clust_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)

  out6 <- capture.output(print(mmm, digits = 6))
  ll_line <- out6[grep("BIC =", out6)]
  expect_match(ll_line, "BIC = -?\\d+\\.\\d{6}")

  res <- withVisible(print(mmm))
  expect_identical(res$value, mmm)
  expect_false(res$visible)
})

# ==============================================================================
# print.net_mmm_clustering
# ==============================================================================

test_that("print.net_mmm_clustering: header + cluster table", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)

  expect_s3_class(grp, "netobject_group")
  cl <- attr(grp, "clustering")
  expect_s3_class(cl, "net_mmm_clustering")

  out <- capture.output(print(cl))
  expect_match(out[1L], "^MMM Clustering \\[k = 2\\]")
  expect_true(any(grepl("Sequences: 40", out)))
  expect_true(any(grepl("AvePP =", out)))
  expect_true(any(grepl("BIC =", out)))

  hdr <- grep("^\\s*Cluster\\b", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Mix%",  out[hdr]))
  expect_true(grepl("AvePP", out[hdr]))
})

test_that("net_mmm_clustering methods reject unsupported dots and bad digits", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  expect_error(print(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(plot(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(print(cl, digits = NA_real_),
               "'digits' must be a single non-negative whole")
  expect_error(print(cl, digits = 2.5),
               "'digits' must be a single non-negative whole")
})

test_that("net_mmm_clustering print sizes match component member data", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  assignment_sizes <- as.integer(tabulate(cl$assignments, nbins = cl$k))
  component_rows <- as.integer(vapply(grp, function(x) nrow(x$data), integer(1L)))
  expect_equal(assignment_sizes, component_rows)
})

test_that("cluster_mmm stashes full sequence data on the clustering attribute", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  # The invariant set by cluster_mmm + .extract_seqplot_input depends on
  # this carrying an N-row data frame, where N matches assignments.
  expect_false(is.null(cl$data))
  expect_equal(NROW(cl$data), length(cl$assignments))
  expect_equal(NROW(cl$data), nrow(d))
})

# ==============================================================================
# print.netobject_group
# ==============================================================================

test_that("print.netobject_group (group_col): legacy header preserved", {
  d <- .make_clust_data()
  d$grp <- rep(c("X", "Y"), each = 20)
  nets <- build_network(d, method = "relative", group = "grp")

  out <- capture.output(print(nets))
  expect_match(out[1L], "^Group Networks \\(2 groups, group_col: grp\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Nodes",   out[hdr]))
  expect_true(grepl("Edges",   out[hdr]))
  expect_true(grepl("Weights", out[hdr]))
  # No N column without a clustering attribute.
  expect_false(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group (cluster_network): clustering source surfaced", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2",
                         dissimilarity = "hamming")

  out <- capture.output(print(grp))
  expect_match(out[1L],
               "^Group Networks \\(2 clusters via ward\\.D2 / hamming\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_length(hdr, 1L)
  expect_true(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group (cluster_mmm): MMM source surfaced", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)

  out <- capture.output(print(grp))
  expect_match(out[1L], "^Group Networks \\(2 clusters from MMM\\)")

  hdr <- grep("^\\s+Group\\s+Nodes", out)
  expect_true(grepl("\\bN\\b", out[hdr]))
})

test_that("print.netobject_group: digits arg + invisible return", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2")

  out6 <- capture.output(print(grp, digits = 6))
  weights_lines <- out6[grep("\\[\\d", out6)]
  # weight ranges shown as [a, b] with `digits` decimals
  expect_true(any(grepl("\\[\\d+\\.\\d{6}, \\d+\\.\\d{6}\\]", weights_lines)))

  res <- withVisible(print(grp))
  expect_identical(res$value, grp)
  expect_false(res$visible)
})

# ==============================================================================
# Regression: .extract_seqplot_input no longer mismatches data/group length
# ==============================================================================

test_that("sequence_plot accepts cluster_network output without length mismatch", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "ward.D2")

  # Previous bug: x[[1]]$data was n1 rows but assignments was 40 ->
  # stopifnot(length(group) == nrow(x)) inside distribution_plot would fire,
  # and the index path would silently take only cluster 1's data.
  expect_silent(p <- sequence_plot(grp, type = "index"))
  expect_silent(p <- distribution_plot(grp))
})

test_that("sequence_plot accepts cluster_mmm output without length mismatch", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)

  expect_silent(p <- sequence_plot(grp, type = "index"))
  expect_silent(p <- distribution_plot(grp))
})

test_that("sequence_plot accepts plain netobject_group without clustering attr", {
  d <- .make_clust_data()
  d$grp <- rep(c("X", "Y"), each = 20)
  nets <- build_network(d, method = "relative", group = "grp")

  # Plain group_col split: helper rebuilds full data via rbind() and labels
  # rows by their group. Should not error.
  expect_silent(p <- sequence_plot(nets, type = "index"))
})

# ==============================================================================
# cluster_network(cluster_by = "mmm") + build_network(net_mmm) interplay
# ==============================================================================

test_that("cluster_network(cluster_by='mmm') forwards MMM args to build_mmm", {
  d <- .make_clust_data()

  # If n_starts is forwarded, two seeds give reproducible-but-different
  # numerical state. If it's silently dropped (the bug), build_mmm runs with
  # the default n_starts and seed; both calls would still be free to differ
  # by chance, but the converged log-likelihood with n_starts = 50 is
  # rarely worse than with n_starts = 1, and we can verify by comparing to
  # the explicit build_mmm call.
  grp_a <- cluster_network(d, k = 2, cluster_by = "mmm",
                           n_starts = 1, max_iter = 5, seed = 1)
  expect_s3_class(grp_a, "netobject_group")

  # If n_starts/max_iter were silently dropped, the LL would match the
  # default n_starts = 50 fit, not this 1-start, 5-iter fit. Compare via
  # the clustering attribute attached to the netobject_group.
  cl_a  <- attr(grp_a, "clustering")
  ref   <- build_mmm(d, k = 2, n_starts = 1, max_iter = 5, seed = 1)
  expect_equal(cl_a$log_likelihood, ref$log_likelihood, tolerance = 1e-6)
})

test_that("cluster_network(cluster_by='mmm') sets attr(, 'clustering') so seqplot panels differ", {
  d <- .make_clust_data()
  grp <- cluster_network(d, k = 2, cluster_by = "mmm",
                         n_starts = 1, max_iter = 20, seed = 1)
  cl  <- attr(grp, "clustering")

  # Without this attribute, .extract_seqplot_input falls back to rbind-ing
  # each member's $data -- which for net_mmm -> netobject_group are 2
  # copies of the full dataset, producing identical panels.
  expect_s3_class(cl, "net_mmm_clustering")
  expect_false(is.null(cl$data))
  expect_equal(NROW(cl$data), length(cl$assignments))
  expect_equal(NROW(cl$data), nrow(d))
  expect_silent(p <- sequence_plot(grp, type = "index"))
})

test_that("cluster_mmm accepts cluster_by + extra ... without erroring", {
  d <- .make_clust_data()

  # API parity with cluster_network(): the unified surface
  # cluster_*(data, k, cluster_by = ..., ...) shouldn't blow up because
  # cluster_mmm() got a name it doesn't use internally.
  expect_silent(grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20,
                                    seed = 1, cluster_by = "mmm"))
  expect_s3_class(grp, "netobject_group")

  # But mistyping cluster_by should produce a clear error, not run
  # something else.
  expect_error(
    cluster_mmm(d, k = 2, n_starts = 1, seed = 1, cluster_by = "pam"),
    "only supports cluster_by"
  )
})

# ==============================================================================
# cluster_data() deprecation alias
# ==============================================================================

# ==============================================================================
# plot.net_mmm_clustering (the attribute object plot method)
# ==============================================================================

test_that("plot.net_mmm_clustering: posterior produces a ggplot", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  p <- plot(cl, type = "posterior")
  expect_s3_class(p, "ggplot")
})

test_that("plot.net_mmm_clustering: predictors is alias for covariates, errors w/o covariates", {
  d <- .make_clust_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  cl <- attr(grp, "clustering")

  # No covariates set on this run -> both predictors and covariates error
  # with the same message.
  expect_error(plot(cl, type = "covariates"), "No covariate analysis")
  expect_error(plot(cl, type = "predictors"), "No covariate analysis")
})

# ==============================================================================
# cluster_data() deprecation alias
# ==============================================================================

test_that("cluster_data() forwards to build_clusters() with a deprecation warning", {
  d <- .make_clust_data()

  # Deprecation warning is raised once per session, so use
  # withCallingHandlers + .Deprecated's "deprecatedWarning" condition class
  # rather than expect_warning() (which can swallow only the first call).
  saw <- FALSE
  withCallingHandlers(
    cl <- cluster_data(d, k = 2, method = "ward.D2"),
    deprecatedWarning = function(w) {
      saw <<- TRUE
      invokeRestart("muffleWarning")
    }
  )
  expect_true(saw)
  expect_s3_class(cl, "net_clustering")
  expect_equal(cl$k, 2L)
})

test_that("cluster_data() is equivalent to direct build_clusters()", {
  d <- .make_clust_data()

  muffled_cluster_data <- function(...) {
    withCallingHandlers(
      cluster_data(...),
      deprecatedWarning = function(w) invokeRestart("muffleWarning")
    )
  }

  alias <- muffled_cluster_data(d, k = 2, method = "ward.D2",
                                dissimilarity = "hamming", seed = 7)
  ref <- build_clusters(d, k = 2, method = "ward.D2",
                        dissimilarity = "hamming", seed = 7)

  expect_identical(class(alias), class(ref))
  expect_identical(unname(alias$assignments), unname(ref$assignments))
  expect_identical(alias$sizes, ref$sizes)
  expect_equal(alias$silhouette, ref$silhouette, tolerance = 1e-12)
  expect_equal(as.matrix(alias$distance), as.matrix(ref$distance),
               tolerance = 1e-12)
})

test_that("cluster_data() inherits build_clusters strict argument rejection", {
  d <- .make_clust_data()

  expect_error(
    withCallingHandlers(
      cluster_data(d, k = 2, typo_arg = TRUE),
      deprecatedWarning = function(w) invokeRestart("muffleWarning")
    ),
    "build_clusters\\(\\) got unsupported argument: typo_arg"
  )
})
