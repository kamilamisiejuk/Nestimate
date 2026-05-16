testthat::skip_on_cran()

# ==============================================================================
# Tests for cluster_diagnostics() and the net_cluster_diagnostics class.
#
# Coverage:
#   1. Each dispatch path (net_clustering, net_mmm, netobject_group with
#      either clustering attribute, net_mmm_clustering directly).
#   2. Per-cluster math: sizes sum to N; percentages sum to 100; sil_mean
#      averaged across observations equals cl$silhouette; class_err_pct
#      averaged (size-weighted) equals quality$classification_error.
#   3. The net_mmm and net_mmm_clustering paths produce identical
#      diagnostics (the second method is a thin shim over the first).
#   4. Print output: header, per-family lines, table rows.
#   5. Plot delegation: the resulting ggplot has the title that
#      plot.net_clustering / plot.net_mmm_clustering would produce on the
#      source object directly.
#   6. as.data.frame() returns the per_cluster table.
#   7. Default method errors with a useful message.
# ==============================================================================

# ---- Synthetic data ----------------------------------------------------------

.make_diag_data <- function(n_per = 25, n_cols = 8, seed = 13) {
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
# Dispatch + structure
# ==============================================================================

test_that("cluster_diagnostics(net_clustering): structure + per_cluster math", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)

  expect_s3_class(diag, "net_cluster_diagnostics")
  expect_identical(diag$family, "distance")
  expect_equal(diag$k, 2L)
  expect_equal(diag$n, 50L)
  expect_equal(sum(diag$sizes), 50L)
  expect_identical(diag$ics, NULL)

  pc <- diag$per_cluster
  expect_equal(nrow(pc), 2L)
  expect_setequal(names(pc),
    c("cluster", "size", "pct", "mean_within_dist", "sil_mean"))
  expect_equal(sum(pc$size), 50L)
  expect_equal(sum(pc$pct),  100, tolerance = 1e-9)
  # Size-weighted sil_mean ought to equal the overall silhouette to
  # within float noise: that's the formula cluster::silhouette uses.
  expect_equal(sum(pc$sil_mean * pc$size) / sum(pc$size),
               diag$overall$silhouette, tolerance = 1e-9)
})

test_that("cluster_diagnostics(net_mmm): structure + class_err_pct decomposition", {
  d <- .make_diag_data()
  mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  diag <- cluster_diagnostics(mmm)

  expect_s3_class(diag, "net_cluster_diagnostics")
  expect_identical(diag$family, "mmm")
  expect_equal(diag$k, 2L)
  expect_equal(diag$n, 50L)
  expect_false(is.null(diag$ics))
  expect_setequal(names(diag$ics),
    c("BIC", "AIC", "ICL", "log_likelihood"))

  pc <- diag$per_cluster
  expect_setequal(names(pc),
    c("cluster", "size", "pct", "mix_pct", "avepp", "class_err_pct"))
  expect_equal(sum(pc$size), 50L)
  expect_equal(sum(pc$mix_pct), 100, tolerance = 1e-6)

  # Size-weighted mean of per-cluster class_err_pct must equal
  # 100 * quality$classification_error -- that's the definition.
  weighted_err <- sum(pc$class_err_pct * pc$size) / sum(pc$size)
  expect_equal(weighted_err / 100,
               diag$overall$classification_error, tolerance = 1e-9)
})

test_that("netobject_group dispatch routes to attr(, 'clustering')", {
  d <- .make_diag_data()

  # Distance-clustering group -> net_clustering attr -> distance family.
  grp_d <- cluster_network(d, k = 2, cluster_by = "ward.D2")
  d_diag <- cluster_diagnostics(grp_d)
  expect_identical(d_diag$family, "distance")
  expect_equal(d_diag$k, 2L)

  # MMM group -> net_mmm_clustering attr -> mmm family.
  grp_m <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  m_diag <- cluster_diagnostics(grp_m)
  expect_identical(m_diag$family, "mmm")
  expect_equal(m_diag$k, 2L)
})

test_that("netobject_group without clustering attribute errors helpfully", {
  d <- .make_diag_data()
  d$grp <- rep(c("X", "Y"), each = 25)
  nets <- build_network(d, method = "relative", group = "grp")

  expect_error(cluster_diagnostics(nets),
               "requires a clustering attribute")
})

test_that("net_mmm_clustering produces same diagnostics as net_mmm", {
  d <- .make_diag_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)
  cl_attr <- attr(grp, "clustering")

  # Compare to a direct build_mmm on the same call. Use the same
  # seed/n_starts -- they should produce identical EM state, so the
  # diagnostics must agree on every numeric except for the source ref.
  ref_mmm <- build_mmm(d, k = 2, n_starts = 2, max_iter = 30, seed = 1)

  d1 <- cluster_diagnostics(cl_attr)
  d2 <- cluster_diagnostics(ref_mmm)

  expect_identical(d1$family, "mmm")
  expect_equal(d1$k, d2$k)
  expect_equal(d1$n, d2$n)
  expect_equal(d1$per_cluster, d2$per_cluster)
  expect_equal(d1$overall, d2$overall)
  expect_equal(d1$ics, d2$ics)
  # source is class-specific (one is net_mmm_clustering, the other net_mmm).
  expect_s3_class(d1$source, "net_mmm_clustering")
  expect_s3_class(d2$source, "net_mmm")
})

test_that("default method errors clearly on unsupported input", {
  expect_error(cluster_diagnostics(list(a = 1)),
               "no method for class")
})

test_that("cluster_diagnostics rejects unsupported extra arguments", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  mmm <- build_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)

  expect_error(cluster_diagnostics(cl, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(cluster_diagnostics(mmm, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(cluster_diagnostics(attr(grp, "clustering"), typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(cluster_diagnostics(grp, typo_arg = TRUE),
               "unsupported argument: typo_arg")
})

# ==============================================================================
# Print
# ==============================================================================

test_that("print(): distance family header + per_cluster table", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)
  out <- capture.output(print(diag))

  expect_match(out[1L], "^Cluster Diagnostics \\(distance\\)")
  expect_true(any(grepl("Sequences: 50", out)))
  expect_true(any(grepl("Clusters: 2", out)))
  expect_true(any(grepl("silhouette = ", out)))

  hdr <- grep("^\\s+Cluster\\s+N\\s", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Mean within-dist", out[hdr]))
  expect_true(grepl("Silhouette",       out[hdr]))
})

test_that("print(): mmm family header + IC line + per_cluster table", {
  d <- .make_diag_data()
  mmm <- build_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  diag <- cluster_diagnostics(mmm)
  out <- capture.output(print(diag))

  expect_match(out[1L], "^Cluster Diagnostics \\(mmm\\)")
  expect_true(any(grepl("BIC = ", out)))
  expect_true(any(grepl("AIC = ", out)))
  expect_true(any(grepl("ICL = ", out)))

  hdr <- grep("^\\s+Cluster\\s+N\\s", out)
  expect_length(hdr, 1L)
  expect_true(grepl("Mix%",       out[hdr]))
  expect_true(grepl("AvePP",      out[hdr]))
  expect_true(grepl("Class.Err%", out[hdr]))
})

test_that("print() returns invisibly with default args (non-breaking)", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)

  res <- withVisible(print(diag))
  expect_identical(res$value, diag)
  expect_false(res$visible)
})

test_that("print and as.data.frame diagnostics reject unsupported dots", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)

  expect_error(print(diag, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(as.data.frame(diag, typo_arg = TRUE),
               "unsupported argument: typo_arg")
})

# ==============================================================================
# Plot delegation
# ==============================================================================

test_that("plot.net_cluster_diagnostics delegates to plot.net_clustering", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)

  for (t in c("silhouette", "mds", "heatmap")) {
    p <- plot(diag, type = t)
    # Pin which `type` failed by including it in the assertion label.
    testthat::expect_true(inherits(p, "ggplot"),
                          info = sprintf("type = %s", t))
  }
  # Default type for distance family is "silhouette".
  expect_s3_class(plot(diag), "ggplot")
})

test_that("plot.net_cluster_diagnostics delegates to plot.net_mmm_clustering for MMM groups", {
  d <- .make_diag_data()
  grp <- cluster_mmm(d, k = 2, n_starts = 1, max_iter = 20, seed = 1)
  diag <- cluster_diagnostics(grp)

  expect_s3_class(plot(diag, type = "posterior"), "ggplot")
  # Default type for mmm family is "posterior".
  expect_s3_class(plot(diag), "ggplot")
})

# ==============================================================================
# as.data.frame
# ==============================================================================

test_that("as.data.frame returns the per_cluster table", {
  d <- .make_diag_data()
  cl <- build_clusters(d, k = 2, method = "ward.D2")
  diag <- cluster_diagnostics(cl)

  df <- as.data.frame(diag)
  expect_s3_class(df, "data.frame")
  expect_identical(df, diag$per_cluster)
  expect_equal(nrow(df), 2L)
})
