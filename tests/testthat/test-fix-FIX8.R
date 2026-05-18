# Regression tests for FIX8 (audit findings A12-F01..F06, A05-F01, A08-F01).
#
# All assertions are structural (layer_data / syntax / object diffs) because
# plot images cannot be inspected -- we assert the argument now changes the
# output in the documented direction, or errors cleanly, or returns the
# documented shape.

# ---- A12-F01: plot_state_frequencies(node_groups=) collapses fill ----------

test_that("node_groups collapses the per-state fill onto semantic groups", {
  skip_if_not_installed("ggplot2")
  data(group_regulation_long, package = "Nestimate")
  tab <- aggregate(Course ~ Actor, group_regulation_long, function(x) x[1])
  sub <- group_regulation_long[group_regulation_long$Actor %in%
           tab$Actor[tab$Course == "A"][1:25], ]
  nw <- build_network(sub, method = "frequency", format = "long",
                       actor = "Actor", action = "Action", order = "Time")
  st <- unique(as.character(nw$nodes$label))
  ng <- setNames(rep(c("Cog", "Soc", "Meta"), length.out = length(st)), st)

  n_fill <- function(p) length(unique(ggplot2::layer_data(p$plot)$fill))
  leg <- function(p) ggplot2::ggplot_build(p$plot)$plot$labels$fill

  for (sty in c("bars", "marimekko")) {
    b0 <- plot_state_frequencies(nw, style = sty, node_groups = NULL)
    b1 <- plot_state_frequencies(nw, style = sty, node_groups = ng)
    # Output must differ NULL vs set (was byte-identical = the bug).
    expect_false(
      isTRUE(all.equal(ggplot2::layer_data(b0$plot),
                       ggplot2::layer_data(b1$plot))),
      info = sty
    )
    # 9 states collapse onto 3 groups.
    expect_equal(n_fill(b0), 9L, info = sty)
    expect_equal(n_fill(b1), 3L, info = sty)
    # Legend relabelled from "State" to "Group".
    expect_equal(leg(b0), "State", info = sty)
    expect_equal(leg(b1), "Group", info = sty)
  }

  # States absent from the mapping keep their own name (none dropped).
  ng_partial <- ng[1:3]
  b2 <- plot_state_frequencies(nw, style = "bars", node_groups = ng_partial)
  expect_equal(n_fill(b2), 9L)  # 3 mapped groups + 6 unmapped states

  # Invalid node_groups (unnamed) errors cleanly.
  expect_error(
    plot_state_frequencies(nw, style = "bars",
                           node_groups = c("Cog", "Soc")),
    "names"
  )
})

# ---- A12-F02: build_gimme(exogenous=) enters the model -----------------------

.fix8_gimme_panel <- function(seed = 11, n_sub = 4, n = 60) {
  set.seed(seed)
  do.call(rbind, lapply(seq_len(n_sub), function(i) {
    v1 <- as.numeric(arima.sim(list(ar = 0.5), n))
    v2 <- 0.6 * c(0, v1[-n]) + rnorm(n, 0, 0.5)
    v3 <- 0.4 * v2 + rnorm(n, 0, 0.6)
    v4 <- rnorm(n)
    data.frame(V1 = v1, V2 = v2, V3 = v3, V4 = v4,
               id = i, t = seq_len(n))
  }))
}

test_that("build_gimme(exogenous=) changes the fitted model", {
  skip_if_not_installed("lavaan")
  p <- .fix8_gimme_panel()
  g0 <- build_gimme(p, vars = paste0("V", 1:4), id = "id", time = "t",
                    exogenous = NULL, seed = 1)
  g1 <- build_gimme(p, vars = paste0("V", 1:4), id = "id", time = "t",
                    exogenous = "V4", seed = 1)

  # coefs / syntax must differ (was byte-identical = the bug).
  expect_false(isTRUE(all.equal(g0$coefs, g1$coefs)))
  expect_false(identical(g0$syntax, g1$syntax))

  bs <- g1$syntax[[1]]
  # V4 is exogenous: no AR self-path, never a regression outcome,
  # row of coefs all zero.
  expect_false(any(bs == "V4~V4lag"))
  expect_true(all(g1$coefs[[1]]["V4", ] == 0))
  # Endogenous->exogenous fixed to zero (nonsense paths present).
  expect_true(any(grepl("^V4~0\\*", bs)))

  # exogenous = NULL must be byte-identical to omitting the arg (no
  # behavioral regression for the default path).
  ga <- build_gimme(p, vars = paste0("V", 1:4), id = "id", time = "t",
                    seed = 1)
  expect_equal(g0$coefs, ga$coefs)
  expect_identical(g0$syntax, ga$syntax)
})

test_that("build_gimme validates exogenous names", {
  skip_if_not_installed("lavaan")
  p <- .fix8_gimme_panel()
  expect_error(
    build_gimme(p, vars = paste0("V", 1:4), id = "id", time = "t",
                exogenous = "V9"),
    "must be among"
  )
  expect_error(
    build_gimme(p, vars = paste0("V", 1:4), id = "id", time = "t",
                exogenous = paste0("V", 1:4)),
    "every variable"
  )
})

# ---- A12-F03: mosaic_plot.netobject recounts from $data ----------------------

test_that("mosaic_plot.netobject recounts non-integer nets from $data", {
  skip_if_not_installed("ggplot2")
  data(group_regulation_long, package = "Nestimate")
  tab <- aggregate(Course ~ Actor, group_regulation_long, function(x) x[1])
  sub <- group_regulation_long[group_regulation_long$Actor %in%
           tab$Actor[tab$Course == "A"][1:25], ]
  nw_rel <- build_network(sub, method = "relative", format = "long",
                          actor = "Actor", action = "Action", order = "Time")
  nw_freq <- build_network(sub, method = "frequency", format = "long",
                           actor = "Actor", action = "Action", order = "Time")
  expect_false(is.null(nw_rel$data))

  p_rel <- mosaic_plot(nw_rel, residuals = "asymptotic")
  expect_s3_class(p_rel, "ggplot")
  # Recount must equal the frequency mosaic on identical data.
  p_freq <- mosaic_plot(nw_freq, residuals = "asymptotic")
  expect_equal(ggplot2::layer_data(p_rel), ggplot2::layer_data(p_freq))

  # Strict integer guard fires only when NO count source is available:
  # non-integer weights, no $data, AND no stored $frequency_matrix. (Codex
  # P2: the estimator's stored count matrix is a valid source on its own,
  # so removing only $data no longer forces an error.)
  no_counts <- nw_rel
  no_counts$data <- NULL
  no_counts$frequency_matrix <- NULL
  no_counts$method <- "glasso"
  expect_error(mosaic_plot(no_counts), "integer-valued")
})

# ---- A12-F04: mosaic_plot data-bearing NULL angle auto-rule reachable -------

test_that("mosaic_plot data-bearing top/left_angle NULL auto-rule reachable", {
  skip_if_not_installed("ggplot2")
  data(group_regulation_long, package = "Nestimate")
  tab <- aggregate(Course ~ Actor, group_regulation_long, function(x) x[1])
  sub <- group_regulation_long[group_regulation_long$Actor %in%
           tab$Actor[tab$Course == "A"][1:25], ]
  nw <- build_network(sub, method = "frequency", format = "long",
                       actor = "Actor", action = "Action", order = "Time")

  # 9 states > 3 -> documented auto rule gives angle 90 (was hard 0).
  p <- mosaic_plot(nw, residuals = "asymptotic")
  expect_equal(ggplot2::ggplot_build(p)$plot$theme$axis.text.y$angle, 90)

  # Explicit override still respected.
  p2 <- mosaic_plot(nw, residuals = "asymptotic", left_angle = 45)
  expect_equal(ggplot2::ggplot_build(p2)$plot$theme$axis.text.y$angle, 45)

  # <= 3 levels -> auto rule gives 0.
  m <- matrix(c(5, 2, 3, 4, 6, 1, 2, 3, 7), 3, 3)
  dimnames(m) <- list(letters[1:3], letters[1:3])
  small <- nw
  small$weights <- m
  small$data <- NULL
  p3 <- mosaic_plot(small, residuals = "asymptotic")
  expect_equal(ggplot2::ggplot_build(p3)$plot$theme$axis.text.y$angle, 0)
})

# ---- A12-F06: compare_model(measures=) warns on invalid names ---------------

test_that("compare_model warns on unknown centrality measures", {
  set.seed(3)
  m1 <- matrix(runif(36, 0, 5), 6, 6)
  m2 <- matrix(runif(36, 0, 5), 6, 6)
  dimnames(m1) <- dimnames(m2) <- list(letters[1:6], letters[1:6])

  expect_warning(
    r <- compare_model(m1, m2, measures = c("Bogus1", "Bogus2"),
                        network = FALSE),
    "unknown centrality"
  )
  expect_equal(nrow(r$centrality_differences), 0L)

  # Partially invalid: warns, keeps the valid one.
  expect_warning(
    r2 <- compare_model(m1, m2, measures = c("Strength", "OutStrength"),
                         network = FALSE),
    "Strength"
  )
  expect_gt(nrow(r2$centrality_differences), 0L)

  # All valid: no warning.
  expect_no_warning(
    r3 <- compare_model(m1, m2,
                        measures = c("OutStrength", "InStrength",
                                     "Betweenness"),
                        network = FALSE)
  )
  expect_gt(nrow(r3$centrality_differences), 0L)
})

# ---- A05-F01: print.net_cluster_diagnostics(digits=) validated -------------

test_that("print.net_cluster_diagnostics validates digits", {
  seqs <- simulate_sequences(n_actors = 30, n_states = 3, seq_length = 6,
                             seed = 1)
  cl <- build_clusters(seqs, k = 2, method = "ward.D2", seed = 1)
  d <- cluster_diagnostics(cl)

  expect_error(print(d, digits = -2),
               "single non-negative whole finite number")
  expect_error(print(d, digits = 2.7),
               "single non-negative whole finite number")
  expect_error(print(d, digits = c(2, 3)),
               "single non-negative whole finite number")
  expect_error(print(d, digits = "x"),
               "single non-negative whole finite number")
  # Valid digits still work.
  expect_output(print(d, digits = 4))
})

# ---- A08-F01: plot.net_sequence_comparison(alpha=) controls display --------

test_that("plot.net_sequence_comparison alpha gates p-value display", {
  skip_if_not_installed("ggplot2")
  set.seed(101)
  g1 <- as.data.frame(matrix(sample(c("A", "B", "C"), 80 * 5, TRUE,
                                     prob = c(.55, .3, .15)), nrow = 80),
                       stringsAsFactors = FALSE)
  g2 <- as.data.frame(matrix(sample(c("A", "B", "C"), 80 * 5, TRUE,
                                     prob = c(.2, .3, .5)), nrow = 80),
                       stringsAsFactors = FALSE)
  names(g1) <- names(g2) <- paste0("T", 1:5)
  seqs <- rbind(g1, g2)
  grp <- rep(c("G1", "G2"), each = 80)
  res <- sequence_compare(seqs, group = grp, sub = 2, min_freq = 4L,
                          test = "chisq", adjust = "none")

  p_strict <- plot(res, style = "pyramid", alpha = 0.001, top_n = 8)
  p_loose <- plot(res, style = "pyramid", alpha = 0.999, top_n = 8)
  b1 <- ggplot2::ggplot_build(p_strict)
  b2 <- ggplot2::ggplot_build(p_loose)

  # Layers must differ between alpha values (was byte-identical = the bug).
  ident <- all(mapply(
    function(d1, d2) isTRUE(all.equal(d1, d2, check.attributes = FALSE)),
    b1$data, b2$data
  ))
  expect_false(ident)

  # The p-label text layer: a non-significant pattern (p ~ 0.34) is starred
  # under alpha = 0.999 but not under alpha = 0.001.
  plabels <- function(b) {
    txt <- which(vapply(b$plot$layers,
                        function(L) inherits(L$geom, "GeomText"),
                        logical(1)))
    b$data[[utils::tail(txt, 1)]]$label
  }
  l_strict <- plabels(b1)
  l_loose <- plabels(b2)
  expect_true(any(grepl("\\*", l_strict)))           # significant ones starred
  expect_gt(sum(grepl("\\*", l_loose)),
            sum(grepl("\\*", l_strict)))              # looser alpha stars more
})
