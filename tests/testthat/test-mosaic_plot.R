# Tests for mosaic_plot() — the tna-equivalent chi-square mosaic for netobjects.

skip_if_not_installed("ggplot2")

# Build a small frequency netobject the test can lean on.
.tiny_freq_net <- function() {
  set.seed(13)
  seqs <- replicate(40, sample(c("A", "B", "C"), size = 6, replace = TRUE),
                    simplify = FALSE)
  df <- do.call(rbind, lapply(seq_along(seqs), function(i) {
    data.frame(id = i, time = seq_along(seqs[[i]]), state = seqs[[i]],
               stringsAsFactors = FALSE)
  }))
  build_network(df, method = "frequency",
                id_col = "id", time_col = "time", action = "state")
}

test_that("mosaic_plot.netobject returns a ggplot for a frequency network", {
  net <- .tiny_freq_net()
  p <- mosaic_plot(net)
  expect_s3_class(p, "ggplot")
})

test_that("mosaic_plot rejects non-integer-weighted networks", {
  set.seed(13)
  seqs <- replicate(20, sample(c("A", "B", "C"), size = 6, replace = TRUE),
                    simplify = FALSE)
  df <- do.call(rbind, lapply(seq_along(seqs), function(i) {
    data.frame(id = i, time = seq_along(seqs[[i]]), state = seqs[[i]],
               stringsAsFactors = FALSE)
  }))
  net_rel <- build_network(df, method = "relative",
                           id_col = "id", time_col = "time", action = "state")
  expect_error(mosaic_plot(net_rel), "integer-valued")
})

test_that("mosaic_plot.netobject_group returns one (group x state) mosaic", {
  net <- .tiny_freq_net()
  grp <- list(A = net, B = net)
  class(grp) <- "netobject_group"
  out <- mosaic_plot(grp, seed = 1L)
  expect_s3_class(out, "ggplot")
  # Single panel containing both groups -- not faceted.
  ld <- ggplot2::layer_data(out, 1L)
  expect_setequal(as.character(unique(ld$PANEL)), "1")
  # 2 groups x 3 states = 6 rectangles in the (group x state) contingency.
  expect_equal(nrow(ld), 6L)
  # Two distinct x-extents (one column per group, group totals identical).
  x_widths <- unique(round(ld$xmax - ld$xmin, 4))
  expect_length(x_widths, 1L)
})

test_that("mosaic_plot.htna defers to the netobject path", {
  # htna inherits from netobject; mosaic of $weights should match exactly
  # what mosaic_plot.netobject would draw on the same matrix.
  net <- .tiny_freq_net()
  ht  <- net
  class(ht) <- c("htna", class(net))
  # Same seed so the permutation residuals match cell-for-cell.
  p_ht  <- mosaic_plot(ht,  seed = 1L)
  p_net <- mosaic_plot(net, seed = 1L)
  expect_s3_class(p_ht, "ggplot")
  expect_equal(ggplot2::layer_data(p_ht), ggplot2::layer_data(p_net))
})

# Fixture: a real mcml from group_regulation_long with a theory-driven
# 3-way clustering of the SRL action codes. Used by the next three tests.
.tiny_mcml <- function() {
  build_mcml(
    group_regulation_long,
    clusters = list(
      Cognition     = c("discuss", "synthesis", "consensus"),
      Metacognition = c("plan", "monitor", "adapt"),
      Social        = c("coregulate", "cohesion", "emotion")
    ),
    method = "sum", type = "frequency",
    actor = "Actor", time = "Time", action = "Action"
  )
}

test_that("mosaic_plot.mcml level='macro' draws one cluster x cluster panel", {
  fit <- .tiny_mcml()
  p <- mosaic_plot(fit, seed = 1L)
  expect_s3_class(p, "ggplot")
})

test_that("mosaic_plot.mcml level='clusters' returns a faceted ggplot", {
  fit <- .tiny_mcml()
  out <- mosaic_plot(fit, level = "clusters", seed = 1L)
  expect_s3_class(out, "ggplot")
  # One facet per cluster.
  ld <- ggplot2::layer_data(out, 1L)
  expect_setequal(as.character(unique(ld$PANEL)),
                  as.character(seq_along(fit$clusters)))
})

test_that("mosaic_plot.mcml level='clusters' errors when $clusters is empty", {
  fit <- structure(list(macro = list(weights = diag(3)), clusters = NULL),
                   class = "mcml")
  expect_error(mosaic_plot(fit, level = "clusters"),
               "x\\$clusters is empty")
})

test_that("mosaic_plot default method names the four data-bearing classes", {
  expect_error(mosaic_plot(1:5),
               "netobject, netobject_group, mcml, htna")
})

test_that("permutation residuals converge to asymptotic stdres on large N", {
  set.seed(101)
  rs <- c(800, 200, 400)
  cs <- c(300, 500, 600)
  probs <- outer(rs / sum(rs), cs / sum(cs))
  counts <- stats::rmultinom(1L, 5000L, as.vector(probs))
  tab <- as.table(matrix(counts, 3, 3,
                         dimnames = list(LETTERS[1:3], letters[1:3])))
  z_perm <- Nestimate:::.mosaic_perm_stdres(tab, n_perm = 2000L, seed = 1L)
  z_asy  <- suppressWarnings(stats::chisq.test(tab))$stdres
  expect_equal(as.numeric(z_perm), as.numeric(z_asy), tolerance = 0.2)
})

test_that("mosaic_plot geometry matches tna::plot_mosaic", {
  skip_if_not_installed("tna")
  set.seed(7)
  seqs <- replicate(80, sample(c("A", "B", "C"), size = 8, replace = TRUE),
                    simplify = FALSE)
  wide <- do.call(rbind, lapply(seqs, function(s) data.frame(t(s))))
  names(wide) <- paste0("T", seq_len(ncol(wide)))

  tna_model <- tna::ftna(wide)
  nest_net <- build_network(wide, method = "frequency", format = "wide")

  p_tna  <- tna::plot_mosaic(tna_model)
  p_nest <- mosaic_plot(nest_net, range = c(-4, 4),
                        residuals = "asymptotic")

  d_tna  <- ggplot2::ggplot_build(p_tna)$data[[1]]
  d_nest <- ggplot2::ggplot_build(p_nest)$data[[1]]

  o <- order(d_tna$xmin, d_tna$ymin)
  d_tna  <- d_tna[o, c("xmin", "xmax", "ymin", "ymax", "fill")]
  o <- order(d_nest$xmin, d_nest$ymin)
  d_nest <- d_nest[o, c("xmin", "xmax", "ymin", "ymax", "fill")]

  expect_equal(nrow(d_tna), nrow(d_nest))
  expect_equal(d_tna$xmin, d_nest$xmin, tolerance = 1e-10)
  expect_equal(d_tna$xmax, d_nest$xmax, tolerance = 1e-10)
  expect_equal(d_tna$ymin, d_nest$ymin, tolerance = 1e-10)
  expect_equal(d_tna$ymax, d_nest$ymax, tolerance = 1e-10)
  expect_equal(d_tna$fill, d_nest$fill)
})
