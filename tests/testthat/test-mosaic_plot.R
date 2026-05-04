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

test_that("mosaic_plot.netobject_group returns one plot per group", {
  net <- .tiny_freq_net()
  grp <- list(A = net, B = net)
  class(grp) <- "netobject_group"
  out <- mosaic_plot(grp)
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    expect_s3_class(out, "gtable")
  } else {
    expect_type(out, "list")
    expect_length(out, 2L)
  }
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
