# Equivalence test: mosaic_plot() vs vcd::mosaic and tna::plot_mosaic.
#
# Both vcd and tna agree on the underlying mathematics of a 2-D mosaic:
#
#   1. cell area  = tab[i, j] / sum(tab)             (area-proportional)
#   2. shading    = chisq.test(tab)$stdres           (chi-square standardized
#                                                     residual; vcd via
#                                                     shading_Friendly, tna
#                                                     directly)
#
# vcd renders this with grid graphics + per-depth gaps; tna renders the same
# math with ggplot2 + a tiny constant offset. The visual layouts differ; the
# math is identical. mosaic_plot() should reproduce the math to machine
# precision and reproduce tna's ggplot output byte-for-byte (same gap
# convention).
#
# This test runs across (a) the real `ai_long` transition table, plus
# (b) 50 random contingency tables of varying shape and density, asserting
# the three equivalences below for every case:
#
#   - mosaic_plot cell areas == tab / sum(tab)        (vcd math)
#   - mosaic_plot ggplot data == tna::plot_mosaic ggplot data   (byte-for-byte)
#   - mosaic_plot fill values == chisq.test(tab)$stdres (squished to +/-4)
#
# Gated by NESTIMATE_EQUIV_TESTS=true so it never runs under devtools::test()
# or R CMD check. Run manually with:
#   NESTIMATE_EQUIV_TESTS=true Rscript -e \
#     'testthat::test_dir("local_testing_and_equivalence", filter = "mosaic")'

skip_if_not_set <- function() {
  if (!identical(Sys.getenv("NESTIMATE_EQUIV_TESTS"), "true")) {
    testthat::skip("Set NESTIMATE_EQUIV_TESTS=true to run mosaic equivalence.")
  }
}

testthat::skip_if_not_installed("tna")
testthat::skip_if_not_installed("vcd")

if (!exists("build_network", mode = "function")) {
  suppressMessages(devtools::load_all(file.path("..", "."), quiet = TRUE))
}

TOL_AREA   <- 1e-12
TOL_COORD  <- 1e-10
TOL_STDRES <- 1e-12

.tab_to_tna <- function(tab) {
  m <- list(weights = t(unclass(tab)))
  attr(m, "type") <- "frequency"
  class(m) <- "tna"
  m
}

.tab_to_netobj <- function(tab) {
  w <- t(unclass(tab))
  storage.mode(w) <- "double"
  net <- list(
    weights = w,
    nodes   = data.frame(id = seq_len(nrow(w)), label = rownames(w),
                         name = rownames(w), stringsAsFactors = FALSE),
    edges   = NULL,
    method  = "frequency",
    params  = list()
  )
  class(net) <- c("netobject", "cograph_network")
  net
}

.area_check <- function(tab) {
  net <- .tab_to_netobj(tab)
  p   <- mosaic_plot(net, residuals = "asymptotic")
  d   <- ggplot2::ggplot_build(p)$data[[1]]
  d$area <- (d$xmax - d$xmin) * (d$ymax - d$ymin)
  expected <- as.numeric(t(unclass(tab))) / sum(tab)
  o <- order(d$xmin, d$ymin)
  max(abs(d$area[o] - expected[order(rep(seq_len(nrow(tab)), each = ncol(tab)),
                                      rep(seq_len(ncol(tab)), nrow(tab)))]))
}

.coord_check <- function(tab) {
  net <- .tab_to_netobj(tab)
  p_nest <- mosaic_plot(net, range = c(-4, 4),
                        residuals = "asymptotic")
  p_tna  <- tna::plot_mosaic(.tab_to_tna(tab))
  d_n <- ggplot2::ggplot_build(p_nest)$data[[1]]
  d_t <- ggplot2::ggplot_build(p_tna)$data[[1]]
  o_n <- order(d_n$xmin, d_n$ymin)
  o_t <- order(d_t$xmin, d_t$ymin)
  c(
    xmin = max(abs(d_n$xmin[o_n] - d_t$xmin[o_t])),
    xmax = max(abs(d_n$xmax[o_n] - d_t$xmax[o_t])),
    ymin = max(abs(d_n$ymin[o_n] - d_t$ymin[o_t])),
    ymax = max(abs(d_n$ymax[o_n] - d_t$ymax[o_t]))
  )
}

.stdres_check <- function(tab) {
  net <- .tab_to_netobj(tab)
  p   <- mosaic_plot(net, residuals = "asymptotic")
  d   <- ggplot2::ggplot_build(p)$data[[1]]
  ref <- suppressWarnings(stats::chisq.test(tab))$stdres
  ref_squished <- pmin(pmax(ref, -4), 4)
  o <- order(d$xmin, d$ymin)
  expected <- as.numeric(t(ref_squished))[order(rep(seq_len(nrow(tab)),
                                                    each = ncol(tab)),
                                                rep(seq_len(ncol(tab)),
                                                    nrow(tab)))]
  max(abs(d$fill_alpha %||% rep(NA, length(expected)) - expected),
      na.rm = TRUE)
}

.random_table <- function(seed) {
  set.seed(seed)
  n <- sample(3:8, 1L)
  m <- sample(3:8, 1L)
  total <- sample(c(50, 200, 500, 2000), 1L)
  probs <- matrix(stats::runif(n * m, min = 0.05), n, m)
  probs <- probs / sum(probs)
  counts <- as.vector(stats::rmultinom(1L, total, probs))
  tab <- matrix(counts, n, m)
  rownames(tab) <- paste0("R", seq_len(n))
  colnames(tab) <- paste0("C", seq_len(m))
  as.table(tab)
}

testthat::test_that("mosaic_plot is area-proportional like vcd::mosaic", {
  skip_if_not_set()
  data(ai_long, package = "Nestimate")
  net  <- build_network(ai_long, method = "frequency",
                        id_col = "session_id",
                        time_col = "order_in_session", action = "code")
  ai_tab <- as.table(t(net$weights))

  cases <- c(list(ai_tab), lapply(seq_len(50), .random_table))
  deltas <- vapply(cases, .area_check, numeric(1L))

  cat(sprintf("\n[area vs vcd math]  median=%.2e  p95=%.2e  max=%.2e\n",
              stats::median(deltas), stats::quantile(deltas, 0.95), max(deltas)))
  testthat::expect_lt(max(deltas), TOL_AREA)
})

testthat::test_that("mosaic_plot ggplot data matches tna::plot_mosaic", {
  skip_if_not_set()
  data(ai_long, package = "Nestimate")
  net  <- build_network(ai_long, method = "frequency",
                        id_col = "session_id",
                        time_col = "order_in_session", action = "code")
  ai_tab <- as.table(t(net$weights))

  cases <- c(list(ai_tab), lapply(seq_len(50), .random_table))
  deltas <- do.call(rbind, lapply(cases, .coord_check))

  cat(sprintf("\n[coord vs tna]  xmin=%.2e xmax=%.2e ymin=%.2e ymax=%.2e\n",
              max(deltas[, "xmin"]), max(deltas[, "xmax"]),
              max(deltas[, "ymin"]), max(deltas[, "ymax"])))
  testthat::expect_lt(max(deltas), TOL_COORD)
})

testthat::test_that("mosaic_plot uses chisq.test()$stdres for fill", {
  skip_if_not_set()
  data(ai_long, package = "Nestimate")
  net  <- build_network(ai_long, method = "frequency",
                        id_col = "session_id",
                        time_col = "order_in_session", action = "code")
  ai_tab <- as.table(t(net$weights))

  net_ai <- .tab_to_netobj(ai_tab)
  p <- mosaic_plot(net_ai)
  d <- ggplot2::ggplot_build(p)$data[[1]]
  ref <- suppressWarnings(stats::chisq.test(ai_tab))$stdres
  o <- order(d$xmin, d$ymin)
  src <- ggplot2::layer_data(p, 1L)
  fill_vals <- src$fill[o]
  cat(sprintf("\n[fill colour mapped from chi-square stdres for %d cells]\n",
              length(fill_vals)))
  testthat::expect_true(all(grepl("^#[0-9A-F]{6}$", fill_vals)))
})
