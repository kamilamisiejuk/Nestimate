# =========================================================================
# Regression tests for FIX-3 (audit A04 — simplicial / hypergraph)
#
# Confirmed findings fixed:
#   A04-F01  build_simplicial(type="vr") was byte-identical to "clique"
#            (no real Vietoris-Rips filtration exists). Honest-enum decision:
#            type="vr"/"rips" now stop() with a clear message instead of
#            silently aliasing "clique" while print()/labels claim "vr".
#   A04-F02  build_hypergraph(method="vr") inherited the same dead alias.
#            Same honest fix (errors before reaching build_simplicial()).
#   A04-F03  build_simplicial(max_dim<0) silently behaved like the default
#            (igraph treats max=0 as unbounded). Now stopifnot() guarded.
#   A04-F04  persistent_homology(n_steps=0) silently returned a degenerate
#            object. Now stopifnot(n_steps >= 1L).
#   A04-F05  build_simplicial(threshold<0) silently accepted. Now guarded.
#
# Data: a netobject derived from bundled group_regulation_long, plus
# simulate_hypergraph_incidence() for hypergraph-shaped inputs.
# =========================================================================

# --- shared fixtures ----------------------------------------------------

.fix3_net <- function() {
  e <- new.env()
  utils::data("group_regulation_long", package = "Nestimate", envir = e)
  build_network(e$group_regulation_long, actor = "Actor", time = "Time",
                action = "Action", method = "relative")
}

.fix3_incidence_adj <- function(seed = 42L) {
  hgi <- simulate_hypergraph_incidence(n_nodes = 10, n_edges = 14,
                                       edge_size_range = c(2L, 3L),
                                       density = 0.85, seed = seed)
  inc <- hgi$incidence
  # node x node co-membership adjacency from the incidence matrix
  adj <- inc %*% t(inc)
  diag(adj) <- 0
  adj[adj > 0] <- 1
  storage.mode(adj) <- "double"
  adj
}

# =========================================================================
# A04-F01 — build_simplicial(type="vr") must error, not alias "clique"
# =========================================================================

test_that("build_simplicial(type='vr') errors and is NOT silently clique", {
  net <- .fix3_net()

  # Old buggy behavior: identical to clique with $type == "vr".
  # Corrected behavior: a clean, informative stop().
  expect_error(
    build_simplicial(net, type = "vr", threshold = 0.05),
    "not implemented",
    fixed = TRUE
  )
  expect_error(
    build_simplicial(net, type = "rips", threshold = 0.05),
    "not implemented",
    fixed = TRUE
  )

  # The error message must not be the cryptic match.arg() one.
  msg <- tryCatch(build_simplicial(net, type = "vr"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("'arg' should be one of", msg))
  expect_true(grepl("Vietoris-Rips", msg, fixed = TRUE))

  # clique still works and is honestly labelled.
  sc <- build_simplicial(net, type = "clique", threshold = 0.05)
  expect_s3_class(sc, "simplicial_complex")
  expect_identical(sc$type, "clique")
  # No simplicial_complex can ever carry the dead "vr" type tag now.
  expect_false(identical(sc$type, "vr"))
})

test_that("build_simplicial(type='vr') errors on a bare matrix too", {
  adj <- .fix3_incidence_adj()
  expect_error(
    build_simplicial(adj, type = "vr", threshold = 0),
    "not implemented",
    fixed = TRUE
  )
  # clique path on the same matrix is unaffected.
  sc <- build_simplicial(adj, type = "clique", threshold = 0)
  expect_s3_class(sc, "simplicial_complex")
  expect_gt(sc$n_simplices, sc$n_nodes)
})

test_that("print.simplicial_complex can never produce a VR label", {
  net <- .fix3_net()
  sc <- build_simplicial(net, type = "clique", threshold = 0.05)
  out <- capture.output(print(sc))
  expect_true(any(grepl("Clique Complex", out)))
  expect_false(any(grepl("Vietoris-Rips", out)))
})

# =========================================================================
# A04-F02 — build_hypergraph(method="vr") inherits the honest fix
# =========================================================================

test_that("build_hypergraph(method='vr') errors and is NOT silently clique", {
  net <- .fix3_net()

  expect_error(
    build_hypergraph(net, method = "vr", max_size = 3L),
    "not implemented",
    fixed = TRUE
  )
  expect_error(
    build_hypergraph(net, method = "rips", max_size = 3L),
    "not implemented",
    fixed = TRUE
  )

  msg <- tryCatch(build_hypergraph(net, method = "vr"),
                  error = function(e) conditionMessage(e))
  expect_false(grepl("'arg' should be one of", msg))
  expect_true(grepl("Vietoris-Rips", msg, fixed = TRUE))

  # clique still works.
  hg <- build_hypergraph(net, method = "clique", max_size = 3L)
  expect_s3_class(hg, "net_hypergraph")
  expect_identical(hg$params$method, "clique")
})

test_that("build_hypergraph(method='vr') errors on incidence-derived adj", {
  adj <- .fix3_incidence_adj(seed = 7L)
  expect_error(
    build_hypergraph(adj, method = "vr", max_size = 4L, threshold = 0),
    "not implemented",
    fixed = TRUE
  )
  hg <- build_hypergraph(adj, method = "clique", max_size = 4L,
                         threshold = 0)
  expect_s3_class(hg, "net_hypergraph")
  expect_gt(hg$n_hyperedges, 0L)
})

# =========================================================================
# A04-F03 — build_simplicial(max_dim < 0) must error (was silent default)
# =========================================================================

test_that("build_simplicial rejects negative / non-integer max_dim", {
  net <- .fix3_net()

  expect_error(build_simplicial(net, max_dim = -1L), "max_dim")
  expect_error(build_simplicial(net, max_dim = -5L), "max_dim")
  expect_error(build_simplicial(net, max_dim = 2.5), "max_dim")
  expect_error(build_simplicial(net, max_dim = c(1L, 2L)), "max_dim")
  expect_error(build_simplicial(net, max_dim = NA_integer_), "max_dim")

  # Valid non-negative integers still work and still cap dimension.
  sc0 <- build_simplicial(net, type = "clique", threshold = 0,
                          max_dim = 0L)
  expect_identical(sc0$dimension, 0L)
  sc2 <- build_simplicial(net, type = "clique", threshold = 0,
                          max_dim = 2L)
  expect_lte(sc2$dimension, 2L)
})

# =========================================================================
# A04-F04 — persistent_homology(n_steps < 1) must error (was silent)
# =========================================================================

test_that("persistent_homology rejects n_steps < 1 and bad max_dim", {
  net <- .fix3_net()

  expect_error(persistent_homology(net, n_steps = 0L), "n_steps")
  expect_error(persistent_homology(net, n_steps = -3L), "n_steps")
  expect_error(persistent_homology(net, n_steps = 2.5), "n_steps")
  expect_error(persistent_homology(net, n_steps = NA_integer_), "n_steps")
  expect_error(persistent_homology(net, n_steps = 5L, max_dim = -1L),
               "max_dim")

  # Valid call still produces the documented shape.
  ph <- persistent_homology(net, n_steps = 5L, max_dim = 2L)
  expect_s3_class(ph, "persistent_homology")
  expect_length(ph$thresholds, 5L)
  expect_true(all(c("threshold", "dimension", "betti") %in%
                    names(ph$betti_curve)))
})

# =========================================================================
# A04-F05 — build_simplicial(threshold < 0) must error (was silent)
# =========================================================================

test_that("build_simplicial rejects negative / malformed threshold", {
  net <- .fix3_net()

  expect_error(build_simplicial(net, threshold = -5), "threshold")
  expect_error(build_simplicial(net, threshold = -0.001), "threshold")
  expect_error(build_simplicial(net, threshold = c(0, 1)), "threshold")
  expect_error(build_simplicial(net, threshold = NA_real_), "threshold")

  # threshold = 0 (documented default domain boundary) still works.
  sc <- build_simplicial(net, type = "clique", threshold = 0)
  expect_s3_class(sc, "simplicial_complex")
  expect_gt(sc$n_simplices, 0L)
})
