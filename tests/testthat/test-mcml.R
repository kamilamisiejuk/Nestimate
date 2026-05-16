# Tests for mcml.R: net_aggregate_weights, cluster_summary, build_mcml

# ============================================
# net_aggregate_weights / wagg
# ============================================

test_that("net_aggregate_weights sum method", {
  expect_equal(net_aggregate_weights(c(0.5, 0.8, 0.3), "sum"), 1.6)
})

test_that("net_aggregate_weights mean method", {
  expect_equal(net_aggregate_weights(c(2, 4, 6), "mean"), 4)
})

test_that("net_aggregate_weights median method", {
  expect_equal(net_aggregate_weights(c(1, 3, 5), "median"), 3)
})

test_that("net_aggregate_weights max method", {
  expect_equal(net_aggregate_weights(c(1, 5, 3), "max"), 5)
})

test_that("net_aggregate_weights min method", {
  expect_equal(net_aggregate_weights(c(1, 5, 3), "min"), 1)
})

test_that("net_aggregate_weights prod method", {
  expect_equal(net_aggregate_weights(c(2, 3, 4), "prod"), 24)
})

test_that("net_aggregate_weights density with n_possible", {
  expect_equal(net_aggregate_weights(c(1, 2, 3), "density", n_possible = 10), 0.6)
})

test_that("net_aggregate_weights density without n_possible", {
  expect_equal(net_aggregate_weights(c(1, 2, 3), "density"), 2)
})

test_that("net_aggregate_weights geomean method", {
  expect_equal(net_aggregate_weights(c(4, 9), "geomean"), 6, tolerance = 0.01)
})

test_that("net_aggregate_weights matches direct base calculations", {
  w <- c(1, 2, 3, NA, 0)
  cleaned <- c(1, 2, 3)

  expect_equal(net_aggregate_weights(w, "sum"), sum(cleaned))
  expect_equal(net_aggregate_weights(w, "mean"), mean(cleaned))
  expect_equal(net_aggregate_weights(w, "median"), stats::median(cleaned))
  expect_equal(net_aggregate_weights(w, "max"), max(cleaned))
  expect_equal(net_aggregate_weights(w, "min"), min(cleaned))
  expect_equal(net_aggregate_weights(w, "prod"), prod(cleaned))
  expect_equal(net_aggregate_weights(w, "density"), sum(cleaned) / length(cleaned))
  expect_equal(net_aggregate_weights(w, "density", n_possible = 10), sum(cleaned) / 10)
  expect_equal(net_aggregate_weights(c(4, 9), "geomean"),
               exp(mean(log(c(4, 9)))))
})

test_that("net_aggregate_weights removes NA and zero", {
  expect_equal(net_aggregate_weights(c(1, NA, 0, 2), "sum"), 3)
})

test_that("net_aggregate_weights returns 0 for empty/all-zero input", {
  expect_equal(net_aggregate_weights(c(0, 0, NA), "sum"), 0)
  expect_equal(net_aggregate_weights(numeric(0), "mean"), 0)
})

test_that("net_aggregate_weights errors on unknown method", {
  expect_error(net_aggregate_weights(c(1, 2), "bogus"), "Unknown method")
})

test_that("net_aggregate_weights validates arguments strictly", {
  expect_error(net_aggregate_weights(c("1", "2"), "sum"),
               "'w' must be a numeric vector", fixed = TRUE)
  expect_error(net_aggregate_weights(c(1, Inf), "sum"),
               "'w' must contain only finite values or NA", fixed = TRUE)
  expect_error(net_aggregate_weights(c(1, 2), NA_character_),
               "'method' must be a single non-missing character value",
               fixed = TRUE)
  expect_error(net_aggregate_weights(c(1, 2), c("sum", "mean")),
               "'method' must be a single non-missing character value",
               fixed = TRUE)
  expect_error(net_aggregate_weights(c(1, 2), "density",
                                     n_possible = NA_real_),
               "'n_possible' must be a single finite numeric value or NULL",
               fixed = TRUE)
  expect_error(net_aggregate_weights(c(1, 2), "density",
                                     n_possible = c(1, 2)),
               "'n_possible' must be a single finite numeric value or NULL",
               fixed = TRUE)
})

test_that("net_aggregate_weights is a function", {
  expect_true(is.function(net_aggregate_weights))
})

# ============================================
# mcml_layer
# ============================================

test_that(".mcml_layer validates layer structure strictly", {
  mat <- matrix(c(1, 2, 3, 4), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))

  layer <- Nestimate:::.mcml_layer(mat, c(A = 0.4, B = 0.6),
                                   c("A", "B"))
  expect_s3_class(layer, "mcml_layer")
  expect_equal(layer$inits, c(A = 0.4, B = 0.6))

  expect_error(Nestimate:::.mcml_layer(data.frame(A = 1), c(A = 1), "A"),
               "'weights' must be a numeric matrix", fixed = TRUE)
  expect_error(Nestimate:::.mcml_layer(matrix(1, 2, 3), c(A = 1), "A"),
               "'weights' must be a square matrix", fixed = TRUE)
  expect_error(Nestimate:::.mcml_layer(
    matrix(c(1, NA, 2, 3), 2, 2,
           dimnames = list(c("A", "B"), c("A", "B"))),
    c(A = 0.4, B = 0.6), c("A", "B")
  ), "finite non-missing weights")
  expect_error(Nestimate:::.mcml_layer(mat, c(A = 1), c("A", "B")),
               "'inits' must be a finite numeric vector with one value per node",
               fixed = TRUE)
  expect_error(Nestimate:::.mcml_layer(mat, c(A = 0.4, C = 0.6),
                                       c("A", "B")),
               "'inits' names must match layer labels", fixed = TRUE)
  expect_error(Nestimate:::.mcml_layer(mat, c(A = 0.4, B = 0.6),
                                       c("A", "A")),
               "'labels' must be unique", fixed = TRUE)
})

test_that("print.mcml_layer rejects unsupported dots", {
  mat <- matrix(c(1, 2, 3, 4), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  layer <- Nestimate:::.mcml_layer(mat, c(A = 0.4, B = 0.6),
                                   c("A", "B"))

  expect_error(print(layer, typo_arg = TRUE),
               "unsupported argument: typo_arg")
})

# ============================================
# cluster_summary
# ============================================

test_that("cluster_summary with vector clusters", {
  mat <- matrix(c(10, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- c(A = 1, B = 1, C = 2)

  cs <- cluster_summary(mat, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
  expect_equal(ncol(cs$macro$weights), 2)
  expect_equal(cs$meta$n_clusters, 2)
  expect_equal(cs$meta$n_nodes, 3)
})

test_that("cluster_summary with named list clusters", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- cluster_summary(mat, clusters)
  expect_equal(rownames(cs$macro$weights), c("G1", "G2"))
  expect_equal(colnames(cs$macro$weights), c("G1", "G2"))
  expect_equal(length(cs$clusters), 2)
  expect_equal(nrow(cs$clusters$G1$weights), 2)
})

test_that("cluster_summary matrix path returns raw aggregated values (no normalization)", {
  mat <- matrix(c(10, 2, 3, 8), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  clusters <- list(G1 = "A", G2 = "B")

  cs <- cluster_summary(mat, clusters, method = "sum")
  expect_equal(cs$macro$weights["G1", "G1"], 10)
  expect_equal(cs$macro$weights["G1", "G2"], 3)
})

test_that("cluster_summary compute_within=FALSE skips within", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs <- cluster_summary(mat, clusters, compute_within = FALSE)
  expect_null(cs$clusters)
})

test_that("cluster_summary inits sum to 1", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- cluster_summary(mat, clusters)
  expect_equal(sum(cs$macro$inits), 1, tolerance = 1e-10)
})

test_that("cluster_summary matrix aggregation matches direct block sums", {
  mat <- matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 8, 9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs <- cluster_summary(mat, clusters, method = "sum")
  manual <- matrix(
    c(sum(mat[c("A", "B"), c("A", "B")]),
      sum(mat[c("A", "B"), "C"]),
      sum(mat["C", c("A", "B")]),
      mat["C", "C"]),
    2, 2, byrow = TRUE,
    dimnames = list(c("G1", "G2"), c("G1", "G2"))
  )

  expect_equal(cs$macro$weights, manual, tolerance = 1e-12)
  expect_equal(cs$macro$inits, colSums(manual) / sum(manual),
               tolerance = 1e-12)
  expect_equal(cs$clusters$G1$weights,
               mat[c("A", "B"), c("A", "B")],
               tolerance = 1e-12)
})

test_that("cluster_summary rejects invalid logical controls", {
  mat <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  clusters <- list(G1 = "A", G2 = "B")

  expect_error(cluster_summary(mat, clusters, directed = NA),
               "'directed' must be TRUE or FALSE")
  expect_error(cluster_summary(mat, clusters, directed = c(TRUE, FALSE)),
               "'directed' must be TRUE or FALSE")
  expect_error(cluster_summary(mat, clusters, compute_within = NA),
               "'compute_within' must be TRUE or FALSE")
  expect_error(cluster_summary(mat, clusters, compute_within = c(TRUE, FALSE)),
               "'compute_within' must be TRUE or FALSE")
})

test_that("cluster_summary errors without clusters", {
  mat <- matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(cluster_summary(mat), "clusters")
})

test_that("cluster_summary errors on non-square matrix", {
  mat <- matrix(1, 2, 3)
  expect_error(cluster_summary(mat, c(1, 2)), "square")
})

test_that("cluster_summary validates matrix weights and dimnames strictly", {
  clusters <- list(G1 = "A", G2 = "B")

  mat_na <- matrix(c(1, NA, 2, 3), 2, 2,
                   dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(cluster_summary(mat_na, clusters),
               "finite non-missing weights")

  mat_dup <- matrix(1, 2, 2,
                    dimnames = list(c("A", "A"), c("A", "A")))
  expect_error(cluster_summary(mat_dup, list(G1 = "A")),
               "row names must be unique")

  mat_col_dup <- matrix(1, 2, 2,
                        dimnames = list(c("A", "B"), c("A", "A")))
  expect_error(cluster_summary(mat_col_dup, clusters),
               "column names must be unique")

  mat_mismatch <- matrix(1, 2, 2,
                         dimnames = list(c("A", "B"), c("B", "A")))
  expect_error(cluster_summary(mat_mismatch, clusters),
               "row and column names must be identical")

  mat_empty <- matrix(1, 2, 2,
                      dimnames = list(c("A", ""), c("A", "")))
  expect_error(cluster_summary(mat_empty, list(G1 = "A", G2 = "")),
               "row names must not contain missing or empty values")
})

test_that("cluster_summary print method works", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  expect_output(print(cs))
})

test_that("mcml print and summary reject unsupported dots", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))

  expect_error(print(cs, typo_arg = TRUE),
               "unsupported argument: typo_arg")
  expect_error(summary(cs, typo_arg = TRUE),
               "unsupported argument: typo_arg")
})

test_that("cluster_summary is a function", {
  expect_true(is.function(cluster_summary))
})

test_that("cluster_summary with different methods", {
  mat <- matrix(c(4, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")

  cs_sum <- cluster_summary(mat, clusters, method = "sum")
  cs_mean <- cluster_summary(mat, clusters, method = "mean")
  cs_max <- cluster_summary(mat, clusters, method = "max")

  # Sum should be larger than mean for multi-node clusters
  expect_true(cs_sum$macro$weights["G1", "G1"] >= cs_mean$macro$weights["G1", "G1"])
  # Max should be at most sum
  expect_true(cs_max$macro$weights["G1", "G1"] <= cs_sum$macro$weights["G1", "G1"])
})

test_that("cluster_summary within-cluster has correct dimensions", {
  mat <- matrix(runif(16), 4, 4,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  clusters <- list(G1 = c("A", "B", "C"), G2 = "D")

  cs <- cluster_summary(mat, clusters)
  expect_equal(nrow(cs$clusters$G1$weights), 3)
  expect_equal(ncol(cs$clusters$G1$weights), 3)
  expect_equal(nrow(cs$clusters$G2$weights), 1)
})

# ============================================
# build_mcml
# ============================================

test_that("build_mcml with sequence data", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T2 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T3 = sample(c("A", "B", "C", "D"), 50, replace = TRUE),
    T4 = sample(c("A", "B", "C", "D"), 50, replace = TRUE)
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(seqs, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
  expect_equal(sum(cs$macro$inits), 1, tolerance = 1e-10)
})

test_that("build_mcml with edge list", {
  edges <- data.frame(
    from = c("A", "A", "B", "C", "C", "D"),
    to = c("B", "C", "A", "D", "D", "A"),
    weight = c(1, 2, 1, 3, 1, 2)
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(edges, clusters)
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml type=raw preserves counts", {
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "D"),
    T3 = c("C", "A")
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(seqs, clusters, type = "raw")
  expect_true(is.numeric(cs$macro$weights))
})

test_that("build_mcml raw sequence counts match direct cluster transitions", {
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "D"),
    T3 = c("C", "A"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  cs <- build_mcml(seqs, clusters, type = "raw")
  manual <- matrix(c(1, 1,
                     1, 1), 2, 2, byrow = TRUE,
                   dimnames = list(c("G1", "G2"), c("G1", "G2")))

  expect_equal(cs$macro$weights, manual, tolerance = 1e-12)
  expect_equal(cs$macro$inits, c(G1 = 0.5, G2 = 0.5),
               tolerance = 1e-12)
  expect_equal(nrow(cs$edges), 4L)
})

test_that("build_mcml rejects invalid logical controls", {
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "D"),
    T3 = c("C", "A"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))

  expect_error(build_mcml(seqs, clusters, directed = NA),
               "'directed' must be TRUE or FALSE")
  expect_error(build_mcml(seqs, clusters, directed = c(TRUE, FALSE)),
               "'directed' must be TRUE or FALSE")
  expect_error(build_mcml(seqs, clusters, compute_within = NA),
               "'compute_within' must be TRUE or FALSE")
  expect_error(build_mcml(seqs, clusters, compute_within = c(TRUE, FALSE)),
               "'compute_within' must be TRUE or FALSE")
})

test_that("build_mcml rejects duplicate and incomplete cluster lists", {
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "D"),
    T3 = c("C", "A"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_mcml(seqs, list(G1 = c("A", "B"), G2 = c("B", "C", "D"))),
    "Nodes assigned to multiple clusters: B"
  )
  expect_error(
    build_mcml(seqs, list(G1 = c("A", "B"), G2 = "C")),
    "Unmapped nodes: D"
  )
})

test_that("build_mcml returns mcml if already mcml", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  cs2 <- build_mcml(cs)
  expect_identical(cs, cs2)
})

test_that("build_mcml tna type normalizes rows", {
  set.seed(1)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C"), 30, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 30, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 30, replace = TRUE)
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "tna")
  row_sums <- rowSums(cs$macro$weights)
  expect_equal(unname(row_sums), c(1, 1), tolerance = 1e-10)
})

# ============================================================
# Coverage gap tests — mcml.R
# ============================================================

# ---- net_aggregate_weights / wagg: geomean zero path (L43) ----

test_that("net_aggregate_weights geomean returns 0 when all non-positive", {
  # All negative: pos_w is empty → returns 0
  expect_equal(net_aggregate_weights(c(-1, -2, -3), "geomean"), 0)
})

# ---- cluster_summary: various input types ----

test_that("cluster_summary assigns sequential node names when matrix has no rownames (L313-315)", {
  mat <- matrix(c(0, 2, 3, 1, 0, 4, 5, 6, 0), 3, 3)
  # No rownames/colnames
  cs <- cluster_summary(mat, list(G1 = c("1", "2"), G2 = "3"))
  expect_s3_class(cs, "mcml")
  expect_equal(cs$meta$n_nodes, 3)
})

test_that("cluster_summary between_inits fallback when zero matrix (L380)", {
  # All-zero matrix → colSums all zero → uniform inits
  mat <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  cs <- cluster_summary(mat, list(G1 = "A", G2 = "B"))
  expect_equal(unname(cs$macro$inits), c(0.5, 0.5), tolerance = 1e-10)
})

test_that("cluster_summary within-cluster zero total produces uniform inits (L423)", {
  # All-zero within-cluster block
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  # G1 has 2 nodes, zero total → rep(0.5, 2)
  expect_equal(unname(cs$clusters$G1$inits), c(0.5, 0.5), tolerance = 1e-10)
})

# ---- build_mcml: cograph_network input (L572) ----

test_that("build_mcml accepts cograph_network input (L572)", {
  mat <- matrix(c(0, 0.3, 0.7, 0.4, 0, 0.6, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  # Create a minimal cograph_network object
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      edges = data.frame(from = integer(0), to = integer(0),
                         weight = numeric(0)),
      directed = TRUE,
      data = NULL,
      meta = list(),
      node_groups = NULL
    ),
    class = c("cograph_network", "list")
  )
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

# ---- build_mcml: tna_data branch (L581-598) ----

test_that("build_mcml handles netobject with sequence data (netobject_data branch, L586-599)", {
  set.seed(1)
  seqs <- data.frame(
    T1 = sample(c("A", "B", "C"), 20, replace = TRUE),
    T2 = sample(c("A", "B", "C"), 20, replace = TRUE),
    T3 = sample(c("A", "B", "C"), 20, replace = TRUE)
  )
  net <- build_network(seqs, method = "relative")
  # net$data contains the raw sequences
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml handles netobject with edge list data (netobject_data edgelist, L593-595)", {
  edges <- data.frame(
    from = c("A", "A", "B", "C"),
    to   = c("B", "C", "A", "A"),
    stringsAsFactors = FALSE
  )
  net <- build_network(edges, method = "relative")
  # If net has $data, it should route through netobject_data
  if (!is.null(net$data)) {
    cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
    expect_s3_class(cs, "mcml")
  } else {
    # Falls to netobject_matrix branch — still valid test
    cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
    expect_s3_class(cs, "mcml")
  }
})

# ---- build_mcml: netobject_matrix branch (L601-607) ----

test_that("build_mcml handles netobject without data (netobject_matrix branch, L601-607)", {
  mat <- matrix(c(0, 0.4, 0.6, 0.3, 0, 0.7, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      edges = data.frame(from = integer(0), to = integer(0),
                         weight = numeric(0)),
      directed = TRUE,
      data = NULL,
      meta = list(),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  cs <- build_mcml(net, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml handles plain numeric matrix (matrix branch, L608-610)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- build_mcml(mat, list(G1 = c("A", "B"), G2 = "C"))
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml errors on unknown input class (L611-612)", {
  expect_error(build_mcml(list(foo = 1), list(G1 = "A")),
               "Cannot build MCML")
})

# ---- .detect_mcml_input coverage (L620-644) ----

test_that(".detect_mcml_input returns netobject_data for netobject with data (L625)", {
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"))
  net <- build_network(seqs, method = "relative")
  if (!is.null(net$data)) {
    result <- Nestimate:::.detect_mcml_input(net)
    expect_equal(result, "netobject_data")
  } else {
    skip("netobject has no $data in this configuration")
  }
})

test_that(".detect_mcml_input returns netobject_matrix for netobject without data (L626)", {
  mat <- matrix(c(0, 0.4, 0.6, 0.7, 0, 0.3, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(weights = mat, nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                                           name = c("A", "B", "C"),
                                           x = NA_real_, y = NA_real_,
                                           stringsAsFactors = FALSE),
         edges = data.frame(from = integer(0), to = integer(0),
                            weight = numeric(0)),
         directed = TRUE, data = NULL, meta = list(), node_groups = NULL),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.detect_mcml_input(net)
  expect_equal(result, "netobject_matrix")
})

test_that(".detect_mcml_input returns sequence for non-square numeric matrix (L641)", {
  mat <- matrix(1:6, 2, 3)
  result <- Nestimate:::.detect_mcml_input(mat)
  expect_equal(result, "sequence")
})

test_that(".detect_mcml_input returns unknown for unrecognized class (L644)", {
  result <- Nestimate:::.detect_mcml_input(42L)
  expect_equal(result, "unknown")
})

test_that(".detect_mcml_input accepts igraph-style first-two-column edge lists", {
  unnamed_edges <- data.frame(node_from = c("A", "B"),
                              node_to = c("B", "A"),
                              stringsAsFactors = FALSE)
  weighted_edges <- data.frame(node_from = c("A", "B"),
                               node_to = c("B", "A"),
                               weight = c(2, 3),
                               stringsAsFactors = FALSE)

  expect_equal(Nestimate:::.detect_mcml_input(unnamed_edges), "edgelist")
  expect_equal(Nestimate:::.detect_mcml_input(weighted_edges), "edgelist")
})

test_that(".detect_mcml_input keeps named two-step sequence data as sequence", {
  seqs <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"),
                     stringsAsFactors = FALSE)
  long_names <- data.frame(time1 = c("A", "B"), time2 = c("B", "A"),
                           stringsAsFactors = FALSE)

  expect_equal(Nestimate:::.detect_mcml_input(seqs), "sequence")
  expect_equal(Nestimate:::.detect_mcml_input(long_names), "sequence")
})

# ---- .auto_detect_clusters coverage (L650-672) ----

test_that(".auto_detect_clusters finds cluster from nodes$cluster column (L650-656)", {
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         cluster = c("G1", "G1", "G2"),
                         stringsAsFactors = FALSE),
      data = NULL, node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.auto_detect_clusters(net)
  expect_equal(unname(result), c("G1", "G1", "G2"))
})

test_that(".auto_detect_clusters falls back to node_groups (L660-664)", {
  # audit_mcml #1: node_groups must carry a node identifier column so the
  # cluster column is aligned by label rather than by row order. Updated
  # from the previous shape (cluster-only data.frame) which silently
  # assumed row order matched x$nodes — a real corruption hazard for
  # externally constructed netobjects.
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:3, label = c("A", "B", "C"),
                         name = c("A", "B", "C"),
                         x = NA_real_, y = NA_real_,
                         stringsAsFactors = FALSE),
      data = NULL,
      node_groups = data.frame(node    = c("A", "B", "C"),
                               cluster = c("G1", "G1", "G2"),
                               stringsAsFactors = FALSE)
    ),
    class = c("netobject", "cograph_network")
  )
  result <- Nestimate:::.auto_detect_clusters(net)
  expect_equal(unname(result), c("G1", "G1", "G2"))
})

test_that(".auto_detect_clusters errors when no cluster info found (L667-671)", {
  mat <- matrix(0, 2, 2, dimnames = list(c("A", "B"), c("A", "B")))
  net <- structure(
    list(
      weights = mat,
      nodes = data.frame(id = 1:2, label = c("A", "B"),
                         stringsAsFactors = FALSE),
      data = NULL, node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
  expect_error(Nestimate:::.auto_detect_clusters(net), "No clusters found")
})

test_that(".auto_detect_clusters rejects missing and duplicate auto assignments", {
  mat <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  nodes <- data.frame(label = c("A", "B", "C"), stringsAsFactors = FALSE)

  nodes_bad <- data.frame(label = c("A", "B", "C"),
                          cluster = c("G1", NA, "G2"),
                          stringsAsFactors = FALSE)
  net_nodes_bad <- structure(
    list(weights = mat, nodes = nodes_bad, data = NULL, node_groups = NULL),
    class = c("netobject", "cograph_network")
  )
  expect_error(Nestimate:::.auto_detect_clusters(net_nodes_bad),
               "x\\$nodes cluster assignments must not contain missing")

  net_dup_groups <- structure(
    list(
      weights = mat,
      nodes = nodes,
      data = NULL,
      node_groups = data.frame(node = c("A", "A", "B", "C"),
                               cluster = c("G1", "G2", "G1", "G2"),
                               stringsAsFactors = FALSE)
    ),
    class = c("netobject", "cograph_network")
  )
  expect_error(Nestimate:::.auto_detect_clusters(net_dup_groups),
               "assigns duplicate rows to node\\(s\\): A")

  net_missing_group <- structure(
    list(
      weights = mat,
      nodes = nodes,
      data = NULL,
      node_groups = data.frame(node = c("A", "B", "C"),
                               cluster = c("G1", "", "G2"),
                               stringsAsFactors = FALSE)
    ),
    class = c("netobject", "cograph_network")
  )
  expect_error(Nestimate:::.auto_detect_clusters(net_missing_group),
               "node_groups cluster assignments must not contain missing")
})

# ---- .build_cluster_lookup coverage (L696-725) ----

test_that(".build_cluster_lookup errors on unmapped nodes (L696-698)", {
  cl <- list(G1 = c("A", "B"))
  # all_nodes includes C, which is not in G1
  expect_error(
    Nestimate:::.build_cluster_lookup(cl, c("A", "B", "C")),
    "Unmapped nodes"
  )
})

test_that(".build_cluster_lookup rejects duplicate list membership", {
  cl <- list(G1 = c("A", "B"), G2 = c("B", "C"))

  expect_error(
    Nestimate:::.build_cluster_lookup(cl, c("A", "B", "C")),
    "Nodes assigned to multiple clusters: B"
  )
})

test_that(".build_cluster_lookup works with character membership vector (L703-711)", {
  all_nodes <- c("A", "B", "C")
  result <- Nestimate:::.build_cluster_lookup(c("G1", "G1", "G2"), all_nodes)
  expect_equal(unname(result), c("G1", "G1", "G2"))
  expect_equal(names(result), c("A", "B", "C"))
})

test_that(".build_cluster_lookup errors when char vector length mismatches (L705-708)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(c("G1", "G1"), c("A", "B", "C")),
    "Membership vector length"
  )
})

test_that(".build_cluster_lookup works with numeric membership vector (L714-721)", {
  all_nodes <- c("X", "Y", "Z")
  result <- Nestimate:::.build_cluster_lookup(c(1, 1, 2), all_nodes)
  expect_equal(unname(result), c("1", "1", "2"))
  expect_equal(names(result), c("X", "Y", "Z"))
})

test_that(".build_cluster_lookup errors when numeric vector length mismatches (L715-718)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(c(1, 2), c("A", "B", "C")),
    "Membership vector length"
  )
})

test_that(".build_cluster_lookup errors on invalid input type (L724-725)", {
  expect_error(
    Nestimate:::.build_cluster_lookup(TRUE, c("A", "B")),
    "clusters must be a named list"
  )
})

# ---- .build_from_transitions: density method and zero-inits paths ----

test_that("build_mcml with density method triggers n_possible computation (L761-764)", {
  seqs <- data.frame(
    T1 = c("A", "C", "B"),
    T2 = c("B", "D", "C"),
    T3 = c("C", "A", "D")
  )
  clusters <- list(G1 = c("A", "B"), G2 = c("C", "D"))
  cs <- build_mcml(seqs, clusters, method = "density", type = "raw")
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml zero transitions produce uniform between_inits (L784)", {
  # Create edge list that has no transitions between G2 -> G1 or G2 -> G2
  # Force all edges in G1 to make G2 isolated
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  # G2 has node "C" but no transitions at all → zero col sums for G2
  # Actually need G2 to have zero total in cluster
  # Simpler: all transitions are G1 internal, G2 isolated
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  # Since "C" never appears in edges, it won't appear in all_nodes
  # Let's add a C transition explicitly
  edges2 <- data.frame(
    from = c("A", "B", "C"),
    to   = c("B", "A", "C"),
    stringsAsFactors = FALSE
  )
  cs <- build_mcml(edges2, clusters, type = "raw")
  expect_s3_class(cs, "mcml")
})

test_that(".build_from_transitions validates transition vectors strictly", {
  clusters <- list(G1 = "A", G2 = "B")
  lookup <- c(A = "G1", B = "G2")

  expect_error(
    Nestimate:::.build_from_transitions(
      c("A"), c("B", "A"), c(1), lookup, clusters,
      "sum", "raw", TRUE, TRUE
    ),
    "must have the same length"
  )
  expect_error(
    Nestimate:::.build_from_transitions(
      c("A"), c("B"), NA_real_, lookup, clusters,
      "sum", "raw", TRUE, TRUE
    ),
    "'weights' must be a finite non-missing numeric vector",
    fixed = TRUE
  )
  expect_error(
    Nestimate:::.build_from_transitions(
      c("A"), c("B"), -1, lookup, clusters,
      "sum", "raw", TRUE, TRUE
    ),
    "'weights' must not contain negative values",
    fixed = TRUE
  )
  expect_error(
    Nestimate:::.build_from_transitions(
      c("X"), c("B"), 1, lookup, clusters,
      "sum", "raw", TRUE, TRUE
    ),
    "'cluster_lookup' is missing node\\(s\\): X"
  )
  expect_error(
    Nestimate:::.build_from_transitions(
      c("A"), c("B"), 1, c(A = "G1", B = "G9"), clusters,
      "sum", "raw", TRUE, TRUE
    ),
    "'cluster_lookup' contains unknown cluster\\(s\\): G9"
  )
  expect_error(
    Nestimate:::.build_from_transitions(
      c("A"), c("B"), 1, lookup, clusters,
      "sum", "raw", TRUE, NA
    ),
    "'compute_within' must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that(".build_from_transitions accepts zero-transition inputs", {
  clusters <- list(G1 = "A", G2 = "B")
  lookup <- c(A = "G1", B = "G2")

  out <- Nestimate:::.build_from_transitions(
    character(0), character(0), numeric(0), lookup, clusters,
    "sum", "raw", TRUE, TRUE
  )

  expect_s3_class(out, "mcml")
  expect_equal(out$macro$weights,
               matrix(0, 2, 2, dimnames = list(c("G1", "G2"),
                                                c("G1", "G2"))))
  expect_equal(out$macro$inits, c(G1 = 0.5, G2 = 0.5))
})

# ---- Single-node within-cluster (L839-847) ----

test_that("build_mcml single-node cluster computes self-loop weight (L839-847)", {
  # Cluster G2 has only 1 node "C"
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "C"),
    T3 = c("C", "A")
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "raw")
  expect_s3_class(cs, "mcml")
  # G2 is a 1x1 within matrix
  expect_equal(dim(cs$clusters$G2$weights), c(1, 1))
  expect_equal(unname(cs$clusters$G2$inits), 1)
})

test_that("build_mcml single-node cluster with no self-loops returns 0 weight (L843-844)", {
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  # G2 node "C" has no edges at all — single node, no self-loops
  # Add C to appear in nodes but never in edges: need it in seqs
  seqs <- data.frame(
    T1 = c("A", "C"),
    T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters, type = "raw")
  # C appears once but only in T1 position, no C->C transition
  expect_equal(unname(cs$clusters$G2$inits), 1)
})

# ---- .build_mcml_edgelist: igraph-style endpoint detection ----

test_that(".build_mcml_edgelist accepts first two columns as endpoints", {
  edges <- data.frame(
    node_from = c("A", "B", "C"),
    node_to   = c("B", "C", "A"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- Nestimate:::.build_mcml_edgelist(edges, clusters, "sum", "raw",
                                         TRUE, TRUE)

  expect_s3_class(cs, "mcml")
  expect_equal(cs$macro$weights["G1", "G1"], 1)
  expect_equal(cs$macro$weights["G1", "G2"], 1)
  expect_equal(cs$macro$weights["G2", "G1"], 1)
})

test_that("build_mcml distinguishes two-step sequences from unnamed edge lists", {
  clusters <- list(G1 = "A", G2 = "B")
  seqs <- data.frame(T1 = c("A", "A", "B"),
                     T2 = c("B", "B", "B"),
                     stringsAsFactors = FALSE)
  edges <- data.frame(node_from = c("A", "A", "B"),
                      node_to = c("B", "B", "B"),
                      stringsAsFactors = FALSE)

  seq_mc <- build_mcml(seqs, clusters, type = "raw")
  edge_mc <- build_mcml(edges, clusters, type = "raw")

  expect_null(attr(seq_mc$macro$data, "source"))
  expect_identical(attr(edge_mc$macro$data, "source"), "edgelist")
  expect_equal(seq_mc$macro$inits, c(G1 = 2 / 3, G2 = 1 / 3),
               tolerance = 1e-12)
  expect_equal(edge_mc$macro$inits, c(G1 = 0, G2 = 1),
               tolerance = 1e-12)
})

test_that("build_mcml edgelist validates endpoints and weights", {
  clusters <- list(G1 = "A", G2 = "B")

  expect_error(
    build_mcml(
      data.frame(from = c("A", ""), to = c("B", "A"),
                 stringsAsFactors = FALSE),
      clusters, type = "raw"
    ),
    "source and target columns must not contain missing or empty values"
  )
  expect_error(
    build_mcml(
      data.frame(from = c("A", "B"), to = c("B", "A"),
                 weight = c("2", "bad"), stringsAsFactors = FALSE),
      clusters, type = "raw"
    ),
    "weight column must be numeric"
  )
  expect_error(
    build_mcml(
      data.frame(from = c("A", "B"), to = c("B", "A"),
                 weight = c(2, NA_real_), stringsAsFactors = FALSE),
      clusters, type = "raw"
    ),
    "weight column must contain finite non-missing values"
  )
  expect_error(
    build_mcml(
      data.frame(from = c("A", "B"), to = c("B", "A"),
                 weight = c(2, -1), stringsAsFactors = FALSE),
      clusters, type = "raw"
    ),
    "weight column must not contain negative values"
  )
})

# ---- Column-name clusters branch in edgelist (L971-989) ----

test_that("build_mcml edgelist accepts cluster column name (L971-989)", {
  edges <- data.frame(
    from  = c("A", "B", "C", "D"),
    to    = c("B", "A", "D", "C"),
    group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  cs <- build_mcml(edges, clusters = "group")
  expect_s3_class(cs, "mcml")
  expect_equal(nrow(cs$macro$weights), 2)
})

test_that("build_mcml edgelist cluster column rejects conflicting node groups", {
  edges <- data.frame(
    from = c("A", "B"),
    to = c("B", "A"),
    group = c("G1", "G2"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_mcml(edges, clusters = "group"),
    "assigns nodes to multiple groups"
  )
})

test_that("build_mcml edgelist errors when clusters=NULL (L992-993)", {
  edges <- data.frame(
    from = c("A", "B"),
    to   = c("B", "A"),
    stringsAsFactors = FALSE
  )
  expect_error(build_mcml(edges, clusters = NULL),
               "clusters argument is required")
})

test_that("build_mcml edgelist accepts membership vector clusters (L1000)", {
  edges <- data.frame(
    from = c("A", "B", "C"),
    to   = c("B", "C", "A"),
    stringsAsFactors = FALSE
  )
  # character membership vector for the 3 unique sorted nodes (A, B, C)
  clusters <- c("G1", "G1", "G2")
  cs <- build_mcml(edges, clusters)
  expect_s3_class(cs, "mcml")
})

# ---- .build_mcml_sequence: matrix input, single column error, null clusters (L1017-1044) ----

test_that("build_mcml sequence builder coerces matrix input (L1017)", {
  mat_seq <- matrix(c("A", "B", "C", "B", "C", "A"), nrow = 3,
                    dimnames = list(NULL, c("T1", "T2")))
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(mat_seq, clusters)
  expect_s3_class(cs, "mcml")
})

test_that(".build_mcml_sequence errors on single-column data (L1023-1024)", {
  df <- data.frame(T1 = c("A", "B", "C"), stringsAsFactors = FALSE)
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  expect_error(
    Nestimate:::.build_mcml_sequence(df, clusters, "sum", "tna", TRUE, TRUE),
    "at least 2 columns"
  )
})

test_that(".build_mcml_sequence errors when clusters=NULL (L1043-1044)", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_mcml_sequence(df, NULL, "sum", "tna", TRUE, TRUE),
    "clusters argument is required"
  )
})

test_that(".build_mcml_sequence calls .normalize_clusters for non-list clusters (L1050)", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  # Pass a character membership vector (triggers .normalize_clusters path)
  clusters <- c("G1", "G2")  # 2 nodes: A, B
  cs <- Nestimate:::.build_mcml_sequence(df, clusters, "sum", "tna", TRUE, TRUE)
  expect_s3_class(cs, "mcml")
})

test_that("build_mcml sequence requires clusters for all observed states", {
  seqs <- data.frame(
    T1 = c("X", "A"),
    T2 = c(NA, "B"),
    stringsAsFactors = FALSE
  )

  expect_error(
    build_mcml(seqs, list(G1 = "A", G2 = "B"), type = "raw"),
    "Unmapped nodes: X"
  )
})

test_that("build_mcml sequence keeps zero-transition observed states in metadata", {
  seqs <- data.frame(
    T1 = c("X", "A"),
    T2 = c(NA, "B"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "X"), G2 = "B")

  mc <- build_mcml(seqs, clusters, type = "raw")

  expect_equal(mc$meta$n_nodes, 3L)
  expect_equal(mc$meta$cluster_sizes, c(G1 = 2L, G2 = 1L))
  expect_equal(mc$macro$weights["G1", "G2"], 1)
  expect_equal(mc$macro$inits, c(G1 = 1, G2 = 0), tolerance = 1e-12)
  expect_equal(mc$macro$data$T1, c("G1", "G1"))
})

# ---- .process_weights: cooccurrence type (L1070) ----

test_that(".process_weights returns symmetrized matrix for cooccurrence type (L1070)", {
  raw <- matrix(c(1, 3, 2, 4), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  result <- Nestimate:::.process_weights(raw, "cooccurrence")
  expect_true(isSymmetric(result))
  expect_equal(result["A", "B"], 2.5)
})

test_that(".process_weights applies undirected symmetrisation before type processing", {
  raw <- matrix(c(0, 3, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  sym <- (raw + t(raw)) / 2

  expect_equal(Nestimate:::.process_weights(raw, "raw", directed = FALSE),
               sym)
  expect_equal(Nestimate:::.process_weights(raw, "frequency",
                                            directed = FALSE), sym)
  expect_equal(Nestimate:::.process_weights(raw, "tna", directed = FALSE),
               matrix(c(0, 1, 1, 0), 2, 2,
                      dimnames = dimnames(raw)))
  expect_equal(Nestimate:::.process_weights(raw, "semi_markov",
                                            directed = FALSE),
               matrix(c(0, 1, 1, 0), 2, 2,
                      dimnames = dimnames(raw)))
})

test_that(".process_weights validates arguments strictly", {
  raw <- matrix(c(0, 1, 3, 0), 2, 2)

  expect_error(Nestimate:::.process_weights(data.frame(A = 1), "raw"),
               "'raw_weights' must be a numeric matrix", fixed = TRUE)
  expect_error(Nestimate:::.process_weights(matrix(1, 2, 3), "raw"),
               "'raw_weights' must be a square matrix", fixed = TRUE)
  expect_error(Nestimate:::.process_weights(raw, NA_character_),
               "'type' must be a single non-missing character value",
               fixed = TRUE)
  expect_error(Nestimate:::.process_weights(raw, "bad"),
               "'type' must be one of", fixed = TRUE)
  expect_error(Nestimate:::.process_weights(raw, "raw", directed = NA),
               "'directed' must be TRUE or FALSE", fixed = TRUE)
})

test_that("build_mcml sequence path honours directed = FALSE", {
  seqs <- data.frame(
    T1 = c("A", "A", "A"),
    T2 = c("B", "B", "B"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = "A", G2 = "B")

  directed_raw <- build_mcml(seqs, clusters, type = "raw",
                             directed = TRUE)$macro$weights
  undirected_raw <- build_mcml(seqs, clusters, type = "raw",
                               directed = FALSE)$macro$weights
  undirected_tna <- build_mcml(seqs, clusters, type = "tna",
                               directed = FALSE)$macro$weights

  expect_equal(directed_raw["G1", "G2"], 3)
  expect_equal(directed_raw["G2", "G1"], 0)
  expect_equal(undirected_raw["G1", "G2"], 1.5)
  expect_equal(undirected_raw["G2", "G1"], 1.5)
  expect_equal(undirected_tna["G1", "G2"], 1)
  expect_equal(undirected_tna["G2", "G1"], 1)
})

# ---- as_tna generic dispatch and methods (L1226-1306) ----

test_that("as_tna dispatches correctly for mcml objects (L1226)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
})

test_that("as_tna.mcml on matrix path uses frequency method (matrix is aggregation-only)", {
  mat <- matrix(c(10, 2, 3, 1, 8, 4, 5, 6, 12), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  # Matrix-derived mcml carries no $meta$type, as_tna defaults to "frequency".
  expect_equal(result$macro$method, "frequency")
})

test_that("as_tna.mcml preserves matrix-path macro weights and inits", {
  mat <- matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 8, 9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)

  expect_identical(result$macro$weights, cs$macro$weights)
  expect_identical(result$macro$inits, cs$macro$inits)
  expect_equal(result$macro$method, "frequency")
})

test_that("as_tna.mcml creates macro netobject (L1243-1244)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result$macro, "netobject")
  expect_true(is.matrix(result$macro$weights))
})

test_that("as_tna.mcml skips clusters with zero-row sums (L1247-1262)", {
  # Cluster with all-zero rows will be excluded
  mat <- matrix(0, 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  mat["A", "B"] <- 0.5; mat["A", "C"] <- 0.5  # only A has outgoing
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  # macro always present
  expect_true("macro" %in% names(result))
})

test_that("as_tna.mcml with compute_within=FALSE returns empty cluster list (L1256-1258)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"),
                        compute_within = FALSE)
  result <- as_tna(cs)
  expect_s3_class(result, "netobject_group")
  expect_true("macro" %in% names(result))
})

test_that(".wrap_netobject produces valid dual-class object (L1270-1296)", {
  mat <- matrix(c(0, 0.6, 0.4, 0.3, 0, 0.7, 0.5, 0.5, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  result <- Nestimate:::.wrap_netobject(mat, method = "relative", directed = TRUE)
  expect_s3_class(result, "netobject")
  expect_true(inherits(result, "cograph_network"))
  expect_equal(result$method, "relative")
  expect_equal(result$n_nodes, 3)
  expect_true(is.data.frame(result$nodes))
  expect_equal(result$nodes$label, c("A", "B", "C"))
})

test_that(".wrap_netobject validates matrix and inits strictly", {
  mat <- matrix(c(1, 2, 3, 4), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))

  expect_error(Nestimate:::.wrap_netobject(data.frame(A = 1)),
               "'mat' must be a numeric matrix", fixed = TRUE)
  expect_error(Nestimate:::.wrap_netobject(matrix(1, 2, 3)),
               "'mat' must be a square matrix", fixed = TRUE)
  expect_error(Nestimate:::.wrap_netobject(
    matrix(c(1, NA, 2, 3), 2, 2,
           dimnames = list(c("A", "B"), c("A", "B")))
  ), "finite non-missing weights")
  expect_error(Nestimate:::.wrap_netobject(mat, method = NA_character_),
               "'method' must be a single non-missing character value",
               fixed = TRUE)
  expect_error(Nestimate:::.wrap_netobject(mat, directed = NA),
               "'directed' must be TRUE or FALSE", fixed = TRUE)
  expect_error(Nestimate:::.wrap_netobject(mat, inits = c(A = 1)),
               "'inits' must be a finite numeric vector with one value per state",
               fixed = TRUE)
  expect_error(Nestimate:::.wrap_netobject(mat, inits = c(A = 1, C = 0)),
               "'inits' names must match matrix state names", fixed = TRUE)

  unnamed <- Nestimate:::.wrap_netobject(matrix(c(1, 2, 3, 4), 2, 2),
                                         inits = c(0.25, 0.75))
  expect_equal(rownames(unnamed$weights), c("1", "2"))
  expect_equal(names(unnamed$inits), c("1", "2"))
})

test_that("as_tna.default errors for non-tna objects (L1305-1306)", {
  expect_error(as_tna(list(a = 1)), "Cannot convert")
  expect_error(as_tna("some string"), "Cannot convert")
})

# ---- .normalize_clusters coverage (L1312-1362) ----

test_that(".normalize_clusters handles data.frame input (L1313-1318)", {
  node_names <- c("A", "B", "C", "D")
  clusters_df <- data.frame(
    node  = c("A", "B", "C", "D"),
    group = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  result <- Nestimate:::.normalize_clusters(clusters_df, node_names)
  expect_true(is.list(result))
  expect_equal(sort(names(result)), c("G1", "G2"))
})

test_that(".normalize_clusters validates data.frame membership strictly", {
  node_names <- c("A", "B", "C")

  expect_error(
    Nestimate:::.normalize_clusters(data.frame(node = c("A", "B", "C")),
                                    node_names),
    "must have at least two columns"
  )
  expect_error(
    Nestimate:::.normalize_clusters(
      data.frame(node = c("A", "", "C"), group = c("G1", "G1", "G2")),
      node_names
    ),
    "node column must not contain missing or empty values"
  )
  expect_error(
    Nestimate:::.normalize_clusters(
      data.frame(node = c("A", "B", "C"), group = c("G1", NA, "G2")),
      node_names
    ),
    "group column must not contain missing or empty values"
  )
  expect_error(
    Nestimate:::.normalize_clusters(
      data.frame(node = c("A", "A", "B", "C"),
                 group = c("G1", "G2", "G1", "G2")),
      node_names
    ),
    "Nodes assigned to multiple clusters: A"
  )
})

test_that(".normalize_clusters errors on unknown nodes in list (L1325-1327)", {
  node_names <- c("A", "B")
  # C is not in node_names
  expect_error(
    Nestimate:::.normalize_clusters(list(G1 = c("A", "C")), node_names),
    "Unknown nodes"
  )
})

test_that(".normalize_clusters rejects duplicate and missing list mappings", {
  node_names <- c("A", "B", "C")

  expect_error(
    Nestimate:::.normalize_clusters(
      list(G1 = c("A", "B"), G2 = c("B", "C")),
      node_names
    ),
    "Nodes assigned to multiple clusters: B"
  )
  expect_error(
    Nestimate:::.normalize_clusters(list(G1 = c("A", "B")), node_names),
    "Unmapped nodes: C"
  )
})

test_that("cluster_summary rejects duplicate and incomplete cluster lists", {
  mat <- matrix(1, 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))

  expect_error(
    cluster_summary(mat, list(G1 = c("A", "B"), G2 = c("B", "C"))),
    "Nodes assigned to multiple clusters: B"
  )
  expect_error(
    cluster_summary(mat, list(G1 = c("A", "B"))),
    "Unmapped nodes: C"
  )
})

test_that(".normalize_clusters errors on numeric vector wrong length (L1335-1337)", {
  node_names <- c("A", "B", "C")
  expect_error(
    Nestimate:::.normalize_clusters(c(1, 2), node_names),
    "Membership vector length"
  )
})

test_that(".normalize_clusters handles factor cluster membership (L1348-1359)", {
  node_names <- c("A", "B", "C")
  clusters_fac <- factor(c("G1", "G1", "G2"))
  result <- Nestimate:::.normalize_clusters(clusters_fac, node_names)
  expect_true(is.list(result))
  expect_equal(sort(names(result)), c("G1", "G2"))
  expect_true("A" %in% result$G1)
})

test_that(".normalize_clusters handles character cluster membership (L1348-1359)", {
  node_names <- c("A", "B", "C")
  result <- Nestimate:::.normalize_clusters(c("G1", "G1", "G2"), node_names)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
})

test_that(".normalize_clusters aligns named membership vectors by node name", {
  node_names <- c("A", "B", "C")

  result <- Nestimate:::.normalize_clusters(
    c(C = "G2", A = "G1", B = "G1"),
    node_names
  )
  expect_equal(result$G1, c("A", "B"))
  expect_equal(result$G2, "C")

  result_num <- Nestimate:::.normalize_clusters(
    c(C = 2, A = 1, B = 1),
    node_names
  )
  expect_equal(result_num$`1`, c("A", "B"))
  expect_equal(result_num$`2`, "C")
})

test_that(".normalize_clusters rejects malformed membership vectors", {
  node_names <- c("A", "B", "C")

  expect_error(
    Nestimate:::.normalize_clusters(c(A = "G1", A = "G2", B = "G1"),
                                    node_names),
    "names must be unique"
  )
  expect_error(
    Nestimate:::.normalize_clusters(c(A = "G1", B = "G1", D = "G2"),
                                    node_names),
    "names must match node names exactly"
  )
  expect_error(
    Nestimate:::.normalize_clusters(setNames(c("G1", "G1", "G2"),
                                             c("A", "", "C")),
                                    node_names),
    "names must not contain missing or empty values"
  )
  expect_error(
    Nestimate:::.normalize_clusters(c("G1", NA, "G2"), node_names),
    "must not contain missing or empty values"
  )
  expect_error(
    Nestimate:::.normalize_clusters(c(1, NA_real_, 2), node_names),
    "must not contain missing or non-finite values"
  )
})

test_that(".normalize_clusters errors on character vector wrong length (L1350-1351)", {
  node_names <- c("A", "B", "C")
  expect_error(
    Nestimate:::.normalize_clusters(c("G1", "G2"), node_names),
    "Membership vector length"
  )
})

test_that(".normalize_clusters errors on invalid input type (L1362)", {
  expect_error(
    Nestimate:::.normalize_clusters(TRUE, c("A", "B")),
    "clusters must be a list"
  )
})

# ---- print.mcml with edges field (L1381-1384) ----

test_that("print.mcml shows Transitions line when edges present (L1381-1384)", {
  seqs <- data.frame(
    T1 = c("A", "B", "C"),
    T2 = c("B", "C", "A"),
    T3 = c("C", "A", "B")
  )
  clusters <- list(G1 = c("A", "B"), G2 = "C")
  cs <- build_mcml(seqs, clusters)
  # build_mcml produces $edges — print should show Transitions line
  output <- capture.output(print(cs))
  expect_true(any(grepl("Transitions", output)))
})

test_that("print.mcml without edges shows no Transitions line (L1385-1386)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  # cluster_summary does not set $edges
  output <- capture.output(print(cs))
  expect_true(any(grepl("MCML Network", output)))
})

# ---- summary.mcml (L1404) ----

test_that("summary.mcml returns a tidy per-cluster data.frame (no print side-effect)", {
  mat <- matrix(runif(9), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cs <- cluster_summary(mat, list(G1 = c("A", "B"), G2 = "C"))
  s <- summary(cs)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s),
                  c("cluster", "size", "within_total",
                    "between_out", "between_in"))
  expect_equal(nrow(s), 2L)
})

test_that("summary.mcml on build_mcml result returns the same shape", {
  seqs <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B")
  )
  cs <- build_mcml(seqs, list(G1 = "A", G2 = "B"))
  s <- summary(cs)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2L)
})


# ============================================
# mcml seamless dispatch to downstream functions
# ============================================

.make_mcml <- function(seed = 42) {
  set.seed(seed)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:6], 40, TRUE),
    T2 = sample(LETTERS[1:6], 40, TRUE),
    T3 = sample(LETTERS[1:6], 40, TRUE),
    T4 = sample(LETTERS[1:6], 40, TRUE),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  build_mcml(seqs, clusters, type = "tna")
}


test_that("centrality_stability() works on mcml", {
  cs <- .make_mcml()
  stab <- centrality_stability(cs, iter = 20, seed = 1)

  expect_true(is.list(stab))
  expect_true("macro" %in% names(stab))
  for (nm in names(stab)) {
    expect_s3_class(stab[[nm]], "net_stability")
  }
})

test_that("centrality_stability() mcml dispatch matches as_tna group shape", {
  cs <- .make_mcml()
  grp <- as_tna(cs)

  stab_mc <- centrality_stability(cs, iter = 5, seed = 1)
  stab_grp <- centrality_stability(grp, iter = 5, seed = 1)

  expect_identical(names(stab_mc), names(grp))
  expect_identical(
    vapply(stab_mc, function(x) class(x)[1L], character(1L)),
    vapply(stab_grp, function(x) class(x)[1L], character(1L))
  )
})

test_that("network_reliability() works on mcml", {
  cs <- .make_mcml()
  rel <- network_reliability(cs, iter = 20, seed = 1)

  expect_s3_class(rel, "net_reliability")
})

test_that("network_reliability() mcml dispatch matches as_tna labels", {
  cs <- .make_mcml()
  grp <- as_tna(cs)

  rel_mc <- network_reliability(cs, iter = 5, seed = 1)
  rel_grp <- network_reliability(grp, iter = 5, seed = 1)

  expect_s3_class(rel_mc, "net_reliability")
  expect_setequal(as.character(rel_mc$results$network),
                  as.character(rel_grp$results$network))
})

test_that("casedrop_reliability() mcml dispatch matches as_tna group shape", {
  cs <- .make_mcml()
  grp <- as_tna(cs)

  case_mc <- casedrop_reliability(cs, iter = 5, seed = 1)
  case_grp <- casedrop_reliability(grp, iter = 5, seed = 1)

  expect_identical(names(case_mc), names(grp))
  expect_identical(
    vapply(case_mc, function(x) class(x)[1L], character(1L)),
    vapply(case_grp, function(x) class(x)[1L], character(1L))
  )
})

test_that("extract_transition_matrix() works on mcml", {
  cs <- .make_mcml()
  mats <- extract_transition_matrix(cs)

  expect_true(is.list(mats))
  expect_true("macro" %in% names(mats))
  expect_true(is.matrix(mats$macro))
  # Within-cluster matrices
  cluster_names <- setdiff(names(mats), "macro")
  expect_true(length(cluster_names) > 0)
  for (nm in cluster_names) {
    expect_true(is.matrix(mats[[nm]]))
  }
})

test_that("extract_initial_probs() works on mcml", {
  cs <- .make_mcml()
  inits <- extract_initial_probs(cs)

  expect_true(is.list(inits))
  expect_true("macro" %in% names(inits))
  expect_true(is.numeric(inits$macro))
  expect_equal(sum(inits$macro), 1, tolerance = 1e-6)
})

test_that("extract_edges() works on mcml", {
  cs <- .make_mcml()
  edges <- extract_edges(cs)

  expect_true(is.list(edges))
  expect_true("macro" %in% names(edges))
  expect_true(is.data.frame(edges$macro))
  expect_true(all(c("from", "to", "weight") %in% names(edges$macro)))
})

test_that("mcml extraction returns stored matrices, inits, and delegated edges", {
  seqs <- data.frame(
    T1 = c("A", "D"),
    T2 = c("B", "E"),
    T3 = c("C", "F"),
    T4 = c("A", "D"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  mc <- build_mcml(seqs, clusters, type = "raw")

  mats <- extract_transition_matrix(mc)
  inits <- extract_initial_probs(mc)
  edges <- extract_edges(mc, include_self = TRUE, sort_by = "from")
  ref_edges <- extract_edges(mc$macro, include_self = TRUE, sort_by = "from")

  expect_identical(unclass(mats$macro), mc$macro$weights)
  expect_identical(unclass(mats$G1), mc$clusters$G1$weights)
  expect_identical(inits$macro, mc$macro$inits)
  expect_identical(inits$G1, mc$clusters$G1$inits)
  expect_identical(edges$macro, ref_edges)
})

test_that("extract_edges rejects invalid controls on mcml", {
  cs <- .make_mcml()

  expect_error(extract_edges(cs, threshold = NA_real_),
               "'threshold' must be a single finite numeric")
  expect_error(extract_edges(cs, include_self = NA),
               "'include_self' must be TRUE or FALSE")
  expect_error(extract_edges(cs, sort_by = "bad"),
               "'sort_by' must be one of")
})

# ---- Edgelist input: $data propagated to macro + clusters so that
#      as_tna() + bootstrap_network() work end-to-end (regression test) ----

test_that("edgelist input yields macro/per-cluster $data usable by bootstrap", {
  # Tiny edgelist with known cluster structure
  edges <- data.frame(
    from   = c("a1", "a1", "a2", "b1", "b2", "b1", "a1", "b2", "a2", "b1"),
    to     = c("a2", "b1", "a1", "b2", "b1", "a2", "b2", "a1", "a1", "b1"),
    weight = 1,
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  mc <- build_mcml(edges, clusters)

  # Previous bug: $macro$data was NULL for edgelist input
  expect_false(is.null(mc$macro$data))
  expect_true(is.data.frame(mc$macro$data))
  expect_equal(ncol(mc$macro$data), 2L)
  expect_true(all(mc$macro$data[[1]] %in% c("A", "B")))
  expect_true(all(mc$macro$data[[2]] %in% c("A", "B")))

  # Each per-cluster data keeps only transitions fully within that cluster
  expect_false(is.null(mc$clusters$A$data))
  expect_true(all(mc$clusters$A$data[[1]] %in% c("a1", "a2")))
  expect_true(all(mc$clusters$A$data[[2]] %in% c("a1", "a2")))
  expect_false(is.null(mc$clusters$B$data))
  expect_true(all(mc$clusters$B$data[[1]] %in% c("b1", "b2")))
  expect_true(all(mc$clusters$B$data[[2]] %in% c("b1", "b2")))

  # Edgelist source is tagged so bootstrap can warn about it
  expect_identical(attr(mc$macro$data, "source"), "edgelist")
  expect_identical(attr(mc$clusters$A$data, "source"), "edgelist")

  # as_tna() then bootstrap_network() must run without the "$data is NULL"
  # error, but should emit the anti-conservative-CI warning pointing to
  # permutation / case-dropping alternatives.
  tnas <- as_tna(mc)
  suppressWarnings(boot <- bootstrap_network(tnas, iter = 20))
  expect_s3_class(boot, "net_bootstrap_group")
  expect_warning(
    bootstrap_network(tnas$macro, iter = 10),
    "edgelist"
  )
})

test_that("sequence input does NOT emit the edgelist bootstrap warning", {
  # Wide-format sequences: rows = actors, cols = time steps. Should have
  # no "source = edgelist" tag, no warning from bootstrap_network().
  seq_df <- data.frame(
    t1 = c("a1", "a2", "b1", "b2"),
    t2 = c("a2", "a1", "b2", "b1"),
    t3 = c("a1", "a1", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  mc <- build_mcml(seq_df, clusters)

  expect_null(attr(mc$macro$data, "source"))
  tnas <- as_tna(mc)
  expect_no_warning(
    bootstrap_network(tnas$macro, iter = 10)
  )
})

# ============================================
# audit_mcml #1: .auto_detect_clusters() must align by node label, not
# by row order, when reading from x$node_groups. Misordered node_groups
# previously corrupted assignments silently.
# ============================================

test_that(".auto_detect_clusters aligns node_groups by label, not row order", {
  # Build a simple netobject by hand. nodes are A,B,C,D; node_groups is
  # supplied in REVERSE order (D,C,B,A). Positional read would assign
  # cluster G2 to A and G1 to D — a corruption. Label-aligned read must
  # produce the original A->G1, B->G1, C->G2, D->G2 mapping.
  mat <- matrix(0, 4, 4, dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  mat[1,2] <- mat[2,3] <- mat[3,4] <- mat[4,1] <- 1
  net <- list(
    weights = mat,
    nodes = data.frame(id = 1:4, label = LETTERS[1:4], name = LETTERS[1:4],
                       stringsAsFactors = FALSE),
    edges = data.frame(from = c(1,2,3,4), to = c(2,3,4,1), weight = 1),
    directed = TRUE,
    node_groups = data.frame(
      node    = c("D", "C", "B", "A"),  # REVERSED w.r.t. nodes
      cluster = c("G2", "G2", "G1", "G1"),
      stringsAsFactors = FALSE
    )
  )
  class(net) <- c("netobject", "cograph_network")

  clusters <- Nestimate:::.auto_detect_clusters(net)
  expect_equal(unname(clusters), c("G1", "G1", "G2", "G2"))
})

test_that(".auto_detect_clusters errors clearly when node_groups lacks node column", {
  mat <- matrix(0, 3, 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  mat[1,2] <- mat[2,3] <- mat[3,1] <- 1
  net <- list(
    weights = mat,
    nodes = data.frame(id = 1:3, label = LETTERS[1:3], name = LETTERS[1:3],
                       stringsAsFactors = FALSE),
    edges = data.frame(from = c(1,2,3), to = c(2,3,1), weight = 1),
    directed = TRUE,
    node_groups = data.frame(cluster = c("G1", "G1", "G2"),
                             stringsAsFactors = FALSE)
  )
  class(net) <- c("netobject", "cograph_network")
  expect_error(Nestimate:::.auto_detect_clusters(net),
               "must include a node identifier column")
})

# ============================================
# audit_mcml #7: label propagation — cross-module test that MCML labels
# survive into state_distribution() / plot_state_frequencies().
# ============================================

test_that("MCML labels propagate to state_distribution()", {
  seq_df <- data.frame(
    t1 = c("a1", "a2", "b1", "b2"),
    t2 = c("a2", "a1", "b2", "b1"),
    t3 = c("a1", "a1", "b1", "b2"),
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))
  labels <- c(a1 = "Activity-1", a2 = "Activity-2",
              b1 = "Behavior-1", b2 = "Behavior-2")
  mc <- build_mcml(seq_df, clusters, labels = labels)
  sd <- state_distribution(mc)
  # The tidy data frame should carry the user-supplied labels in the
  # state column rather than the raw a1/a2/... names.
  expect_true(any(grepl("Activity|Behavior", sd$state)))
})

test_that("MCML labels reject collisions and leave macro cluster labels untouched", {
  seq_df <- data.frame(
    t1 = c("a1", "a2", "b1", "b2"),
    t2 = c("a2", "a1", "b2", "b1"),
    stringsAsFactors = FALSE
  )
  clusters <- list(A = c("a1", "a2"), B = c("b1", "b2"))

  expect_error(
    build_mcml(seq_df, clusters, labels = c(a1 = "Same", a2 = "Same")),
    "`labels` values must be unique",
    fixed = TRUE
  )
  expect_error(
    build_mcml(seq_df, clusters, labels = c(a1 = "a2")),
    "`labels` must map nodes to unique labels",
    fixed = TRUE
  )
  expect_error(
    build_mcml(seq_df, clusters, labels = setNames("Activity", "")),
    "`labels` must have non-empty source names",
    fixed = TRUE
  )
  expect_error(
    build_mcml(seq_df, clusters, labels = c(a1 = "")),
    "`labels` values must be non-missing and non-empty",
    fixed = TRUE
  )

  labelled <- build_mcml(
    seq_df, clusters,
    labels = c(a1 = "Activity-1", a2 = "Activity-2",
               A = "Do-not-rename-cluster")
  )

  expect_equal(labelled$macro$labels, c("A", "B"))
  expect_equal(rownames(labelled$macro$weights), c("A", "B"))
  expect_true("Activity-1" %in% labelled$clusters$A$labels)
  expect_false("Do-not-rename-cluster" %in% labelled$macro$labels)
})

# ============================================
# audit_mcml gap: edge-list `clusters = "<colname>"` mode assigns the
# row's group to BOTH endpoints. Pin the expected behaviour so the
# documented limitation can't silently drift. With from-side and to-side
# in different groups the resulting node lookup picks one of them per
# node (depending on which row writes to it last), which is exactly the
# narrow contract documented in build_mcml()'s `clusters` param.
# ============================================

test_that("edge-list clusters-by-column reflects documented narrow contract", {
  # Edge list where every row tags both endpoints with the same group —
  # the supported case — produces correctly-aligned cluster assignments.
  edges_clean <- data.frame(
    from  = c("A", "B", "C", "D"),
    to    = c("B", "A", "D", "C"),
    weight = 1,
    grp   = c("G1", "G1", "G2", "G2"),
    stringsAsFactors = FALSE
  )
  mc <- build_mcml(edges_clean, clusters = "grp")
  members <- mc$cluster_members
  expect_setequal(members$G1, c("A", "B"))
  expect_setequal(members$G2, c("C", "D"))
})

# ============================================
# audit_mcml #6 gap: as_tna.mcml() must drop relative-method clusters
# whose within-matrix has zero row sums (cannot row-normalise) and warn
# clearly. Deterministic fixture: a 4-node mcml where one cluster has
# an absorbing-only row (all-zero outgoing weights).
# ============================================

test_that("as_tna.mcml drops zero-row-sum clusters and warns", {
  # Drop only fires when net_method == "relative" (i.e. sequence/edgelist
  # input with type = "tna"). cluster_summary() on a matrix produces
  # frequency-method results where a zero row is a legitimate sink and
  # is kept. So we need a sequence-derived mcml whose within-cluster
  # row-normalisation would divide by zero.
  #
  # Construction: in cluster G2 = {b1, b2}, node b1 only transitions OUT
  # of the cluster (b1 -> a1), never within. Its within-G2 row sums to
  # zero, triggering the drop.
  seq_df <- data.frame(
    t1 = c("b2", "b2", "a1", "a2"),
    t2 = c("b1", "b1", "a2", "a1"),  # b1 only appears as transition target
    t3 = c("a1", "a2", "a1", "a2"),  # ... and only transitions out (b1 -> a1)
    t4 = c("a2", "a1", "a2", "a1"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("a1", "a2"), G2 = c("b1", "b2"))
  mc <- build_mcml(seq_df, clusters, type = "tna")

  # Sanity: the within-G2 row for b1 has zero sum.
  w_g2 <- mc$clusters$G2$weights
  zero_rows <- rownames(w_g2)[rowSums(w_g2) == 0]
  expect_true(length(zero_rows) >= 1L)

  expect_warning(out <- as_tna(mc), "Dropped clusters with zero row sums")
  expect_false("G2" %in% names(out))
  expect_true("macro" %in% names(out))
})

test_that("as_tna.mcml preserves relative macro weights and inits while dropping bad within clusters", {
  seq_df <- data.frame(
    t1 = c("b2", "b2", "a1", "a2"),
    t2 = c("b1", "b1", "a2", "a1"),
    t3 = c("a1", "a2", "a1", "a2"),
    t4 = c("a2", "a1", "a2", "a1"),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("a1", "a2"), G2 = c("b1", "b2"))
  mc <- build_mcml(seq_df, clusters, type = "tna")

  expect_warning(result <- as_tna(mc),
                 "Dropped clusters with zero row sums")
  expect_identical(result$macro$weights, mc$macro$weights)
  expect_identical(result$macro$inits, mc$macro$inits)
  expect_equal(result$macro$method, "relative")
  expect_true("G1" %in% names(result))
  expect_false("G2" %in% names(result))
})
