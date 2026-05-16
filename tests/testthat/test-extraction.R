# Tests for R/extraction.R

# -- Shared fixtures --------------------------------------------------------
make_mat <- function(named = TRUE) {
  m <- matrix(c(0, 0.3, 0.7,
                0.4, 0, 0.6,
                0.5, 0.5, 0), nrow = 3, byrow = TRUE)
  if (named) rownames(m) <- colnames(m) <- c("A", "B", "C")
  m
}

make_tna_model <- function(field = "weights", cls = "tna") {
  obj <- list()
  obj[[field]] <- make_mat()
  class(obj) <- cls
  obj
}

# == extract_transition_matrix ==============================================

test_that("extracts from tna-class $weights", {
  result <- extract_transition_matrix(make_tna_model("weights", "tna"))
  expect_equal(result, make_mat())
})

test_that("extracts from ftna/ctna/atna $transition", {
  vapply(c("ftna", "ctna", "atna"), function(cls) {
    result <- extract_transition_matrix(make_tna_model("transition", cls))
    expect_equal(result, make_mat())
    TRUE
  }, logical(1))
})

test_that("extracts from generic list with different field names", {
  fields <- c("weights", "transition_matrix", "transition")
  vapply(fields, function(f) {
    obj <- list()
    obj[[f]] <- make_mat()
    expect_equal(extract_transition_matrix(obj), make_mat())
    TRUE
  }, logical(1))
})

test_that("accepts direct matrix input", {
  expect_equal(extract_transition_matrix(make_mat()), make_mat())
})

test_that("errors when no matrix found", {
  expect_error(extract_transition_matrix(list(x = 1)), "Could not extract")
  expect_error(extract_transition_matrix("bad"), "Could not extract")
})

test_that("type='scaled' row-normalizes", {
  m <- matrix(c(2, 4, 6,
                0, 0, 0,
                1, 2, 3), nrow = 3, byrow = TRUE)
  rownames(m) <- colnames(m) <- c("X", "Y", "Z")
  result <- extract_transition_matrix(m, type = "scaled")
  # Row 1 sums to 1
  expect_equal(sum(result[1, ]), 1)
  expect_equal(result[1, 1], 2 / 12)
  # Row 2 (zero row) stays zero, no NaN
  expect_equal(sum(result[2, ]), 0)
  expect_true(!any(is.nan(result)))
})

test_that("type='raw' returns unmodified values", {
  m <- matrix(c(10, 20, 30, 40), nrow = 2)
  expect_equal(extract_transition_matrix(m, type = "raw"), m)
})

# == extract_initial_probs ==================================================

test_that("extracts $initial from tna-class", {
  obj <- make_tna_model()
  obj$initial <- c(A = 0.5, B = 0.3, C = 0.2)
  expect_equal(extract_initial_probs(obj), obj$initial)
})

test_that("extracts $initial_probs from tna-class", {
  obj <- make_tna_model()
  obj$initial_probs <- c(A = 0.2, B = 0.8)
  expect_equal(extract_initial_probs(obj), c(A = 0.2, B = 0.8))
})

test_that("extracts from list with $initial_probabilities", {
  obj <- list(initial_probabilities = c(X = 0.6, Y = 0.4))
  expect_equal(extract_initial_probs(obj), obj$initial_probabilities)
})

test_that("falls back to uniform with warning when weights exist", {
  obj <- list(weights = make_mat())
  expect_warning(result <- extract_initial_probs(obj), "uniform")
  expect_equal(length(result), 3)
  expect_equal(sum(result), 1)
  expect_equal(as.numeric(result), rep(1 / 3, 3))
})

test_that("normalizes probabilities that do not sum to 1", {
  obj <- list(initial = c(A = 2, B = 3))
  result <- extract_initial_probs(obj)
  expect_equal(sum(result), 1)
  expect_equal(result[["A"]], 0.4)
})

test_that("adds S-names to unnamed initial vector", {
  obj <- list(initial = c(0.5, 0.5))
  result <- extract_initial_probs(obj)
  expect_equal(names(result), c("S1", "S2"))
})

test_that("errors when nothing extractable", {
  expect_error(extract_initial_probs(list(x = 1)), "Could not extract")
})

# == extract_edges ==========================================================

test_that("returns correct edge list from matrix", {
  m <- make_mat()
  edges <- extract_edges(m)
  expect_s3_class(edges, "data.frame")
  expect_named(edges, c("from", "to", "weight"))
  # Default: no self-loops, threshold 0 filters exact zeros
  expect_true(all(edges$from != edges$to))
  expect_true(all(edges$weight >= 0))
})

test_that("include_self adds self-loops", {
  m <- diag(3)
  rownames(m) <- colnames(m) <- c("A", "B", "C")
  edges <- extract_edges(m, include_self = TRUE)
  self_edges <- edges[edges$from == edges$to, ]
  expect_equal(nrow(self_edges), 3)
})

test_that("threshold filters edges", {
  m <- make_mat()
  edges <- extract_edges(m, threshold = 0.5)
  expect_true(all(edges$weight >= 0.5))
})

test_that("sort_by='weight' sorts descending", {
  edges <- extract_edges(make_mat())
  expect_true(all(diff(edges$weight) <= 0))
})

test_that("sort_by='from' sorts alphabetically", {
  edges <- extract_edges(make_mat(), sort_by = "from")
  expect_true(!is.unsorted(edges$from))
})

test_that("sort_by=NULL returns unsorted (no error)", {
  edges <- extract_edges(make_mat(), sort_by = NULL)
  expect_s3_class(edges, "data.frame")
})

test_that("unnamed matrix gets S-prefixed names", {
  m <- matrix(c(0, 1, 0, 0), nrow = 2)
  edges <- extract_edges(m, include_self = TRUE)
  expect_true(all(edges$from %in% c("S1", "S2")))
})

test_that("single-state matrix with self-loop", {
  m <- matrix(1, nrow = 1, dimnames = list("X", "X"))
  edges <- extract_edges(m, include_self = TRUE)
  expect_equal(nrow(edges), 1)
  expect_equal(edges$from, "X")
  expect_equal(edges$to, "X")
  # Without self-loop: empty

  edges_no <- extract_edges(m, include_self = FALSE)
  expect_equal(nrow(edges_no), 0)
})

test_that("extract_edges rejects invalid controls", {
  m <- make_mat()

  expect_error(extract_edges(m, threshold = NA_real_),
               "'threshold' must be a single finite numeric")
  expect_error(extract_edges(m, threshold = c(0, 1)),
               "'threshold' must be a single finite numeric")
  expect_error(extract_edges(m, include_self = NA),
               "'include_self' must be TRUE or FALSE")
  expect_error(extract_edges(m, include_self = c(TRUE, FALSE)),
               "'include_self' must be TRUE or FALSE")
  expect_error(extract_edges(m, sort_by = "bad"),
               "'sort_by' must be one of")
  expect_error(extract_edges(m, sort_by = NA_character_),
               "'sort_by' must be one of")
})


# == extract_initial_probs: additional paths ================================

test_that("extracts $initial_probs from tna-class (L132-133)", {
  # tna/ftna/ctna/atna class with $initial_probs (not $initial)
  obj <- make_tna_model()
  obj$initial_probs <- c(A = 0.4, B = 0.3, C = 0.3)
  class(obj) <- "tna"
  result <- extract_initial_probs(obj)
  expect_equal(result, c(A = 0.4, B = 0.3, C = 0.3))
})

test_that("extracts $initial_probs from tna-subclass ftna (L132-133)", {
  obj <- make_tna_model("transition", "ftna")
  obj$initial_probs <- c(A = 0.5, B = 0.5)
  result <- extract_initial_probs(obj)
  expect_equal(result, c(A = 0.5, B = 0.5))
})

test_that("extracts $initial_probs from generic list (L138-139)", {
  # Non-tna list with $initial_probs (not $initial)
  obj <- list(initial_probs = c(X = 0.6, Y = 0.4))
  result <- extract_initial_probs(obj)
  expect_equal(result, c(X = 0.6, Y = 0.4))
})

test_that("uniform fallback uses S-names when matrix has no rownames (L154)", {
  # Transition matrix with no dimnames -> states inferred as S1, S2
  m <- matrix(c(0, 1, 1, 0), nrow = 2)
  obj <- list(weights = m)
  expect_warning(result <- extract_initial_probs(obj), "uniform")
  expect_equal(names(result), c("S1", "S2"))
  expect_equal(as.numeric(result), c(0.5, 0.5))
})

# == extract_edges: sort_by = "to" ==========================================

test_that("sort_by='to' sorts by to column then from", {
  edges <- extract_edges(make_mat(), sort_by = "to")
  # to column should be non-decreasing
  expect_true(!is.unsorted(edges$to))
})
