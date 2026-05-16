# ---- prepare_onehot() tests ----

test_that("prepare_onehot basic conversion", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 0, 0, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B", "C"))
  expect_true(is.data.frame(result))
  # Should have columns like W0_T1, W0_T2, etc.
  expect_true(any(grepl("^W\\d+_T\\d+$", names(result))))
  # A and B should appear; C should be NA (all zeros)
  vals <- unlist(result[, grepl("^W\\d+_T\\d+$", names(result))])
  expect_true("A" %in% vals)
  expect_true("B" %in% vals)
})

test_that("prepare_onehot with actor grouping", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2),
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"), actor = "actor")
  expect_true(is.data.frame(result))
  # Should have at least 2 rows (one per actor)
  expect_true(nrow(result) >= 2)
})

test_that("prepare_onehot non-overlapping window", {
  df <- data.frame(
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "non-overlapping")
  expect_true(is.data.frame(result))
  expect_true(attr(result, "windowed"))
  expect_equal(attr(result, "window_size"), 2L)
})

test_that("prepare_onehot overlapping window", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "overlapping")
  expect_true(is.data.frame(result))
  expect_true(attr(result, "windowed"))
})

test_that("prepare_onehot aggregate mode", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           window_size = 2, window_type = "non-overlapping",
                           aggregate = TRUE)
  expect_true(is.data.frame(result))
})

test_that("prepare_onehot default actor/session (single group)", {
  df <- data.frame(
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )

  result <- prepare_onehot(df, cols = c("A", "B"))
  expect_true(is.data.frame(result))
  # Should have 1 row (single group)
  expect_equal(nrow(result), 1)
})

test_that("prepare_onehot attributes set correctly", {
  df <- data.frame(
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )

  # Pinned to size 1 (was the implicit default before 2026-05-10).
  result <- prepare_onehot(df, cols = c("A", "B"), window_size = 1L)
  expect_true(attr(result, "windowed"))
  expect_equal(attr(result, "window_size"), 1L)
  expect_equal(attr(result, "window_span"), 2L)
  expect_equal(attr(result, "codes"), c("A", "B"))

  result2 <- prepare_onehot(df, cols = c("A", "B"), window_size = 2)
  expect_true(attr(result2, "windowed"))
})

test_that("prepare_onehot with session", {
  df <- data.frame(
    actor = c(1, 1, 1, 1),
    session = c("s1", "s1", "s2", "s2"),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- prepare_onehot(df, cols = c("A", "B"),
                           actor = "actor", session = "session")
  expect_true(is.data.frame(result))
})

test_that("prepare_onehot output feeds into build_network", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1),
    C = c(0, 0, 0, 0, 0, 0),
    actor = c(1, 1, 1, 2, 2, 2)
  )

  result <- prepare_onehot(df, cols = c("A", "B", "C"), actor = "actor")
  # The wide format should be usable for build_network
  seq_cols <- grep("^W\\d+_T\\d+$", names(result), value = TRUE)
  if (length(seq_cols) >= 2) {
    net <- build_network(result, method = "relative",
                         params = list(format = "wide", cols = seq_cols))
    expect_s3_class(net, "netobject")
  }
})

test_that("prepare_onehot validates inputs", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(prepare_onehot(df, cols = c("X", "Y")))
  expect_error(prepare_onehot(df, cols = "A", actor = "missing_col"))
})

test_that("prepare_onehot preserves multi-hot rows", {
  df <- data.frame(
    A = c(1L, 0L),
    B = c(1L, 1L),
    C = c(0L, 0L)
  )

  # Pinned to size 1 (was the implicit default before 2026-05-10).
  result <- prepare_onehot(df, cols = c("A", "B", "C"), window_size = 1L)
  expect_equal(sum(!is.na(as.matrix(result))), 3L)
  expect_equal(result$W0_T1, "A")
  expect_equal(result$W0_T2, "B")
  expect_equal(result$W1_T2, "B")
})

test_that("prepare_onehot matches tna import_onehot across 100 edge-case data sets", {
  testthat::skip_if_not_installed("tna")
  set.seed(20260507)
  codes <- c("A", "B", "C", "D")

  make_case <- function(i) {
    n <- sample(1:8, 1L)
    mat <- matrix(sample(c(0L, 1L, NA_integer_), n * length(codes),
                         replace = TRUE, prob = c(0.52, 0.42, 0.06)),
                  nrow = n, dimnames = list(NULL, codes))
    if (i %% 10L == 0L) mat[] <- 0L
    if (i %% 10L == 1L) mat[1L, ] <- c(1L, 1L, 0L, NA_integer_)
    df <- as.data.frame(mat)
    if (i %% 4L %in% c(1L, 3L)) {
      df$actor <- rep(letters[1:3], length.out = n)
    }
    if (i %% 4L %in% c(2L, 3L)) {
      df$session <- rep(c("s1", "s2"), length.out = n)
    }
    list(
      data = df,
      window_size = sample(1:min(4L, max(1L, n)), 1L),
      window_type = if (i %% 3L == 0L) "overlapping" else "non-overlapping",
      aggregate = i %% 5L == 0L,
      interval = if (i %% 6L == 0L) sample(1:3, 1L) else NULL,
      actor = if ("actor" %in% names(df)) "actor" else NULL,
      session = if ("session" %in% names(df)) "session" else NULL
    )
  }

  compare_case <- function(case) {
    # Known intentional divergence: tna::import_onehot drops the leading
    # row of every group even when window_size == 1L (no overlap to
    # remove). Nestimate's prepare_onehot was patched on 2026-05-10 to
    # preserve the row, so it intentionally no longer matches tna here.
    # See `tests/testthat/test-equiv-windowed-fixes.R` for the locked
    # post-fix behavior.
    if (case$window_size == 1L && case$window_type == "overlapping") {
      return(TRUE)
    }
    ours_args <- list(
      data = case$data,
      cols = codes,
      window_size = case$window_size,
      window_type = case$window_type,
      aggregate = case$aggregate
    )
    if (!is.null(case$interval)) ours_args$interval <- case$interval
    if (!is.null(case$actor)) ours_args$actor <- case$actor
    if (!is.null(case$session)) ours_args$session <- case$session
    ours <- do.call(prepare_onehot, ours_args)

    tna_args <- list(
      data = case$data,
      cols = quote(c(A, B, C, D)),
      window_size = case$window_size,
      window_type = if (case$window_type == "overlapping") "sliding" else "tumbling",
      aggregate = case$aggregate
    )
    if (!is.null(case$interval)) tna_args$interval <- case$interval
    if (!is.null(case$actor)) tna_args$actor <- quote(actor)
    if (!is.null(case$session)) tna_args$session <- quote(session)
    ref <- as.data.frame(do.call(tna::import_onehot, tna_args))

    expect_equal(names(ours), names(ref))
    expect_equal(unname(as.matrix(ours)), unname(as.matrix(ref)))
    expect_equal(attr(ours, "windowed"), attr(ref, "windowed"))
    expect_equal(attr(ours, "window_size"), attr(ref, "window_size"))
    expect_equal(attr(ours, "window_span"), attr(ref, "window_span"))

    if (nrow(ours) > 0L && ncol(ours) > 1L && any(!is.na(as.matrix(ours)))) {
      seq_cols <- names(ours)
      freq <- build_network(ours, method = "frequency",
                            params = list(format = "wide", cols = seq_cols))
      rel <- build_network(ours, method = "relative",
                           params = list(format = "wide", cols = seq_cols))
      co <- build_network(ours, method = "co_occurrence",
                          params = list(format = "wide", cols = seq_cols))
      expect_equal(freq$weights, tna::ftna(ref)$weights)
      expect_equal(rel$weights, tna::tna(ref)$weights)
      expect_equal(co$weights, tna::ctna(ref)$weights)
    }
    TRUE
  }

  results <- lapply(lapply(seq_len(100L), make_case), compare_case)
  expect_equal(length(results), 100L)
})
