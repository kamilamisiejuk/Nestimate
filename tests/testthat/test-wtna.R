# ---- wtna() tests ----

test_that("wtna transition counts match manual crossprod", {
  # 3-code, 5-row one-hot data
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 0, 0, 0)
  )

  # Pinned to size 1 (was the implicit default before 2026-05-10);
  # the manual crossprod expectation is for size-1 transitions.
  net <- wtna(df, method = "transition", type = "frequency", window_size = 1L)
  expect_s3_class(net, "netobject")
  expect_true(net$directed)

  # Manual: crossprod(X[-5,], X[-1,])
  X <- as.matrix(df)
  expected <- crossprod(X[-5, ], X[-1, ])
  expect_equal(unname(net$weights), unname(expected))
})

test_that("wtna cooccurrence counts match crossprod", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 1, 0, 0)
  )

  # Pinned to size 1 — the manual crossprod is the size-1 reference.
  net <- wtna(df, method = "cooccurrence", type = "frequency", window_size = 1L)
  expect_s3_class(net, "netobject")
  expect_false(net$directed)

  # Manual: crossprod(X)
  X <- as.matrix(df)
  expected <- crossprod(X)
  expect_equal(unname(net$weights), unname(expected))
})

test_that("wtna method='both' returns list of two networks", {
  df <- data.frame(
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1)
  )

  result <- wtna(df, method = "both")
  expect_true(is.list(result))
  expect_true("transition" %in% names(result))
  expect_true("cooccurrence" %in% names(result))
  expect_s3_class(result$transition, "netobject")
  expect_s3_class(result$cooccurrence, "netobject")
  expect_true(result$transition$directed)
  expect_false(result$cooccurrence$directed)
})

test_that("wtna non-overlapping window aggregation", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "frequency",
              window_size = 2, mode = "non-overlapping")
  expect_s3_class(net, "netobject")

  # Non-overlapping: rows 1-2 -> window 0, rows 3-4 -> window 1
  # Window 0: A=1,B=1 (both active); Window 1: A=1,B=1
  # Transitions: (1,1) -> (1,1) = crossprod should be 1 on all pairs
  X <- as.matrix(df)
  wid <- c(0L, 0L, 1L, 1L)
  agg <- rowsum(X, wid)
  agg[agg > 0] <- 1
  expected <- crossprod(agg[-nrow(agg), , drop = FALSE], agg[-1, , drop = FALSE])
  expect_equal(unname(net$weights), unname(expected))
})

test_that("wtna overlapping window aggregation", {
  df <- data.frame(
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "frequency",
              window_size = 2, mode = "overlapping")
  expect_s3_class(net, "netobject")

  # Overlapping with ws=2: windows (1,2), (2,3), (3,4)
  # Pairwise between-window transitions:
  # Window(1,2)->Window(2,3): rows {1,2} x rows {2,3} = 4 pairs
  # Window(2,3)->Window(3,4): rows {2,3} x rows {3,4} = 4 pairs
  X <- as.matrix(df)
  expected <- matrix(0, 2, 2)
  # w1->w2: (1,2)x(2,3)
  for (j in 1:2) for (k in 2:3) expected <- expected + tcrossprod(X[j,], X[k,])
  # w2->w3: (2,3)x(3,4)
  for (j in 2:3) for (k in 3:4) expected <- expected + tcrossprod(X[j,], X[k,])
  expect_equal(unname(net$weights), unname(expected))
})

test_that("wtna per-actor grouping", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2, 2),
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1)
  )

  # Pinned to size 1 — manual per-actor crossprod is the size-1 reference.
  net <- wtna(df, method = "transition", type = "frequency",
              codes = c("A", "B"), actor = "actor", window_size = 1L)
  expect_s3_class(net, "netobject")

  # Manual: compute per actor and sum
  g1 <- as.matrix(df[df$actor == 1, c("A", "B")])
  g2 <- as.matrix(df[df$actor == 2, c("A", "B")])
  t1 <- crossprod(g1[-3, ], g1[-1, ])
  t2 <- crossprod(g2[-3, ], g2[-1, ])
  expected <- t1 + t2
  expect_equal(unname(net$weights), unname(expected))
})

test_that("wtna treats one-hot NA indicators as inactive", {
  df_na <- data.frame(
    actor = c("s1", "s1", "s1", "s2", "s2"),
    A = c(1L, NA, 0L, 1L, 0L),
    B = c(0L, 1L, NA, 0L, 1L),
    C = c(0L, 0L, 1L, NA, 0L)
  )
  df_zero <- df_na
  df_zero[is.na(df_zero)] <- 0L

  trans_na <- wtna(df_na, method = "transition", type = "frequency",
                   codes = c("A", "B", "C"), actor = "actor")
  trans_zero <- wtna(df_zero, method = "transition", type = "frequency",
                     codes = c("A", "B", "C"), actor = "actor")
  co_na <- wtna(df_na, method = "cooccurrence", type = "frequency",
                codes = c("A", "B", "C"), actor = "actor")
  co_zero <- wtna(df_zero, method = "cooccurrence", type = "frequency",
                  codes = c("A", "B", "C"), actor = "actor")

  expect_false(anyNA(trans_na$weights))
  expect_false(anyNA(co_na$weights))
  expect_equal(trans_na$weights, trans_zero$weights)
  expect_equal(co_na$weights, co_zero$weights)
  expect_equal(trans_na$initial, trans_zero$initial)
})

test_that("wtna type='relative' row-normalizes", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0)
  )

  net <- wtna(df, method = "transition", type = "relative")
  rs <- rowSums(net$weights)
  # Non-zero rows should sum to 1
  non_zero <- rs > 0
  if (any(non_zero)) {
    expect_true(all(abs(rs[non_zero] - 1) < 1e-10))
  }
})

test_that("wtna auto-detects one-hot columns", {
  df <- data.frame(
    id = c(1, 2, 3, 4),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1),
    score = c(3.5, 2.1, 4.0, 1.5)
  )

  # id and score are not binary, should be excluded
  # Actually id column has values in 0/1 range... let's use different ids
  df2 <- data.frame(
    name = c("x", "y", "z", "w"),
    A = c(1, 0, 1, 0),
    B = c(0, 1, 0, 1),
    score = c(3.5, 2.1, 4.0, 1.5)
  )

  net <- wtna(df2, method = "transition")
  # Should only use A and B (binary columns)
  expect_equal(sort(net$nodes$label), c("A", "B"))
})

test_that("wtna single row edge case", {
  df <- data.frame(A = 1, B = 0)
  net <- wtna(df, method = "transition", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_true(all(net$weights == 0))
})

test_that("wtna all-zeros edge case", {
  df <- data.frame(
    A = c(0, 0, 0),
    B = c(0, 0, 0)
  )
  net <- wtna(df, method = "transition", type = "frequency")
  expect_s3_class(net, "netobject")
  expect_true(all(net$weights == 0))
})

test_that("wtna single actor with actor param", {
  df <- data.frame(
    actor = c(1, 1, 1),
    A = c(1, 0, 1),
    B = c(0, 1, 0)
  )
  net <- wtna(df, method = "transition", codes = c("A", "B"), actor = "actor")
  expect_s3_class(net, "netobject")
})

test_that("wtna cooccurrence is symmetric", {
  df <- data.frame(
    A = c(1, 0, 1, 1),
    B = c(0, 1, 1, 0),
    C = c(1, 0, 0, 1)
  )

  net <- wtna(df, method = "cooccurrence", type = "frequency")
  expect_true(isSymmetric(unname(net$weights)))
})

test_that("wtna validates inputs", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(wtna(df, codes = c("X", "Y")))
  expect_error(wtna(df, codes = c("A"), method = "transition"))
})

# ---- Registry integration tests ----

test_that("wtna works via build_network registry", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1),
    B = c(0, 1, 0, 1, 0),
    C = c(0, 0, 1, 0, 0)
  )

  net <- build_network(df, method = "wtna")
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  expect_equal(net$method, "wtna")

  net2 <- build_network(df, method = "cna")
  expect_s3_class(net2, "netobject")
  expect_false(net2$directed)
  expect_equal(net2$method, "co_occurrence")
})

test_that("wtna registry with window_size param", {
  df <- data.frame(
    A = c(1, 0, 0, 1, 1, 0),
    B = c(0, 1, 1, 0, 0, 1)
  )

  net <- build_network(df, method = "wtna",
                       params = list(window_size = 2,
                                     mode = "non-overlapping"))
  expect_s3_class(net, "netobject")
})

test_that("wtna registry with actor param", {
  df <- data.frame(
    actor = c(1, 1, 1, 2, 2, 2),
    A = c(1, 0, 1, 0, 1, 0),
    B = c(0, 1, 0, 1, 0, 1)
  )

  net <- build_network(df, method = "wtna",
                       params = list(codes = c("A", "B"), actor = "actor"))
  expect_s3_class(net, "netobject")
})

test_that("wtna registered, wcna is alias for co_occurrence", {
  estimators <- list_estimators()
  expect_true("wtna" %in% estimators$name)
  expect_true(estimators$directed[estimators$name == "wtna"])

  # wcna resolves to co_occurrence
  expect_equal(.resolve_method_alias("wcna"), "co_occurrence")
  expect_equal(.resolve_method_alias("cna"), "co_occurrence")
})

test_that("wtna bootstrap works via wtna()", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  net <- wtna(df, method = "transition")
  boot <- bootstrap_network(net, iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  expect_true(boot$original$directed)
  expect_true(nrow(boot$summary) > 0)
})

test_that("wtna bootstrap via build_network falls back to association path", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  # build_network(method="wtna") stores method as "wtna" (not co_occurrence),

  # so bootstrap dispatches via association path — no error but prefer wtna()
  boot <- bootstrap_network(build_network(df, method = "wtna"),
                            iter = 5, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
})

test_that("wcna bootstrap works via wtna()", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  net <- wtna(df, method = "cooccurrence")
  boot <- bootstrap_network(net, iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "wtna_cooccurrence")
  expect_false(boot$original$directed)
})

test_that("wcna bootstrap works via build_network one-hot data", {
  set.seed(42)
  df <- data.frame(
    A = sample(0:1, 30, replace = TRUE),
    B = sample(0:1, 30, replace = TRUE),
    C = sample(0:1, 30, replace = TRUE)
  )

  boot <- bootstrap_network(build_network(df, method = "cna"),
                            iter = 20, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "co_occurrence")
  expect_true(any(boot$mean > 0))
})

test_that("wtna permutation test works via wtna()", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  net1 <- wtna(df1, method = "transition")
  net2 <- wtna(df2, method = "transition")
  perm <- permutation(net1, net2, iter = 20, seed = 1)

  expect_s3_class(perm, "net_permutation")
  expect_true(is.matrix(perm$p_values))
  expect_true(all(perm$p_values >= 0 & perm$p_values <= 1))
})

test_that("wtna permutation via build_network falls back to association path", {
  set.seed(42)
  df1 <- data.frame(A = sample(0:1, 20, replace = TRUE), B = sample(0:1, 20, replace = TRUE))
  df2 <- data.frame(A = sample(0:1, 20, replace = TRUE), B = sample(0:1, 20, replace = TRUE))

  # build_network(method="wtna") stores method as "wtna" (not co_occurrence),
  # so permutation dispatches via association path
  perm <- permutation(build_network(df1, method = "wtna"),
                      build_network(df2, method = "wtna"),
                      iter = 5, seed = 1)
  expect_s3_class(perm, "net_permutation")
})

test_that("wcna permutation test works via wtna()", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  net1 <- wtna(df1, method = "cooccurrence")
  net2 <- wtna(df2, method = "cooccurrence")
  perm <- permutation(net1, net2, iter = 20, seed = 1)

  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "wtna_cooccurrence")
})

test_that("wcna permutation works via build_network one-hot data", {
  set.seed(42)
  df1 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )
  df2 <- data.frame(
    A = sample(0:1, 20, replace = TRUE),
    B = sample(0:1, 20, replace = TRUE)
  )

  perm <- permutation(build_network(df1, method = "cna"),
                      build_network(df2, method = "cna"),
                      iter = 20, seed = 1)
  expect_s3_class(perm, "net_permutation")
  expect_equal(perm$method, "co_occurrence")
})


# ---- .resolve_codes() coverage ----

test_that(".resolve_codes with numeric indices returns column names (L219-220)", {
  df <- data.frame(id = 1:3, A = c(1, 0, 1), B = c(0, 1, 0))
  result <- Nestimate:::.resolve_codes(df, codes = 2:3)
  expect_equal(result, c("A", "B"))
})

test_that(".resolve_codes with column range string returns range (L225-231)", {
  df <- data.frame(A = c(1, 0), B = c(0, 1), C = c(1, 1), D = c(0, 0))
  result <- Nestimate:::.resolve_codes(df, codes = "B:D")
  expect_equal(result, c("B", "C", "D"))
})

test_that(".resolve_codes range errors on missing start column (L229)", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(
    Nestimate:::.resolve_codes(df, codes = "Z:B"),
    "not found"
  )
})

test_that(".resolve_codes range errors on missing end column (L230)", {
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  expect_error(
    Nestimate:::.resolve_codes(df, codes = "A:Z"),
    "not found"
  )
})

# ---- wtna: overlapping window with fewer rows than window_size ----

test_that("wtna returns zero matrix when n < window_size (overlapping)", {
  # 2-row data, window_size = 5 -> too few rows -> zero weights
  df <- data.frame(A = c(1, 0), B = c(0, 1))
  net <- wtna(df, method = "transition", type = "frequency",
              window_size = 5L, mode = "overlapping")
  expect_true(all(net$weights == 0))
})

# ---- .wtna_compute_by_actor: multiple actor columns (L175-176) ----

test_that(".wtna_compute_by_actor with multi-column actor key (L175-176)", {
  df <- data.frame(
    actor   = c(1, 1, 2, 2),
    group   = c("g1", "g1", "g1", "g1"),
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )
  codes <- c("A", "B")
  result <- Nestimate:::.wtna_compute_by_actor(df, codes = codes,
                                               window_size = 1L,
                                               mode = "non-overlapping",
                                               actor = c("actor", "group"),
                                               method = "transition")
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(2L, 2L))
})

# ---- .wtna_compute_by_actor: method='both' summing (L185-188) ----

test_that(".wtna_compute_by_actor method='both' with multiple actors (L185-188)", {
  df <- data.frame(
    actor = c(1, 1, 2, 2),
    A = c(1, 0, 0, 1),
    B = c(0, 1, 1, 0)
  )
  result <- Nestimate:::.wtna_compute_by_actor(df, codes = c("A", "B"),
                                               window_size = 1L,
                                               mode = "non-overlapping",
                                               actor = "actor",
                                               method = "both")
  expect_true(is.list(result))
  expect_true(all(c("transition", "cooccurrence") %in% names(result)))
  expect_true(is.matrix(result$transition))
  expect_true(is.matrix(result$cooccurrence))
})

# ---- .estimator_wtna: sequence fallback path (L394-397) ----

test_that(".estimator_wtna falls back to frequency for long-format non-binary data (L394-397)", {
  # Long format with action column -> not one-hot -> uses .estimator_frequency
  df <- data.frame(
    id = c(1, 1, 1, 2, 2),
    time = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "A", "B", "A"),
    stringsAsFactors = FALSE
  )
  net <- build_network(df, method = "wtna", format = "long",
                       action = "Action", id = "id", time = "time")
  expect_s3_class(net, "netobject")
  expect_true(net$weights["A", "B"] > 0 || net$weights["B", "A"] > 0)
})


# ---- .wtna_initial_probs edge cases ----

test_that("wtna with all-zero data produces NULL initial (L343)", {
  df_empty <- data.frame(A = c(0,0,0,0), B = c(0,0,0,0))
  net <- wtna(df_empty, method = "transition")
  expect_null(net$initial)
})

test_that("wtna with all-zero actor data produces NULL initial (L366)", {
  df4 <- data.frame(actor = c("a","a","b","b"), A = c(0,0,0,0), B = c(0,0,0,0))
  net <- wtna(df4, method = "transition", actor = "actor")
  expect_null(net$initial)
})

test_that("wtna with mixed actor data produces valid initial (L361)", {
  df3 <- data.frame(actor = c("a","a","b","b"), A = c(0,0,1,0), B = c(0,0,0,1))
  net <- wtna(df3, method = "transition", actor = "actor")
  expect_false(is.null(net$initial))
  expect_equal(sum(net$initial), 1)
})


# ---- .estimator_wtna_core type="relative" (L460-462) ----

test_that(".estimator_wtna_core normalizes to relative (L460-462)", {
  df <- data.frame(
    A = c(1, 0, 1, 0, 1), B = c(0, 1, 0, 1, 0), C = c(0, 0, 0, 0, 0)
  )
  net <- build_network(df, method = "wtna", params = list(type = "relative"))
  expect_s3_class(net, "netobject")
  expect_true(all(rowSums(net$weights) <= 1 + 1e-10))
})


# ---- print.wtna_mixed (L544-551) ----

test_that("print.wtna_mixed shows both components (L544-551)", {
  oh <- data.frame(A = c(1,0,1,0), B = c(0,1,0,1), C = c(1,1,0,0))
  mixed <- wtna(oh, method = "both")
  out <- capture.output(print(mixed))
  expect_true(any(grepl("Mixed Window TNA", out)))
  expect_true(any(grepl("Transition", out)))
  expect_true(any(grepl("Co-occurrence", out)))
})


# ---- Vectorization equivalence tests ----
# Prove colSums+tcrossprod matches naive nested loops exactly.

# Reference (naive) implementations for comparison
.ref_transitions <- function(X, window_size, mode) {
  n <- nrow(X); k <- ncol(X)
  if (n < 2L) return(matrix(0, k, k))
  if (window_size <= 1L) return(crossprod(X[-n, , drop = FALSE], X[-1, , drop = FALSE]))
  weights <- matrix(0, k, k)
  if (mode == "non-overlapping") {
    divides <- n %% window_size == 0L
    q <- n %/% window_size - 1L * divides
    for (i in seq_len(q)) {
      j_idx <- seq((i - 1L) * window_size + 1L, i * window_size)
      k_idx <- seq(i * window_size + 1L, min(n, (i + 1L) * window_size))
      for (j in j_idx) for (ki in k_idx) weights <- weights + tcrossprod(X[j, ], X[ki, ])
    }
  } else {
    n_windows <- n - window_size + 1L
    if (n_windows < 2L) return(weights)
    for (i in seq_len(n_windows - 1L)) {
      j_idx <- seq(i, i + window_size - 1L)
      k_idx <- seq(i + 1L, i + window_size)
      for (j in j_idx) for (ki in k_idx) weights <- weights + tcrossprod(X[j, ], X[ki, ])
    }
  }
  weights
}

.ref_cooccurrence <- function(X, window_size, mode) {
  n <- nrow(X); k <- ncol(X)
  if (window_size <= 1L) return(crossprod(X))
  weights <- matrix(0, k, k)
  if (mode == "non-overlapping") {
    n_windows <- ceiling(n / window_size)
    for (i in seq_len(n_windows)) {
      idx <- seq((i - 1L) * window_size + 1L, min(n, i * window_size))
      for (j in idx) for (ki in idx) weights <- weights + tcrossprod(X[j, ], X[ki, ])
    }
  } else {
    n_windows <- n - window_size + 1L
    if (n_windows < 1L) return(weights)
    for (i in seq_len(n_windows)) {
      idx <- seq(i, i + window_size - 1L)
      for (j in idx) for (ki in idx) weights <- weights + tcrossprod(X[j, ], X[ki, ])
    }
  }
  weights
}

test_that("vectorized transitions match naive loops across 20 random configs", {
  set.seed(999)
  configs <- list(
    list(n = 6, k = 3, ws = 2, mode = "non-overlapping"),
    list(n = 6, k = 3, ws = 2, mode = "overlapping"),
    list(n = 6, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 6, k = 3, ws = 3, mode = "overlapping"),
    list(n = 10, k = 4, ws = 2, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 2, mode = "overlapping"),
    list(n = 10, k = 4, ws = 3, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 3, mode = "overlapping"),
    list(n = 10, k = 4, ws = 5, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 5, mode = "overlapping"),
    list(n = 15, k = 5, ws = 2, mode = "non-overlapping"),
    list(n = 15, k = 5, ws = 2, mode = "overlapping"),
    list(n = 15, k = 5, ws = 4, mode = "non-overlapping"),
    list(n = 15, k = 5, ws = 4, mode = "overlapping"),
    list(n = 20, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 20, k = 3, ws = 3, mode = "overlapping"),
    list(n = 7, k = 3, ws = 2, mode = "non-overlapping"),
    list(n = 7, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 9, k = 4, ws = 4, mode = "overlapping"),
    list(n = 50, k = 6, ws = 5, mode = "non-overlapping")
  )
  for (cfg in configs) {
    X <- matrix(sample(0:1, cfg$n * cfg$k, replace = TRUE), cfg$n, cfg$k)
    ref <- .ref_transitions(X, cfg$ws, cfg$mode)
    vec <- Nestimate:::.wtna_transitions(X, cfg$ws, cfg$mode)
    expect_identical(
      vec, ref,
      info = sprintf("n=%d k=%d ws=%d mode=%s", cfg$n, cfg$k, cfg$ws, cfg$mode)
    )
  }
})

test_that("vectorized cooccurrence matches naive loops across 20 random configs", {
  set.seed(888)
  configs <- list(
    list(n = 6, k = 3, ws = 2, mode = "non-overlapping"),
    list(n = 6, k = 3, ws = 2, mode = "overlapping"),
    list(n = 6, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 6, k = 3, ws = 3, mode = "overlapping"),
    list(n = 10, k = 4, ws = 2, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 2, mode = "overlapping"),
    list(n = 10, k = 4, ws = 3, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 3, mode = "overlapping"),
    list(n = 10, k = 4, ws = 5, mode = "non-overlapping"),
    list(n = 10, k = 4, ws = 5, mode = "overlapping"),
    list(n = 15, k = 5, ws = 2, mode = "non-overlapping"),
    list(n = 15, k = 5, ws = 2, mode = "overlapping"),
    list(n = 15, k = 5, ws = 4, mode = "non-overlapping"),
    list(n = 15, k = 5, ws = 4, mode = "overlapping"),
    list(n = 20, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 20, k = 3, ws = 3, mode = "overlapping"),
    list(n = 7, k = 3, ws = 2, mode = "non-overlapping"),
    list(n = 7, k = 3, ws = 3, mode = "non-overlapping"),
    list(n = 9, k = 4, ws = 4, mode = "overlapping"),
    list(n = 50, k = 6, ws = 5, mode = "non-overlapping")
  )
  for (cfg in configs) {
    X <- matrix(sample(0:1, cfg$n * cfg$k, replace = TRUE), cfg$n, cfg$k)
    ref <- .ref_cooccurrence(X, cfg$ws, cfg$mode)
    vec <- Nestimate:::.wtna_cooccurrence(X, cfg$ws, cfg$mode)
    expect_identical(
      vec, ref,
      info = sprintf("n=%d k=%d ws=%d mode=%s", cfg$n, cfg$k, cfg$ws, cfg$mode)
    )
  }
})

test_that("vectorized wtna end-to-end matches on bundled data", {
  # Use bundled dataset for realistic equivalence
  hw <- head(learning_activities, 30)
  codes <- setdiff(names(hw), "student")
  codes <- codes[vapply(hw[codes], function(x) all(x %in% c(0L, 1L, NA)), logical(1))]
  if (length(codes) < 2) skip("Not enough binary columns")
  X <- as.matrix(hw[, codes, drop = FALSE])
  storage.mode(X) <- "double"
  X[is.na(X)] <- 0

  # Transition ws=3 non-overlapping
  ref <- .ref_transitions(X, 3L, "non-overlapping")
  vec <- Nestimate:::.wtna_transitions(X, 3L, "non-overlapping")
  expect_identical(vec, ref, info = "bundled data transition ws=3")

  # Co-occurrence ws=2 non-overlapping
  ref_co <- .ref_cooccurrence(X, 2L, "non-overlapping")
  vec_co <- Nestimate:::.wtna_cooccurrence(X, 2L, "non-overlapping")
  expect_identical(vec_co, ref_co, info = "bundled data cooccurrence ws=2")
})

# ---- Branch-matrix coverage (task #17) ----
# Crosses method x type x mode x window_size x (actor|NULL). Each combination
# must run without error and produce structurally-valid output. Exists to
# catch regressions where one branch silently starts returning NaN, losing
# dimnames, or producing asymmetric co-occurrence matrices.

test_that("wtna branch matrix: all combinations produce valid output", {
  set.seed(17)
  n_rows <- 30L
  k      <- 4L
  # Small random one-hot-ish matrix (rows pick exactly one code)
  pick <- sample.int(k, n_rows, replace = TRUE)
  X <- matrix(0L, n_rows, k)
  X[cbind(seq_len(n_rows), pick)] <- 1L
  code_names <- paste0("c", seq_len(k))
  colnames(X) <- code_names
  df_single <- as.data.frame(X)

  # Same data with an actor column (two actors, roughly balanced)
  df_multi         <- df_single
  df_multi$actor   <- rep(c("A", "B"), each = n_rows / 2)

  grid <- expand.grid(
    method = c("transition", "cooccurrence", "both"),
    type   = c("frequency", "relative"),
    mode   = c("non-overlapping", "overlapping"),
    ws     = c(1L, 2L, 3L),
    actor  = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  check_one_net <- function(net, expect_symmetric, info) {
    expect_true(is.matrix(net$weights), info = info)
    expect_equal(dim(net$weights), c(k, k), info = info)
    expect_false(any(is.na(net$weights)),  info = paste0(info, " [NA]"))
    expect_false(any(is.nan(net$weights)), info = paste0(info, " [NaN]"))
    expect_false(any(is.infinite(net$weights)), info = paste0(info, " [Inf]"))
    if (expect_symmetric) {
      expect_true(isSymmetric(net$weights), info = paste0(info, " [sym]"))
    }
  }

  for (i in seq_len(nrow(grid))) {
    cfg <- grid[i, ]
    df  <- if (cfg$actor) df_multi else df_single
    info <- sprintf("method=%s type=%s mode=%s ws=%d actor=%s",
                    cfg$method, cfg$type, cfg$mode, cfg$ws, cfg$actor)

    fit <- wtna(
      df,
      method      = cfg$method,
      type        = cfg$type,
      codes       = code_names,
      window_size = cfg$ws,
      mode        = cfg$mode,
      actor       = if (cfg$actor) "actor" else NULL
    )

    # Cooccurrence is only guaranteed symmetric for type="frequency" — the
    # row-normalization in .wtna_finalize (type="relative") deliberately
    # breaks symmetry to produce conditional co-occurrence.
    sym_coocc <- (cfg$type == "frequency")

    if (cfg$method == "both") {
      expect_s3_class(fit, "wtna_mixed")
      check_one_net(fit$transition,   expect_symmetric = FALSE,     info = paste0(info, " [T]"))
      check_one_net(fit$cooccurrence, expect_symmetric = sym_coocc, info = paste0(info, " [C]"))
    } else {
      expect_s3_class(fit, "netobject")
      check_one_net(fit, expect_symmetric = (cfg$method == "cooccurrence" && sym_coocc), info = info)

      # Type-specific invariants
      if (cfg$type == "frequency") {
        expect_true(all(fit$weights >= 0), info = paste0(info, " [nonneg]"))
      } else {
        # Relative: each non-zero row of transitions sums to 1; cooccurrence
        # row-normalizes the same way in .wtna_finalize
        row_sums <- rowSums(fit$weights)
        nz       <- row_sums > 0
        if (any(nz)) {
          expect_equal(
            unname(row_sums[nz]), rep(1, sum(nz)),
            tolerance = 1e-10, info = paste0(info, " [rowsum=1]")
          )
        }
      }
    }
  }
})
