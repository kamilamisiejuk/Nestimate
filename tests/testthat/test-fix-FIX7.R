# Regression tests for FIX7 (audit findings A10-F01 .. A10-F04).
#
# A10-F01: chain_structure(tol=) must NOT double as the absorbing-state
#          threshold; absorbing = P[i,i] == 1 exactly.
# A10-F02: classification and hitting_probabilities must use the same
#          (tol-thresholded) transition support => internally consistent.
# A10-F03: transition_entropy print/summary footer must describe each
#          normalised column truthfully (redundancy_norm = (Hpi-h)/Hpi).
# A10-F04: @return of markov_order_test / transition_entropy must list
#          every slot the real object contains.

test_that("A10-F01: P with diag 0.99 is NOT absorbing at any tol", {
  Pdiag <- matrix(c(0.99, 0.01, 0,
                    0.30, 0.40, 0.30,
                    0.20, 0.30, 0.50), 3, 3, byrow = TRUE,
                  dimnames = list(c("a", "b", "c"), c("a", "b", "c")))

  # Across a wide tol range the 1%-escape recurrent state stays recurrent
  # and absorbing_states stays empty (the bug flipped it at the tol
  # boundary).
  cls <- vapply(c(1e-10, 1e-3, 0.05, 0.5), function(t) {
    cs <- chain_structure(Pdiag, tol = t)
    paste(cs$classification["a"],
          paste(cs$absorbing_states, collapse = ","), sep = "|")
  }, character(1L))
  expect_true(all(cls == "recurrent|"))

  # The explicit bug repro from the audit must now fail to reproduce.
  c2 <- chain_structure(Pdiag, tol = 0.05)
  expect_false(c2$classification["a"] == "absorbing")
  expect_length(c2$absorbing_states, 0L)
})

test_that("A10-F01: a genuine absorbing state (P[i,i]=1) is still detected", {
  Pabs <- matrix(c(1.0, 0.0, 0.0,
                   0.3, 0.4, 0.3,
                   0.2, 0.3, 0.5), 3, 3, byrow = TRUE,
                 dimnames = list(c("d", "e", "f"), c("d", "e", "f")))
  cs <- chain_structure(Pabs, tol = 1e-10)
  expect_equal(unname(cs$classification["d"]), "absorbing")
  expect_equal(cs$absorbing_states, "d")
  # Absorption analysis on the genuine absorbing state still produced.
  expect_false(is.null(cs$absorption_probabilities))
  expect_false(is.null(cs$mean_absorption_time))
})

test_that("A10-F02: classification and hitting_probabilities are consistent", {
  # State X has only a sub-tol (0.001) escape edge. With tol = 1e-2 the
  # support graph zeroes that edge: X becomes its own closed class.
  # hitting_probabilities must agree -> H[X, others] == 0.
  Pt <- matrix(c(0.999, 0.001, 0,
                 0,     0.5,   0.5,
                 0,     0.5,   0.5), 3, 3, byrow = TRUE,
               dimnames = list(c("X", "Y", "Z"), c("X", "Y", "Z")))
  cs <- chain_structure(Pt, tol = 1e-2)

  H <- cs$hitting_probabilities
  cls <- cs$classification

  # General invariant: any state in a CLOSED (recurrent/absorbing) class
  # has zero hitting probability to every state outside that class.
  closed_states <- unlist(cs$recurrent_classes, use.names = FALSE)
  expect_true(length(closed_states) > 0L)
  consistent <- vapply(cs$recurrent_classes, function(cl) {
    others <- setdiff(cs$states, cl)
    if (length(others) == 0L) return(TRUE)
    all(H[cl, others, drop = FALSE] == 0)
  }, logical(1L))
  expect_true(all(consistent))

  # Concrete: X is closed under tol=1e-2 and cannot reach Y or Z.
  expect_true("X" %in% closed_states)
  expect_equal(unname(H["X", "Y"]), 0)
  expect_equal(unname(H["X", "Z"]), 0)
  # No state classified absorbing while reaching others w.p. 1.
  abs_st <- cs$absorbing_states
  if (length(abs_st) > 0L) {
    expect_true(all(H[abs_st, setdiff(cs$states, abs_st), drop = FALSE] == 0))
  }
})

test_that("A10-F02: default-tol full-support chain hitting unchanged", {
  # Regression guard: on a full-support chain (no sub-tol edges) the
  # tol-thresholded reach equals the raw P>0 reach, so hitting is the
  # documented all-reachable result.
  net <- suppressWarnings(
    build_network(as.data.frame(group_regulation_long), method = "relative"))
  cs <- chain_structure(net)
  H <- cs$hitting_probabilities
  expect_false(any(is.na(H)))
  expect_true(all(H >= 0 & H <= 1))
  # Irreducible full-support chain: every off-diagonal hitting prob is 1.
  if (isTRUE(cs$is_irreducible)) {
    od <- H[upper.tri(H) | lower.tri(H)]
    expect_true(all(abs(od - 1) < 1e-8))
  }
})

test_that("A10-F03: redundancy_norm equals (Hpi - h)/Hpi, footer matches", {
  te <- transition_entropy(
    build_network(as.data.frame(trajectories), method = "relative"),
    base = 2)

  # The actually-stored normalised redundancy is the RELATIVE redundancy,
  # not raw / log_b(n) (the old footer claimed the latter).
  rel <- (te$stationary_entropy - te$entropy_rate) / te$stationary_entropy
  expect_equal(te$redundancy_norm, rel)
  expect_false(isTRUE(all.equal(te$redundancy_norm,
                                te$redundancy / te$max_entropy)))

  # The other normalised columns ARE raw / log_b(n).
  expect_equal(unname(te$entropy_rate_norm),
               te$entropy_rate / te$max_entropy)
  expect_equal(unname(te$stationary_entropy_norm),
               te$stationary_entropy / te$max_entropy)

  # Footer text must mention the relative-redundancy formula and must NOT
  # claim every normalised value is raw / log_b(n).
  print_txt <- paste(capture.output(print(te)), collapse = "\n")
  expect_true(grepl("(H(pi) - h(P)) / H(pi)", print_txt, fixed = TRUE))
  expect_false(grepl("Normalised: raw / log_2(n_states); 0 = deterministic, 1 = uniform.",
                      print_txt, fixed = TRUE))

  sum_txt <- paste(capture.output(print(summary(te))), collapse = "\n")
  expect_true(grepl("relative redundancy (H(pi) - h(P)) / H(pi)",
                    sum_txt, fixed = TRUE))
  expect_false(grepl("Normalised values are raw / log_2(n_states), in [0, 1].",
                      sum_txt, fixed = TRUE))
})

test_that("A10-F04: transition_entropy @return lists every real slot", {
  te <- transition_entropy(
    build_network(as.data.frame(trajectories), method = "relative"))
  # Mirrors the \describe{} block in the roxygen @return.
  documented <- c("row_entropy", "row_entropy_norm", "stationary",
                  "stationary_entropy", "stationary_entropy_norm",
                  "entropy_rate", "entropy_rate_norm", "redundancy",
                  "redundancy_norm", "max_entropy", "base", "states")
  expect_setequal(names(te), documented)
})

test_that("A10-F04: markov_order_test @return lists every real slot/column", {
  r <- markov_order_test(simulate_sequences(8, 4, 14, seed = 1),
                         max_order = 3, n_perm = 15, seed = 1)
  documented_slots <- c("optimal_order", "bic_order", "aic_order",
                        "test_table", "permutation_null", "logliks",
                        "layer_dofs", "transition_matrices", "states",
                        "n_sequences", "n_observations", "n_perm",
                        "alpha", "max_order")
  expect_setequal(names(r), documented_slots)

  documented_cols <- c("order", "loglik", "AIC", "BIC", "df", "g2",
                       "p_permutation", "p_asymptotic", "significant")
  expect_setequal(names(r$test_table), documented_cols)
})
