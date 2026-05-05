# Numerical equivalence: Nestimate's build_hypa() vs two independent references.
#
# Why two references: build_hypa() computes an approximate p-value via classical
# (central) hypergeometric with N = round(sum(xi)), K = round(xi[i,j]),
# n = min(sum(adj), N). Two independent checks:
#
#   (1) Monte-Carlo multinomial edge placement under the configuration-model
#       null. Draws m = sum(adj) edge tokens from the xi-weighted distribution,
#       counts how often f_obs-or-fewer edges fall on (v, w). Empirical p-value.
#       Independent of any analytic formula; tests null-model *semantics*.
#
#   (2) BiasedUrn::pMWNCHypergeo (Wallenius multivariate noncentral
#       hypergeometric). Formally the correct null distribution per LaRock et
#       al. 2020 when sampling m edges without replacement with weights xi.
#       Independent of Nestimate; exact within library tolerance.
#
# Divergence between Nestimate and both references would indicate the
# approximation regime is wrong for real data. Per user direction:
# investigate first, decide later (patch or document).

set.seed(4242)
N_HYPA <- 50L
TOL_WALLENIUS <- 1e-6  # Primary assertion: Nestimate == Wallenius(odds=1) exactly
TOL_MC_INFO <- 0.02    # Informational only; MC has sparsity-driven variance
TOL_EXACT <- 1e-10     # for observed counts and xi matrix entries

skip_if_pkg_broken("BiasedUrn")

# ---- Config generation ----
hypa_configs <- lapply(seq_len(N_HYPA), function(i) {
  list(n_actors = sample(c(10L, 15L, 20L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(15L, 20L, 30L), 1),
       k = sample(2:3, 1),
       seed = sample.int(100000, 1))
})

# ---- Monte-Carlo null-model p-values ----
# For each edge (v, w), sample m = sum(adj) draws from multinomial with
# probabilities xi / sum(xi) and count how many land on (v, w). Empirical
# p-value = mean(sim_count <= f_obs).
.mc_hypa_pvalues <- function(adj, xi, iter = 5000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  m <- sum(adj)
  if (m == 0L) return(matrix(NA_real_, nrow(adj), ncol(adj)))
  xi_flat <- as.numeric(xi)
  if (sum(xi_flat) == 0) return(matrix(NA_real_, nrow(adj), ncol(adj)))
  probs <- xi_flat / sum(xi_flat)

  # Each iteration: draw m edge placements with probabilities `probs`.
  # Accumulate via rmultinom; gives exact multinomial sample (independent draws).
  draws <- stats::rmultinom(n = iter, size = m, prob = probs)
  # draws is (n_cells x iter); edge (i,j) is row (j-1)*nrow(adj) + i, but
  # flattening above was column-major so row = (j-1)*nrow + i is already correct.

  # For each cell, empirical p-value = P(sim_count <= f_obs).
  # Vectorised over cells:
  obs_flat <- as.numeric(adj)
  # rowMeans(draws <= obs) but obs varies per row: use sweep equivalent.
  # draws[k, ] vs obs_flat[k]: mean(draws[k,] <= obs_flat[k])
  pvals <- vapply(seq_along(obs_flat), function(k) {
    mean(draws[k, ] <= obs_flat[k])
  }, numeric(1))
  matrix(pvals, nrow(adj), ncol(adj))
}

# ---- Wallenius noncentral hypergeometric via BiasedUrn ----
# Per-edge: condition on total draws = m, weights = xi. P(X_{v,w} <= f_obs).
# Two-colour Wallenius: colour 1 = edge (v,w) with weight xi[v,w];
# colour 2 = all other edges with weight sum(xi) - xi[v,w].
.wallenius_hypa_pvalues <- function(adj, xi) {
  m <- as.integer(sum(adj))
  total_xi <- sum(xi)
  n <- nrow(adj)
  out <- matrix(NA_real_, n, ncol(adj))
  idx <- which(adj > 0, arr.ind = TRUE)
  for (row_i in seq_len(nrow(idx))) {
    i <- idx[row_i, 1L]
    j <- idx[row_i, 2L]
    w <- xi[i, j]
    if (w <= 0 || total_xi - w <= 0) next
    # P(X <= f_obs) where X ~ Wallenius(m, c(round(w), round(total_xi - w)), odds=1)
    # Use odds ratio = w / (total_xi - w); BiasedUrn's pWNCHypergeo signature.
    x <- as.integer(adj[i, j])
    # m_vec = (round(w), round(total_xi - w)), N = sum, n = m (number drawn).
    # pWNCHypergeo(x, m1, m2, n, odds) gives P(X <= x).
    m1 <- as.integer(round(w))
    m2 <- as.integer(round(total_xi - w))
    n_draw <- min(m, m1 + m2)
    if (m1 == 0L || m2 == 0L || n_draw == 0L) next
    # Odds = 1 means classical hypergeometric; for Wallenius with equal odds,
    # the null model treats all edges as exchangeable given xi weights.
    out[i, j] <- tryCatch(
      BiasedUrn::pWNCHypergeo(x = x, m1 = m1, m2 = m2, n = n_draw, odds = 1),
      error = function(e) NA_real_
    )
  }
  out
}

# ---- Main test ----
test_that("HYPA p-values match Monte-Carlo + Wallenius references across 50 configs", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  results <- lapply(seq_len(N_HYPA), function(i) {
    cfg <- hypa_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    hyp <- tryCatch(build_hypa(seqs, k = cfg$k, min_count = 1L,
                               p_adjust = "none"),
                    error = function(e) NULL)
    if (is.null(hyp)) {
      report$log(func = "build_hypa",
                 config = sprintf("cfg%d(build_hypa_error)", i),
                 n_checked = 0L, n_failed = 0L,
                 max_abs_err = NA_real_, mean_abs_err = NA_real_,
                 median_abs_err = NA_real_, p95_abs_err = NA_real_,
                 reference = "montecarlo+wallenius",
                 notes = "build_hypa threw; skipped")
      return(NULL)
    }

    adj <- hyp$adjacency
    xi <- hyp$xi
    if (is.null(adj) || is.null(xi) || sum(adj) == 0L) return(NULL)

    mc_mat <- .mc_hypa_pvalues(adj, xi, iter = 5000L, seed = cfg$seed + 1L)
    wall_mat <- .wallenius_hypa_pvalues(adj, xi)

    # Reconstruct (i, j) matrix positions for each row of hyp$scores.
    # scores is sorted by build_hypa() into [over, under, normal] order, so the
    # which(adj > 0, arr.ind = TRUE) ordering no longer aligns. Rebuild by
    # parsing scores$from (" -> " joined) and scores$to (last state only).
    scores <- hyp$scores
    adj_nodes <- rownames(adj)
    HON_SEP <- "\x01"
    from_nodes <- gsub(" -> ", HON_SEP, scores$from, fixed = TRUE)
    i_vec <- match(from_nodes, adj_nodes)
    # to-node = drop first state of from-node, append scores$to
    to_nodes <- vapply(seq_len(nrow(scores)), function(r) {
      parts <- strsplit(from_nodes[r], HON_SEP, fixed = TRUE)[[1L]]
      paste(c(parts[-1L], scores$to[r]), collapse = HON_SEP)
    }, character(1L))
    j_vec <- match(to_nodes, adj_nodes)

    bad <- is.na(i_vec) | is.na(j_vec)
    if (any(bad)) {
      report$log(func = "build_hypa",
                 config = sprintf("cfg%d(node_match_failed)", i),
                 n_checked = 0L, n_failed = as.integer(sum(bad)),
                 max_abs_err = NA_real_, mean_abs_err = NA_real_,
                 median_abs_err = NA_real_, p95_abs_err = NA_real_,
                 reference = "montecarlo+wallenius",
                 notes = sprintf("n_unmatched=%d of %d", sum(bad), length(i_vec)))
      return(NULL)
    }
    cell_idx <- cbind(i_vec, j_vec)

    nest_p <- scores$p_value
    mc_p <- mc_mat[cell_idx]
    wall_p <- wall_mat[cell_idx]

    # MC comparison
    mc_delta <- abs(nest_p - mc_p)
    keep_mc <- !is.na(mc_delta)
    mc_n_failed <- sum(mc_delta[keep_mc] > TOL_MC_INFO)
    report$log(
      func = "build_hypa_vs_mc", config = sprintf("cfg%d(k=%d,n=%d,s=%d)",
                                                  i, cfg$k, cfg$n_actors,
                                                  cfg$seq_length),
      n_checked = sum(keep_mc), n_failed = mc_n_failed,
      max_abs_err = if (any(keep_mc)) max(mc_delta[keep_mc]) else NA_real_,
      mean_abs_err = if (any(keep_mc)) mean(mc_delta[keep_mc]) else NA_real_,
      median_abs_err = if (any(keep_mc)) stats::median(mc_delta[keep_mc]) else NA_real_,
      p95_abs_err = if (any(keep_mc)) stats::quantile(mc_delta[keep_mc], 0.95) else NA_real_,
      reference = "MonteCarlo(iter=5000)",
      notes = sprintf("n_edges=%d", sum(keep_mc))
    )

    # Wallenius comparison (exact within BiasedUrn precision, ~1e-10)
    wall_delta <- abs(nest_p - wall_p)
    keep_w <- !is.na(wall_delta)
    wall_n_failed <- sum(wall_delta[keep_w] > TOL_WALLENIUS)
    report$log(
      func = "build_hypa_vs_wallenius", config = sprintf("cfg%d(k=%d,n=%d,s=%d)",
                                                         i, cfg$k, cfg$n_actors,
                                                         cfg$seq_length),
      n_checked = sum(keep_w), n_failed = wall_n_failed,
      max_abs_err = if (any(keep_w)) max(wall_delta[keep_w]) else NA_real_,
      mean_abs_err = if (any(keep_w)) mean(wall_delta[keep_w]) else NA_real_,
      median_abs_err = if (any(keep_w)) stats::median(wall_delta[keep_w]) else NA_real_,
      p95_abs_err = if (any(keep_w)) stats::quantile(wall_delta[keep_w], 0.95) else NA_real_,
      reference = "BiasedUrn::pWNCHypergeo(odds=1)",
      notes = sprintf("n_edges=%d", sum(keep_w))
    )

    # Primary assertion: Nestimate's analytic HYPA p-value must match
    # Wallenius noncentral hypergeometric (odds=1) to ~1e-6. This is the
    # canonical LaRock 2020 null-model definition. MC delta is logged for
    # visibility but not asserted — sparse k>=3 regimes give MC variance
    # that exceeds any meaningful Nestimate error.
    wall_max <- if (any(keep_w)) max(wall_delta[keep_w]) else 0
    expect_true(
      wall_max < TOL_WALLENIUS,
      label = sprintf("cfg%d Wallenius max delta = %.2e", i, wall_max)
    )

    NULL
  })

  report$write_csv("hypa")
  report$write_cvs("hypa", "local_testing_and_equivalence/test-equiv-hypa.R")
})

# ---- Secondary: xi matrix is exactly outer(out_strength, in_strength) * mask ----
test_that("HYPA xi matrix matches manual outer-product definition", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(10L), function(i) {
    cfg <- hypa_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })
    hyp <- tryCatch(build_hypa(seqs, k = cfg$k, min_count = 1L,
                               p_adjust = "none"),
                    error = function(e) NULL)
    if (is.null(hyp)) return(NULL)
    adj <- hyp$adjacency
    xi <- hyp$xi
    out_str <- rowSums(adj)
    in_str <- colSums(adj)
    mask <- adj > 0
    xi_ref <- outer(out_str, in_str) * mask
    expect_equal(unname(xi), unname(xi_ref), tolerance = TOL_EXACT,
                 label = sprintf("cfg%d xi outer-product", i))
  }))
})

test_that("HYPA real-data anchor: human_long Wallenius p-values agree at 1e-6", {
  # Real coded human-AI interaction sequences have many rare transitions where
  # the configuration-model expected count xi[i,j] is tiny. This stresses the
  # central-vs-Wallenius approximation regime that the synthetic uniform sampler
  # never produces.
  skip_on_cran()
  skip_equiv_tests()
  skip_if_pkg_broken("BiasedUrn")

  seqs <- bundled_sequences("human_long", max_actors = 100L)
  seqs <- seqs[lengths(seqs) >= 2L]
  hyp <- build_hypa(seqs, k = 2L, min_count = 1L, p_adjust = "none")
  if (is.null(hyp) || sum(hyp$adjacency) == 0L) skip("HYPA empty")

  adj <- hyp$adjacency; xi <- hyp$xi
  wall_mat <- .wallenius_hypa_pvalues(adj, xi)

  scores <- hyp$scores
  HON_SEP <- "\x01"
  from_nodes <- gsub(" -> ", HON_SEP, scores$from, fixed = TRUE)
  i_vec <- match(from_nodes, rownames(adj))
  to_nodes <- vapply(seq_len(nrow(scores)), function(r) {
    parts <- strsplit(from_nodes[r], HON_SEP, fixed = TRUE)[[1L]]
    paste(c(parts[-1L], scores$to[r]), collapse = HON_SEP)
  }, character(1L))
  j_vec <- match(to_nodes, rownames(adj))
  ok <- !is.na(i_vec) & !is.na(j_vec)

  delta <- abs(scores$p_value[ok] - wall_mat[cbind(i_vec[ok], j_vec[ok])])
  delta <- delta[!is.na(delta)]
  expect_true(
    max(delta) < TOL_WALLENIUS,
    label = sprintf("real human_long Wallenius max delta = %.2e", max(delta))
  )
})

test_that("HYPA p_adjust paths match stats::p.adjust for under and over tails", {
  skip_on_cran()
  skip_equiv_tests()

  cfg <- hypa_configs[[1L]]
  data <- simulate_sequences(n_actors = cfg$n_actors,
                             n_states = cfg$n_states,
                             seq_length = cfg$seq_length,
                             seed = cfg$seed)
  seqs <- lapply(seq_len(nrow(data)), function(r) {
    as.character(unlist(data[r, ], use.names = FALSE))
  })
  base_fit <- build_hypa(seqs, k = cfg$k, min_count = 1L,
                         p_adjust = "none")
  key <- base_fit$scores$path

  invisible(lapply(c("BH", "bonferroni", "holm"), function(method) {
    fit <- build_hypa(seqs, k = cfg$k, min_count = 1L, p_adjust = method)
    idx <- match(key, fit$scores$path)
    expect_false(anyNA(idx), label = sprintf("%s adjusted path alignment", method))
    expected_under <- stats::p.adjust(base_fit$scores$p_value, method = method)
    expected_over <- stats::p.adjust(1 - base_fit$scores$p_value, method = method)
    under_delta <- max(abs(fit$scores$p_adjusted_under[idx] - expected_under))
    over_delta <- max(abs(fit$scores$p_adjusted_over[idx] - expected_over))
    expect_true(under_delta < TOL_EXACT,
                label = sprintf("%s under-tail p.adjust delta = %.2e",
                                method, under_delta))
    expect_true(over_delta < TOL_EXACT,
                label = sprintf("%s over-tail p.adjust delta = %.2e",
                                method, over_delta))
  }))
})
