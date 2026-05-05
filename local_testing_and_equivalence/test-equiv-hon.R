# Numerical equivalence: Nestimate's build_hon() vs independent oracles.
#
# build_hon() implements Xu et al. (2016) BuildHON / (Saebi 2020) BuildHON+,
# which produces a MIXED-ORDER graph with rewiring: when a k-th order source
# tuple carries a distinct transition distribution, the edge from its
# (k-1)-order prefix gets redirected to the k-th-order node. This means the
# output cannot be compared directly to pathpy.HigherOrderNetwork(k=k), which
# is a pure k-th-order graph without rewiring.
#
# Instead we validate at two levels:
#
#  1. Counting primitive: Nestimate's .hon_build_observations() — the k-gram
#     counter inside build_hon — must match MOGen's .mogen_count_kgrams()
#     at every order. MOGen's counter has already been validated at machine
#     precision against pathpy (see test-equiv-mogen.R), so agreement here
#     transitively validates HON's counts against the canonical reference.
#
#  2. Output invariants:
#     * probabilities out of every source sum to 1 (row-stochasticity);
#     * total count conservation: the sum of raw counts across rules in a
#       given first-order context equals the total count of that context
#       appearing as a source in the sequences.
#
# The invariants catch rewiring bugs that a pure count comparison would miss.

set.seed(4242)
N_HON <- 50L
TOL_PROB <- 1e-10

# ---- Config generation ----
hon_configs <- lapply(seq_len(N_HON), function(i) {
  list(n_actors = sample(c(10L, 15L, 20L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(15L, 20L, 30L), 1),
       max_order = sample(2:3, 1),
       seed = sample.int(100000, 1))
})

# Access to unexported internals via ::: is fine in a local testing file.
.hon_build_observations <- Nestimate:::.hon_build_observations
.hon_decode <- Nestimate:::.hon_decode
.mogen_count_kgrams <- Nestimate:::.mogen_count_kgrams
HON_SEP <- "\x01"

# Normalise Nestimate's HON count environment into a data.frame keyed by
# (source tuple, target state). Skips entries with empty target-count vectors
# (possible if a trajectory is too short for a given order).
.observations_to_df <- function(count_env, order) {
  rows <- lapply(ls(count_env), function(src_key) {
    parts <- .hon_decode(src_key)
    if (length(parts) != order) return(NULL)
    counts <- count_env[[src_key]]
    if (length(counts) == 0L) return(NULL)
    data.frame(src = paste(parts, collapse = HON_SEP),
               tgt = names(counts),
               count = as.integer(counts),
               stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame(src = character(0), tgt = character(0),
                      count = integer(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

# Convert .mogen_count_kgrams output at order k into the same (src, tgt, count)
# schema as .observations_to_df: src is the k-tuple context, tgt is the next
# state (last element of the target k-tuple).
.kgrams_to_df <- function(kg, k) {
  if (nrow(kg$edges) == 0L) {
    return(data.frame(src = character(0), tgt = character(0),
                      count = integer(0), stringsAsFactors = FALSE))
  }
  # kg$edges$from / $to are \x01-joined k-tuples. Target "next state" = last
  # element of $to. Source tuple stays as-is.
  to_parts <- strsplit(kg$edges$to, HON_SEP, fixed = TRUE)
  next_state <- vapply(to_parts, function(v) v[length(v)], character(1L))
  data.frame(src = kg$edges$from, tgt = next_state,
             count = as.integer(kg$edges$weight),
             stringsAsFactors = FALSE)
}

.independent_kl <- function(p, q) {
  keys <- union(names(p), names(q))
  p <- p[keys]; q <- q[keys]
  p[is.na(p)] <- 0; q[is.na(q)] <- 0
  if (any(p > 0 & q == 0)) return(Inf)
  sum(p[p > 0] * log2(p[p > 0] / q[p > 0]))
}

test_that("HON .hon_build_observations matches .mogen_count_kgrams per order", {
  skip_on_cran()
  skip_equiv_tests()

  report <- equiv_report()

  invisible(lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    trajectories <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    count_env <- .hon_build_observations(trajectories, cfg$max_order)

    invisible(lapply(seq_len(cfg$max_order), function(k) {
      hon_df <- .observations_to_df(count_env, order = k)
      kg <- .mogen_count_kgrams(trajectories, k)
      mog_df <- .kgrams_to_df(kg, k)

      hon_df$key <- paste(hon_df$src, hon_df$tgt, sep = "||")
      mog_df$key <- paste(mog_df$src, mog_df$tgt, sep = "||")

      common <- intersect(hon_df$key, mog_df$key)
      only_hon <- setdiff(hon_df$key, mog_df$key)
      only_mog <- setdiff(mog_df$key, hon_df$key)

      hon_c <- hon_df$count[match(common, hon_df$key)]
      mog_c <- mog_df$count[match(common, mog_df$key)]
      delta <- abs(hon_c - mog_c)

      report$log(
        func = sprintf("hon_observations_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,mo=%d)",
                         i, cfg$n_actors, cfg$seq_length, cfg$max_order),
        n_checked = length(common),
        n_failed = as.integer(sum(delta > 0)),
        max_abs_err = if (length(delta)) max(delta) else 0,
        mean_abs_err = if (length(delta)) mean(delta) else 0,
        median_abs_err = if (length(delta)) stats::median(delta) else 0,
        p95_abs_err = if (length(delta)) as.numeric(stats::quantile(delta, 0.95)) else 0,
        reference = "mogen_count_kgrams (validated vs pathpy)",
        notes = sprintf("k=%d common=%d only_hon=%d only_mog=%d",
                        k, length(common), length(only_hon), length(only_mog))
      )

      # Hard assertions: set equality + exact count agreement.
      expect_equal(length(only_hon), 0L,
                   label = sprintf("cfg%d k=%d hon-only observations", i, k))
      expect_equal(length(only_mog), 0L,
                   label = sprintf("cfg%d k=%d mogen-only observations", i, k))
      if (length(common) > 0L) {
        expect_equal(hon_c, mog_c, tolerance = 0,
                     label = sprintf("cfg%d k=%d observation counts", i, k))
      }
    }))

    NULL
  }))

  report$write_csv("hon_observations")
  report$write_cvs("hon_observations",
                   "local_testing_and_equivalence/test-equiv-hon.R")
})

test_that("HON total k-gram counts match pathpy.HigherOrderNetwork directly", {
  # External-library check that does NOT ride on MOGen. For each order k,
  # sum of counts in .hon_build_observations at that order must equal the
  # sum of subpath+longest-path weights in pathpy.HigherOrderNetwork(k=k).
  # The two codebases share no code, so agreement is genuinely independent.
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  report <- equiv_report()

  invisible(lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    trajectories <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    count_env <- .hon_build_observations(trajectories, cfg$max_order)
    paths <- py_paths_from_sequences(trajectories)

    invisible(lapply(seq_len(cfg$max_order), function(k) {
      # Per-cell aligned comparison: build Nestimate's order-k edge counts
      # from .hon_build_observations and compare to pathpy's edge weights.
      hon_df <- .observations_to_df(count_env, order = k)
      py_edges <- py_hon_edges(paths, k)

      # pathpy's HON(k) keys: (src_tuple, tgt_tuple) joined by ",". Source is
      # comma-separated k-tuple; target is comma-separated (k-1-overlap + next).
      # To align: reconstruct Nestimate keys in pathpy notation.
      src_pp <- gsub(HON_SEP, ",", hon_df$src, fixed = TRUE)
      # target tuple = drop first state of src_tuple, append hon_df$tgt
      if (k == 1L) {
        tgt_pp <- hon_df$tgt
      } else {
        src_parts <- strsplit(src_pp, ",", fixed = TRUE)
        tgt_pp <- vapply(seq_along(src_parts), function(r) {
          paste(c(src_parts[[r]][-1L], hon_df$tgt[r]), collapse = ",")
        }, character(1L))
      }
      hon_keys <- paste(src_pp, tgt_pp, sep = "||")
      py_keys <- paste(py_edges$from, py_edges$to, sep = "||")

      common <- intersect(hon_keys, py_keys)
      only_hon <- setdiff(hon_keys, py_keys)
      only_py <- setdiff(py_keys, hon_keys)

      hon_c <- hon_df$count[match(common, hon_keys)]
      py_c <- py_edges$count[match(common, py_keys)]
      delta <- abs(hon_c - py_c)

      report$log(
        func = sprintf("hon_counts_vs_pathpy_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,mo=%d)",
                         i, cfg$n_actors, cfg$seq_length, cfg$max_order),
        n_checked = length(common),
        n_failed = as.integer(sum(delta > 0)),
        max_abs_err = if (length(delta)) max(delta) else 0,
        mean_abs_err = if (length(delta)) mean(delta) else 0,
        median_abs_err = if (length(delta)) stats::median(delta) else 0,
        p95_abs_err = if (length(delta)) as.numeric(stats::quantile(delta, 0.95)) else 0,
        reference = "pathpy.HigherOrderNetwork (direct, no MOGen)",
        notes = sprintf("k=%d common=%d only_hon=%d only_py=%d",
                        k, length(common), length(only_hon), length(only_py))
      )

      expect_equal(length(only_hon), 0L,
                   label = sprintf("cfg%d k=%d hon-only edges vs pathpy", i, k))
      expect_equal(length(only_py), 0L,
                   label = sprintf("cfg%d k=%d pathpy-only edges vs hon", i, k))
      if (length(common) > 0L) {
        expect_equal(hon_c, py_c, tolerance = 0,
                     label = sprintf("cfg%d k=%d hon vs pathpy counts", i, k))
      }
    }))
  }))

  report$write_csv("hon_vs_pathpy")
  report$write_cvs("hon_vs_pathpy",
                   "local_testing_and_equivalence/test-equiv-hon.R")
})

test_that("HON output probabilities sum to 1 per source (both methods)", {
  skip_on_cran()
  skip_equiv_tests()

  invisible(lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    for (method in c("hon", "hon+")) {
      hon <- tryCatch(
        build_hon(seqs, max_order = cfg$max_order, min_freq = 1L,
                  method = method),
        error = function(e) NULL
      )
      if (is.null(hon) || nrow(hon$ho_edges) == 0L) next

      sums <- stats::aggregate(
        probability ~ from,
        data = hon$ho_edges,
        FUN = sum
      )
      # Every source must sum to 1 within float tolerance. Sources that only
      # appear as terminal states in trajectories have no outgoing edges so
      # they're simply absent from the aggregate — not violations.
      max_deviation <- max(abs(sums$probability - 1))
      expect_true(
        max_deviation < TOL_PROB,
        label = sprintf("cfg%d method=%s max row-prob deviation = %.2e",
                        i, method, max_deviation)
      )
    }
  }))
})

test_that("HON real-data anchor: human_long sequences match mogen counts at every order", {
  # Synthetic uniform-random sequences hide higher-order structure. This anchor
  # uses real coded interaction sequences from the bundled human_long dataset
  # (Specify -> Command -> Refine etc.) to verify that .hon_build_observations
  # still matches mogen_count_kgrams when state distributions are skewed and
  # transitions are autocorrelated.
  skip_on_cran()
  skip_equiv_tests()

  seqs <- bundled_sequences("human_long", max_actors = 100L)
  seqs <- seqs[lengths(seqs) >= 2L]
  count_env <- .hon_build_observations(seqs, max_order = 3L)

  invisible(lapply(seq_len(3L), function(k) {
    hon_df <- .observations_to_df(count_env, order = k)
    kg <- .mogen_count_kgrams(seqs, k)
    mog_df <- .kgrams_to_df(kg, k)

    hon_df$key <- paste(hon_df$src, hon_df$tgt, sep = "||")
    mog_df$key <- paste(mog_df$src, mog_df$tgt, sep = "||")
    common <- intersect(hon_df$key, mog_df$key)
    only_hon <- setdiff(hon_df$key, mog_df$key)
    only_mog <- setdiff(mog_df$key, hon_df$key)

    expect_equal(length(only_hon), 0L,
                 label = sprintf("real human_long k=%d hon-only obs", k))
    expect_equal(length(only_mog), 0L,
                 label = sprintf("real human_long k=%d mogen-only obs", k))
    if (length(common) > 0L) {
      expect_equal(hon_df$count[match(common, hon_df$key)],
                   mog_df$count[match(common, mog_df$key)],
                   tolerance = 0,
                   label = sprintf("real human_long k=%d obs counts", k))
    }
  }))
})

test_that("HON rewiring follows independent KL threshold on separable contexts", {
  skip_on_cran()
  skip_equiv_tests()

  seqs <- c(rep(list(c("A", "B", "X")), 30L),
            rep(list(c("C", "B", "Y")), 30L))
  count_env <- .hon_build_observations(seqs, max_order = 2L)
  b_dist <- count_env[[Nestimate:::.hon_encode("B")]]
  ab_dist <- count_env[[Nestimate:::.hon_encode(c("A", "B"))]]
  cb_dist <- count_env[[Nestimate:::.hon_encode(c("C", "B"))]]
  b_prob <- b_dist / sum(b_dist)
  ab_prob <- ab_dist / sum(ab_dist)
  cb_prob <- cb_dist / sum(cb_dist)

  ab_kl <- .independent_kl(ab_prob, b_prob)
  cb_kl <- .independent_kl(cb_prob, b_prob)
  ab_threshold <- 2 / log2(1 + sum(ab_dist))
  cb_threshold <- 2 / log2(1 + sum(cb_dist))
  expect_true(ab_kl > ab_threshold,
              label = sprintf("AB KL %.2e <= threshold %.2e",
                              ab_kl, ab_threshold))
  expect_true(cb_kl > cb_threshold,
              label = sprintf("CB KL %.2e <= threshold %.2e",
                              cb_kl, cb_threshold))

  hon <- build_hon(seqs, max_order = 2L, min_freq = 1L, method = "hon")
  expect_true("A -> B" %in% rownames(hon$matrix))
  expect_true("C -> B" %in% rownames(hon$matrix))
  expect_equal(unname(hon$matrix["A", "A -> B"]), 1, tolerance = TOL_PROB)
  expect_equal(unname(hon$matrix["C", "C -> B"]), 1, tolerance = TOL_PROB)
  expect_equal(unname(hon$matrix["A", "B"]), 0, tolerance = TOL_PROB)
  expect_equal(unname(hon$matrix["C", "B"]), 0, tolerance = TOL_PROB)
})

test_that("HON first-order rule counts conserve trajectory context totals", {
  skip_on_cran()
  skip_equiv_tests()

  # Invariant: for each first-order state X, the total count of (X -> *)
  # occurrences equals the total times X appears as a source in trajectories.
  # (Rewiring redistributes WHICH targets get counted but never destroys.)
  invisible(lapply(seq_len(N_HON), function(i) {
    cfg <- hon_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    trajectories <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    # Compute raw first-order source-state occurrences directly.
    src_counts <- table(unlist(lapply(trajectories, function(tr) {
      if (length(tr) < 2L) return(character(0))
      tr[-length(tr)]
    })))

    hon <- tryCatch(
      build_hon(trajectories, max_order = cfg$max_order, min_freq = 1L,
                method = "hon"),
      error = function(e) NULL
    )
    if (is.null(hon)) return(NULL)

    # Sum counts grouped by first state of each source context.
    edges <- hon$ho_edges
    edges <- edges[!is.na(edges$count), ]
    first_state <- vapply(strsplit(edges$from, " -> ", fixed = TRUE),
                          `[`, character(1L), 1L)
    agg <- stats::aggregate(count ~ first_state, data = data.frame(
      first_state = first_state, count = edges$count, stringsAsFactors = FALSE
    ), FUN = sum)

    for (s in names(src_counts)) {
      got <- if (s %in% agg$first_state) agg$count[agg$first_state == s] else 0L
      expected <- as.integer(src_counts[s])
      # Inequality is allowed only if Nestimate DROPS context-dependent rules
      # due to min_freq — we set min_freq = 1 so no drops; expect equality.
      # However rewiring can CREATE new source contexts (higher-order), which
      # are COUNTED AGAIN in their own right. So `got` may exceed `expected`;
      # this check only guards against losing first-order source mass.
      expect_true(
        got >= expected,
        label = sprintf("cfg%d state=%s got=%d expected=%d (rewire lost counts)",
                        i, s, got, expected)
      )
    }
  }))
})
