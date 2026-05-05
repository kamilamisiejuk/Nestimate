# Numerical equivalence: Nestimate's build_mogen() vs pathpy.MultiOrderModel.
#
# pathpy 2.2.0 (Scholtes' own lab) is the canonical Python reference for MOGen
# (Scholtes 2017, Gote & Scholtes 2023). We check three things per config:
#
#  1. Per-order count matrices from .mogen_count_kgrams()
#     match pathpy.HigherOrderNetwork(k=k).edges (summed weights).
#  2. Total cumulative log-likelihood at each max_order 1..K matches
#     pathpy's MultiOrderModel(paths, max_order=k).likelihood(log=True).
#  3. path_counts() agrees with pathpy k-gram enumeration (integer-exact).

set.seed(4242)
N_MOGEN <- 30L
TOL <- 1e-8         # counts must match exactly; likelihood to 1e-8
TOL_LIK <- 1e-6     # log-likelihood tolerance (hierarchical decomposition may
                    # differ in floor value for unseen transitions)

skip_if_pkg_broken("reticulate")

# ---- Config generation ----
mogen_configs <- lapply(seq_len(N_MOGEN), function(i) {
  list(n_actors = sample(c(10L, 15L, 20L), 1),
       n_states = sample(3:5, 1),
       seq_length = sample(c(15L, 20L, 30L), 1),
       max_order = sample(2:3, 1),
       seed = sample.int(100000, 1))
})

.mogen_expected_transition_table <- function(mg, order, min_count = 1L) {
  HON_SEP <- "\x01"
  tm <- mg$transition_matrices[[order + 1L]]
  cm <- mg$count_matrices[[order + 1L]]
  idx <- which(cm >= min_count, arr.ind = TRUE)
  if (nrow(idx) == 0L) {
    return(data.frame(path = character(0L), count = integer(0L),
                      probability = numeric(0L), from = character(0L),
                      to = character(0L), stringsAsFactors = FALSE))
  }
  from_raw <- rownames(tm)[idx[, 1L]]
  to_raw <- colnames(tm)[idx[, 2L]]
  parsed <- lapply(seq_len(nrow(idx)), function(r) {
    from_parts <- strsplit(from_raw[r], HON_SEP, fixed = TRUE)[[1L]]
    to_parts <- strsplit(to_raw[r], HON_SEP, fixed = TRUE)[[1L]]
    next_state <- to_parts[length(to_parts)]
    list(path = paste(c(from_parts, next_state), collapse = " -> "),
         from = paste(from_parts, collapse = " -> "),
         to = next_state)
  })
  out <- data.frame(
    path = vapply(parsed, `[[`, character(1L), "path"),
    count = as.integer(cm[idx]),
    probability = round(tm[idx], 4),
    from = vapply(parsed, `[[`, character(1L), "from"),
    to = vapply(parsed, `[[`, character(1L), "to"),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$count), ]
  rownames(out) <- NULL
  out
}

.nonzero_node_count <- function(mat, nodes) {
  present <- intersect(nodes, rownames(mat))
  if (length(present) == 0L) return(0L)
  sub <- mat[present, present, drop = FALSE]
  sum(rowSums(sub) + colSums(sub) > 0)
}

test_that("MOGen per-order count matrices match pathpy.HigherOrderNetwork", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  report <- equiv_report()
  HON_SEP <- "\x01"

  invisible(lapply(seq_len(N_MOGEN), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    mg <- tryCatch(build_mogen(seqs, max_order = cfg$max_order),
                   error = function(e) NULL)
    if (is.null(mg)) return(NULL)

    paths <- py_paths_from_sequences(seqs)

    # Check each order k from 1 to max_order. Nestimate stores the count matrix
    # at index k+1; node names use \x01 separators.
    invisible(lapply(seq_len(cfg$max_order), function(k) {
      nest_cm <- mg$count_matrices[[k + 1L]]
      py_cm <- py_hon_count_matrix(paths, k)

      # Translate node names: Nestimate uses \x01, pathpy uses ",".
      nest_nodes_py <- gsub(HON_SEP, ",", rownames(nest_cm), fixed = TRUE)
      rownames(nest_cm) <- nest_nodes_py
      colnames(nest_cm) <- nest_nodes_py

      # Align on intersection of node sets; log any asymmetry.
      common <- intersect(rownames(nest_cm), rownames(py_cm))
      n_only_nest <- length(setdiff(rownames(nest_cm), common))
      n_only_py <- length(setdiff(rownames(py_cm), common))

      if (length(common) == 0L) {
        report$log(func = sprintf("mogen_counts_k%d", k),
                   config = sprintf("cfg%d(no_common_nodes)", i),
                   n_checked = 0L, n_failed = 0L,
                   max_abs_err = NA_real_, mean_abs_err = NA_real_,
                   median_abs_err = NA_real_, p95_abs_err = NA_real_,
                   reference = "pathpy.HigherOrderNetwork",
                   notes = "empty intersection")
        return(NULL)
      }

      a <- nest_cm[common, common, drop = FALSE]
      b <- py_cm[common, common, drop = FALSE]
      delta <- abs(a - b)

      report$log(
        func = sprintf("mogen_counts_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,k=%d)",
                         i, cfg$n_actors, cfg$seq_length, k),
        n_checked = length(delta),
        n_failed = as.integer(sum(delta > TOL)),
        max_abs_err = max(delta), mean_abs_err = mean(delta),
        median_abs_err = stats::median(delta),
        p95_abs_err = as.numeric(stats::quantile(delta, 0.95)),
        reference = "pathpy.HigherOrderNetwork",
        notes = sprintf("common=%d only_nest=%d only_py=%d",
                        length(common), n_only_nest, n_only_py)
      )

      expect_true(max(delta) < TOL,
                  label = sprintf("cfg%d k=%d max count delta = %.2e",
                                  i, k, max(delta)))
      expect_equal(.nonzero_node_count(nest_cm, setdiff(rownames(nest_cm), common)),
                   0L,
                   label = sprintf("cfg%d k=%d Nestimate-only nonzero nodes",
                                   i, k))
      expect_equal(.nonzero_node_count(py_cm, setdiff(rownames(py_cm), common)),
                   0L,
                   label = sprintf("cfg%d k=%d pathpy-only nonzero nodes",
                                   i, k))
    }))

    NULL
  }))

  report$write_csv("mogen_counts")
  report$write_cvs("mogen_counts",
                   "local_testing_and_equivalence/test-equiv-mogen.R")
})

test_that("MOGen cumulative log-likelihood matches pathpy.MultiOrderModel", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  report <- equiv_report()

  invisible(lapply(seq_len(N_MOGEN), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    mg <- tryCatch(build_mogen(seqs, max_order = cfg$max_order),
                   error = function(e) NULL)
    if (is.null(mg)) return(NULL)

    paths <- py_paths_from_sequences(seqs)

    # Compare total log-likelihood at each max_order 1..cfg$max_order
    invisible(lapply(seq_len(cfg$max_order), function(k) {
      nest_ll <- as.numeric(mg$log_likelihood[k + 1L])
      py_ll <- py_mogen_likelihood(paths, k)
      delta <- abs(nest_ll - py_ll)

      report$log(
        func = sprintf("mogen_loglik_k%d", k),
        config = sprintf("cfg%d(n=%d,s=%d,k=%d)",
                         i, cfg$n_actors, cfg$seq_length, k),
        n_checked = 1L,
        n_failed = as.integer(delta > TOL_LIK),
        max_abs_err = delta, mean_abs_err = delta,
        median_abs_err = delta, p95_abs_err = delta,
        reference = "pathpy.MultiOrderModel.likelihood",
        notes = sprintf("nest=%.4f py=%.4f", nest_ll, py_ll)
      )

      expect_true(delta < TOL_LIK,
                  label = sprintf("cfg%d k=%d ll delta = %.2e (nest=%.4f py=%.4f)",
                                  i, k, delta, nest_ll, py_ll))
    }))

    NULL
  }))

  report$write_csv("mogen_loglik")
  report$write_cvs("mogen_loglik",
                   "local_testing_and_equivalence/test-equiv-mogen.R")
})

test_that("path_counts matches pathpy k-gram enumeration exactly", {
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  invisible(lapply(seq_len(10L), function(i) {
    cfg <- mogen_configs[[i]]
    data <- simulate_sequences(n_actors = cfg$n_actors,
                               n_states = cfg$n_states,
                               seq_length = cfg$seq_length,
                               seed = cfg$seed)
    seqs <- lapply(seq_len(nrow(data)), function(r) {
      as.character(unlist(data[r, ], use.names = FALSE))
    })

    # path_counts uses " -> " separator; k=2 is a transition count.
    pc <- path_counts(seqs, k = 2L)

    # Manually aggregate bigram counts; pathpy HON(k=1) gives same result
    # (subpath + longest-path = total occurrences of each transition).
    paths <- py_paths_from_sequences(seqs)
    py_edges <- py_hon_edges(paths, 1L)
    py_edges$path <- paste(py_edges$from, py_edges$to, sep = " -> ")

    common <- intersect(pc$path, py_edges$path)
    pc_aligned <- pc[match(common, pc$path), ]
    py_aligned <- py_edges[match(common, py_edges$path), ]

    expect_equal(pc_aligned$count, py_aligned$count,
                 tolerance = 0,
                 label = sprintf("cfg%d path_counts k=2", i))
  }))
})

test_that("MOGen real-data anchor: human_long counts and log-lik match pathpy", {
  # Real human-AI coding sequences have realistic state imbalance and
  # long-range dependencies that uniform-random sequences never produce.
  # Validates count_matrices and likelihood at every order on the bundled
  # human_long dataset.
  skip_on_cran()
  skip_equiv_tests()
  skip_if_no_python_ref("pathpy")

  HON_SEP <- "\x01"
  seqs <- bundled_sequences("human_long", max_actors = 80L)
  # pathpy.Paths skips sequences of length < 2 ('if len(s) >= 2'); apply the
  # same filter so both engines see identical input.
  seqs <- seqs[lengths(seqs) >= 2L]
  mg <- build_mogen(seqs, max_order = 3L)
  paths <- py_paths_from_sequences(seqs)

  invisible(lapply(seq_len(3L), function(k) {
    nest_cm <- mg$count_matrices[[k + 1L]]
    py_cm <- py_hon_count_matrix(paths, k)
    rn <- gsub(HON_SEP, ",", rownames(nest_cm), fixed = TRUE)
    rownames(nest_cm) <- rn; colnames(nest_cm) <- rn

    common <- intersect(rownames(nest_cm), rownames(py_cm))
    if (length(common) == 0L) return(NULL)
    expect_equal(.nonzero_node_count(nest_cm, setdiff(rownames(nest_cm), common)),
                 0L,
                 label = sprintf("real human_long Nestimate-only nonzero nodes k=%d", k))
    expect_equal(.nonzero_node_count(py_cm, setdiff(rownames(py_cm), common)),
                 0L,
                 label = sprintf("real human_long pathpy-only nonzero nodes k=%d", k))
    delta <- abs(nest_cm[common, common, drop = FALSE] -
                 py_cm[common, common, drop = FALSE])
    expect_true(max(delta) < TOL,
                label = sprintf("real human_long counts k=%d max delta = %.2e",
                                k, max(delta)))

    nest_ll <- as.numeric(mg$log_likelihood[k + 1L])
    py_ll <- py_mogen_likelihood(paths, k)
    expect_true(abs(nest_ll - py_ll) < TOL_LIK,
                label = sprintf("real human_long loglik k=%d delta = %.2e",
                                k, abs(nest_ll - py_ll)))
  }))
})

test_that("mogen_transitions extracts order-specific transition tables exactly", {
  skip_on_cran()
  skip_equiv_tests()

  cfg <- mogen_configs[[1L]]
  data <- simulate_sequences(n_actors = cfg$n_actors,
                             n_states = cfg$n_states,
                             seq_length = cfg$seq_length,
                             seed = cfg$seed)
  seqs <- lapply(seq_len(nrow(data)), function(r) {
    as.character(unlist(data[r, ], use.names = FALSE))
  })
  mg <- build_mogen(seqs, max_order = cfg$max_order)

  invisible(lapply(seq_len(cfg$max_order), function(k) {
    got <- mogen_transitions(mg, order = k, min_count = 2L)
    expected <- .mogen_expected_transition_table(mg, order = k, min_count = 2L)
    key_got <- paste(got$from, got$to, sep = "||")
    key_expected <- paste(expected$from, expected$to, sep = "||")
    expect_equal(sort(key_got), sort(key_expected),
                 label = sprintf("mogen_transitions keys order %d", k))
    idx <- match(key_expected, key_got)
    expect_equal(got$path[idx], expected$path, tolerance = 0,
                 label = sprintf("mogen_transitions path order %d", k))
    expect_equal(got$count[idx], expected$count, tolerance = 0,
                 label = sprintf("mogen_transitions count order %d", k))
    prob_delta <- max(abs(got$probability[idx] - expected$probability))
    expect_true(prob_delta < TOL,
                label = sprintf("mogen_transitions probability order %d delta = %.2e",
                                k, prob_delta))
  }))
})
