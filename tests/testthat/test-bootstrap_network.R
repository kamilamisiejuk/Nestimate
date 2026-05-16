testthat::skip_on_cran()

# ---- bootstrap_network() Tests ----

# Helper: generate wide sequence data
.make_boot_wide <- function(n = 50, t = 10, states = c("A", "B", "C"),
                            seed = 42) {
  set.seed(seed)
  mat <- matrix(sample(states, n * t, replace = TRUE), nrow = n, ncol = t)
  colnames(mat) <- paste0("T", seq_len(t))
  as.data.frame(mat, stringsAsFactors = FALSE)
}

# Helper: generate frequency-like data for association methods
.make_boot_assoc <- function(n = 100, p = 5, seed = 42) {
  set.seed(seed)
  mat <- matrix(rpois(n * p, lambda = 10), nrow = n, ncol = p)
  colnames(mat) <- paste0("state_", seq_len(p))
  as.data.frame(mat)
}


# ---- Basic functionality: transition methods ----

test_that("bootstrap_network works with method='relative'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 30L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_s3_class(boot$original, "netobject")
  expect_s3_class(boot$model, "netobject")
  expect_equal(boot$method, "relative")
  expect_equal(boot$iter, 30L)
  expect_true(is.matrix(boot$mean))
  expect_true(is.matrix(boot$sd))
  expect_true(is.matrix(boot$p_values))
  expect_true(is.matrix(boot$ci_lower))
  expect_true(is.matrix(boot$ci_upper))
  expect_true(is.matrix(boot$significant))
  expect_equal(nrow(boot$mean), 3)
  expect_equal(ncol(boot$mean), 3)
})

test_that("bootstrap_network works with method='frequency'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "frequency"),
                            iter = 30L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "frequency")
  expect_true(all(boot$mean >= 0))
})

test_that("bootstrap_network warns for one-sequence transition networks", {
  long <- data.frame(action = c("A", "B", "C", "A"),
                     stringsAsFactors = FALSE)
  net <- suppressWarnings(build_tna(long, action = "action"))

  expect_warning(
    boot <- bootstrap_network(net, iter = 2L, seed = 1),
    "one long sequence is not recommended"
  )
  expect_s3_class(boot, "net_bootstrap")
})

test_that("bootstrap_network onehot frequency uses wtna per-sequence counts", {
  df <- data.frame(
    actor = rep(c("s1", "s2", "s3"), each = 4),
    A = c(1L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 1L),
    B = c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 1L, 0L),
    C = c(0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L)
  )
  codes <- c("A", "B", "C")
  net <- suppressWarnings(build_network(
    df, method = "frequency", format = "onehot", codes = codes,
    actor = "actor", window_size = 1L
  ))

  resampling_data <- Nestimate:::.resampling_transition_data(
    net$data, net$metadata, net$params
  )
  pre <- Nestimate:::.precompute_per_sequence(
    resampling_data, net$method, net$params, net$nodes$label
  )
  counts <- matrix(colSums(pre), nrow = 3L, byrow = TRUE,
                   dimnames = list(codes, codes))
  boot <- bootstrap_network(net, iter = 20L, seed = 1)

  expect_equal(counts, net$weights)
  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "frequency")
  expect_true(any(boot$mean > 0))
})

test_that("bootstrap_network works with method='co_occurrence'", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "co_occurrence"),
                            iter = 30L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "co_occurrence")
  expect_false(boot$original$directed)
})


# ---- Basic functionality: association methods ----

test_that("bootstrap_network works with method='cor'", {
  df <- .make_boot_assoc()
  boot <- bootstrap_network(build_network(df, method = "cor"),
                            iter = 20L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "cor")
  expect_false(boot$original$directed)
  expect_true(is.matrix(boot$mean))
  expect_equal(nrow(boot$mean), 5)
})

test_that("bootstrap_network works with method='pcor'", {
  df <- .make_boot_assoc(n = 80, p = 4)
  boot <- bootstrap_network(build_network(df, method = "pcor"),
                            iter = 20L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "pcor")
})

test_that("bootstrap_network works with method='glasso'", {
  df <- .make_boot_assoc(n = 80, p = 4)
  boot <- bootstrap_network(build_network(df, method = "glasso",
                                          params = list(nlambda = 20L)),
                            iter = 20L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "glasso")
})


# ---- Inference ----

test_that("stability inference produces valid p-values", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 50L, seed = 42, inference = "stability")

  expect_true(all(boot$p_values >= 0))
  expect_true(all(boot$p_values <= 1))
  expect_equal(boot$inference, "stability")
  # CR bounds should exist
  expect_true(is.matrix(boot$cr_lower))
  expect_true(is.matrix(boot$cr_upper))
})

test_that("threshold inference produces valid p-values", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 50L, seed = 42, inference = "threshold")

  expect_true(all(boot$p_values >= 0))
  expect_true(all(boot$p_values <= 1))
  expect_equal(boot$inference, "threshold")
  # edge_threshold should be auto-set

  expect_true(is.numeric(boot$edge_threshold))
})


# ---- CI correctness ----

test_that("ci_lower <= ci_upper everywhere", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 50L, seed = 1)

  expect_true(all(boot$ci_lower <= boot$ci_upper + 1e-10))
})


# ---- Summary ----

test_that("summary returns correct columns", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 30L, seed = 1)

  s <- summary(boot)
  expect_true(is.data.frame(s))
  expect_true(all(c("from", "to", "weight", "mean", "sd", "p_value",
                     "sig", "ci_lower", "ci_upper") %in% names(s)))

  # Stability inference should include CR columns
  expect_true(all(c("cr_lower", "cr_upper") %in% names(s)))

  # All edges should have non-zero original weight
  expect_true(all(s$weight != 0))
})

test_that("summary for undirected keeps only upper triangle", {
  df <- .make_boot_assoc()
  boot <- bootstrap_network(build_network(df, method = "cor"),
                            iter = 20L, seed = 1)
  s <- summary(boot)

  if (nrow(s) > 0) {
    expect_true(all(s$from < s$to))
  }
})


# ---- Pruned model ----

test_that("pruned model has <= edges of original", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 50L, seed = 1)

  expect_true(boot$model$n_edges <= boot$original$n_edges)
  # Pruned matrix should be zero where not significant
  non_sig <- boot$p_values >= boot$ci_level
  expect_true(all(boot$model$weights[non_sig] == 0))
})


# ---- Composability ----

test_that("original matches standalone estimate", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative", params = list(format = "wide"))
  boot <- bootstrap_network(net, iter = 20L, seed = 1)

  # Original should match standalone call
  expect_equal(boot$original$weights, net$weights)
})


# ---- Reproducibility ----

test_that("same seed produces identical results", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative")
  boot1 <- bootstrap_network(net, iter = 30L, seed = 99)
  boot2 <- bootstrap_network(net, iter = 30L, seed = 99)

  expect_equal(boot1$p_values, boot2$p_values)
  expect_equal(boot1$mean, boot2$mean)
  expect_equal(boot1$ci_lower, boot2$ci_lower)
})


# ---- Print ----

test_that("print.net_bootstrap produces expected output", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 30L, seed = 1)

  out <- capture.output(print(boot))
  expect_true(any(grepl("Bootstrap Network", out)))
  expect_true(any(grepl("Iterations", out)))
  expect_true(any(grepl("Nodes", out)))
  expect_true(any(grepl("Edges", out)))
})


# ---- Custom estimator ----

test_that("bootstrap works with custom estimator", {
  # Register a simple custom estimator
  custom_fn <- function(data, ...) {
    numeric_cols <- vapply(data, is.numeric, logical(1))
    mat <- as.matrix(data[, numeric_cols, drop = FALSE])
    S <- cor(mat)
    diag(S) <- 0
    list(matrix = S, nodes = colnames(S), directed = FALSE, cleaned_data = data)
  }
  register_estimator("test_custom_boot", custom_fn,
                     "Test custom for bootstrap", directed = FALSE)
  on.exit(remove_estimator("test_custom_boot"), add = TRUE)

  df <- .make_boot_assoc(n = 50, p = 4)
  net <- build_network(df, method = "test_custom_boot")
  boot <- bootstrap_network(net, iter = 20L, seed = 1)

  expect_s3_class(boot, "net_bootstrap")
  expect_equal(boot$method, "test_custom_boot")
})


# ---- Method aliases ----

test_that("bootstrap resolves method aliases via build_network", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "transition"),
                            iter = 20L, seed = 1)
  expect_equal(boot$method, "relative")

  boot2 <- bootstrap_network(build_network(wide, method = "counts"),
                             iter = 20L, seed = 1)
  expect_equal(boot2$method, "frequency")
})


# ---- Validation ----

test_that("invalid inputs error", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative")
  expect_error(bootstrap_network(net, iter = 1L), "iter")
  expect_error(bootstrap_network(net, ci_level = 0), "ci_level")
  expect_error(bootstrap_network(net, ci_level = 1), "ci_level")
  # method = 123 is now a build_network concern, not bootstrap_network
  expect_error(build_network(wide, method = 123))
  expect_error(bootstrap_network(net, inference = "bad"),
               "'arg' should be one of")
})


# ---- Scaling and threshold pass-through ----

test_that("scaling is applied to bootstrap replicates", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative",
                                          scaling = "max"),
                            iter = 20L, seed = 1)
  # Original should have max scaling applied
  expect_true(max(abs(boot$original$weights)) <= 1 + 1e-10)
  # Mean bootstrap values should also be bounded
  expect_true(max(abs(boot$mean)) <= 1 + 1e-10)
})


# ---- netobject_group dispatch (L84-90) ----

test_that("bootstrap_network accepts cograph_network input", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative")
  # Strip the netobject class to get a bare cograph_network
  cograph_net <- net
  class(cograph_net) <- "cograph_network"
  expect_no_error({
    boot <- bootstrap_network(cograph_net, iter = 20L, seed = 1)
  })
  expect_s3_class(boot, "net_bootstrap")
})


# ---- missing $data error (L98-100) ----

test_that("bootstrap_network errors when $data is NULL", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative")
  net$data <- NULL
  expect_error(bootstrap_network(net, iter = 20L),
               "does not contain \\$data")
})


# ---- auto edge_threshold with all-zero weights (L163) ----

test_that("auto edge_threshold defaults to 0 when all weights are zero", {
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative", threshold = 0.99)
  # Force all weights to zero so the auto-threshold logic hits the else branch
  net$weights[] <- 0
  boot <- bootstrap_network(net, iter = 20L, seed = 1, inference = "threshold")
  expect_equal(boot$edge_threshold, 0)
})


# ---- long-format data error in .precompute_per_sequence (L291-295) ----

test_that(".precompute_per_sequence errors on long-format data", {
  # Directly construct a netobject whose stored data still contains an
  # action column so the auto-detect resolves to "long" at bootstrap time.
  states <- c("A", "B", "C")
  set.seed(42)
  wide <- .make_boot_wide()
  net <- build_network(wide, method = "relative")
  # Inject long-format data and params so the fast path triggers the error
  net$data <- data.frame(
    id = rep(1:10, each = 5),
    Action = sample(states, 50, replace = TRUE),
    stringsAsFactors = FALSE
  )
  net$params$format <- "auto"
  net$params$action <- "Action"
  expect_error(
    bootstrap_network(net, iter = 10L),
    "wide-format"
  )
})


# ---- association bootstrap with threshold (L412-414) ----

test_that("bootstrap_network applies threshold in association path", {
  df <- .make_boot_assoc(n = 80, p = 4)
  net <- build_network(df, method = "cor", threshold = 0.01)
  boot <- bootstrap_network(net, iter = 20L, seed = 1)
  expect_s3_class(boot, "net_bootstrap")
  # With a small threshold, significant edges may differ from no-threshold
  expect_true(is.matrix(boot$mean))
})


# ---- print: threshold inference branch (L556) ----

test_that("print shows edge_threshold for threshold inference", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 30L, seed = 1, inference = "threshold")
  out <- capture.output(print(boot))
  expect_true(any(grepl("Threshold", out)))
})


# ---- print: unknown method label (L541 fallback) ----

test_that("print shows generic label for unknown method", {
  wide <- .make_boot_wide()
  boot <- bootstrap_network(build_network(wide, method = "relative"),
                            iter = 20L, seed = 1)
  # Hack the method to something not in the label table
  boot$method <- "custom_unknown_xyz"
  out <- capture.output(print(boot))
  expect_true(any(grepl("custom_unknown_xyz", out)))
})


# ---- Cross-validation against tna::bootstrap ----

test_that("bootstrap_network dispatches for wtna_mixed (L83-95)", {
  set.seed(1)
  oh <- data.frame(
    A = c(1,0,1,0,1,0,1,0), B = c(0,1,0,1,0,1,0,1),
    C = c(1,1,0,0,1,1,0,0)
  )
  mixed <- wtna(oh, method = "both")
  boot_mixed <- bootstrap_network(mixed, iter = 10, seed = 1)
  expect_true(inherits(boot_mixed, "wtna_boot_mixed"))
  expect_true(all(c("transition", "cooccurrence") %in% names(boot_mixed)))
  expect_s3_class(boot_mixed$transition, "net_bootstrap")
  expect_s3_class(boot_mixed$cooccurrence, "net_bootstrap")
})

test_that("print.net_bootstrap_group shows grouped output with shared edges (L660-715)", {
  set.seed(42)
  states <- c("A","B","C")
  n <- 300
  seqs <- data.frame(
    V1 = sample(states, n, TRUE), V2 = sample(states, n, TRUE),
    V3 = sample(states, n, TRUE), V4 = sample(states, n, TRUE),
    grp = rep(c("X","Y"), each = n / 2),
    stringsAsFactors = FALSE
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  boot_grp <- bootstrap_network(nets, iter = 100, seed = 1)
  expect_true(inherits(boot_grp, "net_bootstrap_group"))
  out <- capture.output(print(boot_grp))
  expect_true(any(grepl("Grouped Bootstrap", out)))
  # With enough data, shared significant edges should appear in the table
  expect_true(any(grepl("Edge", out)))
})

test_that("summary.net_bootstrap_group returns combined data frame (L737-741)", {
  set.seed(42)
  states <- c("A","B","C")
  n <- 300
  seqs <- data.frame(
    V1 = sample(states, n, TRUE), V2 = sample(states, n, TRUE),
    V3 = sample(states, n, TRUE), V4 = sample(states, n, TRUE),
    grp = rep(c("X","Y"), each = n / 2),
    stringsAsFactors = FALSE
  )
  nets <- build_network(seqs, method = "relative", group = "grp")
  boot_grp <- bootstrap_network(nets, iter = 100, seed = 1)
  s <- summary(boot_grp)
  expect_true(is.data.frame(s))
  expect_true("group" %in% names(s))
  expect_true(all(c("X", "Y") %in% s$group))
})

test_that("print.wtna_boot_mixed shows both components (L767-772)", {
  set.seed(1)
  oh <- data.frame(
    A = c(1,0,1,0,1,0,1,0), B = c(0,1,0,1,0,1,0,1),
    C = c(1,1,0,0,1,1,0,0)
  )
  mixed <- wtna(oh, method = "both")
  boot_mixed <- bootstrap_network(mixed, iter = 10, seed = 1)
  out <- capture.output(print(boot_mixed))
  expect_true(any(grepl("Mixed Window TNA", out)))
  expect_true(any(grepl("Transition", out)))
  expect_true(any(grepl("Co-occurrence", out)))
})

test_that("summary.wtna_boot_mixed returns list of summaries (L798-801)", {
  set.seed(1)
  oh <- data.frame(
    A = c(1,0,1,0,1,0,1,0), B = c(0,1,0,1,0,1,0,1),
    C = c(1,1,0,0,1,1,0,0)
  )
  mixed <- wtna(oh, method = "both")
  boot_mixed <- bootstrap_network(mixed, iter = 10, seed = 1)
  s <- summary(boot_mixed)
  expect_true(is.list(s))
  expect_true(all(c("transition", "cooccurrence") %in% names(s)))
  expect_true(is.data.frame(s$transition))
  expect_true(is.data.frame(s$cooccurrence))
})

test_that("print.net_bootstrap shows '...and N more' for >5 sig edges (L597)", {
  set.seed(42)
  states <- c("A","B","C")
  wide <- data.frame(
    V1 = sample(states, 500, TRUE), V2 = sample(states, 500, TRUE),
    V3 = sample(states, 500, TRUE), V4 = sample(states, 500, TRUE),
    V5 = sample(states, 500, TRUE), V6 = sample(states, 500, TRUE),
    stringsAsFactors = FALSE
  )
  net <- build_network(wide, method = "relative")
  boot <- bootstrap_network(net, iter = 100, seed = 1)
  out <- capture.output(print(boot))
  expect_true(any(grepl("and.*more significant", out)))
})


# ---- mcml dispatch ----

test_that("bootstrap_network works with mcml objects", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:6], 30, TRUE),
    T2 = sample(LETTERS[1:6], 30, TRUE),
    T3 = sample(LETTERS[1:6], 30, TRUE),
    T4 = sample(LETTERS[1:6], 30, TRUE),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  cs <- build_mcml(seqs, clusters, type = "tna")

  boot <- bootstrap_network(cs, iter = 20, seed = 1)

  expect_s3_class(boot, "net_bootstrap_group")
  expect_true("macro" %in% names(boot))
  expect_s3_class(boot$macro, "net_bootstrap")
  expect_equal(boot$macro$iter, 20L)

  # Within-cluster bootstraps present
  cluster_boots <- setdiff(names(boot), "macro")
  expect_true(length(cluster_boots) > 0)
  for (nm in cluster_boots) {
    expect_s3_class(boot[[nm]], "net_bootstrap")
  }
})

test_that("bootstrap_network mcml dispatch matches explicit as_tna group dispatch", {
  set.seed(42)
  seqs <- data.frame(
    T1 = sample(LETTERS[1:6], 30, TRUE),
    T2 = sample(LETTERS[1:6], 30, TRUE),
    T3 = sample(LETTERS[1:6], 30, TRUE),
    T4 = sample(LETTERS[1:6], 30, TRUE),
    stringsAsFactors = FALSE
  )
  clusters <- list(G1 = c("A", "B", "C"), G2 = c("D", "E", "F"))
  cs <- build_mcml(seqs, clusters, type = "tna")

  boot_mc <- bootstrap_network(cs, iter = 10, seed = 1)
  boot_grp <- bootstrap_network(as_tna(cs), iter = 10, seed = 1)

  expect_identical(class(boot_mc), class(boot_grp))
  expect_identical(names(boot_mc), names(boot_grp))
  expect_identical(
    vapply(boot_mc, function(x) class(x)[1L], character(1L)),
    vapply(boot_grp, function(x) class(x)[1L], character(1L))
  )
  expect_identical(boot_mc$macro$iter, boot_grp$macro$iter)
  expect_equal(boot_mc$macro$summary, boot_grp$macro$summary,
               tolerance = 1e-12)
})

# ---- Branch-matrix coverage (task #17) ----
# Crosses inference x input-class x edge_threshold. Each cell asserts the
# returned object has the expected class and no NA/NaN leaks into the edge
# summary frame. Iter is kept tiny (15) so the full matrix runs in seconds.

test_that("bootstrap_network branch matrix: all combinations produce valid output", {
  set.seed(17)
  wide <- .make_boot_wide(n = 40, t = 8)
  net_tran <- build_network(wide, method = "relative")
  net_freq <- build_network(wide, method = "frequency")
  net_grp  <- {
    wg <- wide
    wg$grp <- rep(c("x", "y"), each = 20)
    build_network(wg, method = "relative", group = "grp")
  }

  grid <- expand.grid(
    inference      = c("stability", "threshold"),
    input          = c("transition", "frequency", "group"),
    edge_threshold = c(NA_real_, 0.05),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    cfg <- grid[i, ]
    obj <- switch(cfg$input,
                  transition = net_tran,
                  frequency  = net_freq,
                  group      = net_grp)
    et  <- if (is.na(cfg$edge_threshold)) NULL else cfg$edge_threshold
    info <- sprintf("inference=%s input=%s edge_threshold=%s",
                    cfg$inference, cfg$input,
                    if (is.null(et)) "NULL" else sprintf("%.2f", et))

    boot <- bootstrap_network(
      obj, iter = 15L, inference = cfg$inference,
      edge_threshold = et, seed = 17L
    )

    if (cfg$input == "group") {
      expect_true(inherits(boot, "net_bootstrap_group"), info = info)
      for (sub in boot) {
        expect_true(inherits(sub, "net_bootstrap"), info = info)
      }
    } else {
      expect_true(inherits(boot, "net_bootstrap"), info = info)
      # Structural guarantees on summary frame + stat matrices
      expect_true(is.data.frame(boot$summary), info = info)
      expect_true(is.matrix(boot$mean),       info = info)
      expect_false(any(is.nan(boot$mean)),    info = paste0(info, " [mean NaN]"))
      expect_false(any(is.nan(boot$sd)),      info = paste0(info, " [sd NaN]"))
      # CI bounds must be ordered lower <= upper elementwise
      expect_true(all(boot$ci_lower <= boot$ci_upper),
                  info = paste0(info, " [ci order]"))
    }
  }
})

test_that("bootstrap_network reproducibility: same seed -> identical stats", {
  set.seed(1)
  net  <- build_network(.make_boot_wide(n = 30, t = 6), method = "relative")
  b1   <- bootstrap_network(net, iter = 20L, inference = "stability", seed = 42L)
  b2   <- bootstrap_network(net, iter = 20L, inference = "stability", seed = 42L)
  # Summary frames + mean matrices must match bit-for-bit under identical seeds
  expect_equal(b1$summary, b2$summary)
  expect_equal(b1$mean,    b2$mean)
  expect_equal(b1$sd,      b2$sd)
})
