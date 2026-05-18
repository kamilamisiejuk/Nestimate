# =========================================================================
# Regression tests for FIX-1 (audit A01 core estimation + A07 mlVAR/MGM)
#
# Confirmed findings fixed:
#   A01-F01 (HIGH)  scaling="minmax" mapped the smallest genuine edge of an
#                   association net (cor/pcor/glasso/ising) to exactly 0,
#                   dropping a real edge and corrupting $n_edges. Fix:
#                   structure-preserving minmax maps the non-zero block into
#                   (eps, 1] so the block minimum becomes eps, not 0.
#   A01-F02 (MED)   level= on a *directed* estimator was a numeric no-op yet
#                   printed "[between-person]" and returned a fake
#                   netobject_ml. Fix: stop() with a clear message.
#   A01-F03 (LOW)   @param mode doc-only: now states it applies to one-hot
#                   input only (verified: still inert for wide, no error).
#   A01-F04 (LOW)   @param method roxygen completeness — verified the 5
#                   omitted estimators + 7 omitted aliases all resolve.
#   A07-F01 (MED)   build_network(method="mgm", level=<vec>) collided with
#                   the multilevel `level` formal. Fix: route a numeric
#                   length-ncol level into params for mgm; informative error
#                   for a wrong-length numeric level.
#   A07-F02 (HIGH)  build_network(method="mgm", threshold="LW"/"none")
#                   collided with the numeric `threshold` formal (dead enum).
#                   Fix: route the mgm threshold enum into params; clear
#                   error for an invalid string; numeric cutoff still works.
#   A07-F03 (MED)   .estimator_mgm(level=) is numerically inert — roxygen now
#                   states it is a length-validated parity-only declaration.
#                   Test asserts the documented "no effect" is true.
#   A07-F04 (MED)   build_mlvar roxygen no longer promises netobject_group
#                   bootstrap iteration. Test asserts bootstrap_network()
#                   still errors (code unchanged — doc was the fix).
#   A07-F06 (LOW)   .estimator_mgm auto-detect crashed on NA via un-na.rm'd
#                   all(col == round(col)). Fix: na.rm + drop NA from unique.
#
# Data: simulate_continuous() (strong rho) for association nets,
# simulate_data("mlvar") for mlVAR, as.factor() mixed data + bundled
# group_regulation_long where a directed sequence network is needed.
# =========================================================================

# --------------------------------------------------------------------- #
# A01-F01 — minmax must not delete the weakest association edge          #
# --------------------------------------------------------------------- #

test_that("A01-F01: scaling='minmax' preserves $n_edges for association nets", {
  skip_if_not_installed("glasso")
  cont <- simulate_continuous(n = 120, p = 6, rho = 0.45, seed = 42)

  for (m in c("cor", "pcor", "glasso")) {
    base   <- build_network(cont, method = m)
    scaled <- build_network(cont, method = m, scaling = "minmax")
    expect_equal(scaled$n_edges, base$n_edges,
                 info = paste("method", m, "edge count must survive minmax"))
    expect_equal(nrow(scaled$edges), nrow(base$edges), info = m)
    # weights stay in [0, 1] and the block minimum is strictly positive
    nz <- scaled$weights[scaled$weights != 0]
    expect_true(all(scaled$weights >= 0 & scaled$weights <= 1), info = m)
    expect_true(min(nz) > 0, info = m)
  }
})

test_that("A01-F01: minmax leaves the transition (include_zeros) path unchanged", {
  # frequency uses include_zeros = TRUE — structural zeros (A->A absent) must
  # still map to 0 and the block max to 1 (this path is NOT touched by the fix).
  seqs <- data.frame(V1 = c("A", "A", "B"),
                     V2 = c("B", "C", "C"),
                     V3 = c("C", "C", NA))
  scaled <- build_network(seqs, method = "frequency", scaling = "minmax")
  expect_equal(scaled$weights["A", "A"], 0)
  expect_equal(scaled$weights["B", "C"], 1)
})

test_that("A01-F01: .apply_scaling minmax edge cases unchanged", {
  # all non-zero equal -> unchanged; all-zero -> unchanged (include_zeros=FALSE)
  expect_equal(.apply_scaling(matrix(c(0, .5, .5, 0), 2), "minmax"),
               matrix(c(0, .5, .5, 0), 2))
  expect_equal(.apply_scaling(matrix(0, 3, 3), "minmax"), matrix(0, 3, 3))
})

# --------------------------------------------------------------------- #
# A01-F02 — level= on a directed estimator must error, not mislabel      #
# --------------------------------------------------------------------- #

test_that("A01-F02: level= errors on directed estimators (no fake netobject_ml)", {
  set.seed(1)
  seqs <- data.frame(id = rep(1:10, each = 4),
                     V1 = sample(LETTERS[1:4], 40, TRUE),
                     V2 = sample(LETTERS[1:4], 40, TRUE),
                     V3 = sample(LETTERS[1:4], 40, TRUE))
  for (m in c("relative", "frequency", "attention")) {
    for (lv in c("between", "within", "both")) {
      expect_error(
        build_network(seqs, method = m, params = list(id = "id"), level = lv),
        "only supported for undirected",
        info = paste(m, lv)
      )
    }
  }
})

test_that("A01-F02: level= still works for undirected association methods", {
  set.seed(2)
  df <- data.frame(id = rep(1:20, each = 5),
                   x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100))
  b  <- build_network(df, method = "cor", params = list(id = "id"),
                       level = "between")
  expect_s3_class(b, "netobject")
  expect_identical(b$level, "between")
  bb <- build_network(df, method = "cor", params = list(id = "id"),
                       level = "both")
  expect_s3_class(bb, "netobject_ml")
})

# --------------------------------------------------------------------- #
# A01-F03 — mode= is inert for wide data (documented as one-hot-only)    #
# --------------------------------------------------------------------- #

test_that("A01-F03: mode= is a documented no-op on wide sequence data", {
  set.seed(1)
  seqs <- data.frame(V1 = sample(LETTERS[1:4], 60, TRUE),
                     V2 = sample(LETTERS[1:4], 60, TRUE),
                     V3 = sample(LETTERS[1:4], 60, TRUE))
  m1 <- build_network(seqs, method = "relative", mode = "non-overlapping")
  m2 <- build_network(seqs, method = "relative", mode = "overlapping")
  expect_identical(m1$weights, m2$weights)  # inert for wide, as documented
  # invalid mode still rejected by match.arg
  expect_error(build_network(seqs, method = "relative", mode = "sliding"))
})

# --------------------------------------------------------------------- #
# A01-F04 — every documented method/alias resolves and builds            #
# --------------------------------------------------------------------- #

test_that("A01-F04: the 5 omitted estimators + 7 omitted aliases all resolve", {
  reg <- list_estimators()$name
  expect_true(all(c("ising", "mgm", "attention", "wtna",
                     "wtna_cooccurrence") %in% reg))
  alias_map <- c(isingfit = "ising", atna = "attention",
                 mixed = "mgm", mixed_graphical = "mgm",
                 wtna_transition = "wtna", wcna = "co_occurrence",
                 # documented already (must remain so)
                 counts = "frequency", partial = "pcor")
  for (a in names(alias_map)) {
    expect_identical(Nestimate:::.resolve_method_alias(a),
                     unname(alias_map[a]), info = a)
  }
})

# --------------------------------------------------------------------- #
# A07-F01 / A07-F02 — mgm level/threshold reachable via build_network    #
# --------------------------------------------------------------------- #

test_that("A07-F02: build_network(method='mgm', threshold='LW'/'none') reachable", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 250
  dd <- data.frame(V1 = rnorm(n))
  dd$V2 <- 0.7 * dd$V1 + rnorm(n, sd = .5)
  dd$V3 <- rnorm(n)

  none_direct <- build_network(dd, method = "mgm", threshold = "none")
  none_params <- build_network(dd, method = "mgm",
                               params = list(threshold = "none"))
  expect_s3_class(none_direct, "netobject")
  # the direct route must be byte-identical to the documented params route
  expect_equal(none_direct$weights, none_params$weights)

  lw_direct <- build_network(dd, method = "mgm", threshold = "LW")
  lw_params <- build_network(dd, method = "mgm",
                             params = list(threshold = "LW"))
  expect_equal(lw_direct$weights, lw_params$weights)

  # both documented enum values are now genuinely reachable as netobjects
  # via the public entry point (the dead-enum is fixed). No numeric claim
  # about LW-vs-none differences is asserted — mgm's LW tau frequently does
  # not bite once glmnet's own L1 shrinkage has zeroed weak coefficients.
  expect_true(is.matrix(lw_direct$weights))
  expect_true(is.matrix(none_direct$weights))

  # invalid string errors clearly; numeric cutoff still works for mgm
  expect_error(build_network(dd, method = "mgm", threshold = "bogus"),
               "must be \"LW\" or \"none\"")
  expect_s3_class(build_network(dd, method = "mgm", threshold = 0.05),
                  "netobject")
})

test_that("A07-F01: build_network(method='mgm', level=<vec>) reachable + clear error", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 250
  dd <- data.frame(V1 = rnorm(n))
  dd$V2 <- 0.7 * dd$V1 + rnorm(n, sd = .5)
  dd$V3 <- rnorm(n)

  lvl_direct <- build_network(dd, method = "mgm", level = c(1L, 1L, 1L))
  lvl_params <- build_network(dd, method = "mgm",
                              params = list(level = c(1L, 1L, 1L)))
  expect_s3_class(lvl_direct, "netobject")
  expect_equal(lvl_direct$weights, lvl_params$weights)

  # wrong-length numeric level -> informative error (mentions mgm + the
  # build_network decomposition enum), not "'arg' must be NULL or a
  # character vector"
  err <- tryCatch(build_network(dd, method = "mgm", level = 2),
                  error = function(e) conditionMessage(e))
  expect_true(grepl("mgm", err) && grepl("integer vector", err))
  expect_false(grepl("must be NULL or a character vector", err))

  # non-mgm methods keep the original enum semantics unchanged
  expect_error(build_network(dd, method = "cor", level = 2))
  expect_error(build_network(dd, method = "cor", threshold = "LW"),
               "is.numeric")
})

# --------------------------------------------------------------------- #
# A07-F03 — mgm level is numerically inert (documented as parity-only)   #
# --------------------------------------------------------------------- #

test_that("A07-F03: .estimator_mgm level= has no numerical effect (as documented)", {
  skip_if_not_installed("glmnet")
  cd <- simulate_continuous(n = 200, p = 5, rho = 0.3, seed = 1)
  L1 <- .estimator_mgm(cd, type = rep("g", 5), level = rep(1L, 5))
  L9 <- .estimator_mgm(cd, type = rep("g", 5), level = rep(99L, 5))
  expect_equal(L1$matrix, L9$matrix)

  # also inert through the reachable build_network params route
  b1 <- build_network(cd, method = "mgm",
                       params = list(type = rep("g", 5), level = rep(1L, 5)))
  b9 <- build_network(cd, method = "mgm",
                       params = list(type = rep("g", 5), level = rep(99L, 5)))
  expect_equal(b1$weights, b9$weights)
})

# --------------------------------------------------------------------- #
# A07-F04 — build_mlvar constituents carry no $data (doc was the fix)    #
# --------------------------------------------------------------------- #

test_that("A07-F04: net_mlvar constituents are data-less; bootstrap still errors", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("corpcor")
  skip_if_not_installed("data.table")
  d   <- simulate_data("mlvar", seed = 3, n_subjects = 12, d = 4, n_obs = 30)
  fit <- suppressWarnings(
    build_mlvar(d, vars = attr(d, "vars"), id = "id", day = "day",
                beep = "beep")
  )
  # the roxygen now correctly states these are data-less wrappers
  expect_null(fit$temporal$data)
  expect_null(fit$contemporaneous$data)
  expect_null(fit$between$data)
  # verbs the doc DOES promise still work
  expect_s3_class(fit, "net_mlvar")
  expect_s3_class(coefs(fit), "data.frame")
  expect_invisible(print(fit))
  # bootstrap_network is NOT promised and still errors (code unchanged)
  expect_error(bootstrap_network(fit, iter = 5),
               "does not contain \\$data")
})

# --------------------------------------------------------------------- #
# A07-F06 — mgm auto-detect must not crash on a numeric NA column        #
# --------------------------------------------------------------------- #

test_that("A07-F06: .estimator_mgm auto-detect no longer crashes on NA && short-circuit", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 300
  k <- sample(1:4, n, TRUE); k[7] <- NA
  dat0 <- data.frame(V1 = rnorm(n), K = as.numeric(k))
  msg <- tryCatch(.estimator_mgm(dat0), error = function(e) conditionMessage(e))
  # the specific cryptic && NA crash must be gone
  expect_false(grepl("missing value where TRUE/FALSE needed", msg))

  # clean data still classifies correctly (no regression in auto-detect)
  r_cont <- .estimator_mgm(
    data.frame(V1 = rnorm(200), V2 = rnorm(200), V3 = rnorm(200))
  )
  expect_identical(unname(r_cont$type), rep("g", 3))
  r_mix <- .estimator_mgm(
    data.frame(V1 = rnorm(200), V2 = rnorm(200),
               G = sample(1:3, 200, TRUE))
  )
  expect_identical(unname(r_mix$type), c("g", "g", "c"))

  # bundled directed sequence network sanity (group_regulation_long) — the
  # A01 changes must not perturb the standard relative-network path.
  e <- new.env()
  utils::data("group_regulation_long", package = "Nestimate", envir = e)
  net <- build_network(e$group_regulation_long, actor = "Actor",
                       time = "Time", action = "Action", method = "relative")
  expect_s3_class(net, "netobject")
  expect_true(net$directed)
  # directed + level still errors cleanly here too (A01-F02 path)
  expect_error(
    build_network(e$group_regulation_long, actor = "Actor", time = "Time",
                  action = "Action", method = "relative", level = "between"),
    "only supported for undirected"
  )
})
