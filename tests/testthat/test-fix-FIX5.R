# Regression tests for FIX5 (audit_folder findings A03-F01..A03-F05)
#
# A03-F01  build_hypa() had no help page; build_hypa's roxygen was
#          mis-attached to the exported internal .hypa_one_order.
# A03-F02  build_hypa(order = c(hi, lo)) cograph slots followed argument
#          order instead of the documented "lowest order".
# A03-F03  summary.net_hon / .net_mogen / .net_honem @return wrongly said
#          "input object, invisibly" (they return a data.frame).
# A03-F04  build_hon @return called $nodes a character vector (data.frame).
# A03-F05  fractional order/max_order/k silently floored in hypa / mogen /
#          path_dependence.
#
# Fixtures: group_regulation_long (Actor / Action), ai_long (session/code).

# --- shared fixtures -------------------------------------------------------

.fix5_grl_seqs <- local({
  s <- split(group_regulation_long$Action, group_regulation_long$Actor)
  s <- s[lengths(s) >= 6 & lengths(s) <= 15]
  s[seq_len(min(80L, length(s)))]
})

.fix5_ai_seqs <- local({
  s <- split(ai_long$code, ai_long$session_id)
  s <- s[lengths(s) >= 5 & lengths(s) <= 20]
  s[seq_len(min(40L, length(s)))]
})

# --- A03-F01: build_hypa documentation contract ----------------------------

test_that("A03-F01: build_hypa has its own complete, correct roxygen", {
  hypa_src <- readLines(test_path("..", "..", "R", "hypa.R"))

  # Locate the build_hypa definition and the roxygen block directly above it.
  def_line <- grep("^build_hypa <- function\\(", hypa_src)
  expect_length(def_line, 1L)

  # Walk upward over the contiguous roxygen (#') lines.
  i <- def_line - 1L
  while (i >= 1L && grepl("^#'", hypa_src[i])) i <- i - 1L
  rox <- hypa_src[(i + 1L):(def_line - 1L)]
  expect_gt(length(rox), 20L)               # a real, full block (not bare @export)

  block <- paste(rox, collapse = "\n")
  expect_match(block, "@export", fixed = TRUE)
  expect_match(block, "Detect Path Anomalies via HYPA", fixed = TRUE)

  # @param must document build_hypa's real formals, not data/k of the helper.
  real_formals <- setdiff(names(formals(build_hypa)), "...")
  expect_setequal(real_formals,
                  c("data", "order", "alpha", "min_count", "p_adjust", "k"))
  for (p in real_formals) {
    expect_match(block, paste0("@param ", p, "\\b"),
                 info = paste("missing @param", p))
  }
  # @return must promise the real S3 class.
  expect_match(block, "@return", fixed = TRUE)
  expect_match(block, "net_hypa", fixed = TRUE)
})

test_that("A03-F01: .hypa_one_order roxygen is internal and accurate", {
  hypa_src <- readLines(test_path("..", "..", "R", "hypa.R"))
  def_line <- grep("^\\.hypa_one_order <- function\\(", hypa_src)
  expect_length(def_line, 1L)

  i <- def_line - 1L
  while (i >= 1L && grepl("^#'", hypa_src[i])) i <- i - 1L
  rox <- paste(hypa_src[(i + 1L):(def_line - 1L)], collapse = "\n")

  # Internal: @noRd, NOT @export.
  expect_match(rox, "@noRd", fixed = TRUE)
  expect_false(grepl("@export", rox, fixed = TRUE))

  # @param documents its real formals (trajectories, ord, ...), not data/k.
  expect_setequal(names(formals(Nestimate:::.hypa_one_order)),
                  c("trajectories", "ord", "alpha", "min_count", "p_adjust"))
  expect_match(rox, "@param trajectories\\b")
  expect_match(rox, "@param ord\\b")
  expect_false(grepl("@param data\\b", rox))
  expect_false(grepl("@param k\\b", rox))

  # @return must say it is a bare list (not an S3 net_hypa object).
  expect_match(rox, "list", fixed = TRUE)
})

test_that("A03-F01: post-document Rd contract holds (skipped pre-document)", {
  # The orchestrator runs document() once after all fixers. If man/ has
  # already been regenerated, assert the real Rd files; otherwise skip.
  rd_dir <- test_path("..", "..", "man")
  build_hypa_rd <- file.path(rd_dir, "build_hypa.Rd")
  skip_if_not(file.exists(build_hypa_rd),
              "man/build_hypa.Rd not yet regenerated (orchestrator step)")

  # Check the regenerated source man/ directly. tools::Rd_db() needs an
  # *installed* package and errors under devtools::load_all(); the source
  # Rd files are the reliable, load_all-safe contract surface.
  expect_true(file.exists(build_hypa_rd))
  # .hypa_one_order is now @noRd: no generated (exported) Rd page.
  expect_false(file.exists(file.path(rd_dir, "dot-hypa_one_order.Rd")))
  rd <- paste(readLines(build_hypa_rd, warn = FALSE), collapse = "\n")
  expect_match(rd, "\\\\title\\{")
  expect_match(rd, "\\\\arguments\\{")
  expect_match(rd, "build_hypa")
})

# --- A03-F02: build_hypa(order=) contract -----------------------------------

test_that("A03-F02: cograph slots follow lowest order, not argument order", {
  y2 <- build_hypa(.fix5_grl_seqs, order = 2L, min_count = 2L)
  y3 <- build_hypa(.fix5_grl_seqs, order = 3L, min_count = 2L)
  yhi_lo <- build_hypa(.fix5_grl_seqs, order = c(3L, 2L), min_count = 2L)
  ylo_hi <- build_hypa(.fix5_grl_seqs, order = c(2L, 3L), min_count = 2L)

  # Primary cograph slot == the LOWEST requested order (order 2), regardless
  # of how the order vector is arranged.
  expect_identical(yhi_lo$weights, y2$weights)
  expect_false(identical(yhi_lo$weights, y3$weights))

  # Arg-order independent: c(3,2) and c(2,3) give identical primary slots.
  expect_identical(yhi_lo$weights, ylo_hi$weights)
  expect_identical(yhi_lo$adjacency, ylo_hi$adjacency)
  expect_identical(yhi_lo$xi, ylo_hi$xi)
  expect_identical(yhi_lo$edges, ylo_hi$edges)
  expect_identical(yhi_lo$nodes, ylo_hi$nodes)

  # $order is stored sorted ascending; $scores still spans every order.
  expect_identical(yhi_lo$order, c(2L, 3L))
  expect_setequal(unique(yhi_lo$scores$order), c(2L, 3L))
})

# --- A03-F03: summary.* @return is a data.frame ----------------------------

test_that("A03-F03: summary.net_hon/.net_mogen/.net_honem return data.frames", {
  hon  <- build_hon(.fix5_ai_seqs, max_order = 2L)
  mog  <- build_mogen(.fix5_ai_seqs, max_order = 2L)
  hem  <- build_honem(hon, dim = 2L)

  sh <- summary(hon)
  expect_s3_class(sh, "data.frame")
  expect_identical(sh, hon$edges)
  expect_false(identical(sh, hon))

  sm <- summary(mog)
  expect_s3_class(sm, "data.frame")
  expect_true(all(c("order", "aic", "bic", "selected") %in% names(sm)))
  expect_false(identical(sm, mog))

  se <- summary(hem)
  expect_s3_class(se, "data.frame")
  expect_true("node" %in% names(se))
  expect_false(identical(se, hem))

  # @return text must no longer claim "input object, invisibly".
  for (f in c("hon.R", "mogen.R", "honem.R")) {
    src <- paste(readLines(test_path("..", "..", "R", f)), collapse = "\n")
    blocks <- regmatches(
      src,
      gregexpr("#' Summary Method for net_[a-z]+.*?@return[^\n]*", src))[[1]]
    for (b in blocks) {
      expect_false(grepl("The input object, invisibly", b, fixed = TRUE),
                   info = paste("stale @return still in", f))
    }
  }
})

# --- A03-F04: build_hon @return $nodes is a data.frame ----------------------

test_that("A03-F04: build_hon $nodes is a data.frame and @return says so", {
  hon <- build_hon(.fix5_ai_seqs, max_order = 2L)
  expect_s3_class(hon$nodes, "data.frame")
  expect_setequal(names(hon$nodes), c("id", "label", "name"))

  hon_src <- paste(readLines(test_path("..", "..", "R", "hon.R")),
                   collapse = "\n")
  expect_false(grepl("Character vector of HON node names", hon_src,
                      fixed = TRUE))
  expect_match(hon_src,
               "\\\\item\\{nodes\\}\\{data\\.frame with columns")
})

# --- A03-F05: non-integer order/max_order/k is a clean error ----------------

test_that("A03-F05: fractional order/max_order/k errors instead of flooring", {
  # hypa
  expect_error(build_hypa(.fix5_grl_seqs, order = 2.7, min_count = 2L),
               "whole number")
  expect_s3_class(build_hypa(.fix5_grl_seqs, order = 2L, min_count = 2L),
                  "net_hypa")
  # double-valued whole numbers are still accepted (2.0 == round(2.0))
  expect_s3_class(build_hypa(.fix5_grl_seqs, order = 2, min_count = 2L),
                  "net_hypa")

  # mogen: build_mogen / mogen_transitions / path_counts
  expect_error(build_mogen(.fix5_grl_seqs, max_order = 2.7), "whole number")
  m <- build_mogen(.fix5_grl_seqs, max_order = 3L)
  expect_s3_class(m, "net_mogen")
  expect_s3_class(build_mogen(.fix5_grl_seqs, max_order = 3), "net_mogen")
  expect_error(mogen_transitions(m, order = 1.5), "whole number")
  expect_true(is.data.frame(mogen_transitions(m, order = 1L)))
  expect_true(is.data.frame(mogen_transitions(m)))   # NULL -> optimal, ok
  expect_error(path_counts(.fix5_grl_seqs, k = 2.7), "whole number")
  expect_true(is.data.frame(path_counts(.fix5_grl_seqs, k = 2L)))
  expect_true(is.data.frame(path_counts(.fix5_grl_seqs, k = 3)))

  # path_dependence
  expect_error(
    path_dependence(as.data.frame(trajectories), order = 2.7,
                    min_count = 2L),
    "whole number")
  expect_s3_class(
    path_dependence(as.data.frame(trajectories), order = 2L,
                    min_count = 2L),
    "net_path_dependence")
})
