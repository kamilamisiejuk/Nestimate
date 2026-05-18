# Regression: Codex P2 on the fix diff -- build_network(method="mgm",
# level=, group=) validated the mgm `level` length against raw ncol(data),
# before the group column is dropped for the per-group recursion. A correct
# level was rejected; a group-padded one was forwarded too long to
# .estimator_mgm(). Fix: validate against the modeled-variable count
# (ncol minus the grouping column(s)).

test_that("mgm + group: a correct level (one per modeled var) is accepted", {
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 160
  df <- data.frame(g1 = rnorm(n), g2 = rnorm(n),
                   c1 = factor(sample(letters[1:3], n, TRUE)),
                   grp = sample(c("X", "Y"), n, TRUE),
                   stringsAsFactors = FALSE)
  # modeled vars = g1,g2,c1 (3); grp is the grouping column and is dropped.
  net <- build_network(df, method = "mgm", level = c(1L, 1L, 3L),
                       group = "grp")
  expect_s3_class(net, "netobject_group")
  expect_setequal(names(net), c("X", "Y"))
})

test_that("mgm + group: wrong-length level errors, naming the modeled count", {
  skip_if_not_installed("glmnet")
  set.seed(2); n <- 120
  df <- data.frame(g1 = rnorm(n), g2 = rnorm(n),
                   c1 = factor(sample(letters[1:3], n, TRUE)),
                   grp = sample(c("X", "Y"), n, TRUE),
                   stringsAsFactors = FALSE)
  # padding to raw ncol (4) is now correctly rejected.
  expect_error(
    build_network(df, method = "mgm", level = c(1L, 1L, 1L, 1L),
                  group = "grp"),
    "length 3"
  )
})

test_that("mgm without group: level length == ncol still works (no regression)", {
  skip_if_not_installed("glmnet")
  set.seed(3); n <- 140
  df <- data.frame(g1 = rnorm(n), g2 = rnorm(n),
                   c1 = factor(sample(letters[1:3], n, TRUE)),
                   stringsAsFactors = FALSE)
  net <- build_network(df, method = "mgm", level = c(1L, 1L, 3L))
  expect_s3_class(net, "netobject")
})
