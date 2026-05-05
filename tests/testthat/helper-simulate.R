# Minimal simulation helpers for Nestimate tests
# These replace Saqrlab simulation functions with standalone versions

simulate_mtna <- function(n_nodes = 5, n_types = 5,
                          type_names = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(type_names)) {
    all_types <- c("Metacognitive", "Cognitive", "Behavioral",
                   "Social", "Motivational", "Affective", "GroupRegulation")
    type_names <- all_types[seq_len(min(n_types, length(all_types)))]
    if (n_types > length(all_types)) {
      type_names <- c(type_names, paste0("Type", (length(all_types) + 1):n_types))
    }
  }
  total <- n_nodes * n_types
  node_names <- character(0)
  node_types_list <- list()
  for (i in seq_len(n_types)) {
    nms <- paste0(type_names[i], "_", seq_len(n_nodes))
    node_names <- c(node_names, nms)
    node_types_list[[type_names[i]]] <- nms
  }
  mat <- matrix(runif(total * total), nrow = total)
  diag(mat) <- 0
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  mat <- mat / row_sums
  rownames(mat) <- colnames(mat) <- node_names
  list(
    matrix = round(mat, 4),
    node_types = node_types_list,
    type_names = type_names,
    n_nodes_per_type = stats::setNames(rep(n_nodes, n_types), type_names)
  )
}

simulate_data <- function(type = "mlvar", seed = NULL, n_subjects = 20,
                          d = 6, n_obs = 40, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (type == "mlvar") {
    vars <- paste0("V", seq_len(d))
    rows <- n_subjects * n_obs

    # Generate known temporal and contemporaneous matrices
    true_B <- matrix(0, d, d)
    diag(true_B) <- runif(d, 0.2, 0.5)
    # Add a few off-diagonal temporal effects
    for (i in seq_len(min(d - 1, 3))) {
      true_B[i, i + 1] <- runif(1, 0.1, 0.3)
    }

    true_Theta <- diag(d)
    for (i in seq_len(min(d - 1, 2))) {
      true_Theta[i, i + 1] <- true_Theta[i + 1, i] <- runif(1, 0.1, 0.3)
    }

    # Generate VAR data
    between_means <- matrix(rnorm(n_subjects * d, 0, 1), n_subjects, d)
    df <- data.frame(
      id = rep(seq_len(n_subjects), each = n_obs),
      day = rep(rep(seq_len(n_obs %/% 10), each = 10), n_subjects),
      beep = rep(rep(seq_len(10), n_obs %/% 10), n_subjects)
    )
    y_mat <- matrix(0, rows, d)
    for (s in seq_len(n_subjects)) {
      idx <- ((s - 1) * n_obs + 1):(s * n_obs)
      y_mat[idx[1], ] <- between_means[s, ] + rnorm(d)
      for (t in 2:n_obs) {
        innovation <- rnorm(d)
        y_mat[idx[t], ] <- between_means[s, ] +
          true_B %*% (y_mat[idx[t - 1], ] - between_means[s, ]) + innovation
      }
    }
    for (j in seq_len(d)) df[[vars[j]]] <- y_mat[, j]

    rownames(true_B) <- colnames(true_B) <- vars
    rownames(true_Theta) <- colnames(true_Theta) <- vars

    attr(df, "vars") <- vars
    attr(df, "type") <- "mlvar"
    attr(df, "true_temporal") <- true_B
    attr(df, "true_contemporaneous") <- true_Theta
    attr(df, "info") <- list()
    return(df)
  }
  stop(sprintf("simulate_data type '%s' not supported in Nestimate test helper", type))
}


# ---- Safe package-load skip ----

#' Skip if a package cannot actually be loaded (catches ABI/linker failures,
#' not just missing installations).
#' @noRd
skip_if_pkg_broken <- function(pkg) {
  ok <- isTRUE(tryCatch(requireNamespace(pkg, quietly = TRUE),
                        error = function(e) FALSE))
  if (!ok) skip(paste0(pkg, " cannot be loaded"))
}

# ---- Equivalence test infrastructure ----

#' Skip equivalence tests unless explicitly enabled
#' @noRd
skip_equiv_tests <- function() {
  run <- Sys.getenv("NESTIMATE_EQUIV_TESTS", unset = "false")
  if (!identical(run, "true")) {
    skip("Equivalence tests skipped (set NESTIMATE_EQUIV_TESTS=true)")
  }
}

#' Generate random sequence data for equivalence testing
#' @noRd
simulate_sequences <- function(n_actors = 10, n_states = 5,
                               seq_length = 20, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  states <- LETTERS[seq_len(n_states)]
  rows <- lapply(seq_len(n_actors), function(i) {
    sample(states, seq_length, replace = TRUE)
  })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(df) <- paste0("T", seq_len(seq_length))
  df
}

#' Real-data sequence anchor for equivalence tests.
#'
#' Returns a list of character sequences extracted from a bundled long-format
#' dataset, ordered within session. Used by higher-order equivalence tests to
#' anchor synthetic-data validation against realistic state imbalance and
#' non-trivial higher-order dependencies that random uniform sampling smooths
#' away.
#'
#' @param dataset One of "human_long", "ai_long", "group_regulation_long".
#' @param max_actors Optional cap on number of sessions returned (NULL = all).
#' @noRd
bundled_sequences <- function(dataset = "human_long", max_actors = NULL) {
  e <- new.env()
  utils::data(list = dataset, package = "Nestimate", envir = e)
  d <- e[[dataset]]
  stopifnot(all(c("session_id", "code", "order_in_session") %in% names(d)))
  d <- d[order(d$session_id, d$order_in_session), ]
  d <- d[!duplicated(d[, c("session_id", "order_in_session")]), ]
  ids <- unique(d$session_id)
  if (!is.null(max_actors)) ids <- ids[seq_len(min(length(ids), max_actors))]
  lapply(ids, function(s) as.character(d$code[d$session_id == s]))
}

#' Generate random continuous data for association method equivalence testing
#' @noRd
simulate_continuous <- function(n = 100, p = 5, rho = 0.3, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Sigma <- diag(p)
  for (i in seq_len(p - 1L)) {
    Sigma[i, i + 1L] <- rho
    Sigma[i + 1L, i] <- rho
  }
  L <- chol(Sigma)
  mat <- matrix(rnorm(n * p), n, p) %*% L
  df <- as.data.frame(mat)
  colnames(df) <- paste0("V", seq_len(p))
  df
}

#' Generate a random binary incidence matrix for hypergraph equivalence tests.
#'
#' Columns are hyperedges, rows are nodes. Each hyperedge's size is drawn
#' uniformly from `edge_size_range`, members are sampled without replacement
#' from the node set. `density` adds per-edge jitter: with probability
#' (1 - density) an edge is skipped (returned as a degenerate size-1 edge, or
#' dropped if too small), ensuring size variability across runs.
#'
#' @param n_nodes Integer >= 3. Number of nodes.
#' @param n_edges Integer >= 1. Number of hyperedges.
#' @param edge_size_range Integer vector length 2. Min/max hyperedge size.
#'   Must satisfy 2 <= min <= max <= n_nodes.
#' @param density Numeric in (0, 1]. Probability that a given edge is kept.
#' @param seed Integer or NULL. Random seed.
#' @return list(incidence = matrix(n_nodes, n_edges_kept),
#'              hyperedges = list of integer node indices,
#'              nodes = character vector of node names)
#' @noRd
simulate_hypergraph_incidence <- function(n_nodes = 15, n_edges = 20,
                                          edge_size_range = c(2L, 4L),
                                          density = 0.8,
                                          seed = NULL) {
  stopifnot(n_nodes >= 3L, n_edges >= 1L,
            length(edge_size_range) == 2L,
            edge_size_range[1] >= 2L,
            edge_size_range[2] >= edge_size_range[1],
            edge_size_range[2] <= n_nodes,
            density > 0, density <= 1)
  if (!is.null(seed)) set.seed(seed)
  nodes <- paste0("n", seq_len(n_nodes))
  sizes <- sample(edge_size_range[1]:edge_size_range[2], n_edges, replace = TRUE)
  keep <- stats::runif(n_edges) < density
  edges <- mapply(function(sz, k) {
    if (!k) return(integer(0))
    sort(sample.int(n_nodes, sz))
  }, sizes, keep, SIMPLIFY = FALSE)
  edges <- edges[vapply(edges, length, integer(1)) >= 2L]
  if (length(edges) == 0L) {
    edges <- list(sort(sample.int(n_nodes, edge_size_range[1])))
  }
  inc <- matrix(0L, nrow = n_nodes, ncol = length(edges),
                dimnames = list(nodes, paste0("e", seq_along(edges))))
  for (j in seq_along(edges)) inc[edges[[j]], j] <- 1L
  list(incidence = inc, hyperedges = edges, nodes = nodes)
}

#' Equivalence report logger
#' @noRd
equiv_report <- function() {
  env <- new.env(parent = emptyenv())
  env$rows <- list()

  find_repo_root <- function(path = getwd()) {
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
      parent <- dirname(path)
      if (identical(parent, path)) return(getwd())
      path <- parent
    }
  }

  env$log <- function(func, config, n_checked, n_failed,
                      max_abs_err, mean_abs_err, median_abs_err,
                      p95_abs_err, reference, notes = "") {
    env$rows[[length(env$rows) + 1L]] <- data.frame(
      func = func, config = config,
      n_checked = n_checked, n_failed = n_failed,
      max_abs_err = max_abs_err, mean_abs_err = mean_abs_err,
      median_abs_err = median_abs_err, p95_abs_err = p95_abs_err,
      reference = reference, notes = notes,
      stringsAsFactors = FALSE
    )
  }

  env$write_csv <- function(module) {
    if (length(env$rows) == 0L) return(invisible(NULL))
    df <- do.call(rbind, env$rows)
    out_dir <- file.path(find_repo_root(), "tmp")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    path <- file.path(out_dir, sprintf("%s_equivalence_report.csv", module))
    write.csv(df, path, row.names = FALSE)
    message(sprintf("Equivalence report: %s (%d checks)", path, sum(df$n_checked)))
  }

  env$write_cvs <- function(module, test_file = NULL) {
    if (length(env$rows) == 0L) return(invisible(NULL))
    df <- do.call(rbind, env$rows)

    # Build vitest-compatible assertion results
    assertions <- lapply(seq_len(nrow(df)), function(i) {
      r <- df[i, ]
      passed <- r$n_failed == 0
      title <- sprintf("%s: %s delta=%.2e", r$func, r$config, r$max_abs_err)
      list(
        ancestorTitles = list(paste0(module, " equivalence")),
        title = title,
        fullName = sprintf("%s equivalence > %s", module, title),
        status = if (passed) "passed" else "failed",
        duration = 0,
        failureMessages = if (passed) list() else list(
          sprintf("max delta %.2e >= tolerance, %d/%d values failed",
                  r$max_abs_err, r$n_failed, r$n_checked)
        ),
        `_cvs` = list(
          delta = r$max_abs_err,
          tolerance = 1e-10,
          rFunction = r$func,
          rPackage = r$reference,
          module = module,
          target = "nestimate"
        )
      )
    })

    n_passed <- sum(df$n_failed == 0)
    n_failed <- sum(df$n_failed > 0)
    result <- list(
      numTotalTestSuites = 1L,
      numPassedTestSuites = if (n_failed == 0) 1L else 0L,
      numFailedTestSuites = if (n_failed > 0) 1L else 0L,
      numTotalTests = nrow(df),
      numPassedTests = n_passed,
      numFailedTests = n_failed,
      testResults = list(list(
        name = if (!is.null(test_file)) test_file
               else sprintf("tests/testthat/test-equiv-%s.R", module),
        assertionResults = assertions
      ))
    )

    inbox_candidates <- c(
      file.path(find_repo_root(), "validation", "data", "inbox"),
      file.path(dirname(find_repo_root()), "validation", "data", "inbox"),
      file.path("..", "..", "..", "validation", "data", "inbox"),
      file.path("..", "..", "validation", "data", "inbox")
    )
    inbox <- inbox_candidates[dir.exists(inbox_candidates) &
                                file.access(inbox_candidates, 2) == 0][1L]
    if (!is.na(inbox)) {
      ts <- format(Sys.time(), "%Y%m%dT%H%M%S")
      path <- file.path(inbox, sprintf("nestimate-%s-%s.json", module, ts))
      writeLines(jsonlite::toJSON(result, auto_unbox = TRUE, pretty = TRUE), path)
      message(sprintf("CVS report: %s", path))
    } else {
      message("CVS report skipped: validation/data/inbox is not writable")
    }
  }

  env$summary <- function() {
    if (length(env$rows) == 0L) return("No results logged.")
    df <- do.call(rbind, env$rows)
    sprintf("Total: %d values checked, %d failed, max delta %.2e",
            sum(df$n_checked), sum(df$n_failed), max(df$max_abs_err))
  }

  env
}
