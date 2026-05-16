# ---- Global Estimator Registry ----

#' @noRd
.estimator_registry <- new.env(parent = emptyenv())


#' Register a Network Estimator
#'
#' @description
#' Register a custom or built-in network estimator function by name.
#' Estimators registered here can be used by \code{\link{estimate_network}}
#' via the \code{method} parameter.
#'
#' @param name Character. Unique name for the estimator (e.g. \code{"relative"},
#'   \code{"glasso"}).
#' @param fn Function. The estimator function. Must accept \code{data} as its
#'   first argument and \code{...} for additional parameters. Must return a list
#'   with at least: \code{matrix} (square numeric matrix), \code{nodes}
#'   (character vector), \code{directed} (logical).
#' @param description Character. Short description of the estimator.
#' @param directed Logical. Whether the estimator produces directed networks.
#'
#' @return Invisible \code{NULL}.
#'
#' @examples
#' my_fn <- function(data, ...) {
#'   m <- cor(data)
#'   diag(m) <- 0
#'   list(matrix = m, nodes = colnames(m), directed = FALSE)
#' }
#' register_estimator("my_cor", my_fn, "Custom correlation", directed = FALSE)
#' df <- data.frame(A = rnorm(20), B = rnorm(20), C = rnorm(20))
#' net <- build_network(df, method = "my_cor")
#' remove_estimator("my_cor")
#'
#' @seealso \code{\link{get_estimator}}, \code{\link{list_estimators}},
#'   \code{\link{remove_estimator}}, \code{\link{estimate_network}}
#'
#' @export
register_estimator <- function(name, fn, description, directed) {
  stopifnot(
    is.character(name), length(name) == 1, nzchar(name),
    is.function(fn),
    is.character(description), length(description) == 1,
    is.logical(directed), length(directed) == 1
  )
  entry <- list(fn = fn, description = description, directed = directed)
  assign(name, entry, envir = .estimator_registry)
  invisible(NULL)
}


#' Retrieve a Registered Estimator
#'
#' @description
#' Retrieve a registered network estimator by name.
#'
#' @param name Character. Name of the estimator to retrieve.
#'
#' @return A list with elements \code{fn}, \code{description}, \code{directed}.
#'
#' @examples
#' est <- get_estimator("relative")
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}}
#'
#' @export
get_estimator <- function(name) {
  stopifnot(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    available <- ls(envir = .estimator_registry)
    stop(
      sprintf("Estimator '%s' not found. Available: %s",
              name, paste(sort(available), collapse = ", ")),
      call. = FALSE
    )
  }
  get(name, envir = .estimator_registry, inherits = FALSE)
}


#' List All Registered Estimators
#'
#' @description
#' Return a data frame summarising all registered network estimators.
#'
#' @return A data frame with columns \code{name}, \code{description},
#'   \code{directed}.
#'
#' @examples
#' list_estimators()
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{get_estimator}}
#'
#' @export
list_estimators <- function() {
  nms <- ls(envir = .estimator_registry)
  if (length(nms) == 0) {
    return(data.frame(
      name = character(0),
      description = character(0),
      directed = logical(0),
      stringsAsFactors = FALSE
    ))
  }
  nms <- sort(nms)
  entries <- lapply(nms, function(nm) {
    e <- get(nm, envir = .estimator_registry, inherits = FALSE)
    data.frame(
      name = nm,
      description = e$description,
      directed = e$directed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, entries)
}


#' Remove a Registered Estimator
#'
#' @description
#' Remove a network estimator from the registry.
#'
#' @param name Character. Name of the estimator to remove.
#'
#' @return Invisible \code{NULL}.
#'
#' @examples
#' register_estimator("test_est", function(data, ...) diag(3),
#'   description = "test", directed = FALSE)
#' remove_estimator("test_est")
#'
#' @seealso \code{\link{register_estimator}}, \code{\link{list_estimators}}
#'
#' @export
remove_estimator <- function(name) {
  stopifnot(is.character(name), length(name) == 1)
  if (!exists(name, envir = .estimator_registry, inherits = FALSE)) {
    stop(
      sprintf("Estimator '%s' not found in registry.", name),
      call. = FALSE
    )
  }
  rm(list = name, envir = .estimator_registry)
  invisible(NULL)
}


#' Register Built-in Estimators
#' @noRd
.register_builtin_estimators <- function() {
  register_estimator("relative", .estimator_relative,
                     "Row-normalized transition probabilities", directed = TRUE)
  register_estimator("frequency", .estimator_frequency,
                     "Raw transition frequency counts", directed = TRUE)
  register_estimator("co_occurrence", .estimator_co_occurrence,
                     "Co-occurrence within sequences", directed = FALSE)
  register_estimator("cor", .estimator_cor,
                     "Pairwise correlation network", directed = FALSE)
  register_estimator("pcor", .estimator_pcor,
                     "Unregularized partial correlations", directed = FALSE)
  register_estimator("glasso", .estimator_glasso,
                     "EBICglasso regularized partial correlations",
                     directed = FALSE)
  register_estimator("ising", .estimator_ising,
                     "Ising model (L1-penalized logistic regression)",
                     directed = FALSE)
  register_estimator("mgm", .estimator_mgm,
                     "Mixed Graphical Model (nodewise lasso, EBIC, LW threshold)",
                     directed = FALSE)
  register_estimator("attention", .estimator_attention,
                     "Decay-weighted attention transitions", directed = TRUE)
  register_estimator("wtna", .estimator_wtna,
                     "Window-based TNA transitions (one-hot)", directed = TRUE)
  register_estimator("wtna_cooccurrence",
                     function(data, codes = NULL, window_size = 3L,
                              mode = "non-overlapping", actor = NULL,
                              type = "frequency", ...) {
                       .estimator_wtna_core(data, codes = codes,
                         window_size = window_size, mode = mode,
                         actor = actor, wtna_method = "cooccurrence",
                         type = type, ...)
                     },
                     "Window-based TNA co-occurrence (one-hot)", directed = FALSE)
}
