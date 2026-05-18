#' Estimate a Network (Deprecated)
#'
#' @description
#' This function is deprecated. Use \code{\link{build_network}} instead.
#'
#' @inheritParams build_network
#' @param method Character. Defaults to \code{"relative"} for backward
#'   compatibility.
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#'
#' @return A \code{netobject} (see \code{\link{build_network}}).
#'
#' @examples
#' data <- data.frame(A = c("x","y","z","x"), B = c("y","x","z","y"))
#' net <- estimate_network(data, method = "relative")
#'
#' @seealso \code{\link{build_network}}
#'
#' @importFrom stats aggregate ave cor complete.cases var
#' @export
estimate_network <- function(data,
                             method = "relative",
                             params = list(),
                             scaling = NULL,
                             threshold = 0,
                             level = NULL,
                             ...) {
  .Deprecated("build_network")
  build_network(
    data = data,
    method = method,
    params = params,
    scaling = scaling,
    threshold = threshold,
    level = level,
    ...
  )
}


# ---- Convenience wrappers for build_network ----

#' Build a Transition Network (TNA)
#'
#' Convenience wrapper for \code{build_network(method = "relative")}.
#' Computes row-normalized transition probabilities from sequence data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' net <- build_tna(seqs)
#' @export
build_tna <- function(data, ...) {
  build_network(data, method = "relative", ...)
}

#' Build a Frequency Transition Network (FTNA)
#'
#' Convenience wrapper for \code{build_network(method = "frequency")}.
#' Computes raw transition counts from sequence data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' net <- build_ftna(seqs)
#' @export
build_ftna <- function(data, ...) {
  build_network(data, method = "frequency", ...)
}

#' Build an Attention-Weighted Transition Network (ATNA)
#'
#' Convenience wrapper for \code{build_network(method = "attention")}.
#' Computes decay-weighted transitions from sequence data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' net <- build_atna(seqs)
#' @export
build_atna <- function(data, ...) {
  build_network(data, method = "attention", ...)
}

#' Build a Co-occurrence Network (CNA)
#'
#' Convenience wrapper for \code{build_network(method = "co_occurrence")}.
#' Computes co-occurrence counts from binary or sequence data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}, \code{\link{cooccurrence}} for
#'   delimited-field, bipartite, and other non-sequence co-occurrence formats.
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C"), V2 = c("B","C","A"))
#' net <- build_cna(seqs)
#' @export
build_cna <- function(data, ...) {
  build_network(data, method = "co_occurrence", ...)
}

#' Build a Correlation Network
#'
#' Convenience wrapper for \code{build_network(method = "cor")}.
#' Computes Pearson correlations from numeric data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_cor(srl_strategies)
#' @export
build_cor <- function(data, ...) {
  build_network(data, method = "cor", ...)
}

#' Build a Partial Correlation Network
#'
#' Convenience wrapper for \code{build_network(method = "pcor")}.
#' Computes partial correlations from numeric data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_pcor(srl_strategies)
#' @export
build_pcor <- function(data, ...) {
  build_network(data, method = "pcor", ...)
}

#' Build a Graphical Lasso Network (EBICglasso)
#'
#' Convenience wrapper for \code{build_network(method = "glasso")}.
#' Computes L1-regularized partial correlations with EBIC model selection.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' data(srl_strategies)
#' net <- build_glasso(srl_strategies)
#' @export
build_glasso <- function(data, ...) {
  build_network(data, method = "glasso", ...)
}

#' Build an Ising Network
#'
#' Convenience wrapper for \code{build_network(method = "ising")}.
#' Computes L1-regularized logistic regression network for binary data.
#'
#' @inheritParams build_network
#' @param ... Additional arguments passed to \code{\link{build_network}}.
#' @return A \code{netobject} (see \code{\link{build_network}}).
#' @seealso \code{\link{build_network}}
#' @examples
#' \donttest{
#' bin_data <- data.frame(matrix(rbinom(200, 1, 0.5), ncol = 5))
#' net <- build_ising(bin_data)
#' }
#' @export
build_ising <- function(data, ...) {
  build_network(data, method = "ising", ...)
}


# ---- Shared internal helpers ----
# These are used by build_network(), bootstrap_network(), and other functions.


# ---- Method alias resolution ----

#' Resolve method aliases to canonical names
#' @noRd
.resolve_method_alias <- function(method) {
  aliases <- c(
    ebicglasso        = "glasso",
    regularized       = "glasso",
    partial           = "pcor",
    correlation       = "cor",
    corr              = "cor",
    transition        = "relative",
    tna               = "relative",
    counts            = "frequency",
    ftna              = "frequency",
    cna               = "co_occurrence",
    wcna              = "co_occurrence",
    wtna_transition   = "wtna",
    wtna_cooccurrence = "wtna_cooccurrence",
    isingfit          = "ising",
    atna              = "attention",
    mixed_graphical   = "mgm",
    mixed             = "mgm"
  )
  if (method %in% names(aliases)) {
    aliases[[method]]
  } else {
    method
  }
}


# ---- Post-estimation scaling ----

#' Apply scaling transformations to a network matrix
#' @noRd
.apply_scaling <- function(mat, scaling, include_zeros = FALSE) {
  for (s in scaling) {
    mat <- switch(s,
      minmax = {
        vals <- if (isTRUE(include_zeros)) as.vector(mat) else mat[mat != 0]
        if (length(vals) == 0) {
          mat
        } else {
          rng <- range(vals, na.rm = TRUE)
          if (rng[1] == rng[2]) mat
          else {
            if (isTRUE(include_zeros)) {
              mat[] <- (mat - rng[1]) / (rng[2] - rng[1])
            } else {
              # Structure-preserving min-max for association networks: the
              # smallest genuine edge must not collide with the structural-
              # zero value (which would drop a real edge and corrupt
              # n_edges). Map the non-zero block into the strictly positive
              # interval (eps, 1] so the block minimum becomes eps, not 0.
              eps <- sqrt(.Machine$double.eps)
              nz <- mat != 0
              mat[nz] <- eps + (1 - eps) *
                (mat[nz] - rng[1]) / (rng[2] - rng[1])
            }
            mat
          }
        }
      },
      max = {
        max_abs <- max(abs(mat))
        if (max_abs > 0) mat / max_abs else mat # nocov
      },
      rank = {
        # Rank scaling always preserves structural zeros: a non-edge must
        # remain a non-edge after scaling. include_zeros = TRUE only affects
        # min-max / range scalings; for ranks it would promote absent
        # transitions to positive weights and densify the network.
        nz <- mat != 0
        if (any(nz)) {
          mat[nz] <- rank(mat[nz], ties.method = "average")
          mat
        } else {
          mat
        }
      },
      normalize = {
        rs <- rowSums(abs(mat))
        nonzero_rows <- rs > 0
        mat[nonzero_rows, ] <- mat[nonzero_rows, ] / rs[nonzero_rows]
        mat
      },
      mat  # default: no change
    )
  }
  mat
}


# ---- Edge extraction ----

#' Extract non-zero edges from a network matrix
#'
#' For undirected networks, uses upper triangle only.
#' For directed networks, uses all non-diagonal non-zero entries.
#'
#' @noRd
.extract_edges_from_matrix <- function(mat, directed = FALSE) {
  if (directed) {
    # Keep self-loops too: every non-zero entry is a real edge.
    idx <- which(mat != 0, arr.ind = TRUE)
  } else {
    # row <= col keeps the upper triangle PLUS the diagonal (one row
    # per self-loop, no double-count for undirected networks).
    idx <- which(mat != 0 & row(mat) <= col(mat), arr.ind = TRUE)
  }

  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from   = as.integer(idx[, 1]),
    to     = as.integer(idx[, 2]),
    weight = mat[idx],
    stringsAsFactors = FALSE
  )
}


# ---- Multilevel decomposition ----

#' Decompose data for multilevel analysis
#'
#' @param data Data frame.
#' @param id_col Character. Grouping variable name.
#' @param level Character: "between" or "within".
#'
#' @return Transformed data frame.
#' @noRd
.decompose_multilevel <- function(data, id_col, level) {
  stopifnot(is.data.frame(data)) # nocov
  grp_var <- id_col[1]

  if (!grp_var %in% names(data)) {
    stop("id_col '", grp_var, "' not found in data.", call. = FALSE)
  }

  # Get numeric columns (exclude id columns and "rid")
  exclude <- c(id_col, "rid")
  numeric_cols <- vapply(data, is.numeric, logical(1))
  keep <- setdiff(names(data)[numeric_cols], exclude)

  if (length(keep) < 2) {
    stop("At least 2 numeric columns are required for multilevel decomposition.")
  }

  mat <- data[, keep, drop = FALSE]
  id_vals <- data[[grp_var]]

  if (level == "between") {
    # Aggregate to person means
    mat$.id <- id_vals
    agg <- aggregate(. ~ .id, data = mat, FUN = mean)
    result <- agg[, names(agg) != ".id", drop = FALSE]
    return(as.data.frame(result))

  } else if (level == "within") {
    # Drop persons with < 2 observations
    tab <- table(id_vals)
    multi <- names(tab[tab >= 2])
    keep_rows <- id_vals %in% multi
    if (any(!keep_rows)) {
      n_single <- sum(!keep_rows)
      message("Dropping ", n_single,
              " single-observation rows (within-person centering).")
      mat <- mat[keep_rows, , drop = FALSE]
      id_vals <- id_vals[keep_rows]
    }

    if (nrow(mat) < 3) {
      stop("Fewer than 3 rows remain after dropping ",
           "single-observation persons.")
    }

    # Person-mean center each variable
    mat_m <- as.matrix(mat)
    for (j in seq_len(ncol(mat_m))) {
      mat_m[, j] <- mat_m[, j] - ave(mat_m[, j], id_vals, FUN = mean)
    }

    return(as.data.frame(mat_m))
  }

  data
}
