# ---- Sequence Clustering ----

# Available metrics and methods (private constants)
.clustering_metrics <- c(
  "hamming", "osa", "lv", "dl", "lcs", "qgram", "cosine", "jaccard", "jw"
)
.clustering_methods <- c(
  "pam", "ward.D2", "ward.D", "complete", "average",
  "single", "mcquitty", "median", "centroid"
)

# ==============================================================================
# 1. Sequence Encoding
# ==============================================================================

#' Encode sequences as integer matrix
#'
#' @param data Data frame or matrix of sequences (rows = sequences, cols = time)
#' @param na_syms Character vector of symbols to treat as NA
#' @return List with int_mat, len, states, n_states
#' @noRd
.encode_sequences <- function(data, na_syms = c("*", "%")) {
  mat <- as.matrix(data)
  storage.mode(mat) <- "character"
  mat[mat %in% na_syms] <- NA
  obs <- !is.na(mat)
  last_obs <- max.col(obs, ties.method = "last")
  all_na <- rowSums(obs) == 0L
  last_obs[all_na] <- 0L
  states <- sort(unique(mat[obs]))
  n_states <- length(states)
  int_mat <- matrix(match(mat, states), nrow = nrow(mat), ncol = ncol(mat))
  list(int_mat = int_mat, len = last_obs, states = states, n_states = n_states)
}

# ==============================================================================
# 2. Distance Functions (R fallback)
# ==============================================================================

#' @noRd
.hamming_dist_r <- function(x, y, n, m, weights, ...) {
  sum((x != y) * weights)
}

#' @noRd
.levenshtein_dist_r <- function(x, y, n, m, ...) {
  d <- matrix(0L, n + 1L, m + 1L)
  d[, 1L] <- 0:n
  d[1L, ] <- 0:m
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      cost <- as.integer(x[i] != y[j])
      d[i + 1L, j + 1L] <- min(
        d[i, j + 1L] + 1L,
        d[i + 1L, j] + 1L,
        d[i, j] + cost
      )
    }
  }
  d[n + 1L, m + 1L]
}

#' @noRd
.osa_dist_r <- function(x, y, n, m, ...) {

  d <- matrix(0L, n + 1L, m + 1L)
  d[, 1L] <- 0:n
  d[1L, ] <- 0:m
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      cost <- as.integer(x[i] != y[j])
      d[i + 1L, j + 1L] <- min(
        d[i, j + 1L] + 1L,
        d[i + 1L, j] + 1L,
        d[i, j] + cost
      )
      if (i > 1L && j > 1L && x[i] == y[j - 1L] && x[i - 1L] == y[j]) {
        d[i + 1L, j + 1L] <- min(d[i + 1L, j + 1L], d[i - 1L, j - 1L] + 1L)
      }
    }
  }
  d[n + 1L, m + 1L]
}

#' @noRd
.dl_dist_r <- function(x, y, n, m, ...) {
  if (n == 0L) return(m)
  if (m == 0L) return(n)
  x_sub <- x[seq_len(n)]
  y_sub <- y[seq_len(m)]
  alphabet <- unique(c(x_sub, y_sub))
  last <- stats::setNames(rep(0L, length(alphabet)), as.character(alphabet))
  d <- matrix(0L, n + 2L, m + 2L)
  maxdist <- n + m
  d[1L, ] <- maxdist
  d[, 1L] <- maxdist
  for (i in 0:n) d[i + 2L, 2L] <- i
  for (j in 0:m) d[2L, j + 2L] <- j
  for (i in seq_len(n)) {
    db <- 0L
    for (j in seq_len(m)) {
      k <- last[[as.character(y[j])]]
      if (is.null(k)) k <- 0L # nocov
      l <- db
      if (x[i] == y[j]) {
        cost <- 0L
        db <- j
      } else {
        cost <- 1L
      }
      d[i + 2L, j + 2L] <- min(
        d[i + 1L, j + 2L] + 1L,
        d[i + 2L, j + 1L] + 1L,
        d[i + 1L, j + 1L] + cost,
        d[k + 1L, l + 1L] + (i - k - 1L) + (j - l - 1L) + 1L
      )
    }
    last[[as.character(x[i])]] <- i
  }
  d[n + 2L, m + 2L]
}

#' @noRd
.lcs_dist_r <- function(x, y, n, m, ...) {
  lcs_len <- matrix(0L, n + 1L, m + 1L)
  for (i in seq_len(n)) {
    for (j in seq_len(m)) {
      if (x[i] == y[j]) {
        lcs_len[i + 1L, j + 1L] <- lcs_len[i, j] + 1L
      } else {
        lcs_len[i + 1L, j + 1L] <- max(lcs_len[i, j + 1L], lcs_len[i + 1L, j])
      }
    }
  }
  n + m - 2L * lcs_len[n + 1L, m + 1L]
}

#' Extract q-gram frequency profile from integer sequence
#' @noRd
.get_qgram_r <- function(x, n, q) {
  if (n < q) return(integer()) # nocov start
  ng <- n - q + 1L
  grams <- character(ng)
  for (i in seq_len(ng)) {
    grams[i] <- paste0(x[i:(i + q - 1L)], collapse = "\x01")
  }
  tab <- table(grams)
  counts <- as.integer(tab)
  names(counts) <- names(tab)
  counts # nocov end
}

#' @noRd
.qgram_dist_r <- function(x, y, qx, qy, ...) {
  grams <- union(names(qx), names(qy))
  x_all <- qx[grams]
  y_all <- qy[grams]
  x_all[is.na(x_all)] <- 0L
  y_all[is.na(y_all)] <- 0L
  sum(abs(x_all - y_all))
}

#' @noRd
.cosine_dist_r <- function(x, y, qx, qy, ...) {
  grams <- union(names(qx), names(qy))
  x_all <- qx[grams]
  y_all <- qy[grams]
  x_all[is.na(x_all)] <- 0L
  y_all[is.na(y_all)] <- 0L
  num <- sum(x_all * y_all)
  den <- sqrt(sum(x_all^2)) * sqrt(sum(y_all^2))
  if (den == 0) return(1)
  1 - num / den
}

#' @noRd
.jaccard_dist_r <- function(qx, qy, ...) {
  qx_names <- names(qx)
  qy_names <- names(qy)
  if (length(qx_names) == 0L && length(qy_names) == 0L) return(0)
  1 - length(intersect(qx_names, qy_names)) / length(union(qx_names, qy_names))
}

#' @noRd
.jaro_dist_r <- function(x, y, n, m, ...) {
  if (n == 0L && m == 0L) return(0) # nocov start
  if (n == 0L || m == 0L) return(1) # nocov end
  match_dist <- max(floor(max(n, m) / 2) - 1L, 0L)
  x_match <- rep(FALSE, n)
  y_match <- rep(FALSE, m)
  matches <- 0L
  for (i in seq_len(n)) {
    start <- max(1L, i - match_dist)
    end <- min(m, i + match_dist)
    if (start > end) next # nocov
    for (j in start:end) {
      if (!y_match[j] && x[i] == y[j]) {
        x_match[i] <- TRUE
        y_match[j] <- TRUE
        matches <- matches + 1L
        break
      }
    }
  }
  if (matches == 0L) return(1) # nocov
  trans <- 0L
  k <- 1L
  for (i in seq_len(n)) {
    if (x_match[i]) {
      while (k <= m && !y_match[k]) k <- k + 1L
      if (k <= m && x[i] != y[k]) trans <- trans + 1L
      k <- k + 1L
    }
  }
  trans <- trans / 2L
  jaro_sim <- (matches / n + matches / m + (matches - trans) / matches) / 3
  1 - jaro_sim
}

#' @noRd
.jw_dist_r <- function(x, y, n, m, p = 0.1, max_l = 4L, ...) {
  jaro_sim <- 1 - .jaro_dist_r(x, y, n, m)
  if (jaro_sim == 0) return(1) # nocov
  pref <- min(n, m, max_l)
  l <- 0L
  for (k in seq_len(pref)) {
    if (x[k] == y[k]) {
      l <- l + 1L
    } else {
      break
    }
  }
  p <- min(p, 0.25)
  jw <- jaro_sim + l * p * (1 - jaro_sim)
  1 - jw
}

# ==============================================================================
# 3. Vectorized Distance Matrix (hamming, qgram-family)
# ==============================================================================

#' Hamming distance matrix via column-wise outer() — no pair loop
#' @noRd
.hamming_matrix_outer <- function(mat_s, weights) {
  n <- nrow(mat_s)
  k <- ncol(mat_s)
  d <- matrix(0, n, n)
  if (length(weights) == 1L) {
    # Unweighted: just count mismatches per column
    for (c in seq_len(k)) {
      d <- d + outer(mat_s[, c], mat_s[, c], "!=")
    }
  } else {
    for (c in seq_len(k)) {
      d <- d + outer(mat_s[, c], mat_s[, c], "!=") * weights[c]
    }
  }
  stats::as.dist(d)
}

#' Build term-document matrix of q-gram counts (rows = seqs, cols = qgrams)
#' @noRd
.build_qgram_tdm <- function(mat_s, len, q) {
  n <- nrow(mat_s)
  # Collect all q-gram profiles
  profiles <- vector("list", n)
  all_grams <- character(0L)
  for (i in seq_len(n)) {
    if (len[i] >= q) {
      ng <- len[i] - q + 1L
      grams <- character(ng)
      for (j in seq_len(ng)) {
        grams[j] <- paste0(mat_s[i, j:(j + q - 1L)], collapse = "\x01")
      }
      tab <- table(grams)
      profiles[[i]] <- tab
      all_grams <- union(all_grams, names(tab))
    }
  }
  # Build matrix: rows = sequences, cols = unique qgrams
  tdm <- matrix(0L, nrow = n, ncol = length(all_grams))
  colnames(tdm) <- all_grams
  for (i in seq_len(n)) {
    if (!is.null(profiles[[i]])) {
      tdm[i, names(profiles[[i]])] <- as.integer(profiles[[i]])
    }
  }
  tdm
}

#' Q-gram distance matrix via term-document matrix subtraction
#' @noRd
.qgram_matrix_tdm <- function(tdm) {
  n <- nrow(tdm)
  # |a_i - b_i| summed over all qgrams = Manhattan distance between rows
  # Expand via outer difference per column, then sum absolute values
  d <- matrix(0, n, n)
  for (c in seq_len(ncol(tdm))) {
    d <- d + abs(outer(tdm[, c], tdm[, c], "-"))
  }
  stats::as.dist(d)
}

#' Cosine distance matrix via matrix multiply
#' @noRd
.cosine_matrix_tdm <- function(tdm) {
  # cosine_sim = (A %*% t(A)) / (norm_a * norm_b)
  norms <- sqrt(rowSums(tdm^2))
  norms[norms == 0] <- 1  # avoid division by zero (returns dist=1 for zero vecs)
  dot <- tcrossprod(tdm)
  sim <- dot / outer(norms, norms)
  # Clamp to [0, 1] for numerical safety
  sim[sim > 1] <- 1
  sim[sim < 0] <- 0
  # Handle zero-norm rows: distance = 1
  zero_rows <- which(rowSums(tdm^2) == 0)
  d <- 1 - sim
  diag(d) <- 0
  if (length(zero_rows) > 0L) {
    d[zero_rows, ] <- 1
    d[, zero_rows] <- 1
    d[cbind(zero_rows, zero_rows)] <- 0
  }
  stats::as.dist(d)
}

#' Jaccard distance matrix via set intersection/union on qgram presence
#' @noRd
.jaccard_matrix_tdm <- function(tdm) {
  # Binary presence matrix
  bin <- (tdm > 0L) * 1L
  # intersection = A %*% t(B), union = |A| + |B| - intersection
  isect <- tcrossprod(bin)
  sizes <- rowSums(bin)
  union_mat <- outer(sizes, sizes, "+") - isect
  # Avoid 0/0 (both empty → distance 0)
  d <- ifelse(union_mat == 0, 0, 1 - isect / union_mat)
  diag(d) <- 0
  stats::as.dist(d)
}

# ==============================================================================
# 4. Pairwise Distance Matrix (DP distances — pair loop unavoidable)
# ==============================================================================

#' Compute pairwise dissimilarity matrix using pure R
#' @noRd
.dissimilarity_matrix_r <- function(enc, dissimilarity, lambda,
                                    q = 2L, p = 0.1) {
  int_mat <- enc$int_mat
  len <- enc$len
  n <- nrow(int_mat)
  k <- ncol(int_mat)
  sentinel <- enc$n_states + 1L

  # Full-width matrix with sentinel for NA
  mat_s <- int_mat
  mat_s[is.na(mat_s)] <- sentinel

  # --- Vectorized fast paths ---
  if (dissimilarity == "hamming") {
    weights <- 1
    if (lambda > 0) {
      weights <- exp(-lambda * seq(0, k - 1))
      weights <- weights / max(weights)
    }
    return(.hamming_matrix_outer(mat_s, weights))
  }

  if (dissimilarity %in% c("qgram", "cosine", "jaccard")) {
    tdm <- .build_qgram_tdm(mat_s, len, q)
    if (dissimilarity == "qgram") return(.qgram_matrix_tdm(tdm))
    if (dissimilarity == "cosine") return(.cosine_matrix_tdm(tdm))
    return(.jaccard_matrix_tdm(tdm))
  }

  # --- Pair loop for DP distances (osa, lv, dl, lcs, jw) ---
  d_fun <- switch(dissimilarity,
    osa = .osa_dist_r,
    lv = .levenshtein_dist_r,
    dl = .dl_dist_r,
    lcs = .lcs_dist_r,
    jw = .jw_dist_r
  )

  d <- matrix(0, n, n)
  pairs <- combn(n, 2L)

  dists <- vapply(seq_len(ncol(pairs)), function(idx) {
    i <- pairs[1L, idx]
    j <- pairs[2L, idx]
    if (dissimilarity == "jw") {
      d_fun(x = mat_s[i, ], y = mat_s[j, ],
            n = len[i], m = len[j], p = p)
    } else {
      d_fun(x = mat_s[i, ], y = mat_s[j, ],
            n = len[i], m = len[j])
    }
  }, numeric(1L))

  d[lower.tri(d)] <- dists
  d <- d + t(d)
  stats::as.dist(d)
}

# ==============================================================================
# 4. stringdist Fast Path
# ==============================================================================

#' Compute pairwise dissimilarity matrix using stringdist (C-level)
#' @noRd
.dissimilarity_matrix_stringdist <- function(enc, dissimilarity, lambda,
                                             q = 2L, p = 0.1) {
  int_mat <- enc$int_mat
  len <- enc$len
  n_states <- enc$n_states
  k <- ncol(int_mat)

  base_chars <- c(letters, LETTERS, as.character(0:9))
  sentinel_char <- base_chars[n_states + 1L]

  # Build character matrix (sentinel for NA)
  chr_mat <- matrix(sentinel_char, nrow = nrow(int_mat), ncol = k)
  obs <- !is.na(int_mat)
  chr_mat[obs] <- base_chars[int_mat[obs]]

  if (dissimilarity == "hamming") {
    # Weighted hamming not supported by stringdist — fall back to R
    if (lambda > 0) {
      return(.dissimilarity_matrix_r(enc, dissimilarity, lambda, q, p)) # nocov
    }
    # Full width strings (same length, NAs are sentinel chars)
    strings <- vapply(seq_len(nrow(chr_mat)), function(i) {
      paste0(chr_mat[i, ], collapse = "")
    }, character(1L))
    stringdist::stringdistmatrix(strings, method = "hamming")
  } else {
    # Truncate to effective length per sequence
    strings <- vapply(seq_len(nrow(chr_mat)), function(i) {
      if (len[i] == 0L) return("") # nocov
      paste0(chr_mat[i, seq_len(len[i])], collapse = "")
    }, character(1L))
    if (dissimilarity == "jw") {
      stringdist::stringdistmatrix(strings, method = "jw", p = p)
    } else if (dissimilarity %in% c("qgram", "cosine", "jaccard")) {
      stringdist::stringdistmatrix(strings, method = dissimilarity, q = q)
    } else {
      stringdist::stringdistmatrix(strings, method = dissimilarity)
    }
  }
}

# ==============================================================================
# 5. Input Extraction
# ==============================================================================

#' Extract wide-format sequence data from various object types
#'
#' Supports: data.frame, matrix, netobject, tna, cograph_network.
#' For tna objects, decodes integer-encoded data using stored labels.
#' Rejects association-method netobjects (cor/pcor/glasso/ising) since
#' those contain numeric observations, not categorical sequences.
#' @param data Input object
#' @return A data frame of wide-format sequences
#' @noRd
.extract_sequence_data <- function(data) {
  # --- netobject (from Nestimate) ---
  if (inherits(data, "netobject")) {
    seq_methods <- c("relative", "frequency", "co_occurrence", "attention")
    if (!is.null(data$method) && !(data$method %in% seq_methods)) {
      stop(sprintf(
        "build_clusters() requires sequence data, but this netobject uses method '%s'. ",
        data$method
      ), "Pass the original sequence data frame instead.", call. = FALSE)
    }
    if (is.null(data$data)) { # nocov start
      stop("netobject has no $data. Pass the original sequence data.", call. = FALSE)
    } # nocov end
    return(as.data.frame(data$data, stringsAsFactors = FALSE))
  }

  # --- tna model (integer-encoded with $labels) ---
  if (inherits(data, "tna")) {
    if (is.null(data$data) || is.null(data$labels)) { # nocov start
      stop("tna object missing $data or $labels.", call. = FALSE)
    } # nocov end
    # Decode integer matrix back to state names
    decoded <- matrix(data$labels[data$data], nrow = nrow(data$data),
                      ncol = ncol(data$data))
    colnames(decoded) <- colnames(data$data)
    return(as.data.frame(decoded, stringsAsFactors = FALSE))
  }

  # --- cograph_network ---
  if (inherits(data, "cograph_network")) {
    if (is.null(data$data)) { # nocov start
      stop("cograph_network has no $data. Pass the original sequence data.",
           call. = FALSE)
    } # nocov end
    return(as.data.frame(data$data, stringsAsFactors = FALSE))
  }

  # --- raw data.frame or matrix ---
  if (is.data.frame(data) || is.matrix(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }

  stop("Unsupported input type: ", paste(class(data), collapse = ", "),
       ". Expected data.frame, matrix, netobject, tna, or cograph_network.",
       call. = FALSE)
}

# ==============================================================================
# 6. Main Function
# ==============================================================================

#' Cluster Sequences by Dissimilarity
#'
#' @description
#' Clusters wide-format sequences using pairwise string dissimilarity and
#' either PAM (Partitioning Around Medoids) or hierarchical clustering.
#' Supports 9 distance metrics including temporal weighting for Hamming
#' distance. When the \pkg{stringdist} package is available, uses C-level
#' distance computation for 100-1000x speedup on edit distances.
#'
#' @param data Input data. Accepts multiple formats:
#' \describe{
#'   \item{data.frame / matrix}{Wide-format sequences (rows = sequences,
#'     columns = time points, values = state names).}
#'   \item{netobject}{A network object from \code{\link{build_network}}.
#'     Extracts the stored sequence data. Only valid for sequence-based
#'     methods (relative, frequency, co_occurrence, attention).}
#'   \item{tna}{A tna model from the tna package. Decodes the
#'     integer-encoded sequence data using stored labels.}
#'   \item{cograph_network}{A cograph network object. Extracts the
#'     stored sequence data.}
#' }
#' @param k Integer. Number of clusters (must be between 2 and
#'   \code{nrow(data) - 1}).
#' @param dissimilarity Character. Distance metric. One of \code{"hamming"},
#'   \code{"osa"} (optimal string alignment), \code{"lv"} (Levenshtein),
#'   \code{"dl"} (Damerau-Levenshtein), \code{"lcs"} (longest common
#'   subsequence), \code{"qgram"}, \code{"cosine"}, \code{"jaccard"},
#'   \code{"jw"} (Jaro-Winkler). Default: \code{"hamming"}.
#' @param method Character. Clustering method. \code{"pam"} for Partitioning
#'   Around Medoids, or a hierarchical method: \code{"ward.D2"},
#'   \code{"ward.D"}, \code{"complete"}, \code{"average"}, \code{"single"},
#'   \code{"mcquitty"}, \code{"median"}, \code{"centroid"}.
#'   Default: \code{"pam"}.
#' @param na_syms Character vector. Symbols treated as missing values.
#'   Default: \code{c("*", "\%")}.
#'
#'   \strong{Missing-value distance rule:} after symbols are converted to
#'   \code{NA}, missing values are encoded as a single comparable sentinel
#'   state -- \emph{not} pairwise-deleted. Two missing values in the same
#'   position match (distance contribution 0); a missing value paired with
#'   any observed state mismatches (distance contribution 1 for Hamming,
#'   etc.). This is the conventional behaviour for aligned sequence
#'   matrices because pairwise deletion would change the effective length
#'   of every pair and break the metric. If you want pairwise deletion or
#'   a different missing-value semantic, drop or recode the missing cells
#'   before passing the data in.
#' @param weighted Logical. Apply exponential decay weighting to Hamming
#'   distance positions? Only valid when \code{dissimilarity = "hamming"}.
#'   Default: \code{FALSE}.
#' @param lambda Numeric. Non-negative decay rate for weighted Hamming.
#'   Higher values weight earlier positions more strongly. Default: 1.
#' @param seed Integer or NULL. Random seed for reproducibility. Default:
#'   \code{NULL}.
#' @param q Integer. Size of q-grams for \code{"qgram"}, \code{"cosine"},
#'   and \code{"jaccard"} distances. Default: \code{2L}.
#' @param p Numeric. Winkler prefix penalty for Jaro-Winkler distance.
#'   Must be between 0 and 0.25. Default: \code{0.1}.
#' @param covariates Optional. Post-hoc covariate analysis of cluster
#'   membership. Accepts:
#'   \describe{
#'     \item{string}{Single column name, e.g. \code{"Age"}. Resolved
#'       against \code{x$metadata} (and \code{x$data}) for
#'       \code{netobject} or \code{cograph_network} input.}
#'     \item{character vector}{\code{c("Age", "Gender")}, same lookup.}
#'     \item{formula}{\code{~ Age + Gender}, same lookup; supports
#'       \code{"Age + Gender"} string form too.}
#'     \item{data.frame}{All columns used as covariates verbatim;
#'       must have one row per sequence.}
#'     \item{NULL}{No covariate analysis (default).}
#'   }
#'   For \code{netobject} or \code{cograph_network} input, names are
#'   resolved against \code{$metadata} first and then non-state columns
#'   of \code{$data}, so a typical call looks like
#'   \code{build_clusters(net, k = 3, covariates = "session_label")}
#'   without pre-extracting a data.frame. \code{tna} input requires the
#'   data.frame form. Results are stored in \code{$covariates}.
#' @param estimator Multinomial logit fitter for the covariate analysis.
#'   \code{"firth"} (default) uses Firth's penalised likelihood via
#'   \code{brglm2::brmultinom} — bias-reduced and finite under
#'   quasi-complete separation. \code{"multinom"} uses classical ML via
#'   \code{nnet::multinom}; emits a warning because rare-cell separation
#'   produces astronomical ORs with degenerate CIs (silent failure).
#'   \code{"chisq"} runs WeightedCluster-style descriptive tests
#'   (chi-square + Cramer's V + standardized adjusted residuals for
#'   factors; Kruskal-Wallis + eta-squared for numerics).
#' @param ... Unsupported. Supplying unused arguments raises an error.
#'
#' @return An object of class \code{"net_clustering"} containing:
#' \describe{
#'   \item{data}{The original input data.}
#'   \item{k}{Number of clusters.}
#'   \item{assignments}{Named integer vector of cluster assignments.}
#'   \item{silhouette}{Overall average silhouette width.}
#'   \item{sizes}{Named integer vector of cluster sizes.}
#'   \item{method}{Clustering method used.}
#'   \item{dissimilarity}{Distance metric used.}
#'   \item{distance}{The computed dissimilarity matrix (\code{dist} object).}
#'   \item{medoids}{Integer vector of medoid row indices (PAM only; NULL
#'     for hierarchical methods).}
#'   \item{seed}{Seed used (or NULL).}
#'   \item{weighted}{Logical, whether weighted Hamming was used.}
#'   \item{lambda}{Lambda value used (0 if not weighted).}
#' }
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' cl <- build_clusters(seqs, k = 2)
#' cl
#' \donttest{
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:3], 20, TRUE), V2 = sample(LETTERS[1:3], 20, TRUE),
#'   V3 = sample(LETTERS[1:3], 20, TRUE), V4 = sample(LETTERS[1:3], 20, TRUE)
#' )
#' cl <- build_clusters(seqs, k = 2)
#' print(cl)
#' summary(cl)
#' }
#'
#' @importFrom utils combn
#' @export
build_clusters <- function(data, k, dissimilarity = "hamming", method = "pam",
                         na_syms = c("*", "%"), weighted = FALSE, lambda = 1,
                         seed = NULL, q = 2L, p = 0.1,
                         covariates = NULL,
                         estimator = c("firth", "multinom", "chisq"),
                         ...) {
  estimator <- match.arg(estimator)
  dots <- list(...)
  if (length(dots) > 0L) {
    dot_names <- names(dots)
    dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
    stop("build_clusters() got unsupported argument",
         if (length(dots) == 1L) ": " else "s: ",
         paste(dot_names, collapse = ", "), call. = FALSE)
  }

  # --- Extract sequence data from supported objects ---
  raw_data <- data
  data <- .extract_sequence_data(data)

  # --- Input validation ---
  # Use named-condition stopifnot for type checks so the failure message
  # names the offending argument instead of dumping the predicate.
  stopifnot(
    "'k' must be a single numeric value"            = is.numeric(k) && length(k) == 1L,
    "'dissimilarity' must be a single character"    = is.character(dissimilarity) && length(dissimilarity) == 1L,
    "'method' must be a single character"           = is.character(method) && length(method) == 1L,
    "'na_syms' must be a character vector"          = is.character(na_syms),
    "'weighted' must be a single logical"           = is.logical(weighted) && length(weighted) == 1L,
    "'lambda' must be a single numeric"             = is.numeric(lambda) && length(lambda) == 1L,
    "'q' must be a single numeric"                  = is.numeric(q) && length(q) == 1L,
    "'p' must be a single numeric"                  = is.numeric(p) && length(p) == 1L
  )
  if (!is.finite(k) || k != floor(k)) {
    stop("'k' must be a whole finite number.", call. = FALSE)
  }
  if (is.na(weighted)) {
    stop("'weighted' must be TRUE or FALSE, not NA.", call. = FALSE)
  }
  if (!is.finite(lambda) || lambda < 0) {
    stop("'lambda' must be a finite non-negative number.", call. = FALSE)
  }
  if (!is.finite(q) || q != floor(q) || q < 1) {
    stop("'q' must be a positive whole number.", call. = FALSE)
  }
  if (!is.finite(p) || p < 0 || p > 0.25) {
    stop("'p' must be a finite number between 0 and 0.25.", call. = FALSE)
  }
  k <- as.integer(k)
  q <- as.integer(q)
  n <- nrow(data)
  if (k < 2L) {
    stop(sprintf(
      "'k' must be at least 2 (got k = %d). Clustering with one group is not meaningful.",
      k
    ), call. = FALSE)
  }
  if (k > n - 1L) {
    stop(sprintf(
      paste0("'k' (= %d) must be <= n - 1 where n = %d sequences. ",
             "Either reduce k or supply more sequences (need at least %d)."),
      k, n, k + 1L
    ), call. = FALSE)
  }

  dissimilarity <- match.arg(dissimilarity, .clustering_metrics)
  method <- match.arg(method, .clustering_methods)

  if (weighted && dissimilarity != "hamming") {
    stop("Weighting is only supported for Hamming distance.", call. = FALSE)
  }
  lambda <- if (weighted) lambda else 0

  if (!is.null(seed)) set.seed(seed)

  # --- Encode sequences ---
  enc <- .encode_sequences(data, na_syms)

  if (enc$n_states == 0L) {
    stop("No observed sequence states were found after applying na_syms = ",
         sprintf("c(%s)", paste(shQuote(na_syms), collapse = ", ")),
         ". Check that the input data contains at least one non-missing ",
         "state cell, or adjust na_syms to exclude symbols that should ",
         "be treated as observed states.", call. = FALSE)
  }

  # --- Compute distance matrix ---
  use_stringdist <- requireNamespace("stringdist", quietly = TRUE) &&
    enc$n_states < 62L &&  # need room for sentinel char
    !(dissimilarity == "hamming" && lambda > 0)  # weighted hamming not in stringdist

  if (use_stringdist) {
    dist_mat <- .dissimilarity_matrix_stringdist(
      enc, dissimilarity, lambda, q, p
    )
  } else {
    dist_mat <- .dissimilarity_matrix_r(enc, dissimilarity, lambda, q, p)
  }

  # --- Cluster ---
  medoids <- NULL
  if (method == "pam") {
    clust_result <- cluster::pam(dist_mat, diss = TRUE, k = k)
    assignments <- clust_result$clustering
    silhouette_score <- clust_result$silinfo$avg.width
    medoids <- clust_result$id.med
  } else {
    hc <- stats::hclust(dist_mat, method = method)
    assignments <- stats::cutree(hc, k = k)
    sil <- cluster::silhouette(assignments, dist = dist_mat)
    silhouette_score <- mean(sil[, 3L])
  }

  sizes <- c(table(assignments))

  # --- Post-hoc covariate analysis ---
  cov_result <- NULL
  cov_resolved <- .resolve_covariates(covariates, raw_data, n)
  if (!is.null(cov_resolved)) {
    cov_result <- .run_covariate_analysis(
      assignments, cov_resolved$cov_df, cov_resolved$rhs, k,
      estimator = estimator
    )
  }

  structure(
    list(
      data = data,
      k = k,
      assignments = assignments,
      silhouette = silhouette_score,
      sizes = sizes,
      method = method,
      dissimilarity = dissimilarity,
      distance = dist_mat,
      medoids = medoids,
      seed = seed,
      weighted = weighted,
      lambda = lambda,
      covariates = cov_result,
      network_method = if (inherits(raw_data, "netobject")) raw_data$method else NULL,
      build_args     = if (inherits(raw_data, "netobject")) raw_data$build_args else NULL
    ),
    class = "net_clustering"
  )
}

#' Cluster sequence data (deprecated alias)
#'
#' Renamed to \code{\link{build_clusters}} in Nestimate 0.4.3. This thin
#' wrapper is preserved so the function name in older tutorials and the
#' historical pkgdown reference continues to work; it issues a one-shot
#' deprecation warning and forwards every argument unchanged.
#'
#' @param ... Passed verbatim to \code{\link{build_clusters}}.
#' @return The \code{net_clustering} object returned by
#'   \code{\link{build_clusters}}.
#' @seealso \code{\link{build_clusters}}, \code{\link{cluster_network}},
#'   \code{\link{cluster_mmm}}.
#' @keywords internal
#' @export
cluster_data <- function(...) {
  .Deprecated("build_clusters",
              package = "Nestimate",
              msg = paste0(
                "cluster_data() was renamed to build_clusters() in ",
                "Nestimate 0.4.3. Please use build_clusters() instead."))
  build_clusters(...)
}


# ==============================================================================
# 6. S3 Methods
# ==============================================================================

# ------------------------------------------------------------------------------
# Shared formatter for cluster summary tables.
#
# `cols` is a named list of vectors (length k each). The first column drives
# row count; all other columns must match that length. Returns a character
# vector of formatted lines, indented and headed by a separator rule. Used by
# print.net_clustering / print.net_mmm / print.net_mmm_clustering /
# print.netobject_group so a netobject_group built by cluster_network() and
# one built by cluster_mmm() print with the same shape.
# ------------------------------------------------------------------------------
.cluster_table_lines <- function(cols, indent = "  ", sep = "  ") {
  stopifnot(is.list(cols), length(cols) >= 1L, !is.null(names(cols)))
  k <- length(cols[[1L]])
  if (k == 0L) return(character(0L))

  names_v <- names(cols)
  # Column widths from header + cell content (cells are already character).
  widths <- vapply(seq_along(cols), function(j) {
    cells <- as.character(cols[[j]])
    max(nchar(names_v[[j]]), if (length(cells)) max(nchar(cells)) else 0L)
  }, integer(1L))

  ncol_v <- length(cols)
  # `print(data.frame)` left-aligns; we keep the same convention but pad
  # only up to (not past) the *last* column so trailing whitespace on the
  # right edge is gone.
  pad <- function(s, w, last) if (last) s else formatC(s, width = w, flag = "-")

  build_row <- function(cells) {
    parts <- vapply(seq_len(ncol_v), function(j) {
      pad(cells[[j]], widths[[j]], last = j == ncol_v)
    }, character(1L))
    paste(parts, collapse = sep)
  }

  header <- build_row(as.list(names_v))
  rows   <- vapply(seq_len(k), function(i) {
    build_row(lapply(cols, function(col) as.character(col[[i]])))
  }, character(1L))

  paste0(indent, c(header, rows))
}

# Helper: format size + percentage as a "30 (50.0%)" string.
.fmt_size_pct <- function(sizes, total = sum(sizes)) {
  pct <- if (isTRUE(total > 0)) sizes / total * 100 else rep(0, length(sizes))
  sprintf("%d (%4.1f%%)", as.integer(sizes), pct)
}

# Per-cluster mean within-cluster distance. Singletons return 0 (no
# within-distance defined). Used by print.net_clustering,
# summary.net_clustering, cluster_diagnostics.net_clustering, and
# cluster_choice's per-row builder.
.per_cluster_within_dist <- function(distance, assignments, k = max(assignments)) {
  if (is.null(distance)) return(rep(NA_real_, k))
  dmat <- as.matrix(distance)
  vapply(seq_len(k), function(cl) {
    members <- which(assignments == cl)
    if (length(members) > 1L) {
      sub <- dmat[members, members]
      mean(sub[lower.tri(sub)])
    } else {
      0
    }
  }, numeric(1L))
}

#' Print Method for net_clustering
#'
#' Compact, fixed-width summary of a sequence-clustering result. The header
#' carries the clustering method and dissimilarity; the per-cluster table
#' carries cluster size (count and percentage) and mean within-cluster
#' distance when available. Optional medoid and covariate lines surface
#' only when those fields are populated.
#'
#' @param x A \code{net_clustering} object.
#' @param digits Integer. Decimal places used for floating-point statistics
#'   in the printout. Default \code{3}. Non-breaking: existing
#'   \code{print(x)} calls keep their previous formatting.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' cl <- build_clusters(seqs, k = 2)
#' print(cl)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 20, TRUE),
#'   V2 = sample(c("A","B","C"), 20, TRUE),
#'   V3 = sample(c("A","B","C"), 20, TRUE)
#' )
#' cl <- build_clusters(seqs, k = 2)
#' print(cl)
#' }
#'
#' @export
print.net_clustering <- function(x, digits = 3L, ...) {
  .net_clustering_check_unused_dots("print.net_clustering", ...)
  .net_clustering_check_digits(digits)
  digits <- as.integer(digits)
  k <- as.integer(x$k)
  sizes <- as.integer(x$sizes)
  n_total <- sum(sizes)

  diss_lab <- if (isTRUE(x$weighted)) {
    sprintf("%s (weighted, lambda = %g)", x$dissimilarity, x$lambda)
  } else {
    x$dissimilarity
  }

  cat(sprintf("Sequence Clustering [%s]\n", x$method))
  cat(sprintf("  Sequences: %d  |  Clusters: %d\n", n_total, k))
  cat(sprintf("  Dissimilarity: %s\n", diss_lab))
  if (!is.null(x$silhouette)) {
    cat(sprintf("  Quality: silhouette = %.*f\n",
                digits, x$silhouette))
  }

  # Per-cluster mean within-cluster distance, when distance matrix is present.
  mean_within <- if (!is.null(x$distance) && !is.null(x$assignments) && k > 0L) {
    .per_cluster_within_dist(x$distance, x$assignments, k)
  } else {
    rep(NA_real_, k)
  }

  cat("\n")
  cols <- list(
    Cluster = sprintf("%d", seq_len(k)),
    N       = .fmt_size_pct(sizes, n_total)
  )
  if (!all(is.na(mean_within))) {
    cols[["Mean within-dist"]] <- sprintf(paste0("%.", digits, "f"),
                                          mean_within)
  }
  if (!is.null(x$medoids) && length(x$medoids) == k) {
    cols[["Medoid"]] <- as.character(x$medoids)
  }
  cat(paste(.cluster_table_lines(cols), collapse = "\n"), "\n", sep = "")

  if (!is.null(x$covariates)) {
    cov_names <- setdiff(
      unique(x$covariates$coefficients$variable), "(Intercept)"
    )
    cat(sprintf("\n  Covariates: %s (post-hoc, %d predictors)\n",
                paste(cov_names, collapse = ", "), length(cov_names)))
  }

  invisible(x)
}

#' Summary Method for net_clustering
#'
#' @param object A \code{net_clustering} object.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' cl <- build_clusters(seqs, k = 2)
#' summary(cl)
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 20, TRUE),
#'   V2 = sample(c("A","B","C"), 20, TRUE),
#'   V3 = sample(c("A","B","C"), 20, TRUE)
#' )
#' cl <- build_clusters(seqs, k = 2)
#' summary(cl)
#' }
#'
#' @export
summary.net_clustering <- function(object, ...) {
  .net_clustering_check_unused_dots("summary.net_clustering", ...)
  k <- object$k
  assignments <- object$assignments
  sizes <- as.integer(tabulate(assignments, nbins = k))
  cluster_stats <- data.frame(
    cluster          = seq_len(k),
    size             = sizes,
    mean_within_dist = .per_cluster_within_dist(object$distance,
                                                 assignments, k),
    stringsAsFactors = FALSE
  )
  rownames(cluster_stats) <- NULL

  cat("Sequence Clustering Summary\n")
  cat("  Method:       ", object$method, "\n")
  cat("  Dissimilarity:", object$dissimilarity, "\n")
  cat("  Silhouette:   ", round(object$silhouette, 4), "\n\n")
  cat("Per-cluster statistics:\n")
  print(cluster_stats, row.names = FALSE)

  if (!is.null(object$covariates)) {
    cat("\n")
    .print_covariate_profiles(object$covariates)
    cat("\nNote: Covariates are post-hoc and do not influence cluster assignments.\n")
    # Return a tidy data.frame of the inferential output. Cluster sizes,
    # fit stats, residuals, and the full untidied list are attached as
    # attributes so power users still have access without spelunking
    # through `$covariates`.
    return(invisible(.tidy_covariates(object$covariates,
                                       cluster_stats = cluster_stats)))
  }

  cluster_stats
}

#' Plot Sequence Clustering Results
#'
#' @param x A \code{net_clustering} object.
#' @param type Character. Plot type: \code{"silhouette"} (per-observation
#'   silhouette bars), \code{"mds"} (2D MDS projection), or
#'   \code{"heatmap"} (distance matrix heatmap ordered by cluster).
#'   Default: \code{"silhouette"}.
#' @param ... Unsupported. Supplying unused arguments raises an error.
#' @return A \code{ggplot} object (invisibly).
#'
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' cl <- build_clusters(seqs, k = 2)
#' plot(cl, type = "silhouette")
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(c("A","B","C"), 20, TRUE),
#'   V2 = sample(c("A","B","C"), 20, TRUE),
#'   V3 = sample(c("A","B","C"), 20, TRUE)
#' )
#' cl <- build_clusters(seqs, k = 2)
#' plot(cl, type = "silhouette")
#' }
#'
#' @import ggplot2
#' @export
plot.net_clustering <- function(x, type = c("silhouette", "mds", "heatmap",
                                             "predictors"),
                                 combined = TRUE, ...) {
  .net_clustering_check_unused_dots("plot.net_clustering", ...)
  type <- match.arg(type)
  stopifnot(is.logical(combined), length(combined) == 1L)

  if (type == "predictors") {
    if (is.null(x$covariates)) {
      stop("No covariate analysis found. Run build_clusters() with covariates.",
           call. = FALSE)
    }
    if (is.null(x$covariates$coefficients)) {
      stop("plot(..., type = \"predictors\") requires odds-ratio coefficients ",
           "from estimator = \"multinom\". The current covariate analysis ",
           "(estimator = \"chisq\") has no coefficients; inspect ",
           "x$covariates$tests and x$covariates$residuals instead.",
           call. = FALSE)
    }
    return(.plot_covariate_forest(
      x$covariates$coefficients,
      sprintf("Covariate Effects on Cluster Membership (ref: Cluster %s)",
              x$covariates$fit$reference_cluster),
      combined = combined
    ))
  }

  if (type == "silhouette") {
    sil <- cluster::silhouette(x$assignments, dist = x$distance)
    sil_df <- data.frame(
      obs = seq_len(nrow(sil)),
      cluster = factor(sil[, 1L]),
      width = sil[, 3L]
    )
    # Order by cluster then width within cluster
    sil_df <- sil_df[order(sil_df$cluster, -sil_df$width), ]
    sil_df$rank <- seq_len(nrow(sil_df))

    p <- ggplot(sil_df, aes(x = .data$rank, y = .data$width,
                             fill = .data$cluster)) +
      geom_col(width = 1) +
      coord_flip() +
      labs(x = "Observation", y = "Silhouette Width",
           title = sprintf("Silhouette Plot (avg = %.3f)", x$silhouette),
           fill = "Cluster") +
      theme_minimal(base_size = 12) +
      geom_hline(yintercept = x$silhouette, linetype = "dashed", alpha = 0.5)

  } else if (type == "mds") {
    mds <- stats::cmdscale(x$distance, k = 2L)
    mds_df <- data.frame(
      dim1 = mds[, 1L],
      dim2 = mds[, 2L],
      cluster = factor(x$assignments)
    )

    p <- ggplot(mds_df, aes(x = .data$dim1, y = .data$dim2,
                             colour = .data$cluster)) +
      geom_point(size = 2.5) +
      labs(x = "MDS Dimension 1", y = "MDS Dimension 2",
           title = sprintf("MDS Projection (k = %d)", x$k),
           colour = "Cluster") +
      theme_minimal(base_size = 12)

  } else {
    # Heatmap — order by cluster assignment
    ord <- order(x$assignments)
    dist_mat <- as.matrix(x$distance)[ord, ord]
    n <- nrow(dist_mat)

    # Build long-format data frame
    idx <- expand.grid(row = seq_len(n), col = seq_len(n))
    idx$distance <- as.vector(dist_mat)
    idx$row <- factor(idx$row, levels = seq_len(n))
    idx$col <- factor(idx$col, levels = rev(seq_len(n)))

    p <- ggplot(idx, aes(x = .data$row, y = .data$col,
                          fill = .data$distance)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      labs(x = "Sequence", y = "Sequence",
           title = "Distance Matrix (ordered by cluster)",
           fill = "Distance") +
      theme_minimal(base_size = 12) +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  }

  print(p)
  invisible(p)
}

.net_clustering_check_unused_dots <- function(method, ...) {
  dots <- list(...)
  if (!length(dots)) {
    return(invisible(TRUE))
  }
  dot_names <- names(dots)
  dot_names[!nzchar(dot_names)] <- paste0("..", which(!nzchar(dot_names)))
  stop(
    method, "() got unsupported argument",
    if (length(dots) == 1L) ": " else "s: ",
    paste(dot_names, collapse = ", "),
    call. = FALSE
  )
}

.net_clustering_check_digits <- function(digits) {
  if (!is.numeric(digits) || length(digits) != 1L || !is.finite(digits) ||
      digits != floor(digits) || digits < 0L) {
    stop("'digits' must be a single non-negative whole finite number.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# ==============================================================================
# 7. Covariate Analysis (post-hoc)
# ==============================================================================

#' Resolve covariates into a data.frame
#'
#' Handles five input forms: formula, character vector, string, data.frame,
#' NULL. Looks up column names in metadata or raw data as needed.
#' @param covariates User-supplied covariates argument
#' @param raw_data Original input to build_clusters (before extraction)
#' @param n Number of sequences
#' @return A list with \code{cov_df} (data.frame) and \code{formula} (formula),
#'   or NULL if covariates is NULL.
#' @noRd
.resolve_covariates <- function(covariates, raw_data, n) {
  if (is.null(covariates)) return(NULL)

  if (!requireNamespace("nnet", quietly = TRUE)) {
    stop("The 'nnet' package is required for covariate analysis. ", # nocov
         "Install it with install.packages('nnet').", call. = FALSE) # nocov
  }

  # --- Covariates supplied as a data.frame: use all columns ---
  if (is.data.frame(covariates)) {
    if (nrow(covariates) != n) {
      stop(sprintf(
        "Covariate data.frame has %d rows but sequence data has %d rows.",
        nrow(covariates), n
      ), call. = FALSE)
    }
    rhs <- paste(names(covariates), collapse = " + ")
    return(list(cov_df = covariates, rhs = rhs))
  }

  # --- Parse into column names ---
  if (inherits(covariates, "formula")) {
    rhs_text <- deparse(covariates[[length(covariates)]])
    cov_names <- all.vars(covariates)
  } else if (is.character(covariates) && length(covariates) == 1L &&
             grepl("+", covariates, fixed = TRUE)) {
    # String like "Age + Gender"
    rhs_text <- covariates
    cov_names <- all.vars(stats::as.formula(paste("~", covariates)))
  } else if (is.character(covariates)) {
    # Character vector like c("Age", "Gender")
    cov_names <- covariates
    rhs_text <- paste(cov_names, collapse = " + ")
  } else {
    stop("covariates must be a formula, character vector, string, or data.frame.",
         call. = FALSE)
  }

  if (length(cov_names) == 0L) {
    stop("No covariate variables specified.", call. = FALSE)
  }

  # --- Look up columns in metadata or raw data ---
  cov_df <- NULL

  # Try netobject / cograph_network: walk $metadata first, then $data columns.
  # Both classes carry the same slot shape, so we resolve from either.
  if (inherits(raw_data, "netobject") ||
      inherits(raw_data, "cograph_network")) {
    avail_df <- NULL
    avail_names <- character(0L)
    if (!is.null(raw_data$metadata) && is.data.frame(raw_data$metadata)) {
      avail_df <- raw_data$metadata
      avail_names <- names(raw_data$metadata)
    }
    if (!is.null(raw_data$data) && is.data.frame(raw_data$data)) {
      extra <- setdiff(names(raw_data$data), avail_names)
      if (length(extra) > 0L) {
        avail_names <- c(avail_names, extra)
        avail_df <- if (is.null(avail_df)) {
          raw_data$data[, extra, drop = FALSE]
        } else {
          cbind(avail_df, raw_data$data[, extra, drop = FALSE])
        }
      }
    }
    if (length(avail_names) == 0L) {
      stop("Network object has no $metadata or $data slot to resolve ",
           "column-name covariates against. Pass covariates as a ",
           "data.frame instead.", call. = FALSE)
    }
    if (all(cov_names %in% avail_names)) {
      cov_df <- avail_df[, cov_names, drop = FALSE]
    } else {
      missing <- setdiff(cov_names, avail_names)
      stop(sprintf(
        "Covariate columns not found in network object: %s. Available: %s",
        paste(missing, collapse = ", "),
        paste(avail_names, collapse = ", ")
      ), call. = FALSE)
    }
  }

  # tna objects: require data.frame covariates (no canonical metadata slot).
  if (is.null(cov_df) && inherits(raw_data, "tna")) {
    stop("Column-name covariates are not supported for tna inputs. ",
         "Pass covariates as a data.frame instead.", call. = FALSE)
  }

  # Try raw data.frame columns
  if (is.null(cov_df) && (is.data.frame(raw_data) || is.matrix(raw_data))) {
    avail <- colnames(raw_data)
    if (all(cov_names %in% avail)) {
      cov_df <- as.data.frame(raw_data[, cov_names, drop = FALSE],
                               stringsAsFactors = FALSE)
    } else {
      missing <- setdiff(cov_names, avail)
      stop(sprintf(
        "Covariate columns not found: %s. Available: %s",
        paste(missing, collapse = ", "),
        paste(avail, collapse = ", ")
      ), call. = FALSE)
    }
  }

  if (is.null(cov_df)) {
    stop("Could not resolve covariate columns from the input data.", # nocov start
         call. = FALSE) # nocov end
  }

  if (nrow(cov_df) != n) { # nocov start
    stop(sprintf(
      "Covariate data has %d rows but sequence data has %d rows.",
      nrow(cov_df), n
    ), call. = FALSE) # nocov end
  }

  list(cov_df = cov_df, rhs = rhs_text)
}

#' Compute per-cluster descriptive profiles
#' @noRd
.compute_cluster_profiles <- function(cov_df, assignments, k) {
  n_total <- length(assignments)
  num_cols <- vapply(cov_df, is.numeric, logical(1L))
  cat_cols <- !num_cols

  # --- Numeric profiles ---
  numeric_profile <- NULL
  if (any(num_cols)) {
    num_names <- names(cov_df)[num_cols]
    rows <- vector("list", k * length(num_names))
    idx <- 0L
    for (cl in seq_len(k)) {
      mask <- assignments == cl
      cl_n <- sum(mask)
      cl_pct <- round(100 * cl_n / n_total, 1)
      for (vname in num_names) {
        idx <- idx + 1L
        vals <- cov_df[[vname]][mask]
        vals <- vals[!is.na(vals)]
        rows[[idx]] <- data.frame(
          cluster = cl, n = cl_n, pct = cl_pct, variable = vname,
          mean = mean(vals), sd = stats::sd(vals),
          median = stats::median(vals),
          stringsAsFactors = FALSE, row.names = NULL
        )
      }
    }
    numeric_profile <- do.call(rbind, rows)
  }

  # --- Categorical profiles ---
  categorical_profile <- NULL
  if (any(cat_cols)) {
    cat_names <- names(cov_df)[cat_cols]
    # Pre-count total rows to avoid O(n^2) list growth
    total_rows <- 0L
    for (cl in seq_len(k)) {
      mask <- assignments == cl
      for (vname in cat_names) {
        total_rows <- total_rows + length(unique(cov_df[[vname]][mask & !is.na(cov_df[[vname]])]))
      }
    }
    rows <- vector("list", total_rows)
    idx <- 0L
    for (cl in seq_len(k)) {
      mask <- assignments == cl
      cl_n <- sum(mask)
      for (vname in cat_names) {
        tab <- table(cov_df[[vname]][mask], useNA = "no")
        for (lev in names(tab)) {
          idx <- idx + 1L
          rows[[idx]] <- data.frame(
            cluster = cl, n = cl_n, variable = vname,
            level = lev, count = as.integer(tab[[lev]]),
            pct = round(100 * as.integer(tab[[lev]]) / cl_n, 1),
            stringsAsFactors = FALSE, row.names = NULL
          )
        }
      }
    }
    categorical_profile <- do.call(rbind, rows)
  }

  list(numeric = numeric_profile, categorical = categorical_profile)
}

#' Run covariate analysis on cluster assignments
#'
#' Dispatches on `estimator`:
#' \itemize{
#'   \item `"firth"` (default) — Firth's penalised multinomial logit via
#'     `brglm2::brmultinom`. Bias-reduced; estimates stay finite under
#'     quasi-complete separation. AIC/BIC/McFadden derived from the same
#'     log-likelihood scale as the ML fit, so they remain comparable.
#'   \item `"multinom"` — classical ML multinomial via `nnet::multinom`.
#'     Identical to the pre-2026-05-10 behaviour. A `warning()` is emitted
#'     because rare-cell pathologies (separation, infinite ORs) are silent
#'     in this estimator.
#'   \item `"chisq"` — descriptive only, WeightedCluster-style: per-covariate
#'     test of cluster association (chi-square + Cramer's V for factors,
#'     Kruskal-Wallis + eta² for numerics). No ORs; a separate residuals
#'     table reports per-cell standardized adjusted residuals so the user
#'     can read which (cluster, level) cells drive the association.
#' }
#'
#' Output shape:
#' \itemize{
#'   \item firth/multinom: `list(profiles, coefficients, fit, model)` —
#'     `coefficients` carries the OR table, `fit` carries AIC/BIC/McFadden.
#'   \item chisq: `list(profiles, tests, residuals, fit)` — `tests` is a
#'     per-covariate row, `residuals` is a tidy frame of standardized
#'     adjusted residuals (factor covariates only).
#' }
#'
#' @noRd
.run_covariate_analysis <- function(assignments, cov_df, rhs, k,
                                     estimator = c("firth", "multinom",
                                                   "chisq")) {
  estimator <- match.arg(estimator)

  # --- Prepare data (shared across all three estimators) ---
  fit_df <- cov_df
  for (nm in names(fit_df)) {
    if (is.character(fit_df[[nm]])) {
      fit_df[[nm]] <- factor(fit_df[[nm]])
    } else if (is.ordered(fit_df[[nm]])) {
      fit_df[[nm]] <- factor(fit_df[[nm]], ordered = FALSE)
    }
  }
  complete <- stats::complete.cases(fit_df)
  n_dropped <- sum(!complete)
  if (n_dropped > 0L) {
    warning(sprintf(
      "Dropped %d rows with NA covariates (%d remaining).",
      n_dropped, sum(complete)
    ), call. = FALSE)
    fit_df <- fit_df[complete, , drop = FALSE]
    assignments <- assignments[complete]
  }
  fit_df$cluster <- factor(assignments)
  for (nm in setdiff(names(fit_df), "cluster")) {
    if (length(unique(fit_df[[nm]][!is.na(fit_df[[nm]])])) < 2L) {
      stop(sprintf("Covariate '%s' is constant. Remove it.", nm), call. = FALSE)
    }
  }
  profiles_df <- if (n_dropped > 0L) cov_df[complete, , drop = FALSE] else cov_df
  profiles <- .compute_cluster_profiles(profiles_df, assignments, k)

  if (estimator == "chisq") {
    return(.covariate_chisq(fit_df, profiles))
  }

  .covariate_logit(fit_df, rhs, profiles, estimator,
                    raw_assignments = assignments)
}


# Factor + numeric per-covariate descriptive tests (WeightedCluster style).
# No multinomial fit; no ORs. Produces a per-covariate table + a tidy
# standardized-residuals frame for the factor covariates so users can
# read which (cluster, level) cells drive each association.
#' @noRd
.covariate_chisq <- function(fit_df, profiles) {
  cov_names <- setdiff(names(fit_df), "cluster")
  cluster <- fit_df$cluster

  test_rows <- list()
  resid_rows <- list()

  for (nm in cov_names) {
    x <- fit_df[[nm]]
    if (is.factor(x)) {
      tab <- table(cluster = cluster, level = x)
      # Suppress chi-sq's small-expected-count warning; we report the
      # statistic regardless and let users read the residuals.
      ct <- suppressWarnings(stats::chisq.test(tab))
      n_total <- sum(tab)
      n_cells <- prod(dim(tab))
      # Cramer's V on a (r x c) table: sqrt(chi^2 / (n * min(r-1, c-1))).
      v <- sqrt(as.numeric(ct$statistic) /
                  (n_total * min(nrow(tab) - 1L, ncol(tab) - 1L)))
      test_rows[[length(test_rows) + 1L]] <- data.frame(
        variable = nm, type = "factor",
        statistic = unname(ct$statistic),
        df = unname(ct$parameter),
        p_value = ct$p.value,
        effect_size = v,
        effect_label = "Cramer's V",
        stringsAsFactors = FALSE
      )
      # stdres is the residual divided by its standard error -> approx N(0,1)
      # under H_0; |stdres| > 2 flags a cell with significant departure.
      sr <- suppressWarnings(ct$stdres)
      sr_df <- as.data.frame(as.table(sr), stringsAsFactors = FALSE)
      names(sr_df) <- c("cluster", "level", "std_residual")
      sr_df$variable <- nm
      sr_df$observed <- as.numeric(as.vector(tab))
      sr_df$expected <- as.numeric(as.vector(ct$expected))
      resid_rows[[length(resid_rows) + 1L]] <-
        sr_df[, c("variable", "cluster", "level", "observed",
                  "expected", "std_residual"),
              drop = FALSE]
    } else if (is.numeric(x)) {
      kw <- suppressWarnings(stats::kruskal.test(x ~ cluster))
      # eta-squared for KW: H / (n - 1) bounded to [0, 1]; sometimes
      # reported as (H - g + 1) / (n - g). Use the simpler, more common
      # form so it matches `effectsize::rank_epsilon_squared`'s scale.
      n_total <- length(x)
      eta_sq <- as.numeric(kw$statistic) / (n_total - 1L)
      eta_sq <- max(0, min(1, eta_sq))
      test_rows[[length(test_rows) + 1L]] <- data.frame(
        variable = nm, type = "numeric",
        statistic = unname(kw$statistic),
        df = unname(kw$parameter),
        p_value = kw$p.value,
        effect_size = eta_sq,
        effect_label = "eta^2 (KW)",
        stringsAsFactors = FALSE
      )
    } else {
      stop(sprintf(
        "Covariate '%s' has unsupported type '%s' for estimator = 'chisq'. ",
        nm, class(x)[1L]),
        "Convert to factor or numeric.", call. = FALSE)
    }
  }

  tests <- do.call(rbind, test_rows)
  if (!is.null(tests)) {
    tests$sig <- ifelse(is.na(tests$p_value), "",
      ifelse(tests$p_value < 0.001, "***",
      ifelse(tests$p_value < 0.01,  "**",
      ifelse(tests$p_value < 0.05,  "*",
      ifelse(tests$p_value < 0.1,   ".", "")))))
  }
  residuals_df <- if (length(resid_rows)) do.call(rbind, resid_rows) else NULL

  list(
    profiles = profiles,
    tests = tests,
    residuals = residuals_df,
    fit = list(method = "chisq",
               n_factor_tests = sum(vapply(cov_names,
                 function(n) is.factor(fit_df[[n]]), logical(1L))),
               n_numeric_tests = sum(vapply(cov_names,
                 function(n) is.numeric(fit_df[[n]]), logical(1L)))),
    model = NULL
  )
}


# Multinomial logit, dispatched on estimator. Firth via brglm2::brmultinom
# (default); ML via nnet::multinom. Output shape is shared so downstream
# print/summary code does not need to branch.
#' @noRd
.covariate_logit <- function(fit_df, rhs, profiles, estimator,
                              raw_assignments) {
  if (estimator == "multinom") {
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' is required for estimator = 'multinom'. ",
           "Install with: install.packages('nnet')", call. = FALSE)
    }
    warning(
      "estimator = 'multinom' uses nnet::multinom (classical ML). ",
      "Coefficient estimates can diverge under quasi-complete separation ",
      "(rare-cell or zero-cell covariate levels), producing astronomical ",
      "ORs with degenerate confidence intervals. The default ",
      "estimator = 'firth' applies Firth's penalty and is robust to this. ",
      "See Heinze & Schemper (2002) and brglm2::brmultinom.",
      call. = FALSE
    )
    fml <- stats::as.formula(paste("cluster ~", rhs))
    fit <- nnet::multinom(fml, data = fit_df, trace = FALSE)
    null_fit <- nnet::multinom(cluster ~ 1, data = fit_df, trace = FALSE)
  } else {
    # estimator == "firth"
    if (!requireNamespace("brglm2", quietly = TRUE)) {
      stop("Package 'brglm2' is required for estimator = 'firth'. ",
           "Install with: install.packages('brglm2'). ",
           "Pass estimator = 'multinom' to use nnet::multinom instead, ",
           "or estimator = 'chisq' for descriptive-only output.",
           call. = FALSE)
    }
    fml <- stats::as.formula(paste("cluster ~", rhs))
    fit <- suppressWarnings(brglm2::brmultinom(fml, data = fit_df,
                                                ref = 1L))
    null_fit <- suppressWarnings(brglm2::brmultinom(cluster ~ 1,
                                                     data = fit_df,
                                                     ref = 1L))
  }

  s <- summary(fit)
  coefs <- s$coefficients
  ses <- s$standard.errors
  if (is.null(ses) && !is.null(s$standard.errors)) ses <- s$standard.errors
  if (is.null(ses)) {
    # brmultinom stores SE on the vcov diagonal
    ses <- matrix(sqrt(diag(stats::vcov(fit))),
                  nrow = nrow(coefs), byrow = FALSE,
                  dimnames = dimnames(coefs))
  }
  if (is.null(dim(coefs))) {
    coefs <- matrix(coefs, nrow = 1L,
                    dimnames = list(levels(fit_df$cluster)[2L],
                                     names(coefs)))
    ses <- matrix(ses, nrow = 1L,
                  dimnames = list(levels(fit_df$cluster)[2L],
                                   names(ses)))
  }
  n_params <- length(coefs)
  min_cl <- min(table(raw_assignments))
  if (min_cl < n_params) {
    warning(sprintf(
      "Smallest cluster has %d observations but model has %d parameters. ",
      min_cl, n_params),
      "Estimates may be unreliable.", call. = FALSE)
  }

  rows <- vector("list", nrow(coefs) * ncol(coefs))
  idx <- 0L
  for (i in seq_len(nrow(coefs))) {
    for (j in seq_len(ncol(coefs))) {
      idx <- idx + 1L
      est <- coefs[i, j]; se <- ses[i, j]
      z_val <- est / se
      p_val <- 2 * (1 - stats::pnorm(abs(z_val)))
      or <- exp(est); ci_lo <- exp(est - 1.96 * se); ci_hi <- exp(est + 1.96 * se)
      sig <- if (is.na(p_val)) "" else if (p_val < 0.001) "***" else
        if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else
        if (p_val < 0.1) "." else ""
      rows[[idx]] <- data.frame(
        cluster = rownames(coefs)[i],
        variable = colnames(coefs)[j],
        estimate = est, std_error = se,
        odds_ratio = or, ci_lower = ci_lo, ci_upper = ci_hi,
        z = z_val, p = p_val, sig = sig,
        stringsAsFactors = FALSE, row.names = NULL
      )
    }
  }
  coef_df <- do.call(rbind, rows)

  null_dev <- tryCatch(-2 * as.numeric(stats::logLik(null_fit)),
                       error = function(e) NA_real_)
  model_dev <- tryCatch(-2 * as.numeric(stats::logLik(fit)),
                        error = function(e) NA_real_)
  mcfadden_r2 <- if (is.na(null_dev) || is.na(model_dev)) NA_real_
                 else 1 - model_dev / null_dev

  list(
    profiles = profiles,
    coefficients = coef_df,
    fit = list(
      method = estimator,
      aic = tryCatch(stats::AIC(fit), error = function(e) NA_real_),
      bic = tryCatch(stats::BIC(fit), error = function(e) NA_real_),
      deviance = model_dev,
      mcfadden_r2 = mcfadden_r2,
      reference_cluster = levels(fit_df$cluster)[1L]
    ),
    model = fit
  )
}

#' Print covariate profiles and coefficient table (shared by build_clusters/mmm)
#' @noRd
.print_covariate_profiles <- function(cov, header = "Post-hoc Covariate Analysis (does not influence cluster membership)") {
  cat(header, "\n\n")

  # Numeric profiles
  if (!is.null(cov$profiles$numeric)) {
    cat("Cluster Profiles (numeric):\n")
    np <- cov$profiles$numeric
    vars <- unique(np$variable)
    rows <- lapply(unique(np$cluster), function(cl) {
      sub <- np[np$cluster == cl, , drop = FALSE]
      row <- data.frame(
        Cluster = cl,
        `N (%)` = sprintf("%d (%.0f%%)", sub$n[1L], sub$pct[1L]),
        stringsAsFactors = FALSE, check.names = FALSE
      )
      for (v in vars) {
        sv <- sub[sub$variable == v, ]
        row[[sprintf("%s Mean (SD)", v)]] <- sprintf("%.2f (%.2f)", sv$mean, sv$sd)
        row[[sprintf("%s Median", v)]] <- sprintf("%.2f", sv$median)
      }
      row
    })
    print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
    cat("\n")
  }

  # Categorical profiles
  if (!is.null(cov$profiles$categorical)) {
    cat("Cluster Profiles (categorical):\n")
    cp <- cov$profiles$categorical
    # Pre-compute all variable=level combos for consistent columns
    all_cols <- unique(paste0(cp$variable, "=", cp$level))
    rows <- lapply(unique(cp$cluster), function(cl) {
      sub <- cp[cp$cluster == cl, , drop = FALSE]
      row <- data.frame(
        Cluster = cl, N = sub$n[1L],
        stringsAsFactors = FALSE, check.names = FALSE
      )
      for (col in all_cols) {
        parts <- strsplit(col, "=", fixed = TRUE)[[1L]]
        v <- parts[1L]
        lev <- paste(parts[-1L], collapse = "=")  # handle = in level names
        slev <- sub[sub$variable == v & sub$level == lev, ]
        if (nrow(slev) > 0L) {
          row[[sprintf("%s N(%%)", col)]] <- sprintf(
            "%d (%.0f%%)", slev$count, slev$pct)
        } else {
          row[[sprintf("%s N(%%)", col)]] <- "0 (0%)"
        }
      }
      row
    })
    print(do.call(rbind, rows), row.names = FALSE, right = FALSE)
    cat("\n")
  }

  # Estimator-specific section: ORs (firth/multinom) or test table (chisq).
  if (identical(cov$fit$method, "chisq")) {
    cat("Per-Covariate Cluster Association Tests:\n")
    tt <- cov$tests
    disp <- data.frame(
      Variable = tt$variable,
      Type = tt$type,
      Statistic = sprintf("%.3f", tt$statistic),
      df = ifelse(is.na(tt$df), "", as.character(tt$df)),
      p = ifelse(tt$p_value < 0.001, "<0.001",
                 sprintf("%.3f", tt$p_value)),
      Effect = paste0(sprintf("%.3f", tt$effect_size), " (", tt$effect_label, ")"),
      Sig = tt$sig,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    print(disp, row.names = FALSE, right = FALSE)

    if (!is.null(cov$residuals) && nrow(cov$residuals) > 0L) {
      cat("\nStandardized Adjusted Residuals (|r| > 2 flags significant cells):\n")
      rs <- cov$residuals
      rs$std_residual <- sprintf("%.2f", rs$std_residual)
      rs$observed <- as.integer(rs$observed)
      rs$expected <- sprintf("%.1f", rs$expected)
      print(rs, row.names = FALSE, right = FALSE)
    }
    cat("\nMethod: chi-square (categorical) / Kruskal-Wallis (numeric); ",
        "descriptive only.\n", sep = "")
  } else {
    cat(sprintf("Predictors of Membership (estimator = %s, reference: Cluster %s):\n",
                cov$fit$method %||% "multinom",
                cov$fit$reference_cluster))
    coef_display <- cov$coefficients[
      cov$coefficients$variable != "(Intercept)", , drop = FALSE]
    disp <- data.frame(
      Cluster = coef_display$cluster,
      Variable = coef_display$variable,
      OR = sprintf("%.2f", coef_display$odds_ratio),
      `95% CI` = sprintf("[%.2f, %.2f]", coef_display$ci_lower,
                          coef_display$ci_upper),
      p = ifelse(coef_display$p < 0.001, "<0.001",
                 sprintf("%.3f", coef_display$p)),
      Sig = coef_display$sig,
      stringsAsFactors = FALSE, check.names = FALSE
    )
    print(disp, row.names = FALSE, right = FALSE)
    cat(sprintf("\nModel: AIC = %.1f | BIC = %.1f | McFadden R-squared = %.2f\n",
                cov$fit$aic, cov$fit$bic, cov$fit$mcfadden_r2))
  }
}

# Tidy a covariate-analysis result list into a single data.frame the user
# can hand to write.csv / kable / dplyr without further reshaping. Returns
# different shapes per estimator (logit vs chisq) — the column set names
# the relevant statistics for each path so the meaning is unambiguous.
#
# Attributes carry the metadata that didn't fit a tidy column (cluster
# sizes, fit stats, residuals frame for chisq). `print.tidy_covariates`
# renders the same human-readable view that `.print_covariate_profiles()`
# does, so the tidy frame remains REPL-friendly.
#' @noRd
.tidy_covariates <- function(cov, cluster_stats = NULL) {
  if (is.null(cov)) return(NULL)
  method <- cov$fit$method %||% "multinom"

  if (method == "chisq") {
    out <- cov$tests
    if (!is.null(out)) {
      out <- out[, c("variable", "type", "statistic", "df", "p_value",
                     "effect_size", "effect_label", "sig"),
                  drop = FALSE]
      names(out)[names(out) == "p_value"] <- "p"
      names(out)[names(out) == "effect_label"] <- "effect_kind"
      names(out)[names(out) == "effect_size"] <- "effect"
      rownames(out) <- NULL
    }
  } else {
    co <- cov$coefficients
    co <- co[co$variable != "(Intercept)", , drop = FALSE]
    out <- data.frame(
      cluster   = co$cluster,
      predictor = co$variable,
      OR        = co$odds_ratio,
      CI_low    = co$ci_lower,
      CI_high   = co$ci_upper,
      p         = co$p,
      sig       = co$sig,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  attr(out, "estimator")       <- method
  attr(out, "fit")             <- cov$fit
  attr(out, "profiles")        <- cov$profiles
  attr(out, "residuals")       <- cov$residuals
  attr(out, "cluster_stats")   <- cluster_stats
  attr(out, "covariates_full") <- cov
  class(out) <- c("tidy_covariates", "data.frame")
  out
}

#' Print method for tidy covariate output
#'
#' Prints a one-line header naming the estimator, then the data.frame.
#' The full human-readable view (per-cluster stats, profiles, OR/test
#' tables) was already printed by `summary()` when this object was
#' produced, so this method intentionally stays minimal to avoid
#' duplication. Auto-prints when the user types the variable at the REPL.
#'
#' @param x A `tidy_covariates` data.frame.
#' @param ... Ignored.
#' @return The input invisibly.
#' @export
print.tidy_covariates <- function(x, ...) {
  cat(sprintf("# Tidy covariate analysis (estimator = %s)\n",
              attr(x, "estimator") %||% "unknown"))
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}


#' Odds ratio forest plot (shared by build_clusters/mmm)
#' @noRd
.plot_covariate_forest <- function(coef_df, title, combined = TRUE) {
  coef_df <- coef_df[coef_df$variable != "(Intercept)", , drop = FALSE]
  coef_df$cluster <- factor(coef_df$cluster)
  coef_df$significant <- coef_df$p < 0.05

  base_p <- function(d, ttl) {
    ggplot(d, aes(
      x = .data$odds_ratio, y = .data$variable,
      xmin = .data$ci_lower, xmax = .data$ci_upper,
      colour = .data$significant
    )) +
      geom_vline(xintercept = 1, linetype = "dashed", alpha = 0.5) +
      geom_pointrange(size = 0.6) +
      scale_colour_manual(
        values = c("TRUE" = "#D62728", "FALSE" = "#7F7F7F"),
        labels = c("TRUE" = "p < 0.05", "FALSE" = "n.s."),
        name = NULL,
        drop = FALSE
      ) +
      scale_x_log10() +
      labs(x = "Odds Ratio (log scale)", y = NULL, title = ttl) +
      theme_minimal(base_size = 12)
  }

  if (!combined) {
    clusters <- levels(coef_df$cluster)
    plots <- lapply(clusters, function(cl) {
      sub <- coef_df[coef_df$cluster == cl, , drop = FALSE]
      base_p(sub, sprintf("%s -- Cluster %s", title, cl))
    })
    names(plots) <- clusters
    return(invisible(plots))
  }

  p <- base_p(coef_df, title) +
    facet_wrap(~ .data$cluster, labeller = label_both)
  print(p)
  invisible(p)
}

# ==============================================================================
# 8. build_network dispatch for net_clustering
# ==============================================================================

#' Build per-cluster networks from clustering results
#'
#' Called by \code{build_network()} when data is a \code{net_clustering}
#' object. Splits sequence data by cluster assignment and builds a network
#' for each cluster.
#'
#' @param x A \code{net_clustering} object.
#' @param method Network estimation method (default: "relative").
#' @param ... Passed to \code{build_network()}.
#' @return A \code{netobject_group}.
#' @noRd
.build_network_clustering <- function(x, method = "relative", ...) {
  assignments <- x$assignments
  seq_data <- x$data
  k <- x$k

  # Merge stored build_args with caller's ...; caller takes precedence
  dots <- list(...)
  build_args <- if (!is.null(x$build_args))
    modifyList(x$build_args, dots) else dots

  nets <- lapply(seq_len(k), function(cl) {
    sub <- seq_data[assignments == cl, , drop = FALSE]
    do.call(build_network, c(list(data = sub, method = method), build_args))
  })

  names(nets) <- paste("Cluster", seq_len(k))
  attr(nets, "group_col") <- "cluster"
  attr(nets, "clustering") <- x
  class(nets) <- "netobject_group"
  nets
}

# ==============================================================================
# 9. cluster_network — one-shot clustering + network estimation
# ==============================================================================

#' Cluster data and build per-cluster networks in one step
#'
#' Combines sequence clustering and network estimation into a single call.
#' Clusters the data using the specified algorithm, then calls
#' \code{\link{build_network}} on each cluster subset.
#'
#' If \code{data} is a \code{netobject} and \code{method} is not provided in
#' \code{...}, the original network method is inherited automatically so the
#' per-cluster networks match the type of the input network.
#'
#' @param data Sequence data. Accepts a data frame, matrix, or
#'   \code{netobject}. See \code{\link{build_clusters}} for supported formats.
#' @param k Integer. Number of clusters.
#' @param cluster_by Character. Clustering algorithm passed to
#'   \code{\link{build_clusters}}'s \code{method} parameter (\code{"pam"},
#'   \code{"ward.D2"}, \code{"ward.D"}, \code{"complete"}, \code{"average"},
#'   \code{"single"}, \code{"mcquitty"}, \code{"median"}, \code{"centroid"}),
#'   or \code{"mmm"} for Mixed Markov Model clustering. Default: \code{"pam"}.
#' @param dissimilarity Character. Distance metric for sequence clustering.
#'   Only valid when \code{cluster_by != "mmm"}. Default: \code{"hamming"}.
#' @param ... Routed to two stages. For distance clustering
#'   (\code{cluster_by != "mmm"}), \code{\link{build_clusters}} arguments
#'   \code{na_syms}, \code{weighted}, \code{lambda}, \code{seed}, \code{q},
#'   \code{p}, and \code{covariates} are intercepted and forwarded to the
#'   clusterer; recognised network-estimation arguments flow to
#'   \code{\link{build_network}}. Unknown arguments error before either stage
#'   is run. When
#'   \code{cluster_by = "mmm"}, the recognised \code{\link{build_mmm}}
#'   arguments (\code{n_starts}, \code{max_iter}, \code{tol}, \code{smooth},
#'   \code{seed}, \code{covariates}) are intercepted instead, and the rest
#'   flows to \code{build_network}. In both modes, when \code{data} is a
#'   \code{netobject}, its \code{build_args} are merged into the
#'   \code{build_network} side (caller's explicit values take precedence).
#' @return A \code{netobject_group}.
#' @seealso \code{\link{build_clusters}}, \code{\link{cluster_mmm}},
#'   \code{\link{build_network}}
#' @examples
#' seqs <- data.frame(V1 = c("A","B","C","A","B"), V2 = c("B","C","A","B","A"),
#'                    V3 = c("C","A","B","C","B"))
#' grp <- cluster_network(seqs, k = 2)
#' grp
#' \donttest{
#' set.seed(1)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:4], 50, TRUE), V2 = sample(LETTERS[1:4], 50, TRUE),
#'   V3 = sample(LETTERS[1:4], 50, TRUE), V4 = sample(LETTERS[1:4], 50, TRUE)
#' )
#' # Default: PAM clustering, relative (transition) networks
#' grp <- cluster_network(seqs, k = 3)
#'
#' # Specify network method (cor requires numeric panel data)
#' \dontrun{
#' panel <- as.data.frame(matrix(rnorm(1500), nrow = 300, ncol = 5))
#' grp <- cluster_network(panel, k = 2, method = "cor")
#' }
#'
#' # MMM-based clustering
#' grp <- cluster_network(seqs, k = 2, cluster_by = "mmm")
#' }
#' @export
cluster_network <- function(data, k, cluster_by = "pam",
                             dissimilarity = "hamming", ...) {
  caller_dots <- list(...)

  # Coerce cograph_network to netobject so the build_args lift below works.
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }

  if (identical(cluster_by, "mmm")) {
    if (!missing(dissimilarity)) {
      stop(
        "`dissimilarity` is only valid for distance clustering. ",
        "Omit it when `cluster_by = \"mmm\"`.",
        call. = FALSE
      )
    }

    # Split caller's `...` into MMM EM hyperparameters vs build_network args
    # so a single cluster_network() call controls both stages without a
    # silently-dropped `n_starts`/`max_iter`/`seed`. Anything else falls
    # through to build_network() unchanged.
    mmm_arg_names  <- c("n_starts", "max_iter", "tol", "smooth",
                        "seed", "covariates", "estimator")
    mmm_args   <- caller_dots[intersect(names(caller_dots), mmm_arg_names)]
    build_dots <- caller_dots[setdiff(names(caller_dots), mmm_arg_names)]
    if (inherits(data, "netobject")) {
      if (!is.null(data$build_args))
        build_dots <- modifyList(data$build_args, build_dots)
      if (is.null(build_dots$method))
        build_dots$method <- data$method
    }
    .validate_cluster_network_build_args(build_dots)
    mmm_fit <- do.call(build_mmm, c(list(data = data, k = k), mmm_args))
    return(do.call(build_network, c(list(data = mmm_fit), build_dots)))
  }

  # Split caller's `...` into build_clusters() args vs build_network() args.
  # Without this split, weighted/lambda/seed/q/p/na_syms/covariates are
  # silently dropped because they reach build_network() instead of the
  # clusterer (audit_clustering finding #1). Splitting on caller_dots only —
  # never on data$build_args — keeps network-history args (e.g. an `atna`
  # net's `lambda` decay) routed to build_network where they belong.
  cluster_arg_names <- c("na_syms", "weighted", "lambda", "seed",
                         "q", "p", "covariates", "estimator")
  cluster_args <- caller_dots[intersect(names(caller_dots), cluster_arg_names)]
  build_dots   <- caller_dots[setdiff(names(caller_dots), cluster_arg_names)]

  if (inherits(data, "netobject")) {
    if (!is.null(data$build_args))
      build_dots <- modifyList(data$build_args, build_dots)
    if (is.null(build_dots$method))
      build_dots$method <- data$method
  }
  .validate_cluster_network_build_args(build_dots)

  cls <- do.call(build_clusters,
                 c(list(data = data, k = k, method = cluster_by,
                        dissimilarity = dissimilarity),
                   cluster_args))
  do.call(build_network, c(list(data = cls), build_dots))
}

.validate_cluster_network_build_args <- function(args) {
  arg_names <- names(args)
  if (is.null(arg_names) || !length(arg_names)) return(invisible(TRUE))

  build_formals <- setdiff(names(formals(build_network)), c("data", "..."))
  estimator_names <- tryCatch(list_estimators()$name, error = function(e) character())
  estimator_formals <- unique(unlist(lapply(estimator_names, function(nm) {
    entry <- tryCatch(get_estimator(nm), error = function(e) NULL)
    if (is.null(entry) || is.null(entry$fn)) return(character())
    setdiff(names(formals(entry$fn)), c("data", "..."))
  }), use.names = FALSE))
  accepted <- unique(c(build_formals, estimator_formals))

  unknown <- setdiff(arg_names, accepted)
  if (length(unknown)) {
    stop(
      "Unknown argument(s) in cluster_network(): ",
      paste(unknown, collapse = ", "), ". ",
      "Pass estimator-specific options through `params = list(...)` ",
      "when using custom estimators.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
