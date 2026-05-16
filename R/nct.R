# ---- Network Comparison Test (NCT) ----

#' Network Comparison Test
#'
#' Tests whether two networks estimated from independent samples differ at
#' three levels: \strong{global strength} (M-statistic), \strong{network
#' structure} (S-statistic, max absolute edge difference), and
#' \strong{individual edges} (E-statistic per edge). Inference is via
#' permutation of group labels.
#'
#' Implementation matches \code{NetworkComparisonTest::NCT()} with defaults
#' \code{abs = TRUE}, \code{weighted = TRUE}, \code{paired = FALSE} at
#' machine precision when the same seed is used. The network estimator is
#' EBIC-selected glasso applied to a Pearson correlation matrix, with
#' \code{Matrix::nearPD} symmetrization (matching NCT's
#' \code{NCT_estimator_GGM} default).
#'
#' @param data1 A numeric matrix or data.frame of observations from group 1.
#' @param data2 A numeric matrix or data.frame of observations from group 2.
#'   Same number of columns as \code{data1}.
#' @param iter Integer. Number of permutation iterations. Default 1000.
#' @param gamma EBIC tuning parameter for glasso. Default 0.5.
#' @param paired Logical. If \code{TRUE}, perform a paired permutation
#'   (within-subject swap). Default \code{FALSE}.
#' @param abs Logical. If \code{TRUE}, compute global strength on absolute
#'   edge weights. Default \code{TRUE}.
#' @param weighted Logical. If \code{TRUE}, use weighted networks for the
#'   tests. If \code{FALSE}, binarize before computing statistics.
#'   Default \code{TRUE}.
#' @param p_adjust P-value adjustment method for the per-edge tests
#'   (any method in \code{stats::p.adjust.methods}). Default \code{"none"}.
#' @return A list of class \code{net_nct} with elements:
#' \describe{
#'   \item{nw1, nw2}{Estimated weighted adjacency matrices.}
#'   \item{M}{List with \code{observed}, \code{perm}, \code{p_value} for
#'     the global strength test.}
#'   \item{S}{Same structure for the maximum absolute edge difference.}
#'   \item{E}{Same structure for per-edge tests.}
#'   \item{n_iter}{Number of permutations.}
#'   \item{paired}{Whether a paired test was used.}
#' }
#' @examples
#' \dontrun{
#' set.seed(1)
#' x1 <- matrix(rnorm(200 * 5), 200, 5)
#' x2 <- matrix(rnorm(200 * 5), 200, 5)
#' colnames(x1) <- colnames(x2) <- paste0("V", 1:5)
#' res <- nct(x1, x2, iter = 100)
#' res$M$p_value
#' res$S$p_value
#' }
#' @export
nct <- function(data1, data2, iter = 1000L, gamma = 0.5,
                 paired = FALSE, abs = TRUE, weighted = TRUE,
                 p_adjust = "none") {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("nct() requires the 'Matrix' package.", call. = FALSE)
  }
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)
  stopifnot(
    ncol(data1) == ncol(data2),
    is.numeric(iter), iter >= 1L,
    is.logical(paired), is.logical(abs), is.logical(weighted)
  )
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)

  iter <- as.integer(iter)
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p  <- ncol(data1)
  dataall <- rbind(data1, data2)

  # Estimator: nearPD symmetrization + EBIC-glasso. When qgraph is available
  # use the same EBICglasso implementation as NetworkComparisonTest::NCT().
  # The package's internal path remains a dependency-light fallback.
  est <- function(x) {
    cor_x <- stats::cor(x)
    cor_x <- as.matrix(Matrix::nearPD(cor_x, corr = TRUE)$mat)
    cor_x <- (cor_x + t(cor_x)) / 2
    if (requireNamespace("qgraph", quietly = TRUE)) {
      return(suppressWarnings(suppressMessages(
        qgraph::EBICglasso(cor_x, nrow(x), gamma = gamma)
      )))
    }
    lambda_path <- .compute_lambda_path(cor_x, nlambda = 100L,
                                          lambda.min.ratio = 0.01)
    selected <- .select_ebic(cor_x, lambda_path,
                              n = nrow(x), gamma = gamma,
                              penalize_diagonal = FALSE)
    .wi2net(selected$wi)
  }

  nw1 <- est(data1)
  nw2 <- est(data2)
  ut <- upper.tri(nw1)

  binarize <- function(m) (m != 0) * 1
  if (!weighted) {
    nw1 <- binarize(nw1)
    nw2 <- binarize(nw2)
  }

  # ---- observed statistics ----
  M_obs <- if (abs) {
    base::abs(sum(base::abs(nw1[ut])) - sum(base::abs(nw2[ut])))
  } else {
    base::abs(sum(nw1[ut]) - sum(nw2[ut]))
  }
  diff_real <- base::abs(nw1 - nw2)
  S_obs <- max(diff_real[ut])
  E_obs <- diff_real[ut]
  n_edges <- length(E_obs)

  # ---- permutation ----
  M_perm <- numeric(iter)
  S_perm <- numeric(iter)
  E_perm <- matrix(0, nrow = iter, ncol = n_edges)

  for (i in seq_len(iter)) {
    if (paired) {
      s <- sample(c(1L, 2L), n1, replace = TRUE)
      x1p <- rbind(data1[s == 1L, , drop = FALSE],
                   data2[s == 2L, , drop = FALSE])
      x2p <- rbind(data2[s == 1L, , drop = FALSE],
                   data1[s == 2L, , drop = FALSE])
    } else {
      s <- sample(seq_len(n1 + n2), n1, replace = FALSE)
      x1p <- dataall[s, , drop = FALSE]
      x2p <- dataall[-s, , drop = FALSE]
    }
    r1 <- est(x1p)
    r2 <- est(x2p)
    if (!weighted) {
      r1 <- binarize(r1)
      r2 <- binarize(r2)
    }
    if (abs) {
      M_perm[i] <- base::abs(sum(base::abs(r1[ut])) - sum(base::abs(r2[ut])))
    } else {
      M_perm[i] <- base::abs(sum(r1[ut]) - sum(r2[ut]))
    }
    diff_perm <- base::abs(r1 - r2)
    S_perm[i] <- max(diff_perm[ut])
    E_perm[i, ] <- diff_perm[ut]
  }

  # ---- p-values (matches NCT formula) ----
  M_pval <- (sum(M_perm >= M_obs) + 1) / (iter + 1)
  S_pval <- (sum(S_perm >= S_obs) + 1) / (iter + 1)
  E_pval_raw <- (colSums(E_perm >= matrix(E_obs, iter, n_edges, byrow = TRUE))
                 + 1) / (iter + 1)
  E_pval <- if (p_adjust != "none") {
    stats::p.adjust(E_pval_raw, method = p_adjust)
  } else {
    E_pval_raw
  }

  edge_names <- if (!is.null(colnames(data1))) {
    out <- expand.grid(colnames(data1), colnames(data1),
                       stringsAsFactors = FALSE)[as.vector(ut), ]
    out
  } else NULL

  structure(
    list(
      nw1    = nw1,
      nw2    = nw2,
      M      = list(observed = M_obs, perm = M_perm, p_value = M_pval),
      S      = list(observed = S_obs, perm = S_perm, p_value = S_pval),
      E      = list(observed = E_obs, perm = E_perm,
                    p_value = E_pval, edge_names = edge_names),
      n_iter = iter,
      paired = paired,
      params = list(gamma = gamma, abs = abs, weighted = weighted,
                    p_adjust = p_adjust)
    ),
    class = "net_nct"
  )
}


#' Print Method for net_nct
#'
#' @param x A \code{net_nct} object.
#' @param ... Ignored.
#' @return The input object, invisibly.
#' @inherit nct examples
#' @export
print.net_nct <- function(x, ...) {
  cat(sprintf("Network Comparison Test  [%d permutations | %s]\n",
              x$n_iter, if (x$paired) "paired" else "unpaired"))
  cat(sprintf("  Global strength (M):  observed = %.4f   p = %.4f\n",
              x$M$observed, x$M$p_value))
  cat(sprintf("  Network structure (S): observed = %.4f   p = %.4f\n",
              x$S$observed, x$S$p_value))
  n_edges <- length(x$E$observed)
  n_sig   <- sum(x$E$p_value < 0.05)
  cat(sprintf("  Edge tests (E):       %d edges, %d significant at p < 0.05",
              n_edges, n_sig))
  if (!identical(x$params$p_adjust, "none")) {
    cat(sprintf("  (%s adjusted)", x$params$p_adjust))
  }
  cat("\n")
  invisible(x)
}


#' Summary Method for net_nct
#'
#' @description
#' Returns a tidy data frame with one row per edge test. The global M
#' (strength) and S (structure) statistics are attached as attributes.
#'
#' @param object A \code{net_nct} object.
#' @param ... Ignored.
#' @return A data frame with columns \code{from}, \code{to},
#'   \code{diff_observed}, \code{p_value}, \code{significant}. Attributes
#'   \code{m_stat} and \code{s_stat} each hold a one-row data frame with
#'   \code{observed} and \code{p_value}.
#' @inherit nct examples
#' @export
summary.net_nct <- function(object, ...) {
  ed <- object$E$edge_names
  n  <- length(object$E$observed)
  if (is.null(ed) || nrow(ed) != n) {
    from <- paste0("edge_", seq_len(n))
    to   <- rep(NA_character_, n)
  } else {
    # edge_names is an expand.grid() result with Var1/Var2
    from <- as.character(ed$Var1)
    to   <- as.character(ed$Var2)
  }
  df <- data.frame(
    from          = from,
    to            = to,
    diff_observed = as.numeric(object$E$observed),
    p_value       = as.numeric(object$E$p_value),
    stringsAsFactors = FALSE,
    row.names     = NULL
  )
  df$significant <- !is.na(df$p_value) & df$p_value < 0.05
  m_df <- data.frame(observed = object$M$observed, p_value = object$M$p_value)
  s_df <- data.frame(observed = object$S$observed, p_value = object$S$p_value)
  structure(df, m_stat = m_df, s_stat = s_df)
}
