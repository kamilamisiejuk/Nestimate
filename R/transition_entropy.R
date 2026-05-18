# ---- Transition Entropy ----

#' @noRd
.te_row_entropy <- function(P, base) {
  ## 0 * log 0 := 0; vectorised, no loops
  L      <- ifelse(P > 0, log(P, base = base), 0)
  -rowSums(P * L)
}

# ---- transition_entropy() ----

#' Transition Entropy of a Markov Chain
#'
#' @description
#' Computes per-state branching entropy, stationary entropy, and the
#' chain-level entropy rate of a Markov transition process. The entropy rate
#' is the Shannon-McMillan-Breiman per-step uncertainty of trajectories under
#' the stationary distribution; it is the canonical information-theoretic
#' summary of a transition matrix.
#'
#' @param x A \code{netobject}, \code{cograph_network}, \code{tna} object,
#'   row-stochastic numeric transition matrix, or a wide sequence data.frame
#'   (rows = actors, columns = time-steps; a relative transition network is
#'   built automatically). Group dispatch on \code{netobject_group}.
#' @param base Numeric. Logarithm base. \code{2} (default) for bits,
#'   \code{exp(1)} for nats, \code{10} for hartleys.
#' @param normalize Logical. If \code{TRUE} (default), rows that do not sum
#'   to 1 are normalised automatically (with a warning).
#'
#' @return An object of class \code{"net_transition_entropy"} with:
#' \describe{
#'   \item{row_entropy}{Named numeric vector, length \eqn{n}. Per-state
#'     branching entropy \eqn{H(P_{i\cdot}) = -\sum_j P_{ij} \log P_{ij}}.}
#'   \item{row_entropy_norm}{Named numeric vector. \code{row_entropy}
#'     divided by the ceiling \eqn{\log_b n} (in \eqn{[0, 1]}; all zeros
#'     when \eqn{n = 1}).}
#'   \item{stationary}{Named numeric vector. Stationary distribution
#'     \eqn{\pi}.}
#'   \item{stationary_entropy}{Scalar. \eqn{H(\pi) = -\sum_i \pi_i \log \pi_i}
#'     - the entropy of \eqn{\pi} treated as an i.i.d. distribution. Upper
#'     bound on the entropy rate.}
#'   \item{stationary_entropy_norm}{Scalar. \code{stationary_entropy}
#'     divided by the ceiling \eqn{\log_b n}.}
#'   \item{entropy_rate}{Scalar. \eqn{h(P) = \sum_i \pi_i H(P_{i\cdot})} -
#'     the Shannon-McMillan-Breiman entropy rate.}
#'   \item{entropy_rate_norm}{Scalar. \code{entropy_rate} divided by the
#'     ceiling \eqn{\log_b n}.}
#'   \item{redundancy}{Scalar. \eqn{H(\pi) - h(P)}, the entropy deficit
#'     attributable to serial dependence; zero for an i.i.d. chain (rows of
#'     \eqn{P} all equal \eqn{\pi}).}
#'   \item{redundancy_norm}{Scalar. The \emph{relative} redundancy
#'     \eqn{(H(\pi) - h(P)) / H(\pi)} (the fraction of the stationary
#'     entropy removed by order-1 memory), \strong{not}
#'     \code{redundancy} divided by \eqn{\log_b n}; \code{0} when
#'     \eqn{H(\pi) = 0}.}
#'   \item{max_entropy}{Scalar. The normalising ceiling \eqn{\log_b n}.}
#'   \item{base}{Logarithm base used.}
#'   \item{states}{Character vector of state names.}
#' }
#'
#' @details
#' Convention \eqn{0 \log 0 := 0} is applied, so absorbing or
#' deterministic rows contribute zero per-row entropy. The chain need not be
#' irreducible; \eqn{\pi} is computed from the eigendecomposition of
#' \eqn{P^\top} as elsewhere in the package. For non-ergodic chains the
#' returned \eqn{\pi} is one stationary distribution among many - interpret
#' with the help of \code{\link{chain_structure}}.
#'
#' The relation \eqn{h(P) \leq H(\pi)} holds with equality iff successive
#' states are independent. The deficit \eqn{H(\pi) - h(P)} is reported as
#' \code{redundancy} - a measure of how much memory the chain has at order 1.
#'
#' @examples
#' \donttest{
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' te  <- transition_entropy(net)
#' print(te)
#' summary(te)
#' plot(te)
#' }
#'
#' @seealso \code{\link{markov_stability}}, \code{\link{passage_time}},
#'   \code{\link{markov_order_test}}, \code{\link{chain_structure}}
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information Theory},
#' 2nd ed., chapter 4. Wiley.
#'
#' Shannon, C.E. (1948). A mathematical theory of communication.
#' \emph{Bell System Technical Journal}, 27, 379-423.
#'
#' @export
transition_entropy <- function(x, base = 2, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      transition_entropy(net, base = base, normalize = normalize)
    })
    class(out) <- c("net_transition_entropy_group", "list")
    return(out)
  }

  stopifnot(is.numeric(base), length(base) == 1L, base > 0, base != 1)

  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  pi <- .mpt_stationary(P)
  names(pi) <- state_names
  if (any(pi <= 0)) {
    warning("Non-positive stationary probabilities; chain may not be ergodic.",
            call. = FALSE)
    pi <- pmax(pi, .Machine$double.eps)
    pi <- pi / sum(pi)
  }

  H_row <- .te_row_entropy(P, base = base)
  names(H_row) <- state_names

  H_pi  <- {
    p     <- pi[pi > 0]
    -sum(p * log(p, base = base))
  }
  h_P   <- sum(pi * H_row)
  H_max <- log(length(state_names), base = base)

  ## Normalised forms (unit-free, in [0, 1]); H_max == 0 only when n == 1
  norm <- function(x) if (H_max > 0) x / H_max else x * 0
  H_row_norm <- norm(H_row)
  H_pi_norm  <- norm(H_pi)
  h_P_norm   <- norm(h_P)
  redund_norm <- if (H_pi > 0) (H_pi - h_P) / H_pi else 0

  structure(
    list(
      row_entropy             = H_row,
      row_entropy_norm        = H_row_norm,
      stationary              = pi,
      stationary_entropy      = H_pi,
      stationary_entropy_norm = H_pi_norm,
      entropy_rate            = h_P,
      entropy_rate_norm       = h_P_norm,
      redundancy              = H_pi - h_P,
      redundancy_norm         = redund_norm,
      max_entropy             = H_max,
      base                    = base,
      states                  = state_names
    ),
    class = "net_transition_entropy"
  )
}

#' Print method for `net_transition_entropy`
#'
#' @param x A `net_transition_entropy` object.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_transition_entropy <- function(x, digits = 3, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Transition Entropy (%d states, %s; ceiling = %.3f)\n\n",
              length(x$states), unit, x$max_entropy))
  cat(sprintf("                          raw            normalised\n"))
  cat(sprintf("  Entropy rate    h(P)  = %.*f %-6s  %.3f\n",
              digits, x$entropy_rate, unit, x$entropy_rate_norm))
  cat(sprintf("  Stationary    H(pi)  = %.*f %-6s  %.3f\n",
              digits, x$stationary_entropy, unit, x$stationary_entropy_norm))
  cat(sprintf("  Redundancy   H(pi)-h = %.*f %-6s  %.3f\n",
              digits, x$redundancy, unit, x$redundancy_norm))
  cat(sprintf(paste0(
    "\nNormalised: h(P) and H(pi) are raw / log_%g(n_states) ",
    "(0 = deterministic, 1 = uniform);\n",
    "  redundancy is the relative redundancy (H(pi) - h(P)) / H(pi), ",
    "not raw / log_%g(n_states).\n"),
    x$base, x$base))
  invisible(x)
}

#' Print method for `net_transition_entropy_group`
#'
#' @param x A `net_transition_entropy_group`.
#' @param ... Forwarded to `print.net_transition_entropy`.
#' @return `x` invisibly.
#' @export
print.net_transition_entropy_group <- function(x, ...) {
  cat(sprintf("Transition Entropy - %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  rates <- vapply(x, function(o) o$entropy_rate, numeric(1))
  cat("Entropy rate per group:\n")
  print(round(rates, 4))
  invisible(x)
}

#' Summary method for `net_transition_entropy`
#'
#' @description
#' Returns a tidy per-state contribution table sorted by share of the
#' chain-level entropy rate (largest first), so the dominant contributors
#' to \eqn{h(P)} are visible at a glance. Each row contains the
#' stationary mass, the raw and normalised row entropy, the additive
#' contribution \eqn{\pi_i H(P_{i\cdot})}, and that contribution as a
#' percentage of \eqn{h(P)}.
#'
#' @param object A `net_transition_entropy` object.
#' @param ... Ignored.
#' @return A `summary.net_transition_entropy` containing
#'   \describe{
#'     \item{table}{tidy per-state data.frame, sorted by `contribution_pct`
#'       descending}
#'     \item{chain}{tidy chain-level data.frame with raw and normalised
#'       \eqn{h(P)}, \eqn{H(\pi)}, redundancy, and ceiling}
#'     \item{base}{logarithm base used}
#'   }
#' @export
summary.net_transition_entropy <- function(object, ...) {
  contrib  <- object$stationary * object$row_entropy
  total    <- object$entropy_rate
  pct      <- if (total > 0) 100 * contrib / total else contrib * 0

  per_state <- data.frame(
    state            = object$states,
    stationary       = unname(object$stationary),
    row_entropy      = unname(object$row_entropy),
    row_entropy_norm = unname(object$row_entropy_norm),
    contribution     = unname(contrib),
    contribution_pct = unname(pct),
    stringsAsFactors = FALSE
  )
  per_state <- per_state[order(-per_state$contribution_pct), ]
  rownames(per_state) <- NULL

  unit <- switch(as.character(object$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", object$base))
  if (isTRUE(all.equal(object$base, exp(1)))) unit <- "nats"

  chain <- data.frame(
    quantity = c("entropy_rate h(P)",
                 "stationary H(pi)",
                 "redundancy H(pi)-h(P)",
                 sprintf("ceiling log_%g(n)", object$base)),
    raw        = c(object$entropy_rate,
                   object$stationary_entropy,
                   object$redundancy,
                   object$max_entropy),
    normalised = c(object$entropy_rate_norm,
                   object$stationary_entropy_norm,
                   object$redundancy_norm,
                   1),
    unit       = unit,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      table = per_state,
      chain = chain,
      base  = object$base
    ),
    class = "summary.net_transition_entropy"
  )
}

#' Print method for `summary.net_transition_entropy`
#'
#' @param x A `summary.net_transition_entropy`.
#' @param digits Integer. Digits to round numeric output. Default `3`.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.summary.net_transition_entropy <- function(x, digits = 3, ...) {
  unit <- x$chain$unit[1L]

  cat(sprintf("Transition Entropy Summary (%s)\n\n", unit))

  cat("Per-state contribution to h(P):\n")
  tbl <- x$table
  tbl$stationary       <- round(tbl$stationary,       digits)
  tbl$row_entropy      <- round(tbl$row_entropy,      digits)
  tbl$row_entropy_norm <- round(tbl$row_entropy_norm, digits)
  tbl$contribution     <- round(tbl$contribution,     digits)
  tbl$contribution_pct <- round(tbl$contribution_pct, 1)
  print(tbl, row.names = FALSE)

  cat("\nChain-level summary:\n")
  ch <- x$chain
  ch$raw        <- round(ch$raw,        digits)
  ch$normalised <- round(ch$normalised, digits)
  ch$unit       <- NULL
  print(ch, row.names = FALSE)

  cat(sprintf(paste0(
    "\nNormalised: row entropy, h(P) and H(pi) are raw / log_%g(n_states), ",
    "in [0, 1];\n",
    "  redundancy is the relative redundancy (H(pi) - h(P)) / H(pi).\n"),
    x$base))
  invisible(x)
}

#' Plot method for `net_transition_entropy`
#'
#' @description
#' Bar chart of per-state row entropy with overlaid horizontal lines at
#' the entropy rate \eqn{h(P)} (chain-level summary) and the maximum row
#' entropy \eqn{\log_b n} (uniform branching). Bar widths are proportional
#' to the stationary probability so the visual area sums to the entropy
#' rate.
#'
#' @param x A `net_transition_entropy` object.
#' @param title Character. Plot title.
#' @param fill Character. Bar fill colour. Default Okabe-Ito blue.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.net_transition_entropy <- function(x,
                                        title = "Transition Entropy",
                                        fill  = "#0072B2",
                                        ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  n      <- length(x$states)
  H_max  <- log(n, base = x$base)

  df <- data.frame(
    state       = factor(x$states, levels = x$states),
    row_entropy = x$row_entropy,
    stationary  = x$stationary,
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$state,
                                   y = .data$row_entropy,
                                   width = .data$stationary /
                                     max(.data$stationary))) +
    ggplot2::geom_col(fill = fill, alpha = 0.85) +
    ggplot2::geom_hline(yintercept = x$entropy_rate,
                        linetype = "dashed", colour = "#D55E00",
                        linewidth = 0.6) +
    ggplot2::geom_hline(yintercept = H_max,
                        linetype = "dotted", colour = "grey30",
                        linewidth = 0.4) +
    ggplot2::annotate("text",
                      x = 0.6, y = x$entropy_rate,
                      label = sprintf("h(P) = %.3f", x$entropy_rate),
                      hjust = 0, vjust = -0.4,
                      colour = "#D55E00", size = 3.4) +
    ggplot2::annotate("text",
                      x = 0.6, y = H_max,
                      label = sprintf("log_%g(n) = %.3f", x$base, H_max),
                      hjust = 0, vjust = -0.4,
                      colour = "grey30", size = 3.2) +
    ggplot2::labs(x = "State",
                  y = sprintf("Row entropy (%s)", unit),
                  title    = title,
                  subtitle = "Bar width proportional to stationary probability") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = 30, hjust = 1),
      plot.title         = ggplot2::element_text(face = "bold")
    )
}
