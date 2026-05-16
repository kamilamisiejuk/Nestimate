# ---- Markov Chain Analysis: Passage Times and Stability ----

#' @noRd
.mpt_extract_P <- function(x) {
  if (is.matrix(x) && is.numeric(x)) return(x)
  if (inherits(x, "netobject") || inherits(x, "cograph_network") ||
      inherits(x, "tna")) {
    P <- x$weights
    if (is.null(P) || !is.matrix(P) || !is.numeric(P))
      stop("Object has no numeric weight matrix.", call. = FALSE)
    return(P)
  }
  if (is.matrix(x) && (is.character(x) || is.logical(x))) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (is.data.frame(x)) {
    net <- build_network(x, method = "relative")
    return(net$weights)
  }
  stop(
    "'x' must be a numeric matrix (transition matrix), a netobject / ",
    "cograph_network / tna object, or a wide sequence data.frame / ",
    "character matrix.",
    call. = FALSE
  )
}

#' @noRd
.mpt_normalize_rows <- function(P, state_names, normalize = TRUE) {
  row_sums <- rowSums(P)

  zero_rows <- which(row_sums <= 0)
  if (length(zero_rows)) {
    bad <- state_names[zero_rows]
    stop(
      "Transition matrix has zero-sum row(s) for state(s): ",
      paste(bad, collapse = ", "),
      ". These states have no outgoing transitions, so the chain is not ",
      "ergodic and mean first passage times are undefined. Remove the ",
      "state(s) or supply a different transition matrix.",
      call. = FALSE
    )
  }

  if (any(abs(row_sums - 1) > 1e-6)) {
    if (!normalize)
      stop("Transition matrix rows must sum to 1.", call. = FALSE)
    warning("Rows do not sum to 1; normalizing.", call. = FALSE)
    P <- P / row_sums
  }

  P
}

#' @noRd
.mpt_stationary <- function(P) {
  ev  <- eigen(t(P))
  idx <- which.min(abs(Re(ev$values) - 1))
  pi  <- abs(Re(ev$vectors[, idx]))
  pi / sum(pi)
}

#' @noRd
.mpt_full <- function(P, pi) {
  n     <- nrow(P)
  PI    <- matrix(pi, nrow = n, ncol = n, byrow = TRUE)
  Z     <- solve(diag(n) - P + PI)
  # Kemeny-Snell: M[i,j] = (Z[j,j] - Z[i,j]) / pi[j]
  Zd    <- matrix(diag(Z), nrow = n, ncol = n, byrow = TRUE)
  pimat <- matrix(pi,      nrow = n, ncol = n, byrow = TRUE)
  M     <- (Zd - Z) / pimat
  diag(M) <- 1 / pi
  M
}

# ---- passage_time() ----

#' Mean First Passage Times
#'
#' @description
#' Computes the full matrix of mean first passage times (MFPT) for a Markov
#' chain. Element \eqn{M_{ij}} is the expected number of steps to travel from
#' state \eqn{i} to state \eqn{j} for the first time. The diagonal equals
#' the mean recurrence time \eqn{1/\pi_i}.
#'
#' @param x A \code{netobject}, \code{cograph_network}, \code{tna} object,
#'   row-stochastic numeric transition matrix, or a wide sequence data.frame
#'   (rows = actors, columns = time-steps; a relative transition network is
#'   built automatically).
#' @param object A \code{net_mpt} object (for \code{summary}).
#' @param states Character vector. Restrict output to these states.
#'   \code{NULL} (default) keeps all states.
#' @param normalize Logical. If \code{TRUE} (default), rows that do not sum
#'   to 1 are normalized automatically (with a warning).
#' @param ... Ignored.
#'
#' @return An object of class \code{"net_mpt"} with:
#' \describe{
#'   \item{matrix}{Full \eqn{n \times n} MFPT matrix. Row \eqn{i}, column
#'     \eqn{j} = expected steps from state \eqn{i} to state \eqn{j}.
#'     Diagonal = mean recurrence time \eqn{1/\pi_i}.}
#'   \item{stationary}{Named numeric vector: stationary distribution \eqn{\pi}.}
#'   \item{return_times}{Named numeric vector: \eqn{1/\pi_i} per state.}
#'   \item{states}{Character vector of state names.}
#' }
#'
#' @details
#' Uses the Kemeny-Snell fundamental matrix formula:
#' \deqn{M_{ij} = \frac{Z_{jj} - Z_{ij}}{\pi_j}, \quad
#'       Z = (I - P + \Pi)^{-1}}
#' where \eqn{\Pi_{ij} = \pi_j}. Requires an ergodic (irreducible,
#' aperiodic) chain.
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' pt  <- passage_time(net)
#' print(pt)
#' \donttest{
#' plot(pt)
#' }
#'
#' @seealso \code{\link{markov_stability}}, \code{\link{build_network}}
#'
#' @references
#' Kemeny, J.G. and Snell, J.L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' @export
passage_time <- function(x, states = NULL, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      passage_time(net, states = states, normalize = normalize)
    })
    class(out) <- c("net_mpt_group", "list")
    return(out)
  }
  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names     <- paste0("S", seq_len(nrow(P)))
    colnames(P)     <- rownames(P) <- state_names
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

  M           <- .mpt_full(P, pi)
  dimnames(M) <- list(state_names, state_names)

  if (!is.null(states)) {
    bad <- setdiff(states, state_names)
    if (length(bad))
      stop("Unknown states: ", paste(bad, collapse = ", "), call. = FALSE)
    M  <- M[states, states, drop = FALSE]
    pi <- pi[states]
  }

  structure(
    list(
      matrix       = M,
      stationary   = pi,
      return_times = 1 / pi,
      states       = rownames(M)
    ),
    class = "net_mpt"
  )
}

#' @export
print.net_mpt <- function(x, digits = 1, ...) {
  n <- length(x$states)
  cat(sprintf("Mean First Passage Times (%d states)\n\n", n))
  print(round(x$matrix, digits))
  cat("\nStationary distribution:\n")
  print(round(x$stationary, 4))
  invisible(x)
}

#' Print method for `net_mpt_group`
#'
#' @param x A `net_mpt_group` (named list of `net_mpt` results).
#' @param ... Forwarded to `print.net_mpt` for each element.
#' @return `x` invisibly.
#' @export
print.net_mpt_group <- function(x, ...) {
  cat(sprintf("Mean First Passage Times -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}

#' @return \code{summary.net_mpt} returns a data frame with one row per state
#'   and columns \code{state}, \code{return_time}, \code{stationary},
#'   \code{mean_out} (mean steps to other states), \code{mean_in} (mean steps
#'   from other states).
#' @rdname passage_time
#' @export
summary.net_mpt <- function(object, ...) {
  M  <- object$matrix
  n  <- nrow(M)
  st <- object$states

  df <- data.frame(
    state       = st,
    return_time = round(object$return_times, 2),
    stationary  = round(object$stationary,   4),
    mean_out    = vapply(seq_len(n),
                         function(i) round(mean(M[i, -i]), 2), numeric(1)),
    mean_in     = vapply(seq_len(n),
                         function(i) round(mean(M[-i, i]), 2), numeric(1)),
    stringsAsFactors = FALSE
  )
  rownames(df) <- NULL
  structure(list(table = df, object = object), class = "summary.net_mpt")
}

#' @export
print.summary.net_mpt <- function(x, ...) {
  cat("Mean First Passage Times - Summary\n\n")
  print(x$table, row.names = FALSE)
  invisible(x)
}

#' @param log_scale Logical. Apply log transform to the fill scale for better
#'   contrast? Default \code{TRUE}.
#' @param digits Integer. Decimal places displayed in cells. Default \code{1}.
#' @param title Character. Plot title.
#' @param low Character. Hex colour for the low end (short passage time).
#'   Default dark green \code{"#004d00"}.
#' @param high Character. Hex colour for the high end (long passage time).
#'   Default pale green \code{"#ccffcc"}.
#' @rdname passage_time
#' @export
plot.net_mpt <- function(x,
                          log_scale = TRUE,
                          digits    = 1,
                          title     = "Mean First Passage Times",
                          low       = "#004d00",
                          high      = "#ccffcc",
                          ...) {
  states <- x$states
  M      <- x$matrix
  n      <- length(states)

  from  <- factor(rep(states, each  = n), levels = rev(states))
  to    <- factor(rep(states, times = n), levels = states)
  vals  <- as.vector(M)
  fill  <- if (log_scale) log(pmax(vals, .Machine$double.eps)) else vals
  label <- formatC(round(vals, digits), format = "f", digits = digits,
                   flag = " ")

  rng <- range(fill, na.rm = TRUE)
  mid <- mean(rng)
  txt_col <- ifelse(fill <= mid, "white", "#003300")

  df <- data.frame(from = from, to = to, fill = fill,
                   label = label, txt_col = txt_col,
                   stringsAsFactors = FALSE)

  ggplot2::ggplot(df, ggplot2::aes(x = to, y = from, fill = fill)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = txt_col),
                       size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient(
      low  = low, high = high,
      name = if (log_scale) "log(Steps)" else "Steps",
      labels = if (log_scale) function(b) round(exp(b), 1) else ggplot2::waiver()
    ) +
    ggplot2::labs(x = "To", y = "From", title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text  = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


# ---- markov_stability() ----

#' Markov Stability Analysis
#'
#' @description
#' Computes per-state stability metrics from a transition network:
#' persistence (self-loop probability), stationary distribution, mean
#' recurrence time, sojourn time, and mean accessibility to/from other states.
#'
#' @param x A \code{netobject}, \code{cograph_network}, \code{tna} object,
#'   row-stochastic numeric transition matrix, or a wide sequence data.frame
#'   (rows = actors, columns = time-steps).
#' @param normalize Logical. Normalize rows to sum to 1? Default \code{TRUE}.
#' @param ... Ignored.
#'
#' @return An object of class \code{"net_markov_stability"} with:
#' \describe{
#'   \item{stability}{Data frame with one row per state and columns:
#'     \code{state}, \code{persistence} (\eqn{P_{ii}}),
#'     \code{stationary_prob} (\eqn{\pi_i}),
#'     \code{return_time} (\eqn{1/\pi_i}),
#'     \code{sojourn_time} (\eqn{1/(1-P_{ii})}),
#'     \code{avg_time_to_others} (mean MFPT leaving state \eqn{i}),
#'     \code{avg_time_from_others} (mean MFPT arriving at state \eqn{i}).}
#'   \item{mpt}{The underlying \code{net_mpt} object.}
#' }
#'
#' @details
#' \strong{Sojourn time} is the expected consecutive time steps spent in a
#' state before leaving: \eqn{1/(1-P_{ii})}. States with
#' \code{persistence = 1} have \code{sojourn_time = Inf}.
#'
#' \strong{avg_time_to_others}: mean passage time from this state to all
#' others; reflects how "sticky" or "isolated" the state is.
#'
#' \strong{avg_time_from_others}: mean passage time from all other states
#' to this one; reflects accessibility (attractor strength).
#'
#' @examples
#' net <- build_network(as.data.frame(trajectories), method = "relative")
#' ms  <- markov_stability(net)
#' print(ms)
#' \donttest{
#' plot(ms)
#' }
#'
#' @seealso \code{\link{passage_time}}
#' @export
markov_stability <- function(x, normalize = TRUE) {
  if (inherits(x, "netobject_group")) {
    out <- lapply(x, function(net) {
      markov_stability(net, normalize = normalize)
    })
    class(out) <- c("net_markov_stability_group", "list")
    return(out)
  }
  P           <- .mpt_extract_P(x)
  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names     <- paste0("S", seq_len(nrow(P)))
    colnames(P)     <- rownames(P) <- state_names
  }

  P <- .mpt_normalize_rows(P, state_names, normalize = normalize)

  mpt  <- passage_time(P, normalize = FALSE)
  M    <- mpt$matrix
  pi   <- mpt$stationary
  n    <- nrow(P)

  persistence <- diag(P)
  sojourn     <- ifelse(persistence >= 1 - .Machine$double.eps,
                        Inf, 1 / (1 - persistence))
  avg_to   <- vapply(seq_len(n), function(i) mean(M[i, -i]), numeric(1))
  avg_from <- vapply(seq_len(n), function(i) mean(M[-i, i]), numeric(1))

  stability_df <- data.frame(
    state                = state_names,
    persistence          = round(persistence,  4),
    stationary_prob      = round(pi,            4),
    return_time          = round(1 / pi,        2),
    sojourn_time         = round(sojourn,       2),
    avg_time_to_others   = round(avg_to,        2),
    avg_time_from_others = round(avg_from,      2),
    stringsAsFactors     = FALSE
  )
  rownames(stability_df) <- NULL

  structure(list(stability = stability_df, mpt = mpt),
            class = "net_markov_stability")
}

#' @export
print.net_markov_stability <- function(x, ...) {
  cat("Markov Stability Analysis\n\n")
  print(x$stability, row.names = FALSE)
  invisible(x)
}

#' Print method for `net_markov_stability_group`
#'
#' @param x A `net_markov_stability_group` (named list of
#'   `net_markov_stability` results).
#' @param ... Forwarded to `print.net_markov_stability` for each element.
#' @return `x` invisibly.
#' @export
print.net_markov_stability_group <- function(x, ...) {
  cat(sprintf("Markov Stability -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}

#' @export
summary.net_markov_stability <- function(object, ...) {
  df      <- object$stability
  attract <- df$state[which.min(df$avg_time_from_others)]
  sticky  <- df$state[which.max(df$sojourn_time)]
  cat(sprintf("Most accessible state (attractor): %s\n", attract))
  cat(sprintf("Most persistent state (stickiest): %s\n\n", sticky))
  df
}

#' @param metrics Character vector. Which metrics to plot. Options:
#'   \code{"persistence"}, \code{"stationary_prob"}, \code{"return_time"},
#'   \code{"sojourn_time"}, \code{"avg_time_to_others"},
#'   \code{"avg_time_from_others"}. Default: all six.
#' @param combined When \code{TRUE} (default), all selected metrics are
#'   shown in one ggplot via \code{facet_wrap(~ metric)}. When \code{FALSE},
#'   returns a named list of single-panel ggplots, one per metric, so each
#'   can be printed, saved, or re-laid-out independently.
#' @rdname markov_stability
#' @export
plot.net_markov_stability <- function(x,
                                       metrics = c("persistence",
                                                   "stationary_prob",
                                                   "return_time",
                                                   "sojourn_time",
                                                   "avg_time_to_others",
                                                   "avg_time_from_others"),
                                       combined = TRUE,
                                       ...) {
  df <- x$stability
  metrics <- match.arg(metrics, several.ok = TRUE)
  stopifnot(is.logical(combined), length(combined) == 1L)

  labels <- c(
    persistence          = "Persistence",
    stationary_prob      = "Stationary Prob.",
    return_time          = "Return Time",
    sojourn_time         = "Sojourn Time",
    avg_time_to_others   = "Avg. Steps to Others",
    avg_time_from_others = "Avg. Steps from Others"
  )

  st_ord <- factor(df$state, levels = df$state[order(df$stationary_prob)])
  n_states <- nrow(df)

  pal <- c("#009E73", "#E69F00", "#56B4E9", "#CC79A7",
           "#D55E00", "#0072B2", "#F0E442", "#999999")
  # Assign colors to states in their original row order, not sorted order
  state_colors <- setNames(pal[seq_len(n_states)], df$state)

  plot_rows <- lapply(metrics, function(m) {
    data.frame(state  = st_ord,
               metric = labels[m],
               value  = df[[m]],
               stringsAsFactors = FALSE,
               row.names = NULL)
  })
  plot_df <- do.call(rbind, plot_rows)
  plot_df$metric <- factor(plot_df$metric, levels = labels[metrics])

  base_plot <- function(d, ttl) {
    ggplot2::ggplot(d, ggplot2::aes(x = state, y = value, fill = state)) +
      ggplot2::geom_col(show.legend = TRUE) +
      ggplot2::scale_fill_manual(values = state_colors, name = NULL) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = NULL, title = ttl) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title       = ggplot2::element_text(face = "bold"),
        legend.position  = "bottom"
      )
  }
  if (!combined) {
    levs <- levels(plot_df$metric)
    plots <- lapply(levs, function(lv) {
      sub <- plot_df[plot_df$metric == lv, , drop = FALSE]
      base_plot(sub, lv)
    })
    names(plots) <- as.character(metrics)
    return(invisible(plots))
  }
  base_plot(plot_df, "Markov Stability Metrics") +
    ggplot2::facet_wrap(~ metric, scales = "free_x", ncol = 2)
}
