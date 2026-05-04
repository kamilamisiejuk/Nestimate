# ---- Association Rule Mining ----

#' Discover Association Rules from Sequential or Transaction Data
#'
#' @description
#' Discovers association rules using the Apriori algorithm with proper
#' candidate pruning. Accepts \code{netobject} (extracts sequences as
#' transactions), data frames, lists, or binary matrices.
#'
#' Support counting is vectorized via \code{crossprod()} for 2-itemsets
#' and logical matrix indexing for k-itemsets.
#'
#' @param x Input data. Accepts:
#'   \describe{
#'     \item{netobject}{Uses \code{$data} sequences — each sequence becomes
#'       a transaction of its unique states.}
#'     \item{list}{Each element is a character vector of items (one transaction).}
#'     \item{data.frame}{Wide format: each row is a transaction, character
#'       columns are item occurrences. Or a binary matrix (0/1).}
#'     \item{matrix}{Binary transaction matrix (rows = transactions,
#'       columns = items).}
#'   }
#' @param min_support Numeric. Minimum support threshold. Default: 0.1.
#' @param min_confidence Numeric. Minimum confidence threshold. Default: 0.5.
#' @param min_lift Numeric. Minimum lift threshold. Default: 1.0.
#' @param max_length Integer. Maximum itemset size. Default: 5.
#'
#' @return An object of class \code{"net_association_rules"} containing:
#' \describe{
#'   \item{rules}{Data frame with columns: antecedent (list), consequent (list),
#'     support, confidence, lift, conviction, count, n_transactions.}
#'   \item{frequent_itemsets}{List of frequent itemsets per level k.}
#'   \item{items}{Character vector of all items.}
#'   \item{n_transactions}{Integer.}
#'   \item{n_rules}{Integer.}
#'   \item{params}{List of min_support, min_confidence, min_lift, max_length.}
#' }
#'
#' @details
#' ## Algorithm
#'
#' Uses level-wise Apriori (Agrawal & Srikant, 1994) with the full pruning
#' step: after the join step generates k-candidates, all (k-1)-subsets are
#' verified as frequent before support counting. This is critical for
#' efficiency at k >= 4.
#'
#' ## Metrics
#'
#' \describe{
#'   \item{support}{P(A and B). Fraction of transactions containing both
#'     antecedent and consequent.}
#'   \item{confidence}{P(B | A). Fraction of antecedent transactions that
#'     also contain the consequent.}
#'   \item{lift}{P(A and B) / (P(A) * P(B)). Values > 1 indicate positive
#'     association; < 1 indicate negative association.}
#'   \item{conviction}{(1 - P(B)) / (1 - confidence). Measures departure
#'     from independence. Higher = stronger implication.}
#' }
#'
#' @references
#' Agrawal, R. & Srikant, R. (1994). Fast algorithms for mining association
#' rules. In \emph{Proc. 20th VLDB Conference}, 487--499.
#'
#' @examples
#' # From a list of transactions
#' trans <- list(
#'   c("plan", "discuss", "execute"),
#'   c("plan", "research", "analyze"),
#'   c("discuss", "execute", "reflect"),
#'   c("plan", "discuss", "execute", "reflect"),
#'   c("research", "analyze", "reflect")
#' )
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5)
#' print(rules)
#'
#' # From a netobject (sequences as transactions)
#' seqs <- data.frame(
#'   V1 = sample(LETTERS[1:5], 50, TRUE),
#'   V2 = sample(LETTERS[1:5], 50, TRUE),
#'   V3 = sample(LETTERS[1:5], 50, TRUE)
#' )
#' net <- build_network(seqs, method = "relative")
#' rules <- association_rules(net, min_support = 0.1)
#'
#' @seealso \code{\link{build_network}}, \code{\link{predict_links}}
#'
#' @export
association_rules <- function(x,
                              min_support = 0.1,
                              min_confidence = 0.5,
                              min_lift = 1.0,
                              max_length = 5L) {

  # ---- Input validation ----
  stopifnot(
    is.numeric(min_support), length(min_support) == 1,
    min_support > 0, min_support <= 1,
    is.numeric(min_confidence), length(min_confidence) == 1,
    min_confidence >= 0, min_confidence <= 1,
    is.numeric(min_lift), length(min_lift) == 1, min_lift >= 0,
    is.numeric(max_length), length(max_length) == 1, max_length >= 2
  )
  max_length <- as.integer(max_length)

  # ---- Parse input to binary transaction matrix ----
  parsed <- .ar_parse_input(x)
  trans_mat <- parsed$matrix  # logical matrix: rows=transactions, cols=items
  items <- parsed$items
  n_trans <- nrow(trans_mat)
  n_items <- length(items)
  min_count <- ceiling(min_support * n_trans)

  # ---- Level 1: frequent 1-itemsets ----
  item_counts <- colSums(trans_mat)
  freq1_mask <- item_counts >= min_count
  freq1_items <- items[freq1_mask]
  if (length(freq1_items) == 0) {
    return(.ar_empty_result(items, n_trans, min_support, min_confidence,
                            min_lift, max_length))
  }
  # Prune transaction matrix to frequent items only
  trans_mat <- trans_mat[, freq1_items, drop = FALSE]
  items <- freq1_items
  n_items <- length(items)
  item_counts <- item_counts[freq1_mask]

  frequent_itemsets <- list()
  frequent_itemsets[[1]] <- lapply(seq_along(items), function(i) {
    list(items = items[i], count = item_counts[i],
         support = item_counts[i] / n_trans)
  })

  # ---- Level 2: frequent 2-itemsets via crossprod ----
  if (max_length >= 2 && n_items >= 2) {
    co_counts <- as.matrix(crossprod(trans_mat * 1L))
    # Upper triangle only (avoid duplicates)
    freq2 <- list()
    for (i in seq_len(n_items - 1L)) {
      for (j in seq(i + 1L, n_items)) {
        if (co_counts[i, j] >= min_count) {
          freq2 <- c(freq2, list(list(
            items = c(items[i], items[j]),
            count = co_counts[i, j],
            support = co_counts[i, j] / n_trans
          )))
        }
      }
    }
    if (length(freq2) > 0) frequent_itemsets[[2]] <- freq2
  }

  # ---- Level k >= 3: join + prune + count ----
  if (max_length >= 3) {
    k <- 3L
    while (k <= max_length) {
      if (length(frequent_itemsets) < k - 1L) break
      prev <- frequent_itemsets[[k - 1L]]
      if (is.null(prev) || length(prev) < 2L) break

      prev_items <- lapply(prev, `[[`, "items")
      # Join step: merge (k-1)-itemsets sharing first k-2 items
      candidates <- .ar_generate_candidates(prev_items, k)
      if (length(candidates) == 0) break

      # Prune step: verify all (k-1)-subsets are frequent
      prev_set <- vapply(prev_items, paste, character(1), collapse = "\t")
      candidates <- Filter(function(cand) {
        subsets <- combn(cand, k - 1L, simplify = FALSE)
        all(vapply(subsets, function(s) paste(s, collapse = "\t"), character(1))
            %in% prev_set)
      }, candidates)
      if (length(candidates) == 0) break

      # Count support via logical AND
      freq_k <- list()
      for (cand in candidates) {
        mask <- Reduce(`&`, lapply(cand, function(it) trans_mat[, it]))
        cnt <- sum(mask)
        if (cnt >= min_count) {
          freq_k <- c(freq_k, list(list(
            items = cand, count = cnt, support = cnt / n_trans
          )))
        }
      }
      if (length(freq_k) == 0) break
      frequent_itemsets[[k]] <- freq_k
      k <- k + 1L
    }
  }

  # ---- Generate rules from frequent itemsets ----
  rules <- .ar_generate_rules(frequent_itemsets, item_counts, items,
                               n_trans, min_confidence, min_lift)

  # ---- Build tidy rules data frame ----
  if (nrow(rules) > 0) {
    tidy_rules <- data.frame(
      antecedent = vapply(rules$antecedent, paste, character(1), collapse = ", "),
      consequent = vapply(rules$consequent, paste, character(1), collapse = ", "),
      support    = rules$support,
      confidence = rules$confidence,
      lift       = rules$lift,
      conviction = rules$conviction,
      count      = rules$count,
      n_transactions = rules$n_transactions,
      stringsAsFactors = FALSE
    )
  } else {
    tidy_rules <- data.frame(
      antecedent = character(0), consequent = character(0),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    )
  }

  # ---- Build tidy frequent itemsets data frame ----
  fi_rows <- lapply(seq_along(frequent_itemsets), function(k) {
    level <- frequent_itemsets[[k]]
    data.frame(
      itemset = vapply(level, function(fi) paste(fi$items, collapse = ", "), character(1)),
      size    = k,
      support = vapply(level, `[[`, numeric(1), "support"),
      count   = vapply(level, `[[`, numeric(1), "count"),
      stringsAsFactors = FALSE
    )
  })
  frequent <- do.call(rbind, fi_rows)
  if (is.null(frequent)) {
    frequent <- data.frame(itemset = character(0), size = integer(0),
                           support = numeric(0), count = numeric(0),
                           stringsAsFactors = FALSE)
  }

  structure(list(
    rules              = tidy_rules,
    frequent           = frequent,
    frequent_itemsets   = frequent_itemsets,
    items              = items,
    n_transactions     = n_trans,
    n_rules            = nrow(tidy_rules),
    params             = list(min_support = min_support,
                              min_confidence = min_confidence,
                              min_lift = min_lift,
                              max_length = max_length)
  ), class = "net_association_rules")
}


# ---- Input Parsing ----

#' @noRd
.ar_parse_input <- function(x) {
  # mcml → netobject_group → use first element
  if (inherits(x, "mcml")) x <- as_tna(x)
  if (inherits(x, "netobject_group")) x <- x[[1]]
  if (inherits(x, "cograph_network") && !inherits(x, "netobject")) {
    x <- .as_netobject(x)
  }

  # tna object: extract $data, decode integer labels if needed
  if (inherits(x, "tna")) {
    if (!is.null(x$data)) {
      df <- as.data.frame(x$data, stringsAsFactors = FALSE)
      labels <- rownames(x$weights)
      # Decode integer-encoded data to state names
      if (!is.null(labels) && length(labels) > 0 &&
          (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
        df[] <- lapply(df, function(col) {
          idx <- as.integer(col)
          ifelse(is.na(idx) | idx < 1L | idx > length(labels),
                 NA_character_, labels[idx])
        })
      }
      transactions <- lapply(seq_len(nrow(df)), function(i) {
        vals <- as.character(unlist(df[i, ], use.names = FALSE))
        unique(vals[!is.na(vals) & vals != ""])
      })
      return(.ar_transactions_to_matrix(transactions))
    }
    stop("tna object has no $data. Build from sequence data.", call. = FALSE)
  }

  # netobject: sequences → transactions (unique states per row)
  if (inherits(x, "netobject")) {
    if (is.null(x$data)) {
      stop("netobject has no $data. Build from sequence data.", call. = FALSE)
    }
    df <- as.data.frame(x$data, stringsAsFactors = FALSE)
    transactions <- lapply(seq_len(nrow(df)), function(i) {
      vals <- unlist(df[i, ], use.names = FALSE)
      unique(vals[!is.na(vals) & vals != ""])
    })
    return(.ar_transactions_to_matrix(transactions))
  }

  # List of character vectors
  if (is.list(x) && !is.data.frame(x)) {
    stopifnot(all(vapply(x, is.character, logical(1))))
    return(.ar_transactions_to_matrix(x))
  }

  # Matrix (assumed binary)
  if (is.matrix(x)) {
    items <- colnames(x) %||% paste0("I", seq_len(ncol(x)))
    colnames(x) <- items
    ord <- order(items)
    return(list(matrix = (x > 0)[, ord, drop = FALSE], items = items[ord]))
  }

  # Data frame
  if (is.data.frame(x)) {
    # Check if binary (0/1 numeric)
    all_numeric <- all(vapply(x, is.numeric, logical(1)))
    if (all_numeric && all(unlist(x) %in% c(0, 1, NA))) {
      items <- colnames(x)
      ord <- order(items)
      mat <- as.matrix(x) > 0
      return(list(matrix = mat[, ord, drop = FALSE], items = items[ord]))
    }
    # Wide character data: each row is a transaction
    transactions <- lapply(seq_len(nrow(x)), function(i) {
      vals <- unlist(x[i, ], use.names = FALSE)
      vals <- as.character(vals)
      unique(vals[!is.na(vals) & vals != ""])
    })
    return(.ar_transactions_to_matrix(transactions))
  }

  stop("x must be a netobject, list, data.frame, or matrix.", call. = FALSE)
}


#' @noRd
.ar_transactions_to_matrix <- function(transactions) {
  items <- sort(unique(unlist(transactions, use.names = FALSE)))
  n_trans <- length(transactions)
  n_items <- length(items)

  # Vectorized construction
  row_idx <- rep(seq_along(transactions), lengths(transactions))
  col_idx <- match(unlist(transactions, use.names = FALSE), items)

  mat <- matrix(FALSE, n_trans, n_items, dimnames = list(NULL, items))
  mat[cbind(row_idx, col_idx)] <- TRUE

  list(matrix = mat, items = items)
}


# ---- Candidate Generation ----

#' @noRd
.ar_generate_candidates <- function(prev_items, k) {
  # Standard Apriori join: two (k-1)-itemsets sharing first k-2 items
  n <- length(prev_items)
  candidates <- list()
  for (i in seq_len(n - 1L)) {
    a <- prev_items[[i]]
    prefix_a <- a[seq_len(k - 2L)]
    for (j in seq(i + 1L, n)) {
      b <- prev_items[[j]]
      prefix_b <- b[seq_len(k - 2L)]
      if (identical(prefix_a, prefix_b) && a[k - 1L] < b[k - 1L]) {
        candidates <- c(candidates, list(c(a, b[k - 1L])))
      }
    }
  }
  candidates
}


# ---- Rule Generation ----

#' @noRd
.ar_generate_rules <- function(frequent_itemsets, item_counts, items,
                                n_trans, min_confidence, min_lift) {
  # Pre-build support lookup for all frequent itemsets
  support_lookup <- new.env(hash = TRUE, parent = emptyenv())
  for (level in frequent_itemsets) {
    for (fi in level) {
      key <- paste(fi$items, collapse = "\t")
      support_lookup[[key]] <- fi$support
    }
  }

  # Item-level support for lift/conviction
  item_support <- setNames(item_counts / n_trans, items)

  rules_list <- list()

  # Generate rules from itemsets of size >= 2
  if (length(frequent_itemsets) < 2L) {
    return(data.frame(
      antecedent = I(list()), consequent = I(list()),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  for (k in seq(2L, length(frequent_itemsets))) {
    level <- frequent_itemsets[[k]]
    if (is.null(level)) next
    for (fi in level) {
      itemset <- fi$items
      sup_ab <- fi$support
      cnt_ab <- fi$count
      n_items_k <- length(itemset)

      # Generate all non-empty proper subsets as antecedents
      for (size in seq_len(n_items_k - 1L)) {
        subsets <- combn(n_items_k, size, simplify = FALSE)
        for (idx in subsets) {
          ante <- itemset[idx]
          cons <- itemset[-idx]

          # Antecedent support
          ante_key <- paste(ante, collapse = "\t")
          sup_a <- support_lookup[[ante_key]]
          if (is.null(sup_a)) next

          # Confidence = P(A,B) / P(A)
          confidence <- sup_ab / sup_a
          if (confidence < min_confidence) next

          # Consequent support (for lift/conviction)
          cons_key <- paste(cons, collapse = "\t")
          sup_b <- support_lookup[[cons_key]]
          if (is.null(sup_b)) {
            # Single-item consequent
            if (length(cons) == 1L) {
              sup_b <- item_support[cons]
            } else {
              next
            }
          }

          # Lift = P(A,B) / (P(A) * P(B))
          lift <- sup_ab / (sup_a * sup_b)
          if (lift < min_lift) next

          # Conviction = (1 - P(B)) / (1 - confidence)
          conviction <- if (confidence >= 1) Inf
                        else (1 - sup_b) / (1 - confidence)

          rules_list <- c(rules_list, list(data.frame(
            antecedent     = I(list(ante)),
            consequent     = I(list(cons)),
            support        = sup_ab,
            confidence     = confidence,
            lift           = lift,
            conviction     = conviction,
            count          = cnt_ab,
            n_transactions = n_trans,
            stringsAsFactors = FALSE
          )))
        }
      }
    }
  }

  if (length(rules_list) == 0) {
    return(data.frame(
      antecedent = I(list()), consequent = I(list()),
      support = numeric(0), confidence = numeric(0),
      lift = numeric(0), conviction = numeric(0),
      count = integer(0), n_transactions = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  rules <- do.call(rbind, rules_list)
  rules <- rules[order(-rules$lift, -rules$confidence), , drop = FALSE]
  rownames(rules) <- NULL
  rules
}


#' @noRd
.ar_empty_result <- function(items, n_trans, min_support, min_confidence,
                              min_lift, max_length) {
  empty_rules <- data.frame(
    antecedent = character(0), consequent = character(0),
    support = numeric(0), confidence = numeric(0),
    lift = numeric(0), conviction = numeric(0),
    count = integer(0), n_transactions = integer(0),
    stringsAsFactors = FALSE
  )
  empty_freq <- data.frame(
    itemset = character(0), size = integer(0),
    support = numeric(0), count = numeric(0),
    stringsAsFactors = FALSE
  )
  structure(list(
    rules = empty_rules,
    frequent = empty_freq,
    frequent_itemsets = list(),
    items = items,
    n_transactions = n_trans,
    n_rules = 0L,
    params = list(min_support = min_support, min_confidence = min_confidence,
                  min_lift = min_lift, max_length = max_length)
  ), class = "net_association_rules")
}


# ---- S3 Methods ----

#' Print Method for net_association_rules
#'
#' @param x A \code{net_association_rules} object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#'
#' @examples
#' trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5,
#'                            min_lift = 0)
#' print(rules)
#'
#' @export
print.net_association_rules <- function(x, ...) {
  cat(sprintf("Association Rules  [%d rules | %d items | %d transactions]\n",
              x$n_rules, length(x$items), x$n_transactions))
  cat(sprintf("  Support >= %.2f  |  Confidence >= %.2f  |  Lift >= %.2f\n",
              x$params$min_support, x$params$min_confidence, x$params$min_lift))

  if (x$n_rules > 0) {
    top <- utils::head(x$rules, 10)
    cat("\n  Top rules (by lift):\n")
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      cat(sprintf("    %d. %s -> %s  (sup=%.3f conf=%.3f lift=%.2f)\n",
                  i, r$antecedent, r$consequent, r$support, r$confidence, r$lift))
    }
    if (x$n_rules > 10) {
      cat(sprintf("    ... and %d more rules\n", x$n_rules - 10))
    }
  }
  invisible(x)
}


#' Summary Method for net_association_rules
#'
#' @param object A \code{net_association_rules} object.
#' @param ... Additional arguments (ignored).
#' @return A data frame summarizing the rules, invisibly.
#'
#' @examples
#' trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"), c("A","C","D"))
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.5,
#'                            min_lift = 0)
#' summary(rules)
#'
#' @export
summary.net_association_rules <- function(object, ...) {
  r <- object$rules
  row.names(r) <- NULL
  r
}


#' Plot Method for net_association_rules
#'
#' @description
#' Scatter plot of association rules: support vs confidence, with point
#' size proportional to lift.
#'
#' @param x A \code{net_association_rules} object.
#' @param ... Additional arguments passed to \code{ggplot2} functions.
#' @return A \code{ggplot} object, invisibly.
#'
#' @examples
#' trans <- list(c("A","B","C"), c("A","B"), c("B","C","D"),
#'               c("A","C","D"), c("A","B","D"), c("B","C"))
#' rules <- association_rules(trans, min_support = 0.3, min_confidence = 0.3,
#'                            min_lift = 0)
#' plot(rules)
#'
#' @import ggplot2
#' @export
plot.net_association_rules <- function(x, ...) {
  r <- x$rules
  if (nrow(r) == 0) {
    message("No rules to plot.")
    return(invisible(NULL))
  }

  df <- data.frame(
    support    = r$support,
    confidence = r$confidence,
    lift       = r$lift,
    rule       = paste0("{", r$antecedent, "} => {", r$consequent, "}"),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = support, y = confidence,
                                         size = lift, color = lift)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::scale_color_gradient(low = "#2E7D32", high = "#C62828") +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(x = "Support", y = "Confidence",
                  title = sprintf("Association Rules (%d rules)", nrow(r)),
                  size = "Lift", color = "Lift") +
    ggplot2::theme_minimal(base_size = 12)

  print(p)
  invisible(p)
}
