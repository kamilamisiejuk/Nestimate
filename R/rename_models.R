# ---- rename_models ----
#
# A `netobject_group` is a named list whose names label each constituent
# network ("Cluster 1", "High", "A", ...). Those labels show up in every
# downstream report (print methods, plot panel titles, bootstrap output,
# permutation tests, etc.), so the group built by `build_network()` may
# need post-hoc renaming once the analyst knows what each cluster
# actually represents.
#
# `rename_models()` is the equivalent of `tna::rename_groups()`. The
# Nestimate name uses "models" because the Nestimate vocabulary refers
# to each constituent of a `netobject_group` as a model (one fitted
# network), reserving "group" for the dispatch class.

#' Rename the models of a `netobject_group`
#'
#' Replaces the names of the constituent networks in a `netobject_group`
#' (or any object inheriting from it). Useful when `build_network()`
#' produced generic labels (e.g. `"Cluster 1"`, `"Cluster 2"`) and you
#' want to substitute meaningful ones (e.g. `"High engagement"`,
#' `"Low engagement"`).
#'
#' @param x A `netobject_group` (or any object inheriting from it, such
#'   as `net_mlvar`).
#' @param new_names A character vector of new names. Must have the same
#'   length as `x`, contain no `NA` or empty strings, and be unique.
#' @return A `netobject_group` of the same class with renamed members.
#' @examples
#' \dontrun{
#'   d   <- tna::group_regulation
#'   grp <- build_network(d, method = "tna",
#'                        group = sample(c("a", "b"), nrow(d), TRUE))
#'   grp <- rename_models(grp, c("High", "Low"))
#'   names(grp)
#' }
#' @export
rename_models <- function(x, new_names) {
  UseMethod("rename_models")
}

#' @rdname rename_models
#' @export
rename_models.netobject_group <- function(x, new_names) {
  stopifnot(
    is.character(new_names),
    length(new_names) == length(x),
    !anyNA(new_names),
    all(nzchar(new_names)),
    !anyDuplicated(new_names)
  )
  names(x) <- new_names
  x
}

#' @rdname rename_models
#' @export
rename_models.default <- function(x, new_names) {
  stop("rename_models() expects a netobject_group (or subclass). Got: ",
       paste(class(x), collapse = ", "), ".", call. = FALSE)
}
