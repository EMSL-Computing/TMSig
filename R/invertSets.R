#' @title Invert a List of Sets, Transposing Sets and Elements
#'
#' @description Invert a list of sets so that elements become set names and set
#'   names become elements.
#'
#' @inheritParams sparseIncidence
#'
#' @note This function is essentially a more limited version of
#'   \code{purrr::transpose_list}.
#'
#' @returns A named list of sets.
#'
#' @export invertSets
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z"),
#'           "D" = c("a", "c", "d"))
#'
#' # Invert sets
#' (y <- invertSets(x))
#'
#' # Jaccard similarity of pairs of elements
#' similarity(y)
#'
#' # Decompose sets into disjoint parts
#' yc <- lapply(y, paste, collapse = ", ")
#' invertSets(yc)

invertSets <- function(x) {
    dt <- .prepare_sets(x)

    split(x = dt[["sets"]], f = dt[["elements"]])
}
