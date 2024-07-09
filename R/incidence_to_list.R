#' @title Convert Incidence Matrix to a Named List of Sets
#'
#' @description `r lifecycle::badge("experimental")`
#'
#'   Converts an incidence matrix to a named list of sets. The inverse of
#'   \code{\link{incidence}}.
#'
#' @param incidence incidence matrix with set names as rows and elements as
#'   columns. Usually, the output of \code{\link{incidence}}.
#'
#' @returns a named list of sets with the same length as \code{nrow(incidence)}.
#'
#' @note Currently, there are no checks to ensure \code{incidence} is a valid
#'   incidence matrix.
#'
#' @export incidence_to_list
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z", "z"), # duplicates
#'           "D" = c("a", NA)) # missing values
#'
#' (imat <- incidence(x)) # incidence matrix
#'
#' incidence_to_list(imat)

incidence_to_list <- function(incidence) {
  # The incidence matrix is assumed to be in the correct form with sets as rows
  # and elements as columns.
  idx <- which(incidence == 1, arr.ind = TRUE, useNames = FALSE)

  # Convert indices to names of sets and elements
  elements <- colnames(incidence)[idx[, 2L]]
  sets <- rownames(incidence)[idx[, 1L]]

  # Convert to factor to prevent ordering by set name when splitting
  sets <- factor(sets, levels = unique(sets))

  x <- split(x = elements, f = sets)

  return(x)
}
