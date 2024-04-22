#' @title Convert incidence matrix to a list of sets
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
#' @note This function does not currently check that the input is valid.
#'
#' @export incidence_to_list
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z", "z"), # duplicates
#'           "D" = c("a", NA)) # missing values
#'
#' (mat <- incidence(x)) # incidence matrix
#'
#' incidence_to_list(mat)

incidence_to_list <- function(incidence) {
  # The incidence matrix is assumed to be in the correct form with sets as rows
  # and elements as columns.
  idx <- which(incidence == 1, arr.ind = TRUE, useNames = FALSE)

  # Convert indices to names of sets and elements
  dt <- data.table(sets = rownames(incidence)[idx[, 1L]],
                   elements = colnames(incidence)[idx[, 2L]],
                   stringsAsFactors = FALSE)

  # Convert to factor to prevent ordering by set name when splitting
  dt[, sets := factor(sets, levels = unique(sets))]

  x <- split(x = dt[["elements"]], f = dt[["sets"]])

  return(x)
}
