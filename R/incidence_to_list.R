#' @title Convert incidence matrix to a named list of sets
#'
#' @description Converts an incidence matrix to a named list. The inverse of
#'   \code{\link{incidence}}.
#'
#' @param incidence incidence matrix with set names as rows and elements as
#'   columns. Usually, the output of \code{\link{incidence}}.
#'
#' @returns a named list of sets with the same length as \code{nrow(incidence)}.
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
  df <- data.frame(sets = rownames(incidence)[idx[, 1L]],
                   elements = colnames(incidence)[idx[, 2L]],
                   row.names = NULL,
                   stringsAsFactors = FALSE)

  o <- unique(df$sets) # prevent ordering by set name
  x <- split(x = df$elements, f = df$sets)[o]

  return(x)
}
