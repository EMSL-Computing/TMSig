#' @title Convert Incidence Matrix to a Named List of Sets
#'
#' @description Converts an incidence matrix to a named list of sets. The
#'   inverse of \code{\link{sparseIncidence}}.
#'
#' @param incidence incidence matrix with set names as rows and elements as
#'   columns. For instance, the output of \code{\link{sparseIncidence}}.
#'
#' @returns a named list of sets with the same length as \code{nrow(incidence)}.
#'
#' @note Currently, there are no checks to ensure \code{incidence} is a valid
#'   incidence matrix.
#'
#' @importFrom Matrix which
#'
#' @export incidenceToList
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z", "z"), # duplicates
#'           "D" = c("a", NA)) # missing values
#'
#' (imat <- sparseIncidence(x)) # incidence matrix
#'
#' incidenceToList(incidence = imat)

incidenceToList <- function(incidence) {
    # The incidence matrix is assumed to be in the correct form with sets as
    # rows and elements as columns.
    idx <- which(incidence == 1, arr.ind = TRUE, useNames = FALSE)

    # Convert indices to names of sets and elements
    elements <- colnames(incidence)[idx[, 2L]]
    sets <- rownames(incidence)[idx[, 1L]]

    # Convert to factor to prevent ordering by set name when splitting
    sets <- factor(sets, levels = unique(sets))

    x <- split(x = elements, f = sets)

    return(x)
}
