#' @title Construct sparse incidence matrix from list of sets
#'
#' @description Given a named list of sets, construct a sparse incidence matrix.
#'
#' @param x a named list of sets. Elements must be of type \code{"character"}.
#'
#' @returns An object of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'   with unique set names as rows and unique elements as columns. A value of 1
#'   indicates the element is a member of the set.
#'
#' @seealso \code{\link{incidence_to_list}}, \code{\link{similarity}},
#'   \code{\link[Matrix]{sparseMatrix}}
#'
#' @import Matrix
#'
#' @export incidence
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z", "z"), # duplicates
#'           "D" = c("a", NA)) # missing values
#'
#' (mat <- incidence(x))
#'
#' # Sizes of sets and their intersections
#' tcrossprod(mat)
#'
#' # Number of elements unique to each set
#' keep <- apply(mat, 2, sum) == 1
#' apply(mat[, keep], 1, sum)

incidence <- function(x) {
  # Validate x, remove missing values, remove duplicate set-element pairs
  dt <- .prepare_sets(x)
  sets <- dt[["sets"]]
  elements <- dt[["elements"]]

  # Unique sets and elements
  sets_unique <- unique(sets)
  elements_unique <- unique(elements)

  # Convert identifiers to positions of 1's
  i <- match(sets, sets_unique)         # row idx
  j <- match(elements, elements_unique) # column idx

  mat_dimnames <- list(sets_unique, elements_unique)

  # Note: sparseMatrix converts integers to numeric
  mat <- sparseMatrix(i = i, j = j, x = 1,
                      dimnames = mat_dimnames)

  return(mat)
}

