#' @title Construct Sparse Incidence Matrix from List of Sets
#'
#' @description Given a named list of sets, construct a sparse incidence matrix.
#'
#' @param x a named list of sets. Elements must be of type \code{"character"}.
#'
#' @details The incidence matrix, \eqn{B}, is defined such that
#'
#' \deqn{B_{ij} = \begin{cases} 1, & \text{if } \text{element } e_j \in
#' \text{set } s_i \\ 0, & \text{otherwise} \end{cases}}
#'
#' @returns An object of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'   with unique set names as rows and unique elements as columns.
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
#' (imat <- incidence(x))
#'
#' # Sizes of sets and their pairwise intersections
#' tcrossprod(imat)
#'
#' # Number of sets in which each element and pair of elements appears
#' crossprod(imat)
#'
#' # Count number of elements unique to each set
#' keep <- apply(imat, 2, sum) == 1
#' apply(imat[, keep], 1, sum)

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

