#' @title Construct sparse incidence matrix from list of sets
#'
#' @description Given a named list of sets, construct a sparse incidence matrix.
#'
#' @param x a named list of sets.
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
#' # Number of unique elements in each set
#' keep <- apply(mat, 2, sum) == 1
#' apply(mat[, keep], 1, sum)

incidence <- function(x) {
  if (!is.list(x) | is.null(names(x)))
    stop("`x` must be a named list of character vectors.")

  sets <- rep(names(x), lengths(x))

  # All genes (may include duplicates from the same set)
  elements <- unlist(x, use.names = FALSE)

  # Remove missing elements and convert type to character
  keep <- !is.na(elements)
  elements <- elements[keep]

  if (length(elements) == 0L)
    stop("All sets in `x` are empty or only contain missing values.")

  if (!is.character(elements))
    stop("`x` must be a named list of character vectors.")

  sets <- sets[keep]

  # Unique genes and sets
  elements_unique <- unique(elements)
  sets_unique <- unique(sets)

  # Convert identifiers to positions of nonzero values
  i <- match(sets, sets_unique)         # row idx
  j <- match(elements, elements_unique) # column idx

  mat.dimnames <- list(sets_unique, elements_unique)

  mat <- sparseMatrix(i = i, j = j, x = 1,
                      dimnames = mat.dimnames,
                      use.last.ij = TRUE) # remove duplicates

  return(mat)
}

