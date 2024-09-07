#' @title Construct a Sparse Incidence Matrix
#'
#' @description Construct a sparse incidence matrix from a named list of sets.
#'
#' @param x a named list of sets. Elements must be of type \code{"character"}.
#' @param ... additional arguments are not currently used.
#'
#' @section Incidence Matrix:
#'
#'   An incidence matrix, \eqn{A}, is defined such that
#'
#'   \deqn{A_{ij} = \begin{cases} 1, & \text{if } \text{element } e_j \in
#'   \text{set } s_i \\ 0, & \text{otherwise} \end{cases}}
#'
#' @details \code{sparseIncidence} differs from
#'   \code{\link[GSEABase:incidence]{GSEABase::incidence}} in that it returns a
#'   sparse matrix, rather than a dense matrix, so it is more memory efficient.
#'   It also removes missing elements, removes empty sets, and combines sets
#'   with duplicate names.
#'
#' @returns An object of class \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'   with unique set names as rows and unique elements as columns.
#'
#' @seealso \code{\link{incidenceToList}}, \code{\link{similarity}},
#'   \code{\link[Matrix]{sparseMatrix}}
#'
#' @examples
#' x <- list("A" = c("a", "b"),
#'           "A" = c("c"), # duplicate sets
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z", "z"), # duplicates
#'           "D" = c("a", NA)) # missing values
#'
#' (imat <- sparseIncidence(x))
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
#'
#' @importFrom Matrix sparseMatrix

sparseIncidence <- function(x, ...) {
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


#' @importFrom methods setGeneric setMethod signature
#'
#' @export
setGeneric(name = "sparseIncidence", def = sparseIncidence)


#' @rdname sparseIncidence
#'
#' @importFrom GSEABase geneIds setName setIdentifier
#'
#' @export
setMethod(f = "sparseIncidence",
          signature = signature(x = "GeneSet"),
          definition = function(x, ...) {
              dots <- names(list(...))
              if (length(dots))
                  warning("Extra arguments disregarded: ", sQuote(dots))

              y <- list(geneIds(x))
              names(y) <- setName(x)

              if (is.null(names(y)))
                  names(y) <- setIdentifier(x)

              sparseIncidence(x = y)
          })


#' @rdname sparseIncidence
#'
#' @importFrom GSEABase geneIds
#'
#' @export
setMethod(f = "sparseIncidence",
          signature = signature(x = "GeneSetCollection"),
          definition = function(x, ...) {
              dots <- names(list(...))
              if (length(dots))
                  warning("Extra arguments disregarded: ", sQuote(dots))

              x <- geneIds(x)

              sparseIncidence(x = x)
          })
