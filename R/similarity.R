#' @title Construct a Matrix of Pairwise Set Similarity Coefficients
#'
#' @description Construct a sparse matrix of Jaccard or overlap similarity
#'   coefficients for each pair of sets in a list.
#'
#' @inheritParams incidence
#' @param type character; the type of similarity measure to use. Either
#'   \code{"jaccard"} or \code{"overlap"}. May be abbreviated.
#'
#' @returns A symmetric \code{\link[Matrix:dgCMatrix-class]{dgCMatrix}}
#'   containing all pairwise set similarity coefficients.
#'
#' @section Set Similarity:
#'
#'   If \eqn{A} and \eqn{B} are sets, we define the Jaccard similarity
#'   coefficient \eqn{J} as the size of their intersection divided by the size
#'   of their union (Jaccard, 1912):
#'
#'   \deqn{J(A, B) = \frac{|A \cap B|}{|A \cup B|}}
#'
#'   The overlap coefficient (also known as the Szymkiewicz–Simpson coefficient)
#'   is defined as the size of the intersection divided by the size of the
#'   smaller set (Szymkiewicz 1934; Simpson, 1943, 1947, 1960; Fallaw 1979):
#'
#'   \deqn{Overlap(A, B) = \frac{|A \cap B|}{min(|A|, |B|)}}
#'
#'   The Jaccard coefficient can identify aliased sets (sets which contain the
#'   same elements, but have different names), while the overlap coefficient can
#'   identify both aliased sets and subsets. Aliases and subsets are not easily
#'   distinguished without also having the matrix of Jaccard coefficients or the
#'   set sizes.
#'
#' @section Optimization:
#'
#'   Calculations are only performed for pairs of sets with nonzero
#'   intersections in the lower triangular part of the matrix. As such,
#'   \code{similarity} is efficient even for large similarity matrices, and it
#'   is especially efficient for sparse similarity matrices.
#'
#' @references Jaccard, P. (1912). The distribution of the flora in the alpine
#'   zone. \emph{The New Phytologist, 11}(2), 37–50.
#'   doi:\href{https://doi.org/10.1111/j.1469-8137.1912.tb05611.x}{10.1111/j.1469-8137.1912.tb05611.x}.
#'   \url{https://www.jstor.org/stable/2427226}
#'
#'   Szymkiewicz, D. (1934). Une contribution statistique à la géographie
#'   floristique. \emph{Acta Societatis Botanicorum Poloniae, 11}(3), 249–265.
#'   doi:\href{https://doi.org/10.5586/asbp.1934.012}{10.5586/asbp.1934.012}.
#'
#'   Simpson, G. G. (1943). Mammals and the nature of continents. \emph{American
#'   Journal of Science, 241}(1), 1–31.
#'
#'   Simpson, G. G. (1947). Holarctic mammalian faunas and continental
#'   relationships during the Cenozoic. \emph{Bulletin of the Geological Society
#'   of America, 58}(7), 613–688.
#'
#'   Simpson, G. G. (1960). Notes on the measurement of faunal resemblance.
#'   \emph{American Journal of Science, 258-A}, 300–311.
#'
#'   Fallaw, W. C. (1979). A Test of the Simpson Coefficient and Other Binary
#'   Coefficients of Faunal Similarity. \emph{Journal of Paleontology, 53}(4),
#'   1029–1034. \url{http://www.jstor.org/stable/1304126}
#'
#' @seealso \code{\link{incidence}}, \code{\link{cluster_sets}}
#'
#' @import Matrix
#'
#' @export similarity
#'
#' @examples
#' x <- list("A" = c("a", "b", "c", "d", "e"),
#'           "B" = c("d", "e", "f", "g"), # overlaps with A
#'           "C" = c("d", "e", "f", "g"), # aliased with B
#'           "D" = c("a", "b", "c")) # subset of A
#'
#' (j <- similarity(x)) # Jaccard coefficients
#' (o <- similarity(x, type = "overlap")) # overlap coefficients

similarity <- function(x,
                       type = c("jaccard", "overlap"))
{
  type <- match.arg(type, choices = c("jaccard", "overlap"))

  # Incidence matrix with set identifiers as rows and elements as columns
  # 1 if the element is a member of the set; 0 otherwise
  incidence_mat <- incidence(x)

  if (nrow(incidence_mat) < 2L)
    stop("`x` must contain 2 or more sets.")

  # Sizes of all pairwise intersections: t(incidence_mat) %*% incidence_mat
  mat <- tcrossprod(incidence_mat)
  set_sizes <- diag(mat)

  # Since mat is symmetric, we only need to perform calculations using the
  # lower triangular part
  mat <- tril(mat, k = -1L)

  # Positions of nonzero intersections
  idx <- which(mat > 0, arr.ind = TRUE, useNames = FALSE)

  # Array of set sizes
  size_array <- array(set_sizes[idx], dim = dim(idx))
  size_intersect <- mat[idx]

  switch(type,
         jaccard = {
           ## |A union B| = |A| + |B| - |A intersect B|
           size_union <- size_array[, 1L] + size_array[, 2L] - size_intersect
           ## |A intersect B| / |A union B|
           sim <- size_intersect / size_union
         },
         overlap = {
           ## |A intersect B| / min(|A|, |B|)
           sim <- size_intersect / pmin(size_array[, 1L], size_array[, 2L])
         })

  # idx only covers the lower triangular part, so we flip the indices to fill
  # in the upper triangular part
  mat[idx] <- mat[idx[, 2:1, drop = FALSE]] <- sim
  diag(mat) <- 1

  return(mat)
}
