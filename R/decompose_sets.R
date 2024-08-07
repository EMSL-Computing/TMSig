#' @title Decompose Pairs of Overlapping Sets Into 3 Disjoint Parts
#'
#' @description `r lifecycle::badge("experimental")`
#'
#'   Decompose all pairs of sufficiently overlapping sets into 3 disjoint parts:
#'   the elements unique to the first set, the elements unique to the second
#'   set, and the elements found in both sets. See the examples section in
#'   \code{\link{invert_sets}} for a method to decompose an entire list of sets.
#'
#' @inheritParams incidence
#' @param overlap integer; only pairs of sets with at least \code{overlap}
#'   elements in common will be decomposed.
#'
#' @section Optimization:
#'
#'   Since the size of the intersection between two sets is at most the size of
#'   the smaller set, any sets with fewer than \code{overlap} elements can be
#'   immediately discarded.
#'
#' @source Decomposition of sets is described by Jiang and Gentleman (2007) in
#'   section 2.3.1 "Overlap among gene sets". It is a method to reduce the
#'   redundancy of significant gene set testing results whereby the decomposed
#'   sets are reanalyzed and the following selections can be made:
#'
#' \itemize{
#'   \item{If the elements unique to set 1 and set 2, elements common to both
#'   sets, or all 3 parts are statistically significant, keep both set 1 and set
#'   2 in the original results. We can not separate their effects.}
#'
#'   \item{If the elements unique to set 1 or the elements unique to set 1 and
#'   common to both sets are statistically significant, only keep set 1 in the
#'   original results. (The same logic can be applied for set 2.)}
#' }
#'
#' @returns A named list of disjoint parts of sets. May contain aliases.
#'
#' @seealso \code{\link{filter_sets}}
#'
#' @references
#'
#' Jiang, Z., & Gentleman, R. (2007). Extensions to gene set enrichment.
#' \emph{Bioinformatics, 23}(3), 306–313.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btl599
#' }{10.1093/bioinformatics/btl599}
#'
#' @import Matrix
#' @importFrom data.table data.table rbindlist `:=`
#'
#' @export decompose_sets
#'
#' @examples
#' x <- list("A" = letters[1:10],
#'           "B" = letters[3:7],
#'           "C" = letters[1:4],
#'           "D" = letters[6:12])
#'
#' decompose_sets(x)
#'
#' decompose_sets(x, overlap = 5L)

decompose_sets <- function(x, overlap = 1L)
{
  lifecycle::signal_stage("experimental", "decompose_sets()")

  if (!is.vector(overlap, mode = "numeric") ||
      isTRUE(is.infinite(overlap)) ||
      length(overlap) != 1L)
    stop("`overlap` must be a single integer specifying the minimum ",
         "intersection size required to decompose pairs of sets.")

  overlap <- max(1L, floor(overlap))

  # This also validates x. Since the size of the intersection between a set and
  # any other set is at most the size of that set, we can pre-filter to sets of
  # size `overlap` or greater.
  x <- filter_sets(x, min_size = overlap)

  if (length(x) < 2L)
    stop("Fewer than 2 sets with at least `overlap` elements.")

  incidence <- incidence(x)

  # Sparse lower triangular matrix of intersection sizes
  imat <- tcrossprod(incidence)
  imat <- tril(imat, k = -1L)

  # Indices of sufficiently overlapping pairs of sets
  idx <- which(imat >= overlap, arr.ind = TRUE, useNames = FALSE)

  if (nrow(idx) == 0L)
    stop("Fewer than 2 sets with at least `overlap` elements in common.")

  # Flip columns to preserve order of sets in x
  idx <- idx[, 2:1, drop = FALSE]

  # Convert indices to set names
  set_pairs <- array(rownames(imat)[idx], dim = dim(idx))

  # Set decomposition ----
  n <- 2L # only pairs of sets are currently supported

  # Coefficients used to convert each disjoint component from binary to int.
  coefs <- matrix(2L ^ seq(n - 1L, 0L), nrow = 1L)
  # e.g., n = 5 --> (16, 8, 4, 2, 1)

  outcomes <- rep(list(0:1), n) # 1 = in set; 0 = not in set
  outcomes <- expand.grid(outcomes)
  outcomes <- as.matrix(outcomes)[-1, ] # convert to matrix, remove null set

  outcome_vec <- as.vector(tcrossprod(coefs, outcomes))

  elements <- colnames(incidence)

  x_decomp <- apply(set_pairs, 1, function(sets_i) {
    ## For a pair of sets A and B, define names of disjoint components:
    # "A NOT B" = elements in A and not in B,
    # "B NOT A" = elements in B and not in A,
    # "A AND B" = elements in A and B
    decomp_labels <- apply(outcomes, 1, function(in_set) {
      label_i <- paste(sort(paste0(ifelse(in_set, "AND ", "NOT "), sets_i)),
                       collapse = " ")
      sub("^AND ", "", label_i)
    })

    # Convert incidence matrix columns from binary to integer
    i_vec <- as.vector(coefs %*% incidence[sets_i, ])

    idx <- match(i_vec, outcome_vec)
    # Reorder indices to preserve order of decomp_labels
    o <- order(idx)
    idx <- idx[o]
    not_null <- which(!is.na(idx)) # NA = null set

    decomp_i <- data.table(set = decomp_labels[idx],
                           element = elements[o],
                           stringsAsFactors = FALSE)[not_null, ]

    return(decomp_i)
  })

  x_decomp <- rbindlist(x_decomp)
  x_decomp[, set := factor(set, levels = unique(set))]

  # Convert data.table to list. set was converted to a factor to maintain the
  # adjacency of disjoint components when using split().
  x_decomp <- split(x = x_decomp[["element"]], f = x_decomp[["set"]])

  return(x_decomp)
}
