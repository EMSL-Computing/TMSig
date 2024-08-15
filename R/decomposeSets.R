#' @title Decompose Pairs of Overlapping Sets Into 3 Disjoint Parts
#'
#' @description Decompose all pairs of sufficiently overlapping sets into 3
#'   disjoint parts: the elements unique to the first set, the elements unique
#'   to the second set, and the elements found in both sets. See the examples
#'   section in \code{\link{invertSets}} for a method to decompose an entire
#'   list of sets.
#'
#' @inheritParams incidence
#' @param overlap integer; only pairs of sets with at least \code{overlap}
#'   elements in common will be decomposed.
#' @param AND character; string used to denote the intersection of two sets.
#'   Defaut is "~AND~", which produces intersections of the form "A ~AND~ B"
#'   (i.e., elements in both A and B).
#' @param MINUS character; string used to denote the difference of two sets.
#'   Defualt is "~MINUS~", which produces differences of the form "A ~MINUS~ B"
#'   (i.e., elements in A and not in B).
#' @param verbose logical; whether to print warnings and messages.
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
#' @seealso \code{\link{filterSets}}
#'
#' @references
#'
#' Jiang, Z., & Gentleman, R. (2007). Extensions to gene set enrichment.
#' \emph{Bioinformatics, 23}(3), 306–313.
#' doi:\href{https://doi.org/10.1093/bioinformatics/btl599
#' }{10.1093/bioinformatics/btl599}
#'
#' @importFrom Matrix tril
#'
#' @export decomposeSets
#'
#' @examples
#' x <- list("A" = letters[1:10],
#'           "B" = letters[3:7],
#'           "C" = letters[1:4],
#'           "D" = letters[6:12])
#'
#' decomposeSets(x)
#'
#' decomposeSets(x, overlap = 5L)

decomposeSets <- function(x,
                          overlap = 1L,
                          AND = "~AND~",
                          MINUS = "~MINUS~",
                          verbose = TRUE) {
    if (!is.character(AND) || !is.character(MINUS))
        stop("`AND` and `MINUS` must be character strings.")

    AND <- paste0(AND, " ") # include trailing space for later
    MINUS <- paste0(MINUS, " ")

    if (!is.vector(overlap, mode = "numeric") ||
        isTRUE(is.infinite(overlap)) || length(overlap) != 1L)
        stop("`overlap` must be a single integer specifying the minimum ",
             "intersection size required to decompose pairs of sets.")

    x <- filterSets(x, min_size = overlap)

    incidence <- incidence(x)
    elements <- colnames(incidence)

    # Sparse lower triangular matrix of intersection sizes
    imat <- tcrossprod(incidence)
    imat <- tril(imat, k = -1L)

    # Indices of sufficiently overlapping pairs of sets
    idx <- which(imat >= overlap, arr.ind = TRUE, useNames = FALSE)
    n_pairs <- nrow(idx)

    if (n_pairs == 0L)
        stop("No pairs of sets with at least `overlap` elements in common.")

    if (verbose)
        message("Decomposing ", n_pairs, " pairs of sets.")

    # Convert indices to set names
    set_pairs <- array(rownames(imat)[idx], dim = dim(idx))

    ## Set decomposition ----
    outcomes <- rep(list(c(MINUS, AND)), 2L) # 1 = in set; 0 = not in set
    outcomes <- expand.grid(outcomes)
    outcomes <- as.matrix(outcomes)[-1, ] # convert to matrix, remove null set

    # Coefficients used to convert each disjoint component from binary to int.
    coefs <- matrix(seq_len(2L), nrow = 1L) # matrix: 2, 1

    x_decomp <- vector(mode = "list", length = n_pairs)
    for (i in seq_len(n_pairs)) { # much faster than apply(set_pairs, 1, ...)
        sets_i <- set_pairs[i, ]
        i_vec <- as.vector(coefs %*% incidence[sets_i, ])

        # Remove elements not in either set
        keep_elements <- which(i_vec != 0L)
        decomp_i <- split(x = elements[keep_elements], f = i_vec[keep_elements])

        ## For a pair of sets A and B, define names of disjoint components:
        # "A ~MINUS~ B" = elements in A and not in B,
        # "B ~MINUS~ A" = elements in B and not in A,
        # "A ~AND~ B" = elements in A and B
        outcomes_i <- outcomes[as.integer(names(decomp_i)), , drop = FALSE]
        names(decomp_i) <- apply(outcomes_i, 1, function(outcome_i) {
            paste(sort(paste0(outcome_i, sets_i)), collapse = " ")
        })

        x_decomp[[i]] <- decomp_i
    }

    x_decomp <- unlist(x_decomp, recursive = FALSE)

    # Remove leading `AND` from names
    names(x_decomp) <- sub(AND, "", names(x_decomp), fixed = TRUE)

    return(x_decomp)
}
