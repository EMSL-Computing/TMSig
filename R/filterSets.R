#' @title Filter a Named List of Sets by Size
#'
#' @description Given a named list of sets, filter to those that contain at
#'   least \code{min_size} and no more than \code{max_size} elements. The sets
#'   are optionally restricted to elements of \code{background} before filtering
#'   by size.
#'
#' @inheritParams sparseIncidence
#' @param background character; optional character vector. \code{x} will be
#'   filtered to only those elements of \code{background}.
#' @param min_size integer (\eqn{\geq 1}); the minimum allowable set size.
#' @param max_size integer (\eqn{\geq 1}); the maximum allowable set size.
#'
#' @returns A named list of sets at most the same size as \code{x}.
#'
#' @importFrom data.table :=
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("a", "a", "d", "e", NA), # duplicates and NA
#'           "C" = c("f", "g"))
#'
#' # All sets have at least 2 elements,
#' # so this just removes duplicates and NA
#' filterSets(x, min_size = 2L)
#'
#' # Limit scope of sets before filtering
#' filterSets(x, min_size = 2L, background = c("a", "c", "e", "z"))

filterSets <- function(x,
                       background = NULL,
                       min_size = 5L,
                       max_size = Inf,
                       ...) {
    if (!is.vector(min_size, mode = "numeric") || length(min_size) != 1L)
        stop("`min_size` must be a single integer.")

    if (!is.vector(max_size, mode = "numeric") || length(max_size) != 1L)
        stop("`max_size` must be a single integer or Inf.")

    if (min_size > max_size)
        stop("`min_size` cannot be greater than `max_size`.")

    min_size <- pmax(1, floor(min_size))
    max_size <- pmax(1, floor(max_size))

    # Convert list of sets to a data.table with columns "sets" and "elements"
    dt <- .prepare_sets(x = x, background = background)
    dt[, sets := factor(sets, levels = unique(sets))] # preserve set order

    x <- split(x = dt[["elements"]], f = dt[["sets"]])
    ss <- lengths(x) # set sizes

    if (all(ss < min_size))
        stop("No sets contain at least `min_size` elements.")

    if (all(ss > max_size))
        stop("All sets contain more than `max_size` elements.")

    keep_set <- which(ss >= min_size & ss <= max_size)

    if (length(keep_set) == 0L)
        stop("No sets satisfy both `min_size` and `max_size` thresholds.")

    x <- x[keep_set]

    return(x)
}


#' @importFrom methods setGeneric setMethod signature
#'
#' @export
setGeneric(name = "filterSets", def = filterSets)


#' @rdname filterSets
#'
#' @importFrom GSEABase geneIds setName setIdentifier
#'
#' @export
setMethod(f = "filterSets",
          signature = signature(x = "GeneSet",
                                background = "character",
                                min_size = "numeric",
                                max_size = "numeric"),
          definition = function(x,
                                background = NULL,
                                min_size = 5L,
                                max_size = Inf,
                                ...) {
              dots <- names(list(...))
              if (length(dots))
                  warning("Extra arguments disregarded: ", sQuote(dots))

              y <- list(geneIds(x))
              names(y) <- setName(x)

              if (is.null(names(y)))
                  names(y) <- setIdentifier(x)

              filterSets(x = y, background = background,
                         min_size = min_size, max_size = max_size)
          })


#' @rdname filterSets
#'
#' @importFrom GSEABase geneIds
#'
#' @export
setMethod(f = "filterSets",
          signature = signature(x = "GeneSetCollection",
                                background = "character",
                                min_size = "numeric",
                                max_size = "numeric"),
          definition = function(x,
                                background = NULL,
                                min_size = 5L,
                                max_size = Inf,
                                ...) {
              dots <- names(list(...))
              if (length(dots))
                  warning("Extra arguments disregarded: ", sQuote(dots))

              x <- geneIds(x)

              filterSets(x = x, background = background,
                         min_size = min_size, max_size = max_size)
          })
