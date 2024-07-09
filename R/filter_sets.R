#' @title Filter a Named List of Sets by Size
#'
#' @description Given a named list of sets, filter to those that contain at
#'   least \code{min_size} and no more than \code{max_size} elements. The sets
#'   are optionally restricted to elements of \code{background} before filtering
#'   by size.
#'
#' @inheritParams incidence
#' @param background character; optional character vector. \code{x} will be
#'   filtered to only those elements of \code{background}. If \code{background}
#'   is any other atomic vector, it will be coerced to type \code{character}.
#' @param min_size integer (\eqn{\geq 1}); the minimum allowable set size.
#' @param max_size integer (\eqn{\geq 1}); the maximum allowable set size.
#'
#' @returns A named list of sets at most the same size as \code{x}.
#'
#' @export filter_sets
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("a", "a", "d", "e", NA), # duplicates and NA
#'           "C" = c("f", "g"))
#'
#' # All sets have at least 2 elements, so this just removes duplicates and NA
#' filter_sets(x, min_size = 2L)
#'
#' # Limit scope of sets before filtering
#' filter_sets(x, min_size = 2L, background = c("a", "c", "e", "z"))

filter_sets <- function(x,
                        background,
                        min_size = 5L,
                        max_size = Inf)
{
  # Relax the type of min_size and max_size to allow doubles
  if (!is.vector(min_size, mode = "numeric") || length(min_size) != 1L)
    stop("`min_size` must be a single integer.")

  if (!is.vector(max_size, mode = "numeric") || length(max_size) != 1L)
    stop("`max_size` must be a single integer or Inf.")

  if (min_size > max_size)
    stop("`min_size` cannot be greater than `max_size`.")

  # Sets must contain at least 1 element
  min_size <- max(1L, floor(min_size))
  max_size <- max(1L, floor(max_size))

  # Convert list of sets to a data.table with columns "sets" and "elements"
  dt <- .prepare_sets(x)
  elements <- dt[["elements"]]
  sets <- dt[["sets"]]

  # Validate background
  if (!missing(background)) {
    if (!is.atomic(background))
      stop("If provided, `background` must be an atomic vector, ",
           "preferably of type \"character\".")

    background <- unique(as.character(background))
    background <- background[!is.na(background)]

    if (length(background) == 0L)
      stop("`background` must contain at least 1 unique, nonmissing element.")

    # Subset to elements in background
    in_bg <- which(elements %in% background)

    if (length(in_bg) == 0L)
      stop("No elements of `x` are present in `background`.")

    if (length(in_bg) != length(elements)) {
      elements <- elements[in_bg]
      sets <- sets[in_bg]
    }
  }

  sets <- factor(sets, levels = unique(sets))

  x <- split(x = elements, f = sets)
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

