#' @title Filter a list of sets by size
#'
#' @description Given a named list of sets, filter to those that contain at
#'   least \code{min_size} and no more than \code{min_size} elements. The sets
#'   are optionally filtered to elements of \code{background} before filtering
#'   by size.
#'
#' @inheritParams incidence
#' @param background character; optional vector of at least 2 elements. \code{x}
#'   will be filtered to only those elements of \code{background}.
#' @param min_size integer (\eqn{\geq 2}); the minimum set size.
#' @param max_size integer (\eqn{\geq 2}); the maximum set size.
#'
#' @details
#'
#' Elements of \code{x} and \code{background} will be coerced to type
#' \code{character}.
#'
#' If either \code{min_size} or \code{max_size} are less than 2, they will be
#' set to 2.
#'
#' @return A named list of sets at most the same size as \code{x}.
#'
#' @seealso \code{\link{find_aliases}}
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
  # Validate input
  if (!is.list(x) | is.null(names(x)))
    stop("`x` must be a named list.")

  # Relax the type of min_size and max_size to allow doubles
  if (!(typeof(min_size) %in% c("double", "integer")) | length(min_size) != 1L)
    stop("`min_size` must be a single integer.")

  if (!(typeof(max_size) %in% c("double", "integer")) | length(max_size) != 1L)
    stop("`min_size` must be a single integer or Inf.")

  if (min_size > max_size)
    stop("`min_size` cannot be greater than `max_size`.")

  # Sets must contain at least 2 elements
  min_size <- max(2L, min_size)
  max_size <- max(2L, max_size)

  # Convert list of sets to a data.frame
  sets <- rep(names(x), lengths(x))

  if (length(sets) == 0L)
    stop("`x` only contains empty sets.")

  elements <- as.character(unlist(x, use.names = FALSE))

  df <- data.frame(sets = sets,
                   elements = elements,
                   row.names = NULL,
                   stringsAsFactors = FALSE)
  df <- unique(df)
  df <- df[!is.na(df$elements), ]

  if (!missing(background)) {
    # Validate background
    if (!typeof(background) %in% c("character", "double", "integer"))
      stop("`background` must be a character vector containing ",
           "2 or more elements.")

    background <- as.character(unique(background))
    background <- background[!is.na(background)]

    if (length(background) < 2L)
      stop("`background` must contain at least 2 unique, nonmissing elements.")

    # Subset df to elements of background
    keep <- df$elements %in% background

    if (all(!keep))
      stop("No elements of `x` are present in `background`.")

    df <- df[keep, ]
  }

  o <- order(unique(df$sets)) # prevent ordering by set name
  x <- split(x = df$elements, f = df$sets)[o]
  s <- lengths(x)

  if (all(s < min_size))
    stop("No sets contain at least `min_size` elements.")

  if (all(s > max_size))
    stop("All sets contain more than `max.size` elements.")

  keep_set <- (s >= min_size) & (s <= max_size)

  if (!any(keep_set))
    stop("No sets satisfy both `min_size` and `max_size` thresholds.")

  x <- x[keep_set]

  return(x)
}

