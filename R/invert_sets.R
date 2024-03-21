#' @title Invert a list of sets, transposing positions of sets and elements
#'
#' @description `r lifecycle::badge("experimental")`
#'
#'   Invert a list of sets so that elements become set names and set names
#'   become elements.
#'
#' @inheritParams incidence
#'
#' @note This function is essentially a simplified version of
#'   \code{purrr::transpose_list}.
#'
#' @returns A named list of sets.
#'
#' @import Matrix
#'
#' @export invert_sets
#'
#' @examples
#' x <- list("A" = c("a", "b", "c"),
#'           "B" = c("c", "d"),
#'           "C" = c("x", "y", "z"),
#'           "D" = c("a", "c", "d"))
#'
#' # Invert sets
#' (y <- invert_sets(x))
#'
#' # Decompose sets into disjoint parts
#' x_dc <- lapply(y, paste, collapse = ", ")
#' (x_dc <- invert_sets(x_dc))

invert_sets <- function(x) {
  ## Parsimonious approach
  # imat <- incidence(x)
  # imat <- t(imat)
  #
  # y <- incidence_to_list(imat)

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

  # The code above is from incidence()
  # TODO wrap reused code like this into an unexported function.

  y <- split(x = sets, f = elements)

  return(y)
}
