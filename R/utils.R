
#' @title Prepare a list of sets for other functions
#'
#' @description Remove missing values, remove duplicate set-element pairs, and
#'   create a \code{data.table}.
#'
#' @inheritParams incidence
#'
#' @returns A \code{data.table} with columns \code{sets} and \code{elements}.
#'
#' @importFrom data.table data.table
#'
#' @noRd

.prepare_sets <- function(x) {
  if (!is.list(x) | is.null(names(x)))
    stop("`x` must be a named list of character vectors.")

  sets <- rep(names(x), lengths(x))

  # All genes (may include duplicates from the same set)
  elements <- unlist(x, use.names = FALSE)

  # Remove missing elements and check type
  keep <- !is.na(elements)
  elements <- elements[keep]

  if (length(elements) == 0L)
    stop("All sets in `x` are empty or only contain missing values.")

  if (!is.character(elements))
    stop("`x` must be a named list of character vectors.")

  sets <- sets[keep]

  dt <- data.table("sets" = sets,
                   "elements" = elements,
                   stringsAsFactors = FALSE)
  dt <- unique(dt) # remove duplicates

  return(dt)
}
