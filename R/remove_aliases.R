#' @title Remove aliased sets
#'
#' @inheritParams incidence
#'
#' @section Function Details:
#'
#'   \code{remove_aliases} is a wrapper for \code{find_aliases}, which itself
#'   uses \code{cluster_sets}. It is recommended to use \code{find_aliases}
#'   first. This way, aliases can be identified and included in the final
#'   analysis results.
#'
#'   All sets are pre-filtered to have at least 2 unique elements. For each
#'   group of aliased sets, the first set is kept according to alphanumeric
#'   order. All non-aliased sets are also kept.
#'
#' @return a named list of sets with aliases removed.
#'
#' @seealso \code{\link{find_aliases}}, \code{\link{filter_sets}}
#'
#' @export remove_aliases
#'
#' @examples
#' # A and B are aliased, C and D are aliased
#' x <- list("C" = 4:6,
#'           "B" = 1:3,
#'           "A" = 1:3,
#'           "D" = 4:6)
#'
#' remove_aliases(x) # B and D removed

remove_aliases <- function(x) {
  # List of aliases. Names are sets to keep
  alias_list <- find_aliases(x)

  x <- x[names(alias_list)]

  return(x)
}
