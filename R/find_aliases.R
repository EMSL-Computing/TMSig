#' @title List aliased sets
#'
#' @description Given a list of sets, collapse groups of aliased sets to the
#'   first set by alphanumeric order and list the rest as aliases.
#'
#' @inheritParams incidence
#'
#' @section Function Details:
#'
#'   \code{find_aliases} uses the fact that a Jaccard similarity coefficient of
#'   1 indicates that sets are aliased. Groups of aliased sets are identified
#'   with \code{\link{cluster_sets}}. For each group, the first set
#'   alphanumerically is taken as the set to keep, while all others are listed
#'   as aliases. All sets are pre-filtered with \code{\link{filter_sets}} to
#'   have at least 2 unique elements.
#'
#' @returns A list of aliases with primary sets as names. If no sets are
#'   aliased, a warning will be issued and the function will return a list with
#'   the same length as \code{x}, but all elements will be \code{character(0)}.
#'
#' @seealso \code{\link{filter_sets}}, \code{\link{similarity}},
#'   \code{\link{cluster_sets}}
#'
#' @examples
#' ## B, C, and D are aliased
#' x <- list("A" = letters[1:5],
#'           "B" = letters[1:4],
#'           "C" = letters[1:4],
#'           "D" = letters[1:4])
#'
#' # Names are sets to keep. List elements are aliases
#' (aliases <- find_aliases(x))
#'
#' # Remove aliases from x
#' (x <- x[names(aliases)])
#'
#' ## No sets are aliased
#' y <- list("A" = letters[1:3],
#'           "B" = letters[2:3])
#'
#' find_aliases(y)
#'
#' @importFrom data.table setDT dcast
#' @importFrom stats setNames
#'
#' @export find_aliases

find_aliases <- function(x) {
  # Suppress "No pair of sets passes the similarity cutoff." message when no
  # sets are aliased.
  suppressMessages(
    d <- cluster_sets(x, method = "jaccard", cutoff = 1, h = 0)
  )
  d$type <- ifelse(duplicated(d$cluster), "alias", "set")
  setDT(d)

  d_wide <- dcast(d, formula = cluster ~ type,
                  fun.aggregate = list,
                  value.var = "set")

  if (!"alias" %in% colnames(d_wide)) {
    message("No sets are aliased.")
    d_wide[["alias"]] <- list(character(0L))
  }

  with(d_wide, setNames(object = alias, nm = set))
}
