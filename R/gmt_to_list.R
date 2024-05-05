#' @title Convert GMT file to named list of sets
#'
#' @description Create a named list of sets from a GMT file, or a file
#'   structured like a GMT.
#'
#' @param path character; path to a GMT file. Files may include one additional
#'   extension after ".gmt", such as ".gmt.gzip".
#' @param check logical; check that \code{path} points to a valid GMT file. If
#'   \code{FALSE}, files with different extensions may be read, so long as they
#'   are in the expected format.
#'
#' @returns A named list of character vectors.
#'
#' @export gmt_to_list
#'
#' @note Similar to \code{fgsea::gmtPathways}.
#'
#' @examples
#' path <- system.file("extdata", "c5.go.v2023.2.Hs.symbols.gmt.gz",
#'                     package = "most")
#'
#' x <- gmt_to_list(path)
#'
#' head(names(x)) # First 6 gene set names
#'
#' x[1] # first set

gmt_to_list <- function(path, check = TRUE) {
  if (check & !grepl("\\.gmt(\\.[^\\.]+)?$", path))
    stop("`path` is not a path to a GMT file.")

  gmt <- readLines(path)
  gmt <- strsplit(gmt, split = "\t")

  # This assumes that the first entry of each line is the set name and the
  # elements begin at the third entry.
  out <- lapply(gmt, function(x) x[-c(1, 2)])
  names(out) <- vapply(gmt, FUN = function(x) x[1],
                       FUN.VALUE = character(1L))

  return(out)
}
