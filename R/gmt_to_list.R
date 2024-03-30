#' @title Convert GMT file to named list
#'
#' @description
#' Create a named list of sets from a GMT file.
#'
#' @param path character; path to a GMT file.
#' @param check logical; check that \code{path} points to a valid GMT file.
#'
#' @returns A named list of character vectors.
#'
#' @export gmt_to_list
#'
#' @note Similar to \code{fgsea::gmtPathways}.
#'
#' @examples
#' path <- system.file("extdata", "c5.go.v2023.2.Hs.symbols.gmt.gz",
#'                     package = "ostRich")
#'
#' x <- gmt_to_list(path)
#'
#' head(names(x)) # First 6 gene set names
#'
#' x[1] # first set

gmt_to_list <- function(path, check = TRUE) {
  if (check & !grepl("\\.gmt(\\..+)?$", path))
    stop("`path` is not a path to a GMT file.")

  gmt <- readLines(path)
  gmt <- strsplit(gmt, split = "\t")

  out <- lapply(gmt, function(x) x[-c(1, 2)])
  names(out) <- vapply(gmt, function(x) x[1], character(1L))

  return(out)
}
