#' @title Extend the Range of Values Out to the Nearest Digit
#'
#' @description Extend the range of values out to the nearest digit.
#'
#' @param x any \code{\link[base]{numeric}} object.
#' @param nearest numeric; the range of \code{x} will be extended out to the
#'   value specified by \code{nearest}. Default is 1, which extends the range
#'   out to the nearest integer.
#'
#' @returns A \code{numeric} vector of length 2 containing the minimum and
#'   maximum values of \code{x} after extending them outward to the value
#'   provided by \code{nearest}.
#'
#' @seealso \code{\link[base]{range}}, \code{\link[grDevices]{extendrange}}
#'
#' @export range_extend
#'
#' @examples
#' set.seed(0)
#' x <- runif(5, min = -10, max = 10)
#' range(x) # -4.689827  8.164156
#'
#' range_extend(x) # -5  9
#' range_extend(x, nearest = 2) # -6  10
#' range_extend(x, nearest = 0.1) # -4.7  8.2

range_extend <- function(x, nearest = 1) {
  r <- range(x, na.rm = TRUE) / nearest

  c(floor(r[1]), ceiling(r[2])) * nearest
}
