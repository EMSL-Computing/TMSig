#' @title Color Functions for enrichmap
#'
#' @description Heatmap color functions for plotting GSEA-like and CAMERA-like
#'   results with \code{\link{enrichmap}}.
#'
#' @param statistics numeric matrix or vector of statistics. Missing values are
#'   removed. Used to compute limits for the color legend.
#' @param colors vector of 2 colors for the most negative and most positive
#'   values of \code{statistics}. Default is "#3366ff" (blue) and "darkred".
#'
#' @returns a named list of breaks and colors for the heatmap legend.
#'
#' @details For \code{gseaColorFun}, the \code{statistics} are expected to be
#'   normalized enrichment scores (NES). Due to how the NES is formulated,
#'   values between -1 and 1 are never significant or otherwise interesting, so
#'   they are given a white fill so as to not appear in the heatmap (see
#'   examples).
#'
#' @seealso \code{\link{extendRangeNum}}
#'
#' @name enrichmapColorFunctions
#'
#' @examples
#' set.seed(0)
#' x <- rnorm(10, mean = 0, sd = 2)
#'
#' cameraColorFun(x)
#' cameraColorFun(x[x >= 0]) # positive only
#'
#' gseaColorFun(x)
#' gseaColorFun(x[x >= 0]) # positive only


#' @rdname enrichmapColorFunctions
#' @export gseaColorFun
gseaColorFun <- function(statistics,
                         colors = c("#3366ff", "darkred"))
{
    r <- range(statistics, na.rm = TRUE)

    # Extend range of values out to the nearest tenth
    extended_range <- extendRangeNum(r, nearest = 0.1)

    # Due to how the NES is formulated, NES between -1 and 1 are never
    # significant or otherwise interesting, so they are given a white color so
    # as to not appear in the heatmap.
    if (all(r >= -1)) {
        # All NES are effectively positive
        breaks <- c(-1, 1, extended_range[2])
        colors <- c("white", "white", colors[2])
    } else if (all(r <= +1)) {
        # All NES are effectively negative
        breaks <- c(extended_range[1], -1, 1)
        colors <- c(colors[1], "white", "white")
    } else {
        # Mix of positive and negative NES
        breaks <- c(extended_range[1], -1, 1, extended_range[2])
        colors <- c(colors[1], rep("white", 2), colors[2])
    }

    return(list(breaks = breaks, colors = colors))
}



#' @rdname enrichmapColorFunctions
#' @export cameraColorFun
cameraColorFun <- function(statistics,
                           colors = c("#3366ff", "darkred"))
{
    r <- range(statistics, na.rm = TRUE)

    # Extend range of values out to the nearest tenth
    extended_range <- extendRangeNum(r, nearest = 0.1)

    if (all(r >= 0)) {
        breaks <- c(0, extended_range[2])
        colors <- c("white", colors[2])
    } else if (all(r < 0)) {
        breaks <- c(extended_range[1], 0)
        colors <- c(colors[1], "white")
    } else {
        breaks <- c(extended_range[1], 0, extended_range[2])
        colors <- c(colors[1], "white", colors[2])
    }

    return(list(breaks = breaks, colors = colors))
}
