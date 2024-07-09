#' @title Prepare a List of Sets for Other Functions
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
  elements <- unlist(x, recursive = FALSE, use.names = FALSE)

  if (!is.vector(elements, mode = "character"))
    stop("`x` must be a named list of character vectors.")

  # Remove missing elements and check type
  keep <- which(!is.na(elements))

  if (length(keep) == 0L)
    stop("All sets in `x` are empty or only contain missing values.")

  if (length(keep) != length(elements)) {
    elements <- elements[keep]
    sets <- sets[keep]
  }

  dt <- data.table(sets = sets,
                   elements = elements,
                   stringsAsFactors = FALSE)

  # Remove duplicate set-element pairs
  dt <- unique(dt)

  return(dt)
}


## enrichmap utility functions -------------------------------------------------

#' @title Format the Cells of the Heatmap
#'
#' @param j integer; column index.
#' @param i integer; row index.
#' @param x numeric; x-coordinate of the cell.
#' @param y numeric; y-coordinate of the cell.
#' @param w numeric; width of the cell.
#' @param h numeric; height of the cell.
#' @param f character; the filled color of the cell.
#'
#' @details See \code{\link[ComplexHeatmap]{Heatmap}} \code{layer_fun} parameter
#'   for details. This function has access to all objects in the
#'   \code{enrichmap} function environment that are created before the
#'   environment is set, including \code{statistic_mat}, \code{padj_mat}, and
#'   all function arguments.
#'
#' @importFrom ComplexHeatmap pindex
#' @importFrom grid grid.rect grid.circle gpar
#'
#' @noRd

.layer_fun <- function(j, i, x, y, w, h, f) {
  # Cell background
  grid.rect(x = x, y = y, width = w, height = h,
            gp = gpar(col = NA,
                      fill = ifelse(
                        pindex(padj_mat, i, j) < padj_cutoff, padj_fill,
                        ifelse(is.na(pindex(padj_mat, i, j)),
                               NA_character_, "white")
                      )
            ))

  # Matrix of radii (optionally scaled to row or column max)
  rmat <- -log10(padj_mat)

  if (scale_by != "max") {
    # Scale bubbles relative to row or column max
    margin <- 1L + (scale_by == "column")
    rmat <- sweep(rmat, MARGIN = margin,
                  apply(rmat, MARGIN = margin, max, na.rm = TRUE),
                  FUN = "/")
  } else {
    # Scale bubbles relative to global max
    rmat <- rmat / max(rmat, na.rm = TRUE)
  }

  # Limits on bubble radii for significant adjusted p-values
  r_min <- 0.20
  r_max <- 0.95 # upper limit because a black border is added to the bubbles
  rmat <- ifelse(padj_mat < padj_cutoff,
                 rmat * (r_max - r_min) + r_min,
                 rmat * r_max)

  # Draw bubbles. col_fun is taken from heatmap_color_fun.
  grid.circle(
    x = x, y = y,
    r = pindex(rmat, i, j) / 2 * cell_size,
    # Significant bubbles get a black outline to separate from padj_fill
    gp = gpar(col = ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                           "black", NA),
              fill = col_fun(pindex(statistic_mat, i, j)))
  )
}



#' @title Save the Heatmap to a File
#'
#' @param filename character; path to the file. While a default is provided, it
#'   will never be used by enrichmap. The file extension will be used to
#'   determine the graphics device.
#' @param width numeric; width of the file in \code{units}. Default is 480
#'   (pixels).
#' @param height numeric; height of the file in \code{units}. Default is 480
#'   (pixels).
#' @param units character; units that define \code{height} and \code{width}.
#'   Default is "px" (pixels).
#' @param res numeric; resolution of file (in ppi: pixels per inch). Default is
#'   500.
#' @param ... additional arguments passed to the graphics device.
#'
#' @importFrom grDevices bmp jpeg png tiff pdf
#' @importFrom utils modifyList
#'
#' @noRd

.save_heatmap <- function(filename = "set_bubble_heatmap%03d.png",
                          width = 480,
                          height = 480,
                          units = "px",
                          res = 500,
                          ...)
{
  # see tools::file_ext
  file_ext <- sub(".*\\.([[:alnum:]]+)$", "\\1", filename)

  save_fun <- switch(
    file_ext,
    "bmp" = bmp,
    "jpg" = jpeg,
    "pdf" = pdf,
    "png" = png,
    "tiff" = tiff,
    stop("`filename` must have one of the following extensions: ",
         "bmp, jpg, pdf, png, or tiff.")
  )

  default_args <- list(filename = filename,
                       file = filename,
                       width = width,
                       height = height,
                       units = units,
                       quality = 100,
                       res = res,
                       compression = "lzw")
  save_args <- modifyList(default_args, val = list(...))
  save_args <- save_args[names(save_args) %in% names(formals(save_fun))]

  do.call(what = save_fun, args = save_args)
}
