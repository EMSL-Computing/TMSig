#' @title Set Enrichment Bubble Heatmap
#'
#' @description Create a bubble heatmap summarizing molecular signature analysis
#'   results, such as those from \code{\link{cameraPR.matrix}}. May also be used
#'   to generate bubble heatmaps of differential analysis results.
#'
#' @param x an object that can be coerced to a \code{data.table} with columns
#'   \code{contrast_column}, \code{set_column}, \code{statistic_column}, and
#'   \code{padj_column}.
#' @param n_top integer; number of sets (rows) to display. Defaults to the top
#'   15 sets with the highest median \eqn{-log_{10}}(adjusted p-values) across
#'   contrasts.
#' @param set_column character; the name of a column in \code{x} containing
#'   unique set identifiers that will be used as the row names in the heatmap.
#'   Default is "GeneSet".
#' @param statistic_column character; the name of a column in \code{x}
#'   containing the statistic for each combination of contrast and molecular
#'   signature. Determines the heatmap body colors.
#' @param contrast_column character; the name of a column in \code{x} containing
#'   contrasts that will be used as columns for the heatmap. Entries of
#'   \code{x[[rownames_colum]]} must be uniquely defined for each contrast
#'   group.
#' @param padj_column character; the name of a column in \code{x} containing the
#'   adjusted p-values. Determines the diameter of each bubble in the heatmap.
#' @param padj_legend_title character; title of the background fill legend.
#'   Defaults to \code{padj_column}.
#' @param padj_aggregate_fun function; a function used to aggregate the adjusted
#'   p-values in \code{x[[pvalue_column]]} across contrasts for each unique
#'   entry in \code{x[[set_column]]}. The default computes the median of the
#'   \eqn{-log_{10}} adjusted p-values.
#' @param padj_cutoff numeric; cutoff for terms to be statistically significant.
#'   If \code{plot_sig_only=TRUE}, only those molecular signatures with at least
#'   one \code{padj_column} value less than this threshold may appear in the
#'   heatmap. Default is 0.05.
#' @param plot_sig_only logical; whether to plot only those \code{n_top} terms
#'   that have at least one \code{padj_column} value less than
#'   \code{padj_cutoff}.
#' @param padj_fill character; the background color used for values in
#'   \code{padj_column} that are less than \code{padj_cutoff}. Default is
#'   "grey".
#' @param colors character; vector of length 2 specifying the colors for the
#'   largest negative and largest positive values of
#'   \code{x[[statistic_column]]}, respectively. Default is "#3366ff" (blue) and
#'   "darkred".
#' @param heatmap_color_fun function; used to create the legend for the heatmap
#'   bubble fill. See \code{\link{enrichmap_color_functions}} for details.
#' @param scale_by character; whether to scale the bubbles such that the term
#'   with the largest \eqn{-log_{10}} adjusted p-value in each row
#'   (\code{scale_by="row"}), column (\code{scale_by="column"}), or overall
#'   (\code{scale_by="max"}) is of maximum diameter. Default is "row" to better
#'   visualize patterns across contrasts. May be abbreviated.
#' @param cell_size \code{\link[grid]{unit}} object; the size of each heatmap
#'   cell (used for both height and width). Default is \code{unit(14,
#'   "points")}. This also controls the default text size, which defaults to 90%
#'   the size of \code{cell_size}.
#' @param filename character; the file name used to save the heatmap. If missing
#'   (default), the heatmap will be displayed instead.
#' @param height numeric; height of the file in \code{units}.
#' @param width numeric; width of the file in \code{units}.
#' @param units character; units that define \code{height} and \code{width}.
#'   Defaults to "in" (inches). See \code{\link[grid]{unit}} for possible units.
#' @param heatmap_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#' @param padj_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Legend}}. Modifies the adjusted p-value legend.
#' @param save_args list; additional arguments passed to the graphics device
#'   determined by the \code{filename} extension. See
#'   \code{\link[grDevices]{png}} and \code{\link[grDevices]{pdf}} for options.
#' @param draw_args list; additional arguments passed to
#'   \code{\link[ComplexHeatmap]{draw-HeatmapList-method}}.
#'
#' @returns Nothing. Displays heatmap or saves the heatmap to a file (if
#'   \code{filename} is provided).
#'
#' @details The diameter of each bubble is determined by the \eqn{-log_{10}}
#'   adjusted p-values. By default, the bubbles are scaled such that the
#'   contrast with the largest \eqn{-log_{10}} adjusted p-value per row
#'   (\code{scale_by="row"}) has a bubble diameter of \code{0.95 * cell_size},
#'   and all other bubbles in that row are scaled relative to this maximum
#'   diameter; this is to better visualize patterns across contrasts. Bubbles
#'   can also be scaled so that largest \eqn{-log_{10}} adjusted p-value by
#'   column (\code{scale_by="column"}) or in the entire heatmap
#'   (\code{scale_by="max"}) has the maximum diameter. If the adjusted p-value
#'   is below \code{padj_cutoff}, the bubble diameter will be no smaller than
#'   \code{0.2 * cell_size}. If the adjusted p-value is greater than or equal to
#'   \code{padj_cutoff}, there is no limit on how small the diameter can be, so
#'   it may seem as if the bubbles have disappeared, leaving a blank cell.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap Legend draw
#' @importFrom grDevices dev.off
#' @importFrom grid gpar unit is.unit
#' @importFrom stats median
#' @importFrom utils modifyList
#'
#' @export enrichmap
#'
#' @seealso \code{\link[ComplexHeatmap]{ComplexHeatmap-package}}
#'
#' @examples
#' ## Simulate results of cameraPR.matrix
#' set.seed(1)
#' df <- 5000L
#' x <- data.frame(
#'     Contrast = rep(paste("Contrast", 1:3), each = 4),
#'     GeneSet = rep(paste("GeneSet", 1:4), times = 3),
#'     TwoSampleT = 5 * rt(n = 12L, df = df)
#' )
#'
#' # Calculate z-statistics, two-sided p-values, and BH adjusted p-values
#' x$ZScore <- limma::zscoreT(x = x$TwoSampleT, df = df)
#' x$PValue <- 2 * pnorm(abs(x$ZScore), lower.tail = FALSE)
#' x$FDR <- p.adjust(x$PValue, method = "BH")
#'
#' ## Plot results
#' # Same as enrichmap(x, statistic_column = "ZScore")
#' enrichmap(x = x,
#'           set_column = "GeneSet",
#'           statistic_column = "ZScore",
#'           contrast_column = "Contrast",
#'           padj_column = "FDR",
#'           padj_cutoff = 0.05)
#'
#' # Include gene sets with adjusted p-values above padj_cutoff (0.05). Also
#' # update adjusted p-value legend title.
#' enrichmap(x = x,
#'           statistic_column = "ZScore",
#'           plot_sig_only = FALSE,
#'           padj_legend_title = "BH Adjusted\nP-Value")

enrichmap <- function(x,
                      n_top = 15L,
                      set_column = "GeneSet",
                      statistic_column = "TwoSampleT",
                      contrast_column = "Contrast",
                      padj_column = "FDR",
                      padj_legend_title = padj_column,
                      padj_aggregate_fun = function(padj) {
                        median(-log10(padj), na.rm = TRUE)
                      },
                      padj_cutoff = 0.05,
                      plot_sig_only = TRUE,
                      padj_fill = "grey",
                      colors = c("#3366ff", "darkred"),
                      heatmap_color_fun = camera_color_fun,
                      scale_by = c("row", "column", "max"),
                      cell_size = unit(14, "points"),
                      filename,
                      height = 5,
                      width = 5,
                      units = "in",
                      heatmap_args = list(),
                      padj_args = list(),
                      save_args = list(),
                      draw_args = list())
{
  scale_by <- match.arg(scale_by,
                        choices = c("row", "column", "max"))

  # Create statistic_mat, padj_mat
  ls <- .enrichmap_prepare_x(
    x = x,
    n_top = n_top,
    set_column = set_column,
    statistic_column = statistic_column,
    contrast_column = contrast_column,
    padj_column = padj_column,
    padj_aggregate_fun = padj_aggregate_fun,
    padj_cutoff = padj_cutoff,
    plot_sig_only = plot_sig_only
  )

  statistic_mat <- ls[["statistic_mat"]]
  padj_mat <- ls[["padj_mat"]]

  colorRamp2_args <- heatmap_color_fun(statistic_mat, colors)

  ## Create heatmap ------------------------------------------------------------
  if (!is.unit(cell_size))
    stop("`cell_size` must be a unit object.")

  fontsize <- 0.9 * cell_size # default font size for all components

  # Arguments that will be passed to ComplexHeatmap::Heatmap
  base_heatmap_args <- list(
    matrix = statistic_mat,
    col = do.call(what = circlize::colorRamp2,
                  args = colorRamp2_args),
    name = statistic_column,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = fontsize, fontface = "bold"),
      at = colorRamp2_args$breaks,
      border = "black",
      legend_height = max(cell_size * 5, unit(21.1, "mm")),
      grid_width = max(cell_size, unit(4, "mm"))
    ),
    border = TRUE,
    row_labels = rownames(statistic_mat),
    column_labels = colnames(statistic_mat),
    row_names_gp = gpar(fontsize = fontsize),
    column_names_gp = gpar(fontsize = fontsize),
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = ifelse(anyNA(x), "average", "complete"),
    na_col = "black",
    height = cell_size * nrow(statistic_mat),
    width = cell_size * ncol(statistic_mat),
    layer_fun = .layer_fun
  )

  heatmap_args <- .update_heatmap_args(base_heatmap_args = base_heatmap_args,
                                       heatmap_args = heatmap_args)

  # Color function for bubbles and statistic legend (used by .layer_fun)
  col_fun <- heatmap_args$col

  # Allow layer_fun to access all objects created before this point.
  if (!is.null(heatmap_args[["layer_fun"]]))
    environment(heatmap_args$layer_fun) <- environment()

  # Create heatmap
  ht <- do.call(what = Heatmap, args = heatmap_args)

  # Legend for background fill ----
  padj_legend_labels <- paste(c("<", "\u2265"), padj_cutoff) # \u2265 = ">="
  padj_legend_fill <- c(padj_fill, heatmap_args[["rect_gp"]][["fill"]])

  if (anyNA(statistic_mat)) {
    padj_legend_labels[3] <- "N/A"
    padj_legend_fill[3] <- heatmap_args[["na_col"]]
  }

  base_lt_args <- list(
    at = seq_along(padj_legend_labels),
    title = padj_legend_title,
    title_gp = gpar(fontsize = fontsize, fontface = "bold"),
    labels = padj_legend_labels, # \u2265 = ">="
    labels_gp = gpar(fontsize = fontsize),
    legend_gp = gpar(fill = padj_legend_fill),
    grid_height = heatmap_args$heatmap_legend_param$grid_width,
    grid_width = heatmap_args$heatmap_legend_param$grid_width,
    border = "black",
    nrow = length(padj_legend_labels),
    direction = "horizontal"
  )
  lt_args <- modifyList(x = base_lt_args, val = padj_args, keep.null = TRUE)

  lt <- do.call(what = Legend, args = lt_args)

  if (!missing(filename)) {
    on.exit(dev.off()) # close the device after writing the heatmap to a file

    base_save_args <- list(filename = filename,
                           height = height, width = width,
                           units = units)
    save_args <- modifyList(x = base_save_args,
                            val = save_args,
                            keep.null = TRUE)

    do.call(what = .save_heatmap, args = save_args)
  }

  # Draw heatmap ----
  draw_args <- modifyList(
    x = list(object = ht,
             annotation_legend_list = list(lt),
             merge_legends = TRUE,
             legend_gap = unit(0.15, "in"),
             align_heatmap_legend = "heatmap_top",
             align_annotation_legend = "heatmap_top"),
    val = draw_args,
    keep.null = TRUE
  )

  do.call(what = draw, args = draw_args)
}
