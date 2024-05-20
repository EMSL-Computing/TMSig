#' @title Set Enrichment Bubble Heatmap
#'
#' @description Create a bubble heatmap summarizing set enrichment analysis
#'   results, such as those from \code{\link{cameraPR.matrix}}.
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
#' @param statistic_column similar to \code{set_column}. The name of a column
#'   containing the statistic for each pathway. Determines the heatmap body
#'   colors.
#' @param contrast_column character; the name of a column in \code{x} containing
#'   contrasts that will be used as columns for the heatmap. Entries of
#'   \code{x[[rownames_colum]]} must be uniquely defined for each contrast
#'   group.
#' @param padj_column similar to \code{set_column}. The name of a column
#'   containing the adjusted p-values that determine the area of each circle in
#'   the heatmap.
#' @param padj_legend_title character; title of the background fill legend.
#'   Defaults to \code{padj_column}.
#' @param padj_aggregate_fun function; a function used to aggregate the adjusted
#'   p-values in \code{x[[pvalue_column]]} across contrasts for each unique
#'   entry in \code{x[[set_column]]}. Passed to \code{\link[base]{match.fun}}.
#'   The default computes the median of the \eqn{-log_{10}} adjusted p-values.
#' @param padj_cutoff numeric; cutoff for terms to be statistically significant.
#'   If \code{plot_sig_only=TRUE}, only those pathways with at least one
#'   \code{padj_column} value less than this threshold may appear in the
#'   heatmap. Default is 0.05.
#' @param plot_sig_only logical; whether to plot only those \code{n_top} terms
#'   that have at least one \code{padj_column} value less than
#'   \code{padj_cutoff}.
#' @param padj_fill character; the background color used for values in
#'   \code{padj_column} that are less than \code{padj_cutoff}. Default is
#'   "grey".
#' @param colors vector of length 2 specifying the colors for the largest
#'   negative and largest positive values of \code{x[[statistic_column]]},
#'   respectively. Default is "#3366ff" (blue) and "darkred".
#' @param heatmap_color_fun function; used to create the legend for the heatmap
#'   bubble fill. See \code{\link{enrichmap_color_functions}} for details.
#' @param scale_by character; whether to scale the circles such that the
#'   most-significant term in each row (\code{scale_by="row"}), column
#'   (\code{scale_by="column"}), or overall (\code{scale_by="max"}) is of
#'   maximum area. Default is "row" to better visualize patterns across
#'   contrasts.
#' @param cell_size \code{\link[grid]{unit}} object; the size of each heatmap
#'   cell (used for both height and width). Default is \code{unit(14,
#'   "points")}.
#' @param filename character; the file name used to save the heatmap. If missing
#'   (default), the heatmap will be displayed instead.
#' @param height numeric; height of the file in \code{units}.
#' @param width numeric; width of the file in \code{units}.
#' @param units character; units that define \code{height} and \code{width}.
#'   Defaults to "in" (inches).
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
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap max_text_width Heatmap Legend draw
#' @importFrom data.table as.data.table `:=` `.N` setorderv dcast
#' @importFrom grDevices dev.off
#' @importFrom grid gpar unit
#' @importFrom stats median
#' @importFrom utils modifyList
#'
#' @export enrichmap
#'
#' @seealso \code{\link[ComplexHeatmap]{ComplexHeatmap-package}}

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
  # Check padj_cutoff
  if (padj_cutoff < 0 | padj_cutoff > 1)
    stop("`padj_cutoff` must be between 0 and 1.")

  scale_by <- match.arg(scale_by, scale_by)

  # n_top does not strictly have to be an integer, since it will be rounded
  if (!is.numeric(n_top))
    stop("`n_top` must be an integer.")

  n_top <- max(round(n_top), 1L)

  # Required columns
  cols_to_keep <- c(set_column, padj_column,
                    statistic_column, contrast_column)

  # Check that all required columns are present
  col_present <- cols_to_keep %in% colnames(x)
  if (any(!col_present)) {
    missing_cols <- cols_to_keep[!col_present]
    stop(sprintf("Missing columns: %s", paste(missing_cols, collapse = ", ")))
  }

  x <- as.data.table(x)

  # Identify sets that are significant in at least 1 contrast
  x[, `:=`(any_sig = any(get(padj_column) < padj_cutoff)),
    by = set_column]

  # Filter to significant sets for plotting
  if (plot_sig_only)
    x <- subset(x, subset = any_sig == TRUE)

  if (nrow(x) == 0L)
    stop("No terms are significant at `padj_cutoff`. ",
         "Consider setting plot_sig_only=FALSE.")

  # Select n_top most significant terms
  x[, `:=`(criteria = padj_aggregate_fun(get(padj_column))),
    by = set_column]

  setorderv(x,
            cols = c(contrast_column, "criteria"),
            order = c(1, -1))

  # Avoid using utils::head - conflicts with Matrix::head
  top_pathways <- unique(x[, get(set_column)])
  top_pathways <- top_pathways[seq_len(min(n_top, length(top_pathways)))]

  x <- subset(x, get(set_column) %in% top_pathways)

  x <- unique(x[, cols_to_keep, with = FALSE]) # remove unnecessary columns

  n <- x[, .N, by = set_column][["N"]] # number of entries by row name

  # Multiple set_column entries per contrast
  if (!all(n == max(n)))
    stop("`set_column` is not uniquely defined for each contrast.")

  # Reshape data to wide format and convert to a matrix
  x <- dcast(x,
             formula = get(set_column) ~ get(contrast_column),
             value.var = c(statistic_column, padj_column))
  x <- as.matrix(x, rownames = 1)

  # Matrices of adjusted p-values and set statistics
  padj_mat <- x[, grepl(paste0("^", padj_column), colnames(x)),
                drop = FALSE]
  statistic_mat <- x[, grepl(paste0("^", statistic_column), colnames(x)),
                     drop = FALSE]

  colnames(padj_mat) <- colnames(statistic_mat) <- contrasts <-
    sub(paste0("^", padj_column, "_"), "", colnames(padj_mat))

  colorRamp2_args <- heatmap_color_fun(statistic_mat, colors)

  # Create heatmap -------------------------------------------------------------
  # Arguments that will be passed to ComplexHeatmap::Heatmap
  base_heatmap_args <- list(
    matrix = statistic_mat,
    col = do.call(what = circlize::colorRamp2,
                  args = colorRamp2_args),
    heatmap_legend_param = list(
      title = statistic_column,
      at = colorRamp2_args$breaks,
      border = "black",
      legend_height = max(cell_size * 5, unit(21.1, "mm")),
      grid_width = max(cell_size, unit(4, "mm"))
    ),
    border = TRUE,
    row_names_max_width = max_text_width(rownames(x)),
    column_names_max_height = max_text_width(contrasts),
    cluster_columns = FALSE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = ifelse(anyNA(x), "average", "complete"),
    na_col = "black",
    height = cell_size * nrow(statistic_mat),
    width = cell_size * ncol(statistic_mat),
    layer_fun = .layer_fun
  )

  # Update with user-supplied arguments
  heatmap_args <-  modifyList(x = base_heatmap_args,
                              val = heatmap_args,
                              keep.null = TRUE)

  # Mark missing values
  heatmap_args$rect_gp <- gpar(col = NA, fill = heatmap_args$na_col)

  # Color function for circles and statistic legend
  col_fun <- heatmap_args$col

  # If layer_fun is specified, set the environment to be the environment of
  # enrichmap. This allows it to access all objects created before this point.
  if (!is.null(heatmap_args$layer_fun))
    environment(heatmap_args$layer_fun) <- environment()

  # Create heatmap
  ht <- do.call(what = Heatmap, args = heatmap_args)

  # Legend for background fill ----
  # base args
  lt_args <- list(
    title = padj_legend_title,
    at = 1:2,
    labels = paste(c("<", "\u2265"), padj_cutoff),
    legend_gp = gpar(fill = c(padj_fill, "white")),
    grid_height = heatmap_args$heatmap_legend_param$grid_width,
    grid_width = heatmap_args$heatmap_legend_param$grid_width,
    border = "black",
    nrow = 2,
    direction = "horizontal"
  )
  lt_args <- modifyList(x = lt_args, val = padj_args, keep.null = TRUE)

  lt <- do.call(what = Legend, args = lt_args)

  if (!missing(filename)) {
    on.exit(dev.off())

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
    val = draw_args, keep.null = TRUE
  )

  do.call(what = draw, args = draw_args)
}
