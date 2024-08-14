#' @title Prepare a List of Sets for Other Functions
#'
#' @description Remove missing values, remove duplicate set-element pairs, and
#'   create a \code{data.table}.
#'
#' @inheritParams incidence
#' @inheritParams filterSets
#'
#' @returns A \code{data.table} with columns \code{sets} and \code{elements}.
#'
#' @importFrom data.table data.table
#'
#' @noRd

.prepare_sets <- function(x, background = NULL) {
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

    if (!is.null(background)) {
        if (!is.atomic(background))
            stop("If provided, `background` must be an atomic vector, ",
                 "preferably of type \"character\".")

        background <- unique(as.character(background))
        background <- background[!is.na(background)]

        if (length(background) == 0L)
            stop("`background` must contain at least 1 unique, ",
                 "nonmissing element.")

        # Subset to elements in background
        in_bg <- which(elements %in% background)

        if (length(in_bg) == 0L)
            stop("No elements of `x` are present in `background`.")

        if (length(in_bg) != length(elements)) {
            elements <- elements[in_bg]
            sets <- sets[in_bg]
        }
    }

    dt <- data.table(sets = sets,
                     elements = elements,
                     stringsAsFactors = FALSE)

    # Remove duplicate set-element pairs
    dt <- unique(dt)

    return(dt)
}


#' @title Validate min_size and max_size for filterSets
#'
#' @inheritParams filterSets
#'
#' @returns A length 2 numeric vector of the set size limits. May include
#'   \code{Inf}.
#'
#' @noRd

.checkSizeRange <- function(min_size = 1L, max_size = Inf) {
    # Relax the type of min_size and max_size to allow doubles
    if (!is.vector(min_size, mode = "numeric") || length(min_size) != 1L)
        stop("`min_size` must be a single integer.")

    if (!is.vector(max_size, mode = "numeric") || length(max_size) != 1L)
        stop("`max_size` must be a single integer or Inf.")

    if (min_size > max_size)
        stop("`min_size` cannot be greater than `max_size`.")

    # Sets must contain at least 1 element
    min_size <- max(1L, floor(min_size))
    max_size <- max(1L, floor(max_size))
    sizeRange <- floor(c(min_size, max_size))
    sizeRange <- pmax(sizeRange, 1)

    return(sizeRange)
}


## enrichmap utility functions -------------------------------------------------

#' @title Prepare data for enrichmap
#'
#' @inheritParams enrichmap
#'
#' @importFrom data.table as.data.table `:=` `.N` setorderv dcast
#'
#' @noRd

.enrichmap_prepare_x <- function(x,
                                 n_top = 15L,
                                 set_column = "GeneSet",
                                 statistic_column = "TwoSampleT",
                                 contrast_column = "Contrast",
                                 padj_column = "FDR",
                                 padj_aggregate_fun = function(padj)
                                     median(-log10(padj), na.rm = TRUE),
                                 padj_cutoff = 0.05,
                                 plot_sig_only = TRUE) {
    if (padj_cutoff < 0 | padj_cutoff > 1)
        stop("`padj_cutoff` must be between 0 and 1.")

    if (!is.vector(n_top, mode = "numeric"))
        stop("`n_top` must be an integer.")

    n_top <- max(floor(n_top), 1L)

    # Required columns
    required_cols <- c(set_column, padj_column,
                       statistic_column, contrast_column)

    x <- as.data.table(x)
    x <- unique(x[, required_cols, with = FALSE, drop = FALSE])

    # Identify sets that are significant in at least 1 contrast
    x[, `:=`(any_sig = any(get(padj_column) < padj_cutoff)),
      by = set_column]

    # Filter to significant sets for plotting
    if (plot_sig_only)
        x <- subset(x, subset = any_sig)

    if (nrow(x) == 0L)
        stop("No terms are significant at `padj_cutoff`. ",
             "Consider setting plot_sig_only=FALSE.")

    # Select n_top most significant terms
    x[, criteria := padj_aggregate_fun(get(padj_column)), by = set_column]

    setorderv(x, cols = c(contrast_column, "criteria"), order = c(1, -1))

    # Filter to top pathways
    top_pathways <- unique(x[[set_column]])
    top_pathways <- top_pathways[seq_len(min(n_top, length(top_pathways)))]
    x <- subset(x, subset = get(set_column) %in% top_pathways)

    # All n should be 1 if there are no duplicates
    n <- x[, .N, by = c(contrast_column, set_column)][["N"]]

    if (any(n != 1L))
        stop("set_column=", sQuote(set_column),
             " is not uniquely defined for each contrast.")

    # Reshape data to wide format and convert to a matrix
    x <- dcast(x,
               formula = get(set_column) ~ get(contrast_column),
               value.var = c(statistic_column, padj_column),
               fill = NA)
    x <- as.matrix(x, rownames = 1)

    # Split into matrices of adjusted p-values and set statistics
    padjColumnIndices <- grep(paste0("^", padj_column, "_"), colnames(x))
    out <- list("statistic_mat" = x[, -padjColumnIndices, drop = FALSE],
                "padj_mat" = x[, padjColumnIndices, drop = FALSE])

    out <- lapply(out, function(mat) {
        colnames(mat) <- sub(paste0("^", padj_column, "_"), "", colnames(mat))
        return(mat)
    })

    return(out)
}


#' @title Update enrichmap heatmap arguments
#'
#' @description Update rect_gp, row_names_max_width, and
#'   column_names_max_height.
#'
#' @param base_heatmap_args list of base heatmap args.
#' @param heatmap_args optional list of user-supplied arguments. Used to update
#'   \code{base_heatmap_args}.
#'
#' @returns A list of arguments that will be passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @importFrom ComplexHeatmap max_text_width
#' @importFrom grid gpar
#' @importFrom utils modifyList
#'
#' @noRd

.update_heatmap_args <- function(base_heatmap_args = list(),
                                 heatmap_args = list()) {
    if (!is.list(heatmap_args))
        stop("`heatmap_args` must be a list of arguments that ",
             "will be passed to ComplexHeatmap::Heatmap.")

    # Update with user-supplied arguments
    heatmap_args <-  modifyList(x = base_heatmap_args,
                                val = heatmap_args,
                                keep.null = TRUE)

    # Default rect_gp
    if (is.null(heatmap_args[["rect_gp"]]))
        heatmap_args[["rect_gp"]] <- gpar(col = NA, fill = "white")

    if (is.null(heatmap_args[["rect_gp"]][["col"]]))
        heatmap_args[["rect_gp"]][["col"]] <- NA

    if (is.null(heatmap_args[["rect_gp"]][["fill"]]))
        heatmap_args[["rect_gp"]][["fill"]] <- "white"

    # If row or column labels were updated, use the new values to calculate max
    # text width and height to avoid overlapping elements or unnecessary spacing
    heatmap_args[["row_names_max_width"]] <-
        max_text_width(heatmap_args[["row_labels"]],
                       gp = heatmap_args[["row_names_gp"]])

    heatmap_args[["column_names_max_height"]] <-
        max_text_width(heatmap_args[["column_labels"]],
                       gp = heatmap_args[["column_names_gp"]])

    return(heatmap_args)
}


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

# NOTE: heatmap_args, statistic_mat, padj_mat, padj_cutoff, padj_fill, and
# cell_size exist in the enrichmap function environment. See enrichmap code for
# how this is handled.
.layer_fun <- function(j, i, x, y, w, h, f) {
    # Cell background
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = heatmap_args$rect_gp$col,
                        fill = ifelse(
                          is.na(pindex(padj_mat, i, j)),
                          heatmap_args[["na_col"]],
                          ifelse(pindex(padj_mat, i, j) < padj_cutoff,
                                 padj_fill, # grey
                                 heatmap_args[["rect_gp"]][["fill"]]) # white
                        )
              ))

    # Matrix of bubble diameters (optionally scaled to row or column max)
    dmat <- -log10(padj_mat)

    if (scale_by != "max") {
        # Scale bubbles relative to row or column max
        margin <- 1L + (scale_by == "column")
        dmat <- sweep(dmat, MARGIN = margin,
                      apply(dmat, MARGIN = margin, max, na.rm = TRUE),
                      FUN = "/")
    } else {
        # Scale bubbles relative to global max
        dmat <- dmat / max(dmat, na.rm = TRUE)
    }

    # Limits on bubble diameters for significant adjusted p-values
    r_min <- 0.20
    r_max <- 0.95 # upper limit because a black border is added to the bubbles
    dmat <- ifelse(padj_mat < padj_cutoff,
                   dmat * (r_max - r_min) + r_min,
                   dmat * r_max)

    # Draw bubbles. col_fun is taken from heatmap_color_fun.
    grid.circle(
        x = x, y = y,
        r = pindex(dmat, i, j) / 2 * cell_size,
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
                          ...) {
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
