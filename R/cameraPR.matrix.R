#' @title Competitive Gene Set Test Accounting for Inter-gene Correlation
#'
#' @description Pre-ranked Correlation-Adjusted MEan RAnk gene set testing
#'   (CAMERA-PR) tests whether a set of genes is highly ranked relative to other
#'   genes in terms of some measure of differential expression, accounting for
#'   inter-gene correlation (Wu & Smyth, 2012). See
#'   \code{\link[limma]{cameraPR}} for details.
#'
#'   While the language is gene-centric, any \emph{a priori} groups of molecules
#'   may be tested.
#'
#' @param statistic a matrix of statistics (e.g., moderated t-statistics;
#'   possibly with \code{\link[base]{NA}}s) with genes/molecules as row names
#'   and one or more contrasts or coefficients as column names.
#' @param index a named list of sets to test. Passed to \code{\link{incidence}}.
#'   \code{index} must be a list of character vectors, not the result of
#'   \code{\link[limma]{ids2indices}}, so it is more restrictive than what
#'   \code{\link[limma]{cameraPR.default}} allows.
#' @param use.ranks logical; whether to perform a parametric test (\code{FALSE};
#'   default) or a rank-based test (\code{TRUE}).
#' @param inter.gene.cor numeric; the inter-gene correlation within tested sets.
#'   May be a single value or a named vector with names matching those of
#'   \code{index}.
#' @param sort logical; should the results of each contrast be sorted by
#'   p-value? Default is \code{TRUE}.
#' @param adjust.globally logical; whether p-values from different contrasts
#'   should be adjusted together. It is recommended to set this to \code{TRUE}
#'   when testing a set of closely related contrasts. See Section 13.3 of the
#'   LIMMA User's Guide (\code{\link[limma]{limmaUsersGuide}}) for details.
#'   Default is \code{FALSE}.
#' @inheritParams limma::cameraPR
#'
#' @returns A \code{data.frame} with the following columns:
#'
#'   \item{Contrast}{factor; the contrast of interest.}
#'   \item{GeneSet}{character; the gene set being tested.}
#'   \item{NGenes}{integer; number of genes in the set with values in the
#'   \code{statistic} matrix for a given contrast.}
#'   \item{Correlation}{numeric; inter-gene correlation (only included if
#'   \code{inter.gene.cor} was not a single value).}
#'   \item{Direction}{character; direction of change ("Up" or "Down").}
#'   \item{TwoSampleT}{numeric; two-sample t-statistic (only included if
#'   \code{use.ranks=FALSE}).}
#'   \item{df}{integer; degrees of freedom (only included if
#'   \code{use.ranks=FALSE}). Two less than the number of non-missing values in
#'   each column of the \code{statistic} matrix.}
#'   \item{PValue}{numeric; two-tailed p-value.}
#'   \item{FDR}{numeric; Benjamini and Hochberg FDR adjusted p-value (not
#'   included when only one set and contrast was tested).}
#'
#' @section Test Assumptions:
#'
#'   If \code{use.ranks=FALSE}, the parametric version of CAMERA-PR will be
#'   used. Since this is a modification of Student's two sample t-test, it is
#'   assumed that the statistics in each column of \code{statistic} are
#'   approximately normally distributed. In \code{\link[limma]{camera.default}},
#'   the moderated t-statistics are converted to z-statistics with
#'   \code{\link[limma]{zscoreT}} and used for the analysis.
#'
#'   If \code{use.ranks=TRUE}, a modified signed rank test will be used, so
#'   there are no distributional assumptions.
#'
#' @references Wu, D., and Smyth, G. K. (2012). Camera: a competitive gene set
#'   test accounting for inter-gene correlation. \emph{Nucleic Acids Research}
#'   40, e133.
#'   \href{https://doi.org/10.1093/nar/gks461}{doi:10.1093/nar/gks461}.
#'
#'   Goeman, J. J., and Bühlmann, P. (2007). Analyzing gene expression data in
#'   terms of gene sets: methodological issues. \emph{Bioinformatics} 23,
#'   980-987.
#'   \href{https://doi.org/10.1093/bioinformatics/btm051}{doi:10.1093/bioinformatics/btm051}.
#'
#' @author Di Wu, Gordon Smyth, and Tyler Sagendorf
#'
#' @seealso \code{\link[limma]{cameraPR}},
#'   \code{\link[limma]{rankSumTestWithCorrelation}}
#'
#' @import Matrix
#' @importFrom data.table data.table `:=` setcolorder setorderv rbindlist setDF
#'   frank
#' @importFrom stats p.adjust pt var
#' @importFrom limma cameraPR
#'
#' @export cameraPR.matrix
#' @exportS3Method limma::cameraPR
#'
#' @examples
#' require(limma)
#' require(stats)
#'
#' # Simulate experimental data with control and treatment groups (3 samples
#' # each)
#' group <- rep(c("control", "treatment"), each = 3)
#' design <- model.matrix(~ 0 + group)
#' contrasts <- makeContrasts(contrasts = "grouptreatment - groupcontrol",
#'                            levels = colnames(design))
#'
#' ngenes <- 1000L
#' nsamples <- length(group)
#'
#' set.seed(0)
#' y <- matrix(data = rnorm(ngenes * nsamples),
#'             nrow = ngenes, ncol = nsamples,
#'             dimnames = list(paste0("gene", seq_len(ngenes)),
#'                             make.unique(group)))
#'
#' # First set of 20 genes are genuinely differentially expressed (trt1 and trt2
#' # are lower than control)
#' index1 <- 1:20
#' y[index1, 1:3] <- y[index1, 1:3] + 1
#'
#' # Second set of 20 genes are not DE
#' index2 <- 21:40
#'
#' # Generate matrix of moderated t-statistics
#' fit <- lmFit(y, design)
#' fit.contr <- contrasts.fit(fit, contrasts = contrasts)
#' fit.smooth <- eBayes(fit.contr)
#'
#' index <- list(set1 = rownames(y)[index1],
#'               set2 = rownames(y)[index2])
#'
#' # Only set1 is DE
#' cameraPR(statistic = fit.smooth$t, index = index)
#'
#' # Non-parametric version
#' cameraPR(statistic = fit.smooth$t, index = index, use.ranks = TRUE)

cameraPR.matrix <- function(statistic,
                            index,
                            use.ranks = FALSE,
                            inter.gene.cor = 0.01,
                            sort = TRUE,
                            adjust.globally = FALSE,
                            ...)
{
  dots <- names(list(...))
  if (length(dots))
    warning("Extra arguments disregarded: ", sQuote(dots))

  genes <- rownames(statistic)
  contrast_names <- colnames(statistic)
  ncontrasts <- length(contrast_names)

  if (!all(is.numeric(statistic)))
    stop("`statistic` must be a numeric matrix.")

  if (is.null(genes) | is.null(contrast_names))
    stop("`statistic` matrix must have row and column names.")

  # Number of non-NA statistics (number of genes) in each contrast
  G <- apply(!is.na(statistic), 2, sum)

  if (any(G < 3L))
    stop("Each column of `statistic` must have at least 3 nonmissing values.")

  # data.table with columns "sets" and "elements"
  dt <- .prepare_sets(index) # this also validates index
  dt <- subset(dt, subset = elements %in% genes) # restrict sets to background

  if (nrow(dt) == 0L)
    stop("No genes in `index` match rownames of `statistic` matrix.")

  # Include genes from statistic matrix that are not members of any set
  other_genes <- setdiff(genes, dt$elements)

  if (length(other_genes)) {
    dt <- rbind(dt,
                data.table(sets = rep("GENES_NOT_IN_SETS",
                                      length(other_genes)),
                           elements = other_genes,
                           stringsAsFactors = FALSE))
    # Preserve order when splitting (affects final results if sort = FALSE)
    dt[, sets := factor(sets, levels = unique(sets))]

    index <- split(x = dt$elements, f = dt$sets)
  }

  ## Sparse incidence matrix
  imat <- incidence(index)
  # Remove row containing genes not in any set (keeps those extra genes)
  imat <- imat[rownames(imat) != "GENES_NOT_IN_SETS", genes, drop = FALSE]

  all_set_names <- rownames(imat)
  nsets <- length(all_set_names)

  m <- as.matrix(imat %*% !is.na(statistic)) # number of genes in each set
  storage.mode(m) <- "integer"

  # Identify sets that are too small or too large in at least one contrast
  extreme_sets <- which(apply(m, 1, function(mi) any(mi < 2L | mi == G)))
  if (length(extreme_sets)) {
    # If all sets will be dropped, throw an error
    if (length(extreme_sets) == nsets)
      stop("No sets in `index` have at least 2 and fewer than",
           "apply(!is.na(statistic), 2, sum) genes with nonmissing ",
           "values in `statistic`.")

    # Drop sets that are too small or too large
    warning("Sets in `index` with fewer than 2 and at most ",
            "apply(!is.na(statistic), 2, sum) genes with nonmissing ",
            "values in at least one contrast will be dropped.")
    m <- m[-extreme_sets, , drop = FALSE]
    imat <- imat[-extreme_sets, , drop = FALSE]

    all_set_names <- rownames(imat)
    nsets <- length(all_set_names)
  }

  # Check inter.gene.cor
  if (anyNA(inter.gene.cor) | is.null(inter.gene.cor))
    stop("NA or NULL `inter.gene.cor` not allowed.")

  if (any(abs(inter.gene.cor) >= 1))
    stop("`inter.gene.cor` must be between -1 and 1.")

  if (length(inter.gene.cor) > 1L) {
    inter.gene.cor <- inter.gene.cor[all_set_names]

    if (anyNA(inter.gene.cor))
      stop("Length of `inter.gene.cor` must be 1 or the same length as ",
           "index. If the latter, names of `inter.gene.cor` should match ",
           "names of `index`.")

    fixed.cor <- FALSE
  } else {
    fixed.cor <- TRUE
  }

  # Number of genes not in each set
  m2 <- t(G - t(m))

  if (use.ranks) {
    ## Based on limma::rankSumTestWithCorrelation

    # Matrix of ranks of statistics by contrast. May include NA's
    rank_mat <- apply(statistic, 2, function(ci) frank(ci, na.last = "keep"))
    rank_mat[is.na(rank_mat)] <- 0L # avoid propagating NA's in the sum

    sumRanksInSet <- as.matrix(imat %*% rank_mat)
    m_prod <- m * m2 # used several times
    U <- m_prod + m * (m + 1L) / 2 - sumRanksInSet
    mu <- 0.5 * m_prod

    sigma2 <- asin(1) +
      (m2 - 1L) * (asin(0.5) + (m - 1L) * asin(inter.gene.cor / 2)) +
      (m - 1L) * asin((inter.gene.cor + 1) / 2)
    sigma2 <- sigma2 * m_prod / 2 / pi

    # Different calculation for sigma2 if the correlation is zero
    zero.cor.idx <- which(inter.gene.cor == 0)
    if (length(zero.cor.idx)) {

      sigma2.zero.cor <- t(t(m_prod) / 12 * (G + 1L))

      # If inter.gene.cor is 0 and scalar, update all sigma2
      if (length(inter.gene.cor) == 1L) {
        sigma2 <- sigma2.zero.cor
      } else {
        sigma2[zero.cor.idx, ] <- sigma2.zero.cor[zero.cor.idx, ]
      }
    }

    # Adjust sigma2 for ties in ranks
    adjustment <- apply(rank_mat, 2, function(r) {
      NTIES <- table(r[r != 0]) # remove 0's (formerly NA's)
      sum(NTIES * (NTIES + 1L) * (NTIES - 1L)) # 0, if there are no ties
    })
    adjustment <- t(adjustment / t(m * (m + 1L) * (m - 1L)))
    sigma2 <- sigma2 * (1 - adjustment)

    # Two-sample z-statistics
    zlowertail <- (U + 0.5 - mu) / sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu) / sqrt(sigma2)

    Down <- pt(zuppertail, df = Inf, lower.tail = FALSE)
    Up <- pt(zlowertail, df = Inf)
  } else {
    # Global mean and variance of statistics in each contrast
    meanStat <- apply(statistic, 2, mean, na.rm = TRUE)
    varStat <- apply(statistic, 2, var, na.rm = TRUE)

    df.camera <- G - 2L
    vif <- 1L + (m - 1L) * inter.gene.cor

    # Replace missing values with 0 to calculate meanStatInSet
    statistic[is.na(statistic)] <- 0

    # imat %*% statistic results in the sum of statistics in each set, so divide
    # by m (set size) to get the mean
    meanStatInSet <- as.matrix(imat %*% statistic / m)

    # Difference between the means of statistics in a set and not in a set
    delta <- t(G / t(m2) * (t(meanStatInSet) - meanStat))

    # The statistics are assumed to have been sampled from Normal distributions
    # with the same variance, but potentially different means.
    varStatPooled <- t(
      ((G - 1L) * varStat - (t(delta^2 * m) * t(m2) / G)) / df.camera
    )

    # Two-sample t-statistics
    two.sample.t <- delta / sqrt(varStatPooled * (vif / m + 1 / m2))

    # Convert to a matrix to have same size as set_stats
    df.camera <- matrix(rep(df.camera, times = nsets),
                        nrow = nsets,
                        ncol = ncontrasts,
                        byrow = TRUE,
                        dimnames = dimnames(two.sample.t))

    Up <- pt(two.sample.t, df = df.camera, lower.tail = FALSE)
    Down <- pt(two.sample.t, df = df.camera)
  }

  D <- which(Down < Up)

  # Direction of change
  Direction <- matrix(data = "Up",
                      nrow = nsets, ncol = ncontrasts,
                      dimnames = list(all_set_names,
                                      contrast_names))
  Direction[D] <- "Down"

  # Two-sided p-values
  TwoSided <- 2 * pmin(Down, Up)

  ## Assemble into list of data tables for each contrast
  tab <- lapply(colnames(statistic), function(contrast_i) {
    tab_i <- data.table(
      GeneSet = all_set_names,
      NGenes = m[, contrast_i],
      Direction = Direction[, contrast_i],
      PValue = TwoSided[, contrast_i],
      stringsAsFactors = FALSE
    )

    if (!use.ranks) {
      tab_i[, `:=`(df = df.camera[, contrast_i],
                   TwoSampleT = two.sample.t[, contrast_i])]
      setcolorder(tab_i, neworder = c("TwoSampleT", "df"), before = "PValue")
    }

    return(tab_i)
  })

  names(tab) <- colnames(statistic)
  tab <- rbindlist(tab, idcol = "Contrast")
  tab[, Contrast := factor(Contrast, levels = contrast_names)]

  # Include column for inter-gene correlation
  if (!fixed.cor) {
    tab[, Correlation := inter.gene.cor[GeneSet]]
    setcolorder(tab, neworder = "Correlation", before = "Direction")
  }

  if (nrow(tab) > 1L) { # more than 1 gene set and/or contrast
    if (adjust.globally) {
      # Adjust p-values across related contrasts
      tab[, FDR := p.adjust(PValue, method = "BH")]
    } else if (nsets > 1L) {
      # Adjust p-values separately by contrast
      tab[, FDR := p.adjust(PValue, method = "BH"), by = Contrast]
    }

    # Sort by p-value
    if (sort)
      setorderv(tab, cols = c("Contrast", "PValue"), order = c(1, 1))
  }

  setDF(tab) # convert to data.frame to match cameraPR.default

  return(tab)
}
