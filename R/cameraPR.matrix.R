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
#' @param index a named list of sets to test. Passed to
#'   \code{\link{sparseIncidence}}. \code{index} must be a list of character
#'   vectors, not the result of \code{\link[limma]{ids2indices}}, so it is more
#'   restrictive than what \code{\link[limma]{cameraPR.default}} allows.
#' @param use.ranks logical; whether to perform a parametric test (\code{FALSE};
#'   default) or a rank-based test (\code{TRUE}).
#' @param inter.gene.cor numeric; the inter-gene correlation within tested sets.
#'   May be a single value or a named vector with names matching those of
#'   \code{index}.
#' @param sort logical; should the results of each contrast be sorted by
#'   p-value? Default is \code{TRUE}.
#' @param alternative character; the alternative hypothesis. Must be one of
#'   "\code{two.sided}" (default), "\code{greater}", or "\code{less}". May be
#'   abbreviated. A warning will be issued if anything other than
#'   "\code{two.sided}" is specified when \code{use.ranks=FALSE}.
#' @param adjust.globally logical; whether p-values from different contrasts
#'   should be adjusted together. It is recommended to set this to \code{TRUE}
#'   when testing a set of closely related contrasts. See Section 13.3 of the
#'   LIMMA User's Guide (\code{\link[limma]{limmaUsersGuide}}) for details.
#'   Default is \code{FALSE}.
#' @param min.size integer; the minimum set size. To be considered for testing,
#'   sets must have at least \code{min.size} elements with non-missing values in
#'   all contrasts. The default value of 2 is the minimum possible set size
#'   required for testing, though a value of 10 or higher is recommended; higher
#'   values tend to produce more robust results.
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
#'   \item{ZScore}{numeric; the z-score equivalent of \code{TwoSampleT}.}
#'   \item{PValue}{numeric; one- or two-sided (if
#'   \code{alternative="two.sided"}) p-value.}
#'   \item{FDR}{numeric; Benjamini and Hochberg FDR adjusted p-value.}
#'
#' @section Test Assumptions:
#'
#'   If \code{use.ranks=FALSE}, the parametric version of CAMERA-PR will be
#'   used. Since this is a modification of Student's two sample t-test, it is
#'   assumed that the statistics in each column of \code{statistic} are
#'   approximately Normally distributed. In \code{\link[limma]{camera.default}},
#'   the moderated t-statistics are converted to z-statistics with
#'   \code{\link[limma]{zscoreT}} and used for the analysis.
#'
#'   If \code{use.ranks=TRUE}, a modified Wilcoxon rank sum test will be used.
#'
#' @references Wu, D., and Smyth, G. K. (2012). Camera: a competitive gene set
#'   test accounting for inter-gene correlation. \emph{Nucleic Acids Research}
#'   40, e133. doi:\href{https://doi.org/10.1093/nar/gks461
#'   }{10.1093/nar/gks461}.
#'
#'   Goeman, J. J., and Bühlmann, P. (2007). Analyzing gene expression data in
#'   terms of gene sets: methodological issues. \emph{Bioinformatics} 23,
#'   980-987. doi:\href{https://doi.org/10.1093/bioinformatics/btm051
#'   }{10.1093/bioinformatics/btm051}.
#'
#' @author Di Wu, Gordon Smyth, and Tyler Sagendorf
#'
#' @seealso \code{\link[limma]{cameraPR}},
#'   \code{\link[limma]{rankSumTestWithCorrelation}}
#'
#' @importFrom data.table frank data.table `:=` rbindlist setcolorder setorderv
#'   setDF
#' @importFrom limma cameraPR zscoreT
#' @importFrom stats p.adjust pt var
#'
#' @export cameraPR.matrix
#' @exportS3Method limma::cameraPR
#'
#' @examples
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
#' # First set of 20 genes are genuinely differentially expressed
#' # (trt1 and trt2 are lower than control)
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
#' # Compute z-score equivalents of moderated t-statistics
#' statistic <- zscoreT(fit.smooth$t, fit.smooth$df.total)
#' head(statistic)
#'
#' # Only set1 is DE
#' cameraPR(statistic = statistic, index = index)
#'
#' # Non-parametric version
#' cameraPR(statistic = statistic, index = index, use.ranks = TRUE)

cameraPR.matrix <- function(statistic,
                            index,
                            use.ranks = FALSE,
                            inter.gene.cor = 0.01,
                            sort = TRUE,
                            alternative = c("two.sided", "less", "greater"),
                            adjust.globally = FALSE,
                            min.size = 2L,
                            ...) {
    dots <- names(list(...))
    if (length(dots))
        warning("Extra arguments disregarded: ", sQuote(dots))

    alternative <- match.arg(alternative, c("two.sided", "less", "greater"))

    if (any(!is.numeric(statistic)))
        stop("`statistic` must be a numeric matrix.")

    background <- rownames(statistic)
    contr.cols <- colnames(statistic)

    if (is.null(background) || is.null(contr.cols))
        stop("`statistic` matrix must have row and column names.")

    G <- colSums(!is.na(statistic))

    if (any(G < 3L))
        stop("Each column of `statistic` must have ",
             "at least 3 nonmissing values.")

    min.size <- max(2L, floor(min.size))

    if (any(min.size >= G))
        stop("`min.size` must be smaller than the number of non-missing ",
             "values in each contrast column of the `statistic` matrix.")

    index <- filterSets(x = index,
                        background = background,
                        min_size = min.size,
                        max_size = length(background) - 1L)

    imat <- sparseIncidence(index) # sparse incidence matrix
    sets <- rownames(imat)
    elements <- colnames(imat)

    # Number of genes in each set
    m <- as.matrix(imat %*% !is.na(statistic[elements, , drop = FALSE]))

    # Remove sets that are too small or too large in at least one contrast
    extreme.sets <-
        which(apply(m, 1, function(mi) any(mi < min.size | mi == G)))
    if (length(extreme.sets)) {
        # If all sets will be dropped, throw an error
        if (length(extreme.sets) == nrow(imat))
            stop("No sets in `index` have at least `min.size` and fewer than ",
                 "min(apply(!is.na(statistic), 2, sum)) genes with nonmissing ",
                 "values in `statistic`.")

        m <- m[-extreme.sets, , drop = FALSE]
        imat <- imat[-extreme.sets, , drop = FALSE]
        sets <- sets[-extreme.sets]
    }

    # Check inter.gene.cor
    if (is.null(inter.gene.cor) || anyNA(inter.gene.cor))
        stop("NA or NULL `inter.gene.cor` not allowed.")

    if (any(abs(inter.gene.cor) >= 1))
        stop("`inter.gene.cor` must be between -1 and 1.")

    fixed.cor <- length(inter.gene.cor) == 1L
    if (!fixed.cor) {
        inter.gene.cor <- inter.gene.cor[sets]

        if (anyNA(inter.gene.cor))
            stop("Length of `inter.gene.cor` must be 1 or the same length as ",
                 "`index`. If the latter, names of `inter.gene.cor` should ",
                 "match names of `index`.")
    }

    m2 <- t(G - t(m)) # number of genes not in each set

    if (use.ranks) {
        rankMat <- apply(statistic, 2, function(ci) frank(ci, na.last = "keep"))
        dimnames(rankMat) <- dimnames(statistic) # apply() removes attributes

        rankMat[is.na(rankMat)] <- 0L # avoid propagating NA's
        sumRanksInSet <-
            as.matrix(imat %*% rankMat[elements, , drop = FALSE])

        sigma2 <- asin(1) + (m - 1L) * asin((inter.gene.cor + 1) / 2) +
            (m2 - 1L) * (asin(0.5) + (m - 1L) * asin(inter.gene.cor / 2))
        sigma2 <- sigma2 / 2 / pi * m * m2

        # Different calculation for sigma2 if the correlation is zero
        zero.cor.idx <- which(inter.gene.cor == 0)
        if (length(zero.cor.idx)) {
            sigma2.zero.cor <- t(t(m / 12 * m2) * (G + 1L))

            if (fixed.cor)
                sigma2 <- sigma2.zero.cor
            else
                sigma2[zero.cor.idx, , drop = FALSE] <-
                    sigma2.zero.cor[zero.cor.idx, , drop = FALSE]
        }

        # Adjust sigma2 for ties in ranks
        adjustment <- apply(rankMat, 2, function(r) {
            NTIES <- table(r[r != 0]) # remove 0's (formerly NA's)
            sum(NTIES * (NTIES + 1L) * (NTIES - 1L)) # 0, if there are no ties
        })
        adjustment <- adjustment / G / (G + 1L) / (G - 1L)
        sigma <- sqrt(t((1 - adjustment) * t(sigma2)))

        U_minus_mu <- t((G + 1L) / 2L * t(m)) - sumRanksInSet
        zlowertail <- (U_minus_mu + 0.5) / sigma
        zuppertail <- (U_minus_mu - 0.5) / sigma

        Down <- pt(zuppertail, df = Inf, lower.tail = FALSE)
        Up <- pt(zlowertail, df = Inf)
    } else {
        # Degrees of freedom for each contrast
        df <- matrix(data = rep(G - 2L, each = nrow(m)),
                     nrow = nrow(m), ncol = ncol(m), dimnames = dimnames(m))

        meanStat <- colMeans(statistic, na.rm = TRUE)
        varStat <- apply(statistic, 2, var, na.rm = TRUE)

        # imat %*% statistic results in the sum of statistics in each set, so
        # divide by m (set size) to get the means
        statistic[is.na(statistic)] <- 0 # avoid propagating NA's
        meanStatInSet <-
            as.matrix(imat %*% statistic[elements, , drop = FALSE] / m)

        # Difference between the means of statistics in a set and not in a set
        delta <- t(G / t(m2) * (t(meanStatInSet) - meanStat))

        varStatPooled <-
            t((G - 1L) * varStat - (t(delta^2 * m) * t(m2) / G)) / df

        vif <- 1L + (m - 1L) * inter.gene.cor # variance inflation factor
        TwoSampleT <- delta / sqrt(varStatPooled * (vif / m + 1 / m2))

        Up <- pt(TwoSampleT, df = df, lower.tail = FALSE)
        Down <- pt(TwoSampleT, df = df)
    }

    Direction <- matrix(data = "Up", nrow = nrow(m), ncol = ncol(m),
                        dimnames = dimnames(m))

    switch(alternative,
           two.sided = {
               PValue <- 2 * pmin(Down, Up)
               Direction[Down < Up] <- "Down"
           },
           less = {
               PValue <- Down
               Direction[] <- "Down"
           },
           greater = {
               PValue <- Up
           })

    tab <- lapply(contr.cols, function(contrast_i) {
        tab_i <- data.table(
            GeneSet = sets,
            NGenes = m[, contrast_i],
            Direction = Direction[, contrast_i],
            stringsAsFactors = FALSE
        )

        if (!use.ranks) {
            tab_i[, `:=`(TwoSampleT = TwoSampleT[, contrast_i],
                         df = df[, contrast_i])]
            tab_i[, ZScore := zscoreT(TwoSampleT, df)]
        }
        tab_i[, PValue := PValue[, contrast_i]]

        return(tab_i)
    })

    names(tab) <- contr.cols
    tab <- rbindlist(tab, idcol = "Contrast")
    tab[, Contrast := factor(Contrast, levels = contr.cols)]

    if (!fixed.cor) {
        tab[, Correlation := inter.gene.cor[GeneSet]]
        setcolorder(tab, "Correlation", before = "Direction")
    }

    if (adjust.globally) # adjust p-values across related contrasts
        tab[, FDR := p.adjust(PValue, method = "BH")]
    else # adjust p-values separately by contrast
        tab[, FDR := p.adjust(PValue, method = "BH"), by = Contrast]

    if (sort)
        setorderv(tab, cols = c("Contrast", "PValue"), order = c(1, 1))

    setDF(tab) # convert to data.frame to match cameraPR.default

    return(tab)
}
