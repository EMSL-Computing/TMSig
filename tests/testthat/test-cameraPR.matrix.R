set.seed(0)
ngenes <- 6000L
ncontrasts <- 4L
statistic <- matrix(data = rnorm(n = ngenes * ncontrasts),
                    nrow = ngenes, ncol = ncontrasts,
                    dimnames = list(paste0("gene", seq_len(ngenes)),
                                    paste0("contrast", seq_len(ncontrasts))))

# Introduce missing values
statistic[1:100, 1] <- NA
statistic[200:400, 2] <- NA


test_that("each contrast must have at least 3 nonmissing values", {
    statistic2 <- statistic
    statistic2[3:ngenes, 3] <- NA # only 2 nonmissing values

    err1 <- expect_error(
        cameraPR.matrix(statistic2, index = list("A" = paste0("gene", 1:10)))
    )$message

    expect_identical(
        err1,
        "Each column of `statistic` must have at least 3 nonmissing values."
    )
})


test_that("min.size is valid", {
    index1 <- list("A" = rownames(statistic)[1:400],
                   "B" = rownames(statistic)[1:100])

    expect_no_warning(
        out <- cameraPR.matrix(statistic = statistic,
                               index = index1,
                               min.size = 101L)
    )

    # "A" is kept only
    expect_true(length(unique(out$GeneSet)) ==  1L)

    err1 <- expect_error(
        cameraPR.matrix(statistic = statistic,
                        index = index1,
                        min.size = 10000)
    )$message

    expect_identical(
        err1,
        paste0("`min.size` must be smaller than the number of non-missing ",
               "values in each contrast column of the `statistic` matrix.")
    )
})


test_that("gene sets with extreme sizes will be dropped", {
    # Warned when at least one set can not be tested
    index1 <- list("A" = "gene1",
                   "B" = c("gene1", "gene2"),
                   "C" = paste0("gene", 101:103))

    # Error when sets have insufficient data, or they include all genes in
    # statistic matrix in at least one contrast
    index2 <- c(index1[1:2], list("C" = paste0("gene", 101:ngenes)))
    err <- expect_error(
        cameraPR.matrix(statistic, index = index2)
    )$message

    expect_identical(
        err,
        paste0("No sets in `index` have at least `min.size` and fewer than ",
               "min(apply(!is.na(statistic), 2, sum)) genes with nonmissing ",
               "values in `statistic`.")
    )
})


test_that("inter.gene.cor is valid", {
    index1 <- list("A" = rownames(statistic)[1:400],
                   "B" = rownames(statistic)[1:200])

    # inter.gene.cor must be between -1 and 1 (non-inclusive)
    err0 <- expect_error(
        cameraPR.matrix(statistic, index1, inter.gene.cor = 1)
    )$message

    expect_identical(
        err0,
        "`inter.gene.cor` must be between -1 and 1."
    )

    # Length of inter.gene.cor must be 1 or the same length as index
    err1 <- expect_error(
        cameraPR.matrix(statistic, index1, inter.gene.cor = rep(0.01, 3))
    )$message

    expect_identical(
        err1,
        paste0("Length of `inter.gene.cor` must be 1 or the same length as ",
               "`index`. If the latter, names of `inter.gene.cor` should ",
               "match names of `index`.")
    )

    # Length of inter.gene.cor is correct, but names do not match names of
    # index1
    err2 <- expect_error(
        cameraPR.matrix(statistic, index1, inter.gene.cor = rep(0.01, 2))
    )$message

    expect_identical(err1, err2)

    # No errors if inter.gene.cor has the same names as index
    expect_no_error(
        cameraPR.matrix(statistic, index1,
                        inter.gene.cor = structure(rep(0.01, 2),
                                                   names = names(index1)))
    )

    # No errors even if sets are removed for being too small or too large
    index2 <- c(index1, list("C" = rownames(statistic)[1:10]))
    suppressWarnings(
        expect_no_error(
            cameraPR.matrix(statistic, index2,
                            inter.gene.cor = structure(
                                rep(0.01, length(index2)),
                                names = names(index2)
                            ))
        )
    )
})


test_that("parametric results are correct", {
    object <- cameraPR.matrix(statistic,
                              index = list("A" = paste0("gene", 1:200)))

    # These values were also checked using output from cameraPR.default to
    # establish prior validity
    expected <- data.frame(
        Contrast = factor(colnames(statistic)),
        GeneSet = "A",
        NGenes = c(100L, 199L, 200L, 200L),
        Direction = c(rep("Down", 3), "Up"),
        TwoSampleT = c(-0.38944590, -0.24192298, -1.11542975, 0.04874584),
        df = c(5898L, 5797L, 5998L, 5998L),
        PValue = NA_real_
    )

    expected$PValue[1:3] <- 2 * pt(expected$TwoSampleT[1:3],
                                   df = expected$df[1:3])
    expected$PValue[4] <- 2 * pt(expected$TwoSampleT[4],
                                 df = expected$df[4],
                                 lower.tail = FALSE)
    expected$FDR <- expected$PValue

    expect_equal(object, expected)
})


test_that("nonparametric results are correct", {
    statistic[2, ] <- statistic[1, ] # introduce ties in ranks

    object <- cameraPR.matrix(statistic,
                              index = list("A" = paste0("gene", 1:200)),
                              use.ranks = TRUE)

    # These values were also checked using output from cameraPR.default to
    # establish prior validity
    expected <- data.frame(
        Contrast = factor(colnames(statistic)),
        GeneSet = "A",
        NGenes = c(100L, 199L, 200L, 200L),
        Direction = c(rep("Down", 3), "Up"),
        PValue = c(0.38398249, 0.77902378, 0.21212502, 0.82771353)
    )
    expected$FDR <- expected$PValue

    expect_equal(object, expected)
})
