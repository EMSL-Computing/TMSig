x <- list("A" = 1:3,
          "D" = 2:3,
          "B" = 2:4,
          "C" = 8:10,
          "A" = 4)
x <- lapply(x, as.character)

incidence <- sparseIncidence(x)


test_that("set order does not change", {
    expect_no_error(
        out <- incidenceToList(incidence)
    )

    expect_identical(
        names(out),
        rownames(incidence)
    )
})


test_that("output is correct", {
    out <- incidenceToList(incidence)

    expect_identical(
        out,
        lapply(list("A" = 1:4, "D" = 2:3, "B" = 2:4, "C" = 8:10),
               as.character)
    )
})

