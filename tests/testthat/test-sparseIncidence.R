test_that("x must be a named list of character vectors", {
    err <- "`x` must be a named list of character vectors."

    expect_error(
        sparseIncidence(list(c("a", "b"))), # no names
        err
    )

    expect_error(
        sparseIncidence(1), # not a list
        err
    )

    expect_error(
        sparseIncidence(list("A" = 1)), # elements are not strings
        err
    )
})


test_that("sparseIncidence works for a single set", {
    x <- list("A" = letters[1:5])

    expect_no_error(
        out <- sparseIncidence(x)
    )

    expect_identical(
        out,
        sparseMatrix(i = rep(1, 5), j = 1:5, x = 1,
                     dimnames = list(c("A"), letters[1:5]))
    )
})


test_that("duplicates are handled properly", {
    x <- list("A" = c("a", "b"),
              "B" = letters[rep(3, 10)])

    out <- sparseIncidence(x)

    expect_true(all(attr(out, "x") == 1))
})


test_that("missing values are discarded", {
    x <- list("A" = c("a", "b"),
              "B" = c("b", NA))

    out <- sparseIncidence(x)

    expect_identical(colnames(out), c("a", "b"))
})


test_that("all values cannot be missing", {
    x <- list("A" = rep(NA_character_, 2))

    expect_error(
        sparseIncidence(x),
        "All sets in `x` are empty or only contain missing values."
    )
})


test_that("duplicate list names combine those unique elements", {
    x <- list("A" = c("a", "b"),
              "A" = c("a", "c"),
              "B" = c("b", "c"))

    expect_no_error(
        out <- sparseIncidence(x)
    )

    expect_identical(
        out,
        sparseMatrix(i = c(1, 1, 1, 2, 2), j = c(1:3, 2:3), x = 1,
                     dimnames = list(c("A", "B"), c("a", "b", "c")))
    )
})

