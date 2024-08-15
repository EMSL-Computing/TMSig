test_that("overlap must be a length 1 integer", {
    x <- list("A" = c("a", "b"),
              "B" = c("a", "b"))

    expect_error(
        decomposeSets(x, overlap = NULL),
        "`overlap` must be a single integer"
    )

    expect_error(
        decomposeSets(x, overlap = rep(2L, 2)),
        "`overlap` must be a single integer"
    )
})


test_that("at least one pair must have overlap elements in common", {
    x <- list("A" = c("a", "b"),
              "B" = c("b", "c"))

    expect_error(
        decomposeSets(x, overlap = 2L),
        "No pairs of sets with at least `overlap` elements in common."
    )
})


test_that("sets are properly decomposed", {
    x <- list("A" = c("a", "b", "c"),
              "B" = c("b", "c", "d", "e"), # overlaps with A
              "C" = c("b", "c", "d"), # subset of A, not B
              "D" = c("a", "e")) # insufficient overlap with other sets

    object <- decomposeSets(x, overlap = 2L)

    expected <- list("B ~MINUS~ A" = c("d", "e"),
                     "A ~MINUS~ B" = "a",
                     "A ~AND~ B" = c("b", "c"),
                     "C ~MINUS~ A" = "d",
                     "A ~MINUS~ C" = "a",
                     "A ~AND~ C" = c("b", "c"),
                     "B ~MINUS~ C" = "e",
                     "B ~AND~ C" = c("b", "c", "d"))

    expect_identical(object, expected)
})

