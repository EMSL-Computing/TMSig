x <- list("A" = letters[1:4],
          "B" = letters[4:5],
          "C" = letters[4:5], # aliased with B
          "D" = letters[1:3], # subset of A
          "E" = c(rep("a", 10), NA)) # duplicates and NA


test_that("sets must be a named list", {

    expect_error(
        filterSets(list(c("a", "b"),
                        c("b", "c"))),
        "`x` must be a named list of character vectors."
    )

    expect_error(
        filterSets(c("a", "b", "c")),
        "`x` must be a named list of character vectors."
    )
})


test_that("an error is thrown if all sets are empty or only contain NA", {
    err1 <- expect_error(
        filterSets(list("A" = character(0L),
                        "B" = NA_character_))
    )$message

    expect_identical(
        err1,
        "All sets in `x` are empty or only contain missing values."
    )
})


test_that("min_size and max_size must each be a single integer", {

    expect_error(
        filterSets(x, min_size = c(1, 2)),
        "`min_size` must be a single integer."
    )

    expect_error(
        filterSets(x, min_size = 1, max_size = c(1, 2)),
        "`max_size` must be a single integer or Inf."
    )
})


test_that("sets must be non-NA character vectors", {

    expect_error(
        filterSets(list("A" = NA_character_, "B" = c())),
        "All sets in `x` are empty or only contain missing values."
    )

    expect_error(
        filterSets(list("A" = c(1, 2))),
        "`x` must be a named list of character vectors."
    )

})


test_that("an error is thrown when min_size is greater than max_size", {
    expect_error(
        filterSets(x, min_size = 3L, max_size = 2L),
        "`min_size` cannot be greater than `max_size`."
    )
})


test_that("an error is thrown when all sets are smaller than min_size", {
    expect_error(
        filterSets(x, min_size = 5L),
        "No sets contain at least `min_size` elements"
    )
})


test_that("an error is thrown when all sets are larger than max_size", {
    x <- list("A" = letters[1:5],
              "B" = letters[4:8])
    expect_error(
        filterSets(x, min_size = 2L, max_size = 2L),
        "All sets contain more than `max.size` elements."
    )
})


test_that("an error is thrown when no sets satisfy both size requirements", {
    x <- list("A" = "a", # satisfies max_size only
              "B" = c("a", "b", "c", "d")) # satisfies min_size only

    # sets satisfy min_size OR max_size
    expect_error(
        filterSets(x, min_size = 2L, max_size = 3L),
        "No sets satisfy both `min_size` and `max_size` thresholds."
    )
})


test_that("smallest allowable set size is 1", {
    x <- list("A" = c("a", "b"),
              "B" = c("a"),
              "C" = c("a", "b", "c"),
              "D" = c())

    expect_no_error(
        out <- filterSets(x, min_size = 1L, max_size = 1L)
    )

    expect_identical(names(out), "B")

    expect_true(all(range(lengths(out)) >= 1L))
})


test_that("duplicates and NA values are removed from the sets", {
    x <- list("A" = c("a", "a", "b", "c", NA),
              "B" = rep("a", 2))
    out <- filterSets(x, min_size = 2L)

    expect_identical(out,
                     list("A" = c("a", "b", "c")))
})


test_that("background filtering works", {
    background <- c("a", "b", "c", "e", "f", "g", NA)
    out <- filterSets(x,
                      background = background,
                      min_size = 2L)

    elements <- unique(unlist(out)) # should not contain NA

    expect_identical(names(out), c("A", "D"))

    expect_identical(elements, c("a", "b", "c"))
})


test_that("at least some elements of sets must be present in background", {
    x <- list("A" = letters[1:3],
              "B" = letters[4:8])

    expect_error(
        filterSets(x, background = letters[20:24], min_size = 2L),
        "No elements of `x` are present in `background`."
    )
})


test_that("duplicate set names will have their elements combined", {
    x <- list("A" = c("a", "b"),
              "A" = c("b", "c"))

    expect_no_error(
        out <- filterSets(x, min_size = 3L)
    )

    expect_identical(out, list("A" = c("a", "b", "c")))
})

