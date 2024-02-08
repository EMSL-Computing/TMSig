x <- list("A" = letters[1:4],
          "B" = letters[4:5],
          "C" = letters[4:5], # aliased with B
          "D" = letters[1:3], # subset of A
          "E" = c(rep("a", 10), NA)) # duplicates and NA


test_that("an error is thrown if all sets are empty", {
  expect_error(
    filter_sets(list("A" = c())),
    "`x` only contains empty sets."
  )
})


test_that("an error is thrown when min_size is greater than max_size", {
  expect_error(
    filter_sets(x, min_size = 3L, max_size = 2L),
    "`min_size` cannot be greater than `max_size`."
  )
})


test_that("an error is thrown when all sets are smaller than min_size", {
  expect_error(
    filter_sets(x),
    "No sets contain at least `min_size` elements"
  )
})


test_that("an error is thrown when all sets are larger than max_size", {
  x <- list("A" = letters[1:5],
            "B" = letters[4:8])
  expect_error(
    filter_sets(x, min_size = 2L, max_size = 2L),
    "All sets contain more than `max.size` elements."
  )
})


test_that("an error is thrown when no sets satisfy both size requirements", {
  x <- list("A" = "a", # satisfies max_size only
            "B" = c("a", "b", "c", "d")) # satisfies min_size only

  # sets satisfy min_size OR max_size
  expect_error(
    filter_sets(x, min_size = 2L, max_size = 3L),
    "No sets satisfy both `min_size` and `max_size` thresholds."
  )
})


test_that("smallest allowable set size is 2", {
  x <- list("A" = c("a", "b"),
            "B" = c("a"),
            "C" = c("a", "b", "c"))

  expect_no_error(
    out <- filter_sets(x, min_size = 1L, max_size = 1L)
  )

  expect_identical(names(out), "A")

  expect_true(all(range(lengths(out)) >= 2L))
})


test_that("duplicates and NA values are removed from the sets", {
  x <- list("A" = c("a", "a", "b", "c", NA),
            "B" = rep("a", 2))
  out <- filter_sets(x, min_size = 2L)

  expect_identical(out,
                   list("A" = c("a", "b", "c")))
})


test_that("background must be a character, numeric, or integer vector", {
  expect_error(
    filter_sets(x, background = list()),
    "`background` must be a character vector containing 2 or more elements."
  )
})


test_that("background must contain at least 2 elements", {
  expect_error(
    filter_sets(x, background = "a"),
    "`background` must contain at least 2 unique, nonmissing elements."
  )
})


test_that("background filtering works", {
  background <- c("a", "b", "c", "e", "f", "g", NA)
  out <- filter_sets(x,
                     background = background,
                     min_size = 2L)

  elements <- unique(unlist(out)) # should not contain NA

  expect_identical(names(out), c("A", "D"))

  expect_identical(elements, c("a", "b", "c"))
})


test_that("at least some elements of sets must be present in background", {
  x <- list("A" = letters[1:3],
            "B" = letters[4:8])

  background <- letters[20:24]

  expect_error(
    filter_sets(x, background = background, min_size = 2L),
    "No elements of `x` are present in `background`."
  )
})


test_that("duplicate set names will have their elements combined", {
  x <- list("A" = c("a", "b"),
            "A" = c("c"))

  expect_no_error(
    out <- filter_sets(x, min_size = 3L)
  )

  expect_identical(
    out,
    list("A" = c("a", "b", "c"))
  )
})

