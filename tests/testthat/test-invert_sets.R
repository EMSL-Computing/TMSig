test_that("input must be a named list", {
  expect_error(
    invert_sets(list(letters[1:3])),
    "`x` must be a named list of character vectors."
  )

  expect_error(
    invert_sets(c()),
    "`x` must be a named list of character vectors."
  )
})


test_that("elements must be characters", {
  expect_error(
    invert_sets(list("A" = 1)),
    "`x` must be a named list of character vectors."
  )
})


test_that("all sets are not empty", {
  expect_error(
    invert_sets(list("A" = c(), "B" = c(NA_character_))),
    "All sets in `x` are empty or only contain missing values."
  )
})


test_that("output is correct", {
  x <- list("A" = c("a", "b", "c"),
            "B" = c("c", "d"),
            "C" = c("x", "y", "z"),
            "D" = c("a", "c", "d"))

  y <- invert_sets(x)

  expect_identical(names(y),
                   letters[c(1:4, 24:26)])

  expect_identical(y[["c"]], c("A", "B", "D"))
})
