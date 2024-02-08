test_that("x must contain at least 2 sets", {
  x <- list("A" = 1:2)

  expect_error(
    decompose_sets(x, overlap = 2L),
    "`x` must contain 2 or more sets."
  )
})


test_that("overlap must be a length 1 integer", {
  x <- list("A" = 1:2, "B" = 1:2)

  expect_error(
    decompose_sets(x, overlap = NULL),
    "`overlap` must be a single integer"
  )

  expect_error(
    decompose_sets(x, overlap = rep(2L, 2)),
    "`overlap` must be a single integer"
  )
})


test_that("at least one pair must have overlap elements in common", {
  x <- list("A" = 1:2, "B" = 2:3)

  expect_error(
    decompose_sets(x, overlap = 2L),
    paste0("No pairs of sets in `x` have at least `overlap` ",
           "elements in common.")
  )
})


test_that("sets are properly decomposed", {
  x <- list("A" = 1:3,
            "B" = 2:5, # overlaps with A
            "C" = 2:3, # subset of A, not B
            "D" = c(1, 5)) # insufficient overlap with other sets

  object <- decompose_sets(x, overlap = 2L)

  expected <- list("A NOT B" = 1,
                   "B NOT A" = 4:5,
                   "A AND B" = 2:3,
                   "A NOT C" = 1,
                   "A AND C" = 2:3,
                   "B NOT C" = 4:5,
                   "B AND C" = 2:3)
  expected <- lapply(expected, as.character)

  expect_identical(object, expected)
})

