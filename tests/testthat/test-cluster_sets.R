library(data.table)

test_that("cutoff must be between 0 and 1", {
  x <- list("A" = c("a", "b"),
            "B" = c("a", "b"))

  expect_error(
    cluster_sets(x, cutoff = -0.01),
    "`cutoff` must be between 0 and 1."
  )

  expect_error(
    cluster_sets(x, cutoff = 1.01),
    "`cutoff` must be between 0 and 1."
  )
})


test_that("x must contain 2 or more sets", {
  x <- list("A" = c("a", "b"))

  expect_error(
    cluster_sets(x),
    "`x` must contain 2 or more sets."
  )
})


test_that("a message is produced when no pairs of sets pass threshold", {
  x <- list("A" = c("a", "b"), "B" = c("c", "d"))

  expect_message(
    object <- cluster_sets(x),
    "No pair of sets passes the similarity cutoff."
  )

  expected <- data.table(set = c("A", "B"),
                         cluster = 1:2,
                         set_size = rep(2L, 2),
                         stringsAsFactors = FALSE)

  expect_identical(object, expected)
})


test_that("results are correct", {
  x <- list("A" = letters[1:5],
            "B" = letters[1:4], # subset of A
            "C" = letters[1:4], # aliased with B
            "D" = letters[1:3], # subset of A, B, C
            "E" = c("a", "a", NA), # duplicates and NA
            "F" = c("x", "y", "z"), # distinct elements
            "G" = letters[3:6]) # overlaps with A-E

  df <- cluster_sets(x, cutoff = 0.5)
  expected <- data.table(set = c("A", "B", "C", "D", "G", "E", "F"),
                         cluster = c(1L, 1L, 1L, 1L, 2L, 3L, 4L),
                         set_size = c(5L, 4L, 4L, 3L, 4L, 1L, 3L),
                         stringsAsFactors = FALSE)

  expect_identical(df, expected)
})

