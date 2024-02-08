test_that("similarity works with only 2 sets", {
  x <- list("A" = 1:2,
            "B" = 1:2)

  expect_identical(
    similarity(x),
    sparseMatrix(i = rep(1:2, times = 2), j = rep(1:2, each = 2), x = 1,
                 dimnames = list(names(x), names(x)))
  )
})


test_that("similarity requires 2 or more sets", {
  x <- list("A" = letters[1:5])

  expect_error(
    similarity(x),
    "`x` must contain 2 or more sets."
  )
})


test_that("similarity matrices are symmetric", {
  x <- list("A" = letters[1:5],
            "B" = letters[4:8],
            "C" = letters[4:8], # aliased with B
            "D" = letters[1:3]) # subset of A

  expect_true(isSymmetric(similarity(x, method = "jaccard")))

  expect_true(isSymmetric(similarity(x, method = "overlap")))
})


test_that("Jaccard coefficients are correct", {
  x <- list("A" = letters[1:5],
            "B" = letters[4:8],
            "C" = letters[4:8], # aliased with B
            "D" = letters[1:3]) # subset of A

  jacc <- similarity(x)

  x <- c(1, 0.25, 0.25, 0.6, 0.25, 1, 1, 0.25, 1, 1, 0.6, 1)

  expect_identical(attr(jacc, "x"), x)
})


test_that("overlap coefficients are correct", {
  x <- list("A" = letters[1:5],
            "B" = letters[4:8],
            "C" = letters[4:8], # aliased with B
            "D" = letters[1:3]) # subset of A

  overlap <- similarity(x, method = "overlap")

  x <- c(1, 0.4, 0.4, 1, 0.4, 1, 1, 0.4, 1, 1, 1, 1)

  expect_identical(attr(overlap, "x"), x)
})

