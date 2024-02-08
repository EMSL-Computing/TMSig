x <- list("C" = 4:6,
          "B" = 1:3, # aliased with A
          "A" = 1:3,
          "F" = 1:3, # aliased with A
          "D" = 4:6, # aliased with C
          "E" = 6:10) # not aliased

test_that("aliases are correctly identified", {
  out <- find_aliases(x)

  expect_identical(list("C" = "D",
                        "A" = c("B", "F"),
                        "E" = character(0)),
                   out)
})


test_that("a message is printed when no sets are aliased", {

  expect_message(
    out <- find_aliases(list("A" = 1:3,
                             "B" = 2:4)),
    "No sets are aliased"
  )

  expect_identical(out,
                   list("A" = character(0L),
                        "B" = character(0L)))
})
