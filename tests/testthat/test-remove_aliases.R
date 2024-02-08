x <- list("C" = 4:6,
          "B" = 1:3,
          "A" = 1:3,
          "D" = 4:6,
          "E" = 6:10) # not aliased

test_that("aliases are removed", {
  out <- remove_aliases(x)

  expect_equal(names(out),
               c("C", "A", "E"))
})
