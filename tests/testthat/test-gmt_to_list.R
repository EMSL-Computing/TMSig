test_that("non-GMT files throw an error", {
  expect_error(
    gmt_to_list("non_gmt_file.txt"),
    "`path` is not a path to a GMT file."
  )

  expect_error(
    gmt_to_list("non_gmt_file.gmtt"),
    "`path` is not a path to a GMT file."
  )
})


test_that("GMT files can be compressed", {
  expect_error(
    suppressWarnings(
      gmt_to_list("fake_file.gmt.gzip")
    ),
    "cannot open the connection" # fake file does not exist
  )
})


test_that("gmt_to_list output is correct", {
  path <- system.file("extdata", "c5.go.v2023.2.Hs.symbols.gmt.gz",
                      package = "ostRich")

  x <- gmt_to_list(path)

  # Correct number of sets
  expect_equal(length(x), 10461L)

  # Sets are named
  expect_true(!is.null(names(x)))

  # Spot check first set
  expect_identical(
    x[["GOBP_10_FORMYLTETRAHYDROFOLATE_METABOLIC_PROCESS"]],
    c("AASDHPPT", "ALDH1L1", "ALDH1L2", "MTHFD1", "MTHFD1L", "MTHFD2L")
  )
})
