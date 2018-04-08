library(MetaGxPancreas)

context("Checking loadPancreasEsets")


test_that("ensure datasets and duplicates are properly loaded from the hub and package", {
  esetsAndDuplicates = MetaGxPancreas::loadPancreasEsets()
  esets = esetsAndDuplicates$esets
  expect_equal(class(esets[[1]])[1], "ExpressionSet")
  expect_equal(class(esets[[1]]@assayData$exprs), "matrix")
})
