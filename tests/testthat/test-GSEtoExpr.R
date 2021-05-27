context("Test GSEtoExpr")

test_that("test GSEtoExpr", {
  expect_equal(GSEtoExpr("GSE128562",destdir="../../tmp"),"Conversion Done!")
})
