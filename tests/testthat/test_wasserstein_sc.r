library("testthat")
library("waddR")


test_that("Correctness of wasserstein single cell output", {
  skip("TODO")
  expect_true(FALSE)
})


test_that("Example run wasserstein single cell", {
  
  # input data: x, y represent the expression levels for 16 individuals in two
  # conditions of one gene
  x <- c( 100, 120, 0, 0, 50, 150, 113, 90)
  y <- c(99, 110, 0, 0, 20, 21, 39, 100)
  ## shorter input z
  z <- c( 100, 120.00, 0, 1, 100, 150, 200)
  
  data <- matrix(c(x,y), nrow=1)
  data2 <- matrix(c(x,z), nrow=1)
  condition1 <- c(rep(0, length(x)), rep(1, length(y)))
  condition2 <- c(rep(0, length(x)), rep(2, length(z)))

  # test for wasserstein.test.sc( ... method="OS")
  expect_known("value", wasserstein.sc(data, condition1, 100, method="OS"),
               file="known.values/testresult_example_run_wass_sc1")
  # test missing for wasserstein.test.sc( ... method="TS")
  expect_known("value", wasserstein.sc(data, condition1, 100, method="TS"),
               file="known.values/testresult_example_run_wass_sc2")
  # test for different length input vectors
  expect_known("value", wasserstein.sc(data2, condition2, 100, method="OS"),
               file="known.values/testresult_example_run_wass_sc3")
})
