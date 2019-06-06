library("testthat")
library("Rcpp")

sourceCpp("/home/julian/Desktop/scrna-diffexpr/test/wasserstein_test.cpp")

test_that("NumericVectorAbs", {
  
  v1 <- c(0,2.3,2,3,4,-3,-134)
  empty <- c()
  
  expect_equal(numericVectorAbs(v1), abs(v1))
  expect_error(numericVectorAbs(empty))
})

test_that("numericVectorMeans", {
  
  set.seed(13)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  
  expect_equal(numericVectorMean(v), mean(v))
  expect_equal(numericVectorMean(v2), mean(v2))
})


test_that("wasserstein metric", {
  # test versus an R implementation
  library("transport")
  
  a <- c(13,21,34,23)
  b <- c(1,1,1,2.3)
  p <- 2
  # case with equally long vectors a and b
  expect_equal(wasserstein_metric(a,b,p), wasserstein1d(a,b,p))
  expect_equal(wasserstein_metric(a,b), wasserstein1d(a,b))
  
})

# test for nonsense
test_that("Bogus tests", {
  
  x <- c(1, 2, 3)
  expect_false( length(x) == 2.7 )
  expect_false( typeof(x) == "data.frame") 
})

