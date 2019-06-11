library("testthat")
library("Rcpp")
library("RcppArmadillo")


sourceCpp("/home/julian/Desktop/scrna-diffexpr/diffexpR/wasserstein_test.cpp")
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

test_that("permutate", {
  a <- c(1:100)
  a.len <- length(a)
  permutation1 = permutate(a)

  expect_equal(length(permutation1), a.len)
  expect_true(all(a %in% permutation1))
})

test_that("permutations", {
  a <- rnorm(100)
  n <- 1000
  df <- sapply(1:n, function(e) sample(a, replace = FALSE))
  
  permutations_set1 = permutations(a, n)
  expect_equal(length(df), length(permutations_set1))
  expect_true(all( 
                  sapply(permutations_set1, function(perm) {a %in% perm})
  ))
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

# check sanity
test_that("Bogus tests", {
  
  x <- c(1, 2, 3)
  expect_false( length(x) == 2.7 )
  expect_false( typeof(x) == "data.frame") 
  expect_true(all(x %in% sample(x)))
})

