library("testthat")
library("Rcpp")
library("RcppArmadillo")


sourceCpp("../src/wasserstein_test.cpp")

##########################################################################
##            NUMERIC VECTOR FUNCTIONS EXPOSED TO R                     ##
##########################################################################

#### NUMERIC VECTOR ABSOLUTE VALUE
test_that("NumericVectorAbs", {
  
  v1 <- c(0,2.3,2,3,4,-3,-134)
  empty <- c()
  
  expect_equal(abs(v1), abs(v1))
  expect_error(abs(empty))
})

#### NUMERIC VECTOR MEAN
test_that("numericVectorMeans", {
  
  set.seed(42)
  v <- rnorm(100)
  v2 <- c(-1,3.5,100)
  
  expect_equal(mean(v), mean(v))
  expect_equal(mean(v2), mean(v2))
})

#### PERMUTATIONS OF NUMERIC VECTOR
test_that("permutations", {
  set.seed(123)
  a <- rnorm(100)
  n <- 1000
  df <- sapply(1:n, function(e) sample(a, replace = FALSE))
  
  permutations_matrix = permutations(a, n)
  expect_equal(length(df), length(permutations_matrix))
  expect_true(all( 
                  apply(permutations_matrix, 2, function(perm) {a %in% perm})
  ))
})

#### CUMULATIVE SUM OF NUMERIC VECTOR
test_that("cumSum", {
  expect_equal(cumSum(c(1,2,3,4,5)), c(1,3,6,10,15))
  expect_equal(cumSum(c(1,2,3,4,5), 3), c(1,3,6))
  expect_equal(cumSum(c(1,2,3,4), 0), c(1,3,6,10))
})

#### INTERVAL TABLE
test_that("Interval table", {
  # setup
  set.seed(42)
  wa<- rnorm(200, 20, 1)
  set.seed(43)
  wb <- rnorm(210, 21, 1) 
  ua <- (wa/sum(wa))[-2000]
  ub <- (wb/sum(wb))[-2010]
  cua <- c(cumsum(ua))
  cub <- c(cumsum(ub))
  
  wa2 <- c(1, 1)
  wb2 <- c(2) 
  ua2 <- (wa/sum(wa))[-2]
  ub2 <- (wb/sum(wb))[-1]
  cua2 <- c(cumsum(ua))
  cub2 <- c(cumsum(ub))
  
  set.seed(123)
  wa3 <- runif(c(1:60))
  set.seed(234)
  wb3 <- runif(c(1:100))
  ua3 <- (wa/sum(wa))[-2]
  ub3 <- (wb/sum(wb))[-1]
  cua3 <- c(cumsum(ua))
  cub3 <- c(cumsum(ub))
  
  # cpp implementation
  cppresult1 <- interval_table(cub, cua,1)
  cppresult2 <- interval_table(cua2, cub2,1)
  cppresult3 <- interval_table(cua3, cub3,1)
  
  # pure R
  temp <- cut(cub, breaks = c(-Inf, cua, Inf))
  rresult1 <- table(temp) + 1
  temp <- cut(cua2, breaks = c(-Inf, cub2, Inf))
  rresult2 <- table(temp) + 1
  temp <- cut(cua3, breaks = c(-Inf, cub3, Inf))
  rresult3 <- table(temp) + 1

  expect_true(all(rresult1 == cppresult1))
  expect_true(length(rresult1) == length(cppresult1))
  expect_true(all(rresult2 == cppresult2))
  expect_true(all(rresult3 == cppresult3))
  expect_type(interval_table(rnorm(1000, 10,0.3), rnorm(834, 20, 10)), "integer")

})

#### Repeat Weighted
test_that("rep weighted", {
  expect_equal(rep_weighted(c(1,2,3),c(1,2,2)), rep(c(1,2,3), times=c(1,2,2)))
  expect_equal(rep_weighted(c(1,2,3,4), c(1,2,2,2)), c(1,2,2,3,3,4,4))
  expect_equal(rep_weighted(c(1,2), c(0,1)), c(2))
  expect_error(rep_weighted(c(1,2), c()))
})

#### Concat NumericVectors
test_that("concat", {
  
  expect_equal(concat(c(1,2,3,4), c(5,6,7,8)), c(1,2,3,4,5,6,7,8))
  expect_equal(concat(c(1,2,3), c(1)), c(1,2,3,1))
  expect_equal(concat(c(1,2,3),c(1,2,3)), c(1,2,3,1,2,3))
  expect_error(concat(c(), c(1,2,3)))
})

#### WASSERSTEIN DISTANCE 
test_that("wasserstein metric", {
  # test versus an R implementation 
  if (!requireNamespace("transport", quietly = TRUE)) {
    # Copy of the transport implementation of wasserstein metric
    wasserstein1d <- function (a, b, p = 1, wa = NULL, wb = NULL) {
      m <- length(a)
      n <- length(b)
      stopifnot(m > 0 && n > 0)
      if (m == n && is.null(wa) && is.null(wb)) {
        return( mean( abs( sort(b) - sort(a) )^p )^(1/p))
      }
      if (is.null(wa)) {
        wa <- rep(1, m)
      }
      if (is.null(wb)) {
        wb <- rep(1, n)
      }
      stopifnot(length(wa) == m && length(wb) == n)
      
      # normalizes values to sum up to 1
      ua <- (wa/sum(wa))[-m]
      ub <- (wb/sum(wb))[-n]
      
      
      cua <- c(cumsum(ua))
      cub <- c(cumsum(ub))
      temp <- cut(cub, breaks = c(-Inf, cua, Inf))
      arep <- table(temp) + 1
      temp <- cut(cua, breaks = c(-Inf, cub, Inf))
      brep <- table(temp) + 1
      
      # repeat each element in a and b as often as the intervals in arep and brep are mentionned
      aa <- rep(sort(a), times = arep)
      bb <- rep(sort(b), times = brep)
      
      # combine ecdf of weights vectors for a and b 
      uu <- sort(c(cua, cub))
      
      ## Quantiles of empirical Fa and Fb
      uu0 <- c(0, uu)
      uu1 <- c(uu, 1)
      areap <- sum((uu1 - uu0) * abs(bb - aa)^p)^(1/p)
      #    print("uu1-uu0 = ")
      #    print(uu1-uu0)
      #    print("bb-aa = ")
      #    print(bb-aa)
      #    print("(uu1 - uu0) * abs(bb - aa) = ")
      #    print((uu1 - uu0) * abs(bb - aa))
      #    print("sum((uu1 - uu0) * abs(bb - aa)^p) = ")
      #    print(sum((uu1 - uu0) * abs(bb - aa)^p))
      return(areap)
    }
  } else {
    library("transport")
  }
  
  a <- c(13,21,34,23)
  b <- c(1,1,1,2.3)
  p <- 2
  # case with equally long vectors a and b
  expect_equal(wasserstein_metric(a,b,p), wasserstein1d(a,b,p))
  expect_equal(wasserstein_metric(a,b), wasserstein1d(a,b))
  
  # vectors of different lengths
  c <- c(34,4343,3090,1309,23.2)
  set.seed(42)
  a2 <- rnorm(100,10,1)
  set.seed(24)
  b2 <- rnorm(102,10.5, 1)
  expect_equal(wasserstein_metric(a2,b2,p), wasserstein1d(a2,b2,p))
  expect_equal(wasserstein_metric(a,c,p), wasserstein1d(a,c,p))
})


##########################################################################
##                    FUNCTIONS AVAILABLE ONLY FROM CPP                 ##
##########################################################################
#### PERMUTATIONS OF ARMADILLO DATA STRUCTURES
test_that("permutations internal datatypes",{
  a <- rnorm(100)
  n <- 100
  df <- sapply(1:n, function(e) sample(a, replace = FALSE))
  
  permutations_matrix = permutations_internal_test_export(a, n)
  expect_equal(length(df), length(permutations_matrix))
  expect_true(all(apply(permutations_matrix, 2, function(perm) {a %in% perm})))
})




#### SANITY CHECK
test_that("Bogus tests", {
  
  x <- c(1, 2, 3)
  expect_false( length(x) == 2.7 )
  expect_false( typeof(x) == "data.frame") 
  expect_true(all(x %in% sample(x)))
})

