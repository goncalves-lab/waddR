library("testthat")
library("diffexpR")

##########################################################################
##                    WASSERSTEIN METRIC CALCULATION                    ##
##########################################################################

# R implementation to test against, from the package transport
if (TRUE) {#!requireNamespace("transport", quietly = TRUE)) {
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
    
    cat("a length = ",length(a))
    cat(" b length = ",length(b), "\n")
    cat("ua length = ", length(ua))
    cat(" ub length = ", length(ub), "\n")
    cat("cua length = ", length(cua))
    cat(" cub length = ", length(cub), "\n")
    cat("arep length = ", length(arep))
    cat(" brep length = ", length(brep), "\n")
    cat("arep  = ", arep)
    cat(" brep  = ", brep, "\n")
    cat("a_weighted length = ", length(aa))
    cat(" b_weighted length= ", length(bb), "\n")
    
    # combine ecdf of weights vectors for a and b 
    uu <- sort(c(cua, cub))
    
    ## Quantiles of empirical Fa and Fb
    uu0 <- c(0, uu)
    uu1 <- c(uu, 1)
    areap <- sum((uu1 - uu0) * abs(bb - aa)^p)^(1/p)
    return(areap)
  }
} else {
  library("transport")
}


##################################################
#       WASSERSTEIN IMPLEMENTATION (NEW)         #
##################################################

# test correctnes of wasserstein versus an R implementation 
test_that("wasserstein correctness", {
  
  a <- c(13, 21, 34, 23)
  b <- c(1,  1,  1,  2.3)
  p <- 2
  # case with equally long vectors a and b
  expect_equal(wasserstein(a,b,p), wasserstein1d(a,b,p), tolerance=1)
  expect_equal(wasserstein(a,b,1), wasserstein1d(a,b,1), tolerance=2)
  
  set.seed(42)
  a2 <- rnorm(100, 10, 1)
  set.seed(24)
  b2 <- rnorm(102, 10.5, 1)
  expect_equal(wasserstein(a2,b2,p), wasserstein1d(a2,b2,p), tolerance=1)
})

# test consistency of results
test_that("wasserstein consistency test", {
  #skip("skip the verbose tests")
  # check against a weird behaviour where NaN's were produced seemingly randomly
  x <- c(2, 1, 3) 
  y <- c(3, 3, 2, 6)
  results <- 1:100000
  for (i in results) { results[i] =  wasserstein(x, y, p=2)}
  first = results[1]
  expect_true(all(results == first))
} )


##################################################
# WASSERSTEIN_METRIC IMPLEMENTATION (DEPRECATED) #
##################################################

# test correctnes of wasserstein_metric versus an R implementation 
test_that("wasserstein_metric correctness", {
  #skip("wasserstein_metric will be deprecated, skipping all test for it")
  a <- c(13, 21, 34, 23)
  b <- c(1,  1,  1,  2.3)
  p <- 2
  # case with equally long vectors a and b
  expect_equal(wasserstein_metric(a,b,p), wasserstein1d(a,b,p))
  expect_equal(wasserstein_metric(a,b), wasserstein1d(a,b))
  
  # vectors of different lengths
  x <- c(34, 4343, 3090, 1309, 23.2)
  set.seed(42)
  a2 <- rnorm(100,10,1)
  set.seed(24)
  b2 <- rnorm(102,10.5, 1)
  expect_equal(wasserstein_metric(a2,b2,p), wasserstein1d(a2,b2,p))
  expect_equal(wasserstein_metric(a,x,p), wasserstein1d(a,x,p))
})

# test consistency of results
test_that("wasserstein_metric consistency test", {
  #skip("wasserstein_metric will be deprecated, skipping all test for it")
  # check against a weird behaviour where NaN's were produced seemingly randomly
  x <- c(2, 1, 3) 
  y <- c(3, 3, 2, 6)
  results <- 1:100000
  for (i in results) { results[i] =  wasserstein_metric(x, y, p=2)}
  first = results[1]
  expect_true(all(results == first))
} )

