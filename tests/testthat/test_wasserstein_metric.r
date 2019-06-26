library("testthat")
library("diffexpR")

##########################################################################
##                    WASSERSTEIN METRIC CALCULATION                    ##
##########################################################################

# R implementation to test against, from the package transport
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

test_that("wasserstein metric", {
  # test versus an R implementation 
  
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

