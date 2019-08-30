library("testthat")
library("waddR")

##########################################################################
##                    WASSERSTEIN METRIC CALCULATION                    ##
##########################################################################

# R implementation to test against, from the package transport
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

  # repeat each element in a and b as often as the intervals in arep and brep
  # are mentionned
  aa <- rep(sort(a), times = arep)
  bb <- rep(sort(b), times = brep)
  
  # combine ecdf of weights vectors for a and b 
  uu <- sort(c(cua, cub))
  
  ## Quantiles of empirical Fa and Fb
  uu0 <- c(0, uu)
  uu1 <- c(uu, 1)
  areap <- sum((uu1 - uu0) * abs(bb - aa)^p)^(1/p)
  
  #cat("a length = ",length(a))
  #cat(" b length = ",length(b), "\n")
  #cat("ua length = ", length(ua))
  #cat(" ub length = ", length(ub), "\n")
  #cat("cua length = ", length(cua))
  #cat(" cub length = ", length(cub), "\n")
  #cat("arep length = ", length(arep))
  #cat(" brep length = ", length(brep), "\n")
  #cat("cua  = ", paste(as.character(cua), collapse=", "), "\n")
  #cat("cub  = ", paste(as.character(cua), collapse=", "), "\n")
  
  #cat("a_rep  = ", paste(as.character(arep), collapse=", "), "\n")
  #cat("b_rep  = ", paste(as.character(brep), collapse=", "), "\n")
  
  #cat("aa  = ", paste(as.character(aa), collapse=", "), "\n")
  #cat("bb  = ", paste(as.character(bb), collapse=", "), "\n")
  
  #cat("uu0  = ", paste(as.character(uu0), collapse=", "), "\n")
  #cat("uu1  = ", paste(as.character(uu1), collapse=", ")+, "\n")
  return(areap)
}


##################################################
#   WASSERSTEIN IMPLEMENTATION (Approximation)   #
##################################################

# test correctness of wasserstein approximation 
test_that("squared_wass_approx correctness", {
  
  a <- c(13, 21, 34, 23)
  b <- c(1,  1,  1,  2.3)
  # case with equally long vectors a and b
  expect_known("value", squared_wass_approx(a,b),
               file="known.values/testresult_squared_wass_approx_correctness1")
  set.seed(42)
  a2 <- rnorm(100, 10, 1)
  set.seed(24)
  b2 <- rnorm(102, 10.5, 1)
  expect_known("value", squared_wass_approx(a2,b2),
               file="known.values/testresult_squared_wass_approx_correctness3")
})


##################################################
#   WASSERSTEIN DECOMPOSITION (Approximation)    #
##################################################

# test correctnes of decomposed wasserstein approximation
test_that("squared_wass_decomp correctness", {
  a <- c(13, 21, 34, 23)
  b <- c(1,  1,  1,  2.3)
  
  # case with equally long vectors a and b
  res <- squared_wass_decomp(a,b)
  location <- (mean(a) - mean(b))**2
  expect_equal(res$location, location)
  size <- (sd(a) - sd(b))**2
  expect_equal(res$size, size)
  a.quant <- quantile(a, probs=(seq(1000)-0.5)/1000, type=1)
  b.quant <- quantile(b, probs=(seq(1000)-0.5)/1000, type=1)
  shape <- (2 * sd(a) * sd(b) * (1 - cor(a.quant, b.quant)))
  expect_equal(res$shape, shape)
  expect_equal(res$distance, location + size + shape)

  # test for special case that caused the shape calculation to fail
  x <- rep(0,107)
  y <- c(rep(0,154), 0.8463086, rep(0, 26))
  res <- squared_wass_decomp(x, y)
  location <- (mean(x) - mean(y))**2
  expect_equal(res$location, location)
  size <- (sd(x) - sd(y))**2
  expect_equal(res$size, size)
  x.quant <- quantile(x, probs=(seq(1000)-0.5)/1000, type=1)
  y.quant <- quantile(y, probs=(seq(1000)-0.5)/1000, type=1)
  # shape is zero, because sd(x) == 0
  shape <- 0
  expect_equal(res$shape, shape)
  expect_equal(res$distance, location + size + shape)
})


##################################################
#       WASSERSTEIN_METRIC IMPLEMENTATION        #
##################################################

# test correctnes of wasserstein_metric versus an R implementation 
test_that("wasserstein_metric correctness", {
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
  
  # Bug when length(a) == 1 and legth(b) > 2
  a3 <- c(23)
  b3 <- c(0,1,0,0)
  expect_equal(wasserstein_metric(a3, b3, p), wasserstein1d(a3,b3,p))
  
  # a simulated permutation procedure as used in the wasserstein.test method
  set.seed(24)
  ctrl <- rnorm(300 ,0 ,1)
  set.seed(24)
  dd1 <- rnorm(300, 1, 1)
  z <- c(ctrl,dd1)
  shuffle <- permutations(z, num_permutations=1000)
  wass_metric.values <- apply(  shuffle, 2,
                                function (k) {
                                    w <- wasserstein_metric(
                                        k[seq_len(length(x))],
                                        k[seq((length(x)+1):length(z))],
                                        p=2) **2
                                    return(w)})
  wass1d.values <- apply(  shuffle, 2,
                           function (k) {
                               w <- wasserstein1d(
                                   k[seq_len(length(x))],
                                   k[seq((length(x)+1):length(z))],
                                   p=2) **2
                               return(w)})
  expect_equal(wass_metric.values, wass1d.values)
})


# test consistency of results
test_that("wasserstein_metric consistency test", {
  x <- c(2, 1, 3) 
  y <- c(3, 3, 2, 6)
  results <- sapply(1:10000, function(i) { wasserstein_metric(x, y, p=2)})
  first = results[1]
  expect_true(all(results == first))
} )

