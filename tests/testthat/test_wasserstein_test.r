library("testthat")
library("diffexpR")


##########################################################################
##                    WASSERSTEIN.TEST FUNCTION                         ##
##########################################################################
test_that("Input Validation for wasserstein test", {
  
  invalid_method <- "TS"
  expect_error(wasserstein.test(rnorm(199, 20, 0.4), rnorm(239, 19, 2), 2, method=invalid_method))
  
})


test_that("Example Run of Wasserstein Test (I)", {
  v<-rnorm(500)
  w<-rnorm(500)
  
  expect_output(wasserstein.test(v,w,method="SP"))
  expect_output(wasserstein.test(v,w,seedex=34,permnum=1000,method="SP"))
  expect_output(wasserstein.test(v,w,method="asy"))
  
})


test_that("Example Run of Wasserstein Test (II)", {
  
  v<-rnorm(500)
  w<-rnorm(500,1,2)
  
  expect_output(wasserstein.test(v,w,method="SP"))
  expect_output(wasserstein.test(v,w,seedex=34,permnum=1000,method="SP"))
  expect_output(wasserstein.test(v,w,method="asy"))
})