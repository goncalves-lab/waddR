library("testthat")

source("/home/julian/Desktop/scrna-diffexpr/wasserstein_test.r")

test_that("Example Run of Wasserstein Test", {
  v<-rnorm(500)
  w<-rnorm(500)
  
  expect_output(wasserstein.test(v,w,method="SP"))
  expect_output(wasserstein.test(v,w,seedex=34,permnum=1000,method="SP"))
  expect_output(wasserstein.test(v,w,method="asy"))
  
})


test_that("Example Run of Wasserstein Test. (II)", {
  
  v<-rnorm(500)
  w<-rnorm(500,1,2)
  
  expect_output(wasserstein.test(v,w,method="SP"))
  expect_output(wasserstein.test(v,w,seedex=34,permnum=1000,method="SP"))
  expect_output(wasserstein.test(v,w,method="asy"))
})