library("testthat")
library("waddR")

##########################################################################
##                    WASSERSTEIN.TEST FUNCTION                         ##
##########################################################################
test_that("Input Validation for wasserstein test", {
  invalid_method <- "TS"
  expect_error(wasserstein.test(rnorm(199, 20, 0.4),
                                rnorm(239, 19, 2),
                                method=invalid_method))
})


test_that("Correctness of wasserstein test", {
  set.seed(42)
  x<- rnorm(100,20,3)
  y <- rnorm(134, 30,10)
  output <- wasserstein.test(x, y, method="ASY", permnum=1000)
  expect_known("value", output,
               file="known.values/testresult_correctness_wasserstein_test_1")
})


test_that("Example Run of Wasserstein Test", {
  v<-rnorm(500)
  w<-rnorm(500,1,2)
  names.sp <- c("d.wass", "d.wass^2", "d.comp^2", "d.comp", "location", "size",
                "shape", "rho", "pval", "p.ad.gpd", "N.exc", "perc.loc",
                "perc.size", "perc.shape", "decomp.error")
  names.asy <- c("d.wass","d.wass^2","d.comp^2","d.comp",
                 "location","size","shape","rho","pval","perc.loc",
                 "perc.size","perc.shape","decomp.error")
  
  expect_named(wasserstein.test(v, w, method="SP", permnum=10),
               expected=names.sp, ignore.order=TRUE)
  expect_named(wasserstein.test(v, w, method="ASY"),
               expected=names.asy, ignore.order=TRUE)
  
  # example from the wasserstein_test vignette that displayed an error
  set.seed(24)
  ctrl <- rnorm(300 ,0 ,1)
  set.seed(24)
  dd1 <- rnorm(300, 1, 1)
  
  expect_named(wasserstein.test(ctrl, dd1, method="SP", permnum=10000),
               expected=names.sp, ignore.order=TRUE)

})
