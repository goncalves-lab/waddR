library("testthat")
library("waddR")

# load the reference data
# Checking for and loading the R object empcdf.ref
if (!exists("empcdf.ref")) {
  load(system.file("data/empcdf_ref.RData", package="waddR"))
}


##########################################################################
##                    WASSERSTEIN.TEST FUNCTION                         ##
##########################################################################
test_that("Input Validation for wasserstein test", {
  invalid_method <- "TS"
  expect_error(wasserstein.test(rnorm(199, 20, 0.4),
                                rnorm(239, 19, 2), 2,
                                method=invalid_method))
})


test_that("Correctness of wasserstein test", {
  set.seed(42)
  x<- rnorm(100,20,3)
  y <- rnorm(134, 30,10)
  expect_known("value", wasserstein.test(x, y, method="SP"),
               file="known.values/testresult_correctness_wasserstein_test_1")
})


test_that("Example Run of Wasserstein Test", {
  v<-rnorm(500)
  w<-rnorm(500,1,2)
  names.sp <-c("d.transport","d.transport^2","d.comp^2","d.comp","location",
            "size", "shape","rho","pval","p.ad.gpd","N.exc","perc.loc",
            "perc.size","perc.shape","decomp.error")
  names.asy <- c("d.transport","d.transport^2","d.comp^2","d.comp",
                 "location","size","shape","rho","pval","perc.loc",
                 "perc.size","perc.shape","decomp.error")
  
  expect_named(wasserstein.test(v,w,method="SP"),
               expected=names.sp, ignore.order=TRUE)
  expect_named(wasserstein.test(v,w,permnum=1000,method="SP"),
               expected=names.sp, ignore.order=TRUE)
  expect_named(wasserstein.test(v,w,method="asy"),
               expected=names.asy, ignore.order=TRUE)

})
