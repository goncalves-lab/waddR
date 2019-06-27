library("testthat")
library("diffexpR")


##########################################################################
##                    WASSERSTEIN.TEST FUNCTION                         ##
##########################################################################
test_that("Input Validation for wasserstein test", {
  
  invalid_method <- "TS"
  expect_error(wasserstein.test(rnorm(199, 20, 0.4), rnorm(239, 19, 2), 2, method=invalid_method))
  
})


test_that("Correctness of wasserstein test", {
  set.seed(42)
  result <- wasserstein.test(rnorm(100,20,3), rnorm(134, 30,10), method="SP")

  expect_true(round(result["d.transport"], 4) == round(12.01397, 4))
  expect_true(round(result["perc.loc"]) == round(62.32))
  expect_true(round(result["perc.size"]) == round(36.56))
  expect_true(round(result["perc.shape"]) == round(1.12))
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
  names.sc.os <- c("d.transport","d.transport^2","d.comp^2","d.comp","location",
                "size","shape","rho","p.nonzero","p.ad.gpd","N.exc","perc.loc",
                "perc.size","perc.shape","decomp.error","p.zero","p.combined",
                "p.adj.nonzero","p.adj.zero","p.adj.combined")
  names.sc.ts <- c("d.transport","d.transport^2","d.comp^2","d.comp","location",
                   "size","shape","rho","pval","p.ad.gpd","N.exc","perc.loc",
                   "perc.size","perc.shape","decomp.error","pval.adj")
 
  
  expect_named(wasserstein.test(v,w,method="SP"), expected=names.sp, ignore.order=TRUE)
  expect_named(wasserstein.test(v,w,seedex=34,permnum=1000,method="SP"), expected=names.sp, ignore.order=TRUE)
  expect_named(wasserstein.test(v,w,method="asy"), expected=names.asy, ignore.order=TRUE)
  
  # test missing for wasserstein.test.sc( ... method="OS")
  # test missing for wasserstein.test.sc( ... method="TS")

})