library("testthat")
library("diffexpR")


test_that("Correctness of wasserstein single cell output", {
  # TODO
  expect_true(FALSE)
})


test_that("Example run wasserstein single cell", {
  
  # input data
  ## Interpretation: x, y are expression levels for 16 individuals in two conditions for one gene
  x <- c( 100, 120, 0, 0, 50, 150, 113, 90)
  y <- c(99, 110, 0, 0, 20, 21, 39, 100)
  ## shorter input z
  z <- c( 100, 120.00, 0, 1, 100, 150, 200)
  
  data <- matrix(c(x,y), nrow=1)
  data2 <- matrix(c(x,z), nrow=1)
  condition1 <- c(rep(0, length(x)), rep(1, length(y)))
  condition2 <- c(rep(0, length(x)), rep(2, length(z)))

  # expected output
  colnames.sc.ts <- c("d.transport","d.transport^2","d.comp^2","d.comp","location",
                   "size","shape","rho","p.nonzero","p.ad.gpd","N.exc","perc.loc",
                   "perc.size","perc.shape","decomp.error","p.zero","p.combined",
                   "p.adj.nonzero","p.adj.zero","p.adj.combined")
  colnames.sc.os <- c("d.transport","d.transport^2","d.comp^2","d.comp","location",
                   "size","shape","rho","pval","p.ad.gpd","N.exc","perc.loc",
                   "perc.size","perc.shape","decomp.error","pval.adj")
  
  # test for wasserstein.test.sc( ... method="OS")
  expect_true(all(colnames(wasserstein.sc(data, condition1, 24, 100, method="OS")) == colnames.sc.os))
  # test missing for wasserstein.test.sc( ... method="TS")
  expect_true(all(colnames(wasserstein.sc(data, condition1, 24, 100, method="TS")) == colnames.sc.ts))
  # test for different length input vectors
  expect_true(all(colnames(wasserstein.sc(data2, condition2, 24, 100, method="OS")) == colnames.sc.os))
})