library("Rcpp")
sourceCpp('~/Desktop/scrna-diffexpr/test/wasserstein_test.cpp')
source("~/Desktop/scrna-diffexpr/test/wasserstein_metric.r")
library("rbenchmark")

 
set.seed(13)
x <- rnorm(200)
set.seed(14)
y <- rnorm(200)

print(benchmark(wasserstein_metric(x,y,p=2), wasserstein1d(x,y,p=2)))
#z <- rnorm(1000000)
#stopifnot(all.equal(cppSort(x), sort(x)))
#print(benchmark(cppSort(z), sort(z), order="relative")[,1:4])

