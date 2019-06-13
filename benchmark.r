library("Rcpp")
#source("~/Desktop/scrna-diffexpr/diffexpR/wasserstein_metric.r")
sourceCpp("/home/julian/Desktop/scrna-diffexpr/diffexpR/wasserstein_test.cpp")

#library("transport")
library("rbenchmark")

 
set.seed(13)
x <- rnorm(20000)
set.seed(14)
y <- rnorm(20000)
z <- rnorm(18000)

#print(benchmark(permutations(x, 1000),  sapply(1:1000, function(j) sample(x, length(x), replace = FALSE))))

print(benchmark(wasserstein_metric(x,z,p=2), wasserstein1d(x,z,p=2)))
#z <- rnorm(1000000)