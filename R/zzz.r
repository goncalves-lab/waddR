#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd
#'@import arm
#'@import eva
#'@import BiocParallel
NULL

# load asymptotic reference distribution
# contains:
#   value.integral  : a distribution
#   empcdf.ref      : an empirical cumulative distribution function
#                     based on the values in value.integral
load(system.file("data/ref_distr.dat", package="waddR"))

# cleanup after our cpp libraries
.onUnload <- function (libpath) {
  library.dynam.unload("waddR", libpath)
}
