#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd
#'@importFrom arm bayesglm
#'@importFrom eva gpdAd gpdFit pgpd
#'@importFrom BiocParallel bplapply
NULL

# Documentation of the RData objects
#' An empirical cumulative distribution function of a simulated Brownian bridge
#' distribution. It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#' 
#' @name empcdf.ref  
#' @docType data
#' @keywords data
empcdf.ref <- NULL

#' A global reference to empcdf.ref
#' @name global.empcdf.ref
#' @docType data
#' @keywords data
global.empcdf.ref <- NULL

# Non-exported definition to check if non-exported functions are available.
# This will cause tests for these functions (causing issues on build systems)
# to be skipped.
NONEXPORTS.AVAILABLE <- TRUE

# cleanup after our cpp libraries
.onUnload <- function (libpath) {
  library.dynam.unload("waddR", libpath)
}
