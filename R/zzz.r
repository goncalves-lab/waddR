#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd
#'@import arm
#'@import eva
#'@import BiocParallel
NULL

# Documentation of the RData objects
#
#' This object contains a distribution
#'
#' @name value.integral
#' @docType data
#' @keywords data
value.integral <- NULL

#' An empirical cumulative distribution function based on the values in
#' value.integral
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
