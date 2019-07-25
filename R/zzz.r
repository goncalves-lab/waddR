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
NULL

#' An empirical cumulative distribution function based on the values in value.integral
#'
#' @name empcdf.ref  
#' @docType data
#' @keywords data
NULL

# Loading of the R objects value.integral and empcdf.ref
load(system.file("data/VALUE_INTEGRAL.RData", package="waddR"))
load(system.file("data/EMPCDF.RData", package="waddR"))



# cleanup after our cpp libraries
.onUnload <- function (libpath) {
  library.dynam.unload("waddR", libpath)
}
