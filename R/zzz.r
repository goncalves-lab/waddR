#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd
#'@importFrom arm bayesglm
#'@importFrom eva gpdAd gpdFit pgpd
#'@importFrom BiocParallel bplapply
#'@importFrom BiocFileCache BiocFileCache bfcadd bfcquery bfcdownload
#'@importFrom BiocFileCache bfcpath bfcrpath bfccount bfcneedsupdate
NULL

# Documentation of the RData objects
#' empcdf.ref
#'
#' An empirical cumulative distribution function of a simulated Brownian bridge
#' distribution. It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#' 
#' @name empcdf.ref  
#' @docType data
#' @keywords data
empcdf.ref <- NULL

#' brownianbridge.empcdf
#'
#' A global reference to empcdf.ref.
#' It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#'
#' @name brownianbridge.empcdf
#' @docType data
#' @keywords data
#' @export
brownianbridge.empcdf <- NULL

# TODO: FTP hosting or data package!!
#' brownianbridge.empcdf.url
#'
#' Url for downloading the simulated Brownian bridge distribution. 
#' It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#' 
#' @name brownianbridge.empcdf.url
#' @docType data
#' @export
brownianbridge.empcdf.url <- "ftp://0.0.0.0/data/empcdf_ref.RData"

# Non-exported definition to check if non-exported functions are available.
# This will cause tests for these functions (causing issues on build systems)
# to be skipped.
NONEXPORTS.AVAILABLE <- TRUE

# cleanup after our cpp libraries
.onUnload <- function (libpath) {
    library.dynam.unload("waddR", libpath)
}
