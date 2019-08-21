#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd na.exclude
#'@importFrom arm bayesglm
#'@importFrom eva gpdAd gpdFit pgpd
#'@importFrom BiocParallel bplapply
#'@importFrom BiocFileCache BiocFileCache bfcadd bfcquery bfcdownload
#'@importFrom BiocFileCache bfcpath bfcrpath bfccount bfcneedsupdate
#'@importFrom SingleCellExperiment SingleCellExperiment counts logcounts
NULL


#' brownianbridge.empcdf
#'
#' An empirical cumulative distribution function of a simulated Brownian bridge
#' distribution. It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#'
#' @name brownianbridge.empcdf
#' 
.brownianBridgeEmpcdf <- NULL


# TODO: upload as data package
#' brownianbridge.empcdf.url
#'
#' Url for downloading the simulated Brownian bridge distribution. 
#' It is used an empirical quantile function to determine 
#' p-values in the asymptotic wasserstein test function wasserstein.test.asy
#' 
#' @name brownianbridge.empcdf.url
#' @docType data
#' 
.brownianBridgeEmpcdf.url <- paste0("https://github.com/goncalves-lab/",
                                    "waddR-data/blob/master/data/",
                                    "empcdf_ref.RData?raw=true")


# Non-exported definition to check if non-exported functions are available.
# This will cause tests for these functions (causing issues on build systems)
# to be skipped.
NONEXPORTS.AVAILABLE <- TRUE


.onLoad <- function (libname, pkgname) {

    # Load the reference distributions from cache
    if (is.null(.brownianBridgeEmpcdf)) {
        brownianBridgeEmpcdf.path <- .cache.getOrDownload(
            url=.brownianBridgeEmpcdf.url,
            rname="empcdf.ref")
        load(brownianBridgeEmpcdf.path)
        assign(".brownianBridgeEmpcdf", empcdf.ref, .GlobalEnv)
    }
}


# cleanup after our cpp libraries
.onUnload <- function (libpath) {
    library.dynam.unload("waddR", libpath)
}
