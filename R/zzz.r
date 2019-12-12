#'@useDynLib waddR
#'@importFrom Rcpp sourceCpp
#'@importFrom parallel nextRNGStream
#'@importFrom methods is
#'@importFrom stats binomial cor ecdf p.adjust pchisq quantile sd na.exclude
#'@importFrom arm bayesglm
#'@importFrom eva gpdAd gpdFit pgpd
#'@importFrom BiocParallel bplapply bpmapply
#'@importFrom BiocFileCache BiocFileCache bfcadd bfcquery bfcdownload
#'@importFrom BiocFileCache bfcpath bfcrpath bfccount bfcneedsupdate
#'@importFrom SingleCellExperiment SingleCellExperiment counts logcounts
NULL


#' Compute value of the asymptotic CDF occuring in the asymptotic theory-based test
#'
#' Computes the values of the cumulative distribution function (CDF) of the integral over the squared standard Brownian bridge in the unit interval, where the computation is based on Monte Carlo simulations.
#' This CDF occurs as an asymptotic distribution in the asymptotic theory-based test using the 2-Wasserstein distance, see Schefzik et al. (2019) for details. 
#' It is used to determine the corresponding p-values in the function \code{.wassersteinTestAsy}.
#'
#' @param v a number
#' @return Value at \code{v} of the asymptotic CDF
#' @name .brownianBridgeEmpcdf
#' 
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
.brownianBridgeEmpcdf <- function(v) NULL # distribution loaded in .onLoad


# TODO: upload as data package
#' URL for downloading the asymptotic CDF occuring in the asymptotic theory-based test
#'
#' URL for downloading the cumulative distribution function (CDF) of the integral over the squared standard Brownian bridge in the unit interval, where the computation is based on Monte Carlo simulations.
#' This CDF occurs as an asymptotic distribution in the asymptotic theory-based test using the 2-Wasserstein distance, see Schefzik et al. (2019) for details. 
#' It is used to determine the corresponding p-values in the function \code{.wassersteinTestAsy}.
#' 
#' @name brownianBridgeEmpcdf.url
#' @docType data
#' 
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
brownianBridgeEmpcdf.url <- paste0( "https://github.com/goncalves-lab/",
                                    "waddR-data/blob/master/data/",
                                    "empcdf_ref.RData?raw=true")


# Non-exported definition to check if non-exported functions are available.
# This will cause tests for these functions (causing issues on build systems)
# to be skipped.
NONEXPORTS.AVAILABLE <- TRUE


.onLoad <- function (libname, pkgname) {

    # Load the reference distributions from cache
    if (is.null(.brownianBridgeEmpcdf(0))) {
        brownianBridgeEmpcdf.path <- .cache.getOrDownload(
            url=brownianBridgeEmpcdf.url,
            rname="empcdf.ref")
        empcdf.ref <- NULL # this variable will be loaded from the file
        load(brownianBridgeEmpcdf.path)
        assign(".brownianBridgeEmpcdf", 
               empcdf.ref,
               environment(.brownianBridgeEmpcdf))
    }
}


# cleanup after our cpp libraries
.onUnload <- function (libpath) {
    library.dynam.unload("waddR", libpath)
}
