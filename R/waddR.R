#' waddR: Statistical Test for Detecting Differential Distributions Based on the
#' Wasserstein Distance
#'
#' Wasserstein distance based statistical test for detecting and describing
#' differential distributions in one-dimensional data. Functions for wasserstein
#' distance calculation, differential distribution testing, and a specialized
#' test for differential expression in scRNA data are provided.
#'
#' The Wasserstein package offers utilities for three distinct use cases:
#' \itemize{ \item Computation of the 2-Wasserstein distance \item Two-sample
#' test to check for differences between two distributions \item Detect
#' differential gene expression distributions in scRNAseq data }
#'
#' @section Wasserstein Distance functions: The 2-Wasserstein distance is a
#'   metric to describe the distance between two distributions, representing two
#'   diferent conditions A and B. This package specifically considers the
#'   squared 2-Wasserstein distance d := W^2 which offers a decomposition into
#'   location, size, and shape terms. It offers three functions to calculate the
#'   2-Wasserstein distance, all of which are implemented in Cpp and exported to
#'   R with Rcpp for better performance. \code{wasserstein_metric} is a Cpp
#'   reimplementation of the wasserstein1d method from the package
#'   \code{transport} and offers the most exact results. The functions
#'   \code{squared_wass_approx} and \code{squared_wass_decomp} compute
#'   approximations of the squared 2-Wasserstein distance with
#'   \code{suared_wass_decomp} also returning the decomosition terms for
#'   location, size, and shape. See \code{?wasserstein_metric},
#'   \code{?squared_wass_aprox}, and \code{?squared_wass_decomp} as well as the
#'   accompanying paper Schefzik and Goncalves 2019.
#'
#' @section Two-Sample Testing: This package provides two testing procedures
#'   using the 2-Wasserstein distance to test whether two distributions F_A and
#'   F_B given in the form of samples are different ba specifically testing the
#'   null hypothesis H0: F_A = F_B against the alternative hypothesis H1: 
#'   F_A != F_B.
#'
#'   The first, semi-parametric (SP), procedure uses a test based on
#'   permutations combined with a generalized pareto distribution approximation
#'   to estimate small pvalues accurately.
#'
#'   The second procedure (ASY) uses a test based on asymptotic theory which is
#'   valid only if the samples can be assumed to come from continuous
#'   distributions.
#'
#'   See the documentation of these functions \code{?wasserstein.test.sp},
#'   \test{?wasserstein.test.asy} for more details.
#'
#' @section Single Cell Test: The waddR package provides an adaptation of the
#'   semi-parametric testing procedure based on the 2-Wasserstein distance which
#'   is specifically tailored to identify differential distributions in
#'   single-cell RNA-seqencing (scRNA-seq) data. In particular, a two-stage (TS)
#'   approach has been implemented that takes account of the specific nature of
#'   scRNA-seq data by separately testing for differential proportions of zero
#'   gene expression (using a logistic regression model) and differences in
#'   non-zero gene expression (using the semi-parametric 2-Wasserstein
#'   distance-based test) between two conditions.
#'
#'   See the documentation of the Single Cell testing function
#'   \code{?wasserstein.sc} and the test for zero expression levels
#'   \code{?testZeroes} for more details.
#' 
"_PACKAGE"