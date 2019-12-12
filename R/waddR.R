#' waddR: Statistical tests for detecting differential distributions based on
#' the 2-Wasserstein distance
#'
#' The package offers statistical tests based on the 2-Wasserstein distance for detecting and characterizing differences between two distributions given in the form of samples. Functions for calculating the 2-Wasserstein distance and testing for differential distributions are provided, as welll as specifically tailored test for differential expression in single-cell RNA sequencing data. 
#'
#' The waddR package provides tools to address the following tasks:
#' \enumerate{ \item Computation of the 2-Wasserstein distance \item Two-sample
#' test to check for differences between two distributions \item Detect
#' differential gene expression distributions in single-cell RNA sequencing data }
#'
#' @section 1. 2-Wasserstein Distance functions: The 2-Wasserstein distance is a
#'   metric to quantify the difference between two distributions, representing e.g. 
#'   two different conditions A and B. The waddR package specifically considers the
#'   squared 2-Wasserstein distance which can be decomposed into
#'   location, size, and shape terms, thus providing a characterization of potential differences. It offers three functions to calculate
#'   the (squared) 2-Wasserstein distance, which are implemented in Cpp and
#'   exported to R with Rcpp for faster computation. \code{wasserstein_metric}
#'   is a Cpp reimplementation of the \code{wasserstein1d} function from the R package
#'   transport. The functions
#'   \code{squared_wass_approx} and \code{squared_wass_decomp} compute
#'   approximations of the squared 2-Wasserstein distance, with
#'   \code{squared_wass_decomp} also returning the decomposition terms for
#'   location, size, and shape. 
#'
#' See \code{?wasserstein_metric},
#'   \code{?squared_wass_aprox}, and \code{?squared_wass_decomp} as well as the
#'   accompanying paper Schefzik et al. (2019).
#'
#' @section 2. Testing for differences between distributions: The waddR package provides two testing procedures
#'   using the 2-Wasserstein distance to test whether two distributions \eqn{F_A} and
#'   \eqn{F_B} given in the form of samples are different by testing the
#'   null hypothesis \eqn{H_0: F_A = F_B} against the alternative hypothesis \eqn{H_1: F_A
#'   \neq F_B}.
#'
#'   The first, semi-parametric (SP), procedure uses a permutation-based test combined with a generalized Pareto distribution approximation
#'   to estimate small p-values accurately.
#'
#'   The second procedure uses a test based on asymptotic theory (ASY) which is
#'   valid only if the samples can be assumed to come from continuous
#'   distributions.
#'
#'   See \code{?wasserstein.test} for more
#'   details.
#'
#' @section 3. Testing for differences between distributions in the context of single-cell RNA sequencing (scRNA-seq) data: The waddR package provides an adaptation of the
#'   semi-parametric testing procedure based on the 2-Wasserstein distance
#'   which is specifically tailored to identify differential distributions in scRNA-seq data. In particular, a two-stage
#'   (TS) approach is implemented that takes account of the specific
#'   nature of scRNA-seq data by separately testing for differential
#'   proportions of zero gene expression (using a logistic regression model)
#'   and differences in non-zero gene expression (using the semi-parametric
#'   2-Wasserstein distance-based test) between two conditions.
#'
#'   See \code{?wasserstein.sc} and \code{?testZeroes} for more details.
#'   
#'@section References: Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
"_PACKAGE"
