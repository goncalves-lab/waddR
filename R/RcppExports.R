# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' vector_factor_addition
#'
#' @param x vector 
#' @param summand numerical
#' @return a vector containing the sum of each element in x and the summand
#'
NULL

#' vector_vector_addition
#'
#' @param x vector 
#' @param y vector
#' @return a vector containing at each position i the sum of x[i] + y[i]
#'
NULL

#' vector_factor_multiplication
#'
#' @param x vector 
#' @param factor numerical
#' @return a vector containing the product of each element in x and the factor
#'
NULL

#' vector_vector_multiplication
#'
#' @param x vector 
#' @param y vector
#' @return a vector containing at each position i the product  x[i] * y[i]
#'
NULL

#' vector_vector_subtract
#'
#' @param x vector 
#' @param y vector
#' @return vector representing the subtraction x - y
#'
NULL

#' vector_factor_division
#'
#' @param x vector 
#' @param divisor numerical
#' @return a vector containing the result of each 
#'   element in x divided by divisor
#'
NULL

#' vector_vector_division
#'
#' @param x vector 
#' @param y vector
#' @return a vector containing the product of each element
#'   in x and the factor
#'
NULL

#' vector_pow
#'
#' @param x vector 
#' @param exp numerical representing the exponent
#' @return a vector containing a copy of x with each 
#'   element raised to the power of exp
#'
NULL

#' vector_sum
#'
#' @param x vector with numerical elements 
#' @return The sum of all elements in x
#'
NULL

#' vector_mean
#'
#' @param x vector with numericals
#' @return the average of all elements in x
#'
NULL

#' vector_standard_deviation
#'
#' @param x vector with numericals
#' @return the standard deviation of all elements in x
#'
NULL

#' vector_absolute
#'
#' @param x vector with numericals
#' @return vector with the absolute values of all elements in x
#'
NULL

#' max
#'
#' @param x unsorted vector with numerals
#' @return the maximum value in x
NULL

#' min
#'
#' @param x unsorted vector with numerals
#' @return the minimum value in x
NULL

#' Repeat weighted
#'
#' Each element x[i] in the given input vector x is repeated according
#' to the weight vector at position i
#'
#' @param x vector with numeric elements
#' @param freq_table vector<int> representing the numer of repeats
#' @return the weight-repeated NumericVector of x
#'
NULL

#' vector_concatenate
#'
#' `concat` returns a vector that represents the concatenation of two input
#' vectors. The elements of the second given vector are appended to a copy of
#' the first vector.
#'
#' @param x vector
#' @param y vector
#' @return concatenation of y after x
#'
NULL

#' vector_vector_correlation
#'
#' @param x vector with numericals
#' @param y vector with numericals
#' @param mean_x optional pre-calculated mean of x
#' @param mean_y optional pre-calculated mean of y
#' @return pearsons correlation of the vectors x and y
#'
NULL

#' vector_cumulative_sum
#' 
#' The cumulative sum x[i] is the sum of all previous elements in x:
#' x[i] = x[i-1] + x[i-2] + ... + x[0]
#'
#' @param x vector of numericals
#' @return a vector containing cumulative sums of the values in x
#'
NULL

#' interval_table
#'
#' Given a vector datavec and a vector with interval breaks,
#' a histogram with the number of elements in datavec that fall into each
#' of the intervals is returned.
#'
#' @param datavec sorted vector with elements to be distributed over the
#'  intervals
#' @param interval_breaks vector with n interval_borders that are
#'  interpreted as interval breaks:\cr 
#' (-Inf, breaks[0]], (breaks[0], breaks[1]), ... , (breaks(n), Inf)
#' @param ini_value default frequency value
#'
#' @return frequency with which elements of datavec fall into each of the 
#'  intervals defined by interval_breaks
#'
NULL

#' Return permutations of a given vector as columns in a matrix
#'
#' Returns permutations of a given vector as columns in a matrix
#' @param x vector that is to be permutated
#' @param num_permutations number of permutations to be performed
#' @return a matrix, where each of the \code{num_permutations} columns represents one permutation of the input vector
#'
#' @examples
#' x <- 1:10
#' set.seed(24)
#' permutations(x, 5)
#'
#' @export
permutations <- function(x, num_permutations) {
    .Call('_waddR_permutations', PACKAGE = 'waddR', x, num_permutations)
}

#' Compute the squared 2-Wasserstein distance based on a decomposition
#'
#' Computes the squared 2-Wasserstein distance between two vectors based on a decomposition into location, size and shape terms.
#' For a detailed description of the (empirical) calculation of the invoved quantities, see Schefzik et al. (2020).
#'
#' @param x sample (vector) representing the distribution of condition \eqn{A}
#' @param y sample (vector) representing the distribution of condition \eqn{B}
#' @return A list of 4:
#' \itemize{
#' \item distance: the sum location+size+shape
#' \item location: location part in the decoposition of the 2-Wasserstein distance
#' \item size: size part in the decoposition of the 2-Wasserstein distance
#' \item shape: shape part in the decoposition of the 2-Wasserstein distance
#'}
#' 
#' @references 
#'Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
#'
#' @seealso See the functions \code{wasserstein_metric} and \code{squared_wass_approx} for
#' alternative implementations of the 2-Wasserstein distance
#'
#' @examples
#' set.seed(24)
#' x<-rnorm(100)
#' y1<-rnorm(150)
#' y2<-rexp(150,3)
#' y3<-rpois(150,2)
#'
#' squared_wass_decomp(x,y1)
#' squared_wass_decomp(x,y2)
#' squared_wass_decomp(x,y3)
#' 
#' @export
squared_wass_decomp <- function(x, y) {
    .Call('_waddR_squared_wass_decomp', PACKAGE = 'waddR', x, y)
}

#' Compute approximated squared 2-Wasserstein distance
#'
#' Calculates an approximated squared 2-Wasserstein distance based on the mean squared difference between 1000 equidistant
#' quantiles corresponding to the empirical distributions of two input vectors \eqn{x} and \eqn{y}
#'
#' @param x sample (vector) representing the distribution of condition \eqn{A}
#' @param y sample (vector) representing the distribution of condition \eqn{B}
#' @return The approximated squared 2-Wasserstein distance between \eqn{x} and \eqn{y}
#'
#' @references Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
#'
#' @seealso See the functions \code{wasserstein_metric} and \code{squared_wass_decomp} for
#' alternative implementations of the 2-Wasserstein distance
#'
#' @examples
#' set.seed(24)
#' x<-rnorm(100)
#' y1<-rnorm(150)
#' y2<-rexp(150,3)
#' y3<-rpois(150,2)
#'
#' squared_wass_approx(x,y1)
#' squared_wass_approx(x,y2)
#' squared_wass_approx(x,y3)
#'
#' @export
squared_wass_approx <- function(x, y) {
    .Call('_waddR_squared_wass_approx', PACKAGE = 'waddR', x, y)
}

#' Calculate the p-Wasserstein distance
#'
#' Calculates the \eqn{p}-Wasserstein distance (metric) between two vectors \eqn{x} and \eqn{y}
#'
#' This implementation of the \eqn{p}-Wasserstein distance is a Rcpp reimplementation of
#' the \code{wasserstein1d} function from the R package \code{transport} by Schuhmacher et al.
#' 
#' @param x sample (vector) representing the distribution of condition \eqn{A}
#' @param y sample (vector) representing the distribution of condition \eqn{B}
#' @param p order of the Wasserstein distance
#' @param wa_ optional vector of weights for \code{x}
#' @param wb_ optional vector of weights for \code{y}
#' @return The \eqn{p}-Wasserstein distance between \eqn{x} and \eqn{y}
#'
#' @references Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
#'
#' @seealso See the functions \code{squared_wass_approx} and \code{squared_wass_decomp} for
#' alternative implementations of the 2-Wasserstein distance.
#'
#' @examples
#' set.seed(24)
#' x<-rnorm(100)
#' y1<-rnorm(150)
#' y2<-rexp(150,3)
#' y3<-rpois(150,2)
#'
#' #calculate 2-Wasserstein distance between x and y1
#' wasserstein_metric(x,y1,p=2)
#' #calculate squared 2-Wasserstein distance between x and y1
#' wasserstein_metric(x,y1,p=2)^2
#'
#' #calculate 2-Wasserstein distance between x and y2
#' wasserstein_metric(x,y2,p=2)
#' #calculate squared 2-Wasserstein distance between x and y2
#' wasserstein_metric(x,y2,p=2)^2
#'
#' #calculate 2-Wasserstein distance between x and y3
#' wasserstein_metric(x,y3,p=2)
#' #calculate squared 2-Wasserstein distance between x and y3
#' wasserstein_metric(x,y3,p=2)^2
#'
#' @export
wasserstein_metric <- function(x, y, p = 1, wa_ = NULL, wb_ = NULL) {
    .Call('_waddR_wasserstein_metric', PACKAGE = 'waddR', x, y, p, wa_, wb_)
}

add_test_export <- function(x_, y_) {
    .Call('_waddR_add_test_export', PACKAGE = 'waddR', x_, y_)
}

add_test_export_sv <- function(x_, summand_) {
    .Call('_waddR_add_test_export_sv', PACKAGE = 'waddR', x_, summand_)
}

multiply_test_export <- function(x_, y_) {
    .Call('_waddR_multiply_test_export', PACKAGE = 'waddR', x_, y_)
}

multiply_test_export_sv <- function(x_, factor_) {
    .Call('_waddR_multiply_test_export_sv', PACKAGE = 'waddR', x_, factor_)
}

pow_test_export <- function(x_, exp) {
    .Call('_waddR_pow_test_export', PACKAGE = 'waddR', x_, exp)
}

abs_test_export <- function(x_) {
    .Call('_waddR_abs_test_export', PACKAGE = 'waddR', x_)
}

sum_test_export <- function(x_) {
    .Call('_waddR_sum_test_export', PACKAGE = 'waddR', x_)
}

subtract_test_export <- function(x_, y_) {
    .Call('_waddR_subtract_test_export', PACKAGE = 'waddR', x_, y_)
}

divide_test_export_sv <- function(x_, y_) {
    .Call('_waddR_divide_test_export_sv', PACKAGE = 'waddR', x_, y_)
}

divide_test_export_vectors <- function(x_, y_) {
    .Call('_waddR_divide_test_export_vectors', PACKAGE = 'waddR', x_, y_)
}

mean_test_export <- function(x_) {
    .Call('_waddR_mean_test_export', PACKAGE = 'waddR', x_)
}

sd_test_export <- function(x_) {
    .Call('_waddR_sd_test_export', PACKAGE = 'waddR', x_)
}

cumSum_test_export <- function(x_, last_index = 0L) {
    .Call('_waddR_cumSum_test_export', PACKAGE = 'waddR', x_, last_index)
}

cor_test_export <- function(x_, y_) {
    .Call('_waddR_cor_test_export', PACKAGE = 'waddR', x_, y_)
}

rep_weighted_test_export <- function(x_, weights_) {
    .Call('_waddR_rep_weighted_test_export', PACKAGE = 'waddR', x_, weights_)
}

concat_test_export <- function(x_, y_) {
    .Call('_waddR_concat_test_export', PACKAGE = 'waddR', x_, y_)
}

interval_table_test_export <- function(data_, breaks_, default_freq = 0L) {
    .Call('_waddR_interval_table_test_export', PACKAGE = 'waddR', data_, breaks_, default_freq)
}

equidist_quantile_test_export <- function(x_, K, d = 0, type = 1L) {
    .Call('_waddR_equidist_quantile_test_export', PACKAGE = 'waddR', x_, K, d, type)
}

quantile_test_export <- function(x_, q_, type = 1L) {
    .Call('_waddR_quantile_test_export', PACKAGE = 'waddR', x_, q_, type)
}

