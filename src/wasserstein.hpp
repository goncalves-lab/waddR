
#ifndef WASSERSTEIN_H
#define WASSERSTEIN_H


#include <iostream>
#include <armadillo>
#include <math.h>
#include <RcppArmadillo.h>
//#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

template <typename T> 
std::vector<T> multiply(const std::vector<T> & x, const std::vector<T> y);

template <typename T>
std::vector<T> multiply(const std::vector<T> & x, const T & factor);

template <typename T>
std::vector<T> pow(const std::vector<T> & x, const T exp);

template <typename T>
T sum(const std::vector<T> & x);

template <typename T>
std::vector<T> subtract(const std::vector<T> & x, const std::vector<T> y);

template <typename T>
std::vector<T> divide(const std::vector<T> & x, const T divisor);

template <typename T>
double mean(const std::vector<T> & x);

template <typename T>
double sd(const std::vector<T> & x);

template <typename T>
std::vector<T> abs(const std::vector<T> & x);

template <typename T>
double cor(const std::vector<T> & x, const std::vector<T> & y, 
	double mean_x = std::numeric_limits<double>::infinity(), 
	double mean_y = std::numeric_limits<double>::infinity());

template <typename T>
std::vector<T> cumSum(const std::vector<T> & x, const int last_index=0);

template <typename T>
std::vector<T> lin_interpolated_quantiles(std::vector<T> & x, int K);

template <typename T>
std::vector<T> emp_equi_quantiles(std::vector<T> & x_, int & K,
	bool INTERPOLATE=false, bool USE_ECDF=true);

std::vector<double> rep_weighted(std::vector<double> x,
						   std::vector<int> freq_table);


std::vector<double> concat(std::vector<double> x, std::vector<double> y);

std::vector<int> interval_table(std::vector<double> datavec,
							std::vector<double> interval_breaks,
							const int init_value=0);


Rcpp::NumericMatrix permutations(Rcpp::NumericVector x, const int num_permutations);

double sq_wasserstein(Rcpp::NumericVector a_, Rcpp::NumericVector b_, double p=1);

Rcpp::List sq_wasserstein_decomp(Rcpp::NumericVector a_, Rcpp::NumericVector b_, double p=1);

double wasserstein_metric(Rcpp::NumericVector a_, 
						  Rcpp::NumericVector b_,
						  Rcpp::Nullable<Rcpp::NumericVector> wa_=R_NilValue, 
						  Rcpp::Nullable<Rcpp::NumericVector> wb_=R_NilValue,
						  double p=1);


#endif