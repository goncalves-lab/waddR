// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <armadillo>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace std;
using namespace Rcpp;

#define END "\n";

template <typename T> 
vector<T> multiply(const vector<T> & x, const vector<T> y);

template <typename T>
vector<T> pow(const vector<T> & x, const T exp);

template <typename T>
T sum(const vector<T> & x);

template <typename T>
vector<T> subtract(const vector<T> & x, const vector<T> y);

template <typename T>
vector<T> divide(const vector<T> & x, const T divisor);

template <typename T>
double mean(const vector<T> & x);

template <typename T>
double sd(const vector<T> & x);

template <typename T>
vector<T> abs(const vector<T> & x);

template <typename T>
double cor(const vector<T> & x, const vector<T> & y, 
	double mean_x = numeric_limits<double>::infinity(), 
	double mean_y = numeric_limits<double>::infinity());

template <typename T>
vector<T> cumSum(const vector<T> & x, const int last_index=0);

template <typename T>
vector<T> lin_interpolated_quantiles(vector<T> & x, int K);

template <typename T>
vector<T> emp_equi_quantiles(vector<T> & x_, int & K,
	bool INTERPOLATE=false, bool USE_ECDF=true);

vector<double> rep_weighted(vector<double> x,
						   vector<int> freq_table);


vector<double> concat(vector<double> x, vector<double> y);

vector<int> interval_table(vector<double> datavec,
							vector<double> interval_breaks,
							const int init_value=0);


NumericMatrix permutations(NumericVector x, const int num_permutations);

double sq_wasserstein(NumericVector a_, NumericVector b_, double p=1);

Rcpp::List sq_wasserstein_decomp(NumericVector a_, NumericVector b_, double p=1)

double wasserstein_metric(NumericVector a_, 
						  NumericVector b_, double p=1, 
						  Nullable<NumericVector> wa_=R_NilValue, 
						  Nullable<NumericVector> wb_=R_NilValue);