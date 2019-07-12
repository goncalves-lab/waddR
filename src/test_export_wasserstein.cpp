#include "./wasserstein.hpp"


using namespace arma;
using namespace std;
using namespace Rcpp;

// [Rcpp::export]
NumericVector multiply_test_export(NumericVector x_, NumericVector y_)
{

	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());
	result = multiply(x, y);
	NumericVector output(result.begin(), result.end());
	return output;
}


// [Rcpp::export]
NumericVector multiply_test_export(NumericVector x_, double factor)
{

	vector<double> 	x(x_.begin(), x_.end()),
					result(x.size());
	result = multiply(x, factor);
	NumericVector output(result.begin(), result.end());
	return output;
}


// [Rcpp::export]
double mean_test_export(NumericVector x_)
{

	vector<double> 	x(x_.begin(), x_.end());
	double result;

	result = mean(x);

	return result;
	
}