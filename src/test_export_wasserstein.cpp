// [[Rcpp::depends(RcppArmadillo)]]
#include "./wasserstein.hpp"

// [Rcpp::export]
NumericVector multiply_test_export(NumericVector x_, NumericVector y_)
{

	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());

	result = multiply(x, y);

	vector<double> output(result.begin(), result.end());

	return output;
	
}