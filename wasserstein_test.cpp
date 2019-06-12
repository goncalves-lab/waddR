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


/*
Returns the size of an R NumericVector object or -1 if it is null.
@param x_ : NumericVector object or R_NilValue
*/
int nullable_size(Nullable<NumericVector> x_ = R_NilValue)
{
    if (x_.isNotNull()) {
        NumericVector x(x_.get());
        return x.size();
    }
    return -1;
}


/*
Returns true if a given NumericVector object contains fewer than 1
elements or is of type R_Nilvalue.
@param x_ : NumericVector object or R_NilValue
*/
bool is_empty(Nullable<NumericVector> x_)
{
	int size = nullable_size(x_);

	if (size < 1) {
		return true;
	} else {
		return false;
	}
}


/*
Returns a NumericVector object of size n, filled with e.
@param n : size of the NumericVector that is to be created
@param e : element to initialize the ne NumericVector object with
// [[Rcpp::export]]
NumericVector rep(int n, double e)
{
	NumericVector x(n);

	for(int i=0; i<x.size(); i++) {
		x[i] = e;
	}

	return x;
}*/

/*
Returns the average of a NumericVector object that is not Null
@param x : NumericVector object
*/
// [[Rcpp::export]]
double mean(NumericVector x)
{	
	int vector_size = x.size();
	double sum = 0;
	double result = 0;

	for (int i = 0; i < vector_size; i++)
	{
		sum += x[i];
	}

	result = sum / vector_size;
	
	return result;
}


/*
Returns a NumericVector object with the absolute values from
a given input vector.
@param x : NumericVector object, not Null
*/
// [[Rcpp::export]]
NumericVector abs(NumericVector x) 
{	
	NumericVector abs_vector(x.size());
	for (int i=0; i<x.size(); i++)
	{	
		abs_vector[i] =  abs(x[i]);
	}
	
	return abs_vector;
}


/*
Returns the sum of all elements in a NumericVector object.
@param x : NumericVector object
// [[Rcpp::export]]
double sum(NumericVector x)
{
	double total = 0;

	NumericVector::iterator it;
	for(it = x.begin(); it != x.end(); ++it) {
		total += *it;
	}

	return total;

}
*/


/*
Returns the sum of all elements in a NumericVector object.
@param x : NumericVector object
// [[Rcpp::export]]
double sum(IntegerVector x)
{
	double total = 0;

	IntegerVector::iterator it;
	for(it = x.begin(); it != x.end(); ++it) {
		total += *it;
	}

	return total;

}*/


/* Returns the cumulative sums of a NumericalVector object.
@param x :	NumericalVector
*/
NumericVector cumSum(NumericVector x)
{
	NumericVector out(x.size());

	for (int i=0; i<x.size(); i++){
		if (i==0){
			out[i] = x[i];
		} else {
			out[i] = x[i] + x[i-1];
		}
	}

	return out;
}


/*
Returns permutations of a given NumericVector as columns in a NumericMatrix
object. 
@param x :	NumericVector representing a vector that is to be permutated
@param num_permutations : 	Int representing the number of permutations
							that are to be performed.
*/
// [[Rcpp::export]]
NumericMatrix permutations(NumericVector x, const int num_permutations) 
{

    // Matrix to store all permutations as columns
    int n_rows = x.size();
    int n_cols = num_permutations;
    NumericMatrix m(n_rows, n_cols);

	// create permutations
    for (int i=0; i<n_cols; i++){

    	m(_, i) = sample(x, n_rows, false);

    }

    return m;
}


/*
Return a given number of permutations of a given Vector object as 
columns in a Matrix.
@param x :	arma::vec oject
*/
mat permutations_internal(vec x, const int num_permutations)
{

	// Matrix to store all permutations as columns
	const int n_rows = x.size();
	const int n_cols = num_permutations;
	mat m(n_rows, n_cols);

	for (int i=0; i<n_cols; i++){
		m.col(i) = shuffle(x);
	}

	return m;
}


/*
Wrapper around permutations_internal for testing in R.
*/
// [[Rcpp::export]]
NumericMatrix permutations_internal_test_export(
	NumericVector x, const int num_permutations)
{
	// convert input NumericVector to arma::vec
	vec input_vec(x.size());
	for (int i=0; i<x.size(); i++) {input_vec[i] = x[i];}

	// convert outpu arma::mat to NumericMatrix
	mat result_mat(x.size(), num_permutations);
	result_mat = permutations_internal(input_vec, num_permutations);

	NumericMatrix output_matrix(x.size(), num_permutations);
	for (int i=0; i<x.size(); i++){
		for (int j=0; j<num_permutations; j++){
			output_matrix(i,j) = result_mat(i,j);
		}
	}
	return output_matrix;
}


/*
Returns the Wasserstein Metric of two input vectors a and b.
Reimplementation in Cpp of the function wasserstein1d in the package transport.


@param a : NumericVector a representing a distribution
@param b : NumericVector b representing a distribution
@param p : int p representing the exponent in the root mean squared difference
@param wa : NumericVector wa representing a weight matrix for vector a
@param wb : NumericVector wb representing a weight matrix for vector b
*/
// [[Rcpp::export]]
double wasserstein_metric(NumericVector a, 
						  NumericVector b, float p=1, 
						  Nullable<NumericVector> wa_=R_NilValue, 
						  Nullable<NumericVector> wb_=R_NilValue) {

	int m = nullable_size(a);
	int n = nullable_size(b);

	if (m < 1) {
		throw "Invalid input! a must not be empty";
	}

	if (n < 1) {
		throw "Invalid input! b must not be empty";
	}

	// No weight vectors are given
	if (m == n && wa_.isNull() && wb_.isNull()) {

		// compute root mean squared absolute difference of a and b
		// in R: mean(abs(sort(b) - sort(a))^p)^(1/p)
		NumericVector sorted_a = a.sort();
		NumericVector sorted_b = b.sort();
		NumericVector sq_abs_diff = pow(
			abs(sorted_b - sorted_a),
			p);
		double mrsad = pow(mean(sq_abs_diff), 1/p);
		return mrsad;
	
	}

	// At least one weight vector is given
	// If only one weight vector is undefined, set all its weights to 1
	int default_weight = 1;
	NumericVector wa;
	NumericVector wb;
	if (wa_.isNull()) {
		wa = rep(m, default_weight);
	} 
	if (wb_.isNull()) {
		wb = rep(n, default_weight);
	}
	
	NumericVector 	ua(m),
					ub(n),
					cua(m),
					cub(n);
	IntegerVector arep(m), brep(n);

	ua = (wa / sum(wa));//.erase(m);
	ub = (wb / sum(wb));//.erase(n);
	
	cua = cumSum(ua);
	cub = cumSum(ub);
	
	arep = table(ua);
	brep = table(ub);

	int len_aa = sum(arep);
	int len_bb = sum(brep);


	return (double) 0.1;

}
