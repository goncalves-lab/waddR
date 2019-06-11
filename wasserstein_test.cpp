// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm>
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

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
*/
// [[Rcpp::export]]
NumericVector numericVectorRep(int n, double e)
{
	NumericVector x(n);

	for(int i=0; i<x.size(); i++) {
		x[i] = e;
	}

	return x;
}


/*
Returns the average of a NumericVector object that is not Null
@param x : NumericVector object
*/
// [[Rcpp::export]]
double numericVectorMean(NumericVector x)
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
NumericVector numericVectorAbs(NumericVector x) 
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
*/
// [[Rcpp::export]]
double numericVectorSum(NumericVector x)
{
	double total = 0;

	NumericVector::iterator it;
	for(it = x.begin(); it != x.end(); ++it) {
		total += *it;
	}

	return total;

}


/*
Returns a random permutation of a given input Vector.
@param x :	Numericvector
@param replace :	boolean representing the use of inplace-permutation
@param prob :	Probability Vector, assigning to each element of input vector x
				a probability of selecting it in a pseudo-random sample
*/
// [[Rcpp::export]]
NumericVector permutate(NumericVector x, 
						bool replace=false,
						NumericVector prob = NumericVector::create())
{

	int size = x.size();
	NumericVector out = RcppArmadillo::sample(x, size, replace, prob);

	return out;
}


/*
Returns a given number of permutations of a given input NumericVector
object. 
@param x :	NumericVector representing a vector that is to be permutated
@param num_permutations : 	Int representing the number of permutations
							that are to be performed.
*/
// [[Rcpp::export]]
List permutations(NumericVector x, int num_permutations) 
{
    
	// Store permutations in List with predefined length
    List permutations_list(num_permutations);

    /* dont use names for perfomance
    // Vector for column names
    CharacterVector namevec;
    string col_name_prefix = "PERM";
    */

    bool replace = false;
    NumericVector prob = NumericVector::create();

    for (int i=0; i<num_permutations; i++){

    	permutations_list[i] = permutate(x, replace, prob);
    	/* dont use column names for performance
    	namevec.push_back(col_name_prefix + string(0,(char) i ));
    	*/
    }

    /*
    permutations_list.attr("names") = namevec;
    */
    return permutations_list;
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
			numericVectorAbs(sorted_b - sorted_a),
			p);
		double mrsad = pow(numericVectorMean(sq_abs_diff), 1/p);
		return mrsad;
	
	}

	// At least one weight vector is given
	// If only one weight vector is undefined, set all its weights to 1
	int default_weight = 1;
	NumericVector wa;
	NumericVector wb;
	if (wa_.isNull()) {

		wa = numericVectorRep(m, default_weight);
		wb = NumericVector(wb_);
	
	} else { // wb_ is null

		wa = NumericVector(wa_);
		wb = numericVectorRep(n, default_weight);
	
	}
	
	NumericVector ua(m);
	ua = (wa / numericVectorSum(wa));

}
