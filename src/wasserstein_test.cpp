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
Returns the average of a NumericVector object that is not Null
@param x : NumericVector object
*/
double mean(Nullable<NumericVector> x_)
{	

	if(x_.isNull()) {
		stop("Input Vector is Null");
	} else {
		NumericVector x (x_.get());
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
}


/*
Wrapper around the Rcpp function `mean` for export into R and testing.
The function `mean` is not exported due to name conflicts with R base.
*/
//' @export
// [[Rcpp::export]]
double mean_export(Nullable<NumericVector> x_)
{
	return mean(x_);
}

/*
Returns a NumericVector object with the absolute values from
a given input vector.
@param x : NumericVector object, not Null
*/
NumericVector abs(Nullable<NumericVector> x_) 
{	
	if(x_.isNull()) {
		stop("Input Vector is Null");
	} else {
		
		NumericVector x (x_.get());
		NumericVector abs_vector(x.size());
		for (int i=0; i<x.size(); i++)
		{	
			abs_vector[i] =  abs(x[i]);
		}
		
		return abs_vector;
	}

}


/*
Wrapper around the function `abs` for export into R and testing.
The function `abs` is not exported due to name conflicts with R base.
*/
//' @export
// [[Rcpp::export]]
NumericVector abs_export(Nullable<NumericVector> x_)
{
	return abs(x_);
}


/* Returns the cumulative sums of a NumericalVector object.
@param x :	NumericalVector
*/
//' @export
//[[Rcpp::export]]
NumericVector cumSum(Nullable<NumericVector> x_, const int last_index=0)
{
	if(x_.isNull()) {
		stop("Input Vector is Null");
	} else {

		NumericVector x (x_.get());
		const int upper = ((last_index > 0) && (last_index <= x.size()) 
							? last_index
							: x.size());

		NumericVector out(upper);

		for (int i=0; i<upper; i++){
			if (i==0){
				out[i] = x[i];
			} else {
				out[i] = x[i] + out[i-1];
			}
		}

		return out;
	}
}


/* Returns a weighted NumericVector of a given input vector.
@param x :	NumericVector object with the elements that are to be weighted
@param cua :	NumericVector object with cumulative distribution values of
				the weights for x
@param cub :	NumericVector object with cumulative distribution values of
				the weights for a reference vector y
*/
//' @export
// [[Rcpp::export]]
NumericVector rep_weighted(NumericVector x,
						   IntegerVector freq_table,
						   int len_x=-1)
{
	// build a new vector x_weighted, that repeats every element at position i
	// in x according to the frequency given at position i in freq_table
	int length = (len_x != -1) ? len_x : sum(freq_table);
	NumericVector x_weighted(length);

	// iterate over all fields of the new weighted vector
	NumericVector::iterator it = x_weighted.begin();
	NumericVector::iterator it_end = x_weighted.end();

	if(it != it_end) {

    	// iterator over all elements in the original vector
    	for(int i=0; i<x.size(); i++) {

    		// iterator over the number of repeats assigned 
    		// to each element in the original vector
    		for (int n=0; n<freq_table[i]; n++) {
    			
    			// copy value from original vector
    			*it = x[i];

    			// increment the iterator
    			++it;
    		}
    	}
    }

    return x_weighted;

}


/* Returns a NumericVector that represents the concatenation of two input vectors
@param x : 	NumericVector
@param y :	NumericVector
*/
//' @export
//[[Rcpp::export]]
NumericVector concat(Nullable<NumericVector> x_, Nullable<NumericVector> y_)
{	
	if (x_.isNull()){
		return y_.get();
	} else if (y_.isNull()){
		return x_.get();
	} else { // both x_ and y_ are not null

		NumericVector x(x_.get());
		NumericVector y(y_.get());

		NumericVector out(x.size() + y.size());

		NumericVector::iterator x_it = x.begin();
		NumericVector::iterator x_it_end = x.end();
		NumericVector::iterator y_it = y.begin();
		NumericVector::iterator y_it_end = y.end();
		NumericVector::iterator out_it = out.begin();
		NumericVector::iterator out_it_end = out.end();

		for(; x_it != x_it_end; ++x_it, ++out_it){
			*out_it = *x_it;
		}
		for(; y_it != y_it_end; ++y_it, ++out_it){
			*out_it = *y_it;
		}

		return out;
	}
}

/* Given a NumericVector x and a NumericVector of interval breaks,
a table with the number of elements in x that fall into each of the
intervals is returned.
@param datavec :	NumericVector with elements
@param interval_breaks : 	NumericVector with n interval_borders that are
					interpreted as interval breaks: 
					(-Inf, breaks[0]], (breaks[0], breaks[1]), ... , (breaks(n), Inf) 
*/
//' @export
// [[Rcpp::export]]
IntegerVector interval_table(NumericVector datavec,
							NumericVector interval_breaks,
							const int init_value=0)
{
	// count the elements in cua that occur
	// in the intervals defined by cub
	IntegerVector freq_table(interval_breaks.size()+1, init_value);

	double lower_bound(- numeric_limits<double>::infinity());
	double upper_bound = interval_breaks[0];
	int data_i = 0, interval_i=0;

	// iteration over the intervals (lower_bond, upper_bound]
	// defined by interval_breaks
	for (; interval_i<interval_breaks.size()+1; interval_i++){


		// iteration over the elements in cua, to count
		// how many elemtents fall into each of the intervals
		// defined by interval_breaks
		while (data_i < datavec.size()) {

			// if value at data_i in cua lies in interval
			// (lower_bound, upper_bound] 
			if ((datavec[data_i] > lower_bound) 
				&& (datavec[data_i] <= upper_bound)){ 

				// => increment the freq_table count for that interval
				// and don't check value at data_i in cua again
				++data_i;
				++freq_table[interval_i];

			} else {
				break;
			}
		}

		// update the interval:
		if (interval_i == interval_breaks.size()-1) {
			upper_bound = - numeric_limits<double>::infinity();
			lower_bound = interval_breaks[interval_breaks.size()];
		}
		else {
			upper_bound = interval_breaks[interval_i + 1];
			lower_bound = interval_breaks[interval_i];
		}
	}

	return freq_table;
}

/*
Returns permutations of a given NumericVector as columns in a NumericMatrix
object. 
@param x :	NumericVector representing a vector that is to be permutated
@param num_permutations : 	Int representing the number of permutations
							that are to be performed.
*/
//' @export
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
//' @export
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
//' @export
// [[Rcpp::export]]
double wasserstein_metric(NumericVector a, 
						  NumericVector b, float p=1, 
						  Nullable<NumericVector> wa_=R_NilValue, 
						  Nullable<NumericVector> wb_=R_NilValue) {

	const int m = nullable_size(a);
	const int n = nullable_size(b);

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
		double mrsad = std::pow((double) mean(sq_abs_diff), (double) 1/p);
		return mrsad;
	
	}

	// At least one weight vector is given
	// If only one weight vector is undefined, set all its weights to 1
	const int default_weight = 1;
	NumericVector wa(a.size(), default_weight);
	NumericVector wb(b.size(), default_weight);
	if (!wa_.isNull()) {
		wa = wa_;
	}
	if (!wb_.isNull()) {
		wb = wb_;
	}

	NumericVector 	ua(m), ub(n), cua(m-1), cub(n-1);
	IntegerVector a_rep, b_rep;

	// normalize the weights to add up to 1
	ua = (wa / sum(wa));
	ub = (wb / sum(wb));
	
	// cumulative distribution without the last value
	// => last value will be considered as an open interval to Infinity
	cua = cumSum(ua, m-1);
	cub = cumSum(ub, n-1);
	
	a_rep = interval_table(cub, cua, 1);
	b_rep = interval_table(cua, cub, 1);

	const int len_a_weighted = sum(a_rep);
	const int len_b_weighted = sum(b_rep);
	NumericVector a_weighted(len_b_weighted);
	NumericVector b_weighted(len_a_weighted);
	a_weighted = rep_weighted(a.sort(), a_rep);
	b_weighted = rep_weighted(b.sort(), b_rep);

	NumericVector uu0(cua.size() + cub.size());
	NumericVector uu1(cua.size() + cub.size());

	uu0 = concat(cua, cub);
	uu0.insert(0,0);
	uu1 = concat(cua, cub);
	uu1.push_back(1);

	double wsum = sum( (uu1 - uu0) * pow(abs(b_weighted - a_weighted), p));
	double areap = std::pow((double) wsum, (double) (1/p));

	//cout << "(uu1 - uu0) = " << NumericVector(uu1 - uu0) << END;
	//cout << "b_weighted - a_weighted = " << NumericVector(b_weighted - a_weighted) << END;
	//cout << "(uu1 - uu0) * abs(b_weighted - a_weighted) = " << NumericVector((uu1 - uu0) * abs(b_weighted - a_weighted)) << END;
	//cout << "sum( (uu1 - uu0) * pow(abs(b_weighted - a_weighted), p) = "<< sum( (uu1 - uu0) * pow(abs(b_weighted - a_weighted), p)) << END;
	return areap;
}
