// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <armadillo>
#include <math.h>
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

#define END "\n";

using namespace arma;
using namespace std;
using namespace Rcpp;



/*=============================================

				VECTOR OPERATIONS

==============================================*/ 

//' vector_factor_addition
//'
//' @param x vector 
//' @param summand numerical
//' @return a vector containing the sum of each element in x and the summand
//'
template <typename T>
vector<T> operator+(vector<T> x, T summand)
{
	vector<T> result(x.begin(), x.end());

	for (T & element : result) {
		element = element + summand;
	}
	return result;
}

//' vector_vector_addition
//'
//' @param x vector 
//' @param y vector
//' @return a vector containing at each position i the sum of x[i] + y[i]
//'
template <typename T>
vector<T> operator+(vector<T> x, vector<T> y)
{
	if (x.size() != y.size()) {
		stop("add: Sizes of vectors x and y are incompatible.");
	} else {
		vector<T> result(x.begin(), x.end());
		for (int i=0; i<x.size(); i++) {
			result[i] = x[i] + y[i];
		}
		return result;
	}
}

//' vector_factor_multiplication
//'
//' @param x vector 
//' @param factor numerical
//' @return a vector containing the product of each element in x and the factor
//'
template <typename T>
vector<T> operator*(const vector<T> & x, const T & factor)
{ 
  vector<T> result(x.begin(), x.end());
  
  for (T & element : result) {
    element = element * factor;
  }

  return result;
}

//' vector_vector_multiplication
//'
//' @param x vector 
//' @param y vector
//' @return a vector containing at each position i the product  x[i] * y[i]
//'
template <typename T>
vector<T> operator*(const vector<T> & x, const vector<T> & y)
{ 
  vector<T> result(x.begin(), x.end());
  if (x.size() != y.size()) {
  	warning("multiply: Sizes of vectors x and y are incompatible. Attempting multiplication by factor with y[0] ...");
  	if (y.size() >=1) {
  		const double factor = y[0];
  		result = x * factor;
  	} else {
  		stop("Invalid vector y");
  	}
  } else {
    for (int i=0; i<x.size(); i++) {
      result[i] = x[i] * y[i];
    }
  }
  return result;
}

//' vector_vector_subtract
//'
//' @param x vector 
//' @param y vector
//' @return vector representing the subtraction x - y
//'
template <typename T>
vector<T> operator-(const vector<T> & x, const vector<T> & y)
{
  if (x.size() != y.size()) {
  	stop("subtract: Sizes of vectors x and y are incompatible.");
  } else {
  	vector<T> result(x.begin(), x.end());
    for (int i=0; i<x.size(); i++) {
      result[i] = x[i] - y[i];
    }
    return result;
  }
}

//' vector_factor_division
//'
//' @param x vector 
//' @param divisor numerical
//' @return a vector containing the result of each 
//'   element in x divided by divisor
//'
template <typename T>
vector<T> operator/(const vector<T> & x, const T & divisor)
{
  vector<T> result(x.begin(), x.end());
  for (T& element : result) {
    element = element / divisor;
  }
  return result;
}

//' vector_vector_division
//'
//' @param x vector 
//' @param y vector
//' @return a vector containing the product of each element
//'   in x and the factor
//'
template <typename T>
vector<T> operator/(const vector<T> & x, const vector<T> & y)
{ 
  vector<T> result(x.begin(), x.end());
  if (x.size() != y.size()) {
  	warning("divide: Sizes of vectors x and y are incompatible. Attempting division of x by y[0] ...");
  	if (y.size() >=1) {
  		const double divisor = y[0];
  		result = x / divisor;
  	} else {
  		stop("Invalid vector y");
  	}
  } else {
    for (int i=0; i<x.size(); i++) {
      result[i] = x[i] / y[i];
    }
  }
  return result;
}


//' vector_pow
//'
//' @param x vector 
//' @param exp numerical representing the exponent
//' @return a vector containing a copy of x with each 
//'   element raised to the power of exp
//'
template <typename T>
vector<T> pow(const vector<T> & x, const T & exp)
{
  vector<T> result(x.begin(), x.end());
  for (T& element : result) {
    element = pow(element, exp);
  }
  return result;
}

//' vector_sum
//'
//' @param x vector with numerical elements 
//' @return The sum of all elements in x
//'
template <typename T>
T sum(const vector<T> & x) 
{
  T result = 0;
  for(const T& element : x) {
    result += element;
  }
  return result;
}


//' vector_mean
//'
//' @param x vector with numericals
//' @return the average of all elements in x
//'
template <typename T>
double mean(const vector<T> & x)
{	

 		double sum = 0.0;
 		double result = 0.0;

 		for (const double& element : x)
 		{
 			sum += element;
 		}

 		result = sum / x.size();
 		return result;
}

//' vector_standard_deviation
//'
//' @param x vector with numericals
//' @return the standard deviation of all elements in x
//'
template <typename T>
double sd(const vector<T> & x)
{
	double mean_x = mean(x);
	double sum = 0.0;

	for (const T& number : x) {
		sum += pow((number - mean_x), 2);
	}

	double result = sqrt(sum / (x.size() - 1));
	return result;
}

//' vector_absolute
//'
//' @param x vector with numericals
//' @return vector with the absolute values of all elements in x
//'
template <typename T>
vector<T> abs(const vector<T> & x) 
{	

	vector<T> out(x.begin(), x.end());
	for (double& element : out)
	{	
		element = abs(element);
	}
	
	return out;

}


/*=============================================

			ALGORITHMS ON VECTORS

==============================================*/ 


//' vector_vector_correlation
//'
//' @param x vector with numericals
//' @param y vector with numericals
//' @param mean_x optional pre-calculated mean of x
//' @param mean_y optional pre-calculated mean of y
//' @return pearsons correlation of the vectors x and y
//'
template <typename T>
double cor(const vector<T> & x, const vector<T> & y, 
	double mean_x = std::numeric_limits<double>::infinity(), 
	double mean_y = std::numeric_limits<double>::infinity())
{	
	if (x.size() != y.size()){
		stop("cor: Vector x and y have incompatible size.");
	}

	// calculate empirical means if needed
	if (mean_x == numeric_limits<double>::infinity()) {
		mean_x = mean(x);
	}
	if (mean_y == numeric_limits<double>::infinity()) {
		mean_y = mean(y);
	}

	double 	r = 0.0, numerator = 0.0, denom_x = 0.0, denom_y = 0.0,
			delta_xval = 0.0, delta_yval = 0.0, denominator = 0.0;

	for (int i=0; i<x.size(); i++){
		delta_xval = x[i] - mean_x;
		delta_yval = y[i] - mean_y;

		numerator += delta_xval * delta_yval;

		denom_x += pow(delta_xval, 2);
		denom_y += pow(delta_yval, 2); 
	}

	denom_y = pow(denom_y, 1.0/2.0);
	denom_x = pow(denom_x, 1.0/2.0);
	denominator = denom_x * denom_y;

	r = numerator / denominator;

	return r;

}


//' vector_cumulative_sum
//' 
//' The cumulative sum x[i] is the sum of all previous elements in x:
//' x[i] = x[i-1] + x[i-2] + ... + x[0]
//'
//' @param x vector of numericals
//' @return a vector containing cumulative sums of the values in x
//'
template <typename T>
vector<T> cumSum(const std::vector<T> & x, const int last_index=0)
{
	
	int upper = ((last_index > 0) && (last_index <= x.size()) 
				? last_index
				: x.size());

	vector<T> out(upper);

	for (int i=0; i<upper; i++){
		if (i==0){
			out[i] = x[i];
		} else {
			out[i] = x[i] + out[i-1];
		}
	}

	return out;
}


template <typename T>
vector<T> lin_interpolated_quantiles(const vector<T> & x, const int K)
{
  warning("NotYetImplemented: Using a function that currently produces dummy output!");
  
  vector<T> FOO(x.begin(), x.end());
  
  return FOO;
}


template <typename T>
vector<T> emp_equi_quantiles(const vector<T> & x_, const int & K=1000,
	const bool USE_ECDF=false, const bool INTERPOLATE=false)
{
	int n = x_.size();
	vector<T> x(n), equi_quantiles(K);

	if (USE_ECDF) {
		x = cumSum(x_);
	} else {
		x = vector<T>(x_.begin(), x_.end());
		std::sort(x.begin(), x.end());
	}

	if (INTERPOLATE){
		//NOT YET IMPLEMENTED
		equi_quantiles = lin_interpolated_quantiles(x, K);

	} else {

		double x_pos = 0.0;
		int q_idx=0, q_num = 1, x_idx;
		for (; q_idx < K; q_num++, q_idx++) {

			x_pos = n * q_num / K ;
			x_idx = (int) min(max((double) 0.0, ceil(x_pos) - 1), (double) n-1);
			
			// check if the data index x_pos is a whole number
			if (x_pos == (double) x_idx) {
				equi_quantiles[q_idx] = 1/2 * (x[x_idx] + x[x_idx + 1]);
			} else { // not a whole number 
				// => x_idx has been rounded down to be used as an index
				equi_quantiles[q_idx] = x[x_idx];
			}
		}

	}

	return equi_quantiles;
}


//' Repeat weighted
//'
//' Each element x[i] in the given input vector x is repeated according
//' to the weight vector at position i
//'
//' @param x vector with numeric elements
//' @param freq_table vector<int> representing the numer of repeats
//' @return the weight-repeated NumericVector of x
//'
vector<double> rep_weighted(vector<double> x,
						   vector<int> freq_table)
{
	// build a new vector x_weighted, that repeats every element at position i
	// in x according to the frequency given at position i in freq_table
	int length = sum(freq_table);
	vector<double> x_weighted(length);

	// iterate over all fields of the new weighted vector
	vector<double>::iterator it = x_weighted.begin();
	vector<double>::iterator it_end = x_weighted.end();

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


//' vector_concatenate
//'
//' `concat` returns a vector that represents the concatenation of two input vectors.
//' The elements of the second given vector are appended to a copy of the first vector.
//'
//' @param x vector
//' @param y vector
//' @return concatenation of y after x
//'
vector<double> concat(vector<double> & x, vector<double> & y)
{	
	vector<double> out(x.size() + y.size());
	vector<double>::iterator x_it = x.begin();
	vector<double>::iterator x_it_end = x.end();
	vector<double>::iterator y_it = y.begin();
	vector<double>::iterator y_it_end = y.end();
	vector<double>::iterator out_it = out.begin();
	vector<double>::iterator out_it_end = out.end();

	for(; x_it != x_it_end; ++x_it, ++out_it){
		*out_it = *x_it;
	}
	for(; y_it != y_it_end; ++y_it, ++out_it){
		*out_it = *y_it;
	}

	return out;

}


//' interval_table
//'
//' Given a vector datavec and a vector with interval breaks,
//' a histogram with the number of elements in datavec that fall into each
//' of the intervals is returned.
//'
//' @param datavec vector with elements to be distributed over the intervals
//' @param interval_breaks vector with n interval_borders that are
//'		interpreted as interval breaks:\cr 
//'		(-Inf, breaks[0]], (breaks[0], breaks[1]), ... , (breaks(n), Inf)
//' @param ini_value default frequency value
//'
//' @return frequency with which elements of datavec fall into each of the 
//'    intervals defined by interval_breaks
//'
vector<int> interval_table(	const vector<double> & datavec,
							const vector<double> & interval_breaks,
							const int init_value=0)
{
	// count the elements in cua that occur
	// in the intervals defined by cub
	vector<int> freq_table(interval_breaks.size()+1, init_value);

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
			upper_bound = numeric_limits<double>::infinity();
			lower_bound = interval_breaks[interval_breaks.size()];
		}
		else {
			upper_bound = interval_breaks[interval_i + 1];
			lower_bound = interval_breaks[interval_i];
		}
	}

	return freq_table;
}


//' permutations
//'
//' Returns permutations of a given NumericVector as columns in a NumericMatrix
//' object. 
//' @param x :	NumericVector representing a vector that is to be permutated
//' @param num_permutations : 	Int representing the number of permutations
//'							that are to be performed.
//' @return a matrix containing in every column one permutations of the input vector
//'
//' @export
// [[Rcpp::export]]
NumericMatrix permutations(const NumericVector x, const int num_permutations) 
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


/*=============================================

		WASSERSTEIN DISTANCE CALCULATION

==============================================*/


//' sq_wasserstein
//'
//' Squared Wasserstein distance between two vectors.
//'
//' The squared wasserstein distance has useful properties for decomposition.
//' @param a Vector representing an empirical distribution under condition A 
//' @param b Vector representing an empirical distribution under condition B
//'	@param p exponent of the wasserstine distance
//' @return The squared wasserstein distance between a and b
//'
//' @export
//[[Rcpp::export]]
double sq_wasserstein(const NumericVector & a_, const NumericVector & b_, const double & p=1)
{

	int NUM_QUANTILES = 1000;
	vector<double> 	a(a_.begin(), a_.end()),
					b(b_.begin(), b_.end()),
					quantiles_a = emp_equi_quantiles(a, NUM_QUANTILES),
					quantiles_b = emp_equi_quantiles(b, NUM_QUANTILES);

	double location, shape, size, d, mean_a = mean(a), mean_b = mean(b),
		sd_a = sd(a), sd_b = sd(b), 
		quantile_cor_ab = cor(quantiles_a, quantiles_b);

	location = pow(mean_a - mean_b, 2);
	size = pow(sd_a - sd_b, 2);
	shape = 2 * sd_a * sd_b * (1 - quantile_cor_ab);
	d = location + size + shape;

	/*
	// DEBUGGING ...
	Rcout << "mean_a = " << mean_a <<END;
	Rcout << "mean_b = " << mean_b <<END;
	Rcout << "sd_a = " << sd_a <<END;
	Rcout << "sd_b = " << sd_b<<END;
	Rcout << "quantile_cor_ab = " << quantile_cor_ab <<END;
	Rcout << "location = " << location <<END;
	Rcout << "size = " << size <<END;
	Rcout << "shape = " << shape <<END;
	Rcout << "d = " << d <<END;
	*/
	return d;
}

//' squared_wass_decomp
//'
//' Approximation of the squared Wasserstein distance between two vectors,
//' decomposed into size, location and shape. 
//' As an approximation of the distribution, 1000 quantiles are computed for each vector.
//'
//' @param a Vector representing an empirical distribution under condition A 
//' @param b Vector representing an empirical distribution under condition B
//'	@param p exponent of the wasserstine distance
//' @return An named Rcpp::List with the wasserstein distance between a and b,
//'    decomposed into terms for size, location, and shape
//'
//' @export
//[[Rcpp::export]]
Rcpp::List squared_wass_decomp(const NumericVector & a_, const NumericVector & b_,const double & p=1)
{

	int NUM_QUANTILES = 1000;
	vector<double> 		a(a_.begin(), a_.end()), b(b_.begin(), b_.end()),
						quantiles_a = emp_equi_quantiles(a, NUM_QUANTILES),
						quantiles_b = emp_equi_quantiles(b, NUM_QUANTILES);

	double 	location, shape, size, d, 
			mean_a = mean(a), mean_b = mean(b),
			sd_a = sd(a), sd_b = sd(b), 
			quantile_cor_ab = cor(quantiles_a, quantiles_b);


	location 	= pow(mean_a - mean_b, 2);
	size 		= pow(sd_a - sd_b, 2);
	shape 		= 2 * sd_a * sd_b * (1 - quantile_cor_ab);
	d 			= location + size + shape;
	
	return Rcpp::List::create(
		Rcpp::Named("distance") = d,
		Rcpp::Named("location") = location,
		Rcpp::Named("size") = size,
		Rcpp::Named("shape") = shape
		);
}

//' squared_wass_approx
//'
//' Approximation of the squared wasserstein distance.
//' Calculation based on the mean squared difference between the equidistant 
//' quantiles of the two input vectors a and b.
//' As an approximation of the distribution, 1000 quantiles are computed for each vector.
//'
//' @param a Vector representing an empirical distribution under condition A 
//' @param b Vector representing an empirical distribution under condition B
//'	@param p exponent of the wasserstine distance
//' @return The approximated squared wasserstein distance between a and b
//'
//' @export
//[[Rcpp::export]]
double squared_wass_approx(const NumericVector & a_, const NumericVector & b_, const double p=1)
{
	const int			NUM_QUANTILES = 1000;
	double				wass_approx;
	const vector<double>	a(a_.begin(), a_.end()),
							b(b_.begin(), b_.end());
	vector<double>		quantiles_a(NUM_QUANTILES), 
						quantiles_b(NUM_QUANTILES),
						squared_quantile_diff(NUM_QUANTILES);

	quantiles_a				= emp_equi_quantiles(a, NUM_QUANTILES);
	quantiles_b				= emp_equi_quantiles(b, NUM_QUANTILES);
	const vector<double> quantile_diff	= quantiles_a - quantiles_b;
	squared_quantile_diff	= pow(quantile_diff, (double) 2.0);
	wass_approx 			= mean(squared_quantile_diff);

	return wass_approx;
}


//' wasserstein_metric
//'
//' @param a_ NumericVector representing an empirical distribution under condition A 
//' @param b_ NumericVector representing an empirical distribution under condition B
//'	@param p exponent of the wasserstine distance
//' @param wa_ NumericVector representing the weights of datapoints (interpreted as clusters) in a
//' @param wb_ NumericVector representing the weights of datapoints (interpreted as clusters) in b
//' @return The wasserstein (transport) distance between a and b
//'
//' @export
//[[Rcpp::export]]
double wasserstein_metric(NumericVector a_, 
						  NumericVector b_,
						  double p=1,
						  Nullable<NumericVector> wa_=R_NilValue, 
						  Nullable<NumericVector> wb_=R_NilValue) {
	//warning("wasserstein_metric will be deprecated in the next update. Use wasserstein() instead");
	vector<double> a(a_.begin(), a_.end());
	vector<double> b(b_.begin(), b_.end());
	sort(a.begin(), a.end());
	sort(b.begin(), b.end());

	// No weight vectors are given
	if (a.size() == b.size() && wa_.isNull() && wb_.isNull()) {

		Rcout << "a.size() = " << a.size() << ", b.size() = " << b.size() <<END;
		// compute root mean squared absolute difference of a and b
		// in R: mean(abs(sort(b) - sort(a))^p)^(1/p)
		vector<double> sq_abs_diff = pow(
		  abs(b - a),
		  p);
		double mrsad = pow((double) mean(sq_abs_diff), (double) 1.0/p);
		return mrsad;

	}

	// At least one weight vector is given
	// If only one weight vector is undefined, set all its weights to 1
	double default_weight = 1.0;
	std::vector<double> wa(a.size(), default_weight);
	std::vector<double> wb(b.size(), default_weight);

	if (!wa_.isNull()) {
	  NumericVector dumpwa = wa_.get();
		wa = vector<double>(dumpwa.begin(), dumpwa.end());
	}
	if (!wb_.isNull()) {
	  NumericVector dumpwb =wb_.get();
		wb = vector<double>(dumpwb.begin(), dumpwb.end());
	}

	// normalize the weights to add up to 1
	vector<double> ua(wa.size());
	ua = wa / sum(wa);
	ua.pop_back();
	vector<double> ub(wb.size());
	ub = wb / sum(wb);
	ub.pop_back();

	// cumulative distribution without the last value
	// => last value will be considered as an open interval to Infinity
	vector<double> cua(ua.size());
	cua = cumSum(ua);
	vector<double> cub(ub.size());
	cub = cumSum(ub);

	vector<int> a_rep = interval_table(cub, cua, 0);
	vector<int> b_rep = interval_table(cua, cub, 0);
	a_rep = a_rep + 1;
	b_rep = b_rep + 1;

	//Rcout << "a_rep = "; for (int ea : a_rep) { Rcout << ea; } Rcout << END;
	//Rcout << "b_rep = "; for (int eb : b_rep) { Rcout << eb; } Rcout << END;

	int len_a_weighted = sum(a_rep);
	int len_b_weighted = sum(b_rep);
	//Rcout << "len_b_weighted = " << sum(b_rep) << END;
	//Rcout << "len_a_weighted = " << sum(a_rep) << END;

	vector<double> a_weighted(len_b_weighted);
	vector<double> b_weighted(len_a_weighted);
	a_weighted = rep_weighted(a, a_rep);
	b_weighted = rep_weighted(b, b_rep);

	vector<double> uu0(cua.size() + cub.size());
	vector<double> uu1(cua.size() + cub.size());
	//Rcout << "AFTER INITIALIZATION: uu1.size() = " << uu1.size() << ", uu0.size() = " << uu0.size() <<END;

	uu0 = concat(cua, cub);
	uu0.insert(uu0.begin(),0);
	uu1 = concat(cua, cub);
	uu1.push_back(1);

	double wsum = 0.0;
	double areap = 0.0;

	//Rcout << "AFTER CONCATENATING SOMETHING: uu1.size() = " << uu1.size() << ", uu0.size() = " << uu0.size() <<END;
	//Rcout << "b_weighted.size() = " << b_weighted.size() << ", a_weighted.size() = " << a_weighted.size() <<END;

	wsum = sum((uu1 - uu0) * pow(abs(b_weighted - a_weighted), p));
	areap = pow((double) wsum, (double) (1/p));

	return areap;

}



/*=============================================

			EXPORTS FOR TESTING IN R

==============================================*/ 

// [[Rcpp::export]]
NumericVector add_test_export(NumericVector & x_, NumericVector & y_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());

	result = x + y;
	NumericVector output(result.begin(), result.end());
	return output;
}


// [[Rcpp::export]]
NumericVector add_test_export_sv(NumericVector & x_, double & summand_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					result(x.size());
	const double summand = summand_;

	result = x + summand;
	NumericVector output(result.begin(), result.end());
	return output;
}


// [[Rcpp::export]]
NumericVector multiply_test_export(NumericVector & x_, NumericVector & y_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());

	result = x * y;
	NumericVector output(result.begin(), result.end());
	return output;
}


// [[Rcpp::export]]
NumericVector multiply_test_export_sv(NumericVector & x_, double & factor_)
{

	vector<double> 	x(x_.begin(), x_.end()),
					result(x.size());
	const double factor = factor_;

	result = x * factor;
	NumericVector output(result.begin(), result.end());
	return output;
}


// [[Rcpp::export]]
NumericVector pow_test_export(NumericVector & x_, const double & exp)
{
	vector<double>	x(x_.begin(), x_.end()),
					result(x.size());
	
	result = pow(x, exp);

	NumericVector 	output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
NumericVector abs_test_export(NumericVector & x_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					result(x.size());
	
	result = abs(x);

	NumericVector 	output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
double sum_test_export(NumericVector & x_)
{
	vector<double> 	x(x_.begin(), x_.end());

	double result = sum(x);

	return result;
}

// [[Rcpp::export]]
NumericVector subtract_test_export(NumericVector & x_, NumericVector & y_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());

	result = x - y;
	
	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
NumericVector divide_test_export_sv(NumericVector & x_, double & y_)
{
	vector<double>	x(x_.begin(), x_.end()),
					result(x.size());
	const double 	y = y_;

	result = x / y;

	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
NumericVector divide_test_export_vectors(NumericVector & x_, NumericVector & y_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size());

	result = x / y;

	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
double mean_test_export(NumericVector & x_)
{
	vector<double> 	x(x_.begin(), x_.end());
	double 			result = mean(x);
	return result;
}

// [[Rcpp::export]]
double sd_test_export(NumericVector & x_)
{
	vector<double> 	x(x_.begin(), x_.end());
	double 			result = sd<double>(x);
	return result;	
}

// [[Rcpp::export]]
NumericVector cumSum_test_export(NumericVector & x_, int last_index=0)
{

	vector<double> 	x(x_.begin(), x_.end()),
					result;
	
	result = cumSum(x, last_index);
	
	NumericVector 	output(result.begin(), result.end());
	return output;	
}

// [[Rcpp::export]]
double cor_test_export(NumericVector x_, NumericVector y_)
{	

	vector<double> 	x(x_.size()),
					y(y_.size());
	
	for (int i=0; i<x.size(); i++) { x[i] = x_[i];}
	for (int j=0; j<y.size(); j++) { y[j] = y_[j];}

	double result = cor(x, y);
	
	return result;
}

// [[Rcpp::export]]
NumericVector emp_equi_quantiles_test_export(NumericVector & x_, const int & K)
{
	const 	vector<double> 	x(x_.begin(), x_.end());
			vector<double>	result(x.size());

	result = emp_equi_quantiles(x,K);
	
	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
NumericVector rep_weighted_test_export(NumericVector & x_, NumericVector & weights_)
{
	vector<int>		weights(weights_.begin(), weights_.end());
	vector<double> 	x(x_.begin(), x_.end()),
					result(sum(weights));

	result = rep_weighted(x, weights);

	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
NumericVector concat_test_export(NumericVector & x_, NumericVector & y_)
{
	vector<double> 	x(x_.begin(), x_.end()),
					y(y_.begin(), y_.end()),
					result(x.size() + y.size());

	result = concat(x, y);

	NumericVector output(result.begin(), result.end());
	return output;
}

// [[Rcpp::export]]
IntegerVector interval_table_test_export(	NumericVector & data_,
											NumericVector & breaks_,
											const int & default_freq=0)
{
	vector<double> 	x(data_.begin(), data_.end()),
					y(breaks_.begin(), breaks_.end());
	vector<int>		result(x.size() + y.size());

	result = interval_table(x, y, default_freq);

	IntegerVector output(result.begin(), result.end());
	return output;
}



/*
	Program entry point for debugging
*/
int main() 
{

	return 0;
}