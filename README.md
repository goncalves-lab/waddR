# Differential Expression Detection in scRNA seq Data

`waddR` is an R package that provides a Wasserstein distance based statistical
test for detecting and describing differential distributions in one-dimensional
data.
Functions for wasserstein distance calculation, differential distribution
testing, and a specialized test for differential expression in scRNA data are
provided.

The Wasserstein package offers utilities for three distinct use cases:
    * Computation of the 2-Wasserstein distance 
    * Check for differences between two distributions
    * Detect differential gene expression distributions in scRNAseq data

# Installation

## Via Package Repository

The package is in development and not yet available via an official package
repository.

## From Github

The latest package version can be installed with the `install_github` utility from the 
`devtools` package.
Simply run `devtools::install_github("goncalves-lab/waddR")` from an R session.

# Using `waddR`

## Wasserstein Distance functions

The 2-Wasserstein distance is a metric to describe the distance between two
distributions, representing two diferent conditions A and B. This package
specifically considers the squared 2-Wasserstein distance d := W^2 which
offers a decomposition into location, size, and shape terms.

The package `waddR` offers three functions to calculate the 2-Wasserstein
distance, all of which are implemented in Cpp and exported to R with Rcpp for
better performance.
The function `wasserstein_metric` is a Cpp reimplementation of the
function `wasserstein1d` from the package `transport` and offers the most exact
results.
The functions `squared_wass_approx` and `squared_wass_decomp` compute
approximations of the squared 2-Wasserstein distance with `squared_wass_decomp`
also returning the decomosition terms for location, size, and shape. 


## Two-Sample Testing

This package provides two testing procedures using the 2-Wasserstein distance
to test whether two distributions F_A and F_B given in the form of samples are
different ba specifically testing the null hypothesis H0: F_A = F_B against the
alternative hypothesis H1: F_A != F_B.

The first, semi-parametric (SP), procedure uses a test based on permutations
combined with a generalized pareto distribution approximation to estimate small
pvalues accurately.

The second procedure (ASY) uses a test based on asymptotic theory which is
valid only if the samples can be assumed to come from continuous
distributions.

See the documentation of these functions `?wasserstein.test`,
`?wasserstein.test.sp`, `?wasserstein.test.asy` for more details.

## Single Cell Test: The waddR package provides an adaptation of the

semi-parametric testing procedure based on the 2-Wasserstein distance
which is specifically tailored to identify differential distributions in
single-cell RNA-seqencing (scRNA-seq) data. In particular, a two-stage
(TS) approach has been implemented that takes account of the specific
nature of scRNA-seq data by separately testing for differential
proportions of zero gene expression (using a logistic regression model)
and differences in non-zero gene expression (using the semi-parametric
2-Wasserstein distance-based test) between two conditions.

See the documentation of the Single Cell testing function `?wasserstein.sc`
and the test for zero expression levels `?testZeroes` for more details.

# Documentation

We have included detailed examples of how to use all functions provided with
`waddR` in our vignettes.
They are available online [here](https://github.com/goncalves-lab/waddR) 
(update this link!) or from an R session with the following command:
`browseVignettes("waddR")`

# Running Tests

Tests can be run by calling `test()` from the `devtools` package.
All tests are implemented using the `testthat` package and reside in `tests/testhat`


# References

