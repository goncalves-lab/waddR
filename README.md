# Statistical tests for detecting differential distributions based on the 2-Wasserstein distance

The `waddR` package offers statistical tests based on the 2-Wasserstein distance for detecting and characterizing differences between two distributions given in the form of samples. Functions for calculating the 2-Wasserstein distance and testing for differential distributions are provided, as well as specifically tailored test for differential expression in single-cell RNA sequencing data.

The package provides tools to address the following tasks:
1. Computation of the 2-Wasserstein distance 
2. Two-sample tests to check for differences between two distributions
3. Detection of differential gene expression distributions in single-cell RNA sequencing data

## Installation

### Requirements

* R >= 3.6.0

### Via Package Repository

Available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/waddR.html):
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("waddR")
```

### From Github

The latest package version can be installed from Github using `BiocManager`:

```
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("goncalves-lab/waddR")
```

## Running Tests

Tests can be run by calling `test()` from the `devtools` package.
All tests are implemented using the `testthat` package and reside in `tests/testhat`



## Using `waddR`

### 2-Wasserstein distance functions

The 2-Wasserstein distance is a metric to quantify the distance between two
distributions, representing two different conditions A and B. The `waddR` package
specifically considers the squared 2-Wasserstein distance, which
offers a decomposition into location, size, and shape terms, thus providing a characterization of potential differences.

The `waddR` package offers three functions to calculate the 2-Wasserstein
distance, which are implemented in C++ and exported to R with Rcpp for
faster computation.
The function `wasserstein_metric` is a C++ reimplementation of the
function `wasserstein1d` from the R package `transport`.
The functions `squared_wass_approx` and `squared_wass_decomp` compute
approximations of the squared 2-Wasserstein distance, with `squared_wass_decomp`
also returning the decomposition terms for location, size, and shape. 

See `?wasserstein_metric`, `?squared_wass_aprox`, and `?squared_wass_decomp`, as well as the accompanying paper Schefzik et al. (2020).

### Testing for differences between two distributions

The `waddR` package provides two testing procedures using the 2-Wasserstein distance
to test whether two distributions given in the form of samples are
different by specifically testing the null hypothesis of no difference against the
alternative hypothesis that the two distributions are different.

The first, semi-parametric (SP), procedure uses a permutation-based test combined with a generalized Pareto distribution approximation to 
estimate small p-values accurately.

The second procedure uses a test based on asymptotic theory (ASY) which is
valid only if the samples can be assumed to come from continuous
distributions.

See `?wasserstein.test` for more details.

### Testing for differences between two distributions in the context of single-cell RNA sequencing data:

The `waddR` package provides an adaptation of the
semi-parametric testing procedure based on the 2-Wasserstein distance
which is specifically tailored to identify differential distributions in
single-cell RNA-seqencing (scRNA-seq) data. In particular, a two-stage
(TS) approach is implemented that takes account of the specific
nature of scRNA-seq data by separately testing for differential
proportions of zero gene expression (using a logistic regression model)
and differences in non-zero gene expression (using the semi-parametric
2-Wasserstein distance-based test) between two conditions.

See `?wasserstein.sc` and `?testZeroes` for more details.

## Documentation

We have included detailed examples of how to use the functions provided with
`waddR` in our vignettes.
They are available online [here](https://github.com/goncalves-lab/waddR) 
*(update this link once it is final)* or from an R session with the
following command: 
`browseVignettes("waddR")`

# References

Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.

