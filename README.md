# Differential Expression Detection in scRNA seq Data

`diffexpR` is an R package for detecting differential expression in scRNA Data from two different conditions using the Wasserstein distance and a differential proportion test.

# Installation

## Enable Bioconductor Repositories
This project imports other packages that are hosted on the Bioconductor Software repository.
It is required to select a Bioconductor mirror, or the installation of `diffexpR` will fail.
Repositories where R queries for packages upon calling `install.packages()`can be set from within an R session can be set by calling: `setRepositories()`.
Entering `1 2` selects the CRAN and Bioconductor Software repositories.

For the best experience when using and testing the package, it is also recommended to have the `devtools` package installed.

## From Github

The package can be installed with the `install_github` utility from the `devtools` package.
Simply run `install_github("goncalves-lab/diffexpR")` from an R session.

## Via Package Repository

The package is not yet available via the official repositories.


# Running `diffexpR`

## Wasserstein Distance

## Wasserstein Test

## Running on single cell RNA sequencing data

# Running Tests

Tests can be run by calling `test()` from the `devtools` package.
All tests are implemented using the `testthat` package and reside in `tests/testhat`

# Documentation

# References

