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

This an example on how to analyse the expression of 18 immune cell in two hypothetical conditions.
For an example dataset we use GSE41265, which can be downloaded from the [conquer database](http://imlspenticton.uzh.ch:3838/conquer/). 
We will assume that the file `GSE41265.rds` is already available in working directory.
To read it, install the packages `MultiAssayExperiment` and `SummarizedExperiment` from BioConductor.

### Loading all required Libraries and Data

For a more detailed description of how to handle datasets from the conquer database, see the [conqer R tutorial](https://github.com/markrobinsonuzh/conquer/blob/master/shiny-download/tutorial.md).
```
# load libraries and the .rds data file
suppressPackageStartupMessages(library(diffexpR))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(MultiAssayExperiment))
gse41265 <- readRDS("GSE41265.rds")

# extract gene abundance from the dataset
gene_abd <- experiments(gse41265)[[1]]

# we are interested in looking at the "Transcripts per Million" matrix
data <- assays(gse41265_gene)[["TPM"]]

# for faster runtime, we randomly select 100 of the 45,686 transcripts
rand_rows <- sort(sample(dim(data)[1], 100))
data <- data[rand_rows, ]

# for this demonstration, we define arbitrary conditions
conditions <- c(rep(1,9), rep(2,9))

# we want a one step test of all transcripts in the two conditions
test_results <- wasserstein.sc(data, conditions, 24, 1000, "TS")

# for every gene, a test result is printed, describing how the distribution differs in size,
# shape and location between the two conditions, while providing p-values and errors
head(results)
                      d.transport d.transport^2     d.comp^2     d.comp
ENSMUSG00000000215.10   0.0000000     0.0000000    0.0000000  0.0000000
ENSMUSG00000008999.7    0.6087344     0.3705575    0.4049929  0.6363905
ENSMUSG00000013523.13  55.1313604  3039.4669034 3371.4203244 58.0639331
ENSMUSG00000017204.4    1.0423733     1.0865420    1.2220874  1.1054806
ENSMUSG00000018401.17   4.5537914    20.7370158   22.8496752  4.7801334
ENSMUSG00000018634.10   0.0000000     0.0000000    0.0000000  0.0000000
                          location         size      shape       rho  pval
ENSMUSG00000000215.10 0.000000e+00    0.0000000 0.00000000 0.0000000 1.000
ENSMUSG00000008999.7  9.524064e-02    0.2988548 0.01089746 0.9837274 0.659
ENSMUSG00000013523.13 3.838363e+02 2984.8860996 2.69795548 0.8202325 0.444
ENSMUSG00000017204.4  1.807983e-03    0.8105763 0.40970319 0.8391326 0.469
ENSMUSG00000018401.17 3.835737e+00   18.9428596 0.07107891 0.9599687 0.296
ENSMUSG00000018634.10 0.000000e+00    0.0000000 0.00000000 0.0000000 1.000
                      p.ad.gpd N.exc perc.loc perc.size perc.shape decomp.error
ENSMUSG00000000215.10       NA    NA      NaN       NaN        NaN   0.00000000
ENSMUSG00000008999.7        NA    NA    23.52     73.79       2.69   0.03443539
ENSMUSG00000013523.13       NA    NA    11.39     88.53       0.08 331.95342096
ENSMUSG00000017204.4        NA    NA     0.15     66.33      33.52   0.13554541
ENSMUSG00000018401.17       NA    NA    16.79     82.90       0.31   2.11265936
ENSMUSG00000018634.10       NA    NA      NaN       NaN        NaN   0.00000000
                      pval.adj
ENSMUSG00000000215.10        1
ENSMUSG00000008999.7         1
ENSMUSG00000013523.13        1
ENSMUSG00000017204.4         1
ENSMUSG00000018401.17        1
ENSMUSG00000018634.10        1

```

# Running Tests

Tests can be run by calling `test()` from the `devtools` package.
All tests are implemented using the `testthat` package and reside in `tests/testhat`

# Documentation

# References

