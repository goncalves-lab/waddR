% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.r
\name{.brownianBridgeEmpcdf}
\alias{.brownianBridgeEmpcdf}
\title{Compute value of the asymptotic CDF occuring in the asymptotic theory-based test}
\usage{
.brownianBridgeEmpcdf(v)
}
\arguments{
\item{v}{a number}
}
\value{
Value at \code{v} of the asymptotic CDF
}
\description{
Computes the values of the cumulative distribution function (CDF) of the integral over the squared standard Brownian bridge in the unit interval, where the computation is based on Monte Carlo simulations.
This CDF occurs as an asymptotic distribution in the asymptotic theory-based test using the 2-Wasserstein distance, see Schefzik et al. (2020) for details. 
It is used to determine the corresponding p-values in the function \code{.wassersteinTestAsy}.
}
\references{
Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
}
