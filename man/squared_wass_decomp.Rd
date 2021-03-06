% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{squared_wass_decomp}
\alias{squared_wass_decomp}
\title{Compute the squared 2-Wasserstein distance based on a decomposition}
\usage{
squared_wass_decomp(x, y)
}
\arguments{
\item{x}{sample (vector) representing the distribution of condition \eqn{A}}

\item{y}{sample (vector) representing the distribution of condition \eqn{B}}
}
\value{
A list of 4:
\itemize{
\item distance: the sum location+size+shape
\item location: location part in the decoposition of the 2-Wasserstein distance
\item size: size part in the decoposition of the 2-Wasserstein distance
\item shape: shape part in the decoposition of the 2-Wasserstein distance
}
}
\description{
Computes the squared 2-Wasserstein distance between two vectors based on a decomposition into location, size and shape terms.
For a detailed description of the (empirical) calculation of the invoved quantities, see Schefzik et al. (2020).
}
\examples{
set.seed(24)
x<-rnorm(100)
y1<-rnorm(150)
y2<-rexp(150,3)
y3<-rpois(150,2)

squared_wass_decomp(x,y1)
squared_wass_decomp(x,y2)
squared_wass_decomp(x,y3)

}
\references{
Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
}
\seealso{
See the functions \code{wasserstein_metric} and \code{squared_wass_approx} for
alternative implementations of the 2-Wasserstein distance
}
