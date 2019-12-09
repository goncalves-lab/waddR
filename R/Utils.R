
#' Compute relative error
#'
#' Computes the relative error between two numbers \eqn{x} and \eqn{y}.
#' 
#' The relative error \eqn{e} is defined as \eqn{e = | 1 -  x/y |}
#' 
#' @param x a number
#' @param y a number
#' 
#' @return The relative error between \eqn{x} and \eqn{y}
#'
.relativeError <- function(x,y) {
    if ((x == y) || ((x == 0) && (y == 0)))
        rel.error <- 0
    else {
        rel.error <- abs(1 - (x / y))  
    }
    return(rel.error)
}


#' Compute the quantile-quantile correlation
#'
#' Computes the quantile-quantile correlation of two samples \eqn{x} and \eqn{y}, using the quantile
#' type 1 implementation in R and the Pearson correlation.
#' 
#' @param x numeric vector 
#' @param y numeric vector 
#' @param pr levels at which the quantiles of \eqn{x} and \eqn{y} are computed; by default, 1000 equidistant quantiles at levels \eqn{\frac{k-0.5}{1000}}, where \eqn{k=1,\ldots,1000}, are used 
#' @return quantile-quantile correlation of \eqn{x} and \eqn{y}
#' 
.quantileCorrelation <- function(x, y, pr=NULL) {
    stopifnot(length(x) != 0, length(y) != 0)
    
    if (is.null(pr)){
        pr <- (seq_len(1000) - 0.5) / 1000
    }
    
    # quantiles of type 1
    quant.x <- quantile(x, probs=pr, type=1)
    quant.y <- quantile(y, probs=pr, type=1)
    if (sd(quant.x) !=0 & sd(quant.y) != 0){
        return(cor(quant.x, quant.y, method="pearson"))
    } else {
        return(0)  
    }   
}
