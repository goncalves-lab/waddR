
#' .relativeError
#'
#' Computes the relative error between two numericals (in the sense of
#' is.numeric) x and y.
#' If x and y are vectors, it is assumed that length(x) == length(y).
#' 
#' The relative error e is defined as: e = | 1 -  x/y |
#' 
#' @param x numerical (in the sense of is.numeric)
#' @param y numerical (in the sense of is.numeric)
#' 
#' @return The relative error between x and y
#'
.relativeError <- function(x,y) {
    if ((x == y) || ((x == 0) && (y == 0)))
        rel.error <- 0
    else {
        rel.error <- abs(1 - (x / y))  
    }
    return(rel.error)
}


#'.quantileCorrelation
#'
#' Computes the quantile-quantile correlation of x and y, using the quantile
#' type 1 implementation in R and Pearson correlation.
#' 
#' @param x numeric vector representing distribution x
#' @param y numeric vector representing distribution y
#' @param pr probabilities for computing the quantiles of x and y. 
#'  By default is set to 1000 equidistant quantiles starting at q_1 = 0.5/1000.
#' @return quantile-quantile correlation of x and y
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