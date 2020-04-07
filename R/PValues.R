
#'Compute combined p-value using Fisher's method
#'
#'Aggregates two p-values into a combined p-value according to Fisher’s method 
#'
#'@details Aggregates two p-values into a combined p-value according to
#' Fisher’s method
#'
#'@param x vector of the two p-values that are to be aggregated
#'
#'@return The combined p-value 
#'
.fishersCombinedPval <- function(x) {
    if(sum(is.na(x)) == 0) {
        p.comb <- pchisq(-2 * sum(log(x)), df=2*length(x), lower.tail=FALSE)
    } else if (sum(is.na(x)) == 1){
        p.comb <- x[!is.na(x)]
    } else {
        p.comb <- NA
    }
    return(p.comb)
}


#'Compute combined p-value using Fisher's method for several pairs of p-values
#'
#'For a given set of \eqn{N} pairs of p-values, aggregates each respective pair of 
#'p-values into a combined p-value according to Fisher’s method 
#'
#'@details For a given set of \eqn{N} pairs of p-values, aggregates each respective 
#' pair of p-values into a combined p-value according to Fisher’s method.
#' Applies the function \code{.fishersCombinedPval} to a whole set of \eqn{N} pairs of
#' p-values. 
#'
#'@param r vector of length \eqn{N} of the p-values corresponding to the first test
#'@param s vector of length \eqn{N} of the p-values corresponding to the second test
#'
#'@return A vector of length \eqn{N} of the combined p-values
#'
.combinePVal <- function(r,s){
    apply(cbind(r,s), 1, function(x)
        .fishersCombinedPval(x))
}


#TODO: Use .gpdFittedPValue again, once the issues with eva are fixed
#' Compute p-value based on generalized Pareto distribution fitting
#'
#' Computes a p-value based on a generalized Pareto distribution (GPD) fitting. This procedure may be used in the semi-parametric, 2-Wasserstein distance-based test to estimate small p-values accurately, instead of obtaining the p-value from a permutation test.
#' 
#' @param val value of a specific test statistic, based on original group labels
#' @param distr.ordered vector of values of the test statistic obtained by repeatedly permuting the original group labels
#'@return A vector of three, see Schefzik et al. (2020) for details:
#' \itemize{
#' \item pvalue.gpd: p-value corresponding to the GPD fitting
#' test
#' \item ad.pval: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' \item N.exc: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit, i.e. p-value of Anderson-Darling test
#' \eqn{\geq 0.05}
#' }
#'
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2020). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
NULL
#.gpdFittedPValue <- function(val, distr.ordered) {
#    
#    # list of possible exceedance thresholds (decreasing)
#    poss.exc.num <- seq(from=250, to=10, by=-10)
#    
#    bsn <- length(distr.ordered)
#    r <- 1
#    repeat {
#        
#        # set threshold for exceedance according to paper
#        N.exc <- poss.exc.num[r]
#        
#        # compute set of N.exc exceedances
#        exceedances <- distr.ordered[seq_len(N.exc)]
#        
#        # check whether the N.exc largest permutation values follow 
#        # a GPD using an Anderson-Darling test
#        gpd.ad.check <- gpdAd(exceedances)
#        ad.pval <- gpd.ad.check$p.value
#        
#        r <- r + 1
#        if (ad.pval > 0.05) {break}
#    }
#    
#    # calculate exceedance threshold for so-obtained N.exc
#    t.exc <- (distr.ordered[N.exc] 
#              + distr.ordered[N.exc+1]) / 2
#    
#    # fit GPD distribution to the exceedances using maximum 
#    # likelihood estimation
#    gpd.fit <- gpdFit(  data=distr.ordered,
#                        threshold=t.exc,
#                        method="mle")
#    
#    # extract fitted parameters
#    fit.scale <- as.numeric(gpd.fit$par.ests[1])
#    fit.shape <- as.numeric(gpd.fit$par.ests[2])
#    
#    # compute GPD p-value (see paper)
#    pvalue.gpd <- (N.exc / bsn) * (1 - pgpd(q=val-t.exc,
#                                            loc=0,
#                                            scale=fit.scale,
#                                            shape=fit.shape))
#    
#    pvalue.gpd <- as.numeric(pvalue.gpd)
#    pvalue.wass <- c("pvalue.gpd"=pvalue.gpd,
#                     "ad.pval"=ad.pval,
#                     "N.exc"=N.exc)
#    return(pvalue.wass)
#}
