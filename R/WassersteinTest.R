
#' Permutation procedure that calculates the squared 2-Wasserstein distances for
#' random shuffles of two input samples representing two conditions
#'
#' Permutation procedure that calculates the squared 2-Wasserstein distances for
#' random shuffles of two input samples representing two conditions A and B, respectively (i.e. the elements of the pooled sample vector are randomly assigned to either condition A or condition B, where the sample sizes are as in the original samples).
#'
#' @param x sample (vector) representing the distribution of
#' condition A
#' @param y sample (vector) representing the distribution of
#' condition B
#' @param permnum number of permutations to be
#'  performed
#'
#' @return Vector with squared 2-Wasserstein distances computed for random
#'  shuffles of the two input samples
#'  
.wassPermProcedure <- function(x, y, permnum) {
    z <- c(x,y)
    
    shuffle <- permutations(z, num_permutations=permnum)
    return(apply(   shuffle, 2, 
                    function (k) {
                        return( wasserstein_metric(
                                    k[seq_len(length(x))],
                                    k[seq((length(x)+1), length(z))],
                                    p=2) **2)
                    }))
}


#' Semi-paramteric test using the 2-Wasserstein distance to check for differential distributions
#' 
#' Two-sample test to check for differences between two distributions
#' using the 2-Wasserstein distance: Semi-parametric
#' implementation using a permutation test with a generalized Pareto
#' distribution (GPD) approximation to estimate small p-values accurately
#' 
#' This is the semi-parametric version of \code{wasserstein.test}, for the
#' asymptotic theory-based procedure see \code{.wassersteinTestAsy}.
#' 
#'@details Details concerning the permutation testing procedure with GPD
#' approximation to estimate small p-values accurately can be found in Schefzik
#' et al. (2019).
#'
#'@param x sample (vector) representing the distribution of
#' condition A
#'@param y sample (vector) representing the distribution of
#' condition B
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'@return A vector of 15, see Schefzik et al. (2019) for details:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed by
#' quantile approximation
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#' computed by quantile approximation
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#' computed by decomposition approximation
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#' decomposition approximation
#' \item location: location term in the decomposition of the squared
#' 2-Wasserstein distance between the two samples
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item rho: correlation coefficient in the quantile-quantile plot
#' \item pval: p-value of the semi-parametric 2-Wasserstein distance-based
#' test
#' \item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' (otherwise NA). NOTE: GPD fitting is currently not supported!
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
<<<<<<< HEAD
#' required to obtain a good GPD fit (i.e. p-value of Anderson-Darling test
#' greater or eqaul to 0.05) (otherwise NA) NOTE: GPD fitting is currently 
#' not supported!
#' \item perc.loc: fraction (in %) of the location part with respect to the
=======
#' required to obtain a good GPD fit, i.e. p-value of Anderson-Darling test
#' \eqn{\geq 0.05} (otherwise NA)
#' \item perc.loc: fraction (in \%) of the location part with respect to the
>>>>>>> 988f3e9401c07aae38f7c082377ae38d93cedd38
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.size: fraction (in \%) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.shape: fraction (in \%) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#' distance obtained by the quantile approximation and the squared
#' 2-Wasserstein distance obtained by the decomposition approximation
#' }
#'
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
.wassersteinTestSp <- function(x, y, permnum=10000){
    stopifnot(permnum>0)
    if (length(x) !=0 & length(y) != 0){

        # wasserstein distance between the samples
        value <- wasserstein_metric(x, y, p=2)
        value.sq <- value**2

        # permutation procedure to calculate the wasserstein distances of
        # random shuffles of x and y
        bsn <- permnum
        wass.values <- .wassPermProcedure(x, y, bsn)
        wass.values.ordered <- sort(wass.values, decreasing=TRUE)

        # computation of an approximative p-value
        num.extr <- sum(wass.values >= value.sq)
        pvalue.ecdf <- num.extr/bsn
        pvalue.ecdf.pseudo <- (1 + num.extr) / (bsn + 1)

        # gdp fitting needed
        pvalue.wass <- pvalue.ecdf
        pvalue.gdpfit <- NA
        N.exc <- NA
        env <- environment()
        if (num.extr < 10) {
            #TODO: Use .gdpFittedPValue again, once the issues with eva are fixed
            #tryCatch({
            #        res <- .gdpFittedPValue(value.sq,
            #                                wass.values.ordered)
            #        assign("pvalue.wass", unname(res["pvalue.gpd"]), env)
            #        assign("pvalue.gdpfit", unname(res["ad.pva"]), env)
            #        assign("N.exc", unname(res["N.exc"]), env)
            #    }, error=function(...) {
            #        assign("pvalue.wass", pvalue.ecdf.pseudo, env)
            #        assign("pvalue.gdpfit", NA, env)
            #        assign("N.exc", NA, env)
            #    })
            
            # For now, just use pseudo pvalues
            assign("pvalue.wass", pvalue.ecdf.pseudo, env)
        }

        # correlation of quantile-quantile plot
        rho.xy <- .quantileCorrelation(x, y)

        # decomposition of wasserstein distance
        wass.comp <- squared_wass_decomp(x, y)
        location <- wass.comp$location
        size <- wass.comp$size
        shape <- wass.comp$shape
        d.comp.sq <- wass.comp$distance
        if (is.na(d.comp.sq)) {
            d.comp.sq <- sum(na.exclude(location, shape, size))
        }
        d.comp <- sqrt(d.comp.sq)
        perc.loc <- round(((location / d.comp.sq) * 100), 2)
        perc.size <- round(((size / d.comp.sq) * 100), 2)
        perc.shape <- round(((shape / d.comp.sq)*100), 2)
        decomp.error <- .relativeError(d.comp.sq, value.sq)

        output <- c("d.wass"=value, "d.wass^2"=value.sq, "d.comp^2"=d.comp.sq,
                    "d.comp"=d.comp, "location"=location, "size"=size,
                    "shape"=shape, "rho"=rho.xy, "pval"=pvalue.wass,
                    "p.ad.gdp"=pvalue.gdpfit, "N.exc"=N.exc,
                    "perc.loc"=perc.loc, "perc.size"=perc.size,
                    "perc.shape"=perc.shape, "decomp.error"=decomp.error)

    } else { output <- c("d.wass"=NA, "d.wass^2"=NA, "d.comp^2"=NA,
                         "d.comp"=NA, "location"=NA, "size"=NA,
                         "shape"=NA, "rho"=NA, "pval"=NA,
                         "p.ad.gdp"=NA, "N.exc"=NA,
                         "perc.loc"=NA, "perc.size"=NA,
                         "perc.shape"=NA, "decomp.error"=NA) }

    return(output)
}


#'Asymptotic theory-based test using the 2-Wasserstein distance to check for differential distributions
#' 
#' Two-sample test to check for differences between two distributions
#' using the 2-Wasserstein distance: Implementation using a test
#' based on asymptotic theory
#'
#' This is the asymptotic version of \code{wasserstein.test}, for the
#' semi-parametric procedure see \code{.wassersteinTestSp}.
#'
#'@details Details concerning the testing procedure based on asymptotic theory
#' can be found in Schefzik et al (2019).
#'
#'@param x sample (vector) representing the distribution of
#' condition A
#'@param y sample (vector) representing the distribution of
#' condition B
#'
#'@return A vector of 13, see Schefzik et al. (2019) for details:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed
#' by quantile approximation
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#' computed by quantile approximation
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#' computed by decomposition approximation
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#' decomposition approximation
#' \item location: location term in the decomposition of the squared
#' 2-Wasserstein distance between the two samples
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item rho: correlation coefficient in the quantile-quantile plot
#' \item pval: p-value of the 2-Wasserstein distance-based test using
#' asymptotic theory
#' \item perc.loc: fraction (in \%) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.size: fraction (in \%) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.shape: fraction (in \%) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#' distance obtained by the quantile approximation and the squared
#' 2-Wasserstein distance obtained by the decomposition approximation
#' }
#'
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
.wassersteinTestAsy <- function(x, y){

    if (length(x) != 0 & length(y) != 0) {

        value <- wasserstein_metric(x, y, p=2)
        value.sq <- value **2 

        # compute p-value based on asymptotoc theory (brownian bridge)
        pr <- seq(from=0, to=1, by=1/10000)
        empir.cdf.y <- ecdf(y)

        parts.trf <- ((empir.cdf.y(quantile(x, probs=pr, type=1))) - pr) **2

        trf.int <- (1 / length(pr)) * sum(parts.trf)
        test.stat <- (length(x) * length(y) 
                    / (length(x) + length(y))) * trf.int

        # p-value
        pvalue.wass <- 1 - .brownianBridgeEmpcdf(test.stat)

        # correlation of quantile-quantile plot
        rho.xy <- .quantileCorrelation(x, y)

        # decomposition of wasserstein distance
        wass.comp <- squared_wass_decomp(x, y)
        location <- wass.comp$location
        size <- wass.comp$size
        shape <- wass.comp$shape
        d.comp.sq <- wass.comp$distance
        if (is.na(d.comp.sq)) {
            d.comp.sq <- sum(na.exclude(location, shape, size))
        }
        d.comp <- sqrt(d.comp.sq)
        perc.loc <- round(((location / d.comp.sq) * 100), 2)
        perc.size <- round(((size / d.comp.sq) * 100), 2)
        perc.shape <- round(((shape / d.comp.sq)*100), 2)

        decomp.error <- .relativeError(value.sq, d.comp.sq)
        
        output <- c("d.wass"=value, "d.wass^2"=value.sq, "d.comp^2"=d.comp.sq,
                    "d.comp"=d.comp, "location"=location, "size"=size,
                    "shape"=shape, "rho"=rho.xy, "pval"=pvalue.wass,
                    "perc.loc"=perc.loc, "perc.size"=perc.size,
                    "perc.shape"=perc.shape, "decomp.error"=decomp.error)
    } else { output <-c("d.wass"=NA, "d.wass^2"=NA, "d.comp^2"=NA,
                        "d.comp"=NA, "location"=NA, "size"=NA,
                        "shape"=NA, "rho"=NA, "pval"=NA,
                        "perc.loc"=NA, "perc.size"=NA,
                        "perc.shape"=NA, "decomp.error"=NA)
    }
    return(output)
}


#'Two-sample test to check for differences between two distributions
#'using the 2-Wasserstein distance
#'
#'Two-sample test to check for differences between two distributions
#'using the 2-Wasserstein distance, either using the
#'semi-parametric permutation testing procedure with a generalized Pareto distribution (GPD) approximation to
#'estimate small p-values accurately or the test based on asymptotic theory
#'
#'@name wasserstein.test
#'@details Details concerning the two testing procedures (i.e. the semi-paarmetric permutation
#' testing procedure with GPD approximation and the test based on asymptotic theory) can be found in
#' Schefzik et al. (2019).
#'
#'@param x sample (vector) representing the distribution of
#' condition A
#'@param y sample (vector) representing the distribution of
#' condition B
#'@param method testing procedure to be employed: "SP" for the semi-parametric
#' permutation testing procedure with GPD approximation, "ASY" for the test based on asymptotic theory; if no method is specified, "SP" will be used by default.
#'@param permnum number of permutations used in the permutation testing
#' procedure (if method=”SP” is performed); default is 10000
#' 
#'@return A vector, see Schefzik et al. (2019) for details:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed by
#' quantile approximation
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#' computed by quantile approximation
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#' computed by decomposition approximation
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#' decomposition approximation
#' \item location: location term in the decomposition of the squared
#' 2-Wasserstein distance between the two samples
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#' \item rho: correlation coefficient in the quantile-quantile plot
#' \item pval: The p-value of the semi-parametric or the asymptotic theory-based test, depending on the specified method
#' \item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' (otherwise NA). This output is only returned when performing the
#' semi-parametric test (method="SP")!
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit, i.e. p-value of Anderson-Darling test
#' \eqn{\geq 0.05} (otherwise NA). This output is only returned when
#' performing the semi-parametric test (method="SP")!
#' \item perc.loc: fraction (in \%) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.size: fraction (in \%) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.shape: fraction (in \%) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#' distance computed by the quantile approximation and the squared
#' 2-Wasserstein distance computed by the decomposition approximation
#' }
#'
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
#'
#'@examples
#'# generate two input distributions
#'x<-rnorm(500)
#'y<-rnorm(500,4,1.5)
#'wasserstein.test(x,y,method="ASY")
#'# Run with default options: method="SP", permnum=10000
#'wasserstein.test(x,y)
#'# Run with a seed for the semi-parametric test ("SP")
#'set.seed(42)
#'wasserstein.test(x,y, method="SP")
#'
#'@export
#'
wasserstein.test <- function(x, y, method=c("SP", "ASY"), permnum=10000){
    method <- match.arg(method)
    switch(method,
           "SP"=.wassersteinTestSp(x, y, permnum),
           "ASY"=.wassersteinTestAsy(x, y))
}
