
#'.wassPermProcedure
#'
#' Permutation Procedure that calculates the squared 2-Wasserstein distance for
#' random shuffles of two input distributions and returns them as a vector.
#'
#' @param x numeric vector representing distribution A
#' @param y numeric vector representing distribution B
#' @param permnum integer interpreted as the number of permutations to be
#'  performed
#' 
#' @return vector with squared 2-Wasserstein distances computed on random
#'  shuffles of the two input vectors
#'  
.wassPermProcedure <- function(x, y, permnum) {
    z <- c(x,y)
    shuffle <- permutations(z, num_permutations=permnum)
    wass.values <- apply(   shuffle, 2, 
                            function (k) {
                                wasserstein_metric(
                                    k[seq_len(length(x))],
                                    k[seq((length(x)+1):length(z))],
                                    p=2) **2
                            })
    return(wass.values)
}


#'wasserstein.test.sp
#' 
#' Two-sample test to check for differences between two distributions
#' (conditions) using the 2-Wasserstein distance: Semi-parametric
#' implementation using a permutation test with a generalized Pareto
#' distribution (GPD) approximation to estimate small p-values accurately
#' 
#' This is the semi-parametric version of wasserstein.test, for the
#' asymptotic procedure see wasserstein.test.asy
#' 
#'@details Details concerning the permutation testing procedure with GPD
#' approximation to estimate small p-values accurately can be found in Schefzik
#' and Goncalves (2019).
#'
#'@param x univariate sample (vector) representing the distribution of
#' condition A
#'@param y univariate sample (vector) representing the distribution of
#' condition B
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'
#'@return A vector concerning the testing results, precisely (see Schefzik and
#' Goncalves (2019) for details)
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
#' (otherwise NA)
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit (i.e. p-value of Anderson-Darling test
#' greater or eqaul to 0.05) (otherwise NA)
#' \item perc.loc: fraction (in %) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.size: fraction (in %) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.shape: fraction (in %) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#' distance computed by the quantile approximation and the squared
#' 2-Wasserstein distance computed by the decomposition approximation
#' }
#'
#'@references Schefzik, R. and Goncalves, A. (2019).
#'
#'@examples
#'set.seed(32)
#'x<-rnorm(500)
#'y<-rnorm(500,2,1.5)
#'wasserstein.test.sp(x,y,10000)
#'
#'
#'@export
#'
wasserstein.test.sp<-function(x,y,permnum=10000){
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
        if (num.extr < 10) {
            tryCatch({
                    res <- .gdpFittedPValue(value.sq,
                                            wass.values.ordered,
                                            pvalue.ecdf)
                    pvalue.wass <<- res$pvalue.gdp
                    pvalue.gdpfit <<- res$ad.pval
                    N.exc <<- res$N.exc
                }, error=function(...) {
                    pvalue.wass <<- pvalue.ecdf.pseudo
                    pvalue.gdpfit <<- NA
                    N.exc <<- NA
                })
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
        decomp.error <- .relativeError(value.sq, d.comp.sq)

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


#'wasserstein.test.asy
#' 
#' Two-sample test to check for differences between two distributions
#' (conditions) using the 2-Wasserstein distance: Implementation using a test
#' based on asymptotic theory
#'
#' This is the asymptotic version of wasserstein.test, for the
#' semi-parametric procedure see wasserstein.test.sp
#'
#'@details Details concerning the testing procedure based on asymptotic theory
#' can be found in Schefzik and Goncalves (2019).
#'
#'@param x univariate sample (vector) representing the distribution of
#' condition A
#'@param y univariate sample (vector) representing the distribution of
#' condition B
#'
#'@return A vector concerning the testing results, precisely (see Schefzik and
#' Goncalves (2019) for details)
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed
#' by quantile approximation
#' \item d.wass^2 squared 2-Wasserstein distance between the two samples
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
#' \item perc.loc: fraction (in %) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.size: fraction (in %) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item perc.shape: fraction (in %) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#' distance computed by the quantile approximation and the squared
#' 2-Wasserstein distance computed by the decomposition approximation
#' }
#'
#'@references Schefzik, R. and Goncalves, A. (2019).
#'
#'@examples
#'set.seed(32)
#'x<-rnorm(500)
#'y<-rnorm(500,2,1.5)
#'wasserstein.test.asy(x,y)
#'
#'
#'@export
#'
wasserstein.test.asy <- function(x, y){

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
        pvalue.wass <- 1 - brownianbridge.empcdf(test.stat)

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


#'wasserstein.test
#'
#'Two-sample test to check for differences between two distributions
#'(conditions) using the 2-Wasserstein distance, either using the
#'semi-parametric permutation testing procedure with GPD approximation to
#'estimate small p-values accurately or the test based on asymptotic theory
#'
#'@name wasserstein.test
#'@details Details concerning the two testing procedures (i.e. the permutation
#' testing procedure with GPD approximation to estimate small p-values
#' accurately and the test based on asymptotic theory) can be found in
#' Schefzik and Goncalves (2019).
#'
#'@param x univariate sample (vector) representing the distribution of
#' condition A
#'@param y univariate sample (vector) representing the distribution of
#' condition B
#'@param permnum number of permutations used in the permutation testing
#' procedure (if method=”SP” is performed); default is 10000
#'@param method testing procedure to be employed: “SP” for the semi-parametric
#' permutation testing procedure with GPD approximation to estimate small
#' p-values accurately; “ASY” for the test based on asymptotic theory
#'
#'@return A vector concerning the testing results (see Schefzik and Goncalves
#' (2019) for details), see the documentations for wasserstein.test.sp and
#' wasserstein.test.asy, respectively, for detailed descriptions.
#'
#'@references Schefzik, R. and Goncalves, A. (2019).
#'
#'@examples
#'set.seed(32)
#'x<-rnorm(500)
#'y<-rnorm(500,2,1.5)
#'wasserstein.test(x,y,10000,method="SP")
#'wasserstein.test(x,y,method="asy")
#'
#'
#'@export
#'
wasserstein.test <- function(x, y, permnum=10000, method){
    if(toupper(method) == "SP"){
        RES <- wasserstein.test.sp(x, y, permnum)
    } else if (toupper(method) == "ASY") {
        RES <- wasserstein.test.asy(x, y)  
    } else {
        stop("Argument 'method' must be one of {SP, ASY} : ", method)
    }
    return(RES) 
}
