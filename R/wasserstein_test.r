
#' relativeError
#'
#' Computes the relative error between two numericals (in the sense of
#' is.numeric) x and y.
#' If x and y are vectors, it is assumed that length(x) == length(y).
#' 
#' The relative error $\delta$ is defined as: 
#' $$ \delta = \bigg| 1 -  \frac{x}{y} \bigg| $$
#' 
#' @param x numerical (in the sense of is.numeric)
#' @param y numerical (in the sense of is.numeric)
#' 
#' @return The relative error between x and y
#' 
#' @example 
#' pi1 <- 3.1415
#' pi2 <- 3.141
#' relativeError(x,y)
#' 
relativeError <- function(x,y) {
    rel.error <- abs(1 - (x / y))  
    return(rel.error)
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
#'\itemize{
#'\item d.wass: 2-Wasserstein distance between the two samples computed by
#' quantile approximation
#'\item d.wass^2: squared 2-Wasserstein distance between the two samples
#' computed by quantile approximation
#'\item d.comp^2: squared 2-Wasserstein distance between the two samples
#' computed by decomposition approximation
#'\item d.comp: 2-Wasserstein distance between the two samples computed by
#' decomposition approximation
#'\item location: location term in the decomposition of the squared
#' 2-Wasserstein distance between the two samples
#'\item size: size term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#'\item shape: shape term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#'\item rho: correlation coefficient in the quantile-quantile plot
#'\item pval: p-value of the semi-parametric 2-Wasserstein distance-based
#' test
#'\item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' (otherwise NA)
#'\item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit (i.e. p-value of Anderson-Darling test
#' greater or eqaul to 0.05) (otherwise NA)
#'\item perc.loc: fraction (in %) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item perc.size: fraction (in %) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item perc.shape: fraction (in %) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item decomp.error: relative error between the squared 2-Wasserstein
#' distance computed by the quantile approximation and the squared
#' 2-Wasserstein distance computed by the decomposition approximation
#'}
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
wasserstein.test.sp<-function(x,y,permnum){

    if (length(x) !=0 & length(y) != 0){

        value <- wasserstein_metric(x, y, p=2)
        value.sq <- value**2

        # computation of an approximative p-value (permutation procedure)
        z <- c(x,y)
        nsample <- length(z)
        bsn <- permnum

        shuffle <- permutations(z, num_permutations=bsn)
        wass.val <- apply(shuffle, 2, 
                            function (k) {
                                wasserstein_metric(
                                    k[seq_len(length(x))],
                                    k[seq((length(x)+1):length(z))],
                                    p=2)
                        })
        wass.val <- wass.val**2

        # list of possible exceedance thresholds (decreasing)
        poss.exc.num <- rev(seq(from=10, to=250, by=10))

        # order the values of the permutation-based test statistics
        wass.val.ordered <- sort(wass.val,decreasing=TRUE)

        # algorithm
        num.extr <- sum(wass.val >= value.sq)
        pvalue.ecdf <- num.extr/bsn
        pvalue.ecdf.pseudo <- (1 + num.extr) / (bsn + 1)

        if (num.extr >= 10){
            pvalue.wass <- c(pvalue.ecdf, NA, NA)  
        } else {

            procedure <- function(x) {
            
                r <- 1
                
                repeat {

                    # set threshold for exceedance according to paper
                    N.exc <- poss.exc.num[r]

                    # compute set of N.exc exceedances
                    exceedances <- wass.val.ordered[seq_len(N.exc)]

                    # check whether the N.exc largest permutation values follow 
                    # a GPD using an Anderson-Darling test
                    gpd.ad.check <- gpdAd(exceedances)
                    ad.pval <- gpd.ad.check$p.value

                    r <- r + 1

                    if (ad.pval > 0.05) {break}
                }

                # calculate exceedance threshold for so-obtained N.exc
                t.exc <- (wass.val.ordered[N.exc] 
                            + wass.val.ordered[N.exc+1]) / 2

                # fit GPD distribution to the exceedances using maximum 
                # likelihood estimation
                gpd.fit <- gpdFit(data=wass.val.ordered,
                                    threshold=t.exc,
                                    method="mle")

                # extract fitted parameters
                fit.scale <- as.numeric(gpd.fit$par.ests[1])
                fit.shape <- as.numeric(gpd.fit$par.ests[2])
                
                # compute GPD p-value (see paper)
                pvalue.gpd <- (N.exc / bsn) * (1 - pgpd(q=value.sq-t.exc,
                                                        loc=0,
                                                        scale=fit.scale,
                                                        shape=fit.shape))

                pvalue.gpd <- as.numeric(pvalue.gpd)
                
                pvalue.wass <- c(pvalue.gpd, ad.pval, N.exc)

                return(pvalue.wass)

            }

            tr <- try(procedure(25), silent=TRUE)

            if (is(tr,"try-error")) {
                pvalue.wass<-c(pvalue.ecdf.pseudo,NA,NA)
            } else {
                pvalue.wass<-tr
            }
        }


        # decomposition of wasserstein distance
        mu.x <- mean(x)
        mu.y <- mean(y)

        sigma.x <- sd(x)
        sigma.y <- sd(y)

        pr <- ((seq_len(1000)) - 0.5) / 1000

        quant.x <- quantile(x, probs=pr, type=1)
        quant.y <- quantile(y, probs=pr, type=1)

        if (sd(quant.x) !=0 & sd(quant.y) != 0){
            rho.xy <- cor(quant.x, quant.y, method="pearson")
        } else {
            rho.xy <- 0  
        }

        wass.comp <- squared_wass_decomp(x, y)
        location <- wass.comp$location
        size <- wass.comp$size
        shape <- wass.comp$shape
        
        d.comp.sq <- wass.comp$distance
        d.comp <- sqrt(d.comp.sq)
        
        perc.loc <- round(((location / d.comp.sq) * 100), 2)
        perc.size <- round(((size / d.comp.sq) * 100), 2)
        perc.shape <- round(((shape / d.comp.sq)*100), 2)
        
        decomp.error <- relativeError(value.sq, d.comp.sq)

    } else {
        value <- NA
        value.sq <- NA
        wass.comp.sq <- NA
        wass.comp <- NA
        location <- NA
        size <- NA
        shape <- NA
        rho.xy <- NA
        pvalue.wass <- c(NA, NA, NA)
        perc.loc <- NA
        perc.size <- NA
        perc.shape <- NA
        decomp.error <- NA
    }

    # create output
    output <- c(value, value.sq, d.comp.sq, d.comp, location, size, 
                shape, rho.xy, pvalue.wass, perc.loc, perc.size, perc.shape,
                decomp.error)
    names(output)<-c("d.wass", "d.wass^2", "d.comp^2", "d.comp",
                    "location", "size","shape", "rho","pval", "p.ad.gpd",
                    "N.exc", "perc.loc", "perc.size", "perc.shape",
                    "decomp.error")
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
#'\itemize{
#'\item d.wass: 2-Wasserstein distance between the two samples computed
#' by quantile approximation
#'\item d.wass^2 squared 2-Wasserstein distance between the two samples
#' computed by quantile approximation
#'\item d.comp^2: squared 2-Wasserstein distance between the two samples
#' computed by decomposition approximation
#'\item d.comp: 2-Wasserstein distance between the two samples computed by
#' decomposition approximation
#'\item location: location term in the decomposition of the squared
#' 2-Wasserstein distance between the two samples
#'\item size: size term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#'\item shape: shape term in the decomposition of the squared 2-Wasserstein
#' distance between the two samples
#'\item rho: correlation coefficient in the quantile-quantile plot
#'\item pval: p-value of the 2-Wasserstein distance-based test using
#' asymptotic theory
#'\item perc.loc: fraction (in %) of the location part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item perc.size: fraction (in %) of the size part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item perc.shape: fraction (in %) of the shape part with respect to the
#' overall squared 2-Wasserstein distance obtained by the decomposition
#' approximation
#'\item decomp.error: relative error between the squared 2-Wasserstein
#' distance computed by the quantile approximation and the squared
#' 2-Wasserstein distance computed by the decomposition approximation
#'}
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

        # Load the reference distributions from cache
        if (is.null(brownianbridge.empcdf)) {
            brownianbridge.empcdf.path <- cache.getOrDownload(
                url = brownianbridge.empcdf.url,
                rname = "empcdf.ref")
            load(brownianbridge.empcdf.path)
            # make the variable empcdf.ref in local scope globally available
            brownianbridge.empcdf <- empcdf.ref
        }

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

        # decomposition of wasserstein distance
        mu.x <- mean(x)
        mu.y <- mean(y)

        sigma.x <- sd(x)
        sigma.y <- sd(y)

        pr <- (seq_len(1000) - 0.5) / 1000

        quant.x <- quantile(x, probs=pr, type=1)
        quant.y <- quantile(y, probs=pr, type=1)

        if (sd(quant.x) != 0 & sd(quant.y) != 0) {
            rho.xy <- cor(quant.x, quant.y, method="pearson")
        } else {
            rho.xy <- 0  
        }

        wass.comp <- squared_wass_decomp(x, y)
        location <- wass.comp$location
        size <- wass.comp$size
        shape <- wass.comp$shape

        d.comp.sq <- wass.comp$distance
        d.comp <- sqrt(d.comp.sq)

        perc.loc <- round(((location / d.comp.sq) * 100), 2)
        perc.size <- round(((size / d.comp.sq) * 100), 2)
        perc.shape <- round(((shape / d.comp.sq)*100), 2)

        decomp.error <- relativeError(value.sq, d.comp.sq)


    } else {
        value <- NA
        value.sq <- NA
        d.comp.sq <- NA
        wass.comp <- NA
        location <- NA
        size <- NA
        shape <- NA
        rho.xy <- NA
        pvalue.wass <- NA
        perc.loc <- NA
        perc.size <- NA
        perc.shape <- NA
        decomp.error <- NA
    }


    ##create output
    output <- c(value, value.sq, d.comp.sq, wass.comp, location, size,
                shape, rho.xy, pvalue.wass, perc.loc, perc.size, perc.shape,
                decomp.error)
    names(output) <- c("d.wass", "d.wass^2", "d.comp^2", "d.comp",
                        "location", "size", "shape", "rho", "pval", "perc.loc",
                        "perc.size", "perc.shape", "decomp.error")

    return(output)
}




# Wasserstein test wrapper function, exposed for this package

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
