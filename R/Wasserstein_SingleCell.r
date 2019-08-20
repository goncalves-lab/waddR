
#'.fishersCombinedPval
#'
#'Aggregates two p-values into a combined p-value according to Fisher’s method 
#'
#'@details Aggregates two p-values into a combined p-value according to
#' Fisher’s method.
#'
#'@param x vector of the two p-values that are to be aggregated
#'
#'@return The combined p-value 
#'
#'@examples
#'x<-c(0.0004,2e-7)
#'.fishersCombinedPval(x)
#'
#'
#'@export
#'
.fishersCombinedPval<-function(x) { 
    if(sum(is.na(x)) == 0) {
        x <- ifelse(x==0,1e-100,x)
        p.comb <- pchisq(-2 * sum(log(x)), df=2*length(x), lower.tail=FALSE)
    } else if (sum(is.na(x)) == 1){
        p.comb <- x[!is.na(x)]
    } else {
        p.comb <- NA
    }
    return(p.comb)
}


#'.combinePVal
#'
#'For a given set of N pairs of p-values, aggregates each respective pair of 
#'p-values into a combined p-value according to Fisher’s method 
#'
#'@details For a given set of pairs of p-values, aggregates each respective 
#' pair of p-values into a combined p-value according to Fisher’s method.
#' Applies the .fishersCombinedPval function to a whole set of N pairs of
#' p-values. 
#'
#'@param r vector of length N of the p-values corresponding to the first test
#'@param s vector of length N of the p-values corresponding to the second test
#'
#'@return A vector of length N of the combined p-values
#'
#'@examples
#'r<-c(0.007,6e-23,0.77,0.23,2e-22)
#'s<-c(1e-6,3e-54,0.55,0.001,4e-5)
#'.combinePVal(r,s)
#'
#'
#'@export
#'
.combinePVal<-function(r,s){
    apply(cbind(r,s), 1, function(x)
        .fishersCombinedPval(x))
}


#'testZeroes
#'
#' Test for differential proportions of zero expression between two conditions
#' for a specified set of genes
#'
#'@details Test for differential proportions of zero expression between two
#' conditions that is not explained by the detection rate using a (Bayesian)
#' logistic regression model. Adapted from the scDD package (Korthauer et al.
#' 2016).
#'
#'@param x matrix of single-cell RNA-sequencing expression data with genes in
#' rows and samples (cells) in columns
#'@param y vector of condition labels
#'@param these vector of row numbers (i.e. gene numbers) employed to test for
#' differential proportions of zero expression. Default is 1:nrow(dat)
#'
#'@return A vector of (unadjusted) p-values 
#'
#'@references Korthauer et al. (2016).
#'
#'@examples
#'dat<-matrix(c(rnorm(100,42,1), rnorm(102,45,3)), nrow=1)
#'condition<-c(rep(1,100), rep(2,102))
#'testZeroes(dat,condition)
#'
#'
#'@export
#'
setGeneric("testZeroes",
    function(x, y, ...) standardGeneric("testZeroes"))


#'@export
setMethod("testZeroes",
    c(x="matrix", y="vector"),
    function(x, y, these=seq_len(nrow(x))) {
        detection <- colSums(x > 0) / nrow(x)
        onegene <- function(j, dat, detection, cond, these) {
            trows = dat[these[j], ]
            if (sum( y== 0) > 0) {
                M1 <- suppressWarnings(
                            bayesglm(   trows > 0 ~ detection + factor(cond),
                                        family=binomial(link="logit"),
                                        Warning=FALSE))
                return(summary(M1)$coefficients[3, 4])
            } else {
                return(NA)
            }
        }
        
        pval <- unlist(bplapply(seq_along(these), onegene, dat=x,
                                detection=detection, cond=y, these=these))
        return(pval)
    })


#'@export
setMethod("testZeroes",
    c(x="SingleCellExperiment", y="SingleCellExperiment"),
    function(x, y) {
        dat <- cbind(counts(x), counts(y))
        condition <- c(rep(1, dim(counts(x))[2]), rep(2, dim(counts(y))[2]))
        return(testZeroes(dat, condition))
    })


#'.testWass
#'
#'Two-sample test for single-cell RNA-sequencing data to check for differences
#'between two distributions (conditions) using the 2-Wasserstein distance:
#'Semi-parametric implementation using a permutation test with a generalized 
#'Pareto distribution (GPD) approximation to estimate small p-values accurately
#'
#'@details Details concerning the permutation testing procedures for
#' single-cell RNA-sequencing data can be found in Schefzik and 
#' Goncalves (2019).
#'
#'@param dat matrix of single-cell RNA-sequencing expression data with genes in
#' rows and samples (cells) in columns
#'@param condition vector of condition labels
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'@param inclZero logical; if TRUE, the one-stage method (i.e. semi-parametric
#' testing applied to all (zero and non-zero) expression values) is performed;
#' if FALSE, the two-stage method (i.e. semi-parametric testing applied to
#' non-zero expression values only, combined with a separate testing for
#' differential proportions of zero expression using logistic regression) is
#' performed. Default is TRUE
#'
#'@return A vector concerning the testing results, precisely (see Schefzik and
#' Goncalves (2019) for details) in case of inclZero=TRUE:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed
#'  by quantile approximation
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#'  computed by quantile approximation
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#'  computed by decomposition approximation
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#'  decomposition approximation
#' \item location: location term in the decomposition of the squared
#'  2-Wasserstein distance between the two samples
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples
#' \item rho: correlation coefficient in the quantile-quantile plot
#' \item pval: p-value of the semi-parametric 2-Wasserstein distance-based
#'  test
#' \item p.ad.gpd in case the GPD fitting is performed: p-value of the
#'  Anderson-Darling test to check whether the GPD actually fits the data well
#'  (otherwise NA)
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#'  (starting with 250 and iteratively decreased by 10 if necessary) that are
#'  required to obtain a good GPD fit (i.e. p-value of Anderson-Darling test
#'  greater or eqaul to 0.05)(otherwise NA)
#' \item perc.loc: fraction (in %) of the location part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.size: fraction (in %) of the size part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.shape: fraction (in %) of the shape part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#'  distance computed by the quantile approximation and the squared
#'  2-Wasserstein distance computed by the decomposition approximation
#' \item pval.adj: adjusted p-value of the semi-parametric 2-Wasserstein
#'  distance-based test according to the method of Benjamini-Hochberg
#' }
#' In case of inclZero=FALSE:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed
#'  by quantile approximation
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#'  computed by quantile approximation
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#'  computed by decomposition approximation
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#'  decomposition approximation
#' \item location: location term in the decomposition of the squared
#'  2-Wasserstein distance between the two samples
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples
#' \item rho: correlation coefficient in the quantile-quantile plot
#' \item p.nonzero: p-value of the semi-parametric 2-Wasserstein distance-based
#'  test (based on non-zero expression only)
#' \item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#'  Anderson-Darling test to check whether the GPD actually fits the data well
#'  (otherwise NA)
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#'  (starting with 250 and iteratively decreased by 10 if necessary) that are
#'  required to obtain a good GPD fit (i.e. p-value of Anderson-Darling test
#'  greater or eqaul to 0.05)(otherwise NA)
#' \item perc.loc: fraction (in %) of the location part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.size: fraction (in %) of the size part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.shape: fraction (in %) of the shape part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#'  distance computed by the quantile approximation and the squared 
#'  2-Wasserstein distance computed by the decomposition approximation
#' \item p.zero: p-value of the test for differential proportions of zero
#'  expression (logistic regression model)
#' \item p.combined: combined p-value of p.nonzero and p.zero obtained by
#'  Fisher’s method
#' \item p.adj.nonzero: adjusted p-value of the semi-parametric 2-Wasserstein
#'  distance-based test (based on non-zero expression only) according to the
#'  method of Benjamini-Hochberg
#' \item p.adj.zero: adjusted p-value of the test for differential proportions
#'  of zero expression (logistic regression model) according to the method of
#'  Benjamini-Hochberg
#' \item p.adj.combined: adjusted combined p-value of p.nonzero and p.zero
#'  obtained by Fisher’s method according to the method of Benjamini-Hochberg
#' }
#' If 
#'
#'@references Schefzik and Goncalves (2019).
#'
#'@examples
#'dat <- matrix(c(rnorm(100, 42, 1), rnorm(102, 45, 3)), nrow=1)
#'condition <- c(rep(1, 100), rep(2, 102))
#'testWass(dat, condition, 10000, inclZero=TRUE)
#'testWass(dat, condition, 10000, inclZero=FALSE)
#'
#'
#'        
.testWass <- function(dat, condition, permnum, inclZero=TRUE){

    if (!inclZero){

        onegene <- function(x, dat, condition){
            x1 <- dat[x,][condition==unique(condition)[1]]
            x2 <- dat[x,][condition==unique(condition)[2]]

            x1 <- (x1[x1>0])
            x2 <- (x2[x2>0])

            suppressWarnings(wasserstein.test.sp(x1,x2,permnum))
        }

        wass.res <- bplapply(seq_len(nrow(dat)), onegene, 
                            condition=condition, dat=dat)

        wass.res1 <- do.call(rbind, wass.res)

        wass.pval.adj <- p.adjust(wass.res1[,9], method="BH")

        pval.zero <- testZeroes(dat,condition)
        pval.adj.zero <- p.adjust(pval.zero,method="BH")  

        pval.combined <- .combinePVal(wass.res1[,9],pval.zero)
        pval.adj.combined <- p.adjust(pval.combined,method="BH")

        RES <- cbind(wass.res1,pval.zero,pval.combined,wass.pval.adj,
                    pval.adj.zero,pval.adj.combined)
        row.names(RES) <- rownames(dat)
        colnames(RES) <- c( colnames(wass.res1), "p.zero", "p.combined",
                            "p.adj.nonzero","p.adj.zero","p.adj.combined")
        return(RES)
    } else {

        onegene <- function(x, dat, condition){
            x1 <- dat[x,][condition==unique(condition)[1]]
            x2 <- dat[x,][condition==unique(condition)[2]]
            
            suppressWarnings(wasserstein.test.sp(x1,x2,permnum))
        }
        
        wass.res <- bplapply(seq_len(nrow(dat)), onegene, 
                        condition=condition, dat=dat)

        wass.res1 <- do.call(rbind,wass.res)

        wass.pval.adj <- p.adjust(wass.res1[,9], method="BH")
        
        RES <- cbind(wass.res1, wass.pval.adj)
        row.names(RES) <- rownames(dat)
        return(RES)
    }
}


#'wasserstein.sc
#'
#' Two-sample test for single-cell RNA-sequencing data to check for differences
#' between two distributions (conditions) using the 2-Wasserstein distance:
#' Semi-parametric implementation using a permutation test with a generalized
#' Pareto distribution (GPD) approximation to estimate small p-values
#' accurately
#'
#'@details Details concerning the permutation testing procedures for
#' single-cell RNA-sequencing data can be found in Schefzik and Goncalves
#' (2019). Corresponds to the function .testWass when identifying the argument
#' inclZero=TRUE in .testWass with the argument method=”OS” and the argument
#' inclZero=FALSE in .testWass with the argument method=”TS”.
#'
#'@param x matrix of single-cell RNA-sequencing expression data with genes in
#' rows and samples (cells) in columns
#'@param y vector of condition labels
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'@param method method employed in the testing procedure: “OS” for the
#' one-stage method (i.e. semi-parametric testing applied to all (zero and
#' non-zero) expression values); “TS” for the two-stage method (i.e.
#' semi-parametric testing applied to non-zero expression values only, combined
#' with a separate testing for differential proportions of zero expression
#' using logistic regression)
#'
#'@return See the corresponding values in the description of the function
#' .testWass, where the argument inclZero=TRUE in .testWass has to be
#' identified with the argument method=”OS”, and the argument inclZero=FALSE in
#' .testWass with the argument method=”TS”.
#'
#'@references Schefzik and Goncalves (2019).
#'
#'@examples
#'
#'# some data in two conditions
#'cond1 <- matrix(rnorm(100, 42, 1), nrow=1)
#'cond2 <- matrix(rnorm(100, 45, 3), nrow=1)
#'
#'# call wasserstein.sc with a matrix
#'# and a vector denoting conditions
#'dat <- cbind(cond1, cond2)
#'condition <- c(rep(1, 100), rep(2, 100))
#'wasserstein.sc(dat, condition, 100, method="TS")
#'
#'# call wasserstein.sc with two SingleCellExperiment objects
#'sce1 <- SingleCellExperiment(
#'            assays=list(counts=cond1, logcounts=log10(cond1)))
#'sce2 <- SingleCellExperiment(
#'            assays=list(counts=cond2, logcounts=log10(cond2)))
#'wasserstein.sc(sce1, sce2, 100, method="TS")
#'
#'@export
setGeneric("wasserstein.sc",
    function(x, y, ...) standardGeneric("wasserstein.sc"))


#'@export
setMethod("wasserstein.sc", 
    c(x="matrix", y="vector"),
    function(x, y, permnum, method) {
        method <- toupper(method)
        stopifnot(method %in% c("OS", "TS"))
        stopifnot(dim(x)[2] == length(y))

        if(method == "OS")
            .testWass(x, y, permnum, inclZero=TRUE)
        if(method == "TS")
            .testWass(x, y, permnum, inclZero=FALSE)
    })


#'@export
setMethod("wasserstein.sc",
    c(x="SingleCellExperiment", y="SingleCellExperiment"),
    function(x, y, permnum, method) {
        method <- toupper(method)
        stopifnot(method %in% c("OS", "TS"))
        stopifnot(dim(counts(x))[1] == dim(counts(y))[1])

        dat <- cbind(counts(x), counts(y))
        condition <- c(rep(1, dim(counts(x))[2]), rep(2, dim(counts(y))[2]))
        if(method == "OS")
            .testWass(dat, condition, permnum, inclZero=TRUE)
        if(method == "TS")
            .testWass(dat, condition, permnum, inclZero=FALSE)
    })

