
#' Test for differential proportions of zero gene expression
#'
#' Test for differential proportions of zero expression between two conditions
#' for a specified set of genes; adapted from the R/Bioconductor package \code{scDD} by Korthauer et al. (2106)
#'
#' Test for differential proportions of zero gene expression between two
#' conditions using a logistic regression model accounting for the cellular detection rate. Adapted from the R/Bioconductor package \code{scDD} by Korthauer et al.(2016).
#'
#' @param x matrix of single-cell RNA-sequencing expression data with genes in
#'   rows and cells (samples) in columns
#' @param y vector of condition labels
#' @param these vector of row numbers (i.e. gene numbers) employed to test for
#'   differential proportions of zero expression; default is seq_len(nrow(dat))
#'
#' @return A vector of (unadjusted) p-values
#'
#' @references K.D. Korthauer, L.-F. Chu, M. A. Newton, Y. Li, J. Thomson, R. Stewart, and C. Kendziorski (2016). A statistical approach for identifying differential distributions in single-cell RNA-seq experiments. Genome Biology, 17:222.
#'
#' @examples
#' x1 <- c(rnorm(100,42,1), rnorm(102,45,3))
#' x2 <- c(rnorm(100,0,1), rnorm(102,0,2))
#' dat <- matrix(c(x1,x2), nrow=2, byrow=TRUE)
#' condition <- c(rep(1,100), rep(2,102))
#' # test over all rows
#' testZeroes(dat, condition)
#' # only consider the second row
#' testZeroes(dat, condition, these=c(2))
#'
#' @name testZeroes
#' @export
#' @docType methods
#' @rdname testZeroes-method
setGeneric("testZeroes",
    function(x, y, these=seq_len(nrow(x))) standardGeneric("testZeroes"))


#'@rdname testZeroes-method
#'@aliases testZeroes,matrix,vector,ANY-method
setMethod("testZeroes",
    c(x="matrix", y="vector"),
    function(x, y, these=seq_len(nrow(x))) {
        detection <- colSums(x > 0) / nrow(x)
        onegene <- function(j, dat, detection, cond, these) {
            trow = dat[these[j], ]
            if (sum(trow == 0) > 0) {
                M1 <- suppressWarnings(
                            bayesglm(   trow > 0 ~ detection + factor(cond),
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


#'@rdname testZeroes-method
#'@aliases testZeroes,SingleCellExperiment,SingleCellExperiment,vector-method
setMethod("testZeroes",
    c(x="SingleCellExperiment", y="SingleCellExperiment", these="vector"),
    function(x, y, these=seq_len(nrow(x))) {
        dat <- cbind(counts(x), counts(y))
        condition <- c(rep(1, dim(counts(x))[2]), rep(2, dim(counts(y))[2]))
        return(testZeroes(dat, condition, these))
    })


#'Check for differential distributions in single-cell RNA sequencing data via a semi-paramteric test using the 2-Wasserstein distance
#'           
#'Two-sample test for single-cell RNA-sequencing data to check for differences
#'between two distributions using the 2-Wasserstein distance:
#'Semi-parametric implementation using a permutation test with a generalized 
#'Pareto distribution (GPD) approximation to estimate small p-values accurately
#'
#'@details Details concerning the permutation testing procedures for
#' single-cell RNA-sequencing data can be found in Schefzik et al. (2019).
#'
#'@param dat matrix of single-cell RNA-sequencing expression data, with rowas corresponding to genes and columns corresponding to cells (samples)
#'@param condition vector of condition labels
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'@param inclZero logical; if TRUE, a one-stage method is performed, i.e. the semi-parametric
#' test based on the 2-Wasserstein distance is applied to all (zero and non-zero) expression values;
#' if FALSE, a two-stage method is performed, i.e. the semi-parametric test based on the 2-Wasserstein distance is applied to
#' non-zero expression values only, and a separate test for
#' differential proportions of zero expression using logistic regression is conducted; default is TRUE
#'@param seed Number to be used as a L'Ecuyer-CMRG seed, which itself
#' seeds the generation of an nextRNGStream() for each gene. Internally, when
#' this argument is given, a seed is specified by calling
#' `RNGkind("L'Ecuyer-CMRG")` followed by `set.seed(seed)`.
#' The `RNGkind` and `.Random.seed` will be reset on termination of this
#' function. Default is NULL, and no seed is set.
#'@return Matrix, where each row contains the testing results of the respective gene from \code{dat}.
#'  For the corresponding values of each row (gene), see the description of the function
#' \code{wasserstein.sc}, where the argument \code{inclZero=TRUE} in \code{.testWass} has to be
#' identified with the argument \code{method=”OS”}, and the argument \code{inclZero=FALSE} with the argument \code{method=”TS”}.
#' 
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
#'
.testWass <- function(dat, condition, permnum, inclZero=TRUE, seed=NULL){
    ngenes <- nrow(dat)
    seeds <- NULL
    
    if (!is.null(seed)) {
        if (exists(".Random.seed")) {
            oseed <- .Random.seed
        } else {
            oseed <- NULL
        }
        
        # set a reproducible L'Ecuyer seed based on the previous seed.
        # This way, the user doesn't have to worry about the RNGkind.
        okind <- RNGkind("L'Ecuyer-CMRG")[1]
        set.seed(seed)
        
        # Create seeds from RNGStream for each gene for Reproducibility
        seed <- .Random.seed
        seeds <- replicate(ngenes,
                           { seed <<- nextRNGStream(seed) },
                           simplify=FALSE)
        on.exit({
            if (is.null(oseed)) {
                RNGkind(okind)
                rm(.Random.seed)
            } else {
                # reset the RNGkind and seed, so the user environment is
                # unaffected
                RNGkind(okind)
                .Random.seed <<- oseed
            }
        })
    }
    
    # parallel worker 
    onegene <- function(x, seed=NULL){
        x1 <- dat[x,][condition==unique(condition)[1]]
        x2 <- dat[x,][condition==unique(condition)[2]]
        
        if (!inclZero) {
            x1 <- (x1[x1>0])
            x2 <- (x2[x2>0])
        }
        
        if (!is.null(seed)) {
            .Random.seed <<- seed
        }
        
        suppressWarnings(.wassersteinTestSp(x1, x2, permnum))
    }
    
    # run worker
    if (!is.null(seeds)) {
        wass.res <- t(simplify2array({
                            bpmapply(onegene, seq(ngenes), seeds)
                        }))
    } else {
        wass.res <- t(simplify2array({
                            bpmapply(onegene, seq(ngenes))
                        }))
    }

    #wass.res1 <- do.call(rbind, wass.res)
    wass.pval.adj <- p.adjust(wass.res[,9], method="BH")
    
    if (!inclZero){
        # zeroes were excluded => test them separately now
        pval.zero <- testZeroes(dat, condition)
        pval.adj.zero <- p.adjust(pval.zero, method="BH")
        pval.combined <- .combinePVal(wass.res[,9],pval.zero)
        pval.adj.combined <- p.adjust(pval.combined,method="BH")

        RES <- cbind(wass.res,pval.zero,pval.combined,wass.pval.adj,
                    pval.adj.zero,pval.adj.combined)
        row.names(RES) <- rownames(dat)
        colnames(RES) <- c( colnames(wass.res), "p.zero", "p.combined",
                            "p.adj.nonzero","p.adj.zero","p.adj.combined")
        return(RES)
    
    } else {

        RES <- cbind(wass.res, wass.pval.adj)
        row.names(RES) <- rownames(dat)
        colnames(RES) <- c( colnames(wass.res), "pval.adj")
        return(RES)
    }
}


#'Two-sample semi-parametric test for single-cell RNA-sequencing data to check for differences between two distributions using the 2-Wasserstein distance
#'
#' Two-sample test for single-cell RNA-sequencing data to check for differences
#' between two distributions using the 2-Wasserstein distance:
#' Semi-parametric implementation using a permutation test with a generalized
#' Pareto distribution (GPD) approximation to estimate small p-values
#' accurately
#'
#' Details concerning the permutation testing procedures for
#' single-cell RNA-sequencing data can be found in Schefzik et al.
#' (2019). Corresponds to the function \code{.testWass} when identifying the argument
#' \code{inclZero=TRUE} in \code{.testWass} with the argument \code{method=”OS”} and the argument
#' \code{inclZero=FALSE} with the argument \code{method=”TS”}.
#'
#'@param x matrix of single-cell RNA-sequencing expression data with genes in
#' rows and cells (samples) in columns
#'@param y vector of condition labels
#'@param method method employed in the testing procedure: if “OS” , a one-stage test is perofrmed, i.e. the semi-parametric test is applied to all (zero and
#' non-zero) expression values; if “TS”, a two-stage test is performed, i.e.
#' the semi-parametric test is applied to non-zero expression values only and combined
#' with a separate test for differential proportions of zero expression
#' using logistic regression
#'@param permnum number of permutations used in the permutation testing
#' procedure
#'@param seed number to be used to generate a L'Ecuyer-CMRG seed, which itself
#' seeds the generation of an nextRNGStream() for each gene to achieve
#' reproducibility; default is NULL, and no seed is set
#'@return Matrix, where each row contains the testing results of the respective gene from \code{dat}. The corresponding values of each row (gene) are as follows, see Schefzik et al. (2019) for details.     
#' In case of \code{inclZero=TRUE}:
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
#'  test when applied to all (zero and non-zero) respective gene expression values
#' \item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' (otherwise NA)
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit, i.e. p-value of Anderson-Darling test
#' \eqn{\geq 0.05} (otherwise NA)
#' \item perc.loc: fraction (in \%) of the location part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.size: fraction (in \%) of the size part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item perc.shape: fraction (in \%) of the shape part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation
#' \item decomp.error: relative error between the squared 2-Wasserstein
#'  distance computed by the quantile approximation and the squared
#'  2-Wasserstein distance computed by the decomposition approximation
#' \item pval.adj: adjusted p-value of the semi-parametric 2-Wasserstein
#'  distance-based test according to the method of Benjamini-Hochberg (i.e. adjusted p-value corresponding to pval)
#' }
#' In case of \code{inclZero=FALSE}:
#' \itemize{
#' \item d.wass: 2-Wasserstein distance between the two samples computed
#'  by quantile approximation (based on non-zero expression only)
#' \item d.wass^2: squared 2-Wasserstein distance between the two samples
#'  computed by quantile approximation (based on non-zero expression only)
#' \item d.comp^2: squared 2-Wasserstein distance between the two samples
#'  computed by decomposition approximation (based on non-zero expression only)
#' \item d.comp: 2-Wasserstein distance between the two samples computed by
#'  decomposition approximation (based on non-zero expression only)
#' \item location: location term in the decomposition of the squared
#'  2-Wasserstein distance between the two samples (based on non-zero expression only)
#' \item size: size term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples (based on non-zero expression only)
#' \item shape: shape term in the decomposition of the squared 2-Wasserstein
#'  distance between the two samples (based on non-zero expression only)
#' \item rho: correlation coefficient in the quantile-quantile plot (based on non-zero expression only)
#' \item p.nonzero: p-value of the semi-parametric 2-Wasserstein distance-based
#'  test when applied to non-zero respective gene expression values
#' \item p.ad.gpd: in case the GPD fitting is performed: p-value of the
#' Anderson-Darling test to check whether the GPD actually fits the data well
#' (otherwise NA)
#' \item N.exc: in case the GPD fitting is performed: number of exceedances
#' (starting with 250 and iteratively decreased by 10 if necessary) that are
#' required to obtain a good GPD fit, i.e. p-value of Anderson-Darling test
#' \eqn{\geq 0.05} (otherwise NA)
#' \item perc.loc: fraction (in \%) of the location part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation (based on non-zero expression only)
#' \item perc.size: fraction (in \%) of the size part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation (based on non-zero expression only)
#' \item perc.shape: fraction (in \%) of the shape part with respect to the
#'  overall squared 2-Wasserstein distance obtained by the decomposition
#'  approximation (based on non-zero expression only)
#' \item decomp.error: relative error between the squared 2-Wasserstein
#'  distance computed by the quantile approximation and the squared 
#'  2-Wasserstein distance computed by the decomposition approximation (based on non-zero expression only)
#' \item p.zero: p-value of the test for differential proportions of zero
#'  expression (logistic regression model)
#' \item p.combined: combined p-value of p.nonzero and p.zero obtained by
#'  Fisher’s method
#' \item p.adj.nonzero: adjusted p-value of the semi-parametric 2-Wasserstein
#'  distance-based test based on non-zero expression only according to the
#'  method of Benjamini-Hochberg (i.e. adjusted p-value corresponding to p.nonzero)
#' \item p.adj.zero: adjusted p-value of the test for differential proportions
#'  of zero expression (logistic regression model) according to the method of
#'  Benjamini-Hochberg (i.e. adjusted p-value corresponding to p.zero)
#' \item p.adj.combined: adjusted combined p-value of p.nonzero and p.zero
#'  obtained by Fisher’s method according to the method of Benjamini-Hochberg (i.e. adjusted p-value corresponding to p.combined)
#' }
#'
#'@references Schefzik, R., Flesch, J., and Goncalves, A. (2019). waddR: Using the 2-Wasserstein distance to identify differences between distributions in two-sample testing, with application to single-cell RNA-sequencing data.
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
#'wasserstein.sc(dat, condition, method="TS", 100)
#'
#'# call wasserstein.sc with two SingleCellExperiment objects
#'sce1 <- SingleCellExperiment::SingleCellExperiment(
#'            assays=list(counts=cond1, logcounts=log10(cond1)))
#'sce2 <- SingleCellExperiment::SingleCellExperiment(
#'            assays=list(counts=cond2, logcounts=log10(cond2)))
#'wasserstein.sc(sce1, sce2, method="TS", 100)
#'
#'# for reproducible p-values
#'wasserstein.sc(sce1, sce2, seed=123)
#'
#' @name wasserstein.sc
#' @export
#' @docType methods
#' @rdname wasserstein.sc-method
setGeneric("wasserstein.sc",
    function(x, y, method=c("TS", "OS"), permnum=10000, seed=NULL)
        standardGeneric("wasserstein.sc"))


#'@rdname wasserstein.sc-method
#'@aliases wasserstein.sc-method,matrix,vector,ANY,ANY,ANY-method
setMethod("wasserstein.sc", 
    c(x="matrix", y="vector"),
    function(x, y, method=c("TS", "OS"), permnum=10000, seed=NULL) {
        stopifnot(length(unique(y)) == 2)
        stopifnot(dim(x)[2] == length(y))
        
        method <- match.arg(method)
        switch(method,
               "TS"=.testWass(x, y, permnum, inclZero=FALSE, seed=seed),
               "OS"=.testWass(x, y, permnum, inclZero=TRUE, seed=seed))
    })


#'@rdname wasserstein.sc-method
#'@aliases
#'  wasserstein.sc,SingleCellExperiment,SingleCellExperiment,ANY,ANY,ANY-method
setMethod("wasserstein.sc",
    c(x="SingleCellExperiment", y="SingleCellExperiment"),
    function(x, y, method=c("TS", "OS"), permnum=10000, seed=NULL) {
        stopifnot(dim(counts(x))[1] == dim(counts(y))[1])
        
        
        dat <- cbind(counts(x), counts(y))
        condition <- c(rep(1, dim(counts(x))[2]), rep(2, dim(counts(y))[2]))
        method <- match.arg(method)
        switch(method,
               "TS"=.testWass(dat, condition, permnum, 
                              inclZero=FALSE, seed=seed),
               "OS"=.testWass(dat, condition, permnum, 
                              inclZero=TRUE, seed=seed))
    })

