library(BiocParallel)
library(eva)


######
###test for differential proportions of zeroes; adapted from scDD package

testZeroes <- function(dat, cond, these=1:nrow(dat)){
  detection <- colSums(dat>0)/nrow(dat)
  
  onegene <- function(j, dat, detection, cond, these){
    y=dat[these[j],]
    if (sum(y==0) > 0){
      M1 <- suppressWarnings(arm::bayesglm(y>0 ~ detection + factor(cond), 
                                           family=binomial(link="logit"),
                                           Warning=FALSE))
      return(summary(M1)$coefficients[3,4])
    }else{
      return(NA)
    }
  }
  
  pval <- unlist(bplapply(seq_along(these), onegene, dat=dat, 
                          detection=detection, cond=cond, these=these))
  
  return(pval)
}


###########
####compute Fisher's combined p-value (zero&non-zero)

fishersCombinedPval<-function(x){ 
  if(sum(is.na(x)) == 0){
    ifelse(x==0,1e-100,x)
    p.comb<-pchisq(-2 * sum(log(x)), df=2*length(x),lower.tail=FALSE)
  }else if(sum(is.na(x)) == 1){
    p.comb<-x[!is.na(x)]
  }else{
    p.comb<-NA
  }
  return(p.comb)
}


### compute combined p-value
CombinePVal<-function(r,s){
  apply(cbind(r,s), 1, function(x)
    fishersCombinedPval(x))
}


#######################
#####apply semi-parametric Wasserstein test; either include zero expression values (one-stage approach) or not (two-stage approach with separate test for differential proportions of zeroes)
#@param inclZero  : boolean representing whether or not to check the distribution of
#                   the expression value 0 (zero) independently
testWass<-function(dat, condition,seedex,permnum, inclZero=TRUE){
  
  if (!inclZero){
    
    onegene <- function(x, dat, condition){
      x1 <- dat[x,][condition==unique(condition)[1]]
      x2 <- dat[x,][condition==unique(condition)[2]]
      
      x1 <- (x1[x1>0])
      x2 <- (x2[x2>0])
      
      suppressWarnings(wasserstein.test.sp(x1,x2,seedex,permnum))
    }
    
    wass.res<- bplapply(seq_len(nrow(dat)), onegene, 
                        condition=condition, dat=dat)
    
    wass.res1<-do.call(rbind,wass.res)
    
    wass.pval.adj <- p.adjust(wass.res1[,9], method="BH")
    
    pval.zero<-testZeroes(dat,condition)
    pval.adj.zero<-p.adjust(pval.zero,method="BH")  
    
    pval.combined<-CombinePVal(wass.res1[,9],pval.zero)
    pval.adj.combined<-p.adjust(pval.combined,method="BH")
    
    RES<-cbind(wass.res1,pval.zero,pval.combined,wass.pval.adj,pval.adj.zero,pval.adj.combined)
    row.names(RES)<-rownames(dat)
    colnames(RES)<-c("d.transport","d.transport^2","d.comp^2","d.comp","location","size","shape","rho","p.nonzero","p.ad.gpd","N.exc","perc.loc","perc.size","perc.shape","decomp.error","p.zero","p.combined","p.adj.nonzero","p.adj.zero","p.adj.combined")
    
  }
  
  if(inclZero){
    
    onegene <- function(x, dat, condition){
      x1 <- dat[x,][condition==unique(condition)[1]]
      x2 <- dat[x,][condition==unique(condition)[2]]
      
      suppressWarnings(wasserstein.test.sp(x1,x2,seedex,permnum))
    }
    
    wass.res<- bplapply(seq_len(nrow(dat)), onegene, 
                        condition=condition, dat=dat)
    
    wass.res1<-do.call(rbind,wass.res)
    
    wass.pval.adj <- p.adjust(wass.res1[,9], method="BH")
    
    RES<-cbind(wass.res1,wass.pval.adj)
    row.names(RES)<-rownames(dat)
    colnames(RES)<-c("d.transport","d.transport^2","d.comp^2","d.comp","location","size","shape","rho","pval","p.ad.gpd","N.exc","perc.loc","perc.size","perc.shape","decomp.error","pval.adj")
    
    
  }
  
  return(RES)
  
}



#' Testing procedure for single-cell data based on the wasserstein distance
#'
#' wasserstein.sc wraps around the testWass function to provide a interface for calling
#' wasserstein distance based test of differential expression in sigle cell RNA data ...
#'
#' @param dat single cell expression values for two conditions
#' @param condition vector of the same length as dat, describing which \cr
#'  condition every datapoint is associated to
#' @param seedex seed for generating reproducable distributions
#' @param permnum number of permutations that should be produced in calculating the p-values
#' @param method one of "OS" or "TS", describing wether a one-stage ("OS") or
#'  two-stage ("TS") test should be performed. \cr
#'  TS - the distribution of the zero values is
#'  independently tested for differential proportions of zeroes \cr
#'  OS - all datapoints are tested with the wasserstein distance derived model
#'
#' @return Test statistics 
#'
#' @examples 
#' dat <- c(100, 120, 100, 90, 0, 0, 90, 80, 10, 50))
#' condition <- c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2)
#' wasserstein.sc(dat, condition, 24, 10000, "OS")
#'
#'
#' @export
wasserstein.sc<-function(dat,condition,seedex,permnum,method){
  if(method=="OS")
    RES<-testWass(dat, condition,seedex,permnum, inclZero=TRUE)
  if(method=="TS")  
    RES<-testWass(dat, condition,seedex,permnum, inclZero=FALSE)  
  return(RES) 
}

