library(BiocParallel)
library(eva)
sourceCpp("src/wasserstein_test.cpp")


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









########################################
########semi-parametric Wasserstein approach to check differential expression for non-zero expression

wasserstein.test.sp<-function(x,y,seedex,permnum){
  
  set.seed(seedex)  
  
  
  #######Wasserstein part
  
  if (length(x)!=0&length(y)!=0){
    
    value<-wasserstein_metric(x,y,p=2)
    value.sq<-value^2
    
    
    ##### computation of an approximative p-value (permutation procedure)
    z<-c(x,y)
    nsample<-length(z)
    bsn<-permnum
    
    shuffle <- permutations(z, n = bsn)
    wass.val <- apply(shuffle, 2, function (k) {wasserstein_metric(k[1:length(x)], k[(length(x)+1):5], p=2)})
    wass.val<-wass.val^2
    
    
    
    ###list of possible exceedance thresholds (decreasing)
    poss.exc.num<-rev(seq(from=10,to=250,by=10))
    
    
    
    
    ##order the values of the permutation-based test statistics
    wass.val.ordered<-sort(wass.val,decreasing=TRUE)
    
    
    ###algorithm
    
    num.extr<-sum(wass.val>=value.sq)
    pvalue.ecdf<-num.extr/bsn
    pvalue.ecdf.pseudo<-(1+num.extr)/(bsn+1)
    
    
    if (num.extr>=10){
      pvalue.wass<-c(pvalue.ecdf,NA,NA)  
    } else {
      
      procedure<-function(x){
        
        r<-1
        
        repeat {
          
          
          ##set threshold for exceedance according to paper
          N.exc<-poss.exc.num[r]
          
          ##compute set of N.exc exceedances
          exceedances<-wass.val.ordered[1:N.exc]
          
          ##check whether the N.exc largest permutation values follow a GPD using an Anderson-Darling test
          gpd.ad.check<-gpdAd(exceedances)
          ad.pval<-gpd.ad.check$p.value
          
          r<-r+1
          
          if(ad.pval>0.05){break}
          
        }
        
        
        ###calculate exceedance threshold for so-obtained N.exc
        t.exc<-(wass.val.ordered[N.exc]+wass.val.ordered[N.exc+1])/2
        
        ###fit GPD distribution to the exceedances using maximum likelihood estimation
        #gpd.fit<-gpdFit(x=wass.val.ordered,u=t.exc,type="mle")
        gpd.fit<-gpdFit(data=wass.val.ordered,threshold=t.exc,method="mle")
        ##extract fitted parameters
        fit.scale<-as.numeric(gpd.fit$par.ests[1])
        fit.shape<-as.numeric(gpd.fit$par.ests[2])
        
        ##check goodness-of-fit of the fitted GPD distribution via Anderson-Darling test 
        #gof.test<-ad.test(exceedances, "pgpd", xi=fit.shape,beta=fit.scale)
        #pvalue.gof<-gof.test$p.value
        #goodfit<-pvalue.gof>0.05
        
        ###compute GPD p-value (see paper)
        pvalue.gpd<-(N.exc/bsn)*(1-pgpd(q=value.sq-t.exc,loc=0,scale=fit.scale,shape=fit.shape))
        pvalue.gpd<-as.numeric(pvalue.gpd)
        
        pvalue.wass<-c(pvalue.gpd,ad.pval,N.exc)
        
        return(pvalue.wass)
        
      }
      
      tr<-try(procedure(25),silent=TRUE)
      
      if(is(tr,"try-error")) {pvalue.wass<-c(pvalue.ecdf.pseudo,NA,NA)}
      else {pvalue.wass<-tr}
      
    }
    
    
    
    ####decomposition of wasserstein distance
    mu.x<-mean(x)
    mu.y<-mean(y)
    
    sigma.x<-sd(x)
    sigma.y<-sd(y)
    
    
    pr<-((1:1000)-0.5)/1000
    
    quant.x<-quantile(x,probs=pr,type=1)
    quant.y<-quantile(y,probs=pr,type=1)
    
    if(sd(quant.x)!=0&sd(quant.y)!=0){
      rho.xy<-cor(quant.x,quant.y,method="pearson")
    } else{
      rho.xy<-0  
    }
    
    
    
    location<-(mu.x-mu.y)^2
    size<-(sigma.x-sigma.y)^2
    shape<-2*sigma.x*sigma.y*(1-rho.xy)
    
    wass.comp.sq<-location+size+shape
    wass.comp<-sqrt(wass.comp.sq)
    
    
    perc.loc<-round(((location/wass.comp.sq)*100),2)
    perc.size<-round(((size/wass.comp.sq)*100),2)
    perc.shape<-round(((shape/wass.comp.sq)*100),2)
    
    
    decomp.error<-abs(value.sq-wass.comp.sq)
    
    
    
  } else {
    value<-NA
    value.sq<-NA
    wass.comp.sq<-NA
    wass.comp<-NA
    location<-NA
    size<-NA
    shape<-NA
    rho.xy<-NA
    pvalue.wass<-c(NA,NA,NA)
    perc.loc<-NA
    perc.size<-NA
    perc.shape<-NA
    decomp.error<-NA
  }
  
  
  ##create output
  output<-c(value,value.sq,wass.comp.sq,wass.comp,location,size,shape,rho.xy,pvalue.wass,perc.loc,perc.size,perc.shape,decomp.error)
  return(output)
  
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



################################
#######overall testing procedure for single-cell data based on wasserstein distance:
##either one-stage (OS) approach : Wasserstein for zero and non-zero expression data
##or two-stage (TS) approach: Wasserstein for non-zero expression combined with test for differential proportions of zeroes

wasserstein.sc<-function(dat,condition,seedex,permnum,method){
  if(method=="OS")
    RES<-testWass(dat, condition,seedex,permnum, inclZero=TRUE)
  if(method=="TS")  
    RES<-testWass(dat, condition,seedex,permnum, inclZero=FALSE)  
  return(RES) 
}

