
library(RLRsim)
library(nlme)
library(CompQuadForm)
library(MASS)

ReLRT=function(y,X,G,a1,a2,nsim,nmlp,method="REML")
{
  y  = y            ##the response
  X  = as.matrix(X) ##interpret + covariates
  G  = as.matrix(G) ##genotypes
  a1 = a1           ##a1 and a2 are the parameters for the  
  a2 = a2           ##probability density function of beta distribution
  nsim = nsim       ## the number of sampling for the exact null distribution of LRT and ReLRT
  nmlp = nmlp      ## the number of MLP mlpoximation of the exact null distribution, nmlp needs to be equal to or less than nsim, and can be a vector
  ### method can be ML and REML
  q0 = dim(G)[2]
  p0 = dim(X)[2]
  n  = length(y)
  colnames(X)=paste("X",seq(1:p0),sep="")
  colnames(G)=paste("G",seq(1:q0),sep="")
  propmaf   = colMeans(G)/2
  G[, which(propmaf>0.5)]=2-G[,which(propmaf>0.5)]
  maf        = ifelse(propmaf<0.5,propmaf,1-propmaf)
  maf[(maf==0)]=1e-5
  W  =dbeta(maf,a1,a2)/max(dbeta(maf,a1,a2))
  GW=G%*%diag(W,ncol=length(W),nrow=length(W))
  ONE=as.factor(rep(1,n))
  dataset=data.frame(y,X,ONE)
  #####################
  if (p0==1) 
  {
    formu1 = paste(1,collapse = "+")
  }
  else 
  {
    xnames = paste0("X", 2:p0)
    formu1 = paste(xnames,collapse = "+")
  }
  formu= as.formula(paste("y~",formu1))
  dataset$RE=GW
  ## fit the null model
  mH0= lm(formu,data=dataset)
  mui=(svd(t(qr.resid(qr(X), GW)), nu = 0, nv = 0)$d)^2 ###eigenvalue: mu
  xi =(svd(t(GW), nu = 0, nv = 0)$d)^2                  ###eigenvalue: epislon 
  
  ##########RELRT
  if (method=="REML")
  {
    mH1=lme(formu,data=dataset,random=list(ONE=pdIdent(~RE-1)),method="REML") 
    simLRT=RLRTSim(X=X,Z=GW,qrX=qr(X),sqrt.Sigma=diag(q0),nsim=nsim)
    obsLRT=2*(logLik(mH1,REML=T)-logLik(mH0,REML=T))
    obsLRT=ifelse(obsLRT<0,0,obsLRT)
    pvalue.mlp=matrix(NA,length(nmlp),1)
    AUDresult=matrix(NA,length(nmlp),2)
    for (i in 1:length(nmlp))
    {
      AUD=ReLRTn.MLP(simLRT[1:nmlp[i]],ReLRT.test=obsLRT,mui)
      pvalue.mlp[i,]=c(AUD$pvalue)
      AUDresult[i,]=c(AUD$a,AUD$p)
    }
  }
  
  ##########LRT
  else 
  {
    mH1=lme(formu,data=dataset,random=list(ONE=pdIdent(~RE-1)),method="ML") 
    simLRT=LRTSim(X=X,Z=GW,q=0,sqrt.Sigma=diag(q0),nsim=nsim)
    obsLRT=2*(logLik(mH1,REML=F)-logLik(mH0,REML=F))
    obsLRT=ifelse(obsLRT<0,0,obsLRT)
    pvalue.mlp=matrix(NA,length(nmlp),1)
    AUDresult=matrix(NA,length(nmlp),2)
    for ( i in 1:length(nmlp))
    {
      AUD=LRTn.MLP(simLRT[1:nmlp[i]],LRT.test=obsLRT,n,p0,mui,xi)
      pvalue.mlp[i,]=c(AUD$pvalue)
      AUDresult[i,]  =c(AUD$a,AUD$p)
    }
  }
  
  pvalue=mean(simLRT>=obsLRT)
  pvalue.mlp=pvalue.mlp
  AUDresult=AUDresult
  fixbeta=mH1$coef$fixed
  ranbeta=mH1$coef$random$ONE
  covvar=as.numeric(VarCorr(mH1)[q0:(q0+1),1])
  object = list(obsLRT=obsLRT,pvalue=pvalue,pvalue.mlp=pvalue.mlp,
                covvar=covvar,AUDresult=AUDresult,fixbeta=fixbeta,
                ranbeta=ranbeta)
}


########### cumulative probability
pchibarsqAUD=function(p,df=1,mix=0.5,a=1,lower.tail=TRUE,log.p=FALSE) 
{
  p=p/a
  df=rep(df,length.out=length(p))
  mix=rep(mix,length.out=length(p))
  c1=ifelse(df==1, if (lower.tail) 
    1
    else 0, pchisq(p,df-1,lower.tail=lower.tail))
  c2=pchisq(p,df,lower.tail=lower.tail)
  r=mix*c1+(1-mix)*c2
  if (log.p) 
    log(r)
  else r
}

################## MLP mlpoximation for LRT
LRTn.MLP=function(LRT,LRT.test,n,p0,mui,xi)
{
  x=LRT
  pQ=1-davies((n-p0)*sum(xi)/n,lambda=mui)$Qq
  apQ=mean(x)/(1-pQ)
  p.pQ =1-pchibarsqAUD(LRT.test,df=1,mix=pQ,a=apQ)
  object = list(pvalue=p.pQ,a=apQ,p=pQ)
}

################## MLP mlpoximation for ReLRT
ReLRTn.MLP=function(ReLRT,ReLRT.test,mui)
{
  x=ReLRT
  pQ=1-davies(sum(mui),lambda=mui)$Qq
  apQ=mean(x)/(1-pQ)
  p.pQ =1-pchibarsqAUD(ReLRT.test,df=1,mix=pQ, a=apQ)
  object = list(pvalue=p.pQ,a=apQ,p=pQ)
}

