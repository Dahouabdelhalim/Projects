##### Appendix 2: This script uses demographic data for red gorgonian coral
##### to compare IPMs with different growth models to PPMs
# Note: this script was run with R v. 3.6.0

rm(list=ls(all=TRUE))
library(dplyr)
library(popbio)
library(MuMIn)
library(binr)
library(quantreg)
library(betareg)
library(sn)
library(gamlss)
library(abind)
library(reshape2)
library(zoo)

########################################################################
##### Functions 

# This function transforms continuous size data to a (0,1) interval
# based on size-dependent min and max values
betaFn<-function(x,min,max){
  y=(x-min)/(max-min)
}

# This function back-transforms data from a (0,1) interval
# to the original data scale based on size-dependent min and max values
backbeta<-function(x,min,max){
  (x*(max-min))+min
}

# This function back-transforms variance from a (0,1) interval
# to the original data scale based on size-dependent min and max values
backvar<-function(sigma,min,max){
  sigma*((max-min)^2)
}

# This function calculates the alpha and beta shape parameters of a 
# beta distribution based on the mean and variance
shapebeta<-function(mu,sigma){
  alpha<-((1-mu)/sigma-1/mu)*mu^2
  beta<-alpha*(1/mu-1)
  return(params=list(alpha=alpha,beta=beta))
}

# This function calculates the mean and variance of a beta distribution
# based on the alpha and beta shape parameters
revshapebeta<-function(alpha,beta){
  mu<-alpha/(alpha+beta)
  sigma<-(alpha*beta)/(((alpha+beta)^2)*(alpha+beta+1))
  return(params=list(mu=mu,sigma=sigma))
}

# This function calculates the skewness of a beta distribution
# based on the alpha and beta shape parameters
skewbeta<-function(alpha,beta){
  (2*(beta-alpha)*sqrt(alpha+beta+1))/((alpha+beta+2)*sqrt(alpha*beta))
}

# This function calculates the lifespan as the time for a new recruit
# to have less than 1% probability of still being alive
lifespan <- function(nx){
  nclasses=dim(nx)[1]
  vec=c(100,rep(0,(nclasses-1)))
  nx[1,]=0
  jj=1
  while (sum(vec)>1){
    vec=nx%*%vec
    jj=jj+1
    print(sum(vec))
  }
  return(jj)
}

# This function calculates the skewness of a skewed normal
# distribution as a function of the mu (xi), sigma (omega), and nu (alpha) parameters
skewsn<-function(mu,sigma,nu){
  c=nu/sqrt(1+nu^2)
  ((4-pi)/2)*(((c*sqrt(2/pi))^3)/((1-(2*c^2)/pi)^(3/2)))
}


#########################################################################
### Load, manipulate, and sort data

# reproductive rates
sizes <- c(15,25,35,45) # 45 is actually >40, and I replaced 0's with 0.001 
gonads <- c(1100,31734,135473,291871) # these are actually oocyte #s, from como et al 1995 table 5.  
correctionfromoocytes = 2.77774E-06

# constant reproduction model
nlsmod=nls(gonads~b*(sizes^c),start=list(b=4557,c=1))
nlpred=predict(nlsmod)

# function to estimate reproduction for a given size
Gf_i <- function(x) {
  dataf <- as.data.frame(x)
  names(dataf) <- c("sizes")
  repall <- correctionfromoocytes *predict(nlsmod, dataf)
  return(repall)
}

# demographic data
dat <- read.csv("Gorgonian raw data.csv")

# note all surviving size-class-1 polyps recruit into size class 2
recruitsurv = mean(0.667, 0.636,0.750,0.769, 0.667) # first-class polyps survive at the same rate as does size class 2, 
recruitsurv=recruitsurv # the mean of size class 2 survival rates from Linares et al. 2007 
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$repro <- Gf_i(size$t0) # number of successful recruits that get to first size class
size$t1[which(size$survival==0)] <- NA 
minsize <- 0 
maxsize <- 75 # this is set manually to not have a large category at end that has no individuals.
size$t0[which(size$t0>maxsize)]=maxsize-0.1
size$t1[which(size$t1>maxsize)]=maxsize-0.1
size <- size[-which(is.na(size$t0)),]
size$t0sq <- size$t0^2 

#############################################
### Construct IPMs and PPMs for 200 bootstrapped datasets

Nboot=200
N=nrow(size)

# store matrices
BetaIPMs=SkewIPMs=NormIPMs=DiscreteMXs=c() 
# store moments of the estimated growth distributions
bmean=nmean=snmean=dmean=list() 
bvar=nvar=snvar=dvar=list()
bskew=snskew=dskew=list()
# store reactivities
breac=nreac=snreac=dreac=list()

# PPM bins
n.bin=50
stpsize=(60-minsize)*(1/(n.bin-1)) # 49 evenly divided size classes plus 1 larger size class (60-75) to ensure adequate sample sizes

for (i in 1:Nboot){
  # bootstrap the data set
  df<-size[sample(1:N,N,replace=TRUE),]
  # fit the quantile regression
  mqa<-rq(t1~t0,data=df,tau=c(0.001,0.999))
  mqb<-rq(t1~t0+I(t0^2),data=df,tau=c(0.001,0.999))
  mqc<-rq(t1~t0+I(t0^2)+I(t0^3),data=df,tau=c(0.001,0.999))
  mqAIC<-lapply(list(mqa,mqb,mqc),AIC)
  mqMin<-list(mqa,mqb,mqc)[[which.min(lapply(mqAIC,'[',1))]] # best-supported model for minimum size bound
  mqMax<-list(mqa,mqb,mqc)[[which.min(lapply(mqAIC,'[',2))]] # best-supported model for maximum size bound
  # transform to (0,1) interval 
  MaxSize<-predict(mqMax,df)[,2]
  MinSize<-predict(mqMin,df)[,1]
  df$t1b<-betaFn(df$t1,min=MinSize,max=MaxSize)
  df$t1b[df$t1b>0.99]<-0.99 # set any values outside bounds to just within bounds
  df$t1b[df$t1b<0.01]<-0.01
  ## Model selection for growth models for each distribution
  ## Comparing both parametric and non-parametric (cubic spline) models
  # beta regression model selection
  mba<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),data=df)
  mbb<-betareg(t1b~t0|t0+I(t0^2),data=df)
  mbc<-betareg(t1b~t0+I(t0^2)|t0,data=df)
  mbd<-betareg(t1b~t0|t0,data=df)
  mbe<-betareg(t1b~t0+I(t0^2),data=df)
  mbf<-betareg(t1b~t0,data=df)
  mbg<-gamlss(t1b~cs(t0,3),sigma.formula=~cs(t0,3),data=df[!is.na(df$t1),],family=BE,method=RS(250))
  mbh<-gamlss(t1b~cs(t0,3),data=df[!is.na(df$t1),],family=BE,method=RS(250))
  mbAIC<-AIC(mba,mbb,mbc,mbd,mbe,mbf,mbg,mbh)
  mb<-list(mba,mbb,mbc,mbd,mbe,mbf,mbg,mbh)[[which.min(mbAIC$AIC)]]
  # skewed normal regression model selection
  msna<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),nu.formula=~t0+I(t0^2),data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnb<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),nu.formula=~t0,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnc<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,nu.formula=~t0,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnd<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msne<-gamlss(t1~t0+I(t0^2),sigma.formula=~1,nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnf<-gamlss(t1~t0,sigma.formula=~t0+I(t0^2),nu.formula=~t0,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msng<-gamlss(t1~t0,sigma.formula=~t0,nu.formula=~t0,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnh<-gamlss(t1~t0,sigma.formula=~t0,nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msni<-gamlss(t1~t0,sigma.formula=~1,nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnj<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),nu.formula=~cs(t0,3),data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnk<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnl<-gamlss(t1~cs(t0,3),sigma.formula=~1,nu.formula=~1,data=df[!is.na(df$t1),],family=SN1,method=RS(250))
  msnAIC<-unlist(lapply(list(msna,msnb,msnc,msnd,msne,msnf,msng,msnh,msni,msnj,msnk,msnl),AIC))
  msn<-list(msna,msnb,msnc,msnd,msne,msnf,msng,msnh,msni,msnj,msnk,msnl)[[which.min(msnAIC)]]
  # normal linear regression model selection
  mna<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),data=df[!is.na(df$t1),],family='NO')
  mnb<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,data=df[!is.na(df$t1),],family='NO')
  mnc<-gamlss(t1~t0+I(t0^2),sigma.formula=~1,data=df[!is.na(df$t1),],family='NO')
  mnd<-gamlss(t1~t0,sigma.formula=~t0+I(t0^2),data=df[!is.na(df$t1),],family='NO')
  mne<-gamlss(t1~t0,sigma.formula=~t0,data=df[!is.na(df$t1),],family='NO')
  mnf<-gamlss(t1~t0,sigma.formula=~1,data=df[!is.na(df$t1),],family='NO')
  mng<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),data=df[!is.na(df$t1),],family='NO')
  mnh<-gamlss(t1~cs(t0,3),data=df[!is.na(df$t1),],family='NO')
  mnAIC<-unlist(lapply(list(mna,mnb,mnc,mnd,mne,mnf,mng,mnh),AIC))
  mn<-list(mna,mnb,mnc,mnd,mne,mnf,mng,mnh)[[which.min(mnAIC)]]
  # fit other vital rate models 
  sur_models <- list(glm(survival~ t0, family= "binomial",data=df),
                     glm(survival~ t0 + t0sq , family= "binomial",data=df))
  min_AIC <- min(AIC(sur_models[[1]], sur_models[[2]])$AIC)
  bestsur <- sur_models[[which(min_AIC == AIC(sur_models[[1]], sur_models[[2]])$AIC)]] 
  bestrep= nlsmod
  # create IPMs
  vec.bin = c(seq(minsize, 60, stpsize/2),seq(60,maxsize,length.out=round((maxsize-60)/(stpsize/2)))[-1]) # size bins that evenly divide the bins used by the PPM
  bin.num = length(vec.bin)-1
  vec.bin[length(vec.bin)] = maxsize
  binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)]) # mesh points
  indata <- as.data.frame(cbind(binmids, binmids^2))
  names(indata) <- c("t0", "t0sq")
  inreprodata=as.data.frame(binmids)
  names(inreprodata)=c("sizes")
  # predict survival and reproduction rates
  sur_vals <- predict(bestsur,indata, type='response')
  rep_vals <- predict(nlsmod, inreprodata, type='response')*correctionfromoocytes 
  rep_vals[rep_vals<0] <- 0
  ## predict growth transition rates for each distribution using best-supported models
  ## Note: 'safe prediction' warnings from gamlss caused by predicting for wider range of 
  ## values than used to fit the models (i.e., extrapolation)
  # beta growth model
  ifelse(class(mb)[1]=='gamlss',
    beta_params<-list(mu=predict(mb,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('mu')),
                sigma=predict(mb,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('sigma')),
                mins=predict(mqMin,indata)[,1],
                maxs=predict(mqMax,indata)[,2]),
    beta_params<-c(shapebeta(mu=predict(mb,newdata=indata,type='response'),sigma=predict(mb,newdata=indata,type='variance')),
                list(mins=predict(mqMin,indata)[,1],
                maxs=predict(mqMax,indata)[,2])))
  # skewed normal growth model
  sn_params<-list(mu=predict(msn,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('mu')),
                  sigma=predict(msn,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('sigma')),
                  nu=predict(msn,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('nu')))
  # normal growth model
  n_params<-list(mu=predict(mn,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('mu')),
                 sigma=predict(mn,newdata=indata,data=df[!is.na(df$t1),],type='response',what=c('sigma')))
  # plot fit to the data
  pn<-centiles.pred(mn,xname='t0',xvalues=indata$t0,data=df[!is.na(df$t1),],cent=c(1,50,99))
  psn<-centiles.pred(msn,xname='t0',xvalues=indata$t0,data=df[!is.na(df$t1),],cent=c(1,50,99))
  ifelse(class(mb)[1]=='gamlss',
    pb<-centiles.pred(mb,xname='t0',xvalues=indata$t0,data=df[!is.na(df$t1),],cent=c(1,50,99)),
    pb<-predict(mb,indata,type='quantile',c(0.01,0.5,0.99)))
  pb[,2:4]<-apply(pb[,2:4],2,backbeta,min=beta_params$mins,max=beta_params$maxs)
  plot(df$t0,df$t1,col='grey')
  lines(indata$t0,pn[,3],col='blue');lines(indata$t0,pn[,2],lty=2,col='blue');lines(indata$t0,pn[,4],lty=2,col='blue')
  lines(indata$t0,psn[,3],col='green');lines(indata$t0,psn[,2],lty=2,col='green');lines(indata$t0,psn[,4],lty=2,col='green')
  lines(indata$t0,pb[,3],col='red');lines(indata$t0,pb[,2],lty=2,col='red');lines(indata$t0,pb[,4],lty=2,col='red')
  ## Create growth transition matrices for each distribution
  gmx_beta=gmx_sn=gmx_n=matrix(data=NA, nrow=bin.num, ncol=bin.num)
  n<-length(vec.bin)
  if(class(mb)[1]=='gamlss') { # get shape parameters for dbeta from gamlss parameters
      beta_params$alpha=beta_params$mu*(1-beta_params$sigma^2)/(beta_params$sigma^2)
      beta_params$beta=beta_params$alpha*(1-beta_params$mu)/beta_params$mu
  }
  for (ss in 1:(bin.num)) {
    zbeta=betaFn(x=vec.bin,min=beta_params$mins[ss],max=beta_params$maxs[ss])
    # use the difference in the cdf for each bin edge to get probabilities of growing to each bin
    # and ensure that growth probabilities sum to 1
    grow_beta<-pbeta(zbeta,shape1=beta_params$alpha[ss],shape2=beta_params$beta[ss])
    grow_beta<-grow_beta[2:n]-grow_beta[1:(n-1)]
    ifelse(sum(grow_beta)>0,gmx_beta[,ss]<-grow_beta/sum(grow_beta),gmx_beta[,ss]<-NA)
    grow_skew<-psn(vec.bin,xi=sn_params$mu[ss],omega=sn_params$sigma[ss],alpha=sn_params$nu[ss])
    grow_skew<-grow_skew[2:n]-grow_skew[1:(n-1)]
    ifelse(sum(grow_skew)>0,gmx_sn[,ss]<-grow_skew/sum(grow_skew),gmx_sn[,ss]<-NA)
    grow_norm<-pnorm(vec.bin,mean=n_params$mu[ss],sd=n_params$sigma[ss])
    grow_norm<-grow_norm[2:n]-grow_norm[1:(n-1)]
    ifelse(sum(grow_norm)>0,gmx_n[,ss]<-grow_norm/sum(grow_norm), gmx_n[,ss]<-NA)
  } # end ss loop
  
  # make the survival/growth and reproduction matrices
  survmx <- t(matrix(rep(sur_vals,(bin.num)),(bin.num)))
  reprow <- rep_vals

  mx <- matrix(0, (bin.num+1), (bin.num+1))
  mx[2,1] = recruitsurv
  mx[1,2:(bin.num+1)] = reprow 
  mxbeta=mxsn=mxn=mx
  # separate survival/growth matrices for each distribution
  mxbeta[2:(bin.num+1), 2:(bin.num+1)] = gmx_beta*survmx
  mxsn[2:(bin.num+1), 2:(bin.num+1)] = gmx_sn*survmx
  mxn[2:(bin.num+1), 2:(bin.num+1)] = gmx_n*survmx
  # store the IPM matrices
  BetaIPMs<-abind(BetaIPMs,mxbeta,along=3)
  SkewIPMs<-abind(SkewIPMs,mxsn,along=3)
  NormIPMs<-abind(NormIPMs,mxn,along=3)
  
  ## get the estimated moments for each growth distribution
  # beta growth model
  if(class(mb)[1]=='gamlss'){
    sigma=(beta_params$sigma^2)*beta_params$mu*(1-beta_params$mu) # this is the variance
    params=shapebeta(mu=beta_params$mu,sigma=sigma)
    # get the moments on the original data scale, not the (0,1) interval
    bmean[[i]]<-backbeta(x=beta_params$mu,min=beta_params$mins,max=beta_params$maxs)
    bvar[[i]]<-backvar(sigma=sigma,min=beta_params$mins,max=beta_params$maxs)
    bskew[[i]]<-skewbeta(alpha=params$alpha,beta=params$beta)
  } else if(class(mb)[1]=='betareg'){
    params=shapebeta(mu=beta_params$mu,sigma=beta_params$sigma)
    # get the moments on the original data scale, not the (0,1) interval
    bmean[[i]]<-backbeta(x=beta_params$mu,min=beta_params$mins,max=beta_params$maxs)
    bvar[[i]]<-backvar(sigma=beta_params$sigma,min=beta_params$mins,max=beta_params$maxs)
    bskew[[i]]<-skewbeta(alpha=params$alpha,beta=params$beta)
  }
  # skewed normal growth model
  snmean[[i]]<-sn_params$mu+sn_params$sigma*sign(sn_params$nu)*sqrt((2*sn_params$nu^2)/(pi*(1+sn_params$nu^2))) # this is the mean
  snvar[[i]]<-(sn_params$sigma^2)*(1-((2*sn_params$nu^2)/(pi*(1+sn_params$nu^2)))) # this is the variance
  snskew[[i]]<-skewsn(mu=sn_params$mu,sigma=sn_params$sigma,nu=sn_params$nu)
  # normal growth model
  nmean[[i]]<-n_params$mu
  nvar[[i]]<-n_params$sigma^2
  
  ## make the discrete PPM
  ss<-df$t0
  n.bin=50
  vec.bin = seq(minsize, 60, stpsize)
  vec.bin = c(vec.bin,maxsize) # 1 large bin at the end to ensure adequate sample sizes with bootstrapping
  nums=hist(ss,breaks=vec.bin, plot=FALSE)$counts 
  if (min(nums)>0) {
    surv <- rep(NA, n.bin)                      # survivorship for each class
    grow <- matrix(NA, n.bin, n.bin)            # store growth probabilites for each class
    reproduction <- rep(NA, n.bin)
    
    for(j in 1:(length(vec.bin)-1)){
      # set limits for subset according to bin edges
      bounds <- c(vec.bin[j], vec.bin[j+1])
      # subset data according to bounds
      subset <- df[df$t0 > bounds[1] & df$t0 <= bounds[2],]
      # calculate survivorship for this class
      surv[j] <- sum(subset$survival, na.rm=TRUE) / length(subset$t0[which(!is.na(subset$survival))]) 
      # store histo as object, to access counts per bin
      histo <- hist(subset$t1, breaks = vec.bin, plot = FALSE)
      # $counts returns the number of individuals of a certain size class
      grow[,j] <- histo$counts/length(subset$t0[subset$survival==1]) 
      reproduction[j] <- mean(subset$repro)
    }
    
    M1 <- matrix(NA, n.bin, n.bin)   # initiate projection matrix
    M <- matrix(0, (n.bin+1), (n.bin+1))
    # populate survival/growth matrix
    for(j in 1:length(surv)) M1[,j] <- surv[j] * grow[,j]
    # add the creation of babies and their transition to first size class
    M[2:(n.bin+1), 2:(n.bin+1)] = M1
    M[2,1] = recruitsurv
    M[1,2:(n.bin+1)] = reproduction 
  }  
  DiscreteMXs<-abind(DiscreteMXs,M,along=3)

  ## Get the reactivity corrected for matrix dimension
  ## using an initial population vector of the largest PPM size class
  ## or the proportional representation of the narrower IPM size classes within it
  wts<-hist(df$t0[df$t0>60],breaks=seq(60,maxsize,length.out=round((maxsize-60)/(stpsize/2))), plot=FALSE)$counts
  v<-c(rep(0,bin.num-length(wts)+1),wts/sum(wts))
  # one-time step population growth relative to lambda
  dreac[[i]]<-sum(M%*%c(rep(0,n.bin),1))/lambda(M)
  breac[[i]]<-sum(mxbeta%*%v)/lambda(mxbeta)
  snreac[[i]]<-sum(mxsn%*%v)/lambda(mxsn)
  nreac[[i]]<-sum(mxn%*%v)/lambda(mxn)
}

save.image("~Gorgonian bootstrap.RData")

#######################################################
########### Look at the results

# lambda
Dlams<-apply(DiscreteMXs,3,lambda)
Blams<-apply(BetaIPMs,3,lambda)
Slams<-apply(SkewIPMs,3,lambda)
Nlams<-apply(NormIPMs,3,lambda)
Lams<-melt(cbind(Dlams,Blams,Slams,Nlams))
boxplot(Lams$value~Lams$Var2,names=c('Discrete','Beta','Skewed','Normal'),ylab='Lambda',main='Gorgonian')

# lifespan
Dls<-apply(DiscreteMXs,3,lifespan)
Bls<-apply(BetaIPMs,3,lifespan)
Sls<-apply(SkewIPMs,3,lifespan)
Nls<-apply(NormIPMs,3,lifespan)
Life<-melt(cbind(Dls,Bls,Sls,Nls))
boxplot(Life$value~Life$Var2,names=c('Discrete','Beta','Skewed','Normal'),ylab='Life span',main='Gorgonian')

# Reactivity
Reac<-melt(cbind(unlist(dreac),unlist(breac),unlist(snreac),unlist(nreac)))
boxplot(Reac$value~Reac$Var2,names=c('Discrete','Beta','Skewed','Normal'),ylab='Reactivity',main='Gorgonian')

### Get actual moments of the dataset
### using a rolling window of 50 data points
t0<-size$t0[order(size$t0)]
t1<-size$t1[order(size$t0)]
startvals<-rollapply(t0,width=50,mean)
realmeans<-rollapply(t1,width=50,mean,na.rm=T)
realvars<-rollapply(t1,width=50,var,na.rm=T)
realskews<-rollapply(t1,width=50,skewness,na.rm=T)

# compare estimated moments for the same size classes
vec.bin = c(seq(minsize, 60, stpsize/2),seq(60,maxsize,length.out=round((maxsize-60)/(stpsize/2)))[-1])  
vec.bin[length(vec.bin)] = maxsize
binmids <- 0.5*(vec.bin[2:length(vec.bin)] + vec.bin[1:(length(vec.bin)-1)])

means<-melt(cbind(unlist(bmean),unlist(snmean),unlist(nmean)))
means$Var1<-binmids[means$Var1]
means$Var2<-c('beta','skew','norm')[means$Var2]

vars<-melt(cbind(unlist(bvar),unlist(snvar),unlist(nvar)))
vars$Var1<-binmids[vars$Var1]
vars$Var2<-c('beta','skew','norm')[vars$Var2]

skews<-melt(cbind(unlist(bskew),unlist(snskew)))
skews$Var1<-rep(binmids,200)[skews$Var1]
skews$Var2<-c('beta','skew')[skews$Var2]

boxplot(value~Var2*Var1,data=means,outline=FALSE,whisklty=0,staplelty=0,
        col=c('red','green','blue'),ylim=c(0,81),xaxt='n',pars=list(boxwex=1.5),
        ylab='Mean',xlab='Size at time t')
axis(side=1,at=seq(1,length(binmids)*3,3),round(binmids,0))
lines(seq(1,length(binmids)*3,3),predict(smooth.spline(startvals,realmeans,spar=1),x=binmids)$y,lwd=2,col='grey',type='l')
legend('topleft',c('Beta','Skewed Normal','Normal'),col=c('red','green','blue'),pch=15)

boxplot(value~Var2*Var1,data=vars,outline=FALSE,whisklty=0,staplelty=0,
        col=c('red','green','blue'),ylim=c(0,170),xaxt='n',pars=list(boxwex=1.5),
        ylab='Variance',xlab='Size at time t')
axis(side=1,at=seq(1,length(binmids)*3,3),round(binmids,0))
lines(seq(1,length(binmids)*3,3),predict(smooth.spline(startvals,realvars,spar=1),x=binmids)$y,lwd=2,col='grey',type='l')

boxplot(value~Var2*Var1,data=skews,outline=FALSE,whisklty=0,staplelty=0,
        col=c('red','green'),ylim=c(-3,1.5),xaxt='n',pars=list(boxwex=1.5),
        ylab='Skew',xlab='Size at time t')
axis(side=1,at=seq(1,length(binmids)*2,2),round(binmids,0))
lines(seq(1,length(binmids)*2,2),predict(smooth.spline(startvals,realskews,spar=1),x=binmids)$y,lwd=2,col='grey',type='l')
lines(abline(h=0,col='blue',lwd=2))
