#### Appendix 3: This script uses red gorgonian data to illustrate one method
#### to simultaneously fit the minimum, maximum, mean, and precision of the beta model
# Note: this script was run with R v. 3.6.0

## Note: we use a constrained search to prevent biologically implausible bounds
## and to ensure comparison to the same data set (e.g., not excluding additional data points)

rm(list=ls())
require(quantreg)
require(betareg)
require(dplyr)
require(MuMIn)
require(MASS)
require(MBESS)
require(gamlss)

########################################################################
##### Functions related to beta approach

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

#########################################################################
### Load, manipulate, and sort data
dat <- read.csv("Gorgonian raw data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA 
minsize=0 
maxsize=75 # this is set manually to not have a large category at end that has no plants.
size$t0[which(size$t0>maxsize)]=maxsize-0.1
size$t1[which(size$t1>maxsize)]=maxsize-0.1
size <- size[-which(is.na(size$t0)),]
size <- size[-which(is.na(size$t1)),]
size$t02<-size$t0^2
size$t03<-size$t0^3

#############################################################

## Step 1: get covariance matrix to use for sampling parameters 
## for the minimum and maximum sizes
quant<-rq(t1~t0+I(t0^2)+I(t0^3),data=size,tau=c(0.001,0.999))
quantcovmin<-summary.rq(rq(t1~t0+I(t0^2)+I(t0^3),data=size,tau=c(0.001)),covariance=T)$cov 
quantcovmax<-summary.rq(rq(t1~t0+I(t0^2)+I(t0^3),data=size,tau=c(0.999)),covariance=T)$cov
minsizes=predict(quant,size)[,1] # size dependent minimum
maxsizes=predict(quant,size)[,2] # size dependent maximum

plot(size$t0,size$t1,col='grey',ylab='Size at time t+1',xlab='Size at time t')
lines(0:100,predict(quant,data.frame(t0=0:100))[,1])
lines(0:100,predict(quant,data.frame(t0=0:100))[,2])

## Step 2: transform to (0,1) interval 
size$t1b<-betaFn(size$t1,min=minsizes,max=maxsizes)
# adjust any data points outside the bounds to just within the bounds
size$t1b[size$t1b>0.99]<-0.99
size$t1b[size$t1b<0.01]<-0.01
# back-transform with adjusted data points; this will be used to ensure
# new parameters don't exclude any additional data points
size$t1n<-backbeta(size$t1b,min=minsizes,max=maxsizes)

## Step 3: sample parameter values for max and min sizes 
## from their covariance matrices 
# Note: this constrains the parameter search space 
ndraw<-1000000 
# choose a large number of draws because many parameter sets will be excluded
# if they exclude additional data points
plo<-mvrnorm(ndraw,mu=coef(quant)[,1],Sigma=quantcovmin)
phi<-mvrnorm(ndraw,mu=coef(quant)[,2],Sigma=quantcovmax)

## Step 4: for each set of parameter values, transform data to (0,1) interval
## and, if not excluding any more data points, fit the beta-regression
## then get the log likelihood of the entire model to the data
LL<-numeric(ndraw)
for(i in 1:ndraw){
  mins<-plo[i,1]+plo[i,2]*size$t0+plo[i,3]*size$t02+plo[i,4]*size$t03 # new min sizes
  maxs<-phi[i,1]+phi[i,2]*size$t0+phi[i,3]*size$t02+phi[i,4]*size$t03 # new max sizes
  nout<-length(c(size$t1n[size$t1n>=maxs],size$t1n[size$t1n<=mins])) # number of data points outside new interval
  if(nout==0){ # only consider parameters that don't exclude additional data points
    size$t1bfit<-betaFn(size$t1n,mins,maxs) # transform to (0,1) interval using new bounds
    # fit and compare beta models
    b1<-gamlss(t1bfit~cs(t0,3),sigma.formula=~cs(t0,3),data=size,family=BE,method=RS(250))
    b2<-betareg(t1bfit~t0+I(t0^2)|t0+I(t0^2),data=size)
    b3<-betareg(t1bfit~t0|t0+I(t0^2),data=size)
    b4<-betareg(t1bfit~t0+I(t0^2)|t0,data=size)
    b5<-betareg(t1bfit~t0+I(t0^2),data=size)
    mbAIC<-unlist(lapply(list(b1,b2,b3,b4,b5),AIC))
    betamod<-list(b1,b2,b3,b4,b5)[[which.min(mbAIC)]]
    LL[i]<-logLik(betamod)-sum(log(maxs-mins)) # likelihood of the new model
  } else {
    LL[i]<-NA
  }
}

length(LL[!is.na(LL)]) # number of parameter combinations tested

## Step 5: find the parameters of the sample that maximize the log likelihood
# Note: these are not the maximum likelihood estimates because the search is constrained
best<-which.max(LL)
plo[best,];phi[best,]
mins<-plo[best,1]+plo[best,2]*size$t0+plo[best,3]*size$t02+plo[best,4]*size$t03
maxs<-phi[best,1]+phi[best,2]*size$t0+phi[best,3]*size$t02+phi[best,4]*size$t03

plot(size$t0,size$t1,col='grey',ylab='Size at time t+1',xlab='Size at time t',ylim=c(-3,75))
lines(0:100,predict(quant,data.frame(t0=0:100))[,1])
lines(0:100,predict(quant,data.frame(t0=0:100))[,2])
lines(0:100,plo[best,1]+plo[best,2]*0:100+plo[best,3]*(0:100)^2+plo[best,4]*(0:100)^3,lty=2)
lines(0:100,phi[best,1]+phi[best,2]*0:100+phi[best,3]*(0:100)^2+phi[best,4]*(0:100)^3,lty=2)
legend('topleft',lty=c(1,2),c('Original','Updated'))

## Step 6: Compare the AIC of the old and new models
# fit new model
size$t1bfit<-betaFn(size$t1n,mins,maxs)
b1<-gamlss(t1bfit~cs(t0,3),sigma.formula=~cs(t0,3),data=size,family=BE,method=RS(250))
b2<-betareg(t1bfit~t0+I(t0^2)|t0+I(t0^2),data=size)
b3<-betareg(t1bfit~t0|t0+I(t0^2),data=size)
b4<-betareg(t1bfit~t0+I(t0^2)|t0,data=size)
b5<-betareg(t1bfit~t0+I(t0^2),data=size)
mbAIC<-unlist(lapply(list(b1,b2,b3,b4,b5),AIC))
betamod<-list(b1,b2,b3,b4,b5)[[which.min(mbAIC)]]
LLNew<-logLik(betamod)-sum(log(maxs-mins)) # log likelihood on the original data scale
AICNew<- -2*LLNew + 2*(length(coef(quant))+betamod$df.fit)

# Original model
b1<-gamlss(t1b~cs(t0,3),sigma.formula=~cs(t0,3),data=size,family=BE,method=RS(250))
b2<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),data=size)
b3<-betareg(t1b~t0|t0+I(t0^2),data=size)
b4<-betareg(t1b~t0+I(t0^2)|t0,data=size)
b5<-betareg(t1b~t0+I(t0^2),data=size)
mbAIC<-unlist(lapply(list(b1,b2,b3,b4,b5),AIC))
betamod<-list(b1,b2,b3,b4,b5)[[which.min(mbAIC)]]
LLOrig<-logLik(betamod)-sum(log(maxsizes-minsizes)) # log likelihood on the original data scale
AICOrig<- -2*LLOrig + 2*(length(coef(quant))+betamod$df.fit)

