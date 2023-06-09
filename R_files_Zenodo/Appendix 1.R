#### Appendix 1: This script uses demographic data on the growth of red gorgonian coral 
#### to compare various approaches to modeling growth data.
# Note: this script was run with R v. 3.6.0

rm(list=ls(all=TRUE))
library(quantreg)
library(betareg)
library(moments)
library(zoo)
library(dplyr)
library(MuMIn)
library(gamlss)
library(sn)
library(ggplot2)

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

#########################################################################
### Load, manipulate, and sort data
dat <- read.csv("Gorgonian raw data.csv")
dat.size <- as.data.frame(dat)
colnames(dat.size)[6:8] <- c("Size"="t0", "Sizenext"="t1", "Survnext"= "survival") # rename column header
size <- arrange(dat.size, t0) # arrange initial column in ascending orders
size$t1[which(size$survival==0)] <- NA # remove individuals that died from growth analyses
minsize=0 
maxsize=75 # this is set manually to not have a large category at end that has no plants.
size$t0[which(size$t0>maxsize)]=maxsize-0.1
size$t1[which(size$t1>maxsize)]=maxsize-0.1
size <- size[-which(is.na(size$t0)),] # remove individuals with missing size data
size <- size[-which(is.na(size$t1)),]

## Look for departures from normality for a given starting size
# sliding window of the mean and skew for 50 data points
rollmean<-rollapply(size$t0,width=50,mean,na.rm=T)
rollskew<-rollapply(size$t1,width=50,skewness,na.rm=T)
plot(rollmean,rollskew,xlab='Size at time t',ylab='Skewness in size t+1',type='l')
abline(h=0,col='grey') # skew shifts from positive to negative with size

##################################################################
#### Start with the beta approach 

## Fit size-dependent minimum and maximum values
# Note: the choice of quantiles will depend on data availability
# and potential for outliers and/or measurement error
quant1<-rq(t1~t0,data=size,tau=c(0.001,0.999)) # fit quantile regression
quant2<-rq(t1~t0+I(t0^2),data=size,tau=c(0.001,0.999)) # fit quantile regression
quant3<-rq(t1~t0+I(t0^2)+I(t0^3),data=size,tau=c(0.001,0.999)) # fit quantile regression
lapply(list(quant1,quant2,quant3),AIC) # cubic model is best-supported for both quantiles
minsizes=predict(quant3,size)[,1] # size dependent minimum
maxsizes=predict(quant3,size)[,2] # size dependent maximum

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,predict(quant3,data.frame(t0=0:100))[,1],col='black')
lines(0:100,predict(quant3,data.frame(t0=0:100))[,2],col='black')

## Transform size at time t+1 to (0,1) interval
size$t1b<-betaFn(x=size$t1,min=minsizes,max=maxsizes)
range(size$t1b,na.rm=T) # produces some values outside (0,1) interval
# setting these values to just within the interval
# alternatives would be to exclude as measurement error or 
# adjust min and max bounds slightly to accomodate all points
size$t1b[size$t1b>0.99]=0.99
size$t1b[size$t1b<0.01]=0.01

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')

## Fit and compare beta regression models
# parametric models
m1<-betareg(t1b~t0+I(t0^2)|t0+I(t0^2),link='loglog',data=size) # size and size2 effects on mean and precision
m2<-betareg(t1b~t0|t0+I(t0^2),link='loglog',data=size) # size effects on mean and precision, size2 effect on precision
m3<-betareg(t1b~t0+I(t0^2)|t0,link='loglog',data=size) # size effects on mean and precision, size2 effect on mean
m4<-betareg(t1b~t0|t0,link='loglog',data=size) # size effects on mean and precision
m5<-betareg(t1b~t0+I(t0^2),link='loglog',data=size) # size and size2 effects on mean, constant precision
m6<-betareg(t1b~t0,link='loglog',data=size) # size effects on mean, constant precision
# non-parametric models
m7<-gamlss(t1b~cs(t0,3),sigma.formula=~cs(t0,3),data=size,family=BE) # size effects on mean and scale parameter
m8<-gamlss(t1b~cs(t0,3),sigma.formula=~1,data=size,family=BE) # size effects on mean, constant scale parameter 
AIC(m1,m2,m3,m4,m5,m6,m7,m8) # m7 is best supported (m1 is best-supported parametric model)
bestgrowth=m7
bestpargrowth=m1

plot(size$t0,size$t1b,xlab='Size at time t',ylab='Beta-transformed size t+1',col='grey')
lines(0:100,centiles.pred(bestgrowth,xname='t0',xvalues=0:100,cent=c(50))[,2],col='darkred',lty=1,lwd=2)
lines(0:100,centiles.pred(bestgrowth,xname='t0',xvalues=0:100,cent=c(1))[,2],col='darkred',lty=2,lwd=2)
lines(0:100,centiles.pred(bestgrowth,xname='t0',xvalues=0:100,cent=c(99))[,2],col='darkred',lty=2,lwd=2)
# show best-supported parametric model for comparison
lines(0:100,predict(bestpargrowth,data.frame(t0=0:100),type='response'),col='orange',lty=1,lwd=2)
lines(0:100,predict(bestpargrowth,data.frame(t0=0:100),type='quantile',at=c(0.01)),col='orange',lty=2,lwd=2) 
lines(0:100,predict(bestpargrowth,data.frame(t0=0:100),type='quantile',at=c(0.99)),col='orange',lty=2,lwd=2) 

## Get shape parameters and pdf for a given starting size
newSize<-predict(quant3,data.frame(t0=0:100))
# mu and sigma shape parameters for gamlss BE parameterization of the beta distribution
# note: mean = mu, variance = sigma*mu*(1-mu)
betaMu<-predict(bestgrowth,newdata=data.frame(t0=0:100),data=size,type='response',what=c('mu')) 
betaSigma<-predict(bestgrowth,newdata=data.frame(t0=0:100),data=size,type='response',what=c('sigma')) 
# alpha and beta shape parameters for dbeta parameterization
# for the best-supported parametric model for comparison
betaMean<-predict(bestpargrowth,data.frame(t0=0:100),type='response') # mean on (0,1) scale
betaVar<-predict(bestpargrowth,data.frame(t0=0:100),type='variance') # variance on (0,1) scale
params<-shapebeta(mu=betaMean,sigma=betaVar) # alpha and beta parameters for a given starting size
# get PDFs for starting size=45
zb<-seq(0.0001,0.9999,length.out=100)
db<-dBE(zb,mu=betaMu[46],sigma=betaSigma[46]) # nonparametric density
dbparam<-dbeta(zb,shape1=params$alpha[46],shape2=params$beta[46]) # parametric density
zn<-backbeta(x=zb,min=newSize[46,1],max=newSize[46,2]) # back-transform to original data scale
db<-db/(zn[100]-zn[1]) # back-transform to original data scale
dbparam<-dbparam/(zn[100]-zn[1]) # back-transform to original data scale
sub<-size[which(size$t0>43&size$t0<47),] # get actual size at t+1 for starting size between 44 and 48

hist(sub$t1,breaks=20,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(zn,db,type='l',lwd=2,col='darkred')
lines(zn,dbparam,type='l',lwd=2,col='orange')

## Back-transform model predictions to original data scale
betacent<-centiles.pred(bestgrowth,xname='t0',xvalues=0:100,cent=c(1,50,99))
# best-supported parametric model for comparison
betaquant<-predict(bestpargrowth,data.frame(t0=0:100),type='quantile',at=c(0.01,0.99))

plot(size$t0,size$t1,xlab='Size time t',ylab='Size time t+1',col='grey')
lines(0:100,backbeta(betacent[,3],min=newSize[,1],max=newSize[,2]),lty=1,lwd=2,col='darkred')
lines(0:100,backbeta(betacent[,2],min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='darkred')
lines(0:100,backbeta(betacent[,4],min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='darkred')
# show best-supported parametric model for comparison
lines(0:100,backbeta(betaMean,min=newSize[,1],max=newSize[,2]),lty=1,lwd=2,col='orange')
lines(0:100,backbeta(betaquant[,1],min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='orange')
lines(0:100,backbeta(betaquant[,2],min=newSize[,1],max=newSize[,2]),lty=2,lwd=2,col='orange')

## Get AIC on the original data scale
mu<-predict(bestgrowth,data=size,type='response',what=c('mu'))
sigma<-predict(bestgrowth,data=size,type='response',what=c('sigma'))
LikBeta<-log(dBE(size$t1b,mu=mu,sigma=sigma))-log(maxsizes-minsizes) # transform to original datascale
LLBeta<-sum(LikBeta) 
KBeta<-length(coef(quant3))+bestgrowth$df.fit
AICBeta<- -2*LLBeta + 2*KBeta
# compare to AIC for best-supported parametric model
bmean<-predict(bestpargrowth,data=size,type='response')
bvar<-predict(bestpargrowth,data=size,type='variance')
bparams<-shapebeta(bmean,bvar)
LikBetaParam<-log(dbeta(size$t1b,shape1=bparams$alpha,shape2=bparams$beta))-log(maxsizes-minsizes) # transform to original datascale
LLBetaParam<-sum(LikBetaParam)
KBetaParam<-length(coef(quant3))+length(coef(bestpargrowth))
AICBetaParam<- -2*LLBetaParam + 2*KBetaParam

#####################################################################
#### Now compare to the normal approach fit by separate linear
#### regressions of the mean and variance

## Fit and compare linear regressions for mean growth
m1<-lm(t1~t0+I(t0^2),data=size) # size and size2 on mean 
m2<-lm(t1~t0,data=size) # size on mean
m3<-lm(t1~1,data=size) # no size effect
m4<-nls(t1 ~ a + b*(t0^c), data= size, start=list(a= 1, b=1,c=1)) # power function of size on mean

AIC(m1,m2,m3,m4) # m4 is best-supported
bestmean<-m4

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,predict(bestmean,data.frame(t0=0:100),type='response'),col='darkblue',lty=1,lwd=2)

## Get residuals
size$resid<-size$t1-predict(bestmean,size)
size$resid2<-size$resid^2

plot(size$t0,size$resid,xlab='Size at time t',ylab='Residuals in size t+1',col='grey')

## Fit and compare linear regressions for variance in growth
m1<-lm(resid2~t0+I(t0^2),data=size) # size and size2 on variance 
m2<-lm(resid2~t0,data=size) # size on variance
m3<-lm(resid2~1,data=size) # no size effect

AIC(m1,m2,m3) # m1 is best-supported
bestvar<-m1

plot(size$t0,size$resid2,xlab='Size at time t',ylab='Squared residuals in size t+1',col='grey')
lines(0:100,predict(bestvar,data.frame(t0=0:100),type='response'),col='darkblue',lty=1,lwd=2)

## Get parameters and pdf for a given starting size
normMean<-predict(bestmean,data.frame(t0=0:100)) 
normVar<-predict(bestvar,data.frame(t0=0:100)) 

# get pdf for starting size = 45
dn<-dnorm(0:100,mean=normMean[46],sd=sqrt(normVar[46])) 
sub<-size[which(size$t0>43&size$t0<47),] # get actual size at t+1 for starting size between 42 and 48

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(0:100,dn,type='l',lwd=2,col='darkblue')

## Visualize model predictions
normLow<-qnorm(0.01,mean=normMean,sd=sqrt(normVar))
normHi<-qnorm(0.99,mean=normMean,sd=sqrt(normVar))

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,normMean,col='darkblue',lwd=2,lty=1)
lines(0:100,normLow,col='darkblue',lwd=2,lty=2)
lines(0:100,normHi,col='darkblue',lwd=2,lty=2)

## Compare AIC to the exact same data as the beta approach
size$t1n<-backbeta(size$t1b,min=minsizes,max=maxsizes) # this accounts for the adjustment of some data points to inside the size bounds
mn<-predict(bestmean,size)
mv<-predict(bestvar,size)
LikNorm<-log(dnorm(size$t1n,mean=mn,sd=sqrt(mv)))
LLNorm<-sum(LikNorm)
KNorm<-length(c(coef(bestmean),coef(bestvar)))
AICNorm<- -2*LLNorm + 2*KNorm

####################################################################
#### Now compare to the normal approach fit by a single
#### model to jointly estimate the mean and variance

## Fit and compare linear regressions
# parametric models
m1<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),data=size,family=NO) # size and size2 on mean and variance
m2<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,data=size,family=NO) # size and size2 on mean, size on variance 
m3<-gamlss(t1~t0+I(t0^2),sigma.formula=~1,data=size,family=NO) # size and size2 on mean, constant variance 
m4<-gamlss(t1~t0,sigma.formula=~t0+I(t0^2),data=size,family=NO) # size and size2 on variance, size on mean 
m5<-gamlss(t1~t0,sigma.formula=~t0,data=size,family=NO) # size on mean and variance  
m6<-gamlss(t1~t0,sigma.formula=~1,data=size,family=NO) # size on mean, constant variance 
# non-parametric
m7<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),data=size,family=NO)
m8<-gamlss(t1~cs(t0,3),sigma.formula=~1,data=size,family=NO)
AIC(m1,m2,m3,m4,m5,m6,m7,m8) # m7 is best-supported
bestnorm<-m7

## Get parameters and pdf for a given starting size
# get pdf for starting size = 45
normMean2<-predict(bestnorm,newdata=data.frame(t0=0:100),data=size,type='response',what=c('mu'))
normSD<-predict(bestnorm,newdata=data.frame(t0=0:100),data=size,type='response',what=c('sigma'))

dn<-dNO(0:100,mu=normMean2[46],sigma=normSD[46]) 
sub<-size[which(size$t0>44&size$t0<46),] # get actual size at t+1 for starting size between 42 and 48

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(0:100,dn,type='l',lwd=2,col='blue')

## Visualize model predictions
normLow2<-qnorm(0.01,mean=normMean2,sd=normSD)
normHi2<-qnorm(0.99,mean=normMean2,sd=normSD)

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,normMean2,col='blue',lwd=2,lty=1)
lines(0:100,normLow2,col='blue',lwd=2,lty=2)
lines(0:100,normHi2,col='blue',lwd=2,lty=2)

## Compare AIC to the exact same data as the beta approach
normMu<-predict(bestnorm,data=size,type='response',what=c('mu'))
normSigma<-predict(bestnorm,data=size,type='response',what=c('sigma'))
size$t1n<-backbeta(size$t1b,min=minsizes,max=maxsizes) 
LikNorm2<-log(dnorm(size$t1n,mean=normMu,sd=normSigma))
LLNorm2<-sum(LikNorm2)
KNorm2<-bestnorm$df.fit
AICNorm2<- -2*LLNorm2 + 2*KNorm2

################################################################
### Now compare to a skewed normal
## Fit and compare regressions
# parametric models
m1<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),nu.formula=~t0+I(t0^2),data=size,family=SN1,method=RS(250)) # size and size2 on mean, variance, and skew
m2<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0+I(t0^2),nu.formula=~t0,data=size,family=SN1,method=RS(250)) # size and size2 on mean and variance, size on skew
m3<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,nu.formula=~t0,data=size,family=SN1,method=RS(250)) # size and size2 on mean, size on variance and skew
m4<-gamlss(t1~t0+I(t0^2),sigma.formula=~t0,nu.formula=~1,data=size,family=SN1,method=RS(250)) # size and size2 on mean, size on variance, constant skew
m5<-gamlss(t1~t0+I(t0^2),sigma.formula=~1,nu.formula=~1,data=size,family=SN1,method=RS(250)) # size and size2 on mean, constant variance and skew
m6<-gamlss(t1~t0,sigma.formula=~t0+I(t0^2),nu.formula=~t0,data=size,family=SN1,method=RS(250)) # size on mean and skew, size and size2 on variance
m7<-gamlss(t1~t0,sigma.formula=~t0,nu.formula=~t0,data=size,family=SN1,method=RS(250)) # size on mean, variance, and skew
m8<-gamlss(t1~t0,sigma.formula=~t0,nu.formula=~1,data=size,family=SN1,method=RS(250)) # size on mean and variance, constant skew
m9<-gamlss(t1~t0,sigma.formula=~1,nu.formula=~1,data=size,family=SN1,method=RS(250)) # size on mean, constant variance and skew
# non-parametric
m10<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),nu.formula=~cs(t0,3),data=size,family=SN1,method=RS(250))
m11<-gamlss(t1~cs(t0,3),sigma.formula=~cs(t0,3),nu.formula=~1,data=size,family=SN1,method=RS(250))
m12<-gamlss(t1~cs(t0,3),sigma.formula=~1,nu.formula=~1,data=size,family=SN1,method=RS(250))

AIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12) # m10 is best-supported
bestskew<-m10

## Get parameters and pdf for a given starting size
# get pdf for starting size = 45
skewmu<-predict(bestskew,newdata=data.frame(t0=0:100),data=size,type='response',what=c('mu'))
skewsigma<-predict(bestskew,newdata=data.frame(t0=0:100),data=size,type='response',what=c('sigma'))
skewnu<-predict(bestskew,newdata=data.frame(t0=0:100),data=size,type='response',what=c('nu'))

# Note: the dSN1 function in the gamlss package can sometimes return nonsensical 
# values for x values far outside the distribution, so using dsn here to avoid this issue
dskew<-dsn(0:100,xi=skewmu[46],omega=skewsigma[46],alpha=skewnu[46]) 
sub<-size[which(size$t0>44&size$t0<46),] # get actual size at t+1 for starting size between 42 and 48

hist(sub$t1,breaks=10,freq=F,xlim=c(newSize[46,1]*0.8,newSize[46,2]*1.1),main='',xlab="Size at time t+1")
lines(0:100,dskew,type='l',lwd=2,col='forestgreen')

## Visualize model predictions
skewcent<-centiles.pred(bestskew,xname='t0',xvalues=0:100,cent=c(1,50,99))

plot(size$t0,size$t1,xlab='Size at time t',ylab='Size at time t+1',col='grey')
lines(0:100,skewcent[,3],col='forestgreen',lwd=2,lty=1)
lines(0:100,skewcent[,2],col='forestgreen',lwd=2,lty=2)
lines(0:100,skewcent[,4],col='forestgreen',lwd=2,lty=2)

## Compare AIC to the exact same data as the beta approach
skewMu<-predict(bestskew,data=size,type='response',what=c('mu'))
skewSigma<-predict(bestskew,data=size,type='response',what=c('sigma'))
skewNu<-predict(bestskew,data=size,type='response',what=c('nu'))
size$t1n<-backbeta(size$t1b,min=minsizes,max=maxsizes) 
LikSkew<-log(dsn(size$t1n,xi=skewMu,omega=skewSigma,alpha=skewNu))
LLSkew<-sum(LikSkew)
KSkew<-bestskew$df.fit
AICSkew<- -2*LLSkew + 2*KSkew

### Other pdfs for starting size x-1
x=11
zb<-seq(0.0001,0.9999,length.out=101)
zn<-backbeta(x=zb,min=newSize[x,1],max=newSize[x,2]) # back-transform to original data scale
db<-dBE(zb,mu=betaMu[x],sigma=betaSigma[x])
db<-db/(zn[101]-zn[1]) # back-transform to original data scale
dsn<-dsn(seq(0,20,length.out=100),xi=skewmu[x],omega=skewsigma[x],alpha=skewnu[x]) 
dn<-dNO(seq(0,20,length.out=100),mu=normMean2[x],sigma=normSD[x]) 
sub<-size[which(size$t0>(x-2)&size$t0<x),] # get actual size at t+1 for starting size between 42 and 48

hist(sub$t1,breaks=20,freq=F,xlim=c(newSize[x,1]*0.5,newSize[x,2]*1.5),main='',xlab="Size at time t+1")
lines(seq(0,20,length.out=100),dsn,type='l',lwd=2,col='forestgreen')
lines(seq(0,20,length.out=100),dn,type='l',lwd=2,col='purple')
lines(zn,db,type='l',lwd=2,col='darkred')

## Compare which data points are predicted well or poorly
## across methods
ggplot(size, aes(x=t0, y=t1, color=LikNorm2-LikBeta)) + geom_point() +
      scale_color_gradient2(midpoint=0, low="red", mid="white", high="blue", space ="Lab" )
ggplot(size, aes(x=t0, y=t1, color=LikNorm2-LikSkew)) + geom_point() +
  scale_color_gradient2(midpoint=0, low="forestgreen", mid="white", high="blue", space ="Lab" )
ggplot(size, aes(x=t0, y=t1, color=LikSkew-LikBeta)) + geom_point() +
  scale_color_gradient2(midpoint=0, low="red", mid="white", high="forestgreen", space ="Lab" )


