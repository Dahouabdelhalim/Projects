## This file contains example code for estimating age-at-maturity and other life history
## parameters from length-at-age using a Lester biphasic growth model. The model is fit
## using a profile likelihood approach. This code corresponds with the Ecological
## Applications article by Honsey et al. entitled "Accurate estimates of age-at-maturity 
## from the growth trajectories of fishes and other ectotherms".

# Load required library
library(boot) # use install.packages('boot') if package isn't already installed

## Step 1: Read in your data and name it 'data'.


## Step 2: Store age and length data as vectors named 'age' and 'len', respectively (see examples)
age<-data$Age ## example of storing age data as vector named 'age'
len<-data$TotalLength ## example of storing length data as vector named 'len'


## Step 3: Store Lester model likelihood function as 'Biphas.Lik.MA' excluding age-at-maturity parameter 
## Optional: include marginal likelihoods for immature growth slope and intercept
###### likelihood function
Biphas.Lik.MA = function(parms) { 
  b0 = parms[1]
  h1 = parms[2]
  g = inv.logit(parms[3])   
  sighat = sqrt(parms[4]) 
  
  age.i = age[age<=mat.age]
  len.i = len[age<=mat.age]
  age.m = age[age>mat.age]
  len.m = len[age>mat.age]
  
  t1 = -b0/h1                                          
  Linf = 3*h1/g
  k = log(1 + g/3)
  t0 = mat.age + 
    suppressWarnings(log(1-(g*(mat.age-t1)/3)))/log(1+g/3)                                  
  mn.i = b0 + h1*age.i
  mn.m = Linf*(1-exp(-k*(age.m-t0)))
  
  b0.lik = dnorm(b0,mean=b0est,sd=25,log=T) #optional (also, distribution can be adjusted if needed)
  h1.lik = dnorm(h1,mean=h1est,sd=5, log=T) #optional (also, distribution can be adjusted if needed)
  L.i = dnorm(len.i,mean=mn.i,sd=sighat,log=T)
  L.m = dnorm(len.m,mean=mn.m,sd=sighat,log=T)
  return(sum(c(L.i,L.m, b0.lik,h1.lik)))
}


## OPTIONAL Step 4: Use linear model of first few years to inform marginal likelihoods
## for early growth slope & intercept
immdata<- data[ which(age <= (min(age)+3)), ] #choose data within first four ages -- number of ages can be changed
immout<-lm(TotalLength~Age, data=immdata) #linear regression on "immature" data -- change formulation as needed to match your data column names
b0est<-immout$coefficients[[1]] # store intercept estimate, used for prior likelihood
h1est<-immout$coefficients[[2]] # store slope estimate, used for prior likelihood


## Step 5: List starting values for each parameter
b0 = b0est # early growth intercept (if you skipped Step 4, you should put a number here)
h1 = h1est # early growth slope (if you skipped Step 4, you should put a number here)
g = 0.15  # cost to somatic growth of maturity
sighat = 25 # standard deviation
parms=c(b0,h1,logit(g),sighat^2)  #compile parameters



## Step 6: Create a vector of potential values for age-at-maturity
## and a matrix for storing parameter estimates
Mat.age = seq(min(age),max(age),by=0.025) #  range of mat.age values for profile likelihood calculation -- adjust as needed
lik<-b0<-h1<-g<-var<-rep(NA,length(Mat.age)) # create empty vectors for parameters
mat.age.Lik = cbind(Mat.age,lik,b0,h1,g,var) # create matrix for storing parameter estimates



## Step 7: Optimize likelihood function for each potential age-at-maturity value
## and store parameter estimates
for(j in 1:length(Mat.age)) {
  mat.age = Mat.age[j] # fix age-at-maturity at a given value
  L.out = try(optim(par=parms,fn=Biphas.Lik.MA, 
                    control=list(fnscale=-1,reltol=1e-8)), silent=T) # optimize likelihood function
  check<-is.numeric(L.out[[1]]) # check to see if model converged
  
  ## store values only if model worked
  if (check[[1]] == "TRUE"){
    
    #Store parameter values (back-transform g)
    mat.age.Lik[j,2] <- L.out$value
    mat.age.Lik[j,3] <- L.out$par[[1]]
    mat.age.Lik[j,4] <- L.out$par[[2]]
    mat.age.Lik[j,5] <- inv.logit(L.out$par[[3]])
    mat.age.Lik[j,6] <- L.out$par[[4]]
  }}


## Step 8: Analyze results
mat.age.Lik<-as.data.frame(mat.age.Lik) #convert to data frame for easier referencing
mat.age.Lik<-mat.age.Lik[which(mat.age.Lik$lik != "NA"),] #remove failed runs
mle = max(mat.age.Lik$lik) ## find maximum likelihood
MLE<-mat.age.Lik[which(mat.age.Lik$lik == mle),] ## maximum likelihood estimates for all parameters
MLE ## print maximum likelihood estimates


## Additional steps: Plotting likelihood ratio profile, calculating confidence intervals,
## and plotting curves onto data

## Plot likelihood profile
rlike = exp(mat.age.Lik$lik-mle)
plot(mat.age.Lik$Mat.age,rlike,xlab="Age-at-Maturity", 
     ylab="Likelihood Ratio")


## Confidence interval in terms of chi-squared (~ 95% CI)
ndx1 = which(mat.age.Lik$lik>(mle-1.92)) # change '1.92' to 0.228 for 50% CI, 1.36 for 90% CI
points(mat.age.Lik$Mat.age[ndx1],rep(0,length(ndx1)),col='red',lwd=6)
CI = c(min(mat.age.Lik$Mat.age[ndx1]),max(mat.age.Lik$Mat.age[ndx1]))
CI # print confidence interval

## Plot biphasic growth curves onto data
plot(age,len,xlim=c(0,max(age)),ylim=c(0,max(len)))  
g. = MLE[[5]]
h1. = MLE[[4]]
mT = MLE[[1]]
b0. = MLE[[3]]
t1. = -b0./h1.
Linf = 3*h1./g.
k = log(1 + g./3)
t0 = mT + log(1-(g.*(mT-t1.)/3))/log(1+g./3)
abline(b0.,h1.,lwd=3)  
matX = seq(mT,max(age),length.out=25)
matY = Linf*(1-exp(-k*(matX-t0)))
lines(matX,matY,col='red',type='l',lwd=6,lty=1)