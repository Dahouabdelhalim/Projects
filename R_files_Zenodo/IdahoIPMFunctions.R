########################################################################
# Parameters and functions for single-species Idaho IPMs.
# Code is patterned on size-quality code in Ch. 6 of The Book 
#
# NOTATION: x is size, W is competition pressure.
# At the moment, PSSP is the only 100% complete parameter set. 
# As noted in the paper, the size-scaling of fecundity estimated
# for PSSP is applied to ARTR, with checks that the results are
# robust to changes in this assumption. 
#
# BIG WARNING: what is called "W" in the IPM kernel is actually   
# log(lambda+W), with lambda= 1% of a seedling, so that lambda+W > 0. 
# Back-transformation to "real W" is used in all the 
# demographic rate functions fitted to data. 
#
# Last update to original script July 5 2016, SPE 
# 
# Extended to include Group effects by SPE, Feb 2020 
#######################################################################

################################################################## 
#  Student t-distribution, for use in demog functions
################################################################## 

# Student-t distribution as parameterized in BUGS/JAGS 
jagsdt <- function(x,tau,df){
       fac1 = gamma( (df+1)/2 )/gamma(df/2);
       fac2 = sqrt(tau/(df*pi)); 
       fac3 = (1 + {tau*x^2}/df)^{-(df+1)/2}
       return(fac1*fac2*fac3)
}

# Student-t distribution as parameterized in gamlss 
dT <- function(x,mu,sigma,nu) {
    a=1/2; b=nu/2; B = gamma(a)*gamma(b)/gamma(a+b); 
    fac1 = 1/(sigma*B*sqrt(nu)); 
    fac2 = 1 + ((x-mu)^2)/(nu*sigma^2);
    pwr = -(1+nu)/2; 
    return(fac1*(fac2^pwr))    
}

#######################################################################
# Parameters. The code here creates a data frame "params" with 4 rows,
# one for each species. Parameters are named, and accessed by name in 
# the demog functions. To make this work, extract a row and give it
# the names from the data frame, as in:
# pars <- params[4,]; names(pars) <- names(params); 
#######################################################################

################## Survival Parameters 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival parameters a0, a1, b for NON-seedlings. 
# logit(s) = a0 + a1*x + b*W 
# Source: Fit_Survival_Older.R 

survPars = read.csv("SurvPars.csv") 
parNames=names(survPars)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Survival parameters a0, bW for SEEDLINGS 
# logit(s) = a0 + b*W 
# Source: Fit_Survival_Seedlings.R 

survSeedPars = read.csv("SeedlingSurvPars.csv");
parNames=c(parNames,c("a0Sseed","bWSseed"));   

############### Growth Parameters 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Growth parameters a0 a1 b tau0 c df, for NON-seedlings 
# Mean size = a0 + a1*x + b*W 
# Residuals are t(df) with precision tau0*exp(c*z) 
# using the JAGS parameterization of t distribution.  
# OK for ARTR and PSSP, others still need refitting with new
# classification of seedlings versus older individuals 
# Source: JAGS scripts in SingleSppIPM/IPM_v3 

growPars = read.csv("growPars.csv"); 
parNames=c(parNames,names(growPars))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters for Gaussian mixture model of size for age-2 individuals 
# Parameters f, mu, sigma: fraction too small to measure, and mean,sd  
# of larger individuals in the second component of the mixture 
# Source: Fit_Growth_Seedlings.R 

age2SizePars = read.csv("age2SizePars.csv"); 
parNames = c(parNames,names(age2SizePars))

############### Recruitment Parameters 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Seed production as a function of parent size 
# Mean seed production C*(area^pwr) = C*exp(pwr*x)  
# To make fecundity=lifespan, we give each plant 1 offspring
# per year, guaranteed, regardless of age or size 

fecPars = matrix(c(1, 0, 1, 0, 1, 0, 1, 0), 4,2,byrow=TRUE); 
fecPars <- data.frame(fecPars); 
names(fecPars) <- c("Cseed","pwrseed"); 
parNames = c(parNames,names(fecPars));

############### Competition Parameters  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mean and StdDev of initial log(W+w0) for seedlings (Age = 1)
# Note, the lognormal model is only good for PSSP. 
# We'll need something more general for other species. 
# Source: SeedlingW.R 
wRecPars = matrix(c( 
-2, 2, 
NA, NA,
NA, NA,
-1.14, 1.03),
4,2,byrow=TRUE);   
parNames=c(parNames,c("muW","sdW"));   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# AR(1) coefficients for log-transformed W
# Parameters a0 a1 tau0 c df W0 
# Source: ModelingW.R 

wPars = read.csv("wPars.csv");  
parNames = c(parNames,names(wPars))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Ceiling on growth function to eliminate eviction out the top 
# Floor at the bottom is set by minimum size, log(0.25) 
Epars = matrix(c(8,NA,NA,6.0),4,1); Epars=data.frame(Epars); 
parNames=c(parNames,"Ulim"); 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Floor and ceiling on transitions for W, to eliminate eviction 
WEpars = matrix(c(
-6, 2,
NA, NA,
NA, NA,
-5, 2),
4, 2.5,byrow=TRUE); 
parNames=c(parNames,c("wLlim","wUlim")); 


################## Quadrat Group Parameters #############################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Group effect parameters, seedling survival 
SSG = matrix(c(0.07354396, -0.23961535,  0.22278200,  0.29098337, -0.05889084, -0.28880314),1,6);
SeedSurvGparams = rbind(SSG, matrix(0,3,6));   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
# Group effect parameters, adult survival 
AdultSurvGparams = matrix(c(-0.323598954,  0.163705362, -0.109165424, -0.718369219,  0.989199680, -0.001771444,
0,0,0,0,0,0,
0,0,0,0,0,0,
-0.005277883,  0.074223526, -0.132938796, -0.164246152,  0.000389793,  0.227849513),
4,6,byrow=TRUE);  


#################################################################
# Put it all together 
params=cbind(survPars,survSeedPars,growPars,age2SizePars,fecPars,wRecPars,wPars,Epars,WEpars); 
names(params) <- parNames; 

# [1] "a0S"     "a1S"     "bWS"     "a0Sseed" "bWSseed" "a0"      "a1"      "bW"      "tau"    
#[10] "c"       "df"      "f"       "mu2"     "sd2"     "Cseed"   "pwrseed" "muW"     "sdW"    
#[19] "a0W"     "a1W"     "tau0W"   "cW"      "dfW"     "W0"      "Ulim"    "wLlim"   "wUlim"   
 
################################################################
# Demographic rate functions 
################################################################   

########################### Survival  ##########################
#~~~~~~~~~~~~~ Non-seedlings 
s_x <- function(x,W,pars,Group) {
	with(pars,{
    realW = exp(W)-W0; 
	u = exp(a0S + a1S*x + bWS*realW + AGpars[Group]); 
	return(u/(1+u))
	}) 
} 

#~~~~~~~~~~~~~~~~Seedlings 
s_seed <- function(W,pars,Group){ 
	with(pars,{
    realW = exp(W)-W0; 
	u = exp(a0Sseed + bWSseed*realW + SGpars[Group]); 
	return(u/(1+u))
	}) 
} 

########################### Growth ##########################
# ~~~~~~~~~~~~~~~~~~ Non-seedlings
g_x1x <- function(x1,x,W,pars) {
    with(pars,{
    Wlim = pmax(wLlim,pmin(W,wUlim)); #floor and ceiling
    realW = exp(Wlim)-W0; 
    xlim=pmax(-1.38,pmin(x,Ulim));  # floor and ceiling 
    mu <- a0 + a1*xlim + bW*realW 
    tau <- tau*exp(c*xlim); 
    return(jagsdt(x1-mu,tau,df))
    }) 
}    

# ~~~~~~~~~~~~~~~ Seedlings growing to age 2: Gaussian mixture model 
c2_x1 <- function(x1,pars) {
    with(pars,{
    return(f*dnorm(x1,mean=-1.38,sd=0.2)+(1-f)*dnorm(x1,mean=mu2,sd=sd2))
    }) 
}   

########################### Recruitment #########################
#~~~~~~~~~~~~~~~~ Mean seed production 
b_x <- function(x,pars) {
  with(pars,{
    return(Cseed*exp(x*pwrseed))
  } )
}

#~~~~~~~~~~~~~~~~ Variance in seed production: var=mean for poisson
sigma2.b_x  <- function(x,pars) {
	with(pars, {
    return(Cseed*exp(x*pwrseed))
    }) 
}

########################### Competition #########################
#~~~~~~~~~~~~~~~~~~~~~ Non-seedlings 
# log(W0 + W) is AR(1) with t-distributed innovations and 
# nonconstant variance. The IPM kernel is done in terms of 
# log(W0 + W), so this function computes the AR(1) density
# WITHOUT any transformation to/from log scale ("realW").    
g_w1w <- function(W1,W,pars) {
	with(pars,{
    Wlim = pmax(wLlim,pmin(W,wUlim)); #floor and ceiling
  	mu <- a0W + a1W*Wlim 
    sigma <- tau0W + cW*Wlim   
    return(dT(W1,mu=mu,sigma=sigma,nu=dfW))
    }) 
        
}

#~~~~~~~~~~~~~~~~~~~~~ Seedlings: initial competition  
# Normal distribution for log(W0 + W). The IPM is done in
# terms of log(W0+W) so this function does not use any
# transformations to/from log scale ("realW"). 
c1_w1 = function(w1,pars) {
    with(pars,dnorm(w1,mean=muW,sd=sdW))
}

#################################################################
# Functions to make IPM kernels by midpoint rule
# Variable names and code structure are patterned 
# after the size-quality code in sec. 6.6.2 of The Book  
# using the 'expand.grid' method. 
# RECALL: what's called W in the IPM kernels is log(W0  + real W) 
#################################################################

########### Transition density P(x', W', x, W) for non-seedlings 
p_x1W1 <- function(x1,W1,x,W,pars,Group) {
    return( s_x(x,W,pars,Group)*g_x1x(x1,x,W,pars)*g_w1w(W1,W,pars) ) 
}


mk_K <- function(mx, mW, Lx, Ux, LW, UW, pars,Group) {
	hx <- (Ux - Lx)/mx; yx <- Lx + ((1:mx) - 1/2) * hx
	hW <- (UW - LW)/mW; yW <- LW + ((1:mW) - 1/2) * hW
	
    # make P2, the iteration matrix for non-seedlings 
	Kd<- expand.grid(xp=yx,Wp=yW,x=yx,W=yW)
	P2 <- with(Kd, p_x1W1(xp,Wp,x,W,pars,Group))
    dim(P2) <- c(mx*mW,mx*mW); # collapse it straight to 2D 
    
    # make P1, survival and growth of seedlings
    P1 <- matrix(0,mx*mW,mW); 
    Pd <- expand.grid(xp=yx,Wp=yW);
    for(j in 1:mW) {
        P1[,j] <- with(Pd, g_w1w(Wp,yW[j],pars)*c2_x1(xp,pars)*s_seed(yW[j],pars,Group)) 
    } 
    
    # make F2, fecundity of adults 
    Fmat <- matrix(0,mW,mx); 
    W.init <- c1_w1(yW,pars);
    for(j in 1:mx) Fmat[,j]=W.init*b_x(yx[j],pars) 
    F2 = Fmat; 
    for(k in 2:mW) F2=cbind(F2,Fmat)    
     
    
	return(list(P2 = hx*hW*P2, P1 = hx*hW*P1, F2 = hW*F2, yx=yx, yW=yW, hx=hx, hW=hW) )
}	


#####################################################
# Utility functions  
#####################################################   

# size distribution for ages 2 and above 
sizeDist=function(nt,mx,mW) {
    nmat=matrix(nt,mx,mW)
    return(apply(nmat,1,sum))
}

# W distribution for ages 2 and above 
# as usual this is log(W0+W), the IPM's scale for W
wDist=function(nt,mx,mW) {
    nmat=matrix(nt,mx,mW)
    return(apply(nmat,2,sum))
}

# Population mean and variance of a set of values in vector vals, 
# having probabilities in the vector probs. 
pmean=function(vals,probs) {
    sum(vals*probs)
}

pvar=function(vals,probs) {
    mu=pmean(vals,probs)
    pmean( (vals-mu)^2, probs)
}    
    
        	
