###
### Gardner et al. (2022) Integrated animal movement and spatial capture-recapture models: simulation, implementation, and inference
### Appendix: Random Walk movement and SCR
###
###
### This script simulates and analyzes SCR-movement models
### using continuous-space discrete-time random walk movement processes
###
### Current settings reflect Scenario 7a in Gardner et al 2022. (correlated random walk with telemetry)
###
### Reducing to a Bivariate Normal Random walk without directional persistence requires fixing gamma=0 in the simulation and model
###
### MCMC requires SCR_RandomWalks_Distributions_Samplers_and_Functions.R to run.
###  
###
### 
###
### A few modifications from Gardner et al. 2022 indexing
###   t indexed as tt. Time steps.
###   T indexed as TT. Total number of occasions
###   h = 1,2,... telemHits is the number of telemetry hits per occasion.
###

set.seed(1234567)

## required packages
library(nimble)
library(coda)
library(msm)
library(abind)
  
## source file with nimble functions, distributions, and samplers
 source("SCR_RandomWalks_Distributions_Samplers_and_Functions.R")

## Parameters to monitor
params <- c("sigmaMOVE", "gamma", "sigmaDet","lambda0", "N", "psi")
  
  # save s and z if needed. otherwise set to NA
  parameters2 <-NA   #<- c("s", "z")	#		
  nt2 <- 10	         # if saving s and z, set high thinning value for parameters2 to prevent memory issues

## MCMC settings 
# test run
nburnin = 1; niter = 5000+nburnin; nchains = 2; adaptInterval = 200; nthin = 1
#nburnin = 15000; niter = 55000+nburnin; nchains = 3; adaptInterval = 1000; nthin = 2

# custom full trajectory sampler does not currently adapt. Thus a reasonable scaling parameter must be selected for proposed locations  
SamplerScale <- 0.01 	# proposal scale (sd) for custom trajectory sampler. Scale for proposal of new locations for augmented individuals when z==1. 
  

##
##
##
## Section 1. Data generation and formatting
##
##
##
  
## Define study
N <- 100          ## Abundance
M <- 250          ## data augmentation limit
TT <- 25          ## occasions
lambda0 <- 5      ## Local baseline encounter rate. 
sigmaDet <- 0.05  ## Scale parameter of local detection function
sigmaMOVE <- 0.10 ## standard deviation for movement 
gamma <- 0.50     ## directional persistence 0<=gamma<=1. gamma==0 is a BVN random walk. i.e. no persistence
ylim <- xlim <- c(0,10) ## state-space 

# possible telemetry data
nTelem <- 10      ## number of individuals with telemetry. Set to 0 or some positive integer <=N
telemHits <- 2    ## number of telemetry locations per occasion (e.g., GPS collar settings)


## Section 1b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SCR Trap coordinates
x0 <- seq(3, 7, length.out=10)
X <- cbind(rep(x0, each=10),rep(x0, times=10))
J <- nrow(X)    ## number of traps


## Section 1c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate initial locations and subsequent movement (truncated at state-space)
sall<-  array(NA, c(N,2,TT)) # placeholder

## Location at occasion 1
sall[,1:2,1] <- cbind(runif(N, xlim[1], xlim[2]), runif(N, ylim[1], ylim[2]))

 # Location at occasion 2 (Gardner et al. eqs. 2 and 3)
    sall[,1:2,2] <- cbind(rtnorm(N, sall[,1,1], sigmaMOVE, lower=xlim[1], upper=xlim[2]), 
                          rtnorm(N, sall[,2,1], sigmaMOVE, lower=ylim[1], upper=ylim[2]))

  
 # Location at occasion 3:TT (Gardner et al. eqs. 2 and 3)
   for(tt in 3:TT){
    # expected location in x and y
    e.x <- sall[,1,tt-1] + (sall[,1,tt-1]-sall[,1,tt-2])*gamma 
    e.y <- sall[,2,tt-1] + (sall[,2,tt-1]-sall[,2,tt-2])*gamma 

    # realized location is expected location and noise    
    sall[,1:2,tt] <- cbind(rtnorm(N, e.x, sigmaMOVE, lower=xlim[1], upper=xlim[2]), 
                           rtnorm(N, e.y, sigmaMOVE, lower=ylim[1], upper=ylim[2]) )
  }
  

## Section 1d ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SCR observation process (Gardner et al. eq. 1)
d <- lambda <- yall <- array(NA,dim=c(N,J,TT)) # placeholders
for(i in 1:N) {
  for(tt in 1:TT) { 
    d[i,,tt] <- sqrt((sall[i,1,tt] - X[,1])^2 +(sall[i,2,tt] - X[,2])^2)
    lambda[i,,tt] <- lambda0*exp(-d[i,,tt]^2 /(2*sigmaDet^2))
    yall[i,,tt] <- rpois(J, lambda[i,,tt])
  }
}
  
detected <- apply(yall>0, 1, any, na.rm=TRUE)
y <- yall[detected,,,drop=FALSE]
scrID <- which(detected)
nind <- length(scrID)
if(nTelem==0) nTotal <- nind # total number of unique individuals with spatial data (SCR or telemetry)



## Section 1e ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Add in telemetry data if nTelem > 0
## here, any in individual can have telemetry 
## telemetry is active from the start of the study to the end
## ID of telemetry individuals are known in SCR observations 
## Gardner et al. section 2.3
  
if(nTelem >0){
    telemID <- sample(N, nTelem)              # nTelem individuals had telemetry from start of study
    inBoth <- telemID[telemID%in%scrID]       # telem individual observed in y
    telemOnly <- telemID[!(telemID%in%scrID)] # telem individual not observed in y
    nTelemBoth <- length(inBoth)
    nTelemOnly <- length(telemOnly)
    nTotal <- nind + nTelemOnly  # total number of unique individuals with spatial data (SCR or telemetry)

    # match individuals across datasets. This can be bit tricky.
    telemSCRrowID <- NA
    for(i in 1:nTelem){
      # if in both telem and scr
      if(telemID[i]%in%scrID) telemSCRrowID[i] <- which(telemID[i]==scrID)
      # if in telem only
      if(!(telemID[i]%in%scrID) ) telemSCRrowID[i] <- nind + which(telemID[i]==telemOnly)
    }
    
    ## generate location data for telemetered individuals, restricted to state-space 
    ## currently assumes telemetry locations are random about s[i,1:2,tt]       
    muTelem <- array(NA, c(nTelem, 2, TT, telemHits)) # telemetry locations. observed without error.
    for(i in 1:nTelem){
     for(tt in 1:TT){
      muTelem[i,1,tt, 1:telemHits] <- rtnorm(telemHits, sall[telemID[i],1,tt], sd=sigmaDet, lower=xlim[1], upper=xlim[2])	# x-direction
      muTelem[i,2,tt, 1:telemHits] <- rtnorm(telemHits, sall[telemID[i],2,tt], sd=sigmaDet, lower=ylim[1], upper=ylim[2])	# y-direction
     }
    }
 
    # telemOnly have observed (known) all zero encounter histories
    yTotal <-abind(y, array(0, c(length(telemOnly), J, TT)), along=1)

  }#if(nTelem >0)
  


##
## simple plot for example
##
  
  if(1==0){
    #all non-detected trajectories
    matplot(t(sall[!detected,1,]), t(sall[!detected,2,]), pch=16, col=rgb(1,0,0,0.2),
      xlim=xlim, ylim=ylim, xlab="X", ylab="Y") #all non-detected trajectories
    matlines(t(sall[!detected,1,]), t(sall[!detected,2,]), lty=1, col=rgb(1,0,0,0.5))
    points(X[,1], X[,2], pch="x", col="black", cex=1)
    # all detected individuals
    matpoints(t(sall[detected,1,]), t(sall[detected,2,]), pch=16, col=rgb(0,0,1,0.5)) 
    matlines(t(sall[detected,1,]), t(sall[detected,2,]), lty=1, col=rgb(0,0,1,0.5)) 

   # telemetry locations
   if(nTelem>0){
    for(i in 1:nTelem){
     matpoints(t(muTelem[i,1,,]), t(muTelem[i,2,,]), pch=15, col=rgb(0,0,0,0.5), cex=.5) 
    }
   }  
    
  }#if(1==0)
##END PLOT
  
  
  
   
##
##
##
## Section 2. NIMBLE model using data augmentation.
##            For clarity, we provide separate models for when nTelem = 0 of nTelem >0.
##            However, the Telemetry model would work when nTelem = 0.
##
##            Model:
##            1. Data augmentation    
##            2. Continuous-space discrete-time random walk movement
##            3. Custom distribution to evaluate full trajectory of augmented individuals
##            4. Full-trajectory sampler for augmented individuals that proposes full trajectory conditional on z[i]
##
##
 
  
  
 
 
  
## Section 2a ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SCR-movement without telemetry

if(nTelem == 0){ 

## NIMBLE model
code <- nimbleCode({

 ## Process for observed individuals 
 for(i in 1:nTotal){
  z[i] ~ dbern(psi)

  ## location at tt=1
   # starting location in continuous space
   s[i,1,1] ~ dunif(xlim[1],xlim[2])  # x-location
   s[i,2,1] ~ dunif(ylim[1],ylim[2])  # y-location

   # detection process conditional on location s and state z
   y[i,1:J,1]~dpoisVec(s=s[i,1:2,1], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

  ## for tt = 2, movement and detection
   # location. For efficiency, we jointly evaluate the x- and y- locations
   s[i,1:2,2] ~ dRW1(s = s[i,1:2,1], sigmaMOVE = sigmaMOVE, gamma = gamma, xlim=xlim[1:2], ylim=ylim[1:2])

   # detection process conditional on location s and state z
    y[i,1:J,2]~dpoisVec(s=s[i,1:2,2], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

  ## for tt = 3:TT, movement and detection
  for(tt in 3:TT){
   # location. For efficiency, we jointly evaluate the x- and y- locations
   s[i,1:2,tt] ~ dRW2(s = s[i,1:2,(tt-2):(tt-1)], sigmaMOVE = sigmaMOVE, gamma = gamma, xlim=xlim[1:2], ylim=ylim[1:2])

   # detection process conditional on location s and state z
    y[i,1:J,tt]~dpoisVec(s=s[i,1:2,tt], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])
  }# tt
 }#i

 ## Processes for augmented individuals 
 for(i in (nTotal+1):M){
  z[i] ~ dbern(psi) 

  # joint distribution for initial location and movement
  s[i,1:2,1:TT] ~ dRWTrajectory(sigmaMOVE=sigmaMOVE,gamma=gamma,
                                   xlim=xlim[1:2],ylim=ylim[1:2], TT=TT, z=z[i])
        
  # lambdaStar = sum of trap and occasion specific lambdas, given trajectory and z[i]
  lambdaStar[i] <- GetLambdaStar(s=s[i,1:2,1:TT], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, TT=TT, z=z[i]) 
  zeros[i]~dpois(lambdaStar[i]) # augmented individuals were not detected on all occasions.

 }#i


 ## abundance
 N <- sum(z[1:M])

 ## Priors
 psi ~ dbeta(.000001,1)			# scale prior for data augmentation
 sigmaMOVE ~ dunif(0,5)		# movement between occaions
 sigmaDet ~ dunif(0,5) 		# detection (movement within an occasion)
 lambda0 ~ dunif(0,10)		# detection rate at distance zero
 gamma ~ dbeta(1,1)         # directional persistence

})#NimModel

### END model

}#if(nTelem == 0) 




## Section 2b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## SCR-movement with telemetry

if(nTelem >0){ 

## NIMBLE model
code <- nimbleCode({

 ## Process for observed individuals 
 for(i in 1:nTotal){
  z[i] ~ dbern(psi*(1-q[i]) + q[i]) # Gardner et al. eq. 9
   # where q[i] = 1 if individual i has telemetry data and 0 otherwise.

  ## location at tt=1
   # starting location in continuous space
   s[i,1,1] ~ dunif(xlim[1],xlim[2])  # x-location
   s[i,2,1] ~ dunif(ylim[1],ylim[2])  # y-location

   # detection process conditional on location s and state z
   y[i,1:J,1]~dpoisVec(s=s[i,1:2,1], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

  ## for tt = 2, movement and detection
   # location. For efficiency, we jointly evaluate the x- and y- locations
   s[i,1:2,2] ~ dRW1(s = s[i,1:2,1], sigmaMOVE = sigmaMOVE, gamma = gamma, xlim=xlim[1:2], ylim=ylim[1:2])

   # detection process conditional on location s and state z
    y[i,1:J,2]~dpoisVec(s=s[i,1:2,2], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

  ## for tt = 3:TT, movement and detection
  for(tt in 3:TT){
   # location. For efficiency, we jointly evaluate the x- and y- locations
   s[i,1:2,tt] ~ dRW2(s = s[i,1:2,(tt-2):(tt-1)], sigmaMOVE = sigmaMOVE, gamma = gamma, xlim=xlim[1:2], ylim=ylim[1:2])

   # detection process conditional on location s and state z
    y[i,1:J,tt]~dpoisVec(s=s[i,1:2,tt], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])
  }# tt
 }#i


 ## Processes for augmented individuals 
 for(i in (nTotal+1):M){
  z[i] ~ dbern(psi*(1-q[i]) + q[i]) # Gardner et al. eq. 9
   # where q[i] = 1 if individual i has telemetry data and 0 otherwise.

  # joint distribution for initial location and movement
  s[i,1:2,1:TT] ~ dRWTrajectory(sigmaMOVE=sigmaMOVE,gamma=gamma,
                                     xlim=xlim[1:2],ylim=ylim[1:2], TT=TT, z=z[i])
        
  # lambdaStar = sum of trap and occasion specific lambdas, given trajectory and z[i]
  lambdaStar[i] <- GetLambdaStar(s=s[i,1:2,1:TT], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, TT=TT, z=z[i]) 
  zeros[i]~dpois(lambdaStar[i]) # augmented individuals were not detected on all occasions.

 }#i


 ## Additional data from telemetry locations. See Gardner et al. section 2.3 
  for(i in 1:nTelem){         # total number of telemetry individuals
   for(tt in 1:TT){           # assumes telemetry available for all occasions
    for(h in 1:telemHits){    # assumes constant number of locations per occasion
     # telemSCRrowID is given as data to link i in nTelem to correct row in SCR data
      muTelem[i,1,tt,h] ~ T(dnorm(s[telemSCRrowID[i],1,tt], sd=sigmaDet), xlim[1], xlim[2]) # x-location
      muTelem[i,2,tt,h] ~ T(dnorm(s[telemSCRrowID[i],2,tt], sd=sigmaDet), ylim[1], ylim[2]) # x-location
    }#h
   }#tt
  }#nTelem

 ## abundance
 N <- sum(z[1:M])

 ## Priors
 psi ~ dbeta(.000001,1)		# scale prior for data augmentation
 sigmaMOVE ~ dunif(0,5)		# movement between occaions
 sigmaDet ~ dunif(0,5) 		# detection (movement within an occasion)
 lambda0 ~ dunif(0,10)		# detection rate at distance zero
 gamma ~ dbeta(1,1)       # directional persistence


})#NimModel

### END model

}#if(nTelem > 0) 

  
## Section 2c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### organize data
 # augmentation
 nz <- M-nTotal
 zeros <- rep(0, M)
 zAug <- c(rep(1, nTotal), rep(NA, nz) )  # z is known for observed individuals (both trapped and telemetry)

 ## stuff to nimble if no telemetry data
 if(nTelem==0){
    constants<-list(nTotal=nTotal, M=M, TT=TT, X=X, J=J) 
                   
    nim.data<-list(z=zAug, y = y, zeros=zeros, 
                   xlim = xlim, ylim = ylim)  
    
 }#if(nTelem==0)



 ## stuff to nimble if telemetry data
 if(nTelem>0){
  # indicator for telemetry individuals so that they do not contirubute to z
  qAug <- rep(0, M)
  qAug[telemSCRrowID] <- 1

  constants<-list(nTotal=nTotal, M=M, TT=TT, X=X, J=J, 
                  nTelem=nTelem, telemSCRrowID=telemSCRrowID, telemHits=telemHits) 
    
  nim.data<-list(z=zAug, q=qAug, y = yTotal, zeros=zeros, 
                 muTelem=muTelem,
                 xlim = xlim, ylim = ylim)  
 
 }#if(nTelem>0)



 # Inits function
 # simplistic location inits 
 init.s <- array(NA, c(M, 2, TT) )
  for(i in 1:nind){
    occdet <- which(apply(y[i,,], 2, sum)>0)
    for(tt in 1:TT){
      init.s[i, 1:2, tt] <-  c(mean( X[which(y[i,,occdet[which.min(abs(occdet  - tt))]]>0),1] ),
                               mean( X[which(y[i,,occdet[which.min(abs(occdet  - tt))]]>0),2] ))
    }#tt
  }#nind
  # for augmented individuals. 
  for(i in (nTotal+1):M){
    init.s[i, 1:2, 1:TT] <-  c(runif(1, xlim[1],xlim[2]), runif(1, ylim[1],ylim[2]) )
  }
  if(nTelem>0){    # for telemetry individuals pick a location near mean (if exactly at mean if has difficulty moving)
    for(i in 1:nTelem){
     for(tt in 1:TT){
       init.s[telemSCRrowID[i], 1:2, tt] <-  c(rtnorm(1, muTelem[i,1,tt,],sd=.1, lower=xlim[1], upper=xlim[2] ),
                                               rtnorm(1, muTelem[i,2,tt,],sd=.1, lower=ylim[1], upper=ylim[2] ))
    }#tt
   }#nTelem
  }#if(nTelem>0)   
    
 ## Here we simply set values close to true values. 
 inits <- list(lambda0=rnorm(1, lambda0,sd=.05),  sigmaMOVE = runif(1, 0.05, 0.15), sigmaDet=runif(1, 0.05, 0.15),  
               s = init.s, z = c(rep(NA, nTotal), rep(0, nz)),
               gamma=runif(1, .4, .6) )
    
 ## Compile and run using NIMBLE
 rm(Rmodel,Cmodel,Rmcmc,Cmcmc) # just to be safe
 Rmodel <- nimbleModel(code=code, constants=constants, data=nim.data,check=FALSE, calculate=F, inits=inits)
 conf <- configureMCMC(Rmodel,monitors=params, control = list(adaptInterval = adaptInterval), thin=nthin) # start with small scale and adapt from there
  ## change samplers here
  #conf$printSamplers()
    
  ## for observed individuals, jointly propose s[i, 1:2, tt]
  conf$removeSampler(paste0("s[",1:nTotal,", 1:2, 1:",TT,"]"))
  for(i in 1:nTotal){
    for(tt in 1:TT){
        conf$addSampler(target = paste0("s[",i,", 1:2, ",tt,"]"),
                        type = 'RW_block', silent=TRUE,
                        control = list(adaptInterval = adaptInterval, propCov='identity', scale=SamplerScale, adaptScaleOnly=TRUE))
    }
  }
  ## for augmented individuals, use custom full trajectory proposal
  conf$removeSampler(paste0("s[",(nTotal+1),":",M,", 1:2, 1:",TT,"]"))
  for(i in (nTotal+1):M){
    conf$addSampler(target = paste0("s[",i,", 1:2, 1:",TT,"]"),
                    type = 'myRWtrajectorySampler',
                    control = list(TT=TT, ind=i, scale=SamplerScale)) # set scale parameter for candidate locations if z==1
  }
    
  # setting smaller initial scale parameters then adapating reduces chances of getting stuck at incorrect local maximums due to poor initial values
  conf$removeSampler(c("gamma", "sigmaMOVE","sigmaDet","lambda0"))
  conf$addSampler(target = c("sigmaMOVE"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
  conf$addSampler(target = c("gamma"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
  conf$addSampler(target = c("sigmaDet"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
  conf$addSampler(target = c("lambda0"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
    
    
  ## 'sprinkle' sigmaMOVE sampler in a few times since it is the worst performing parameter
  for(extra in 1:5){
     conf$addSampler(target = c("sigmaMOVE"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
     conf$addSampler(target = c("gamma"), type="RW", control = list(adaptInterval = adaptInterval, scale=.01), thin=nthin)
  }
   nsamplers <- length(conf$getSamplers())
   SampBreaks <- round(seq(1,nsamplers-(5*2), length = (5*2)+2))[2:((5*2)+1)]   # roughly evenly spaced
   SampOrder <- order(rank( c( 1:(nsamplers-(5*2)), SampBreaks), ties.method = "first")) # re-ordered rank
      
   conf$setSamplerExecutionOrder(SampOrder)
   #conf$printSamplers(executionOrder = TRUE)
    
   ## be sure to thin if saving s and/or z 
   if(!is.na(parameters2[1])) {
     conf$addMonitors2(parameters2)
     conf$setThin2(nt2)
   }
    
 Rmcmc <- buildMCMC(conf)  
 Cmodel <- compileNimble(Rmodel)
 Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
    
 ##Run multiple chains
 out <- runMCMC(Cmcmc, niter = niter , nburnin = nburnin , nchains = nchains, inits=inits,
                setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  

    
 summary(out)
 plot(out[,"N"])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### END 
  
