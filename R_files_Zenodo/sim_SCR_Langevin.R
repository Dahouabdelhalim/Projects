###
### Gardner et al. (2022) Integrated animal movement and spatial capture-recapture models: simulation, implementation, and inference
### Appendix: Langevin movement and SCR
###
###
### This script simulates and analyzes SCR-movement models
### using a Langevin diffusion process
###
### Current settings reflect Scenario 13 in Gardner et al. 2022.
###
### Requires SCR_Langevin_Mvmt_Distributions_Samplers_and_Functions.R to run.
###  
### A few modifications from Gardner et al. 2022 indexing
###   G indexed as nPix for clarity. Number of grid cells
###   t indexed as tt. Time steps.
###   T indexed as TT. Total number of occasions


# required packages
library(raster)
library(ctmcmove)
library(msm)
library(nimble)
library(coda)
library(reshape2)
library(fields)

source("SCR_Langevin_Mvmt_Distributions_Samplers_and_Functions.R") # should throw "Cannot deregister ..." which is okay
set.seed(123456)  	# required to generate identical raster for all scenarios 

## Parameters to monitor
params <- c("sigmaMOVE", "sigmaDet", "lambda0", "delta", "N", "psi")

## MCMC settings 
# test run
nburnin = 1; niter = 1000+nburnin; nchains = 2; adaptInterval = 200; nthin = 1
#settings used in the paper
#nburnin = 15000; niter = 55000+nburnin; nchains = 3; adaptInterval = 1000; nthin = 2

SamplerScale <- 0.05 	# proposal scale (sd) for custom trajectory sampler. Scale for proposal of new locations if z==1. 

##
##
##
## Section 1. Data generation and formatting
##
##
##

## Define study
N <- 100                   	## Abundance
M <- 250					          ## data augmentation limit
TT <- 25					          ## sampling occasions
ylim <- xlim <- c(0,10)  	  ## state-space. Verify this matches extent of raster used in analysis. 
delta <- 0.7                ## influence of covariate on initial distribution and movement. 
sigmaMOVE <- 0.10           ## movement (BVN RW)
sigmaDet <- 0.10            ## Scale parameter of local detection function
lambda0 <- 2.0              ## Local baseline encounter rate. Lots of photos if near camera


## Section 1b ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gradient calculation uses neighbors. 
# Thus we buffer the state-space to derive gradients, then restrict analyses to the state-space. 
# be careful with set-up. We want state-space to be exactly xlim and ylim.
# easiest to have each grid cell fall completely in the state-space rather than adjust for grid cell area.
ylimB <- xlimB <- c(-1,11) 	  ## buffered state-space for potential function.  
B <- 60                       ## buffered number of rows and columns
pixLen <- diff(xlimB)/B       ## length of grid cell edge
co <- seq(xlimB[1]+pixLen/2, xlimB[2]-(pixLen/2), length=B) ## centroids
Z <- cbind(rep(co, each=B), rep(co, times=B)) ## placeholder


## Function to generate spatial covariate (modified from Andy Royle's code)
spcov <- function(R,alpha) {
    v <- sqrt(nrow(R))
    D <- as.matrix(dist(R))
    V <- exp(-D/alpha)
    cov1 <- t(chol(V)) %*% rnorm(nrow(R))
    Rd <- as.data.frame(R)
    colnames(Rd) <- c("x", "y")
    Rd$C <- as.numeric((cov1 - mean(cov1)) / sd(cov1))
    return(Rd)
}

# create raster covariate
cov1 <- spcov(Z, alpha=1)
r1 <- rasterFromXYZ(cov1, crs=CRS("+proj=utm"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset to state-space
Scoords <- as(extent(xlim[1], xlim[2], ylim[1], ylim[2]), 'SpatialPolygons')
r <- crop(r1, Scoords )


#make small grid to do gradient calculations
smallgrid<-seq(-.85, 10.85, by=.1)
gradgrid<-expand.grid(smallgrid, smallgrid)
gradgrid<-gradgrid[,2:1]


##create a neighborhood matrix
distmat<-rdist(gradgrid, cov1[,1:2])
neighMat<-matrix(NA, dim(gradgrid)[1], 4)
for(i in 1:dim(gradgrid)[1]){
neighMat[i,]<-which(distmat[i,] %in% sort(distmat[i,])[1:4])  # this gives you an index vector
}

neighMat<-neighMat[,c(1,2,4,3)]   #so the values go lower left, upper left, lower right, upper right

#make a raster with the "vals" as the grid ID to relate back to the neighborhood matrix
gradgrid$ID<-1:dim(gradgrid)[1]
colnames(gradgrid)<-c("x", "y", "vals")
rastergradgrid <- rasterFromXYZ(gradgrid, crs=CRS("+proj=utm"))

## grab original covariate values
habitat <- as.data.frame(r, xy=TRUE) # let's call it habitat for simplicity. 

## 
gcoords <- as.matrix(habitat[,c("x", "y")])
nPix <- nrow(gcoords)
CellArea <- prod(res(r))    # constant in this simulation 
grid.dist <- sqrt(CellArea) # distance between centroids



## Section 1c ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Generate starting locations (RSF) and movement

sall <-  array(NA, c(N,2,TT))
pi <- exp(delta*habitat$C)/sum(exp(delta*habitat$C))    ## expected distribution prob of N individuals at occasion 1
dfdx<-NULL
dfdy<-NULL
## Starting location. Select grid cell, then unif() within grid cell
 # starting grid cell
 start_pix <- rep(1:nPix,rmultinom(1, N, pi))   
 #sall[,1:2,1] <- gcoords[start_pix,] # if simply assigning centroid.
 # uniform location within the grid cell
 sall[,1:2,1] <- cbind( runif(N, gcoords[start_pix,1]-.5*grid.dist, gcoords[start_pix,1]+.5*grid.dist), 
                        runif(N, gcoords[start_pix,2]-.5*grid.dist, gcoords[start_pix,2]+.5*grid.dist) )

## subsequent movement
for(tt in 2:TT){
  # grab small grid ID associated with location at tt-1
   gradID <-extract(rastergradgrid, sall[,1:2,tt-1])

  for(ii in 1:N){
   dfdx[ii] <-((cov1[neighMat[gradID[ii],2],2] - sall[ii,2,tt-1])*(cov1[neighMat[gradID[ii],4],3] - cov1[neighMat[gradID[ii],1],3]) + (sall[ii,2,tt-1] - cov1[neighMat[gradID[ii],1],2])*
              (cov1[neighMat[gradID[ii],3],3]- cov1[neighMat[gradID[ii],2],3]))/CellArea

   dfdy[ii] <-((cov1[neighMat[gradID[ii],4],1] - sall[ii,1,tt-1])*(cov1[neighMat[gradID[ii],2],3] - cov1[neighMat[gradID[ii],1],3]) + (sall[ii,1,tt-1] - cov1[neighMat[gradID[ii],1],1])*
              (cov1[neighMat[gradID[ii],3],3]- cov1[neighMat[gradID[ii],4],3]))/CellArea
   } 

     mu_x <-  sall[,1,tt-1] + (sigmaMOVE^2/2)*(delta*dfdx)  # expected x location
     mu_y <-  sall[,2,tt-1] + (sigmaMOVE^2/2)*(delta*dfdy)  # expected y location
  # realized locations
  sall[,1:2,tt] <- cbind(rtnorm(N,mu_x, sd=sigmaMOVE, lower=xlim[1], upper=xlim[2]),  # x loc
                         rtnorm(N,mu_y, sd=sigmaMOVE, lower=ylim[1], upper=ylim[2]) ) # y loc
}#tt


## SCR capture-recapture process

## Trap coordinates
x0 <- seq(3, 7, length.out=10)
X <- cbind(rep(x0, each=10),
           rep(x0, times=10))
J <- nrow(X)    ## number of traps

## Capture-recapture process
d <- lambda <- yall <- array(NA,dim=c(N,J,TT))
for(i in 1:N) {
  for(tt in 1:TT) { 
    d[i,,tt] <- sqrt((sall[i,1,tt] - X[,1])^2 +
                    (sall[i,2,tt] - X[,2])^2)
    lambda[i,,tt] <- lambda0*exp(-d[i,,tt]^2 /(2*sigmaDet^2))
    yall[i,,tt] <- rpois(J, lambda[i,,tt])
  }
}

detected <- apply(yall>0, 1, any, na.rm=TRUE)
y <- yall[detected,,,drop=FALSE]
scrID <- which(detected)
(nind <- length(scrID))



##
##
##
## Section 2. NIMBLE model using data augmentation  
##            1. Habitat influenced BVN movement in continuous space 
##            2. custom distribution to keep starting locations in continuous space
##            3. custom distribution to evaluate full trajectory of augmented individuals
##            4. Full-trajectory sampler for augmented individuals that proposes full trajectory conditional on z[i]
##
##

## Nimble Model
code <- nimbleCode({

## Process for observed individuals 
for(i in 1:nind){
 z[i] ~ dbern(psi)

 # location at time tt=1
 # starting location in continuous space
 s[i, 1:2, 1] ~ dPointProcess(pi=pi[1:nPix], nPix=nPix, grid.dist = grid.dist, 
                              xlim=xlim[1:2], ylim=ylim[1:2],
                              pixMat = pixMat[1:nPix_x, 1:nPix_y],
                              gcoords=gcoords[1:nPix, 1:2])
   # We do not use dcat() as large nPix causes model to be extremely slow
   # the common approach is
    #start_pix[i] ~ dcat(pi[1:nPix])
    #s[i,1:2,1] <- gcoords[start_pix[i],1:2]  

 #detection process conditional on location s and state z
 y[i,1:J,1]~dpoisVec(s=s[i,1:2,1], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

 #for tt>1, movement and detection
 for(tt in 2:TT){
  #movement process as a function of habitat gradients
  # get expected location via nimble function
   mu.xy[i,1:2,tt-1] <- GetLangExpectLoc(s=s[i,1:2,tt-1], offsetX = offsetX, offsetY = offsetY, grad.dist = grad.dist,
                             cov1 = cov1[1:3600,1:3], neighMat = neighMat[1:13924,1:4], grad_mat=grad_mat[1:118,1:118],
                             CellArea=CellArea,
                             delta = delta, sigmaMOVE=sigmaMOVE)

   # BVN RW but expected location is mu.xy and locations restricted to state-space
   s[i,1:2,tt] ~ dRW(s = mu.xy[i,1:2,tt-1], sigmaMOVE = sigmaMOVE, xlim=xlim[1:2], ylim=ylim[1:2])

   #detection process conditional on location s and state z
   y[i,1:J,tt]~dpoisVec(s=s[i,1:2,tt], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, z=z[i])

 }#tt
}#i


## Processes for augmented individuals 
 for(i in (nind+1):M){
  z[i] ~ dbern(psi) 

  # joint distribution for initial location (at centroids) and movement (continuous space)
  s[i,1:2,1:TT] ~ dRWTrajectoryLangevin(sigmaMOVE=sigmaMOVE,
                 xlim=xlim[1:2],ylim=ylim[1:2],
                 pi=pi[1:nPix],nPix=nPix,offsetX = offsetX , offsetY = offsetX,  grad.dist = grad.dist,
                 cov1 = cov1[1:3600,1:3], neighMat = neighMat[1:13924,1:4], grad_mat=grad_mat[1:118,1:118],
                 grid.dist = grid.dist,gcoords=gcoords[1:nPix, 1:2],
                 pixMat = pixMat[1:nPix_x, 1:nPix_y],
                 delta = delta, 
                 TT=TT, CellArea=CellArea)

  # lambdaStar = sum of trap and occasion specific lambdas, given trajectory and z[i]
  lambdaStar[i] <- GetLambdaStar(s=s[i,1:2,1:TT], lambda0=lambda0, sigmaDet=sigmaDet, X=X[1:J, 1:2], J=J, TT=TT, z=z[i]) 
  zeros[i]~dpois(lambdaStar[i])

 }#nind


## abundance
N <- sum(z[1:M])

## initial distribution
pi[1:nPix] <- exp(delta*habitat[1:nPix])/sum(exp(delta*habitat[1:nPix]))

## Priors
 psi ~ dbeta(.000001,1)			#scale prior for data augmentation
 sigmaMOVE ~ dunif(0,5)		# movement between occasions
 sigmaDet ~ dunif(0,5) 		# detection (movement within an occasion)
 lambda0 ~ dunif(0,10)		# detection rate at distance zero
 delta ~ dnorm(0, sd=10)      # resource selection (initial and movement)

})#NimModel

### END model


### organize data

# number of rows and columns, square grid for simplicity 
nPix_x  <- nrow(r)
nPix_y  <- ncol(r)

pix_grad <- nrow(rastergradgrid)
# format gradient covariates as look-up style matrix
grad_mat <- as.matrix(dcast(gradgrid, x~y, fun.aggregate = sum, value.var ="vals")[,1:pix_grad+1])

# assign grid cell number in same manner
pixMat <- as.matrix(dcast(cbind(habitat, 1:nPix), x~y, fun.aggregate = sum, value.var ="1:nPix")[,1:nPix_x+1])

# augmentation
nz <- M-nind
sAug <- array(NA, c(M, 2, TT)) # latent locations
zeros <- rep(0, M)
zAug <- c(rep(1, nind), rep(NA, nz) )

##nimble requires constants, data, intial values 
constants<-list(nind=nind, M=M, TT=TT, X=X, J=J, habitat=habitat$C, nPix=nPix, nPix_x=nPix_x, nPix_y=nPix_y) 

nim.data<-list(z=zAug, y = y, s = sAug, zeros=zeros, 
		       xlim = xlim, ylim = ylim, offsetY=-.9, offsetX=-.9,grad.dist=.1,
               grid.dist=grid.dist, pixMat=pixMat, gcoords=gcoords, 
               grad_mat=grad_mat, cov1=cov1, neighMat=neighMat, CellArea=CellArea)  

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

 for(i in (nind+1):M){
    init.s[i, 1:2, 1:TT] <-  c(runif(1, xlim[1],xlim[2]), runif(1, ylim[1],ylim[2]) )
 }


## Here we simply set close to true values. 
inits <- list(lambda0=rnorm(1, lambda0,sd=.1),  sigmaMOVE = runif(1, 0.05, 0.15), sigmaDet=runif(1, 0.05, 0.15),  
             s = init.s, z = c(rep(NA, nind), rep(0, nz)),
             delta=rnorm(1, delta,sd=.1)) 


# Compile and run using NIMBLE
rm(Rmodel,conf,Cmodel,Rmcmc,Cmcmc) # just to be safe
Rmodel <- nimbleModel(code=code, constants=constants, data=nim.data,check=T, calculate=T, inits=inits)

conf <- configureMCMC(Rmodel,monitors=params, control = list(adaptInterval = adaptInterval, scale=.1), thin=nthin) # start with small scale and adapt from there
 ## change samplers here
 ## for observed individuals, jointly propose s[i, 1:2, tt] 
 conf$removeSampler(paste0("s[",1:nind,", 1:2, 1:",TT,"]"))
 for(i in 1:nind){
  for(tt in 1:TT){
   conf$addSampler(target = paste0("s[",i,", 1:2, ",tt,"]"),
                  type = 'RW_block', silent=TRUE, 
                  control = list(adaptInterval = adaptInterval, propCov=diag(2), scale=SamplerScale, adaptScaleOnly=TRUE)) 
  }
 }

 ## use custom full trajectory proposal for augmented individuals 
 conf$removeSampler(paste0("s[",(nind+1),":",M,", 1:2, 1:",TT,"]"))
 for(i in (nind+1):M){
  conf$addSampler(target = paste0("s[",i,", 1:2, 1:",TT,"]"),
                  type = 'myRWtrajectoryLangevinSampler',
                  control = list(TT=TT, ind=i, scale=SamplerScale, nPix=nPix)) # set scale parameter for candidate locations if z==1
 }

Rmcmc <- buildMCMC(conf)   
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel, showCompilerOutput = FALSE)

##Run 
 out <- runMCMC(Cmcmc, niter = niter , nburnin = nburnin , nchains = nchains, inits=inits,
  setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  

summary(out)
plot(out[,"N"])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### END