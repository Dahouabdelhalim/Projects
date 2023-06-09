###
### Gardner et al. (2021) Integrated animal movement and spatial capture-recapture models: simulation, implementation, and inference
### Appendix: Langevin movement and SCR
###
### Nimble functions, distributions, and samplers for Langevin movement SCR
###


########################################################################
##
## (function)
## Get expected location based on Langevin movement and potential gradients
## Gardner et al. equations 4 and 5
## 
########################################################################

GetLangExpectLoc <- nimbleFunction(
  run = function(s = double(1), offsetX=double(0), offsetY=double(0), 
                 grad.dist=double(0), cov1 =double(2),neighMat=double(2),grad_mat=double(2),CellArea=double(0),
                 delta = double(0), sigmaMOVE=double(0) ){
    returnType(double(1))
    ans<-rep(0,2) #placeholder
    xL<-trunc(((s[1]+abs(offsetX))/grad.dist)) +1         #look-up table row
    yL<-trunc(((s[2]+abs(offsetY))/grad.dist)) +1         #look-up table column
    
    gradID<-grad_mat[xL,yL]    #get the gradient matrix ID for s
 
     # Gardner et al. eq 5
     dfdx <-((cov1[neighMat[gradID,2],2] - s[2])*(cov1[neighMat[gradID,4],3] - cov1[neighMat[gradID,1],3]) + (s[2] - cov1[neighMat[gradID,1],2])*
          (cov1[neighMat[gradID,3],3]- cov1[neighMat[gradID,2],3]))/(CellArea)

     dfdy <-((cov1[neighMat[gradID,4],1] - s[1])*(cov1[neighMat[gradID,2],3] - cov1[neighMat[gradID,1],3]) + (s[1] - cov1[neighMat[gradID,1],1])*
          (cov1[neighMat[gradID,3],3]- cov1[neighMat[gradID,4],3]))/(CellArea)

    # Gardner et al. eq 4 (expected value)
    ans[1]<- s[1] + (sigmaMOVE^2/2)* (delta*dfdx)  # expected x location
    ans[2]<- s[2] + (sigmaMOVE^2/2)* (delta*dfdy)  # expected y location
    return(ans)
  }
)




########################################################################
###  
### (function)
### Get sum of lambdas across all traps and occasions. 
### distance calculations ignored if z==0.
### Gardner et al. equation 1
###  							   			 
########################################################################

GetLambdaStar <- nimbleFunction(
  run = function(s = double(2), lambda0=double(0), sigmaDet=double(0), 
                 X=double(2), J=double(0), TT=double(0), z=double(0)){
    returnType(double(0))
    if(z==0) return(0)
    if(z==1){
     lambda <- matrix(0, J, TT)
     d2 <- matrix(0, J, TT)
     for(tt in 1:TT){
      d2[1:J,tt] <- ((s[1,tt]-X[1:J,1])^2 + (s[2,tt]-X[1:J,2])^2) # Euclidean distance
      lambda[1:J,tt] <- lambda0*exp(-d2[1:J,tt]/(2*sigmaDet^2))  # encounter rate
     }#tt
     return(sum(lambda[1:J, 1:TT]))
   }#if z==1
  }
)




########################################################################
###  
###  (distribution)
###  Vector Poisson
###  Gardner et al. section 2.1 
###  							   			 
########################################################################


deregisterDistributions("dpoisVec")  # forced to reload in case any edits
dpoisVec <- nimbleFunction(
  run = function(	x = double(1),s = double(1), X=double(2),  
			sigmaDet=double(0), lambda0=double(0),
                 	J=double(0), z=double(0),
			log = integer(0, default = 0)) {
    returnType(double(0))

    # if z==0, detection == 0
    if(z==0){ 
 	 if(log) return(0)
  	 else return(1)
    }#z==0

    # if z==1, detection is a function of distance to traps
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     lam <- lambda0*exp(-d2/(2*sigmaDet^2))
	 loglike <- sum(dpois(x[1:J], lam[1:J], log=TRUE)) 
	 if(log) return(loglike)
	 else return(exp(loglike))
    }#z==1
 }
)

#random value generator
rpoisVec <- nimbleFunction(
 run = function(n = integer(0), s = double(1), X=double(2),  
			sigmaDet=double(0), lambda0=double(0),
                 	J=double(0), z=double(0)) {
 returnType(double(1))

    # if z==0, detection == 0. currently only uses z=1 individuals. but keep for flexibility.
    if(z==0) return(rep(0, J))

    # if z==1, detection is a function of distance to traps
    if(z==1){
     d2 <- ((s[1]-X[1:J,1])^2 + (s[2]-X[1:J,2])^2)
     lam <- lambda0*exp(-d2/(2*sigmaDet^2))
     rp <- rpois(J, lam[1:J])
     return(rp[1:J])
    }#z==1

})





########################################################################
###  
###  (distribution)
###  Inhomogeneous point process where points are in continuous space and:
###    1. points restricted to state-space, 
###    2. continuous-space points assigned to correct grid cell
###
###  This is used for locations at occasion 1. Gardner et al. equations 6 and 7
###  							   			 
########################################################################


deregisterDistributions("dPointProcess")  # forced to reload in case any edits
dPointProcess <- nimbleFunction(
  run = function( x = double(1),pi=double(1), nPix=double(0), grid.dist=double(0), 
                  xlim=double(1), ylim=double(1),
                  pixMat=double(2), gcoords=double(2), 
			log = integer(0, default = 0)) {
   returnType(double(0))

    # assign continuous space location to grid cell
     sx <- trunc(((x[1]+abs(xlim[1]))/grid.dist)) +1
     sy <- trunc(((x[2]+abs(ylim[1]))/grid.dist)) +1
     pix <- pixMat[sx, sy]

     loglike <- dcat(pix, pi[1:nPix], log=TRUE)
 
     if(log) return(loglike)
     else return(exp(loglike))
 }
)

#random value generator. 
rPointProcess <- nimbleFunction(
 run = function(n = integer(0), pi=double(1), nPix=double(0), grid.dist=double(0), 
                  xlim=double(1), ylim=double(1),
                  pixMat=double(2), gcoords=double(2) ) {
 returnType(double(1))
    # select grid cell
    pix <- rcat(1, pi[1:nPix])
    # xy coordinates are uniform within grid cell
    x <- c(runif(1, max(xlim[1],gcoords[pix,1]-.5*grid.dist), min(xlim[2], gcoords[pix,1]+.5*grid.dist)), 	# x-location
	       runif(1, max(ylim[1],gcoords[pix,2]-.5*grid.dist), min(ylim[2], gcoords[pix,2]+.5*grid.dist)))	# y-location
    
	return(x[1:2])


})


########################################################################
###  
###  (distribution)
###  BVN random walk truncated at state space 
###  							   			 
########################################################################

deregisterDistributions("dRW")  # forced to reload in case any edits
dRW <- nimbleFunction(
  run = function(	x = double(1),s = double(1), sigmaMOVE=double(0),  
			xlim=double(1), ylim=double(1),
			log = integer(0, default = 0)) {
    returnType(double(0))

     # Truncated normal with the limit set as the boundary of the state space (xlim, ylim).
      loglike <- dnorm(x[1], s[1], sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], s[1], sd=sigmaMOVE) - pnorm(xlim[1], s[1], sd=sigmaMOVE))) + 	# x-direction
                 dnorm(x[2], s[2], sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], s[2], sd=sigmaMOVE) - pnorm(ylim[1], s[2], sd=sigmaMOVE)))		# y-direction

    if(log) return(loglike)
    else return(exp(loglike))
 }
)

rRW <- nimbleFunction(
  run = function(n = integer(0), s = double(1), sigmaMOVE=double(0),  
			xlim=double(1), ylim=double(1)) {
    returnType(double(1))
    ans <- rep(0, 2)

    # nimble functions don't seem to like rtnorm. Write out truncated normal instead.
    tmp.x <- runif(1, pnorm(xlim[1], s[1], sd=sigmaMOVE), pnorm(xlim[2], s[1], sd=sigmaMOVE))
    ans[1] <- qnorm(tmp.x, s[1], sigmaMOVE)

    tmp.y <- runif(1, pnorm(ylim[1], s[2], sd=sigmaMOVE), pnorm(ylim[2], s[2], sd=sigmaMOVE))
    ans[2] <- qnorm(tmp.y, s[2], sigmaMOVE)

    return(ans[1:2])
 }
)



########################################################################
###  
###  (distribtion)
###  Full trajectory of augmented individuals
###  Langevin movement, potential gradients, locations restricted to state space
###  							   			 
########################################################################

deregisterDistributions("dRWTrajectoryLangevin") # forced to reload in case any edits
dRWTrajectoryLangevin <- nimbleFunction(
  run = function(x = double(2), sigmaMOVE=double(0),
                 pi = double(1),nPix=double(0),offsetX=double(0), offsetY=double(0), 
                 xlim=double(1),ylim=double(1),cov1=double(2),neighMat=double(2),
                 grid.dist=double(0), pixMat=double(2), gcoords=double(2),
                 grad.dist=double(0),grad_mat=double(2),
                 delta=double(0), TT=double(0), CellArea=double(0), 
				 log = integer(0, default = 0)) {
    returnType(double(0))

    loglike<-matrix(0, nr=2, nc=TT-1)
    mu<-matrix(0, nr=2, nc=(TT-1))
    xL <- rep(0, (TT-1) )
    yL <- rep(0, (TT-1) )
    gradID<- rep(0, (TT-1) )
    dfdx<- rep(0, (TT-1) )
    dfdy<- rep(0, (TT-1) )

    # initial continuous space location as grid cell number
     sx <- trunc(((x[1,1]+abs(xlim[1]))/grid.dist)) +1
     sy <- trunc(((x[2,1]+abs(ylim[1]))/grid.dist)) +1
     start_pix <- pixMat[sx, sy]
     # loglike
     loglikePix <- dcat(start_pix, pi, log=TRUE)

    # locations for occasions >1

    for(tt in 2:TT){

    xL[tt-1]<-trunc(((x[1,tt-1]+abs(offsetX))/grad.dist)) +1       #look-up table row for gradient
    yL[tt-1]<-trunc(((x[2,tt-1]+abs(offsetY))/grad.dist)) +1       #look-up table column for gradient
  
    gradID[tt-1]<-grad_mat[xL[tt-1],yL[tt-1]]    #get the gradient matrix ID for x
 
     dfdx[tt-1] <-((cov1[neighMat[gradID[tt-1],2],2] - x[2,tt-1])*(cov1[neighMat[gradID[tt-1],4],3] - cov1[neighMat[gradID[tt-1],1],3]) + (x[2,tt-1] - cov1[neighMat[gradID[tt-1],1],2])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],2],3]))/CellArea

     dfdy[tt-1] <-((cov1[neighMat[gradID[tt-1],4],1] - x[1,tt-1])*(cov1[neighMat[gradID[tt-1],2],3] - cov1[neighMat[gradID[tt-1],1],3]) + (x[1,tt-1] - cov1[neighMat[gradID[tt-1],1],1])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],4],3]))/CellArea

      mu[1,tt-1]<- x[1,tt-1] + (sigmaMOVE^2/2) * (dfdx[tt-1]*delta)           #gradient in X
      mu[2,tt-1]<- x[2,tt-1] + (sigmaMOVE^2/2) * (dfdy[tt-1]*delta)           #gradient in Y

      # truncated normal
      loglike[1, tt-1] <- dnorm(x[1,tt], mu[1,tt-1], sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], mu[1,tt-1], sd=sigmaMOVE) - pnorm(xlim[1], mu[1,tt-1], sd=sigmaMOVE)))
      loglike[2, tt-1] <- dnorm(x[2,tt], mu[2,tt-1], sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], mu[2,tt-1], sd=sigmaMOVE) - pnorm(ylim[1], mu[2,tt-1], sd=sigmaMOVE)))
   }#tt

 ans <- loglikePix  + sum(loglike[1:2,1:(TT-1)]) 
 if(log) return(ans)
 else return(exp(ans))

 }
)


#random value generator
rRWTrajectoryLangevin  <- nimbleFunction(
 run = function(n = integer(0), sigmaMOVE=double(0),
                 pi = double(1),nPix=double(0),offsetX=double(0), offsetY=double(0), 
                 xlim=double(1),ylim=double(1),cov1=double(2),neighMat=double(2),
                 grid.dist=double(0), pixMat=double(2), gcoords=double(2),
                 grad.dist=double(0),grad_mat=double(2),
                 delta=double(0), TT=double(0), CellArea=double(0)) {

 returnType(double(2))
    x<-matrix(0, nr=2, nc=TT)
    mu<-matrix(0, nr=2, nc=(TT-1))
    xL <- rep(0, (TT-1) )
    yL <- rep(0, (TT-1) )
    gradID<- rep(0, (TT-1) )
    dfdx<- rep(0, (TT-1) )
    dfdy<- rep(0, (TT-1) )

    # initial location as grid cell number
     # starting cell
      start_pix <- rcat(1, pi[1:nPix])
     # starting coordinates
	  # assume location is unif within selected grid cell
      x[1,1] <- runif(1, gcoords[start_pix,1]-.5*grid.dist, gcoords[start_pix,1]+.5*grid.dist) 	# x-location
	  x[2,1] <- runif(1, gcoords[start_pix,2]-.5*grid.dist, gcoords[start_pix,2]+.5*grid.dist)	# y-location


    # occasions >1
    for(tt in 2:TT){
      xL[tt-1]<-trunc(((x[1,tt-1]+abs(offsetX))/grad.dist)) +1       #look-up table row for gradient
      yL[tt-1]<-trunc(((x[2,tt-1]+abs(offsetY))/grad.dist)) +1       #look-up table column for gradient
  
     gradID[tt-1]<-grad_mat[xL[tt-1],yL[tt-1]]    #get the gradient matrix ID for x
 
     dfdx[tt-1] <-((cov1[neighMat[gradID[tt-1],2],2] - x[2,tt-1])*(cov1[neighMat[gradID[tt-1],4],3] - cov1[neighMat[gradID[tt-1],1],3]) + (x[2,tt-1] - cov1[neighMat[gradID[tt-1],1],2])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],2],3]))/CellArea

     dfdy[tt-1] <-((cov1[neighMat[gradID[tt-1],4],1] - x[1,tt-1])*(cov1[neighMat[gradID[tt-1],2],3] - cov1[neighMat[gradID[tt-1],1],3]) + (x[1,tt-1] - cov1[neighMat[gradID[tt-1],1],1])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],4],3]))/CellArea

      mu[1,tt-1]<- x[1,tt-1] + (sigmaMOVE^2/2) * (dfdx[tt-1]*delta)           #gradient in X
      mu[2,tt-1]<- x[2,tt-1] + (sigmaMOVE^2/2) * (dfdy[tt-1]*delta)           #gradient in Y

      # nimble functions don't seem to like rtnorm. Write out truncated normal instead.
      tmp.x <- runif(1, pnorm(xlim[1], mu[1,tt-1], sd=sigmaMOVE), pnorm(xlim[2], mu[1,tt-1], sd=sigmaMOVE))
      x[1,tt] <- qnorm(tmp.x,  mu[1,tt-1], sigmaMOVE)

      tmp.y <- runif(1, pnorm(ylim[1], mu[2,tt-1], sd=sigmaMOVE), pnorm(ylim[2],  mu[2,tt-1], sd=sigmaMOVE))
      x[2,tt] <- qnorm(tmp.y, mu[2,tt-1], sigmaMOVE)

   }#tt

 return(x)
})




########################################################################
## 	
## 	(sampler)
##  Proposes candidate trajectory using:
##   Langevin movement, potential gradients, locations restricted to state space
## 	 conditional on z[i]                         			   
## 	
########################################################################

myRWtrajectoryLangevinSampler <- nimbleFunction(
 contains = sampler_BASE,
 setup = function(model, mvSaved, target, control) {
   # Defined stuff
   ind<-control$ind
   TT<-control$TT    
   scale <-control$scale
   nPix <- control$nPix
   grad.dist<-model$grad.dist
   offsetX<-model$offsetX
   offsetY<-model$offsetY
   xlim<-model$xlim
   ylim<-model$ylim
   gcoords<-model$gcoords
   pixMat<-model$pixMat
   grid.dist<-model$grid.dist
   grad_mat<-model$grad_mat
   cov1<-model$cov1
   neighMat<-model$neighMat
   CellArea<-model$CellArea
   calcNodes <- model$getDependencies(target)
 },

 run = function() {
   # grab current values
   z.curr <- model$z[ind]

   ## If z==0 generate candidate trajectory that is conditional on mvmt parameters but not current trajectory
   ## Accept with probablity 1.0 (i.e., if z==0 Gibbs update)
   if(z.curr==0){
    sigmaMOVE<-model[["sigmaMOVE"]]
    pi<-model[["pi"]]
    delta<-model[["delta"]]

    s.cand<-matrix(0, nrow=2, ncol=TT)
    mu<-matrix(0, nr=2, nc=(TT-1))
    xL <- rep(0, (TT-1) )
    yL <- rep(0, (TT-1) )
    gradID<- rep(0, (TT-1) )
    dfdx<- rep(0, (TT-1) )
    dfdy<- rep(0, (TT-1) )


     # starting cell
      start_pix.cand <- rcat(1, pi[1:nPix])
     # starting coordinates
       s.cand[1,1] <- runif(1, max(xlim[1],gcoords[start_pix.cand,1]-.5*grid.dist), min(xlim[2], gcoords[start_pix.cand,1]+.5*grid.dist)) 	# x-location
 	   s.cand[2,1] <- runif(1, max(ylim[1],gcoords[start_pix.cand,2]-.5*grid.dist), min(ylim[2], gcoords[start_pix.cand,2]+.5*grid.dist))	# y-location

     for(tt in 2:TT){
      
      xL[tt-1]<-trunc(((s.cand[1,tt-1]+abs(offsetX))/grad.dist)) +1       #look-up table row for gradient
      yL[tt-1]<-trunc(((s.cand[2,tt-1]+abs(offsetY))/grad.dist)) +1       #look-up table column for gradient
  
     gradID[tt-1]<-grad_mat[xL[tt-1],yL[tt-1]]    #get the gradient matrix ID for x
 
     dfdx[tt-1] <-((cov1[neighMat[gradID[tt-1],2],2] - s.cand[2,tt-1])*(cov1[neighMat[gradID[tt-1],4],3] - cov1[neighMat[gradID[tt-1],1],3]) + (s.cand[2,tt-1] - cov1[neighMat[gradID[tt-1],1],2])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],2],3]))/CellArea

     dfdy[tt-1] <-((cov1[neighMat[gradID[tt-1],4],1] - s.cand[1,tt-1])*(cov1[neighMat[gradID[tt-1],2],3] - cov1[neighMat[gradID[tt-1],1],3]) + (s.cand[1,tt-1] - cov1[neighMat[gradID[tt-1],1],1])*
          (cov1[neighMat[gradID[tt-1],3],3]- cov1[neighMat[gradID[tt-1],4],3]))/CellArea

      mu[1,tt-1]<- s.cand[1,tt-1] + (sigmaMOVE^2/2) * (dfdx[tt-1]*delta)           #gradient in X
      mu[2,tt-1]<- s.cand[2,tt-1] + (sigmaMOVE^2/2) * (dfdy[tt-1]*delta)           #gradient in Y

      # Write out truncated normal.
      tmp.x <- runif(1, pnorm(xlim[1], mu[1,tt-1], sd=sigmaMOVE), pnorm(xlim[2], mu[1,tt-1], sd=sigmaMOVE))
      s.cand[1,tt] <- qnorm(tmp.x,  mu[1,tt-1], sigmaMOVE)

      tmp.y <- runif(1, pnorm(ylim[1], mu[2,tt-1], sd=sigmaMOVE), pnorm(ylim[2],  mu[2,tt-1], sd=sigmaMOVE))
      s.cand[2,tt] <- qnorm(tmp.y, mu[2,tt-1], sigmaMOVE)

     }#tt

    # store proposal into model
    model$s[ind, 1:2, 1:TT] <<- s.cand[1:2, 1:TT]

    # update dependencies
    model$calculate(calcNodes)

    # Accept candidate trajectory with probabliy 1.0 and keep the model and mvSaved objects consistent
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)

   }#if z.curr==0

   ## If z==1 generate candidate trajectory conditional on current trajectory and use MH update to accept/reject 
   if(z.curr==1){

     # initial model logProb
     model_lp_initial <- model$getLogProb(calcNodes)

     # generate candidate trajectory conditional on current trajectory
     ## this could be more efficient, including adapting the proposal scale parameter
     ## may not matter much since z[i] regularly switches 
     s.curr<-model$s[ind, 1:2, 1:TT]
     s.cand<-matrix(0, nrow=2, ncol=TT)

     ## for tt in 1:TT
      ## truncated normal
      sz.x <- runif(TT, pnorm(xlim[1], s.curr[1,1:TT], sd=scale), pnorm(xlim[2], s.curr[1,1:TT], sd=scale))
      s.cand[1,1:TT] <- qnorm(sz.x[1:TT], s.curr[1,1:TT], scale)

      sz.y <- runif(TT, pnorm(ylim[1], s.curr[2,1:TT], sd=scale), pnorm(ylim[2], s.curr[2,1:TT], sd=scale))
      s.cand[2,1:TT] <- qnorm(sz.y[1:TT], s.curr[2,1:TT], scale)

     # store proposal into model
     model$s[ind, 1:2, 1:TT] <<- s.cand[1:2, 1:TT]

     # proposal model logProb
     model_lp_proposed <- model$calculate(calcNodes)

     # log-Metropolis-Hastings ratio
     log_MH_ratio <- model_lp_proposed - model_lp_initial
     # Metropolis-Hastings step: determine whether or
     # not to accept the newly proposed value
     accept <- decide(log_MH_ratio)
     if(accept) {
        copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
     } else {
        copy(from = mvSaved, to = model, row = 1, nodes = calcNodes, logProb = TRUE)
     }
   }#z==1

 },
 methods = list( reset = function () {} )
)

## END sampler