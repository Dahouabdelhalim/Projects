###
### Gardner et al. (2021) Integrated animal movement and spatial capture-recapture models: simulation, implementation, and inference
### Appendix: Random walk movement models and SCR
###
### Nimble functions, distributions, and samplers for Random Walk movement SCR
###




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
      d2[1:J,tt] <- ((s[1,tt]-X[1:J,1])^2 + (s[2,tt]-X[1:J,2])^2) # euclidian distance
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
###  Random walk truncated at state space using only previous location
###  Gardner et al. equation 2							   			 
########################################################################

deregisterDistributions("dRW1")  # forced to reload in case any edits
dRW1 <- nimbleFunction(
  run = function(	x = double(1),s = double(1), 
                  sigmaMOVE=double(0), gamma=double(0), 
			xlim=double(1), ylim=double(1),
			log = integer(0, default = 0)) {
    returnType(double(0))

     # Truncated normal in x- and y- directions
      loglike <- dnorm(x[1], s[1], sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], s[1], sd=sigmaMOVE) - pnorm(xlim[1], s[1], sd=sigmaMOVE))) + 	# x-direction
                 dnorm(x[2], s[2], sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], s[2], sd=sigmaMOVE) - pnorm(ylim[1], s[2], sd=sigmaMOVE)))	# y-direction

    if(log) return(loglike)
    else return(exp(loglike))
 }
)

rRW1 <- nimbleFunction(
  run = function(n = integer(0), s = double(1), 
                 sigmaMOVE=double(0), gamma=double(0), 
                 xlim=double(1), ylim=double(1)
                ) {

    returnType(double(1))
    ans <- rep(0, 2)

     # Truncated normal in x- and y- directions
      tmp.x <- runif(1, pnorm(xlim[1], s[1], sd=sigmaMOVE), pnorm(xlim[2], s[1], sd=sigmaMOVE))
       ans[1] <- qnorm(tmp.x, s[1], sigmaMOVE)
      tmp.y <- runif(1, pnorm(ylim[1], s[2], sd=sigmaMOVE), pnorm(ylim[2], s[2], sd=sigmaMOVE))
       ans[2] <- qnorm(tmp.y, s[2], sigmaMOVE)


   return(ans[1:2])
 }
)


########################################################################
###  
###  (distribution)
###  Random walk truncated at state space using two previous locations
###  Gardner et al. equation 3							   			 
########################################################################

deregisterDistributions("dRW2")  # forced to reload in case any edits
dRW2 <- nimbleFunction(
  run = function(	x = double(1),s = double(2), 
                  sigmaMOVE=double(0), gamma=double(0), 
			xlim=double(1), ylim=double(1),
			log = integer(0, default = 0)) {
    returnType(double(0))

     # expected location
      e.xy <- c(s[1,2] + gamma*(s[1,2] - s[1,1]) , s[2,2] + gamma*(s[2,2] - s[2,1]) )

     # Truncated normal in x- and y- directions
	  loglike <- dnorm(x[1], e.xy[1], sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], e.xy[1], sd=sigmaMOVE) - pnorm(xlim[1], e.xy[1], sd=sigmaMOVE))) +
                   dnorm(x[2], e.xy[2], sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], e.xy[2], sd=sigmaMOVE) - pnorm(ylim[1], e.xy[2], sd=sigmaMOVE)))

    if(log) return(loglike)
    else return(exp(loglike))
 }
)

rRW2 <- nimbleFunction(
  run = function(n = integer(0), s = double(2), 
                 sigmaMOVE=double(0), gamma=double(0), 
                 xlim=double(1), ylim=double(1)
                ) {

    returnType(double(1))
    ans <- rep(0, 2)

     # expected location
     e.xy <- c(s[1,2] + gamma*(s[1,2] - s[1,1]) , s[2,2] + gamma*(s[2,2] - s[2,1]) )

     # Truncated normal in x- and y- directions
     tmp.x <- runif(1, pnorm(xlim[1], e.xy[1], sd=sigmaMOVE), pnorm(xlim[2], e.xy[1], sd=sigmaMOVE))
     ans[1] <- qnorm(tmp.x, e.xy[1], sigmaMOVE)

     tmp.y <- runif(1, pnorm(ylim[1], e.xy[2], sd=sigmaMOVE), pnorm(ylim[2], e.xy[2], sd=sigmaMOVE))
     ans[2] <- qnorm(tmp.y, e.xy[2], sigmaMOVE)

   return(ans[1:2])
 }
)





########################################################################
###  
### (distribution)
###  Custom distribution for full Random Walk trajectory of augmented individuals 
###  truncated at state-space							   			 
###  
########################################################################

deregisterDistributions("dRWTrajectory")  # just reloading after edits
dRWTrajectory <- nimbleFunction(
  run = function(x = double(2), sigmaMOVE=double(0), gamma=double(0),
                 xlim=double(1),ylim=double(1), TT=double(0), z=double(0),
                 log = integer(0, default = 0)) {

    returnType(double(0))
    loglike<-matrix(0, nr=2, nc=TT)

    # initial location
    loglike[1,1] <- dunif(x[1,1], xlim[1], xlim[2], log=TRUE)
    loglike[2,1] <- dunif(x[2,1], ylim[1], ylim[2], log=TRUE)

    # movement and locations truncated at state-space
    for(tt in 2:2){
     # Truncated normal in x- and y- directions
     loglike[1,tt] <- dnorm(x[1,tt], x[1,tt-1], sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], x[1,tt-1], sd=sigmaMOVE) - pnorm(xlim[1], x[1,tt-1], sd=sigmaMOVE)))
     loglike[2,tt] <- dnorm(x[2,tt], x[2,tt-1], sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], x[2,tt-1], sd=sigmaMOVE) - pnorm(ylim[1], x[2,tt-1], sd=sigmaMOVE)))
    }#k

    # movement and associated locations truncated at state-space
    for(tt in 3:TT){
     e.x <- x[1,tt-1] + gamma*(x[1,tt-1] - x[1,tt-2]) 
     e.y <- x[2,tt-1] + gamma*(x[2,tt-1] - x[2,tt-2])  

     # Truncated normal in x- and y- directions
     loglike[1,tt] <- dnorm(x[1,tt], e.x, sd=sigmaMOVE, log=TRUE)- log((pnorm(xlim[2], e.x, sd=sigmaMOVE) - pnorm(xlim[1], e.x, sd=sigmaMOVE)))
     loglike[2,tt] <- dnorm(x[2,tt], e.y, sd=sigmaMOVE, log=TRUE)- log((pnorm(ylim[2], e.y, sd=sigmaMOVE) - pnorm(ylim[1], e.y, sd=sigmaMOVE)))
    }#k

  ans <- sum(loglike[1:2,1:TT])
  if(log) return(ans)
  else return(exp(ans))
 }
)



#random value generator
rRWTrajectory <- nimbleFunction(
 run = function(n = integer(0), sigmaMOVE=double(0), gamma=double(0),
                 xlim=double(1),ylim=double(1), TT=double(0), z=double(0)) {
 returnType(double(2))
    x<-matrix(0, nr=2, nc=TT)

    # initial location
    x[1,1]<-runif(1, xlim[1], xlim[2])
    x[2,1]<-runif(1, ylim[1], ylim[2])

    # movement
    for(tt in 2:2){
     # Truncated normal
     tmp.x <- runif(1, pnorm(xlim[1], x[1,tt-1], sd=sigmaMOVE), pnorm(xlim[2], x[1,tt-1], sd=sigmaMOVE))
      x[1,tt] <- qnorm(tmp.x, x[1,tt-1], sigmaMOVE)
     tmp.y <- runif(1, pnorm(ylim[1], x[2,tt-1], sd=sigmaMOVE), pnorm(ylim[2], x[2,tt-1], sd=sigmaMOVE))
      x[2,tt] <- qnorm(tmp.y, x[2,tt-1], sigmaMOVE)
    }#tt

    for(tt in 3:TT){
     e.x <- x[1,tt-1] + gamma*(x[1,tt-1] - x[1,tt-2])  
     e.y <- x[2,tt-1] + gamma*(x[2,tt-1] - x[2,tt-2])  

     # Truncated normal
     tmp.x <- runif(1, pnorm(xlim[1], e.x, sd=sigmaMOVE), pnorm(xlim[2], e.x, sd=sigmaMOVE))
      x[1,tt] <- qnorm(tmp.x, e.x, sigmaMOVE)
     tmp.y <- runif(1, pnorm(ylim[1], e.y, sd=sigmaMOVE), pnorm(ylim[2], e.y, sd=sigmaMOVE))
      x[2,tt] <- qnorm(tmp.y, e.y, sigmaMOVE)
    }#tt

 return(x)
})






########################################################################
##
##   (sampler)			   
##   Proposes candidate trajectory using:
##   Random Walk where locations are restricted to state space
##   conditional on z[i]     
##
########################################################################

myRWtrajectorySampler <- nimbleFunction(
 contains = sampler_BASE,
 setup = function(model, mvSaved, target, control) {
   # Defined stuff
   ind<-control$ind
   TT<-control$TT    
   scale <-control$scale
   xlim<-model$xlim
   ylim<-model$ylim
   calcNodes <- model$getDependencies(target)
 },

 run = function() {
   # grab current values
   z.curr <- model$z[ind]

   ## If z==0 generate candidate trajectory that is not conditional on current trajectory and accept with probablity 1.0
   if(z.curr==0){
    sigmaMOVE<-model[["sigmaMOVE"]]
    gamma<-model[["gamma"]]


    s.cand<-matrix(0, nrow=2, ncol=TT) # placeholder

     s.cand[1,1]<-runif(1, xlim[1], xlim[2])
     s.cand[2,1]<-runif(1, ylim[1], ylim[2])

    # movement
    for(tt in 2:2){
     # Truncated normal
     tmp.x <- runif(1, pnorm(xlim[1], s.cand[1,tt-1], sd=sigmaMOVE), pnorm(xlim[2], s.cand[1,tt-1], sd=sigmaMOVE))
      s.cand[1,tt] <- qnorm(tmp.x, s.cand[1,tt-1], sigmaMOVE)
     tmp.y <- runif(1, pnorm(ylim[1], s.cand[2,tt-1], sd=sigmaMOVE), pnorm(ylim[2], s.cand[2,tt-1], sd=sigmaMOVE))
      s.cand[2,tt] <- qnorm(tmp.y, s.cand[2,tt-1], sigmaMOVE)
    }#tt

    for(tt in 3:TT){
     e.x <- s.cand[1,tt-1] + gamma*(s.cand[1,tt-1] - s.cand[1,tt-2])  
     e.y <- s.cand[2,tt-1] + gamma*(s.cand[2,tt-1] - s.cand[2,tt-2])  

     # Truncated normal
     tmp.x <- runif(1, pnorm(xlim[1], e.x, sd=sigmaMOVE), pnorm(xlim[2], e.x, sd=sigmaMOVE))
      s.cand[1,tt] <- qnorm(tmp.x, e.x, sigmaMOVE)
     tmp.y <- runif(1, pnorm(ylim[1], e.y, sd=sigmaMOVE), pnorm(ylim[2], e.y, sd=sigmaMOVE))
      s.cand[2,tt] <- qnorm(tmp.y, e.y, sigmaMOVE)
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

     ## Truncated normal to keep proposal location in state-space
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










