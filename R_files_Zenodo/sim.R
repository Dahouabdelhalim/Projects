#########################################################################
###########   HIERARCHICAL DISEASE MODEL - TWO PATHOGENS   ##############
#########################################################################

library(deSolve) # for simulating ODEs

#########################################################################
#########################   ANALYTICAL SUB-MODELS   #####################
#########################################################################

local.sim <- function(t,y,params){
  R <- y[1]; P1 <- y[2]; P2 <- y[3]
  with(as.list(params), { 
    dR <- S*(1-(R/Rmax)) - (u1*R*P1*Q1) - (u2*R*P2*Q2)
    dP1 <- (P1*u1*R) - (m1*P1)
    dP2 <- (P2*u2*R) - (m2*P2)
    res <- c(dR, dP1, dP2)
    list(res)})}

# R: resource (concentration in host)
# Rmax: maximum homeostatic balance of resources in host
# S: resource supply
# P1: pathogen 1 abundance in host
# P2: pathogen 2 abundance in host
# Q's: resource quota of pathogens
# u's: resource-dependent birth rates of pathogens
# d's: death rates of pathogens


#########################################################################
###################    STOCHASTIC HIERARCHICAL MODEL    #################
#########################################################################

disease.metacomm.tr.freq <- function(n, t.ext, t.int, m, u1, u2, d, S, R.max, Q, beta, v1, v2, P1.start, P2.start){
  
  # 1) initialize. 
  hosts <- seq(1:n); cycles <- seq(1:t.ext) # how many hosts? how many transmission cycles?
  meta.com.sum <- data.frame(cycle=cycles) # create meta.com.sum for end of each cycle & count
  R <- rep(R.max, n); P1 <- rbinom(n, size=1, prob=P1.start); P2 <- rbinom(n, size=1, prob=P2.start) # starting conditions
  t <- seq(from=0, to=t.int, by=1) # number time steps to simulate local dynamics 
  params <- c(S=S, Q1=Q, Q2=Q, u1=u1, u2=u2, m1=m, m2=m, Rmax=R.max) # parameters for within-host dynamics
  last.host <- matrix(ncol=4) # random last host; will keep track of all internal dynamics to plot
  R.star.1 <- m/u1; R.star.2 <- m/u2 # R* for pathogens 1 & 2
  slope1 <- (v2-v1)/(R.star.2-R.star.1); int1 <- d+v2-(slope1*R.star.2) # slope & intercept for steep segment of resource-death function
  slope2 <- -v2/(R.max-R.star.2); int2 <- d-(slope2*R.max) # slope & intercept for shallow segment of resource-death function
  meta.com.sum <- data.frame() # empty data frame to fill in with summary statistics from metacommunity 
  
  # 2) loop through transmission cycles
  for(j in 1:length(cycles)){ # start to loop through cycles
    print(j) # keep track of where we are in simulation
    
    # 2A) simulate the within-host dynamics
    comp.mod <- function(R, P1, P2){
      sim <- lsoda(y=c(R=R, P1=P1, P2=P2), times=t, func=local.sim, parms=params, rtol=1, atol=1)} # simulate within hosts
    mat.sum <- t(matrix(mapply(comp.mod, R=R, P1=P1, P2=P2), ncol=n)) # vectorize across all hosts
    R <- mat.sum[, 2*(t.int+1)]; P1 <- mat.sum[, 3*(t.int+1)]; P2 <- mat.sum[, 4*(t.int+1)] # pull out results 
    last.host <- rbind(last.host, matrix(mat.sum[n,], ncol=4)) # and the random last host
    P1.hosts <- which(P1>1); P2.hosts <- which(P2>1) # define infections as hosts where pathogen abundance >1
    P12.hosts <- which(P1>1 & P2>1); P0.hosts <- which(P1<1 & P2<1) 
    sumlist <- list(P1.prev=length(P1.hosts)/n, P1.intens=mean(P1[P1.hosts]), P1.R=mean(R[P1.hosts]),
                    P2.prev=length(P2.hosts)/n, P2.intens=mean(P2[P2.hosts]), P2.R=mean(R[P2.hosts]),
                    P12.prev=length(P12.hosts)/n, co.P1.intens=mean(P1[P12.hosts]), co.P2.intens=mean(P2[P12.hosts]), 
                    P0.prev=length(P0.hosts)/n, mean.R=mean(R), co.R=mean(R[P12.hosts])) # calculate infection prevalence, etc.
    
    # 2B) stochastic death and transmission. 
    d.prob <- ifelse(R>=R.star.2, 1-exp(-(slope2*R + int2)), 1-exp(-(slope1*R + int1))) # list of n; depends on R in each host
    died <- rbinom(n=n, size=1, prob=d.prob) # binomial death probabilities
    trP1s <- sum(rbinom(n=length(which(P1>1)), size=1, prob=1 - exp(-beta * (1-sumlist[["P1.prev"]])))) # transmission probabilities
    trP2s <- sum(rbinom(n=length(which(P2>1)), size=1, prob=1 - exp(-beta * (1-sumlist[["P2.prev"]])))) # transmission probabilities
    P1 <- ifelse(died==1, 0, P1); P2 <- ifelse(died==1, 0, P2); R <- ifelse(died==1, R.max, R) # kill (reset) hosts that died
    if(trP1s<length(which(P1<1))){newP1 <- sample(which(P1<1), size=trP1s, replace=F)} else {newP1 <- which(P1<1)} # pick hosts to get infected
    if(trP2s<length(which(P2<1))){newP2 <- sample(which(P2<1), size=trP2s, replace=F)} else {newP2 <- which(P2<1)} # pick hosts to get infected
    P1[newP1] <- 1; P2[newP2] <- 1 # distribute new infections
    
    # 2C) save meta.com.sum statistics for each step of competition, death and transmission
    meta.com.sum <- rbind(meta.com.sum, data.frame(sumlist)) 
  } # close the cycle loop
  
  meta.com.sum$cycle <- cycles; rownames(meta.com.sum) <- cycles
  meta.com.sum <<- meta.com.sum # save to global environment
  
  # 3) plot the simulation. 
  par(mfrow=c(2,1),mar=c(1,2.5,0,0), oma=c(0,0,0.5,0.5)) 
  
  colnames(last.host) <- c("time", "R", "P1", "P2")
  last.host <- as.data.frame(last.host)
  last.host <- last.host[last.host$time!=t.int , ]
  last.host$time2=seq(1:length(last.host$time))
  plot(last.host$time2, last.host$R, type="l", col="dark green", ylim=c(0,22), lwd=5,
       xaxt='n', yaxt='n', ann=FALSE, lty=3)
  lines(last.host$time2, last.host$P2, col="blue", lty=2, lwd=3)
  lines(last.host$time2, last.host$P1, col="red", lwd=3)
  axis(side=2, at=c(0,10,20))
  
  plot(meta.com.sum$cycle, meta.com.sum$P1.prev, ylim=c(0,1), col="red", xlab='n', type='l', lwd=3,
       xaxt='n', yaxt='n', ann=FALSE)
  lines(meta.com.sum$cycle, meta.com.sum$P2.prev, col="blue", lwd=3, lty=2)
  axis(side=2, at=c(0 ,0.5, 1))
  axis(side=1, at=c(1, t.ext/2, t.ext), labels=F)
  
  
} # end function


#########################################################################
########################   EXAMPLE SIMULATION   #########################
#########################################################################

disease.metacomm.tr.freq(n=100, S=5, R.max=10, m=0.1, u1=0.1, u2=0.02, Q=2, t.int=10, t.ext=20,
                 beta=0.5, d=0.0, v1=0.3, v2=0.01, P1.start=0.5, P2.start=0.5)

