# load libraries
library(R2OpenBUGS)
library(coda)
library(data.table)
library(matlab)


mymodel <- function(){
  for(i in 1:N) # rooms x floors loop
  {
    for (t in 1:Tmax){ # weeks loop
      # ZIP model â€“ not executed
      #w[i, t] ~ dbern(psi) 
      #Y[i, t] ~ dpois(eff.mu[i, t]) 
      #eff.mu[i, t]<-w[i, t]*mu[i, t] 
      
      Y[i, t] ~ dpois(mu[i, t]) # mean prevalence for week t in room i of floor j
      
      log(mu[i, t]) <- log(PAR[i, t]) # offset
      +beta[floor[i]]   # intercept per floor k
      +U[i] # unstructured spatial random effect, irrespective of floor
      +V[floor[i]] # structured spatial random effect on floor k
      +theta[t] # structured temporal effect RW(1) for *all* floors
      #+delta[t] # unstructured temporal effect
      
      #risk ratio
      RR[i, t] <- exp(#beta0
        beta[floor[i]]
        +U[i] 
        +V[i] 
        +theta[t]
        #+delta[t]
      )
    } 
    U[i] ~ dnorm(0, tau.U) # unstructured spatial random effect for each room
  } 
  
  ## PRIORS
  
  # SPATIAL
  for(i in 1:Nfloor)
  {
    V[(((i-1)*22)+1):(22*i)] ~ car.normal(adj[]
                                          , weights[] # Adj: ID numbers of the rooms which are adjacent to room i
                                          #	CHECK: A room cannot be specified as its own neighbor
                                          # CHECK symmetric so w_ij=w_ji
                                          , num[] # vector of length N, number of neighbors of room i, i.e. count the polygons in adj for each room
                                          , tauomega.V # precision parameter of the Gaussian prior
    )
  }
  
  for(k in 1:sumNumNeigh) {
    weights[k] <- 1
  }
  
  #------------------------------
  # NEW: prior for floor random intercept: mean = 0
  for(j in 1:Nfloor) {
    beta[j]~dnorm(beta0, prec.tau2) 
  }
  
  prec.tau2~dgamma(0.001, 0.001) # for floors loop
  
  #---------------------
  # U+V
  # weakly informative prior on total tolerance T = 1 / (total variance U + V)
  tau.T ~ dgamma(0.5, 0.005) # tolerance of the total random effects variance U + V
  
  # prior to make sure U+V=T
  p ~ dbeta(1,1) 
  sigma.Z <- sqrt(p/tau.T) 
  omega.V <- sigma.Z
  sigma.U <- sqrt((1-p)/tau.T)
  
  # definitions for tolerance
  tauomega.V <- 1/(omega.V*omega.V) # in ICAR: tolerance for the structured part of the variance: square of the estimate of sd(omega)
  tau.U <- 1/(sigma.U*sigma.U) 
  
  # TEMPORAL
  delta[1]<-0 # to start RW
  theta[1]<-0 # always starts at 0, which isn't totally correct but works
  
  for(t in 2:Tmax){
    delta[t]~dnorm(0, tau.delta)
    theta[t]~dnorm(theta[t-1], tau.theta) # recursive definition of RW - why we need theta[1] defined above
  }
  
  #--------------------------------------------------------
  
  #delta0 ~ dnorm(0, 0.001)
  tau.delta~dgamma(0.1, 0.1)
  tau.theta~dgamma(0.1, 0.1)
  #tau.theta1~dgamma(0.1, 0.1)
  
  
  # FIXED EFFECTS
  beta1 ~ dnorm(0, 0.0001)
  
  # Functions of interest
  #sigma.V <- sqrt(1/tau.V) # standard deviation of non-spatial
  #RRl95 <- exp(-1.96*sigma.V)
  #RRu95 <- exp(1.96*sigma.V) 
  # V and U variance components
  mean.V<-mean(V[1:N])
  sd.V <- sd(V[1:N])       
  mean.U<-mean(U[1:N])
  sd.U<-sd(U[1:N])
  vratio <- sd.V*sd.V/(sd.V*sd.V+sigma.U*sigma.U) 
  mean.theta<-mean(theta[1:Tmax])
  sd.theta<-sd(theta[1:Tmax])
  mean.delta<-mean(delta[1:Tmax])
  sd.delta<-mean(delta[1:Tmax])
  #mean.eps<-mean(eps[1:N, 1:Tmax])
  #sd.eps<-sd(eps[1:N, 1:Tmax])
}


# data
mydata <-list(Tmax = nweeks # weeks
              ,N   = nroom*nfloor
              ,Nfloor=nfloor
              ,Nroom=nroom
              ,Y   = structure(.Data=    Ylong,      .Dim=c(nroom*nfloor, nweeks)) 
              ,PAR = structure(.Data=    PARlong,    .Dim=c(nroom*nfloor, nweeks)) 
              ,gender=structure(.Data=   genderlong, .Dim=c(nroom*nfloor, nweeks))
              ,m_age_ind=structure(.Data= m_age_indlong,.Dim=c(nroom*nfloor, nweeks))
              ,icu_ind=structure(.Data=icu_indlong, .Dim=c(nroom*nfloor, nweeks))
              ,floor=floor # nested indexing in model
              # ICAR
              ,num = adjdata$num # number of neighbors for each area
              ,adj =  adjdata$adj # adjacency matrix 
              ,sumNumNeigh = adjdata$sumNumNeigh # sum of num
)

lineinits <- function() {
  list(  beta0 = 0 # fixed
         #,beta1=0
         #,beta2=0 
         #,beta3=0 
         , tau.T = 1 
         , p=0.5 # split U:V starting value
         , U=rep(0,nroom*nfloor) # rooms and floors
         #, V=rep(0,nroom) # rooms on a floor have spatial dependence structure
         #, theta=rep(0.1, 104) # time structured part of model
         ,prec.tau2=0.001
         #,alpha0=0 # imputation
         #,tau.PAR=0.1 # imputation
  )
}

params <- c(
  "beta"
  ,"RR"
  ,"mean.U"
  ,"sd.U"
  ,"mean.V"
  ,"sd.V"
  ,"vratio"
  ,"mean.theta"
  ,"mean.delta"
  ,"sd.delta"
)

tic()
lineout <- bugs(data = mydata
                , inits = lineinits
                , parameters.to.save = params
                , model.file = mymodel
                , codaPkg = T
                , n.chains = 1
                , n.iter =  10000 # 10000
                , n.burnin =4000 # 4000
                , n.thin=1
                , debug=T
)
toc() 

#--------------------------------------------------------------------------------------------
# summarise results
# using coda
line.coda <- read.bugs(lineout)
summary(line.coda)

