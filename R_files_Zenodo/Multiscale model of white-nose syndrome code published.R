##Code for simulations and figures in "Multi-scale model fo regional population decline in little brown bats due to white-nose syndrome"
##Authors: Andrew M. Kramer, Claire S. Teitelbaum, Ashton Griffin, and John M. Drake

#Load necessary pacakges
library(adaptivetau) #For stochastic simulations via approximation of Gillespie's direct method
library(doParallel) #for parallel processing
library(lhs) #for Latin hyper cube parameter estimation
library(logspline) # For approximating population size distribution

#initial conditions and params-----------
#function for generating cave population sizes based on NY data in Turner et al 2011
#Parameters from auxiliary script "WNSBats_PowerLaw.R"
# randCavePops<-function(n,distrib="spline"){
#   if(distrib=="weibull") cave.pops<-rweibull(n, shape = 0.3155586, scale = 1508.5363225)
#   if(distrib=="lnorm") cave.pops<-rlnorm(n,meanlog = 5.596804,sdlog=3.522228)
#   if(distrib=="spline"){
#     pops<-readRDS("NY_Ml_pops.rds")
#     spline.est<-logspline(pops, lbound=1, ubound=200000)
#     cave.pops<-rlogspline(n, spline.est)
#   }
#   return(ceiling(cave.pops))
# }
#


# #Set the number of simulations
n.sims <- 100
# #Set the number of caves
nCaves <- 54 #From NY WNS cave data, O'Regan et al 2015

# #Generate cave population sizes
# set.seed(10)
# S_0 <-list()
# for (i in 1:n.sims){
#   S_0[[i]]<-randCavePops(nCaves)
# }
# saveRDS(S_0,"ListOfCavePops500.rds")
S_0 <- readRDS("ListOfCavePops500.rds")
starting_metapops<-sapply(S_0,sum)

#between cave distances
#For current case the distances are equal based on approximately uniform dispersal probability over region of interest (NY state), based on Norquay et al 2013
distances <- matrix(1, nrow = nCaves , ncol = nCaves)
diag(distances) <- 0

#population growth rates independent of WNS
##Frick *et al.* Science Table S2
geom_growth <- c(1.084,1.200,1.177
                 ,0.983,0.978,1.077
                 ,1.118,1.116,1.119,
                 1.000,0.977,1.082,
                 1.024,1.155,1.109,
                 1.131,1.049,1.008,
                 1.114,1.034,1.071,
                 1.054)

growthRange <- range(geom_growth)


#are caves infected
infected <- rep(F,nCaves)
nInfected <- 4 #Start with 4 initially infected caves as in NY dataset
infected2 <- infected
infected2[c(1:nInfected)] <- T
infected <- rbind(as.logical(infected),as.logical(infected2))

#time
winter_length <- 212 #rough estimation of average hibernating days, Meyer 2016
n_years <- 10

#dispersal/fidelity based on Norquay et al 2013, this is probability of hibernating in same cave in subsequent year
fidelity <- 0.96 #Main value is 0.96, did LowFidelity with 0.92 and HighFidelity with 1.0

#parameter ranges to test
phi <- c(12.5,200) # Meyer et al 2016
epsilon <- c(0,1) # Meyer et al 2016
psi <- c(0,0.4) # Reynolds et al. 2015
kappa <- 0.1 #This is fixed because it is not-identifiable apart from psi, and we can save lhs computation time
beta1 <- log10(c(10^(-13),10^(-11))) # Meyer et al 2016
beta2_FD <- c(0.00014,0.21) # Meyer et al 2016
beta2_DD <- c(0.00001,0.18) #beta2 for density dependent transmission
beta2 <- log10(range(c(beta2_DD , beta2_FD))) #sample across whole range, log transform
gamma1 <- c(1/100,1/50) #Lorch et al 2011
gamma2 <- c(1/80,1/10) #Lorch et al 2011
mu <- c(1/80,1/10) #1/gamma1+ 1/gamma2+ 1/mu ~ 150 : Lorch et al. 2011
alpha <- log10(c(0.9995,0.999999)) #fraction of times visit to infected cave DOES NOT result in bat infection
V <- c(1,14) #Cave visits during swarming, rough guess
#d = c(0,1) #in this version, we will not fit this with LHS, but will use 0 and 1 separately

#number of combinations to be generated in LHS, affects number of simulations and sampling of parameter ranges.
ncombs <- 1000

#use lhs to get parameter sets
# set.seed(457)
# # #---LHS routine, takes time, cannot be parallelized (1 day for ncombs=1000)
#  lhsPts = optimumLHS(ncombs,10)
#  colnames(lhsPts) = c("phi","epsilon","psi",
#                       "beta1","beta2","gamma1" , "gamma2","mu",
#                       "alpha","V")
#  saveRDS(lhsPts,"lhsPts1000_15p.rds")
lhsPts<-readRDS("lhsPts1000_10p.rds")

# #get parameter sets using the lhs
# lhsVals <- sapply(colnames(lhsPts),
#                   function(x){
#                     parmRange = get(x)
#                     (max(parmRange)-min(parmRange))*lhsPts[,x]+min(parmRange)
#                   }
# )
# lhsVals[,4] = 10^lhsVals[,4] #This back-transforms the environment to bat transmission probability to original scale
# lhsVals[,5] = 10^lhsVals[,5] #This back-transforms the bat-to-bat transmission probability to the original scale if needed
# lhsVals[,9] = 10^lhsVals[,9] #This back-transforms alpha to the original scale
# 
# #Include fixed value for kappa to satisfy mechanistic model
# lhsVals <- cbind(lhsVals[,1:3],
#                  kappa=rep(kappa,dim(lhsVals)[1]),
#                  lhsVals[,4:ncol(lhsVals)]) #add back in kappa
# 
# #now add "d" back in (because it is no longer parts of lhs):
# lhsVals <- rbind(lhsVals,
#                  lhsVals)
# 
# lhsVals <- cbind(lhsVals,
#                  d=c(rep(0,ncombs)
#                      ,rep(1,ncombs)))
# write.csv(lhsVals, 'lhsVals1000_16.csv',row.names=F)

# Read in pre-generated lhs parameters sets
lhsVals <- read.csv("lhsVals1000_16.csv") 

##FUNCTION FOR SIMULATIONS
multi_model <<- function(S_0,
                         n_years,
                         nInfected,
                         distances,
                         geom_growth,
                         winter_length,
                         allparms,
                         fidelity){
  
  #Set our cave populations
  cavePops <- matrix(S_0 , byrow=T , nrow = 1)
  nCaves <- length(S_0)
  
  #this institutes a carrying capacity
  max_pop <- sum(cavePops[1,]) 
  winterPops <- NULL
  annGrowthRates <- NULL
  
  #initialize infective bats entering caves
  swarming_infections <- matrix(c(rep(1,nInfected),
                                  rep(0,length(S_0)-nInfected)),
                                ncol=length(S_0)) 
  
  #intialize the recorder of infected caves
  infected <- matrix(c(rep(FALSE,length(S_0)),rep(TRUE,nInfected),rep(FALSE,length(S_0)-nInfected)),
                     nrow=2,
                     byrow=TRUE) 
  
  
  #Event Matrix
  nu <- matrix(c(1, -1,  0 , 0 , 0 , 0 , 0 , #what happens to E?
                 0,  0, -1 , -1, 0 , 0 , 0 , #what happens to S?
                 0,  0,  1 , 1 , -1, 0 , 0 , #what happens to I?
                 0,  0,  0 , 0 , 1 , -1, 0 , #what happens to F1?
                 0,  0,  0 , 0 , 0 , 1 , -1), #what happens to F2
               nrow = 5 , ncol = 7 , byrow = T) #change in number given each event
  
  #Rate function
  a <- function(x, parms, t){
    
    return(c(parms["phi"]*(parms["epsilon"]*x["I"]+x["F1"]+x["F2"]) + #rate of spores growing from bat shedding
               ifelse(x["E"]<10^10 , parms["psi"]*x["E"] , 0), #rate of spore growth in environment
             parms["kappa"]*x["E"], #rate of spore death in environemtn
             parms["beta1"]*x["S"]*x["E"], #rate of new infections
             ifelse((x["S"]+x["I"]+x["F1"]+x["F2"])>0,
                    parms["beta2"]*x["S"]*(parms["epsilon"]*x["I"]+x["F1"]+x["F2"])/
                      (x["S"]+x["I"]+x["F1"]+x["F2"])^parms["d"], 0), #This corrects for division by zero when all bats die during winter
             parms["gamma1"]*x["I"],
             parms["gamma2"]*x["F1"],
             parms["mu"]*x["F2"])) #propensity vector (functions of vars)
  }
  
  
  withinWinter <- NULL
  
  for(i in 2:(n_years+1)){
    
    #within cave model-----------
    gillespieResult <- NULL
    withinWinterYear <- NULL
    
    #gillespie simulation independently for each cave
    for(cave in 1:nCaves){
      
      parms = sapply(allparms , '[' , cave) #extract parameters for focal cave
      names(parms) = c("phi","epsilon","psi","kappa","beta1","beta2",
                       "gamma1","gamma2","mu","alpha","V","d")
      
      if(infected[i,cave] & !infected[i-1,cave]){ #Newly infected
        E = 0
        I = swarming_infections[i-1,cave] # row i-1 is the infected bats coming out of fall swarming prior to hibernation, 1 bat in first year
      } else if(infected[i,cave] & infected[i-1,cave]) { #Long term infection
        E = round(as.numeric(winterPops[[i-2]][cave,"E"]*exp((parms["psi"]-parms["kappa"])*(365-winter_length)))) #row i-2 in "winterpops" corresponds to row i-1 in "infected"
        I = swarming_infections[i-1,cave] # row i-1 is the infected bats coming out of fall swarming prior to hibernation
      } else if(!infected[i,cave]) {
        E = 0
        I = 0
      }
      
      if(infected[i,cave]){ #only run gillespie simulation for infected caves
        x0 <- c(E = E,
                S = cavePops[nrow(cavePops), cave][[1]] ,
                I = I , F1 = 0 , F2 = 0) #extract initial conditions for focal cave
        
        ##Approximate solution to birth-death process using 'adaptivetau' algorithm
        out<-ssa.adaptivetau(x0  ,nu , a , parms  ,tf= winter_length ,
                             deterministic = c(TRUE,TRUE,FALSE,FALSE,FALSE))
        
        gillespieResult = rbind(gillespieResult , out[nrow(out),]) #save results for cave
        
        #save a few time points from within the winter
        timePoints = seq(0  ,winter_length , 30)
        toDraw = sapply(timePoints , function(x) which.min(abs( out[,"time"] - x )))
        withinWinterYear[[cave]] = out[toDraw , ]
        withinWinterYear[[cave]] = cbind(withinWinterYear[[cave]], contactTransRate = parms["beta2"] *
                                           out[toDraw,"S"] * (parms["epsilon"]*out[toDraw,"I"] + out[toDraw,"F1"] + out[toDraw,"F2"]))
        withinWinterYear[[cave]] = cbind(withinWinterYear[[cave]], enviroTransRate = parms["beta1"] *
                                           out[toDraw,"S"] *out[toDraw,"E"])
        # }
      } else { #otherwise assume 100% winter survival
        gillespieResult = rbind(gillespieResult , cbind(NA , E = 0 , S = cavePops[nrow(cavePops),cave] ,
                                                        I = 0 , F1 = 0 , F2 = 0))
        withinWinterYear[[cave]] =  cbind(NA , E = 0 , S = cavePops[nrow(cavePops),cave] ,
                                          I = 0 , F1 = 0 , F2 = 0 , contractTransRate = NA ,
                                          enviroTransRate = NA)
      }
      
    }
    
    withinWinter[[i-1]] = withinWinterYear
    
    gillespieResult = as.data.frame(gillespieResult[,c("E","S","I","F1" , "F2")])
    winterPops[[i-1]] = gillespieResult
    gillespieResult$N = apply(gillespieResult[,c("S","I","F1" , "F2")] , 1 , sum) #get end population size for each cave
    
    
    
    #between cave model -------------
    #population growth
    N_0 = gillespieResult$N #initial number for population growth model
    if(sum(N_0) * min(geom_growth) <= max_pop){
      #if the population size is larger than initially (before WNS, assumed to be at carrying capacity):
      #then just choose a different growth rate until you get a result that's below carrying capacity
      #first, choose whether to go into while loop or just bump down to carrying capacity
      #only go into loop if there is a growth rate that will get you below carrying capacity
      growthRatesRand = sample(geom_growth) #randomize order of selecting growth rate from known range
      count = 1
      growthRate = growthRatesRand[count]
      N_disp = N_0 * growthRate
      while(sum(N_disp) > max_pop){
        growthRate = growthRatesRand[count]
        N_disp = round(N_0 * growthRate)
        count = count+1
      }
    } else {
      N_disp = (max_pop/sum(N_0))*N_0
      growthRate = sum(N_disp) / sum(N_0)
    }
    annGrowthRates = c(annGrowthRates , growthRate)
    
    #dispersal
    #bats disperse with probability "fidelity"
    sameCave = sapply(N_disp , function(n){
      disp = rbinom(n=n,size=1,prob=.96)
      return(sum(disp))
    })
    dispersers = N_disp - sameCave
    #choose destination cave randomly
    destinations = sapply(1:nCaves , function(x){
      destinations = sample(c(1:nCaves)[-x] , dispersers[x] , replace = T)
    } )
    destinations = unlist(destinations)
    destinations = sapply(c(1:nCaves),function(x) length(which(destinations==x)))
    
    new_N = sameCave + destinations
    
    # incoming infected bats at beginning of winter
    alpha <- parms["alpha"]
    V <- parms["V"]
    infectedEndWinter <- ifelse(winterPops[[i-1]][,"E"]>=1, infected[i,] , F) #These are the caves with potential for transmitting spores, currently does not depend on number of spores in cave
    nInfected <- length(which(infectedEndWinter))
    nUninfected <- length(which(!infectedEndWinter))
    infectedbats <- sapply(new_N , function(n){
      temp <- rbinom(n=n,size=1,prob= 1-(nUninfected/nCaves + (1 - nUninfected/nCaves)*alpha)^V[1])
      return(sum(temp))
    })
    swarming_infections <- rbind(swarming_infections, as.vector(infectedbats)) # Keeps track of the infected bats after each swarm period
    infected <- rbind(infected, infectedEndWinter | (infectedbats > 0)) # Adds new record of infected with holdovers and newly infected
    S_0 <- new_N
    cavePops <- rbind(cavePops , as.numeric(S_0))
  }
  
  
  return(list(cavePops = cavePops , winterPops = winterPops ,
              annGrowthRates = annGrowthRates , infected = infected[2:nrow(infected),],
              withinWinter = withinWinter, swarming_infections=swarming_infections))
}



##################################Base simulations############################
lappend <<- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

#These are run in parallel
#Output is written to file and then aggregated/distilled afterwards
n.cores <- 30
cl <- makeCluster(n.cores) #Sets the number of cores
registerDoParallel(cl) #Registers the cores with the server.
options<-list(preschedule=FALSE)
start.time <- Sys.time()
lhs.arc <- lhsVals
set.seed(160, kind = "L'Ecuyer-CMRG")
output <- foreach (z = 1:length(lhsVals[,1]), .options.snow=options) %dopar% {
  
  #library(GillespieSSA)
  library(adaptivetau)
  allparms = lapply(colnames(lhsVals) , function(x){
    rep(lhsVals[z,x],nCaves)
  })
  
  wns.sims <- list()
  #n.sims = 100
  
  #results_base = list()
  for(i in 1:n.sims){
    wns.sims[[i]] = multi_model(S_0[[i]] , n_years , nInfected , distances , geom_growth , winter_length , allparms , fidelity)
    #try(saveRDS(list(lhsval=z,results=wns.sims[[i]]),paste("WNS_sims",z,i,".rds",sep="_")))
  }
  
  output <- list(lhsval = z, results = wns.sims)
  saveRDS(output,paste("outputProgress16Base/wns_output_nc30_lhs1000_100sims",z,".rds",sep=""))
}

stopCluster(cl)


#High Fidelity simulations
#dispersal/fidelity
fidelity <- 1 #Main value is 0.96, did LowFidelity with 0.92 and HighFidelity with 1.0

n.cores <- 30
cl <- makeCluster(n.cores) #Sets the number of cores
registerDoParallel(cl) #Registers the cores with the server.
options<-list(preschedule=FALSE)
start.time <- Sys.time()
lhs.arc <- lhsVals
set.seed(160, kind = "L'Ecuyer-CMRG")
output <- foreach (z = 1:length(lhsVals[,1]), .options.snow=options) %dopar% {
  
  #library(GillespieSSA)
  library(adaptivetau)
  allparms = lapply(colnames(lhsVals) , function(x){
    rep(lhsVals[z,x],nCaves)
  })

  wns.sims <- list()
  #n.sims = 100
  
  #results_base = list()
  for(i in 1:n.sims){
    wns.sims[[i]] = multi_model(S_0[[i]] , n_years , nInfected , distances , geom_growth , winter_length , allparms , fidelity)
    #try(saveRDS(list(lhsval=z,results=wns.sims[[i]]),paste("WNS_sims",z,i,".rds",sep="_")))
  }
  
  output <- list(lhsval = z, results = wns.sims)
  saveRDS(output,paste("outputProgress16HighFidelity/wns_output_nc30_lhs1000_100sims",z,".rds",sep=""))
}

stopCluster(cl)


##Low Fidelity
fidelity <- 0.92 #Main value is 0.96, did LowFidelity with 0.92 and HighFidelity with 1.0

n.cores <- 30
cl <- makeCluster(n.cores) #Sets the number of cores
registerDoParallel(cl) #Registers the cores with the server.
options<-list(preschedule=FALSE)
start.time <- Sys.time()
lhs.arc <- lhsVals
set.seed(160, kind = "L'Ecuyer-CMRG")
output <- foreach (z = 1:length(lhsVals[,1]), .options.snow=options) %dopar% {
  
  #library(GillespieSSA)
  library(adaptivetau)
  allparms = lapply(colnames(lhsVals) , function(x){
    rep(lhsVals[z,x],nCaves)
  })
 
  wns.sims <- list()
  #n.sims = 100
  
  #results_base = list()
  for(i in 1:n.sims){
    wns.sims[[i]] = multi_model(S_0[[i]] , n_years , nInfected , distances , geom_growth , winter_length , allparms , fidelity)
    #try(saveRDS(list(lhsval=z,results=wns.sims[[i]]),paste("WNS_sims",z,i,".rds",sep="_")))
  }
  
  output <- list(lhsval = z, results = wns.sims)
  saveRDS(output,paste("outputProgress16LowFidelity/wns_output_nc30_lhs1000_100sims",z,".rds",sep=""))
}

stopCluster(cl)



########################Deterministic Version#############
##These are run with within-cave dynamics determinstic
library(doParallel) #for parallel processing
library(lhs) #for Latin hyper cube parameter estimation
library(deSolve) #differential equation solver

#initial conditions and params-----------
#initial conditions
S_0 = readRDS("ListOfCavePops500.rds")
S_0_all = S_0
nCaves = length(S_0[[1]])
#between cave distances
# A--10--B---15---C
distances = matrix(1, nrow = nCaves , ncol = nCaves)
diag(distances) = 0
#population growth rates independent of WN
geom_growth = c(1.084,1.200,1.177,0.983,0.978,1.077,1.118,1.116,1.119,1.000,0.977,1.082,1.024,1.155,1.109,1.131,1.049,1.008,1.114,1.034,1.071,1.054)
growthRange = range(geom_growth)
winter_length <<- 212 #rough estimation of average hibernating days, Meyer 2016
n_years <<- 10
nInfectedInit=4


##Function
determ_model <<- function(S_0 , n_years , nInfectedInit , distances , geom_growth , winter_length,
                          allparms , fidelity=0.96){
  
  cavePops = matrix(S_0 , byrow=T , nrow = 1)
  nCaves = length(S_0)
  max_pop = sum(cavePops[1,]) #this institutes a carrying capacity
  swarming_infections = matrix(c(rep(1,nInfectedInit),rep(0,length(S_0)-nInfectedInit)),
                               ncol=length(S_0)) #initialize infective bats entering caves
  infected = matrix(c(rep(FALSE,length(S_0)),rep(TRUE,nInfectedInit),rep(FALSE,length(S_0)-nInfectedInit))
                    ,nrow=2,byrow=TRUE) #intialize the recorder of infected caves
  
  
  withinCaveModel = function(t,x,parameters){
    E = ifelse(x[1]>0,x[1],0)
    S = ifelse(x[2]>0,x[2],0)
    I = ifelse(x[3]>0,x[3],0)
    F1 = ifelse(x[4]>0,x[4],0)
    F2 = ifelse(x[5]>0,x[5],0)
    
    with(as.list(parameters),{
      dE = phi*(epsilon*I+F1+F2) - kappa*E + psi*E
      dS = (-1)*beta1*S*E - beta2*S*(epsilon*I+F1+F2)/(S+I+F1+F2)^d
      dI = ifelse((S+I+F1+F2)>0,
                  beta1*S*E + (beta2*S*(epsilon*I+F1+F2))/(S+I+F1+F2)^d - gamma1*I,
                  0)
      dF1 = gamma1*I - gamma2*F1
      dF2 = gamma2*F1 - mu*F2
      dx = c(dE,dS,dI,dF1,dF2)
      list(dx)
    })
    
  }
  
  withinWinter = NULL
  winterPops = NULL
  annGrowthRates = NULL
  for(i in 2:(n_years+1)){
    
    #within cave model-----------
    determResult = NULL
    withinWinterYear = vector(mode="list",length(S_0))
    
    # simulation independently for each cave
    for(cave in 1:nCaves){
      
      if(infected[i,cave] & !infected[i-1,cave]){ #Newly infected
        E = 0 
        I = swarming_infections[i-1,cave] # row i-1 is the infected bats coming out of fall swarming prior to hibernation, 1 bat in first year
      } else if(infected[i,cave] & infected[i-1,cave]) { #Long term infection
        E = winterPops[[i-2]][cave,"E"]*exp((parms["psi"]-parms["kappa"])*(365-winter_length)) #row i-2 in "winterpops" corresponds to row i-1 in "infected"
        E = round(as.numeric(E))
        I = swarming_infections[i-1,cave] # row i-1 is the infected bats coming out of fall swarming prior to hibernation
      } else if(!infected[i,cave]) {
        E = 0
        I = 0
      } 
      
      x0 <- round(c(E = E,
                    S = cavePops[nrow(cavePops), cave][[1]] , 
                    I = I , F1 = 0 , F2 = 0))
      
      parms = sapply(allparms , '[' , cave) #extract parameters for focal cave
      names(parms) = c("phi","epsilon","psi","kappa",
                       "beta1","beta2","gamma1","gamma2","mu",
                       "alpha","V","d")
      
      if(infected[i,cave] & (x0["S"] > 0 | x0["I"] > 0  | x0["F1"] > 0 | x0["F2"] > 0)){
        out <- tryCatch(
          as.data.frame(lsoda(x0,seq(1,winter_length,1),withinCaveModel,parms,
                              rtol=10^(-4) , atol=10^(-4))),
          error = function(x){
            out <- as.data.frame(matrix(rep(0,12),nrow=2))
            names(out) = c("time","E","S","I","F1","F2")
            return(out)
            # withinWinter[[i-1]] = NA
          } )
        if(any(na.omit(unlist(out))<0)){
          out <- apply(out,2,function(x){
            orig = x
            new=x
            new[orig<0]=0
            return(new)
          }) #change any negative vals to zero
        }
        
        timePoints = seq(0  ,winter_length , 30)
        toDraw = sapply(timePoints , function(x) which.min(abs( out[,"time"] - x )))
        withinWinterYear[[cave]] = out[toDraw , ]
        withinWinterYear[[cave]] = cbind(withinWinterYear[[cave]], contactTransRate =
                                           parms[["beta2"]] *
                                           out[toDraw,"S"] * (parms["epsilon"]*out[toDraw,"I"] + out[toDraw,"F1"]+ out[toDraw,"F2"]))
        withinWinterYear[[cave]] = cbind(withinWinterYear[[cave]], enviroTransRate = parms["beta1"] *
                                           out[toDraw,"S"] *out[toDraw,"E"])
        
      } else {
        out <- as.data.frame(matrix(c(1,1,sapply(x0,rep,2)),nrow=2))
        names(out) = c("time","E","S","I","F1","F2")
        withinWinterYear[[cave]] = NA
      }
      
      determResult = rbind(determResult , out[max(which(!is.na(out[,"E"]))),]) #save results for cave
      withinWinter[[i-1]] = withinWinterYear #save results for cave
      
    }
    
    determResult = as.data.frame(determResult[,c("E","S","I","F1","F2")])
    winterPops[[i-1]] = determResult
    determResult$N = apply(determResult[,c("S","I","F1","F2")] , 1 , sum) #get end population size for each cave
    
    
    
    #between cave model -------------
    #population growth
    N_0 = determResult$N #initial number for population growth model
    if(sum(N_0) * min(geom_growth) <= max_pop){
      #if the population size is larger than initially (before WNS, assumed to be at carrying capacity):
      #then just choose a different growth rate until you get a result that's below carrying capacity
      #first, choose whether to go into while loop or just bump down to carrying capacity
      #only go into loop if there is a growth rate that will get you below carrying capacity
      growthRatesRand = sample(geom_growth) #randomize order of selecting growth rate from known range
      count = 1
      growthRate = growthRatesRand[count]
      N_disp = N_0 * growthRate
      while(sum(N_disp) > max_pop){
        growthRate = growthRatesRand[count]
        N_disp = round(N_0 * growthRate)
        count = count+1
      }
    } else {
      N_disp = (max_pop/sum(N_0))*N_0
      growthRate = sum(N_disp) / sum(N_0)
    }
    annGrowthRates = c(annGrowthRates , growthRate)
    
    #dispersal
    #bats disperse with probability "fidelity"
    sameCave = sapply(N_disp , function(n){
      disp = rbinom(n=n,size=1,prob=.96)
      return(sum(disp))
    })
    dispersers = N_disp - sameCave
    #choose destination cave randomly
    destinations = sapply(1:nCaves , function(x){
      destinations = sample(c(1:nCaves)[-x] , dispersers[x] , replace = T)
    } )
    destinations = unlist(destinations)
    destinations = sapply(c(1:nCaves),function(x) length(which(destinations==x)))
    
    new_N = sameCave + destinations
    
    # incoming infected bats at beginning of winter---- NEW VERSION
    alpha = parms["alpha"]
    V = parms["V"]
    infectedEndWinter <- ifelse(winterPops[[i-1]][,"E"]>=1, infected[i,] , F) #These are the caves with potential for transmitting spores, currently does not depend on number of spores in cave
    nInfected <- length(which(infectedEndWinter))
    nUninfected <- length(which(!infectedEndWinter))
    
    infectedbats <- sapply(new_N , function(n){
      temp = rbinom(n=n,size=1,prob= 1-(nUninfected/nCaves + (1 - nUninfected/nCaves)*alpha)^V[1])
      return(sum(temp))
    })
    
    #attempt at making this deterministic. not working right now because it produces non-integers
    # infectedbats <- sapply(new_N , function(n){
    #   temp = n*(1-(nUninfected/nCaves + (1 - nUninfected/nCaves)*alpha)^V[1])
    #   return(sum(temp))
    # })
    
    swarming_infections <- rbind(swarming_infections, as.vector(infectedbats)) # Keeps track of the infected bats after each swarm period
    infected <- rbind(infected, infectedEndWinter | (infectedbats > 0)) # Adds new record of infected with holdovers and newly infected
    S_0 = new_N
    cavePops = rbind(cavePops , as.numeric(S_0))
  }
  
  
  return(list(cavePops = cavePops , winterPops = winterPops , 
              annGrowthRates = annGrowthRates , infected = infected[2:nrow(infected),],
              withinWinter = withinWinter, swarming_infections=swarming_infections))
}

##########################Determinstic simulations
cl = makeCluster(30)
registerDoParallel(cl) #Registers the cores with the server.


output = vector(mode="list",length=nrow(lhsVals))

output <- foreach (z = 1:nrow(lhsVals)) %dopar% {
  
  library(deSolve)
  
  out1 = vector(mode="list",length=length(S_0_all))
  
  for(sim in 1:length(S_0_all)){
    
    allparms = lapply(lhsVals[z,],function(x) rep(x,nCaves))
    
    modelOut = determ_model(S_0_all[[sim]] , n_years , nInfectedInit , distances , geom_growth ,
                            winter_length,allparms , fidelity=0.96 )
    
    out1[[sim]] = modelOut
    
  }
  
  fileName = paste0("determ_sims16/output_deterministic_1000LHS_100sim_", z,".rds")
  
  saveRDS(out1,fileName)
  out1
}

saveRDS(output,"output_deterministic_1000LHS_100sim_16.rds")


stopCluster(cl)

output<-list()
for (i in 1:2000){
  fileName = paste0("determ_sims16/output_deterministic_1000LHS_100sim_", i,".rds")
  output[[i]]<-readRDS(fileName)
}

#Evaluation of deterministic output

prevalenceY1=NULL
prevalenceY2=NULL
propSurviveY1=NULL
propSurviveY2=NULL
yearsToFullInfection=NULL
nInfectedY2 = NULL
nInfectedY3 = NULL
nInfectedY4 = NULL
nInfectedY5 = NULL
sporeCount = NULL
batPopY5 = NULL

#this is what we're getting as our output for each metric
statToCompute = function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}

##Read through all the simulation output and extract summary stats
for(i in 1:2000){
  #dats =  readRDS(paste0("determ_sims/output_deterministic_1000LHS_100sim_", i,".rds"))
  dats= output[[i]]
  
  #proportion surviving at cave scale in year 1
  propSurviveY1_temp = sapply(dats , function(y){
    startPops = y[["cavePops"]]
    winterPops = y[["winterPops"]]
    endPops = t(sapply(winterPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum)))
    infected = y[["infected"]] ; infected = infected[2:nrow(infected),]
    firstInfection = apply(infected,2,function(cave) which(cave)[1])
    propSurvive = sapply(1:length(firstInfection) , function(x){
      endPops[firstInfection[x] , x]/
        startPops[firstInfection[x] , x]}) 
    mean(na.omit(propSurvive))
  })
  propSurviveY1 = cbind(propSurviveY1 , statToCompute(propSurviveY1_temp))
  
  
  #proportion surviving at cave scale in year 2
  propSurviveY2_temp = sapply(dats , function(y){
    startPops = y[["cavePops"]]
    winterPops = y[["winterPops"]]
    endPops = t(sapply(winterPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum)))
    infected = y[["infected"]] ; infected = infected[2:nrow(infected),]
    firstInfection = apply(infected,2,function(cave) which(cave)[1])
    firstInfection = ifelse(firstInfection >= nrow(endPops) , NA , firstInfection)
    
    propSurvive = sapply(1:length(firstInfection) , function(x){
      endPops[firstInfection[x]+1 , x]/
        startPops[firstInfection[x]+1 , x]}) 
    mean(na.omit(propSurvive))
  })
  propSurviveY2 = cbind(propSurviveY2 , statToCompute(propSurviveY2_temp))
  
  
  #prevalence at end of winter in year 1
  prevalenceY1_temp = sapply(dats , function(y){
    winterPops = y[["winterPops"]]
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    firstInfection = apply(infected,2,function(cave) which(cave)[1])
    
    winterPopYearInfected = sapply(1:length(firstInfection) ,
                                   function(x) winterPops[[firstInfection[x]]][x,])
    prev = sapply(winterPopYearInfected , function(cave){
      sum(cave[c("I","F1","F2")])/sum(cave[c("S","I","F1","F2")])
    })
    mean(na.omit(prev))
  })
  prevalenceY1 = cbind(prevalenceY1 , statToCompute(prevalenceY1_temp))
  
  
  #prevalence at end of dats in and after year 2
  prevalenceY2_temp = sapply(dats , function(y){
    winterPops = y[["winterPops"]]
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    firstInfection = apply(infected,2,function(cave) which(cave)[1])
    firstInfection = ifelse(firstInfection >= length(winterPops) , NA , firstInfection)
    
    winterPopYearInfected = sapply(1:length(firstInfection) ,
                                   function(x) winterPops[[firstInfection[x]+1]][x,])
    prev = sapply(winterPopYearInfected , function(cave){
      sum(cave[c("I","F1","F2")])/sum(cave[c("S","I","F1","F2")])
    })
    mean(na.omit(prev))
  })
  prevalenceY2 = cbind(prevalenceY2 , statToCompute(prevalenceY2_temp))
  
  #Time until all caves infected
  yearsToFullInfection_temp = sapply(dats , function(y){
    infected = y[["infected"]]
    nInfected = rowSums(infected)
    yearsToInfected = which(nInfected==ncol(infected))[1]
    if (is.na(yearsToInfected)) yearsToInfected <- 12
    yearsToInfected
  })
  yearsToFullInfection = cbind(yearsToFullInfection , statToCompute(yearsToFullInfection_temp))
  
  
  #spore count
  sporeCount_temp = sapply(dats , function(y){
    endPops = y[["winterPops"]]
    endE = t(sapply(endPops,function(X) X$E)) #number spores at end of winter in all caves
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    infectedE = endE[which(infected)] #number of spores at end of winter in infected caves
    mean(infectedE)
  })
  sporeCount = cbind(sporeCount , statToCompute(sporeCount_temp))
  
  
  #new infections in each year
  nInfectedY2_temp = sapply(dats , function(y){
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    length(which(!infected[1,] & infected[2,]))
  })
  nInfectedY2 = cbind(nInfectedY2 , statToCompute(nInfectedY2_temp))
  
  nInfectedY3_temp = sapply(dats , function(y){
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    length(which(!infected[2,] & infected[3,]))
  })
  nInfectedY3 = cbind(nInfectedY3 , statToCompute(nInfectedY3_temp))
  
  nInfectedY4_temp = sapply(dats , function(y){
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    length(which(!infected[3,] & infected[4,]))
  })
  nInfectedY4 = cbind(nInfectedY4 , statToCompute(nInfectedY4_temp))
  
  nInfectedY5_temp = sapply(dats , function(y){
    infected = y[["infected"]]  ; infected = infected[2:nrow(infected),]
    length(which(!infected[4,] & infected[5,]))
  })
  nInfectedY5 = cbind(nInfectedY5 , statToCompute(nInfectedY5_temp))
  
  batPopY5_temp<-sapply(dats, function(y){
    cavePops=sum(y[["cavePops"]][5,])
  })
  batPopY5=cbind(batPopY5,statToCompute(batPopY5_temp))
  
  rm(dats)
  
  if(i %in% seq(1,2000,20)){
    print(paste((i/2000)*100 , "% complete at" , Sys.time()))
  }
}

evaluationMetricsDeterm = list(
  prevalenceY1=prevalenceY1,
  prevalenceY2=prevalenceY2,
  propSurviveY1=propSurviveY1,
  propSurviveY2=propSurviveY2,
  yearsToFullInfection=yearsToFullInfection,
  nInfectedY2 = nInfectedY2,
  nInfectedY3 = nInfectedY3,
  nInfectedY4 = nInfectedY4,
  nInfectedY5 = nInfectedY5,
  sporeCount = sporeCount,
  batPopY5 = batPopY5)


saveRDS(evaluationMetricsDeterm,
        "evaluationMetricsDeterm_16.rds")
rm(output)
gc()


#####################Evaluation for base simulations
##Read through all the simulation output and extract summary stats
##This is rerun later in the file for Low and High fidelity)
output <- NULL
for(i in 1:2000){
  output[[i]] = readRDS(paste("outputProgress16Base/wns_output_nc30_lhs1000_100sims",i,".rds",sep=""))
}
saveRDS(output,file='wns_output_nc30_lhs1000_100sims_16.rds')
#output<-readRDS('wns_output_nc30_lhs1000_100sims_16.rds')

length(output) == nrow(lhsVals) #Consistency check

allResults = lapply(1:length(output) , function(x){
  results = output[[x]][['results']]
  a = lapply(names(results[[1]]),function(y){
    lapply(1:length(results),function(z){
      new = results[[z]][[y]]
    })
  })
  names(a) = names(results[[1]])
  return(a)
})

#this is what we're getting as our output for each metric
statToCompute <- function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}


#proportion surviving at cave scale in year 1
propSurviveY1 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year1Survival = sapply(which(!is.na(yearInfected) & yearInfected <= length(endPops)) ,
                           function(cave){
                             endPops[[yearInfected[cave]]][cave]/
                               startPops[yearInfected[cave],cave]
                           })
    mean(year1Survival) #This is the average over all the infected caves in the simulation
  })
  statToCompute(propSurvive) #The appropriate quantiles of the mean survival from each set of simulations
})


#proportion surviving at cave scale in year 2
#this is beginning to end of winter, not through entire year
propSurviveY2 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year2Survival = sapply(which(!is.na(yearInfected) & yearInfected <= (length(endPops)-1)) , function(cave){
      endPops[[yearInfected[cave]+1]][cave]/
        startPops[yearInfected[cave]+1,cave]
    })
    mean(year2Survival)
  })
  statToCompute(propSurvive)
})


#proportion surviving at cave scale after year 2
propSurviveY3 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    allSurvival = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-2)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],"[",cave)
      winterStart = cavePops[[sim]][allOtherYears,cave]
      winterEnd/winterStart
    })
    mean(unlist(allSurvival))
  })
  statToCompute(propSurvive)
})


#prevalence at end of winter in year 1
prevalenceY1 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    year1Prev = sapply(which(!is.na(yearInfected)& yearInfected <= length(endPops)) , function(cave){
      sum(endPops[[yearInfected[cave]]][cave,c("I","F1","F2")])/
        sum(endPops[[yearInfected[cave]]][cave,c("S","I","F1","F2")])
    })
    mean(year1Prev)
  })
  statToCompute(prev)
})


#prevalence at end of winter in and after year 2
prevalenceY2 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    allPrev = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-1)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],function(x){
        sum(x[cave,c("S","I","F1","F2")])/sum(x[cave,c("S","I","F1","F2")])
      })
    })
    mean(unlist(allPrev))
  })
  statToCompute(prev)
})


#survival in year 5
survivalY5 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  survive = sapply(1:length(cavePops),function(sim){
    popSize = cavePops[[sim]]
    sum(popSize[5,]) / sum(popSize[1,])
  })
  statToCompute(survive)
})

yearsToFullInfection = sapply(allResults , function(y){
  infected = y[["infected"]]
  survive = sapply(1:length(infected),function(sim){
    infectedSim = infected[[sim]]
    nInfected = rowSums(infectedSim)
    yearsToInfected = which(nInfected==ncol(infectedSim))[1]
    if (is.na(yearsToInfected)) yearsToInfected <- 12
    return (yearsToInfected)
  })
  statToCompute(survive)
})


#spore count
sporeCount = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  spores = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    endE = t(sapply(endPops,function(X) X$E)) #number spores at end of winter in all caves
    infectedSim = infected[[sim]][1:10,] #which caves are infected? (all years)
    infectedE = endE[which(infectedSim)] #number of spores at end of winter in infected caves
    spores = median(infectedE) #mean spore count in infected caves for this simulation
  })
  statToCompute(spores)
})


#new infections in each year
nInfectedY2 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][1,] &
                                 infected[[sim]][2,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY3 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][2,] &
                                 infected[[sim]][3,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY4 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][3,] &
                                 infected[[sim]][4,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY5 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][4,] &
                                 infected[[sim]][5,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#Total bat population after 5 years
batPopY5<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY5=sapply(1:length(cavePops),function(sim){
    popSumY5<-sum(cavePops[[sim]][5,])
  })
  statToCompute(popY5)
})

#Total bat populations after 10 years
batPopY10<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY10=sapply(1:length(cavePops),function(sim){
    popSumY10<-sum(cavePops[[sim]][10,])
  })
  statToCompute(popY10)
})

#Proportion of caves without bats
propExtinct<-sapply(allResults, function(y){
  winterPops = y[["winterPops"]]
  popSizes = lapply(winterPops,function(x) rowSums(x[[10]][,-1]))
  extinctionProp = sapply(popSizes, function(x){
    length(which(x==0))/length(x)
  })
  statToCompute(extinctionProp)
})

saveRDS(list(prevalenceY1=prevalenceY1,
             prevalenceY2=prevalenceY2,
             propSurviveY1=propSurviveY1,
             propSurviveY2=propSurviveY2,
             survivalY5=survivalY5,
             yearsToFullInfection=yearsToFullInfection,
             nInfectedY2 = nInfectedY2,
             nInfectedY3 = nInfectedY3,
             nInfectedY4 = nInfectedY4,
             nInfectedY5 = nInfectedY5,
             sporeCount = sporeCount,
             batPopY5 = batPopY5,
             batPopY10 = batPopY10,
             propExtinct = propExtinct),
        "evaluationMetrics1000_16.rds")

allResults = lapply(1:length(output) , function(x){
  results = output[[x]][['results']]
  a = lapply(names(results[[1]]),function(y){
    lapply(1:length(results),function(z){
      new = results[[z]][[y]]
    })
  })
  names(a) = names(results[[1]])
  return(a)
})
statToCompute = function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}

#Choose for appropriate experiment (baseline, gammas, HighFidelity, LowFidelity)
evaluationMetrics = readRDS("evaluationMetrics1000_16.rds")
for(i in names(evaluationMetrics)){
  assign(i,evaluationMetrics[[i]])
}

###########Comparison to goodness of fit measures##########
##Draw various plots

#Read in the lhsVals
lhsVals <- read.csv('lhsVals1000_16.csv')

# Read in data on time series of infection in NY
observed <- read.csv("NYPA_caves_infections_DK.csv") #This is data from O'Regan et al 2015 with state labels added
observed.NY <- observed[observed$State=="NY",]
time.series.NY <- colSums(observed.NY[,2:7])
time.series.NY.c = cumsum(time.series.NY)

#what are the expected values of each of these metrics?
expectedvals =list(propSurviveY2=0.15, # Langwig et al. 2015 EID Pg. 1024 & Fig. 2
                   prevalenceY1=c(0.4,0.88), #0.88: Langwig et al. 2015 EID Fig 1A, Verant et al 2018
                   #prevalenceY1 = 0.5: Frick et al 2017 Ecology Fig 1A
                   #prevalenceY2=1, #Langwig et al 2015 EID Fig 1B
                   #prevalenceY2=0.9, #Frick et al 2017 Ecology Fig 1A
                   #survivalY5=0.1,
                   yearsToFullInfection=c(3,6), #O'Regan et al 2014
                   sporeCount = c(10^6,10^12), #Bracketing Reynolds et al 2015
                   nInfectedY2 = time.series.NY[3], #O'Regan et al 2014
                   nInfectedY3 = time.series.NY[4], #O'Regan et al 2014
                   nInfectedY4 = time.series.NY[5], #O'Regan et al 2014
                   nInfectedY5 = time.series.NY[6]) #O'Regan et al 2014


#plot each parameter and parameter set for combined frequency and density dependence
png("parameterMatching1000_revised.png",7,10,units="in",res=500)
#Ignore kappa for the figure
fig.params<-c(1:3,5:11)
#colnames(lhsVals) = c("setNum",colnames(lhsVals)[-1])
x.labels<-c(expression(phi),expression(epsilon),expression(psi),expression(kappa),expression(beta[1]),
            expression(beta[2]),expression(gamma[1]),expression(gamma[2]),expression(mu),
            expression(alpha),"V")
par(mfrow=c(length(expectedvals),length(fig.params)) , mar = c(0,0,0,0) ,
    oma = c(5,5,0.2,0.2))
for(metric in names(expectedvals)){
  for(param in fig.params){
    ymax<-max(na.omit(get(metric)["97.5%",]))+max(na.omit(get(metric)["97.5%",]))/16
    ymin<-0
    if(metric=="yearsToFullInfection"){
      ylab = "All infected (yrs)"
      ymin<-2
    }
    if(metric=="prevalenceY1"){
      ylab = "Year 1 prevalence"
    }
    if(metric=="propSurviveY2"){
      ylab = "Year 2 survival"
    } 
    if(metric=="sporeCount"){ 
      ymax<-1e15+1e15/15
      ylab = bquote(paste("Environ. ",italic("Pd")))
    }
    if(metric=="nInfectedY2"){
      ylab="Year 2 infections"
    }
    if(metric=="nInfectedY3"){
      ylab="Year 3 infections"
    }
    if(metric=="nInfectedY4"){
      ylab="Year 4 infections"
    }
    if(metric=="nInfectedY5"){
      ylab="Year 5 infections"
    }
    xmin = min(lhsVals[,param])-(max(lhsVals[,param])-min(lhsVals[,param]))/10
    xmax = max(lhsVals[,param])+(max(lhsVals[,param])-min(lhsVals[,param]))/10
    plot(x=lhsVals[,param] , y=get(metric)["50%",] ,
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),
         pch=20,axes=FALSE,xpd=TRUE)
    box()
    segments(x0=lhsVals[,param] , y0 = get(metric)["2.5%",] , y1 = get(metric)["97.5%",])
    expectedValueMetric = expectedvals[[metric]]
    if(length(expectedValueMetric)==1){
      abline(h = expectedValueMetric , col = "red")
    } else {
      polygon(c(par()$usr[1],par()$usr[2],par()$usr[2],par()$usr[1]) ,
              c(min(expectedValueMetric),min(expectedValueMetric),
                max(expectedValueMetric),max(expectedValueMetric)),
              border="red",col=rgb(1,0,0,.2))
    }
    
    if(param==1){
      mtext(ylab,side=2,line=4,cex=.8)
      axis(2,las=1)
    }
    if(metric == names(expectedvals)[length(expectedvals)]){
      axis(1,labels=FALSE)
      text(x=(axTicks(side=1)+((max(axTicks(side=1)-min(axTicks(side=1)))*0.1)))[c(1,3,length(axTicks(side=1)))],
           y=par("usr")[3]-2,
           labels=axTicks(side=1)[c(1,3,length(axTicks(side=1)))],srt=45,pos=2,xpd=NA)
      #mtext(colnames(lhsVals)[param],side=1,line=3,cex=.8)
      mtext(text=x.labels[param],side=1,line=3,cex=.8)
    }
  }
}
dev.off()


##Split plot to show the frequency and density dependent independently
##First density dependent
png("parameterMatching1000_16_DDrevised.png",7,10,units="in",res=500)
#Ignore kappa for the figure
fig.params<-c(1:3,5:11)
#colnames(lhsVals) = c("setNum",colnames(lhsVals)[-1])
par(mfrow=c(length(expectedvals),length(fig.params)) , mar = c(0,0,0,0) ,
    oma = c(5,5,0.2,0.2))
for(metric in names(expectedvals)){
  for(param in fig.params){
    ymax<-max(na.omit(get(metric)["97.5%",]))+max(na.omit(get(metric)["97.5%",]))/16
    ymin<-0
    if(metric=="yearsToFullInfection"){
      ylab = "All infected (yrs)"
      ymin = 2
    }
    if(metric=="prevalenceY1"){
      ylab = "Year 1 prevalence"
    }
    if(metric=="propSurviveY2"){
      ylab = "Year 2 survival"
    } 
    if(metric=="sporeCount"){ 
      ymax<-1e15+1e15/15
      ylab = bquote(paste("Environ. ",italic("Pd")))
    }
    if(metric=="nInfectedY2"){
      ylab="Year 2 infections"
    }
    if(metric=="nInfectedY3"){
      ylab="Year 3 infections"
    }
    if(metric=="nInfectedY4"){
      ylab="Year 4 infections"
    }
    if(metric=="nInfectedY5"){
      ylab="Year 5 infections"
    }
    xmin = min(lhsVals[,param])-(max(lhsVals[,param])-min(lhsVals[,param]))/10
    xmax = max(lhsVals[,param])+(max(lhsVals[,param])-min(lhsVals[,param]))/10
    plot(x=lhsVals[,param][1:1000] , y=get(metric)["50%",][1:1000] , #The first 1000 are d=0, i.e. density-dependent
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),
         pch=20,axes=FALSE,xpd=TRUE)
    box()
    segments(x0=lhsVals[,param][1:1000] , y0 = get(metric)["2.5%",][1:1000] , y1 = get(metric)["97.5%",][1:1000])
    expectedValueMetric = expectedvals[[metric]]
    if(length(expectedValueMetric)==1){
      abline(h = expectedValueMetric , col = "red")
    } else {
      polygon(c(par()$usr[1],par()$usr[2],par()$usr[2],par()$usr[1]) ,
              c(min(expectedValueMetric),min(expectedValueMetric),
                max(expectedValueMetric),max(expectedValueMetric)),
              border="red",col=rgb(1,0,0,.2))
    }
    
    if(param==1){
      mtext(ylab,side=2,line=4,cex=.8)
      axis(2,las=1)
    }
    if(metric == names(expectedvals)[length(expectedvals)]){
      axis(1,labels=FALSE)
      text(x=(axTicks(side=1)+((max(axTicks(side=1)-min(axTicks(side=1)))*0.1)))[c(1,3,length(axTicks(side=1)))],
           y=par("usr")[3]-2,
           labels=axTicks(side=1)[c(1,3,length(axTicks(side=1)))],srt=45,pos=2,xpd=NA)
      mtext(text=x.labels[param],side=1,line=3,cex=.8)
    }
  }
}
dev.off()


#Plot for frequency-dependent transmission
png("parameterMatching1000_16_FDrevised.png",7,10,units="in",res=500)
#Ignore kappa for the figure
fig.params<-c(1:3,5:11)
#colnames(lhsVals) = c("setNum",colnames(lhsVals)[-1])
par(mfrow=c(length(expectedvals),length(fig.params)) , mar = c(0,0,0,0) ,
    oma = c(5,5,0.2,0.2))
for(metric in names(expectedvals)){
  for(param in fig.params){
    ymax<-max(na.omit(get(metric)["97.5%",]))+max(na.omit(get(metric)["97.5%",]))/16
    ymin<-0
    if(metric=="yearsToFullInfection"){
      ylab = "All infected (yrs)"
      ymin = 2
    }
    if(metric=="prevalenceY1"){
      ylab = "Year 1 prevalence"
    }
    if(metric=="propSurviveY2"){
      ylab = "Year 2 survival"
    } 
    if(metric=="sporeCount"){ 
      ymax<-1e15+1e15/15
      ylab = bquote(paste("Environ. ",italic("Pd")))
    }
    if(metric=="nInfectedY2"){
      ylab="Year 2 infections"
    }
    if(metric=="nInfectedY3"){
      ylab="Year 3 infections"
    }
    if(metric=="nInfectedY4"){
      ylab="Year 4 infections"
    }
    if(metric=="nInfectedY5"){
      ylab="Year 5 infections"
    }
    xmin = min(lhsVals[,param])-(max(lhsVals[,param])-min(lhsVals[,param]))/10
    xmax = max(lhsVals[,param])+(max(lhsVals[,param])-min(lhsVals[,param]))/10
    plot(x=lhsVals[,param][1001:2000] , y=get(metric)["50%",][1001:2000] , #The second 1000 are d=1, i.e. frequency-dependent
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),
         pch=20,axes=FALSE,xpd=TRUE)
    box()
    segments(x0=lhsVals[,param][1001:2000] , y0 = get(metric)["2.5%",][1001:2000] , y1 = get(metric)["97.5%",][1001:2000])
    expectedValueMetric = expectedvals[[metric]]
    if(length(expectedValueMetric)==1){
      abline(h = expectedValueMetric , col = "red")
    } else {
      polygon(c(par()$usr[1],par()$usr[2],par()$usr[2],par()$usr[1]) ,
              c(min(expectedValueMetric),min(expectedValueMetric),
                max(expectedValueMetric),max(expectedValueMetric)),
              border="red",col=rgb(1,0,0,.2))
    }
    
    if(param==1){
      mtext(ylab,side=2,line=4,cex=.8)
      axis(2,las=1)
    }
    if(metric == names(expectedvals)[length(expectedvals)]){
      axis(1,labels=FALSE)
      text(x=(axTicks(side=1)+((max(axTicks(side=1)-min(axTicks(side=1)))*0.1)))[c(1,3,length(axTicks(side=1)))],
           y=par("usr")[3]-2,
           labels=axTicks(side=1)[c(1,3,length(axTicks(side=1)))],srt=45,pos=2,xpd=NA)
      mtext(text=x.labels[param],side=1,line=3,cex=.8)
    }
  }
}
dev.off()

###Look for best parameter sets = those that match most goodness of fit measures
#does each expected value fall within the quantile?
isMatch = sapply(names(expectedvals) ,
                 function(x){
                   if(length(expectedvals[[x]])==1){
                     get(x)["2.5%",]<=expectedvals[[x]] & get(x)["97.5%",]>=expectedvals[[x]]
                   } else{
                     get(x)["2.5%",]<=max(expectedvals[[x]]) &
                       get(x)["97.5%",]>=min(expectedvals[[x]])
                   }
                 } )
#how many metrics fall within  the quantile for each set?
nMatch = apply(isMatch,1,function(x) length(which(x)))
#which are the best parameter sets?
rankSets = order(nMatch,decreasing=T)
bestSets = cbind(lhsVals[rankSets,] , nMatch = nMatch[rankSets], isMatch[rankSets,])

write.csv(bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch']) |  bestSets[,'nMatch']==max(bestSets[,'nMatch'])-1 | bestSets[,'nMatch']==max(bestSets[,'nMatch'])-2),],"bestParameterSets1000_16.csv",row.names = F)
#bestSets<-read.csv("bestParameterSets1000_16.csv")

write.csv(row.names(bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch']) |  bestSets[,'nMatch']==max(bestSets[,'nMatch'])-1 | bestSets[,'nMatch']==max(bestSets[,'nMatch'])-2),]),"bestParameterSetsRowNames_16.csv",row.names = F)

if (max(bestSets[,'nMatch'])==7){
  bestParamsLarge = bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch']) |  bestSets[,'nMatch']==max(bestSets[,'nMatch'])-1 | bestSets[,'nMatch']==max(bestSets[,'nMatch'])-2),]
  bPL.rows<-as.numeric(row.names(bestParamsLarge[order(bestParamsLarge$d),])) #
  color.bPL<-c(rep("dodgerblue3",sum(bestParamsLarge$d==0)),rep("darkorange",sum(bestParamsLarge$d==1)))
}else{
  bestParamsLarge = bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch']) |  bestSets[,'nMatch']==max(bestSets[,'nMatch'])-1),]
  bPL.rows<-as.numeric(row.names(bestParamsLarge[order(bestParamsLarge$d),])) #
  color.bPL<-c(rep("dodgerblue3",sum(bestParamsLarge$d==0)),rep("darkorange",sum(bestParamsLarge$d==1)))
}

##Overlap plot for the large set of best parameters (7,6,5 overlaps)
png("parameterMatching1000_16_bestParamsLargerevised.png",7,10,units="in",res=500)
#Ignore kappa for the figure
fig.params<-c(1:3,5:11)
#colnames(lhsVals) = c("setNum",colnames(lhsVals)[-1])
par(mfrow=c(length(expectedvals),length(fig.params)) , mar = c(0,0,0,0) ,
    oma = c(5,5,0.2,0.2))
for(metric in names(expectedvals)){
  for(param in fig.params){
    ymax<-max(na.omit(get(metric)["97.5%",]))+max(na.omit(get(metric)["97.5%",]))/16
    ymin<-0
    if(metric=="yearsToFullInfection"){
      ylab = "All infected (yrs)"
      ymin = 2
    }
    if(metric=="prevalenceY1"){
      ylab = "Year 1 prevalence"
    }
    if(metric=="propSurviveY2"){
      ylab = "Year 2 survival"
    } 
    if(metric=="sporeCount"){ 
      ymax<-1e15+1e15/15
      ylab = bquote(paste("Environ. ",italic("Pd")))
    }
    if(metric=="nInfectedY2"){
      ylab="Year 2 infections"
    }
    if(metric=="nInfectedY3"){
      ylab="Year 3 infections"
    }
    if(metric=="nInfectedY4"){
      ylab="Year 4 infections"
    }
    if(metric=="nInfectedY5"){
      ylab="Year 5 infections"
    }
    xmin = min(lhsVals[,param])-(max(lhsVals[,param])-min(lhsVals[,param]))/10
    xmax = max(lhsVals[,param])+(max(lhsVals[,param])-min(lhsVals[,param]))/10
    plot(x=lhsVals[,param][bPL.rows] , y=get(metric)["50%",][bPL.rows] , #These are the best parameter sets and ordered by density dependence
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=color.bPL,
         pch=20,axes=FALSE,xpd=TRUE)
    box()
    segments(x0=lhsVals[,param][bPL.rows] , y0 = get(metric)["2.5%",][bPL.rows] , y1 = get(metric)["97.5%",][bPL.rows],col=color.bPL)
    expectedValueMetric = expectedvals[[metric]]
    if(length(expectedValueMetric)==1){
      abline(h = expectedValueMetric , col = "green4")
    } else {
      polygon(c(par()$usr[1],par()$usr[2],par()$usr[2],par()$usr[1]) ,
              c(min(expectedValueMetric),min(expectedValueMetric),
                max(expectedValueMetric),max(expectedValueMetric)),
              border="green4",col=rgb(0,.55,0,.2))
    }
    
    if(param==1){
      mtext(ylab,side=2,line=4,cex=.8)
      axis(2,las=1)
    }
    if(metric == names(expectedvals)[length(expectedvals)]){
      axis(1,labels=FALSE)
      text(x=(axTicks(side=1)+((max(axTicks(side=1)-min(axTicks(side=1)))*0.1)))[c(1,3,length(axTicks(side=1)))],
           y=par("usr")[3]-2,
           labels=axTicks(side=1)[c(1,3,length(axTicks(side=1)))],srt=45,pos=2,xpd=NA)
      mtext(text=x.labels[param],side=1,line=3,cex=.8)
    }
  }
}
dev.off()

if(max(bestSets[,'nMatch'])==7){
  bestParamsSmall = bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch']) |  bestSets[,'nMatch']==max(bestSets[,'nMatch'])-1),]
  bPS.rows<-as.numeric(row.names(bestParamsSmall[order(bestParamsSmall$d),]))
  color.bPS<-c(rep("dodgerblue3",sum(bestParamsSmall$d==0)),rep("darkorange",sum(bestParamsSmall$d==1)))
}else{
  bestParamsSmall = bestSets[which(bestSets[,'nMatch']==max(bestSets[,'nMatch'])),]
  bPS.rows<-as.numeric(row.names(bestParamsSmall[order(bestParamsSmall$d),])) 
  color.bPS<-c(rep("dodgerblue3",sum(bestParamsSmall$d==0)),rep("darkorange",sum(bestParamsSmall$d==1)))
}

##Overlap plot for the large set of best parameters (7,6,5 overlaps)
png("parameterMatching1000_16_bestParamsSmallrevisedT.png",7,10,units="in",res=500)
#Ignore kappa for the figure
fig.params<-c(1:3,5:11)
#colnames(lhsVals) = c("setNum",colnames(lhsVals)[-1])
par(mfrow=c(length(expectedvals),length(fig.params)) , mar = c(0,0,0,0) ,
    oma = c(5,5,0.2,0.2))
for(metric in names(expectedvals)){
  for(param in fig.params){
    ymax<-max(na.omit(get(metric)["97.5%",]))+max(na.omit(get(metric)["97.5%",]))/16
    ymin<-0
    if(metric=="yearsToFullInfection"){
      ylab = "All infected (yrs)"
      ymin = 2
    }
    if(metric=="prevalenceY1"){
      ylab = "Year 1 prevalence"
    }
    if(metric=="propSurviveY2"){
      ylab = "Year 2 survival"
    } 
    if(metric=="sporeCount"){ 
      ymax<-1e15+1e15/15
      ylab = bquote(paste("Environ. ",italic("Pd")))
    }
    if(metric=="nInfectedY2"){
      ylab="Year 2 infections"
    }
    if(metric=="nInfectedY3"){
      ylab="Year 3 infections"
    }
    if(metric=="nInfectedY4"){
      ylab="Year 4 infections"
    }
    if(metric=="nInfectedY5"){
      ylab="Year 5 infections"
    }
    xmin = min(lhsVals[,param])-(max(lhsVals[,param])-min(lhsVals[,param]))/10
    xmax = max(lhsVals[,param])+(max(lhsVals[,param])-min(lhsVals[,param]))/10
    plot(x=lhsVals[,param][bPS.rows] , y=get(metric)["50%",][bPS.rows] , #These are the best parameter sets and ordered by density dependence
         ylim=c(ymin,ymax),xlim=c(xmin,xmax),col=color.bPS,
         pch=20,axes=FALSE,xpd=TRUE)
    box()
    segments(x0=lhsVals[,param][bPS.rows] , y0 = get(metric)["2.5%",][bPS.rows] , y1 = get(metric)["97.5%",][bPS.rows],col=color.bPS)
    expectedValueMetric = expectedvals[[metric]]
    if(length(expectedValueMetric)==1){
      abline(h = expectedValueMetric , col = "green4")
    } else {
      polygon(c(par()$usr[1],par()$usr[2],par()$usr[2],par()$usr[1]) ,
              c(min(expectedValueMetric),min(expectedValueMetric),
                max(expectedValueMetric),max(expectedValueMetric)),
              border="green4",col=rgb(0,0.55,0,.2))
    }
    
    if(param==1){
      mtext(ylab,side=2,line=4,cex=.8)
      axis(2,las=1)
    }
    if(metric == names(expectedvals)[length(expectedvals)]){
      axis(1,labels=FALSE)
      text(x=(axTicks(side=1)+((max(axTicks(side=1)-min(axTicks(side=1)))*0.1)))[c(1,3,length(axTicks(side=1)))],
           y=par("usr")[3]-2,
           labels=axTicks(side=1)[c(1,3,length(axTicks(side=1)))],srt=45,pos=2,xpd=NA)
      mtext(text=x.labels[param],side=1,line=3,cex=.8)
    }
  }
}
dev.off()



fig.params<-c(1:3,5:11)
paramsToPlot = c("epsilon","psi","beta1","beta2","gamma1","gamma2","mu","d","alpha","V")

##Histogram plots to look at the distribution of best parameters
#First is for the large set
png("parameterHistoPlot_bPL_16revised.png",w=7,h=3,units="in",res=300)
par(mfrow=c(2,5),mar=c(3,0,0.5,1),oma=c(0,3.5,0,0),las=1)
for (i in fig.params){
  ymax=325
  if (i %in% c(1,7)){
    ylab="Frequency"
    labels=TRUE
  }else{
    ylab=NA
    labels=FALSE
  }
  hist(bestParamsLarge[,i],breaks=6,ylim=c(0,ymax),#xlim=c(min(lhsVals[,i]),max(lhsVals[,i])),
              main=NA,xlab=colnames(bestParamsLarge)[i],axes=F)
  segments(x0=median(bestParamsLarge[,i]),y0=0,y1=ymax,col="red")
  axis(1)
  mtext(side=1,text=x.labels[i],cex=0.7,line=2)
  axis(2,labels=labels)
  mtext(side=2,text=ylab,las=3,cex=0.7,line=2.5)
}
dev.off()

png("parameterHistoPlot_bPS_16revised.png",w=7,h=3,units="in",res=300)
par(mfrow=c(2,5),mar=c(3,0,0.5,1),oma=c(0,3.5,0,0),las=1)
for (i in fig.params){
  ymax=25
  if (i %in% c(1,7)){
    ylab="Frequency"
    labels=TRUE
  }else{
    ylab=NA
    labels=FALSE
  }
  hist(bestParamsSmall[,i],breaks=6,ylim=c(0,ymax),#xlim=c(min(lhsVals[,i]),max(lhsVals[,i])),
       main=NA,xlab=colnames(bestParamsSmall)[i],axes=F)
  segments(x0=median(bestParamsSmall[,i]),y0=0,y1=ymax,col="red")
  axis(1)
  mtext(side=1,text=x.labels[i],cex=0.7,line=2)
  axis(2,labels=labels)
  mtext(side=2,text=ylab,las=3,cex=0.7,line=2.5)
}
dev.off()


##Correlation of parameter sets
x.labels1<-c(expression(epsilon),expression(psi),expression(beta[1]),
            expression(beta[2]),expression(gamma[1]),expression(gamma[2]),expression(mu),
            "d",expression(alpha),"V")
png("parameterCorPlot1000_bPS_16ddrevised.png",20,20,units="in",res=300)
par(mfrow=c(length(paramsToPlot),length(paramsToPlot)),mar=c(2,2,0.5,0.5))
for(i1 in 1:length(paramsToPlot)){
  i=paramsToPlot[i1]
  for(j1 in 1:length(paramsToPlot)){
    j=paramsToPlot[j1]
    if(i1<j1){
      plot(y=bestParamsSmall[,i],x=bestParamsSmall[,j],
           xlab = NA , ylab = NA , pch=19 , cex = 2)
    }
    if(i1==j1){
      plot(1,1,type="n",axes=F,xlab="",ylab="")
      text(x=1,y=1,x.labels1[i1],cex=3)
      box()
    }
    if(j1<i1){
      plot(1,y=1,type="n",axes=F,xlab="",ylab="")
      text(x=1,y=1,
           round(cor(bestParamsSmall[,c(i,j)])[1,2],2),
           cex=3,font=ifelse(cor.test(bestParamsSmall[,i],bestParamsSmall[,j])[3]<=(0.002),2,1))
      box()
    }
  }
}
dev.off()

png("parameterCorPlot1000_bPL_16revised.png",20,20,units="in",res=300)
par(mfrow=c(length(paramsToPlot),length(paramsToPlot)),mar=c(2,2,0.5,0.5))
for(i1 in 1:length(paramsToPlot)){
  i=paramsToPlot[i1]
  for(j1 in 1:length(paramsToPlot)){
    j=paramsToPlot[j1]
    if(i1<j1){
      plot(y=bestParamsLarge[,i],x=bestParamsLarge[,j],
           xlab = NA , ylab = NA , pch=19 , cex = 2)
    }
    if(i1==j1){
      plot(1,1,type="n",axes=F,xlab="",ylab="")
      text(x=1,y=1,x.labels1[i1],cex=3)
      box()
    }
    if(j1<i1){
      plot(1,y=1,type="n",axes=F,xlab="",ylab="")
      text(x=1,y=1,
           round(cor(bestParamsLarge[,c(i,j)])[1,2],2),
           cex=3,font=ifelse(cor.test(bestParamsLarge[,i],bestParamsLarge[,j])[3]<=(0.005),2,1))
      box()
    }
  }
}
dev.off()

###Plots to show difference between deterministic and stochastic simulations
determ = readRDS("evaluationMetricsDeterm_16.rds")

png(width=7,height=7,units="in",res=300, "Stochasticity effect on survival 16.png")
plot(c(0,1),c(0,1),type="n",las=1,
     xlab = "Survival year after infection (Stochastic)" , ylab = "Survival year after infection (Deterministic)")
s = evaluationMetrics[["propSurviveY2"]]
d = determ[["propSurviveY2"]]
points(s['50%',],d['50%',] , pch = 19 , col = rgb(0,0,0,.5))
segments(x0 = s['2.5%',] , x1 = s['97.5%',] , y0 = d['50%',], col = rgb(0,0,0,.5))
segments(y0 = d['2.5%',] , y1 = d['97.5%',] , x0 = s['50%',], col = rgb(0,0,0,.5))
abline(0,1 , lty=2)
dev.off()

d = determ[["nInfectedY2"]]
s = evaluationMetrics[["nInfectedY2"]]

lims = max(c(d,s))
png(width=7,height=7,units="in",res=300, "Stochasticity effect on cave infection 16.png")
plot(c(0,lims),c(0,lims),type="n",las=1,
     xlab = "Caves infected year 2 (Stochastic)" , ylab = "Caves infected year 2 (Deterministic)")

points(s['50%',],d['50%',] , pch = 19, col = rgb(0,0,0,.5))
segments(x0 = s['2.5%',] , x1 = s['97.5%',] , y0 = d['50%',], col = rgb(0,0,0,.5))
segments(y0 = d['2.5%',] , y1 = d['97.5%',] , x0 = s['50%',], col = rgb(0,0,0,.5))

abline(0,1 , lty=2)
dev.off()

png(width=4,height=4,units="in",res=300, "Stochasticity effect on prevalence 16.png")
par(mar=c(3,3,0.5,0.5),mgp=c(2,0.5,0))
plot(c(0,1),c(0,1),type="n",las=1,
     xlab = "Prevalence (Stochastic)" , ylab = "Prevalence (Deterministic)")
d = determ[["prevalenceY1"]]
s = evaluationMetrics[["prevalenceY1"]]
points(s['50%',],d['50%',] , pch = 19, col = rgb(0,0,0,.5))
segments(x0 = s['2.5%',] , x1 = s['97.5%',] , y0 = d['50%',], col = rgb(0,0,0,.5))
segments(y0 = d['2.5%',] , y1 = d['97.5%',] , x0 = s['50%',], col = rgb(0,0,0,.5))

abline(0,1 , lty=2)
dev.off()


##Gamma experiment to see whether cave quality matters
##Requires new simulations, but only run for best parameter sets
bestParamsSmall = read.csv("bestParameterSets1000_16.csv")
bestParamsSmall<-bestParamsSmall[which(bestParamsSmall$nMatch>5),]

#run model with new gamma1 values
cl = makeCluster(30) #Sets the number of cores
registerDoParallel(cl) #Registers the cores with the server.
options<-list(preschedule=FALSE)

start.time = Sys.time()
set.seed(159)
nGammas = 100
output <- foreach (z = 1:(nGammas*nrow(bestParamsSmall)), .options.snow=options) %dopar% {
  
  bestGammas = bestParamsSmall[,"gamma1"]
  
  gammaMeans = bestGammas
  gammaVars = sd(bestGammas)
  nGammas = 100
  
  gammaCoV = gammaVars/mean(gammaMeans)
  #gammaCoV = 1 #Verant et al.: fungal colony size can vary from 0-1 much over 0-10 C, temps at which LBBs hibernate (Webb et al.)
  
  
  paramIndex = ceiling(z / nGammas) # get the "base" parameter set
  
  #simulate new values of gamma for all caves
  originalGamma = bestParamsSmall[paramIndex , "gamma1"]
  newGammas = rnorm(nCaves , mean = originalGamma , sd = originalGamma*gammaCoV)
  while(any(newGammas<0)){
    newGammas[newGammas<0] = rnorm(length(which(newGammas<0)) ,
                                   mean = originalGamma , sd = originalGamma*gammaCoV)
  }
  
  library(adaptivetau)
  newParms = bestParamsSmall[paramIndex,] #get corresponding parameter values
  allparms = lapply(newParms , function(x){
    rep(x,nCaves)
  })
  allparms[["gamma1"]] = newGammas #assign already-simulated vals of gamma
  
  
  
  wns.sims = list()
  for(i in 1:n.sims){
    wns.sims[[i]] = multi_model(S_0[[i]] , n_years , nInfected , distances , geom_growth , winter_length , allparms , fidelity)
  }
  output = list(parameterSet = paramIndex , gammas = newGammas ,
                otherParms = newParms[-which(names(newParms)=="gamma1")] , results = wns.sims)
  fileName = paste0("gamma_sims16/output_gamma_1000LHS_100sim_", z,".rds")
  
  saveRDS(output,fileName)
}

stopCluster(cl)

end.time = Sys.time()
total.time = end.time - start.time

#saveRDS(output,file='gammaExpt_output_nc28_lhs500_100sims_100gammas_16.rds')

#Read in output files from gamma experiment
output <- NULL
for(i in 1:2000){
  output[[i]] = readRDS(paste("gamma_sims16/output_gamma_1000LHS_100sim_",i,".rds",sep=""))
}
saveRDS(output,file='gammaExpt_output_nc28_lhs500_100sims_100gammas_16.rds')


#Evaluate gamma experiment
length(output) == nrow(lhsVals)
#test whether this is going to work.


allResults = lapply(1:length(output) , function(x){
  results = output[[x]][['results']]
  a = lapply(names(results[[1]]),function(y){
    lapply(1:length(results),function(z){
      new = results[[z]][[y]]
    })
  })
  names(a) = names(results[[1]])
  return(a)
})


#this is what we're getting as our output for each metric
statToCompute <- function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}


#proportion surviving at cave scale in year 1
propSurviveY1 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year1Survival = sapply(which(!is.na(yearInfected) & yearInfected <= length(endPops)) ,
                           function(cave){
                             endPops[[yearInfected[cave]]][cave]/
                               startPops[yearInfected[cave],cave]
                           })
    mean(year1Survival) #This is the average over all the infected caves in the simulation
  })
  statToCompute(propSurvive) #The appropriate quantiles of the mean survival from each set of simulations
})


#proportion surviving at cave scale in year 2
#this is beginning to end of winter, not through entire year
propSurviveY2 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year2Survival = sapply(which(!is.na(yearInfected) & yearInfected <= (length(endPops)-1)) , function(cave){
      endPops[[yearInfected[cave]+1]][cave]/
        startPops[yearInfected[cave]+1,cave]
    })
    mean(year2Survival)
  })
  statToCompute(propSurvive)
})


#proportion surviving at cave scale after year 2
propSurviveY3 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    allSurvival = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-2)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],"[",cave)
      winterStart = cavePops[[sim]][allOtherYears,cave]
      winterEnd/winterStart
    })
    mean(unlist(allSurvival))
  })
  statToCompute(propSurvive)
})


#prevalence at end of winter in year 1
prevalenceY1 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    year1Prev = sapply(which(!is.na(yearInfected)& yearInfected <= length(endPops)) , function(cave){
      sum(endPops[[yearInfected[cave]]][cave,c("I","F1","F2")])/
        sum(endPops[[yearInfected[cave]]][cave,c("S","I","F1","F2")])
    })
    mean(year1Prev)
  })
  statToCompute(prev)
})


#prevalence at end of winter in and after year 2
prevalenceY2 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    allPrev = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-1)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],function(x){
        sum(x[cave,c("S","I","F1","F2")])/sum(x[cave,c("S","I","F1","F2")])
      })
    })
    mean(unlist(allPrev))
  })
  statToCompute(prev)
})


#survival in year 5
survivalY5 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  survive = sapply(1:length(cavePops),function(sim){
    popSize = cavePops[[sim]]
    sum(popSize[5,]) / sum(popSize[1,])
  })
  statToCompute(survive)
})

yearsToFullInfection = sapply(allResults , function(y){
  infected = y[["infected"]]
  survive = sapply(1:length(infected),function(sim){
    infectedSim = infected[[sim]]
    nInfected = rowSums(infectedSim)
    yearsToInfected = which(nInfected==ncol(infectedSim))[1]
    if (is.na(yearsToInfected)) yearsToInfected <- 12
    return (yearsToInfected)
  })
  statToCompute(survive)
})


#spore count
sporeCount = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  spores = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    endE = t(sapply(endPops,function(X) X$E)) #number spores at end of winter in all caves
    infectedSim = infected[[sim]][1:10,] #which caves are infected? (all years)
    infectedE = endE[which(infectedSim)] #number of spores at end of winter in infected caves
    spores = median(infectedE) #mean spore count in infected caves for this simulation
  })
  statToCompute(spores)
})


#new infections in each year
nInfectedY2 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][1,] &
                                 infected[[sim]][2,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY3 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][2,] &
                                 infected[[sim]][3,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY4 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][3,] &
                                 infected[[sim]][4,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY5 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][4,] &
                                 infected[[sim]][5,])) #which caves are infected?
  })
  statToCompute(newInfections)
})


#Total bat population after 5 years
batPopY5<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY5=sapply(1:length(cavePops),function(sim){
    popSumY5<-sum(cavePops[[sim]][5,])
  })
  statToCompute(popY5)
})


#Total bat populations after 10 years
batPopY10<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY10=sapply(1:length(cavePops),function(sim){
    popSumY10<-sum(cavePops[[sim]][10,])
  })
  statToCompute(popY10)
})


#Proportion of caves without bats
propExtinct<-sapply(allResults, function(y){
  winterPops = y[["winterPops"]]
  popSizes = lapply(winterPops,function(x) rowSums(x[[10]][,-1]))
  extinctionProp = sapply(popSizes, function(x){
    length(which(x==0))/length(x)
  })
  statToCompute(extinctionProp)
})

saveRDS(list(prevalenceY1=prevalenceY1,
             prevalenceY2=prevalenceY2,
             propSurviveY1=propSurviveY1,
             propSurviveY2=propSurviveY2,
             survivalY5=survivalY5,
             yearsToFullInfection=yearsToFullInfection,
             nInfectedY2 = nInfectedY2,
             nInfectedY3 = nInfectedY3,
             nInfectedY4 = nInfectedY4,
             nInfectedY5 = nInfectedY5,
             sporeCount = sporeCount,
             batPopY5 = batPopY5,
             batPopY10 = batPopY10,
             propExtinct = propExtinct),
        "evaluationMetrics1000_GammaExp_16.rds")
rm(output)
gc()


output <- NULL
for(i in 1:2000){
  output[[i]] = readRDS(paste("outputProgress16HighFidelity/wns_output_nc30_lhs1000_100sims",i,".rds",sep=""))
}
saveRDS(output,file='wns_output_nc30_lhs1000_100sims_HighFidelity_16.rds')

length(output) == nrow(lhsVals)
#test whether this is going to work.


allResults = lapply(1:length(output) , function(x){
  results = output[[x]][['results']]
  a = lapply(names(results[[1]]),function(y){
    lapply(1:length(results),function(z){
      new = results[[z]][[y]]
    })
  })
  names(a) = names(results[[1]])
  return(a)
})


#this is what we're getting as our output for each metric
statToCompute <- function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}

#proportion surviving at cave scale in year 1
propSurviveY1 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year1Survival = sapply(which(!is.na(yearInfected) & yearInfected <= length(endPops)) ,
                           function(cave){
                             endPops[[yearInfected[cave]]][cave]/
                               startPops[yearInfected[cave],cave]
                           })
    mean(year1Survival) #This is the average over all the infected caves in the simulation
  })
  statToCompute(propSurvive) #The appropriate quantiles of the mean survival from each set of simulations
})


#proportion surviving at cave scale in year 2
#this is beginning to end of winter, not through entire year
propSurviveY2 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year2Survival = sapply(which(!is.na(yearInfected) & yearInfected <= (length(endPops)-1)) , function(cave){
      endPops[[yearInfected[cave]+1]][cave]/
        startPops[yearInfected[cave]+1,cave]
    })
    mean(year2Survival)
  })
  statToCompute(propSurvive)
})


#proportion surviving at cave scale after year 2
propSurviveY3 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    allSurvival = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-2)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],"[",cave)
      winterStart = cavePops[[sim]][allOtherYears,cave]
      winterEnd/winterStart
    })
    mean(unlist(allSurvival))
  })
  statToCompute(propSurvive)
})


#prevalence at end of winter in year 1
prevalenceY1 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    year1Prev = sapply(which(!is.na(yearInfected)& yearInfected <= length(endPops)) , function(cave){
      sum(endPops[[yearInfected[cave]]][cave,c("I","F1","F2")])/
        sum(endPops[[yearInfected[cave]]][cave,c("S","I","F1","F2")])
    })
    mean(year1Prev)
  })
  statToCompute(prev)
})


#prevalence at end of winter in and after year 2
prevalenceY2 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    allPrev = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-1)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],function(x){
        sum(x[cave,c("S","I","F1","F2")])/sum(x[cave,c("S","I","F1","F2")])
      })
    })
    mean(unlist(allPrev))
  })
  statToCompute(prev)
})


#survival in year 5
survivalY5 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  survive = sapply(1:length(cavePops),function(sim){
    popSize = cavePops[[sim]]
    sum(popSize[5,]) / sum(popSize[1,])
  })
  statToCompute(survive)
})

yearsToFullInfection = sapply(allResults , function(y){
  infected = y[["infected"]]
  survive = sapply(1:length(infected),function(sim){
    infectedSim = infected[[sim]]
    nInfected = rowSums(infectedSim)
    yearsToInfected = which(nInfected==ncol(infectedSim))[1]
    if (is.na(yearsToInfected)) yearsToInfected <- 12
    return (yearsToInfected)
  })
  statToCompute(survive)
})


#spore count
sporeCount = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  spores = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    endE = t(sapply(endPops,function(X) X$E)) #number spores at end of winter in all caves
    infectedSim = infected[[sim]][1:10,] #which caves are infected? (all years)
    infectedE = endE[which(infectedSim)] #number of spores at end of winter in infected caves
    spores = median(infectedE) #mean spore count in infected caves for this simulation
  })
  statToCompute(spores)
})


#new infections in each year
nInfectedY2 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][1,] &
                                 infected[[sim]][2,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#newinfections in each year
nInfectedY3 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][2,] &
                                 infected[[sim]][3,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY4 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][3,] &
                                 infected[[sim]][4,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY5 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][4,] &
                                 infected[[sim]][5,])) #which caves are infected?
  })
  statToCompute(newInfections)
})


#Total bat population after 5 years
batPopY5<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY5=sapply(1:length(cavePops),function(sim){
    popSumY5<-sum(cavePops[[sim]][5,])
  })
  statToCompute(popY5)
})


#Total bat populations after 10 years
batPopY10<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY10=sapply(1:length(cavePops),function(sim){
    popSumY10<-sum(cavePops[[sim]][10,])
  })
  statToCompute(popY10)
})


#Proportion of caves without bats
propExtinct<-sapply(allResults, function(y){
  winterPops = y[["winterPops"]]
  popSizes = lapply(winterPops,function(x) rowSums(x[[10]][,-1]))
  extinctionProp = sapply(popSizes, function(x){
    length(which(x==0))/length(x)
  })
  statToCompute(extinctionProp)
})

saveRDS(list(prevalenceY1=prevalenceY1,
             prevalenceY2=prevalenceY2,
             propSurviveY1=propSurviveY1,
             propSurviveY2=propSurviveY2,
             survivalY5=survivalY5,
             yearsToFullInfection=yearsToFullInfection,
             nInfectedY2 = nInfectedY2,
             nInfectedY3 = nInfectedY3,
             nInfectedY4 = nInfectedY4,
             nInfectedY5 = nInfectedY5,
             sporeCount = sporeCount,
             batPopY5 = batPopY5,
             batPopY10 = batPopY10,
             propExtinct = propExtinct),
        "evaluationMetrics1000_HighFidelity_16.rds")
rm(output)
gc()


output <- NULL
for(i in 1:2000){
  output[[i]] = readRDS(paste("outputProgress16LowFidelity/wns_output_nc30_lhs1000_100sims",i,".rds",sep=""))
}
saveRDS(output,file='wns_output_nc30_lhs1000_100sims_LowFidelity_16.rds')

length(output) == nrow(lhsVals)
#test whether this is going to work.

allResults = lapply(1:length(output) , function(x){
  results = output[[x]][['results']]
  a = lapply(names(results[[1]]),function(y){
    lapply(1:length(results),function(z){
      new = results[[z]][[y]]
    })
  })
  names(a) = names(results[[1]])
  return(a)
})


#this is what we're getting as our output for each metric
statToCompute <- function(x){
  quantile(x,c(0.025,0.5,0.975),na.rm=T)
}

#proportion surviving at cave scale in year 1
propSurviveY1 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year1Survival = sapply(which(!is.na(yearInfected) & yearInfected <= length(endPops)) ,
                           function(cave){
                             endPops[[yearInfected[cave]]][cave]/
                               startPops[yearInfected[cave],cave]
                           })
    mean(year1Survival) #This is the average over all the infected caves in the simulation
  })
  statToCompute(propSurvive) #The appropriate quantiles of the mean survival from each set of simulations
})


#proportion surviving at cave scale in year 2
#this is beginning to end of winter, not through entire year
propSurviveY2 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    year2Survival = sapply(which(!is.na(yearInfected) & yearInfected <= (length(endPops)-1)) , function(cave){
      endPops[[yearInfected[cave]+1]][cave]/
        startPops[yearInfected[cave]+1,cave]
    })
    mean(year2Survival)
  })
  statToCompute(propSurvive)
})



#proportion surviving at cave scale after year 2
propSurviveY3 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  propSurvive = sapply(1:length(cavePops),function(sim){
    startPops = cavePops[[sim]]
    endPops = winterPops[[sim]]
    endPops = lapply(endPops , function(x) apply(x[,c("S","I","F1","F2")],1,sum))
    yearInfected = firstInfection[[sim]]
    allSurvival = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-2)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],"[",cave)
      winterStart = cavePops[[sim]][allOtherYears,cave]
      winterEnd/winterStart
    })
    mean(unlist(allSurvival))
  })
  statToCompute(propSurvive)
})


#prevalence at end of winter in year 1
prevalenceY1 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    year1Prev = sapply(which(!is.na(yearInfected)& yearInfected <= length(endPops)) , function(cave){
      sum(endPops[[yearInfected[cave]]][cave,c("I","F1","F2")])/
        sum(endPops[[yearInfected[cave]]][cave,c("S","I","F1","F2")])
    })
    mean(year1Prev)
  })
  statToCompute(prev)
})


#prevalence at end of winter in and after year 2
prevalenceY2 = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  firstInfection = lapply(infected,function(x){
    apply(x,2,function(cave) which(cave)[1])
  })
  prev = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    yearInfected = firstInfection[[sim]]
    allPrev = sapply(which(!is.na(yearInfected)& yearInfected <= (length(endPops)-1)) , function(cave){
      allOtherYears  = (yearInfected[cave]+2):length(endPops)
      winterEnd = sapply(endPops[allOtherYears],function(x){
        sum(x[cave,c("S","I","F1","F2")])/sum(x[cave,c("S","I","F1","F2")])
      })
    })
    mean(unlist(allPrev))
  })
  statToCompute(prev)
})


#survival in year 5
survivalY5 = sapply(allResults , function(y){
  cavePops = y[["cavePops"]]
  survive = sapply(1:length(cavePops),function(sim){
    popSize = cavePops[[sim]]
    sum(popSize[5,]) / sum(popSize[1,])
  })
  statToCompute(survive)
})


#Years until all caves infected
yearsToFullInfection = sapply(allResults , function(y){
  infected = y[["infected"]]
  survive = sapply(1:length(infected),function(sim){
    infectedSim = infected[[sim]]
    nInfected = rowSums(infectedSim)
    yearsToInfected = which(nInfected==ncol(infectedSim))[1]
    if (is.na(yearsToInfected)) yearsToInfected <- 12
    return (yearsToInfected)
  })
  statToCompute(survive)
})


#spore count
sporeCount = sapply(allResults , function(y){
  winterPops = y[["winterPops"]]
  infected = y[["infected"]]
  spores = sapply(1:length(winterPops),function(sim){
    endPops = winterPops[[sim]]
    endE = t(sapply(endPops,function(X) X$E)) #number spores at end of winter in all caves
    infectedSim = infected[[sim]][1:10,] #which caves are infected? (all years)
    infectedE = endE[which(infectedSim)] #number of spores at end of winter in infected caves
    spores = median(infectedE) #mean spore count in infected caves for this simulation
  })
  statToCompute(spores)
})


#new infections in each year
nInfectedY2 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][1,] &
                                 infected[[sim]][2,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY3 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][2,] &
                                 infected[[sim]][3,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY4 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][3,] &
                                 infected[[sim]][4,])) #which caves are infected?
  })
  statToCompute(newInfections)
})

#new infections in each year
nInfectedY5 = sapply(allResults , function(y){
  infected = y[["infected"]]
  newInfections = sapply(1:length(infected),function(sim){
    infectedSim = length(which(!infected[[sim]][4,] &
                                 infected[[sim]][5,])) #which caves are infected?
  })
  statToCompute(newInfections)
})


#Total bat population after 5 years
batPopY5<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY5=sapply(1:length(cavePops),function(sim){
    popSumY5<-sum(cavePops[[sim]][5,])
  })
  statToCompute(popY5)
})


#Total bat populations after 10 years
batPopY10<-sapply(allResults, function(y){
  cavePops=y[["cavePops"]]
  popY10=sapply(1:length(cavePops),function(sim){
    popSumY10<-sum(cavePops[[sim]][10,])
  })
  statToCompute(popY10)
})


#Proportion of caves without bats
propExtinct<-sapply(allResults, function(y){
  winterPops = y[["winterPops"]]
  popSizes = lapply(winterPops,function(x) rowSums(x[[10]][,-1]))
  extinctionProp = sapply(popSizes, function(x){
    length(which(x==0))/length(x)
  })
  statToCompute(extinctionProp)
})

saveRDS(list(prevalenceY1=prevalenceY1,
             prevalenceY2=prevalenceY2,
             propSurviveY1=propSurviveY1,
             propSurviveY2=propSurviveY2,
             survivalY5=survivalY5,
             yearsToFullInfection=yearsToFullInfection,
             nInfectedY2 = nInfectedY2,
             nInfectedY3 = nInfectedY3,
             nInfectedY4 = nInfectedY4,
             nInfectedY5 = nInfectedY5,
             sporeCount = sporeCount,
             batPopY5 = batPopY5,
             batPopY10 = batPopY10,
             propExtinct = propExtinct),
        "evaluationMetrics1000_LowFidelity_16.rds")


####Plots comparing the experiments
bestParamsRowsSmall<-read.csv("bestParameterSetsRowNames_16.csv")[1:dim(bestParamsSmall)[1],]


evaluationMetricsBase<-readRDS("evaluationMetrics1000_16.rds")
evaluationMetricsGamma<-readRDS("evaluationMetrics1000_GammaExp_16.rds") #This structure differs, it is 100 random simulations using the mean and variance of the best parameter sets that are referred in the others
evaluationMetricsLowFid<-readRDS("evaluationMetrics1000_LowFidelity_16.rds")
evaluationMetricsHighFid<-readRDS("evaluationMetrics1000_HighFidelity_16.rds")

out.cols<-which(names(evaluationMetricsBase)=="batPopY10" | names(evaluationMetricsBase)=="propExtinct")

#Plot distribution of median values from simulations using best parameter sets (see note about Gammas above)
labels<-c("Final Population Size","Cave extinction")
png(width=4,height=4,units="in",res=400, "CompareExperiments_16revised.png")
par(mfrow=c(1,length(out.cols)),las=1,mar=c(0,4.25,0,0),oma=c(6,0,0.3,0.3))
for (i in 1:length(out.cols)){
  if (i==1){
    data<-rbind(cbind(log10(evaluationMetricsBase[[out.cols[i]]][2,bestParamsRowsSmall]+1),0),
                cbind(log10(evaluationMetricsLowFid[[out.cols[i]]][2,bestParamsRowsSmall]+1),1),
                cbind(log10(evaluationMetricsHighFid[[out.cols[i]]][2,bestParamsRowsSmall]+1),2),
                cbind(log10(evaluationMetricsGamma[[out.cols[i]]][2,]+1),3))
  }else{
    data<-rbind(cbind(evaluationMetricsBase[[out.cols[i]]][2,bestParamsRowsSmall],0),
                cbind(evaluationMetricsLowFid[[out.cols[i]]][2,bestParamsRowsSmall],1),
                cbind(evaluationMetricsHighFid[[out.cols[i]]][2,bestParamsRowsSmall],2),
                cbind(evaluationMetricsGamma[[out.cols[i]]][2,],3))
  }
  boxplot(data[,1]~data[,2],xlab=NA,
          col=c("dodgerblue4","dodgerblue3","skyblue3","skyblue1"),xaxt="n")
  axis(1,labels=c("Base Model","Low Fidel.","High Fidel.","Hibern. Qual."),at=c(1,2,3,4),las=2)
  #mtext(line=4,side=1,text=labels[i])
  if (i==1) mtext(line=2.5,side=2,text=expression("log"[10]*"(Final population size)"),las=0)
  if (i==2) mtext(line=2.5,side=2,text="Proportion colonies extinct",las=0)
}
dev.off()


##Sample trajectory plots to show what output looks like
evaluationMetrics<-readRDS('evaluationMetrics1000_16.rds')
bestResults = allResults[bestParamsRowsSmall]
n.dens.dep<-length(which(bestParamsRowsSmall<=(length(allResults)/2)))
n.freq.dep<-length(which(bestParamsRowsSmall>(length(allResults)/2)))

#colors<-c(rainbow(6,alpha=0.1)[c(2,3,5,6)])
library(RColorBrewer)
#palette<-brewer.pal(8,"Set2")
#rgbs<-col2rgb(c("yellow2","forestgreen","royalblue3","slateblue3"))
#rgbs<-col2rgb(palette[c(3,1)])
palette<-c("dodgerblue3","darkorange")
rgbs<-col2rgb(palette[c(1,2)])
colors<-rgb(rgbs[1,],rgbs[2,],rgbs[3,],max=255,alpha=30)

#Population trajectories
png(width=6,height=3.5,units="in",res=400,"Trajectories16revised.png")
par(las=1,mar=c(3,4.5,0.5,0.5),mgp=c(3,0.5,0))
plot(c(0,10),c(0,max(starting_metapops)),xlab=NA,ylab=NA,type="n",yaxt="n")
axis(side=2,at=c(0,5e5,1e6,1.5e6),labels=expression(0,5%*%10^5,1%*%10^6,1.5%*%10^6))
for(set in 1:length(bestResults)){
  cavePops = bestResults[[set]][["cavePops"]]
  popSizes = lapply(cavePops,rowSums)
  if (set <= n.dens.dep) color<-colors[1]
  if (set > n.dens.dep) color<-colors[2]
  for(j in 1:length(popSizes)){
    lines(popSizes[[j]]~c(0:10) , col = color)
  }
}
legend("topright",col=palette[c(1,2)],legend=c("Dens. dep. matches","Freq. dep. matches"),lwd=1.5)
mtext(side=1,line=2,"Year")
mtext(side=2,line=3.5,"Metapopulation size",las=3)
dev.off()

#Estimate the percent population decline for plausible parameter sets
#Requires the full output as above
median_decline<-NULL
extinct<-rep(FALSE,42)
for(set in 1:length(bestResults)){
  metapops<-sapply(bestResults[[set]]["cavePops"][[1]],FUN=rowSums)
  median_decline[set]<-1-median(metapops[10,]/metapops[1,])
  if (median_decline[set]>0.9999) extinct[set]<-TRUE
}
#

#The proportion of caves extinct
#rgbs<-col2rgb(palette[c(3,1)])
colors<-rgb(rgbs[1,],rgbs[2,],rgbs[3,],max=255,alpha=70)
png(width=6,height=3.5,units="in",res=400,"CavesExtinct16revised.png")
par(las=1,mar=c(3,4,0.5,0.5),mgp=c(3,0.5,0))
plot(c(0,1),c(0,20),xlab=NA,ylab="Density",type="n")
for(set in 1:length(bestResults)){
  winterPops = bestResults[[set]][["winterPops"]]
  popSizes = lapply(winterPops,function(x) rowSums(x[[10]][,-1]))
  extinctionProp = sapply(popSizes, function(x) length(which(x==0))/length(x))
  if (set <= n.dens.dep) color<-colors[1]
  if (set > n.dens.dep) color<-colors[2]
  lines(density(extinctionProp,bw=0.02),xlim=c(0,1) , col = color)
}
legend("topleft",col=palette[c(1,2)],legend=c("Dens. dep. matches","Freq. dep. matches"),lwd=1.5)
mtext(side=1,line=2,"Proportion of colonies extinct")
dev.off()


#The number of new cave infections in each year
#rgbs<-col2rgb(palette[c(3,1)])
colors<-rgb(rgbs[1,],rgbs[2,],rgbs[3,],max=255,alpha=30)
png(width=6,height=3.5,units="in",res=400,"CavesInfect16revised.png")
par(las=1,mar=c(3,3,0.5,0.5),mgp=c(2,0.5,0))
plot(c(0,10),c(0,35),xlab=NA,ylab="Number of new infections",type="n")
for(set in 1:length(bestResults)){
  infections = bestResults[[set]][["infected"]]
  if (set <= n.dens.dep) color<-colors[1]
  if (set > n.dens.dep) color<-colors[2]
  for(i in 1:length(infections)){
    newInfections = diff(rowSums(infections[[i]]))
    newInfections = ifelse(newInfections<0,0,newInfections) #This makes new infections zero when there are more recovered caves than newly infected (i.e. negative infections)
    newInfections = c(4,newInfections)
    lines(newInfections ~ c(0:10) , col = color)
  }
}
legend("topright",col=palette[c(1,2)],legend=c("Dens. dep. matches","Freq. dep. matches"),lwd=1.5)
mtext(side=1,line=2,"Year")
dev.off()
