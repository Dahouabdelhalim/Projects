#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#         **Details**         #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# use alt + o to collapse sections


#-#-# AUTHOR
# J K Wilson-Aggarwal


#-#-# DATE
# 26/06/19


#-#-# LIBRARIES
# epimodel v1.3.0 
# igraph v1.2.2
# network v1.13.0.1
# networkdynamic v0.9.0
# intergraph v2.0.2


#-#-# DESCRIPTION
# This code was used to run simulations for the disease model described in:

# Wilson-Aggarwal et al. 2019. High-resolution contact networks of free-ranging
# domestic dogs Canis familiaris and implications for transmission of infection


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#      **Disease model**      #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Epimodel provides a framework/code for simulating epidemics on dynamic networks and for recording
# the results at each time step. The package is very flexible and provides functions/modules
# that are initiated at each time step e.g. infection module, progress module...

# We will modify the code in the epimodel modules to simulate disease through time on a static
# binary network, weighted binary network and random networks. The model will use equations and
# parameters outlined in the paper.


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#              |              #----
#           Modules           #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

library("EpiModel")

# tutorial for setting up modules for modelling transmission through the observed network:
# http://www.epimodel.org/tut.html


#-#-#-#-#-#-# INITIALISATION MODULE

# This module is only called at time point 1 and sets up the model.
# The default module does lots of work, such as simulating an initial network 
# from the model fit, that we do not need to do with an observed network. Here, 
# there are three key steps for the initialization module: 
# - set up the master data list (dat) including its core data structures
# - initialize infection among the nodes in the network
# - and use the get_prev.net function to record summary epidemiological statistics.

# we will set x, param, init and control before running the models.
# x = network dynamic object
# param = object with parameter details (infection probability and act rate)
# control = object with module settings
# s = simulation number (which is recorded in an internal module of epimodel)
# seed = a modification to keep track of the seeded individual during the simulations

# getAnywhere(initialize.net) # get the code for the default initialization function


# initialisation module modified for seeded individuals & a static network
net.init <- function(x, param, init, control, s, seed = i) {
  
  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control
  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()
  
  # Network Parameters
  nw <- activate.vertices(x, onset = 1, terminus = Inf)
  
  dat$nw <- nw
  dat$param$modes <- 1
  
  # Initialization
  
  ## Infection Status and Time Modules
  n <- network.size(dat$nw)
  dat$attr$status <- rep("s", n) # seed susceptible
  #dat$attr$status[sample(1:n, init$i.num)] <- "i" # seed random infectious
  dat$attr$status[seed] <- "i" # set seed of infectious
  
  dat$attr$active <- rep(1, n)
  dat$attr$entrTime <- rep(1, n)
  dat$attr$exitTime <- rep(NA, n)
  
  dat$attr$expTime <- rep(0, n)
  dat$attr$expPeriod <- rep(0, n)
  dat$attr$infTime <- rep(0, n)
  dat$attr$infPeriod <- rep(0, n)
  
  # assign parameters for infectious individual
  dat$attr$infTime[dat$attr$status == "i"] <- 2 # as timestep one is used to run the initialisation module
  dat$attr$expDur[dat$attr$status == "i"] <- round(rgamma(n = init$i.num, shape = 1.1, scale = 20.1)) # incubation time
  dat$attr$infDur[dat$attr$status == "i"] <- round(rgamma(n = init$i.num, shape = 3.0, scale = 0.9)) # infectious time
  
  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)
  
  dat$epi$se.flow[1] <- 0
  dat$epi$ei.flow[1] <- 0  # number of infected in the time step
  dat$epi$ir.flow[1] <- 0  # number of recovered/removed in the time step
  dat$epi$e.num[1] <- 0 # total number of exposed
  dat$epi$r.num[1] <- 0 # total number of removed/recovered
  
  dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  
  
  return(dat)
}


# further modification for random networks
net.init.random <- function(x, param, init, control, s, seed = i) {
  
  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control
  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()
  
  
  # check that s is automatically updated at the end of each simulation
  #print(paste0("sim = ", s))
  #Sys.sleep(5)
  
  # get random network s
  nw <- x[[s]]
  
  # Network Parameters
  nw <- activate.vertices(nw, onset = 1, terminus = Inf)
  
  dat$nw <- nw
  dat$param$modes <- 1
  
  # Initialization
  
  ## Infection Status and Time Modules
  n <- network.size(dat$nw)
  dat$attr$status <- rep("s", n) # seed susceptible
  #dat$attr$status[sample(1:n, init$i.num)] <- "i" # seed random infectious
  dat$attr$status[seed] <- "i" # set seed of infectious
  
  dat$attr$active <- rep(1, n)
  dat$attr$entrTime <- rep(1, n)
  dat$attr$exitTime <- rep(NA, n)
  
  dat$attr$expTime <- rep(0, n)
  dat$attr$expPeriod <- rep(0, n)
  dat$attr$infTime <- rep(0, n)
  dat$attr$infPeriod <- rep(0, n)
  
  # randomly seed infectious individuals
  dat$attr$infTime[dat$attr$status == "i"] <- 2
  dat$attr$expDur[dat$attr$status == "i"] <- round(rgamma(n = init$i.num, shape = 1.1, scale = 20.1)) # incubation time
  dat$attr$infDur[dat$attr$status == "i"] <- round(rgamma(n = init$i.num, shape = 3.0, scale = 0.9)) # infectious time
  
  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)
  
  dat$epi$se.flow[1] <- 0
  dat$epi$ei.flow[1] <- 0  # number of infected in the time step
  dat$epi$ir.flow[1] <- 0  # number of recovered/removed in the time step
  dat$epi$e.num[1] <- 0 # total number of exposed
  dat$epi$r.num[1] <- 0 # total number of removed/recovered
  
  dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  
  
  return(dat)
}


#-#-#-#-#-#-# INFECTION MODULE 

# This module determines who is infected at each time step and
# determines the incubation period and infectious period for 
# the newly infected.

# getAnywhere(infection.net) # original code in epimodel package 


# modified infection module using parameters and equations from paper
infect <- function (dat, at) {
  
  ## Variables ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate
  
  nw <- dat$nw
  tea.status <- dat$control$tea.status
  
  # Vector of infected and susceptible IDs
  idsSus <- which(active == 1 & status == "s") 
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # set initial TEAs
  if(at == 2 & tea.status == TRUE){
    
    nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "s", onset = 1, terminus = Inf, v = idsSus)
    nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "i", onset = 1, terminus = Inf, v = idsInf)
    
  }
  
  
  # Initialize vectors
  nExp <- 0
  
  
  ## Processes ##
  # If some infected AND some susceptible then ...
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, at)
    
    # If some discordant edges then ...
    if (!(is.null(del))) {
      
      # calculate final probability
      del$transProb <- inf.prob # Infection probabilities
      del$act.rate <- act.rate # act rates
      del$finalProb <- 1 - (1 - del$transProb)^del$act.rate # final probability
      
      
      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb) # DECIDES IF THE VIRUS IS TRANSMITTED
      del <- del[which(transmit == 1), ]
      
      # Set new infections vector
      idsNewInf <- unique(del$sus) # IDs of newly infected
      nExp <- length(idsNewInf) # n newly infected
      
      
      # Update attributes
      if (nExp > 0) {
        
        if (tea.status == TRUE) {
          
          nw <- activate.vertex.attribute(nw, prefix = "testatus", 
                                          value = "e", onset = at, terminus = Inf, 
                                          v = idsNewInf)
        }
        
        dat$attr$status[idsNewInf] <- "e" # change from susceptible to exposed
        dat$attr$expTime[idsNewInf] <- at # add time of infection
        dat$attr$expPeriod[idsNewInf] <- round(rgamma(n = nExp, shape = 1.1, scale = 20.1)) # assign incubation time
        dat$attr$infPeriod[idsNewInf] <- round(rgamma(n = nExp, shape = 3.0, scale = 0.9)) # assign infectious time
        
      }
      
      if (any(names(nw$gal) %in% "vertex.pid")) {
        
        del$sus <- get.vertex.pid(nw, del$sus)
        del$exp <- get.vertex.pid(nw, del$exp)
        
      }
    }
  }
  
  if (nExp > 0) {
    
    del <- del[!duplicated(del$sus), ]
    
    if (at == 2) {
      
      dat$stats$transmat <- del
    } else {
      
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
      
    }
  }
  
  if (at == 2) {
    
    dat$epi$se.flow[at] <- 0
    
  } else {
    
    dat$epi$se.flow[at] <- nExp
    
  }
  
  dat$nw <- nw
  #dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  
  return(dat)
}


# infection module further modified for weights
infect.weighted <- function (dat, at) {
  
  ## Variables ##
  active <- dat$attr$active
  status <- dat$attr$status
  
  inf.prob <- dat$param$inf.prob
  act.rate <- dat$param$act.rate
  
  nw <- dat$nw
  tea.status <- dat$control$tea.status
  
  # Vector of infected and susceptible IDs
  idsSus <- which(active == 1 & status == "s") 
  idsInf <- which(active == 1 & status == "i")
  nActive <- sum(active == 1)
  nElig <- length(idsInf)
  
  # set initial TEAs
  if(at == 2 & tea.status == TRUE){
    
    nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "s", onset = 1, terminus = Inf, v = idsSus)
    nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "i", onset = 1, terminus = Inf, v = idsInf)
    
  }
  
  
  # Initialize vectors
  nExp <- 0
  
  
  ## Processes ##
  # If some infected AND some susceptible then ...
  if (nElig > 0 && nElig < nActive) {
    
    # Get discordant edgelist
    del <- discord_edgelist(dat, idsInf, idsSus, at)
    
    
    # If some discordant edges then ...
    if (!(is.null(del))) {
      
      # for each dyad get the edge weight
      for(i in 1:nrow(del)){
        del$w[i] <- unlist(get.edges(x = net, v = del$inf[i], alter = del$sus[i], na.omit = T))[4]
      } 
      
      # calculate final probability
      del$transProb <- inf.prob # Infection probabilities
      del$act.rate <- act.rate*(del$w / (1 + del$w)) # weighted act rates
      del$finalProb <- (1 - (1 - del$transProb)^del$act.rate)*2 # final probability multiplied by 2 to make mean the same as the unweighted sims
      
      # Randomize transmissions and subset df
      transmit <- rbinom(nrow(del), 1, del$finalProb) # DECIDES IF THE VIRUS IS TRANSMITTED
      del <- del[which(transmit == 1), ]
      
      # Set new infections vector
      idsNewInf <- unique(del$sus)
      nExp <- length(idsNewInf)
      
      
      # Update attributes
      if (nExp > 0) {
        
        if (tea.status == TRUE) {
          
          nw <- activate.vertex.attribute(nw, prefix = "testatus", 
                                          value = "e", onset = at, terminus = Inf, 
                                          v = idsNewInf)
        }
        
        dat$attr$status[idsNewInf] <- "e" # change s to e
        dat$attr$expTime[idsNewInf] <- at # changed infTime to expTime
        dat$attr$expPeriod[idsNewInf] <- round(rgamma(n = nExp, shape = 1.1, scale = 20.1)) # incubation time
        dat$attr$infPeriod[idsNewInf] <- round(rgamma(n = nExp, shape = 3.0, scale = 0.9)) # infectious time
        
      }
      
      if (any(names(nw$gal) %in% "vertex.pid")) {
        
        del$sus <- get.vertex.pid(nw, del$sus)
        del$exp <- get.vertex.pid(nw, del$exp)
        
      }
    }
  }
  
  if (nExp > 0) {
    
    del <- del[!duplicated(del$sus), ]
    
    if (at == 2) {
      
      dat$stats$transmat <- del
    } else {
      
      dat$stats$transmat <- rbind(dat$stats$transmat, del)
      
    }
  }
  
  if (at == 2) {
    
    dat$epi$se.flow[at] <- 0
    
  } else {
    
    dat$epi$se.flow[at] <- nExp
    
  }
  
  dat$nw <- nw
  #dat$nw <- set.vertex.attribute(dat$nw, "status", dat$attr$status)
  
  return(dat)
}


#-#-#-#-#-#-# PROGRESS MODULE

# This module tracks and records changes at each time point.

# getAnywhere(recovery.net) # original code in epimodel package 


# modified progress module
progress <- function (dat, at) {
  
  ## VARIABLES ##
  tea.status <- dat$control$tea.status
  nw <- dat$nw
  
  active <- dat$attr$active
  status <- dat$attr$status 
  infTime <- dat$attr$infTime
  infPeriod <- dat$attr$infPeriod
  
  infDur <- ifelse(active == 1 & status == "i", at - infTime, 0)
  
  
  # INFECTED -> REMOVED 
  nRec <- 0
  
  vecDeaths <- which(active == 1 & status == "i" & infDur > infPeriod) 
  nRec <- length(vecDeaths) 
  
  if (nRec > 0){
    
    dat$attr$active[vecDeaths] <- 0
    dat$attr$status[vecDeaths] <- "r"
    dat$attr$exitTime[vecDeaths] <- at 
    
    
    if (tea.status == TRUE) {
      
      nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "r", 
                                      onset = at, terminus = Inf, v = vecDeaths)
      
      nw <- deactivate.vertices(nw, onset = at, terminus = Inf, 
                                v = vecDeaths, deactivate.edges = TRUE) 
    }
    
    
  }
  
  
  # update parameters
  active <- dat$attr$active
  status <- dat$attr$status
  
  
  # EXPOSED -> INFECTIOUS
  nInf <- 0
  expTime <- dat$attr$expTime
  expPeriod <- dat$attr$expPeriod
  
  expDur <- ifelse(active == 1 & status == "e", at - expTime, 0) # exposed to incubation rate
  
  vecInf <- which(active == 1 & status == "e" & expDur > expPeriod) # vector of infectious ind
  nInf <- length(vecInf)
  
  if (nInf > 0) {
    
    dat$attr$status[vecInf] <- "i"
    dat$attr$infTime[vecInf] <- at
    
    if (tea.status == TRUE) {
      
      nw <- activate.vertex.attribute(nw, prefix = "testatus", value = "i", 
                                      onset = at, terminus = Inf, v = vecInf)
    }
    
  }
  
  
  # update parameters
  active <- dat$attr$active
  status <- dat$attr$status
  
  # summary statistics - FINISH
  dat$epi$ei.flow[at] <- nInf  # number of new infected in the time step
  dat$epi$ir.flow[at] <- nRec  # number of new recovered/removed in the time step
  dat$epi$e.num[at] <- sum(active == 1 & status == "e") # total number of exposed
  dat$epi$i.num[at] <- sum(active == 1 & status == "i") # total number of infectious
  dat$epi$s.num[at] <- sum(active == 1 & status == "s") # total number of susceptible
  dat$epi$r.num[at] <- sum(active == 0 & status == "r") # total number of removed/recovered
  
  #nw <- set.vertex.attribute(nw, "status", status)
  dat$nw <- nw
  
  return(dat)
} 



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#          Act rate           #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# calculate and set the act rate for simulations. See paper for specific details and references.

#b <- 1.8 * (3.1)^-1 * mean(degree)/(mean(degree^2)-mean(degree))
#a <- log(1-b)/log(1-0.49)
#rm(a, b)

# act rate when R0 = 1.2: 
# Kakale 
# act = 0.052
# Magrao 
# act = 0.080

# act rate when R0 = 1.8:
# Kakale 
# act = 0.123
# Magrao
# act  = 0.132

# act rate when R0 = 2.4:
# Kakale
# act = 0.166
# Magrao
# act = 0.179


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#        Set up networks      #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-#-#-#-#-#-#-#-# CREATE OBSERVED DYNAMIC NETWORK OBJECT

# N.B. code for adding vertex attributes from meta data is not included

library("igraph")
library("intergraph")
library("network")
library("networkDynamic")


matrix <- as.matrix(read.csv("network.csv", row.names = 1))# import static network
graph <- graph.adjacency(matrix,mode="undirected", weighted=T, diag=FALSE)
net <- asNetwork(graph)

detach("package:igraph")
detach("package:intergraph")


# create networkDynamic object from static network
DYnet <- list() # list object


# create list with the static networks replicated for n time steps
for(i in 1:300){
  
  DYnet[[i]] <- net # set static network as time step 1
  
}

# turn network lists into a network dynamic object
DYNET <- networkDynamic(network.list = DYnet, create.TEAs = F)

set.edge.attribute(DYNET, "weight", net%e%"weight")

n.ind <- network.size(net)
n.edges <- network.edgecount(net)



#-#-#-#-#-#-#-#-#-# CREATE RANDOM NETWORKS

# Generate 100 random networks and for each create a networkdynamic object with 300 time steps 

library("EpiModel")
library("networkDynamic")
library("intergraph")
library("igraph")


DYNET.r <- list() # list object

set.seed(3157)

for(i in 1:100){
  
  DYnet.r.ind <- list() # list object
  
  rnet <- erdos.renyi.game(n.ind, n.edges, type = "gnm") # random graph with same no. of nodes and edges as observed
  rnet <- as.matrix(get.adjacency(rnet, sparse = F)) # convert igraph object to matrix
  rnet <- as.network(rnet, directed = F) # convert matrix to a network object
  
  for(j in 1:300){
    
    DYnet.r.ind[[j]] <- rnet # set static network as time step j
    
  }
  
  DYNET.r[[i]] <- networkDynamic(network.list = DYnet.r.ind, create.TEAs = F)
  
}
rm(rnet, i, j, DYnet.r.ind)


detach("package:igraph")
detach("package:intergraph")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#    Binomial - Simulations   #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

library("EpiModel")
library("networkDynamic")

#-#-#-#-#-#-#-#-#-# SET UP MODEL

sims.B <- list() # list to store results for simulations
  
init <- init.net(i.num = 1) # number of initially infected individuals
  
param <- param.net(inf.prob = 0.49, # infection probability
                    act.rate = act) # number of acts that occur within a partnership each time unit

control <- control.net(type = "SI",
                       nsims = 100, nsteps = 300, 
                       #nwstats.formula = formation, 
                       #delete.nodes = TRUE,
                       initialize.FUN = net.init,
                       recovery.FUN = NULL,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       module.order = c("infection.FUN", "progress.FUN", "get_prev.FUN"),
                       skip.check = TRUE,
                       #depend = FALSE,
                       tea.status=TRUE,
                       save.nwstats = FALSE, 
                       save.network = TRUE)
  

#-#-#-#-#-#-#-#-#-# SIMULATIONS
set.seed(7551)
for(i in 1:n.ind){
  
  sims.B[[i]] <- netsim(DYNET, param, init, control)
  
}

# saveRDS(sim.B, "sim_B.rds") # save r object, change name to reflect settlement and R0


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#    Weighted - Simulations   #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# reset sims
sims.W <- list()

# update control to use infect.weighted module
control <- control.net(type = "SI",
                       nsims = 100, nsteps = 300, 
                       #nwstats.formula = formation, 
                       #delete.nodes = TRUE,
                       initialize.FUN = net.init,
                       recovery.FUN = NULL,
                       infection.FUN = infect.weighted,
                       progress.FUN = progress,
                       module.order = c("infection.FUN", "progress.FUN", "get_prev.FUN"),
                       skip.check = TRUE,
                       #depend = FALSE,
                       tea.status=TRUE,
                       save.nwstats = FALSE, 
                       save.network = TRUE)


#-#-#-#-#-#-#-#-#-# SIMULATIONS

set.seed(7551)

for(i in 1:n.ind){
  
  # run simulations
  sims.W[[i]] <- netsim(DYNET, param, init, control)
  
}

# saveRDS(sims.W, "sims_W.rds") # save r object, change name to reflect settlement and R0


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#     Random - Simulations    #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

#-#-#-#-#-#-#-#-#-# SET UP MODEL

sims.R <- list()

control <- control.net(type = "SI",
                       nsims = 100, nsteps = 300, 
                       #nwstats.formula = formation, 
                       #delete.nodes = TRUE,
                       initialize.FUN = net.init.random,
                       recovery.FUN = NULL,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       module.order = c("infection.FUN", "progress.FUN", "get_prev.FUN"),
                       skip.check = TRUE,
                       #depend = FALSE,
                       tea.status=TRUE,
                       save.nwstats = FALSE, 
                       save.network = TRUE)


#-#-#-#-#-#-#-#-#-# SIMULATIONS
set.seed(759)

for(i in 1:n.ind){
  
  sims.R[[i]] <- netsim(DYNET.r, param, init, control)
  
}


#saveRDS(sims.R, "sims_R.rds") # save r object, change name to reflect settlement and R0



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#    Get metrics from sims    #----
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


# create empty master dataframe
sim.df <- data.frame(id = factor(), x = numeric(), x.hh = numeric(), x.time = numeric())
sims <- sims.B # change to simulation of interest


# for each seeded individual...
for(i in 1:length(sims)){
  
  # get the number removed by end of each sim
  x.r <- unlist(sims[[i]]$epi$r.num[300,]) # epidemic size
  y.r <- rep(i, length(x.r)) # id for individual
  
  # create empty vectors to store info
  x.t <- vector() # time point epidemic ended
  x.hh <- vector() # no. of infections outside of seeded individuals household
  
  index <- which(net%v%"hh" == c(net%v%"hh")[i]) # houseold of seeded individual
  
  # for each simulation ...
  for(j in 1:100){
    
    # get the time step the epidemic stops
    x.t[j] <- which(sims[[i]]$epi$i.num[,j] < 1 & sims[[i]]$epi$e.num[,j] < 1)[1] # which timepoint are there no infected or exposed
    
    # get number of infections outside of household
    transmat <- get_transmat(sims[[i]], sim = j) # data on who was infected in sim
    
    x.hh[j] <- nrow(transmat[!transmat$inf %in% index, ]) # number of infections outside of household

  }
  
  # store data in new dataframe
  new.r <- data.frame(id = y.r, x = x.r, x.hh = x.hh, x.time = x.t)
  
  # rbind new dataframe with master
  sim.df <- rbind(sim.df, new.r) 
  rownames(sim.df) <- seq(length=nrow(sim.df)) # reset row names
  
}
rm(index, transmat, i, j, x.r, y.r, new.r, x.t, x.hh)


# add individual and household IDs to master dataframe
sim.df$hh <- c(net%v%"hh")[sim.df$id]
sim.df$name <- c(net%v%"vertex.names")[sim.df$id]


sims.b.df <- sims.df # store dataframe before running again for another set of simulations

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#----