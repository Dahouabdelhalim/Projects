library(igraph)
library(sna)

### BASIC SIMULATION FUNCTIONS

### function to create a network (matrix n x n) from data collection (obs) and focal id using simple ratio index (sri) 
make_network <- function(obs, focal.id) {
  N <- ncol(obs)
  network <- matrix(0,nrow=N, ncol=N)
  
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      x <- sum(obs[which(focal.id %in% c(i,j)),c(i,j)]>0)
      fa <- sum(focal.id == i)
      fb <- sum(focal.id == j)
      sri <- x/(fa+fb)
      if(!is.nan(sri)){
        network[i,j] <- sri
        network[j,i] <- sri
      }else{
        network[i,j] <- 0
        network[j,i] <- 0
        }
    }
  }
  return(network)
}

### function to generate pre-network permutations (swaps of individuals between focals)
rand_network2 <- function(obs.p, focal.id, n.perm,n_focals) {
  N <- ncol(obs.p)
  networks_rand <- array(0, c(n.perm,N,N))
  for (i in 1:n.perm) {
    # first randomly select two focal observations
    repeat {
      o <- 1:n_focals
      a <- sample(o,1)
      b <- sample(o[-a],1)
      
      # check if these are different individuals and they have associates
      if ((focal.id[a] != focal.id[b]) & (sum(obs.p[a,])>0) & (sum(obs.p[b,])>0)) {
        # next select two associates to swap
        d <- sample(which(obs.p[a,] > 0),1)
        e <- sample(which(obs.p[b,] > 0),1)
        
        # check they do not occur in the other focal
        if ((obs.p[a,e] == 0) & obs.p[b,d] == 0) {
          
          # now check we have 4 distinct individuals, otherwise repeat this process
          if (!(d %in% c(focal.id[a], focal.id[b], e)) & !(e %in% c(focal.id[a], focal.id[b], d))) {
            break;
          }
        }
      }
    }
    
    # swap individuals
    obs.p[a,d] <- 0
    obs.p[b,d] <- 1
    obs.p[b,e] <- 0
    obs.p[a,e] <- 1
    # caculate network
    networks_rand[i,,] <- make_network(obs.p,focal.id)
  }
  return(networks_rand)
}

### Function to allocate number of observations to groups
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

### MAIN SIMULATION FUNCTION ###
### Arguments ###
# GS --> numeric argument indicating group size
# ObsBia --> numeric argument indicating the degree of observation bias [0.5-1.0]
# FemPhenotypeBias --> boolean argument indicating whether a phenotype bias is present among females
# nfocals --> numeric argument indicating number of focal samples
# N.perm --> numeric argument indicating number of permutations 
Simulation<-function(GS,ObsBias,FemSexRatio,FemPhenotypeBias,nfocals,N.Perm)
{
  # Set parameters
  N <- GS
  n_focals <- nfocals
  # Generate nodes
  NumFem<-round(GS * FemSexRatio)
  NumMal<-GS - NumFem
  Sex<-c(rep("F",NumFem),rep("M",NumMal))
  Sex<-sample(Sex,GS,replace=F)
  ids <- data.frame(ID=1:(N),SEX=Sex)
  # Generate a distribution of group sizes
  group_size <- sample(c(1:(N/2)),n_focals,replace=TRUE)
  # Create blank observation matrix
  obs <- matrix(0,nrow=n_focals,ncol=N)
  ## set number of observations of an individual in a group per individual
  ids$OBS <- rand_vect(N,sum(group_size),pos.only=TRUE)
  ## Variables to Allocate individuals to groups, 
  GroupID<-c(1:n_focals)
  group_size.tmp <- group_size
  # IF Fem phenotype is stronger than males, start with males so that they end up in smaller groups
  if(FemPhenotypeBias == T)
  {
    which.males <- which(ids$SEX=="M")
    which.females <- which(ids$SEX=="F")
    for (i in which.males) 
    {
      g <- sample(GroupID[which(group_size.tmp>0)],ids$OBS[i])
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,i] <- 1
    }
    for (i in which.females) 
    {
      if ((sum(group_size.tmp>0) < ids$OBS[i])) 
      {
        Needed<-ids$OBS[i]-(sum(group_size.tmp>0))
        group.tmp<-group_size
        group.tmp[group_size.tmp>0]=0
        BiggestGroups<-sort(group.tmp,decreasing = T,index.return=T)$ix
        ExtraGroups<-BiggestGroups[1:Needed]
        g<-c(GroupID[which(group_size.tmp>0)],ExtraGroups)
      }else 
      {
        g<-sample(GroupID[which(group_size.tmp>0)],ids$OBS[i])
      }
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,i] <- 1
    }
  }else # IF Fem phenotype is equal to males, allocate indivdiuals to groups at random
  {
    Inds<-c(1:GS)
    for (. in 1:GS) 
    {
      id<-Inds[1]
      if(length(Inds)>1){id<-sample(Inds,1)}
      Inds<-Inds[-which(Inds==id)]
      if ((sum(group_size.tmp>0) < ids$OBS[id])) 
      {
        Needed<-ids$OBS[id]-(sum(group_size.tmp>0))
        Fullgroups<-which(group_size.tmp==0)
        ExtraGroups<-sample(Fullgroups,Needed,replace=F) 
        g<-c(GroupID[which(group_size.tmp>0)],ExtraGroups)
      }else 
      {
        g<-sample(GroupID[which(group_size.tmp>0)],ids$OBS[id])
      }
      group_size.tmp[g] <- group_size.tmp[g]-1
      obs[g,id] <- 1
    }
    
  }
  # Select a focal individual from each group
  focal.id <- apply(obs,1,function(x) { sample(which(x==1),1)})

  # Now remove cases where individuals occur in a group for which they are focal
  obs[cbind(1:n_focals,focal.id)] <- 0

  ## NOW DO NETWORK ANALYSIS ON THESE DATA
  # Calculate network
  Net.Ori <- make_network(obs,focal.id)

  # Remove some observations according to the degre of observation bias ObsBias
  # Generate probability of being observed (males=1,females=ObsBias)
  ids$OBS_PROB <- ObsBias
  ids$OBS_PROB[which(ids$SEX=="M")] <- 1

  # Remove observations from GBI
  obs.Bias <- obs
  for (i in 1:N) {
    obs.Bias[which(obs.Bias[,i] > 0),i] <- sample(c(0,1),sum(obs.Bias[,i]),replace=TRUE,prob=c(1-ids$OBS_PROB[i],ids$OBS_PROB[i]))
  }
  # Calculate new network
  Net.Bias <- make_network(obs.Bias,focal.id)
  # Calculate Degrees
  ids$DEGREE <- rowSums(Net.Ori)
  ids$DEGREE.Bias <- rowSums(Net.Bias)
  # Calculate effects
  coef.Ori <- coefficients(lm(DEGREE~SEX,data=ids))[2]
  coef.Bias <- coefficients(lm(DEGREE.Bias~SEX,data=ids))[2]
  ### Data permutations
  # Create random networks with pre-network permutations
  n.perm <- N.Perm
  networks_Perm <- rand_network2(obs.Bias, focal.id, n.perm,n_focals)
  # Calculate degree distribution for each network
  deg_Perm <- apply(networks_Perm,1,function(x) { rowSums(x)})
  # Get coefficients for each randomisation
  coefs_Perm <- apply(deg_Perm,2,function(x) { coefficients(lm(x~SEX,data=ids))[2] })
  ## Create random networks with Node permutations
  deg_Perm.Nodes <- matrix(0,nrow=N,ncol=N.Perm)
  coefs.Perm_Nodes <- rep(0,N.Perm)
  for (i in 1:N.Perm) {
    s <- sample(1:N)
    net <- Net.Bias[s,s]
    deg_Perm.Nodes <- rowSums(net)
    coefs.Perm_Nodes[i] <- coefficients(lm(deg_Perm.Nodes~SEX,data=ids))[2]
  }
  # Collect data in a list
  Result<-list()
  Result[[1]]<-sum(coef.Bias>coefs_Perm) / n.perm
  Result[[2]]<-sum(coef.Bias>coefs.Perm_Nodes) / n.perm
  Result[[3]]<-ids
  Result[[4]]<-obs
  Result[[5]]<-obs.Bias
  names(Result)<-c("P-value_PreNetwork_Perm","P-value_Node_Perm","IdsDataFrame","ObsDataFrame_Original","ObsDataFrame_Bias")
  Result
}
####################################################################################
### MAIN CODE TO RUN SIMULATIONS IN SCENARIOS WITH OBSERVATION BIAS [0.5-1.0]
### Sample parameter space with latin hypercube sampling
library(lhs)
NumCombinations<-500
VariablesToSample<-4
VarNames<-c("GroupSize",    ## Range 10-100
            "FEM.REMOVAL",  ## Range 0.5-1
            "FemSexRatio",  ## Range 0.2-0.8
            "FOCALS.NUM")   ## Range 100 - 2000
LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=4)
Mat[,1]<-round((10 + (LHS[,1]*(100-10))),0)
Mat[,2]<-round(0.5 + (LHS[,2]*(1-0.5)),2)
Mat[,3]<-round(0.2 + (LHS[,3]*(0.8-0.2)),2)
Mat[,4]<-round((100 + (LHS[,4]*(2000-100))),0)
FemPhenotypeBias<-c(TRUE,FALSE)
## number of repetitions per combination of parameters
nSim<-1 
#### set directories
setwd("C:/My_Directory/...")
dir<-"C:/My_Directory/.../SimulationFiles/"
### Run simulations
for (a in 1:length(FemPhenotypeBias))
{
  ### vectors to collect simulated data...
  GROUP.SIZE<-c()
  FEM.REMOVAL<-c()
  FEM.SEXRATIO<-c()
  FEM.BIAS<-c()
  FOCALS.NUM<-c()
  MEDIAN.DEGREE.MALES<-c()
  MEDIAN.DEGREE.MALES.BIAS<-c()
  MEDIAN.DEGREE.FEMALES<-c()
  MEDIAN.DEGREE.FEMALES.BIAS<-c()
  P.VALUE_PRE<-c()
  P.VALUE_NODES<-c()
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    for(c in 1:nSim)
    {
      filename <- paste(as.character(FemPhenotypeBias[a]),"_",
                        as.character(Mat[b,1]),"_",
                        as.character(Mat[b,2]),"_",
                        as.character(Mat[b,3]),"_",
                        as.character(Mat[b,4]),"_",
                        as.character(c),".RData",sep="")
      cat("processing Sim: ",b,"\\n")
      Result[[c]]<-Simulation(Mat[b,1],Mat[b,2],Mat[b,3],FemPhenotypeBias[a],Mat[b,4],N.Perm = 1000)
      GROUP.SIZE<-c(GROUP.SIZE,Mat[b,1])
      FEM.REMOVAL<-c(FEM.REMOVAL,Mat[b,2])
      FEM.SEXRATIO<-c(FEM.SEXRATIO,Mat[b,3])
      FOCALS.NUM<-c(FOCALS.NUM,Mat[b,4])
      FEM.BIAS<-c(FEM.BIAS,FemPhenotypeBias[a])
      DF<-Result[[c]][[3]]
      MEDIAN.DEGREE.MALES<-c(MEDIAN.DEGREE.MALES,median(DF$DEGREE[DF$SEX=="M"]))
      MEDIAN.DEGREE.MALES.BIAS<-c(MEDIAN.DEGREE.MALES.BIAS,median(DF$DEGREE.Bias[DF$SEX=="M"]))
      MEDIAN.DEGREE.FEMALES<-c(MEDIAN.DEGREE.FEMALES,median(DF$DEGREE[DF$SEX=="F"]))
      MEDIAN.DEGREE.FEMALES.BIAS<-c(MEDIAN.DEGREE.FEMALES.BIAS,median(DF$DEGREE.Bias[DF$SEX=="F"]))
      P.VALUE_PRE<-c(c(P.VALUE_PRE,Result[[c]][[1]]))
      P.VALUE_NODES<-c(P.VALUE_NODES,Result[[c]][[2]])
      ### save data from a simulation
      save(Result,file=paste(dir,filename,sep=""))
    }
  }
### save data from all simulations in scenarios with or without female phenotype bias
result<-data.frame(GROUP.SIZE,FEM.REMOVAL,FEM.SEXRATIO,FEM.BIAS,FOCALS.NUM,MEDIAN.DEGREE.MALES,MEDIAN.DEGREE.MALES.BIAS,MEDIAN.DEGREE.FEMALES,MEDIAN.DEGREE.FEMALES.BIAS,P.VALUE_PRE,P.VALUE_NODES)
save(result,file=paste("SimResults_ObsBias_FemPheno_",as.character(FemPhenotypeBias[a]),"_1000Perm.RData",sep = ""))
}

####################################################################################
### MAIN CODE TO RUN SIMULATIONS IN SCENARIOS WITH NO OBSERVATION BIAS
### Sample parameter space with latin hypercube sampling
library(lhs)
NumCombinations<-500
VariablesToSample<-3
VarNames<-c("GroupSize",    ## Range 10-100
            "FEM.REMOVAL",  ## Range 0.5-1
            "FemSexRatio",  ## Range 0.2-0.8
            "FOCALS.NUM")   ## Range 100 - 2000
LHS<-randomLHS(NumCombinations,VariablesToSample)
Mat<-matrix(NA,nrow=NumCombinations,ncol=4)
Mat[,1]<-round((10 + (LHS[,1]*(100-10))),0)
Mat[,2]<-rep(1,NumCombinations)              # no obs bias keep constant to 1
Mat[,3]<-round(0.2 + (LHS[,2]*(0.8-0.2)),2)
Mat[,4]<-round((100 + (LHS[,3]*(2000-100))),0)
FemPhenotypeBias<-c(TRUE,FALSE)
## number of repetitions per combination of parameters
nSim<-1 
#### set directories
setwd("C:/My_Directory/...")
dir<-"C:/My_Directory/.../SimulationFiles/"
### Run simulations
for (a in 1:length(FemPhenotypeBias))
{
  ### vectors to collect simulated data...
  GROUP.SIZE<-c()
  FEM.REMOVAL<-c()
  FEM.SEXRATIO<-c()
  FEM.BIAS<-c()
  FOCALS.NUM<-c()
  MEDIAN.DEGREE.MALES<-c()
  MEDIAN.DEGREE.MALES.BIAS<-c()
  MEDIAN.DEGREE.FEMALES<-c()
  MEDIAN.DEGREE.FEMALES.BIAS<-c()
  P.VALUE_PRE<-c()
  P.VALUE_NODES<-c()
  for(b in 1:nrow(Mat))
  {
    Result<-list()
    for(c in 1:nSim)
    {
      filename <- paste(as.character(FemPhenotypeBias[a]),"_",
                        as.character(Mat[b,1]),"_",
                        as.character(Mat[b,2]),"_",
                        as.character(Mat[b,3]),"_",
                        as.character(Mat[b,4]),"_",
                        as.character(c),".RData",sep="")
      cat("processing Sim: ",b,"\\n")
      Result[[c]]<-Simulation(Mat[b,1],Mat[b,2],Mat[b,3],FemPhenotypeBias[a],Mat[b,4],N.Perm = 1000)
      GROUP.SIZE<-c(GROUP.SIZE,Mat[b,1])
      FEM.REMOVAL<-c(FEM.REMOVAL,Mat[b,2])
      FEM.SEXRATIO<-c(FEM.SEXRATIO,Mat[b,3])
      FOCALS.NUM<-c(FOCALS.NUM,Mat[b,4])
      FEM.BIAS<-c(FEM.BIAS,FemPhenotypeBias[a])
      DF<-Result[[c]][[3]]
      MEDIAN.DEGREE.MALES<-c(MEDIAN.DEGREE.MALES,median(DF$DEGREE[DF$SEX=="M"]))
      MEDIAN.DEGREE.MALES.BIAS<-c(MEDIAN.DEGREE.MALES.BIAS,median(DF$DEGREE.Bias[DF$SEX=="M"]))
      MEDIAN.DEGREE.FEMALES<-c(MEDIAN.DEGREE.FEMALES,median(DF$DEGREE[DF$SEX=="F"]))
      MEDIAN.DEGREE.FEMALES.BIAS<-c(MEDIAN.DEGREE.FEMALES.BIAS,median(DF$DEGREE.Bias[DF$SEX=="F"]))
      P.VALUE_PRE<-c(c(P.VALUE_PRE,Result[[c]][[1]]))
      P.VALUE_NODES<-c(P.VALUE_NODES,Result[[c]][[2]])
      ### save data from a simulation
      save(Result,file=paste(dir,filename,sep=""))
    }
  }
  ### save data from all simulations in scenarios with or without female phenotype bias
  result<-data.frame(GROUP.SIZE,FEM.REMOVAL,FEM.SEXRATIO,FEM.BIAS,FOCALS.NUM,MEDIAN.DEGREE.MALES,MEDIAN.DEGREE.MALES.BIAS,MEDIAN.DEGREE.FEMALES,MEDIAN.DEGREE.FEMALES.BIAS,P.VALUE_PRE,P.VALUE_NODES)
  save(result,file=paste("SimResults_NoObsBias_FemPheno_",as.character(FemPhenotypeBias[a]),"_1000Perm.RData",sep = ""))
}




