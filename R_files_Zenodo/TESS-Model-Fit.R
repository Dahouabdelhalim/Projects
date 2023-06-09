library(TESS)
library(parallel)
library(digest)

#Parameters
setwd("PATH/TO/FOLDER")
treefilename <- "trees.trees"
n.trees <- 30 #how many trees to use
exclusion <- as.numeric(c("70", "84", "71", "82", "43", "12", "27", "46", "63", "83", "2", "11", "60", "1", "89", "37", "75", "72", "81", "65", "29", "21", "96", "88", "42", "22", "92", "67", "100", "24", "80", "79", "9", "30", "14", "73", "87", "76", "99", "45", "18", "32", "49", "68", "91", "6", "16", "66", "3", "35", "41", "34", "86", "17", "78", "38", "28", "98", "56", "50", "31", "53", "85", "95", "36", "77", "97", "15", "47", "26")) #List of trees not to sample (ie if already analyzed on previous run)
seed <- 12345
set.seed(seed) #Comment to randomize tree selection
samplingProb <-  1 #For the likelihood functions
samplingStrat <- "diversified" #You can also a random sampling strategy
pRateChangeTime <- 10 #Percent of tree age at which rate change occurs
marginalIterations <- 1000 #number of iterations for calculating marginal Likelihood
marginalBurnin <- marginalIterations/10 #number of burn ins for calculating marginal likelihood
no_cores <- detectCores() - 1 #Uses one less than all the cores. Modify as necessary

#Attempt to close previous clusters in case previously abandoned
tryCatch({
  cl <- makeCluster(no_cores, type="FORK")
  parSapply(cl,1,length)
  stopCluster(cl)
},error = function(e){cat("Attempting to close previous clusters...")})

tryCatch(file.remove("tess_Output/treesDone.csv"),warning= function(w)NULL)


#Import data
mytree <- read.nexus(treefilename)


#Subset the chosen number of trees from full list
trees.2.use <- sample(setdiff(seq(1, length(mytree)),exclusion), min(n.trees,length(mytree)-length(exclusion)), replace = F); ## use trees after 25% burnin
times <- sapply(1:length(trees.2.use), function(x){as.numeric(branching.times(mytree[[trees.2.use[x]]]))})
colnames(times) <- trees.2.use
rm(mytree)#clear up some memory, all analyses in this script are based solely off the vector times

#Output list of trees used this run
if(!dir.exists("tess_Output"))
  dir.create("tess_Output")
write.table(matrix(trees.2.use,ncol=(length(trees.2.use))), file = "tess_Output/treesUsed.csv", row.names=FALSE, col.names=FALSE, sep=",")

#-----------------------------------------------------------#

#Set up functions to get Marginal Likelihoods (pass in index for times / trees.2.use)
getMarginalConstBD <-function(x){
  err <- NULL;
  withCallingHandlers(tryCatch({
    #Set likelihood function for Constant BD
    likelihoodConstBD <- function(params) {
      speciation <- params[1] + params[2]
      extinction <- params[2]
      lnl <- tess.likelihood(times[,x],
                             lambda = speciation,
                             mu = extinction,
                             samplingProbability = samplingProb,
                             samplingStrategy = samplingStrat,
                             log = TRUE)
      return (lnl)
    }
    
    #Calculate and retun Marginal Likelihood for Constant BD
    return(
      tess.steppingStoneSampling(
        likelihoodFunction = likelihoodConstBD,
        priors = priorsConstBD,
        parameters = runif(2,0,1),
        logTransforms = c(TRUE,TRUE),
        iterations = marginalIterations,
        burnin = marginalBurnin,
        K = 50)
    )
  },
  error=function(e) {
    err <<- conditionMessage(e);
    if(!dir.exists("debug"))
      dir.create("debug")
    cat(dput(err), file = paste0("debug/debug_ConstBD_tree", trees.2.use[x], ".txt"))
    return(NA)
  }
  ))
}

getMarginalDecrBD <- function(x){
  err <- NULL;
  withCallingHandlers(tryCatch({
    #Set likelihood function for Decreasing BD
  likelihoodDecrBD <- function(params) {
    speciation <- function(t) params[1] + params[2] * exp(-params[3]*t)
    extinction <- function(t) params[1]
    lnl <- tess.likelihood(times[,x],
                           lambda = speciation,
                           mu = extinction,
                           samplingProbability = samplingProb,
                           samplingStrategy = samplingStrat,
                           log = TRUE)
    return (lnl)
  }
  
  #Calculate and retun Marginal Likelihood for Deccreasing BD
  return(
    tess.steppingStoneSampling(
      likelihoodFunction = likelihoodDecrBD,
      priors = priorsDecrBD,
      parameters = runif(3,0,1),
      logTransforms = c(TRUE,TRUE,TRUE),
      iterations = marginalIterations,
      burnin = marginalBurnin,
      K = 50)
  )
  },
  error=function(e) {
    err <<- conditionMessage(e);
    if(!dir.exists("debug"))
      dir.create("debug")
    cat(dput(err), file = paste0("debug/debug_DecrBD_tree", trees.2.use[x], ".txt"))
    return(NA)
  }
  ))
}

getMarginalEpisodicBD <- function(x){
  err <- NULL;
  withCallingHandlers(tryCatch({
    #Set likelihood function for Episodic BD
  likelihoodEpisodicBD <- function(params) {
    speciation <- c(params[1]+params[2],params[3]+params[4])
    extinction <- c(params[2],params[4])
    lnl <- tess.likelihood.rateshift(times[,x],
                                     lambda = speciation,
                                     mu = extinction,
                                     rateChangeTimesLambda = max(times[,x])/100*pRateChangeTime,
                                     rateChangeTimesMu = max(times[,x])/100*pRateChangeTime,
                                     samplingProbability = samplingProb,
                                     samplingStrategy = samplingStrat,
                                     log = TRUE)
    return (lnl)
  }
  
  #Calculate and retun Marginal Likelihood for Episodic BD
  return(
    tess.steppingStoneSampling(
      likelihoodFunction = likelihoodEpisodicBD,
      priors = priorsEpisodicBD,
      parameters = runif(4,0,1),
      logTransforms = c(TRUE,TRUE,TRUE,TRUE),
      iterations = marginalIterations,
      burnin = marginalBurnin,
      K = 50)
  )
  },
  error=function(e) {
    err <<- conditionMessage(e);
    if(!dir.exists("debug"))
      dir.create("debug")
    cat(dput(err), file = paste0("debug/debug_EpisodicBD_tree", trees.2.use[x], ".txt"))
    return(NA)
  }
  ))
}

#Set up a Cacheing wrapper function, depends on setting the global variabl "which.function" as the function you want to run
cacheWrapTess <- function(i){
  err <- NULL;
  withCallingHandlers(tryCatch({
    #dg <- digest(list(which.function, i))
    cache_filename <- sprintf("cache/CacheTess_%s_%s.Rdata",which.debug,as.character(i))
    if (file.exists(cache_filename))
      load(cache_filename)
    else{
      i <- which.function(i);
      if(!dir.exists("cache"))
        dir.create("cache")
      if(!is.na(i)) save(i, file = cache_filename)
      
      #Create and/or update progress file
      if(!dir.exists("tess_Output"))
        dir.create("tess_Output")
      treesDone <- as.matrix(t(c(0,0,0)))
      colnames(treesDone) <- c("Const","Decr","Episodic")
      tryCatch(treesDone <- read.csv("tess_Output/treesDone.csv", header=T, row.names = 1),
               warning = function(w) treesDone <- as.matrix(t(c(0,0,0))))
      if(which.debug == "ConstBD")
        treesDone[1] <- treesDone[1] +1
      if(which.debug == "DecrBD")
        treesDone[2] <- treesDone[2] +1
      if (which.debug == "EpisodicBD")
        treesDone[3] <- treesDone[3] +1
      write.csv(treesDone, file = "tess_Output/treesDone.csv")
    }
    return(i)
  },
  error=function(e) {
    err <<- conditionMessage(e);
    if(!dir.exists("debug"))
      dir.create("debug")
    cat(dput(err), file = paste0("debug/debug_cache",which.debug,"_tree", trees.2.use[i], ".txt"))
    return(NA)
  }
  ))
}

##Set priors for likelihood functions

#Birth-death processes with constant rates, 
prior_delta <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsConstBD <- c("diversification"=prior_delta,
                   "turnover"=prior_tau)

#Birth-death processes with continuously varying rates 
prior_delta <- function(x) { dexp(x,rate=0.1,log=TRUE) }
prior_lambda <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_alpha <- function(x) { dexp(x,rate=0.1,log=TRUE) }
priorsDecrBD <- c("turnover"=prior_delta,
                  "initial speciation"=prior_lambda,
                  "speciation decay"=prior_alpha)

#Birth-death processes with episodically varying rates
prior_delta_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_before <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_delta_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
prior_tau_after <- function(x) { dexp(x,rate=10.0,log=TRUE) }
priorsEpisodicBD <- c("diversification before"=prior_delta_before,
                      "turnover before"=prior_tau_before,
                      "diversification after"=prior_delta_after,
                      "turnover after"=prior_tau_after)

rm(list=ls(pattern="prior_")) #Clear out the temporary function used to build the three priors lists

##Calculate Marginal Likelihoods for all the trees
#ConstantBD
which.function <- getMarginalConstBD
which.debug <- "ConstBD" #Used to help name debug files in case of error
cl <- makeCluster(no_cores, type="FORK")
system.time(MLConstBD <- parSapply(cl, 1:n.trees,cacheWrapTess))
stopCluster(cl)

#DecrBD
which.function <- getMarginalDecrBD
which.debug <- "DecrBD" #Used to help name debug files in case of error
cl <- makeCluster(no_cores, type="FORK")
system.time(MLDecrBD <- parSapply(cl, 1:n.trees,cacheWrapTess))
stopCluster(cl)

#EpisodicBD
which.function <- getMarginalEpisodicBD
which.debug <- "EpisodicBD" #Used to help name debug files in case of error
cl <- makeCluster(no_cores, type="FORK")
system.time(MLEpisodicBD <- parSapply(cl, 1:length(trees.2.use),cacheWrapTess))
stopCluster(cl)

#For testing
# MLConstBD <- sample(1:100, n.trees)
# MLDecrBD <- sample(1:100, n.trees)
# MLEpisodicBD <- sample(1:100, n.trees)

##Calculate Bayes Factors
#Organize marginal likelihoods from each model into matrix marginal.data
marginal.data <- sapply(1:n.trees, function(i) c("ConstBD"=MLConstBD[i],"DecrBD"=MLDecrBD[i],"EpisodicBD"=MLEpisodicBD[i]))
colnames(marginal.data) <- trees.2.use  

#Create a two column matrix where each row has each combination of models to compare (used to index when calculating BF)
modelGrid <- expand.grid(M0=names(marginal.data[,1]),M1=names(marginal.data[,1]))

#Calculate BF, rename columns, rows
BF.data <- sapply(1:n.trees, function(i) 2 * (marginal.data[modelGrid$M0,i] - marginal.data[modelGrid$M1,i]))
colnames(BF.data) <- trees.2.use
rownames(BF.data) <- NULL #Prevents conflict when cbind is called below

#add the modelGrid to first the beginning as a quasi-rowname
BF.data <- as.data.frame(BF.data)
BF.data <- cbind(modelGrid, BF.data)

#Create a single row dataframe of the most supported model for each tree
supportedModel <- as.data.frame(t(sapply(1:(n.trees+2), function(x){
  tryCatch({
    if(any(which(BF.data[x] == max(BF.data[x])))) 
      as.character(BF.data$M0[which(BF.data[x] == max(BF.data[x]))])
    else NA
  }
  ,error = function(e)NA)#used for first two columns where max() throws an error
})))
names(supportedModel) <- names(BF.data) #rename col to allow the rbind below
#BF.data <- rbind(BF.data,supportedModel) #Add supportedModel to the last row of BF.data

#Output Data
if(!dir.exists("tess_Output"))
  dir.create("tess_Output")

write.csv(marginal.data, file = "tess_Output/tess-MarginalLikelihood.csv")
write.csv(BF.data, file = "tess_Output/tess-BF.csv", row.names = F)
write.csv(supportedModel, file = "tess_Output/tess-SupportedModel.csv", row.names = F)

#Option to clear cache after all things are done running
tryCatch({
  if (as.logical(readline(prompt = "\\n\\nDon't forget to rename files in tess_Output folder to keep your results! \\nClear cache and debug? (T/F): ")))
    file.remove(dir("cache",full.names = T), dir("debug",full.names = T))
}
,error = function(e){print("Unrecognized input! Cache *not* cleared!")}
)

