########################################################################
## Van Irsel et al. 2021. State-dependent environmental sensitivity of reproductive success and survival in a shorebird. Ibis.
## 
## Authors: Jurrian van Irsel & Andrew M. Allen
## Date: 13-12-2021
##
## Title: R-Script for running multi-state live-dead recovery capture-mark-recapture models using the package RMark
##
######################################################################
##
## 1. MULTI-STATE LIVE-DEAD RECOVERY MODEL RUNNING WITHOUT TERRITORY QUALITY
##  1.1. Date preparation
##  1.2. Run RMark models
## 2. MULTI-STATE LIVE-DEAD RECOVERY MODEL RUNNING WITH TERRITORY QUALITY
##  2.1. Date preparation
##  2.2. Run RMark models
##
######################################################################
##
## This R-script is used to run the multi-state live-dead capture-mark-recapture models described in Van Irsel et al. 2021. State-dependent environmental sensitivity of reproductive success and survival in a shorebird. For computational matters, the script replicates only the top performing models (Table 2 in Van Irsel et al. 2021). However, the script can quickly be extended to include the complete model set. The script runs a multi-state live-dead recovery model with three states included, namely Non-breeders (N), Successful breeders (S) and Failed breeders (F). To avoid near boundery estimates we fixed the resighting probability of the Succesful and Failed breeder states (S and F) to 1. Hence, the resighting probability is only estimated for the Non-breeder state (N). 
#In the second part of the script, the multi-state live-dead capture-mark-recapture models are specified to run the territory quality analysis in the supplementary materials. In these models territory quality is used as an individual covariate to explain additional variation between individuals besides their previous breeding state. For readability we refer to territory quality as strategy from here on.
#Please be aware that we used MARK version 8.2 with RMark version 2.2.7. Using different versions of the package RMark and/or MARK may affect the results due to different model specifications in newer versions of RMark.
##
#####################################################################
## libraries
library(RMark) #RMark version 2.2.7 and MARK version 8.2

######################################################################
rm(list = ls()) # remove all files in the R workspace environment

##############################
# 1.1. Date preparation ---------------------------------------------
## Get data
ms = read.csv("Mark data.csv", check.names = TRUE, sep = ";")
colnames(ms)[2] <- "ch" #change colname to useful name in RMark
ms$ch <- as.character(ms$ch) #change the recapture data to character

summary(ms) #Check the data

#Get covariate data
covariates <- read.csv("covariates.csv", check.names = TRUE, sep = ";") #read covariate data of environmental covariates

##############################
#Prepare process data
ms.pr=process.data(data = ms, begin.time = c(2002), 
                   model = "MSLiveDead", 
                   strata.labels = c("S", "N", "F")) #create process data and specify strata, Successful breeder (S), Failed breeder (F) and Non-breeder (N)

ms.pr #check process data for any misspecifications

###############################
#Make design data and fix resighting probabilities
mscombi.ddl=make.design.data(ms.pr, parameters=list(Psi=list(pim.type="time")))

mscombi.ddl$p$fix=NA #create empty pim-chart
mscombi.ddl$p$fix[mscombi.ddl$p$stratum == "F" | mscombi.ddl$p$stratum == "S"] = 1 #Fix the resighting probability of the Failed (F) and Successful breeding state (S) to 1

summary(mscombi.ddl) #Check design data

##############################
#Add environmental covariates to the design data ------------------
mscovariate.ddl <- mscombi.ddl
mscovariate.ddl$S = merge_design.covariates(mscovariate.ddl$S, covariates) #Add covariates to the survival parameter
mscovariate.ddl$Psi = merge_design.covariates(mscovariate.ddl$Psi, covariates) #Add covariates to the transition parameter
mscovariate.ddl$S$wsinext #Check merging for survival
mscovariate.ddl$Psi$wsinext #Check merging for transition

#####################################################################
##
## 1.2. Run RMark models ------------------------------------
##
#####################################################################
#Specify MARK model structures
model.FSND_covariates = function()
{
  ##
  #Specify the model structure for the survival parameter
  S.stratum.time=list(formula = ~stratum + stratum:time)

  S.winter3 = list(formula=~ stratum + stratum:musselschier + stratum:cockle + stratum:wsinext)
  
  S.stratum.mussel_windchill = list(formula=~ stratum + stratum:musselschier + stratum:windchillnext)
  
  S.stratum.cockle_windchill = list(formula=~ stratum + stratum:cockle + stratum:windchillnext)
 
  ##
  #Specify model structure for the resighting parameter
  p.stratum.time=list(formula=~ -1 + stratum:time)
  
  ##
  #Specify model structure for the transition parameter
  Psi.stratum.time = list(formula=~ stratum:tostratum + stratum:tostratum:time, remove.intercept = TRUE)
  
  Psi.incubation1 = list(formula=~ stratum:tostratum + stratum:tostratum:surprecipincubationnext + stratum:tostratum:avgtempincubationnext + stratum:tostratum:tidalmaxnext, remove.intercept = TRUE)
  
  Psi.incubation2 = list(formula=~ stratum:tostratum + stratum:tostratum:ragwormnext+ stratum:tostratum:avgtempincubationnext + stratum:tostratum:tidalmaxnext, remove.intercept = TRUE)
  
  ##
  #model structure for the dead recovery parameter
  r.time=list(formula=~time)
  
  ##############################
  #Create model list as specify process data and design data
  ms.model.list=create.model.list("MSLiveDead")
  ms.results=mark.wrapper(ms.model.list,
                          data=ms.pr,
                          ddl=mscovariate.ddl, 
                          output=FALSE, adjust = TRUE, 
                          default.fixed=T)
  return(ms.results)
}

#Run model in MARK and store results
ms.results <- model.FSND_covariates()

#Extract top performing model in this case covariate model for Survival and Combination covariate model for the transition parameter (PSI)
top_perfrom_model = ms.results[["#model number#"]]

#Extract real values
RealRes <- top_perfrom_model$results$real
write.csv(RealRes, "Real_values.csv")

#Extract beta values
BetaRes <- top_perfrom_model$results$beta
write.csv(BetaRes, "Beta_values.csv")

#####################################################################
##
################ END OF MARK MODELS WITHOUT TERRITORY ###############
##
#####################################################################

#####################################################################
##
## 2.MULTI-STATE LIVE-DEAD RECOVERY MODEL RUNNING WITH TERRITORY QUALITY ------
##
#####################################################################
## This part of 
## libraries
library(RMark)

######################################################################
rm(list = ls()) # remove all files in the R workspace environment

##############################
# 2.1. Date preparation ---------------------------------------------
## Get data
ms = read.csv("Mark data_strategy.csv", check.names = TRUE, sep = ";")
colnames(ms)[3] <- "ch" #change colname to useful name in RMark
ms$ch <- as.character(ms$ch) #change the recapture data to character
ms$strategy <- as.factor(ms$strategy) #change the territory quality data to factor

summary(ms) #Check the data

##############################
#Prepare process data
ms.pr = process.data(data = ms, begin.time = c(2002), 
                   model = "MSLiveDead", 
                   strata.labels = c("S", "N", "F"),
                   groups = "strategy")

ms.pr #check process data for any misspecifications

###############################
#Make design data and fix resighting probabilities
mscombi.ddl = make.design.data(ms.pr, parameters=list(Psi=list(pim.type="time")))

mscombi.ddl$p$fix=NA #create empty pim-chart
mscombi.ddl$p$fix[mscombi.ddl$p$stratum == "F" | mscombi.ddl$p$stratum == "S"] = 1 #Fix the resighting probability of the states F and S to 1

summary(mscombi.ddl) #Check design data

##############################
#Add environmental covariates to the design data ------------------
mscovariate.ddl <- mscombi.ddl
mscovariate.ddl$S = merge_design.covariates(mscovariate.ddl$S, covariates) #Add covariates to the survival parameter
mscovariate.ddl$Psi = merge_design.covariates(mscovariate.ddl$Psi, covariates) #Add covariates to the transition parameter
mscovariate.ddl$S$wsinext #Check merging for survival
mscovariate.ddl$Psi$wsinext #Check merging for transition

#####################################################################
##
## 2.2. Run RMark models ------------------------------------
##
#####################################################################
#Specify MARK model structures
model.FSND_strategy_covariates = function()
{
  ##
  #Specify the model structure for the survival parameter
  S.stratum.time=list(formula = ~stratum + stratum:time)
  
  S.stratum.strategy.time=list(formula = ~stratum:strategy + stratum:strategy:time)
  
  S.strategy.winter3=list(formula=~ stratum:strategy + stratum:strategy:musselschier + stratum:strategy:cockle + stratum:strategy:wsinext)
  
  S.strategy.winter2=list(formula=~ stratum:strategy + stratum:strategy:surprecipwinternext + stratum:strategy:cockle + stratum:strategy:windchillnext)
  
  S.stratum.strategy.cockle_windchill = list(formula = ~stratum:strategy + stratum:strategy:cockle + stratum:strategy:windchillnext)
  
  S.stratum.strategy.mussel_windchill = list(formula = ~stratum:strategy + stratum:strategy:musselschier + stratum:strategy:windchillnext)
  
  ##
  #Specify model structure for the resighting parameter
  p.stratum.time=list(formula=~ -1 + stratum:time)
  
  ##
  #Specify model structure for the transition parameter
  Psi.stratum.time = list(formula=~ stratum:tostratum:strategy + stratum:tostratum:strategy:time, remove.intercept = TRUE)
  
  Psi.strategy.incubation1 = list(formula=~ stratum:tostratum:strategy + stratum:tostratum:strategy:surprecipincubationnext + stratum:tostratum:strategy:avgtempincubationnext + stratum:tostratum:strategy:tidalmaxnext, remove.intercept = TRUE)
  
  Psi.strategy.incubation2 = list(formula=~ stratum:tostratum:strategy + stratum:tostratum:strategy:ragwormnext+ stratum:tostratum:strategy:avgtempincubationnext + stratum:tostratum:strategy:tidalmaxnext, remove.intercept = TRUE)
  
  ##
  #model structure for the dead recovery parameter
  r.time=list(formula=~time)
  
  ##############################
  #Create model list as specify process data and design data
  ms.model.list = create.model.list("MSLiveDead")
  ms.results = mark.wrapper(ms.model.list,
                          data=ms.pr,
                          ddl=mscovariate.ddl, 
                          output=FALSE, 
                          adjust = TRUE, 
                          default.fixed=T)
  return(ms.results)
}

#Run model in MARK and store results
ms.results_strategy <- model.FSND_strategy_covariates()

#Extract top performing model in this case covariate model for Survival and Combination covariate model for the transition parameter (PSI)
top_perfrom_strategy_model = ms.results_strategy[["#model number#"]]

#Extract real values
RealRes_strategy <- top_perfrom_strategy_model$results$real
write.csv(RealRes_strategy, "Real_values_strategy.csv")

#Extract beta values
BetaRes_strategy <- top_perfrom_strategy_model$results$beta
write.csv(BetaRes_strategy, "Beta_values_strategy.csv")

#####################################################################
##
################ END OF MARK MODELS WITH TERRITORY QUALITY ##########
##
#####################################################################
