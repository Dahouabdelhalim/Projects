#-------------------------------------------------------------------------------------------#
#                             OCCUPANCY MODELS FOR PHILORIA                                 #
#-------------------------------------------------------------------------------------------#

##### REQUIRED LIBRARIES #####
library(tuple)
library(lubridate)
library(unmarked)
library(AICcmodavg)
library(dplyr)
library(readxl)

#### IMPORT DATA AND PREPARE FOR ANALYSIS ####
setwd('') # Add own

### P. kundagungan
site.covs.kund <- as.data.frame(read_excel('P_kundagungan_2020-2021_data.xlsx', sheet=1, na="NA", col_types=c(rep("text", 2), rep("numeric", 10))))[, 3:12]
obs.covs.kund <- as.data.frame(read_excel('P_kundagungan_2020-2021_data.xlsx', sheet=2, na="NA", col_types=c("text", rep("numeric", 5))))[, 2:6]
dets.kund <- as.data.frame(read_excel('P_kundagungan_2020-2021_data.xlsx', sheet=3, na="NA", col_types=c("text", rep("numeric", 5))))[, 2:6]

## Data adjustments
dets.kund <- ifelse(is.na(dets.kund), NA, 
                    ifelse(dets.kund > 0,1,0)) # convert counts to detection and non-detection data 
obs.covs.kund <- (as.matrix(obs.covs.kund) - mean(as.matrix(obs.covs.kund), na.rm=T))/(2*sd(as.matrix(obs.covs.kund), na.rm=T)) # centre observation covariates
obs.covs.kund <- ifelse(is.na(obs.covs.kund), 0, obs.covs.kund) # replace missing covariate values with mean (0)
site.burnt.kund <- site.covs.kund$BurntOrUnburnt
for(i in 1:ncol(site.covs.kund)){
  site.covs.kund[, i] <- (site.covs.kund[, i] - mean(site.covs.kund[, i], na.rm=T)) / (2*sd(site.covs.kund[, i], na.rm=T)) # centre site covariates
}
site.covs.kund$BurntOrUnburnt <- site.burnt.kund # put binary burnt or not burn variable back in
head(site.covs.kund)

### P. richmondensis
site.covs.rich <- as.data.frame(read_excel('P_richmondensis_2020-2021_data.xlsx', sheet=1, na="NA", col_types=c(rep("text", 2), rep("numeric", 10))))[, 3:12]
obs.covs.rich <- as.data.frame(read_excel('P_richmondensis_2020-2021_data.xlsx', sheet=2, na="NA", col_types=c("text", rep("numeric", 3))))[, 2:4]
dets.rich <- as.data.frame(read_excel('P_richmondensis_2020-2021_data.xlsx', sheet=3, na="NA", col_types=c("text", rep("numeric", 3))))[, 2:4]

## Data adjustments
dets.rich <- ifelse(is.na(dets.rich), NA, 
                    ifelse(dets.rich > 0,1,0)) # convert counts to detection and non-detection data 
obs.covs.rich <- (as.matrix(obs.covs.rich) - mean(as.matrix(obs.covs.rich), na.rm=T))/(2*sd(as.matrix(obs.covs.rich), na.rm=T)) # centre observation covariates
obs.covs.rich <- ifelse(is.na(obs.covs.rich), 0, obs.covs.rich) # replace missing covariate values with mean (0)
site.burnt.rich <- site.covs.rich$BurntOrUnburnt
for(i in 1:ncol(site.covs.rich)){
  site.covs.rich[, i] <- (site.covs.rich[, i] - mean(site.covs.rich[, i], na.rm=T)) / (2*sd(site.covs.rich[, i], na.rm=T)) # centre site covariates
}
site.covs.rich$BurntOrUnburnt <- site.burnt.rich # put binary burnt or not burnt variable back in
head(site.covs.rich)

#### COMPLETE MODEL SELECTION IN UNMARKED #####

### Convert datasets to unmarked dataframes
kund.dat <- unmarkedFrameOccu(y=dets.kund, siteCovs=site.covs.kund, obsCovs=list(AirTemp=obs.covs.kund))   # organize data
summary(kund.dat) # take a look
rich.dat <- unmarkedFrameOccu(y=dets.rich, siteCovs=site.covs.rich, obsCovs=list(AirTemp=obs.covs.rich))   # organize data
summary(rich.dat) # take a look

### FIT MODELS ###

## First define model structures to include in the analysis
names(site.covs.kund)
mod.vars <- c("Elevation", "BurntOrUnburnt", "PercentBedBurnt", "AverageBankBurnSeverity", "GeebamFireSeverity",     
"SpinachCover", "AverageMuddiness", "MinimumMuddiness", "DroughtSeverityNd2", "DroughtSeverityVARI")

## Define all the possible combinations of variables
first.coms <- combn(mod.vars, 1)[1, ]
second.coms <- combn(mod.vars, 2)
second.coms.vect <- vector(length=ncol(second.coms))
for(i in 1:length(second.coms.vect)){
  second.coms.vect[i] <- paste(second.coms[, i], collapse="+")
}
third.coms <- combn(mod.vars, 3)
third.coms.vect <- vector(length=ncol(third.coms))
for(i in 1:length(third.coms.vect)){
  third.coms.vect[i] <- paste(third.coms[, i], collapse="+")
}
fourth.coms <- combn(mod.vars, 4)
fourth.coms.vect <- vector(length=ncol(fourth.coms))
for(i in 1:length(fourth.coms.vect)){
  fourth.coms.vect[i] <- paste(fourth.coms[, i], collapse="+")
}
fifth.coms <- combn(mod.vars, 5)
fifth.coms.vect <- vector(length=ncol(fifth.coms))
for(i in 1:length(fifth.coms.vect)){
  fifth.coms.vect[i] <- paste(fifth.coms[, i], collapse="+")
}
mod.struc <- c(first.coms, second.coms.vect, third.coms.vect, fourth.coms.vect, fifth.coms.vect)

## Take out models that include correlated variables
x <- unique(c(grep("BurntOrUnburnt.*PercentBedBurnt", mod.struc),
              grep("BurntOrUnburnt.*AverageBankBurnSeverity", mod.struc),
              grep("BurntOrUnburnt.*GeebamFireSeverity", mod.struc),
              grep("PercentBedBurnt.*AverageBankBurnSeverity", mod.struc),
              grep("PercentBedBurnt.*GeebamFireSeverity", mod.struc),
              grep("AverageBankBurnSeverity.*GeebamFireSeverity", mod.struc),
              grep("AverageMuddiness.*MinimumMuddiness", mod.struc),
              grep("AverageMuddiness.*DroughtSeverityNd2", mod.struc),
              grep("AverageMuddiness.*DroughtSeverityVARI", mod.struc),
              grep("MinimumMuddiness.*DroughtSeverityNd2", mod.struc),
              grep("MinimumMuddiness.*DroughtSeverityVARI", mod.struc),
              grep("DroughtSeverityNd2.*DroughtSeverityVARI", mod.struc)))
length(x)
mod.struc <- mod.struc[-x]
length(mod.struc)

## Take out models that don't include an effect of elevation (elevation in all models a priori)
x <- grep("Elevation", mod.struc)
mod.struc <- mod.struc[x]
length(mod.struc)
mod.struc

## Check for duplicate models
dupes <- duplicated(mod.struc) 
sum(dupes) # none

## Model function (to run all of the models in 'mod.struc' sequentially)
mod.func <- function(mod.num, dataset){
  # Specify the model format
  mod.str <- paste(paste('~', 'AirTemp'), paste('~', mod.struc[mod.num])) 
  # Fit the model
  mod <- occu(as.formula(mod.str), data=dataset)
  # Gather model selection statistics
  aic <- aictab(list(mod), modnames=FALSE)
  # Gather any error messages
  err <- ifelse(mod@opt$convergence==0, 'None', "Failed to converge")
  # Gather up information for the model selection results
  return(c(mod.num, mod.struc[mod.num], aic$K, aic$LL, aic$AICc, err))
}

### FIT MODELS FOR P. KUNDAGUNGAN

## Fit each of the models in the set and return model selection statistics
kund.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(kund.mod.mat) <- c("Model", "OccupancyVariables", "Nparameters", "LogLik", "AIC", "Messages")
for(i in 1:nrow(kund.mod.mat)){
  kund.mod.mat[i, ] <- mod.func(mod.num=i, dataset=kund.dat)
}

# Add calculation of model weight for each model
mod.mat.dframe <- as.data.frame(kund.mod.mat)
mod.mat.dframe$deltaAIC <- as.numeric(as.character(mod.mat.dframe$AIC)) - min(as.numeric(as.character(mod.mat.dframe$AIC)))
mod.mat.dframe$weight <- exp(-0.5 * mod.mat.dframe$deltaAIC)
mod.mat.dframe$sum.weight <- sum(mod.mat.dframe$weight)
mod.mat.dframe$ModelWeight <- mod.mat.dframe$weight / mod.mat.dframe$sum.weight

# Some rejigging before export to get into shape
mod.mat.dframe <- mod.mat.dframe[, -c(8:9)]
mod.mat.dframe <- mod.mat.dframe[order(as.numeric(as.character(mod.mat.dframe$AIC))), ] # reorder by AIC
mod.mat.dframe <- mod.mat.dframe[, c(1:5,7:8,6)] # reorder the columns so that 'Messages' is last

# Export model selection statistics (and any error messages)
setwd('C:\\\\Users\\\\User\\\\Dropbox\\\\Research\\\\Papers\\\\Active\\\\Philoria_post_fire\\\\EcologyEvolution_submission\\\\For_submission\\\\Zenodo_archive')
write.csv(mod.mat.dframe, file='Model_selection_statistics_P_kundagungan_occupancy.csv')

# Refit top model to look at regression coefficients
top.mod.Pk <- occu(as.formula(paste('~ ', 'AirTemp', '~ ', mod.mat.dframe$OccupancyVariables[1])), data=kund.dat)
summary(top.mod.Pk)

# Confidence intervals for parameters of top model
confint(top.mod.Pk, type='state')
confint(top.mod.Pk, type='det')

### FIT MODELS FOR P. RICHMONDENSIS

## Fit each of the models in the set and return model selection statistics
rich.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(rich.mod.mat) <- c("Model", "OccupancyVariables", "Nparameters", "LogLik", "AIC", "Messages")
for(i in 1:nrow(rich.mod.mat)){
  rich.mod.mat[i, ] <- mod.func(mod.num=i, dataset=rich.dat)
}

# Add calculation of model weight for each model
mod.mat.dframe <- as.data.frame(rich.mod.mat)
mod.mat.dframe$deltaAIC <- as.numeric(as.character(mod.mat.dframe$AIC)) - min(as.numeric(as.character(mod.mat.dframe$AIC)))
mod.mat.dframe$weight <- exp(-0.5 * mod.mat.dframe$deltaAIC)
mod.mat.dframe$sum.weight <- sum(mod.mat.dframe$weight)
mod.mat.dframe$ModelWeight <- mod.mat.dframe$weight / mod.mat.dframe$sum.weight

# Some rejigging before export to get into shape
mod.mat.dframe <- mod.mat.dframe[, -c(8:9)]
mod.mat.dframe <- mod.mat.dframe[order(as.numeric(as.character(mod.mat.dframe$AIC))), ] # reorder by AIC
mod.mat.dframe <- mod.mat.dframe[, c(1:5,7:8,6)] # reorder the columns so that 'Messages' is last

# Export model selection statistics (and any error messages)
setwd('C:\\\\Users\\\\User\\\\Dropbox\\\\Research\\\\Papers\\\\Active\\\\Philoria_post_fire\\\\EcologyEvolution_submission\\\\For_submission\\\\Zenodo_archive')
write.csv(mod.mat.dframe, file='Model_selection_statistics_P_richmondensis_occupancy.csv')

# Refit top model to look at regression coefficients
top.mod.Pr <- occu(as.formula(paste('~ ', 'AirTemp', '~ ', mod.mat.dframe$OccupancyVariables[1])), data=rich.dat)
summary(top.mod.Pr)

# Confidence intervals for parameters of top model
confint(top.mod.Pr, type='state')
confint(top.mod.Pr, type='det')

### COVARIATE RELATIONSHIPS FROM THE TOP MODELS ####

## P. kundagungan

# Data to predict to
elevation.unburnt <- data.frame(BurntOrUnburnt=rep(0,21), AverageMuddiness=rep(0, 21), Elevation=(seq(600, 1000, 20)-876.0775)/(2*106.3115)) # range of elevation for unburnt sites at mean muddiness (with standardisation of elevation, as per analysis) 
elevation.burnt <- data.frame(BurntOrUnburnt=rep(1,21), AverageMuddiness=rep(0, 21), Elevation=(seq(600, 1000, 20)-876.0775)/(2*106.3115)) # range of elevation for burnt sites at mean muddiness (with standardisation of elevation, as per analysis) 
muddy.unburnt <- data.frame(BurntOrUnburnt=rep(0,26), AverageMuddiness=(seq(0, 100, 4)-62.02431)/(2*30.37526), Elevation=rep(0, 26)) # range of average muddiness for unburnt sites at mean elevation (with standardisation of muddiness, as per analysis) 
muddy.burnt <- data.frame(BurntOrUnburnt=rep(1,26), AverageMuddiness=(seq(0, 100, 4)-62.02431)/(2*30.37526), Elevation=rep(0, 26)) # range of average muddiness for burnt sites at mean elevation (with standardisation of muddiness, as per analysis) 

# Predictions and export
setwd('C:\\\\Users\\\\User\\\\Dropbox\\\\Research\\\\Papers\\\\Active\\\\Philoria_post_fire\\\\EcologyEvolution_submission\\\\For_submission\\\\Zenodo_archive')
elev.unburnt.pred <- predict(top.mod.Pk, type='state', newdata=elevation.unburnt)
write.csv(elev.unburnt.pred, "P_kundagungan_occ_vs_elevation_unburnt.csv")
elev.burnt.pred <- predict(top.mod.Pk, type='state', newdata=elevation.burnt)
write.csv(elev.burnt.pred, "P_kundagungan_occ_vs_elevation_burnt.csv")
muddy.unburnt.pred <- predict(top.mod.Pk, type='state', newdata=muddy.unburnt)
write.csv(muddy.unburnt.pred, "P_kundagungan_occ_vs_muddiness_unburnt.csv")
muddy.burnt.pred <- predict(top.mod.Pk, type='state', newdata=muddy.burnt)
write.csv(muddy.burnt.pred, "P_kundagungan_occ_vs_muddiness_burnt.csv")

## P. richmondensis

# Data to predict to
elev.pred.dat <- data.frame(AverageMuddiness=rep(0, 21), Elevation=(seq(400, 800, 20)-571.8822)/(2*87.83562)) # range of elevation at mean muddiness (with standardisation of elevation, as per analysis) 
muddy.pred.dat <- data.frame(AverageMuddiness=(seq(0, 100, 4)-68.53667)/(2*29.26409), Elevation=rep(0, 26)) # range of average muddiness at mean elevation (with standardisation of muddiness, as per analysis) 

# Predictions and export
setwd('C:\\\\Users\\\\User\\\\Dropbox\\\\Research\\\\Papers\\\\Active\\\\Philoria_post_fire\\\\EcologyEvolution_submission\\\\For_submission\\\\Zenodo_archive')
elev.pred <- predict(top.mod.Pr, type='state', newdata=elev.pred.dat)
write.csv(elev.pred, "P_richmondensis_occ_vs_elevation.csv")
muddy.pred <- predict(top.mod.Pr, type='state', newdata=muddy.pred.dat)
write.csv(muddy.unburnt.pred, "P_richmondensis_occ_vs_muddiness.csv")