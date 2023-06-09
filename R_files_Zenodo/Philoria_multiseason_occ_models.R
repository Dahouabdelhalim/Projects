#-------------------------------------------------------------------------------------------#
#                     MULTISEASON OCCUPANCY MODELS FOR PHILORIA                             #
#-------------------------------------------------------------------------------------------#

##### REQUIRED LIBRARIES #####
library(tuple)
library(dplyr)
library(unmarked)
library(AICcmodavg)
library(readxl)

### IMPORT DATA AND PREPARE FOR ANALYSIS ####
setwd('') # Add own

### P. kundagungan
site.covs.kund <- as.data.frame(read_excel('P_kundagungan_multiseason_data.xlsx', sheet=1, na="NA", col_types=c(rep("text", 1), rep("numeric", 10))))[, 2:11]
obs.covs.kund <- as.data.frame(read_excel('P_kundagungan_multiseason_data.xlsx', sheet=2, na="NA", col_types=c("text", rep("numeric", 10))))[, 2:11]
dets.kund <- as.data.frame(read_excel('P_kundagungan_multiseason_data.xlsx', sheet=3, na="NA", col_types=c("text", rep("numeric", 10))))[, 2:11]

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
site.covs.rich <- as.data.frame(read_excel('P_richmondensis_multiseason_data.xlsx', sheet=1, na="NA", col_types=c(rep("text", 1), rep("numeric", 10))))[, 2:11]
obs.covs.rich <- as.data.frame(read_excel('P_richmondensis_multiseason_data.xlsx', sheet=2, na="NA", col_types=c("text", rep("numeric", 9))))[, 2:10]
dets.rich <- as.data.frame(read_excel('P_richmondensis_multiseason_data.xlsx', sheet=3, na="NA", col_types=c("text", rep("numeric", 9))))[, 2:10]

## Data adjustments
dets.rich <- ifelse(is.na(dets.rich), NA, 
                    ifelse(dets.rich > 0,1,0)) # convert counts to detection and non-detection data 
obs.covs.rich <- (as.matrix(obs.covs.rich) - mean(as.matrix(obs.covs.rich), na.rm=T))/(2*sd(as.matrix(obs.covs.rich), na.rm=T)) # centre observation covariates
obs.covs.rich <- ifelse(is.na(obs.covs.rich), 0, obs.covs.rich) # replace missing covariate values with mean (0)
site.burnt.rich <- site.covs.rich$BurntOrUnburnt
for(i in 1:ncol(site.covs.rich)){
  site.covs.rich[, i] <- (site.covs.rich[, i] - mean(site.covs.rich[, i], na.rm=T)) / (2*sd(site.covs.rich[, i], na.rm=T)) # centre site covariates
}
site.covs.rich$BurntOrUnburnt <- site.burnt.rich # put binary burnt or not burn variable back in
head(site.covs.rich)

#### COMPILE DATA FOR FITTING MODELS IN UNMARKED #####

## P. kundagungan
kund.dat = unmarkedMultFrame(y=dets.kund, siteCovs=site.covs.kund, obsCovs=list(AirTemp=data.frame(obs.covs.kund)),
                              numPrimary=2)
summary(kund.dat)

### P. richmondensis
rich.dat <- unmarkedMultFrame(y=dets.rich[1:37, ], siteCovs=site.covs.rich[1:37, ], obsCovs=list(AirTemp=data.frame(obs.covs.rich[1:37, ])),
                                  numPrimary=3)
summary(rich.dat)

#### FIT MODELS FOR P. KUNDAGUNGAN #####

## Define the model structures to try (up to 3 vars for extinction)
mod.vars <- c('Elevation', 'BurntOrUnburnt', 'GeebamFireSeverity', 'AverageMuddiness', 'DroughtSeverityNd2', 'DroughtSeverityVARI')
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
mod.struc <- c(first.coms, second.coms.vect, third.coms.vect)

## Take out models that include correlated variables
x <- unique(c(grep("BurntOrUnburnt+GeebamFireSeverity", mod.struc, fixed=TRUE),
              grep("AverageMuddiness+DroughtSeverityVARI", mod.struc, fixed=TRUE),
              grep("AverageMuddiness+DroughtSeverityNd2", mod.struc, fixed=TRUE),
              grep("DroughtSeverityNd2+DroughtSeverityVARI", mod.struc, fixed=TRUE)))
mod.struc <- mod.struc[-x]

## Take out a models that doesn't converge
mod.struc <- mod.struc[-which(mod.struc=='Elevation+DroughtSeverityVARI')]
mod.struc <- mod.struc[-which(mod.struc=='Elevation+BurntOrUnburnt+DroughtSeverityVARI')]
mod.struc <- mod.struc[-which(mod.struc=='DroughtSeverityVARI')]
length(mod.struc) # 20 models remaining

## Model function (to run all of the models in 'mod.struc' at once)
mod.func <- function(mod.num, dataset){
  # Now identify the variables that are in each model and specify the model format accordingly
  ext.str <- mod.struc[mod.num]
  # Fit the model
  mod <- colext(psiformula = ~Elevation, 
                gammaformula = ~1,
                epsilonformula = as.formula(paste('~',ext.str, sep="")),
                pformula = ~AirTemp,
                data=dataset)
  # Gather model selection statistics
  aic <- aictab(list(mod), modnames=FALSE)
  # Gather any error messages
  err <- ifelse(mod@opt$convergence==0, 'None', "Failed to converge")
  # Gather up information for the model selection results
  return(c(mod.num, ext.str, aic$K, aic$LL, aic$AICc, err))
}

### FIT MODELS FOR P. KUNDAGUNGAN ###

## Fit each of the models in the set and return model selection statistics
kund.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(kund.mod.mat) <- c("Model", "ExtinctionVariables", "Nparameters", "LogLik", "AIC", "Messages")
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
write.csv(mod.mat.dframe, file='Model_selection_statistics_P_kundagungan_multiseason_occupancy.csv')

# Refit top model
mod <- colext(psiformula = ~Elevation, 
              gammaformula = ~1,
              epsilonformula = ~BurntOrUnburnt+AverageMuddiness,
              pformula = ~AirTemp,
              data=kund.dat)
summary(mod)

# Get confidence intervals for the parameters
confint(mod, type='psi')
confint(mod, type='ext')
confint(mod, type='col')
confint(mod, type='det')

# Estimate occupancy each year
m1 <- nonparboot(mod, B = 100)
occ.pred <- data.frame(Year = c(1:2), MeanOcc = smoothed(mod)[2,], SE = m1@smoothed.mean.bsse[2,])
occ.pred
write.csv(occ.pred, file='P_kundagungan_occ_preds.csv')

### Ext and col predictions for covariate figures

## Extinction vs burn status
burnt <- data.frame(BurntOrUnburnt = c(0,1), AverageMuddiness = rep(mean(kund.dat@siteCovs$AverageMuddiness),2))
avmuddy.unburnt <- data.frame(BurntOrUnburnt=rep(0,26), AverageMuddiness=(seq(0, 100, 4)-62)/(2*30.37526)) 
avmuddy.burnt <- data.frame(BurntOrUnburnt=rep(1,26), AverageMuddiness=(seq(0, 100, 4)-62)/(2*30.37526)) 
burnt.pred <- predict(mod, type='ext', newdata=burnt)
write.csv(burnt.pred, "P_kundagungan_ext_vs_burnt.csv")
burnt.avmuddy <- predict(mod, type='ext', newdata=avmuddy.burnt)
write.csv(burnt.avmuddy, "P_kundagungan_ext_vs_avMuddy_burnt.csv")
unburnt.avmuddy <- predict(mod, type='ext', newdata=avmuddy.unburnt)
write.csv(unburnt.avmuddy, "P_kundagungan_ext_vs_avMuddy_unburnt.csv")

#### FIT MODELS FOR P. RICHMONDENSIS #####

## Define the model structures to try (up to 3 vars for extinction)
mod.vars <- c('Elevation', 'DroughtSeverityNd2', 'DroughtSeverityVARI')
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
mod.struc <- c(first.coms, second.coms.vect, third.coms.vect)

## Take out models that include correlated variables
x <- unique(c(grep("AverageMuddiness+DroughtSeverityVARI", mod.struc, fixed=TRUE),
              grep("AverageMuddiness+DroughtSeverityNd2", mod.struc, fixed=TRUE),
              grep("DroughtSeverityNd2+DroughtSeverityVARI", mod.struc, fixed=TRUE)))
mod.struc <- mod.struc[-x]
length(mod.struc) # 5 models remaining

## Model function (to run all of the models in 'mod.struc' at once)
mod.func <- function(mod.num, dataset){
  # Now identify the variables that are in each model and specify the model format accordingly
  ext.str <- mod.struc[mod.num]
  # Fit the model
  mod <- colext(psiformula = ~Elevation, 
                gammaformula = ~AverageMuddiness,
                epsilonformula = as.formula(paste('~',ext.str, sep="")),
                pformula = ~AirTemp,
                data=dataset)
  # Gather model selection statistics
  aic <- aictab(list(mod), modnames=FALSE)
  # Gather any error messages
  err <- ifelse(mod@opt$convergence==0, 'None', "Failed to converge")
  # Gather up information for the model selection results
  return(c(mod.num, ext.str, aic$K, aic$LL, aic$AICc, err))
}

## Fit each of the models in the set and return model selection statistics
rich.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(rich.mod.mat) <- c("Model", "ExtinctionVariables", "Nparameters", "LogLik", "AIC", "Messages")
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
write.csv(mod.mat.dframe, file='Model_selection_statistics_P_richmondensis_multiseason_occupancy.csv')

# Refit top model
mod <- colext(psiformula = ~Elevation, 
              gammaformula = ~AverageMuddiness,
              epsilonformula = ~DroughtSeverityVARI,
              pformula = ~AirTemp,
              data=rich.dat)
summary(mod)

# Get confidence intervals for the parameters
confint(mod, type='psi')
confint(mod, type='ext')
confint(mod, type='col')
confint(mod, type='det')

# Estimate occupancy each year
m1 <- nonparboot(mod, B = 100)
occ.pred <- data.frame(Year = c(1:3), MeanOcc = smoothed(mod)[2,], SE = m1@smoothed.mean.bsse[2,])
occ.pred
write.csv(occ.pred, file='P_richmondensis_occ_preds.csv')

### Ext and col predictions for covariate figures

## Extinction vs VARI
vari <- data.frame(DroughtSeverityVARI = (seq(0.1, 0.60, 0.02)-0.3762899)/(2*0.1206517)) 
vari.pred <- predict(mod, type='ext', newdata=vari)
write.csv(vari.pred, "P_richmondensis_ext_vs_vari.csv")

## Colonisation vs AvMuddy
avmuddy <- data.frame(AverageMuddiness = (seq(0,100,4)-68.53667)/(2*29.26409)) 
avmuddy.pred <- predict(mod, type='col', newdata=avmuddy)
write.csv(avmuddy.pred, "P_richmondensis_col_vs_avmuddy.csv")