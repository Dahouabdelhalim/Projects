#--------------------------------------------------------------------------#
#                     COUNT MODELS FOR PHILORIA                            #
#--------------------------------------------------------------------------#

##### REQUIRED LIBRARIES #####
library(tuple)
library(lme4)
library(AICcmodavg)
library(dplyr)
library(readxl)

#### IMPORT AND COMPILE SURVEY DATA #####
setwd('') # Add own
glm.dat <- as.data.frame(read_excel('Philoria_count_data_2021_2022.xlsx', sheet=1, na="NA", col_types=c(rep("text", 1), rep("numeric", 13))))
str(glm.dat)
glm.dat$SiteID <- as.factor(glm.dat$SiteID)
str(glm.dat)

#### COMPLETE MODEL SELECTION IN LME4 #####

### Model function (to run all of the models in 'mod.struc' at once)
mod.func <- function(mod.num, dataset){
  # Now identify the variables that are in each model and specify the model format accordingly
  vars <- mod.struc[mod.num]
  mod.str <- paste('Count ~ ', paste(vars[1:length(vars)], collapse=" + "), "+ (1|SiteID)")
  # Fit the model
  mod <- glmer(mod.str, family = poisson, data=dataset)
  # Gather model selection statistics
  aic <- aictab(list(mod), modnames=FALSE)
  # Gather any error messages
  err <- ifelse(is.null(mod@optinfo$conv$lme4$messages)==TRUE, 'None', mod@optinfo$conv$lme4$messages)
  # Gather up information for the model selection results
  return(c(mod.num, paste(mod.str[1:length(mod.str)], collapse="+"), aic$K, aic$LL, aic$AICc, err))
}

### First define model structures to include in the analysis
names(glm.dat)
mod.vars <- c("Elevation", "BurntOrUnburnt", "PercentBedBurnt", "AverageBankBurnSeverity", "GeebamFireSeverity",     
              "SpinachCover", "AverageMuddiness", "MinimumMuddiness", "DroughtSeverityNd2", "DroughtSeverityVARI", "AirTemp")

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

## Take out models that don't include an effect of air temperature (air temp in all models a priori)
x <- grep("AirTemp", mod.struc)
mod.struc <- mod.struc[x]
length(mod.struc)
mod.struc

## Check for duplicate models
dupes <- duplicated(mod.struc) 
sum(dupes) # none

### P.KUNDAGUNGAN ###

## Compile data
kund.sites <- c("KUN_1",  "KUN_2",  "KUN_3",  "KUN_4",  "KUN_5",  "KUN_6",  "KUN_7",  "KUN_8",  "KUN_9",  "KUN_10", "KUN_11", "KUN_12", "KUN_13", "KUN_14",
                "KUN_15", "KUN_16", "KUN_17", "KUN_18", "KUN_19", "KUN_20", "KUN_21", "KUN_22", "KUN_23", "KUN_24", "KUN_25", "KUN_26", "KUN_27", "KUN_28",
                "KUN_29", "KUN_30", "KUN_31", "KUN_32", "KUN_33", "KUN_34", "KUN_35", "KUN_F1", "KUN_F2", "KUN_F3", "KUN_E1", "9907", "9908", "9909",  
                "9920", "9921", "9922", "9924", "9926", "9927")
kund.match <- matchAll(as.factor(kund.sites), glm.dat$SiteID)
kund.dat <- glm.dat[kund.match, ]
kund.dat$SiteID
head(kund.dat)

## Standardise covariates
scale.func <- function(x){
  z <- (x - mean(x, na.rm=T)) / (2 * sd(x, na.rm=T))
  return(z)
}
names(kund.dat)
kund.dat$AirTemp <- scale.func(kund.dat$AirTemp)
kund.dat$Elevation <- scale.func(kund.dat$Elevation)
kund.dat$PercentBedBurnt <- scale.func(kund.dat$PercentBedBurnt)
kund.dat$AverageBankBurnSeverity <- scale.func(kund.dat$AverageBankBurnSeverity)
kund.dat$GeebamFireSeverity <- scale.func(kund.dat$GeebamFireSeverity)
kund.dat$SpinachCover <- scale.func(kund.dat$SpinachCover)
kund.dat$AverageMuddiness <- scale.func(kund.dat$AverageMuddiness)
kund.dat$MinimumMuddiness <- scale.func(kund.dat$MinimumMuddiness)
kund.dat$DroughtSeverityNd2 <- scale.func(kund.dat$DroughtSeverityNd2)
kund.dat$DroughtSeverityVARI <- scale.func(kund.dat$DroughtSeverityVARI)
kund.dat$Rain2020 <- scale.func(kund.dat$Rain2020)

## Fit each of the models in the set and return model selection statistics
kund.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(kund.mod.mat) <- c("Model", "Variables", "Nparameters", "LogLik", "AIC", "ErrorMessages")
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
write.csv(mod.mat.dframe, file='Count_model_selection_statistics_P_kundagungan_counts.csv')

## Refit top model to look at regression coefficients
mod <- glmer(mod.mat.dframe$Variables[1], family = poisson, data=kund.dat)
summary(mod)

# Confidence intervals for parameters of top model
confint(mod)

### P.RICHMONDENSIS ###

## Compile data
rich.sites <- c("RIC_1", "RIC_2", "RIC_3", "RIC_4", "RIC_5", "RIC_6", "RIC_7", "RIC_8", "RIC_9", "RIC_10", "RIC_11", "RIC_12", "RIC_13", 
                "RIC_14", "RIC_15", "RIC_16", "RIC_17", "RIC_18", "RIC_19", "RIC_20", "RIC_21", "RIC_22", "RIC_23", "RIC_24", "RIC_25", "RIC_26", 
                "RIC_27", "RIC_28", "RIC_29", "RIC_30", "RIC_31", "RIC_32", "RIC_33", "RIC_34", "RIC_35", "RIC_36", "RIC_37", "RIC_F1", "RIC_F2", 
                "RIC_F3", "RIC_F4", "RIC_F5", "RIC_F6", "RIC_F7", "RIC_F8", "RIC_F9", "RIC_F10", "RIC_F11", "RIC_F12", "RIC_F13")
rich.match <- matchAll(as.factor(rich.sites), glm.dat$SiteID)
rich.dat <- glm.dat[rich.match, ]
rich.dat$SiteID
head(rich.dat)

## Standardise covariates
names(rich.dat)
rich.dat$AirTemp <- scale.func(rich.dat$AirTemp)
rich.dat$Elevation <- scale.func(rich.dat$Elevation)
rich.dat$PercentBedBurnt <- scale.func(rich.dat$PercentBedBurnt)
rich.dat$AverageBankBurnSeverity <- scale.func(rich.dat$AverageBankBurnSeverity)
rich.dat$GeebamFireSeverity <- scale.func(rich.dat$GeebamFireSeverity)
rich.dat$SpinachCover <- scale.func(rich.dat$SpinachCover)
rich.dat$AverageMuddiness <- scale.func(rich.dat$AverageMuddiness)
rich.dat$MinimumMuddiness <- scale.func(rich.dat$MinimumMuddiness)
rich.dat$DroughtSeverityNd2 <- scale.func(rich.dat$DroughtSeverityNd2)
rich.dat$DroughtSeverityVARI <- scale.func(rich.dat$DroughtSeverityVARI)
rich.dat$Rain2020 <- scale.func(rich.dat$Rain2020)

## Fit each of the models in the set and return model selection statistics
rich.mod.mat <- matrix(nrow=length(mod.struc), ncol=6)
colnames(rich.mod.mat) <- c("Model", "Variables", "Nparameters", "LogLik", "AIC", "ErrorMessages")
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
write.csv(mod.mat.dframe, file='Count_model_selection_statistics_P_richmondensis_counts.csv')

## Refit top model to look at regression coefficients
mod <- glmer(mod.mat.dframe$Variables[1], family = poisson, data=rich.dat)
summary(mod)

# Confidence intervals for parameters of top model
confint(mod)