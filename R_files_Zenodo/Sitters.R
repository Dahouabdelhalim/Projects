# Sitters et al. 2020. Nutrient availability controls the impact of mammalian herbivores
# on soil carbon and nitrogen pools in grasslands. In press, Global Change Biology.

# First analysis - MuMIn on local controls - requires the file 'soil.RR.CNpool_local_env.txt'
# Second analysis - MuMIn on environmental drivers - requires the file 'soil.RR.CNpool_local_env.txt'

# load libraries
library(lme4)
library(MuMIn)
library(arm)
library(car)

# set working directory
setwd(".")

# read in data
soil <- read.delim('soil.RR.CNpool_local_env.txt')
soil$Fertilization <- factor(soil$Fertilization,levels=c('Unfertilized','Fertilized'))

# In this first analysis, multi-model inference is used to examine which local controls over soil C and N 
# were responsible for changes in soil C and N pools due to herbivore exclusion. The effects of the local 
# controls on the C and N response ratios were modelled as a full LMM with site ID as a random effect. The 
# models include fertilization as a fixed factor to observe any significant interactions between fertilization 
# and the local controls. Multi-model inference uses model averaging based on Akaike???s information criterion 
# (AIC) to arrive at consistent parameter estimates of the most important explanatory variables in the full LMM, 
# by averaging a set of top models which share similarly high levels of parsimony. Top models were defined as 
# those that fell within 4 AIC units of the model with the lowest AIC value. The parameter estimates of this 
# analysis are depicted in Fig. 4.

# MuMIn for C pool withou microbial data
model.Cpool <- lmer(RR_Cpool~Fertilization + RR_live.biom + RR_dead.biom + RR_root.biom +
                      Fertilization*RR_live.biom + Fertilization*RR_dead.biom + Fertilization*RR_root.biom 
                           +(1|site_code),data=soil,REML=FALSE,na.action='na.fail')
stdz.model.Cpool <- standardize(model.Cpool,standardize.y=FALSE) # standardize regression predictors 
model.set.Cpool <- dredge(stdz.model.Cpool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Cpool4 <- get.models(model.set.Cpool,subset=delta<4) # get the top models within 4 AIC units of 
                                                                # the model with the lowest AIC value
Cpool4 <- model.avg(top.models.Cpool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Cpool4)) # get the parameter estimates of the average model

## MuMIn for N pool without microbial data
model.Npool <- lmer(RR_Npool~Fertilization + RR_live.biom + RR_dead.biom + RR_root.biom +
                      Fertilization*RR_live.biom + Fertilization*RR_dead.biom + Fertilization*RR_root.biom 
                    +(1|site_code),data=soil,REML=FALSE,na.action='na.fail')
stdz.model.Npool <- standardize(model.Npool,standardize.y=FALSE) # standardize regression predictors 
model.set.Npool <- dredge(stdz.model.Npool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Npool4 <- get.models(model.set.Npool,subset=delta<4) # get the top models within 4 AIC units of 
                                                                # the model with the lowest AIC value
Npool4 <- model.avg(top.models.Npool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Npool4)) # get the parameter estimates of the average model

# Create dataset including microbial data and no missing values
completeFun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}
soil.microbe <- completeFun(soil, "RR_mic.biom")

# MuMIn for C pool with microbial data
model.Cpool <- lmer(RR_Cpool~Fertilization + RR_live.biom + RR_dead.biom + RR_root.biom + RR_mic.act + RR_mic.biom +
                      Fertilization*RR_live.biom + Fertilization*RR_dead.biom + Fertilization*RR_root.biom + 
                        Fertilization*RR_mic.act + Fertilization*RR_mic.biom 
                          +(1|site_code),data=soil.microbe,REML=FALSE,na.action='na.fail')
stdz.model.Cpool <- standardize(model.Cpool,standardize.y=FALSE) # standardize regression predictors 
model.set.Cpool <- dredge(stdz.model.Cpool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Cpool4 <- get.models(model.set.Cpool,subset=delta<4) # get the top models within 4 AIC units of 
                                                                # the model with the lowest AIC value
Cpool4 <- model.avg(top.models.Cpool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Cpool4)) # get the parameter estimates of the average model

# MuMIn for N pool with microbial data
model.Npool <- lmer(RR_Npool~Fertilization + RR_live.biom + RR_dead.biom + RR_root.biom + RR_mic.act + RR_mic.biom +
                      Fertilization*RR_live.biom + Fertilization*RR_dead.biom + Fertilization*RR_root.biom + 
                      Fertilization*RR_mic.act + Fertilization*RR_mic.biom 
                    +(1|site_code),data=soil.microbe,REML=FALSE,na.action='na.fail')
stdz.model.Npool <- standardize(model.Npool,standardize.y=FALSE) # standardize regression predictors 
model.set.Npool <- dredge(stdz.model.Npool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Npool4 <- get.models(model.set.Npool,subset=delta<4) # get the top models within 4 AIC units of 
                                                                # the model with the lowest AIC value
Npool4 <- model.avg(top.models.Npool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Npool4)) # get the parameter estimates of the average model

####################################################################################################################

# In this first analysis, multi-model inference is used to examine which across-site environmental drivers 
# affect the impact of herbivore exclusion on soil C and N pools. The effects of the envrironmental drivers 
# on the C and N response ratios were modelled as a full LMM with site ID as a random effect. The 
# models include fertilization as a fixed factor to observe any significant interactions between fertilization 
# and the local controls. Multi-model inference uses model averaging based on Akaike???s information criterion 
# (AIC) to arrive at consistent parameter estimates of the most important explanatory variables in the full LMM, 
# by averaging a set of top models which share similarly high levels of parsimony. Top models were defined as 
# those that fell within 4 AIC units of the model with the lowest AIC value. The parameter estimates of this 
# analysis are depicted in Fig. 5.

## MuMIn for C pool
model.Cpool <- lmer(RR_Cpool~Fertilization + MAT + TEMP_VAR + TEMP_WET_Q + MAP + MAP_VAR +
                             N_Dep + plant.biom_yr0 + pct_N_y0 + Fertilization*MAT +
                             Fertilization*TEMP_VAR + Fertilization*TEMP_WET_Q + Fertilization*MAP +
                             Fertilization*MAP_VAR + Fertilization*N_Dep + Fertilization*plant.biom_yr0 + 
                             Fertilization*pct_N_y0 + (1|site_code),data=soil,REML=FALSE,na.action='na.fail')

stdz.model.Cpool <- standardize(model.Cpool,standardize.y=FALSE) # standardize regression predictors 
model.set.Cpool <- dredge(stdz.model.Cpool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Cpool4 <- get.models(model.set.Cpool,subset=delta<4) # get the top models within 4 AIC units of 
# the model with the lowest AIC value
Cpool4 <- model.avg(top.models.Cpool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Cpool4)) # get the parameter estimates of the average model

## MuMIn for N pool
model.Npool <- lmer(RR_Npool~Fertilization + MAT + TEMP_VAR + TEMP_WET_Q + MAP + MAP_VAR +
                      N_Dep + plant.biom_yr0 + pct_N_y0 + Fertilization*MAT +
                      Fertilization*TEMP_VAR + Fertilization*TEMP_WET_Q + Fertilization*MAP +
                      Fertilization*MAP_VAR + Fertilization*N_Dep + Fertilization*plant.biom_yr0 + 
                      Fertilization*pct_N_y0 + (1|site_code),data=soil,REML=FALSE,na.action='na.fail')

stdz.model.Npool <- standardize(model.Npool,standardize.y=FALSE) # standardize regression predictors 
model.set.Npool <- dredge(stdz.model.Npool,extra='R^2') # generate set of models and calculate R2 of each
top.models.Npool4 <- get.models(model.set.Npool,subset=delta<4) # get the top models within 4 AIC units of 
# the model with the lowest AIC value
Npool4 <- model.avg(top.models.Npool4,revised.var = TRUE) # average this set of top models
summary(model.avg(top.models.Npool4)) # get the parameter estimates of the average model
