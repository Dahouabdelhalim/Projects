aCTI=read.csv(file.choose()) # CTI_Annual.csv

library(MCMCglmm)


# standarising data and calculating weights
aCTI$Weight=(1/aCTI$se^2)
aCTI$Weight=aCTI$Weight/max(aCTI$Weight)
aCTI$nContinent=as.numeric(aCTI$Continent)-mean(as.numeric(aCTI$Continent))

# Analyses

# GLMM ANALYSES

# winter
aWm2.0=blmer(dCTI ~ CTI_1+dTemp+Continent +(1|Country), weights=Weight,data = aCTI[which(aCTI$Season=="Wi"),],REML=FALSE,control = lmerControl(optCtrl = list(maxfun = 20000)))
aWm2.1=blmer(dCTI ~ CTI_1+dTemp+Continent +(1+dTemp|Country), weights=Weight,data = aCTI[which(aCTI$Season=="Wi"),],REML=FALSE,control = lmerControl(optCtrl = list(maxfun = 20000)))

# AIC comparison

library(MuMIn)
library(parameters)

outputctiall<-model.sel(aWm2.0, aWm2.1)
outputctiall

summary(aWm2.1)
parameters::p_value(aWm2.1)

# summer

aBm2.0=blmer(dCTI ~ CTI_1+dTemp+nContinent +(1|Country), weights=Weight,data = aCTI[which(aCTI$Season=="Br"),],REML=FALSE,control = lmerControl(optCtrl = list(maxfun = 20000)))
aBm2.1=blmer(dCTI ~ CTI_1+dTemp+nContinent +(1+dTemp|Country), weights=Weight,data = aCTI[which(aCTI$Season=="Br"),],REML=FALSE,control = lmerControl(optCtrl = list(maxfun = 20000)))

# AIC comparison

outputctiall<-model.sel(aBm2.0, aBm2.1)
outputctiall

summary(aBm2.1)
parameters::p_value(aBm2.1)



