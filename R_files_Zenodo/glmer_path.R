setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

pathUR<-read.csv("pathogens_2015UR.csv")
pathR<-read.csv("pathogens_2015R.csv")

library(piecewiseSEM)
library(lme4)
library(nlme)
library(olsrr)

#CORRELATIONS AMONG LANDCOVER TYPES
pathUR$CROP_500<-pathUR$AG_500-pathUR$PASTURE_500
pathUR$CROP_1250<-pathUR$AG_1250-pathUR$PASTURE_1250
xUR<-pathUR[1:11,1:44] #same for unresolved network

#Croplands correlated with agricultural cover
summary(lm(xUR$CROP_500~xUR$AG_500))
summary(lm(xUR$CROP_750~xUR$AG_750))
summary(lm(xUR$CROP_1250~xUR$AG_1250))

#Natural area correlated with agricultural cover
summary(lm(xUR$NATURAL_500~xUR$AG_500))
summary(lm(xUR$NATURAL_750~xUR$AG_750))
summary(lm(xUR$NATURAL_1250~xUR$AG_1250))

#Pasture lands correlated with agricultural cover
summary(lm(xUR$PASTURE_500~xUR$AG_500))
summary(lm(xUR$PASTURE_750~xUR$AG_750))
summary(lm(xUR$PASTURE_1250~xUR$AG_1250))

#The network and landcover data are the same across each species within a site
#the unique site level values are in the first 11 rows so those are called in the first two models below

## UNRESOLVED NETWORK
#CROPLAND
## 500m
mod.list=psem(lm(connectance~AG_500,data = pathUR[1:11,1:42] ),
              lm(bimpdegree.z~AG_500,data = pathUR[1:11,1:42] ), 
              lm(modularity.z~AG_500,data = pathUR[1:11,1:42] ), 
              lm(WNODFZ~AG_500,data = pathUR[1:11,1:42] ), 
              lm(ABN~AG_500,data = pathUR[1:11,1:42] ), 
              lm(bee.species~AG_500,data = pathUR[1:11,1:42] ), 
              glmer(pathogen~AG_500 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathUR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

##750m
mod.list=psem(lm(connectance~AG_750,data = pathUR[1:11,1:42] ),
              lm(bimpdegree.z~AG_750,data = pathUR[1:11,1:42] ), 
              lm(modularity.z~AG_750,data = pathUR[1:11,1:42] ), 
              lm(WNODFZ~AG_750,data = pathUR[1:11,1:42] ), 
              lm(ABN~AG_750,data = pathUR[1:11,1:42] ), 
              lm(bee.species~AG_750,data = pathUR[1:11,1:42] ), 
              glmer(pathogen~AG_750 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathUR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

##1250
mod.list=psem(lm(connectance~AG_1250,data = pathUR[1:11,1:42] ),
              lm(bimpdegree.z~AG_1250,data = pathUR[1:11,1:42] ), 
              lm(modularity.z~AG_1250,data = pathUR[1:11,1:42] ), 
              lm(WNODFZ~AG_1250,data = pathUR[1:11,1:42] ), 
              lm(ABN~AG_1250,data = pathUR[1:11,1:42] ), 
              lm(bee.species~AG_1250,data = pathUR[1:11,1:42] ), 
              glmer(pathogen~AG_1250 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathUR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

## only the 500m scale predicted network metrics (though all models had the same AIC value, because the GLMER is the same across the models and has 575 data points compared to 11, so differences are obscured)

### model simplification for cropland at 500 m (removed no-significant terms)
mod.list=psem(lm(bimpdegree.z~AG_500,data = pathUR[1:11,1:42]), 
              lm(ABN~AG_500,data = pathUR[1:11,1:42]), 
              lm(connectance~AG_500,data = pathUR[1:11,1:42]),
              glmer(pathogen~AG_500 + ABN + connectance+bimpdegree.z+ (1|site) + (1|Gen.sp), data=pathUR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% ABN,
              bimpdegree.z %~~% ABN)

summary(mod.list) ## not including abundance or AG_500 resulted in a significant d-sep test, so were included. Also, if path between connectance and landscape not evaluated model assumed that connectance was not a site-level metric and correlations with B.imp degree and ABN had 575 df instead of 11, so was also included.

#Does excluding honey bees influence results?
nohbUR<-subset(pathUR, !Gen.sp=="Apis.mellifera")
mod.list=psem(lm(bimpdegree.z~AG_500,data = nohbUR[1:11,1:42]), 
              lm(ABN~AG_500,data = nohbUR[1:11,1:42]), 
              lm(connectance~AG_500,data = nohbUR[1:11,1:42]),
              glmer(pathogen~AG_500 + ABN + connectance+bimpdegree.z+ (1|site) + (1|Gen.sp), data=nohbUR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% ABN,
              bimpdegree.z %~~% ABN)
summary(mod.list)

#TESTING MODEL ASSUMPTIONS 

model1 <- lm(bimpdegree.z~AG_500, data=xUR)
model2 <- lm(connectance~AG_500, data=xUR)
model3 <- lm(ABN~AG_500, data=xUR)

#normality
ols_plot_resid_qq(model1)
ols_plot_resid_qq(model2)
ols_plot_resid_qq(model3)

ols_test_normality(model1)
ols_test_normality(model2)
ols_test_normality(model3)

#Homogeneity of variance
ols_test_breusch_pagan(model1)
ols_test_breusch_pagan(model2)
ols_test_breusch_pagan(model3)

## RESOLVED NETWORK
#CROPLAND
## 500m
mod.list=psem(lm(connectance~AG_500,data = pathR[1:11,1:42] ),
              lm(bimpdegree.z~AG_500,data = pathR[1:11,1:42] ), 
              lm(modularity.z~AG_500,data = pathR[1:11,1:42] ), 
              lm(WNODFZ~AG_500,data = pathR[1:11,1:42] ), 
              lm(ABN~AG_500,data = pathR[1:11,1:42] ), 
              lm(bee.species~AG_500,data = pathR[1:11,1:42] ), 
              glmer(pathogen~AG_500 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

##750m
mod.list=psem(lm(connectance~AG_750,data = pathR[1:11,1:42] ),
              lm(bimpdegree.z~AG_750,data = pathR[1:11,1:42] ), 
              lm(modularity.z~AG_750,data = pathR[1:11,1:42] ), 
              lm(WNODFZ~AG_750,data = pathR[1:11,1:42] ), 
              lm(ABN~AG_750,data = pathR[1:11,1:42] ), 
              lm(bee.species~AG_750,data = pathR[1:11,1:42] ), 
              glmer(pathogen~AG_750 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

##1250
mod.list=psem(lm(connectance~AG_1250,data = pathR[1:11,1:42] ),
              lm(bimpdegree.z~AG_1250,data = pathR[1:11,1:42] ), 
              lm(modularity.z~AG_1250,data = pathR[1:11,1:42] ), 
              lm(WNODFZ~AG_1250,data = pathR[1:11,1:42] ), 
              lm(ABN~AG_1250,data = pathR[1:11,1:42] ), 
              lm(bee.species~AG_1250,data = pathR[1:11,1:42] ), 
              glmer(pathogen~AG_1250 + connectance+bimpdegree.z+modularity.z+ WNODFZ + ABN +bee.species +(1|site) + (1|Gen.sp), data=pathR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% modularity.z,
              connectance %~~% WNODFZ,
              connectance %~~% ABN,
              connectance %~~% bee.species,
              bimpdegree.z %~~% modularity.z,
              bimpdegree.z %~~% WNODFZ,
              bimpdegree.z %~~% ABN,
              bimpdegree.z %~~% bee.species,
              modularity.z %~~% WNODFZ,
              modularity.z %~~% ABN,
              modularity.z %~~% bee.species,
              WNODFZ %~~% ABN,
              WNODFZ %~~% bee.species,
              ABN %~~% bee.species)

summary(mod.list)

## only the 500m scale predicted network metrics (though all models had the same AIC value, because the GLMER is the same across the models and has 575 data points compared to 11, so differences are obscured)

### model simplification for cropland at 500 m (removed no-significant terms)
mod.list=psem(lm(bimpdegree.z~AG_500,data = pathR[1:11,1:42]), 
              lm(ABN~AG_500,data = pathR[1:11,1:42]), 
              lm(connectance~AG_500,data = pathR[1:11,1:42]),
              glmer(pathogen~AG_500 + ABN + connectance+bimpdegree.z+(1|site) + (1|Gen.sp), data=pathR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% ABN,
              bimpdegree.z %~~% ABN)

summary(mod.list) ## not including abundance resulted in a significant d-sep test, so was included. Also, if path between connectance and landscape not evaluated model assumed that connectance was not a site-level metric and correlations with B.imp degree and ABN had 575 df instead of 11, so was also included.

#Does excluding honey bees influence results?
nohbR<-subset(pathR, !Gen.sp=="Apis.mellifera")
mod.list=psem(lm(bimpdegree.z~AG_500,data = nohbR[1:11,1:42]), 
              lm(ABN~AG_500,data = nohbR[1:11,1:42]), 
              lm(connectance~AG_500,data = nohbR[1:11,1:42]),
              glmer(pathogen~AG_500 + ABN + connectance+bimpdegree.z+ (1|site) + (1|Gen.sp), data=nohbR, family = binomial),
              connectance %~~% bimpdegree.z,
              connectance %~~% ABN,
              bimpdegree.z %~~% ABN)
summary(mod.list)

#TESTING MODEL ASSUMPTIONS 
xR<-pathR[1:11,1:42]

model1 <- lm(bimpdegree.z~AG_500, data=xR)
model2 <- lm(connectance~AG_500, data=xR)
model3 <- lm(ABN~AG_500, data=xR)

#normality
ols_plot_resid_qq(model1)
ols_plot_resid_qq(model2)
ols_plot_resid_qq(model3)

ols_test_normality(model1)
ols_test_normality(model2)
ols_test_normality(model3)

#Homogeneity of variance
ols_test_breusch_pagan(model1)
ols_test_breusch_pagan(model2)
ols_test_breusch_pagan(model3)

#H2
##resolved
mod.list=psem(lm(bimpdegree.z~AG_500,data = pathR[1:11,1:42]), 
              lm(H2z~AG_500,data = pathR[1:11,1:42]), 
              glmer(pathogen~AG_500 + H2z+bimpdegree.z+(1|site) + (1|Gen.sp), data=pathR, family = binomial),
              H2z %~~% bimpdegree.z)
summary(mod.list)

summary(lm(xR$WNODFZ~xR$H2z))

##unresolved
mod.list=psem(lm(bimpdegree.z~AG_500,data = pathUR[1:11,1:42]), 
              lm(H2z~AG_500,data = pathUR[1:11,1:42]), 
              glmer(pathogen~AG_500 + H2z+bimpdegree.z+(1|site) + (1|Gen.sp), data=pathUR, family = binomial),
              H2z %~~% bimpdegree.z)
summary(mod.list)

summary(lm(xUR$WNODFZ~xUR$H2z))

