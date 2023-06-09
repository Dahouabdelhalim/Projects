####################################################

## CODE FOR FITTING GENERALIZED LINEAR MIXED MODELS TO ASSESS THE IMPACT OF HUMAN DISTURBACE ON 
## NON-BREEDING SHOREBIRDS IN LA PAZ, BAJA CALIFORNIA SUR, MEXICO

## INTIAL OCTOBER 2019
## REVISED FEBRUARY 2022
## CICESE UNIDAD LA PAZ
## LA PAZ, BCS

## AUTHORS - MAC and MER

####

# read in packages that will be used / cargue los programas

library(lme4)
library(glmmTMB)
library(broom)
library(sp)
library(lubridate)
library(reshape2)
library(ggplot2)
library(emmeans)
library(parameters)

## set the working directory  / Establezca el directorio de trabajo
## need to specify
setwd("C:/pfss/msp/mexico_disturbance/") #En mi caso

## read in data set // Lea los datos de la base de SHoreBirD Ensenada de La Paz

ds3<-read.csv("ds3_fin.csv", header=T)


### start analyzing the data

## make the % flooded into proportion flooded
ds3$pinund <- ds3$INUNDADO/100
## make habitat a factor
ds3$habfac <- as.factor(ds3$HABITAT)
## make month a date format
ds3$mes <- month(dmy(ds3$FECHA))
## create an index for season; win=winter and mig=migration
ds3$seas<-ifelse(ds3$mes%in%c(11,12,1,2), "win", "mig")


## make offset to account for different sized sampling units

ds3$lha<-log(ds3$ha)

## scale disturbances to max value

ds3$sraptor<-ds3$Raptorxm/max(ds3$Raptorxm)
ds3$shuman<-ds3$EvntAntrpxm/max(ds3$EvntAntrpxm)


## remove species not used in the analysis (<0.02 prob of occurrence)

ds3<-subset(ds3, ds3$ESPECIE%in%c("AMAV","AMOY","BBPL","BNST","DOWI","DUNL","GRYE","KILL","LBCU", "LESA", "LEYE", "MAGO", "RUTU", "SAND",
                                  "SEPL", "SNPL","SPSA", "WESA", "WHIM", "WILL", "WIPL"))
## all species - migration
## subset migration data
migrac <- subset(ds3,ds3$mes == 8 | ds3$mes == 9 | ds3$mes == 10 | ds3$mes == 3 | ds3$mes == 4)
alshbdm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), data = migrac, ziformula = ~1, family="nbinom2")
summary(alshbdm)
plot(residuals(alshbdm)~fitted(alshbdm))

## looked at the deviations in the random effects
r1<-ranef(alshbdm)
str(r1)
reff<-data.frame(r1$cond$ESPECIE)
names(reff)<-c("intM","shumanM")
plot(reff$shuman~reff$int, xlab = "Deviation in Abundance", 
               ylab = "Deviation in Effect of Disturbance", main = "Migration")

cor.test(reff$intM,reff$shumanM, method="spearman")


## get 95% CI 
alshbdm_mp<-model_parameters(alshbdm)
alshbdm_mp1<-data.frame(cbind(alshbdm_mp$Parameter,alshbdm_mp$Coefficient,alshbdm_mp$CI_low,alshbdm_mp$CI_high))
alshbdm_mp1<-alshbdm_mp1[1:7,]

## create a data frame with all parameters
names(alshbdm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save to working directory
write.csv(alshbdm_mp1, "alshbdm_mp1.csv")

## all species winter
## subset the winter data
invrn <- subset(ds3,ds3$mes == 11 | ds3$mes == 12 | ds3$mes == 1 | ds3$mes == 2)
## run all species model
alshbdi <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = invrn, family="nbinom2")
summary(alshbdi)
plot(residuals(alshbdi)~fitted(alshbdi))

## look at random effect deviations
r2<-ranef(alshbdi)
str(r2)
reff2<-data.frame(r2$cond$ESPECIE)
names(reff2)<-c("intW","shumanW")

## look at effects
reff2

## plot effects
plot(reff2$shuman~reff2$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Winter")

## get 95% CI 
alshbdi_mp<-model_parameters(alshbdi)
alshbdi_mp1<-data.frame(cbind(alshbdi_mp$Parameter,alshbdi_mp$Coefficient,alshbdi_mp$CI_low,alshbdi_mp$CI_high))
alshbdi_mp1<-alshbdi_mp1[1:7,]

## make data frame with parameters
names(alshbdi_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save the file
write.csv(alshbdi_mp1, "alshbdinv_mp1.csv")

# Western Sandpiper in migration
# WESA EN MIGRACION
## subset winter WESA
wesam <- subset(migrac,migrac$ESPECIE == "WESA")
## fit model - note: no need for species random effects
wesamm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD), ziformula = ~1, data = wesam, family="nbinom2")
summary(wesamm)
plot(residuals(wesamm)~fitted(wesamm))

## get 95% CI 
wesamm_mp<-model_parameters(wesamm)
wesamm_mp1<-data.frame(cbind(wesamm_mp$Parameter,wesamm_mp$Coefficient,wesamm_mp$CI_low,wesamm_mp$CI_high))
wesamm_mp1<-wesamm_mp1[1:7,]
## make data frame of parameters
names(wesamm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save the file
write.csv(wesamm_mp1, "wesamm_mp1.csv")


# Western Sandpiper in Winter // WESA EN INVIERNO
## subset winter WESA
wesai <- subset(invrn,invrn$ESPECIE == "WESA")
## fit model
wesaim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD), ziformula = ~1, data = wesai, family="nbinom2")
summary(wesaim)

wesaim_mp<-model_parameters(wesaim)
wesaim_mp1<-data.frame(cbind(wesaim_mp$Parameter,wesaim_mp$Coefficient,wesaim_mp$CI_low,wesaim_mp$CI_high))
wesaim_mp1<-wesaim_mp1[1:7,]
## make data frame for parameters
names(wesaim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save the file
write.csv(wesaim_mp1, "wesaim_mp1.csv")

# Large shorebirds in migration
# GRANDES EN MIGRACION
## subset large species in winter
grndsm <- subset(migrac,migrac$ESPECIE == "AMAV" | migrac$ESPECIE == "AMOY" | migrac$ESPECIE == "BNST" | migrac$ESPECIE == "LBCU" | migrac$ESPECIE == "MAGO" | migrac$ESPECIE == "WHIM" | migrac$ESPECIE == "WILL" )
grndsmm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = grndsm, family="nbinom2")
summary(grndsmm)

##look at random effect deviations
r3<-ranef(grndsmm)
str(r3)
reff3<-data.frame(r3$cond$ESPECIE)
names(reff3)<-c("intM","shumanM")

## look at deviations
reff3

## get 95% CI 
grndsmm_mp<-model_parameters(grndsmm)
grndsmm_mp1<-data.frame(cbind(grndsmm_mp$Parameter,grndsmm_mp$Coefficient,grndsmm_mp$CI_low,grndsmm_mp$CI_high))
grndsmm_mp1<-grndsmm_mp1[1:7,]
## make data frame of parameters
names(grndsmm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(grndsmm_mp1, "grndsmmig_mp1.csv")

# Larger shorebirds in winter
# GRANDES EN INVIERNO
## subset for large body birds in winter
grndsi <- subset(invrn,invrn$ESPECIE == "AMAV" | invrn$ESPECIE == "AMOY" | invrn$ESPECIE == "BNST" | invrn$ESPECIE == "LBCU" | invrn$ESPECIE == "MAGO" | invrn$ESPECIE == "WHIM" | invrn$ESPECIE == "WILL" )
## fit model
grndsim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = grndsi, family="nbinom2")
summary(grndsim)

## look at deviations of random effects
r4<-ranef(grndsim)
str(r4)
reff4<-data.frame(r4$cond$ESPECIE)
names(reff4)<-c("intW","shumanW")

reff4

## get 95% CI 
grndsim_mp<-model_parameters(grndsim)
grndsim_mp1<-data.frame(cbind(grndsim_mp$Parameter,grndsim_mp$Coefficient,grndsim_mp$CI_low,grndsim_mp$CI_high))
grndsim_mp1<-grndsim_mp1[1:7,]
## make data frame of parameters
names(grndsim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(grndsim_mp1, "grndsminv_mp1.csv")

## Medium shorebirds in migration
#MEDIANAS EN MIGRACION
## subset for medium shorebirds in migration
mdnsm <- subset(migrac, migrac$ESPECIE == "BBPL" | migrac$ESPECIE == "DOWI" | migrac$ESPECIE == "GRYE" | migrac$ESPECIE == "KILL" | migrac$ESPECIE == "LEYE" | migrac$ESPECIE == "RUTU")
## fit model
mdnsmm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(1|ESPECIE), ziformula = ~1, data = mdnsm, family="nbinom2")
summary(mdnsmm)

## look at deviations of the random effects
r5<-ranef(mdnsmm)
str(r5)
reff5<-data.frame(r5$cond$ESPECIE)
names(reff5)<-c("intM")
reff5

## get 95% CI 
mdnsm_mp<-model_parameters(mdnsmm)
mdnsm_mp1<-data.frame(cbind(mdnsm_mp$Parameter,mdnsm_mp$Coefficient,mdnsm_mp$CI_low,mdnsm_mp$CI_high))
mdnsm_mp1<-mdnsm_mp1[1:7,]
## make data frame of parameters
names(mdnsm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(mdnsm_mp1,"mdnsmig_mp1.csv")

## Medium shorebirds in winter
# MEDIANAS EN INVIERNO
## subset medium birds in winter
mdnsi <- subset(invrn,invrn$ESPECIE == "BBPL" | invrn$ESPECIE == "DOWI" | invrn$ESPECIE == "GRYE" | invrn$ESPECIE == "KILL" | invrn$ESPECIE == "LEYE" | invrn$ESPECIE == "RUTU")
## fit model
mdnsim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = mdnsi, family="nbinom2")
summary(mdnsim)

## look at deviations of random effects
r6<-ranef(mdnsim)
str(r6)
reff6<-data.frame(r6$cond$ESPECIE)
names(reff6)<-c("intW","shumanW")
reff6

plot(reff6$shuman~reff6$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Medium - Winter")

## get 95% CI 
mdnsim_mp<-model_parameters(mdnsim)
mdnsim_mp1<-data.frame(cbind(mdnsim_mp$Parameter,mdnsim_mp$Coefficient,mdnsim_mp$CI_low,mdnsim_mp$CI_high))
mdnsim_mp1<-mdnsim_mp1[1:7,]
## make data frame of parameters
names(mdnsim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(mdnsim_mp1,"mdnsminv_mp1.csv")

# Small shorebirds in migration
# PEQUEÑAS EN MIGRACION
## subset small birds in migration
pqnsm <- subset(migrac,migrac$ESPECIE == "DUNL" | migrac$ESPECIE == "LESA" | migrac$ESPECIE == "SAND" | migrac$ESPECIE == "SEPL" | migrac$ESPECIE == "SNPL" | migrac$ESPECIE == "SPSA"  | migrac$ESPECIE == "WESA" | migrac$ESPECIE == "WIPL")
## fit model
pqnsmm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = pqnsm, family="nbinom2")
summary(pqnsmm)

## look at deviations of random effects
r7<-ranef(pqnsmm)
str(r7)
reff7<-data.frame(r7$cond$ESPECIE)
names(reff7)<-c("intM","shumanM")
reff7

plot(reff7$shuman~reff7$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Small - Migration")

## get 95% CI 
pqnsmm_mp<-model_parameters(pqnsmm)
pqnsmm_mp1<-data.frame(cbind(pqnsmm_mp$Parameter,pqnsmm_mp$Coefficient,pqnsmm_mp$CI_low,pqnsmm_mp$CI_high))
pqnsmm_mp1<-pqnsmm_mp1[1:7,]
## make data frame of parameters
names(pqnsmm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(pqnsmm_mp1, "pqnsmmig_mp1.csv")

# Small shorebirds in winter
# PEQUEÑAS EN INVIERNO
## subset small birds in winter
pqnsi <- subset(invrn,invrn$ESPECIE == "DUNL" | invrn$ESPECIE == "LESA" | invrn$ESPECIE == "SAND" | invrn$ESPECIE == "SEPL" | invrn$ESPECIE == "SNPL" | invrn$ESPECIE == "SPSA"  | invrn$ESPECIE == "WESA" | invrn$ESPECIE == "WIPL")
## fit model
pqnsim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = pqnsi, family="nbinom2")
summary(pqnsim)

## look at deviations of random effects
r8<-ranef(pqnsim)
str(r8)
reff8<-data.frame(r8$cond$ESPECIE)
names(reff8)<-c("intW","shumanW")
reff8
plot(reff8$shuman~reff8$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Small - Winter")

## get 95% CI 
pqnsim_mp<-model_parameters(pqnsim)
pqnsim_mp1<-data.frame(cbind(pqnsim_mp$Parameter,pqnsim_mp$Coefficient,pqnsim_mp$CI_low,pqnsim_mp$CI_high))
pqnsim_mp1<-pqnsim_mp1[1:7,]
## make data frame of parameters
names(pqnsim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(pqnsim_mp1, "pqnsim_mp1.csv")

# Probing shorebirds in migration
# SONDEADORAS EN MIGRACION
## subset probing birds in migration
snddrsm <- subset(migrac,migrac$ESPECIE == "AMAV" | migrac$ESPECIE == "DOWI" | migrac$ESPECIE == "DUNL" | migrac$ESPECIE == "LBCU" |  migrac$ESPECIE == "LESA" | migrac$ESPECIE == "MAGO" | migrac$ESPECIE == "SAND" | migrac$ESPECIE == "WESA" | migrac$ESPECIE == "WHIM"  | migrac$ESPECIE == "WILL")
## fit model
snddrsmm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = snddrsm, family="nbinom2")
summary(snddrsmm)

## look at the deviations of the random effects
r9<-ranef(snddrsmm)
str(r9)
reff9<-data.frame(r9$cond$ESPECIE)
names(reff9)<-c("intM","shumanM")
reff9

plot(reff9$shuman~reff9$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Probing - Migration")

## get 95% CI 
snddrsmm_mp<-model_parameters(snddrsmm)
snddrsmm_mp1<-data.frame(cbind(snddrsmm_mp$Parameter,snddrsmm_mp$Coefficient,snddrsmm_mp$CI_low,snddrsmm_mp$CI_high))
snddrsmm_mp1<-snddrsmm_mp1[1:7,]
## make a data frame of parameters
names(snddrsmm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(snddrsmm_mp1, "snddrsmm_mp1.csv")

# Probing shorebirds in winter
# SONDEADORAS EN INVIERNO
## subset probing birds in winter
snddrsi <- subset(invrn,invrn$ESPECIE == "AMAV" | invrn$ESPECIE == "DOWI" | invrn$ESPECIE == "DUNL" | invrn$ESPECIE == "LBCU" |  invrn$ESPECIE == "LESA" | invrn$ESPECIE == "MAGO" | invrn$ESPECIE == "SAND" | invrn$ESPECIE == "WESA" | invrn$ESPECIE == "WHIM"  | invrn$ESPECIE == "WILL")
## fit model
snddrsim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = snddrsi, family="nbinom2")
summary(snddrsim)

## look at deviations of random effects
r10<-ranef(snddrsim)
str(r10)
reff10<-data.frame(r10$cond$ESPECIE)
names(reff10)<-c("intW","shumanW")
reff10
plot(reff10$shuman~reff10$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Probing - Winter")

## get 95% CI 
snddrsim_mp<-model_parameters(snddrsim)
snddrsim_mp1<-data.frame(cbind(snddrsim_mp$Parameter,snddrsim_mp$Coefficient,snddrsim_mp$CI_low,snddrsim_mp$CI_high))
snddrsim_mp1<-snddrsim_mp1[1:7,]
## make data frame of parameters
names(snddrsim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(snddrsim_mp1, "snddrsim_mp1.csv")

# Visual foraging shorebirds in migration
# VISUALES EN MIGRACION
## subset visual foraging birds in migration
vslsm <- subset(migrac,migrac$ESPECIE == "AMOY" | migrac$ESPECIE == "BBPL" | migrac$ESPECIE == "BNST" | migrac$ESPECIE == "GRYE" | migrac$ESPECIE == "KILL" | migrac$ESPECIE == "LEYE" | migrac$ESPECIE == "RUTU" | migrac$ESPECIE == "SEPL" | migrac$ESPECIE == "SNPL" | migrac$ESPECIE == "SPSA"  | migrac$ESPECIE == "WIPL")
## fit model
vslsmm <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = vslsm, family="nbinom2")
summary(vslsmm)

## look at deviations of the random effects
r11<-ranef(vslsmm)
str(r11)
reff11<-data.frame(r11$cond$ESPECIE)
names(reff11)<-c("intM","shumanM")
reff11
plot(reff11$shuman~reff11$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Visual - Migration")

## get 95% CI 
vslsmm_mp<-model_parameters(vslsmm)
vslsmm_mp1<-data.frame(cbind(vslsmm_mp$Parameter,vslsmm_mp$Coefficient,vslsmm_mp$CI_low,vslsmm_mp$CI_high))
vslsmm_mp1<-vslsmm_mp1[1:7,]
## make data frame of parameters
names(vslsmm_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
## save file
write.csv(vslsmm_mp1, "C:/pfss/msp/mexico_disturbance/vslsmm_mp1.csv")

# Visual foraging shorebirds in winter
# VISUALES EN INVIERNO
## subset visual foraing birds in winter
vslsi <- subset(invrn,invrn$ESPECIE == "AMOY" | invrn$ESPECIE == "BBPL" | invrn$ESPECIE == "BNST" | invrn$ESPECIE == "GRYE" | invrn$ESPECIE == "KILL" | invrn$ESPECIE == "LEYE" | invrn$ESPECIE == "RUTU" | invrn$ESPECIE == "SEPL" | invrn$ESPECIE == "SNPL" | invrn$ESPECIE == "SPSA"  | invrn$ESPECIE == "WIPL")
## fit the model
vslsim <- glmmTMB(Conteo ~ lha + habfac + pinund + shuman  + sraptor + (1|UNIDAD)+(shuman|ESPECIE), ziformula = ~1, data = vslsi, family="nbinom2")
summary(vslsim)

## look at deviations of the random effects
r12<-ranef(vslsim)
str(r12)
reff12<-data.frame(r12$cond$ESPECIE)
names(reff12)<-c("intW","shumanW")
reff12
plot(reff12$shuman~reff12$int, xlab = "Deviation in Abundance", 
     ylab = "Deviation in Effect of Disturbance", main = "Visual - Winter")

## get 95% CI 
vslsim_mp<-model_parameters(vslsim)
vslsim_mp1<-data.frame(cbind(vslsim_mp$Parameter,vslsim_mp$Coefficient,vslsim_mp$CI_low,vslsim_mp$CI_high))
vslsim_mp1<-vslsim_mp1[1:7,]
## make data frame of parameters
names(vslsim_mp1)<-c("Parameter","Estimate","95CI Lower","95CI Upper")
write.csv(vslsim_mp1, "vslsim_mp1.csv")
