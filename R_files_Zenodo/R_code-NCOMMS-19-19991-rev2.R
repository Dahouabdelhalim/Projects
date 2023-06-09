# All numerical analyses were performed using R open source software (version 3.6.0).

# Analysis of Mortality vs. Temperature

#Datasets:
  # 3 general datasets:
    # (1)   Dataset with tropical and sub-tropical host species infected with bacteria
    # (2)   Dataset with temperate host species infected with bacteria
    # (3)   Dataset with hosts infected with virus

  # 14 host(phylum) and pathogen (family) specific datasets related to the host habitat (tropical-subtropical or temperate):
    # (1)   Fish-Aeromonadaceae (Aeromonas)(tropical-subtropical)
    # (2)   Fish-Aeromonadaceae (Aeromonas) (temperate)
    # (3)   Fish-Enterobacteriaceae (Edwardsiella only) (tropical-subtropical)
    # (4)   Fish-Enterobacteriaceae (Edwardsiella only) (temperate)
    # (5)   Fish-Flavobacteriaceae (Flavobacterium columnare) (tropical-subtropical)
    # (6)   Fish-Flavobacteriaceae (F. columnare) (temperate)
    # (7)   Fish - Strepotoccocae (Streptococcus, Lactococcus) (tropical-subtropical)
    # (8)   Arthropods-Vibrionaceae (Vibrio) (tropical-subtropical)
    # (9)   Fish - Vibrionaceae (Vibrio) (tropical-subtropical)
    # (10)  Molluscs - Vibrionaceae (Vibrio) (tropical-subtropical)
    # (11)  Molluscs - Vibrionaceae (Vibrio) (temperate)
    # (12)  Fish - Yersiniaceae (Yersinia ruckerii) (temperate)
    # (13)  Fish - Koi herpesvirus
    # (14)  Molluscs - Ostreid herpesvirus

#reading dataset:
data<-read.csv("All-temperate-bacteria.csv", header=T, sep=";") #example on reading the general dataset (2)

#packages required: lme4, MuMIn, 
library(lme4)
library(MuMIn)

# 1.1. Selection of random effects (re)
# for both general and specific datasets
lm1<-lm(Mortality~T,data=data) #without re

lmm1<-lmer(Mortality~T+(1 | Reference),data=data,REML=FALSE) #without host or pathogen re
lmm1a<-lmer(Mortality~T+(1 | Pathogen.species)+(1 | Reference),data=data, REML=FALSE) #Pathogen species as re
lmm1b<-lmer(Mortality~T+(1 | Host.family)+(1 | Reference),data=data, REML=FALSE) #Host family as re
lmm1c<-lmer(Mortality~T+(1 | Pathogen.species)+(1 | Host.family)+(1 | Reference),data=data, REML=FALSE) #Both pathogen species and host family as re

#only in general datasets where several host phylums and pathogen families were combined
lmm1d<-lmer(Mortality~T+(1 | Host.phylum)+(1 | Reference),data=data, REML=FALSE) #Host phylum as re
lmm1e<-lmer(Mortality~T+(1 | Pathogen.family)+(1 | Reference),data=data, REML=FALSE) #Pathogen family as re
lmm1f<-lmer(Mortality~T+(1 | Host.phylum/Host.family)+(1 | Reference),data=data, REML=FALSE) 
lmm1g<-lmer(Mortality~T+(1 | Pathogen.family/Pathogen.species)+(1 | Reference),data=data, REML=FALSE)
lmm1h<-lmer(Mortality~T+(1 | Pathogen.family/Pathogen.species)+(1 | Host.family)+(1 | Reference),data=data, REML=FALSE)
lmm1i<-lmer(Mortality~T+(1 | Pathogen.species)+(1 | Host.phylum)+(1 | Reference),data=data, REML=FALSE)

#function for selecting the best model (random effects) using the Akaike weights
model.sel(lmm1,lmm1a,lmm1b,lmm1c)
model.sel(lmm1,lmm1a,lmm1b,lmm1c,lmm1d,lmm1e,lmm1f,lmm1g,lmm1h,lmm1i)

#1.2. Selection of fixed effects (nested linear or mixed effect model)

#re need to be included in following models according to the previously selected models
#example provided for selected model with (1 | Reference) as random effect

lmm1<-lmer(Mortality~T+(1 | Reference),data=data, REML=FALSE)
lmm2a<-lmer(Mortality~T+Host.life.stage+(1 | Reference),data=data, REML=FALSE)
lmm2b<-lmer(Mortality~T+Mode.of.infection+(1 | Reference),data=data, REML=FALSE)
lmm2c<-lmer(Mortality~T+log(Dose)+(1 | Reference),data=data, REML=FALSE)
lmm3a<-lmer(Mortality~T+Host.life.stage+Mode.of.infection+(1 | Reference),data=data, REML=FALSE)
lmm3b<-lmer(Mortality~T+Host.life.stage+log(Dose)+(1 | Reference),data=data, REML=FALSE)
lmm3c<-lmer(Mortality~T+Mode.of.infection+log(Dose)+(1 | Reference),data=data, REML=FALSE)
lmm3d<-lmer(Mortality~T+Mode.of.infection*log(Dose)+(1 | Reference),data=data, REML=FALSE)
lmm4a<-lmer(Mortality~T+Mode.of.infection+log(Dose)+Host.life.stage+(1 | Reference),data=data, REML=FALSE)
lmm4b<-lmer(Mortality~T+Mode.of.infection*log(Dose)+Host.life.stage+(1 | Reference),data=data, REML=FALSE)

model.sel(lmm1,lmm2a,lmm2b,lmm2c,lmm3a,lmm3b,lmm3c,lmm3d,lmm4a,lmm4b)

#if no randonm-effects were necessary:
lm1<-lm(Mortality~T,data=data)
lm2a<-lm(Mortality~T+Host.life.stage,data=data)
lm2b<-lm(Mortality~T+Mode.of.infection,data=data)
lm2c<-lm(Mortality~T+log(Dose),data=data)
lm3a<-lm(Mortality~T+Host.life.stage+Mode.of.infection,data=data)
lm3b<-lm(Mortality~T+Host.life.stage+log(Dose),data=data)
lm3c<-lm(Mortality~T+Mode.of.infection+log(Dose),data=data)
lm3d<-lm(Mortality~T+Mode.of.infection*log(Dose),data=data)
lm4a<-lm(Mortality~T+Mode.of.infection+log(Dose)+Host.life.stage,data=data)
lm4b<-lm(Mortality~T+Mode.of.infection*log(Dose)+Host.life.stage,data=data)
model.sel(lm1,lm2a,lm2b,lm2c,lm3a,lm3b,lm3c,lm3d,lm4a,lm4b)


#check normality of model residues
qqnorm(resid(lmm3d))
#check homoscedascity of variance of model- plot residuals vs fitted
plot(lmm2b)


#Details of selected model:
summary(lmm2c)

#Confidence intervals (95%) of parameters considered in selected:
confint(lmm2c)
#R-squared of selected model
r.squaredGLMM(lmm2c)

#Predicted changes in mortality in response to T° (from previously selected models for tropical and subtropical/temperate species)

predicted<-predict(lmm1,newdata=data)
write.csv(predicted,file="predicted-M-aero-warm-refre.csv")

#before the next step data needs to be re-arranged so each predicted value matches the T from which was predicted
#we create another file containing both of the informations (predicted mortality and Temperature)
data<-read.csv("vibrio-mol-warm-corr-refre.csv", header=T, sep=";")
data1<-read.csv("yers-temp-refre.csv", header=T, sep=";")

library(ggplot2)
#creating graphs of predicted mortality vs T
b<-ggplot(data=data,aes(x=Predicted.mortality,y=T))+
  geom_smooth(method="lm",formula=y~x,se=T,size=0.4,data=data,color="coral2")+
  geom_smooth(method="lm",formula=y~x,se=T,size=0.4,data=data1,color="cornflowerblue")+
  geom_point(data=data,aes(x=Mortality,y=T),alpha=0.27,size=3.5,color="coral2")+
  geom_point(data=data1,aes(x=Mortality,y=T),alpha=0.25,size=4, color="cornflowerblue",pch=18)
  
(b = b + theme_bw()+theme(panel.grid.minor = element_blank(),panel.grid.major = element_blank()))
b + labs (x="Mortality (%)",y="Temperature (°C)")
  

# Analysis of aquaculture-derived MAR and socioeconomic, environmental and health variables
#Dataset Table 2 - Figure 3

# correlation tests between MAR aquaculture and each of the variables:

cor.test(data$MAR.aquaculture,data$Temperature,method="pearson")

# comparison between nested multiple regression models 

lm<-lm(MAR.aquaculture~Temperature+GDP.per.capita+MAR.clinical,data=data)
lm<-lm(MAR.aquaculture~Temperature+GDP.per.capita+MAR.clinical+Temperature*GDP.per.capita,data=data)
lm<-lm(MAR.aquaculture~Temperature+GDP.per.capita+MAR.clinical+Temperature*MAR.clinical,data=data)
lm<-lm(MAR.aquaculture~Temperature+GDP.per.capita+MAR.clinical+GDP.per.capita*MAR.clinical,data=data)
lm<-lm(MAR.aquaculture~Temperature+GDP.per.capita+MAR.clin+Temperature*GDP.per.capita*MAR.clinical,data=data)

summary(lm)

# Lindeman, Merenda and Gold (LMG) metric to calculate variable relative importance 

library(relaimpo)

metrics <- calc.relimp(lm, type = ("lmg"))
metrics
randomlmg <- booteval.relimp(boot.relimp(lm, b = 1000), bty = "perc", level = 0.95)

