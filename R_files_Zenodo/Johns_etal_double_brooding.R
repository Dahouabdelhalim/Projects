# Mike Johns
# Double brooding analysis
# RMark multistate mark recapture 
# Generalized linear modeling
# January 2018
# Final working script, used for final analysis and writeup

setwd("ENTER PATH INFORMATION")

require(reshape)
require(RMark)
require(devtools)
require(scales)
require(dplyr)
require(gridExtra)
require(lme4)
library(AICcmodavg)

### Survival Analysis ####
#==========================================================================================================
### read in capture histories 
capt.hist<-read.csv("db.capture.histories.csv",colClasses="character")

#=============================================================================================
# state Codes: 
# 1=Pre-breeder, 2=Single Brooded, 3=Double Brooded. 0=not detected/dead

# process data
closed.proc=process.data(data = capt.hist, begin.time = 1984, model = "Multistrata", groups="db")

# Create design data
ms.ddl=make.design.data(closed.proc)

### add mean summer upwelling index ###
upwell<-read.csv("CSV/upwelling.csv")
ms.ddl$S=merge_design.covariates(ms.ddl$S,upwell)
ms.ddl$p=merge_design.covariates(ms.ddl$p,upwell)
ms.ddl$Psi=merge_design.covariates(ms.ddl$Psi,upwell)

#### Define Constraints ####
ms.ddl$Psi$fix=NA

# birds can't survive past age 23
ms.ddl$S$fix[ms.ddl$S$Age == 24]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 25]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 26]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 27]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 28]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 29]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 30]= 0
ms.ddl$S$fix[ms.ddl$S$Age == 31]= 0

ms.ddl$S$fix[ms.ddl$S$age == 24]= 0
ms.ddl$S$fix[ms.ddl$S$age == 25]= 0
ms.ddl$S$fix[ms.ddl$S$age == 26]= 0
ms.ddl$S$fix[ms.ddl$S$age == 27]= 0
ms.ddl$S$fix[ms.ddl$S$age == 28]= 0
ms.ddl$S$fix[ms.ddl$S$age == 29]= 0
ms.ddl$S$fix[ms.ddl$S$age == 30]= 0
ms.ddl$S$fix[ms.ddl$S$age == 31]= 0

# birds can't be detected past age 23
ms.ddl$p$fix[ms.ddl$p$age == 24]= 0
ms.ddl$p$fix[ms.ddl$p$age == 25]= 0
ms.ddl$p$fix[ms.ddl$p$age == 26]= 0
ms.ddl$p$fix[ms.ddl$p$age == 27]= 0
ms.ddl$p$fix[ms.ddl$p$age == 28]= 0
ms.ddl$p$fix[ms.ddl$p$age == 29]= 0
ms.ddl$p$fix[ms.ddl$p$age == 30]= 0
ms.ddl$p$fix[ms.ddl$p$age == 31]= 0

# Examine parameteres to confirm structure
head(ms.ddl$S)
head(ms.ddl$p)
head(ms.ddl$Psi)

#============================================================================================
# Define models and write into function:
run.ms=function()
{
  # S - modeling survival:
  S.1=list(formula= ~ 1)
  S.2=list(formula= ~ stratum)
  S.3=list(formula= ~ stratum + Age)
  S.4=list(formula= ~ stratum + Age + stratum:Age)
  S.5=list(formula= ~ stratum + time)
  S.6=list(formula= ~ stratum + time + Age)
  S.7=list(formula= ~ stratum + time + Age + Age:stratum)
  S.8=list(formula= ~ stratum + time + Age + upall)
  S.9=list(formula= ~ stratum + time + Age + upall + Age:stratum)
  S.10=list(formula= ~ stratum + time + Age + upall + upall:stratum)
  S.11=list(formula= ~ stratum + time + Age + upall + upall:stratum + Age:stratum)
  S.12=list(formula= ~ stratum + time + upall)
  S.13=list(formula= ~ stratum + time + upall + upall:stratum)
  S.19=list(formula= ~ Age)
  S.20=list(formula= ~ time)
  S.21=list(formula= ~ db)
  S.22=list(formula= ~ db + Age)
  S.23=list(formula= ~ db + Age + db:Age)
  S.24=list(formula= ~ db + time)
  S.25=list(formula= ~ db + time + Age)
  S.26=list(formula= ~ db + time + Age + Age:db)
  S.27=list(formula= ~ db + time + Age + upall)
  S.28=list(formula= ~ db + time + Age + upall + Age:db)
  S.29=list(formula= ~ db + time + Age + upall + upall:db)
  S.30=list(formula= ~ db + time + Age + upall + upall:db + Age:db)
  S.31=list(formula= ~ db + time + upall)
  S.32=list(formula= ~ db + time + upall + upall:db)
  S.33=list(formula= ~ time + Age)
  S.34=list(formula= ~ db + stratum)
  S.35=list(formula= ~ db + stratum + Age)
  S.36=list(formula= ~ db + stratum + Age + db:Age)
  S.37=list(formula= ~ db + stratum + time)
  S.38=list(formula= ~ db + stratum + time + Age)
  S.39=list(formula= ~ db + stratum + time + Age + Age:db)
  S.40=list(formula= ~ db + stratum + time + Age + upall)
  S.41=list(formula= ~ db + stratum + time + Age + upall + Age:db)
  S.42=list(formula= ~ db + stratum + time + Age + upall + upall:db)
  S.43=list(formula= ~ db + stratum + time + Age + upall + upall:db + Age:db)
  S.44=list(formula= ~ db + stratum + time + upall)
  S.45=list(formula= ~ db + stratum + time + upall + upall:db) 
  
  # encounter probability: more likely to encounter young birds, in certain years, and 
  # states (pre vs post breeders) - see Townsend & Anderson 2007 Evolution
  p.all=list(formula=~ stratum + Age + time)
  
  # Psi - transition probability:
  # transiton probability depends on interaction of state i and i+1
  Psi.1=list(formula= ~ stratum:tostratum)
  
  # Create model list and run set of models
  ms.model.list=create.model.list("Multistrata")
  ms.results=mark.wrapper(ms.model.list, data=closed.proc,ddl=ms.ddl)

}

# run models: THIS WILL TAKE A WHILE!!!! 
ms.results=run.ms()

#===========================================================================================
# Print AIC model selection results
# with chat correction for fully saturated model
adjust.chat(chat=1.49,ms.results) 

# without chat to check how dispersion effects selection resutls
ms.results

# names of models
names(ms.results)

# store top model in object with short name 
top=ms.results$S.39.p.all.Psi.1

# examine the output from top-ranked model 
PIMS(top,"S",simplified=F)
tail(extract.indices(top,"S"))
top$results$real[500:600,]

betas=summary(top)$beta
betas

#### Lifetime Reproductive Success Analysis ####
#===============================================================================================
# read in data
lrs<-read.csv("CSV/LRS.csv")
lrs$sex<-as.factor(lrs$sex)
lrs$DB_bin<-as.factor(lrs$DB_bin)
lrs$band<-as.factor(lrs$band)

# subset for females ##############################################################
lrs_f<-subset(lrs,sex==5, select=c(hy,LRS,DB,DB_bin,TB,LS,band))

# summary stats
tapply(lrs_f$TB, lrs_f$DB_bin, mean)
#mean 0 = 1.17, 1 = 3.13, 2 = 7
#sd 0 = 1.2, 1 = 1.7, 2 = 3.1

# models
lrs1<-glmer(LRS ~ DB_bin + TB + (1|hy), data=lrs_f,family="poisson")
lrs2<-glmer(LRS ~ DB_bin*TB + (1|hy), data=lrs_f,family="poisson")
lrs3<-glmer(LRS ~ DB_bin + (1|hy), data=lrs_f,family="poisson")
lrs4<-glmer(LRS ~ TB + (1|hy), data=lrs_f,family="poisson")
lrs5<-glmer(LRS ~ DB_bin + LS + (1|hy), data=lrs_f,family="poisson")
lrs6<-glmer(LRS ~ DB_bin*LS + (1|hy), data=lrs_f,family="poisson")
lrs7<-glmer(LRS ~ LS + (1|hy), data=lrs_f,family="poisson")

aictab(list(lrs1,lrs2,lrs3,lrs4,lrs5,lrs6,lrs7))

# top model is DB*TB
summary(lrs2)
plot(lrs2)

# subset for males ##############################################################
lrs_m<-subset(lrs,sex==4, select=c(hy,LRS,DB,DB_bin,TB,LS,band))
lrs_m$DB_bin<-as.factor(lrs_m$DB_bin)

# summary stats 
tapply(lrs_m$LS, lrs_m$DB_bin, sd)
#LRS
#mean 0 = 1.4, 1 = 3.9, 2 = 8.7
#sd 0 = 1.4, 1 = 1.9, 2 = 3.5

# models
lrs1m<-glmer(LRS ~ DB_bin + TB + (1|hy), data=lrs_m,family="poisson")
lrs2m<-glmer(LRS ~ DB_bin*TB + (1|hy), data=lrs_m,family="poisson")
lrs3m<-glmer(LRS ~ DB_bin + (1|hy), data=lrs_m,family="poisson")
lrs4m<-glmer(LRS ~ TB + (1|hy), data=lrs_m,family="poisson")
lrs5m<-glmer(LRS ~ DB_bin + LS + (1|hy), data=lrs_m,family="poisson")
lrs6m<-glmer(LRS ~ DB_bin*LS + (1|hy), data=lrs_m,family="poisson")
lrs7m<-glmer(LRS ~ LS + (1|hy), data=lrs_m,family="poisson")

aictab(list(lrs1m,lrs2m,lrs3m,lrs4m,lrs5m,lrs6m,lrs7m))

summary(lrs2m)
plot(lrs2m)


