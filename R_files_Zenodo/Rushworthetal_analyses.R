### R script for GLMMs and LRTs from Rushworth et al. Evolution Letters
### "Identifying the fitness consequences of sex in complex natural environments" 


### Load in packages

library(glmmTMB)
library(car) #for Levene Tests
library(ggplot2)
library(ggeffects)
library(robustbase) # for identifying outliers
#library(MASS) # for exploring distributions
#library(fitdistrplus) # for exploring distributions
library(pbkrtest)
library(tidyverse)


### Note that p-value corrections were applied using the "p.adjust" command with method="holm" ###


# read in data after converting Excel documents to .csv
data <- read.csv("Rushworthetal_fielddata.csv")

# remove plants that did not survive to be planted in the field (WinGH == NA)
data <- data[which(data$WinGH != "NA"),]

# remove the genotype that turned out to be B. stricta
data <- data[which(data$Geno != "S11"),]

# to replicate the authors' results, gardens that failed to reproduce must also be removed; see publication methods
data <- data[!(data$Cohort == 1 & data$Loc == "MIL"),]
data <- data[!(data$Cohort == 2 & data$Loc == "CUM"),]


# change Cohort (Experimental Year in manuscript) to a factor with values 1 and 2
data$Cohort <- as.factor(data$Cohort)

# add homozygosity to the dataset
hetdata <- read.csv("Rushworthetal_heterozygosity.csv")
data$hom = hetdata[match(data$Geno, hetdata$Code),"homozygosity"]
data$het = 1-data$hom







###############################################################################################################
############################################## Total fitness ##################################################
###############################################################################################################


##############################################  sex vs asex  ##################################################



### GLMM, zero-inflated negative binomial ###

### OUTLIERS ### 

adjbox(data$SeedsmaxY1, range=1.5) # this shows anything above 3000 as an outlier
fitness_outliers <- data[which(data$SeedsmaxY1>3000),] # N=60 obs

# remove outliers
datains <- data[which(data$SeedsmaxY1<3000),] # N=3058 obs



# GLMM
totalw <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|RS:Geno) + (1|Loc:Block),
                  data=datains,
                  ziformula=~.,
                  family="nbinom2")
summary(totalw)
ggemmeans(totalw,~ RS | Cohort)
ggemmeans(totalw,~ RS | Loc)




### LRTs ###

# (1|RS:Geno)
TWran1 <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|Loc:Block),
                   data=datains,
                   zeroInflation=~.,
                   family="nbinom")
anova(TWran1,totalw)


# (1|Loc:Block)
TWran2 <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|RS:Geno),
                   data=datains,
                   zeroInflation=~.,
                   family="nbinom")
anova(TWran1,totalw)


# RS:Loc
### Model 1: RS*Loc (totalw)
### Model 2: RS + Loc (totalwB)

totalwB <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + Cohort + RS:Cohort + (1|RS:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(totalwB, totalw)


# RS:Cohort
### Model 1: RS*Loc + Cohort + RS:Cohort (totalw)
### Model 2: RS*Loc + Cohort (totalwC)

totalwC <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(totalwC, totalw)


# Cohort
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: RS + Loc (totalwD)

totalwFix <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                      data=datains,
                      zeroInflation=~.,
                      family="nbinom")

totalwD <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + (1|RS:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")

anova(totalwD, totalwFix)


# Loc
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: RS + Cohort (totalwE)

totalwE <- glmmTMB(round(SeedsmaxY1) ~ RS + Cohort + (1|RS:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(totalwE, totalwFix)


# RS
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: Loc + Cohort (totalwF)

totalwF <- glmmTMB(round(SeedsmaxY1) ~ Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(totalwF, totalwFix)



########################################  hyb vs non vs sex ("group" term)  ########################################



### GLMM, zero-inflated negative binomial ###


# GLMM
Gtw <- glmmTMB(round(SeedsmaxY1) ~ Group*Loc + Cohort + Group:Cohort + (1|Group:Geno) + (1|Loc:Block),
                  data=datains,
                  ziformula=~.,
                  family="nbinom2")

summary(Gtw)
ggemmeans(Gtw,~ Group | Cohort)
ggemmeans(Gtw,~ Group | Loc)
ggemmeans(Gsurv,terms="Group")
ggemmeans(Gsurv,terms="Loc")





### LRTs ###

# (1|Group:Geno)
GTWr1 <- glmmTMB(round(SeedsmaxY1) ~ Group*Loc + Cohort + Group:Cohort + (1|Loc:Block),
                   data=datains,
                   zeroInflation=~.,
                   family="nbinom")
anova(GTWr1,Gtw)


# (1|Loc:Block)
GTWr2 <- glmmTMB(round(SeedsmaxY1) ~ Group*Loc + Cohort + Group:Cohort + (1|Group:Geno),
                   data=datains,
                   zeroInflation=~.,
                   family="nbinom")
anova(GTWr2,Gtw)


# Group:Loc
### Model 1: Group*Loc (Gtw)
### Model 2: Group + Loc (GtwB)

GtwB <- glmmTMB(round(SeedsmaxY1) ~ Group + Loc + Cohort + Group:Cohort + (1|Group:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(GtwB, Gtw)


# Group:Cohort
### Model 1: Group*Loc + Cohort + Group:Cohort (Gtw)
### Model 2: Group*Loc + Cohort (GtwC)

GtwC <- glmmTMB(round(SeedsmaxY1) ~ Group*Loc + Cohort + (1|Group:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(GtwC, Gtw)


# Cohort
### Model 1: Group + Loc + Cohort (fixed terms only)
### Model 2: Group + Loc (GtwD)

GtwFix <- glmmTMB(round(SeedsmaxY1) ~ Group + Loc + Cohort + (1|Group:Geno) + (1|Loc:Block),
                      data=datains,
                      zeroInflation=~.,
                      family="nbinom")

GtwD <- glmmTMB(round(SeedsmaxY1) ~ Group + Loc + (1|Group:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")

anova(GtwD, GtwFix)


# Loc
### Model 1: Group + Loc + Cohort (fixed terms only)
### Model 2: Group + Cohort (GtwE)

GtwE <- glmmTMB(round(SeedsmaxY1) ~ Group + Cohort + (1|Group:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(GtwE, GtwFix)


# Group
### Model 1: Group + Loc + Cohort (fixed terms only)
### Model 2: Loc + Cohort (GtwF)

GtwF <- glmmTMB(round(SeedsmaxY1) ~ Loc + Cohort + (1|Group:Geno) + (1|Loc:Block),
                    data=datains,
                    zeroInflation=~.,
                    family="nbinom")
anova(GtwF, GtwFix)







###############################################################################################################
###############################################   SURVIVAL   #################################################
###############################################################################################################


##############################################  sex vs asex  ##################################################


surv <- glmmTMB(survY1 ~ RS*Loc + Cohort + RS:Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                data=data,
                zeroInflation=FALSE,
                family="binomial")
summary(surv)

ggemmeans(mort,terms="Loc") 
ggemmeans(mort,pairwise ~ "RS" | "Loc")



### LRTs ###


# (1|RS:Geno)
ran1 <- glmmTMB(survY1 ~ RS*Loc + Cohort + RS:Cohort + scale(WinGH) + (1|Loc:Block),
                 data=data,
                 zeroInflation=FALSE,
                 family="binomial")
anova(ran1,surv)

# (1|Loc:Block)
ran2 <- glmmTMB(mortY1 ~ RS*Loc + Cohort + RS:Cohort + scale(WinGH) + (1|RS:Geno),
                 data=data,
                 zeroInflation=FALSE,
                 family="binomial")
anova(ran2,surv)



# RS:Cohort
### Model 1: RS*Loc + Cohort + RS:Cohort + WinGH (mort) 
### Model 2: RS*Loc + Cohort + WinGH (mortB)
survRC <- glmmTMB(survY1 ~ RS*Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")
anova(survRC,surv)


# RS:Loc
### Model 1: RS*Loc + Cohort + RS:Cohort + WinGH (surv)
### Model 2: RS + Loc + Cohort + RS:Cohort + WinGH (survRC)
survRL <- glmmTMB(survY1 ~ RS + Loc + Cohort + RS:Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")
anova(survRL, surv)


# Cohort
### Model 1: RS + Loc + Cohort + WinGH (survFix) (fixed terms only)
### Model 2: RS + Loc + WinGH (survC)

survFix <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                    data=data,
                    zeroInflation=FALSE,
                    family="binomial")

survC <- glmmTMB(survY1 ~ RS + Loc + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")

anova(survC, survFix)


# Loc 
### Model 1: RS + Loc + Cohort + WinGH (survFix) (fixed terms only)
### Model 2: RS + Cohort + WinGH (survL)

survL <- glmmTMB(survY1 ~ RS + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")

anova(survL, survFix)

# RS 
### Model 1: RS + Loc + Cohort + WinGH (survFix)
### Model 2: Loc + Cohort + WinGH (survR)

survR <- glmmTMB(survY1 ~ Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")

anova(survR, survFix)


# WinGH 
### Model 1: RS + Loc + Cohort + WinGH (survFix)
### Model 2: RS + Loc + Cohort (survW)

survW <- glmmTMB(survY1 ~ RS + Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")
anova(survW,survFix)



########################################  hyb vs non vs sex ("group" term)  ########################################


Gsurv <- glmmTMB(survY1 ~ Group*Loc + Cohort + Group:Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                  data=data,
                  zeroInflation=FALSE,
                  family="binomial")
summary(Gsurv)


ggemmeans(Gsurv,terms="Group")
ggemmeans(Gsurv,terms="Loc")
ggemmeans(Gsurv,~ Group | Loc)
ggemmeans(Gsurv,~ Group | Cohort)




### LRTs ###


### random effects ###
# (1|Group:Geno)
GSran1 <- glmmTMB(survY1 ~ Group*Loc + Cohort + Group:Cohort + scale(WinGH) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GSran1,Gsurv)

# (1|Loc:Block)
GSran2 <- glmmTMB(survY1 ~ Group*Loc + Cohort + Group:Cohort + scale(WinGH) + (1|Group:Geno),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GSran2,Gsurv)


# Group:Cohort
### Model 1: Group*Loc + Cohort + scale(WinGH) (GsurvA)
### Model 2: Group*Loc + Cohort + Group:Cohort + scale(WinGH) (Gsurv)
GsurvA <- glmmTMB(survY1 ~ Group*Loc + Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GsurvA,Gsurv)


# Group:Loc
### Model 1: Group*Loc + Cohort + Group:Cohort + scale(WinGH) (GsurvB)
### Model 2: Group + Loc + Cohort + Group:Cohort + scale(WinGH) (Gsurv)
GsurvB <- glmmTMB(survY1 ~ Group + Loc + Cohort + Group:Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GsurvB,Gsurv)


# WinGH
### Model 1: Group + Loc + Cohort + scale(WinGH) (GsurvFix)
### Model 2: Group + Loc + Cohort (GmortC)
GsurvFix <- glmmTMB(survY1 ~ Group + Loc + Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                       data=data,
                       zeroInflation=FALSE,
                       family="binomial")
GsurvC <- glmmTMB(survY1 ~ Group + Loc + Cohort + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")

anova(GsurvC,GsurvFix)


# Cohort
### Model 1: Group + Loc + Cohort + scale(WinGH) (GsurvFix)
### Model 2: Group + Loc + scale(WinGH) (GsurvD)
GsurvD <- glmmTMB(survY1 ~ Group + Loc + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GsurvD,GsurvFix)


# Loc
### Model 1: Group + Loc + Cohort + scale(WinGH) (GsurvFix)
### Model 2: Group + Cohort + scale(WinGH) (GsurvE)
GsurvE <- glmmTMB(survY1 ~ Group + Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GmortE,GmortFixed)


# Group
### Model 1: Group + Loc + Cohort + scale(WinGH) (GsurvFix)
### Model 2: Loc + Cohort + scale(WinGH) (GsurvF)
GsurvF <- glmmTMB(survY1 ~ Loc + Cohort + scale(WinGH) + (1|Group:Geno) + (1|Loc:Block),
                   data=data,
                   zeroInflation=FALSE,
                   family="binomial")
anova(GsurvF,GsurvFix)








###############################################################################################################
################################################## Fecundity ##################################################
###############################################################################################################



##############################################  sex vs asex  ##################################################


fecdata <- subset(data, Relseedsmax > 0)
length(which(fecdata$Relseedsmax>0)) # 995 obs

# subset to remove outliers 
fecdata <- fecdata[which(fecdata$SeedsmaxY1 <8000),] #989 obs

### GLMM


fecmodel <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|RS:Geno) + (1|Loc:Block),
                data=fecdata,
                family="nbinom2")
summary(fecmodel) 

emmeans(fecmodel, ~ RS | Cohort)


### LRTs


# (1|RS:Geno)
Fran1 <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|Loc:Block),
                 data=fecdata,
                 family="nbinom2")
anova(Fran1, fecmodel)

# (1|Loc:Block)
Fran2 <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + RS:Cohort + (1|RS:Geno),
                 data=fecdata,
                 family="nbinom2")
anova(Fran2, fecmodel)


# RS:Loc
### Model 1: RS*Loc + Cohort + RS:Cohort (fecmodel)
### Model 2: RS + Loc + Cohort + RS:Cohort (fecRL)

fecRL <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + Cohort + RS:Cohort + (1|RS:Geno) + (1|Loc:Block),
                 data=fecdata,
                 family="nbinom2")
anova(fecRL, fecmodel)


# RS:Cohort
### Model 1: RS*Loc + Cohort + RS:Cohort (fec)
### Model 2: RS*Loc + Cohort (fecRC)

fecRC <- glmmTMB(round(SeedsmaxY1) ~ RS*Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                 data=fecdata,
                 family="nbinom2")
anova(fecRC, fecmodel)


# Cohort
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: RS + Loc (fecC)

fecFix <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                  data=fecdata,
                  family="nbinom2")

fecC <- glmmTMB(round(SeedsmaxY1) ~ RS + Loc + (1|RS:Geno) + (1|Loc:Block),
                data=fecdata,
                family="nbinom2")

anova(fecC, fecFix)

# Loc
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: RS + Cohort (fecL)

fecL <- glmmTMB(round(SeedsmaxY1) ~ RS + Cohort + (1|RS:Geno) + (1|Loc:Block),
                data=fecdata,
                family="nbinom2")
anova(fecL, fecFix)

# RS
### Model 1: RS + Loc + Cohort (fixed terms only)
### Model 2: Loc + Cohort (fecR)

fecR <- glmmTMB(round(SeedsmaxY1) ~ Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                data=fecdata,
                family="nbinom2")
anova(fecR, fecFix)







###############################################################################################################
############################################### LEAF HERBIVORY ################################################
###############################################################################################################



##############################################  sex vs asex  ##################################################


library(lme4)
library(lsmeans)
library(pbkrtest)

herbdata <- subset(data, LogDamPC != -999) # N = 1379


# linear model
damage <- lme4::lmer(LogDamPC ~ scale(H) + RS*Loc + Cohort + RS:Cohort + (1|RS:Geno) + (1|Loc:Block), data=herbdata, REML=FALSE)
summary(damage) 

# updating the models to run F tests
damage_noran1 <- update(damage, . ~ . -(1|RS:Geno))
damage_noran2 <- update(damage, . ~ . -(1|Loc:Block))
damage_noRL <- update(damage, . ~ . -RS:Loc)
damage_noRC <- update(damage, . ~ . -RS:Cohort)
damage_fixed <- update(damage, . ~ . -RS:Loc -RS:Cohort)
summary(damage_fixed)
damage_noL <- update(damage_fixed, . ~ . -Loc)
damage_noC <- update(damage_fixed, . ~ . -Cohort)
damage_noR <- update(damage_fixed, . ~ . -RS)
damage_noH <- update(damage_fixed, . ~ . -scale(H))
summary(damage_noR)

# LRT for random effects
anova(damage_noran1, damage) # test effect of (1|RS:Geno)
anova(damage_noran2, damage) # test effect of (1|Loc:Block)


# F test with Kenward-Roger approximation 

(damage.kr.rl <- KRmodcomp(damage,damage_noRL)) # test effect of RS*Loc
(damage.kr.rc <- KRmodcomp(damage,damage_noRC)) # test effect of RS*Cohort
(damage.kr.c <- KRmodcomp(damage_noC, damage_fixed)) # test effect of Cohort
(damage.kr.l <- KRmodcomp(damage_noL, damage_fixed)) # test effect of Loc
(damage.kr.r <- KRmodcomp(damage_noR, damage_fixed)) # test effect of RS
(damage.kr.h <- KRmodcomp(damage_noH, damage_fixed)) # test effect of height! 




########################################  hyb vs non vs sex ("group" term)  ########################################


# linear model 

damG <- lme4::lmer(LogDamPC ~ scale(H) + Group*Loc + Cohort + Group:Cohort + (1|Group:Geno) + (1|Loc:Block), data=herbdata, REML=FALSE)
summary(damG) 

# updating the models to run F tests/LRTs
damG_noran1 <- update(damG, . ~ . -(1|Group:Geno))
damG_noran2 <- update(damG, . ~ . -(1|Loc:Block))
damG_noGL <- update(damG, . ~ . -Group:Loc)
damG_noGC <- update(damG, . ~ . -Group:Cohort)
damG_fixed <- update(damG_noGC, . ~ . -Group:Loc)
summary(damG_fixed)
damG_noL <- update(damG_fixed, . ~ . -Loc)
damG_noC <- update(damG_fixed, . ~ . -Cohort)
damG_noG <- update(damG_fixed, . ~ . -Group)
damG_noH <- update(damG_fixed, . ~ . -scale(H))
summary(damG_noG)


# LRTs for random effects
anova(damG_noran1, damG) # test effect of (1|Group:Geno)
anova(damG_noran2, damG) # test effect of (1|Loc:Block)

# F test with Kenward-Roger approximation 

(damG.kr.gl <- KRmodcomp(damG,damG_noGL)) # test effect of Group*Loc
(damG.kr.gc <- KRmodcomp(damG,damG_noGC)) # test effect of Group*Cohort
(damG.kr.c <- KRmodcomp(damG_noC, damG_fixed)) # test effect of Cohort
(damG.kr.l <- KRmodcomp(damG_noL, damG_fixed)) # test effect of Loc
(damG.kr.g <- KRmodcomp(damG_noG, damG_fixed)) # test effect of Group
(damG.kr.h <- KRmodcomp(damG_noH, damG_fixed)) # test effect of height




###################################################################################################################
######################################## Year 2 survival with herbivory ###########################################
###################################################################################################################



# outliers
adjbox(herbdata$LogDamPC, range=1.5) #anything > -0.5 is an outlier
herbdata3 <- herbdata[which(herbdata$LogDamPC < -0.5),] # N=1350 obs


# GLMM
herbsurv2 <- glmmTMB(survY2 ~ Group*Loc + Cohort + Group:Cohort + LogDamPC + Group:LogDamPC + Loc:LogDamPC + 
                       scale(H) + (1|Group:Geno) + (1|Loc:Block),
                     data=herbdata3,
                     family="binomial")

summary(herbsurv2)


### EMMs ###

gherbemms_survY2 <- ggemmeans(herbsurv2,~ LogDamPC|Group)




### LRTs ###


# (1|Group:Geno)
HS2ran1 <- glmmTMB(survY2 ~ Group*Loc + Cohort + Group:Cohort + LogDamPC + Group:LogDamPC + Loc:LogDamPC + scale(H) + (1|Loc:Block),
                   data=herbdata3,
                   family="binomial")
anova(HS2ran1,herbsurv2)


# (1|Loc:Block)
HS2ran2 <- glmmTMB(survY2 ~ Group*Loc + Cohort + Group:Cohort + LogDamPC + Group:LogDamPC + Loc:LogDamPC + scale(H) + (1|Group:Geno),
                   data=herbdata3,
                   family="binomial")
anova(HS2ran2,herbsurv2)


# Group:LogDamPC
HS2_GD <- glmmTMB(survY2 ~ Group*Loc + Cohort + Group:Cohort + LogDamPC + Loc:LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                  data=herbdata3,
                  family="binomial")
anova(HS2_GD,herbsurv2)


# Group:Cohort
HS2_GC <- glmmTMB(survY2 ~ Group*Loc + Cohort + LogDamPC + Group:LogDamPC + Loc:LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                  data=herbdata3,
                  family="binomial")
anova(HS2_GC,herbsurv2)


# Group:Loc
HS2_GL <- glmmTMB(survY2 ~ Group + Loc + Cohort + Group:Cohort + LogDamPC + Group:LogDamPC + Loc:LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                  data=herbdata3,
                  family="binomial")
anova(HS2_GL,herbsurv2)


# Loc:LogDamPC
HS2_LD <- glmmTMB(survY2 ~ Group*Loc + Cohort + Group:Cohort + LogDamPC + Group:LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                  data=herbdata3,
                  family="binomial")
anova(HS2_LD,herbsurv2)


## fixed effects only
HS2_fix <- glmmTMB(survY2 ~ Group + Loc + Cohort + LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                   data=herbdata3,
                   family="binomial")
## LogDamPC
HS2_D <- glmmTMB(survY2 ~ Group + Loc + Cohort + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                 data=herbdata3,
                 family="binomial")
anova(HS2_D,HS2_fix)


## Group
HS2_G <- glmmTMB(survY2 ~ Loc + Cohort + LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                 data=herbdata3,
                 family="binomial")
anova(HS2_G,HS2_fix)


## Loc
HS2_L <- glmmTMB(survY2 ~ Group + Cohort + LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                 data=herbdata3,
                 family="binomial")
anova(HS2_L,HS2_fix)


## Cohort
HS2_C <- glmmTMB(survY2 ~ Group + Loc + LogDamPC + scale(H) + (1|Group:Geno) + (1|Loc:Block),
                 data=herbdata3,
                 family="binomial")
anova(HS2_C,HS2_fix)


## Height
HS2_H <- glmmTMB(survY2 ~ Group + Loc + Cohort + LogDamPC + (1|Group:Geno) + (1|Loc:Block),
                 data=herbdata3,
                 family="binomial")
anova(HS2_H,HS2_fix)




###################################################################################################################
####################################   Impact of heterozygosity on fitness   ######################################
###################################################################################################################




#####################################   Heterozygous vs. homozygous sexuals   #####################################


# subset data to be sexuals-only 
sexdata <- data[which(data$BS=="S"),]
# be sure "het" and "survY1" are both factors

# make het a categorical variable; hetcat=1 is het=0; hetcat=2 is het>0
sexdata$hetcat <- cut(data$het, c(-1,0.07,1), labels=F)

sexdata$hetcat <- as.factor(sexdata$hetcat)



### GLMM ###
survFixcat <- glmmTMB(survY1 ~ hetcat + Loc + Cohort + scale(WinGH) + (1|hetcat:Geno) + (1|Loc:Block),
                      data=sexdata,
                      family="binomial")

summary(survFixcat)


### LRTs ###


# (1|hetcat:Geno)
ran1cat <- glmmTMB(survY1 ~ hetcat + Loc + Cohort + scale(WinGH) + (1|Loc:Block),
                   data=sexdata,
                   family="binomial")
anova(ran1cat,survFixcat)

# (1|Loc:Block)
ran2cat <- glmmTMB(survY1 ~ hetcat + Loc + Cohort + scale(WinGH) + (1|hetcat:Geno),
                   data=sexdata,
                   family="binomial")
anova(ran2cat,survFixcat)



# Cohort
### Model 1: hetcat + Loc + Cohort + WinGH (survFixcat) 
### Model 2: hetcat + Loc + WinGH (survD)

survDcat <- glmmTMB(survY1 ~ hetcat + Loc + scale(WinGH) + (1|hetcat:Geno) + (1|Loc:Block),
                    data=sexdata,
                    family="binomial")

anova(survDcat, survFixcat)


# Loc 
### Model 1: hetcat + Loc + Cohort + WinGH (survFixcat)
### Model 2: hetcat + Cohort + WinGH (survE)

survEcat <- glmmTMB(survY1 ~ hetcat + Cohort + scale(WinGH) + (1|hetcat:Geno) + (1|Loc:Block),
                    data=sexdata,
                    family="binomial")

anova(survEcat, survFixcat)


### Model 1: hetcat + Loc + Cohort + WinGH (survFixcat)
### Model 2: Loc + Cohort + WinGH (survF)

survFcat <- glmmTMB(survY1 ~ Loc + Cohort + scale(WinGH) + (1|hetcat:Geno) + (1|Loc:Block),
                    data=sexdata,
                    family="binomial")

anova(survFcat, survFixcat)


# WinGH 
### Model 1: hetcat + Loc + Cohort + WinGH (survFixcat)
### Model 2: hetcat + Loc + Cohort (survGcat)

survGcat <- glmmTMB(mortY1 ~ hetcat + Loc + Cohort + (1|hetcat:Geno) + (1|Loc:Block),
                    data=sexdata,
                    family="binomial")
anova(survGcat,survFixcat)



#######################################   Asexuals vs. heterozygous sexuals   ######################################



asexhetsex <- filter(data,hetcat==2) # this will get all heterozygous genotypes, and then can separate by RS
# N = 1459 A, 198 S; very imbalanced!



### GLMM

survFixSA <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                     data=asexhetsex,
                     family="binomial")

summary(survFixSA)


plot(ggeffect(survFixSA, terms="RS"))
ggeffect(survFixSA, terms=c("RS"))




### LRTs ###


# (1|RS:Geno)
ran1SA <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|Loc:Block),
                  data=asexhetsex,
                  family="binomial")
anova(ran1SA,survFixSA)

# (1|Loc:Block)
ran2SA <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|RS:Geno),
                  data=asexhetsex,
                  family="binomial")
anova(ran2SA,survFixSA)



# Cohort
### Model 1: RS + Loc + Cohort + WinGH (survFixSA) 
### Model 2: RS + Loc + WinGH (survDSA)

survDSA <- glmmTMB(survY1 ~ RS + Loc + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                   data=asexhetsex,
                   family="binomial")

anova(survDSA, survFixSA)


# Loc 
### Model 1: RS + Loc + Cohort + WinGH (survFixSA)
### Model 2: RS + Cohort + WinGH (survESA)

survESA <- glmmTMB(survY1 ~ RS + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                   data=asexhetsex,
                   family="binomial")

anova(survESA, survFixSA)

# RS  0.4948
### Model 1: RS + Loc + Cohort + WinGH (survFixSA)
### Model 2: Loc + Cohort + WinGH (survFSA)

survFSA <- glmmTMB(survY1 ~ Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                   data=asexhetsex,
                   family="binomial")

anova(survFSA, survFixSA)


# WinGH 
### Model 1: RS + Loc + Cohort + WinGH (survFixSA)
### Model 2: RS + Loc + Cohort (survGSA)

survGSA <- glmmTMB(mortY1 ~ RS + Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                   data=asexhetsex,
                   family="binomial")
anova(survGSA,survFixSA)




########################################   Asexuals vs. homozygous sexuals   #######################################



# a new variable of RS+hetcat
data$RShetcat <- paste0(data$RS, sexdata$hetcat)

# filter by certain values to obtain all asexuals and only homozygous sexuals; can then use RS term in model
AShom <- filter(data, RShetcat == "A2" | RShetcat == "S1")
# N=2864, 1418 A2 and 1446 S1 


### GLMM
survAhom <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                       data=AShom,
                       family="binomial")

summary(survAhom)

plot(ggeffect(survAhom, terms="RS"))
ggeffect(survAhom, terms=c("RS"))

### LRTs ###


# (1|RS:Geno)
ran1Ahom <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|Loc:Block),
                    data=asexhomsex,
                    family="binomial")
anova(ran1Ahom,survAhom)

# (1|Loc:Block)
ran2Ahom <- glmmTMB(survY1 ~ RS + Loc + Cohort + scale(WinGH) + (1|RS:Geno),
                    data=asexhomsex,
                    family="binomial")
anova(ran2Ahom,survAhom)



# Cohort
### Model 1: RS + Loc + Cohort + WinGH (survAhom) 
### Model 2: RS + Loc + WinGH (survDAhom)

survDAhom <- glmmTMB(survY1 ~ RS + Loc + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                     data=asexhomsex,
                     family="binomial")

anova(survDAhom, survAhom)


# Loc 
### Model 1: RS + Loc + Cohort + WinGH (survAhom) 
### Model 2: RS + Cohort + WinGH (survEAhom)

survEAhom <- glmmTMB(survY1 ~ RS + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                     data=asexhomsex,
                     family="binomial")

anova(survEAhom, survAhom)

# RS  
### Model 1: RS + Loc + Cohort + WinGH (survAhom)
### Model 2: Loc + Cohort + WinGH (survFAhom)

survFAhom <- glmmTMB(survY1 ~ Loc + Cohort + scale(WinGH) + (1|RS:Geno) + (1|Loc:Block),
                     data=asexhomsex,
                     family="binomial")

anova(survFAhom, survAhom)


# WinGH 
### Model 1: RS + Loc + Cohort + WinGH (survAhom)
### Model 2: RS + Loc + Cohort (survGAhom)

survGAhom <- glmmTMB(survY1 ~ RS + Loc + Cohort + (1|RS:Geno) + (1|Loc:Block),
                     data=asexhomsex,
                     family="binomial")

anova(survGAhom,survFixAhom)


### COMPARE ALL 3 P-VALUES FOR RS ###

# APO VS. HOM SEX: 7.257e-08
# HET SEX VS. HOM SEX: 0.007839
# APO VS. HET SEX: 0.4948

p.adjust(c(7.257e-08, 0.007839, 0.4948),method="holm")
# 2.1771e-07 0.015678 0.4948
