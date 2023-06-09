#01 Apr 2022
#V4.1.1
#---------


# Cleaned up for publication 


### Paired with

#Tank means
# "Larval_end-point_tank_totals.csv"

#Mass data
# "Larval_end-point_size.csv"



#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting
library(multcomp)   # needed for contrasts
library(car)        # needed for type II ANOVA command 'Anova()'


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this R file is saved in


Larv.tanks <- read.csv("Larval_end-point_tank_totals.csv", header = TRUE)

Density.initial <- Larv.tanks$Trial.Stocking/50 # no. larvae per liter (tanks were 50 L)
Treatment <- Larv.tanks$Treatment




##################### GLM for Survival ###########################


Survivorship <- Larv.tanks$Trial.Survival

## First, checking for effect of interaction
survivorship.glm.int <- glm(Survivorship ~ Treatment * Density.initial, family = quasibinomial(link = "logit"))
summary(survivorship.glm.int)

surv1.int <- Anova(survivorship.glm.int)
surv1.int

## No effect of interaction - remove



survivorship.glm <- glm(Survivorship ~ Treatment + Density.initial, family = quasibinomial(link = "logit"))
summary(survivorship.glm)


surv1 <- Anova(survivorship.glm)
surv1


#Collecting P-values for treatment and density separately to later use Holm's correction for all analyses in this file
p.all.trt <- vector(mode = "numeric", length = 4)
p.all.den <- vector(mode = "numeric", length = 4)
p.all.trt[1] <- surv1$`Pr(>Chisq)`[1]
p.all.den[1] <- surv1$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(survivorship.glm) 
qqnorm(resid(survivorship.glm, type = "pearson"))
qqline(resid(survivorship.glm, type = "pearson"))





##################### GLM for Swim bladder inflation ###########################


Inflated <- Larv.tanks$SB.inflation


## First, checking for effect of interaction

inflated.glm.int <- glm(Inflated ~ Treatment * Density.initial, family = quasibinomial(link = "logit"))
summary(inflated.glm.int)

inf1.int <- Anova(inflated.glm.int)
inf1.int

## No effect of interaction - remove



inflated.glm <- glm(Inflated ~ Treatment + Density.initial, family = quasibinomial(link = "logit"))
summary(inflated.glm)

inf1 <- Anova(inflated.glm)
inf1

p.all.trt[2] <- inf1$`Pr(>Chisq)`[1]
p.all.den[2] <- inf1$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(inflated.glm) 
qqnorm(resid(inflated.glm, type = "pearson"))
qqline(resid(inflated.glm, type = "pearson"))




















##################################################################################################################
#......................... Looking at mean tank sizes, SGR, and densities for larval period ..........................


Larv <- read.csv("Larval_end-point_size.csv", header = TRUE)

Larv.2 <- Larv %>%
  mutate(Treatment = Tank)

Larv.2$Treatment <- as.character(Larv.2$Treatment)
Larv.2$Treatment[grepl("1", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("2", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("3", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("4", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("5", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("6", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("7", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("8", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("9", Larv.2$Treatment)] <- "C"


#Calculating means and standard deviations for mass and length of inflated fish
Larv.size.inf.means <- Larv.2 %>%
  filter(SB == "1") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.inf.weight = mean(Weight..mg.), Avg.larv.inf.TL = mean(length..mm.),
            SD.larv.inf.weight = sd(Weight..mg.), SD.larv.inf.TL = sd(length..mm.))

#Calculating means and standard deviations for mass and length of uninflated fish
Larv.size.uninf.means <- Larv.2 %>%
  filter(SB == "0") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.uninf.weight = mean(Weight..mg.), Avg.larv.uninf.TL = mean(length..mm.),
            SD.larv.uninf.weight = sd(Weight..mg.), SD.larv.uninf.TL = sd(length..mm.))


#Calculating out beginning and end values for error bars that would be +/- 1 standard deviation
Larv.size.means.plus.dens <- full_join(Larv.tanks, Larv.size.inf.means, by = "Tank", "Treatment")
Larv.size.means.plus.dens2 <- full_join(Larv.size.means.plus.dens, Larv.size.uninf.means, by = "Tank") %>%
  mutate(weight.larv.inf.min = Avg.larv.inf.weight - SD.larv.inf.weight) %>%
  mutate(weight.larv.inf.max = Avg.larv.inf.weight + SD.larv.inf.weight) %>%
  mutate(weight.larv.uninf.min = Avg.larv.uninf.weight - SD.larv.uninf.weight) %>%
  mutate(weight.larv.uninf.max = Avg.larv.uninf.weight + SD.larv.uninf.weight) %>%
  mutate(length.larv.inf.min = Avg.larv.inf.TL - SD.larv.inf.TL) %>%
  mutate(length.larv.inf.max = Avg.larv.inf.TL + SD.larv.inf.TL) %>%
  mutate(length.larv.uninf.min = Avg.larv.uninf.TL - SD.larv.uninf.TL) %>%
  mutate(length.larv.uninf.max = Avg.larv.uninf.TL + SD.larv.uninf.TL) %>%
  mutate(Density.larv = Trial.Stocking/50) %>%
  mutate(SGR.P1 = (((Avg.larv.inf.weight/2.55)^(1/10)) - 1))




#------------------- Only using inflated fish ----------------------

## First, checking for effect of interaction

Inf.mass.glm.int <- glm(Avg.larv.inf.weight ~ Treatment * Density.larv, data = Larv.size.means.plus.dens2, family = gaussian)
summary(Inf.mass.glm.int)

mass1.int <- Anova(Inf.mass.glm.int)
mass1.int

## No effect of interaction - remove




Inf.mass.glm <- glm(Avg.larv.inf.weight ~ Treatment + Density.larv, data = Larv.size.means.plus.dens2, family = gaussian)
summary(Inf.mass.glm)

mass1 <- Anova(Inf.mass.glm)
mass1

p.all.trt[3] <- mass1$`Pr(>Chisq)`[1]
p.all.den[3] <- mass1$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(Inf.mass.glm) 
qqnorm(resid(Inf.mass.glm, type = "pearson"))
qqline(resid(Inf.mass.glm, type = "pearson"))










#----------- Contrasts: 
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept

Inf.mass.glm.cont <- glm(Avg.larv.inf.weight ~ Treatment + Density.larv - 1, data = Larv.size.means.plus.dens2, family = gaussian)
summary(Inf.mass.glm.cont)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(Inf.mass.glm.cont, linfct = Cont1)
summary(cont.test1)

Cont2 <- matrix(c(0, -1, 1, 0), 1) #PUFA vs. PUFA + a-t
cont.test2 <- glht(Inf.mass.glm.cont, linfct = Cont2)
summary(cont.test2)









############# SGR ################

#------------------- Only using inflated fish ----------------------

## First, checking for effect of interaction

Inf.SGR.glm.int <- glm(SGR.P1 ~ Treatment * Density.larv, data = Larv.size.means.plus.dens2, family = quasibinomial(link = "logit"))
summary(Inf.SGR.glm.int)

SGR1.int <- Anova(Inf.SGR.glm.int)
SGR1.int

## No effect of interaction - remove




Inf.SGR.glm <- glm(SGR.P1 ~ Treatment + Density.larv, data = Larv.size.means.plus.dens2, family = quasibinomial(link = "logit"))
summary(Inf.SGR.glm)

SGR1 <- Anova(Inf.SGR.glm)
SGR1
          
p.all.trt[4] <- SGR1$`Pr(>Chisq)`[1]
p.all.den[4] <- SGR1$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(Inf.SGR.glm) 
qqnorm(resid(Inf.SGR.glm, type = "pearson"))
qqline(resid(Inf.SGR.glm, type = "pearson"))










#----------- Contrasts: 
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept

Inf.SGR.glm.cont <- glm(SGR.P1 ~ Treatment + Density.larv - 1, data = Larv.size.means.plus.dens2, family = quasibinomial(link = "logit"))
summary(Inf.SGR.glm.cont)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(Inf.SGR.glm.cont, linfct = Cont1)
summary(cont.test1)

Cont2 <- matrix(c(0, -1, 1, 0), 1) #PUFA vs. PUFA + a-t
cont.test2 <- glht(Inf.SGR.glm.cont, linfct = Cont2)
summary(cont.test2)































#%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holm's Correction for multiple comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.new.trt <- p.adjust(p.all.trt, method = "holm", n = length(p.all.trt))
p.new.trt
#Order is: survival, swim bladder inflation, mass, SGR

p.new.den <- p.adjust(p.all.den, method = "holm", n = length(p.all.den))
p.new.den
#Order is: survival, swim bladder inflation, mass, SGR


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















#$$$$$$$$$$$$$$$ Quick check to see if mass results are same if inflated and uninflated fish not separated




#Calculating means and standard deviations for mass and length of inflated fish
Larv.size.means <- Larv.2 %>%
  group_by(Tank) %>%
  summarize(Avg.larv.weight = mean(Weight..mg.), Avg.larv.TL = mean(length..mm.),
            SD.larv.weight = sd(Weight..mg.), SD.larv.TL = sd(length..mm.))

#Calculating out beginning and end values for error bars that would be +/- 1 standard deviation
Larv.size.means.plus.dens3 <- full_join(Larv.tanks, Larv.size.means, by = "Tank", "Treatment") %>%
  mutate(weight.larv.min = Avg.larv.weight - SD.larv.weight) %>%
  mutate(weight.larv.max = Avg.larv.weight + SD.larv.weight) %>%
  mutate(length.larv.min = Avg.larv.TL - SD.larv.TL) %>%
  mutate(length.larv.max = Avg.larv.TL + SD.larv.TL) %>%
  mutate(Density.larv = Trial.Stocking/50) %>%
  mutate(SGR.P1 = (((Avg.larv.weight/2.55)^(1/10)) - 1))


#---- Mass ----
mass.glm <- glm(Avg.larv.weight ~ Treatment + Density.larv, data = Larv.size.means.plus.dens3, family = gaussian)
summary(mass.glm)

mass2 <- Anova(mass.glm)
mass2





#Contrasts: 

mass.glm.cont <- glm(Avg.larv.weight ~ Treatment + Density.larv - 1, data = Larv.size.means.plus.dens3, family = gaussian)
summary(mass.glm.cont)

Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(mass.glm.cont, linfct = Cont1)
summary(cont.test1)
Cont2 <- matrix(c(0, -1, 1, 0), 1) #PUFA vs. PUFA + a-t
cont.test2 <- glht(mass.glm.cont, linfct = Cont2)
summary(cont.test2)








qqp(resid(mass.glm))





#----- SGR -----

SGR.glm <- glm(SGR.P1 ~ Treatment + Density.larv, data = Larv.size.means.plus.dens3, family = quasibinomial(link = "logit"))
summary(SGR.glm)

SGR2 <- Anova(SGR.glm)
SGR2



#Contrasts:


SGR.glm.cont <- glm(SGR.P1 ~ Treatment + Density.larv - 1, data = Larv.size.means.plus.dens3, family = quasibinomial(link = "logit"))
summary(SGR.glm.cont)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(SGR.glm.cont, linfct = Cont1)
summary(cont.test1)

Cont2 <- matrix(c(0, -1, 1, 0), 1) #PUFA vs. PUFA + a-t
cont.test2 <- glht(SGR.glm.cont, linfct = Cont2)
summary(cont.test2)





qqp(resid(SGR.glm))







