#!/usr/bin/env R
# Author: Barbara Diez Rodriguez (2021/07/01)
# Model testing


library(lme4)
library(piecewiseSEM)


#datasets needed ####
    #aphid_fitness
    #aphid_development
    #tannin_concentrations
    #tannin_induction
    #aphid_life_tables

#models ####

#effects of treatment on aphid performance


aphids <- read.csv("path/to/aphid_fitnes/dataset", header = TRUE) #aphid fitness data
aphids <- subset(aphids, aphids$value == "real" & aphids$aphid_treatment == "infested" ) #subset dataset and keep only real observations and infested values


aphids$N_treatment <- as.factor(aphids$N_treatment)

aov1<- aov(data=aphids, total_aphids ~ N_treatment) # ns. Remove treatment from models.
summary.aov(aov1)


#differences in aphid fitness between SwaspGTs

aphids$GT <- as.factor(aphids$GT)

aphids$env <- as.factor(aphids$env)

aov2 <- aov(data = aphids, total_aphids ~ GT + env + GT:env)
summary(aov2)


#differences in cummulative number of nymphs between SwaspGTs


x <- read.csv("path/to/aphid_development/dataset", header = TRUE) #growth curves dataset

x$GT <- as.factor(x$GT)
x$age <- as.factor(x$age)
summary(x)
rm.anova <- aov(data=x, cum_nymphs ~  GT*age + Error(plantID/age))
summary(rm.anova)



#correlation between total number of aphids and tannin induction

tannins <- read.csv2("path/to/tannin_concentrations/dataset", header = TRUE) #tannin induction data with aphid numbers

tannins$env <-as.factor(tannins$env)
tannins$GT <-as.factor(tannins$GT)
local <- subset(tannins, tannins$induction == "local")

glm1 <- glmer(data=local, family = poisson,  total_aphids ~ induced_tannins + (1|env/GT))
summary(glm1)
rsquared(glm1)




ind <- read.csv("path/to/tannin_induction/dataset", header = TRUE)

ind$GT <- as.factor(ind$GT)


glm2 <- glmer(data = ind, family = poisson, total_aphids ~ control_CTs + (1|GT))
summary(glm2)
rsquared(glm2)

glm3 <- glmer(data = ind, family = poisson, total_aphids ~ local_CTs + (1|GT))
summary(glm3)
rsquared(glm3)

glm4 <- glmer(data = ind, family = poisson, total_aphids ~  systemic_CTs + (1|GT))
summary(glm4)
rsquared(glm4)



#differences in life table parameters between SwaspGTs

life <- read.csv("path/to/Aphid_life_tables/dataset", check.names = F)
colnames(life)
life$GT <- as.factor(life$GT)
life1 <- subset(life, life$month == "november")
life2 <- subset(life, life$month == "march")

shapiro.test(life1$rm)
aov1 <- aov(data=life1, rm ~ GT) 
summary(aov1)

shapiro.test(life1$DT)
kruskal.test(data=life1, DT ~ GT)


shapiro.test(life2$rm)
kruskal.test(data=life2, rm ~ GT)

shapiro.test(life2$DT)
kruskal.test(data=life2, DT ~ GT)

colnames(life)
shapiro.test(life2$p)
kruskal.test(data = life2, p ~ GT)
