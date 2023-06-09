mdata <- read.csv("~/Documents/Hybrid migration/Hybrid migration data/hybmigr-submit.csv", header = TRUE, sep = ",", row.names = NULL)
mdata$length <- as.numeric(mdata$length)
mdata$pred_by_cormorant <- factor(mdata$pred_by_cormorant)
mdata$survived_to_spring <- factor(mdata$survived_to_spring)
mdata$migr_trips <- as.numeric(mdata$migr_trips)
mdata$obs_period_days <- as.numeric(mdata$obs_period_days)
mdata$first_migr_numeric <- as.numeric(mdata$first_migr_numeric)
mdata$first_migr <- as.Date(mdata$first_migr,"%m/%d/%Y")
mdata$first_migr_noNA <- as.numeric(mdata$first_migr_noNA)

library(car)
library(multcomp)
library(MuMIn)

#ANALYSIS 1: Migration trips between species and hybrids (GLM)

#Limit to individuals that showed activity after the winter
springsurv <- subset(mdata, mdata$survived_to_spring =="1")

# Poisson GLM for species differences in migration trips
trips.model <- glm(migr_trips~Species, data=springsurv, family=poisson)
summary(trips.model)
anova(trips.model, test = "Chisq")

# Tukey Post Hoc
tuk.species <- glht(trips.model, linfct = mcp(Species = "Tukey"))
summary(tuk.species)

#ANALYSIS 2: Date of lake departure between species and hybrids

kruskal.test(first_migr_numeric ~ Species, data = mdata)
pairwise.wilcox.test(mdata$first_migr_numeric, mdata$Species, p.adjust.method = "BH")

#ANALYSIS 3: Model predicting cormorant predation

# Limit used variables to avoid NA's in the dataset
mdata.log <- subset(mdata, select = c(Species,length, pred_by_cormorant, obs_period_days, migr_trips, migr_frequency, first_migr_noNA))

#Full model
full.logmodel <- glm(pred_by_cormorant~Species+obs_period_days+Species*obs_period_days
                     +length+migr_frequency+first_migr_noNA, data=mdata.log, family=binomial)

summary(full.logmodel)
Anova(full.logmodel)

options(na.action = "na.fail")  
output <- dredge(full.logmodel, rank=AIC)
output
options(na.action = "na.omit") 

#Final model
log.model <- glm(pred_by_cormorant~Species+migr_frequency+obs_period_days, data=mdata.log, family=binomial)
summary(log.model)
Anova(log.model)
