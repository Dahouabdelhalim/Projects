##########################################
###   SOCIAL HABITAT CHOICE ANALYSES   ###
##########################################

library(lme4)
library(MuMIn)
library(glmmTMB)
library(dplyr)
library(lmtest)

#first reading in the data and ensuring it's suitably restricted

data <- read.csv("Social_hab_choice_daily_scale.csv")

#removing individuals with few pre-experiment observations (100 obs is the point at which a correlation between number of obs and degree is no longer present)

data <- subset(data, no.obs.pre.exp > 99)

data <- subset(data, seen.after.end.feb == "Yes" | max.day.seen > 27)

data$species <- as.factor(as.character(data$species))

data$density.restriction <- as.factor(as.character(data$density.restriction))

# Probability of being seen on home feeder --------------------------------

on.home.feeder.dat <- subset(data, !is.na(on.home.feeder) & !is.na(on.other.feeder.in.pair))

on.home.feeder.dat$experimental.day.factor <- as.factor(as.character(on.home.feeder.dat$experimental.day))

prob.home.feeder.1.autocorr <- glmmTMB(on.home.feeder ~ species + experimental.day + deg.weighted + density.restriction + (1|feeder.site) + (1|tag) + ar1(experimental.day.factor +0|tag), data = on.home.feeder.dat, family = "binomial")

summary(prob.home.feeder.1.autocorr)

prob.home.feeder.2.autocorr <- glmmTMB(on.home.feeder ~ species + experimental.day + deg.weighted * density.restriction + (1|feeder.site) + (1|tag) + ar1(experimental.day.factor +0|tag), data = on.home.feeder.dat, family = "binomial")

summary(prob.home.feeder.2.autocorr)

#comparing models

AICc(prob.home.feeder.1.autocorr, prob.home.feeder.2.autocorr)

lrtest(prob.home.feeder.2.autocorr, prob.home.feeder.1.autocorr)

#no indication that birds with different social phenotypes respond differently to the 2 treatments

# Probability of being seen on the other feeder in the pair ---------------

on.other.feeder.dat <- subset(data, !is.na(on.other.feeder.in.pair) & !is.na(on.home.feeder))

on.other.feeder.dat$experimental.day.factor <- as.factor(as.character(on.other.feeder.dat$experimental.day))

prob.other.feeder.1.autocorr <- glmmTMB(on.other.feeder.in.pair ~ species + experimental.day + deg.weighted + density.restriction + (1|feeder.site) + (1|tag) + ar1(experimental.day.factor +0|tag), data = on.other.feeder.dat, family = "binomial")

summary(prob.other.feeder.1.autocorr)

prob.other.feeder.2.autocorr <- glmmTMB(on.other.feeder.in.pair ~ species + experimental.day + deg.weighted * density.restriction + (1|feeder.site) + (1|tag) + ar1(experimental.day.factor +0|tag), data = on.other.feeder.dat, family = "binomial")

summary(prob.other.feeder.2.autocorr)

lrtest(prob.other.feeder.2.autocorr, prob.other.feeder.1.autocorr)

# Looking at feeding behaviour at a site ----------------------------------

#even if birds stay at their feeding site, they may adjust their feeding behaviour in response to the social environment
#to look at this, I've used information on how long their feeding bouts were and how long intervals between them were for any given day

visit.length.1 <- glmer(mean.visit.length ~ species + scale(experimental.day) + scale(deg.weighted) + density.restriction + (1|feeder.site) + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log))

#feeder.site singular so remove and re-run

visit.length.1 <- glmer(mean.visit.length ~ species + scale(experimental.day) + scale(deg.weighted) + density.restriction + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log))

summary(visit.length.1)

visit.length.2 <- glmer(mean.visit.length ~ species + scale(experimental.day) + scale(deg.weighted) * density.restriction + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log), glmerControl(optimizer = "bobyqa"))

summary(visit.length.2)

lrtest(visit.length.2, visit.length.1)

#NOW LOOKING AT INTER VISIT INTERVALS

interval.length.1382.cutoff.1 <- glmer(mean.interval.length.1382.cutoff ~  species + scale(experimental.day) + scale(deg.weighted) + density.restriction + (1|feeder.site) + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log))

#feeder.site singular so remove and re-run

interval.length.1382.cutoff.1 <- glmer(mean.interval.length.1382.cutoff ~  species + scale(experimental.day) + scale(deg.weighted) + density.restriction + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log), glmerControl(optimizer = "bobyqa"))

summary(interval.length.1382.cutoff.1)

interval.length.1382.cutoff.2 <- glmer(mean.interval.length.1382.cutoff ~  species + scale(experimental.day) + scale(deg.weighted) * density.restriction + (1|tag), data = on.home.feeder.dat, family = Gamma(link = log), glmerControl(optimizer = "bobyqa"))

summary(interval.length.1382.cutoff.2)

lrtest(interval.length.1382.cutoff.2, interval.length.1382.cutoff.1)
