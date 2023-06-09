# bamm post-analysis

setwd("~/Documents/HYDIV")
rm(list = ls())

library(tidyverse)
library(ape)
library(BAMMtools)
library(coda)
library(dplyr)


# read in the family-level phylogeny
phy <- read.tree("data/full_tree.tre")
a
# get total event  data (take a long time)
edata_25 <- getEventData(phy, eventdata = "data/event_data_prior25.txt", burnin = 0.2)
edata_50 <- getEventData(phy, eventdata = "data/event_data_prior50.txt", burnin = 0.2)
edata_100 <- getEventData(phy, eventdata = "data/event_data_prior100.txt", burnin = 0.2)
edata_200 <- getEventData(phy, eventdata = "data/event_data_prior200.txt", burnin = 0.2)

# get out files
mcmcout_25 <- read.csv("data/mcmc_out_prior25.txt", header = T)
mcmcout_50 <- read.csv("data/mcmc_out_prior50.txt", header = T)
mcmcout_100 <- read.csv("data/mcmc_out_prior100.txt", header = T)
mcmcout_250 <- read.csv("data/mcmc_out_prior200.txt", header = T)

# discard 10% of samples as burnin
burnstart_25 <- floor(0.2 * nrow(mcmcout_25))
postburn_25 <- mcmcout[burnstart_25:nrow(mcmcout_25), ]
plot(postburn_25$logLik ~ postburn_25$generation)

burnstart_50 <- floor(0.2 * nrow(mcmcout_50))
postburn_50 <- mcmcout[burnstart_50:nrow(mcmcout_50), ]
plot(postburn_50$logLik ~ postburn_50$generation)

burnstart_100 <- floor(0.2 * nrow(mcmcout_100))
postburn_100 <- mcmcout[burnstart_100:nrow(mcmcout_100), ]
plot(postburn_100$logLik ~ postburn_100$generation)

burnstart_200 <- floor(0.2 * nrow(mcmcout_200))
postburn_200 <- mcmcout[burnstart_200:nrow(mcmcout_200), ]
plot(postburn_200$logLik ~ postburn_200$generation)

# want these to be AT LEAST 200
effectiveSize(postburn_25$N_shifts)
effectiveSize(postburn_25$logLik)

effectiveSize(postburn_50$N_shifts)
effectiveSize(postburn_50$logLik)

effectiveSize(psotburn_100$N_shifts)
effectiveSize(postburn_100$logLik)

effectiveSize(postburn_200$logLik)
effectiveSize(postburn_200$logLik)

# create objects for average tip diversification rates 
bamm_data_25 <- data.frame(Family = phy$tip.label, speciation = edata_25$meanTipLambda, extinction = edata_25$meanTipMu, net = edata_25$meanTipLambda - edata_25$meanTipMu)
bamm_data_50 <- data.frame(Family = phy$tip.label, speciation = edata_50$meanTipLambda, extinction = edata_50$meanTipMu, net = edata_50$meanTipLambda - edata_50$meanTipMu)
bamm_data_100 <- data.frame(Family = phy$tip.label, speciation = edata_100$meanTipLambda, extinction = edata_100$meanTipMu, net = edata_100$meanTipLambda - edata_100$meanTipMu)
bamm_data_200 <- data.frame(Family = phy$tip.label, speciation = edata_200$meanTipLambda, extinction = edata_200$meanTipMu, net = edata_200$meanTipLambda - edata_200$meanTipMu)

