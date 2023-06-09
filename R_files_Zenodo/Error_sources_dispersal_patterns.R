### R script for the manuscript: Dispersal in a house sparrow metapopulation: an integrative 
### case study of genetic assignment calibrated with ecological data and pedigree information ###

# Dilan Saatoglu, dilansaatoglu@gmail.com
# 2021.07.05

# Packages used - install if it's uninstalled.
# install.packages(c("tidyverse","dplyr","glmmTMB","bbmle","DescTools"))

##############################################################################
### Pairwise number of misassignments vs pairwise genetic differentiation ###
#############################################################################

Pairwise.Fst = read.table("Pairwise_Fst.txt", header = T, sep = " ")
Pairwise.Fst$island.1 = as.factor(Pairwise.Fst$island.1)
Pairwise.Fst$island.2 = as.factor(Pairwise.Fst$island.2)
Pairwise.Fst$year = as.factor(Pairwise.Fst$year)

library(glmmTMB)

Pairwise.Fst$island.combination = paste0(Pairwise.Fst$island.1, Pairwise.Fst$island.2)
Pairwise.Fst$island.combination = as.factor(Pairwise.Fst$island.combination)

mod.pairwiseFst = glmmTMB(nr.mismatch ~ pw.Fst + (1 | island.combination), ziformula = ~1,
                   family=poisson,
                   data=Pairwise.Fst)

summary(mod.pairwiseFst)

######################################################################################
### Number of misassignments vs observed heterozygosity, sample size of baseline ###
### populations, and proportion of the population included in the baseline       ###
######################################################################################

# 1. Number of misassignments given away
# Data retrieval

WA.givenaway = read.table("Misassignments_givenaway.txt", header = T, sep = " ")
WA.givenaway$island = as.factor(WA.givenaway$island)
WA.givenaway$year = as.factor(WA.givenaway$year)

WA.received = read.table("Misassignments_received.txt", header = T, sep = " ")
WA.received$island = as.factor(WA.received$island)
WA.received$year = as.factor(WA.received$year)

  # i. Number of misassignments given away vs observed heterozygosity
library(glmmTMB)

mod.Ho_givenaway = glmmTMB(misassignment ~ Ho, ziformula = ~1, 
                       family=poisson,
                       data=WA.givenaway)

summary(mod.Ho_givenaway)

  # ii. Number of misassignments given away vs log(e)(sample size)
library(glmmTMB)

WA.givenaway$log.samplesize = log(WA.givenaway$sample.size)

mod.samplesize_givenaway = glmmTMB(misassignment ~ log.samplesize, ziformula = ~1, 
                           family=poisson,
                           data=WA.givenaway)

summary(mod.samplesize_givenaway)

  # iii. Number of misassignments given away vs proportion of population sampled
library(glmmTMB)

WA.givenaway$proportion.sampled = WA.givenaway$sample.size/WA.givenaway$Population.size

mod.prosampled_givenaway = glmmTMB(misassignment ~ proportion.sampled, ziformula = ~1, 
                                   family=poisson,
                                   data=WA.givenaway)

summary(mod.prosampled_givenaway)

# 2. Number of misassignments received
  # i. Number of misassignments received vs observed heterozygosity
library(glmmTMB)

mod.Ho_received = glmmTMB(misassignment ~ Ho, ziformula = ~1, 
                           family=poisson,
                           data=WA.received)

summary(mod.Ho_received)

  # ii. Number of misassignments received vs log(e)(sample size)
library(glmmTMB)

WA.received$log.samplesize = log(WA.received$sample.size)

mod.samplesize_received = glmmTMB(misassignment ~ log.samplesize, ziformula = ~1, 
                                   family=poisson,
                                   data=WA.received)

summary(mod.samplesize_received)

  # iii. Number of misassignments received vs proportion of population sampled
library(glmmTMB)

WA.received$proportion.sampled = WA.received$sample.size/WA.received$Population.size

mod.prosampled_received = glmmTMB(misassignment ~ proportion.sampled, ziformula = ~1, 
                                   family=poisson,
                                   data=WA.received)

summary(mod.prosampled_received)

#######################################################################
### General dispersal patterns in the house sparrow metapopulation ###
#######################################################################

Dispersal.data = read.table("Dispersaldata.txt", header = T, sep = " ")
Dispersal.data[1:6] = lapply(Dispersal.data, as.factor)

##############################################
### Spatio-temporal variation in dispersal ###
##############################################

# Data management

library(tidyverse)
library(dplyr)

Genpattern.data = 
  Dispersal.data[(Dispersal.data$natal.island %in% c("20","22","23","24","26","27","28","38","77","88")),]

Genpattern.data = Genpattern.data %>%
  group_by(hatch.year,natal.island,sex,dispersal.status) %>%
  tally() 

Genpattern.data = as.data.frame(Genpattern.data)

Genpattern.data$sex = as.numeric(Genpattern.data$sex)
Genpattern.data$dispersal.status = as.numeric(Genpattern.data$dispersal.status)

Genpattern.data[3:4] = lapply(names(Genpattern.data)[3:4], 
                              function(x) paste(x, Genpattern.data[, x], sep = "."))
Genpattern.data = Genpattern.data %>%
  unite(key, sex, dispersal.status, sep = ".") %>%
  spread(key, n, fill = 0)

colnames(Genpattern.data)[3:6] = c("Number.maleresidents","Number.maledispersers",
                                   "Number.femaleresidents","Number.femaledispersers")

Genpattern.data$Number.recruits = rowSums(Genpattern.data[3:6])
Genpattern.data$Number.residents = rowSums(Genpattern.data[,c(3,5)])
Genpattern.data$Number.dispersers = rowSums(Genpattern.data[,c(4,6)])

# Models for spatio-temporal variation in dispersal

library(glmmTMB)
library(bbmle)

mod.tempvar.1 = glmmTMB(cbind(Number.dispersers,Number.residents) ~ 1, family=binomial, 
               data=Genpattern.data)
summary(mod.tempvar.1)

mod.tempvar.2 = glmmTMB(cbind(Number.dispersers,Number.residents) ~ 1 + hatch.year, family=binomial,
               data=Genpattern.data)
summary(mod.tempvar.2)

mod.tempvar.3 = glmmTMB(cbind(Number.dispersers,Number.residents) ~ 1 + natal.island, family=binomial,
               data=Genpattern.data)
summary(mod.tempvar.3)

mod.tempvar.4 = glmmTMB(cbind(Number.dispersers,Number.residents) ~ 1 + 
                        hatch.year + natal.island, family=binomial,
                        data=Genpattern.data)
summary(mod.tempvar.4)

anova(mod.tempvar.1,mod.tempvar.2,mod.tempvar.3,mod.tempvar.4, test = "LRT")
AICtab(mod.tempvar.1,mod.tempvar.2,mod.tempvar.3,mod.tempvar.4)

  # i.Temporal variation in number of dispersers
anova(mod.tempvar.1,mod.tempvar.2, test="LRT") ### LRT for the year
# x2(14) = 34.666, p=0.001647 ***
# χ2(df)= chisq, p 

  # ii. Spatial variation in number of dispersers
anova(mod.tempvar.1,mod.tempvar.3, test="LRT") ### LRT for the natal island
# x2(9) = 510.25, p=2.2e-16 ***
# χ2(df)= chisq, p


######################################################
### Variation between-within habitats in dispersal ###
#####################################################

# Data management
# Because the sampling period for non-farm islands did not include years before 2004, 
# only data from years 2004–2012 was used.
# Inviduals hatched and ended up on 8 SNP islands were used.
# (Island codes: 20,22,23,24,26,27,28,38)

Habitat.data = Dispersal.data[!(Dispersal.data$hatch.year %in% 
                                     c("1998","1999","2000","2001","2002","2003")) &
                                   Dispersaldata$natal.island %in% 
                                   c("20","22","23","24","26","27","28","38") &
                                   Dispersaldata$adult.island %in% 
                                   c("20","22","23","24","26","27","28","38"),]
                                 
Habitat.data = Habitat.data[,c(1,4,2,3)]

for (i in 1:nrow(Habitat.data)){
  if(Habitat.data[i,3] %in% c("23","24","22")) {
    Habitat.data[i,5] = "Non-farm"
  } else {
    Habitat.data[i,5] = "Farm"
  }
}

for (i in 1:nrow(Habitat.data)){
  if(Habitat.data[i,4] %in% c("23","24","22")) {
    Habitat.data[i,6] = "Non-farm"
  } else {
    Habitat.data[i,6] = "Farm"
  }
}

for (i in 1:nrow(Habitat.data)){
  if(Habitat.data[i,2] == 1 & (Habitat.data[i,5] == Habitat.data[i,6])) {
    Habitat.data[i,7] = "1"
  } else {
    Habitat.data[i,7] = "0"
  }
}

for (i in 1:nrow(Habitat.data)){
  if(Habitat.data[i,2] == 1 & (Habitat.data[i,5] != Habitat.data[i,6])) {
    Habitat.data[i,8] = "1"
  } else {
    Habitat.data[i,8] = "0"
  }
}

colnames(Habitat.data)[5:8] = c("natal.habitat","adult.habitat",
                                 "withinhab.disperser","betweenhab.disperser")
Habitat.data[1:8] = lapply(Habitat.data, factor)

  # i. Variation in dispersal within habitat types

library(glmmTMB)

habitat.within = glmmTMB(withinhab.disperser ~ natal.habitat, 
                       data = Habitat.data, family = binomial)
summary(habitat.within)

  # ii. Variation in dispersal between habitat types

habitat.between = glmmTMB(betweenhab.disperser ~ natal.habitat, 
               data = Habitat.data, family = binomial)
summary(habitat.between)

#######################################################################
### Variation in dispersal between the sexes - Sex biased dispersal ###
#######################################################################

library(glmmTMB)

sex.mod.0 = glmmTMB(dispersal.status ~ 1 + (1| natal.island) + (1 | hatch.year),
              data = Dispersal.data, family = binomial)
summary(sex.mod.0)

sex.mod.1 = glmmTMB(dispersal.status ~ sex + (1| natal.island) + (1 | hatch.year),
                    data = Dispersal.data, family = binomial)
summary(sex.mod.1)

anova(sex.mod.0,sex.mod.1, test = "LRT")


###############################################################################################
### Sex ratio difference between the datasets (ecological & integrative case study dataset) ###
###############################################################################################

# Data management
# Numbers from the Dispersal.data dataframe
# Proportion for female dispersers for integrative case data = 356/607*100 = 58.6
# Proportion for female dispersers for only ecological data = 160/291*100 = 55
# Proportion for male dispersers for integrative case data = 251/607*100 = 41.4
# Proportion for male dispersers for only ecological data = 131/291*100 = 45

Sexdifferproportions.data = data.frame( Status = c("Female.dispersers","Male.dispersers"),
                         Integrative = c(58.6,41.4),
                         Ecological = c(55,45))
library(DescTools)

GTest(Sexdifferproportions.data[1,2:3])

# G = 0.114, df = 1, P > 0.05