#R code for the analyses on group social network metrics in the manuscript "Resource distribution alters individual social network positions but not group network structures in experimental populations of forked fungus beetles"
#This script measures the effect of resource distribution on group social network metrics (average shortest path length, density, global clustering coefficient). This script uses paired t-tests to answer how resource distribution affects population-level social network metrics. The population-level social network metrics are built from networks that include all interactions and both sexes.
#Data are available at **add dryad link**

#set working directory to where you saved the data file
setwd("")

#load libraries
library(tidyverse)
library(car)

#Change R's defult to contrasts 
options(contrasts=c("contr.sum", "contr.poly")) 

#upload data
ObsPopSN <- read.csv("groupSN.csv")

#restructure group SN dataframe for paired t-tests####
ObsPopSN %>%
  filter(Treatment == "Clumped") %>%
  dplyr::select(avg.shortest.path, global.cc, densityw, Condo) %>%
  dplyr::rename(avg.shortest.path.clump=avg.shortest.path, global.cc.clump=global.cc, densityw.clump=densityw) -> PopSNClump

ObsPopSN %>%
  filter(Treatment == "Dispersed") %>%
  dplyr::select(avg.shortest.path, global.cc, densityw, Condo) %>%
  dplyr::rename( avg.shortest.path.disperse=avg.shortest.path, global.cc.disperse=global.cc, densityw.disperse=densityw) -> PopSNDisperse

PopSNLong <- merge(PopSNClump, PopSNDisperse, by="Condo", all.x=TRUE, all.y=FALSE)

#Do group social network metrics differ across resource distribution treatments?####
##weighted network density####

###paired t-test
t.test(PopSNLong$densityw.clump, PopSNLong$densityw.disperse, paired=TRUE, conf.level=0.95)

###check assumptions
densityw.difference = PopSNLong$densityw.clump - PopSNLong$densityw.disperse
hist(densityw.difference)
qqp(densityw.difference)

##average shortest path length####

###paired t-test
t.test(PopSNLong$avg.shortest.path.clump, PopSNLong$avg.shortest.path.disperse, paired=TRUE, conf.level=0.95)

###check assumptions
avg.shortest.path.difference = PopSNLong$avg.shortest.path.clump - PopSNLong$avg.shortest.path.disperse
hist(avg.shortest.path.difference)
qqp(avg.shortest.path.difference)

##global clustering coefficient####

###paired t-test
t.test(PopSNLong$global.cc.clump, PopSNLong$global.cc.disperse, paired=TRUE, conf.level=0.95)

###check assumptions
global.cc.difference = PopSNLong$global.cc.clump - PopSNLong$global.cc.disperse
hist(global.cc.difference)
qqp(global.cc.difference)