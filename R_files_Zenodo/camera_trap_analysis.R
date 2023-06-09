
#### Small mammals reduce distance-dependence and increase seed predation risk in tropical rainforest fragments ####
#### Biotropica
#### Authors: Aparna KRISHNAN*, Anand M OSURI, Meghna KRISHNADAS

#### Description: This file contains the analysis for section 2.3 (Mammal visitation), Table 2, and Table S3 in the manuscript and the supplementary material. The code summarises mammal seed predator visitation rates and proportion of seed removal vs. predation by mammal seed predators.

library(dplyr)
library(ggplot2)

#setwd("~/Documents/NCF/Seed Predation Manuscript/Data Submission/Camera Trap")

########### MAIN DOCUMENT ##########

#### TABLE 2 ####
cmtp <- read.csv(file = "camera_trap_data.csv", header = T, stringsAsFactors = F)

# Vists  of each predator species group to camera traps
  # 1 if the animal was passing by or predating
  # 0 if the animal did not visit the camera trap at all
cmtp$temp <- ifelse(cmtp$Activity == "None", 0,1)

cmtp_sum <- cmtp %>% group_by(Tree_Species, Landscape_Location, Plot_No, Plot_Location, Habitat_Status, Species_Group) %>% summarise(
  visits = sum(temp),
  days = unique(Camera_Days))

# Visitation Rate of each predator species group to each camera trap
cmtp_sum$vrate <- 100*cmtp_sum$visits / cmtp_sum$days

# Table 2
# Average Visitation Rate for each species group
cmtp_sum %>% group_by(Species_Group, Habitat_Status) %>% summarise(
  av_vrate = mean(vrate),
  se_vrate = sd(vrate)/sqrt(length(vrate)),
  ss = length(vrate))

#### SEED REMOVAL for Section 2.3 ####

# Removed by rodents

#total seed fates observed
tot.obs <- sum(cmtp$No_Removed, na.rm = T) + sum(cmtp$No_Predated, na.rm = T )
# 181

# % seeds removed by Rodents
rod.remv <- sum(cmtp$No_Removed[which(cmtp$Species_Group == "Rodent")], na.rm = T)
100*rod.remv/tot.obs
# 6.629834

# % seeds removed by malabar spiny dormouse
spn.remv <- sum(cmtp$No_Removed[which(cmtp$Visitor_Species == "Spiny Tailed Dormouse")], na.rm = T)
100*spn.remv/tot.obs
# 2.209945

########### SUPPORTING INFORMATION ##########

#### TABLE S2 ####
cmpt_ss <- cmtp %>% group_by(Tree_Species, Habitat_Status, Plot_No, Plot_Location, Landscape_Location) %>% summarise(days = unique(Camera_Days))

# Number of camera traps and trap nights in each fragment and contiguous forest location
SupTable2 <- cmpt_ss %>% group_by(Tree_Species, Habitat_Status, Landscape_Location) %>% summarise(totdays = sum(days), ss = length(days))

#### TABLE S3 ####
cmtp_sum <- cmtp %>% group_by(Tree_Species, Landscape_Location, Plot_No, Plot_Location, Habitat_Status, Visitor_Species) %>% summarise(
  visits = sum(temp),
  days = unique(Camera_Days),
  pred = sum(No_Predated, na.rm = T),
  rmvd = sum(No_Removed, na.rm = T))

cmtp_sum$vrate <- 100*cmtp_sum$visits / cmtp_sum$days

# visitation rate
SupTable3_vrate <- cmtp_sum %>% group_by(Visitor_Species, Habitat_Status) %>% summarise(
  av_vrate = round(mean(vrate), 2),
  se_vrate = round(sd(vrate)/sqrt(length(vrate)),2),
  ss = length(vrate))

# Proportion removed and predated
SupTable3_rmvl <- cmtp_sum %>% group_by(Visitor_Species) %>% summarise(
  pct_pred = round(100*sum(pred)/tot.obs,2),
  pct_rmvd= round(100*sum(rmvd)/tot.obs,2),
  ss = length(vrate))