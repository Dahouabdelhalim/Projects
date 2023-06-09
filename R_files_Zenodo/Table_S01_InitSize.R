# NSF FX data presentation

# Table S01
# Initial sizes of trees

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

##############################################
################## New York ##################
##############################################

dat <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]

trees_to_use <- dat$Use_growth_01==1

# Initial Robinia diameter and height
mean(dat[dat$Species=="ROPS" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="ROPS" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Betula diameter and height
mean(dat[dat$Species=="BENI" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="BENI" & trees_to_use,]$Height_cm_01,na.rm=TRUE)

##############################################
################### Oregon ###################
##############################################

dat <- read.csv("OR_FX_Size_FoliarCNIsotope_Data.csv")[1:64,1:122]

trees_to_use <- dat$Use_growth_01==1

# Initial Alnus diameter and height
mean(dat[dat$Species=="ALRU" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="ALRU" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Pseudotsuga diameter and height
mean(dat[dat$Species=="PSME" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="PSME" & trees_to_use,]$Height_cm_01,na.rm=TRUE)

###############################################
################### Waiakea ###################
###############################################

dat <- read.csv("HI_W_FX_Size_FoliarCNIsotope_Data.csv")[1:108,1:111]

trees_to_use <- dat$Use_growth_01==1

# Initial Gliricidia diameter and height
mean(dat[dat$Species=="GLSE" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="GLSE" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Casuarina diameter and height
mean(dat[dat$Species=="CAEQ" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="CAEQ" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Psidium diameter and height
mean(dat[dat$Species=="PSCA" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="PSCA" & trees_to_use,]$Height_cm_01,na.rm=TRUE)

###############################################
################### Volcano ###################
###############################################

dat <- read.csv("HI_V_FX_Size_FoliarCNIsotope_Data.csv")[1:96,1:111]

trees_to_use <- dat$Use_growth_01==1

# Initial Acacia diameter and height
mean(dat[dat$Species=="ACKO" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="ACKO" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Morella diameter and height
mean(dat[dat$Species=="MOFA" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="MOFA" & trees_to_use,]$Height_cm_01,na.rm=TRUE)
# Initial Dodonaea diameter and height
mean(dat[dat$Species=="DOVI" & trees_to_use,]$Diam_mm_01,na.rm=TRUE)
mean(dat[dat$Species=="DOVI" & trees_to_use,]$Height_cm_01,na.rm=TRUE)

#############################################################
#############################################################
#############################################################