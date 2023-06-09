# ***************************************************************************************
# GSEARIOModel: Version 2.0 by Yaxin Zhang
# GSEARIOModel is extended based on Wang,D (2020),DaopingW/economic-impact-model: Disaster Footprint Model (Version 1.0) [Source code].https://doi.org/10.5281/zenodo.4290117.
# ***************************************************************************************

# Appendix: Scenario Set
# this file will read the manual profile data from "./Data/ScenarioSet.xlsx"
# to create a full scenario set
# called from GSEARIOModel_Main
# coding with UTF-8

# Load data -------------------------------------------------------------

library(readxl)

ScenarioSet22 <- read_excel('./Data/ScenarioSet.xlsx',sheet = 'test',"B1:J300",col_names = T)
ScenarioSetBAU <- read_excel('./Data/ScenarioSet.xlsx',sheet = 'BAU',"B1:J300",col_names = T)

ScenarioSet22 <- ScenarioSet22[-which(is.na(ScenarioSet22[,1])),]
ScenarioSetBAU <- ScenarioSetBAU[-which(is.na(ScenarioSetBAU[,1])),]

ScenarioSet22_ad <- paste(ScenarioSet22$Sce,"_LR",ScenarioSet22$LR,"_D",ScenarioSet22$D,"_S",ScenarioSet22$SNum,"_T",ScenarioSet22$T,
                          "_M",ScenarioSet22$M,"_E",ScenarioSet22$E,"_",ScenarioSet22$GG,sep = "")
ScenarioSetBAU_ad <- paste(ScenarioSetBAU$Sce,"_LR",ScenarioSetBAU$LR,"_D",ScenarioSetBAU$D,"_S",ScenarioSetBAU$SNum,"_T",ScenarioSetBAU$T,
                           "_M",ScenarioSetBAU$M,"_E",ScenarioSetBAU$E,"_",ScenarioSetBAU$GG,sep = "")

ScenarioSet <- c(ScenarioSetBAU_ad,ScenarioSet22_ad)



