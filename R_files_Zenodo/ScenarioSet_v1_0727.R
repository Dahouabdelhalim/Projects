library(readxl)

ScenarioSetTES <- read_excel('./Data/202107/ScenarioSet202107.xlsx',sheet = 'TES',"B1:J200",col_names = T)
ScenarioSetLCS <- read_excel('./Data/202107/ScenarioSet202107.xlsx',sheet = 'LCS',"B1:J200",col_names = T)
ScenarioSetREA <- read_excel('./Data/202107/ScenarioSet202107.xlsx',sheet = 'REA',"B1:J200",col_names = T)
ScenarioSetBAU <- read_excel('./Data/202107/ScenarioSet202107.xlsx',sheet = 'BAU',"B1:J200",col_names = T)

ScenarioSetTES <- ScenarioSetTES[-which(is.na(ScenarioSetTES[,1])),]
ScenarioSetLCS <- ScenarioSetLCS[-which(is.na(ScenarioSetLCS[,1])),]
ScenarioSetREA <- ScenarioSetREA[-which(is.na(ScenarioSetREA[,1])),]
ScenarioSetBAU <- ScenarioSetBAU[-which(is.na(ScenarioSetBAU[,1])),]

ScenarioSetTES_ad <- paste(ScenarioSetTES$Sce,"_LR",ScenarioSetTES$LR,"_D",ScenarioSetTES$D,"_S",ScenarioSetTES$SNum,"_T",ScenarioSetTES$T,
                           "_M",ScenarioSetTES$M,"_E",ScenarioSetTES$E,"_",ScenarioSetTES$GG,sep = "")
ScenarioSetLCS_ad <- paste(ScenarioSetLCS$Sce,"_LR",ScenarioSetLCS$LR,"_D",ScenarioSetLCS$D,"_S",ScenarioSetLCS$SNum,"_T",ScenarioSetLCS$T,
                           "_M",ScenarioSetLCS$M,"_E",ScenarioSetLCS$E,"_",ScenarioSetLCS$GG,sep = "")
ScenarioSetREA_ad <- paste(ScenarioSetREA$Sce,"_LR",ScenarioSetREA$LR,"_D",ScenarioSetREA$D,"_S",ScenarioSetREA$SNum,"_T",ScenarioSetREA$T,
                           "_M",ScenarioSetREA$M,"_E",ScenarioSetREA$E,"_",ScenarioSetREA$GG,sep = "")
ScenarioSetBAU_ad <- paste(ScenarioSetBAU$Sce,"_LR",ScenarioSetBAU$LR,"_D",ScenarioSetBAU$D,"_S",ScenarioSetBAU$SNum,"_T",ScenarioSetBAU$T,
                           "_M",ScenarioSetBAU$M,"_E",ScenarioSetBAU$E,"_",ScenarioSetBAU$GG,sep = "")

ScenarioSet <- c(ScenarioSetBAU_ad,ScenarioSetTES_ad,ScenarioSetLCS_ad,ScenarioSetREA_ad)



TT = 78
Recover_rate_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)

Recover_CN_rate_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)
Recover_US_rate_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)
Recover_EU_rate_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)

XT_CN_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)
XT_US_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)
XT_EU_sum = matrix(seq(1,TT,length = TT),nrow = 1, ncol = TT)

sce_output <- "22-32-senana1" #Server AL
temp_sw <- as.character(Sys.Date())
Scefile_name <- paste('./Output/202107/Sce-',sce_output,"-",temp_sw,sep = "")
ScefileXT_name <- paste('./Output/202107/Sce-',sce_output,"-",temp_sw,"/XT",sep = "")

RemainVariable <- c("Scenario","ScenarioSet","Scefile_name","ScefileXT_name",
  "Recover_rate_sum","Recover_CN_rate_sum","Recover_US_rate_sum","Recover_EU_rate_sum","Recover_RR_rate_sum",
  "XT_CN_sum","XT_US_sum","XT_EU_sum",
  "Demand_p","Demand_r_RR","Demand_r_PK","Demand_r_TS","Demand_r_RS","Demand_r_OT_CN",
  "Labor_p","Labor_r","Ef","Sect_Fra","Sect_Key","Sect_Fra_self","Sect_Key_self")

