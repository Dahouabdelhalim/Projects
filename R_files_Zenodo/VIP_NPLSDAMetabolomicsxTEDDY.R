###########################################################
##########   VIP_NPLSDAMetabolomicsxTEDDY.R     ###########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the metabolomics dataset to calculate the NPLSDA and VIP selection

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()
###########################################################
# Functions:
source ("/media/data/leobalzano/ScriptsForTEDDY/Scripts/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data:
# Metabolomics
mp149_GCTOF <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/mp149_GCTOF.csv")
MP149_NEGATIVE_LIPIDOMICS <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/MP149_NEGATIVE_LIPIDOMICS.csv")
MP149_POSITIVE_LIPIDOMICS <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/MP149_POSITIVE_LIPIDOMICS.csv")

# Checkpoint 1:
POSLIPXw <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/POSLIPXw.txt")
NEGLIPXw <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/NEGLIPXw.txt")
GCTOFXw <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/GCTOFXw.txt")

# Checkpoint 2: 
load("/media/data/leobalzano/ScriptsForTEDDY/Data/arrayGCTOF.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/arrayNegLip.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/arrayPosLip.RData")

# Checkpoint 3: 
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayGCTOFsizeGE.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayPosLipsizeGE.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayNegLipsizeGE.RData")

# Checkpoint 4: 
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayGCTOF136.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayNegLip136.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/fullarrayPosLip136.RData")

# Checkpoint 5: 
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayUNIONGCTOFSelVars.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayINTERSECTNegLipSelVars.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayUNIONPosLipSelVars.RData")

# Response Variable
load("/media/data/leobalzano/ScriptsForTEDDY/Data/Outcomedummyarray.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/outcomedummyarray136.RData")

# Misc
metabolomicsdatatotalpieceresolved <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/metabolomicsdatatotalpieceresolved.csv")
CaseIAnew <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/CaseIAnew")
DEMOGRAPHIC160424 <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/DEMOGRAPHIC160424.csv")
HLA_DIAGNOSIS20160418 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/HLA_DIAGNOSIS20160418")

# List of Cases with at least 3 out of 5 time points with data
patients3tps <- read.table("/media/data/leobalzano/ScriptsForTEDDY/Data/patients3tps.txt", quote="\\"", comment.char="")

###########################################################
# Libraries:
require("abind")
library(data.table)
library("gdata")
require("VennDiagram")

###########################################################
mp149_GCTOF <- read.csv("~/Documents/TEDDY data/MP149/Diet/Metabolomics/mp149_GCTOF.csv")
GCTOF<-mp149_GCTOF
dim (GCTOF)    #11560 x 366
GCTOF[1:10,1:10]
GCTOF[1:10,365:366]
GCTOF2<-as.data.frame(table(GCTOF$Vial_Barcode_Number))

GCTOF3<-GCTOF2[GCTOF2$Freq>1,]
GCTOF2[GCTOF2$Var1==1441709,]
colnames(GCTOF)[1:2]<-c("Mask.Id","Sample.Mask.Id")



MP149_NEGATIVE_LIPIDOMICS <- read.csv("~/Documents/TEDDY data/MP149/Diet/Metabolomics/MP149_NEGATIVE_LIPIDOMICS.csv")
dataneg<-MP149_NEGATIVE_LIPIDOMICS
dim (dataneg)    #11587 x 445 Ahora 11560 x 445
dataneg[1:10,1:10]
dataneg[1:10,444:445]

chirulito<-as.data.frame(table(dataneg$sample_mask_id))
chirulito2<- chirulito[chirulito$Freq>1,]

MP149_POSITIVE_LIPIDOMICS <- read.csv("~/Documents/TEDDY data/MP149/Diet/Metabolomics/MP149_POSITIVE_LIPIDOMICS.csv")
dataPOS<-MP149_POSITIVE_LIPIDOMICS
dim (dataPOS)    #11581 x 516 Ahora 11560 x 516
dataPOS[1:10,1:10]
dataPOS[1:10,515:516]

chirulo<-as.data.frame(table(dataPOS$sample_mask_id))
chirulo2<-chirulo[chirulo$Freq>1,]
dataPOS[dataPOS$sample_mask_id==1441709,1:10] # Ahora no se repite
dataPOS[dataPOS$mask_id==411401,1:10]


MP149_BARCODE_MAPPING_2016.08.29_14.05.14 <- read.csv("~/Documents/TEDDY data/MP149/Diet/Metabolomics/MP149_BARCODE_MAPPING_2016-08-29_14-05-14.csv")
head(MP149_BARCODE_MAPPING_2016.08.29_14.05.14)
colnames(dataneg)[1:2]<-c("Mask.Id","Sample.Mask.Id")
dataneg[1:10,1:4]


colnames(dataPOS)[1:2]<-c("Mask.Id","Sample.Mask.Id")
dataPOS[1:10,1:4]
data2POS<- merge (MP149_BARCODE_MAPPING_2016.08.29_14.05.14, dataPOS, by=c("Sample.Mask.Id","Mask.Id"))
dim(dataPOS); dim(MP149_BARCODE_MAPPING_2016.08.29_14.05.14);dim(data2POS)
summary(data2POS$Draw.Dte.Agedys)
sum(is.na(data2POS$Draw.Dte.Agedys)); sum(is.na(data2POS$Mask.Id)); sum(is.na(data2POS$Sample.Mask.Id))
colnames(data2POS)
dim(data2POS)
data2POS[1:10,1:7]

data2neg<- merge (MP149_BARCODE_MAPPING_2016.08.29_14.05.14, dataneg, by=c("Sample.Mask.Id","Mask.Id"))
dim(dataneg); dim(MP149_BARCODE_MAPPING_2016.08.29_14.05.14); dim(data2neg)
sum(is.na(data2neg$Draw.Dte.Agedys)); sum(is.na(data2neg$Mask.Id)); sum(is.na(data2neg$Sample.Mask.Id))
colnames(data2neg)
dim(data2neg)
data2neg[1:10,1:7]

data2GCTOF<- merge (MP149_BARCODE_MAPPING_2016.08.29_14.05.14, GCTOF, by=c("Sample.Mask.Id","Mask.Id"))
dim(GCTOF); dim(MP149_BARCODE_MAPPING_2016.08.29_14.05.14); dim(data2GCTOF)
sum(is.na(data2GCTOF$Draw.Dte.Agedys)); sum(is.na(data2GCTOF$Mask.Id)); sum(is.na(data2GCTOF$Sample.Mask.Id))
colnames(data2GCTOF)
dim(data2GCTOF)
data2GCTOF[1:10,1:7]

###########################################################
#############     Creation of Agemonth     ################
###########################################################
# Once resolved the Agemonth binning calculation we pull it and merge it with our original data
colnames(metabolomicsdatatotalpieceresolved)
agemonthtable<-metabolomicsdatatotalpieceresolved[,c(1,2,6,7)]
head(agemonthtable)

# Positive lipidomics
data3POS<-merge(agemonthtable,data2POS, by= c("Mask.Id","Draw.Dte.Agedys","Sample.Mask.Id"))
dim(agemonthtable);dim(data2POS);dim(data3POS)
colnames(data3POS[1:10])
colnames(data3POS)[4]<- "Agemonth"
data3POS[1:15,1:5]

# Negative lipidomics
data3neg<-merge(agemonthtable,data2neg, by= c("Mask.Id","Draw.Dte.Agedys","Sample.Mask.Id"))
dim(agemonthtable);dim(data2neg);dim(data3neg)
colnames(data3neg[1:10])
colnames(data3neg)[4]<- "Agemonth"
data3neg[1:15,1:5]

# GCTOF
data3GCTOF
data3GCTOF<-merge(agemonthtable,data2GCTOF, by= c("Mask.Id","Draw.Dte.Agedys","Sample.Mask.Id"))
dim(agemonthtable);dim(data2GCTOF);dim(data3GCTOF)
colnames(data3GCTOF[1:10])
colnames(data3GCTOF)[4]<- "Agemonth"
data3GCTOF[1:15,1:5]

###########################################################
################      Merge with CasesIA     ##############
###########################################################
CaseIAnew<-cbind(CaseIAnew, Diagnosis = "IA")
head(CaseIAnew)
dim(CaseIAnew)  #1671 x 5

# Positive lipidomics
data3POS[data3POS$Mask.Id== 202138,1:10]

data4POS<-merge(CaseIAnew,data3POS, by="Mask.Id")#, all.y = T)   202138
dim(CaseIAnew);dim(data3POS);dim(data4POS)  #10530
sum(is.na(data4POS$Mask.Id))
sum(is.na(data4POS$Outcome))
data4POS[1:100,c(1:10,50)]
data4POS[data4POS$Mask.Id== 202138,1:10]


# Negative lipidomics
data4neg<-merge(CaseIAnew,data3neg, by="Mask.Id")
dim(CaseIAnew);dim(data3neg);dim(data4neg)  #10530
sum(is.na(data4neg$Mask.Id))
sum(is.na(data4neg$Outcome))
data4neg[1:100,c(1:9,50)]
data4neg[data4neg$Mask.Id== 202138,1:10]   # No debe estar presente



# GCTOF
data4GCTOF<-merge(CaseIAnew,data3GCTOF, by="Mask.Id")
dim(CaseIAnew);dim(data3GCTOF);dim(data4GCTOF)  #10530
sum(is.na(data4GCTOF$Mask.Id))
sum(is.na(data4GCTOF$Outcome))
data4GCTOF[1:100,c(1:9,50)]
data4GCTOF[data4GCTOF$Mask.Id== 202138,1:10]   # No debe estar presente

###########################################################
##########     Merge with demographics data     ###########
###########################################################
dim(DEMOGRAPHIC160424)
dimnames(DEMOGRAPHIC160424)

# Positive lipidomics
data5POS<-merge(DEMOGRAPHIC160424,data4POS, by="Mask.Id")
dim(DEMOGRAPHIC160424);dim(data4POS);dim(data5POS)  #10530
data5POS[1:50,c(1:10,449)]
dim(data5POS)
summary(data5POS$Gender)
sum(is.na(data5POS$Gender))
sum(is.na(data5POS$Mask.Id))
sum(is.na(data5POS$Outcome))

# Negative lipidomics
data5neg<-merge(DEMOGRAPHIC160424,data4neg, by="Mask.Id")
dim(DEMOGRAPHIC160424);dim(data4neg);dim(data5neg)  #10530
data5neg[1:50,c(1:10,449)]
dim(data5neg)
summary(data5neg$Gender)
sum(is.na(data5neg$Gender))
sum(is.na(data5neg$Mask.Id))
sum(is.na(data5neg$Outcome))


# GCTOF
data5GCTOF<-merge(DEMOGRAPHIC160424,data4GCTOF, by="Mask.Id")
dim(DEMOGRAPHIC160424);dim(data4GCTOF);dim(data5GCTOF)  #10530
data5GCTOF[1:50,c(1:10,449)]
dim(data5GCTOF)
summary(data5GCTOF$Gender)
sum(is.na(data5GCTOF$Gender))
sum(is.na(data5GCTOF$Mask.Id))
sum(is.na(data5GCTOF$Outcome))

###########################################################
############     Merging with HLA Diagnosis     ###########
###########################################################

# Positive lipidomics
data6POS<-merge(HLA_DIAGNOSIS20160418,data5POS,by="Mask.Id", all.y = T)
dim(data5POS); dim(HLA_DIAGNOSIS20160418);dim(data6POS)
sum(is.na(data6POS$Persist.Conf.Gad))
data6POS[1:50,c(1:15,22)]
sum(is.na(data6POS$Case.Endptage))

# Negative lipidomics
data6neg<-merge(HLA_DIAGNOSIS20160418,data5neg,by="Mask.Id", all.y = T)
dim(data5neg); dim(HLA_DIAGNOSIS20160418);dim(data6neg)
sum(is.na(data6neg$Persist.Conf.Gad))
data6neg[1:50,c(1:15,22)]
sum(is.na(data6neg$Case.Endptage))

# GCTOF
data6GCTOF<-merge(HLA_DIAGNOSIS20160418,data5GCTOF,by="Mask.Id", all.y = T)
dim(data5GCTOF); dim(HLA_DIAGNOSIS20160418);dim(data6GCTOF)
sum(is.na(data6GCTOF$Persist.Conf.Gad))
data6GCTOF[1:50,c(1:15,22)]
sum(is.na(data6GCTOF$Case.Endptage))

#########################################################################################
########     Converting Case.Endptage(22) to monthly, dividing by 30     ################
#########################################################################################
# Positive lipidomics
sum(is.na(data6POS$Case.Endptage))
colnames (data6POS[1:45])
data6POS[,22] <- round(data6POS[,22]/30)
data6POS[1:30,c(2,22:26,29)]

# Negative lipidomics
sum(is.na(data6neg$Case.Endptage))
colnames (data6neg[1:45])
data6neg[,22] <- round(data6neg[,22]/30)
data6neg[1:30,c(2,22:26,29)]

# GCTOF
sum(is.na(data6GCTOF$Case.Endptage))
colnames (data6GCTOF[1:45])
data6GCTOF[,22] <- round(data6GCTOF[,22]/30)
data6GCTOF[1:30,c(2,22:26,29)]


#########################################################################################
################     Order tables by Outcome (23) and Agemonth(27)     ##################
#########################################################################################
# Positive lipidomics
data6POS<-data6POS[with(data6POS, order(-data6POS[,23],data6POS[,27])),]  
data6POS[1:6,c(2,22:27)]

# Negative lipidomics
data6neg<-data6neg[with(data6neg, order(-data6neg[,23],data6neg[,27])),]  
data6neg[1:6,c(2,22:27)]

# GCTOF
data6GCTOF<-data6GCTOF[with(data6GCTOF, order(-data6GCTOF[,23],data6GCTOF[,27])),]  
data6GCTOF[1:6,c(2,22:27)]



###########################################################################################
###     Create column that is the substraction of Agemonth(27) - enpointStage (22),     ###
###     now i have the 0;-3; -6; -9; -12                                                ###
###########################################################################################
# Seroconversion column
# Positive lipidomics
sum(is.na(data6POS$Agemonth))   #To check for any error Value have to be 0
sum(is.na(data6POS$Mask.Id))
sum(is.na(data6POS$Case.Endptage)) #have to be 0
sum(is.na(data6POS$CASE.IND))

colnames(data6POS[1:40])
Sconv <- data6POS[,27]-data6POS[,22]
Sconv
data7POS<-cbind(Sconv,data6POS)
data7POS[1:5,1:30]
dim(data7POS)     # 947700  x  1351
sum(is.na(data7POS$Sconv))
sum(is.na(data7POS$Agemonth))
sum(is.na(data7POS$Case.Endptage))   # 0

completeIAPositivelipidomics<- data7POS

# Negative lipidomics
sum(is.na(data6neg$Agemonth))   #To check for any error Value have to be 0
sum(is.na(data6neg$Mask.Id))
sum(is.na(data6neg$Case.Endptage)) #have to be 0
sum(is.na(data6neg$CASE.IND))

colnames(data6neg[1:40])
Sconv <- data6neg[,27]-data6neg[,22]
Sconv
data7neg<-cbind(Sconv,data6neg)
data7neg[1:5,1:30]
dim(data7neg)     # 947700  x  1351
sum(is.na(data7neg$Sconv))
sum(is.na(data7neg$Agemonth))
sum(is.na(data7neg$Case.Endptage))   # 0

completeIANegativelipidomics<-data7neg

# GCTOF
sum(is.na(data6GCTOF$Agemonth))   #To check for any error Value have to be 0
sum(is.na(data6GCTOF$Mask.Id))
sum(is.na(data6GCTOF$Case.Endptage)) #have to be 0
sum(is.na(data6GCTOF$CASE.IND))

colnames(data6GCTOF[1:40])
Sconv <- data6GCTOF[,27]-data6GCTOF[,22]
Sconv
data7GCTOF<-cbind(Sconv,data6GCTOF)
data7GCTOF[1:5,1:30]
dim(data7GCTOF)     # 947700  x  1351
sum(is.na(data7GCTOF$Sconv))
sum(is.na(data7GCTOF$Agemonth))
sum(is.na(data7GCTOF$Case.Endptage))   # 0
completeIAGCTOF<- data7GCTOF

###########################################################################################
###########     Order table by Sconv(1), CASE.IND(22) and -Outcome(24)     ################
###########################################################################################
# Positive lipidomics
colnames(completeIAPositivelipidomics[1:40])
completeIAPositivelipidomics<-completeIAPositivelipidomics[with(completeIAPositivelipidomics, order(completeIAPositivelipidomics[,1],completeIAPositivelipidomics[,22],-completeIAPositivelipidomics[,24])),]  #Ordered by Sconv, CASE.IND and Outcome
completeIAPositivelipidomics[1:20,c(1:2,22:30)]

# Negative lipidomics
colnames(completeIANegativelipidomics[1:40])
completeIANegativelipidomics<-completeIANegativelipidomics[with(completeIANegativelipidomics, order(completeIANegativelipidomics[,1],completeIANegativelipidomics[,22],-completeIANegativelipidomics[,24])),]  #Ordered by Sconv, CASE.IND and Outcome
completeIANegativelipidomics[1:20,c(1:2,22:30)]

# GCTOF
colnames(completeIAGCTOF[1:40])
completeIAGCTOF<-completeIAGCTOF[with(completeIAGCTOF, order(completeIAGCTOF[,1],completeIAGCTOF[,22],-completeIAGCTOF[,24])),]  #Ordered by Sconv, CASE.IND and Outcome
completeIAGCTOF[1:20,c(1:2,22:30)]


###########################################################################################
######################     Eliminate just cases with Sconv <-13,     ######################
###########################################################################################
# Positive lipidomics
POSLIP12m<-completeIAPositivelipidomics[!(completeIAPositivelipidomics$Sconv<(-13)),]
POSLIP12m[1:50,1:5]
dim(POSLIP12m)  #6382 for 544 months

# Negative lipidomics
NEGLIP12m<-completeIANegativelipidomics[!(completeIANegativelipidomics$Sconv<(-13)),]
NEGLIP12m[1:50,1:5]
dim(NEGLIP12m)  #6382 for 473 months

# GCTOF
GCTOF12m<-completeIAGCTOF[!(completeIAGCTOF$Sconv<(-13)),]
GCTOF12m[1:50,1:5]
dim(GCTOF12m)  #6382 for 394 months

###########################################################################################
##########################     Eliminate cases with Sconv >0,     #########################
###########################################################################################
# Positive lipidomics
POSLIP12m2<-POSLIP12m[!(POSLIP12m$Sconv>(1)),]
POSLIP12m2[1:50,1:5]
dim(POSLIP12m2)  #5300 for 544 months

# Negative lipidomics
NEGLIP12m2<-NEGLIP12m[!(NEGLIP12m$Sconv>(1)),]
NEGLIP12m2[1:50,1:5]
dim(NEGLIP12m2)  #5300 for 473 months

# GCTOF
GCTOF12m2<-GCTOF12m[!(GCTOF12m$Sconv>(1)),]
GCTOF12m2[1:50,1:5]
dim(GCTOF12m2)  #5300 for 394 months

###########################################################################################
########################     Regroup data in the 3 time points     ########################
###########################################################################################
# Positive lipidomics
POSLIP12m2$Sconv[POSLIP12m2$Sconv==-13] <- -12
POSLIP12m2$Sconv[POSLIP12m2$Sconv==-11] <- -12

POSLIP12m2$Sconv[POSLIP12m2$Sconv==-10] <- -9
POSLIP12m2$Sconv[POSLIP12m2$Sconv==-8] <- -9

POSLIP12m2$Sconv[POSLIP12m2$Sconv==-7] <- -6
POSLIP12m2$Sconv[POSLIP12m2$Sconv==-5] <- -6

POSLIP12m2$Sconv[POSLIP12m2$Sconv==-4] <- -3
POSLIP12m2$Sconv[POSLIP12m2$Sconv==-2] <- -3

POSLIP12m2$Sconv[POSLIP12m2$Sconv==-1] <- 0
POSLIP12m2$Sconv[POSLIP12m2$Sconv== 1] <- 0


# Negative lipidomics
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-13] <- -12
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-11] <- -12

NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-10] <- -9
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-8] <- -9

NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-7] <- -6
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-5] <- -6

NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-4] <- -3
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-2] <- -3

NEGLIP12m2$Sconv[NEGLIP12m2$Sconv==-1] <- 0
NEGLIP12m2$Sconv[NEGLIP12m2$Sconv== 1] <- 0


# GCTOF
GCTOF12m2$Sconv[GCTOF12m2$Sconv==-13] <- -12
GCTOF12m2$Sconv[GCTOF12m2$Sconv==-11] <- -12

GCTOF12m2$Sconv[GCTOF12m2$Sconv==-10] <- -9
GCTOF12m2$Sconv[GCTOF12m2$Sconv==-8] <- -9

GCTOF12m2$Sconv[GCTOF12m2$Sconv==-7] <- -6
GCTOF12m2$Sconv[GCTOF12m2$Sconv==-5] <- -6

GCTOF12m2$Sconv[GCTOF12m2$Sconv==-4] <- -3
GCTOF12m2$Sconv[GCTOF12m2$Sconv==-2] <- -3

GCTOF12m2$Sconv[GCTOF12m2$Sconv==-1] <- 0
GCTOF12m2$Sconv[GCTOF12m2$Sconv== 1] <- 0


GCTOF12m2[GCTOF12m2$Mask.Id==714949,c(1,2,21:40)] # Examples of Cases used as controls
GCTOF12m2[GCTOF12m2$Mask.Id==918646,c(1,2,21:40)]

####################################################################################################################################
##################     Create the Case-Control Groups based on the separation due to Agemonth and CASE.IND      ####################
####################################################################################################################################


colnames(POSLIP12m2[1:40])
POSLIP12m2<-POSLIP12m2[with(POSLIP12m2, order(POSLIP12m2[,22], POSLIP12m2[,28], -POSLIP12m2[,24])),]  #Ordered by CASE.IND(22), Agemonth(28) and -Outcome (24)
NEGLIP12m2<-NEGLIP12m2[with(NEGLIP12m2, order(NEGLIP12m2[,22], NEGLIP12m2[,28], -NEGLIP12m2[,24])),]  #Ordered by CASE.IND(22), Agemonth(28) and -Outcome (24)
GCTOF12m2<-GCTOF12m2[with(GCTOF12m2, order(GCTOF12m2[,22], GCTOF12m2[,28], -GCTOF12m2[,24])),]  #Ordered by CASE.IND(22), Agemonth(28) and -Outcome (24)

colnames(POSLIP12m2[1:35])
POSLIP12m2[1:30,c(1,2,22,24:27,29,30,31)]
NEGLIP12m2[1:30,c(1,2,22,24:27,29,30,31)]
GCTOF12m2[1:30,c(1,2,22,24:27,29,30,31)]

identical(POSLIP12m2$Mask.Id,NEGLIP12m2$Mask.Id)
identical(POSLIP12m2$Mask.Id,GCTOF12m2$Mask.Id)
if(all(POSLIP12m2$Mask.Id == NEGLIP12m2$Mask.Id) & all(POSLIP12m2$Mask.Id == GCTOF12m2$Mask.Id))
{
  print("All Mask.ID's are identical")
}

# Since they are identical, we are going to create just one NewGroupsMetabolomics and then paste it in all three tables
# Create just one group set, in this case due to the variations of Agemonth, the best way to do it is by
# using CASE.IND (22) and Sconv (1) nevertheless it is not like previous cases, this is very variable
colnames(POSLIP12m2[1:35])

gruposMetabolomics <-apply(POSLIP12m2[,c(22,1)], 1, paste, collapse ="_")
head(gruposMetabolomics)
unigruposMetabolomics <- unique (gruposMetabolomics)
length(unigruposMetabolomics)   # 1560 groups

# Positive lipidomics
POSLIP12m3 <- cbind (gruposMetabolomics, POSLIP12m2)

# Negative lipidomics
NEGLIP12m3 <- cbind (gruposMetabolomics, NEGLIP12m2)

# GCTOF
GCTOF12m3 <- cbind (gruposMetabolomics, GCTOF12m2)

POSLIP12m3[1:50,c(1,6:7,23,25,30,31)]   #datos(1212 rows)Agemonth (30)

#######
NewGroupsDB=NULL
n=0
Nn ="word"
grupin=NULL
for (i in 1:nrow(POSLIP12m3)) {
  if (POSLIP12m3[i,1] == Nn) {
    grupin = n
  }
  if (POSLIP12m3[i,1] != Nn) {
    Nn = POSLIP12m3[i,1]
    n=n+1
    grupin = n
  }
  NewGroupsDB=rbind(NewGroupsDB,as.vector(grupin))
  
}
as.vector(NewGroupsDB)
dim(NewGroupsDB)    #5300
NewGroupsMetabolomics<-NewGroupsDB

# Positive lipidomics
POSLIP12m4 <- cbind (NewGroupsMetabolomics, POSLIP12m3)
POSLIP12m4[1:10,c(1,2,3,4,5,6:8,24:26,31,32)]   
POSLIP12m4[1:100,c(1,2,3,4,24:26,32)]   

# Negative lipidomics
NEGLIP12m4 <- cbind (NewGroupsMetabolomics, NEGLIP12m3)
NEGLIP12m4[1:10,c(1,2,3,4,5,6:8,24:26,31,32)]   
NEGLIP12m4[1:100,c(1,2,3,4,24:26,32)]   

# GCTOF
GCTOF12m3
GCTOF12m4 <- cbind (NewGroupsMetabolomics, GCTOF12m3)
GCTOF12m4[1:10,c(1,2,3,4,5,6:8,24:26,31,32)]   
GCTOF12m4[1:100,c(1,2,3,4,24:26,32)]   

##############################################################################################################################
###  Creation of AgeGroups which is the variable that groups individuals in terms of age to try to determine if there is
###  any difference of disease appearance due to the age of the kid
###  We decide to use the following frames: Year1<=13,Year2<=25,Year3<=37,Year4>37
##############################################################################################################################
identical(POSLIP12m4$Case.Endptage,NEGLIP12m4$Case.Endptage)
identical(POSLIP12m4$Case.Endptage,GCTOF12m4$Case.Endptage)

# Since the Case.Endptage columns are identical, we are going to create just one AgeGroup and then paste it in all three tables
POSLIP12m4[3085:3094,1:30]


#### Creation of the variable AgeGroup
AgeGroup = NULL

for(i in 1:nrow(POSLIP12m4)) {
  
  grupin = NULL
  if(POSLIP12m4$Case.Endptage[i]>0 & POSLIP12m4$Case.Endptage[i]<=13){
    grupin = 1
  }
  if(POSLIP12m4$Case.Endptage[i]>13 & POSLIP12m4$Case.Endptage[i]<=25){ 
    grupin = 2
  }
  if(POSLIP12m4$Case.Endptage[i]>25 & POSLIP12m4$Case.Endptage[i]<=37){ 
    grupin = 3
  }
  if(POSLIP12m4$Case.Endptage[i]>37){ 
    grupin = 4
  }
  AgeGroup=c(AgeGroup,grupin)
}

dim(AgeGroup)
as.vector(AgeGroup)

# Positive lipidomics
POSLIP12m5<-cbind(AgeGroup,POSLIP12m4)
POSLIP12m5[1:200,c(1:6,25,26)]
dim(POSLIP12m5)  #53007 x 1354     for 13 months

# Negative lipidomics
NEGLIP12m5<-cbind(AgeGroup,NEGLIP12m4)
NEGLIP12m5[1:200,c(1:6,25,26)]
dim(NEGLIP12m5)  #53007 x 1354     for 13 months

# GCTOF
GCTOF12m5<-cbind(AgeGroup,GCTOF12m4)
GCTOF12m5[1:200,c(1:6,25,26)]
dim(GCTOF12m5)  #53007 x 1354     for 13 months

#############################################################################################################################
#######     Creating First AAb that appears as a categorical variable related to case controls groups     ###################
#############################################################################################################################
identical(POSLIP12m4$Days.To.Persist.Conf.Gad,NEGLIP12m4$Days.To.Persist.Conf.Gad)
identical(POSLIP12m4$Days.To.Persist.Conf.Gad,GCTOF12m4$Days.To.Persist.Conf.Gad)
# Similarly the columns are the same so i am going to create the column and then paste ir to each table

POSLIP12m5[1:10,c(1:6,10:12)]
ata8<-POSLIP12m5

ata8[,10:12][is.na(ata8[,10:12])] <- 9999
ata8[1:10,c(1:2,10:12,26)]

#############################################################################################################################
###########   Creation of the column with a for loop
# Controls must have the value of the first AAb that the case have
# So this have to be done in two steps, first, the AAb in the cases and 0 in controls, 
# Second: the controls assume the first AAb of the case they are attached to.
# Control = 0; GADA = 1; IAA = 2; IA2A = 3; IAA/GADA/IA2A = 4; IAA/GADA = 5; IAA/IA2A = 6; GADA/IA2A = 7.

firstAAb = NULL

for(i in 1:nrow(ata8)) {
  AAb = NULL
  
  if((ata8[i,10]==ata8[i,11]) & (ata8[i,10]==ata8[i,12]) & (ata8[i,10] == 9999)){
    AAb = 0
  }
  if((ata8[i,10]<ata8[i,11]) & (ata8[i,10]<ata8[i,12])){
    AAb = 1
  }
  if((ata8[i,11]<ata8[i,10]) & (ata8[i,11]<ata8[i,12])){ 
    AAb = 2
  }
  if((ata8[i,12]<ata8[i,10]) & (ata8[i,12]<ata8[i,11])){ 
    AAb = 3
  }
  if((ata8[i,10]==ata8[i,11]) & (ata8[i,10]==ata8[i,12]) & (ata8[i,10] != 9999)){ 
    AAb = 4 
  }
  if((ata8[i,10]==ata8[i,11]) & (ata8[i,10]<ata8[i,12])){ 
    AAb = 5
  }
  if((ata8[i,11]==ata8[i,12]) & (ata8[i,11]<ata8[i,10])){ 
    AAb = 6
  }
  if((ata8[i,10]==ata8[i,12]) & (ata8[i,10]<ata8[i,11])){ 
    AAb = 7
  }
  firstAAb=rbind(firstAAb, AAb)
  
}

dim(firstAAb)
as.vector(firstAAb)

ata8[,10:12][ata8[,10:12]==9999] <- NA

# Positive lipidomics
POSLIP12m6<-cbind(as.vector(firstAAb), POSLIP12m5)
colnames(POSLIP12m6)[1]<-"FirstAAb"
POSLIP12m6[3235:3255,c(1:6,11:13,27,28,29,30)]
dim(POSLIP12m6)   #5300   x  1355
POSLIP12m6[1:50,c(1:6,12:14,27,28,29,30)]
POSLIP12m6[POSLIP12m6$Mask.Id==756775,c(1:6,12:14,27,28,29,30)] # Example of control that after a while becomes a case
colnames(POSLIP12m6[1:35])

# Negative lipidomics
NEGLIP12m6<-cbind(as.vector(firstAAb), NEGLIP12m5)
colnames(NEGLIP12m6)[1]<-"FirstAAb"
NEGLIP12m6[3235:3255,c(1:6,11:13,27,28,29,30)]
colnames(NEGLIP12m6[1:35])

# GCTOF
GCTOF12m6<-cbind(as.vector(firstAAb), GCTOF12m5)
colnames(GCTOF12m6)[1]<-"FirstAAb"
GCTOF12m6[3235:3255,c(1:6,11:13,27,28,29,30)]
colnames(GCTOF12m6[1:35])

#########################################################################################################################
# Second step (It has to be done with respect to Outcome)
# For the characteristics of the data, the tables must be ordered again by NewGroupsMetabolomics(3), Outcome (28)
ata9<-POSLIP12m6

ata9<-ata9[with(ata9, order(ata9[,3],-ata9[,28])),]  
ata9[1:50,c(1:6,11:13,28,29,30)]
ata9[1:50,c(1:6,28,29,30)]


firstAAbCC = NULL
AAbCC = NULL
truquin=NULL
#for(i in 1:length(i)) {
for(i in 1:length(ata9$NewGroupsMetabolomics)) {
  
  if(ata9[i,28] != 0) {               # Outcome
    AAbCC = ata9[i,1]
    truquin=ata9[i,1]
  }
  if(ata9[i,28] == 0 ) {              # Outcome
    AAbCC =truquin
  }
  firstAAbCC=rbind(firstAAbCC, AAbCC)
}
dim(firstAAbCC)
as.vector(firstAAbCC)



# Positive lipidomics
POSLIP12m7<-cbind(as.vector(firstAAbCC), POSLIP12m6)
colnames(POSLIP12m7)[1]<-"FirstAAbCC"
POSLIP12m7[1:60,c(1:5,12:14,27,28,29)]
dim(POSLIP12m7)   #5467   x  993
summary(POSLIP12m7$FirstAAbCC)
POSLIP12m7[POSLIP12m7$CASE.IND==47,c(1:4,6,29)]
POSLIP12m7[1:40,c(1,2,3,4,6,29,31,7)]

# Negative lipidomics
NEGLIP12m7<-cbind(as.vector(firstAAbCC), NEGLIP12m6)
colnames(NEGLIP12m7)[1]<-"FirstAAbCC"
NEGLIP12m7[1:60,c(1:5,12:14,27,28,29)]
dim(NEGLIP12m7)   #5467   x  993


# GCTOF

GCTOF12m7<-cbind(as.vector(firstAAbCC), GCTOF12m6)
colnames(GCTOF12m7)[1]<-"FirstAAbCC"
GCTOF12m7[1:60,c(1:5,12:14,27,28,29)]
dim(GCTOF12m7)   #5467   x  993


################################################################################################################################
# In this particular approach, we must delete all controls that eventually become cases, we can detect them because they have
# 0 value in Outcome but a value in FirstAAb
# So the better to do is eliminate them 
################################################################################################################################
# Positive lipidomics
colnames (POSLIP12m7[1:40])
dim(POSLIP12m7) #5300 x 1358
sum(POSLIP12m7$Outcome==0 & POSLIP12m7$FirstAAb!= 0) # 180
eliminatingdatos<-POSLIP12m7[POSLIP12m7$Outcome==0 & POSLIP12m7$FirstAAb!= 0, ]
dim(eliminatingdatos)

POSLIP12m8<- subset(POSLIP12m7,!(POSLIP12m7$Outcome==0 & POSLIP12m7$FirstAAb != 0))
dim(POSLIP12m7);dim(POSLIP12m8)     # 5120 x 1358
IAm12Positivelipidomics<-POSLIP12m8

# Negative lipidomics
colnames (NEGLIP12m7[1:40])
dim(NEGLIP12m7) #5300 x 1358
sum(NEGLIP12m7$Outcome==0 & NEGLIP12m7$FirstAAb!= 0) # 180
eliminatingdatos<-NEGLIP12m7[NEGLIP12m7$Outcome==0 & NEGLIP12m7$FirstAAb!= 0, ]
dim(eliminatingdatos)

NEGLIP12m8<- subset(NEGLIP12m7,!(NEGLIP12m7$Outcome==0 & NEGLIP12m7$FirstAAb != 0))
dim(NEGLIP12m7);dim(NEGLIP12m8)     # 5120 x 1358
IAm12Negativelipidomics<-NEGLIP12m8

# GCTOF
colnames (GCTOF12m7[1:40])
dim(GCTOF12m7) #5300 x 1358
sum(GCTOF12m7$Outcome==0 & GCTOF12m7$FirstAAb!= 0) # 180
eliminatingdatos<-GCTOF12m7[GCTOF12m7$Outcome==0 & GCTOF12m7$FirstAAb!= 0, ]
dim(eliminatingdatos)

GCTOF12m8<- subset(GCTOF12m7,!(GCTOF12m7$Outcome==0 & GCTOF12m7$FirstAAb != 0))
dim(GCTOF12m7);dim(GCTOF12m8)     # 5120 x 1358

IAm12GCTOF<-GCTOF12m8

###########################################################
###     Eliminating the groups with just one member     ###
###########################################################

# Positive lipidomics
setDT(IAm12Positivelipidomics)
POSLIP12m9 <- IAm12Positivelipidomics[,n:=.N,NewGroupsMetabolomics][n>1,,][,n:=NULL]
sort(table(POSLIP12m9$NewGroupsMetabolomics))
sort(table(IAm12Positivelipidomics$NewGroupsMetabolomics))
dim(POSLIP12m9)    #5010 x 1364
class(POSLIP12m9)
POSLIP12m9<-data.frame(POSLIP12m9)
POSLIP12m9[1:10,1:10]

# Negative lipidomics
setDT(IAm12Negativelipidomics)
NEGLIP12m9 <- IAm12Negativelipidomics[,n:=.N,NewGroupsMetabolomics][n>1,,][,n:=NULL]
sort(table(NEGLIP12m9$NewGroupsMetabolomics))
sort(table(IAm12Negativelipidomics$NewGroupsMetabolomics))
dim(NEGLIP12m9)    #5010 x 1364
class(NEGLIP12m9)
NEGLIP12m9<-data.frame(NEGLIP12m9)
NEGLIP12m9[1:10,1:10]


# GCTOF
setDT(IAm12GCTOF)
GCTOF12m9 <- IAm12GCTOF[,n:=.N,NewGroupsMetabolomics][n>1,,][,n:=NULL]
sort(table(GCTOF12m9$NewGroupsMetabolomics))
sort(table(IAm12GCTOF$NewGroupsMetabolomics))
dim(GCTOF12m9)    #5010 x 1364
class(GCTOF12m9)
GCTOF12m9<-data.frame(GCTOF12m9)
GCTOF12m9[1:10,1:10]

###########################################################
####     Eliminating the groups with just Controls     ####
###########################################################
# Positive lipidomics
data<-POSLIP12m9
colnames(data[1:45])
annots <- data[,c(1:35)]  # move annotations


POSLIP12m10<-NULL
new.groups <- unique(annots[,4])
for (j in 1:length(new.groups)){  
  #for (j in 1:10){
  myannots <- annots[annots[,4] == new.groups[j],]
  mydata <- data[annots[,4] == new.groups[j],]  
  if(sum(myannots[,29])>0) {
    misdatos <- mydata
  } else {
    misdatos <- NULL
  }
  
  POSLIP12m10<-rbind(POSLIP12m10,misdatos)
}

dim(POSLIP12m9)
dim(POSLIP12m10)   # 4840 x 549
POSLIP12m10[1:100,c(1:4,6,7,29,514)]
POSLIP12m9[1:100,c(1:4,6,7,29,514)]

IAm12Positivelipidomics<-POSLIP12m10

# Negative lipidomics
data<-NEGLIP12m9
colnames(data[1:45])
annots <- data[,c(1:35)]  # move annotations


NEGLIP12m10<-NULL
new.groups <- unique(annots[,4])
for (j in 1:length(new.groups)){  
  #for (j in 1:4){
  myannots <- annots[annots[,4] == new.groups[j],]
  mydata <- data[annots[,4] == new.groups[j],]  
  if(sum(myannots[,29])>0) {
    misdatos <- mydata
  }else {
    misdatos <- NULL
  }
  NEGLIP12m10<-rbind(NEGLIP12m10,misdatos)
}

dim(NEGLIP12m9)
dim(NEGLIP12m10)
NEGLIP12m10[,c(1:10)]
NEGLIP12m9[1:16,c(1:10)]

IAm12Negativelipidomics<-NEGLIP12m10
dim(IAm12Negativelipidomics)
IAm12Negativelipidomics[1:10,1:10]
# GCTOF
data<-GCTOF12m9
colnames(data[1:45])
annots <- data[,c(1:35)]  # move annotations


GCTOF12m10<-NULL
new.groups <- unique(annots[,4])
for (j in 1:length(new.groups)){  
  #for (j in 1:4){
  myannots <- annots[annots[,4] == new.groups[j],]
  mydata <- data[annots[,4] == new.groups[j],]  
  if(sum(myannots[,29])>0) {
    misdatos <- mydata
  }else {
    misdatos <- NULL
  }
  GCTOF12m10<-rbind(GCTOF12m10,misdatos)
}

dim(GCTOF12m9)
dim(GCTOF12m10)
GCTOF12m10[,c(1:10)]
GCTOF12m9[1:16,c(1:10)]

IAm12GCTOF<-GCTOF12m10


###########################################################
#     Function To calculate the internal substraction     #
###########################################################
# FUNCTIONS
###########
matchAge2 <- function (ages, case, h = 4) {
  ages.1 <- ages[case == 1]
  n = length(ages.1)
  ages.0 <- ages[case == 0]
  order <- c(1: length(ages))
  for ( i in 1: length(ages.1)) {
    abs.dif <- abs(ages.1[i] - ages.0)
    if (min(abs.dif) < h) {
      match <- which(abs.dif == min(abs.dif))
      order[n+match] <- n+i
    }
  }
  order
}

matchAge <- function (ages, case, h = 4) {
  ages.1 <- ages[case == 1]
  ages.0 <- ages[case == 0]
  pairs <- data.frame(case.month = 1:length(ages.1), controls = rep(0, length (ages.1)))
  for ( i in 1: length(ages.1)) {
    abs.dif <- abs(ages.1[i] - ages.0)
    if (min(abs.dif) < h) {
      match <- which(abs.dif == min(abs.dif))[1]
      pairs[i,2] <- match
    }
  }
  pairs <- pairs[pairs[,2] > 0,]
  pairs
}


###########################################################
#  DATA

# Positive lipidomics
IAm12Positivelipidomics[1:100,c(1:4,6,7,29,514)]
data<-IAm12Positivelipidomics
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)

# Calculate Ratios
new.groups <- unique(annots[,4])
tablecita = NULL
tablecitapromediada=NULL
for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #  for (j in 1:5){
  print(paste(j, "out of", length(new.groups)))
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  for ( a in "IA") {
    if (is.element (a, unique(myannots[,30]))){     # Diagnosis(30)
      myannots2 <- myannots[myannots[,30] == a,]    # Diagnosis(30)
      mydata2 <- mydata[myannots[,30] == a,]        # Diagnosis(30)
      myages <- myannots2[,6]                       # Sconv(6)
      mycase <- myannots2[,29]                      # Outcome(29)
      if (any(mycase==1) & any(mycase == 0)) {
        myorder <- matchAge2(myages, mycase, h = 1)  # to average controls
        data3 <- apply(mydata2,2,function (x) tapply(x, myorder, mean))
        data4 <- data.frame(apply(data3,2,mean))
        colnames(data4)<-new.groups[j]
        data5<-data.frame(t(data4))
        myages2 <- tapply(myages, myorder, mean)
        mycase2 <- tapply(mycase, myorder, mean)
        mypairs <- matchAge(myages2, mycase2 , h = 4)  # to get index for computing ratios
        print (paste(a, nrow(myannots), nrow(mypairs)))
      }
      tablecita<- rbind(tablecita,data3)
      tablecitapromediada <-rbind(tablecitapromediada,data5)
      
    }
  }
}


dim(tablecita)
tablecita[,1:30]
dim(tablecitapromediada)
tablecitapromediada[,1:30]

# Function
# Positive lipidomics
dim(tablecitapromediada)
length(unique(IAm12Positivelipidomics[,4]))

data<-IAm12Positivelipidomics
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)
new.groups <- unique(annots[,4])
miu<-rownames(tablecitapromediada)
tablecitasub = NULL

for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #for (j in 1:50){
  print(paste(j, "out of", length(new.groups)))
  print(new.groups[j])
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  miuvalue<-as.matrix(tablecitapromediada[j,])
  sustraccion <-sweep(mydata,2,miuvalue[1,])
  #sustraccion2<-cbind (new.groups[j],sustraccion)
  tablecitasub<- rbind(tablecitasub,sustraccion)
}      

POSLIPXw<-cbind(annots,tablecitasub)
dim(annots)
dim(tablecitasub); dim(POSLIPXw)
POSLIPXw[1:100,c(1:4,6,7,29,514)]
#write.table(POSLIPXw, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/POSLIPXw.txt", sep="\\t")
#POSLIPXw <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/POSLIPXw.txt")
##############################################################################################################################
# Negative lipidomics

data<-IAm12Negativelipidomics
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)

# Calculate Ratios
new.groups <- unique(annots[,4])
tablecita = NULL
tablecitapromediada=NULL
for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #  for (j in 1:5){
  print(paste(j, "out of", length(new.groups)))
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  for ( a in "IA") {
    if (is.element (a, unique(myannots[,30]))){     # Diagnosis(30)
      myannots2 <- myannots[myannots[,30] == a,]    # Diagnosis(30)
      mydata2 <- mydata[myannots[,30] == a,]        # Diagnosis(30)
      myages <- myannots2[,6]                       # Sconv(6)
      mycase <- myannots2[,29]                      # Outcome(29)
      if (any(mycase==1) & any(mycase == 0)) {
        myorder <- matchAge2(myages, mycase, h = 1)  # to average controls
        data3 <- apply(mydata2,2,function (x) tapply(x, myorder, mean))
        data4 <- data.frame(apply(data3,2,mean))
        colnames(data4)<-new.groups[j]
        data5<-data.frame(t(data4))
        myages2 <- tapply(myages, myorder, mean)
        mycase2 <- tapply(mycase, myorder, mean)
        mypairs <- matchAge(myages2, mycase2 , h = 4)  # to get index for computing ratios
        print (paste(a, nrow(myannots), nrow(mypairs)))
      }
      tablecita<- rbind(tablecita,data3)
      tablecitapromediada <-rbind(tablecitapromediada,data5)
      
    }
  }
}


dim(tablecita)
tablecita[,1:30]
dim(tablecitapromediada)
tablecitapromediada[,1:30]
#######
# Negative lipidomics
data<-IAm12Negativelipidomics
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)
new.groups <- unique(annots[,4])
miu<-rownames(tablecitapromediada)
tablecitasub = NULL

for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #for (j in 1:50){
  print(paste(j, "out of", length(new.groups)))
  print(new.groups[j])
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  miuvalue<-as.matrix(tablecitapromediada[j,])
  sustraccion <-sweep(mydata,2,miuvalue[1,])
  #sustraccion2<-cbind (new.groups[j],sustraccion)
  tablecitasub<- rbind(tablecitasub,sustraccion)
}      

NEGLIPXw<-cbind(annots,tablecitasub)
dim(annots)
dim(tablecitasub); dim(NEGLIPXw)
NEGLIPXw[1:100,c(1:4,6,7,29,514)]
#write.table(NEGLIPXw, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NEGLIPXw.txt", sep="\\t")
#NEGLIPXw <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NEGLIPXw.txt")


##############################################################################################################################
# GCTOF
data<-IAm12GCTOF
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)

# Calculate Ratios
new.groups <- unique(annots[,4])
tablecita = NULL
tablecitapromediada=NULL
for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #  for (j in 1:5){
  print(paste(j, "out of", length(new.groups)))
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  for ( a in "IA") {
    if (is.element (a, unique(myannots[,30]))){     # Diagnosis(30)
      myannots2 <- myannots[myannots[,30] == a,]    # Diagnosis(30)
      mydata2 <- mydata[myannots[,30] == a,]        # Diagnosis(30)
      myages <- myannots2[,6]                       # Sconv(6)
      mycase <- myannots2[,29]                      # Outcome(29)
      if (any(mycase==1) & any(mycase == 0)) {
        myorder <- matchAge2(myages, mycase, h = 1)  # to average controls
        data3 <- apply(mydata2,2,function (x) tapply(x, myorder, mean))
        data4 <- data.frame(apply(data3,2,mean))
        colnames(data4)<-new.groups[j]
        data5<-data.frame(t(data4))
        myages2 <- tapply(myages, myorder, mean)
        mycase2 <- tapply(mycase, myorder, mean)
        mypairs <- matchAge(myages2, mycase2 , h = 4)  # to get index for computing ratios
        print (paste(a, nrow(myannots), nrow(mypairs)))
      }
      tablecita<- rbind(tablecita,data3)
      tablecitapromediada <-rbind(tablecitapromediada,data5)
      
    }
  }
}


dim(tablecita)
tablecita[,1:30]
dim(tablecitapromediada)
tablecitapromediada[,1:30]
#######
# Negative lipidomics
data<-IAm12GCTOF
annots <- data[,c(1:35)]  # move annotations
data <- data[,-c(1:35)]  # keep only gene expression
#data <- data.frame(annots,data)
new.groups <- unique(annots[,4])
miu<-rownames(tablecitapromediada)
tablecitasub = NULL

for (j in 1:length(new.groups)){    # the length of the groups is 1358 
  #for (j in 1:50){
  print(paste(j, "out of", length(new.groups)))
  print(new.groups[j])
  myannots <- annots[annots[,4] == new.groups[j],]  # NewGroupsMetabolomics(4)
  mydata <- data[annots[,4] == new.groups[j],]      # NewGroupsMetabolomics(4)
  miuvalue<-as.matrix(tablecitapromediada[j,])
  sustraccion <-sweep(mydata,2,miuvalue[1,])
  #sustraccion2<-cbind (new.groups[j],sustraccion)
  tablecitasub<- rbind(tablecitasub,sustraccion)
}      

GCTOFXw<-cbind(annots,tablecitasub)
dim(annots)
dim(tablecitasub); dim(GCTOFXw)
GCTOFXw[1:100,c(1:4,6,7,29,514)]
#write.table(GCTOFXw, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/GCTOFXw.txt", sep="\\t")
#GCTOFXw <- read.delim("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/GCTOFXw.txt")

###########################################################
###########################################################
#################     Checkpoint 1     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing POSLIPXw, NEGLIPXw and GCTOFXw


colnames(POSLIPXw[1:37])
outcomemetabol<-POSLIPXw[,c(7,29)]  # Mask.ID(7) and Outcome(29)
outcomemetabol2<-POSLIPXw[,c(7,29)]  # Mask.ID(7) and Outcome(29)
head(outcomemetabol)
dim(outcomemetabol)
duplicated(outcomemetabol)
outcomemetabol<-outcomemetabol[!duplicated(outcomemetabol),]
dim(outcomemetabol)   # 1556
duplicated(outcomemetabol2);outcomemetabol2<-outcomemetabol2[!duplicated(outcomemetabol2),];dim(outcomemetabol2)   # 1556


####OutcomeNPLS
rownames(outcomemetabol)<-outcomemetabol[,1]; rownames(outcomemetabol2)<-outcomemetabol2[,1]
head(outcomemetabol);head(outcomemetabol2)
outcomemetabol<-outcomemetabol[,-1]

####
# Array
outcomemetabol<-as.matrix(outcomemetabol)
OutcomeMetabolsdummyarray <- array(data = NA, dim = c(1556,1,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))
dim(OutcomeMetabolsdummyarray)

OutcomeMetabolsdummyarray
OutcomeMetabolsdummyarray[,,1] <- outcomemetabol
OutcomeMetabolsdummyarray[,,2] <- outcomemetabol
OutcomeMetabolsdummyarray[,,3] <- outcomemetabol
OutcomeMetabolsdummyarray[,,4] <- outcomemetabol
OutcomeMetabolsdummyarray[,,5] <- outcomemetabol

rownames(OutcomeMetabolsdummyarray)<-rownames(outcomemetabol2)
colnames(OutcomeMetabolsdummyarray)<-colnames (outcomemetabol2[2])
head(OutcomeMetabolsdummyarray)
#as.data.frame(OutcomeMetabolsdummyarray)
####################################################################################################################
# Making the 3D arrays
individuals<-unique(POSLIPXw$Mask.Id)
length(individuals)
individuals<-data.frame(individuals)
colnames(individuals)[1]<-"Mask.Id"
head(individuals)

####################################################################################################################
# Positive lipidomics
colnames(POSLIPXw[1:40])
dim(POSLIPXw) # 4840 x 549
POSLIP<-POSLIPXw[,c(7,36:549)]
PosLip12<-POSLIP[POSLIPXw$Sconv=="-12",]
dim(PosLip12)
PosLip9<-POSLIP[POSLIPXw$Sconv=="-9",]
PosLip6<-POSLIP[POSLIPXw$Sconv=="-6",]
PosLip3<-POSLIP[POSLIPXw$Sconv=="-3",]
PosLip0<-POSLIP[POSLIPXw$Sconv=="0",]
dim(PosLip12)+dim(PosLip9)+dim(PosLip6)+dim(PosLip3)+dim(PosLip0)
####################################################################################################################
# Merging with all cases
# PosLip
PosLip12total<-merge(individuals,PosLip12, by="Mask.Id", all.x = T);dim(PosLip12total);dim(individuals);dim(PosLip12)
PosLip12total[1:5,1:5]
rownames(PosLip12total)<- PosLip12total[,1]
PosLip12total<- as.matrix(PosLip12total[,c(-1)])
dim(PosLip12total)

PosLip9total<-merge(individuals,PosLip9, by="Mask.Id", all.x = T);dim(PosLip9total);dim(individuals);dim(PosLip9)
PosLip9total[1:5,1:5]
rownames(PosLip9total)<- PosLip9total[,1]
PosLip9total<- as.matrix(PosLip9total[,c(-1)])
dim(PosLip9total)

PosLip6total<-merge(individuals,PosLip6, by="Mask.Id", all.x = T);dim(PosLip6total);dim(individuals);dim(PosLip6)
PosLip6total[1:5,1:5]
rownames(PosLip6total)<- PosLip6total[,1]
PosLip6total<- as.matrix(PosLip6total[,c(-1)])
dim(PosLip6total)

PosLip3total<-merge(individuals,PosLip3, by="Mask.Id", all.x = T);dim(PosLip3total);dim(individuals);dim(PosLip3)
PosLip3total[1:5,1:5]
rownames(PosLip3total)<- PosLip3total[,1]
PosLip3total<- as.matrix(PosLip3total[,c(-1)])
dim(PosLip3total)

PosLip0total<-merge(individuals,PosLip0, by="Mask.Id", all.x = T);dim(PosLip0total);dim(individuals);dim(PosLip0)
PosLip0total[1:5,1:5]
rownames(PosLip0total)<- PosLip0total[,1]
PosLip0total<- as.matrix(PosLip0total[,c(-1)])
dim(PosLip0total)

# Dimensions are 1556*514*5
PosLip0total[1:10,1:4]
#################################################################################
arrayPosLip <- array(data = NA, dim = c(1556,514,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayPosLip
arrayPosLip[,,1] <- PosLip12total
arrayPosLip[,,2] <- PosLip9total
arrayPosLip[,,3] <- PosLip6total
arrayPosLip[,,4] <- PosLip3total
arrayPosLip[,,5] <- PosLip0total

rownames(arrayPosLip)<-rownames(PosLip12total)
colnames(arrayPosLip)<-colnames (PosLip12total)
arrayPosLip
dim(arrayPosLip)   # 1556 * 514 * 5


####################################################################################################################
# Negative lipidomics
colnames(NEGLIPXw[1:40])
dim(NEGLIPXw) # 4840 x 478
NEGLIP<-NEGLIPXw[,c(7,36:478)]

NegLip12<-NEGLIP[NEGLIPXw$Sconv=="-12",]
dim(NegLip12);dim(NEGLIP)
NegLip9<-NEGLIP[NEGLIPXw$Sconv=="-9",]
NegLip6<-NEGLIP[NEGLIPXw$Sconv=="-6",]
NegLip3<-NEGLIP[NEGLIPXw$Sconv=="-3",]
NegLip0<-NEGLIP[NEGLIPXw$Sconv=="0",]
dim(NegLip12)+dim(NegLip9)+dim(NegLip6)+dim(NegLip3)+dim(NegLip0)
####################################################################################################################
# Merging with all cases
# NegLip
NegLip12total<-merge(individuals,NegLip12, by="Mask.Id", all.x = T);dim(NegLip12total);dim(individuals);dim(NegLip12)
NegLip12total[1:5,1:5]
rownames(NegLip12total)<- NegLip12total[,1]
NegLip12total<- as.matrix(NegLip12total[,c(-1)])
dim(NegLip12total)

NegLip9total<-merge(individuals,NegLip9, by="Mask.Id", all.x = T);dim(NegLip9total);dim(individuals);dim(NegLip9)
NegLip9total[1:5,1:5]
rownames(NegLip9total)<- NegLip9total[,1]
NegLip9total<- as.matrix(NegLip9total[,c(-1)])
dim(NegLip9total)

NegLip6total<-merge(individuals,NegLip6, by="Mask.Id", all.x = T);dim(NegLip6total);dim(individuals);dim(NegLip6)
NegLip6total[1:5,1:5]
rownames(NegLip6total)<- NegLip6total[,1]
NegLip6total<- as.matrix(NegLip6total[,c(-1)])
dim(NegLip6total)

NegLip3total<-merge(individuals,NegLip3, by="Mask.Id", all.x = T);dim(NegLip3total);dim(individuals);dim(NegLip3)
NegLip3total[1:5,1:5]
rownames(NegLip3total)<- NegLip3total[,1]
NegLip3total<- as.matrix(NegLip3total[,c(-1)])
dim(NegLip3total)

NegLip0total<-merge(individuals,NegLip0, by="Mask.Id", all.x = T);dim(NegLip0total);dim(individuals);dim(NegLip0)
NegLip0total[1:5,1:5]
rownames(NegLip0total)<- NegLip0total[,1]
NegLip0total<- as.matrix(NegLip0total[,c(-1)])
dim(NegLip0total)

# Dimensions are 1556*443*5
NegLip0total[1:10,1:4]
#################################################################################
arrayNegLip <- array(data = NA, dim = c(1556,443,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayNegLip
arrayNegLip[,,1] <- NegLip12total
arrayNegLip[,,2] <- NegLip9total
arrayNegLip[,,3] <- NegLip6total
arrayNegLip[,,4] <- NegLip3total
arrayNegLip[,,5] <- NegLip0total

rownames(arrayNegLip)<-rownames(NegLip12total)
colnames(arrayNegLip)<-colnames (NegLip12total)
arrayNegLip
dim(arrayNegLip)   # 1556 * 443 * 5



####################################################################################################################
# GCTOF

colnames(GCTOFXw[1:40])
dim(GCTOFXw) # 4840 x 399
GCTOF<-GCTOFXw[,c(7,36:399)]

GCTOF12<-GCTOF[GCTOFXw$Sconv=="-12",]
dim(GCTOF12);dim(GCTOF)
GCTOF9<-GCTOF[GCTOFXw$Sconv=="-9",]
GCTOF6<-GCTOF[GCTOFXw$Sconv=="-6",]
GCTOF3<-GCTOF[GCTOFXw$Sconv=="-3",]
GCTOF0<-GCTOF[GCTOFXw$Sconv=="0",]
dim(GCTOF12)+dim(GCTOF9)+dim(GCTOF6)+dim(GCTOF3)+dim(GCTOF0)
####################################################################################################################
# Merging with all cases
# GCTOF
Grt12total<-merge(individuals,GCTOF12, by="Mask.Id", all.x = T);dim(Grt12total);dim(individuals);dim(GCTOF12)
Grt12total[1:5,1:5]
rownames(Grt12total)<- Grt12total[,1]
Grt12total<- as.matrix(Grt12total[,c(-1)])
dim(Grt12total)

Grt9total<-merge(individuals,GCTOF9, by="Mask.Id", all.x = T);dim(Grt9total)
rownames(Grt9total)<- Grt9total[,1]
Grt9total<- as.matrix(Grt9total[,c(-1)])
dim(Grt9total)

Grt6total<-merge(individuals,GCTOF6, by="Mask.Id", all.x = T);dim(Grt6total)
rownames(Grt6total)<- Grt6total[,1]
Grt6total<- as.matrix(Grt6total[,c(-1)])
dim(Grt6total)

Grt3total<-merge(individuals,GCTOF3, by="Mask.Id", all.x = T);dim(Grt3total)
rownames(Grt3total)<- Grt3total[,1]
Grt3total<- as.matrix(Grt3total[,c(-1)])
dim(Grt3total)

Grt0total<-merge(individuals,GCTOF0, by="Mask.Id", all.x = T);dim(Grt0total)
rownames(Grt0total)<- Grt0total[,1]
Grt0total<- as.matrix(Grt0total[,c(-1)])
dim(Grt0total)
# Dimensions are 1556*364*5
Grt0total[1:10,1:4]
#################################################################################
arrayGCTOF <- array(data = NA, dim = c(1556,364,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))

arrayGCTOF
arrayGCTOF[,,1] <- Grt12total
arrayGCTOF[,,2] <- Grt9total
arrayGCTOF[,,3] <- Grt6total
arrayGCTOF[,,4] <- Grt3total
arrayGCTOF[,,5] <- Grt0total

rownames(arrayGCTOF)<-rownames(Grt12total)
colnames(arrayGCTOF)<-colnames (Grt12total)
arrayGCTOF
dim(arrayGCTOF)   # 1556 * 364 * 5


###########################################################
###########################################################
#################     Checkpoint 2     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing arrayGCTOF, arrayPosLip and arrayNegLip
###########################################################
# subsetting to the size of GE data
dim(arrayGCTOF)

arrayGCTOF308<- arrayGCTOF[rownames(arrayGCTOF) %in% rownames(Outcomedummyarray),,]
dim(arrayGCTOF308)

arrayNegLip308<- arrayGCTOF[rownames(arrayNegLip) %in% rownames(Outcomedummyarray),,]
dim(arrayNegLip308)

arrayPosLip308<- arrayPosLip[rownames(arrayPosLip) %in% rownames(Outcomedummyarray),,]
dim(arrayPosLip308)

###########################################################
# NPLSDA

# Step1: Determining the best fitted model
# GCTOF
modelGCTOF308<-bestfittedmodel (X=arrayGCTOF308,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
#Best model 434

# NegLip
modelNegLip308<-bestfittedmodel (X=arrayNegLip308,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by TimemodelGCTOFsizeGE<-bestfittedmodel (X=arrayGCTOFsizeGE,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
#Best model 434

# PosLip
modelPosLip308<-bestfittedmodel (X=arrayPosLip308,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by TimemodelGCTOFsizeGE<-bestfittedmodel (X=arrayGCTOFsizeGE,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
#Best model 423

# Step2: Imputing the best fitted model data
# GCTOF
fullarrayGCTOFsizeGE<-Imputemethod(X=arrayGCTOF308,fac=c(4, 3, 4), conver = 1e-07, max.iter = 1000)
summary(fullarrayGCTOFsizeGE)
dim(fullarrayGCTOFsizeGE)

# NegLip
fullarrayNegLipsizeGE<-Imputemethod(X=arrayNegLip308,fac=c(4, 3, 4), conver = 1e-07, max.iter = 1000)
summary(fullarrayNegLipsizeGE)
dim(fullarrayNegLipsizeGE)

# PosLip
fullarrayPosLipsizeGE<-Imputemethod(X=arrayPosLip308,fac=c(4, 2, 3), conver = 1e-07, max.iter = 1000)
summary(fullarrayPosLipsizeGE)
dim(fullarrayPosLipsizeGE)


###########################################################
###########################################################
#################     Checkpoint 3     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing fullarrayPosLipsizeGE, fullarrayNegLipsizeGE and fullarrayGCTOFsizeGE
###########################################################
# Setp3: NPLSDA
NPLSDAGCTOFsizeGEBestModel<-NPLSDAmod(XN=fullarrayGCTOFsizeGE, YN=Outcomedummyarray, outcome.Y=NULL, factors=3, centering=2) 
summary(NPLSDAGCTOFsizeGEBestModel)

NPLSDANegLipsizeGEBestModel<-NPLSDAmod(XN=fullarrayNegLipsizeGE, YN=Outcomedummyarray, outcome.Y=NULL, factors=3, centering=2) 
summary(NPLSDANegLipsizeGEBestModel)

NPLSDAPosLipsizeGEBestModel<-NPLSDAmod(XN=fullarrayPosLipsizeGE, YN=Outcomedummyarray, outcome.Y=NULL, factors=3, centering=2) 
summary(NPLSDAPosLipsizeGEBestModel)

#Step4: Plotting
ploteoNPLSDAFullarrayGCTOFsizeGE<- plotNPLSDAmod (X=NPLSDAGCTOFsizeGEBestModel, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                  cutoff = 20, factors=3, penalty=1) 

ploteoNPLSDAFullarrayNegLipsizeGE<- plotNPLSDAmod (X=NPLSDANegLipsizeGEBestModel, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                  cutoff = 20, factors=3, penalty=1) 

ploteoNPLSDAFullarrayPosLipsizeGE<- plotNPLSDAmod (X=NPLSDAPosLipsizeGEBestModel, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                  cutoff = 20, factors=3, penalty=1) 


###########################################################
###############     Variable Selection     ################
###########################################################
# Positive Lipidomics # Based on the results we are going to select variables by the most varying
# in each of the three components separatelly
dim(arrayGCTOF);dim(arrayNegLip);dim(arrayPosLip)

arrayGCTOF136<-arrayGCTOF[is.element(rownames(arrayGCTOF), patients3tps[,1]),,]
dim(arrayGCTOF136)   # 136 x 364 x 5 

arrayNegLip136<-arrayNegLip[is.element(rownames(arrayNegLip), patients3tps[,1]),,]
dim(arrayNegLip136)   # 136 x 443 x 5 

arrayPosLip136<-arrayPosLip[is.element(rownames(arrayPosLip), patients3tps[,1]),,]
dim(arrayPosLip136)   # 136 x 514 x 5 

###########################################################
# GCTOF
modelsGCTOF136<-bestfittedmodel(arrayGCTOF136, centering = 0)   # Best model: 424

fullarrayGCTOF136<-Imputemethod(X=arrayGCTOF136, fac=c(4, 2, 4), conver = 1e-07, max.iter = 1000)

# NegLip
modelsNegLip136<-bestfittedmodel(arrayNegLip136, centering = 0)   # Best model: 424 also

fullarrayNegLip136<-Imputemethod(X=arrayNegLip136, fac=c(4, 2, 4), conver = 1e-07, max.iter = 1000)

# PosLip
modelsPosLip136<-bestfittedmodel(arrayPosLip136, centering = 0)   # Best model: 444

fullarrayPosLip136<-Imputemethod(X=arrayPosLip136, fac=c(4, 4, 4), conver = 1e-07, max.iter = 1000)

###########################################################
###########################################################
#################     Checkpoint 4     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing fullarrayPosLip136, fullarrayNegLip136 and fullarrayGCTOF136
###########################################################
#   GCTOF:
NPLSDA136GCTOF<-NPLSDAmod (XN=fullarrayGCTOF136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                           centering=0)

ploteoNPLSDAGCTOFtotal136indvs<- plotNPLSDAmod (X=NPLSDA136GCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=1) 


###########################################################
#   NegLip:
NPLSDA136NegLip<-NPLSDAmod (XN=fullarrayNegLip136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                            centering=0)

ploteoNPLSDANegLiptotal136indvs<- plotNPLSDAmod (X=NPLSDA136NegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

NPLSDA136NegLip$NPLSDAQ22D

round(NPLSDA136NegLip$NPLSDAQ22D[1,6], digits=3)
###########################################################
#   PosLip:
NPLSDA136PosLip<-NPLSDAmod (XN=fullarrayPosLip136, YN=outcomedummyarray136,factors = 2,COMP= c(4,2,4), conver = 1e-16, max.iteration = 10000,
                            centering=0)

ploteoNPLSDAPosLiptotal136indvs<- plotNPLSDAmod (X=NPLSDA136PosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                 cutoff = 20, factors=2, penalty=1) 

###########################################################
############    Variable Selection GCTOF     ##############
###########################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136GCTOF$VIP3Dmodel2)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayGCTOF136)
PosLipselVars
#selectedGCTOFbydifVIPs<-list(VIP3Dmodel2.95p.85v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 3), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP3Dmodel1
dim(NPLSDA136GCTOF$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136GCTOF$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayGCTOF136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=3, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=3, penalty=2) 
#dev.off()


#summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.1comp.213g=PosLipselVars))
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.12comp.336g=PosLipselVars))


######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136GCTOF)
NPLSDA136GCTOF$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136GCTOF$VIP2D)

colp1<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136GCTOF$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayGCTOF136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                        genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedGCTOFbydifVIPs)
#selectedGCTOFbydifVIPs<-append(selectedGCTOFbydifVIPs,list(VIP2D.95p.comp12.82v= as.vector(genes)))
#save(selectedGCTOFbydifVIPs,file="selectedGCTOFbydifVIPs.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedGCTOFbydifVIPs.RData")
###########################################################################################################
# Venn Diagram of the four genelists #

summary(selectedGCTOFbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedGCTOFbydifVIPs, removeEMPTYstrings)

summary(geneLS2)



lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")

VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDAGCTOFLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

#pdf("VenndiagramsselectedGCTOF.pdf")
grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunion<- union(selectedGCTOFbydifVIPs[[1]],selectedGCTOFbydifVIPs[[2]])
length(theunion) #91
theintersection<- intersect(selectedGCTOFbydifVIPs[[1]],selectedGCTOFbydifVIPs[[2]])
length(theintersection) # 76

####
FullarrayUNIONGCTOFSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                          theunion),]
dim(FullarrayUNIONGCTOFSelVars)

### Now NPLSDA and Graph
NPLSDAUnionGCTOF<-NPLSDAmod(XN=FullarrayUNIONGCTOFSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONGCTOF.pdf")
ploteoNPLSDAUnionGCTOF<- plotNPLSDAmod (X=NPLSDAUnionGCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                        cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTGCTOFSelVars<-fullarrayGCTOF136[,is.element(colnames(fullarrayGCTOF136),
                                                              theintersection),]
dim(FullarrayINTERSECTGCTOFSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionGCTOF<-NPLSDAmod(XN=FullarrayINTERSECTGCTOFSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionGCTOF.pdf")
ploteoNPLSDAIntersectionGCTOF<- plotNPLSDAmod (X=NPLSDAIntersectionGCTOF, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                               cutoff = 20, factors=2, penalty=2) 
#dev.off()


#############################################################################################################
###################################    Variable Selection NegLip     ########################################
#############################################################################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136NegLip$VIP3Dmodel2)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayNegLip136)
PosLipselVars
#selectedNegLipbydifVIPs<-list(VIP3Dmodel2.95p.98v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP3Dmodel1
dim(NPLSDA136NegLip$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136NegLip$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayNegLip136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()


summary(selectedgenesbydifVIPs)
#selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.1comp.213g=PosLipselVars))
selectedgenesbydifVIPs<-append(selectedgenesbydifVIPs,list(VIP3Dmodel1.99p.12comp.336g=PosLipselVars))


######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136NegLip)
NPLSDA136NegLip$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136NegLip$VIP2D)

colp1<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136NegLip$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayNegLip136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                         genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedNegLipbydifVIPs)
#selectedNegLipbydifVIPs<-append(selectedNegLipbydifVIPs,list(VIP2D.95p.comp12.153v= as.vector(genes)))
#save(selectedNegLipbydifVIPs,file="selectedNegLipbydifVIPs.RData")
###########################################################################################################
# We selected VIP3D/2 95% and VIP2D 95% Comps 1&2
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedNegLipbydifVIPs.RData")

# Venn Diagram of the four genelists #
summary(selectedNegLipbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedNegLipbydifVIPs, removeEMPTYstrings)

summary(geneLS2)



lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:
require("VennDiagram")

VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDANegLipLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
plot.new()
grid.draw(venn.plot)

#pdf("VenndiagramsselectedNegLip.pdf")
#grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunionNegLip<- union(selectedNegLipbydifVIPs[[1]],selectedNegLipbydifVIPs[[2]])
length(theunionNegLip) #160
theintersectionNegLip<- intersect(selectedNegLipbydifVIPs[[1]],selectedNegLipbydifVIPs[[2]])
length(theintersectionNegLip) # 91

####
FullarrayUNIONNegLipSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                            theunionNegLip),]
dim(FullarrayUNIONNegLipSelVars)

### Now NPLSDA and Graph
NPLSDAUnionNegLip<-NPLSDAmod(XN=FullarrayUNIONNegLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONNegLip.pdf")
ploteoNPLSDAUnionNegLip<- plotNPLSDAmod (X=NPLSDAUnionNegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                         cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTNegLipSelVars<-fullarrayNegLip136[,is.element(colnames(fullarrayNegLip136),
                                                                theintersectionNegLip),]
dim(FullarrayINTERSECTNegLipSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionNegLip<-NPLSDAmod(XN=FullarrayINTERSECTNegLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionNegLip.pdf")
ploteoNPLSDAIntersectionNegLip<- plotNPLSDAmod (X=NPLSDAIntersectionNegLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=2) 
#dev.off()
#summary(selectedNegLipbydifVIPs)
#save(selectedNegLipbydifVIPs, file="selectedNegLipbydifVIPs.RData")
################################################################################################################

#############################################################################################################
###################################    Variable Selection PosLip     ########################################
#############################################################################################################
# Variable selection by VIP3Dmodel2
summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP3Dmodel2

### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136PosLip$VIP3Dmodel2)

thrs<-99
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:length(vipsoutcomemet)){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
######################################################################################################
# Retain just these variables in GCTOF 
dim(fullarrayPosLip136)
PosLipselVars
#selectedPosLipbydifVIPs<-list(VIP3Dmodel2.99p.28v=PosLipselVars)

#summary(selectedGCTOFbydifVIPs)
#selectedgenesbydifVIPs<-c(selectedgenesbydifVIPs,tulo="pelotas")

FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

#FullarrayGenesVIPSelVarsVIP3D.2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP3D.2,file="FullarrayGenesVIPSelVarsVIP3D.2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

#NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3
#NPLSDAarraygenesVIPselected$NPLSDAvariates$X
#dim(NPLSDAarraygenesVIPselected$NPLSDAvariatesperMode3)
# Plotting
#pdf("VIP3Dmodel2.99p.1013g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

ploteoNPLSDAVIPselectedgenes

######################################################################################################################
# Variable selection by VIP3Dmodel1

summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP3Dmodel1
dim(NPLSDA136PosLip$VIP3Dmodel1)



### Desde aca
vipsoutcomemet<-data.frame(NPLSDA136PosLip$VIP3Dmodel1)

thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)  # 2050
######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayPosLip136)
PosLipselVars

FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         PosLipselVars),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP1.comp1<-FullarrayGenesVIPSelVars
#FullarrayGenesVIPSelVarsVIP1.comps12<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP1.comp1,file="FullarrayGenesVIPSelVarsVIP1.comp1.RData")
#save(FullarrayGenesVIPSelVarsVIP1.comps12,file="FullarrayGenesVIPSelVarsVIP1.comps12.RData")


### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf("VIP3Dmodel1.99p.1comp.213g.pdf")
#pdf("VIP3Dmodel1.99p.12comp.336g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()


#summary(selectedPosLipbydifVIPs)
#selectedPosLipbydifVIPs<-append(selectedPosLipbydifVIPs,list(VIP3Dmodel1.95p.comp12.46v=PosLipselVars))

#save(selectedPosLipbydifVIPs, file="selectedPosLipbydifVIPs.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/selectedPosLipbydifVIPs.RData")
######################################################################################################
# Gene Expression USING VIP2D

summary(NPLSDA136PosLip)
NPLSDA136PosLip$VIP2D


### Desde aca
vipsoutcome2D<-data.frame(NPLSDA136PosLip$VIP2D)

colp1<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-12", sep="!_____")
colp2<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-9", sep="!_____")
colp3<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-6", sep="!_____")
colp4<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"-3", sep="!_____")
colp5<-paste(rownames(NPLSDA136PosLip$VIP3Dmodel1),"0", sep="!_____")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

#apply(vipsoutcomemet, 2, function(x) is.numeric(x))
vipsoutcomemet<-vipsoutcome2D
thrs<-95
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
for (i in 1:2){
  a<-quantile(vipsoutcomemet[,i], probs = c(50, 95,99, thrs)/100)
  print(i)
  print(a)
  temp<-vipsoutcomemet[vipsoutcomemet[,i]>=a[4],]
  namesubset<-paste('t',i,sep='')
  vips<-list(temp)
  vipsoutcomemetsubset[[namesubset]]<-vips
}
vipsoutcomemetsubset
rm(vaya)
#rownames(vipsoutcomemetsubset$t1[[1]])
vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"!_____")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

######################################################################################################
# Retain just these variables in Gene Expression 2D
dim(fullarrayPosLip136)
genes
#as.vector(genes)
FullarrayGenesVIPSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                         genes),]
dim(FullarrayGenesVIPSelVars)
#FullarrayGenesVIPSelVarsVIP2Dcomp2<-FullarrayGenesVIPSelVars
#save(FullarrayGenesVIPSelVarsVIP2Dcomp2,file="FullarrayGenesVIPSelVarsVIP2Dcomp2.RData")

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("VIP2D.99p.comp2.1022g.pdf")
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 
#dev.off()

#summary(selectedNegLipbydifVIPs)
#selectedNegLipbydifVIPs<-append(selectedNegLipbydifVIPs,list(VIP2D.95p.comp12.153v= as.vector(genes)))
#save(selectedNegLipbydifVIPs,file="selectedNegLipbydifVIPs.RData")
###########################################################################################################
# We selected for PosLip: VIP3D/2 99% and VIP3D/1 95%

# Venn Diagram of the four genelists #

summary(selectedPosLipbydifVIPs)

# If this is a bit confusing you can also write a function and then use it in lapply() 
removeEMPTYstrings <- function(x) {
  
  newVectorWOstrings <- x[x != ""]
  return(newVectorWOstrings)
  
}
geneLS2 <- lapply(selectedPosLipbydifVIPs, removeEMPTYstrings)

summary(geneLS2)
lapply(geneLS2, tail) # Both methods return the same results

# Now we can plot a Venn diagram with the VennDiagram R package, as follows:


VENN.LIST <- geneLS2
venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"),#,"darkseagreen4","goldenrod2"), 
                          alpha=c(0.5,0.5), cex = 2, cat.fontface=4, 
                          category.names=c("VIP3D.2", "VIP2D.comp12"), main="NPLSDAPosLipLists",
                          cex.main=2
)

# To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
grid.draw(venn.plot)

#pdf("VenndiagramsselectedPosLip.pdf")
#grid.draw(venn.plot)
#dev.off()
####################################################################################################################
# Test the Union and the intersection between the two metabolites sets

theunionPosLip<- union(selectedPosLipbydifVIPs[[1]],selectedPosLipbydifVIPs[[2]])
length(theunionPosLip) #63
theintersectionPosLip<- intersect(selectedPosLipbydifVIPs[[1]],selectedPosLipbydifVIPs[[2]])
length(theintersectionPosLip) # 11

####
FullarrayUNIONPosLipSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                            theunionPosLip),]
dim(FullarrayUNIONPosLipSelVars)

### Now NPLSDA and Graph
NPLSDAUnionPosLip<-NPLSDAmod(XN=FullarrayUNIONPosLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("UNIONPosLip.pdf")
ploteoNPLSDAUnionPosLip<- plotNPLSDAmod (X=NPLSDAUnionPosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                         cutoff = 20, factors=2, penalty=2) 
#dev.off()

################################################################################################################
FullarrayINTERSECTPosLipSelVars<-fullarrayPosLip136[,is.element(colnames(fullarrayPosLip136),
                                                                theintersectionPosLip),]
dim(FullarrayINTERSECTPosLipSelVars)

### Now NPLSDA and Graph
NPLSDAIntersectionPosLip<-NPLSDAmod(XN=FullarrayINTERSECTPosLipSelVars, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
#pdf ("IntersectionPosLip.pdf")
ploteoNPLSDAIntersectionPosLip<- plotNPLSDAmod (X=NPLSDAIntersectionPosLip, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=2) 
#dev.off()





#####################################################
######     Merge of all Metabolomics Data     #######
#####################################################
#Merge winners

###########################################################
###########################################################
#################     Checkpoint 5     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing FullarrayUNIONGCTOFSelVars, 
# FullarrayINTERSECTNegLipSelVars and FullarrayUNIONPosLipSelVars
###########################################################
dim(FullarrayUNIONGCTOFSelVars)          # 91
#save(FullarrayUNIONGCTOFSelVars, file= "FullarrayUNIONGCTOFSelVars.RData")
dim(FullarrayINTERSECTNegLipSelVars)     # 91
sort(table(colnames(FullarrayINTERSECTNegLipSelVars)))
#save(FullarrayINTERSECTNegLipSelVars, file= "FullarrayINTERSECTNegLipSelVars.RData")
dim(FullarrayUNIONPosLipSelVars)         # 63
sort(table(colnames(FullarrayUNIONPosLipSelVars)))
#save(FullarrayUNIONPosLipSelVars, file="FullarrayUNIONPosLipSelVars.RData")

###########################################################
# There are repeated names in PosLip and NegLip, so we are adding a little pieace of name
FullarrayINTERSECTNegLipSelVars
colp1<-paste(colnames(FullarrayINTERSECTNegLipSelVars),"NegLip", sep="!_")
colnames(FullarrayINTERSECTNegLipSelVars)<-colp1

FullarrayUNIONPosLipSelVars
colp2<-paste(colnames(FullarrayUNIONPosLipSelVars),"PosLip", sep="!_")
colnames(FullarrayUNIONPosLipSelVars)<-colp2


posneglip<-abind(FullarrayUNIONPosLipSelVars,FullarrayINTERSECTNegLipSelVars,along=2);dim(FullarrayUNIONPosLipSelVars);dim(FullarrayINTERSECTNegLipSelVars);dim(posneglip)
Allmetabolomics136<-abind(posneglip,FullarrayUNIONGCTOFSelVars,along=2);dim(posneglip);dim(FullarrayUNIONGCTOFSelVars);dim(Allmetabolomics136)
#
#save(Allmetabolomics136,file="Allmetabolomics136.RData")
DefinitiveMetabolomicsListFeb2<-as.data.frame(colnames(Allmetabolomics136))
#write.table(DefinitiveMetabolomicsListFeb2, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/DefinitiveMetabolomicsListFeb2.txt", sep="\\t")
#DefinitiveMetabolomicsListFeb2 <- read.csv("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/DefinitiveMetabolomicsListFeb2.txt", sep="")

###########################################################
# NPLSDA and Plotting
dim(outcomedummyarray136)
dim(Allmetabolomics136)
#colnames(Allmetabolomics136)
#sort(table(colnames(Allmetabolomics136)))

NPLSDAallmetabolomicsSelVars<-NPLSDAmod(XN=Allmetabolomics136, YN=outcomedummyarray136, outcome.Y=NULL, factors=2, centering=0) 


#save(NPLSDAallmetabolomicsSelVars,file="NPLSDAallmetabolomicsSelVars.RData")
# Plotting
#pdf ("NPLSDAallmetabolomicsSelVars.pdf")
ploteo<- plotNPLSDAmod (X=NPLSDAallmetabolomicsSelVars, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                        cutoff = 20, factors=2, penalty=2)
#dev.off()

###########################################################