###########################################################
########   VIP_NPLSDADietaryBiomarkersxTEDDY.R     ########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the Dietary Biomarkers dataset to calculate the NPLSDA and VIP selection

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()
###########################################################
# Functions:
source ("/media/data/leobalzano/ScriptsForTEDDY/Scripts/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

###########################################################
# Data:
# Dietary biomarkers
load("/media/data/leobalzano/ScriptsForTEDDY/Data/DSIA.RData")
vitC <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/vitC.csv")
fatty <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/mp149_fatty_acids.csv")
toco <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/mp149_tocopherols.csv")
metabvitD <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/mp149_vit_d.csv")
metabCarot <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/mp149_carotenoids.csv")

# Checkpoint 1:
load("/media/data/leobalzano/ScriptsForTEDDY/Data/metabyDS12m.RData")

# Checkpoint 2: 
DietaryBiomarkersOctoberXw <- read.delim ("/media/data/leobalzano/ScriptsForTEDDY/Data/DietaryBiomarkersOctoberXw.txt", sep= "\\t")

# Checkpoint 3: 
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayDB.RData")

# Response Variable
load("/media/data/leobalzano/ScriptsForTEDDY/Data/Outcomedummyarray.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/outcomedummyarray136.RData")

# Misc
DEMOGRAPHIC160424 <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/DEMOGRAPHIC160424.csv")
HLA_DIAGNOSIS20160418 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/HLA_DIAGNOSIS20160418")


# List of Cases with at least 3 out of 5 time points with data
patients3tps <- read.table("/media/data/leobalzano/ScriptsForTEDDY/Data/patients3tps.txt", quote="\\"", comment.char="")

###########################################################
# Libraries:
library("outliers")
require("abind")
library(data.table)

###########################################################
########     Merge with demographics data     #############
dim(DEMOGRAPHIC160424)

# Merge
dim(DSIA)
DemoSia<-merge (DEMOGRAPHIC160424,DSIA, by="Mask.Id")
dim(DemoSia)

###########################################################
##########     Merging with HLA Diagnosis     #############
HLA_DIAGNOSIS20160418[1:10,1:10]
dim(HLA_DIAGNOSIS20160418)
dim(DemoSia)
HLADemoSia<- merge (HLA_DIAGNOSIS20160418,DemoSia, by ="Mask.Id")
dim(HLADemoSia)

###########################################################
#######     Merging with metabolic Ascorbic Acid     ######

dim(vitC)
head(vitC)
vitC<-vitC[with(vitC, order(vitC[,2],vitC[,1])),]  #Ordered by Mask.Id, Agemonth
vitC<-vitC[,c(2,1,4)]
colnames(vitC)<- c("Mask.Id", "Agemonth", "metabvitC")
vitC<-as.data.frame(vitC)

HLADemoSia2<- merge (vitC,HLADemoSia, by =c("Mask.Id", "Agemonth"), all=T)
dim(HLADemoSia)
dim(HLADemoSia2)
HLADemoSia2[20:40,1:30]

sum(is.na(HLADemoSia2$Mask.Id))
sum(is.na(HLADemoSia2$Agemonth))
as.factor(HLADemoSia2$Agemonth)

###########################################################
#######     Merging with metabolic Fatty Acids     ########

head(fatty)
dim(fatty)
dim(HLADemoSia2)

fatty[1:10,]
colnames(fatty)
# Erasing unnecessary columns: sumoffattyacids, unit
fatty2<-fatty[,-c(25,27)]
fatty2[1:10,]
dim(fatty2)
colnames(fatty2)
# Ordering the columns, first Mask.ID then Agemonth, then the others
fatty3<-fatty2[,c(24,10, 1:9,11:23,25:27)]
fatty3[1:10,]
#Change names
colnames(fatty3)[1:2]<- c("Mask.Id", "Agemonth")
#Order rows 
fatty4<-fatty3[with(fatty3, order(fatty3[,1],fatty3[,2])),]  #Ordered by Mask.Id, Agemonth
fatty4[1:10,]

dim(fatty4)
dim(HLADemoSia2)

HLADemoSia3<- merge (fatty4,HLADemoSia2, by =c("Mask.Id", "Agemonth"), all=TRUE)
dim(HLADemoSia2)
dim(HLADemoSia3)
HLADemoSia3[1:20,1:30]

sum(is.na(HLADemoSia3$Mask.Id)) # has to be 0
sum(is.na(HLADemoSia3$Agemonth)) # has to be 0
as.factor(HLADemoSia3$Agemonth) # from 3 to 90

###########################################################
########     Merging with metabolic Tocopherols     #######

dim(toco)
toco[1:10,]

# Erasing unnecessary columns: unit
toco2<-toco[,-4]
toco2[1:10,]
dim(toco2)
# Ordering the columns, first Mask.ID then Agemonth, then the others
toco3<-toco2[,c(3,1,4,2)]
toco3[1:10,]
#Change names
colnames(toco3)[1:2]<- c("Mask.Id", "Agemonth")
#Order rows 
toco4<-toco3[with(toco3, order(toco3[,1],toco3[,2])),]  #Ordered by Mask.Id, Agemonth
toco4[1:10,]


dim(toco4)
dim(HLADemoSia3)

HLADemoSia4<- merge (toco4,HLADemoSia3, by =c("Mask.Id", "Agemonth"), all=TRUE)
dim(HLADemoSia3);dim(HLADemoSia4)
HLADemoSia4[1:20,1:30]

sum(is.na(HLADemoSia4$Mask.Id))
sum(is.na(HLADemoSia4$Agemonth))    ####Aca estaba la vaina habia antes 3 NA... los cambiaron
dim(HLADemoSia4)  # before 29061x109   rows, now  29058x108    rows
#HLADemoSia4<-HLADemoSia4[!(is.na(HLADemoSia4$Agemonth)),]

###########################################################
#######     Merging with metabolic Vitamin D     ##########

dim(metabvitD)
metabvitD[1:100,]

sum(is.na(metabvitD$num_result))
sum(is.na(metabvitD$result))
class(metabvitD$num_result)
class(metabvitD$result)
as.factor(metabvitD$num_result)
as.factor(metabvitD$result)

metabvitD2<-metabvitD[,c(-3,-5)]
metabvitD2[1:10,]
dim(metabvitD2)
#Change names
colnames(metabvitD2)[1:3]<- c("Mask.Id", "Agemonth", "metabvitD")
#Order rows 
metabvitD3<-metabvitD2[with(metabvitD2, order(metabvitD2[,1],metabvitD2[,2])),]  #Ordered by Mask.Id, Agemonth
metabvitD3[1:10,]

dim(metabvitD3)
dim(HLADemoSia4)

HLADemoSia5<- merge (metabvitD3,HLADemoSia4, by =c("Mask.Id", "Agemonth"), all=TRUE)
dim(HLADemoSia4)
dim(HLADemoSia5)
HLADemoSia5[1:20,1:30]


sum(is.na(HLADemoSia5$Mask.Id))
sum(is.na(HLADemoSia5$Agemonth))

###########################################################
#######     Merging with metabolic Carotenoids     ########

dim(metabCarot)
summary(metabCarot)
metabCarot[1:10,4:5]

# Erasing unnecessary column: unit
metabCarot2<-metabCarot[,-5]
metabCarot2[1:10,]
dim(metabCarot2)
#Change mask_id name
colnames(metabCarot2)[4]<- "Mask.Id"
#Change Agemonth name
colnames(metabCarot2)[16]<- "Agemonth"
#Order columns 
metabCarot3<-metabCarot2[,c(4,16,1:3,5:15)]
metabCarot3[1:10,]
#Change the other names
colnames(metabCarot3)[10:16]<- c("Carotc1", "Carotc2", "Carotc4","Carotc5","Carotc6","Carotc7","Carotc8b")
#Order rows 
metabCarot4<-metabCarot3[with(metabCarot3, order(metabCarot3[,1],metabCarot3[,2])),]  #Ordered by Mask.Id, Agemonth
metabCarot4[1:10,]

dim(metabCarot4)
dim(HLADemoSia5)
summary(HLADemoSia5)

HLADemoSia6<- merge (metabCarot4,HLADemoSia5, by =c("Mask.Id", "Agemonth"), all=TRUE)
dim(HLADemoSia5)
dim(HLADemoSia6)
HLADemoSia5[1:20,1:30]
HLADemoSia6[1:20,1:30]

sum(is.na(HLADemoSia5$Mask.Id))
sum(is.na(HLADemoSia5$Agemonth))

DSMetab<-HLADemoSia6

# Ordering the columns, first Mask.ID then Agemonth, then all descriptive ones and then the others
colnames (DSMetab)
dim(DSMetab)   # 29230x123
DSMetab<-DSMetab[,c(1:2,46:68,3:45,69:123)]
DSMetab[1:10,]

sum(is.na(DSMetab$Mask.Id))
sum(is.na(DSMetab$Agemonth))    #antes 1 Agemonth que es NA ahora 0.

dim(DSMetab)  # before 29230x122   rows, now  29218x124    rows
sum(is.na(DSMetab$metabvitD))
sum(!is.na(DSMetab$metabvitD))

###########################################################
# Eliminating individuals that do not belongs to IA (so T1D)
sum(is.na(DSMetab$Diagnosis))
as.factor(DSMetab$Diagnosis)
dim(DSMetab) # before eliminating =29218, after eliminating the T1D= 28407
DSMetab<-DSMetab[!(is.na(DSMetab$Diagnosis)),]

###########################################################

DSMetab[1:54,1:30]
dim(DSMetab)
DSMetab[1:5,c(2,23:25)]

####1) Converting Case.Endptage to monthly, dividing by 30
ata<-DSMetab
summary(DSMetab)
ata[,23] <- round(ata[,23]/30)
ata[1:30,c(2,23:26)]

colnames(ata)

####2)Order table
ata<-ata[with(ata, order(-ata[,24],ata[,25],ata[,2])),]  #Ordered by outcome, Diagnosis and agemonth
head(ata)
dim (ata)

ata[ata$Mask.Id==754319,1:60]
#ata[1:20,1:50]
#toco4[toco4$Mask.Id==754319,]
#toco[toco$mask_id==754319,]
#toco[1:10,]

####3)Create column that is the substraction of Agemonth - enpointStage, now i have the 0;-3 and -6
#Seroconversion column
sum(is.na(ata$Agemonth))   #To check for any error
sum(is.na(ata$Case.Endptage)) #822
sum(is.na(ata$Diagnosis))
colnames(ata)
Sconv<- ata[,2]-ata[,23]
Sconv
ata2<-cbind(Sconv,ata)
ata2[1:5,1:30]
dim(ata2)     # 28407  x  124
sum(is.na(ata2$Sconv))
sum(is.na(ata2$Agemonth))
sum(is.na(ata2$Case.Endptage))   # 811
dim(ata2)    # 28407 x 124     
table(as.factor(ata2$Sconv))
table(as.factor(is.na(ata2[27:69])))   # 1.047.763 NAs in Dietary Biomarker data


####5)Eliminate just cases and controls with Sconv <-13, 
ata4<-ata2[!(ata2$Sconv<(-13)),]
dim(ata4)  #23888 rows  
ata4[1:10,1:10]
colnames(ata4)
table(as.factor(ata4$Sconv))
table(as.factor(is.na(ata4[27:69])))   # 930.335 NAs in Dietary Biomarker data

####6)Eliminate cases and controls with Sconv >0, 
ata5<-ata4[!(ata4$Sconv>(1)),]
#ata5[1:50,1:5]
dim(ata5)  #4845 rows
table(as.factor(ata5$Sconv))
table(as.factor(is.na(ata5[27:69])))   # 123.261 NAs in Dietary Biomarker data
sum(is.na(ata5$Sconv))  #Just to check for errors
sum(is.na(ata5$Agemonth))   #Just to check for errors
sum(is.na(ata5$Outcome))   #Just to check for errors
sum(is.na(ata5$Mask.Id))   #Just to check for errors
sum(is.na(ata5$Case.Endptage))   #Just to check for errors


####9) Regroup data in the 3 time points
ata7<-ata5
ata7$Sconv[ata7$Sconv==-13] <- -12
ata7$Sconv[ata7$Sconv==-11] <- -12

ata7$Sconv[ata7$Sconv==-10] <- -9
ata7$Sconv[ata7$Sconv==-8] <- -9

ata7$Sconv[ata7$Sconv==-7] <- -6
ata7$Sconv[ata7$Sconv==-5] <- -6

ata7$Sconv[ata7$Sconv==-4] <- -3
ata7$Sconv[ata7$Sconv==-2] <- -3

ata7$Sconv[ata7$Sconv==-1] <- 0
ata7$Sconv[ata7$Sconv== 1] <- 0
table(as.factor(ata7$Sconv))
summary(ata7)
###########################################################
###  Creation of AgeGroups which is the variable that groups individuals in terms of age to try to determine if there is
###  any difference of disease appearance due to the age of the kid
###  We decide to use the following frames: Year1<=13,Year2<=25,Year3<=37,Year4>37

dim(ata7)  #4845 for 13 months
#View(ata7)
#### Creation of the variable AgeGroup
AgeGroup = NULL

for(i in 1:nrow(ata7)) {
  
  grupin = NULL
  if(ata7$Case.Endptage[i]>0 & ata7$Case.Endptage[i]<=13){
    grupin = 1
  }
  if(ata7$Case.Endptage[i]>13 & ata7$Case.Endptage[i]<=25){ 
    grupin = 2
  }
  if(ata7$Case.Endptage[i]>25 & ata7$Case.Endptage[i]<=37){ 
    grupin = 3
  }
  if(ata7$Case.Endptage[i]>37){ 
    grupin = 4
  }
  AgeGroup=c(AgeGroup,grupin)
}

dim(AgeGroup)
as.vector(AgeGroup)

ata8<-cbind(AgeGroup,ata7)
ata8[1:200,c(1:10)]
dim(ata8)  #4845 for 13 months

###How many NA do we have and how much this amount of missing values they represent
ata8[1:10,11:28]
colnames (ata8)
sum(is.na(ata8[,28:124]))     #Total  antes de cambiar 175346 despues 172804
sum(!is.na(ata8[,28:124]))     #Total antes de cambiar 294619 despues 297161
172804 + 297161      #Total Total=262990
#  for 13 months 469965 (antes igual que ahora, lo que indica que cambiaron los valores de las celdas pero no el numero de ellas)

# Antes
17534600/469965   ## 37.31 % of missing values for 13 months 
#Ahora
17280400/469965   ## 36.8 % of missing values, (mejoro un pelin)

metabyDS12m<-ata8

###########################################################
###########################################################
#################     Checkpoint 1     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing metabyDS12m

colnames(metabyDS12m)
DB12m<-metabyDS12m[,1:70]
colnames(DB12m)

##################################################################################################################
# Create just one group set, in this case due to the variations of Agemonth, the best way to do it is by
# using CASE.IND (24) and Sconv (2) nevertheless it is not like previous cases, this is very variable
##################################################################################################################
colnames(DB12m)
ata7<-DB12m
colnames(ata7)
ata7[1:3,c(1,2,22,24:27,29,30,31,32)]

# First we have to order the data by CASE.IND (24), Agemonth (4) and -Outcome (26)
ata7<-ata7[with(ata7, order(ata7[,24], ata7[,4], -ata7[,26])),]  #Ordered by CASE.IND(24), Agemonth(4) and -Outcome (26)
ata7[1:10,c(1,2,4,24:27,28,30,31,32)]

# Second Create the groups
gruposDB <-apply(ata7[,c(24,2)], 1, paste, collapse ="_")   # CASE.IND (24) and Sconv (2)
head(gruposDB)
unigroupDB <- unique (gruposDB)
length(unigroupDB)   # 1212 groups

datos <- cbind (gruposDB, ata7)
datos[1:50,c(1,6:7,23,25,30,31)]   #datos(1212 rows)Agemonth (30)

NewGroupsDB=NULL
n=0
Nn ="word"
grupin=NULL
for (i in 1:nrow(datos)) {
  if (datos[i,1] == Nn) {
    grupin = n
  }
  if (datos[i,1] != Nn) {
    Nn = datos[i,1]
    n=n+1
    grupin = n
  }
  NewGroupsDB=rbind(NewGroupsDB,as.vector(grupin))
  
}
as.vector(NewGroupsDB)
dim(NewGroupsDB)    #2978
#NewGroupsMetabolomics<-NewGroupsDB

ata10 <- cbind (NewGroupsDB, datos)
ata10[1:10,c(1,2,3,4,5,6:8,24:26,31,32)]   
ata10[1:100,c(1,2,3,4,24:26,32)]   
###################################################################################################################
###########        Creation of the column of the first AAB that appears with a for loop       #####################
###################################################################################################################
### For Facility, we call ata9<-DB12m and perform the for loop 
#Transform to months the appearance of the AAb
colnames (ata10)
ata9<-ata10
ata9[1:5,11:13]

#Transform to months the appearance of the AAb
ata9[,11:13] <- round(ata9[,11:13]/30)


ata9[1:14,11:15]
ata9[,11:13][is.na(ata9[,11:13])] <- 9999


firstAAb = NULL
for(i in 1:nrow(ata9)) {
  AAb = NULL
  if((ata9[i,11]<ata9[i,12]) & (ata9[i,11]<ata9[i,13])){
    AAb = 1#"GADA"
  }
  if((ata9[i,12]<ata9[i,11]) & (ata9[i,12]<ata9[i,13])){ 
    AAb = 2#"IAA"
  }
  if((ata9[i,13]<ata9[i,11]) & (ata9[i,13]<ata9[i,12])){ 
    AAb = 3#"IA2A"
  }
  if((ata9[i,11]==ata9[i,12]) & (ata9[i,11]==ata9[i,13])){ 
    AAb = 4#"IAA/GADA/IA2A"
  }
  if((ata9[i,11]==ata9[i,12]) & (ata9[i,11]<ata9[i,13])){ 
    AAb = 5#"IAA/GADA"
  }
  if((ata9[i,12]==ata9[i,13]) & (ata9[i,12]<ata9[i,11])){ 
    AAb = 6#"IAA/IA2A"
  }
  if((ata9[i,11]==ata9[i,13]) & (ata9[i,11]<ata9[i,12])){ 
    AAb = 7#"GADA/IA2A"
  }
  if((ata9[i,11]==ata9[i,12]) & (ata9[i,11]==ata9[i,13]) & (ata9[i,11]==9999)){
    AAb = 0 #Controls
  }
  firstAAb=rbind(firstAAb, AAb)
}

dim(firstAAb)
as.vector(firstAAb)

ata9[,11:13][ata9[,11:13]==9999] <- NA
ata10<-cbind(as.vector(firstAAb),ata9)
colnames(ata10)[1]<-"FirstAAb"
ata10[1:20,c(1,2,12:14,22:30)]
dim(ata10)   #4845   x  73
head(ata10)

########################################################################################################################
#    Removing cases that are used as controls 
ata10[1:10,c(1:4,9:11,27,29:30)]
aabs<-rowSums(ata10[,9:11]) # From Persist.Conf.Gad to Persist.Conf.Ia2a (9:11)
ata11<-cbind(aabs,ata10)
head(ata11)
dim(ata11)

casesusedascontrols<-ata11[ata11$Outcome==0 & ata11$aabs!=0,]
dim(casesusedascontrols) #150
casesusedascontrols[,c(1:5,7:10,27,28,30)]

ata12<-ata11[!(ata11$Outcome==0 & ata11$aabs!=0),]; dim(ata11);dim(casesusedascontrols);dim(ata12)
#
sum(ata11$Outcome==0 & ata11$aabs!=0);sum(ata12$Outcome==0 & ata12$aabs!=0)

############################################################
#########      Representation by the First AAb      ########
DB12mnoimputednocasesusedascontrols<-ata12
length(unique(DB12mnoimputednocasesusedascontrols$FirstAAb))  # 8
table(DB12mnoimputednocasesusedascontrols$FirstAAb)
#   0    1    2    3    4    5    6    7 
#3483  381  609   15   17  172   12    6 

##################################################################################################################
#########     Create the Case-Control Groups based on the separation due to Agemonth and CASE.IND      ###########
##################################################################################################################
ata7<-DB12mnoimputednocasesusedascontrols
# Eliminate aabs column used to eliminate all cases that were used as controls
ata7[1:10,1:2]
ata7<-ata7[,-1]
colnames(ata7)
ata7<-ata7[with(ata7, order(ata7[,27], ata7[,7], -ata7[,29])),]  #Ordered by CASE.IND(27), Agemonth(7) and -Outcome (29)
ata7[1:10,c(1,2,22,24:27,28,30,31,32)]

##################################################################################################################
# Second step (It has to be done with respect to Outcome)
# For the characteristics of the data, the tables must be ordered again by NewGroupsDB(2), Outcome (29)
ata9<-ata7
ata9<-ata9[with(ata9, order(ata9[,2],-ata9[,29])),]  
ata9[1:50,c(1:6,11:13,28,29,30)]
ata9[1:50,c(1:6,29,30)]

firstAAbCC = NULL
AAbCC = NULL
truquin=NULL
#for(i in 1:length(i)) {
for(i in 1:length(ata9$NewGroupsDB)) {
  
  if(ata9[i,29] != 0) {               # Outcome
    AAbCC = ata9[i,1]
    truquin=ata9[i,1]
  }
  if(ata9[i,29] == 0 ) {              # Outcome
    AAbCC =truquin
  }
  firstAAbCC=rbind(firstAAbCC, AAbCC)
}
dim(firstAAbCC)
as.vector(firstAAbCC)

ata91<-cbind(as.vector(firstAAbCC),ata9)
colnames(ata91)[1]<-"FirstAAbCC"
ata91[1:60,c(1:5,12:14,27,28,29)]
dim(ata91)   #5467   x  993
summary(ata91$FirstAAbCC)
ata91[ata91$CASE.IND==47,c(1:4,6,29,30)]

colnames(ata91)
##################################################################################################################
# Creation of AgeGroupCC (It has to be done with respect to Outcome)
ata8<-ata91
ata8[1:20,c(1,5,6,28,30)]

firstAAbCC = NULL
AAbCC = NULL
truquin=NULL
#for(i in 1:length(i)) {
for(i in 1:length(ata8$NewGroupsDB)) {
  
  if(ata8[i,30] != 0) {            # Outcome = 30
    AAbCC = ata8[i,5]              # AgeGroup = 5
    truquin=ata8[i,5]              # AgeGroup = 5
  }
  if(ata8[i,30] == 0 ) {
    AAbCC =truquin
  }
  firstAAbCC=rbind(firstAAbCC, AAbCC)
}
as.vector(firstAAbCC)
ata8[,30]
dim(firstAAbCC)
dim(ata8)   #4695   x  133

ata81<-cbind(as.vector(firstAAbCC),ata8)
colnames(ata81)[1]<-"AgeGroupCC"
ata81[1:30,c(1:3,6,7,29,31)]
dim(ata81)    #5283 x 999

##################################################################################################################
#### Now creating the First&SecAAb column that has the information about the first two aab that appears.
# Controls must have the value of the first AAb that the case have
# So this have to be done in two steps, first, the AAb in the cases and NA in controls, 
# Second: the controls assume the first AAb of the case they are attached to.
# Control = 0 (for the moment); GADAIAA = 1; GADAIA2A = 2; IAAGADA = 3; IAAIA2A = 4; IA2AGADA = 5; IA2AIAA = 6; 
#                               IAAGADAsametime = 7;IAAIA2Asametime = 8;GADAIA2Asametime = 9;allsametime = 10; GADAothertwo = 11; 
#                               IAAothertwo = 12; IA2Aothertwo = 13.
# The AgeGroup also changes! We have to create a new agegroup
##################################################################################################################

ata8<-ata81
dim(ata8)   #4695 x 75
ata8<- ata8[,c(1:13,18,14:75)]
colnames (ata8)
ata8[1:10,c(1:6,15:17)]   # 15 es GADA, 16 es IAA, 17 es IA2A

ata8[,15:17][is.na(ata8[,15:17])] <- 9999
ata8[1:10,c(1:2,15:17,26)]


firstAAb = NULL
for(i in 1:nrow(ata8)) {
  AAb = NULL
  
  if((ata8[i,15]==ata8[i,16]) & (ata8[i,15]==ata8[i,17]) & (ata8[i,15] == 9999)){
    AAb = 0
  }
  if((ata8[i,15]<ata8[i,16]) & (ata8[i,15]<ata8[i,17]) & (ata8[i,16] == 9999) & (ata8[i,17] == 9999)){
    AAb = 0
  }
  if((ata8[i,16]<ata8[i,15]) & (ata8[i,16]<ata8[i,17]) & (ata8[i,15] == 9999) & (ata8[i,17] == 9999)){
    AAb = 0
  }
  if((ata8[i,17]<ata8[i,15]) & (ata8[i,17]<ata8[i,16]) & (ata8[i,15] == 9999) & (ata8[i,16] == 9999)){
    AAb = 0
  }
  if((ata8[i,15]<ata8[i,16]) & (ata8[i,16]<ata8[i,17]) & (ata8[i,16] != 9999)){
    AAb = 1
  }
  if((ata8[i,15]<ata8[i,17]) & (ata8[i,17]<ata8[i,16])){ 
    AAb = 2
  }
  if((ata8[i,16]<ata8[i,15]) & (ata8[i,15]<ata8[i,17])){ 
    AAb = 3
  }
  if((ata8[i,16]<ata8[i,17]) & (ata8[i,17]<ata8[i,15])){ 
    AAb = 4 
  }
  if((ata8[i,17]<ata8[i,15]) & (ata8[i,15]<ata8[i,16])){ 
    AAb = 5
  }
  if((ata8[i,17]<ata8[i,16]) & (ata8[i,16]<ata8[i,15])){ 
    AAb = 6
  }
  if((ata8[i,15]==ata8[i,16]) & (ata8[i,15]<ata8[i,17])){ 
    AAb = 7
  }
  if((ata8[i,16]==ata8[i,17]) & (ata8[i,16]<ata8[i,15])){ 
    AAb = 8
  }
  if((ata8[i,15]==ata8[i,17]) & (ata8[i,15]<ata8[i,16])){ 
    AAb = 9
  }
  if((ata8[i,15]==ata8[i,16]) & (ata8[i,15]==ata8[i,17]) & (ata8[i,15] != 9999)){ 
    AAb = 10
  }
  if((ata8[i,15]<ata8[i,16]) & (ata8[i,16]==ata8[i,17]) & (ata8[i,16] != 9999)){ 
    AAb = 11
  }
  if((ata8[i,16]<ata8[i,15]) & (ata8[i,15]==ata8[i,17]) & (ata8[i,15] != 9999)){ 
    AAb = 12
  }
  if((ata8[i,17]<ata8[i,15]) & (ata8[i,15]==ata8[i,16]) & (ata8[i,15] != 9999)){ 
    AAb = 13
  }
  firstAAb=rbind(firstAAb, AAb)
}

dim(firstAAb)
as.vector(firstAAb)



head(ata8)
ata8[,15:17][ata8[,15:17]==9999] <- NA
ata9<-cbind(as.vector(firstAAb),ata8)
colnames(ata9)[1]<-"FirstandSecondAAb"
ata9[3235:3255,c(1,13:18,27,28,29,30)]
dim(ata9)   #4845   x  125
summary(ata9$FirstandSecondAAb)
ata9[1:150,c(1,16:18,27,28,29,30)]

ata9[ata9$Mask.Id==756775,] # Example of control that after a while becomes a case not anymore
ata9[1:15,c(1:6,16:18,31,33)]

#############################################################################################
# Second step (It has to be done with respect to Outcome, column )
colnames(ata9)
firstAAbCC = NULL
AAbCC = NULL
truquin=NULL
#for(i in 1:length(i)) {
for(i in 1:length(ata9$NewGroupsDB)) {
  
  if(ata9[i,33] != 0) {       # Outcome = 33
    AAbCC = ata9[i,1]         # FirstandSecondAAb = 1
    truquin=ata9[i,1]
  }
  if(ata9[i,33] == 0 ) {
    AAbCC =truquin
  }
  firstAAbCC=rbind(firstAAbCC, AAbCC)
}
dim(firstAAbCC)
as.vector(firstAAbCC)


ata91<-cbind(as.vector(firstAAbCC),ata9)
colnames(ata91)[1]<-"FirstandsecondAAbCC"
ata91[1:100,c(1:3,12:19)]
dim(ata91)   #4845   x  128
summary(ata91$FirstandsecondAAbCC)
ata91[ata91$CASE.IND==47,]
##################################################################################################################################
#### Creation of the variable AgeGroup2AAbs
# Rounding and dividing by 30 the variables of AAb appearance
ata91[1:10,c(1:3,12:19)]
ata91[,17:19][is.na(ata91[,17:19])] <- 9999

#For easyness, ata7<- ata91
ata7<-ata91

AgeGroup = NULL

for(i in 1:nrow(ata7)) {
  
  grupin = NULL
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]==ata7[i,19]) & (ata7[i,17] == 9999)){
    grupin = 0
  }
  
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,17] == 9999) & (ata7[i,19]!= 9999)){
    grupin = 5
  }
  if((ata7[i,19]==ata7[i,17]) & (ata7[i,19] == 9999) & (ata7[i,18]!= 9999)){
    grupin = 5
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18] == 9999) & (ata7[i,17]!= 9999)){
    grupin = 5
  }
  
  if((ata7[i,17]<ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] > 37)){
    grupin = 4
  }
  if((ata7[i,17]<ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] > 37)){
    grupin = 4
  }
  
  if((ata7[i,18]<ata7[i,19]) & (ata7[i,19]<ata7[i,17]) & (ata7[i,19] > 37)){
    grupin = 4
  }
  if((ata7[i,18]<ata7[i,17]) & (ata7[i,17]<ata7[i,19]) & (ata7[i,17] > 37)){
    grupin = 4
  }
  
  if((ata7[i,19]<ata7[i,17]) & (ata7[i,17]<ata7[i,18]) & (ata7[i,17] > 37)){
    grupin = 4
  }
  if((ata7[i,19]<ata7[i,18]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] > 37)){
    grupin = 4
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] > 37)){
    grupin = 4
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] > 37)){
    grupin = 4
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] > 37)){
    grupin = 4
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]>ata7[i,19]) & (ata7[i,18] > 37) & (ata7[i,18] < 9999)){
    grupin = 4
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]>ata7[i,18]) & (ata7[i,19] > 37) & (ata7[i,19] < 9999)){
    grupin = 4
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]>ata7[i,17]) & (ata7[i,18] > 37) & (ata7[i,18] < 9999)){
    grupin = 4
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]==ata7[i,19]) & (ata7[i,17] > 37) & (ata7[i,17] != 9999)){
    grupin = 4
  }
  
  
  
  
  
  
  
  if((ata7[i,17]<ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 37)){
    grupin = 3
  }
  if((ata7[i,17]<ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 37)){
    grupin = 3
  }
  if((ata7[i,18]<ata7[i,19]) & (ata7[i,19]<ata7[i,17]) & (ata7[i,19] <= 37)){
    grupin = 3
  }
  if((ata7[i,18]<ata7[i,17]) & (ata7[i,17]<ata7[i,19]) & (ata7[i,17] <= 37)){
    grupin = 3
  }
  if((ata7[i,19]<ata7[i,17]) & (ata7[i,17]<ata7[i,18]) & (ata7[i,17] <= 37)){
    grupin = 3
  }
  if((ata7[i,19]<ata7[i,18]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 37)){
    grupin = 3
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 37)){
    grupin = 3
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 37)){
    grupin = 3
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 37)){
    grupin = 3
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]>ata7[i,19]) & (ata7[i,18] <= 37) & (ata7[i,18] < 9999)){
    grupin = 3
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]>ata7[i,18]) & (ata7[i,19] <= 37) & (ata7[i,19] < 9999)){
    grupin = 3
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]>ata7[i,17]) & (ata7[i,18] <= 37)& (ata7[i,18] < 9999)){
    grupin = 3
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]==ata7[i,19]) & (ata7[i,17] <= 37)){
    grupin = 3
  }
  
  
  
  
  if((ata7[i,17]<ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 25)){
    grupin = 2
  }
  if((ata7[i,17]<ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 25)){
    grupin = 2
  }
  if((ata7[i,18]<ata7[i,19]) & (ata7[i,19]<ata7[i,17]) & (ata7[i,19] <= 25)){
    grupin = 2
  }
  if((ata7[i,18]<ata7[i,17]) & (ata7[i,17]<ata7[i,19]) & (ata7[i,17] <= 25)){
    grupin = 2
  }
  if((ata7[i,19]<ata7[i,17]) & (ata7[i,17]<ata7[i,18]) & (ata7[i,17] <= 25)){
    grupin = 2
  }
  if((ata7[i,19]<ata7[i,18]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 25)){
    grupin = 2
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 25)){
    grupin = 2
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 25)){
    grupin = 2
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 25)){
    grupin = 2
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]>ata7[i,19]) & (ata7[i,18] <= 25) & (ata7[i,18] < 9999)){
    grupin = 2
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]>ata7[i,18]) & (ata7[i,19] <= 25) & (ata7[i,19] < 9999)){
    grupin = 2
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]>ata7[i,17]) & (ata7[i,18] <= 25)& (ata7[i,18] < 9999)){
    grupin = 2
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]==ata7[i,19]) & (ata7[i,17] <= 25)){
    grupin = 2
  }
  
  
  if((ata7[i,17]<ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 13)){
    grupin = 1
  }
  if((ata7[i,17]<ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 13)){
    grupin = 1
  }
  if((ata7[i,18]<ata7[i,19]) & (ata7[i,19]<ata7[i,17]) & (ata7[i,19] <= 13)){
    grupin = 1
  }
  if((ata7[i,18]<ata7[i,17]) & (ata7[i,17]<ata7[i,19]) & (ata7[i,17] <= 13)){
    grupin = 1
  }
  if((ata7[i,19]<ata7[i,17]) & (ata7[i,17]<ata7[i,18]) & (ata7[i,17] <= 13)){
    grupin = 1
  }
  if((ata7[i,19]<ata7[i,18]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 13)){
    grupin = 1
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]<ata7[i,19]) & (ata7[i,18] <= 13)){
    grupin = 1
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]<ata7[i,18]) & (ata7[i,19] <= 13)){
    grupin = 1
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]<ata7[i,17]) & (ata7[i,18] <= 13)){
    grupin = 1
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]>ata7[i,19]) & (ata7[i,18] <= 13) & (ata7[i,18] < 9999)){
    grupin = 1
  }
  if((ata7[i,17]==ata7[i,19]) & (ata7[i,19]>ata7[i,18]) & (ata7[i,19] <= 13) & (ata7[i,19] < 9999)){
    grupin = 1
  }
  if((ata7[i,18]==ata7[i,19]) & (ata7[i,18]>ata7[i,17]) & (ata7[i,18] <= 13)& (ata7[i,18] < 9999)){
    grupin = 1
  }
  if((ata7[i,17]==ata7[i,18]) & (ata7[i,18]==ata7[i,19]) & (ata7[i,17] <=13)){
    grupin = 1
  }
  
  AgeGroup=c(AgeGroup,grupin)
}

length(AgeGroup)
dim(ata7)
ata7[370:380,c(1:2,17:19)]
as.vector(AgeGroup)

ata8<-cbind(AgeGroup,ata7)
ata8[1:200,c(1:10)]
dim(ata8)  #3094 x 112     4845 for 13 months

colnames(ata8)[1]<-"AgeGroup2AAbs"
ata8[3235:3255,c(1:7,9,10)]
dim(ata8)   #4695   x  133

########################################################################################################################
# Second step (It has to be done with respect to Outcome)
colnames(ata8)
ata8[1:20,c(1:7,9,10)]
firstAAbCC = NULL
AAbCC = NULL
truquin=NULL
#for(i in 1:length(i)) {
for(i in 1:length(ata8$NewGroupsDB)) {
  
  if(ata8[i,35] != 0) {            # Outcome = 35
    AAbCC = ata8[i,1]
    truquin=ata8[i,1]
  }
  if(ata8[i,35] == 0 ) {
    AAbCC =truquin
  }
  firstAAbCC=rbind(firstAAbCC, AAbCC)
}
as.vector(firstAAbCC)
ata8[,35]
dim(firstAAbCC)
dim(ata8)   #4695   x  133

ata81<-cbind(as.vector(firstAAbCC),ata8)
colnames(ata81)[1]<-"AgeGroup2AAbsCC"
ata81[1:30,c(1:3,36)]
dim(ata81) # 4695 x 80
colnames(ata81)

##################################################################################################################################
# Creating a set of groups of Sconv and Outcome so we can be able to graph the canonical biplot
##################################################################################################################################
colnames(ata81)

# First we have to order the data
ata92<-ata81
# Order by Sconv (11), -Outcome (36)
ata92<-ata92[with(ata92, order(ata92[,11],-ata92[,36])),]  #Ordered by SConv(11) and Outcome(36)
ata92[1:40,c(7:11,36)]


gruposSconvOutcome <-apply(ata92[,c(11,36)], 1, paste, collapse ="_")    #Sconv; Outcome
head(gruposSconvOutcome)
UNIgruposSconvOutcome <- unique (gruposSconvOutcome)
length(UNIgruposSconvOutcome)     #10 groups

ata93 <- cbind (gruposSconvOutcome, ata92)
ata93[1:10,c(1,2,3,19:31)]   #annots (3094 x 28)
dim(ata93)

# For Loop to create the values of the 10 groups
NewGroupsmetabol=NULL
n=0
Nn ="word"
grupin=NULL
for (i in 1:nrow(ata93)) {
  if (ata93[i,1] == Nn) {
    grupin = n
  }
  if (ata93[i,1] != Nn) {
    Nn = ata93[i,1]
    n=n+1
    grupin = n
  }
  NewGroupsmetabol=rbind(NewGroupsmetabol,as.vector(grupin))
  
}
as.vector(NewGroupsmetabol)
dim(NewGroupsmetabol)    #5449
GroupsSconvOutcome<-NewGroupsmetabol

ata93 <- cbind (GroupsSconvOutcome, ata93)
ata93[2000:2020,c(1:3)]   
dim(ata93)  # 4695 x 82
ata93$GroupsSconvOutcome
colnames(ata93)
# Again we reorder the table by Sconv (13), CASE.IND (36) and Outcome (38)
ata93<-ata93[with(ata93, order(ata93[,13],ata93[,36],-ata93[,38])),]  #Ordered by Sconv, CASE.IND and Outcome
ata93[1:10,c(1:2,8:11,32:38)]

DB12mnoimputednocasesusedascontrolsoctober<-ata93
#write.table(DB12mnoimputednocasesusedascontrolsoctober,  "/home/leobalzano/Documents/TEDDY data/MP149/Diet/metabolites/DB12mnoimputednocasesusedascontrolsoctober.txt", sep="\\t")
#save(file = "DB12mnoimputednocasesusedascontrolsoctober.RData", DB12mnoimputednocasesusedascontrolsoctober)
##########################################################################################################################
# Eliminating all the all-NA rows 
who.is.NA <- as.vector(apply(DB12mnoimputednocasesusedascontrolsoctober, 1, function (x) all(is.na(x[40:82]))))  # identify which individuals have all DB values at NA
datos4 <- DB12mnoimputednocasesusedascontrolsoctober[which(!who.is.NA),] ; dim(DB12mnoimputednocasesusedascontrolsoctober); dim(datos4)  # remove these individuals
sort(table(DB12mnoimputednocasesusedascontrolsoctober$NewGroupsDB))
sort(table(datos4$NewGroupsDB))
#sort(data.frame(table(datos4$NewGroupsDB)))
colnames(datos4)
datos4order<-datos4[with(datos4, order(datos4[,10],-datos4[,38])),]  #Ordered by NewGroupsDB(10) and -Outcome (38)
datos4order[1:40,]
length(unique(datos4order$NewGroupsDB))
tabla<-data.frame(table(datos4order$NewGroupsDB))
head(tabla)
dim(tabla)
colnames(tabla)[1]<- "NewGroupsDB"

datos5<-merge(tabla,datos4order, by= "NewGroupsDB");dim(tabla);dim(datos4order);dim(datos5)
datos5[1:100,1:5]
colnames(datos5)
datos5<-datos5[with(datos5, order(datos5[,1],-datos5[,39])),]  #Ordered by NewGroupsDB(1) and -Outcome (39)
dim(DB12mnoimputednocasesusedascontrolsoctober);dim(datos5)
dim(datos5[datos5$Outcome==1,])


totvalue=NULL
value=NULL
for(i in 1:nrow(datos5)) {
  if (datos5$Freq[i] >1) {
    value=1
  }
  if (datos5$Freq[i] ==1) {
    value=0.5
  }
  totvalue<-rbind(totvalue,value)  
}
dim(totvalue)
summary (totvalue)
head(totvalue)
colnames(totvalue)<-"heatmapfactor"

datos6<- cbind(datos5,totvalue)
datos6[1:10,1:10]
datos6<-datos6[,c(84,1:83)]
colnames(datos6)
datos61<-datos6[,c(1:3,16)]
datos7<-merge (DB12mnoimputednocasesusedascontrolsoctober, datos61,by = c("NewGroupsDB","Mask.Id"), all.x = TRUE)
dim(datos61);dim(DB12mnoimputednocasesusedascontrolsoctober);dim(datos7)
#
colnames(datos7)
datos7[1:200,c(1,14,38,83,84)]

#Transform NA in heatmapfactor to 0
datos7[,83][is.na(datos7[,83])] <- 0

DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor<-datos7
#write.table(DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor,  "/home/leobalzano/Documents/TEDDY data/MP149/Diet/metabolites/DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor.txt", sep="\\t")
datos7[datos7$Mask.Id==204614,]


head(DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor)

how.many.NA <- apply(datos4[40:82], 2, function (x) length(which(is.na(x))))
how.many.NA # counts the numbers of NAs per variable. Less is VitaminD with 38 NAs and most is lutein and others with 728

length(unique(datos4$Mask.Id))  # 1453 individuals with Dietary biomarker data  Igual luego de los cambios de octubre
length(unique(datos4$Mask.Id[datos4$Outcome ==1])) # 378 cases with Dietary biomarker data Igual luego de los cambios de octubre
length(unique(datos4$Mask.Id[datos4$Outcome ==0])) # 1075 cases with Dietary biomarker data Igual luego de los cambios de octubre

###########################################################################################################
############################     Function To calculate the internal substraction     ######################
###########################################################################################################
#Eliminating the groups with just one member
sort(table(DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor$NewGroupsDB))
DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor[1:5,c(10,13,14,38)]
datos4[datos4$NewGroupsDB== 9,c(10,13,14,38)]
datos4[datos4$Mask.Id== 881791,c(10,13,14,38)]

# This is the data to perform the heatmaps.
DBoctobernoimputed<-datos4
getwd()

# Function
substraction<- function (X, design) 
{
  X = as.matrix(X)
  rep.measures = factor(design[, 1])
  factors = design[, -1, drop = FALSE]
  if ((ncol(factors) == 0) | (ncol(factors) == 1)) {
    indiv.names = rownames(X)
    rownames(X) = as.character(rep.measures)
    X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, 
                                mean, na.rm = TRUE), nrow = length(unique(rep.measures)), 
                          ncol = dim(X)[2], dimnames = list(levels(as.factor(rep.measures)), 
                                                            colnames(X)))
    Xb = X.mean.indiv[as.character(rep.measures), ]
    Xw = X - Xb
    dimnames(Xw) = list(indiv.names, colnames(X))
  }
  return(Xw)
}

# X=datos6, design = NewGroupsDB
colnames(DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor)
DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor[1:3,40:82]

X<-DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor[,40:82]
colnames(X)

design<- data.frame(sample = DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor$NewGroupsDB)

Xw<- substraction (X=X, design=design) 


X[25:28,1:2]
mean(X[25:28,1])

Xw [25:28,1:2]
head(X)
head(Xw)
data.frame(design[25:28,])
head(data.frame(design))
dim(Xw)

colnames(datos6)
datos7<-cbind(DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor[,1:41],Xw)

dim(datos6)
dim(datos7)

#DietaryBiomarkersOctoberX <- DB12mnoimputednocasesusedascontrolsoctoberwithheatmapfactor
#write.table(DietaryBiomarkersOctoberX,  "/home/leobalzano/Documents/TEDDY data/MP149/Diet/Fortucker/DietaryBiomarkersOctoberX.txt", sep="\\t")
#DietaryBiomarkersOctoberXw <- datos7
#write.table(DietaryBiomarkersOctoberXw,  "/home/leobalzano/Documents/TEDDY data/MP149/Diet/Fortucker/DietaryBiomarkersOctoberXw.txt", sep="\\t")
###########################################################################################################

###########################################################
###########################################################
#################     Checkpoint 2     ####################
###########################################################
###########################################################
# From this point on you can continue just by importing DietaryBiomarkersOctoberXw
DietaryBiomarkersOctoberXw <- read.delim ("/media/data/leobalzano/ScriptsForTEDDY/Data/DietaryBiomarkersOctoberXw.txt", sep= "\\t")
dim(DietaryBiomarkersOctoberXw)
DietaryBiomarkersOctoberXw [50:100,74:84]

###########################################################
#############     Creating the array     ##################
###########################################################
laXwDB<-DietaryBiomarkersOctoberXw
colnames(laXwDB)
laXwDB<-laXwDB[,-c(42,43)]
laXwDBout1<-laXwDB
#laXwDBout1<- subset(laXwDB, laXwDB$Outcome ==1)

dim(laXwDBout1)  #1212 x 82 
colnames(laXwDBout1)

# Order by Mask.Id Variables has Mask.Id(2), Sconv(14), FirstAAb(11), AgeGroup(13) and the rest
laXwDBout1<-laXwDBout1[with(laXwDBout1,order(laXwDBout1[,2])),]

variablesdb<-laXwDBout1
head(variablesdb)  

#####################################################################################
dim(variablesdb)
variablesdb[1:10,1:40]
allcasesdb<-unique(variablesdb$Mask.Id)
length(allcasesdb)  # 1621
#####################################################################################
colnames(variablesdb)
VariablesIADBDS<-variablesdb[,1:39]

#####################################################################################
DietSummaryOctoberXw <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/DietSummaryOctoberXw.txt", sep="\\t")
laXwDB<-DietSummaryOctoberXw
colnames(laXwDB)
laXwDBout1<- laXwDB
dim(laXwDBout1)  #1078 x 85 
colnames(laXwDBout1)

# Order by Mask.Id Variables has Mask.Id(2), Sconv(9), FirstAAb(6), AgeGroup (8) and the rest
laXwDBout1<-laXwDBout1[with(laXwDBout1,order(laXwDBout1[,2])),]
variablesds<-laXwDBout1
head(variablesds)  

#####################################################################################
dim(variablesds)
variablesds[1:10,1:40]
allcasesds<-unique(variablesds$Mask.Id)
length(allcasesds)  # 1580
#####################################################################################
allcasesdb
allcasesds

setdiff(allcasesdb,allcasesds)
allcasesunion<-union(allcasesdb,allcasesds)
length(allcasesunion)   # 1621
length(allcasesdb)
length(allcasesds)

is.element(allcasesdb,allcasesds)
is.element(allcasesds,allcasesdb)

#####################################################################################
allcasesunion<-data.frame(allcasesunion)
colnames(allcasesunion)[1]<-"Mask.Id"

#####################################################################################
#Groups by timepoint
Grt12db<- subset (variablesdb,variablesdb[,14] == "-12")
Grt9db<- subset (variablesdb,variablesdb[,14] == "-9")
Grt6db<- subset (variablesdb,variablesdb[,14] == "-6")
Grt3db<- subset (variablesdb,variablesdb[,14] == "-3")
Grt0db<- subset (variablesdb,variablesdb[,14] == "0")

# Merging with all cases
Grt12dbtotal<-merge(allcasesunion,Grt12db, by="Mask.Id", all.x = T);dim(Grt12dbtotal);dim(allcasesunion);dim(Grt12db)
Grt12dbtotal[1:5,1:5]
rownames(Grt12dbtotal)<- Grt12dbtotal[,1]
colnames(Grt12dbtotal)
Grt12dbtotal<- as.matrix(Grt12dbtotal[,c(-1:-39)])
dim(Grt12dbtotal)
head(Grt12dbtotal)

Grt9dbtotal<-merge(allcasesunion,Grt9db, by="Mask.Id", all.x = T);dim(Grt9dbtotal);dim(allcasesunion);dim(Grt9db)
rownames(Grt9dbtotal)<- Grt9dbtotal[,1]
Grt9dbtotal<- as.matrix(Grt9dbtotal[,c(-1:-39)])
dim(Grt9dbtotal)
head(Grt9dbtotal)

Grt6dbtotal<-merge(allcasesunion,Grt6db, by="Mask.Id", all.x = T);dim(Grt6dbtotal)
rownames(Grt6dbtotal)<- Grt6dbtotal[,1]
Grt6dbtotal<- as.matrix(Grt6dbtotal[,c(-1:-39)])
dim(Grt6dbtotal)

Grt3dbtotal<-merge(allcasesunion,Grt3db, by="Mask.Id", all.x = T);dim(Grt3dbtotal)
rownames(Grt3dbtotal)<- Grt3dbtotal[,1]
Grt3dbtotal<- as.matrix(Grt3dbtotal[,c(-1:-39)])
dim(Grt3dbtotal)

Grt0dbtotal<-merge(allcasesunion,Grt0db, by="Mask.Id", all.x = T);dim(Grt0dbtotal)
rownames(Grt0dbtotal)<- Grt0dbtotal[,1]
Grt0dbtotal<- as.matrix(Grt0dbtotal[,c(-1:-39)])
dim(Grt0dbtotal)

# Dimensions are 418*43*5
colnames(Grt0dbtotal)
Grt0dbtotal[1:10,1:4]
dim(Grt0dbtotal)

#################################################################################
###########################     3D array     ####################################
#################################################################################
# Dietary Biomarkers
ArrayDietBiomindiv <- array(data = NA, dim = c(1621,43,5),dimnames = list(NULL, NULL, c("-12","-9","-6", "-3", "0")))


ArrayDietBiomindiv
ArrayDietBiomindiv[,,1] <- Grt12dbtotal
ArrayDietBiomindiv[,,2] <- Grt9dbtotal
ArrayDietBiomindiv[,,3] <- Grt6dbtotal
ArrayDietBiomindiv[,,4] <- Grt3dbtotal
ArrayDietBiomindiv[,,5] <- Grt0dbtotal

rownames(ArrayDietBiomindiv)<-rownames(Grt12dbtotal)
colnames(ArrayDietBiomindiv)<-colnames (Grt12dbtotal)
ArrayDietBiomindiv
dim(ArrayDietBiomindiv)   # 1621 * 43 * 5
getwd()
# save(ArrayDietBiomindiv, file = "ArrayDietBiomindiv.RData")

####################################################################################################################
# Determining the best fitted model
modelDB<-bestfittedmodel (X=ArrayDietBiomindiv,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
modelDB
# the best minimal fitted model would be 4,2,3

# Impute method
FullarrayDB<-Imputemethod(X=ArrayDietBiomindiv,fac=c(4, 2, 3), conver = 1e-07, max.iter = 1000)
summary(FullarrayDB)

###########################################################
###########################################################
#################     Checkpoint 3     ####################
###########################################################
###########################################################
dim(FullarrayDB)

FullarrayDB136<- FullarrayDB[rownames(FullarrayDB) %in% rownames(outcomedummyarray136),,]
dim(FullarrayDB136)

##################################################################################
###   NPLSDA
NPLSDAFullarrayDB136<-NPLSDAmod(XN=FullarrayDB136, YN=outcomedummyarray136,centering = 0)
summary(NPLSDAFullarrayDB136)
NPLSDAFullarrayDB136$NPLSDAQ22D

plotNPLSDAmod(X=NPLSDAFullarrayDB136)
