###########################################################
####     PartialCorrelationNetworksGADAIAAxTEDDY.R     ####
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to perform the Partial Correlation networks for GADA first and IAA first datasets
# We Perform ParCor to the cases only, then subset the correlation matrix and plot the network
# Also, this analysis, we are selecting the central core of correlation by a Shrinkage procedure...
# ... based on the lowest root mean square error of prediction (RMSEP) value from...
# ... Ridge, Lasso or ElasticNet methods to guarantee the best fitted model.
# This script is son of PartialCorrelationNetworksxTEDDY.R on which we already made the ENselected vars
# using GADA vs IAA, so the selection is the same, what is different is that we need to subset the individuals to 
# GADA first or IAA first, then perform ParCor and then subset again to the pre-selected (315) features

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()
###########################################################
# Data: 
wholedescriptivevars <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/wholedescriptivevars.txt")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/DeltaMatrices_FullArray.Rdata")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayDB.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/outcomedummyarray136.RData")
definitivemetabolitesSelectedConvertedmarch9 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/definitivemetabolitesSelectedConvertedmarch9.txt")
AgeGroupdummycasesfrom136 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/AgeGroupdummycasesfrom136.txt", row.names=1)
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayGEMETABDBBestModel.RData")
ENSelVarswithDBGADAIAA <- read.table("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/NetworksAnalyzes/ENSelVarswithDBGADAIAA.txt")


###########################################################
# Libraries:
require("ropls")
library(reshape2)
library("glmnet", lib.loc="/usr/local/lib/R/site-library")
library("qgraph")
library("parcor")
require(reshape)

###########################################################
class(DeltaMatrices_FullArray)
summary(DeltaMatrices_FullArray)
dim(DeltaMatrices_FullArray$`Period 12-9`)   # 136 x 1107
DeltaMatrices_FullArray$`Period 12-9`[1:10,1:10]
###########################################################
# First we have to add to each table the deltas from the three DB variables
dim(FullarrayDB)

FullarrayDB136<- FullarrayDB[rownames(FullarrayDB) %in% rownames(outcomedummyarray136),,]
dim(FullarrayDB136)
colnames(FullarrayDB136) # Retain just vitD (15),VitC (43) and alphaTocopherol(16)
PieceFullarrayDB136<-FullarrayDB136[,c(15,16,43),]
dim(PieceFullarrayDB136)
PieceFullarrayDB136[1:10,,1]
###########################################################
# Changing the colnames of the metabolomics variables
colnames(DeltaMatrices_FullArray$`Period 12-9`)
###  Time period -12 -9
tulo<- DeltaMatrices_FullArray$`Period 12-9`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym12m9<-merge(DeltaMatrices_FullArray$`Period 12-9`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 12-9`); dim (tulo);dim(DeltaMatrices_FullArraym12m9)
rownames(DeltaMatrices_FullArraym12m9)<-DeltaMatrices_FullArraym12m9[,1]
DeltaMatrices_FullArraym12m9 <-DeltaMatrices_FullArraym12m9[,-1]
dim(DeltaMatrices_FullArraym12m9)
DeltaMatrices_FullArraym12m9[1:10,900:910]

###  Time period -9 -6
tulo<- DeltaMatrices_FullArray$`Period 9-6`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym9m6<-merge(DeltaMatrices_FullArray$`Period 9-6`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 9-6`); dim (tulo);dim(DeltaMatrices_FullArraym9m6)
rownames(DeltaMatrices_FullArraym9m6)<-DeltaMatrices_FullArraym9m6[,1]
DeltaMatrices_FullArraym9m6 <-DeltaMatrices_FullArraym9m6[,-1]
dim(DeltaMatrices_FullArraym9m6)
DeltaMatrices_FullArraym9m6[1:10,1100:1107]

###  Time period -6 -3
tulo<- DeltaMatrices_FullArray$`Period 6-3`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym6m3<-merge(DeltaMatrices_FullArray$`Period 6-3`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 6-3`); dim (tulo);dim(DeltaMatrices_FullArraym6m3)
rownames(DeltaMatrices_FullArraym6m3)<-DeltaMatrices_FullArraym6m3[,1]
DeltaMatrices_FullArraym6m3 <-DeltaMatrices_FullArraym6m3[,-1]
dim(DeltaMatrices_FullArraym6m3)
DeltaMatrices_FullArraym6m3[1:10,1100:1107]

###  Time period -3 to 0
tulo<- DeltaMatrices_FullArray$`Period 3-0`[,863:1107]
dim(tulo)
colnames(tulo)
colnames(tulo)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( tulo ) ,definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                               "FunCat" ]

DeltaMatrices_FullArraym3m0<-merge(DeltaMatrices_FullArray$`Period 3-0`[,1:862],tulo, by=0)
dim(DeltaMatrices_FullArray$`Period 3-0`); dim (tulo);dim(DeltaMatrices_FullArraym3m0)
rownames(DeltaMatrices_FullArraym3m0)<-DeltaMatrices_FullArraym3m0[,1]
DeltaMatrices_FullArraym3m0 <-DeltaMatrices_FullArraym3m0[,-1]
dim(DeltaMatrices_FullArraym3m0)
DeltaMatrices_FullArraym3m0[1:10,1100:1107]

###########################################################
# Substraction of tf-to per Dietary biomarker variable
tablem12m9<-PieceFullarrayDB136[,,2] - PieceFullarrayDB136[,,1]
PieceFullarrayDB136[1:10,,1]
PieceFullarrayDB136[1:10,,2]
tablem12m9 [1:10,] # perfect
dim(tablem12m9)
anyNA(tablem12m9)
tablem9m6<-PieceFullarrayDB136[,,3] - PieceFullarrayDB136[,,2]
tablem6m3<-PieceFullarrayDB136[,,4] - PieceFullarrayDB136[,,3]
tablem3m0<-PieceFullarrayDB136[,,5] - PieceFullarrayDB136[,,4]
dim(tablem12m9);dim(tablem9m6);dim(tablem6m3);dim(tablem3m0)

###########################################################
# Merging this tables to the array
### Time period m9 -m12
Deltaallvarsm12m9<- merge (DeltaMatrices_FullArraym12m9,tablem12m9, by= 0)
dim(DeltaMatrices_FullArraym12m9);dim(tablem12m9); dim(Deltaallvarsm12m9)
Deltaallvarsm12m9[1:10,1:10]
rownames(Deltaallvarsm12m9)<-Deltaallvarsm12m9[,1]
Deltaallvarsm12m9 <-Deltaallvarsm12m9[,-1]
Deltaallvarsm12m9[1:10,1:10]
dim(Deltaallvarsm12m9)
colnames(Deltaallvarsm12m9)
Deltaallvarsm12m9[1:10,1100:1110]

### Time period m9 -m6
Deltaallvarsm9m6<- merge (DeltaMatrices_FullArraym9m6,tablem9m6, by= 0)
dim(DeltaMatrices_FullArraym9m6);dim(tablem9m6); dim(Deltaallvarsm9m6)
Deltaallvarsm9m6[1:10,1:10]
rownames(Deltaallvarsm9m6)<-Deltaallvarsm9m6[,1]
Deltaallvarsm9m6 <-Deltaallvarsm9m6[,-1]
Deltaallvarsm9m6[1:10,1:10]
dim(Deltaallvarsm9m6)
colnames(Deltaallvarsm9m6)
### Time period m6 -m3
Deltaallvarsm6m3<- merge (DeltaMatrices_FullArraym6m3,tablem6m3, by= 0)
dim(DeltaMatrices_FullArraym6m3);dim(tablem6m3); dim(Deltaallvarsm6m3)
Deltaallvarsm6m3[1:10,1:10]
rownames(Deltaallvarsm6m3)<-Deltaallvarsm6m3[,1]
Deltaallvarsm6m3 <-Deltaallvarsm6m3[,-1]
Deltaallvarsm6m3[1:10,1:10]
dim(Deltaallvarsm6m3)
colnames(Deltaallvarsm6m3)
### Time period m3 - m0
Deltaallvarsm3m0<- merge (DeltaMatrices_FullArraym3m0,tablem3m0, by= 0)
dim(DeltaMatrices_FullArraym3m0);dim(tablem3m0); dim(Deltaallvarsm3m0)
Deltaallvarsm3m0[1:10,1:10]
rownames(Deltaallvarsm3m0)<-Deltaallvarsm3m0[,1]
Deltaallvarsm3m0 <-Deltaallvarsm3m0[,-1]
Deltaallvarsm3m0[1:10,1:10]
dim(Deltaallvarsm3m0)
colnames(Deltaallvarsm3m0)

###########################################################
###########################################################
#################     Shrinkage method     ################
###########################################################
###########################################################
# This have to be done with the CASES ONLY SEPARATING GADA AND IAA
dim(AgeGroupdummycasesfrom136)

Deltaallvarsm12m9cases<-Deltaallvarsm12m9 [is.element(rownames (Deltaallvarsm12m9),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm12m9cases)

Deltaallvarsm9m6cases<-Deltaallvarsm9m6 [is.element(rownames (Deltaallvarsm9m6),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm9m6cases)
Deltaallvarsm6m3cases<-Deltaallvarsm6m3 [is.element(rownames (Deltaallvarsm6m3),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm6m3cases)
Deltaallvarsm3m0cases<-Deltaallvarsm3m0 [is.element(rownames (Deltaallvarsm3m0),rownames(AgeGroupdummycasesfrom136)),];dim(Deltaallvarsm3m0cases)

###########################################################
# Creating the FAAb vector
# GADA is 1, IAA is 2
colnames(wholedescriptivevars)
preFAAbcases<-wholedescriptivevars[wholedescriptivevars$Mask.Id %in% rownames(Deltaallvarsm12m9cases),]
dim(preFAAbcases)
FAAbcases<-unique(preFAAbcases[,c(1,4)])
dim(FAAbcases)
dim(Deltaallvarsm12m9cases)

FAAbcasesGADAIAA<-FAAbcases[FAAbcases$FirstAAb != 4 & FAAbcases$FirstAAb != 5 & FAAbcases$FirstAAb != 7,]
dim(FAAbcasesGADAIAA)
table(FAAbcasesGADAIAA$FirstAAb)
###########################################################
FAAbcasesGADA<-FAAbcases[FAAbcases$FirstAAb == 1,]
FAAbcasesIAA<-FAAbcases[FAAbcases$FirstAAb == 2,]

###########################################################
# GADA #
Deltaallvarsm12m9GADA<-Deltaallvarsm12m9cases[is.element(rownames(Deltaallvarsm12m9cases),FAAbcasesGADA$Mask.Id),];dim(Deltaallvarsm12m9cases);dim(FAAbcasesGADA);dim(Deltaallvarsm12m9GADA)
Deltaallvarsm9m6GADA<-Deltaallvarsm9m6cases[is.element(rownames(Deltaallvarsm9m6cases),FAAbcasesGADA$Mask.Id),];dim(Deltaallvarsm9m6cases);dim(FAAbcasesGADA);dim(Deltaallvarsm9m6GADA)
Deltaallvarsm6m3GADA<-Deltaallvarsm6m3cases[is.element(rownames(Deltaallvarsm6m3cases),FAAbcasesGADA$Mask.Id),];dim(Deltaallvarsm6m3cases);dim(FAAbcasesGADA);dim(Deltaallvarsm6m3GADA)
Deltaallvarsm3m0GADA<-Deltaallvarsm3m0cases[is.element(rownames(Deltaallvarsm3m0cases),FAAbcasesGADA$Mask.Id),];dim(Deltaallvarsm3m0cases);dim(FAAbcasesGADA);dim(Deltaallvarsm3m0GADA)


# IAA #
Deltaallvarsm12m9IAA<-Deltaallvarsm12m9cases[is.element(rownames(Deltaallvarsm12m9cases),FAAbcasesIAA$Mask.Id),];dim(Deltaallvarsm12m9cases);dim(FAAbcasesIAA);dim(Deltaallvarsm12m9IAA)
Deltaallvarsm9m6IAA<-Deltaallvarsm9m6cases[is.element(rownames(Deltaallvarsm9m6cases),FAAbcasesIAA$Mask.Id),];dim(Deltaallvarsm9m6cases);dim(FAAbcasesIAA);dim(Deltaallvarsm9m6IAA)
Deltaallvarsm6m3IAA<-Deltaallvarsm6m3cases[is.element(rownames(Deltaallvarsm6m3cases),FAAbcasesIAA$Mask.Id),];dim(Deltaallvarsm6m3cases);dim(FAAbcasesIAA);dim(Deltaallvarsm6m3IAA)
Deltaallvarsm3m0IAA<-Deltaallvarsm3m0cases[is.element(rownames(Deltaallvarsm3m0cases),FAAbcasesIAA$Mask.Id),];dim(Deltaallvarsm3m0cases);dim(FAAbcasesIAA);dim(Deltaallvarsm3m0IAA)



###########################################################
# Subsetting both datasets.
dim(FullarrayGEMETABDBBestModel)
FullarrayGEMETABDBBestModel [1:10,1:10,1]
# GADA
FullarrayGEMETABDBBestModelCASESGADA<-FullarrayGEMETABDBBestModel[rownames(FullarrayGEMETABDBBestModel) %in% FAAbcasesGADA$Mask.Id,,]
dim(FullarrayGEMETABDBBestModelCASESGADA)
# IAA
FullarrayGEMETABDBBestModelCASESIAA<-FullarrayGEMETABDBBestModel[rownames(FullarrayGEMETABDBBestModel) %in% FAAbcasesIAA$Mask.Id,,]
dim(FullarrayGEMETABDBBestModelCASESIAA)

###########################################################
###########################################################
####################     NETWORKS     #####################
###########################################################
###########################################################
# 1.- Perform ParCor to the cases only, then subset the correlation matrix and plot the network
# Time Period -12 to -9

#################################################
#############     CORRELATIONS     ##############
#################################################
# Compute correlations:
##############################
### Time Period -12 to -9 ####
##############################
##############################
###          GADA         ####
##############################

Correlation12to9GADA <- cor(Deltaallvarsm12m9GADA, use = "pairwise")  

class(Correlation12to9GADA)
isSymmetric(Correlation12to9GADA)
PCorMat12to9GADA <- cor2pcor(Correlation12to9GADA)
class(PCorMat12to9GADA)
isSymmetric(PCorMat12to9GADA)
PCorMat12to9GADA<-round(PCorMat12to9GADA,digits = 6)
isSymmetric(PCorMat12to9GADA)
colnames(PCorMat12to9GADA)<-colnames(Correlation12to9GADA)
rownames(PCorMat12to9GADA)<-rownames(Correlation12to9GADA)
dim(PCorMat12to9GADA)
min(PCorMat12to9GADA)
max(PCorMat12to9GADA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat12to9GADA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat12to9GADA))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat12to9GADA2<-PCorMat12to9GADA[rownames(PCorMat12to9GADA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat12to9GADA2)
PCorMat12to9GADASizeEN<-PCorMat12to9GADA2[,colnames(PCorMat12to9GADA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat12to9GADASizeEN)
PCorMat12to9GADASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat12to9GADASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat12to9GADASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat12to9GADASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat12to9GADASizeEN_melted)
head(PCorMat12to9GADASizeEN_melted)
colnames(PCorMat12to9GADASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat12to9GADASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat12to9GADASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9GADA[1:10,1:10]
dim(Deltaallvarsm12m9GADA)
vectitodeGE<-colnames(Deltaallvarsm12m9GADA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9GADA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9GADA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable12to9GADA<-PCorMat12to9GADASizeEN_melted[,-3]

NodeTable12to9GADA2<- data.frame(Var=c(as.vector(NodeTable12to9GADA[,1]),as.vector(NodeTable12to9GADA[,2])))
length(as.vector(NodeTable12to9GADA[,1])); length(as.vector(NodeTable12to9GADA[,2]));dim(NodeTable12to9GADA2)
head(NodeTable12to9GADA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable12to9GADA2$Omic<-""
head(NodeTable12to9GADA2)
NodeTable12to9GADA2[1:50,]

NodeTable12to9GADA2$Omic [NodeTable12to9GADA2$Var %in% vectitodeGE] <- 1
NodeTable12to9GADA2$Omic [NodeTable12to9GADA2$Var %in% vectitodeMetab] <- 2
NodeTable12to9GADA2$Omic [NodeTable12to9GADA2$Var %in% vectitodeDB] <- 3
table(NodeTable12to9GADA2$Omic)
#sum(table(NodeTable12to9GADA2$Omic))
NodeTable12to9GADA2[850:1065,]
dim(NodeTable12to9GADA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat12to9GADASizeEN_melted$Correlations)>=thrs) # 294 true
dim(PCorMat12to9GADASizeEN_melted)
PCorMat12to9GADASizeEN_meltedCO0.7<-PCorMat12to9GADASizeEN_melted[abs(PCorMat12to9GADASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat12to9GADASizeEN_meltedCO0.7)
head(PCorMat12to9GADASizeEN_meltedCO0.7)
PCorMat12to9GADASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat12to9GADASizeEN_meltedCO0.7[,1]),as.vector(PCorMat12to9GADASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat12to9GADASizeEN_meltedCO0.7[,2]));dim(PCorMat12to9GADASizeEN_meltedCO0.72)
dim(unique(PCorMat12to9GADASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat12to9GADASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable12to9GADA2$Weight2<-tulo$Freq [match(NodeTable12to9GADA2$Var,tulo$Var1)]
NodeTable12to9GADA2[1:500,]
dim(NodeTable12to9GADA2)

NodeTable12to9GADA2[NodeTable12to9GADA2$Var=="PLEKHA8",]
dim(NodeTable12to9GADA2[NodeTable12to9GADA2$Var=="PLEKHA8",])

NodeTable12to9GADA2[NodeTable12to9GADA2$Var=="unknown272",]
dim(NodeTable12to9GADA2[NodeTable12to9GADA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable12to9GADA2$Weight2[is.na(NodeTable12to9GADA2$Weight2)]<-0
NodeTable12to9GADA2$Weight<-NodeTable12to9GADA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable12to9GADA2$Nameito<- abbreviate(NodeTable12to9GADA2$Var, 7, method = "both")
NodeTable12to9GADA2[NodeTable12to9GADA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm12m9GADA[1:10,1:10]
dim(Deltaallvarsm12m9GADA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm12m9GADA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)

NodeTable12to9GADA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable12to9GADA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable12to9GADA<-unique (NodeTable12to9GADA2)
length(unique (NodeTable12to9GADA2$Var)) # 406
dim(NodeTable12to9GADA)
#write.table(NodeTable12to9GADA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable12to9GADA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -12 to -9 ####
##############################
##############################
###           IAA         ####
##############################

Correlation12to9IAA <- cor(Deltaallvarsm12m9IAA, use = "pairwise")  

class(Correlation12to9IAA)
isSymmetric(Correlation12to9IAA)
PCorMat12to9IAA <- cor2pcor(Correlation12to9IAA)
class(PCorMat12to9IAA)
isSymmetric(PCorMat12to9IAA)
PCorMat12to9IAA<-round(PCorMat12to9IAA,digits = 5)
isSymmetric(PCorMat12to9IAA)
colnames(PCorMat12to9IAA)<-colnames(Correlation12to9IAA)
rownames(PCorMat12to9IAA)<-rownames(Correlation12to9IAA)
dim(PCorMat12to9IAA)
min(PCorMat12to9IAA)
max(PCorMat12to9IAA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat12to9IAA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat12to9IAA))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat12to9IAA2<-PCorMat12to9IAA[rownames(PCorMat12to9IAA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat12to9IAA2)
PCorMat12to9IAASizeEN<-PCorMat12to9IAA2[,colnames(PCorMat12to9IAA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat12to9IAASizeEN)
PCorMat12to9IAASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat12to9IAASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat12to9IAASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat12to9IAASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat12to9IAASizeEN_melted)
head(PCorMat12to9IAASizeEN_melted)
colnames(PCorMat12to9IAASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat12to9IAASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat12to9IAASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm12m9IAA[1:10,1:10]
dim(Deltaallvarsm12m9IAA)
vectitodeGE<-colnames(Deltaallvarsm12m9IAA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm12m9IAA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm12m9IAA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable12to9IAA<-PCorMat12to9IAASizeEN_melted[,-3]

NodeTable12to9IAA2<- data.frame(Var=c(as.vector(NodeTable12to9IAA[,1]),as.vector(NodeTable12to9IAA[,2])))
length(as.vector(NodeTable12to9IAA[,1])); length(as.vector(NodeTable12to9IAA[,2]));dim(NodeTable12to9IAA2)
head(NodeTable12to9IAA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable12to9IAA2$Omic<-""
head(NodeTable12to9IAA2)
NodeTable12to9IAA2[1:50,]

NodeTable12to9IAA2$Omic [NodeTable12to9IAA2$Var %in% vectitodeGE] <- 1
NodeTable12to9IAA2$Omic [NodeTable12to9IAA2$Var %in% vectitodeMetab] <- 2
NodeTable12to9IAA2$Omic [NodeTable12to9IAA2$Var %in% vectitodeDB] <- 3
table(NodeTable12to9IAA2$Omic)
#sum(table(NodeTable12to9IAA2$Omic))
NodeTable12to9IAA2[850:1065,]
dim(NodeTable12to9IAA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat12to9IAASizeEN_melted$Correlations)>=thrs) # 1753 true
dim(PCorMat12to9IAASizeEN_melted)
PCorMat12to9IAASizeEN_meltedCO0.7<-PCorMat12to9IAASizeEN_melted[abs(PCorMat12to9IAASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat12to9IAASizeEN_meltedCO0.7)
head(PCorMat12to9IAASizeEN_meltedCO0.7)
PCorMat12to9IAASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat12to9IAASizeEN_meltedCO0.7[,1]),as.vector(PCorMat12to9IAASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat12to9IAASizeEN_meltedCO0.7[,2]));dim(PCorMat12to9IAASizeEN_meltedCO0.72)
dim(unique(PCorMat12to9IAASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat12to9IAASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable12to9IAA2$Weight2<-tulo$Freq [match(NodeTable12to9IAA2$Var,tulo$Var1)]
NodeTable12to9IAA2[1:500,]
dim(NodeTable12to9IAA2)

NodeTable12to9IAA2[NodeTable12to9IAA2$Var=="PLEKHA8",]
dim(NodeTable12to9IAA2[NodeTable12to9IAA2$Var=="PLEKHA8",])

NodeTable12to9IAA2[NodeTable12to9IAA2$Var=="unknown272",]
dim(NodeTable12to9IAA2[NodeTable12to9IAA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable12to9IAA2$Weight2[is.na(NodeTable12to9IAA2$Weight2)]<-0
NodeTable12to9IAA2$Weight<-NodeTable12to9IAA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable12to9IAA2$Nameito<- abbreviate(NodeTable12to9IAA2$Var, 7, method = "both")
NodeTable12to9IAA2[NodeTable12to9IAA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm12m9IAA[1:10,1:10]
dim(Deltaallvarsm12m9IAA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm12m9IAA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable12to9IAA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable12to9IAA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable12to9IAA<-unique (NodeTable12to9IAA2)
length(unique (NodeTable12to9IAA$Var)) # 315
dim(NodeTable12to9IAA)
#write.table(NodeTable12to9IAA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable12to9IAA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -9 to -6 #####
##############################
##############################
###          GADA         ####
##############################

Correlation9to6GADA <- cor(Deltaallvarsm9m6GADA, use = "pairwise")  

class(Correlation9to6GADA)
isSymmetric(Correlation9to6GADA)
PCorMat9to6GADA <- cor2pcor(Correlation9to6GADA)
class(PCorMat9to6GADA)
isSymmetric(PCorMat9to6GADA)
PCorMat9to6GADA<-round(PCorMat9to6GADA,digits = 6)
isSymmetric(PCorMat9to6GADA)
colnames(PCorMat9to6GADA)<-colnames(Correlation9to6GADA)
rownames(PCorMat9to6GADA)<-rownames(Correlation9to6GADA)
dim(PCorMat9to6GADA)
min(PCorMat9to6GADA)
max(PCorMat9to6GADA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to6GADA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to6GADA))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat9to6GADA2<-PCorMat9to6GADA[rownames(PCorMat9to6GADA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat9to6GADA2)
PCorMat9to6GADASizeEN<-PCorMat9to6GADA2[,colnames(PCorMat9to6GADA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat9to6GADASizeEN)
PCorMat9to6GADASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat9to6GADASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat9to6GADASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat9to6GADASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat9to6GADASizeEN_melted)
head(PCorMat9to6GADASizeEN_melted)
colnames(PCorMat9to6GADASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat9to6GADASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat9to6GADASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm9m6GADA[1:10,1:10]
dim(Deltaallvarsm9m6GADA)
vectitodeGE<-colnames(Deltaallvarsm9m6GADA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm9m6GADA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm9m6GADA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable9to6GADA<-PCorMat9to6GADASizeEN_melted[,-3]

NodeTable9to6GADA2<- data.frame(Var=c(as.vector(NodeTable9to6GADA[,1]),as.vector(NodeTable9to6GADA[,2])))
length(as.vector(NodeTable9to6GADA[,1])); length(as.vector(NodeTable9to6GADA[,2]));dim(NodeTable9to6GADA2)
head(NodeTable9to6GADA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable9to6GADA2$Omic<-""
head(NodeTable9to6GADA2)
NodeTable9to6GADA2[1:50,]

NodeTable9to6GADA2$Omic [NodeTable9to6GADA2$Var %in% vectitodeGE] <- 1
NodeTable9to6GADA2$Omic [NodeTable9to6GADA2$Var %in% vectitodeMetab] <- 2
NodeTable9to6GADA2$Omic [NodeTable9to6GADA2$Var %in% vectitodeDB] <- 3
table(NodeTable9to6GADA2$Omic)
#sum(table(NodeTable9to6GADA2$Omic))
NodeTable9to6GADA2[850:1065,]
dim(NodeTable9to6GADA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat9to6GADASizeEN_melted$Correlations)>=thrs) # 110 true
dim(PCorMat9to6GADASizeEN_melted)
PCorMat9to6GADASizeEN_meltedCO0.7<-PCorMat9to6GADASizeEN_melted[abs(PCorMat9to6GADASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat9to6GADASizeEN_meltedCO0.7)
head(PCorMat9to6GADASizeEN_meltedCO0.7)
PCorMat9to6GADASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat9to6GADASizeEN_meltedCO0.7[,1]),as.vector(PCorMat9to6GADASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat9to6GADASizeEN_meltedCO0.7[,2]));dim(PCorMat9to6GADASizeEN_meltedCO0.72)
dim(unique(PCorMat9to6GADASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat9to6GADASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable9to6GADA2$Weight2<-tulo$Freq [match(NodeTable9to6GADA2$Var,tulo$Var1)]
NodeTable9to6GADA2[1:500,]
dim(NodeTable9to6GADA2)

NodeTable9to6GADA2[NodeTable9to6GADA2$Var=="PLEKHA8",]
dim(NodeTable9to6GADA2[NodeTable9to6GADA2$Var=="PLEKHA8",])

NodeTable9to6GADA2[NodeTable9to6GADA2$Var=="unknown272",]
dim(NodeTable9to6GADA2[NodeTable9to6GADA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable9to6GADA2$Weight2[is.na(NodeTable9to6GADA2$Weight2)]<-0
NodeTable9to6GADA2$Weight<-NodeTable9to6GADA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable9to6GADA2$Nameito<- abbreviate(NodeTable9to6GADA2$Var, 7, method = "both")
NodeTable9to6GADA2[NodeTable9to6GADA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm9m6GADA[1:10,1:10]
dim(Deltaallvarsm9m6GADA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm9m6GADA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)

NodeTable9to6GADA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable9to6GADA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable9to6GADA<-unique (NodeTable9to6GADA2)
length(unique (NodeTable9to6GADA2$Var)) # 315
dim(NodeTable9to6GADA) # 315
#write.table(NodeTable9to6GADA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable9to6GADA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
#### Time Period -9 to -6 ####
##############################
##############################
###           IAA         ####
##############################

Correlation9to6IAA <- cor(Deltaallvarsm9m6IAA, use = "pairwise")  

class(Correlation9to6IAA)
isSymmetric(Correlation9to6IAA)
PCorMat9to6IAA <- cor2pcor(Correlation9to6IAA)
class(PCorMat9to6IAA)
isSymmetric(PCorMat9to6IAA)
PCorMat9to6IAA<-round(PCorMat9to6IAA,digits = 5)
isSymmetric(PCorMat9to6IAA)
colnames(PCorMat9to6IAA)<-colnames(Correlation9to6IAA)
rownames(PCorMat9to6IAA)<-rownames(Correlation9to6IAA)
dim(PCorMat9to6IAA)
min(PCorMat9to6IAA)
max(PCorMat9to6IAA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat9to6IAA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat9to6IAA))
#####PAVPAVPAVPAPAVPAVPAVPAVPAV
################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat9to6IAA2<-PCorMat9to6IAA[rownames(PCorMat9to6IAA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat9to6IAA2)
PCorMat9to6IAASizeEN<-PCorMat9to6IAA2[,colnames(PCorMat9to6IAA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat9to6IAASizeEN)
PCorMat9to6IAASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat9to6IAASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat9to6IAASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat9to6IAASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat9to6IAASizeEN_melted)
head(PCorMat9to6IAASizeEN_melted)
colnames(PCorMat9to6IAASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat9to6IAASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat9to6IAASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm9m6IAA[1:10,1:10]
dim(Deltaallvarsm9m6IAA)
vectitodeGE<-colnames(Deltaallvarsm9m6IAA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm9m6IAA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm9m6IAA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable9to6IAA<-PCorMat9to6IAASizeEN_melted[,-3]

NodeTable9to6IAA2<- data.frame(Var=c(as.vector(NodeTable9to6IAA[,1]),as.vector(NodeTable9to6IAA[,2])))
length(as.vector(NodeTable9to6IAA[,1])); length(as.vector(NodeTable9to6IAA[,2]));dim(NodeTable9to6IAA2)
head(NodeTable9to6IAA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable9to6IAA2$Omic<-""
head(NodeTable9to6IAA2)
NodeTable9to6IAA2[1:50,]

NodeTable9to6IAA2$Omic [NodeTable9to6IAA2$Var %in% vectitodeGE] <- 1
NodeTable9to6IAA2$Omic [NodeTable9to6IAA2$Var %in% vectitodeMetab] <- 2
NodeTable9to6IAA2$Omic [NodeTable9to6IAA2$Var %in% vectitodeDB] <- 3
table(NodeTable9to6IAA2$Omic)
#sum(table(NodeTable9to6IAA2$Omic))
NodeTable9to6IAA2[850:1065,]
dim(NodeTable9to6IAA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat9to6IAASizeEN_melted$Correlations)>=thrs) # 81 true
dim(PCorMat9to6IAASizeEN_melted)
PCorMat9to6IAASizeEN_meltedCO0.7<-PCorMat9to6IAASizeEN_melted[abs(PCorMat9to6IAASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat9to6IAASizeEN_meltedCO0.7)
head(PCorMat9to6IAASizeEN_meltedCO0.7)
PCorMat9to6IAASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat9to6IAASizeEN_meltedCO0.7[,1]),as.vector(PCorMat9to6IAASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat9to6IAASizeEN_meltedCO0.7[,2]));dim(PCorMat9to6IAASizeEN_meltedCO0.72)
dim(unique(PCorMat9to6IAASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat9to6IAASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable9to6IAA2$Weight2<-tulo$Freq [match(NodeTable9to6IAA2$Var,tulo$Var1)]
NodeTable9to6IAA2[1:500,]
dim(NodeTable9to6IAA2)

NodeTable9to6IAA2[NodeTable9to6IAA2$Var=="PLEKHA8",]
dim(NodeTable9to6IAA2[NodeTable9to6IAA2$Var=="PLEKHA8",])

NodeTable9to6IAA2[NodeTable9to6IAA2$Var=="unknown272",]
dim(NodeTable9to6IAA2[NodeTable9to6IAA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable9to6IAA2$Weight2[is.na(NodeTable9to6IAA2$Weight2)]<-0
NodeTable9to6IAA2$Weight<-NodeTable9to6IAA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable9to6IAA2$Nameito<- abbreviate(NodeTable9to6IAA2$Var, 7, method = "both")
NodeTable9to6IAA2[NodeTable9to6IAA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm9m6IAA[1:10,1:10]
dim(Deltaallvarsm9m6IAA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm9m6IAA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable9to6IAA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable9to6IAA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable9to6IAA<-unique (NodeTable9to6IAA2)
length(unique (NodeTable9to6IAA$Var)) # 315
dim(NodeTable9to6IAA)
#write.table(NodeTable9to6IAA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable9to6IAA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->
##############################
### Time Period -6 to -3 #####
##############################
##############################
###          GADA         ####
##############################

Correlation6to3GADA <- cor(Deltaallvarsm6m3GADA, use = "pairwise")  

class(Correlation6to3GADA)
isSymmetric(Correlation6to3GADA)
PCorMat6to3GADA <- cor2pcor(Correlation6to3GADA)
class(PCorMat6to3GADA)
isSymmetric(PCorMat6to3GADA)
PCorMat6to3GADA<-round(PCorMat6to3GADA,digits = 6)
isSymmetric(PCorMat6to3GADA)
colnames(PCorMat6to3GADA)<-colnames(Correlation6to3GADA)
rownames(PCorMat6to3GADA)<-rownames(Correlation6to3GADA)
dim(PCorMat6to3GADA)
min(PCorMat6to3GADA)
max(PCorMat6to3GADA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat6to3GADA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat6to3GADA))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat6to3GADA2<-PCorMat6to3GADA[rownames(PCorMat6to3GADA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat6to3GADA2)
PCorMat6to3GADASizeEN<-PCorMat6to3GADA2[,colnames(PCorMat6to3GADA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat6to3GADASizeEN)
PCorMat6to3GADASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat6to3GADASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat6to3GADASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat6to3GADASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat6to3GADASizeEN_melted)
head(PCorMat6to3GADASizeEN_melted)
colnames(PCorMat6to3GADASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat6to3GADASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat6to3GADASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm6m3GADA[1:10,1:10]
dim(Deltaallvarsm6m3GADA)
vectitodeGE<-colnames(Deltaallvarsm6m3GADA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm6m3GADA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm6m3GADA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable6to3GADA<-PCorMat6to3GADASizeEN_melted[,-3]

NodeTable6to3GADA2<- data.frame(Var=c(as.vector(NodeTable6to3GADA[,1]),as.vector(NodeTable6to3GADA[,2])))
length(as.vector(NodeTable6to3GADA[,1])); length(as.vector(NodeTable6to3GADA[,2]));dim(NodeTable6to3GADA2)
head(NodeTable6to3GADA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable6to3GADA2$Omic<-""
head(NodeTable6to3GADA2)
NodeTable6to3GADA2[1:50,]

NodeTable6to3GADA2$Omic [NodeTable6to3GADA2$Var %in% vectitodeGE] <- 1
NodeTable6to3GADA2$Omic [NodeTable6to3GADA2$Var %in% vectitodeMetab] <- 2
NodeTable6to3GADA2$Omic [NodeTable6to3GADA2$Var %in% vectitodeDB] <- 3
table(NodeTable6to3GADA2$Omic)
#sum(table(NodeTable6to3GADA2$Omic))
NodeTable6to3GADA2[850:1065,]
dim(NodeTable6to3GADA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat6to3GADASizeEN_melted$Correlations)>=thrs) # 336 true
dim(PCorMat6to3GADASizeEN_melted)
PCorMat6to3GADASizeEN_meltedCO0.7<-PCorMat6to3GADASizeEN_melted[abs(PCorMat6to3GADASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat6to3GADASizeEN_meltedCO0.7)
head(PCorMat6to3GADASizeEN_meltedCO0.7)
PCorMat6to3GADASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat6to3GADASizeEN_meltedCO0.7[,1]),as.vector(PCorMat6to3GADASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat6to3GADASizeEN_meltedCO0.7[,2]));dim(PCorMat6to3GADASizeEN_meltedCO0.72)
dim(unique(PCorMat6to3GADASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat6to3GADASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable6to3GADA2$Weight2<-tulo$Freq [match(NodeTable6to3GADA2$Var,tulo$Var1)]
NodeTable6to3GADA2[1:500,]
dim(NodeTable6to3GADA2)

NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="PLEKHA8",]
dim(NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="PLEKHA8",])

NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="unknown272",]
dim(NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="unknown272",])

NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="sphingomyelin35",]
dim(NodeTable6to3GADA2[NodeTable6to3GADA2$Var=="sphingomyelin35",])

##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable6to3GADA2$Weight2[is.na(NodeTable6to3GADA2$Weight2)]<-0
NodeTable6to3GADA2$Weight<-NodeTable6to3GADA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable6to3GADA2$Nameito<- abbreviate(NodeTable6to3GADA2$Var, 7, method = "both")
NodeTable6to3GADA2[NodeTable6to3GADA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm6m3GADA[1:10,1:10]
dim(Deltaallvarsm6m3GADA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm6m3GADA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)

NodeTable6to3GADA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable6to3GADA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable6to3GADA<-unique (NodeTable6to3GADA2)
length(unique (NodeTable6to3GADA2$Var)) # 315
dim(NodeTable6to3GADA) # 315
#write.table(NodeTable6to3GADA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable6to3GADA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
#### Time Period -6 to -3 ####
##############################
##############################
###           IAA         ####
##############################

Correlation6to3IAA <- cor(Deltaallvarsm6m3IAA, use = "pairwise")  

class(Correlation6to3IAA)
isSymmetric(Correlation6to3IAA)
PCorMat6to3IAA <- cor2pcor(Correlation6to3IAA)
class(PCorMat6to3IAA)
isSymmetric(PCorMat6to3IAA)
PCorMat6to3IAA<-round(PCorMat6to3IAA,digits = 5)
isSymmetric(PCorMat6to3IAA)
colnames(PCorMat6to3IAA)<-colnames(Correlation6to3IAA)
rownames(PCorMat6to3IAA)<-rownames(Correlation6to3IAA)
dim(PCorMat6to3IAA)
min(PCorMat6to3IAA)
max(PCorMat6to3IAA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat6to3IAA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat6to3IAA))
#####PAVPAVPAVPAPAVPAVPAVPAVPAV
################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat6to3IAA2<-PCorMat6to3IAA[rownames(PCorMat6to3IAA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat6to3IAA2)
PCorMat6to3IAASizeEN<-PCorMat6to3IAA2[,colnames(PCorMat6to3IAA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat6to3IAASizeEN)
PCorMat6to3IAASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat6to3IAASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat6to3IAASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat6to3IAASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat6to3IAASizeEN_melted)
head(PCorMat6to3IAASizeEN_melted)
colnames(PCorMat6to3IAASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat6to3IAASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat6to3IAASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm6m3IAA[1:10,1:10]
dim(Deltaallvarsm6m3IAA)
vectitodeGE<-colnames(Deltaallvarsm6m3IAA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm6m3IAA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm6m3IAA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable6to3IAA<-PCorMat6to3IAASizeEN_melted[,-3]

NodeTable6to3IAA2<- data.frame(Var=c(as.vector(NodeTable6to3IAA[,1]),as.vector(NodeTable6to3IAA[,2])))
length(as.vector(NodeTable6to3IAA[,1])); length(as.vector(NodeTable6to3IAA[,2]));dim(NodeTable6to3IAA2)
head(NodeTable6to3IAA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable6to3IAA2$Omic<-""
head(NodeTable6to3IAA2)
NodeTable6to3IAA2[1:50,]

NodeTable6to3IAA2$Omic [NodeTable6to3IAA2$Var %in% vectitodeGE] <- 1
NodeTable6to3IAA2$Omic [NodeTable6to3IAA2$Var %in% vectitodeMetab] <- 2
NodeTable6to3IAA2$Omic [NodeTable6to3IAA2$Var %in% vectitodeDB] <- 3
table(NodeTable6to3IAA2$Omic)
#sum(table(NodeTable6to3IAA2$Omic))
NodeTable6to3IAA2[850:1065,]
dim(NodeTable6to3IAA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat6to3IAASizeEN_melted$Correlations)>=thrs) # 23 true
dim(PCorMat6to3IAASizeEN_melted)
PCorMat6to3IAASizeEN_meltedCO0.7<-PCorMat6to3IAASizeEN_melted[abs(PCorMat6to3IAASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat6to3IAASizeEN_meltedCO0.7)
head(PCorMat6to3IAASizeEN_meltedCO0.7)
PCorMat6to3IAASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat6to3IAASizeEN_meltedCO0.7[,1]),as.vector(PCorMat6to3IAASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat6to3IAASizeEN_meltedCO0.7[,2]));dim(PCorMat6to3IAASizeEN_meltedCO0.72)
dim(unique(PCorMat6to3IAASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat6to3IAASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable6to3IAA2$Weight2<-tulo$Freq [match(NodeTable6to3IAA2$Var,tulo$Var1)]
NodeTable6to3IAA2[1:500,]
dim(NodeTable6to3IAA2)

NodeTable6to3IAA2[NodeTable6to3IAA2$Var=="PLEKHA8",]
dim(NodeTable6to3IAA2[NodeTable6to3IAA2$Var=="PLEKHA8",])

NodeTable6to3IAA2[NodeTable6to3IAA2$Var=="unknown272",]
dim(NodeTable6to3IAA2[NodeTable6to3IAA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable6to3IAA2$Weight2[is.na(NodeTable6to3IAA2$Weight2)]<-0
NodeTable6to3IAA2$Weight<-NodeTable6to3IAA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable6to3IAA2$Nameito<- abbreviate(NodeTable6to3IAA2$Var, 7, method = "both")
NodeTable6to3IAA2[NodeTable6to3IAA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm6m3IAA[1:10,1:10]
dim(Deltaallvarsm6m3IAA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm6m3IAA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable6to3IAA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable6to3IAA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable6to3IAA<-unique (NodeTable6to3IAA2)
length(unique (NodeTable6to3IAA$Var)) # 315
dim(NodeTable6to3IAA)
#write.table(NodeTable6to3IAA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable6to3IAA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
### Time Period -3 to  0 #####
##############################
##############################
###          GADA         ####
##############################

Correlation3to0GADA <- cor(Deltaallvarsm3m0GADA, use = "pairwise")  

class(Correlation3to0GADA)
isSymmetric(Correlation3to0GADA)
PCorMat3to0GADA <- cor2pcor(Correlation3to0GADA)
class(PCorMat3to0GADA)
isSymmetric(PCorMat3to0GADA)
PCorMat3to0GADA<-round(PCorMat3to0GADA,digits = 6)
isSymmetric(PCorMat3to0GADA)
colnames(PCorMat3to0GADA)<-colnames(Correlation3to0GADA)
rownames(PCorMat3to0GADA)<-rownames(Correlation3to0GADA)
dim(PCorMat3to0GADA)
min(PCorMat3to0GADA)
max(PCorMat3to0GADA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat3to0GADA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat3to0GADA))

################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat3to0GADA2<-PCorMat3to0GADA[rownames(PCorMat3to0GADA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat3to0GADA2)
PCorMat3to0GADASizeEN<-PCorMat3to0GADA2[,colnames(PCorMat3to0GADA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat3to0GADASizeEN)
PCorMat3to0GADASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat3to0GADASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat3to0GADASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat3to0GADASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat3to0GADASizeEN_melted)
head(PCorMat3to0GADASizeEN_melted)
colnames(PCorMat3to0GADASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat3to0GADASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat3to0GADASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm3m0GADA[1:10,1:10]
dim(Deltaallvarsm3m0GADA)
vectitodeGE<-colnames(Deltaallvarsm3m0GADA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm3m0GADA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm3m0GADA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable3to0GADA<-PCorMat3to0GADASizeEN_melted[,-3]

NodeTable3to0GADA2<- data.frame(Var=c(as.vector(NodeTable3to0GADA[,1]),as.vector(NodeTable3to0GADA[,2])))
length(as.vector(NodeTable3to0GADA[,1])); length(as.vector(NodeTable3to0GADA[,2]));dim(NodeTable3to0GADA2)
head(NodeTable3to0GADA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable3to0GADA2$Omic<-""
head(NodeTable3to0GADA2)
NodeTable3to0GADA2[1:50,]

NodeTable3to0GADA2$Omic [NodeTable3to0GADA2$Var %in% vectitodeGE] <- 1
NodeTable3to0GADA2$Omic [NodeTable3to0GADA2$Var %in% vectitodeMetab] <- 2
NodeTable3to0GADA2$Omic [NodeTable3to0GADA2$Var %in% vectitodeDB] <- 3
table(NodeTable3to0GADA2$Omic)
#sum(table(NodeTable3to0GADA2$Omic))
NodeTable3to0GADA2[850:1065,]
dim(NodeTable3to0GADA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat3to0GADASizeEN_melted$Correlations)>=thrs) # 131 true
dim(PCorMat3to0GADASizeEN_melted)
PCorMat3to0GADASizeEN_meltedCO0.7<-PCorMat3to0GADASizeEN_melted[abs(PCorMat3to0GADASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat3to0GADASizeEN_meltedCO0.7)
head(PCorMat3to0GADASizeEN_meltedCO0.7)
PCorMat3to0GADASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat3to0GADASizeEN_meltedCO0.7[,1]),as.vector(PCorMat3to0GADASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat3to0GADASizeEN_meltedCO0.7[,2]));dim(PCorMat3to0GADASizeEN_meltedCO0.72)
dim(unique(PCorMat3to0GADASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat3to0GADASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable3to0GADA2$Weight2<-tulo$Freq [match(NodeTable3to0GADA2$Var,tulo$Var1)]
NodeTable3to0GADA2[1:500,]
dim(NodeTable3to0GADA2)

NodeTable3to0GADA2[NodeTable3to0GADA2$Var=="PLEKHA8",]
dim(NodeTable3to0GADA2[NodeTable3to0GADA2$Var=="PLEKHA8",])

NodeTable3to0GADA2[NodeTable3to0GADA2$Var=="unknown272",]
dim(NodeTable3to0GADA2[NodeTable3to0GADA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable3to0GADA2$Weight2[is.na(NodeTable3to0GADA2$Weight2)]<-0
NodeTable3to0GADA2$Weight<-NodeTable3to0GADA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable3to0GADA2$Nameito<- abbreviate(NodeTable3to0GADA2$Var, 7, method = "both")
NodeTable3to0GADA2[NodeTable3to0GADA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm3m0GADA[1:10,1:10]
dim(Deltaallvarsm3m0GADA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm3m0GADA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)

NodeTable3to0GADA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable3to0GADA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable3to0GADA<-unique (NodeTable3to0GADA2)
length(unique (NodeTable3to0GADA2$Var)) # 315
dim(NodeTable3to0GADA) # 315
#write.table(NodeTable3to0GADA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable3to0GADA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

##############################
#### Time Period -3 to  0 ####
##############################
##############################
###           IAA         ####
##############################

Correlation3to0IAA <- cor(Deltaallvarsm3m0IAA, use = "pairwise")  

class(Correlation3to0IAA)
isSymmetric(Correlation3to0IAA)
PCorMat3to0IAA <- cor2pcor(Correlation3to0IAA)
class(PCorMat3to0IAA)
isSymmetric(PCorMat3to0IAA)
PCorMat3to0IAA<-round(PCorMat3to0IAA,digits = 5)
isSymmetric(PCorMat3to0IAA)
colnames(PCorMat3to0IAA)<-colnames(Correlation3to0IAA)
rownames(PCorMat3to0IAA)<-rownames(Correlation3to0IAA)
dim(PCorMat3to0IAA)
min(PCorMat3to0IAA)
max(PCorMat3to0IAA)
########################################################
# Genes= 862, Metabolites= 245, DB= 3
colgroup<-c(rep("dodgerblue2",862),rep("gold2",245), rep("deeppink1",3))
length(colgroup)==length(colnames(PCorMat3to0IAA))
namegroup<-c(rep("Genes",862),rep("Metabolites",245), rep("deeppink1",3))
length(namegroup)==length(colnames(PCorMat3to0IAA))
#####PAVPAVPAVPAPAVPAVPAVPAVPAV
################################################
#  Subsetting the ParCor matrix to the mentioned actors
dim(ENSelVarswithDBGADAIAA) # 315

PCorMat3to0IAA2<-PCorMat3to0IAA[rownames(PCorMat3to0IAA) %in% ENSelVarswithDBGADAIAA$x,]
dim(PCorMat3to0IAA2)
PCorMat3to0IAASizeEN<-PCorMat3to0IAA2[,colnames(PCorMat3to0IAA2) %in% ENSelVarswithDBGADAIAA$x]
dim(PCorMat3to0IAASizeEN)
PCorMat3to0IAASizeEN[1:10,1:10]
################################################
# Melting the matrix
upperTriangle = upper.tri(PCorMat3to0IAASizeEN, diag=F) #turn into a upper triangle
correlations.upperTriangle = PCorMat3to0IAASizeEN #take a copy of the original cor-mat
correlations.upperTriangle[!upperTriangle]<-NA #set everything not in upper triangle o NA
PCorMat3to0IAASizeEN_melted<-na.omit(melt(correlations.upperTriangle, value.name ="Correlation")) 
dim(PCorMat3to0IAASizeEN_melted)
head(PCorMat3to0IAASizeEN_melted)
colnames(PCorMat3to0IAASizeEN_melted)<- c("Var1", "Var2","Correlations")
#write.table(PCorMat3to0IAASizeEN_melted,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/PCorMat3to0IAASizeEN_melted.txt", sep="\\t")

##########################################################################################
# Three vectors to assign the category GE (1),MEtab (2) or DB (3)
Deltaallvarsm3m0IAA[1:10,1:10]
dim(Deltaallvarsm3m0IAA)
vectitodeGE<-colnames(Deltaallvarsm3m0IAA)[1:862]
vectitodeMetab<-colnames(Deltaallvarsm3m0IAA)[863:1107]
vectitodeDB<-colnames(Deltaallvarsm3m0IAA)[1108:1110]
#########################################################################
# Creation of the Node Table

NodeTable3to0IAA<-PCorMat3to0IAASizeEN_melted[,-3]

NodeTable3to0IAA2<- data.frame(Var=c(as.vector(NodeTable3to0IAA[,1]),as.vector(NodeTable3to0IAA[,2])))
length(as.vector(NodeTable3to0IAA[,1])); length(as.vector(NodeTable3to0IAA[,2]));dim(NodeTable3to0IAA2)
head(NodeTable3to0IAA2)
#########################################################################
# Creation of clasifier of what type of dataset, the feature is
NodeTable3to0IAA2$Omic<-""
head(NodeTable3to0IAA2)
NodeTable3to0IAA2[1:50,]

NodeTable3to0IAA2$Omic [NodeTable3to0IAA2$Var %in% vectitodeGE] <- 1
NodeTable3to0IAA2$Omic [NodeTable3to0IAA2$Var %in% vectitodeMetab] <- 2
NodeTable3to0IAA2$Omic [NodeTable3to0IAA2$Var %in% vectitodeDB] <- 3
table(NodeTable3to0IAA2$Omic)
#sum(table(NodeTable3to0IAA2$Omic))
NodeTable3to0IAA2[850:1065,]
dim(NodeTable3to0IAA2)

##########################################################################################
# Calculate the weight of each variable

thrs<-0.7
table(abs(PCorMat3to0IAASizeEN_melted$Correlations)>=thrs) # 171 true
dim(PCorMat3to0IAASizeEN_melted)
PCorMat3to0IAASizeEN_meltedCO0.7<-PCorMat3to0IAASizeEN_melted[abs(PCorMat3to0IAASizeEN_melted$Correlations)>=thrs,]
dim(PCorMat3to0IAASizeEN_meltedCO0.7)
head(PCorMat3to0IAASizeEN_meltedCO0.7)
PCorMat3to0IAASizeEN_meltedCO0.72<- data.frame(Var=c(as.vector(PCorMat3to0IAASizeEN_meltedCO0.7[,1]),as.vector(PCorMat3to0IAASizeEN_meltedCO0.7[,2])))
length(as.vector(PCorMat3to0IAASizeEN_meltedCO0.7[,2]));dim(PCorMat3to0IAASizeEN_meltedCO0.72)
dim(unique(PCorMat3to0IAASizeEN_meltedCO0.72))

tulo<-as.data.frame(sort(table(PCorMat3to0IAASizeEN_meltedCO0.72$Var),decreasing = TRUE))

NodeTable3to0IAA2$Weight2<-tulo$Freq [match(NodeTable3to0IAA2$Var,tulo$Var1)]
NodeTable3to0IAA2[1:500,]
dim(NodeTable3to0IAA2)

NodeTable3to0IAA2[NodeTable3to0IAA2$Var=="PLEKHA8",]
dim(NodeTable3to0IAA2[NodeTable3to0IAA2$Var=="PLEKHA8",])

NodeTable3to0IAA2[NodeTable3to0IAA2$Var=="unknown272",]
dim(NodeTable3to0IAA2[NodeTable3to0IAA2$Var=="unknown272",])
##########################################################################################
# Adjusting of Weights to remove NA's
NodeTable3to0IAA2$Weight2[is.na(NodeTable3to0IAA2$Weight2)]<-0
NodeTable3to0IAA2$Weight<-NodeTable3to0IAA2$Weight2 +1

##########################################################################################
# Creating a abbreviated name to plot
NodeTable3to0IAA2$Nameito<- abbreviate(NodeTable3to0IAA2$Var, 7, method = "both")
NodeTable3to0IAA2[NodeTable3to0IAA2$Omic==2,]
##########################################################################################
# Creating Expression Level column to know if the feature is up or downregulated
Deltaallvarsm3m0IAA[1:10,1:10]
dim(Deltaallvarsm3m0IAA)
vectitodebehaviormedio<-data.frame(ExpressionLevel=colMeans(Deltaallvarsm3m0IAA))
head(vectitodebehaviormedio)
dim(vectitodebehaviormedio)


NodeTable3to0IAA2$ExpressionLevel<-vectitodebehaviormedio$ExpressionLevel[match(NodeTable3to0IAA2$Var,rownames(vectitodebehaviormedio))]

##########################################################################################
# Definitive table #
NodeTable3to0IAA<-unique (NodeTable3to0IAA2)
length(unique (NodeTable3to0IAA$Var)) # 315
dim(NodeTable3to0IAA)
#write.table(NodeTable3to0IAA,file="/media/data/leobalzano/EnrichmentFaab/ParCorFaab/NodeTable3to0IAA.txt", sep="\\t")
##########################################################################################
# To Cytoscape
# 1) give a tab at the beginning of the first row (the colnames row)
# 2) Erase all "" in both datasets
# 3) Import PCorMat9to12caseswithDB_melted_CO0.7 as Network specifying...
#     3.1) In Network Collection: Create a new network collection
#     3.2) Then hit select none and then...
#     3.3) Assign Var 1 as the Source Node
#     3.4) Assign Var2 as the target Node
#     3.5) Assign Correlations as the Edge Attribute
# 4) To the node table add NodeTable9to12caseswithDB as follows:
#     4.1) File -> Import -> Table -> File -> NodeTable9to12caseswithDB
#     4.2) In Where to Import Table Data: Select To selected networks only
#     4.3) In Network List: Select PCorMat9to12caseswithDB_melted_CO0.7.txt
#     4.4) In Import Data as: Select Node Table Columns
#     4.5) Then hit select none and then...
#     4.6) Assign Var as the Key
#     4.7) Assign Omic, Weight and Nameito as Attribute
# 5) Go to Tools-> NetworkAnalyzer -> Network Analysis -> Analyze network->

