###########################################################
########   DatasetsasinputinPaintomicsxTEDDY.R     ########
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This script is to create the tables to include in paintomics for pathways analyzes

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()
###########################################################
# Functions:
TEDDYtools = "/media/data/leobalzano/ScriptsForTEDDY/Scripts/TEDDYtools2.R" # Location of TEDDYtools
source(TEDDYtools)

PaintomicsSubsetter<- function (X=LM_globalTargets, Y=LM_globalEXPMAT, Z= metabolomicstargets, ZZ=globalmetabolomicsMatrix,
                                ageg=2, FAAb=0, timepoints=3, timeframe=c(12,9,6,3,0) ) {
  if (ageg==0 & FAAb==0) {
    print("Both matrices are exactly the same, No Subset was created")
    targetsGE = X
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = Z
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg!=0 & FAAb==0) {
    targetsGE = filter(X,agegroup==ageg)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,AgeGroup==ageg)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg==0 & FAAb!=0) {
    targetsGE = filter(X,FirstAAb==FAAb)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,FirstAAbCC==FAAb)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
  }
  if (ageg!=0 & FAAb!=0) {
    print("Warning: Subsetting by two variables could reduce too much the number of Invividuals")
    targetsGE = filter(X,agegroup==ageg,FirstAAb==FAAb)
    datasetGE = Y[,as.character(targetsGE$sample_mask_id)]
    targetsMET = filter(Z,AgeGroup==ageg,FirstAAbCC==FAAb)
    datasetMET = ZZ[,as.character(targetsMET$Sample.Mask.Id)]
    
  }
  # Get Patients with at least 3 time points  (N=136)
  #GeneExpression
  patientTable = table(targetsGE$mask_id)
  patients = names(patientTable)[patientTable>=timepoints]
  #length(patients)
  Xsmall = filter(X,mask_id %in% patients)
  
  #Metabolomics
  patientTableMET = table(targetsMET$Mask.Id)
  patients = names(patientTableMET)[patientTableMET>=timepoints]
  #length(patients)
  Zsmall = filter(Z,Mask.Id %in% patients)
  
  
  # Paintomics Genes tables
  V2_pltDF = c()
  V2_pltDFMET = c()
  for(timepoint in timeframe){
    #GeneExpression
    testTargets = filter(Xsmall,time==timepoint)
    meansDFv2 = getGroupMean(groupcolumn = "outcome",stargets = testTargets,sexprmat = Y)
    LogRv2 = as.numeric(meansDFv2[,3]) - as.numeric(meansDFv2[,2])
    names(LogRv2) = rownames(meansDFv2)
    V2_pltDF = cbind(V2_pltDF, LogRv2)
    
    #Metabolomics
    testTargetsMET = filter(Zsmall,Sconv==timepoint*-1)
    meansDFv2MET = getGroupMeanMET(groupcolumn = "Outcome",stargets = testTargetsMET,sexprmat = ZZ)
    #head(meansDFv2MET)
    LogRv2MET = as.numeric(meansDFv2MET[,3]) - as.numeric(meansDFv2MET[,2])
    names(LogRv2MET) = rownames(meansDFv2MET)
    V2_pltDFMET = cbind(V2_pltDFMET, LogRv2MET)
    
  }
  columnnames=NULL
  for (i in 1:length(timeframe)) {
    names<-paste("Time",timeframe[i], sep="_-")
    columnnames<-c(columnnames,names)
  }
  colnames(V2_pltDF) = columnnames
  colnames(V2_pltDFMET) = columnnames
  result<-list(targets=targets,datasetGE=datasetGE,datasetMET=datasetMET, paintomicsGETable= V2_pltDF,paintomicsMETTable= V2_pltDFMET)
  return(result)
  
}

###########################################################
# Data:
wholedescriptivevars <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/wholedescriptivevars.txt")
Allmetabolomicstogether <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/Allmetabolomicstogether.txt")
AllmetabolitesConverted <- read.csv("/media/data/leobalzano/ScriptsForTEDDY/Data/AllmetabolitesConverted.txt", sep= "\\t")
load( '/media/data/leobalzano/ScriptsForTEDDY/Data/rawTEDDY_GEX_inversevst.ro' )


###########################################################
# Libraries:
#require("abind")
#library(gplots)
#library(data.table)
#library("gdata")
#require("VennDiagram")

###########################################################
#Retrieve Expression Data

# Global targets
LM_globalTargets = loadTargets_v2(SELdisease = "IA",time_point = c(12,9,6,3,0),complete_case = "flexible",targets = targets) #Uses all samples

# Process Expression Matrix
RAW_EXPMAT = log2(RAW_EXPMATv2)

# Annotation
LM_globalEXPMAT = getExprmat(RAWEXPRMAT = RAW_EXPMAT,STARGETS = LM_globalTargets, ciNORM = FALSE, gNORM = FALSE)

# Add information from FAAB
FAAb = na.omit(read.delim(file = "/media/data/leobalzano/ScriptsForTEDDY/Data/wholedescriptivevars.txt", sep="\\t",header = TRUE,stringsAsFactors = F)[,c("FirstAAbCC","sample_mask_id")])
colnames(FAAb) = c("FirstAAb","sample_mask_id")

LM_globalTargets = addFeatures(newDFdata = FAAb, TARGETS = LM_globalTargets, selCol = "sample_mask_id")

LM_globalTargets$FirstAAb[(is.na(LM_globalTargets$FirstAAb))] = 132
LM_globalEXPMAT = LM_globalEXPMAT[,as.character(LM_globalTargets$sample_mask_id)]

#If you are planning to subset by age group or any other annotation, please do it after annotation

###########################################################
# SUBSET HERE
dim(LM_globalTargets)# 742 x 20
LM_globalTargets[1:3,1:10]
LM_globalEXPMAT[1:3,1:10]
colnames(Allmetabolomicstogether[1:40])
Allmetabolomicstogether[1:10,1:10]
dim(Allmetabolomicstogether)
Allmetabolomicstogether[1:5,50:60]
rownames(Allmetabolomicstogether)<-Allmetabolomicstogether[,1]
metabolomicstargets<-Allmetabolomicstogether[,1:35]
dim(metabolomicstargets)
metabolomicstargets[1:10,1:10]
globalmetabolomicsMatrix<-t(Allmetabolomicstogether[,-(1:35)])
dim(globalmetabolomicsMatrix)
globalmetabolomicsMatrix[1:3,1:10]
globalmetabolomicsMatrix<-log2(globalmetabolomicsMatrix)
globalmetabolomicsMatrix[1:3,1:10]
namesrows<-gsub("._neglip", "!_neglip",rownames(globalmetabolomicsMatrix))
namesrows2<-gsub("._poslip", "!_poslip",namesrows)
rownames(globalmetabolomicsMatrix)<-namesrows2

###########################################################
# For AgeGroup:
test<-PaintomicsSubsetter(X=LM_globalTargets, Y=LM_globalEXPMAT,  Z= metabolomicstargets, ZZ=globalmetabolomicsMatrix,
                          ageg=0, FAAb=0,timepoints=3, timeframe=c(12,9,6,3,0) )

test$paintomicsGETable
test$paintomicsMETTable
test$targets
test$datasetGE
###########################################################
# Merging with AllmetabolitesConverted
test$paintomicsMETTable
AllmetabolitesConverted <- read.delim2("~/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/AllmetabolitesConverted.txt")
head(AllmetabolitesConverted)
rownames(AllmetabolitesConverted)<-AllmetabolitesConverted[,1];AllmetabolitesConverted<-AllmetabolitesConverted[-1]
log2ratMetabolsm12tom0originalnames<-merge(AllmetabolitesConverted,test$paintomicsMETTable, by=0)
dim(AllmetabolitesConverted);dim(test$paintomicsMETTable);dim(log2ratMetabolsm12tom0originalnames);
head(log2ratMetabolsm12tom0originalnames)
#write.table(log2ratMetabolsm12tom0originalnames, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/log2ratMetabolsm12tom0originalnames.txt", sep="\\t")

###########################################################
# Rename with FunCatsamename all dataset
#colnames(log2ratMetabolsm12tom0originalnames)
funcatsn<-log2ratMetabolsm12tom0originalnames[,c(5,7:11)]
head(funcatsn)

#  Make Unique based on the names
UniqueSamemetnamesm12tom0from1321<-aggregate(.~FunCatsameName,funcatsn,mean)
dim(UniqueSamemetnamesm12tom0from1321)
head(UniqueSamemetnamesm12tom0from1321)
rownames(UniqueSamemetnamesm12tom0from1321)<-UniqueSamemetnamesm12tom0from1321[,1]

UniqueSamemetnamesm12tom0from1321<-UniqueSamemetnamesm12tom0from1321[,-1]
#write.table(UniqueSamemetnamesm12tom0from1321, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/UniqueSamemetnamesm12tom0from1321.txt", sep="\\t")
#a<-grep("rol", UniqueSamemetnamesm12tom0from1321[,1], value=TRUE)
#fila<-UniqueSamemetnamesm12tom0from1321[is.element(UniqueSamemetnamesm12tom0from1321[,1],a),]
#fila

###########################################################
#####       So for Paintomics we have to create      ######
###########################################################
# 1) PaintomicsGlobal
# GeneExpression:
# test$paintomicsGETable
# 1) PaintGlobalGenes<- test$paintomicsGETable
#write.table(PaintGlobalGenes, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/Stratification/PaintGlobalGenes.txt", sep="\\t")

# 2) PaintDefinitiveGeneListNewData.txt # This is the list of 862 genes selected by NPLSDA

# Metabolomics:
# UniqueSamemetnamesm12tom0from1321
# 3) PaintGlobalMetabs<-UniqueSamemetnamesm12tom0from1321
#write.table(PaintGlobalMetabs, "/home/leobalzano/Documents/TEDDY data/MP149/Diet/NPLS/NPLSDAmetabolomicsDecember/Paintomics/Stratification/PaintGlobalMetabs.txt", sep="\\t")
# 4) PaintlistUniqueMetabolites.txt  # This is the list of Unique names of the selected 245 metabolites

# Use the four datasets as input in Paintomics
