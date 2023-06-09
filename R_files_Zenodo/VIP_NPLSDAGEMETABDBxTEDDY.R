###########################################################
############   VIP_NPLSDAGEMETABDBxTEDDY.R     ############
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# Tis script is to perform an NPLSDA for selected features from GE (862), Metabolomics(245)
# And Dietary Biomarkers (3)

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
setwd("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY") # Laptop
getwd()
###########################################################
# Functions:
#source ("/media/data/leobalzano/ScriptsForTEDDY/Scripts/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA
# Laptop:
source ("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Scripts/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA


###########################################################
# Data:
# Desktop:
# Gene Expression
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayGE862x136.RData")

# Metabolomics
load("/media/data/leobalzano/ScriptsForTEDDY/Data/Allmetabolomics136.RData")

# Dietary biomarkers
load("/media/data/leobalzano/ScriptsForTEDDY/Data/FullarrayDB.RData")

# Response Variable
load("/media/data/leobalzano/ScriptsForTEDDY/Data/Outcomedummyarray.RData")
load("/media/data/leobalzano/ScriptsForTEDDY/Data/outcomedummyarray136.RData")

# Misc
definitivemetabolitesSelectedConvertedmarch9 <- read.delim("/media/data/leobalzano/ScriptsForTEDDY/Data/definitivemetabolitesSelectedConvertedmarch9.txt", sep="\\t")

# Laptop
# Gene Expression
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/FullarrayGE862x136.RData")

# Metabolomics
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/Allmetabolomics136.RData")

# Dietary biomarkers
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/FullarrayDB.RData")

# Response Variable
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/Outcomedummyarray.RData")
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/outcomedummyarray136.RData")

# Misc
definitivemetabolitesSelectedConvertedmarch9 <- read.delim("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/definitivemetabolitesSelectedConvertedmarch9.txt", sep="\\t")



###########################################################
# Libraries:
library("outliers")
require("abind")

###########################################################
#     Changing the names of the selected metabolites     ##
AllmetabolomicsNewName136<-Allmetabolomics136
colnames(AllmetabolomicsNewName136)<- definitivemetabolitesSelectedConvertedmarch9[ match( colnames( AllmetabolomicsNewName136 ) ,
                                        definitivemetabolitesSelectedConvertedmarch9[ , "ID_Var" ]),
                                                                 "FunCat" ]
colnames(AllmetabolomicsNewName136)
###########################################################
#############     Merge GE and Metab     ##################
FullarrayGEMETABenX<-abind(FullarrayGE862x136,AllmetabolomicsNewName136, along=2)
#
dim(FullarrayGEMETABenX)
#save(FullarrayGEMETABenX, file="FullarrayGEMETABenX.RData")

###########################################################
# Putting DB in the same size as the other previous to merge
dim(FullarrayGEMETABenX)
dim(FullarrayDB)

FullarrayDB136<- FullarrayDB[rownames(FullarrayDB) %in% rownames(outcomedummyarray136),,]
dim(FullarrayDB136)
colnames(FullarrayDB136) # Retain just vitD (15),VitC (43) and alphaTocopherol(16)
PieceFullarrayDB136<-FullarrayDB136[,c(15,16,43),]
dim(PieceFullarrayDB136)
PieceFullarrayDB136[1:10,,1]
colnames(PieceFullarrayDB136)
#############################################################################
# Merge with GE and METAB
dim(FullarrayGEMETABenX)
FullarrayGEMETABenX[1:10,1:5,1]
dim(PieceFullarrayDB136)

Finalarray<- abind (FullarrayGEMETABenX,PieceFullarrayDB136, along=2);dim(FullarrayGEMETABenX);dim(PieceFullarrayDB136);dim(Finalarray)
anyNA(Finalarray)
#############################################################################
# 1) Calculate the bestfittedmodel

model<-bestfittedmodel (X=Finalarray,centering=0) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
model 
# Result: The best model is 4,2,3 according to the analysis.
#############################################################################
# 2) Calculate the bestfittedmodel
NPLSDAGEMETABDBBestModel<-NPLSDAmod(XN=Finalarray,YN= outcomedummyarray136, outcome.Y=NULL, COMP= c(4,2,3), factors=2, centering=0) # This is the same as the previous one

summary(NPLSDAGEMETABDBBestModel)
NPLSDAGEMETABDBBestModel$FactorsX
NPLSDAGEMETABDBBestModel$NPLSDAQ22D
NPLSDAGEMETABDBBestModel$NPLSDAQ22Dcum
NPLSDAGEMETABDBBestModel$NPLSDAQ2mean3D
NPLSDAGEMETABDBBestModel$VIP
NPLSDAGEMETABDBBestModel$Ypred
NPLSDAGEMETABDBBestModel$residuals
NPLSDAGEMETABDBBestModelv2020<-NPLSDAGEMETABDBBestModel
#save(NPLSDAGEMETABDBBestModelv2020, file= "/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/NPLSDAGEMETABDBBestModelv2020.RData")
#load("/media/data/leobalzano/ScriptsForTEDDY/Data/NPLSDAGEMETABDBBestModel.RData")
load("/Users/leobalzano/Dropbox (UFL)/ScriptsForTEDDY/Data/NPLSDAGEMETABDBBestModelv2020.RData")

#############################################################################
# 3) Plotting NPLSDA
#pdf("NPLSDAGEMETABDBBestModel.pdf")
ploteoNPLSDAGEMETABDBBestModel<- plotNPLSDAmod (X=NPLSDAGEMETABDBBestModel, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                                cutoff = 20, factors=2, penalty=1) 
#dev.off()


