###########################################################
###############     DEMONPLSDAxTEDDY.R     ################
###########################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# This is a demo script for NPLSDA analysis
# In this script we show how to use the NPLSDAfunctionsApr11.R
# You can evaluate VIP selection strategy 

###########################################################
setwd("/media/data/leobalzano/ScriptsForTEDDY")
getwd()
###########################################################
# Functions:
source ("/media/data/leobalzano/ScriptsForTEDDY/Scripts/NPLSDAfunctionsApr11.R") # These are the functions created to perform the NPLSDA

# Data: 
# Gene Expression
load("/media/data/leobalzano/ScriptsForTEDDY/Data/arraygeneexpressionsexample.RData")

# Response Variable
load("/media/data/leobalzano/ScriptsForTEDDY/Data/Outcomedummyarray.RData")

###########################################################
##########     Example with a piece of data     ###########
###########################################################
# 1.- Determining the best fitted model
dim(arraygeneexpressionsexample)
arraygeneexpressionsexample    # Note the presence of NA missing values.
model<-bestfittedmodel (X=arraygeneexpressionsexample,centering=2) # 0= No centering; 1= centering by Individuals; 2= centering by Variables;3= centering by Time
model
# If btm>CriVal, that is the best minimal fitted model, in the example, the best minimal fitted model would be 3,2,2

###########################################################
# 2.- Imputing data based on the model that gathers the maximum ammount of information
Fullarray<-Imputemethod(X=arraygeneexpressionsexample,fac=c(3, 2, 2), conver = 1e-07, max.iter = 1000)
# Compare before and after the imputation:
arraygeneexpressionsexample [1:10,1:10,1]
Fullarray[1:10,1:10,1]

###########################################################
# 3.- Tucker3 Analysis
Tuck3<-tucker3mod2(X=Fullarray,COMP=c(3, 2, 2), conver = 1e-07, max.iter = 1000)
summary (Tuck3)
Tuck3

###########################################################
# 4.- NPLSDA
# Best Fitted model was 3,3,2

# Here you can go with:
# 1) The no imputed data and using the adequate set of components to calculate
NPLSDAwithBestModel<-NPLSDAmod(XN=arraygeneexpressionsexample,YN= Outcomedummyarray, COMP=c(3, 2, 2),outcome.Y=NULL, factors=2, centering=0)

# 2) The no imputed data imputing a conventional 2x2x2 components
NPLSDAwithBestModel<-NPLSDAmod(XN=arraygeneexpressionsexample,YN= Outcomedummyarray, outcome.Y=NULL, factors=2, centering=0)

# 3) Use the data that you already impute in the example 2
NPLSDAwithBestModel<-NPLSDAmod(XN=Fullarray,YN= Outcomedummyarray, outcome.Y=NULL, factors=2, centering=0) # This is the same as the previous one

summary(NPLSDAwithBestModel)
NPLSDAwithBestModel$FactorsX
NPLSDAwithBestModel$VIP2D
NPLSDAwithBestModel$Ypred
NPLSDAwithBestModel$residuals

###########################################################
# 5.- Plotting NPLSDA
# Here you can go to the original source and you can notice that there are a lot more graphs turned off just for convenience
ploteoNPLSDAGENEXpiece<- plotNPLSDAmod (X=NPLSDAwithBestModel, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                        cutoff = 20, factors=2, penalty=1) 


###########################################################
###############     Variable Selection     ################
###########################################################
# Variable selection was performed by VIP calculated by three different manners
# Variable selection by VIP3Dmodel2
summary(NPLSDAwithBestModel)
NPLSDAwithBestModel$VIP3Dmodel2


vipsoutcomemet<-data.frame(NPLSDAwithBestModel$VIP3Dmodel2)

thrs<-99 # threshold applied
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

vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t4[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t5[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)

# Retain just these variables in Gene Expression 
dim(Fullarray)
PosLipselVars

FullarrayGenesVIPSelVars<-Fullarray[,is.element(colnames(Fullarray),
                                                  PosLipselVars),]
dim(FullarrayGenesVIPSelVars)

### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=Outcomedummyarray, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 



###########################################################
# Variable selection by VIP3Dmodel1
summary(NPLSDAwithBestModel)
NPLSDAwithBestModel$VIP3Dmodel1
dim(NPLSDAwithBestModel$VIP3Dmodel1)


### Start From here
vipsoutcomemet<-data.frame(NPLSDAwithBestModel$VIP3Dmodel1)
#apply(vipsoutcomemet, 2, function(x) is.numeric(x))

thrs<-95 # You can play with the threshold
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
#for (i in 1){      # Here you can play with the number of components 
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

# Retain just these variables in Gene Expression 2D
dim(Fullarray)
PosLipselVars



FullarrayGenesVIPSelVars<-Fullarray[,is.element(colnames(Fullarray),
                                                          PosLipselVars),]
dim(FullarrayGenesVIPSelVars)



### Now NPLSDA and Graph
NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=Outcomedummyarray, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 


###########################################################
# Gene Expression USING VIP2D
summary(NPLSDAwithBestModel)
NPLSDAwithBestModel$VIP2D


vipsoutcome2D<-data.frame(NPLSDAwithBestModel$VIP2D)

colp1<-paste(rownames(NPLSDAwithBestModel$VIP3Dmodel1),"-12", sep="_")
colp2<-paste(rownames(NPLSDAwithBestModel$VIP3Dmodel1),"-9", sep="_")
colp3<-paste(rownames(NPLSDAwithBestModel$VIP3Dmodel1),"-6", sep="_")
colp4<-paste(rownames(NPLSDAwithBestModel$VIP3Dmodel1),"-3", sep="_")
colp5<-paste(rownames(NPLSDAwithBestModel$VIP3Dmodel1),"0", sep="_")


rownames(vipsoutcome2D)<-c(colp1,colp2,colp3,colp4,colp5)

vipsoutcomemet<-vipsoutcome2D
thrs<-95    # Here you can play with the threshold
vipsoutcomemetsubset95<-NULL
temp<-NULL
vipsoutcomemetsubset<-list()
#for (i in 2){      # Here you can ply with the number of components
for (i in 1:2) {
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

vaya<-c (as.vector(rownames(vipsoutcomemetsubset$t1[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t2[[1]])),
         as.vector(rownames(vipsoutcomemetsubset$t3[[1]]))) 
duplicated(vaya)
PosLipselVars<-vaya[!duplicated(vaya)]
length(PosLipselVars)
string<-strsplit (PosLipselVars,"_")
string<-listTOdata.frame(string)
genelist<- string[,1]
duplicated (genelist)
genes<-genelist[!duplicated(genelist)]
length(genes)

# Retain just these variables in Gene Expression 2D
dim(Fullarray)
genes

FullarrayGenesVIPSelVars<-Fullarray[,is.element(colnames(Fullarray),
                                                  genes),]
dim(FullarrayGenesVIPSelVars)


NPLSDAarraygenesVIPselected<-NPLSDAmod(XN=FullarrayGenesVIPSelVars, YN=Outcomedummyarray, outcome.Y=NULL, factors=2, centering=0) 

# Plotting
ploteoNPLSDAVIPselectedgenes<- plotNPLSDAmod (X=NPLSDAarraygenesVIPselected, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                                              cutoff = 20, factors=2, penalty=2) 



