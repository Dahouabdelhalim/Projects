#--------------------------------------------------------------------------------------
#
#  July 2022
#  Code associated with publication:
#  "Kinetic modulation of bacterial hydrolases by microbial community structure in coastal waters" 
#   Authors: N. Abad, A. Uranga, B. Ayo, J.M. Arrieta,Z. Baña, I. Artolozaga, I. Azúa, J. Iriberri, Santos J. González-Rojí and M. Unanue
#
#
#  Creative Commons Licence: Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)
#--------------------------------------------------------------------------------------

### Adapted from Ruan et al (2006): doi:10.1093/bioinformatics/btl417

### The modified module from the original source is displayed below. The explanations of the original code 
### have been retained for the sake of better comprehension. Specific modifications related to a function are 
### indicated by the comment "#New: modified function".


### References of the modified code 

# Storey JD, Bass AJ, Dabney A, Robinson D (2022). 
#     qvalue: Q-value estimation for false discovery rate control. 
#     R package version 2.28.0, http://github.com/jdstorey/qvalue.

# Zeileis, A. & Grothendieck, G. (2005). 
#     zoo: S3 Infrastructure for Regular and Irregular Time Series. 
#     Journal of Statistical Software, 14(6), 1-27. doi:10.18637/jss.v014.i06


#--------------------------------------------------------------------------------------

######################################################################

# LSAscript.R:

# This script Computes LSA scores and PCC scores for the 
# ARISA time series as well as the env factors

############## 0. setup   ############################################

# 0_RAW DATA FRAME: 

# "m x n" (rowsxcols) data matrix, indicate missing values (NA)

## NOTE: set working directory in the same folder where all files provided 
## by Ruan et al 2006 package are stored.


setwd("C:\\\\Users\\\\...\\\\code\\\\") #New: change me

source("LSAfunctions.R")

datos<-read.table(file.choose(),header=TRUE, sep="\\t", na.strings= "",row.names = 1)
dim (datos)
str (datos)
head(datos)

#New: Linear interpolation of missing values (NAs) (Zeleis et al 2021)

library(zoo) 

idx <- colSums(!is.na(datos)) > 1
idx
datos[ , idx] <- na.approx(datos[ , idx])
dataTe<-datos
dataTe
str(dataTe)


############## 1. Preparing data  ##############################################

# normalization: this is not necessary for analysis but rather in the result plotting 
# since what matters in the next step is their ranks

dataXe <- normalization(dataTe)
dataXe <- as.data.frame(dataXe)
dataXe.nt <- normalTransform(dataXe)


############# 2. Computing Local Similarity scores ####################

NumPermu <- 1000

LocalSimilarity3Nai<- function(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F)  #New: modified function
  
{
  
  if(!is.array(lsTS1))
  {
    ls.TS1x <- t(lsTS1)
    ls.TS2x <- t(lsTS2)
  }
  
  ls.TS1 <- ls.TS1x
  ls.TS2 <- ls.TS2x
  
  
  if(scale==T)
  {
    lsTS1 <- (lsTS1-mean(lsTS1))/sd(lsTS1)  # Splus: stdev(lsTS1)
    lsTS2 <- (lsTS2-mean(lsTS2))/sd(lsTS2)  # Splus: stdev(lsTS2)	
  }
  
  scoreMatrixPos <- matrix(0, numTimePoints+1,numTimePoints+1);
  scoreMatrixNeg <- matrix(0, numTimePoints+1,numTimePoints+1);
  
  scoreMax <- 0.;
  PosOrNeg <- 0.;
  startX <- 0; # start of sub seq in dt1 from back to front
  startY <- 0; # start of sub seq in dt2 from back to front
  
  Thresh2 <- 0.000001
  Tabla=rep(NA, 5)
  TablaFinal=rep(NA, 5)
  
  for(i in 2:(numTimePoints+1)){
    for(j in 2:(numTimePoints+1)){
      if(abs(i-j) > maxDelay)
        next;
      
      scoreMatrixPos[i,j] <- scoreMatrixPos[i-1,j-1] + lsTS1[i-1] * lsTS2[j-1];
      if(scoreMatrixPos[i,j] < 0)
        scoreMatrixPos[i,j] <- 0;
      
      scoreMatrixNeg[i,j] <- scoreMatrixNeg[i-1,j-1] - lsTS1[i-1] * lsTS2[j-1];
      if(scoreMatrixNeg[i,j] < 0)
        scoreMatrixNeg[i,j] <- 0;
      
      if(scoreMatrixPos[i,j] > 0)
      {
        scoreMax <- scoreMatrixPos[i,j];
        startX <- i;
        startY <- j;
        PosOrNeg <- 1;	
      }	
      
      if(scoreMatrixNeg[i,j] > 0)
      {
        scoreMax <- scoreMatrixNeg[i,j];
        startX <- i;
        startY <- j;
        PosOrNeg <- 0;	
      }
      
      if(PosOrNeg == 1) {
        for(k in 1:numTimePoints)
          if(scoreMatrixPos[startX-k,startY-k]<=Thresh2){    		break;
          }
      } else {
        for(k in 1:numTimePoints)
          if(scoreMatrixNeg[startX-k, startY-k]<=Thresh2)
          {     		break;
          }
        
      }
      length = k;
      
   
      Exp=c(scoreMax, (startX-length), (startY-length), length, PosOrNeg)
      Tabla=rbind(Tabla,Exp)
      
    }
    
    if(!is.null(dim(Tabla))){
      Tabla=Tabla[-1,]
    }
    
    TablaFinal=rbind(TablaFinal,Tabla)
    Tabla=rep(NA, 5)
  }
  
  TablaFinal=TablaFinal[-1,]
  
  
  return(TablaFinal);
}

LocalSimilarity3Nailong<- function(lsTS1, lsTS2, maxDelay=1, numTimePoints=1, scale=F)  #New: modified function

  {
  
  if(!is.array(lsTS1))
  {
    ls.TS1x <- t(lsTS1)
    ls.TS2x <- t(lsTS2)
  }
  
  ls.TS1 <- ls.TS1x
  ls.TS2 <- ls.TS2x
  
  
  if(scale==T)
  {
    lsTS1 <- (lsTS1-mean(lsTS1))/sd(lsTS1)  # Splus: stdev(lsTS1)
    lsTS2 <- (lsTS2-mean(lsTS2))/sd(lsTS2)  # Splus: stdev(lsTS2)	
  }
  
  
  scoreMatrixPos <- matrix(0, numTimePoints+1,numTimePoints+1);
  scoreMatrixNeg <- matrix(0, numTimePoints+1,numTimePoints+1);
  
  scoreMax <- 0.;
  PosOrNeg <- 0.;
  startX <- 0; # start of sub seq in dt1 from back to front
  startY <- 0; # start of sub seq in dt2 from back to front
  
  Thresh2 <- 0.000001
  Tabla=rep(NA, 7)
  TablaFinal=rep(NA, 7)
  
  for(i in 2:(numTimePoints+1)){
    for(j in 2:(numTimePoints+1)){
      if(abs(i-j) > maxDelay)
        next;
      
      scoreMatrixPos[i,j] <- scoreMatrixPos[i-1,j-1] + lsTS1[i-1] * lsTS2[j-1];
      if(scoreMatrixPos[i,j] < 0)
        scoreMatrixPos[i,j] <- 0;
      
      scoreMatrixNeg[i,j] <- scoreMatrixNeg[i-1,j-1] - lsTS1[i-1] * lsTS2[j-1];
      if(scoreMatrixNeg[i,j] < 0)
        scoreMatrixNeg[i,j] <- 0;
      
      if(scoreMatrixPos[i,j] > 0)
      {
        scoreMax <- scoreMatrixPos[i,j];
        startX <- i;
        startY <- j;
        PosOrNeg <- 1;	
      }	
      
      if(scoreMatrixNeg[i,j] > 0)
      {
        scoreMax <- scoreMatrixNeg[i,j];
        startX <- i;
        startY <- j;
        PosOrNeg <- 0;	
      }
      
      if(PosOrNeg == 1) {
        for(k in 1:numTimePoints)
          if(scoreMatrixPos[startX-k,startY-k]<=Thresh2){    		break;
          }
      } else {
        for(k in 1:numTimePoints)
          if(scoreMatrixNeg[startX-k, startY-k]<=Thresh2)
          {     		break;
          }
        
      }
      length = k;
      

      if (PosOrNeg==1) {
        TS1 <-lsTS1[(startX-length):( (startX-length)+(length-1) )] 
        TS2 <-lsTS2[(startY-length):( (startY-length)+(length-1) )] 
      } else {
        TS1 <-lsTS1[(startX-length):( (startX-length)+(length-1) )] 
        TS2 <-rev(lsTS2[(startY-length):( (startY-length)+(length-1) )] )
      }
      corTmp <- cor(TS1, TS2)
      
      #calculate the pvalue
      if( is.na(corTmp) ){
        CorPval <- NA
      } else {
        CorPval <- 0.5 + sign(corTmp) * (0.5 - pt(corTmp * sqrt((dim1[[1]]-1)/(1-corTmp^2)), df=(dim1[[1]]-1)) )
      }
      Exp=c(scoreMax, (startX-length), (startY-length), length, PosOrNeg, corTmp, CorPval)
      Tabla=rbind(Tabla,Exp)
      
    }
    
    if(!is.null(dim(Tabla))){
      Tabla=Tabla[-1,]
    }
    
    TablaFinal=rbind(TablaFinal,Tabla)
    Tabla=rep(NA, 7)
  }
  
  TablaFinal=TablaFinal[-1,]
  

  
  return(TablaFinal);
}


data1=dataXe.nt
N=dim(data1)[[2]]
delay=3   #New: change me according to the established criteria
permu=NumPermu


dim1 <- dim(data1)

if(is.null(dim1))
  return (0)

rIdx <- 0;
SuperTabla <- rep(NA,10)
for(i in 1:(N-1)){
  for(j in (i+1):N){
    vectorPvalue <- 0
    lsTmp <- LocalSimilarity3Nailong(lsTS1 = data1[,i], lsTS2 = data1[,j], maxDelay=delay, numTimePoints=length(data1[,1]), scale=F); #New: modified function
    

    nrepeat=length(lsTmp[,1])	  # 	first seq index
    column1=rep(i,nrepeat)      #  	first seq index
    column2=rep(j,nrepeat)      #  	second seq index
    
    for( ncount in 1:nrepeat){
      pvalue <- sigTesting3(stOtu1 = data1[,i], stOtu2 = data1[,j], numPermu=permu, maxDelay=delay, numTimePoints=length(data1[,1]), scale=F,indexNai = ncount) #New: modified function
      vectorPvalue <- cbind(vectorPvalue, pvalue)
    }
    column10=vectorPvalue[-1]
    
    Exp=cbind(column1,column2,lsTmp,column10)
    
    SuperTabla=rbind(SuperTabla,Exp)
  }
}

SuperTabla=SuperTabla[-1,]


# 'normalize' the LS scores
SuperTabla[,3] <- SuperTabla[,3]/dim1[[1]]

dimnames(SuperTabla)[[2]] <- c("rIdx", "cIdx", "LSscore", "startR", "startC", "length", "PorN",  "cor", "corpVal","pValue")

# remove the empty rows
SuperTabla <- SuperTabla[SuperTabla[,1]>=1,]

fileName <- paste("SuperTabla.AllDelays", delay, 
                  ".ran", sample(10000:999999,1), ".txt", sep="") 

write.table(SuperTabla, paste(fileName, sep=""))

dim(SuperTabla)


# reading LS scores
ls.res.all <- SuperTabla
  
dim1 <- dim(ls.res.all) # n rows x 10
ls.res.all[1:3,]

dimnames2 <- dimnames(ls.res.all)[[2]]

tmp <- ls.res.all[,c(1:10, 10)]  # extend the result data frame to hold q-values
tmp[,11] <- 0

dimnames(tmp)[[2]] <- c(dimnames2, "q-val")


##### 3.Computing q-values by QVALUE (Storey et al 2022, http://github.com/jdstorey/qvalue) ########

library("devtools")

library(qvalue)  # need to install it in R first

p <- ls.res.all[,8]
qobj <- qvalue(p)

plot(qobj)
hist(qobj)
rownames(ls.res.all)<-NULL
ls.res.all <-as.data.frame(ls.res.all)
ls.res.all[,11] <- qobj$qvalues
head(ls.res.all)

ls.res.all2<-round(ls.res.all,4) 
ls.res.all2

head(ls.res.all2)

# Note: This need to be run only once 

write.table(ls.res.all2, file=paste("Set Path", "ls.res.all.txt", sep="")) # change me

####data ready:
# ls.res.all: an N by 11 matrix, with LSA and PCC score and their significance
# ls.dataTe.txt: data without normalization
# ls.dataXe.txt: normalized data
# ls.dataXe.nt.txt: data with Normal score transformation


############# 4. Analyzing results ###########################################

ls.res.all <- ls.res.all2

tmp <- ls.res.all
tmp
D<-(tmp[,4])-(tmp[,5])
D
tmp[,12] <- D
head(tmp)

## Filtering of the correlations

tmp <- tmp[tmp[,1]<=58,]
tmp <- tmp[tmp[,2]>58,]
dim(tmp) # n rows x 12

tmp1 <- tmp[tmp[,8]<=0.05,] # Filtering correlations with LS score p-value <0.05
dim(tmp1) 

tmp2 <- tmp1[tmp1[,11]<=0.05,] # Filtering correlations with Q-value <0.05
dim(tmp2)
 
tmp3 <- tmp2[tmp2[,12]<1,] # Filtering correlations with Delay <1
dim(tmp3)

tmp4 <- tmp3[tmp3[,12]>-1,] # Filtering correlations with Delay >-1
dim(tmp4)

tmp5 <- tmp4[tmp4[,6]>16,] # Filtering correlations with length >16
dim(tmp5)
head(tmp5)

tmp6 <- tmp5[tmp5[,3]>0.10,] # Filtering correlations with LS score value > 0.10
dim(tmp6)
head(tmp6)

# Matrix ordination and naming rIdx, cIdx, lenght y Delay

tmp6ord<-tmp6[order(tmp6[,1], tmp6[,2],tmp6[,6],tmp6[,12],decreasing=TRUE),]
head(tmp6ord)

write.csv(tmp6ord, file=paste("Set Path", "Set Name.csv", sep=","))


#Variable names

for (i in 1:58){
  tmp6ord$rIdx[tmp6ord$rIdx==i]<- colnames(dataTe)[i]
  tmp6ord$cIdx[tmp6ord$cIdx==i]<- colnames(dataTe)[i]
}

head(tmp6ord)

write.csv(tmp6ord, file=paste("Set Path", "Set Name.csv", sep=",")) #change me

