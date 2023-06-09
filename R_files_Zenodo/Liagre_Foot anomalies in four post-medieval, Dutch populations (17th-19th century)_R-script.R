#############################################################
## Elle Liagre (she/her)                                   ##
## Contact: Laboratory for Human Osteoarchaeology,         ##
##          Faculty of Archaeology, Leiden University,     ##
##          Einsteinweg 2, 2333 CC Leiden, The Netherlands ##
##          elleliagre@gmail.com                           ##
## Affiliation: Leiden University. Faculty of Archaeology, ##
##              Leiden, the Netherlands                    ##
##                                                         ##
## Date: 15/06/2021                                        ##
#############################################################

#############################################################################
##                           Description of code                           ##
## Data processing and statistical analyses of the data published in the   ##
## article "Kinship analysis using foot anomalies in a rural Dutch         ##
## cemetery (17th-19th c.)".                                               ##
#############################################################################


## Please clear environment `rm(list = ls())` and run the whole script ##


## Start of Analysis ##


# Dependencies ------------------------------------------------------------

library(tidyverse)
library(tibble)
library(janitor)
library(readr)

# Functions ---------------------------------------------------------------

###---###---###---###---###---###---###---###---###---###---###---
#
# Spatial Analysis Routines for Archaeology in R
# SARA-R
#
# V 0.1 Feb 17, 2014 (best considered an alpha version)
#
# Created by Jim Keron, Western University jkeron@uwo.ca
#
# This script is a series of R functions which must be executed before using 
# any of the routines.They will eventually be built into an R Library but not 
# yet there.
#
###---###---###---###---###---###---###---###---###---###---###---


###---###---###---###---###---###---###---###---###---###---###---
#
# This function calculates the distance between two points in Cartesion Space
# given by (x,y) (ie two dimensional space
#
###---###---###---###---###---###---###---###---###---###---###---
abDistance <- function(x1,y1,x2,y2)
{
  dist <- sqrt((x1-x2)^2+(y1-y2)^2)
  return(dist)
}


###---###---###---###---###---###---###---###---###---###---###---
#
# This function selects a random set of points of size n from an input data 
# frameisk list without replacement
# The result is a data frame of points (x,y) of size n
# All remaining points not selected are put into a second dataframe
#
###---###---###---###---###---###---###---###---###---###---###---
randABSplit <- function(indata,n)
{
  ###############Select sample for Type A
  #create index list to randomly select as sample the same size as the actuals file
  
  indexList <- 1:length(indata$x)
  indexSample <- sample(indexList, size=n, replace=FALSE)
  # set up x and y list of size from input n
  x1 <- numeric(n)
  y1 <- numeric(n)
  for (i in 1:n) {
    x1[i] <- indata$x[indexSample[i]]
    y1[i] <- indata$y[indexSample[i]]
  }
  typeASamp <- data.frame(x1,y1)
  #Assemble the left overs###
  m <- length(indata$x)- n
  x2 <- numeric(m)
  y2 <- numeric(m)
  o <- m+n
  p <- 1
  for(i in 1:o) {
    alreadySelected <- 0
    for (j in 1:n) {
      if (indexSample[j]==i) {
        alreadySelected <- 1
      }
    }
    if(alreadySelected==0) {
      x2[p] <- infile$x[i]
      y2[p] <- infile$y[i]
      p <- p+1
    }
  }
  typeBSamp <- data.frame(x2,y2)
  # done
  return(list("a"=typeASamp,"b"=typeBSamp))
}


###---###---###---###---###---###---###---###---###---###---###---
#
# This function calculates the p value from a set of statistics from a
# in a vector that were calculated from a set of randomizations
#
###---###---###---###---###---###---###---###---###---###---###---
pValueLessThan <- function(resultsList,actualStat)
{
  randomizations <- length(resultsList)
  countR <-0
  for (i in 1:randomizations)
  {
    if(resultsList[i]< actualStat)
    {
      countR <- countR+1
    }
  }
  probabilty <- countR/randomizations
  return(probabilty)
}
# endfunction


###---###---###---###---###---###---###---###---###---###---###---
#
# This counts the number of pairs of points withing the distance specified
# 
# Input listArray is the number of points to be checked
# mDist is the within-which distance
#
###---###---###---###---###---###---###---###---###---###---###---
proxCount <- function(listArray,mDist)
{
  nPoints <- length(listArray$x)
  count <- 0
  k <- nPoints-1
  for (i in 1:k) {
    l <- i+1
    for (j in l:nPoints) {
      aDist <-
        abDistance(listArray$x[i],listArray$y[i],listArray$x[j],listArray$y[j])
      if (aDist < mDist) {
        count <- count +1
      }
    }
    
  }
  return(count)
} #


###---###---###---###---###---###---###---###---###---###---###---
#
# This function calculates Hodder and Okells A statistc
#
###---###---###---###---###---###---###---###---###---###---###---
hoddersA <- function(a,b)
{
  # Calculate RAA####################
  n <- length(a$x)
  sumd <- 0
  k <- n-1
  for (i in 1:k) {
    l <- i+1
    for (j in l:n) {
      d1 <- abDistance(a$x[i],a$y[i],a$x[j],a$y[j])
      sumd <- sumd + d1
    }
  }
  raa <- sumd/((n^2-n)/2)
  # Calculate RBB####################
  n <- length(b$x)
  sumd <- 0
  k <- n-1
  for (i in 1:k) {
    l <- i+1
    for (j in l:n) {
      d1 <- abDistance(b$x[i],b$y[i],b$x[j],b$y[j])
      sumd <- sumd + d1
    }
  }
  rbb <- sumd/((n^2-n)/2)
  # Calculate RAB####################
  n <- length(a$x)
  m <- length(b$x)
  sumd <- 0
  for (i in 1:n) {
    for (j in 1:m) {
      d1 <- abDistance(a$x[i],a$y[i],b$x[j],b$y[j])
      sumd <- sumd + d1
    }
  }
  rab <- sumd/(n*m)
  # Calculate A-Staistic
  a <- (raa*rbb)/rab^2
  return(a)
}
#end function

###---###---###---###---###---###---###---###---###---###---###---
#
# This function calculates the p value from a set of statistics from a
# in a vector that were calculated from a set of randomizations
#
###---###---###---###---###---###---###---###---###---###---###---
pValueGreaterThan <- function(resultsList,actualStat)
{
  randomizations <- length(resultsList)
  countR <-0
  cv <- unlist(actualStat)
  for (i in 1:randomizations)
  {
    if(resultsList[i]>=cv)
    {
      countR <- countR+1
    }
  }
  probabilty <- countR/randomizations
  return(probabilty)
}
# endfunction


###---###---###---###---###---###---###---###---###---###---###---
#
# This function selects a random set of points of size n from an At-Risk list
# without replacement
# The result is a data frame of points (x,y) of size n
#
###---###---###---###---###---###---###---###---###---###---###---
randSampleXY <- function(atRiskList,n)
{
  x <- numeric(n)
  y <- numeric(n)
  indexList <- 1:length(atRiskList$x)

  #create index list to randomly select as sample the same size as the actuals file
  indexSample <- sample(indexList, size=n, replace=FALSE)
  for (i in 1:n) {
    x[i] <- atRiskList$x[indexSample[i]]
    y[i] <- atRiskList$y[indexSample[i]]
  }
  samplePoints <- data.frame(x,y)
  return(samplePoints)
}


###---###---###---###---###---###---###---###---###---###---###---
#
#This counts the number of events one type within a fixed distance of a 
# second (fixed) type
#Input The fixed points The varied points
# mDist is the within-which distance
#
###---###---###---###---###---###---###---###---###---###---###---
crossProxCount <- function(fixedSex,variedSex,mDist)
{
  count <- 0
  fixedCount <- length(fixedSex$x)
  variedCount <- length(variedSex$x)
  for (i in 1:fixedCount) {
    for (j in 1:variedCount) {
      aDist <-
        abDistance(fixedSex$x[i],fixedSex$y[i],variedSex$x[j],variedSex$y[j])
      if (aDist < mDist) {
        count <- count +1
      }
    }
  }
  return(count)
}


###---###---###---###---###---###---###---###---###---###---###---
#
# Nearest Neighbour calculations
#
###---###---###---###---###---###---###---###---###---###---###---
#
# Two sets of data are passed, from and to.
# fromData is the set of point from which the nearest neighbour is calculated.
# toData is the set of points which can form the nearest neighbour
# They can be the same et of points for a traditional NN analysis
# Or they can be different for a Cross NN analysis
#
# The value returned is the average nearest neighbour between from and to.
###---###---###---###---###---###---###---###---###---###---###---
avgNNDist <- function(fromData,toData)
{
  totDist <- 0
  nFrom <- length(fromData$x)
  nTo <- length(toData$x)
  avg <- 0
  for (i in 1:nFrom) {
    nnDist <- 999999999
    for (j in 1:nTo) {
      aDist <-
        abDistance(fromData$x[i],fromData$y[i],toData$x[j],toData$y[j])
      if (aDist < nnDist & aDist>0) {
        nnDist <- aDist
      }
    }
    totDist <- totDist + nnDist
  }
  avg <- totDist/nFrom
  return(avg)
} #
#### End Function
  


# Load data ---------------------------------------------------------------


# Please load the .csv file named 'Appendix S1'
Data <<- read.csv2(file.choose(), header = T, stringsAsFactors = T, 
                   fileEncoding = "UTF-8-BOM")

class(Data)
head(Data)
names(Data)

## Source of Data ##
      # This data was collected by Elle Liagre and available as supporting
      # information (Table S1) with the article "Kinship analysis using foot 
      # anomalies in a rural Dutch cemetery (17th-19th c.)". For more 
      # information on the dataset, see the 'Materials and methods' section 
      # of the article.


# Preparing data: feet2individual ------------------------------------------

#                                                                            #
##   The data collection happened for every foot. However, it is needed     ##
###  for every individual instead. Thus, for every trait both feet need    ###
###  to be combined. For more information, see the 'Materials and methods' ###
##   section of the article.                                                ##
#                                                                            #


Data_Ind <- Data    # new data frame for data per individuals

# AccessNav #
Data_Ind$AccessNav_L <- as.numeric(Data_Ind$AccessNav_L)
Data_Ind$AccessNav_R <- as.numeric(Data_Ind$AccessNav_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$AccessNav <- Data_Ind$AccessNav_L + Data_Ind$AccessNav_R
Data_Ind$AccessNav[Data_Ind$AccessNav == 2] <- 1
Data_Ind$AccessNav[Data_Ind$AccessNav == 7] <- 1
Data_Ind$AccessNav[Data_Ind$AccessNav == 10] <- 1
Data_Ind$AccessNav[Data_Ind$AccessNav == 12] <- 6
Data_Ind$AccessNav[Data_Ind$AccessNav == 15] <- 6
Data_Ind$AccessNav[Data_Ind$AccessNav == 18] <- 9
  # removing columns for left and right
Data_Ind$AccessNav_L <- NULL
Data_Ind$AccessNav_R <- NULL

# BrachyD #
Data_Ind$BrachyD_L <- as.numeric(Data_Ind$BrachyD_L)
Data_Ind$BrachyD_R <- as.numeric(Data_Ind$BrachyD_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$BrachyD <- Data_Ind$BrachyD_L + Data_Ind$BrachyD_R
Data_Ind$BrachyD[Data_Ind$BrachyD == 2] <- 1
Data_Ind$BrachyD[Data_Ind$BrachyD == 7] <- 1
Data_Ind$BrachyD[Data_Ind$BrachyD == 10] <- 1
Data_Ind$BrachyD[Data_Ind$BrachyD == 12] <- 6
Data_Ind$BrachyD[Data_Ind$BrachyD == 15] <- 6
Data_Ind$BrachyD[Data_Ind$BrachyD == 18] <- 9
  # removing columns for left and right
Data_Ind$BrachyD_L <- NULL
Data_Ind$BrachyD_R <- NULL

# BrachyMT1 #
Data_Ind$BrachyMT1_L <- as.numeric(Data_Ind$BrachyMT1_L)
Data_Ind$BrachyMT1_R <- as.numeric(Data_Ind$BrachyMT1_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$BrachyMT1 <- Data_Ind$BrachyMT1_L + Data_Ind$BrachyMT1_R
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 2] <- 1
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 7] <- 1
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 10] <- 1
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 12] <- 6
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 15] <- 6
Data_Ind$BrachyMT1[Data_Ind$BrachyMT1 == 18] <- 9
  # removing columns for left and right
Data_Ind$BrachyMT1_L <- NULL
Data_Ind$BrachyMT1_R <- NULL

# BrachyMT4 #
Data_Ind$BrachyMT4_L <- as.numeric(Data_Ind$BrachyMT4_L)
Data_Ind$BrachyMT4_R <- as.numeric(Data_Ind$BrachyMT4_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$BrachyMT4 <- Data_Ind$BrachyMT4_L + Data_Ind$BrachyMT4_R
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 2] <- 1
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 7] <- 1
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 10] <- 1
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 12] <- 6
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 15] <- 6
Data_Ind$BrachyMT4[Data_Ind$BrachyMT4 == 18] <- 9
  # removing columns for left and right
Data_Ind$BrachyMT4_L <- NULL
Data_Ind$BrachyMT4_R <- NULL

# BrachyPP1 #
Data_Ind$BrachyPP1_L <- as.numeric(Data_Ind$BrachyPP1_L)
Data_Ind$BrachyPP1_R <- as.numeric(Data_Ind$BrachyPP1_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$BrachyPP1 <- Data_Ind$BrachyPP1_L + Data_Ind$BrachyPP1_R
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 2] <- 1
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 7] <- 1
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 10] <- 1
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 12] <- 6
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 15] <- 6
Data_Ind$BrachyPP1[Data_Ind$BrachyPP1 == 18] <- 9
  # removing columns for left and right
Data_Ind$BrachyPP1_L <- NULL
Data_Ind$BrachyPP1_R <- NULL

# CalcCubCoal #
Data_Ind$CalcCubCoal_L <- as.numeric(Data_Ind$CalcCubCoal_L)
Data_Ind$CalcCubCoal_R <- as.numeric(Data_Ind$CalcCubCoal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$CalcCubCoal <- Data_Ind$CalcCubCoal_L + Data_Ind$CalcCubCoal_R
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 2] <- 1
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 7] <- 1
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 10] <- 1
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 12] <- 6
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 15] <- 6
Data_Ind$CalcCubCoal[Data_Ind$CalcCubCoal == 18] <- 9
  # removing columns for left and right
Data_Ind$CalcCubCoal_L <- NULL
Data_Ind$CalcCubCoal_R <- NULL

# CalcNavCoal #
Data_Ind$CalcNavCoal_L <- as.numeric(Data_Ind$CalcNavCoal_L)
Data_Ind$CalcNavCoal_R <- as.numeric(Data_Ind$CalcNavCoal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$CalcNavCoal <- Data_Ind$CalcNavCoal_L + Data_Ind$CalcNavCoal_R
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 2] <- 1
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 7] <- 1
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 10] <- 1
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 12] <- 6
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 15] <- 6
Data_Ind$CalcNavCoal[Data_Ind$CalcNavCoal == 18] <- 9
  # removing columns for left and right
Data_Ind$CalcNavCoal_L <- NULL
Data_Ind$CalcNavCoal_R <- NULL

# CF2CF3Coal #
Data_Ind$CF2CF3Coal_L <- as.numeric(Data_Ind$CF2CF3Coal_L)
Data_Ind$CF2CF3Coal_R <- as.numeric(Data_Ind$CF2CF3Coal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$CF2CF3Coal <- Data_Ind$CF2CF3Coal_L + Data_Ind$CF2CF3Coal_R
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 2] <- 1
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 7] <- 1
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 10] <- 1
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 12] <- 6
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 15] <- 6
Data_Ind$CF2CF3Coal[Data_Ind$CF2CF3Coal == 18] <- 9
  # removing columns for left and right
Data_Ind$CF2CF3Coal_L <- NULL
Data_Ind$CF2CF3Coal_R <- NULL

# CF3MT3Coal #
Data_Ind$CF3MT3Coal_L <- as.numeric(Data_Ind$CF3MT3Coal_L)
Data_Ind$CF3MT3Coal_R <- as.numeric(Data_Ind$CF3MT3Coal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$CF3MT3Coal <- Data_Ind$CF3MT3Coal_L + Data_Ind$CF3MT3Coal_R
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 2] <- 1
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 7] <- 1
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 10] <- 1
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 12] <- 6
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 15] <- 6
Data_Ind$CF3MT3Coal[Data_Ind$CF3MT3Coal == 18] <- 9
  # removing columns for left and right
Data_Ind$CF3MT3Coal_L <- NULL
Data_Ind$CF3MT3Coal_R <- NULL

# CF1Intermet #
Data_Ind$CF1Intermet_L <- as.numeric(Data_Ind$CF1Intermet_L)
Data_Ind$CF1Intermet_R <- as.numeric(Data_Ind$CF1Intermet_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$CF1Intermet <- Data_Ind$CF1Intermet_L + Data_Ind$CF1Intermet_R
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 2] <- 1
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 7] <- 1
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 10] <- 1
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 12] <- 6
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 15] <- 6
Data_Ind$CF1Intermet[Data_Ind$CF1Intermet == 18] <- 9
  # removing columns for left and right
Data_Ind$CF1Intermet_L <- NULL
Data_Ind$CF1Intermet_R <- NULL

# MT1Intermet #
Data_Ind$MT1Intermet_L <- as.numeric(Data_Ind$MT1Intermet_L)
Data_Ind$MT1Intermet_R <- as.numeric(Data_Ind$MT1Intermet_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$MT1Intermet <- Data_Ind$MT1Intermet_L + Data_Ind$MT1Intermet_R
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 2] <- 1
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 7] <- 1
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 10] <- 1
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 12] <- 6
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 15] <- 6
Data_Ind$MT1Intermet[Data_Ind$MT1Intermet == 18] <- 9
  # removing columns for left and right
Data_Ind$MT1Intermet_L <- NULL
Data_Ind$MT1Intermet_R <- NULL

# MT2Intermet #
Data_Ind$MT2Intermet_L <- as.numeric(Data_Ind$MT2Intermet_L)
Data_Ind$MT2Intermet_R <- as.numeric(Data_Ind$MT2Intermet_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$MT2Intermet <- Data_Ind$MT2Intermet_L + Data_Ind$MT2Intermet_R
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 2] <- 1
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 7] <- 1
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 10] <- 1
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 12] <- 6
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 15] <- 6
Data_Ind$MT2Intermet[Data_Ind$MT2Intermet == 18] <- 9
  # removing columns for left and right
Data_Ind$MT2Intermet_L <- NULL
Data_Ind$MT2Intermet_R <- NULL

# TaloCalcCoal #
Data_Ind$TaloCalcCoal_L <- as.numeric(Data_Ind$TaloCalcCoal_L)
Data_Ind$TaloCalcCoal_R <- as.numeric(Data_Ind$TaloCalcCoal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$TaloCalcCoal <- Data_Ind$TaloCalcCoal_L + Data_Ind$TaloCalcCoal_R
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 2] <- 1
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 7] <- 1
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 10] <- 1
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 12] <- 6
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 15] <- 6
Data_Ind$TaloCalcCoal[Data_Ind$TaloCalcCoal == 18] <- 9
  # removing columns for left and right
Data_Ind$TaloCalcCoal_L <- NULL
Data_Ind$TaloCalcCoal_R <- NULL

# TaloNavCoal #
Data_Ind$TaloNavCoal_L <- as.numeric(Data_Ind$TaloNavCoal_L)
Data_Ind$TaloNavCoal_R <- as.numeric(Data_Ind$TaloNavCoal_R)
  # making new column with sum of columns for left and right
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$TaloNavCoal <- Data_Ind$TaloNavCoal_L + Data_Ind$TaloNavCoal_R
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 2] <- 1
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 7] <- 1
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 10] <- 1
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 12] <- 6
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 15] <- 6
Data_Ind$TaloNavCoal[Data_Ind$TaloNavCoal == 18] <- 9
  # removing columns for left and right
Data_Ind$TaloNavCoal_L <- NULL
Data_Ind$TaloNavCoal_R <- NULL


# Preparing data: clustering ----------------------------------------------

#                                                                         #
##  The traits 'Os Intermetatarseum (IntermetCF1, IntermetMT1 and        ##
### IntermetMT2)' and 'Calcaneonavicular and Talocalcaneal coalition    ### 
### (CalcNavCoal and TaloCalcCoal)' need to be clustered.For more       ###
##  information, see the 'Materials and methods' section of the article. ##
#                                                                         #

## Os Intermetatarseum ##
# CF1Intermet and MT1Intermet #
Data_Ind$CF1Intermet <- as.numeric(Data_Ind$CF1Intermet)
Data_Ind$MT1Intermet <- as.numeric(Data_Ind$MT1Intermet)
  # making new column with sum of columns for 2 out of 3 traits
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$OsIntermet1 <- Data_Ind$CF1Intermet + Data_Ind$MT1Intermet
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 2] <- 1
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 7] <- 1
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 10] <- 1
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 12] <- 6
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 15] <- 6
Data_Ind$OsIntermet1[Data_Ind$OsIntermet1 == 18] <- 9


# CF1Intermet, MT1Intermet and MT2Intermet #
Data_Ind$OsIntermet1 <- as.numeric(Data_Ind$OsIntermet1)
Data_Ind$MT2Intermet <- as.numeric(Data_Ind$MT2Intermet)
  # making new column with sum of columns for 3rd trait with sum of 1st and 2nd 
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$OsIntermet <- Data_Ind$OsIntermet1 + Data_Ind$MT2Intermet
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 2] <- 1
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 7] <- 1
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 10] <- 1
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 12] <- 6
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 15] <- 6
Data_Ind$OsIntermet[Data_Ind$OsIntermet == 18] <- 9
  # removing columns for left and right
Data_Ind$OsIntermet1 <- NULL


## Ankle coalition ##
Data_Ind$TaloCalcCoal <- as.numeric(Data_Ind$TaloCalcCoal)
Data_Ind$CalcNavCoal <- as.numeric(Data_Ind$CalcNavCoal)
  # making new column with sum of columns for both traits,
  # but keeping same values as in data collection (0, 1, 6, 9)
Data_Ind$AnkleCoal <- Data_Ind$TaloCalcCoal + Data_Ind$CalcNavCoal
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 2] <- 1
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 7] <- 1
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 10] <- 1
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 12] <- 6
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 15] <- 6
Data_Ind$AnkleCoal[Data_Ind$AnkleCoal == 18] <- 9



# Subsets of data ---------------------------------------------------------

#                                                      #
## To compare the data, the dataframe shoudl be split ##
## up according to their respective sites and sexes.  ##
#                                                      #

## Preparing Middenbeemster dataset ##
  # the Middenbeemster dataset (MB) #
MB <<- Data_Ind[Data_Ind$Site == "MB",]         
MB_NumberInd <- length(MB$Key.ID)
paste("The number of studied individuals (#) in the Middenbeemster skeletal", 
      "collection is", MB_NumberInd, "( N =",MB_NumberInd,").", 
      "The selected portion of the collection is", 
      round(MB_NumberInd/452,2)*100, "%.")

  # MB according to sex #
      # removing individuals with 'indeterminate (I)' sex in a new data set
MB_wo_probables <- MB[MB$Sex != "I",]    
      # gathering probable and determined sexes in same category
MB_wo_probables$Sex[MB_wo_probables$Sex == "PF"] <- "F"   
MB_wo_probables$Sex[MB_wo_probables$Sex == "PM"] <- "M"
MB_F <<- MB_wo_probables[MB_wo_probables$Sex == "F",]    #Females in MB
MB_M <<- MB_wo_probables[MB_wo_probables$Sex == "M",]    #Males in MB


## Preparing reference dataset ##
  # the Arnhem dataset (AR) #
AR <<- Data_Ind[Data_Ind$Site == "AR",]         
AR_NumberInd <- length(AR$Key.ID)
paste("The number of studied individuals (#) in the Arnhem skeletal", 
      "collection is", AR_NumberInd, "( N =",AR_NumberInd,").", 
      "The selected portion of the collection is", 
      round(AR_NumberInd/750,3)*100, "%.")

  # AR according to sex #
      # removing individuals with 'indeterminate (I)' sex in a new data set
AR_wo_probables <- AR[AR$Sex != "I",]    
      # gathering probable and determined sexes in same category
AR_wo_probables$Sex[AR_wo_probables$Sex == "PF"] <- "F"   
AR_wo_probables$Sex[AR_wo_probables$Sex == "PM"] <- "M"
AR_F <<- AR_wo_probables[AR_wo_probables$Sex == "F",]    #Females in AR
AR_M <<- AR_wo_probables[AR_wo_probables$Sex == "M",]    #Males in AR


  # the Eindhoven dataset (EH) #
EH <<- Data_Ind[Data_Ind$Site == "EH",]         
EH_NumberInd <- length(EH$Key.ID)
paste("The number of studied individuals (#) in the Eindhoven skeletal", 
      "collection is", EH_NumberInd, "( N =",EH_NumberInd,").", 
      "The selected portion of the collection is", 
      round(EH_NumberInd/752,3)*100, "%.")

  # EH according to sex #
      # removing individuals with 'indeterminate (I)' sex in a new data set
EH_wo_probables <- EH[EH$Sex != "I",]    
      # gathering probable and determined sexes in same category
EH_wo_probables$Sex[EH_wo_probables$Sex == "PF"] <- "F"   
EH_wo_probables$Sex[EH_wo_probables$Sex == "PM"] <- "M"
EH_F <<- EH_wo_probables[EH_wo_probables$Sex == "F",]    #Females in EH
EH_M <<- EH_wo_probables[EH_wo_probables$Sex == "M",]    #Males in EH


  # the Zwolle dataset (ZW) #
ZW <<- Data_Ind[Data_Ind$Site == "ZW",]         
ZW_NumberInd <- length(ZW$Key.ID)
paste("The number of studied individuals (#) in the Zwolle skeletal", 
      "collection is", ZW_NumberInd, "( N =",ZW_NumberInd,").", 
      "The selected portion of the collection is", 
      round(ZW_NumberInd/529,3)*100, "%.")

  # ZW according to sex #
      # removing individuals with 'indeterminate (I)' sex in a new data set
ZW_wo_probables <- ZW[ZW$Sex != "I",]
      # gathering probable and determined sexes in same category
ZW_wo_probables$Sex[ZW_wo_probables$Sex == "PF"] <- "F"   
ZW_wo_probables$Sex[ZW_wo_probables$Sex == "PM"] <- "M"
ZW_F <<- ZW_wo_probables[ZW_wo_probables$Sex == "F",]    #Females in ZW
ZW_M <<- ZW_wo_probables[ZW_wo_probables$Sex == "M",]    #Males in ZW


  # the reference sample #
      # combining all skeletal collections in reference samples
Ref_Overall <<- rbind(AR_wo_probables, EH_wo_probables, ZW_wo_probables) 
      #creating reference sample according to sex
Ref_Overall_M <<- rbind(AR_M, EH_M, ZW_M) 
Ref_Overall_F <<- rbind(AR_F, EH_F, ZW_F)

# Frequency calculations: MB ----------------------------------------------

#                                                             #
## Calculating the frequency of the traits in the MB dataset ##
#                                                             #

## AccessNav ##
# Overall N and frequency #
MB_AccessNav <- MB_wo_probables %>% count(AccessNav)
MB_AccessNav <- 
  MB_AccessNav[(MB_AccessNav$AccessNav != 6 & MB_AccessNav$AccessNav != 9),]
MB_AccessNav_0 <- MB_AccessNav[1,2]
MB_AccessNav_0[is.na(MB_AccessNav_0)] <-0
MB_AccessNav_1 <- MB_AccessNav[2,2]
MB_AccessNav_1[is.na(MB_AccessNav_1)] <-0
MB_AccessNav_Total <- MB_AccessNav_0 + MB_AccessNav_1
MB_AccessNav_Affected <- round((MB_AccessNav_1/MB_AccessNav_Total)*100, 2)

# Male N and frequency
MB_M_AccessNav <- MB_M %>% count(AccessNav)
MB_M_AccessNav <- 
  MB_M_AccessNav[(MB_M_AccessNav$AccessNav != 6 & MB_M_AccessNav$AccessNav != 9),]
MB_M_AccessNav_0 <- MB_M_AccessNav[1,2]
MB_M_AccessNav_0[is.na(MB_M_AccessNav_0)] <-0
MB_M_AccessNav_1 <- MB_M_AccessNav[2,2]
MB_M_AccessNav_1[is.na(MB_M_AccessNav_1)] <-0
MB_M_AccessNav_Total <- MB_M_AccessNav_0 + MB_M_AccessNav_1
MB_M_AccessNav_Affected <- round((MB_M_AccessNav_1/MB_M_AccessNav_Total)*100, 2)

# Female N and frequency
MB_F_AccessNav <- MB_F %>% count(AccessNav)
MB_F_AccessNav <- 
  MB_F_AccessNav[(MB_F_AccessNav$AccessNav != 6 & MB_F_AccessNav$AccessNav != 9),]
MB_F_AccessNav_0 <- MB_F_AccessNav[1,2]
MB_F_AccessNav_0[is.na(MB_F_AccessNav_0)] <-0
MB_F_AccessNav_1 <- MB_F_AccessNav[2,2]
MB_F_AccessNav_1[is.na(MB_F_AccessNav_1)] <-0
MB_F_AccessNav_Total <- MB_F_AccessNav_0 + MB_F_AccessNav_1
MB_F_AccessNav_Affected <- round((MB_F_AccessNav_1/MB_F_AccessNav_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_AccessNav_row <- c("AccessNav", MB_AccessNav_Total, MB_AccessNav_1, 
                      MB_AccessNav_Affected, MB_M_AccessNav_Total, 
                      MB_M_AccessNav_1, MB_M_AccessNav_Affected, 
                      MB_F_AccessNav_Total, MB_F_AccessNav_1, 
                      MB_F_AccessNav_Affected)

## BrachyD ##
# Overall N and frequency #
MB_BrachyD <- MB_wo_probables %>% count(BrachyD)
MB_BrachyD <- MB_BrachyD[(MB_BrachyD$BrachyD != 6 & MB_BrachyD$BrachyD != 9),]
MB_BrachyD_0 <- MB_BrachyD[1,2]; MB_BrachyD_0[is.na(MB_BrachyD_0)] <-0
MB_BrachyD_1 <- MB_BrachyD[2,2]; MB_BrachyD_1[is.na(MB_BrachyD_1)] <-0
MB_BrachyD_Total <- na.omit(MB_BrachyD_0) + na.omit(MB_BrachyD_1)
MB_BrachyD_Affected <- round((MB_BrachyD_1/MB_BrachyD_Total)*100, 2)

# Male N and frequency
MB_M_BrachyD <- MB_M %>% count(BrachyD)
MB_M_BrachyD <- 
  MB_M_BrachyD[(MB_M_BrachyD$BrachyD != 6 & MB_M_BrachyD$BrachyD != 9),]
MB_M_BrachyD_0 <- MB_M_BrachyD[1,2]; MB_M_BrachyD_0[is.na(MB_M_BrachyD_0)] <-0
MB_M_BrachyD_1 <- MB_M_BrachyD[2,2]; MB_M_BrachyD_1[is.na(MB_M_BrachyD_1)] <-0
MB_M_BrachyD_Total <- MB_M_BrachyD_0 + MB_M_BrachyD_1
MB_M_BrachyD_Affected <- round((MB_M_BrachyD_1/MB_M_BrachyD_Total)*100, 2)

# Female N and frequency
MB_F_BrachyD <- MB_F %>% count(BrachyD)
MB_F_BrachyD <- 
  MB_F_BrachyD[(MB_F_BrachyD$BrachyD != 6 & MB_F_BrachyD$BrachyD != 9),]
MB_F_BrachyD_0 <- MB_F_BrachyD[1,2]; MB_F_BrachyD_0[is.na(MB_F_BrachyD_0)] <-0
MB_F_BrachyD_1 <- MB_F_BrachyD[2,2]; MB_F_BrachyD_1[is.na(MB_F_BrachyD_1)] <-0
MB_F_BrachyD_Total <- MB_F_BrachyD_0 + MB_F_BrachyD_1
MB_F_BrachyD_Affected <- round((MB_F_BrachyD_1/MB_F_BrachyD_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_BrachyD_row <- c("BrachyD", MB_BrachyD_Total, MB_BrachyD_1, 
                    MB_BrachyD_Affected, MB_M_BrachyD_Total, 
                    MB_M_BrachyD_1, MB_M_BrachyD_Affected, 
                    MB_F_BrachyD_Total, MB_F_BrachyD_1, 
                    MB_F_BrachyD_Affected)


## BrachyMT1 ##
# Overall N and frequency #
MB_BrachyMT1 <- MB_wo_probables %>% count(BrachyMT1)
MB_BrachyMT1 <- 
  MB_BrachyMT1[(MB_BrachyMT1$BrachyMT1 != 6 & MB_BrachyMT1$BrachyMT1 != 9),]
MB_BrachyMT1_0 <- MB_BrachyMT1[1,2]; MB_BrachyMT1_0[is.na(MB_BrachyMT1_0)] <-0
MB_BrachyMT1_1 <- MB_BrachyMT1[2,2]; MB_BrachyMT1_1[is.na(MB_BrachyMT1_1)] <-0
MB_BrachyMT1_Total <- na.omit(MB_BrachyMT1_0) + na.omit(MB_BrachyMT1_1)
MB_BrachyMT1_Affected <- round((MB_BrachyMT1_1/MB_BrachyMT1_Total)*100, 2)

# Male N and frequency
MB_M_BrachyMT1 <- MB_M %>% count(BrachyMT1)
MB_M_BrachyMT1 <- 
  MB_M_BrachyMT1[(MB_M_BrachyMT1$BrachyMT1 != 6 & MB_M_BrachyMT1$BrachyMT1 != 9),]
MB_M_BrachyMT1_0 <- MB_M_BrachyMT1[1,2]
MB_M_BrachyMT1_0[is.na(MB_M_BrachyMT1_0)] <-0
MB_M_BrachyMT1_1 <- MB_M_BrachyMT1[2,2]
MB_M_BrachyMT1_1[is.na(MB_M_BrachyMT1_1)] <-0
MB_M_BrachyMT1_Total <- MB_M_BrachyMT1_0 + MB_M_BrachyMT1_1
MB_M_BrachyMT1_Affected <- round((MB_M_BrachyMT1_1/MB_M_BrachyMT1_Total)*100, 2)

# Female N and frequency
MB_F_BrachyMT1 <- MB_F %>% count(BrachyMT1)
MB_F_BrachyMT1 <- 
  MB_F_BrachyMT1[(MB_F_BrachyMT1$BrachyMT1 != 6 & MB_F_BrachyMT1$BrachyMT1 != 9),]
MB_F_BrachyMT1_0 <- MB_F_BrachyMT1[1,2]
MB_F_BrachyMT1_0[is.na(MB_F_BrachyMT1_0)] <-0
MB_F_BrachyMT1_1 <- MB_F_BrachyMT1[2,2]
MB_F_BrachyMT1_1[is.na(MB_F_BrachyMT1_1)] <-0
MB_F_BrachyMT1_Total <- MB_F_BrachyMT1_0 + MB_F_BrachyMT1_1
MB_F_BrachyMT1_Affected <- round((MB_F_BrachyMT1_1/MB_F_BrachyMT1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_BrachyMT1_row <- c("BrachyMT1", MB_BrachyMT1_Total, MB_BrachyMT1_1, 
                      MB_BrachyMT1_Affected, MB_M_BrachyMT1_Total, 
                      MB_M_BrachyMT1_1, MB_M_BrachyMT1_Affected, 
                      MB_F_BrachyMT1_Total, MB_F_BrachyMT1_1, 
                      MB_F_BrachyMT1_Affected)


## BrachyMT4 ##
# Overall N and frequency #
MB_BrachyMT4 <- MB_wo_probables %>% count(BrachyMT4)
MB_BrachyMT4 <- 
  MB_BrachyMT4[(MB_BrachyMT4$BrachyMT4 != 6 & MB_BrachyMT4$BrachyMT4 != 9),]
MB_BrachyMT4_0 <- MB_BrachyMT4[1,2]
MB_BrachyMT4_0[is.na(MB_BrachyMT4_0)] <-0
MB_BrachyMT4_1 <- MB_BrachyMT4[2,2]
MB_BrachyMT4_1[is.na(MB_BrachyMT4_1)] <-0
MB_BrachyMT4_Total <- na.omit(MB_BrachyMT4_0) + na.omit(MB_BrachyMT4_1)
MB_BrachyMT4_Affected <- round((MB_BrachyMT4_1/MB_BrachyMT4_Total)*100, 2)

# Male N and frequency
MB_M_BrachyMT4 <- MB_M %>% count(BrachyMT4)
MB_M_BrachyMT4 <- 
  MB_M_BrachyMT4[(MB_M_BrachyMT4$BrachyMT4 != 6 & MB_M_BrachyMT4$BrachyMT4 != 9),]
MB_M_BrachyMT4_0 <- MB_M_BrachyMT4[1,2]
MB_M_BrachyMT4_0[is.na(MB_M_BrachyMT4_0)] <-0
MB_M_BrachyMT4_1 <- MB_M_BrachyMT4[2,2]
MB_M_BrachyMT4_1[is.na(MB_M_BrachyMT4_1)] <-0
MB_M_BrachyMT4_Total <- MB_M_BrachyMT4_0 + MB_M_BrachyMT4_1
MB_M_BrachyMT4_Affected <- round((MB_M_BrachyMT4_1/MB_M_BrachyMT4_Total)*100, 2)

# Female N and frequency
MB_F_BrachyMT4 <- MB_F %>% count(BrachyMT4)
MB_F_BrachyMT4 <- 
  MB_F_BrachyMT4[(MB_F_BrachyMT4$BrachyMT4 != 6 & MB_F_BrachyMT4$BrachyMT4 != 9),]
MB_F_BrachyMT4_0 <- MB_F_BrachyMT4[1,2]
MB_F_BrachyMT4_0[is.na(MB_F_BrachyMT4_0)] <-0
MB_F_BrachyMT4_1 <- MB_F_BrachyMT4[2,2]
MB_F_BrachyMT4_1[is.na(MB_F_BrachyMT4_1)] <-0
MB_F_BrachyMT4_Total <- MB_F_BrachyMT4_0 + MB_F_BrachyMT4_1
MB_F_BrachyMT4_Affected <- round((MB_F_BrachyMT4_1/MB_F_BrachyMT4_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_BrachyMT4_row <- c("BrachyMT4", MB_BrachyMT4_Total, MB_BrachyMT4_1, 
                      MB_BrachyMT4_Affected, MB_M_BrachyMT4_Total, 
                      MB_M_BrachyMT4_1, MB_M_BrachyMT4_Affected, 
                      MB_F_BrachyMT4_Total, MB_F_BrachyMT4_1, 
                      MB_F_BrachyMT4_Affected)

## BrachyPP1 ##
# Overall N and frequency #
MB_BrachyPP1 <- MB_wo_probables %>% count(BrachyPP1)
MB_BrachyPP1 <- 
  MB_BrachyPP1[(MB_BrachyPP1$BrachyPP1 != 6 & MB_BrachyPP1$BrachyPP1 != 9),]
MB_BrachyPP1_0 <- MB_BrachyPP1[1,2]; MB_BrachyPP1_0[is.na(MB_BrachyPP1_0)] <-0
MB_BrachyPP1_1 <- MB_BrachyPP1[2,2]; MB_BrachyPP1_1[is.na(MB_BrachyPP1_1)] <-0
MB_BrachyPP1_Total <- na.omit(MB_BrachyPP1_0) + na.omit(MB_BrachyPP1_1)
MB_BrachyPP1_Affected <- round((MB_BrachyPP1_1/MB_BrachyPP1_Total)*100, 2)

# Male N and frequency
MB_M_BrachyPP1 <- MB_M %>% count(BrachyPP1)
MB_M_BrachyPP1 <- 
  MB_M_BrachyPP1[(MB_M_BrachyPP1$BrachyPP1 != 6 & MB_M_BrachyPP1$BrachyPP1 != 9),]
MB_M_BrachyPP1_0 <- MB_M_BrachyPP1[1,2]
MB_M_BrachyPP1_0[is.na(MB_M_BrachyPP1_0)] <-0
MB_M_BrachyPP1_1 <- MB_M_BrachyPP1[2,2]
MB_M_BrachyPP1_1[is.na(MB_M_BrachyPP1_1)] <-0
MB_M_BrachyPP1_Total <- MB_M_BrachyPP1_0 + MB_M_BrachyPP1_1
MB_M_BrachyPP1_Affected <- round((MB_M_BrachyPP1_1/MB_M_BrachyPP1_Total)*100, 2)

# Female N and frequency
MB_F_BrachyPP1 <- MB_F %>% count(BrachyPP1)
MB_F_BrachyPP1 <- 
  MB_F_BrachyPP1[(MB_F_BrachyPP1$BrachyPP1 != 6 & MB_F_BrachyPP1$BrachyPP1 != 9),]
MB_F_BrachyPP1_0 <- MB_F_BrachyPP1[1,2]
MB_F_BrachyPP1_0[is.na(MB_F_BrachyPP1_0)] <-0
MB_F_BrachyPP1_1 <- MB_F_BrachyPP1[2,2]
MB_F_BrachyPP1_1[is.na(MB_F_BrachyPP1_1)] <-0
MB_F_BrachyPP1_Total <- MB_F_BrachyPP1_0 + MB_F_BrachyPP1_1
MB_F_BrachyPP1_Affected <- round((MB_F_BrachyPP1_1/MB_F_BrachyPP1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_BrachyPP1_row <- c("BrachyPP1", MB_BrachyPP1_Total, MB_BrachyPP1_1, 
                      MB_BrachyPP1_Affected, MB_M_BrachyPP1_Total, 
                      MB_M_BrachyPP1_1, MB_M_BrachyPP1_Affected, 
                      MB_F_BrachyPP1_Total, MB_F_BrachyPP1_1, 
                      MB_F_BrachyPP1_Affected)

## CalcCubCoal ##
# Overall N and frequency #
MB_CalcCubCoal <- MB_wo_probables %>% count(CalcCubCoal)
MB_CalcCubCoal <- 
  MB_CalcCubCoal[(MB_CalcCubCoal$CalcCubCoal != 6 & MB_CalcCubCoal$CalcCubCoal != 9),]
MB_CalcCubCoal_0 <- MB_CalcCubCoal[1,2]
MB_CalcCubCoal_0[is.na(MB_CalcCubCoal_0)] <-0
MB_CalcCubCoal_1 <- MB_CalcCubCoal[2,2]
MB_CalcCubCoal_1[is.na(MB_CalcCubCoal_1)] <-0
MB_CalcCubCoal_Total <- na.omit(MB_CalcCubCoal_0) + na.omit(MB_CalcCubCoal_1)
MB_CalcCubCoal_Affected <- round((MB_CalcCubCoal_1/MB_CalcCubCoal_Total)*100, 2)

# Male N and frequency
MB_M_CalcCubCoal <- MB_M %>% count(CalcCubCoal)
MB_M_CalcCubCoal <- 
  MB_M_CalcCubCoal[(MB_M_CalcCubCoal$CalcCubCoal != 6 & MB_M_CalcCubCoal$CalcCubCoal != 9),]
MB_M_CalcCubCoal_0 <- MB_M_CalcCubCoal[1,2]
MB_M_CalcCubCoal_0[is.na(MB_M_CalcCubCoal_0)] <-0
MB_M_CalcCubCoal_1 <- MB_M_CalcCubCoal[2,2]
MB_M_CalcCubCoal_1[is.na(MB_M_CalcCubCoal_1)] <-0
MB_M_CalcCubCoal_Total <- MB_M_CalcCubCoal_0 + MB_M_CalcCubCoal_1
MB_M_CalcCubCoal_Affected <- round((MB_M_CalcCubCoal_1/MB_M_CalcCubCoal_Total)*100, 2)

# Female N and frequency
MB_F_CalcCubCoal <- MB_F %>% count(CalcCubCoal)
MB_F_CalcCubCoal <- 
  MB_F_CalcCubCoal[(MB_F_CalcCubCoal$CalcCubCoal != 6 & MB_F_CalcCubCoal$CalcCubCoal != 9),]
MB_F_CalcCubCoal_0 <- MB_F_CalcCubCoal[1,2]
MB_F_CalcCubCoal_0[is.na(MB_F_CalcCubCoal_0)] <-0
MB_F_CalcCubCoal_1 <- MB_F_CalcCubCoal[2,2]
MB_F_CalcCubCoal_1[is.na(MB_F_CalcCubCoal_1)] <-0
MB_F_CalcCubCoal_Total <- MB_F_CalcCubCoal_0 + MB_F_CalcCubCoal_1
MB_F_CalcCubCoal_Affected <- 
  round((MB_F_CalcCubCoal_1/MB_F_CalcCubCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_CalcCubCoal_row <- c("CalcCubCoal", MB_CalcCubCoal_Total, MB_CalcCubCoal_1, 
                        MB_CalcCubCoal_Affected, MB_M_CalcCubCoal_Total, 
                        MB_M_CalcCubCoal_1, MB_M_CalcCubCoal_Affected, 
                        MB_F_CalcCubCoal_Total, MB_F_CalcCubCoal_1, 
                        MB_F_CalcCubCoal_Affected)

## TaloNavCoal ##
# Overall N and frequency #
MB_TaloNavCoal <- MB_wo_probables %>% count(TaloNavCoal)
MB_TaloNavCoal <- 
  MB_TaloNavCoal[(MB_TaloNavCoal$TaloNavCoal != 6 & MB_TaloNavCoal$TaloNavCoal != 9),]
MB_TaloNavCoal_0 <- MB_TaloNavCoal[1,2]
MB_TaloNavCoal_0[is.na(MB_TaloNavCoal_0)] <-0
MB_TaloNavCoal_1 <- MB_TaloNavCoal[2,2]
MB_TaloNavCoal_1[is.na(MB_TaloNavCoal_1)] <-0
MB_TaloNavCoal_Total <- na.omit(MB_TaloNavCoal_0) + na.omit(MB_TaloNavCoal_1)
MB_TaloNavCoal_Affected <- round((MB_TaloNavCoal_1/MB_TaloNavCoal_Total)*100, 2)

# Male N and frequency
MB_M_TaloNavCoal <- MB_M %>% count(TaloNavCoal)
MB_M_TaloNavCoal <- 
  MB_M_TaloNavCoal[(MB_M_TaloNavCoal$TaloNavCoal != 6 & MB_M_TaloNavCoal$TaloNavCoal != 9),]
MB_M_TaloNavCoal_0 <- MB_M_TaloNavCoal[1,2]
MB_M_TaloNavCoal_0[is.na(MB_M_TaloNavCoal_0)] <-0
MB_M_TaloNavCoal_1 <- MB_M_TaloNavCoal[2,2]
MB_M_TaloNavCoal_1[is.na(MB_M_TaloNavCoal_1)] <-0
MB_M_TaloNavCoal_Total <- MB_M_TaloNavCoal_0 + MB_M_TaloNavCoal_1
MB_M_TaloNavCoal_Affected <- 
  round((MB_M_TaloNavCoal_1/MB_M_TaloNavCoal_Total)*100, 2)

# Female N and frequency
MB_F_TaloNavCoal <- MB_F %>% count(TaloNavCoal)
MB_F_TaloNavCoal <- 
  MB_F_TaloNavCoal[(MB_F_TaloNavCoal$TaloNavCoal != 6 & MB_F_TaloNavCoal$TaloNavCoal != 9),]
MB_F_TaloNavCoal_0 <- MB_F_TaloNavCoal[1,2]
MB_F_TaloNavCoal_0[is.na(MB_F_TaloNavCoal_0)] <-0
MB_F_TaloNavCoal_1 <- MB_F_TaloNavCoal[2,2]
MB_F_TaloNavCoal_1[is.na(MB_F_TaloNavCoal_1)] <-0
MB_F_TaloNavCoal_Total <- MB_F_TaloNavCoal_0 + MB_F_TaloNavCoal_1
MB_F_TaloNavCoal_Affected <- 
  round((MB_F_TaloNavCoal_1/MB_F_TaloNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_TaloNavCoal_row <- c("TaloNavCoal", MB_TaloNavCoal_Total, MB_TaloNavCoal_1,
                        MB_TaloNavCoal_Affected, MB_M_TaloNavCoal_Total, 
                        MB_M_TaloNavCoal_1, MB_M_TaloNavCoal_Affected, 
                        MB_F_TaloNavCoal_Total, MB_F_TaloNavCoal_1, 
                        MB_F_TaloNavCoal_Affected)

## CF2CF3Coal ##
# Overall N and frequency #
MB_CF2CF3Coal <- MB_wo_probables %>% count(CF2CF3Coal)
MB_CF2CF3Coal <- 
  MB_CF2CF3Coal[(MB_CF2CF3Coal$CF2CF3Coal != 6 & MB_CF2CF3Coal$CF2CF3Coal != 9),]
MB_CF2CF3Coal_0 <- MB_CF2CF3Coal[1,2]
MB_CF2CF3Coal_0[is.na(MB_CF2CF3Coal_0)] <-0
MB_CF2CF3Coal_1 <- MB_CF2CF3Coal[2,2]
MB_CF2CF3Coal_1[is.na(MB_CF2CF3Coal_1)] <-0
MB_CF2CF3Coal_Total <- na.omit(MB_CF2CF3Coal_0) + na.omit(MB_CF2CF3Coal_1)
MB_CF2CF3Coal_Affected <- round((MB_CF2CF3Coal_1/MB_CF2CF3Coal_Total)*100, 2)

# Male N and frequency
MB_M_CF2CF3Coal <- MB_M %>% count(CF2CF3Coal)
MB_M_CF2CF3Coal <- 
  MB_M_CF2CF3Coal[(MB_M_CF2CF3Coal$CF2CF3Coal != 6 & MB_M_CF2CF3Coal$CF2CF3Coal != 9),]
MB_M_CF2CF3Coal_0 <- MB_M_CF2CF3Coal[1,2]
MB_M_CF2CF3Coal_0[is.na(MB_M_CF2CF3Coal_0)] <-0
MB_M_CF2CF3Coal_1 <- MB_M_CF2CF3Coal[2,2]
MB_M_CF2CF3Coal_1[is.na(MB_M_CF2CF3Coal_1)] <-0
MB_M_CF2CF3Coal_Total <- MB_M_CF2CF3Coal_0 + MB_M_CF2CF3Coal_1
MB_M_CF2CF3Coal_Affected <- 
  round((MB_M_CF2CF3Coal_1/MB_M_CF2CF3Coal_Total)*100, 2)

# Female N and frequency
MB_F_CF2CF3Coal <- MB_F %>% count(CF2CF3Coal)
MB_F_CF2CF3Coal <- 
  MB_F_CF2CF3Coal[(MB_F_CF2CF3Coal$CF2CF3Coal != 6 & MB_F_CF2CF3Coal$CF2CF3Coal != 9),]
MB_F_CF2CF3Coal_0 <- MB_F_CF2CF3Coal[1,2]
MB_F_CF2CF3Coal_0[is.na(MB_F_CF2CF3Coal_0)] <-0
MB_F_CF2CF3Coal_1 <- MB_F_CF2CF3Coal[2,2]
MB_F_CF2CF3Coal_1[is.na(MB_F_CF2CF3Coal_1)] <-0
MB_F_CF2CF3Coal_Total <- MB_F_CF2CF3Coal_0 + MB_F_CF2CF3Coal_1
MB_F_CF2CF3Coal_Affected <- 
  round((MB_F_CF2CF3Coal_1/MB_F_CF2CF3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_CF2CF3Coal_row <- c("CF2CF3Coal", MB_CF2CF3Coal_Total, MB_CF2CF3Coal_1,
                       MB_CF2CF3Coal_Affected, MB_M_CF2CF3Coal_Total, 
                       MB_M_CF2CF3Coal_1, MB_M_CF2CF3Coal_Affected, 
                       MB_F_CF2CF3Coal_Total, MB_F_CF2CF3Coal_1, 
                       MB_F_CF2CF3Coal_Affected)

## CF3MT3Coal ##
# Overall N and frequency #
MB_CF3MT3Coal <- MB_wo_probables %>% count(CF3MT3Coal)
MB_CF3MT3Coal <- 
  MB_CF3MT3Coal[(MB_CF3MT3Coal$CF3MT3Coal != 6 & MB_CF3MT3Coal$CF3MT3Coal != 9),]
MB_CF3MT3Coal_0 <- MB_CF3MT3Coal[1,2]
MB_CF3MT3Coal_0[is.na(MB_CF3MT3Coal_0)] <-0
MB_CF3MT3Coal_1 <- MB_CF3MT3Coal[2,2]
MB_CF3MT3Coal_1[is.na(MB_CF3MT3Coal_1)] <-0
MB_CF3MT3Coal_Total <- na.omit(MB_CF3MT3Coal_0) + na.omit(MB_CF3MT3Coal_1)
MB_CF3MT3Coal_Affected <- round((MB_CF3MT3Coal_1/MB_CF3MT3Coal_Total)*100, 2)

# Male N and frequency
MB_M_CF3MT3Coal <- MB_M %>% count(CF3MT3Coal)
MB_M_CF3MT3Coal <- 
  MB_M_CF3MT3Coal[(MB_M_CF3MT3Coal$CF3MT3Coal != 6 & MB_M_CF3MT3Coal$CF3MT3Coal != 9),]
MB_M_CF3MT3Coal_0 <- MB_M_CF3MT3Coal[1,2]
MB_M_CF3MT3Coal_0[is.na(MB_M_CF3MT3Coal_0)] <-0
MB_M_CF3MT3Coal_1 <- MB_M_CF3MT3Coal[2,2]
MB_M_CF3MT3Coal_1[is.na(MB_M_CF3MT3Coal_1)] <-0
MB_M_CF3MT3Coal_Total <- MB_M_CF3MT3Coal_0 + MB_M_CF3MT3Coal_1
MB_M_CF3MT3Coal_Affected <- 
  round((MB_M_CF3MT3Coal_1/MB_M_CF3MT3Coal_Total)*100, 2)

# Female N and frequency
MB_F_CF3MT3Coal <- MB_F %>% count(CF3MT3Coal)
MB_F_CF3MT3Coal <- 
  MB_F_CF3MT3Coal[(MB_F_CF3MT3Coal$CF3MT3Coal != 6 & MB_F_CF3MT3Coal$CF3MT3Coal != 9),]
MB_F_CF3MT3Coal_0 <- MB_F_CF3MT3Coal[1,2]
MB_F_CF3MT3Coal_0[is.na(MB_F_CF3MT3Coal_0)] <-0
MB_F_CF3MT3Coal_1 <- MB_F_CF3MT3Coal[2,2]
MB_F_CF3MT3Coal_1[is.na(MB_F_CF3MT3Coal_1)] <-0
MB_F_CF3MT3Coal_Total <- MB_F_CF3MT3Coal_0 + MB_F_CF3MT3Coal_1
MB_F_CF3MT3Coal_Affected <- 
  round((MB_F_CF3MT3Coal_1/MB_F_CF3MT3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_CF3MT3Coal_row <- c("CF3MT3Coal", MB_CF3MT3Coal_Total, MB_CF3MT3Coal_1, 
                       MB_CF3MT3Coal_Affected, MB_M_CF3MT3Coal_Total, 
                       MB_M_CF3MT3Coal_1, MB_M_CF3MT3Coal_Affected, 
                       MB_F_CF3MT3Coal_Total, MB_F_CF3MT3Coal_1, 
                       MB_F_CF3MT3Coal_Affected)

## OsIntermet ##
# Overall N and frequency #
MB_OsIntermet <- MB_wo_probables %>% count(OsIntermet)
MB_OsIntermet <- 
  MB_OsIntermet[(MB_OsIntermet$OsIntermet != 6 & MB_OsIntermet$OsIntermet != 9),]
MB_OsIntermet_0 <- MB_OsIntermet[1,2]
MB_OsIntermet_0[is.na(MB_OsIntermet_0)] <-0
MB_OsIntermet_1 <- MB_OsIntermet[2,2]
MB_OsIntermet_1[is.na(MB_OsIntermet_1)] <-0
MB_OsIntermet_Total <- na.omit(MB_OsIntermet_0) + na.omit(MB_OsIntermet_1)
MB_OsIntermet_Affected <- round((MB_OsIntermet_1/MB_OsIntermet_Total)*100, 2)

# Male N and frequency
MB_M_OsIntermet <- MB_M %>% count(OsIntermet)
MB_M_OsIntermet <- 
  MB_M_OsIntermet[(MB_M_OsIntermet$OsIntermet != 6 & MB_M_OsIntermet$OsIntermet != 9),]
MB_M_OsIntermet_0 <- MB_M_OsIntermet[1,2]
MB_M_OsIntermet_0[is.na(MB_M_OsIntermet_0)] <-0
MB_M_OsIntermet_1 <- MB_M_OsIntermet[2,2]
MB_M_OsIntermet_1[is.na(MB_M_OsIntermet_1)] <-0
MB_M_OsIntermet_Total <- MB_M_OsIntermet_0 + MB_M_OsIntermet_1
MB_M_OsIntermet_Affected <- 
  round((MB_M_OsIntermet_1/MB_M_OsIntermet_Total)*100, 2)

# Female N and frequency
MB_F_OsIntermet <- MB_F %>% count(OsIntermet)
MB_F_OsIntermet <- 
  MB_F_OsIntermet[(MB_F_OsIntermet$OsIntermet != 6 & MB_F_OsIntermet$OsIntermet != 9),]
MB_F_OsIntermet_0 <- MB_F_OsIntermet[1,2]
MB_F_OsIntermet_0[is.na(MB_F_OsIntermet_0)] <-0
MB_F_OsIntermet_1 <- MB_F_OsIntermet[2,2]
MB_F_OsIntermet_1[is.na(MB_F_OsIntermet_1)] <-0
MB_F_OsIntermet_Total <- MB_F_OsIntermet_0 + MB_F_OsIntermet_1
MB_F_OsIntermet_Affected <- 
  round((MB_F_OsIntermet_1/MB_F_OsIntermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_OsIntermet_row <- c("OsIntermet", MB_OsIntermet_Total, MB_OsIntermet_1, 
                       MB_OsIntermet_Affected, MB_M_OsIntermet_Total, 
                       MB_M_OsIntermet_1, MB_M_OsIntermet_Affected, 
                       MB_F_OsIntermet_Total, MB_F_OsIntermet_1, 
                       MB_F_OsIntermet_Affected)

# OsIntermet unclustered: CF1Intermet #
# Overall N and frequency #
MB_CF1Intermet <- MB_wo_probables %>% count(CF1Intermet)
MB_CF1Intermet <- 
  MB_CF1Intermet[(MB_CF1Intermet$CF1Intermet != 6 & MB_CF1Intermet$CF1Intermet != 9),]
MB_CF1Intermet_0 <- MB_CF1Intermet[1,2]
MB_CF1Intermet_0[is.na(MB_CF1Intermet_0)] <-0
MB_CF1Intermet_1 <- MB_CF1Intermet[2,2]
MB_CF1Intermet_1[is.na(MB_CF1Intermet_1)] <-0
MB_CF1Intermet_Total <- na.omit(MB_CF1Intermet_0) + na.omit(MB_CF1Intermet_1)
MB_CF1Intermet_Affected <- 
  round((MB_CF1Intermet_1/MB_CF1Intermet_Total)*100, 2)


# Male N and frequency
MB_M_CF1Intermet <- MB_M %>% count(CF1Intermet)
MB_M_CF1Intermet <- 
  MB_M_CF1Intermet[(MB_M_CF1Intermet$CF1Intermet != 6 & MB_M_CF1Intermet$CF1Intermet != 9),]
MB_M_CF1Intermet_0 <- MB_M_CF1Intermet[1,2]
MB_M_CF1Intermet_0[is.na(MB_M_CF1Intermet_0)] <-0
MB_M_CF1Intermet_1 <- MB_M_CF1Intermet[2,2]
MB_M_CF1Intermet_1[is.na(MB_M_CF1Intermet_1)] <-0
MB_M_CF1Intermet_Total <- MB_M_CF1Intermet_0 + MB_M_CF1Intermet_1
MB_M_CF1Intermet_Affected <- 
  round((MB_M_CF1Intermet_1/MB_M_CF1Intermet_Total)*100, 2)

# Female N and frequency
MB_F_CF1Intermet <- MB_F %>% count(CF1Intermet)
MB_F_CF1Intermet <- 
  MB_F_CF1Intermet[(MB_F_CF1Intermet$CF1Intermet != 6 & MB_F_CF1Intermet$CF1Intermet != 9),]
MB_F_CF1Intermet_0 <- MB_F_CF1Intermet[1,2]
MB_F_CF1Intermet_0[is.na(MB_F_CF1Intermet_0)] <-0
MB_F_CF1Intermet_1 <- MB_F_CF1Intermet[2,2]
MB_F_CF1Intermet_1[is.na(MB_F_CF1Intermet_1)] <-0
MB_F_CF1Intermet_Total <- MB_F_CF1Intermet_0 + MB_F_CF1Intermet_1
MB_F_CF1Intermet_Affected <- 
  round((MB_F_CF1Intermet_1/MB_F_CF1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_CF1Intermet_row <- c("CF1Intermet", MB_CF1Intermet_Total, MB_CF1Intermet_1, 
                        MB_CF1Intermet_Affected, MB_M_CF1Intermet_Total, 
                        MB_M_CF1Intermet_1, MB_M_CF1Intermet_Affected, 
                        MB_F_CF1Intermet_Total, MB_F_CF1Intermet_1, 
                        MB_F_CF1Intermet_Affected)

# OsIntermet unclustered: MT1Intermet #
# Overall N and frequency #
MB_MT1Intermet <- MB_wo_probables %>% count(MT1Intermet)
MB_MT1Intermet <- 
  MB_MT1Intermet[(MB_MT1Intermet$MT1Intermet != 6 & MB_MT1Intermet$MT1Intermet != 9),]
MB_MT1Intermet_0 <- MB_MT1Intermet[1,2]
MB_MT1Intermet_0[is.na(MB_MT1Intermet_0)] <-0
MB_MT1Intermet_1 <- MB_MT1Intermet[2,2]
MB_MT1Intermet_1[is.na(MB_MT1Intermet_1)] <-0
MB_MT1Intermet_Total <- na.omit(MB_MT1Intermet_0) + na.omit(MB_MT1Intermet_1)
MB_MT1Intermet_Affected <- round((MB_MT1Intermet_1/MB_MT1Intermet_Total)*100, 2)

# Male N and frequency
MB_M_MT1Intermet <- MB_M %>% count(MT1Intermet)
MB_M_MT1Intermet <- 
  MB_M_MT1Intermet[(MB_M_MT1Intermet$MT1Intermet != 6 & MB_M_MT1Intermet$MT1Intermet != 9),]
MB_M_MT1Intermet_0 <- MB_M_MT1Intermet[1,2]
MB_M_MT1Intermet_0[is.na(MB_M_MT1Intermet_0)] <-0
MB_M_MT1Intermet_1 <- MB_M_MT1Intermet[2,2]
MB_M_MT1Intermet_1[is.na(MB_M_MT1Intermet_1)] <-0
MB_M_MT1Intermet_Total <- MB_M_MT1Intermet_0 + MB_M_MT1Intermet_1
MB_M_MT1Intermet_Affected <- 
  round((MB_M_MT1Intermet_1/MB_M_MT1Intermet_Total)*100, 2)

# Female N and frequency
MB_F_MT1Intermet <- MB_F %>% count(MT1Intermet)
MB_F_MT1Intermet <- 
  MB_F_MT1Intermet[(MB_F_MT1Intermet$MT1Intermet != 6 & MB_F_MT1Intermet$MT1Intermet != 9),]
MB_F_MT1Intermet_0 <- MB_F_MT1Intermet[1,2]
MB_F_MT1Intermet_0[is.na(MB_F_MT1Intermet_0)] <-0
MB_F_MT1Intermet_1 <- MB_F_MT1Intermet[2,2]
MB_F_MT1Intermet_1[is.na(MB_F_MT1Intermet_1)] <-0
MB_F_MT1Intermet_Total <- MB_F_MT1Intermet_0 + MB_F_MT1Intermet_1
MB_F_MT1Intermet_Affected <- 
  round((MB_F_MT1Intermet_1/MB_F_MT1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_MT1Intermet_row <- c("MT1Intermet", MB_MT1Intermet_Total, MB_MT1Intermet_1, 
                        MB_MT1Intermet_Affected, MB_M_MT1Intermet_Total, 
                        MB_M_MT1Intermet_1, MB_M_MT1Intermet_Affected, 
                        MB_F_MT1Intermet_Total, MB_F_MT1Intermet_1, 
                        MB_F_MT1Intermet_Affected)

# OsIntermet unclustered: MT2Intermet #
# Overall N and frequency #
MB_MT2Intermet <- MB_wo_probables %>% count(MT2Intermet)
MB_MT2Intermet <- 
  MB_MT2Intermet[(MB_MT2Intermet$MT2Intermet != 6 & MB_MT2Intermet$MT2Intermet != 9),]
MB_MT2Intermet_0 <- MB_MT2Intermet[1,2]
MB_MT2Intermet_0[is.na(MB_MT2Intermet_0)] <-0
MB_MT2Intermet_1 <- MB_MT2Intermet[2,2]
MB_MT2Intermet_1[is.na(MB_MT2Intermet_1)] <-0
MB_MT2Intermet_Total <- na.omit(MB_MT2Intermet_0) + na.omit(MB_MT2Intermet_1)
MB_MT2Intermet_Affected <- round((MB_MT2Intermet_1/MB_MT2Intermet_Total)*100, 2)

# Male N and frequency
MB_M_MT2Intermet <- MB_M %>% count(MT2Intermet)
MB_M_MT2Intermet <- 
  MB_M_MT2Intermet[(MB_M_MT2Intermet$MT2Intermet != 6 & MB_M_MT2Intermet$MT2Intermet != 9),]
MB_M_MT2Intermet_0 <- MB_M_MT2Intermet[1,2]
MB_M_MT2Intermet_0[is.na(MB_M_MT2Intermet_0)] <-0
MB_M_MT2Intermet_1 <- MB_M_MT2Intermet[2,2]
MB_M_MT2Intermet_1[is.na(MB_M_MT2Intermet_1)] <-0
MB_M_MT2Intermet_Total <- MB_M_MT2Intermet_0 + MB_M_MT2Intermet_1
MB_M_MT2Intermet_Affected <- 
  round((MB_M_MT2Intermet_1/MB_M_MT2Intermet_Total)*100, 2)

# Female N and frequency
MB_F_MT2Intermet <- MB_F %>% count(MT2Intermet)
MB_F_MT2Intermet <- 
  MB_F_MT2Intermet[(MB_F_MT2Intermet$MT2Intermet != 6 & MB_F_MT2Intermet$MT2Intermet != 9),]
MB_F_MT2Intermet_0 <- MB_F_MT2Intermet[1,2]
MB_F_MT2Intermet_0[is.na(MB_F_MT2Intermet_0)] <-0
MB_F_MT2Intermet_1 <- MB_F_MT2Intermet[2,2]
MB_F_MT2Intermet_1[is.na(MB_F_MT2Intermet_1)] <-0
MB_F_MT2Intermet_Total <- MB_F_MT2Intermet_0 + MB_F_MT2Intermet_1
MB_F_MT2Intermet_Affected <- 
  round((MB_F_MT2Intermet_1/MB_F_MT2Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_MT2Intermet_row <- c("MT2Intermet", MB_MT2Intermet_Total, MB_MT2Intermet_1, 
                        MB_MT2Intermet_Affected, MB_M_MT2Intermet_Total, 
                        MB_M_MT2Intermet_1, MB_M_MT2Intermet_Affected, 
                        MB_F_MT2Intermet_Total, MB_F_MT2Intermet_1, 
                        MB_F_MT2Intermet_Affected)


## AnkleCoal ##
# Overall N and frequency #
MB_AnkleCoal <- MB_wo_probables %>% count(AnkleCoal)
MB_AnkleCoal <- 
  MB_AnkleCoal[(MB_AnkleCoal$AnkleCoal != 6 & MB_AnkleCoal$AnkleCoal != 9),]
MB_AnkleCoal_0 <- MB_AnkleCoal[1,2]; MB_AnkleCoal_0[is.na(MB_AnkleCoal_0)] <-0
MB_AnkleCoal_1 <- MB_AnkleCoal[2,2]; MB_AnkleCoal_1[is.na(MB_AnkleCoal_1)] <-0
MB_AnkleCoal_Total <- na.omit(MB_AnkleCoal_0) + na.omit(MB_AnkleCoal_1)
MB_AnkleCoal_Affected <- round((MB_AnkleCoal_1/MB_AnkleCoal_Total)*100, 2)

# Male N and frequency
MB_M_AnkleCoal <- MB_M %>% count(AnkleCoal)
MB_M_AnkleCoal <- 
  MB_M_AnkleCoal[(MB_M_AnkleCoal$AnkleCoal != 6 & MB_M_AnkleCoal$AnkleCoal != 9),]
MB_M_AnkleCoal_0 <- MB_M_AnkleCoal[1,2]
MB_M_AnkleCoal_0[is.na(MB_M_AnkleCoal_0)] <-0
MB_M_AnkleCoal_1 <- MB_M_AnkleCoal[2,2]
MB_M_AnkleCoal_1[is.na(MB_M_AnkleCoal_1)] <-0
MB_M_AnkleCoal_Total <- MB_M_AnkleCoal_0 + MB_M_AnkleCoal_1
MB_M_AnkleCoal_Affected <- 
  round((MB_M_AnkleCoal_1/MB_M_AnkleCoal_Total)*100, 2)

# Female N and frequency
MB_F_AnkleCoal <- MB_F %>% count(AnkleCoal)
MB_F_AnkleCoal <- 
  MB_F_AnkleCoal[(MB_F_AnkleCoal$AnkleCoal != 6 & MB_F_AnkleCoal$AnkleCoal != 9),]
MB_F_AnkleCoal_0 <- MB_F_AnkleCoal[1,2]
MB_F_AnkleCoal_0[is.na(MB_F_AnkleCoal_0)] <-0
MB_F_AnkleCoal_1 <- MB_F_AnkleCoal[2,2]
MB_F_AnkleCoal_1[is.na(MB_F_AnkleCoal_1)] <-0
MB_F_AnkleCoal_Total <- MB_F_AnkleCoal_0 + MB_F_AnkleCoal_1
MB_F_AnkleCoal_Affected <- 
  round((MB_F_AnkleCoal_1/MB_F_AnkleCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_AnkleCoal_row <- c("AnkleCoal", MB_AnkleCoal_Total, MB_AnkleCoal_1,
                      MB_AnkleCoal_Affected, MB_M_AnkleCoal_Total, 
                      MB_M_AnkleCoal_1, MB_M_AnkleCoal_Affected,  
                      MB_F_AnkleCoal_Total, MB_F_AnkleCoal_1, 
                      MB_F_AnkleCoal_Affected)

# AnkleCoal unclustered: CalcNavCoal #
# Overall N and frequency #
MB_CalcNavCoal <- MB_wo_probables %>% count(CalcNavCoal)
MB_CalcNavCoal <- 
  MB_CalcNavCoal[(MB_CalcNavCoal$CalcNavCoal != 6 & MB_CalcNavCoal$CalcNavCoal != 9),]
MB_CalcNavCoal_0 <- MB_CalcNavCoal[1,2]
MB_CalcNavCoal_0[is.na(MB_CalcNavCoal_0)] <-0
MB_CalcNavCoal_1 <- MB_CalcNavCoal[2,2]
MB_CalcNavCoal_1[is.na(MB_CalcNavCoal_1)] <-0
MB_CalcNavCoal_Total <- na.omit(MB_CalcNavCoal_0) + na.omit(MB_CalcNavCoal_1)
MB_CalcNavCoal_Affected <- round((MB_CalcNavCoal_1/MB_CalcNavCoal_Total)*100, 2)

# Male N and frequency
MB_M_CalcNavCoal <- MB_M %>% count(CalcNavCoal)
MB_M_CalcNavCoal <- 
  MB_M_CalcNavCoal[(MB_M_CalcNavCoal$CalcNavCoal != 6 & MB_M_CalcNavCoal$CalcNavCoal != 9),]
MB_M_CalcNavCoal_0 <- MB_M_CalcNavCoal[1,2]
MB_M_CalcNavCoal_0[is.na(MB_M_CalcNavCoal_0)] <-0
MB_M_CalcNavCoal_1 <- MB_M_CalcNavCoal[2,2]
MB_M_CalcNavCoal_1[is.na(MB_M_CalcNavCoal_1)] <-0
MB_M_CalcNavCoal_Total <- MB_M_CalcNavCoal_0 + MB_M_CalcNavCoal_1
MB_M_CalcNavCoal_Affected <- 
  round((MB_M_CalcNavCoal_1/MB_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
MB_F_CalcNavCoal <- MB_F %>% count(CalcNavCoal)
MB_F_CalcNavCoal <- 
  MB_F_CalcNavCoal[(MB_F_CalcNavCoal$CalcNavCoal != 6 & MB_F_CalcNavCoal$CalcNavCoal != 9),]
MB_F_CalcNavCoal_0 <- MB_F_CalcNavCoal[1,2]
MB_F_CalcNavCoal_0[is.na(MB_F_CalcNavCoal_0)] <-0
MB_F_CalcNavCoal_1 <- MB_F_CalcNavCoal[2,2]
MB_F_CalcNavCoal_1[is.na(MB_F_CalcNavCoal_1)] <-0
MB_F_CalcNavCoal_Total <- MB_F_CalcNavCoal_0 + MB_F_CalcNavCoal_1
MB_F_CalcNavCoal_Affected <- 
  round((MB_F_CalcNavCoal_1/MB_F_CalcNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_CalcNavCoal_row <- c("CalcNavCoal", MB_CalcNavCoal_Total, MB_CalcNavCoal_1, 
                        MB_CalcNavCoal_Affected, MB_M_CalcNavCoal_Total, 
                        MB_M_CalcNavCoal_1, MB_M_CalcNavCoal_Affected, 
                        MB_F_CalcNavCoal_Total, MB_F_CalcNavCoal_1, 
                        MB_F_CalcNavCoal_Affected)


# AnkleCoal unclustered: TaloCalcCoal #
# Overall N and frequency #
MB_TaloCalcCoal <- MB_wo_probables %>% count(TaloCalcCoal)
MB_TaloCalcCoal <- 
  MB_TaloCalcCoal[(MB_TaloCalcCoal$TaloCalcCoal != 6 & MB_TaloCalcCoal$TaloCalcCoal != 9),]
MB_TaloCalcCoal_0 <- MB_TaloCalcCoal[1,2]
MB_TaloCalcCoal_0[is.na(MB_TaloCalcCoal_0)] <-0
MB_TaloCalcCoal_1 <- MB_TaloCalcCoal[2,2]
MB_TaloCalcCoal_1[is.na(MB_TaloCalcCoal_1)] <-0
MB_TaloCalcCoal_Total <- na.omit(MB_TaloCalcCoal_0) + na.omit(MB_TaloCalcCoal_1)
MB_TaloCalcCoal_Affected <- 
  round((MB_TaloCalcCoal_1/MB_TaloCalcCoal_Total)*100, 2)

# Male N and frequency
MB_M_TaloCalcCoal <- MB_M %>% count(TaloCalcCoal)
MB_M_TaloCalcCoal <- 
  MB_M_TaloCalcCoal[(MB_M_TaloCalcCoal$TaloCalcCoal != 6 & MB_M_TaloCalcCoal$TaloCalcCoal != 9),]
MB_M_TaloCalcCoal_0 <- MB_M_TaloCalcCoal[1,2]
MB_M_TaloCalcCoal_0[is.na(MB_M_TaloCalcCoal_0)] <-0
MB_M_TaloCalcCoal_1 <- MB_M_TaloCalcCoal[2,2]
MB_M_TaloCalcCoal_1[is.na(MB_M_TaloCalcCoal_1)] <-0
MB_M_TaloCalcCoal_Total <- MB_M_TaloCalcCoal_0 + MB_M_TaloCalcCoal_1
MB_M_TaloCalcCoal_Affected <- 
  round((MB_M_TaloCalcCoal_1/MB_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
MB_F_TaloCalcCoal <- MB_F %>% count(TaloCalcCoal)
MB_F_TaloCalcCoal <- 
  MB_F_TaloCalcCoal[(MB_F_TaloCalcCoal$TaloCalcCoal != 6 & MB_F_TaloCalcCoal$TaloCalcCoal != 9),]
MB_F_TaloCalcCoal_0 <- MB_F_TaloCalcCoal[1,2]
MB_F_TaloCalcCoal_0[is.na(MB_F_TaloCalcCoal_0)] <-0
MB_F_TaloCalcCoal_1 <- MB_F_TaloCalcCoal[2,2]
MB_F_TaloCalcCoal_1[is.na(MB_F_TaloCalcCoal_1)] <-0
MB_F_TaloCalcCoal_Total <- MB_F_TaloCalcCoal_0 + MB_F_TaloCalcCoal_1
MB_F_TaloCalcCoal_Affected <- 
  round((MB_F_TaloCalcCoal_1/MB_F_TaloCalcCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
MB_TaloCalcCoal_row <- c("TaloCalcCoal", MB_TaloCalcCoal_Total, MB_TaloCalcCoal_1, 
                         MB_TaloCalcCoal_Affected, MB_M_TaloCalcCoal_Total, 
                         MB_M_TaloCalcCoal_1, MB_M_TaloCalcCoal_Affected, 
                         MB_F_TaloCalcCoal_Total, MB_F_TaloCalcCoal_1, 
                         MB_F_TaloCalcCoal_Affected)


## Table with overview of trait frequencies for MB ##
MB_freq_head <- c("Anomaly", "Overall N", "Overall Affected", "Overall %", 
                  "Male N", "Male affected", "Male %",
                  "Female N", "Female affected", "Female %")
MB_freq <- rbind(MB_AccessNav_row, MB_BrachyD_row, MB_BrachyMT1_row,
                 MB_BrachyMT4_row, MB_BrachyPP1_row, MB_CalcCubCoal_row,
                 MB_TaloNavCoal_row, MB_CF2CF3Coal_row, MB_CF3MT3Coal_row,
                 MB_OsIntermet_row, MB_CF1Intermet_row, MB_MT1Intermet_row,
                 MB_MT2Intermet_row, MB_AnkleCoal_row, MB_CalcNavCoal_row, 
                 MB_TaloCalcCoal_row)
MB_freq <- data.frame(rbind(MB_freq_head, MB_freq))
MB_freq <- row_to_names(MB_freq, 1, remove_row = T)


# Frequency calculations: reference sample --------------------------------


#                                                                    #
## Calculating the frequency of the traits in the reference sample  ##
#                                                                    #

### Arnhem collection ###

## AccessNav ##
# Overall N and frequency #
AR_AccessNav <- AR_wo_probables %>% count(AccessNav)
AR_AccessNav <- 
  AR_AccessNav[(AR_AccessNav$AccessNav != 6 & AR_AccessNav$AccessNav != 9),]
AR_AccessNav_0 <- AR_AccessNav[1,2]
AR_AccessNav_0[is.na(AR_AccessNav_0)] <-0
AR_AccessNav_1 <- AR_AccessNav[2,2]
AR_AccessNav_1[is.na(AR_AccessNav_1)] <-0
AR_AccessNav_Total <- AR_AccessNav_0 + AR_AccessNav_1
AR_AccessNav_Affected <- round((AR_AccessNav_1/AR_AccessNav_Total)*100, 2)

# Male N and frequency
AR_M_AccessNav <- AR_M %>% count(AccessNav)
AR_M_AccessNav <- 
  AR_M_AccessNav[(AR_M_AccessNav$AccessNav != 6 & AR_M_AccessNav$AccessNav != 9),]
AR_M_AccessNav_0 <- AR_M_AccessNav[1,2]
AR_M_AccessNav_0[is.na(AR_M_AccessNav_0)] <-0
AR_M_AccessNav_1 <- AR_M_AccessNav[2,2]
AR_M_AccessNav_1[is.na(AR_M_AccessNav_1)] <-0
AR_M_AccessNav_Total <- AR_M_AccessNav_0 + AR_M_AccessNav_1
AR_M_AccessNav_Affected <- round((AR_M_AccessNav_1/AR_M_AccessNav_Total)*100, 2)

# Female N and frequency
AR_F_AccessNav <- AR_F %>% count(AccessNav)
AR_F_AccessNav <- 
  AR_F_AccessNav[(AR_F_AccessNav$AccessNav != 6 & AR_F_AccessNav$AccessNav != 9),]
AR_F_AccessNav_0 <- AR_F_AccessNav[1,2]
AR_F_AccessNav_0[is.na(AR_F_AccessNav_0)] <-0
AR_F_AccessNav_1 <- AR_F_AccessNav[2,2]
AR_F_AccessNav_1[is.na(AR_F_AccessNav_1)] <-0
AR_F_AccessNav_Total <- AR_F_AccessNav_0 + AR_F_AccessNav_1
AR_F_AccessNav_Affected <- round((AR_F_AccessNav_1/AR_F_AccessNav_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_AccessNav_row <- c("AccessNav", AR_AccessNav_Total, AR_AccessNav_1, 
                      AR_AccessNav_Affected, AR_M_AccessNav_Total, 
                      AR_M_AccessNav_1, AR_M_AccessNav_Affected, 
                      AR_F_AccessNav_Total, AR_F_AccessNav_1, 
                      AR_F_AccessNav_Affected)


## BrachyD ##
# Overall N and frequency #
AR_BrachyD <- AR_wo_probables %>% count(BrachyD)
AR_BrachyD <- AR_BrachyD[(AR_BrachyD$BrachyD != 6 & AR_BrachyD$BrachyD != 9),]
AR_BrachyD_0 <- AR_BrachyD[1,2]; AR_BrachyD_0[is.na(AR_BrachyD_0)] <-0
AR_BrachyD_1 <- AR_BrachyD[2,2]; AR_BrachyD_1[is.na(AR_BrachyD_1)] <-0
AR_BrachyD_Total <- AR_BrachyD_0 + AR_BrachyD_1
AR_BrachyD_Affected <- round((AR_BrachyD_1/AR_BrachyD_Total)*100, 2)

# Male N and frequency
AR_M_BrachyD <- AR_M %>% count(BrachyD)
AR_M_BrachyD <- 
  AR_M_BrachyD[(AR_M_BrachyD$BrachyD != 6 & AR_M_BrachyD$BrachyD != 9),]
AR_M_BrachyD_0 <- AR_M_BrachyD[1,2]; AR_M_BrachyD_0[is.na(AR_M_BrachyD_0)] <-0
AR_M_BrachyD_1 <- AR_M_BrachyD[2,2]; AR_M_BrachyD_1[is.na(AR_M_BrachyD_1)] <-0
AR_M_BrachyD_Total <- AR_M_BrachyD_0 + AR_M_BrachyD_1
AR_M_BrachyD_Affected <- round((AR_M_BrachyD_1/AR_M_BrachyD_Total)*100, 2)

# Female N and frequency
AR_F_BrachyD <- AR_F %>% count(BrachyD)
AR_F_BrachyD <- 
  AR_F_BrachyD[(AR_F_BrachyD$BrachyD != 6 & AR_F_BrachyD$BrachyD != 9),]
AR_F_BrachyD_0 <- AR_F_BrachyD[1,2]; AR_F_BrachyD_0[is.na(AR_F_BrachyD_0)] <-0
AR_F_BrachyD_1 <- AR_F_BrachyD[2,2]; AR_F_BrachyD_1[is.na(AR_F_BrachyD_1)] <-0
AR_F_BrachyD_Total <- AR_F_BrachyD_0 + AR_F_BrachyD_1
AR_F_BrachyD_Affected <- round((AR_F_BrachyD_1/AR_F_BrachyD_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_BrachyD_row <- c("BrachyD", AR_BrachyD_Total, AR_BrachyD_1, 
                    AR_BrachyD_Affected, AR_M_BrachyD_Total, 
                    AR_M_BrachyD_1, AR_M_BrachyD_Affected, 
                    AR_F_BrachyD_Total, AR_F_BrachyD_1, 
                    AR_F_BrachyD_Affected)


## BrachyMT1 ##
# Overall N and frequency #
AR_BrachyMT1 <- AR_wo_probables %>% count(BrachyMT1)
AR_BrachyMT1 <- 
  AR_BrachyMT1[(AR_BrachyMT1$BrachyMT1 != 6 & AR_BrachyMT1$BrachyMT1 != 9),]
AR_BrachyMT1_0 <- AR_BrachyMT1[1,2]; AR_BrachyMT1_0[is.na(AR_BrachyMT1_0)] <-0
AR_BrachyMT1_1 <- AR_BrachyMT1[2,2]; AR_BrachyMT1_1[is.na(AR_BrachyMT1_1)] <-0
AR_BrachyMT1_Total <- AR_BrachyMT1_0 + AR_BrachyMT1_1
AR_BrachyMT1_Affected <- round((AR_BrachyMT1_1/AR_BrachyMT1_Total)*100, 2)

# Male N and frequency
AR_M_BrachyMT1 <- AR_M %>% count(BrachyMT1)
AR_M_BrachyMT1 <- 
  AR_M_BrachyMT1[(AR_M_BrachyMT1$BrachyMT1 != 6 & AR_M_BrachyMT1$BrachyMT1 != 9),]
AR_M_BrachyMT1_0 <- AR_M_BrachyMT1[1,2]
AR_M_BrachyMT1_0[is.na(AR_M_BrachyMT1_0)] <-0
AR_M_BrachyMT1_1 <- AR_M_BrachyMT1[2,2]
AR_M_BrachyMT1_1[is.na(AR_M_BrachyMT1_1)] <-0
AR_M_BrachyMT1_Total <- AR_M_BrachyMT1_0 + AR_M_BrachyMT1_1
AR_M_BrachyMT1_Affected <- round((AR_M_BrachyMT1_1/AR_M_BrachyMT1_Total)*100, 2)

# Female N and frequency
AR_F_BrachyMT1 <- AR_F %>% count(BrachyMT1)
AR_F_BrachyMT1 <- 
  AR_F_BrachyMT1[(AR_F_BrachyMT1$BrachyMT1 != 6 & AR_F_BrachyMT1$BrachyMT1 != 9),]
AR_F_BrachyMT1_0 <- AR_F_BrachyMT1[1,2]
AR_F_BrachyMT1_0[is.na(AR_F_BrachyMT1_0)] <-0
AR_F_BrachyMT1_1 <- AR_F_BrachyMT1[2,2]
AR_F_BrachyMT1_1[is.na(AR_F_BrachyMT1_1)] <-0
AR_F_BrachyMT1_Total <- AR_F_BrachyMT1_0 + AR_F_BrachyMT1_1
AR_F_BrachyMT1_Affected <- round((AR_F_BrachyMT1_1/AR_F_BrachyMT1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_BrachyMT1_row <- c("BrachyMT1", AR_BrachyMT1_Total, AR_BrachyMT1_1, 
                      AR_BrachyMT1_Affected, AR_M_BrachyMT1_Total, 
                      AR_M_BrachyMT1_1, AR_M_BrachyMT1_Affected, 
                      AR_F_BrachyMT1_Total, AR_F_BrachyMT1_1, 
                      AR_F_BrachyMT1_Affected)

## BrachyMT4 ##
# Overall N and frequency #
AR_BrachyMT4 <- AR_wo_probables %>% count(BrachyMT4)
AR_BrachyMT4 <- 
  AR_BrachyMT4[(AR_BrachyMT4$BrachyMT4 != 6 & AR_BrachyMT4$BrachyMT4 != 9),]
AR_BrachyMT4_0 <- AR_BrachyMT4[1,2]
AR_BrachyMT4_0[is.na(AR_BrachyMT4_0)] <-0
AR_BrachyMT4_1 <- AR_BrachyMT4[2,2]
AR_BrachyMT4_1[is.na(AR_BrachyMT4_1)] <-0
AR_BrachyMT4_Total <- AR_BrachyMT4_0 + AR_BrachyMT4_1
AR_BrachyMT4_Affected <- round((AR_BrachyMT4_1/AR_BrachyMT4_Total)*100, 2)

# Male N and frequency
AR_M_BrachyMT4 <- AR_M %>% count(BrachyMT4)
AR_M_BrachyMT4 <- 
  AR_M_BrachyMT4[(AR_M_BrachyMT4$BrachyMT4 != 6 & AR_M_BrachyMT4$BrachyMT4 != 9),]
AR_M_BrachyMT4_0 <- AR_M_BrachyMT4[1,2]
AR_M_BrachyMT4_0[is.na(AR_M_BrachyMT4_0)] <-0
AR_M_BrachyMT4_1 <- AR_M_BrachyMT4[2,2]
AR_M_BrachyMT4_1[is.na(AR_M_BrachyMT4_1)] <-0
AR_M_BrachyMT4_Total <- AR_M_BrachyMT4_0 + AR_M_BrachyMT4_1
AR_M_BrachyMT4_Affected <- round((AR_M_BrachyMT4_1/AR_M_BrachyMT4_Total)*100, 2)

# Female N and frequency
AR_F_BrachyMT4 <- AR_F %>% count(BrachyMT4)
AR_F_BrachyMT4 <- 
  AR_F_BrachyMT4[(AR_F_BrachyMT4$BrachyMT4 != 6 & AR_F_BrachyMT4$BrachyMT4 != 9),]
AR_F_BrachyMT4_0 <- AR_F_BrachyMT4[1,2]
AR_F_BrachyMT4_0[is.na(AR_F_BrachyMT4_0)] <-0
AR_F_BrachyMT4_1 <- AR_F_BrachyMT4[2,2]
AR_F_BrachyMT4_1[is.na(AR_F_BrachyMT4_1)] <-0
AR_F_BrachyMT4_Total <- AR_F_BrachyMT4_0 + AR_F_BrachyMT4_1
AR_F_BrachyMT4_Affected <- round((AR_F_BrachyMT4_1/AR_F_BrachyMT4_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_BrachyMT4_row <- c("BrachyMT4", AR_BrachyMT4_Total, AR_BrachyMT4_1, 
                      AR_BrachyMT4_Affected, AR_M_BrachyMT4_Total, 
                      AR_M_BrachyMT4_1, AR_M_BrachyMT4_Affected, 
                      AR_F_BrachyMT4_Total, AR_F_BrachyMT4_1, 
                      AR_F_BrachyMT4_Affected)

## BrachyPP1 ##
# Overall N and frequency #
AR_BrachyPP1 <- AR_wo_probables %>% count(BrachyPP1)
AR_BrachyPP1 <- 
  AR_BrachyPP1[(AR_BrachyPP1$BrachyPP1 != 6 & AR_BrachyPP1$BrachyPP1 != 9),]
AR_BrachyPP1_0 <- AR_BrachyPP1[1,2]; AR_BrachyPP1_0[is.na(AR_BrachyPP1_0)] <-0
AR_BrachyPP1_1 <- AR_BrachyPP1[2,2]; AR_BrachyPP1_1[is.na(AR_BrachyPP1_1)] <-0
AR_BrachyPP1_Total <- AR_BrachyPP1_0 + AR_BrachyPP1_1
AR_BrachyPP1_Affected <- round((AR_BrachyPP1_1/AR_BrachyPP1_Total)*100, 2)

# Male N and frequency
AR_M_BrachyPP1 <- AR_M %>% count(BrachyPP1)
AR_M_BrachyPP1 <- 
  AR_M_BrachyPP1[(AR_M_BrachyPP1$BrachyPP1 != 6 & AR_M_BrachyPP1$BrachyPP1 != 9),]
AR_M_BrachyPP1_0 <- AR_M_BrachyPP1[1,2]
AR_M_BrachyPP1_0[is.na(AR_M_BrachyPP1_0)] <-0
AR_M_BrachyPP1_1 <- AR_M_BrachyPP1[2,2]
AR_M_BrachyPP1_1[is.na(AR_M_BrachyPP1_1)] <-0
AR_M_BrachyPP1_Total <- AR_M_BrachyPP1_0 + AR_M_BrachyPP1_1
AR_M_BrachyPP1_Affected <- round((AR_M_BrachyPP1_1/AR_M_BrachyPP1_Total)*100, 2)

# Female N and frequency
AR_F_BrachyPP1 <- AR_F %>% count(BrachyPP1)
AR_F_BrachyPP1 <- 
  AR_F_BrachyPP1[(AR_F_BrachyPP1$BrachyPP1 != 6 & AR_F_BrachyPP1$BrachyPP1 != 9),]
AR_F_BrachyPP1_0 <- AR_F_BrachyPP1[1,2]
AR_F_BrachyPP1_0[is.na(AR_F_BrachyPP1_0)] <-0
AR_F_BrachyPP1_1 <- AR_F_BrachyPP1[2,2]
AR_F_BrachyPP1_1[is.na(AR_F_BrachyPP1_1)] <-0
AR_F_BrachyPP1_Total <- AR_F_BrachyPP1_0 + AR_F_BrachyPP1_1
AR_F_BrachyPP1_Affected <- round((AR_F_BrachyPP1_1/AR_F_BrachyPP1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_BrachyPP1_row <- c("BrachyPP1", AR_BrachyPP1_Total, AR_BrachyPP1_1, 
                      AR_BrachyPP1_Affected, AR_M_BrachyPP1_Total, 
                      AR_M_BrachyPP1_1, AR_M_BrachyPP1_Affected, 
                      AR_F_BrachyPP1_Total, AR_F_BrachyPP1_1, 
                      AR_F_BrachyPP1_Affected)

## CalcCubCoal ##
# Overall N and frequency #
AR_CalcCubCoal <- AR_wo_probables %>% count(CalcCubCoal)
AR_CalcCubCoal <- 
  AR_CalcCubCoal[(AR_CalcCubCoal$CalcCubCoal != 6 & AR_CalcCubCoal$CalcCubCoal != 9),]
AR_CalcCubCoal_0 <- AR_CalcCubCoal[1,2]
AR_CalcCubCoal_0[is.na(AR_CalcCubCoal_0)] <-0
AR_CalcCubCoal_1 <- AR_CalcCubCoal[2,2]
AR_CalcCubCoal_1[is.na(AR_CalcCubCoal_1)] <-0
AR_CalcCubCoal_Total <- AR_CalcCubCoal_0 + AR_CalcCubCoal_1
AR_CalcCubCoal_Affected <- 
  round((AR_CalcCubCoal_1/AR_CalcCubCoal_Total)*100, 2)

# Male N and frequency
AR_M_CalcCubCoal <- AR_M %>% count(CalcCubCoal)
AR_M_CalcCubCoal <- 
  AR_M_CalcCubCoal[(AR_M_CalcCubCoal$CalcCubCoal != 6 & AR_M_CalcCubCoal$CalcCubCoal != 9),]
AR_M_CalcCubCoal_0 <- AR_M_CalcCubCoal[1,2]
AR_M_CalcCubCoal_0[is.na(AR_M_CalcCubCoal_0)] <-0
AR_M_CalcCubCoal_1 <- AR_M_CalcCubCoal[2,2]
AR_M_CalcCubCoal_1[is.na(AR_M_CalcCubCoal_1)] <-0
AR_M_CalcCubCoal_Total <- AR_M_CalcCubCoal_0 + AR_M_CalcCubCoal_1
AR_M_CalcCubCoal_Affected <- 
  round((AR_M_CalcCubCoal_1/AR_M_CalcCubCoal_Total)*100, 2)

# Female N and frequency
AR_F_CalcCubCoal <- AR_F %>% count(CalcCubCoal)
AR_F_CalcCubCoal <- 
  AR_F_CalcCubCoal[(AR_F_CalcCubCoal$CalcCubCoal != 6 & AR_F_CalcCubCoal$CalcCubCoal != 9),]
AR_F_CalcCubCoal_0 <- AR_F_CalcCubCoal[1,2]
AR_F_CalcCubCoal_0[is.na(AR_F_CalcCubCoal_0)] <-0
AR_F_CalcCubCoal_1 <- AR_F_CalcCubCoal[2,2]
AR_F_CalcCubCoal_1[is.na(AR_F_CalcCubCoal_1)] <-0
AR_F_CalcCubCoal_Total <- AR_F_CalcCubCoal_0 + AR_F_CalcCubCoal_1
AR_F_CalcCubCoal_Affected <- 
  round((AR_F_CalcCubCoal_1/AR_F_CalcCubCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_CalcCubCoal_row <- c("CalcCubCoal", AR_CalcCubCoal_Total, AR_CalcCubCoal_1, 
                        AR_CalcCubCoal_Affected, AR_M_CalcCubCoal_Total, 
                        AR_M_CalcCubCoal_1, AR_M_CalcCubCoal_Affected, 
                        AR_F_CalcCubCoal_Total, AR_F_CalcCubCoal_1, 
                        AR_F_CalcCubCoal_Affected)

## TaloNavCoal ##
# Overall N and frequency #
AR_TaloNavCoal <- AR_wo_probables %>% count(TaloNavCoal)
AR_TaloNavCoal <- 
  AR_TaloNavCoal[(AR_TaloNavCoal$TaloNavCoal != 6 & AR_TaloNavCoal$TaloNavCoal != 9),]
AR_TaloNavCoal_0 <- AR_TaloNavCoal[1,2]
AR_TaloNavCoal_0[is.na(AR_TaloNavCoal_0)] <-0
AR_TaloNavCoal_1 <- AR_TaloNavCoal[2,2]
AR_TaloNavCoal_1[is.na(AR_TaloNavCoal_1)] <-0
AR_TaloNavCoal_Total <- AR_TaloNavCoal_0 + AR_TaloNavCoal_1
AR_TaloNavCoal_Affected <- round((AR_TaloNavCoal_1/AR_TaloNavCoal_Total)*100, 2)

# Male N and frequency
AR_M_TaloNavCoal <- AR_M %>% count(TaloNavCoal)
AR_M_TaloNavCoal <- 
  AR_M_TaloNavCoal[(AR_M_TaloNavCoal$TaloNavCoal != 6 & AR_M_TaloNavCoal$TaloNavCoal != 9),]
AR_M_TaloNavCoal_0 <- AR_M_TaloNavCoal[1,2]
AR_M_TaloNavCoal_0[is.na(AR_M_TaloNavCoal_0)] <-0
AR_M_TaloNavCoal_1 <- AR_M_TaloNavCoal[2,2]
AR_M_TaloNavCoal_1[is.na(AR_M_TaloNavCoal_1)] <-0
AR_M_TaloNavCoal_Total <- AR_M_TaloNavCoal_0 + AR_M_TaloNavCoal_1
AR_M_TaloNavCoal_Affected <- 
  round((AR_M_TaloNavCoal_1/AR_M_TaloNavCoal_Total)*100, 2)

# Female N and frequency
AR_F_TaloNavCoal <- AR_F %>% count(TaloNavCoal)
AR_F_TaloNavCoal <- 
  AR_F_TaloNavCoal[(AR_F_TaloNavCoal$TaloNavCoal != 6 & AR_F_TaloNavCoal$TaloNavCoal != 9),]
AR_F_TaloNavCoal_0 <- AR_F_TaloNavCoal[1,2]
AR_F_TaloNavCoal_0[is.na(AR_F_TaloNavCoal_0)] <-0
AR_F_TaloNavCoal_1 <- AR_F_TaloNavCoal[2,2]
AR_F_TaloNavCoal_1[is.na(AR_F_TaloNavCoal_1)] <-0
AR_F_TaloNavCoal_Total <- AR_F_TaloNavCoal_0 + AR_F_TaloNavCoal_1
AR_F_TaloNavCoal_Affected <- 
  round((AR_F_TaloNavCoal_1/AR_F_TaloNavCoal_Total)*100, 2)

# overview of Taloulated frequencies for table (see below) #
AR_TaloNavCoal_row <- c("TaloNavCoal", AR_TaloNavCoal_Total, AR_TaloNavCoal_1, 
                        AR_TaloNavCoal_Affected, AR_M_TaloNavCoal_Total, 
                        AR_M_TaloNavCoal_1, AR_M_TaloNavCoal_Affected, 
                        AR_F_TaloNavCoal_Total, AR_F_TaloNavCoal_1, 
                        AR_F_TaloNavCoal_Affected)

## CF2CF3Coal ##
# Overall N and frequency #
AR_CF2CF3Coal <- AR_wo_probables %>% count(CF2CF3Coal)
AR_CF2CF3Coal <- 
  AR_CF2CF3Coal[(AR_CF2CF3Coal$CF2CF3Coal != 6 & AR_CF2CF3Coal$CF2CF3Coal != 9),]
AR_CF2CF3Coal_0 <- AR_CF2CF3Coal[1,2]
AR_CF2CF3Coal_0[is.na(AR_CF2CF3Coal_0)] <-0
AR_CF2CF3Coal_1 <- AR_CF2CF3Coal[2,2]
AR_CF2CF3Coal_1[is.na(AR_CF2CF3Coal_1)] <-0
AR_CF2CF3Coal_Total <- AR_CF2CF3Coal_0 + AR_CF2CF3Coal_1
AR_CF2CF3Coal_Affected <- round((AR_CF2CF3Coal_1/AR_CF2CF3Coal_Total)*100, 2)

# Male N and frequency
AR_M_CF2CF3Coal <- AR_M %>% count(CF2CF3Coal)
AR_M_CF2CF3Coal <- 
  AR_M_CF2CF3Coal[(AR_M_CF2CF3Coal$CF2CF3Coal != 6 & AR_M_CF2CF3Coal$CF2CF3Coal != 9),]
AR_M_CF2CF3Coal_0 <- AR_M_CF2CF3Coal[1,2]
AR_M_CF2CF3Coal_0[is.na(AR_M_CF2CF3Coal_0)] <-0
AR_M_CF2CF3Coal_1 <- AR_M_CF2CF3Coal[2,2]
AR_M_CF2CF3Coal_1[is.na(AR_M_CF2CF3Coal_1)] <-0
AR_M_CF2CF3Coal_Total <- AR_M_CF2CF3Coal_0 + AR_M_CF2CF3Coal_1
AR_M_CF2CF3Coal_Affected <- 
  round((AR_M_CF2CF3Coal_1/AR_M_CF2CF3Coal_Total)*100, 2)

# Female N and frequency
AR_F_CF2CF3Coal <- AR_F %>% count(CF2CF3Coal)
AR_F_CF2CF3Coal <- 
  AR_F_CF2CF3Coal[(AR_F_CF2CF3Coal$CF2CF3Coal != 6 & AR_F_CF2CF3Coal$CF2CF3Coal != 9),]
AR_F_CF2CF3Coal_0 <- AR_F_CF2CF3Coal[1,2]
AR_F_CF2CF3Coal_0[is.na(AR_F_CF2CF3Coal_0)] <-0
AR_F_CF2CF3Coal_1 <- AR_F_CF2CF3Coal[2,2]
AR_F_CF2CF3Coal_1[is.na(AR_F_CF2CF3Coal_1)] <-0
AR_F_CF2CF3Coal_Total <- AR_F_CF2CF3Coal_0 + AR_F_CF2CF3Coal_1
AR_F_CF2CF3Coal_Affected <- 
  round((AR_F_CF2CF3Coal_1/AR_F_CF2CF3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_CF2CF3Coal_row <- c("CF2CF3Coal", AR_CF2CF3Coal_Total, AR_CF2CF3Coal_1, 
                       AR_CF2CF3Coal_Affected, AR_M_CF2CF3Coal_Total, 
                       AR_M_CF2CF3Coal_1, AR_M_CF2CF3Coal_Affected, 
                       AR_F_CF2CF3Coal_Total, AR_F_CF2CF3Coal_1, 
                       AR_F_CF2CF3Coal_Affected)

## CF3MT3Coal ##
# Overall N and frequency #
AR_CF3MT3Coal <- AR_wo_probables %>% count(CF3MT3Coal)
AR_CF3MT3Coal <- 
  AR_CF3MT3Coal[(AR_CF3MT3Coal$CF3MT3Coal != 6 & AR_CF3MT3Coal$CF3MT3Coal != 9),]
AR_CF3MT3Coal_0 <- AR_CF3MT3Coal[1,2]
AR_CF3MT3Coal_0[is.na(AR_CF3MT3Coal_0)] <-0
AR_CF3MT3Coal_1 <- AR_CF3MT3Coal[2,2]
AR_CF3MT3Coal_1[is.na(AR_CF3MT3Coal_1)] <-0
AR_CF3MT3Coal_Total <- AR_CF3MT3Coal_0 + AR_CF3MT3Coal_1
AR_CF3MT3Coal_Affected <- round((AR_CF3MT3Coal_1/AR_CF3MT3Coal_Total)*100, 2)

# Male N and frequency
AR_M_CF3MT3Coal <- AR_M %>% count(CF3MT3Coal)
AR_M_CF3MT3Coal <- 
  AR_M_CF3MT3Coal[(AR_M_CF3MT3Coal$CF3MT3Coal != 6 & AR_M_CF3MT3Coal$CF3MT3Coal != 9),]
AR_M_CF3MT3Coal_0 <- AR_M_CF3MT3Coal[1,2]
AR_M_CF3MT3Coal_0[is.na(AR_M_CF3MT3Coal_0)] <-0
AR_M_CF3MT3Coal_1 <- AR_M_CF3MT3Coal[2,2]
AR_M_CF3MT3Coal_1[is.na(AR_M_CF3MT3Coal_1)] <-0
AR_M_CF3MT3Coal_Total <- AR_M_CF3MT3Coal_0 + AR_M_CF3MT3Coal_1
AR_M_CF3MT3Coal_Affected <- 
  round((AR_M_CF3MT3Coal_1/AR_M_CF3MT3Coal_Total)*100, 2)

# Female N and frequency
AR_F_CF3MT3Coal <- AR_F %>% count(CF3MT3Coal)
AR_F_CF3MT3Coal <- 
  AR_F_CF3MT3Coal[(AR_F_CF3MT3Coal$CF3MT3Coal != 6 & AR_F_CF3MT3Coal$CF3MT3Coal != 9),]
AR_F_CF3MT3Coal_0 <- AR_F_CF3MT3Coal[1,2]
AR_F_CF3MT3Coal_0[is.na(AR_F_CF3MT3Coal_0)] <-0
AR_F_CF3MT3Coal_1 <- AR_F_CF3MT3Coal[2,2]
AR_F_CF3MT3Coal_1[is.na(AR_F_CF3MT3Coal_1)] <-0
AR_F_CF3MT3Coal_Total <- AR_F_CF3MT3Coal_0 + AR_F_CF3MT3Coal_1
AR_F_CF3MT3Coal_Affected <- 
  round((AR_F_CF3MT3Coal_1/AR_F_CF3MT3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_CF3MT3Coal_row <- c("CF3MT3Coal", AR_CF3MT3Coal_Total, AR_CF3MT3Coal_1, 
                       AR_CF3MT3Coal_Affected, AR_M_CF3MT3Coal_Total, 
                       AR_M_CF3MT3Coal_1, AR_M_CF3MT3Coal_Affected, 
                       AR_F_CF3MT3Coal_Total, AR_F_CF3MT3Coal_1, 
                       AR_F_CF3MT3Coal_Affected)
## OsIntermet ##
# Overall N and frequency #
AR_OsIntermet <- AR_wo_probables %>% count(OsIntermet)
AR_OsIntermet <- 
  AR_OsIntermet[(AR_OsIntermet$OsIntermet != 6 & AR_OsIntermet$OsIntermet != 9),]
AR_OsIntermet_0 <- AR_OsIntermet[1,2]
AR_OsIntermet_0[is.na(AR_OsIntermet_0)] <-0
AR_OsIntermet_1 <- AR_OsIntermet[2,2]
AR_OsIntermet_1[is.na(AR_OsIntermet_1)] <-0
AR_OsIntermet_Total <- AR_OsIntermet_0 + AR_OsIntermet_1
AR_OsIntermet_Affected <- round((AR_OsIntermet_1/AR_OsIntermet_Total)*100, 2)

# Male N and frequency
AR_M_OsIntermet <- AR_M %>% count(OsIntermet)
AR_M_OsIntermet <- 
  AR_M_OsIntermet[(AR_M_OsIntermet$OsIntermet != 6 & AR_M_OsIntermet$OsIntermet != 9),]
AR_M_OsIntermet_0 <- AR_M_OsIntermet[1,2]
AR_M_OsIntermet_0[is.na(AR_M_OsIntermet_0)] <-0
AR_M_OsIntermet_1 <- AR_M_OsIntermet[2,2]
AR_M_OsIntermet_1[is.na(AR_M_OsIntermet_1)] <-0
AR_M_OsIntermet_Total <- AR_M_OsIntermet_0 + AR_M_OsIntermet_1
AR_M_OsIntermet_Affected <- 
  round((AR_M_OsIntermet_1/AR_M_OsIntermet_Total)*100, 2)

# Female N and frequency
AR_F_OsIntermet <- AR_F %>% count(OsIntermet)
AR_F_OsIntermet <- 
  AR_F_OsIntermet[(AR_F_OsIntermet$OsIntermet != 6 & AR_F_OsIntermet$OsIntermet != 9),]
AR_F_OsIntermet_0 <- AR_F_OsIntermet[1,2]
AR_F_OsIntermet_0[is.na(AR_F_OsIntermet_0)] <-0
AR_F_OsIntermet_1 <- AR_F_OsIntermet[2,2]
AR_F_OsIntermet_1[is.na(AR_F_OsIntermet_1)] <-0
AR_F_OsIntermet_Total <- AR_F_OsIntermet_0 + AR_F_OsIntermet_1
AR_F_OsIntermet_Affected <- 
  round((AR_F_OsIntermet_1/AR_F_OsIntermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_OsIntermet_row <- c("OsIntermet", AR_OsIntermet_Total, AR_OsIntermet_1, 
                       AR_OsIntermet_Affected, AR_M_OsIntermet_Total, 
                       AR_M_OsIntermet_1, AR_M_OsIntermet_Affected, 
                       AR_F_OsIntermet_Total, AR_F_OsIntermet_1, 
                       AR_F_OsIntermet_Affected)
# OsIntermet unclustered: CF1Intermet #
# Overall N and frequency #
AR_CF1Intermet <- AR_wo_probables %>% count(CF1Intermet)
AR_CF1Intermet <- 
  AR_CF1Intermet[(AR_CF1Intermet$CF1Intermet != 6 & AR_CF1Intermet$CF1Intermet != 9),]
AR_CF1Intermet_0 <- AR_CF1Intermet[1,2]
AR_CF1Intermet_0[is.na(AR_CF1Intermet_0)] <-0
AR_CF1Intermet_1 <- AR_CF1Intermet[2,2]
AR_CF1Intermet_1[is.na(AR_CF1Intermet_1)] <-0
AR_CF1Intermet_Total <- na.omit(AR_CF1Intermet_0) + na.omit(AR_CF1Intermet_1)
AR_CF1Intermet_Affected <- round((AR_CF1Intermet_1/AR_CF1Intermet_Total)*100, 2)

# Male N and frequency
AR_M_CF1Intermet <- AR_M %>% count(CF1Intermet)
AR_M_CF1Intermet <- 
  AR_M_CF1Intermet[(AR_M_CF1Intermet$CF1Intermet != 6 & AR_M_CF1Intermet$CF1Intermet != 9),]
AR_M_CF1Intermet_0 <- AR_M_CF1Intermet[1,2]
AR_M_CF1Intermet_0[is.na(AR_M_CF1Intermet_0)] <-0
AR_M_CF1Intermet_1 <- AR_M_CF1Intermet[2,2]
AR_M_CF1Intermet_1[is.na(AR_M_CF1Intermet_1)] <-0
AR_M_CF1Intermet_Total <- AR_M_CF1Intermet_0 + AR_M_CF1Intermet_1
AR_M_CF1Intermet_Affected <- 
  round((AR_M_CF1Intermet_1/AR_M_CF1Intermet_Total)*100, 2)

# Female N and frequency
AR_F_CF1Intermet <- AR_F %>% count(CF1Intermet)
AR_F_CF1Intermet <- 
  AR_F_CF1Intermet[(AR_F_CF1Intermet$CF1Intermet != 6 & AR_F_CF1Intermet$CF1Intermet != 9),]
AR_F_CF1Intermet_0 <- AR_F_CF1Intermet[1,2]
AR_F_CF1Intermet_0[is.na(AR_F_CF1Intermet_0)] <-0
AR_F_CF1Intermet_1 <- AR_F_CF1Intermet[2,2]
AR_F_CF1Intermet_1[is.na(AR_F_CF1Intermet_1)] <-0
AR_F_CF1Intermet_Total <- AR_F_CF1Intermet_0 + AR_F_CF1Intermet_1
AR_F_CF1Intermet_Affected <- 
  round((AR_F_CF1Intermet_1/AR_F_CF1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_CF1Intermet_row <- c("CF1Intermet", AR_CF1Intermet_Total, AR_CF1Intermet_1, 
                        AR_CF1Intermet_Affected, AR_M_CF1Intermet_Total, 
                        AR_M_CF1Intermet_1, AR_M_CF1Intermet_Affected, 
                        AR_F_CF1Intermet_Total, AR_F_CF1Intermet_1, 
                        AR_F_CF1Intermet_Affected)

# OsIntermet unclustered: MT1Intermet #
# Overall N and frequency #
AR_MT1Intermet <- AR_wo_probables %>% count(MT1Intermet)
AR_MT1Intermet <- 
  AR_MT1Intermet[(AR_MT1Intermet$MT1Intermet != 6 & AR_MT1Intermet$MT1Intermet != 9),]
AR_MT1Intermet_0 <- AR_MT1Intermet[1,2]
AR_MT1Intermet_0[is.na(AR_MT1Intermet_0)] <-0
AR_MT1Intermet_1 <- AR_MT1Intermet[2,2]
AR_MT1Intermet_1[is.na(AR_MT1Intermet_1)] <-0
AR_MT1Intermet_Total <- na.omit(AR_MT1Intermet_0) + na.omit(AR_MT1Intermet_1)
AR_MT1Intermet_Affected <- round((AR_MT1Intermet_1/AR_MT1Intermet_Total)*100, 2)

# Male N and frequency
AR_M_MT1Intermet <- AR_M %>% count(MT1Intermet)
AR_M_MT1Intermet <- 
  AR_M_MT1Intermet[(AR_M_MT1Intermet$MT1Intermet != 6 & AR_M_MT1Intermet$MT1Intermet != 9),]
AR_M_MT1Intermet_0 <- AR_M_MT1Intermet[1,2]
AR_M_MT1Intermet_0[is.na(AR_M_MT1Intermet_0)] <-0
AR_M_MT1Intermet_1 <- AR_M_MT1Intermet[2,2]
AR_M_MT1Intermet_1[is.na(AR_M_MT1Intermet_1)] <-0
AR_M_MT1Intermet_Total <- AR_M_MT1Intermet_0 + AR_M_MT1Intermet_1
AR_M_MT1Intermet_Affected <- 
  round((AR_M_MT1Intermet_1/AR_M_MT1Intermet_Total)*100, 2)

# Female N and frequency
AR_F_MT1Intermet <- AR_F %>% count(MT1Intermet)
AR_F_MT1Intermet <- 
  AR_F_MT1Intermet[(AR_F_MT1Intermet$MT1Intermet != 6 & AR_F_MT1Intermet$MT1Intermet != 9),]
AR_F_MT1Intermet_0 <- AR_F_MT1Intermet[1,2]
AR_F_MT1Intermet_0[is.na(AR_F_MT1Intermet_0)] <-0
AR_F_MT1Intermet_1 <- AR_F_MT1Intermet[2,2]
AR_F_MT1Intermet_1[is.na(AR_F_MT1Intermet_1)] <-0
AR_F_MT1Intermet_Total <- AR_F_MT1Intermet_0 + AR_F_MT1Intermet_1
AR_F_MT1Intermet_Affected <- 
  round((AR_F_MT1Intermet_1/AR_F_MT1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_MT1Intermet_row <- c("MT1Intermet", AR_MT1Intermet_Total, AR_MT1Intermet_1, 
                        AR_MT1Intermet_Affected, AR_M_MT1Intermet_Total, 
                        AR_M_MT1Intermet_1, AR_M_MT1Intermet_Affected, 
                        AR_F_MT1Intermet_Total, AR_F_MT1Intermet_1, 
                        AR_F_MT1Intermet_Affected)

# OsIntermet unclustered: MT2Intermet #
# Overall N and frequency #
AR_MT2Intermet <- AR_wo_probables %>% count(MT2Intermet)
AR_MT2Intermet <- 
  AR_MT2Intermet[(AR_MT2Intermet$MT2Intermet != 6 & AR_MT2Intermet$MT2Intermet != 9),]
AR_MT2Intermet_0 <- AR_MT2Intermet[1,2]
AR_MT2Intermet_0[is.na(AR_MT2Intermet_0)] <-0
AR_MT2Intermet_1 <- AR_MT2Intermet[2,2]
AR_MT2Intermet_1[is.na(AR_MT2Intermet_1)] <-0
AR_MT2Intermet_Total <- na.omit(AR_MT2Intermet_0) + na.omit(AR_MT2Intermet_1)
AR_MT2Intermet_Affected <- round((AR_MT2Intermet_1/AR_MT2Intermet_Total)*100, 2)

# Male N and frequency
AR_M_MT2Intermet <- AR_M %>% count(MT2Intermet)
AR_M_MT2Intermet <- 
  AR_M_MT2Intermet[(AR_M_MT2Intermet$MT2Intermet != 6 & AR_M_MT2Intermet$MT2Intermet != 9),]
AR_M_MT2Intermet_0 <- AR_M_MT2Intermet[1,2]
AR_M_MT2Intermet_0[is.na(AR_M_MT2Intermet_0)] <-0
AR_M_MT2Intermet_1 <- AR_M_MT2Intermet[2,2]
AR_M_MT2Intermet_1[is.na(AR_M_MT2Intermet_1)] <-0
AR_M_MT2Intermet_Total <- AR_M_MT2Intermet_0 + AR_M_MT2Intermet_1
AR_M_MT2Intermet_Affected <- 
  round((AR_M_MT2Intermet_1/AR_M_MT2Intermet_Total)*100, 2)

# Female N and frequency
AR_F_MT2Intermet <- AR_F %>% count(MT2Intermet)
AR_F_MT2Intermet <- 
  AR_F_MT2Intermet[(AR_F_MT2Intermet$MT2Intermet != 6 & AR_F_MT2Intermet$MT2Intermet != 9),]
AR_F_MT2Intermet_0 <- AR_F_MT2Intermet[1,2]
AR_F_MT2Intermet_0[is.na(AR_F_MT2Intermet_0)] <-0
AR_F_MT2Intermet_1 <- AR_F_MT2Intermet[2,2]
AR_F_MT2Intermet_1[is.na(AR_F_MT2Intermet_1)] <-0
AR_F_MT2Intermet_Total <- AR_F_MT2Intermet_0 + AR_F_MT2Intermet_1
AR_F_MT2Intermet_Affected <- 
  round((AR_F_MT2Intermet_1/AR_F_MT2Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_MT2Intermet_row <- c("MT2Intermet", AR_MT2Intermet_Total, AR_MT2Intermet_1, 
                        AR_MT2Intermet_Affected, AR_M_MT2Intermet_Total, 
                        AR_M_MT2Intermet_1, AR_M_MT2Intermet_Affected, 
                        AR_F_MT2Intermet_Total, AR_F_MT2Intermet_1, 
                        AR_F_MT2Intermet_Affected)


## AnkleCoal ##
# Overall N and frequency #
AR_AnkleCoal <- AR_wo_probables %>% count(AnkleCoal)
AR_AnkleCoal <- 
  AR_AnkleCoal[(AR_AnkleCoal$AnkleCoal != 6 & AR_AnkleCoal$AnkleCoal != 9),]
AR_AnkleCoal_0 <- AR_AnkleCoal[1,2]; AR_AnkleCoal_0[is.na(AR_AnkleCoal_0)] <-0
AR_AnkleCoal_1 <- AR_AnkleCoal[2,2]; AR_AnkleCoal_1[is.na(AR_AnkleCoal_1)] <-0
AR_AnkleCoal_Total <- na.omit(AR_AnkleCoal_0) + na.omit(AR_AnkleCoal_1)
AR_AnkleCoal_Affected <- round((AR_AnkleCoal_1/AR_AnkleCoal_Total)*100, 2)

# Male N and frequency
AR_M_AnkleCoal <- AR_M %>% count(AnkleCoal)
AR_M_AnkleCoal <- 
  AR_M_AnkleCoal[(AR_M_AnkleCoal$AnkleCoal != 6 & AR_M_AnkleCoal$AnkleCoal != 9),]
AR_M_AnkleCoal_0 <- AR_M_AnkleCoal[1,2]
AR_M_AnkleCoal_0[is.na(AR_M_AnkleCoal_0)] <-0
AR_M_AnkleCoal_1 <- AR_M_AnkleCoal[2,2]
AR_M_AnkleCoal_1[is.na(AR_M_AnkleCoal_1)] <-0
AR_M_AnkleCoal_Total <- AR_M_AnkleCoal_0 + AR_M_AnkleCoal_1
AR_M_AnkleCoal_Affected <- round((AR_M_AnkleCoal_1/AR_M_AnkleCoal_Total)*100, 2)

# Female N and frequency
AR_F_AnkleCoal <- AR_F %>% count(AnkleCoal)
AR_F_AnkleCoal <- 
  AR_F_AnkleCoal[(AR_F_AnkleCoal$AnkleCoal != 6 & AR_F_AnkleCoal$AnkleCoal != 9),]
AR_F_AnkleCoal_0 <- AR_F_AnkleCoal[1,2]
AR_F_AnkleCoal_0[is.na(AR_F_AnkleCoal_0)] <-0
AR_F_AnkleCoal_1 <- AR_F_AnkleCoal[2,2]
AR_F_AnkleCoal_1[is.na(AR_F_AnkleCoal_1)] <-0
AR_F_AnkleCoal_Total <- AR_F_AnkleCoal_0 + AR_F_AnkleCoal_1
AR_F_AnkleCoal_Affected <- round((AR_F_AnkleCoal_1/AR_F_AnkleCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_AnkleCoal_row <- c("AnkleCoal", AR_AnkleCoal_Total, AR_AnkleCoal_1,
                      AR_AnkleCoal_Affected, AR_M_AnkleCoal_Total, 
                      AR_M_AnkleCoal_1, AR_M_AnkleCoal_Affected, AR_F_AnkleCoal_Total, 
                      AR_F_AnkleCoal_1, AR_F_AnkleCoal_Affected)

# AnkleCoal unclustered: CalcNavCoal #
# Overall N and frequency #
AR_CalcNavCoal <- AR_wo_probables %>% count(CalcNavCoal)
AR_CalcNavCoal <- 
  AR_CalcNavCoal[(AR_CalcNavCoal$CalcNavCoal != 6 & AR_CalcNavCoal$CalcNavCoal != 9),]
AR_CalcNavCoal_0 <- AR_CalcNavCoal[1,2]
AR_CalcNavCoal_0[is.na(AR_CalcNavCoal_0)] <-0
AR_CalcNavCoal_1 <- AR_CalcNavCoal[2,2]
AR_CalcNavCoal_1[is.na(AR_CalcNavCoal_1)] <-0
AR_CalcNavCoal_Total <- na.omit(AR_CalcNavCoal_0) + na.omit(AR_CalcNavCoal_1)
AR_CalcNavCoal_Affected <- round((AR_CalcNavCoal_1/AR_CalcNavCoal_Total)*100, 2)

# Male N and frequency
AR_M_CalcNavCoal <- AR_M %>% count(CalcNavCoal)
AR_M_CalcNavCoal <- 
  AR_M_CalcNavCoal[(AR_M_CalcNavCoal$CalcNavCoal != 6 & AR_M_CalcNavCoal$CalcNavCoal != 9),]
AR_M_CalcNavCoal_0 <- AR_M_CalcNavCoal[1,2]
AR_M_CalcNavCoal_0[is.na(AR_M_CalcNavCoal_0)] <-0
AR_M_CalcNavCoal_1 <- AR_M_CalcNavCoal[2,2]
AR_M_CalcNavCoal_1[is.na(AR_M_CalcNavCoal_1)] <-0
AR_M_CalcNavCoal_Total <- AR_M_CalcNavCoal_0 + AR_M_CalcNavCoal_1
AR_M_CalcNavCoal_Affected <- 
  round((AR_M_CalcNavCoal_1/AR_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
AR_F_CalcNavCoal <- AR_F %>% count(CalcNavCoal)
AR_F_CalcNavCoal <- 
  AR_F_CalcNavCoal[(AR_F_CalcNavCoal$CalcNavCoal != 6 & AR_F_CalcNavCoal$CalcNavCoal != 9),]
AR_F_CalcNavCoal_0 <- AR_F_CalcNavCoal[1,2]
AR_F_CalcNavCoal_0[is.na(AR_F_CalcNavCoal_0)] <-0
AR_F_CalcNavCoal_1 <- AR_F_CalcNavCoal[2,2]
AR_F_CalcNavCoal_1[is.na(AR_F_CalcNavCoal_1)] <-0
AR_F_CalcNavCoal_Total <- AR_F_CalcNavCoal_0 + AR_F_CalcNavCoal_1
AR_F_CalcNavCoal_Affected <- 
  round((AR_F_CalcNavCoal_1/AR_F_CalcNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_CalcNavCoal_row <- c("CalcNavCoal", AR_CalcNavCoal_Total, AR_CalcNavCoal_1, 
                        AR_CalcNavCoal_Affected, AR_M_CalcNavCoal_Total, 
                        AR_M_CalcNavCoal_1, AR_M_CalcNavCoal_Affected, 
                        AR_F_CalcNavCoal_Total, AR_F_CalcNavCoal_1, 
                        AR_F_CalcNavCoal_Affected)

# AnkleCoal unclustered: TaloCalcCoal #
# Overall N and frequency #
AR_TaloCalcCoal <- AR_wo_probables %>% count(TaloCalcCoal)
AR_TaloCalcCoal <- 
  AR_TaloCalcCoal[(AR_TaloCalcCoal$TaloCalcCoal != 6 & AR_TaloCalcCoal$TaloCalcCoal != 9),]
AR_TaloCalcCoal_0 <- AR_TaloCalcCoal[1,2]
AR_TaloCalcCoal_0[is.na(AR_TaloCalcCoal_0)] <-0
AR_TaloCalcCoal_1 <- AR_TaloCalcCoal[2,2]
AR_TaloCalcCoal_1[is.na(AR_TaloCalcCoal_1)] <-0
AR_TaloCalcCoal_Total <- na.omit(AR_TaloCalcCoal_0) + na.omit(AR_TaloCalcCoal_1)
AR_TaloCalcCoal_Affected <- 
  round((AR_TaloCalcCoal_1/AR_TaloCalcCoal_Total)*100, 2)

# Male N and frequency
AR_M_TaloCalcCoal <- AR_M %>% count(TaloCalcCoal)
AR_M_TaloCalcCoal <- 
  AR_M_TaloCalcCoal[(AR_M_TaloCalcCoal$TaloCalcCoal != 6 & AR_M_TaloCalcCoal$TaloCalcCoal != 9),]
AR_M_TaloCalcCoal_0 <- AR_M_TaloCalcCoal[1,2]
AR_M_TaloCalcCoal_0[is.na(AR_M_TaloCalcCoal_0)] <-0
AR_M_TaloCalcCoal_1 <- AR_M_TaloCalcCoal[2,2]
AR_M_TaloCalcCoal_1[is.na(AR_M_TaloCalcCoal_1)] <-0
AR_M_TaloCalcCoal_Total <- AR_M_TaloCalcCoal_0 + AR_M_TaloCalcCoal_1
AR_M_TaloCalcCoal_Affected <- 
  round((AR_M_TaloCalcCoal_1/AR_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
AR_F_TaloCalcCoal <- AR_F %>% count(TaloCalcCoal)
AR_F_TaloCalcCoal <- 
  AR_F_TaloCalcCoal[(AR_F_TaloCalcCoal$TaloCalcCoal != 6 & AR_F_TaloCalcCoal$TaloCalcCoal != 9),]
AR_F_TaloCalcCoal_0 <- AR_F_TaloCalcCoal[1,2]
AR_F_TaloCalcCoal_0[is.na(AR_F_TaloCalcCoal_0)] <-0
AR_F_TaloCalcCoal_1 <- AR_F_TaloCalcCoal[2,2]
AR_F_TaloCalcCoal_1[is.na(AR_F_TaloCalcCoal_1)] <-0
AR_F_TaloCalcCoal_Total <- AR_F_TaloCalcCoal_0 + AR_F_TaloCalcCoal_1
AR_F_TaloCalcCoal_Affected <- 
  round((AR_F_TaloCalcCoal_1/AR_F_TaloCalcCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
AR_TaloCalcCoal_row <- c("TaloCalcCoal", AR_TaloCalcCoal_Total, AR_TaloCalcCoal_1, 
                         AR_TaloCalcCoal_Affected, AR_M_TaloCalcCoal_Total, 
                         AR_M_TaloCalcCoal_1, AR_M_TaloCalcCoal_Affected, 
                         AR_F_TaloCalcCoal_Total, AR_F_TaloCalcCoal_1, 
                         AR_F_TaloCalcCoal_Affected)


## Table with overview of trait frequencies for AR ##
AR_freq <- data.frame(rbind(AR_AccessNav_row, AR_BrachyD_row, AR_BrachyMT1_row,
                            AR_BrachyMT4_row, AR_BrachyPP1_row, 
                            AR_CalcCubCoal_row, AR_TaloNavCoal_row, 
                            AR_CF2CF3Coal_row, AR_CF3MT3Coal_row, 
                            AR_OsIntermet_row, AR_CF1Intermet_row,  
                            AR_MT1Intermet_row, AR_MT2Intermet_row, 
                            AR_AnkleCoal_row, AR_CalcNavCoal_row, 
                            AR_TaloCalcCoal_row))
names(AR_freq) <- c("Anomaly", "Overall N (AR)", "Overall affected (AR)", 
                    "Overall % (AR)", "Male N (AR)", "Male affected (AR)",
                    "Male % (AR)", "Female N (AR)", "Female affected (AR)", 
                    "Female % (AR)")

### the Eindhoven collection ###


## AccessNav ##
# Overall N and frequency #
EH_AccessNav <- EH_wo_probables %>% count(AccessNav)
EH_AccessNav <- 
  EH_AccessNav[(EH_AccessNav$AccessNav != 6 & EH_AccessNav$AccessNav != 9),]
EH_AccessNav_0 <- EH_AccessNav[1,2]
EH_AccessNav_0[is.na(EH_AccessNav_0)] <-0
EH_AccessNav_1 <- EH_AccessNav[2,2]
EH_AccessNav_1[is.na(EH_AccessNav_1)] <-0
EH_AccessNav_Total <- EH_AccessNav_0 + EH_AccessNav_1
EH_AccessNav_Affected <- round((EH_AccessNav_1/EH_AccessNav_Total)*100, 2)

# Male N and frequency
EH_M_AccessNav <- EH_M %>% count(AccessNav)
EH_M_AccessNav <- 
  EH_M_AccessNav[(EH_M_AccessNav$AccessNav != 6 & EH_M_AccessNav$AccessNav != 9),]
EH_M_AccessNav_0 <- EH_M_AccessNav[1,2]
EH_M_AccessNav_0[is.na(EH_M_AccessNav_0)] <-0
EH_M_AccessNav_1 <- EH_M_AccessNav[2,2]
EH_M_AccessNav_1[is.na(EH_M_AccessNav_1)] <-0
EH_M_AccessNav_Total <- EH_M_AccessNav_0 + EH_M_AccessNav_1
EH_M_AccessNav_Affected <- round((EH_M_AccessNav_1/EH_M_AccessNav_Total)*100, 2)

# Female N and frequency
EH_F_AccessNav <- EH_F %>% count(AccessNav)
EH_F_AccessNav <- 
  EH_F_AccessNav[(EH_F_AccessNav$AccessNav != 6 & EH_F_AccessNav$AccessNav != 9),]
EH_F_AccessNav_0 <- EH_F_AccessNav[1,2]
EH_F_AccessNav_0[is.na(EH_F_AccessNav_0)] <-0
EH_F_AccessNav_1 <- EH_F_AccessNav[2,2]
EH_F_AccessNav_1[is.na(EH_F_AccessNav_1)] <-0
EH_F_AccessNav_Total <- EH_F_AccessNav_0 + EH_F_AccessNav_1
EH_F_AccessNav_Affected <- round((EH_F_AccessNav_1/EH_F_AccessNav_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_AccessNav_row <- c("AccessNav", EH_AccessNav_Total, EH_AccessNav_1, 
                      EH_AccessNav_Affected, EH_M_AccessNav_Total, 
                      EH_M_AccessNav_1, EH_M_AccessNav_Affected, 
                      EH_F_AccessNav_Total, EH_F_AccessNav_1, 
                      EH_F_AccessNav_Affected)


## BrachyD ##
# Overall N and frequency #
EH_BrachyD <- EH_wo_probables %>% count(BrachyD)
EH_BrachyD <- EH_BrachyD[(EH_BrachyD$BrachyD != 6 & EH_BrachyD$BrachyD != 9),]
EH_BrachyD_0 <- EH_BrachyD[1,2]; EH_BrachyD_0[is.na(EH_BrachyD_0)] <-0
EH_BrachyD_1 <- EH_BrachyD[2,2]; EH_BrachyD_1[is.na(EH_BrachyD_1)] <-0
EH_BrachyD_Total <- EH_BrachyD_0 + EH_BrachyD_1
EH_BrachyD_Affected <- round((EH_BrachyD_1/EH_BrachyD_Total)*100, 2)

# Male N and frequency
EH_M_BrachyD <- EH_M %>% count(BrachyD)
EH_M_BrachyD <- 
  EH_M_BrachyD[(EH_M_BrachyD$BrachyD != 6 & EH_M_BrachyD$BrachyD != 9),]
EH_M_BrachyD_0 <- EH_M_BrachyD[1,2]; EH_M_BrachyD_0[is.na(EH_M_BrachyD_0)] <-0
EH_M_BrachyD_1 <- EH_M_BrachyD[2,2]; EH_M_BrachyD_1[is.na(EH_M_BrachyD_1)] <-0
EH_M_BrachyD_Total <- EH_M_BrachyD_0 + EH_M_BrachyD_1
EH_M_BrachyD_Affected <- round((EH_M_BrachyD_1/EH_M_BrachyD_Total)*100, 2)

# Female N and frequency
EH_F_BrachyD <- EH_F %>% count(BrachyD)
EH_F_BrachyD <- 
  EH_F_BrachyD[(EH_F_BrachyD$BrachyD != 6 & EH_F_BrachyD$BrachyD != 9),]
EH_F_BrachyD_0 <- EH_F_BrachyD[1,2]; EH_F_BrachyD_0[is.na(EH_F_BrachyD_0)] <-0
EH_F_BrachyD_1 <- EH_F_BrachyD[2,2]; EH_F_BrachyD_1[is.na(EH_F_BrachyD_1)] <-0
EH_F_BrachyD_Total <- EH_F_BrachyD_0 + EH_F_BrachyD_1
EH_F_BrachyD_Affected <- round((EH_F_BrachyD_1/EH_F_BrachyD_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_BrachyD_row <- c("BrachyD", EH_BrachyD_Total, EH_BrachyD_1, 
                    EH_BrachyD_Affected, EH_M_BrachyD_Total, 
                    EH_M_BrachyD_1, EH_M_BrachyD_Affected, 
                    EH_F_BrachyD_Total, EH_F_BrachyD_1, 
                    EH_F_BrachyD_Affected)


## BrachyMT1 ##
# Overall N and frequency #
EH_BrachyMT1 <- EH_wo_probables %>% count(BrachyMT1)
EH_BrachyMT1 <- 
  EH_BrachyMT1[(EH_BrachyMT1$BrachyMT1 != 6 & EH_BrachyMT1$BrachyMT1 != 9),]
EH_BrachyMT1_0 <- EH_BrachyMT1[1,2]; EH_BrachyMT1_0[is.na(EH_BrachyMT1_0)] <-0
EH_BrachyMT1_1 <- EH_BrachyMT1[2,2]; EH_BrachyMT1_1[is.na(EH_BrachyMT1_1)] <-0
EH_BrachyMT1_Total <- EH_BrachyMT1_0 + EH_BrachyMT1_1
EH_BrachyMT1_Affected <- round((EH_BrachyMT1_1/EH_BrachyMT1_Total)*100, 2)

# Male N and frequency
EH_M_BrachyMT1 <- EH_M %>% count(BrachyMT1)
EH_M_BrachyMT1 <- 
  EH_M_BrachyMT1[(EH_M_BrachyMT1$BrachyMT1 != 6 & EH_M_BrachyMT1$BrachyMT1 != 9),]
EH_M_BrachyMT1_0 <- EH_M_BrachyMT1[1,2]
EH_M_BrachyMT1_0[is.na(EH_M_BrachyMT1_0)] <-0
EH_M_BrachyMT1_1 <- EH_M_BrachyMT1[2,2]
EH_M_BrachyMT1_1[is.na(EH_M_BrachyMT1_1)] <-0
EH_M_BrachyMT1_Total <- EH_M_BrachyMT1_0 + EH_M_BrachyMT1_1
EH_M_BrachyMT1_Affected <- round((EH_M_BrachyMT1_1/EH_M_BrachyMT1_Total)*100, 2)

# Female N and frequency
EH_F_BrachyMT1 <- EH_F %>% count(BrachyMT1)
EH_F_BrachyMT1 <- 
  EH_F_BrachyMT1[(EH_F_BrachyMT1$BrachyMT1 != 6 & EH_F_BrachyMT1$BrachyMT1 != 9),]
EH_F_BrachyMT1_0 <- EH_F_BrachyMT1[1,2]
EH_F_BrachyMT1_0[is.na(EH_F_BrachyMT1_0)] <-0
EH_F_BrachyMT1_1 <- EH_F_BrachyMT1[2,2]
EH_F_BrachyMT1_1[is.na(EH_F_BrachyMT1_1)] <-0
EH_F_BrachyMT1_Total <- EH_F_BrachyMT1_0 + EH_F_BrachyMT1_1
EH_F_BrachyMT1_Affected <- round((EH_F_BrachyMT1_1/EH_F_BrachyMT1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_BrachyMT1_row <- c("BrachyMT1", EH_BrachyMT1_Total, EH_BrachyMT1_1, 
                      EH_BrachyMT1_Affected, EH_M_BrachyMT1_Total, 
                      EH_M_BrachyMT1_1, EH_M_BrachyMT1_Affected, 
                      EH_F_BrachyMT1_Total, EH_F_BrachyMT1_1, 
                      EH_F_BrachyMT1_Affected)
## BrachyMT4 ##
# Overall N and frequency #
EH_BrachyMT4 <- EH_wo_probables %>% count(BrachyMT4)
EH_BrachyMT4 <- 
  EH_BrachyMT4[(EH_BrachyMT4$BrachyMT4 != 6 & EH_BrachyMT4$BrachyMT4 != 9),]
EH_BrachyMT4_0 <- EH_BrachyMT4[1,2]
EH_BrachyMT4_0[is.na(EH_BrachyMT4_0)] <-0
EH_BrachyMT4_1 <- EH_BrachyMT4[2,2]
EH_BrachyMT4_1[is.na(EH_BrachyMT4_1)] <-0
EH_BrachyMT4_Total <- EH_BrachyMT4_0 + EH_BrachyMT4_1
EH_BrachyMT4_Affected <- round((EH_BrachyMT4_1/EH_BrachyMT4_Total)*100, 2)

# Male N and frequency
EH_M_BrachyMT4 <- EH_M %>% count(BrachyMT4)
EH_M_BrachyMT4 <- 
  EH_M_BrachyMT4[(EH_M_BrachyMT4$BrachyMT4 != 6 & EH_M_BrachyMT4$BrachyMT4 != 9),]
EH_M_BrachyMT4_0 <- EH_M_BrachyMT4[1,2]
EH_M_BrachyMT4_0[is.na(EH_M_BrachyMT4_0)] <-0
EH_M_BrachyMT4_1 <- EH_M_BrachyMT4[2,2]
EH_M_BrachyMT4_1[is.na(EH_M_BrachyMT4_1)] <-0
EH_M_BrachyMT4_Total <- EH_M_BrachyMT4_0 + EH_M_BrachyMT4_1
EH_M_BrachyMT4_Affected <- round((EH_M_BrachyMT4_1/EH_M_BrachyMT4_Total)*100, 2)

# Female N and frequency
EH_F_BrachyMT4 <- EH_F %>% count(BrachyMT4)
EH_F_BrachyMT4 <- 
  EH_F_BrachyMT4[(EH_F_BrachyMT4$BrachyMT4 != 6 & EH_F_BrachyMT4$BrachyMT4 != 9),]
EH_F_BrachyMT4_0 <- EH_F_BrachyMT4[1,2]
EH_F_BrachyMT4_0[is.na(EH_F_BrachyMT4_0)] <-0
EH_F_BrachyMT4_1 <- EH_F_BrachyMT4[2,2]
EH_F_BrachyMT4_1[is.na(EH_F_BrachyMT4_1)] <-0
EH_F_BrachyMT4_Total <- EH_F_BrachyMT4_0 + EH_F_BrachyMT4_1
EH_F_BrachyMT4_Affected <- round((EH_F_BrachyMT4_1/EH_F_BrachyMT4_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_BrachyMT4_row <- c("BrachyMT4", EH_BrachyMT4_Total, EH_BrachyMT4_1, 
                      EH_BrachyMT4_Affected, EH_M_BrachyMT4_Total, 
                      EH_M_BrachyMT4_1, EH_M_BrachyMT4_Affected, 
                      EH_F_BrachyMT4_Total, EH_F_BrachyMT4_1, 
                      EH_F_BrachyMT4_Affected)

## BrachyPP1 ##
# Overall N and frequency #
EH_BrachyPP1 <- EH_wo_probables %>% count(BrachyPP1)
EH_BrachyPP1 <- 
  EH_BrachyPP1[(EH_BrachyPP1$BrachyPP1 != 6 & EH_BrachyPP1$BrachyPP1 != 9),]
EH_BrachyPP1_0 <- EH_BrachyPP1[1,2]; EH_BrachyPP1_0[is.na(EH_BrachyPP1_0)] <-0
EH_BrachyPP1_1 <- EH_BrachyPP1[2,2]; EH_BrachyPP1_1[is.na(EH_BrachyPP1_1)] <-0
EH_BrachyPP1_Total <- EH_BrachyPP1_0 + EH_BrachyPP1_1
EH_BrachyPP1_Affected <- round((EH_BrachyPP1_1/EH_BrachyPP1_Total)*100, 2)

# Male N and frequency
EH_M_BrachyPP1 <- EH_M %>% count(BrachyPP1)
EH_M_BrachyPP1 <- 
  EH_M_BrachyPP1[(EH_M_BrachyPP1$BrachyPP1 != 6 & EH_M_BrachyPP1$BrachyPP1 != 9),]
EH_M_BrachyPP1_0 <- EH_M_BrachyPP1[1,2]
EH_M_BrachyPP1_0[is.na(EH_M_BrachyPP1_0)] <-0
EH_M_BrachyPP1_1 <- EH_M_BrachyPP1[2,2]
EH_M_BrachyPP1_1[is.na(EH_M_BrachyPP1_1)] <-0
EH_M_BrachyPP1_Total <- EH_M_BrachyPP1_0 + EH_M_BrachyPP1_1
EH_M_BrachyPP1_Affected <- round((EH_M_BrachyPP1_1/EH_M_BrachyPP1_Total)*100, 2)

# Female N and frequency
EH_F_BrachyPP1 <- EH_F %>% count(BrachyPP1)
EH_F_BrachyPP1 <- 
  EH_F_BrachyPP1[(EH_F_BrachyPP1$BrachyPP1 != 6 & EH_F_BrachyPP1$BrachyPP1 != 9),]
EH_F_BrachyPP1_0 <- EH_F_BrachyPP1[1,2]
EH_F_BrachyPP1_0[is.na(EH_F_BrachyPP1_0)] <-0
EH_F_BrachyPP1_1 <- EH_F_BrachyPP1[2,2]
EH_F_BrachyPP1_1[is.na(EH_F_BrachyPP1_1)] <-0
EH_F_BrachyPP1_Total <- EH_F_BrachyPP1_0 + EH_F_BrachyPP1_1
EH_F_BrachyPP1_Affected <- round((EH_F_BrachyPP1_1/EH_F_BrachyPP1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_BrachyPP1_row <- c("BrachyPP1", EH_BrachyPP1_Total, EH_BrachyPP1_1, 
                      EH_BrachyPP1_Affected, EH_M_BrachyPP1_Total, 
                      EH_M_BrachyPP1_1, EH_M_BrachyPP1_Affected, 
                      EH_F_BrachyPP1_Total, EH_F_BrachyPP1_1, 
                      EH_F_BrachyPP1_Affected)

## CalcCubCoal ##
# Overall N and frequency #
EH_CalcCubCoal <- EH_wo_probables %>% count(CalcCubCoal)
EH_CalcCubCoal <- 
  EH_CalcCubCoal[(EH_CalcCubCoal$CalcCubCoal != 6 & EH_CalcCubCoal$CalcCubCoal != 9),]
EH_CalcCubCoal_0 <- EH_CalcCubCoal[1,2]
EH_CalcCubCoal_0[is.na(EH_CalcCubCoal_0)] <-0
EH_CalcCubCoal_1 <- EH_CalcCubCoal[2,2]
EH_CalcCubCoal_1[is.na(EH_CalcCubCoal_1)] <-0
EH_CalcCubCoal_Total <- EH_CalcCubCoal_0 + EH_CalcCubCoal_1
EH_CalcCubCoal_Affected <- 
  round((EH_CalcCubCoal_1/EH_CalcCubCoal_Total)*100, 2)

# Male N and frequency
EH_M_CalcCubCoal <- EH_M %>% count(CalcCubCoal)
EH_M_CalcCubCoal <- 
  EH_M_CalcCubCoal[(EH_M_CalcCubCoal$CalcCubCoal != 6 & EH_M_CalcCubCoal$CalcCubCoal != 9),]
EH_M_CalcCubCoal_0 <- EH_M_CalcCubCoal[1,2]
EH_M_CalcCubCoal_0[is.na(EH_M_CalcCubCoal_0)] <-0
EH_M_CalcCubCoal_1 <- EH_M_CalcCubCoal[2,2]
EH_M_CalcCubCoal_1[is.na(EH_M_CalcCubCoal_1)] <-0
EH_M_CalcCubCoal_Total <- EH_M_CalcCubCoal_0 + EH_M_CalcCubCoal_1
EH_M_CalcCubCoal_Affected <- 
  round((EH_M_CalcCubCoal_1/EH_M_CalcCubCoal_Total)*100, 2)

# Female N and frequency
EH_F_CalcCubCoal <- EH_F %>% count(CalcCubCoal)
EH_F_CalcCubCoal <- 
  EH_F_CalcCubCoal[(EH_F_CalcCubCoal$CalcCubCoal != 6 & EH_F_CalcCubCoal$CalcCubCoal != 9),]
EH_F_CalcCubCoal_0 <- EH_F_CalcCubCoal[1,2]
EH_F_CalcCubCoal_0[is.na(EH_F_CalcCubCoal_0)] <-0
EH_F_CalcCubCoal_1 <- EH_F_CalcCubCoal[2,2]
EH_F_CalcCubCoal_1[is.na(EH_F_CalcCubCoal_1)] <-0
EH_F_CalcCubCoal_Total <- EH_F_CalcCubCoal_0 + EH_F_CalcCubCoal_1
EH_F_CalcCubCoal_Affected <- 
  round((EH_F_CalcCubCoal_1/EH_F_CalcCubCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_CalcCubCoal_row <- c("CalcCubCoal", EH_CalcCubCoal_Total, EH_CalcCubCoal_1, 
                        EH_CalcCubCoal_Affected, EH_M_CalcCubCoal_Total, 
                        EH_M_CalcCubCoal_1, EH_M_CalcCubCoal_Affected, 
                        EH_F_CalcCubCoal_Total, EH_F_CalcCubCoal_1, 
                        EH_F_CalcCubCoal_Affected)

## TaloNavCoal ##
# Overall N and frequency #
EH_TaloNavCoal <- EH_wo_probables %>% count(TaloNavCoal)
EH_TaloNavCoal <- 
  EH_TaloNavCoal[(EH_TaloNavCoal$TaloNavCoal != 6 & EH_TaloNavCoal$TaloNavCoal != 9),]
EH_TaloNavCoal_0 <- EH_TaloNavCoal[1,2]
EH_TaloNavCoal_0[is.na(EH_TaloNavCoal_0)] <-0
EH_TaloNavCoal_1 <- EH_TaloNavCoal[2,2]
EH_TaloNavCoal_1[is.na(EH_TaloNavCoal_1)] <-0
EH_TaloNavCoal_Total <- EH_TaloNavCoal_0 + EH_TaloNavCoal_1
EH_TaloNavCoal_Affected <- round((EH_TaloNavCoal_1/EH_TaloNavCoal_Total)*100, 2)

# Male N and frequency
EH_M_TaloNavCoal <- EH_M %>% count(TaloNavCoal)
EH_M_TaloNavCoal <- 
  EH_M_TaloNavCoal[(EH_M_TaloNavCoal$TaloNavCoal != 6 & EH_M_TaloNavCoal$TaloNavCoal != 9),]
EH_M_TaloNavCoal_0 <- EH_M_TaloNavCoal[1,2]
EH_M_TaloNavCoal_0[is.na(EH_M_TaloNavCoal_0)] <-0
EH_M_TaloNavCoal_1 <- EH_M_TaloNavCoal[2,2]
EH_M_TaloNavCoal_1[is.na(EH_M_TaloNavCoal_1)] <-0
EH_M_TaloNavCoal_Total <- EH_M_TaloNavCoal_0 + EH_M_TaloNavCoal_1
EH_M_TaloNavCoal_Affected <- 
  round((EH_M_TaloNavCoal_1/EH_M_TaloNavCoal_Total)*100, 2)

# Female N and frequency
EH_F_TaloNavCoal <- EH_F %>% count(TaloNavCoal)
EH_F_TaloNavCoal <- 
  EH_F_TaloNavCoal[(EH_F_TaloNavCoal$TaloNavCoal != 6 & EH_F_TaloNavCoal$TaloNavCoal != 9),]
EH_F_TaloNavCoal_0 <- EH_F_TaloNavCoal[1,2]
EH_F_TaloNavCoal_0[is.na(EH_F_TaloNavCoal_0)] <-0
EH_F_TaloNavCoal_1 <- EH_F_TaloNavCoal[2,2]
EH_F_TaloNavCoal_1[is.na(EH_F_TaloNavCoal_1)] <-0
EH_F_TaloNavCoal_Total <- EH_F_TaloNavCoal_0 + EH_F_TaloNavCoal_1
EH_F_TaloNavCoal_Affected <- 
  round((EH_F_TaloNavCoal_1/EH_F_TaloNavCoal_Total)*100, 2)

# overview of Taloulated frequencies for table (see below) #
EH_TaloNavCoal_row <- c("TaloNavCoal", EH_TaloNavCoal_Total, EH_TaloNavCoal_1, 
                        EH_TaloNavCoal_Affected, EH_M_TaloNavCoal_Total, 
                        EH_M_TaloNavCoal_1, EH_M_TaloNavCoal_Affected, 
                        EH_F_TaloNavCoal_Total, EH_F_TaloNavCoal_1, 
                        EH_F_TaloNavCoal_Affected)

## CF2CF3Coal ##
# Overall N and frequency #
EH_CF2CF3Coal <- EH_wo_probables %>% count(CF2CF3Coal)
EH_CF2CF3Coal <- 
  EH_CF2CF3Coal[(EH_CF2CF3Coal$CF2CF3Coal != 6 & EH_CF2CF3Coal$CF2CF3Coal != 9),]
EH_CF2CF3Coal_0 <- EH_CF2CF3Coal[1,2]
EH_CF2CF3Coal_0[is.na(EH_CF2CF3Coal_0)] <-0
EH_CF2CF3Coal_1 <- EH_CF2CF3Coal[2,2]
EH_CF2CF3Coal_1[is.na(EH_CF2CF3Coal_1)] <-0
EH_CF2CF3Coal_Total <- EH_CF2CF3Coal_0 + EH_CF2CF3Coal_1
EH_CF2CF3Coal_Affected <- round((EH_CF2CF3Coal_1/EH_CF2CF3Coal_Total)*100, 2)

# Male N and frequency
EH_M_CF2CF3Coal <- EH_M %>% count(CF2CF3Coal)
EH_M_CF2CF3Coal <- 
  EH_M_CF2CF3Coal[(EH_M_CF2CF3Coal$CF2CF3Coal != 6 & EH_M_CF2CF3Coal$CF2CF3Coal != 9),]
EH_M_CF2CF3Coal_0 <- EH_M_CF2CF3Coal[1,2]
EH_M_CF2CF3Coal_0[is.na(EH_M_CF2CF3Coal_0)] <-0
EH_M_CF2CF3Coal_1 <- EH_M_CF2CF3Coal[2,2]
EH_M_CF2CF3Coal_1[is.na(EH_M_CF2CF3Coal_1)] <-0
EH_M_CF2CF3Coal_Total <- EH_M_CF2CF3Coal_0 + EH_M_CF2CF3Coal_1
EH_M_CF2CF3Coal_Affected <- 
  round((EH_M_CF2CF3Coal_1/EH_M_CF2CF3Coal_Total)*100, 2)

# Female N and frequency
EH_F_CF2CF3Coal <- EH_F %>% count(CF2CF3Coal)
EH_F_CF2CF3Coal <- 
  EH_F_CF2CF3Coal[(EH_F_CF2CF3Coal$CF2CF3Coal != 6 & EH_F_CF2CF3Coal$CF2CF3Coal != 9),]
EH_F_CF2CF3Coal_0 <- EH_F_CF2CF3Coal[1,2]
EH_F_CF2CF3Coal_0[is.na(EH_F_CF2CF3Coal_0)] <-0
EH_F_CF2CF3Coal_1 <- EH_F_CF2CF3Coal[2,2]
EH_F_CF2CF3Coal_1[is.na(EH_F_CF2CF3Coal_1)] <-0
EH_F_CF2CF3Coal_Total <- EH_F_CF2CF3Coal_0 + EH_F_CF2CF3Coal_1
EH_F_CF2CF3Coal_Affected <- 
  round((EH_F_CF2CF3Coal_1/EH_F_CF2CF3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_CF2CF3Coal_row <- c("CF2CF3Coal", EH_CF2CF3Coal_Total, EH_CF2CF3Coal_1, 
                       EH_CF2CF3Coal_Affected, EH_M_CF2CF3Coal_Total, 
                       EH_M_CF2CF3Coal_1, EH_M_CF2CF3Coal_Affected, 
                       EH_F_CF2CF3Coal_Total, EH_F_CF2CF3Coal_1, 
                       EH_F_CF2CF3Coal_Affected)

## CF3MT3Coal ##
# Overall N and frequency #
EH_CF3MT3Coal <- EH_wo_probables %>% count(CF3MT3Coal)
EH_CF3MT3Coal <- 
  EH_CF3MT3Coal[(EH_CF3MT3Coal$CF3MT3Coal != 6 & EH_CF3MT3Coal$CF3MT3Coal != 9),]
EH_CF3MT3Coal_0 <- EH_CF3MT3Coal[1,2]
EH_CF3MT3Coal_0[is.na(EH_CF3MT3Coal_0)] <-0
EH_CF3MT3Coal_1 <- EH_CF3MT3Coal[2,2]
EH_CF3MT3Coal_1[is.na(EH_CF3MT3Coal_1)] <-0
EH_CF3MT3Coal_Total <- EH_CF3MT3Coal_0 + EH_CF3MT3Coal_1
EH_CF3MT3Coal_Affected <- round((EH_CF3MT3Coal_1/EH_CF3MT3Coal_Total)*100, 2)

# Male N and frequency
EH_M_CF3MT3Coal <- EH_M %>% count(CF3MT3Coal)
EH_M_CF3MT3Coal <- 
  EH_M_CF3MT3Coal[(EH_M_CF3MT3Coal$CF3MT3Coal != 6 & EH_M_CF3MT3Coal$CF3MT3Coal != 9),]
EH_M_CF3MT3Coal_0 <- EH_M_CF3MT3Coal[1,2]
EH_M_CF3MT3Coal_0[is.na(EH_M_CF3MT3Coal_0)] <-0
EH_M_CF3MT3Coal_1 <- EH_M_CF3MT3Coal[2,2]
EH_M_CF3MT3Coal_1[is.na(EH_M_CF3MT3Coal_1)] <-0
EH_M_CF3MT3Coal_Total <- EH_M_CF3MT3Coal_0 + EH_M_CF3MT3Coal_1
EH_M_CF3MT3Coal_Affected <- 
  round((EH_M_CF3MT3Coal_1/EH_M_CF3MT3Coal_Total)*100, 2)

# Female N and frequency
EH_F_CF3MT3Coal <- EH_F %>% count(CF3MT3Coal)
EH_F_CF3MT3Coal <- 
  EH_F_CF3MT3Coal[(EH_F_CF3MT3Coal$CF3MT3Coal != 6 & EH_F_CF3MT3Coal$CF3MT3Coal != 9),]
EH_F_CF3MT3Coal_0 <- EH_F_CF3MT3Coal[1,2]
EH_F_CF3MT3Coal_0[is.na(EH_F_CF3MT3Coal_0)] <-0
EH_F_CF3MT3Coal_1 <- EH_F_CF3MT3Coal[2,2]
EH_F_CF3MT3Coal_1[is.na(EH_F_CF3MT3Coal_1)] <-0
EH_F_CF3MT3Coal_Total <- EH_F_CF3MT3Coal_0 + EH_F_CF3MT3Coal_1
EH_F_CF3MT3Coal_Affected <- 
  round((EH_F_CF3MT3Coal_1/EH_F_CF3MT3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_CF3MT3Coal_row <- c("CF3MT3Coal", EH_CF3MT3Coal_Total, EH_CF3MT3Coal_1, 
                       EH_CF3MT3Coal_Affected, EH_M_CF3MT3Coal_Total, 
                       EH_M_CF3MT3Coal_1, EH_M_CF3MT3Coal_Affected, 
                       EH_F_CF3MT3Coal_Total, EH_F_CF3MT3Coal_1, 
                       EH_F_CF3MT3Coal_Affected)
## OsIntermet ##
# Overall N and frequency #
EH_OsIntermet <- EH_wo_probables %>% count(OsIntermet)
EH_OsIntermet <- 
  EH_OsIntermet[(EH_OsIntermet$OsIntermet != 6 & EH_OsIntermet$OsIntermet != 9),]
EH_OsIntermet_0 <- EH_OsIntermet[1,2]
EH_OsIntermet_0[is.na(EH_OsIntermet_0)] <-0
EH_OsIntermet_1 <- EH_OsIntermet[2,2]
EH_OsIntermet_1[is.na(EH_OsIntermet_1)] <-0
EH_OsIntermet_Total <- EH_OsIntermet_0 + EH_OsIntermet_1
EH_OsIntermet_Affected <- round((EH_OsIntermet_1/EH_OsIntermet_Total)*100, 2)

# Male N and frequency
EH_M_OsIntermet <- EH_M %>% count(OsIntermet)
EH_M_OsIntermet <- 
  EH_M_OsIntermet[(EH_M_OsIntermet$OsIntermet != 6 & EH_M_OsIntermet$OsIntermet != 9),]
EH_M_OsIntermet_0 <- EH_M_OsIntermet[1,2]
EH_M_OsIntermet_0[is.na(EH_M_OsIntermet_0)] <-0
EH_M_OsIntermet_1 <- EH_M_OsIntermet[2,2]
EH_M_OsIntermet_1[is.na(EH_M_OsIntermet_1)] <-0
EH_M_OsIntermet_Total <- EH_M_OsIntermet_0 + EH_M_OsIntermet_1
EH_M_OsIntermet_Affected <- 
  round((EH_M_OsIntermet_1/EH_M_OsIntermet_Total)*100, 2)

# Female N and frequency
EH_F_OsIntermet <- EH_F %>% count(OsIntermet)
EH_F_OsIntermet <- 
  EH_F_OsIntermet[(EH_F_OsIntermet$OsIntermet != 6 & EH_F_OsIntermet$OsIntermet != 9),]
EH_F_OsIntermet_0 <- EH_F_OsIntermet[1,2]
EH_F_OsIntermet_0[is.na(EH_F_OsIntermet_0)] <-0
EH_F_OsIntermet_1 <- EH_F_OsIntermet[2,2]
EH_F_OsIntermet_1[is.na(EH_F_OsIntermet_1)] <-0
EH_F_OsIntermet_Total <- EH_F_OsIntermet_0 + EH_F_OsIntermet_1
EH_F_OsIntermet_Affected <- 
  round((EH_F_OsIntermet_1/EH_F_OsIntermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_OsIntermet_row <- c("OsIntermet", EH_OsIntermet_Total, EH_OsIntermet_1, 
                       EH_OsIntermet_Affected, EH_M_OsIntermet_Total, 
                       EH_M_OsIntermet_1, EH_M_OsIntermet_Affected, 
                       EH_F_OsIntermet_Total, EH_F_OsIntermet_1, 
                       EH_F_OsIntermet_Affected)
# OsIntermet unclustered: CF1Intermet #
# Overall N and frequency #
EH_CF1Intermet <- EH_wo_probables %>% count(CF1Intermet)
EH_CF1Intermet <- 
  EH_CF1Intermet[(EH_CF1Intermet$CF1Intermet != 6 & EH_CF1Intermet$CF1Intermet != 9),]
EH_CF1Intermet_0 <- EH_CF1Intermet[1,2]
EH_CF1Intermet_0[is.na(EH_CF1Intermet_0)] <-0
EH_CF1Intermet_1 <- EH_CF1Intermet[2,2]
EH_CF1Intermet_1[is.na(EH_CF1Intermet_1)] <-0
EH_CF1Intermet_Total <- na.omit(EH_CF1Intermet_0) + na.omit(EH_CF1Intermet_1)
EH_CF1Intermet_Affected <- round((EH_CF1Intermet_1/EH_CF1Intermet_Total)*100, 2)

# Male N and frequency
EH_M_CF1Intermet <- EH_M %>% count(CF1Intermet)
EH_M_CF1Intermet <- 
  EH_M_CF1Intermet[(EH_M_CF1Intermet$CF1Intermet != 6 & EH_M_CF1Intermet$CF1Intermet != 9),]
EH_M_CF1Intermet_0 <- EH_M_CF1Intermet[1,2]
EH_M_CF1Intermet_0[is.na(EH_M_CF1Intermet_0)] <-0
EH_M_CF1Intermet_1 <- EH_M_CF1Intermet[2,2]
EH_M_CF1Intermet_1[is.na(EH_M_CF1Intermet_1)] <-0
EH_M_CF1Intermet_Total <- EH_M_CF1Intermet_0 + EH_M_CF1Intermet_1
EH_M_CF1Intermet_Affected <- 
  round((EH_M_CF1Intermet_1/EH_M_CF1Intermet_Total)*100, 2)

# Female N and frequency
EH_F_CF1Intermet <- EH_F %>% count(CF1Intermet)
EH_F_CF1Intermet <- 
  EH_F_CF1Intermet[(EH_F_CF1Intermet$CF1Intermet != 6 & EH_F_CF1Intermet$CF1Intermet != 9),]
EH_F_CF1Intermet_0 <- EH_F_CF1Intermet[1,2]
EH_F_CF1Intermet_0[is.na(EH_F_CF1Intermet_0)] <-0
EH_F_CF1Intermet_1 <- EH_F_CF1Intermet[2,2]
EH_F_CF1Intermet_1[is.na(EH_F_CF1Intermet_1)] <-0
EH_F_CF1Intermet_Total <- EH_F_CF1Intermet_0 + EH_F_CF1Intermet_1
EH_F_CF1Intermet_Affected <- 
  round((EH_F_CF1Intermet_1/EH_F_CF1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_CF1Intermet_row <- c("CF1Intermet", EH_CF1Intermet_Total, EH_CF1Intermet_1, 
                        EH_CF1Intermet_Affected, EH_M_CF1Intermet_Total, 
                        EH_M_CF1Intermet_1, EH_M_CF1Intermet_Affected, 
                        EH_F_CF1Intermet_Total, EH_F_CF1Intermet_1, 
                        EH_F_CF1Intermet_Affected)

# OsIntermet unclustered: MT1Intermet #
# Overall N and frequency #
EH_MT1Intermet <- EH_wo_probables %>% count(MT1Intermet)
EH_MT1Intermet <- 
  EH_MT1Intermet[(EH_MT1Intermet$MT1Intermet != 6 & EH_MT1Intermet$MT1Intermet != 9),]
EH_MT1Intermet_0 <- EH_MT1Intermet[1,2]
EH_MT1Intermet_0[is.na(EH_MT1Intermet_0)] <-0
EH_MT1Intermet_1 <- EH_MT1Intermet[2,2]
EH_MT1Intermet_1[is.na(EH_MT1Intermet_1)] <-0
EH_MT1Intermet_Total <- na.omit(EH_MT1Intermet_0) + na.omit(EH_MT1Intermet_1)
EH_MT1Intermet_Affected <- round((EH_MT1Intermet_1/EH_MT1Intermet_Total)*100, 2)

# Male N and frequency
EH_M_MT1Intermet <- EH_M %>% count(MT1Intermet)
EH_M_MT1Intermet <- 
  EH_M_MT1Intermet[(EH_M_MT1Intermet$MT1Intermet != 6 & EH_M_MT1Intermet$MT1Intermet != 9),]
EH_M_MT1Intermet_0 <- EH_M_MT1Intermet[1,2]
EH_M_MT1Intermet_0[is.na(EH_M_MT1Intermet_0)] <-0
EH_M_MT1Intermet_1 <- EH_M_MT1Intermet[2,2]
EH_M_MT1Intermet_1[is.na(EH_M_MT1Intermet_1)] <-0
EH_M_MT1Intermet_Total <- EH_M_MT1Intermet_0 + EH_M_MT1Intermet_1
EH_M_MT1Intermet_Affected <- 
  round((EH_M_MT1Intermet_1/EH_M_MT1Intermet_Total)*100, 2)

# Female N and frequency
EH_F_MT1Intermet <- EH_F %>% count(MT1Intermet)
EH_F_MT1Intermet <- 
  EH_F_MT1Intermet[(EH_F_MT1Intermet$MT1Intermet != 6 & EH_F_MT1Intermet$MT1Intermet != 9),]
EH_F_MT1Intermet_0 <- EH_F_MT1Intermet[1,2]
EH_F_MT1Intermet_0[is.na(EH_F_MT1Intermet_0)] <-0
EH_F_MT1Intermet_1 <- EH_F_MT1Intermet[2,2]
EH_F_MT1Intermet_1[is.na(EH_F_MT1Intermet_1)] <-0
EH_F_MT1Intermet_Total <- EH_F_MT1Intermet_0 + EH_F_MT1Intermet_1
EH_F_MT1Intermet_Affected <- 
  round((EH_F_MT1Intermet_1/EH_F_MT1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_MT1Intermet_row <- c("MT1Intermet", EH_MT1Intermet_Total, EH_MT1Intermet_1, 
                        EH_MT1Intermet_Affected, EH_M_MT1Intermet_Total, 
                        EH_M_MT1Intermet_1, EH_M_MT1Intermet_Affected, 
                        EH_F_MT1Intermet_Total, EH_F_MT1Intermet_1, 
                        EH_F_MT1Intermet_Affected)

# OsIntermet unclustered: MT2Intermet #
# Overall N and frequency #
EH_MT2Intermet <- EH_wo_probables %>% count(MT2Intermet)
EH_MT2Intermet <- 
  EH_MT2Intermet[(EH_MT2Intermet$MT2Intermet != 6 & EH_MT2Intermet$MT2Intermet != 9),]
EH_MT2Intermet_0 <- EH_MT2Intermet[1,2]
EH_MT2Intermet_0[is.na(EH_MT2Intermet_0)] <-0
EH_MT2Intermet_1 <- EH_MT2Intermet[2,2]
EH_MT2Intermet_1[is.na(EH_MT2Intermet_1)] <-0
EH_MT2Intermet_Total <- na.omit(EH_MT2Intermet_0) + na.omit(EH_MT2Intermet_1)
EH_MT2Intermet_Affected <- round((EH_MT2Intermet_1/EH_MT2Intermet_Total)*100, 2)

# Male N and frequency
EH_M_MT2Intermet <- EH_M %>% count(MT2Intermet)
EH_M_MT2Intermet <- 
  EH_M_MT2Intermet[(EH_M_MT2Intermet$MT2Intermet != 6 & EH_M_MT2Intermet$MT2Intermet != 9),]
EH_M_MT2Intermet_0 <- EH_M_MT2Intermet[1,2]
EH_M_MT2Intermet_0[is.na(EH_M_MT2Intermet_0)] <-0
EH_M_MT2Intermet_1 <- EH_M_MT2Intermet[2,2]
EH_M_MT2Intermet_1[is.na(EH_M_MT2Intermet_1)] <-0
EH_M_MT2Intermet_Total <- EH_M_MT2Intermet_0 + EH_M_MT2Intermet_1
EH_M_MT2Intermet_Affected <- 
  round((EH_M_MT2Intermet_1/EH_M_MT2Intermet_Total)*100, 2)

# Female N and frequency
EH_F_MT2Intermet <- EH_F %>% count(MT2Intermet)
EH_F_MT2Intermet <- 
  EH_F_MT2Intermet[(EH_F_MT2Intermet$MT2Intermet != 6 & EH_F_MT2Intermet$MT2Intermet != 9),]
EH_F_MT2Intermet_0 <- EH_F_MT2Intermet[1,2]
EH_F_MT2Intermet_0[is.na(EH_F_MT2Intermet_0)] <-0
EH_F_MT2Intermet_1 <- EH_F_MT2Intermet[2,2]
EH_F_MT2Intermet_1[is.na(EH_F_MT2Intermet_1)] <-0
EH_F_MT2Intermet_Total <- EH_F_MT2Intermet_0 + EH_F_MT2Intermet_1
EH_F_MT2Intermet_Affected <- 
  round((EH_F_MT2Intermet_1/EH_F_MT2Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_MT2Intermet_row <- c("MT2Intermet", EH_MT2Intermet_Total, EH_MT2Intermet_1, 
                        EH_MT2Intermet_Affected, EH_M_MT2Intermet_Total, 
                        EH_M_MT2Intermet_1, EH_M_MT2Intermet_Affected, 
                        EH_F_MT2Intermet_Total, EH_F_MT2Intermet_1, 
                        EH_F_MT2Intermet_Affected)


## AnkleCoal ##
# Overall N and frequency #
EH_AnkleCoal <- EH_wo_probables %>% count(AnkleCoal)
EH_AnkleCoal <- 
  EH_AnkleCoal[(EH_AnkleCoal$AnkleCoal != 6 & EH_AnkleCoal$AnkleCoal != 9),]
EH_AnkleCoal_0 <- EH_AnkleCoal[1,2]; EH_AnkleCoal_0[is.na(EH_AnkleCoal_0)] <-0
EH_AnkleCoal_1 <- EH_AnkleCoal[2,2]; EH_AnkleCoal_1[is.na(EH_AnkleCoal_1)] <-0
EH_AnkleCoal_Total <- na.omit(EH_AnkleCoal_0) + na.omit(EH_AnkleCoal_1)
EH_AnkleCoal_Affected <- round((EH_AnkleCoal_1/EH_AnkleCoal_Total)*100, 2)

# Male N and frequency
EH_M_AnkleCoal <- EH_M %>% count(AnkleCoal)
EH_M_AnkleCoal <- 
  EH_M_AnkleCoal[(EH_M_AnkleCoal$AnkleCoal != 6 & EH_M_AnkleCoal$AnkleCoal != 9),]
EH_M_AnkleCoal_0 <- EH_M_AnkleCoal[1,2]
EH_M_AnkleCoal_0[is.na(EH_M_AnkleCoal_0)] <-0
EH_M_AnkleCoal_1 <- EH_M_AnkleCoal[2,2]
EH_M_AnkleCoal_1[is.na(EH_M_AnkleCoal_1)] <-0
EH_M_AnkleCoal_Total <- EH_M_AnkleCoal_0 + EH_M_AnkleCoal_1
EH_M_AnkleCoal_Affected <- round((EH_M_AnkleCoal_1/EH_M_AnkleCoal_Total)*100, 2)

# Female N and frequency
EH_F_AnkleCoal <- EH_F %>% count(AnkleCoal)
EH_F_AnkleCoal <- 
  EH_F_AnkleCoal[(EH_F_AnkleCoal$AnkleCoal != 6 & EH_F_AnkleCoal$AnkleCoal != 9),]
EH_F_AnkleCoal_0 <- EH_F_AnkleCoal[1,2]
EH_F_AnkleCoal_0[is.na(EH_F_AnkleCoal_0)] <-0
EH_F_AnkleCoal_1 <- EH_F_AnkleCoal[2,2]
EH_F_AnkleCoal_1[is.na(EH_F_AnkleCoal_1)] <-0
EH_F_AnkleCoal_Total <- EH_F_AnkleCoal_0 + EH_F_AnkleCoal_1
EH_F_AnkleCoal_Affected <- round((EH_F_AnkleCoal_1/EH_F_AnkleCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_AnkleCoal_row <- c("AnkleCoal", EH_AnkleCoal_Total, EH_AnkleCoal_1,
                      EH_AnkleCoal_Affected, EH_M_AnkleCoal_Total, 
                      EH_M_AnkleCoal_1, EH_M_AnkleCoal_Affected, EH_F_AnkleCoal_Total, 
                      EH_F_AnkleCoal_1, EH_F_AnkleCoal_Affected)

# AnkleCoal unclustered: CalcNavCoal #
# Overall N and frequency #
EH_CalcNavCoal <- EH_wo_probables %>% count(CalcNavCoal)
EH_CalcNavCoal <- 
  EH_CalcNavCoal[(EH_CalcNavCoal$CalcNavCoal != 6 & EH_CalcNavCoal$CalcNavCoal != 9),]
EH_CalcNavCoal_0 <- EH_CalcNavCoal[1,2]
EH_CalcNavCoal_0[is.na(EH_CalcNavCoal_0)] <-0
EH_CalcNavCoal_1 <- EH_CalcNavCoal[2,2]
EH_CalcNavCoal_1[is.na(EH_CalcNavCoal_1)] <-0
EH_CalcNavCoal_Total <- na.omit(EH_CalcNavCoal_0) + na.omit(EH_CalcNavCoal_1)
EH_CalcNavCoal_Affected <- round((EH_CalcNavCoal_1/EH_CalcNavCoal_Total)*100, 2)

# Male N and frequency
EH_M_CalcNavCoal <- EH_M %>% count(CalcNavCoal)
EH_M_CalcNavCoal <- 
  EH_M_CalcNavCoal[(EH_M_CalcNavCoal$CalcNavCoal != 6 & EH_M_CalcNavCoal$CalcNavCoal != 9),]
EH_M_CalcNavCoal_0 <- EH_M_CalcNavCoal[1,2]
EH_M_CalcNavCoal_0[is.na(EH_M_CalcNavCoal_0)] <-0
EH_M_CalcNavCoal_1 <- EH_M_CalcNavCoal[2,2]
EH_M_CalcNavCoal_1[is.na(EH_M_CalcNavCoal_1)] <-0
EH_M_CalcNavCoal_Total <- EH_M_CalcNavCoal_0 + EH_M_CalcNavCoal_1
EH_M_CalcNavCoal_Affected <- 
  round((EH_M_CalcNavCoal_1/EH_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
EH_F_CalcNavCoal <- EH_F %>% count(CalcNavCoal)
EH_F_CalcNavCoal <- 
  EH_F_CalcNavCoal[(EH_F_CalcNavCoal$CalcNavCoal != 6 & EH_F_CalcNavCoal$CalcNavCoal != 9),]
EH_F_CalcNavCoal_0 <- EH_F_CalcNavCoal[1,2]
EH_F_CalcNavCoal_0[is.na(EH_F_CalcNavCoal_0)] <-0
EH_F_CalcNavCoal_1 <- EH_F_CalcNavCoal[2,2]
EH_F_CalcNavCoal_1[is.na(EH_F_CalcNavCoal_1)] <-0
EH_F_CalcNavCoal_Total <- EH_F_CalcNavCoal_0 + EH_F_CalcNavCoal_1
EH_F_CalcNavCoal_Affected <- 
  round((EH_F_CalcNavCoal_1/EH_F_CalcNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_CalcNavCoal_row <- c("CalcNavCoal", EH_CalcNavCoal_Total, EH_CalcNavCoal_1, 
                        EH_CalcNavCoal_Affected, EH_M_CalcNavCoal_Total, 
                        EH_M_CalcNavCoal_1, EH_M_CalcNavCoal_Affected, 
                        EH_F_CalcNavCoal_Total, EH_F_CalcNavCoal_1, 
                        EH_F_CalcNavCoal_Affected)

# AnkleCoal unclustered: TaloCalcCoal #
# Overall N and frequency #
EH_TaloCalcCoal <- EH_wo_probables %>% count(TaloCalcCoal)
EH_TaloCalcCoal <- 
  EH_TaloCalcCoal[(EH_TaloCalcCoal$TaloCalcCoal != 6 & EH_TaloCalcCoal$TaloCalcCoal != 9),]
EH_TaloCalcCoal_0 <- EH_TaloCalcCoal[1,2]
EH_TaloCalcCoal_0[is.na(EH_TaloCalcCoal_0)] <-0
EH_TaloCalcCoal_1 <- EH_TaloCalcCoal[2,2]
EH_TaloCalcCoal_1[is.na(EH_TaloCalcCoal_1)] <-0
EH_TaloCalcCoal_Total <- na.omit(EH_TaloCalcCoal_0) + na.omit(EH_TaloCalcCoal_1)
EH_TaloCalcCoal_Affected <- 
  round((EH_TaloCalcCoal_1/EH_TaloCalcCoal_Total)*100, 2)

# Male N and frequency
EH_M_TaloCalcCoal <- EH_M %>% count(TaloCalcCoal)
EH_M_TaloCalcCoal <- 
  EH_M_TaloCalcCoal[(EH_M_TaloCalcCoal$TaloCalcCoal != 6 & EH_M_TaloCalcCoal$TaloCalcCoal != 9),]
EH_M_TaloCalcCoal_0 <- EH_M_TaloCalcCoal[1,2]
EH_M_TaloCalcCoal_0[is.na(EH_M_TaloCalcCoal_0)] <-0
EH_M_TaloCalcCoal_1 <- EH_M_TaloCalcCoal[2,2]
EH_M_TaloCalcCoal_1[is.na(EH_M_TaloCalcCoal_1)] <-0
EH_M_TaloCalcCoal_Total <- EH_M_TaloCalcCoal_0 + EH_M_TaloCalcCoal_1
EH_M_TaloCalcCoal_Affected <- 
  round((EH_M_TaloCalcCoal_1/EH_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
EH_F_TaloCalcCoal <- EH_F %>% count(TaloCalcCoal)
EH_F_TaloCalcCoal <- 
  EH_F_TaloCalcCoal[(EH_F_TaloCalcCoal$TaloCalcCoal != 6 & EH_F_TaloCalcCoal$TaloCalcCoal != 9),]
EH_F_TaloCalcCoal_0 <- EH_F_TaloCalcCoal[1,2]
EH_F_TaloCalcCoal_0[is.na(EH_F_TaloCalcCoal_0)] <-0
EH_F_TaloCalcCoal_1 <- EH_F_TaloCalcCoal[2,2]
EH_F_TaloCalcCoal_1[is.na(EH_F_TaloCalcCoal_1)] <-0
EH_F_TaloCalcCoal_Total <- EH_F_TaloCalcCoal_0 + EH_F_TaloCalcCoal_1
EH_F_TaloCalcCoal_Affected <- 
  round((EH_F_TaloCalcCoal_1/EH_F_TaloCalcCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
EH_TaloCalcCoal_row <- c("TaloCalcCoal", EH_TaloCalcCoal_Total, EH_TaloCalcCoal_1, 
                         EH_TaloCalcCoal_Affected, EH_M_TaloCalcCoal_Total, 
                         EH_M_TaloCalcCoal_1, EH_M_TaloCalcCoal_Affected, 
                         EH_F_TaloCalcCoal_Total, EH_F_TaloCalcCoal_1, 
                         EH_F_TaloCalcCoal_Affected)


## Table with overview of trait frequencies for EH ##
EH_freq <- data.frame(rbind(EH_AccessNav_row, EH_BrachyD_row, EH_BrachyMT1_row,
                            EH_BrachyMT4_row, EH_BrachyPP1_row, EH_CalcCubCoal_row,
                            EH_TaloNavCoal_row, EH_CF2CF3Coal_row, EH_CF3MT3Coal_row,
                            EH_OsIntermet_row, EH_CF1Intermet_row, EH_MT1Intermet_row,
                            EH_MT2Intermet_row, EH_AnkleCoal_row, EH_CalcNavCoal_row, 
                            EH_TaloCalcCoal_row))
names(EH_freq) <- c("Anomaly", "Overall N (EH)", "Overall affected (EH)", 
                    "Overall % (EH)", "Male N (EH)", "Male affected (EH)",
                    "Male % (EH)", "Female N (EH)", "Female affected (EH)", 
                    "Female % (EH)")

### The Zwolle collection ###


## AccessNav ##
# Overall N and frequency #
ZW_AccessNav <- ZW_wo_probables %>% count(AccessNav)
ZW_AccessNav <- 
  ZW_AccessNav[(ZW_AccessNav$AccessNav != 6 & ZW_AccessNav$AccessNav != 9),]
ZW_AccessNav_0 <- ZW_AccessNav[1,2]
ZW_AccessNav_0[is.na(ZW_AccessNav_0)] <-0
ZW_AccessNav_1 <- ZW_AccessNav[2,2]
ZW_AccessNav_1[is.na(ZW_AccessNav_1)] <-0
ZW_AccessNav_Total <- ZW_AccessNav_0 + ZW_AccessNav_1
ZW_AccessNav_Affected <- round((ZW_AccessNav_1/ZW_AccessNav_Total)*100, 2)

# Male N and frequency
ZW_M_AccessNav <- ZW_M %>% count(AccessNav)
ZW_M_AccessNav <- 
  ZW_M_AccessNav[(ZW_M_AccessNav$AccessNav != 6 & ZW_M_AccessNav$AccessNav != 9),]
ZW_M_AccessNav_0 <- ZW_M_AccessNav[1,2]
ZW_M_AccessNav_0[is.na(ZW_M_AccessNav_0)] <-0
ZW_M_AccessNav_1 <- ZW_M_AccessNav[2,2]
ZW_M_AccessNav_1[is.na(ZW_M_AccessNav_1)] <-0
ZW_M_AccessNav_Total <- ZW_M_AccessNav_0 + ZW_M_AccessNav_1
ZW_M_AccessNav_Affected <- round((ZW_M_AccessNav_1/ZW_M_AccessNav_Total)*100, 2)

# Female N and frequency
ZW_F_AccessNav <- ZW_F %>% count(AccessNav)
ZW_F_AccessNav <- 
  ZW_F_AccessNav[(ZW_F_AccessNav$AccessNav != 6 & ZW_F_AccessNav$AccessNav != 9),]
ZW_F_AccessNav_0 <- ZW_F_AccessNav[1,2]
ZW_F_AccessNav_0[is.na(ZW_F_AccessNav_0)] <-0
ZW_F_AccessNav_1 <- ZW_F_AccessNav[2,2]
ZW_F_AccessNav_1[is.na(ZW_F_AccessNav_1)] <-0
ZW_F_AccessNav_Total <- ZW_F_AccessNav_0 + ZW_F_AccessNav_1
ZW_F_AccessNav_Affected <- round((ZW_F_AccessNav_1/ZW_F_AccessNav_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_AccessNav_row <- c("AccessNav", ZW_AccessNav_Total, ZW_AccessNav_1, 
                      ZW_AccessNav_Affected, ZW_M_AccessNav_Total, 
                      ZW_M_AccessNav_1, ZW_M_AccessNav_Affected, 
                      ZW_F_AccessNav_Total, ZW_F_AccessNav_1, 
                      ZW_F_AccessNav_Affected)


## BrachyD ##
# Overall N and frequency #
ZW_BrachyD <- ZW_wo_probables %>% count(BrachyD)
ZW_BrachyD <- ZW_BrachyD[(ZW_BrachyD$BrachyD != 6 & ZW_BrachyD$BrachyD != 9),]
ZW_BrachyD_0 <- ZW_BrachyD[1,2]; ZW_BrachyD_0[is.na(ZW_BrachyD_0)] <-0
ZW_BrachyD_1 <- ZW_BrachyD[2,2]; ZW_BrachyD_1[is.na(ZW_BrachyD_1)] <-0
ZW_BrachyD_Total <- ZW_BrachyD_0 + ZW_BrachyD_1
ZW_BrachyD_Affected <- round((ZW_BrachyD_1/ZW_BrachyD_Total)*100, 2)

# Male N and frequency
ZW_M_BrachyD <- ZW_M %>% count(BrachyD)
ZW_M_BrachyD <- 
  ZW_M_BrachyD[(ZW_M_BrachyD$BrachyD != 6 & ZW_M_BrachyD$BrachyD != 9),]
ZW_M_BrachyD_0 <- ZW_M_BrachyD[1,2]; ZW_M_BrachyD_0[is.na(ZW_M_BrachyD_0)] <-0
ZW_M_BrachyD_1 <- ZW_M_BrachyD[2,2]; ZW_M_BrachyD_1[is.na(ZW_M_BrachyD_1)] <-0
ZW_M_BrachyD_Total <- ZW_M_BrachyD_0 + ZW_M_BrachyD_1
ZW_M_BrachyD_Affected <- round((ZW_M_BrachyD_1/ZW_M_BrachyD_Total)*100, 2)

# Female N and frequency
ZW_F_BrachyD <- ZW_F %>% count(BrachyD)
ZW_F_BrachyD <- 
  ZW_F_BrachyD[(ZW_F_BrachyD$BrachyD != 6 & ZW_F_BrachyD$BrachyD != 9),]
ZW_F_BrachyD_0 <- ZW_F_BrachyD[1,2]; ZW_F_BrachyD_0[is.na(ZW_F_BrachyD_0)] <-0
ZW_F_BrachyD_1 <- ZW_F_BrachyD[2,2]; ZW_F_BrachyD_1[is.na(ZW_F_BrachyD_1)] <-0
ZW_F_BrachyD_Total <- ZW_F_BrachyD_0 + ZW_F_BrachyD_1
ZW_F_BrachyD_Affected <- round((ZW_F_BrachyD_1/ZW_F_BrachyD_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_BrachyD_row <- c("BrachyD", ZW_BrachyD_Total, ZW_BrachyD_1, 
                    ZW_BrachyD_Affected, ZW_M_BrachyD_Total, 
                    ZW_M_BrachyD_1, ZW_M_BrachyD_Affected, 
                    ZW_F_BrachyD_Total, ZW_F_BrachyD_1, 
                    ZW_F_BrachyD_Affected)


## BrachyMT1 ##
# Overall N and frequency #
ZW_BrachyMT1 <- ZW_wo_probables %>% count(BrachyMT1)
ZW_BrachyMT1 <- 
  ZW_BrachyMT1[(ZW_BrachyMT1$BrachyMT1 != 6 & ZW_BrachyMT1$BrachyMT1 != 9),]
ZW_BrachyMT1_0 <- ZW_BrachyMT1[1,2]; ZW_BrachyMT1_0[is.na(ZW_BrachyMT1_0)] <-0
ZW_BrachyMT1_1 <- ZW_BrachyMT1[2,2]; ZW_BrachyMT1_1[is.na(ZW_BrachyMT1_1)] <-0
ZW_BrachyMT1_Total <- ZW_BrachyMT1_0 + ZW_BrachyMT1_1
ZW_BrachyMT1_Affected <- round((ZW_BrachyMT1_1/ZW_BrachyMT1_Total)*100, 2)

# Male N and frequency
ZW_M_BrachyMT1 <- ZW_M %>% count(BrachyMT1)
ZW_M_BrachyMT1 <- 
  ZW_M_BrachyMT1[(ZW_M_BrachyMT1$BrachyMT1 != 6 & ZW_M_BrachyMT1$BrachyMT1 != 9),]
ZW_M_BrachyMT1_0 <- ZW_M_BrachyMT1[1,2]
ZW_M_BrachyMT1_0[is.na(ZW_M_BrachyMT1_0)] <-0
ZW_M_BrachyMT1_1 <- ZW_M_BrachyMT1[2,2]
ZW_M_BrachyMT1_1[is.na(ZW_M_BrachyMT1_1)] <-0
ZW_M_BrachyMT1_Total <- ZW_M_BrachyMT1_0 + ZW_M_BrachyMT1_1
ZW_M_BrachyMT1_Affected <- round((ZW_M_BrachyMT1_1/ZW_M_BrachyMT1_Total)*100, 2)

# Female N and frequency
ZW_F_BrachyMT1 <- ZW_F %>% count(BrachyMT1)
ZW_F_BrachyMT1 <- 
  ZW_F_BrachyMT1[(ZW_F_BrachyMT1$BrachyMT1 != 6 & ZW_F_BrachyMT1$BrachyMT1 != 9),]
ZW_F_BrachyMT1_0 <- ZW_F_BrachyMT1[1,2]
ZW_F_BrachyMT1_0[is.na(ZW_F_BrachyMT1_0)] <-0
ZW_F_BrachyMT1_1 <- ZW_F_BrachyMT1[2,2]
ZW_F_BrachyMT1_1[is.na(ZW_F_BrachyMT1_1)] <-0
ZW_F_BrachyMT1_Total <- ZW_F_BrachyMT1_0 + ZW_F_BrachyMT1_1
ZW_F_BrachyMT1_Affected <- round((ZW_F_BrachyMT1_1/ZW_F_BrachyMT1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_BrachyMT1_row <- c("BrachyMT1", ZW_BrachyMT1_Total, ZW_BrachyMT1_1, 
                      ZW_BrachyMT1_Affected, ZW_M_BrachyMT1_Total, 
                      ZW_M_BrachyMT1_1, ZW_M_BrachyMT1_Affected, 
                      ZW_F_BrachyMT1_Total, ZW_F_BrachyMT1_1, 
                      ZW_F_BrachyMT1_Affected)
## BrachyMT4 ##
# Overall N and frequency #
ZW_BrachyMT4 <- ZW_wo_probables %>% count(BrachyMT4)
ZW_BrachyMT4 <- 
  ZW_BrachyMT4[(ZW_BrachyMT4$BrachyMT4 != 6 & ZW_BrachyMT4$BrachyMT4 != 9),]
ZW_BrachyMT4_0 <- ZW_BrachyMT4[1,2]
ZW_BrachyMT4_0[is.na(ZW_BrachyMT4_0)] <-0
ZW_BrachyMT4_1 <- ZW_BrachyMT4[2,2]
ZW_BrachyMT4_1[is.na(ZW_BrachyMT4_1)] <-0
ZW_BrachyMT4_Total <- ZW_BrachyMT4_0 + ZW_BrachyMT4_1
ZW_BrachyMT4_Affected <- round((ZW_BrachyMT4_1/ZW_BrachyMT4_Total)*100, 2)

# Male N and frequency
ZW_M_BrachyMT4 <- ZW_M %>% count(BrachyMT4)
ZW_M_BrachyMT4 <- 
  ZW_M_BrachyMT4[(ZW_M_BrachyMT4$BrachyMT4 != 6 & ZW_M_BrachyMT4$BrachyMT4 != 9),]
ZW_M_BrachyMT4_0 <- ZW_M_BrachyMT4[1,2]
ZW_M_BrachyMT4_0[is.na(ZW_M_BrachyMT4_0)] <-0
ZW_M_BrachyMT4_1 <- ZW_M_BrachyMT4[2,2]
ZW_M_BrachyMT4_1[is.na(ZW_M_BrachyMT4_1)] <-0
ZW_M_BrachyMT4_Total <- ZW_M_BrachyMT4_0 + ZW_M_BrachyMT4_1
ZW_M_BrachyMT4_Affected <- round((ZW_M_BrachyMT4_1/ZW_M_BrachyMT4_Total)*100, 2)

# Female N and frequency
ZW_F_BrachyMT4 <- ZW_F %>% count(BrachyMT4)
ZW_F_BrachyMT4 <- 
  ZW_F_BrachyMT4[(ZW_F_BrachyMT4$BrachyMT4 != 6 & ZW_F_BrachyMT4$BrachyMT4 != 9),]
ZW_F_BrachyMT4_0 <- ZW_F_BrachyMT4[1,2]
ZW_F_BrachyMT4_0[is.na(ZW_F_BrachyMT4_0)] <-0
ZW_F_BrachyMT4_1 <- ZW_F_BrachyMT4[2,2]
ZW_F_BrachyMT4_1[is.na(ZW_F_BrachyMT4_1)] <-0
ZW_F_BrachyMT4_Total <- ZW_F_BrachyMT4_0 + ZW_F_BrachyMT4_1
ZW_F_BrachyMT4_Affected <- round((ZW_F_BrachyMT4_1/ZW_F_BrachyMT4_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_BrachyMT4_row <- c("BrachyMT4", ZW_BrachyMT4_Total, ZW_BrachyMT4_1, 
                      ZW_BrachyMT4_Affected, ZW_M_BrachyMT4_Total, 
                      ZW_M_BrachyMT4_1, ZW_M_BrachyMT4_Affected, 
                      ZW_F_BrachyMT4_Total, ZW_F_BrachyMT4_1, 
                      ZW_F_BrachyMT4_Affected)

## BrachyPP1 ##
# Overall N and frequency #
ZW_BrachyPP1 <- ZW_wo_probables %>% count(BrachyPP1)
ZW_BrachyPP1 <- 
  ZW_BrachyPP1[(ZW_BrachyPP1$BrachyPP1 != 6 & ZW_BrachyPP1$BrachyPP1 != 9),]
ZW_BrachyPP1_0 <- ZW_BrachyPP1[1,2]; ZW_BrachyPP1_0[is.na(ZW_BrachyPP1_0)] <-0
ZW_BrachyPP1_1 <- ZW_BrachyPP1[2,2]; ZW_BrachyPP1_1[is.na(ZW_BrachyPP1_1)] <-0
ZW_BrachyPP1_Total <- ZW_BrachyPP1_0 + ZW_BrachyPP1_1
ZW_BrachyPP1_Affected <- round((ZW_BrachyPP1_1/ZW_BrachyPP1_Total)*100, 2)

# Male N and frequency
ZW_M_BrachyPP1 <- ZW_M %>% count(BrachyPP1)
ZW_M_BrachyPP1 <- 
  ZW_M_BrachyPP1[(ZW_M_BrachyPP1$BrachyPP1 != 6 & ZW_M_BrachyPP1$BrachyPP1 != 9),]
ZW_M_BrachyPP1_0 <- ZW_M_BrachyPP1[1,2]
ZW_M_BrachyPP1_0[is.na(ZW_M_BrachyPP1_0)] <-0
ZW_M_BrachyPP1_1 <- ZW_M_BrachyPP1[2,2]
ZW_M_BrachyPP1_1[is.na(ZW_M_BrachyPP1_1)] <-0
ZW_M_BrachyPP1_Total <- ZW_M_BrachyPP1_0 + ZW_M_BrachyPP1_1
ZW_M_BrachyPP1_Affected <- round((ZW_M_BrachyPP1_1/ZW_M_BrachyPP1_Total)*100, 2)

# Female N and frequency
ZW_F_BrachyPP1 <- ZW_F %>% count(BrachyPP1)
ZW_F_BrachyPP1 <- 
  ZW_F_BrachyPP1[(ZW_F_BrachyPP1$BrachyPP1 != 6 & ZW_F_BrachyPP1$BrachyPP1 != 9),]
ZW_F_BrachyPP1_0 <- ZW_F_BrachyPP1[1,2]
ZW_F_BrachyPP1_0[is.na(ZW_F_BrachyPP1_0)] <-0
ZW_F_BrachyPP1_1 <- ZW_F_BrachyPP1[2,2]
ZW_F_BrachyPP1_1[is.na(ZW_F_BrachyPP1_1)] <-0
ZW_F_BrachyPP1_Total <- ZW_F_BrachyPP1_0 + ZW_F_BrachyPP1_1
ZW_F_BrachyPP1_Affected <- round((ZW_F_BrachyPP1_1/ZW_F_BrachyPP1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_BrachyPP1_row <- c("BrachyPP1", ZW_BrachyPP1_Total, ZW_BrachyPP1_1, 
                      ZW_BrachyPP1_Affected, ZW_M_BrachyPP1_Total, 
                      ZW_M_BrachyPP1_1, ZW_M_BrachyPP1_Affected, 
                      ZW_F_BrachyPP1_Total, ZW_F_BrachyPP1_1, 
                      ZW_F_BrachyPP1_Affected)

## CalcCubCoal ##
# Overall N and frequency #
ZW_CalcCubCoal <- ZW_wo_probables %>% count(CalcCubCoal)
ZW_CalcCubCoal <- 
  ZW_CalcCubCoal[(ZW_CalcCubCoal$CalcCubCoal != 6 & ZW_CalcCubCoal$CalcCubCoal != 9),]
ZW_CalcCubCoal_0 <- ZW_CalcCubCoal[1,2]
ZW_CalcCubCoal_0[is.na(ZW_CalcCubCoal_0)] <-0
ZW_CalcCubCoal_1 <- ZW_CalcCubCoal[2,2]
ZW_CalcCubCoal_1[is.na(ZW_CalcCubCoal_1)] <-0
ZW_CalcCubCoal_Total <- ZW_CalcCubCoal_0 + ZW_CalcCubCoal_1
ZW_CalcCubCoal_Affected <- 
  round((ZW_CalcCubCoal_1/ZW_CalcCubCoal_Total)*100, 2)

# Male N and frequency
ZW_M_CalcCubCoal <- ZW_M %>% count(CalcCubCoal)
ZW_M_CalcCubCoal <- 
  ZW_M_CalcCubCoal[(ZW_M_CalcCubCoal$CalcCubCoal != 6 & ZW_M_CalcCubCoal$CalcCubCoal != 9),]
ZW_M_CalcCubCoal_0 <- ZW_M_CalcCubCoal[1,2]
ZW_M_CalcCubCoal_0[is.na(ZW_M_CalcCubCoal_0)] <-0
ZW_M_CalcCubCoal_1 <- ZW_M_CalcCubCoal[2,2]
ZW_M_CalcCubCoal_1[is.na(ZW_M_CalcCubCoal_1)] <-0
ZW_M_CalcCubCoal_Total <- ZW_M_CalcCubCoal_0 + ZW_M_CalcCubCoal_1
ZW_M_CalcCubCoal_Affected <- 
  round((ZW_M_CalcCubCoal_1/ZW_M_CalcCubCoal_Total)*100, 2)

# Female N and frequency
ZW_F_CalcCubCoal <- ZW_F %>% count(CalcCubCoal)
ZW_F_CalcCubCoal <- 
  ZW_F_CalcCubCoal[(ZW_F_CalcCubCoal$CalcCubCoal != 6 & ZW_F_CalcCubCoal$CalcCubCoal != 9),]
ZW_F_CalcCubCoal_0 <- ZW_F_CalcCubCoal[1,2]
ZW_F_CalcCubCoal_0[is.na(ZW_F_CalcCubCoal_0)] <-0
ZW_F_CalcCubCoal_1 <- ZW_F_CalcCubCoal[2,2]
ZW_F_CalcCubCoal_1[is.na(ZW_F_CalcCubCoal_1)] <-0
ZW_F_CalcCubCoal_Total <- ZW_F_CalcCubCoal_0 + ZW_F_CalcCubCoal_1
ZW_F_CalcCubCoal_Affected <- 
  round((ZW_F_CalcCubCoal_1/ZW_F_CalcCubCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_CalcCubCoal_row <- c("CalcCubCoal", ZW_CalcCubCoal_Total, ZW_CalcCubCoal_1, 
                        ZW_CalcCubCoal_Affected, ZW_M_CalcCubCoal_Total, 
                        ZW_M_CalcCubCoal_1, ZW_M_CalcCubCoal_Affected, 
                        ZW_F_CalcCubCoal_Total, ZW_F_CalcCubCoal_1, 
                        ZW_F_CalcCubCoal_Affected)

## TaloNavCoal ##
# Overall N and frequency #
ZW_TaloNavCoal <- ZW_wo_probables %>% count(TaloNavCoal)
ZW_TaloNavCoal <- 
  ZW_TaloNavCoal[(ZW_TaloNavCoal$TaloNavCoal != 6 & ZW_TaloNavCoal$TaloNavCoal != 9),]
ZW_TaloNavCoal_0 <- ZW_TaloNavCoal[1,2]
ZW_TaloNavCoal_0[is.na(ZW_TaloNavCoal_0)] <-0
ZW_TaloNavCoal_1 <- ZW_TaloNavCoal[2,2]
ZW_TaloNavCoal_1[is.na(ZW_TaloNavCoal_1)] <-0
ZW_TaloNavCoal_Total <- ZW_TaloNavCoal_0 + ZW_TaloNavCoal_1
ZW_TaloNavCoal_Affected <- round((ZW_TaloNavCoal_1/ZW_TaloNavCoal_Total)*100, 2)

# Male N and frequency
ZW_M_TaloNavCoal <- ZW_M %>% count(TaloNavCoal)
ZW_M_TaloNavCoal <- 
  ZW_M_TaloNavCoal[(ZW_M_TaloNavCoal$TaloNavCoal != 6 & ZW_M_TaloNavCoal$TaloNavCoal != 9),]
ZW_M_TaloNavCoal_0 <- ZW_M_TaloNavCoal[1,2]
ZW_M_TaloNavCoal_0[is.na(ZW_M_TaloNavCoal_0)] <-0
ZW_M_TaloNavCoal_1 <- ZW_M_TaloNavCoal[2,2]
ZW_M_TaloNavCoal_1[is.na(ZW_M_TaloNavCoal_1)] <-0
ZW_M_TaloNavCoal_Total <- ZW_M_TaloNavCoal_0 + ZW_M_TaloNavCoal_1
ZW_M_TaloNavCoal_Affected <- 
  round((ZW_M_TaloNavCoal_1/ZW_M_TaloNavCoal_Total)*100, 2)

# Female N and frequency
ZW_F_TaloNavCoal <- ZW_F %>% count(TaloNavCoal)
ZW_F_TaloNavCoal <- 
  ZW_F_TaloNavCoal[(ZW_F_TaloNavCoal$TaloNavCoal != 6 & ZW_F_TaloNavCoal$TaloNavCoal != 9),]
ZW_F_TaloNavCoal_0 <- ZW_F_TaloNavCoal[1,2]
ZW_F_TaloNavCoal_0[is.na(ZW_F_TaloNavCoal_0)] <-0
ZW_F_TaloNavCoal_1 <- ZW_F_TaloNavCoal[2,2]
ZW_F_TaloNavCoal_1[is.na(ZW_F_TaloNavCoal_1)] <-0
ZW_F_TaloNavCoal_Total <- ZW_F_TaloNavCoal_0 + ZW_F_TaloNavCoal_1
ZW_F_TaloNavCoal_Affected <- 
  round((ZW_F_TaloNavCoal_1/ZW_F_TaloNavCoal_Total)*100, 2)

# overview of Taloulated frequencies for table (see below) #
ZW_TaloNavCoal_row <- c("TaloNavCoal", ZW_TaloNavCoal_Total, ZW_TaloNavCoal_1, 
                        ZW_TaloNavCoal_Affected, ZW_M_TaloNavCoal_Total, 
                        ZW_M_TaloNavCoal_1, ZW_M_TaloNavCoal_Affected, 
                        ZW_F_TaloNavCoal_Total, ZW_F_TaloNavCoal_1, 
                        ZW_F_TaloNavCoal_Affected)

## CF2CF3Coal ##
# Overall N and frequency #
ZW_CF2CF3Coal <- ZW_wo_probables %>% count(CF2CF3Coal)
ZW_CF2CF3Coal <- 
  ZW_CF2CF3Coal[(ZW_CF2CF3Coal$CF2CF3Coal != 6 & ZW_CF2CF3Coal$CF2CF3Coal != 9),]
ZW_CF2CF3Coal_0 <- ZW_CF2CF3Coal[1,2]
ZW_CF2CF3Coal_0[is.na(ZW_CF2CF3Coal_0)] <-0
ZW_CF2CF3Coal_1 <- ZW_CF2CF3Coal[2,2]
ZW_CF2CF3Coal_1[is.na(ZW_CF2CF3Coal_1)] <-0
ZW_CF2CF3Coal_Total <- ZW_CF2CF3Coal_0 + ZW_CF2CF3Coal_1
ZW_CF2CF3Coal_Affected <- round((ZW_CF2CF3Coal_1/ZW_CF2CF3Coal_Total)*100, 2)

# Male N and frequency
ZW_M_CF2CF3Coal <- ZW_M %>% count(CF2CF3Coal)
ZW_M_CF2CF3Coal <- 
  ZW_M_CF2CF3Coal[(ZW_M_CF2CF3Coal$CF2CF3Coal != 6 & ZW_M_CF2CF3Coal$CF2CF3Coal != 9),]
ZW_M_CF2CF3Coal_0 <- ZW_M_CF2CF3Coal[1,2]
ZW_M_CF2CF3Coal_0[is.na(ZW_M_CF2CF3Coal_0)] <-0
ZW_M_CF2CF3Coal_1 <- ZW_M_CF2CF3Coal[2,2]
ZW_M_CF2CF3Coal_1[is.na(ZW_M_CF2CF3Coal_1)] <-0
ZW_M_CF2CF3Coal_Total <- ZW_M_CF2CF3Coal_0 + ZW_M_CF2CF3Coal_1
ZW_M_CF2CF3Coal_Affected <- 
  round((ZW_M_CF2CF3Coal_1/ZW_M_CF2CF3Coal_Total)*100, 2)

# Female N and frequency
ZW_F_CF2CF3Coal <- ZW_F %>% count(CF2CF3Coal)
ZW_F_CF2CF3Coal <- 
  ZW_F_CF2CF3Coal[(ZW_F_CF2CF3Coal$CF2CF3Coal != 6 & ZW_F_CF2CF3Coal$CF2CF3Coal != 9),]
ZW_F_CF2CF3Coal_0 <- ZW_F_CF2CF3Coal[1,2]
ZW_F_CF2CF3Coal_0[is.na(ZW_F_CF2CF3Coal_0)] <-0
ZW_F_CF2CF3Coal_1 <- ZW_F_CF2CF3Coal[2,2]
ZW_F_CF2CF3Coal_1[is.na(ZW_F_CF2CF3Coal_1)] <-0
ZW_F_CF2CF3Coal_Total <- ZW_F_CF2CF3Coal_0 + ZW_F_CF2CF3Coal_1
ZW_F_CF2CF3Coal_Affected <- 
  round((ZW_F_CF2CF3Coal_1/ZW_F_CF2CF3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_CF2CF3Coal_row <- c("CF2CF3Coal", ZW_CF2CF3Coal_Total, ZW_CF2CF3Coal_1, 
                       ZW_CF2CF3Coal_Affected, ZW_M_CF2CF3Coal_Total, 
                       ZW_M_CF2CF3Coal_1, ZW_M_CF2CF3Coal_Affected, 
                       ZW_F_CF2CF3Coal_Total, ZW_F_CF2CF3Coal_1, 
                       ZW_F_CF2CF3Coal_Affected)

## CF3MT3Coal ##
# Overall N and frequency #
ZW_CF3MT3Coal <- ZW_wo_probables %>% count(CF3MT3Coal)
ZW_CF3MT3Coal <- 
  ZW_CF3MT3Coal[(ZW_CF3MT3Coal$CF3MT3Coal != 6 & ZW_CF3MT3Coal$CF3MT3Coal != 9),]
ZW_CF3MT3Coal_0 <- ZW_CF3MT3Coal[1,2]
ZW_CF3MT3Coal_0[is.na(ZW_CF3MT3Coal_0)] <-0
ZW_CF3MT3Coal_1 <- ZW_CF3MT3Coal[2,2]
ZW_CF3MT3Coal_1[is.na(ZW_CF3MT3Coal_1)] <-0
ZW_CF3MT3Coal_Total <- ZW_CF3MT3Coal_0 + ZW_CF3MT3Coal_1
ZW_CF3MT3Coal_Affected <- round((ZW_CF3MT3Coal_1/ZW_CF3MT3Coal_Total)*100, 2)

# Male N and frequency
ZW_M_CF3MT3Coal <- ZW_M %>% count(CF3MT3Coal)
ZW_M_CF3MT3Coal <- 
  ZW_M_CF3MT3Coal[(ZW_M_CF3MT3Coal$CF3MT3Coal != 6 & ZW_M_CF3MT3Coal$CF3MT3Coal != 9),]
ZW_M_CF3MT3Coal_0 <- ZW_M_CF3MT3Coal[1,2]
ZW_M_CF3MT3Coal_0[is.na(ZW_M_CF3MT3Coal_0)] <-0
ZW_M_CF3MT3Coal_1 <- ZW_M_CF3MT3Coal[2,2]
ZW_M_CF3MT3Coal_1[is.na(ZW_M_CF3MT3Coal_1)] <-0
ZW_M_CF3MT3Coal_Total <- ZW_M_CF3MT3Coal_0 + ZW_M_CF3MT3Coal_1
ZW_M_CF3MT3Coal_Affected <- 
  round((ZW_M_CF3MT3Coal_1/ZW_M_CF3MT3Coal_Total)*100, 2)

# Female N and frequency
ZW_F_CF3MT3Coal <- ZW_F %>% count(CF3MT3Coal)
ZW_F_CF3MT3Coal <- 
  ZW_F_CF3MT3Coal[(ZW_F_CF3MT3Coal$CF3MT3Coal != 6 & ZW_F_CF3MT3Coal$CF3MT3Coal != 9),]
ZW_F_CF3MT3Coal_0 <- ZW_F_CF3MT3Coal[1,2]
ZW_F_CF3MT3Coal_0[is.na(ZW_F_CF3MT3Coal_0)] <-0
ZW_F_CF3MT3Coal_1 <- ZW_F_CF3MT3Coal[2,2]
ZW_F_CF3MT3Coal_1[is.na(ZW_F_CF3MT3Coal_1)] <-0
ZW_F_CF3MT3Coal_Total <- ZW_F_CF3MT3Coal_0 + ZW_F_CF3MT3Coal_1
ZW_F_CF3MT3Coal_Affected <- 
  round((ZW_F_CF3MT3Coal_1/ZW_F_CF3MT3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_CF3MT3Coal_row <- c("CF3MT3Coal", ZW_CF3MT3Coal_Total, ZW_CF3MT3Coal_1, 
                       ZW_CF3MT3Coal_Affected, ZW_M_CF3MT3Coal_Total, 
                       ZW_M_CF3MT3Coal_1, ZW_M_CF3MT3Coal_Affected, 
                       ZW_F_CF3MT3Coal_Total, ZW_F_CF3MT3Coal_1, 
                       ZW_F_CF3MT3Coal_Affected)
## OsIntermet ##
# Overall N and frequency #
ZW_OsIntermet <- ZW_wo_probables %>% count(OsIntermet)
ZW_OsIntermet <- 
  ZW_OsIntermet[(ZW_OsIntermet$OsIntermet != 6 & ZW_OsIntermet$OsIntermet != 9),]
ZW_OsIntermet_0 <- ZW_OsIntermet[1,2]
ZW_OsIntermet_0[is.na(ZW_OsIntermet_0)] <-0
ZW_OsIntermet_1 <- ZW_OsIntermet[2,2]
ZW_OsIntermet_1[is.na(ZW_OsIntermet_1)] <-0
ZW_OsIntermet_Total <- ZW_OsIntermet_0 + ZW_OsIntermet_1
ZW_OsIntermet_Affected <- round((ZW_OsIntermet_1/ZW_OsIntermet_Total)*100, 2)

# Male N and frequency
ZW_M_OsIntermet <- ZW_M %>% count(OsIntermet)
ZW_M_OsIntermet <- 
  ZW_M_OsIntermet[(ZW_M_OsIntermet$OsIntermet != 6 & ZW_M_OsIntermet$OsIntermet != 9),]
ZW_M_OsIntermet_0 <- ZW_M_OsIntermet[1,2]
ZW_M_OsIntermet_0[is.na(ZW_M_OsIntermet_0)] <-0
ZW_M_OsIntermet_1 <- ZW_M_OsIntermet[2,2]
ZW_M_OsIntermet_1[is.na(ZW_M_OsIntermet_1)] <-0
ZW_M_OsIntermet_Total <- ZW_M_OsIntermet_0 + ZW_M_OsIntermet_1
ZW_M_OsIntermet_Affected <- 
  round((ZW_M_OsIntermet_1/ZW_M_OsIntermet_Total)*100, 2)

# Female N and frequency
ZW_F_OsIntermet <- ZW_F %>% count(OsIntermet)
ZW_F_OsIntermet <- 
  ZW_F_OsIntermet[(ZW_F_OsIntermet$OsIntermet != 6 & ZW_F_OsIntermet$OsIntermet != 9),]
ZW_F_OsIntermet_0 <- ZW_F_OsIntermet[1,2]
ZW_F_OsIntermet_0[is.na(ZW_F_OsIntermet_0)] <-0
ZW_F_OsIntermet_1 <- ZW_F_OsIntermet[2,2]
ZW_F_OsIntermet_1[is.na(ZW_F_OsIntermet_1)] <-0
ZW_F_OsIntermet_Total <- ZW_F_OsIntermet_0 + ZW_F_OsIntermet_1
ZW_F_OsIntermet_Affected <- 
  round((ZW_F_OsIntermet_1/ZW_F_OsIntermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_OsIntermet_row <- c("OsIntermet", ZW_OsIntermet_Total, ZW_OsIntermet_1, 
                       ZW_OsIntermet_Affected, ZW_M_OsIntermet_Total, 
                       ZW_M_OsIntermet_1, ZW_M_OsIntermet_Affected, 
                       ZW_F_OsIntermet_Total, ZW_F_OsIntermet_1, 
                       ZW_F_OsIntermet_Affected)
# OsIntermet unclustered: CF1Intermet #
# Overall N and frequency #
ZW_CF1Intermet <- ZW_wo_probables %>% count(CF1Intermet)
ZW_CF1Intermet <- 
  ZW_CF1Intermet[(ZW_CF1Intermet$CF1Intermet != 6 & ZW_CF1Intermet$CF1Intermet != 9),]
ZW_CF1Intermet_0 <- ZW_CF1Intermet[1,2]
ZW_CF1Intermet_0[is.na(ZW_CF1Intermet_0)] <-0
ZW_CF1Intermet_1 <- ZW_CF1Intermet[2,2]
ZW_CF1Intermet_1[is.na(ZW_CF1Intermet_1)] <-0
ZW_CF1Intermet_Total <- na.omit(ZW_CF1Intermet_0) + na.omit(ZW_CF1Intermet_1)
ZW_CF1Intermet_Affected <- round((ZW_CF1Intermet_1/ZW_CF1Intermet_Total)*100, 2)

# Male N and frequency
ZW_M_CF1Intermet <- ZW_M %>% count(CF1Intermet)
ZW_M_CF1Intermet <- 
  ZW_M_CF1Intermet[(ZW_M_CF1Intermet$CF1Intermet != 6 & ZW_M_CF1Intermet$CF1Intermet != 9),]
ZW_M_CF1Intermet_0 <- ZW_M_CF1Intermet[1,2]
ZW_M_CF1Intermet_0[is.na(ZW_M_CF1Intermet_0)] <-0
ZW_M_CF1Intermet_1 <- ZW_M_CF1Intermet[2,2]
ZW_M_CF1Intermet_1[is.na(ZW_M_CF1Intermet_1)] <-0
ZW_M_CF1Intermet_Total <- ZW_M_CF1Intermet_0 + ZW_M_CF1Intermet_1
ZW_M_CF1Intermet_Affected <- 
  round((ZW_M_CF1Intermet_1/ZW_M_CF1Intermet_Total)*100, 2)

# Female N and frequency
ZW_F_CF1Intermet <- ZW_F %>% count(CF1Intermet)
ZW_F_CF1Intermet <- 
  ZW_F_CF1Intermet[(ZW_F_CF1Intermet$CF1Intermet != 6 & ZW_F_CF1Intermet$CF1Intermet != 9),]
ZW_F_CF1Intermet_0 <- ZW_F_CF1Intermet[1,2]
ZW_F_CF1Intermet_0[is.na(ZW_F_CF1Intermet_0)] <-0
ZW_F_CF1Intermet_1 <- ZW_F_CF1Intermet[2,2]
ZW_F_CF1Intermet_1[is.na(ZW_F_CF1Intermet_1)] <-0
ZW_F_CF1Intermet_Total <- ZW_F_CF1Intermet_0 + ZW_F_CF1Intermet_1
ZW_F_CF1Intermet_Affected <- 
  round((ZW_F_CF1Intermet_1/ZW_F_CF1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_CF1Intermet_row <- c("CF1Intermet", ZW_CF1Intermet_Total, ZW_CF1Intermet_1, 
                        ZW_CF1Intermet_Affected, ZW_M_CF1Intermet_Total, 
                        ZW_M_CF1Intermet_1, ZW_M_CF1Intermet_Affected, 
                        ZW_F_CF1Intermet_Total, ZW_F_CF1Intermet_1, 
                        ZW_F_CF1Intermet_Affected)

# OsIntermet unclustered: MT1Intermet #
# Overall N and frequency #
ZW_MT1Intermet <- ZW_wo_probables %>% count(MT1Intermet)
ZW_MT1Intermet <- 
  ZW_MT1Intermet[(ZW_MT1Intermet$MT1Intermet != 6 & ZW_MT1Intermet$MT1Intermet != 9),]
ZW_MT1Intermet_0 <- ZW_MT1Intermet[1,2]
ZW_MT1Intermet_0[is.na(ZW_MT1Intermet_0)] <-0
ZW_MT1Intermet_1 <- ZW_MT1Intermet[2,2]
ZW_MT1Intermet_1[is.na(ZW_MT1Intermet_1)] <-0
ZW_MT1Intermet_Total <- na.omit(ZW_MT1Intermet_0) + na.omit(ZW_MT1Intermet_1)
ZW_MT1Intermet_Affected <- round((ZW_MT1Intermet_1/ZW_MT1Intermet_Total)*100, 2)

# Male N and frequency
ZW_M_MT1Intermet <- ZW_M %>% count(MT1Intermet)
ZW_M_MT1Intermet <- 
  ZW_M_MT1Intermet[(ZW_M_MT1Intermet$MT1Intermet != 6 & ZW_M_MT1Intermet$MT1Intermet != 9),]
ZW_M_MT1Intermet_0 <- ZW_M_MT1Intermet[1,2]
ZW_M_MT1Intermet_0[is.na(ZW_M_MT1Intermet_0)] <-0
ZW_M_MT1Intermet_1 <- ZW_M_MT1Intermet[2,2]
ZW_M_MT1Intermet_1[is.na(ZW_M_MT1Intermet_1)] <-0
ZW_M_MT1Intermet_Total <- ZW_M_MT1Intermet_0 + ZW_M_MT1Intermet_1
ZW_M_MT1Intermet_Affected <- 
  round((ZW_M_MT1Intermet_1/ZW_M_MT1Intermet_Total)*100, 2)

# Female N and frequency
ZW_F_MT1Intermet <- ZW_F %>% count(MT1Intermet)
ZW_F_MT1Intermet <- 
  ZW_F_MT1Intermet[(ZW_F_MT1Intermet$MT1Intermet != 6 & ZW_F_MT1Intermet$MT1Intermet != 9),]
ZW_F_MT1Intermet_0 <- ZW_F_MT1Intermet[1,2]
ZW_F_MT1Intermet_0[is.na(ZW_F_MT1Intermet_0)] <-0
ZW_F_MT1Intermet_1 <- ZW_F_MT1Intermet[2,2]
ZW_F_MT1Intermet_1[is.na(ZW_F_MT1Intermet_1)] <-0
ZW_F_MT1Intermet_Total <- ZW_F_MT1Intermet_0 + ZW_F_MT1Intermet_1
ZW_F_MT1Intermet_Affected <- 
  round((ZW_F_MT1Intermet_1/ZW_F_MT1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_MT1Intermet_row <- c("MT1Intermet", ZW_MT1Intermet_Total, ZW_MT1Intermet_1, 
                        ZW_MT1Intermet_Affected, ZW_M_MT1Intermet_Total, 
                        ZW_M_MT1Intermet_1, ZW_M_MT1Intermet_Affected, 
                        ZW_F_MT1Intermet_Total, ZW_F_MT1Intermet_1, 
                        ZW_F_MT1Intermet_Affected)

# OsIntermet unclustered: MT2Intermet #
# Overall N and frequency #
ZW_MT2Intermet <- ZW_wo_probables %>% count(MT2Intermet)
ZW_MT2Intermet <- 
  ZW_MT2Intermet[(ZW_MT2Intermet$MT2Intermet != 6 & ZW_MT2Intermet$MT2Intermet != 9),]
ZW_MT2Intermet_0 <- ZW_MT2Intermet[1,2]
ZW_MT2Intermet_0[is.na(ZW_MT2Intermet_0)] <-0
ZW_MT2Intermet_1 <- ZW_MT2Intermet[2,2]
ZW_MT2Intermet_1[is.na(ZW_MT2Intermet_1)] <-0
ZW_MT2Intermet_Total <- na.omit(ZW_MT2Intermet_0) + na.omit(ZW_MT2Intermet_1)
ZW_MT2Intermet_Affected <- round((ZW_MT2Intermet_1/ZW_MT2Intermet_Total)*100, 2)

# Male N and frequency
ZW_M_MT2Intermet <- ZW_M %>% count(MT2Intermet)
ZW_M_MT2Intermet <- 
  ZW_M_MT2Intermet[(ZW_M_MT2Intermet$MT2Intermet != 6 & ZW_M_MT2Intermet$MT2Intermet != 9),]
ZW_M_MT2Intermet_0 <- ZW_M_MT2Intermet[1,2]
ZW_M_MT2Intermet_0[is.na(ZW_M_MT2Intermet_0)] <-0
ZW_M_MT2Intermet_1 <- ZW_M_MT2Intermet[2,2]
ZW_M_MT2Intermet_1[is.na(ZW_M_MT2Intermet_1)] <-0
ZW_M_MT2Intermet_Total <- ZW_M_MT2Intermet_0 + ZW_M_MT2Intermet_1
ZW_M_MT2Intermet_Affected <- 
  round((ZW_M_MT2Intermet_1/ZW_M_MT2Intermet_Total)*100, 2)

# Female N and frequency
ZW_F_MT2Intermet <- ZW_F %>% count(MT2Intermet)
ZW_F_MT2Intermet <- 
  ZW_F_MT2Intermet[(ZW_F_MT2Intermet$MT2Intermet != 6 & ZW_F_MT2Intermet$MT2Intermet != 9),]
ZW_F_MT2Intermet_0 <- ZW_F_MT2Intermet[1,2]
ZW_F_MT2Intermet_0[is.na(ZW_F_MT2Intermet_0)] <-0
ZW_F_MT2Intermet_1 <- ZW_F_MT2Intermet[2,2]
ZW_F_MT2Intermet_1[is.na(ZW_F_MT2Intermet_1)] <-0
ZW_F_MT2Intermet_Total <- ZW_F_MT2Intermet_0 + ZW_F_MT2Intermet_1
ZW_F_MT2Intermet_Affected <- 
  round((ZW_F_MT2Intermet_1/ZW_F_MT2Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_MT2Intermet_row <- c("MT2Intermet", ZW_MT2Intermet_Total, ZW_MT2Intermet_1, 
                        ZW_MT2Intermet_Affected, ZW_M_MT2Intermet_Total, 
                        ZW_M_MT2Intermet_1, ZW_M_MT2Intermet_Affected, 
                        ZW_F_MT2Intermet_Total, ZW_F_MT2Intermet_1, 
                        ZW_F_MT2Intermet_Affected)


## AnkleCoal ##
# Overall N and frequency #
ZW_AnkleCoal <- ZW_wo_probables %>% count(AnkleCoal)
ZW_AnkleCoal <- 
  ZW_AnkleCoal[(ZW_AnkleCoal$AnkleCoal != 6 & ZW_AnkleCoal$AnkleCoal != 9),]
ZW_AnkleCoal_0 <- ZW_AnkleCoal[1,2]; ZW_AnkleCoal_0[is.na(ZW_AnkleCoal_0)] <-0
ZW_AnkleCoal_1 <- ZW_AnkleCoal[2,2]; ZW_AnkleCoal_1[is.na(ZW_AnkleCoal_1)] <-0
ZW_AnkleCoal_Total <- na.omit(ZW_AnkleCoal_0) + na.omit(ZW_AnkleCoal_1)
ZW_AnkleCoal_Affected <- round((ZW_AnkleCoal_1/ZW_AnkleCoal_Total)*100, 2)

# Male N and frequency
ZW_M_AnkleCoal <- ZW_M %>% count(AnkleCoal)
ZW_M_AnkleCoal <- 
  ZW_M_AnkleCoal[(ZW_M_AnkleCoal$AnkleCoal != 6 & ZW_M_AnkleCoal$AnkleCoal != 9),]
ZW_M_AnkleCoal_0 <- ZW_M_AnkleCoal[1,2]
ZW_M_AnkleCoal_0[is.na(ZW_M_AnkleCoal_0)] <-0
ZW_M_AnkleCoal_1 <- ZW_M_AnkleCoal[2,2]
ZW_M_AnkleCoal_1[is.na(ZW_M_AnkleCoal_1)] <-0
ZW_M_AnkleCoal_Total <- ZW_M_AnkleCoal_0 + ZW_M_AnkleCoal_1
ZW_M_AnkleCoal_Affected <- round((ZW_M_AnkleCoal_1/ZW_M_AnkleCoal_Total)*100, 2)

# Female N and frequency
ZW_F_AnkleCoal <- ZW_F %>% count(AnkleCoal)
ZW_F_AnkleCoal <- 
  ZW_F_AnkleCoal[(ZW_F_AnkleCoal$AnkleCoal != 6 & ZW_F_AnkleCoal$AnkleCoal != 9),]
ZW_F_AnkleCoal_0 <- ZW_F_AnkleCoal[1,2]
ZW_F_AnkleCoal_0[is.na(ZW_F_AnkleCoal_0)] <-0
ZW_F_AnkleCoal_1 <- ZW_F_AnkleCoal[2,2]
ZW_F_AnkleCoal_1[is.na(ZW_F_AnkleCoal_1)] <-0
ZW_F_AnkleCoal_Total <- ZW_F_AnkleCoal_0 + ZW_F_AnkleCoal_1
ZW_F_AnkleCoal_Affected <- round((ZW_F_AnkleCoal_1/ZW_F_AnkleCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_AnkleCoal_row <- c("AnkleCoal", ZW_AnkleCoal_Total, ZW_AnkleCoal_1,
                      ZW_AnkleCoal_Affected, ZW_M_AnkleCoal_Total, 
                      ZW_M_AnkleCoal_1, ZW_M_AnkleCoal_Affected, ZW_F_AnkleCoal_Total, 
                      ZW_F_AnkleCoal_1, ZW_F_AnkleCoal_Affected)

# AnkleCoal unclustered: CalcNavCoal #
# Overall N and frequency #
ZW_CalcNavCoal <- ZW_wo_probables %>% count(CalcNavCoal)
ZW_CalcNavCoal <- 
  ZW_CalcNavCoal[(ZW_CalcNavCoal$CalcNavCoal != 6 & ZW_CalcNavCoal$CalcNavCoal != 9),]
ZW_CalcNavCoal_0 <- ZW_CalcNavCoal[1,2]
ZW_CalcNavCoal_0[is.na(ZW_CalcNavCoal_0)] <-0
ZW_CalcNavCoal_1 <- ZW_CalcNavCoal[2,2]
ZW_CalcNavCoal_1[is.na(ZW_CalcNavCoal_1)] <-0
ZW_CalcNavCoal_Total <- na.omit(ZW_CalcNavCoal_0) + na.omit(ZW_CalcNavCoal_1)
ZW_CalcNavCoal_Affected <- round((ZW_CalcNavCoal_1/ZW_CalcNavCoal_Total)*100, 2)

# Male N and frequency
ZW_M_CalcNavCoal <- ZW_M %>% count(CalcNavCoal)
ZW_M_CalcNavCoal <- 
  ZW_M_CalcNavCoal[(ZW_M_CalcNavCoal$CalcNavCoal != 6 & ZW_M_CalcNavCoal$CalcNavCoal != 9),]
ZW_M_CalcNavCoal_0 <- ZW_M_CalcNavCoal[1,2]
ZW_M_CalcNavCoal_0[is.na(ZW_M_CalcNavCoal_0)] <-0
ZW_M_CalcNavCoal_1 <- ZW_M_CalcNavCoal[2,2]
ZW_M_CalcNavCoal_1[is.na(ZW_M_CalcNavCoal_1)] <-0
ZW_M_CalcNavCoal_Total <- ZW_M_CalcNavCoal_0 + ZW_M_CalcNavCoal_1
ZW_M_CalcNavCoal_Affected <- 
  round((ZW_M_CalcNavCoal_1/ZW_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
ZW_F_CalcNavCoal <- ZW_F %>% count(CalcNavCoal)
ZW_F_CalcNavCoal <- 
  ZW_F_CalcNavCoal[(ZW_F_CalcNavCoal$CalcNavCoal != 6 & ZW_F_CalcNavCoal$CalcNavCoal != 9),]
ZW_F_CalcNavCoal_0 <- ZW_F_CalcNavCoal[1,2]
ZW_F_CalcNavCoal_0[is.na(ZW_F_CalcNavCoal_0)] <-0
ZW_F_CalcNavCoal_1 <- ZW_F_CalcNavCoal[2,2]
ZW_F_CalcNavCoal_1[is.na(ZW_F_CalcNavCoal_1)] <-0
ZW_F_CalcNavCoal_Total <- ZW_F_CalcNavCoal_0 + ZW_F_CalcNavCoal_1
ZW_F_CalcNavCoal_Affected <- 
  round((ZW_F_CalcNavCoal_1/ZW_F_CalcNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_CalcNavCoal_row <- c("CalcNavCoal", ZW_CalcNavCoal_Total, ZW_CalcNavCoal_1, 
                        ZW_CalcNavCoal_Affected, ZW_M_CalcNavCoal_Total, 
                        ZW_M_CalcNavCoal_1, ZW_M_CalcNavCoal_Affected, 
                        ZW_F_CalcNavCoal_Total, ZW_F_CalcNavCoal_1, 
                        ZW_F_CalcNavCoal_Affected)

# AnkleCoal unclustered: TaloCalcCoal #
# Overall N and frequency #
ZW_TaloCalcCoal <- ZW_wo_probables %>% count(TaloCalcCoal)
ZW_TaloCalcCoal <- 
  ZW_TaloCalcCoal[(ZW_TaloCalcCoal$TaloCalcCoal != 6 & ZW_TaloCalcCoal$TaloCalcCoal != 9),]
ZW_TaloCalcCoal_0 <- ZW_TaloCalcCoal[1,2]
ZW_TaloCalcCoal_0[is.na(ZW_TaloCalcCoal_0)] <-0
ZW_TaloCalcCoal_1 <- ZW_TaloCalcCoal[2,2]
ZW_TaloCalcCoal_1[is.na(ZW_TaloCalcCoal_1)] <-0
ZW_TaloCalcCoal_Total <- na.omit(ZW_TaloCalcCoal_0) + na.omit(ZW_TaloCalcCoal_1)
ZW_TaloCalcCoal_Affected <- 
  round((ZW_TaloCalcCoal_1/ZW_TaloCalcCoal_Total)*100, 2)

# Male N and frequency
ZW_M_TaloCalcCoal <- ZW_M %>% count(TaloCalcCoal)
ZW_M_TaloCalcCoal <- 
  ZW_M_TaloCalcCoal[(ZW_M_TaloCalcCoal$TaloCalcCoal != 6 & ZW_M_TaloCalcCoal$TaloCalcCoal != 9),]
ZW_M_TaloCalcCoal_0 <- ZW_M_TaloCalcCoal[1,2]
ZW_M_TaloCalcCoal_0[is.na(ZW_M_TaloCalcCoal_0)] <-0
ZW_M_TaloCalcCoal_1 <- ZW_M_TaloCalcCoal[2,2]
ZW_M_TaloCalcCoal_1[is.na(ZW_M_TaloCalcCoal_1)] <-0
ZW_M_TaloCalcCoal_Total <- ZW_M_TaloCalcCoal_0 + ZW_M_TaloCalcCoal_1
ZW_M_TaloCalcCoal_Affected <- 
  round((ZW_M_TaloCalcCoal_1/ZW_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
ZW_F_TaloCalcCoal <- ZW_F %>% count(TaloCalcCoal)
ZW_F_TaloCalcCoal <- 
  ZW_F_TaloCalcCoal[(ZW_F_TaloCalcCoal$TaloCalcCoal != 6 & ZW_F_TaloCalcCoal$TaloCalcCoal != 9),]
ZW_F_TaloCalcCoal_0 <- ZW_F_TaloCalcCoal[1,2]
ZW_F_TaloCalcCoal_0[is.na(ZW_F_TaloCalcCoal_0)] <-0
ZW_F_TaloCalcCoal_1 <- ZW_F_TaloCalcCoal[2,2]
ZW_F_TaloCalcCoal_1[is.na(ZW_F_TaloCalcCoal_1)] <-0
ZW_F_TaloCalcCoal_Total <- ZW_F_TaloCalcCoal_0 + ZW_F_TaloCalcCoal_1
ZW_F_TaloCalcCoal_Affected <- 
  round((ZW_F_TaloCalcCoal_1/ZW_F_TaloCalcCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
ZW_TaloCalcCoal_row <- c("TaloCalcCoal", ZW_TaloCalcCoal_Total, ZW_TaloCalcCoal_1, 
                         ZW_TaloCalcCoal_Affected, ZW_M_TaloCalcCoal_Total, 
                         ZW_M_TaloCalcCoal_1, ZW_M_TaloCalcCoal_Affected, 
                         ZW_F_TaloCalcCoal_Total, ZW_F_TaloCalcCoal_1, 
                         ZW_F_TaloCalcCoal_Affected)


## Table with overview of trait frequencies for ZW ##
ZW_freq <- data.frame(rbind(ZW_AccessNav_row, ZW_BrachyD_row, ZW_BrachyMT1_row,
                            ZW_BrachyMT4_row, ZW_BrachyPP1_row, ZW_CalcCubCoal_row,
                            ZW_TaloNavCoal_row, ZW_CF2CF3Coal_row, ZW_CF3MT3Coal_row,
                            ZW_OsIntermet_row, ZW_CF1Intermet_row, ZW_MT1Intermet_row,
                            ZW_MT2Intermet_row, ZW_AnkleCoal_row, ZW_CalcNavCoal_row, 
                            ZW_TaloCalcCoal_row))
names(ZW_freq) <- c("Anomaly", "Overall N (ZW)", "Overall affected (ZW)", 
                    "Overall % (ZW)", "Male N (ZW)", "Male affected (ZW)",
                    "Male % (ZW)", "Female N (ZW)", "Female affected (ZW)", 
                    "Female % (ZW)")


### The reference sample ###


## AccessNav ##
# Overall N and frequency #
Ref_Overall_AccessNav <- Ref_Overall %>% count(AccessNav)
Ref_Overall_AccessNav <- 
  Ref_Overall_AccessNav[(Ref_Overall_AccessNav$AccessNav != 6 & Ref_Overall_AccessNav$AccessNav != 9),]
Ref_Overall_AccessNav_0 <- Ref_Overall_AccessNav[1,2]
Ref_Overall_AccessNav_0[is.na(Ref_Overall_AccessNav_0)] <-0
Ref_Overall_AccessNav_1 <- Ref_Overall_AccessNav[2,2]
Ref_Overall_AccessNav_1[is.na(Ref_Overall_AccessNav_1)] <-0
Ref_Overall_AccessNav_Total <- Ref_Overall_AccessNav_0 + Ref_Overall_AccessNav_1
Ref_Overall_AccessNav_Affected <- 
  round((Ref_Overall_AccessNav_1/Ref_Overall_AccessNav_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_AccessNav <- Ref_Overall_M %>% count(AccessNav)
Ref_Overall_M_AccessNav <- 
  Ref_Overall_M_AccessNav[(Ref_Overall_M_AccessNav$AccessNav != 6 & Ref_Overall_M_AccessNav$AccessNav != 9),]
Ref_Overall_M_AccessNav_0 <- Ref_Overall_M_AccessNav[1,2]
Ref_Overall_M_AccessNav_0[is.na(Ref_Overall_M_AccessNav_0)] <-0
Ref_Overall_M_AccessNav_1 <- Ref_Overall_M_AccessNav[2,2]
Ref_Overall_M_AccessNav_1[is.na(Ref_Overall_M_AccessNav_1)] <-0
Ref_Overall_M_AccessNav_Total <- 
  Ref_Overall_M_AccessNav_0 + Ref_Overall_M_AccessNav_1
Ref_Overall_M_AccessNav_Affected <- 
  round((Ref_Overall_M_AccessNav_1/Ref_Overall_M_AccessNav_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_AccessNav <- Ref_Overall_F %>% count(AccessNav)
Ref_Overall_F_AccessNav <- 
  Ref_Overall_F_AccessNav[(Ref_Overall_F_AccessNav$AccessNav != 6 & Ref_Overall_F_AccessNav$AccessNav != 9),]
Ref_Overall_F_AccessNav_0 <- Ref_Overall_F_AccessNav[1,2]
Ref_Overall_F_AccessNav_0[is.na(Ref_Overall_F_AccessNav_0)] <-0
Ref_Overall_F_AccessNav_1 <- Ref_Overall_F_AccessNav[2,2]
Ref_Overall_F_AccessNav_1[is.na(Ref_Overall_F_AccessNav_1)] <-0
Ref_Overall_F_AccessNav_Total <- 
  Ref_Overall_F_AccessNav_0 + Ref_Overall_F_AccessNav_1
Ref_Overall_F_AccessNav_Affected <- 
  round((Ref_Overall_F_AccessNav_1/Ref_Overall_F_AccessNav_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_AccessNav_row <- 
  c("AccessNav", Ref_Overall_AccessNav_Total, Ref_Overall_AccessNav_1, 
     Ref_Overall_AccessNav_Affected, Ref_Overall_M_AccessNav_Total, 
     Ref_Overall_M_AccessNav_1, Ref_Overall_M_AccessNav_Affected, 
     Ref_Overall_F_AccessNav_Total, Ref_Overall_F_AccessNav_1, 
     Ref_Overall_F_AccessNav_Affected)

## BrachyD ##
# Overall N and frequency #
Ref_Overall_BrachyD <- Ref_Overall %>% count(BrachyD)
Ref_Overall_BrachyD <- Ref_Overall_BrachyD[(Ref_Overall_BrachyD$BrachyD != 6 & Ref_Overall_BrachyD$BrachyD != 9),]
Ref_Overall_BrachyD_0 <- Ref_Overall_BrachyD[1,2]
Ref_Overall_BrachyD_0[is.na(Ref_Overall_BrachyD_0)] <-0
Ref_Overall_BrachyD_1 <- Ref_Overall_BrachyD[2,2]
Ref_Overall_BrachyD_1[is.na(Ref_Overall_BrachyD_1)] <-0
Ref_Overall_BrachyD_Total <- Ref_Overall_BrachyD_0 + Ref_Overall_BrachyD_1
Ref_Overall_BrachyD_Affected <- 
  round((Ref_Overall_BrachyD_1/Ref_Overall_BrachyD_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_BrachyD <- Ref_Overall_M %>% count(BrachyD)
Ref_Overall_M_BrachyD <- 
  Ref_Overall_M_BrachyD[(Ref_Overall_M_BrachyD$BrachyD != 6 & Ref_Overall_M_BrachyD$BrachyD != 9),]
Ref_Overall_M_BrachyD_0 <- Ref_Overall_M_BrachyD[1,2]
Ref_Overall_M_BrachyD_0[is.na(Ref_Overall_M_BrachyD_0)] <-0
Ref_Overall_M_BrachyD_1 <- Ref_Overall_M_BrachyD[2,2]
Ref_Overall_M_BrachyD_1[is.na(Ref_Overall_M_BrachyD_1)] <-0
Ref_Overall_M_BrachyD_Total <- 
  Ref_Overall_M_BrachyD_0 + Ref_Overall_M_BrachyD_1
Ref_Overall_M_BrachyD_Affected <- 
  round((Ref_Overall_M_BrachyD_1/Ref_Overall_M_BrachyD_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_BrachyD <- Ref_Overall_F %>% count(BrachyD)
Ref_Overall_F_BrachyD <- 
  Ref_Overall_F_BrachyD[(Ref_Overall_F_BrachyD$BrachyD != 6 & Ref_Overall_F_BrachyD$BrachyD != 9),]
Ref_Overall_F_BrachyD_0 <- Ref_Overall_F_BrachyD[1,2]
Ref_Overall_F_BrachyD_0[is.na(Ref_Overall_F_BrachyD_0)] <-0
Ref_Overall_F_BrachyD_1 <- Ref_Overall_F_BrachyD[2,2]
Ref_Overall_F_BrachyD_1[is.na(Ref_Overall_F_BrachyD_1)] <-0
Ref_Overall_F_BrachyD_Total <- 
  Ref_Overall_F_BrachyD_0 + Ref_Overall_F_BrachyD_1
Ref_Overall_F_BrachyD_Affected <- 
  round((Ref_Overall_F_BrachyD_1/Ref_Overall_F_BrachyD_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_BrachyD_row <- 
  c("BrachyD", Ref_Overall_BrachyD_Total, Ref_Overall_BrachyD_1, 
     Ref_Overall_BrachyD_Affected, Ref_Overall_M_BrachyD_Total, 
     Ref_Overall_M_BrachyD_1, Ref_Overall_M_BrachyD_Affected, 
     Ref_Overall_F_BrachyD_Total, Ref_Overall_F_BrachyD_1, 
     Ref_Overall_F_BrachyD_Affected)


## BrachyMT1 ##
# Overall N and frequency #
Ref_Overall_BrachyMT1 <- Ref_Overall %>% count(BrachyMT1)
Ref_Overall_BrachyMT1 <- 
  Ref_Overall_BrachyMT1[(Ref_Overall_BrachyMT1$BrachyMT1 != 6 & Ref_Overall_BrachyMT1$BrachyMT1 != 9),]
Ref_Overall_BrachyMT1_0 <- Ref_Overall_BrachyMT1[1,2]
Ref_Overall_BrachyMT1_0[is.na(Ref_Overall_BrachyMT1_0)] <-0
Ref_Overall_BrachyMT1_1 <- Ref_Overall_BrachyMT1[2,2]
Ref_Overall_BrachyMT1_1[is.na(Ref_Overall_BrachyMT1_1)] <-0
Ref_Overall_BrachyMT1_Total <- Ref_Overall_BrachyMT1_0 + Ref_Overall_BrachyMT1_1
Ref_Overall_BrachyMT1_Affected <- 
  round((Ref_Overall_BrachyMT1_1/Ref_Overall_BrachyMT1_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_BrachyMT1 <- Ref_Overall_M %>% count(BrachyMT1)
Ref_Overall_M_BrachyMT1 <- 
  Ref_Overall_M_BrachyMT1[(Ref_Overall_M_BrachyMT1$BrachyMT1 != 6 & Ref_Overall_M_BrachyMT1$BrachyMT1 != 9),]
Ref_Overall_M_BrachyMT1_0 <- Ref_Overall_M_BrachyMT1[1,2]
Ref_Overall_M_BrachyMT1_0[is.na(Ref_Overall_M_BrachyMT1_0)] <-0
Ref_Overall_M_BrachyMT1_1 <- Ref_Overall_M_BrachyMT1[2,2]
Ref_Overall_M_BrachyMT1_1[is.na(Ref_Overall_M_BrachyMT1_1)] <-0
Ref_Overall_M_BrachyMT1_Total <- 
  Ref_Overall_M_BrachyMT1_0 + Ref_Overall_M_BrachyMT1_1
Ref_Overall_M_BrachyMT1_Affected <- 
  round((Ref_Overall_M_BrachyMT1_1/Ref_Overall_M_BrachyMT1_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_BrachyMT1 <- Ref_Overall_F %>% count(BrachyMT1)
Ref_Overall_F_BrachyMT1 <- 
  Ref_Overall_F_BrachyMT1[(Ref_Overall_F_BrachyMT1$BrachyMT1 != 6 & Ref_Overall_F_BrachyMT1$BrachyMT1 != 9),]
Ref_Overall_F_BrachyMT1_0 <- Ref_Overall_F_BrachyMT1[1,2]
Ref_Overall_F_BrachyMT1_0[is.na(Ref_Overall_F_BrachyMT1_0)] <-0
Ref_Overall_F_BrachyMT1_1 <- Ref_Overall_F_BrachyMT1[2,2]
Ref_Overall_F_BrachyMT1_1[is.na(Ref_Overall_F_BrachyMT1_1)] <-0
Ref_Overall_F_BrachyMT1_Total <- 
  Ref_Overall_F_BrachyMT1_0 + Ref_Overall_F_BrachyMT1_1
Ref_Overall_F_BrachyMT1_Affected <- 
  round((Ref_Overall_F_BrachyMT1_1/Ref_Overall_F_BrachyMT1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_BrachyMT1_row <- 
  c("BrachyMT1", Ref_Overall_BrachyMT1_Total, Ref_Overall_BrachyMT1_1, 
     Ref_Overall_BrachyMT1_Affected, Ref_Overall_M_BrachyMT1_Total, 
     Ref_Overall_M_BrachyMT1_1, Ref_Overall_M_BrachyMT1_Affected, 
     Ref_Overall_F_BrachyMT1_Total, Ref_Overall_F_BrachyMT1_1, 
     Ref_Overall_F_BrachyMT1_Affected)

## BrachyMT4 ##
# Overall N and frequency #
Ref_Overall_BrachyMT4 <- Ref_Overall %>% count(BrachyMT4)
Ref_Overall_BrachyMT4 <- 
  Ref_Overall_BrachyMT4[(Ref_Overall_BrachyMT4$BrachyMT4 != 6 & Ref_Overall_BrachyMT4$BrachyMT4 != 9),]
Ref_Overall_BrachyMT4_0 <- Ref_Overall_BrachyMT4[1,2]
Ref_Overall_BrachyMT4_0[is.na(Ref_Overall_BrachyMT4_0)] <-0
Ref_Overall_BrachyMT4_1 <- Ref_Overall_BrachyMT4[2,2]
Ref_Overall_BrachyMT4_1[is.na(Ref_Overall_BrachyMT4_1)] <-0
Ref_Overall_BrachyMT4_Total <- Ref_Overall_BrachyMT4_0 + Ref_Overall_BrachyMT4_1
Ref_Overall_BrachyMT4_Affected <- 
  round((Ref_Overall_BrachyMT4_1/Ref_Overall_BrachyMT4_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_BrachyMT4 <- Ref_Overall_M %>% count(BrachyMT4)
Ref_Overall_M_BrachyMT4 <- 
  Ref_Overall_M_BrachyMT4[(Ref_Overall_M_BrachyMT4$BrachyMT4 != 6 & Ref_Overall_M_BrachyMT4$BrachyMT4 != 9),]
Ref_Overall_M_BrachyMT4_0 <- Ref_Overall_M_BrachyMT4[1,2]
Ref_Overall_M_BrachyMT4_0[is.na(Ref_Overall_M_BrachyMT4_0)] <-0
Ref_Overall_M_BrachyMT4_1 <- Ref_Overall_M_BrachyMT4[2,2]
Ref_Overall_M_BrachyMT4_1[is.na(Ref_Overall_M_BrachyMT4_1)] <-0
Ref_Overall_M_BrachyMT4_Total <- 
  Ref_Overall_M_BrachyMT4_0 + Ref_Overall_M_BrachyMT4_1
Ref_Overall_M_BrachyMT4_Affected <- 
  round((Ref_Overall_M_BrachyMT4_1/Ref_Overall_M_BrachyMT4_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_BrachyMT4 <- Ref_Overall_F %>% count(BrachyMT4)
Ref_Overall_F_BrachyMT4 <- 
  Ref_Overall_F_BrachyMT4[(Ref_Overall_F_BrachyMT4$BrachyMT4 != 6 & Ref_Overall_F_BrachyMT4$BrachyMT4 != 9),]
Ref_Overall_F_BrachyMT4_0 <- Ref_Overall_F_BrachyMT4[1,2]
Ref_Overall_F_BrachyMT4_0[is.na(Ref_Overall_F_BrachyMT4_0)] <-0
Ref_Overall_F_BrachyMT4_1 <- Ref_Overall_F_BrachyMT4[2,2]
Ref_Overall_F_BrachyMT4_1[is.na(Ref_Overall_F_BrachyMT4_1)] <-0
Ref_Overall_F_BrachyMT4_Total <- 
  Ref_Overall_F_BrachyMT4_0 + Ref_Overall_F_BrachyMT4_1
Ref_Overall_F_BrachyMT4_Affected <- 
  round((Ref_Overall_F_BrachyMT4_1/Ref_Overall_F_BrachyMT4_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_BrachyMT4_row <- 
  c("BrachyMT4", Ref_Overall_BrachyMT4_Total, Ref_Overall_BrachyMT4_1, 
     Ref_Overall_BrachyMT4_Affected, Ref_Overall_M_BrachyMT4_Total, 
     Ref_Overall_M_BrachyMT4_1, Ref_Overall_M_BrachyMT4_Affected, 
     Ref_Overall_F_BrachyMT4_Total, Ref_Overall_F_BrachyMT4_1, 
     Ref_Overall_F_BrachyMT4_Affected)

## BrachyPP1 ##
# Overall N and frequency #
Ref_Overall_BrachyPP1 <- Ref_Overall %>% count(BrachyPP1)
Ref_Overall_BrachyPP1 <- 
  Ref_Overall_BrachyPP1[(Ref_Overall_BrachyPP1$BrachyPP1 != 6 & Ref_Overall_BrachyPP1$BrachyPP1 != 9),]
Ref_Overall_BrachyPP1_0 <- Ref_Overall_BrachyPP1[1,2]
Ref_Overall_BrachyPP1_0[is.na(Ref_Overall_BrachyPP1_0)] <-0
Ref_Overall_BrachyPP1_1 <- Ref_Overall_BrachyPP1[2,2]
Ref_Overall_BrachyPP1_1[is.na(Ref_Overall_BrachyPP1_1)] <-0
Ref_Overall_BrachyPP1_Total <- Ref_Overall_BrachyPP1_0 + Ref_Overall_BrachyPP1_1
Ref_Overall_BrachyPP1_Affected <- 
  round((Ref_Overall_BrachyPP1_1/Ref_Overall_BrachyPP1_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_BrachyPP1 <- Ref_Overall_M %>% count(BrachyPP1)
Ref_Overall_M_BrachyPP1 <- 
  Ref_Overall_M_BrachyPP1[(Ref_Overall_M_BrachyPP1$BrachyPP1 != 6 & Ref_Overall_M_BrachyPP1$BrachyPP1 != 9),]
Ref_Overall_M_BrachyPP1_0 <- Ref_Overall_M_BrachyPP1[1,2]
Ref_Overall_M_BrachyPP1_0[is.na(Ref_Overall_M_BrachyPP1_0)] <-0
Ref_Overall_M_BrachyPP1_1 <- Ref_Overall_M_BrachyPP1[2,2]
Ref_Overall_M_BrachyPP1_1[is.na(Ref_Overall_M_BrachyPP1_1)] <-0
Ref_Overall_M_BrachyPP1_Total <- 
  Ref_Overall_M_BrachyPP1_0 + Ref_Overall_M_BrachyPP1_1
Ref_Overall_M_BrachyPP1_Affected <- 
  round((Ref_Overall_M_BrachyPP1_1/Ref_Overall_M_BrachyPP1_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_BrachyPP1 <- Ref_Overall_F %>% count(BrachyPP1)
Ref_Overall_F_BrachyPP1 <- 
  Ref_Overall_F_BrachyPP1[(Ref_Overall_F_BrachyPP1$BrachyPP1 != 6 & Ref_Overall_F_BrachyPP1$BrachyPP1 != 9),]
Ref_Overall_F_BrachyPP1_0 <- Ref_Overall_F_BrachyPP1[1,2]
Ref_Overall_F_BrachyPP1_0[is.na(Ref_Overall_F_BrachyPP1_0)] <-0
Ref_Overall_F_BrachyPP1_1 <- Ref_Overall_F_BrachyPP1[2,2]
Ref_Overall_F_BrachyPP1_1[is.na(Ref_Overall_F_BrachyPP1_1)] <-0
Ref_Overall_F_BrachyPP1_Total <- 
  Ref_Overall_F_BrachyPP1_0 + Ref_Overall_F_BrachyPP1_1
Ref_Overall_F_BrachyPP1_Affected <- 
  round((Ref_Overall_F_BrachyPP1_1/Ref_Overall_F_BrachyPP1_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_BrachyPP1_row <- 
  c("BrachyPP1", Ref_Overall_BrachyPP1_Total, Ref_Overall_BrachyPP1_1, 
     Ref_Overall_BrachyPP1_Affected, Ref_Overall_M_BrachyPP1_Total, 
     Ref_Overall_M_BrachyPP1_1, Ref_Overall_M_BrachyPP1_Affected, 
     Ref_Overall_F_BrachyPP1_Total, Ref_Overall_F_BrachyPP1_1, 
     Ref_Overall_F_BrachyPP1_Affected)

## CalcCubCoal ##
# Overall N and frequency #
Ref_Overall_CalcCubCoal <- Ref_Overall %>% count(CalcCubCoal)
Ref_Overall_CalcCubCoal <- 
  Ref_Overall_CalcCubCoal[(Ref_Overall_CalcCubCoal$CalcCubCoal != 6 & Ref_Overall_CalcCubCoal$CalcCubCoal != 9),]
Ref_Overall_CalcCubCoal_0 <- Ref_Overall_CalcCubCoal[1,2]
Ref_Overall_CalcCubCoal_0[is.na(Ref_Overall_CalcCubCoal_0)] <-0
Ref_Overall_CalcCubCoal_1 <- Ref_Overall_CalcCubCoal[2,2]
Ref_Overall_CalcCubCoal_1[is.na(Ref_Overall_CalcCubCoal_1)] <-0
Ref_Overall_CalcCubCoal_Total <- 
  Ref_Overall_CalcCubCoal_0 + Ref_Overall_CalcCubCoal_1
Ref_Overall_CalcCubCoal_Affected <- 
  round((Ref_Overall_CalcCubCoal_1/Ref_Overall_CalcCubCoal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_CalcCubCoal <- Ref_Overall_M %>% count(CalcCubCoal)
Ref_Overall_M_CalcCubCoal <- 
  Ref_Overall_M_CalcCubCoal[(Ref_Overall_M_CalcCubCoal$CalcCubCoal != 6 & Ref_Overall_M_CalcCubCoal$CalcCubCoal != 9),]
Ref_Overall_M_CalcCubCoal_0 <- Ref_Overall_M_CalcCubCoal[1,2]
Ref_Overall_M_CalcCubCoal_0[is.na(Ref_Overall_M_CalcCubCoal_0)] <-0
Ref_Overall_M_CalcCubCoal_1 <- Ref_Overall_M_CalcCubCoal[2,2]
Ref_Overall_M_CalcCubCoal_1[is.na(Ref_Overall_M_CalcCubCoal_1)] <-0
Ref_Overall_M_CalcCubCoal_Total <- 
  Ref_Overall_M_CalcCubCoal_0 + Ref_Overall_M_CalcCubCoal_1
Ref_Overall_M_CalcCubCoal_Affected <- 
  round((Ref_Overall_M_CalcCubCoal_1/Ref_Overall_M_CalcCubCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CalcCubCoal <- Ref_Overall_F %>% count(CalcCubCoal)
Ref_Overall_F_CalcCubCoal <- 
  Ref_Overall_F_CalcCubCoal[(Ref_Overall_F_CalcCubCoal$CalcCubCoal != 6 & Ref_Overall_F_CalcCubCoal$CalcCubCoal != 9),]
Ref_Overall_F_CalcCubCoal_0 <- Ref_Overall_F_CalcCubCoal[1,2]
Ref_Overall_F_CalcCubCoal_0[is.na(Ref_Overall_F_CalcCubCoal_0)] <-0
Ref_Overall_F_CalcCubCoal_1 <- Ref_Overall_F_CalcCubCoal[2,2]
Ref_Overall_F_CalcCubCoal_1[is.na(Ref_Overall_F_CalcCubCoal_1)] <-0
Ref_Overall_F_CalcCubCoal_Total <- 
  Ref_Overall_F_CalcCubCoal_0 + Ref_Overall_F_CalcCubCoal_1
Ref_Overall_F_CalcCubCoal_Affected <- 
  round((Ref_Overall_F_CalcCubCoal_1/Ref_Overall_F_CalcCubCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_CalcCubCoal_row <- 
  c("CalcCubCoal", Ref_Overall_CalcCubCoal_Total, Ref_Overall_CalcCubCoal_1, 
     Ref_Overall_CalcCubCoal_Affected, Ref_Overall_M_CalcCubCoal_Total, 
     Ref_Overall_M_CalcCubCoal_1, Ref_Overall_M_CalcCubCoal_Affected, 
     Ref_Overall_F_CalcCubCoal_Total, Ref_Overall_F_CalcCubCoal_1, 
     Ref_Overall_F_CalcCubCoal_Affected)

## TaloNavCoal ##
# Overall N and frequency #
Ref_Overall_TaloNavCoal <- Ref_Overall %>% count(TaloNavCoal)
Ref_Overall_TaloNavCoal <- 
  Ref_Overall_TaloNavCoal[(Ref_Overall_TaloNavCoal$TaloNavCoal != 6 & Ref_Overall_TaloNavCoal$TaloNavCoal != 9),]
Ref_Overall_TaloNavCoal_0 <- Ref_Overall_TaloNavCoal[1,2]
Ref_Overall_TaloNavCoal_0[is.na(Ref_Overall_TaloNavCoal_0)] <-0
Ref_Overall_TaloNavCoal_1 <- Ref_Overall_TaloNavCoal[2,2]
Ref_Overall_TaloNavCoal_1[is.na(Ref_Overall_TaloNavCoal_1)] <-0
Ref_Overall_TaloNavCoal_Total <- 
  Ref_Overall_TaloNavCoal_0 + Ref_Overall_TaloNavCoal_1
Ref_Overall_TaloNavCoal_Affected <- 
  round((Ref_Overall_TaloNavCoal_1/Ref_Overall_TaloNavCoal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_TaloNavCoal <- Ref_Overall_M %>% count(TaloNavCoal)
Ref_Overall_M_TaloNavCoal <- 
  Ref_Overall_M_TaloNavCoal[(Ref_Overall_M_TaloNavCoal$TaloNavCoal != 6 & Ref_Overall_M_TaloNavCoal$TaloNavCoal != 9),]
Ref_Overall_M_TaloNavCoal_0 <- Ref_Overall_M_TaloNavCoal[1,2]
Ref_Overall_M_TaloNavCoal_0[is.na(Ref_Overall_M_TaloNavCoal_0)] <-0
Ref_Overall_M_TaloNavCoal_1 <- Ref_Overall_M_TaloNavCoal[2,2]
Ref_Overall_M_TaloNavCoal_1[is.na(Ref_Overall_M_TaloNavCoal_1)] <-0
Ref_Overall_M_TaloNavCoal_Total <- 
  Ref_Overall_M_TaloNavCoal_0 + Ref_Overall_M_TaloNavCoal_1
Ref_Overall_M_TaloNavCoal_Affected <- 
  round((Ref_Overall_M_TaloNavCoal_1/Ref_Overall_M_TaloNavCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_TaloNavCoal <- Ref_Overall_F %>% count(TaloNavCoal)
Ref_Overall_F_TaloNavCoal <- 
  Ref_Overall_F_TaloNavCoal[(Ref_Overall_F_TaloNavCoal$TaloNavCoal != 6 & Ref_Overall_F_TaloNavCoal$TaloNavCoal != 9),]
Ref_Overall_F_TaloNavCoal_0 <- Ref_Overall_F_TaloNavCoal[1,2]
Ref_Overall_F_TaloNavCoal_0[is.na(Ref_Overall_F_TaloNavCoal_0)] <-0
Ref_Overall_F_TaloNavCoal_1 <- Ref_Overall_F_TaloNavCoal[2,2]
Ref_Overall_F_TaloNavCoal_1[is.na(Ref_Overall_F_TaloNavCoal_1)] <-0
Ref_Overall_F_TaloNavCoal_Total <- 
  Ref_Overall_F_TaloNavCoal_0 + Ref_Overall_F_TaloNavCoal_1
Ref_Overall_F_TaloNavCoal_Affected <- 
  round((Ref_Overall_F_TaloNavCoal_1/Ref_Overall_F_TaloNavCoal_Total)*100, 2)

# overview of Taloulated frequencies for table (see below) #
Ref_Overall_TaloNavCoal_row <- 
  c("TaloNavCoal", Ref_Overall_TaloNavCoal_Total, Ref_Overall_TaloNavCoal_1, 
    Ref_Overall_TaloNavCoal_Affected, Ref_Overall_M_TaloNavCoal_Total, 
    Ref_Overall_M_TaloNavCoal_1, Ref_Overall_M_TaloNavCoal_Affected, 
    Ref_Overall_F_TaloNavCoal_Total, Ref_Overall_F_TaloNavCoal_1, 
    Ref_Overall_F_TaloNavCoal_Affected)

## CF2CF3Coal ##
# Overall N and frequency #
Ref_Overall_CF2CF3Coal <- Ref_Overall %>% count(CF2CF3Coal)
Ref_Overall_CF2CF3Coal <- 
  Ref_Overall_CF2CF3Coal[(Ref_Overall_CF2CF3Coal$CF2CF3Coal != 6 & Ref_Overall_CF2CF3Coal$CF2CF3Coal != 9),]
Ref_Overall_CF2CF3Coal_0 <- Ref_Overall_CF2CF3Coal[1,2]
Ref_Overall_CF2CF3Coal_0[is.na(Ref_Overall_CF2CF3Coal_0)] <-0
Ref_Overall_CF2CF3Coal_1 <- Ref_Overall_CF2CF3Coal[2,2]
Ref_Overall_CF2CF3Coal_1[is.na(Ref_Overall_CF2CF3Coal_1)] <-0
Ref_Overall_CF2CF3Coal_Total <- 
  Ref_Overall_CF2CF3Coal_0 + Ref_Overall_CF2CF3Coal_1
Ref_Overall_CF2CF3Coal_Affected <- 
  round((Ref_Overall_CF2CF3Coal_1/Ref_Overall_CF2CF3Coal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_CF2CF3Coal <- Ref_Overall_M %>% count(CF2CF3Coal)
Ref_Overall_M_CF2CF3Coal <- 
  Ref_Overall_M_CF2CF3Coal[(Ref_Overall_M_CF2CF3Coal$CF2CF3Coal != 6 & Ref_Overall_M_CF2CF3Coal$CF2CF3Coal != 9),]
Ref_Overall_M_CF2CF3Coal_0 <- Ref_Overall_M_CF2CF3Coal[1,2]
Ref_Overall_M_CF2CF3Coal_0[is.na(Ref_Overall_M_CF2CF3Coal_0)] <-0
Ref_Overall_M_CF2CF3Coal_1 <- Ref_Overall_M_CF2CF3Coal[2,2]
Ref_Overall_M_CF2CF3Coal_1[is.na(Ref_Overall_M_CF2CF3Coal_1)] <-0
Ref_Overall_M_CF2CF3Coal_Total <- 
  Ref_Overall_M_CF2CF3Coal_0 + Ref_Overall_M_CF2CF3Coal_1
Ref_Overall_M_CF2CF3Coal_Affected <- 
  round((Ref_Overall_M_CF2CF3Coal_1/Ref_Overall_M_CF2CF3Coal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CF2CF3Coal <- Ref_Overall_F %>% count(CF2CF3Coal)
Ref_Overall_F_CF2CF3Coal <- 
  Ref_Overall_F_CF2CF3Coal[(Ref_Overall_F_CF2CF3Coal$CF2CF3Coal != 6 & Ref_Overall_F_CF2CF3Coal$CF2CF3Coal != 9),]
Ref_Overall_F_CF2CF3Coal_0 <- Ref_Overall_F_CF2CF3Coal[1,2]
Ref_Overall_F_CF2CF3Coal_0[is.na(Ref_Overall_F_CF2CF3Coal_0)] <-0
Ref_Overall_F_CF2CF3Coal_1 <- Ref_Overall_F_CF2CF3Coal[2,2]
Ref_Overall_F_CF2CF3Coal_1[is.na(Ref_Overall_F_CF2CF3Coal_1)] <-0
Ref_Overall_F_CF2CF3Coal_Total <- 
  Ref_Overall_F_CF2CF3Coal_0 + Ref_Overall_F_CF2CF3Coal_1
Ref_Overall_F_CF2CF3Coal_Affected <- 
  round((Ref_Overall_F_CF2CF3Coal_1/Ref_Overall_F_CF2CF3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_CF2CF3Coal_row <- 
  c("CF2CF3Coal", Ref_Overall_CF2CF3Coal_Total, Ref_Overall_CF2CF3Coal_1, 
    Ref_Overall_CF2CF3Coal_Affected, Ref_Overall_M_CF2CF3Coal_Total, 
    Ref_Overall_M_CF2CF3Coal_1, Ref_Overall_M_CF2CF3Coal_Affected, 
    Ref_Overall_F_CF2CF3Coal_Total, Ref_Overall_F_CF2CF3Coal_1, 
    Ref_Overall_F_CF2CF3Coal_Affected)

## CF3MT3Coal ##
# Overall N and frequency #
Ref_Overall_CF3MT3Coal <- Ref_Overall %>% count(CF3MT3Coal)
Ref_Overall_CF3MT3Coal <- 
  Ref_Overall_CF3MT3Coal[(Ref_Overall_CF3MT3Coal$CF3MT3Coal != 6 & Ref_Overall_CF3MT3Coal$CF3MT3Coal != 9),]
Ref_Overall_CF3MT3Coal_0 <- Ref_Overall_CF3MT3Coal[1,2]
Ref_Overall_CF3MT3Coal_0[is.na(Ref_Overall_CF3MT3Coal_0)] <-0
Ref_Overall_CF3MT3Coal_1 <- Ref_Overall_CF3MT3Coal[2,2]
Ref_Overall_CF3MT3Coal_1[is.na(Ref_Overall_CF3MT3Coal_1)] <-0
Ref_Overall_CF3MT3Coal_Total <- 
  Ref_Overall_CF3MT3Coal_0 + Ref_Overall_CF3MT3Coal_1
Ref_Overall_CF3MT3Coal_Affected <- 
  round((Ref_Overall_CF3MT3Coal_1/Ref_Overall_CF3MT3Coal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_CF3MT3Coal <- Ref_Overall_M %>% count(CF3MT3Coal)
Ref_Overall_M_CF3MT3Coal <- 
  Ref_Overall_M_CF3MT3Coal[(Ref_Overall_M_CF3MT3Coal$CF3MT3Coal != 6 & Ref_Overall_M_CF3MT3Coal$CF3MT3Coal != 9),]
Ref_Overall_M_CF3MT3Coal_0 <- Ref_Overall_M_CF3MT3Coal[1,2]
Ref_Overall_M_CF3MT3Coal_0[is.na(Ref_Overall_M_CF3MT3Coal_0)] <-0
Ref_Overall_M_CF3MT3Coal_1 <- Ref_Overall_M_CF3MT3Coal[2,2]
Ref_Overall_M_CF3MT3Coal_1[is.na(Ref_Overall_M_CF3MT3Coal_1)] <-0
Ref_Overall_M_CF3MT3Coal_Total <- 
  Ref_Overall_M_CF3MT3Coal_0 + Ref_Overall_M_CF3MT3Coal_1
Ref_Overall_M_CF3MT3Coal_Affected <- 
  round((Ref_Overall_M_CF3MT3Coal_1/Ref_Overall_M_CF3MT3Coal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CF3MT3Coal <- Ref_Overall_F %>% count(CF3MT3Coal)
Ref_Overall_F_CF3MT3Coal <- 
  Ref_Overall_F_CF3MT3Coal[(Ref_Overall_F_CF3MT3Coal$CF3MT3Coal != 6 & Ref_Overall_F_CF3MT3Coal$CF3MT3Coal != 9),]
Ref_Overall_F_CF3MT3Coal_0 <- Ref_Overall_F_CF3MT3Coal[1,2]
Ref_Overall_F_CF3MT3Coal_0[is.na(Ref_Overall_F_CF3MT3Coal_0)] <-0
Ref_Overall_F_CF3MT3Coal_1 <- Ref_Overall_F_CF3MT3Coal[2,2]
Ref_Overall_F_CF3MT3Coal_1[is.na(Ref_Overall_F_CF3MT3Coal_1)] <-0
Ref_Overall_F_CF3MT3Coal_Total <- 
  Ref_Overall_F_CF3MT3Coal_0 + Ref_Overall_F_CF3MT3Coal_1
Ref_Overall_F_CF3MT3Coal_Affected <- 
  round((Ref_Overall_F_CF3MT3Coal_1/Ref_Overall_F_CF3MT3Coal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_CF3MT3Coal_row <- 
  c("CF3MT3Coal", Ref_Overall_CF3MT3Coal_Total, Ref_Overall_CF3MT3Coal_1, 
    Ref_Overall_CF3MT3Coal_Affected, Ref_Overall_M_CF3MT3Coal_Total, 
    Ref_Overall_M_CF3MT3Coal_1, Ref_Overall_M_CF3MT3Coal_Affected, 
    Ref_Overall_F_CF3MT3Coal_Total, Ref_Overall_F_CF3MT3Coal_1, 
    Ref_Overall_F_CF3MT3Coal_Affected)

## OsIntermet ##
# Overall N and frequency #
Ref_Overall_OsIntermet <- Ref_Overall %>% count(OsIntermet)
Ref_Overall_OsIntermet <- 
  Ref_Overall_OsIntermet[(Ref_Overall_OsIntermet$OsIntermet != 6 & Ref_Overall_OsIntermet$OsIntermet != 9),]
Ref_Overall_OsIntermet_0 <- Ref_Overall_OsIntermet[1,2]
Ref_Overall_OsIntermet_0[is.na(Ref_Overall_OsIntermet_0)] <-0
Ref_Overall_OsIntermet_1 <- Ref_Overall_OsIntermet[2,2]
Ref_Overall_OsIntermet_1[is.na(Ref_Overall_OsIntermet_1)] <-0
Ref_Overall_OsIntermet_Total <- 
  Ref_Overall_OsIntermet_0 + Ref_Overall_OsIntermet_1
Ref_Overall_OsIntermet_Affected <- 
  round((Ref_Overall_OsIntermet_1/Ref_Overall_OsIntermet_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_OsIntermet <- Ref_Overall_M %>% count(OsIntermet)
Ref_Overall_M_OsIntermet <- 
  Ref_Overall_M_OsIntermet[(Ref_Overall_M_OsIntermet$OsIntermet != 6 & Ref_Overall_M_OsIntermet$OsIntermet != 9),]
Ref_Overall_M_OsIntermet_0 <- Ref_Overall_M_OsIntermet[1,2]
Ref_Overall_M_OsIntermet_0[is.na(Ref_Overall_M_OsIntermet_0)] <-0
Ref_Overall_M_OsIntermet_1 <- Ref_Overall_M_OsIntermet[2,2]
Ref_Overall_M_OsIntermet_1[is.na(Ref_Overall_M_OsIntermet_1)] <-0
Ref_Overall_M_OsIntermet_Total <- 
  Ref_Overall_M_OsIntermet_0 + Ref_Overall_M_OsIntermet_1
Ref_Overall_M_OsIntermet_Affected <- 
  round((Ref_Overall_M_OsIntermet_1/Ref_Overall_M_OsIntermet_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_OsIntermet <- Ref_Overall_F %>% count(OsIntermet)
Ref_Overall_F_OsIntermet <- 
  Ref_Overall_F_OsIntermet[(Ref_Overall_F_OsIntermet$OsIntermet != 6 & Ref_Overall_F_OsIntermet$OsIntermet != 9),]
Ref_Overall_F_OsIntermet_0 <- Ref_Overall_F_OsIntermet[1,2]
Ref_Overall_F_OsIntermet_0[is.na(Ref_Overall_F_OsIntermet_0)] <-0
Ref_Overall_F_OsIntermet_1 <- Ref_Overall_F_OsIntermet[2,2]
Ref_Overall_F_OsIntermet_1[is.na(Ref_Overall_F_OsIntermet_1)] <-0
Ref_Overall_F_OsIntermet_Total <- 
  Ref_Overall_F_OsIntermet_0 + Ref_Overall_F_OsIntermet_1
Ref_Overall_F_OsIntermet_Affected <- 
  round((Ref_Overall_F_OsIntermet_1/Ref_Overall_F_OsIntermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_OsIntermet_row <- 
  c("OsIntermet", Ref_Overall_OsIntermet_Total, Ref_Overall_OsIntermet_1, 
    Ref_Overall_OsIntermet_Affected, Ref_Overall_M_OsIntermet_Total, 
    Ref_Overall_M_OsIntermet_1, Ref_Overall_M_OsIntermet_Affected, 
    Ref_Overall_F_OsIntermet_Total, Ref_Overall_F_OsIntermet_1, 
    Ref_Overall_F_OsIntermet_Affected)

# OsIntermet unclustered: CF1Intermet #
# Overall N and frequency #
Ref_Overall_CF1Intermet <- Ref_Overall %>% count(CF1Intermet)
Ref_Overall_CF1Intermet <- 
  Ref_Overall_CF1Intermet[(Ref_Overall_CF1Intermet$CF1Intermet != 6 & Ref_Overall_CF1Intermet$CF1Intermet != 9),]
Ref_Overall_CF1Intermet_0 <- Ref_Overall_CF1Intermet[1,2]
Ref_Overall_CF1Intermet_0[is.na(Ref_Overall_CF1Intermet_0)] <-0
Ref_Overall_CF1Intermet_1 <- Ref_Overall_CF1Intermet[2,2]
Ref_Overall_CF1Intermet_1[is.na(Ref_Overall_CF1Intermet_1)] <-0
Ref_Overall_CF1Intermet_Total <- 
  na.omit(Ref_Overall_CF1Intermet_0) + na.omit(Ref_Overall_CF1Intermet_1)
Ref_Overall_CF1Intermet_Affected <- 
  round((Ref_Overall_CF1Intermet_1/Ref_Overall_CF1Intermet_Total)*100, 2)


# Male N and frequency
Ref_Overall_M_CF1Intermet <- Ref_Overall_M %>% count(CF1Intermet)
Ref_Overall_M_CF1Intermet <- 
  Ref_Overall_M_CF1Intermet[(Ref_Overall_M_CF1Intermet$CF1Intermet != 6 & Ref_Overall_M_CF1Intermet$CF1Intermet != 9),]
Ref_Overall_M_CF1Intermet_0 <- Ref_Overall_M_CF1Intermet[1,2]
Ref_Overall_M_CF1Intermet_0[is.na(Ref_Overall_M_CF1Intermet_0)] <-0
Ref_Overall_M_CF1Intermet_1 <- Ref_Overall_M_CF1Intermet[2,2]
Ref_Overall_M_CF1Intermet_1[is.na(Ref_Overall_M_CF1Intermet_1)] <-0
Ref_Overall_M_CF1Intermet_Total <- 
  Ref_Overall_M_CF1Intermet_0 + Ref_Overall_M_CF1Intermet_1
Ref_Overall_M_CF1Intermet_Affected <- 
  round((Ref_Overall_M_CF1Intermet_1/Ref_Overall_M_CF1Intermet_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CF1Intermet <- Ref_Overall_F %>% count(CF1Intermet)
Ref_Overall_F_CF1Intermet <- 
  Ref_Overall_F_CF1Intermet[(Ref_Overall_F_CF1Intermet$CF1Intermet != 6 & Ref_Overall_F_CF1Intermet$CF1Intermet != 9),]
Ref_Overall_F_CF1Intermet_0 <- Ref_Overall_F_CF1Intermet[1,2]
Ref_Overall_F_CF1Intermet_0[is.na(Ref_Overall_F_CF1Intermet_0)] <-0
Ref_Overall_F_CF1Intermet_1 <- Ref_Overall_F_CF1Intermet[2,2]
Ref_Overall_F_CF1Intermet_1[is.na(Ref_Overall_F_CF1Intermet_1)] <-0
Ref_Overall_F_CF1Intermet_Total <- 
  Ref_Overall_F_CF1Intermet_0 + Ref_Overall_F_CF1Intermet_1
Ref_Overall_F_CF1Intermet_Affected <- 
  round((Ref_Overall_F_CF1Intermet_1/Ref_Overall_F_CF1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_CF1Intermet_row <- c("CF1Intermet", Ref_Overall_CF1Intermet_Total, 
                                 Ref_Overall_CF1Intermet_1, 
                                 Ref_Overall_CF1Intermet_Affected, 
                                 Ref_Overall_M_CF1Intermet_Total, 
                                 Ref_Overall_M_CF1Intermet_1, 
                                 Ref_Overall_M_CF1Intermet_Affected, 
                                 Ref_Overall_F_CF1Intermet_Total, 
                                 Ref_Overall_F_CF1Intermet_1, 
                                 Ref_Overall_F_CF1Intermet_Affected)

# OsIntermet unclustered: MT1Intermet #
# Overall N and frequency #
Ref_Overall_MT1Intermet <- Ref_Overall %>% count(MT1Intermet)
Ref_Overall_MT1Intermet <- 
  Ref_Overall_MT1Intermet[(Ref_Overall_MT1Intermet$MT1Intermet != 6 & Ref_Overall_MT1Intermet$MT1Intermet != 9),]
Ref_Overall_MT1Intermet_0 <- Ref_Overall_MT1Intermet[1,2]
Ref_Overall_MT1Intermet_0[is.na(Ref_Overall_MT1Intermet_0)] <-0
Ref_Overall_MT1Intermet_1 <- Ref_Overall_MT1Intermet[2,2]
Ref_Overall_MT1Intermet_1[is.na(Ref_Overall_MT1Intermet_1)] <-0
Ref_Overall_MT1Intermet_Total <- 
  na.omit(Ref_Overall_MT1Intermet_0) + na.omit(Ref_Overall_MT1Intermet_1)
Ref_Overall_MT1Intermet_Affected <- 
  round((Ref_Overall_MT1Intermet_1/Ref_Overall_MT1Intermet_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_MT1Intermet <- Ref_Overall_M %>% count(MT1Intermet)
Ref_Overall_M_MT1Intermet <- 
  Ref_Overall_M_MT1Intermet[(Ref_Overall_M_MT1Intermet$MT1Intermet != 6 & Ref_Overall_M_MT1Intermet$MT1Intermet != 9),]
Ref_Overall_M_MT1Intermet_0 <- Ref_Overall_M_MT1Intermet[1,2]
Ref_Overall_M_MT1Intermet_0[is.na(Ref_Overall_M_MT1Intermet_0)] <-0
Ref_Overall_M_MT1Intermet_1 <- Ref_Overall_M_MT1Intermet[2,2]
Ref_Overall_M_MT1Intermet_1[is.na(Ref_Overall_M_MT1Intermet_1)] <-0
Ref_Overall_M_MT1Intermet_Total <- 
  Ref_Overall_M_MT1Intermet_0 + Ref_Overall_M_MT1Intermet_1
Ref_Overall_M_MT1Intermet_Affected <- 
  round((Ref_Overall_M_MT1Intermet_1/Ref_Overall_M_MT1Intermet_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_MT1Intermet <- Ref_Overall_F %>% count(MT1Intermet)
Ref_Overall_F_MT1Intermet <- 
  Ref_Overall_F_MT1Intermet[(Ref_Overall_F_MT1Intermet$MT1Intermet != 6 & Ref_Overall_F_MT1Intermet$MT1Intermet != 9),]
Ref_Overall_F_MT1Intermet_0 <- Ref_Overall_F_MT1Intermet[1,2]
Ref_Overall_F_MT1Intermet_0[is.na(Ref_Overall_F_MT1Intermet_0)] <-0
Ref_Overall_F_MT1Intermet_1 <- Ref_Overall_F_MT1Intermet[2,2]
Ref_Overall_F_MT1Intermet_1[is.na(Ref_Overall_F_MT1Intermet_1)] <-0
Ref_Overall_F_MT1Intermet_Total <- 
  Ref_Overall_F_MT1Intermet_0 + Ref_Overall_F_MT1Intermet_1
Ref_Overall_F_MT1Intermet_Affected <- 
  round((Ref_Overall_F_MT1Intermet_1/Ref_Overall_F_MT1Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_MT1Intermet_row <- c("MT1Intermet", Ref_Overall_MT1Intermet_Total, 
                                 Ref_Overall_MT1Intermet_1, 
                                 Ref_Overall_MT1Intermet_Affected, 
                                 Ref_Overall_M_MT1Intermet_Total, 
                                 Ref_Overall_M_MT1Intermet_1, 
                                 Ref_Overall_M_MT1Intermet_Affected, 
                                 Ref_Overall_F_MT1Intermet_Total, 
                                 Ref_Overall_F_MT1Intermet_1, 
                                 Ref_Overall_F_MT1Intermet_Affected)

# OsIntermet unclustered: MT2Intermet #
# Overall N and frequency #
Ref_Overall_MT2Intermet <- Ref_Overall %>% count(MT2Intermet)
Ref_Overall_MT2Intermet <- 
  Ref_Overall_MT2Intermet[(Ref_Overall_MT2Intermet$MT2Intermet != 6 & Ref_Overall_MT2Intermet$MT2Intermet != 9),]
Ref_Overall_MT2Intermet_0 <- Ref_Overall_MT2Intermet[1,2]
Ref_Overall_MT2Intermet_0[is.na(Ref_Overall_MT2Intermet_0)] <-0
Ref_Overall_MT2Intermet_1 <- Ref_Overall_MT2Intermet[2,2]
Ref_Overall_MT2Intermet_1[is.na(Ref_Overall_MT2Intermet_1)] <-0
Ref_Overall_MT2Intermet_Total <- 
  na.omit(Ref_Overall_MT2Intermet_0) + na.omit(Ref_Overall_MT2Intermet_1)
Ref_Overall_MT2Intermet_Affected <- 
  round((Ref_Overall_MT2Intermet_1/Ref_Overall_MT2Intermet_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_MT2Intermet <- Ref_Overall_M %>% count(MT2Intermet)
Ref_Overall_M_MT2Intermet <- 
  Ref_Overall_M_MT2Intermet[(Ref_Overall_M_MT2Intermet$MT2Intermet != 6 & Ref_Overall_M_MT2Intermet$MT2Intermet != 9),]
Ref_Overall_M_MT2Intermet_0 <- Ref_Overall_M_MT2Intermet[1,2]
Ref_Overall_M_MT2Intermet_0[is.na(Ref_Overall_M_MT2Intermet_0)] <-0
Ref_Overall_M_MT2Intermet_1 <- Ref_Overall_M_MT2Intermet[2,2]
Ref_Overall_M_MT2Intermet_1[is.na(Ref_Overall_M_MT2Intermet_1)] <-0
Ref_Overall_M_MT2Intermet_Total <- 
  Ref_Overall_M_MT2Intermet_0 + Ref_Overall_M_MT2Intermet_1
Ref_Overall_M_MT2Intermet_Affected <- 
  round((Ref_Overall_M_MT2Intermet_1/Ref_Overall_M_MT2Intermet_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_MT2Intermet <- Ref_Overall_F %>% count(MT2Intermet)
Ref_Overall_F_MT2Intermet <- 
  Ref_Overall_F_MT2Intermet[(Ref_Overall_F_MT2Intermet$MT2Intermet != 6 & Ref_Overall_F_MT2Intermet$MT2Intermet != 9),]
Ref_Overall_F_MT2Intermet_0 <- Ref_Overall_F_MT2Intermet[1,2]
Ref_Overall_F_MT2Intermet_0[is.na(Ref_Overall_F_MT2Intermet_0)] <-0
Ref_Overall_F_MT2Intermet_1 <- Ref_Overall_F_MT2Intermet[2,2]
Ref_Overall_F_MT2Intermet_1[is.na(Ref_Overall_F_MT2Intermet_1)] <-0
Ref_Overall_F_MT2Intermet_Total <- 
  Ref_Overall_F_MT2Intermet_0 + Ref_Overall_F_MT2Intermet_1
Ref_Overall_F_MT2Intermet_Affected <- 
  round((Ref_Overall_F_MT2Intermet_1/Ref_Overall_F_MT2Intermet_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_MT2Intermet_row <- c("MT2Intermet", Ref_Overall_MT2Intermet_Total, 
                                 Ref_Overall_MT2Intermet_1, 
                                 Ref_Overall_MT2Intermet_Affected, 
                                 Ref_Overall_M_MT2Intermet_Total, 
                                 Ref_Overall_M_MT2Intermet_1, 
                                 Ref_Overall_M_MT2Intermet_Affected, 
                                 Ref_Overall_F_MT2Intermet_Total, 
                                 Ref_Overall_F_MT2Intermet_1, 
                                 Ref_Overall_F_MT2Intermet_Affected)


## AnkleCoal ##
# Overall N and frequency #
Ref_Overall_AnkleCoal <- Ref_Overall %>% count(AnkleCoal)
Ref_Overall_AnkleCoal <- 
  Ref_Overall_AnkleCoal[(Ref_Overall_AnkleCoal$AnkleCoal != 6 & Ref_Overall_AnkleCoal$AnkleCoal != 9),]
Ref_Overall_AnkleCoal_0 <- Ref_Overall_AnkleCoal[1,2]; 
Ref_Overall_AnkleCoal_0[is.na(Ref_Overall_AnkleCoal_0)] <-0
Ref_Overall_AnkleCoal_1 <- Ref_Overall_AnkleCoal[2,2]; 
Ref_Overall_AnkleCoal_1[is.na(Ref_Overall_AnkleCoal_1)] <-0
Ref_Overall_AnkleCoal_Total <- 
  na.omit(Ref_Overall_AnkleCoal_0) + na.omit(Ref_Overall_AnkleCoal_1)
Ref_Overall_AnkleCoal_Affected <- 
  round((Ref_Overall_AnkleCoal_1/Ref_Overall_AnkleCoal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_AnkleCoal <- Ref_Overall_M %>% count(AnkleCoal)
Ref_Overall_M_AnkleCoal <- 
  Ref_Overall_M_AnkleCoal[(Ref_Overall_M_AnkleCoal$AnkleCoal != 6 & Ref_Overall_M_AnkleCoal$AnkleCoal != 9),]
Ref_Overall_M_AnkleCoal_0 <- Ref_Overall_M_AnkleCoal[1,2]
Ref_Overall_M_AnkleCoal_0[is.na(Ref_Overall_M_AnkleCoal_0)] <-0
Ref_Overall_M_AnkleCoal_1 <- Ref_Overall_M_AnkleCoal[2,2]
Ref_Overall_M_AnkleCoal_1[is.na(Ref_Overall_M_AnkleCoal_1)] <-0
Ref_Overall_M_AnkleCoal_Total <- 
  Ref_Overall_M_AnkleCoal_0 + Ref_Overall_M_AnkleCoal_1
Ref_Overall_M_AnkleCoal_Affected <- 
  round((Ref_Overall_M_AnkleCoal_1/Ref_Overall_M_AnkleCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_AnkleCoal <- Ref_Overall_F %>% count(AnkleCoal)
Ref_Overall_F_AnkleCoal <- 
  Ref_Overall_F_AnkleCoal[(Ref_Overall_F_AnkleCoal$AnkleCoal != 6 & Ref_Overall_F_AnkleCoal$AnkleCoal != 9),]
Ref_Overall_F_AnkleCoal_0 <- Ref_Overall_F_AnkleCoal[1,2]
Ref_Overall_F_AnkleCoal_0[is.na(Ref_Overall_F_AnkleCoal_0)] <-0
Ref_Overall_F_AnkleCoal_1 <- Ref_Overall_F_AnkleCoal[2,2]
Ref_Overall_F_AnkleCoal_1[is.na(Ref_Overall_F_AnkleCoal_1)] <-0
Ref_Overall_F_AnkleCoal_Total <- 
  Ref_Overall_F_AnkleCoal_0 + Ref_Overall_F_AnkleCoal_1
Ref_Overall_F_AnkleCoal_Affected <- 
  round((Ref_Overall_F_AnkleCoal_1/Ref_Overall_F_AnkleCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_AnkleCoal_row <- c("AnkleCoal", Ref_Overall_AnkleCoal_Total, 
                               Ref_Overall_AnkleCoal_1,
                               Ref_Overall_AnkleCoal_Affected, 
                               Ref_Overall_M_AnkleCoal_Total, 
                               Ref_Overall_M_AnkleCoal_1, 
                               Ref_Overall_M_AnkleCoal_Affected, 
                               Ref_Overall_F_AnkleCoal_Total, 
                               Ref_Overall_F_AnkleCoal_1, 
                               Ref_Overall_F_AnkleCoal_Affected)

# AnkleCoal unclustered: CalcNavCoal #
# Overall N and frequency #
Ref_Overall_CalcNavCoal <- Ref_Overall %>% count(CalcNavCoal)
Ref_Overall_CalcNavCoal <- 
  Ref_Overall_CalcNavCoal[(Ref_Overall_CalcNavCoal$CalcNavCoal != 6 & Ref_Overall_CalcNavCoal$CalcNavCoal != 9),]
Ref_Overall_CalcNavCoal_0 <- Ref_Overall_CalcNavCoal[1,2]
Ref_Overall_CalcNavCoal_0[is.na(Ref_Overall_CalcNavCoal_0)] <-0
Ref_Overall_CalcNavCoal_1 <- Ref_Overall_CalcNavCoal[2,2]
Ref_Overall_CalcNavCoal_1[is.na(Ref_Overall_CalcNavCoal_1)] <-0
Ref_Overall_CalcNavCoal_Total <- 
  na.omit(Ref_Overall_CalcNavCoal_0) + na.omit(Ref_Overall_CalcNavCoal_1)
Ref_Overall_CalcNavCoal_Affected <- 
  round((Ref_Overall_CalcNavCoal_1/Ref_Overall_CalcNavCoal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_CalcNavCoal <- Ref_Overall_M %>% count(CalcNavCoal)
Ref_Overall_M_CalcNavCoal <- 
  Ref_Overall_M_CalcNavCoal[(Ref_Overall_M_CalcNavCoal$CalcNavCoal != 6 & Ref_Overall_M_CalcNavCoal$CalcNavCoal != 9),]
Ref_Overall_M_CalcNavCoal_0 <- Ref_Overall_M_CalcNavCoal[1,2]
Ref_Overall_M_CalcNavCoal_0[is.na(Ref_Overall_M_CalcNavCoal_0)] <-0
Ref_Overall_M_CalcNavCoal_1 <- Ref_Overall_M_CalcNavCoal[2,2]
Ref_Overall_M_CalcNavCoal_1[is.na(Ref_Overall_M_CalcNavCoal_1)] <-0
Ref_Overall_M_CalcNavCoal_Total <- 
  Ref_Overall_M_CalcNavCoal_0 + Ref_Overall_M_CalcNavCoal_1
Ref_Overall_M_CalcNavCoal_Affected <- 
  round((Ref_Overall_M_CalcNavCoal_1/Ref_Overall_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CalcNavCoal <- Ref_Overall_F %>% count(CalcNavCoal)
Ref_Overall_F_CalcNavCoal <- 
  Ref_Overall_F_CalcNavCoal[(Ref_Overall_F_CalcNavCoal$CalcNavCoal != 6 & Ref_Overall_F_CalcNavCoal$CalcNavCoal != 9),]
Ref_Overall_F_CalcNavCoal_0 <- Ref_Overall_F_CalcNavCoal[1,2]
Ref_Overall_F_CalcNavCoal_0[is.na(Ref_Overall_F_CalcNavCoal_0)] <-0
Ref_Overall_F_CalcNavCoal_1 <- Ref_Overall_F_CalcNavCoal[2,2]
Ref_Overall_F_CalcNavCoal_1[is.na(Ref_Overall_F_CalcNavCoal_1)] <-0
Ref_Overall_F_CalcNavCoal_Total <- 
  Ref_Overall_F_CalcNavCoal_0 + Ref_Overall_F_CalcNavCoal_1
Ref_Overall_F_CalcNavCoal_Affected <- 
  round((Ref_Overall_F_CalcNavCoal_1/Ref_Overall_F_CalcNavCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_CalcNavCoal_row <- c("CalcNavCoal", Ref_Overall_CalcNavCoal_Total, 
                                 Ref_Overall_CalcNavCoal_1, 
                                 Ref_Overall_CalcNavCoal_Affected, 
                                 Ref_Overall_M_CalcNavCoal_Total, 
                                 Ref_Overall_M_CalcNavCoal_1, 
                                 Ref_Overall_M_CalcNavCoal_Affected, 
                                 Ref_Overall_F_CalcNavCoal_Total, 
                                 Ref_Overall_F_CalcNavCoal_1, 
                                 Ref_Overall_F_CalcNavCoal_Affected)

# AnkleCoal unclustered: TaloCalcCoal #
# Overall N and frequency #
Ref_Overall_TaloCalcCoal <- Ref_Overall %>% count(TaloCalcCoal)
Ref_Overall_TaloCalcCoal <- 
  Ref_Overall_TaloCalcCoal[(Ref_Overall_TaloCalcCoal$TaloCalcCoal != 6 & Ref_Overall_TaloCalcCoal$TaloCalcCoal != 9),]
Ref_Overall_TaloCalcCoal_0 <- Ref_Overall_TaloCalcCoal[1,2]
Ref_Overall_TaloCalcCoal_0[is.na(Ref_Overall_TaloCalcCoal_0)] <-0
Ref_Overall_TaloCalcCoal_1 <- Ref_Overall_TaloCalcCoal[2,2]
Ref_Overall_TaloCalcCoal_1[is.na(Ref_Overall_TaloCalcCoal_1)] <-0
Ref_Overall_TaloCalcCoal_Total <- 
  na.omit(Ref_Overall_TaloCalcCoal_0) + na.omit(Ref_Overall_TaloCalcCoal_1)
Ref_Overall_TaloCalcCoal_Affected <- 
  round((Ref_Overall_TaloCalcCoal_1/Ref_Overall_TaloCalcCoal_Total)*100, 2)

# Male N and frequency
Ref_Overall_M_TaloCalcCoal <- Ref_Overall_M %>% count(TaloCalcCoal)
Ref_Overall_M_TaloCalcCoal <- 
  Ref_Overall_M_TaloCalcCoal[(Ref_Overall_M_TaloCalcCoal$TaloCalcCoal != 6 & Ref_Overall_M_TaloCalcCoal$TaloCalcCoal != 9),]
Ref_Overall_M_TaloCalcCoal_0 <- Ref_Overall_M_TaloCalcCoal[1,2]
Ref_Overall_M_TaloCalcCoal_0[is.na(Ref_Overall_M_TaloCalcCoal_0)] <-0
Ref_Overall_M_TaloCalcCoal_1 <- Ref_Overall_M_TaloCalcCoal[2,2]
Ref_Overall_M_TaloCalcCoal_1[is.na(Ref_Overall_M_TaloCalcCoal_1)] <-0
Ref_Overall_M_TaloCalcCoal_Total <- 
  Ref_Overall_M_TaloCalcCoal_0 + Ref_Overall_M_TaloCalcCoal_1
Ref_Overall_M_TaloCalcCoal_Affected <- 
  round((Ref_Overall_M_TaloCalcCoal_1/Ref_Overall_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_TaloCalcCoal <- Ref_Overall_F %>% count(TaloCalcCoal)
Ref_Overall_F_TaloCalcCoal <- 
  Ref_Overall_F_TaloCalcCoal[(Ref_Overall_F_TaloCalcCoal$TaloCalcCoal != 6 & Ref_Overall_F_TaloCalcCoal$TaloCalcCoal != 9),]
Ref_Overall_F_TaloCalcCoal_0 <- Ref_Overall_F_TaloCalcCoal[1,2]
Ref_Overall_F_TaloCalcCoal_0[is.na(Ref_Overall_F_TaloCalcCoal_0)] <-0
Ref_Overall_F_TaloCalcCoal_1 <- Ref_Overall_F_TaloCalcCoal[2,2]
Ref_Overall_F_TaloCalcCoal_1[is.na(Ref_Overall_F_TaloCalcCoal_1)] <-0
Ref_Overall_F_TaloCalcCoal_Total <- 
  Ref_Overall_F_TaloCalcCoal_0 + Ref_Overall_F_TaloCalcCoal_1
Ref_Overall_F_TaloCalcCoal_Affected <- 
  round((Ref_Overall_F_TaloCalcCoal_1/Ref_Overall_F_TaloCalcCoal_Total)*100, 2)

# overview of calculated frequencies for table (see below) #
Ref_Overall_TaloCalcCoal_row <- c("TaloCalcCoal", Ref_Overall_TaloCalcCoal_Total, 
                                  Ref_Overall_TaloCalcCoal_1, 
                                  Ref_Overall_TaloCalcCoal_Affected, 
                                  Ref_Overall_M_TaloCalcCoal_Total, 
                                  Ref_Overall_M_TaloCalcCoal_1, 
                                  Ref_Overall_M_TaloCalcCoal_Affected, 
                                  Ref_Overall_F_TaloCalcCoal_Total, 
                                  Ref_Overall_F_TaloCalcCoal_1, 
                                  Ref_Overall_F_TaloCalcCoal_Affected)


## Table with overview of trait frequencies for Ref_Overall ##
# Table with overall frequencies #
Ref_Overall_freq <- data.frame(rbind(Ref_Overall_AccessNav_row, 
                                     Ref_Overall_BrachyD_row, 
                                     Ref_Overall_BrachyMT1_row,
                                     Ref_Overall_BrachyMT4_row, 
                                     Ref_Overall_BrachyPP1_row, 
                                     Ref_Overall_CalcCubCoal_row,
                                     Ref_Overall_TaloNavCoal_row, 
                                     Ref_Overall_CF2CF3Coal_row, 
                                     Ref_Overall_CF3MT3Coal_row,
                                     Ref_Overall_OsIntermet_row, 
                                     Ref_Overall_CF1Intermet_row, 
                                     Ref_Overall_MT1Intermet_row,
                                     Ref_Overall_MT2Intermet_row, 
                                     Ref_Overall_AnkleCoal_row, 
                                     Ref_Overall_CalcNavCoal_row, 
                                     Ref_Overall_TaloCalcCoal_row))
names(Ref_Overall_freq) <- c("Anomaly", "Overall N (Ref_Overall)", 
                             "Overall Affected (Ref_Overall)", 
                             "Overall % (Ref_Overall)", "Male N (Ref_Overall)", 
                             "Male Affected (Ref_Overall)", 
                             "Male % (Ref_Overall)", "Female N (Ref_Overall)", 
                             "Female affected (Ref_Overall)", 
                             "Female % (Ref_Overall)")


# Table with overview of trait frequencies for reference sample #
Ref_freq <- data.frame(cbind(Ref_Overall_freq, AR_freq[,c(2:10)], 
                             EH_freq[,c(2:10)], ZW_freq[,c(2:10)]))


# Comparing study sample and reference sample -----------------------------

#                                                                       #
## Comparing the study sample (MB) with the reference sample for males ##
#                                                                       #

## Fisher's exact test for AccessNav ##
# Contingency table #
Tab_M_AccessNav <- data.frame(
  "MB" = c(MB_M_AccessNav_0, MB_M_AccessNav_1),
  "Reference" = c(Ref_Overall_AccessNav_0, Ref_Overall_AccessNav_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_AccessNav <- stats::fisher.test(Tab_M_AccessNav, alternative = "two.sided")


## Fisher's exact test for BrachyD ##
# Contingency table #
Tab_M_BrachyD <- data.frame(
  "MB" = c(MB_M_BrachyD_0, MB_M_BrachyD_1),
  "Reference" = c(Ref_Overall_BrachyD_0, Ref_Overall_BrachyD_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_BrachyD <- stats::fisher.test(Tab_M_BrachyD)


## Fisher's exact test for BrachyMT1 ##
# Contingency table #
Tab_M_BrachyMT1 <- data.frame(
  "MB" = c(MB_M_BrachyMT1_0, MB_M_BrachyMT1_1),
  "Reference" = c(Ref_Overall_BrachyMT1_0, Ref_Overall_BrachyMT1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_BrachyMT1 <- stats::fisher.test(Tab_M_BrachyMT1)


## Fisher's exact test for BrachyMT4 ##
# Contingency table #
Tab_M_BrachyMT4 <- data.frame(
  "MB" = c(MB_M_BrachyMT4_0, MB_M_BrachyMT4_1),
  "Reference" = c(Ref_Overall_BrachyMT4_0, Ref_Overall_BrachyMT4_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_BrachyMT4 <- stats::fisher.test(Tab_M_BrachyMT4)


## Fisher's exact test for BrachyPP1 ##
# Contingency table #
Tab_M_BrachyPP1 <- data.frame(
  "MB" = c(MB_M_BrachyPP1_0, MB_M_BrachyPP1_1),
  "Reference" = c(Ref_Overall_BrachyPP1_0, Ref_Overall_BrachyPP1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_BrachyPP1 <- stats::fisher.test(Tab_M_BrachyPP1)


## Fisher's exact test for CalcCubCoal ##
# Contingency table #
Tab_M_CalcCubCoal <- data.frame(
  "MB" = c(MB_M_CalcCubCoal_0, MB_M_CalcCubCoal_1),
  "Reference" = c(Ref_Overall_CalcCubCoal_0, Ref_Overall_CalcCubCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_CalcCubCoal <- stats::fisher.test(Tab_M_CalcCubCoal)


## Fisher's exact test for TaloNavCoal ##
# Contingency table #
Tab_M_TaloNavCoal <- data.frame(
  "MB" = c(MB_M_TaloNavCoal_0, MB_M_TaloNavCoal_1),
  "Reference" = c(Ref_Overall_TaloNavCoal_0, Ref_Overall_TaloNavCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_TaloNavCoal <- stats::fisher.test(Tab_M_TaloNavCoal)


## Fisher's exact test for CF2CF3Coal ##
# Contingency table #
Tab_M_CF2CF3Coal <- data.frame(
  "MB" = c(MB_M_CF2CF3Coal_0, MB_M_CF2CF3Coal_1),
  "Reference" = c(Ref_Overall_CF2CF3Coal_0, Ref_Overall_CF2CF3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_CF2CF3Coal <- stats::fisher.test(Tab_M_CF2CF3Coal)


## Fisher's exact test for CF3MT3Coal ##
# Contingency table #
Tab_M_CF3MT3Coal <- data.frame(
  "MB" = c(MB_M_CF3MT3Coal_0, MB_M_CF3MT3Coal_1),
  "Reference" = c(Ref_Overall_CF3MT3Coal_0, Ref_Overall_CF3MT3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_CF3MT3Coal <- stats::fisher.test(Tab_M_CF3MT3Coal)


## Fisher's exact test for OsIntermet ##
# Contingency table #
Tab_M_OsIntermet <- data.frame(
  "MB" = c(MB_M_OsIntermet_0, MB_M_OsIntermet_1),
  "Reference" = c(Ref_Overall_OsIntermet_0, Ref_Overall_OsIntermet_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_OsIntermet <- stats::fisher.test(Tab_M_OsIntermet)


## Fisher's exact test for AnkleCoal ##
# Contingency table #
Tab_M_AnkleCoal <- data.frame(
  "MB" = c(MB_M_AnkleCoal_0, MB_M_AnkleCoal_1),
  "Reference" = c(Ref_Overall_AnkleCoal_0, Ref_Overall_AnkleCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
M_AnkleCoal <- stats::fisher.test(Tab_M_AnkleCoal)


## Overview of comparison values for males in a table ##
# Converting p.values into a single data frame #
M_pvalues <- data.frame(
  "p.value" = c(round(M_AccessNav$p.value,3), round(M_BrachyD$p.value,3), 
                round(M_BrachyMT1$p.value,3), round(M_BrachyMT4$p.value,3), 
                round(M_BrachyPP1$p.value,3), round(M_CalcCubCoal$p.value, 3),
                round(M_TaloNavCoal$p.value, 3), round(M_CF2CF3Coal$p.value, 3), 
                round(M_CF3MT3Coal$p.value, 3), round(M_OsIntermet$p.value, 3), 
                round(M_AnkleCoal$p.value, 3)),
  row.names = c("AccessNav", "BrachyD", "BrachyMT1", "BrachyMT4", 
                "BrachyPP1", "CalcCubCoal", "TaloNavCoal", "CF2CF3Coal", 
                "CF3MT3Intermet", "OsIntermet", "AnkleCoal"),
  stringsAsFactors = FALSE)
# combining frequency tables and p.values for males in one data frame #
Freq_M <- data.frame(cbind(Ref_Overall_freq[c(1:10, 14),c(1:4)], 
                           MB_freq[c(1:10, 14),c(5:7)], M_pvalues))
# renaming into significant names #
Freq_M <- Freq_M %>% rename(
  Ref.N = Overall.N..Ref_Overall., 
  Ref.Aff = Overall.Affected..Ref_Overall.,
  Ref.perc = Overall....Ref_Overall.,
  MB.Male.N = Male.N, 
  MB.Male.Aff = Male.affected,
  MB.Male.perc = Male..)


#                                                                         #
## Comparing the study sample (MB) with the reference sample for females ##
#                                                                         #


## Fisher's exact test for AccessNav ##
# Contingency table #
Tab_F_AccessNav <- data.frame(
  "MB" = c(MB_F_AccessNav_0, MB_F_AccessNav_1),
  "Reference" = c(Ref_Overall_AccessNav_0, Ref_Overall_AccessNav_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_AccessNav <- stats::fisher.test(Tab_F_AccessNav)
Tab_M_AccessNav

## Fisher's exact test for BrachyD ##
# Contingency table #
Tab_F_BrachyD <- data.frame(
  "MB" = c(MB_F_BrachyD_0, MB_F_BrachyD_1),
  "Reference" = c(Ref_Overall_BrachyD_0, Ref_Overall_BrachyD_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_BrachyD <- stats::fisher.test(Tab_F_BrachyD)


## Fisher's exact test for BrachyMT1 ##
# Contingency table #
Tab_F_BrachyMT1 <- data.frame(
  "MB" = c(MB_F_BrachyMT1_0, MB_F_BrachyMT1_1),
  "Reference" = c(Ref_Overall_BrachyMT1_0, Ref_Overall_BrachyMT1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_BrachyMT1 <- stats::fisher.test(Tab_F_BrachyMT1)


## Fisher's exact test for BrachyMT4 ##
# Contingency table #
Tab_F_BrachyMT4 <- data.frame(
  "MB" = c(MB_F_BrachyMT4_0, MB_F_BrachyMT4_1),
  "Reference" = c(Ref_Overall_BrachyMT4_0, Ref_Overall_BrachyMT4_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_BrachyMT4 <- stats::fisher.test(Tab_F_BrachyMT4)


## Fisher's exact test for BrachyPP1 ##
# Contingency table #
Tab_F_BrachyPP1 <- data.frame(
  "MB" = c(MB_F_BrachyPP1_0, MB_F_BrachyPP1_1),
  "Reference" = c(Ref_Overall_BrachyPP1_0, Ref_Overall_BrachyPP1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_BrachyPP1 <- stats::fisher.test(Tab_F_BrachyPP1)


## Fisher's exact test for CalcCubCoal ##
# Contingency table #
Tab_F_CalcCubCoal <- data.frame(
  "MB" = c(MB_F_CalcCubCoal_0, MB_F_CalcCubCoal_1),
  "Reference" = c(Ref_Overall_CalcCubCoal_0, Ref_Overall_CalcCubCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_CalcCubCoal <- stats::fisher.test(Tab_F_CalcCubCoal)


## Fisher's exact test for TaloNavCoal ##
# Contingency table #
Tab_F_TaloNavCoal <- data.frame(
  "MB" = c(MB_F_TaloNavCoal_0, MB_F_TaloNavCoal_1),
  "Reference" = c(Ref_Overall_TaloNavCoal_0, Ref_Overall_TaloNavCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_TaloNavCoal <- stats::fisher.test(Tab_F_TaloNavCoal)


## Fisher's exact test for CF2CF3Coal ##
# Contingency table #
Tab_F_CF2CF3Coal <- data.frame(
  "MB" = c(MB_F_CF2CF3Coal_0, MB_F_CF2CF3Coal_1),
  "Reference" = c(Ref_Overall_CF2CF3Coal_0, Ref_Overall_CF2CF3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_CF2CF3Coal <- stats::fisher.test(Tab_F_CF2CF3Coal)


## Fisher's exact test for CF3MT3Coal ##
# Contingency table #
Tab_F_CF3MT3Coal <- data.frame(
  "MB" = c(MB_F_CF3MT3Coal_0, MB_F_CF3MT3Coal_1),
  "Reference" = c(Ref_Overall_CF3MT3Coal_0, Ref_Overall_CF3MT3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_CF3MT3Coal <- stats::fisher.test(Tab_F_CF3MT3Coal)


## Fisher's exact test for OsIntermet ##
# Contingency table #
Tab_F_OsIntermet <- data.frame(
  "MB" = c(MB_F_OsIntermet_0, MB_F_OsIntermet_1),
  "Reference" = c(Ref_Overall_OsIntermet_0, Ref_Overall_OsIntermet_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_OsIntermet <- stats::fisher.test(Tab_F_OsIntermet)


## Fisher's exact test for AnkleCoal ##
# Contingency table #
Tab_F_AnkleCoal <- data.frame(
  "MB" = c(MB_F_AnkleCoal_0, MB_F_AnkleCoal_1),
  "Reference" = c(Ref_Overall_AnkleCoal_0, Ref_Overall_AnkleCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
F_AnkleCoal <- stats::fisher.test(Tab_F_AnkleCoal)



# converting p.values into a single data frame #
F_pvalues <- data.frame(
  "p.value" = c(round(F_AccessNav$p.value,3), round(F_BrachyD$p.value,3), 
                round(F_BrachyMT1$p.value,3), round(F_BrachyMT4$p.value,3),
                round(F_BrachyPP1$p.value,3), round(F_CalcCubCoal$p.value, 3),
                round(F_TaloNavCoal$p.value, 3), round(F_CF2CF3Coal$p.value, 3), 
                round(F_CF3MT3Coal$p.value, 3), round(F_OsIntermet$p.value, 3),
                round(F_AnkleCoal$p.value, 3)),
  row.names = c("AccessNav", "BrachyD", "BrachyMT1", "BrachyMT4", 
                "BrachyPP1", "CalcCubCoal", "TaloNavCoal", "CF2CF3Coal", 
                "CF3MT3Intermet", "OsIntermet", "AnkleCoal"),
  stringsAsFactors = FALSE)
# combining frequency tables and p.values for females in one data frame #
Freq_F <- data.frame(cbind(Ref_Overall_freq[c(1:10, 14),c(1:4)], 
                           MB_freq[c(1:10, 14),c(8:10)], F_pvalues))
# renaming into significant names
Freq_F <- Freq_F %>% rename(
  Ref.N = Overall.N..Ref_Overall., 
  Ref.Aff = Overall.Affected..Ref_Overall.,
  Ref.perc = Overall....Ref_Overall.,
  MB.Female.N = Female.N, 
  MB.Female.Aff = Female.affected,
  MB.Female.perc = Female..)


# Comparing complete reference sample with MB -----------------------------

## Fisher's exact test for AccessNav ##
# Contingency table #
Tab_AccessNav <- data.frame(
  "MB" = c(MB_AccessNav_0, MB_AccessNav_1),
  "Reference" = c(Ref_Overall_AccessNav_0, Ref_Overall_AccessNav_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_AccessNav <- stats::fisher.test(Tab_AccessNav)


## Fisher's exact test for BrachyD ##
# Contingency table #
Tab_BrachyD <- data.frame(
  "MB" = c(MB_BrachyD_0, MB_BrachyD_1),
  "Reference" = c(Ref_Overall_BrachyD_0, Ref_Overall_BrachyD_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_BrachyD <- stats::fisher.test(Tab_BrachyD)


## Fisher's exact test for BrachyMT1 ##
# Contingency table #
Tab_BrachyMT1 <- data.frame(
  "MB" = c(MB_BrachyMT1_0, MB_BrachyMT1_1),
  "Reference" = c(Ref_Overall_BrachyMT1_0, Ref_Overall_BrachyMT1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_BrachyMT1 <- stats::fisher.test(Tab_BrachyMT1)


## Fisher's exact test for BrachyMT4 ##
# Contingency table #
Tab_BrachyMT4 <- data.frame(
  "MB" = c(MB_BrachyMT4_0, MB_BrachyMT4_1),
  "Reference" = c(Ref_Overall_BrachyMT4_0, Ref_Overall_BrachyMT4_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_BrachyMT4 <- stats::fisher.test(Tab_BrachyMT4)


## Fisher's exact test for BrachyPP1 ##
# Contingency table #
Tab_BrachyPP1 <- data.frame(
  "MB" = c(MB_BrachyPP1_0, MB_BrachyPP1_1),
  "Reference" = c(Ref_Overall_BrachyPP1_0, Ref_Overall_BrachyPP1_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_BrachyPP1 <- stats::fisher.test(Tab_BrachyPP1)


## Fisher's exact test for CalcCubCoal ##
# Contingency table #
Tab_CalcCubCoal <- data.frame(
  "MB" = c(MB_CalcCubCoal_0, MB_CalcCubCoal_1),
  "Reference" = c(Ref_Overall_CalcCubCoal_0, Ref_Overall_CalcCubCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_CalcCubCoal <- stats::fisher.test(Tab_CalcCubCoal)


## Fisher's exact test for TaloNavCoal ##
# Contingency table #
Tab_TaloNavCoal <- data.frame(
  "MB" = c(MB_TaloNavCoal_0, MB_TaloNavCoal_1),
  "Reference" = c(Ref_Overall_TaloNavCoal_0, Ref_Overall_TaloNavCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_TaloNavCoal <- stats::fisher.test(Tab_TaloNavCoal)


## Fisher's exact test for CF2CF3Coal ##
# Contingency table #
Tab_CF2CF3Coal <- data.frame(
  "MB" = c(MB_CF2CF3Coal_0, MB_CF2CF3Coal_1),
  "Reference" = c(Ref_Overall_CF2CF3Coal_0, Ref_Overall_CF2CF3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_CF2CF3Coal <- stats::fisher.test(Tab_CF2CF3Coal)


## Fisher's exact test for CF3MT3Coal ##
# Contingency table #
Tab_CF3MT3Coal <- data.frame(
  "MB" = c(MB_CF3MT3Coal_0, MB_CF3MT3Coal_1),
  "Reference" = c(Ref_Overall_CF3MT3Coal_0, Ref_Overall_CF3MT3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_CF3MT3Coal <- stats::fisher.test(Tab_CF3MT3Coal)


## Fisher's exact test for OsIntermet ##
# Contingency table #
Tab_OsIntermet <- data.frame(
  "MB" = c(MB_OsIntermet_0, MB_OsIntermet_1),
  "Reference" = c(Ref_Overall_OsIntermet_0, Ref_Overall_OsIntermet_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_OsIntermet <- stats::fisher.test(Tab_OsIntermet)


## Fisher's exact test for AnkleCoal ##
# Contingency table #
Tab_AnkleCoal <- data.frame(
  "MB" = c(MB_AnkleCoal_0, MB_AnkleCoal_1),
  "Reference" = c(Ref_Overall_AnkleCoal_0, Ref_Overall_AnkleCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_AnkleCoal <- stats::fisher.test(Tab_AnkleCoal)



# converting p.values into a single data frame #
Ref_pvalues <- data.frame(
  "p.value" = c(round(Ref_AccessNav$p.value,3), round(Ref_BrachyD$p.value,3), 
                round(Ref_BrachyMT1$p.value,3), round(Ref_BrachyMT4$p.value,3), 
                round(Ref_BrachyPP1$p.value,3), round(Ref_CalcCubCoal$p.value, 3),
                round(Ref_TaloNavCoal$p.value, 3), round(Ref_CF2CF3Coal$p.value, 3), 
                round(Ref_CF3MT3Coal$p.value, 3), round(Ref_OsIntermet$p.value, 3), 
                round(Ref_AnkleCoal$p.value, 3)),
  row.names = c("AccessNav", "BrachyD", "BrachyMT1", "BrachyMT4", 
                "BrachyPP1", "CalcCubCoal", "TaloNavCoal", "CF2CF3Coal", 
                "CF3MT3Intermet", "OsIntermet", "AnkleCoal"),
  stringsAsFactors = FALSE)

# combining frequency tables and p.values in one data frame #
Freq_Ref <- data.frame(cbind(Ref_Overall_freq[c(1:10, 14),c(1:4)], 
                             MB_freq[c(1:10, 14),c(2:4)], Ref_pvalues))
# renaming into significant names
Freq_Ref <- Freq_Ref %>% rename(
  Ref.N = Overall.N..Ref_Overall., 
  Ref.Aff = Overall.Affected..Ref_Overall.,
  Ref.perc = Overall....Ref_Overall., 
  MB.N = Overall.N, 
  MB.Aff = Overall.Affected,
  MB.perc = Overall..)


# Export of created files ----------------------------------------------


# Frequency table of MB #
write.csv(MB_freq, "MB_freq.csv", row.names = T)
# Frequency table of reference sample
write.csv(Ref_freq, "Ref_freq.csv", row.names = T)
# Comparison table males #
write.csv(Freq_M, "Freq_M.csv", row.names = T)
# Comparison tables females #
write.csv(Freq_F, "Freq_F.csv", row.names = T)
# Comparison tables overall reference #
write.csv(Freq_Ref, "Freq_Ref.csv", row.names = T)


# Differences in frequencies between sexes --------------------------------

#                              #
##  Middenbeemster collection ##
#                              #

### The examined traits are limited to the present ones in the previous steps ###

## Fisher's exact test for AccessNav ##
# Contingency table #
Tab_AccessNav_sex <- data.frame(
  "MB_M" = c(MB_M_AccessNav_0, MB_M_AccessNav_1),
  "MB_F" = c(MB_F_AccessNav_0, MB_F_AccessNav_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_AccessNav_sex <- fisher.test(Tab_AccessNav_sex)


## Fisher's exact test for CF3MT3Coal ##
# Contingency table #
Tab_CF3MT3Coal_sex <- data.frame(
  "MB_M" = c(MB_M_CF3MT3Coal_0, MB_M_CF3MT3Coal_1),
  "MB_F" = c(MB_F_CF3MT3Coal_0, MB_F_CF3MT3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_CF3MT3Coal_sex <- fisher.test(Tab_CF3MT3Coal_sex)


## Fisher's exact test for OsIntermet ##
# Contingency table #
Tab_OsIntermet_sex <- data.frame(
  "MB_M" = c(MB_M_OsIntermet_0, MB_M_OsIntermet_1),
  "MB_F" = c(MB_F_OsIntermet_0, MB_F_OsIntermet_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_OsIntermet_sex <- fisher.test(Tab_OsIntermet_sex)

## Fisher's exact test for AnkleCoal ##
# Contingency table #
Tab_AnkleCoal_sex <- data.frame(
  "MB_M" = c(MB_M_AnkleCoal_0, MB_M_AnkleCoal_1),
  "MB_F" = c(MB_F_AnkleCoal_0, MB_F_AnkleCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_AnkleCoal_sex <- fisher.test(Tab_AnkleCoal_sex)

MB_sex_diff <- data.frame(
  "AccessNav" = MB_AccessNav_sex$p.value,
  "CF3MT3Coal" = MB_CF3MT3Coal_sex$p.value,
  "OsIntermet" = MB_OsIntermet_sex$p.value,
  "AnkleCoal" = MB_AnkleCoal_sex$p.value,
  row.names = "MB sex p-values")
print(MB_sex_diff)


#                     #
##  Reference sample ##
#                     #

### The examined traits are limited to the present ones in the previous steps ###

## Fisher's exact test for AccessNav ##
# Contingency table #
Tab_AccessNav_sex_ref <- data.frame(
  "Ref_Overall_M" = c(Ref_Overall_M_AccessNav_0, Ref_Overall_M_AccessNav_1),
  "Ref_Overall_F" = c(Ref_Overall_F_AccessNav_0, Ref_Overall_F_AccessNav_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_Overall_AccessNav_sex <- fisher.test(Tab_AccessNav_sex_ref)


## Fisher's exact test for CF3MT3Coal ##
# Contingency table #
Tab_CF3MT3Coal_sex_ref <- data.frame(
  "Ref_Overall_M" = c(Ref_Overall_M_CF3MT3Coal_0, Ref_Overall_M_CF3MT3Coal_1),
  "Ref_Overall_F" = c(Ref_Overall_F_CF3MT3Coal_0, Ref_Overall_F_CF3MT3Coal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_Overall_CF3MT3Coal_sex <- fisher.test(Tab_CF3MT3Coal_sex_ref)


## Fisher's exact test for OsIntermet ##
# Contingency table #
Tab_OsIntermet_sex_ref <- data.frame(
  "Ref_Overall_M" = c(Ref_Overall_M_OsIntermet_0, Ref_Overall_M_OsIntermet_1),
  "Ref_Overall_F" = c(Ref_Overall_F_OsIntermet_0, Ref_Overall_F_OsIntermet_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_Overall_OsIntermet_sex <- fisher.test(Tab_OsIntermet_sex_ref)

## Fisher's exact test for AnkleCoal ##
# Contingency table #
Tab_AnkleCoal_sex_ref <- data.frame(
  "Ref_Overall_M" = c(Ref_Overall_M_AnkleCoal_0, Ref_Overall_M_AnkleCoal_1),
  "Ref_Overall_F" = c(Ref_Overall_F_AnkleCoal_0, Ref_Overall_F_AnkleCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
Ref_Overall_AnkleCoal_sex <- fisher.test(Tab_AnkleCoal_sex_ref)

Ref_Overall_sex_diff <- data.frame(
  "AccessNav" = Ref_Overall_AccessNav_sex$p.value,
  "CF3MT3Coal" = Ref_Overall_CF3MT3Coal_sex$p.value,
  "OsIntermet" = Ref_Overall_OsIntermet_sex$p.value,
  "AnkleCoal" = Ref_Overall_AnkleCoal_sex$p.value,
  row.names = "Ref_Overall sex p-values")
print(Ref_Overall_sex_diff)



# Biases in kinship group -------------------------------------------------

## Starting from before clustering ##


#                          #
## Biases for CalcNavCoal ##
#                          #


## CalcNavCoal in MB ##
# Male N and frequency
MB_M_CalcNavCoal <- MB_M %>% count(CalcNavCoal)
MB_M_CalcNavCoal <- 
  MB_M_CalcNavCoal[(MB_M_CalcNavCoal$CalcNavCoal != 6 & MB_M_CalcNavCoal$CalcNavCoal != 9),]
MB_M_CalcNavCoal_0 <- MB_M_CalcNavCoal[1,2]
MB_M_CalcNavCoal_0[is.na(MB_M_CalcNavCoal_0)] <-0
MB_M_CalcNavCoal_1 <- MB_M_CalcNavCoal[2,2]
MB_M_CalcNavCoal_1[is.na(MB_M_CalcNavCoal_1)] <-0


# Female N and frequency
MB_F_CalcNavCoal <- MB_F %>% count(CalcNavCoal)
MB_F_CalcNavCoal <- 
  MB_F_CalcNavCoal[(MB_F_CalcNavCoal$CalcNavCoal != 6 & MB_F_CalcNavCoal$CalcNavCoal != 9),]
MB_F_CalcNavCoal_0 <- MB_F_CalcNavCoal[1,2]
MB_F_CalcNavCoal_0[is.na(MB_F_CalcNavCoal_0)] <-0
MB_F_CalcNavCoal_1 <- MB_F_CalcNavCoal[2,2]
MB_F_CalcNavCoal_1[is.na(MB_F_CalcNavCoal_1)] <-0

## Fisher's exact test for CalcNavCoal ##
# Contingency table #
Tab_CalcNavCoal_sex <- data.frame(
  "MB_M" = c(MB_M_CalcNavCoal_0, MB_M_CalcNavCoal_1),
  "MB_F" = c(MB_F_CalcNavCoal_0, MB_F_CalcNavCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_CalcNavCoal_sex <- fisher.test(Tab_CalcNavCoal_sex)
print(MB_CalcNavCoal_sex)


## CalcNavCoal for full sample ##

# Male N and frequency
Ref_Overall_M_CalcNavCoal <- Ref_Overall_M %>% count(CalcNavCoal)
Ref_Overall_M_CalcNavCoal <- 
  Ref_Overall_M_CalcNavCoal[(Ref_Overall_M_CalcNavCoal$CalcNavCoal != 6 & Ref_Overall_M_CalcNavCoal$CalcNavCoal != 9),]
Ref_Overall_M_CalcNavCoal_0 <- Ref_Overall_M_CalcNavCoal[1,2]
Ref_Overall_M_CalcNavCoal_0[is.na(Ref_Overall_M_CalcNavCoal_0)] <-0
Ref_Overall_M_CalcNavCoal_1 <- Ref_Overall_M_CalcNavCoal[2,2]
Ref_Overall_M_CalcNavCoal_1[is.na(Ref_Overall_M_CalcNavCoal_1)] <-0
Ref_Overall_M_CalcNavCoal_Total <- 
  Ref_Overall_M_CalcNavCoal_0 + Ref_Overall_M_CalcNavCoal_1
Ref_Overall_M_CalcNavCoal_Affected <- 
  round((Ref_Overall_M_CalcNavCoal_1/Ref_Overall_M_CalcNavCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_CalcNavCoal <- Ref_Overall_F %>% count(CalcNavCoal)
Ref_Overall_F_CalcNavCoal <- 
  Ref_Overall_F_CalcNavCoal[(Ref_Overall_F_CalcNavCoal$CalcNavCoal != 6 & Ref_Overall_F_CalcNavCoal$CalcNavCoal != 9),]
Ref_Overall_F_CalcNavCoal_0 <- Ref_Overall_F_CalcNavCoal[1,2]
Ref_Overall_F_CalcNavCoal_0[is.na(Ref_Overall_F_CalcNavCoal_0)] <-0
Ref_Overall_F_CalcNavCoal_1 <- Ref_Overall_F_CalcNavCoal[2,2]
Ref_Overall_F_CalcNavCoal_1[is.na(Ref_Overall_F_CalcNavCoal_1)] <-0
Ref_Overall_F_CalcNavCoal_Total <- 
  Ref_Overall_F_CalcNavCoal_0 + Ref_Overall_F_CalcNavCoal_1
Ref_Overall_F_CalcNavCoal_Affected <- 
  round((Ref_Overall_F_CalcNavCoal_1/Ref_Overall_F_CalcNavCoal_Total)*100, 2)


## Fisher's exact test for CalcNavCoal ##
# Contingency table #
Tab_CalcNavCoal_sex_overall <- data.frame(
  "Overall_M" = c((MB_M_CalcNavCoal_0 + Ref_Overall_M_CalcNavCoal_0), 
             (MB_M_CalcNavCoal_1 + Ref_Overall_M_CalcNavCoal_1)),
  "Overall_F" = c((MB_F_CalcNavCoal_0 + Ref_Overall_F_CalcNavCoal_0), 
             (MB_F_CalcNavCoal_1 + Ref_Overall_F_CalcNavCoal_1)),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
CalcNavCoal_sex_overall <- fisher.test(Tab_CalcNavCoal_sex_overall)
print(CalcNavCoal_sex_overall)


#                           #
## Biases for TaloCalcCoal ##
#                           #

## TaloCalcCoal for MB ##
# Male N and frequency
MB_M_TaloCalcCoal <- MB_M %>% count(TaloCalcCoal)
MB_M_TaloCalcCoal <- 
  MB_M_TaloCalcCoal[(MB_M_TaloCalcCoal$TaloCalcCoal != 6 & MB_M_TaloCalcCoal$TaloCalcCoal != 9),]
MB_M_TaloCalcCoal_0 <- MB_M_TaloCalcCoal[1,2]
MB_M_TaloCalcCoal_0[is.na(MB_M_TaloCalcCoal_0)] <-0
MB_M_TaloCalcCoal_1 <- MB_M_TaloCalcCoal[2,2]
MB_M_TaloCalcCoal_1[is.na(MB_M_TaloCalcCoal_1)] <-0


# Female N and frequency
MB_F_TaloCalcCoal <- MB_F %>% count(TaloCalcCoal)
MB_F_TaloCalcCoal <- 
  MB_F_TaloCalcCoal[(MB_F_TaloCalcCoal$TaloCalcCoal != 6 & MB_F_TaloCalcCoal$TaloCalcCoal != 9),]
MB_F_TaloCalcCoal_0 <- MB_F_TaloCalcCoal[1,2]
MB_F_TaloCalcCoal_0[is.na(MB_F_TaloCalcCoal_0)] <-0
MB_F_TaloCalcCoal_1 <- MB_F_TaloCalcCoal[2,2]
MB_F_TaloCalcCoal_1[is.na(MB_F_TaloCalcCoal_1)] <-0


## Fisher's exact test for TaloCalcCoal ##
# Contingency table #
Tab_TaloCalcCoal_sex <- data.frame(
  "MB_M" = c(MB_M_TaloCalcCoal_0, MB_M_TaloCalcCoal_1),
  "MB_F" = c(MB_F_TaloCalcCoal_0, MB_F_TaloCalcCoal_1),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
MB_TaloCalcCoal_sex <- fisher.test(Tab_TaloCalcCoal_sex)
print(MB_TaloCalcCoal_sex)



## TaloCalcCoal for full sample ##

# Male N and frequency
Ref_Overall_M_TaloCalcCoal <- Ref_Overall_M %>% count(TaloCalcCoal)
Ref_Overall_M_TaloCalcCoal <- 
  Ref_Overall_M_TaloCalcCoal[(Ref_Overall_M_TaloCalcCoal$TaloCalcCoal != 6 & Ref_Overall_M_TaloCalcCoal$TaloCalcCoal != 9),]
Ref_Overall_M_TaloCalcCoal_0 <- Ref_Overall_M_TaloCalcCoal[1,2]
Ref_Overall_M_TaloCalcCoal_0[is.na(Ref_Overall_M_TaloCalcCoal_0)] <-0
Ref_Overall_M_TaloCalcCoal_1 <- Ref_Overall_M_TaloCalcCoal[2,2]
Ref_Overall_M_TaloCalcCoal_1[is.na(Ref_Overall_M_TaloCalcCoal_1)] <-0
Ref_Overall_M_TaloCalcCoal_Total <- 
  Ref_Overall_M_TaloCalcCoal_0 + Ref_Overall_M_TaloCalcCoal_1
Ref_Overall_M_TaloCalcCoal_Affected <- 
  round((Ref_Overall_M_TaloCalcCoal_1/Ref_Overall_M_TaloCalcCoal_Total)*100, 2)

# Female N and frequency
Ref_Overall_F_TaloCalcCoal <- Ref_Overall_F %>% count(TaloCalcCoal)
Ref_Overall_F_TaloCalcCoal <- 
  Ref_Overall_F_TaloCalcCoal[(Ref_Overall_F_TaloCalcCoal$TaloCalcCoal != 6 & Ref_Overall_F_TaloCalcCoal$TaloCalcCoal != 9),]
Ref_Overall_F_TaloCalcCoal_0 <- Ref_Overall_F_TaloCalcCoal[1,2]
Ref_Overall_F_TaloCalcCoal_0[is.na(Ref_Overall_F_TaloCalcCoal_0)] <-0
Ref_Overall_F_TaloCalcCoal_1 <- Ref_Overall_F_TaloCalcCoal[2,2]
Ref_Overall_F_TaloCalcCoal_1[is.na(Ref_Overall_F_TaloCalcCoal_1)] <-0
Ref_Overall_F_TaloCalcCoal_Total <- 
  Ref_Overall_F_TaloCalcCoal_0 + Ref_Overall_F_TaloCalcCoal_1
Ref_Overall_F_TaloCalcCoal_Affected <- 
  round((Ref_Overall_F_TaloCalcCoal_1/Ref_Overall_F_TaloCalcCoal_Total)*100, 2)


## Fisher's exact test for TaloCalcCoal ##
# Contingency table #
Tab_TaloCalcCoal_sex_overall <- data.frame(
  "Overall_M" = c((MB_M_TaloCalcCoal_0 + Ref_Overall_M_TaloCalcCoal_0), 
             (MB_M_TaloCalcCoal_1 + Ref_Overall_M_TaloCalcCoal_1)),
  "Overall_F" = c((MB_F_TaloCalcCoal_0 + Ref_Overall_F_TaloCalcCoal_0), 
             (MB_F_TaloCalcCoal_1 + Ref_Overall_F_TaloCalcCoal_1)),
  row.names = c("Absent", "Present"),
  stringsAsFactors = FALSE
)
# Fisher's exact test with significance level 95% (alpha 0.05) #
TaloCalcCoal_sex_overall <- fisher.test(Tab_TaloCalcCoal_sex_overall)
print(TaloCalcCoal_sex_overall)



# Spatial distribution ----------------------------------------------------

### In the following section, the spatial distribution of the graves of the  ###
### hypothetical kinship group formed by the Middenbeemster individuals with ###
### the ankle coalition traits will be analysed by spatial statistics.       ###

## Preparing the dataset for the following statistics ##
  # Creating dataframe with columns: trait, x-coordinate and y-coordinate #
MB_AnkleCoal_spatial <- MB_wo_probables[, c("AnkleCoal", "x.coord", "y.coord")]
  # Removing individuals with indiscernible (6) or missing (9) values #
MB_AnkleCoal_spatial <- 
  MB_AnkleCoal_spatial[MB_AnkleCoal_spatial$AnkleCoal != 6 & MB_AnkleCoal_spatial$AnkleCoal != 9,]
  # Removing rows with no coordinates #
MB_AnkleCoal_spatial <- na.omit(MB_AnkleCoal_spatial)
  # Renaming the columns to enter in statistical analysis #
MB_AnkleCoal_spatial <- MB_AnkleCoal_spatial %>% rename(
  type = AnkleCoal,
  x = x.coord,
  y = y.coord
)



#                                                                       #
##  The following lines of codes were written by Keron (2015)          ##
### Reference:  Keron, J.R., 2015. The Use of Point Pattern Analysis  ###
### in Archaeology: Some Methods and Applications. Electronic Thesis  ###
##  and Dissertation Repository 3137.                                  ##
#                                                                       #


### Hodder and Okells A ###
#
# Note the input file must have 3 columns labeled x,y,and type all lower case
#
## Set variables ##
infilename <- "SPBF.txt"
typeA <- 1
typeB <- 0
randNum = 99
#
###---###---###---###
# Read the Input File
infile <- MB_AnkleCoal_spatial
n <- length(infile$type[infile$type == typeA])
x <- numeric(n)
y <- numeric(n)
j <- 1
for (i in 1:length(infile$x))
{
  if (infile$type[i]==typeA)
  {
    x[j] <- infile$x[i]
    y[j] <- infile$y[i]
    j <- j+1
  }
}
type1 <- data.frame(x,y)
m <- length(infile$type[infile$type == typeB])
x <- numeric(m)
y <- numeric(m)
j <- 1
for (i in 1:length(infile$x))
{
  if (infile$type[i]==typeB)
  {
    x[j] <- infile$x[i]
    y[j] <- infile$y[i]
    j <- j+1
  }
}
type2 <- data.frame(x,y)
###---###---###---###---###---###---###---###---###---###---###
# Calculate Hodder and Okells A of actual type 1 vs actual type2
aStat <- hoddersA(type1,type2)
xtitle <- "Hodders A "
xtitle
aStat
###---###---###---###---###---###---###---###---###---###---###
# Set up to run Monte Carlo on the at risk points.
resultsList <- numeric(randNum)
set.seed(14541)
for (j in 1:randNum) {
  myList <- randABSplit(infile,n)
  s1 <- myList$a
  s2 <- myList$b
  resultsList[j] <- hoddersA(s1,s2)
}
###---###---###---###---###---###---###---###---###---###---###
probLessThan <- pValueLessThan(resultsList, aStat)
probLessThan
# #
####End




### Proximity count ###
## Script Start ##
#
# Note that the input file must have 3 columns labeled "x, y, and type" (lower case)
#
### Set critical variables ###
#Set the Radius Parameter
mDist <- c(3, 5, 7)
#
# set the type to be counted
intype <- 1
#
# Set the input file name
infilename <- "PPOPFR.txt"
#
# set number of randomizations
randNum <- 999
###---###---###---###---###---###---###---###---###---###---###---
#
prob <- numeric(length(mDist))
actCount <- numeric(length(mDist))
for (p in 1:length(mDist))
{
  ###########################Regular oop#####################
  n <- length(infile$type[infile$type == intype])
  x <- numeric(n)
  y <- numeric(n)
  j <- 1
  for (i in 1:length(infile$x))
  {
    if (infile$type[i]==intype)
    {
      x[j] <- infile$x[i]
      y[j] <- infile$y[i]
      j <- j+1
    }
  }
  actuals <- data.frame(x,y)
  ###---###---###---###---###---###---###---###---###---###---###---
  # Calculate the actual count of points withing the specified distance
  #
  actCount[p] <- proxCount(actuals,mDist[p])
  ###---###---###---###---###---###---###---###---###---###---###---
  # Set up to run Monte Carlo on the at risk points
  #
  resultsList <- numeric(randNum)
  set.seed(14541)
  for (j in 1:randNum) {
    samplePoints <- randSampleXY(infile,n)
    resultsList[j] <- proxCount(samplePoints,mDist[p])
  }
  ###---###---###---###---###---###---###---###---###---###---###---
  prob[p] <- pValueGreaterThan(resultsList, actCount[p])
}
results <- data.frame(actCount,mDist,prob)
infilename
results
#
#
####End




###---###---###---###---###---###---###---###---###---###---###---
# Nearest Neighbour Random Label Script
###---###---###---###---###---###---###---###---###---###---###---
# This script calculates the average nearest neighbour for two given types.
#
#Note that the input file must have 3 columns labeled "x, y, and type" (lower case)
#
# type should have two values ;0' and presence '9'
MB_AnkleCoal_spatial$type[MB_AnkleCoal_spatial$type == 1] <- 9
# Note this implementation has not been generalized to allow multiple types ie 0,1 and 9)
#
# Note that all of these are case sensitive
###---###---###---###---###---###---###---###---###---###---###---
#
### Set key variables ###
# set the type to be counted Mostly 9
traitPresent <- 9
#
# Set the input file name
infilename <- "Sq.txt"
#
# set number of randomizations
randNum <- 999
###---###---###---###---###---###---###---###---###---###---###---
#
### select all males ###
infile <- MB_AnkleCoal_spatial
n <- length(infile$type[infile$type == 9])
x <- numeric(n)
y <- numeric(n)
j <- 1
for (i in 1:length(infile$x))
{
  if (infile$type[i]==9)
  {
    x[j] <- infile$x[i]
    
    y[j] <- infile$y[i]
    j <- j+1
  }
}
traitn <- data.frame(x,y)
###---###---###---###---###---###---###---###---###---###---###---
### ### ### ### ###
Title <- "AnkleCoal"
ActualAvgNN <- avgNNDist(traitn, traitn)
###---###---###---###---###---###---###---###---###---###---###---
# Set up to run Monte Carlo on the at risk points.
###---###---###---###---###---###---###---###---###---###---###---
set.seed(14541)
### Male to Male ###
resultsList <- numeric(randNum)
for (j in 1:randNum) {
  randList <- randSampleXY(infile,n)
  resultsList[j] <- avgNNDist(randList,randList)
}
RandomAvgNN <- mean(resultsList)
prob <- pValueLessThan(resultsList, ActualAvgNN)
nnratio <- ActualAvgNN / RandomAvgNN
#
###---###---###---###---###---###---###---###---###---###---###---
#
results <- data.frame(Title,traitPresent, ActualAvgNN,RandomAvgNN,nnratio,prob)
"Nearest Neighbour - Random Labeling"
"randomizations"
randNum
results
#
####End###

