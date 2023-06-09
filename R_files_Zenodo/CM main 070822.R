#############################################
##Classical Model R code#####################
#############################################

setwd("path of the working folder")

library(dplyr)          
library(tidyr)        
library(purrr)

# load other functions
source("R code/Functions/calibrationScore1.R")
source("R code/Functions/informationScore1.R")
source("R code/Functions/constructDM1.R")
source("R code/Functions/globalWeights_opt.R")

#specify the number of quantiles - for now this is done manually
Nquantiles = 5
quantilesPercent = c("5%","25%", "50%","75%","95%")
quantiles = c(0.05, 0.25, 0.50,0.75, 0.95)

#import realizations file, which contain info about the questions
realizations_file = read.csv("realizations.csv")
#get the vector of realizations
realizations = realizations_file[,2]
#get the list of questions
questions = realizations_file[,1]

#get state information about the variables
state_info_file = read.csv("Question_List.csv")

#create a list of lists with all the experts assessments
m = list() # create empty list of matrices for assessment data

# read experts assessments
badReads = 0
goodReads = 0
count_equal_quant = 0
equalQuantiles = 0

#select which data to choose
path = "path info"
#read the expert data - needs to have Exp in the name of the csv file
pattern = "Exp.*\\\\.csv$"      
csvFiles = list.files(path = path, pattern = pattern)
numbers = as.numeric(regmatches(csvFiles, regexpr("[0-9]+", csvFiles)))
csvFiles = csvFiles[order(numbers)]

expertMappingRev = array()

for (file in csvFiles) {
  tryCatch(
    {
      
      filename = sprintf("%s/%s", path, file)
      dataAll = read.csv(filename)
      
      # get assessments 
        data = dataAll[3:(Nquantiles+2)][,]
      # get expert ID
      expertID = as.integer(dataAll[1][1,])
      
      expertMappingRev[goodReads+1] = expertID
      
      # check if assessments are strictly increasing
      for (row in 1:NROW(data)) {
      
        if (any(diff(t(data[row,]))==0)) {
      
          count_equal_quant = count_equal_quant +1
          
          k = which(diff(t(data[row,]))==0)
          
          if (length(k)==1)
          {
            if (k==4)
            {
              data[row,5]= data[row,5]+0.0001
            }
           else { 
          if (data[row,k+1]+0.0001< data[row,k+2])
          {
            data[row,k+1] = data[row,k+1]+0.0001  
          }
          else if (data[row,k+1]+0.0001>= data[row,k+2])
          {
            data[row,k+1] = data[row,k+2]/10
          }
           }
          }
        
        else if (length(k)==2)#this is covering all the cases when three quantiles are equal
        {
          if (k[2]==4)
          {
            data[row,5]= data[row,5]+0.01
          }
          else{
          
          if (data[row,k[2]+1]+0.001< data[row,k[2]+2])
          {
            data[row,k[2]+1] = data[row,k[2]+1]+0.001  
          }
          else if (data[row,k[2]+1]+0.001>= data[row,k[2]+2])
          {
            data[row,k[2]+1] = data[row,k[2]+2]/10
          }
          }
          if (data[row,k[1]+1]+0.0001< data[row,k[1]+2])
          {
            data[row,k[1]+1] = data[row,k[1]+1]+0.0001  
          }
          else if (data[row,k[1]+1]+0.0001>= data[row,k[1]+2])
          {
            data[row,k[1]+1] = data[row,k[1]+2]/10
          }
        }
          
          else if (length(k)==3)#this is covering both the case when first three or last four quantiles are equal 
          {
            if (data[row,1]+0.0001>= data[row,2])
            {
              data[row,2] = data[row,1]+0.0001  
            }
            
            if (data[row,2]+0.001>= data[row,3])
            {
              data[row,3] = data[row,2]+0.001  
            }
            
            if (data[row,3]+0.01>= data[row,4])
            {
              data[row,4] = data[row,3]+0.01-0.001  
            }
            
            if (data[row,4]+0.1>= data[row,5])
            {
              data[row,5] = data[row,4]+0.1-0.01
            }
          }
          
          else if (length(k)==4)
          {
            data[row,2] = data[row,1]+0.0001
            data[row,3] = data[row,1]+0.001
            data[row,4] = data[row,1]+0.01
            data[row,5] = data[row,1]+0.1
          }
        }
        if (any(diff(t(data[row,]))<0)) {
          errorString = sprintf("Error assessments for %s question %d are strictly decreasing", filename, questions[row])
          stop(errorString)
        }
      }
      equalQuantiles = c(equalQuantiles, count_equal_quant)
      
      mm = as.matrix(data, nrow=NROW(data), ncol=Nquantiles)
      if (typeof(mm) == "double" || typeof(mm) == "integer") {
        m[[goodReads+1]] = mm
        goodReads = goodReads + 1
        count_equal_quant = 0
      } else {
        print("Error cannot interpret data as numeric")
        badReads = badReads + 1
      }
    },
    warning=function(cond) {
      print(cond)
      message("Warning cannot read csv file")
      badReads = badReads + 1
    }
  )
}

cat("Found", length(csvFiles), "csv files\\n")
numExperts = goodReads      # number of experts
cat("Read", numExperts, "experts\\n")
cat("Cannot read", badReads, "experts\\n")
cat("Equal quantiles", equalQuantiles[-1], "\\n")
N = NROW(m[[1]])   # total number of questions

# create labels for all experts
labelsExperts = c()
for (e in 1:numExperts) {
  labelsExperts = c(labelsExperts, sprintf("Expert %d", expertMappingRev[[e]]))
}

Nexperts = length(m)
Nstates =   length((unique(state_info_file$location)))

##########################################################################################################################################
#############performance of experts#######################################################################################################
##########################################################################################################################################
#################################
#store cal scores for all states#
#################################
calscores_all = matrix(0, nrow=Nstates, ncol=Nexperts+1)

for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
#  location_name = state_info_file$location_name[which(state_info_file$location==location)][1]
  ids = which(state_info_file[,3] == location)
  numb_cal = length(ids)
  calscores_all[i,] = c(numb_cal,calculateCalibrationScore(purrr::map(m, ~ .x[ids,]),realizations[ids]))
}

#############################
#info scores for all states##
#############################
infoscores_all = matrix(0, nrow=Nstates, ncol=Nexperts)

for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
  location_name = state_info_file$location_name[which(state_info_file$location==location)][1]
  ids = which(state_info_file[,3] == location)
  infoscores_all[i,] = colMeans(calculateInformationScore(purrr::map(m, ~ .x[ids,]),realizations[ids])[[1]][which(!is.na(realizations[ids])),]) 
}

##########################################################################################################################################
############costruction of DM ############################################################################################################
##########################################################################################################################################
# Equal Weights Decision Maker (EWDM)
# Performance Weights Decision Maker (PWDM)
equalWeights = 1/numExperts

#####################################################################
###performance scores of DMs#########################################
#####################################################################
calscores_all_DM = matrix(0, nrow=Nstates, ncol=2) #we're considering 2 DMs

for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
  ids = which(state_info_file[,3] == location)
  tmp = calculateInformationScore(purrr::map(m, ~ .x[ids,]),realizations[ids], k=0.1, bounds=NULL)
  L_i = tmp[[2]]
  U_i = tmp[[3]]
  
  EWDM_i = constructDM(purrr::map(m, ~ .x[ids,]),rep(equalWeights,numExperts), NULL, L_i, U_i, quantiles)
  
  perf_opt = perfWeights_opt(purrr::map(m, ~ .x[ids,]),realizations[ids])
  PWDM_opt_i = constructDM(purrr::map(m, ~ .x[ids,]), perf_opt, NULL, L_i, U_i, quantiles)
  
  calscores_all_DM[i,1] = calculateCalibrationScore(list(PWDM_opt_i),realizations[ids])
  calscores_all_DM[i,2] = calculateCalibrationScore(list(EWDM_i),realizations[ids])
}

colnames(calscores_all_DM) = c("GWDM opt","EWDM")


#################
##info score DMs#
#################
infoscores_all_DM = matrix(0, nrow=Nstates, ncol=2)

for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
  ids = which(state_info_file[,3] == location)
  tmp = calculateInformationScore(purrr::map(m, ~ .x[ids,]),realizations[ids], k=0.1, bounds=NULL)
  L_i = tmp[[2]]
  U_i = tmp[[3]]
  
  EWDM_i = constructDM(purrr::map(m, ~ .x[ids,]),rep(equalWeights,numExperts), NULL, L_i, U_i, quantiles)

  perf_opt = perfWeights_opt(purrr::map(m, ~ .x[ids,]),realizations[ids])
  PWDM_opt_i = constructDM(purrr::map(m, ~ .x[ids,]), perf_opt, NULL, L_i, U_i, quantiles)
  
  infoscores_all_DM[i,1] = colMeans(calculateInformationScore(list(PWDM_opt_i),realizations[ids])[[1]])
  infoscores_all_DM[i,2] = colMeans(calculateInformationScore(list(EWDM_i),realizations[ids])[[1]])
}

colnames(infoscores_all_DM) = c("GWDM opt","EWDM")

#################################################################################
########accuracy and precision###################################################
#################################################################################
accuracy <- function(my_data,my_realizations)
{
  accuracy_ratio = my_realizations/my_data
  
  accuracy = exp(mean(log(accuracy_ratio)))
  
  return(accuracy)
}


precision <- function(my_data,my_realizations)
{
  accuracy_ratio = my_realizations/my_data
  
  precision = exp(sd(log(accuracy_ratio)))
  
  return(precision)
}

###############################################
#models' performance in accuracy and precision#
###############################################
accuracy_models = matrix(0, nrow=Nstates, ncol=Nexperts)
for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
  ids = which(state_info_file[,3] == location)
  m_id = purrr::map(m, ~ .x[ids,])
  
  for (j in 1:Nexperts)
  {
    accuracy_models[i,j]= accuracy(m_id[[j]][,3],realizations[ids])
  }
}

precision_models = matrix(0, nrow=Nstates, ncol=Nexperts)
for (i in 1:Nstates)
{
  location = sort(unique(state_info_file$location))[i]
  ids = which(state_info_file[,3] == location)
  m_id = purrr::map(m, ~ .x[ids,])
  
  for (j in 1:Nexperts)
  {
    precision_models[i,j]= precision(m_id[[j]][,3],realizations[ids])
  }
}
