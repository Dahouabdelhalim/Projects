# VIRA Annual Turnover Analysis
#*** This script uses the Delta method to estimate the standard error for Virginia Rail (VIRA) colonization and extinction annual covariates graphed in Figure 3. ***
#    It uses the package msm to do the calculations and the script was written by Steve Beissinger.
#    All occupancy model estimates (parameter means and covariance matrices) were taken from models implemented in Program Presence, so are pasted into this script
#  The approach follows the same flow and command structure as in Deltamethod_annualcovars_BLRA.R

#Delta method to estimate the standard error for VIRA colonization and extinction covariates, so we could calculate the uncertainty and display it by adding 1SE or 1.96*SE to the mean values
#Results are from Present models in C:\\Users\\beis\\Box Sync\\Black_Rails\\Extinction Dynamics\\2021\\Analysis\\VIRA\\AnnCovs  on the Dell Desktop

#Patterned after scripts for BLRA in file "Deltamethod_turnovercovs.R in "C:/Users/beis/Box Sync/Black_Rails/Extinction Dynamics/2021/Analysis/Results")
# The code is in 4 parts with subsections (e.g., 3.1, 3.2, etc.)
# First section loads some packages. The second applies the Delta Method using the msm package to do the calculation of the standard error for back-transformed estimates
# for the model with only main effects.
# The third section applies it to the best VIRA Model with Main Effects for extinction  and a single two-way interaction between precipitation and West Nile virus index for extinction.
# All covariates were centered and standardized [(x - mean)/SD] and in the case of WNVvi (West Nile Virus vector index) were logged first (ln + 0.1).  
# The fourth section combines the variables into data files for graphing with another script.


#1.0 Load some packages
require(msm)  #does delta method 
library(dplyr)
library(readxl)
library(magrittr)
library(tidyr)

#Set Working Directory
setwd("")   


#*** 2.0 ALL PARAMETERS FROM THE MODEL WITH ONLY MAIN EFFECTS FOR PRECIPITATION ***

#2.1.1 Precip Colonization SE
covmat.colprecip = matrix(c( 0.006412, -0.000806, -0.000806,  0.006178), nrow=2, ncol=2, byrow=TRUE)
means.colprecip <- c(-1.867665,  0.301336)
stdvalue <- c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
              -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)
deltamethod(~ 1/(1+exp(-(x1 + x2*2))), means.colprecip, covmat.colprecip) 

#2.1.2 Need to write a function to pass stdvalue through it
# see "as.formula" explanation in the msm help file
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f)))", z)
form

#2.1.3 Write a loop to get it to run across all the stdvalues - this works
colprecip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  colprecip.se[[i]] <- deltamethod(as.formula(form[i]), means.colprecip, covmat.colprecip)
}


#2.1.4 Calculate Precip Extinction SE for each standard value from -2 to +2
covmat.extprecip = matrix(c(0.008747, 0.001757, 0.001757, 0.007951), nrow=2, ncol=2, byrow=TRUE)
means.extprecip <- c(-1.434815, -0.134830)    #slope and intercept are incredibly similar to BLRA values!
deltamethod(~ 1/(1+exp(-(x1 + x2))), means.extprecip, covmat.extprecip)  #works! 
#writing a loop to get it to run across all the stdvalues - this works
extprecip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  extprecip.se[[i]] <- deltamethod(as.formula(form[i]), means.extprecip, covmat.extprecip)
}

#2.1.5 Precip MEANS for col and ext 
colprecip.mean <-  1/(1+exp(-(-1.867665 + 0.301336*stdvalue)))   
extprecip.mean <-  1/(1+exp(-(-1.434815 + -0.134830*stdvalue)))   

#2.1.6 Now calculate Precip means + 1.96SE  or 1 SE across all standard values
colprecip.plus2se <- colprecip.mean + 1.96*colprecip.se
colprecip.minus2se <- colprecip.mean - 1.96*colprecip.se
extprecip.plus2se <- extprecip.mean + 1.96*extprecip.se
extprecip.minus2se <- extprecip.mean - 1.96*extprecip.se
colprecip.plus1se <- colprecip.mean + 1*colprecip.se
colprecip.minus1se <- colprecip.mean - 1*colprecip.se
extprecip.plus1se <- extprecip.mean + 1*extprecip.se
extprecip.minus1se <- extprecip.mean - 1*extprecip.se

#2.1.7 Create a dataframe and write the file if needed
VIRA_precip_all <- data.frame(stdvalue, colprecip.mean, colprecip.se, colprecip.plus2se, colprecip.minus2se, extprecip.mean, extprecip.se, extprecip.plus2se, extprecip.minus2se,
                         colprecip.plus1se, colprecip.minus1se, extprecip.plus1se, extprecip.minus1se)
str(VIRA_precip_all)
#saveRDS(VIRA_precip_all, file = "VIRA_precip_all.RDS")



#****** 3.0  BEST MODEL with two-way interaction for Colonization and main effects for Extinction. Do extinction first ********
#WNVcur is the vector index for the current year, unlike BLRA which uses the previous year

#3.1 Enter the untransformed Parameter values from the Best Model
InitialPsi <- -0.635065
Col.Intercept <-	-1.788575
Col.Precip <-	0.292327
Col.WNVcur <-	0.021578
Col.WNVcur_Precip <-	0.166967
Ext.Intercept <-	-1.445917
Ext.Precip <-	-0.053697
Ext.WNVcur <-	0.176981

# vector of standard values
stdvalue <- c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
              -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)

# 3.2 Estimating standard errors for Extinction 
# 3.2.1   deltamethod(~ 1/(1+exp(-(x1 + x2*0 + x3*stdvalue))), Means.Ext, Covmat.Ext) # precip is x2, WNV is x3, intercept is x1
Covmat.Ext <- matrix(c(0.009012,	0.001083,	-0.001393, 0.001083, 0.009106,	0.003635, -0.001393,	0.003635,	0.009999),nrow=3)
Means.Ext <- c(Ext.Intercept[1], Ext.Precip[1], Ext.WNVcur[1])  

# 3.2.2 Writing a loop to get it to run across all the stdvalues for Ext precip
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3*0)))", z)
form
Ext.Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.Precip.se[[i]] <- deltamethod(as.formula(form[i]), Means.Ext, Covmat.Ext)
}

# 3.2.3 writing a loop to get it to run across all the stdvalues for Ext WNV
form <- sprintf("~ 1/(1+exp(-(x1 + x2*0 + x3*%f)))", z)
form
Ext.WNVcur.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNVcur.se[[i]] <- deltamethod(as.formula(form[i]), Means.Ext, Covmat.Ext)
}

#3.3 Calculating Extinction means and SEs
Ext.Precip.mean <- 1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue)))
Ext.WNVcur.mean <- 1/(1+exp(-(Ext.Intercept + Ext.WNVcur*stdvalue)))
#Means plus 1 SE
Ext.Precip.plus1se <- Ext.Precip.mean + Ext.Precip.se
Ext.Precip.minus1se <- Ext.Precip.mean - Ext.Precip.se
Ext.WNVcur.plus1se <- Ext.WNVcur.mean + Ext.WNVcur.se
Ext.WNVcur.minus1se <- Ext.WNVcur.mean - Ext.WNVcur.se

#3.4 Putting the Extinction pieces into a dataframe
EXT_ALL <-  data.frame(stdvalue, Ext.Precip.mean, Ext.Precip.plus1se, Ext.Precip.minus1se, Ext.WNVcur.mean, 
                       Ext.WNVcur.plus1se, Ext.WNVcur.minus1se)
str(EXT_ALL)
#saveRDS(EXT_ALL, file = "VIRA_EXT_ALL.RDS")

#3.5 Colonization calculations with the interaction term
# Step 1 first put together the variance covariance matrix elements in covmat  and the mean estimates.

# 3.5.1 Estimating standard errors for Colonization
##    deltamethod(~ 1/(1+exp(-(x1 + x2*0 + x3*stdvalue + x4*stdvalue))), Means.Col, Covmat.Col) # precip is x2, WNV is x3, interaction is x4, intercept is x1
Covmat.Col <- matrix(c( 0.00691,	-0.002232,	-0.001852,	0.00214, -0.002232,	0.007153,	0.003242,	-0.001777,
                        -0.001852,	0.003242,	0.009019,	-0.002154, 0.00214,	-0.001777,	-0.002154,	0.006523),nrow=4)
Means.Col <- c(Col.Intercept[1], Col.Precip[1], Col.WNVcur[1], Col.WNVcur_Precip[1])  

#step 2 write the delta method statement
# (~ 1/(1+exp(-(intercept + precip*stdvalue + WNVcur + Interaction*stdvalue)))
##deltamethod(~ 1/(1+exp(-(x1 + x2*stdvalue + x3*0 + x4*stdvalue*0))), means.extwnv_precip, covmat.extwnv_precip) # for +1SD  so here WNV is the x2
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f*1 + x3*0 + x4*%f*0)))", z, z)
form
# step 3 run a loop to get it to run across all the stdvalues 
Col.Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.Precip.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}


#3.5.2 Do formula steps again for the +1SD for WNV, since this changes the WNV and interaction terms, or multiply each by 1
##deltamethod(~ 1/(1+exp(-(x1 - x2*stdvalue + x3*1 - x4*stdvalue*1))),  # for +1SD (last term), precip is x2, WNV is x3, x4 is interaction, intercept is x1
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f*1 + x3*1 + x4*%f*1)))", z, z)
form
Col.WNVcur_Precip_plus1SD.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.WNVcur_Precip_plus1SD.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}

#3.5.3 Do formula steps again for the -1SD, since this changes the sign for the precip and interaction terms, or multiply each by -1
##deltamethod(~ 1/(1+exp(-(x1 - x2*stdvalue + x3*-1 - x4*stdvalue*-1))),  # for +1SD (last term), precip is x2, WNV is x3, x4 is interaction, intercept is x1
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f*1 + x3*-1 + x4*%f*-1)))", z, z)
form
Col.WNVcur_Precip_minus1SD.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.WNVcur_Precip_minus1SD.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}


##3.5.4 Calculating the 1 SE confidence intervals around each Extinction parameter and the 1 standard unit interaction
#EXT Means with 1 se
Col.Precip.mean <- 1/(1+exp(-(Col.Intercept + Col.Precip*stdvalue)))
Col.Precip.meanplus1SE <- Col.Precip.mean + Col.Precip.se
Col.Precip.meanminus1SE <- Col.Precip.mean - Col.Precip.se

Col.WNVcur_Precip.plus1SD <- 1/(1+exp(-(Col.Intercept + Col.Precip*stdvalue + Col.WNVcur*1 + Col.WNVcur_Precip*stdvalue*1)))
Col.WNVcur_Precip.plus1SD.plus1SE    <- Col.WNVcur_Precip.plus1SD + Col.WNVcur_Precip_plus1SD.se
Col.WNVcur_Precip.plus1SD.minus1SE    <- Col.WNVcur_Precip.plus1SD - Col.WNVcur_Precip_plus1SD.se

Col.WNVcur_Precip.minus1SD <- 1/(1+exp(-(Col.Intercept + Col.Precip*stdvalue + Col.WNVcur*-1 + Col.WNVcur_Precip*stdvalue*-1)))
Col.WNVcur_Precip.minus1SD.plus1SE    <- Col.WNVcur_Precip.minus1SD + Col.WNVcur_Precip_minus1SD.se
Col.WNVcur_Precip.minus1SD.minus1SE    <- Col.WNVcur_Precip.minus1SD - Col.WNVcur_Precip_minus1SD.se

#3.5 Putting the colonization pieces into a dataframe
COL_ALL <-  data.frame(stdvalue, Col.Precip.mean, Col.Precip.meanplus1SE, Col.Precip.meanminus1SE, Col.WNVcur_Precip.plus1SD,
                       Col.WNVcur_Precip.plus1SD.plus1SE, Col.WNVcur_Precip.plus1SD.minus1SE, Col.WNVcur_Precip.minus1SD,
                       Col.WNVcur_Precip.minus1SD.plus1SE, Col.WNVcur_Precip.minus1SD.minus1SE)
str(COL_ALL)
#saveRDS(COL_ALL, file = "VIRA_COL_ALL.RDS")




