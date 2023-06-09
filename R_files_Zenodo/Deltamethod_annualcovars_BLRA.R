#*** This script uses the Delta method to estimate the standard error for Black Rail (BLRA)  colonization and extinction annual covariates graphed in Figure 3. ***
#    It uses the package msm to do the calculations and was written by Steve Beissinger.
#    All occupancy model estimates (parameter means and covariance matrices) were taken from models implemented in Program Presence, so are pasted into this script
# Key references to consult include https://rdrr.io/cran/msm/man/deltamethod.html   as well as Olivier Gimenez's blog cited below and Appendix B in Gentle Intro to Mark manual
# Also Larkin Powell's article https://bioone.org/journals/The-Condor/volume-109/issue-4/0010-5422%282007%29109[949:AVODPU]2.0.CO;2/APPROXIMATING-VARIANCE-OF-DEMOGRAPHIC-PARAMETERS-USING-THE-DELTA-METHOD/10.1650/0010-5422%282007%29109[949:AVODPU]2.0.CO;2.full

# The code is in 4 parts with subsections (e.g., 3.1, 3.2, etc.)
# First section loads some packages, the second explores how to apply the Delta Method using the msm package to do the calculation of the standard error for back-transformed estimates
# The third section applies it to the best BLRA Model with both Main Effects for colonization and all two-way interactions for extinction.
# All covariates were centered and standardized [(x - mean)/SD] and in the case of WNVvi (West Nile Virus vector index) were logged first (ln + 0.1).  
# the fourth section combines the variables into data files for graphing with another script.

#*******************************

#1.  Load some packages
require(msm)  #does delta method 
# library(ggplot2)
library(dplyr)
library(readxl)
library(magrittr)
library(tidyr)



#Set Working Directory
setwd("C:/Users/beis/Box Sync/Black_Rails/Extinction Dynamics/2021/Analysis/Results")  #for Dell desktop

#2. TRY SOME EXAMPLES of the Delta Method
#Info and example provided with the msm package 
example(deltamethod)

#Practice examples of the Delta Method from https://oliviergimenez.github.io/blog/delta-method/
deltamethod(~ 1/(1+exp(-x1)), -0.4473122, 0.3362757^2)  #Olivier's example

#Trying with our notation gives same answer
deltamethod(~ 1/(1+exp(-(x1))), -0.4473122, 0.3362757^2)

#Example from MacKenzie et al. 2018 Occupancy book page 99
covmat = matrix(c(0.30, -0.05, -0.05, 0.40), nrow=2, ncol=2, byrow=TRUE)
x1 <- 1.1
x2 <- 1.8
means <- c(1.1,  1.8)
deltamethod(~ 1/(1+exp(-(x1 + x2))), means, covmat)  #works




#****** 3.0  BEST BLRA MODEL with main effects for Colonization and all three two-way interactions for Extinction. Do extinction first ********

#3.1 Enter the untransformed Parameter values from the Best BLRA Model
InitialPsi <- 0.634254
Col.Intercept <-	-2.476518
Col.Precip <-	0.473476
Col.WNV <-	-0.425366
Col.Freezedays <-	0.187805
Ext.Intercept <-	-1.637213
Ext.Precip <-	-0.251026
Ext.WNV <-	0.029965
Ext.Freezedays <-	0.006246
Ext.WNV_Freezedays <-	0.289374
Ext.WNV_Precip <-	0.332088
Ext.Precip_Freezedays	<- -0.261967

# vector of standard values
stdvalue <- c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
              -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)


#3.2 Some back-transform calculations of Extinction functions with interactions (terms drop out when set to mean of zero). 1 SD refers to one standard unit, as all of the data are standardized
Ext.WNV.mean <- 1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue)))
Ext.WNV_Precip.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + Ext.Precip + Ext.WNV_Precip*stdvalue)))
Ext.WNV_Precip.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + -1*Ext.Precip + Ext.WNV_Precip*stdvalue*-1)))

Ext.Freezedays.mean <- 1/(1+exp(-(Ext.Intercept + Ext.Freezedays*stdvalue)))
Ext.WNV_Freezedays.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + Ext.Freezedays + Ext.WNV_Freezedays*stdvalue)))
Ext.WNV_Freezedays.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + -1*Ext.Freezedays + Ext.WNV_Freezedays*stdvalue*-1)))

Ext.Precip.mean <- 1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue)))
Ext.Precip_Freezedays.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue + Ext.Freezedays + Ext.Precip_Freezedays*stdvalue)))
Ext.Precip_Freezedays.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue - Ext.Freezedays - Ext.Precip_Freezedays*stdvalue)))


#3.3 Bind together the extinction data into a data.frame
Extinction_all <- data.frame(stdvalue, Ext.WNV.mean, Ext.WNV_Precip.plus1SD, Ext.WNV_Precip.minus1SD, Ext.WNV_Freezedays.plus1SD,Ext.WNV_Freezedays.minus1SD,   
                             Ext.Precip.mean, Ext.Freezedays.mean, Ext.Precip_Freezedays.plus1SD, Ext.Precip_Freezedays.minus1SD)
str(Extinction_all)
head(Extinction_all)



# 3.5 Estimating standard errors for Full model using the Delta method. Sections 3.5.1 and 3.5.2 are no longer needed


#3.5.3 Parameter values of the untransformed estimates from BEST MODEL with main effects for Colonization and all two-way interactions for Extinction
InitialPsi <- 0.634254
Col.Intercept <-	-2.476518
Col.Precip <-	0.473476
Col.WNV <-	-0.425366
Col.Freezedays <-	0.187805
Ext.Intercept <-	-1.637213
Ext.Precip <-	-0.251026
Ext.WNV <-	0.029965
Ext.Freezedays <-	0.006246
Ext.WNV_Freezedays <-	0.289374
Ext.WNV_Precip <-	0.332088
Ext.Precip_Freezedays	<- -0.261967

# Some back-transform calculations of extinction functions with interactions (terms drop out when set to mean of zero). 1 SD refers to one standard unit, as all of the data are standardized
Ext.WNV.mean <- 1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue)))
Ext.Precip.mean <- 1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue)))
Ext.Freezedays.mean <- 1/(1+exp(-(Ext.Intercept + Ext.Freezedays*stdvalue)))

Ext.WNV_Precip.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + Ext.Precip + Ext.WNV_Precip*stdvalue)))
Ext.WNV_Precip.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + -1*Ext.Precip + Ext.WNV_Precip*stdvalue*-1)))

Ext.WNV_Freezedays.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + Ext.Freezedays + Ext.WNV_Freezedays*stdvalue)))
Ext.WNV_Freezedays.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.WNV*stdvalue + -1*Ext.Freezedays + Ext.WNV_Freezedays*stdvalue*-1)))

Ext.Precip_Freezedays.plus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue + Ext.Freezedays + Ext.Precip_Freezedays*stdvalue)))
Ext.Precip_Freezedays.minus1SD <-  1/(1+exp(-(Ext.Intercept + Ext.Precip*stdvalue - Ext.Freezedays - Ext.Precip_Freezedays*stdvalue)))



#3.5.4 Variance covariance matrix for the Best Model with rows and columns in order of the covariates above
vcmatrix <- matrix(c(0.038686	,	-0.000733	,	0.000233	,	0.001146	,	0.000475	,	0.000202	,	-0.000621	,	-0.000239	,	-0.001103	,	0.00157	,	0.000777	,	0.000052	,
                       -0.000733	,	0.008828	,	-0.003504	,	0.000216	,	-0.000064	,	0.00234	,	-0.000485	,	-0.000339	,	0.000273	,	-0.000391	,	-0.000312	,	0.000635	,
                       0.000233	,	-0.003504	,	0.008684	,	0.000187	,	0.003994	,	-0.000939	,	0.000311	,	0.000203	,	-0.000498	,	0.000103	,	0.000182	,	-0.000384	,
                       0.001146	,	0.000216	,	0.000187	,	0.009387	,	-0.005894	,	0.000122	,	-0.000162	,	0.001451	,	-0.000614	,	-0.000469	,	-0.000377	,	-0.000196	,
                       0.000475	,	-0.000064	,	0.003994	,	-0.005894	,	0.017257	,	-0.000252	,	0.000306	,	-0.000783	,	0.000681	,	0.000058	,	0.000184	,	0.000055	,
                       0.000202	,	0.00234	,	-0.000939	,	0.000122	,	-0.000252	,	0.0126	,	0.0016	,	-0.005095	,	0.004108	,	-0.007471	,	0.000075	,	0.002819	,
                       -0.000621	,	-0.000485	,	0.000311	,	-0.000162	,	0.000306	,	0.0016	,	0.012597	,	0.00113	,	0.007289	,	-0.002562	,	-0.00535	,	0.002462	,
                       -0.000239	,	-0.000339	,	0.000203	,	0.001451	,	-0.000783	,	-0.005095	,	0.00113	,	0.020607	,	-0.011188	,	0.013525	,	-0.001252	,	0.002512	,
                       -0.001103	,	0.000273	,	-0.000498	,	-0.000614	,	0.000681	,	0.004108	,	0.007289	,	-0.011188	,	0.018764	,	-0.007583	,	-0.006117	,	0.00692	,
                       0.00157	,	-0.000391	,	0.000103	,	-0.000469	,	0.000058	,	-0.007471	,	-0.002562	,	0.013525	,	-0.007583	,	0.022141	,	0.00432	,	0.007116	,
                       0.000777	,	-0.000312	,	0.000182	,	-0.000377	,	0.000184	,	0.000075	,	-0.00535	,	-0.001252	,	-0.006117	,	0.00432	,	0.016726	,	-0.007745	,
                       0.000052	,	0.000635	,	-0.000384	,	-0.000196	,	0.000055	,	0.002819	,	0.002462	,	0.002512	,	0.00692	,	0.007116	,	-0.007745	,	0.018843), nrow=12)
# standard values
stdvalue <- c(-2,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,
              -0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2)

## 3.5.5 Code for Ext WNV*Precip calculation of SE using the Delta method
#3.5.5.1 Step 1 first put together the variance covariance matrix elements in covmat.extwmv_precip and the mean estimates. Did this all the hard way!
covmat.extwnv_precip <- matrix(vcmatrix[6:8,6:8],nrow=3)
covmat.extwnv_precip <- rbind(covmat.extwnv_precip,vcmatrix[11,6:8])
addcol1 <- cbind(vcmatrix[6:8,11])
addcol2 <- rbind(addcol1, vcmatrix[11,11])
covmat.extwnv_precip <- cbind(covmat.extwnv_precip,addcol2)
means.extwnv_precip <- c(Ext.Intercept[1], Ext.Precip [1], Ext.WNV [1], Ext.WNV_Precip[1])

#3.5.5.2 Step 2 write the delta method statement for +1 standard unit (SD) and get for all std values
##deltamethod(~ 1/(1+exp(-(x1 + x2 + x3*stdvalue + x4*stdvalue))), means.extwnv_precip, covmat.extwnv_precip) # for +1SD  precip is x2 and WNV is x3, intercept is x1
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2 + x3*%f + x4*%f)))", z, z)
form
Ext.WNV_Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNV_Precip.se[[i]] <- deltamethod(as.formula(form[i]), means.extwnv_precip, covmat.extwnv_precip)
}

#3.5.5.3 Do formula steps again for the -1SD, since this changes the sign for the precip and interaction terms, or multiply each by -1
##deltamethod(~ 1/(1+exp(-(x1 - x2* + x3*stdvalue - x4*stdvalue))), means.extwnv_precip, covmat.extwnv_precip) # for -1SD (last term), precip is x2, WNV is x3, intercept is x1
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*-1 + x3*%f + x4*%f*-1)))", z, z)
form
Ext.WNV_Precip_minus.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNV_Precip_minus.se[[i]] <- deltamethod(as.formula(form[i]), means.extwnv_precip, covmat.extwnv_precip)
}


#3.5.5.4 for Ext WNV*Freezedays  (here I pasted the portions of variance covariance matrix above output by Program Presence and grabbed the means from above)
covmat.extwnv_freezedays <- matrix(c(0.0126,-0.005095,0.004108,-0.007471, -0.005095,0.020607,	-0.011188,	0.013525, 0.004108,	-0.011188,
         0.018764,	-0.007583, -0.007471,	0.013525,	-0.007583,	0.022141), nrow = 4)
means.extwnv_freezedays <- c(Ext.Intercept[1], Ext.WNV[1], Ext.Freezedays[1], Ext.WNV_Freezedays[1])
#step 2 write the delta method statement
# (~ 1/(1+exp(-(intercept + WNV*stdvalue + Freezedays + Interaction*stdvalue)))
##deltamethod(~ 1/(1+exp(-(x1 + x2*stdvalue + x3 + x4*stdvalue))), means.extwnv_precip, covmat.extwnv_precip) # for +1SD  so here WNV is the x2
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3 + x4*%f)))", z, z)
form
# step 3 run a loop to get it to run across all the stdvalues - this works
Ext.WNV_Freezedays.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNV_Freezedays.se[[i]] <- deltamethod(as.formula(form[i]), means.extwnv_freezedays, covmat.extwnv_freezedays)
}
# Do formula steps again for the -1SD, which changes the delta method formula and results in a different se to subtract from mean. x2 is WNV, x3 is Freezedays
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3*-1 - x4*%f*-1)))", z, z)
form
Ext.WNV_Freezedays_minus.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNV_Freezedays_minus.se[[i]] <- deltamethod(as.formula(form[i]), means.extwnv_freezedays, covmat.extwnv_freezedays)
}
###


# 3.5.5.5 Do this for Ext Precip*Freezedays +1SD  (here I pasted the portions of variance covariance matrix above output by Program Presence and grabbed the means from above)
covmat.extprecip_freezedays <- matrix(c(0.0126,	0.0016,	0.004108,	0.002819, 0.0016,	0.012597,	0.007289,	0.002462, 
                                        0.004108,	0.007289,	0.018764,	0.00692, 0.002819,	0.002462,	0.00692,	0.018843), nrow=4)
means.extprecip_freezedays <- c(Ext.Intercept[1], Ext.Precip[1], Ext.Freezedays[1], Ext.Precip_Freezedays[1])              
#step 2 write the delta method statement
# (~ 1/(1+exp(-(intercept + Precip*stdvalue + Freezedays + Interaction*stdvalue)))
##deltamethod(~ 1/(1+exp(-(x1 + x2*stdvalue + x3 + x4*stdvalue))), means.extprecip_freezedays, covmat.extprecip_freezedays) # for +1SD  so here Precip is the x2
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3 + x4*%f)))", z, z)
form
# step 3 run a loop to get it to run across all the stdvalues - this works
Ext.Precip_Freezedays.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.Precip_Freezedays.se[[i]] <- deltamethod(as.formula(form[i]), means.extprecip_freezedays, covmat.extprecip_freezedays)
}
# Do formula steps again for the -1SD, which changes the delta method formula and results in a different se to subtract from mean. x2 is Precip, x3 is Freezedays
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3*-1 + x4*%f*-1)))", z, z)
form
# step 3 run a loop to get it to run across all the stdvalues - this works
Ext.Precip_Freezedays_minus.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.Precip_Freezedays_minus.se[[i]] <- deltamethod(as.formula(form[i]), means.extprecip_freezedays, covmat.extprecip_freezedays)
}


#3.5.5.6 Delta method estimates for standard errors of each extinction main effect only 
# Precipitation se's
covmat.extprecip = matrix(c(0.0126,	0.0016, 0.0016, 0.0126), nrow=2, ncol=2, byrow=TRUE)
means.extprecip <- c(-1.637213, -0.251026)
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f)))", z)
form
#writing a loop to get it to run across all the stdvalues - this works
Ext.Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.Precip.se[[i]] <- deltamethod(as.formula(form[i]), means.extprecip, covmat.extprecip)
}

# WNV se's
covmat.extwnv = matrix(c(0.0126,	-0.005095, -0.005095, 0.0126), nrow=2, ncol=2, byrow=TRUE)
means.extwnv <- c(-1.637213, 0.029965)
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f)))", z)
form
#writing a loop to get it to run across all the stdvalues - this works
Ext.WNV.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.WNV.se[[i]] <- deltamethod(as.formula(form[i]), means.extwnv, covmat.extwnv)
}

# Freezedays se's
covmat.extfreezedays = matrix(c(0.0126, 0.004108, 0.004108, 0.0126), nrow=2, ncol=2, byrow=TRUE)
means.extfreezedays <- c(-1.637213, 0.006246)
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f)))", z)
form
#writing a loop to get it to run across all the stdvalues - this works
Ext.Freezedays.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Ext.Freezedays.se[[i]] <- deltamethod(as.formula(form[i]), means.extfreezedays, covmat.extfreezedays)
}

##3.5.5.7 Calculating the 1  SE confidence intervals around each Extinction parameter and interaction
#EXT Means with 1 se
Ext.Precip.meanplus1SE <- Ext.Precip.mean + Ext.Precip.se
Ext.Precip.meanminus1SE <- Ext.Precip.mean - Ext.Precip.se
Ext.WNV.meanplus1SE    <- Ext.WNV.mean +  Ext.WNV.se
Ext.WNV.meanminus1SE    <- Ext.WNV.mean -  Ext.WNV.se
Ext.Freezedays.meanplus1SE  <- Ext.Freezedays.mean +  Ext.Freezedays.se
Ext.Freezedays.meanminus1SE  <- Ext.Freezedays.mean -  Ext.Freezedays.se

#Interaction means plus 1 standard unit with 1 se
Ext.WNV_Precip.plus1SD_plusSE   <-    Ext.WNV_Precip.plus1SD +  Ext.WNV_Precip.se
Ext.WNV_Precip.plus1SD_minusSE  <-    Ext.WNV_Precip.plus1SD -  Ext.WNV_Precip_minus.se
Ext.WNV_Precip.minus1SD_plusSE  <-    Ext.WNV_Precip.minus1SD +  Ext.WNV_Precip.se
Ext.WNV_Precip.minus1SD_minusSE  <-    Ext.WNV_Precip.minus1SD -  Ext.WNV_Precip_minus.se

Ext.WNV_Freezedays.plus1SD_plusSE   <- Ext.WNV_Freezedays.plus1SD +  Ext.WNV_Freezedays.se
Ext.WNV_Freezedays.plus1SD_minusSE  <- Ext.WNV_Freezedays.plus1SD -  Ext.WNV_Freezedays_minus.se
Ext.WNV_Freezedays.minus1SD_plusSE  <- Ext.WNV_Freezedays.minus1SD +  Ext.WNV_Freezedays.se
Ext.WNV_Freezedays.minus1SD_minusSE  <- Ext.WNV_Freezedays.minus1SD -  Ext.WNV_Freezedays_minus.se

Ext.Precip_Freezedays.plus1SD_plusSE   <- Ext.Precip_Freezedays.plus1SD +  Ext.Precip_Freezedays.se
Ext.Precip_Freezedays.plus1SD_minusSE   <- Ext.Precip_Freezedays.plus1SD -  Ext.Precip_Freezedays_minus.se
Ext.Precip_Freezedays.minus1SD_plusSE   <- Ext.Precip_Freezedays.minus1SD +  Ext.Precip_Freezedays.se
Ext.Precip_Freezedays.minus1SD_minusSE   <- Ext.Precip_Freezedays.minus1SD - Ext.Precip_Freezedays_minus.se


# 3.5.6 Estimating standard errors for Colonization - using parameters in 3.5.3 using Col.Intercept, Col.Precip, Col.WNV, and Col.Freezedays
## 3.5.6.1   deltamethod(~ 1/(1+exp(-(x1 + x2*0 + x3*stdvalue + x4*0))), means.extwnv_precip, covmat.extwnv_precip) # for +1SD  precip is x2 and WNV is x3, intercept is x1
Covmat.Col <- matrix(c(0.008828,	-0.003504,	0.000216,	-0.000064, -0.003504,	0.008684,	0.000187,	0.003994, 0.000216,	0.000187,	0.009387,	-0.005894, -0.000064,	
                       0.003994,	-0.005894,	0.017257),nrow=4)
Means.Col <- c(Col.Intercept[1], Col.Precip[1], Col.WNV[1], Col.Freezedays[1])  

# 3.5.6.2 Writing a loop to get it to run across all the stdvalues for col precip
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f + x3*0 + x4*0)))", z)
form
Col.Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.Precip.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}

# 3.5.6.3 writing a loop to get it to run across all the stdvalues for col WNV
form <- sprintf("~ 1/(1+exp(-(x1 + x2*0 + x3*%f + x4*0)))", z)
form
Col.WNV.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.WNV.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}

# 3.5.6.4 writing a loop to get it to run across all the stdvalues for col Freezedays
form <- sprintf("~ 1/(1+exp(-(x1 + x2*0 + x3*0 + x4*%f)))", z)
form
Col.Freezedays.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  Col.Freezedays.se[[i]] <- deltamethod(as.formula(form[i]), Means.Col, Covmat.Col)
}

# 3.5.6.5 Test if the previous formulation for Col.Precip.se (3.5.6.2) is correct and it indeed works!
TestCovmat.Col <- matrix(c(0.008828,	-0.003504, -0.003504,	0.008684),nrow=2, ncol=2, byrow=TRUE)
TestMeans.Col <- c(-2.476518,  0.473476)  
z <- stdvalue
form <- sprintf("~ 1/(1+exp(-(x1 + x2*%f)))", z)
form
TestCol.Precip.se <- array(dim=c(41))  #to hold the se's
for(i in 1:41)  {
  TestCol.Precip.se[[i]] <- deltamethod(as.formula(form[i]), TestMeans.Col, TestCovmat.Col)
}

#3.5.6.6 Calculating Colonization means and SEs
#Means
Col.Precip.mean <- 1/(1+exp(-(Col.Intercept + Col.Precip*stdvalue)))
Col.WNV.mean <- 1/(1+exp(-(Col.Intercept + Col.WNV*stdvalue)))
Col.Freezedays.mean <- 1/(1+exp(-(Col.Intercept + Col.Freezedays*stdvalue)))

#Means plus 1 or 2 SE
Col.Precip.plus1se <- Col.Precip.mean + Col.Precip.se
Col.Precip.minus1se <- Col.Precip.mean - Col.Precip.se
Col.WNV.plus1se <- Col.WNV.mean + Col.WNV.se
Col.WNV.minus1se <- Col.WNV.mean - Col.WNV.se
Col.Freezedays.plus1se <- Col.Freezedays.mean + Col.Freezedays.se
Col.Freezedays.minus1se <- Col.Freezedays.mean - Col.Freezedays.se
Col.Precip.plus2se <- Col.Precip.mean + 2*Col.Precip.se
Col.Precip.minus2se <- Col.Precip.mean - 2*Col.Precip.se
Col.WNV.plus2se <- Col.WNV.mean + 2*Col.WNV.se
Col.WNV.minus2se <- Col.WNV.mean - 2*Col.WNV.se
Col.Freezedays.plus2se <- Col.Freezedays.mean + 2*Col.Freezedays.se
Col.Freezedays.minus2se <- Col.Freezedays.mean - 2*Col.Freezedays.se



#******4.    NOW PUT THE PIECES TOGETHER INTO DATA FRAMES FOR GRAPHING. This output is used in graphs for Fig. 3 (scripts in VIRA&BLRA_Annual_Turnover_Fig3.R)

# Extinction pieces starting with the mean values for each interaction +- 1 Standard unit:   Ext.WNV.mean,Ext.WNV_Precip.plus1SD, Ext.WNV_Precip.minus1SD, Extwnv_precip.se
EXT_ALL <-  data.frame(stdvalue, Ext.WNV.mean, Ext.WNV.meanplus1SE, Ext.WNV.meanminus1SE,Ext.WNV_Precip.plus1SD, Ext.WNV_Precip.plus1SD_plusSE, Ext.WNV_Precip.plus1SD_minusSE,
                Ext.WNV_Precip.minus1SD, Ext.WNV_Precip.minus1SD_plusSE, 	Ext.WNV_Precip.plus1SD_minusSE,Ext.WNV_Freezedays.plus1SD, Ext.WNV_Freezedays.plus1SD_plusSE, 	
                Ext.WNV_Freezedays.plus1SD_minusSE, Ext.WNV_Freezedays.minus1SD, Ext.WNV_Freezedays.minus1SD_plusSE,	Ext.WNV_Freezedays.minus1SD_minusSE,
                Ext.Precip.mean, Ext.Precip.meanplus1SE, Ext.Precip.meanminus1SE, Ext.Precip_Freezedays.plus1SD, Ext.Precip_Freezedays.plus1SD_plusSE, 
                Ext.Precip_Freezedays.plus1SD_minusSE, Ext.Precip_Freezedays.minus1SD, Ext.Precip_Freezedays.minus1SD_plusSE,	Ext.Precip_Freezedays.minus1SD_minusSE,
                Ext.Freezedays.mean, Ext.Freezedays.meanplus1SE, Ext.Freezedays.meanminus1SE) 
str(EXT_ALL)
#saveRDS(EXT_ALL, file = "BLRA_EXT_ALL.RDS") # optional save as an RDS file that can be read in using the readRDS("EXT_ALL.RDS") command

# Colonization pieces (Col.Precip.se, Col.WNV.se, Col.Freezedays.se)
COL_ALL <-  data.frame(stdvalue, Col.Precip.mean, Col.WNV.mean, Col.Freezedays.mean, Col.Precip.plus1se, Col.Precip.plus2se, Col.Precip.minus1se, Col.Precip.minus2se,
                       Col.WNV.plus1se, Col.WNV.plus2se, Col.WNV.minus1se, Col.WNV.minus2se, Col.Freezedays.plus1se, Col.Freezedays.plus2se, Col.Freezedays.minus1se,Col.Freezedays.minus2se)
str(COL_ALL)
#saveRDS(COL_ALL, file = "BLRA_COL_ALL.RDS")






