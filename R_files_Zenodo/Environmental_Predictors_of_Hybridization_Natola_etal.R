# Environmental Analyses for Natola et al.
# 
# Date: June 2022
#
# Created by: Shawn Billerman and Libby Natola
#
#####################################################################################################

#setwd('/Users/smb223/Dropbox/Manuscript_Drafts/Sapsuckers/Natola_Environmental_Data')
setwd("~/Documents/UBC/Bioinformatics/rnrb/clines/")

#####################################################################################################

## ISOCLINE DATA ##

Oregon <- read.csv("caor_trans_worldclim.csv", header = T)
#	coordinates(Oregon) <- ~Long + Lat

Manning <- read.csv("man_trans_worldclim.csv", header = T)
#	coordinates(Manning) <- ~Long + Lat

North <- read.csv("nbc_trans_worldclim.csv", header = T)
#	coordinates(North) <- ~Long + Lat

require(sp)
require(rgdal)
require(raster)
require(randomForest)
require(rfUtilities)
require(maptools)
require(mapproj)
require(dismo)
require(maps)
require(proj4)
require(geoR)
require(spatialEco)
data(stateMapEnv)

# ITTERATIVE TEST FOR MULTI-COLINEARITY
# Determines which variables show multicolinearity

# For Oregon/California Transect

xdata.Oregon <- Oregon[,c(4:22)]

cl.Oregon <- multi.collinear(xdata.Oregon, p=0.05) 
  for(l in cl.Oregon) {
    cl.test <- xdata.Oregon[,-which(names(xdata.Oregon)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: Temperature_Seasonality"
# [1] "REMOVE VARIABLE: Max_Temperature_of_Warmest_Month"
# [1] "REMOVE VARIABLE: Min_Temperature_of_Coldest_Month"
# [1] "REMOVE VARIABLE: Temperature_Annual_Range"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Driest_Quarter"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Warmest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Wettest_Month"
# [1] "REMOVE VARIABLE: Precipitation_of_Wettest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Driest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Warmest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Coldest_Quarter"

Oregon.trim <- Oregon[,-which(names(Oregon) %in% cl.Oregon)]

# For Manning Transect

Manning <- na.omit(Manning)

xdata.Manning <- Manning[,c(4:22)]

cl.Manning <- multi.collinear(xdata.Manning, p=0.05) 
  for(l in cl.Manning) {
    cl.test <- xdata.Manning[,-which(names(xdata.Manning)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: Temperature_Seasonality"
# [1] "REMOVE VARIABLE: Max_Temperature_of_Warmest_Month"
# [1] "REMOVE VARIABLE: Min_Temperature_of_Coldest_Month"
# [1] "REMOVE VARIABLE: Temperature_Annual_Range"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Wettest_Quarter"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Warmest_Quarter"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Coldest_Quarter"
# [1] "REMOVE VARIABLE: Annual_Precipitation"
# [1] "REMOVE VARIABLE: Precipitation_of_Driest_Month"
# [1] "REMOVE VARIABLE: Precipitation_Seasonality"
# [1] "REMOVE VARIABLE: Precipitation_of_Wettest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Driest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Warmest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Coldest_Quarter"

Manning.trim <- Manning[,-which(names(Manning) %in% cl.Manning)]

# For Northern BC Transect

North <- na.omit(North)

xdata.North <- North[,c(4:22)]

cl.North <- multi.collinear(xdata.North, p=0.05) 
  for(l in cl.North) {
    cl.test <- xdata.North[,-which(names(xdata.North)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: Max_Temperature_of_Warmest_Month"
# [1] "REMOVE VARIABLE: Min_Temperature_of_Coldest_Month"
# [1] "REMOVE VARIABLE: Temperature_Annual_Range"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Warmest_Quarter"
# [1] "REMOVE VARIABLE: Mean_Temperature_of_Coldest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Wettest_Month"
# [1] "REMOVE VARIABLE: Precipitation_of_Wettest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Driest_Quarter"
# [1] "REMOVE VARIABLE: Precipitation_of_Coldest_Quarter"

North.trim <- North[,-which(names(North) %in% cl.North)]


#####################################################################################################

# Random Forests models using regression mode to get a sense of variable importance in each of the three transects
# Using the trimmed response variable set (after removing multicollinear variables) - interesting that each transect had different sets of multicollinear variables, with the northern transect having the least collinearity

# Oregon/California transect

b = 1001 # number of bootstrap replicates

Oregon.rf <- ( rf.model.Oregon <- rf.modelSel(x = Oregon.trim[,c(4:12)], y = Oregon.trim[,14], imp.scale = "mir", ntree = b, final = TRUE))

# Variable importance for selected parameters: 
#                                            imp
# Annual_Mean_Temperature             0.07502136
# Mean_Diurnal_Range                  0.62065420
# Isothermality                       0.14640623
# Mean_Temperature_of_Wettest_Quarter 0.11824121
# Mean_Temperature_of_Coldest_Quarter 0.14011660
# Annual_Precipitation                0.76664621
# Precipitation_of_Driest_Month       0.17618534
# Precipitation_Seasonality           1.00000000
# caor_trans_pts_elev_values          0.14667006
#
#   THRESHOLD    VAREXP          MSE NPARAMETERS
# 1 0.0000000 0.9996672 4.487881e-05           9 
# 2 0.1401166 0.9993175 1.098472e-04           7 
# 3 0.1466701 0.9979542 2.603120e-04           5 
# 4 0.6206542 0.9975664 3.079773e-04 		   3
# 
# While the program chose the model with 9 varibles, the model with only 3 variables (Precipitation Seasonality, Annual Precipitation, and Mean Diurnal Temperature Range) still explained more than 99% of the variation
# I think using the RF model is a good approach for then further exploring these top variables and understanding *how* they impact hybrid index along the isocline

# Manning transect

# Model looking at all climate and elevation variables

Manning.lm <- lmer(caor_loess_envirocline.iso_ancestry ~ caor_trans_pts_elev_values, data = Manning.trim)

b = 1001 # number of bootstrap replicates

Manning.rf <- ( rf.model.Manning <- rf.modelSel(x = Manning.trim[,c(4:9)], y = Manning.trim[,11], imp.scale = "se", ntree = b, final = TRUE))

# Variable importance for selected parameters: 
#                                           imp
# Annual_Mean_Temperature            0.12861007
# Mean_Diurnal_Range                 0.64547276
# Isothermality                      0.01923105
# Mean_Temperature_of_Driest_Quarter 0.32352699
# Precipitation_of_Wettest_Month     1.00000000
# man_trans_pts_elev_values          0.08084320

# Variable importance test for selected parameters:
# 
#    THRESHOLD    VAREXP          MSE NPARAMETERS 
# 3 0.22606853 0.9977912 0.0005825510           3 
# 1 0.00000000 0.9977876 0.0005502580           6 
# 4 0.56498632 0.9974872 0.0006525216           2 
# 2 0.09278491 0.9971812 0.0007181877           4 
#
# The model with 3 variables (Precipitation of Wettest Month, Mean Diurnal Temperature Range, and Mean Temperature of the Driest Quarter) explained more than 99% of the variation

# Northern BC transect

# Model looking at all climate and elevation variables

North.lm <- lmer(caor_loess_envirocline.iso_ancestry ~ caor_trans_pts_elev_values, data = North.trim)

# For the hell of it, trying an RF regression model, to see what happens

b = 1001 # number of bootstrap replicates

North.rf <- ( rf.model.North <- rf.modelSel(x = North.trim[,c(4:14)], y = North.trim[,16], imp.scale = "mir", ntree = b, final = TRUE))
# Variable importance for selected parameters: 
#                                           imp
# Annual_Mean_Temperature             0.05623997
# Mean_Diurnal_Range                  0.24730816
# Isothermality                       0.11348576
# Temperature_Seasonality             1.00000000
# Mean_Temperature_of_Wettest_Quarter 0.05635967
# Mean_Temperature_of_Driest_Quarter  0.44869961
# Annual_Precipitation                0.27619303
# Precipitation_of_Driest_Month       0.10397082
# Precipitation_Seasonality           0.24837400
# Precipitation_of_Warmest_Quarter    0.10603055
# nbc_trans_pts_elev_values           0.04085301
#
# Variable importance test for selected parameters:
# 
#    THRESHOLD    VAREXP          MSE NPARAMETERS 
# 2 0.08016524 0.9994734 3.688646e-05           8 
# 3 0.11348576 0.9994171 4.164173e-05           6 
# 1 0.00000000 0.9993977 4.509674e-05          11 
# 4 0.26228352 0.9990680 6.676678e-05           3 
#
# Similar to the Oregon model, while the program chose the model with 8 varibles, the model with only 3 variables (Temperature Seasonality, Mean Temperature of the Driest Quarter, and Mean Annual Precipitation) still explained more than 99% of the variation
# Interesting that in the northern transect, temperature was much more important compared to both Oregon/California an Manning, in which precipitation was much more important.

# In all cases, all models did well and explained over 99% of variance in ancestry, so for simiplicity, taking the model with only 3 variable for ease of explanation.

#####################################################################################################

# Investigation of specific variables

# In loess_envirocline.iso_ancestry, 0 is "pure" RBSA, and 1 is "pure" RNSA

# Oregon
# top 3 variables - Precipitation Seasonality, Annual Precipitation, and Mean Diurnal Temperature Range

Oregon.lm <- lm(caor_loess_envirocline.iso_ancestry ~ Precipitation_Seasonality + Annual_Precipitation + Mean_Diurnal_Range, data = Oregon.trim)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.28742 -0.09152 -0.01678  0.08705  0.36266 
# 
# Coefficients:
#                             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                4.453e+00  1.483e-01   30.03   <2e-16 ***
# Precipitation_Seasonality  2.128e-02  6.785e-04   31.37   <2e-16 ***
# Annual_Precipitation      -1.068e-03  4.687e-05  -22.79   <2e-16 ***
# Mean_Diurnal_Range        -2.631e-02  8.704e-04  -30.23   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1177 on 509 degrees of freedom
# Multiple R-squared:  0.8877,	Adjusted R-squared:  0.887 
# F-statistic:  1341 on 3 and 509 DF,  p-value: < 2.2e-16
#
# In this model, mean annual precipitation and mean diurnal temperature range had a strong negative relationship with hybrid ancestry; mean annual precipitation makes the most sense, as RNSA would be in areas overall that are drier. Similarly, precipitation seasonality had a strong positive relationship with ancestry, with the environment of RNSA experiencing greater seasonality in precipitation (most precip in winter, versus more evenly distributed along the coast)

# Manning
# top 3 variables - Precipitation of Wettest Month, Mean Diurnal Temperature Range, and Mean Temperature of the Driest Quarter

Manning.lm <- lm(man_loess_envirocline.iso_ancestry ~ Precipitation_of_Wettest_Month + Mean_Diurnal_Range + Mean_Temperature_of_Driest_Quarter, data = Manning.trim)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.36494 -0.19547  0.02806  0.19904  0.33152 
#
# Coefficients:
#                                      Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         2.5669150  0.8055783   3.186 0.001819 ** 
# Precipitation_of_Wettest_Month      0.0011980  0.0011876   1.009 0.315017    
# Mean_Diurnal_Range                 -0.0233071  0.0066310  -3.515 0.000613 ***
# Mean_Temperature_of_Driest_Quarter  0.0028093  0.0005238   5.363 3.82e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2161 on 125 degrees of freedom
# Multiple R-squared:  0.8145,	Adjusted R-squared:  0.8101 
# F-statistic:   183 on 3 and 125 DF,  p-value: < 2.2e-16
#
# Weirdly, precipitation of the wettest month did not come out as signficant here despite being the top ranked predictor variable in the RF model. But, of the two significant variables, mean diurnal temperature range had a negative relationship with hybrid ancestry, while mean temperature of the driest quarter had a positive relationship (hotter summers in interior west where RNSA is)

# Northern BC
# top 3 variables - Temperature Seasonality, Mean Temperature of the Driest Quarter, and Mean Annual Precipitation
North.lm <- lm(nbc_loess_envirocline.iso_ancestry ~ Temperature_Seasonality + Mean_Temperature_of_Driest_Quarter + Annual_Precipitation, data = North.trim)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.26448 -0.09245  0.03573  0.07449  0.29823 
# 
# Coefficients:
                                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                         1.426e+00  1.007e-01   14.17  < 2e-16 ***
# Temperature_Seasonality            -1.573e-04  1.103e-05  -14.26  < 2e-16 ***
# Mean_Temperature_of_Driest_Quarter  3.973e-03  2.096e-04   18.95  < 2e-16 ***
# Annual_Precipitation               -1.726e-04  3.488e-05   -4.95 1.01e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.1098 on 503 degrees of freedom
# Multiple R-squared:  0.8235,	Adjusted R-squared:  0.8225 
# F-statistic: 782.3 on 3 and 503 DF,  p-value: < 2.2e-16
#
# In the Northern BC isocline, all 3 variables are significant, with temperature seasonality and annual precipitation both having a negative relationship with hybrid ancestry, which makes sense: the more highly seasonal precipitation of the Rocky Mountain ecoregion, with greater differences between summer and winter, as well as the overall drier conditions, match with RNSA, while the higher summer temperatures also matches with RNSA 

# I realized a potential problem with using the isoclines themselves. In mapping out the isocline for the California and Oregon birds, I realized that the isocline passes through a very large area of unsuitable habitat, especially for the eastern half of the isocline, so the relationships found in these tests may not be very accurate and may be strongly biased. I think in this case, it may be necessary to try

#####################################################################################################
#####################################################################################################


## INDIVIDUAL LOCALITY DATA ##

# download the relevant bioclim and elevation raster layers (need tiles 12 and 11)
# using raster:: to specify that we want the getData function from the raster package

bioclim_12 <- raster::getData('worldclim', var = 'bio', res = 0.5, lon = -117, lat = 45)
bioclim_11 <- raster::getData('worldclim', var = 'bio', res = 0.5, lon = -123, lat = 42)
alt_12 <- raster::getData('worldclim', var = 'alt', res = 0.5, lon = -117, lat = 45)
alt_11 <- raster::getData('worldclim', var = 'alt', res = 0.5, lon = -123, lat = 42)

# combine the two tiles for each both bioclim and elevation
mergedbio <- merge(bioclim_12, bioclim_11)
mergedalt <- merge(alt_12, alt_11)

# just to be sure, make sure that bioclim and elevation are in the same projection, then add the elevation raster layer to the bioclim raster layer
projection(mergedalt) <- projection(mergedbio)
bioclim <- addLayer(mergedbio, mergedalt)

# break up the large raster layer into separate raster objects
lrs<-numeric(0)
for (i in 1:20){
	lrs[i]<-paste("bc", i, sep = "") 
	assign(lrs[i],raster(bioclim, layer = i))
}

# layer.1 - Annual Mean Temperature
# layer.2 - Mean Diurnal Range (mean of monthly (max temp - min temp))
# layer.3 - Isothermality (layer.2/layer.7)(x100)
# layer.4 - Temperature Seasonality (standard deviation x 100)
# layer.5 - Max Temperature of the Warmest Month
# layer.6 - Min Temperature of the Coldest Month
# layer.7 - Temperature Annual Range (layer.5 - layer.6)
# layer.8 - Mean Temperature of Wettest Quarter
# layer.9 - Mean Temperature of Driest Quarter
# layer.10 - Mean Temperature of Warmest Quarter
# layer.11 - Mean Temperature of Coldest Quarter
# layer.12 - Annual Precipitation
# layer.13 - Precipitation of Wettest Month
# layer.14 - Precipitation of Driest Month
# layer.15 - Precipitation Seasonality (coefficient of variation)
# layer.16 - Precipitation of Wettest Quarter
# layer.17 - Precipitation of Driest Quarter
# layer.18 - Precipitation of Warmest Quarter
# layer.19 - Precipitation of Coldest Quarter
# layer - elevation

# OREGON and CALIFORNIA

Oregon <- read.csv("Locality_Data/CAOR_Locality_Data.csv", header = T)

# define the spatial coordinates
latpr <- Oregon$Lat
lonpr <- Oregon$Long
sites.oregon <- SpatialPoints(cbind(lonpr, latpr))

# using the spatial coordinates, create spatial points data frame
Oregon.spatial <- SpatialPointsDataFrame(sites.oregon, Oregon, match.ID = T)

# extract the climate data from each point
pred.oregon <- raster::extract(bioclim, Oregon.spatial)

Oregon.spatial@data <- data.frame(Oregon.spatial@data, pred.oregon)
Oregon.spatial@data <- na.omit(Oregon.spatial@data)

# ITTERATIVE TEST FOR MULTI-COLINEARITY
# Determines which variables show multicolinearity

# For Oregon/California Transect

xdata.Oregon.spatial <- Oregon.spatial@data[,7:ncol(Oregon.spatial@data)]

cl.Oregon.spatial <- multi.collinear(xdata.Oregon.spatial, p=0.05) 
  for(l in cl.Oregon.spatial) {
    cl.test <- xdata.Oregon.spatial[,-which(names(xdata.Oregon.spatial)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: layer.4"
# [1] "REMOVE VARIABLE: layer.5"
# [1] "REMOVE VARIABLE: layer.6"
# [1] "REMOVE VARIABLE: layer.7"
# [1] "REMOVE VARIABLE: layer.9"
# [1] "REMOVE VARIABLE: layer.10"
# [1] "REMOVE VARIABLE: layer.11"
# [1] "REMOVE VARIABLE: layer.15"
# [1] "REMOVE VARIABLE: layer.16"
# [1] "REMOVE VARIABLE: layer.17"
# [1] "REMOVE VARIABLE: layer.19"

Oregon.spatial@data <- Oregon.spatial@data[,-which(names(Oregon.spatial@data) %in% cl.Oregon.spatial)]

# Random Forests models using regression mode to get a sense of variable importance in each of the three transects
# Using the trimmed response variable set (after removing multicollinear variables) - interesting that each transect had different sets of multicollinear variables, with the northern transect having the least collinearity

# Oregon/California transect

b = 1001 # number of bootstrap replicates

Oregon.spatial.rf <- ( rf.model.Oregon.spatial <- rf.modelSel(x = Oregon.spatial@data[,c(7:15)], y = Oregon.spatial@data[,5], imp.scale = "mir", ntree = b, final = TRUE))

# Variable importance for selected parameters: 
#                 imp
# layer.1  0.30433858
# layer.2  0.14517921
# layer.3  0.00160644
# layer.8  0.66322543
# layer.12 0.04283698
# layer.13 0.04289120
# layer.14 1.00000000
# layer.18 0.90201217
# layer    0.19107382
#
# Variable importance test for selected parameters: 
#
#   THRESHOLD    VAREXP        MSE NPARAMETERS  
# 4 0.6632254 0.7358155 0.03329873           3 
# 3 0.1910738 0.7308428 0.03393683           5 
# 2 0.0428912 0.7220043 0.03469869           7 
# 1 0.0000000 0.7219975 0.03489639           9 
# 
# The top ranked model had 3 variables, which included layer.14 (precipitation of driest month), layer.18 (precipitation of warmest quarter), and layer.8 (mean temperature of wettest quarter)
# This model seems a bit more realistic, with it explaining 73.5% of the variance in ancestry
# I think using the RF model is a good approach for then further exploring these top variables and understanding *how* they impact hybrid index along the isocline

sel.vars.Oregon <- Oregon.spatial.rf$parameters[[4]]
rf.data.Oregon <- data.frame(y = Oregon.spatial@data[,5], Oregon.spatial@data[,sel.vars.Oregon])

rf.final.Oregon <- randomForest(rf.data.Oregon[,2:ncol(rf.data.Oregon)], rf.data.Oregon[,"y"], ntree = b, importance = T, norm.votes = T)

Oregon.rf.ci <- rf.partial.ci(m = rf.final.Oregon, x = rf.data.Oregon[,2:ncol(rf.data.Oregon)])

Oregon.rf.fit <- rf.regression.fit(rf.final.Oregon)

# $fit
#                     [,1]
# RMSE               0.182
# R.squared          0.734
# Cohen.f2           2.758
# Accuracy           0.922
# Overfitting.ratio 36.000
# 
# $message
# [1] "Model is not overfit"

# MANNING PARK

Manning <- read.csv("Locality_Data/MAN_Locality_Data.csv", header = T)

# define the spatial coordinates
latpr <- Manning $Lat
lonpr <- Manning $Long
sites.manning <- SpatialPoints(cbind(lonpr, latpr))

# using the spatial coordinates, create spatial points data frame
Manning.spatial <- SpatialPointsDataFrame(sites.manning, Manning, match.ID = T)

# extract the climate data from each point
pred.manning <- raster::extract(bioclim, Manning.spatial)

Manning.spatial@data <- data.frame(Manning.spatial@data, pred.manning)
Manning.spatial@data <- na.omit(Manning.spatial@data)

# ITTERATIVE TEST FOR MULTI-COLINEARITY
# Determines which variables show multicolinearity

# For Manning Park Transect

xdata.Manning.spatial <- Manning.spatial@data[,8:ncol(Manning.spatial@data)]

cl.Manning.spatial <- multi.collinear(xdata.Manning.spatial, p=0.05) 
  for(l in cl.Manning.spatial) {
    cl.test <- xdata.Manning.spatial[,-which(names(xdata.Manning.spatial)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: layer.3"
# [1] "REMOVE VARIABLE: layer.4"
# [1] "REMOVE VARIABLE: layer.5"
# [1] "REMOVE VARIABLE: layer.6"
# [1] "REMOVE VARIABLE: layer.7"
# [1] "REMOVE VARIABLE: layer.10"
# [1] "REMOVE VARIABLE: layer.11"
# [1] "REMOVE VARIABLE: layer.13"
# [1] "REMOVE VARIABLE: layer.14"
# [1] "REMOVE VARIABLE: layer.16"
# [1] "REMOVE VARIABLE: layer.17"
# [1] "REMOVE VARIABLE: layer.18"
# [1] "REMOVE VARIABLE: layer.19"
# [1] "REMOVE VARIABLE: layer"

Manning.spatial@data <- Manning.spatial@data[,-which(names(Manning.spatial@data) %in% cl.Manning.spatial)]

# Using the trimmed response variable set (after removing multicollinear variables) - interesting that each transect had different sets of multicollinear variables, with the northern transect having the least collinearity

# Manning Park transect

b = 1001 # number of bootstrap replicates

Manning.spatial.rf <- ( rf.model.Manning.spatial <- rf.modelSel(x = Manning.spatial@data[,c(8:13)], y = Manning.spatial@data[,5], imp.scale = "mir", ntree = b, final = TRUE))

# Variable importance for selected parameters: 
#                 imp
# layer.1  0.3761416
# layer.2  1.0000000
# layer.8  0.3651918
# layer.9  0.2899968
# layer.12 0.7226368
# layer.15 0.5354407
#
# Variable importance test for selected parameters: 
#
#   THRESHOLD    VAREXP        MSE NPARAMETERS 
# 2 0.3679293 0.6388847 0.08529924           4 
# 3 0.4557911 0.6296073 0.08628765           3 
# 4 0.6758378 0.6280045 0.08718375           2 
# 1 0.0000000 0.6249151 0.08731896           6 
# 
# The top ranked model had 4 variables, which included layer.2 (mean diurnal range), layer.12 (mean annual precipitation), layer.15 (precipitation seasonality), and layer.1 (mean annual temperature)
# Like the Oregon/California transect model, this model seems a bit more realistic, with it explaining 63.9% of the variance in ancestry

sel.vars.Manning <- Manning.spatial.rf$parameters[[2]]
rf.data.Manning <- data.frame(y = Manning.spatial@data[,5], Manning.spatial@data[,sel.vars.Manning])

rf.final.Manning <- randomForest(rf.data.Manning[,2:ncol(rf.data.Manning)], rf.data.Manning[,"y"], ntree = b, importance = T, norm.votes = T)

# Manning.rf.ci <- rf.partial.ci(m = rf.final.Manning, x = rf.data.Manning[,2:ncol(rf.data.Manning)])

Manning.rf.fit <- rf.regression.fit(rf.final.Manning)

# $fit
#                     [,1]
# RMSE               0.290
# R.squared          0.635
# Cohen.f2           1.739
# Accuracy           0.972
# Overfitting.ratio 11.000
# 
# $message
# [1] "Model is not overfit"


# NORTHERN BRITISH COLUMBIA

North <- read.csv("Locality_Data/NBC_Locality_Data.csv", header = T)

# define the spatial coordinates
latpr <- North$Lat
lonpr <- North$Long
sites.north <- SpatialPoints(cbind(lonpr, latpr))

# using the spatial coordinates, create spatial points data frame
North.spatial <- SpatialPointsDataFrame(sites.north, North, match.ID = T)

# extract the climate data from each point
pred.north <- raster::extract(bioclim, North.spatial)

North.spatial@data <- data.frame(North.spatial@data, pred.north)
North.spatial@data <- na.omit(North.spatial@data)

# ITTERATIVE TEST FOR MULTI-COLINEARITY
# Determines which variables show multicolinearity

# For Northern British Columbia Transect

xdata.North.spatial <- North.spatial@data[,7:ncol(North.spatial@data)]

cl.North.spatial <- multi.collinear(xdata.North.spatial, p=0.05) 
  for(l in cl.North.spatial) {
    cl.test <- xdata.North.spatial[,-which(names(xdata.North.spatial)==l)]
    print(paste("REMOVE VARIABLE", l, sep=": "))
    multi.collinear(cl.test, p=0.05)	
  }
  
# [1] "REMOVE VARIABLE: layer.5"
# [1] "REMOVE VARIABLE: layer.6"
# [1] "REMOVE VARIABLE: layer.7"
# [1] "REMOVE VARIABLE: layer.10"
# [1] "REMOVE VARIABLE: layer.11"
# [1] "REMOVE VARIABLE: layer.13"
# [1] "REMOVE VARIABLE: layer.14"
# [1] "REMOVE VARIABLE: layer.16"
# [1] "REMOVE VARIABLE: layer.17"
# [1] "REMOVE VARIABLE: layer.18"
# [1] "REMOVE VARIABLE: layer.19"
# [1] "REMOVE VARIABLE: layer"

North.spatial@data <- North.spatial@data[,-which(names(North.spatial@data) %in% cl.North.spatial)]

# Using the trimmed response variable set (after removing multicollinear variables) - interesting that each transect had different sets of multicollinear variables, with the northern transect having the least collinearity

# Northern British Columbia transect

b = 1001 # number of bootstrap replicates

North.spatial.rf <- ( rf.model.North.spatial <- rf.modelSel(x = North.spatial@data[,c(7:14)], y = North.spatial@data[,5], imp.scale = "mir", ntree = b, final = TRUE))

# Variable importance for selected parameters: 
#                 imp
# layer.1  0.12622527
# layer.2  0.49424747
# layer.3  0.08219880
# layer.4  1.00000000
# layer.8  0.07938766
# layer.9  0.25562608
# layer.12 0.77655684
# layer.15 0.06653883
#
# Variable importance test for selected parameters: 
#
#    THRESHOLD    VAREXP        MSE NPARAMETERS 
# 4 0.56482481 0.7596907 0.04404313           2 
# 3 0.19092567 0.7253486 0.04997528           4 
# 2 0.08149601 0.7226542 0.05203585           6 
# 1 0.00000000 0.7203192 0.05180149           8 
# 
# The top ranked model had only 2 variables, which included layer.4 (temperature seasonality) and layer.12 (mean annual precipitation). In the next top-ranked model with 4 variables, layer.2 (mean diurnal range) and layer.9 (mean temperature of driest quarter) came out as important as well
# The top-ranked model explained about 76% of variation in ancestry, while the next ranked model explained 72.5% of variation

sel.vars.North <- North.spatial.rf$parameters[[4]]
rf.data.North <- data.frame(y = North.spatial@data[,5], North.spatial@data[,sel.vars.North])

rf.final.North <- randomForest(rf.data.North[,2:ncol(rf.data.North)], rf.data.North[,"y"], ntree = b, importance = T, norm.votes = T)

# North.rf.ci <- rf.partial.ci(m = rf.final.North, x = rf.data.North[,2:ncol(rf.data.North)])

North.rf.fit <- rf.regression.fit(rf.final.North)

# $fit
#                     [,1]
# RMSE               0.210
# R.squared          0.754
# Cohen.f2           3.071
# Accuracy           0.248
# Overfitting.ratio 38.000
#
# $message
# [1] "Model is not overfit"

#####################################################################################################

# Investigation of specific variables from the actual hybrid locality data

# In the files, a q-value of 1 is "pure" RBSA, and 0 is "pure" RNSA

# Oregon
# top 3 variables - layer.14 (precipitation of driest month), layer.18 (precipitation of warmest quarter), and layer.8 (mean temperature of wettest quarter)

Oregon.lm <- lm(q ~ layer.14 + layer.18 + layer.8, data = Oregon.spatial@data)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.92691 -0.09730  0.01233  0.17665  0.54672 
#
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.916910   0.197905   9.686 3.39e-16 ***
# layer.14    -0.051659   0.020457  -2.525   0.0131 *  
# layer.18    -0.006658   0.004627  -1.439   0.1532    
# layer.8      0.004787   0.003224   1.485   0.1406    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2513 on 104 degrees of freedom
# Multiple R-squared:  0.5144,	Adjusted R-squared:  0.5004 
# F-statistic: 36.73 on 3 and 104 DF,  p-value: 2.871e-16
#
# In this model, precipitation of the driest month had a significant negative relationship with hybrid ancestry. Both precipitatin of the warmest quarter and mean temperature of the wettest quarter were not significant in the model, but temperature of the wettest quarter had a positive relationship, suggesting that warmer winter temperatures (if that's when its wettest) correlates with RBSA ancestry

# Manning
# top 4 variables - layer.2 (mean diurnal range), layer.12 (mean annual precipitation), layer.15 (precipitation seasonality), and layer.1 (mean annual temperature)

Manning.lm <- lm(q ~ layer.2 + layer.12 + layer.15 + layer.1, data = Manning.spatial@data)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.46366 -0.09430 -0.02331  0.10055  0.65145 
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)  
# (Intercept)  5.9842918  2.6701123   2.241   0.0309 *
# layer.2     -0.0524288  0.0215479  -2.433   0.0198 *
# layer.12    -0.0010590  0.0008603  -1.231   0.2259  
# layer.15     0.0223903  0.0103443   2.165   0.0368 *
# layer.1      0.0020540  0.0053801   0.382   0.7048  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2846 on 38 degrees of freedom
# Multiple R-squared:  0.6893,	Adjusted R-squared:  0.6566 
# F-statistic: 21.08 on 4 and 38 DF,  p-value: 3.183e-09
#
# In this case, mean diurnal range had a significant negative relationship with ancestry, while precipitiation seasonality had a significant positive relationship with ancestry. The positive precipitation seasonality relationship does not make sense in this case, as I would expect RBSA to be in areas that received more consistent precipitation, at least relative to RNSA

# Northern BC
# top 2 variables - layer.4 (temperature seasonality) and layer.12 (mean annual precipitation)
North.lm <- lm(q ~ layer.4 + layer.12, data = North.spatial@data)
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.36208 -0.19003 -0.07938  0.12565  0.71456 
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  8.079e-01  6.743e-01   1.198 0.234740    
# layer.4     -1.435e-04  7.391e-05  -1.941 0.056089 .  
# layer.12     8.979e-04  2.297e-04   3.909 0.000205 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.2506 on 73 degrees of freedom
# Multiple R-squared:  0.6689,	Adjusted R-squared:  0.6598 
# F-statistic: 73.74 on 2 and 73 DF,  p-value: < 2.2e-16
#
# In the Northern BC transect, only mean annual precipitatin was significant, with a positive effect on ancestry, which makes sense if RBSA occupies wetter areas. Temperature seasonality was very close to signficant, and had a negative effect on ancestry, which also makes sense if RNSA occur in the Rocky Mountain region, which should experience hotter summers and colder winters