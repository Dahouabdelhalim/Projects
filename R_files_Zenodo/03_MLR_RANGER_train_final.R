######################################################################
# Woody Cover Estimation using Sentinel-1 time series and LiDAR data
# 
# This script belongs to the Publication: 
#   Urban et. al (2020): Woody Cover Classification in the Savanna Ecosystem 
#         of the Kruger National Park Using Sentinel-1 Time Series. Koedoe, X,X.
#
#
# This script was prepared by Dr. Marcel Urban (marcel.urban@uni-jena.de), Kai Heckel (kai.heckel@uni-jena.de) and
# Patrick Schratz (p.schratz@lmu.de)
# 
# Step 3/4
# Training
#
######################################################################

######################################################################
# LOADING LIBRARIES
######################################################################

require(devtools)
pack_list <- c("splitstackshape","raster","ranger","data.table","parallelMap","mlr", "parallel")
require(pacman)
p_load(pack_list, character.only = TRUE)

######################################################################

######################################################################
# DATA PREPARATION
######################################################################

#import .rda file saved in 01_MLR_tune_spatial.R
homedir = dirname(rstudioapi::getSourceEditorContext()$path)
data1 <- readRDS(paste(homedir,'/S1_A_VH_VV_16_17_lidar.rda',sep=''))

######################################################################

######################################################################
# TRAINING USING RANGER
######################################################################

#definition of formular (fo) which are the layernames
fo <- lidar ~ VH_20160102 + VH_20160114 + VH_20160126 + VH_20160207 + VH_20160302 +
   VH_20160314 + VH_20160326 + VH_20160407 + VH_20160419 + VH_20160501 + VH_20160513 + 
   VH_20160525 + VH_20160606 + VH_20160630 + VH_20160712 + VH_20160724 + VH_20160805 + 
   VH_20160817 + VH_20160829 + VH_20160910 + VH_20160922 + VH_20161004 + VH_20161016 + 
   VH_20161028 + VH_20161109 + VH_20161121 + VH_20161203 + VH_20161215 + VH_20161227 + 
   VH_20170108 + VH_20170120 + VH_20170201 + VH_20170213 + VH_20170225 + VH_20170309 + 
   VH_20170321 + VH_20170402 + VH_20170414 + VH_20170426 + VV_20160102 + VV_20160114 + 
   VV_20160126 + VV_20160207 + VV_20160302 + VV_20160314 + VV_20160326 + VV_20160407 + 
   VV_20160419 + VV_20160501 + VV_20160513 + VV_20160525 + VV_20160606 + VV_20160630 + 
   VV_20160712 + VV_20160724 + VV_20160805 + VV_20160817 + VV_20160910 + VV_20160922 + 
   VV_20161004 + VV_20161016 + VV_20161028 + VV_20161109 + VV_20161121 + VV_20161203 + 
   VV_20161215 + VV_20161227 + VV_20170108 + VV_20170120 + VV_20170201 + VV_20170213 + 
   VV_20170225 + VV_20170309 + VV_20170321 + VV_20170402 + VV_20170414 + VV_20170426 

#training the random forest using ranger
# optionally you can extract the optimal processing parameters from the initial tuning phase to use for training
#rf_tune <- readRDS(paste(homedir,'/rf_tuneGrid_results_mtry_1_to_4_ntrees_10_50_100_300_700_S1_16_17.rda',sep=''))
out <- ranger(data1, formula = fo, mtry = 1, num.trees = 300, num.threads = 3, importance="permutation", verbose = TRUE)
saveRDS(out, paste(homedir,'/rf_train_ranger_S1_16_17.rda',sep=''))

##### TRAINING END

######################################################################
######################################################################
######################################################################





