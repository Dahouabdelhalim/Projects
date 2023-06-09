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
# Step 4/4
# Prediction
#
######################################################################

######################################################################
# LOADING LIBRARIES
######################################################################

require(devtools)
pack_list <- c("splitstackshape","raster","ranger","data.table","parallelMap","mlr", "parallel", "rgdal")
require(pacman)
p_load(pack_list, character.only = TRUE)

######################################################################

######################################################################
# DATA PREPARATION
######################################################################

# data input
homedir = dirname(rstudioapi::getSourceEditorContext()$path)
mymodel <- readRDS(paste(homedir,'/rf_train_ranger_S1_16_17.rda',sep='')) #RandomForestModel: output from from 03_RANGER_train.R 
data_dir <- '/.../.../'                             		#Folder containing the data for prediction 
out_dir <- '/.../.../'                  						#Folder where to store the prediction

# input of data stack
setwd(homedir)
files <- "S1_A_VH_VV_16_17_subset_example_lower_Sabie"
data_input <- brick(files)

# define nodata values (e.g. -99) and loop through data stack
for(i in 1:nlayers(data_input)){
  print(names(data_input)[i])
  data_input[[i]][data_input[[i]] == -99] <- NA
}

# create data frame
data1 <- as.data.frame(data_input, xy = TRUE)

# change variable names
setnames(data1, c("x", "y", "VH_20160102", "VH_20160114", "VH_20160126",
  "VH_20160207", "VH_20160302", "VH_20160314", "VH_20160326",
  "VH_20160407", "VH_20160419", "VH_20160501", "VH_20160513",
  "VH_20160525", "VH_20160606", "VH_20160630", "VH_20160712",
  "VH_20160724", "VH_20160805", "VH_20160817", "VH_20160829",
  "VH_20160910", "VH_20160922", "VH_20161004", "VH_20161016",
  "VH_20161028", "VH_20161109", "VH_20161121", "VH_20161203",
  "VH_20161215", "VH_20161227", "VH_20170108", "VH_20170120",
  "VH_20170201", "VH_20170213", "VH_20170225", "VH_20170309",
  "VH_20170321", "VH_20170402", "VH_20170414", "VH_20170426",
  "VV_20160102", "VV_20160114", "VV_20160126", "VV_20160207",
  "VV_20160302", "VV_20160314", "VV_20160326", "VV_20160407",
  "VV_20160419", "VV_20160501", "VV_20160513", "VV_20160525",
  "VV_20160606", "VV_20160630", "VV_20160712", "VV_20160724",
  "VV_20160805", "VV_20160817", "VV_20160910", "VV_20160922",
  "VV_20161004", "VV_20161016", "VV_20161028", "VV_20161109",
  "VV_20161121", "VV_20161203", "VV_20161215", "VV_20161227",
  "VV_20170108", "VV_20170120", "VV_20170201", "VV_20170213",
  "VV_20170225", "VV_20170309", "VV_20170321", "VV_20170402",
  "VV_20170414", "VV_20170426"))

# define a variable "coords", which stores the geographic coordinates
coords <- as.data.frame(data1[c("x", "y")])

# remove the x and y coords from the data frame as they should not be used as a
# feature, only for partitioning
data1$x <- NULL
data1$y <- NULL

# (optional) save data frame into file
# saveRDS(data1, "/.../.../S1_A_VH_VV_16_17_subset_example.rda")

# INFO: you can read the saved file with the following command: 
# data1 = readRDS("/.../.../S1_A_VH_VV_16_17_subset_example.rda")

###############################################################################################
# PREDICTION
######################################################################
result <- predict(mymodel, data = data1, verbose=TRUE)

# preparation for output
result_xy = cbind(as.vector(coords$x), as.vector(coords$y), result$predictions)
result_xy = as.data.frame(result_xy)
setnames(result_xy, c("x", "y","wcover"))
coordinates(result_xy) <- ~x+y
gridded(result_xy) <- TRUE

# define projection
projection(result_xy) <- CRS("+proj=utm +zone=36 +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
#projection(result_xy) <- CRS("+proj=longlat +datum=WGS84")

# save woody cover dataset to file
outfile <- raster(result_xy) 
writeRaster(outfile, paste(homedir, '/S1_A_VH_VV_16_17_subset_example_lower_Sabie_woody_cover.tif',sep=''), format="GTiff", datatype='INT2U', overwrite=TRUE, na.rm=TRUE)

##### PREDICTION END 

######################################################################
######################################################################
######################################################################