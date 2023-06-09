

# This script contains the code used in Thornhill et al 2018
# to generate modeled range predictions for each species.




###########################################################################

# STAGE 1: FIT MAXENT MODELS


# set directories, based on parameters file containing these variables
source("user_parameters.r") 
cdir <- filled_climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- maxent_output_dir
splist <- paste0(spdir, 'combined/0_Species_list_v2.rdata')
bgfile <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')
version <- maxent_run_version


# load libraries
options(java.parameters = "-Xmx1g" )
require(dismo)
require(rJava)
require(raster) 
library(doParallel)
library(rgdal)


# climate data setup
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd", "djf", "jja", "ppt"))
mxModelType <- paste(climnames,collapse="-")
predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep="/"))))
names(predictors) = climnames 


# load background data
bg <- readRDS(bgfile)


# load species list
allSpecies <- readRDS(splist)


# define parameter arguments for maxent models
mxArgs <- c("-a", "-z", "outputformat=raw", 
            "maximumbackground=10000", 
            "nothreshold", "nohinge")


# set up cluster for parallel processing
cl <- makeCluster(nodes)
registerDoParallel(cl)


# loop through species in parallel
results <- foreach(i=1:length(allSpecies)) %dopar% {
        
        mySpecies <- allSpecies[i]
        
        # load libraries
        options(java.parameters = "-Xmx1g" )
        require(dismo)
        require(rJava)
        require(raster) 
        
        # load species occurrence data
        pres <- readRDS(paste0(spdir, "/atomic/", mySpecies, ".rdata"))
        coordinates(pres) <- ~longitude + latitude
        orig.project <- "+proj=longlat +ellps=WGS84"
        projection(pres) <- CRS(orig.project)
        occur <- spTransform(pres, projection(predictors))
        
        # directory to write files to
        mx.dir <- paste(odir, version, mySpecies, sep="/")
        if(file.exists(mx.dir) == F) {dir.create(mx.dir, recursive=T)}
        
        # drop occurrences that fall in the water
        valid <- !is.na(extract(predictors[[1]], occur))
        occur <- occur[valid,]
        
        # fit and save maxent model
        mx <- try(maxent(predictors, p=occur, a=bg, path=mx.dir, args=mxArgs))
        saveRDS(mx, file = paste0(mx.dir, "/ModelObject.rdata"))
}

# close computing cluster
stopCluster(cl)






###########################################################################

# STAGE 2: MAKE RANGE PREDICTIONS


# set directories, as above
source("user_parameters.r") 
clim_dir <- filled_climate_data_dir
model_dir <- paste0(project_stem_dir, "/Output/Maxent/v5")
richness_dir <- paste0(project_stem_dir, "/Output/Richness/V5")
pr_dir <- paste0(project_stem_dir, "/Data/Richness/point_richness_rasters")
spdir <- occurrence_data_dir_processed
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')
outdir <- predicted_ranges_dir
spp_dirs <- list.dirs(model_dir, full.names=T, recursive=F)


# load occurrences
allocc <- readRDS(paste0(spdir, "/occurrences_clean.rds"))
species <- unique(allocc$current_name_binomial)
background <- readRDS(maxent_background)


# load climate data
files <- list.files(path=clim_dir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt"))
predictors <- stack(lapply(files, function(x) readRDS(paste(clim_dir, x, sep="/"))))
names(predictors) <- climnames


# california boundary, for use as mask
cali <- rgdal::readOGR("Data/Shapefiles/states", "cb_2014_us_state_500k")
cali <- cali[cali$NAME=="California",]
cali <- spTransform(cali, crs(readRDS(paste0(clim_dir, "/", files[1]))))
cali_mask <- mask(predictors[[1]], cali)


# this function produces a thresholded (or continuous) range map
# from a fitted maxent model and a set of climate rasters
# and optionally also a set of occurrences if distance is to be considered 
binary_range <- function(points=NULL, background=NULL, maxmod=NULL, predictors=NULL, threshold=NULL, 
                         constraint=NULL, constraint_method="none", sigma=50, continuous=FALSE){
        
        # continuous suitability prediction
        pred <- predict(maxmod, predictors)
        if(continuous) return(pred)
        
        # create a raster template matching non-NA climate data
        calrst <- predictors[[1]]
        values(calrst)[!is.na(values(calrst))] <- 0
        
        # if distance-constrained prediction is requested
        if(constraint_method=="distance"){
                
                # create raster of distances to occurences and mask to study area
                dst <- distanceFromPoints(calrst, constraint)
                dst <- mask(dst, calrst)
                
                # transform linear distances to gaussian decay surface
                gauss <- function(x, sigma) exp(-.5 * (x/sigma) ^ 2)
                dst <- calc(dst, function(x) gauss(x, sigma*1000))
                
                # multiply climatic suitability by distance surface
                pred <- pred * dst
        }
        
        # identify and apply a threshold, converting continuous prediction to binary
        extract <- raster::extract
        eval <- evaluate(p=extract(pred, points), a=extract(pred, background))
        thresh <- threshold(eval, stat=threshold)
        pred <- reclassify(pred, c(0, thresh, 0, thresh, 1, 1))
        
        return(pred)
}


# set up cluster for parallel processing
cl <- makeCluster(nodes)
registerDoParallel(cl)


# loop through species in parallel
results <- foreach(spp = spp_dirs,
                   .packages=c("raster", "dismo", "ggplot2", "tidyr", 
                               "dplyr", "grid", "gridExtra", "rgdal")) %dopar% {
                                       
                                       # load species-specific data
                                       points <- allocc[allocc$current_name_binomial==basename(spp),]
                                       maxmod <- readRDS(paste0(spp, "/ModelObject.rdata"))
                                       if(class(maxmod)=="try-error") return("aborted -- no maxent model")
                                       
                                       # parameters for range predictions
                                       threshold <- "spec_sens" # threshold statistic: equal sensitivity-specificity
                                       sigma <- 50 # gaussian distance decay bandwidth, in km
                                       
                                       # standard maxent range prediction
                                       mxm <- binary_range(points, background, maxmod, predictors, threshold)
                                       
                                       # hybrid range prediction, maxent * distance
                                       ptd <- binary_range(points, background, maxmod, predictors, threshold,
                                                           constraint=points, constraint_method="distance")
                                       
                                       # combine and mask
                                       s <- stack(mxm, ptd)
                                       s <- mask(s, cali_mask)
                                       
                                       # save individual layers
                                       writeRaster(s[[1]], paste0(outdir, "/rasters/binary_maxent/", basename(spp), ".tif"), overwrite=T)
                                       writeRaster(s[[2]], paste0(outdir, "/rasters/distance_hybrid/", basename(spp), ".tif"), overwrite=T)
                               }

# close computing cluster
stopCluster(cl)

