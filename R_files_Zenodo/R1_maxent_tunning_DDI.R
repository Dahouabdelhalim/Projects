## R version 4.0.2

T1<-Sys.time()
library(ENMeval)
library(raster)
library(dismo)
library(dplyr)

################################# distribution data #################################
occs.e <- read.csv("H:/SDMs_2020/Cambaroides_japonicus/trial_10th_LULC_OSM/bg_random_Jamie/East_lineage/presence_spThin_5km.csv")
occs.e[,1] <- "C_japonicus_east"

occs.w <- read.csv("H:/SDMs_2020/Cambaroides_japonicus/trial_10th_LULC_OSM/bg_random_Jamie/West_lineage/presence_spThin_5km.csv")
occs.w[,1] <- "C_japonicus_west"

occs.both <- read.csv("H:/SDMs_2020/Cambaroides_japonicus/trial_10th_LULC_OSM/bg_random_Jamie/two_together/presence_spThin_5km.csv")
occs.both[,1] <- "C_japonicus"

################################# east and west polygons #################################
## freshwater regions: OSM
OSM_fw<-raster("H:/SDMs_2020/Cambaroides_japonicus/trial_10th_LULC_OSM/OMS_freshwater/OSM_freshwater.tif")

# east polygon
slope<-raster("H:/SDMs/slope_elevation_30s/chelsa/slope.tif")

japanRange<- extent(138, 146, 40, 46) %>% as("SpatialPolygons") 
crs(japanRange) <- crs(slope)

slope.crop<-crop(slope, japanRange)

east.shp <-matrix(c(
  142.35, 42.3,
  142.5, 42.6,
  142.84, 43.1,
  142.91, 43.2,
  142.89, 43.3,
  143.1, 43.62,
  143.38, 43.75,
  143.51, 43.75,
  143.74, 43.9,
  143.74, 44.2,
  143.74, 46,
  146, 46,
  146, 40,
  142.35, 40,
  142.35, 42.3
), byrow=T, ncol=2)

east.shp <- SpatialPolygons(list(Polygons(list(Polygon(east.shp)), ID=1)))
crs(east.shp) <- crs(slope)

# west polygon
west.shp = rgeos::gDifference(japanRange, east.shp)

################################# environmental layers #################################
files<- list.files(path="H:/SDMs/chelsa_30s/present_30s", pattern='tif',full.names=TRUE)
wdClim <- stack(files)
wdClim.crop<-crop(wdClim, japanRange)

## slope
slope.crop<-crop(slope, japanRange)

## LULC
LULC.files<- list.files(path="H:/SDMs/LULC/RCP45_2010", pattern='tif',full.names=TRUE)
TBLULC <- stack(LULC.files)
names(TBLULC)
names(TBLULC)<-c(
  "barren", "cropland", "forest", "grass", "impervious", 
  "shrub", "snowIce", "urbanGreenSpace", "water", "wetland")

LULC_crop<-crop(TBLULC, japanRange)
LULC_crop <-stack(LULC_crop)

################# only remain freshwater regions ####################
env_stack <-stack(wdClim.crop, LULC_crop, slope.crop, OSM_fw)

## subset
final.var <- c(
  "bio5", "bio6", "bio13", "bio14",
  "slope",
  "forest", "water", "wetland")

env.selected <-subset(env_stack, final.var)
envs<-stack(env.selected)
envs.e = mask(envs, east.shp)
envs.w = mask(envs, west.shp)

############################ randomly generate background data #################################
bg.e <- rasterToPoints(envs.e)[, 1:2] %>% as.data.frame()
names(bg.e) <- c("lon", "lat")

bg.w <- rasterToPoints(envs.w)[, 1:2] %>% as.data.frame()
names(bg.w) <- c("lon", "lat")

bg <- rasterToPoints(envs)[, 1:2] %>% as.data.frame()
names(bg) <- c("lon", "lat")

## tune Maxent ENMeval
eval.e <- ENMevaluate(
  occs=occs.e[,2:3], 
  envs=envs.e, 
  bg=bg.e, 
  partitions = 'block', orientation = 'lon_lon', 
  tune.args =list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5)), 
  algorithm = "maxent.jar", clamp = TRUE)

eval.w <- ENMevaluate(
  occs=occs.w[,2:3], 
  envs=envs.w, 
  bg=bg.w, 
  partitions = 'block', orientation = "lat_lat", 
  tune.args =list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5)), 
  algorithm = "maxent.jar", clamp = TRUE)

eval.both <- ENMevaluate(
  occs=occs.both[,2:3], 
  envs=envs, 
  bg=bg, 
  partitions = 'block', 
  tune.args =list(fc = c("L", "LQ", "H", "LQH"), rm = seq(0.5, 4, 0.5)), 
  algorithm = "maxent.jar", clamp = TRUE)

Sys.time()-T1 ##