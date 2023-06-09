## Alex Liu, Marcus Griffiths (2022)
## Data pre-processing for 3D gel imager
## generalized version 6.21.2022

## 1) Set working directory --------------------------------------------------------------------------------------------
getwd()
setwd("C:/Users/USERNAME/FILEPATH") #PC
setwd("C:\\\\Users\\\\USERNAME\\\\FILEPATH") #PC
setwd("/Users/USERNAME/FILEPATH") #macOS

## 2) install & load following packages --------------------------------------------------------------------------------
library(tidyverse)      #includes packages ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats

## 3) Load data files --------------------------------------------------------------------------------------------------
## enter experiment parameter for filenames ----------------------------------------------------------------------------
exp <- "GAA"
species <- "Ta"
expname <- paste0(species,exp)

## enter filenames for datafiles ---------------------------------------------------------------------------------------
Gia3D_csv <- "TaGAA_Gia3D_10-11-2022_gia3d_v2.csv"
DynamicRoots_csv <- "TaGAA_DynamicRoots_2022-10-11_23-31-28.csv"
Biomassdata_csv <- "TaGAA_Biomass_data.csv"
ShootImagedata_csv <- "TaGAA_PlantCVshootdata.csv"

## 4) Load PlantID file for barcode mapping ----------------------------------------------------------------------------
plant_treatment_mapping_csv <- "TaGAA_metadata_12_7_2021.csv"

plantid <- read_csv(plant_treatment_mapping_csv)
plantid <- plantid %>% select(PlantID, Barcode)


## delete?
## 5) Get working directory
# wd <- getwd()
# setwd(paste0(wd,expname))


## 6) Define parse functions -------------------------------------------------------------------------------------------
## Pull out the Gia3D data that we want and make the proper plant ID ---------------------------------------------------
parse_Gia3D <- function(csv){
  data <- read_csv(csv,na = c("NA", "na", "n.a.", ""))
  # get barcode out of 4 columns: Species	Experiment	Plant
  # ex. pennycress	GAA	p0422
  # want TaGAAp0422
  data <- data %>% mutate(Barcode = paste0(species,data$Experiment,data$Plant))
  data <- data %>% select(Barcode
                          ,ImagingDay
                          ,'ConvexHullVolume3D(mm^3)' # RootConvexHullVolume3D_mm3,
                          ,Solidity3D # RootSolidity3D,
                          ,Bushiness3D # RootBushiness3D,
                          ,'Depth3D(mm)' # RootDepth3D_mm,
                          ,'MaximumNetworkWidth3D(mm)' # RootMaximumNetworkWidth3D_mm,
                          ,'LengthDistribution3D(mm)' # RootLengthDistribution3D_mm,
                          ,WidthDepthRatio3D # RootWidthDepth_Ratio
                          )
  return(data)
}

## Parse DynamicRoots data file ----------------------------------------------------------------------------------------
parse_DynamicRoots <- function(csv){
  data <- read_csv(csv) %>%
  # get barcode out of FileNames
  # ex. TaGAAp0101d05_aliu_2021-12-07_12-45-05_rootwork
  # want TaGAAp0101
  separate(FileNames, into = c('Barcode_ImagingDay'), sep = '_') %>% 
  # remove d## from end
  separate(Barcode_ImagingDay, into = c('Barcode','ImagingDay'), sep = -3)
  # data <- data %>% select()
  # data <- data %>% mutate(gsub('.{0,3}$', '', barcode))
  
  ## delete??
  # PlantID	RootTotalCount	RootVolume_mm3	RootLength_mm	RootMeanVolume_mm3	
  # RootMeanLength_mm	RootMeanTortuosity	RootMeanRadius_mm	RootMeanSoilAngle_degrees	RootMeanBranchingAngle_degrees	
  # RootMedianVolume_mm3	RootMedianLength_mm	RootMedianTortuosity	RootMedianRadius_mm	RootMedianSoilAngle_degrees	RootMedianBranchingAngle_degrees	RootCountLateral	RootVolumeLateral_mm3	RootLengthLateral_mm	RootMeanVolumeLateral_mm3	RootMeanLengthLateral_mm	RootMeanTortuosityLateral	RootMeanRadiusLateral_mm	RootMeanSoilAngleLateral_degrees	RootMeanBranchingAngleLateral_degrees	RootMedianVolumeLateral_mm3	RootMedianLengthLateral_mm	RootMedianTortuosityLateral	RootMedianRadiusLateral_mm	RootMedianSoilAngleLateral_degrees	RootMediaBranchingAngle_degrees	RootCountFirstOrderLateral	RootVolumeFirstOrderLateral_mm3	RootLengthFirstOrderLateral_mm	RootMeanVolumeFirstOrderLateral_mm3	RootMeanLengthFirstOrderLateral_mm	RootMeanTortuosityFirstOrderLateral	RootMeanRadiusFirstOrderLateral_mm	RootMeanSoilAngleFirstOrderLateral_degrees	RootMeanBranchingAngleFirstOrderLateral_degrees	RootMedianVolumeFirstOrderLateral_mm	RootMedianLengthFirstOrderLateral_mm	RootMedianTortuosityFirstOrderLateral	RootMedianRadiusFirstOrderLateral_mm	RootMedianSoilAngleFirstOrderLateral_degrees	RootMedianBranchingAngleFirstOrderLateral_degrees	RootDensityFirstOrderLateral_TL	RootDensityFirstOrderLateral_BRTL	InterbranchMeanDistance_mm	InterbranchMedianDistance_mm	RootVolumePrimary_mm3	RootLengthPrimary_mm	RootTortuosityPrimary	RootRadiusPrimary_mm	RootSoilAnglePrimary_degrees
  
  return(data)
}

## Parse biomass data file ---------------------------------------------------------------------------------------------
parse_Biomassdata <- function(csv){
  data <- read_csv(csv)
  return(data)
}

## Parse shoot image data file -----------------------------------------------------------------------------------------
parse_ShootImagedata <- function(csv){
  data <- read_csv(csv)
  return(data)
}

## 7) Load & prepare PlantID file for barcode mapping ------------------------------------------------------------------
plantid <- read_csv(plant_treatment_mapping_csv)

plant <- function(x) {
  len <- 4-nchar(x)
  paste0("p",paste0(integer(len), collapse=""),x)
}
# format of plant should be p00xx
plantid$Plant <- lapply(plantid$Plant,plant )

## Make plantID column  ------------------------------------------------------------------------------------------------
# ex 11122021_TaGAA_132_p0108_tHighN_r1
plantid <- plantid %>% mutate(PlantID = paste(plantid$Date
                                              ,plantid$ExperimentNumber
                                              ,plantid$Geno
                                              ,plantid$Plant
                                              ,paste0("t",plantid$Treatment)
                                              ,paste0("r",plantid$Block)
                                              ,sep = '_'
                                              ))

plantid <- plantid %>% select(PlantID, Barcode)
#metadata is weird, need to add a 4th digit otherwise rsa-Gia will fail
#pattern: regex looking for 'p' followed immediately by a numeric digit
#replace: \\\\1 (group 1), 0 , \\\\2 (group 2)
#plantid$PlantID <- gsub('(p)([0-9])', '\\\\10\\\\2', plantid$PlantID)
#plantid$Barcode <- gsub('(p)([0-9])', '\\\\10\\\\2', plantid$Barcode)

RootGIA3Ddata <- parse_Gia3D(Gia3D_csv)
RootDynamicdata <- parse_DynamicRoots(DynamicRoots_csv)
Biomassdata <- parse_Biomassdata(Biomassdata_csv)
ShootImagedata <- parse_ShootImagedata(ShootImagedata_csv)

## 8) mapping of datafiles with PlantID  -------------------------------------------------------------------------------

##Delete?
# mapping 
#RootGIA3Ddata <- RootGIA3Ddata %>% mutate_at(c('Barcode'), funs(ifelse(. %in% plantid$Barcode, plantid$PlantID[match(., plantid$Barcode)], .)))
#RootDynamicdata <- RootDynamicdata %>% mutate_at(c('Barcode'), funs(ifelse(. %in% plantid$Barcode, plantid$PlantID[match(., plantid$Barcode)], .)))

RootGIA3Ddata <- RootGIA3Ddata %>% mutate(PlantID = ifelse(RootGIA3Ddata$Barcode %in% plantid$Barcode, plantid$PlantID[match(RootGIA3Ddata$Barcode, plantid$Barcode)], RootGIA3Ddata$Barcode))
RootGIA3Ddata$PlantID <- paste(RootGIA3Ddata$PlantID, RootGIA3Ddata$ImagingDay, sep = '_')
RootGIA3Ddata <- RootGIA3Ddata %>% select(PlantID
                                          ,everything()
                                          ,-Barcode
                                          ,-ImagingDay
                                          )

RootDynamicdata <- RootDynamicdata %>% mutate(PlantID = ifelse(RootDynamicdata$Barcode %in% plantid$Barcode, plantid$PlantID[match(RootDynamicdata$Barcode, plantid$Barcode)], RootDynamicdata$Barcode))
RootDynamicdata$PlantID <- paste(RootDynamicdata$PlantID, RootDynamicdata$ImagingDay, sep = '_')
RootDynamicdata <- RootDynamicdata %>% select(PlantID
                                              ,everything()
                                              ,-Barcode
                                              ,-ImagingDay
                                              ,-'Scale(mm)'
                                              ,-Resolution
                                              ,-Threshold
                                              )

ShootImagedata <- ShootImagedata %>% mutate(PlantID = ifelse(ShootImagedata$Barcode %in% plantid$Barcode, plantid$PlantID[match(ShootImagedata$Barcode, plantid$Barcode)], ShootImagedata$Barcode))
ShootImagedata <- ShootImagedata %>% select(PlantID
                                            ,everything()
                                            ,-Barcode
                                            )

Biomassdata <- Biomassdata %>% mutate(PlantID = ifelse(Biomassdata$Barcode %in% plantid$Barcode, plantid$PlantID[match(Biomassdata$Barcode, plantid$Barcode)], Biomassdata$Barcode))
Biomassdata <- Biomassdata %>% select(PlantID
                                      ,everything()
                                      ,-Barcode
                                      )

## 9) Save preprocessed dataframes -------------------------------------------------------------------------------------
dir.create(paste(expname,"_dataprocessing", sep=""), showWarnings = TRUE)

write_csv(RootGIA3Ddata, paste0(expname,"_dataprocessing/",expname,"_GIA3D_data.csv"))
write_csv(RootDynamicdata, paste0(expname,"_dataprocessing/",expname,"_DynamicRoots_data.csv"))
write_csv(Biomassdata, paste0(expname,"_dataprocessing/",expname,"_Biomass_data.csv"))
write_csv(ShootImagedata, paste0(expname,"_dataprocessing/",expname,"_ShootImage_data.csv"))
