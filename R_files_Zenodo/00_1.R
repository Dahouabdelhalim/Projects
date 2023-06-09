#load libraries
library(sp)
library(rredlist)
library(raster)

#load directories of modern and historic ranges as tif files

  #set the directory where modern range data is kept
  modDir <- "C:/Users/JoeG/Documents/IUCN/species_range_split/split_pres_filt_only/raster_100k/reclass/"
  
  #load all tif files from this directory
  modFiles <- list.files(modDir,pattern = ".*tif$")

  
  #set the directory where natural range data is kept
  natDir <- "C:/Users/JoeG/OneDrive - WCMC/WILMA2/work_in_progress/input_data_R/natural_ranges_some_modified/mollweide/species_select/reclass/"
  
  #load all tif files from this directory
  natFiles <- list.files(natDir,pattern = ".*tif$")
  
  #remove species which are missing data
  natFiles <- natFiles[natFiles %in% modFiles]
  modFiles <- modFiles[modFiles %in% natFiles]
  
  #change file names into paths
  modFiles <- paste0(modDir, modFiles)
  natFiles <- paste0(natDir,natFiles)
  
#load land cover and habitat information
  lc <- raster("C:/Users/JoeG/OneDrive - WCMC/WILMA2/raw_data/ESA_CCA_10km_moll.tif")
  cw <- read.csv(("C:/Users/JoeG/OneDrive - WCMC/WILMA2/raw_data/ESA_crosswalk.csv"))
  thisCrs <- crs(lc)
  thisRes <- res(lc)
  thisExt <- extent(lc)

  #load a test range
  modCheck <- raster(modFiles[1])
  natCheck <- raster(natFiles[1])  
  