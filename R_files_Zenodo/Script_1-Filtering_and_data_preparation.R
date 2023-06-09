# load all libraries
library(tidyverse)
library(rgdal)
library(rgeos)
library(sp)
library(MASS)
library(patchwork)
library(raster)
library(ENMeval)
library(geosphere)
library(dismo)
library(sf)
library(plotrix)
library(rgdal)
library(MASS)
library(raster)
library(ecospat)
library(ks)
library(R2jags)
library(runjags)
library(ENMTools)

### INITIAL DATA PREPARATION ###

# read in species information - species determined by the Frogs of Australia application,
# if a species occurs overseas, determined by the IUCN, and threatened species determined from
# Gillespie et al, 2020
species_information <- read.csv('Species_Information.csv')

# read in raw data from GBIF, ALA, FrogID and sensitive data
gbif_raw <- read.table('gbif_raw.csv',
                       sep = '\\t',
                       fill = T,
                       header = T)
ala_raw <- read.csv('ala_raw.csv')
frogID_raw <- read.table('frogID_public_raw.csv',
                         sep = '\\t',
                         fill = T,
                         header = T)
sensitive <- read.csv('Sensitive_data.csv')

# combine data
combined_raw_data <- data.frame(Species = c(gbif_raw$verbatimScientificName,
                                            ala_raw$scientificName,
                                            frogID_raw$scientificName,
                                            sensitive$scientificName),
                                Lat = c(gbif_raw$decimalLatitude,
                                        ala_raw$decimalLatitude,
                                        frogID_raw$decimalLatitude,
                                        sensitive$decimalLatitude),
                                Lon = c(gbif_raw$decimalLongitude,
                                        ala_raw$decimalLongitude,
                                        frogID_raw$decimalLongitude,
                                        sensitive$decimalLongitude),
                                Year = c(gbif_raw$year,
                                         ala_raw$year,
                                         frogID_raw$year,
                                         rep('2018', times = nrow(sensitive))), # year not available for sensitive records, though all from after FrogId launch.
                                Cood_Uncert = c(gbif_raw$coordinateUncertaintyInMeters,
                                                ala_raw$coordinateUncertaintyInMeters,
                                                frogID_raw$coordinateUncertaintyInMeters,
                                                sensitive$coordinateUncertaintyInMeters))
# remove large data files no longer needed
rm(ala_raw, frogID_raw, gbif_raw, sensitive)


### FILTERING ###
# must have all above information
combined_raw_data$Lat <- as.numeric(combined_raw_data$Lat)
combined_raw_data$Lon <- as.numeric(combined_raw_data$Lon)
combined_raw_data$Cood_Uncert <- as.numeric(combined_raw_data$Cood_Uncert)
combined_raw_data <- combined_raw_data[complete.cases(combined_raw_data),]
# coordinate_uncertainty < 1000 
combined_raw_data <- subset(combined_raw_data, Cood_Uncert <= 1000)
# coordinates must have at least 2 decimal places
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
combined_raw_data$DP_Lat <- decimalplaces(combined_raw_data$Lat)
combined_raw_data$DP_Lon <- decimalplaces(combined_raw_data$Lon)
combined_raw_data <- subset(combined_raw_data, DP_Lat >= 2 & DP_Lon >= 2)
# remove duplicates
combined_raw_data <- combined_raw_data[!duplicated(combined_raw_data),]

# Need to get occurrence records for all frogs that occur in the regions where Australian frogs occur overseas
overseas_frogs <- subset(species_information, Overseas_occurrences == 'Yes (IUCN)')$Species
# crop all records to where IUCN designates overseas boundaries, in order to build bias file
# load in the IUCN shapefiles
iucn_shapefile <- readOGR(dsn = 'IUCN_shapefiles',
                          layer = 'data_0')
# convert occurrence data to spatial object for easy cropping
coordinates(combined_raw_data) <- ~Lon+Lat
# complete loop to get occs
occs_for_overseas_bias <- list()
for(i in 1:length(overseas_frogs)) {
  # get the individual frog's shapefile
  frog_shapefile <- iucn_shapefile[iucn_shapefile@data$BINOMIAL==overseas_frogs[i],]
  # buffer by 1 degree - because ranges are broken by ocean here, method used later not applicable
  frog_shapefile_buf <- gBuffer(frog_shapefile, width = 1)
  # crop out occs
  occs <- as.data.frame(crop(combined_raw_data, frog_shapefile_buf))
  # add a column to say which species these occs are for
  occs$Species_for_Bias <- overseas_frogs[i]
  # add to list
  occs_for_overseas_bias[[i]] <- occs
}

# Now, for the remainder, crop to Australia
# first, read in Australian shapefile
world <- readOGR(dsn = 'Shapefiles',
                 layer = 'TM_WORLD_BORDERS-0.3')
aus <- world[world@data$name=='Australia',]
australia_frog_occurrences <- as.data.frame(crop(combined_raw_data, aus))

# clean up environment
rm(combined_raw_data,
   frog_shapefile,
   frog_shapefile_buf,
   occs)

# save image to complete next script
save.image(file = 'Script_1_image.RData')

## FILTERING COMPLETE ##

