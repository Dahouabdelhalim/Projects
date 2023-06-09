###### This script is for automated cleaning of GBIF and Atlantic Forest Butterflies databases ######

library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)
library(countrycode)
library(raster)
library(CoordinateCleaner)

dat <- read.csv("Infomap_Brassolini_input.csv")

#plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")

ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.5)+
  theme_bw()

# Reference data
getData("countries")->world

# Renaming the column. If you use the latest version of CoordinateCleaner (2.0-6), just rename the column with the ISO-3 values by "iso_a3_eh"
names(world)[names(world) == "ISO"] <- "iso_a3_eh"
  
#convert country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

#flag problems
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude",
                          countries = "countryCode", 
                          species = "species",
                          tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                    "zeros", "seas"),
									country_ref=world)
summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

#Exclude problematic records
dat_cl <- dat[flags$.summary,]
nrow(dat_cl)

#The flagged records
dat_fl <- dat[!flags$.summary,]
nrow(dat_fl)

write.csv(dat_cl, file = "Infomap_Brassolini_input_cleaned.csv")
