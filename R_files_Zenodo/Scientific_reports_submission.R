##
##   Project Name:  Galápagos fish disease monitoring and ecology
##                   
##   Objective:     Combine data from targeted surveys of fish disease prevalence with 
##                  ancillary data to determine the extent of infection and severity 
##                  in the central Galápagos Archipelago.
##
##   Approach:      Prepare elements for the preparation of maps and collate data from 
##                  individual site spreadsheets. Analyses follow progressions for
##                  basic understanding of taxonomic distribution and abundance and 
##                  relationships with geophysical factors.
##
##                  
##
##   Authors:       Robert Lamb and Franz Smith,  Department of Ecology & Evolutionary 
##                  Biology, Brown University with contributions from Anaide Wrublevski 
##                  Aued (Universidade Federal de Santa Catarina), Jenifer Suarez 
##                  (Direccion de Parque Nacional Galápagos), and Jon Witman (Department 
##                  of Ecology & Evolutionary Biology, Brown University)
##
##   Date:          2018-01-14
##
##   Notes:         This file is intended to provide a guide to the different steps necessary to 
##                  conduct the analyses & visual outputs

##
##  1. Set up the core functionality
##
# clear the decks
rm(list=ls())

# call to core packages for data connectivity & manipulation
library(dplyr)
library(tidyr)
library(magrittr)
library(lubridate)
library(forcats)
library(readr)
library(readxl)
library(ncf)
library(ape)
library(car)

# libraries for visualisations
library(ggplot2)
library(GGally)
library(Cairo)
library(extrafont)
library(RColorBrewer)

# required for mapping, spatial analyses and ancillary visuals
library(raster)
library(scales)
library(maptools)
library(gdistance)
library(parallel)
library(cluster)
library(rgdal)
library(ggsn)
library(stringi)

## point to working directory (ensure all data files associated with this manuscript are also in working directory folder)

setwd("~/Desktop/glps_fishDisease")

# call to custom themes and functions

# A function to rename levels of a factor

rename.levels <- function (x, orig, new) 
{
  if (!is.factor(x)) 
    stop("x must be a factor")
  if (length(orig) != length(new)) 
    stop("Number of new labels must equal number of old labels.")
  for (i in 1:length(orig)) {
    levels(x)[levels(x) == orig[i]] <- new[i]
  }
  return(x)
}

## A function to remove plotting attributes

theme_nothing <- function (base_size = 12, base_family = "sans") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(rect              = element_blank(), 
          line              = element_blank(), 
          axis.ticks.length = unit(0.01, "cm"), 
          # axis.ticks.margin = unit(0, "lines"),   ## depreciated
          axis.text         = element_blank(), 
          axis.ticks        = element_blank(), 
          axis.title        = element_blank())
}

## 1) Fig. 1a: Plot ONI over time to show coincidence of 
## extreme departures with fish disease outbreak. Dataset
## is from the Niño 3.4 zone.

# read in the data
oni <- 
  read_excel("ONI_3.4.xlsx") %>%
  tidyr::gather(-Year, key = "Month", value = "ONI") %>%
  mutate(Date = as.Date(paste(Year, Month, 1), format="%Y %b %d")) 


# create plot
oni %>%
  ggplot(aes(x = Date, y = ONI)) +
  geom_ribbon(aes(x = Date, ymin = 0, ymax = ONI * (ONI > 0)),
              fill = alpha('red', 0.7)) +
  geom_ribbon(aes(x = Date, ymax = 0, ymin = ONI * (ONI < 0)),
              fill = alpha('blue', 0.7)) +
  geom_ribbon(aes(x = Date, ymax = 0.5, ymin = -0.5),
              fill = alpha('white', 1)) +
  geom_line() +
  scale_x_date(date_breaks = "5 years", date_label = "%Y") +
  xlab("") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.4)


## 2) Fig. 1b: Creation of Chlorophyll a anomaly plot

## NOTE: skip to second ## NOTE line 218 to load anomaly values from intermediate data folder as a .rda file if access to original MODIS images is not available on local copy

## Set up for mapping and site positions
##
# import the galapagos coastline shapefile
gal   <- shapefile("galcosta.shp")
gal.f <- fortify(gal)



# call to the site locations and other attributes
sites      <- read.csv("siteLocations.csv", header           = TRUE,
                       stringsAsFactors = FALSE,
                       strip.white      = TRUE)

# identify the site positions and transform them to UTM for mapping
sites_xy      <- sites[, c("easting","northing")]
sitePositions <- SpatialPoints(sites_xy, proj4string=CRS("+proj=longlat +ellps=WGS84"))
sitePositions <- spTransform(sitePositions, 
                             CRS = CRS(paste0("+proj=utm +zone=15 +south +datum=WGS84 +units=m",
                                              " +no_defs +ellps=WGS84 +towgs84=0,0,0")))

sites_xy.utm <- coordinates(sitePositions)
sites <- cbind(sites, sites_xy.utm)
sites <- sites[,c(3:ncol(sites))]

# What is the study area?
box.hpts <- chull(x = sites[,5], y = sites[,6])
box.hpts <- c(box.hpts, box.hpts[1])
box.chull.coords <- sites[box.hpts,]
hull <- Polygon(box.chull.coords[5:6], hole = F)
hull@area/1000000

# create the bounding box for the zoom to match the sampling frame
fdzoom <- as.data.frame(bbox(sites_xy.utm)) 

fdzoom$min = fdzoom$min-20000
fdzoom$max = fdzoom$max+20000

fdzoom = as.matrix(fdzoom)



# call to chlorophyll data
load("predictors_chla.rda")

##
##  Create anomaly calculations
##
# set list of months

anomaly_date <- c("2016.01.05", "2016.01.21", "2016.01.29", "2015.12.29")

# get month from anomaly date
months_for_anomaly <- anomaly_date %>% stri_sub(from = 6, to = 7)

# identify positions for month
months_to_select <- which(names(predictors_chla) %>% stri_sub(from = 7, to = 8) %in% months_for_anomaly)

# subset individual month
anomaly_stack <- predictors_chla[[months_to_select]]

# calculate mean
anomaly_month_raster <- anomaly_stack %>% mean(na.rm = TRUE)

# transform date to select
date_to_select <- paste0("X", anomaly_date) %>% gsub("-", ".", .)
# date_to_select <- anomaly_date %>% gsub("-", ".", .)

# get position of month for calculation
date_position <- which(names(predictors_chla) %in% date_to_select)
# get average across multiple dates for anomaly

date_raster = predictors_chla[[date_position]]

anomaly_date_raster = date_raster %>% mean(na.rm = TRUE)

# calculate anomaly
anomaly_raster <- (anomaly_date_raster - anomaly_month_raster) 

## -- project to utm -- ##
# set extent from source
anomaly_raster_extent <-
  anomaly_raster %>% 
  projectExtent(crs = CRS(paste0("+proj=utm +zone=15 +south +datum=WGS84",
                                 " +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")))

# modify resolution to 'square' to original
res(anomaly_raster_extent) <- anomaly_raster_extent %>% res() %>% min() -100

# reproject
anomaly_raster %<>%
  projectRaster(from = anomaly_raster, to = anomaly_raster_extent)

# NOTE: load raster values directly if no access to original MODIS images on local copy

load("chla_anomaly_raster.rda")  

# create a function to overlay the sites
addCoastline <- function() {
  plot(galcosta, col = "grey70", border = "black", add = TRUE)
}

# set colour palette
# c_palette <- colorRampPalette(c("purple","blue","cyan","green","yellow","red"))(255)
c_palette <- viridis_pal()(255)

# plot figure to screen

chla_anomaly <- anomaly_raster %>%
  plot(col    = c_palette,
       addfun = addCoastline,
       #main   = paste0("modis chla anomaly ", anomaly_date),
       main = "",
       zlim = c(-1, .2))


# Return anomaly values for each site     
chla <- data.frame(coordinates(sites[,5:6]),
                   sites$siteName, 
                   extract(anomaly_raster, sites[,5:6])) %>%
  mutate(Type = "chla")
names(chla) <- c("easting", "northing", "siteName", "value", "type")

## 3) Fig.1c: Creation of Sea Surface Temperature anomaly plot


## NOTE: skip to second ## NOTE line 296 to load anomaly values from intermediate data folder as a .rda file if access to original MODIS images is not available on local copy

# call to temperature data
load("predictors_sst.rda")

### identify positions for month
months_to_select <-  which(names(predictors_sst) %>% stri_sub(from = 7, to = 8) %in% months_for_anomaly)

# subset individual month
anomaly_stack <- predictors_sst[[months_to_select]]

# calculate mean
anomaly_month_raster <- anomaly_stack %>% mean(na.rm = TRUE)

# transform date to select
date_to_select <- paste0("X", anomaly_date) %>% gsub("-", ".", .)

# get position of month for calculation
date_position <- which(names(predictors_sst) %in% date_to_select)

# get average across multiple dates for anomaly

date_raster = predictors_sst[[date_position]]

anomaly_date_raster = date_raster %>% mean(na.rm = TRUE)

# calculate anomaly
anomaly_raster <- anomaly_date_raster - anomaly_month_raster


## -- project to utm -- ##
# set extent from source
anomaly_raster_extent <-
  anomaly_raster %>% 
  projectExtent(crs = CRS(paste0("+proj=utm +zone=15 +south +datum=WGS84",
                                 " +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")))

# modify resolution to 'square' to original
res(anomaly_raster_extent) <- anomaly_raster_extent %>% res() %>% min() -100

# reproject
anomaly_raster %<>%
  projectRaster(from = anomaly_raster, to = anomaly_raster_extent)

# NOTE: load raster values directly if no access to original MODIS images on local copy

load("sst_anomaly_raster.rda")  

# create a function to overlay the sites
addCoastline <- function() {
  plot(galcosta, col = "grey70", border = "black", add = TRUE)
}

# set colour palette
# c_palette <- colorRampPalette(c("purple","blue","cyan","green","yellow","red"))(255)
c_palette <- viridis_pal(option = "magma")(255)

# plot figure to screen

anomaly_raster %>%
  plot(col    = c_palette,
       addfun = addCoastline,
       main   = "")

# Return anomaly values for each site     
sst <- data.frame(coordinates(sites[,5:6]),
                  sites$siteName, 
                  extract(anomaly_raster, sites[,5:6])) %>%
  mutate(Type = "sst")
names(sst) <- c("easting", "northing", "siteName", "value", "type") 

## 4) Fig. 3a: Creation of map displaying study sites and number of species affected by ulcerative skin disease


# call to the fish IDs
speciesCodes <- read.csv("fishID.csv", header = TRUE)

# call to the historic fish transect data
fishComm <- read_excel("community_transects.xlsx") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  tidyr::gather(sizeClass, Value, -Site, -Date, -Observer, -Depth, -Temperature, -Water.movement,
                -Site.Exposure, -Substrate, -Transect, -Block, -Transect.Length,
                -Block.Length, -Species, -Sick) %>%
  group_by(Site, Date, Depth, Species, Sick) %>%
  summarise(Abundance = sum(Value, na.rm = TRUE)) %>%
  left_join(speciesCodes) %>%                      
  data.frame()

# Bring in community-wide sick fish transects

focus_community = read_excel("Sick-Fish-Full-Transects.xlsx") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  tidyr::gather(sizeClass, Value, -Site, -Date, -Observer, -Depth, -Temperature, -Water.movement,
                -Site.Exposure, -Substrate, -Transect, -Block, -Transect.Length,
                -Block.Length, -Species, -Sick) %>%
  group_by(Site, Date, Depth, Species, Sick) %>%
  summarise(Abundance = sum(Value, na.rm = TRUE)) %>%
  left_join(speciesCodes) %>%                      
  data.frame()

# Bring in focus species transects

focus_species = read_csv("focus_species_transects.csv") %>%
  rename(Date = date, Site = siteName, Depth = depth, Block = block, Block.Length = block.length, Transect.Length = transect.length, Observer = diver, Scientific.Name = genusSpecies) %>%
  mutate(Sick = as.character(as.factor(infectionCoverage))) %>%
  mutate(Sick = fct_recode(Sick, Y = "0.05", Y ="0.15", Y ="0.1", Y ="0.4", Y ="0.5", Y ="0.2", Y ="0.25", Y ="0.3", Y ="0.8", Y ="0.75", Y ="1", Y ="0.6", Y ="0.7", N = "0")) %>%
  dplyr::select(-c(infectionCoverage, commonName, transect.width)) %>%
  tidyr::gather(sizeClass, Value, -Site, -Date, -Observer, -Depth, -Block, -Transect.Length,
                -Block.Length, -Sick, -Scientific.Name) %>%
  replace_na(replace = list(Value = 0)) %>%
  mutate(Value = as.integer(Value)) %>%
  group_by(Site, Date, Depth, Scientific.Name, Sick) %>%
  summarise(Abundance = sum(Value, na.rm = TRUE)) %>%
  left_join(speciesCodes) %>%                      
  data.frame()

# Merge

focus = rbind(focus_community, focus_species)

# bring in roving counts
roving   <- read_excel("roving_transects.xlsx") %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  tidyr::gather(sizeClass, Value, -Site, -Date, -Depth, -Patch, -Observer,
                -Scientific.Name, -commonName, -infectionCoverage, -Sick, -Notes) %>%
  group_by(Site, Date, Depth, Scientific.Name, Sick) %>%
  summarise(Abundance = sum(Value, na.rm=TRUE)) %>%                           
  data.frame()

# Select columns of interest from transects and re-order
trans  <- fishComm[, c(1:3, 5, 6, 8)]
trans1 <- trans[,c(1:3, 6, 4, 5)]

trans2  <- focus[, c(1:3, 5, 6, 8)]
trans3 <- trans2[,c(1:3, 6, 4, 5)]

# Merge datasets

trans4 = rbind(trans1, trans3)
allfish <- rbind(trans4, roving) %>%
  mutate(Sick = as.factor(Sick)) %>%
  mutate(Sick = fct_recode(Sick, Y = "0.05", Y = "0.1", N = "Isopod.2", N = "undefined")) %>%
  filter(Scientific.Name != "NA") %>%
  mutate(Scientific.Name = fct_recode(Scientific.Name, `Anisotremus interruptus` = "Anistremus interruptus", `Myripristis leiognathus` = "Myripristis leigognathos", `Paranthias colonus` = "Paranthius colonus", `Pareques perissa` = "Pareques pareques", `Prionurus laticlavius` = "Prionurus laticlavus", `Microspathodon dorsalis` = "Microspathadon dorsalis", `Johnrandallia nigrirostris` = "johnrandallia nigrirostris"))


# Site levels we do not want
no_site = c("Baltra Protected", "Camano_Protected", "Camano.Protected", "Punta Nunez", "Cerro Gallina", "palmas.protected", "crater")

# Get diversity and disease data
disease.dist <- allfish %>%
  mutate(Site = fct_recode(Site, norte.baltra = "baltra", norte.baltra = "Baltra", rocas.beagle = "beagle", rocas.beagle = "Beagle", cousins = "Cousins", guy.fawkes = "Guy Fawkes", islote.gardner = "Gardner", rocas.gordon = "gordon", pinzon = "Pinzon", camano.exposed = "Camano.Exposed", ch.pequeno = "4 Hermanos", daphne.menor = "Daphne", daphne.menor = "daphne", la.botella = "La Botella", las.cuevas = "Las Cuevas", palmas.exposed = "Palmas Exposed", palmas.protected = "Palmas Protected", islote.tortuga = "tortuga", champion = "Champion", ch.gem.west = "cuarto.hermano", ch.gem.east = "tercero.hermano")) %>%
  rename(siteName = Site) %>%
  mutate(Scientific.Name = as.factor(Scientific.Name)) %>%
  group_by(siteName) %>%
  summarise(Richness = length(unique(Scientific.Name)), 
            Healthy  = length(unique(Scientific.Name[which(Sick!="Y")])), 
            Sick     = length(unique(Scientific.Name[which(Sick=="Y")]))) %>%
  right_join(sites[,c("siteName", "easting", "northing")])
# mutate(site.abb = rename.levels(as.factor(siteName), orig = c("cousins", "daphne.menor", "norte.baltra", "rocas.beagle", "guy.fawkes", "rocas.gordon", "camano.exposed", "ch.gem.west", "ch.pequeno", "ch.gem.east", "islote.tortuga", "corona.diablo", "champion", "las.cuevas"), new = c("RC", "DM", "BAL", "RB", "GF", "GOR", "CM", "CGW", "CHE", "CGE", "IT", "CD", "CHA", "LC"))) %>%

# Plot
disease_map <- disease.dist %>%
  filter(!siteName %in% no_site) %>%
  ggplot(aes(easting, northing)) +
  geom_polygon(data=gal.f, aes(x=long, y=lat, group=group), fill="grey75", colour="grey75", size=0.2) +
  geom_point(aes(colour=Sick), alpha = .7, size = 8) +
  geom_point(colour="black", alpha = .7, size = 8, shape = 1) +
  theme_nothing() +
  coord_equal(ratio = 1, xlim = c(fdzoom[1], fdzoom[3]), ylim = c(fdzoom[2], fdzoom[4])) +
  scale_colour_gradientn(colours = rev(viridis_pal()(11))) +
  scalebar(gal.f, dist = 15, location = "bottomleft", anchor = c(x=fdzoom[1], y = fdzoom[2]+17000), height = .01, st.dist = .01, st.bottom = F) +
  north(gal.f, location = "bottomleft", anchor = c(x=fdzoom[1]+650, y = fdzoom[2]+31000), symbol = 3, scale = .04)

## How many fish species were observed, and how many of these had signs of infection?

total_richness <- length(unique(allfish$Scientific.Name)) - 3 #for Chelonia mydas (green sea turtle) and Zalophus californianus (sea lion), Arctocephalus galapagoensis (fur seal))

## 5) Fig. 3b: Creation of barplot showing mean and standard error of abundance for focal target species

# create a list of sites to exclude
sites_to_exclude <-
  c("palmas.protected", "islote.gardner", "la.botella", "pinzon", "crater")  

# Bring in file with total number of species infected for plotting purposes

load("community_prevalence.rda")

# call to the survey data
prevalence <- read_csv("focus_species_transects.csv") %>%
  mutate(date = as.POSIXct(date, tz = "UTC")) %>%
  filter(siteName != "ch.pequeno")

# Call to ch.pequeno on Jan 13

ch.pequeno <- 
  read_csv("focus_species_transects.csv") %>%
  mutate(date = as.POSIXct(date, tz = "UTC")) %>%
  filter(siteName == "ch.pequeno",
         date == "2016-01-12 19:00:00 EST")

# data from november 2016
nov_16_prevalence <- 
  read_csv("Hermanos_nov_2016_focus_species_transects.csv") %>%
  mutate(date = as.POSIXct(date, tz = "UTC")) 

# data from july 2017
jul_17_prevalence <- 
  read_excel("Hermanos_july_2017_focus_species_transects.xlsx")

# Combine and curate
prevalence %<>%
  bind_rows(nov_16_prevalence, jul_17_prevalence, ch.pequeno) %>%
  dplyr::filter(genusSpecies != "") %>%
  tidyr::gather(sizeClass, abundance, -siteName, -date, -depth, 
                -block, -block.length, 
                -transect.length, -transect.width,
                -diver, -genusSpecies, -commonName, 
                -infectionCoverage) %>%
  mutate(sizeClass = gsub("X", "", sizeClass)) %>%
  mutate(genusSpecies = gsub ("Myripristis leigognathos", "Myripristis leiognathus", genusSpecies)) %>%
  mutate(genusSpecies = gsub ("Microspathadon dorsalis", "Microspathodon dorsalis", genusSpecies)) %>%
  mutate(abundance = as.numeric(as.character(abundance))) %>%
  mutate(siteName = as.factor(siteName)) %>%
  replace_na(list(abundance = 0)) %>%
  right_join(sites[,c("siteName", "easting", "northing")])               


# modify the siteNames to match the current convention
prevalence$siteName <- rename.levels(factor(prevalence$siteName), 
                                     orig = c("camano.exposed", "cousins",
                                              "hermanos.west", "hermanos.east",
                                              "daphne.menor", "gordon",
                                              "guy.fawkes", "islote.tortuga",
                                              "norte.baltra", "rocas.beagle",
                                              "hermanos.south", "hermanos.2",
                                              "las.cuevas", "champion",
                                              "corona.diablo", "palmas.exposed",
                                              "palmas.protected"),
                                     new  = c("camano.exposed", "cousins",
                                              "hermanos.gemelo.oeste", "hermanos.pequeno",
                                              "daphne.menor", "rocas.gordon",
                                              "guy.fawkes", "islote.tortuga",
                                              "norte.baltra", "rocas.beagle",
                                              "hermanos.grande", "hermanos.gemelo.este",
                                              "las.cuevas", "champion",
                                              "corona.diablo", "palmas.exposed",
                                              "palmas.protected"))

# convert the date
prevalence$sampleDate <- 
  as.Date(as.character(prevalence$date, format="%Y-%m-%d %H:%M:%S"))

##
##  Modify the data structrue & summarise
##
# create a copy of the data
dat <- prevalence

# put the different depths into  strata
sort(unique(dat$depth))
dat$strata <- 
  rename.levels(as.factor(dat$depth), 
                orig = c("5", "6", "7", "10", "11", "12", "13", "14", "15", "16", NA),
                new  = c("shallow", "shallow", "shallow", "deep", "deep", "deep", "deep", "deep", "deep", "deep", NA))

# List key taxa
keyTaxa <- 
  c("Holacanthus passer", 
    "Myripristis leiognathus",
    "Stegastes beebei",
    "Microspathodon dorsalis")

# summarise the prevalence data by depth strata for total and calculate overall proportion infected during peak outbreak
prev <- 
  dat %>%
  mutate(period = substr(date, start = 1, stop = 7)) %>%
  filter(period == "2016-01") %>%
  filter(genusSpecies %in% keyTaxa) %>% 
  group_by(siteName, genusSpecies, strata, block, date) %>%
  summarise(total      = sum(abundance, na.rm = TRUE),
            infected   = sum(abundance[infectionCoverage > 0], na.rm = TRUE)) %>%
  mutate(prevalence = infected / total,
         healthy    = total - infected) %>%
  #group_by(siteName, genusSpecies) %>%
  tidyr::gather(healthy, infected, key = "disease", value = "abundance") %>%
  #mutate(abundance = log10(abundance+1)) %>%
  ungroup() %>%
  mutate(genusSpecies = gsub(" ", "\\n", genusSpecies)) 


##  Summarize peak prevalence at 4 Hermanos east (pequeno)

peak_prev <- prev %>%
  filter(siteName == "ch.pequeno") %>%
  group_by(genusSpecies) %>%
  summarise(total = sum(total),
            prevalence = mean(prevalence, na.rm = T))

## Arrange sites by distance to focal disease center at 4 Hermanos East (pequeno)

disease.dists <- as.matrix(dist(cbind(sites$easting, sites$northing)))

epi_sites <- sites %>%
  mutate(dist_4hermp = disease.dists[,11]) %>%
  dplyr::select(siteName, dist_4hermp) %>%
  arrange(dist_4hermp, siteName)

##
## summarize data for plotting
histo_2 <- prev %>%
  filter(!siteName %in% sites_to_exclude) %>%
  group_by(siteName, genusSpecies, disease) %>%
  dplyr::summarize(mabund      = mean(abundance), 
                   seabund     = sd(abundance) / sqrt(length(abundance)), 
                   mprevalence = mean(prevalence, na.rm = T),
                   n = length(prevalence)) %>%
  right_join(sites, by = "siteName") %>%
  ungroup() %>%
  right_join(epi_sites, by = "siteName") %>%
  arrange(dist_4hermp, siteName) %>%
  # arrange(mprevalence, siteName) %>%
  mutate(siteName = factor(siteName, unique(siteName))) %>%
  mutate(site.abb = rename.levels(siteName, orig = c("cousins", "daphne.menor", "norte.baltra", "rocas.beagle", "guy.fawkes", "rocas.gordon", "palmas.exposed", "camano.exposed", "ch.gem.west", "ch.pequeno", "ch.gem.east", "ch.grande", "islote.tortuga", "corona.diablo", "champion", "las.cuevas"), new = c("RC", "DM", "BAL", "RB", "GF", "GOR", "PAL", "CM", "CGW", "CHE", "CGE", "CHG", "IT", "CD", "CHA", "LC"))) %>%
  mutate(upper = mabund + seabund, 
         lower = mabund - seabund) %>%
  filter(genusSpecies != "NA") 

# plot in order of distance from disease epicenter


disease_histograms_2 <- histo_2 %>%
  ggplot(aes(site.abb, mabund, fill = factor(disease), ymax = upper, ymin = lower)) +
  geom_bar(stat = "identity", position = "dodge", width = 1) +
  geom_linerange(stat = "identity", position = position_dodge(1)) +
  facet_grid(genusSpecies~., scales = "free_y") +
  theme_bw() +
  scale_fill_manual(values = c("grey", "grey50")) +
  scale_color_manual(values = c("grey", "grey50")) +
  # scale_fill_gradient(low = "blue", high = "yellow") +
  #scale_y_continuous(labels = c("2", "1", "0", "1", "2")) +
  #geom_label(data = comm.prev, aes(label = Sick, x = 1.2, y = 1.75, family = "Helvetica", fontface = "bold", size = 8), color = "black", fill = "white", label.r = unit(0.02, "lines"), label.size = 0) +
  #scale_color_gradientn(colours = rev(viridis_pal()(11))) +
  labs(x = "", y = expression("Abundance • 50m"^"-2")) +
  theme(axis.text.x = element_text(colour="black", size=13, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 11),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold")) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(),
        #strip.text = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

## 6) Fig. 4) Spatial auto-correlation analysis of prevalence for M. leiognathus and S. beebei

#Read in data

load("geographic_prevalence.rda")
load("species_prevalence.rda")

# Wrangle species-specific prevalence data

spec_prev <- prev %>%
  group_by(siteName, genusSpecies) %>%
  summarise(mean.prev = mean(prevalence, na.rm = T)) %>%
  filter(mean.prev != 'NaN') %>%
  spread(genusSpecies, mean.prev)

disease.dist <- merge(disease.dist, spec_prev, by = "siteName")

# Estimate number of distance classes appropriate using Sturge's rule 
1+3.3*(log(length(unique(sites$siteName))))

#Reformat data
spec_prev = disease.dist %>%
  tidyr::gather(-siteName, -Richness, -Healthy, -easting, -northing, key = Species, value = Prevalence) %>%
  mutate(Species = fct_recode(Species, Community = "Sick"))

spec_prev$Prevalence[which(is.na(spec_prev$Prevalence))] =0

# Create vectors for storing correlation data

species <- unique(spec_prev$Species)

storage = vector("list", length(species))
list(1:length(species))

for(i in 1:length(species)) {
  data = spec_prev[which(spec_prev$Species == species[i]),] %>%
    filter(Prevalence!="NA")
  spec_fit <- correlog(x=data$easting, y=data$northing, z=data$Prevalence, increment=12000, resamp=1000, na.rm = T)
  corr.plot <- as_tibble(cbind(as.numeric(spec_fit$mean.of.class), spec_fit$correlation, spec_fit$p, spec_fit$n)) %>%
    rename(Distance = V1, Correlation = V2, p = V3) %>%
    mutate(Color = if_else(p < .05, "black", "white")) %>%
    mutate(km = Distance/1000,
           Species = species[i])
  storage[[i]] = corr.plot
}

all.corr <- rbind(storage[[4]], storage[[5]])

# Visualize

corr.plot <- all.corr %>%
  ggplot(aes(km, Correlation)) +
  #geom_point(aes(color = V4, shape = Color, size = 3)) +
  geom_point(aes(color = Color), size = 3) +
  geom_point(shape = 21, size = 3) +
  geom_line(aes(linetype = as.factor(Species)), lwd = .35) +
  #scale_color_viridis()+
  scale_color_grey() +
  theme_bw() +
  #facet_wrap(~Species) +
  theme(legend.position = "none",
        axis.text.x  = element_text(colour = "black", size = 12),
        axis.text.y  = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14)) +
  labs(x = "Distance class (km)",
       y = "Correlation (Moran's I)") +
  theme(strip.text.x = element_text(size = 14),
        strip.background = element_blank())

## 7) Extended Data Fig. 2: 

load("sst.anomaly.rda")
load("chl.a.anomaly.rda")



# Wrangle species-specific prevalence data

spec_prev <- prev %>%
  group_by(siteName, genusSpecies) %>%
  summarise(mean.prev = mean(prevalence, na.rm = T)) %>%
  filter(mean.prev != 'NaN') %>%
  spread(genusSpecies, mean.prev)

disease.dist <- merge(disease.dist, spec_prev, by = "siteName")



# load summary satellite temp dat

load("summary_satellite_temps.rda")

# Test whether satellite measures of temp predict prevalence?

sat_predictors <- sat_temps %>%
  filter(sst != "NA") %>%
  filter(siteName != "corona.diablo") %>%
  rename(mean.Temperature = sst) %>%
  group_by(siteName) %>%
  summarise(deg.weeks.26 = length(mean.Temperature[which(mean.Temperature>26)]),
            deg.weeks.25 = length(mean.Temperature[which(mean.Temperature>25)]),
            mean_temp = mean(mean.Temperature, na.rm = T),
            var_temp = var(mean.Temperature, na.rm = T),
            sd_temp = sd(mean.Temperature, na.rm = T),
            max_temp = max(mean.Temperature, na.rm = T)) %>%
  right_join(disease.dist) %>%
  tidyr::gather(`Holacanthus\\npasser`, `Microspathodon\\ndorsalis`, `Myripristis\\nleiognathus`, `Stegastes\\nbeebei`, Sick, key = Species, value = Prevalence)

# Linear regression analysis
sick_regression <- sat_predictors %>%
  filter(Species == "Sick")
sick_mod <- lm(log10(Prevalence+1) ~ deg.weeks.26, data = sick_regression) %>%
  summary()

# Plot Extended Data Fig. 2

sick_regression_plot <- sat_predictors %>%
  filter(Species == "Sick") %>%
  ggplot(aes(x = deg.weeks.26, y = log10(Prevalence+1))) +
  geom_jitter(size = 3, width = .1, height = .1) +
  geom_smooth(method = "lm", se = F, col = "black") +
  theme_bw() +
  labs(x = "Weeks above 26 C",
       y = "Number species affected (log [x+1]") +
  theme(axis.text.x  = element_text(colour = "black", size = 12),
        axis.text.y  = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 14, face = "bold"))

## 8) Figure 5: Test for correlation between abundance and prevalence for each target species across sites

monitoring_prevalence <- prev %>%
  group_by(siteName, genusSpecies) %>%
  summarise(mean_abund = mean(total, na.rm = T),
            prevalence = mean(prevalence, na.rm = T))

monitoring_plot <- monitoring_prevalence %>%
  mutate(significance = as.character(fct_recode(genusSpecies, black = "Holacanthus\\npasser", black = "Myripristis\\nleiognathus", black = "Stegastes\\nbeebei", grey10 = "Microspathodon\\ndorsalis"))) %>%
  ggplot(aes(mean_abund, prevalence)) +
  geom_point(color = "black") +
  geom_smooth(aes(color = significance, linetype = significance2), method = "lm", se = F) +
  facet_wrap(~genusSpecies, scales = ("free"), ncol = 2) +
  theme_bw() + 
  scale_color_manual(values = c("black", NA)) +
  theme(legend.position = "none",
        axis.text.x  = element_text(colour = "black", size = 12),
        axis.text.y  = element_text(colour = "black", size = 12),
        axis.title.y = element_text(colour = "black", size = 14),
        axis.title.x = element_text(colour = "black", size = 14),
        strip.text.x = element_blank(),
        strip.background = element_blank()) +
  labs(x = expression("Mean abundance • 50m"^"-2"),
       y = "Disease prevalence")



## models

HOPA_mod = lm(prevalence~mean_abund, data = monitoring_prevalence[which(monitoring_prevalence$genusSpecies=="Holacanthus\\npasser"),])
summary(HOPA_mod)

STBE_mod = lm(prevalence~mean_abund, data = monitoring_prevalence[which(monitoring_prevalence$genusSpecies=="Stegastes\\nbeebei"),])
summary(STBE_mod)

MYLE_mod = lm(prevalence~mean_abund, data = monitoring_prevalence[which(monitoring_prevalence$genusSpecies=="Myripristis\\nleiognathus"),])
summary(MYLE_mod)

MIDO_mod = lm(prevalence~mean_abund, data = monitoring_prevalence[which(monitoring_prevalence$genusSpecies=="Microspathodon\\ndorsalis"),])
summary(MIDO_mod)

## 9) Fig. 6: Comparison among years of abundance for H. passer and S. beebei using data from stationary videos at Cuatro Hermanos East

# call to data
HOPA.abund <- 
  read_excel("passer_beebei_fishcams.xlsx") %>%
  mutate(Site = "ch.pequeno", 
         Camera     = as.factor(Camera), 
         Segment    = as.factor(Segment), 
         Prevalence = Sick_abundance / MaxN) %>%
  mutate(rel.abund = MaxN*10/visibility) %>%
  filter(Species == "HOPA") %>%
  filter(Date > "2003-06-01") %>%
  mutate(MaxN = rel.abund)

STBE.abund <- 
  read_excel("passer_beebei_fishcams.xlsx") %>%
  mutate(Site    = "ch.pequeno", 
         Camera  = as.factor(Camera), 
         Segment = as.factor(Segment), 
         Prevalence = Sick_abundance / MaxN) %>%
  mutate(rel.abund = MaxN*10/benthic_area) %>%
  filter(Species == "STBE") %>%
  filter(Date > "2003-06-01") %>%
  mutate(MaxN = rel.abund)

data_to_plot <-
  bind_rows(STBE.abund, HOPA.abund) %>%
  mutate(Species = fct_recode(Species, `Holacanthus passer` = "HOPA", 
                              `Stegastes beebei`   = "STBE")) %>%
  tidyr::gather(MaxN, Prevalence, key = Variable, value = Value) %>%
  group_by(Variable, Species, Date) %>%
  summarize(mean = na.omit(mean(Value)), 
            se   = sd(na.omit((Value)) / sqrt(length(na.omit(Value))))) %>%
  mutate(upper = mean + se, 
         lower = mean - se) 

##
##     
# create plot
prevalence_abundance <- data_to_plot %>%
  mutate(Date = as.Date(Date)) %>%
  ggplot(aes(Date, mean)) +
  geom_path(aes(Date, mean, linetype = Species)) +
  geom_point(aes(Date, mean), size = 3) +
  theme_bw() +
  geom_linerange(aes(ymax  = upper, 
                     ymin  = lower),
                 size = 1.4) +
  geom_point(shape = 1, size = 3) +
  facet_wrap(~Variable, nrow = 2, scales = "free_y") +
  #scale_colour_brewer(palette = "Paired") +
  labs(title = "", 
       x     = "", 
       #y      ="Prevalence +/- SE\\nAbundance (MaxN +/- SE)") +
       y      ="        Prevalence                                             Abundance (MaxN)") +
  scale_x_date(date_breaks = "6 months", date_labels = "%b-%y", limits = c(as.Date("2012-05-01"), as.Date("2017-09-01"))) +          
  theme(legend.position  = c(.29,.42), 
        legend.title     = element_blank(),
        legend.text = element_text(size = 11, face = "bold"),
        strip.background = element_blank(), 
        strip.text       = element_blank(),
        axis.text.x = element_text(colour="black", size=10, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(colour = "black", size = 14, face = "bold"))

## 10) Extended Data Table 3: Statistical analyses for differences in abundance for each species by year

#Test for normality and homogeneity of variance

#Stegastes beebei
leveneTest(STBE.abund$rel.abund, STBE.abund$Date)
qqnorm(STBE.abund$rel.abund)
qqline(STBE.abund$rel.abund)
STBE_mod = lm(rel.abund~as.factor(Date), data = STBE.abund)
summary(STBE_mod)
anova(STBE_mod)
TukeyHSD(aov(STBE_mod))


#Holacanthus passer

leveneTest(log(HOPA.abund$rel.abund), HOPA.abund$Date)
qqnorm(log(HOPA.abund$rel.abund))
qqline(log(HOPA.abund$rel.abund))
HOPA_mod = lm(log(rel.abund)~as.factor(Date), data = HOPA.abund)
summary(HOPA_mod)
anova(HOPA_mod)
TukeyHSD(aov(HOPA_mod))


## How much did populations drop during the outbreak?

ENSO <- data_to_plot %>%
  filter(Variable == "MaxN") %>%        
  filter(year(Date) == "2016") %>%
  group_by(Species, Date) %>%
  summarize(ENSO = mean(mean)) %>%
  group_by(Species) %>%
  summarize(ENSO = min(ENSO))

baseline <- data_to_plot %>%
  filter(Variable == "MaxN") %>%
  filter(year(Date) != "2016") %>%
  group_by(Species) %>%
  summarize(baseline = mean(mean))

mortality <- merge(ENSO, baseline, by = "Species") %>%
  mutate(baseline = as.numeric(baseline),
         ENSO = as.numeric(ENSO)) %>%
  mutate(mortality = (baseline - ENSO) / baseline)

