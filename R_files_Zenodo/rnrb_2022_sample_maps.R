# Introduction ----
# rnrb_2022_sample_maps.R
# Started by Libby Natola on 28 Feb 2022.

## snatched bits of https://slcladal.github.io/maps.html and https://gis.stackexchange.com/questions/224035/how-to-create-a-crisp-topographical-terrain-map-with-ggplot2 

# set working directory
setwd("~/Documents/UBC/Bioinformatics/rnrb/maps")

# # install packages (if necessary)
# install.packages("OpenStreetMap")
# install.packages("DT")
# install.packages("RColorBrewer")
# install.packages("mapproj")
# install.packages("sf")
# install.packages("RgoogleMaps")
# install.packages("scales")
# install.packages("rworldmap")
# install.packages("maps")
# install.packages("tidyverse")
# install.packages("rnaturalearth")
# install.packages("rnaturalearthdata")
# install.packages("rgeos")
# install.packages("ggspatial")
# install.packages("maptools")
# install.packages("leaflet")
# install.packages("sf")
# install.packages("tmap")
# install.packages("here")
# install.packages("rgdal")
# install.packages("scales")
# install.packages("flextable")
# # install package from github
# devtools::install_github("dkahle/ggmap", ref = "tidyup")
# # install klippy for copy-to-clipboard button in code chunks
# remotes::install_github("rlesur/klippy")
install.packages('ggnewscale')


# set options
options(stringsAsFactors = F)         # no automatic data transformation
options("scipen" = 100, "digits" = 4) # suppress math annotation
op <- options(gvis.plot.tag='chart')  # set gViz options

# load packages
library(OpenStreetMap)
library(DT)
library(RColorBrewer)
library(mapproj)
library(sf)
library(RgoogleMaps)
library(scales)
library(rworldmap)
library(maps)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(maptools)
library(leaflet)
library(sf)
library(tmap)
library(here)
library(rgdal)
library(scales)
library(flextable)
library(ggplot2)
library(cowplot)
library(rcartocolor)
# activate klippy for copy-to-clipboard button
klippy::klippy()

# add central data
manning_loc_Data <- read.table("~/Documents/UBC/Bioinformatics/rnrb/sample_info_man.txt", header=TRUE)
manning_loc_Data <- as.data.frame(manning_loc_Data)
manning_loc_Data$trans <- print("man")

# add north data
nbc_loc_Data <- read.table("~/Documents/UBC/Bioinformatics/rnrb/sample_info_nbc.txt", header = TRUE)
nbc_loc_Data <- as.data.frame(nbc_loc_Data)
nbc_loc_Data$trans <- print("nbc")

# add south data
caor_loc_Data <- read.table("~/Documents/UBC/Bioinformatics/rnrb/sample_info_caor.txt", header = TRUE)
caor_loc_Data <- as.data.frame(caor_loc_Data)
caor_loc_Data$trans <- print("caor")

# combined 
combined <- rbind(manning_loc_Data, nbc_loc_Data, caor_loc_Data)

# bounding boxes
man_bb <- as.data.frame(rbind(c(-122.0, 49.63),
               c(-119.90, 49.63),
               c(-119.90, 48.95),
               c(-122.0, 48.95),
               c(-122.0, 49.63)))
nbc_bb <- as.data.frame(rbind(c(-127.3,50.5),
                              c(-127.3,53),
                              c(-119.2,53),
                              c(-119.2,50.5),
                              c(-127.3,50.5)))
caor_bb <- as.data.frame(rbind(c(-123.85,40.1),
                              c(-123.85,45.9),
                              c(-117,45.9),
                              c(-117,40.1),
                              c(-123.85,40.1)))

levels(combined$trans)
combined$trans <- factor(combined$trans, levels = c("nbc","man","caor"))
                           
# make world map                  
world <- ne_countries(scale = "medium", returnclass = "sf")

# plot it
pdf("three_trans_map.pdf")
ggplot() +
  geom_sf(data = world, fill="gray96") + 
  geom_polygon(data = man_bb, aes(x=V1, y=V2), fill = NA, color = "#440154", size = 1.2) +
  geom_polygon(data = nbc_bb, aes(x=V1, y=V2), fill = NA, color = "#FDE725", size = 1.2) +
  geom_polygon(data = caor_bb, aes(x=V1, y=V2), fill = NA, color = "#21918C", size = 1.2) +
  labs( x = "Longitude", y = "Latitude", color = "Transect") +
  coord_sf(xlim = c(-127,-117), ylim = c(40,53), expand = T) +
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw() + theme(legend.position=c(-0.3, 0.8)) +
  geom_point(shape = 20,alpha = 0.7, size = 2.5, data=combined, aes(x=long, y=lat, color=trans)) + scale_color_manual(values = c("#FDE725", "#440154", "#21918C"), labels = c("North", "Central", "South"))
dev.off()

# get the data  
library(rasterVis)
library(raster)
can1 <- getData('GADM', country='CAN', level=1)
# subset it to BC
can1 <- subset(can1,NAME_1 %in% "British Columbia")

# get the centroid location and build out tiles from there?
can1c <- gCentroid(can1) %>% coordinates()
nbc1 <- getData("SRTM",lat=can1c[2],lon=can1c[1])

# Get tiles around that to fill out the whole plot
nbc2 <- getData("SRTM", lat = 54.75318, lon = -127.3) # the left side i think
nbc3 <- getData("SRTM",lat = 54.75318, lon = -119.2) # the right side what am i doing

# merge tiles, make the files smaller
nbc_topo <- merge(nbc1,nbc2,nbc3)
nbc_topo <- aggregate(nbc_topo, fact = 16) 

# turn raster to points in dataframe, rename columns to long lat and alt
nbc_topo.p  <-  rasterToPoints(nbc_topo)
nbc_topo_df <-  data.frame(nbc_topo.p)
colnames(nbc_topo_df) = c("lon", "lat", "alt")

# plot it
library(ggnewscale)
nbc_inset <- ggplot(nbc_topo_df, aes(lon,lat)) +
  geom_raster(aes(fill = alt)) + 
  scale_fill_gradientn(colours=c("grey45","white"), limits=c(0,4000)) +
  coord_sf(xlim = c(-127.3,-119.2), ylim = c(50.5,53), expand = T) + 
  new_scale("fill") +
  geom_point(data=nbc_loc_Data, colour = "black", alpha = 0.7, size = 3.5, aes(x=long, y=lat, shape = spp, fill=spp)) + 
  scale_fill_manual(name = "Phenotypic species", values = c("#DE73FF", "#FF4343", "#4343FF"), 
                    labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  scale_shape_manual(name = "Phenotypic species", values = c(23, 22, 21), labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  labs( x = "Longitude", y = "Latitude", color = "Phenotypic Species ID", fill="Altitude") +
  theme(panel.border = element_rect(colour = "#FDE725", fill=NA, size=2))

pdf("nbc_inset.pdf")
nbc_inset
dev.off()

ggplot() +
  geom_sf(data = world, fill="gray96") + 
  coord_sf(xlim = c(-127.3,-119.2), ylim = c(50.5,53), expand = T) + 
  geom_point(shape = 20,alpha = 0.7, size = 2.5, data=nbc_loc_Data, aes(x=long, y=lat, color=spp)) + 
  theme(panel.border = element_rect(colour = "#FDE725", fill=NA, size=2))

# man data wrangling 
man1 <- getData("SRTM", lat = 49.63, lon = -122) # the left side i think
man2 <- getData("SRTM",lat = 49.63, lon = -120) # the right side what am i doing
man_topo <- merge(man1, man2)
man_topo <- aggregate(man_topo, fact = 5) 
man_topo.p  <-  rasterToPoints(man_topo)
man_topo_df <-  data.frame(man_topo.p)
colnames(man_topo_df) = c("lon", "lat", "alt")

man_inset <- ggplot(man_topo_df, aes(lon,lat)) +
  geom_raster(aes(fill = alt)) + 
  scale_fill_gradientn(colours=c("grey45","white"), limits=c(0,4000)) +
  coord_sf(xlim = c(-122,-119.9), ylim = c(49.63,48.95), expand = T) + 
  new_scale("fill") +
  geom_point(data=manning_loc_Data, colour = "black", alpha = 0.7, size = 3.5, aes(x=long, y=lat, shape = spp, fill=spp)) + 
  scale_fill_manual(name = "Phenotypic species", values = c("#DE73FF", "#FF4343", "#4343FF"), 
                     labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  scale_shape_manual(name = "Phenotypic species", values = c(23, 22, 21), labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  labs( x = "Longitude", y = "Latitude", color = "Phenotypic Species ID", fill="Altitude") +
  theme(panel.border = element_rect(colour = "#440154", fill=NA, size=2))

pdf("man_inset.pdf")
man_inset
dev.off()

ggplot() +
  geom_sf(data = world, fill="gray96") + 
  coord_sf(xlim = c(-122,-119.9), ylim = c(49.63,48.95), expand = T) + 
  geom_point(shape = 20,alpha = 0.7, size = 2.5, data=manning_loc_Data, aes(x=long, y=lat, color=spp)) + 
  theme(panel.border = element_rect(colour = "#440154", fill=NA, size=2))

# california and oregon data wrangling
us1 <- getData('GADM', country='USA', level=1)
us1 <- subset(us1,NAME_1 %in% c("California", "Oregon"))
caor1 <- getData("SRTM", lat = 40, lon = -124) # left very bottom
caor2 <- getData("SRTM",lat = 40, lon = -120) # right very bottom
caor3 <- getData("SRTM",lat = 45, lon = -124) # left center
caor4 <- getData("SRTM",lat = 45, lon = -120) # right center
caor5 <- getData("SRTM",lat = 46, lon = -124) # left top
caor6 <- getData("SRTM",lat = 46, lon = -120) # right top
caor_topo <- merge(caor1, caor2, caor3, caor4, caor5, caor6)
caor_topo <- aggregate(caor_topo, fact = 16) 
caor_topo.p  <-  rasterToPoints(caor_topo)
caor_topo_df <-  data.frame(caor_topo.p)
colnames(caor_topo_df) = c("lon", "lat", "alt")

caor_inset <- ggplot(caor_topo_df, aes(lon,lat)) +
  geom_raster(aes(fill = alt)) + scale_fill_gradientn(colours=c("grey45","white"), limits=c(0,4000)) +
  coord_sf(xlim = c(-123.85,-117), ylim = c(45.9,40.1), expand = T) + 
  new_scale("fill") +
  geom_point(data=caor_loc_Data, colour = "black", alpha = 0.7, size = 3.5, aes(x=long, y=lat, shape = spp, fill=spp)) + 
  scale_fill_manual(name = "Phenotypic species", values = c("#DE73FF", "#FF4343", "#4343FF"), 
                    labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  scale_shape_manual(name = "Phenotypic species", values = c(23, 22, 21), labels = c("Hybrid", "Red-breasted", "Red-naped")) +
  labs( x = "Longitude", y = "Latitude", color = "Phenotypic Species ID", fill="Altitude") +
  theme(panel.border = element_rect(colour = "#21918C", fill=NA, size=2))

pdf("caor_inset.pdf")
caor_inset
dev.off()

ggplot() +
  geom_sf(data = world, fill="gray96") + 
  coord_sf(xlim = c(-123.85,-117), ylim = c(45.9,40.1), expand = T) + 
  geom_point(shape = 20,alpha = 0.7, size = 2.5, data=caor_loc_Data, aes(x=long, y=lat, color=spp)) + scale_color_manual(values = c("#DE73FF", "#4343FF", "#FF4343"), labels = c("Hybrid", "Red-naped", "Red-breasted")) +
  theme(panel.border = element_rect(colour = "#21918C", fill=NA, size=2))
