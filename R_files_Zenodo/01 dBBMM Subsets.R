#### TITLE ####
# How do King Cobras move across a major highway? Unintentional
# wildlife crossing structures may facilitate movement.
#### AUTHORS ####
# Max Dolton Jones, Benjamin Michael Marshall, Samantha Nicole Smith, 
# Matt Crane, InÃªs Silva, Taksin Artchawakom, Pongthep Suwanwaree, 
# Surachit Waengsothorn, Wolfgang Wüster, Matt Goode, Colin Thomas Strine
#### YEAR ####
# 2021

## 01 dBBMM SCRIPT

## Subsetting dBBMM analysis for road crossing events

# Libraries -------------------------------------------------------------------------

library(ggplot2) 
library(scales) 
library(dplyr) 
library(move) 
library(adehabitatHR) 
library(ggspatial) 
library(rgeos) 
library(stringr)
library(ggmap)
library(readr)

# Create Folders ----------------------------------------------------------

loc.data <- paste0("./DATA/")
loc.new <- paste0("./Figures/New_dBBMM")
loc.shape <- paste0("./Shapefiles")

# Read Resources for dBBMM ----------------------------------------------

# Replace the file name below with the snake and time you are 
# interested in
all.data <- read_csv(file = paste0(loc.data, "AM054_20190504.csv"),
                     locale = locale(tz = "Asia/Bangkok"))

# Set dBBMM parameters --------------------------------------------------

set_grid.ext <- 5
set_dimsize <- 1000

set_tiff.width <- 2000
set_tiff.height <- 2000

move.obj <- move(x = all.data$x, y = all.data$y, proj = CRS("+init=epsg:32647"), 
                 time = all.data$datetime)

move.obj

# setting window size and margin. Both are odd, margin must be smaller than window.
ws <- 15
mrg <- 3

set_loc.error <- 5


# Pulls the current start time
ind.start <- Sys.time()

dbbmm <- brownian.bridge.dyn(object = move.obj, 
                             location.error = set_loc.error, 
                             margin = mrg, 
                             window.size = ws,  
                             ext = set_grid.ext, 
                             dimSize = set_dimsize,
                             verbose = FALSE)

print(paste("------ dBBMM computation time:", 
            round(difftime(Sys.time(), ind.start, units = "min"), 3), "mins"))

dbbmm

# Extract the motion variance from the subset
all.data$var <- getMotionVariance(dbbmm)

# Run the line below to save the motion variance data from subset
#write_csv(path = "./Data/AM054_20190504_var.csv", x = all.data)

# The next lines will create a simple motion variance plot
# for the dBBMM subset
var.plot <- ggplot(all.data) +
  geom_area(aes(x = datetime, y = var), alpha = 0.5, fill = "black") +
  geom_line(aes(x = datetime, y = var)) +
  theme_bw() +
  scale_size_manual(values = 0.6) +
  scale_x_datetime(labels = date_format("%Y-%m-%d"),
                   breaks = date_breaks("day")) +
  labs(x = "Date (YYYY-MM-DD)", y = (expression(paste(sigma^2,"m")))) +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = 2))

var.plot

ggsave(file = paste0("./Figures/New_dBBMM/Variance_AM054_20190504.png"), width = 200, height = 120,
       dpi = 600, units = "mm")
#ggsave(file = paste0("./Figures/New_dBBMM/Variance_AF010_20150601.pdf"), width = 200, height = 120,
 #      units = "mm")


# Now prepare the dBBMM output for plotting
dbbmm.sp <- as(dbbmm, "SpatialPixelsDataFrame")
dbbmm.sp.ud <- new("estUD", dbbmm.sp)
# @ is
dbbmm.sp.ud@vol = FALSE
dbbmm.sp.ud@h$meth = "dBBMM"
dbbmm.ud <- getvolumeUD(dbbmm.sp.ud, standardize = TRUE) 

poly.090 <- getverticeshr(dbbmm.ud, percent = 90)
poly.095 <- getverticeshr(dbbmm.ud, percent = 95)
poly.099 <- getverticeshr(dbbmm.ud, percent = 99)

fort.90 <- fortify(poly.090)
fort.95 <- fortify(poly.095)
fort.99 <- fortify(poly.099)

# Read in culvert data for plotting
cul.north <- read_csv(file = paste0(loc.data, "culvert_north.csv"))
cul.south <- read_csv(file = paste0(loc.data, "culvert_south.csv"))
culvert.data <- read_csv(file = paste0(loc.data, 
                                       "crossing_structure_data.csv"))

# specify the coordinate reference system in use
crs.proj <- CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Read in shapefiles for plotting
# Highway 304
hw.shp <- readOGR(paste0(loc.shape, "/304 cut.shp"))
hw.shp <- spTransform(hw.shp, crs.proj)

# Irrigation canal
klong <- readOGR(paste0(loc.shape, "/Klong.shp"))
klong <- spTransform(klong, crs.proj)

klong2 <- readOGR(paste0(loc.shape, "/Klong_2.shp"))
klong2 <- spTransform(klong2, crs.proj)

##cross map ----

# Now we have everything, we can plot the dBBMM subset
dBBMM.plot <- ggplot() +
  geom_path(data = hw.shp, aes(x = long, y = lat, group = group), linetype = 2, size = 1) +
  geom_polygon(data = poly.099,
                       aes(x = long, y = lat, group = group, fill = "orange", colour = "orange"),
                       alpha = 0.3) +
  geom_polygon(data = poly.095,
               aes(x = long, y = lat, group = group, fill = "orange", colour = "orange"),
               alpha = 0.5) +
  geom_polygon(data = poly.090,
               aes(x = long, y = lat, group = group, fill = "orange", colour = "orange"),
               alpha = 0.8) +
  geom_polygon(data = klong, aes(x = long, y = lat, group = group),
               alpha = 0.4, 
               fill = "blue", size = 0.5) +
  geom_polygon(data = klong2, aes(x = long, y = lat, group = group),
               alpha = 0.4, 
               fill = "blue", size = 0.5) +
  geom_point(data = all.data, aes(x = x, y = y), alpha = 0.65, size = 0.75) +
  geom_path(data = all.data, aes(x = x, y = y), alpha = 0.65, size = 1) +
  geom_point(data = cul.north, aes(x = x, y = y), size = 1) +
  geom_point(data = cul.south, aes(x = x, y = y), size = 1) +
  coord_equal(xlim = range(fort.99$long), y = range(fort.99$lat),
              expand = TRUE) +
  geom_text(data = culvert.data, aes(x = North_Easting, y = North_Northing,
                label = ID), nudge_y = 100, nudge_x = -50, hjust = 1,
            size = 5) +
  labs(title = paste0(all.data$id, ": ",  all.data$datetime[1], " - ", all.data$datetime[25]), x = "Easting", y = "Northing") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(hjust = 1, face = 2),
        plot.title = element_text(face = 4))

dBBMM.plot

ggsave(file = paste0("./Figures/New_dBBMM/dBBMM_AM054_20190504_New.png"), width = 200, height = 120,
       dpi = 600, units = "mm")
#ggsave(file = paste0("./Figures/New_dBBMM/dBBMM_AF010_20150601.pdf"), width = 200, height = 120,
 #      units = "mm")

#end