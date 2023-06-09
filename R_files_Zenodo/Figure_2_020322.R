# -------------------------------------------------------------------------
# Code for generating Figure 2
# -------------------------------------------------------------------------
setwd("~/Data Repository")

library(lme4)
library(raster)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(rgdal)
library(ggnewscale)
library(scales)
library(sf)



# ---------------------------------
# Read in data
# ---------------------------------
eBird_Filtered <- read_csv("MAD_Output_All_Filter_1.csv")

# Get grid-specific shifts
Grid_coefficients <- function (gridid){
  eBird_Filtered <- eBird_Filtered[eBird_Filtered$GridID == gridid ,]
  regression <- lmer(MAD ~ Year + (1|Common_Name), data = eBird_Filtered)
  MAD.Shift <- data.frame(fixef(regression))
  grid.id <- unique(eBird_Filtered$GridID)
  data.frame(grid.id, MAD.Shift[2,])
}

#Create dataframe of slopes == MAD Shifts
Grids <- unique(eBird_Filtered$GridID)
Grid_Shifts_lmer <- lapply(Grids, Grid_coefficients)   
Grid_Shifts_lmer <- data.frame(do.call(rbind, Grid_Shifts_lmer))
colnames(Grid_Shifts_lmer) <- c("GridID", "MAD_Shift_3")



# ----------------------------------------
# MAPPPING CODE
# ----------------------------------------
# Extract lat long from the grids
eBird_Filtered$lat <- str_extract_all(eBird_Filtered$GridID, "([-]*[0-9]+[x])")
eBird_Filtered$lat <- as.numeric(str_remove_all(eBird_Filtered$lat, "[x]"))
eBird_Filtered$long <- as.numeric(str_extract_all(eBird_Filtered$GridID, "[0-9]+$"))
eBird_coords <- data.frame(x = eBird_Filtered$long, y = eBird_Filtered$lat)
# Make spatialpointsdataframe
eBird_sp <- SpatialPointsDataFrame(coords = eBird_coords, data = eBird_Filtered)
proj4string(eBird_sp) <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
# Set radius
radius = 100000
# define the plot edges based upon the plot radius.
unique_Grids <- eBird_Filtered %>%
  select(GridID, lat, long) %>%
  distinct()

yPlus <- unique_Grids$lat+radius
xPlus <- unique_Grids$long+radius
yMinus <- unique_Grids$lat-radius
xMinus <- unique_Grids$long-radius
# calculate polygon coordinates for each plot centroid.
square = cbind(xMinus, yPlus,  # NW corner
               xPlus, yPlus,  # NE corner
               xPlus, yMinus,  # SE corner
               xMinus, yMinus, # SW corner
               xMinus, yPlus)  # NW corner again - close ploygon
# Extract the plot ID information
ID = unique_Grids$GridID

# create spatial polygons from coordinates
polys <- SpatialPolygons(mapply(function(poly, id)
{
  xy <- matrix(poly, ncol=2, byrow=TRUE)
  Polygons(list(Polygon(xy)), ID=id)
},
split(square, row(square)), ID),
proj4string = CRS(as.character("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")))
# Make into spdf
polys <- SpatialPolygonsDataFrame(Sr = polys, data = Grid_Shifts_lmer, match.ID = F)



# ---------------------------------------
# Make a map
# ---------------------------------------
# Convert polys
polys_sf <- st_as_sf(polys)
# Get Canada and USA
can1 <- getData('GADM', country="CAN", level=1) # provinces
us1 <- getData('GADM', country="USA", level=1)

# Lakes
lakes <- readOGR("Figure 2/ne_10m_lakes/ne_10m_lakes.shp") # file path to stored lake shapefile
lakes <- lakes[lakes$scalerank == 0,] # subset largest lakes - scalerank ranks lake size 0 = largest
# Just get the 5 great lakes
lakes <- lakes[lakes$name %in% c("Lake Erie", "Lake Michigan", "Lake Huron", "Lake Ontario", "Lake Superior"),]

# Specify a geographic extent for the map by defining the top-left and bottom-right geographic coordinates
# This can be modified
mapExtent <- rbind(c(-156, 80), c(-68, 19))

# Specify the required projection using a proj4 string
# Use http://www.spatialreference.org/ to find the required string
# ESRI 102008
newProj <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

# Project the map extent (first need to specify that it is longlat)
mapExtentPr <- spTransform(SpatialPoints(mapExtent,
                                         proj4string=CRS("+proj=longlat")),
                           newProj)

# Project other layers
can1Pr <- spTransform(can1, newProj)
us1Pr <- spTransform(us1, newProj)
lakesPr <- spTransform(lakes, newProj)

# Convert to sf
can1Pr_sf <- st_as_sf(can1Pr)
us1Pr_sf <- st_as_sf(us1Pr)
lakesPr_sf <- st_as_sf(lakesPr)

# Highlight provinces and states of interest
theseJurisdictions <- c("Ontario",
                        "QuÃ©bec",
                        "New Brunswick",
                        "Nova Scotia",
                        "Prince Edward Island",
                        "Connecticut",
                        "Delaware",
                        "District of Columbia",
                        "Illinois",
                        "Indiana",
                        "Kentucky",
                        "Maine",
                        "Maryland",
                        "Massachusetts",
                        "Michigan",
                        "New Hampshire",
                        "New Jersey",
                        "New York",
                        "Ohio",
                        "Pennsylvania",
                        "Rhode Island",
                        "Vermont",
                        "Virginia",
                        "West Virginia",
                        "Wisconsin")
# Add
can1Pr_sf$Highlight <- ifelse(can1Pr_sf$NAME_1 %in% theseJurisdictions, "black", "white")
us1Pr_sf$Highlight <- ifelse(us1Pr_sf$NAME_1 %in% theseJurisdictions, "black", "white")
lakesPr_sf$Highlight <- "white"

# Add a column symbolizing the shifts
polys_sf <- polys_sf %>% 
  mutate(MAD_Shift_Sign = ifelse(MAD_Shift_3 < 0, "Negative", "Positive"))

# Plot 
ggplot() +
  geom_sf(data = can1Pr_sf, aes(fill = Highlight), show.legend = FALSE) +
  geom_sf(data = us1Pr_sf, aes(fill = Highlight), show.legend = FALSE) +
  geom_sf(data = lakesPr_sf, aes(fill = Highlight), show.legend = FALSE) +
  scale_fill_manual(values = c("black", "white")) +
  new_scale_fill() +
  geom_sf(data = polys_sf, aes(fill = MAD_Shift_3), alpha = 0.9) +
  stat_sf_coordinates(data = filter(polys_sf, MAD_Shift_Sign == "Positive"), 
                      geom = "point", 
                      shape = "+", 
                      colour = "black",
                      size = 4) +
  stat_sf_coordinates(data = filter(polys_sf, MAD_Shift_Sign == "Negative"), 
                      geom = "point", 
                      shape = "-", 
                      colour = "black",
                      size = 4) +
  scale_fill_gradientn(colours = c("#66ccff", "#b3e6ff", "white", "#ff9999", "#ff3333"),
  values = rescale(c(max(polys_sf$MAD_Shift_3), 0.05, 0,  -0.37, min(polys_sf$MAD_Shift_3)))) +
  labs(fill = paste("Mean arrival date shift\\n(days/year)", sep = "")) +
  coord_sf(xlim = c(st_bbox(polys_sf)[1]*0.3, st_bbox(polys_sf)[3]*1.2),
           ylim = c(st_bbox(polys_sf)[2]*0.75, st_bbox(polys_sf)[4])*2) +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        text = element_text(family = "serif"))
ggsave(filename = "Figure_2.png",
       height = 15,
       width = 20,
       units = "cm",
       dpi = 300)
