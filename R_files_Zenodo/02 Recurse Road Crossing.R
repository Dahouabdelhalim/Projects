#### TITLE ####
# How do King Cobras move across a major highway? Unintentional
# wildlife crossing structures may facilitate movement.
#### AUTHORS ####
# Max Dolton Jones, Benjamin Michael Marshall, Samantha Nicole Smith, 
# Matt Crane, InÃªs Silva, Taksin Artchawakom, Pongthep Suwanwaree, 
# Surachit Waengsothorn, Wolfgang Wüster, Matt Goode, Colin Thomas Strine
#### YEAR ####
# 2021

## 02 RECURSE SCRIPT

# Libraries -------------------------------------------------------------------------

library(readr)
library(rgdal)
library(sp)
library(move)
library(recurse)
library(ggplot2)
library(scico)
library(dplyr)
library(ggspatial)
library(stringr)

# Create Folders ----------------------------------------------------------

loc.data <- paste0("./Data/")
loc.fig <- paste0("./Figures/")
loc.shape <- paste0("./Shapefiles")

## Highway 304 crossings ---- 

# Read Resources for Recurse ----------------------------------------------

# read in the csv file with our telemetry data

full.recurse.data <- read.csv(file = paste0(loc.data, "Ophiophagus hannah 2014_2020 recurse new.csv"),
                     header = TRUE, stringsAsFactors = FALSE, sep = ",")



# just make sure that R is reading the datetime column (t) as a date and time

full.recurse.data$t <- as.POSIXct(full.recurse.data$t, 
                                format = "%Y-%m-%d %H:%M:%S",
                                tz = "Asia/Bangkok")

# the class now should be POSIXct
class(full.recurse.data$t)

# specify the coordinate reference system in use
crs.proj <- CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# read in our shape files
# 304northside is used for determining the crossings
# across the Highway 304
northside304 <- readOGR(paste0(loc.shape, "/304northside.shp"))
northside304 <- spTransform(northside304, crs.proj)
northside304 <- SpatialPolygons(northside304@polygons,
                                proj4string = northside304@proj4string)

# South_304NW-TLC is used for determining crossings across
# the 304NW-TLC road
southNW_TLC <- readOGR(paste0(loc.shape, "/South_304NW-TLC.shp"))
southNW_TLC <- spTransform(southNW_TLC, crs.proj)
southNW_TLC <- SpatialPolygons(southNW_TLC@polygons,
                                proj4string = southNW_TLC@proj4string)

# these three csv files contain dates for when three snakes were not tracked
# we will need these later for plotting
miss.data.6 <- read_csv(file = paste0(loc.data, "missing_data_6.csv"))
miss.data.7 <- read_csv(file = paste0(loc.data, "missing_data_7.csv"))
miss.data.10 <- read_csv(file = paste0(loc.data, "missing_data_10.csv"))

# make sure R understands the date and time for the missing data
miss.data.6$begin <- as.POSIXct(miss.data.6$begin, 
                                   format = "%Y-%m-%dT%H:%M:%SZ",
                                   tz = "Asia/Bangkok")

miss.data.6$end <- as.POSIXct(miss.data.6$end, 
                                 format = "%Y-%m-%dT%H:%M:%SZ",
                                 tz = "Asia/Bangkok")

miss.data.7$begin <- as.POSIXct(miss.data.7$begin, 
                                   format = "%Y-%m-%dT%H:%M:%SZ",
                                   tz = "Asia/Bangkok")

miss.data.7$end <- as.POSIXct(miss.data.7$end, 
                                 format = "%Y-%m-%dT%H:%M:%SZ",
                                 tz = "Asia/Bangkok")

miss.data.10$begin <- as.POSIXct(miss.data.10$begin, 
                                    format = "%Y-%m-%dT%H:%M:%SZ",
                                    tz = "Asia/Bangkok")

miss.data.10$end <- as.POSIXct(miss.data.10$end, 
                                  format = "%Y-%m-%dT%H:%M:%SZ",
                                  tz = "Asia/Bangkok")

### Reformat data ---- 

names(full.recurse.data)
full.recurse.data <- full.recurse.data[,c(2,3,4,5)]
names(full.recurse.data) <- c("x", "y", "t", "id")

full.recurse.data$id <- as.factor(full.recurse.data$id)

## Resident time in polygon ------
all.visits.north <- getRecursionsInPolygon(full.recurse.data, northside304, timeunits = "days")

# for a summary of the entrance and exit times for our animals:
all.visits.north$revisitStats

# turn these stats into an object for saving
visits.north <- all.visits.north$revisitStats

# save revisit stats
write_csv(path = paste0(loc.data, "highway 304 crossing data.csv"), x = visits.north)

# create an object with the snake IDs

snakes <- c("AM006", "AM007",
            "AF010", "JM013", "AM015", "AF017", 
            "AM018", "JM019", "AM024", "JM025",
            "AM026", "JF027", "JM034", "AM054", "JF055", 
            "AF056", "AF058", "AM059", "AF086", "AF096", "AF099")

## Full track time from data ---- 
snake.track.time <- NULL
for(snk in snakes){
  snk.id <- full.recurse.data[full.recurse.data$id == snk,]
  track.time <- cbind(as.character(snk.id$id[1]), as.character(min(snk.id$t)),
                      as.character(max(snk.id$t)))
  snake.track.time <- rbind(snake.track.time, track.time)
}
snake.track.time <- as.data.frame(snake.track.time)
names(snake.track.time) <- c("id", "startdate", "enddate")
snake.track.time$startdate <- as.POSIXct(strptime(snake.track.time$startdate,
                                                  format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
snake.track.time$enddate <- as.POSIXct(strptime(snake.track.time$enddate,
                                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

write_csv(path = paste0(loc.data, "visit.north.new.plot.csv"), x = visits.north)
write_csv(path = paste0(loc.data, "snake.track.time.csv"), x = snake.track.time)
visits.north.new <- read.csv(file = paste0(loc.data, "visit.north.new.plot.csv"),
                              header = TRUE, stringsAsFactors = FALSE, sep = ",")
snake.track.time.new <- read.csv(file = paste0(loc.data, "snake.track.time.csv"),
                              header = TRUE, stringsAsFactors = FALSE, sep = ",")

class(snake.track.time.new$startdate)

visits.north.new$entranceTime <- as.POSIXct(visits.north.new$entranceTime, 
                                         format = "%Y-%m-%dT%H:%M:%SZ",
                                         tz = "Asia/Bangkok")

visits.north.new$exitTime <- as.POSIXct(visits.north.new$exitTime, 
                                            format = "%Y-%m-%dT%H:%M:%SZ",
                                            tz = "Asia/Bangkok")

snake.track.time.new$startdate <- as.POSIXct(snake.track.time.new$startdate, 
                                            format = "%Y-%m-%dT%H:%M:%SZ",
                                            tz = "Asia/Bangkok")

snake.track.time.new$enddate <- as.POSIXct(snake.track.time.new$enddate, 
                                        format = "%Y-%m-%dT%H:%M:%SZ",
                                        tz = "Asia/Bangkok")


color.plot.1 <- c("#333333", "#333333", "#333333", "#333333")

# Plotting all the road crossing events ---- 
visit.new <- visits.north.new %>%
  mutate(sexage = str_extract(id, "^.."))

track.new <- snake.track.time.new %>%
  mutate(sexage = str_extract(id, "^.."))

# now we have our visiting times (road crossings) and the full
# tracking data for each snake, we can plot the crossings along
# each snake's tracking period
  
 cross.graph <- ggplot() +
  geom_segment(data = visit.new,
               aes(x = entranceTime, xend = exitTime, y = id, yend = id, colour = sexage),
               size = 5) +
  geom_segment(data = track.new,
               aes(x = startdate, xend = enddate, y = id, yend = id, colour = sexage),
               alpha = 0.4, size = 5) +
  scale_color_manual(values = color.plot.1) +
  geom_rect(data = miss.data.6, aes(xmin = begin, xmax = end,
                                       ymin = 0.74, ymax = 1.24),
                                       fill = "#FF3300", colour = "#FF3300") +
  geom_rect(data = miss.data.7, aes(xmin = begin, xmax = end,
                                       ymin = 1.74, ymax = 2.24),
                                       fill = "#FF3300", colour = "#FF3300") +
  geom_rect(data = miss.data.10, aes(xmin = begin, xmax = end,
                                        ymin = 2.74, ymax = 3.24),
                                        fill = "#FF3300", colour = "#FF3300") +
  coord_cartesian() +
  scale_x_datetime(date_breaks = "3 months", date_labels = "%Y-%m") +
  scale_y_discrete(limits = c("AM006", "AM007",
                              "AF010", "JM013", "AM015", "AF017", 
                              "AM018", "JM019", "AM024", "JM025",
                              "AM026", "JF027", "JM034", "AM054", "JF055", 
                              "AF056", "AF058", "AM059", "AF086", "AF096", "AF099")) +
  theme_bw() +
  labs(title = "Road crossing events", y = "Snake", x = "Date (Year - Month)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

cross.graph

ggsave(file = paste0("./Figures/304 cross graph new 2.png"), width = 200, height = 180,
       dpi = 600, units = "mm")
ggsave(file = paste0("./Figures/Figure - highway road crossing.pdf"), width = 200, height = 180,
       units = "mm")

## 304NW-TLC crossings and nesting ---- 

# The above sections of script can be used to determine
# crossings across the 304NW-TLC road by using the specified
# shapefile

# Below is the data and code used for creating Fig.4

# Movement data for just adult females used in plot
data.var <- read_csv(file = paste0(loc.data, "data with movvar_AF_2020.csv"),
                     locale = locale(tz = "Asia/Bangkok"))

# Data for before the snakes interacted with the road
bef.crossing <- read_csv(file = paste0(loc.data, "/before_crossing.csv"))

# Data for after the snakes have first crossed the road
south.TLC <- read_csv(file = paste0(loc.data, "/south_of_TLC.csv"))

# Data for when the snakes are actively nesting
nesting <- read_csv(file = paste0(loc.data, "/female_nesting.csv"))

# This is the data for when the snakes moved back over the
# road following nesting
north.TLC <- read_csv(file = paste0(loc.data, "/north_of_TLC.csv"))

library(scales)

movevarplot <- ggplot() +
  geom_rect(data = nesting,
            aes(xmin = entranceTime, xmax = exitTime, ymin = -100, ymax = -60),
            fill = "#FF9900", colour = NA) +
  geom_rect(data = north.TLC,
            aes(xmin = begin, xmax = end, ymin = -180, ymax = -140),
            fill = "#33CCFF", colour = NA, alpha = 0.4) +
  geom_rect(data = bef.crossing,
            aes(xmin = begin, xmax = end, ymin = -180, ymax = -140),
            fill = "#33CCFF", colour = NA, alpha = 0.4) +
  geom_rect(data = south.TLC,
            aes(xmin = begin, xmax = end, ymin = -180, ymax = -140),
            fill = "#33CCFF", colour = NA) +
  #geom_line(data = data.var, aes(x = datetime, y = var), size = 1) +
  facet_grid(id ~., scales = "free_y", space = "free_y") +
  scale_y_continuous(breaks = seq(0, 600, 100),
                     minor_breaks = seq(0, 600, 100)) +
  scale_x_datetime(labels = date_format("%b"), 
                   breaks = date_breaks("month"),
                   # expand = c(0.1, 0),
                   expand = c(0, 0),
                   position = "top") +  
  geom_text(data = data.var, aes(label = id), x = Inf, y = Inf,
            hjust = 1.1, vjust = 1.3, fontface = 4, size = 6) +
  coord_cartesian(xlim = as.POSIXct(c("2020-03-01 00:00:00",
                                      "2021-01-15 00:00:00"),
                                    format = "%Y-%m-%d %H:%M:%S")) +
  labs(x="Date (month)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(hjust = 1, face = 2, size = 14),
        axis.title.y = element_text(vjust = 1, angle = 0, face = 2),
        legend.position = "none",        
        strip.background = element_blank(),
        strip.text.y = element_blank())

movevarplot

ggsave(file = paste0("./Figures/South_TLC_AF.png"),
       width = 300, height = 140,
       dpi = 600, units = "mm")
ggsave(file = paste0("./Figures/South_TLC_AF.pdf"),
       width = 300, height = 140,
       units = "mm")


#end
