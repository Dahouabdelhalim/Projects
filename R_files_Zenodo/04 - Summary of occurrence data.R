# 04 - Summary of occurrence data

# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sp)
library(pracma)
library(rgeos)
library(scico)
library(lubridate)

# Read data and make folders ----------------------------------------------

set.seed(1)

loc.gbifdata <- "./Data/GBIFData/"
# loc.figure <- "./Figures/"

gbif.data <- read.csv(file = paste0(loc.gbifdata, "7_compiled_GBIFData_coordcleanComplete.csv"),
                      stringsAsFactors = FALSE)

flr.data <- read.csv(file = "./Data/Flickr_results_clean.csv",
                     stringsAsFactors = FALSE)

map <- map_data("world")

# Generate crop areas for continental summary -----------------------------

x_coord <- c(-180, 0, -40, -75,-90, -180)
y_coord <- c(90, 90, 15, 15, 0, 0)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

namerica <- sps


x_coord <- c(-180, 0, -40, -75, -90, -180)
y_coord <- c(-60, -60, 15, 15, 0, 0)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

samerica <- sps

x_coord <- c(-25, -5, 10, 18, 29, 45, 55, 60, 90, 0, -40)
y_coord <- c(25, 36, 38, 35, 35, 11, 15, 15, -60, -60, 20)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

africa <- sps

x_coord <- c(-35, 0, 55, 50, 30, 25, 25, 18, 10, -5)
y_coord <- c(20, 90, 90, 42, 42, 40, 35, 35, 38, 36)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

europe <- sps

x_coord <- c(180, 180, 55, 50, 30, 25, 25, 29, 45, 55, 60, 70, 130, 130)
y_coord <- c(20, 90, 90, 42, 42, 40, 35, 35, 11, 15, 15, -15, -10, 0)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

asia <- sps

x_coord <- c(130, 130, 180, 180, 90, 70)
y_coord <- c(-10, 0, 20, -60, -60, -15)
xym <- cbind(x_coord, y_coord)
p <- Polygon(xym)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))
# sps <- rgeos::gBuffer(sps, width = 3)
# plot(sps)
proj4string(sps) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

oceania <- sps

ggplot(map)+
  geom_polygon(aes(x = long, y = lat, group = group)) +
  geom_polygon(data = namerica, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "red") +
  geom_polygon(data = samerica, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "orange") +
  geom_polygon(data = africa, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "lightblue") +
  geom_polygon(data = europe, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "blue") +
  geom_polygon(data = asia, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "purple") +
  geom_polygon(data = oceania, aes(x = long, y = lat, group = group),
               alpha = 0.25, fill = "green") +
  coord_quickmap() +
  scale_x_continuous(breaks = seq(-180, 180, 10)) +
  scale_y_continuous(breaks = seq(-90, 90, 10)) +
  theme_bw()

#
spl <- rbind(namerica, samerica, africa, europe, asia, oceania, makeUniqueIDs = TRUE)
polyData <- data.frame("con" = c("NAmerica", "SAmerica", "Africa", "Europe", "Asia", "Oceania"))
row.names(polyData) <- c("1", "11", "12", "13", "14", "15")

contintents <- SpatialPolygonsDataFrame(Sr = spl,
                         data = polyData)


# Generating summary table ------------------------------------------------

sp <- SpatialPoints(coords = 
                cbind(gbif.data$lon, gbif.data$lat), 
              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

loc_cont <- over(x = sp, y = contintents)

gbif.data$continent <- loc_cont$con

occ_summaryTABLE <- gbif.data %>% 
  group_by(continent, species) %>% 
  count() %>% 
  ungroup() %>%
  group_by(continent) %>%
  summarise(num_sp = n(),
            occ = sum(n),
            m_occ_sp = mean(n, na.rm = TRUE),
            se_occ_sp = std_err(n)) %>% 
  filter(!is.na(continent))

occ_summaryTABLE

sum(occ_summaryTABLE$num_sp)
sum(occ_summaryTABLE$occ)

land <- readRDS(file = "landmass.Rds")

# plot(land)

# Projecting continents for area estimations ------------------------------

# Albers Equal Area Conic
prj_Namerica <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
prj_Samerica <- CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs ")
prj_Nasia <- CRS("+proj=aea +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
prj_Africa <- CRS("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
prj_Europe <- CRS("+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs ")
prj_Sasia <- CRS("+proj=aea +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

land_namerica <- gIntersection(land, namerica)
land_samerica <- gIntersection(land, samerica)
land_asia <- gIntersection(land, asia)
land_africa <- gIntersection(land, africa)
land_europe <- gIntersection(land, europe)
land_oceania <- gIntersection(land, oceania)

proj_namerica <- spTransform(x = land_namerica, CRSobj = prj_Namerica)
proj_samerica <- spTransform(x = land_samerica, CRSobj = prj_Samerica)
proj_asia <- spTransform(x = land_asia, CRSobj = prj_Nasia)
proj_africa <- spTransform(x = land_africa, CRSobj = prj_Africa)
proj_europe <- spTransform(x = land_europe, CRSobj = prj_Europe)
proj_oceania <- spTransform(x = land_oceania, CRSobj = prj_Sasia)

# km^2
areas <- c(
  gArea(proj_africa)/1000,
  gArea(proj_asia)/1000,
  gArea(proj_europe)/1000,
  gArea(proj_namerica)/1000,
  gArea(proj_oceania)/1000,
  gArea(proj_samerica)/1000)

occ_summaryTABLE$cont_area_km2 <- areas

occ_summaryTABLE2 <- occ_summaryTABLE %>% 
  mutate(area_millkm2 = cont_area_km2/1000000000,
         occ_100km2 = occ/cont_area_km2 * 100*100) %>% 
  select(-cont_area_km2)

FINALTABLE <- occ_summaryTABLE2 %>%
  mutate(continent = c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")) %>% 
  rename("Continent" = continent, "# Species" = num_sp, "# Occurrences" = occ,
         "Mean occurrences per species" = m_occ_sp, "SE" = se_occ_sp,
         "Area (mill km2)" = area_millkm2, "Occurrences per 100 km2" = occ_100km2)

write.csv(x = FINALTABLE, row.names = FALSE, file = "./Tables/TABLE - SampleSizeSummary.csv")

# Map of GBIF records -----------------------------------------------------

(world_gbif_cont <- 
    ggplot(map)+
    geom_polygon(aes(x = long, y = lat, group = group), alpha = 0) +
    geom_polygon(data = land_namerica, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[1]) +
    geom_polygon(data = land_samerica, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[2]) +
    geom_polygon(data = land_africa, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[3]) +
    geom_polygon(data = land_europe, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[4]) +
    geom_polygon(data = land_asia, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[5]) +
    geom_polygon(data = land_oceania, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[6]) +
    geom_text(data = data.frame(x = seq(-180,180,30), y = rep(0,180*2/30+1)),
              aes(x = x, y = y),
              label = seq(-180,180,30), vjust = 0.5, hjust = 0.5,
              size = 3, colour = "grey25") +
    geom_text(data = data.frame(y = seq(-80, 80, 20)[-5], x = rep(0,80*2/20+1)[-5]),
              aes(x = x, y = y),
              label = seq(-80, 80, 20)[-5], vjust = 0.5, hjust = 0.5,
              size = 3, colour = "grey25") +
    geom_point(data = gbif.data, aes(x =  lon, y = lat),
               pch = '.', colour = "black") +
    coord_map("mollweide") +
    scale_x_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180)) +
    scale_y_continuous(breaks = seq(-80, 80, 20)) +
    theme_bw() +
    labs(x = "", y = "") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())
)

ggsave(file = "./Figures/FIGURE - Map of GBIF and Contintents.png",
       plot = world_gbif_cont, width = 18, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(file = "./Figures/FIGURE - Map of GBIF and Contintents.pdf",
       plot = world_gbif_cont, width = 18, height = 12, units = "cm",
       device = "pdf")


# Comparable Flickr data map ----------------------------------------------

(world_flickr_cont <- 
    ggplot(map)+
    geom_polygon(aes(x = long, y = lat, group = group), alpha = 0) +
    geom_polygon(data = land_namerica, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[1]) +
    geom_polygon(data = land_samerica, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[2]) +
    geom_polygon(data = land_africa, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[3]) +
    geom_polygon(data = land_europe, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[4]) +
    geom_polygon(data = land_asia, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[5]) +
    geom_polygon(data = land_oceania, aes(x = long, y = lat, group = group),
                 alpha = 1, fill = scico(n = 6, palette = "roma")[6]) +
    geom_text(data = data.frame(x = seq(-180,180,30), y = rep(0,180*2/30+1)),
              aes(x = x, y = y),
              label = seq(-180,180,30), vjust = 0.5, hjust = 0.5,
              size = 3, colour = "grey25") +
    geom_text(data = data.frame(y = seq(-80, 80, 20)[-5], x = rep(0,80*2/20+1)[-5]),
              aes(x = x, y = y),
              label = seq(-80, 80, 20)[-5], vjust = 0.5, hjust = 0.5,
              size = 3, colour = "grey25") +
    geom_point(data = flr.data, aes(x =  longitude, y = latitude),
               pch = '.', colour = "black") +
    coord_map("mollweide") +
    scale_x_continuous(breaks = seq(-180, 180, 30), limits = c(-180, 180)) +
    scale_y_continuous(breaks = seq(-80, 80, 20)) +
    theme_bw() +
    labs(x = "", y = "") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
)

ggsave(file = "./Figures/FIGURE - Map of Flickr and Contintents.png",
       plot = world_flickr_cont, width = 18, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(file = "./Figures/FIGURE - Map of Flickr and Contintents.pdf",
       plot = world_flickr_cont, width = 18, height = 12, units = "cm",
       device = "pdf")

# Temporal trends in data ------------------------------------------------

gbif_trend <- gbif.data %>% 
  group_by(year) %>% 
  count() %>% 
  mutate(source = "GBIF")

flr_trend <- flr.data %>% 
  group_by(year(as.POSIXct(datetaken))) %>% 
  count() %>% 
  mutate(source = "Flickr") %>% 
  rename("year" = `year(as.POSIXct(datetaken))`)

trend_data <- rbind(gbif_trend, flr_trend)

(trend_plot <- 
    ggplot(trend_data) +
    geom_point(aes(x = year, y = n, colour = source, shape = source)) +
    # scale_y_sqrt() +
    scale_y_log10() +
    scale_x_continuous(breaks = seq(1700, 2020, 20)) +
    theme_bw() +
    labs(colour = "Source", shape = "Source",
         x = "Year", y = "Number of Records (log)") +
    theme(legend.position = c(0.08,0.9),
          legend.background = element_blank(),
    ) +
    scale_colour_scico_d(palette = "roma", begin = 0.075, end = 0.9,
                         direction = -1)
)
  
ggsave(file = "./Figures/SUPPLEMENTARY FIGURE 4 - Trend of Records.png",
       plot = trend_plot, width = 18, height = 12, units = "cm",
       dpi = 300, device = "png")
ggsave(file = "./Figures/SUPPLEMENTARY FIGURE 4 - Trend of Records.pdf",
       plot = trend_plot, width = 18, height = 12, units = "cm",
       device = "pdf")
