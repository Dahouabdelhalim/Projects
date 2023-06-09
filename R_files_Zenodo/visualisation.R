# Load packages
library(rgeos)
library(ggplot2)
library(dplyr)
library(raster)
library(maptools)
library(foreach)
library(gridExtra)
library(htmlTable)
library(tableHTML)

# Create folder for figures
if(!dir.exists("figs")){
  dir.create("figs")
}

# Source functions:
source("custom_functions.R")

# Load in spatial layers:
load(file = "spatial_rasters")

# Create maps and stacked plots:
boundary <- shapefile("Kapalga_outline_9453.shp")
crs(boundary) <- "+proj=utm +zone=53 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
boundary <- boundary[ , -(1:26)]
boundary$id <- as.character(1)

fire_ts15 <- rasterToPolygons(Early_fire[[15]])
fire_ts15 <- fire_ts15[fire_ts15$layer.15 == 1, ]
fire_ts15 <- fire_ts15[ , -1]
fire_ts15$id <- as.character(2:(nrow(fire_ts15) + 1))

fire_ts17 <- rasterToPolygons(Late_fire[[17]])
fire_ts17 <- fire_ts17[fire_ts17$layer.17 == 1, ]
fire_ts17 <- fire_ts17[ , -1]
fire_ts17$id <- as.character(2:(nrow(fire_ts17) + 1))

fire_ts1517 <- spRbind(fire_ts15, fire_ts17)

shift_matrix <- matrix(c(1.5, 1.3, 0, 0.5), 2, 2)

no_fire <- transform_polygons(boundary, shift_matrix)
early_fire <- transform_polygons(fire_ts15, shift_matrix)
late_fire <- transform_polygons(fire_ts17, shift_matrix)
all_late_fires <- transform_polygons(fire_ts1517, shift_matrix)


png('figs/burn_diagrams.png',
    pointsize = 8,
    res = 300,
    width = 1600,
    height = 1000)

# plot early/late fires
p_early_late <- ggplot() + 
  geom_polygon(data = no_fire, aes(x = x, y = y, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = late_fire, aes(x = x, y = y + 12000, group = id), fill = "black", color = NA) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 12000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 24000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = early_fire, aes(x = x, y = y + 36000, group = id), fill = "black", color = NA) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 36000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 48000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 60000, group = id), fill = NA, color = "black", size = 0.1) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.margin=unit(c(0, -1.0, 0, -0.8), "cm")) +
  annotate("text", x = 11530000, y = 4311000, size = 2, color = "black", hjust = 0, label = "November/December") +
  annotate("text", x = 11530000, y = 4323000, size = 2, color = "black", hjust = 0, label = "September/October") +
  annotate("text", x = 11530000, y = 4335000, size = 2, color = "black", hjust = 0, label = "July/August") +
  annotate("text", x = 11530000, y = 4347000, size = 2, color = "black", hjust = 0, label = "May/June") +
  annotate("text", x = 11530000, y = 4359000, size = 2, color = "black", hjust = 0, label = "March/April") +
  annotate("text", x = 11530000, y = 4371000, size = 2, color = "black", hjust = 0, label = "January/February") +
  annotate("text", x = 11525000, y = 4298000, size = 3, color = "black", hjust = 0, label = "Early/Late Fires") +
  annotate("segment", x = 11543000, xend = 11552000, y = 4352000, yend = 4352000, colour = "black", size = 0.2) +
  annotate("segment", x = 11543000, xend = 11552000, y = 4328000, yend = 4328000, colour = "black", size = 0.2) +
  annotate("segment", x = 11552000, xend = 11552000, y = 4328000, yend = 4352000, colour = "black", size = 0.2) +
  annotate("segment", x = 11552000, xend = 11560000, y = 4335000, yend = 4335000, colour = "black", size = 0.2, arrow = arrow(length = unit(0.1, "cm"))) +
  labs(x = "", y = "") +
  scale_x_continuous(limits = c(11488000, 11560000)) +
  coord_fixed()

# plot all early fires
p_all_early <- ggplot() + 
  geom_polygon(data = no_fire, aes(x = x, y = y, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 12000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 24000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 36000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = all_late_fires, aes(x = x, y = y + 36000, group = id), fill = "black", color = NA) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 48000, group = id), fill = NA, color = "black", size = 0.1) +
  geom_polygon(data = no_fire, aes(x = x, y = y + 60000, group = id), fill = NA, color = "black", size = 0.1) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.margin=unit(c(0, -0.8, 0, -1.0), "cm")) +
  annotate("text", x = 11530000, y = 4311000, size = 2, color = "black", hjust = 0, label = "November/December") +
  annotate("text", x = 11530000, y = 4323000, size = 2, color = "black", hjust = 0, label = "September/October") +
  annotate("text", x = 11530000, y = 4335000, size = 2, color = "black", hjust = 0, label = "July/August") +
  annotate("text", x = 11530000, y = 4347000, size = 2, color = "black", hjust = 0, label = "May/June") +
  annotate("text", x = 11530000, y = 4359000, size = 2, color = "black", hjust = 0, label = "March/April") +
  annotate("text", x = 11530000, y = 4371000, size = 2, color = "black", hjust = 0, label = "January/February") +
  annotate("text", x = 11525000, y = 4298000, size = 3, color = "black", hjust = 0, label = "All Early Fires") +
  labs(x = "", y = "") +
  scale_x_continuous(limits = c(11488000, 11560000)) +
  coord_fixed()

grid.arrange(p_early_late, p_all_early, ncol = 2, padding = 0)

dev.off()


# Create EMA plot & trend plots
load("possum_pop_data.RData")
load("possum_AF_pop_data.RData")
load("possum_AEF_pop_data.RData")
load("possum_AF_AEF_pop_data.RData")

load("melomys_pop_data.RData")
load("melomys_AF_pop_data.RData")
load("melomys_AEF_pop_data.RData")
load("melomys_AF_AEF_pop_data.RData")

load("bandicoot_pop_data.RData")
load("bandicoot_AF_pop_data.RData")
load("bandicoot_AEF_pop_data.RData")
load("bandicoot_AF_AEF_pop_data.RData")

pop_data_list <- list("Common brushtail possum" = list("Observed fire history - unadjusted intensity" = possum_pop_data,
                                                       "Observed fire history - adjusted intensity" = possum_AF_pop_data,
                                                       "EDS prescribed burning scenario - unadjusted intensity" = possum_AEF_pop_data,
                                                       "EDS prescribed burning scenario - adjusted intensity" = possum_AF_AEF_pop_data),
                      "Grassland melomys" = list("Observed fire history - unadjusted intensity" = melomys_pop_data,
                                                 "Observed fire history - adjusted intensity" = melomys_AF_pop_data,
                                                 "EDS prescribed burning scenario - unadjusted intensity" = melomys_AEF_pop_data,
                                                 "EDS prescribed burning scenario - adjusted intensity" = melomys_AF_AEF_pop_data),
                      "Northern brown bandicoot" = list("Observed fire history - unadjusted intensity" = bandicoot_pop_data,
                                                        "Observed fire history - adjusted intensity" = bandicoot_AF_pop_data,
                                                        "EDS prescribed burning scenario - unadjusted intensity" = bandicoot_AEF_pop_data,
                                                        "EDS prescribed burning scenario - adjusted intensity" = bandicoot_AF_AEF_pop_data))


png('figs/ema_comparison.png',
    pointsize = 6,
    res = 300,
    width = 2400,
    height = 1200)

plot_ema(pop_data_list,
         all_points = FALSE,
         p = 0.01)

dev.off()


png('figs/possum_pop_trends.png',
    pointsize = 8,
    res = 300,
    width = 1600,
    height = 1600)

pop_data_species <- pop_data_list[["Common brushtail possum"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(2, 2))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()


png('figs/melomys_pop_trends.png',
    pointsize = 8,
    res = 300,
    width = 1600,
    height = 1600)

pop_data_species <- pop_data_list[["Grassland melomys"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(2, 2))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()


png('figs/bandicoot_pop_trends.png',
    pointsize = 8,
    res = 300,
    width = 1600,
    height = 1600)

pop_data_species <- pop_data_list[["Northern brown bandicoot"]]
max_pop <- ceiling(max(vapply(pop_data_species, max, 1)))

par(mfrow = c(2, 2))

for (i in 1:length(pop_data_species)) {
  pop_data_scenario <- pop_data_species[[i]]
  plot_pop(pop_data = pop_data_scenario,
           upper_limit = max_pop,
           title = names(pop_data_species)[i])
}

dev.off()

pop_data <- pop_data_list[["Northern brown bandicoot"]]
lapply(pop_data, function(x) round(mean(apply(x, 1, function(y) min(y))), 0))

####### New table to compare percentage changes in final mean estimates of population

perc_change <- foreach(i = 1:length(pop_data_list), .combine = cbind) %do% {
  
  mean_last_pops <- unlist(lapply(pop_data_list[[i]], function(x) mean(apply(x, 1, function(y) tail(y, 1)))))
  outer(1:4, 1:4, FUN = function(x, y) (mean_last_pops[x] - mean_last_pops[y]) / mean_last_pops[x] * 100)
}

htmlTable(
  x = round(perc_change, 2),
  caption = paste("Table X. Percentage change of mean final populations",
                   "between scenarios for each species. Table is read as",
                   "change from column to row. Positive values indicate",
                   "increases and negative values indicate decreases."),
  cgroup = names(pop_data_list),
  n.cgroup=c(4, 4, 4),
  rnames = 1:4)

tableHTML(round(perc_change, 2),
          second_header = list(c(1, 4, 4, 4), c("", names(pop_data_list))),
          )

          