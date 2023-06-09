# Figure 1 - plot map of Mimulus guttatus collection locations -----------------

# import data
setwd("~/Desktop/DRYAD/")
dat <- read.csv(file = "01_Rotter_Mimgut_locations.csv", header = T, stringsAsFactors = F )

# plot study locations on a world map using ggplot
library(maps)
library(ggmap)
library(ggplot2)
library(ggrepel)

# import world map
world <- map_data("world")


# plot study site locations ----------------------------------------------------
dev.new(width = 8.5, height = 5, noRStudioGD = TRUE)

mp <- ggplot(dat, aes(x = lon, y = lat, fill = Subregion))

mp + 
    geom_point(shape = 21, size = 2.5) +
    xlim(-165,10) +
    ylim(27.5, 70) +
    labs(x = "Longitude", y = "Latitude", font = "bold") +
    theme(axis.title = element_text(size = 14, face = "bold")) +
    theme(legend.text = element_text(size = 14)) +
    theme(legend.title = element_text(size = 14, face = "bold")) +
    geom_text_repel(aes(label = Population, fontface = "bold"), max.overlaps = 25) +
    scale_fill_discrete(name = "Biogeographic Region") +
    theme(legend.position = "top") +
    
    geom_map(
        data = world, map = world,
        aes(long, lat, map_id = region),
        color = "black", fill = "gray60", size = 0.1, alpha = 0.1) 

# quartz.save(file = "Figure1_population_map.jpg", type = "jpg", dpi = 300)

