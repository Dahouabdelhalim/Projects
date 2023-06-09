# A global map of the fraction of derived woody species on islands

# libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(sp)
library(viridis)
library(cowplot)

# load data
load(file = "output/final_island_dataset_plots.rda")

polys <- st_read(dsn = "input", layer = "geoentities_simple")# spatial data

# prepare for plotting
dat <- iw_plot %>% 
  dplyr::select(geo_entity, 
         entity_ID, 
         entity_class,
         phylowood_iwspnum,
         GIFT_list_spnum,
         GIFT_list_spnum_no_monocots,
         geology_simple,
         phylowood_dwspnum) %>% 
  mutate(prop_iw = phylowood_iwspnum / GIFT_list_spnum * 100)%>% 
  mutate(prop_iw_nomonocots = phylowood_iwspnum / GIFT_list_spnum_no_monocots * 100)

## Island polygons
st_geometry(polys) <-  NULL
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
wgs1984 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

bg <-  ne_countries(scale = 110, returnclass = "sf")%>% st_transform(rob)

plo<- polys %>% 
  dplyr::select(entt_ID, point_x, point_y, area, entty_c) %>% 
  right_join(dat, by = c("entt_ID" = "entity_ID")) %>% 
  mutate(point_x = ifelse(is.na(point_x), longitude, point_x))%>% 
  mutate(point_y = ifelse(is.na(point_y), latitude, point_y))%>% 
  # dplyr::select(-point_x, -point_y) %>% 
  mutate(species = log(phylowood_iwspnum)) %>% 
  mutate(species = ifelse(species < 0, -1, species)) %>% 
  arrange(species)

# re-project
repro <- SpatialPoints(plo[, c("point_x", "point_y")], proj4string = wgs1984) %>%
  spTransform(rob) %>%
  coordinates() %>% 
  as_tibble()

plo <- plo %>% 
  dplyr::select(-point_x, -point_y) %>% 
  bind_cols(repro)

# The plot
ggplot()+
  geom_sf(data = bg, fill = "grey90", color = "grey90")+
  geom_point(data = plo, aes(x = point_x, y = point_y, color = species, shape = geology_simple), 
             size = 4 )+
  scale_color_viridis(name = "Insular woody species",
                      breaks = log(c(1, 10, 50, 90)),
                      labels = c(1, 10, 50, 90))+
  scale_shape_manual(values = c(17,16), name = "")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(2, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18))

ggsave("output/figure_number_of_IW_species_per_island.pdf")

# Supplementary plot fraction
plo <- plo %>% arrange(prop_iw)

p1 <- ggplot()+
  geom_sf(data = bg, fill = "grey90", color = "grey90")+
  geom_point(data = plo, aes(x = point_x, y = point_y, color = prop_iw, shape = geology_simple), 
             size = 3 )+
  scale_color_viridis(name = "Proportion of\\nInsular woody\\nspecies")+
  scale_shape_manual(values = c(17,16), name = "")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# Supplementary plot insular woodiness
plo <- plo %>% arrange(prop_iw_nomonocots)

p2 <- ggplot()+
  geom_sf(data = bg, fill = "grey90", color = "grey90")+
  geom_point(data = plo, aes(x = point_x, y = point_y, color = prop_iw_nomonocots, shape = geology_simple), 
             size = 3 )+
  scale_color_viridis(name = "Proportion of \\nInsular woody\\nof non-monocot species")+
  scale_shape_manual(values = c(17,16), name = "")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

props <- plot_grid(p1, p2, labels = c('A', 'B'), ncol = 1)

ggsave("output/supplementary_figure_proportions_IW.jpg", height = 10, width = 8)


# Supplementary plot derived woody not insular woody
plo <- plo %>% 
  mutate(dw_ex = phylowood_dwspnum - phylowood_iwspnum) %>% 
  arrange(dw_ex)

ggplot()+
  geom_sf(data = bg, fill = "grey90", color = "grey90")+
  geom_point(data = plo, aes(x = point_x, y = point_y, color = log(dw_ex +1), shape = geology_simple), 
             size = 4 )+
  scale_color_viridis(name = "Number of\\nDW not IW\\nspecies [log]",
                      breaks = c(log(1), log(10), log(50), log(100)),
                      labels = c(1, 10, 50, 100))+
  scale_shape_manual(values = c(17,16), name = "")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave("output/supplementary_figure_number_of_DW_species_not_IW.jpg", height = 5, width = 8)
