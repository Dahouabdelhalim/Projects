library(phylowood) 
library(tidyverse)
library(viridis)
library(readxl)
library(sf)
library(rnaturalearth)
library(writexl)

# get the numbers from phylowood_sp and phylowood_isl
load(file = "output/GIFT_archipelago_ID.rda")
load("output/GIFT_island_characteristica.rda")

iw_species <- phylowood_sp %>% 
  filter(geo_insularwoody) %>%
  dplyr::select(geo_insularwoody, geo_archipelago) %>% 
  separate_rows(geo_archipelago, sep = ",") %>% 
  group_by(geo_archipelago) %>% 
  count()

shifts <- phylowood_isl %>% 
  group_by(geo_archipelago) %>% 
  dplyr::summarize(shifts = sum(trt_minshift_nr)) %>% 
  filter(shifts > 0) %>% 
  arrange(desc(shifts)) %>% 
  left_join(iw_species) %>% 
  mutate(geo_archipelago = recode(geo_archipelago, 
                                 `Galapagos Islands` = "GalÃ¡pagos Islands",
                                 `Mozambique Channel Islands` = "Mozambique Channel Islands: Bassas da India, Europa, Glorioso Islands, Juan de Nova")) %>% #also in the archipelago file
  mutate(ID = substr(geo_archipelago, start = 1, stop = 3) %>% toupper())

write_xlsx(shifts %>%  dplyr::select(geo_archipelago, ID), path = "output/archipelagos_roberto.xlsx")

# get the coordinates for the archipelagos from GIFT
polys <- st_read(dsn = "input", layer = "geoentities_simple") %>% 
  dplyr::select(entt_ID, ge_ntty, point_x, point_y, area)

archp <- geoentities_env_misc%>% 
  left_join(arch_gift %>% filter(Arch_level == 1) %>% dplyr::select(Arch_ID, Arch_name), by = c("arch_lvl_1" = "Arch_ID")) %>% 
  dplyr::select(entity_ID, Arch_name, latitude, longitude) %>% 
  st_as_sf(coords = c("longitude","latitude")) %>% 
  right_join(shifts, by = c("Arch_name" = "geo_archipelago"))

circs <- archp%>% 
  group_by(Arch_name) %>% 
  summarise() %>% 
  st_convex_hull() %>% 
  st_centroid() %>% 
  left_join(shifts, by = c("Arch_name" = "geo_archipelago"))

add <- circs %>% 
  filter(sf::st_is_empty(.))

st_geometry(add) <- NULL

add <- add %>% 
  filter(ID != "AMB") %>% 
  mutate(X = c(-79.8880, 58.31)) %>% 
  mutate(Y = c(-26.34, -48.4)) %>% 
  st_as_sf(coords = 5:6)

circs <- circs %>% 
  filter(!st_is_empty(.)) %>% 
  bind_rows(add)%>%
  mutate(n= ifelse(is.na(n), 0, n))

shifts$geo_archipelago[!shifts$geo_archipelago %in% archp$Arch_name]

sort(unique(archp$Arch_name))

# Plot map
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
wgs1984 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

bg <-  ne_countries(scale = 110, returnclass = "sf")%>% st_transform(rob)

st_crs(circs) <-  wgs1984
circs <- circs  %>%  st_transform(rob) %>% 
  mutate(labl = sprintf("%s: %s/%s", ID, n, shifts)) %>% 
  bind_cols(data.frame(st_coordinates(.)))


ggplot()+
  geom_sf(data = bg, fill = "grey80", color = "grey80", size = 0.1)+
  geom_point(data = circs, aes(x = X, y = Y, fill = shifts, size = n), alpha=0.5, shape=21, color="black") +
  scale_size(range = c(2, 25), name="Number of insular\\nwoody species") +
  scale_fill_viridis(discrete=FALSE, name = "Evolutionary\\nshifts to\\nwoodiness", alpha = 0.5) +
  geom_text(data = data.frame(x = 5153441.0, y = -5060865.043), aes(x = x, y = y), label = "*", size = 5)+
  #geom_sf(data = convxh, aes(color = Arch_name), fill = "transparent")+ 
  #geom_sf_label(data = convxh, aes(label = labl, color = Arch_name), size = 4)+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1, "cm"),
        legend.key.height = unit(1, "cm"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16))

ggsave("output/figure_map_shifts_per_archipelago.pdf")
