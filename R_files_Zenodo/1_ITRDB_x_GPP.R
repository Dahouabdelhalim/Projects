##### 
##### Identify RW sites that are nearby flux tower sites and overlap GPP time series
#####

library(sf)
library(geosphere)

### flux sites  
flux_species <- read_csv("Data/flux_species.csv")
flux_metadata <- read_csv("Data/flux_metadata.csv") 
flux_points <- flux_metadata %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = "+proj=longlat +datum=WGS84")

### RW sites 
rw_metadata <- read_csv("Data/rw_metadata.csv")
rw_points <- rw_metadata %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = "+proj=longlat +datum=WGS84")

rw_onsite_metadata <- read_csv("Data/rw_onsite_metadata.csv")
rw_onsite_sites <- unique(rw_onsite_metadata$Site)

### Distance RW-Flux towers
dist <- distm(flux_points %>% as_Spatial,
              rw_points %>% as_Spatial) # res in m

dist_x <- function(x = 1:1000, dist){
  res <- c()
  for(i in x){
    dist_i <- map(1:nrow(dist), function(y) which(dist[y,]<i*1e3)) # in km
    res_i <- sum(map(dist_i, function(y) length(y) != 0) %>% unlist)
    res <- c(res, res_i)
  }
  return(res)
}

x <- seq(0,5000,by = 50)
dist_0_1000 <- dist_x(x, dist)
plot(x, dist_0_1000, type = "l")
points(1000, dist_0_1000[x == 1000], col = "red")
points(500, dist_0_1000[x == 500], col = "red")
points(200, dist_0_1000[x == 200], col = "red")
points(100, dist_0_1000[x == 100], col = "red")
points(50, dist_0_1000[x == 50], col = "red")

## Create table containing flux sites and names and distance of corresponding closest RW sites
dist_obj <- map(1:nrow(dist), function(x) which(dist[x,]<1000*1e3))

stdMAT <- function(x) (x-mean(rw_metadata$MAT, na.rm = T))/sd(rw_metadata$MAT, na.rm = T)
stdMAP <- function(x) (log(x)-mean(log(rw_metadata$MAP), na.rm = T))/sd(log(rw_metadata$MAP), na.rm = T)

flux_rw_overlap <- flux_metadata %>% 
  ungroup %>%
  mutate(dist_obj = dist_obj,
         RW_site = map(dist_obj, function(x) rw_metadata$Site[x]),
         RW_geodist = map2(seq_along(Site), dist_obj, function(x,y) dist[x,unlist(y)]/1e3)) %>%
  unnest(-dist_obj) %>% select(-dist_obj) %>%
  filter(!is.na(RW_site)) %>%
  mutate(RW_geodist = round(RW_geodist,1)) %>%
  left_join(flux_species) %>% 
  filter(!is.na(Group)) %>% ## NA-filled group correspond to those flux sites that do not have sufficient overlap with RW data or species classified as "other"
  left_join(rw_metadata %>% select(RW_site = Site, RW_species = Species, RW_genera = Genera, RW_family = Family, RW_group = Group, RW_evergreen = Evergreen, 
                                   RW_start = Start, RW_end = End, RW_MAT = MAT, RW_MAP = MAP, RW_elev = Elev, RW_source = Source)) %>%
  mutate(RW_climdist = sqrt(
    (stdMAP(MAP)-stdMAP(RW_MAP))^2+
      (stdMAT(MAT)-stdMAT(RW_MAT))^2
  )) %>% 
  group_by(Site, Species, RW_site) %>% 
  mutate(RW_overlap = sum(Start:End %in% RW_start:RW_end)) %>%
  filter(RW_overlap >= 4) %>%
  mutate(RW_phylodist = (RW_species != Species) + (RW_genera != Genera) + 
                            (RW_family != Family) + (RW_group != Group) +
                            (RW_evergreen != Evergreen)) %>% 
  group_by(Site,RW_site) %>% 
  summarise(RW_phylodist = weighted.mean(RW_phylodist, Proportion),
            across(c(Country:Source, contains("RW")), unique)) %>% 
  filter(!is.na(RW_climdist))

write_csv(flux_rw_overlap %>% select(Site, RW_site, RW_phylodist, RW_geodist, RW_climdist, RW_overlap), "Data/flux_rw_overlap.csv")

map <- ggplot()+
  geom_sf(data = world)+
  geom_point(data = rw_metadata %>%
               filter(Site %in% flux_rw_overlap$RW_site),
             aes(x = Lon, y = Lat, shape = "1"),
             color = "black", size = 2, alpha = 0.6)+
  # geom_point(data = flux_rw_overlap %>% 
  #              filter(!Site %in% rw_onsite_sites) %>% 
  #              group_by(Lon, Lat, End, Start) %>% summarise(RW_overlap = sum(RW_overlap)),
  #            aes(x = Lon, y = Lat, size = RW_overlap),
  #            color = "red", alpha = 0.8, show.legend = F)+
  geom_point(data = flux_rw_overlap %>%
            # filter(Site %in% rw_onsite_sites) %>% 
            group_by(Lon, Lat, End, Start, Site) %>% summarise(RW_overlap = sum(RW_overlap)), 
            aes(x = Lon, y = Lat, size = RW_overlap, fill = Site %in% rw_onsite_sites, shape = "2"),
          color = "black", alpha = 0.8)+
  scale_size(range = c(1,5), guide = F)+
  scale_fill_manual(name = "Onsite RW data", values = c("red", "yellow"), labels = c("Has onsite RW", "No onsite RW"), guide = F)+
  scale_shape_manual(name = NULL, values = c(3,21), guide = F)+
  coord_sf(expand = F, ylim = c(-60, 80), xlim = c(-165, 165))+
  xlab("Longitude")+ylab("Latitude")

map_plot <- map+
  coord_sf(expand = F, ylim = c(10, 72), xlim = c(-155, 45))+
  annotation_custom(ggplotGrob(map+theme(text = element_blank(), axis.ticks = element_blank(), plot.margin = NULL)), 
                    xmin = -32, ymax = 35)+
  scale_shape_manual(name = "Site", values = c(3,21), labels = c(expression(RW[region], "GPP")))+
  theme(legend.position = c(0.01,0.01), legend.justification = c(0,0), legend.background = element_blank())
save(map_plot, file = "Graphs/Final/Figure 1A.RData")

group_plot <- flux_rw_overlap %>% 
  group_by(x = RW_group) %>% 
  summarise(n = sum(RW_overlap)) %>% 
  ungroup %>% 
  arrange(n) %>% 
  mutate(x = factor(x, levels = unique(x))) %>% 
  ggplot(aes(x = x, y = n))+
  geom_col(fill = "grey")+
  xlab(NULL)+ylab(expression(GPP-RW[region]~site~x~yr))+
  coord_cartesian(expand = T, ylim = c(470*8, 80000))+
  theme(axis.text.x = element_text(angle = 25, hjust = 0.9))
save(group_plot, file = "Graphs/Final/Figure 1B.RData")
