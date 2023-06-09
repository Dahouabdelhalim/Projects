#####
##### Map Flux and RW sites in climatic space
#####

library(plotbiomes)

### Metadata 
flux_rw_overlap <- read_csv("Data/flux_rw_overlap.csv")
flux_metadata <- read_csv("Data/flux_metadata.csv") %>% 
  filter(Site %in% flux_rw_overlap$Site) %>% 
  mutate(temp_c = MAT, )
rw_metadata <- read_csv("Data/rw_metadata.csv") %>% filter(Site %in% flux_rw_overlap$RW_site)
rw_onsite_metadata <- read_csv("Data/rw_onsite_metadata.csv") 
rw_onsite_sites <- unique(rw_onsite_metadata$Site)

### Plot
whit_plot <- ggplot() + 
  geom_polygon(data = Whittaker_biomes %>% mutate(MAT = temp_c, MAP = precp_cm),
               aes(x = MAT, y = MAP, fill = biome),
               colour = "gray98", size = 1, show.legend = F) +
  stat_density_2d(data = rw_metadata, aes(x = MAT, y = MAP/10), color = "black", alpha = 1, bins = 10, h = c(9,90), size = 0.5)+
  geom_point(data = flux_rw_overlap %>% group_by(Site) %>% summarise(Span=sum(RW_overlap)) %>% left_join(flux_metadata %>% select(Site, MAT, MAP)),
             aes(x = MAT, y = MAP/10, size = Span, color = Site %in% rw_onsite_sites), 
             shape = 16, alpha = 0.8, show.legend = T)+
  geom_point(data = flux_rw_overlap %>% group_by(Site) %>% summarise(Span=sum(RW_overlap)) %>% left_join(flux_metadata %>% select(Site, MAT, MAP)),
             aes(x = MAT, y = MAP/10, size = Span), 
             shape = 1, alpha = 0.8, color = "black",show.legend = T)+
  # scale_fill_brewer(type = "qual", palette = 3)
  # geom_abline(intercept = 150, slope = 10)
  coord_cartesian(expand = F)+
  scale_size(name = expression(atop(Paired~GPP-RW[region], obs.~(site~x~year))), breaks = c(1000,3000,6000))+
  scale_color_manual(name = expression(RW[on-site]~available), values = c("red", "yellow"), labels = c("No", "Yes"))+
  scale_y_continuous(name = "MAP (mm)", labels = function(x) as.numeric(x)*10)+
  scale_x_continuous(name = "MAT (°C)")+
  scale_fill_grey(guide = F, start = 0.4, end = 0.8)+
  theme(legend.position = c(0.01,0.99), legend.justification = c(0,1), legend.box = "horizontal", legend.background = element_blank())
whit_plot
save(whit_plot, file = "Graphs/Final/Figure 1C.RData")
# ggsave("Graphs/Whittaker bioclimatic space.pdf", width = 6, height = 6)


whit_plot_template <- ggplot() + 
  geom_polygon(data = Whittaker_biomes,
               aes(x = temp_c, y = precp_cm, fill = biome),
               colour = "gray98", size = 1) +
  scale_fill_manual(name = "Whittaker biomes",
                    breaks = names(Ricklefs_colors), labels = names(Ricklefs_colors),
                    values = Ricklefs_colors)+
  scale_y_continuous(name = "MAP (mm)", labels = function(x) as.numeric(x)*10)+
  scale_x_continuous(name = "MAT (°C)")+
  theme(legend.position = c(0.01,0.99), legend.justification = c(0,1), legend.box = "horizontal", legend.background = element_blank())
whit_plot_template
ggsave("Graphs/Whittaker template.pdf", width = 6, height = 6)
