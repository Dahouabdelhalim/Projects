# zoom in plots for four archipelagos
# A global map of the fraction of derived woody species on islands
# libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(sp)
library(viridis)

# load data
load(file = "output/final_island_dataset_plots.rda")
polys <- st_read(dsn = "input", layer = "geoentities_simple")# spatial data

iw_plot <- iw_plot %>% 
  mutate(prop_iw = phylowood_iwspnum / GIFT_list_spnum * 100)%>% 
  mutate(prop_iw_nomonocots = phylowood_iwspnum / GIFT_list_spnum_no_monocots * 100)

#Canaries
dat <- filter(iw_plot, Arch_name == "Canary Islands")
plo_ma <- polys %>% 
  right_join(dat %>% dplyr::select(entity_ID, geo_entity, prop_iw_nomonocots), 
             by = c( "entt_ID" = "entity_ID")) %>% 
  mutate(prop_iw_nomonocots = cut(x = as.numeric(prop_iw_nomonocots ), 
                                  breaks = c(0, 5,10,15,20, 1000),
                                  labels = c("0-5", "5-10", "10-15", "15-20", ">20")))%>% 
  mutate(prop_iw_nomonocots = factor(prop_iw_nomonocots,
                                     levels = c("0-5", "5-10", "10-15", "15-20", ">20")))

ggplot()+
  geom_sf(data = plo_ma, aes(fill = prop_iw_nomonocots), size = 0.01)+
  scale_fill_viridis(name = "Fraction of\\nderived\\nwoody species\\n[%]",
                     discrete = TRUE, 
                     drop = FALSE, 
                     option = "B")+
  ggtitle("Canary Islands")+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = c(1,0),
        legend.justification= c(1,0),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 28))

ggsave(file = "output/figure_canaries_proportion_IW.pdf")

# Hawaii
dat <- filter(iw_plot , Arch_name == "Hawaiian Islands")
plo_ma <- polys %>% 
  right_join(dat %>% dplyr::select(entity_ID, geo_entity,  prop_iw_nomonocots), 
             by = c( "entt_ID" = "entity_ID")) %>% 
  mutate(prop_iw_nomonocots = cut(x = as.numeric(prop_iw_nomonocots ), 
                                  breaks = c(0, 5,10,15,20, 1000),
                                  labels = c("0-5", "5-10", "10-15", "15-20", ">20")))%>% 
  mutate(prop_iw_nomonocots = factor(prop_iw_nomonocots,
                                     levels = c("0-5", "5-10", "10-15", "15-20", ">20")))


ggplot()+
  geom_sf(data = plo_ma, aes(fill =  prop_iw_nomonocots), size = 0.01)+
  scale_fill_viridis(name = "Fraction of\\nderived\\nwoody species\\n[%]",
                     discrete = TRUE, 
                     drop = FALSE,
                     option = "B")+
  xlim(-160.5, -154.5)+
  ylim(18.5, 22.5)+
  ggtitle("Hawaii (main islands)")+
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0,0),
        legend.justification= c(-0.2,-0.1),
        legend.key.height = unit(1.2, "cm"),
        legend.key.width = unit(1.2, "cm"),
        legend.text = element_text(size = 20),
        plot.title = element_text(size = 28))

ggsave(file = "output/figure_hawaii_fraction_dw_species.pdf")

# #New Zealand
# dat <- filter(iw_plot, Arch_name == "New Zealand") 
# 
# plo_ma <- polys %>% 
#   right_join(dat %>% dplyr::select(entity_ID, geo_entity, prop_iw_nomonocots), 
#              by = c( "entt_ID" = "entity_ID"))  %>% 
#   mutate(prop_iw_nomonocots = cut(x = as.numeric(prop_iw_nomonocots ), 
#                                   breaks = c(-1, 5,10,15,20, 1000),
#                                   labels = c("0-5", "5-10", "10-15", "15-20", ">20")))%>% 
#   mutate(prop_iw_nomonocots = factor(prop_iw_nomonocots,
#                                      levels = c("0-5", "5-10", "10-15", "15-20", ">20"))) 
# 
# 
# ggplot()+
#   geom_sf(data = plo_ma, aes(fill = prop_iw_nomonocots), size = 0.01)+
#   scale_fill_viridis(name = "Fraction of\\nderived\\nwoody species\\n[%]",
#                      discrete = TRUE, 
#                      drop = FALSE)+
#   ggtitle("New Zealand")+
#   theme_bw()+
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title = element_blank(),
#         legend.position = c(1,0),
#         legend.justification= c(1.1,-0.1),
#         legend.key.height = unit(1.2, "cm"),
#         legend.key.width = unit(1.2, "cm"),
#         legend.text = element_text(size = 20),
#         plot.title = element_text(size = 28))
# 
# ggsave(file = "output/figure_new_zealand_fraction_dw_species.pdf")

# # Mascarenes
# dat <- filter(iw_plot, Arch_name  %in% c("Mascarene Islands", "Rodrigues Island"))
# plo_ma <- polys %>% 
#   right_join(dat %>% dplyr::select(entity_ID, geo_entity, prop_iw_nomonocots), 
#              by = c( "entt_ID" = "entity_ID")) %>% 
#   mutate(prop_iw_nomonocots = cut(x = as.numeric(prop_iw_nomonocots ), 
#                                   breaks = c(0, 5,10,15,20, 1000),
#                                   labels = c("0-5", "5-10", "10-15", "15-20", ">20")))%>% 
#   mutate(prop_iw_nomonocots = factor(prop_iw_nomonocots,
#                                      levels = c("0-5", "5-10", "10-15", "15-20", ">20")))
# 
# 
# ggplot()+
#   geom_sf(data = plo_ma, aes(fill = prop_iw_nomonocots), size = 0.01)+
#   scale_fill_viridis(name = "Fraction of\\nderived\\nwoody species\\n[%]",
#                      discrete = TRUE, 
#                      drop = FALSE)+
#   ggtitle("Mascarenes")+
#   theme_bw()+
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title = element_blank(),
#         legend.key.height = unit(1.2, "cm"),
#         legend.key.width = unit(1.2, "cm"),
#         legend.text = element_text(size = 20),
#         plot.title = element_text(size = 28),
#         legend.position = "bottom")
# 
# ggsave(file = "output/figure_mascarenes_fraction_dw_species.pdf")


# phylowood_sp %>% filter(geo_archipelago == "Mascarene Islands") %>%  filter(geo_insularwoody) %>% pull(tax_genus) %>% unique()

