library(tidyverse)
library(car)
library(olsrr)
library(rnaturalearth)
library(sp)
library(sf)

load(file = "output/final_island_dataset_models.rda")

#  pairs plot
pairs(iw_model[,13:26])

# test correlation with vartiance inflation factor
model <- lm(phylowood_iwspnum ~ area + 
              dist + 
              mam_total +
              mean_CHELSA_NFD +
              mean_ai_yr + 
              mean_b12_v0 +
              mean_b1_v0 +
              mean_elev +
              mean_wc2.0_bio_30s_01 + 
              mean_wc2.0_bio_30s_04 + 
              mean_wc2.0_bio_30s_12 + 
              mean_wc2.0_bio_30s_15 +
              mean_wc2.0_bio_30s_17 +
              mean_wc2.0_bio_30s_18,
            data = iw_model)

ols_vif_tol(model)

# remove MAP and MAT and precipitation of the driest quarter
model <- lm(phylowood_iwspnum ~ area + 
              dist + 
              mam_total +
              mean_CHELSA_NFD +
              mean_ai_yr + 
              mean_b12_v0 +
              mean_b1_v0 +
              mean_elev +
              mean_wc2.0_bio_30s_04 + 
              mean_wc2.0_bio_30s_15 +
              mean_wc2.0_bio_30s_18,
            data = iw_model)

ols_vif_tol(model)
pairs(iw_model[,c(13:20,22,24,26)])

# write to disk for SEM/SAR modelling
iw_model <- iw_model %>% 
  dplyr::select(geo_entity, 
                Arch_name,
                GIFT_list_spnum,
                phylowood_iwspnum,
                geology_simple,
                latitude:mean_elev,
                mean_wc2.0_bio_30s_04,
                mean_wc2.0_bio_30s_15,
                mean_wc2.0_bio_30s_18)

save(iw_model, file = "output/data_for_island_model_final.rda")
write_csv(iw_model, file = "output/data_for_island_model_final.csv")

# plot number of IW species against all variables
plo <- iw_model %>% 
  mutate(fraction_iw = phylowood_iwspnum / GIFT_list_spnum * 100) %>% 
  dplyr::select(fraction_iw, geology_simple, area:mean_wc2.0_bio_30s_18 ) %>% 
    rename("Island area" = area,
           "Distance to the nearest continent" = dist,
           "Time-integrated species richness of large mammal herbivores" = mam_total,
           "Mean aridity index" = mean_ai_yr,
           "Past climate change velocity in temperature" = mean_b1_v0,
           "Past climate change velocity in precipitation" = mean_b12_v0,
           "Number of frost days" = mean_CHELSA_NFD,
           "Mean elevation" = mean_elev,
           "Temperature Seasonality" = mean_wc2.0_bio_30s_04,
           "Precipitation Seasonality" = mean_wc2.0_bio_30s_15,
           "Precipitation of Warmest Quarter" = mean_wc2.0_bio_30s_18) %>% 
  pivot_longer(cols = `Island area`:`Precipitation of Warmest Quarter`)

ggplot(data = plo, aes(x= value, y = fraction_iw, colour = geology_simple))+
  geom_point()+
  #geom_smooth(method = "glm")+
  facet_wrap(name~., ncol =  3, scales = "free_x")+
  ylab("")+
  scale_color_discrete(name = "Archipelago type")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

ggsave("output/supplementary_figure_island_level_data.pdf", height = 10, width = 8)

plo <- iw_model %>% 
  mutate(fraction_iw = phylowood_iwspnum / GIFT_list_spnum * 100) %>% 
  dplyr::select(fraction_iw, geology_simple, area:mean_wc2.0_bio_30s_18) %>% 
  pivot_longer(cols = area:mean_wc2.0_bio_30s_18) %>% 
  filter(geology_simple == "oceanic")

ggplot(data = plo, aes(x= value, y = fraction_iw))+
  geom_point()+
  #geom_smooth(method = "glm")+
  facet_wrap(name~., ncol =  3, scales = "free_x")+
  ylab("")+
  scale_color_discrete(name = "Archipelago type")+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title.x = element_blank())

ggsave("output/supplementary_figure_island_level_data_oceanic.jpeg", height = 10, width = 8)


# plot the islands included in the analysis
rob <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
wgs1984 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

bg <-  ne_countries(scale = 110, returnclass = "sf")%>% st_transform(rob)

trop <- bind_rows(data.frame(longitude = seq(-180, 180, 1),
                             latitude = 23.44, 
                             type = "cancer"),
                  data.frame(longitude = seq(-180, 180, 1),
                             latitude = -23.44, 
                             type = "capricorn"))


trop <- SpatialPoints(trop[, c("longitude", "latitude")], proj4string = wgs1984) %>%
  spTransform(rob) %>%
  coordinates() %>% 
  as_tibble()%>% 
  bind_cols(trop %>% dplyr::select(-longitude, -latitude)) 

repro <- SpatialPoints(iw_model[, c("longitude", "latitude")], proj4string = wgs1984) %>%
  spTransform(rob) %>%
  coordinates() %>% 
  as_tibble()

plo <- iw_model %>% 
  dplyr::select(-longitude, -latitude) %>% 
  bind_cols(repro) 

ggplot()+
  geom_sf(data = bg, fill = "grey90", color = "grey90")+
  geom_line(data = trop, aes(x = longitude, y = latitude, group = type),
            color = "grey50", alpha = 0.5)+
  geom_point(data = plo, aes(x = longitude,
                             y = latitude, 
                             color = geology_simple, shape = geology_simple), 
             size = 4)+
  scale_color_manual(values = c("black", "red"))+
  theme_bw()+
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(1.5, "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_blank()) #+ 
 # guides(color=FALSE)

ggsave("output/supplementary_figure_islands_included_island_analysis.pdf")
