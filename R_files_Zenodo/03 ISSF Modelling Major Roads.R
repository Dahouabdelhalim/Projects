#### TITLE ####
# How do King Cobras move across a major highway? Unintentional
# wildlife crossing structures may facilitate movement.
#### AUTHORS ####
# Max Dolton Jones, Benjamin Michael Marshall, Samantha Nicole Smith, 
# Matt Crane, InÃªs Silva, Taksin Artchawakom, Pongthep Suwanwaree, 
# Surachit Waengsothorn, Wolfgang Wüster, Matt Goode, Colin Thomas Strine
#### YEAR ####
# 2021

## 03 Individual ISSF SCRIPT

# Libraries ---------------------------------------------------------------

library(rgdal)
library(dplyr)
library(raster)
library(readr)
library(stringr)
library(amt)
library(broom)
library(ggplot2)
library(here)
library(scico)

# Create Folders ----------------------------------------------------------

loc.data <- paste0("./Data/")
loc.fig <- paste0("./Figures/")
loc.shape <- paste0("./Shapefiles")

# Read in Data and Resources ----------------------------------------------
landuse <- readOGR(paste0(loc.shape, "/Landuse_2018-12-11.shp"))

list.files(loc.data)

core <- readOGR(paste0(loc.shape, "/SBRcore.shp")) %>%
  crop(extent(landuse))

data.var <- read_csv(file = paste0(loc.data, "Ophiophagus hannah 2013_2020 ISSF_adults.csv"),
                     locale = locale(tz = "Asia/Bangkok"))


crs.proj <- CRS("+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Function to fit iSSF model
fitted_ssf <- function(data, model){
  mod <- fit_issf(model, data = data)
  return(mod)
}

# Create Landuse Rasters --------------------------------

(extent(landuse)[1] - extent(landuse)[2]) / 10
(extent(landuse)[3] - extent(landuse)[4]) / 10

###### major roads #########

mj.roads <- readOGR(dsn = paste0(loc.shape, "/Major_Road.shp"))
r <- raster(x = extent(landuse), nrows = 2039, ncols = 1819)
proj4string(r) <- CRS(proj4string(landuse))
r[] <- NA
mj.roads.r <- rasterize(mj.roads, r, field=1)
plot(mj.roads.r)
mj.roaddist.r <- distance(mj.roads.r)
class(mj.roaddist.r)
# Check:
plot(mj.roaddist.r)
lines(mj.roads)
writeRaster(x = mj.roaddist.r, filename = paste0(loc.shape, "/Major_roads_dist.tif"),
            format = "GTiff", overwrite = TRUE)

# We need to invert the raster for ease of interpretation
dist.mj.roads <- raster(paste0(loc.shape, "/Major_roads_dist.tif"))
dist.mj.roads[] <- abs(dist.mj.roads[] - max(dist.mj.roads[]))
plot(dist.mj.roads)
writeRaster(x = dist.mj.roads, filename = paste0(loc.shape, "/Major_Roads_dist_INVERT.tif"),
            format = "GTiff", overwrite = TRUE)


# Read Habitat Raster and Set SSF Parameters ------------------------------

# You can enter multiple covariates as .tif files here
covariates <-  list(paste0(loc.shape, "/Major_Roads_dist_INVERT.tif"))
                    
rp <- stack(x = covariates)

cov.names <- c("dist_major_road")
names(rp) <- cov.names 

## Plot the variables to check they are all present
plot(rp)

# Run again if you are doing playing with values in the loop
all.ssfResults <- NULL
ind.ssfResults <- NULL
# rm(data.1, ssf1, tr1)

# Specify the snakes to use. All adults. 
snake_ids <- c("AM006", "AM007", "AM015", "AM024", "AM026",
               "AM054", "AM059", "AF010", "AF017", "AF056", "AF058",
               "AF086", "AF096", "AF099")


## Create model sets to run for the iSSF analysis
mods <- list()

mods[[1]] <- case_ ~  log_sl * cos_ta + strata(step_id_)
mods[[2]] <- case_ ~  dist_major_road + dist_major_road:log_sl + dist_major_road:cos_ta +
  log_sl * cos_ta + strata(step_id_)


# SSF Loop per Individual -------------------------------------------------

i <- 0
ssfResults <- vector(mode = "list", length = length(snake_ids))
for(snake in snake_ids) {
  i <- i+1
   #snake <- "OPHA007"
  data.2 <- dplyr::filter(data.var, id == snake)
  
  # Creates the track data using package 'amt'
  tr1 <- mk_track(data.2, x, y, datetime, crs = crs.proj, id = id)
  
  # Calculates step and angle variables, then filters out non-movement (including due to GPS error)
  ssf1 <- tr1 %>%
    steps() %>%
    filter(sl_ > 5)  
  
  summary(ssf1)

  # Create random steps and extract covariates
  ssf1 <- ssf1 %>% 
    random_steps(n = 200) %>% 
    extract_covariates(rp, where = "end") 
  
  # Creates the transformed variables used in the models
  ssf1 <- ssf1 %>% 
    mutate(dist_major_road = log(dist_major_road),
           cos_ta = cos(ta_), 
           log_sl = log(sl_))
  
  # Set the individual ID for the current data set
  ssf1$id <- snake 
  
  # Runs a loop and fits an iSSF for all model sets
  ind.ssfResults <- NULL
  iter <- 0
  for(model in mods){
    iter <- iter + 1
    ssffit <-  fit_issf(ssf1, formula = model)
    results_ssf <- tidy(ssffit$model)
    results_ssf$id <- snake
    results_ssf$model <- paste0('model',iter)
    results_ssf$AIC <- AIC(ssffit$model)
    ind.ssfResults <- rbind(ind.ssfResults, results_ssf)
  }
  ssfResults[[i]] <- ind.ssfResults
} # end of indi loop

# Review Results ----------------------------------------------------------

all.ssfResults.2 <- do.call(rbind, ssfResults)

write_csv(path = "./Data/Secondary Table - ISSF_major_roads.csv", x = all.ssfResults)

# all.ssfResults <- read_csv(file = paste0(loc.output, "ssf_ALLresults.csv"))

# Create AIC table for all individuals
ssf.AIC.results.2 <- all.ssfResults.2 %>%
  group_by(model, id) %>%
  summarise(AIC = mean(AIC)) %>%
  arrange(AIC)

write_csv(path = paste0(loc.output, "ISSF_major_roads.csv"), x = ssf.AIC.results)

ssf.AIC.results.t2.2 <- ssf.AIC.results.2 %>% 
  group_by(id) %>% 
  mutate(AIC = ifelse(AIC < min(AIC)+2, paste0("*",
                                              round(AIC, digits = 2),
                                              "*"),
                                              round(AIC, digits = 2)))
# ssf.AIC.results.t2$AIC <- round(ssf.AIC.results$AIC)
table.3 <- reshape2::dcast(ssf.AIC.results.t2.2, model ~ id)

table.3$formula <-
c("log_sl*cos_ta+strata(step_id_)",
  "log_sl*cos_ta+strata(step_id_)+dist_major_road+dist_major_road:log_sl+dist_major_road:cos_ta")

table.3 <- table.3 %>% 
  dplyr::select(model, formula, contains("0")) %>% 
  rename("Model #" = model,
         "Model formula" = formula)

table.3

write_csv(path = "./Data/Table - ISSF models.csv", x = table.2)

# Return the top model for each individual
top.models <- ssf.AIC.results.2 %>%
  group_by(id) %>%
  filter(AIC == min(AIC)) %>%
  arrange(model)

top.models

write_csv(path = paste0(loc.output, "ssf_AICresults_ISSF_major_roads.csv"), x = top.models)

# Compare AIC values for each model by individual
ssf.AIC.results.2 %>% 
  group_by(id) %>% 
  mutate(best = ifelse(AIC == min(AIC), "*", NA)) %>% 
  ungroup() %>% 
  ggplot() +
  geom_col(aes(x = id, y = AIC, fill = model), position = position_dodge()) +
  geom_text(aes(x = id, y = AIC +100, label = best, group = model),
            position = position_dodge(width = 0.9),
            fontface = 2, size = 8) +
  scale_fill_manual(values = scico(n = 9, palette = "roma"),
                    labels = paste0("Model ", seq(1,9,1))) +
  theme_bw() +
  labs(x = "ID", fill = "Model") +
  theme(legend.position = c(0.9, 0.70),
        legend.background = element_blank(),
        axis.text.x = element_text(hjust = 0.5),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(hjust = 1, face = 2),
        plot.title = element_text(face = 4))

ggsave(file = paste0("./Figures/Secondary Figure - individual ISSF.png"), width = 200, height = 120,
       dpi = 600, units = "mm")
ggsave(file = paste0("./Figures/Secondary Figure - individual ISSF.pdf"), width = 200, height = 120,
       units = "mm")

all.results <- rbind(all.ssfResults, all.ssfResults.2)

write_csv(path = "./Data/all snakes road results ISSF.csv", x = all.results)

# Load in this new csv for plotting, it has two new columns
# season and sex
new.data <- read_csv(file = paste0(loc.data, "all snakes road results ISSF.csv"),
                       locale = locale(tz = "Asia/Bangkok"))

# Filter the data for plotting
issf.fil <- new.data %>%
  dplyr::filter((model == 'model1' & id %in% c("AM006", "AM007", "AM015", "AM024", "AM026",
                                               "AM054", "AM059", "AF010", "AF017", "AF056", "AF058",
                                               "AF086", "AF096", "AF099")) |
                  (model == 'model2' & id %in% c("AM006", "AM007", "AM015", "AM024", "AM026",
                                                 "AM054", "AM059", "AF010", "AF017", "AF056", "AF058",
                                                 "AF086", "AF096", "AF099")) |
                  term %in% cov.names) %>% 
  mutate(term = str_to_title(sub("dist_", "", term)))

road.mod <- issf.fil %>%
  dplyr::filter((term == 'Major_road' & model %in% c("model2")))

write_csv(habitats, "Selection by adults.csv")

colors1 <- c("#00aedb", "#f37735")

road.plot <- ggplot(road.mod) +
  geom_hline(yintercept = 0, linetype = 2, alpha = 0.5) +
  geom_point(aes(x = id, y = estimate, colour = season, group = season), size = 2, 
             position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(x = id, ymin = conf.low, ymax = conf.high, colour = season, group = season),
                alpha = 0.85, linetype = 7, size = 1, width = 0.35, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = colors1) +
  labs(y = expression(beta), x = "Snake") +
  facet_wrap(sex ~ ., scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 55, hjust = 0.5, vjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = 4, hjust = 0),
        legend.position = "none",
        legend.title = element_text(face = 2),
        legend.text = element_text(lineheight = 1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 1, face = 2),
        axis.title.x = element_text(angle = 0, hjust = 1, face = 2))

road.plot

ggsave(file = paste0("./Figures/Individual season ISSF plot.png"), width = 200, height = 120,
       dpi = 600, units = "mm")
ggsave(file = paste0("./Figures/Individual season ISSF plots.pdf"), width = 200, height = 120,
       units = "mm")

#end