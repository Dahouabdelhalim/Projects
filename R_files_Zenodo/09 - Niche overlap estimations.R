## 09 - Niche overlap estimations

# Libraries ---------------------------------------------------------------

library(pracma)
library(stringr)
library(raster)
library(dplyr)
library(ENMeval)
library(ggplot2)
library(scico)
library(reshape2)

# Loop parameters ---------------------------------------------------------

runScens <- c("GBIF",
              "FlickrSupp",
              "GBIFran",
              "FlickrSuppRan")

list.files(list.files("./Outputs/", pattern = "_mod", full.names = TRUE)[1], pattern = "Mods.Rdata")

species_list <- unique(str_extract(list.files(path = "./Outputs/GBIF_mod/",
                                              pattern = "evalMods.Rdata"), "^...."))


# Loop to run niche overlaps ----------------------------------------------

# an inefficient loop that calculates the niche overalp for all models for all
# scenarios for all species.

# comparisons between models with different regulisations are removed after
# niche overlap calulations.

nc.i <- 0
i <- 0
niche_compare <- vector(mode = "list", length = length(species_list))
for(sp in species_list){
  print(paste0("------ ", sp, " ------"))
  
  nc.i <- nc.i + 1
  
  rasters <- raster::stack(paste0("./Outputs/", sp, "_rasters.grd"))
  
  for(run in runScens){
    print(paste0("--- ", run, " ---"))
    
    if(run == "GBIF"){
      loc.model <- "./Outputs/GBIF_mod/"
      species_files <- list.files(loc.model, pattern = sp)
      load(paste0(loc.model,
                  grep(pattern = "evalMods.Rdata", x = species_files, value = TRUE)))
      evalMods_GBIF <- evalMods
    } else if(run == "FlickrSupp"){
      loc.model <- "./Outputs/FlickrSupp_mod/"
      species_files <- list.files(loc.model, pattern = sp)
      load(paste0(loc.model,
                  grep(pattern = "evalMods.Rdata", x = species_files, value = TRUE)))
      evalMods_FlickrSupp <- evalMods
    } else if(run == "GBIFran"){
      loc.model <- "./Outputs/GBIFran_mod/"
      species_files <- list.files(loc.model, pattern = sp)
      load(paste0(loc.model,
                  grep(pattern = "evalMods.Rdata", x = species_files, value = TRUE)))
      evalMods_GBIFran <- evalMods
    } else if(run == "FlickrSuppRan"){
      loc.model <- "./Outputs/FlickrSuppRan_mod/"
      species_files <- list.files(loc.model, pattern = sp)
      load(paste0(loc.model,
                  grep(pattern = "evalMods.Rdata", x = species_files, value = TRUE)))
      evalMods_FlickrSuppRan <- evalMods
    }
    
    try(rm(list=ls(pattern = "_Modprediction")))
    
    mod_names <- NULL
    for(mod in names(evalMods)){
      
      cat(mod, sep = "| ")
      
      i <- i+1
      
      model <- evalMods[names(evalMods) == mod]
      mod_names <- c(mod_names, mod)
      prediction <- maxnet.predictRaster(mod = model[[1]], env = rasters,
                                   type = "logistic", clamp = FALSE)
      assign(paste0(mod, "_", run, "_Modprediction"), prediction)
    } # mod loop
      
    pred_list <- mget(ls(pattern = "_Modprediction"))
    pred_stack <- stack(pred_list)
    assign(paste0(run, "_Stack"), pred_stack)
      
  } # run loop
  
  GBIF_test <- stack(FlickrSupp_Stack, GBIF_Stack)
  
  GBIF_res <- calc.niche.overlap(GBIF_test)
  mods_rown <- row.names(GBIF_res)
  GBIF_res <- as.data.frame(GBIF_res)
  GBIF_res$mod <- mods_rown
  
  GBIF_trim_res <- GBIF_res %>% 
    dplyr::filter(str_detect(mod, pattern = "GBIF")) %>% 
    dplyr::select(-contains("GBIF"))
  
  geo_niche <- diag(as.matrix(GBIF_trim_res[,-17]))  
  
  
  GBIFran_test <- stack(FlickrSuppRan_Stack, GBIFran_Stack)
  
  GBIFran_res <- calc.niche.overlap(GBIFran_test)
  mods_rown <- row.names(GBIFran_res)
  GBIFran_res <- as.data.frame(GBIFran_res)
  GBIFran_res$mod <- mods_rown
  
  GBIFran_trim_res <- GBIFran_res %>% 
    dplyr::filter(str_detect(mod, pattern = "GBIF")) %>% 
    dplyr::select(-contains("GBIF"))
  
  ran_niche <- diag(as.matrix(GBIFran_trim_res[,-17]))  
  
  niche_compare[[nc.i]] <- data.frame("species" = sp,
                                      "mod" = mod_names,
                                      "geo" = geo_niche,
                                      "ran" = ran_niche)
  
} # species loop

# saveRDS(file = "./Outputs/niche_results.Rds", object = niche_compare)

# Summarise and plot niche overlap values ---------------------------------

niche_compare <- readRDS("./Outputs/niche_results.Rds")

niche_df <- bind_rows(niche_compare)

summary.labels <- niche_df %>% 
  melt() %>%
  group_by(species) %>% 
  mutate(mean.n = format(nsmall = 3, round(mean(value), digits = 3)),
         se.n = format(nsmall = 3, round(std_err(value), digits = 3))) %>% 
  slice(1)

summary.labels$mean.n
summary.labels$se.n

(niche.plot <- niche_df %>% 
  reshape2::melt() %>% 
  arrange(desc(value)) %>% 
  mutate(species.fac = factor(species, levels = unique(species))) %>% 
  ggplot() +
  geom_point(aes(x = species.fac, y = value, colour = variable),
             position = position_dodge(0.5)) +
  geom_text(data = summary.labels,
            aes(x = species, y = 1.03,
                label = paste0(mean.n, " Â± ", se.n)),
             position = position_dodge(0.5), hjust = 1) +
  scale_colour_scico_d(palette = "roma", begin = 0.075, end = 0.9,
                       direction = -1, labels = c("Geo", "Random")) +
  scale_y_continuous(limits = c(0.825, 1.035), breaks = seq(0.85, 1, 0.05),
                     expand = c(0,0)) +
  coord_flip() +
  theme_bw() +
  labs(x = "Species", y = "Schoener's D Niche Overlap", colour = "GBIF Datasets") +
  theme(legend.position = c(0.12, 0.12)) +
  NULL)

ggsave(filename = paste0("./Figures/FIGURE - Niche Overlap.png"),
       plot = niche.plot, width = 24, height = 16, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0("./Figures/FIGURE - Niche Overlap.pdf"),
       plot = niche.plot, width = 24, height = 16, units = "cm",
       device = "pdf")
