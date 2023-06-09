## 07 - SDM evaluation and sample summaries

# Libraries ---------------------------------------------------------------

library(PRROC)
library(stringr)
library(raster)
library(ENMeval)
library(scico)
library(ggplot2)
library(ggrepel)
library(dplyr)

# PRRC and ROC Evaluation -------------------------------------------------

# a nested loop that goes into every model training scenario, and for every
# species, looks at every model to calculate ROC and PRRC AUC values against all
# testing datasets.
# inappropriate tests are removed after the loop
# (eg. Flickr test data testing a Flickr supplemented model)

runScens <- c("GBIF",
              "FlickrSupp",
              "GBIFran",
              "FlickrSuppRan")

test_list <- c("GBIF", "Flickr", "GBIFran")
i <- 0
i.inter <- 0
internal_results <- vector(mode = "list", length = 18*3)
PRresults <- list()
PRcurves <- list()
ROCcurves <- list()
dir.create("./Outputs/PRRC/")
dir.create("./Outputs/ROC/")
for(run in runScens){
  
  print(paste0("--- ", run, " ---"))
  
  if(run == "GBIF"){
    loc.model <- "./Outputs/GBIF_mod/"
    
  } else if(run == "FlickrSupp"){
    loc.model <- "./Outputs/FlickrSupp_mod/"
    
  } else if(run == "GBIFran"){
    loc.model <- "./Outputs/GBIFran_mod/"
    
  } else if(run == "FlickrSuppRan"){
    loc.model <- "./Outputs/FlickrSuppRan_mod/"
  }
  
  species_list <- unique(str_extract(list.files(path = loc.model,
                                                pattern = "evalMods.Rdata"), "^...."))
  
  for(sp in species_list){
    
    print(paste0("------ ", sp, " ------"))
    
    species_files <- list.files(loc.model, pattern = sp)
    
    rasters <- raster::stack(paste0("./Outputs/", sp, "_rasters.grd"))
    
    test_flickr <- readRDS(paste0("./Outputs/", sp, "_testFlickr.Rds"))
    test_gbif <- readRDS(paste0("./Outputs/", sp, "_testGBIF.Rds"))
    test_gbifran <- readRDS(paste0("./Outputs/", sp, "_testGBIFran.Rds"))
    bg <- readRDS(paste0("./Outputs/", sp, "_bg.Rds"))
    
    train_flickr <- readRDS(paste0("./Outputs/", sp, "_trainFlickr.Rds"))
    train_gbif <- readRDS(paste0("./Outputs/", sp, "_trainGBIF.Rds"))
    train_gbifran <- readRDS(paste0("./Outputs/", sp, "_trainGBIFran.Rds"))
    train_FlickrRan <- readRDS(paste0("./Outputs/", sp, "_trainFlickrRan.Rds"))

    load(paste0(loc.model,
                grep(pattern = "evalMods.Rdata", x = species_files, value = TRUE)))
    
    evalTbl <- readRDS(paste0(loc.model,
                grep(pattern = "evalTbl.Rds", x = species_files, value = TRUE)))
    
    internal.AUC <- mean(evalTbl$avg.test.AUC)
    internal.varAUC <- mean(evalTbl$var.test.AUC)
    internal <- data.frame("species" = sp,
                           "train" = run,
                           "internalAUC" = internal.AUC, 
                           "internalAUCvar" = internal.varAUC)
    
    i.inter <- i.inter+1
    
    internal_results[[i.inter]] <- internal
    
    for(mod in names(evalMods)){
      cat(mod, sep = "| ")
      
      i <- i+1
      
      model <- evalMods[names(evalMods) == mod]
      
      pred <- maxnet.predictRaster(mod = model[[1]], env = rasters,
                                   type = "logistic", clamp = FALSE)
      
      prAUC <- NULL
      rocAUC <- NULL
      prdf_comb <- NULL
      rocdf_comb <- NULL
      for(test in test_list){
        
        if(test == "GBIF"){
          test_vals <- extract(x = pred, y = test_gbif[,2:1])
        } else if(test == "Flickr"){
          test_vals <- extract(x = pred, y = test_flickr[,2:1])
        } else if(test == "GBIFran"){
          test_vals <- extract(x = pred, y = test_gbifran[,2:1])
        }
        bg.vals <- extract(x = pred, y = bg)
        
        try(
          rm(pr)
        )
        
        try(
          pr <- pr.curve(scores.class0 = test_vals, weights.class0 = test_vals,
                         curve = TRUE)
        )
        
        if(exists("pr")){
          prdf <- as.data.frame(pr$curve)
          names(prdf) <- c("Recall", "Precision", "Threshold")
          prdf$train <- run
          prdf$test <- test
          prdf$model <- mod
          prdf$species <- sp
          
          prAUC <- c(prAUC, pr$auc.integral)
          # used for plot
          prdf_comb <- rbind(prdf_comb, prdf)
        } else {
          prdf <- data.frame("Recall" = NA, "Precision" = NA,
                             "Threshold" = NA)
          prdf$train <- run
          prdf$test <- test
          prdf$model <- mod
          prdf$species <- sp
          
          prAUC <- c(prAUC, NA)
          
          # used for plot
          prdf_comb <- rbind(prdf_comb, prdf)
        }
        
        roc <- roc.curve(scores.class0 = test_vals, scores.class1 = bg.vals, curve = TRUE)
        
        rocdf <- as.data.frame(roc$curve)
        names(rocdf) <- c("FPR", "Sensitivity", "Threshold")
        rocdf$train <- run
        rocdf$test <- test
        rocdf$model <- mod
        rocdf$species <- sp
        
        rocAUC <- c(rocAUC, roc$auc)
        # used for plot
        rocdf_comb <- rbind(rocdf_comb, rocdf)
        
      } # end of test data loop
      
      PRcurves[[i]] <- prdf_comb
      ROCcurves[[i]] <- rocdf_comb
      
      spPRresults <- data.frame(test_list, prAUC, rocAUC)
      spPRresults$train <- run
      spPRresults$model <- mod
      spPRresults$species <- sp
      
      PRresults[[i]] <- spPRresults
      
    } # end of model loop
  } # end of species loop
} # end of run loop
beepr::beep(1)

# saveRDS(PRresults, file = "./Outputs/PRROCresults.Rds")
# saveRDS(PRcurves, file = "./Outputs/PRRCcurves.Rds")
# saveRDS(ROCcurves, file = "./Outputs/ROCcurves.Rds")

# Summary plots for ROC and PRRC ---------------------------------------------------
PRresults <- readRDS(file = "./Outputs/PRROCresults.Rds")

prres <- bind_rows(PRresults)
prres$train[prres$train == "GBIF"] <- "GBIF geographic space"
prres$train[prres$train == "FlickrSupp"] <- "Geo. Flickr supplemented"
prres$train[prres$train == "FlickrSuppRan"] <- "Random Flickr supplemented"
prres$train[prres$train == "GBIFran"] <- "GBIF random"
prres$train <- str_wrap(prres$train, width = 10)

prres$test_list <- as.character(prres$test_list)
prres$test_list[prres$test_list == "GBIF"] <- "GBIF geo-independent"
prres$test_list[prres$test_list == "FlickrSupp"] <- "Flickr data"
prres$test_list[prres$test_list == "GBIFran"] <- "GBIF random"

prres %>%
  group_by(species, test_list, train) %>%
  summarise(mPRRC = mean(prAUC, na.rm = TRUE),
            mROC = mean(rocAUC, na.rm = TRUE)) %>%
  print(n = 800)

unique(prres$train)
unique(prres$test)

means_PRRC_ROC <- prres %>%
  filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
  filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
  filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "GBIF random") %>%
  filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "GBIF geo-independent") %>%
  filter(!train == "GBIF\\nrandom" | !test_list == "GBIF geo-independent") %>%
  filter(!train == "GBIF\\ngeographic\\nspace" | !test_list == "GBIF random") %>% 
  group_by(test_list, train) %>%
  summarise(mPRRC = mean(prAUC, na.rm = TRUE),
            mROC = mean(rocAUC, na.rm = TRUE)) %>% 
  mutate(test_list_f = factor(test_list,
                              levels = c("GBIF random", "GBIF geo-independent", "Flickr")))

strip_names <- c(
  "Test:\\nFlickr data",
  "Test:\\nGBIF geo-independent sample",
  "Test:\\nGBIF random sample")

names(strip_names) <- c("Flickr",
                        "GBIF geo-independent",
                        "GBIF random")

(ROCaucSummary_plot <- prres %>%
    filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
    filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
    filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "GBIF random") %>%
    filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "GBIF geo-independent") %>%
    filter(!train == "GBIF\\nrandom" | !test_list == "GBIF geo-independent") %>%
    filter(!train == "GBIF\\ngeographic\\nspace" | !test_list == "GBIF random") %>%  
    mutate(test_list_f = factor(test_list,
                                levels = c("GBIF random", "GBIF geo-independent", "Flickr"))) %>% 
    ggplot() +
    geom_violin(aes(x = train, y = rocAUC), 
                colour = NA, fill = "black", alpha = 0.2) +
    geom_point(aes(x = train, y = rocAUC, colour = species),
               position = position_dodge(width = 0.5)
    ) +
    geom_text(data = means_PRRC_ROC,
              aes(x = train, y = 1.1, label = round(mROC, digits = 3))) +
    facet_wrap(.~test_list_f, scales = "free_x",
               labeller = as_labeller(strip_names)
    ) +
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, face = 4),
          strip.background = element_blank(),
          axis.text.x = element_text(hjust = 0.5, angle = 0)
    ) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_colour_scico_d(palette = "roma") +
    labs(x = "Training dataset", y = "ROC AUC", colour = "Species")
)

ggsave(filename = paste0("./Figures/FIGURE - ROC Summary.png"),
       plot = ROCaucSummary_plot, width = 24, height = 16, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0("./Figures/FIGURE - ROC Summary.pdf"),
       plot = ROCaucSummary_plot, width = 24, height = 16, units = "cm",
       device = "pdf")

(PRRCaucSummary_plot <- 
    prres %>%
    filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
    filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "Flickr") %>%
    filter(!train == "Geo.\\nFlickr\\nsupplemented" | !test_list == "GBIF random") %>%
    filter(!train == "Random\\nFlickr\\nsupplemented" | !test_list == "GBIF geo-independent") %>%
    filter(!train == "GBIF\\nrandom" | !test_list == "GBIF geo-independent") %>%
    filter(!train == "GBIF\\ngeographic\\nspace" | !test_list == "GBIF random") %>%  
    mutate(test_list_f = factor(test_list,
                                levels = c("GBIF random", "GBIF geo-independent", "Flickr"))) %>% 
    ggplot() +
    geom_violin(aes(x = train, y = prAUC), 
                colour = NA, fill = "black", alpha = 0.2) +
    geom_point(aes(x = train, y = prAUC, colour = species),
               position = position_dodge(width = 0.5)
    ) +
    geom_text(data = means_PRRC_ROC,
              aes(x = train, y = 1.1, label = round(mPRRC, digits = 3))) +
    facet_wrap(.~test_list_f, scales = "free_x",
               labeller = as_labeller(strip_names)
    ) +
    theme_bw() +
    theme(strip.text = element_text(hjust = 0, face = 4),
          strip.background = element_blank(),
          axis.text.x = element_text(hjust = 0.5, angle = 0),
    ) +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.2)) +
    scale_colour_scico_d(palette = "roma") +
    labs(x = "Training dataset", y = "PR-RC AUC", colour = "Species")
)

ggsave(filename = paste0("./Figures/FIGURE - PRRC Summary.png"),
       plot = PRRCaucSummary_plot, width = 24, height = 16, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0("./Figures/FIGURE - PRRC Summary.pdf"),
       plot = PRRCaucSummary_plot, width = 24, height = 16, units = "cm",
       device = "pdf")

# Sample size extraction --------------------------------------------------

mcp <- function(xy) {
  xy <- as.data.frame(sp::coordinates(xy))
  coords.t <- chull(xy[, 1], xy[, 2])
  xy.bord <- xy[coords.t, ]
  xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
  return(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(as.matrix(xy.bord))), 1))))
}

land <- readRDS(file = "landmass.Rds")
prj_Nasia <- CRS("+proj=aea +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
prj_Africa <- CRS("+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
prj_Samerica <- CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs ")

i.samp <- 0
samples <- vector(mode = "list", length = 18)
internal_results <- vector(mode = "list", length = 3*18)

species_list <- unique(str_extract(list.files(path = "./Outputs/GBIFran_mod/",
                                              pattern = "evalMods.Rdata"), "^...."))

# loop that pulls out the sample sizes of training and testing datasets as well
# as calculating estimated area covered for each species

for(sp in species_list){
  
  i.samp <- i.samp + 1 
  
  print(paste0("------ ", sp, " ------"))
  
  test_flickr <- readRDS(paste0("./Outputs/", sp, "_testFlickr.Rds"))
  test_gbif <- readRDS(paste0("./Outputs/", sp, "_testGBIF.Rds"))
  test_gbifran <- readRDS(paste0("./Outputs/", sp, "_testGBIFran.Rds"))
  
  train_flickr <- readRDS(paste0("./Outputs/", sp, "_trainFlickr.Rds"))
  train_gbif <- readRDS(paste0("./Outputs/", sp, "_trainGBIF.Rds"))
  train_gbifran <- readRDS(paste0("./Outputs/", sp, "_trainGBIFran.Rds"))
  train_FlickrRan <- readRDS(paste0("./Outputs/", sp, "_trainFlickrRan.Rds"))
  
  occs <- rbind(train_flickr, train_gbif)
  
  area_mcp <- mcp(occs[,2:1])
  proj4string(area_mcp) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  
  inter_mcp <- gIntersection(land, area_mcp)
  
  if(mean(area_mcp@bbox[c(1,3)]) > 50 & mean(area_mcp@bbox[c(1,3)]) < 150){
    proj_mcp <- spTransform(inter_mcp, prj_Nasia)
  } else if(mean(area_mcp@bbox[c(1,3)]) > -40 & mean(area_mcp@bbox[c(1,3)]) < 40){
    proj_mcp <- spTransform(inter_mcp, prj_Africa)
  } else if(mean(area_mcp@bbox[c(1,3)]) > -120 & mean(area_mcp@bbox[c(1,3)]) < -25){
    proj_mcp <- spTransform(inter_mcp, prj_Samerica)
  }
  
  occ_area <- gArea(proj_mcp)/1000000
  
  samplesizes <- data.frame("species" = sp,
                      "Train_FlickrSuppSize" = dim(train_flickr)[1],
                      "Train_FlickrRanSuppSize" = dim(train_FlickrRan)[1],
                      "Train_GBIFGeoSize" = dim(train_gbif)[1],
                      "Train_GBIFRanSize" = dim(train_gbifran)[1],
                      "Test_FlickrSuppSize" = dim(test_flickr)[1],
                      "Test_GBIFGeoSize" = dim(test_gbif)[1],
                      "Test_GBIFRanSize" = dim(test_gbifran)[1],
                      "Total_GBIF" = dim(test_gbifran)[1] + dim(train_gbifran)[1],
                      "Occ_Area_km2" = occ_area)
  
  samplesizes <- samplesizes %>% 
    mutate(All_occ = Total_GBIF + abs(Train_GBIFGeoSize-Train_FlickrSuppSize))
  
  samples[[i.samp]] <- samplesizes
  
} # sp loop

# saveRDS(samples, file = "./Outputs/SampleSummaryData.Rds")
# samples <- readRDS(file = "./Outputs/SampleSummaryData.Rds")

# Sample summary and plot ------------------------------------------------------------

samples <- readRDS(file = "./Outputs/SampleSummaryData.Rds")

ss_full <- bind_rows(samples)

(sample_plot <- ss_full %>% 
  mutate(neg_flickr = Train_GBIFGeoSize-Train_FlickrSuppSize,
         percent_flickr  = Test_FlickrSuppSize / All_occ *100 ) %>%
  arrange(Total_GBIF) %>%
  mutate(species = factor(species, levels = (species))) %>% 
  ggplot() +
  geom_segment(aes(x = species, xend = species,
                   y = Total_GBIF, yend = 0, colour = "GBIF"),
               size = 5) +
  geom_segment(aes(x = species, xend = species,
                   y = neg_flickr, yend = 0, colour = "Flickr"),
               size = 5) +
  geom_text(aes(x = species, y = -100, label = All_occ),
            size = 4) +
  coord_flip() +
  theme_bw() +
  scale_colour_scico_d(palette = "roma", begin = 0.075, end = 0.9,
                       direction = -1) +
  scale_y_continuous(breaks = seq(-100, 600, 100), 
                     labels = c(100, seq(0, 600, 100))) +
  theme(legend.position = c(0.9, 0.1),
        legend.background = element_blank(),
        axis.title.y = element_text(angle = 0, hjust = 0),
        plot.subtitle = element_text(face = 3)
        ) +
  labs(x = "Species", y = "Number of locations", colour = "Data source",
       subtitle = "  Total")
)

ggsave(filename = paste0("./Figures/FIGURE - Sample Summary.png"),
       plot = sample_plot, width = 24, height = 16, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0("./Figures/FIGURE - Sample Summary.pdf"),
       plot = sample_plot, width = 24, height = 16, units = "cm",
       device = "pdf")

qqnorm(ss_full$Occ_Area_km2)
qqnorm(ss_full$All_occ)
shapiro.test(ss_full$Occ_Area_km2)
shapiro.test(ss_full$All_occ)

sp_results <- cor.test(ss_full$Occ_Area_km2, ss_full$Total_GBIF, method = "spearman")

breaks <- c(250000, 500000, 750000, 1000000, 2500000, 5000000, 7500000,
  10000000, 25000000)

(area_sample_plot <- 
    ggplot(ss_full, aes(x = Occ_Area_km2, y = All_occ)) +
    geom_point(aes(colour = species)) +
    geom_text_repel(aes(label = species)) +
    theme_bw() +
    scale_colour_scico_d(palette = "roma") +
    scale_y_log10(breaks = c(1, 10, 25, 50, 75, 100, 250, 500, 750, 1000), limits = c(1, 1000),
                  minor_breaks = NULL) +
    scale_x_log10(breaks = breaks,
                  labels = breaks/1000000,
                  minor_breaks = NULL,
                  limits = c(250000, 28000000),
                  expand = c(0,0)) +
    theme(legend.position = "none") +
    labs(x = expression(paste("Land mass covered by occurrence MCP (Mill. ", km^2, ")")),
         y = "Total number of occurrences", colour = "Species",
         title = "Spearman's correlation: ",
         subtitle = paste0("S = ", round(digits = 3, sp_results$statistic),
                           ", p-value = ", round(digits = 2, sp_results$p.value),
                           ", rho = ", round(digits = 2, sp_results$estimate)))
)

ggsave(filename = paste0("./Figures/FIGURE - Sample Area Plot.png"),
       plot = area_sample_plot, width = 24, height = 16, units = "cm",
       dpi = 300, device = "png")
ggsave(filename = paste0("./Figures/FIGURE - Sample Area Plot.pdf"),
       plot = area_sample_plot, width = 24, height = 16, units = "cm",
       device = "pdf")

# PRRC and ROC Curves - NOT RUN -----------------------------------------------------

# loop to generate individual curve plots for ROC and PRRC values. Not needed
# and takes a long time.

# PRcurves <- readRDS(file = "./Outputs/PRRCcurves.Rds")
# ROCcurves <- readRDS(file = "./Outputs/ROCcurves.Rds")
# 
# PRcomplete <- bind_rows(PRcurves)
# 
# PRcomplete <- PRcomplete %>%
#   filter(!train == "FlickrSupp" | !test == "Flickr") %>%
#   filter(!train == "GBIFran" | !test == "GBIF") %>%
#   filter(!train == "GBIF" | !test == "GBIFran") 
# 
# label_loc <- PRcomplete %>%
#   group_by(test, model) %>%
#   top_n(Recall, n = 1) %>%
#   slice(1)
# 
# PRcomplete$test[PRcomplete$test == "GBIF"] <- "Test: GBIF Env. Space"
# PRcomplete$test[PRcomplete$test == "GBIFran"] <- "Test: GBIF Random"
# PRcomplete$test[PRcomplete$test == "Flickr"] <- "Test: Flickr"
# 
# PRcomplete$train[PRcomplete$train == "GBIF"] <- "Train: GBIF Env. Space"
# PRcomplete$train[PRcomplete$train == "GBIFran"] <- "Train: GBIF Random"
# PRcomplete$train[PRcomplete$train == "FlickrSupp"] <- "Train: Flickr Supp."
# 
# ROCcomplete <- bind_rows(ROCcurves)
# 
# ROCcomplete <- ROCcomplete %>%
#   filter(!train == "FlickrSupp" | !test == "Flickr") %>%
#   filter(!train == "GBIFran" | !test == "GBIF") %>%
#   filter(!train == "GBIF" | !test == "GBIFran") 
# 
# label_loc <- ROCcomplete %>%
#   group_by(test, model) %>%
#   top_n(FPR, n = 1) %>%
#   slice(1)
# 
# ROCcomplete$test[ROCcomplete$test == "GBIF"] <- "Test: GBIF Env. Space"
# ROCcomplete$test[ROCcomplete$test == "GBIFran"] <- "Test: GBIF Random"
# ROCcomplete$test[ROCcomplete$test == "Flickr"] <- "Test: Flickr"
# 
# ROCcomplete$train[ROCcomplete$train == "GBIF"] <- "Train: GBIF Env. Space"
# ROCcomplete$train[ROCcomplete$train == "GBIFran"] <- "Train: GBIF Random"
# ROCcomplete$train[ROCcomplete$train == "FlickrSupp"] <- "Train: Flickr Supp."
# 
# for(sp in unique(ROCcomplete$species)){
#   
#   for(run in unique(PRcomplete$train)){
#     
#     prrc_plot <-
#       PRcomplete %>% 
#       filter(species == sp, train == run) %>% 
#       ggplot(aes(x = Recall, y = Precision, colour = Threshold)) +
#       geom_line(aes(group = model), size = 1.25, alpha = 0.5) +
#       facet_grid(train~test) +
#       theme_bw() +
#       scale_colour_scico(palette = "roma") +
#       scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
#       scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.05)) +
#       theme(
#         legend.position = c(0.1, 0.25),
#         legend.background = element_blank(),
#         plot.title = element_text(face = 2),
#         plot.subtitle = element_text(face = 3),
#         strip.background = element_blank(),
#         strip.text = element_text(face = 4, hjust = 0),
#         strip.text.y = element_blank()
#       ) +
#       labs(title = paste("Precision-Recall Curves:", sp),
#            subtitle = run,
#            size = "Test Set")
#     NULL
#     
#     prrc_plot
#     
#     ggsave(filename = paste0("./Figures/", sp, "_", sub("Train: ", "", run), "PRRC.png"),
#            plot = prrc_plot, width = 18, height = 12, units = "cm",
#            dpi = 300, device = "png")
#     
#     roc_plot <- 
#       ROCcomplete %>% 
#       filter(species == sp, train == run) %>% 
#       ggplot(aes(x = FPR, y = Sensitivity, colour = Threshold)) +
#       geom_line(aes(group = model), size = 1.25, alpha = 0.5) +
#       facet_grid(test~train) +
#       theme_bw() +
#       scale_colour_scico(palette = "roma") +
#       scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
#       scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.05)) +
#       theme(
#         legend.position = c(0.1, 0.25),
#         legend.background = element_blank(),
#         plot.title = element_text(face = 2),
#         plot.subtitle = element_text(face = 3),
#         strip.background = element_blank(),
#         strip.text = element_text(face = 4, hjust = 0),
#         strip.text.y = element_blank()
#       ) +
#       labs(title = paste("Receiver Operating Characteristic Curve:", sp),
#            subtitle = run,
#            size = "Test Set")
#     
#     ggsave(filename = paste0("./Figures/", sp, "_", sub("Train: ", "", run), "_ROC.png"),
#            plot = roc_plot, width = 18, height = 12, units = "cm",
#            dpi = 300, device = "png")
#     
#   } # run/train loop
# } # species loop
