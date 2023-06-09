## 06 - Sub-sampling and SDModelling

# Libraries ---------------------------------------------------------------

library(blockCV)
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(raster)
library(here)
library(scico)
library(dplyr)
library(stringr)


# Folder locations --------------------------------------------------------

loc.raster <- "./Data/Raster Data/"
loc.figure <- "./Figures/"

# Maxent modelling loop ---------------------------------------------------

# loop to run for selected and random species, separated because of different
# data locations
random <- c("non-random", "random")

for(ran in random){
  
  if(ran == "non-random"){
    loc.species <- "./Data/SpeciesData/"
    species_list <- unique(str_extract(list.files(path = loc.species,
                                                  pattern = "all.csv"), "^...."))
  } else {
    loc.species <- "./Data/RandomSpeciesData/"
    species_list <- unique(str_extract(list.files(path = loc.species,
                                                  pattern = "all.csv"), "^...."))
  }
  
  # Species loop ----
  
  for(sp in species_list){
    set.seed(1)
    
    # Load occurrencs ----
    
    file <- grep(list.files(path = loc.species, pattern = sp), pattern = "all.csv", value = TRUE)
    df <- read.csv(file = paste0(loc.species, file), stringsAsFactors = FALSE)
    df <- df %>%
      group_by(source) %>%
      distinct(latitude, longitude, .keep_all = TRUE)
    
    df <- df[!duplicated(df$latitude) & !duplicated(df$longitude),]
    df <- df[!df$latitude == 0 | !df$longitude == 0,]
    df <- df[!is.na(df$latitude) | !is.na(df$longitude),]

    # Load rasters ----
    
    # get centre to select which rasters to use
    if(mean(range(df$longitude)) > 50 & mean(range(df$longitude)) < 150){
      # read in ASIA
      rast.data <- stack(paste0(loc.raster, list.files(loc.raster, pattern = "_005.tif")))
    } else if(mean(range(df$longitude)) > -25 & mean(range(df$longitude)) < 75){
      # read Africa
      rast.data <- stack(paste0(loc.raster, list.files(loc.raster, pattern = "_africa005.tif")))
    } else if(mean(range(df$longitude)) > -120 & mean(range(df$longitude)) < -25){
      # read S America
      rast.data <- stack(paste0(loc.raster, list.files(loc.raster, pattern = "_samerica005.tif")))
    }
    # rename them to something consistent
    names(rast.data) <- c("bio1_005", "bio12_005", "bio15_005",
                          "bio2_005", "bio7_005", "dem_005", "hfootprint_005")
    # remove layer to be used for bias file
    model.rasters <- rast.data[[names(rast.data)[!names(rast.data) == "hfootprint_005"]]]
    
    # Spatially thin occurrences ----
    
    # thinning par is the same as the res of the climate data
    output <- spThin::thin(df, 'latitude', 'longitude', 'name', thin.par = 4.65,
                           reps = 100, locs.thinned.list.return = TRUE,
                           write.files = FALSE, verbose = FALSE)
    # find the iteration that returns the max number of occurrences
    maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
    # if there's more than one max, pick the first one
    maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
    # subset occs to match only thinned occs
    df <- df[as.numeric(rownames(maxThin)),]

    # Generate extent area ----
    occs.xy <- df[,2:1]
    
    # generate MCP around the occureneces
    mcp <- function(xy) {
      xy <- as.data.frame(sp::coordinates(xy))
      coords.t <- chull(xy[, 1], xy[, 2])
      xy.bord <- xy[coords.t, ]
      xy.bord <- rbind(xy.bord[nrow(xy.bord), ], xy.bord)
      return(sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(as.matrix(xy.bord))), 1))))
    }
    bgExt <- mcp(occs.xy)
    
    # and a buffer around that MCP - now based on half mean distance betweem points
    bgExt <- rgeos::gBuffer(bgExt, width = mean(dist(occs.xy))/2 )
    proj4string(bgExt) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    # crop the environmental rasters by the background extent shape
    envsBgCrop <- raster::crop(model.rasters, bgExt)
    # writeRaster(x = envsBgCrop,
    #             filename = paste0(loc.model, "/envsBgCrop.grd"),
    #             format = "raster", overwrite = TRUE)
    
    # mask the background extent shape from the cropped raster
    envsBgMsk <- raster::mask(envsBgCrop, bgExt)
    
    # Ensure occurrences have environmental data ----
    # extract environmental values at occ grid cells
    locs.vals <- raster::extract(envsBgMsk, df[, c('longitude', 'latitude')])
    # remove occs without environmental values
    df <- df[complete.cases(locs.vals),]
    
    # Generating background points bias file -----------------
    
    fp <- rast.data[["hfootprint_005"]]
    fp.BgCrop <- raster::crop(fp, bgExt)
    mask.fp <- raster::mask(fp.BgCrop, bgExt)
    mask.fp <- mask(x = mask.fp, mask = envsBgMsk[[1]], updatevalue = 0, maskvalue = NA,
                    inverse = FALSE)
    mask.fp[is.na(mask.fp[])] <- 0
    
    bg.xy <- data.frame(xyFromCell(mask.fp,
                                   sample(ncell(mask.fp),
                                          10000,
                                          prob = values(mask.fp))))
    
    # convert matrix output to data frame
    bg.xy <- as.data.frame(bg.xy)
    names(bg.xy) <- c("longitude", "latitude")
    # write.csv(x = bg.xy, file = paste0("./Data/BGData/", sp, "_bgxy.csv"), row.names = FALSE)
    saveRDS(object = bg.xy, file = paste0("./Outputs/", sp, "_bg.Rds"))
    
    # Create datasets for training and testing ----
    
    # flickr only test dataset
    test_Flickr <- df[df$source == "Flickr",]
    
    # GBIF only dataframe from geo and random subsampling
    df_bCV <- df[!df$source == "Flickr",]
    
    # blockCV used to generate geo subset based on environmental space
    spdf <- SpatialPointsDataFrame(coords = df_bCV[,2:1],
                                   data = df_bCV,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    spauto.rast <- spatialAutoRange(rasterLayer = envsBgCrop, sampleNumber = 10000,
                                    doParallel = TRUE, border = bgExt)
    block.size <- spauto.rast$range
    
    try(rm(sb))
    try(
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size, # size of the blocks
                         k = 5,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    )
    
    # while loop to change block size if the sb fails to be generated
    while(!exists("sb")){
      block.size <- block.size*0.95
      try(
        sb <- spatialBlock(speciesData = spdf,
                           border = bgExt,
                           rasterLayer = envsBgCrop,
                           theRange = block.size,
                           k = 5,
                           selection = "random",
                           iteration = 500, # find evenly dispersed folds
                           biomod2Format = TRUE,
                           xOffset = 0, 
                           yOffset = 0)
      )
    }
    
    # second while loop that makes sure all groups have occurrences in them
    while(length(unique(sb$foldID)) < 5){
      block.size <- block.size*0.95
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size,
                         k = 5,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    }
    
    df_bCV$group <- sb$foldID
    # make sure bg points are assigned to a 0 group
    bg.grp <- rep(0, dim(bg.xy)[1])
    
    blockcv.map <- sb$plots +
      # geom_point(data = df_bCV, aes(x = longitude, y = latitude, colour = as.factor(group))) +
      scale_fill_scico(palette = "roma", direction = 1) +
      scale_colour_viridis_d(begin = 0, end = 1) +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) +
      theme_bw()
    
    png(filename = paste0(loc.figure, "/", sp, "_BlockCV splits.png"), width = 18,
        height = 12, units = "cm", res = 1200)
    print(blockcv.map)
    dev.off()

    png(filename = paste0(loc.figure, "/", sp, "_BlockCV Autocorrel range.png"), width = 18,
        height = 12, units = "cm", res = 1200)
    print(spauto.rast$plots$barchart)
    dev.off()
    
    # select median group size to be GBIF geo test dataset
    group_test <- names(sort(table(df_bCV$group))[3])
    test_gbif <- df_bCV[df_bCV$group == group_test,]
    
    # and the remaining as the training dataset for GBIF geo
    # create modelling train dataset for only gbif
    train_gbif <- df_bCV[!df_bCV$group == group_test,]
    
    # make sure the groups start at 1 and are consecutive
    if(!group_test == 5){
      train_gbif$group[train_gbif$group == 5] <- 
      unique(c(1:5)[!c(1:5) %in% train_gbif$group])
    }
    
    # subset completely randomly to create the GBIF random test and train datasets
    test_gbifran <- sample_n(df_bCV, dim(test_gbif)[1])
    train_gbifran <- df_bCV[!df_bCV$id %in% test_gbifran$id,]
    
    # re-run blockCV to create groups within the subset training datasets for
    # internal cross validation of the model
    spdf <- SpatialPointsDataFrame(coords = train_gbifran[,2:1],
                                   data = train_gbifran,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    block.size <- spauto.rast$range
    try(rm(sb))
    try(
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size, # size of the blocks
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    )
    
    while(!exists("sb")){
      block.size <- block.size*0.95
      try(
        sb <- spatialBlock(speciesData = spdf,
                           border = bgExt,
                           rasterLayer = envsBgCrop,
                           theRange = block.size,
                           k = 4,
                           selection = "random",
                           iteration = 500, # find evenly dispersed folds
                           biomod2Format = TRUE,
                           xOffset = 0, 
                           yOffset = 0)
      )
    }
    
    while(length(unique(sb$foldID)) < 4){
      block.size <- block.size*0.95
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size,
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    }
    train_gbifran$group <- sb$foldID
    
    # create model train datasets and internal groupings for Flickr supplemented data
    test_Flickr$group <- NA
    train_Flickr <- rbind(train_gbif, test_Flickr)
    
    spdf <- SpatialPointsDataFrame(coords = train_Flickr[,2:1],
                                   data = train_Flickr,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    block.size <- spauto.rast$range
    try(rm(sb))
    try(
      sb <- spatialBlock(speciesData = spdf,
                         # species = "name",
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size, # size of the blocks
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    )
    
    while(!exists("sb")){
      block.size <- block.size*0.95
      try(
        sb <- spatialBlock(speciesData = spdf,
                           # species = "name",
                           border = bgExt,
                           rasterLayer = envsBgCrop,
                           theRange = block.size,
                           k = 4,
                           selection = "random",
                           iteration = 500, # find evenly dispersed folds
                           biomod2Format = TRUE,
                           xOffset = 0, 
                           yOffset = 0)
      )
    }
    
    while(length(unique(sb$foldID)) < 4){
      block.size <- block.size*0.95
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size,
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    }
    train_Flickr$group <- sb$foldID
    
    # repeat the above for GBIF random data now supplemented by Flickr data 
    train_FlickrRan <- rbind(train_gbifran, test_Flickr)
    spdf <- SpatialPointsDataFrame(coords = train_FlickrRan[,2:1],
                                   data = train_FlickrRan,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    block.size <- spauto.rast$range
    try(rm(sb))
    try(
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size, # size of the blocks
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    )
    
    while(!exists("sb")){
      block.size <- block.size*0.95
      try(
        sb <- spatialBlock(speciesData = spdf,
                           border = bgExt,
                           rasterLayer = envsBgCrop,
                           theRange = block.size,
                           k = 4,
                           selection = "random",
                           iteration = 500, # find evenly dispersed folds
                           biomod2Format = TRUE,
                           xOffset = 0,
                           yOffset = 0)
      )
    }
    
    while(length(unique(sb$foldID)) < 4){
      block.size <- block.size*0.95
      sb <- spatialBlock(speciesData = spdf,
                         border = bgExt,
                         rasterLayer = envsBgCrop,
                         theRange = block.size,
                         k = 4,
                         selection = "random",
                         iteration = 500, # find evenly dispersed folds
                         biomod2Format = TRUE,
                         xOffset = 0, 
                         yOffset = 0)
    }
    train_FlickrRan$group <- sb$foldID
    
    
    # Save all training and test datasets ----
    
    # three tests, as flickr only test is run against both GBIF models
    # four training sets, two GBIF only, two Flickr supplemented
    # GBIF only differ in methods used to subset out test datasets
    # (environmentally-independent, pure random)
        
    # test dataset of purely gbif data, environmentally independent
    # (to test against gbif only and gbif+flickr model)
    test_gbif
    saveRDS(object = test_gbif, file = paste0("./Outputs/", sp, "_testGBIF.Rds"))
    # test dataset that is the same size as test_gbif but randomly selected
    test_gbifran
    saveRDS(object = test_gbifran, file = paste0("./Outputs/", sp, "_testGBIFran.Rds"))
    # test dataset of flickr data only (only to be tested against any GBIF trained model)
    test_Flickr
    saveRDS(object = test_Flickr, file = paste0("./Outputs/", sp, "_testFlickr.Rds"))
    # training dataset of gbif only - geographic space
    train_gbif
    saveRDS(object = train_gbif, file = paste0("./Outputs/", sp, "_trainGBIF.Rds"))
    # training dataset of gbif only - completely random
    train_gbifran
    saveRDS(object = train_gbifran, file = paste0("./Outputs/", sp, "_trainGBIFran.Rds"))
    # training dataset of gbif and flickr (near entire dataset)
    train_Flickr
    saveRDS(object = train_Flickr, file = paste0("./Outputs/", sp, "_trainFlickr.Rds"))
    # training dataset of gbif and flickr (near entire dataset)
    train_FlickrRan
    saveRDS(object = train_FlickrRan, file = paste0("./Outputs/", sp, "_trainFlickrRan.Rds"))
    
    # save raster layers for later niche testing
    writeRaster(file = paste0("./Outputs/", sp, "_rasters.grd"), x = envsBgCrop,
                overwrite = TRUE)
    
    # Run Maxent with varying levels of regularization and model types ----
    
    # define the vector of regularization multipliers to test
    rms <- seq(1, 8, 1)
    # iterate model building over all chosen parameter settings
    model.types <- c("L", "LQ")
    # the four training dataset scenarios
    runScens <- c("GBIF",
                  "FlickrSupp",
                  "GBIFran",
                  "FlickrSuppRan")
    
    for(run in runScens){
      
      if(run == "GBIF"){
        loc.model <- "./Outputs/GBIF_mod/"
        dir.create(loc.model)
        train_data <- train_gbif
        
      } else if(run == "FlickrSupp"){
        loc.model <- "./Outputs/FlickrSupp_mod/"
        dir.create(loc.model)
        train_data <- train_Flickr
        
      } else if(run == "GBIFran"){
        loc.model <- "./Outputs/GBIFran_mod/"
        dir.create(loc.model)
        train_data <- train_gbifran
        
      } else if(run == "FlickrSuppRan"){
        loc.model <- "./Outputs/FlickrSuppRan_mod/"
        dir.create(loc.model)
        train_data <- train_FlickrRan
        
      }
      
      occs.xy <- train_data[,2:1]
      occs.grp <- train_data$group
      
      # remove older raster data to free memory
      rm(list = ls()[ls() %in% c("envsBgCrop", "fp.BgCrop", "mask.fp",
                                 "sb", "spauto.rast", "blockcv.map")])
      gc()
      
      start <- Sys.time()
      e <- ENMeval::ENMevaluate(occ = as.matrix(occs.xy), 
                                env = envsBgMsk, 
                                bg.coords = as.matrix(bg.xy),
                                RMvalues = rms, 
                                fc = model.types, 
                                method = 'user', 
                                occ.grp = occs.grp, 
                                bg.grp = bg.grp,
                                bin.output = TRUE, 
                                parallel = FALSE, 
                                numCores = 6)
      end <- Sys.time()
      print(difftime(end, start, units = "mins"))
      
      # unpack the results data frame, the list of models, and the RasterStack of raw predictions
      evalTbl <- e@results
      saveRDS(object = evalTbl, file = paste0(loc.model, "/", sp, "_evalTbl.Rds"))
      
      evalMods <- e@models
      names(evalMods) <- paste0(evalTbl$features, evalTbl$rm)
      save(object = evalMods, file = paste0(loc.model, "/", sp, "_evalMods.Rdata"))
      
      model.list <- names(evalMods)
      write.csv(x = model.list, file = paste0(loc.model, "/", sp, "_modellist.txt"),
                row.names = FALSE)
      
      evalPreds <- e@predictions
      saveRDS(object = evalPreds, file = paste0(loc.model, "/", sp, "_evalPreds.Rds"))
      
    } # end of runScens loop
  } # end of species loop
} # end of random/nonrandom loop
