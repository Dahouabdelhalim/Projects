
# load previous image
load('Script_2_image.RData')

### COMPLETE MAXENT PIPELINE ###

# we will be using the ENMeval library to trial a number of model parameters, and choose the best.
# first, let's load the current climatic conditions
current_climate <- stack(raster('Climate/Current/wc2.1_2.5m_bio_1.tif'),  # annual mean temperature
                         raster('Climate/Current/wc2.1_2.5m_bio_2.tif'),  # mean diurnal range
                         raster('Climate/Current/wc2.1_2.5m_bio_12.tif'), # annual precipitation
                         raster('Climate/Current/wc2.1_2.5m_bio_16.tif'), # precipitation in the wet quarter
                         raster('Climate/Current/wc2.1_2.5m_bio_17.tif'), # precipitation in the dry quarter
                         raster('Climate/wc2.1_2.5m_elev.tif'))           # elevation

# make a vector of climate models we will investigate for future time periods
climate_models <- list.files(path = 'Climate/Future',
                             pattern = 'ssp126_2021-2040.tif')
climate_models <- gsub('wc2.1_2.5m_bioc_',
                       '',
                       climate_models)
climate_models <- gsub('_ssp126_2021-2040.tif',
                       '',
                       climate_models)
# make a vector of time scenarios
time <- c('2021-2040',
          '2041-2060',
          '2061-2080',
          '2081-2100')
# and emissions scenarios
emissions <- c('ssp126', 'ssp370')
# read in elevation for later
elevation <- raster('Climate/wc2.1_2.5m_elev.tif')

# For the pipeline, we need to define the occurrences, the background, and the study area
# we will do this for the overseas species first, because they will be slightly different
for(i in 1:length(overseas_frogs)) {
  tryCatch({ 
    print(paste('Starting',
                i,
                'of',
                length(overseas_frogs),
                'at',
                Sys.time()))
    # first, define occurrences
    occurrences <- subset(australia_frog_occurrences, Species == overseas_frogs[i])
    occs_extra <- subset(occs_for_overseas_bias[[i]], Species == overseas_frogs[i])
    final_occs <- rbind(occurrences[,1:3],
                        occs_extra[,1:3])
    # then, define the study area
    # australian territory first - Frogs of Australia App
    aus_territory <- readOGR(dsn = 'Shapefiles',
                             layer = gsub(' ', '_', overseas_frogs[i]))
    # buffer the centroid by ~10% or 1 degree, whatever is bigger
    {
      poly_centroid <- gCentroid(aus_territory)
      bbox_poly <- as.data.frame(bbox(aus_territory))
      point_1 <- c(bbox_poly[1,1],
                   bbox_poly[2,2])
      point_2 <- c(bbox_poly[1,2],
                   bbox_poly[2,2])
      point_3 <- c(bbox_poly[1,1],
                   bbox_poly[2,1])
      point_4 <- c(bbox_poly[1,2],
                   bbox_poly[2,1])
      point_5 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,2])
      point_6 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,1])
      point_7 <- c(bbox_poly[1,1],
                   poly_centroid@coords[1,2])
      point_8 <- c(bbox_poly[1,2],
                   poly_centroid@coords[1,2])
      dis_1 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_1, r = 6378.137)
      dis_2 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_2, r = 6378.137)
      dis_3 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_3, r = 6378.137)
      dis_4 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_4, r = 6378.137)
      dis_5 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_5, r = 6378.137)
      dis_6 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_6, r = 6378.137)
      dis_7 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_7, r = 6378.137)
      dis_8 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_8, r = 6378.137)
      mean_dis <- mean(dis_1,
                       dis_2,
                       dis_3,
                       dis_4,
                       dis_5,
                       dis_6,
                       dis_7,
                       dis_8)
      mean_dis_10_percent <- mean_dis * 0.1 # mean distance from centroid to 8 points around bbox, in km
      # convert to degrees
      mean_dis_degree_aus <- mean_dis_10_percent/distHaversine(c(1, poly_centroid@coords[2]),
                                                               c(2, poly_centroid@coords[2]),
                                                               r = 6378.137) 
      # if less than 0.5 degrees, buffer by this distance, other wise buffer by 0.5 degrees
      poly_buf <- gBuffer(aus_territory, width = ifelse(mean_dis_degree_aus < 0.5,
                                                        mean_dis_degree_aus,
                                                        0.5))
      # crop to aus
      poly_buf_aus_territory <- crop(poly_buf, aus)
      
    }
    # now, need to repeat this process for the overseas territory  aus_territory <- readOGR(dsn = 'Shapefiles',
    # buffer the centroid by ~10% or 1 degree, whatever is bigger
    {
      # get the shapefile
      frog_shapefile <- iucn_shapefile[iucn_shapefile@data$BINOMIAL==overseas_frogs[i],]
      # keep only the part outside of australia
      outside <- gDifference(frog_shapefile, aus)
      outside <- crop(outside, world)
      poly_centroid <- gCentroid(outside)
      bbox_poly <- as.data.frame(bbox(outside))
      point_1 <- c(bbox_poly[1,1],
                   bbox_poly[2,2])
      point_2 <- c(bbox_poly[1,2],
                   bbox_poly[2,2])
      point_3 <- c(bbox_poly[1,1],
                   bbox_poly[2,1])
      point_4 <- c(bbox_poly[1,2],
                   bbox_poly[2,1])
      point_5 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,2])
      point_6 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,1])
      point_7 <- c(bbox_poly[1,1],
                   poly_centroid@coords[1,2])
      point_8 <- c(bbox_poly[1,2],
                   poly_centroid@coords[1,2])
      dis_1 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_1, r = 6378.137)
      dis_2 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_2, r = 6378.137)
      dis_3 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_3, r = 6378.137)
      dis_4 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_4, r = 6378.137)
      dis_5 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_5, r = 6378.137)
      dis_6 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_6, r = 6378.137)
      dis_7 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_7, r = 6378.137)
      dis_8 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_8, r = 6378.137)
      mean_dis <- mean(dis_1,
                       dis_2,
                       dis_3,
                       dis_4,
                       dis_5,
                       dis_6,
                       dis_7,
                       dis_8)
      mean_dis_10_percent <- mean_dis * 0.1 # mean distance from centroid to 8 points around bbox, in km
      # convert to degrees
      mean_dis_degree <- mean_dis_10_percent/distHaversine(c(1, poly_centroid@coords[2]),
                                                           c(2, poly_centroid@coords[2]),
                                                           r = 6378.137) 
      # if less than 0.5 degrees, buffer by this distance, other wise buffer by 0.5 degrees
      poly_buf <- gBuffer(outside, width = ifelse(mean_dis_degree < 0.5,
                                                  mean_dis_degree,
                                                  0.5))
      # crop to aus
      poly_buf_not_au <- crop(poly_buf, world)
    }
    # now, need to crop occurrences outside the buffered zone and thin to 1 occ/2.5m cell
    coordinates(final_occs) <- ~Lon+Lat
    buffer <- bind(poly_buf_aus_territory, poly_buf_not_aus)
    final_occs <- crop(final_occs, buffer)  
    cells <- cellFromXY(current_climate, final_occs)
    dups <- duplicated(cells)
    final_occs_thin <- as.data.frame(final_occs[!dups,])
    # buffer background by same distance again to make study area
    study_area <- gBuffer(buffer, width = mean_dis_degree_aus)
    study_area <- crop(study_area, world)
    # crop predictors to this area
    predictors_study_area <- crop(current_climate, study_area)
    predictors_study_area <- mask(predictors_study_area, study_area)
    # crop to 5% of background if necessary
    cell_background <- ncell(predictors_study_area)
    if(nrow(final_occs_thin) > 0.05 * cell_background) {
      final_occs_thin <- final_occs_thin[sample(nrow(final_occs_thin), size = 0.05*cell_background),]
    }
    # sample the background points (20% of available cells)  
    # get background occurrences
    # first, crop the bias layer, and adjust scaled values to above 0 for sampling
    bias_crop <- crop(overseas_bias[[i]], predictors_study_area)
    values(bias_crop) <- values(bias_crop) + abs(min(values(bias_crop), na.rm = T))
    values(bias_crop)[is.na(values(bias_crop))] <- 0
    # select more than needed
    n.bg <- 0.25*cell_background
    biased.bg <- as.data.frame(xyFromCell(bias_crop, sample(ncell(bias_crop), n.bg, prob = values(bias_crop), replace = TRUE)))
    # thin
    coordinates(biased.bg) <- ~x+y
    cells <- cellFromXY(bias_crop, biased.bg)
    dups <- duplicated(cells)
    biased.bg <- as.data.frame(biased.bg[!dups,])
    # finally, select 20% of background (if needed)
    if(nrow(biased.bg) > 0.2*cell_background) {
      biased.bg <- biased.bg[sample(nrow(biased.bg), size = 0.2*cell_background),]
    }
    ### READY TO RUN ENM EVAL
    occs_for_test <- final_occs_thin[,c(3,2)]
    colnames(occs_for_test) <- c('x','y')
    enm_eval_mod <- ENMevaluate(occs = occs_for_test,
                                envs = predictors_study_area,
                                bg = biased.bg,
                                tune.args = list(fc = c('L', 'LQ', 'H', 'LQH'),
                                                 rm = c(0.5,1,1.5,2)),
                                partitions = 'block',
                                algorithm = 'maxent.jar',
                                parallel = TRUE,
                                numCores = 8,
                                parallelType = 'doParallel',
                                quiet = TRUE,
                                updateProgress = TRUE)
    
    res <- eval.results(enm_eval_mod)
    # select the best model as per ENMEval 2.0 Vignette
    opt.aicc <- res %>% filter(delta.AICc == 0)
    opt.seq <- res %>% 
      filter(or.10p.avg == min(or.10p.avg)) %>% 
      filter(auc.val.avg == max(auc.val.avg))
    mod.seq <- eval.models(enm_eval_mod)[[opt.seq$tune.args]]
    # get maxSSS
    maxsss <- as.numeric(mod.seq@results['Maximum.training.sensitivity.plus.specificity.area',])
    # add info to summary
    opt.seq$occurrences_in_model <- nrow(final_occs_thin)
    opt.seq$maxsss <- maxsss
    opt.seq$species <- overseas_frogs[i]
    write.csv(opt.seq,
              paste('Statistical_results/',
                    gsub(' ', '_', overseas_frogs[i]),
                    '_enmevaluate_summary.csv',
                    sep = ''))
    # run a null model with the selected parameters
    mod.null <- ENMnulls(enm_eval_mod, 
                         mod.settings = list(fc = as.character(opt.seq$fc), rm = as.numeric(as.character(opt.seq$rm))), 
                         no.iter = 100,
                         quiet = TRUE)
    stats <- null.emp.results(mod.null)
    # write out stats results
    write.csv(stats,
              paste('Statistical_results/',
                    gsub(' ', '_', overseas_frogs[i]),
                    '_null_model_stats.csv',
                    sep = ''))
    # predict model to current environment - note that we are only interested in Australia now
    predictors_aus <- mask(predictors_study_area, aus)
    current_projection <- predict(mod.seq,
                                  predictors_aus,
                                  predictors_aus)
    current_projection_extend <- extend(current_projection, aus)
    writeRaster(current_projection_extend,
                overwrite =T,
                filename = paste('Result_rasters/',
                                 gsub(' ', '_', overseas_frogs[i]), '_',
                                 'current_projection.tif',
                                 sep = ''))

    # make current projection binary
    binary <- current_projection_extend
    binary[binary < maxsss] <- NA
    binary[binary >= maxsss] <- 1
    writeRaster(current_projection_extend,
                overwrite =T,
                filename = paste('Result_rasters/',
                                 gsub(' ', '_', overseas_frogs[i]), '_',
                                 'current_projection_binary.tif',
                                 sep = ''))
    # make into polygon for start
    binary_border <- rasterToPolygons(binary, dissolve = TRUE)
    # calculate current statistics
    current_area_raster <- area(binary, na.rm = TRUE, weights = FALSE)
    current_area_raster <- current_area_raster[!is.na(current_area_raster)]
    current_area <- length(current_area_raster) * median(current_area_raster)
    centroid <- gCentroid(binary_border)
    
    # make a dataframe for the results
    results <- data.frame(Species = overseas_frogs[i],
                          Area = current_area,
                          Model = 'Current',
                          Emissions = 'Current',
                          Year = 'Present',
                          Centroid_X = centroid@coords[1],
                          Centroid_Y = centroid@coords[2])
    
    # get a minimum elevation the species is at now
    min_elevation_raster <- crop(elevation, binary_border)
    min_elevation_raster <- mask(elevation, binary_border)
    minimum_elevation <- min(values(min_elevation_raster), na.rm = TRUE)
    elev <- crop(min_elevation_raster, gBuffer(binary_border, width = 0.5))
    elev[elev<minimum_elevation] <- NA
    elev_border <- rasterToPolygons(elev, dissolve = T)
    # now, make predictions for the future, under a range of climate scenarios and models
    for(j in 1:length(climate_models)) {
      # keep track of time
      print(paste('Starting Model',
                  j,
                  '-',
                  climate_models[j],
                  'at',
                  Sys.time()))
      for(k in 1:length(emissions)) {
        for(l in 1:length(time)) {
          # read in relevant predictors
          predictors <- stack(raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 1),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 2),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 12),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 16),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 17),
                              raster('Climate/wc2.1_2.5m_elev.tif'))
          names(predictors) <- names(predictors_aus)
          # get an extent estimate - allow the species to dispers 0.1 degrees, but no decrease in elevation
          if(l == 1) {
            border <- gBuffer(binary_border, width = 0.1)
            border <- crop(border, elev_border)
          } else {
            border <- gBuffer(border, width = 0.1)
            border <- crop(border, elev_border)
          }
          # make predictions
          prediction <- predict(mod.seq,
                                predictors,
                                ext = border)
          prediction <- mask(prediction, aus)
          prediction_aus <- extend(prediction, aus)
          prediction_aus_binary <- prediction_aus
          prediction_aus_binary[prediction_aus_binary<maxsss] <- NA
          prediction_aus_binary[prediction_aus_binary>=maxsss] <- 1
          # calculate current statistics
          prediction_aus_binary_area <- area(prediction_aus_binary, na.rm = TRUE, weights = FALSE)
          prediction_aus_binary_area <- prediction_aus_binary_area[!is.na(prediction_aus_binary_area)]
          prediction_aus_binary_area <- length(prediction_aus_binary_area) * median(prediction_aus_binary_area)
          border <- rasterToPolygons(prediction_aus_binary, dissolve = T)
          centroid <- gCentroid(border)
          
          # make a dataframe for the results
          results_ <- data.frame(Species = overseas_frogs[i],
                                 Area = prediction_aus_binary_area,
                                 Model = climate_models[j],
                                 Emissions = emissions[k],
                                 Year = time[l],
                                 Centroid_X = centroid@coords[1],
                                 Centroid_Y = centroid@coords[2])
          results <- rbind(results, results_)
          write.csv(results,
                    paste('Statistical_results/',
                          gsub(' ', '_', overseas_frogs[i]),
                          '_prediction_statistics.csv',
                          sep = ''))
          # store results
          writeRaster(prediction_aus,
                      overwrite =T,
                      filename = paste('Result_rasters/',
                                       gsub(' ', '_', overseas_frogs[i]), '_',
                                       climate_models[j], '_',
                                       emissions[k], '_',
                                       time[l],
                                       '.tif',
                                       sep = ''))
          writeRaster(prediction_aus_binary,
                      overwrite =T,
                      filename = paste('Result_rasters/',
                                       gsub(' ', '_', overseas_frogs[i]), '_',
                                       climate_models[j], '_',
                                       emissions[k], '_',
                                       time[l],
                                       '.tif',
                                       sep = ''))
        }
      }
    }
  }, error=function(e){cat("ERROR :", conditionMessage(e),"\\n")})
}
# now, for the remaining species
remaining_frogs <- species_information$Species
remaining_frogs <- remaining_frogs[!remaining_frogs %in% overseas_frogs]
for(i in 1:length(remaining_frogs)) {
  tryCatch({ 
    print(paste('Starting',
                i,
                'of',
                length(remaining_frogs),
                'at',
                Sys.time()))
    # first, define occurrences
    final_occs <- subset(australia_frog_occurrences, Species == remaining_frogs[i])
    
    # then, define the study area
    # australian territory first - Frogs of Australia App
    aus_territory <- readOGR(dsn = 'Shapefiles',
                             layer = gsub(' ', '_', remaining_frogs[i]))
    # buffer the centroid by ~10% or 1 degree, whatever is bigger
    {
      poly_centroid <- gCentroid(aus_territory)
      bbox_poly <- as.data.frame(bbox(aus_territory))
      point_1 <- c(bbox_poly[1,1],
                   bbox_poly[2,2])
      point_2 <- c(bbox_poly[1,2],
                   bbox_poly[2,2])
      point_3 <- c(bbox_poly[1,1],
                   bbox_poly[2,1])
      point_4 <- c(bbox_poly[1,2],
                   bbox_poly[2,1])
      point_5 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,2])
      point_6 <- c(poly_centroid@coords[1,1],
                   bbox_poly[2,1])
      point_7 <- c(bbox_poly[1,1],
                   poly_centroid@coords[1,2])
      point_8 <- c(bbox_poly[1,2],
                   poly_centroid@coords[1,2])
      dis_1 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_1, r = 6378.137)
      dis_2 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_2, r = 6378.137)
      dis_3 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_3, r = 6378.137)
      dis_4 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_4, r = 6378.137)
      dis_5 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_5, r = 6378.137)
      dis_6 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_6, r = 6378.137)
      dis_7 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_7, r = 6378.137)
      dis_8 <- distHaversine(c(poly_centroid@coords[1,1],
                               poly_centroid@coords[1,2]), 
                             point_8, r = 6378.137)
      mean_dis <- mean(dis_1,
                       dis_2,
                       dis_3,
                       dis_4,
                       dis_5,
                       dis_6,
                       dis_7,
                       dis_8)
      mean_dis_10_percent <- mean_dis * 0.1 # mean distance from centroid to 8 points around bbox, in km
      # convert to degrees
      mean_dis_degree_aus <- mean_dis_10_percent/distHaversine(c(1, poly_centroid@coords[2]),
                                                               c(2, poly_centroid@coords[2]),
                                                               r = 6378.137) 
      # if less than 0.5 degrees, buffer by this distance, other wise buffer by 0.5 degrees
      poly_buf <- gBuffer(aus_territory, width = ifelse(mean_dis_degree_aus < 0.5,
                                                        mean_dis_degree_aus,
                                                        0.5))
      # crop to aus
      buffer <- crop(poly_buf, aus)
      
    }
    
    # now, need to crop occurrences outside the buffered zone and thin to 1 occ/2.5m cell
    coordinates(final_occs) <- ~Lon+Lat
    final_occs <- crop(final_occs, buffer)  
    cells <- cellFromXY(current_climate, final_occs)
    dups <- duplicated(cells)
    final_occs_thin <- as.data.frame(final_occs[!dups,])
    # buffer background by same distance again to make study area
    study_area <- gBuffer(buffer, width = mean_dis_degree_aus)
    study_area <- crop(study_area, aus)
    # crop predictors to this area
    predictors_study_area <- crop(current_climate, study_area)
    predictors_study_area <- mask(predictors_study_area, study_area)
    # crop to 5% of background if necessary
    cell_background <- ncell(predictors_study_area)
    if(nrow(final_occs_thin) > 0.05 * cell_background) {
      final_occs_thin <- final_occs_thin[sample(nrow(final_occs_thin), size = 0.05*cell_background),]
    }
    # sample the background points (20% of available cells)  
    # get background occurrences
    # first, crop the bias layer, and adjust scaled values to above 0 for sampling
    bias_crop <- crop(bias_layer_australia, predictors_study_area)
    values(bias_crop) <- values(bias_crop) + abs(min(values(bias_crop), na.rm = T))
    values(bias_crop)[is.na(values(bias_crop))] <- 0
    # select more than needed
    n.bg <- 0.25*cell_background
    biased.bg <- as.data.frame(xyFromCell(bias_crop, sample(ncell(bias_crop), n.bg, prob = values(bias_crop), replace = TRUE)))
    # thin
    coordinates(biased.bg) <- ~x+y
    cells <- cellFromXY(bias_crop, biased.bg)
    dups <- duplicated(cells)
    biased.bg <- as.data.frame(biased.bg[!dups,])
    # finally, select 20% of background (if needed)
    if(nrow(biased.bg) > 0.2*cell_background) {
      biased.bg <- biased.bg[sample(nrow(biased.bg), size = 0.2*cell_background),]
    }
    ### READY TO RUN ENM EVAL
    occs_for_test <- final_occs_thin[,c(3,2)]
    colnames(occs_for_test) <- c('x','y')
    enm_eval_mod <- ENMevaluate(occs = occs_for_test,
                                envs = predictors_study_area,
                                bg = biased.bg,
                                tune.args = list(fc = c('L', 'LQ', 'H', 'LQH'),
                                                 rm = c(0.5,1,1.5,2)),
                                partitions = 'block',
                                algorithm = 'maxent.jar',
                                parallel = TRUE,
                                numCores = 8,
                                parallelType = 'doParallel',
                                quiet = TRUE,
                                updateProgress = TRUE)
    res <- eval.results(enm_eval_mod)
    # select the best model as per ENMEval 2.0 Vignette
    opt.aicc <- res %>% filter(delta.AICc == 0)
    opt.seq <- res %>% 
      filter(or.10p.avg == min(or.10p.avg)) %>% 
      filter(auc.val.avg == max(auc.val.avg))
    mod.seq <- eval.models(enm_eval_mod)[[opt.seq$tune.args]]
    # get maxSSS
    maxsss <- as.numeric(mod.seq@results['Maximum.training.sensitivity.plus.specificity.training.omission',])
    opt.seq$occurrences_in_model <- nrow(final_occs_thin)
    opt.seq$maxsss <- maxsss
    opt.seq$species <- remaining_frogs[i]
    write.csv(opt.seq,
              paste('Statistical_results/',
                    gsub(' ', '_', remaining_frogs[i]),
                    '_enmevaluate_summary.csv',
                    sep = ''))
    # run a null model with the selected parameters
    mod.null <- ENMnulls(enm_eval_mod, 
                         mod.settings = list(fc = as.character(opt.seq$fc), rm = as.numeric(as.character(opt.seq$rm))), 
                         no.iter = 100,
                         quiet = TRUE)
    stats <- null.emp.results(mod.null)
    # write out stats results
    write.csv(stats,
              paste('Statistical_results/',
                    gsub(' ', '_', remaining_frogs[i]),
                    '_null_model_stats.csv',
                    sep = ''))
    # predict model to current environment - note that we are only interested in Australia now
    predictors_aus <- mask(predictors_study_area, aus)
    current_projection <- predict(mod.seq,
                                  predictors_aus,
                                  predictors_aus)
    current_projection_extend <- extend(current_projection, aus)
    writeRaster(current_projection_extend,
                overwrite =T,
                filename = paste('Result_rasters/',
                                 gsub(' ', '_', remaining_frogs[i]), '_',
                                 'current_projection.tif',
                                 sep = ''))
    # make current projection binary
    binary <- current_projection_extend
    binary[binary < maxsss] <- NA
    binary[binary >= maxsss] <- 1
    writeRaster(current_projection_extend,
                overwrite =T,
                filename = paste('Result_rasters/',
                                 gsub(' ', '_', remaining_frogs[i]), '_',
                                 'current_projection_binary.tif',
                                 sep = ''))
    # make into polygon for start
    binary_border <- rasterToPolygons(binary, dissolve = TRUE)
    # calculate current statistics
    current_area_raster <- area(binary, na.rm = TRUE, weights = FALSE)
    current_area_raster <- current_area_raster[!is.na(current_area_raster)]
    current_area <- length(current_area_raster) * median(current_area_raster)
    centroid <- gCentroid(binary_border)
    
    # make a dataframe for the results
    results <- data.frame(Species = remaining_frogs[i],
                          Area = current_area,
                          Model = 'Current',
                          Emissions = 'Current',
                          Year = 'Present',
                          Centroid_X = centroid@coords[1],
                          Centroid_Y = centroid@coords[2])
    
    # get a minimum elevation the species is at now
    min_elevation_raster <- crop(elevation, binary_border)
    min_elevation_raster <- mask(elevation, binary_border)
    minimum_elevation <- min(values(min_elevation_raster), na.rm = TRUE)
    elev <- crop(min_elevation_raster, gBuffer(binary_border, width = 0.5))
    elev[elev<minimum_elevation] <- NA
    elev_border <- rasterToPolygons(elev, dissolve = T)
    # now, make predictions for the future, under a range of climate scenarios and models
    for(j in 1:length(climate_models)) {
      # keep track of time
      print(paste('Starting Model',
                  j,
                  '-',
                  climate_models[j],
                  'at',
                  Sys.time()))
      for(k in 1:length(emissions)) {
        for(l in 1:length(time)) {
          # read in relevant predictors
          predictors <- stack(raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 1),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 2),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 12),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 16),
                              raster(paste('Climate/Future/wc2.1_2.5m_bioc_', 
                                           climate_models[j], '_',
                                           emissions[k], '_',
                                           time[l],
                                           '.tif',
                                           sep = ''),
                                     band = 17),
                              raster('Climate/wc2.1_2.5m_elev.tif'))
          names(predictors) <- names(predictors_aus)
          # get an extent estimate - allow the species to dispers 0.1 degrees, but no decrease in elevation
          if(l == 1) {
            border <- gBuffer(binary_border, width = 0.1)
            border <- crop(border, elev_border)
          } else {
            border <- gBuffer(border, width = 0.1)
            border <- crop(border, elev_border)
          }
          # make predictions
          prediction <- predict(mod.seq,
                                predictors,
                                ext = border)
          prediction <- mask(prediction, aus)
          prediction_aus <- extend(prediction, aus)
          prediction_aus_binary <- prediction_aus
          prediction_aus_binary[prediction_aus_binary<maxsss] <- NA
          prediction_aus_binary[prediction_aus_binary>=maxsss] <- 1
          # calculate current statistics
          prediction_aus_binary_area <- area(prediction_aus_binary, na.rm = TRUE, weights = FALSE)
          prediction_aus_binary_area <- prediction_aus_binary_area[!is.na(prediction_aus_binary_area)]
          prediction_aus_binary_area <- length(prediction_aus_binary_area) * median(prediction_aus_binary_area)
          border <- rasterToPolygons(prediction_aus_binary, dissolve = T)
          centroid <- gCentroid(border)
          
          # make a dataframe for the results
          results_ <- data.frame(Species = remaining_frogs[i],
                                 Area = prediction_aus_binary_area,
                                 Model = climate_models[j],
                                 Emissions = emissions[k],
                                 Year = time[l],
                                 Centroid_X = centroid@coords[1],
                                 Centroid_Y = centroid@coords[2])
          results <- rbind(results, results_)
          write.csv(results,
                    paste('Statistical_results/',
                          gsub(' ', '_', remaining_frogs[i]),
                          '_prediction_statistics.csv',
                          sep = ''))
          # store results
          writeRaster(prediction_aus,
                      overwrite =T,
                      filename = paste('Result_rasters/',
                                       gsub(' ', '_', remaining_frogs[i]), '_',
                                       climate_models[j], '_',
                                       emissions[k], '_',
                                       time[l],
                                       '.tif',
                                       sep = ''))
          writeRaster(prediction_aus_binary,
                      overwrite =T,
                      filename = paste('Result_rasters/',
                                       gsub(' ', '_', remaining_frogs[i]), '_',
                                       climate_models[j], '_',
                                       emissions[k], '_',
                                       time[l],
                                       '.tif',
                                       sep = ''))
        }
      }
    }
  }, error=function(e){cat("ERROR :", conditionMessage(e),"\\n")})
}

# save image
save.image(file = 'Script_3_image.RData')