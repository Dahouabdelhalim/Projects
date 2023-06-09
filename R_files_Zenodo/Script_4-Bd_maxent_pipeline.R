# load image
load('Script_3_image.RData')

# load in Bd occurrence data
records <- read.csv('Bd_records_simple.csv')

# thin
coordinates(records) <- ~Lon+Lat
cells <- cellFromXY(current_climate, records)
dups <- duplicated(cells)
final_occs_thin <- records[!dups,]

# create a bias layer
# turn into a raster
occs <- rasterize(final_occs_thin, current_climate, 1)
# calculate density raster
presences <- which(values(occs) == 1)
pres.locs <- coordinates(occs)[presences, ]
dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occs), ncol(occs)))
dens.ras <- raster(dens)
# number of points to select
cell_background <- 100000
biased.bg <- as.data.frame(xyFromCell(dens.ras, sample(ncell(dens.ras), n.bg, prob = values(dens.ras), replace = TRUE)))
# thin
coordinates(biased.bg) <- ~x+y
cells <- cellFromXY(dens.ras, biased.bg)
dups <- duplicated(cells)
biased.bg <- as.data.frame(biased.bg[!dups,])

### READY TO RUN ENM EVAL
occs_for_test <- final_occs_thin
colnames(occs_for_test) <- c('x','y')
enm_eval_mod <- ENMevaluate(occs = occs_for_test,
                            envs = current_climate,
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
opt.seq$species <- 'Bd'
write.csv(opt.seq,
          paste('Statistical_results/',
                gsub(' ', '_', 'Bd'),
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
                gsub(' ', '_', 'Bd'),
                '_null_model_stats.csv',
                sep = ''))
# predict model to current environment - note that we are only interested in Australia now
predictors_aus <- mask(current_climate, aus)
current_projection <- predict(mod.seq,
                              predictors_aus,
                              predictors_aus)
current_projection_extend <- extend(current_projection, aus)
writeRaster(current_projection_extend,
            overwrite =T,
            filename = paste('Result_rasters/',
                             gsub(' ', '_', 'Bd'), '_',
                             'current_projection.tif',
                             sep = ''))

# make current projection binary
binary <- current_projection_extend
binary[binary < maxsss] <- NA
binary[binary >= maxsss] <- 1
writeRaster(current_projection_extend,
            overwrite =T,
            filename = paste('Result_rasters/',
                             gsub(' ', '_', 'Bd'), '_',
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
results <- data.frame(Species = 'Bd',
                      Area = current_area,
                      Model = 'Current',
                      Emissions = 'Current',
                      Year = 'Present',
                      Centroid_X = centroid@coords[1],
                      Centroid_Y = centroid@coords[2])

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

      # make predictions
      prediction <- predict(mod.seq,
                            predictors,
                            ext = predictors)
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
      results_ <- data.frame(Species = 'Bd',
                             Area = prediction_aus_binary_area,
                             Model = climate_models[j],
                             Emissions = emissions[k],
                             Year = time[l],
                             Centroid_X = centroid@coords[1],
                             Centroid_Y = centroid@coords[2])
      results <- rbind(results, results_)
      write.csv(results,
                paste('Statistical_results/',
                      gsub(' ', '_', 'Bd'),
                      '_prediction_statistics.csv',
                      sep = ''))
      # store results
      writeRaster(prediction_aus,
                  overwrite =T,
                  filename = paste('Result_rasters/',
                                   gsub(' ', '_', 'Bd'), '_',
                                   climate_models[j], '_',
                                   emissions[k], '_',
                                   time[l],
                                   '.tif',
                                   sep = ''))
      writeRaster(prediction_aus_binary,
                  overwrite =T,
                  filename = paste('Result_rasters/',
                                   gsub(' ', '_', 'Bd'), '_',
                                   climate_models[j], '_',
                                   emissions[k], '_',
                                   time[l],
                                   '.tif',
                                   sep = ''))
    }
  }
}

save.image(file = 'Script_4_image.RData')
