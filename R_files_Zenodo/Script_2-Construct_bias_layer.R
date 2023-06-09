
# load previous image
load('Script_1_image.RData')

### CREATE A BIAS LAYER FOR MAXENT ###

## Will need to create a general bias layer for all of Australia, as well as for each of the overseas territories for
## species occurring overseas
# read in a 2.5m raster (resolution used) for simplicity
ras <- raster('Climate/Current/wc2.1_2.5m_bio_1.tif')
# convert to spatial
coordinates(australia_frog_occurrences) <- ~Lon+Lat
# turn into a raster
occs <- rasterize(australia_frog_occurrences, ras, 1)
# crop to Australia
occs <- crop(occs, aus)
# calculate density raster
presences <- which(values(occs) == 1)
pres.locs <- coordinates(occs)[presences, ]
dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occs), ncol(occs)))
dens.ras <- raster(dens)
# mask to Australia for ease
bias_layer_australia <- mask(dens.ras, aus)

# now, need to do the same for each of the overseas territories as well. 
# calculate for the overseas territory. Then Australia to match.
overseas_bias <- list()
for(i in 1:length(overseas_frogs)) {
  # isolate occs
  occs <- occs_for_overseas_bias[[i]]
  # convert to spatial
  coordinates(occs) <- ~Lon+Lat
  # get the shapefile
  frog_shapefile <- iucn_shapefile[iucn_shapefile@data$BINOMIAL==overseas_frogs[i],]
  # buffer
  frog_shapefile <- gBuffer(frog_shapefile, width = 1)
  # keep only the part outside of australia
  outside <- gDifference(frog_shapefile, aus)
  outside <- crop(outside, world)
  # crop to occurrences outside here
  occs <- crop(occs, outside)
  # turn into a raster
  occs <- rasterize(occs, ras, 1)
  # crop to buffer zone
  occs <- crop(occs, frog_shapefile)
  # calculate density raster
  presences <- which(values(occs) == 1)
  pres.locs <- coordinates(occs)[presences, ]
  dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occs), ncol(occs)))
  dens.ras <- raster(dens)
  # mask to Australia for ease
  bias_layer_frog <- mask(dens.ras, outside)
  # scale and connect with Australia
  aus_scaled <- scale(bias_layer_australia)
  bias_frog_scaled <- scale(bias_layer_frog)
  # get bbox to crop to speed up
  bbox_aus <- extent(aus_scaled)
  bbox_frog <- extent(bias_frog_scaled)
  # extend
  bbox <- c(min(bbox_aus@xmin,
                bbox_frog@xmin),
            max(bbox_aus@xmax,
                bbox_frog@xmax),
            min(bbox_aus@ymin,
                bbox_frog@ymin),
            max(bbox_aus@ymax,
                bbox_frog@ymax))
  aus_scaled_ <- extend(aus_scaled, bbox)
  bias_frog_scaled_ <- extend(bias_frog_scaled, bbox)
  aus_scaled_ <- resample(aus_scaled_, bias_frog_scaled_, method = 'bilinear')
  # merge and add to list
  layer <- merge(aus_scaled_, bias_frog_scaled_)
  overseas_bias[[i]] <- layer
}
# put occs back to a data frame
australia_frog_occurrences <- as.data.frame(australia_frog_occurrences)

# save image
save.image(file = 'Script_2_image.RData')
