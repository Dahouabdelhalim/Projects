## PEI data set
## He, Juanxia (AAFC/AAC) <juanxia.he@canada.ca>

library(rgdal)
library(maptools)
library(sf)
library(XML)
library(raster)

som = readOGR("one_example_of_training_data_from_random_sampling_DSS_polygons/SOM_PEI.shp", "SOM_PEI")
## 624 points
## SOM is the soil organic matter (unit: percent, g/100g )  in topsoil (from the soil surface to a 23 cm depth) according to the attached paper. No bulk density values are available in the dataset.
npdb = readOGR("validationPoints/PEI_NPDB_Judith.shp", "PEI_NPDB_Judith")
summary(as.factor(npdb$SOIL_TYPE))
## 672

r30m = raster("covariates/DEM.tif")
proj4string(r30m) = proj4string(som)
# 4131, 6242

## GlobalForestWatch data ----
system(paste0('gdalwarp /mnt/DATA/GlobalForestChange2000-2014/first.vrt /data/Canada/tmp/DataForTom/covariates/PEI_Landsat_2000_20m.tif -t_srs \\"', proj4string(r30m), '\\" -co \\"BIGTIFF=YES\\" -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -co \\"COMPRESS=DEFLATE\\" -tr ', res(r30m)[1], ' ', res(r30m)[1], ' -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))
system(paste0('gdalwarp /mnt/DATA/GlobalForestChange2000-2014/last.vrt /data/Canada/tmp/DataForTom/covariates/PEI_Landsat_2018_20m.tif -t_srs \\"', proj4string(r30m), '\\" -co \\"BIGTIFF=YES\\" -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -co \\"COMPRESS=DEFLATE\\" -tr ', res(r30m)[1], ' ', res(r30m)[1], ' -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))
system(paste0('gdalwarp /mnt/DATA/GlobalForestChange2000-2014/treecover2000.vrt /data/Canada/tmp/DataForTom/covariates/PEI_treecover2000_20m.tif -t_srs \\"', proj4string(r30m), '\\" -co \\"BIGTIFF=YES\\" -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -co \\"COMPRESS=DEFLATE\\" -tr ', res(r30m)[1], ' ', res(r30m)[1], ' -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))

## soil map:
library(foreign)
geo = "/data/Canada/tmp/DataForTom/DSS_soil_polygons/Soil_region_Original_adding_group_ID_2018.dbf"
shpF.db <- read.dbf(geo)
str(shpF.db)
## 25 levels
summary(shpF.db$soilname)
shpF.db$TYPE_INT = as.integer(as.factor(shpF.db$soilname))
geo_leg = data.frame(Value=1:length(levels(shpF.db$soilname)), Classes=levels(shpF.db$soilname))
write.dbf(shpF.db, "/data/Canada/tmp/DataForTom/DSS_soil_polygons/Soil_region_Original_adding_group_ID_2018.dbf")
cellsize = 30
te = as.vector(extent(r30m))[c(1,3,2,4)]
system(paste0('saga_cmd -c=64 grid_gridding 0 -INPUT \\"/data/Canada/tmp/DataForTom/DSS_soil_polygons/Soil_region_Original_adding_group_ID_2018.shp\\" -FIELD \\"TYPE_INT\\" -GRID \\"/data/Canada/tmp/DataForTom/DSS_soil_polygons_30m.sgrd\\" -GRID_TYPE 2 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', te[1]+cellsize/2,' -TARGET_USER_XMAX ', te[3]-cellsize/2, ' -TARGET_USER_YMIN ', te[2]+cellsize/2,' -TARGET_USER_YMAX ', te[4]-cellsize/2))
system(paste0('gdal_translate /data/Canada/tmp/DataForTom/DSS_soil_polygons_30m.sdat /data/Canada/tmp/DataForTom/covariates/DSS_soil_polygons_30m.tif -ot \\"Byte\\" -co \\"COMPRESS=DEFLATE\\"'))
unlink("DSS_soil_polygons_30m.*")

## Climatic layers:
system(paste0('gdalwarp /mnt/DATA/GlobalForestChange2000-2014/treecover2000.vrt /data/Canada/tmp/DataForTom/covariates/PEI_treecover2000_250m.tif -r \\"average\\" -t_srs \\"', proj4string(r30m), '\\" -co \\"BIGTIFF=YES\\" -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -co \\"COMPRESS=DEFLATE\\" -tr 250 250 -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))
r250m = raster("covariates/PEI_treecover2000_250m.tif")

clm.lst = c(list.files("/mnt/DATA/LandGIS/layers1km", pattern = glob2rx("clm_lst_mod11a2.*.night_m_1km_s0..0cm_2000..2017_v1.0.tif$"), full.names=TRUE), list.files("/mnt/DATA/LandGIS/layers1km", pattern = glob2rx("clm_precipitation_sm2rain.*_m_1km_s0..0cm_2007..2018_v0.2.tif$"), full.names=TRUE))
x = parallel::mclapply(clm.lst, function(j){  system(paste0('gdalwarp ', j,' /data/Canada/tmp/DataForTom/covariates/PEI_', basename(j), '_250m.tif -t_srs \\"', proj4string(r250m), '\\" -co \\"BIGTIFF=YES\\" -multi -wo \\"NUM_THREADS=4\\" -r \\"cubicspline\\" -wm 2000 -co \\"COMPRESS=DEFLATE\\" -tr 250 250 -te ', paste(as.vector(extent(r250m))[c(1,3,2,4)], collapse=" "))) }, mc.cores = length(clm.lst)/2 )
#plot(raster("covariates/PEI_clm_lst_mod11a2.aug.night_m_1km_s0..0cm_2000..2017_v1.0.tif_250m.tif"))

## stack everything together:
tif30 = list.files("/data/Canada/tmp/DataForTom/covariates", pattern = glob2rx("*.tif$"), full.names = TRUE)
zip()
tif30 = tif30[-grep("250m", tif30)]
r30m.sp = stack(tif30)
r30m.sp = as(r30m.sp, "SpatialGridDataFrame")
summary(as.factor(r30m.sp$DSS_soil_polygons_30m))
names(r30m.sp)
r30m.sp$SurficialGeology_code = as.factor(r30m.sp$SurficialGeology_code)
r30m.sp$DSS_soil_polygons_30m = plyr::join(data.frame(Value=r30m.sp$DSS_soil_polygons_30m), geo_leg)$Classes
r30m.sp$landcover_2016_reclassify = as.factor(r30m.sp$landcover_2016_reclassify)
names(r30m.sp)
plot(r30m.sp[5])
x30m.lst = lapply(r30m.sp@data, function(i){sd(as.numeric(i), na.rm=TRUE)})
r30m.sp = r30m.sp[!x30m.lst==0]

tif250 = list.files("/data/Canada/tmp/DataForTom/covariates", pattern = glob2rx("*250m.tif$"), full.names = TRUE)
r250m.sp = stack(tif250)
r250m.sp = as(r250m.sp, "SpatialGridDataFrame")
names(r250m.sp)
plot(r250m.sp[12])
for(i in 1:ncol(r250m.sp@data)){   r250m.sp@data[,i] <- ifelse(is.na(r250m.sp@data[,1]), NA, r250m.sp@data[,i]) }
plot(raster(r250m.sp[14]))
points(som, pch="+")

pei.soil = list(som.points=som, type.points=npdb, grids30m = r30m.sp, grids250m = r250m.sp)
saveRDS(pei.soil, "PEI.soil.rds")

## large file - try 100 m:
for(i in 1:length(tif30)){
  if(length(grep("polygons", tif30[i]))>0 | length(grep("Geology", tif30[i]))>0 | length(grep("landcover", tif30[i]))>0){
    system(paste0('gdalwarp ', tif30[i], ' ', gsub("covariates", "covariates100m", tif30[i]), ' -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -r \\"mode\\" -co \\"COMPRESS=DEFLATE\\" -tr 100 100 -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))
  } else {
    system(paste0('gdalwarp ', tif30[i], ' ', gsub("covariates", "covariates100m", tif30[i]), ' -multi -wo \\"NUM_THREADS=ALL_CPUS\\" -wm 2000 -r \\"average\\" -co \\"COMPRESS=DEFLATE\\" -tr 100 100 -te ', paste(as.vector(extent(r30m))[c(1,3,2,4)], collapse=" ")))
  }
}
tif100 = list.files("/data/Canada/tmp/DataForTom/covariates100m", pattern = glob2rx("*.tif$"), full.names = TRUE)
r100m.sp = stack(tif100)
r100m.sp = as(r100m.sp, "SpatialGridDataFrame")
summary(as.factor(r100m.sp$DSS_soil_polygons_30m))
r100m.sp$SurficialGeology_code = as.factor(r100m.sp$SurficialGeology_code)
r100m.sp$DSS_soil_polygons_30m = plyr::join(data.frame(Value=r100m.sp$DSS_soil_polygons_30m), geo_leg)$Classes
r100m.sp$landcover_2016_reclassify = as.factor(r100m.sp$landcover_2016_reclassify)
#plot(r100m.sp[7])
x100m.lst = lapply(r100m.sp@data, function(i){sd(as.numeric(i), na.rm=TRUE)})
r100m.sp = r100m.sp[!x100m.lst==0]
names(r100m.sp)
plot(r100m.sp[14])
pei100m.soil = list(som.points=som, type.points=npdb, grids100m = r100m.sp, grids250m = r250m.sp)
saveRDS(pei100m.soil, "PEI100m.soil.rds")
