###############################################################################
###                                                                         ###   
###        R-script to translate SMAP HDF5 files into GeoTIFF format        ###
###                                                                         ###
###############################################################################

# 11 January, 2022

# Marian Schönauer
# Department  of  Forest  Work  Science  and  Engineering 
# University  of  Göttingen,  Göttingen
# E-mail address: marian.schoenauer@uni-goettingen.de

# This code is free for use and can be cited as:

# Schönauer, M., 2022. R-script to translate SMAP HDF5 files into GeoTIFF format. Zenodo. doi:10.5281/zenodo.5835820



#install.packages("BiocManager")
#BiocManager::install("rhdf5")
library("rhdf5")
library("raster")
library("rgdal")
library("gdalUtils")
library("smapr")
library("dplyr")

rm(list=ls())
path <- paste("D:/SMAP_L4/")

# list of h5 files
files <- list.files(path,
                    pattern = "\\\\.h5$",
                    recursive = TRUE)

# overview of the hdf5 file
rhdf5::h5ls(paste0(path, files[1]))

# select long-lat values
latitude <- values(raster(rhdf5::h5read(paste0(path, files[1]), name = "/cell_lat")))
longitude <- values(raster(rhdf5::h5read(paste0(path, files[1]), name = "/cell_lon")))
dat<-data.frame(lon = c(min(longitude, na.rm=T),max(longitude, na.rm=T)), # upper left extent
                lat = c(max(latitude, na.rm=T),min(latitude, na.rm=T))) # lower right extent
cord.dec = SpatialPoints(cbind(dat$lon, dat$lat), proj4string = CRS("+proj=longlat"))
cord<-spTransform(cord.dec, CRS("+init=epsg:6933"))

# function to translate every file contained in the list 'files'
doit<-function(x) {
  f <- paste0(path, x)
  v1<-rhdf5::h5read(f, name = "/Analysis_Data/sm_surface_analysis")
  out1 <- t(v1)
  outr <- raster(out1)
  crs(outr) <- "+proj=longlat"
  values(outr) <- ifelse(values(outr)<0, NA, values(outr))
  pathTiff<-file.path(paste0(path, "0tiff2/sm_surface_", substr(f, 45, 55)))
  writeRaster(outr, pathTiff, overwrite = T, format = "GTiff")
  gdal_translate(src_dataset = paste0(pathTiff, ".tif"),
                 dst_dataset = file.path(paste0(pathTiff, "_ref.tif")), # set filename here, (1) _am.tif (2) _pm.tif
                 of = "GTiff",
                 a_ullr = c(xmin(cord), ymax(cord), xmax(cord), ymin(cord)),
                 a_srs = "+proj=cea +lon_0=0 +lat_ts=30 +ellps=WGS84 +units=m",
                 co="TILED=YES",
                 verbose=TRUE)
  file.remove(paste0(pathTiff, ".tif"))
}

# apply the function to the list of input files
sapply(files, doit)
