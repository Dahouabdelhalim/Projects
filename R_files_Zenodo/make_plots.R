### feed multiple clipped rasters in to match Curtis' polygons
### change distance from 200 to 150m because some of these polygons are small
### mod by DSG Nov 6-11 2013, 6 Dec 2013, 11 Dec 2013, 30 Jan 2014
### new *.R created using Rominger's revisions 4 Feb 2014 to
### find_pixels.R  &  find_rng.R  &  gnerate_plots.R

rm(list=ls())

setwd('~/Dropbox/plot_selection')
source('find_rng.R')      # function for finding desired canopy height range
source('find_pixels.R')   # function for finding random pixels with desired canopy range
source('generate_plots.R') # function to make candidate plots


## create base directory for appending output directories
base <- getwd()

##  Use on Curtis' polygons with clip of CAO raster (mean height)
##  List of files in ./SitePolygons_Feb2014

##   RUN FOR EACH SITE
all.file <- list.files('cao', pattern='.tif')
all.file <- all.file[!grepl('\\\\.xml', all.file)]

for(i in all.file) { # first one ran
	print(i)
	
	base.name <- gsub('_CAOMN.tif', '', i)
	
	OUT.name <- paste(base.name, 'PTS200', sep='_')
	RAST.name <- paste('cao/', i, sep='')
	POLY.name <- paste('SitePolygons_Feb2014/', base.name, '.shp', sep='')
	
	# create filename for the output
	OUT <- paste(base, 'run2014-04-10', OUT.name, sep="/")
	
	# load CAO raster
	RAST <- raster(RAST.name)
	
	# load one of Curtis's site polygons
	POLY <- readShapeSpatial(POLY.name,
	    proj4string=CRS('+proj=utm +zone=5 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))
	
	# change projection as needed
	if(proj4string(POLY) != proj4string(RAST)) {
		POLY <- spTransform(POLY, CRS(proj4string(RAST)))
	}
	
	set.seed(103)
	x <- generate.plots(n=20, d=200, q=c(0.6,0.95), r=0, site.poly=POLY, cao=RAST,
    	                flow.poly=NULL, file=OUT)
}
