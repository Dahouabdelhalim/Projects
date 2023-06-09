##  function to generate randomized plots within pre-selected site
##  polygons. This function creates `n' plots that are minimum `d'
##  distance appart. Each plot is placed such that it falls within a
##  specified rang of canopy height values. That range is specified as
##  a quantile `q' [e.g. c(0.9,1) for the upper 10% quantile] of
##  canopy heights recorded in the CAO LiDAR data `cao.' The LiDAR
##  data are clipped to the shape of the lava flow we're interested
##  in, the lava flow polygon being given by `flow.poly.' Writes the
##  plot points to a file if argument `file' is provided

generate.plots <- function(n, d, q, r, site.poly, cao, flow.poly, file) {
    if(missing(r)) r <- 0
    
    if(is.null(flow.poly)) {
    		flow2use <- NULL
    } else {
    		flow2use <- findFlow(site.poly,flow.poly)
    }
    # browser()
    rng <- find.rng(x=cao, site=site.poly, q=q, r=r,
                    flow=flow2use, return.r=TRUE)
    plts <- rand.pxl(x=rng$r, n=n, d=d, rng=rng$rng)
    
    cat('\\n', attr(plts,'note'), ';\\n', sprintf('found %s plots',
                                                length(plts)), '\\n', sep='')
    
    vals <- values(rng$r)
    vals <- vals[!is.na(vals)]
    
    if(!missing('file')) {
        plts <- SpatialPointsDataFrame(plts,
                                       data=data.frame(type=rep('plot_centroid',
                                                           length(plts))),
                                       proj4string=CRS(proj4string(plts)))
        
        ## write spatial points to ESRI shapefile
        writeOGR(plts, file, gsub('.*/','',file), driver='ESRI Shapefile', overwrite_layer=TRUE)
        
        ## write raster values to text file
        write.table(vals, file=paste(file, '/', gsub('.*/','',file), '.txt', sep=''),
                    sep='\\t', row.names=FALSE, col.names=FALSE)
                    
        ## write range of values considered
        write.table(rng$rng, file=paste(file, '/', gsub('.*/','',file), '_rng.txt', sep=''),
                    sep='\\t', row.names=FALSE, col.names=FALSE)
    }
    
    out <- list(rng=rng$rng, vals=vals, plts=plts)
    return(out)
}


##  test it on simulated data
# setwd('~/Dropbox/hawaiidimensions/SiteSelection')
# source('find_rng.R')
# source('find_pixels.R')

# # load simulated CAO data
# laup.r <- raster('cao/laup_sim.tif')

# # load one of Curtis's site polygons
# laup.p <- readShapeSpatial('PreliminarySitePolygons/LaupahoehoeHipNet_4-14K.shp',
                           # proj4string=CRS('+proj=utm +zone=5 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

# # load flow polygons clipped by elevation
# elevGeo <- readShapeSpatial('elevFlow/geo3to5k.shp',
                            # proj4string=CRS('+proj=utm +zone=5 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'))

# ##  specifying the file will make the function write
# ##  the plots to an ESRI shapefile, you can also omit 
# ##  the file argument and then the function will return
# ##  the points in memory
# x <- generate.plots(n=20, d=200, q=c(0.6,0.95), site.poly=laup.p, cao=laup.r,
                    # flow.poly=elevGeo, file='./laupTest')
# x <- generate.plots(n=20, d=200, q=c(0.6,0.95), site.poly=laup.p, cao=laup.r,
                    # flow.poly=NULL, file='./laupTest')

# ##  I've also made the function return the range of
# ##  canopy heights being used (i.e. across the entire flow)
# ##  and the range of heights found within the flow
# x
