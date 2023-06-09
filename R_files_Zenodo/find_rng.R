library(sp)
library(rgdal)
library(raster)
library(maptools)

##  function to determine acceptable rng of canopy
##  heights based on criteria of:
##    x: the CAO data
##    site: the site polygon
##    q: quantile range we're aiming for (e.g. [0.9,1])
##    r: evaluate that quantile range within `r' distance from site
##    flow: consider all pixels within the flow that contains the site
##          (if flow is NULL, other constraints on pixles are used (e.g. site polygon))
##    all.val: determines if all values are returned or just the quantile specified by `q'
##    return.r: logical, should masked raster be returned
##
##  i'm keeping this function sepparate from the function that selects
##  candidate pixels so that we can evaluate how `r' changes our
##  acceptable range

find.rng <- function(x,site,q,r,flow=NULL,all.vals=FALSE,return.r=FALSE) {
	## determine which polygon to use as a mask for the CAO data
	# browser()
	if(!is.null(flow)) { # use the lava flow as mask
		poly2mask <- flow
		cat('using flow for mask \\n')
	} else if(max(apply(bbox(site),1,diff)) < 2*r) { # use a radius as mask
		cPoly <- coordinates(site)[1,]
		rPoly <- cbind(cos(seq(0,2,0.02)*pi)*r+cPoly[1], sin(seq(0,2,0.02)*pi)*r+cPoly[2])
		poly2mask <- SpatialPolygons(
			list(Polygons(list(rPoly=Polygon(rPoly)),ID=0)),
			proj4string=CRS(proj4string(site))
		)
		cat('using radius for mask \\n')
	} else { # use the site polygon as mask
		poly2mask <- site
		cat('using site for mask \\n')
	}
	# browser()
	## masking is a time-intensive step! cropping speeds it up some
	x <- crop(x,extent(bbox(poly2mask)))
	## do the actual masking
	x <- mask(x,poly2mask)
	
	## if raster already masked to site polygon no need for
	## further masking below
	if(max(apply(bbox(site),1,diff)) > 2*r & is.null(flow)) {
		r.out <- x
	}
	
	# the `!is.na' part is very important because all the raster
	# subsetting above is actually just setting out-of-bounds pixles
	# to NA
	vals <- values(x)
	vals <- vals[!is.na(vals)]
	
	if(all.vals) {
		out <- vals
	} else {
		## `rm.outliers' removes potential erronous pixels at
		## extremes of distribuitons before calculating quantile
		rng <- quantile(rm.outliers(vals), prob=q)
		names(rng) <- NULL
		
		out <- rng
	}
	
	## should the maksed raster be returned?
	if(return.r) {
		## further clip raster to site polygon if haven't already
		if(!exists('r.out')) {
			x <- crop(x,extent(bbox(site)))
			cat('further masking to site for output \\n')
			r.out <- mask(x,site)
		}
		out <- list(rng=rng,r=r.out)
	}
	
	return(out)
}


##  helper function to remove anamolys LiDAR data. Looks for point mass at 0
##  and also trailing right hand tail. Right had tail criteria is to cut
##  canopy height values that are > 0.9 quantile AND whoes frequencies deviate
##  from an exponential decay from the mode.

rm.outliers <- function(x) {
	# browser()
	x <- x[!is.na(x)]
	
	## remove potential point mass near 0
	tooSmall <- 0.045 * diff(range(x))
	if(sum(x < tooSmall) > 1.25*(sum(x < 2*tooSmall) - sum(x < tooSmall))) {
		x <- x[x > tooSmall]
	}
	
	## calculate frequencies excluding top 10% quantile
	brks <- seq(0, max(x)*1.1, by=1.5*tooSmall)
	brks <- c(brks, max(brks) + diff(brks[1:2]))
	x4tailCalc <- x[x < brks[which.min(abs(brks - quantile(x[x > mean(x)], 0.95)))]]
	xhist <- data.frame(hist(x4tailCalc, brks, 
                             plot=FALSE)[c('mids', 'counts')])
	
	## define `tail' as region of distribuiton right of mode
	xtail <- xhist[which.max(xhist[,2]) <= 1:nrow(xhist) & xhist[,2] > 0, ]
	xtail[,2] <- log(xtail[,2])
	
	## fit exponential to tail
	tail.mod <- lm(counts~mids, data=xtail)
	
	## calculate frequency for all data
	xhist.all <- data.frame(hist(x, brks, 
                                 plot=FALSE)[c('mids', 'counts')])
	xhist.all <- xhist.all[xhist.all$counts > 0, ]
	xhist.all[,2] <- log(xhist.all[,2])
	
	## calculate which canopy height classes are outliers
	tail.line <- predict(tail.mod, newdata=xhist.all[, 1, drop=FALSE])
	if(length(tail.mod$residuals) > 2) {
		tail.line <- tail.line + 2*sqrt(sum(tail.mod$residuals^2)/(length(tail.mod$residuals)-2))
	}
	
	bad.mid <- xhist.all[xhist.all[,2] > tail.line & xhist.all[,1] > quantile(x, 0.9), 1]
	
	## remove outliers from `x'
	if(length(bad.mid) > 0) {
		x <- x[x < min(bad.mid)]
	}
	
	return(x)
}

##  helper function to identify which polygon in large shapefile of flows
##  corresponds to the actual site polygon
findFlow <- function(site,flows) {
	ID <- over(SpatialPoints(coordinates(site),proj4string=CRS(proj4string(site))),flows)
	
	return(flows[flows$FID_geol_f == ID$FID_geol_f, ])
}
