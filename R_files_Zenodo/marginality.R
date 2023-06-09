autocrop <-
function(dmap,pops,margin=1) 
{
# Crops a distribution map by removing all outside areas where the species is
# not present
# dmap: distribution map of the species
# pops: geographical coordinates of the marginal populations (optional)
# margin: size (in pixels) of the outer margin left around the cropping window

x <- as.matrix(dmap)
xr <- range((1:dim(x)[2])[as.logical(colSums(x,na.rm=TRUE))])
yr <- range((1:dim(x)[1])[as.logical(rowSums(x,na.rm=TRUE))])
ext <- extent(dmap,yr[1],yr[2],xr[1],xr[2])
if(!missing(pops)) {
	if((p<-proj4string(dmap))!=proj4string(pops)) pops <- spTransform(pops,CRS(p))
	r <- res(dmap)
	ext <- extent(min(xmin(ext),xmin(pops)-margin*r[1]),
		max(xmax(ext),xmax(pops)+margin*r[1]),
		min(ymin(ext),ymin(pops)-margin*r[2]),
		max(ymax(ext),ymax(pops)+margin*r[2])) }
crop(dmap,ext)
}
dist.area <-
function(mspa,method=c("grid","direct"),...) 
{
# Area (in ha) of the nearest core
# mspa: raster with the MSPA of a distribution map that classifies pixels as:
#       0: background (no forest)
#       2: edge
#       3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#       10, 11, 12... : labels of core patches
#       NA: no data
# method: method for computing distances, either as paths going through the
#       centers of pixels, or directly as the Euclidean distance
# ...: additional parameters for the function 'dist.cores'
# OUTPUT: raster giving the map of the area of the nearest core

dist.cores(mspa,method=method,what="area.nearest.core",...)$area
}
dist.centroid <-
function(dmap,conductance=c(forest=100,land=50,sea=10),id=NULL,bbmask=NULL,
	method=c("grass","R")) 
{
# Least-cost distance to the centroid 
# dmap: raster with
#       NA: no land (sea)
#       0: non-forest land
#       1: forest
# conductance: conductance values for forest, non-forest land, and sea.
# bbmask: optional mask for the bounding box of map
# OUTPUT: raster giving the map of distances from the centroid.

# calculating the centroid of the distribution map
centroid <- SpatialPoints(matrix(colMeans(xyFromCell(dmap,Which(dmap==1,
	cells=TRUE,na.rm=TRUE))),ncol=2),proj4string=crs(dmap))

# preparing the conductance raster 
# dmap <- subs(dmap,data.frame(id=0:1,v=conductance[c("land","forest")]))
# This computation would be the preferred one because it prevents any
# confusion between the different conductance values, but it is more
# memory consuming and can crash for large maps. Therefore, the following
# computation is used:
if(missing(bbmask)) bbmask <- dmap
dmap[dmap==1] <- conductance["forest"]
dmap[dmap==0] <- conductance["land"] 
dmap[is.na(dmap)] <- conductance["sea"]

# doing computations externally with GRASS GIS
dist.centroid.grass <- function(dmap,id) {
	# directory where are GRASS files are located
	GRASS.rep <- "C:/Program Files/GRASS GIS 7.8" # change to your own settings
	if(is.null(id)) id <- basename(tempfile("species"))
	library(rgrass7)
	writeRaster(1/dmap,"tmp.tif",overwrite=TRUE) # cost = 1/conductance
	initGRASS(GRASS.rep,home=tempdir(),gisDbase=tempdir(),override=TRUE)
	execGRASS("r.in.gdal",flags="overwrite",parameters=list(input="tmp.tif",
		output=id,location=id))
	initGRASS(GRASS.rep,home=tempdir(),gisDbase=tempdir(),location=id,mapset=
		"PERMANENT",override=TRUE)
	execGRASS("r.cost",flags="k",parameters=list(input=id,output="dcentroid",
		start_coordinates=c(centroid@coords)))
	execGRASS("r.out.gdal",flags="overwrite",parameters=list(input="dcentroid",
		output="tmp-centroid.tif"))
	return(raster("tmp-centroid.tif")) }

# doing computation internally in R
dist.centroid.R <- function(dmap,bbmask) {
	# preparing the conductance raster 
	# dmap <- subs(dmap,data.frame(id=0:1,v=conductance[c("land","forest")]))
	# This computation would be the preferred one because it prevents any
	# confusion between the different conductance values, but it is more
	# memory consuming and can crash for large maps. Therefore, the following
	# computation is used:
	dmap[dmap==1] <- conductance["forest"]
	dmap[dmap==0] <- conductance["land"] 
	dmap[is.na(dmap)] <- conductance["sea"]
	library(gdistance)
	while(TRUE) {
		# preparing the transition matrix
		Tr <- try(transition(dmap,transitionFunction=mean,directions=8),silent=TRUE)
		if(class(Tr)!="try-error") Tr <- try(geoCorrection(Tr,type="c"),silent=TRUE) 
		if(class(Tr)!="try-error") break
		cat("Conductance matrix too large to be allocated in memory. Aggregating pixels 2 by 2.\\n"); flush.console()
		dmap <- raster::aggregate(dmap,fact=2) 
		if(!is.null(bbmask)) bbmask <- raster::aggregate(bbmask,fact=2) }
	return(accCost(Tr,centroid)) }

if("grass" %in% method) {
	x <- try(dist.centroid.grass(dmap,id),silent=TRUE)
	if(class(x)=="try-error") x <- dist.centroid.R(dmap,bbmask) }
else x <- dist.centroid.R(dmap,bbmask)
if(is.null(bbmask)) return(x) else return(mask(x,bbmask))
}
dist.cores <-
function(mspa,sizemin=100,dmin=50,power=2,method=c("proximity","grid","direct"),id,workdir,
	what=c("area.nearest.core","dist.nearest.core","gravity","isolation"),
	keep=10) 
{
# File prefixes:
# d: distance to the nearest core
# a: area of the nearest core
# l: distance to the nearest core greater than 'sizemin' and further than 'dmin'
# w: sum of weight 1/d^power
# m: sum of weight 1/d^power times area of the core 
# f: distance to nearest core greater than 'sizemin'
# s: distance to second nearest core greater than 'sizemin'
# keep: number of steps to be kept as files before being cleaned

method <- match.arg(method)
if((method=="proximity") & (projsys(mspa)!="laea")) stop("Cannot use proximity method with non-metric projection system")
ext <- c(if("area.nearest.core" %in% what) c("d","a"),
	if("dist.nearest.core" %in% what) "l",
	if("gravity" %in% what) c("w","m"),
	if("isolation" %in% what) c("f","s"))
if(missing(id)) id <- basename(tempfile("species"))
if(missing(workdir)) workdir <- getOption("rasterTmpDir")
coredir <- paste0(workdir,id)
if(!dir.exists(coredir)) dir.create(coredir)
namefile <- function(no) structure(paste0(coredir,"/",ext,"core",no,".img"),names=ext)
rasterOptions(progress="text")
rasterOptions(todisk=TRUE)
rasterOptions(memfrac=0.1)

# Checking where computations stopped (in case computations were not achieved in a single slot)
done <- list.files(coredir)
done <- done[grep("core\\\\d+\\\\.img",done)]
done <- if(length(done)==0) 0 else max(as.numeric(gsub("\\\\D","",done)))
if(done) {
	filedone <- namefile(done)
	cat("Computations already done till core",done,"\\n")
	if("area.nearest.core" %in% what) {
		dcur <- raster(filedone["d"])   # map of the distance to the nearest core
		acur <- raster(filedone["a"]) } # area of the nearest core
	if("dist.nearest.core" %in% what) {
		lcur <- raster(filedone["l"]) } # map of the distance to the nearest core greater than 'sizemin' and further than 'dmin'
	if("gravity" %in% what) {
		wcur <- raster(filedone["w"])   # sum of weight 1/d^power
		mcur <- raster(filedone["m"]) } # sum of weight 1/d^power times area of the core
	if("isolation" %in% what) {
		fcur <- raster(filedone["f"])   # distance to nearest large core
		scur <- raster(filedone["s"]) }}# distance to second nearest large core
else { # starting computations from scratch
	 # the value of 50,000 km (more than the distance between any two points 
	 # on Earth) is used rather than 'Inf' because 'Inf' values in rasters 
	 # are confused with NA values
	if("area.nearest.core" %in% what) {
		dcur <- acur <- mspa
		dcur[] <- 50e3
		acur[] <- NA }
	if("dist.nearest.core" %in% what) {
		lcur <- mspa 
		lcur[] <- 50e3 } 	
	if("gravity" %in% what) {
		wcur <- mcur <- mspa
		wcur[] <- mcur[] <- 0 }
	if("isolation" %in% what) {
		fcur <- scur <- mspa
		fcur[] <- scur[] <- 50e3  }}

# Area of the cores in hectares
mspa[mspa < 10] <- NA # removing everything but cores
a <- switch(projsys(mspa),
	longlat = {
		b <- as.data.frame(zonal(area(mspa,na.rm=TRUE),mspa,fun="sum",na.rm=TRUE)) # area in km2
		b$sum <- b$sum*100 # conversion km2 -> ha 
		b }, 
	laea = {
		b <- table(mspa[])*prod(res(mspa))/1e4 # conversion m2 -> ha
		data.frame(zone=as.numeric(names(b)),sum=c(b)) })
a <- a[order(a$zone),]	

# Main loop on cores
cat(nrow(a),"cores to process\\n")
for(i in 1:nrow(a)) {
	core <- a$zone[i]
	if(core <= done) next
	cat("Computing distance to core",core,"[",i,"/",nrow(a),"]\\n"); flush.console()
	dcore <- mspa
	dcore[dcore != core] <- NA # removing everything but the current core
	d <- switch(method, # map of distances to the current core
		grid = gridDistance(boundaries(dcore),origin=1,omit=0),
		direct = distance(boundaries(dcore)), 
		proximity = proximity(dcore,core))
	d <- d/1000 # conversion m -> km

	# updating the rasters
	if("gravity" %in% what) {
		p <- 1/d^power
		wcur <- wcur+p # updating the sum of weigths
		mcur <- mcur+p*a$sum[i] } # updating the sum of wieghts times the areas of cores
	d[is.na(d)] <- 0 # setting distance to 0 inside the core
	if("isolation" %in% what) if(a$sum[i] >= sizemin) {
		# if the core is too small, no further action is needed
		# updating the maps of distance to the first and second nearest large core
		# when the current core is the nearest
		fi <- (d < fcur)
		o <- try(scur[fi] <- fcur[fi],silent=TRUE)
		if(class(o)=="try-error") {
			cat(o)
			scur <- mask(scur,fi,maskvalue=1)
			scur <- cover(scur,fcur) }
		o <- try(fcur[fi] <- d[fi],silent=TRUE)
		if(class(o)=="try-error") {
			cat(o)
			fcur <- mask(fcur,fi,maskvalue=1) 
			fcur <- cover(fcur,d) }
		# updating the maps of distance to the second nearest large core when 
		# the current core is the second nearest
		fi <- !fi & (d < scur)
		o <- try(scur[fi] <- d[fi],silent=TRUE)
		if(class(o)=="try-error") {
			cat(o)
			scur <- mask(scur,fi,maskvalue=1) 
			scur <- cover(scur,d) }}
	if("area.nearest.core" %in% what) {
		acur[d < dcur] <- a$sum[i] # updating the area of the nearest core
		dcur[] <- pmin(dcur[],d[]) } # updating the map of distances to the nearest core
	if("dist.nearest.core" %in% what) if(a$sum[i] >= sizemin) {
		# if the core is too small, no further action is needed
		d[d < dmin] <- 50e3 # excluding points closer than 'dmin'
			# the value of 50,000 km (more than the distance between any
			# two points on Earth) is used rather than 'Inf' because
			# 'Inf' values in rasters are confused with NA 
		lcur[] <- pmin(lcur[],d[]) }

	# saving the updated rasters and cleaning the previous ones
	fich <- namefile(core)
	if("area.nearest.core" %in% what) {
		writeRaster(dcur,file=fich["d"],overwrite=TRUE)
		writeRaster(acur,file=fich["a"],overwrite=TRUE) }
	if("dist.nearest.core" %in% what) {
		writeRaster(lcur,file=fich["l"],overwrite=TRUE) }
	if("gravity" %in% what) {
		writeRaster(wcur,file=fich["w"],overwrite=TRUE) 
		writeRaster(mcur,file=fich["m"],overwrite=TRUE) }
	if("isolation" %in% what) {
		writeRaster(fcur,file=fich["f"],overwrite=TRUE)
		writeRaster(scur,file=fich["s"],overwrite=TRUE) }
	if(any(file.exists(cls <- namefile(a$zone[max(0,i-keep)])))) unlink(cls) }

filedone <- namefile(a$zone[nrow(a)])
list(area=if("area.nearest.core" %in% what) raster(filedone["a"]) else NULL,
	nearest=if("dist.nearest.core" %in% what) raster(filedone["l"]) else NULL,
	gravity=if("gravity" %in% what) {
		mspa <- subs(mspa,a,by="zone",which="sum") # raster map of the area of the cores
		p <- raster(filedone["m"])/raster(filedone["w"])
		isna <- is.na(p)
		p[isna] <- mspa[isna] 
		p } else NULL,
	isolation1=if("isolation" %in% what) raster(filedone["f"]) else NULL,
	isolation2=if("isolation" %in% what) raster(filedone["s"]) else NULL)
}
dist.edge <-
function(mspa,method=c("proximity","grid","direct")) 
{
# Distance to the border
# mspa: raster with the MSPA of a distribution map that classifies pixels as:
#	0: background (no forest)
#	2: edge
#	3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#	10, 11, 12... : labels of core patches
#	NA: no data
# method: method for computing distances, either as paths going through the
#	centers of pixels, or directly as the Euclidean distance
# OUTPUT: raster giving the map of distances from the border as a positive 
#	number if the pixel is inside a core and as a negative number if the
#	pixel is outside a core.

method <- match.arg(method)
d <- mspa
d[d!=2] <- NA
d <- switch(method,
	grid = gridDistance(d,origin=2),
	direct = distance(d),
	proximity = proximity(d,origin=2))
mspa[(mspa<10) | is.na(mspa)] <- 0
overlay(d,mspa,fun=function(x,y) ifelse(y>0,x,-x)) # distance inside
	# cores will be positive, distance outside cores will be negative
}
dist.gravity <-
function(mspa,power=2,method=c("grid","direct"),...) 
{
# Gravity index = weighted mean of the areas of the cores where the weights at
# each location are defined as the inverse of a power of the distances from
# that location to the cores
# mspa: raster with the MSPA of a distribution map that classifies pixels as:
#       0: background (no forest)
#       2: edge
#       3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#       10, 11, 12... : labels of core patches
#       NA: no data
# power: power to raise the distances to get the inverse of the weights
# method: method for computing distances, either as paths going through the
#       centers of pixels, or directly as the Euclidean distance
# ...: additional parameters for the function 'dist.cores'
# OUTPUT: raster giving the gravity index

dist.cores(mspa,power=power,method=method,what="gravity",...)$gravity
}
dist.history <-
function(dmap,what=c("lat","lon"),rproj) 
{
# Computes the North/South edge index and the West/East edge index
# dmap: distribution map, i.e. raster with cell values as follows:
#       NA: no land (sea)
#       0: non-forest land
#       1: forest

if(projsys(dmap)!="longlat") stop("use a longitude-latitude projection system")
what <- match.arg(what)
library(raster)
rasterOptions(progress="text")
rasterOptions(todisk=TRUE)
rasterOptions(memfrac=0.1)
p <- xyFromCell(dmap,Which(dmap==1,cells=TRUE,na.rm=TRUE))[,c(lon=1,lat=2)[what]]
x <- init(dmap,switch(what,lon='x',lat='y'))
library(EnvStats)
x <- calc(x,function(x) pemp(x,p))
if(!missing(rproj)) x <- projectRaster(x,rproj)
x
}
dist.isolation <-
function(mspa,sizemin=100,method=c("proximity","grid","direct","extern"),...) 
{
# Distance to the second nearest core greater than 'sizemin'
# mspa: raster with the MSPA of a distribution map that classifies pixels as:
#       0: background (no forest)
#       2: edge
#       3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#       10, 11, 12... : labels of core patches
#       NA: no data
# method: method for computing distances, either as paths going through the
#       centers of pixels ('grid'), or directly as the Euclidean distance
#	  ('direct'), or using the GDAL proximity function ('proximity'), or 
#	  using the GDAL proximity externally to R ('extern')
# ...: additional parameters for the function 'dist.cores'
# OUTPUT: raster map giving the distance to the second nearest core

if(match.arg(method)!="extern") return(dist.cores(mspa,sizemin=sizemin,method=
	method,what="isolation",...)$isolation2)

# Interface to rgdal_proximity to do computations externally to R
block <- if(hasArg(block)) list(...)$block else 10
workdir <- if(hasArg(workdir)) list(...)$workdir else paste0(getOption(
	"rasterTmpDir"),basename(tempfile("species")))
workdir <- normalizePath(paste0(getwd(),"/",workdir),mustWork=FALSE)
if(!dir.exists(workdir)) dir.create(workdir)
rasterOptions(progress="text")
# area of the cores in hectares
mspa[mspa < 10] <- NA # removing everything but cores
if(hasArg(a)) a <- list(...)$a else {
	a <- switch(projsys(mspa),
	        longlat = {
	                b <- as.data.frame(zonal(area(mspa,na.rm=TRUE),mspa,fun="sum",na.rm=TRUE)) # area in km2
	                b$sum <- b$sum*100 # conversion km2 -> ha 
	                b }, 
	        laea = {
	                b <- table(mspa[])*prod(res(mspa))/1e4 # conversion m2 -> ha
	                data.frame(zone=as.numeric(names(b)),sum=c(b)) })
	a <- a[a$sum >= sizemin,]
	a <- a[order(a$zone),]   }
N <- nrow(a)
# initializing variables
cores.done <- if(file.exists(fich <- paste0(workdir,"/done.txt"))) scan(fich,
	quiet=TRUE) else NULL
if(file.exists(fich <- paste0(workdir,"/fcur.img"))) {
	if(is.null(cores.done)) stop("Inconsistency in the cores already processed")
	fcur <- raster(fich) }
else {
	if(!is.null(cores.done)) stop("Inconsistency in the cores already processed")
	# the value of 50,000,000 m (more than the distance between any two points 
	# on Earth) is used rather than 'Inf' because 'Inf' values in rasters 
	# are confused with NA values
	fcur <- mspa; fcur[] <- 50e6 }
if(file.exists(fich <- paste0(workdir,"/scur.img"))) {
	if(is.null(cores.done)) stop("Inconsistency in the cores already processed")
	scur <- raster(fich) }
else {
	if(!is.null(cores.done)) stop("Inconsistency in the cores already processed")
	scur <- mspa; scur[] <- 50e6 }
cat(length(cores.done),"cores already processed out of",nrow(a),"\\n")
a <- a[!(a$zone %in% cores.done),]
while(nrow(a)) {
	k <- min(block,nrow(a))
	# processing one block of k cores
	proximity(mspa,origin=a$zone[1:k],REP=workdir,outfile=paste0("core",
		a$zone[1:k],".tif"))
	for(i in 1:k) {
		core <- a$zone[i]
		n <- N-nrow(a)+i
		cat("Processing core",core,"[",i,"/",k,"][",n,"/",N,"]\\n"); flush.console()
		d <- raster(paste0(workdir,"/core",core,".tif"))
		# updating the maps of distance to the first and second nearest large core
		# when the current core is the nearest
		fi <- (d < fcur)
		if(n >= 2) {
			o <- try(scur[fi] <- fcur[fi],silent=TRUE)
			if(class(o)=="try-error") {
				cat(o)
			scur <- mask(scur,fi,maskvalue=1)
			scur <- cover(scur,fcur) }}
		o <- try(fcur[fi] <- d[fi],silent=TRUE)
		if(class(o)=="try-error") {
                cat(o)
			fcur <- mask(fcur,fi,maskvalue=1)
			fcur <- cover(fcur,d) }
		# updating the maps of distance to the second nearest large core when 
		# the current core is the second nearest
		if(n >= 2) {
			fi <- !fi & (d < scur)
			o <- try(scur[fi] <- d[fi],silent=TRUE)
			if(class(o)=="try-error") {
				cat(o)
				scur <- mask(scur,fi,maskvalue=1)
				scur <- cover(scur,d) }}}
	# updating the files
	cat("/!\\\\ Writing files: Do not interrupt now!\\n"); flush.console()
	writeRaster(fcur,paste0(workdir,"/fcur.img"),overwrite=TRUE)
	writeRaster(scur,paste0(workdir,"/scur.img"),overwrite=TRUE)
	write(a$zone[1:k],paste0(workdir,"/done.txt"),ncol=1,append=TRUE)
	cat("Files written, you can interrupt now...\\n"); flush.console() 
	file.remove(paste0(workdir,"/core",a$zone[1:k],".tif"))
	a <- a[-(1:k),,drop=FALSE] }
# conversion m -> km
fcur <- fcur/1000
scur <- scur/1000
list(isolation1=fcur,isolation2=scur)
}
dist.nearest <-
function(mspa,sizemin=100,dmin=50,method=c("grid","direct"),...) 
{
# Distance to the nearest large isolated core
# mspa: raster with the MSPA of a distribution map that classifies pixels as:
#       0: background (no forest)
#       2: edge
#       3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#       10, 11, 12... : labels of core patches
#       NA: no data
# sizemin: minimum area (in ha) of the nearest core
# dmin: minimum distance (in km) to the nearest core
# method: method for computing distances, either as paths going through the
#       centers of pixels, or directly as the Euclidean distance
# ...: additional parameters for the function 'dist.cores'
# OUTPUT: raster giving the map of the distance to the nearest core greater 
#	    than 'sizemin' and further than 'dmin'

dist.cores(mspa,sizemin=sizemin,dmin=dmin,method=method,what="dist.nearest.core",...)$nearest
}
is.valid.dmap <-
function(x) 
{
# Controls whether the distribution map of a species is in a valid format
# x: raster object with the distribution map

v <- unique(x[])
v <- v[!is.na(v)]
ok <- v %in% c(0,1)
if(!all(ok)) cat("Non-valid values:",v[!ok],"\\n")
return(all(ok))
}
mspa <-
function(x,neigh=8,npix=1,label.patches=TRUE,method=c("R","Guidos"),plot=TRUE,
	imp=FALSE,consider.water.as.background=TRUE,inner.edge=TRUE) 
{
# x: matrix or object of class 'raster' with the following pixel values:
#	0: land cell where the species is absent
#	1: land cell where the species is present
#	NA: non-land cell (i.e. water)
# neigh: type of neighbourhood. Either 4 (Von Neumann neighbourhood) or
#	8 (Moore neighbourhood)
# npix: width of the edge in pixels
# label.patches: if TRUE, all cores are labelled (sequentially starting from 10);
#	otherwise, the value of 1 is returned wherever the cell is part of a core
# consider.water.as.background: if TRUE, non-land cells (i.e. water) will be
#	considered as cells where the species is absent (i.e. background) rather
#	than as missing data. It makes a difference for the identification of
#	edges.
# inner.edge: if TRUE, inner edges (i.e. perforation) will be considered as
#	edges. The R implementation of MSPA only considers outer edges as
#	edges, so the Guidos implementation of MSPA is always used when 
#	'inner.edge' is TRUE.
# OUTPUT: matrix or object of class 'raster' that classifies the pixels as:
#	0: background (no forest)
#	1 or 10, 11, 12... : core or label of core patches
#	2: edge
#	3: neither core nor edge, i.e. islet, loop, bridge, perforation or branch
#	NA: no data
# Cores and edges are defined as in the MSPA of GuidosToolbox

library(raster)
if(!require(EBImage)) { 
	# See https://doi.org/10.18129/B9.bioc.EBImage 
	install.packages("BiocManager")
	BiocManager::install("EBImage") 
	library(EBImage) }
if(!(neigh %in% c(4,8))) stop("Neighbourhood must be 4 or 8")
method <- if(inner.edge) "Guidos" else match.arg(method)
what <- class(x)
crs <- try(proj4string(x),silent=TRUE)
ext <- try(extent(x),silent=TRUE)
if(consider.water.as.background) { water<-x; x[is.na(x)]<-0 }
switch(method,
	Guidos = {
		# Re-encoding map using Guidos input coding, i.e.
		# 0: missing data
		# 1: background (cell where the species is absent)
		# 2: foreground (cell where the species is present)
		x <- reclassify(x,rbind(c(0,1),c(1,2)))
		if(!consider.water.as.background) x[is.na(x)] <- 0 
		# GuidosToolbox MSPA can be downloaded at https://forest.jrc.ec.europa.eu/en/activities/lpa/mspa
		GuidosREP <- c("C:/GuidosToolbox/MSPAstandalone",paste0(.libPaths(),
			"/RMSPA/resources")) # adjust to your own settings
		GuidosREP <- GuidosREP[dir.exists(GuidosREP)][1]
		if(length(GuidosREP)==0) stop("could not find where GuidosToolbox MSPA is installed")
		com <- paste("mspa_win64.exe -eew",npix,"-graphfg",neigh,"-transition 0 -internal 0 -i input.tif -o output.tif")
		cat("Writing raster...\\n"); flush.console()
		writeRaster(x,paste0(GuidosREP,"/input.tif"),datatype="INT1U",overwrite=TRUE)
		wd <- getwd()
		setwd(GuidosREP)
		file.remove("output.tif")
		system(com,intern=TRUE)
		setwd(wd)
		cat("Reading raster...\\n")
		x <- raster(paste0(GuidosREP,"/output.tif"))
		# Re-encoding Guidos output where
		# 0: land cell where the species is absent (background)
		# 1: branch
		# 3: edge
		# 5: perforation
		# 9: islet
		# 17: core
		# 33: bridge
		# 35: bridge in edge
		# 37: bridge in perforation
		# 65: loop
		# 67: loop in edge
		# 69: loop in perforation
		# 129: non-land cell (missing data)
		x[x==129] <- NA
		x <- reclassify(x,rbind(c(1,3),c(3,2),c(5,if(inner.edge) 2 else 3),
			c(9,3),c(17,1),c(33,3),c(35,2),c(37,3),c(65,3),c(67,2),c(69,3)))
		if(label.patches) {
			core <- x <- as.matrix(x)
			core[core >= 2] <- 0
			core <- num.cores(core)
			core[x==2] <- 2
			core[x==3] <- 3 }},
	R = {
		if(method=="R") neigh <- switch(as.character(neigh),
			"4" = matrix(c(0,1,0,1,1,1,0,1,0),3),
			"8" = matrix(1,3,3))
		core <- x <- as.matrix(x)

		# erosion to identify core zone
		for(k in 1:npix) core <- erode(core,neigh)
		core[core > 1] <- NA
		dilated <- core
		for(k in 1:npix) dilated <- dilate(dilated,neigh)
		dilated[dilated > 1] <- NA
	
		# flood-fill to identify edges
		edge <- core
		edge[is.na(edge)] <- 0
		n <- nrow(edge)
		p <- ncol(edge)
		i <- rbind(cbind(1:n,1),cbind(1:n,p))
		if(p > 2) i <- rbind(i,cbind(1,2:(p-1)),cbind(n,2:(p-1)))
		while(any(o <- (edge[i]==0))) edge <- floodFill(edge,i[which(o)[1],],2)
		edge[edge==1] <- 0
		edge[edge==2] <- 1
		edge <- edge*(x-core)*dilated
		
		# other (neither core, nor edge), i.e. islet, loop, bridge, perforation or branch
		other <- x-core-edge
		if(label.patches) core <- num.cores(core)
		core[as.logical(edge)] <- 2
		core[as.logical(other)] <- 3 })
		
if(plot | imp) { 
	if(imp) png("mspa.png",width=3000,height=3000,res=600)
	col <- c("grey80","forestgreen","red","blue")
	if(label.patches) plot(raster(core))
	else {
		plot(raster(core),legend=FALSE,col=col)
		par(xpd=NA)
		legend("right",c("background","core","edge","other","no data"),fill=
			c(col,"white"),inset=-0.25,bty="n")
		par(xpd=FALSE) }
	if(imp) dev.off() }
if(what=="matrix") return(core)
core <- raster(core,ext@xmin,ext@xmax,ext@ymin,ext@ymax,crs=CRS(crs))
if(consider.water.as.background) core <- mask(core,water)
return(core)
}
num.cores <-
function(x,start=10,method=c("auto","clump","flood"),verbose=TRUE) 
{
# Gives numbers (starting from 'start') to the connected forest components
# x: raster or matrix with the following pixel values:
#	0: land (no forest)
#	1: forest
#	NA: no land (sea)
# OUTPUT: matrix with the following pixel values:
#	0: land (no forest)
#	10, 11, 12...: labels of the cores
#	NA: no mand (sea)

method <- match.arg(method)
if(method=="auto") {
	y <- try(Recall(x,start,method="clump"),silent=TRUE)
	if(class(y)!="try-error") return(y)
	return(Recall(x,start,method="flood")) }
x <- as.matrix(x)
x <- switch(method,
	clump = {
		fi <- x==0
		fi[is.na(fi)] <- FALSE
		x <- as.matrix(clump(raster(x),directions=4,gaps=FALSE))+start-1
		x[fi] <- 0 
		x },
	flood = {
		fi <- is.na(x)
		x[fi] <- 0
		while(nrow(i <- which(x==1,arr.ind=TRUE))) {
			if(verbose) { cat(start,":",length(i),"\\n"); flush.console() }
			x <- floodFill(x,i[1,],start)
			start <- start+1 }
		x[fi] <- NA
		x })
}
projsys <-
function(x) 
{
sub("\\\\+proj=(\\\\w+?) .+","\\\\1",sp::proj4string(x))
}
proximity <-
function(x,origin=2,omit,...) 
{
# Produces a raster proximity map using 'gdal_proximity' (distributed with 
# QGIS). An alternative would be to use the 'proximity' function of the
# RQGIS package (see https://rdrr.io/github/juoe/spatialtools/man/proximity.html)
# but this package is no longer maintained and distributed with R.
# x: raster layer
# origin: value(s) of the cells from which the distance is calculated
# omit: value(s) of the cells for which the distance does not have to be computed
# ...: potential additional arguments such as
#	REP: directory where to run QGIS computations
#	infile: name(s) of the input raster file(s) to process
#	outfile: name(s) of the processed output raster file(s)
#	do.run: if FALSE, write the command file but does not call it
# The path to QGIS will have to be adjusted

qgispath <- c("C:\\\\Program Files\\\\QGIS 3.16","C:\\\\OSGeo4W64") # adjust the path to QGIS depending on your own settings
qgispath <- qgispath[dir.exists(qgispath)]
REP <- if(hasArg(REP)) list(...)$REP else tempdir()
infile <- if(hasArg(infile)) list(...)$infile else paste0(deparse(substitute(x)),".tif")
outfile <- if(hasArg(outfile)) list(...)$outfile else "output.tif"
do.run <- if(hasArg(do.run)) list(...)$do.run else TRUE
if(length(outfile)!=length(origin)) origin <- paste(origin,collapse=",")
com <- paste0(REP,"/proximity.bat")
write(paste0("CHDIR \\"",qgispath,"\\\\bin\\""),com)
# gdal_proximity does not overwrite output file, so be careful to remove
# any pre-existing output file
write(paste("DEL",paste0(REP,"\\\\",outfile)),com,append=TRUE)
# write("SET PYTHONPATH=",com,append=TRUE)
write(paste0("SET PYTHONPATH=",qgispath,"\\\\apps\\\\Python37\\\\Scripts"),com,append=TRUE)
write(paste0("SET PYTHONHOME=",qgispath,"\\\\apps\\\\Python37"),com,append=TRUE)
write(paste0("PATH ",qgispath,"\\\\apps\\\\Python37;",qgispath,"\\\\apps\\\\Python37\\\\Scripts;%PATH%"),com,append=TRUE)
write(paste("IF NOT EXIST",paste0(REP,"\\\\",outfile),
	"( python3 -m gdal_proximity -srcband 1 -distunits GEO -values",origin,
	"-nodata 0.0 -ot Float32 -of GTiff",paste0(REP,"\\\\",infile),paste0(REP,
	"\\\\",outfile),")"),com,append=TRUE)
cat("Writing input file...\\n"); flush.console()
writeRaster(x,paste0(REP,"/",infile),overwrite=TRUE)
if(!do.run) return()
system(paste0(REP,"/proximity.bat"),intern=TRUE)
if(length(outfile)==1) {
	cat("Reading output file...\\n"); flush.console() 
	d <- raster(paste0(REP,"/",outfile)) 
	if(!missing(omit)) d <- mask(d,x,maskvalue=omit,updatevalue=0)
	return(d) }
}
read.pops <-
function(file=file.choose(),plot=TRUE,imp=FALSE) 
{
# Reads the geographical coordinates of marginal populations
# file: ASCII file with comma as decimal sign and semicolon as column separator.
#	The file must include column labels on its first line with a column
#	labelled "lon" that gives the longitudes and a column labelled "lat"
#	that gives the latitudes.

x <- read.table(file,dec=",",sep=";",quote="\\"",header=TRUE)
if(any(fi <- duplicated(x))) {
	warning("Removing ",sum(fi)," duplicated populations at row ",paste(which(fi),collapse=", "))  
	x <- x[-fi,] }
v <- tolower(names(x))
i <- grep("lon",v)
if((length(i)==0) | (length(i)>1)) stop("unable to identify the column corresponding to longitudes")
j <- grep("lat",v)
if((length(j)==0) | (length(j)>1)) stop("unable to identify the column corresponding to latitudes")
library(sp)
z <- x[,-c(i,j),drop=FALSE]
# By default, the population name is supposed to be given in the first column
# that is neither longitude nor latitude
names(z)[1] <- "pop.name"
z$pop.name <- singularize(z$pop.name)
x <- SpatialPointsDataFrame(x[,c(i,j)],z,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
if(plot | imp) {
	library(maps)
	library(mapdata)
	if(imp) png("pops.png",width=3000,height=3000,res=600)
	plot(x,col="red") 
	box()
	axis(1)
	axis(2)
	map("world",add=TRUE)
	if(imp) dev.off() }
x
}
smooth <-
function(x,npix=2,plot=TRUE) 
{
# Smooths a distribution map by successively applying dilations and erosions

kern <- matrix(1,3,3)
crs <- proj4string(x)
ext <- extent(x)
if(plot) { par(mfrow=c(1,2)); plot(x) }
x <- as.matrix(x)
library(EBImage)
for(k in 1:npix) {
	x <- dilate(x,kern)
	x[abs(x) > 1] <- NA }
for(k in 1:npix) {
	x <- erode(x,kern)
	x[abs(x) > 1] <- NA }
x <- raster(x,ext@xmin,ext@xmax,ext@ymin,ext@ymax,crs=CRS(crs))
if(plot) plot(x)
x
}
