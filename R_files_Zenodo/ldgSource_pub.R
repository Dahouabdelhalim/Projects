# supporting source code for LDG analyses
require(rgdal)
require(rgeos) # method for over
require(maptools) # for proj4string

m <- readOGR(dsn="/Users/jmarcot/Dropbox/code/R/ldg/dat/mammal_ranges-areas", layer="mammal_ranges-areas")

makeIntervals <- function(startDate, endDate, intervalLength=1) {
	# function begins intervals at start date, and includes an interval >= end date
	# this function works back in time, so if the start date is greater than the end date, it swaps them
	if (startDate>endDate) {
		holder<-startDate
		startDate<-endDate
		endDate<-holder
	}
	tops<-seq(from=startDate, to=endDate, by=intervalLength)
	mat<-data.frame(tops, tops+intervalLength)
	colnames(mat)<-c("ageTop", "ageBase")
	rownames(mat)<-paste(rowMeans(mat), "Ma")
	mat
}

appendTaxonNames <- function(occs, taxonomic.level=c("family", "subfamily", "tribe", "genus", "species"), keep.indet=FALSE, use.original=FALSE, thisSep=" ") {
	taxonomic.level <- match.arg(taxonomic.level)
	if (tolower(taxonomic.level)=="species") {
		if (!keep.indet) occs <- occs[occs$occurrence.species_name != "sp.",]
		occs$taxon <- factor(paste(occs$occurrence.genus_name, occs$occurrence.species_name, sep=thisSep))
	} else if (tolower(taxonomic.level)=="genus") { occs$taxon <- factor(occs$occurrence.genus_name) # for genus level 
	} else if (tolower(taxonomic.level)=="tribe") { occs$taxon <- factor(occs$Tribe) # for genus level 
	} else if (tolower(taxonomic.level)=="subfamily") { occs$taxon <- factor(occs$Subfamily) # for genus level 
	} else if (tolower(taxonomic.level)=="family") occs$taxon <- factor(occs$family_name) # for family level
	if (keep.indet) {
		levels(occs$taxon) <- c(levels(occs$taxon), paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
		occs$taxon[is.na(occs$taxon)] <- factor(paste("[taxon indet. ", seq_along(occs$taxon[occs$taxon==""]), "]", sep=""))
	} else occs <- occs[!is.na(occs$taxon),]
	occs
}

getOccurrenceDatesFromOccs <- function(occs, age.determination=c("midpoint", "random")) {
	age.determination <- match.arg(age.determination)
	thisCols <- getCollectionDatesFromOccs(occs=occs[, c("collection_no", "ma_max", "ma_min")], age.determination=age.determination)
	thisCols[match(occs$collection_no, as.numeric(rownames(thisCols)))]
}

getCollectionDatesFromOccs <- function(occs, age.determination=c("midpoint", "random")) {
	age.determination <- match.arg(age.determination)
	cols <- unique(occs[, c("collection_no", "ma_max", "ma_min")])
	if (age.determination=="random") thisCols <- cbind(cols[,1], apply(cols[,2:3], 1, function(x) { runif(1, max=x[1], min=x[2]) })) else thisCols <- cbind(cols[,1], rowMeans(cols[,2:3]))
	array(thisCols[,2], dimnames=list(thisCols[,1]))
}

buildLatBands<-function(bandWidth, northLat, southLat) {
	x <- seq(from=northLat, to=southLat, by=-bandWidth)
	cbind(bottom=x[2:length(x)], top=x[1:(length(x)-1)])
}

doOneBox<-function(band, thisOccs, thisInterval, thisVec, do.count=TRUE, thisQuota=0.0, use.paleolat=TRUE) {
	if (use.paleolat) { z <- thisVec[thisOccs$ma_rand >= thisInterval["ageTop"] & thisOccs$ma_rand < thisInterval["ageBase"] & thisOccs$gp_mid_lat > band[1] & thisOccs$gp_mid_lat <=band[2]]
	} else z <- thisVec[thisOccs$ma_rand >= thisInterval["ageTop"] & thisOccs$ma_rand < thisInterval["ageBase"] & thisOccs$latdec > band[1] & thisOccs$latdec<=band[2]]
	if (thisQuota > 0) {
		if (thisQuota > length(z)) { return(NA)	# below quota, so not subsampled, no value returned
		} else z <- sample(z, thisQuota)
	}
	if (do.count) return(length(unique(z))) else return(unique(z))
}

doOneInterval<-function(intv, thisQuota=0.0, thisOccs, intervals, bands, thisVec, do.count, use.paleolat=TRUE) {
	apply(X=bands, MARGIN=1, FUN=doOneBox, thisOccs=thisOccs, thisInterval=data.matrix(intervals)[intv,], thisVec=thisVec, do.count=do.count, thisQuota=thisQuota, use.paleolat=use.paleolat)
}

fillBookendsOneTaxon<-function(thisTaxon, thisSample) {
	z <- sapply(thisSample, function(x) thisTaxon %in% x)
	z[min(which(z)):max(which(z))]<-TRUE
	z
}

fillBookendsOneIntv <- function(thisSample) {
	unsampledIntvs <- sapply(thisSample, function(x) all(is.na(x)))
	if (sum(!unsampledIntvs) > 1) {
		# taxList <- unique(unlist(thisSample)[!is.na(unlist(thisSample))])
		taxList <- unique(unlist(lapply(thisSample, function(x) lapply(x, as.character))))
		taxList <- taxList[!is.na(taxList)]
		richness <- rowSums(sapply(taxList, fillBookendsOneTaxon, thisSample))
		richness[unsampledIntvs] <- NA
		richness
	} else array(NA, dim=length(thisSample))
}

getSample <- function(thisOccs, intervals, bands, column=c("cols", "occs", "taxa"), do.count=TRUE, quota=rep(0.0, nrow(intervals)), doBookend=TRUE, use.paleolat=TRUE, do.parallel=FALSE) {
	column <- match.arg(column)
	if (column=="cols") { thisVec <- thisOccs$collection_no
	} else if (column=="occs") { thisVec <- thisOccs$occurrence_no
	} else if (column=="taxa") {
		thisVec <- thisOccs$taxon
		if (doBookend) do.count <- FALSE
	}

	if (do.parallel) { subsample <- mcmapply(FUN=doOneInterval, intv=seq_len(nrow(intervals)), thisQuota=quota, MoreArgs = list(thisOccs=thisOccs, intervals=intervals, bands=bands, thisVec=thisVec, do.count=do.count, use.paleolat=use.paleolat), SIMPLIFY=do.count, mc.cores=detectCores()-2)
	} else subsample <- mapply(FUN=doOneInterval, intv=seq_len(nrow(intervals)), thisQuota=quota, MoreArgs = list(thisOccs=thisOccs, intervals=intervals, bands=bands, thisVec=thisVec, do.count=do.count, use.paleolat=use.paleolat), SIMPLIFY=do.count)
	if (column=="taxa") {
		if (doBookend) {
			if (do.parallel) { mclapply(subsample, fillBookendsOneIntv, mc.cores=detectCores()-2)
			} else vapply(X=subsample, FUN=fillBookendsOneIntv, FUN.VALUE=rep(1, nrow(bands)))
		} else if (do.count) sapply(subsample, function(x) sapply(x, length))
	} else subsample
}

getLongsOneInterval<-function(intOccs, thisOccs) {
	thisLongs <- sapply(intOccs, function(x) thisOccs$lngdec[thisOccs$occurrence_no %in% x], USE.NAMES=TRUE)
	thisPaleoLongs <- sapply(intOccs, function(x) thisOccs$gp_mid_lng[thisOccs$occurrence_no %in% x], USE.NAMES=TRUE)
	cbind( t(sapply(thisLongs, function(x) if (length(x)>0) c(long_min=min(x, na.rm=TRUE), long_max=max(x, na.rm=TRUE)) else c(long_min=NA, long_max=NA) )),
	t(sapply(thisPaleoLongs, function(x) if (length(x)>0) c(paleolong_min =min(x, na.rm=TRUE), paleolong_max =max(x, na.rm=TRUE)) else c(paleolong_min =NA, paleolong_max =NA) )))
}

makeBoxFromBands <- function(thisBand, thisLongBoundary, scalar=0.1, thisProj=proj4string(m)) {
	if (any(!is.finite(c(thisBand, thisLongBoundary)))) return(NA)
	west <- thisLongBoundary[1]
	east <- thisLongBoundary[2]
	north <- thisBand[2]
	south <- thisBand[1]

	ew <- seq(from=west, to=east, by=scalar)
	ns <- seq(from=south, to=north, by=scalar)

	thisBox <- cbind(c(ew, rep(east, length(ns)), rev(ew), rep(west, length(ns))), c(rep(south, length(ew)), ns, rep(north, length(ew)), rev(ns)))
	thisBox <- project(thisBox, proj=thisProj)
	SpatialPolygons(list(Polygons(list(Polygon(thisBox, hole=FALSE)), 1)), proj4string=CRS(thisProj))
}

getExtantRichnessFromOneBoxSet <- function(boxSet) {
	rez <- sapply(boxSet, function(x) if (class(x)=="SpatialPolygons") over(m, x) else rep(NA, nrow(m)))
	spVec <- vector()
	for (i in seq_len(nrow(m))) spVec[i] <- as.character(m[i,]@data$SPP_NAME)
	rownames(rez) <- spVec
	colSums(is.finite(rez))
}

getExtantRichness <- function(longBoundaryMat, bands) {
	boxSet <- list()
	for (i in seq_len(nrow(bands))) boxSet[[i]] <- makeBoxFromBands(thisBand=bands[i,], thisLongBoundary=longBoundaryMat[i,c("paleolong_min", "paleolong_max")])
	getExtantRichnessFromOneBoxSet(boxSet)
}

getFossilSlopeFromOneBandWidth <- function(thisBandWidth, northLat, southLat, intervals, min.lat.samples=3, minLatSpan=10.0, doBookend=TRUE, globalQuota=FALSE, use.paleolat=TRUE, write.Box.Files=FALSE, do.subsample=TRUE, do.extant.richness=TRUE, do.parallel=FALSE) {
	cat("\\t",intervals[1,2]-intervals[1,1]," Myr intervals\\t", round(thisBandWidth, digits=2),"ยบ latitudinal bands\\r", sep="") #"\\tRep ",r,"of",reps,
	bands <- buildLatBands(thisBandWidth, northLat, southLat)

	thisOccs <- cbind(occs, ma_rand=getOccurrenceDatesFromOccs(occs, age.determination))
	thisOccs <- thisOccs[is.finite(thisOccs$ma_rand) & thisOccs$ma_rand < max(intervals) & thisOccs$ma_rand > min(intervals),]


	if (do.subsample) {
		taxSample <- getSample(thisOccs=thisOccs, intervals=intervals, bands=bands, column="taxa", do.count=TRUE, doBookend=doBookend, use.paleolat=use.paleolat, do.parallel=do.parallel)
		if (globalQuota) { quota <- max(taxSample, na.rm=TRUE)
		} else quota <- apply(taxSample, 2, max, na.rm=TRUE)
		# taxSample <- getSample(thisOccs=thisOccs, intervals=intervals, bands=bands, column="taxa", do.count=TRUE, quota=quota, doBookend, use.paleolat, do.parallel=do.parallel)
	} else quota <- rep(0.0, nrow(intervals))

	occSample <- getSample(thisOccs, intervals, bands, column="occs", do.count=FALSE, quota=quota, doBookend, use.paleolat= use.paleolat, do.parallel=do.parallel)
	fossilRichness <- sapply(occSample, function(y) sapply(y, function(x) length(unique(thisOccs$taxon[thisOccs$occurrence_no %in% x]))))
	dimnames(fossilRichness) <- list(rowMeans(bands), rowMeans(intervals))
	fossilRichness[fossilRichness==0] <- NA

	extantRichness <- NA
	if (do.extant.richness) {
		bandLongs <- lapply(occSample, getLongsOneInterval, thisOccs)

		for (i in seq_along(bandLongs)) getExtantRichness(bandLongs[[i]], bands)
		extantRichness <- sapply(bandLongs, getExtantRichness, bands=bands)
		dimnames(extantRichness) <- list(rowMeans(bands), rowMeans(intervals))
		extantRichness[extantRichness==0] <- NA
	}

	list(fossilRichness=fossilRichness, extantRichness=extantRichness)

}

getOneBoxSetFromBandMat <- function(bandMat, longBoundary, scalar=0.1, thisProj=proj4string(m)) {
	lapply(listifyMatrixByRow(bandMat), makeBoxFromBands, thisLongBoundary=longBoundary, scalar=scalar, thisProj=thisProj)
}

getRegFromOneInterval <- function(richVec, lats, min.lat.samples=3, min.lat.span=10) {
	if (sum(is.finite(richVec)) >= min.lat.samples & diff(range(lats[is.finite(richVec)])) >= min.lat.span) {
		glm(log(richVec[is.finite(richVec)]) ~ lats[is.finite(richVec)])
	} else NA
}

getRegFromOneMatrix <- function(richMat, min.lat.samples=3, min.lat.span=10) {
	apply(richMat, 2, getRegFromOneInterval, lats=as.numeric(rownames(richMat)))
}


getModernRegList <- function(thisReps=100, northLat, southLat, westLong, eastLong) {
	bandMatList <- lapply(seq(from=1, to=5, length.out=thisReps), buildLatBands, northLat, southLat)
	boxSetList <- lapply(bandMatList, getOneBoxSetFromBandMat, longBoundary =c(westLong, eastLong))

	richList <- sapply(boxSetList, getExtantRichnessFromOneBoxSet)
	for (i in seq_along(richList)) names(richList[[i]]) <- rowMeans(bandMatList[[i]])
	regList <- lapply(richList, function(x) glm(log(x) ~ as.numeric(names(x))))
	
	regList
}

getExantSlopeVecOneRep <- function(x) {
	sapply(x$extant, function(y) if (!is.na(y)) y$coefficients[2] else NA)
}

makeFossilCIMatOneRep <- function(x) {
	t(sapply(x$fossil, function(y) if (!is.na(y)) { if (var(y$residuals) > 0) confint(y)[2,] else rep(NA, 2) } else rep(NA, 2)))
}

getSigOneRep <- function(x) {
	t(apply(X=cbind(getExantSlopeVecOneRep(x), makeFossilCIMatOneRep(x)), MARGIN=1, FUN=function(z) c(extant= z[1] < z[2] | z[1] > z[3], zero=0.0 < z[2] | 0.0 > z[3])))
}

listifyMatrixByRow <- function(m) {
	intList<-list()
	for (i in 1:nrow(m)) intList[[i]] <- data.matrix(m)[i,]
	intList
}

getContiguousIntv <- function(intv) {
  thisList <- list()
  i <- 1
  thisVec <- vector()
  while (i <= length(intv)) {
    if (length(thisVec) == 0) { thisVec <- intv[i]
    } else {
      if (thisVec[length(thisVec)] == intv[i] - 1) { thisVec <- c(thisVec, intv[i])
      } else {
        thisList[[length(thisList) + 1]] <- thisVec
        thisVec <- intv[i]
      }
    }
    i = i+1
  }
  thisList[[length(thisList) + 1]] <- thisVec
  thisList
}


overlayCzTimescale <- function(do.subepochs=FALSE, color=TRUE, thisAlpha=0.33, borderCol="white", invertTime=FALSE, scale.cex=0.75) {
	textCol <- rgb(0,0,0,thisAlpha)
	# textShadowCol<-"gray50"
	old.cex<-par()$cex
	par(cex=old.cex * scale.cex)
	epochs=data.frame(name=c("Camb",
							"eO",
							"mO",
							"lO",
							"eS",
							"mS",
							"lS",
							"eD",
							"mD",
							"lD",
							"Miss",
							"Penn",
							"eP",
							"mP",
							"lP",
							"eT", 
							"mT",
							"lT",
							"eJ",
							"mJ",
							"lJ",
							"eK",
							"lK",
							"Paleo",
							"Eo",
							"Oligo",
							"Mio",
							"Plio",
							"Pl.",
							"Recent"),
					ageBase =c(542,
								486,
								472,
								461,
								444,
								428,
								423,
								416,
								398,
								385,
								359,
								318,
								299,
								271,
								260,

								251.0,
								245.0,
								235.0,
								201.6,
								176.0,
								161.0,
								145.5,
								99.6,
								65.5,
								55.8,
								33.9,
								23.03,
								5.33,
								2.58,
								0))
	if (color) { epochs<-data.frame(epochs, rgb =c(
		rgb(0.533, 0.671, 0.494, thisAlpha),
		rgb(0.0, 0.686, 0.565, thisAlpha),
		rgb(0.118, 0.737, 0.624, thisAlpha),
		rgb(0.459, 0.796, 0.690, thisAlpha),
		rgb(0.565, 0.835, 0.788, thisAlpha),
		rgb(0.675, 0.871, 0.831, thisAlpha),
		rgb(0.725, 0.894, 0.867, thisAlpha),
		rgb(0.906, 0.698, 0.471, thisAlpha),
		rgb(0.953, 0.788, 0.557, thisAlpha),
		rgb(0.953, 0.875, 0.710, thisAlpha),
		rgb(0.427, 0.624, 0.533, thisAlpha),
		rgb(0.584, 0.769, 0.780, thisAlpha),
		rgb(0.937, 0.463, 0.404, thisAlpha),
		rgb(0.988, 0.549, 0.478, thisAlpha),
		rgb(0.996, 0.702, 0.647, thisAlpha),

		rgb(0.643, 0.365, 0.627, thisAlpha),
		rgb(0.718, 0.510, 0.71, thisAlpha),
		rgb(0.745, 0.616, 0.776, thisAlpha),
		rgb(0.0, 0.718, 0.906, thisAlpha),
		rgb(0.392, 0.816, 0.918, thisAlpha),
		rgb(0.647, 0.882, 0.973, thisAlpha),
		rgb(0.5803922, 0.7960784, 0.4745098, thisAlpha),
		rgb(0.7803922, 0.8784314, 0.6156863, thisAlpha),
		rgb(0.9803922, 0.6980392, 0.4862745, thisAlpha),
		rgb(0.9843137, 0.7372549, 0.5294118, thisAlpha),
		rgb(0.9960784, 0.8588235, 0.6745098, thisAlpha),
		rgb(1, 0.945098, 0, thisAlpha),
		rgb(1, 0.9764706, 0.6823529, thisAlpha),
		rgb(1, 0.9411765, 0.7490196, thisAlpha),
		rgb(1, 0.9529412, 0.9333333, thisAlpha)), stringsAsFactors = FALSE) 
	} else { epochs<-data.frame(epochs, rgb =c(
		rgb(0.57, 0.57, 0.57, thisAlpha),
		rgb(0.51, 0.51, 0.51, thisAlpha),
		rgb(0.59, 0.59, 0.59, thisAlpha),
		rgb(0.68, 0.68, 0.68, thisAlpha),
		rgb(0.79, 0.79, 0.79, thisAlpha),
		rgb(0.79, 0.79, 0.79, thisAlpha),
		rgb(0.79, 0.79, 0.79, thisAlpha),
		rgb(0.69, 0.69, 0.69, thisAlpha),
		rgb(0.77, 0.77, 0.77, thisAlpha),
		rgb(0.85, 0.85, 0.85, thisAlpha),
		rgb(0.51, 0.51, 0.51, thisAlpha),
		rgb(0.68, 0.68, 0.68, thisAlpha),
		rgb(0.53, 0.53, 0.53, thisAlpha),
		rgb(0.61, 0.61, 0.61, thisAlpha),
		rgb(0.73, 0.73, 0.73, thisAlpha),

		rgb(0.39, 0.39, 0.39, thisAlpha),
		rgb(0.51, 0.51, 0.51, thisAlpha),
		rgb(0.59, 0.59, 0.59, thisAlpha),
		rgb(0.58, 0.58, 0.58, thisAlpha),
		rgb(0.71, 0.71, 0.71, thisAlpha),
		rgb(0.81, 0.81, 0.81, thisAlpha),
		rgb(0.69, 0.69, 0.69, thisAlpha),
		rgb(0.81, 0.81, 0.81, thisAlpha),
		rgb(0.71, 0.71, 0.71, thisAlpha),
		rgb(0.75, 0.75, 0.75, thisAlpha),
		rgb(0.85, 0.85, 0.85, thisAlpha),
		rgb(0.93, 0.93, 0.93, thisAlpha),
		rgb(0.96, 0.96, 0.96, thisAlpha),
		rgb(0.93, 0.93, 0.93, thisAlpha),
		rgb(1, 1, 1, thisAlpha)), stringsAsFactors = FALSE) }
	# } else { epochs <- data.frame(epochs, rgb = rep("000000", times=nrow(epochs)), stringsAsFactors = FALSE) }
	
	if (do.subepochs) subepochs=data.frame(modifier=c("e", "m", "l", "e", "m", "l", "e", "l", "e", "m", "l", "PP"), ageBase=c(65.5, 61.7, 58.7, 55.8, 48.6, 37.2, 33.9, 28.4, 23.03, 15.97, 11.61, 5.33), stringsAsFactors = FALSE)

	if (invertTime) {
		epochs$ageBase<- par()$usr[2]/1.137059-epochs$ageBase
		if (do.subepochs) subepochs$ageBase<- par()$usr[2]/1.137059-subepochs$ageBase
	}

	top=par()$usr[4]
	bottom=par()$usr[3]

#lays down colors
	for (i in 1:(nrow(epochs)-1)) {
		polygon(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom), col=epochs$rgb[i], border=borderCol, lty=0)
	}

#lays down Era boundaries
	abline(v=251,lwd=1.5, col=borderCol)
	abline(v=65.5,lwd=1.5, col=borderCol)

	if (do.subepochs) {
		top=(0.95*(top - bottom)) + bottom
		for (i in 1:(nrow(subepochs)-1)) {
			polygon(cbind(c(subepochs[i,2], subepochs[i,2], subepochs[(i+1),2], subepochs[(i+1),2]), c(bottom, top, top, bottom)), lwd=0.25, border=borderCol)
			text (x=mean(c(as.numeric(subepochs[i,2]), as.numeric(subepochs[(i+1),2]))), y=((0.975*(top - bottom)) + bottom), label=subepochs[i,1], col=textCol)
		}
	}

#lays down borders and text
	top=par()$usr[4]
	for (i in 1:(nrow(epochs)-1)) {
		polygon(cbind(c(epochs[i,2], epochs[i,2], epochs[(i+1),2], epochs[(i+1),2]), c(bottom, top, top, bottom)), col=NA, border=borderCol, lwd=0.5)
		# text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, adj=c(0.475,0.525), font=2) #cex=par()$cex+0.05
		# text (mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), (0.975*top), label=epochs[i,1], col=textShadowCol, font=2) #cex=par()$cex+0.05
		text (x=mean(c(as.numeric(epochs[i,2]), as.numeric(epochs[i+1,2]))), y=((0.975*(top - bottom)) + bottom), label=epochs[i,1], col=textCol, font=2)
	}
	abline(h=((0.95*(top - bottom)) + bottom), lwd=1.0, col=borderCol)

	par(cex=old.cex)
}

getAlroyStatisticsOneInterval <- function(thisInt, zachos2001) {
	intTopes <- which(is.finite(zachos2001$Age) & is.finite(zachos2001$d18Oadj) & zachos2001$Age >= thisInt["ageTop"] & zachos2001$Age < thisInt["ageBase"])
	thisLM <- lm(formula=zachos2001$d18Oadj[intTopes] ~ zachos2001$Age[intTopes])
	zachos2001[intTopes,]
	# plot(x=zachos2001$Age[intTopes], y=zachos2001$d18Oadj[intTopes])
	# abline(thisLM)

	# mean(zachos2001$d18Oadj[intTopes])
	c(n=length(intTopes), midpoint=(thisLM$coefficients[2] * mean(thisInt)) + thisLM$coefficients[1], stdev=sd(thisLM$residuals))
}

getAlroyStatistics <- function(intervals) {
	dat <- t(apply(intervals, 1, getAlroyStatisticsOneInterval,zachos2001))
  data.frame(n=dat[,1], midpoint=dat[,2], change=c(diff(dat[,2]), NA), volatility=c(abs(diff(dat[,2])), NA), stdev=dat[,3]) #det.stdev=c(diff(dat[,3]), NA)
}
