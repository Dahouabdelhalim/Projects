library(sp)
library(raster)

##  function to generate a random sample of size `n' of pixels
##  for potential use as actual plots
##    x: site raster with values outsite the site polygon masked
##    n: number of candidate plots to generate
##    d: minimum distance between plots
##    rng: range of acceptable canopy heights

rand.pxl <- function(x,n,d,rng) {
	out.note <- 'all good'
	# browser()
	vals <- values(x)
	good.cells <- which(vals >= min(rng,na.rm=TRUE) & vals <= max(rng,na.rm=TRUE) & !is.na(vals))
	
	if(length(good.cells) < n) {
		n <- length(good.cels)
		out.note <- 'too few cells in canopy height range'
	}
	
	good.pnts <- xyFromCell(x, good.cells, spatial=TRUE)
	good.pnts.dist <- spDists(good.pnts)
	
	##  pick points >= d distance appart
	pos.pnts <- replicate(100,brut.seq.samp(good.pnts.dist,n,d))
	
	##  see if we got any sets w/ >= n points
	# browser()
	npnts <- sapply(pos.pnts,length)
	if(any(npnts >= n)) {
		out.pnts <- sample(sample(pos.pnts[npnts >= n], 1)[[1]], n)
	} else {
		out.pnts <- pos.pnts[[which.max(npnts)]]
		if(out.note != 'all good') {
			out.note <- paste(out.note, sprintf('could not find enough pixels >= %s m appart', d), sep='; ')
		} else {
			out.note <- sprintf('could not find enough pixels >= %s m appart', d)
		}
		
	}
	
	out <- good.pnts[out.pnts,]
	attr(out,'note') <- out.note
	return(out)
}

##  functions to generate random samples >= d distance appart

##  best option so far--generates truely random sample
##  then purges those points that are too close then tries
##  to add back as many points as possible
brut.seq.samp <- function(pdist,n,d) {
	# browser()
	try1 <- brut.samp(pdist, n*3, d) # was: n*ceiling(nrow(pdist)/n/2)
	
	if(length(try1) < n) {
		pos <- which(apply(pdist[try1,] >= d, 2, all))
		try2 <- c(try1,sample(pos,size=ifelse(length(pos) > n-length(try1), n-length(try1), length(pos))))
		
		return(purge2close(pdist < d, try2))
	} else {
		return(try1)
	}
}

##  just tries to get lucky and find n points >=d appart by chance
brut.samp <- function(pdist,n,d) {
	pnts <- sample(nrow(pdist),n)
	
	purge2close(pdist < d, pnts)
}

##  gets rid of points in too close of proximity
purge2close <- function(plogic, pnts) {
	# browser()
	this.logic <- plogic[pnts,pnts]
	
	if(any(this.logic[lower.tri(this.logic)])) {
		pnts <- pnts[-which.max(rowSums(this.logic))]
		return(purge2close(plogic,pnts))
	} else {
		return(pnts)
	}
}

##  seq.samp is still buggy
seq.samp <- function(pdist,n,d) {
	pnts <- sample(nrow(pdist),1)
	ok2add <- pdist[pnts,] >= d
	for(i in 1:(n-1)) {
		if(any(ok2add)) {
			new.pnt <- sample(which(ok2add),1)
			pnts <- c(pnts,new.pnt)
			ok2add <- pdist[new.pnt,] >= d & ok2add
			print(min(pdist[pnts,pnts][lower.tri(pdist[pnts,pnts])]))
		} else {
			break
		}
	}
	
	return(pnts)
}
