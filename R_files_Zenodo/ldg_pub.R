# i did it myself

require(sp)

source("~/Dropbox/code/R/ldg/src/ldgSource_pub.R")

reps <- 10000
repsPerBandwidth <- 1
intvLength <- 2
minAge <- 0
maxAge <- 65.5
northLat <-70			
southLat <-25			
eastLong <- -90			
westLong <- -180		
min.lat.samples <- 3
min.lat.span <- 10
globalQuota <- FALSE
doBookend <- TRUE
includeUnstandardizedIntervals <- FALSE
age.determination <- "random"
taxonomic.level <- "species"
intervals <- makeIntervals(minAge, maxAge, intvLength)

bandWidths <- sort(rep(seq(from=1, to=5, length.out=reps), times=repsPerBandwidth))

do.new.analysis <- TRUE
do.parallel <- TRUE
if (do.parallel) require(parallel)

paleoLatLongs <- TRUE	#TRUE = use Paleo Lat/Longs, FALSE = use Modern Lat/Longs

filename_occs <- "allMammals-20150521-occs.csv"
filename_occs <- "~/Dropbox/code/common_dat/allMammals-20150521-occs.csv"
cat("Getting Occurrences from File (", filename_occs, ")\\n\\n", sep="")
	occs <- read.csv(filename_occs, stringsAsFactors=FALSE)

	occs <- occs[is.finite(occs$ma_max) & is.finite(occs$ma_min),]
	switched <- which(occs$ma_max < occs$ma_min)
	if (length(switched > 0)) {
		holder <- occs[switched,"ma_max"]
		occs[switched,"ma_max"] <- occs[switched,"ma_min"]
		occs[switched,"ma_min"] <- holder
	}
	occs <- occs[occs$ma_min < maxAge & occs$ma_max > minAge,]
	occs <- occs[!(occs$order_name %in% c("Cetacea", "Sirenia", "Demostylia")), ]
	occs <- occs[!(occs$family_name %in% c("Odobenidae", "Otariidae",  "Phocidae")), ]
	occs <- occs[occs$latdec<=northLat & occs$latdec>southLat & occs$lngdec<=eastLong & occs$lngdec>=westLong,]
	occs <- appendTaxonNames(occs, taxonomic.level= taxonomic.level, keep.indet=FALSE)

	if (globalQuota) mes1<-"\\tGlobal\\n" else mes1<-"\\tper-Interval\\n"
	if (doBookend) {
		if (includeUnstandardizedIntervals) mes2 <- " and including unstandardized intervals to estimate bookend richness\\n" else mes2<-", but only from subsampled intervals\\n"
		mes2 <- paste("\\tAdding 'bookend' richness", mes2)
	} else mes2<-("\\tUsing only (sub)sampled richness within intervals\\n")
	cat("\\tRequested Reps:\\t\\t",reps,"\\n\\tInterval Length:\\t",intvLength,"\\n\\tLatitude range:\\t\\t",southLat,"º - ",northLat,"º (obs: ", min(occs$gp_mid_lat),"º - ", max(occs$gp_mid_lat),"º)\\n\\tLongitude range:\\t",westLong,"º - ",eastLong,"º (obs: ", min(occs$lngdec),"º - ", max(occs$lngdec),"º)\\n\\tQuota: ", mes1, mes2, sep="")
	cat("\\tNo. Occs:\\t\\t", nrow(occs), "\\n\\tLatitudinal bands range from ", min(bandWidths)," to ", max(bandWidths), "\\t mean: ", mean(bandWidths), "\\tmedian: ", median(bandWidths), "\\n\\ttaxonomic level:\\t", taxonomic.level,"\\n\\n", sep="")

if (do.new.analysis){
	if (do.parallel) { rez <- mclapply(bandWidths, getFossilSlopeFromOneBandWidth, northLat=northLat, southLat=southLat, intervals=intervals, globalQuota=globalQuota, doBookend=doBookend, use.paleolat=TRUE, write.Box.Files=FALSE, do.subsample=TRUE, do.extant.richness=TRUE, do.parallel=FALSE, mc.cores=detectCores()-2)
	} else rez <- 			    lapply(bandWidths, getFossilSlopeFromOneBandWidth, northLat=northLat, southLat=southLat, intervals=intervals, globalQuota=globalQuota, doBookend=doBookend, use.paleolat=TRUE, write.Box.Files=FALSE, do.subsample=TRUE, do.extant.richness=TRUE, do.parallel=FALSE)
} else load("ldg_rez_2My_10000_20151023.R")

#modern LDG
regList <- getModernRegList(thisReps=100, northLat, southLat, westLong, eastLong)

if (do.parallel) { rez.fv <- mclapply(sort(rep(seq(from=1, to=5, length.out=reps/10), times=repsPerBandwidth)), getFossilSlopeFromOneBandWidth, northLat=northLat, southLat=southLat, intervals=intervals, globalQuota=globalQuota, doBookend=doBookend, use.paleolat=TRUE, write.Box.Files=TRUE, do.subsample=FALSE,  do.extant.richness=FALSE, do.parallel=FALSE, mc.cores=detectCores()-2)
} else rez.fv <- 			   lapply(sort(rep(seq(from=1, to=5, length.out=reps/10), times=repsPerBandwidth)), getFossilSlopeFromOneBandWidth, northLat=northLat, southLat=southLat, intervals=intervals, globalQuota=globalQuota, doBookend=doBookend, use.paleolat=TRUE, write.Box.Files=TRUE, do.subsample=FALSE,  do.extant.richness=FALSE, do.parallel=FALSE)

reg.fv <- lapply(rez.fv, function(x) getRegFromOneMatrix(x[[1]], min.lat.samples=min.lat.samples, min.lat.span=min.lat.span))
reg.fossil <- lapply(rez, function(x) getRegFromOneMatrix(x[[1]], min.lat.samples=min.lat.samples, min.lat.span=min.lat.span))
worked.mat <- t(!sapply(reg.fossil, is.na))
reg.extant <- list()
for (i in seq_len(nrow(worked.mat))) {
	thisList <- getRegFromOneMatrix(rez[[i]][[2]], min.lat.samples=min.lat.samples, min.lat.span=min.lat.span)
	for (j in seq_len(ncol(worked.mat))) {
		if (!worked.mat[i,j]) thisList[[j]] <- NA
	}
	reg.extant[[i]] <- thisList
}

reg <- list()
for (i in seq_along(reg.fossil)) {
	reg[[i]]	<- list(fossil=reg.fossil[[i]], extant=reg.extant[[i]])
}

slopes <- sapply(reg, function(x) sapply(x, function(y) sapply(y, function(z) if (!is.na(z)) z$coefficients[2] else NA), simplify="array"), simplify="array")
slopes.q <- apply(slopes, c(1,2), quantile, probs=c(0,0.025, 0.25, 0.5, 0.75, 0.975, 1.0), na.rm=TRUE)
slopes.fv <- sapply(reg.fv, function(y) sapply(y, function(x) if (!is.na(x)) x$coefficients[2] else NA))
slopes.fv.q <- apply(slopes.fv, 1, quantile, probs=c(0,0.025, 0.25, 0.5, 0.75, 0.975, 1.0), na.rm=TRUE)

sig <- sapply(reg, FUN= getSigOneRep, simplify="array")

c.1 <- apply(slopes[,"fossil",], 1, function(x) sum(is.finite(x)))
c.4 <- rowSums(apply(slopes, 3, function(thisRep) thisRep[,1] > thisRep[,2]), na.rm=TRUE)/c.1

table1 <- data.frame(c.1, round(cbind(slopes.q["50%",,], c.4, apply(sig, c(1,2), sum, na.rm=TRUE)/apply(sig, c(1,2), function(x) sum(!is.na(x)))), digits=3))
rownames(table1) <- gsub(" ", "_", rownames(intervals))
colnames(table1) <- c("n", "med.fossil", "med.adj.Recent", "p.fossil.weaker", "p.sig", "p.sig.0")
table1


# ### Figure 1
	thisBandMat <- buildLatBands(bandWidth=3.0, northLat=northLat, southLat=southLat)
	thisBoxSet <- getOneBoxSetFromBandMat(bandMat=thisBandMat, longBoundary=c(westLong, eastLong))
	thisRichness <- getExtantRichnessFromOneBoxSet(thisBoxSet)

	pdf(file = "Fig1_extantSlopeBiplot_3º_20160413.pdf", width = 3.42, height = 3.42, pointsize=10)

	par(mar=c(3, 3.5, 0.5, 0.5), mgp=c(2,1,0))
	plot(rowMeans(thisBandMat), log(thisRichness), ylim=c(3.5, 6.0), xlab="Latitude (ºN)", ylab="log-Richness", type="n")
	abline(lm(log(thisRichness) ~ rowMeans(thisBandMat)), lwd=2)
	points(rowMeans(thisBandMat), log(thisRichness), pch=21, col="black", bg="gray50", cex=1.25)
	dev.off()
	
### Figure 2
	pdf(file = "Fig2_fossil_vs_extant_log_10000_20160413.pdf", width = 7, height = 4, pointsize=9)
	par(mar=c(3, 3.5, 0.5, 0.5), mgp=c(2,1,0))
	intv <- which(is.finite(slopes.q[4,,1]))[-1]
	names(intv) <- NULL
	intvList <- getContiguousIntv(intv)
	plot(rowMeans(intervals), slopes.q[4,,1], type="n", xlim=c(65,0), ylim=c(-0.1, 0.1), xlab="Time (Ma)", ylab="Estimated Slope of LDG")
	overlayCzTimescale(do.subepochs=TRUE)
	abline(h=0, lty=3, col="gray50")

	# #range of face-value modern slopes
	polygon(c(-10,90,90,-10), sort(rep(range(sapply(regList, function(x) x$coefficients[2])), 2)), border=adjustcolor("black", alpha.f=0.5), col=adjustcolor("black", alpha.f=0.25))

	# # fossil slope confidence envelope
	for (i in seq_along(intvList))
	  if (length(intvList[[i]]) > 1) {
	    polygon(c(rowMeans(intervals)[intvList[[i]]], rev(rowMeans(intervals)[intvList[[i]]])), c(slopes.q[2,intvList[[i]],1], rev(slopes.q[6,intvList[[i]],1])), border=NA, col=adjustcolor(col="firebrick4", alpha.f=0.25))
	  } else {
	    stitch <- diff(data.matrix(intervals)[intvList[[i]],])/4
	    thisSlopes <- c(slopes.q[2,intvList[[i]],1], slopes.q[6,intvList[[i]],1])
	    polygon(c(rep(rowMeans(intervals)[intvList[[i]]] + stitch,2), rep(rowMeans(intervals)[intvList[[i]]] - stitch, 2)), c(thisSlopes,rev(thisSlopes)), border=NA, col=adjustcolor(col="firebrick4", alpha.f=0.25))
	  }

	# # adjusted modern slope confidence envelope
	for (i in seq_along(intvList))
	  if (length(intvList[[i]]) > 1) {
	    polygon(c(rowMeans(intervals)[intvList[[i]]], rev(rowMeans(intervals)[intvList[[i]]])), c(slopes.q[2,intvList[[i]],2], rev(slopes.q[6,intvList[[i]],2])), border=NA, col=adjustcolor(col="dodgerblue4", alpha.f=0.25))
	  } else {
	    stitch <- diff(data.matrix(intervals)[intvList[[i]],])/4
	    thisSlopes <- c(slopes.q[2,intvList[[i]],2], slopes.q[6,intvList[[i]],2])
	    polygon(c(rep(rowMeans(intervals)[intvList[[i]]] + stitch,2), rep(rowMeans(intervals)[intvList[[i]]] - stitch, 2)), c(thisSlopes,rev(thisSlopes)), border=NA, col=adjustcolor(col="dodgerblue4", alpha.f=0.25))
	  }

	# face-value fossil slopes
	for (i in seq_along(intvList)) points(rowMeans(intervals)[intvList[[i]]], slopes.fv.q[4,intvList[[i]]], type="l", lty=2, lwd=1)
	for (i in seq_along(intvList)) points(rowMeans(intervals)[intvList[[i]]], slopes.fv.q[4,intvList[[i]]], type="p", cex=1, pch=21, col="black")

	for (i in seq_along(intvList)){
		points(rowMeans(intervals)[intvList[[i]]], slopes.q[4,intvList[[i]],1], type="o", lwd=1.5, cex=1, pch=21, col="firebrick4", bg="firebrick1")
		points(rowMeans(intervals)[intvList[[i]]], slopes.q[4,intvList[[i]],2], type="o", lwd=1.5, cex=1, pch=22, col="dodgerblue4", bg="dodgerblue1")
	  }
	box()
	dev.off()

#### Figure 3
	zachos2001 <- read.csv ("zachos2001.csv")
	topes <- getAlroyStatistics(intervals)

	pdf(file = "Fig3_slopeVsTopes_20160413.pdf", width = 3.42, height = 3.42, pointsize=8)
	
	par(mar=c(3, 3.5, 0.5, 0.5), mgp=c(2,1,0))
	x <- -diff(topes$midpoint)[-1]
	y <- -diff(slopes.q["50%",,"fossil"])[-1]
	plot(x, y, xlim=c(-0.2, 0.6), ylim=c(-0.02, 0.06), xlab= expression(paste(Delta, " ", delta^18,plain(O))), ylab=expression(paste(Delta, " median fossil slope")), type="n")
	abline(h=0.0, lty=2, col="gray50")
	abline(v=0.0, lty=2, col="gray50")
	thisReg <- lm(y ~ x)
	abline(thisReg, lwd=2)
	points(x, y, pch=21, bg="gray50", cex=1.5)
	
	thisCor <- cor.test(x, y, method="spearman")
	text(x=0.25, y=0.05, adj=c(0, 0), labels=expression(paste("Spearman's ", rho, " = -0.525")))
	text(x=0.25, y=0.0475, adj=c(0, 0), paste("\\np = ", round(thisCor$p.value, 3)))
	dev.off()

# ###Suppemental Figure 1
	pdf(file = "FigS1_correlationFigure.pdf", width = 9.45, height = 7, pointsize=8)
	x <- -diff(topes$midpoint)[-1]
	y <- -diff(slopes.q["50%",,"fossil"])[-1]
	names(y) <- seq(3, 63, 2)

	quartz(height=6, width=12)
	par(mfrow=c(1,2), mar=c(4,4,1,1))
	plot(x, y, xlim=range(x[!is.na(y)]), xlab= expression(paste(Delta, " ", delta^18,plain(O))), ylab=expression(paste(Delta, " median fossil slope")), pch=21, bg="gray50", cex=1.5)
	abline(h=0.0, lty=2, col="gray50")
	abline(v=0.0, lty=2, col="gray50")
	text(x, y, labels=names(y), cex=0.5, pos=4)
	thisReg <- lm(y ~ x)
	abline(thisReg, lwd=2)
	thisCor <- cor.test(x, y, method="spearman")
	# text(x=min(x[!is.na(y)]), y=min(y, na.rm=TRUE), adj=c(0, 0), labels=paste("Spearman's rho = ", round(thisCor$estimate,3), "\\np = ", round(thisCor$p.value, 3)))
	text(x=0.3, y=0.05, adj=c(0, 0), labels=expression(paste("Spearman's ", rho, " = -0.525")))
	text(x=0.3, y=0.0475, adj=c(0, 0), paste("\\np = ", round(thisCor$p.value, 3)))

	x <- x[is.finite(y)]
	y <- y[is.finite(y)]
	
	corList <- list()
	for (i in seq_along(x)) {
		corList[[i]] <- cor.test(x[-i], y[-i], method="spearman")
	}
	
	t(sapply(corList, function(x) c(x$estimate, x$p.value)))
	barplot(sapply(corList, function(x) x$p.value), names.arg=names(y), main="p-value when intervals are removed", xlab="interval", ylab="p.value")
	abline(h=0.05, lty=3, col="red")
	dev.off()

# ###Suppemental Figure S2
	pdf(file = "FigS2_bandWidthSplatter2_log_10000_20160412.pdf", width = 9.45, height = 7, pointsize=8)
	par(mfrow=rep(ceiling(sqrt(nrow(intervals))), 2), mar=c(4,4,1,1), mgp=c(2,1,0))
	for (i in seq_len(nrow(intervals))) {
		colorVec <- array(adjustcolor(col="black", alpha.f=0.1), dim=reps)
		colorVec[sig[i,1,]] <- adjustcolor(col="dodgerblue4", alpha.f=0.5)
		# symVec <- array(, dim=length(m[i,]))
		plot(bandWidths, slopes[i,2,], ylim=c(min(slopes, na.rm=TRUE), max(slopes, na.rm=TRUE)), type="p", pch=22, cex=0.5, col=colorVec, ylab="slope", main=rownames(intervals)[i])
		colorVec[sig[i,1,]] <- adjustcolor(col="firebrick4", alpha.f=0.5)
		if (all(is.na(slopes[i,1,]))) { rect(100,-100,-100,100, col="gray50", fg="gray50", bg="gray50")
		} else points(bandWidths, slopes[i,1,], type="p", pch=21, col=colorVec, cex=0.5)
	}
	dev.off()
	
#################

