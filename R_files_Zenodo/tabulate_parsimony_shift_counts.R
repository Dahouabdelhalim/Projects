### make empty shift chart
shift.chart <- data.frame("method"=vector(mode="character",length=5),"tissue"=vector(mode="character",length=5),"H_old"=vector(mode="integer",length=5),"H_young"=vector(mode="integer",length=5),"G_old"=vector(mode="integer",length=5),"G_young"=vector(mode="integer",length=5),stringsAsFactors=FALSE)
shift.chart[,"method"] <- c("poly.DEL","poly.DEL","poly.ACC","pred.DEL","pred.ACC")
shift.chart[,"tissue"] <- c("none","poly.DEL","poly.ACC","pred.DEL","pred.ACC")

### fill empty shift chart with parsimony shift data
breakpoint <- 3 # 3 million years breakpoint
min_age <- 1 # ignore splits < this age (set at 1 Ma)
for (row_n in 1:5) {
	hcol <- paste("host.",shift.chart[row_n,"method"],sep="")
	gcol <- paste("geo.",shift.chart[row_n,"method"],sep="")
	tcol <- paste("tissue.",shift.chart[row_n,"tissue"],sep="")
	shift.chart[row_n,"H_old"] <- length(which(pars.shift.counts[,"age"] > breakpoint & (pars.shift.counts[,hcol] | pars.shift.counts[,tcol])))/length(which(pars.shift.counts[,"age"] > breakpoint & !is.na(pars.shift.counts[,hcol])))
	shift.chart[row_n,"H_young"] <- length(which((pars.shift.counts[,"age"] < breakpoint & pars.shift.counts[,"age"] > min_age) & (pars.shift.counts[,hcol] | pars.shift.counts[,tcol])))/length(which((pars.shift.counts[,"age"] < breakpoint & pars.shift.counts[,"age"] > min_age) & !is.na(pars.shift.counts[,hcol])))
	shift.chart[row_n,"G_old"] <- length(which(pars.shift.counts[,"age"] > breakpoint & pars.shift.counts[,gcol]))/length(which(pars.shift.counts[,"age"] > breakpoint))
	shift.chart[row_n,"G_young"] <- length(which((pars.shift.counts[,"age"] < breakpoint & pars.shift.counts[,"age"] > min_age) & pars.shift.counts[,gcol]))/length(which(pars.shift.counts[,"age"] < breakpoint & pars.shift.counts[,"age"] > min_age))
	}

row_n <- 2  ##  use this comparison (row 2 = polymorphic coding, DELTRAN optimization)
freq.shifts <- matrix(c(shift.chart[row_n,"H_old"],shift.chart[row_n,"H_young"],shift.chart[row_n,"G_old"],shift.chart[row_n,"G_young"]),ncol=2,nrow=2)
dimnames(freq.shifts) <- list(c("old","young"),c("host","geography"))

barplot(freq.shifts,beside=TRUE,names.arg=c(paste(">",breakpoint,"Ma,host",sep=""),paste(min_age,"-",breakpoint,"Ma,host",sep=""),paste(">",breakpoint,"Ma,region",sep=""),paste(min_age,"-",breakpoint,"Ma,region",sep="")),main="Proportion of Splits Associated With Host and Region Differences")

### repeat shift chart by clade
nall <- 116:229
nA <- 183:229
nB <- 151:182
nC <- 124:147
nodes <- list(nall,nA,nB,nC)
shift.chart <- data.frame(clade=c("all","A","B","C"),"H_old"=vector(mode="integer",length=4),"H_young"=vector(mode="integer",length=4),"G_old"=vector(mode="integer",length=4),"G_young"=vector(mode="integer",length=4),stringsAsFactors=FALSE)
breakpoint <- 3 # 3 million years breakpoint

for (row_n in 1:nrow(shift.chart)) {
	hcol <- "host.poly.DEL"
	gcol <- "geo.poly.DEL"
	tcol <- "tissue.poly.DEL"
	shift.subset <- pars.shift.counts[which(pars.shift.counts[,"node"] %in% nodes[[row_n]]),c("node","age",hcol,gcol,tcol)]
	names(shift.subset) <- c("node","age","host","geo","tissue")
	shift.subset.old <- shift.subset[which(shift.subset[,"age"] > breakpoint),]
	shift.subset.young <- shift.subset[which(shift.subset[,"age"] < breakpoint & shift.subset[,"age"] > min_age),]
	
	shift.chart[row_n,"H_old"] <- length(which(shift.subset.old[,"host"] | shift.subset.old[,"tissue"])) / length(which(!is.na(shift.subset.old[,"host"])))
	shift.chart[row_n,"H_young"] <- length(which(shift.subset.young[,"host"] | shift.subset.young[,"tissue"])) / length(which(!is.na(shift.subset.young[,"host"])))
	shift.chart[row_n,"G_old"] <- length(which(shift.subset.old[,"geo"])) / nrow(shift.subset.old)
	shift.chart[row_n,"G_young"] <- length(which(shift.subset.young[,"geo"])) / nrow(shift.subset.young)
	}

freq.shifts <- list()
freq.shifts$all <- matrix(c(shift.chart[1,"H_old"],shift.chart[1,"H_young"],shift.chart[1,"G_old"],shift.chart[1,"G_young"]),ncol=2,nrow=2)
freq.shifts$A <- matrix(c(shift.chart[2,"H_old"],shift.chart[2,"H_young"],shift.chart[2,"G_old"],shift.chart[2,"G_young"]),ncol=2,nrow=2)
freq.shifts$B <- matrix(c(shift.chart[3,"H_old"],shift.chart[3,"H_young"],shift.chart[3,"G_old"],shift.chart[3,"G_young"]),ncol=2,nrow=2)
freq.shifts$C <- matrix(c(shift.chart[4,"H_old"],shift.chart[3,"H_young"],shift.chart[4,"G_old"],shift.chart[4,"G_young"]),ncol=2,nrow=2)
dimnames(freq.shifts$all) <- list(c("old","young"),c("host","geography"))
dimnames(freq.shifts$A) <- list(c("old","young"),c("host","geography"))
dimnames(freq.shifts$B) <- list(c("old","young"),c("host","geography"))
dimnames(freq.shifts$C) <- list(c("old","young"),c("host","geography"))

layout(matrix(c(1,2,3),nrow=1,ncol=3))
barplot(freq.shifts$A,beside=TRUE,names.arg=c(paste(">",breakpoint,"Ma,host",sep=""),paste(min_age,"-",breakpoint,"Ma,host",sep=""),paste(">",breakpoint,"Ma,region",sep=""),paste("<",breakpoint,"Ma,region",sep="")),main="Clade A",ylim=c(0,1))
barplot(freq.shifts$B,beside=TRUE,names.arg=c(paste(">",breakpoint,"Ma,host",sep=""),paste(min_age,"-",breakpoint,"Ma,host",sep=""),paste(">",breakpoint,"Ma,region",sep=""),paste("<",breakpoint,"Ma,region",sep="")),main="Clade B",ylim=c(0,1))
barplot(freq.shifts$C,beside=TRUE,names.arg=c(paste(">",breakpoint,"Ma,host",sep=""),paste(min_age,"-",breakpoint,"Ma,host",sep=""),paste(">",breakpoint,"Ma,region",sep=""),paste("<",breakpoint,"Ma,region",sep="")),main="Clade C",ylim=c(0,1))


