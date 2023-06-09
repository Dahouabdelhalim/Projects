## Function to plot reaction norms
require(viridis)
require(colormap)

#########################################################################################
## Function to plot the community composition 
library('mapplots')
ComComp.plot <- function(dat=mean.PTS09.Pau, CG=TRUE, show.biomass=FALSE){

	if(CG==TRUE){
	
	sal.lvl <- unique(c(dat$Sal.ori, dat$Sal.des))
	sal.lvl <- sal.lvl[order(sal.lvl)]

	for(i in 1:length(sal.lvl)){
		#print(i)

		plot(0:6, rep(-100,7), ylim=c(0,4), main=paste('salinity origin - ', sal.lvl[i]))
		
		for(j in 1:length(sal.lvl)){
			#print(j)

			dat.tmp <- dat[dat$Sal.ori==sal.lvl[i] & dat$Sal.des==sal.lvl[j],]

			if(show.biomass==TRUE){
				frac.tmp <- tapply(dat.tmp$Bio.Fraction, dat.tmp$ID.ori, mean)
			} else {
				frac.tmp <- tapply(dat.tmp$Fraction, dat.tmp$ID.ori, mean)
			}

			if(length(frac.tmp)!=0){
			for(k in 1:length(frac.tmp)){
				add.pie(z=c(frac.tmp[k],1-frac.tmp[k]), x = j, y = k, radius = 0.25)
			}}
		}
	}
	}

	if(CG==FALSE){
	
	sal.lvl <- unique(dat[[1]]$Salinity)
	sal.lvl <- sal.lvl[order(sal.lvl)]

	nn <- length(dat)

	plot(0:6, rep(-100,7), ylim=c(0,4))

	if(nn==1){
	for(i in 1:length(sal.lvl)){
		dat.tmp <- dat[[1]][dat[[1]]$Salinity==sal.lvl[i],]
		if(show.biomass==TRUE){
			frac.tmp <- tapply(dat.tmp$Bio.Fraction, dat.tmp$ID, mean)
		} else {
			frac.tmp <- tapply(dat.tmp$Fraction, dat.tmp$ID, mean)
		}
		if(length(frac.tmp)!=0){
			for(k in 1:length(frac.tmp)){
				add.pie(z=c(frac.tmp[k],1-frac.tmp[k]), x = i, y = k, radius = 0.25)
			}
		}
	}}

	if(nn==2){
	for(i in 1:length(sal.lvl)){
		dat1.tmp <- dat[[1]][dat[[1]]$Salinity==sal.lvl[i],]
		dat2.tmp <- dat[[2]][dat[[2]]$Salinity==sal.lvl[i],]
		if(show.biomass==TRUE){
			frac1.tmp <- tapply(dat1.tmp$Bio.Fraction, dat1.tmp$ID, mean)
			frac2.tmp <- tapply(dat2.tmp$Bio.Fraction, dat2.tmp$ID, mean)
		} else {	
			frac1.tmp <- tapply(dat1.tmp$Fraction, dat1.tmp$ID, mean)
			frac2.tmp <- tapply(dat2.tmp$Fraction, dat2.tmp$ID, mean)
		}
		if(length(frac1.tmp)!=0 & length(frac2.tmp)!=0){
			names.dat1 <- names(frac1.tmp)
			names.dat2 <- names(frac2.tmp)
			names.tmp <- unique(c(names.dat1, names.dat2))

			for(k in 1:length(names.tmp)){
				if(length(which(names.dat1==names.tmp[k]))!=0){z1 <- frac1.tmp[which(names.dat1==names.tmp[k])]} else {z1 <- 0}
				if(length(which(names.dat2==names.tmp[k]))!=0){z2 <- frac2.tmp[which(names.dat2==names.tmp[k])]} else {z2 <- 0}
				z3 <- 1-z1-z2
				add.pie(z=c(z1,z2,z3), x = i, y = k, radius = 0.25)
			}
		}
		if(length(frac1.tmp)!=0 & length(frac2.tmp)==0){
			names.dat1 <- names(frac1.tmp)
			#names.dat2 <- names(frac2.tmp)
			names.tmp <- unique(names.dat1)

			for(k in 1:length(names.tmp)){
				if(length(which(names.dat1==names.tmp[k]))!=0){z1 <- frac1.tmp[k]} else {z1 <- 0}
				#if(length(which(names.dat2==names.tmp[k]))!=0){z2 <- frac2.tmp[k]} else {z2 <- 0}
				z3 <- 1-z1
				add.pie(z=c(z1,0,z3), x = i, y = k, radius = 0.25)
			}
		}
	}}		
	}
}

			











