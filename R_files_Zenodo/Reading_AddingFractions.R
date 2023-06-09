## Adding the fractions to the data set

mean.PTS09.Pau$Fraction <- 0
mean.PTS09.Spite$Fraction <- 0

## Fraction calculated on the biomass
mean.PTS09.Pau$Bio.Fraction <- 0
mean.PTS09.Spite$Bio.Fraction <- 0

id.lvl <- unique(c(mean.PTS09.Pau$ID.new, mean.PTS09.Spite$ID.new))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean.PTS09.Pau$ID.new == id.lvl[i])
	pos.spite <- which(mean.PTS09.Spite$ID.new == id.lvl[i])
	if(length(pos.pau)==1 & length(pos.spite)==1){
		mean.PTS09.Pau$Fraction[pos.pau] <- mean.PTS09.Pau$Density[pos.pau]/(mean.PTS09.Pau$Density[pos.pau]+mean.PTS09.Spite$Density[pos.spite])
		mean.PTS09.Spite$Fraction[pos.spite] <- mean.PTS09.Spite$Density[pos.spite]/(mean.PTS09.Pau$Density[pos.pau]+mean.PTS09.Spite$Density[pos.spite])
	}
	if(length(pos.pau)==1 & length(pos.spite)==0){
		mean.PTS09.Pau$Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1){
		mean.PTS09.Spite$Fraction[pos.spite] <- 1
	}
}

id.lvl <- unique(c(mean.PTS09.Pau$ID.new, mean.PTS09.Spite$ID.new))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean.PTS09.Pau$ID.new == id.lvl[i])
	pos.spite <- which(mean.PTS09.Spite$ID.new == id.lvl[i])	

	pau.tmp <- dat09.PTS.Pau[dat09.PTS.Pau$ID_new == id.lvl[i],]
	pau.biomass.tmp <- sum(pau.tmp$mean_area)
	spite.tmp <- dat09.PTS.Spite[dat09.PTS.Spite$ID_new == id.lvl[i],]
	spite.biomass.tmp <- sum(spite.tmp$mean_area)

	if(length(pos.pau)==1 & length(pos.spite)==1){
		mean.PTS09.Pau$Bio.Fraction[pos.pau] <- pau.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
		mean.PTS09.Spite$Bio.Fraction[pos.spite] <- spite.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
	}
	if(length(pos.pau)==1 & length(pos.spite)==0){
		mean.PTS09.Pau$Bio.Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1){
		mean.PTS09.Spite$Bio.Fraction[pos.spite] <- 1
	}
}

##-----------------------------------------------------------------------------
## Calculate the fraction for the ancestral data at the beginning and end
mean05.PTS.Pau$Fraction <- 0
mean05.PTS.Spite$Fraction <- 0
mean05.PTS.Pau$Bio.Fraction <- 0
mean05.PTS.Spite$Bio.Fraction <- 0

id.lvl <- unique(c(mean05.PTS.Pau$ID, mean05.PTS.Spite$ID))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean05.PTS.Pau$ID == id.lvl[i])
	pos.spite <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.pau)==1 & length(pos.spite)==1){
		mean05.PTS.Pau$Fraction[pos.pau] <- mean05.PTS.Pau$Density[pos.pau]/(mean05.PTS.Pau$Density[pos.pau]+mean05.PTS.Spite$Density[pos.spite])
		mean05.PTS.Spite$Fraction[pos.spite] <- mean05.PTS.Spite$Density[pos.spite]/(mean05.PTS.Pau$Density[pos.pau]+mean05.PTS.Spite$Density[pos.spite])
	}
	if(length(pos.pau)==1 & length(pos.spite)==0){
		mean05.PTS.Pau$Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1){
		mean05.PTS.Spite$Fraction[pos.spite] <- 1
	}
}

id.lvl <- unique(c(mean05.PTS.Pau$ID, mean05.PTS.Spite$ID))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean05.PTS.Pau$ID == id.lvl[i])
	pos.spite <- which(mean05.PTS.Spite$ID == id.lvl[i])	

	pau.tmp <- dat05.PTS.Pau[dat05.PTS.Pau$Sample_ID == id.lvl[i],]
	pau.biomass.tmp <- sum(pau.tmp$mean_area)
	spite.tmp <- dat05.PTS.Spite[dat05.PTS.Spite$Sample_ID == id.lvl[i],]
	spite.biomass.tmp <- sum(spite.tmp$mean_area)

	if(length(pos.pau)==1 & length(pos.spite)==1){
		mean05.PTS.Pau$Bio.Fraction[pos.pau] <- pau.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
		mean05.PTS.Spite$Bio.Fraction[pos.spite] <- spite.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
	}
	if(length(pos.pau)==1 & length(pos.spite)==0){
		mean05.PTS.Pau$Bio.Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1){
		mean05.PTS.Spite$Bio.Fraction[pos.spite] <- 1
	}
}

##-----------------------------------------------------------------------------
mean22.PTS.Pau.CG$Fraction <- 0
mean22.PTS.Spite.CG$Fraction <- 0
mean22.PTS.Tet.CG$Fraction <- 0

id.lvl <- unique(c(mean22.PTS.Pau.CG$ID, mean22.PTS.Spite.CG$ID, mean22.PTS.Tet.CG$ID))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean22.PTS.Pau.CG$ID == id.lvl[i])
	pos.spite <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	pos.tet <- which(mean22.PTS.Tet.CG$ID == id.lvl[i])
	if(length(pos.pau)==1 & length(pos.spite)==1 & length(pos.tet)==1){
		mean22.PTS.Pau.CG$Fraction[pos.pau] <- mean22.PTS.Pau.CG$Density[pos.pau]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Spite.CG$Density[pos.spite]+mean22.PTS.Tet.CG$Density[pos.tet])
		mean22.PTS.Spite.CG$Fraction[pos.spite] <- mean22.PTS.Spite.CG$Density[pos.spite]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Spite.CG$Density[pos.spite]+mean22.PTS.Tet.CG$Density[pos.tet])
		mean22.PTS.Tet.CG$Fraction[pos.tet] <- mean22.PTS.Tet.CG$Density[pos.tet]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Spite.CG$Density[pos.spite]+mean22.PTS.Tet.CG$Density[pos.tet])
	}
	if(length(pos.pau)==1 & length(pos.spite)==0 & length(pos.tet)==1){
		mean22.PTS.Pau.CG$Fraction[pos.pau] <- mean22.PTS.Pau.CG$Density[pos.pau]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Tet.CG$Density[pos.tet])
		mean22.PTS.Tet.CG$Fraction[pos.tet] <- mean22.PTS.Tet.CG$Density[pos.tet]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Tet.CG$Density[pos.tet])
	}
	if(length(pos.pau)==1 & length(pos.spite)==1 & length(pos.tet)==0){
		mean22.PTS.Pau.CG$Fraction[pos.pau] <- mean22.PTS.Pau.CG$Density[pos.pau]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Spite.CG$Density[pos.spite])
		mean22.PTS.Spite.CG$Fraction[pos.spite] <- mean22.PTS.Spite.CG$Density[pos.spite]/(mean22.PTS.Pau.CG$Density[pos.pau]+mean22.PTS.Spite.CG$Density[pos.spite])
	}
	if(length(pos.pau)==0 & length(pos.spite)==1 & length(pos.tet)==1){
		mean22.PTS.Spite.CG$Fraction[pos.spite] <- mean22.PTS.Spite.CG$Density[pos.spite]/(mean22.PTS.Spite.CG$Density[pos.spite]+mean22.PTS.Tet.CG$Density[pos.tet])
		mean22.PTS.Tet.CG$Fraction[pos.tet] <- mean22.PTS.Tet.CG$Density[pos.tet]/(mean22.PTS.Spite.CG$Density[pos.spite]+mean22.PTS.Tet.CG$Density[pos.tet])
	}
	if(length(pos.pau)==1 & length(pos.spite)==0 & length(pos.tet)==0){
		mean22.PTS.Pau.CG$Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1 & length(pos.tet)==0){
		mean22.PTS.Spite.CG$Fraction[pos.spite] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==0 & length(pos.tet)==1){
		mean22.PTS.Tet.CG$Fraction[pos.tet] <- 1
	}
}

mean22.PTS.Pau.CG$Bio.Fraction <- 0
mean22.PTS.Spite.CG$Bio.Fraction <- 0
mean22.PTS.Tet.CG$Bio.Fraction <- 0

id.lvl <- unique(c(mean22.PTS.Pau.CG$ID, mean22.PTS.Spite.CG$ID, mean22.PTS.Tet.CG$ID))
for(i in 1:length(id.lvl)){
	pos.pau <- which(mean22.PTS.Pau.CG$ID == id.lvl[i])
	pos.spite <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	pos.tet <- which(mean22.PTS.Tet.CG$ID == id.lvl[i])

	pau.tmp <- dat22.PTS.Pau.CG[dat22.PTS.Pau.CG$ID == id.lvl[i],]
	pau.biomass.tmp <- sum(pau.tmp$mean_area)
	spite.tmp <- dat22.PTS.Spite.CG[dat22.PTS.Spite.CG$ID == id.lvl[i],]
	spite.biomass.tmp <- sum(spite.tmp$mean_area)
	tet.tmp <- dat22.PTS.Tet.CG[dat22.PTS.Tet.CG$ID == id.lvl[i],]
	tet.biomass.tmp <- sum(tet.tmp$mean_area)

	if(length(pos.pau)==1 & length(pos.spite)==1 & length(pos.tet)==1){
		mean22.PTS.Pau.CG$Bio.Fraction[pos.pau] <- pau.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp+tet.biomass.tmp)
		mean22.PTS.Spite.CG$Bio.Fraction[pos.spite] <- spite.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp+tet.biomass.tmp)
		mean22.PTS.Tet.CG$Bio.Fraction[pos.tet] <- tet.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp+tet.biomass.tmp)
	}
	if(length(pos.pau)==1 & length(pos.spite)==0 & length(pos.tet)==1){
		mean22.PTS.Pau.CG$Bio.Fraction[pos.pau] <- pau.biomass.tmp/(pau.biomass.tmp+tet.biomass.tmp)
		mean22.PTS.Tet.CG$Bio.Fraction[pos.tet] <- tet.biomass.tmp/(pau.biomass.tmp+tet.biomass.tmp)
	}
	if(length(pos.pau)==1 & length(pos.spite)==1 & length(pos.tet)==0){
		mean22.PTS.Pau.CG$Bio.Fraction[pos.pau] <- pau.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
		mean22.PTS.Spite.CG$Bio.Fraction[pos.spite] <- spite.biomass.tmp/(pau.biomass.tmp+spite.biomass.tmp)
	}
	if(length(pos.pau)==0 & length(pos.spite)==1 & length(pos.tet)==1){
		mean22.PTS.Spite.CG$Bio.Fraction[pos.spite] <- spite.biomass.tmp/(spite.biomass.tmp+tet.biomass.tmp)
		mean22.PTS.Tet.CG$Bio.Fraction[pos.tet] <- tet.biomass.tmp/(spite.biomass.tmp+tet.biomass.tmp)
	}
	if(length(pos.pau)==1 & length(pos.spite)==0 & length(pos.tet)==0){
		mean22.PTS.Pau.CG$Bio.Fraction[pos.pau] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==1 & length(pos.tet)==0){
		mean22.PTS.Spite.CG$Bio.Fraction[pos.spite] <- 1
	}
	if(length(pos.pau)==0 & length(pos.spite)==0 & length(pos.tet)==1){
		mean22.PTS.Tet.CG$Bio.Fraction[pos.tet] <- 1
	}
}

##-----------------------------------------------------------------------------
## Now add the biomass fractions to the individual data to then be used in the regression analysis 
dat09.Pau$Bio.Fraction <- 1
dat09.PTS.Pau$Bio.Fraction <- 0

id.lvl <- unique(dat09.PTS.Pau$ID_new)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat09.PTS.Pau$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Pau$ID.new == id.lvl[i])
	dat09.PTS.Pau$Bio.Fraction[pos.ind] <- mean.PTS09.Pau$Bio.Fraction[pos.mean]
}

dat09.Spite$Bio.Fraction <- 1
dat09.PTS.Spite$Bio.Fraction <- 0

id.lvl <- unique(dat09.PTS.Spite$ID_new)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat09.PTS.Spite$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Spite$ID.new == id.lvl[i])
	dat09.PTS.Spite$Bio.Fraction[pos.ind] <- mean.PTS09.Spite$Bio.Fraction[pos.mean]
}

##-----------------------------------------------------------------------------
## Now add the biomass fractions and densities of the other species to the individual data to then be used in the regression analysis 
dat09.Pau$Bio.Fraction.other <- 0
dat09.PTS.Pau$Bio.Fraction.other <- 0

id.lvl <- unique(dat09.PTS.Pau$ID_new)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat09.PTS.Pau$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Spite$ID.new == id.lvl[i])
	if(length(pos.mean)!=0){
		dat09.PTS.Pau$Bio.Fraction.other[pos.ind] <- mean.PTS09.Spite$Bio.Fraction[pos.mean]
	}
}

dat09.Pau$Density.other <- 0
dat09.PTS.Pau$Density.other <- 0

id.lvl <- unique(dat09.PTS.Pau$ID_new)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat09.PTS.Pau$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Spite$ID.new == id.lvl[i])
	if(length(pos.mean)!=0){
		dat09.PTS.Pau$Density.other[pos.ind] <- mean.PTS09.Spite$Density[pos.mean]
	}
}

dat09.Spite$Bio.Fraction.other <- 0
dat09.PTS.Spite$Bio.Fraction.other <- 0

id.lvl <- unique(dat09.PTS.Spite$ID_new)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat09.PTS.Spite$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Pau$ID.new == id.lvl[i])
	if(length(pos.mean)!=0){
		dat09.PTS.Spite$Bio.Fraction.other[pos.ind] <- mean.PTS09.Pau$Bio.Fraction[pos.mean]
	}
}

dat09.Spite$Density.other <- 0
dat09.PTS.Spite$Density.other <- 0

id.lvl <- unique(dat09.PTS.Spite$ID_new)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat09.PTS.Spite$ID_new == id.lvl[i])
	pos.mean <- which(mean.PTS09.Pau$ID.new == id.lvl[i])
	if(length(pos.mean)!=0){
		dat09.PTS.Spite$Density.other[pos.ind] <- mean.PTS09.Pau$Density[pos.mean]
	}
}


dat22.Pau.CG$Bio.Fraction.other <- 0
dat22.PTS.Pau.CG$Bio.Fraction.other <- 0

id.lvl <- unique(dat22.PTS.Pau.CG$ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat22.PTS.Pau.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Pau.CG$Bio.Fraction.other[pos.ind] <- mean22.PTS.Spite.CG$Bio.Fraction[pos.mean]
	}
}

dat22.Pau.CG$Density.other <- 0
dat22.PTS.Pau.CG$Density.other <- 0

id.lvl <- unique(dat22.PTS.Pau.CG$ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat22.PTS.Pau.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Pau.CG$Density.other[pos.ind] <- mean22.PTS.Spite.CG$Density[pos.mean]
	}
}

dat22.Pau.CG$Density <- NA
id.lvl <- unique(dat22.Pau.CG$ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat22.Pau.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.anc.Pau.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.Pau.CG$Density[pos.ind] <- mean22.anc.Pau.CG$Density[pos.mean]
	}
}

dat22.PTS.Pau.CG$Density <- NA
id.lvl <- unique(dat22.PTS.Pau.CG$ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat22.PTS.Pau.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Pau.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Pau.CG$Density[pos.ind] <- mean22.PTS.Pau.CG$Density[pos.mean]
	}
}

dat05.Pau$Bio.Fraction.other <- 0
dat05.PTS.Pau$Bio.Fraction.other <- 0

id.lvl <- unique(dat05.PTS.Pau$Sample_ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat05.PTS.Pau$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Pau$Bio.Fraction.other[pos.ind] <- mean05.PTS.Spite$Bio.Fraction[pos.mean]
	}
}

dat05.Pau$Density.other <- 0
dat05.PTS.Pau$Density.other <- 0

id.lvl <- unique(dat05.PTS.Pau$Sample_ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat05.PTS.Pau$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Pau$Density.other[pos.ind] <- mean05.PTS.Spite$Density[pos.mean]
	}
}

dat05.Pau$Density <- NA
id.lvl <- unique(dat05.Pau$Sample_ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat05.Pau$Sample_ID == id.lvl[i])
	pos.mean <- which(mean.05.Pau$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.Pau$Density[pos.ind] <- mean.05.Pau$Density[pos.mean]
	}
}

dat05.PTS.Pau$Density <- NA
id.lvl <- unique(dat05.PTS.Pau$Sample_ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat05.PTS.Pau$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Pau$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Pau$Density[pos.ind] <- mean05.PTS.Pau$Density[pos.mean]
	}
}


## Spite ------

dat22.Spite.CG$Bio.Fraction.other <- 0
dat22.PTS.Spite.CG$Bio.Fraction.other <- 0

id.lvl <- unique(dat22.PTS.Spite.CG$ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat22.PTS.Spite.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Spite.CG$Bio.Fraction.other[pos.ind] <- mean22.PTS.Spite.CG$Bio.Fraction[pos.mean]
	}
}

dat22.Spite.CG$Density.other <- 0
dat22.PTS.Spite.CG$Density.other <- 0

id.lvl <- unique(dat22.PTS.Spite.CG$ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat22.PTS.Spite.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Spite.CG$Density.other[pos.ind] <- mean22.PTS.Spite.CG$Density[pos.mean]
	}
}

dat22.Spite.CG$Density <- NA
id.lvl <- unique(dat22.Spite.CG$ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat22.Spite.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.anc.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.Spite.CG$Density[pos.ind] <- mean22.anc.Spite.CG$Density[pos.mean]
	}
}

dat22.PTS.Spite.CG$Density <- NA
id.lvl <- unique(dat22.PTS.Spite.CG$ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat22.PTS.Spite.CG$ID == id.lvl[i])
	pos.mean <- which(mean22.PTS.Spite.CG$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat22.PTS.Spite.CG$Density[pos.ind] <- mean22.PTS.Spite.CG$Density[pos.mean]
	}
}

dat05.Spite$Bio.Fraction.other <- 0
dat05.PTS.Spite$Bio.Fraction.other <- 0

id.lvl <- unique(dat05.PTS.Spite$Sample_ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat05.PTS.Spite$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Spite$Bio.Fraction.other[pos.ind] <- mean05.PTS.Spite$Bio.Fraction[pos.mean]
	}
}

dat05.Spite$Density.other <- 0
dat05.PTS.Spite$Density.other <- 0

id.lvl <- unique(dat05.PTS.Spite$Sample_ID)
for(i in 1:length(id.lvl)){
	#print(i)
	pos.ind <- which(dat05.PTS.Spite$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Spite$Density.other[pos.ind] <- mean05.PTS.Spite$Density[pos.mean]
	}
}

dat05.Spite$Density <- NA
id.lvl <- unique(dat05.Spite$Sample_ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat05.Spite$Sample_ID == id.lvl[i])
	pos.mean <- which(mean.05.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.Spite$Density[pos.ind] <- mean.05.Spite$Density[pos.mean]
	}
}

dat05.PTS.Spite$Density <- NA
id.lvl <- unique(dat05.PTS.Spite$Sample_ID)
for(i in 1:length(id.lvl)){
	pos.ind <- which(dat05.PTS.Spite$Sample_ID == id.lvl[i])
	pos.mean <- which(mean05.PTS.Spite$ID == id.lvl[i])
	if(length(pos.mean)!=0){
		dat05.PTS.Spite$Density[pos.ind] <- mean05.PTS.Spite$Density[pos.mean]
	}
}


