## Read in the data -- change pathways
## GET THE CODE FOR THE PATHWAY SET AUTOMATICALLY

## Data at day 4 selection experiment
### Community
dat22.PTS.Pau.CG <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat22_PTS_Pau_CG.txt', sep=""))
dat22.PTS.Spite.CG <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat22_PTS_Spite_CG.txt', sep=""))
dat22.PTS.Tet.CG <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat22_PTS_Tet_CG.txt', sep=""), header=T)

### Isolation
dat22.Pau.CG <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat22_Pau_CG.txt', sep=""))
dat22.Spite.CG <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat22_Spite_CG.txt', sep=""))

## Data at day 78 selection experiment
### Community
dat05.PTS.Pau <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat05_PTS_Pau_CG.txt', sep=""))
dat05.PTS.Spite <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat05_PTS_Spite_CG.txt', sep=""))

### Isolation
dat05.Pau <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat05_Pau_CG.txt', sep=""))
dat05.Spite <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat05_Spite_CG.txt', sep=""))

## Data of common garden
### Community
dat09.PTS.Pau <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat09_PTS_Pau.txt', sep=""))
dat09.PTS.Spite <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat09_PTS_Spite.txt', sep=""))

### Isolation
dat09.Pau <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat09_Pau.txt', sep=""))
dat09.Spite <- read.table(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Data\\\\dat09_Spite.txt', sep=""))

## --------------------------------------------------------------------------
## Calculate densities and mean data sets

## Day 4 -- selection phase
cat22.anc.Pau.CG <- paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity)
mean22.anc.Pau.CG <- apply(dat22.Pau.CG[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat22.anc.Pau.CG, function(y) mean(y, na.rm=T)))
sd22.anc.Pau.CG <- apply(dat22.Pau.CG[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat22.anc.Pau.CG, function(y) sd(y, na.rm=T)))

names22.anc.Pau.CG <- RTD(mean22.anc.Pau.CG)
colnames(names22.anc.Pau.CG) <- c('ID', 'Salinity')
mean22.anc.Pau.CG <- cbind(names22.anc.Pau.CG, mean22.anc.Pau.CG)
sd22.anc.Pau.CG <- cbind(names22.anc.Pau.CG, sd22.anc.Pau.CG)

cat22.Pau.dens <- paste(dat22.Pau.CG$ID,dat22.Pau.CG$Salinity)
dens22.Pau <- tapply(dat22.Pau.CG$N_frames, cat22.Pau.dens, function(x) sum(x)/500/34.4*100)
speed22.Pau <- tapply(dat22.Pau.CG$gross_speed, cat22.Pau.dens, function(x) mean(x, na.rm=T))
mean22.anc.Pau.CG$Speed <- mean22.anc.Pau.CG$Density <- sd22.anc.Pau.CG$Density <- NA
for(i in 1:length(dens22.Pau)){
	pos.dens <- which(rownames(mean22.anc.Pau.CG) == names(dens22.Pau)[i])
	mean22.anc.Pau.CG$Density[pos.dens] <- dens22.Pau[[i]]
	sd22.anc.Pau.CG$Density[pos.dens] <- dens22.Pau[[i]]
	mean22.anc.Pau.CG$Speed[pos.dens] <- speed22.Pau[[i]]
}

cat22.anc.Spite.CG <- paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity)
mean22.anc.Spite.CG <- apply(dat22.Spite.CG[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat22.anc.Spite.CG, function(y) mean(y, na.rm=T)))
sd22.anc.Spite.CG <- apply(dat22.Spite.CG[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat22.anc.Spite.CG, function(y) sd(y, na.rm=T)))

names22.anc.Spite.CG <- RTD(mean22.anc.Spite.CG)
colnames(names22.anc.Spite.CG) <- c('ID', 'Salinity')
mean22.anc.Spite.CG <- cbind(names22.anc.Spite.CG, mean22.anc.Spite.CG)
sd22.anc.Spite.CG <- cbind(names22.anc.Spite.CG, sd22.anc.Spite.CG)

cat22.Spite.dens <- paste(dat22.Spite.CG$ID,dat22.Spite.CG$Salinity)
dens22.Spite <- tapply(dat22.Spite.CG$N_frames, cat22.Spite.dens, function(x) sum(x)/500/34.4*100)
speed22.Spite <- tapply(dat22.Spite.CG$gross_speed, cat22.Spite.dens, function(x) mean(x, na.rm=T))
mean22.anc.Spite.CG$Speed <- mean22.anc.Spite.CG$Density <- sd22.anc.Spite.CG$Density <- NA
for(i in 1:length(dens22.Spite)){
	pos.dens <- which(rownames(mean22.anc.Spite.CG) == names(dens22.Spite)[i])
	mean22.anc.Spite.CG$Density[pos.dens] <- dens22.Spite[[i]]
	sd22.anc.Spite.CG$Density[pos.dens] <- dens22.Spite[[i]]
	mean22.anc.Spite.CG$Speed[pos.dens] <- speed22.Spite[[i]]
}


cat22.PTS.Pau.CG <- paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity)
mean22.PTS.Pau.CG <- apply(dat22.PTS.Pau.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Pau.CG, function(y) mean(y, na.rm=T)))
sd22.PTS.Pau.CG <- apply(dat22.PTS.Pau.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Pau.CG, function(y) sd(y, na.rm=T)))

names22.PTS.Pau.CG <- RTD(mean22.PTS.Pau.CG)
colnames(names22.PTS.Pau.CG) <- c('ID', 'Salinity')
mean22.PTS.Pau.CG <- cbind(names22.PTS.Pau.CG, mean22.PTS.Pau.CG)
sd22.PTS.Pau.CG <- cbind(names22.PTS.Pau.CG, sd22.PTS.Pau.CG)

## Calculate the density (number of individuals for each ID)
cat22.PTS.Pau.dens <- paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity)
dens22.PTS.Pau <- tapply(dat22.PTS.Pau.CG$N_frames, cat22.PTS.Pau.dens, function(x) sum(x)/500/34.4*100)
speed22.PTS.Pau <- tapply(dat22.PTS.Pau.CG$gross_speed, cat22.PTS.Pau.dens, function(x) mean(x, na.rm=T))

mean22.PTS.Pau.CG$Speed <- mean22.PTS.Pau.CG$Density <- sd22.PTS.Pau.CG$Density <- NA
for(i in 1:length(dens22.PTS.Pau)){
	pos.dens <- which(rownames(mean22.PTS.Pau.CG) == names(dens22.PTS.Pau)[i])
	mean22.PTS.Pau.CG$Density[pos.dens] <- dens22.PTS.Pau[[i]]
	sd22.PTS.Pau.CG$Density[pos.dens] <- dens22.PTS.Pau[[i]]
	mean22.PTS.Pau.CG$Speed[pos.dens] <- speed22.PTS.Pau[[i]]
}

cat22.PTS.Spite.CG <- paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity)
mean22.PTS.Spite.CG <- apply(dat22.PTS.Spite.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Spite.CG, function(y) mean(y, na.rm=T)))
sd22.PTS.Spite.CG <- apply(dat22.PTS.Spite.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Spite.CG, function(y) sd(y, na.rm=T)))

names22.PTS.Spite.CG <- RTD(mean22.PTS.Spite.CG)
colnames(names22.PTS.Spite.CG) <- c('ID', 'Salinity')
mean22.PTS.Spite.CG <- cbind(names22.PTS.Spite.CG, mean22.PTS.Spite.CG)
sd22.PTS.Spite.CG <- cbind(names22.PTS.Spite.CG, sd22.PTS.Spite.CG)

## Calculate the density (number of individuals for each ID)
cat22.PTS.Spite.dens <- paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity)
dens22.PTS.Spite <- tapply(dat22.PTS.Spite.CG$N_frames, cat22.PTS.Spite.dens, function(x) sum(x)/500/34.4*100)
speed22.PTS.Spite <- tapply(dat22.PTS.Spite.CG$gross_speed, cat22.PTS.Spite.dens, function(x) mean(x, na.rm=T))

mean22.PTS.Spite.CG$Speed <- mean22.PTS.Spite.CG$Density <- sd22.PTS.Spite.CG$Density <- NA
for(i in 1:length(dens22.PTS.Spite)){
	pos.dens <- which(rownames(mean22.PTS.Spite.CG) == names(dens22.PTS.Spite)[i])
	mean22.PTS.Spite.CG$Density[pos.dens] <- dens22.PTS.Spite[[i]]
	sd22.PTS.Spite.CG$Density[pos.dens] <- dens22.PTS.Spite[[i]]
	mean22.PTS.Spite.CG$Speed[pos.dens] <- speed22.PTS.Spite[[i]]
}

cat22.PTS.Tet.CG <- paste(dat22.PTS.Tet.CG$ID, dat22.PTS.Tet.CG$Salinity)
mean22.PTS.Tet.CG <- apply(dat22.PTS.Tet.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Tet.CG, function(y) mean(y, na.rm=T)))
sd22.PTS.Tet.CG <- apply(dat22.PTS.Tet.CG[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat22.PTS.Tet.CG, function(y) sd(y, na.rm=T)))

names22.PTS.Tet.CG <- RTD(mean22.PTS.Tet.CG)
colnames(names22.PTS.Tet.CG) <- c('ID', 'Salinity')
mean22.PTS.Tet.CG <- cbind(names22.PTS.Tet.CG, mean22.PTS.Tet.CG)
sd22.PTS.Tet.CG <- cbind(names22.PTS.Tet.CG, sd22.PTS.Tet.CG)

## Calculate the density (number of individuals for each ID)
cat22.PTS.Tet.dens <- paste(dat22.PTS.Tet.CG$ID, dat22.PTS.Tet.CG$Salinity)
dens22.PTS.Tet <- tapply(dat22.PTS.Tet.CG$N_frames, cat22.PTS.Tet.dens, function(x) sum(x)/500/34.4*100)
speed22.PTS.Tet <- tapply(dat22.PTS.Tet.CG$gross_speed, cat22.PTS.Tet.dens, function(x) mean(x, na.rm=T))

mean22.PTS.Tet.CG$Speed <- mean22.PTS.Tet.CG$Density <- sd22.PTS.Tet.CG$Density <- NA
for(i in 1:length(dens22.PTS.Tet)){
	pos.dens <- which(rownames(mean22.PTS.Tet.CG) == names(dens22.PTS.Tet)[i])
	mean22.PTS.Tet.CG$Density[pos.dens] <- dens22.PTS.Tet[[i]]
	sd22.PTS.Tet.CG$Density[pos.dens] <- dens22.PTS.Tet[[i]]
	mean22.PTS.Tet.CG$Speed[pos.dens] <- speed22.PTS.Tet[[i]]
}

## ---------------------------------------------------------------------
## Day 78 -- selection phase
dat05.Pau$eccentricity <- dat05.Pau$mean_major/dat05.Pau$mean_minor

cat.05.Pau <- paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity)
mean.05.Pau <- apply(dat05.Pau[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat.05.Pau, function(y) mean(y, na.rm=T)))
sd.05.Pau <- apply(dat05.Pau[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat.05.Pau, function(y) sd(y, na.rm=T)))

names.05.Pau <- RTD(mean.05.Pau)
colnames(names.05.Pau) <- c('ID', 'Salinity')
mean.05.Pau <- cbind(names.05.Pau, mean.05.Pau)
sd.05.Pau <- cbind(names.05.Pau, sd.05.Pau)

cat05.Pau.dens <- paste(dat05.Pau$Sample_ID,dat05.Pau$Salinity)
dens05.Pau <- tapply(dat05.Pau$N_frames, cat05.Pau.dens, function(x) sum(x)/500/34.4*100)
speed05.Pau <- tapply(dat05.Pau$gross_speed, cat05.Pau.dens, function(x) mean(x, na.rm=T))
mean.05.Pau$Speed <- mean.05.Pau$Density <- sd.05.Pau$Density <- NA
for(i in 1:length(dens05.Pau)){
	pos.dens <- which(rownames(mean.05.Pau) == names(dens05.Pau)[i])
	mean.05.Pau$Density[pos.dens] <- dens05.Pau[[i]]
	sd.05.Pau$Density[pos.dens] <- dens05.Pau[[i]]
	mean.05.Pau$Speed[pos.dens] <- speed05.Pau[[i]]
}

dat05.Spite$eccentricity <- dat05.Spite$mean_major/dat05.Spite$mean_minor

cat.05.Spite <- paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity)
mean.05.Spite <- apply(dat05.Spite[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat.05.Spite, function(y) mean(y, na.rm=T)))
sd.05.Spite <- apply(dat05.Spite[,c(4,6,8,10,12,14,19:22,33)],2,function(x) tapply(x, cat.05.Spite, function(y) sd(y, na.rm=T)))

names.05.Spite <- RTD(mean.05.Spite)
colnames(names.05.Spite) <- c('ID', 'Salinity')
mean.05.Spite <- cbind(names.05.Spite, mean.05.Spite)
sd.05.Spite <- cbind(names.05.Spite, sd.05.Spite)

cat05.Spite.dens <- paste(dat05.Spite$Sample_ID,dat05.Spite$Salinity)
dens05.Spite <- tapply(dat05.Spite$N_frames, cat05.Spite.dens, function(x) sum(x)/500/34.4*100)
speed05.Spite <- tapply(dat05.Spite$gross_speed, cat05.Spite.dens, function(x) mean(x, na.rm=T))
mean.05.Spite$Speed <- mean.05.Spite$Density <- sd.05.Spite$Density <- NA
for(i in 1:length(dens05.Spite)){
	pos.dens <- which(rownames(mean.05.Spite) == names(dens05.Spite)[i])
	mean.05.Spite$Density[pos.dens] <- dens05.Spite[[i]]
	sd.05.Spite$Density[pos.dens] <- dens05.Spite[[i]]
	mean.05.Spite$Speed[pos.dens] <- speed05.Spite[[i]]
}

cat05.PTS.Pau <- paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity)
mean05.PTS.Pau <- apply(dat05.PTS.Pau[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat05.PTS.Pau, function(y) mean(y, na.rm=T)))
sd05.PTS.Pau <- apply(dat05.PTS.Pau[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat05.PTS.Pau, function(y) sd(y, na.rm=T)))

names05.PTS.Pau <- RTD(mean05.PTS.Pau)
colnames(names05.PTS.Pau) <- c('ID', 'Salinity')
mean05.PTS.Pau <- cbind(names05.PTS.Pau, mean05.PTS.Pau)
sd05.PTS.Pau <- cbind(names05.PTS.Pau, sd05.PTS.Pau)

## Calculate the density (number of individuals for each ID)
cat05.PTS.Pau.dens <- paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity)
dens05.PTS.Pau <- tapply(dat05.PTS.Pau$N_frames, cat05.PTS.Pau.dens, function(x) sum(x)/500/34.4*100)
speed05.PTS.Pau <- tapply(dat05.PTS.Pau$gross_speed, cat05.PTS.Pau.dens, function(x) mean(x, na.rm=T))

mean05.PTS.Pau$Speed <- mean05.PTS.Pau$Density <- sd05.PTS.Pau$Density <- NA
for(i in 1:length(dens05.PTS.Pau)){
	pos.dens <- which(rownames(mean05.PTS.Pau) == names(dens05.PTS.Pau)[i])
	mean05.PTS.Pau$Density[pos.dens] <- dens05.PTS.Pau[[i]]
	sd05.PTS.Pau$Density[pos.dens] <- dens05.PTS.Pau[[i]]
	mean05.PTS.Pau$Speed[pos.dens] <- speed05.PTS.Pau[[i]]
}

cat05.PTS.Spite <- paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity)
mean05.PTS.Spite <- apply(dat05.PTS.Spite[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat05.PTS.Spite, function(y) mean(y, na.rm=T)))
sd05.PTS.Spite <- apply(dat05.PTS.Spite[,c(4,6,8,10,12,14,19:22,34)],2,function(x) tapply(x, cat05.PTS.Spite, function(y) sd(y, na.rm=T)))

names05.PTS.Spite <- RTD(mean05.PTS.Spite)
colnames(names05.PTS.Spite) <- c('ID', 'Salinity')
mean05.PTS.Spite <- cbind(names05.PTS.Spite, mean05.PTS.Spite)
sd05.PTS.Spite <- cbind(names05.PTS.Spite, sd05.PTS.Spite)

## Calculate the density (number of individuals for each ID)
cat05.PTS.Spite.dens <- paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity)
dens05.PTS.Spite <- tapply(dat05.PTS.Spite$N_frames, cat05.PTS.Spite.dens, function(x) sum(x)/500/34.4*100)
speed05.PTS.Spite <- tapply(dat05.PTS.Spite$gross_speed, cat05.PTS.Spite.dens, function(x) mean(x, na.rm=T))

mean05.PTS.Spite$Speed <- mean05.PTS.Spite$Density <- sd05.PTS.Spite$Density <- NA
for(i in 1:length(dens05.PTS.Spite)){
	pos.dens <- which(rownames(mean05.PTS.Spite) == names(dens05.PTS.Spite)[i])
	mean05.PTS.Spite$Density[pos.dens] <- dens05.PTS.Spite[[i]]
	sd05.PTS.Spite$Density[pos.dens] <- dens05.PTS.Spite[[i]]
	mean05.PTS.Spite$Speed[pos.dens] <- speed05.PTS.Spite[[i]]
}

## ---------------------------------------------------------------------
## Day 82  -- common garden
dat09.Pau$eccentricity <- dat09.Pau$mean_major/dat09.Pau$mean_minor

## The replicate variable indicates which replicate it is - but mean can be taken across the samples (samples can be merged). 
### Calculating mean and sd
cat.09.Pau <- paste(dat09.Pau$ID_original, dat09.Pau$ID_new, dat09.Pau$Salinity_Origin, dat09.Pau$Salinity_Destination, dat09.Pau$replicate)
mean.09.Pau <- apply(dat09.Pau[,c(4,6,8,10,12,14,19:22,38)],2,function(x) tapply(x, cat.09.Pau, function(y) mean(y, na.rm=T))) #,38
sd.09.Pau <- apply(dat09.Pau[,c(4,6,8,10,12,14,19:22,38)],2,function(x) tapply(x, cat.09.Pau, function(y) sd(y, na.rm=T))) #,38

names.09.Pau <- RTD(mean.09.Pau)
colnames(names.09.Pau) <- c('ID.ori', 'ID.new','Sal.ori', 'Sal.des', 'replicate')
mean.09.Pau <- cbind(names.09.Pau, mean.09.Pau)
sd.09.Pau <- cbind(names.09.Pau, sd.09.Pau)

cat09.Pau.dens2 <- paste(dat09.Pau$ID_new, dat09.Pau$Sample)
dens09.Pau2 <- tapply(dat09.Pau$N_frames, cat09.Pau.dens2, function(x) sum(x)/500/34.4*100)
id.lvl <- unique(dat09.Pau$ID_new)
dens09.sample <- data.frame(CTD(dens09.Pau2))
colnames(dens09.sample) <- c('ID.new', 'Sample')
dens09.sample$Density <- dens09.Pau2

mean.dens.09.Pau <- tapply(dens09.sample$Density, dens09.sample$ID.new, mean)

cat09.Pau.dens <- paste(dat09.Pau$ID_original,dat09.Pau$ID_new,dat09.Pau$Salinity_Origin,dat09.Pau$Salinity_Destination,dat09.Pau$replicate)
dens09.Pau <- tapply(dat09.Pau$N_frames, cat09.Pau.dens, function(x) sum(x)/500/34.4*100)
speed09.Pau <- tapply(dat09.Pau$gross_speed, cat09.Pau.dens, function(x) mean(x, na.rm=T))
mean.09.Pau$Speed <- mean.09.Pau$Density <- sd.09.Pau$Density <- NA
for(i in 1:length(dens09.Pau)){
	pos.dens <- which(rownames(mean.09.Pau) == names(dens09.Pau)[i])
	mean.09.Pau$Density[pos.dens] <- dens09.Pau[[i]]
	sd.09.Pau$Density[pos.dens] <- dens09.Pau[[i]]
	mean.09.Pau$Speed[pos.dens] <- speed09.Pau[[i]]
}
for(i in 1:length(id.lvl)){
	mean.09.Pau$Density[mean.09.Pau$ID.new == id.lvl[i]] <- mean.dens.09.Pau[names(mean.dens.09.Pau) == id.lvl[i]][[1]]
}	

dat09.Spite$eccentricity <- dat09.Spite$mean_major/dat09.Spite$mean_minor

## The replicate variable indicates which replicate it is - but mean can be taken across the samples (samples can be merged). 
### Calculating mean and sd
cat.09.Spite <- paste(dat09.Spite$ID_original, dat09.Spite$ID_new, dat09.Spite$Salinity_Origin, dat09.Spite$Salinity_Destination, dat09.Spite$replicate)
mean.09.Spite <- apply(dat09.Spite[,c(4,6,8,10,12,14,19:22,38)],2,function(x) tapply(x, cat.09.Spite, function(y) mean(y, na.rm=T))) #,38
sd.09.Spite <- apply(dat09.Spite[,c(4,6,8,10,12,14,19:22,38)],2,function(x) tapply(x, cat.09.Spite, function(y) sd(y, na.rm=T))) #,38

names.09.Spite <- RTD(mean.09.Spite)
colnames(names.09.Spite) <- c('ID.ori', 'ID.new','Sal.ori', 'Sal.des', 'replicate')
mean.09.Spite <- cbind(names.09.Spite, mean.09.Spite)
sd.09.Spite <- cbind(names.09.Spite, sd.09.Spite)

cat09.Spite.dens2 <- paste(dat09.Spite$ID_new, dat09.Spite$Sample)
dens09.Spite2 <- tapply(dat09.Spite$N_frames, cat09.Spite.dens2, function(x) sum(x)/500/34.4*100)
id.lvl <- unique(dat09.Spite$ID_new)
dens09.sample <- data.frame(CTD(dens09.Spite2))
colnames(dens09.sample) <- c('ID.new', 'Sample')
dens09.sample$Density <- dens09.Spite2

mean.dens.09.Spite <- tapply(dens09.sample$Density, dens09.sample$ID.new, mean)

cat09.Spite.dens <- paste(dat09.Spite$ID_original,dat09.Spite$ID_new,dat09.Spite$Salinity_Origin,dat09.Spite$Salinity_Destination,dat09.Spite$replicate)
dens09.Spite <- tapply(dat09.Spite$N_frames, cat09.Spite.dens, function(x) sum(x)/500/34.4*100)
speed09.Spite <- tapply(dat09.Spite$gross_speed, cat09.Spite.dens, function(x) mean(x, na.rm=T))
mean.09.Spite$Speed <- mean.09.Spite$Density <- sd.09.Spite$Density <- NA
for(i in 1:length(dens09.Spite)){
	pos.dens <- which(rownames(mean.09.Spite) == names(dens09.Spite)[i])
	mean.09.Spite$Density[pos.dens] <- dens09.Spite[[i]]
	sd.09.Spite$Density[pos.dens] <- dens09.Spite[[i]]
	mean.09.Spite$Speed[pos.dens] <- speed09.Spite[[i]]
}
for(i in 1:length(id.lvl)){
	mean.09.Spite$Density[mean.09.Spite$ID.new == id.lvl[i]] <- mean.dens.09.Spite[names(mean.dens.09.Spite) == id.lvl[i]][[1]]
}	

### Calculating mean and sd
cat.PTS09.Pau <- paste(dat09.PTS.Pau$ID_original, dat09.PTS.Pau$ID_new, dat09.PTS.Pau$Salinity_Origin, dat09.PTS.Pau$Salinity_Destination, dat09.PTS.Pau$replicate)
mean.PTS09.Pau <- apply(dat09.PTS.Pau[,c(4,6,8,10,12,14,19:22,39)],2,function(x) tapply(x, cat.PTS09.Pau, function(y) mean(y, na.rm=T)))
sd.PTS09.Pau <- apply(dat09.PTS.Pau[,c(4,6,8,10,12,14,19:22,39)],2,function(x) tapply(x, cat.PTS09.Pau, function(y) sd(y, na.rm=T)))

names.PTS09.Pau <- RTD(mean.PTS09.Pau)
colnames(names.PTS09.Pau) <- c('ID.ori', 'ID.new','Sal.ori', 'Sal.des', 'replicate')
mean.PTS09.Pau <- cbind(names.PTS09.Pau, mean.PTS09.Pau)
sd.PTS09.Pau <- cbind(names.PTS09.Pau, sd.PTS09.Pau)

cat09.PTS.Pau.dens2 <- paste(dat09.PTS.Pau$ID_new, dat09.PTS.Pau$Sample)
dens09.PTS.Pau2 <- tapply(dat09.PTS.Pau$N_frames, cat09.PTS.Pau.dens2, function(x) sum(x)/500/34.4*100)
id.lvl <- unique(dat09.PTS.Pau$ID_new)
dens09.sample <- data.frame(CTD(dens09.PTS.Pau2))
colnames(dens09.sample) <- c('ID.new', 'Sample')
dens09.sample$Density <- dens09.PTS.Pau2

mean.dens.PTS09.Pau <- tapply(dens09.sample$Density, dens09.sample$ID.new, mean)

## Calculate the density (number of individuals for each ID and replicate)
cat09.PTS.Pau.dens <- paste(dat09.PTS.Pau$ID_original,dat09.PTS.Pau$ID_new,dat09.PTS.Pau$Salinity_Origin,dat09.PTS.Pau$Salinity_Destination,dat09.PTS.Pau$replicate)
dens09.PTS.Pau <- tapply(dat09.PTS.Pau$N_frames, cat09.PTS.Pau.dens, function(x) sum(x)/500/34.4*100)
speed09.PTS.Pau <- tapply(dat09.PTS.Pau$gross_speed, cat09.PTS.Pau.dens, function(x) mean(x, na.rm=T))
mean.PTS09.Pau$Speed <- mean.PTS09.Pau$Density <- sd.PTS09.Pau$Density <- NA
for(i in 1:length(dens09.PTS.Pau)){
	pos.dens <- which(rownames(mean.PTS09.Pau) == names(dens09.PTS.Pau)[i])
	mean.PTS09.Pau$Density[pos.dens] <- dens09.PTS.Pau[[i]]
	sd.PTS09.Pau$Density[pos.dens] <- dens09.PTS.Pau[[i]]
	mean.PTS09.Pau$Speed[pos.dens] <- speed09.PTS.Pau[[i]]
}
for(i in 1:length(id.lvl)){
	mean.PTS09.Pau$Density[mean.PTS09.Pau$ID.new == id.lvl[i]] <- mean.dens.PTS09.Pau[names(mean.dens.PTS09.Pau) == id.lvl[i]][[1]]
}

### Calculating mean and sd -- mixed cultures
cat.PTS09.Spite <- paste(dat09.PTS.Spite$ID_original, dat09.PTS.Spite$ID_new, dat09.PTS.Spite$Salinity_Origin, dat09.PTS.Spite$Salinity_Destination, dat09.PTS.Spite$replicate)
mean.PTS09.Spite <- apply(dat09.PTS.Spite[,c(4,6,8,10,12,14,19:22,39)],2,function(x) tapply(x, cat.PTS09.Spite, function(y) mean(y, na.rm=T)))
sd.PTS09.Spite <- apply(dat09.PTS.Spite[,c(4,6,8,10,12,14,19:22,39)],2,function(x) tapply(x, cat.PTS09.Spite, function(y) sd(y, na.rm=T)))

names.PTS09.Spite <- RTD(mean.PTS09.Spite)
colnames(names.PTS09.Spite) <- c('ID.ori', 'ID.new','Sal.ori', 'Sal.des', 'replicate')
mean.PTS09.Spite <- cbind(names.PTS09.Spite, mean.PTS09.Spite)
sd.PTS09.Spite <- cbind(names.PTS09.Spite, sd.PTS09.Spite)

## Calculate the density (number of individuals for each ID and replicate)
cat09.PTS.Spite.dens <- paste(dat09.PTS.Spite$ID_original,dat09.PTS.Spite$ID_new,dat09.PTS.Spite$Salinity_Origin,dat09.PTS.Spite$Salinity_Destination,dat09.PTS.Spite$replicate)
dens09.PTS.Spite <- tapply(dat09.PTS.Spite$N_frames, cat09.PTS.Spite.dens, function(x) sum(x)/500/34.4*100)
speed09.PTS.Spite <- tapply(dat09.PTS.Spite$gross_speed, cat09.PTS.Spite.dens, function(x) mean(x, na.rm=T))

cat09.PTS.Spite.dens2 <- paste(dat09.PTS.Spite$ID_new, dat09.PTS.Spite$Sample)
dens09.PTS.Spite2 <- tapply(dat09.PTS.Spite$N_frames, cat09.PTS.Spite.dens2, function(x) sum(x)/500/34.4*100)
id.lvl <- unique(dat09.PTS.Spite$ID_new)
dens09.sample <- data.frame(CTD(dens09.PTS.Spite2))
colnames(dens09.sample) <- c('ID.new', 'Sample')
dens09.sample$Density <- dens09.PTS.Spite2

mean.dens.PTS09.Spite <- tapply(dens09.sample$Density, dens09.sample$ID.new, mean)

## Calculate the density (number of individuals for each ID and replicate)
cat09.PTS.Spite.dens <- paste(dat09.PTS.Spite$ID_original,dat09.PTS.Spite$ID_new,dat09.PTS.Spite$Salinity_Origin,dat09.PTS.Spite$Salinity_Destination,dat09.PTS.Spite$replicate)
dens09.PTS.Spite <- tapply(dat09.PTS.Spite$N_frames, cat09.PTS.Spite.dens, function(x) sum(x)/500/34.4*100)
speed09.PTS.Spite <- tapply(dat09.PTS.Spite$gross_speed, cat09.PTS.Spite.dens, function(x) mean(x, na.rm=T))
mean.PTS09.Spite$Speed <- mean.PTS09.Spite$Density <- sd.PTS09.Spite$Density <- NA
for(i in 1:length(dens09.PTS.Spite)){
	pos.dens <- which(rownames(mean.PTS09.Spite) == names(dens09.PTS.Spite)[i])
	mean.PTS09.Spite$Density[pos.dens] <- dens09.PTS.Spite[[i]]
	sd.PTS09.Spite$Density[pos.dens] <- dens09.PTS.Spite[[i]]
	mean.PTS09.Spite$Speed[pos.dens] <- speed09.PTS.Spite[[i]]
}
for(i in 1:length(id.lvl)){
	mean.PTS09.Spite$Density[mean.PTS09.Spite$ID.new == id.lvl[i]] <- mean.dens.PTS09.Spite[names(mean.dens.PTS09.Spite) == id.lvl[i]][[1]]
}

## -----------------------------------------------------------------------------
## Adding some variables to the data needed for statistical analysis

dat09.Pau$Community <- 'mono'
dat09.PTS.Pau$Community <- 'mixed'

### Add the density variable
dat09.Pau$Density <- NA
id09 <- unique(mean.09.Pau$ID.new)
for(i in 1:length(id09)){
	pos.dat <- which(dat09.Pau$ID_new == id09[i])
	pos.mean <- which(mean.09.Pau$ID.new == id09[i])
	dat09.Pau$Density[pos.dat] <- mean.09.Pau$Density[pos.mean]
}
dat09.PTS.Pau$Density <- NA
idPTS09 <- unique(mean.PTS09.Pau$ID.new)
for(i in 1:length(idPTS09)){
	pos.dat <- which(dat09.PTS.Pau$ID_new == idPTS09[i])
	pos.mean <- which(mean.PTS09.Pau$ID.new == idPTS09[i])
	dat09.PTS.Pau$Density[pos.dat] <- mean.PTS09.Pau$Density[pos.mean]
}

## Regression on mixed and monocultures on the individual data
dat09.Spite$Community <- 'mono'
dat09.PTS.Spite$Community <- 'mixed'

### Add the density variable
dat09.Spite$Density <- NA
id09 <- unique(mean.09.Spite$ID.new)
for(i in 1:length(id09)){
	pos.dat <- which(dat09.Spite$ID_new == id09[i])
	pos.mean <- which(mean.09.Spite$ID.new == id09[i])
	dat09.Spite$Density[pos.dat] <- mean.09.Spite$Density[pos.mean]
}
dat09.PTS.Spite$Density <- NA
idPTS09 <- unique(mean.PTS09.Spite$ID.new)
for(i in 1:length(idPTS09)){
	pos.dat <- which(dat09.PTS.Spite$ID_new == idPTS09[i])
	pos.mean <- which(mean.PTS09.Spite$ID.new == idPTS09[i])
	dat09.PTS.Spite$Density[pos.dat] <- mean.PTS09.Spite$Density[pos.mean]
}

dat22.Pau.CG$Time <- 0
dat05.Pau$Time <- 1
dat22.PTS.Pau.CG$Time <- 0
dat05.PTS.Pau$Time <- 1

dat22.Pau.CG$Community <- 'mono'
dat05.Pau$Community <- 'mono'
dat22.PTS.Pau.CG$Community <- 'mixed'
dat05.PTS.Pau$Community <- 'mixed'

dat22.Spite.CG$Time <- 0
dat05.Spite$Time <- 1
dat22.PTS.Spite.CG$Time <- 0
dat05.PTS.Spite$Time <- 1

dat22.Spite.CG$Community <- 'mono'
dat05.Spite$Community <- 'mono'
dat22.PTS.Spite.CG$Community <- 'mixed'
dat05.PTS.Spite$Community <- 'mixed'




