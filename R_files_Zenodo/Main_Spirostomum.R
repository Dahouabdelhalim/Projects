## Overview

filestem <- 'C:\\\\Users\\\\Lynn\\\\'

## Libraries needed
library(lmerTest)
library(MuMIn)
library(car)

## Functions needed
source(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Function_RTD.R', sep=""))
source(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Function_RNplot.R', sep=""))

## Reading the data
source(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Reading_Data.R', sep=""))
source(paste(filestem, 'switchdrive\\\\ExperimentLuis_CommunityGarden_CommunityAssembly\\\\October2018\\\\Scripts_ForPublication\\\\Reading_AddingFractions.R', sep=""))


## ----------------------------------------------------------------------------------
## SPIROSTOMUM TERES -- SELECTION EXPERIMENT
## ----------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------
## Temporal trait shift -- selection -- Figure 2 Main text 
## ----------------------------------------------------------------------------------

par(mfrow = c(1,3))

m22.area <- tapply(dat22.Spite.CG$mean_area, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), mean)
sd22.area <- tapply(dat22.Spite.CG$mean_area, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), sd)
m22.area <- data.frame(cbind(RTD(m22.area), m22.area, sd22.area))
mPTS22.area <- tapply(dat22.PTS.Spite.CG$mean_area, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), mean)
sdPTS22.area <- tapply(dat22.PTS.Spite.CG$mean_area, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), sd)
mPTS22.area <- data.frame(cbind(RTD(mPTS22.area), mPTS22.area, sdPTS22.area))

m05.area <- tapply(dat05.Spite$mean_area, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), mean)
sd05.area <- tapply(dat05.Spite$mean_area, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), sd)
m05.area <- data.frame(cbind(RTD(m05.area), m05.area, sd05.area))
mPTS05.area <- tapply(dat05.PTS.Spite$mean_area, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), mean)
sdPTS05.area <- tapply(dat05.PTS.Spite$mean_area, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), sd)
mPTS05.area <- data.frame(cbind(RTD(mPTS05.area), mPTS05.area, sdPTS05.area))

colnames(m22.area) <- colnames(mPTS22.area) <- colnames(m05.area) <- colnames(mPTS05.area) <- c('ID', 'sal', 'area', 'sd.area')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(1700,10800), type = 'n', main = 'area')
for(i in 1:5){
	tmp.22 <- m22.area[m22.area$sal == sal.lvl[i],]
	tmp.PTS22 <- mPTS22.area[mPTS22.area$sal == sal.lvl[i],]
	tmp.05 <- m05.area[m05.area$sal == sal.lvl[i],]
	tmp.PTS05 <- mPTS05.area[mPTS05.area$sal == sal.lvl[i],]

	id.mo <- unique(c(tmp.22$ID, tmp.05$ID))
	id.mi <- unique(c(tmp.PTS22$ID, tmp.PTS05$ID))

	for(rr in 1:length(id.mo)){
		xx. <- tmp.22$area[tmp.22$ID == id.mo[rr]]
		yy. <- tmp.05$area[tmp.05$ID == id.mo[rr]]
		xx.sd <- tmp.22$sd.area[tmp.22$ID == id.mo[rr]]
		yy.sd <- tmp.05$sd.area[tmp.05$ID == id.mo[rr]]

		if(length(xx.) != 0){
			points(0.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(0.2+0.05*(i-1), xx.-xx.sd, 0.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(0.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(0.8+0.05*(i-1), yy.-yy.sd, 0.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(0.2+0.05*(i-1), xx., 0.8+0.05*(i-1), yy.,)}
	}

	for(rr in 1:length(id.mi)){
		xx. <- tmp.PTS22$area[tmp.PTS22$ID == id.mi[rr]]
		yy. <- tmp.PTS05$area[tmp.PTS05$ID == id.mi[rr]]
		xx.sd <- tmp.PTS22$sd.area[tmp.PTS22$ID == id.mi[rr]]
		yy.sd <- tmp.PTS05$sd.area[tmp.PTS05$ID == id.mi[rr]]

		if(length(xx.) != 0){
			points(1.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(1.2+0.05*(i-1), xx.-xx.sd, 1.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(1.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(1.8+0.05*(i-1), yy.-yy.sd, 1.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(1.2+0.05*(i-1), xx., 1.8+0.05*(i-1), yy.,)}
	}
}



m22.ar <- tapply(dat22.Spite.CG$mean_ar, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), mean)
sd22.ar <- tapply(dat22.Spite.CG$mean_ar, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), sd)
m22.ar <- data.frame(cbind(RTD(m22.ar), m22.ar, sd22.ar))
mPTS22.ar <- tapply(dat22.PTS.Spite.CG$mean_ar, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), mean)
sdPTS22.ar <- tapply(dat22.PTS.Spite.CG$mean_ar, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), sd)
mPTS22.ar <- data.frame(cbind(RTD(mPTS22.ar), mPTS22.ar, sdPTS22.ar))

m05.ar <- tapply(dat05.Spite$mean_ar, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), mean)
sd05.ar <- tapply(dat05.Spite$mean_ar, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), sd)
m05.ar <- data.frame(cbind(RTD(m05.ar), m05.ar, sd05.ar))
mPTS05.ar <- tapply(dat05.PTS.Spite$mean_ar, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), mean)
sdPTS05.ar <- tapply(dat05.PTS.Spite$mean_ar, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), sd)
mPTS05.ar <- data.frame(cbind(RTD(mPTS05.ar), mPTS05.ar, sdPTS05.ar))

colnames(m22.ar) <- colnames(mPTS22.ar) <- colnames(m05.ar) <- colnames(mPTS05.ar) <- c('ID', 'sal', 'ar', 'sd.ar')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(2.7,11), type = 'n', main = 'ar')
for(i in 1:5){
	tmp.22 <- m22.ar[m22.ar$sal == sal.lvl[i],]
	tmp.PTS22 <- mPTS22.ar[mPTS22.ar$sal == sal.lvl[i],]
	tmp.05 <- m05.ar[m05.ar$sal == sal.lvl[i],]
	tmp.PTS05 <- mPTS05.ar[mPTS05.ar$sal == sal.lvl[i],]

	id.mo <- unique(c(tmp.22$ID, tmp.05$ID))
	id.mi <- unique(c(tmp.PTS22$ID, tmp.PTS05$ID))

	for(rr in 1:length(id.mo)){
		xx. <- tmp.22$ar[tmp.22$ID == id.mo[rr]]
		yy. <- tmp.05$ar[tmp.05$ID == id.mo[rr]]
		xx.sd <- tmp.22$sd.ar[tmp.22$ID == id.mo[rr]]
		yy.sd <- tmp.05$sd.ar[tmp.05$ID == id.mo[rr]]

		if(length(xx.) != 0){
			points(0.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(0.2+0.05*(i-1), xx.-xx.sd, 0.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(0.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(0.8+0.05*(i-1), yy.-yy.sd, 0.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(0.2+0.05*(i-1), xx., 0.8+0.05*(i-1), yy.,)}
	}

	for(rr in 1:length(id.mi)){
		xx. <- tmp.PTS22$ar[tmp.PTS22$ID == id.mi[rr]]
		yy. <- tmp.PTS05$ar[tmp.PTS05$ID == id.mi[rr]]
		xx.sd <- tmp.PTS22$sd.ar[tmp.PTS22$ID == id.mi[rr]]
		yy.sd <- tmp.PTS05$sd.ar[tmp.PTS05$ID == id.mi[rr]]

		if(length(xx.) != 0){
			points(1.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(1.2+0.05*(i-1), xx.-xx.sd, 1.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(1.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(1.8+0.05*(i-1), yy.-yy.sd, 1.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(1.2+0.05*(i-1), xx., 1.8+0.05*(i-1), yy.,)}
	}
}


m22.speed <- tapply(dat22.Spite.CG$gross_speed, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), mean)
sd22.speed <- tapply(dat22.Spite.CG$gross_speed, paste(dat22.Spite.CG$ID, dat22.Spite.CG$Salinity), sd)
m22.speed <- data.frame(cbind(RTD(m22.speed), m22.speed, sd22.speed))
mPTS22.speed <- tapply(dat22.PTS.Spite.CG$gross_speed, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), mean)
sdPTS22.speed <- tapply(dat22.PTS.Spite.CG$gross_speed, paste(dat22.PTS.Spite.CG$ID, dat22.PTS.Spite.CG$Salinity), sd)
mPTS22.speed <- data.frame(cbind(RTD(mPTS22.speed), mPTS22.speed, sdPTS22.speed))

m05.speed <- tapply(dat05.Spite$gross_speed, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), mean)
sd05.speed <- tapply(dat05.Spite$gross_speed, paste(dat05.Spite$Sample_ID, dat05.Spite$Salinity), sd)
m05.speed <- data.frame(cbind(RTD(m05.speed), m05.speed, sd05.speed))
mPTS05.speed <- tapply(dat05.PTS.Spite$gross_speed, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), mean)
sdPTS05.speed <- tapply(dat05.PTS.Spite$gross_speed, paste(dat05.PTS.Spite$Sample_ID, dat05.PTS.Spite$Salinity), sd)
mPTS05.speed <- data.frame(cbind(RTD(mPTS05.speed), mPTS05.speed, sdPTS05.speed))

colnames(m22.speed) <- colnames(mPTS22.speed) <- colnames(m05.speed) <- colnames(mPTS05.speed) <- c('ID', 'sal', 'speed', 'sd.speed')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(100,1000), type = 'n', main = 'speed')
for(i in 1:5){
	tmp.22 <- m22.speed[m22.speed$sal == sal.lvl[i],]
	tmp.PTS22 <- mPTS22.speed[mPTS22.speed$sal == sal.lvl[i],]
	tmp.05 <- m05.speed[m05.speed$sal == sal.lvl[i],]
	tmp.PTS05 <- mPTS05.speed[mPTS05.speed$sal == sal.lvl[i],]

	id.mo <- unique(c(tmp.22$ID, tmp.05$ID))
	id.mi <- unique(c(tmp.PTS22$ID, tmp.PTS05$ID))

	for(rr in 1:length(id.mo)){
		xx. <- tmp.22$speed[tmp.22$ID == id.mo[rr]]
		yy. <- tmp.05$speed[tmp.05$ID == id.mo[rr]]
		xx.sd <- tmp.22$sd.speed[tmp.22$ID == id.mo[rr]]
		yy.sd <- tmp.05$sd.speed[tmp.05$ID == id.mo[rr]]

		if(length(xx.) != 0){
			points(0.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(0.2+0.05*(i-1), xx.-xx.sd, 0.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(0.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(0.8+0.05*(i-1), yy.-yy.sd, 0.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(0.2+0.05*(i-1), xx., 0.8+0.05*(i-1), yy.,)}
	}

	for(rr in 1:length(id.mi)){
		xx. <- tmp.PTS22$speed[tmp.PTS22$ID == id.mi[rr]]
		yy. <- tmp.PTS05$speed[tmp.PTS05$ID == id.mi[rr]]
		xx.sd <- tmp.PTS22$sd.speed[tmp.PTS22$ID == id.mi[rr]]
		yy.sd <- tmp.PTS05$sd.speed[tmp.PTS05$ID == id.mi[rr]]

		if(length(xx.) != 0){
			points(1.2+0.05*(i-1), xx., pch = 21, cex = 2)
			segments(1.2+0.05*(i-1), xx.-xx.sd, 1.2+0.05*(i-1), xx.+xx.sd)}
		if(length(yy.) != 0){
			points(1.8+0.05*(i-1), yy., pch = 21, cex = 2)
			segments(1.8+0.05*(i-1), yy.-yy.sd, 1.8+0.05*(i-1), yy.+yy.sd)}
		if(length(xx.) != 0 & length(yy.) != 0){
			segments(1.2+0.05*(i-1), xx., 1.8+0.05*(i-1), yy.,)}
	}
}

### -------------------------------------------------------------------------------------------------
## Temporal trait shift to salinity -- selection -- Supplementary Table S2
### -------------------------------------------------------------------------------------------------

tt.sal <- c(dat22.Spite.CG$Salinity, dat05.Spite$Salinity, dat22.PTS.Spite.CG$Salinity, dat05.PTS.Spite$Salinity)
tt.time <- c(dat22.Spite.CG$Time, dat05.Spite$Time, dat22.PTS.Spite.CG$Time, dat05.PTS.Spite$Time)
tt.id <- c(dat22.Spite.CG$ID, dat05.Spite$Sample_ID, dat22.PTS.Spite.CG$ID, dat05.PTS.Spite$Sample_ID)
tt.com <- factor(c(dat22.Spite.CG$Community, dat05.Spite$Community, dat22.PTS.Spite.CG$Community, dat05.PTS.Spite$Community), levels=c('mono','mixed'))
tt.bf <- c(dat22.Spite.CG$Bio.Fraction.other, dat05.Spite$Bio.Fraction.other, dat22.PTS.Spite.CG$Bio.Fraction.other, dat05.PTS.Spite$Bio.Fraction.other)
tt.dens <- c(dat22.Spite.CG$Density, dat05.Spite$Density, dat22.PTS.Spite.CG$Density, dat05.PTS.Spite$Density)
tt.area <- c(dat22.Spite.CG$mean_area, dat05.Spite$mean_area, dat22.PTS.Spite.CG$mean_area, dat05.PTS.Spite$mean_area)
tt.ar <- c(dat22.Spite.CG$mean_ar, dat05.Spite$mean_ar, dat22.PTS.Spite.CG$mean_ar, dat05.PTS.Spite$mean_ar)
tt.speed <- c(dat22.Spite.CG$gross_speed, dat05.Spite$gross_speed, dat22.PTS.Spite.CG$gross_speed, dat05.PTS.Spite$gross_speed)

fit.area <- lmer(tt.area ~ tt.sal*tt.com*tt.time + tt.bf + tt.dens + (1|tt.id))
summary(fit.area)
r.squaredGLMM(fit.area)

plot(fit.area)
plot(cooks.distance(fit.area))
qqnorm(resid(fit.area))
qqline(resid(fit.area))
hist(resid(fit.area), breaks=50)

fit.ar <- lmer(tt.ar ~ tt.sal*tt.com*tt.time + tt.bf + tt.dens + (1|tt.id))
summary(fit.ar)
r.squaredGLMM(fit.ar)

plot(fit.ar)
plot(cooks.distance(fit.ar))
qqnorm(resid(fit.ar))
qqline(resid(fit.ar))
hist(resid(fit.ar), breaks=50)

fit.speed <- lmer(tt.speed ~ tt.sal*tt.com*tt.time + tt.bf + tt.dens + (1|tt.id))
summary(fit.speed)
r.squaredGLMM(fit.speed)

plot(fit.speed)
plot(cooks.distance(fit.speed))
qqnorm(resid(fit.speed))
qqline(resid(fit.speed))
hist(resid(fit.speed), breaks=50)

### -------------------------------------------------------------------------------------------------
## Phenotypic trait response to salinity and competition -- common garden -- Figure 3 main text
### -------------------------------------------------------------------------------------------------

dat09.combined.Spite <- rbind(dat09.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Spite$Salinity_Origin==0 & dat09.combined.Spite$Salinity_Destination==0)
pos.X42 <- which(dat09.combined.Spite$Salinity_Origin==4 | dat09.combined.Spite$Salinity_Origin==2)
w0.dat09.combined.Spite <- dat09.combined.Spite[-c(pos.0, pos.X42),]
w0.dat09.combined.Spite$Community <- factor(w0.dat09.combined.Spite$Community, levels=c('mono','mixed'))
dim(w0.dat09.combined.Spite)

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)

par(mfrow=c(1,3))
c.area <- summary(fit.area)$coefficients

plot(0:5, rep(-10000,6), type='n', ylim=c(-3500,2000))
abline(h=0)
rect(0,0,0.5,c.area[2,1])
segments(0.25,c.area[2,1]-c.area[2,2],0.25,c.area[2,1]+c.area[2,2])
rect(1,0,1.5,c.area[3,1])
segments(1.25,c.area[3,1]-c.area[3,2],1.25,c.area[3,1]+c.area[3,2])

rect(2,0,2.5,c.area[4,1])
segments(2.25,c.area[4,1]-c.area[4,2],2.25,c.area[4,1]+c.area[4,2])
rect(3,0,3.5,c.area[8,1])
segments(3.25,c.area[8,1]-c.area[8,2],3.25,c.area[8,1]+c.area[8,2])
rect(4,0,4.5,c.area[9,1])
segments(4.25,c.area[9,1]-c.area[9,2],4.25,c.area[9,1]+c.area[9,2])

c.ar <- summary(fit.ar)$coefficients

plot(0:5, rep(-10000,6), type='n', ylim=c(-4,2))
abline(h=0)
rect(0,0,0.5,c.ar[2,1])
segments(0.25,c.ar[2,1]-c.ar[2,2],0.25,c.ar[2,1]+c.ar[2,2])
rect(1,0,1.5,c.ar[3,1])
segments(1.25,c.ar[3,1]-c.ar[3,2],1.25,c.ar[3,1]+c.ar[3,2])

rect(2,0,2.5,c.ar[4,1])
segments(2.25,c.ar[4,1]-c.ar[4,2],2.25,c.ar[4,1]+c.ar[4,2])
rect(3,0,3.5,c.ar[8,1])
segments(3.25,c.ar[8,1]-c.ar[8,2],3.25,c.ar[8,1]+c.ar[8,2])
rect(4,0,4.5,c.ar[9,1])
segments(4.25,c.ar[9,1]-c.ar[9,2],4.25,c.ar[9,1]+c.ar[9,2])


c.speed <- summary(fit.speed)$coefficients

plot(0:5, rep(-10000,6), type='n', ylim=c(-200,70))
abline(h=0)
rect(0,0,0.5,c.speed[2,1])
segments(0.25,c.speed[2,1]-c.speed[2,2],0.25,c.speed[2,1]+c.speed[2,2])
rect(1,0,1.5,c.speed[3,1])
segments(1.25,c.speed[3,1]-c.speed[3,2],1.25,c.speed[3,1]+c.speed[3,2])

rect(2,0,2.5,c.speed[4,1])
segments(2.25,c.speed[4,1]-c.speed[4,2],2.25,c.speed[4,1]+c.speed[4,2])
rect(3,0,3.5,c.speed[8,1])
segments(3.25,c.speed[8,1]-c.speed[8,2],3.25,c.speed[8,1]+c.speed[8,2])
rect(4,0,4.5,c.speed[9,1])
segments(4.25,c.speed[9,1]-c.speed[9,2],4.25,c.speed[9,1]+c.speed[9,2])


### -------------------------------------------------------------------------------------------------
## Phenotypic trait response to salinity and competition -- common garden -- Supplementary Table S4
### -------------------------------------------------------------------------------------------------

dat09.combined.Spite <- rbind(dat09.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Spite$Salinity_Origin==0 & dat09.combined.Spite$Salinity_Destination==0)
pos.X42 <- which(dat09.combined.Spite$Salinity_Origin==4 | dat09.combined.Spite$Salinity_Origin==2)
w0.dat09.combined.Spite <- dat09.combined.Spite[-c(pos.0, pos.X42),]
w0.dat09.combined.Spite$Community <- factor(w0.dat09.combined.Spite$Community, levels=c('mono','mixed'))
dim(w0.dat09.combined.Spite)

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.area)
r.squaredGLMM(fit.area)

plot(fit.area)
plot(cooks.distance(fit.area))
qqnorm(resid(fit.area))
qqline(resid(fit.area))
hist(resid(fit.area), breaks=50)

fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.ar)
r.squaredGLMM(fit.ar)

plot(fit.ar)
plot(cooks.distance(fit.ar))
qqnorm(resid(fit.ar))
qqline(resid(fit.ar))
hist(resid(fit.ar), breaks=50)

fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.speed)
r.squaredGLMM(fit.speed)

plot(fit.speed)
plot(cooks.distance(fit.speed))
qqnorm(resid(fit.speed))
qqline(resid(fit.speed))
hist(resid(fit.speed), breaks=50)

### -------------------------------------------------------------------------------------------------
## Reaction norms of common garden results -- Figure S2
### -------------------------------------------------------------------------------------------------

## area
## Preparatory calculations
cat. <- paste(dat09.combined.Spite$ID_original, dat09.combined.Spite$Salinity_Origin, dat09.combined.Spite$Salinity_Destination, dat09.combined.Spite$replicate,
		dat09.combined.Spite$Community)
des.area <- tapply(dat09.combined.Spite$mean_area, cat., mean)
des.area.sd <- tapply(dat09.combined.Spite$mean_area, cat., sd)
des.area <- data.frame(cbind(RTD(des.area), des.area, des.area.sd))
colnames(des.area) <- c('ID_ori', 'sal.ori', 'sal.des', 'replicate', 'com', 'm.area', 'sd.area')

des.mo <- des.area[des.area$com == 'mono',]
des.mi <- des.area[des.area$com == 'mixed',]

mm.des.mo.id <- tapply(des.mo$m.area, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), mean)
sd.des.mo.id <- tapply(des.mo$m.area, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), sd)
mm.des.mo.id <- data.frame(cbind(RTD(mm.des.mo.id), mm.des.mo.id, sd.des.mo.id))
colnames(mm.des.mo.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.area', 'sd.area')

mm.des.mi.id <- tapply(des.mi$m.area, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), mean)
sd.des.mi.id <- tapply(des.mi$m.area, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), sd)
mm.des.mi.id <- data.frame(cbind(RTD(mm.des.mi.id), mm.des.mi.id, sd.des.mi.id))
colnames(mm.des.mi.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.area', 'sd.area')

mm.des.mo <- tapply(des.mo$m.area, paste(des.mo$sal.ori, des.mo$sal.des), mean)
sd.des.mo <- tapply(des.mo$m.area, paste(des.mo$sal.ori, des.mo$sal.des), sd)
mm.des.mo <- data.frame(cbind(RTD(mm.des.mo), mm.des.mo, sd.des.mo))
colnames(mm.des.mo) <- c('sal.ori', 'sal.des', 'm.area', 'sd.area')

mm.des.mi <- tapply(des.mi$m.area, paste(des.mi$sal.ori, des.mi$sal.des), mean)
sd.des.mi <- tapply(des.mi$m.area, paste(des.mi$sal.ori, des.mi$sal.des), sd)
mm.des.mi <- data.frame(cbind(RTD(mm.des.mi), mm.des.mi, sd.des.mi))
colnames(mm.des.mi) <- c('sal.ori', 'sal.des', 'm.area', 'sd.area')

## Plot monocultures
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(1500, 11500), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - area') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.area[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.area[id.tmp$sal.des == 0]-id.tmp$sd.area[id.tmp$sal.des == 0],-0.2, id.tmp$m.area[id.tmp$sal.des == 0]+id.tmp$sd.area[id.tmp$sal.des == 0])
		}
		x05 <- id.tmp$m.area[id.tmp$sal.des == 0.5]
		if(length(x05) > 0){
			points(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]-id.tmp$sd.area[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]+id.tmp$sd.area[id.tmp$sal.des == 0.5])
		}
		x1 <- id.tmp$m.area[id.tmp$sal.des == 1]
		if(length(x1) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1], pch = 21, cex = 1)
			segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]-id.tmp$sd.area[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]+id.tmp$sd.area[id.tmp$sal.des == 1])
		}
		x2 <- id.tmp$m.area[id.tmp$sal.des == 2]
		if(length(x2) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2], pch = 22, cex = 1)
			segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]-id.tmp$sd.area[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]+id.tmp$sd.area[id.tmp$sal.des == 2])
		}	
		if(length(x05) > 0){
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1])
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.area[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]-id.tmp$sd.area[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]+id.tmp$sd.area[id.tmp$sal.des == 4])
			if(length(x05) > 0){			
				segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mo[mm.des.mo$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.area[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.area[tmp.$sal.des == 0]-tmp.$sd.area[tmp.$sal.des == 0],-0.2, tmp.$m.area[tmp.$sal.des == 0]+tmp.$sd.area[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]-tmp.$sd.area[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]-tmp.$sd.area[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]-tmp.$sd.area[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]-tmp.$sd.area[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]-tmp.$sd.area[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]-tmp.$sd.area[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.area[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]-tmp.$sd.area[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]-tmp.$sd.area[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], lwd = 2)
	}
}

## Plot competition
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(1500, 11500), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - area') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.area[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.area[id.tmp$sal.des == 0]-id.tmp$sd.area[id.tmp$sal.des == 0],-0.2, id.tmp$m.area[id.tmp$sal.des == 0]+id.tmp$sd.area[id.tmp$sal.des == 0])
		}
		x05 <- id.tmp$m.area[id.tmp$sal.des == 0.5]
		if(length(x05) > 0){
			points(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]-id.tmp$sd.area[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]+id.tmp$sd.area[id.tmp$sal.des == 0.5])
		}
		x1 <- id.tmp$m.area[id.tmp$sal.des == 1]
		if(length(x1) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1], pch = 21, cex = 1)
			segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]-id.tmp$sd.area[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]+id.tmp$sd.area[id.tmp$sal.des == 1])
		}
		x2 <- id.tmp$m.area[id.tmp$sal.des == 2]
		if(length(x2) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2], pch = 22, cex = 1)
			segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]-id.tmp$sd.area[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]+id.tmp$sd.area[id.tmp$sal.des == 2])
		}	
		if(length(x05) > 0 & length(x1) > 0){
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1])
		}
		if(length(x05) > 0 & length(x2) > 0){
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.area[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]-id.tmp$sd.area[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]+id.tmp$sd.area[id.tmp$sal.des == 4])
			if(length(x05) > 0){			
				segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.area[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.area[tmp.$sal.des == 0]-tmp.$sd.area[tmp.$sal.des == 0],-0.2, tmp.$m.area[tmp.$sal.des == 0]+tmp.$sd.area[tmp.$sal.des == 0])
	}
	x05 <- tmp.$m.area[tmp.$sal.des == 0.5]
	if(length(x05) > 0){
		points(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5], pch = 21, cex = 2)
		segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]-tmp.$sd.area[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]-tmp.$sd.area[tmp.$sal.des == 0.5])
	}
	x1 <- tmp.$m.area[tmp.$sal.des == 1]
	if(length(x1) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], pch = 21, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]-tmp.$sd.area[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]-tmp.$sd.area[tmp.$sal.des == 1])
	}
	x2 <- tmp.$m.area[tmp.$sal.des == 2]
	if(length(x2) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], pch = 22, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]-tmp.$sd.area[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]-tmp.$sd.area[tmp.$sal.des == 2])
	}
	if(length(x05) > 0 & length(x1) > 0){
		segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], lwd = 2)
	}
	if(length(x05) > 0 & length(x2) > 0){
		segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], lwd = 2)
	}

	x <- tmp.$m.area[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]-tmp.$sd.area[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]-tmp.$sd.area[tmp.$sal.des == 4])
		if(length(x05) > 0){		
			segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], lwd = 2)
		}
	}
}


## Cell shape 
cat. <- paste(dat09.combined.Spite$ID_original, dat09.combined.Spite$Salinity_Origin, dat09.combined.Spite$Salinity_Destination, dat09.combined.Spite$replicate,
		dat09.combined.Spite$Community)
des.ar <- tapply(dat09.combined.Spite$mean_ar, cat., mean)
des.ar.sd <- tapply(dat09.combined.Spite$mean_ar, cat., sd)
des.ar <- data.frame(cbind(RTD(des.ar), des.ar, des.ar.sd))
colnames(des.ar) <- c('ID_ori', 'sal.ori', 'sal.des', 'replicate', 'com', 'm.ar', 'sd.ar')

des.mo <- des.ar[des.ar$com == 'mono',]
des.mi <- des.ar[des.ar$com == 'mixed',]

mm.des.mo.id <- tapply(des.mo$m.ar, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), mean)
sd.des.mo.id <- tapply(des.mo$m.ar, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), sd)
mm.des.mo.id <- data.frame(cbind(RTD(mm.des.mo.id), mm.des.mo.id, sd.des.mo.id))
colnames(mm.des.mo.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.ar', 'sd.ar')

mm.des.mi.id <- tapply(des.mi$m.ar, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), mean)
sd.des.mi.id <- tapply(des.mi$m.ar, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), sd)
mm.des.mi.id <- data.frame(cbind(RTD(mm.des.mi.id), mm.des.mi.id, sd.des.mi.id))
colnames(mm.des.mi.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.ar', 'sd.ar')

mm.des.mo <- tapply(des.mo$m.ar, paste(des.mo$sal.ori, des.mo$sal.des), mean)
sd.des.mo <- tapply(des.mo$m.ar, paste(des.mo$sal.ori, des.mo$sal.des), sd)
mm.des.mo <- data.frame(cbind(RTD(mm.des.mo), mm.des.mo, sd.des.mo))
colnames(mm.des.mo) <- c('sal.ori', 'sal.des', 'm.ar', 'sd.ar')

mm.des.mi <- tapply(des.mi$m.ar, paste(des.mi$sal.ori, des.mi$sal.des), mean)
sd.des.mi <- tapply(des.mi$m.ar, paste(des.mi$sal.ori, des.mi$sal.des), sd)
mm.des.mi <- data.frame(cbind(RTD(mm.des.mi), mm.des.mi, sd.des.mi))
colnames(mm.des.mi) <- c('sal.ori', 'sal.des', 'm.ar', 'sd.ar')

## Linear reaction norm figure
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(2, 7.7), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - ar') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 0])){
			segments(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]-id.tmp$sd.ar[id.tmp$sal.des == 0],-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]+id.tmp$sd.ar[id.tmp$sal.des == 0])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			points(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 0.5])){
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]-id.tmp$sd.ar[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]+id.tmp$sd.ar[id.tmp$sal.des == 0.5])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 1]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 1])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]-id.tmp$sd.ar[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]+id.tmp$sd.ar[id.tmp$sal.des == 1])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 2]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2], pch = 22, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 2])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]-id.tmp$sd.ar[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]+id.tmp$sd.ar[id.tmp$sal.des == 2])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1])
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.ar[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4], pch = 24, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 4])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4]-id.tmp$sd.ar[id.tmp$sal.des == 4], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4]+id.tmp$sd.ar[id.tmp$sal.des == 4])
			}
			x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
			if(length(x) > 0){	
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mo[mm.des.mo$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.ar[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.ar[tmp.$sal.des == 0]-tmp.$sd.ar[tmp.$sal.des == 0],-0.2, tmp.$m.ar[tmp.$sal.des == 0]+tmp.$sd.ar[tmp.$sal.des == 0], col= 'red')
	}
	points(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]+tmp.$sd.ar[tmp.$sal.des == 0.5], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]+tmp.$sd.ar[tmp.$sal.des == 1], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]+tmp.$sd.ar[tmp.$sal.des == 2], col = 'red')
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.ar[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]+tmp.$sd.ar[tmp.$sal.des == 4], col = 'red')
		segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], lwd = 2)
	}
}

## MAKE THE ONE FOR THE COMMUNITY
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(2, 7.7), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - ar') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 0])){
			segments(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]-id.tmp$sd.ar[id.tmp$sal.des == 0],-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]+id.tmp$sd.ar[id.tmp$sal.des == 0])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			points(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 0.5])){
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]-id.tmp$sd.ar[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]+id.tmp$sd.ar[id.tmp$sal.des == 0.5])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 1]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 1])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]-id.tmp$sd.ar[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]+id.tmp$sd.ar[id.tmp$sal.des == 1])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 2]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2], pch = 22, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 2])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]-id.tmp$sd.ar[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]+id.tmp$sd.ar[id.tmp$sal.des == 2])
			}
		}
		x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
		y <- id.tmp$m.ar[id.tmp$sal.des == 1]
		z <- id.tmp$m.ar[id.tmp$sal.des == 2]
		if(length(x) > 0 & length(y) > 0){	
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1])
		}
		if(length(x) > 0 & length(z) > 0){	
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.ar[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4], pch = 24, cex = 1)
			if(!is.na(id.tmp$sd.ar[id.tmp$sal.des == 4])){
			segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4]-id.tmp$sd.ar[id.tmp$sal.des == 4], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4]+id.tmp$sd.ar[id.tmp$sal.des == 4])
			}
			x <- id.tmp$m.ar[id.tmp$sal.des == 0.5]
			y <- id.tmp$m.ar[id.tmp$sal.des == 4]
			if(length(x) > 0 & length(y) > 0){	
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:3){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.ar[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.ar[tmp.$sal.des == 0]-tmp.$sd.ar[tmp.$sal.des == 0],-0.2, tmp.$m.ar[tmp.$sal.des == 0]+tmp.$sd.ar[tmp.$sal.des == 0], col= 'red')
	}
	points(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]+tmp.$sd.ar[tmp.$sal.des == 0.5], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]+tmp.$sd.ar[tmp.$sal.des == 1], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]+tmp.$sd.ar[tmp.$sal.des == 2], col = 'red')
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.ar[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]+tmp.$sd.ar[tmp.$sal.des == 4], col = 'red')
		segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], lwd = 2)
	}
}

## speed
cat. <- paste(dat09.combined.Spite$ID_original, dat09.combined.Spite$Salinity_Origin, dat09.combined.Spite$Salinity_Destination, dat09.combined.Spite$replicate,
		dat09.combined.Spite$Community)
des.speed <- tapply(dat09.combined.Spite$gross_speed, cat., mean)
des.speed.sd <- tapply(dat09.combined.Spite$gross_speed, cat., sd)
des.speed <- data.frame(cbind(RTD(des.speed), des.speed, des.speed.sd))
colnames(des.speed) <- c('ID_ori', 'sal.ori', 'sal.des', 'replicate', 'com', 'm.speed', 'sd.speed')

des.mo <- des.speed[des.speed$com == 'mono',]
des.mi <- des.speed[des.speed$com == 'mixed',]

mm.des.mo.id <- tapply(des.mo$m.speed, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), mean)
sd.des.mo.id <- tapply(des.mo$m.speed, paste(des.mo$sal.ori, des.mo$sal.des, des.mo$ID_ori), sd)
mm.des.mo.id <- data.frame(cbind(RTD(mm.des.mo.id), mm.des.mo.id, sd.des.mo.id))
colnames(mm.des.mo.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.speed', 'sd.speed')

mm.des.mi.id <- tapply(des.mi$m.speed, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), mean)
sd.des.mi.id <- tapply(des.mi$m.speed, paste(des.mi$sal.ori, des.mi$sal.des, des.mi$ID_ori), sd)
mm.des.mi.id <- data.frame(cbind(RTD(mm.des.mi.id), mm.des.mi.id, sd.des.mi.id))
colnames(mm.des.mi.id) <- c('sal.ori', 'sal.des', 'ID_ori', 'm.speed', 'sd.speed')

mm.des.mo <- tapply(des.mo$m.speed, paste(des.mo$sal.ori, des.mo$sal.des), mean)
sd.des.mo <- tapply(des.mo$m.speed, paste(des.mo$sal.ori, des.mo$sal.des), sd)
mm.des.mo <- data.frame(cbind(RTD(mm.des.mo), mm.des.mo, sd.des.mo))
colnames(mm.des.mo) <- c('sal.ori', 'sal.des', 'm.speed', 'sd.speed')

mm.des.mi <- tapply(des.mi$m.speed, paste(des.mi$sal.ori, des.mi$sal.des), mean)
sd.des.mi <- tapply(des.mi$m.speed, paste(des.mi$sal.ori, des.mi$sal.des), sd)
mm.des.mi <- data.frame(cbind(RTD(mm.des.mi), mm.des.mi, sd.des.mi))
colnames(mm.des.mi) <- c('sal.ori', 'sal.des', 'm.speed', 'sd.speed')

## Linear reaction norm figure
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(150, 820), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - speed') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 0])){
			segments(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]-id.tmp$sd.speed[id.tmp$sal.des == 0],-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]+id.tmp$sd.speed[id.tmp$sal.des == 0])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			points(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 0.5])){
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]-id.tmp$sd.speed[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]+id.tmp$sd.speed[id.tmp$sal.des == 0.5])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 1]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 1])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]-id.tmp$sd.speed[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]+id.tmp$sd.speed[id.tmp$sal.des == 1])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 2]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2], pch = 22, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 2])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]-id.tmp$sd.speed[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]+id.tmp$sd.speed[id.tmp$sal.des == 2])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1])
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.speed[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4], pch = 24, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 4])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4]-id.tmp$sd.speed[id.tmp$sal.des == 4], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4]+id.tmp$sd.speed[id.tmp$sal.des == 4])
			}
			x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
			if(length(x) > 0){	
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mo[mm.des.mo$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.speed[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.speed[tmp.$sal.des == 0]-tmp.$sd.speed[tmp.$sal.des == 0],-0.2, tmp.$m.speed[tmp.$sal.des == 0]+tmp.$sd.speed[tmp.$sal.des == 0], col= 'red')
	}
	points(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]+tmp.$sd.speed[tmp.$sal.des == 0.5], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]+tmp.$sd.speed[tmp.$sal.des == 1], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]+tmp.$sd.speed[tmp.$sal.des == 2], col = 'red')
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.speed[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]+tmp.$sd.speed[tmp.$sal.des == 4], col = 'red')
		segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], lwd = 2)
	}
}

## competition
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(150, 820), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Spite - speed') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 0])){
			segments(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]-id.tmp$sd.speed[id.tmp$sal.des == 0],-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]+id.tmp$sd.speed[id.tmp$sal.des == 0])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
		if(length(x) > 0){	
			points(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 0.5])){
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]-id.tmp$sd.speed[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]+id.tmp$sd.speed[id.tmp$sal.des == 0.5])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 1]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1], pch = 21, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 1])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]-id.tmp$sd.speed[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]+id.tmp$sd.speed[id.tmp$sal.des == 1])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 2]
		if(length(x) > 0){	
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2], pch = 22, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 2])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]-id.tmp$sd.speed[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]+id.tmp$sd.speed[id.tmp$sal.des == 2])
			}
		}
		x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
		y <- id.tmp$m.speed[id.tmp$sal.des == 1]
		z <- id.tmp$m.speed[id.tmp$sal.des == 2]
		if(length(x) > 0 & length(y) > 0){	
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1])
		}
		if(length(x) > 0 & length(z) > 0){	
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2])
		}

		x <- id.tmp$m.speed[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4], pch = 24, cex = 1)
			if(!is.na(id.tmp$sd.speed[id.tmp$sal.des == 4])){
			segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4]-id.tmp$sd.speed[id.tmp$sal.des == 4], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4]+id.tmp$sd.speed[id.tmp$sal.des == 4])
			}
			x <- id.tmp$m.speed[id.tmp$sal.des == 0.5]
			y <- id.tmp$m.speed[id.tmp$sal.des == 4]
			if(length(x) > 0 & length(y) > 0){	
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4])
			}
		}
	}
}
for(i in 1:3){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.speed[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.speed[tmp.$sal.des == 0]-tmp.$sd.speed[tmp.$sal.des == 0],-0.2, tmp.$m.speed[tmp.$sal.des == 0]+tmp.$sd.speed[tmp.$sal.des == 0], col= 'red')
	}
	points(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]+tmp.$sd.speed[tmp.$sal.des == 0.5], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]+tmp.$sd.speed[tmp.$sal.des == 1], col = 'red')
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]+tmp.$sd.speed[tmp.$sal.des == 2], col = 'red')
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.speed[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]+tmp.$sd.speed[tmp.$sal.des == 4], col = 'red')
		segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], lwd = 2)
	}
}

### -------------------------------------------------------------------------------------------------
## Bootstrap robustness analysis -- Figure S4
### -------------------------------------------------------------------------------------------------

dat09.combined.Spite <- rbind(dat09.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Spite$Salinity_Origin==0 & dat09.combined.Spite$Salinity_Destination==0)
pos.X42 <- which(dat09.combined.Spite$Salinity_Origin==4 | dat09.combined.Spite$Salinity_Origin==2)
w0.dat09.combined.Spite <- dat09.combined.Spite[-c(pos.0, pos.X42),]
w0.dat09.combined.Spite$Community <- factor(w0.dat09.combined.Spite$Community, levels=c('mono','mixed'))
dim(w0.dat09.combined.Spite)

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.area)

fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.ar)

fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Spite)
summary(fit.speed)

## make the data on which we will do the bootstrapping
boot.data <- w0.dat09.combined.Spite

## create a variable that contains the bootstrap categories
tt0 <- paste(w0.dat09.combined.Spite$Community,w0.dat09.combined.Spite$Salinity_Origin,w0.dat09.combined.Spite$Salinity_Destination,w0.dat09.combined.Spite$ID_original,w0.dat09.combined.Spite$replicate)
boot.data$boot.cat <- tt0

## what are the unique categories of the bootstrap
boot.lvl <- unique(tt0)

## how much do you want to select
bb <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

## how much do you want to bootstrap
perm <- 1000

boot.area.coef <- boot.area.p <- array(NA, dim=c(perm,10,length(bb)))
dimnames(boot.area.coef) <- dimnames(boot.area.p) <- list(NULL, c('Intercept', 'Salinity_Origin', 'Salinity_Destination', 'Community', 'Bio.Fraction', 'Density', 'Sal_oriXSal_des', 'Sal_oriXCom', 'Sal_DesXCom', 'Sal_OriXSal_DesXCom'), NULL)
boot.ar.coef <- boot.ar.p <- array(NA, dim=c(perm,10,length(bb)))
dimnames(boot.ar.coef) <- dimnames(boot.ar.p) <- list(NULL, c('Intercept', 'Salinity_Origin', 'Salinity_Destination', 'Community', 'Bio.Fraction', 'Density', 'Sal_oriXSal_des', 'Sal_oriXCom', 'Sal_DesXCom', 'Sal_OriXSal_DesXCom'), NULL)
boot.speed.coef <- boot.speed.p <- array(NA, dim=c(perm,10,length(bb)))
dimnames(boot.speed.coef) <- dimnames(boot.speed.p) <- list(NULL, c('Intercept', 'Salinity_Origin', 'Salinity_Destination', 'Community', 'Bio.Fraction', 'Density', 'Sal_oriXSal_des', 'Sal_oriXCom', 'Sal_DesXCom', 'Sal_OriXSal_DesXCom'), NULL)

for(xx in 1:length(bb)){

	print(paste('data used -', bb[xx]*100, '%'))

	pb <- txtProgressBar(min = 0, max = perm, style = 3)
	for(kk in 1:perm){
		boot.tmp <- vector()
		for(i in 1:length(boot.lvl)){
			pos.cat <- which(boot.data$boot.cat == boot.lvl[i])

			n.boot <- ceiling(length(pos.cat)*bb[xx])
			boot.pos.tmp <- sample(pos.cat, n.boot, replace=F)
			boot.tmp <- rbind(boot.tmp, boot.data[boot.pos.tmp,])
		}

		boot.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=boot.tmp)
		sum.area <- summary(boot.area)
		coef.area <- sum.area$coefficients

		boot.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=boot.tmp)
		sum.ar <- summary(boot.ar)
		coef.ar <- sum.ar$coefficients

		boot.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=boot.tmp)
		sum.speed <- summary(boot.speed)
		coef.speed <- sum.speed$coefficients
	
		for(j in 1:nrow(coef.area)){
			boot.area.coef[kk, j, xx] <- coef.area[j, 1]
			boot.area.p[kk, j, xx] <- coef.area[j, 5]
			boot.ar.coef[kk, j, xx] <- coef.ar[j, 1]
			boot.ar.p[kk, j, xx] <- coef.ar[j, 5]
			boot.speed.coef[kk, j, xx] <- coef.speed[j, 1]
			boot.speed.p[kk, j, xx] <- coef.speed[j, 5]
		}

		Sys.sleep(0.1)
		# update progress bar
		setTxtProgressBar(pb, kk)	
	}
	close(pb)
}

boot.area.coef

## Plot the effect sizes along the different percentage data used
## Also plot the p-values

## area
pathway <- #"Give in Pathway"
pdf(paste(pathway, 'EffectSize_Spite_Coef_area.pdf', sep=""))

for(kk in 2:10){
	x1 <- boot.area.coef[,kk,]
	m.x1 <- apply(x1, 2, mean)
	q.x1 <- apply(x1, 2, function(x) quantile(x, c(0.025, 0.975)))

	x2 <- boot.area.p[,kk,]
	m.x2 <- apply(x2, 2, mean)
	q.x2 <- apply(x2, 2, function(x) quantile(x, c(0.025, 0.975)))

	if(min(q.x1) < 0){
		ymin.x1 <- min(q.x1)+0.1*min(q.x1)
	} else {
		ymin.x1 <- min(q.x1)-0.1*min(q.x1)
	}
	if(max(q.x1) < 0){
		ymax.x1 <- max(q.x1)-0.1*max(q.x1)
	} else {
		ymax.x1 <- max(q.x1)+0.1*max(q.x1)
	}

	par(mar=c(5.1,4.1,4.1,5.1))

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x1, ymax.x1), type = 'n', main=colnames(boot.area.coef)[kk])
	abline(h=0)
	for(i in 1:length(bb)){
		points(bb[i], m.x1[i], pch=22, cex=2.5)
		segments(bb[i], q.x1[1,i], bb[i], q.x1[2,i])
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i], m.x1[i], bb[i+1], m.x1[i+1])
	}
	points(1, summary(fit.area)$coefficients[kk,1], pch = 22, cex=2.5)
	segments(0.9, m.x1[9], 1, summary(fit.area)$coefficients[kk,1])

	if(min(q.x2) < 0){
		ymin.x2 <- min(q.x2)+0.1*min(q.x2)
	} else {
		ymin.x2 <- min(q.x2)-0.1*min(q.x2)
	}
	if(max(q.x2) < 0){
		ymax.x2 <- max(q.x2)-0.1*max(q.x2)
	} else {
		ymax.x2 <- max(q.x2)+0.1*max(q.x2)
	}
	if(ymin.x2 < 0.05 & ymax.x2 < 0.05){
		ymin.x2 <- 0
		ymax.x2 <- 0.051
	}

	## Allow a second plot on the same graph
	par(new=TRUE)

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x2, ymax.x2), type = 'n', axes='FALSE', xlab="", ylab="")
	for(i in 1:length(bb)){
		points(bb[i]+0.01, m.x2[i], pch = 21, cex =2.5, bg='red')
		segments(bb[i]+0.01, q.x2[1,i], bb[i]+0.01, q.x2[2,i], col='red')
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i]+0.01, m.x2[i], bb[i+1]+0.01, m.x2[i+1], col='red')
	}		
	mtext("p-value",side=4,col="red",line=4) 
	axis(4, seq(ymin.x2, ymax.x2, length.out=5), col="red",col.axis="red",las=1)
	abline(h=0.05, col='red')
	points(1+0.01, summary(fit.area)$coefficients[kk,5], pch = 21, bg='red', cex=2.5)
	segments(0.9+0.01, m.x2[9], 1+0.01, summary(fit.area)$coefficients[kk,5], col='red')

	axis(1, seq(0,1, by=0.1)) 
}

dev.off()

## ar
pathway <- #"Give in Pathway"
pdf(paste(pathway, 'EffectSize_Spite_Coef_ar.pdf', sep=""))

for(kk in 2:10){
	x1 <- boot.ar.coef[,kk,]
	m.x1 <- apply(x1, 2, mean)
	q.x1 <- apply(x1, 2, function(x) quantile(x, c(0.025, 0.975)))

	x2 <- boot.ar.p[,kk,]
	m.x2 <- apply(x2, 2, mean)
	q.x2 <- apply(x2, 2, function(x) quantile(x, c(0.025, 0.975)))

	if(min(q.x1) < 0){
		ymin.x1 <- min(q.x1)+0.1*min(q.x1)
	} else {
		ymin.x1 <- min(q.x1)-0.1*min(q.x1)
	}
	if(max(q.x1) < 0){
		ymax.x1 <- max(q.x1)-0.1*max(q.x1)
	} else {
		ymax.x1 <- max(q.x1)+0.1*max(q.x1)
	}

	par(mar=c(5.1,4.1,4.1,5.1))

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x1, ymax.x1), type = 'n', main=colnames(boot.ar.coef)[kk])
	abline(h=0)
	for(i in 1:length(bb)){
		points(bb[i], m.x1[i], pch=22, cex=2.5)
		segments(bb[i], q.x1[1,i], bb[i], q.x1[2,i])
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i], m.x1[i], bb[i+1], m.x1[i+1])
	}
	points(1, summary(fit.ar)$coefficients[kk,1], pch = 22, cex=2.5)
	segments(0.9, m.x1[9], 1, summary(fit.ar)$coefficients[kk,1])

	if(min(q.x2) < 0){
		ymin.x2 <- min(q.x2)+0.1*min(q.x2)
	} else {
		ymin.x2 <- min(q.x2)-0.1*min(q.x2)
	}
	if(max(q.x2) < 0){
		ymax.x2 <- max(q.x2)-0.1*max(q.x2)
	} else {
		ymax.x2 <- max(q.x2)+0.1*max(q.x2)
	}
	if(ymin.x2 < 0.05 & ymax.x2 < 0.05){
		ymin.x2 <- 0
		ymax.x2 <- 0.051
	}

	## Allow a second plot on the same graph
	par(new=TRUE)

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x2, ymax.x2), type = 'n', axes='FALSE', xlab="", ylab="")
	for(i in 1:length(bb)){
		points(bb[i]+0.01, m.x2[i], pch = 21, cex =2.5, bg='red')
		segments(bb[i]+0.01, q.x2[1,i], bb[i]+0.01, q.x2[2,i], col='red')
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i]+0.01, m.x2[i], bb[i+1]+0.01, m.x2[i+1], col='red')
	}		
	mtext("p-value",side=4,col="red",line=4) 
	axis(4, seq(ymin.x2, ymax.x2, length.out=5), col="red",col.axis="red",las=1)
	abline(h=0.05, col='red')
	points(1+0.01, summary(fit.ar)$coefficients[kk,5], pch = 21, bg='red', cex=2.5)
	segments(0.9+0.01, m.x2[9], 1+0.01, summary(fit.ar)$coefficients[kk,5], col='red')

	axis(1, seq(0,1, by=0.1)) 
}

dev.off()


## speed
pathway <- #"Give in Pathway"
pdf(paste(pathway, 'EffectSize_Spite_Coef_speed.pdf', sep=""))

for(kk in 2:10){
	x1 <- boot.speed.coef[,kk,]
	m.x1 <- apply(x1, 2, mean)
	q.x1 <- apply(x1, 2, function(x) quantile(x, c(0.025, 0.975)))

	x2 <- boot.speed.p[,kk,]
	m.x2 <- apply(x2, 2, mean)
	q.x2 <- apply(x2, 2, function(x) quantile(x, c(0.025, 0.975)))

	if(min(q.x1) < 0){
		ymin.x1 <- min(q.x1)+0.1*min(q.x1)
	} else {
		ymin.x1 <- min(q.x1)-0.1*min(q.x1)
	}
	if(max(q.x1) < 0){
		ymax.x1 <- max(q.x1)-0.1*max(q.x1)
	} else {
		ymax.x1 <- max(q.x1)+0.1*max(q.x1)
	}

	par(mar=c(5.1,4.1,4.1,5.1))

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x1, ymax.x1), type = 'n', main=colnames(boot.speed.coef)[kk])
	abline(h=0)
	for(i in 1:length(bb)){
		points(bb[i], m.x1[i], pch=22, cex=2.5)
		segments(bb[i], q.x1[1,i], bb[i], q.x1[2,i])
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i], m.x1[i], bb[i+1], m.x1[i+1])
	}
	points(1, summary(fit.speed)$coefficients[kk,1], pch = 22, cex=2.5)
	segments(0.9, m.x1[9], 1, summary(fit.speed)$coefficients[kk,1])

	if(min(q.x2) < 0){
		ymin.x2 <- min(q.x2)+0.1*min(q.x2)
	} else {
		ymin.x2 <- min(q.x2)-0.1*min(q.x2)
	}
	if(max(q.x2) < 0){
		ymax.x2 <- max(q.x2)-0.1*max(q.x2)
	} else {
		ymax.x2 <- max(q.x2)+0.1*max(q.x2)
	}
	if(ymin.x2 < 0.05 & ymax.x2 < 0.05){
		ymin.x2 <- 0
		ymax.x2 <- 0.051
	}

	## Allow a second plot on the same graph
	par(new=TRUE)

	plot(bb, rep(-10000, length(bb)), xlim=c(bb[1]-0.1, bb[length(bb)]+0.1), ylim= c(ymin.x2, ymax.x2), type = 'n', axes='FALSE', xlab="", ylab="")
	for(i in 1:length(bb)){
		points(bb[i]+0.01, m.x2[i], pch = 21, cex =2.5, bg='red')
		segments(bb[i]+0.01, q.x2[1,i], bb[i]+0.01, q.x2[2,i], col='red')
	}
	for(i in 1:(length(bb)-1)){
		segments(bb[i]+0.01, m.x2[i], bb[i+1]+0.01, m.x2[i+1], col='red')
	}		
	mtext("p-value",side=4,col="red",line=4) 
	axis(4, seq(ymin.x2, ymax.x2, length.out=5), col="red",col.axis="red",las=1)
	abline(h=0.05, col='red')
	points(1+0.01, summary(fit.speed)$coefficients[kk,5], pch = 21, bg='red', cex=2.5)
	segments(0.9+0.01, m.x2[9], 1+0.01, summary(fit.speed)$coefficients[kk,5], col='red')

	axis(1, seq(0,1, by=0.1)) 
}

dev.off()



### -------------------------------------------------------------------------------------------------
## Reaction norms analysis on selected populations -- Table S7-S8
### -------------------------------------------------------------------------------------------------

## area
## 0. Try to fit the reaction norm approach to each selection population
sal.lvl <- c(0, 0.5, 1, 2, 4)

### Anc 0.5-1 vs All Des (0.5-1)
dat22.Spite.1 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 1,])

lmer.area.mo.X.1 <- lmer.area.mi.X.1 <- aov.area.mo.X.1 <- aov.area.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Spite.X.1 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 1 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.area.mo <- c(dat22.Spite.1$mean_area, dat09.Spite.X.1$mean_area)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.1$mean_area)), dat09.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.1$Salinity, dat09.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Spite.1$Density, dat09.Spite.X.1$Density)
	tt.id.mo <- c(dat22.Spite.1$ID, dat09.Spite.X.1$ID_new)

	lmer.mo <- lmer(tt.area.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.area.mo.X.1 <- rbind(lmer.area.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.area.mo.X.1 <- rbind(aov.area.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.1 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 1 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.area.mi <- c(dat22.Spite.1$mean_area, dat09.PTS.Spite.X.1$mean_area)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.1$mean_area)), dat09.PTS.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.1$Salinity, dat09.PTS.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Spite.1$Density, dat09.PTS.Spite.X.1$Density)
	tt.id.mi <- c(dat22.Spite.1$ID, dat09.PTS.Spite.X.1$ID_new)

	lmer.mi <- lmer(tt.area.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.area.mi.X.1 <- rbind(lmer.area.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.area.mi.X.1 <- rbind(aov.area.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.area.mo.X.1 <- data.frame(lmer.area.mo.X.1) ## output monocultures Table S7
lmer.area.mi.X.1 <- data.frame(lmer.area.mi.X.1) ## output competition Table S8
aov.area.mo.X.1 <- data.frame(aov.area.mo.X.1)
aov.area.mi.X.1 <- data.frame(aov.area.mi.X.1)
colnames(lmer.area.mo.X.1) <- colnames(lmer.area.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.area.mo.X.1) <- colnames(aov.area.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')


### Anc 0.5-2 vs All Des (0.5-2)
dat22.Spite.2 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 2,])

lmer.area.mo.X.2 <- lmer.area.mi.X.2 <- aov.area.mo.X.2 <- aov.area.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Spite.X.2 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 2 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.area.mo <- c(dat22.Spite.2$mean_area, dat09.Spite.X.2$mean_area)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.2$mean_area)), dat09.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.2$Salinity, dat09.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Spite.2$Density, dat09.Spite.X.2$Density)
	tt.id.mo <- c(dat22.Spite.2$ID, dat09.Spite.X.2$ID_new)

	lmer.mo <- lmer(tt.area.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.area.mo.X.2 <- rbind(lmer.area.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.area.mo.X.2 <- rbind(aov.area.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.2 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 2 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.area.mi <- c(dat22.Spite.2$mean_area, dat09.PTS.Spite.X.2$mean_area)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.2$mean_area)), dat09.PTS.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.2$Salinity, dat09.PTS.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Spite.2$Density, dat09.PTS.Spite.X.2$Density)
	tt.id.mi <- c(dat22.Spite.2$ID, dat09.PTS.Spite.X.2$ID_new)

	lmer.mi <- lmer(tt.area.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.area.mi.X.2 <- rbind(lmer.area.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.area.mi.X.2 <- rbind(aov.area.mi.X.2, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.area.mo.X.2 <- data.frame(lmer.area.mo.X.2)
lmer.area.mi.X.2 <- data.frame(lmer.area.mi.X.2)
aov.area.mo.X.2 <- data.frame(aov.area.mo.X.2)
aov.area.mi.X.2 <- data.frame(aov.area.mi.X.2)
colnames(lmer.area.mo.X.2) <- colnames(lmer.area.mi.X.2) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.area.mo.X.2) <- colnames(aov.area.mi.X.2) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')

## Cell shape 
## 0. Try to fit the reaction norm approach to each selection population
sal.lvl <- c(0, 0.5, 1, 2, 4)

### Anc 0.5-1 vs All Des (0.5-1)
dat22.Spite.1 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 1,])

lmer.ar.mo.X.1 <- lmer.ar.mi.X.1 <- aov.ar.mo.X.1 <- aov.ar.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Spite.X.1 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 1 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.ar.mo <- c(dat22.Spite.1$mean_ar, dat09.Spite.X.1$mean_ar)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.1$mean_ar)), dat09.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.1$Salinity, dat09.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Spite.1$Density, dat09.Spite.X.1$Density)
	tt.id.mo <- c(dat22.Spite.1$ID, dat09.Spite.X.1$ID_new)

	lmer.mo <- lmer(tt.ar.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.ar.mo.X.1 <- rbind(lmer.ar.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,3], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.ar.mo.X.1 <- rbind(aov.ar.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.1 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 1 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.ar.mi <- c(dat22.Spite.1$mean_ar, dat09.PTS.Spite.X.1$mean_ar)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.1$mean_ar)), dat09.PTS.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.1$Salinity, dat09.PTS.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Spite.1$Density, dat09.PTS.Spite.X.1$Density)
	tt.id.mi <- c(dat22.Spite.1$ID, dat09.PTS.Spite.X.1$ID_new)

	lmer.mi <- lmer(tt.ar.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.ar.mi.X.1 <- rbind(lmer.ar.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,3], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.ar.mi.X.1 <- rbind(aov.area.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.ar.mo.X.1 <- data.frame(lmer.ar.mo.X.1)
lmer.ar.mi.X.1 <- data.frame(lmer.ar.mi.X.1)
aov.ar.mo.X.1 <- data.frame(aov.ar.mo.X.1)
aov.ar.mi.X.1 <- data.frame(aov.ar.mi.X.1)
colnames(lmer.ar.mo.X.1) <- colnames(lmer.ar.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.ar.mo.X.1) <- colnames(aov.ar.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')


### Anc 0.5-2 vs All Des (0.5-2)
dat22.Spite.2 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 2,])

lmer.ar.mo.X.2 <- lmer.ar.mi.X.2 <- aov.ar.mo.X.2 <- aov.ar.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Spite.X.2 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 2 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.ar.mo <- c(dat22.Spite.2$mean_ar, dat09.Spite.X.2$mean_ar)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.2$mean_ar)), dat09.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.2$Salinity, dat09.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Spite.2$Density, dat09.Spite.X.2$Density)
	tt.id.mo <- c(dat22.Spite.2$ID, dat09.Spite.X.2$ID_new)

	lmer.mo <- lmer(tt.ar.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.ar.mo.X.2 <- rbind(lmer.ar.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,3], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.ar.mo.X.2 <- rbind(aov.ar.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.2 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 2 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.ar.mi <- c(dat22.Spite.2$mean_ar, dat09.PTS.Spite.X.2$mean_ar)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.2$mean_ar)), dat09.PTS.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.2$Salinity, dat09.PTS.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Spite.2$Density, dat09.PTS.Spite.X.2$Density)
	tt.id.mi <- c(dat22.Spite.2$ID, dat09.PTS.Spite.X.2$ID_new)

	lmer.mi <- lmer(tt.ar.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.ar.mi.X.2 <- rbind(lmer.ar.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,3], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.ar.mi.X.2 <- rbind(aov.ar.mi.X.2, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.ar.mo.X.2 <- data.frame(lmer.ar.mo.X.2)
lmer.ar.mi.X.2 <- data.frame(lmer.ar.mi.X.2)
aov.ar.mo.X.2 <- data.frame(aov.ar.mo.X.2)
aov.ar.mi.X.2 <- data.frame(aov.ar.mi.X.2)
colnames(lmer.ar.mo.X.2) <- colnames(lmer.ar.mi.X.2) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.ar.mo.X.2) <- colnames(aov.ar.mi.X.2) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')

## speed
## 0. Try to fit the reaction norm approach to each selection population
sal.lvl <- c(0, 0.5, 1, 2, 4)

### Anc 0.5-1 vs All Des (0.5-1)
dat22.Spite.1 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 1,])

lmer.speed.mo.X.1 <- lmer.speed.mi.X.1 <- aov.speed.mo.X.1 <- aov.speed.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Spite.X.1 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 1 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.speed.mo <- c(dat22.Spite.1$gross_speed, dat09.Spite.X.1$gross_speed)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.1$gross_speed)), dat09.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.1$Salinity, dat09.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Spite.1$Density, dat09.Spite.X.1$Density)
	tt.id.mo <- c(dat22.Spite.1$ID, dat09.Spite.X.1$ID_new)

	lmer.mo <- lmer(tt.speed.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.speed.mo.X.1 <- rbind(lmer.speed.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.speed.mo.X.1 <- rbind(aov.ar.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.1 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 1 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.speed.mi <- c(dat22.Spite.1$gross_speed, dat09.PTS.Spite.X.1$gross_speed)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.1$gross_speed)), dat09.PTS.Spite.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.1$Salinity, dat09.PTS.Spite.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Spite.1$Density, dat09.PTS.Spite.X.1$Density)
	tt.id.mi <- c(dat22.Spite.1$ID, dat09.PTS.Spite.X.1$ID_new)

	lmer.mi <- lmer(tt.speed.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.speed.mi.X.1 <- rbind(lmer.speed.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.speed.mi.X.1 <- rbind(aov.speed.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.speed.mo.X.1 <- data.frame(lmer.speed.mo.X.1)
lmer.speed.mi.X.1 <- data.frame(lmer.speed.mi.X.1)
aov.speed.mo.X.1 <- data.frame(aov.speed.mo.X.1)
aov.speed.mi.X.1 <- data.frame(aov.speed.mi.X.1)
colnames(lmer.speed.mo.X.1) <- colnames(lmer.speed.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.speed.mo.X.1) <- colnames(aov.speed.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')


### Anc 0.5-2 vs All Des (0.5-2)
dat22.Spite.2 <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == 2,])

lmer.speed.mo.X.2 <- lmer.speed.mi.X.2 <- aov.speed.mo.X.2 <- aov.speed.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Spite.X.2 <- rbind(dat09.Spite[dat09.Spite$Salinity_Destination == 0.5 & dat09.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.Spite[dat09.Spite$Salinity_Destination == 2 & dat09.Spite$Salinity_Origin == sal.lvl[sal],])

	tt.speed.mo <- c(dat22.Spite.2$gross_speed, dat09.Spite.X.2$gross_speed)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Spite.2$gross_speed)), dat09.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Spite.2$Salinity, dat09.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Spite.2$Density, dat09.Spite.X.2$Density)
	tt.id.mo <- c(dat22.Spite.2$ID, dat09.Spite.X.2$ID_new)

	lmer.mo <- lmer(tt.speed.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.speed.mo.X.2 <- rbind(lmer.speed.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.speed.mo.X.2 <- rbind(aov.speed.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
}
for(sal in 1:3){
	dat09.PTS.Spite.X.2 <- rbind(dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 0.5 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Destination == 2 & dat09.PTS.Spite$Salinity_Origin == sal.lvl[sal],])
	
	tt.speed.mi <- c(dat22.Spite.2$gross_speed, dat09.PTS.Spite.X.2$gross_speed)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Spite.2$gross_speed)), dat09.PTS.Spite.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Spite.2$Salinity, dat09.PTS.Spite.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Spite.2$Density, dat09.PTS.Spite.X.2$Density)
	tt.id.mi <- c(dat22.Spite.2$ID, dat09.PTS.Spite.X.2$ID_new)

	lmer.mi <- lmer(tt.speed.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.speed.mi.X.2 <- rbind(lmer.speed.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.speed.mi.X.2 <- rbind(aov.speed.mi.X.2, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.speed.mo.X.2 <- data.frame(lmer.speed.mo.X.2)
lmer.speed.mi.X.2 <- data.frame(lmer.speed.mi.X.2)
aov.speed.mo.X.2 <- data.frame(aov.speed.mo.X.2)
aov.speed.mi.X.2 <- data.frame(aov.speed.mi.X.2)
colnames(lmer.speed.mo.X.2) <- colnames(lmer.speed.mi.X.2) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.speed.mo.X.2) <- colnames(aov.speed.mi.X.2) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')



### -------------------------------------------------------------------------------------------------
## Reaction norms analysis on selected populations -- Figure S6 --> mean trait components are use for Figure 4 Main text
### -------------------------------------------------------------------------------------------------

## plot the estimates!
## Using barplots

## area
dev.new(height = 14, width = 21)
plast. <- c(lmer.area.mo.X.1$plast, lmer.area.mo.X.2$plast, lmer.area.mi.X.1$plast, lmer.area.mi.X.2$plast)
gen. <- c(lmer.area.mo.X.1$evo, lmer.area.mo.X.2$evo, lmer.area.mi.X.1$evo, lmer.area.mi.X.2$evo)
ep. <- c(lmer.area.mo.X.1$evoplast, lmer.area.mo.X.2$evoplast, lmer.area.mi.X.1$evoplast, lmer.area.mi.X.2$evoplast)
se.plast. <- c(lmer.area.mo.X.1$se.plast, lmer.area.mo.X.2$se.plast, lmer.area.mi.X.1$se.plast, lmer.area.mi.X.2$se.plast)
se.gen. <- c(lmer.area.mo.X.1$se.evo, lmer.area.mo.X.2$se.evo, lmer.area.mi.X.1$se.evo, lmer.area.mi.X.2$se.evo)
se.ep. <- c(lmer.area.mo.X.1$se.evoplast, lmer.area.mo.X.2$se.evoplast, lmer.area.mi.X.1$se.evoplast, lmer.area.mi.X.2$se.evoplast)
ymin <- min(c(plast., gen., ep.)-c(se.plast., se.gen., se.ep.))
ymax <- max(c(plast., gen., ep.)+c(se.plast., se.gen., se.ep.))
par(mfrow = c(4,5))
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.area.mo.X.1$plast[sal])
	segments(0.15, lmer.area.mo.X.1$plast[sal]-lmer.area.mo.X.1$se.plast[sal], 0.15, lmer.area.mo.X.1$plast[sal]+lmer.area.mo.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.area.mo.X.1$evo[sal])
	segments(0.5, lmer.area.mo.X.1$evo[sal]-lmer.area.mo.X.1$se.evo[sal], 0.5, lmer.area.mo.X.1$evo[sal]+lmer.area.mo.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.area.mo.X.1$evoplast[sal])
	segments(0.85, lmer.area.mo.X.1$evoplast[sal]-lmer.area.mo.X.1$se.evoplast[sal], 0.85, lmer.area.mo.X.1$evoplast[sal]+lmer.area.mo.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.area.mo.X.2$plast[sal])
	segments(0.15, lmer.area.mo.X.2$plast[sal]-lmer.area.mo.X.2$se.plast[sal], 0.15, lmer.area.mo.X.2$plast[sal]+lmer.area.mo.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.area.mo.X.2$evo[sal])
	segments(0.5, lmer.area.mo.X.2$evo[sal]-lmer.area.mo.X.2$se.evo[sal], 0.5, lmer.area.mo.X.2$evo[sal]+lmer.area.mo.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.area.mo.X.2$evoplast[sal])
	segments(0.85, lmer.area.mo.X.2$evoplast[sal]-lmer.area.mo.X.2$se.evoplast[sal], 0.85, lmer.area.mo.X.2$evoplast[sal]+lmer.area.mo.X.2$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.area.mi.X.1$plast[sal])
	segments(0.15, lmer.area.mi.X.1$plast[sal]-lmer.area.mi.X.1$se.plast[sal], 0.15, lmer.area.mi.X.1$plast[sal]+lmer.area.mi.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.area.mi.X.1$evo[sal])
	segments(0.5, lmer.area.mi.X.1$evo[sal]-lmer.area.mi.X.1$se.evo[sal], 0.5, lmer.area.mi.X.1$evo[sal]+lmer.area.mi.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.area.mi.X.1$evoplast[sal])
	segments(0.85, lmer.area.mi.X.1$evoplast[sal]-lmer.area.mi.X.1$se.evoplast[sal], 0.85, lmer.area.mi.X.1$evoplast[sal]+lmer.area.mi.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.area.mi.X.2$plast[sal])
	segments(0.15, lmer.area.mi.X.2$plast[sal]-lmer.area.mi.X.2$se.plast[sal], 0.15, lmer.area.mi.X.2$plast[sal]+lmer.area.mi.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.area.mi.X.2$evo[sal])
	segments(0.5, lmer.area.mi.X.2$evo[sal]-lmer.area.mi.X.2$se.evo[sal], 0.5, lmer.area.mi.X.2$evo[sal]+lmer.area.mi.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.area.mi.X.2$evoplast[sal])
	segments(0.85, lmer.area.mi.X.2$evoplast[sal]-lmer.area.mi.X.2$se.evoplast[sal], 0.85, lmer.area.mi.X.2$evoplast[sal]+lmer.area.mi.X.2$se.evoplast[sal])
}

## cell shape
dev.new(height = 14, width = 21)
plast. <- c(lmer.ar.mo.X.1$plast, lmer.ar.mo.X.2$plast, lmer.ar.mi.X.1$plast, lmer.ar.mi.X.2$plast)
gen. <- c(lmer.ar.mo.X.1$evo, lmer.ar.mo.X.2$evo, lmer.ar.mi.X.1$evo, lmer.ar.mi.X.2$evo)
ep. <- c(lmer.ar.mo.X.1$evoplast, lmer.ar.mo.X.2$evoplast, lmer.ar.mi.X.1$evoplast, lmer.ar.mi.X.2$evoplast)
se.plast. <- c(lmer.ar.mo.X.1$se.plast, lmer.ar.mo.X.2$se.plast, lmer.ar.mi.X.1$se.plast, lmer.ar.mi.X.2$se.plast)
se.gen. <- c(lmer.ar.mo.X.1$se.evo, lmer.ar.mo.X.2$se.evo, lmer.ar.mi.X.1$se.evo, lmer.ar.mi.X.2$se.evo)
se.ep. <- c(lmer.ar.mo.X.1$se.evoplast, lmer.ar.mo.X.2$se.evoplast, lmer.ar.mi.X.1$se.evoplast, lmer.ar.mi.X.2$se.evoplast)
ymin <- min(c(plast., gen., ep.)-c(se.plast., se.gen., se.ep.))
ymax <- max(c(plast., gen., ep.)+c(se.plast., se.gen., se.ep.))
par(mfrow = c(4,5))
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(-3.5, 2.5), main = 'comparison ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mo.X.1$plast[sal])
	segments(0.15, lmer.ar.mo.X.1$plast[sal]-lmer.ar.mo.X.1$se.plast[sal], 0.15, lmer.ar.mo.X.1$plast[sal]+lmer.ar.mo.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mo.X.1$evo[sal])
	segments(0.5, lmer.ar.mo.X.1$evo[sal]-lmer.ar.mo.X.1$se.evo[sal], 0.5, lmer.ar.mo.X.1$evo[sal]+lmer.ar.mo.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mo.X.1$evoplast[sal])
	segments(0.85, lmer.ar.mo.X.1$evoplast[sal]-lmer.ar.mo.X.1$se.evoplast[sal], 0.85, lmer.ar.mo.X.1$evoplast[sal]+lmer.ar.mo.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(-3.5, 2.5), main = 'comparison of ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mo.X.2$plast[sal])
	segments(0.15, lmer.ar.mo.X.2$plast[sal]-lmer.ar.mo.X.2$se.plast[sal], 0.15, lmer.ar.mo.X.2$plast[sal]+lmer.ar.mo.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mo.X.2$evo[sal])
	segments(0.5, lmer.ar.mo.X.2$evo[sal]-lmer.ar.mo.X.2$se.evo[sal], 0.5, lmer.ar.mo.X.2$evo[sal]+lmer.ar.mo.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mo.X.2$evoplast[sal])
	segments(0.85, lmer.ar.mo.X.2$evoplast[sal]-lmer.ar.mo.X.2$se.evoplast[sal], 0.85, lmer.ar.mo.X.2$evoplast[sal]+lmer.ar.mo.X.2$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(-3.5, 2.5), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mi.X.1$plast[sal])
	segments(0.15, lmer.ar.mi.X.1$plast[sal]-lmer.ar.mi.X.1$se.plast[sal], 0.15, lmer.ar.mi.X.1$plast[sal]+lmer.ar.mi.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mi.X.1$evo[sal])
	segments(0.5, lmer.ar.mi.X.1$evo[sal]-lmer.ar.mi.X.1$se.evo[sal], 0.5, lmer.ar.mi.X.1$evo[sal]+lmer.ar.mi.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mi.X.1$evoplast[sal])
	segments(0.85, lmer.ar.mi.X.1$evoplast[sal]-lmer.ar.mi.X.1$se.evoplast[sal], 0.85, lmer.ar.mi.X.1$evoplast[sal]+lmer.ar.mi.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(-3.5, 2.5), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mi.X.2$plast[sal])
	segments(0.15, lmer.ar.mi.X.2$plast[sal]-lmer.ar.mi.X.2$se.plast[sal], 0.15, lmer.ar.mi.X.2$plast[sal]+lmer.ar.mi.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mi.X.2$evo[sal])
	segments(0.5, lmer.ar.mi.X.2$evo[sal]-lmer.ar.mi.X.2$se.evo[sal], 0.5, lmer.ar.mi.X.2$evo[sal]+lmer.ar.mi.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mi.X.2$evoplast[sal])
	segments(0.85, lmer.ar.mi.X.2$evoplast[sal]-lmer.ar.mi.X.2$se.evoplast[sal], 0.85, lmer.ar.mi.X.2$evoplast[sal]+lmer.ar.mi.X.2$se.evoplast[sal])
}



## speed
dev.new(height = 14, width = 21)
plast. <- c(lmer.speed.mo.X.1$plast, lmer.speed.mo.X.2$plast, lmer.speed.mi.X.1$plast, lmer.speed.mi.X.2$plast)
gen. <- c(lmer.speed.mo.X.1$evo, lmer.speed.mo.X.2$evo, lmer.speed.mi.X.1$evo, lmer.speed.mi.X.2$evo)
ep. <- c(lmer.speed.mo.X.1$evoplast, lmer.speed.mo.X.2$evoplast, lmer.speed.mi.X.1$evoplast, lmer.speed.mi.X.2$evoplast)
se.plast. <- c(lmer.speed.mo.X.1$se.plast, lmer.speed.mo.X.2$se.plast, lmer.speed.mi.X.1$se.plast, lmer.speed.mi.X.2$se.plast)
se.gen. <- c(lmer.speed.mo.X.1$se.evo, lmer.speed.mo.X.2$se.evo, lmer.speed.mi.X.1$se.evo, lmer.speed.mi.X.2$se.evo)
se.ep. <- c(lmer.speed.mo.X.1$se.evoplast, lmer.speed.mo.X.2$se.evoplast, lmer.speed.mi.X.1$se.evoplast, lmer.speed.mi.X.2$se.evoplast)
ymin <- min(c(plast., gen., ep.)-c(se.plast., se.gen., se.ep.))
ymax <- max(c(plast., gen., ep.)+c(se.plast., se.gen., se.ep.))
par(mfrow = c(4,5))
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.speed.mo.X.1$plast[sal])
	segments(0.15, lmer.speed.mo.X.1$plast[sal]-lmer.speed.mo.X.1$se.plast[sal], 0.15, lmer.speed.mo.X.1$plast[sal]+lmer.speed.mo.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.speed.mo.X.1$evo[sal])
	segments(0.5, lmer.speed.mo.X.1$evo[sal]-lmer.speed.mo.X.1$se.evo[sal], 0.5, lmer.speed.mo.X.1$evo[sal]+lmer.speed.mo.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.speed.mo.X.1$evoplast[sal])
	segments(0.85, lmer.speed.mo.X.1$evoplast[sal]-lmer.speed.mo.X.1$se.evoplast[sal], 0.85, lmer.speed.mo.X.1$evoplast[sal]+lmer.speed.mo.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.speed.mo.X.2$plast[sal])
	segments(0.15, lmer.speed.mo.X.2$plast[sal]-lmer.speed.mo.X.2$se.plast[sal], 0.15, lmer.speed.mo.X.2$plast[sal]+lmer.speed.mo.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.speed.mo.X.2$evo[sal])
	segments(0.5, lmer.speed.mo.X.2$evo[sal]-lmer.speed.mo.X.2$se.evo[sal], 0.5, lmer.speed.mo.X.2$evo[sal]+lmer.speed.mo.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.speed.mo.X.2$evoplast[sal])
	segments(0.85, lmer.speed.mo.X.2$evoplast[sal]-lmer.speed.mo.X.2$se.evoplast[sal], 0.85, lmer.speed.mo.X.2$evoplast[sal]+lmer.speed.mo.X.2$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.speed.mi.X.1$plast[sal])
	segments(0.15, lmer.speed.mi.X.1$plast[sal]-lmer.speed.mi.X.1$se.plast[sal], 0.15, lmer.speed.mi.X.1$plast[sal]+lmer.speed.mi.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.speed.mi.X.1$evo[sal])
	segments(0.5, lmer.speed.mi.X.1$evo[sal]-lmer.speed.mi.X.1$se.evo[sal], 0.5, lmer.speed.mi.X.1$evo[sal]+lmer.speed.mi.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.speed.mi.X.1$evoplast[sal])
	segments(0.85, lmer.speed.mi.X.1$evoplast[sal]-lmer.speed.mi.X.1$se.evoplast[sal], 0.85, lmer.speed.mi.X.1$evoplast[sal]+lmer.speed.mi.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.speed.mi.X.2$plast[sal])
	segments(0.15, lmer.speed.mi.X.2$plast[sal]-lmer.speed.mi.X.2$se.plast[sal], 0.15, lmer.speed.mi.X.2$plast[sal]+lmer.speed.mi.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.speed.mi.X.2$evo[sal])
	segments(0.5, lmer.speed.mi.X.2$evo[sal]-lmer.speed.mi.X.2$se.evo[sal], 0.5, lmer.speed.mi.X.2$evo[sal]+lmer.speed.mi.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.speed.mi.X.2$evoplast[sal])
	segments(0.85, lmer.speed.mi.X.2$evoplast[sal]-lmer.speed.mi.X.2$se.evoplast[sal], 0.85, lmer.speed.mi.X.2$evoplast[sal]+lmer.speed.mi.X.2$se.evoplast[sal])
}



### -------------------------------------------------------------------------------------------------
## Phenotypic plasticity response to salinity -- Figure S8 and Table S11-S12
### -------------------------------------------------------------------------------------------------


## Calculating plasticity responses for the descendant population 
## Using regression with density

sal.lvl <- c(1, 2)
plast.anc.spite.mo.area.lm <- plast.anc.spite.mo.ar.lm <- plast.anc.spite.mo.speed.lm <- vector()
for(sal in 1:2){
	dat.tmp <- rbind(dat22.Spite.CG[dat22.Spite.CG$Salinity == 0.5,], dat22.Spite.CG[dat22.Spite.CG$Salinity == sal.lvl[sal],])
	Sal. <- factor(dat.tmp$Salinity, levels = c(0.5, sal.lvl[sal]))

	lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.area <- summary(lmer.area)$coefficients

	lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.ar <- summary(lmer.ar)$coefficients

	lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.speed <- summary(lmer.speed)$coefficients

	plast.anc.spite.mo.area.lm <- rbind(plast.anc.spite.mo.area.lm, cbind(sal.lvl[sal], sum.area[2,1], sum.area[2,2], sum.area[2,5], 
							sum.area[3,1], sum.area[3,2], sum.area[3,5]))
	plast.anc.spite.mo.ar.lm <- rbind(plast.anc.spite.mo.ar.lm, cbind(sal.lvl[sal], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5], 
							sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
	plast.anc.spite.mo.speed.lm <- rbind(plast.anc.spite.mo.speed.lm, cbind(sal.lvl[sal], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5], 
							sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
}
plast.anc.spite.mo.area.lm <- data.frame(plast.anc.spite.mo.area.lm)
plast.anc.spite.mo.ar.lm <- data.frame(plast.anc.spite.mo.ar.lm)
plast.anc.spite.mo.speed.lm <- data.frame(plast.anc.spite.mo.speed.lm)
colnames(plast.anc.spite.mo.area.lm) <- c('salinity', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.spite.mo.ar.lm) <- c('salinity', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.spite.mo.speed.lm) <- c('salinity', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')


sal.lvl <- c(1, 2)
plast.anc.spite.mi.area.lm <- plast.anc.spite.mi.ar.lm <- plast.anc.spite.mi.speed.lm <- vector()
for(sal in 1:2){
	dat.tmp <- rbind(dat22.PTS.Spite.CG[dat22.PTS.Spite.CG$Salinity == 0.5,], dat22.PTS.Spite.CG[dat22.PTS.Spite.CG$Salinity == sal.lvl[sal],])
	Sal. <- factor(dat.tmp$Salinity, levels = c(0.5, sal.lvl[sal]))

	lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.area <- summary(lmer.area)$coefficients

	lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.ar <- summary(lmer.ar)$coefficients

	lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.speed <- summary(lmer.speed)$coefficients

	plast.anc.spite.mi.area.lm <- rbind(plast.anc.spite.mi.area.lm, cbind(sal.lvl[sal], sum.area[2,1], sum.area[2,2], sum.area[2,5], 
							sum.area[3,1], sum.area[3,2], sum.area[3,5]))
	plast.anc.spite.mi.ar.lm <- rbind(plast.anc.spite.mi.ar.lm, cbind(sal.lvl[sal], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5], 
							sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
	plast.anc.spite.mi.speed.lm <- rbind(plast.anc.spite.mi.speed.lm, cbind(sal.lvl[sal], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5], 
							sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
}
plast.anc.spite.mi.area.lm <- data.frame(plast.anc.spite.mi.area.lm)
plast.anc.spite.mi.ar.lm <- data.frame(plast.anc.spite.mi.ar.lm)
plast.anc.spite.mi.speed.lm <- data.frame(plast.anc.spite.mi.speed.lm)
colnames(plast.anc.spite.mi.area.lm) <- c('salinity', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.spite.mi.ar.lm) <- c('salinity', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.spite.mi.speed.lm) <- c('salinity', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')



dat09.combined.Spite <- rbind(dat09.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

sal.ori <- c(0, 0.5, 1, 2, 4)
sal.des <- c(1, 2, 4)

plast.desc.spite.mo.area.lm <- plast.desc.spite.mo.ar.lm <- plast.desc.spite.mo.speed.lm <- vector()
for(ori in 1:5){
	dat.tmp <- dat09.Spite[dat09.Spite$Salinity_Origin == sal.ori[ori],]

	for(des in 1:3){
		nr <- nrow(dat.tmp[dat.tmp$Salinity_Destination == sal.des[des],])
		
		if(nr > 0){
			dat.tmp2 <- rbind(dat.tmp[dat.tmp$Salinity_Destination == 0.5,], dat.tmp[dat.tmp$Salinity_Destination == sal.des[des],])
			Sal. <- factor(dat.tmp2$Salinity_Destination, levels = c(0.5, sal.des[des]))

			lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.area <- summary(lmer.area)$coefficients

			
			lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.ar <- summary(lmer.ar)$coefficients

			lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.speed <- summary(lmer.speed)$coefficients	

			nr <- nrow(sum.area)

			if(nr == 3){
				plast.desc.spite.mo.area.lm <- rbind(plast.desc.spite.mo.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									sum.area[3,1], sum.area[3,2], sum.area[3,5]))
				plast.desc.spite.mo.ar.lm <- rbind(plast.desc.spite.mo.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
				plast.desc.spite.mo.speed.lm <- rbind(plast.desc.spite.mo.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
			}

			if(nr == 2){
				plast.desc.spite.mo.area.lm <- rbind(plast.desc.spite.mo.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									'NA', 'NA', 'NA'))
				plast.desc.spite.mo.ar.lm <- rbind(plast.desc.spite.mo.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									'NA', 'NA', 'NA'))
				plast.desc.spite.mo.speed.lm <- rbind(plast.desc.spite.mo.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									'NA', 'NA', 'NA'))
			}
		}
	}
}
plast.desc.spite.mo.area.lm <- data.frame(plast.desc.spite.mo.area.lm)
plast.desc.spite.mo.ar.lm <- data.frame(plast.desc.spite.mo.ar.lm)
plast.desc.spite.mo.speed.lm <- data.frame(plast.desc.spite.mo.speed.lm)
colnames(plast.desc.spite.mo.area.lm) <- c('sal.ori', 'sal.des', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.spite.mo.ar.lm) <- c('sal.ori', 'sal.des', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.spite.mo.speed.lm) <- c('sal.ori', 'sal.des', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')


plast.desc.spite.mi.area.lm <- plast.desc.spite.mi.ar.lm <- plast.desc.spite.mi.speed.lm <- vector()
for(ori in 1:5){
	dat.tmp <- dat09.PTS.Spite[dat09.PTS.Spite$Salinity_Origin == sal.ori[ori],]

	for(des in 1:3){
		nr <- nrow(dat.tmp[dat.tmp$Salinity_Destination == sal.des[des],])
		
		if(nr > 0){
			dat.tmp2 <- rbind(dat.tmp[dat.tmp$Salinity_Destination == 0.5,], dat.tmp[dat.tmp$Salinity_Destination == sal.des[des],])
			Sal. <- factor(dat.tmp2$Salinity_Destination, levels = c(0.5, sal.des[des]))

			lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.area <- summary(lmer.area)$coefficients

			
			lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.ar <- summary(lmer.ar)$coefficients

			lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID_original) + (1|replicate:ID_original), data = dat.tmp2)
			sum.speed <- summary(lmer.speed)$coefficients	

			nr <- nrow(sum.area)

			if(nr == 3){
				plast.desc.spite.mi.area.lm <- rbind(plast.desc.spite.mi.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									sum.area[3,1], sum.area[3,2], sum.area[3,5]))
				plast.desc.spite.mi.ar.lm <- rbind(plast.desc.spite.mi.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
				plast.desc.spite.mi.speed.lm <- rbind(plast.desc.spite.mi.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
			}
			if(nr == 2){
				plast.desc.spite.mi.area.lm <- rbind(plast.desc.spite.mi.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									'NA', 'NA', 'NA'))
				plast.desc.spite.mi.ar.lm <- rbind(plast.desc.spite.mi.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									'NA', 'NA', 'NA'))
				plast.desc.spite.mi.speed.lm <- rbind(plast.desc.spite.mi.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									'NA', 'NA', 'NA'))
			}
		}
	}
}
plast.desc.spite.mi.area.lm <- data.frame(plast.desc.spite.mi.area.lm)
plast.desc.spite.mi.ar.lm <- data.frame(plast.desc.spite.mi.ar.lm)
plast.desc.spite.mi.speed.lm <- data.frame(plast.desc.spite.mi.speed.lm)
colnames(plast.desc.spite.mi.area.lm) <- c('sal.ori', 'sal.des', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.spite.mi.ar.lm) <- c('sal.ori', 'sal.des', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.spite.mi.speed.lm) <- c('sal.ori', 'sal.des', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')

### -------------------------------------------------------------------------------------------------
## Phenotypic plasticity response to salinity -- Figure S8 
### -------------------------------------------------------------------------------------------------

## You need the code of Table S11-S12 to make this figure

## Make figure -- biomass
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))
plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-7000, 4300))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.spite.mo.area.lm[i,2])
	segments(i-0.5, plast.anc.spite.mo.area.lm[i,2]-plast.anc.spite.mo.area.lm[i,3], i-0.5, plast.anc.spite.mo.area.lm[i,2]+plast.anc.spite.mo.area.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.spite.mi.area.lm[i,2])
	segments(i+3, plast.anc.spite.mi.area.lm[i,2]-plast.anc.spite.mi.area.lm[i,3], i+3, plast.anc.spite.mi.area.lm[i,2]+plast.anc.spite.mi.area.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.spite.mo.area.lm[plast.desc.spite.mo.area.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.spite.mi.area.lm[plast.desc.spite.mi.area.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-7000, 4300))
	abline(h = 0)
	
	for(des in 1:3){
		x <- as.numeric(as.character(plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]))
		y <- plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, as.numeric(as.character(plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]])))
			segments(des-0.5, as.numeric(as.character(plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]))-as.numeric(as.character(plast.mo.tmp$se.plast.area[plast.mo.tmp$sal.des == sal.des[des]])),
					des-0.5, as.numeric(as.character(plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]))+as.numeric(as.character(plast.mo.tmp$se.plast.area[plast.mo.tmp$sal.des == sal.des[des]])))
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.area[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.area[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}


## Make figure -- cell shape
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))
plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-1.5, 3))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.spite.mo.ar.lm[i,2])
	segments(i-0.5, plast.anc.spite.mo.ar.lm[i,2]-plast.anc.spite.mo.ar.lm[i,3], i-0.5, plast.anc.spite.mo.ar.lm[i,2]+plast.anc.spite.mo.ar.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.spite.mi.ar.lm[i,2])
	segments(i+3, plast.anc.spite.mi.ar.lm[i,2]-plast.anc.spite.mi.ar.lm[i,3], i+3, plast.anc.spite.mi.ar.lm[i,2]+plast.anc.spite.mi.ar.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.spite.mo.ar.lm[plast.desc.spite.mo.ar.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.spite.mi.ar.lm[plast.desc.spite.mi.ar.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-1.5, 3))
	abline(h = 0)
	
	for(des in 1:3){
		x <- as.numeric(as.character(plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]))
		y <- plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, as.numeric(as.character(plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]])))
			segments(des-0.5, as.numeric(as.character(plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]))-as.numeric(as.character(plast.mo.tmp$se.plast.ar[plast.mo.tmp$sal.des == sal.des[des]])),
					des-0.5, as.numeric(as.character(plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]))+as.numeric(as.character(plast.mo.tmp$se.plast.ar[plast.mo.tmp$sal.des == sal.des[des]])))
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.ar[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.ar[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}


## Make figure - dispersal ability
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))

plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-300, 400))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.spite.mo.speed.lm[i,2])
	segments(i-0.5, plast.anc.spite.mo.speed.lm[i,2]-plast.anc.spite.mo.speed.lm[i,3], i-0.5, plast.anc.spite.mo.speed.lm[i,2]+plast.anc.spite.mo.speed.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.spite.mi.speed.lm[i,2])
	segments(i+3, plast.anc.spite.mi.speed.lm[i,2]-plast.anc.spite.mi.speed.lm[i,3], i+3, plast.anc.spite.mi.speed.lm[i,2]+plast.anc.spite.mi.speed.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.spite.mo.speed.lm[plast.desc.spite.mo.speed.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.spite.mi.speed.lm[plast.desc.spite.mi.speed.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-300, 400))
	abline(h = 0)
	
	for(des in 1:3){
		x <- as.numeric(as.character(plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]))
		y <- plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, as.numeric(as.character(plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]])))
			segments(des-0.5, as.numeric(as.character(plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]))-as.numeric(as.character(plast.mo.tmp$se.plast.speed[plast.mo.tmp$sal.des == sal.des[des]])),
					des-0.5, as.numeric(as.character(plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]))+as.numeric(as.character(plast.mo.tmp$se.plast.speed[plast.mo.tmp$sal.des == sal.des[des]])))
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.speed[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.speed[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}












