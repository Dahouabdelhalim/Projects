## R code for figures ans statistical analysis of manuscript 
## "XXX" 
## R code for figures provide raw draft for figure which is then styled in Illustrator

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
## PARAMECIUM AURELIA 
## ----------------------------------------------------------------------------------

## ----------------------------------------------------------------------------------
## Temporal trait shift -- selection -- Figure 2 Main text 
## ----------------------------------------------------------------------------------

par(mfrow = c(1,3))

m22.area <- tapply(dat22.Pau.CG$mean_area, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), mean)
sd22.area <- tapply(dat22.Pau.CG$mean_area, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), sd)
m22.area <- data.frame(cbind(RTD(m22.area), m22.area, sd22.area))
mPTS22.area <- tapply(dat22.PTS.Pau.CG$mean_area, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), mean)
sdPTS22.area <- tapply(dat22.PTS.Pau.CG$mean_area, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), sd)
mPTS22.area <- data.frame(cbind(RTD(mPTS22.area), mPTS22.area, sdPTS22.area))

m05.area <- tapply(dat05.Pau$mean_area, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), mean)
sd05.area <- tapply(dat05.Pau$mean_area, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), sd)
m05.area <- data.frame(cbind(RTD(m05.area), m05.area, sd05.area))
mPTS05.area <- tapply(dat05.PTS.Pau$mean_area, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), mean)
sdPTS05.area <- tapply(dat05.PTS.Pau$mean_area, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), sd)
mPTS05.area <- data.frame(cbind(RTD(mPTS05.area), mPTS05.area, sdPTS05.area))

colnames(m22.area) <- colnames(mPTS22.area) <- colnames(m05.area) <- colnames(mPTS05.area) <- c('ID', 'sal', 'area', 'sd.area')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(1000,6500), type = 'n', main = 'area')
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



m22.ar <- tapply(dat22.Pau.CG$mean_ar, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), mean)
sd22.ar <- tapply(dat22.Pau.CG$mean_ar, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), sd)
m22.ar <- data.frame(cbind(RTD(m22.ar), m22.ar, sd22.ar))
mPTS22.ar <- tapply(dat22.PTS.Pau.CG$mean_ar, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), mean)
sdPTS22.ar <- tapply(dat22.PTS.Pau.CG$mean_ar, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), sd)
mPTS22.ar <- data.frame(cbind(RTD(mPTS22.ar), mPTS22.ar, sdPTS22.ar))

m05.ar <- tapply(dat05.Pau$mean_ar, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), mean)
sd05.ar <- tapply(dat05.Pau$mean_ar, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), sd)
m05.ar <- data.frame(cbind(RTD(m05.ar), m05.ar, sd05.ar))
mPTS05.ar <- tapply(dat05.PTS.Pau$mean_ar, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), mean)
sdPTS05.ar <- tapply(dat05.PTS.Pau$mean_ar, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), sd)
mPTS05.ar <- data.frame(cbind(RTD(mPTS05.ar), mPTS05.ar, sdPTS05.ar))

colnames(m22.ar) <- colnames(mPTS22.ar) <- colnames(m05.ar) <- colnames(mPTS05.ar) <- c('ID', 'sal', 'ar', 'sd.ar')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(1.7,3.5), type = 'n', main = 'ar')
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


m22.speed <- tapply(dat22.Pau.CG$gross_speed, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), mean)
sd22.speed <- tapply(dat22.Pau.CG$gross_speed, paste(dat22.Pau.CG$ID, dat22.Pau.CG$Salinity), sd)
m22.speed <- data.frame(cbind(RTD(m22.speed), m22.speed, sd22.speed))
mPTS22.speed <- tapply(dat22.PTS.Pau.CG$gross_speed, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), mean)
sdPTS22.speed <- tapply(dat22.PTS.Pau.CG$gross_speed, paste(dat22.PTS.Pau.CG$ID, dat22.PTS.Pau.CG$Salinity), sd)
mPTS22.speed <- data.frame(cbind(RTD(mPTS22.speed), mPTS22.speed, sdPTS22.speed))

m05.speed <- tapply(dat05.Pau$gross_speed, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), mean)
sd05.speed <- tapply(dat05.Pau$gross_speed, paste(dat05.Pau$Sample_ID, dat05.Pau$Salinity), sd)
m05.speed <- data.frame(cbind(RTD(m05.speed), m05.speed, sd05.speed))
mPTS05.speed <- tapply(dat05.PTS.Pau$gross_speed, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), mean)
sdPTS05.speed <- tapply(dat05.PTS.Pau$gross_speed, paste(dat05.PTS.Pau$Sample_ID, dat05.PTS.Pau$Salinity), sd)
mPTS05.speed <- data.frame(cbind(RTD(mPTS05.speed), mPTS05.speed, sdPTS05.speed))

colnames(m22.speed) <- colnames(mPTS22.speed) <- colnames(m05.speed) <- colnames(mPTS05.speed) <- c('ID', 'sal', 'speed', 'sd.speed')

sal.lvl <- c(0, 0.5, 1, 2, 4)
plot(-100, -100, xlim = c(0,2), ylim = c(200,1000), type = 'n', main = 'speed')
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
## Temporal trait shift to salinity -- selection -- Supplementary Table S1
### -------------------------------------------------------------------------------------------------

tt.sal <- c(dat22.Pau.CG$Salinity, dat05.Pau$Salinity, dat22.PTS.Pau.CG$Salinity, dat05.PTS.Pau$Salinity)
tt.time <- c(dat22.Pau.CG$Time, dat05.Pau$Time, dat22.PTS.Pau.CG$Time, dat05.PTS.Pau$Time)
tt.id <- c(dat22.Pau.CG$ID, dat05.Pau$Sample_ID, dat22.PTS.Pau.CG$ID, dat05.PTS.Pau$Sample_ID)
tt.com <- factor(c(dat22.Pau.CG$Community, dat05.Pau$Community, dat22.PTS.Pau.CG$Community, dat05.PTS.Pau$Community), levels=c('mono','mixed'))
tt.bf <- c(dat22.Pau.CG$Bio.Fraction.other, dat05.Pau$Bio.Fraction.other, dat22.PTS.Pau.CG$Bio.Fraction.other, dat05.PTS.Pau$Bio.Fraction.other)
tt.dens <- c(dat22.Pau.CG$Density, dat05.Pau$Density, dat22.PTS.Pau.CG$Density, dat05.PTS.Pau$Density)
tt.area <- c(dat22.Pau.CG$mean_area, dat05.Pau$mean_area, dat22.PTS.Pau.CG$mean_area, dat05.PTS.Pau$mean_area)
tt.ar <- c(dat22.Pau.CG$mean_ar, dat05.Pau$mean_ar, dat22.PTS.Pau.CG$mean_ar, dat05.PTS.Pau$mean_ar)
tt.speed <- c(dat22.Pau.CG$gross_speed, dat05.Pau$gross_speed, dat22.PTS.Pau.CG$gross_speed, dat05.PTS.Pau$gross_speed)

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
## Phenotypic trait response to salinity and competition -- common garden -- Figure 3 Main text
### -------------------------------------------------------------------------------------------------

dat09.combined.Pau <- rbind(dat09.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Pau$Salinity_Origin==0 & dat09.combined.Pau$Salinity_Destination==0)
w0.dat09.combined.Pau <- dat09.combined.Pau[-pos.0,]
w0.dat09.combined.Pau$Community <- factor(w0.dat09.combined.Pau$Community, levels=c('mono','mixed'))

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)

par(mfrow=c(1,3))
c.area <- summary(fit.area)$coefficients

plot(0:5, rep(-10000,6), type='n', ylim=c(-300,600))
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

plot(0:5, rep(-10000,6), type='n', ylim=c(-0.6,0.2))
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

plot(0:5, rep(-10000,6), type='n', ylim=c(-300,200))
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
## Phenotypic trait response to salinity and competition -- common garden -- Supplementary Table S3
### -------------------------------------------------------------------------------------------------

dat09.combined.Pau <- rbind(dat09.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Pau$Salinity_Origin==0 & dat09.combined.Pau$Salinity_Destination==0)
w0.dat09.combined.Pau <- dat09.combined.Pau[-pos.0,]
w0.dat09.combined.Pau$Community <- factor(w0.dat09.combined.Pau$Community, levels=c('mono','mixed'))

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
summary(fit.area)
r.squaredGLMM(fit.area)

plot(fit.area)
plot(cooks.distance(fit.area))
qqnorm(resid(fit.area))
qqline(resid(fit.area))
hist(resid(fit.area), breaks=50)

fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
summary(fit.ar)
r.squaredGLMM(fit.ar)

plot(fit.ar)
plot(cooks.distance(fit.ar))
qqnorm(resid(fit.ar))
qqline(resid(fit.ar))
hist(resid(fit.ar), breaks=50)

fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
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

## Area
## Preparatory calculations
cat. <- paste(dat09.combined.Pau$ID_original, dat09.combined.Pau$Salinity_Origin, dat09.combined.Pau$Salinity_Destination, dat09.combined.Pau$replicate,
		dat09.combined.Pau$Community)
des.area <- tapply(dat09.combined.Pau$mean_area, cat., mean)
des.area.sd <- tapply(dat09.combined.Pau$mean_area, cat., sd)
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

## Plot for monocultures
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(1700, 4200), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - area') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.area[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.area[id.tmp$sal.des == 0]-id.tmp$sd.area[id.tmp$sal.des == 0],-0.2, id.tmp$m.area[id.tmp$sal.des == 0]+id.tmp$sd.area[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]-id.tmp$sd.area[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]+id.tmp$sd.area[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]-id.tmp$sd.area[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]+id.tmp$sd.area[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]-id.tmp$sd.area[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]+id.tmp$sd.area[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2])

		x <- id.tmp$m.area[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]-id.tmp$sd.area[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]+id.tmp$sd.area[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4])
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

## Plot for competition
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(1700, 4200), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - area') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.area[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.area[id.tmp$sal.des == 0]-id.tmp$sd.area[id.tmp$sal.des == 0],-0.2, id.tmp$m.area[id.tmp$sal.des == 0]+id.tmp$sd.area[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]-id.tmp$sd.area[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5]+id.tmp$sd.area[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]-id.tmp$sd.area[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1]+id.tmp$sd.area[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]-id.tmp$sd.area[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2]+id.tmp$sd.area[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 2])

		x <- id.tmp$m.area[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]-id.tmp$sd.area[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 4]+id.tmp$sd.area[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.area[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.area[id.tmp$sal.des == 4])
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.area[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.area[tmp.$sal.des == 0]-tmp.$sd.area[tmp.$sal.des == 0],-0.2, tmp.$m.area[tmp.$sal.des == 0]+tmp.$sd.area[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]-tmp.$sd.area[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5]+tmp.$sd.area[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]-tmp.$sd.area[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1]+tmp.$sd.area[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]-tmp.$sd.area[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2]+tmp.$sd.area[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.area[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]-tmp.$sd.area[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4]+tmp.$sd.area[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.area[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.area[tmp.$sal.des == 4], lwd = 2)
	}
}

## Cell shape 
cat. <- paste(dat09.combined.Pau$ID_original, dat09.combined.Pau$Salinity_Origin, dat09.combined.Pau$Salinity_Destination, dat09.combined.Pau$replicate,
		dat09.combined.Pau$Community)
des.ar <- tapply(dat09.combined.Pau$mean_ar, cat., mean)
des.ar.sd <- tapply(dat09.combined.Pau$mean_ar, cat., sd)
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

## reaction norms monocultures
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(2.2, 4.3), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - ar') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]-id.tmp$sd.ar[id.tmp$sal.des == 0],-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]+id.tmp$sd.ar[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]-id.tmp$sd.ar[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]+id.tmp$sd.ar[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]-id.tmp$sd.ar[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]+id.tmp$sd.ar[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]-id.tmp$sd.ar[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]+id.tmp$sd.ar[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2])

		x <- id.tmp$m.ar[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 4]-id.tmp$sd.ar[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 4]+id.tmp$sd.ar[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4])
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mo[mm.des.mo$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.ar[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.ar[tmp.$sal.des == 0]-tmp.$sd.ar[tmp.$sal.des == 0],-0.2, tmp.$m.ar[tmp.$sal.des == 0]+tmp.$sd.ar[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.ar[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], lwd = 2)
	}
}

## reaction norm competition
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(2.2, 4.3), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - ar') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]-id.tmp$sd.ar[id.tmp$sal.des == 0],-0.2, id.tmp$m.ar[id.tmp$sal.des == 0]+id.tmp$sd.ar[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]-id.tmp$sd.ar[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5]+id.tmp$sd.ar[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]-id.tmp$sd.ar[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1]+id.tmp$sd.ar[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]-id.tmp$sd.ar[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2]+id.tmp$sd.ar[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 2])

		x <- id.tmp$m.ar[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 4]-id.tmp$sd.ar[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 4]+id.tmp$sd.ar[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.ar[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.ar[id.tmp$sal.des == 4])
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.ar[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.ar[tmp.$sal.des == 0]-tmp.$sd.ar[tmp.$sal.des == 0],-0.2, tmp.$m.ar[tmp.$sal.des == 0]+tmp.$sd.ar[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5]-tmp.$sd.ar[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1]-tmp.$sd.ar[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2]-tmp.$sd.ar[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.ar[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4]-tmp.$sd.ar[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.ar[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.ar[tmp.$sal.des == 4], lwd = 2)
	}
}

## Dispersal ability
## Preparatory calculations
cat. <- paste(dat09.combined.Pau$ID_original, dat09.combined.Pau$Salinity_Origin, dat09.combined.Pau$Salinity_Destination, dat09.combined.Pau$replicate,
		dat09.combined.Pau$Community)
des.speed <- tapply(dat09.combined.Pau$gross_speed, cat., mean)
des.speed.sd <- tapply(dat09.combined.Pau$gross_speed, cat., sd)
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

## plot monocultures
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(450, 1900), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - speed') 
for(i in 1:5){
	tmp. <- mm.des.mo.id[mm.des.mo.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]-id.tmp$sd.speed[id.tmp$sal.des == 0],-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]+id.tmp$sd.speed[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]-id.tmp$sd.speed[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]+id.tmp$sd.speed[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]-id.tmp$sd.speed[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]+id.tmp$sd.speed[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]-id.tmp$sd.speed[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]+id.tmp$sd.speed[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2])

		x <- id.tmp$m.speed[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 4]-id.tmp$sd.speed[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 4]+id.tmp$sd.speed[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4])
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mo[mm.des.mo$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.speed[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.speed[tmp.$sal.des == 0]-tmp.$sd.speed[tmp.$sal.des == 0],-0.2, tmp.$m.speed[tmp.$sal.des == 0]+tmp.$sd.speed[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.speed[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], lwd = 2)
	}
}

## plot competition
sal. <- c(0,0.5,1,2,4)
plot(-10000, -10000, xlim = c(-0.5,4), ylim = c(450, 1900), type = 'n', xlab = 'evolved salinity', ylab = 'plasticity', main = 'Pau - speed') 
for(i in 1:5){
	tmp. <- mm.des.mi.id[mm.des.mi.id$sal.ori == sal.[i],]
	id. <- unique(tmp.$ID_ori)

	for(id in 1:length(id.)){
		id.tmp <- tmp.[tmp.$ID_ori == id.[id],]

		if(i == 1){
			points(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0], pch = 21, cex = 1)
			segments(-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]-id.tmp$sd.speed[id.tmp$sal.des == 0],-0.2, id.tmp$m.speed[id.tmp$sal.des == 0]+id.tmp$sd.speed[id.tmp$sal.des == 0])
		}
		points(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5], pch = 21, cex = 1)
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]-id.tmp$sd.speed[id.tmp$sal.des == 0.5], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5]+id.tmp$sd.speed[id.tmp$sal.des == 0.5])
		points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1], pch = 21, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]-id.tmp$sd.speed[id.tmp$sal.des == 1], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1]+id.tmp$sd.speed[id.tmp$sal.des == 1])
		points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2], pch = 22, cex = 1)
		segments(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]-id.tmp$sd.speed[id.tmp$sal.des == 2], 0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2]+id.tmp$sd.speed[id.tmp$sal.des == 2])
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 1])
		segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 2])

		x <- id.tmp$m.speed[id.tmp$sal.des == 4]
		if(length(x) > 0){
			points(0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4], pch = 24, cex = 1)
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 4]-id.tmp$sd.speed[id.tmp$sal.des == 4], 0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 4]+id.tmp$sd.speed[id.tmp$sal.des == 4])
			segments(0.5*(i-1), id.tmp$m.speed[id.tmp$sal.des == 0.5],0.5*(i-1)+0.4, id.tmp$m.speed[id.tmp$sal.des == 4])
		}
	}
}
for(i in 1:5){
	tmp. <- mm.des.mi[mm.des.mi$sal.ori == sal.[i],]

	if(i == 1){
		points(-0.2, tmp.$m.speed[tmp.$sal.des == 0], pch = 21, cex = 2)
		segments(-0.2, tmp.$m.speed[tmp.$sal.des == 0]-tmp.$sd.speed[tmp.$sal.des == 0],-0.2, tmp.$m.speed[tmp.$sal.des == 0]+tmp.$sd.speed[tmp.$sal.des == 0])
	}
	points(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5], pch = 21, cex = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5], 0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5]-tmp.$sd.speed[tmp.$sal.des == 0.5])
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], pch = 21, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1]-tmp.$sd.speed[tmp.$sal.des == 1])
	points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], pch = 22, cex = 2)
	segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2]-tmp.$sd.speed[tmp.$sal.des == 2])
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 1], lwd = 2)
	segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 2], lwd = 2)

	x <- tmp.$m.speed[tmp.$sal.des == 4]
	if(length(x) > 0){
		points(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], pch = 24, cex = 2)
		segments(0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4], 0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4]-tmp.$sd.speed[tmp.$sal.des == 4])
		segments(0.5*(i-1), tmp.$m.speed[tmp.$sal.des == 0.5],0.5*(i-1)+0.4, tmp.$m.speed[tmp.$sal.des == 4], lwd = 2)
	}
}


### -------------------------------------------------------------------------------------------------
## Bootstrap robustness analysis -- Figure S3
### -------------------------------------------------------------------------------------------------

dat09.combined.Pau <- rbind(dat09.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

pos.0 <- which(dat09.combined.Pau$Salinity_Origin==0 & dat09.combined.Pau$Salinity_Destination==0)
w0.dat09.combined.Pau <- dat09.combined.Pau[-pos.0,]
w0.dat09.combined.Pau$Community <- factor(w0.dat09.combined.Pau$Community, levels=c('mono','mixed'))

fit.area <- lmer(mean_area ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
summary(fit.area)
r.squaredGLMM(fit.area)

fit.ar <- lmer(mean_ar ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
summary(fit.ar)

fit.speed <- lmer(gross_speed ~ Salinity_Origin*Salinity_Destination*Community + Bio.Fraction.other + Density + (1|Salinity_Origin:Salinity_Destination:replicate) + (1|ID_original) , data=w0.dat09.combined.Pau)
summary(fit.speed)

## make the data on which we will do the bootstrapping
boot.data <- w0.dat09.combined.Pau

## create a variable that contains the bootstrap categories
tt0 <- paste(w0.dat09.combined.Pau$Community,w0.dat09.combined.Pau$Salinity_Origin,w0.dat09.combined.Pau$Salinity_Destination,w0.dat09.combined.Pau$ID_original,w0.dat09.combined.Pau$replicate)
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

## Making the plots
## choose kk as a value or use pdf to create a pdf 

## area
pathway <- #"Give in Pathway"
pdf(paste(pathway, 'EffectSizeWithCI_Coef_area.pdf', sep=""))

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


## cell shape
pathway <- #"Give in Pathway"
pdf(paste(pathway, 'EffectSize_Coef_ar.pdf', sep=""))

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
pdf(paste(pathway, 'EffectSize_Coef_speed.pdf', sep=""))

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
## Reaction norms analysis on selected populations -- Table S5
### -------------------------------------------------------------------------------------------------

#remove ID 120
# pos.120 <- which(dat09.PTS.Pau$ID_original == '120')
# dat09.PTS.Pau <- dat09.PTS.Pau[-pos.120,]

## 0. Try to fit the reaction norm approach to each selection population
sal.lvl <- c(0, 0.5, 1, 2, 4)

### Anc 0.5-1 vs All Des (0.5-1)
dat22.Pau.1 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 1,])

## area
lmer.area.mo.X.1 <- lmer.area.mi.X.1 <- aov.area.mo.X.1 <- aov.area.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Pau.X.1 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 1 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.1 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 1 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.area.mo <- c(dat22.Pau.1$mean_area, dat09.Pau.X.1$mean_area)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.1$mean_area)), dat09.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.1$Salinity, dat09.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Pau.1$Density, dat09.Pau.X.1$Density)
	tt.id.mo <- c(dat22.Pau.1$ID, dat09.Pau.X.1$ID_new)

	tt.area.mi <- c(dat22.Pau.1$mean_area, dat09.PTS.Pau.X.1$mean_area)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.1$mean_area)), dat09.PTS.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.1$Salinity, dat09.PTS.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Pau.1$Density, dat09.PTS.Pau.X.1$Density)
	tt.id.mi <- c(dat22.Pau.1$ID, dat09.PTS.Pau.X.1$ID_new)

	lmer.mo <- lmer(tt.area.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.area.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.area.mo.X.1 <- rbind(lmer.area.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.area.mo.X.1 <- rbind(aov.area.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.area.mi.X.1 <- rbind(lmer.area.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.area.mi.X.1 <- rbind(aov.area.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.area.mo.X.1 <- data.frame(lmer.area.mo.X.1) ## gives output for Table S5
lmer.area.mi.X.1 <- data.frame(lmer.area.mi.X.1) ## gives output for Table S6
aov.area.mo.X.1 <- data.frame(aov.area.mo.X.1)
aov.area.mi.X.1 <- data.frame(aov.area.mi.X.1)
colnames(lmer.area.mo.X.1) <- colnames(lmer.area.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.area.mo.X.1) <- colnames(aov.area.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')

### Anc 0.5-2 vs All Des (0.5-2)
dat22.Pau.2 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 2,])

lmer.area.mo.X.2 <- lmer.area.mi.X.2 <- aov.area.mo.X.2 <- aov.area.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Pau.X.2 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 2 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.2 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 2 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.area.mo <- c(dat22.Pau.2$mean_area, dat09.Pau.X.2$mean_area)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.2$mean_area)), dat09.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.2$Salinity, dat09.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Pau.2$Density, dat09.Pau.X.2$Density)
	tt.id.mo <- c(dat22.Pau.2$ID, dat09.Pau.X.2$ID_new)

	tt.area.mi <- c(dat22.Pau.2$mean_area, dat09.PTS.Pau.X.2$mean_area)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.2$mean_area)), dat09.PTS.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.2$Salinity, dat09.PTS.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Pau.2$Density, dat09.PTS.Pau.X.2$Density)
	tt.id.mi <- c(dat22.Pau.2$ID, dat09.PTS.Pau.X.2$ID_new)

	lmer.mo <- lmer(tt.area.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.area.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.area.mo.X.2 <- rbind(lmer.area.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.area.mo.X.2 <- rbind(aov.area.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.area.mi.X.2 <- rbind(lmer.area.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], 	sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
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

## ar
## 0. Try to fit the reaction norm approach to each selection population
sal.lvl <- c(0, 0.5, 1, 2, 4)

### Anc 0.5-1 vs All Des (0.5-1)
dat22.Pau.1 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 1,])

lmer.ar.mo.X.1 <- lmer.ar.mi.X.1 <- aov.ar.mo.X.1 <- aov.ar.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Pau.X.1 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 1 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.1 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 1 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.ar.mo <- c(dat22.Pau.1$mean_ar, dat09.Pau.X.1$mean_ar)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.1$mean_ar)), dat09.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.1$Salinity, dat09.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Pau.1$Density, dat09.Pau.X.1$Density)
	tt.id.mo <- c(dat22.Pau.1$ID, dat09.Pau.X.1$ID_new)

	tt.ar.mi <- c(dat22.Pau.1$mean_ar, dat09.PTS.Pau.X.1$mean_ar)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.1$mean_ar)), dat09.PTS.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.1$Salinity, dat09.PTS.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Pau.1$Density, dat09.PTS.Pau.X.1$Density)
	tt.id.mi <- c(dat22.Pau.1$ID, dat09.PTS.Pau.X.1$ID_new)

	lmer.mo <- lmer(tt.ar.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.ar.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.ar.mo.X.1 <- rbind(lmer.ar.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.ar.mo.X.1 <- rbind(aov.ar.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.ar.mi.X.1 <- rbind(lmer.ar.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.ar.mi.X.1 <- rbind(aov.ar.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.ar.mo.X.1 <- data.frame(lmer.ar.mo.X.1) ## gives output for Table S5
lmer.ar.mi.X.1 <- data.frame(lmer.ar.mi.X.1) ## gives output for Table S6
aov.ar.mo.X.1 <- data.frame(aov.ar.mo.X.1)
aov.ar.mi.X.1 <- data.frame(aov.ar.mi.X.1)
colnames(lmer.ar.mo.X.1) <- colnames(lmer.ar.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.ar.mo.X.1) <- colnames(aov.ar.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')

### Anc 0.5-2 vs All Des (0.5-2)
dat22.Pau.2 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 2,])

lmer.ar.mo.X.2 <- lmer.ar.mi.X.2 <- aov.ar.mo.X.2 <- aov.ar.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Pau.X.2 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 2 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.2 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 2 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.ar.mo <- c(dat22.Pau.2$mean_ar, dat09.Pau.X.2$mean_ar)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.2$mean_ar)), dat09.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.2$Salinity, dat09.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Pau.2$Density, dat09.Pau.X.2$Density)
	tt.id.mo <- c(dat22.Pau.2$ID, dat09.Pau.X.2$ID_new)

	tt.ar.mi <- c(dat22.Pau.2$mean_ar, dat09.PTS.Pau.X.2$mean_ar)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.2$mean_ar)), dat09.PTS.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.2$Salinity, dat09.PTS.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Pau.2$Density, dat09.PTS.Pau.X.2$Density)
	tt.id.mi <- c(dat22.Pau.2$ID, dat09.PTS.Pau.X.2$ID_new)

	lmer.mo <- lmer(tt.ar.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.ar.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.ar.mo.X.2 <- rbind(lmer.ar.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.ar.mo.X.2 <- rbind(aov.ar.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.ar.mi.X.2 <- rbind(lmer.ar.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], 	sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
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
dat22.Pau.1 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 1,])

lmer.speed.mo.X.1 <- lmer.speed.mi.X.1 <- aov.speed.mo.X.1 <- aov.speed.mi.X.1 <- vector()
for(sal in 1:5){
	dat09.Pau.X.1 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 1 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.1 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 1 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.speed.mo <- c(dat22.Pau.1$gross_speed, dat09.Pau.X.1$gross_speed)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.1$gross_speed)), dat09.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.1$Salinity, dat09.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '1'))
	tt.dens.mo <- c(dat22.Pau.1$Density, dat09.Pau.X.1$Density)
	tt.id.mo <- c(dat22.Pau.1$ID, dat09.Pau.X.1$ID_new)

	tt.speed.mi <- c(dat22.Pau.1$gross_speed, dat09.PTS.Pau.X.1$gross_speed)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.1$gross_speed)), dat09.PTS.Pau.X.1$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.1$Salinity, dat09.PTS.Pau.X.1$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '1'))
	tt.dens.mi <- c(dat22.Pau.1$Density, dat09.PTS.Pau.X.1$Density)
	tt.id.mi <- c(dat22.Pau.1$ID, dat09.PTS.Pau.X.1$ID_new)

	lmer.mo <- lmer(tt.speed.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.speed.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.speed.mo.X.1 <- rbind(lmer.speed.mo.X.1, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.speed.mo.X.1 <- rbind(aov.speed.mo.X.1, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.speed.mi.X.1 <- rbind(lmer.speed.mi.X.1, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
					sum.mi$coefficients[5,2], sum.mi$coefficients[5,5], sum.mi$coefficients[4,1], sum.mi$coefficients[4,2],
					sum.mi$coefficients[4,5]))

	aov.speed.mi.X.1 <- rbind(aov.speed.mi.X.1, c(sal.lvl[sal], aov.mi$"Pr(>F)"[c(1,2,4,3)]))
}
lmer.speed.mo.X.1 <- data.frame(lmer.speed.mo.X.1) ## gives output for Table S5
lmer.speed.mi.X.1 <- data.frame(lmer.speed.mi.X.1) ## gives output for Table S6
aov.speed.mo.X.1 <- data.frame(aov.speed.mo.X.1)
aov.speed.mi.X.1 <- data.frame(aov.speed.mi.X.1)
colnames(lmer.speed.mo.X.1) <- colnames(lmer.speed.mi.X.1) <- c('sal.ori', 'plast', 'se.plast', 'p.plast', 'evo', 'se.evo', 'p.evo', 'evoplast', 'se.evoplast', 'p.evoplast', 'dens', 'se.dens', 'p.dens')
colnames(aov.speed.mo.X.1) <- colnames(aov.speed.mi.X.1) <- c('sal.ori', 'p.plast', 'p.evo', 'p.evoplast', 'p.dens')

### Anc 0.5-2 vs All Des (0.5-2)
dat22.Pau.2 <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == 2,])

lmer.speed.mo.X.2 <- lmer.speed.mi.X.2 <- aov.speed.mo.X.2 <- aov.speed.mi.X.2 <- vector()
for(sal in 1:5){
	dat09.Pau.X.2 <- rbind(dat09.Pau[dat09.Pau$Salinity_Destination == 0.5 & dat09.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.Pau[dat09.Pau$Salinity_Destination == 2 & dat09.Pau$Salinity_Origin == sal.lvl[sal],])
	dat09.PTS.Pau.X.2 <- rbind(dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 0.5 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],], 
				dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Destination == 2 & dat09.PTS.Pau$Salinity_Origin == sal.lvl[sal],])
	
	tt.speed.mo <- c(dat22.Pau.2$gross_speed, dat09.Pau.X.2$gross_speed)
	tt.sal.gen.mo <- c(rep('A', length(dat22.Pau.2$gross_speed)), dat09.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mo <- factor(tt.sal.gen.mo, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mo <- c(dat22.Pau.2$Salinity, dat09.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mo <- factor(tt.sal.plast.mo, levels = c('0.5', '2'))
	tt.dens.mo <- c(dat22.Pau.2$Density, dat09.Pau.X.2$Density)
	tt.id.mo <- c(dat22.Pau.2$ID, dat09.Pau.X.2$ID_new)

	tt.speed.mi <- c(dat22.Pau.2$gross_speed, dat09.PTS.Pau.X.2$gross_speed)
	tt.sal.gen.mi <- c(rep('A', length(dat22.Pau.2$gross_speed)), dat09.PTS.Pau.X.2$Salinity_Origin)
	tt.sal.gen.mi <- factor(tt.sal.gen.mi, levels = c('A', sal.lvl[sal]))
	tt.sal.plast.mi <- c(dat22.Pau.2$Salinity, dat09.PTS.Pau.X.2$Salinity_Destination)
	tt.sal.plast.mi <- factor(tt.sal.plast.mi, levels = c('0.5', '2'))
	tt.dens.mi <- c(dat22.Pau.2$Density, dat09.PTS.Pau.X.2$Density)
	tt.id.mi <- c(dat22.Pau.2$ID, dat09.PTS.Pau.X.2$ID_new)

	lmer.mo <- lmer(tt.speed.mo ~ tt.sal.plast.mo*tt.sal.gen.mo + tt.dens.mo + (1|tt.id.mo))
	sum.mo <- summary(lmer.mo)
	aov.mo <- anova(lmer.mo)

	lmer.mi <- lmer(tt.speed.mi ~ tt.sal.plast.mi*tt.sal.gen.mi + tt.dens.mi + (1|tt.id.mi))
	sum.mi <- summary(lmer.mi)
	aov.mi <- anova(lmer.mi)

	lmer.speed.mo.X.2 <- rbind(lmer.speed.mo.X.2, c(sal.lvl[sal], sum.mo$coefficients[2,1], sum.mo$coefficients[2,2], sum.mo$coefficients[2,5], 
					sum.mo$coefficients[3,1], sum.mo$coefficients[3,2], 	sum.mo$coefficients[3,5], sum.mo$coefficients[5,1], 
					sum.mo$coefficients[5,2], sum.mo$coefficients[5,5], sum.mo$coefficients[4,1], sum.mo$coefficients[4,2],
					sum.mo$coefficients[4,5]))

	aov.speed.mo.X.2 <- rbind(aov.speed.mo.X.2, c(sal.lvl[sal], aov.mo$"Pr(>F)"[c(1,2,4,3)]))
	
	lmer.speed.mi.X.2 <- rbind(lmer.speed.mi.X.2, c(sal.lvl[sal], sum.mi$coefficients[2,1], sum.mi$coefficients[2,2], sum.mi$coefficients[2,5], 
					sum.mi$coefficients[3,1], sum.mi$coefficients[3,2], 	sum.mi$coefficients[3,5], sum.mi$coefficients[5,1], 
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
## Reaction norms analysis on selected populations -- Figure S5 --> mean trait components are use for Figure 4 Main text
### -------------------------------------------------------------------------------------------------

## You need to run the previous code of Table S5 in order to make this figure

## Plot the estimates!

## area
dev.new(height = 14, width = 18)
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

## ar
dev.new(height = 14, width = 18)
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
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mo.X.1$plast[sal])
	segments(0.15, lmer.ar.mo.X.1$plast[sal]-lmer.ar.mo.X.1$se.plast[sal], 0.15, lmer.ar.mo.X.1$plast[sal]+lmer.ar.mo.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mo.X.1$evo[sal])
	segments(0.5, lmer.ar.mo.X.1$evo[sal]-lmer.ar.mo.X.1$se.evo[sal], 0.5, lmer.ar.mo.X.1$evo[sal]+lmer.ar.mo.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mo.X.1$evoplast[sal])
	segments(0.85, lmer.ar.mo.X.1$evoplast[sal]-lmer.ar.mo.X.1$se.evoplast[sal], 0.85, lmer.ar.mo.X.1$evoplast[sal]+lmer.ar.mo.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mono')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mo.X.2$plast[sal])
	segments(0.15, lmer.ar.mo.X.2$plast[sal]-lmer.ar.mo.X.2$se.plast[sal], 0.15, lmer.ar.mo.X.2$plast[sal]+lmer.ar.mo.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mo.X.2$evo[sal])
	segments(0.5, lmer.ar.mo.X.2$evo[sal]-lmer.ar.mo.X.2$se.evo[sal], 0.5, lmer.ar.mo.X.2$evo[sal]+lmer.ar.mo.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mo.X.2$evoplast[sal])
	segments(0.85, lmer.ar.mo.X.2$evoplast[sal]-lmer.ar.mo.X.2$se.evoplast[sal], 0.85, lmer.ar.mo.X.2$evoplast[sal]+lmer.ar.mo.X.2$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mi.X.1$plast[sal])
	segments(0.15, lmer.ar.mi.X.1$plast[sal]-lmer.ar.mi.X.1$se.plast[sal], 0.15, lmer.ar.mi.X.1$plast[sal]+lmer.ar.mi.X.1$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mi.X.1$evo[sal])
	segments(0.5, lmer.ar.mi.X.1$evo[sal]-lmer.ar.mi.X.1$se.evo[sal], 0.5, lmer.ar.mi.X.1$evo[sal]+lmer.ar.mi.X.1$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mi.X.1$evoplast[sal])
	segments(0.85, lmer.ar.mi.X.1$evoplast[sal]-lmer.ar.mi.X.1$se.evoplast[sal], 0.85, lmer.ar.mi.X.1$evoplast[sal]+lmer.ar.mi.X.1$se.evoplast[sal])
}
for(sal in 1:5){
	plot(-100, -100, type = 'n', xlim = c(0,1), ylim = c(ymin, ymax), main = 'comparison of ancestral mixed')
	abline(h=0)
	rect(0.15-0.12, 0, 0.15+0.12, lmer.ar.mi.X.2$plast[sal])
	segments(0.15, lmer.ar.mi.X.2$plast[sal]-lmer.ar.mi.X.2$se.plast[sal], 0.15, lmer.ar.mi.X.2$plast[sal]+lmer.ar.mi.X.2$se.plast[sal])
	rect(0.5-0.12, 0, 0.5+0.12, lmer.ar.mi.X.2$evo[sal])
	segments(0.5, lmer.ar.mi.X.2$evo[sal]-lmer.ar.mi.X.2$se.evo[sal], 0.5, lmer.ar.mi.X.2$evo[sal]+lmer.ar.mi.X.2$se.evo[sal])
	rect(0.85-0.12, 0, 0.85+0.12, lmer.ar.mi.X.2$evoplast[sal])
	segments(0.85, lmer.ar.mi.X.2$evoplast[sal]-lmer.ar.mi.X.2$se.evoplast[sal], 0.85, lmer.ar.mi.X.2$evoplast[sal]+lmer.ar.mi.X.2$se.evoplast[sal])
}

## speed
dev.new(height = 14, width = 18)
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
## Phenotypic plasticity response to salinity -- Table S9-S10
### -------------------------------------------------------------------------------------------------

#remove ID 120
# pos.120 <- which(dat09.PTS.Pau$ID_original == '120')
# dat09.PTS.Pau <- dat09.PTS.Pau[-pos.120,]

## Calculating plasticity responses for the ancestral and descendant population 
## Using regression with density

## Monocultures ancestral
sal.lvl <- c(1, 2)
plast.anc.pau.mo.area.lm <- plast.anc.pau.mo.ar.lm <- plast.anc.pau.mo.speed.lm <- vector()
for(sal in 1:2){
	dat.tmp <- rbind(dat22.Pau.CG[dat22.Pau.CG$Salinity == 0.5,], dat22.Pau.CG[dat22.Pau.CG$Salinity == sal.lvl[sal],])
	Sal. <- factor(dat.tmp$Salinity, levels = c(0.5, sal.lvl[sal]))

	lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.area <- summary(lmer.area)$coefficients

	lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.ar <- summary(lmer.ar)$coefficients

	lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.speed <- summary(lmer.speed)$coefficients

	plast.anc.pau.mo.area.lm <- rbind(plast.anc.pau.mo.area.lm, cbind(sal.lvl[sal], sum.area[2,1], sum.area[2,2], sum.area[2,5], 
							sum.area[3,1], sum.area[3,2], sum.area[3,5]))
	plast.anc.pau.mo.ar.lm <- rbind(plast.anc.pau.mo.ar.lm, cbind(sal.lvl[sal], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5], 
							sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
	plast.anc.pau.mo.speed.lm <- rbind(plast.anc.pau.mo.speed.lm, cbind(sal.lvl[sal], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5], 
							sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
}
plast.anc.pau.mo.area.lm <- data.frame(plast.anc.pau.mo.area.lm) ## contains the output of Table S9 and S10
plast.anc.pau.mo.ar.lm <- data.frame(plast.anc.pau.mo.ar.lm)
plast.anc.pau.mo.speed.lm <- data.frame(plast.anc.pau.mo.speed.lm)
colnames(plast.anc.pau.mo.area.lm) <- c('salinity', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.pau.mo.ar.lm) <- c('salinity', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.pau.mo.speed.lm) <- c('salinity', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')

## competition ancestral
sal.lvl <- c(1, 2, 4)
plast.anc.pau.mi.area.lm <- plast.anc.pau.mi.ar.lm <- plast.anc.pau.mi.speed.lm <- vector()
for(sal in 1:3){
	dat.tmp <- rbind(dat22.PTS.Pau.CG[dat22.PTS.Pau.CG$Salinity == 0.5,], dat22.PTS.Pau.CG[dat22.PTS.Pau.CG$Salinity == sal.lvl[sal],])
	Sal. <- factor(dat.tmp$Salinity, levels = c(0.5, sal.lvl[sal]))

	lmer.area <- lmer(mean_area ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.area <- summary(lmer.area)$coefficients

	lmer.ar <- lmer(mean_ar ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.ar <- summary(lmer.ar)$coefficients

	lmer.speed <- lmer(gross_speed ~ Sal. + Density + (1|ID), data = dat.tmp)
	sum.speed <- summary(lmer.speed)$coefficients

	plast.anc.pau.mi.area.lm <- rbind(plast.anc.pau.mi.area.lm, cbind(sal.lvl[sal], sum.area[2,1], sum.area[2,2], sum.area[2,5], 
							sum.area[3,1], sum.area[3,2], sum.area[3,5]))
	plast.anc.pau.mi.ar.lm <- rbind(plast.anc.pau.mi.ar.lm, cbind(sal.lvl[sal], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5], 
							sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
	plast.anc.pau.mi.speed.lm <- rbind(plast.anc.pau.mi.speed.lm, cbind(sal.lvl[sal], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5], 
							sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
}
plast.anc.pau.mi.area.lm <- data.frame(plast.anc.pau.mi.area.lm)
plast.anc.pau.mi.ar.lm <- data.frame(plast.anc.pau.mi.ar.lm)
plast.anc.pau.mi.speed.lm <- data.frame(plast.anc.pau.mi.speed.lm)
colnames(plast.anc.pau.mi.area.lm) <- c('salinity', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.pau.mi.ar.lm) <- c('salinity', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.anc.pau.mi.speed.lm) <- c('salinity', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')


## Descendant
dat09.combined.Pau <- rbind(dat09.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

sal.ori <- c(0, 0.5, 1, 2, 4)
sal.des <- c(1, 2, 4)

## Monocultures descendants
plast.desc.pau.mo.area.lm <- plast.desc.pau.mo.ar.lm <- plast.desc.pau.mo.speed.lm <- vector()
for(ori in 1:5){
	dat.tmp <- dat09.Pau[dat09.Pau$Salinity_Origin == sal.ori[ori],]

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

			plast.desc.pau.mo.area.lm <- rbind(plast.desc.pau.mo.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									sum.area[3,1], sum.area[3,2], sum.area[3,5]))
			plast.desc.pau.mo.ar.lm <- rbind(plast.desc.pau.mo.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
			plast.desc.pau.mo.speed.lm <- rbind(plast.desc.pau.mo.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
		}
	}
}
plast.desc.pau.mo.area.lm <- data.frame(plast.desc.pau.mo.area.lm)
plast.desc.pau.mo.ar.lm <- data.frame(plast.desc.pau.mo.ar.lm)
plast.desc.pau.mo.speed.lm <- data.frame(plast.desc.pau.mo.speed.lm)
colnames(plast.desc.pau.mo.area.lm) <- c('sal.ori', 'sal.des', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.pau.mo.ar.lm) <- c('sal.ori', 'sal.des', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.pau.mo.speed.lm) <- c('sal.ori', 'sal.des', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')

## competition descendants
pos.120 <- which(dat09.PTS.Pau$ID_original == 120)
#dat09.PTS.Pau <- dat09.PTS.Pau[-pos.120,]

plast.desc.pau.mi.area.lm <- plast.desc.pau.mi.ar.lm <- plast.desc.pau.mi.speed.lm <- vector()
for(ori in 1:5){
	dat.tmp <- dat09.PTS.Pau[dat09.PTS.Pau$Salinity_Origin == sal.ori[ori],]

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

			plast.desc.pau.mi.area.lm <- rbind(plast.desc.pau.mi.area.lm, cbind(sal.ori[ori], sal.des[des], sum.area[2,1], sum.area[2,2], sum.area[2,5],
									sum.area[3,1], sum.area[3,2], sum.area[3,5]))
			plast.desc.pau.mi.ar.lm <- rbind(plast.desc.pau.mi.ar.lm, cbind(sal.ori[ori], sal.des[des], sum.ar[2,1], sum.ar[2,2], sum.ar[2,5],
									sum.ar[3,1], sum.ar[3,2], sum.ar[3,5]))
			plast.desc.pau.mi.speed.lm <- rbind(plast.desc.pau.mi.speed.lm, cbind(sal.ori[ori], sal.des[des], sum.speed[2,1], sum.speed[2,2], sum.speed[2,5],
									sum.speed[3,1], sum.speed[3,2], sum.speed[3,5]))
		}
	}
}
plast.desc.pau.mi.area.lm <- data.frame(plast.desc.pau.mi.area.lm)
plast.desc.pau.mi.ar.lm <- data.frame(plast.desc.pau.mi.ar.lm)
plast.desc.pau.mi.speed.lm <- data.frame(plast.desc.pau.mi.speed.lm)
colnames(plast.desc.pau.mi.area.lm) <- c('sal.ori', 'sal.des', 'plast.area', 'se.plast.area', 'p.plast.area', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.pau.mi.ar.lm) <- c('sal.ori', 'sal.des', 'plast.ar', 'se.plast.ar', 'p.plast.ar', 'dens', 'se.dens', 'p.dens')
colnames(plast.desc.pau.mi.speed.lm) <- c('sal.ori', 'sal.des', 'plast.speed', 'se.plast.speed', 'p.plast.speed', 'dens', 'se.dens', 'p.dens')

### -------------------------------------------------------------------------------------------------
## Phenotypic plasticity response to salinity -- Figure S7 
### -------------------------------------------------------------------------------------------------

## You need to run previous code of Table S9-S10 in order to make this figure
 
## area
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))
plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-700, 2000))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.pau.mo.area.lm[i,2])
	segments(i-0.5, plast.anc.pau.mo.area.lm[i,2]-plast.anc.pau.mo.area.lm[i,3], i-0.5, plast.anc.pau.mo.area.lm[i,2]+plast.anc.pau.mo.area.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.pau.mi.area.lm[i,2])
	segments(i+3, plast.anc.pau.mi.area.lm[i,2]-plast.anc.pau.mi.area.lm[i,3], i+3, plast.anc.pau.mi.area.lm[i,2]+plast.anc.pau.mi.area.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.pau.mo.area.lm[plast.desc.pau.mo.area.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.pau.mi.area.lm[plast.desc.pau.mi.area.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-700, 2000))
	abline(h = 0)
	
	for(des in 1:3){
		x <- plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]
		y <- plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]])
			segments(des-0.5, plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]-plast.mo.tmp$se.plast.area[plast.mo.tmp$sal.des == sal.des[des]],
					des-0.5, plast.mo.tmp$plast.area[plast.mo.tmp$sal.des == sal.des[des]]+plast.mo.tmp$se.plast.area[plast.mo.tmp$sal.des == sal.des[des]])
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.area[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.area[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.area[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}

## cell shape
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))
plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-1, 0.4))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.pau.mo.ar.lm[i,2])
	segments(i-0.5, plast.anc.pau.mo.ar.lm[i,2]-plast.anc.pau.mo.ar.lm[i,3], i-0.5, plast.anc.pau.mo.ar.lm[i,2]+plast.anc.pau.mo.ar.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.pau.mi.ar.lm[i,2])
	segments(i+3, plast.anc.pau.mi.ar.lm[i,2]-plast.anc.pau.mi.ar.lm[i,3], i+3, plast.anc.pau.mi.ar.lm[i,2]+plast.anc.pau.mi.ar.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.pau.mo.ar.lm[plast.desc.pau.mo.ar.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.pau.mi.ar.lm[plast.desc.pau.mi.ar.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-1, 0.4))
	abline(h = 0)
	
	for(des in 1:3){
		x <- plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]
		y <- plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]])
			segments(des-0.5, plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]-plast.mo.tmp$se.plast.ar[plast.mo.tmp$sal.des == sal.des[des]],
					des-0.5, plast.mo.tmp$plast.ar[plast.mo.tmp$sal.des == sal.des[des]]+plast.mo.tmp$se.plast.ar[plast.mo.tmp$sal.des == sal.des[des]])
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.ar[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.ar[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.ar[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}

## speed
dev.new(width = 12, height = 8)
par(mfrow = c(2,3))

plot(-1000, -1000, type = 'n', xlim = c(0,6.5), ylim = c(-900, 150))
abline(h = 0)
for(i in 1:2){
	rect(i-0.5-0.2, 0, i-0.5+0.2, plast.anc.pau.mo.speed.lm[i,2])
	segments(i-0.5, plast.anc.pau.mo.speed.lm[i,2]-plast.anc.pau.mo.speed.lm[i,3], i-0.5, plast.anc.pau.mo.speed.lm[i,2]+plast.anc.pau.mo.speed.lm[i,3])
}
for(i in 1:3){
	rect(i+3-0.2, 0, i+3+0.2, plast.anc.pau.mi.speed.lm[i,2])
	segments(i+3, plast.anc.pau.mi.speed.lm[i,2]-plast.anc.pau.mi.speed.lm[i,3], i+3, plast.anc.pau.mi.speed.lm[i,2]+plast.anc.pau.mi.speed.lm[i,3])
}

for(ori in 1:5){
	plast.mo.tmp <- plast.desc.pau.mo.speed.lm[plast.desc.pau.mo.speed.lm$sal.ori == sal.ori[ori],]
	plast.mi.tmp <- plast.desc.pau.mi.speed.lm[plast.desc.pau.mi.speed.lm$sal.ori == sal.ori[ori],]
	
	plot(-1000, -1000,  type = 'n', xlim = c(0,6.5), ylim = c(-900, 150))
	abline(h = 0)
	
	for(des in 1:3){
		x <- plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]
		y <- plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]

		if(length(x) > 0){
			rect(des-0.5-0.2, 0, des-0.5+0.2, plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]])
			segments(des-0.5, plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]-plast.mo.tmp$se.plast.speed[plast.mo.tmp$sal.des == sal.des[des]],
					des-0.5, plast.mo.tmp$plast.speed[plast.mo.tmp$sal.des == sal.des[des]]+plast.mo.tmp$se.plast.speed[plast.mo.tmp$sal.des == sal.des[des]])
		}

		if(length(y) > 0){
			rect(des+3-0.2, 0, des+3+0.2, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]])
			segments(des+3, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]-plast.mi.tmp$se.plast.speed[plast.mi.tmp$sal.des == sal.des[des]],
					des+3, plast.mi.tmp$plast.speed[plast.mi.tmp$sal.des == sal.des[des]]+plast.mi.tmp$se.plast.speed[plast.mi.tmp$sal.des == sal.des[des]])
		}
	}
}

### -------------------------------------------------------------------------------------------------
## Genetic trait difference for high salinity selected populations -- Table S13 
### -------------------------------------------------------------------------------------------------

## In the main text we show results excluding ID 120, run the commented code to get those results
## In Supplementary Figure S12 we show the results including ID 120, do not run the commented code to get those results
## Analysis genetic differences
#pos.120 <- which(dat09.combined.Pau$ID_original == 120)
#dat09.combined.Pau <- dat09.combined.Pau[-pos.120,]

tmp.sal2.05 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '2' & dat09.combined.Pau$Salinity_Destination == '0.5',]
Com <- factor(tmp.sal2.05$Community, level = c('mono', 'mixed'))
fit.area.205 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.05)
summary(fit.area.205)
fit.ar.205 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.05)
summary(fit.ar.205)
fit.speed.205 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.05)
summary(fit.speed.205)

tmp.sal2.1 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '2' & dat09.combined.Pau$Salinity_Destination == '1',]
Com <- factor(tmp.sal2.1$Community, level = c('mono', 'mixed'))
fit.area.21 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.1)
summary(fit.area.21)
fit.ar.21 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.1)
summary(fit.ar.21)
fit.speed.21 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.1)
summary(fit.speed.21)

tmp.sal2.2 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '2' & dat09.combined.Pau$Salinity_Destination == '2',]
Com <- factor(tmp.sal2.2$Community, level = c('mono', 'mixed'))
fit.area.22 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.2)
summary(fit.area.22)
fit.ar.22 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.2)
summary(fit.ar.22)
fit.speed.22 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.2)
summary(fit.speed.22)

tmp.sal2.4 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '2' & dat09.combined.Pau$Salinity_Destination == '4',]
Com <- factor(tmp.sal2.4$Community, level = c('mono', 'mixed'))
fit.area.24 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.4)
summary(fit.area.24)
fit.ar.24 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.4)
summary(fit.ar.24)
fit.speed.24 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal2.4)
summary(fit.speed.24)

tmp.sal4.05 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '4' & dat09.combined.Pau$Salinity_Destination == '0.5',]
Com <- factor(tmp.sal4.05$Community, level = c('mono', 'mixed'))
fit.area.405 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.05)
summary(fit.area.405)
fit.ar.405 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.05)
summary(fit.ar.405)
fit.speed.405 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.05)
summary(fit.speed.405)

tmp.sal4.1 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '4' & dat09.combined.Pau$Salinity_Destination == '1',]
Com <- factor(tmp.sal4.1$Community, level = c('mono', 'mixed'))
fit.area.41 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.1)
summary(fit.area.41)
fit.ar.41 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.1)
summary(fit.ar.41)
fit.speed.41 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.1)
summary(fit.speed.41)

tmp.sal4.2 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '4' & dat09.combined.Pau$Salinity_Destination == '2',]
Com <- factor(tmp.sal4.2$Community, level = c('mono', 'mixed'))
fit.area.42 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.2)
summary(fit.area.42)
fit.ar.42 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.2)
summary(fit.ar.42)
fit.speed.42 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.2)
summary(fit.speed.42)

tmp.sal4.4 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '4' & dat09.combined.Pau$Salinity_Destination == '4',]
Com <- factor(tmp.sal4.4$Community, level = c('mono', 'mixed'))
fit.area.44 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.4)
summary(fit.area.44)
fit.ar.44 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.4)
summary(fit.ar.44)
fit.speed.44 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal4.4)
summary(fit.speed.44)


### -------------------------------------------------------------------------------------------------
## Genetic trait difference for high salinity selected populations -- Figure 5 main text --> Supplementary Figure S12
### -------------------------------------------------------------------------------------------------

## Plot the effect size
par(mfrow = c(1,2))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-1100, 1400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.area.205)$coefficients[2,1])
segments(1, summary(fit.area.205)$coefficients[2,1]-summary(fit.area.205)$coefficients[2,2], 1, summary(fit.area.205)$coefficients[2,1]+summary(fit.area.205)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.area.21)$coefficients[2,1])
segments(2, summary(fit.area.21)$coefficients[2,1]-summary(fit.area.21)$coefficients[2,2], 2, summary(fit.area.21)$coefficients[2,1]+summary(fit.area.21)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.area.22)$coefficients[2,1])
segments(3, summary(fit.area.22)$coefficients[2,1]-summary(fit.area.22)$coefficients[2,2], 3, summary(fit.area.22)$coefficients[2,1]+summary(fit.area.22)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.area.24)$coefficients[2,1])
segments(4, summary(fit.area.24)$coefficients[2,1]-summary(fit.area.24)$coefficients[2,2], 4, summary(fit.area.24)$coefficients[2,1]+summary(fit.area.24)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-1100, 1400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.area.405)$coefficients[2,1])
segments(1, summary(fit.area.405)$coefficients[2,1]-summary(fit.area.405)$coefficients[2,2], 1, summary(fit.area.405)$coefficients[2,1]+summary(fit.area.405)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.area.41)$coefficients[2,1])
segments(2, summary(fit.area.41)$coefficients[2,1]-summary(fit.area.41)$coefficients[2,2], 2, summary(fit.area.41)$coefficients[2,1]+summary(fit.area.41)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.area.42)$coefficients[2,1])
segments(3, summary(fit.area.42)$coefficients[2,1]-summary(fit.area.42)$coefficients[2,2], 3, summary(fit.area.42)$coefficients[2,1]+summary(fit.area.42)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.area.44)$coefficients[2,1])
segments(4, summary(fit.area.44)$coefficients[2,1]-summary(fit.area.44)$coefficients[2,2], 4, summary(fit.area.44)$coefficients[2,1]+summary(fit.area.44)$coefficients[2,2])


par(mfrow = c(1,2))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-0.6, 0.5))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.ar.205)$coefficients[2,1])
segments(1, summary(fit.ar.205)$coefficients[2,1]-summary(fit.ar.205)$coefficients[2,2], 1, summary(fit.ar.205)$coefficients[2,1]+summary(fit.ar.205)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.ar.21)$coefficients[2,1])
segments(2, summary(fit.ar.21)$coefficients[2,1]-summary(fit.ar.21)$coefficients[2,2], 2, summary(fit.ar.21)$coefficients[2,1]+summary(fit.ar.21)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.ar.22)$coefficients[2,1])
segments(3, summary(fit.ar.22)$coefficients[2,1]-summary(fit.ar.22)$coefficients[2,2], 3, summary(fit.ar.22)$coefficients[2,1]+summary(fit.ar.22)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.ar.24)$coefficients[2,1])
segments(4, summary(fit.ar.24)$coefficients[2,1]-summary(fit.ar.24)$coefficients[2,2], 4, summary(fit.ar.24)$coefficients[2,1]+summary(fit.ar.24)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-0.6, 0.5))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.ar.405)$coefficients[2,1])
segments(1, summary(fit.ar.405)$coefficients[2,1]-summary(fit.ar.405)$coefficients[2,2], 1, summary(fit.ar.405)$coefficients[2,1]+summary(fit.ar.405)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.ar.41)$coefficients[2,1])
segments(2, summary(fit.ar.41)$coefficients[2,1]-summary(fit.ar.41)$coefficients[2,2], 2, summary(fit.ar.41)$coefficients[2,1]+summary(fit.ar.41)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.ar.42)$coefficients[2,1])
segments(3, summary(fit.ar.42)$coefficients[2,1]-summary(fit.ar.42)$coefficients[2,2], 3, summary(fit.ar.42)$coefficients[2,1]+summary(fit.ar.42)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.ar.44)$coefficients[2,1])
segments(4, summary(fit.ar.44)$coefficients[2,1]-summary(fit.ar.44)$coefficients[2,2], 4, summary(fit.ar.44)$coefficients[2,1]+summary(fit.ar.44)$coefficients[2,2])

par(mfrow = c(1,2))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-800, 400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.speed.205)$coefficients[2,1])
segments(1, summary(fit.speed.205)$coefficients[2,1]-summary(fit.speed.205)$coefficients[2,2], 1, summary(fit.speed.205)$coefficients[2,1]+summary(fit.speed.205)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.speed.21)$coefficients[2,1])
segments(2, summary(fit.speed.21)$coefficients[2,1]-summary(fit.speed.21)$coefficients[2,2], 2, summary(fit.speed.21)$coefficients[2,1]+summary(fit.speed.21)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.speed.22)$coefficients[2,1])
segments(3, summary(fit.speed.22)$coefficients[2,1]-summary(fit.speed.22)$coefficients[2,2], 3, summary(fit.speed.22)$coefficients[2,1]+summary(fit.speed.22)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.speed.24)$coefficients[2,1])
segments(4, summary(fit.speed.24)$coefficients[2,1]-summary(fit.speed.24)$coefficients[2,2], 4, summary(fit.speed.24)$coefficients[2,1]+summary(fit.speed.24)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-800, 400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.speed.405)$coefficients[2,1])
segments(1, summary(fit.speed.405)$coefficients[2,1]-summary(fit.speed.405)$coefficients[2,2], 1, summary(fit.speed.405)$coefficients[2,1]+summary(fit.speed.405)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.speed.41)$coefficients[2,1])
segments(2, summary(fit.speed.41)$coefficients[2,1]-summary(fit.speed.41)$coefficients[2,2], 2, summary(fit.speed.41)$coefficients[2,1]+summary(fit.speed.41)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.speed.42)$coefficients[2,1])
segments(3, summary(fit.speed.42)$coefficients[2,1]-summary(fit.speed.42)$coefficients[2,2], 3, summary(fit.speed.42)$coefficients[2,1]+summary(fit.speed.42)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.speed.44)$coefficients[2,1])
segments(4, summary(fit.speed.44)$coefficients[2,1]-summary(fit.speed.44)$coefficients[2,2], 4, summary(fit.speed.44)$coefficients[2,1]+summary(fit.speed.44)$coefficients[2,2])


### -------------------------------------------------------------------------------------------------
## Trait difference for lower salinity selected populations -- Table S14 
### -------------------------------------------------------------------------------------------------

tmp.sal0.0 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0' & dat09.combined.Pau$Salinity_Destination == '0',]
Com <- factor(tmp.sal0.0$Community, level = c('mono', 'mixed'))
fit.area.00 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.0)
summary(fit.area.00)
fit.ar.00 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.0)
summary(fit.ar.00)
fit.speed.00 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.0)
summary(fit.speed.00)

tmp.sal0.05 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0' & dat09.combined.Pau$Salinity_Destination == '0.5',]
Com <- factor(tmp.sal0.05$Community, level = c('mono', 'mixed'))
fit.area.005 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.05)
summary(fit.area.005)
fit.ar.005 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.05)
summary(fit.ar.005)
fit.speed.005 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.05)
summary(fit.speed.005)

tmp.sal0.1 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0' & dat09.combined.Pau$Salinity_Destination == '1',]
Com <- factor(tmp.sal0.1$Community, level = c('mono', 'mixed'))
fit.area.01 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.1)
summary(fit.area.01)
fit.ar.01 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.1)
summary(fit.ar.01)
fit.speed.01 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.1)
summary(fit.speed.01)

tmp.sal0.2 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0' & dat09.combined.Pau$Salinity_Destination == '2',]
Com <- factor(tmp.sal0.2$Community, level = c('mono', 'mixed'))
fit.area.02 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.2)
summary(fit.area.02)
fit.ar.02 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.2)
summary(fit.ar.02)
fit.speed.02 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.2)
summary(fit.speed.02)

tmp.sal0.4 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0' & dat09.combined.Pau$Salinity_Destination == '4',]
Com <- factor(tmp.sal0.4$Community, level = c('mono', 'mixed'))
fit.area.04 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.4)
summary(fit.area.04)
fit.ar.04 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.4)
summary(fit.ar.04)
fit.speed.04 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal0.4)
summary(fit.speed.04)

## Salinity 0.5
tmp.sal05.05 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0.5' & dat09.combined.Pau$Salinity_Destination == '0.5',]
Com <- factor(tmp.sal05.05$Community, level = c('mono', 'mixed'))
fit.area.0505 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.05)
summary(fit.area.0505)
fit.ar.0505 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.05)
summary(fit.ar.0505)
fit.speed.0505 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.05)
summary(fit.speed.0505)

tmp.sal05.1 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0.5' & dat09.combined.Pau$Salinity_Destination == '1',]
Com <- factor(tmp.sal05.1$Community, level = c('mono', 'mixed'))
fit.area.051 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.1)
summary(fit.area.051)
fit.ar.051 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.1)
summary(fit.ar.051)
fit.speed.051 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.1)
summary(fit.speed.051)

tmp.sal05.2 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0.5' & dat09.combined.Pau$Salinity_Destination == '2',]
Com <- factor(tmp.sal05.2$Community, level = c('mono', 'mixed'))
fit.area.052 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.2)
summary(fit.area.052)
fit.ar.052 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.2)
summary(fit.ar.052)
fit.speed.052 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.2)
summary(fit.speed.052)

tmp.sal05.4 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '0.5' & dat09.combined.Pau$Salinity_Destination == '4',]
Com <- factor(tmp.sal05.4$Community, level = c('mono', 'mixed'))
fit.area.054 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.4)
summary(fit.area.054)
fit.ar.054 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.4)
summary(fit.ar.054)
fit.speed.054 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal05.4)
summary(fit.speed.054)

## Salinity 1
tmp.sal1.05 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '1' & dat09.combined.Pau$Salinity_Destination == '0.5',]
Com <- factor(tmp.sal1.05$Community, level = c('mono', 'mixed'))
fit.area.105 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.05)
summary(fit.area.105)
fit.ar.105 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.05)
summary(fit.ar.105)
fit.speed.105 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.05)
summary(fit.speed.105)

tmp.sal1.1 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '1' & dat09.combined.Pau$Salinity_Destination == '1',]
Com <- factor(tmp.sal1.1$Community, level = c('mono', 'mixed'))
fit.area.11 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.1)
summary(fit.area.11)
fit.ar.11 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.1)
summary(fit.ar.11)
fit.speed.11 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.1)
summary(fit.speed.11)

tmp.sal1.2 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '1' & dat09.combined.Pau$Salinity_Destination == '2',]
Com <- factor(tmp.sal1.2$Community, level = c('mono', 'mixed'))
fit.area.12 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.2)
summary(fit.area.12)
fit.ar.12 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.2)
summary(fit.ar.12)
fit.speed.12 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.2)
summary(fit.speed.12)

tmp.sal1.4 <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Origin == '1' & dat09.combined.Pau$Salinity_Destination == '4',]
Com <- factor(tmp.sal1.4$Community, level = c('mono', 'mixed'))
fit.area.14 <- lmer(mean_area ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.4)
summary(fit.area.14)
fit.ar.14 <- lmer(mean_ar ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.4)
summary(fit.ar.14)
fit.speed.14 <- lmer(gross_speed ~ Com + Density + (1|ID_original) + (1|ID_original:replicate), data = tmp.sal1.4)
summary(fit.speed.14)

### -------------------------------------------------------------------------------------------------
## Trait difference for lower salinity selected populations -- Figure S9 
### -------------------------------------------------------------------------------------------------

## You need to run the code of Table S14 in order to make this figure

## Plot the effect size
par(mfrow = c(1,3))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-1100, 1400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.area.00)$coefficients[2,1])
segments(1, summary(fit.area.00)$coefficients[2,1]-summary(fit.area.00)$coefficients[2,2], 1, summary(fit.area.00)$coefficients[2,1]+summary(fit.area.00)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.area.005)$coefficients[2,1])
segments(2, summary(fit.area.005)$coefficients[2,1]-summary(fit.area.005)$coefficients[2,2], 2, summary(fit.area.005)$coefficients[2,1]+summary(fit.area.005)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.area.01)$coefficients[2,1])
segments(3, summary(fit.area.01)$coefficients[2,1]-summary(fit.area.01)$coefficients[2,2], 3, summary(fit.area.01)$coefficients[2,1]+summary(fit.area.01)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.area.02)$coefficients[2,1])
segments(4, summary(fit.area.02)$coefficients[2,1]-summary(fit.area.02)$coefficients[2,2], 4, summary(fit.area.02)$coefficients[2,1]+summary(fit.area.02)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-1100, 1400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.area.0505)$coefficients[2,1])
segments(1, summary(fit.area.0505)$coefficients[2,1]-summary(fit.area.0505)$coefficients[2,2], 1, summary(fit.area.0505)$coefficients[2,1]+summary(fit.area.0505)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.area.051)$coefficients[2,1])
segments(2, summary(fit.area.051)$coefficients[2,1]-summary(fit.area.051)$coefficients[2,2], 2, summary(fit.area.051)$coefficients[2,1]+summary(fit.area.051)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.area.052)$coefficients[2,1])
segments(3, summary(fit.area.052)$coefficients[2,1]-summary(fit.area.052)$coefficients[2,2], 3, summary(fit.area.052)$coefficients[2,1]+summary(fit.area.052)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-1100, 1400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.area.105)$coefficients[2,1])
segments(1, summary(fit.area.105)$coefficients[2,1]-summary(fit.area.105)$coefficients[2,2], 1, summary(fit.area.105)$coefficients[2,1]+summary(fit.area.105)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.area.11)$coefficients[2,1])
segments(2, summary(fit.area.11)$coefficients[2,1]-summary(fit.area.11)$coefficients[2,2], 2, summary(fit.area.11)$coefficients[2,1]+summary(fit.area.11)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.area.12)$coefficients[2,1])
segments(3, summary(fit.area.12)$coefficients[2,1]-summary(fit.area.12)$coefficients[2,2], 3, summary(fit.area.12)$coefficients[2,1]+summary(fit.area.12)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.area.14)$coefficients[2,1])
segments(4, summary(fit.area.14)$coefficients[2,1]-summary(fit.area.14)$coefficients[2,2], 4, summary(fit.area.14)$coefficients[2,1]+summary(fit.area.14)$coefficients[2,2])


par(mfrow = c(1,3))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-0.7, 0.5))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.ar.00)$coefficients[2,1])
segments(1, summary(fit.ar.00)$coefficients[2,1]-summary(fit.ar.00)$coefficients[2,2], 1, summary(fit.ar.00)$coefficients[2,1]+summary(fit.ar.00)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.ar.005)$coefficients[2,1])
segments(2, summary(fit.ar.005)$coefficients[2,1]-summary(fit.ar.005)$coefficients[2,2], 2, summary(fit.ar.005)$coefficients[2,1]+summary(fit.ar.005)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.ar.01)$coefficients[2,1])
segments(3, summary(fit.ar.01)$coefficients[2,1]-summary(fit.ar.01)$coefficients[2,2], 3, summary(fit.ar.01)$coefficients[2,1]+summary(fit.ar.01)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.ar.02)$coefficients[2,1])
segments(4, summary(fit.ar.02)$coefficients[2,1]-summary(fit.ar.02)$coefficients[2,2], 4, summary(fit.ar.02)$coefficients[2,1]+summary(fit.ar.02)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-0.7, 0.5))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.ar.0505)$coefficients[2,1])
segments(1, summary(fit.ar.0505)$coefficients[2,1]-summary(fit.ar.0505)$coefficients[2,2], 1, summary(fit.ar.0505)$coefficients[2,1]+summary(fit.ar.0505)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.ar.051)$coefficients[2,1])
segments(2, summary(fit.ar.051)$coefficients[2,1]-summary(fit.ar.051)$coefficients[2,2], 2, summary(fit.ar.051)$coefficients[2,1]+summary(fit.ar.051)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.ar.052)$coefficients[2,1])
segments(3, summary(fit.ar.052)$coefficients[2,1]-summary(fit.ar.052)$coefficients[2,2], 3, summary(fit.ar.052)$coefficients[2,1]+summary(fit.ar.052)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-0.7, 0.5))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.ar.105)$coefficients[2,1])
segments(1, summary(fit.ar.105)$coefficients[2,1]-summary(fit.ar.105)$coefficients[2,2], 1, summary(fit.ar.105)$coefficients[2,1]+summary(fit.ar.105)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.ar.11)$coefficients[2,1])
segments(2, summary(fit.ar.11)$coefficients[2,1]-summary(fit.ar.11)$coefficients[2,2], 2, summary(fit.ar.11)$coefficients[2,1]+summary(fit.ar.11)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.ar.12)$coefficients[2,1])
segments(3, summary(fit.ar.12)$coefficients[2,1]-summary(fit.ar.12)$coefficients[2,2], 3, summary(fit.ar.12)$coefficients[2,1]+summary(fit.ar.12)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.ar.14)$coefficients[2,1])
segments(4, summary(fit.ar.14)$coefficients[2,1]-summary(fit.ar.14)$coefficients[2,2], 4, summary(fit.ar.14)$coefficients[2,1]+summary(fit.ar.14)$coefficients[2,2])



par(mfrow = c(1,3))
plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-800, 400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.speed.00)$coefficients[2,1])
segments(1, summary(fit.speed.00)$coefficients[2,1]-summary(fit.speed.00)$coefficients[2,2], 1, summary(fit.speed.00)$coefficients[2,1]+summary(fit.speed.00)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.speed.005)$coefficients[2,1])
segments(2, summary(fit.speed.005)$coefficients[2,1]-summary(fit.speed.005)$coefficients[2,2], 2, summary(fit.speed.005)$coefficients[2,1]+summary(fit.speed.005)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.speed.01)$coefficients[2,1])
segments(3, summary(fit.speed.01)$coefficients[2,1]-summary(fit.speed.01)$coefficients[2,2], 3, summary(fit.speed.01)$coefficients[2,1]+summary(fit.speed.01)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.speed.02)$coefficients[2,1])
segments(4, summary(fit.speed.02)$coefficients[2,1]-summary(fit.speed.02)$coefficients[2,2], 4, summary(fit.speed.02)$coefficients[2,1]+summary(fit.speed.02)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-800, 400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.speed.0505)$coefficients[2,1])
segments(1, summary(fit.speed.0505)$coefficients[2,1]-summary(fit.speed.0505)$coefficients[2,2], 1, summary(fit.speed.0505)$coefficients[2,1]+summary(fit.speed.0505)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.speed.051)$coefficients[2,1])
segments(2, summary(fit.speed.051)$coefficients[2,1]-summary(fit.speed.051)$coefficients[2,2], 2, summary(fit.speed.051)$coefficients[2,1]+summary(fit.speed.051)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.speed.052)$coefficients[2,1])
segments(3, summary(fit.speed.052)$coefficients[2,1]-summary(fit.speed.052)$coefficients[2,2], 3, summary(fit.speed.052)$coefficients[2,1]+summary(fit.speed.052)$coefficients[2,2])

plot(-2000, -2000, type = 'n',  xlim = c(0,5), ylim = c(-800, 400))
abline(h=0)
rect(1-0.4, 0, 1+0.4, summary(fit.speed.105)$coefficients[2,1])
segments(1, summary(fit.speed.105)$coefficients[2,1]-summary(fit.speed.105)$coefficients[2,2], 1, summary(fit.speed.105)$coefficients[2,1]+summary(fit.speed.105)$coefficients[2,2])
rect(2-0.4, 0, 2+0.4, summary(fit.speed.11)$coefficients[2,1])
segments(2, summary(fit.speed.11)$coefficients[2,1]-summary(fit.speed.11)$coefficients[2,2], 2, summary(fit.speed.11)$coefficients[2,1]+summary(fit.speed.11)$coefficients[2,2])
rect(3-0.4, 0, 3+0.4, summary(fit.speed.12)$coefficients[2,1])
segments(3, summary(fit.speed.12)$coefficients[2,1]-summary(fit.speed.12)$coefficients[2,2], 3, summary(fit.speed.12)$coefficients[2,1]+summary(fit.speed.12)$coefficients[2,2])
rect(4-0.4, 0, 4+0.4, summary(fit.speed.14)$coefficients[2,1])
segments(4, summary(fit.speed.14)$coefficients[2,1]-summary(fit.speed.14)$coefficients[2,2], 4, summary(fit.speed.14)$coefficients[2,1]+summary(fit.speed.14)$coefficients[2,2])


### -------------------------------------------------------------------------------------------------
## Checking character displacement -- Figure S10
### -------------------------------------------------------------------------------------------------

dat09.combined.Pau <- rbind(dat09.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Pau[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])
dat09.combined.Spite <- rbind(dat09.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,39:43)], dat09.PTS.Spite[,c(30,32:33,35,4,6,8,10,12,14,19:22,40:44)])

dat.pau <- dat09.combined.Pau[dat09.combined.Pau$Salinity_Destination == 0.5,]
dat.spite <- dat09.combined.Spite[dat09.combined.Spite$Salinity_Destination == 0.5,]


sal.lvl <- c(0,0.5,1,2,4)
plot(-100, -100, ylim = c(1000,10100), xlim = c(0,4), main = 'area')
for(ss in 1:5){
	pau.mo <- dat.pau[dat.pau$Community == 'mono' & dat.pau$Salinity_Origin == sal.lvl[ss],]
	spite.mo <- dat.spite[dat.spite$Community == 'mono' & dat.spite$Salinity_Origin == sal.lvl[ss],]

	id.pau.mo <- unique(pau.mo$ID_original)
	id.spite.mo <- unique(spite.mo$ID_original)

	for(id in 1:length(id.pau.mo)){
		tt. <- mean(pau.mo$mean_area[pau.mo$ID_original == id.pau.mo[id]])
		sd. <- sd(pau.mo$mean_area[pau.mo$ID_original== id.pau.mo[id]])
		points(sal.lvl[ss], tt., pch = 21, cex = 2)
		segments(sal.lvl[ss], tt.-sd., sal.lvl[ss], tt.+sd.)
	}
	for(id in 1:length(id.spite.mo)){
		tt. <- mean(spite.mo$mean_area[spite.mo$ID_original== id.spite.mo[id]])
		sd. <- sd(spite.mo$mean_area[spite.mo$ID_original== id.spite.mo[id]])
		points(sal.lvl[ss], tt., pch = 22, cex = 2)
		segments(sal.lvl[ss], tt.-sd., sal.lvl[ss], tt.+sd.)
	}
}

sal.lvl <- c(0,0.5,1,2,4) 
plot(-100, -100, ylim = c(2,8.2), xlim = c(0,4), main = 'ar')
for(ss in 1:5){
	pau.mo <- dat.pau[dat.pau$Community == 'mono' & dat.pau$Salinity_Origin == sal.lvl[ss],]
	spite.mo <- dat.spite[dat.spite$Community == 'mono' & dat.spite$Salinity_Origin == sal.lvl[ss],]

	id.pau.mo <- unique(pau.mo$ID_original)
	id.spite.mo <- unique(spite.mo$ID_original)

	for(id in 1:length(id.pau.mo)){
		tt. <- mean(pau.mo$mean_ar[pau.mo$ID_original == id.pau.mo[id]])
		sd. <- sd(pau.mo$mean_ar[pau.mo$ID_original== id.pau.mo[id]])
		points(sal.lvl[ss], tt., pch = 21, cex = 2)
		segments(sal.lvl[ss], tt.-sd., sal.lvl[ss], tt.+sd.)
	}
	for(id in 1:length(id.spite.mo)){
		tt. <- mean(spite.mo$mean_ar[spite.mo$ID_original== id.spite.mo[id]])
		sd. <- sd(spite.mo$mean_ar[spite.mo$ID_original== id.spite.mo[id]])
		points(sal.lvl[ss], tt., pch = 22, cex = 2)
		segments(sal.lvl[ss], tt.-sd., sal.lvl[ss], tt.+sd.)
	}
}

### -------------------------------------------------------------------------------------------------
## Plotting the community composition -- Figure S11
### -------------------------------------------------------------------------------------------------

#### day 4
ComComp.plot(dat=list(mean22.PTS.Pau.CG, mean22.PTS.Spite.CG), show.biomass=TRUE, CG=FALSE)
#### day 78
ComComp.plot(dat=list(mean05.PTS.Pau), show.biomass=TRUE, CG=FALSE)
#### day 82
ComComp.plot(dat=mean.PTS09.Pau, show.biomass=TRUE)






