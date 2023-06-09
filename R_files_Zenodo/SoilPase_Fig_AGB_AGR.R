# NSF FX data plotting

# Figure 1
# Aboveground Biomass absolute growth rate from 2016-2018 for Pase paper

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

library(colorBlindness)

#########################################################################
############################# Set up figure #############################
#########################################################################

#pdf(file="Pase_Fig01.pdf",
#	width=5,height=4)
setEPS()
postscript(file="Pase_Fig01.eps",
           width=5,height=4)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,4),omi=c(.4,.25,.2,.2),mai=c(0,.3,0,0))

##############################################
################## New York ##################
##############################################

# Run code to get the data to plot
source("NY_growth_analysis.R")
# Run code for plotting specifics
source("Fig_3-6_setup.R")

nonfixcol <- Blue2DarkRed12Steps[1]
rhizcol <- Blue2DarkRed12Steps[11]
actcol <- Blue2DarkRed12Steps[9]

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_02 == 0,]$AGR_AGB_2015_2016 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2015_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_AGB_2015_2018 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_AGB_2017_2018 <- NA

dat[dat$Use_growth_01 == 0,]$AGB_est_kg_01 <- NA
dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_07 == 0,]$AGB_est_kg_07 <- NA

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016),]$AGR_AGB_2015_2016,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.002,2),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5,cex.axis=0.8)
axis(side=2, at=c(0.01,1),labels=c("0.01","1"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=c(2),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016),]$AGR_AGB_2015_2016,
	col=nonfixcol,pch=1,cex=0.8)	
mtext("a",at=-0.35,padj=1.5,cex=0.8)
mtext("2015-2016",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.005,15),log="y")
axis(side=2, at=c(0.01,1,10),labels=c("0.01","1","10"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.005,0.009,0.001),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=nonfixcol,pch=1,cex=0.8)	
mtext("b",at=-0.35,padj=1.5,cex=0.8)
mtext("2016-2017",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,50),log="y")
axis(side=2, at=c(0.01,1,10),labels=c("0.01","1","10"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=nonfixcol,pch=1,cex=0.8)	
mtext("c",at=-0.35,padj=1.5,cex=0.8)
mtext("2017-2018",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2018),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2018),]$AGR_AGB_2015_2018,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.03,50),log="y")
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"),cex.axis=0.8)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd.ticks = 0.25)
ct <- summary(lme_ROPS_AGB_AGR_2015678)[[20]] # Coefficient table
ROPSmean <- mean(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ROPSmean,rhizcol,17,TRUE)
BENImean <- mean(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_BENI_AGB_AGR_20156789)[[20]],
	BENImean,nonfixcol,16,TRUE)
plot_siglets(-0.1,c(6,18,18,18),c("a","b","b","ab"),nonfixcol)
mtext("d",at=-0.35,padj=1.5,cex=0.8)
mtext("All years",at=1.5,cex=0.8)
mtext("New York",side=4,padj=0.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

##############################################
################### Oregon ###################
##############################################

source("OR_growth_analysis.R")
source("Fig_3-6_setup.R")

nonfixcol <- Blue2DarkRed12Steps[1]
rhizcol <- Blue2DarkRed12Steps[11]
actcol <- Blue2DarkRed12Steps[9]

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA

dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA

plot.new()

points(0.2+c(0.1,0.1),c(.6,.45),pch=c(2,1),col=c(rhizcol,nonfixcol),cex=0.8,xpd=TRUE)
text(0.2+0.1,0.525,"a-d",pos=2,font=1,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.6,"Robinia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.45,"Betula",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=0.8)

points(0.2+c(0.1,0.1),c(.2,.05),pch=c(0,1),col=c(actcol,nonfixcol),cex=0.8,xpd=TRUE)
text(0.2+0.1,0.125,"e-g",pos=2,font=1,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.2,"Alnus",pos=4,col=actcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.05,"Pseudotsuga",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=0.8)

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=actcol,pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.05,3),log="y")
axis(side=2, at=c(0.1,1),labels=c("0.1","1"),cex.axis=0.8)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,3,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=nonfixcol,pch=1,cex=0.8)	
mtext("e",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=actcol,pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.3,10),log="y")
axis(side=2, at=c(1,2,5),labels=c("1","2","5"),cex.axis=0.8)
axis(side=2, at=seq(0.3,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=nonfixcol,pch=1,cex=0.8)	
mtext("f",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2018),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2018),]$AGR_AGB_2016_2018,
	col="white",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.8,3),log="y")
axis(side=2, at=c(1,2,3,4),labels=TRUE,cex.axis=0.8)
axis(side=2, at=seq(1,5,1),labels=FALSE,lwd.ticks = 0.25)
ct <- summary(lme_ALRU_AGB_AGR_201678)[[20]] # Coefficient table
ALRUmean <- mean(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ALRUmean,actcol,15,TRUE)
PSMEmean <- mean(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSME_AGB_AGR_20167890)[[20]],
	PSMEmean,nonfixcol,16,TRUE)
mtext("Oregon",side=4,padj=0.5,cex=0.8)
mtext("g",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Waiakea ###################
###############################################

source("HI_W_growth_analysis.R")
source("Fig_3-6_setup.R")

nonfixcol <- Blue2DarkRed12Steps[1]
rhizcol <- Blue2DarkRed12Steps[11]
actcol <- Blue2DarkRed12Steps[9]

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_03 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA

dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA

plot.new()

points(0.2+c(0.1,0.1,0.1),c(.65,.5,.35),pch=c(2,0,1),col=c(rhizcol,actcol,nonfixcol),cex=0.8,xpd=TRUE)
text(0.2+0.1,0.5,"h-j",pos=2,font=1,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.65,"Gliricidia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.5,"Casuarina",pos=4,col=actcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.35,"Psidium",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=0.8)

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.002,5),log="y")
axis(side=2, at=c(0.01,1),labels=c("0.01","1"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,5,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=nonfixcol,pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=actcol,pch=0,cex=0.8)	
mtext("h",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.003,14),log="y")
axis(side=2, at=c(0.01,1,10),labels=c("0.01","1","10"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.003,0.009,0.001),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=nonfixcol,pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=actcol,pch=0,cex=0.8)	
mtext("i",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2018),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2018),]$AGR_AGB_2016_2018,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.04,3),log="y")
axis(side=2, at=c(0.3,3),labels=c("0.3","3"),cex.axis=0.8)
axis(side=2, at=c(1),labels=c("1"),cex.axis=0.8)
axis(side=2, at=seq(0.1,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd.ticks = 0.25)
ct <- summary(lme_GLSE_AGB_AGR_201678)[[20]] # Coefficient table
GLSEmean <- mean(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(0.1,ct,GLSEmean,rhizcol,17,TRUE)
PSCAmean <- mean(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSCA_AGB_AGR_2016789)[[20]],
	PSCAmean,nonfixcol,16,TRUE)
CAEQmean <- mean(dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(0,summary(lme_CAEQ_AGB_AGR_2016789)[[20]],
	CAEQmean,actcol,15,TRUE)
mtext("j",at=-0.25,padj=1.5,cex=0.8)
mtext("Waiakea",side=4,padj=0.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Volcano ###################
###############################################

source("HI_V_growth_analysis.R")
source("Fig_3-6_setup.R")

nonfixcol <- Blue2DarkRed12Steps[1]
rhizcol <- Blue2DarkRed12Steps[11]
actcol <- Blue2DarkRed12Steps[9]

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_03 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA

dat[dat$Use_growth_02 == 0 & dat$Use_growth_01 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA

plot.new()

points(0.2+c(0.1,0.1,0.1),c(.65,.5,.35),pch=c(2,0,1),col=c(rhizcol,actcol,nonfixcol),cex=0.8,xpd=TRUE)
text(0.2+0.1,0.5,"k-m",pos=2,font=1,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.65,"Acacia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.5,"Morella",pos=4,col=actcol,font=3,xpd=TRUE,cex=0.8)
text(0.2+0.1,0.35,"Dodonaea",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=0.8)

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.001,10),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5,cex.axis=0.8)
axis(side=2, at=c(0.001,0.1,10),labels=c("0.001","0.1","10"),cex.axis=0.8)
axis(side=2, at=c(1),labels=c("1"),cex.axis=0.8)
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=nonfixcol,pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col=actcol,pch=0,cex=0.8)	
mtext("k",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=rhizcol,pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,30),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5,cex.axis=0.8)
axis(side=2, at=c(0.01,1,10),labels=c("0.01","1","10"),cex.axis=0.8)
axis(side=2, at=c(0.1),labels=c("0.1"),cex.axis=0.8)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(20,30,10),labels=FALSE,lwd.ticks = 0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=nonfixcol,pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col=actcol,pch=0,cex=0.8)	
mtext("l",at=-0.35,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2018),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2018),]$AGR_AGB_2016_2018,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.06,8),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5,cex.axis=0.8)
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"),cex.axis=0.8)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd.ticks = 0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd.ticks = 0.25)
ct <- summary(lme_ACKO_AGB_AGR_201678)[[20]] # Coefficient table
ACKOmean <- mean(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ACKOmean,rhizcol,17,TRUE)
DOVImean <- mean(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_DOVI_AGB_AGR_2016789)[[20]],
	DOVImean,nonfixcol,16,TRUE)
MOFAmean <- mean(dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,na.rm=TRUE)
plot_u_se_lme(0,summary(lme_MOFA_AGB_AGR_2016789)[[20]],
	MOFAmean,actcol,15,TRUE)
mtext("m",at=-0.15,padj=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

mtext("Treatment", side=1, outer=T, at=0.52, padj=2.75,cex=0.8)
mtext(expression(Absolute ~ growth ~ rate ~ of ~ AGB ~ (kg ~ y^-1)),
	side=2, outer=T, at=0.5, padj=0,cex=0.8)
mtext("Volcano",side=4,padj=0.5,cex=0.8)


dev.off()

#############################################################
#############################################################
#############################################################