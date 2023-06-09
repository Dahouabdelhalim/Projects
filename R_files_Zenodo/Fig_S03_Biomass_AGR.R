# NSF FX data plotting

# Figure S3
# Total Biomass absolute growth rate

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS03.pdf",
	width=7.05,height=3)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(2,5),omi=c(.4,.25,.2,.2),mai=c(0,.3,0,0))

##############################################
################## New York ##################
##############################################

# Run code to get the data to plot
source("NY_growth_analysis.R")
# Run code for plotting specifics
source("Fig_3-6_setup.R")

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_02 == 0,]$AGR_Biomass_2015_2016 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_Biomass_2015_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_Biomass_2016_2017 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_Biomass_2015_2018 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_Biomass_2017_2018 <- NA
dat[dat$Use_growth_08 == 0,]$AGR_Biomass_2015_2019 <- NA
dat[dat$Use_growth_08 == 0,]$AGR_Biomass_2018_2019 <- NA

dat[dat$Use_growth_01 == 0,]$Biomass_est_kg_01 <- NA
dat[dat$Use_growth_02 == 0,]$Biomass_est_kg_02 <- NA
dat[dat$Use_growth_04 == 0,]$Biomass_est_kg_04 <- NA
dat[dat$Use_growth_07 == 0,]$Biomass_est_kg_07 <- NA
dat[dat$Use_growth_08 == 0,]$Biomass_est_kg_08 <- NA

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2016),]$AGR_Biomass_2015_2016,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2015_2016,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2015_2016,na.rm=TRUE)),
	log="y")
axis(side=2, at=c(0.01,0.1,1),labels=c("0.01","0.1","1"))
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(2),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2015_2016),]$AGR_Biomass_2015_2016,
	col="blue",pch=1,cex=0.8)	
mtext("a",at=-0.35,padj=1.5)
mtext("Years 1-2",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2016_2017),]$AGR_Biomass_2016_2017,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2016_2017,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2016_2017,na.rm=TRUE)),
	log="y")
axis(side=2, at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.005,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2016_2017),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2016_2017),]$AGR_Biomass_2016_2017,
	col="blue",pch=1,cex=0.8)	
mtext("b",at=-0.35,padj=1.5)
mtext("Years 2-3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2017_2018),]$AGR_Biomass_2017_2018,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2017_2018,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2017_2018,na.rm=TRUE)),
	log="y")
axis(side=2, at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2017_2018),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2017_2018),]$AGR_Biomass_2017_2018,
	col="blue",pch=1,cex=0.8)	
mtext("c",at=-0.35,padj=1.5)
mtext("Years 3-4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2018_2019),]$AGR_Biomass_2018_2019,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2018_2019,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2018_2019,na.rm=TRUE)),
	log="y")
axis(side=2, at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2018_2019),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_Biomass_2018_2019),]$AGR_Biomass_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("d",at=-0.35,padj=1.5)
mtext("Years 4-5",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_Biomass_2015_2019),]$AGR_Biomass_2015_2019,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.1,60),log="y")
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"))
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,60,10),labels=FALSE,lwd=0.25)
ct <- summary(lme_ROPS_Biomass_AGR_20156789)[[20]] # Coefficient table
ROPSmean <- mean(dat[dat$Species=="ROPS" & !is.na(dat$Biomass_est_kg_07),]$Biomass_est_kg_07,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ROPSmean,"red",17,TRUE)
BENImean <- mean(dat[dat$Species=="BENI" & !is.na(dat$Biomass_est_kg_07),]$Biomass_est_kg_07,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_BENI_Biomass_AGR_20156789)[[20]],
	BENImean,"blue",16,TRUE)
plot_siglets(-0.1,c(1,5,5,5),c("a","b","b","ab"),"blue")
mtext("e",at=-0.35,padj=1.5)
mtext("All years",at=1.5)
mtext("New York",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

##############################################
################### Oregon ###################
##############################################

source("OR_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_04 == 0,]$AGR_Biomass_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_Biomass_2017_2018 <- NA
dat[dat$Use_growth_05 == 0,]$AGR_Biomass_2016_2019 <- NA
dat[dat$Use_growth_05 == 0,]$AGR_Biomass_2018_2019 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_Biomass_2016_2020 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_Biomass_2019_2020 <- NA

dat[dat$Use_growth_04 == 0,]$Biomass_est_kg_04 <- NA
dat[dat$Use_growth_05 == 0,]$Biomass_est_kg_05 <- NA
dat[dat$Use_growth_06 == 0,]$Biomass_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2017),]$AGR_Biomass_2016_2017,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2016_2017,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2016_2017,na.rm=TRUE)),
	log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.1,0.3,1),labels=c("0.1","0.3","1"))
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,3,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2016_2017),]$AGR_Biomass_2016_2017,
	col="blue",pch=1,cex=0.8)	
mtext("f",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2017_2018),]$AGR_Biomass_2017_2018,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2017_2018,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2017_2018,na.rm=TRUE)),
	log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(1,2,5),labels=c("1","2","5"))
axis(side=2, at=seq(0.3,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2017_2018),]$AGR_Biomass_2017_2018,
	col="blue",pch=1,cex=0.8)	
mtext("g",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2018_2019),]$AGR_Biomass_2018_2019,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2018_2019,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2018_2019,na.rm=TRUE)),
	log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(1,3,10),labels=c("1","3","10"))
axis(side=2, at=seq(0.5,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(20),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019),]$AGR_Biomass_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("h",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2019_2020),]$AGR_Biomass_2019_2020,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.9*min(dat$AGR_Biomass_2019_2020,na.rm=TRUE),
	1.1*max(dat$AGR_Biomass_2019_2020,na.rm=TRUE)),
	log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(1,3,10),labels=c("1","3","10"))
axis(side=2, at=seq(0.5,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(20),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_Biomass_2018_2019),]$AGR_Biomass_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("i",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_Biomass_2016_2020),]$AGR_Biomass_2016_2020,
	col="white",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(1,6),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(1,2,3,4),labels=TRUE)
axis(side=2, at=seq(1,5,1),labels=FALSE,lwd=0.25)
ct <- summary(lme_ALRU_Biomass_AGR_20167890)[[20]] # Coefficient table
ALRUmean <- mean(dat[dat$Species=="ALRU" & !is.na(dat$Biomass_est_kg_04),]$Biomass_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ALRUmean,"orange",15,TRUE)
PSMEmean <- mean(dat[dat$Species=="PSME" & !is.na(dat$Biomass_est_kg_04),]$Biomass_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSME_Biomass_AGR_20167890)[[20]],
	PSMEmean,"blue",16,TRUE)
mtext("Oregon",side=4,padj=0.5)
mtext("j",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")


mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Biomass ~ AGR ~ (kg ~ y^-1)),
	side=2, outer=T, at=0.5, padj=0)


dev.off()

#############################################################
#############################################################
#############################################################