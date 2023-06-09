# NSF FX data plotting

# Figure S4
# Total Biomass relative growth rate

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS04.pdf",
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

dat[dat$Use_growth_02 == 0,]$RGR_Biomass_2015_2016 <- NA
dat[dat$Use_growth_04 == 0,]$RGR_Biomass_2015_2017 <- NA
dat[dat$Use_growth_04 == 0,]$RGR_Biomass_2016_2017 <- NA
dat[dat$Use_growth_07 == 0,]$RGR_Biomass_2015_2018 <- NA
dat[dat$Use_growth_07 == 0,]$RGR_Biomass_2017_2018 <- NA
dat[dat$Use_growth_08 == 0,]$RGR_Biomass_2015_2019 <- NA
dat[dat$Use_growth_08 == 0,]$RGR_Biomass_2018_2019 <- NA

dat[dat$Use_growth_01 == 0,]$Biomass_est_kg_01 <- NA
dat[dat$Use_growth_02 == 0,]$Biomass_est_kg_02 <- NA
dat[dat$Use_growth_04 == 0,]$Biomass_est_kg_04 <- NA
dat[dat$Use_growth_07 == 0,]$Biomass_est_kg_07 <- NA
dat[dat$Use_growth_08 == 0,]$Biomass_est_kg_08 <- NA

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2016),]$RGR_Biomass_2015_2016,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2015_2016,na.rm=TRUE)))
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2015_2016),]$RGR_Biomass_2015_2016,
	col="blue",pch=1,cex=0.8)	
mtext("a",at=-0.35,padj=1.5)
mtext("Years 1-2",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2016_2017),]$RGR_Biomass_2016_2017,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2016_2017,na.rm=TRUE)))
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2016_2017),]$RGR_Biomass_2016_2017,
	col="blue",pch=1,cex=0.8)	
axis(side=2, at=0:3, labels=c("0","1","2","3"))
axis(side=2, at=seq(0.5,2.5,0.5), labels=FALSE)
mtext("b",at=-0.35,padj=1.5)
mtext("Years 2-3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2017_2018),]$RGR_Biomass_2017_2018,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2017_2018,na.rm=TRUE)))
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2017_2018),]$RGR_Biomass_2017_2018,
	col="blue",pch=1,cex=0.8)	
axis(side=2, at=0:2, labels=c("0","1","2"))
axis(side=2, at=seq(0.5,1.5,0.5), labels=FALSE)
mtext("c",at=-0.35,padj=1.5)
mtext("Years 3-4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2018_2019),]$RGR_Biomass_2018_2019,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2018_2019,na.rm=TRUE)))
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$RGR_Biomass_2018_2019),]$RGR_Biomass_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("d",at=-0.35,padj=1.5)
mtext("Years 4-5",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$RGR_Biomass_2015_2019),]$RGR_Biomass_2015_2019,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,3))
axis(side=2, at=0:3, labels=c("0","1","2","3"))
axis(side=2, at=seq(0.5,2.5,0.5), labels=FALSE)
ct <- summary(lme_ROPS_Biomass_RGR_20156789)[[20]] # Coefficient table
ROPSmean <- mean(dat[dat$Species=="ROPS" & !is.na(dat$Biomass_est_kg_07),]$Biomass_est_kg_07,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ROPSmean,"red",17,FALSE)
BENImean <- mean(dat[dat$Species=="BENI" & !is.na(dat$Biomass_est_kg_07),]$Biomass_est_kg_07,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_BENI_Biomass_RGR_20156789)[[20]],
	BENImean,"blue",16,FALSE)
plot_siglets(-0.1,2.3,c("a","b","b","b"),"blue")
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

dat[dat$Use_growth_04 == 0,]$RGR_Biomass_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$RGR_Biomass_2017_2018 <- NA
dat[dat$Use_growth_05 == 0,]$RGR_Biomass_2016_2019 <- NA
dat[dat$Use_growth_05 == 0,]$RGR_Biomass_2018_2019 <- NA
dat[dat$Use_growth_06 == 0,]$RGR_Biomass_2016_2020 <- NA
dat[dat$Use_growth_06 == 0,]$RGR_Biomass_2019_2020 <- NA

dat[dat$Use_growth_04 == 0,]$Biomass_est_kg_04 <- NA
dat[dat$Use_growth_05 == 0,]$Biomass_est_kg_05 <- NA
dat[dat$Use_growth_06 == 0,]$Biomass_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2017),]$RGR_Biomass_2016_2017,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2016_2017,na.rm=TRUE)))
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2016_2017),]$RGR_Biomass_2016_2017,
	col="blue",pch=1,cex=0.8)	
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
mtext("f",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2017_2018),]$RGR_Biomass_2017_2018,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2017_2018,na.rm=TRUE)))
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2017_2018),]$RGR_Biomass_2017_2018,
	col="blue",pch=1,cex=0.8)	
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=0:1, labels=c("0","1"))
axis(side=2, at=seq(0.5,1.5,0.5), labels=FALSE)
mtext("g",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2018_2019),]$RGR_Biomass_2018_2019,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2018_2019,na.rm=TRUE)))
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2018_2019),]$RGR_Biomass_2018_2019,
	col="blue",pch=1,cex=0.8)	
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0,0.4,0.8), labels=c("0","0.4","0.8"))
axis(side=2, at=seq(0,1.2,0.2), labels=FALSE)
mtext("h",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2019_2020),]$RGR_Biomass_2019_2020,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,max(dat$RGR_Biomass_2019_2020,na.rm=TRUE)))
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$RGR_Biomass_2019_2020),]$RGR_Biomass_2019_2020,
	col="blue",pch=1,cex=0.8)	
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0,0.4,0.8,1.2,1.6), labels=c("0","0.4","0.8","1.2","1.6"))
axis(side=2, at=seq(0,1.8,0.2), labels=FALSE)
mtext("i",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$RGR_Biomass_2016_2020),]$RGR_Biomass_2016_2020,
	col="white",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0,2.5))
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=0:2, labels=c("0","1","2"))
axis(side=2, at=seq(0.5,1.5,0.5), labels=FALSE)
ct <- summary(lme_ALRU_Biomass_RGR_20167890)[[20]] # Coefficient table
ALRUmean <- mean(dat[dat$Species=="ALRU" & !is.na(dat$Biomass_est_kg_04),]$Biomass_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ALRUmean,"orange",15,FALSE)
PSMEmean <- mean(dat[dat$Species=="PSME" & !is.na(dat$Biomass_est_kg_04),]$Biomass_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSME_Biomass_RGR_20167890)[[20]],
	PSMEmean,"blue",16,FALSE)
mtext("Oregon",side=4,padj=0.5)
mtext("j",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Biomass ~ RGR ~ (y^-1)),
	side=2, outer=T, at=0.4, padj=0)

dev.off()

#############################################################
#############################################################
#############################################################