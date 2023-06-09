# NSF FX data plotting

# Figure 3
# Aboveground Biomass absolute growth rate

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="Fig03.pdf",
	width=7.05,height=5)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,5),omi=c(.4,.25,.2,.2),mai=c(0,.3,0,0))

##############################################
################## New York ##################
##############################################

# Run code to get the data to plot
source("NY_growth_analysis.R")
# Run code for plotting specifics
source("Fig_3-6_setup.R")

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_02 == 0,]$AGR_AGB_2015_2016 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2015_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_AGB_2015_2018 <- NA
dat[dat$Use_growth_07 == 0,]$AGR_AGB_2017_2018 <- NA
dat[dat$Use_growth_08 == 0,]$AGR_AGB_2015_2019 <- NA
dat[dat$Use_growth_08 == 0,]$AGR_AGB_2018_2019 <- NA

dat[dat$Use_growth_01 == 0,]$AGB_est_kg_01 <- NA
dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_07 == 0,]$AGB_est_kg_07 <- NA
dat[dat$Use_growth_08 == 0,]$AGB_est_kg_08 <- NA

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2016),]$AGR_AGB_2015_2016,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.002,2),log="y")
axis(side=2, at=c(0.01,0.1,1),labels=c("0.01","0.1","1"))
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(2),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2015_2016),]$AGR_AGB_2015_2016,
	col="blue",pch=1,cex=0.8)	
mtext("a",at=-0.35,padj=1.5)
mtext("Years 1-2",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.005,15),log="y")
axis(side=2, at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.005,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="blue",pch=1,cex=0.8)	
mtext("b",at=-0.35,padj=1.5)
mtext("Years 2-3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,50),log="y")
axis(side=2, at=c(0.01,1),labels=c("0.01","1"))
axis(side=2, at=c(0.1,10),labels=c("0.1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="blue",pch=1,cex=0.8)	
mtext("c",at=-0.35,padj=1.5)
mtext("Years 3-4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,50),log="y")
axis(side=2, at=c(0.01,1),labels=c("0.01","1"))
axis(side=2, at=c(0.1,10),labels=c("0.1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2018_2019),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("d",at=-0.35,padj=1.5)
mtext("Years 4-5",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2019),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGR_AGB_2015_2019),]$AGR_AGB_2015_2019,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.1,50),log="y")
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"))
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,40,10),labels=FALSE,lwd=0.25)
ct <- summary(lme_ROPS_AGB_AGR_20156789)[[20]] # Coefficient table
ROPSmean <- mean(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_07),]$AGB_est_kg_07,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ROPSmean,"red",17,TRUE)
BENImean <- mean(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_07),]$AGB_est_kg_07,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_BENI_AGB_AGR_20156789)[[20]],
	BENImean,"blue",16,TRUE)
plot_siglets(-0.1,c(1,4,4,4),c("a","b","b","ab"),"blue")
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

dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA
dat[dat$Use_growth_05 == 0,]$AGR_AGB_2016_2019 <- NA
dat[dat$Use_growth_05 == 0,]$AGR_AGB_2018_2019 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_AGB_2019_2020 <- NA

dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_05 == 0,]$AGB_est_kg_05 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.05,3),log="y")
axis(side=2, at=c(0.1,0.3,1),labels=c("0.1","0.3","1"))
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,3,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="blue",pch=1,cex=0.8)	
mtext("f",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.3,10),log="y")
axis(side=2, at=c(1,2,5),labels=c("1","2","5"))
axis(side=2, at=seq(0.3,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="blue",pch=1,cex=0.8)	
mtext("g",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.5,20),log="y")
axis(side=2, at=c(1,3,10),labels=c("1","3","10"))
axis(side=2, at=seq(0.5,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(20),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="blue",pch=1,cex=0.8)	
mtext("h",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2019_2020),]$AGR_AGB_2019_2020,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(1,30),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(1,3,10),labels=c("1","3","10"))
axis(side=2, at=seq(0.5,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(20),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGR_AGB_2019_2020),]$AGR_AGB_2019_2020,
	col="blue",pch=1,cex=0.8)	
mtext("i",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2020),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGR_AGB_2016_2020),]$AGR_AGB_2016_2020,
	col="white",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(1,5),log="y")
axis(side=2, at=c(1,2,3,4),labels=TRUE)
axis(side=2, at=seq(1,5,1),labels=FALSE,lwd=0.25)
ct <- summary(lme_ALRU_AGB_AGR_20167890)[[20]] # Coefficient table
ALRUmean <- mean(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ALRUmean,"orange",15,TRUE)
PSMEmean <- mean(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSME_AGB_AGR_20167890)[[20]],
	PSMEmean,"blue",16,TRUE)
mtext("Oregon",side=4,padj=0.5)
mtext("j",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Waiakea ###################
###############################################

source("HI_W_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_03 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_AGB_2016_2019 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_AGB_2018_2019 <- NA

dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.002,5),log="y")
axis(side=2, at=c(0.01,0.1,1),labels=c("0.01","0.1","1"))
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,5,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="orange",pch=0,cex=0.8)	
mtext("k",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.003,14),log="y")
axis(side=2, at=c(0.01,1),labels=c("0.01","1"))
axis(side=2, at=c(0.1,10),labels=c("0.1","10"))
axis(side=2, at=seq(0.003,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="orange",pch=0,cex=0.8)	
mtext("l",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2018_2019),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.04,25),log="y")
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"))
axis(side=2, at=seq(0.04,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=c(20),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2018_2019),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2018_2019),]$TRT),
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="orange",pch=0,cex=0.8)	
mtext("m",at=-0.3,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()
points(c(0,0),c(.6,.45),pch=c(2,1),col=c("red","blue"),cex=0.8,xpd=TRUE)
text(0,0.525,"a-e",pos=2,font=1,xpd=TRUE)
text(0,0.6,"Robinia",pos=4,col="red",font=3,xpd=TRUE)
text(0,0.45,"Betula",pos=4,col="blue",font=3,xpd=TRUE)

points(c(0,0),c(.2,.05),pch=c(0,1),col=c("orange","blue"),cex=0.8,xpd=TRUE)
text(0,0.125,"f-j",pos=2,font=1,xpd=TRUE)
text(0,0.2,"Alnus",pos=4,col="orange",font=3,xpd=TRUE)
text(0,0.05,"Pseudotsuga",pos=4,col="blue",font=3,xpd=TRUE)

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2019),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2019),]$AGR_AGB_2016_2019,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.1,10),log="y")
axis(side=2, at=c(0.3,3),labels=c("0.3","3"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.1,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
ct <- summary(lme_GLSE_AGB_AGR_2016789)[[20]] # Coefficient table
GLSEmean <- mean(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,GLSEmean,"red",17,TRUE)
PSCAmean <- mean(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_PSCA_AGB_AGR_2016789)[[20]],
	PSCAmean,"blue",16,TRUE)
CAEQmean <- mean(dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0,summary(lme_CAEQ_AGB_AGR_2016789)[[20]],
	CAEQmean,"orange",15,TRUE)
mtext("n",at=-0.25,padj=1.5)
mtext("Waiakea",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Volcano ###################
###############################################

source("HI_V_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out growth and size metrics from dead/unhealthy trees

dat[dat$Use_growth_03 == 0,]$AGR_AGB_2016_2017 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2016_2018 <- NA
dat[dat$Use_growth_04 == 0,]$AGR_AGB_2017_2018 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_AGB_2016_2019 <- NA
dat[dat$Use_growth_06 == 0,]$AGR_AGB_2018_2019 <- NA

dat[dat$Use_growth_02 == 0 & dat$Use_growth_01 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.001,10),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.001,0.1,1,10),labels=c("0.001","0.1","1","10"))
axis(side=2, at=seq(0.002,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017),]$AGR_AGB_2016_2017,
	col="orange",pch=0,cex=0.8)	
mtext("o",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,30),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,0.1,1,10),labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,30,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018),]$AGR_AGB_2017_2018,
	col="orange",pch=0,cex=0.8)	
mtext("p",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2018_2019),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.01,60),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,1),labels=c("0.01","1"))
axis(side=2, at=c(0.1,10),labels=c("0.1","10"))
axis(side=2, at=seq(0.02,0.09,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,60,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2018_2019),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2018_2019),]$TRT),
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2018_2019),]$AGR_AGB_2018_2019,
	col="orange",pch=0,cex=0.8)	
mtext("q",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()
points(c(0,0,0),c(.9,.75,.6),pch=c(2,0,1),col=c("red","orange","blue"),cex=0.8,xpd=TRUE)
text(0,0.75,"k-n",pos=2,font=1,xpd=TRUE)
text(0,0.9,"Gliricidia",pos=4,col="red",font=3,xpd=TRUE)
text(0,0.75,"Casuarina",pos=4,col="orange",font=3,xpd=TRUE)
text(0,0.6,"Psidium",pos=4,col="blue",font=3,xpd=TRUE)

points(c(0,0,0),c(.35,.2,.05),pch=c(2,0,1),col=c("red","orange","blue"),cex=0.8,xpd=TRUE)
text(0,0.2,"o-r",pos=2,font=1,xpd=TRUE)
text(0,0.35,"Acacia",pos=4,col="red",font=3,xpd=TRUE)
text(0,0.2,"Morella",pos=4,col="orange",font=3,xpd=TRUE)
text(0,0.05,"Dodonaea",pos=4,col="blue",font=3,xpd=TRUE)

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2019),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2019),]$AGR_AGB_2016_2019,
	col="white",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(0.1,10),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.1,1,10),labels=c("0.1","1","10"))
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,9,1),labels=FALSE,lwd=0.25)
ct <- summary(lme_ACKO_AGB_AGR_2016789)[[20]] # Coefficient table
ACKOmean <- mean(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0.1,ct,ACKOmean,"red",17,TRUE)
DOVImean <- mean(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(-0.1,summary(lme_DOVI_AGB_AGR_2016789)[[20]],
	DOVImean,"blue",16,TRUE)
MOFAmean <- mean(dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,na.rm=TRUE)
plot_u_se_lme(0,summary(lme_MOFA_AGB_AGR_2016789)[[20]],
	MOFAmean,"orange",15,TRUE)
mtext("r",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Absolute ~ growth ~ rate ~ of ~ AGB ~ (kg ~ y^-1)),
	side=2, outer=T, at=0.5, padj=0)
mtext("Volcano",side=4,padj=0.5)


dev.off()

#############################################################
#############################################################
#############################################################