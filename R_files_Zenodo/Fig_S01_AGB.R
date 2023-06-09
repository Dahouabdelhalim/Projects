# NSF FX data plotting

# Figure S01
# Aboveground Biomass

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS01.pdf",
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

# Filter out size metrics from dead/unhealthy trees

dat[dat$Use_growth_01 == 0,]$AGB_est_kg_01 <- NA
dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_07 == 0,]$AGB_est_kg_07 <- NA
dat[dat$Use_growth_08 == 0,]$AGB_est_kg_08 <- NA

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_01),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_01),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01,
	col="blue",pch=1,cex=0.8)	
mtext("a",at=-0.35,padj=1.5)
mtext("Year 1",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_02),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_02),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="blue",pch=1,cex=0.8)	
mtext("b",at=-0.35,padj=1.5)
mtext("Year 2",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_04),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_04),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="blue",pch=1,cex=0.8)	
mtext("c",at=-0.35,padj=1.5)
mtext("Year 3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_07),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_07),]$AGB_est_kg_07,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_07),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_07),]$AGB_est_kg_07,
	col="blue",pch=1,cex=0.8)
mtext("d",at=-0.35,padj=1.5)
mtext("Year 4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_08),]$TRT+0.1),
	dat[dat$Species=="ROPS" & !is.na(dat$AGB_est_kg_08),]$AGB_est_kg_08,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_08),]$TRT-0.1),
	dat[dat$Species=="BENI" & !is.na(dat$AGB_est_kg_08),]$AGB_est_kg_08,
	col="blue",pch=1,cex=0.8)
mtext("e",at=-0.35,padj=1.5)
mtext("Year 5",at=1.5)
mtext("New York",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

##############################################
################### Oregon ###################
##############################################

source("OR_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out size metrics from dead/unhealthy trees

dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_05 == 0,]$AGB_est_kg_05 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_01),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_01),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01,
	col="blue",pch=1,cex=0.8)	
mtext("f",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_03),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_03),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="blue",pch=1,cex=0.8)	
mtext("g",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_04),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_04),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="blue",pch=1,cex=0.8)	
mtext("h",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_05),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_05),]$AGB_est_kg_05,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_05),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_05),]$AGB_est_kg_05,
	col="blue",pch=1,cex=0.8)	
mtext("i",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_06),]$TRT+0.1),
	dat[dat$Species=="ALRU" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="orange",pch=0,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_01),]$AGB_est_kg_01)*0.95,
	100),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,1,100),labels=c("0.01","1","100"))
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_06),]$TRT-0.1),
	dat[dat$Species=="PSME" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="blue",pch=1,cex=0.8)	
mtext("Oregon",side=4,padj=0.5)
mtext("j",at=-0.4,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Waiakea ###################
###############################################

source("HI_W_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out size metrics from dead/unhealthy trees

dat[dat$Use_growth_02 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_02),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_02),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_02),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="orange",pch=0,cex=0.8)	
mtext("k",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_03),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_03),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_03),]$TRT),
	dat[dat$Species=="CAEQ" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="orange",pch=0,cex=0.8)	
mtext("l",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_04),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_04),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_04),]$TRT),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="orange",pch=0,cex=0.8)	
mtext("m",at=-0.3,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_06),]$TRT+0.1),
	dat[dat$Species=="GLSE" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_06),]$TRT-0.1),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_06),]$TRT),
	dat[dat$Species=="PSCA" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="orange",pch=0,cex=0.8)	
mtext("n",at=-0.3,padj=1.5)
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

mtext("Waiakea",side=4,padj=0.5)

###############################################
################### Volcano ###################
###############################################

source("HI_V_growth_analysis.R")
source("Fig_3-6_setup.R")

# Filter out size metrics from dead/unhealthy trees

dat[dat$Use_growth_02 == 0 & dat$Use_growth_01 == 0,]$AGB_est_kg_02 <- NA
dat[dat$Use_growth_03 == 0,]$AGB_est_kg_03 <- NA
dat[dat$Use_growth_04 == 0,]$AGB_est_kg_04 <- NA
dat[dat$Use_growth_06 == 0,]$AGB_est_kg_06 <- NA

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_02),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_02),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_02),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02,
	col="orange",pch=0,cex=0.8)	
mtext("o",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_03),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_03),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_03),]$TRT),
	dat[dat$Species=="MOFA" & !is.na(dat$AGB_est_kg_03),]$AGB_est_kg_03,
	col="orange",pch=0,cex=0.8)	
mtext("p",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_04),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_04),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_04),]$TRT),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_04),]$AGB_est_kg_04,
	col="orange",pch=0,cex=0.8)	
mtext("q",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot(jitter(dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_06),]$TRT+0.1),
	dat[dat$Species=="ACKO" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="red",pch=2,cex=0.8,xlim=c(-0.5,3.5),xaxt="n",yaxt="n",
	ylim=c(min(dat[!is.na(dat$AGB_est_kg_02),]$AGB_est_kg_02)*0.95,
	100),log="y")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","   15+P"), padj=-.5)
axis(side=2, at=c(0.01,100),labels=c("0.01","100"))
axis(side=2, at=c(1),labels=c("1"))
axis(side=2, at=seq(0.0001,0.0009,0.0001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.001,0.009,0.001),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.02,0.1,0.01),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.2,0.9,0.1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(2,10,1),labels=FALSE,lwd=0.25)
axis(side=2, at=seq(20,90,10),labels=FALSE,lwd=0.25)
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_06),]$TRT-0.1),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="blue",pch=1,cex=0.8)	
points(jitter(dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_06),]$TRT),
	dat[dat$Species=="DOVI" & !is.na(dat$AGB_est_kg_06),]$AGB_est_kg_06,
	col="orange",pch=0,cex=0.8)	
mtext("r",at=-0.35,padj=1.5)
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

mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Aboveground ~ biomass ~ (kg)),
	side=2, outer=T, at=0.5, padj=0)
mtext("Volcano",side=4,padj=0.5)


dev.off()

#############################################################
#############################################################
#############################################################