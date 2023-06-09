# NSF FX data plotting

# Figure 2
# Soil phosphatase by treatment and species

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

library(colorBlindness)

#########################################################################
############################# Set up figure #############################
#########################################################################

#pdf(file="Pase_Fig02.pdf",
#	width=6,height=4)
setEPS()
postscript(file="Pase_Fig02.eps",
           width=6,height=4)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(2,3),omi=c(.4,.25,.2,.2),mai=c(.2,.3,0,0))

N <- read.csv("NY_FX_Soil_Pase.csv")
O <- read.csv("OR_FX_Soil_Pase.csv")
W <- read.csv("HI_W_FX_Soil_Pase.csv")
V <- read.csv("HI_V_FX_Soil_Pase.csv")

N_TRT <- rep(-1,nrow(N))
N_TRT[N$Treatment=="LN"] <- 0
N_TRT[N$Treatment=="MN"] <- 1
N_TRT[N$Treatment=="HN"] <- 2
N_TRT[N$Treatment=="PHN"] <- 3
N_jTRT <- jitter(N_TRT)

O_TRT <- rep(-1,nrow(O))
O_TRT[O$Treatment=="LN"] <- 0
O_TRT[O$Treatment=="MN"] <- 1
O_TRT[O$Treatment=="HN"] <- 2
O_TRT[O$Treatment=="PHN"] <- 3
O_jTRT <- jitter(O_TRT)

W_TRT <- rep(-1,nrow(W))
W_TRT[W$Treatment=="LN"] <- 0
W_TRT[W$Treatment=="MN"] <- 1
W_TRT[W$Treatment=="HN"] <- 2
W_TRT[W$Treatment=="PHN"] <- 3
W_jTRT <- jitter(W_TRT,amount=0.1)

V_TRT <- rep(-1,nrow(V))
V_TRT[V$Treatment=="LN"] <- 0
V_TRT[V$Treatment=="MN"] <- 1
V_TRT[V$Treatment=="HN"] <- 2
V_TRT[V$Treatment=="PHN"] <- 3
V_jTRT <- jitter(V_TRT,amount=0.1)

nonfixcol <- Blue2DarkRed12Steps[1]
rhizcol <- Blue2DarkRed12Steps[11]
actcol <- Blue2DarkRed12Steps[9]

##############################################
################## New York ##################
##############################################

ylx <- max(c(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04),na.rm=TRUE)*1.16
plot(N_jTRT[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI"]-0.25,
	N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI",]$Soil_Pase_umol_g_hr_04,
	xlim=c(-0.5,3.5),ylim=c(0,ylx),
	pch=1,xaxt="n",col=nonfixcol,cex=1.5,cex.axis=0.8)
abline(h=0,col="gray")
points(N_jTRT[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS"]+0.25,
	N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS",]$Soil_Pase_umol_g_hr_04,
	col=rhizcol,pch=2,cex=1.5)

BENI_04_LN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
BENI_04_LN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(-0.25,BENI_04_LN_u,pch=16,col=nonfixcol,cex=2)
segments(-0.25,BENI_04_LN_u - BENI_04_LN_se,-0.25,BENI_04_LN_u + BENI_04_LN_se,col=nonfixcol)
BENI_04_MN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
BENI_04_MN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(0.75,BENI_04_MN_u,pch=16,col=nonfixcol,cex=2)
segments(0.75,BENI_04_MN_u - BENI_04_MN_se,0.75,BENI_04_MN_u + BENI_04_MN_se,col=nonfixcol)
BENI_04_HN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
BENI_04_HN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(1.75,BENI_04_HN_u,pch=16,col=nonfixcol,cex=2)
segments(1.75,BENI_04_HN_u - BENI_04_HN_se,1.75,BENI_04_HN_u + BENI_04_HN_se,col=nonfixcol)
BENI_04_PHN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
BENI_04_PHN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(2.75,BENI_04_PHN_u,pch=16,col=nonfixcol,cex=2)
segments(2.75,BENI_04_PHN_u - BENI_04_PHN_se,2.75,BENI_04_PHN_u + BENI_04_PHN_se,col=nonfixcol)

ROPS_04_LN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ROPS_04_LN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0.25,ROPS_04_LN_u,pch=17,col=rhizcol,cex=2)
segments(0.25,ROPS_04_LN_u - ROPS_04_LN_se,0.25,ROPS_04_LN_u + ROPS_04_LN_se,col=rhizcol)
ROPS_04_MN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ROPS_04_MN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1.25,ROPS_04_MN_u,pch=17,col=rhizcol,cex=2)
segments(1.25,ROPS_04_MN_u - ROPS_04_MN_se,1.25,ROPS_04_MN_u + ROPS_04_MN_se,col=rhizcol)
ROPS_04_HN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ROPS_04_HN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2.25,ROPS_04_HN_u,pch=17,col=rhizcol,cex=2)
segments(2.25,ROPS_04_HN_u - ROPS_04_HN_se,2.25,ROPS_04_HN_u + ROPS_04_HN_se,col=rhizcol)
ROPS_04_PHN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ROPS_04_PHN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3.25,ROPS_04_PHN_u,pch=17,col=rhizcol,cex=2)
segments(3.25,ROPS_04_PHN_u - ROPS_04_PHN_se,3.25,ROPS_04_PHN_u + ROPS_04_PHN_se,col=rhizcol)

mtext("a",at=-0.35,padj=1.5,cex=0.8)
mtext("New York 2017",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

summary(NY_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
	data=N[N$Use_growth_04==1,]))
anova(NY_FX_Pase_04_lm)

# a a a a a a a a
#text(c(-0.2,0.2,0.8,1.2,1.8,2.2,2.8,3.2),40000,"A",cex=1.25)

ylx <- max(c(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1,]$Soil_Pase_umol_g_hr_07),na.rm=TRUE)*1.06
plot(N_jTRT[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI"]-0.25,
	N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI",]$Soil_Pase_umol_g_hr_07,
	xlim=c(-0.5,3.5),ylim=c(0,ylx),
	pch=1,xaxt="n",col=nonfixcol,cex=1.5,cex.axis=0.8)
abline(h=0,col="gray")
points(N_jTRT[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS"]+0.25,
	N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS",]$Soil_Pase_umol_g_hr_07,
	col=rhizcol,pch=2,cex=1.5)

BENI_07_LN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
BENI_07_LN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07))
points(-0.25,BENI_07_LN_u,pch=16,col=nonfixcol,cex=2)
segments(-0.25,BENI_07_LN_u - BENI_07_LN_se,-0.25,BENI_07_LN_u + BENI_07_LN_se,col=nonfixcol)
BENI_07_MN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
BENI_07_MN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07))
points(0.75,BENI_07_MN_u,pch=16,col=nonfixcol,cex=2)
segments(0.75,BENI_07_MN_u - BENI_07_MN_se,0.75,BENI_07_MN_u + BENI_07_MN_se,col=nonfixcol)
BENI_07_HN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
BENI_07_HN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07))
points(1.75,BENI_07_HN_u,pch=16,col=nonfixcol,cex=2)
segments(1.75,BENI_07_HN_u - BENI_07_HN_se,1.75,BENI_07_HN_u + BENI_07_HN_se,col=nonfixcol)
BENI_07_PHN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
BENI_07_PHN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07))
points(2.75,BENI_07_PHN_u,pch=16,col=nonfixcol,cex=2)
segments(2.75,BENI_07_PHN_u - BENI_07_PHN_se,2.75,BENI_07_PHN_u + BENI_07_PHN_se,col=nonfixcol)

ROPS_07_LN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
ROPS_07_LN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07))
points(0.25,ROPS_07_LN_u,pch=17,col=rhizcol,cex=2)
segments(0.25,ROPS_07_LN_u - ROPS_07_LN_se,0.25,ROPS_07_LN_u + ROPS_07_LN_se,col=rhizcol)
ROPS_07_MN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
ROPS_07_MN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07))
points(1.25,ROPS_07_MN_u,pch=17,col=rhizcol,cex=2)
segments(1.25,ROPS_07_MN_u - ROPS_07_MN_se,1.25,ROPS_07_MN_u + ROPS_07_MN_se,col=rhizcol)
ROPS_07_HN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
ROPS_07_HN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07))
points(2.25,ROPS_07_HN_u,pch=17,col=rhizcol,cex=2)
segments(2.25,ROPS_07_HN_u - ROPS_07_HN_se,2.25,ROPS_07_HN_u + ROPS_07_HN_se,col=rhizcol)
ROPS_07_PHN_u <- mean(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)
ROPS_07_PHN_se <- sd(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)/
  sqrt(length(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07))
points(3.25,ROPS_07_PHN_u,pch=17,col=rhizcol,cex=2)
segments(3.25,ROPS_07_PHN_u - ROPS_07_PHN_se,3.25,ROPS_07_PHN_u + ROPS_07_PHN_se,col=rhizcol)

mtext("b",at=-0.35,padj=1.5,cex=0.8)
mtext("New York 2018",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5,cex.axis=0.8)

summary(NY_FX_Pase_07_lm <- lm(Soil_Pase_umol_g_hr_07 ~ Species * Treatment,
	data=N[N$Use_growth_07==1,]))
anova(NY_FX_Pase_07_lm)

# a a a a a a a a
#text(c(-0.2,0.2,0.8,1.2,1.8,2.2,2.8,3.2),23500,"A",cex=1.25)

##############################################
################### Oregon ###################
##############################################

ylx <- max(c(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04),na.rm=TRUE)*1.06
plot(O_jTRT[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME"]-0.25,
	O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME",]$Soil_Pase_umol_g_hr_04,
	xlim=c(-0.5,3.5),ylim=c(0,ylx),
	pch=1,xaxt="n",col=nonfixcol,cex=1.5,cex.axis=0.8)
abline(h=0,col="gray")
points(O_jTRT[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU"]+0.25,
	O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU",]$Soil_Pase_umol_g_hr_04,
	col=actcol,pch=0,cex=1.5)

PSME_04_LN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSME_04_LN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(-0.25,PSME_04_LN_u,pch=16,col=nonfixcol,cex=2)
segments(-0.25,PSME_04_LN_u - PSME_04_LN_se,-0.25,PSME_04_LN_u + PSME_04_LN_se,col=nonfixcol)
PSME_04_MN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSME_04_MN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(0.75,PSME_04_MN_u,pch=16,col=nonfixcol,cex=2)
segments(0.75,PSME_04_MN_u - PSME_04_MN_se,0.75,PSME_04_MN_u + PSME_04_MN_se,col=nonfixcol)
PSME_04_HN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSME_04_HN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(1.75,PSME_04_HN_u,pch=16,col=nonfixcol,cex=2)
segments(1.75,PSME_04_HN_u - PSME_04_HN_se,1.75,PSME_04_HN_u + PSME_04_HN_se,col=nonfixcol)
PSME_04_PHN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSME_04_PHN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(2.75,PSME_04_PHN_u,pch=16,col=nonfixcol,cex=2)
segments(2.75,PSME_04_PHN_u - PSME_04_PHN_se,2.75,PSME_04_PHN_u + PSME_04_PHN_se,col=nonfixcol)

ALRU_04_LN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ALRU_04_LN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0.25,ALRU_04_LN_u,pch=15,col=actcol,cex=2)
segments(0.25,ALRU_04_LN_u - ALRU_04_LN_se,0.25,ALRU_04_LN_u + ALRU_04_LN_se,col=actcol)
ALRU_04_MN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ALRU_04_MN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1.25,ALRU_04_MN_u,pch=15,col=actcol,cex=2)
segments(1.25,ALRU_04_MN_u - ALRU_04_MN_se,1.25,ALRU_04_MN_u + ALRU_04_MN_se,col=actcol)
ALRU_04_HN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ALRU_04_HN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2.25,ALRU_04_HN_u,pch=15,col=actcol,cex=2)
segments(2.25,ALRU_04_HN_u - ALRU_04_HN_se,2.25,ALRU_04_HN_u + ALRU_04_HN_se,col=actcol)
ALRU_04_PHN_u <- mean(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ALRU_04_PHN_se <- sd(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3.25,ALRU_04_PHN_u,pch=15,col=actcol,cex=2)
segments(3.25,ALRU_04_PHN_u - ALRU_04_PHN_se,3.25,ALRU_04_PHN_u + ALRU_04_PHN_se,col=actcol)

mtext("c",at=-0.35,padj=1.5,cex=0.8)
mtext("Oregon 2018",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

summary(OR_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
	data=O[O$Use_growth_04==1,]))
anova(OR_FX_Pase_04_lm)
TukeyHSD(aov(OR_FX_Pase_04_lm))

# Treatment significant but not species or treatment x species
# LN > HN, MN, PHN
text(0:3,30,c("A","B","B","B"),cex=1)

###############################################
################### Waiakea ###################
###############################################

ylx <- max(c(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5,]$Soil_Pase_umol_g_hr_04),na.rm=TRUE)*1.22
plot(W_jTRT[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="PSCA"]-0.2,
	W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="PSCA",]$Soil_Pase_umol_g_hr_04,
	xlim=c(-0.5,3.5),ylim=c(0,ylx),
	pch=1,xaxt="n",col=nonfixcol,cex=1.5,cex.axis=0.8)
abline(h=0,col="gray")
points(W_jTRT[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE"],
	W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE",]$Soil_Pase_umol_g_hr_04,
	col=rhizcol,pch=2,cex=1.5)
points(W_jTRT[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$Soil_Pase_umol_g_hr_04,
	col=actcol,pch=0,cex=1.5)

PSCA_04_LN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSCA_04_LN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(-0.2,PSCA_04_LN_u,pch=16,col=nonfixcol,cex=2)
segments(-0.2,PSCA_04_LN_u - PSCA_04_LN_se,-0.2,PSCA_04_LN_u + PSCA_04_LN_se,col=nonfixcol)
PSCA_04_MN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSCA_04_MN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(0.8,PSCA_04_MN_u,pch=16,col=nonfixcol,cex=2)
segments(0.8,PSCA_04_MN_u - PSCA_04_MN_se,0.8,PSCA_04_MN_u + PSCA_04_MN_se,col=nonfixcol)
PSCA_04_HN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSCA_04_HN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(1.8,PSCA_04_HN_u,pch=16,col=nonfixcol,cex=2)
segments(1.8,PSCA_04_HN_u - PSCA_04_HN_se,1.8,PSCA_04_HN_u + PSCA_04_HN_se,col=nonfixcol)
PSCA_04_PHN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
PSCA_04_PHN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(2.8,PSCA_04_PHN_u,pch=16,col=nonfixcol,cex=2)
segments(2.8,PSCA_04_PHN_u - PSCA_04_PHN_se,2.8,PSCA_04_PHN_u + PSCA_04_PHN_se,col=nonfixcol)

GLSE_04_LN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
GLSE_04_LN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0,GLSE_04_LN_u,pch=17,col=rhizcol,cex=2)
segments(0,GLSE_04_LN_u - GLSE_04_LN_se,0,GLSE_04_LN_u + GLSE_04_LN_se,col=rhizcol)
GLSE_04_MN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
GLSE_04_MN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1,GLSE_04_MN_u,pch=17,col=rhizcol,cex=2)
segments(1,GLSE_04_MN_u - GLSE_04_MN_se,1,GLSE_04_MN_u + GLSE_04_MN_se,col=rhizcol)
GLSE_04_HN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
GLSE_04_HN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2,GLSE_04_HN_u,pch=17,col=rhizcol,cex=2)
segments(2,GLSE_04_HN_u - GLSE_04_HN_se,2,GLSE_04_HN_u + GLSE_04_HN_se,col=rhizcol)
GLSE_04_PHN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
GLSE_04_PHN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3,GLSE_04_PHN_u,pch=17,col=rhizcol,cex=2)
segments(3,GLSE_04_PHN_u - GLSE_04_PHN_se,3,GLSE_04_PHN_u + GLSE_04_PHN_se,col=rhizcol)

CAEQ_04_LN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
CAEQ_04_LN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0.2,CAEQ_04_LN_u,pch=15,col=actcol,cex=2)
segments(0.2,CAEQ_04_LN_u - CAEQ_04_LN_se,0.2,CAEQ_04_LN_u + CAEQ_04_LN_se,col=actcol)
CAEQ_04_MN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
CAEQ_04_MN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1.2,CAEQ_04_MN_u,pch=15,col=actcol,cex=2)
segments(1.2,CAEQ_04_MN_u - CAEQ_04_MN_se,1.2,CAEQ_04_MN_u + CAEQ_04_MN_se,col=actcol)
CAEQ_04_HN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
CAEQ_04_HN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2.2,CAEQ_04_HN_u,pch=15,col=actcol,cex=2)
segments(2.2,CAEQ_04_HN_u - CAEQ_04_HN_se,2.2,CAEQ_04_HN_u + CAEQ_04_HN_se,col=actcol)
CAEQ_04_PHN_u <- mean(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
CAEQ_04_PHN_se <- sd(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3.2,CAEQ_04_PHN_u,pch=15,col=actcol,cex=2)
segments(3.2,CAEQ_04_PHN_u - CAEQ_04_PHN_se,3.2,CAEQ_04_PHN_u + CAEQ_04_PHN_se,col=actcol)

mtext("d",at=-0.35,padj=1.5,cex=0.8)
mtext("Waiakea 2018",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5,cex.axis=0.8)

summary(HI_W_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
	data=W[W$Use_growth_04==1,]))
anova(HI_W_FX_Pase_04_lm)
TukeyHSD(aov(HI_W_FX_Pase_04_lm))

# Treatment significant but not species or treatment x species
# LN > MN, PHN
text(0:3,28,c("A","B","AB","B"),cex=1)

###############################################
################### Legend ####################
###############################################

plot.new()

points(c(0.3,0.3),c(.95,.85),pch=c(2,1),col=c(rhizcol,nonfixcol),cex=1,xpd=TRUE)
text(0.3,0.95,"Robinia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.3,0.85,"Betula",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.28,0.9,"a-b",pos=2,xpd=TRUE,cex=1)

points(c(0.3,0.3),c(.7,.6),pch=c(0,1),col=c(actcol,nonfixcol),cex=1,xpd=TRUE)
text(0.3,0.7,"Alnus",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.3,0.6,"Pseudotsuga",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.28,0.65,"c",pos=2,xpd=TRUE,cex=1)

points(c(0.3,0.3,0.3),c(.45,.35,.25),pch=c(2,0,1),col=c(rhizcol,actcol,nonfixcol),cex=1,xpd=TRUE)
text(0.3,0.45,"Gliricidia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.3,0.35,"Casuarina",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.3,0.25,"Psidium",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.28,0.35,"d",pos=2,xpd=TRUE,cex=1)

points(c(0.3,0.3,0.3),c(.1,.0,-.1),pch=c(2,0,1),col=c(rhizcol,actcol,nonfixcol),cex=1,xpd=TRUE)
text(0.3,0.1,"Acacia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.3,0.0,"Morella",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.3,-0.1,"Dodonaea",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.28,0.0,"e",pos=2,xpd=TRUE,cex=1)

###############################################
################### Volcano ###################
###############################################

ylx <- max(c(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5,]$Soil_Pase_umol_g_hr_04),na.rm=TRUE)*1.06
plot(V_jTRT[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="DOVI"]-0.2,
	V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="DOVI",]$Soil_Pase_umol_g_hr_04,
	xlim=c(-0.5,3.5),ylim=c(0,ylx),
	pch=1,xaxt="n",col=nonfixcol,cex=1.5,cex.axis=0.8)
abline(h=0,col="gray")
points(V_jTRT[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="ACKO"],
	V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="ACKO",]$Soil_Pase_umol_g_hr_04,
	col=rhizcol,pch=2,cex=1.5)
points(V_jTRT[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04>=0.5 & V$Species=="MOFA",]$Soil_Pase_umol_g_hr_04,
	col=actcol,pch=0,cex=1.5)

DOVI_04_LN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
DOVI_04_LN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(-0.2,DOVI_04_LN_u,pch=16,col=nonfixcol,cex=2)
segments(-0.2,DOVI_04_LN_u - DOVI_04_LN_se,-0.2,DOVI_04_LN_u + DOVI_04_LN_se,col=nonfixcol)
DOVI_04_MN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
DOVI_04_MN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(0.8,DOVI_04_MN_u,pch=16,col=nonfixcol,cex=2)
segments(0.8,DOVI_04_MN_u - DOVI_04_MN_se,0.8,DOVI_04_MN_u + DOVI_04_MN_se,col=nonfixcol)
DOVI_04_HN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
DOVI_04_HN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(1.8,DOVI_04_HN_u,pch=16,col=nonfixcol,cex=2)
segments(1.8,DOVI_04_HN_u - DOVI_04_HN_se,1.8,DOVI_04_HN_u + DOVI_04_HN_se,col=nonfixcol)
DOVI_04_PHN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
DOVI_04_PHN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(2.8,DOVI_04_PHN_u,pch=16,col=nonfixcol,cex=2)
segments(2.8,DOVI_04_PHN_u - DOVI_04_PHN_se,2.8,DOVI_04_PHN_u + DOVI_04_PHN_se,col=nonfixcol)

ACKO_04_LN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ACKO_04_LN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0,ACKO_04_LN_u,pch=17,col=rhizcol,cex=2)
segments(0,ACKO_04_LN_u - ACKO_04_LN_se,0,ACKO_04_LN_u + ACKO_04_LN_se,col=rhizcol)
ACKO_04_MN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ACKO_04_MN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1,ACKO_04_MN_u,pch=17,col=rhizcol,cex=2)
segments(1,ACKO_04_MN_u - ACKO_04_MN_se,1,ACKO_04_MN_u + ACKO_04_MN_se,col=rhizcol)
ACKO_04_HN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ACKO_04_HN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2,ACKO_04_HN_u,pch=17,col=rhizcol,cex=2)
segments(2,ACKO_04_HN_u - ACKO_04_HN_se,2,ACKO_04_HN_u + ACKO_04_HN_se,col=rhizcol)
ACKO_04_PHN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
ACKO_04_PHN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3,ACKO_04_PHN_u,pch=17,col=rhizcol,cex=2)
segments(3,ACKO_04_PHN_u - ACKO_04_PHN_se,3,ACKO_04_PHN_u + ACKO_04_PHN_se,col=rhizcol)

MOFA_04_LN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
MOFA_04_LN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04))
points(0.2,MOFA_04_LN_u,pch=15,col=actcol,cex=2)
segments(0.2,MOFA_04_LN_u - MOFA_04_LN_se,0.2,MOFA_04_LN_u + MOFA_04_LN_se,col=actcol)
MOFA_04_MN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
MOFA_04_MN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04))
points(1.2,MOFA_04_MN_u,pch=15,col=actcol,cex=2)
segments(1.2,MOFA_04_MN_u - MOFA_04_MN_se,1.2,MOFA_04_MN_u + MOFA_04_MN_se,col=actcol)
MOFA_04_HN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
MOFA_04_HN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04))
points(2.2,MOFA_04_HN_u,pch=15,col=actcol,cex=2)
segments(2.2,MOFA_04_HN_u - MOFA_04_HN_se,2.2,MOFA_04_HN_u + MOFA_04_HN_se,col=actcol)
MOFA_04_PHN_u <- mean(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)
MOFA_04_PHN_se <- sd(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)/
  sqrt(length(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04))
points(3.2,MOFA_04_PHN_u,pch=15,col=actcol,cex=2)
segments(3.2,MOFA_04_PHN_u - MOFA_04_PHN_se,3.2,MOFA_04_PHN_u + MOFA_04_PHN_se,col=actcol)

mtext("e",at=-0.35,padj=1.5,cex=0.8)
mtext("Volcano 2018",at=1.5,cex=0.8)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5,cex.axis=0.8)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5,cex.axis=0.8)

summary(HI_V_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
	data=V[V$Use_growth_04==1,]))
anova(HI_V_FX_Pase_04_lm)

# a a a a a a a a
#text(c(-0.2,0.2,0.8,1.2,1.8,2.2,2.8,3.2),80,"A",cex=1.25)

mtext("Treatment", side=1, outer=T, at=0.52, padj=2, cex=0.8)
mtext(expression(Phosphatase ~ activity ~ "(umol" ~ g ~ ""^"-1" ~ dry ~ soil ~ hr ~ ""^"-1" ~ ")"),
	side=2, outer=T, at=0.5, padj=0, cex=0.8)

dev.off()

#############################################################
#############################################################
#############################################################