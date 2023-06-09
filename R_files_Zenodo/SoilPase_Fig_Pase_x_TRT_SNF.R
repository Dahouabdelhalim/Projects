# NSF FX data plotting

# Figure 4
# Soil phosphatase by treatment and species and symbiotic N fixation

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

library(colorBlindness)

#########################################################################
############################# Set up figure #############################
#########################################################################

#pdf(file="Pase_Fig04.pdf",
#	width=6,height=4)
setEPS()
postscript(file="Pase_Fig04.eps",
           width=6,height=4)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(2,3),omi=c(.5,.25,.2,.2),mai=c(.2,.3,0.2,0))

N <- read.csv("NY_FX_Soil_Pase.csv")
N_SNF <- read.csv("NY_FX_SNF.csv")
N <- cbind(N,N_SNF[,c(8,11,14:51)])

O <- read.csv("OR_FX_Soil_Pase.csv")
O_SNF <- read.csv("OR_FX_SNF.csv")
O <- cbind(O,O_SNF[,c(8,11,14:44)])

W <- read.csv("HI_W_FX_Soil_Pase.csv")
W_SNF <- read.csv("HI_W_FX_SNF.csv")
W <- cbind(W,W_SNF[,c(8,11:29)])

V <- read.csv("HI_V_FX_Soil_Pase.csv")
V_SNF <- read.csv("HI_V_FX_SNF.csv")
V <- cbind(V,V_SNF[,c(8,11:29)])

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

xln <- min(c(0,N[N$Use_growth_04==1,]$Nfix_g_N_yr_u_04),na.rm=TRUE)
xlx <- max(N[N$Use_growth_04==1,]$Nfix_g_N_yr_u_04,na.rm=TRUE)
ylx <- max(N[!is.na(N$Soil_Pase_umol_g_hr_04) & N$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)*1.02

plot(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	xlim=c(xln,xlx),ylim=c(0,ylx),#xaxt="n",#yaxt="n",
	pch=2,col=rhizcol,cex=1.5,cex.axis=0.8)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=rhizcol,cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=rhizcol,cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=rhizcol,cex=1.5)
abline(v=0,col="gray")
abline(h=0,col="gray")
points(rep(0,nrow(N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="LN",])),
	N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="MN",])),
	N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="HN",])),
	N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="PHN",])),
	N[N$Use_growth_04==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=nonfixcol,cex=1.5)
mtext("a",at=xlx/10,padj=1.5,cex=0.8)
mtext("New York 2017",at=mean(c(xln,xlx)),cex=0.8)

summary(NY_FX_Pase_Nfix_g_N_yr_u_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ 
	Nfix_g_N_yr_u_04 * Treatment,data=N[N$Use_growth_04==1 & N$Species=="ROPS",]))
anova(NY_FX_Pase_Nfix_g_N_yr_u_04_lm)
# Nothing significant

xln <- min(c(0,N[N$Use_growth_07==1,]$Nfix_g_N_yr_u_07),na.rm=TRUE)
xlx <- max(N[N$Use_growth_07==1,]$Nfix_g_N_yr_u_07,na.rm=TRUE)
ylx <- max(N[!is.na(N$Soil_Pase_umol_g_hr_07) & N$Use_growth_07==1,]$Soil_Pase_umol_g_hr_07,na.rm=TRUE)*1.02

plot(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Nfix_g_N_yr_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,
	xlim=c(xln,xlx),ylim=c(0,ylx),#xaxt="n",yaxt="n",
	pch=2,col=rhizcol,cex=1.5,cex.axis=0.8)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Nfix_g_N_yr_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,
	pch=3,col=rhizcol,cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Nfix_g_N_yr_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,
	pch=4,col=rhizcol,cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Nfix_g_N_yr_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,
	pch=5,col=rhizcol,cex=1.5)
abline(v=0,col="gray")
abline(h=0,col="gray")
points(rep(0,nrow(N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="LN",])),
	N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="LN",]$Soil_Pase_umol_g_hr_07,
	pch=2,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="MN",])),
	N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="MN",]$Soil_Pase_umol_g_hr_07,
	pch=3,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="HN",])),
	N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="HN",]$Soil_Pase_umol_g_hr_07,
	pch=4,col=nonfixcol,cex=1.5)
points(rep(0,nrow(N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="PHN",])),
	N[N$Use_growth_07==1 & N$Species=="BENI" & N$Treatment=="PHN",]$Soil_Pase_umol_g_hr_07,
	pch=5,col=nonfixcol,cex=1.5)
mtext("b",at=xlx/10,padj=1.5,cex=0.8)
mtext("New York 2018",at=mean(c(xln,xlx)),cex=0.8)

summary(NY_FX_Pase_Nfix_g_N_yr_u_07_lm <- lm(Soil_Pase_umol_g_hr_07 ~ 
	Nfix_g_N_yr_u_07 * Treatment,data=N[N$Use_growth_07==1 & N$Species=="ROPS",]))
anova(NY_FX_Pase_Nfix_g_N_yr_u_07_lm)
# Nothing significant

##############################################
################### Oregon ###################
##############################################

xln <- min(c(0,O[O$Use_growth_04==1,]$Nfix_g_N_yr_u_04),na.rm=TRUE)
xlx <- max(O[O$Use_growth_04==1,]$Nfix_g_N_yr_u_04,na.rm=TRUE)
ylx <- max(O[!is.na(O$Soil_Pase_umol_g_hr_04) & O$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)*1.02

plot(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	xlim=c(xln,xlx),ylim=c(0,ylx),#xaxt="n",yaxt="n",
	pch=2,col=actcol,cex=1.5,cex.axis=0.8)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=actcol,cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=actcol,cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=actcol,cex=1.5)
abline(v=0,col="gray")
abline(h=0,col="gray")
points(rep(0,nrow(O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="LN",])),
	O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=nonfixcol,cex=1.5)
points(rep(0,nrow(O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="MN",])),
	O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=nonfixcol,cex=1.5)
points(rep(0,nrow(O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="HN",])),
	O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=nonfixcol,cex=1.5)
points(rep(0,nrow(O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="PHN",])),
	O[O$Use_growth_04==1 & O$Species=="PSME" & O$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=nonfixcol,cex=1.5)

mtext("c",at=xlx/10,padj=1.5,cex=0.8)
mtext("Oregon 2018",at=mean(c(xln,xlx)),cex=0.8)

summary(OR_FX_Pase_Nfix_g_N_yr_u_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ 
	Nfix_g_N_yr_u_04 * Treatment,data=O[O$Use_growth_04==1 & O$Species=="ALRU",]))
anova(OR_FX_Pase_Nfix_g_N_yr_u_04_lm)
# Treatment significant, but lm not significant

###############################################
################### Waiakea ###################
###############################################

xln <- min(c(0,W[W$Use_growth_04==1,]$Nfix_g_N_yr_u_04),na.rm=TRUE)
xlx <- max(W[W$Use_growth_04==1,]$Nfix_g_N_yr_u_04,na.rm=TRUE)
ylx <- max(W[!is.na(W$Soil_Pase_umol_g_hr_04) & W$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)*1.02

plot(W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	xlim=c(xln,xlx),ylim=c(0,ylx),#xaxt="n",yaxt="n",
	pch=2,col=rhizcol,cex=1.5,cex.axis=0.8)
points(W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=rhizcol,cex=1.5)
points(W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=rhizcol,cex=1.5)
points(W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="GLSE" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=rhizcol,cex=1.5)
abline(v=0,col="gray")
abline(h=0,col="gray")
points(W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=actcol,cex=1.5)
points(W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=actcol,cex=1.5)
points(W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=actcol,cex=1.5)
points(W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	W[W$Use_growth_04==1 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=actcol,cex=1.5)

points(rep(0,nrow(W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="LN",])),
	W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=nonfixcol,cex=1.5)
points(rep(0,nrow(W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="MN",])),
	W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=nonfixcol,cex=1.5)
points(rep(0,nrow(W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="HN",])),
	W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=nonfixcol,cex=1.5)
points(rep(0,nrow(W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="PHN",])),
	W[W$Use_growth_04==1 & W$Species=="PSCA" & W$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=nonfixcol,cex=1.5)

mtext("d",at=xlx/10,padj=1.5,cex=0.8)
mtext("Waiakea 2018",at=mean(c(xln,xlx)),cex=0.8)

summary(HI_W_FX_Pase_Nfix_g_N_yr_u_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ 
	Nfix_g_N_yr_u_04 * Species * Treatment,data=W[W$Use_growth_04==1 & W$Species!="PSCA",]))
anova(HI_W_FX_Pase_Nfix_g_N_yr_u_04_lm)
# Says that treatment is signficant but the lm itself isn't significant

###############################################
################### Legend ####################
###############################################

plot.new()

text(0.4,0.95,"Robinia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.4,0.85,"Betula",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.44,0.9,"a-b",pos=2,xpd=TRUE,cex=1)

text(0.4,0.7,"Alnus",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.4,0.6,"Pseudotsuga",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.44,0.65,"c",pos=2,xpd=TRUE,cex=1)

text(0.4,0.45,"Gliricidia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.4,0.35,"Casuarina",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.4,0.25,"Psidium",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.44,0.35,"d",pos=2,xpd=TRUE,cex=1)

text(0.4,0.1,"Acacia",pos=4,col=rhizcol,font=3,xpd=TRUE,cex=1)
text(0.4,0.0,"Morella",pos=4,col=actcol,font=3,xpd=TRUE,cex=1)
text(0.4,-0.1,"Dodonaea",pos=4,col=nonfixcol,font=3,xpd=TRUE,cex=1)
text(0.44,0.0,"e",pos=2,xpd=TRUE,cex=1)

points(c(0.0,0.0,0.0,0.0)-0.1,c(0.2,0.1,0.0,-0.1),pch=c(2,3,4,5),col="black",cex=1,xpd=TRUE)
text(c(0.025,0.025,0.025,0.025)-0.1,c(0.2,0.1,0.0,-0.1),c("Control","+10","+15","+15+P"),
	pos=4,col="black",xpd=TRUE,cex=1)

###############################################
################### Volcano ###################
###############################################

xln <- min(c(0,V[V$Use_growth_04==1,]$Nfix_g_N_yr_u_04),na.rm=TRUE)
xlx <- max(V[V$Use_growth_04==1,]$Nfix_g_N_yr_u_04,na.rm=TRUE)
ylx <- max(V[!is.na(V$Soil_Pase_umol_g_hr_04) & V$Use_growth_04==1,]$Soil_Pase_umol_g_hr_04,na.rm=TRUE)*1.02

plot(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	xlim=c(xln,xlx),ylim=c(0,ylx),#xaxt="n",yaxt="n",
	pch=2,col=rhizcol,cex=1.5,cex.axis=0.8)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=rhizcol,cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=rhizcol,cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=rhizcol,cex=1.5)
abline(v=0,col="gray")
abline(h=0,col="gray")
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=actcol,cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=actcol,cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=actcol,cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Nfix_g_N_yr_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=actcol,cex=1.5)

points(rep(0,nrow(V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="LN",])),
	V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="LN",]$Soil_Pase_umol_g_hr_04,
	pch=2,col=nonfixcol,cex=1.5)
points(rep(0,nrow(V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="MN",])),
	V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="MN",]$Soil_Pase_umol_g_hr_04,
	pch=3,col=nonfixcol,cex=1.5)
points(rep(0,nrow(V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="HN",])),
	V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="HN",]$Soil_Pase_umol_g_hr_04,
	pch=4,col=nonfixcol,cex=1.5)
points(rep(0,nrow(V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="PHN",])),
	V[V$Use_growth_04==1 & V$Species=="DOVI" & V$Treatment=="PHN",]$Soil_Pase_umol_g_hr_04,
	pch=5,col=nonfixcol,cex=1.5)

mtext("e",at=xlx/10,padj=1.5,cex=0.8)
mtext("Volcano 2018",at=mean(c(xln,xlx)),cex=0.8)

summary(HI_V_FX_Pase_Nfix_g_N_yr_u_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ 
	Nfix_g_N_yr_u_04 * Species * Treatment,data=V[V$Use_growth_04==1 & V$Species!="DOVI",]))
anova(HI_V_FX_Pase_Nfix_g_N_yr_u_04_lm)
# Nothing significant



mtext(expression(Nitrogen ~ fixation ~ (g ~ N ~ tree^-1 ~ y^-1)), side=1, outer=T, at=0.52, padj=1.625, cex=0.8)
mtext(expression(Phosphatase ~ activity ~ "(umol" ~ g ~ ""^"-1" ~ dry ~ soil ~ hr ~ ""^"-1" ~ ")"),
	side=2, outer=T, at=0.5, padj=0, cex=0.8)

dev.off()

#############################################################
#############################################################
#############################################################