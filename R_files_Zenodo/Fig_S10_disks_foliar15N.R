# NSF FX data plotting

# Figure S10
# Cellulose disks vs. foliar atom % 15N

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS10.pdf",
	width=7.05,height=2)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(1,4),omi=c(.4,.2,.2,.01),mai=c(.1,.3,0,0))

N <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]
O <- read.csv("OR_FX_Size_FoliarCNIsotope_Data.csv")[1:64,1:101]
W <- read.csv("HI_W_FX_Size_FoliarCNIsotope_Data.csv")[1:108,1:111]
V <- read.csv("HI_V_FX_Size_FoliarCNIsotope_Data.csv")[1:96,1:111]

##############################################
################## New York ##################
##############################################

xl <- c(0,1.1*max(N$disk_d15N_07,na.rm=TRUE))
yl <- c(0,1.1*max(N[N$Species=="BENI" & !is.na(N$disk_d15N_07),]$foliar_d15N_07,na.rm=TRUE))
plot(N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1,]$disk_d15N_07,
	N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1,]$foliar_d15N_07,
	xlim=xl,ylim=yl,
	pch=1,col="blue",cex=1.5)
points(N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07<1,]$disk_d15N_07,
	N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07<1,]$foliar_d15N_07,
	pch=19,col="blue",cex=0.5)
abline(0,1,lty=3)
mtext("a",at=xl[2]*0.05,padj=1.5)
mtext("New York",at=xl[2]*0.5)

summary(N_lm <- lm(foliar_d15N_07 ~ disk_d15N_07,data=N[N$Use_growth_07==1,]))
abline(N_lm,lty=2)
print(N_r2 <- round(summary(N_lm)$adj.r.squared,2))
text(xl[2]*0.5,yl[2]*0.1,expression(paste("R",""^2,"< 0.01")),pos=4)

##############################################
################### Oregon ###################
##############################################

xl <- c(0,1.1*max(O$disk_d15N_04,na.rm=TRUE))
yl <- c(0,1.1*max(O[O$Species=="PSME",]$foliar_d15N_04,na.rm=TRUE))
plot(O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1,]$disk_d15N_04,
	O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1,]$foliar_d15N_04,
	xlim=xl,ylim=yl,
	pch=1,col="blue",cex=1.5)
abline(0,1,lty=3)
mtext("b",at=xl[2]*0.05,padj=1.5)
mtext("Oregon",at=xl[2]*0.5)

summary(O_lm <- lm(foliar_d15N_04 ~ disk_d15N_04,data=O[O$Use_growth_04==1,]))
abline(O_lm)
print(O_r2 <- round(summary(O_lm)$adj.r.squared,2))
text(xl[2]*0.5,yl[2]*0.1,expression(paste("R",""^2,"= 0.69")),pos=4)

###############################################
################### Waiakea ###################
###############################################

xl <- c(0,1.1*max(W$disk_d15N_04,na.rm=TRUE))
yl <- c(0,1.1*max(W[W$Species=="PSCA",]$foliar_d15N_04,na.rm=TRUE))
plot(W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5,]$disk_d15N_04,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5,]$foliar_d15N_04,
	xlim=xl,ylim=yl,
	pch=1,col="blue",cex=1.5)
points(W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04==0,]$disk_d15N_04,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04==0,]$foliar_d15N_04,
	pch=19,col="blue",cex=0.5)
abline(0,1,lty=3)
mtext("c",at=xl[2]*0.05,padj=1.5)
mtext("Waiakea",at=xl[2]*0.5)

summary(W_lm <- lm(foliar_d15N_04 ~ disk_d15N_04,data=W[W$Use_growth_04==1,]))
abline(W_lm,lty=2)
print(W_r2 <- round(summary(W_lm)$adj.r.squared,2))
text(xl[2]*0.5,yl[2]*0.15,expression(paste("R",""^2,"= 0.07")),pos=4)

summary(W_lm_nooutlier <- lm(foliar_d15N_04 ~ disk_d15N_04,
	data=W[W$Use_growth_04==1 & W$disk_d15N_04 < 8000,]))
abline(W_lm_nooutlier,lty=1)
print(W_r2_nooutlier <- round(summary(W_lm_nooutlier)$adj.r.squared,2))
text(xl[2]*0.5,yl[2]*0.3,expression(paste("R",""^2,"= 0.12")),pos=4)

###############################################
################### Volcano ###################
###############################################

xl <- c(0,1.1*max(V$disk_d15N_04,na.rm=TRUE))
yl <- c(0,1.1*max(V[V$Species=="DOVI",]$foliar_d15N_04,na.rm=TRUE))
plot(V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1,]$disk_d15N_04,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1,]$foliar_d15N_04,
	xlim=xl,ylim=yl,
	pch=1,col="blue",cex=1.5)
points(V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04<1,]$disk_d15N_04,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04<1,]$foliar_d15N_04,
	pch=19,col="blue",cex=0.5)
abline(0,1,lty=3)
mtext("d",at=xl[2]*0.05,padj=1.5)
mtext("Volcano",at=xl[2]*0.5)

summary(V_lm <- lm(foliar_d15N_04 ~ disk_d15N_04,data=V[V$Use_growth_04==1,]))
abline(V_lm)
print(V_r2 <- round(summary(V_lm)$adj.r.squared,2))
text(xl[2]*0.5,yl[2]*0.9,expression(paste("R",""^2,"= 0.44")),pos=4)

mtext(expression(paste("Cellulose disk ",delta,""^15*N)), side=1, outer=T, at=0.52, padj=1.5)
mtext(expression(paste("Foliar ",delta,""^15*N)),side=2, outer=T, at=0.5, padj=0)


dev.off()

#############################################################
#############################################################
#############################################################