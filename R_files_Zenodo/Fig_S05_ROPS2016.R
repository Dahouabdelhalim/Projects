# NSF FX data plotting

# Figure S5
# Foliar atom % 15N and %Ndfa for Robinia 2016

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS05.pdf",width=3.3,height=2)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(1,2),omi=c(.4,0,.2,0.02),mai=c(0,0.5,0,0))

purefix_15N_AP <- 0.3663

N <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]
N_TRT <- rep(-1,nrow(N))
N_TRT[N$Treatment=="LN"] <- 0
N_TRT[N$Treatment=="MN"] <- 1
N_TRT[N$Treatment=="HN"] <- 2
N_TRT[N$Treatment=="PHN"] <- 3
N_jTRT <- jitter(N_TRT)

# Set unhealthy BENIs foliar d15N to NA
N[N$Species=="BENI" & N$Use_growth_02==0,]$foliar_15N_AP_02 <- NA

# Set up columns for all years and set the reference to the paired tree
N$ref_15N_AP_u_02 <- NA
N[N$Species=="ROPS",]$ref_15N_AP_u_02 <- N[N$Species=="BENI",]$foliar_15N_AP_02

# Foliar 15N
ylx <- max(c(N[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1 & N$Species=="ROPS",]$ref_15N_AP_u_02,
	N[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1 & N$Species=="ROPS",]$foliar_15N_AP_02),na.rm=TRUE)*1.02
plot(N_jTRT[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1],
	N[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1,]$ref_15N_AP_u_02,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),ylab="",
	pch=1,xaxt="n",yaxt="n",col="blue",cex=1)
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-2,cex.axis=0.6)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-2,cex.axis=0.6)
axis(side=2, at=seq(0.4,1.6,0.2), labels=TRUE, padj=1,cex.axis=0.6)
abline(h=purefix_15N_AP,col="gray")
points(N_jTRT[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1 & N$Species=="ROPS"],
	N[!is.na(N$foliar_15N_AP_02) & N$Use_growth_02==1 & N$Species=="ROPS",]$foliar_15N_AP_02,
	col="red",pch=2,cex=1)
mtext("a",at=-0.35,padj=1.5)
mtext("Year 2",at=1.5,cex.lab=0.7)
mtext(expression(Atom ~ "%" ~ ""^"15"*N), side=2, outer=FALSE, padj=-1.75, cex.lab=0.7)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# %Ndfa

source("NY_SNF_analysis_pNdfa_forpaper_base.R")
N <- dat
N_TRT <- rep(-1,nrow(N))
N_TRT[N$Treatment=="LN"] <- 0
N_TRT[N$Treatment=="MN"] <- 1
N_TRT[N$Treatment=="HN"] <- 2
N_TRT[N$Treatment=="PHN"] <- 3
N_jTRT <- jitter(N_TRT)

plot(N_jTRT[!is.na(N$Ndfa_u_02) & N$Use_growth_02==1],
	N[!is.na(N$Ndfa_u_02) & N$Use_growth_02==1,]$Ndfa_u_02,
	xlim=c(-0.5,3.5),ylim=c(-53,103),ylab="",
	pch=2,xaxt="n",col="red",cex=1,yaxt="n")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-2,cex.axis=0.6)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-2,cex.axis=0.6)
axis(side=2, at=c(-50,0,50,100), labels=TRUE, padj=1,cex.axis=0.6)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("b",at=-0.45,padj=1.5)
mtext("Year 2",at=1.5,cex.lab=0.7)
mtext(expression("%" ~ N[dfa]), side=2, outer=FALSE, padj=-1.75, cex.lab=0.7)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")


mtext("Treatment", side=1, outer=T, at=0.6, padj=1.75,cex.lab=0.8)


dev.off()

#############################################################
#############################################################
#############################################################