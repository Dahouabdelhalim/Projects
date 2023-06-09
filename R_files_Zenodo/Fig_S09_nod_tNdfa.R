# NSF FX data plotting

# Figure S9
# Nodule biomass against total Ndfa

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS09.pdf",
	width=3.3,height=3.3)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(1,1),omi=c(.4,.4,.1,.1))
par(mai=c(.2,.2,0,0))

source("NY_SNF_analysis_tNdfa_forpaper_base.R")
N <- dat
xl <- c(min(c(0,N[N$Use_growth_08==1 & !is.na(N$Nfix_g_N_yr_u_08),]$Nodule_drymass_mg_08),na.rm=TRUE),
	max(N[N$Use_growth_08==1 & !is.na(N$Nfix_g_N_yr_u_08),]$Nodule_drymass_mg_08,na.rm=TRUE))
yl <- c(min(c(0,N[N$Use_growth_08==1 & !is.na(N$Nfix_g_N_yr_u_08),]$Nfix_g_N_yr_u_08)),
	max(N[N$Use_growth_08==1 & !is.na(N$Nfix_g_N_yr_u_08),]$Nfix_g_N_yr_u_08))

# Analysis section: 

summary(lm_tNdfa_nod <- lm(Nfix_g_N_yr_u_08 ~ Nodule_drymass_mg_08,
	data=N[N$Use_growth_08==1,]))
tab <- summary(lm_tNdfa_nod)[[4]]

# Plotting

# 2019
plot(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Nodule_drymass_mg_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Nfix_g_N_yr_u_08,
	xlim=xl,ylim=yl,#xaxt="n",yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Nodule_drymass_mg_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Nfix_g_N_yr_u_08,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Nodule_drymass_mg_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Nfix_g_N_yr_u_08,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Nodule_drymass_mg_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Nfix_g_N_yr_u_08,
	pch=5,col="red",cex=1.5)

points(xl,tab[1,1] + tab[2,1]*xl,
	type="l",col="red",lty=1)
text(125,yl[2]*0.9,c("ns"),col="red",cex=1.3)

mtext(expression(Nodule ~ biomass ~ (mg)), side=1, outer=T, at=0.52, padj=1.2)
mtext(expression(Nitrogen ~ fixation ~ (g ~ N ~ tree^-1 ~ y^-1)),side=2, outer=T, at=0.5, padj=-0.8)

dev.off()

#############################################################
#############################################################
#############################################################