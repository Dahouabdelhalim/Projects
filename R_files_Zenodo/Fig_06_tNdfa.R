# NSF FX data plotting

# Figure 6
# N fixation per biomass (average across years)

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="Fig06.pdf",
	width=7.05,height=6)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,4),omi=c(.4,.5,.2,.2))
par(mai=c(0,0,0,0))

##############################################
################## New York ##################
##############################################

source("NY_SNF_analysis_tNdfa_forpaper_base.R")
N <- dat
N_TRT <- rep(-1,nrow(N))
N_TRT[N$Treatment=="LN"] <- 0
N_TRT[N$Treatment=="MN"] <- 1
N_TRT[N$Treatment=="HN"] <- 2
N_TRT[N$Treatment=="PHN"] <- 3
N_jTRT <- jitter(N_TRT)

# 2017
plot(N_jTRT[!is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & dat$Use_growth_04==1],
	dat[!is.na(dat$Nfix_g_N_kg_biomass_yr_u_04) & dat$Use_growth_04==1,]$Nfix_g_N_kg_biomass_yr_u_04,
	xaxt="n",xlim=c(-0.5,3.5),ylim=c(0,16.4),
	pch=2,col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("a",at=-0.45,padj=1.5)
mtext("Year 3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2018

plot(N_jTRT[!is.na(N$Nfix_g_N_kg_biomass_yr_u_07) & N$Use_growth_07==1],
	N[!is.na(N$Nfix_g_N_kg_biomass_yr_u_07) & N$Use_growth_07==1,]$Nfix_g_N_kg_biomass_yr_u_07,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),ylim=c(0,16.4),
	pch=2,col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("b",at=-0.45,padj=1.5)
mtext("Year 4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019

plot(N_jTRT[!is.na(N$Nfix_g_N_kg_biomass_yr_u_08) & N$Use_growth_08==1],
	N[!is.na(N$Nfix_g_N_kg_biomass_yr_u_08) & N$Use_growth_08==1,]$Nfix_g_N_kg_biomass_yr_u_08,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),ylim=c(0,16.4),
	pch=2,col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("c",at=-0.45,padj=1.5)
mtext("Year 5",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# All years
plot(0:3,c(tNdfa_kg_y_base_LN[1],tNdfa_kg_y_base_MN[1],
	tNdfa_kg_y_base_HN[1],tNdfa_kg_y_base_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(0,16.4),
	pch=17,xaxt="n",yaxt="n",col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(0,tNdfa_kg_y_base_LN[2],0,tNdfa_kg_y_base_LN[3],col="red")
segments(1,tNdfa_kg_y_base_MN[2],1,tNdfa_kg_y_base_MN[3],col="red")
segments(2,tNdfa_kg_y_base_HN[2],2,tNdfa_kg_y_base_HN[3],col="red")
segments(3,tNdfa_kg_y_base_PHN[2],3,tNdfa_kg_y_base_PHN[3],col="red")
text(0:3,rep(10,4),c("a","b","b","b"),col="red",cex=1.3)
mtext("d",at=-0.45,padj=1.5)
mtext("All years",at=1.5)
mtext("New York",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")


##############################################
################### Oregon ###################
##############################################

source("OR_SNF_analysis_tNdfa_forpaper_base.R")
O <- dat
O_TRT <- rep(-1,nrow(O))
O_TRT[O$Treatment=="LN"] <- 0
O_TRT[O$Treatment=="MN"] <- 1
O_TRT[O$Treatment=="HN"] <- 2
O_TRT[O$Treatment=="PHN"] <- 3
O_jTRT <- jitter(O_TRT)

# 2018 
plot(O_jTRT[!is.na(O$Nfix_g_N_kg_biomass_yr_u_04) & O$Use_growth_04==1],
	O[!is.na(O$Nfix_g_N_kg_biomass_yr_u_04) & O$Use_growth_04==1,]$Nfix_g_N_kg_biomass_yr_u_04,
	xaxt="n",xlim=c(-0.5,3.5),
	ylim=c(0,5),
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("e",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(O_jTRT[!is.na(O$Nfix_g_N_kg_biomass_yr_u_05) & O$Use_growth_05==1],
	O[!is.na(O$Nfix_g_N_kg_biomass_yr_u_05) & O$Use_growth_05==1,]$Nfix_g_N_kg_biomass_yr_u_05,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),
	ylim=c(0,5),
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("f",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2020
plot(O_jTRT[!is.na(O$Nfix_g_N_kg_biomass_yr_u_06) & O$Use_growth_06==1],
	O[!is.na(O$Nfix_g_N_kg_biomass_yr_u_06) & O$Use_growth_06==1,]$Nfix_g_N_kg_biomass_yr_u_06,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),
	ylim=c(0,5),
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
mtext("g",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# All years
plot(0:3,c(tNdfa_kg_y_base_LN[1],tNdfa_kg_y_base_MN[1],
	tNdfa_kg_y_base_HN[1],tNdfa_kg_y_base_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(0,5),
	pch=15,xaxt="n",yaxt="n",col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(0,tNdfa_kg_y_base_LN[2],0,tNdfa_kg_y_base_LN[3],col="orange")
segments(1,tNdfa_kg_y_base_MN[2],1,tNdfa_kg_y_base_MN[3],col="orange")
segments(2,tNdfa_kg_y_base_HN[2],2,tNdfa_kg_y_base_HN[3],col="orange")
segments(3,tNdfa_kg_y_base_PHN[2],3,tNdfa_kg_y_base_PHN[3],col="orange")
text(0:3,rep(4,4),c("a","b","b","b"),col="orange",cex=1.3)
mtext("h",at=-0.45,padj=1.5)
mtext("Oregon",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Waiakea ###################
###############################################

source("HI_W_SNF_analysis_tNdfa_forpaper_base.R")
W <- dat
W_TRT <- rep(-1,nrow(W))
W_TRT[W$Treatment=="LN"] <- 0
W_TRT[W$Treatment=="MN"] <- 1
W_TRT[W$Treatment=="HN"] <- 2
W_TRT[W$Treatment=="PHN"] <- 3
W_jTRT <- jitter(W_TRT,amount=0.1)

# 2018
plot(W_jTRT[!is.na(W$Nfix_g_N_kg_biomass_yr_u_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$Nfix_g_N_kg_biomass_yr_u_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE",]$Nfix_g_N_kg_biomass_yr_u_04,
	xaxt="n",xlim=c(-0.5,3.5),
	ylim=c(-0.5,14.3),
	pch=2,col="red",cex=1.5)
points(W_jTRT[!is.na(W$Nfix_g_N_kg_biomass_yr_u_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$Nfix_g_N_kg_biomass_yr_u_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$Nfix_g_N_kg_biomass_yr_u_04,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("i",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(W_jTRT[!is.na(W$Nfix_g_N_kg_biomass_yr_u_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$Nfix_g_N_kg_biomass_yr_u_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE",]$Nfix_g_N_kg_biomass_yr_u_06,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),
	ylim=c(-0.5,14.3),
	pch=2,col="red",cex=1.5)
points(W_jTRT[!is.na(W$Nfix_g_N_kg_biomass_yr_u_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$Nfix_g_N_kg_biomass_yr_u_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ",]$Nfix_g_N_kg_biomass_yr_u_06,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
mtext("j",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()

points(.3,.75,pch=2,col="red",cex=1,xpd=TRUE)
text(.3,0.75,"Robinia",pos=4,col="red",font=3,xpd=TRUE,cex=1.3)
text(.3,0.75,"a-d",pos=2,font=1,xpd=TRUE,cex=1.3)

points(.3,.5,pch=0,col="orange",cex=1,xpd=TRUE)
text(.3,0.5,"Alnus",pos=4,col="orange",font=3,xpd=TRUE,cex=1.3)
text(.3,0.5,"e-h",pos=2,font=1,xpd=TRUE,cex=1.3)

points(c(.3,.3),c(.25,.15),pch=c(2,0),col=c("red","orange"),cex=1,xpd=TRUE)
text(.3,0.25,"Gliricidia",pos=4,col="red",font=3,xpd=TRUE,cex=1.3)
text(.3,0.15,"Casuarina",pos=4,col="orange",font=3,xpd=TRUE,cex=1.3)
text(.3,0.2,"i-k",pos=2,font=1,xpd=TRUE,cex=1.3)

# All years
plot(0:3-0.15,c(tNdfa_kg_y_base_GLSE_LN[1],tNdfa_kg_y_base_GLSE_MN[1],
	tNdfa_kg_y_base_GLSE_HN[1],tNdfa_kg_y_base_GLSE_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(-0.5,14.3),xaxt="n",yaxt="n",
	pch=17,col="red",cex=1.5)
points(0:3+0.15,c(tNdfa_kg_y_base_CAEQ_LN[1],tNdfa_kg_y_base_CAEQ_MN[1],
	tNdfa_kg_y_base_CAEQ_HN[1],tNdfa_kg_y_base_CAEQ_PHN[1]),
	pch=15,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(0-0.15,tNdfa_kg_y_base_GLSE_LN[2],0-0.15,tNdfa_kg_y_base_GLSE_LN[3],col="red")
segments(1-0.15,tNdfa_kg_y_base_GLSE_MN[2],1-0.15,tNdfa_kg_y_base_GLSE_MN[3],col="red")
segments(2-0.15,tNdfa_kg_y_base_GLSE_HN[2],2-0.15,tNdfa_kg_y_base_GLSE_HN[3],col="red")
segments(3-0.15,tNdfa_kg_y_base_GLSE_PHN[2],3-0.15,tNdfa_kg_y_base_GLSE_PHN[3],col="red")
segments(0+0.15,tNdfa_kg_y_base_CAEQ_LN[2],0+0.15,tNdfa_kg_y_base_CAEQ_LN[3],col="orange")
segments(1+0.15,tNdfa_kg_y_base_CAEQ_MN[2],1+0.15,tNdfa_kg_y_base_CAEQ_MN[3],col="orange")
segments(2+0.15,tNdfa_kg_y_base_CAEQ_HN[2],2+0.15,tNdfa_kg_y_base_CAEQ_HN[3],col="orange")
segments(3+0.15,tNdfa_kg_y_base_CAEQ_PHN[2],3+0.15,tNdfa_kg_y_base_CAEQ_PHN[3],col="orange")
text((0:3)-0.2,rep(8,4),c("a","a","a","a"),col="red",cex=1.3)
text((0:3)+0.2,rep(3,4),c("a","a","a","a"),col="orange",cex=1.3)
mtext("k",at=-0.45,padj=1.5)
mtext("Waiakea",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Volcano ###################
###############################################

source("HI_V_SNF_analysis_tNdfa_forpaper_base.R")
V <- dat
V_TRT <- rep(-1,nrow(V))
V_TRT[V$Treatment=="LN"] <- 0
V_TRT[V$Treatment=="MN"] <- 1
V_TRT[V$Treatment=="HN"] <- 2
V_TRT[V$Treatment=="PHN"] <- 3
V_jTRT <- jitter(V_TRT,amount=0.1)

# 2018
plot(V_jTRT[!is.na(V$Nfix_g_N_kg_biomass_yr_u_04) & V$Use_growth_04==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$Nfix_g_N_kg_biomass_yr_u_04) & V$Use_growth_04==1 & V$Species=="ACKO",]$Nfix_g_N_kg_biomass_yr_u_04,
	xaxt="n",xlim=c(-0.5,3.5),
	ylim=c(-1.5,4.7),
	pch=2,col="red",cex=1.5)
points(V_jTRT[!is.na(V$Nfix_g_N_kg_biomass_yr_u_04) & V$Use_growth_04==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$Nfix_g_N_kg_biomass_yr_u_04) & V$Use_growth_04==1 & V$Species=="MOFA",]$Nfix_g_N_kg_biomass_yr_u_04,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
mtext("l",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(V_jTRT[!is.na(V$Nfix_g_N_kg_biomass_yr_u_06) & V$Use_growth_06==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$Nfix_g_N_kg_biomass_yr_u_06) & V$Use_growth_06==1 & V$Species=="ACKO",]$Nfix_g_N_kg_biomass_yr_u_06,
	xaxt="n",yaxt="n",xlim=c(-0.5,3.5),
	ylim=c(-1.5,4.7),
	pch=2,col="red",cex=1.5)
points(V_jTRT[!is.na(V$Nfix_g_N_kg_biomass_yr_u_06) & V$Use_growth_06==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$Nfix_g_N_kg_biomass_yr_u_06) & V$Use_growth_06==1 & V$Species=="MOFA",]$Nfix_g_N_kg_biomass_yr_u_06,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
mtext("m",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()
points(c(.3,.3),c(.55,.45),pch=c(2,0),col=c("red","orange"),cex=1,xpd=TRUE)
text(.3,0.55,"Acacia",pos=4,col="red",font=3,xpd=TRUE,cex=1.3)
text(.3,0.45,"Morella",pos=4,col="orange",font=3,xpd=TRUE,cex=1.3)
text(.3,0.5,"l-n",pos=2,font=1,xpd=TRUE,cex=1.3)

# All years
plot(0:3-0.15,c(tNdfa_kg_y_base_ACKO_LN[1],tNdfa_kg_y_base_ACKO_MN[1],
	tNdfa_kg_y_base_ACKO_HN[1],tNdfa_kg_y_base_ACKO_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(-1.5,4.7),xaxt="n",yaxt="n",
	pch=17,col="red",cex=1.5)
points(0:3+0.15,c(tNdfa_kg_y_base_MOFA_LN[1],tNdfa_kg_y_base_MOFA_MN[1],
	tNdfa_kg_y_base_MOFA_HN[1],tNdfa_kg_y_base_MOFA_PHN[1]),
	pch=15,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
segments(0-0.15,tNdfa_kg_y_base_ACKO_LN[2],0-0.15,tNdfa_kg_y_base_ACKO_LN[3],col="red")
segments(1-0.15,tNdfa_kg_y_base_ACKO_MN[2],1-0.15,tNdfa_kg_y_base_ACKO_MN[3],col="red")
segments(2-0.15,tNdfa_kg_y_base_ACKO_HN[2],2-0.15,tNdfa_kg_y_base_ACKO_HN[3],col="red")
segments(3-0.15,tNdfa_kg_y_base_ACKO_PHN[2],3-0.15,tNdfa_kg_y_base_ACKO_PHN[3],col="red")
segments(0+0.15,tNdfa_kg_y_base_MOFA_LN[2],0+0.15,tNdfa_kg_y_base_MOFA_LN[3],col="orange")
segments(1+0.15,tNdfa_kg_y_base_MOFA_MN[2],1+0.15,tNdfa_kg_y_base_MOFA_MN[3],col="orange")
segments(2+0.15,tNdfa_kg_y_base_MOFA_HN[2],2+0.15,tNdfa_kg_y_base_MOFA_HN[3],col="orange")
segments(3+0.15,tNdfa_kg_y_base_MOFA_PHN[2],3+0.15,tNdfa_kg_y_base_MOFA_PHN[3],col="orange")
text((0:3)-0.2,rep(3.7,4),c("a","a","a","a"),col="red",cex=1.3)
text((0:3)+0.2,rep(2.5,4),c("a","a","a","a"),col="orange",cex=1.3)
mtext("n",at=-0.45,padj=1.5)
mtext("Volcano",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")


mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Nitrogen ~ fixation ~ (g ~ N ~ kg ~ biomass^-1 ~ yr^-1)),side=2, outer=T, at=0.5, padj=-1.4)


dev.off()

#############################################################
#############################################################
#############################################################