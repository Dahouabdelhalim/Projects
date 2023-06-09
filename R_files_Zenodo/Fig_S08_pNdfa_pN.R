# NSF FX data plotting

# Figure S8
# % Ndfa against foliar %N

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS08.pdf",
	width=7.05,height=6)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,4),omi=c(.4,.5,.2,.2))
par(mai=c(0,0,0,0))

##############################################
################## New York ##################
##############################################

source("NY_SNF_analysis_pNdfa_forpaper_base.R")
N <- dat
yl <- c(min(c(N[N$Use_growth_04==1 & N$Species=="ROPS",]$foliar_N_mg_g_04,
	N[N$Use_growth_07==1 & N$Species=="ROPS",]$foliar_N_mg_g_07,
	N[N$Use_growth_08==1 & N$Species=="ROPS",]$foliar_N_mg_g_08),na.rm=TRUE),
	max(c(N[N$Use_growth_04==1 & N$Species=="ROPS",]$foliar_N_mg_g_04,
	N[N$Use_growth_07==1 & N$Species=="ROPS",]$foliar_N_mg_g_07,
	N[N$Use_growth_08==1 & N$Species=="ROPS",]$foliar_N_mg_g_08),na.rm=TRUE))

pN_ROPS_789 <- 
	c(
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_07) & 
		dat$Use_growth_07 == 1,]$foliar_N_mg_g_07,
	dat[dat$Species=="ROPS" & !is.na(dat$Ndfa_u_08) & 
		dat$Use_growth_08 == 1,]$foliar_N_mg_g_08
	)

# Analysis section: 
# We tried two different models. The main one is uncommented. The other 
# adds treatment as a fixed effect.

# None of them show significant effects of pNdfa on ln AGR of AGB for any
# species. The same models are coded for the other sites and species below, 
# though these notes aren't repeated.

#summary(lme_pNdfa_pN_TRT_base <- lme(pN_ROPS_789 ~ pNdfa_u_789 * Treatment_ROPS_789,
#	random=~1 | Tree_ROPS_789))
#tab <- summary(lme_pNdfa_pN_TRT_base)[[20]]
summary(lme_pNdfa_pN_base <- lme(pN_ROPS_789 ~ pNdfa_u_789,
	random=~1 | Tree_ROPS_789))
tab <- summary(lme_pNdfa_pN_base)[[20]]

# 2017
plot(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$foliar_N_mg_g_04,
	xaxt="n",xlim=c(-60,103),ylim=yl,#yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="red",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("a",at=-50,padj=1.5)
mtext("Year 3",at=mean(c(-60,103)))

# 2018
plot(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$foliar_N_mg_g_07,
	xaxt="n",xlim=c(-60,103),ylim=yl,yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$foliar_N_mg_g_07,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$foliar_N_mg_g_07,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$foliar_N_mg_g_07,
	pch=5,col="red",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("b",at=-50,padj=1.5)
mtext("Year 4",at=mean(c(-60,103)))

# 2019
plot(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$foliar_N_mg_g_08,
	xaxt="n",xlim=c(-60,103),ylim=yl,yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$foliar_N_mg_g_08,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$foliar_N_mg_g_08,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$foliar_N_mg_g_08,
	pch=5,col="red",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("c",at=-50,padj=1.5)
mtext("Year 5",at=mean(c(-60,103)))

# All years
plot(0:100,tab[1,1] + tab[2,1]*(0:100),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",
	type="l",col="red",lty=1)
abline(v=c(0,100),col="gray")
text(50,yl[2]*0.9,c("ns"),col="red",cex=1.3)
mtext("d",at=-50,padj=1.5)
mtext("All years",at=mean(c(-60,103)))
mtext("New York",side=4,padj=0.5)

##############################################
################### Oregon ###################
##############################################

source("OR_SNF_analysis_pNdfa_forpaper_base.R")
O <- dat
yl <- c(min(c(O[O$Use_growth_04==1 & O$Species=="ALRU",]$foliar_N_mg_g_04,
	O[O$Use_growth_05==1 & O$Species=="ALRU",]$foliar_N_mg_g_05,
	O[O$Use_growth_06==1 & O$Species=="ALRU",]$foliar_N_mg_g_06),na.rm=TRUE),
	max(c(O[O$Use_growth_04==1 & O$Species=="ALRU",]$foliar_N_mg_g_04,
	O[O$Use_growth_05==1 & O$Species=="ALRU",]$foliar_N_mg_g_05,
	O[O$Use_growth_06==1 & O$Species=="ALRU",]$foliar_N_mg_g_06),na.rm=TRUE))

pN_ALRU_890 <- 
	c(
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_05) & 
		dat$Use_growth_05 == 1,]$foliar_N_mg_g_05,
	dat[dat$Species=="ALRU" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
#summary(lme_pNdfa_pN_TRT_base <- lme(pN_ALRU_890 ~ pNdfa_u_890 * Treatment_ALRU_890,
#	random=~1 | Tree_ALRU_890))
#tab <- summary(lme_pNdfa_pN_TRT_base)[[20]]
summary(lme_pNdfa_pN_base <- lme(pN_ALRU_890 ~ pNdfa_u_890,
	random=~1 | Tree_ALRU_890))
tab <- summary(lme_pNdfa_pN_base)[[20]]

# 2018
plot(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$foliar_N_mg_g_04,
	xaxt="n",xlim=c(-60,103),ylim=yl,#yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("e",at=-50,padj=1.5)

# 2019
plot(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="LN",]$foliar_N_mg_g_05,
	xaxt="n",xlim=c(-60,103),ylim=yl,yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="MN",]$foliar_N_mg_g_05,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="HN",]$foliar_N_mg_g_05,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$foliar_N_mg_g_05,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("f",at=-50,padj=1.5)

# 2020
plot(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="LN",]$foliar_N_mg_g_06,
	xaxt="n",xlim=c(-60,103),ylim=yl,yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="MN",]$foliar_N_mg_g_06,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="HN",]$foliar_N_mg_g_06,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$foliar_N_mg_g_06,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
axis(side=1, c(0,50,100), labels=c("0","50","100  "), padj=-.5)
mtext("g",at=-50,padj=1.5)

# All years
plot(0:100,tab[1,1] + tab[2,1]*(0:100),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",
	type="l",col="orange",lty=1)
abline(v=c(0,100),col="gray")
text(50,yl[2]*0.9,c("ns"),col="orange",cex=1.3)
mtext("h",at=-50,padj=1.5)
mtext("Oregon",side=4,padj=0.5)

###############################################
################### Waiakea ###################
###############################################

source("HI_W_SNF_analysis_pNdfa_forpaper_base.R")
W <- dat
yl <- c(min(c(W[W$Use_growth_04>=0.5 & W$Species!="PSCA",]$foliar_N_mg_g_04,
	W[W$Use_growth_06>=0.5 & W$Species!="PSCA",]$foliar_N_mg_g_06),na.rm=TRUE),
	max(c(W[W$Use_growth_04>=0.5 & W$Species!="PSCA",]$foliar_N_mg_g_04,
	W[W$Use_growth_06>=0.5 & W$Species!="PSCA",]$foliar_N_mg_g_06),na.rm=TRUE))

pN_GLSE_89 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 >=0.5,]$foliar_N_mg_g_04,
	dat[dat$Species=="GLSE" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 >=0.5,]$foliar_N_mg_g_06
	)
#summary(lme_pNdfa_pN_TRT_GLSE_base <- lme(pN_GLSE_89 ~ pNdfa_u_GLSE_89 * Treatment_GLSE_89,
#	random=~1 | Tree_GLSE_89))
#tab_GLSE <- summary(lme_pNdfa_pN_TRT_GLSE_base)[[20]]
summary(lme_pNdfa_pN_GLSE_base <- lme(pN_GLSE_89 ~ pNdfa_u_GLSE_89,
	random=~1 | Tree_GLSE_89))
tab_GLSE <- summary(lme_pNdfa_pN_GLSE_base)[[20]]

pN_CAEQ_89 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 >=0.5,]$foliar_N_mg_g_04,
	dat[dat$Species=="CAEQ" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 >=0.5,]$foliar_N_mg_g_06
	)
#summary(lme_pNdfa_pN_TRT_CAEQ_base <- lme(pN_CAEQ_89 ~ pNdfa_u_CAEQ_89 * Treatment_CAEQ_89,
#	random=~1 | Tree_CAEQ_89))
#tab_CAEQ <- summary(lme_pNdfa_pN_TRT_CAEQ_base)[[20]]
summary(lme_pNdfa_pN_CAEQ_base <- lme(pN_CAEQ_89 ~ pNdfa_u_CAEQ_89,
	random=~1 | Tree_CAEQ_89))
tab_CAEQ <- summary(lme_pNdfa_pN_CAEQ_base)[[20]]

# 2018
plot(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$foliar_N_mg_g_04,
	xaxt="n",#yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$foliar_N_mg_g_04,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("i",at=-50,padj=1.5)

# 2019
plot(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$foliar_N_mg_g_06,
	xaxt="n",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$foliar_N_mg_g_06,
	pch=3,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$foliar_N_mg_g_06,
	pch=4,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$foliar_N_mg_g_06,
	pch=5,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$foliar_N_mg_g_06,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$foliar_N_mg_g_06,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$foliar_N_mg_g_06,
	pch=4,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$foliar_N_mg_g_06,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("j",at=-50,padj=1.5)

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
plot(0:100,tab_GLSE[1,1] + tab_GLSE[2,1]*(0:100),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",
	type="l",col="red",lty=1)
points(0:100,tab_CAEQ[1,1] + tab_CAEQ[2,1]*(0:100),
	type="l",col="orange",lty=1)
abline(v=c(0,100),col="gray")
text(40,yl[2]*0.9,c("ns"),col="red",cex=1.3)
text(60,yl[2]*0.9,c("ns"),col="orange",cex=1.3)
mtext("k",at=-50,padj=1.5)
mtext("Waiakea",side=4,padj=0.5)

###############################################
################### Volcano ###################
###############################################

source("HI_V_SNF_analysis_pNdfa_forpaper_base.R")
V <- dat
yl <- c(min(c(V[V$Use_growth_04>=0.5 & V$Species!="DOVI",]$foliar_N_mg_g_04,
	V[V$Use_growth_06>=0.5 & V$Species!="DOVI",]$foliar_N_mg_g_06),na.rm=TRUE),
	max(c(V[V$Use_growth_04>=0.5 & V$Species!="DOVI",]$foliar_N_mg_g_04,
	V[V$Use_growth_06>=0.5 & V$Species!="DOVI",]$foliar_N_mg_g_06),na.rm=TRUE))

pN_ACKO_89 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	dat[dat$Species=="ACKO" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
#summary(lme_pNdfa_pN_TRT_ACKO_base <- lme(pN_ACKO_89 ~ pNdfa_u_ACKO_89 * Treatment_ACKO_89,
#	random=~1 | Tree_ACKO_89))
#tab_ACKO <- summary(lme_pNdfa_pN_TRT_ACKO_base)[[20]]
summary(lme_pNdfa_pN_ACKO_base <- lme(pN_ACKO_89 ~ pNdfa_u_ACKO_89,
	random=~1 | Tree_ACKO_89))
tab_ACKO <- summary(lme_pNdfa_pN_ACKO_base)[[20]]

pN_MOFA_89 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_04) & 
		dat$Use_growth_04 == 1,]$foliar_N_mg_g_04,
	dat[dat$Species=="MOFA" & !is.na(dat$Ndfa_u_06) & 
		dat$Use_growth_06 == 1,]$foliar_N_mg_g_06
	)
#summary(lme_pNdfa_pN_TRT_MOFA_base <- lme(pN_MOFA_89 ~ pNdfa_u_MOFA_89 * Treatment_MOFA_89,
#	random=~1 | Tree_MOFA_89))
#tab_MOFA <- summary(lme_pNdfa_pN_TRT_MOFA_base)[[20]]
summary(lme_pNdfa_pN_MOFA_base <- lme(pN_MOFA_89 ~ pNdfa_u_MOFA_89,
	random=~1 | Tree_MOFA_89))
tab_MOFA <- summary(lme_pNdfa_pN_MOFA_base)[[20]]

# 2018
plot(V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="LN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="LN",]$foliar_N_mg_g_04,
	xaxt="n",#yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="MN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="red",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="HN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="red",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="PHN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="ACKO" & V$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="red",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="LN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="LN",]$foliar_N_mg_g_04,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="MN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="MN",]$foliar_N_mg_g_04,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="HN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="HN",]$foliar_N_mg_g_04,
	pch=4,col="orange",cex=1.5)
points(V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="PHN",]$Ndfa_u_04,
	V[V$Use_growth_04>=0.5 & V$Species=="MOFA" & V$Treatment=="PHN",]$foliar_N_mg_g_04,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
axis(side=1, c(0,50,100), labels=TRUE, padj=-.5)
mtext("l",at=-50,padj=1.5)

# 2019
plot(V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="LN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="LN",]$foliar_N_mg_g_06,
	xaxt="n",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="MN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="MN",]$foliar_N_mg_g_06,
	pch=3,col="red",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="HN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="HN",]$foliar_N_mg_g_06,
	pch=4,col="red",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="PHN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="ACKO" & V$Treatment=="PHN",]$foliar_N_mg_g_06,
	pch=5,col="red",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="LN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="LN",]$foliar_N_mg_g_06,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="MN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="MN",]$foliar_N_mg_g_06,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="HN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="HN",]$foliar_N_mg_g_06,
	pch=4,col="orange",cex=1.5)
points(V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="PHN",]$Ndfa_u_06,
	V[V$Use_growth_06>=0.5 & V$Species=="MOFA" & V$Treatment=="PHN",]$foliar_N_mg_g_06,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
axis(side=1, c(0,50,100), labels=TRUE, padj=-.5)
mtext("m",at=-50,padj=1.5)

plot.new()
points(c(.3,.3),c(.7,.6),pch=c(2,0),col=c("red","orange"),cex=1,xpd=TRUE)
text(.3,0.7,"Acacia",pos=4,col="red",font=3,xpd=TRUE,cex=1.3)
text(.3,0.6,"Morella",pos=4,col="orange",font=3,xpd=TRUE,cex=1.3)
text(.3,0.65,"l-n",pos=2,font=1,xpd=TRUE,cex=1.3)

points(c(0.1,0.1,0.1,0.1),c(0.3,0.2,0.1,0),pch=c(2,3,4,5),col="black",cex=1,xpd=TRUE)
text(c(0.2,0.2,0.2,0.2),c(0.3,0.2,0.1,0),c("Control","+10","+15","+15+P"),
	pos=4,col="black",xpd=TRUE,cex=1.3)

# All years
plot(0:100,tab_ACKO[1,1] + tab_ACKO[2,1]*(0:100),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",
	type="l",col="red",lty=1)
points(0:100,tab_MOFA[1,1] + tab_MOFA[2,1]*(0:100),
	type="l",col="orange",lty=1)
abline(v=c(0,100),col="gray")
text(40,yl[2]*0.9,c("ns"),col="red",cex=1.3)
text(60,yl[2]*0.9,c("ns"),col="orange",cex=1.3)
axis(side=1, c(0,50,100), labels=TRUE, padj=-.5)
mtext("n",at=-50,padj=1.5)
mtext("Volcano",side=4,padj=0.5)

mtext(expression("%" ~ of ~ nitrogen ~ from ~ fixation ~ ("%" ~ N[dfa])), side=1, outer=T, at=0.52, padj=1.625)
mtext(expression(Foliar ~ N ~ (mg ~ g^-1)),side=2, outer=T, at=0.5, padj=-1.4)

dev.off()

#############################################################
#############################################################
#############################################################