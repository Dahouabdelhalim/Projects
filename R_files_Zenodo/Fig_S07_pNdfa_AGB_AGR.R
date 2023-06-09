# NSF FX data plotting

# Figure S7
# % Ndfa against aboveground biomass growth rate

rm(list=ls())
library(emmeans)

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="FigS07.pdf",
	width=7.05,height=6)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,4),omi=c(.4,.5,.2,.2))
par(mai=c(0,0,0,0))

##############################################
################## New York ##################
##############################################

source("NY_SNF_analysis_pNdfa_forpaper_base.R")
N <- dat
yl <- c(min(c(N[N$Use_growth_07==1 & N$Species=="ROPS",]$AGR_AGB_2016_2017,
	N[N$Use_growth_07==1 & N$Species=="ROPS",]$AGR_AGB_2017_2018,
	N[N$Use_growth_08==1 & N$Species=="ROPS",]$AGR_AGB_2018_2019)),
	max(c(N[N$Use_growth_07==1 & N$Species=="ROPS",]$AGR_AGB_2016_2017,
	N[N$Use_growth_07==1 & N$Species=="ROPS",]$AGR_AGB_2017_2018,
	N[N$Use_growth_08==1 & N$Species=="ROPS",]$AGR_AGB_2018_2019)))

AGR_ROPS_AGB_789 <- AGR_ROPS_AGB_56789[(length(Tree_ROPS_56789)-length(Tree_ROPS_789)+1):length(Tree_ROPS_56789)]
AGB_ROPS_789 <- AGB_ROPS_56789[(length(Tree_ROPS_56789)-length(Tree_ROPS_789)+1):length(Tree_ROPS_56789)]

# Analysis section: 
# We tried two different models. The main one, which includes aboveground
# biomass as a covariate, is uncommented. The other adds treatment as a 
# fixed effect along with aboveground biomass as a covariate.

# None of them show significant effects of pNdfa on ln AGR of AGB for any
# species. The same models are coded for the other sites and species below, 
# though these notes aren't repeated.

#summary(lme_pNdfa_lnAGRAGB_base <- lme(log(AGR_ROPS_AGB_789) ~ pNdfa_u_789 * Treatment_ROPS_789 + AGB_ROPS_789,
#	random=~1 | Tree_ROPS_789))
#tab <- summary(lme_pNdfa_lnAGRAGB_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_base <- lme(log(AGR_ROPS_AGB_789) ~ pNdfa_u_789 + AGB_ROPS_789,
	random=~1 | Tree_ROPS_789))
tab <- summary(lme_pNdfa_lnAGRAGB_base)[[20]]

# 2017
plot(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="LN",]$AGR_AGB_2016_2017,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="MN",]$AGR_AGB_2016_2017,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="HN",]$AGR_AGB_2016_2017,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_04,
	N[N$Use_growth_04==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$AGR_AGB_2016_2017,
	pch=5,col="red",cex=1.5)
axis(side=2, c(1,3,10), labels=TRUE)
axis(side=2, c(20), labels=FALSE,lwd=0.25)
axis(side=2, seq(0.1,0.9,0.1), labels=FALSE,lwd=0.25)
axis(side=2, seq(1,9,1), labels=FALSE,lwd=0.25)
abline(v=c(0,100),col="gray")
mtext("a",at=-50,padj=1.5)
mtext("Year 3",at=mean(c(-60,103)))

# 2018
plot(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="LN",]$AGR_AGB_2017_2018,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_07,
	N[N$Use_growth_07==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="red",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("b",at=-50,padj=1.5)
mtext("Year 4",at=mean(c(-60,103)))

# 2019
plot(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="LN",]$AGR_AGB_2018_2019,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="red",cex=1.5)
points(N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$Ndfa_u_08,
	N[N$Use_growth_08==1 & N$Species=="ROPS" & N$Treatment=="PHN",]$AGR_AGB_2018_2019,
	pch=5,col="red",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("c",at=-50,padj=1.5)
mtext("Year 5",at=mean(c(-60,103)))

# All years
plot(0:100,exp(tab[1,1] + tab[2,1]*(0:100) + tab[3,1]*
	mean(N[N$Use_growth_07==1 & N$Species=="ROPS",]$AGB_est_kg_07,na.rm=TRUE)),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",log="y",
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
yl <- c(min(c(O[O$Use_growth_04==1 & O$Species=="ALRU",]$AGR_AGB_2017_2018,
	O[O$Use_growth_05==1 & O$Species=="ALRU",]$AGR_AGB_2018_2019,
	O[O$Use_growth_06==1 & O$Species=="ALRU",]$AGR_AGB_2019_2020)),
	max(c(O[O$Use_growth_04==1 & O$Species=="ALRU",]$AGR_AGB_2017_2018,
	O[O$Use_growth_05==1 & O$Species=="ALRU",]$AGR_AGB_2018_2019,
	O[O$Use_growth_06==1 & O$Species=="ALRU",]$AGR_AGB_2019_2020)))

AGR_ALRU_AGB_890 <- AGR_ALRU_AGB_6789[(length(Tree_ALRU_6789)-length(Tree_ALRU_890)+1):length(Tree_ALRU_6789)]
AGB_ALRU_890 <- AGB_ALRU_6789[(length(Tree_ALRU_6789)-length(Tree_ALRU_890)+1):length(Tree_ALRU_6789)]
#summary(lme_pNdfa_lnAGRAGB_TRT_base <- lme(log(AGR_ALRU_AGB_890) ~ pNdfa_u_890 * Treatment_ALRU_890 + AGB_ALRU_890,
#	random=~1 | Tree_ALRU_890))
#tab <- summary(lme_pNdfa_lnAGRAGB_TRT_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_base <- lme(log(AGR_ALRU_AGB_890) ~ pNdfa_u_890 + AGB_ALRU_890,
	random=~1 | Tree_ALRU_890))
tab <- summary(lme_pNdfa_lnAGRAGB_base)[[20]]

# 2018
plot(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="LN",]$AGR_AGB_2017_2018,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_04,
	O[O$Use_growth_04==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="orange",cex=1.5)
axis(side=2, c(0.5,1,2,5), labels=c("0.5","1","2","5"))
axis(side=2, seq(0.3,0.9,0.1), labels=FALSE,lwd=0.25)
axis(side=2, seq(1,9,1), labels=FALSE,lwd=0.25)
abline(v=c(0,100),col="gray")
mtext("e",at=-50,padj=1.5)

# 2019
plot(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="LN",]$AGR_AGB_2018_2019,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_05,
	O[O$Use_growth_05==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$AGR_AGB_2018_2019,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
mtext("f",at=-50,padj=1.5)

# 2020
plot(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="LN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="LN",]$AGR_AGB_2019_2020,
	xaxt="n",xlim=c(-60,103),log="y",ylim=yl,yaxt="n",
	pch=2,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="MN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="MN",]$AGR_AGB_2019_2020,
	pch=3,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="HN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="HN",]$AGR_AGB_2019_2020,
	pch=4,col="orange",cex=1.5)
points(O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$Ndfa_u_06,
	O[O$Use_growth_06==1 & O$Species=="ALRU" & O$Treatment=="PHN",]$AGR_AGB_2019_2020,
	pch=5,col="orange",cex=1.5)
abline(v=c(0,100),col="gray")
axis(side=1, c(0,50,100), labels=c("0","50","100  "), padj=-.5)
mtext("g",at=-50,padj=1.5)

# All years
plot(0:100,exp(tab[1,1] + tab[2,1]*(0:100) + tab[3,1]*
	mean(O[O$Use_growth_04==1 & O$Species=="ALRU",]$AGB_est_kg_04,na.rm=TRUE)),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",log="y",
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
yl <- c(min(c(W[W$Use_growth_04>=0.5 & W$Species!="PSCA",]$AGR_AGB_2017_2018,
	W[W$Use_growth_06>=0.5 & W$Species!="PSCA",]$AGR_AGB_2018_2019)),
	max(c(W[W$Use_growth_04>=0.5 & W$Species!="PSCA",]$AGR_AGB_2017_2018,
	W[W$Use_growth_06>=0.5 & W$Species!="PSCA",]$AGR_AGB_2018_2019)))

AGR_GLSE_AGB_89 <- AGR_GLSE_AGB_6789[(length(Tree_GLSE_6789)-length(Tree_GLSE_89)+1):length(Tree_GLSE_6789)]
AGB_GLSE_89 <- AGB_GLSE_6789[(length(Tree_GLSE_6789)-length(Tree_GLSE_89)+1):length(Tree_GLSE_6789)]
#summary(lme_pNdfa_lnAGRAGB_TRT_GLSE_base <- lme(log(AGR_GLSE_AGB_89) ~ pNdfa_u_GLSE_89 * Treatment_GLSE_89 + AGB_GLSE_89,
#	random=~1 | Tree_GLSE_89))
#tab_GLSE <- summary(lme_pNdfa_lnAGRAGB_TRT_GLSE_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_GLSE_base <- lme(log(AGR_GLSE_AGB_89) ~ pNdfa_u_GLSE_89 + AGB_GLSE_89,
	random=~1 | Tree_GLSE_89))
tab_GLSE <- summary(lme_pNdfa_lnAGRAGB_GLSE_base)[[20]]

AGR_CAEQ_AGB_89 <- AGR_CAEQ_AGB_6789[(length(Tree_CAEQ_6789)-length(Tree_CAEQ_89)+1):length(Tree_CAEQ_6789)]
AGB_CAEQ_89 <- AGB_CAEQ_6789[(length(Tree_CAEQ_6789)-length(Tree_CAEQ_89)+1):length(Tree_CAEQ_6789)]
#summary(lme_pNdfa_lnAGRAGB_TRT_CAEQ_base <- lme(log(AGR_CAEQ_AGB_89) ~ pNdfa_u_CAEQ_89 * Treatment_CAEQ_89 + AGB_CAEQ_89,
#	random=~1 | Tree_CAEQ_89))
#tab_CAEQ <- summary(lme_pNdfa_lnAGRAGB_TRT_CAEQ_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_CAEQ_base <- lme(log(AGR_CAEQ_AGB_89) ~ pNdfa_u_CAEQ_89 + AGB_CAEQ_89,
	random=~1 | Tree_CAEQ_89))
tab_CAEQ <- summary(lme_pNdfa_lnAGRAGB_CAEQ_base)[[20]]

# 2018
plot(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$AGR_AGB_2017_2018,
	xaxt="n",log="y",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="red",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$AGR_AGB_2017_2018,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="orange",cex=1.5)
points(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Ndfa_u_04,
	W[W$Use_growth_04>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="orange",cex=1.5)
axis(side=2, at=c(0.01,0.1,1,10), labels=c("0.01","0.1","1","10"))
axis(side=2, at=seq(0.003,0.01,0.001), labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.01,0.1,0.01), labels=FALSE,lwd=0.25)
axis(side=2, at=seq(0.1,1,0.1), labels=FALSE,lwd=0.25)
axis(side=2, at=seq(1,10,1), labels=FALSE,lwd=0.25)
axis(side=2, at=seq(10,20,10), labels=FALSE,lwd=0.25)
abline(v=c(0,100),col="gray")
mtext("i",at=-50,padj=1.5)

# 2019
plot(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="LN",]$AGR_AGB_2018_2019,
	xaxt="n",log="y",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="GLSE" & W$Treatment=="PHN",]$AGR_AGB_2018_2019,
	pch=5,col="red",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="LN",]$AGR_AGB_2018_2019,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="orange",cex=1.5)
points(W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$Ndfa_u_06,
	W[W$Use_growth_06>=0.5 & W$Species=="CAEQ" & W$Treatment=="PHN",]$AGR_AGB_2018_2019,
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
plot(0:100,exp(tab_GLSE[1,1] + tab_GLSE[2,1]*(0:100) + tab_GLSE[3,1]*
	mean(W[W$Use_growth_04>=0.5 & W$Species=="GLSE",]$AGB_est_kg_04,na.rm=TRUE)),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",log="y",
	type="l",col="red",lty=1)
points(0:100,exp(tab_CAEQ[1,1] + tab_CAEQ[2,1]*(0:100) + tab_CAEQ[3,1]*
	mean(W[W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$AGB_est_kg_04,na.rm=TRUE)),
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
yl <- c(min(c(V[V$Use_growth_04==1 & V$Species!="DOVI",]$AGR_AGB_2017_2018,
	V[V$Use_growth_06==1 & V$Species!="DOVI",]$AGR_AGB_2018_2019)),
	max(c(V[V$Use_growth_04==1 & V$Species!="DOVI",]$AGR_AGB_2017_2018,
	V[V$Use_growth_06==1 & V$Species!="DOVI",]$AGR_AGB_2018_2019)))

AGR_ACKO_AGB_89 <- AGR_ACKO_AGB_6789[(length(Tree_ACKO_6789)-length(Tree_ACKO_89)+1):length(Tree_ACKO_6789)]
AGB_ACKO_89 <- AGB_ACKO_6789[(length(Tree_ACKO_6789)-length(Tree_ACKO_89)+1):length(Tree_ACKO_6789)]
#summary(lme_pNdfa_lnAGRAGB_TRT_ACKO_base <- lme(log(AGR_ACKO_AGB_89) ~ pNdfa_u_ACKO_89 * Treatment_ACKO_89 + AGB_ACKO_89,
#	random=~1 | Tree_ACKO_89))
#tab_ACKO <- summary(lme_pNdfa_lnAGRAGB_TRT_ACKO_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_ACKO_base <- lme(log(AGR_ACKO_AGB_89) ~ pNdfa_u_ACKO_89 + AGB_ACKO_89,
	random=~1 | Tree_ACKO_89))
tab_ACKO <- summary(lme_pNdfa_lnAGRAGB_ACKO_base)[[20]]

AGR_MOFA_AGB_89 <- AGR_MOFA_AGB_6789[(length(Tree_MOFA_6789)-length(Tree_MOFA_89)+1):length(Tree_MOFA_6789)]
AGB_MOFA_89 <- AGB_MOFA_6789[(length(Tree_MOFA_6789)-length(Tree_MOFA_89)+1):length(Tree_MOFA_6789)]
#summary(lme_pNdfa_lnAGRAGB_TRT_MOFA_base <- lme(log(AGR_MOFA_AGB_89) ~ pNdfa_u_MOFA_89 * Treatment_MOFA_89 + AGB_MOFA_89,
#	random=~1 | Tree_MOFA_89))
#tab_MOFA <- summary(lme_pNdfa_lnAGRAGB_TRT_MOFA_base)[[20]]
summary(lme_pNdfa_lnAGRAGB_MOFA_base <- lme(log(AGR_MOFA_AGB_89) ~ pNdfa_u_MOFA_89 + AGB_MOFA_89,
	random=~1 | Tree_MOFA_89))
tab_MOFA <- summary(lme_pNdfa_lnAGRAGB_MOFA_base)[[20]]

# 2018
plot(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="LN",]$AGR_AGB_2017_2018,
	xaxt="n",log="y",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="red",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="red",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="red",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="LN",]$AGR_AGB_2017_2018,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="MN",]$AGR_AGB_2017_2018,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="HN",]$AGR_AGB_2017_2018,
	pch=4,col="orange",cex=1.5)
points(V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Ndfa_u_04,
	V[V$Use_growth_04==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$AGR_AGB_2017_2018,
	pch=5,col="orange",cex=1.5)
axis(side=2, c(.02,.2,2,20), labels=c("0.02","0.2","2","20"))
axis(side=2, seq(0.03,0.1,0.01), labels=FALSE,lwd=0.25)
axis(side=2, seq(0.3,1,0.1), labels=FALSE,lwd=0.25)
axis(side=2, seq(3,10,1), labels=FALSE,lwd=0.25)
axis(side=2, c(30), labels=FALSE,lwd=0.25)
abline(v=c(0,100),col="gray")
axis(side=1, c(0,50,100), labels=TRUE, padj=-.5)
mtext("l",at=-50,padj=1.5)

# 2019
plot(V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="LN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="LN",]$AGR_AGB_2018_2019,
	xaxt="n",log="y",yaxt="n",
	xlim=c(-60,103),
	ylim=yl,
	pch=2,col="red",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="MN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="red",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="HN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="red",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="ACKO" & V$Treatment=="PHN",]$AGR_AGB_2018_2019,
	pch=5,col="red",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="LN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="LN",]$AGR_AGB_2018_2019,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="MN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="MN",]$AGR_AGB_2018_2019,
	pch=3,col="orange",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="HN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="HN",]$AGR_AGB_2018_2019,
	pch=4,col="orange",cex=1.5)
points(V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$Ndfa_u_06,
	V[V$Use_growth_06==1 & V$Species=="MOFA" & V$Treatment=="PHN",]$AGR_AGB_2018_2019,
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
plot(0:100,exp(tab_ACKO[1,1] + tab_ACKO[2,1]*(0:100) + tab_ACKO[3,1]*
	mean(V[V$Use_growth_04==1 & V$Species=="ACKO",]$AGB_est_kg_04,na.rm=TRUE)),
	xlim=c(-60,103),ylim=yl,xaxt="n",yaxt="n",log="y",
	type="l",col="red",lty=1)
points(0:100,exp(tab_MOFA[1,1] + tab_MOFA[2,1]*(0:100) + tab_MOFA[3,1]*
	mean(V[V$Use_growth_04==1 & V$Species=="MOFA",]$AGB_est_kg_04,na.rm=TRUE)),
	type="l",col="orange",lty=1)
abline(v=c(0,100),col="gray")
text(40,yl[2]*0.9,c("ns"),col="red",cex=1.3)
text(60,yl[2]*0.9,c("ns"),col="orange",cex=1.3)
axis(side=1, c(0,50,100), labels=TRUE, padj=-.5)
mtext("n",at=-50,padj=1.5)
mtext("Volcano",side=4,padj=0.5)

mtext(expression("%" ~ of ~ nitrogen ~ from ~ fixation ~ ("%" ~ N[dfa])), side=1, outer=T, at=0.52, padj=1.625)
mtext(expression(Absolute ~ growth ~ rate ~ of ~ AGB ~ (kg ~ y^-1)),side=2, outer=T, at=0.5, padj=-1.4)

dev.off()

#############################################################
#############################################################
#############################################################