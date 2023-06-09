# NSF FX data plotting

# Figure 5
# % Ndfa

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="Fig05.pdf",
	width=7.05,height=6)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,4),omi=c(.4,.5,.2,.2))
par(mai=c(0,0,0,0))

purefix_15N_AP <- 0.3663

##############################################
################## New York ##################
##############################################

source("NY_SNF_analysis_pNdfa_forpaper_base.R")
N <- dat
N_TRT <- rep(-1,nrow(N))
N_TRT[N$Treatment=="LN"] <- 0
N_TRT[N$Treatment=="MN"] <- 1
N_TRT[N$Treatment=="HN"] <- 2
N_TRT[N$Treatment=="PHN"] <- 3
N_jTRT <- jitter(N_TRT)

# 2017
plot(N_jTRT[!is.na(N$Ndfa_u_04) & N$Use_growth_04==1],
	N[!is.na(N$Ndfa_u_04) & N$Use_growth_04==1,]$Ndfa_u_04,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=2,xaxt="n",col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("a",at=-0.45,padj=1.5)
mtext("Year 3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2018
plot(N_jTRT[!is.na(N$Ndfa_u_07) & N$Use_growth_07==1],
	N[!is.na(N$Ndfa_u_07) & N$Use_growth_07==1,]$Ndfa_u_07,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=2,xaxt="n",yaxt="n",col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("b",at=-0.45,padj=1.5)
mtext("Year 4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019

plot(N_jTRT[!is.na(N$Ndfa_u_08) & N$Use_growth_08==1],
	N[!is.na(N$Ndfa_u_08) & N$Use_growth_08==1,]$Ndfa_u_08,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=2,xaxt="n",yaxt="n",col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("c",at=-0.45,padj=1.5)
mtext("Year 5",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# All years
plot(0:3,c(pNdfa_base_LN[1],pNdfa_base_MN[1],pNdfa_base_HN[1],pNdfa_base_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=17,xaxt="n",yaxt="n",col="red",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
segments(0,pNdfa_base_LN[2],0,pNdfa_base_LN[3],col="red")
segments(1,pNdfa_base_MN[2],1,pNdfa_base_MN[3],col="red")
segments(2,pNdfa_base_HN[2],2,pNdfa_base_HN[3],col="red")
segments(3,pNdfa_base_PHN[2],3,pNdfa_base_PHN[3],col="red")
text(0:3,rep(40,4),c("a","b","b","b"),col="red",cex=1.3)
mtext("d",at=-0.45,padj=1.5)
mtext("All years",at=1.5)
mtext("New York",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")


##############################################
################### Oregon ###################
##############################################

source("OR_SNF_analysis_pNdfa_forpaper_base.R")
O <- dat
O_TRT <- rep(-1,nrow(O))
O_TRT[O$Treatment=="LN"] <- 0
O_TRT[O$Treatment=="MN"] <- 1
O_TRT[O$Treatment=="HN"] <- 2
O_TRT[O$Treatment=="PHN"] <- 3
O_jTRT <- jitter(O_TRT)

# 2018 
plot(O_jTRT[!is.na(O$Ndfa_u_04) & O$Use_growth_04==1],
	O[!is.na(O$Ndfa_u_04) & O$Use_growth_04==1,]$Ndfa_u_04,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=0,xaxt="n",col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("e",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(O_jTRT[!is.na(O$Ndfa_u_05) & O$Use_growth_05==1],
	O[!is.na(O$Ndfa_u_05) & O$Use_growth_05==1,]$Ndfa_u_05,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=0,xaxt="n",yaxt="n",col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("f",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2020
plot(O_jTRT[!is.na(O$Ndfa_u_06) & O$Use_growth_06==1],
	O[!is.na(O$Ndfa_u_06) & O$Use_growth_06==1,]$Ndfa_u_06,
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=0,xaxt="n",yaxt="n",col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
mtext("g",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# All years
plot(0:3,c(pNdfa_base_LN[1],pNdfa_base_MN[1],pNdfa_base_HN[1],pNdfa_base_PHN[1]),
	xlim=c(-0.5,3.5),ylim=c(-3,103),
	pch=15,xaxt="n",yaxt="n",col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
segments(0,pNdfa_base_LN[2],0,pNdfa_base_LN[3],col="orange")
segments(1,pNdfa_base_MN[2],1,pNdfa_base_MN[3],col="orange")
segments(2,pNdfa_base_HN[2],2,pNdfa_base_HN[3],col="orange")
segments(3,pNdfa_base_PHN[2],3,pNdfa_base_PHN[3],col="orange")
text(0:3,rep(40,4),c("a","b","b","b"),col="orange",cex=1.3)
mtext("h",at=-0.45,padj=1.5)
mtext("Oregon",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Waiakea ###################
###############################################

source("HI_W_SNF_analysis_pNdfa_forpaper_base.R")
W <- dat
W_TRT <- rep(-1,nrow(W))
W_TRT[W$Treatment=="LN"] <- 0
W_TRT[W$Treatment=="MN"] <- 1
W_TRT[W$Treatment=="HN"] <- 2
W_TRT[W$Treatment=="PHN"] <- 3
W_jTRT <- jitter(W_TRT,amount=0.1)

# 2018
plot(W_jTRT[!is.na(W$Ndfa_u_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$Ndfa_u_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE",]$Ndfa_u_04,
	xlim=c(-0.5,3.5),ylim=c(-20,103),
	pch=2,xaxt="n",yaxt="n",col="red",cex=1.5)
points(W_jTRT[!is.na(W$Ndfa_u_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$Ndfa_u_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$Ndfa_u_04,
	pch=0,col="orange",cex=1.5)
axis(side=2,at=c(-20,0,20,40,60,80,100),labels=c("","0","","40","","80",""))
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
mtext("i",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(W_jTRT[!is.na(W$Ndfa_u_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$Ndfa_u_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE",]$Ndfa_u_06,
	xlim=c(-0.5,3.5),ylim=c(-20,103),
	pch=2,xaxt="n",yaxt="n",col="red",cex=1.5)
points(W_jTRT[!is.na(W$Ndfa_u_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$Ndfa_u_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ",]$Ndfa_u_06,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
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
plot((0:3)-0.15,c(pNdfa_base_LN_GLSE[1],pNdfa_base_MN_GLSE[1],
	pNdfa_base_HN_GLSE[1],pNdfa_base_PHN_GLSE[1]),
	xlim=c(-0.5,3.5),ylim=c(-20,103),
	pch=17,xaxt="n",yaxt="n",col="red",cex=1.5)
points((0:3)+0.15,c(pNdfa_base_LN_CAEQ[1],pNdfa_base_MN_CAEQ[1],
	pNdfa_base_HN_CAEQ[1],pNdfa_base_PHN_CAEQ[1]),
	pch=15,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
segments(0-0.15,pNdfa_base_LN_GLSE[2],0-0.15,pNdfa_base_LN_GLSE[3],col="red")
segments(1-0.15,pNdfa_base_MN_GLSE[2],1-0.15,pNdfa_base_MN_GLSE[3],col="red")
segments(2-0.15,pNdfa_base_HN_GLSE[2],2-0.15,pNdfa_base_HN_GLSE[3],col="red")
segments(3-0.15,pNdfa_base_PHN_GLSE[2],3-0.15,pNdfa_base_PHN_GLSE[3],col="red")
segments(0+0.15,pNdfa_base_LN_CAEQ[2],0+0.15,pNdfa_base_LN_CAEQ[3],col="orange")
segments(1+0.15,pNdfa_base_MN_CAEQ[2],1+0.15,pNdfa_base_MN_CAEQ[3],col="orange")
segments(2+0.15,pNdfa_base_HN_CAEQ[2],2+0.15,pNdfa_base_HN_CAEQ[3],col="orange")
segments(3+0.15,pNdfa_base_PHN_CAEQ[2],3+0.15,pNdfa_base_PHN_CAEQ[3],col="orange")
text((0:3)-0.2,rep(40,4),c("a","a","a","a"),col="red",cex=1.3)
text((0:3)+0.2,rep(40,4),c("a","ab","ab","b"),col="orange",cex=1.3)
mtext("k",at=-0.45,padj=1.5)
mtext("Waiakea",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

###############################################
################### Volcano ###################
###############################################

source("HI_V_SNF_analysis_pNdfa_forpaper_base.R")
V <- dat
V_TRT <- rep(-1,nrow(V))
V_TRT[V$Treatment=="LN"] <- 0
V_TRT[V$Treatment=="MN"] <- 1
V_TRT[V$Treatment=="HN"] <- 2
V_TRT[V$Treatment=="PHN"] <- 3
V_jTRT <- jitter(V_TRT,amount=0.1)

# 2018
plot(V_jTRT[!is.na(V$Ndfa_u_04) & V$Use_growth_04==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$Ndfa_u_04) & V$Use_growth_04==1 & V$Species=="ACKO",]$Ndfa_u_04,
	xlim=c(-0.5,3.5),ylim=c(-60,103),
	pch=2,xaxt="n",col="red",cex=1.5)
points(V_jTRT[!is.na(V$Ndfa_u_04) & V$Use_growth_04==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$Ndfa_u_04) & V$Use_growth_04==1 & V$Species=="MOFA",]$Ndfa_u_04,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
mtext("l",at=-0.45,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

# 2019
plot(V_jTRT[!is.na(V$Ndfa_u_06) & V$Use_growth_06==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$Ndfa_u_06) & V$Use_growth_06==1 & V$Species=="ACKO",]$Ndfa_u_06,
	xlim=c(-0.5,3.5),ylim=c(-60,103),
	pch=2,xaxt="n",yaxt="n",col="red",cex=1.5)
points(V_jTRT[!is.na(V$Ndfa_u_06) & V$Use_growth_06==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$Ndfa_u_06) & V$Use_growth_06==1 & V$Species=="MOFA",]$Ndfa_u_06,
	pch=0,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
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
plot((0:3)-0.15,c(pNdfa_base_LN_ACKO[1],pNdfa_base_MN_ACKO[1],
	pNdfa_base_HN_ACKO[1],pNdfa_base_PHN_ACKO[1]),
	xlim=c(-0.5,3.5),ylim=c(-60,103),
	pch=17,xaxt="n",yaxt="n",col="red",cex=1.5)
points((0:3)+0.15,c(pNdfa_base_LN_MOFA[1],pNdfa_base_MN_MOFA[1],
	pNdfa_base_HN_MOFA[1],pNdfa_base_PHN_MOFA[1]),
	pch=15,col="orange",cex=1.5)
segments(-.2,0,3.2,0,col="gray")
segments(-.2,100,3.2,100,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1,3), labels=c("10","15+P"), padj=-.5)
segments(0-0.15,pNdfa_base_LN_ACKO[2],0-0.15,pNdfa_base_LN_ACKO[3],col="red")
segments(1-0.15,pNdfa_base_MN_ACKO[2],1-0.15,pNdfa_base_MN_ACKO[3],col="red")
segments(2-0.15,pNdfa_base_HN_ACKO[2],2-0.15,pNdfa_base_HN_ACKO[3],col="red")
segments(3-0.15,pNdfa_base_PHN_ACKO[2],3-0.15,pNdfa_base_PHN_ACKO[3],col="red")
segments(0+0.15,pNdfa_base_LN_MOFA[2],0+0.15,pNdfa_base_LN_MOFA[3],col="orange")
segments(1+0.15,pNdfa_base_MN_MOFA[2],1+0.15,pNdfa_base_MN_MOFA[3],col="orange")
segments(2+0.15,pNdfa_base_HN_MOFA[2],2+0.15,pNdfa_base_HN_MOFA[3],col="orange")
segments(3+0.15,pNdfa_base_PHN_MOFA[2],3+0.15,pNdfa_base_PHN_MOFA[3],col="orange")
text((0:3)-0.2,rep(-10,4),c("a","b","ab","ab"),col="red",cex=1.3)
text((0:3)+0.2,rep(-10,4),c("a","a","a","a"),col="orange",cex=1.3)
mtext("n",at=-0.45,padj=1.5)
mtext("Volcano",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression("%" ~ of ~ nitrogen ~ from ~ fixation ~ ("%" ~ N[dfa])),side=2, outer=T, at=0.5, padj=-1.5)


dev.off()

#############################################################
#############################################################
#############################################################