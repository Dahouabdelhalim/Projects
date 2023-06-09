# NSF FX data plotting

# Figure 4
# Foliar atom % 15N

rm(list=ls())

setwd("/Users/duncanmenge//Documents/Academia/Pubs/2022_Menge_etal_NSF_FX_main_fixation_paper/Data_and_code/")

#########################################################################
############################# Set up figure #############################
#########################################################################

pdf(file="Fig04.pdf",
  width=7.05,height=6)
# omi and mai coordinates are bottom, left, top, right
par(mfrow=c(4,3),omi=c(.4,.2,.2,.2),mai=c(0,.3,0,0))

purefix_15N_AP <- 0.3663

N <- read.csv("NY_FX_Size_FoliarCNIsotope_Data.csv")[1:66,1:166]
O <- read.csv("OR_FX_Size_FoliarCNIsotope_Data.csv")[1:64,1:101]
W <- read.csv("HI_W_FX_Size_FoliarCNIsotope_Data.csv")[1:108,1:111]
V <- read.csv("HI_V_FX_Size_FoliarCNIsotope_Data.csv")[1:96,1:111]

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

##############################################
################## New York ##################
##############################################

# Set unhealthy BENIs foliar d15N to NA
N[N$Species=="BENI" & N$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
N[N$Species=="BENI" & N$Use_growth_07==0,]$foliar_15N_AP_07 <- NA
N[N$Species=="BENI" & N$Use_growth_08==0,]$foliar_15N_AP_08 <- NA

# Set up columns for all years and set the reference to the paired tree
N$ref_15N_AP_u_04 <- NA
N[N$Species=="ROPS",]$ref_15N_AP_u_04 <- N[N$Species=="BENI",]$foliar_15N_AP_04
N$ref_15N_AP_u_07 <- NA
N[N$Species=="ROPS",]$ref_15N_AP_u_07 <- N[N$Species=="BENI",]$foliar_15N_AP_07
N$ref_15N_AP_u_08 <- NA
N[N$Species=="ROPS",]$ref_15N_AP_u_08 <- N[N$Species=="BENI",]$foliar_15N_AP_08

ylx <- max(c(N[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1 & N$Species=="ROPS",]$ref_15N_AP_u_04,
	N[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1 & N$Species=="ROPS",]$foliar_15N_AP_04),na.rm=TRUE)*1.02
plot(N_jTRT[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1],
	N[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1,]$ref_15N_AP_u_04,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(N_jTRT[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1 & N$Species=="ROPS"],
	N[!is.na(N$foliar_15N_AP_04) & N$Use_growth_04==1 & N$Species=="ROPS",]$foliar_15N_AP_04,
	col="red",pch=2,cex=1.5)
mtext("a",at=-0.35,padj=1.5)
mtext("Year 3",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1 & N$Species=="ROPS",]$ref_15N_AP_u_07,
	N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1 & N$Species=="ROPS",]$foliar_15N_AP_07),na.rm=TRUE)*1.02
plot(N_jTRT[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1],
	N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1,]$ref_15N_AP_u_07,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(N_jTRT[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1 & N$Species=="ROPS"],
	N[!is.na(N$foliar_15N_AP_07) & N$Use_growth_07==1 & N$Species=="ROPS",]$foliar_15N_AP_07,
	col="red",pch=2,cex=1.5)
mtext("b",at=-0.35,padj=1.5)
mtext("Year 4",at=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(N[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1 & N$Species=="ROPS",]$ref_15N_AP_u_08,
	N[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1 & N$Species=="ROPS",]$foliar_15N_AP_08),na.rm=TRUE)*1.02
plot(N_jTRT[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1],
	N[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1,]$ref_15N_AP_u_08,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(N_jTRT[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1 & N$Species=="ROPS"],
	N[!is.na(N$foliar_15N_AP_08) & N$Use_growth_08==1 & N$Species=="ROPS",]$foliar_15N_AP_08,
	col="red",pch=2,cex=1.5)
mtext("c",at=-0.35,padj=1.5)
mtext("Year 5",at=1.5)
mtext("New York",side=4,padj=0.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

##############################################
################### Oregon ###################
##############################################

# Set up columns for all years and set the reference to the paired tree
O$ref_15N_AP_u_04 <- NA
O[O$Species=="ALRU",]$ref_15N_AP_u_04 <- O[O$Species=="PSME",]$foliar_15N_AP_04
O$ref_15N_AP_u_05 <- NA
O[O$Species=="ALRU",]$ref_15N_AP_u_05 <- O[O$Species=="PSME",]$foliar_15N_AP_05
O$ref_15N_AP_u_06 <- NA
O[O$Species=="ALRU",]$ref_15N_AP_u_06 <- O[O$Species=="PSME",]$foliar_15N_AP_06

ylx <- max(c(O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1 & O$Species=="ALRU",]$ref_15N_AP_u_04,
	O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1 & O$Species=="ALRU",]$foliar_15N_AP_04),na.rm=TRUE)*1.02
plot(O_jTRT[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1],
	O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1,]$ref_15N_AP_u_04,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(O_jTRT[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1 & O$Species=="ALRU"],
	O[!is.na(O$foliar_15N_AP_04) & O$Use_growth_04==1 & O$Species=="ALRU",]$foliar_15N_AP_04,
	col="orange",pch=0,cex=1.5)
mtext("d",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(O[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1 & O$Species=="ALRU",]$ref_15N_AP_u_05,
	O[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1 & O$Species=="ALRU",]$foliar_15N_AP_05),na.rm=TRUE)*1.02
plot(O_jTRT[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1],
	O[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1,]$ref_15N_AP_u_05,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(O_jTRT[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1 & O$Species=="ALRU"],
	O[!is.na(O$foliar_15N_AP_05) & O$Use_growth_05==1 & O$Species=="ALRU",]$foliar_15N_AP_05,
	col="orange",pch=0,cex=1.5)
mtext("e",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(O[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1 & O$Species=="ALRU",]$ref_15N_AP_u_06,
	O[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1 & O$Species=="ALRU",]$foliar_15N_AP_06),na.rm=TRUE)*1.02
plot(O_jTRT[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1],
	O[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1,]$ref_15N_AP_u_06,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(O_jTRT[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1 & O$Species=="ALRU"],
	O[!is.na(O$foliar_15N_AP_06) & O$Use_growth_06==1 & O$Species=="ALRU",]$foliar_15N_AP_06,
	col="orange",pch=0,cex=1.5)
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1), labels=c("10"), padj=-.5)
axis(side=1, at=c(3), labels=c("15+P"),padj=-.5)
mtext("f",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")
mtext("Oregon",side=4,padj=0.5)

###############################################
################### Waiakea ###################
###############################################

# Set unhealthy PSCAs foliar d15N to NA
W[W$Species=="PSCA" & W$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
W[W$Species=="PSCA" & W$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Set up columns for all years and set the reference to the paired tree
W$ref_15N_AP_u_04 <- NA
W[W$Species=="GLSE",]$ref_15N_AP_u_04 <- W[W$Species=="PSCA",]$foliar_15N_AP_04
W[W$Species=="CAEQ",]$ref_15N_AP_u_04 <- W[W$Species=="PSCA",]$foliar_15N_AP_04
W$ref_15N_AP_u_06 <- NA
W[W$Species=="GLSE",]$ref_15N_AP_u_06 <- W[W$Species=="PSCA",]$foliar_15N_AP_06
W[W$Species=="CAEQ",]$ref_15N_AP_u_06 <- W[W$Species=="PSCA",]$foliar_15N_AP_06

# Waiakea only: Filter out the ones that weren't isotopically labeled
W[W$Species=="GLSE",]$ref_15N_AP_u_04[c(9,18,27)] <- NA
W[W$Species=="GLSE",]$ref_15N_AP_u_06[c(9,18,27)] <- NA
W[W$Species=="GLSE",]$foliar_15N_AP_04[c(9,18,27)] <- NA
W[W$Species=="GLSE",]$foliar_15N_AP_06[c(9,18,27)] <- NA
W[W$Species=="CAEQ",]$ref_15N_AP_u_04[c(9,18,27)] <- NA
W[W$Species=="CAEQ",]$ref_15N_AP_u_06[c(9,18,27)] <- NA
W[W$Species=="CAEQ",]$foliar_15N_AP_04[c(9,18,27)] <- NA
W[W$Species=="CAEQ",]$foliar_15N_AP_06[c(9,18,27)] <- NA

ylx <- max(c(W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & (W$Species=="GLSE" | W$Species=="CAEQ"),]$ref_15N_AP_u_04,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & (W$Species=="GLSE" | W$Species=="CAEQ"),]$foliar_15N_AP_04),na.rm=TRUE)*1.05
plot(W_jTRT[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE",]$ref_15N_AP_u_04,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(W_jTRT[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$ref_15N_AP_u_04,
	pch=1,col="blue",cex=1.5)
points(W_jTRT[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="GLSE",]$foliar_15N_AP_04,
	col="red",pch=2,cex=1.5)
points(W_jTRT[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$foliar_15N_AP_04) & W$Use_growth_04>=0.5 & W$Species=="CAEQ",]$foliar_15N_AP_04,
	col="orange",pch=0,cex=1.5)
mtext("g",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & (W$Species=="GLSE" | W$Species=="CAEQ"),]$ref_15N_AP_u_06,
	W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & (W$Species=="GLSE" | W$Species=="CAEQ"),]$foliar_15N_AP_06),na.rm=TRUE)*1.02
plot(W_jTRT[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE",]$ref_15N_AP_u_06,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
points(W_jTRT[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ",]$ref_15N_AP_u_06,
	pch=1,col="blue",cex=1.5)
points(W_jTRT[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE"]-0.2,
	W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="GLSE",]$foliar_15N_AP_06,
	col="red",pch=2,cex=1.5)
points(W_jTRT[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ"]+0.2,
	W[!is.na(W$foliar_15N_AP_06) & W$Use_growth_06>=0.5 & W$Species=="CAEQ",]$foliar_15N_AP_06,
	col="orange",pch=0,cex=1.5)
mtext("h",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()

lshift <- 0 #0.25
fcex <- 1.3 #0.8

points(c(0.3,0.3)-lshift,c(.8,.7),pch=c(2,1),col=c("red","blue"),cex=1,xpd=TRUE)
text(0.3-lshift,0.8,"Robinia",pos=4,col="red",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.7,"Betula",pos=4,col="blue",font=3,xpd=TRUE,cex=fcex)
text(0.28-lshift,0.75,"a-c",pos=2,xpd=TRUE,cex=1.3)

points(c(0.3,0.3)-lshift,c(.55,.45),pch=c(0,1),col=c("orange","blue"),cex=1,xpd=TRUE)
text(0.3-lshift,0.55,"Alnus",pos=4,col="orange",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.45,"Pseudotsuga",pos=4,col="blue",font=3,xpd=TRUE,cex=fcex)
text(0.28-lshift,0.5,"d-f",pos=2,xpd=TRUE,cex=1.3)

points(c(0.3,0.3,0.3)-lshift,c(.3,.2,.1),pch=c(2,0,1),col=c("red","orange","blue"),cex=1,xpd=TRUE)
text(0.3-lshift,0.3,"Gliricidia",pos=4,col="red",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.2,"Casuarina",pos=4,col="orange",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.1,"Psidium",pos=4,col="blue",font=3,xpd=TRUE,cex=fcex)
text(0.28-lshift,0.2,"g-h",pos=2,xpd=TRUE,cex=1.3)
mtext("Waiakea",side=4,padj=0.5)

###############################################
################### Volcano ###################
###############################################

# Set unhealthy DOVIs foliar d15N to NA
V[V$Species=="DOVI" & V$Use_growth_04==0,]$foliar_15N_AP_04 <- NA
V[V$Species=="DOVI" & V$Use_growth_06==0,]$foliar_15N_AP_06 <- NA

# Set up columns for all years and set the reference to the paired tree
V$ref_15N_AP_u_04 <- NA
V[V$Species=="ACKO",]$ref_15N_AP_u_04 <- V[V$Species=="DOVI",]$foliar_15N_AP_04
V[V$Species=="MOFA",]$ref_15N_AP_u_04 <- V[V$Species=="DOVI",]$foliar_15N_AP_04
V$ref_15N_AP_u_06 <- NA
V[V$Species=="ACKO",]$ref_15N_AP_u_06 <- V[V$Species=="DOVI",]$foliar_15N_AP_06
V[V$Species=="MOFA",]$ref_15N_AP_u_06 <- V[V$Species=="DOVI",]$foliar_15N_AP_06

ylx <- max(c(V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & (V$Species=="ACKO" | V$Species=="MOFA"),]$ref_15N_AP_u_04,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & (V$Species=="ACKO" | V$Species=="MOFA"),]$foliar_15N_AP_04),na.rm=TRUE)*1.02
plot(V_jTRT[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="ACKO",]$ref_15N_AP_u_04,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
axis(side=1, c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1), labels=c("10"), padj=-.5)
axis(side=1, at=c(3), labels=c("15+P"),padj=-.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="MOFA",]$ref_15N_AP_u_04,
	pch=1,col="blue",cex=1.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="ACKO",]$foliar_15N_AP_04,
	col="red",pch=2,cex=1.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$foliar_15N_AP_04) & V$Use_growth_04==1 & V$Species=="MOFA",]$foliar_15N_AP_04,
	col="orange",pch=0,cex=1.5)
mtext("i",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

ylx <- max(c(V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & (V$Species=="ACKO" | V$Species=="MOFA"),]$ref_15N_AP_u_06,
	V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & (V$Species=="ACKO" | V$Species=="MOFA"),]$foliar_15N_AP_06),na.rm=TRUE)*1.02
plot(V_jTRT[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="ACKO",]$ref_15N_AP_u_06,
	xlim=c(-0.5,3.5),ylim=c(purefix_15N_AP,ylx),
	pch=1,xaxt="n",col="blue",cex=1.5)
abline(h=purefix_15N_AP,col="gray")
axis(side=1, at=c(0,2), labels=c("C","15"), padj=-.5)
axis(side=1, at=c(1), labels=c("10"), padj=-.5)
axis(side=1, at=c(3), labels=c("15+P"),padj=-.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="MOFA",]$ref_15N_AP_u_06,
	pch=1,col="blue",cex=1.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="ACKO"]-0.2,
	V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="ACKO",]$foliar_15N_AP_06,
	col="red",pch=2,cex=1.5)
points(V_jTRT[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="MOFA"]+0.2,
	V[!is.na(V$foliar_15N_AP_06) & V$Use_growth_06==1 & V$Species=="MOFA",]$foliar_15N_AP_06,
	col="orange",pch=0,cex=1.5)
mtext("j",at=-0.35,padj=1.5)
abline(v=c(0.5,1.5,2.5),lty=c(2,2,1),col="gray")

plot.new()
points(c(0.3,0.3,0.3)-lshift,c(.6,.5,.4),pch=c(2,0,1),col=c("red","orange","blue"),cex=1,xpd=TRUE)
text(0.3-lshift,0.6,"Acacia",pos=4,col="red",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.5,"Morella",pos=4,col="orange",font=3,xpd=TRUE,cex=fcex)
text(0.3-lshift,0.4,"Dodonaea",pos=4,col="blue",font=3,xpd=TRUE,cex=fcex)
text(0.28-lshift,0.5,"i-j",pos=2,xpd=TRUE,cex=1.3)
mtext("Volcano",side=4,padj=0.5)


mtext("Treatment", side=1, outer=T, at=0.52, padj=2)
mtext(expression(Atom ~ "%" ~ ""^"15"*N),side=2, outer=T, at=0.5, padj=0)


dev.off()

#############################################################
#############################################################
#############################################################