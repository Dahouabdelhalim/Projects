###
# Scripts for Fig 3 in Cooke et al. 2017 AREES
###
# 3 April 2017
# author: NS Upham
###

# Part (a) phylogeny
###
setwd("/set/to/working/dir") 
source("DR_circleTree_functions.R")
library(ape); library(phytools); library(phyloch); library(viridis)

# subset data from mammal phylo to caribbean mammals
mamPhy<-read.beast("MamPhy_fullPosterior_BDvr_pcsFIXED_NDexp_MCC_target.tre")

toKeep<-read.table("caribbeanMam_spList_65sp.txt")
colnames(toKeep)<-c("tiplabel","sp","fam","ord")
toDrop<-setdiff(mamPhy$tip.label,as.vector(toKeep$tiplabel))
caribPhy<-ladderize(drop.tip2(mamPhy,toDrop))

# get ES tip data
####
#all mammals
mam_ES<-read.table("MamPhy_5911sp_tipES_NDexp.txt")
ES<-mam_ES$ES
names(ES)<-mam_ES$tiplabel

# caribbean mammals
cES_table<-read.table("caribbeanMam_65sp_ESvalues.txt")
cES<-cES_table$x
names(cES)<-rownames(cES_table)


# get the CARIBBEAN MAMMAL branch colors... relatve to those of ALL mammals....
cES_ordered <- cES[match(caribPhy$tip.label,names(cES))]
cES_anc <- ace(cES_ordered, caribPhy, method="pic", scaled=TRUE)$ace

# Match ancestors totree edges
    match.anc <- match(as.numeric(names(cES_anc)), caribPhy$edge[,2])[-1]

    # Assign rates to each internal node
    reconRate <- vector(mode="numeric", length=length(caribPhy$edge[,2]))
    reconRate[match.anc] <- cES_anc[2:64] #[2:7237]

    # Assign rates to tips
    tip.idx <- sort(caribPhy$tip.label, index.return=TRUE)$ix

    reconRate[match(tip.idx, caribPhy$edge[,2])] <- cES[sort(names(cES))]

    # Create colour palette
    reconColors <- reconRate
	
	cols<-viridis(100)
	range <- quantile(ES, seq(0,1, 0.01))[2:101]
    for (i in 1:100) {
        if (i==100) {range[i] <- range[i]+0.1}
        if (i==1) {reconColors[which(reconRate <= range[i])] <- cols[i] } 
        if (i > 1) {reconColors[which((reconRate <= range[i]) & (reconRate > range[i-1]))] <- cols[i]}
        }

# which are missing DNA?
missing_all<-read.table("MamPhy_FIN4_1813sp_missing_LIST.txt", header=FALSE)
sampled_carib<-as.vector(setdiff(caribPhy$tip.label,missing_all$V1))
missing_carib<-as.vector(setdiff(caribPhy$tip.label,sampled_carib))
tipCols<-rep("black",length(caribPhy$tip.label))
tipCols[match(missing_carib,caribPhy$tip.label)]<-"grey"

# change tip labels
caribPhy_plot<-caribPhy
match.tip<-match(caribPhy$tip.label, as.vector(toKeep$tiplabel))
caribPhy_plot$tip.label<-as.vector(toKeep$sp[match.tip])


# plot caribbean mammal tree
pdf(file="ES_phylo_caribbeanMammals_65species_viridis_tipCols.pdf",onefile=TRUE)
XX1=-120
XX2=120
plot2.phylo(caribPhy_plot, show.tip.label=TRUE, cex=0.4, label.offset=0.4, type="p", edge.width=3, no.margin=FALSE, root.edge=TRUE, edge.color=as.matrix(reconColors), tip.color=tipCols,x.lim=c(-12.62546, 114.96015))#,x.lim=c(XX1,XX2), y.lim=c(XX1,XX2))
HPDbars(caribPhy_plot,label="height_95%_HPD", col="black",lwd=1)
axisPhylo()
dev.off()

# ES density plot
x.tick <- quantile(ES, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(NA,round(x.tick,1)), side=1, line=1.3, cex=1, lwd=1, tck=-0.05, cex.axis=1.2, mgp=c(1,1,0))
dens.rate <- density(ES)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1.2, tck=-0.05, mgp=c(1,1,0))

pdf("ES_distribution_caribMammals_65sp_relativeToAllMammals_viridis.pdf",onefile=TRUE)

plot(density(cES), col="dark grey", main="", bty="n", xlab="Evolutionary distinctiveness (Mya)", ylab="",axes=F, xlim=range(ES))
polygon(density(cES), col="light grey", border="black", bty="n",main="")
x.tick <- quantile(ES, c(0.01,0.5,0.99,1))
axis(at=c(0,x.tick), labels=c(NA,round(x.tick,1)), side=1, line=1.3, cex=1, lwd=1, tck=-0.05, cex.axis=1.2, mgp=c(1,1,0))
dens.rate <- density(cES)$y
axis(at=c(min(dens.rate),0.45*max(dens.rate),0.9*max(dens.rate)), labels=c(0,0.45,0.9), side=2, cex=1, las=1, lwd=1, cex.axis=1.2, tck=-0.05, mgp=c(1,1,0))
seg.tick <- quantile(ES, c(0.01,0.5,0.99))
segments(seg.tick[[1]],0,seg.tick[[1]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")
segments(seg.tick[[2]],0,seg.tick[[2]],max(dens.rate)*1, lty=2, lwd=2,col="black")
segments(seg.tick[[3]],0,seg.tick[[3]],max(dens.rate)*0.5, lty=2, lwd=2,col="black")

color.legend(1.82, -0.01, 23.3, 0, legend=NULL, rect.col= cols, gradient="x", align="lt", cex=0.75, col="black", lwd=0.1, border="grey") #xl, yb, xr, yt
#       1%  50%  99% 100% 
#  NA  1.8  4.7 23.3 77.7 
dev.off()


# part (b) body mass histogram
#####
BM_table<-read.csv("bmTable_withAgeCats.csv", header=TRUE)

CubanInvasives<-read.table("Invasive_Caribbean_mamSp.txt", header=TRUE)
invasMass<-CubanInvasives$Mass
names(invasMass)<-CubanInvasives$Species

###
# SUBSET the data:

# extant bats
extantBats<-BM_table[which(BM_table$Order=="Chiroptera" & BM_table$Extirpated_WI==0),1:6]

# extinct bats
extinctBats<-BM_table[which(BM_table$Order=="Chiroptera" & BM_table$Extirpated_WI==1),1:6]

# extant non-volant
extantNV<-BM_table[which(BM_table$Order!="Chiroptera" & BM_table$Extirpated_WI==0),1:6]

# extinct non-volant
extinctNV<-BM_table[which(BM_table$Order!="Chiroptera" & BM_table$Extirpated_WI==1),1:6]


# colors
library(viridis)
cols<-c(grey(0.8,alpha=0.7),grey(0,alpha=1),grey(0.8,alpha=0.7))
colHerb<-"chartreuse4"
cexLeg<-1
ageCat_Cols<-c("deepskyblue2","darkorchid3","darkgoldenrod2")


# with AGE CATEGORIES
pdf(file="bodyMass_Fig_bothLog_ageCats_withHerb.pdf",onefile=TRUE, width=6,height=6)

#quartz(width=6,height=6)
layout(matrix(c(1:4), 2, 2, byrow = FALSE))
par(oma = c(5,4,5,0) + 0.1, mar = rep(0.5,4) + 0.1)

dat<-extantBats
hist(log(dat$Mass_FINAL), breaks=seq(0,5,(1/6)), xlab="",ylab="", col=cols[1], main="",ylim=c(0,8),xlim=c(0,5), yaxt="n",xaxt="n", border="black")#, xlim=c(0,65))#, xaxt="n",yaxt=)
mtext(side=3, "Bat species", font=2, line=1.5)
mtext(side=2, "Extant", font=2, line=2)
axis(side=2, at=c(0,2,4,6,8), labels=c(0,NA,4,NA,8))
axis(side=1, at=c(0:5), labels=FALSE)
text(x=1,y=7,labels=paste("n = ",length(na.omit(dat$Mass_FINAL))," / ",length(dat$Mass_FINAL),sep=""), font=3)
datHerb<-extantBats[which(extantBats$bat_diet==1),]
den<-10
hist(log(datHerb$Mass_FINAL), add=TRUE, breaks=seq(0,5,(1/6)), col=colHerb, xlab="",ylab="",border="black")#, density=den, angle=45,col=colHerb)#, xlim=c(0,65))#, xaxt="n",yaxt=)
legend(x=3,y=7.7,legend=c("native","herbivore"),fill=c(grey(0.8,alpha=0.7),colHerb),border="black", cex=cexLeg, bty="n")


dat<-extinctBats
preHol<-extinctBats[which(extinctBats$whenExtinct=="preHolocene"),"Mass_FINAL"] #4
Hol<-extinctBats[which(extinctBats$whenExtinct=="Holocene"),"Mass_FINAL"] # 4
Unk<-extinctBats[which(extinctBats$whenExtinct=="unknown"),"Mass_FINAL"] #4, but no mass data for any!!

hist(log(dat$Mass_FINAL), breaks=seq(0,5,(1/6)), xlab="Body mass (g)", col=ageCat_Cols[2], ylab="", main="",ylim=c(0,8),xlim=c(0,5), yaxt="n",border="black")#, xlim=c(0,65))
hist(log(preHol), add=TRUE, breaks=seq(0,5,(1/6)), col=ageCat_Cols[1], ylab="", main="",ylim=c(0,8),xlim=c(0,5), yaxt="n",border="black")#, xlim=c(0,65))

mtext(side=2, "Extinct", font=2, line=2)
axis(side=2, at=c(0,2,4,6,8), labels=c(0,NA,4,NA,8))
text(x=1,y=7,labels=paste("n = ",length(na.omit(dat$Mass_FINAL))," / ",length(dat$Mass_FINAL),sep=""), font=3)
legend(x=2.5,y=7.7,legend=c("pre-Holocene","Holocene","unknown"),fill=ageCat_Cols,border="black", cex=cexLeg, bty="n")

datHerb<-extinctBats[which(extinctBats$bat_diet==1),]
den<-10
hist(log(datHerb$Mass_FINAL), add=TRUE, breaks=seq(0,5,(1/6)), xlab="",ylab="", density=den, angle=45,col=colHerb, main="",ylim=c(0,6),xlim=c(0,5), yaxt="n",border="black")#, xlim=c(0,65))#, xaxt="n",yaxt=)



dat<-c(extantNV$Mass_FINAL,invasMass)
hist(log(dat),breaks=seq(0,15,0.5), xlab="",ylab="", col=cols[3],main="", ylim=c(0,8),xlim=c(0,15), xaxt="n", yaxt="n",border="black")
hist(log(invasMass), breaks=seq(0,15,0.5), add=TRUE, col=cols[2], xlab="",ylab="", main="",border=grey(0.3))

mtext(side=3, "Non-volant species", font=2, line=1.5)
axis(side=1, at=c(0,5,10,15), labels=FALSE)
axis(side=2, at=c(0,2,4,6,8), labels=FALSE)
text(x=3,y=7,labels=paste("n = ",length(na.omit(extantNV$Mass_FINAL))," / ",length(extantNV$Mass_FINAL),sep=""), font=3)
legend(x=9,y=7.7,legend=c("native","invasive"),fill=c(grey(0.8,alpha=0.7),grey(0)),border="black", cex=cexLeg, bty="n")


dat<-extinctNV
preHol<-extinctNV[which(extinctNV$whenExtinct=="preHolocene"),"Mass_FINAL"] #9
Hol<-extinctNV[which(extinctNV$whenExtinct=="Holocene"),"Mass_FINAL"] # 21, minus 1
Unk<-extinctNV[which(extinctNV$whenExtinct=="unknown"),"Mass_FINAL"] #29 minus 8

hist(log(dat$Mass_FINAL), breaks=seq(0,15,0.5), col=ageCat_Cols[3], xlab="",ylab="", main="", ylim=c(0,8), xlim=c(0,15), xaxt="n", yaxt="n",border="black")
hist(log(Hol), add=TRUE,breaks=seq(0,15,0.5), col=ageCat_Cols[2], xlab="",ylab="", main="", ylim=c(0,8), xlim=c(0,15), xaxt="n", yaxt="n",border="black")


axis(side=1, at=c(0,5,10,15), labels=TRUE)
axis(side=2, at=c(0,2,4,6,8), labels=FALSE)#c(0,NA,4,NA,8))
text(x=3,y=7,labels=paste("n = ",length(na.omit(dat$Mass_FINAL))," / ",length(dat$Mass_FINAL),sep=""), font=3)

title(main="", 
	  xlab = expression(bold('Mean body mass (ln grams)')),
      ylab = "Frequency",
      outer = TRUE, line = 2,cex.main=1.5,font.main=1,cex.axis=1.5,cex.lab=1.3,font.axis=1, font.lab=2)

dev.off()






