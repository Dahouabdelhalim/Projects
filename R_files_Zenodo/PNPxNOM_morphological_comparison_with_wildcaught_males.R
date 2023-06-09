# comparisons with (wild-caught) parental species
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# load packages and function
library(here)
library(beeswarm)
library(ape)
library("Hmisc")
library(corrplot)
# http://www.sthda.com/english/wiki/correlation-matrix-formatting-and-visualization
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


# Python datasheet 
#@@@@@@@@@@@@@@@@@
Python<-read.table("Pundamilia Python.txt", sep="\\t", header=T)
# from van Rijssel et al. 2018 PRSB: https://doi.org/10.5061/dryad.vk23f
# subset to only nyerereis 
Nyerereis<-Python[Python$Species=="nyererei",]
# and relevant phenotypes
NyrereiPhenos<-Nyerereis[,c(1,2,6:19,52:53,55)]
# recode gaps 0=S, 0.5=M, 1=L
NyrereiPhenos$Gaps<-as.character(NyrereiPhenos$Gaps)
NyrereiPhenos$Gaps[NyrereiPhenos$Gaps=="M"]<-"0.5"
NyrereiPhenos$Gaps[NyrereiPhenos$Gaps=="L"]<-"1"
NyrereiPhenos$Gaps[NyrereiPhenos$Gaps=="S"]<-"0"
NyrereiPhenos$Gaps<-as.numeric(NyrereiPhenos$Gaps)
# adjust colnames
colnames(NyrereiPhenos)[1]<-"id"
colnames(NyrereiPhenos)[2]<-"Type"
colnames(NyrereiPhenos)[18]<-"Gap"
colnames(NyrereiPhenos)[19]<-"Shape"
# do within PNP size correction
NyrereiPhenos$BDcorr<-residuals(lm(log10(NyrereiPhenos$BD)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$HLcorr<-residuals(lm(log10(NyrereiPhenos$HL)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$HWcorr<-residuals(lm(log10(NyrereiPhenos$HW)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$LJLcorr<-residuals(lm(log10(NyrereiPhenos$LJL)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$LJWcorr<-residuals(lm(log10(NyrereiPhenos$LJW)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$SnLcorr<-residuals(lm(log10(NyrereiPhenos$SnL)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$SnWcorr<-residuals(lm(log10(NyrereiPhenos$SnW)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$ChDcorr<-residuals(lm(log10(NyrereiPhenos$ChD)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$PODcorr<-residuals(lm(log10(NyrereiPhenos$POD)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$IOWcorr<-residuals(lm(log10(NyrereiPhenos$IOW)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$EyLcorr<-residuals(lm(log10(NyrereiPhenos$EyL)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$EyDcorr<-residuals(lm(log10(NyrereiPhenos$EyD)~log10(NyrereiPhenos$SL), na.action=na.exclude))
NyrereiPhenos$POWcorr<-residuals(lm(log10(NyrereiPhenos$POW)~log10(NyrereiPhenos$SL), na.action=na.exclude))
# produce correlation matrix
Nyrereicors <- rcorr(as.matrix(NyrereiPhenos[,c(20:26,30,27:29,32,17:19)]), type ="pearson")
# rename columns
colnames(Nyrereicors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape" ) 
rownames(Nyrereicors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape" ) 
colnames(Nyrereicors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape" ) 
rownames(Nyrereicors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape" ) 
# plot it
{
  pdf("PNP_traitcorrelations.pdf", width=5, height=5, pointsize=11)
  par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(3,1,0))
  corrplot(Nyrereicors$r, type="upper", order="original", p.mat = Nyrereicors$P, sig.level = 0.05, 
           insig = "blank", tl.col = "black", tl.srt = 45, cl.ratio = 0.15, cl.align.text="r", cl.offset=-0.75, 
           title = "P. sp. 'nyererei-like' males", method="color", mar=c(0.5,0,0.75,0), cl.cex = 1)
  dev.off()
}


# Makobe datasheet
#@@@@@@@@@@@@@@@@@
Makobe<-read.table("Neochromis Makobe.txt", sep="\\t", header=T)
# from van Rijssel et al. 2018 PRSB: https://doi.org/10.5061/dryad.vk23f
# subset to only omnis 
Omnis<-Makobe[Makobe$Species=="omnicaeruleus",]
# subsample to same number as Nyrereis
OmniSubset<-sample(Omnis$Fishec, 98, replace = F)
Omnis<-subset(Omnis, Fishec%in%OmniSubset)
# and relevant phenotypes
OmniPhenos<-Omnis[,c(1,2,6:19,24:25,27)]
# recode gaps 0=S, 0.5=M, 1=L
OmniPhenos$Gaps<-as.character(OmniPhenos$Gaps)
OmniPhenos$Gaps[OmniPhenos$Gaps=="M"]<-"0.5"
OmniPhenos$Gaps[OmniPhenos$Gaps=="L"]<-"1"
OmniPhenos$Gaps[OmniPhenos$Gaps=="S"]<-"0"
OmniPhenos$Gaps<-as.numeric(OmniPhenos$Gaps)
# adjust colnames
colnames(OmniPhenos)[1]<-"id"
colnames(OmniPhenos)[2]<-"Type"
colnames(OmniPhenos)[18]<-"Gap"
colnames(OmniPhenos)[19]<-"Shape"
# do within PNP size correction
OmniPhenos$BDcorr<-residuals(lm(log10(OmniPhenos$BD)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$HLcorr<-residuals(lm(log10(OmniPhenos$HL)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$HWcorr<-residuals(lm(log10(OmniPhenos$HW)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$LJLcorr<-residuals(lm(log10(OmniPhenos$LJL)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$LJWcorr<-residuals(lm(log10(OmniPhenos$LJW)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$SnLcorr<-residuals(lm(log10(OmniPhenos$SnL)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$SnWcorr<-residuals(lm(log10(OmniPhenos$SnW)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$ChDcorr<-residuals(lm(log10(OmniPhenos$ChD)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$PODcorr<-residuals(lm(log10(OmniPhenos$POD)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$IOWcorr<-residuals(lm(log10(OmniPhenos$IOW)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$EyLcorr<-residuals(lm(log10(OmniPhenos$EyL)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$EyDcorr<-residuals(lm(log10(OmniPhenos$EyD)~log10(OmniPhenos$SL), na.action=na.exclude))
OmniPhenos$POWcorr<-residuals(lm(log10(OmniPhenos$POW)~log10(OmniPhenos$SL), na.action=na.exclude))
# produce correlation matrix
Omnicors <- rcorr(as.matrix(OmniPhenos[,c(20:26,30,27:29,32,17:19)]), type ="pearson")
# rename columns
colnames(Omnicors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape") 
rownames(Omnicors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape") 
colnames(Omnicors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape") 
rownames(Omnicors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows","Gap","T-Shape") 
# plot it
{
  pdf("NOM_traitcorrelations.pdf", width=5, height=5, pointsize=11)
  par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(3,1,0))
  corrplot(Omnicors$r, type="upper", order="original", p.mat = Omnicors$P, sig.level = 0.05, 
           insig = "blank", tl.col = "black", tl.srt = 45, cl.ratio = 0.15, cl.align.text="r", cl.offset=-0.75, 
           title = "N. omicaeruleus males", method="color", mar=c(0.5,0,0.75,0), cl.cex = 1)
  dev.off()
}


# read in F2 male data and subset to traits available in wild
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
allphenos<-read.table("allPhenos.csv", header=T, sep=",") # already size-corrected
# subset to males
F2phenosMales<-allphenos[allphenos$sex==1,]
str(F2phenosMales)
F2phenosMales$InnerRows<-as.numeric(F2phenosMales$InnerRows)
# produce correlation matrix
F2cors <- rcorr(as.matrix(F2phenosMales[,c(4:11,13:16,26,25,21)]), type ="pearson")
# rename columns
colnames(F2cors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD", "POD","IOW","POW", "T-Rows","Gap", "T-Shape") 
rownames(F2cors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD", "POD","IOW","POW", "T-Rows","Gap", "T-Shape") 
colnames(F2cors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD", "POD","IOW","POW", "T-Rows","Gap", "T-Shape") 
rownames(F2cors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD", "POD","IOW","POW", "T-Rows","Gap", "T-Shape") 
# plot it
{
  pdf("F2s_traitcorrelations.pdf", width=5, height=5, pointsize=11)
  par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(3,1,0))
  corrplot(F2cors$r, type="upper", order="original", p.mat = F2cors$P, sig.level = 0.05, 
           insig = "blank", tl.col = "black", tl.srt = 45, cl.ratio = 0.15, cl.align.text="r", cl.offset=-0.75, 
           title = "F2 males", method="color", mar=c(0.5,0,0.75,0), cl.cex = 1)
  dev.off()
}



# test the three matrices against each other (each with and without teeth)
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# Fs vs PNP
mantel.test(m1=F2cors$r, m2=Nyrereicors$r, nperm = 999, graph = FALSE, alternative = "two.sided")
#$z.stat 4.389994, p 0.001
mantel.test(m1=F2cors$r[1:13, 1:13], m2=Nyrereicors$r[1:13, 1:13], nperm = 999, graph = FALSE, alternative = "two.sided")
# $z.stat 4.190753, p 0.001

# F2s vs NOM
mantel.test(m1=F2cors$r, m2=Omnicors$r, nperm = 999, graph = FALSE, alternative = "two.sided")
#$z.stat 6.535013, p 0.001
mantel.test(m1=F2cors$r[1:13, 1:13], m2=Omnicors$r[1:13, 1:13], nperm = 999, graph = FALSE, alternative = "two.sided")
# $z.stat 6.261003, p 0.003

# the two parental species
mantel.test(m1=Nyrereicors$r, m2=Omnicors$r, nperm = 999, graph = FALSE, alternative = "two.sided")
#$z.stat 6.122101, p 0.001
mantel.test(m1=Nyrereicors$r[1:13, 1:13], m2=Omnicors$r[1:13, 1:13], nperm = 999, graph = FALSE, alternative = "two.sided")
#$z.stat  6.07037, p 0.001



# combine the two parental species for a correlation matrix  
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
Parentals<-rbind(OmniPhenos[,1:19], NyrereiPhenos[,1:19])
# size correction combined
Parentals$BDcorr<-residuals(lm(log10(Parentals$BD)~log10(Parentals$SL), na.action=na.exclude))
Parentals$HLcorr<-residuals(lm(log10(Parentals$HL)~log10(Parentals$SL), na.action=na.exclude))
Parentals$HWcorr<-residuals(lm(log10(Parentals$HW)~log10(Parentals$SL), na.action=na.exclude))
Parentals$LJLcorr<-residuals(lm(log10(Parentals$LJL)~log10(Parentals$SL), na.action=na.exclude))
Parentals$LJWcorr<-residuals(lm(log10(Parentals$LJW)~log10(Parentals$SL), na.action=na.exclude))
Parentals$SnLcorr<-residuals(lm(log10(Parentals$SnL)~log10(Parentals$SL), na.action=na.exclude))
Parentals$SnWcorr<-residuals(lm(log10(Parentals$SnW)~log10(Parentals$SL), na.action=na.exclude))
Parentals$ChDcorr<-residuals(lm(log10(Parentals$ChD)~log10(Parentals$SL), na.action=na.exclude))
Parentals$PODcorr<-residuals(lm(log10(Parentals$POD)~log10(Parentals$SL), na.action=na.exclude))
Parentals$IOWcorr<-residuals(lm(log10(Parentals$IOW)~log10(Parentals$SL), na.action=na.exclude))
Parentals$EyLcorr<-residuals(lm(log10(Parentals$EyL)~log10(Parentals$SL), na.action=na.exclude))
Parentals$EyDcorr<-residuals(lm(log10(Parentals$EyD)~log10(Parentals$SL), na.action=na.exclude))
Parentals$POWcorr<-residuals(lm(log10(Parentals$POW)~log10(Parentals$SL), na.action=na.exclude))
# produce correlation matrix
Parentalcors <- rcorr(as.matrix(Parentals[,c(20:26, 30, 27:29, 32, 17:19)]), type ="pearson")
# rename columns
colnames(Parentalcors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows", "Gap","T-Shape") 
rownames(Parentalcors$r)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows", "Gap","T-Shape") 
colnames(Parentalcors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows", "Gap","T-Shape") 
rownames(Parentalcors$P)<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","ChD","POD","IOW","POW","T-Rows", "Gap","T-Shape") 
# plot it
{
  pdf("Parentals_traitcorrelations.pdf", width=5, height=5, pointsize=11)
  par(mfrow=c(1,1), mar=c(0,0,0,0), oma=c(0,0,0,0), mgp=c(3,1,0))
  corrplot(Parentalcors$r, type="upper", order="original", p.mat = Parentalcors$P, sig.level = 0.1, 
           insig = "blank", tl.col = "black", tl.srt = 45, cl.ratio = 0.15, cl.align.text="r", cl.offset=-0.75, 
           title = "Parental species males", method="color", mar=c(0.5,0,0.75,0), cl.cex = 1)
  dev.off()
}



# now combine them to plot trait distributions  including F2s (to also look for transgressive trait values)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# versions needed to do combined size correction
Morpho<-read.table("allPhenoMeans.csv", header=T, sep=",")
teeth<-read.table("allPhenos.csv", header = T, sep=",")[c(1,21:26)]
Allphenos<-merge(Morpho, teeth, by="id", all=F) # only those with both morpo and teeth data

# only F2 males
F2malephenos<-Allphenos[Allphenos$Sex=="male",]
F2malephenos$InnerRows<-as.numeric(F2malephenos$InnerRows)

# only traits available in wild sp
F2malephenos<-F2malephenos[,c(1,20,3:16,26,25,21)]
# adjust colnames
colnames(F2malephenos)<- c("id","Type","SL","BD","HL","HW","LJL","LJW","SnL","SnW","EyL","EyD","ChD","POD","IOW","POW", "Rows","Gap","Shape")

ALLthree<-rbind(NyrereiPhenos[,1:19], OmniPhenos[,1:19], F2malephenos)

# define colours and point symbols
factor(ALLthree$Type)
cols<-c("firebrick2", "dodgerblue", "black")[factor(ALLthree$Type)]

# plot (log) traits against (log) SL
{
pdf("PNPxNOM_phenos_morphology_wild_andF2_males_traitvsSL_check.pdf", width=16, height=12, pointsize=10)
par(mfrow=c(4,4))
plot(log10(ALLthree$SL), log10(ALLthree$BD), pch=19, col=cols, xlab="log10(SL)", ylab="log10(BD)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$HL), pch=19, col=cols, xlab="log10(SL)", ylab="log10(HL)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$HW), pch=19, col=cols, xlab="log10(SL)", ylab="log10(HW)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$LJL), pch=19, col=cols, xlab="log10(SL)", ylab="log10(LJL)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$LJW), pch=19, col=cols, xlab="log10(SL)", ylab="log10(LJW)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$SnL), pch=19, col=cols, xlab="log10(SL)", ylab="log10(SnL)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$SnW), pch=19, col=cols, xlab="log10(SL)", ylab="log10(SnW)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$EyL), pch=19, col=cols, xlab="log10(SL)", ylab="log10(EyL)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$EyD), pch=19, col=cols, xlab="log10(SL)", ylab="log10(EyD)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$ChD), pch=19, col=cols, xlab="log10(SL)", ylab="log10(ChD)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$POD), pch=19, col=cols, xlab="log10(SL)", ylab="log10(POD)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$IOW), pch=19, col=cols, xlab="log10(SL)", ylab="log10(IOW)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
plot(log10(ALLthree$SL), log10(ALLthree$POW), pch=19, col=cols, xlab="log10(SL)", ylab="log10(POW)", cex.lab=2)
legend("bottomright", legend=c("PNP", "NOM", "F2"), col=c("firebrick2", "dodgerblue", "black"), pch=19, bty = "o")
dev.off()
}

# perform correlation tests and homogeneity of slopes
summary(aov(log10(ALLthree$BD)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$HL)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$HW)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$LJL)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$LJW)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$SnL)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$SnW)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$EyL)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$EyD)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$ChD)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$POD)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$IOW)~log10(ALLthree$SL)*ALLthree$Type))
summary(aov(log10(ALLthree$POW)~log10(ALLthree$SL)*ALLthree$Type))

# almost all not only significant effect of Type but also interaction, so slopes mostly non-homogeneous
# to compare them directly, I still need tp combine them for size correction
# (but definitely keep only within F2 size correction for any analyses only within F2s (and especially QTLmapping analyses))

ALLthree$BDcorr<-residuals(lm(log10(ALLthree$BD)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$HLcorr<-residuals(lm(log10(ALLthree$HL)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$HWcorr<-residuals(lm(log10(ALLthree$HW)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$LJLcorr<-residuals(lm(log10(ALLthree$LJL)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$LJWcorr<-residuals(lm(log10(ALLthree$LJW)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$SnLcorr<-residuals(lm(log10(ALLthree$SnL)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$SnWcorr<-residuals(lm(log10(ALLthree$SnW)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$ChDcorr<-residuals(lm(log10(ALLthree$ChD)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$PODcorr<-residuals(lm(log10(ALLthree$POD)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$IOWcorr<-residuals(lm(log10(ALLthree$IOW)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$EyLcorr<-residuals(lm(log10(ALLthree$EyL)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$EyDcorr<-residuals(lm(log10(ALLthree$EyD)~log10(ALLthree$SL), na.action=na.exclude))
ALLthree$POWcorr<-residuals(lm(log10(ALLthree$POW)~log10(ALLthree$SL), na.action=na.exclude))

ALLthree<-droplevels(ALLthree)
ALLthree$Type<-as.character(ALLthree$Type)
ALLthree$Type[ALLthree$Type=="omnicaeruleus"]<-"NOM"
ALLthree$Type[ALLthree$Type=="nyererei"]<-"PNP"

ALLthree$Type <- factor(ALLthree$Type,levels=c("PNP","F2", "NOM"))


cols3<-c("firebrick2", "black", "dodgerblue")[factor(ALLthree$Type)]
{
  pdf("PNPxNOM_phenos_morphology_wild_andF2_males_byType.pdf", width=14, height=10, pointsize=10)
  par(mfrow=c(4,4), mar=c(3,5,2,1))#, oma=c(3,3,0,1))
  
  beeswarm(ALLthree$BDcorr~ALLthree$Type, ylab="BD", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$BDcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$BDcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$BDcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$BDcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$BDcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$BDcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$HLcorr~ALLthree$Type, ylab="HL", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$HLcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$HLcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$HLcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$HLcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$HLcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$HLcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$HWcorr~ALLthree$Type, ylab="HW", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$HWcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$HWcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$HWcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$HWcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$HWcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$HWcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$LJLcorr~ALLthree$Type, ylab="LJL", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$LJLcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$LJLcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$LJLcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$LJLcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$LJLcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$LJLcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$LJWcorr~ALLthree$Type, ylab="LJW", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$LJWcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$LJWcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$LJWcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$LJWcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$LJWcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$LJWcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$SnLcorr~ALLthree$Type, ylab="SnL", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$SnLcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$SnLcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$SnLcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$SnLcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$SnLcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$SnLcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$SnWcorr~ALLthree$Type, ylab="SnW", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$SnWcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$SnWcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$SnWcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$SnWcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$SnWcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$SnWcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$ChDcorr~ALLthree$Type, ylab="ChD", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$ChDcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$ChDcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$ChDcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$ChDcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$ChDcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$ChDcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$PODcorr~ALLthree$Type, ylab="POD", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$PODcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$PODcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$PODcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$PODcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$PODcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$PODcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$IOWcorr~ALLthree$Type, ylab="IOW", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$IOWcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$IOWcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$IOWcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$IOWcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$IOWcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$IOWcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$POWcorr~ALLthree$Type, ylab="POW", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$POWcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$POWcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$POWcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$POWcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$POWcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$POWcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$EyLcorr~ALLthree$Type, ylab="EyL", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$EyLcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$EyLcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$EyLcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$EyLcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$EyLcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$EyLcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  beeswarm(ALLthree$EyDcorr~ALLthree$Type, ylab="EyD", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$EyDcorr[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$EyDcorr[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$EyDcorr[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$EyDcorr[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$EyDcorr[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$EyDcorr[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  boxplot(ALLthree$Rows~ALLthree$Type, ylab="T. Rows", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$Rows[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$Rows[ALLthree$Type=="PNP"],na.rm = T), col="firebrick2", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$Rows[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$Rows[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$Rows[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$Rows[ALLthree$Type=="NOM"],na.rm = T), col="dodgerblue", lwd=2)
  
  boxplot(ALLthree$Gap~ALLthree$Type, ylab="Gap", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$Gap[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$Gap[ALLthree$Type=="PNP"],na.rm = T), col="firebrick2", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$Gap[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$Gap[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$Gap[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$Gap[ALLthree$Type=="NOM"],na.rm = T), col="dodgerblue", lwd=2)
  
  beeswarm(ALLthree$Shape~ALLthree$Type, ylab="T. Shape", xlab="", cex=0.75, cex.lab=2, cex.axis=2, pwcol=cols3, pch = 19)
  segments(x0 = 0.75, y0 = mean(ALLthree$Shape[ALLthree$Type=="PNP"], na.rm = T), x1 =  1.25, y1 = mean(ALLthree$Shape[ALLthree$Type=="PNP"],na.rm = T), col="black", lwd=2)
  segments(x0 = 1.75, y0 = mean(ALLthree$Shape[ALLthree$Type=="F2"], na.rm = T), x1 =  2.25, y1 = mean(ALLthree$Shape[ALLthree$Type=="F2"],na.rm = T), col="black", lwd=2)
  segments(x0 = 2.75, y0 = mean(ALLthree$Shape[ALLthree$Type=="NOM"], na.rm = T), x1 =  3.25, y1 = mean(ALLthree$Shape[ALLthree$Type=="NOM"],na.rm = T), col="black", lwd=2)
  
  dev.off()
}


