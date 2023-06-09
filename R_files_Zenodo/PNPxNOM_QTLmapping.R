# PNPxNOM - QTL mapping 
#@@@@@@@@@@@@@@@@@@@@@@

# load libraries
library(here)
library(qtl)

# read in data
PNPxNOMdat <- read.cross("csv", "", "GenoPheno.csv", map.function="kosambi")

# apply jittermap
PNPxNOMdat<-jittermap(PNPxNOMdat)


# data checking 
#@@@@@@@@@@@@@@

# overall summary
summary(PNPxNOMdat) 
# kind of cross
class(PNPxNOMdat) 
# number of individuals
nind(PNPxNOMdat) 
# number of phenotypes
nphe(PNPxNOMdat) 
# number of chromosomes
nchr(PNPxNOMdat) 
# total number of markers
totmar(PNPxNOMdat) 
# number of markers per chromosome
nmar(PNPxNOMdat) 
# mean, max, min number of markers and sd
mean(nmar(PNPxNOMdat)) 
max(nmar(PNPxNOMdat)) 
min(nmar(PNPxNOMdat)) 
sd(nmar(PNPxNOMdat))

# plot the linkage map
plotMap(PNPxNOMdat, alternate.chrid = T)

# map summary statistics
map<-pull.map(PNPxNOMdat)
summary.map(map)
# 1295 makers, 1237.1 cM, average spacing 1.0 cM

# plot phenotypic data distributions
par(mfrow=c(7,6))
for(i in 1:42)  plotPheno(PNPxNOMdat, pheno.col=i) 

# normality (in numeric traits)? (males+females)
# shapiro.test(PNPxNOMdat$pheno$[trait]) # mostly ok

# plot missing data (in genotypes, per individual)
par(mfrow=c(1,1))
plotMissing(PNPxNOMdat, reorder = T) 

# estimate the sex-averaged recombination fraction between all pairs of genetic markers.
PNPxNOMdat<-est.rf(PNPxNOMdat)
checkAlleles(PNPxNOMdat)
# no apparent problems

# plot the recombination fraction
plotRF(PNPxNOMdat, alternate.chrid=TRUE)
# purple:large LOD score or small recombination fraction, while yellow is the reverse 

# identification of genotyping errors
# large scores indicate likely genotyping errors; 'error LOD scores < 4 can probably be ignored'
PNPxNOMdat<- calc.errorlod(PNPxNOMdat, map.function="kosambi") 
top.errorlod(PNPxNOMdat, cutoff=5)

# plot a measure of the proportion of missing information in the genotype data.
plotInfo(PNPxNOMdat)
hist(nmissing(PNPxNOMdat, what="mar"), breaks=50)

# create table showing the observed numbers of individuals with each genotype at each marker,
# including P-values from chi-square tests for Mendelian segregation.
gt<-geno.table(PNPxNOMdat)
gt[ gt$P.value < 0.05, ]
# so still some distored loci in there, but not too many at 0.1 level


# preparing for mapping 
#@@@@@@@@@@@@@@@@@@@@@@

# calculate QTL genotype probabilities conditional on the available marker data. 
PNPxNOMdat <- calc.genoprob(PNPxNOMdat, step=1, error.prob=0.05, map.function="kosambi")

# subset to males and females (e.g. for colour mapping I need only males)
sex<-PNPxNOMdat$pheno$sex # 161 individuals
PNPxNOMdatM<-subset(PNPxNOMdat,ind=(sex==1)) #  132 males
PNPxNOMdatF<-subset(PNPxNOMdat,ind=(sex==0)) #  29 females


# QTL mapping of mapping sex 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@

# 161 individuals; 132 males, 29 females

# single QTL model (only em, hk and mr methods available for binary traits)
outSex<- scanone(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$sex, model="binary", method="em")
save(outSex, file=as.character(paste("outSex", sep="/")))
# display the maximum LODscore on each chromosome
summary(outSex)
# returns just the highest peak from output of scanone
max(outSex)  
# get a genome-wide LOD significance threshold.
opermSex<-scanone(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$sex, model="binary", method="em", n.perm=1000, verbose=T)
save(opermSex, file=as.character(paste("opermSex", sep="/")))
summary(outSex, perms=opermSex, alpha=0.1, pvalues=TRUE)

# plot it
plot(outSex[c(1,2,3)], xlab="", ylab="", main="", bandcol="lightgrey", 
     ylim=c(0,6), col="black", xaxt="n", yaxt="n", alternate.chrid=F, bgrect = "grey95")
abline(h=summary(opermSex, alpha=0.10)[1], lty=2, lwd=1, col="black")
abline(h=summary(opermSex, alpha=0.05)[1], lty=3, lwd=1, col="black")
axis(2, las=2, at=seq(0,6,2))
mtext("Linkage Group", 1, line=2.5, cex=1)
mtext("LOD score", 2, at=grconvertY(0.85,"ndc","user"), line=2, cex=1)

# effect plots
plotPXG(PNPxNOMdat, marker="chr6_3411618", pheno.col=PNPxNOMdat$pheno$sex, main = "", ylab = "sex")
effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$sex, mname1="chr6_3411618", main = "", ylab = "sex") 

# calculate PVE
length(which(!is.na(PNPxNOMdat$pheno$sex)))
round((1-10^((-2* 5.1)/161))*100,2) # 13.57


# check all other traits for effects of sex and family (covariates) 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# visually inspect relationship between covariate(s) and phenotype(s) 
ftext<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","EyD","ChD","POD","IOW","POW","CPL","CPD","gutlength","fattytissue",
         "toothshapeL","toothShapeU","toothDensityL","toothDensityU","Gap","InnerRows","tricuspid.tg.")
# plot it 
  par(mfrow=c(6,4))
  for(i in c(1:16,18:24)){  
    boxplot(PNPxNOMdat$pheno[,i] ~ PNPxNOMdat$pheno$sex, horizontal=TRUE, xlab=ftext[i], col=c("red","blue"), ylab="Sex (0=female, 1=male)")
  }
  
  par(mfrow=c(6,4))
  for(i in c(1:16,18:24)){  
    boxplot(PNPxNOMdat$pheno[,i] ~ PNPxNOMdat$pheno$sex*PNPxNOMdat$pheno$family, horizontal=TRUE, xlab=ftext[i], col=c("red","blue"), ylab="Sex*Fam")
  }

# statistical analysis
result<-list()
for(i in c(1:16,18:24)){  
  summary<-anova(aov(PNPxNOMdat$pheno[,i] ~ PNPxNOMdat$pheno$sex*PNPxNOMdat$pheno$family))
  result[[i]]<-summary
}
names(result)<-ftext

# significant (p<0.05) effect of sex and family:  HL, HW, LJL, LJW, SnL, SnW, EyL, EyD, POW, gap, innerrows, tricuspid
# significant (p<0.05) effect of sex only:        BD, ChD, POD, IOW, 
# significant (p<0.05) effect of family only:     CPL, CPD, gutlength, toothdensityL
# no significant effects of family or sex:        toothshapeL, toothshapeU, toothdensityU
# (only EyD signifiant interaction sex* family)


# check effect of fam in males only 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
ftext2<-c("BD","HL","HW","LJL","LJW","SnL","SnW","EyL","EyD","ChD","POD","IOW","POW","CPL","CPD","gutlength","fattytissue",
         "toothshapeL","toothShapeU","toothDensityL","toothDensityU","Gap","InnerRows","tricuspid.tg.", 
         "Rdorsalfin1","Rdorsalfin2","Rhead","Rdorsum1","Rdorsum2","Rflanks1.tg.","Rcheek.tg.", "Rgillcover.tg.","Rnose",
         "Yflanks","Yupperlip.tg.","Ycheek.tg.","Ygillcover", "Ypelvicfin.tg.","Ynose", "No_stripes")
# plot it
  par(mfrow=c(8,5))
  for(i in 1:40){  
    boxplot(PNPxNOMdatM$pheno[,i] ~ PNPxNOMdatM$pheno$family, horizontal=TRUE, xlab=ftext2[i], col=c("lightblue","blue"), ylab="Family")
  }

# statistical analysis
result2<-list()
for(i in 1:40){  
  summary<-anova(aov(PNPxNOMdatM$pheno[,i] ~ PNPxNOMdatM$pheno$family))
  result2[[i]]<-summary
}
names(result2)<-ftext2

# significant (p<0.05) family effects in: HL, LJL, LJW, SnL, SnW, POW, CPD, 
#                                         gutlength, gap, innerrows, tricuspid, 
#                                         Rdorsalfin1, Rdorsalfin2, Rhead, Rdorsum1, Rdorsum2, Yflanks, Ygillcover, No_stripes


# QTL mapping 
#@@@@@@@@@@@@

# external morphology; males & females together, with additive covariates where significant effect 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# sex and family: HL, HW, LJL, LJW, SnL, SnW, EyL, EyD, POW
#**********************************************************
outMorphoSexFam<- scanone(PNPxNOMdat, pheno.col=c("HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "POW"), model="normal", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
save(outMorphoSexFam, file=as.character(paste("outMorphoSexFam", sep="/")))
opermMorphoSexFam<-scanone(PNPxNOMdat, pheno.col=c("HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "POW"), model="normal", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family),  n.perm=1000, verbose=T)
save(opermMorphoSexFam, file=as.character(paste("opermMorphoSexFam", sep="/")))
summary(outMorphoSexFam, perms=opermMorphoSexFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_morphology_SexFam<-summary(outMorphoSexFam, perms=opermMorphoSexFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_morphology_SexFam, file = "PNPxNOM_QTLs_morphology_SexFam.csv")

# plot it
  ftext=c("HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "POW")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:11){
    x=0; # if(i==12){x=7}
    plot(outMorphoSexFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=11){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermMorphoSexFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermMorphoSexFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.70,"ndc","user"), line=2, cex=0.8)


# sex only: BD, ChD, POD, IOW
#*****************************
outMorphoSex<- scanone(PNPxNOMdat, pheno.col=c("BD", "ChD", "POD", "IOW"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$sex)
save(outMorphoSex, file=as.character(paste("outMorphoSex", sep="/")))
opermMorphoSex<-scanone(PNPxNOMdat, pheno.col=c("BD", "ChD", "POD", "IOW"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$sex,  n.perm=1000, verbose=T)
save(opermMorphoSex, file=as.character(paste("opermMorphoSex", sep="/")))
summary(outMorphoSex, perms=opermMorphoSex, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_morphology_Sex<-summary(outMorphoSex, perms=opermMorphoSex, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_morphology_Sex, file = "PNPxNOM_QTLs_morphology_Sex.csv")

# plot it
  ftext=c("BD", "ChD", "POD", "IOW")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:6){
    x=0; # if(i==12){x=7}
    plot(outMorphoSex[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=6){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermMorphoSex, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermMorphoSex, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.85,"ndc","user"), line=2, cex=0.8)


# family only: CPL, CPD
#***********************
outMorphoFam<- scanone(PNPxNOMdat, pheno.col=c("CPL", "CPD"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family)
save(outMorphoFam, file=as.character(paste("outMorphoFam", sep="/")))
opermMorphoFam<-scanone(PNPxNOMdat, pheno.col=c("CPL", "CPD"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family,  n.perm=1000, verbose=T)
save(opermMorphoFam, file=as.character(paste("opermMorphoFam", sep="/")))
summary(outMorphoFam, perms=opermMorphoFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_morphology_Fam<-summary(outMorphoFam, perms=opermMorphoFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_morphology_Fam, file = "PNPxNOM_QTLs_morphology_Fam.csv")

# plot it
  ftext=c("CPL", "CPD")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:4){
    x=0; # if(i==12){x=7}
    plot(outMorphoFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=4){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermMorphoFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermMorphoFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)


# intestines (guts); males & females together, with family as additive covariate 
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
outgutFam<- scanone(PNPxNOMdat, pheno.col=c("gutlength"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family)
save(outgutFam, file=as.character(paste("outgutFam", sep="/")))
opermgutFam<-scanone(PNPxNOMdat, pheno.col=c("gutlength"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family,  n.perm=1000, verbose=T)
save(opermgutFam, file=as.character(paste("opermgutFam", sep="/")))
summary(outgutFam, perms=opermgutFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_gut_Fam<-summary(outgutFam, perms=opermgutFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_gut_Fam, file = "PNPxNOM_QTLs_gut_Fam.csv")

# plot it
  ftext=c("gutlength")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:3){
    x=0; # if(i==12){x=7}
    plot(outgutFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=3){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermgutFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermgutFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)


# teeth; males & females together, with sex/family as additive covariate where significant 
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# no significant effects of family or sex:  toothShapeL, toothShapeU, toothdensityU
#**********************************************************************************
outteeth<- scanone(PNPxNOMdat, pheno.col=c("toothshapeL", "toothShapeU", "toothDensityU"), model="normal", method="em")
save(outteeth, file=as.character(paste("outteeth", sep="/")))
opermteeth<-scanone(PNPxNOMdat, pheno.col=c("toothshapeL", "toothShapeU", "toothDensityU"), model="normal", method="em",  n.perm=1000, verbose=T)
save(opermteeth, file=as.character(paste("opermteeth", sep="/")))
summary(outteeth, perms=opermteeth, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_teeth<-summary(outteeth, perms=opermteeth, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_teeth, file = "PNPxNOM_QTLs_teeth.csv")
# plot it
  ftext=c("toothShapeL", "toothShapeU", "toothDensityU")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:5){
    x=0; # if(i==12){x=7}
    plot(outteeth[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=5){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermteeth, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermteeth, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topleft", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)


# significant (p<0.05) effect of sex and family:  Gap, InnerRows (tricuspid.tg. separately because binary!)
#**********************************************************************************
outteethSexFam<- scanone(PNPxNOMdat, pheno.col=c("InnerRows", "Gap"), model="normal", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
save(outteethSexFam, file=as.character(paste("outteethSexFam", sep="/")))
opermteethSexFam<-scanone(PNPxNOMdat, pheno.col=c("InnerRows", "Gap"), model="normal", method="em",  addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), n.perm=1000, verbose=T)
save(opermteethSexFam, file=as.character(paste("opermteethSexFam", sep="/")))
summary(outteethSexFam, perms=opermteethSexFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_teeth_SexFam<-summary(outteethSexFam, perms=opermteethSexFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_teeth_SexFam, file = "PNPxNOM_QTLs_teeth_SexFam.csv")

# plot it
  ftext=c("InnerRows", "Gap")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:4){
    x=0; # if(i==12){x=7}
    plot(outteethSexFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=4){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermteethSexFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermteethSexFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topleft", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)

# significant (p<0.05) effect of family only: toothdensityL
#**********************************************************************************
outteethFam<- scanone(PNPxNOMdat, pheno.col=c("toothDensityL"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family)
save(outteethFam, file=as.character(paste("outteethFam", sep="/")))
opermteethFam<-scanone(PNPxNOMdat, pheno.col=c("toothDensityL"), model="normal", method="em",  addcovar = PNPxNOMdat$pheno$family, n.perm=1000, verbose=T)
save(opermteethFam, file=as.character(paste("opermteethFam", sep="/")))
summary(outteethFam, perms=opermteethFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_teeth_Fam<-summary(outteethFam, perms=opermteethFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_teeth_Fam, file = "PNPxNOM_QTLs_teeth_Fam.csv")

# plot it
  ftext=c("toothDensityL")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:3){
    x=0; # if(i==12){x=7}
    plot(outteethFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=3){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermteethFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermteethFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topleft", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)
  

# binary model for tricuspid (both covars) 
#*****************************************
outtricuspid<- scanone(PNPxNOMdat, pheno.col=c("tricuspid.tg."), model="binary", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
save(outtricuspid, file=as.character(paste("outtricuspid", sep="/")))
opermtricuspid<-scanone(PNPxNOMdat, pheno.col=c("tricuspid.tg."), model="binary", method="em",  addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), n.perm=1000, verbose=T)
save(opermtricuspid, file=as.character(paste("opermtricuspid", sep="/")))
summary(outtricuspid, perms=opermtricuspid, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_tricuspid<-summary(outtricuspid, perms=opermtricuspid, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_tricuspid, file = "PNPxNOM_QTLs_tricuspid.csv")

# plot it
  ftext=c("tricuspid.tg.")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:3){
    x=0; # if(i==12){x=7}
    plot(outtricuspid[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=3){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermtricuspid, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermtricuspid, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)


# colour (binary!); only males, with family as additive covariate where significant 
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# significant (p<0.05) family effects in: Rdorsalfin1, Rdorsalfin2, Rhead, Rdorsum1, Rdorsum2, Yflanks, Ygillcover
#*****************************************************************************************
outcolourFam<- scanone(PNPxNOMdatM, pheno.col=c("Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1", "Rdorsum2", "Yflanks", "Ygillcover"), model="binary", method="em", addcovar = PNPxNOMdatM$pheno$family)
save(outcolourFam, file=as.character(paste("outcolourFam", sep="/")))
opermcolourFam<-scanone(PNPxNOMdatM, pheno.col=c("Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1", "Rdorsum2", "Yflanks", "Ygillcover"), model="binary", method="em",  addcovar = PNPxNOMdatM$pheno$family, n.perm=1000, verbose=T)
save(opermcolourFam, file=as.character(paste("opermcolourFam", sep="/")))
summary(outcolourFam, perms=opermcolourFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_colour_Fam<-summary(outcolourFam, perms=opermcolourFam, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_colour_Fam, file = "PNPxNOM_QTLs_colour_Fam.csv")

# plot it
  ftext=c("Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1", "Rdorsum2", "Yflanks", "Ygillcover")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:9){
    x=0; # if(i==12){x=7}
    plot(outcolourFam[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=9){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermcolourFam, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermcolourFam, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)


# no family effects in: Rflanks1.tg., Rcheek.tg., Rgillcover.tg., Rnose, Yupperlip.tg., Ycheek.tg., Ypelvicfin.tg., Ynose
#*****************************************************************************************
outcolour<- scanone(PNPxNOMdatM, pheno.col=c("Rflanks1.tg.", "Rcheek.tg.", "Rgillcover.tg.", "Rnose", "Yupperlip.tg.", "Ycheek.tg.", "Ypelvicfin.tg.", "Ynose"), model="binary", method="em")
save(outcolour, file=as.character(paste("outcolour", sep="/")))
opermcolour<-scanone(PNPxNOMdatM, pheno.col=c("Rflanks1.tg.", "Rcheek.tg.", "Rgillcover.tg.", "Rnose", "Yupperlip.tg.", "Ycheek.tg.", "Ypelvicfin.tg.", "Ynose"), model="binary", method="em", n.perm=1000, verbose=T)
save(opermcolour, file=as.character(paste("opermcolour", sep="/")))
summary(outcolour, perms=opermcolour, alpha=0.1, pvalues=TRUE, format = "allpheno")
PNPxNOM_QTLs_colour<-summary(outcolour, perms=opermcolour, alpha=0.1, pvalues=TRUE, format = "allpheno")
write.csv(x = PNPxNOM_QTLs_colour, file = "PNPxNOM_QTLs_colour.csv")
# plot it
  ftext=c("Rflanks1.tg.", "Rcheek.tg.", "Rgillcover.tg.", "Rnose", "Yupperlip.tg.", "Ycheek.tg.", "Ypelvicfin.tg.", "Ynose")
  par(mfrow=c(15,1), mar=c(1,0,0,0), oma=c(3.25,3.25,0.5,0.25), mgp=c(1.5,0.7,0))
  for(i in 3:10){
    x=0; # if(i==12){x=7}
    plot(outcolour[c(1,2,i)], ylab="", main="", xlab="", bandcol = "grey96", 
         col="black", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "white", ylim=c(0,6+x))
    if(i!=10){rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)}
    #box()
    abline(h=summary(opermcolour, alpha=0.10)[i-2], lty=2, lwd=1)
    abline(h=summary(opermcolour, alpha=0.05)[i-2], lty=3, lwd=1, col="red")
    legend("topright", ftext[i-2], bty="n", x.intersp = 0)
    axis(2, las=2, at=seq(0,6,2))
  }
  mtext("Linkage Group", 1, line=2.5, cex=0.8)
  mtext("LOD score", 2, at=grconvertY(0.925,"ndc","user"), line=2, cex=0.8)

