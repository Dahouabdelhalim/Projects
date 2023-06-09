# PNPxNOM QTLmapping - calculate PVE and plot Figures
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# load R/qtl and here
library(here)
library(qtl)

# read in data
# *************
PNPxNOMdat <- read.cross("csv", "", here("Analysis/GenoPheno_PNPxNOM_7more.csv"), map.function="kosambi")
# apply jittermap
PNPxNOMdat<-jittermap(PNPxNOMdat)
# calculate QTL genotype probabilities conditional on the available marker data. 
PNPxNOMdat <- calc.genoprob(PNPxNOMdat, step=1, error.prob=0.05, map.function="kosambi")
# subset to males and females (e.g. for colour mapping I need only males)
sex<-PNPxNOMdat$pheno$sex # 161 individuals
PNPxNOMdatM<-subset(PNPxNOMdat,ind=(sex==1)) #  132 males
PNPxNOMdatF<-subset(PNPxNOMdat,ind=(sex==0)) #  29 females

# load QTLmapping results (generated in "PNPxNOM_QTLmapping.R")
# ************************
# morphology; males & females together, with sex/family as additive covariate where significant
# sex and family: HL, HW, LJL, LJW, SnL, SnW, EyL, EyD, POW
load(file="outMorphoSexFam")
load(file="opermMorphoSexFam")
# sex only: BD, ChD, POD, IOW
load(file="outMorphoSex")
load(file="opermMorphoSex")
# family only: CPL, CPD
load(file="outMorphoFam")
load(file="opermMorphoFam")
# guts; males & females together, with family as additive covariate
load(file="outgutFam")
load(file="opermgutFam")
# teeth; males & females together, with sex/family as additive covariate where significant
# no significant effects of family or sex:  toothShapeL, toothShapeU, toothdensityU
load(file="outteeth")
load(file="opermteeth")
# significant (p<0.05) effect of sex and family:  Gap, InnerRows (tricuspid.tg. separately because binary!)
load(file="outteethSexFam")
load(file="opermteethSexFam")
# significant (p<0.05) effect of family only: toothdensityL
load(file="outteethFam")
load(file="opermteethFam")
# colour (binary!); only males, with family as additive covariate where significant #
# significant (p<0.05) family effects in: Rdorsalfin1, Rdorsalfin2, Rhead, Rdorsum1, Rdorsum2, Yflanks, Ygillcover, (No_stripes separately because (normal model!))
load(file="outcolourFam")
load(file="opermcolourFam")
# no family effects in: Rflanks1.tg., Rcheek.tg., Rgillcover.tg., Rnose, Yupperlip.tg., Ycheek.tg., Ypelvicfin.tg., Ynose
load(file="outcolour")
load(file="opermcolour")
# sex
load(file="outSex")
load(file="opermSex")



# combined summaries of all analyses for all traits
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

outALL<-c(outMorphoSexFam, outMorphoSex, outMorphoFam, outgutFam, outteeth, outteethSexFam, outteethFam, outcolourFam, outcolour,
          labels=c("SexFam", "Sex", "Fam", "Fam,", "none", "SexFam", "Fam", "Fam,", "none"))
opermALL<-cbind(opermMorphoSexFam,opermMorphoSex, opermMorphoFam, opermgutFam, opermteeth, opermteethSexFam, opermteethFam, opermcolourFam, opermcolour,
            labels=c("SexFam", "Sex", "Fam", "Fam,", "none", "SexFam", "Fam", "Fam,", "none"))

summaryALL_tabByChr<-summary(outALL, perms=opermALL, alpha=0.1, pvalues=TRUE, format="tabByChr", ci.function = "bayesint")



# get Bayesian confidence intervals and calculate PVE
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# bayesint already given in summary:
#                                       chr pos ci.low ci.high  lod   pval
#Rdorsalfin1.Fam, :     c2.loc43        2   43      34    56.0  4.55  0.009
#Rdorsalfin2.Fam, :     c2.loc43        2   43      34    63.3  4.55  0.006
#Rnose.none :           chr3_21045798   3   36.4    23    43    3.52  0.088
#LJL.SexFam :           chr4_4634552    4   13.7    6     36.4  4.18  0.019
#HL.SexFam :            chr13_11189376  13  24.9    14    28    4.19  0.028
#HL.SexFam :            c15.loc10       15  10      0     35    3.86  0.057
#lod.Fam, :             c16.loc35       16  35.0    12    41.0  5.04  0.006 #intestinelength.Fam
#BD.Sex :               chr20_8108045   20  32.2    16    42.0  4.81  0.002
#POD.Sex :              chr20_3616200   20  16.6    5     35.2  3.99  0.033
#toothDensityU.none :   chr22_22167437  22  43      22    50.0  4.04  0.041
#InnerRows.SexFam :     c22.loc48       22  48      33    53.8  4.97  0.001
#lod.Fam :              c22.loc28       22  28      0     32.0  4.83  0.008 #toothDensityL.Fam

# PVE =(1-10^((-2*LOD)/n))*100

#                      LOD  n

#Rdorsalfin1.Fam, :   4.55  132
length(which(!is.na(PNPxNOMdat$pheno$Rdorsalfin1)))
round((1-10^((-2*4.55)/132))*100,2) # 14.68

#Rdorsalfin2.Fam, :   4.55  132
length(which(!is.na(PNPxNOMdat$pheno$Rdorsalfin2)))
round((1-10^((-2*4.55)/132))*100,2) # 14.68

#Rnose.none :         3.52  132
length(which(!is.na(PNPxNOMdat$pheno$Rnose)))
round((1-10^((-2*3.52)/132))*100,2) # 11.56

#LJL.SexFam :         4.18  158
length(which(!is.na(PNPxNOMdat$pheno$LJL)))
round((1-10^((-2* 4.18)/158))*100,2) # 11.47

#HL.SexFam :          4.19  158
#HL.SexFam :          3.86  158
length(which(!is.na(PNPxNOMdat$pheno$HL)))
round((1-10^((-2* 4.19)/158))*100,2) # 11.5
round((1-10^((-2* 3.86)/158))*100,2) # 10.64

#lod.Fam, :           5.04  158   #intestine length
length(which(!is.na(PNPxNOMdat$pheno$gutlength)))
round((1-10^((-2* 5.04)/158))*100,2) # 13.66

#BD.Sex :             4.81  157
length(which(!is.na(PNPxNOMdat$pheno$BD)))
round((1-10^((-2* 4.81)/157))*100,2) # 13.16

#POD.Sex :            3.99 158
length(which(!is.na(PNPxNOMdat$pheno$POD)))
round((1-10^((-2* 3.99)/158))*100,2) # 10.98

#toothDensityU.none : 4.04  161
length(which(!is.na(PNPxNOMdat$pheno$toothDensityU)))
round((1-10^((-2* 4.04)/161))*100,2) # 10.91

#InnerRows.SexFam :   4.97  160
length(which(!is.na(PNPxNOMdat$pheno$InnerRows)))
round((1-10^((-2* 4.97)/160))*100,2) # 13.33

#lod.Fam :            4.83  161   #toothDensityL.Fam
length(which(!is.na(PNPxNOMdat$pheno$toothDensityL)))
round((1-10^((-2* 4.83)/161))*100,2) # 12.9



# QTL plot(Figure 2)
# @@@@@@@@@@@@@@@@@@

{
pdf("Output/PNPxNOM_QTLs.pdf", width=7, height=5, pointsize=14, bg = "white")
par(mfrow=c(6,1), mar=c(1,1,0,1), oma=c(3.25,3.25,0.5,3.25), mgp=c(1.5,0.7,0))

# MORPHO (*with sign. QTL)
# sex and family: HL*, HW, LJL*, LJW, SnL, SnW, EyL, EyD, POW
fcols<-c("#78c679", "grey60", "#31a354", "grey60", "grey60", "grey60", "grey60", "grey60", "grey60")
plot(outMorphoSexFam, lodcolumn = 1, col=fcols[1], main="", xlab="", bandcol = "grey100", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "grey95", ylim=c(0,6) )
abline(h=summary(opermMorphoSexFam, alpha=0.05)[1], lty=2, lwd=1, col=fcols[1])
rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)
axis(2, las=2, at=seq(0,6,2), cex.axis=1)
for(i in 2:9){
  plot(outMorphoSexFam, lodcolumn = i,col=fcols[i], add=T)
  abline(h=summary(opermMorphoSexFam, alpha=0.05)[3], lty=2, lwd=1, col=fcols[3])    
}
# sex only: BD*, ChD, POD*, IOW
fcols<-c("#c2e699", "grey60", "#006837", "grey60")
for(i in 1:4){
  plot(outMorphoSex, lodcolumn = i,col=fcols[i], add=T)
  abline(h=summary(opermMorphoSex, alpha=0.05)[c(1,3)], lty=2, lwd=1, col=fcols[c(1,3)])    
}
# family only: CPL, CPD
fcols<-c("grey60", "grey60")
for(i in 1:2){
  plot(outMorphoFam, lodcolumn = i, col=fcols[i],add=T)
}

# TEETH (*with sign. QTL + also plot anway) (ignore tricuspid)
# no significant effects of family or sex:  toothShapeL+, toothShapeU+, toothdensityU*
fcols<-c("grey60", "grey60", "#fbb4b9")
plot(outteeth, lodcolumn = 1, col=fcols[1], main="", xlab="", bandcol = "grey100", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect = "grey95", ylim=c(0,6))
#abline(h=summary(opermteeth, alpha=0.05)[1], lty=2, lwd=1, col=fcols[1])
rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)
axis(2, las=2, at=seq(0,6,2), cex.axis=1)
for(i in 2:3){
  plot(outteeth, lodcolumn = i, col=fcols[i], add=T)
  abline(h=summary(opermteeth, alpha=0.05)[3], lty=2, lwd=1, col=fcols[3])    
}
# significant (p<0.05) effect of sex and family: InnerRows*, Gap (tricuspid.tg. separately because binary!)
fcols<-c("#f768a1","grey60")
for(i in 1:2){
  plot(outteethSexFam, lodcolumn = i, col=fcols[i], add=T)
  abline(h=summary(opermteethSexFam, alpha=0.05)[1], lty=2, lwd=1, col=fcols[1])    
}
# significant (p<0.05) effect of family only: toothdensityL*
fcols<-c("#ae017e")
for(i in 1:1){
  plot(outteethFam, lodcolumn = i,col=fcols[i], add=T)
  abline(h=summary(opermteethFam, alpha=0.05)[2], lty=2, lwd=1, col=fcols[i])    
}

# INTESTINES, with family as additive covariate
fcols<-c("#2171b5")
plot(outgutFam, lodcolumn = 1, main="", col=fcols[1], xlab="", bandcol = "grey100", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect ="grey95", ylim=c(0,6))
abline(h=summary(opermgutFam, alpha=0.05)[1], lty=2, lwd=1, col=fcols[1])
rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)
axis(2, las=2, at=seq(0,6,2), cex.axis=1)

# COLOUR
  # no fam effect
  fcols<-c("grey60","grey60","grey60","#fdcc8a","grey60","grey60","grey60","grey60")
  plot(outcolour, lodcolumn = 1, main="", col=fcols[1], xlab="", bandcol = "grey100", lwd=2, xaxt="n", yaxt="n", alternate.chrid=T, bgrect ="grey95", ylim=c(0,6))
  rect(par("usr")[1], par("usr")[3], grconvertX(1,"ndc","user"), grconvertY(0,"ndc","user"), col="white", xpd=NA, border=NA)
  axis(2, las=2, at=seq(0,6,2), cex.axis=1)
   for(i in 2:8){
    plot(outcolour, lodcolumn = i,col=fcols[i], add=T)
    abline(h=summary(opermcolour, alpha=0.05)[4], lty=2, lwd=1, col=fcols[4])    
  }
  # with fameffect
  fcols<-c("#fc8d59","#e34a33","grey60","grey60","grey60","grey60","grey60")
  for(i in 1:7){
    plot(outcolourFam, lodcolumn = i,col=fcols[i], add=T)
    abline(h=summary(opermcolourFam, alpha=0.05)[c(1,2)], lty=2, lwd=1, col=fcols[c(1,2)])    
  }
  
  # SEX
  fcols<-c("#1c9099")
  plot(outSex, lodcolumn = 1, main="", col=fcols[1], xlab="", bandcol = "grey100", lwd=2, xaxt="n", yaxt="n", alternate.chrid=F, bgrect ="grey95", ylim=c(0,6))
  abline(h=summary(opermSex, alpha=0.05)[1], lty=2, lwd=1, col=fcols[1])
  axis(2, las=2, at=seq(0,6,2), cex.axis=1)
  
mtext("LG", 1, line=2.5, cex=0.75)
mtext("LOD score", 2, at=grconvertY(0.625,"ndc","user"), line=2, cex=0.75)
dev.off()
}



# effect plots (Figure 3) 
#@@@@@@@@@@@@@@@@@@@@@@@@

{
  pdf("PNPxNOM_QTLs_morphology_effectplots.pdf", width=7, height=10, pointsize=12)
  par(mfrow=c(7,4),mar=c(4,4,1,1),oma=c(0,0,0,0),mgp=c(3,1,0))
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$BD, ylab="BD", main="", marker="chr20_8108045", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$BD, ylab="BD", main="", mname1="chr20_8108045") 
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$LJL, ylab="LJL", main="", marker="chr4_4634552", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$LJL, ylab="LJL", main="", mname1="chr4_4634552") 
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$HL, ylab="HL", main="", marker="chr13_11189376", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$HL, ylab="HL", main="", mname1="chr13_11189376") 
  
  find.marker(cross=PNPxNOMdat, chr = 15, pos = 10 )# chr15_7913411
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$HL, ylab="HL", main="", marker="chr15_7913411", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$HL, ylab="HL", main="", mname1="chr15_7913411") 
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$POD, ylab="POD", main="", marker="chr20_3616200", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$POD, ylab="POD", main="", mname1="chr20_3616200") 
  
  find.marker(cross=PNPxNOMdat, chr = 16, pos = 35) # chr16_15654960
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$gutlength, ylab="Intestine Length", main="", marker="chr16_15654960", pch=19, col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$gutlength, ylab="Intestine Length", main="", mname1="chr16_15654960") 
  
  find.marker(cross=PNPxNOMdat, chr = 22, pos = 28) # chr22_7741136
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$toothDensityL, ylab="Tooth Density LJ", main="", marker="chr22_7741136", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  axis(2, las=2, at=seq(0,1,0.5))
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$toothDensityL, ylab="Tooth Density LJ", main="", mname1="chr22_7741136") 
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$toothDensityU, ylab="Tooth Density UJ", main="", marker="chr22_22167437", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  axis(2, las=2, at=seq(0,1,0.5))
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$toothDensityU, ylab="Tooth Density UJ", main="", mname1="chr22_22167437") 
  
  find.marker(cross=PNPxNOMdat, chr = 22, pos = 48) # "chr22_23367185"
  #plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$InnerRows, ylab="InnerRows", main="", marker="chr22_23367185", pch=19)
  #effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$InnerRows, ylab="InnerRows", main="", mname1="chr22_23367185") 
  # doesn't work; probably because there are no Zeros (only values between 2-5)
  PNPxNOMdat$pheno$newInnerRows<-PNPxNOMdat$pheno$InnerRows-2
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$newInnerRows, ylab="No. Inner Tooth Rows", main="", marker="chr22_23367185", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$newInnerRows, ylab="No. Inner Tooth Rows", main="", mname1="chr22_23367185") 
  # that solved it, but the y-axes are now wrong - will have to label them manually

  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$sex, ylab="sex", main="", marker="chr6_3411618", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  axis(2, las=2, at=seq(0,1,1)) 
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$sex, ylab="sex", main="", mname1="chr6_3411618")

  find.marker(cross=PNPxNOMdat, chr = 2, pos = 43) # chr2_18179005
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$Rdorsalfin1, ylab="Red dorsal fin 1", main="", marker="chr2_18179005", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  axis(2, las=2, at=seq(0,1,1))
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$Rdorsalfin1, ylab="Red dorsal fin 1", main="", mname1="chr2_18179005") 
  
  plotPXG(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$Rnose, ylab="Red nose", main="", marker="chr3_21045798", pch=19, yaxt="n", col = c("firebrick2", "purple", "dodgerblue"), infer = F)
  axis(2, las=2, at=seq(0,1,1))
  effectplot(PNPxNOMdat, pheno.col=PNPxNOMdat$pheno$Rnose, ylab="Red nose", main="", mname1="chr3_21045798") 
  
  dev.off()
}


