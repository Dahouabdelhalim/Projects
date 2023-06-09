# PNPxNOM QTL mapping: muliple QTL models
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# load libraries
library(here)
library(qtl)

# read in data 
PNPxNOMdat <- read.cross("csv", "", here("Analysis/GenoPheno_PNPxNOM_7more.csv"), map.function="kosambi")
# apply jittermap
PNPxNOMdat<-jittermap(PNPxNOMdat)
# calculate QTL genotype probabilities conditional on the available marker data. 
PNPxNOMdat <- calc.genoprob(PNPxNOMdat, step=1, error.prob=0.05, map.function="kosambi")
# subset to males and females (e.g. for colour mapping I need only males)
sex<-PNPxNOMdat$pheno$sex # 161 individuals
PNPxNOMdatM<-subset(PNPxNOMdat,ind=(sex==1)) #  132 males
PNPxNOMdatF<-subset(PNPxNOMdat,ind=(sex==0)) #  29 females

# load scantwo permutations downloaded from cluster and combine the four permutation batches
# and calculate penalties 

# traits with covars sex and fam
permsHL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHL1")
permsHL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHL2")
permsHL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHL3")
permsHL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHL4")
permsHL1all<-c(permsHL1, permsHL2, permsHL3, permsHL4)
penHL<-calc.penalties(permsHL1all, alpha=c(0.05, 0.1))

permsHW1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHW1")
permsHW2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHW2")
permsHW3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHW3")
permsHW4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsHW4")
permsHW1all<-c(permsHW1, permsHW2, permsHW3, permsHW4)
penHW<-calc.penalties(permsHW1all, alpha=c(0.05, 0.1))

permsLJL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJL1")
permsLJL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJL2")
permsLJL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJL3")
permsLJL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJL4")
permsLJL1all<-c(permsLJL1, permsLJL2, permsLJL3, permsLJL4)
penLJL<-calc.penalties(permsLJL1all, alpha=c(0.05, 0.1))

permsLJW1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJW1")
permsLJW2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJW2")
permsLJW3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJW3")
permsLJW4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsLJW4")
permsLJW1all<-c(permsLJW1, permsLJW2, permsLJW3, permsLJW4)
penLJW<-calc.penalties(permsLJW1all, alpha=c(0.05, 0.1))

permsSnL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnL1")
permsSnL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnL2")
permsSnL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnL3")
permsSnL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnL4")
permsSnL1all<-c(permsSnL1, permsSnL2, permsSnL3, permsSnL4)
penSnL<-calc.penalties(permsSnL1all, alpha=c(0.05, 0.1))

permsSnW1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnW1")
permsSnW2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnW2")
permsSnW3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnW3")
permsSnW4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsSnW4")
permsSnW1all<-c(permsSnW1, permsSnW2, permsSnW3, permsSnW4)
penSnW<-calc.penalties(permsSnW1all, alpha=c(0.05, 0.1))

permsEyL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyL1")
permsEyL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyL2")
permsEyL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyL3")
permsEyL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyL4")
permsEyL1all<-c(permsEyL1, permsEyL2, permsEyL3, permsEyL4)
penEyL<-calc.penalties(permsEyL1all, alpha=c(0.05, 0.1))

permsEyD1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyD1")
permsEyD2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyD2")
permsEyD3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyD3")
permsEyD4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsEyD4")
permsEyD1all<-c(permsEyD1, permsEyD2, permsEyD3, permsEyD4)
penEyD<-calc.penalties(permsEyD1all, alpha=c(0.05, 0.1))

permsPOW1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOW1")
permsPOW2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOW2")
permsPOW3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOW3")
permsPOW4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOW4")
permsPOW1all<-c(permsPOW1, permsPOW2, permsPOW3, permsPOW4)
penPOW<-calc.penalties(permsPOW1all, alpha=c(0.05, 0.1))

permsInnerRows1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsInnerRows1")
permsInnerRows2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsInnerRows2")
permsInnerRows3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsInnerRows3")
permsInnerRows4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsInnerRows4")
permsInnerRows1all<-c(permsInnerRows1, permsInnerRows2, permsInnerRows3, permsInnerRows4)
penInnerRows<-calc.penalties(permsInnerRows1all, alpha=c(0.05, 0.1))

permsGap1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsGap1")
permsGap2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsGap2")
permsGap3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsGap3")
permsGap4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsGap4")
permsGap1all<-c(permsGap1, permsGap2, permsGap3, permsGap4)
penGap<-calc.penalties(permsGap1all, alpha=c(0.05, 0.1))

# traits with covar sex
permsBD1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsBD1")
permsBD2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsBD2")
permsBD3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsBD3")
permsBD4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsBD4")
permsBD1all<-c(permsBD1, permsBD2, permsBD3, permsBD4)
penBD<-calc.penalties(permsBD1all, alpha=c(0.05, 0.1))

permsChD1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsChD1")
permsChD2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsChD2")
permsChD3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsChD3")
permsChD4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsChD4")
permsChD1all<-c(permsChD1, permsChD2, permsChD3, permsChD4)
penChD<-calc.penalties(permsChD1all, alpha=c(0.05, 0.1))

permsPOD1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOD1")
permsPOD2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOD2")
permsPOD3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOD3")
permsPOD4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsPOD4")
permsPOD1all<-c(permsPOD1, permsPOD2, permsPOD3, permsPOD4)
penPOD<-calc.penalties(permsPOD1all, alpha=c(0.05, 0.1))

permsIOW1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsIOW1")
permsIOW2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsIOW2")
permsIOW3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsIOW3")
permsIOW4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsIOW4")
permsIOW1all<-c(permsIOW1, permsIOW2, permsIOW3, permsIOW4)
penIOW<-calc.penalties(permsIOW1all, alpha=c(0.05, 0.1))

# traits with covar fam
permsCPL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPL1")
permsCPL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPL2")
permsCPL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPL3")
permsCPL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPL4")
permsCPL1all<-c(permsCPL1, permsCPL2, permsCPL3, permsCPL4)
penCPL<-calc.penalties(permsCPL1all, alpha=c(0.05, 0.1))

permsCPD1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPD1")
permsCPD2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPD2")
permsCPD3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPD3")
permsCPD4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsCPD4")
permsCPD1all<-c(permsCPD1, permsCPD2, permsCPD3, permsCPD4)
penCPD<-calc.penalties(permsCPD1all, alpha=c(0.05, 0.1))

permsgutlength1<-readRDS("R_on_cluster/for_multiple_qtl_model/permsgutlength1")
permsgutlength2<-readRDS("R_on_cluster/for_multiple_qtl_model/permsgutlength2")
permsgutlength3<-readRDS("R_on_cluster/for_multiple_qtl_model/permsgutlength3")
permsgutlength4<-readRDS("R_on_cluster/for_multiple_qtl_model/permsgutlength4")
permsgutlength1all<-c(permsgutlength1, permsgutlength2, permsgutlength3, permsgutlength4)
pengutlength<-calc.penalties(permsgutlength1all, alpha=c(0.05, 0.1))

permstoothDensityL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityL1")
permstoothDensityL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityL2")
permstoothDensityL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityL3")
permstoothDensityL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityL4")
permstoothDensityL1all<-c(permstoothDensityL1, permstoothDensityL2, permstoothDensityL3, permstoothDensityL4)
penToothDensityL<-calc.penalties(permstoothDensityL1all, alpha=c(0.05, 0.1))

# traits with no covars
permstoothshapeL1<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothshapeL1")
permstoothshapeL2<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothshapeL2")
permstoothshapeL3<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothshapeL3")
permstoothshapeL4<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothshapeL4")
permstoothshapeL1all<-c(permstoothshapeL1, permstoothshapeL2, permstoothshapeL3, permstoothshapeL4)
penTootshapeL<-calc.penalties(permstoothshapeL1all, alpha=c(0.05, 0.1))

permstoothShapeU1<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothShapeU1")
permstoothShapeU2<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothShapeU2")
permstoothShapeU3<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothShapeU3")
permstoothShapeU4<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothShapeU4")
permstoothShapeU1all<-c(permstoothShapeU1, permstoothShapeU2, permstoothShapeU3, permstoothShapeU4)
penTootshapeU<-calc.penalties(permstoothShapeU1all, alpha=c(0.05, 0.1))

permstoothDensityU1<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityU1")
permstoothDensityU2<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityU2")
permstoothDensityU3<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityU3")
permstoothDensityU4<-readRDS("R_on_cluster/for_multiple_qtl_model/permstoothDensityU4")
permstoothDensityU1all<-c(permstoothDensityU1, permstoothDensityU2, permstoothDensityU3, permstoothDensityU4)
penToothDensityU<-calc.penalties(permstoothDensityU1all, alpha=c(0.05, 0.1))


#now run stepwise models 
HLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$HL, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penHL[1,])
# plots...
par(mfrow=c(1,2))
plotModel(HLmodel)
plotLodProfile(HLmodel, main="HL")
# path through model space
thetrace <- attr(HLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))


HWmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$HW, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penHW[1,])
par(mfrow=c(1,2))
plotModel(HWmodel)
plotLodProfile(HWmodel, main="HW")
thetrace <- attr(HWmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

LJLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$LJL, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penLJL[1,])
par(mfrow=c(1,2))
plotModel(LJLmodel)
plotLodProfile(LJLmodel, main="LJL")
thetrace <- attr(LJLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

LJWmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$LJW, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penLJW[1,])
par(mfrow=c(1,2))
plotModel(LJWmodel)
plotLodProfile(LJWmodel, main="LJW")
thetrace <- attr(LJWmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

SnLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$SnL, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penSnL[1,])
par(mfrow=c(1,2))
plotModel(SnLmodel)
plotLodProfile(SnLmodel, main="SnL")
thetrace <- attr(SnLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

SnWmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$SnW, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penSnW[1,])
par(mfrow=c(1,2))
plotModel(SnWmodel)
plotLodProfile(SnWmodel, main="SnW")
thetrace <- attr(SnWmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

EyLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$EyL, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penEyL[1,])
par(mfrow=c(1,2))
plotModel(EyLmodel)
plotLodProfile(EyLmodel, main="EyL")
thetrace <- attr(EyLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

EyDmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$EyD, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penEyD[1,])
par(mfrow=c(1,2))
plotModel(EyDmodel)
plotLodProfile(EyDmodel, main="EyD")
thetrace <- attr(EyDmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

POWmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$POW, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penPOW[1,])
par(mfrow=c(1,2))
plotModel(POWmodel)
plotLodProfile(POWmodel, main="POW")
thetrace <- attr(POWmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

InnerRowsmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$InnerRows-2, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penInnerRows[1,])
par(mfrow=c(1,2))
plotModel(InnerRowsmodel)
plotLodProfile(InnerRowsmodel, main="InnerRows")
thetrace <- attr(InnerRowsmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))
# only worked when adjusting counts (starting at 0 instead of 3)

Gapmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$Gap, max.qtl = 5, covar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                            incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penGap[1,])
par(mfrow=c(1,2))
plotModel(Gapmodel)
plotLodProfile(Gapmodel, main="Gap")
thetrace <- attr(Gapmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

BDmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$BD, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$sex), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penBD[1,])
par(mfrow=c(1,2))
plotModel(BDmodel)
plotLodProfile(BDmodel, main="BD")
thetrace <- attr(BDmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

ChDmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$ChD, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$sex), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penChD[1,])
par(mfrow=c(1,2))
plotModel(ChDmodel)
plotLodProfile(ChDmodel, main="ChD")
thetrace <- attr(ChDmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

PODmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$POD, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$sex), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penPOD[1,])
par(mfrow=c(1,2))
plotModel(PODmodel)
plotLodProfile(PODmodel, main="POD")
thetrace <- attr(PODmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

IOWmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$IOW, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$sex), method = "hk", model = "normal", 
                     incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penIOW[1,])
par(mfrow=c(1,2))
plotModel(IOWmodel)
plotLodProfile(IOWmodel, main="IOW")
thetrace <- attr(IOWmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

CPLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$CPL, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penCPL[1,])
par(mfrow=c(1,2))
plotModel(CPLmodel)
plotLodProfile(CPLmodel, main="CPL")
thetrace <- attr(CPLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

CPDmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$CPD, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penCPD[1,])
par(mfrow=c(1,2))
plotModel(CPDmodel)
plotLodProfile(CPDmodel, main="CPD")
thetrace <- attr(CPDmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

gutlengthmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$gutlength, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=pengutlength[1,])
par(mfrow=c(1,2))
plotModel(gutlengthmodel)
plotLodProfile(gutlengthmodel, main="gutlength")
thetrace <- attr(gutlengthmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

toothDensityLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$toothDensityL, max.qtl = 5, covar = as.data.frame(PNPxNOMdat$pheno$family), method = "hk", model = "normal", 
                      incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penToothDensityL[1,])
par(mfrow=c(1,2))
plotModel(toothDensityLmodel)
plotLodProfile(toothDensityLmodel, main="toothDensityL")
thetrace <- attr(toothDensityLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

toothshapeLmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$toothshapeL, max.qtl = 5, method = "hk", model = "normal", 
                                incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penTootshapeL[1,])
par(mfrow=c(1,2))
plotModel(toothshapeLmodel)
plotLodProfile(toothshapeLmodel, main="tootshapeL")
thetrace <- attr(toothshapeLmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

toothShapeUmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$toothShapeU, max.qtl = 5, method = "hk", model = "normal", 
                              incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penTootshapeU[1,])
par(mfrow=c(1,2))
plotModel(toothShapeUmodel)
plotLodProfile(toothShapeUmodel, main="tootShapeU")
thetrace <- attr(toothShapeUmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

toothDensityUmodel<-stepwiseqtl(cross = PNPxNOMdat, pheno.col = PNPxNOMdat$pheno$toothDensityU, max.qtl = 5, method = "hk", model = "normal", 
                              incl.markers = T, refine.locations = F, additive.only = F, scan.pairs = F, keeplodprofile = T, keeptrace = T, verbose = T, pen=penToothDensityU[1,])
par(mfrow=c(1,2))
plotModel(toothDensityUmodel)
plotLodProfile(toothDensityUmodel, main="toothDensityU")
thetrace <- attr(toothDensityUmodel, "trace")
par(mfrow=c(3,3))
for(i in seq(along=thetrace)) plotModel(thetrace[[i]], main=paste("pLOD =",round(attr(thetrace[[i]],"pLOD"), 2)))

# write output
{
sink("PNPxNOM_QTLmapping_results_MQM.txt")
cat("HLmodel","\\n")
print(HLmodel)
cat("\\n")
cat("HWmodel","\\n")
print(HWmodel)
cat("\\n")
cat("BDmodel","\\n")
print(BDmodel)
cat("\\n")
cat("LJLmodel","\\n")
print(LJLmodel)
cat("\\n")
cat("LJWmodel","\\n")
print(LJWmodel)
cat("\\n")
cat("SnLmodel","\\n")
print(SnLmodel)
cat("\\n")
cat("SnWmodel","\\n")
print(SnWmodel)
cat("\\n")
cat("EyLmodel","\\n")
print(EyLmodel)
cat("\\n")
cat("EyDmodel","\\n")
print(EyDmodel)
cat("\\n")
cat("ChDmodel","\\n")
print(ChDmodel)
cat("\\n")
cat("PODmodel","\\n")
print(PODmodel)
cat("\\n")
cat("IOWmodel","\\n")
print(IOWmodel)
cat("\\n")
cat("POWmodel","\\n")
print(POWmodel)
cat("\\n")
cat("CPLmodel","\\n")
print(CPLmodel)
cat("\\n")
cat("CPLmodel","\\n")
print(CPDmodel)
cat("\\n")
cat("gutlengthmodel","\\n")
print(gutlengthmodel)
cat("\\n")
cat("toothDensityLmodel","\\n")
print(toothDensityLmodel)
cat("\\n")
cat("toothDensityUmodel","\\n")
print(toothDensityUmodel)
cat("\\n")
cat("toothshapeLmodel","\\n")
print(toothshapeLmodel)
cat("\\n")
cat("toothShapeUmodel","\\n")
print(toothShapeUmodel)
cat("\\n")
cat("InnerRowsmodel","\\n")
print(InnerRowsmodel)
cat("\\n")
cat("Gapmodel","\\n")
print(Gapmodel)
cat("\\n")
sink()
}

