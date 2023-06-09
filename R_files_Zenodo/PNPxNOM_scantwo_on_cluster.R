#!/bin/R

# load library
library("qtl")

# read in data (previously copied it to scratch)
PNPxNOMdat <- read.cross("csv", "", "/cluster/scratch/annafell/GenoPheno.csv", map.function="kosambi")
PNPxNOMdat<-jittermap(PNPxNOMdat)

# prepare for mapping (calc.genoprob with stepsize 3 to make it faster)
PNPxNOMdat <- calc.genoprob(PNPxNOMdat, step=3, error.prob=0.05, map.function="kosambi")

# run two QTL models

# normal, covars sex and family: c("HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "POW", "InnerRows", "Gap")
PNPxNOMoutscantwo1 <- scantwo(PNPxNOMdat, pheno.col=c("HL", "HW", "LJL", "LJW", "SnL", "SnW", "EyL", "EyD", "POW", "InnerRows", "Gap"), model="normal", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
# save output
save(PNPxNOMoutscantwo1, file="PNPxNOMoutscantwo1")

# normal, covar sex: c("BD", "ChD", "POD", "IOW")
PNPxNOMoutscantwo2 <- scantwo(PNPxNOMdat, pheno.col=c("BD", "ChD", "POD", "IOW"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$sex)
# save output
save(PNPxNOMoutscantwo2, file="PNPxNOMoutscantwo2")

# normal, covar family: c("CPL", "CPD", "gutlength", "toothDensityL")
PNPxNOMoutscantwo3 <- scantwo(PNPxNOMdat, pheno.col=c("CPL", "CPD", "gutlength", "toothDensityL"), model="normal", method="em", addcovar = PNPxNOMdat$pheno$family)
# save output
save(PNPxNOMoutscantwo3, file="PNPxNOMoutscantwo3")

# normal, no covars: c("toothshapeL", "toothShapeU", "toothDensityU")
PNPxNOMoutscantwo4 <- scantwo(PNPxNOMdat, pheno.col=c("toothshapeL", "toothShapeU", "toothDensityU"), model="normal", method="em")
# save output
save(PNPxNOMoutscantwo4, file="PNPxNOMoutscantwo4")

# binary, covars sex and family: "tricuspid.tg."
PNPxNOMoutscantwo5 <- scantwo(PNPxNOMdat, pheno.col=c("tricuspid.tg."), model="binary", method="em", addcovar = cbind(PNPxNOMdat$pheno$sex, PNPxNOMdat$pheno$family))
# save output
save(PNPxNOMoutscantwo5, file="PNPxNOMoutscantwo5")


# subset to males for colour mapping
sex<-PNPxNOMdat$pheno$sex 
PNPxNOMdatM<-subset(PNPxNOMdat,ind=(sex==1))

# only males, binary, covar family: c("Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1", "Rdorsum2", "Yflanks", "Ygillcover")
PNPxNOMoutscantwo6 <- scantwo(PNPxNOMdatM, pheno.col=c("Rdorsalfin1", "Rdorsalfin2", "Rhead", "Rdorsum1", "Rdorsum2", "Yflanks", "Ygillcover"), model="binary", method="em", addcovar = PNPxNOMdatM$pheno$family)
# save output
save(PNPxNOMoutscantwo6, file="PNPxNOMoutscantwo6")

# only males, binary, no covar: c("Rflanks1.tg.", "Rcheek.tg.", "Rgillcover.tg.", "Rnose", "Yupperlip.tg.", "Ycheek.tg.", "Ypelvicfin.tg.", "Ynose")
PNPxNOMoutscantwo7 <- scantwo(PNPxNOMdatM, pheno.col=c("Rflanks1.tg.", "Rcheek.tg.", "Rgillcover.tg.", "Rnose", "Yupperlip.tg.", "Ycheek.tg.", "Ypelvicfin.tg.", "Ynose"), model="binary", method="em")
# save output
save(PNPxNOMoutscantwo7, file="PNPxNOMoutscantwo7")

# only males, normal, covar family: c("No_stripes")
PNPxNOMoutscantwo8 <- scantwo(PNPxNOMdatM, pheno.col=c("No_stripes"), model="normal", method="em", addcovar = PNPxNOMdatM$pheno$family)
# save output
save(PNPxNOMoutscantwo8, file="PNPxNOMoutscantwo8")
