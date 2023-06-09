rm(list=ls())

library(qtl)
library(qtl2)

# read data 
f2 <- read.cross("csvs", "", "Geno.QTL.csv", "Pheno.QTL.csv")

#flip  order of some LG to be consistent with previous linkage map (based on crossing Crab-ecotypes)
f2 <- flip.order(f2, 6)#crab lg6
f2 <- flip.order(f2, 2)#crab lg 2
f2 <- flip.order(f2, 7)#crab lg 4
f2 <- flip.order(f2, 8)#crab lg 7

f2 <- convert2cross2(f2)


#### QTL analysis ##############
# insert pseudomarkers into the genetic map
map <- insert_pseudomarkers(f2$gmap, step=1)

#calculate the QTL genotype probabilities
pr <- calc_genoprob(f2, map, error_prob=0.002)

#### rename linkage group to be consisten with previous map
map2 <- map
map2$`12` <- map$`1`
map2$`1` <- map$`3`
map2$`3` <- map$`4`
map2$`17` <- map$`5`
map2$`4` <- map$`7`
map2$`7` <- map$`8`
map2$`14` <- map$`10`
map2$`10` <- map$`11`
map2$`5` <- map$`12`
map2$`15` <- map$`14`
map2$`11` <- map$`15`
map2$`8` <- map$`16`
map2$`16` <- map$`17`


#batch effect
covar <- f2$pheno[,"batch"]
## batch and sex as covariate
covar2 <- cbind(covar, f2$covar[,"sex"])
# correcting for length
covar3 <- cbind(covar2, f2$pheno[,"shellLength"])

#### find QTLs: THICKNESS, WEIGHT, LENGTH
scan.size <- scan1(pr,f2$pheno[,c("thickness","weight", "shellLength")], addcovar=covar2)
# find QTLs, LOD score, position
find_peaks(scan.size, map2)
# Confidence intervals of QTLs
lod_int(scan.size, map2, lodcolumn = 1, chr=6)
lod_int(scan.size, map2, lodcolumn = 2, chr=6)
lod_int(scan.size, map2, lodcolumn = 3, chr=6)

## QTLs for shell shape and aperture
scan.shape <- scan1(pr, f2$pheno[,c(6:10)], addcovar=covar3)
find_peaks(scan.shape, map2)

## QTL for sex
sex.scan <- scan1(pr, f2$pheno[,"sex.trait"], model="binary")
find_peaks(sex.scan, map2)
lod_int(sex.scan, lodcolumn = 1, map2, chr = 12)


