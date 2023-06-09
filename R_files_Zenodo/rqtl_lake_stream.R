
# load R/qtl data
library(qtl)

#read file for QTL mapping
#QTL mapping done for each F1 family separately
#this script goes through steps for Boot family 1, but the same steps were repeated for each F1 family
boot_fam1 <- read.cross("csv","~/Desktop","Boot_fam1_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#boot_fam2 <- read.cross("csv","~/Desktop","Boot_fam2_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#misty_fam1 <- read.cross("csv","~/Desktop","Misty_fam1_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#misty_fam2 <- read.cross("csv","~/Desktop","Misty_fam2_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#pye_fam1 <- read.cross("csv","~/Desktop","Pye_fam1_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#pye_fam2 <- read.cross("csv","~/Desktop","Pye_fam2_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#pye_fam4 <- read.cross("csv","~/Desktop","Pye_fam4_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#rob_fam3 <- read.cross("csv","~/Desktop","Roberts_fam3_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")
#rob_fam5 <- read.cross("csv","~/Desktop","Roberts_fam5_qtl_final_usedTraits.csv" , genotypes = c("A", "H", "B"), na.strings = "NA")

#summary of .csv
summary(boot_fam1)
#summary(boot_fam2)
#summary(misty_fam1)
#summary(misty_fam2)
#summary(pye_fam1)
#summary(pye_fam2)
#summary(pye_fam4)
#summary(rob_fam3)
#summary(rob_fam5)

#show the linkage map
plot.map(boot_fam1)
#plot.map(boot_fam2)
#plot.map(misty_fam1)
#plot.map(misty_fam2)
#plot.map(pye_fam1)
#plot.map(pye_fam2)
#plot.map(pye_fam4)
#plot.map(rob_fam3)
#plot.map(rob_fam5)

#check genotyping data. check frequency
gt<-geno.table(boot_fam1)
summary(gt)

#get list of markers with genotypes that deviate significantly from expected 1:2:1 ratio
gt[gt$P.value<0.1,]

# use jittermap because markers in same position (see warning message)
boot_fam1 <- calc.genoprob(boot_fam1, step = 1, error.prob = 0.001)
boot_fam1 <- jittermap(boot_fam1)
sex <- pull.pheno(boot_fam1, pheno.col="sex")

#now do QTL mapping separately for each trait, standard length is given here as an example
#standard length (pheno.col=9)

# check for distribution of trait values
plotPheno(boot_fam1, pheno.col=9) #standard length

#do QTL mapping
out.std_length_mm <- scanone(boot_fam1, pheno.col = 9, addcovar = sex)
plot(out.std_length_mm)
title(main = "standard length")

summary(out.std_length_mm) #marker with highest LOD on each linkage group
summary(out.std_length_mm, threshold = 2.7) #marker with highest LOD on each linkage group that exceeds the threshold 

#permutation testing to identify significance threshold for each trait
operm <- scanone(boot_fam1, pheno.col = 9, n.perm=1000, verbose=FALSE)
plot(operm)
summary(operm, alpha=0.05)

