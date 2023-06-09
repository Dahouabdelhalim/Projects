
### Rqtl ###

rm(list=ls(all=TRUE))

library(qtl2)
library(qtl)
library(qtlcharts)




mapthis<-read.cross(format="csv",dir="C:/Users/OWNER/Desktop",
                    file="Genotype&Phenotype.SNP.csv",genotypes=c(AA = "1", BB = "3"),
                    na.string="NA",map.function="morgan",estimate.map = TRUE)

jittermap(mapthis)

mapthis <- calc.genoprob(mapthis)

scan1<-scanone(mapthis,pheno.col = 5,model = "normal", method = "hk")
plot(scan1)


lodint(scan1, chr = 3)


scan1_perm<- scanone(mapthis, pheno.col = 5, model = "normal", method = "hk", n.perm = 1000)
summary(scan1)



scan2<- scantwo(mapthis, pheno.col =5 , method = "mr")

plot(scan2, lower="fv1")


write.csv(z,"C:/Users/OWNER/Desktop/EPI_BILH.Bx.csv",quote=FALSE,row.names=FALSE)
a<-summary(scan2)
z<-a[a$lod.fv1>3,]
z





#######################################################################################################################

newmap<-(qtl2::convert2cross2(mapthis))
map <- insert_pseudomarkers(newmap$gmap)
pr <- calc_genoprob(newmap, map, error_prob=0.002)
pr <- calc_genoprob(newmap, map, error_prob=0.002, cores=4)
apr <- genoprob_to_alleleprob(pr)


kinship <- calc_kinship(pr)
grid <- calc_grid(newmap$gmap)
pr_grid <- probs_to_grid(pr, grid)
kinship_grid <- calc_kinship(pr_grid)
kinship_loco <- calc_kinship(pr, "loco")


out <- scan1(pr, newmap$pheno, cores = 4)

par(mar=c(5.1, 4.1, 1.1, 1.1))
ymx <- maxlod(out) # overall maximum LOD score
plot(out, map, lodcolumn=5, col="gray8", lwd)



lod_int(out, map, lodcolumn = 1, chr = 1)

bayes_int(out, map, lodcolumn=1, chr=1, prob=0.95)

find_peaks(out, map, threshold=4, peakdrop=1.8, drop=1.5)




