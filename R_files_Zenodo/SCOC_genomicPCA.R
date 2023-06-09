#these are packages you'll need to install using "install.packages"
library(vcfR); library(poppr); library(maps); library(mapdata); library(plotrix); 
library(adegenet); library(calibrate); library(mapproj); library(ape)

setwd("<directory>")

#Real data from VCF
my.loci <- read.vcfR("SCOC_50missing.vcf")

snps <- vcfR2genind(my.loci) 
#need to pick number of axes to retain

gl<-vcfR2genlight(my.loci)

pca1 <- glPca(gl)
scatter(pca1)

#loadingplot(pca1)

myCol <- colorplot(pca1$scores, pca1$scores, transp = TRUE, cex = 4)

abline(h=0, v=0, col = "grey")
text(pca1$scores[, 1], pca1$scores[, 2], gl$ind.names, cex=0.5, pos = 1)

#find the optimal number of clusters, then you MUST select the value. It prompts for a value

#Enter your optimal number of clusters in "n.clust= x" below:
myclust <- find.clusters(snps, n.pca = 100, n.clust = 4)
mydapc <- dapc(snps, pop = myclust$grp, n.pca = 100, n.da = 2)
piecolors <- topo.colors(4, alpha = 1) #set colors to match number of clusters

# Grab posterior as data matrix
mypost <- mydapc$posterior
#mypost <- pca1$scores
