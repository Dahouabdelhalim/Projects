## PACKAGES
## KINDISPERSE
install.packages("devtools")
devtools::install_github("moshejasper/kindisperse")
library(kindisperse)
library(ggplot2)

## dbRDA
require(vegan)
require(ecodist)

## sPCA
require(adespatial)
require(calibrate)
require(vcfR)
require(ade4)
require(sp)
install.packages("adegenet",dependencies = TRUE)
library(adegenet)


########################################
## KINDISPERSE

fs = read.delim("full-sibs.txt",header = F)
hs = read.delim("half-sibs.txt",header = F)
cs = read.delim("cousins.txt",header = F)

fullsibs=as.vector(fs)
halfsibs=as.vector(hs)
cousins=as.vector(cs)

## sigma ovi (ovipositional) 
axpermute(fullsibs$V1, nreps = 10000, composite = 2)

## sigma PO (parent-offspring)
axpermute_standard(avect = cousins$V1, bvect = fullsibs$V1, acat = "1C", bcat = "FS", amix = TRUE, amixcat = "H1C", 
                   bcomp = T, bcompvect = halfsibs$V1, bcompcat = "HS",nreps = 10000)



########################################
## dbRDA

## 
#testing road effects across western road 

#road_west is list of all individuals in a dyad that is separated by the western road, and which side of the road each individual is on  
road_west  = read.delim("road_west.txt",header = TRUE, row.names = 1)

#gen_west is nxn matrix of pairwise Rousset's a scores between all dyads separated by the western road  
gen_west = read.delim("gen_west.txt",header = TRUE, row.names = 1)

#geo_west is nxn matrix of pairwise geographical distances between all dyads separated by the western road  
geo_west = read.delim("geo_west.txt",header = TRUE, row.names = 1)
geo_west_dist = as.dist(geo_west,diag = TRUE,upper = TRUE)
geo_west_pc = pcnm(geo_west_dist)

#this determines which, if any, PCs of geographical distance should be included in the model
west_sigdist = capscale(gen_west ~ geo_west_pc$vectors)
anova(west_sigdist, by="margin", permutations=999)

#as above but for each PC 
west_sigdist = capscale(gen_west ~ geo_west_pc$vectors[,1] + geo_west_pc$vectors[,2] + geo_west_pc$vectors[,3] +
                        geo_west_pc$vectors[,4] + geo_west_pc$vectors[,5] + geo_west_pc$vectors[,6] + 
                        geo_west_pc$vectors[,7] + geo_west_pc$vectors[,8] + geo_west_pc$vectors[,9])
anova(west_sigdist, by="margin", permutations=9999)

#the model, including no PCs ofdistance as none were significant
dbrda_west = capscale(gen_west ~ road_west$road)
anova.cca(dbrda_west, by="margin", permutations=9999)

## 
#testing road effects across southern road 

road_south  = read.delim("road_south.txt",header = TRUE, row.names = 1)
geo_south = read.delim("geo_south.txt",header = TRUE, row.names = 1)
gen_south = read.delim("gen_south.txt",header = TRUE, row.names = 1)

geo_south_dist = as.dist(geo_south,diag = TRUE,upper = TRUE)
geo_south_pc = pcnm(geo_south_dist)

south_sigdist = capscale(gen_south ~ geo_south_pc$vectors)
anova(south_sigdist, by="margin", permutations=999)

south_sigdist = capscale(gen_south ~ geo_south_pc$vectors[,1] + geo_south_pc$vectors[,2] + geo_south_pc$vectors[,3] +
                        geo_south_pc$vectors[,4] + geo_south_pc$vectors[,5] + geo_south_pc$vectors[,6] + 
                        geo_south_pc$vectors[,7] + geo_south_pc$vectors[,8] + geo_south_pc$vectors[,9])
anova(south_sigdist, by="margin", permutations=9999)


dbrda_south = capscale(gen_south ~ road_south$road)
anova.cca(dbrda_south, by="margin", permutations=9999)


########################################
## sPCA

vcf<-read.vcfR("safa.vcf")
safa_gen<-vcfR2genind(vcf)
safa_gen_pop=read.delim("safa_gen_pops.txt",row.names = 1)
safa_gen$pop=as.factor(safa_gen_pop$pop)
safa_xy=read.delim("safa_gen_xy.txt")
safa_gen@other$xy = safa_xy

mySpca <- spca(safa_gen, ask=FALSE, type=5, scannf=FALSE)

