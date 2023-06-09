library(geomorph)
library(Morpho)
library(shapes)
library(ggplot2)
library(GGally)
library(lmerTest)
library(lme4)
library(emmeans)
library(multcomp)
library(lattice)
library(Momocs)

#setwd

# Load landmark data
dat.treehop<-scan("TreehopperTPS.Dorsal.GMdata.TPS", what="character") # 7880 items

# Load classifiers
classifiersAll.dorsal<-read.csv("Treehopper.GM.Classifiers.csv", header=T)
dim(classifiersAll.dorsal) # 197 rows and 6 columns

# Separate out classifiers
stage.d<-classifiersAll.dorsal$Stage
stage.d<-as.factor(stage.d)
genus.d<-classifiersAll.dorsal$Genus
genus.d<-as.factor(genus.d)
genstage.d<-classifiersAll.dorsal$GenusStage
genstage.d<-as.factor(genstage.d)


## Organize TPS landmark data and make it to scale
length(dat.treehop) # length is 7880, 197 specimens
dat.treehop[1:100]
# identify which lines of the vector belong to scale
c(1:197)*40
#Using regular sequence, we can get image names and coordinates.
Scale.All<-dat.treehop[c(1:197)*40]
Scale.All
Scale.Alld<-sub("SCALE=", "", Scale.All)
Scale.Alld
mode(Scale.Alld)<-"numeric"
is.numeric(Scale.Alld)

# generate an empty array to store the 197 configurations (specimens).
n<-197
All.Scaled<-array(NA, dim=c(18,2,197)) 
for (i in 1:n){All.Scaled[,,i]<-matrix(dat.treehop[(i*40-38):(i*40-3)],18,2,byrow=T)}

All.Scaled[,,12]
mode(All.Scaled)<-"numeric"

# want to use this loop to set the scale, scale is 1mm=0.001 pixels or so
# set the scale
A.InsAll.d<-All.Scaled
for (i in 1:41)
{A.InsAll.d[,,i]<-(All.Scaled[,,i]*(Scale.Alld[i]))}
A.InsAll.d<-A.InsAll.d[1:18,,]
A.InsAll.d
dim(A.InsAll.d)

dev.off

## Procrustes Superimposition
#define.sliders(A.InsAll.d, nsliders=16)
slidersDorsal<-read.csv("curveslideDORSAL.csv")
hop.super.Dorsal<-gpagen(A.InsAll.d, curves = slidersDorsal)
hop.super.Dorsal<-rotate.coords(hop.super.Dorsal, type = "rotateC")
hop.super.Dorsal<-rotate.coords(hop.super.Dorsal, type = "rotateC")
plot(hop.super.Dorsal)


## Preliminary PCA figure of the Procrustes superimposed data
PCA.dorsal<-gm.prcomp(hop.super.Dorsal$coords)


plot(PCA.dorsal, col=(3:8)[genus.d], pch=(3:5)[stage.d], xlab="PC1", ylab="PC2")
legend(locator(1), pch=22, pt.bg=3:8,  legend=c("Aetalion","Entylia", "Ennya",
                                                "Membracis", "Metheisa", "Polyglypta"))
legend(locator(1), pch=(3:5), pt.bg=3:8, cex=1, legend=c("Ins4", "Ins5", "Ad"))

plot(PCA.dorsal$x[,2], PCA.dorsal$x[,3], col=(3:8)[genus.d], pch=(3:5)[stage.d], xlab="PC2", ylab="PC3")
legend(locator(1), pch=22, pt.bg=3:8,  legend=c("Aetalion","Entylia", "Ennya",
                                                "Membracis", "Metheisa", "Polyglypta"))
legend(locator(1), pch=(3:5), pt.bg=3:8, cex=1, legend=c("Ins4", "Ins5", "Ad"))


## Extremes of PC1 and PC2
dorsalshapes<-PCA.dorsal$shapes
dorsalshape1<-dorsalshapes$shapes.comp1
dorsalshape2<-dorsalshapes$shapes.comp2
plot(dorsalshape1$min)
plot(dorsalshape1$max)
plot(dorsalshape2$min)
plot(dorsalshape2$max)


## Testing for the influence of size on shape, we expect that there is one since these are ontogenetic data
allometry.test.d.lm <- procD.lm(hop.super.Dorsal$coords ~ log(hop.super.Dorsal$Csize) * genus.d, iter=999)
summary(allometry.test.d.lm)

# Residual shape score -- the variation that remains when the effect of size is removed
RSC.dorsal <-plotAllometry(allometry.test.d.lm, hop.super.Dorsal$Csize, method="CAC", pch=19, col=(3:8)[genus.d])$RSC
cor.test(test.CAC.dorsal, PCA.dorsal$x[,1])
RSC.dorsal2d<-arrayspecs(RSC.dorsal, 18, 2)

PC1dorsal<-PCA.dorsal$x[,1]
PC2dorsal<-PCA.dorsal$x[,2]

df.Dorsal<-data.frame(PC1dorsal, PC2dorsal)
colors.for.dplot<-c("darkslateblue", "lightskyblue3", "plum4", "goldenrod2", "sienna3", "maroon")
dplot<-ggplot(df.Dorsal, aes(PC1dorsal,PC2dorsal, colour=(genus.d)))  + geom_point(aes(shape=stage.d), size=3)
dplot+ theme_classic() + scale_color_manual(values=colors.for.dplot)


## PCA plot using RSC1 and RSC2
plot(RSC.dorsal[,1], RSC.dorsal[,2], col=(3:8)[genus.d], pch=(3:5)[stage.d], xlab="RSC1", ylab="RSC2")
legend(locator(1), pch=22, pt.bg=3:8,  legend=c("Aetalion","Entylia", "Ennya",
                                                "Membracis", "Metheisa", "Polyglypta"))
legend(locator(1), pch=(3:5), pt.bg=3:8, cex=1, legend=c("Ins4", "Ins5", "Ad"))


## To test for significance in the mean shape of the different stages, subset the data by genera
Genus_list.d<- unique(genus.d)
test_genusd_mean_shape_LMs <- list()
test_genusd_mean_shape_LMs_pairwise <- list()

i=1
for (i in 1:length(Genus_list.d)){
  temp_genus.d <- Genus_list.d[i]
  temp_grep2.d <- grep(temp_genus.d, genus.d)
  x<-hop.super.Dorsal$coords[,,temp_grep2.d]
  #x1<-PC2[temp_grep2]
  y<-stage.d[temp_grep2.d]
  temp_lm2.d<-procD.lm(x ~ y, iter = 9999)
  #temp_lm2 <- procD.lm(hop.super.All$coords[,,temp_grep2] ~ stage[temp_grep2])
  
  test_genusd_mean_shape_LMs[[temp_genus.d]] <- temp_lm2.d
  test_genusd_mean_shape_LMs_pairwise[[temp_genus.d]] <- summary(pairwise(temp_lm2.d, groups = factor(stage.d[temp_grep2.d])))
  
}

lapply(test_genusd_mean_shape_LMs, "summary") #Summary of differences in mean shape among developmental stages of each genus

test_genusd_mean_shape_LMs_pairwise

# Aetalion, 4 and 5 are more different 
p.adjust(c(0.0030, 0.1508, 0.1741)) # [1] 0.0090 0.3016 0.3016
# Entylia, 4 and 5 are more different
p.adjust(c(0.0001, 0.0001, 0.0034)) # [1] 0.0003 0.0003 0.0034
# Ennya, 4 and 5 are more different
p.adjust(c(0.0016, 0.0001, 0.0012)) # [1] 0.0024 0.0003 0.0024
# Membracis, 5 and adult are more different
p.adjust(c(0.0941, 0.0001, 0.0001)) # [1] 0.0941 0.0003 0.0003
# Metheisa, 5 and adult are more different
p.adjust(c(0.1438, 0.0001, 0.0003)) #[1] 0.1438 0.0003 0.0006
# Polyglypta, 4 and 5 are more different
p.adjust(c(0.0024, 0.001, 0.0218)) # [1] 0.0048 0.0030 0.0218


## To test for significance in the mean shape of the different stages, subset the data by genera

Genus_list.d<- unique(genus.d)
test_genusd_mean_shape_LMsRSC <- list()
test_genusd_mean_shape_LMs_pairwiseRSC <- list()

i=1
for (i in 1:length(Genus_list.d)){
  temp_genus.d <- Genus_list.d[i]
  temp_grep2.d <- grep(temp_genus.d, genus.d)
  #x<-hop.super.Dorsal$coords[,,temp_grep2.d]
  x1<-RSC.dorsal2d[,,temp_grep2.d]
  y<-stage.d[temp_grep2.d]
  temp_lm2.d<-procD.lm(x1 ~ y, iter = 9999)
  #temp_lm2 <- procD.lm(hop.super.All$coords[,,temp_grep2] ~ stage[temp_grep2])
  
  test_genusd_mean_shape_LMsRSC[[temp_genus.d]] <- temp_lm2.d
  test_genusd_mean_shape_LMs_pairwiseRSC[[temp_genus.d]] <- summary(pairwise(temp_lm2.d, groups = factor(stage.d[temp_grep2.d])))
  
}

lapply(test_genusd_mean_shape_LMsRSC, "summary") #Summary of differences in mean shape among developmental stages of each genus

test_genusd_mean_shape_LMs_pairwiseRSC


### Analysis included in the paper ### 
## Remove Aetalion since it's not of interest in this analysis and could be influencing the Procrustes Fit
genus.noAet<-subset(genus.d, genus.d !="Aet")
stage.noAet<-subset(stage.d, genus.d != "Aet")
genstage.noAet2<-subset(genstage.d, genstage.d != "Aet4")
genstage.noAet1<-subset(genstage.noAet2, genstage.noAet2 != "Aet5")
genstage.noAet<-subset(genstage.noAet1, genstage.noAet1 != "AetAd")

A.InsAll.d2d<-two.d.array(A.InsAll.d)
DorsalnoAet<-subset(A.InsAll.d2d, genus.d !="Aet")
DorsalnoAet.matrix<-arrayspecs(DorsalnoAet, 18, 2)

slidersDorsal<-read.csv("curveslideDORSAL.csv")
noAet.hop.super.Dorsal<-gpagen(DorsalnoAet.matrix, curves = slidersDorsal)
noAet.hop.super.Dorsal<-rotate.coords(noAet.hop.super.Dorsal, type = "rotateC")
noAet.hop.super.Dorsal<-rotate.coords(noAet.hop.super.Dorsal, type = "rotateC")
plot(noAet.hop.super.Dorsal)

PCA.dorsal.noAet<-gm.prcomp(noAet.hop.super.Dorsal$coords)

## Significance Testing
## Pairwise comparisons of the shapes of membracid species across developmental stages
Genus_list.d.noAet<- unique(genus.noAet)
test_genusdNoAet_mean_shape_LMs <- list()
test_genusdNoAet_mean_shape_LMs_pairwise <- list()

i=1
for (i in 1:length(Genus_list.d.noAet)){
  temp_genus.dnoAet <- Genus_list.d.noAet[i]
  temp_grep2.dnoAet <- grep(temp_genus.dnoAet, genus.noAet)
  x<-noAet.hop.super.Dorsal$coords[,,temp_grep2.dnoAet]
  #x1<-PC2[temp_grep2]
  y<-stage.d[temp_grep2.dnoAet]
  temp_lm2.dnoAet<-procD.lm(x ~ y, iter = 9999)
  #temp_lm2 <- procD.lm(hop.super.All$coords[,,temp_grep2] ~ stage[temp_grep2])
  
  test_genusdNoAet_mean_shape_LMs[[temp_genus.dnoAet]] <- temp_lm2.dnoAet
  test_genusdNoAet_mean_shape_LMs_pairwise[[temp_genus.dnoAet]] <- summary(pairwise(temp_lm2.dnoAet, groups = factor(stage.d[temp_grep2.dnoAet])))
  
}

lapply(test_genusdNoAet_mean_shape_LMs, "summary") #Summary of differences in mean shape among developmental stages of each genus

test_genusdNoAet_mean_shape_LMs_pairwise

## Significance of Procrustes Difference within Species, between Stages
# Entylia, 4 and 5 are more different
p.adjust(c(0.0119, 0.0001, 0.0181)) # [1] 0.0238 0.0003 0.0238
# Ennya, 4 and 5 are more different
p.adjust(c(0.0041, 0.0002, 0.5236)) # [1] 0.0082 0.0006 0.5236
# Membracis, 5 and adult are more different
p.adjust(c(0.0001, 0.9270, 0.0002)) #[1] 0.0003 0.9270 0.0004
# Metheisa, 5 and adult are more different
p.adjust(c(0.4606, 0.0002, 0.0023)) #[1] 0.4606 0.0006 0.0046
# Polyglypta, 4 and 5 are more different
p.adjust(c(0.0021, 0.0002, 0.5845)) # [1] 0.0042 0.0006 0.5845


dorsalshapesnoAet<-PCA.dorsal.noAet$shapes
dorsalshape1noAet<-dorsalshapesnoAet$shapes.comp1
dorsalshape2noAet<-dorsalshapesnoAet$shapes.comp2
dorsalshape3noAet<-dorsalshapesnoAet$shapes.comp3
plot(dorsalshape1noAet$min)
plot(dorsalshape1noAet$max)
plot(dorsalshape2noAet$min)
plot(dorsalshape2noAet$max)
plot(dorsalshape3noAet$min)
plot(dorsalshape3noAet$max)

## PCs 1, 2, and 3 account for 89% of the total variance with PC4 accounting for only 4.6% of the variance, not a significant amount with a 
# cutoff of 0.05
PC1dorsalnoAet<-PCA.dorsal.noAet$x[,1]
PC2dorsalnoAet<-PCA.dorsal.noAet$x[,2]
PC3dorsalnoAet<-PCA.dorsal.noAet$x[,3]

df.DorsalnoAet<-data.frame(PC1dorsalnoAet, PC2dorsalnoAet, PC3dorsalnoAet)
colors.for.dplot.noAet<-c("lightskyblue3", "plum4", "goldenrod2", "sienna3", "maroon")

Species<-c("Entylia carinata","Ennya chrysura", "Membracis mexicana","Metheisa lucillodes", "Polyglypta costata")

## Final Figure of PC1 and PC2 with Species and Stages
#tiff("PC1 and PC2 of Dorsal GM2.tiff", units="in", width=8, height=4, res=250)
dplot.noAet<-ggplot(df.DorsalnoAet, aes(PC1dorsalnoAet,PC2dorsalnoAet, colour=(genus.noAet))) + geom_point(aes(shape=stage.noAet), size=3)
dplot.noAet + theme_classic() + scale_color_manual(values=colors.for.dplot.noAet, labels = Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) + labs(x="PC1 (59.4%)", y="PC2 (19.5%)", colour="Species", shape="Stage")

tiff("New.PC1 and PC2 of Dorsal GM2.tiff", units="in", width=8, height=4, res=300)
New.dplot.noAet<-ggplot(df.DorsalnoAet, aes(PC1dorsalnoAet,PC2dorsalnoAet, colour=(genus.noAet))) + geom_point(aes(shape=genus.noAet, color= genus.noAet, size=stage.noAet)) 
New.dplot.noAet + theme_classic() + 
  scale_size_manual(values=c(1,2.5,5)) +
  scale_color_manual(name  ="Species",
                     breaks=c("Eca", "Ech", "Mem", "Met", "Pol"),
                     values=colors.for.dplot.noAet,
                     labels=Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
  scale_shape_manual(name  ="Species",
                     breaks=c("Eca", "Ech", "Mem", "Met", "Pol"),
                     values=c(8, 15, 16, 17, 18),
                     labels=Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
                                                         #, labels = Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) + 
labs(x="PC1 (59.4%)", y="PC2 (19.5%)") 
dev.off()
  # scale_size_manual(values=c(1.5,2,3.5)) +
+ 
 # scale_colour_discrete(name  ="Species",
  #                      breaks=c("Eca", "Ech", "Mem", "Met", "Pol"),
   #                     values=colors.for.dplot.noAet,
    #                    labels=Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
#  scale_shape_discrete(name  ="Species",
 #                      breaks=c("Eca", "Ech", "Mem", "Met", "Pol"),
  #                     values=colors.for.dplot.noAet,
   #                    labels=Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) #, colour="Species", shape="Species")
?scale_colour_discrete

tiff("PC2 and PC3 of Dorsal GM.tiff", units="in", width=7, height=5, res=200)
dplot.noAet2<-ggplot(df.DorsalnoAet, aes(PC2dorsalnoAet,PC3dorsalnoAet, colour=(genus.noAet)))  + geom_point(aes(shape=stage.noAet), size=3)
dplot.noAet2+ theme_classic() + scale_color_manual(values=colors.for.dplot.noAet, labels = Species, guide =guide_legend(label.theme = element_text(angle = 0, face = "italic"))) + labs(y="PC3 (10.1%)", x="PC2 (19.5%)", colour="Species", shape="Stage")
dev.off()


### Use the RSC to remove shape that corresponds with size
allometry.test.d.lmnoAet <- procD.lm(noAet.hop.super.Dorsal$coords ~ log(noAet.hop.super.Dorsal$Csize) * genus.noAet, iter=999)
summary(allometry.test.d.lmnoAet)
noAetRSC.dorsal <-plotAllometry(allometry.test.d.lmnoAet, noAet.hop.super.Dorsal$Csize, method="CAC", pch=19, col=(3:8)[genus.noAet])$RSC
cor.test(test.CAC.dorsal, PCA.dorsal$x[,1])
noAetRSC.dorsal2d<-arrayspecs(noAetRSC.dorsal, 18, 2)


## Pairwise comparisons with RSC and group with no Aetalion
Genus_list.d.noAetRSC<- unique(genus.noAet)
test_genusdNoAetRSC_mean_shape_LMs <- list()
test_genusdNoAetRSC_mean_shape_LMs_pairwise <- list()

i=1
for (i in 1:length(Genus_list.d.noAetRSC)){
  temp_genus.dnoAetRSC <- Genus_list.d.noAetRSC[i]
  temp_grep2.dnoAetRSC <- grep(temp_genus.dnoAetRSC, genus.noAet)
  x<-noAetRSC.dorsal2d[,,temp_grep2.dnoAetRSC]
  #x1<-PC2[temp_grep2]
  y<-stage.noAet[temp_grep2.dnoAetRSC]
  temp_lm2.dnoAetRSC<-procD.lm(x ~ y, iter = 9999)
  #temp_lm2 <- procD.lm(hop.super.All$coords[,,temp_grep2] ~ stage[temp_grep2])
  
  test_genusdNoAet_mean_shape_LMs[[temp_genus.dnoAetRSC]] <- temp_lm2.dnoAetRSC
  test_genusdNoAet_mean_shape_LMs_pairwise[[temp_genus.dnoAetRSC]] <- summary(pairwise(temp_lm2.dnoAetRSC, groups = factor(stage.noAet[temp_grep2.dnoAetRSC])))
  
}

lapply(test_genusdNoAetRSC_mean_shape_LMs, "summary") #Summary of differences in mean shape among developmental stages of each genus

test_genusdNoAetRSC_mean_shape_LMs_pairwise



### Mean shapes of all species by stage 
## Putting the gpagen aligned data into a 2D matrix
nAll.d<-length(genstage.noAet)
A.InsAll.d2d<-matrix(NA, nAll.d, 36)
for (i in 1:nAll.d){A.InsAll.d2d[i,]<-noAet.hop.super.Dorsal$coords[,,i]}


## Entylia
EcaIns4.sub<-subset(A.InsAll.d2d, genstage.noAet=="Eca4")
EcaIns4<-subset(genstage.noAet, genstage.noAet=="Eca4")
EcaIns4.mean<-rowsum(EcaIns4.sub, EcaIns4)/length(EcaIns4)
EcaIns4.yall<-matrix(EcaIns4.mean, 18, 2)
plot(EcaIns4.yall, col=("lightskyblue3"), pch=16, cex=4)

EcaIns5.sub<-subset(A.InsAll.d2d, genstage.noAet=="Eca5")
EcaIns5<-subset(genstage.noAet, genstage.noAet=="Eca5")
EcaIns5.mean<-rowsum(EcaIns5.sub, EcaIns5)/length(EcaIns5)
EcaIns5.yall<-matrix(EcaIns5.mean, 18, 2)
plot(EcaIns5.yall, col=("lightskyblue3"), pch=16, cex=4)

EcaAd.sub<-subset(A.InsAll.d2d, genstage.noAet=="EcaAd")
EcaAd<-subset(genstage.noAet, genstage.noAet=="EcaAd")
EcaAd.mean<-rowsum(EcaAd.sub, EcaAd)/length(EcaAd)
EcaAd.yall<-matrix(EcaAd.mean, 18, 2)
plot(EcaAd.yall, col=("lightskyblue3"), pch=16, cex=4)


## Ennya
EchIns4.sub<-subset(A.InsAll.d2d, genstage.noAet=="Ech4")
EchIns4<-subset(genstage.noAet, genstage.noAet=="Ech4")
EchIns4.mean<-rowsum(EchIns4.sub, EchIns4)/length(EchIns4)
EchIns4.yall<-matrix(EchIns4.mean, 18, 2)
plot(EchIns4.yall, col=("plum4"), pch=16, cex=4)

EchIns5.sub<-subset(A.InsAll.d2d, genstage.noAet=="Ech5")
EchIns5<-subset(genstage.noAet, genstage.noAet=="Ech5")
EchIns5.mean<-rowsum(EchIns5.sub, EchIns5)/length(EchIns5)
EchIns5.yall<-matrix(EchIns5.mean, 18, 2)
plot(EchIns5.yall, col=("plum4"), pch=16, cex=4)

EchAd.sub<-subset(A.InsAll.d2d, genstage.noAet=="EchAd")
EchAd<-subset(genstage.noAet, genstage.noAet=="EchAd")
EchAd.mean<-rowsum(EchAd.sub, EchAd)/length(EchAd)
EchAd.yall<-matrix(EchAd.mean, 18, 2)
plot(EchAd.yall, col=("plum4"), pch=16, cex=4)

## Membracis
MemIns4.sub<-subset(A.InsAll.d2d, genstage.noAet=="Mem4")
MemIns4<-subset(genstage.noAet, genstage.noAet=="Mem4")
MemIns4.mean<-rowsum(MemIns4.sub, MemIns4)/length(MemIns4)
MemIns4.yall<-matrix(MemIns4.mean, 18, 2)
plot(MemIns4.yall, col=("goldenrod2"), pch=16, cex=4)

MemIns5.sub<-subset(A.InsAll.d2d, genstage.noAet=="Mem5")
MemIns5<-subset(genstage.noAet, genstage.noAet=="Mem5")
MemIns5.mean<-rowsum(MemIns5.sub, MemIns5)/length(MemIns5)
MemIns5.yall<-matrix(MemIns5.mean, 18, 2)
plot(MemIns5.yall, col=("goldenrod2"), pch=16, cex=4)

MemAd.sub<-subset(A.InsAll.d2d, genstage.noAet=="MemAd")
MemAd<-subset(genstage.noAet, genstage.noAet=="MemAd")
MemAd.mean<-rowsum(MemAd.sub, MemAd)/length(MemAd)
MemAd.yall<-matrix(MemAd.mean, 18, 2)
plot(MemAd.yall, col=("goldenrod2"), pch=16, cex=4)


## Metheisa
MetIns4.sub<-subset(A.InsAll.d2d, genstage.noAet=="Met4")
MetIns4<-subset(genstage.noAet, genstage.noAet=="Met4")
MetIns4.mean<-rowsum(MetIns4.sub, MetIns4)/length(MetIns4)
MetIns4.yall<-matrix(MetIns4.mean, 18, 2)
plot(MetIns4.yall, col=("sienna3"), pch=16, cex=4)

MetIns5.sub<-subset(A.InsAll.d2d, genstage.noAet=="Met5")
MetIns5<-subset(genstage.noAet, genstage.noAet=="Met5")
MetIns5.mean<-rowsum(MetIns5.sub, MetIns5)/length(MetIns5)
MetIns5.yall<-matrix(MetIns5.mean, 18, 2)
plot(MetIns5.yall, col=("sienna3"), pch=16, cex=4)

MetAd.sub<-subset(A.InsAll.d2d, genstage.noAet=="MetAd")
MetAd<-subset(genstage.noAet, genstage.noAet=="MetAd")
MetAd.mean<-rowsum(MetAd.sub, MetAd)/length(MetAd)
MetAd.yall<-matrix(MetAd.mean, 18, 2)
plot(MetAd.yall, col=("sienna3"), pch=16, cex=4)


## Polyglypta
PolIns4.sub<-subset(A.InsAll.d2d, genstage.noAet=="Pol4")
PolIns4<-subset(genstage.noAet, genstage.noAet=="Pol4")
PolIns4.mean<-rowsum(PolIns4.sub, PolIns4)/length(PolIns4)
PolIns4.yall<-matrix(PolIns4.mean, 18, 2)
plot(PolIns4.yall, col=("maroon"), pch=16, cex=4)

PolIns5.sub<-subset(A.InsAll.d2d, genstage.noAet=="Pol5")
PolIns5<-subset(genstage.noAet, genstage.noAet=="Pol5")
PolIns5.mean<-rowsum(PolIns5.sub, PolIns5)/length(PolIns5)
PolIns5.yall<-matrix(PolIns5.mean, 18, 2)
plot(PolIns5.yall, col=("maroon"), pch=16, cex=4)

PolAd.sub<-subset(A.InsAll.d2d, genstage.noAet=="PolAd")
PolAd<-subset(genstage.noAet, genstage.noAet=="PolAd")
PolAd.mean<-rowsum(PolAd.sub, PolAd)/length(PolAd)
PolAd.yall<-matrix(PolAd.mean, 18, 2)
plot(PolAd.yall, col=("maroon"), pch=16, cex=4)

Ins4.18<-c("Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4", "Ins4")
Ins5.18<-c("Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5", "Ins5")
Ad.18<-c("Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad", "Ad")

Eca.18<-c("Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca", "Eca")
Ech.18<-c("Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech", "Ech")
Mem.18<-c("Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem", "Mem")
Met.18<-c("Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met", "Met")
Pol.18<-c("Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol", "Pol")

# Pol
Pol4dat<-cbind(PolIns4.yall, Pol.18, Ins4.18)
#Pol4dat<-as.data.frame(Pol4dat)
Pol5dat<-cbind(PolIns5.yall, Pol.18, Ins5.18)
#Pol5dat<-as.data.frame(Pol4dat)
PolAddat<-cbind(PolAd.yall, Pol.18, Ad.18)
#PolAddat<-as.data.frame(PolAddat)

# Met
Met4dat<-cbind(MetIns4.yall, Met.18, Ins4.18)
#Met4dat<-as.data.frame(Met4dat)
Met5dat<-cbind(MetIns5.yall, Met.18, Ins5.18)
#Met5dat<-as.data.frame(Met5dat)
MetAddat<-cbind(MetAd.yall, Met.18, Ad.18)
#MetAddat<-as.data.frame(MetAddat)

# Mem
Mem4dat<-cbind(MemIns4.yall, Mem.18, Ins4.18)
#Mem4dat<-as.data.frame(Mem4dat)
Mem5dat<-cbind(MemIns5.yall, Mem.18, Ins5.18)
#Mem5dat<-as.data.frame(Mem5dat)
MemAddat<-cbind(MemAd.yall, Mem.18, Ad.18)
#MemAddat<-as.data.frame(MemAddat)

# Ech
Ech4dat<-cbind(EchIns4.yall, Ech.18, Ins4.18)
#Ech4dat<-as.data.frame(Ech4dat)
Ech5dat<-cbind(EchIns5.yall, Ech.18, Ins5.18)
#Ech5dat<-as.data.frame(Ech5dat)
EchAddat<-cbind(EchAd.yall, Ech.18, Ad.18)
#EchAddat<-as.data.frame(EchAddat)

# Eca
Eca4dat<-cbind(EcaIns4.yall, Eca.18, Ins4.18)
#Eca4dat<-as.data.frame(Eca4dat)
Eca5dat<-cbind(EcaIns5.yall, Eca.18, Ins5.18)
#Eca5dat<-as.data.frame(Eca5dat)
EcaAddat<-cbind(EcaAd.yall, Eca.18, Ad.18)
#EcaAddat<-as.data.frame(EcaAddat)

Fullmeandat<-rbind(Eca4dat, Eca5dat, EcaAddat, Ech4dat, Ech5dat, EchAddat, Mem4dat, Mem5dat, MemAddat, Met4dat, Met5dat, MetAddat, Pol4dat, Pol5dat, PolAddat)

Fullmeandat<-as.data.frame(Fullmeandat)
Fullmeandat$V1<-as.numeric(Fullmeandat$V1)
Fullmeandat$V2<-as.numeric(Fullmeandat$V2)
Fullmeandat$Ins4.18<-as.factor(Fullmeandat$Ins4.18)
colors.for.dplot.noAet<-c("lightskyblue3", "plum4", "goldenrod2", "sienna3", "maroon")

#Ins4.18<-c("Ins4", "Ins5", "Ad", levels = c("Ins4", "Ins5", "Ad"))

Species<-c("Entylia carinata","Ennya chrysura", "Membracis mexicana","Metheisa lucillodes", "Polyglypta costata")
names(Species) <- c("Eca", "Ech", "Mem", "Met", "Pol")
Stages<-c("Adult", "Instar 5", "Instar 4")
names(Stages)<-c("Ad", "Ins5",  "Ins4")
Fullmeandat$Ins4.18<-factor(Fullmeandat$Ins4.18, levels = c("Ad", "Ins5",  "Ins4"))

tiff("Dorsal Mean Shapes Plot.tiff", units="in", width=8, height=3.5, res=300)
mean.shapes.plot<-ggplot(data=Fullmeandat, aes(x=V1, y=V2, colour=Eca.18))+geom_point()
mean.shapes.plot +  
  facet_grid(Ins4.18~Eca.18, labeller = labeller(Eca.18=Species, Ins4.18=Stages)) + theme_bw()+theme(axis.title.x=element_blank(), axis.title.y=element_blank(), strip.text.x = element_text(face = "italic"), legend.position = "none")+ scale_color_manual(values=colors.for.dplot.noAet)
dev.off()


