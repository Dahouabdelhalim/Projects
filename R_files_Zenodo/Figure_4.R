
### FIGURE 4A, file: FD.csv ###
library(vegan)

FD.compound<-oenotheraFD[,3:65] 
FD.categories<-oenotheraFD[,1:2] 

FD.mds <- metaMDS(oenotheraFD.compound, distance = "bray", autotransform = FALSE)
FD.mds

stressplot(FD.mds)

# ANOSIM
anosim.FD<-anosim(FD.compound, FD$species, distance = "bray", permutations = 999)
anosim.FD


### FIGURE 4B, file: DC.csv ###
DC.compound<-oenotheraDC[,3:65] 
DC.categories<-oenotheraDC[,1:2] 

DC.mds <- metaMDS(oenotheraDC.compound, distance = "bray", autotransform = FALSE)
DC.mds

stressplot(DC.mds)

# ANOSIM
anosim.DC<-anosim(DC.compound, DC$treatment, distance = "bray", permutations = 999)
anosim.DC
