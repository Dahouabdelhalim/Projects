library(venn)
library(VennDiagram)
library(ggplot2)
library(ggpolypath)

#setwd
DeData<-read.csv("DE Genes Full Subset for Venn- Final.csv", header=T)

# Compared to Leg
AL<-DeData$Anterior.Leg[1:200]
PL<-DeData$Posterior.Leg[1:529]
LW<-DeData$Wing.Leg[1:493]

# Compare to Wing
AW<-DeData$Anterior.Wing[1:191]
PW<-DeData$Posterior.Wing[1:123]
LW<-DeData$Wing.Leg[1:493]

# comparison with legs
AL.new<-unique(AL)
AL.newdf<-as.data.frame(AL.new)
PL.new<-unique(PL)
PL.newdf<-as.data.frame(PL.new)
LW.new<-unique(LW)
LW.newdf<-as.data.frame(LW.new)
leg.list<-list(AL.new, PL.new, LW.new)

# comparison with wings
AW.new<-unique(AW)
AW.newdf<-as.data.frame(AW.new)
PW.new<-unique(PW)
PW.newdf<-as.data.frame(PW.new)
wing.list<-list(AW.new, PW.new, LW.new)


colors.venn<-c("coral1", "lavenderblush", "paletourquoise1")
## All tissues compared to wings
tiff("Relative to Wing DE Tissues Venn Full Dataset.tiff", units="in", width=8, height=8, res=250)
venn(wing.list, ilcs = 2, sncs = 2, box=FALSE, zcolor="#ffdd77, #bb2020, #1188cc", snames = "Anterior-Wing, Posterior-Wing, Leg-Wing")
dev.off()

tiff("Relative to Leg DE Tissues Venn Full Dataset.tiff", units="in", width=8, height=8, res=250)
venn(leg.list, ilcs = 2, sncs = 2, box=FALSE, zcolor = "#ffdd77, #bb2020, #1188cc", snames = "Anterior-Leg, Posterior-Leg, Wing-Leg")
dev.off()

