library(picante)
library(phytools)
library(phangorn)
library(phylometrics)
library(geiger)
library(plyr)

source("TippyMetrics_functions.R")

phy<-read.tree("BEAST_MCC_SolanaceaeTree.tre")
phy<-drop.tip(Tree,tip=c("Convolvulus_arvensis","Evolvulus_glomeratus","Dinetus_truncatus","Montinia_caryophyllacea"))

# Read in character state information
States <- read.csv("RedStates.csv", header=FALSE, as.is=TRUE)

#read vector of tip labels with trait state 1
RedStates <- read.csv("RedStates_only.csv", header=FALSE, as.is=TRUE)

#resorting character state data to match tip order
MatchedTips <- match(Tree$tip.label, States[,1])
SortedData <- States[,][MatchedTips,]

SortedStates <- SortedData[,2]
names(SortedStates) <- SortedData[,1]

### Testing significance following Bromham et al. 2016 ###

		ParsimonyScore(SortedStates, phy) -> ParsimonyScorestat
		treestat(phy, RedStates[,1], SortedStates, func=ParsimonyScore, traitevol="TBM", alternative="greater", a=1000, simplify=F) 

		MNTD(SortedStates, phy) -> MNTDstat
		treestat(phy, RedStates[,1], SortedStates, func=MNTD, traitevol="TBM", alternative="greater", a=1000, simplify=F) 

		fpd(SortedStates, phy) -> FritzPurvisD #0.5546156=red, 0.5509395=nonred
		treestat(phy, RedStates[,1], SortedStates, func=fpd, traitevol="TBM", alternative="greater", a=1000) 

		meanBL(SortedStates, phy) -> meanBLstat
		treestat(phy, RedStates[,1], SortedStates, func=meanBL, traitevol="TBM", alternative="less", a=1000, simplify=F) 

		SlopeASR(SortedStates, phy) -> SlopeASRstat
		treestat(phy, RedStates[,1], SortedStates, func=SlopeASR, traitevol="TBM", alternative="greater", a=1000, simplify=F) 

		meanCS(SortedStates, phy) -> meanCSstat
		treestat(phy, RedStates[,1], SortedStates, func=meanCS, traitevol="TBM", alternative="less", a=1000, simplify=F) 

		noto(SortedStates, phy) -> NtipsOrigin #1.259259=red, 262.4=nonred
		treestat(phy, RedStates[,1], SortedStates, func=noto, traitevol="TBM", alternative="less", a=1000, simplify=F)




		

