#load packages
library(ape)
library(phytools)
library(geiger)
library(nlme)
library(diversitree)
library(hisse)
library(dplyr)
library(caper)

### Note that dryad data was uploaded as .xlsx files, so to run this code these files need to be converted back to .csv files or the code needs to be slightly altered to call .xlsx files in
### Note that when calling .tre files, comments may need to be removed from the .tre file to load the tree. Comments are denoted with [ ].

### call Scarabeoidea data and tree
mytree.s <- read.nexus("~/Desktop/Data/Dataset_S3.tre")
mydata.s <- read.csv("~/Desktop/Data/Dataset_S4.csv",row.names=1)
### remove tips that we do not need
comparison.s <- name.check(phy=mytree.s,data=mydata.s)
mytree.s <- drop.tip(mytree.s,comparison.s$tree_not_data)
###align data with phylogeny
mydata.s <- mydata.s[mytree.s$tip.label,]
### visually confirm that subclades are monophyletic
cols<-setNames(c("black","blue","red", "green", "orange","yellow","purple","aliceblue","bisque","brown","cadetblue"),c("Trogidae","Ochodaeidae","Lucanidae", "Hybosoridae","Glaresidae","Glaphyridae", "Geotrupidae","Eremazidae","Bolboracertidae & Passalidae", "Aphodiinae & Aegialiinae & Scarabaeinae", "Rutelinae & Dynastinae & Melolonthinae & Cetoniinae & Orphninae & Pachypodidae &Euchiridae"))
plotTree(mytree.s, ftype="off", type="fan")
tiplabels(pie=to.matrix(mydata.s$Subclade, sort(unique(mydata.s$Subclade))), piecol = cols, cex=0.5)

### call Phasmatodea data and tree
mytree.p <- read.nexus("~/Desktop/Data/Dataset_S5.tre")
mydata.p <- read.csv("~/Desktop/Data/Dataset_S6.csv",row.names=1)
### remove tips that we do not need
comparison.p <- name.check(phy=mytree.p,data=mydata.p)
mytree.p <- drop.tip(mytree.p,comparison.p$tree_not_data)
###align data with phylogeny
mydata.p <- mydata.p[mytree.p$tip.label,]
### visually confirm that subclades are monophyletic, white=unassigned
cols<-setNames(c("white", "black","blue","darkseagreen4", "green","yellow","purple","aquamarine","darkorchid4","brown","cadetblue", "azure4","brown1", "goldenrod4","coral","cyan", "darkolivegreen1", "darkgreen","deeppink","deepskyblue4"),c("","Achriopterini", "Agathemeridae", "Anisacanthidae", "Antongiliinae", "Aschiphasmatinae", "Cladomorphinae", "Clitumninae", "Damasippoididae", "Diapheromerinae", "Heteropteryginae", "Lanceocercata", "Lonchodinae", "Necrosciinae", "Palophinae", "Pharnaciini", "Phylliinae", "Pseudophasmatinae", "Stephanacridini", "Timematidae"))
plotTree(mytree.p, ftype="off", type="fan")
tiplabels(pie=to.matrix(mydata.p$Subclade, sort(unique(mydata.p$Subclade))), piecol = cols, cex=0.5)
### white in this phylogeny means unassigned subclade

### call Coreidae data and tree
mytree.c <- read.nexus("~/Desktop/Data/Dataset_S1.tre")
mydata.c <- read.csv("~/Desktop/Data/Dataset_S2.csv",row.names=1)
### remove tips that we do not need
comparison.c <- name.check(phy=mytree.c,data=mydata.c)
mytree.c <- drop.tip(mytree.c,comparison.c$tree_not_data)
###align data with phylogeny
mydata.c <- mydata.c[mytree.c$tip.label,]

####################
### SSE Analyses ###
####################

#reorganize data for analysis
mydata.s.fix <- tibble::rownames_to_column(mydata.s, "Tip")
mydata.p.fix <- tibble::rownames_to_column(mydata.p, "Tip")
mydata.c.fix <- tibble::rownames_to_column(mydata.c, "Tip")
mydata.s.sse <- mydata.s.fix[,c(1,5)]
mydata.p.sse <- mydata.p.fix[,c(1,5)]
mydata.c.sse <- mydata.c.fix[,c(1,5)]

### Assgin X = 1 for liberal coding scheme and X = 0 for conservative coding scheme
mydata.s.sse$SSW[mydata.s.sse$SSW=="X"] <- 1
mydata.p.sse$SSW[mydata.p.sse$SSW=="X"] <- 1
mydata.c.sse$SSW[mydata.c.sse$SSW=="X"] <- 1
mydata.s.sse$SSW <- droplevels(mydata.s.sse$SSW)
mydata.p.sse$SSW <- droplevels(mydata.p.sse$SSW)
mydata.c.sse$SSW <- droplevels(mydata.c.sse$SSW)

### Phasmatodea ###

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.p <-  hisse(mytree.p, mydata.p.sse, f=c(0.075,0.075), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.p


#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.p <-  hisse(mytree.p, mydata.p.sse, f=c(0.075,0.075), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.p

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.p <-  hisse(mytree.p, mydata.p.sse, f=c(0.075,0.075), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.p

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.p <-  hisse(mytree.p, mydata.p.sse, f=c(0.075,0.075), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.p

#CID-4 model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.p <- hisse(mytree.p, mydata.p.sse, f=c(0.075,0.075), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.p

### Lanceocerata ###
### call in Lanceocerata tree and reorganize Phasmatodea data for analysis 
mytree.p.Lan <- read.newick("~/Desktop/Data/Dataset_S7.tre")
mydata.p.Lan.sse <- mydata.p.sse[match(mytree.p.Lan$tip.label, mydata.p.sse$Tip),]
mydata.p.Lan.sse$SSW <- as.factor(mydata.p.Lan.sse$SSW)

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.p.Lan <-  hisse(mytree.p.Lan, mydata.p.Lan.sse, f=c(0.23,0.23), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.p.Lan

#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.p.Lan <-  hisse(mytree.p.Lan, mydata.p.Lan.sse, f=c(0.23,0.23), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.p.Lan

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.p.Lan <-  hisse(mytree.p.Lan, mydata.p.Lan.sse, f=c(0.23,0.23), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.p.Lan

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.p.Lan <-  hisse(mytree.p.Lan, mydata.p.Lan.sse, f=c(0.23,0.23), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.p.Lan

#CID-4 Model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.p.Lan <- hisse(mytree.p.Lan, mydata.p.Lan.sse, f=c(0.23,0.23), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.p.Lan

### Lonchodinae ###
### call in Lonchodinae tree and reorganize Phasmatodea data for analysis 
mytree.p.Lonch <- read.newick("~/Desktop/Data/Dataset_S8.tre")
mydata.p.Lonch.sse <- mydata.p.sse[match(mytree.p.Lonch$tip.label, mydata.p.sse$Tip),]
mydata.p.Lonch.sse$SSW <- as.factor(mydata.p.Lonch.sse$SSW)

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.p.Lonch <-  hisse(mytree.p.Lonch, mydata.p.Lonch.sse, f=c(0.09,0.09), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.p.Lonch

#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.p.Lonch <-  hisse(mytree.p.Lonch, mydata.p.Lonch.sse, f=c(0.09,0.09), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.p.Lonch

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.p.Lonch <-  hisse(mytree.p.Lonch, mydata.p.Lonch.sse, f=c(0.09,0.09), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.p.Lonch

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.p.Lonch <-  hisse(mytree.p.Lonch, mydata.p.Lonch.sse, f=c(0.09,0.09), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.p.Lonch

#CID-4 Model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.p.Lonch <- hisse(mytree.p.Lonch, mydata.p.Lonch.sse, f=c(0.09,0.09), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.p.Lonch

### Lonchodinae Simulated Data ###
### call in Lonchodinae tree and simulated data for analysis 
mytree.p.Lonch.sim.sse <- read.newick("~/Desktop/Data/Dataset_S8.tre")
mydata.p.Lonch.sim.sse <- read.csv("~/Desktop/Data/Dataset_S12.csv")
mydata.p.Lonch.sim.sse$SSW <- as.factor(as.numeric((mydata.p.Lonch.sim.sse$SSW)))

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.p.Lonch.sim <-  hisse(mytree.p.Lonch.sim.sse, mydata.p.Lonch.sim.sse, f=c(0.09,0.09), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "raw")
BiSSE.p.Lonch.sim

#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.p.Lonch.sim <-  hisse(mytree.p.Lonch.sim.sse, mydata.p.Lonch.sim.sse, f=c(0.09,0.09), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.p.Lonch.sim

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.p.Lonch.sim <-  hisse(mytree.p.Lonch.sim.sse, mydata.p.Lonch.sim.sse, f=c(0.09,0.09), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.p.Lonch.sim

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.p.Lonch.sim <-  hisse(mytree.p.Lonch.sim.sse, mydata.p.Lonch.sim.sse, f=c(0.09,0.09), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.p.Lonch.sim

#CID-4 Model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.p.Lonch.sim <- hisse(mytree.p.Lonch.sim.sse, mydata.p.Lonch.sim.sse, f=c(0.09,0.09), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.p.Lonch.sim


### Scarabeoidea ###

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.s <-  hisse(mytree.s, mydata.s.sse, f=c(0.005,0.005), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.s


#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.s <-  hisse(mytree.s, mydata.s.sse, f=c(0.005,0.005), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.s

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.s <-  hisse(mytree.s, mydata.s.sse, f=c(0.005,0.005), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.s

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.s <-  hisse(mytree.s, mydata.s.sse, f=c(0.005,0.005), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.s

#CID-4 Model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.s <- hisse(mytree.s, mydata.s.sse, f=c(0.005,0.005), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.s

### Coreidae ###

#BiSSE Model (implemented in HiSSE)
turnover.anc = c(1,2,0,0)
eps.anc = c(1,2,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.c <-  hisse(mytree.c, mydata.c.sse, f=c(0.02,0.02), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.c

#BiSSE Null Model (implemented in HiSSE)
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)
trans.rates.bisse = TransMatMaker(hidden.states=FALSE)
BiSSE.null.c <-  hisse(mytree.c, mydata.c.sse, f=c(0.02,0.02), hidden.states=FALSE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.bisse, output.type = "net.div")
BiSSE.null.c

#CID-2 Model
turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)
#full 8 transition model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID2.c <-  hisse(mytree.c, mydata.c.sse, f=c(0.02,0.02), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal,output.type="net.div")
CID2.c

#Full HiSSE Model
turnover.anc = c(1,2,3,4)
eps.anc = c(1,2,3,4)
#setting up transition rate matrix
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates
#removing transitions from obvserved to hidden traits
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual
HiSSE.c <-  hisse(mytree.c, mydata.c.sse, f=c(0.02,0.02), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual,output.type="net.div")
HiSSE.c

#CID-4 Model
trans.rates = TransMatMaker(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
#setting all transition rates equal
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal
CID4.c <- hisse(mytree.c, mydata.c.sse, f=c(0.02,0.02), turnover.anc=rep(c(1,2,3,4),2),eps.anc=rep(c(1,2,3,4),2), trans.rate=trans.rates.nodual.allequal)
CID4.c


####################
### MS Analyses ###
####################

### call Scarabeoidea MS data and tree
data.ms.s <-  read.csv("~/Desktop/Data/Dataset_S10.csv")
tree.ms.s <-  read.nexus("~/Desktop/Data/Dataset_S3.tre")
speciation.ms.s <- comparative.data(phy = tree.ms.s, data = data.ms.s, names.col = Species_ID.3, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

### check for taxonomic bias ###
summary(lm(Described.speces~Species.in.tree, data=data.ms.s ))

#PGLS regressions testing relationships between diversification rates and frequency of sexually selected weapons among species 
summary(pgls(E.0 ~ proportion.weapons.1, data = speciation.ms.s, lambda ="ML"))
summary(pgls(E.5 ~ proportion.weapons.1, data = speciation.ms.s, lambda ="ML"))
summary(pgls(E.9 ~ proportion.weapons.1, data = speciation.ms.s, lambda ="ML"))
#to get results for conservative coding scheme simply replace proportion.weapons.1 with proportion.weapons.0 

###Scarabeoidea simulated replicates
model.1 <- pgls(E.5 ~ sim.1, data = speciation.ms.s, lambda =0.001)
model.2 <- pgls(E.5 ~ sim.2, data = speciation.ms.s, lambda =0.001)
model.3 <- pgls(E.5 ~ sim.3, data = speciation.ms.s, lambda =0.001)
model.4 <- pgls(E.5 ~ sim.4, data = speciation.ms.s, lambda =0.001)
model.5 <- pgls(E.5 ~ sim.5, data = speciation.ms.s, lambda =0.001)
model.6 <- pgls(E.5 ~ sim.6, data = speciation.ms.s, lambda =0.001)
model.7 <- pgls(E.5 ~ sim.7, data = speciation.ms.s, lambda =0.001)
model.8 <- pgls(E.5 ~ sim.8, data = speciation.ms.s, lambda =0.001)
model.9 <- pgls(E.5 ~ sim.9, data = speciation.ms.s, lambda =0.001)
model.10 <- pgls(E.5 ~ sim.10, data = speciation.ms.s, lambda =0.001)
model.11 <- pgls(E.5 ~ sim.11, data = speciation.ms.s, lambda =0.001)
model.12 <- pgls(E.5 ~ sim.12, data = speciation.ms.s, lambda =0.001)
model.13 <- pgls(E.5 ~ sim.13, data = speciation.ms.s, lambda =0.001)
model.14 <- pgls(E.5 ~ sim.14, data = speciation.ms.s, lambda =0.001)
model.15 <- pgls(E.5 ~ sim.15, data = speciation.ms.s, lambda =0.001)
model.16 <- pgls(E.5 ~ sim.16, data = speciation.ms.s, lambda =0.001)
model.17 <- pgls(E.5 ~ sim.17, data = speciation.ms.s, lambda =0.001)
model.18 <- pgls(E.5 ~ sim.18, data = speciation.ms.s, lambda =0.001)
model.19 <- pgls(E.5 ~ sim.19, data = speciation.ms.s, lambda =0.001)
model.20 <- pgls(E.5 ~ sim.20, data = speciation.ms.s, lambda =0.001)
model.21 <- pgls(E.5 ~ sim.21, data = speciation.ms.s, lambda =0.001)
model.22 <- pgls(E.5 ~ sim.22, data = speciation.ms.s, lambda =0.001)
model.23 <- pgls(E.5 ~ sim.23, data = speciation.ms.s, lambda =0.001)
model.24 <- pgls(E.5 ~ sim.24, data = speciation.ms.s, lambda =0.001)
model.25 <- pgls(E.5 ~ sim.25, data = speciation.ms.s, lambda =0.001)
model.26 <- pgls(E.5 ~ sim.26, data = speciation.ms.s, lambda =0.001)
model.27 <- pgls(E.5 ~ sim.27, data = speciation.ms.s, lambda =0.001)
model.28 <- pgls(E.5 ~ sim.28, data = speciation.ms.s, lambda =0.001)
model.29 <- pgls(E.5 ~ sim.29, data = speciation.ms.s, lambda =0.001)
model.30 <- pgls(E.5 ~ sim.30, data = speciation.ms.s, lambda =0.001)
model.31 <- pgls(E.5 ~ sim.31, data = speciation.ms.s, lambda =0.001)
model.32 <- pgls(E.5 ~ sim.32, data = speciation.ms.s, lambda =0.001)
model.33 <- pgls(E.5 ~ sim.33, data = speciation.ms.s, lambda =0.001)
model.34 <- pgls(E.5 ~ sim.34, data = speciation.ms.s, lambda =0.001)
model.35 <- pgls(E.5 ~ sim.35, data = speciation.ms.s, lambda =0.001)
model.36 <- pgls(E.5 ~ sim.36, data = speciation.ms.s, lambda =0.001)
model.37 <- pgls(E.5 ~ sim.37, data = speciation.ms.s, lambda =0.001)
model.38 <- pgls(E.5 ~ sim.38, data = speciation.ms.s, lambda =0.001)
model.39 <- pgls(E.5 ~ sim.39, data = speciation.ms.s, lambda =0.001)
model.40 <- pgls(E.5 ~ sim.40, data = speciation.ms.s, lambda =0.001)
model.41 <- pgls(E.5 ~ sim.41, data = speciation.ms.s, lambda =0.001)
model.42 <- pgls(E.5 ~ sim.42, data = speciation.ms.s, lambda =0.001)
model.43 <- pgls(E.5 ~ sim.43, data = speciation.ms.s, lambda =0.001)
model.44 <- pgls(E.5 ~ sim.44, data = speciation.ms.s, lambda =0.001)
model.45 <- pgls(E.5 ~ sim.45, data = speciation.ms.s, lambda =0.001)
model.46 <- pgls(E.5 ~ sim.46, data = speciation.ms.s, lambda =0.001)
model.47 <- pgls(E.5 ~ sim.47, data = speciation.ms.s, lambda =0.001)
model.48 <- pgls(E.5 ~ sim.48, data = speciation.ms.s, lambda =0.001)
model.49 <- pgls(E.5 ~ sim.49, data = speciation.ms.s, lambda =0.001)
model.50 <- pgls(E.5 ~ sim.50, data = speciation.ms.s, lambda =0.001)
model.51 <- pgls(E.5 ~ sim.51, data = speciation.ms.s, lambda =0.001)
model.52 <- pgls(E.5 ~ sim.52, data = speciation.ms.s, lambda =0.001)
model.53 <- pgls(E.5 ~ sim.53, data = speciation.ms.s, lambda =0.001)
model.54 <- pgls(E.5 ~ sim.54, data = speciation.ms.s, lambda =0.001)
model.55 <- pgls(E.5 ~ sim.55, data = speciation.ms.s, lambda =0.001)
model.56 <- pgls(E.5 ~ sim.56, data = speciation.ms.s, lambda =0.001)
model.57 <- pgls(E.5 ~ sim.57, data = speciation.ms.s, lambda =0.001)
model.58 <- pgls(E.5 ~ sim.58, data = speciation.ms.s, lambda =0.001)
model.59 <- pgls(E.5 ~ sim.59, data = speciation.ms.s, lambda =0.001)
model.60 <- pgls(E.5 ~ sim.60, data = speciation.ms.s, lambda =0.001)
model.61 <- pgls(E.5 ~ sim.61, data = speciation.ms.s, lambda =0.001)
model.62 <- pgls(E.5 ~ sim.62, data = speciation.ms.s, lambda =0.001)
model.63 <- pgls(E.5 ~ sim.63, data = speciation.ms.s, lambda =0.001)
model.64 <- pgls(E.5 ~ sim.64, data = speciation.ms.s, lambda =0.001)
model.65 <- pgls(E.5 ~ sim.65, data = speciation.ms.s, lambda =0.001)
model.66 <- pgls(E.5 ~ sim.66, data = speciation.ms.s, lambda =0.001)
model.67 <- pgls(E.5 ~ sim.67, data = speciation.ms.s, lambda =0.001)
model.68 <- pgls(E.5 ~ sim.68, data = speciation.ms.s, lambda =0.001)
model.69 <- pgls(E.5 ~ sim.69, data = speciation.ms.s, lambda =0.001)
model.70 <- pgls(E.5 ~ sim.70, data = speciation.ms.s, lambda =0.001)
model.71 <- pgls(E.5 ~ sim.71, data = speciation.ms.s, lambda =0.001)
model.72 <- pgls(E.5 ~ sim.72, data = speciation.ms.s, lambda =0.001)
model.73 <- pgls(E.5 ~ sim.73, data = speciation.ms.s, lambda =0.001)
model.74 <- pgls(E.5 ~ sim.74, data = speciation.ms.s, lambda =0.001)
model.75 <- pgls(E.5 ~ sim.75, data = speciation.ms.s, lambda =0.001)
model.76 <- pgls(E.5 ~ sim.76, data = speciation.ms.s, lambda =0.001)
model.77 <- pgls(E.5 ~ sim.77, data = speciation.ms.s, lambda =0.001)
model.78 <- pgls(E.5 ~ sim.78, data = speciation.ms.s, lambda =0.001)
model.79 <- pgls(E.5 ~ sim.79, data = speciation.ms.s, lambda =0.001)
model.80 <- pgls(E.5 ~ sim.80, data = speciation.ms.s, lambda =0.001)
model.81 <- pgls(E.5 ~ sim.81, data = speciation.ms.s, lambda =0.001)
model.82 <- pgls(E.5 ~ sim.82, data = speciation.ms.s, lambda =0.001)
model.83 <- pgls(E.5 ~ sim.83, data = speciation.ms.s, lambda =0.001)
model.84 <- pgls(E.5 ~ sim.84, data = speciation.ms.s, lambda =0.001)
model.85 <- pgls(E.5 ~ sim.85, data = speciation.ms.s, lambda =0.001)
model.86 <- pgls(E.5 ~ sim.86, data = speciation.ms.s, lambda =0.001)
model.87 <- pgls(E.5 ~ sim.87, data = speciation.ms.s, lambda =0.001)
model.88 <- pgls(E.5 ~ sim.88, data = speciation.ms.s, lambda =0.001)
model.89 <- pgls(E.5 ~ sim.89, data = speciation.ms.s, lambda =0.001)
model.90 <- pgls(E.5 ~ sim.90, data = speciation.ms.s, lambda =0.001)
model.91 <- pgls(E.5 ~ sim.91, data = speciation.ms.s, lambda =0.001)
model.92 <- pgls(E.5 ~ sim.92, data = speciation.ms.s, lambda =0.001)
model.93 <- pgls(E.5 ~ sim.93, data = speciation.ms.s, lambda =0.001)
model.94 <- pgls(E.5 ~ sim.94, data = speciation.ms.s, lambda =0.001)
model.95 <- pgls(E.5 ~ sim.95, data = speciation.ms.s, lambda =0.001)
model.96 <- pgls(E.5 ~ sim.96, data = speciation.ms.s, lambda =0.001)
model.97 <- pgls(E.5 ~ sim.97, data = speciation.ms.s, lambda =0.001)
model.98 <- pgls(E.5 ~ sim.98, data = speciation.ms.s, lambda =0.001)
model.99 <- pgls(E.5 ~ sim.99, data = speciation.ms.s, lambda =0.001)
model.100 <- pgls(E.5 ~ sim.100, data = speciation.ms.s, lambda =0.001)

p.s <- c(summary(model.1)$coefficients[2,4],summary(model.2)$coefficients[2,4],summary(model.3)$coefficients[2,4],summary(model.4)$coefficients[2,4],summary(model.5)$coefficients[2,4],summary(model.6)$coefficients[2,4],summary(model.7)$coefficients[2,4],summary(model.8)$coefficients[2,4],summary(model.9)$coefficients[2,4],summary(model.10)$coefficients[2,4],summary(model.11)$coefficients[2,4],summary(model.12)$coefficients[2,4],summary(model.13)$coefficients[2,4],summary(model.14)$coefficients[2,4],summary(model.15)$coefficients[2,4],summary(model.16)$coefficients[2,4],summary(model.17)$coefficients[2,4],summary(model.18)$coefficients[2,4],summary(model.19)$coefficients[2,4],summary(model.20)$coefficients[2,4],summary(model.21)$coefficients[2,4],summary(model.22)$coefficients[2,4],summary(model.23)$coefficients[2,4],summary(model.24)$coefficients[2,4],summary(model.25)$coefficients[2,4],summary(model.26)$coefficients[2,4],summary(model.27)$coefficients[2,4],summary(model.28)$coefficients[2,4],summary(model.29)$coefficients[2,4],summary(model.30)$coefficients[2,4],summary(model.31)$coefficients[2,4],summary(model.32)$coefficients[2,4],summary(model.33)$coefficients[2,4],summary(model.34)$coefficients[2,4],summary(model.35)$coefficients[2,4],summary(model.36)$coefficients[2,4],summary(model.37)$coefficients[2,4],summary(model.38)$coefficients[2,4],summary(model.39)$coefficients[2,4],summary(model.40)$coefficients[2,4],summary(model.41)$coefficients[2,4],summary(model.42)$coefficients[2,4],summary(model.43)$coefficients[2,4],summary(model.44)$coefficients[2,4],summary(model.45)$coefficients[2,4],summary(model.46)$coefficients[2,4],summary(model.47)$coefficients[2,4],summary(model.48)$coefficients[2,4],summary(model.49)$coefficients[2,4],summary(model.50)$coefficients[2,4],summary(model.51)$coefficients[2,4],summary(model.52)$coefficients[2,4],summary(model.53)$coefficients[2,4],summary(model.54)$coefficients[2,4],summary(model.55)$coefficients[2,4],summary(model.56)$coefficients[2,4],summary(model.57)$coefficients[2,4],summary(model.58)$coefficients[2,4],summary(model.59)$coefficients[2,4],summary(model.60)$coefficients[2,4],summary(model.61)$coefficients[2,4],summary(model.62)$coefficients[2,4],summary(model.63)$coefficients[2,4],summary(model.64)$coefficients[2,4],summary(model.65)$coefficients[2,4],summary(model.66)$coefficients[2,4],summary(model.67)$coefficients[2,4],summary(model.68)$coefficients[2,4],summary(model.69)$coefficients[2,4],summary(model.70)$coefficients[2,4],summary(model.71)$coefficients[2,4],summary(model.72)$coefficients[2,4],summary(model.73)$coefficients[2,4],summary(model.74)$coefficients[2,4],summary(model.75)$coefficients[2,4],summary(model.76)$coefficients[2,4],summary(model.77)$coefficients[2,4],summary(model.78)$coefficients[2,4],summary(model.79)$coefficients[2,4],summary(model.80)$coefficients[2,4],summary(model.81)$coefficients[2,4],summary(model.82)$coefficients[2,4],summary(model.83)$coefficients[2,4],summary(model.84)$coefficients[2,4],summary(model.85)$coefficients[2,4],summary(model.86)$coefficients[2,4],summary(model.87)$coefficients[2,4],summary(model.88)$coefficients[2,4],summary(model.89)$coefficients[2,4],summary(model.90)$coefficients[2,4],summary(model.91)$coefficients[2,4],summary(model.92)$coefficients[2,4],summary(model.93)$coefficients[2,4],summary(model.94)$coefficients[2,4],summary(model.95)$coefficients[2,4],summary(model.96)$coefficients[2,4],summary(model.97)$coefficients[2,4],summary(model.98)$coefficients[2,4],summary(model.99)$coefficients[2,4],summary(model.100)$coefficients[2,4])
slope.s <- c(summary(model.1)$coefficients[2,1],summary(model.2)$coefficients[2,1],summary(model.3)$coefficients[2,1],summary(model.4)$coefficients[2,1],summary(model.5)$coefficients[2,1],summary(model.6)$coefficients[2,1],summary(model.7)$coefficients[2,1],summary(model.8)$coefficients[2,1],summary(model.9)$coefficients[2,1],summary(model.10)$coefficients[2,1],summary(model.11)$coefficients[2,1],summary(model.12)$coefficients[2,1],summary(model.13)$coefficients[2,1],summary(model.14)$coefficients[2,1],summary(model.15)$coefficients[2,1],summary(model.16)$coefficients[2,1],summary(model.17)$coefficients[2,1],summary(model.18)$coefficients[2,1],summary(model.19)$coefficients[2,1],summary(model.20)$coefficients[2,1],summary(model.21)$coefficients[2,1],summary(model.22)$coefficients[2,1],summary(model.23)$coefficients[2,1],summary(model.24)$coefficients[2,1],summary(model.25)$coefficients[2,1],summary(model.26)$coefficients[2,1],summary(model.27)$coefficients[2,1],summary(model.28)$coefficients[2,1],summary(model.29)$coefficients[2,1],summary(model.30)$coefficients[2,1],summary(model.31)$coefficients[2,1],summary(model.32)$coefficients[2,1],summary(model.33)$coefficients[2,1],summary(model.34)$coefficients[2,1],summary(model.35)$coefficients[2,1],summary(model.36)$coefficients[2,1],summary(model.37)$coefficients[2,1],summary(model.38)$coefficients[2,1],summary(model.39)$coefficients[2,1],summary(model.40)$coefficients[2,1],summary(model.41)$coefficients[2,1],summary(model.42)$coefficients[2,1],summary(model.43)$coefficients[2,1],summary(model.44)$coefficients[2,1],summary(model.45)$coefficients[2,1],summary(model.46)$coefficients[2,1],summary(model.47)$coefficients[2,1],summary(model.48)$coefficients[2,1],summary(model.49)$coefficients[2,1],summary(model.50)$coefficients[2,1],summary(model.51)$coefficients[2,1],summary(model.52)$coefficients[2,1],summary(model.53)$coefficients[2,1],summary(model.54)$coefficients[2,1],summary(model.55)$coefficients[2,1],summary(model.56)$coefficients[2,1],summary(model.57)$coefficients[2,1],summary(model.58)$coefficients[2,1],summary(model.59)$coefficients[2,1],summary(model.60)$coefficients[2,1],summary(model.61)$coefficients[2,1],summary(model.62)$coefficients[2,1],summary(model.63)$coefficients[2,1],summary(model.64)$coefficients[2,1],summary(model.65)$coefficients[2,1],summary(model.66)$coefficients[2,1],summary(model.67)$coefficients[2,1],summary(model.68)$coefficients[2,1],summary(model.69)$coefficients[2,1],summary(model.70)$coefficients[2,1],summary(model.71)$coefficients[2,1],summary(model.72)$coefficients[2,1],summary(model.73)$coefficients[2,1],summary(model.74)$coefficients[2,1],summary(model.75)$coefficients[2,1],summary(model.76)$coefficients[2,1],summary(model.77)$coefficients[2,1],summary(model.78)$coefficients[2,1],summary(model.79)$coefficients[2,1],summary(model.80)$coefficients[2,1],summary(model.81)$coefficients[2,1],summary(model.82)$coefficients[2,1],summary(model.83)$coefficients[2,1],summary(model.84)$coefficients[2,1],summary(model.85)$coefficients[2,1],summary(model.86)$coefficients[2,1],summary(model.87)$coefficients[2,1],summary(model.88)$coefficients[2,1],summary(model.89)$coefficients[2,1],summary(model.90)$coefficients[2,1],summary(model.91)$coefficients[2,1],summary(model.92)$coefficients[2,1],summary(model.93)$coefficients[2,1],summary(model.94)$coefficients[2,1],summary(model.95)$coefficients[2,1],summary(model.96)$coefficients[2,1],summary(model.97)$coefficients[2,1],summary(model.98)$coefficients[2,1],summary(model.99)$coefficients[2,1],summary(model.100)$coefficients[2,1])
replicate.s <- 1:100
overview.s <- data.frame(cbind(replicate.s, slope.s, p.s))
sig.s <- subset(overview.s, p.s < 0.05) 
nrow(sig.s) 


### call Phasmatodea MS data and tree
data.ms.p <-  read.csv("~/Desktop/Data/Dataset_S9.csv")
tree.ms.p <-  read.nexus("~/Desktop/Data/Dataset_S5.tre")
speciation.ms.p <- comparative.data(phy = tree.ms.p, data = data.ms.p, names.col = Species_ID.3, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)

### check for taxonomic bias ###
summary(lm(Described.species~Species.in.tree, data=data.ms.p))

#PGLS regressions testing relationships between diversification rates and frequency of sexually selected weapons among species 
summary(pgls(E.0 ~ proportion.weapons.1, data = speciation.ms.p, lambda ="ML"))
summary(pgls(E.5 ~ proportion.weapons.1, data = speciation.ms.p, lambda ="ML"))
summary(pgls(E.9 ~ proportion.weapons.1, data = speciation.ms.p, lambda ="ML"))
#to get results for conservative coding scheme simply replace proportion.weapons.1 with proportion.weapons.0 

###Phasmatodea simulated replicates
model.1 <- pgls(E.5 ~ sim.1, data = speciation.ms.p, lambda =0.001)
model.2 <- pgls(E.5 ~ sim.2, data = speciation.ms.p, lambda =0.001)
model.3 <- pgls(E.5 ~ sim.3, data = speciation.ms.p, lambda =0.001)
model.4 <- pgls(E.5 ~ sim.4, data = speciation.ms.p, lambda =0.001)
model.5 <- pgls(E.5 ~ sim.5, data = speciation.ms.p, lambda =0.001)
model.6 <- pgls(E.5 ~ sim.6, data = speciation.ms.p, lambda =0.001)
model.7 <- pgls(E.5 ~ sim.7, data = speciation.ms.p, lambda =0.001)
model.8 <- pgls(E.5 ~ sim.8, data = speciation.ms.p, lambda =0.001)
model.9 <- pgls(E.5 ~ sim.9, data = speciation.ms.p, lambda =0.001)
model.10 <- pgls(E.5 ~ sim.10, data = speciation.ms.p, lambda =0.001)
model.11 <- pgls(E.5 ~ sim.11, data = speciation.ms.p, lambda =0.001)
model.12 <- pgls(E.5 ~ sim.12, data = speciation.ms.p, lambda =0.001)
model.13 <- pgls(E.5 ~ sim.13, data = speciation.ms.p, lambda =0.001)
model.14 <- pgls(E.5 ~ sim.14, data = speciation.ms.p, lambda =0.001)
model.15 <- pgls(E.5 ~ sim.15, data = speciation.ms.p, lambda =0.001)
model.16 <- pgls(E.5 ~ sim.16, data = speciation.ms.p, lambda =0.001)
model.17 <- pgls(E.5 ~ sim.17, data = speciation.ms.p, lambda =0.001)
model.18 <- pgls(E.5 ~ sim.18, data = speciation.ms.p, lambda =0.001)
model.19 <- pgls(E.5 ~ sim.19, data = speciation.ms.p, lambda =0.001)
model.20 <- pgls(E.5 ~ sim.20, data = speciation.ms.p, lambda =0.001)
model.21 <- pgls(E.5 ~ sim.21, data = speciation.ms.p, lambda =0.001)
model.22 <- pgls(E.5 ~ sim.22, data = speciation.ms.p, lambda =0.001)
model.23 <- pgls(E.5 ~ sim.23, data = speciation.ms.p, lambda =0.001)
model.24 <- pgls(E.5 ~ sim.24, data = speciation.ms.p, lambda =0.001)
model.25 <- pgls(E.5 ~ sim.25, data = speciation.ms.p, lambda =0.001)
model.26 <- pgls(E.5 ~ sim.26, data = speciation.ms.p, lambda =0.001)
model.27 <- pgls(E.5 ~ sim.27, data = speciation.ms.p, lambda =0.001)
model.28 <- pgls(E.5 ~ sim.28, data = speciation.ms.p, lambda =0.001)
model.29 <- pgls(E.5 ~ sim.29, data = speciation.ms.p, lambda =0.001)
model.30 <- pgls(E.5 ~ sim.30, data = speciation.ms.p, lambda =0.001)
model.31 <- pgls(E.5 ~ sim.31, data = speciation.ms.p, lambda =0.001)
model.32 <- pgls(E.5 ~ sim.32, data = speciation.ms.p, lambda =0.001)
model.33 <- pgls(E.5 ~ sim.33, data = speciation.ms.p, lambda =0.001)
model.34 <- pgls(E.5 ~ sim.34, data = speciation.ms.p, lambda =0.001)
model.35 <- pgls(E.5 ~ sim.35, data = speciation.ms.p, lambda =0.001)
model.36 <- pgls(E.5 ~ sim.36, data = speciation.ms.p, lambda =0.001)
model.37 <- pgls(E.5 ~ sim.37, data = speciation.ms.p, lambda =0.001)
model.38 <- pgls(E.5 ~ sim.38, data = speciation.ms.p, lambda =0.001)
model.39 <- pgls(E.5 ~ sim.39, data = speciation.ms.p, lambda =0.001)
model.40 <- pgls(E.5 ~ sim.40, data = speciation.ms.p, lambda =0.001)
model.41 <- pgls(E.5 ~ sim.41, data = speciation.ms.p, lambda =0.001)
model.42 <- pgls(E.5 ~ sim.42, data = speciation.ms.p, lambda =0.001)
model.43 <- pgls(E.5 ~ sim.43, data = speciation.ms.p, lambda =0.001)
model.44 <- pgls(E.5 ~ sim.44, data = speciation.ms.p, lambda =0.001)
model.45 <- pgls(E.5 ~ sim.45, data = speciation.ms.p, lambda =0.001)
model.46 <- pgls(E.5 ~ sim.46, data = speciation.ms.p, lambda =0.001)
model.47 <- pgls(E.5 ~ sim.47, data = speciation.ms.p, lambda =0.001)
model.48 <- pgls(E.5 ~ sim.48, data = speciation.ms.p, lambda =0.001)
model.49 <- pgls(E.5 ~ sim.49, data = speciation.ms.p, lambda =0.001)
model.50 <- pgls(E.5 ~ sim.50, data = speciation.ms.p, lambda =0.001)
model.51 <- pgls(E.5 ~ sim.51, data = speciation.ms.p, lambda =0.001)
model.52 <- pgls(E.5 ~ sim.52, data = speciation.ms.p, lambda =0.001)
model.53 <- pgls(E.5 ~ sim.53, data = speciation.ms.p, lambda =0.001)
model.54 <- pgls(E.5 ~ sim.54, data = speciation.ms.p, lambda =0.001)
model.55 <- pgls(E.5 ~ sim.55, data = speciation.ms.p, lambda =0.001)
model.56 <- pgls(E.5 ~ sim.56, data = speciation.ms.p, lambda =0.001)
model.57 <- pgls(E.5 ~ sim.57, data = speciation.ms.p, lambda =0.001)
model.58 <- pgls(E.5 ~ sim.58, data = speciation.ms.p, lambda =0.001)
model.59 <- pgls(E.5 ~ sim.59, data = speciation.ms.p, lambda =0.001)
model.60 <- pgls(E.5 ~ sim.60, data = speciation.ms.p, lambda =0.001)
model.61 <- pgls(E.5 ~ sim.61, data = speciation.ms.p, lambda =0.001)
model.62 <- pgls(E.5 ~ sim.62, data = speciation.ms.p, lambda =0.001)
model.63 <- pgls(E.5 ~ sim.63, data = speciation.ms.p, lambda =0.001)
model.64 <- pgls(E.5 ~ sim.64, data = speciation.ms.p, lambda =0.001)
model.65 <- pgls(E.5 ~ sim.65, data = speciation.ms.p, lambda =0.001)
model.66 <- pgls(E.5 ~ sim.66, data = speciation.ms.p, lambda =0.001)
model.67 <- pgls(E.5 ~ sim.67, data = speciation.ms.p, lambda =0.001)
model.68 <- pgls(E.5 ~ sim.68, data = speciation.ms.p, lambda =0.001)
model.69 <- pgls(E.5 ~ sim.69, data = speciation.ms.p, lambda =0.001)
model.70 <- pgls(E.5 ~ sim.70, data = speciation.ms.p, lambda =0.001)
model.71 <- pgls(E.5 ~ sim.71, data = speciation.ms.p, lambda =0.001)
model.72 <- pgls(E.5 ~ sim.72, data = speciation.ms.p, lambda =0.001)
model.73 <- pgls(E.5 ~ sim.73, data = speciation.ms.p, lambda =0.001)
model.74 <- pgls(E.5 ~ sim.74, data = speciation.ms.p, lambda =0.001)
model.75 <- pgls(E.5 ~ sim.75, data = speciation.ms.p, lambda =0.001)
model.76 <- pgls(E.5 ~ sim.76, data = speciation.ms.p, lambda =0.001)
model.77 <- pgls(E.5 ~ sim.77, data = speciation.ms.p, lambda =0.001)
model.78 <- pgls(E.5 ~ sim.78, data = speciation.ms.p, lambda =0.001)
model.79 <- pgls(E.5 ~ sim.79, data = speciation.ms.p, lambda =0.001)
model.80 <- pgls(E.5 ~ sim.80, data = speciation.ms.p, lambda =0.001)
model.81 <- pgls(E.5 ~ sim.81, data = speciation.ms.p, lambda =0.001)
model.82 <- pgls(E.5 ~ sim.82, data = speciation.ms.p, lambda =0.001)
model.83 <- pgls(E.5 ~ sim.83, data = speciation.ms.p, lambda =0.001)
model.84 <- pgls(E.5 ~ sim.84, data = speciation.ms.p, lambda =0.001)
model.85 <- pgls(E.5 ~ sim.85, data = speciation.ms.p, lambda =0.001)
model.86 <- pgls(E.5 ~ sim.86, data = speciation.ms.p, lambda =0.001)
model.87 <- pgls(E.5 ~ sim.87, data = speciation.ms.p, lambda =0.001)
model.88 <- pgls(E.5 ~ sim.88, data = speciation.ms.p, lambda =0.001)
model.89 <- pgls(E.5 ~ sim.89, data = speciation.ms.p, lambda =0.001)
model.90 <- pgls(E.5 ~ sim.90, data = speciation.ms.p, lambda =0.001)
model.91 <- pgls(E.5 ~ sim.91, data = speciation.ms.p, lambda =0.001)
model.92 <- pgls(E.5 ~ sim.92, data = speciation.ms.p, lambda =0.001)
model.93 <- pgls(E.5 ~ sim.93, data = speciation.ms.p, lambda =0.001)
model.94 <- pgls(E.5 ~ sim.94, data = speciation.ms.p, lambda =0.001)
model.95 <- pgls(E.5 ~ sim.95, data = speciation.ms.p, lambda =0.001)
model.96 <- pgls(E.5 ~ sim.96, data = speciation.ms.p, lambda =0.001)
model.97 <- pgls(E.5 ~ sim.97, data = speciation.ms.p, lambda =0.001)
model.98 <- pgls(E.5 ~ sim.98, data = speciation.ms.p, lambda =0.001)
model.99 <- pgls(E.5 ~ sim.99, data = speciation.ms.p, lambda =0.001)
model.100 <- pgls(E.5 ~ sim.100, data = speciation.ms.p, lambda =0.001)

p.p <- c(summary(model.1)$coefficients[2,4],summary(model.2)$coefficients[2,4],summary(model.3)$coefficients[2,4],summary(model.4)$coefficients[2,4],summary(model.5)$coefficients[2,4],summary(model.6)$coefficients[2,4],summary(model.7)$coefficients[2,4],summary(model.8)$coefficients[2,4],summary(model.9)$coefficients[2,4],summary(model.10)$coefficients[2,4],summary(model.11)$coefficients[2,4],summary(model.12)$coefficients[2,4],summary(model.13)$coefficients[2,4],summary(model.14)$coefficients[2,4],summary(model.15)$coefficients[2,4],summary(model.16)$coefficients[2,4],summary(model.17)$coefficients[2,4],summary(model.18)$coefficients[2,4],summary(model.19)$coefficients[2,4],summary(model.20)$coefficients[2,4],summary(model.21)$coefficients[2,4],summary(model.22)$coefficients[2,4],summary(model.23)$coefficients[2,4],summary(model.24)$coefficients[2,4],summary(model.25)$coefficients[2,4],summary(model.26)$coefficients[2,4],summary(model.27)$coefficients[2,4],summary(model.28)$coefficients[2,4],summary(model.29)$coefficients[2,4],summary(model.30)$coefficients[2,4],summary(model.31)$coefficients[2,4],summary(model.32)$coefficients[2,4],summary(model.33)$coefficients[2,4],summary(model.34)$coefficients[2,4],summary(model.35)$coefficients[2,4],summary(model.36)$coefficients[2,4],summary(model.37)$coefficients[2,4],summary(model.38)$coefficients[2,4],summary(model.39)$coefficients[2,4],summary(model.40)$coefficients[2,4],summary(model.41)$coefficients[2,4],summary(model.42)$coefficients[2,4],summary(model.43)$coefficients[2,4],summary(model.44)$coefficients[2,4],summary(model.45)$coefficients[2,4],summary(model.46)$coefficients[2,4],summary(model.47)$coefficients[2,4],summary(model.48)$coefficients[2,4],summary(model.49)$coefficients[2,4],summary(model.50)$coefficients[2,4],summary(model.51)$coefficients[2,4],summary(model.52)$coefficients[2,4],summary(model.53)$coefficients[2,4],summary(model.54)$coefficients[2,4],summary(model.55)$coefficients[2,4],summary(model.56)$coefficients[2,4],summary(model.57)$coefficients[2,4],summary(model.58)$coefficients[2,4],summary(model.59)$coefficients[2,4],summary(model.60)$coefficients[2,4],summary(model.61)$coefficients[2,4],summary(model.62)$coefficients[2,4],summary(model.63)$coefficients[2,4],summary(model.64)$coefficients[2,4],summary(model.65)$coefficients[2,4],summary(model.66)$coefficients[2,4],summary(model.67)$coefficients[2,4],summary(model.68)$coefficients[2,4],summary(model.69)$coefficients[2,4],summary(model.70)$coefficients[2,4],summary(model.71)$coefficients[2,4],summary(model.72)$coefficients[2,4],summary(model.73)$coefficients[2,4],summary(model.74)$coefficients[2,4],summary(model.75)$coefficients[2,4],summary(model.76)$coefficients[2,4],summary(model.77)$coefficients[2,4],summary(model.78)$coefficients[2,4],summary(model.79)$coefficients[2,4],summary(model.80)$coefficients[2,4],summary(model.81)$coefficients[2,4],summary(model.82)$coefficients[2,4],summary(model.83)$coefficients[2,4],summary(model.84)$coefficients[2,4],summary(model.85)$coefficients[2,4],summary(model.86)$coefficients[2,4],summary(model.87)$coefficients[2,4],summary(model.88)$coefficients[2,4],summary(model.89)$coefficients[2,4],summary(model.90)$coefficients[2,4],summary(model.91)$coefficients[2,4],summary(model.92)$coefficients[2,4],summary(model.93)$coefficients[2,4],summary(model.94)$coefficients[2,4],summary(model.95)$coefficients[2,4],summary(model.96)$coefficients[2,4],summary(model.97)$coefficients[2,4],summary(model.98)$coefficients[2,4],summary(model.99)$coefficients[2,4],summary(model.100)$coefficients[2,4])
slope.p <- c(summary(model.1)$coefficients[2,1],summary(model.2)$coefficients[2,1],summary(model.3)$coefficients[2,1],summary(model.4)$coefficients[2,1],summary(model.5)$coefficients[2,1],summary(model.6)$coefficients[2,1],summary(model.7)$coefficients[2,1],summary(model.8)$coefficients[2,1],summary(model.9)$coefficients[2,1],summary(model.10)$coefficients[2,1],summary(model.11)$coefficients[2,1],summary(model.12)$coefficients[2,1],summary(model.13)$coefficients[2,1],summary(model.14)$coefficients[2,1],summary(model.15)$coefficients[2,1],summary(model.16)$coefficients[2,1],summary(model.17)$coefficients[2,1],summary(model.18)$coefficients[2,1],summary(model.19)$coefficients[2,1],summary(model.20)$coefficients[2,1],summary(model.21)$coefficients[2,1],summary(model.22)$coefficients[2,1],summary(model.23)$coefficients[2,1],summary(model.24)$coefficients[2,1],summary(model.25)$coefficients[2,1],summary(model.26)$coefficients[2,1],summary(model.27)$coefficients[2,1],summary(model.28)$coefficients[2,1],summary(model.29)$coefficients[2,1],summary(model.30)$coefficients[2,1],summary(model.31)$coefficients[2,1],summary(model.32)$coefficients[2,1],summary(model.33)$coefficients[2,1],summary(model.34)$coefficients[2,1],summary(model.35)$coefficients[2,1],summary(model.36)$coefficients[2,1],summary(model.37)$coefficients[2,1],summary(model.38)$coefficients[2,1],summary(model.39)$coefficients[2,1],summary(model.40)$coefficients[2,1],summary(model.41)$coefficients[2,1],summary(model.42)$coefficients[2,1],summary(model.43)$coefficients[2,1],summary(model.44)$coefficients[2,1],summary(model.45)$coefficients[2,1],summary(model.46)$coefficients[2,1],summary(model.47)$coefficients[2,1],summary(model.48)$coefficients[2,1],summary(model.49)$coefficients[2,1],summary(model.50)$coefficients[2,1],summary(model.51)$coefficients[2,1],summary(model.52)$coefficients[2,1],summary(model.53)$coefficients[2,1],summary(model.54)$coefficients[2,1],summary(model.55)$coefficients[2,1],summary(model.56)$coefficients[2,1],summary(model.57)$coefficients[2,1],summary(model.58)$coefficients[2,1],summary(model.59)$coefficients[2,1],summary(model.60)$coefficients[2,1],summary(model.61)$coefficients[2,1],summary(model.62)$coefficients[2,1],summary(model.63)$coefficients[2,1],summary(model.64)$coefficients[2,1],summary(model.65)$coefficients[2,1],summary(model.66)$coefficients[2,1],summary(model.67)$coefficients[2,1],summary(model.68)$coefficients[2,1],summary(model.69)$coefficients[2,1],summary(model.70)$coefficients[2,1],summary(model.71)$coefficients[2,1],summary(model.72)$coefficients[2,1],summary(model.73)$coefficients[2,1],summary(model.74)$coefficients[2,1],summary(model.75)$coefficients[2,1],summary(model.76)$coefficients[2,1],summary(model.77)$coefficients[2,1],summary(model.78)$coefficients[2,1],summary(model.79)$coefficients[2,1],summary(model.80)$coefficients[2,1],summary(model.81)$coefficients[2,1],summary(model.82)$coefficients[2,1],summary(model.83)$coefficients[2,1],summary(model.84)$coefficients[2,1],summary(model.85)$coefficients[2,1],summary(model.86)$coefficients[2,1],summary(model.87)$coefficients[2,1],summary(model.88)$coefficients[2,1],summary(model.89)$coefficients[2,1],summary(model.90)$coefficients[2,1],summary(model.91)$coefficients[2,1],summary(model.92)$coefficients[2,1],summary(model.93)$coefficients[2,1],summary(model.94)$coefficients[2,1],summary(model.95)$coefficients[2,1],summary(model.96)$coefficients[2,1],summary(model.97)$coefficients[2,1],summary(model.98)$coefficients[2,1],summary(model.99)$coefficients[2,1],summary(model.100)$coefficients[2,1])
replicate.p <- 1:100
overview.p <- data.frame(cbind(replicate.p, slope.p, p.p))
sig.p <- subset(overview.p, p.p < 0.05) 
nrow(sig.p)






