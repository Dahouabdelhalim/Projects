#########################################
#          Load Data
#########################################

require(geomorph)
require(ape)

source("Buseretal2017_customfunctions.R") #code to take averages code from https://doi.org/10.5061/dryad.2p4k0

coords <- readland.tps(file = "Morphological Data.TPS", specID = "ID") #TPS file with landmark data

classifier <- read.csv(file = "Classifiers.csv") # Classifier used to take averages

species<-read.csv("species.csv", row.names=1) # list of species names

phy<-read.nexus("tree_Burns et al 2019.nexus")

####################################################
# Seperate procrustes alignments: Two modules
####################################################

Y.coords<-two.d.array(coords)

hd <- Y.coords[,c(1:8,29:42)] #cranial module

bod<- Y.coords[,c(9:28,43:48)] # post-cranial module


hd.array<-arrayspecs(hd,11,2) 
bod.array<-arrayspecs(bod,13,2)

hd.aligned<-gpagen(hd.array) #procrustes alignment
bod.aligned<-gpagen(bod.array) #procrustes alignment


avg.aligned.hd<-my.landmark.species.average(classifier = classifier, aligned.coords = hd.aligned$coords) #calculating average landmark position across all individuals per species

avg.aligned.bod<-my.landmark.species.average(classifier = classifier, aligned.coords = bod.aligned$coords) #calculating average landmark position across all individuals per species


####################################################
# Seperate procrustes alignments: Three modules
####################################################

Y.coords<-two.d.array(coords)

hd<-Y.coords[,c(1:8,29:42)] #cranial module
abd<-Y.coords[,c(9:12,23:28,43:48)] #abdominal module
tl<-Y.coords[,c(13:22)] #caudal module

hd.array<-arrayspecs(hd,11,2)
abd.array<-arrayspecs(abd,8,2)
tl.array<-arrayspecs(tl,5,2)

hd.aligned<-gpagen(hd.array)
abd.aligned<-gpagen(abd.array)
tl.aligned<-gpagen(tl.array)

avg.aligned.hd<-my.landmark.species.average(classifier = classifier, aligned.coords = hd.aligned$coords)

avg.aligned.abd<-my.landmark.species.average(classifier = classifier, aligned.coords = abd.aligned$coords)

avg.aligned.tl<-my.landmark.species.average(classifier = classifier, aligned.coords = tl.aligned$coords)

###############################################################
# Removing species from the Citharinoidei from the alingments
################################################################
cith<- c("13_Citharinus_sp", "13_Distichodus_decemmaculatus", "13_Distichodus_fasciolatus", "13_Hemigrammocharax","13_Ichthyborus_sp", "13_Neolebias_trilineatus")

d<-rownames(species)

intersection<-d %in% cith

f<-which(intersection== TRUE)

e<-c(31,42,43,51,59,77)

avg.aligned.abd.characoid<-avg.aligned.abd[,,-e]

avg.aligned.hd.characoid<-avg.aligned.hd[,,-e]

avg.aligned.tl.characoid<-avg.aligned.tl[,,-e]

avg.aligned.bod.characoid<-avg.aligned.bod[,,-e]


species.cith<-species[-f,]


phy.characoid <- drop.tip(phy, cith)

################################
# Specify module configuration
################################
land.gp2<-c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "A", "A", "A", "A", "A", 
	   "A", "A", "B", "B","B" ) # two module hypothesis

land.gp3 <- c("A", "A", "A", "A", "B", "B", "C", "C", "C", "C", "C", "B", "B", "B", "A", "A", "A", "A", 
            "A", "A", "A", "B", "B", "B") #three module hypothesis

############################################################################################################
# Testing Phylogenetic Modularity and comparing the strength of modularity between landmark configurations
############################################################################################################
pm3<-phylo.modularity(avg.aligned.hd,land.gp3,phy=phy, iter=9999, print.progress = FALSE)
pm3
pm2<-phylo.modularity(avg.aligned,land.gp2, phy=phy, iter=9999, print.progress = FALSE)
pm2


pmodel.Z<-compare.CR(pm2,pm3, CR.null = TRUE) #compare strength of modularity between modular hypotheses

summary(pmodel.Z)

##################
#Characoidei only
##################
pm3.characoid<-phylo.modularity(avg.aligned.characoid,land.gp3,phy=phy.characoid, iter=9999, print.progress = FALSE)
pm3.characoid
pm2.characoid<-phylo.modularity(avg.aligned.characoid,land.gp2, phy=phy.characoid, iter=9999, print.progress = FALSE)
pm2.characoid


pmodel.Z.characoid<-compare.CR(pm2,pm3, CR.null = TRUE) #compare strength of modularity between modular hypotheses

summary(pmodel.Z.characoid)

######################################################################################################################
## Phylogenetic integration tests between all pairwise comparisons of modules 

#We need to perform pairwise PLS analysis to determine whether all of the modules exhibit integrated evolution.
#####################################################################################################################

################################################################### 
#Integration between cranial module and post-cranial body module
###################################################################

pls.hd.abd<-phylo.integration(avg.aligned.hd, avg.aligned.bod, phy=phy, iter=9999, print.progress = FALSE) #all specimens
summary(pls.hd.abd)
plot(pls.hd.abd)


pls.hd.abd.characoid<-phylo.integration(avg.aligned.hd.characoid, avg.aligned.bod.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE) #characoid only

############################################################### 
#Integration between cranial, abdominal, and caudal modules
###############################################################

pls.hd.abd<-phylo.integration(avg.aligned.hd, avg.aligned.abd, phy=phy, iter=9999, print.progress = FALSE)

pls.hd.tl<-phylo.integration(avg.aligned.hd, avg.aligned.tl, phy=phy, iter=9999, print.progress = FALSE)

pls.abd.tl<-phylo.integration(avg.aligned.abd, avg.aligned.tl, phy=phy, iter=9999, print.progress = FALSE)


pls.hd.abd.characoid<-phylo.integration(avg.aligned.hd.characoid, avg.aligned.abd.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)

pls.hd.tl.characoid<-phylo.integration(avg.aligned.hd.characoid, avg.aligned.tl.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)

pls.abd.tl.characoid<-phylo.integration(avg.aligned.abd.characoid, avg.aligned.tl.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)


#########################################################################
# Comparing the strength of integration between the different modules
#########################################################################
compare.pls(pls.abd.tl, pls.hd.abd, pls.hd.tl, two.tailed = TRUE)
compare.pls(pls.abd.tl.characoid, pls.hd.abd.characoid, pls.hd.tl.characoid, two.tailed = TRUE)


