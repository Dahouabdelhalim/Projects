##First, we need to load the packages to read the phylogeny & run the analyses.
library(ape)
library(nlme)
library(geiger)

##load the data & tree
phy <- read.tree("C_tree.nwk")
dat <- read.csv("Continuous.csv", header = T, row.names = 1)

##Now the data & the tree should coincide in both number of species and names. We can verify this using name.check.
name.check(phy,dat)

##Now that we are sure that the tree and dataset match, we can start to explore the phylogenetic GLS.
##There are two primary packages that can be used to conduct PGLS: ape (with nlme) and caper. 
##There are differences between the two packages in how they work and the information each package and corresponding function returns.
##First, we will explore phylogenetic GLS in ape.
##Phylogenetic GLS is basically a linear model in which the covariance (correlation) structure between species is permitted to match that expected under a Brownian motion process* of evolution on the tree. (*Or other processes.) Consequently, the first step is to define this covariance structure. We do this as follows:
bm<-corBrownian(1, phy)
bm

colnames(dat)

modelo1<-gls(Lablength~Labwidth, data=dat, correlation=bm)
summary(modelo1)

modelo2<-gls(Lablength~Labwidth, data=dat, correlation=corPagel(1,phy))
summary(modelo2)

phy$tip.label
rownames(dat)


Lablength <- data$Lablength
names(Lablength) <- data$Species
class(Lablength)
Lablength <- as.factor(Lablength)
class(Lablength)
Lablength

Labwidth <- data$Labwidth
names(Labwidth) <- data$Species
class(Labwidth)
Labwidth <- as.factor(Labwidth)
class(Labwidth)
Labwidth



library(caper)
data <-read.csv("Continuous.csv",header=TRUE)
comp.data<-comparative.data(phy, data, names.col="Species", vcv.dim=2, warn.dropped=TRUE)
modelo4<-pgls(Lablength~Labwidth, data=comp.data)
summary(modelo4)
