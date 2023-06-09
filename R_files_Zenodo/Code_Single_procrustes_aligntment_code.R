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

###############################################################
# Removing species from the Citharinoidei from the alingments
################################################################
cith<- c("13_Citharinus_sp", "13_Distichodus_decemmaculatus", "13_Distichodus_fasciolatus", "13_Hemigrammocharax","13_Ichthyborus_sp", "13_Neolebias_trilineatus")

d<-rownames(species)

intersection<-d %in% cith

f<-which(intersection== TRUE)

e<-c(31,42,43,51,59,77)

avg.aligned.characoid<-avg.aligned[,,-e]

species.cith<-species[-f,]

phy.characoid <- drop.tip(phy, cith)

################################################################################
# Specify module configuration
###################################################################################
land.gp2<-c("A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "A", "A", "A", "A", "A", 
	   "A", "A", "B", "B","B" ) # two module hypothesis

land.gp3 <- c("A", "A", "A", "A", "B", "B", "C", "C", "C", "C", "C", "B", "B", "B", "A", "A", "A", "A", 
            "A", "A", "A", "B", "B", "B") #three module hypothesis

############################################################################################################
# Testing Phylogenetic Modularity and comparing the strength of modularity between landmark configurations
############################################################################################################

pm3<-phylo.modularity(avg.aligned,land.gp3,phy=phy, iter=9999, print.progress = FALSE)
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
#All species
#subset landmarks into modules
Y.avg<-two.d.array(avg.aligned)
hd<-Y.avg[,c(1:8,29:42)]
bod<-Y.avg[,c(9:22,23:28,43:48)]


hd.array<-arrayspecs(hd,11,2)
bod.array<-arrayspecs(bod,13,2)

# Integration tests
pls.hd.bod<-phylo.integration(hd.array, bod.array, phy=phy, iter=9999, print.progress = FALSE)
summary(pls.hd.bod)
plot(pls.hd.bod)

#Characoidei only
#subset landmarks into modules
Y.avg.characoid<-two.d.array(avg.aligned.characoid)
hd.characoid<-Y.avg.characoid[,c(1:8,29:42)]
bod.characoid<-Y.avg.characoid[,c(9:22,23:28,43:48)]

hd.array.characoid<-arrayspecs(hd.characoid,11,2)
bod.array.characoid<-arrayspecs(bod.characoid,13,2)

# Integration tests
pls.hd.bod.characoid<-phylo.integration(hd.array.characoid, bod.array.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.hd.bod.characoid)
plot(pls.hd.bod.characoid)

############################################################## 
#Integration between cranial, abdominal, and caudal modules
################################################################
#All species
#subset landmarks into modules
Y.avg<-two.d.array(avg.aligned)
hd<-Y.avg[,c(1:8,29:42)]
abd<-Y.avg[,c(9:12,23:28,43:48)]
tl<-Y.avg[,c(13:22)]

hd.array<-arrayspecs(hd,11,2)
abd.array<-arrayspecs(abd,8,2)
tl.array<-arrayspecs(tl,5,2)

# Integration tests

pls.hd.abd<-phylo.integration(hd.array, abd.array, phy=phy, iter=9999, print.progress = FALSE)
summary(pls.hd.abd)
plot(pls.hd.abd)

pls.hd.tl<-phylo.integration(hd.array, tl.array, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.hd.tl)
plot(pls.hd.tl)

pls.abd.tl<-phylo.integration(abd.array, tl.array, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.abd.tl)
plot(pls.abd.tl)

#Characoidei only
#subset landmarks into modules
Y.avg.characoid<-two.d.array(avg.aligned.characoid)
hd.characoid<-Y.avg.characoid[,c(1:8,29:42)]
abd.characoid<-Y.avg.characoid[,c(9:12,23:28,43:48)]
tl.characoid<-Y.avg.characoid[,c(13:22)]

hd.array.characoid<-arrayspecs(hd.characoid,11,2)
abd.array.characoid<-arrayspecs(abd.characoid,8,2)
tl.array.characoid<-arrayspecs(tl.characoid,5,2)

# Integration tests
pls.hd.abd.characoid<-phylo.integration(hd.array.characoid, abd.array.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.hd.abd.characoid)
plot(pls.hd.abd.characoid)

pls.hd.tl.characoid<-phylo.integration(hd.array.characoid, tl.array.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.hd.tl.characoid)
plot(pls.hd.tl.characoid)

pls.abd.tl.characoid<-phylo.integration(abd.array.characoid, tl.array.characoid, phy=phy.characoid, iter=9999, print.progress = FALSE)
summary(pls.abd.tl.characoid)
plot(pls.abd.tl.characoid)

#########################################################################
# Comparing the strength of integration between the different modules
#########################################################################
compare.pls(pls.abd.tl, pls.hd.abd, pls.hd.tl, two.tailed = TRUE)
compare.pls(pls.abd.tl.characoid, pls.hd.abd.characoid, pls.hd.tl.characoid, two.tailed = TRUE)

#####################################################################################
# Shape change orthogonal to the axis of greatest integration between each module
#####################################################################################

###########################
# Two module hypothesis
###########################
#All species
# cranial module

X <- model.matrix( ~ two.d.array(hd.array) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(hd.array)[[3]]

D <- two.b.pls(I - Hx, hd.array)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

# Post cranial module

X <- model.matrix( ~ two.d.array(bod.array) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(bod.array)[[3]]

D <- two.b.pls(I - Hx, bod.array)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

#Characoidei Only
# cranial module

X <- model.matrix( ~ two.d.array(hd.array.characoidei) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(hd.array.characoidei)[[3]]

D <- two.b.pls(I - Hx, hd.array.characoidei)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

# Post cranial module

X <- model.matrix( ~ two.d.array(bod.array.characoidei) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(bod.array.characoidei)[[3]]

D <- two.b.pls(I - Hx, bod.array.characoidei)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)


###########################
# Three module hypothesis
###########################
#All species
# cranial module
X <- model.matrix( ~ two.d.array(hd.array) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(hd.array)[[3]]

D <- two.b.pls(I - Hx, hd.array)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

#abdominal module
X <- model.matrix( ~ two.d.array(abd.array) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(abd.array)[[3]]

D <- two.b.pls(I - Hx, abd.array)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

#caudal module
X <- model.matrix( ~ two.d.array(tl.array) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(tl.array)[[3]]

D <- two.b.pls(I - Hx, tl.array)
DD <- plot(D, label = rownames(X)) 
#picknplot.shape(DD)

#Characoidei only
# cranial module
X <- model.matrix( ~ two.d.array(hd.array.characoidei) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(hd.array.characoidei)[[3]]

D <- two.b.pls(I - Hx, hd.array.characoidei)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

#abdominal module
X <- model.matrix( ~ two.d.array(abd.array.characoidei) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(abd.array.characoidei)[[3]]

D <- two.b.pls(I - Hx, abd.array.characoidei)
DD <- plot(D, label = rownames(X)) 
picknplot.shape(DD)

#caudal module
X <- model.matrix( ~ two.d.array(tl.array.characoidei) + 0)
U <- qr.Q(qr(X))
Hx <- tcrossprod(U) # same as X %*% solve(crossprod(X)) %*% t(X)
I <- diag(nrow(X))
rownames(Hx) <- rownames(I) <- rownames(X) <- dimnames(tl.array.characoidei)[[3]]

D <- two.b.pls(I - Hx, tl.array.characoidei)
DD <- plot(D, label = rownames(X)) 
#picknplot.shape(DD)

########################################################
#     Module Eigenvalue analysis performed in geomorph
########################################################
#All species

# Two module hypothesis

mod.eig.two<-module.eigen(avg.aligned, partition.gp=land.gp2,phy=phy, transform=TRUE)
summary(mod.eig.two)
plot(mod.eig.two)

# Three module hypothesis

mod.eig.three<-module.eigen(avg.aligned, partition.gp=land.gp3,phy=phy, transform=TRUE)
summary(mod.eig.three)
plot(mod.eig.three)

#Characoidei only

# Two module hypothesis

mod.eig.two<-module.eigen(avg.aligned.characoid, partition.gp=land.gp2,phy=phy.characoid, transform=TRUE)
summary(mod.eig.two)
plot(mod.eig.two)

# Three module hypothesis

mod.eig.three<-module.eigen(avg.aligned.characoid, partition.gp=land.gp3,phy=phy.characoid, transform=TRUE)
summary(mod.eig.three)
plot(mod.eig.three)

##########################################################
#    Module heatmap performed in geomorph
##########################################################

# All species

# Two module hypothesis
module.map(avg.aligned, partition.gp = land.gp2, phy=phy, transform. = FALSE)

# Three module hypothesis
module.map(avg.aligned, partition.gp = land.gp3, phy=phy, transform. = FALSE)

# Characoidei Only

# Two module hypothesis
module.map(avg.aligned.characoidei, partition.gp = land.gp2, phy=phy, transform. = FALSE)

# Three module hypothesis
module.map(avg.aligned.characoidei, partition.gp = land.gp3, phy=phy, transform. = FALSE)

