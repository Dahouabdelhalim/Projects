########################################################################################
# This script includes code used to perform phylogenetic generalized least squares 
# regression of body size and ecological variables on Procrustes body shape variables
# using the R package geomorph (Baken et al. 2021; Adams et al. 2021).
########################################################################################

#load necessary libraries
library(geomorph)
library(MASS)
library(ape)
library(phytools)
library(geiger)

#read in landmark coordinates from tps file for icefishes
chan_data<-readland.tps(file="~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/all_chan.tps", specID="ID")

#read in landmark coordinates from tps file for notoperches
trem_data<-readland.tps(file="~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/all_trem.tps", specID="ID")

#read in file which lists the position of sliding semilandmarks
slide<-read.delim(file="~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/noto_sliders.txt")
slide<-as.matrix(slide)

#read in phylogeny
noto.tree<-read.tree("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/notothenioid_timetree.tre")

#load in file which includes data on body size and ecological variables (buoyancy, depth, and first 3 PC axes of diet variation) for icefishes and notoperches
dat<-read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/ice-trem_datafile.csv")


############################################
# First complete analyses for icefishes
############################################

#prune phylogeny to include only icefish species for which we have morphometric data
dimnames(chan_data)[[3]]->chan_species
unique(chan_species)->chan_tax
chan.tree<-keep.tip(noto.tree, chan_tax)

#plot tree and check that it is ultrametric
plot(chan.tree)
is.ultrametric(chan.tree)

#perform Generalized Procrustes analysis to remove effects of size, translation, and rotation on landmark data
gpagen(chan_data, curves=slide)->chann_coords

#next we need to get species-averaged Procrustes shape variables
#(1) first import a "classifiers" file, which is an n x 2 matrix containing information on species
#identity for each specimen included in our dataset
c.classifier<-read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/chan_classifiers.csv")

#(2) next set up an empty array to hold species-averaged coordinate data - for this we need number of species, number of coordinates, and number of dimensions
#create object to hold names of all unique species in our dataset
c.species<-as.character(unique(chan_tax))
#create dummy object to get number of rows and columns for array
c.dim<-chann_coords$coords[,,1]
#create empty array to hold species-averaged data  
c.array<-array(0, dim = c(nrow(c.dim), ncol(c.dim), length(c.species)))
dimnames(c.array) <- list(c(), c(), c(c.species))

#(3) set up for loop to get mean Procrustes shape coordinates for each species
for(i in 1:length(c.species)){
  sp <-c.classifier[c.classifier$Species==c.species[i],]  ##
  if(length(sp$ID)==1){
    rownum <- as.numeric(row.names(sp))
    species.coords <- chann_coords$coords[,,c(rownum)]
    c.array[,,c.species[i]] <- species.coords
  }
  else{
    rownum <- as.numeric(row.names(sp))
    species.coords <- chann_coords$coords[,,c(rownum)]
    mean.species.coords <- mshape(species.coords)
    c.array[,,c.species[i]] <- mean.species.coords
  }
}

#check array of species-averaged variables
print(c.array)

#isolate body size and ecological traits for just icefishes
c.dat<-dat[1:14,]

#create unique variables for each of the traits we'll regress on body shape and omit species with missing data for traits
c.size<-as.matrix(c.dat$MaxSize)
row.names(c.size)<-c.dat$Species
na.omit(c.size)->c.size

c.depth<-as.matrix(c.dat$MeanDepth)
row.names(c.depth)<-c.dat$Species
na.omit(c.depth)->c.depth

c.buoy<-as.matrix(c.dat$Buoyancy)
row.names(c.buoy)<-c.dat$Species
na.omit(c.buoy)->c.buoy

c.dietPC1<-as.matrix(c.dat$DietPC1)
row.names(c.dietPC1)<-c.dat$Species
na.omit(c.dietPC1)->c.dietPC1

c.dietPC2<-as.matrix(c.dat$DietPC2)
row.names(c.dietPC2)<-c.dat$Species
na.omit(c.dietPC2)->c.dietPC2

c.dietPC3<-as.matrix(c.dat$DietPC3)
row.names(c.dietPC3)<-c.dat$Species
na.omit(c.dietPC3)->c.dietPC3

##regress maximum body size against species-averaged Procrustes shape variables
#check that we have same set of species sampled for both traits - we do, so proceed with PGLS
name.check(chan.tree,c.size)

#create geomorph data frame for procD.pgls including species-averaged shape coordinates, phylogeny, and data on maximum body size
shape.size <- geomorph.data.frame(coords = c.array, phy=chan.tree, size = c.size)
#perform pgls regression of size on procrustes shape coordinates
summary(procD.pgls(coords ~ size, phy = phy, data = shape.size))

##regress mean depth against species-averaged Procrustes shape variables
#check that we have same set of species sampled for both traits
name.check(chan.tree,c.depth)
#we do not have same set of species sampled for both traits here, so need to prune phylogeny and shape datasets to match depth dataset
row.names(c.depth)->lab
dep <- c.array[,,lab]
dep_phy<-keep.tip(phy = chan.tree, tip = lab)

shape.depth <- geomorph.data.frame(coords = dep, phy=dep_phy, depth = c.depth)
summary(procD.pgls(coords ~ depth, phy=phy, data=shape.depth))

# regress mean percentage buoyancy against species-averaged Procrustes shape variables
name.check(chan.tree,c.buoy)
row.names(c.buoy)->lab
buo <- c.array[,,lab]
buo_phy<-keep.tip(phy = chan.tree, tip = lab)

shape.buoy <- geomorph.data.frame(coords = buo, phy=buo_phy, buoy = c.buoy)
summary(procD.pgls(coords ~ buoy, phy=phy, data=shape.buoy))

# regress diet PC1 against species-averaged Procrustes shape variables
name.check(chan.tree,c.dietPC1)
lab<-row.names(c.dietPC1)
diet1_phy<-keep.tip(chan.tree, lab)

shape.diet1 <- geomorph.data.frame(coords = c.array, phy=diet1_phy, dietPC1 = c.dietPC1)
summary(procD.pgls(coords ~ dietPC1, phy=phy, data=shape.diet1))

# regress diet PC2 against species-averaged Procrustes shape variables
name.check(chan.tree,c.dietPC2)
lab<-row.names(c.dietPC2)
diet2_phy<-keep.tip(chan.tree, lab)

shape.diet2 <- geomorph.data.frame(coords = c.array, phy=diet2_phy, dietPC2 = c.dietPC2)
summary(procD.pgls(coords ~ dietPC2, phy=phy, data=shape.diet2))

# regress diet PC3 against species-averaged Procrustes shape variables

name.check(chan.tree,c.dietPC3)
lab<-row.names(c.dietPC3)
diet3_phy<-keep.tip(chan.tree, lab)

shape.diet3 <- geomorph.data.frame(coords = c.array, phy=diet3_phy, dietPC3 = c.dietPC3)
summary(procD.pgls(coords ~ dietPC3, phy=phy, data=shape.diet3))

############################################
# Next complete analyses for notoperches
############################################

#prune phylogeny to include only notoperch species for which we have morphometric data
dimnames(trem_data)[[3]]->trem_species
unique(trem_species)->trem_tax
trem.tree<-keep.tip(noto.tree, trem_tax)

#plot tree and check that it is ultrametric
plot(trem.tree)
is.ultrametric(trem.tree)

#perform Generalized Procrustes analysis to remove effects of size, translation, and rotation from landmark data
gpagen(trem_data, curves=slide)->trem_coords

#next we need to get species-averaged Procrustes shape variables
#(1) first import a "classifiers" file, which is an n x 2 matrix containing information on species
#identity for each specimen included in our dataset
t.classifier<-read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/trem_classifiers.csv")

#(2) next set up an empty array to hold species-averaged coordinate data - for this we need number of species, number of coordinates, and number of dimensions
#create object to hold names of all unique species in our dataset
t.species<-as.character(unique(trem_tax))
#create dummy object to get number of rows and columns for array
t.dim<-trem_coords$coords[,,1]
#create empty array to hold species-averaged data  
t.array<-array(0, dim = c(nrow(t.dim), ncol(t.dim), length(t.species)))
dimnames(t.array) <- list(c(), c(), c(t.species))

#(3) set up for loop to get mean Procrustes shape coordinates for each species
for(i in 1:length(t.species)){
  sp <-t.classifier[t.classifier$Species==t.species[i],]  ##
  if(length(sp$ID)==1){
    rownum <- as.numeric(row.names(sp))
    species.coords <- trem_coords$coords[,,c(rownum)]
    t.array[,,t.species[i]] <- species.coords
  }
  else{
    rownum <- as.numeric(row.names(sp))
    species.coords <- trem_coords$coords[,,c(rownum)]
    mean.species.coords <- mshape(species.coords)
    t.array[,,t.species[i]] <- mean.species.coords
  }
}

#check array of species-averaged variables
print(t.array)


#isolate body size and ecological traits for icefishes
t.dat<-dat[15:31,]

t.size<-as.matrix(t.dat$MaxSize)
row.names(t.size)<-t.dat$Tree_label
na.omit(t.size)->t.size

t.depth<-as.matrix(t.dat$MeanDepth)
row.names(t.depth)<-t.dat$Tree_label
na.omit(t.depth)->t.depth

t.buoy<-as.matrix(t.dat$Buoyancy)
row.names(t.buoy)<-t.dat$Tree_label
na.omit(t.buoy)->t.buoy

t.dietPC1<-as.matrix(t.dat$DietPC1)
row.names(t.dietPC1)<-t.dat$Tree_label
na.omit(t.dietPC1)->t.dietPC1

t.dietPC2<-as.matrix(t.dat$DietPC2)
row.names(t.dietPC2)<-t.dat$Tree_label
na.omit(t.dietPC2)->t.dietPC2

t.dietPC3<-as.matrix(t.dat$DietPC3)
row.names(t.dietPC3)<-t.dat$Tree_label
na.omit(t.dietPC3)->t.dietPC3

##regress maximum body size against species-averaged Procrustes shape variables

name.check(trem.tree,t.size)
t.size <- t.size[trem.tree$tip.label,]
names(t.size) == trem.tree$tip.label

shape.size <- geomorph.data.frame(coords = t.array, phy=trem.tree, size = t.size)
summary(procD.pgls(coords ~ size, phy = phy, data = shape.size))

##regress mean depth against species-averaged Procrustes shape variables

name.check(trem.tree,t.depth)
row.names(t.depth)->lab
dep <- t.array[,,lab]
dep_phy<-keep.tip(phy = trem.tree, tip = lab)

shape.depth <- geomorph.data.frame(coords = dep, phy=dep_phy, depth = t.depth)
summary(procD.pgls(coords ~ depth, phy=phy, data=shape.depth))

##regress mean percentage buoyancy against species-averaged Procrustes shape variables

name.check(trem.tree,t.buoy)
row.names(t.buoy)->lab
buo <- t.array[,,lab]
buo_phy<-keep.tip(phy = trem.tree, tip = lab)

shape.buoy <- geomorph.data.frame(coords = buo, phy=buo_phy, buoy = t.buoy)
summary(procD.pgls(coords ~ buoy, phy=phy, data=shape.buoy))

##regress dietPC1 against species-averaged Procrustes shape variables

name.check(trem.tree,t.dietPC1)
diet1_phy<-drop.tip(trem.tree, "2340_Trematomus_tokarevi_2003")
t.dietPC1<-t.dietPC1[diet1_phy$tip.label,]
lab<-names(t.dietPC1)
diet1 <- t.array[,,lab]

shape.diet1 <- geomorph.data.frame(coords = diet1, phy=diet1_phy, dietPC1 = t.dietPC1)
summary(shape.dietPC1<-procD.pgls(coords ~ dietPC1, phy=phy, data=shape.diet1))

##regress dietPC2 against species-averaged Procrustes shape variables

name.check(trem.tree,t.dietPC2)
diet2_phy<-drop.tip(trem.tree, "2340_Trematomus_tokarevi_2003")
t.dietPC2<-t.dietPC2[diet2_phy$tip.label,]
lab<-names(t.dietPC2)
diet2 <- t.array[,,lab]

shape.diet2 <- geomorph.data.frame(coords = diet2, phy=diet2_phy, dietPC2 = t.dietPC2)
summary(procD.pgls(coords ~ dietPC2, phy=phy, data=shape.diet2))

##regress dietPC3 against species-averaged Procrustes shape variables

name.check(trem.tree,t.dietPC3)
diet3_phy<-drop.tip(trem.tree, "2340_Trematomus_tokarevi_2003")
t.dietPC3<-t.dietPC3[diet3_phy$tip.label,]
lab<-names(t.dietPC3)
diet3 <- t.array[,,lab]

shape.diet3 <- geomorph.data.frame(coords = diet3, phy=diet3_phy, dietPC3 = t.dietPC3)
summary(procD.pgls(coords ~ dietPC3, phy=phy, data=shape.diet3))
