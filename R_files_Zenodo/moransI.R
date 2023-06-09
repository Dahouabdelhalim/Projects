setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

library(ape)

sns<-read.csv("snsUR.csv")

#Calculate inverse distance matrix
dist=as.matrix(dist(cbind(sns$northing, sns$easting)))
dist.inv=1/dist
diag(dist.inv)<-0

#Calclulate Moran's I

#PATHOGEN PREVALENCE
## Trypanosome prevalence
Moran.I(sns$true.tryp, dist.inv) # in this case we are hoping for p-values above 0.05
## Neogregarine prevalence
Moran.I(sns$true.neo, dist.inv) # in this case we are hoping for p-values above 0.05
## Nosema prevalence
Moran.I(sns$true.Nc, dist.inv) # in this case we are hoping for p-values above 0.05
## general prevalence
Moran.I(sns$true.prev, dist.inv) # in this case we are hoping for p-values above 0.05
