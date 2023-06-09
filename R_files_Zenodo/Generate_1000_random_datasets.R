# This code allows you to create 1000 random datasets from the Dolphin location datasets. Each time it
# will randomly allocate species to a different location

# Libraries
library(plyr)
library(DataCombine)

setwd("C:/Users/Jonathan/Desktop/Ale-Resource Partitioning/HomeRange")
rm(list=ls()) # clear data list
# Get everything loaded
files<-dir(pattern="DolphinLocs",full.names=T)
list(files)

locs <-read.csv(files[1],sep=",",header=T) # Full data set containing all the individuals for each sex in each year
head(locs)

tail(locs)
levels(locs$Species)

  
Replicated <- replicate(1000,sample(locs$Species))

Total <- cbind(locs,Replicated) 
head(locs)
head(Total)

write.csv(Total, "RandomisedSpecies.csv", row.names=F)

dta <-read.csv("RandomisedSpecies.csv",sep=",",header=T)
head(dta)
dta$Species.y <- NULL
#dta$SS1 <- NULL
#dta$SS2 <- NULL
#dta$SS3 <- NULL
#dta$SS4 <- NULL

nX <- ncol(dta)-4

listofdf <- lapply(1:nX, function(x) NULL)

#for (i in 1:nX) {
#  listofdf[[i]] <- data.frame(dta$Track_id, dta$Sex.x, dta$date, dta$x, dta$y, dta[i+5])
#}


for (i in 1:nX) {
  tempdf <- data.frame(dta$Date, dta$Species, dta$Lat, dta$Long, dta[i+4])
  write.csv(tempdf, paste0("HR_Random", i, ".csv"),row.names=F)
}




