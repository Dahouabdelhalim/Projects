#This script will conduct a RDA analysis of phylogeo data to test whether genetic variation is partitioned because of 
#climatic differences across a spp distribution (IBE), because of geographic distance (IBD), or a putative barrier.
#You will need to make sure you have the WorldClim data in a folder called 'cur_climate', your genetic data as a fasta file
#and your locality data in a csv file (see the getula examples). 
#works with R version 3.2.1


#Load necessary libraries
library(raster)
library(rworldmap)
library(rgdal)

#Input your spp locality data
DataSpecies=read.csv(file="lampropeltis_points.csv")
#Plot the extent of your study region
plot(getMap(), xlim = c(-126, -90), ylim = c(22, 50), asp=1)
#Plot your species localities
points(DataSpecies$long,DataSpecies$lat,pch=20,col="red")

#Load the WorldClim data (or whatever climatic data you have)
bio1=raster("./cur_climate/bio1.bil")
bio2=raster("./cur_climate/bio2.bil") 
bio3=raster("./cur_climate/bio3.bil") 
bio4=raster("./cur_climate/bio4.bil")
bio5=raster("./cur_climate/bio5.bil")
bio6=raster("./cur_climate/bio6.bil")
bio7=raster("./cur_climate/bio7.bil")
bio8=raster("./cur_climate/bio8.bil") 
bio9=raster("./cur_climate/bio9.bil")
bio10=raster("./cur_climate/bio10.bil")
bio11=raster("./cur_climate/bio11.bil")
bio12=raster("./cur_climate/bio12.bil")
bio13=raster("./cur_climate/bio13.bil")
bio14=raster("./cur_climate/bio14.bil")
bio15=raster("./cur_climate/bio15.bil")
bio16=raster("./cur_climate/bio16.bil")
bio17=raster("./cur_climate/bio17.bil")
bio18=raster("./cur_climate/bio18.bil")
bio19=raster("./cur_climate/bio19.bil")

#Define the extent of your study region, change this to wherever you're study region is.
ext_user=extent(c(-126,-90,22,50))
#Crop the climate layers
AnMeanTemp=crop(bio1,ext_user)
MeanDRange=crop(bio2,ext_user)
Iso=crop(bio3,ext_user)
TempSeaso=crop(bio4,ext_user)
MaxTempWarmMonth=crop(bio5,ext_user)
MintempColdMonth=crop(bio6,ext_user)
TempeAnRange=crop(bio7,ext_user)
MeanTempWetQuarter=crop(bio8,ext_user)
MeanTempDriQuarter=crop(bio9,ext_user)
MeanTempWarmQuarter=crop(bio10,ext_user)
MeanTempColdQuarter=crop(bio11,ext_user)
AnPreci=crop(bio12,ext_user)
PreciWetMonth=crop(bio13,ext_user)
PreciDriMonth=crop(bio14,ext_user)
PreciSeaso=crop(bio15,ext_user)
PreciWetQuarter=crop(bio16,ext_user)
PreciDriQuarter=crop(bio17,ext_user)
PreciWarmQuarter=crop(bio18,ext_user)
PreciColdQuarter=crop(bio19,ext_user)

#Make sure you did it right
plot(PreciSeaso)
#Add your taxon data if you want
points(DataSpecies$long,DataSpecies$lat,pch=20,col="red")

#Remove unnecessary files from your workspace
rm(bio1, bio2, bio3, bio4, bio5, bio6, bio7, bio8, bio9, bio10, bio11, bio12, bio13, bio14, bio15, bio17, bio18, bio19)

#Extract the clim data from your locality data
newbio1 <- subset((extract(AnMeanTemp, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio1)
newbio2 <- subset((extract(MeanDRange, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio2)
newbio3 <- subset((extract(Iso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio3)
newbio4 <- subset((extract(TempSeaso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio4)
newbio5 <- subset((extract(MaxTempWarmMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio5)
newbio6 <- subset((extract(MintempColdMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio6)
newbio7 <- subset((extract(TempeAnRange, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio7)
newbio8 <- subset((extract(MeanTempWetQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio8)
newbio9 <- subset((extract(MeanTempDriQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio9)
newbio10 <- subset((extract(MeanTempWarmQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio10)
newbio11 <- subset((extract(MeanTempColdQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio11)
newbio12 <- subset((extract(AnPreci, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio12)
newbio13 <- subset((extract(PreciWetMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio13)
newbio14 <- subset((extract(PreciDriMonth, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio14)
newbio15 <- subset((extract(PreciSeaso, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio15)
newbio16 <- subset((extract(PreciWetQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio16)
newbio17 <- subset((extract(PreciDriQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio17)
newbio18 <- subset((extract(PreciWarmQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio18)
newbio19 <- subset((extract(PreciColdQuarter, DataSpecies[,2:3], method='simple', buffer=5000, fun=mean, df=TRUE)), select=bio19)
#View the first few rows of the vector you just made
head(newbio1)

#Remove correlated variables
rm(newbio3, newbio7)

#Bind these together with your original locality file
lampropeltis_clim_data_extracted<-cbind(DataSpecies, newbio1, newbio2, newbio4, newbio5, newbio6, newbio8, newbio9, newbio10, newbio11, newbio12, newbio13, newbio14, newbio15, newbio16, newbio17, newbio18, newbio19)
head(lampropeltis_clim_data_extracted)

#Write your finished table to save these data!
write.table(lampropeltis_clim_data_extracted, file= "lampropeltis_clim_data_extracted.txt")



#Now we can do the RDA analysis
#First call 'ape', 'vegan', and get your genetic data
library(ape)
library(vegan)
#Load DNA
dna <- read.dna(file = "lampropeltis.fasta", format = "fasta")
#Generate a genetic distance matrix, may want to use a different model of substitution
lampropeltisDist <- dist.dna(dna, model = "GG95")
#PCoA for gen dist matrix
pcoa((as.dist(scale(lampropeltisDist, F))), correction = "none", rn= NULL)->PCOA_lampropeltis

#If you need to recall your clim and locality data:
#lampropeltis_clim_data_extracted<-read.table("getula_clim_data_extracted.txt",na.string="NA", row.names=1, header=T)



#Run Full Model with all variables
rdaFullModel <- rda(PCOA_lampropeltis$vectors~long+lat+bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+locale, data=lampropeltis_clim_data_extracted)
#Summarize this run
#summary(rdaFullModel)
#Get Rsquared ajusted values
RsquareAdj(rdaFullModel)
#Conduct the rda anova
anova(rdaFullModel)
#Plot your results
#plot(rdaFullModel)

#Above I've blocked the functions 'summary' and 'plot' bc you really just need the R sq., adjusted R sq., and the significance from the anova

#Condition out Geography and Locale so that you are only looking at climate data
rdaClim <- rda(PCOA_lampropeltis$vectors~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19  + Condition(lat + long + locale), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaClim)
#Get Rsquared ajusted values
RsquareAdj(rdaClim)
#Conduct the rda anova
anova(rdaClim)
#Plot your results
#plot(rdaClim)


#Condition out Climate Variables and Locale so you are only looking at the distance between samples
rdaDistance <- rda(PCOA_lampropeltis$vectors~ lat + long + Condition(bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19+locale), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaDistance)
#Get Rsquared ajusted values
RsquareAdj(rdaDistance)
#Conduct the rda anova
anova(rdaDistance)
#Plot your results
#plot(rdaDistance)

#Condition out Geography and Climate Variables so you are only looking at locality with respect to phylogeo barrier
rdaLocale <- rda(PCOA_lampropeltis$vectors~ locale + Condition(lat+long+bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaLocale)
#Get Rsquared ajusted values
RsquareAdj(rdaLocale)
#Conduct the rda anova
anova(rdaLocale)
#Plot your results
#plot(rdaLocale)

#Condition out Climate so that you are only looking at both Geography and Locale
rdaGeoLoc <- rda(PCOA_lampropeltis$vectors~lat + long + locale + Condition(bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaGeoLoc)
#Get Rsquared ajusted values
RsquareAdj(rdaGeoLoc)
#Conduct the rda anova
anova(rdaGeoLoc)
#Plot your results
#plot(rdaGeoLoc)

#Condition out Locale so that you are only looking at both Geography and Climate
rdaGeoClim <- rda(PCOA_lampropeltis$vectors~lat + long + bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19 + Condition(locale), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaGeoClim)
#Get Rsquared ajusted values
RsquareAdj(rdaGeoClim)
#Conduct the rda anova
anova(rdaGeoClim)
#Plot your results
#plot(rdaGeoClim)

#Condition out Geography so that you are only looking at both Locale and Climate
rdaLocClim <- rda(PCOA_lampropeltis$vectors~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19 + locale + Condition(lat + long), data=lampropeltis_clim_data_extracted)

#Summarize this run
#summary(rdaLocClim)
#Get Rsquared ajusted values
RsquareAdj(rdaLocClim)
#Conduct the rda anova
anova(rdaLocClim)
#Plot your results
#plot(rdaLocClim)

#To get the % varaince explained
var <- varpart(PCOA_lampropeltis$vectors,~bio1+bio2+bio4+bio5+bio6+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19, ~lat + long, ~locale, data=lampropeltis_clim_data_extracted)

pdf("lampropeltis_varpartsprop.pdf",paper="a4r")
plot(var)
dev.off()
