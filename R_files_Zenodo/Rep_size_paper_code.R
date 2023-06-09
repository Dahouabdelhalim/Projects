### Code for paper on repertoire size and habitat availability ###

setwd("") #folder containing files

#### Habitat metrics ####

## Import ALA records
# Need to download records from ALA as csv
#https://biocache.ala.org.au/occurrences/search?q=lsid:https://biodiversity.org.au/afd/taxa/a225c76e-6f2f-40e5-af45-71c34d44150e#tab_recordsView
#Contact author if any issues
ALA.records <- read.csv("") #Name of downloaded file
View(ALA.records)

#remove points outside known extent (0.5 degrees outside known region)
ALA.adjusted <- subset(ALA.records, Latitude < -27.4154 & Latitude > -29.3937
                       & Longitude < 153.966 & Longitude > 151.8759)

View(ALA.adjusted)

write.csv(ALA.adjusted, "ALA.records.adjusted.csv")

#### Species Distribution Model ####
# Code is a combination of rspatial.org (Hijmans) and Julia Ryeland's code
#https://rspatial.org/raster/sdm/2_sdm_occdata.html

Unfiltered_data <- ALA.adjusted

library(raster)
library(rgdal)
library(dismo)
library(rJava)
#There may be an issue with rJava. This is needed to create the SDM. If there are issues, skip ahead to "Vegetation Suitability"

#plot points
library(maptools)
data(wrld_simpl)
plot(wrld_simpl, xlim=c(150, 155), ylim=c(-30, -26), axes=TRUE, col="light yellow")
box()
points(Unfiltered_data$Longitude, Unfiltered_data$Latitude, col="red", pch=20, cex=0.75)

#turn occurrence dataset to spatial dataframe
library(sp)
coordinates(Unfiltered_data) <- ~Longitude+Latitude
crs(Unfiltered_data) <- crs(wrld_simpl)

#create background points
mask <- raster(Unfiltered_data)
#change resolution of raster
res(mask)
res(mask) <- res(mask)/207.5961 #matches resolution of Bioclim variables
res(mask)
mask
area(mask)

set.seed(123)
bg <- randomPoints(mask, 10000)

#plot results
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)

## Add environmental data
#create extent
ALB.ext <- c(151.8759, 153.966, -29.3937, -27.4154)
extent(ALB.ext)

#get bioclimatic variables
climate <- raster::getData('worldclim', var='bio', res=0.5, lon=151, lat=-27)
crs(climate)
crs(climate)<- "+proj=longlat +datum=WGS84 +ellps=WGS84"
climate.aust <- crop(climate, ALB.ext)
res(climate.aust)
extent(climate.aust)
names(climate.aust) <- c("Bio1","Bio2","Bio3","Bio4","Bio5","Bio6","Bio7","Bio8","Bio9", "Bio10",
                         "Bio11","Bio12","Bio13","Bio14","Bio15","Bio16","Bio17","Bio18","Bio19")
plot(climate.aust, colNA = "blue")
saveRDS(climate.aust, 'climate.aust.rds')

#get values of climate variables for presence and background points
presvals <- extract(climate.aust, Unfiltered_data)
absvals <- extract(climate.aust, bg)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
pb
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
head(sdmdata)
tail(sdmdata)
summary(sdmdata)

pairs(sdmdata[,2:20], cex=0.1)

#Checking colinearlity and standardising predictor layers
library(usdm)
All_var_vif <- vif(sdmdata[,2:20])
All_var_vif
All_var_vifcor <- vifcor(sdmdata[,2:20])
All_var_vifcor

#new set of predictor variables with low colinearity
climate.aust1 <- stack(climate.aust$Bio3, climate.aust$Bio4, climate.aust$Bio6, climate.aust$Bio8, 
                       climate.aust$Bio9, climate.aust$Bio14, climate.aust$Bio15, climate.aust$Bio19)
presvals1 <- extract(climate.aust1, Unfiltered_data)
absvals1 <- extract(climate.aust1, bg)

pb1 <- c(rep(1, nrow(presvals1)), rep(0, nrow(absvals1)))
sdmdata1 <- data.frame(cbind(pb1, rbind(presvals1, absvals1)))
head(sdmdata1)
tail(sdmdata1)
summary(sdmdata1)

pairs(sdmdata1[, 2:9], cex=0.1)

saveRDS(sdmdata1, "sdm.Rds")
saveRDS(presvals1, "pvals.Rds")

## Run model!

#create training and testing dataset
set.seed(1)
group <- kfold(Unfiltered_data, 5)
pres_train <- Unfiltered_data[group !=1,]
pres_test <- Unfiltered_data[group == 1,]

#training and testing background data
set.seed(2)
group <- kfold(bg, 5)
bg_train <- bg[group !=1, ]
bg_test <- bg[group == 1,]

r <- raster(climate.aust1, 1)
plot(!is.na(r), col=c("white", "light grey"), legend = FALSE)
points(bg_train, pch="-", cex=0.5, col="yellow")
points(bg_test, pch="-", cex=0.5, col="black")
points(pres_train, pch="+", col="green")
points(pres_test, pch="+", col="blue")

#Run model
maxent()
xm <- maxent(climate.aust1, pres_train)
plot(xm)
response(xm)
e <- evaluate(pres_test, bg_test, xm, climate.aust1)
e

#plot results
par(mfrow=c(1,2))
density(e)
boxplot(e, col=c('blue', 'red'))
plot(e, 'ROC')
px <- predict(climate.aust1, xm, ext=ALB.ext, progress='')
par(mfrow=c(1,2))
plot(px, main="Maxent, raw values")
tr <- threshold(e, "spec_sens") #0.4750297
tr
tr2 <- threshold(e, "equal_sens_spec") #0.4156999
tr2
plot(px > tr, main = "presence/absence")
points(pres_train, pch="+")

par(mfrow=c(1,1))
plot(px)
plot(px > tr, main = "presence/absence")

writeRaster(px, "SDM_R_ALB.asc", overwrite=TRUE)

##SDM is saved in paper files as "SDM_R_ALB.asc"
#Further processing was undertaken in ArcGIS (see methods)


#### Vegetation suitability ####

library(raster)
library(sp)
library(rgdal)

## exported ALA shapefile with dataframe projection #
#Download data from NVIS and save as tif, with mercator projection to match other spatial data
#http://environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7B3F8AD12F-8300-45EC-A41A-469519A94039%7D
#Contact author if any issues
VEG <- raster("")#name of saved file
crs(VEG)
proj4string(VEG)
proj4string(VEG) <-CRS('+proj=merc +zone=56 +south +datum=WGS84 +units=m +no_defs')
par(mfrow=c(1,1))
plot(VEG)
#now import ALA records that are also in WGS84
#used ArcGIS to convert ALA records to a shapefile
#contact author if any issues
ALA_adj <- readOGR("") #name of ALA shp file
crs(ALA_adj)
ALA_adj
proj4string(ALA_adj)
proj4string(ALA_adj) <-CRS('+proj=merc +zone=56 +south +datum=WGS84 +units=m +no_defs')
plot(ALA_adj, col = "red", pch=20, cex=0.75, add=T)

#import extent
ALB_ext2 <-readOGR("extent_05.shp")
crs(ALB_ext2)
ALB_ext2
proj4string(ALB_ext2)
proj4string(ALB_ext2) <-CRS('+proj=merc +zone=56 +south +datum=WGS84 +units=m +no_defs')
plot(ALB_ext2, col = "yellow", pch=20, cex=0.75, add=T)

#plot all together
plot(VEG)
points(ALA_adj, col = "red", pch=20, cex=0.75)
points(ALB_ext2, col = "yellow", pch=20, cex=0.75)

#Now need to reduce extent further to extent of points only
extent(ALA_adj)
extent(VEG)
extent(ALB_ext2)
ALB.extent2 <- c(16906748, 17139417, -3425851, -3175467)
extent(ALB.extent2)
VEG2<- crop(VEG, ALB.extent2)
extent(VEG2)
plot(VEG2)
points(ALA_adj, col = "red", pch=20, cex=0.75)
points(ALB_ext2, col = "yellow", pch=20, cex=0.75)

## Now need to find number of points in each veg type
class(ALA_adj)
real.veg <- extract(VEG2, ALA_adj, method = 'simple', df=TRUE)
View(real.veg)
#use ddply to get number in each category
library(plyr)
veg.cat <- ddply(real.veg, .(NVIS_MVS3), summarise,
                 count = length(unique(ID)))
View(veg.cat)
write.csv(veg.cat, "veg.cat.csv")


veg.rand <- function(data, n){
  #create random points
  bg <- randomPoints(data, n)#same number of points as real data
  #extract data
  exp.veg <- extract(data, bg, method = 'simple', df=TRUE)
  #summarise by veg type
  exp.veg.cat <- ddply(exp.veg, .(NVIS_MVS3), summarise,
                       count = length(unique(ID)))
  #return results
  exp.veg.cat <- na.omit(exp.veg.cat)
  return(exp.veg.cat = exp.veg.cat)
}

#Get random points based on VEG raster - this ensures same resolution as raster layer
veg.p2 <- rdply(.n = 1000, .expr = veg.rand(VEG2, 3314), .progress = "text") #3314 point to exclude those in cleared areas (see methods)
View(veg.p2)
write.csv(veg.p2, "veg.perm.nocl2.csv")

#Z tests
#1 - Cool temperate rainforest
true.1 <- veg.cat[veg.cat$NVIS_MVS3==1,]
true.1.val <- true.1$count; true.1.val
exp.1 <- veg.p2[veg.p2$NVIS_MVS3==1,]
exp.1.mean <- mean(exp.1$count); exp.1.mean

par(mar = c(2,2,2,2))
hist(exp.1$count, xlim = c(0, 50), col = "black", breaks = 100)
abline(v = true.1$count, col = "blue", lwd = 2)

z.1 <- (true.1.val - exp.1.mean)/sd(exp.1$count); z.1
p.1  <- pnorm(z.1, lower.tail=TRUE); p.1

#2 - Tropical or subtropical rainforest
true.2 <- veg.cat[veg.cat$NVIS_MVS3==2,]
true.2.val <- true.2$count; true.2.val
exp.2 <- veg.p2[veg.p2$NVIS_MVS3==2,]
exp.2.mean <- mean(exp.2$count); exp.2.mean

par(mar = c(2,2,2,2))
hist(exp.2$count, xlim = c(0, 500), col = "black", breaks = 100)
abline(v = true.2$count, col = "blue", lwd = 2)

z.2 <- (true.2.val - exp.2.mean)/sd(exp.2$count); z.2
p.2  <- pnorm(z.2, lower.tail=TRUE); p.2

#3 - Wet sclerophyll
true.3 <- veg.cat[veg.cat$NVIS_MVS3==3,]
true.3.val <- true.3$count; true.3.val
exp.3 <- veg.p2[veg.p2$NVIS_MVS3==3,]
exp.3.mean <- mean(exp.3$count); exp.3.mean

par(mar = c(2,2,2,2))
hist(exp.3$count, xlim = c(100, 650), col = "black", breaks = 100)
abline(v = true.3$count, col = "blue", lwd = 2)

z.3 <- (true.3.val - exp.3.mean)/sd(exp.3$count); z.3
p.3  <- pnorm(z.3, lower.tail=TRUE); p.3

#4 - Open eucalyptus, shrubby understorey
true.4 <- veg.cat[veg.cat$NVIS_MVS3==4,]
true.4.val <- true.4$count; true.4.val
exp.4 <- veg.p2[veg.p2$NVIS_MVS3==4,]
exp.4.mean <- mean(exp.4$count); exp.4.mean

par(mar = c(2,2,2,2))
hist(exp.4$count, xlim = c(0, 400), col = "black", breaks = 100)
abline(v = true.4$count, col = "blue", lwd = 2)

z.4 <- (true.4.val - exp.4.mean)/sd(exp.4$count); z.4
p.4  <- pnorm(z.4, lower.tail=TRUE); p.4

#5 - open eucalyptus, grassy understorey
true.5 <- veg.cat[veg.cat$NVIS_MVS3==5,]
true.5.val <- true.5$count; true.5.val
exp.5 <- veg.p2[veg.p2$NVIS_MVS3==5,]
exp.5.mean <- mean(exp.5$count); exp.5.mean

par(mar = c(2,2,2,2))
hist(exp.5$count, xlim = c(300, 900), col = "black", breaks = 100)
abline(v = true.5$count, col = "blue", lwd = 2)

z.5 <- (true.5.val - exp.5.mean)/sd(exp.5$count); z.5
p.5  <- pnorm(z.5, lower.tail=TRUE); p.5

#6 - warm temperate rainforest
true.6 <- veg.cat[veg.cat$NVIS_MVS3==6,]
true.6.val <- true.6$count; true.6.val
exp.6 <- veg.p2[veg.p2$NVIS_MVS3==6,]
exp.6.mean <- mean(exp.6$count); exp.6.mean

par(mar = c(2,2,2,2))
hist(exp.6$count, xlim = c(0, 1700), col = "black", breaks = 100)
abline(v = true.6$count, col = "blue", lwd = 2)

z.6 <- (true.6.val - exp.6.mean)/sd(exp.6$count); z.6
p.6  <- pnorm(z.6, lower.tail=TRUE); p.6

#9 - Eucalyptus woodlands with tussock grass understorey
true.9 <- veg.cat[veg.cat$NVIS_MVS3==9,]
true.9.val <- true.9$count; true.9.val
exp.9 <- veg.p2[veg.p2$NVIS_MVS3==9,]
exp.9.mean <- mean(exp.9$count); exp.9.mean

par(mar = c(2,2,2,2))
hist(exp.9$count, xlim = c(0, 400), col = "black", breaks = 100)
abline(v = true.9$count, col = "blue", lwd = 2)

z.9 <- (true.9.val - exp.9.mean)/sd(exp.9$count); z.9
p.9  <- pnorm(z.9, lower.tail=TRUE); p.9

#14 - Other acacia forests and woodlands
true.14 <- veg.cat[veg.cat$NVIS_MVS3==14,]
true.14.val <- true.14$count; true.14.val
exp.14 <- veg.p2[veg.p2$NVIS_MVS3==14,]
exp.14.mean <- mean(exp.14$count); exp.14.mean

par(mar = c(2,2,2,2))
hist(exp.14$count, xlim = c(0, 10), col = "black", breaks = 100)
abline(v = true.14$count, col = "blue", lwd = 2)

z.14 <- (true.14.val - exp.14.mean)/sd(exp.14$count); z.14
p.14  <- pnorm(z.14, lower.tail=TRUE); p.14

#15 - Melaleuca open forests and woodlands
true.15 <- veg.cat[veg.cat$NVIS_MVS3==15,]
true.15.val <- true.15$count; true.15.val
exp.15 <- veg.p2[veg.p2$NVIS_MVS3==15,]
exp.15.mean <- mean(exp.15$count); exp.15.mean

par(mar = c(2,2,2,2))
hist(exp.15$count, xlim = c(0, 60), col = "black", breaks = 100)
abline(v = true.15$count, col = "blue", lwd = 2)

z.15 <- (true.15.val - exp.15.mean)/sd(exp.15$count); z.15
p.15  <- pnorm(z.15, lower.tail=TRUE); p.15

#16 - Other forests and woodlands
true.16 <- veg.cat[veg.cat$NVIS_MVS3==16,]
true.16.val <- true.16$count; true.16.val
exp.16 <- veg.p2[veg.p2$NVIS_MVS3==16,]
exp.16.mean <- mean(exp.16$count); exp.16.mean

par(mar = c(2,2,2,2))
hist(exp.16$count, xlim = c(0, 100), col = "black", breaks = 100)
abline(v = true.16$count, col = "blue", lwd = 2)

z.16 <- (true.16.val - exp.16.mean)/sd(exp.16$count); z.16
p.16  <- pnorm(z.16, lower.tail=TRUE); p.16

#42 - Naturally bare, sand, rock, claypan, mudflat
true.42 <- veg.cat[veg.cat$NVIS_MVS3==42,]
true.42.val <- true.42$count; true.42.val
exp.42 <- veg.p2[veg.p2$NVIS_MVS3==42,]
exp.42.mean <- mean(exp.42$count); exp.42.mean

par(mar = c(2,2,2,2))
hist(exp.42$count, xlim = c(0, 50), col = "black", breaks = 100)
abline(v = true.42$count, col = "blue", lwd = 2)

z.42 <- (true.42.val - exp.42.mean)/sd(exp.42$count); z.42
p.42  <- pnorm(z.42, lower.tail=TRUE); p.42

#44 - Freshwater
true.44 <- veg.cat[veg.cat$NVIS_MVS3==44,]
true.44.val <- true.44$count; true.44.val
exp.44 <- veg.p2[veg.p2$NVIS_MVS3==44,]
exp.44.mean <- mean(exp.44$count); exp.44.mean

par(mar = c(2,2,2,2))
hist(exp.44$count, xlim = c(0, 100), col = "black", breaks = 100)
abline(v = true.44$count, col = "blue", lwd = 2)

z.44 <- (true.44.val - exp.44.mean)/sd(exp.44$count); z.44
p.44  <- pnorm(z.44, lower.tail=TRUE); p.44

#49 - Melaleuca shrublands and open shrublands
true.49 <- veg.cat[veg.cat$NVIS_MVS3==49,]
true.49.val <- true.49$count; true.49.val
exp.49 <- veg.p2[veg.p2$NVIS_MVS3==49,]
exp.49.mean <- mean(exp.49$count); exp.49.mean

par(mar = c(2,2,2,2))
hist(exp.49$count, xlim = c(0, 15), col = "black", breaks = 100)
abline(v = true.49$count, col = "blue", lwd = 2)

z.49 <- (true.49.val - exp.49.mean)/sd(exp.49$count); z.49
p.49  <- pnorm(z.49, lower.tail=TRUE); p.49

#54 - Eucalyptus tall open forest with a fine-leaved shrubby understorey
true.54 <- veg.cat[veg.cat$NVIS_MVS3==54,]
true.54.val <- true.54$count; true.54.val
exp.54 <- veg.p2[veg.p2$NVIS_MVS3==54,]
exp.54.mean <- mean(exp.54$count); exp.54.mean

par(mar = c(2,2,2,2))
hist(exp.54$count, xlim = c(0, 100), col = "black", breaks = 100)
abline(v = true.54$count, col = "blue", lwd = 2)

z.54 <- (true.54.val - exp.54.mean)/sd(exp.54$count); z.54
p.54  <- pnorm(z.54, lower.tail=TRUE); p.54

#60 - Eucalyptus tall open forests and open forests with ferns, herbs, sedges, rushes or wet tussock grasses
true.60 <- veg.cat[veg.cat$NVIS_MVS3==60,]
true.60.val <- true.60$count; true.60.val
exp.60 <- veg.p2[veg.p2$NVIS_MVS3==60,]
exp.60.mean <- mean(exp.60$count); exp.60.mean

par(mar = c(2,2,2,2))
hist(exp.60$count, xlim = c(0, 200), col = "black", breaks = 100)
abline(v = true.60$count, col = "blue", lwd = 2)

z.60 <- (true.60.val - exp.60.mean)/sd(exp.60$count); z.60
p.60  <- pnorm(z.60, lower.tail=TRUE); p.60

#62 - Dry rainforest or vine thickets
true.62 <- veg.cat[veg.cat$NVIS_MVS3==62,]
true.62.val <- true.62$count; true.62.val
exp.62 <- veg.p2[veg.p2$NVIS_MVS3==62,]
exp.62.mean <- mean(exp.62$count); exp.62.mean

par(mar = c(2,2,2,2))
hist(exp.62$count, xlim = c(0, 30), col = "black", breaks = 100)
abline(v = true.62$count, col = "blue", lwd = 2)

z.62 <- (true.62.val - exp.62.mean)/sd(exp.62$count); z.62
p.62  <- pnorm(z.62, lower.tail=TRUE); p.62

#63 - Sedgelands, rushes or reeds
true.63 <- veg.cat[veg.cat$NVIS_MVS3==63,]
true.63.val <- true.63$count; true.63.val
exp.63 <- veg.p2[veg.p2$NVIS_MVS3==63,]
exp.63.mean <- mean(exp.63$count); exp.63.mean

par(mar = c(2,2,2,2))
hist(exp.63$count, xlim = c(0, 30), col = "black", breaks = 100)
abline(v = true.63$count, col = "blue", lwd = 2)

z.63 <- (true.63.val - exp.63.mean)/sd(exp.63$count); z.63
p.63  <- pnorm(z.63, lower.tail=TRUE); p.63

#80 - Other sparse shrublands and sparse heathlands
true.80 <- veg.cat[veg.cat$NVIS_MVS3==80,]
true.80.val <- true.80$count; true.80.val
exp.80 <- veg.p2[veg.p2$NVIS_MVS3==80,]
exp.80.mean <- mean(exp.80$count); exp.80.mean

par(mar = c(2,2,2,2))
hist(exp.80$count, xlim = c(0, 50), col = "black", breaks = 100)
abline(v = true.80$count, col = "blue", lwd = 2)

z.80 <- (true.80.val - exp.80.mean)/sd(exp.80$count); z.80
p.80  <- pnorm(z.80, lower.tail=TRUE); p.80

#SDM and suitable vegetation types were combined in ArcMap.
#Patch size and local habitat availability for each individual were calculated in ArcMap.

# Preparing Repertoire Data ####
#Preparing data for analysis

#Find Sample sizes ####
all.mimicry <- read.csv("all.mimicry.cl.csv") #Main datasheet
View(all.mimicry)
table(all.mimicry$Cat)
table(all.mimicry$Begin.File)
table(all.mimicry$Species)
table(all.mimicry$Type.new)
table(all.mimicry$sp2)

all.mimicry <- all.mimicry[order(all.mimicry$Bird.ID, all.mimicry$Begin.File, all.mimicry$Begin.Time..s.),]
View(all.mimicry)

#Give sequential number to units within birds for sample size selection
all.mimicry$X <- 1:10531
all.mimicry$unit.no <- ave(all.mimicry$X,
                           all.mimicry$Bird.ID,
                           FUN = seq_along)
View(all.mimicry)

#get full number for each bird
library(plyr)
samp.size <- ddply(all.mimicry, .(Bird.ID, Population), summarise,
                   N = length(Bird.ID))
View(samp.size)

#remove unknown vocalisations
all.mimicry.v <- all.mimicry[!all.mimicry$Type.new=="UN",]

#Get total possible sample size
rep.size.all <- ddply(all.mimicry.v, .(Population, Bird.ID), summarise,
                      m.rep = length(unique(Type.new)))
View(rep.size.all)

samp.1 <- merge(samp.size, rep.size.all, by = c("Bird.ID", "Population"))
View(samp.1)

#now get different samples for all unit types

#strategy - get smallest sample size, removing smallest sample one at a time
mim.108 <- all.mimicry.v[all.mimicry.v$unit.no <= 108,]
mim.115 <- all.mimicry.v[all.mimicry.v$unit.no <= 115,]
mim.149 <- all.mimicry.v[all.mimicry.v$unit.no <= 149,]
mim.175 <- all.mimicry.v[all.mimicry.v$unit.no <= 175,]
mim.182 <- all.mimicry.v[all.mimicry.v$unit.no <= 182,]
mim.196 <- all.mimicry.v[all.mimicry.v$unit.no <= 196,]
mim.226 <- all.mimicry.v[all.mimicry.v$unit.no <= 226,] #this level keeps 80% of data


rep.size.108 <- ddply(mim.108, .(Population, Bird.ID), summarise,
                      m.rep.108 = length(unique(Type.new)))
View(rep.size.108)

rep.size.115 <- ddply(mim.115, .(Population, Bird.ID), summarise,
                      m.rep.115 = length(unique(Type.new)))
View(rep.size.115)
rep.size.115 <- rep.size.115[!rep.size.115$Bird.ID == "KLRBG",]

rep.size.149 <- ddply(mim.149, .(Population, Bird.ID), summarise,
                      m.rep.149 = length(unique(Type.new)))
View(rep.size.149)
#remove birds that no longer have a large enough sample
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLRBG",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "MWLT3",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLARH",]

rep.size.175 <- ddply(mim.175, .(Population, Bird.ID), summarise,
                      m.rep.175 = length(unique(Type.new)))
View(rep.size.175)
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLRBG",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "MWLT3",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLARH",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLLBG",]

rep.size.182 <- ddply(mim.182, .(Population, Bird.ID), summarise,
                      m.rep.182 = length(unique(Type.new)))
View(rep.size.182)
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLRBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "MWLT3",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLARH",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLLBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "CGMC4",]

rep.size.196 <- ddply(mim.196, .(Population, Bird.ID), summarise,
                      m.rep.196 = length(unique(Type.new)))
View(rep.size.196)
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLRBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "MWLT3",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLARH",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLLBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "CGMC4",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLCB",]

rep.size.226 <- ddply(mim.226, .(Population, Bird.ID), summarise,
                      m.rep.226 = length(unique(Type.new)))
View(rep.size.226)
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLRBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "MWLT3",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLARH",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLLBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "CGMC4",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLCB",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "TMCGB",]

#merge all together
rep.test <- merge(samp.1, rep.size.108, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.115, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.149, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.175, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.182, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.196, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.226, by = c("Population", "Bird.ID"), all=TRUE)
View(rep.test)


#Test for full vocal repertoire

#n = 108, 100% of birds
plot(m.rep.108 ~ m.rep, rep.test)
summary(lm(m.rep.108 ~ m.rep, rep.test))
rep.test$diff.108 <- 1 - (rep.test$m.rep - rep.test$m.rep.108)/rep.test$m.rep
View(rep.test)
mean(rep.test$diff.108) 
min(rep.test$diff.108)
plot(diff.108 ~ N, rep.test)
summary(lm(diff.108 ~ N, rep.test))

#n = 115, 97% of birds
plot(m.rep.115 ~ m.rep, rep.test)
summary(lm(m.rep.115 ~ m.rep, rep.test))
rep.test$diff.115 <- 1 - (rep.test$m.rep - rep.test$m.rep.115)/rep.test$m.rep
mean(rep.test$diff.115, na.rm=TRUE) #89.7%
min(rep.test$diff.115, na.rm=TRUE)
plot(diff.115 ~ N, rep.test)
summary(lm(diff.115 ~ N, rep.test))

#n = 149, 91% of birds
plot(m.rep.149 ~ m.rep, rep.test)
summary(lm(m.rep.149 ~ m.rep, rep.test)) #R-squared 0.8814
rep.test$diff.149 <- 1 - (rep.test$m.rep - rep.test$m.rep.149)/rep.test$m.rep
mean(rep.test$diff.149, na.rm=TRUE) #91%
min(rep.test$diff.149, na.rm=TRUE) #76.9%
plot(diff.149 ~ N, rep.test)
summary(lm(diff.149 ~ N, rep.test))
#This gets close to full repertoire for most birds without leaving out too many birds

#n = 175, 88.6% of birds
plot(m.rep.175 ~ m.rep, rep.test)
summary(lm(m.rep.175 ~ m.rep, rep.test))
rep.test$diff.175 <- 1 - (rep.test$m.rep - rep.test$m.rep.175)/rep.test$m.rep
mean(rep.test$diff.175, na.rm=TRUE) #92.5%
min(rep.test$diff.175, na.rm=TRUE) #76.9%

#n = 182, 85.7% of birds
plot(m.rep.182 ~ m.rep, rep.test)
summary(lm(m.rep.182 ~ m.rep, rep.test))
rep.test$diff.182 <- 1 - (rep.test$m.rep - rep.test$m.rep.182)/rep.test$m.rep
mean(rep.test$diff.182, na.rm=TRUE) #92.8%
min(rep.test$diff.182, na.rm=TRUE) #76.9%

#n = 196, 82.9% of birds
plot(m.rep.196 ~ m.rep, rep.test)
summary(lm(m.rep.196 ~ m.rep, rep.test))
rep.test$diff.196 <- 1 - (rep.test$m.rep - rep.test$m.rep.196)/rep.test$m.rep
mean(rep.test$diff.196, na.rm=TRUE) #93.6%
min(rep.test$diff.196, na.rm=TRUE) #80.5%

#n = 226, 80% of birds
plot(m.rep.226 ~ m.rep, rep.test)
summary(lm(m.rep.226 ~ m.rep, rep.test))
rep.test$diff.226 <- 1 - (rep.test$m.rep - rep.test$m.rep.226)/rep.test$m.rep
mean(rep.test$diff.226, na.rm=TRUE) #95.7%
min(rep.test$diff.226, na.rm=TRUE) #83.3%


### Vocal units
voc.mim <- all.mimicry.v[all.mimicry.v$Cat=="Vocal",]

voc.mim.rep <- ddply(voc.mim, .(Population, Bird.ID), summarise,
                     voc.rep = length(unique(Type.new)))
View(voc.mim.rep)

samp.2 <- merge(samp.size, voc.mim.rep, by = c("Bird.ID", "Population"))
View(samp.2)

mim.108 <- voc.mim[voc.mim$unit.no <= 108,]
mim.115 <- voc.mim[voc.mim$unit.no <= 115,]
mim.149 <- voc.mim[voc.mim$unit.no <= 149,]
mim.175 <- voc.mim[voc.mim$unit.no <= 175,]
mim.182 <- voc.mim[voc.mim$unit.no <= 182,]
mim.196 <- voc.mim[voc.mim$unit.no <= 196,]
mim.226 <- voc.mim[voc.mim$unit.no <= 226,]

rep.size.108 <- ddply(mim.108, .(Population, Bird.ID), summarise,
                      m.rep.108 = length(unique(Type.new)))
View(rep.size.108)

rep.size.115 <- ddply(mim.115, .(Population, Bird.ID), summarise,
                      m.rep.115 = length(unique(Type.new)))
View(rep.size.115)
rep.size.115 <- rep.size.115[!rep.size.115$Bird.ID == "KLRBG",]

rep.size.149 <- ddply(mim.149, .(Population, Bird.ID), summarise,
                      m.rep.149 = length(unique(Type.new)))
View(rep.size.149)
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLRBG",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "MWLT3",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLARH",]

rep.size.175 <- ddply(mim.175, .(Population, Bird.ID), summarise,
                      m.rep.175 = length(unique(Type.new)))
View(rep.size.175)
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLRBG",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "MWLT3",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLARH",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLLBG",]

rep.size.182 <- ddply(mim.182, .(Population, Bird.ID), summarise,
                      m.rep.182 = length(unique(Type.new)))
View(rep.size.182)
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLRBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "MWLT3",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLARH",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLLBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "CGMC4",]

rep.size.196 <- ddply(mim.196, .(Population, Bird.ID), summarise,
                      m.rep.196 = length(unique(Type.new)))
View(rep.size.196)
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLRBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "MWLT3",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLARH",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLLBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "CGMC4",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLCB",]

rep.size.226 <- ddply(mim.226, .(Population, Bird.ID), summarise,
                      m.rep.226 = length(unique(Type.new)))
View(rep.size.226)
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLRBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "MWLT3",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLARH",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLLBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "CGMC4",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLCB",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "TMCGB",]


rep.test <- merge(samp.2, rep.size.108, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.115, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.149, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.175, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.182, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.196, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.226, by = c("Population", "Bird.ID"), all=TRUE)
View(rep.test)


#Test for vocal unit repertoires

#n = 108, 100% of birds
plot(m.rep.108 ~ voc.rep, rep.test)
summary(lm(m.rep.108 ~ voc.rep, rep.test))
rep.test$diff.108 <- 1 - (rep.test$voc.rep - rep.test$m.rep.108)/rep.test$voc.rep
mean(rep.test$diff.108) #90.6%
min(rep.test$diff.108) #68%

#n = 115, 97% of birds
plot(m.rep.115 ~ voc.rep, rep.test)
summary(lm(m.rep.115 ~ voc.rep, rep.test))
rep.test$diff.115 <- 1 - (rep.test$voc.rep - rep.test$m.rep.115)/rep.test$voc.rep
mean(rep.test$diff.115, na.rm=TRUE) #90.9%
min(rep.test$diff.115, na.rm=TRUE) #68%

#n = 149, 91% of birds
plot(m.rep.149 ~ voc.rep, rep.test)
summary(lm(m.rep.149 ~ voc.rep, rep.test)) #R-squared 0.8151
rep.test$diff.149 <- 1 - (rep.test$voc.rep - rep.test$m.rep.149)/rep.test$voc.rep
mean(rep.test$diff.149, na.rm=TRUE) #91.9%
min(rep.test$diff.149, na.rm=TRUE) #76%
#Use this level

#n = 175, 88.6% of birds
plot(m.rep.175 ~ voc.rep, rep.test)
summary(lm(m.rep.175 ~ voc.rep, rep.test))
rep.test$diff.175 <- 1 - (rep.test$voc.rep - rep.test$m.rep.175)/rep.test$voc.rep
mean(rep.test$diff.175, na.rm=TRUE) #92.7%
min(rep.test$diff.175, na.rm=TRUE) #76%

#n = 182, 85.7% of birds
plot(m.rep.182 ~ voc.rep, rep.test)
summary(lm(m.rep.182 ~ voc.rep, rep.test))
rep.test$diff.182 <- 1 - (rep.test$voc.rep - rep.test$m.rep.182)/rep.test$voc.rep
mean(rep.test$diff.182, na.rm=TRUE) #92.9%
min(rep.test$diff.182, na.rm=TRUE) #76%

#n = 196, 82.9% of birds
plot(m.rep.196 ~ voc.rep, rep.test)
summary(lm(m.rep.196 ~ voc.rep, rep.test))
rep.test$diff.196 <- 1 - (rep.test$voc.rep - rep.test$m.rep.196)/rep.test$voc.rep
mean(rep.test$diff.196, na.rm=TRUE) #93.5%
min(rep.test$diff.196, na.rm=TRUE) #76%

#n = 226, 80% of birds
plot(m.rep.226 ~ voc.rep, rep.test)
summary(lm(m.rep.226 ~ voc.rep, rep.test))
rep.test$diff.226 <- 1 - (rep.test$voc.rep - rep.test$m.rep.226)/rep.test$voc.rep
mean(rep.test$diff.226, na.rm=TRUE) #95.1%
min(rep.test$diff.226, na.rm=TRUE) #80%

### Non-vocal units
nonvoc.mim <- all.mimicry.v[all.mimicry.v$Cat=="Nonvocal",]

nonvoc.mim.rep <- ddply(nonvoc.mim, .(Population, Bird.ID), summarise,
                        nonvoc.rep = length(unique(Type.new)))
View(nonvoc.mim.rep)

samp.3 <- merge(samp.size, nonvoc.mim.rep, by = c("Bird.ID", "Population"))
View(samp.3)

#Not altering sample size for non-vocal units as not analysing separately

### Species repertoire
sp.mim <- voc.mim[!voc.mim$Species=="UN",]
View(sp.mim)

sp.rep <- ddply(sp.mim, .(Population, Bird.ID), summarise,
                sp.rep = length(unique(sp2)))
View(sp.rep)

samp.4 <- merge(samp.size, sp.rep, by = c("Bird.ID", "Population"))
View(samp.4)

mim.108 <- sp.mim[sp.mim$unit.no <= 108,]
mim.115 <- sp.mim[sp.mim$unit.no <= 115,]
mim.149 <- sp.mim[sp.mim$unit.no <= 149,]
mim.175 <- sp.mim[sp.mim$unit.no <= 175,]
mim.182 <- sp.mim[sp.mim$unit.no <= 182,]
mim.196 <- sp.mim[sp.mim$unit.no <= 196,]
mim.226 <- sp.mim[sp.mim$unit.no <= 226,]

rep.size.108 <- ddply(mim.108, .(Population, Bird.ID), summarise,
                      m.rep.108 = length(unique(sp2)))
View(rep.size.108)

rep.size.115 <- ddply(mim.115, .(Population, Bird.ID), summarise,
                      m.rep.115 = length(unique(sp2)))
View(rep.size.115)
rep.size.115 <- rep.size.115[!rep.size.115$Bird.ID == "KLRBG",]

rep.size.149 <- ddply(mim.149, .(Population, Bird.ID), summarise,
                      m.rep.149 = length(unique(sp2)))
View(rep.size.149)
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLRBG",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "MWLT3",]
rep.size.149 <- rep.size.149[!rep.size.149$Bird.ID == "KLARH",]

rep.size.175 <- ddply(mim.175, .(Population, Bird.ID), summarise,
                      m.rep.175 = length(unique(sp2)))
View(rep.size.175)
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLRBG",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "MWLT3",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLARH",]
rep.size.175 <- rep.size.175[!rep.size.175$Bird.ID == "KLLBG",]

rep.size.182 <- ddply(mim.182, .(Population, Bird.ID), summarise,
                      m.rep.182 = length(unique(sp2)))
View(rep.size.182)
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLRBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "MWLT3",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLARH",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "KLLBG",]
rep.size.182 <- rep.size.182[!rep.size.182$Bird.ID == "CGMC4",]

rep.size.196 <- ddply(mim.196, .(Population, Bird.ID), summarise,
                      m.rep.196 = length(unique(sp2)))
View(rep.size.196)
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLRBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "MWLT3",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLARH",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLLBG",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "CGMC4",]
rep.size.196 <- rep.size.196[!rep.size.196$Bird.ID == "KLCB",]

rep.size.226 <- ddply(mim.226, .(Population, Bird.ID), summarise,
                      m.rep.226 = length(unique(sp2)))
View(rep.size.226)
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLRBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "MWLT3",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLARH",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLLBG",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "CGMC4",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "KLCB",]
rep.size.226 <- rep.size.226[!rep.size.226$Bird.ID == "TMCGB",]


rep.test <- merge(samp.4, rep.size.108, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.115, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.149, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.175, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.182, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.196, by = c("Population", "Bird.ID"), all=TRUE)
rep.test <- merge(rep.test, rep.size.226, by = c("Population", "Bird.ID"), all=TRUE)
View(rep.test)


#Test for vocal unit repertoires

#n = 108, 100% of birds
plot(m.rep.108 ~ sp.rep, rep.test)
summary(lm(m.rep.108 ~ sp.rep, rep.test)) #R-squared 0.9595
rep.test$diff.108 <- 1 - (rep.test$sp.rep - rep.test$m.rep.108)/rep.test$sp.rep
mean(rep.test$diff.108) #97.1%
min(rep.test$diff.108) #75%

#n = 115, 97% of birds
plot(m.rep.115 ~ sp.rep, rep.test)
summary(lm(m.rep.115 ~ sp.rep, rep.test))
rep.test$diff.115 <- 1 - (rep.test$sp.rep - rep.test$m.rep.115)/rep.test$sp.rep
mean(rep.test$diff.115, na.rm=TRUE) #97.0%
min(rep.test$diff.115, na.rm=TRUE) #75%

#n = 149, 91% of birds
plot(m.rep.149 ~ sp.rep, rep.test)
summary(lm(m.rep.149 ~ sp.rep, rep.test))
rep.test$diff.149 <- 1 - (rep.test$sp.rep - rep.test$m.rep.149)/rep.test$sp.rep
mean(rep.test$diff.149, na.rm=TRUE) #97.1%
min(rep.test$diff.149, na.rm=TRUE) #75%

#n = 175, 88.6% of birds
plot(m.rep.175 ~ sp.rep, rep.test)
summary(lm(m.rep.175 ~ sp.rep, rep.test))
rep.test$diff.175 <- 1 - (rep.test$sp.rep - rep.test$m.rep.175)/rep.test$sp.rep
mean(rep.test$diff.175, na.rm=TRUE) #97.9%
min(rep.test$diff.175, na.rm=TRUE) #75%

#n = 182, 85.7% of birds
plot(m.rep.182 ~ sp.rep, rep.test)
summary(lm(m.rep.182 ~ sp.rep, rep.test))
rep.test$diff.182 <- 1 - (rep.test$sp.rep - rep.test$m.rep.182)/rep.test$sp.rep
mean(rep.test$diff.182, na.rm=TRUE) #97.9%
min(rep.test$diff.182, na.rm=TRUE) #75%

#n = 196, 82.9% of birds
plot(m.rep.196 ~ sp.rep, rep.test)
summary(lm(m.rep.196 ~ sp.rep, rep.test))
rep.test$diff.196 <- 1 - (rep.test$sp.rep - rep.test$m.rep.196)/rep.test$sp.rep
mean(rep.test$diff.196, na.rm=TRUE) #98.6%
min(rep.test$diff.196, na.rm=TRUE) #83.3%

#n = 226, 80% of birds
plot(m.rep.226 ~ sp.rep, rep.test)
summary(lm(m.rep.226 ~ sp.rep, rep.test))
rep.test$diff.226 <- 1 - (rep.test$sp.rep - rep.test$m.rep.226)/rep.test$sp.rep
mean(rep.test$diff.226, na.rm=TRUE) #99%
min(rep.test$diff.226, na.rm=TRUE) #83.3%

#create repertoire accumulation curves ####
#turn into function
require(data.table)
rep.accum <- function(x){
  dt <- as.data.table(unique(x))
  setkey(dt, "unit.no")
  dt[, count := as.numeric(factor(Type.new, levels = unique(Type.new)))]
  setkey(dt, "unit.no", "count")
  dt.out <- dt[J(unique(unit.no)), mult="last"]
  dt.out[, count := cummax(count)]
  result <- as.data.frame(dt.out)
  return(result = result)
}

sp.accum <- function(x){
  dt <- as.data.table(unique(x))
  setkey(dt, "unit.no")
  dt[, count := as.numeric(factor(sp2, levels = unique(sp2)))]
  setkey(dt, "unit.no", "count")
  dt.out <- dt[J(unique(unit.no)), mult="last"]
  dt.out[, count := cummax(count)]
  result <- as.data.frame(dt.out)
  return(result = result)
}


#now with all birds
rep.split <- split(all.mimicry.v, with(all.mimicry.v, Bird.ID), drop = TRUE)
View(rep.split)
rep.split.a <- lapply(rep.split, function(x) rep.accum(x))
View(rep.split.a)
all.rep <- do.call(rbind, rep.split.a)
View(all.rep)
rownames(all.rep) <- NULL

#vocal sounds only
voc.split <- split(voc.mim, with(voc.mim, Bird.ID), drop = TRUE)
voc.split.a <- lapply(voc.split, function(x) rep.accum(x)) #If this step doesn't work, make sure you run line 656
voc.rep <- do.call(rbind, voc.split.a)
rownames(voc.rep) <- NULL
View(voc.rep)

#Species
sp.split <- split(sp.mim, with(sp.mim, Bird.ID), drop = TRUE)
sp.split.a <- lapply(sp.split, function(x) sp.accum(x)) #If this step doesn't work, make sure you run line 803
sp.rep <- do.call(rbind, sp.split.a)
rownames(sp.rep) <- NULL
View(sp.rep)


#Figure 3 ####
#plot it up
library(ggplot2)

#All vocalisations
ggplot(all.rep, aes(unit.no, count, group = Bird.ID, colour=Population)) + 
  geom_line(size = 0.75) +
  theme_bw() + 
  scale_color_manual(values=c("#D95F02", "#7570B3", "#E7298A", 
                              "#66A61E", "#1B9E77", "#E6AB02", "#A6761D")) +
  geom_vline(xintercept = 149, linetype="dashed", 
             color = "black", size=0.5) + #add line for chosen sample size
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(vjust = 0.5, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14)) +
  xlab("Mimetic units sampled") +
  ylab("Total mimetic units")

#Vocal units
ggplot(voc.rep, aes(unit.no, count, group = Bird.ID, colour=Population)) + 
  geom_line(size = 0.75) +
  theme_bw() +
  geom_vline(xintercept = 149, linetype="dashed", 
             color = "black", size=0.5) +
  scale_color_manual(values=c("#D95F02", "#7570B3", "#E7298A", 
                              "#66A61E", "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(vjust = 0.5, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14)) +
  xlab("Mimetic units sampled") +
  ylab("Vocalisations mimicked")


#Species
ggplot(sp.rep, aes(unit.no, count, group = Bird.ID, colour=Population)) + 
  geom_line(size = 0.75) +
  theme_bw() +
  geom_vline(xintercept = 108, linetype="dashed", 
             color = "black", size=0.5)  +
  scale_color_manual(values=c("#D95F02", "#7570B3", "#E7298A", 
                              "#66A61E", "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(vjust = 0.5, size = 14),
        axis.title.y = element_text(vjust = 1, size = 14)) +
  xlab("Mimetic units sampled") +
  ylab("Species mimicked") +
  scale_y_continuous(name = "Species mimicked", limits = c(0, 11.5))


#Put all data together ####
mimicry.samp <- merge(samp.1, samp.2, by = c("Bird.ID", "Population"))
View(mimicry.samp)
mimicry.samp <- merge(mimicry.samp, samp.3, by = c("Bird.ID", "Population"))
mimicry.samp <- merge(mimicry.samp, samp.4, by = c("Bird.ID", "Population"))
mimicry.samp$N <- mimicry.samp$N.x
mimicry.samp$N.x <- NULL
mimicry.samp$N.y <- NULL
mimicry.samp$N.x <- NULL
mimicry.samp$N.y <- NULL
View(mimicry.samp)

#All vocalisations - 149 units
mim.149 <- all.mimicry.v[all.mimicry.v$unit.no <= 149,]
m.rep.149 <- ddply(mim.149, .(Population, Bird.ID), summarise,
                   m.rep.149 = length(unique(Type.new)))
m.rep.149 <- m.rep.149[!m.rep.149$Bird.ID == "KLRBG",]
m.rep.149 <- m.rep.149[!m.rep.149$Bird.ID == "MWLT3",]
m.rep.149 <- m.rep.149[!m.rep.149$Bird.ID == "KLARH",]

mimicry.samp <- merge(mimicry.samp, m.rep.149, 
                      by = c("Bird.ID", "Population"), all = TRUE)

#Vocal units only
mim.149 <- voc.mim[voc.mim$unit.no <= 149,]
v.rep.149 <- ddply(mim.149, .(Population, Bird.ID), summarise,
                   voc.rep.149 = length(unique(Type.new)))
v.rep.149 <- v.rep.149[!v.rep.149$Bird.ID == "KLRBG",]
v.rep.149 <- v.rep.149[!v.rep.149$Bird.ID == "MWLT3",]
v.rep.149 <- v.rep.149[!v.rep.149$Bird.ID == "KLARH",]

mimicry.samp <- merge(mimicry.samp, v.rep.149, 
                      by = c("Bird.ID", "Population"), all = TRUE)

#Nonvocal units
mim.196 <- nonvoc.mim[nonvoc.mim$unit.no <= 196,]
nv.rep.196 <- ddply(mim.196, .(Population, Bird.ID), summarise,
                    nonvoc.rep.196 = length(unique(Type.new)))
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "KLRBG",]
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "MWLT3",]
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "KLARH",]
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "KLLBG",]
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "CGMC4",]
nv.rep.196 <- nv.rep.196[!nv.rep.196$Bird.ID == "KLCB",]

mimicry.samp <- merge(mimicry.samp, nv.rep.196, 
                      by = c("Bird.ID", "Population"), all = TRUE)

#Species
mim.108 <- sp.mim[sp.mim$unit.no <= 108,]
sp.rep.108 <- ddply(mim.108, .(Population, Bird.ID), summarise,
                    sp.rep.108 = length(unique(sp2)))

mimicry.samp <- merge(mimicry.samp, sp.rep.108, 
                      by = c("Bird.ID", "Population"), all = TRUE)

#Ratios
mimicry.samp$all2sp <- mimicry.samp$m.rep/mimicry.samp$sp.rep
mimicry.samp$voc2sp <- mimicry.samp$voc.rep/mimicry.samp$sp.rep
mimicry.samp$nonvoc2sp <- mimicry.samp$nonvoc.rep/mimicry.samp$sp.rep

#Add spatial data
mimicry.samp$Bird.ID <- as.factor(mimicry.samp$Bird.ID)
mimicry.samp$Population <- as.factor(mimicry.samp$Population)
levels(mimicry.samp$Bird.ID)
levels(mimicry.samp$Population)

spatial.dat <- read.csv("spatial.dat.csv")
#merge
mimicry.sample <- merge(mimicry.samp, spatial.dat, by = c("Bird.ID", "Population"))
View(mimicry.sample)

#Get proportions
mimicry.sample$frag.500 <- (1 - mimicry.sample$out_area_500m/mimicry.sample$buff_area_500m)*100
mimicry.sample$frag.1 <- (1 - mimicry.sample$out_area_1km/mimicry.sample$buff_area_1km)*100
mimicry.sample$frag.2 <- (1 - mimicry.sample$out_area_2km/mimicry.sample$buff_area_2km)*100
mimicry.sample$frag.5 <- (1 - mimicry.sample$out_area_5km/mimicry.sample$buff_area_5km)*100
mimicry.sample$frag.10 <- (1 - mimicry.sample$out_area_10km/mimicry.sample$buff_area_10km)*100

#transform variables
library(bestNormalize)
#NOTE - sometimes different transformations are recommended
#OrderNorm was the most consistently recommended transformation for all variables

(allRbn <- bestNormalize(mimicry.sample$frag.500))
orderNorm_ob <- orderNorm(mimicry.sample$frag.500)
mimicry.sample$frag.500.ON <- predict(orderNorm_ob, newdata=mimicry.sample$frag.500)
hist(mimicry.sample$frag.500)
hist(mimicry.sample$frag.500.ON)

(allRbn <- bestNormalize(mimicry.sample$frag.1))
orderNorm_ob <- orderNorm(mimicry.sample$frag.1)
mimicry.sample$frag.1.ON <- predict(orderNorm_ob, newdata=mimicry.sample$frag.1)
hist(mimicry.sample$frag.1)
hist(mimicry.sample$frag.1.ON)

(allRbn <- bestNormalize(mimicry.sample$frag.2))
orderNorm_ob <- orderNorm(mimicry.sample$frag.2)
mimicry.sample$frag.2.ON <- predict(orderNorm_ob, newdata=mimicry.sample$frag.2)
hist(mimicry.sample$frag.2)
hist(mimicry.sample$frag.2.ON)

(allRbn <- bestNormalize(mimicry.sample$frag.5))
orderNorm_ob <- orderNorm(mimicry.sample$frag.5)
mimicry.sample$frag.5.ON <- predict(orderNorm_ob, newdata=mimicry.sample$frag.5)
hist(mimicry.sample$frag.5)
hist(mimicry.sample$frag.5.ON)

(allRbn <- bestNormalize(mimicry.sample$frag.10))
orderNorm_ob <- orderNorm(mimicry.sample$frag.10)
mimicry.sample$frag.10.ON <- predict(orderNorm_ob, newdata=mimicry.sample$frag.10)
hist(mimicry.sample$frag.10)
hist(mimicry.sample$frag.10.ON)

(allRbn <- bestNormalize(mimicry.sample$patch3))
orderNorm_ob <- orderNorm(mimicry.sample$patch3)
mimicry.sample$patch3.ON <- predict(orderNorm_ob, newdata=mimicry.sample$patch3)
hist(mimicry.sample$patch3)
hist(mimicry.sample$patch3.ON)

### Save final dataset ###
write.csv(mimicry.sample, "full.dataset.csv")

##Get averages
library(plyr)
mim.ave <- ddply(mimicry.sample, .(Population), summarise,
                 patch.size = mean(patch3),
                 frag.500 = mean(frag.500),
                 frag.1 = mean(frag.1),
                 frag.2 = mean(frag.2),
                 frag.5 = mean(frag.5),
                 frag.10 = mean(frag.10))
View(mim.ave)

#Run Models ####

#If jumping to this step, read in full dataset
mimicry.sample <- read.csv("full.dataset.csv")

#Check variable correlations
# include only numeric values here
corr_df <- mimicry.sample[,c(34:39)]
#check structure
str(corr_df)

# make correlation table
corrtab <- cor(corr_df, method = "pearson", use = "complete.obs")
round(corrtab, 2)

#                 patch.size.ON frag.500.ON frag.1.ON frag.2.ON frag.5.ON frag.10.ON patch2.ON patch3.ON
#patch.size.ON          1.00        0.70      0.64      0.52      0.46       0.48      0.56      0.64
#frag.500.ON            0.70        1.00      0.88      0.71      0.35       0.38      0.77      0.80
#frag.1.ON              0.64        0.88      1.00      0.91      0.58       0.54      0.94      0.94
#frag.2.ON              0.52        0.71      0.91      1.00      0.72       0.60      0.91      0.90
#frag.5.ON              0.46        0.35      0.58      0.72      1.00       0.92      0.64      0.65
#frag.10.ON             0.48        0.38      0.54      0.60      0.92       1.00      0.60      0.64
#patch2.ON              0.56        0.77      0.94      0.91      0.64       0.60      1.00      0.97
#patch3.ON              0.64        0.80      0.94      0.90      0.65       0.64      0.97      1.00

#Correlations between population and habitat metrics
aov.patch <- aov(patch3.ON ~ Population, mimicry.sample)
summary(aov.patch)
TukeyHSD(aov.patch)

aov.500m <- aov(frag.500.ON ~ Population, mimicry.sample)
summary(aov.500m)
TukeyHSD(aov.500m)

aov.1km <- aov(frag.1.ON ~ Population, mimicry.sample)
summary(aov.1km)
TukeyHSD(aov.1km)

aov.2km <- aov(frag.2.ON ~ Population, mimicry.sample)
summary(aov.2km)
TukeyHSD(aov.2km)

aov.5km <- aov(frag.5.ON ~ Population, mimicry.sample)
summary(aov.5km)
TukeyHSD(aov.5km)

aov.10km <- aov(frag.10.ON ~ Population, mimicry.sample)
summary(aov.10km)
TukeyHSD(aov.10km)


library(lme4)
library(lmerTest)

#Patch size ####
#Counts
M.patch.all <- glmer(m.rep.149 ~ patch3.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(M.patch.all)

V.patch <- glmer(voc.rep.149 ~ patch3.ON + (1|Population), 
                 data=mimicry.sample, family=poisson())
summary(V.patch)
V.patch.glm <- glm(voc.rep.149 ~ patch3.ON, 
                   data=mimicry.sample, family=poisson())
summary(V.patch.glm)
plot(V.patch.glm)
plot(residuals(V.patch.glm,type = "response")) ## response residuals
plot(residuals(V.patch.glm,type = "pearson")) ## pearson residuals
plot(residuals(V.patch.glm,type = "deviance")) ## deviance residuals
plot(residuals(V.patch.glm) ~ predict(V.patch.glm,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.patch.glm) ~ predict(V.patch.glm,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.patch.glm))
qqline(resid(V.patch.glm))

SP.patch <- glmer(sp.rep.108 ~ patch3.ON + (1|Population), 
                  data=mimicry.sample, family=poisson())
summary(SP.patch)
SP.patch.glm <- glm(sp.rep.108 ~ patch3.ON, 
                    data=mimicry.sample, family=poisson())
summary(SP.patch.glm)
plot(SP.patch.glm)
plot(residuals(SP.patch.glm,type = "response")) ## response residuals
plot(residuals(SP.patch.glm,type = "pearson")) ## pearson residuals
plot(residuals(SP.patch.glm,type = "deviance")) ## deviance residuals
plot(residuals(SP.patch.glm) ~ predict(SP.patch.glm,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.patch.glm) ~ predict(SP.patch.glm,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.patch.glm))
qqline(resid(SP.patch.glm))

#Ratios
all2sp.patch <- lmer(all2sp ~ patch3.ON + (1|Population), 
                     data=mimicry.sample)
summary(all2sp.patch)

v2sp.patch <- lmer(voc2sp ~ patch3.ON + (1|Population), 
                   data=mimicry.sample)
summary(v2sp.patch)

nv2sp.patch <- lmer(nonvoc2sp ~ patch3.ON + (1|Population), 
                    data=mimicry.sample)
summary(nv2sp.patch)


### 500 metres ####
#M.500.full <- glmer(m.rep.149 ~ patch3.ON + frag.500.ON + (1|Population), 
#                   data=mimicry.sample, family=poisson())
#summary(M.500.full)
M.500.frag <- glmer(m.rep.149 ~ frag.500.ON + (1|Population), 
                    data=mimicry.sample, family=poisson())
summary(M.500.frag)
plot(M.500.frag)
plot(residuals(M.500.frag,type = "response")) ## response residuals
plot(residuals(M.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(M.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(M.500.frag) ~ predict(M.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(M.500.frag) ~ predict(M.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(M.500.frag))
qqline(resid(M.500.frag))

#V.500.full <- glmer(voc.rep.149 ~ patch3.ON + frag.500.ON + (1|Population), 
#                    data=mimicry.sample, family=poisson())
#summary(V.500.full)
V.500.frag <- glmer(voc.rep.149 ~ frag.500.ON + (1|Population), 
                    data=mimicry.sample, family=poisson())
summary(V.500.frag)
plot(V.500.frag)
plot(residuals(V.500.frag,type = "response")) ## response residuals
plot(residuals(V.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(V.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(V.500.frag) ~ predict(V.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.500.frag) ~ predict(V.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.500.frag))
qqline(resid(V.500.frag))

#SP.500.full <- glmer(sp.rep.108 ~ patch3.ON + frag.500.ON + (1|Population), 
#                    data=mimicry.sample, family=poisson())
#summary(SP.500.full)
SP.500.frag <- glmer(sp.rep.108 ~ frag.500.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(SP.500.frag)
plot(SP.500.frag)
plot(residuals(SP.500.frag,type = "response")) ## response residuals
plot(residuals(SP.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(SP.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(SP.500.frag) ~ predict(SP.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.500.frag) ~ predict(SP.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.500.frag))
qqline(resid(SP.500.frag))

#AR.500.full <- lmer(all2sp ~ patch3.ON + frag.500.ON + (1|Population), 
#                    data=mimicry.sample)
#summary(AR.500.full)
AR.500.frag <- lmer(all2sp ~ frag.500.ON + (1|Population), 
                    data=mimicry.sample)
summary(AR.500.frag)
plot(AR.500.frag)
plot(residuals(AR.500.frag,type = "response")) ## response residuals
plot(residuals(AR.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(AR.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(AR.500.frag) ~ predict(AR.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(AR.500.frag) ~ predict(AR.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(AR.500.frag))
qqline(resid(AR.500.frag))

#VR.500.full <- lmer(voc2sp ~ patch3.ON + frag.500.ON + (1|Population), 
#                    data=mimicry.sample)
#summary(VR.500.full)
VR.500.frag <- lmer(voc2sp ~ frag.500.ON + (1|Population), 
                    data=mimicry.sample)
summary(VR.500.frag)
plot(VR.500.frag)
plot(residuals(VR.500.frag,type = "response")) ## response residuals
plot(residuals(VR.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(VR.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(VR.500.frag) ~ predict(VR.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(VR.500.frag) ~ predict(VR.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(VR.500.frag))
qqline(resid(VR.500.frag))

#NVR.500.full <- lmer(nonvoc2sp ~ patch.size.ON + frag.500.ON + (1|Population), 
                     #data=mimicry.sample)
#summary(NVR.500.full)
plot(NVR.500.full)
plot(residuals(NVR.500.full,type = "response")) ## response residuals
plot(residuals(NVR.500.full,type = "pearson")) ## pearson residuals
plot(residuals(NVR.500.full,type = "deviance")) ## deviance residuals
plot(residuals(NVR.500.full) ~ predict(NVR.500.full,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(NVR.500.full) ~ predict(NVR.500.full,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(NVR.500.full))
qqline(resid(NVR.500.full))

NVR.500.frag <- lmer(nonvoc2sp ~ frag.500.ON + (1|Population), 
                     data=mimicry.sample)
summary(NVR.500.frag)
plot(NVR.500.frag)
plot(residuals(NVR.500.frag,type = "response")) ## response residuals
plot(residuals(NVR.500.frag,type = "pearson")) ## pearson residuals
plot(residuals(NVR.500.frag,type = "deviance")) ## deviance residuals
plot(residuals(NVR.500.frag) ~ predict(NVR.500.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(NVR.500.frag) ~ predict(NVR.500.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(NVR.500.frag))
qqline(resid(NVR.500.frag))

### 1000m ####
#M.1000.full <- glmer(m.rep.149 ~ frag.1.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample, family=poisson())
#summary(M.1000.full)
M.1000.frag <- glmer(m.rep.149 ~ frag.1.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(M.1000.frag)
plot(M.1000.frag)
plot(residuals(M.1000.frag,type = "response")) ## response residuals
plot(residuals(M.1000.frag,type = "pearson")) ## pearson residuals
plot(residuals(M.1000.frag,type = "deviance")) ## deviance residuals
plot(residuals(M.1000.frag) ~ predict(M.1000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(M.1000.frag) ~ predict(M.1000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(M.1000.frag))
qqline(resid(M.1000.frag))

#V.1000.full <- glmer(voc.rep.149 ~ frag.1.ON + patch3.ON + (1|Population), 
#                    data=mimicry.sample, family=poisson())
#summary(V.1000.full)
V.1000.frag <- glmer(voc.rep.149 ~ frag.1.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(V.1000.frag)
plot(V.1000.frag)
plot(residuals(V.1000.frag,type = "response")) ## response residuals
plot(residuals(V.1000.frag,type = "pearson")) ## pearson residuals
plot(residuals(V.1000.frag,type = "deviance")) ## deviance residuals
plot(residuals(V.1000.frag) ~ predict(V.1000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.1000.frag) ~ predict(V.1000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.1000.frag))
qqline(resid(V.1000.frag))

#SP.1000.full <- glmer(sp.rep.108 ~ frag.1.ON + patch3.ON + (1|Population), 
#                      data=mimicry.sample, family=poisson())
#summary(SP.1000.full)
SP.1000.frag <- glmer(sp.rep.108 ~ frag.1.ON + (1|Population), 
                      data=mimicry.sample, family=poisson())
summary(SP.1000.frag)
plot(SP.1000.frag)
plot(residuals(SP.1000.frag,type = "response")) ## response residuals
plot(residuals(SP.1000.frag,type = "pearson")) ## pearson residuals
plot(residuals(SP.1000.frag,type = "deviance")) ## deviance residuals
plot(residuals(SP.1000.frag) ~ predict(SP.1000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.1000.frag) ~ predict(SP.1000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.1000.frag))
qqline(resid(SP.1000.frag))

#AR.1000.full <- lmer(all2sp ~ frag.1.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample)
#summary(AR.1000.full)
AR.1000.a <- lmer(all2sp ~ frag.1.ON + (1|Population), 
                  data=mimicry.sample)
summary(AR.1000.a)
plot(AR.1000.a)
plot(residuals(AR.1000.a,type = "response")) ## response residuals
plot(residuals(AR.1000.a,type = "pearson")) ## pearson residuals
plot(residuals(AR.1000.a,type = "deviance")) ## deviance residuals
plot(residuals(AR.1000.a) ~ predict(AR.1000.a,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(AR.1000.a) ~ predict(AR.1000.a,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(AR.1000.a))
qqline(resid(AR.1000.a))

#VR.1000.full <- lmer(voc2sp ~ frag.1.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample)
#summary(VR.1000.full)
VR.1000.frag <- lmer(voc2sp ~ frag.1.ON + (1|Population), 
                     data=mimicry.sample)
summary(VR.1000.frag)
plot(VR.1000.frag)
plot(residuals(VR.1000.frag,type = "response")) ## response residuals
plot(residuals(VR.1000.frag,type = "pearson")) ## pearson residuals
plot(residuals(VR.1000.frag,type = "deviance")) ## deviance residuals
plot(residuals(VR.1000.frag) ~ predict(VR.1000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(VR.1000.frag) ~ predict(VR.1000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(VR.1000.frag))
qqline(resid(VR.1000.frag))

#NVR.1000.full <- lmer(nonvoc2sp ~ frag.1.ON + patch3.ON + (1|Population), 
#                      data=mimicry.sample)
#summary(NVR.1000.full)
NVR.1000.frag <- lmer(nonvoc2sp ~ frag.1.ON + (1|Population), 
                      data=mimicry.sample)
summary(NVR.1000.frag)
plot(NVR.1000.frag)
plot(residuals(NVR.1000.frag,type = "response")) ## response residuals
plot(residuals(NVR.1000.frag,type = "pearson")) ## pearson residuals
plot(residuals(NVR.1000.frag,type = "deviance")) ## deviance residuals
plot(residuals(NVR.1000.frag) ~ predict(NVR.1000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(NVR.1000.frag) ~ predict(NVR.1000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(NVR.1000.frag))
qqline(resid(NVR.1000.frag))

### 2000 metres ####
#M.2000.full <- glmer(m.rep.149 ~ frag.2.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample, family=poisson())
#summary(M.2000.full)
M.2000.frag <- glmer(m.rep.149 ~ frag.2.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(M.2000.frag)

#V.2000.full <- glmer(voc.rep.149 ~ frag.2.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample, family=poisson())
#summary(V.2000.full)
V.2000.frag <- glmer(voc.rep.149 ~ frag.2.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(V.2000.frag)
plot(V.2000.frag)
plot(residuals(V.2000.frag,type = "response")) ## response residuals
plot(residuals(V.2000.frag,type = "pearson")) ## pearson residuals
plot(residuals(V.2000.frag,type = "deviance")) ## deviance residuals
plot(residuals(V.2000.frag) ~ predict(V.2000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.2000.frag) ~ predict(V.2000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.2000.frag))
qqline(resid(V.2000.frag))

#SP.2000.full <- glmer(sp.rep.108 ~ frag.2.ON + patch3.ON + (1|Population), 
#                      data=mimicry.sample, family=poisson())
#summary(SP.2000.full)
SP.2000.frag <- glmer(sp.rep.108 ~ frag.2.ON + (1|Population), 
                      data=mimicry.sample, family=poisson())
summary(SP.2000.frag)
plot(SP.2000.frag)
plot(residuals(SP.2000.frag,type = "response")) ## response residuals
plot(residuals(SP.2000.frag,type = "pearson")) ## pearson residuals
plot(residuals(SP.2000.frag,type = "deviance")) ## deviance residuals
plot(residuals(SP.2000.frag) ~ predict(SP.2000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.2000.frag) ~ predict(SP.2000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.2000.frag))
qqline(resid(SP.2000.frag))

#AR.2000.full <- lmer(all2sp ~ frag.2.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample)
#summary(AR.2000.full)
AR.2000.frag <- lmer(all2sp ~ frag.2.ON + (1|Population), 
                     data=mimicry.sample)
summary(AR.2000.frag)
plot(AR.2000.frag)
plot(residuals(AR.2000.frag,type = "response")) ## response residuals
plot(residuals(AR.2000.frag,type = "pearson")) ## pearson residuals
plot(residuals(AR.2000.frag,type = "deviance")) ## deviance residuals
plot(residuals(AR.2000.frag) ~ predict(AR.2000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(AR.2000.frag) ~ predict(AR.2000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(AR.2000.frag))
qqline(resid(AR.2000.frag))

#VR.2000.full <- lmer(voc2sp ~ frag.2.ON + patch3.ON + (1|Population), 
#                     data=mimicry.sample)
#summary(VR.2000.full)
VR.2000.frag <- lmer(voc2sp ~ frag.2.ON + (1|Population), 
                     data=mimicry.sample)
summary(VR.2000.frag)
plot(VR.2000.frag)
plot(residuals(VR.2000.frag,type = "response")) ## response residuals
plot(residuals(VR.2000.frag,type = "pearson")) ## pearson residuals
plot(residuals(VR.2000.frag,type = "deviance")) ## deviance residuals
plot(residuals(VR.2000.frag) ~ predict(VR.2000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(VR.2000.frag) ~ predict(VR.2000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(VR.2000.frag))
qqline(resid(VR.2000.frag))

#NVR.2000.full <- lmer(nonvoc2sp ~ frag.2.ON + patch3.ON + (1|Population), 
#                      data=mimicry.sample)
#summary(NVR.2000.full)
NVR.2000.frag <- lmer(nonvoc2sp ~ frag.2.ON + (1|Population), 
                      data=mimicry.sample)
summary(NVR.2000.frag)

### 5000 metres ####
M.5000.full <- glmer(m.rep.149 ~ frag.5.ON + patch3.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(M.5000.full)
M.5000.frag <- glmer(m.rep.149 ~ frag.5.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(M.5000.frag)

V.5000.full <- glmer(voc.rep.149 ~ frag.5.ON + patch3.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(V.5000.full)
V.5000.full2 <- glm(voc.rep.149 ~ frag.5.ON + patch3.ON, 
                    data=mimicry.sample, family=poisson())
summary(V.5000.full2)
V.5000.frag <- glmer(voc.rep.149 ~ frag.5.ON + (1|Population), 
                     data=mimicry.sample, family=poisson())
summary(V.5000.frag)
plot(V.5000.frag)
plot(residuals(V.5000.frag,type = "response")) ## response residuals
plot(residuals(V.5000.frag,type = "pearson")) ## pearson residuals
plot(residuals(V.5000.frag,type = "deviance")) ## deviance residuals
plot(residuals(V.5000.frag) ~ predict(V.5000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.5000.frag) ~ predict(V.5000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.5000.frag))
qqline(resid(V.5000.frag))

SP.5000.full <- glmer(sp.rep.108 ~ frag.5.ON + patch3.ON + (1|Population), 
                      data=mimicry.sample, family=poisson())
summary(SP.5000.full)
SP.5000.full2 <- glm(sp.rep.108 ~ frag.5.ON + patch3.ON, 
                     data=mimicry.sample, family=poisson())
summary(SP.5000.full2)
SP.5000.frag <- glmer(sp.rep.108 ~ frag.5.ON + (1|Population), 
                      data=mimicry.sample, family=poisson())
summary(SP.5000.frag)
plot(SP.5000.frag)
plot(residuals(SP.5000.frag,type = "response")) ## response residuals
plot(residuals(SP.5000.frag,type = "pearson")) ## pearson residuals
plot(residuals(SP.5000.frag,type = "deviance")) ## deviance residuals
plot(residuals(SP.5000.frag) ~ predict(SP.5000.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.5000.frag) ~ predict(SP.5000.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.5000.frag))
qqline(resid(SP.5000.frag))

AR.5000.full <- lmer(all2sp ~ frag.5.ON + patch3.ON + (1|Population), 
                     data=mimicry.sample)
summary(AR.5000.full)
AR.5000.frag <- lmer(all2sp ~ frag.5.ON + (1|Population), 
                     data=mimicry.sample)
summary(AR.5000.frag)

VR.5000.full <- lmer(voc2sp ~ frag.5.ON + patch3.ON + (1|Population), 
                     data=mimicry.sample)
summary(VR.5000.full)
VR.5000.frag <- lmer(voc2sp ~ frag.5.ON + (1|Population), 
                     data=mimicry.sample)
summary(VR.5000.frag)

NVR.5000.full <- lmer(nonvoc2sp ~ frag.5.ON + patch3.ON + (1|Population), 
                      data=mimicry.sample)
summary(NVR.5000.full)
NVR.5000.frag <- lmer(nonvoc2sp ~ frag.5.ON + (1|Population), 
                      data=mimicry.sample)
summary(NVR.5000.frag)
NVR.5000.frag.a <- lm(nonvoc2sp ~ frag.5.ON, 
                      data=mimicry.sample)
summary(NVR.5000.frag.a)

### 10 km ####
M.10.full <- glmer(m.rep.149 ~ frag.10.ON + patch3.ON + (1|Population), 
                   data=mimicry.sample, family=poisson())
summary(M.10.full)
M.10.frag <- glmer(m.rep.149 ~ frag.10.ON + (1|Population), 
                   data=mimicry.sample, family=poisson())
summary(M.10.frag)

V.10.full <- glmer(voc.rep.149 ~ frag.10.ON + patch3.ON + (1|Population), 
                   data=mimicry.sample, family=poisson())
summary(V.10.full)
V.10.full2 <- glm(voc.rep.149 ~ frag.10.ON + patch3.ON, 
                  data=mimicry.sample, family=poisson())
summary(V.10.full2)
V.10.frag <- glmer(voc.rep.149 ~ frag.10.ON + (1|Population), 
                   data=mimicry.sample, family=poisson())
summary(V.10.frag)
plot(V.10.frag)
plot(residuals(V.10.frag,type = "response")) ## response residuals
plot(residuals(V.10.frag,type = "pearson")) ## pearson residuals
plot(residuals(V.10.frag,type = "deviance")) ## deviance residuals
plot(residuals(V.10.frag) ~ predict(V.10.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(V.10.frag) ~ predict(V.10.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(V.10.frag))
qqline(resid(V.10.frag))

V.10.frag.a <- glm(voc.rep.149 ~ frag.10.ON, 
                   data=mimicry.sample, family=poisson())
summary(V.10.frag.a)

SP.10.full <- glmer(sp.rep.108 ~ frag.10.ON + patch3.ON + (1|Population), 
                    data=mimicry.sample, family=poisson())
summary(SP.10.full)
SP.10.full2 <- glm(sp.rep.108 ~ frag.10.ON + patch3.ON, 
                   data=mimicry.sample, family=poisson())
summary(SP.10.full2)
SP.10.frag <- glmer(sp.rep.108 ~ frag.10.ON + (1|Population), 
                    data=mimicry.sample, family=poisson())
summary(SP.10.frag)
plot(SP.10.frag)
plot(residuals(SP.10.frag,type = "response")) ## response residuals
plot(residuals(SP.10.frag,type = "pearson")) ## pearson residuals
plot(residuals(SP.10.frag,type = "deviance")) ## deviance residuals
plot(residuals(SP.10.frag) ~ predict(SP.10.frag,type="response"),
     xlab=expression(hat(mu)),ylab="Deviance residuals",pch=20,col="red")
plot(residuals(SP.10.frag) ~ predict(SP.10.frag,type="link"),
     xlab=expression(hat(eta)),ylab="Deviance residuals",pch=20,col="blue")
qqnorm(resid(SP.10.frag))
qqline(resid(SP.10.frag))

AR.10.full <- lmer(all2sp ~ frag.10.ON + patch3.ON + (1|Population), 
                   data=mimicry.sample)
summary(AR.10.full)
AR.10.frag <- lmer(all2sp ~ frag.10.ON + (1|Population), 
                   data=mimicry.sample)
summary(AR.10.frag)

VR.10.full <- lmer(voc2sp ~ frag.10.ON + patch3.ON + (1|Population), 
                   data=mimicry.sample)
summary(VR.10.full)
VR.10.frag <- lmer(voc2sp ~ frag.10.ON + (1|Population), 
                   data=mimicry.sample)
summary(VR.10.frag)

NVR.10.full <- lmer(nonvoc2sp ~ frag.10.ON + patch3.ON + (1|Population), 
                    data=mimicry.sample)
summary(NVR.10.full)
NVR.10.frag <- lmer(nonvoc2sp ~ frag.10.ON + (1|Population), 
                    data=mimicry.sample)
summary(NVR.10.frag)
NVR.10.frag.a <- lm(nonvoc2sp ~ frag.10.ON, 
                    data=mimicry.sample)
summary(NVR.10.frag.a)

#Figure 4 ####
#Plot results from models

#create simple models
V.patch2 <- glmer(voc.rep.149 ~ patch3 + (1|Population), 
                  data=mimicry.sample, family=poisson())
summary(V.patch2)
SP.patch2 <- glmer(sp.rep.108 ~ patch3 + (1|Population), 
                   data=mimicry.sample, family=poisson())
summary(SP.patch2)
v2sp.patch2 <- lmer(voc2sp ~ patch3 + (1|Population), 
                    data=mimicry.sample)
summary(v2sp.patch2)
V.1000.frag2 <- glmer(voc.rep.149 ~ frag.1 + (1|Population), 
                      data=mimicry.sample, family=poisson())
summary(V.1000.frag2)
SP.1000.frag2 <- glmer(sp.rep.108 ~ frag.1 + (1|Population), 
                       data=mimicry.sample, family=poisson())
summary(SP.1000.frag2)
VR.1000.frag2 <- lmer(voc2sp ~ frag.1 + (1|Population), 
                      data=mimicry.sample)
summary(VR.1000.frag2)

library(jtools)
library(ggplot2)

View(mimicry.sample)
#vocalisations ~ patch size
p1 <- effect_plot(V.patch2, pred = patch3, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = patch3, y = voc.rep.149, colour = Population),
             position = position_jitter(w = 15, h = 0.2),
             size = 2) +
  ylim(9, 23) +
  xlim(0, 530) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 2, size = 10)) +
  xlab("") +
  ylab("Vocalisations mimicked") +
  theme(plot.margin = unit(c(1,0.1,-0.5,0.5), "lines"))
p1
#species ~ patch size
p2 <- effect_plot(SP.patch2, pred = patch3, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = patch3, y = sp.rep.108, colour = Population),
             position = position_jitter(w = 15, h = 0.2),
             size = 2) +
  ylim(2, 13) +
  xlim(0, 530) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 2, size = 10)) +
  xlab("") +
  ylab("Species mimicked")+
  theme(plot.margin = unit(c(0.25,0.1,0.25,0.5), "lines")) #0.25,0.1,0.25,1
p2
#ratio ~ patch size
p3 <- effect_plot(v2sp.patch2, pred = patch3, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = patch3, y = voc2sp, colour = Population),
             position = position_jitter(w = 15, h = 0.1),
             size = 2) +
  ylim(1.5, 4.5) +
  xlim(0, 530) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 5, size = 10)) +
  xlab(expression(Patch~size~(km^2))) +
  ylab("Vocalisations per species mimicked")+
  theme(plot.margin = unit(c(-0.5,0.1,1,1), "lines"))
p3
#vocalisations ~ frag
p4 <- effect_plot(V.1000.frag2, pred = frag.1, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = frag.1, y = voc.rep.149, colour = Population),
             position = position_jitter(w = 0, h = 0.2),
             size = 2) +
  ylim(9, 23) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 1, size = 10)) +
  xlab("") +
  ylab("")+
  theme(plot.margin = unit(c(1,1.1,-0.5,-0.5), "lines"))
p4
#species ~ frag
p5 <- effect_plot(SP.1000.frag2, pred = frag.1, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = frag.1, y = sp.rep.108, colour = Population),
             position = position_jitter(w = 0, h = 0.2),
             size = 2) +
  ylim(2, 13) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 1, size = 10)) +
  xlab("") +
  ylab("")+
  theme(plot.margin = unit(c(0.25,1.1,0.25,-0.5), "lines"))
p5
#ratio ~ frag
p6 <- effect_plot(VR.1000.frag2, pred = frag.1, interval = TRUE, plot.points = FALSE) +
  geom_point(data = mimicry.sample, 
             mapping = aes(x = frag.1, y = voc2sp, colour = Population),
             position = position_jitter(w = 0, h = 0.1),
             size = 2) +
  ylim(1.5, 4.5) +
  theme_bw() + 
  scale_color_manual(
    values=c("#D95F02", "#7570B3", "#E7298A","#66A61E", 
             "#1B9E77", "#E6AB02", "#A6761D")) +
  theme(#legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(vjust = 0.5, size = 10),
    axis.title.y = element_text(vjust = 1, size = 10)) +
  xlab("Proportion suitable habitat within 1km") +
  ylab("") +
  theme(plot.margin = unit(c(-0.5,1.1,1,0), "lines"))
p6
#plot all together
library(ggpubr)
ggarrange(p1, p4, p2, p5, p3, p6, ncol=2, nrow=3, common.legend=TRUE)

# END OF CODE #