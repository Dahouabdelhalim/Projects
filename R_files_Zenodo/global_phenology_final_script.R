setwd("C:/Users/ak222730/Desktop/RWDLAB/Phenology")
data = read.csv('Mills_data_final_2018.csv', header = T, sep = ',')

data = subset(data, data$colorSymbol == 1 | data$colorSymbol == 2)
levels(data$colorSymbol) <- c(levels(data$colorSymbol), 0) 
data$colorSymbol[data$colorSymbol==2]  <- 0
#make 2=0 keep 1 an delete 3 and 0
#1=white, 0=brown 

names(data)[names(data)=="snow"] <- "snow.cover"
names(data)[names(data)=="humanfp"] <- "h.footprint"
names(data)[names(data)=="dist_to_coast"] <- "dist.to.shore"
data$dist.to.shore.km<-data$dist.to.shore/1000
#rename columns

data$alt.sc <- (data$alt-mean(data$alt))/sd(data$alt)
data$bio_2.sc <- (data$bio_2-mean(data$bio_2))/sd(data$bio_2)
data$bio_3.sc <- (data$bio_3-mean(data$bio_3))/sd(data$bio_3)
data$bio_15.sc <- (data$bio_15-mean(data$bio_15))/sd(data$bio_15)
data$snow.cover.sc <- (data$snow.cover-mean(data$snow.cover))/sd(data$snow.cover)
data$h.footprint.sc <- (data$h.footprint-mean(data$h.footprint))/sd(data$h.footprint)
data$dist.to.shore.km.sc <- (data$dist.to.shore.km-mean(data$dist.to.shore.km))/sd(data$dist.to.shore.km)
#standardize covariates 

table(data$spsCode)

data.8spp = data[data$spsCode == 'LEPAM' |
   data$spsCode == 'LEPBR' | data$spsCode == 'LEPTI' | data$spsCode == 'LEPTO' |
   data$spsCode == 'MUSER' | data$spsCode == 'MUSFR' | data$spsCode == 'MUSNI' | data$spsCode == 'VULLA',]
#subset to 8 species

install.packages("lme4")
library("lme4")
install.packages("raster")
library("raster")
install.packages("rgdal")
library("rgdal")
install.packages("RArcInfo")
library("RArcInfo")
install.packages("MuMIn")
library("MuMIn")
install.packages("MASS")
library('MASS')
#install packages
####################################################################################
#
###############################MAKE FULL MODEL WITH ALL SPECIES#####################
#
####################################################################################


model1 = glmer(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc + (1|spsCode), family = binomial, data = data.8spp, na.action = na.pass)
summary(model1)

exp(coef(model1)[[1]] [1,])
#get odds ratio
cc1<-confint(model1, parm="beta_")
ctab1<-cbind(est=fixef(model1), cc1)
exp(ctab1)
###calc confidence intervals (Bolkers method)


##################################AIC WEIGHTS################################
results <- dredge(model1)
importance(results)
#find support for covariates using AIC 


model3 = glmer(colorSymbol ~ bio_2.sc + bio_3.sc + snow.cover.sc + (1|spsCode), family = binomial, data = data.8spp, na.action = na.pass)
summary(model3)
###make model with only 3 predictors
exp(coef(model3)[[1]] [1,])

cc<-confint(model3, parm="beta_")
ctab<-cbind(est=fixef(model3), cc)
exp(ctab)
###calc confidence intervals (Bolkers methods)


###############################################################################
#
###############################################################################
#
###################SPECIES SPECIFIC PREDICTION MAPS############################
#
###############################################################################
#
###############################################################################


#####################PREDICTION #####################################
###need to download species range files from IUCN red list website
###http://www.iucnredlist.org/

################### making sure everything is in the same projection and resolution ############
newproj <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"


rasstack <- stack("rasterstack4cov.tif")
#read back in the rasstack
names(rasstack)<-c("snow.cover", "bio_2", "bio_3","alt")
#rename files in rasstack


#################BEGIN SPECIES SPECIFIC MODELS################################

###############################LEPUS BRACHYURUS################################
data.lepbr = subset(data.8spp, data.8spp$spsCode =='LEPBR')
#subset data

four.lepbr = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepbr, na.action = na.pass) 
summary(four.lepbr)
#make model

fourlepbr <- predict(rasstack, four.lepbr, type='response', progress='window')
plot(fourlepbr)
#make raster

###### cropping to a range
lepbr.shp<-readOGR("LEPBR_range.shp", layer = "LEPBR_range")
###need to download species range file from IUCN website and read it in
projection(lepbr.shp) <- newproj
lepbr1 <- crop(fourlepbr, extent(lepbr.shp))
plot(lepbr1)
lepbr2 = mask(lepbr1, lepbr.shp) 
plot(lepbr2)

#######exporting raster to read in ArcGIS
writeRaster(lepbr2,  filename='FINAL3lepbr', format='GTiff',overwrite=TRUE)
writeRaster(fourlepbr,  filename='FINAL3lepbr_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################


###############################LEPUS TIMIDUS################################
data.lepti = subset(data.8spp, data.8spp$spsCode =='LEPTI')

four.lepti = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepti, na.action = na.pass) 
summary(four.lepti)

fourlepti <- predict(rasstack, four.lepti, type='response', progress='window')
plot(fourlepti)

###### cropping to a range
lepti.shp<-readOGR("LEPTI_range.shp", layer = "LEPTI_range")
projection(lepti.shp) <- newproj
lepti1 <- crop(fourlepti, extent(lepti.shp))
plot(lepti1)
lepti2 = mask(lepti1, lepti.shp) 
plot(lepti2)
writeRaster(lepti2,  filename='FINAL3lepti', format='GTiff',overwrite=TRUE)
writeRaster(fourlepti,  filename='FINAL3lepti_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################LEPUS AMERICANUS################################
data.lepam = subset(data.8spp, data.8spp$spsCode =='LEPAM')

four.lepam = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepam, na.action = na.pass) 
summary(four.lepam)

fourlepam <- predict(rasstack, four.lepam, type='response', progress='window')
plot(fourlepam)

###### cropping to a range
lepam.shp<-readOGR("LEPAM_range.shp", layer = "LEPAM_range")
projection(lepam.shp) <- newproj
lepam1 <- crop(fourlepam, extent(lepam.shp))
plot(lepam1)
lepam2 = mask(lepam1, lepam.shp) 
plot(lepam2)
writeRaster(lepam2,  filename='FINAL3lepam', format='GTiff',overwrite=TRUE)
writeRaster(fourlepam,  filename='FINAL3lepam_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################LEPUS TOWNSENDII################################
data.lepto = subset(data.8spp, data.8spp$spsCode =='LEPTO')

four.lepto = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepto, na.action = na.pass) 
summary(four.lepto)

fourlepto <- predict(rasstack, four.lepto, type='response', progress='window')
plot(fourlepto)

###### cropping to a range
lepto.shp<-readOGR("LEPTO_range.shp", layer = "LEPTO_range")
projection(lepto.shp) <- newproj
lepto1 <- crop(fourlepto, extent(lepto.shp))
plot(lepto1)
lepto2 = mask(lepto1, lepto.shp) 
plot(lepto2)
writeRaster(lepto2,  filename='FINAL3lepto', format='GTiff',overwrite=TRUE)
writeRaster(fourlepto,  filename='FINAL3lepto_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################MUSTELA NIVALIS################################
data.musni = subset(data.8spp, data.8spp$spsCode =='MUSNI')


four.musni = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.musni, na.action = na.pass) 
summary(four.musni)

fourmusni <- predict(rasstack, four.musni, type='response', progress='window')
plot(fourmusni)

###### cropping to a range
musni.shp<-readOGR("MUSNI_range.shp", layer = "MUSNI_range")
projection(musni.shp) <- newproj
musni1 <- crop(fourmusni, extent(musni.shp))
plot(musni1)
musni2 = mask(musni1, musni.shp) 
plot(musni2)
writeRaster(musni2,  filename='FINAL3musni', format='GTiff',overwrite=TRUE)
writeRaster(fourmusni,  filename='FINAL3musni_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################MUSTELA ERMINEA################################
data.muser = subset(data.8spp, data.8spp$spsCode =='MUSER')
#data.muser = read.csv('muser_edited.csv', header = T, sep = ',')

four.muser = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.muser, na.action = na.pass) 
summary(four.muser)

fourmuser <- predict(rasstack, four.muser, type='response', progress='window')
plot(fourmuser)

###### cropping to a range
muser.shp<-readOGR("MUSER_range.shp", layer = "MUSER_range")
projection(muser.shp) <- newproj
muser1 <- crop(fourmuser, extent(muser.shp))
plot(muser1)
muser2 = mask(muser1, muser.shp) 
plot(muser2)
writeRaster(muser2,  filename='FINAL3muser', format='GTiff',overwrite=TRUE)
writeRaster(fourmuser,  filename='FINAL3muser_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################MUSTELA FRENATA################################
data.musfr = subset(data.8spp, data.8spp$spsCode =='MUSFR')

four.musfr = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.musfr, na.action = na.pass) 
summary(four.musfr)

fourmusfr <- predict(rasstack, four.musfr, type='response', progress='window')
plot(fourmusfr)

###### cropping to a range
musfr.shp<-readOGR("MUSFR_range.shp", layer = "MUSFR_range")
projection(musfr.shp) <- newproj
musfr1 <- crop(fourmusfr, extent(musfr.shp))
plot(musfr1)
musfr2 = mask(musfr1, musfr.shp) 
plot(musfr2)
writeRaster(musfr2,  filename='FINAL3musfr', format='GTiff',overwrite=TRUE)
writeRaster(fourmusfr,  filename='FINAL3musfr_uncropped', format='GTiff',overwrite=TRUE)
###################################END PREDICTION################################

###############################VULPES LAGOPUS################################
data.vulla = subset(data.8spp, data.8spp$spsCode =='VULLA')

four.vulla = glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.vulla, na.action = na.pass) 
summary(four.vulla)

fourvulla <- predict(rasstack, four.vulla, type='response', progress='window')
plot(fourvulla)

###### cropping to a range
vulla.shp<-readOGR("VULLA_range.shp", layer = "VULLA_range")
projection(vulla.shp) <- newproj
vulla1 <- crop(fourvulla, extent(vulla.shp))
plot(vulla1)
vulla2 = mask(vulla1, vulla.shp) 
plot(vulla2)
writeRaster(vulla2,  filename='FINAL3vulla', format='GTiff',overwrite=TRUE)
writeRaster(fourvulla,  filename='FINAL3vulla_uncropped', format='GTiff',overwrite=TRUE)

###################################END PREDICTION################################



###############################################################################
#
###############################################################################
#
###########################GRAPH FOR %WHITE VS SNOW#############################
#
###############################################################################
#
###############################################################################


###############################LEPUS BRACHYURUS################################
data.lepbr <- subset(data.8spp, data.8spp$spsCode =='LEPBR')
four.lepbr <- glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepbr, na.action = na.pass)
###############################LEPUS TOWNSENDII################################
data.lepto <- subset(data.8spp, data.8spp$spsCode =='LEPTO')
four.lepto <- glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.lepto, na.action = na.pass) 
one.lepto <- glm(colorSymbol ~ snow.cover, family = binomial, data = data.lepto) 
###############################MUSTELA NIVALIS################################
data.musni <- subset(data.8spp, data.8spp$spsCode =='MUSNI')
four.musni <- glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.musni, na.action = na.pass) 
one.musni <- glm(colorSymbol ~ snow.cover, family = binomial, data = data.musni) 
###############################MUSTELA FRENATA################################
data.musfr <- subset(data.8spp, data.8spp$spsCode =='MUSFR')
four.musfr <- glm(colorSymbol ~ snow.cover + bio_2 + bio_3, family = binomial, data = data.musfr, na.action = na.pass) 
one.musfr <- glm(colorSymbol ~ snow.cover, family = binomial, data = data.musfr) 
#########make models for the 4 species used in the graph

#########################################################################################
############################plot response curves#########################################
#########################################################################################

range(data.lepbr$snow.cover)
range(data.8spp$snow.cover)


xsnow<-seq(0,375,1)
pred.lepbr<-predict.glm(four.lepbr, list(snow.cover=xsnow, bio_2=rep(mean(data.lepbr$bio_2), length(xsnow)), bio_3=rep(mean(data.lepbr$bio_3), length(xsnow))), type='response')
pred.lepto<-predict.glm(four.lepto, list(snow.cover=xsnow, bio_2=rep(mean(data.lepto$bio_2), length(xsnow)), bio_3=rep(mean(data.lepto$bio_3), length(xsnow))), type='response')
pred.musni<-predict.glm(four.musni, list(snow.cover=xsnow, bio_2=rep(mean(data.musni$bio_2), length(xsnow)), bio_3=rep(mean(data.musni$bio_3), length(xsnow))), type='response')
pred.musfr<-predict.glm(four.musfr, list(snow.cover=xsnow, bio_2=rep(mean(data.musfr$bio_2), length(xsnow)), bio_3=rep(mean(data.musfr$bio_3), length(xsnow))), type='response')

plot(xsnow, pred.lepbr, type="lines", xlab="Number of snow covered days", ylab=("P(white)"), main="")
polygon(x=c(0,0,375,375),y=c(0,.2,.2,0),col="sandybrown")
polygon(x=c(0,0,375,375),y=c(0.8,1,1,0.8),col=gray(0.97))
lines(xsnow, pred.lepbr, lwd=2, col="blue")
lines(xsnow, pred.musfr, lwd=2, col="red")
lines(xsnow, pred.musni, lwd=2, col="green")
lines(xsnow, pred.lepto, lwd=2, col="purple")
legend(x="topright", legend=c("Japanese Hare","White-tailed Jackrabbit", "Long-tailed Weasel", "Least Weasel"),col=c("blue","purple", "red", "green"), lwd=2, cex=1,lty=1,pt.cex=2,inset=c(.04,0.23),bty="n")  



###############################################################################
#
###############################################################################
#
########################### 5 FOLD CROSS VALIDATION############################
#
######################(run seperately for each species)########################
#
###############################################################################




###########################################LEPAM#################################################

data.lepam = subset(data.8spp, data.8spp$spsCode =='LEPAM')
#subset data
k = 5
#pick k (number of subsets to divide data into)

data.lepam$id <- sample(1:k, nrow(data.lepam), replace = TRUE)# adds column to data indicating which of 5 groups
list <- 1:k

prediction.lepam <- data.frame()
testsetCopy.lepam <- data.frame()
# prediction and test set data frames that we add to with each iteration over


for(i in 1:k){
  trainingset.lepam <- subset(data.lepam, id %in% list[-i])
  testset.lepam <- subset(data.lepam, id %in% c(i))
  # remove rows with id i from dataframe to create training set
  # select rows with id i to create test set
  mymodel.lepam <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc + dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.lepam, na.action = na.pass) 
  ###full model
  ###to run reduced model change to: mymodel.lepam <- glm(colorSymbol ~ bio_2.sc + bio_3.sc+ snow.cover.sc, family = binomial, data = trainingset.lepam, na.action = na.pass) 
  
  temp.lepam <- as.data.frame(predict.glm(mymodel.lepam,testset.lepam, type="response"))
  #remove response column 1
  
  prediction.lepam <- rbind(prediction.lepam, temp.lepam)
  # append this iteration's predictions to the end of the prediction data frame
  
  testsetCopy.lepam <- rbind(testsetCopy.lepam, as.data.frame(testset.lepam$colorSymbol))
  # append this iteration's test set to the test set copy data frame
  
}

############function for k fold#############

# add predictions and actual values
result.lepam <- cbind(prediction.lepam, testsetCopy.lepam[, 1])
names(result.lepam) <- c("Predicted", "Actual")
result.lepam$Difference <- abs(result.lepam$Actual - result.lepam$Predicted)
summary(result.lepam$Difference)
#summarize results from the difference between the predicted and actual values

dim(result.lepam)[[1]]
table(result.lepam$Difference<.5)["TRUE"]/dim(result.lepam)[[1]]
# calculate proportion of times model predicted actual value 

###########################################LEPTI#################################################

data.lepti = subset(data.8spp, data.8spp$spsCode =='LEPTI')
k = 5

data.lepti$id <- sample(1:k, nrow(data.lepti), replace = TRUE)
list <- 1:k

prediction.lepti <- data.frame()
testsetCopy.lepti <- data.frame()

for(i in 1:k){
  
  trainingset.lepti <- subset(data.lepti, id %in% list[-i])
  testset.lepti <- subset(data.lepti, id %in% c(i))
  
  mymodel.lepti <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.lepti, na.action = na.pass) 
  
  temp.lepti <- as.data.frame(predict.glm(mymodel.lepti,testset.lepti, type="response"))
  
  prediction.lepti <- rbind(prediction.lepti, temp.lepti)
  
  testsetCopy.lepti <- rbind(testsetCopy.lepti, as.data.frame(testset.lepti$colorSymbol))
  
}

result.lepti <- cbind(prediction.lepti, testsetCopy.lepti[, 1])
names(result.lepti) <- c("Predicted", "Actual")
result.lepti$Difference <- abs(result.lepti$Actual - result.lepti$Predicted)
summary(result.lepti$Difference)

dim(result.lepti)[[1]]
table(result.lepti$Difference<.5)["TRUE"]/dim(result.lepti)[[1]]

###########################################LEPBR#################################################

data.lepbr = subset(data.8spp, data.8spp$spsCode =='LEPBR')
k = 5

data.lepbr$id <- sample(1:k, nrow(data.lepbr), replace = TRUE)
list <- 1:k

prediction.lepbr <- data.frame()
testsetCopy.lepbr <- data.frame()

for(i in 1:k){
  trainingset.lepbr <- subset(data.lepbr, id %in% list[-i])
  testset.lepbr <- subset(data.lepbr, id %in% c(i))
  
  mymodel.lepbr <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.lepbr, na.action = na.pass) 
  
  temp.lepbr <- as.data.frame(predict.glm(mymodel.lepbr,testset.lepbr, type="response"))
  
  prediction.lepbr <- rbind(prediction.lepbr, temp.lepbr)
  
  testsetCopy.lepbr <- rbind(testsetCopy.lepbr, as.data.frame(testset.lepbr$colorSymbol))
  
}
result.lepbr <- cbind(prediction.lepbr, testsetCopy.lepbr[, 1])
names(result.lepbr) <- c("Predicted", "Actual")
result.lepbr$Difference <- abs(result.lepbr$Actual - result.lepbr$Predicted)
summary(result.lepbr$Difference)

dim(result.lepbr)[[1]]
table(result.lepbr$Difference<.5)["TRUE"]/dim(result.lepbr)[[1]]

###########################################LEPTO#################################################

data.lepto = subset(data.8spp, data.8spp$spsCode =='LEPTO')
k = 5

data.lepto$id <- sample(1:k, nrow(data.lepto), replace = TRUE)
list <- 1:k
prediction.lepto <- data.frame()
testsetCopy.lepto <- data.frame()

for(i in 1:k){
  trainingset.lepto <- subset(data.lepto, id %in% list[-i])
  testset.lepto <- subset(data.lepto, id %in% c(i))
  
  mymodel.lepto <- glm(colorSymbol ~  alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.lepto, na.action = na.pass) 
  
  temp.lepto <- as.data.frame(predict.glm(mymodel.lepto,testset.lepto, type="response"))
  
    prediction.lepto <- rbind(prediction.lepto, temp.lepto)
  
  testsetCopy.lepto <- rbind(testsetCopy.lepto, as.data.frame(testset.lepto$colorSymbol))
  
}

result.lepto <- cbind(prediction.lepto, testsetCopy.lepto[, 1])
names(result.lepto) <- c("Predicted", "Actual")
result.lepto$Difference <- abs(result.lepto$Actual - result.lepto$Predicted)
summary(result.lepto$Difference)

dim(result.lepto)[[1]]
table(result.lepto$Difference<.5)["TRUE"]/dim(result.lepto)[[1]]

###########################################MUSNI#################################################

data.musni = subset(data.8spp, data.8spp$spsCode =='MUSNI')
k = 5

data.musni$id <- sample(1:k, nrow(data.musni), replace = TRUE)
list <- 1:k

prediction.musni <- data.frame()
testsetCopy.musni <- data.frame()


for(i in 1:k){
  trainingset.musni <- subset(data.musni, id %in% list[-i])
  testset.musni <- subset(data.musni, id %in% c(i))
  
  mymodel.musni <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.musni, na.action = na.pass) 
  
  temp.musni <- as.data.frame(predict.glm(mymodel.musni,testset.musni, type="response"))
  
    prediction.musni <- rbind(prediction.musni, temp.musni)
  
  testsetCopy.musni <- rbind(testsetCopy.musni, as.data.frame(testset.musni$colorSymbol))
  
}
result.musni <- cbind(prediction.musni, testsetCopy.musni[, 1])
names(result.musni) <- c("Predicted", "Actual")
result.musni$Difference <- abs(result.musni$Actual - result.musni$Predicted)
summary(result.musni$Difference)

dim(result.musni)[[1]]
table(result.musni$Difference<.5)["TRUE"]/dim(result.musni)[[1]]

###########################################MUSER#################################################

data.muser = subset(data.8spp, data.8spp$spsCode =='MUSER')
k = 5

data.muser$id <- sample(1:k, nrow(data.muser), replace = TRUE)
list <- 1:k
prediction.muser <- data.frame()
testsetCopy.muser <- data.frame()

for(i in 1:k){
  trainingset.muser <- subset(data.muser, id %in% list[-i])
  testset.muser <- subset(data.muser, id %in% c(i))
  
  mymodel.muser <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.muser, na.action = na.pass) 
  
  temp.muser <- as.data.frame(predict.glm(mymodel.muser,testset.muser, type="response"))
  
  prediction.muser <- rbind(prediction.muser, temp.muser)
  
  testsetCopy.muser <- rbind(testsetCopy.muser, as.data.frame(testset.muser$colorSymbol))
  
}
result.muser <- cbind(prediction.muser, testsetCopy.muser[, 1])
names(result.muser) <- c("Predicted", "Actual")
result.muser$Difference <- abs(result.muser$Actual - result.muser$Predicted)
summary(result.muser$Difference)

dim(result.muser)[[1]]
table(result.muser$Difference<.5)["TRUE"]/dim(result.muser)[[1]]

###########################################MUSFR#################################################

data.musfr = subset(data.8spp, data.8spp$spsCode =='MUSFR')
k = 5

data.musfr$id <- sample(1:k, nrow(data.musfr), replace = TRUE)
list <- 1:k
prediction.musfr <- data.frame()
testsetCopy.musfr <- data.frame()

for(i in 1:k){
  trainingset.musfr <- subset(data.musfr, id %in% list[-i])
  testset.musfr <- subset(data.musfr, id %in% c(i))
  
  mymodel.musfr <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.musfr, na.action = na.pass) 
  
  temp.musfr <- as.data.frame(predict.glm(mymodel.musfr,testset.musfr, type="response"))
  
  prediction.musfr <- rbind(prediction.musfr, temp.musfr)
  
  testsetCopy.musfr <- rbind(testsetCopy.musfr, as.data.frame(testset.musfr$colorSymbol))
  
}
result.musfr <- cbind(prediction.musfr, testsetCopy.musfr[, 1])
names(result.musfr) <- c("Predicted", "Actual")
result.musfr$Difference <- abs(result.musfr$Actual - result.musfr$Predicted)
summary(result.musfr$Difference)

dim(result.musfr)[[1]]
table(result.musfr$Difference<.5)["TRUE"]/dim(result.musfr)[[1]]

###########################################VULLA#################################################

data.vulla = subset(data.8spp, data.8spp$spsCode =='VULLA')
k = 5

data.vulla$id <- sample(1:k, nrow(data.vulla), replace = TRUE)
list <- 1:k
prediction.vulla <- data.frame()
testsetCopy.vulla <- data.frame()


for(i in 1:k){
  trainingset.vulla <- subset(data.vulla, id %in% list[-i])
  testset.vulla <- subset(data.vulla, id %in% c(i))
  
  mymodel.vulla <- glm(colorSymbol ~ alt.sc + bio_2.sc + bio_3.sc+bio_15.sc + snow.cover.sc +dist.to.shore.km.sc + h.footprint.sc, family = binomial, data = trainingset.vulla, na.action = na.pass) 
  
    temp.vulla <- as.data.frame(predict.glm(mymodel.vulla,testset.vulla, type="response"))
  
  prediction.vulla <- rbind(prediction.vulla, temp.vulla)

  testsetCopy.vulla <- rbind(testsetCopy.vulla, as.data.frame(testset.vulla$colorSymbol))
  
}
result.vulla <- cbind(prediction.vulla, testsetCopy.vulla[, 1])
names(result.vulla) <- c("Predicted", "Actual")
result.vulla$Difference <- abs(result.vulla$Actual - result.vulla$Predicted)

dim(result.vulla)[[1]]
table(result.vulla$Difference<.5)["TRUE"]/dim(result.vulla)[[1]]