install.packages("raster")
install.packages("phyloclim")
install.packages("maptools")
install.packages("dismo")
install.packages("ape")
install.packages("geiger")
install.packages("rgdal")
install.packages("rgeos")
install.packages("GISTools")
install.packages("rJava")
install.packages("ENMeval")
install.packages("spocc")
install.packages("rJava")
install.packages("tiff")

library(raster)
library(phyloclim)
library(maptools)
library(dismo)
library(ape)
library(geiger)
library(rgdal)
library(rgeos)
library(GISTools)
library(rJava)
library(ENMeval)
library(spocc)
library(rJava)
library(tiff)

###############################################################################
##################### Worldclim variables PRESENT #############################
###############################################################################

#Read files
files<-stack(list.files(path = "PathTo/Bioclim",pattern='bil', full.names=T)) ### stack all rasters in Bioclim folder
e<-extent(-78.352, -34.099, -20.812, 12.839) #restrict space to South America
files<-crop(files,e) #Crop to South America
plot(files[[1]])
dev.off()

projection(files)<- CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0') # coodinate reference system for WGS84
plot(files[[1]])

####remove highly correlated variables 
test<-getValues(files)
cor.matrix<-cor(test, use="complete.obs")
write.csv(cor.matrix,'cor.matrix.csv') #open csv to check for correlation among variables
files.sub<- dropLayer(files, c(1,2,5,8,9,10,11,13,16,17,18,19)) #### remove selected layers
test2<-getValues(files.sub)
cor.matrix2<- cor(test2, use="complete.obs") 
write.csv(cor.matrix2,'cor.matrix2.csv') #check correlation among remaining variables

names(files.sub)
plot(files.sub[[2]])

predictors<-files.sub

###############################################################################
##################### Worldclim variables 6k ##################################
###############################################################################

#Read files
files.6k <- stack(list.files(path = "PathTo/6k",pattern='tif', full.names=T))
plot(files.6k[[1]])
projection(files.6k)<-CRS('+proj=longlat  +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs' )
e<-extent(-78.352, -34.099, -20.812, 12.839)
pred6k<-crop(files.6k,e)
plot(pred6k[[7]])
dev.off()

#remove high correlated variables
pred6k<- dropLayer(pred6k, c(1,2,5,8,9,10,11,13,16,17,18,19))
names(pred6k)

###############################################################################
##################### Worldclim variables 21k #################################
###############################################################################

#Read files
files.21k <- stack(list.files(path = "PathTo/21k", pattern='tif', full.names=T))
plot(files.21k[[1]])
e<-extent(-78.352, -34.099, -20.812, 12.839)
files.21k<-crop(files.21k,e)
plot(files.21k[[1]])

pred21k<- dropLayer(files.21k,  c(1,2,5,8,9,10,11,13,16,17,18,19))
projection(pred21k)<-CRS('+proj=longlat  +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs' )
names(pred21k) 

###############################################################################
##################### Worldclim variables 130k ################################
###############################################################################

#Read files
files.130k <- stack(list.files(path = "PathTo/130k", pattern='bil', full.names=T))
plot(files.130k[[1]])
e<-extent(-78.352, -34.099, -20.812, 12.839)
files.130K<-crop(files.130k,e)
plot(files.130K[[1]])

pred130k<- dropLayer(files.130K, c(1,2,5,8,9,10,11,13,16,17,18,19))#aqui não tem a altitude, os valores serão diferentes dos anteriores
pred130k <- aggregate(pred130k,fact=5)
pred130k <- resample(pred130k, predictors)
names(pred130k)


###############################################################################
############################## OCCURRENCES ####################################
###############################################################################

occ.sps <- list.files('PathTo/occ_thin',pattern=".csv") #set path to the directory with occurrences after thinning
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\\\.csv")))


###############################################################################
################ RUN MODEL FOR PRESENT AND PAST CONDITIONS ####################
###############################################################################

for (i in 1:length(splist)){
  sp.file <- read.csv(paste('PathTo/occ_thin/', occ.sps[i],sep=""),h=T) ### read sp occurrence
  occs <- sp.file[,2:3] ## select lat long columns  
  
  envs.backg<-predictors
  
  # Randomly sample 10,000 background points from one background extent raster (only one per cell without replacement). Note: Since the raster has <10,000 pixels, you'll get a warning and all pixels will be used for background. We will be sampling from the biome variable because it is missing some grid cells, and we are trying to avoid getting background points with NA.
  bg <- randomPoints(envs.backg[[1]], n=10000)
  
  eval2.par <- ENMevaluate(occs, envs.backg, bg, 
                           method='checkerboard1',
                           RMvalues= seq(0.5, 4, 0.5),
                           fc = c("L","LQ", "H", "LQH", "LQHP", "LQHPT"),
                           parallel=FALSE,
                           algorithm='maxent.jar',
                           progbar = TRUE, 
                           updateProgress = T,
                           overlap=T)
  
  #Transform from Raw to Logistic format ('best' model)
  x<-data.frame(names(eval2.par@predictions))
  colnames(x)<-'names_models'
  x$models_n<-seq(1:48)
  z<-x[x$names_models == names(eval2.par@predictions[[which(eval2.par@results$delta.AICc==0)]]),2]
  z<-z[1]
  p <- predict(eval2.par@models[[z]], envs.backg,
               args=c("outputformat=logistic"))
  
  plot(p)  # logistic output
  
  #write results in logistic format
  writeRaster(p,
              filename=paste(getwd(),"/PathTo/maxent_log/",
                             splist[i],
                             ".asc", sep=""),overwrite=T)
  
  res_df<-eval2.par@results
  write.csv(res_df, file=paste("PathTo/res_df/",
                               splist[i],".csv",sep=''),row.names=F)
  
  #write results('best' model) in raw format
  writeRaster(eval2.par@predictions[[z]],
              filename=paste(getwd(),"/PathTo/maxent_raw/",
                             splist[i],".asc", sep=""),overwrite=T)
  
  # Model figures ('best' model)
  data(wrld_simpl)
  crs(wrld_simpl)<-'+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0' 
  
  pdf(file=paste(getwd(),'/PathTo/Models_figures/',
                 splist[i], ".pdf", sep=""),
      width=8,height=6) ### make pdf figure to explore results
  plot(p,main= splist[i])
  plot(wrld_simpl, add=T, lwd=0.5)
  points(occs,pch=20, cex=0.3, bg = "dimgray",col = "red",lwd="0.2")
  scalebar(1000, xy = c(-85,-25), type = 'bar', divs = 2, below = "km")
  dev.off()
  
  # Let's look at the model object for our "AICc optimal" model:
  aic.opt <- eval2.par@models[[z]]
  maxent_res_a<-aic.opt@results 
  colnames(maxent_res_a)<-'value'
  Parameters<-data.frame(rownames(maxent_res_a))
  colnames(Parameters)<-'Parameters'
  
  maxent_resuts<-cbind(Parameters,maxent_res_a)
  row.names(maxent_resuts)<-NULL
  maxent_resuts$Species<- rep(splist[i],nrow(maxent_resuts))
  
  write.csv(maxent_resuts[,c(3,1,2)],
            file=paste("PathTo/maxent_results/",
                       splist[i],".csv",sep=''),row.names=F)
  
  
  # plot model plus evaluations
  pdf(file=paste(getwd(),'/PathTo/Figure/',
                 splist[i], ".pdf", sep=""),
      width=12,height=7)
  par(mfrow=c(2,3))
  plot(p,main= splist[i])
  plot(wrld_simpl, add=T, lwd=0.5)
  points(occs,pch=20, cex=0.3, bg = "dimgray",col = "red",lwd="0.2")
  scalebar(1000, xy = c(-85,-25), type = 'bar', divs = 2, below = "km")  
  eval.plot(eval2.par@results,
            legend.position ='topright')
  eval.plot(eval2.par@results, 'avg.test.AUC',
            legend =F)
  eval.plot(eval2.par@results, 'avg.diff.AUC',
            legend =F)
  eval.plot(eval2.par@results, 'avg.test.orMTP',
            legend=F)
  eval.plot(eval2.par@results, 'avg.test.or10pct',
            legend =F)
  dev.off()
  
  #Project model to the past and create raster files
  ##0k##
  palmipes.0k <- predict(eval2.par@models[[z]],predictors)
  plot(palmipes.0k)
  writeRaster(palmipes.0k,
              filename=paste(getwd(),"/PathTo/result/0k/",
                             splist[i],
                             ".asc", sep=""),overwrite=T)  
  ##6k##
  palmipes.6k <- predict(eval2.par@models[[z]],pred6k)
  plot(palmipes.6k)
  writeRaster(palmipes.6k,
              filename=paste(getwd(),"/PathTo/result/6k/",
                             splist[i],
                             ".asc", sep=""),overwrite=T)    
  ##21k##
  palmipes.21k <- predict(eval2.par@models[[z]],pred21k)
  plot(palmipes.21k)
  writeRaster(palmipes.21k,
              filename=paste(getwd(),"/PathTo/result/21k/",
                             splist[i],
                             ".asc", sep=""),overwrite=T)    
  ##130k##
  palmipes.130k <- predict(eval2.par@models[[z]],pred130k)
  plot(palmipes.130k)
  writeRaster(palmipes.130k,
              filename=paste(getwd(),"/PathTo/result/130k/",
                             splist[i],
                             ".asc", sep=""),overwrite=T)    
}

############Stable areas##############

#read asc files with model from current and past conditions
palmipes0k <- raster("PathTo/result/0k/filename.asc") 
palmipes6k <- raster('PathTo/result/6k/filename.asc') 
palmipes21k <- raster('PathTo/result/21k/filename.asc')
palmipes130k <- raster('PathTo/result/130k/filename.asc') 

#Set a threshold to create binary maps
palmipes0kthr = palmipes0k >= 0.1709 ###### in "maxentResults", search for "Equal training sensitivity and specificity logistic threshold" value
plot(palmipes0kthr)

palmipes6kthr = palmipes6k >= 0.1709
plot(palmipes6kthr)

palmipes21kthr = palmipes21k >= 0.1709
plot(palmipes21kthr)

palmipes130kthr = palmipes130k >= 0.1709
plot(palmipes130kthr)

#Merge all files
refugiototal <- palmipes0kthr+palmipes6kthr+palmipes21kthr+palmipes130kthr
plot(refugiototal)

writeRaster(refugiototal,filename='08ENMeval/resultado/refugio_atf.asc', overwrite=T)
