###upadated January 2021
options(java.parameters = "-Xmx14000m") ##Check this in every computer and set accoridng to available RAM
library(ENMeval)
library(spThin)
library(rgeos)
library(sp)
library(rgdal)
library(raster)
#Setwd to the project folder
#setwd("/Users/andreapaz/Dropbox/**Tesis_PHD/Ch2/SDM/")
setwd("D:/Andrea/Ch2/")
#Load and stack predictor variables 
predictors <- stack(list.files(path="./Bioclim2_SA/",pattern='*.tif$*', full.names=TRUE))

#Load and prepare occurence data (long,lat for maxent) 


species_localities<-read.csv("All_points_Celio_dbs_jan2021.csv",h=T)
species_list<-unique(species_localities$Genus_species)


for (i in species_list){
Points<-subset(species_localities,Genus_species==i)

ID<-c(1:length(Points$Latitude))
Points_occ<-data.frame(ID,Points[,c(9,8)])
dir<-getwd()
coordinates(Points_occ)<-c('Longitude','Latitude')
writeOGR(Points_occ, dir,layer=paste("Points",i,sep="_"), driver="ESRI Shapefile")
occ<-readOGR(dir, paste("Points",i,sep="_"))
occ<-as.data.frame(occ)
Species<-rep(i,length(occ[,1]))
occ<-data.frame(Species,occ[,c(2,3)])
#occ1<-occ[,c(2,3)]
#Locality data thinning using spThin package thin.par is distance in km
presencia_thinned<-thin(occ, lat.col = "coords.x2", long.col = "coords.x1", spec.col = "Species",
                       thin.par=5, reps=5, locs.thinned.list.return = TRUE, write.files = TRUE,
                      max.files = 5, out.dir=paste(i,"thinned_data_test",sep="_"), out.base = "thinned_data",
                       write.log.file = TRUE, log.file = "spatial_thin_log.txt",
                      verbose = TRUE)
occ1<-presencia_thinned[[1]]
if (length(occ1$Longitude)>=5){ ##for les than 5 points no model is created just a file with the name of the species
  tryCatch({
  #define background as MCP or MCP with buffer
##MCP only  
backgd<-convHull(occ1)@polygons
#with buffer first project
		##Define projection and project
crs(backgd)<- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
new_coord<-"+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 +ellps=GRS67 +units=m +no_defs  "
backgd1<-spTransform(backgd, CRS(new_coord))
		##then buffer at 100km
backgd1 <- gBuffer(backgd1, width=100000)
		##then reproject
backgd<-spTransform(backgd1, crs(backgd))
bkg <- rasterize(backgd, predictors, field=1)
#plot(bkg,add=T)
bkg1<-mask(bkg,predictors$bio_1)
bkg2<-crop(bkg1,backgd)
#plot(bkg2,add=T)
#Create background points
backg <- randomPoints(bkg2, n=10000) 
#plot(backg,add=T) ##just to look at points and slows the script
####Crop predictors to species specific background area (MCP + 100km buffer)
predictors1<-crop(predictors,backgd)
predictors1<-mask(predictors1,backgd)

if (length(occ1$Longitude)>=20){ #this number can vary in this case for more than 20 does randomkflod and less jackknife

  k<-5
modelo<-ENMevaluate(occ1,predictors1 , bg.coords = backg, occ.grp = NULL,
            bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
            fc = c("L", "LQ", "H", "LQH"),
            categoricals = NULL, n.bg = NULL, method = "randomkfold",
            overlap = FALSE, aggregation.factor = c(2, 2),
            kfolds = k, bin.output = FALSE, clamp = TRUE,
            rasterPreds = T, parallel = T, numCores =18,algorithm="maxent.jar")
save(modelo,file=paste('ENM_eval_blocks_R',i,sep="_"))
}
else {
  modelo<-ENMevaluate(occ1,predictors1 , bg.coords = backg, occ.grp = NULL,
                      bg.grp = NULL, RMvalues = seq(0.5, 4, 0.5),
                      fc = c("L", "LQ", "H", "LQH"),
                      categoricals = NULL, n.bg = NULL, method = "jackknife",
                      overlap = FALSE, aggregation.factor = c(2, 2),
                      kfolds = k, bin.output = FALSE, clamp = TRUE,
                      rasterPreds = T, parallel = T, numCores = 18,algorithm="maxent.jar")
  save(modelo,file=paste('ENM_eval_jakk_R',i,sep="_"))
}


###pick model based on 

##OR then AUC

results_ordered<-modelo@results[order(modelo@results$avg.test.orMTP,-modelo@results$avg.test.AUC),]
best_modelOR<-as.integer(rownames(results_ordered[1,]))
model_predOR<-predict(modelo@models[[best_modelOR]],predictors1,args=c("outputformat=cloglog"))
plot(model_predOR)
points(occ1,pch=16,cex=0.6)
##AIC
best_model<-which (modelo@results$delta.AICc == 0)
best_model<-best_model[1]
model_pred<-predict(modelo@models[[best_model]],predictors1,args=c("outputformat=cloglog"))
plot(model_pred)
points(occ1,pch=16,cex=0.6)
    ###Write model to file
writeRaster(model_pred,paste(i,"AIC.asc",sep="_"),format="ascii")
writeRaster(model_predOR,paste(i,"OR_AUC.asc",sep="_"),format="ascii")


####Thresholding with T10
thresholds<-modelo@models[[best_modelOR]]@results
t10<-thresholds['X10.percentile.training.presence.Cloglog.threshold',]
reclass_mat_t10<-matrix(c(0,t10,0,t10,1,1),ncol=3,byrow=T)
pred_thresholded_t10<-reclassify(model_predOR,reclass_mat_t10)
    ###Write thresholded model to file 
writeRaster(pred_thresholded_t10,paste(i,"OR_AUC_T10.asc",sep="_"),format="ascii")

####Thresholding with MPT
thresholds<-modelo@models[[best_modelOR]]@results
mpt<-thresholds['Minimum.training.presence.Cloglog.threshold',]
reclass_mat_mpt<-matrix(c(0,mpt,0,mpt,1,1),ncol=3,byrow=T)
pred_thresholded_mpt<-reclassify(model_predOR,reclass_mat_mpt)
###Write thresholded model to file 
writeRaster(pred_thresholded_mpt,paste(i,"OR_AUC_MPT.asc",sep="_"),format="ascii")

  },error=function(e){cat("ERROR:",i,conditionMessage(e),"\\n")})
}
else{
  print (paste(i,"cannot be modeled less than 5 points"))
  bla<-paste(i,"cannot be modeled less than 5 points")
  write.table(bla,"errors.txt",append=T)
  }
}
pdf(paste("all_regions","curved.pdf", sep="_"))
response(modelo@models[[best_modelOR]])
dev.off()



