#####-----Modeling the dispersal of Pikas-----#####
path<-"E:/Pika/WD_0517/"

####Import libraries####
library(MigClim)
library(rgdal)
library(ecospat)
####Species to model###
spe_list<-c("Ochotona cansus","Ochotona curzoniae","Ochotona dauurica","Ochotona erythrotis","Ochotona macrotis")

####Import IUCN polygons####
Pikas_IUCN <-readOGR(paste(path,"data/Occurrence/Pikas_IUCN/data_0.shp",sep=""))

rcp="rcp60"
setwd(paste("E:/Pika/WD_0517/dispersal_",rcp,sep=""))
for(i in 1:5){
  species=spe_list[i]
  mod.EM<-readRDS(paste(path,"models/",species,"_mod.EM.rds",sep=""))
  
  ##Introduce the binary converting threshold
  threshold<-ecospat.ESM.threshold(mod.EM)
  threshold<-threshold$TSS.th[4]*1000
  
  ####1. Prepare data####
  
  ###1.1 Current occupied habitats (InitialDist)
  #Current suitable habitats
  EFproj_current<-readRDS(paste(path,"EFprojections/",species,"_mod.EFproj.current.rds",sep=""))
  sh_current<-ecospat.binary.model(EFproj_current$EF, threshold)
  
  #Load IUCN polygon for this species
  spe_IUCN<-subset(Pikas_IUCN, BINOMIAL==species) #IUCN range of i species
  
  allsites<-as.data.frame(xyFromCell(sh_current, 1:ncell(sh_current)))
  sh_current<-cbind(allsites, values(sh_current))
  colnames(allsites)<-c("Longitude","Latitude")
  coordinates(allsites) <- ~ Longitude + Latitude
  proj4string(allsites) <- proj4string(Pikas_IUCN)
  sh_current$IUCN<-complete.cases(over(allsites,spe_IUCN))
  #Current occupied habitats
  oh_current<-sh_current
  oh_current$`values(sh_current)`[is.na(oh_current$`values(sh_current)`)]<-0
  oh_current[oh_current$`values(sh_current)`==1&oh_current$IUCN==F,3]<-0 
  
  #add a barrier column
  water<-raster(paste(path,"data/Permanent_inland_water/permanent_inland_water.tif",sep=""))
  water=crop(water,c(61,130,17,63))
  water<-data.frame(raster::extract(water,allsites, sp = T))
  water<-data.frame('permanent_inland_water'=water$permanent_inland_water)
  water$permanent_inland_water[is.na(water$permanent_inland_water)]<-0
  water$permanent_inland_water[water$permanent_inland_water>0]<-1
  oh_current<-cbind(oh_current,water)
  
  ###1.2 Habitat suitability in the future
  
  ##considering no change
  combined_nochange<-oh_current
  hsmap0<-ecospat.binary.model(EFproj_current$EF, threshold)
  hsmap0<-values(hsmap0)
  hsmap0[is.na(hsmap0)]<-0
  hsmap0[hsmap0==1]<-1000
  combined_nochange<-cbind(combined_nochange,hsmap0)
  
  for(scenario in c("wc","lu","wclu")){
    for(period in c("2021_2040","2041_2060","2061_2080")){
      EFproj_future<-readRDS(paste(path,"EFprojections/",species,"_",rcp,"_mod.EFproj_",period,"_",scenario,".rds",sep=""))
      sh_future<-ecospat.binary.model(EFproj_future$EF, threshold)
      sh_future<-values(sh_future)
      sh_future[is.na(sh_future)]<-0
      sh_future[sh_future==1]<-1000
      assign(paste("hsmap_",period,"_",scenario,sep=""),sh_future)
    }
  }
  combined_wc<-do.call(cbind,list(oh_current,hsmap_2021_2040_wc,hsmap_2041_2060_wc,hsmap_2061_2080_wc))
  colnames(combined_wc)<-c("X","Y","oh_current","IUCN","barrier","hsmap1","hsmap2","hsmap3")
  
  combined_lu<-do.call(cbind,list(oh_current,hsmap_2021_2040_lu,hsmap_2041_2060_lu,hsmap_2061_2080_lu))
  colnames(combined_lu)<-c("X","Y","oh_current","IUCN","barrier","hsmap1","hsmap2","hsmap3")
  
  combined_wclu<-do.call(cbind,list(oh_current,hsmap_2021_2040_wclu,hsmap_2041_2060_wclu,hsmap_2061_2080_wclu))
  colnames(combined_wclu)<-c("X","Y","oh_current","IUCN","barrier","hsmap1","hsmap2","hsmap3")
  
  
  ####2. Modeling dispersal####
  #plot(1:10,exp(-(1:10)/3),type="l",lwd=2,xlab="Distance (km)", ylab='Dispersal probability',ylim=c(0,1))
  dist<-1:10
  dneg<-exp(-dist/3)
  spe_name<-sub("Ochotona ","",species)
  
  MigClim.migrate(iniDist=combined_nochange[,1:3],
                  hsMap=combined_nochange[,6],
                  envChgSteps=1, 
                  dispSteps=60, 
                  dispKernel=dneg,
                  barrier=combined_nochange[,5],
                  iniMatAge=1,
                  propaguleProd=c(1.0),
                  simulName=paste(spe_name, "_nochange",sep=""),
                  replicateNb=3, 
                  overWrite=TRUE, 
                  testMode=FALSE, 
                  fullOutput=FALSE,
                  keepTempFiles=FALSE)
  
  MigClim.migrate(iniDist=combined_wc[,1:3],
                  hsMap=combined_wc[,6:8],
                  envChgSteps=3, 
                  dispSteps=20, 
                  dispKernel=dneg,
                  barrier=combined_wc[,5],
                  iniMatAge=1,
                  propaguleProd=c(1.0),
                  simulName=paste(spe_name, "_wc_",rcp,sep=""),
                  replicateNb=3, 
                  overWrite=TRUE, 
                  testMode=FALSE, 
                  fullOutput=FALSE,
                  keepTempFiles=FALSE)
  
  MigClim.migrate(iniDist=combined_lu[,1:3],
                  hsMap=combined_lu[,6:8],
                  envChgSteps=3, 
                  dispSteps=20, 
                  dispKernel=dneg,
                  barrier=combined_lu[,5],
                  iniMatAge=1,
                  propaguleProd=c(1.0),
                  simulName=paste(spe_name, "_lu_",rcp,sep=""),
                  replicateNb=3, 
                  overWrite=TRUE, 
                  testMode=FALSE, 
                  fullOutput=FALSE,
                  keepTempFiles=FALSE)
  
  MigClim.migrate(iniDist=combined_wclu[,1:3],
                  hsMap=combined_wclu[,6:8],
                  envChgSteps=3, 
                  dispSteps=20, 
                  dispKernel=dneg,
                  barrier=combined_wclu[,5],
                  iniMatAge=1,
                  propaguleProd=c(1.0),
                  simulName=paste(spe_name, "_wclu_",rcp,sep=""),
                  replicateNb=3, 
                  overWrite=TRUE, 
                  testMode=FALSE, 
                  fullOutput=FALSE,
                  keepTempFiles=FALSE)
}



