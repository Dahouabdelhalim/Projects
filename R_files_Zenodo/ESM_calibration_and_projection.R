#####-----Run Species Distribution Models for Each Pika Species-----#####
#####-----May 17th 2021-----#####

setwd("E:/Pika/WD_0517/")

####Import packages and functions####
library(rgdal)
library(biomod2)
library(raster)
library(usdm)
library(ecospat)
memory.limit(size=31500000)
rasterOptions(maxmemory=2.5e+10,tmpdir="F:/TempR/",tmptime=12)

####Import IUCN polygons####
Pikas_IUCN <-readOGR("E:/Pika/WD_0517/data/Occurrence/Pikas_IUCN/data_0.shp")

####Species to model###
spe_list<-c("Ochotona cansus","Ochotona curzoniae","Ochotona dauurica","Ochotona erythrotis","Ochotona macrotis")

####Import environmental data####
path="E:/Pika/WD_0517/data/Ensembled/"
env_1979_2013<-raster::stack(paste(path,"env_1979_2013.tif",sep="")) #raster stack contains all environmental variables
varnames<-read.csv("E:/Pika/WD_0517/data/Ensembled/varnames.csv",head=F)$V1 #names of environmental variables
names(env_1979_2013)<-varnames

####Load variable selection result####
env_1979_2013.NoCor<-readRDS(file = "env_1979_2013.NoCor.rds") #selected variables using VIF<2.5

#####select variables
env_1979_2013.Reduced <- exclude(env_1979_2013, env_1979_2013.NoCor)
env_1979_2013<-env_1979_2013.Reduced

####2. A forloop for running SDMs for each Pika species####
for(i in 1:length(spe_list)){
      species=spe_list[i]
      print(paste("Now runing ",species,sep=""))
      
      ####3.Prepare occurrence data####
      
      ###3.1 Load in all sources of Occurrence data
      ##Load IUCN polygon
      spe_IUCN<-subset(Pikas_IUCN, BINOMIAL==species) #IUCN range of i species
    
      ##Load Occurrence data (csv files)
      path_occ="E:/Pika/WD_0517/data/Occurrence/"
      occ<-as.data.frame(read.csv(paste(path_occ,species,".csv",sep=""),header=T))[,1:2]
      colnames(occ)<-c("Longitude","Latitude")
      coordinates(occ) <- ~ Longitude + Latitude
      proj4string(occ) <- proj4string(spe_IUCN)
    
      ####4.Formating the data with the BIOMOD_FormatingData() function####
      mod.dat<-BIOMOD_FormatingData(resp.var=occ,
                                    expl.var=env_1979_2013,
                                    resp.name=species,
                                    PA.nb.absences=1000,
                                    PA.strategy="random")
      
      Biomod.tuning<-BIOMOD_tuning(mod.dat,
                                   models=c('GLM','ANN','MAXENT.Phillips'),
                                   env.ME = env_1979_2013)
      saveRDS(Biomod.tuning,paste(species,"_Biomod_tuning.rds",sep=""))
      Biomod.tuning<-readRDS(file = paste("tuning/",species,"_Biomod_tuning.rds",sep=""))
      
      ####5.Calibration of simple bivariate models####
      mod<-ecospat.ESM.Modeling(mod.dat,
                                models=c('GLM','ANN','MAXENT.Phillips'),
                                models.options = Biomod.tuning$models.options,
                                NbRunEval=10,
                                DataSplit=80,
                                weighting.score=c('SomersD'))

      ####6.Evaluation and average of simple bivariate models to ESMs####
      mod.EM<-ecospat.ESM.EnsembleModeling(mod,weighting.score="SomersD")
      saveRDS(mod.EM, paste(species,"_mod.EM.rds",sep=""))
        
      EM_Eval<-mod.EM$ESM.evaluations
      write.csv(EM_Eval,paste(species,"_EM_Eval_rep10.csv",spe=""),row.names=F)
      
      EM_Contrib<-ecospat.ESM.VarContrib(mod,mod.EM)
      write.csv(EM_Contrib,paste(species,"_EM_Contrib.csv",spe=""),row.names=F)
    
      ####7.Make projections under current condition####
      mod.proj.current<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                               new.env=env_1979_2013)
      saveRDS(mod.proj.current, paste(species,"_mod.proj.current.rds",sep=""))
      
      mod.EFproj.current<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj.current,
                                                         ESM.EnsembleModeling.output=mod.EM)
      
      saveRDS(mod.EFproj.current, paste(species,"_mod.EFproj.current.rds",sep=""))
      
      ####9.Make projections under future condition####
      for(rcp in c("rcp45","rcp60","rcp85")){
          ####Import environmental data for future####
          path=paste("E:/Pika/WD_0517/data/Ensembled/",rcp,"/",sep="")
          env_2021_2040_wc<-raster::stack(paste(path,"env_2021-2040_wc.tif",sep=""))
          env_2041_2060_wc<-raster::stack(paste(path,"env_2041-2060_wc.tif",sep=""))
          env_2061_2080_wc<-raster::stack(paste(path,"env_2061-2080_wc.tif",sep=""))
          
          env_2021_2040_lu<-raster::stack(paste(path,"env_2021-2040_lu.tif",sep=""))
          env_2041_2060_lu<-raster::stack(paste(path,"env_2041-2060_lu.tif",sep=""))
          env_2061_2080_lu<-raster::stack(paste(path,"env_2061-2080_lu.tif",sep=""))
          
          env_2021_2040_wclu<-raster::stack(paste(path,"env_2021-2040_wclu.tif",sep=""))
          env_2041_2060_wclu<-raster::stack(paste(path,"env_2041-2060_wclu.tif",sep=""))
          env_2061_2080_wclu<-raster::stack(paste(path,"env_2061-2080_wclu.tif",sep=""))
          
          names(env_2021_2040_wc)<-varnames
          names(env_2041_2060_wc)<-varnames
          names(env_2061_2080_wc)<-varnames
          
          names(env_2021_2040_lu)<-varnames
          names(env_2041_2060_lu)<-varnames
          names(env_2061_2080_lu)<-varnames
          
          names(env_2021_2040_wclu)<-varnames
          names(env_2041_2060_wclu)<-varnames
          names(env_2061_2080_wclu)<-varnames
          
          env_2021_2040_wc<-exclude(env_2021_2040_wc, env_1979_2013.NoCor)
          env_2041_2060_wc<-exclude(env_2041_2060_wc, env_1979_2013.NoCor)
          env_2061_2080_wc<-exclude(env_2061_2080_wc, env_1979_2013.NoCor)
          
          env_2021_2040_lu<-exclude(env_2021_2040_lu, env_1979_2013.NoCor)
          env_2041_2060_lu<-exclude(env_2041_2060_lu, env_1979_2013.NoCor)
          env_2061_2080_lu<-exclude(env_2061_2080_lu, env_1979_2013.NoCor)
          
          env_2021_2040_wclu<-exclude(env_2021_2040_wclu, env_1979_2013.NoCor)
          env_2041_2060_wclu<-exclude(env_2041_2060_wclu, env_1979_2013.NoCor)
          env_2061_2080_wclu<-exclude(env_2061_2080_wclu, env_1979_2013.NoCor)
      
          for(period in c("2021_2040","2041_2060","2061_2080")){
            ###only climate change
            mod.proj_wc<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                                new.env=eval(parse(text = paste("env_",period,"_wc",sep=""))))
            saveRDS(mod.proj_wc, paste(species,"_",rcp,"_mod.proj_",period,"_wc.rds",sep=""))
            
            mod.EFproj_wc<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj_wc,
                                                                    ESM.EnsembleModeling.output=mod.EM)
            saveRDS(mod.EFproj_wc, paste(species,"_",rcp,"_mod.EFproj_",period,"_wc.rds",sep=""))
            
            ###only land-use change
            mod.proj_lu<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                                new.env=eval(parse(text = paste("env_",period,"_lu",sep=""))))
            saveRDS(mod.proj_lu, paste(species,"_",rcp,"_mod.proj_",period,"_lu.rds",sep=""))
            
            mod.EFproj_lu<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj_lu,
                                                          ESM.EnsembleModeling.output=mod.EM)
            saveRDS(mod.EFproj_lu, paste(species,"_",rcp,"_mod.EFproj_",period,"_lu.rds",sep=""))
            
            ###both climate-change and land-use change
            mod.proj_wclu<-ecospat.ESM.Projection(ESM.modeling.output=mod,
                                                new.env=eval(parse(text = paste("env_",period,"_wclu",sep=""))))
            saveRDS(mod.proj_wclu, paste(species,"_",rcp,"_mod.proj_",period,"_wclu.rds",sep=""))
            
            mod.EFproj_wclu<-ecospat.ESM.EnsembleProjection(ESM.prediction.output=mod.proj_wclu,
                                                          ESM.EnsembleModeling.output=mod.EM)
            saveRDS(mod.EFproj_wclu, paste(species,"_",rcp,"_mod.EFproj_",period,"_wclu.rds",sep=""))
          }
      }
}


