############################### Adam et al., 2021 - Diminishing potential for tropical reefs to function as coral diversity strongholds under climate change conditions #####################
#### R script describes calibration and validation of the coral species distribution models. Furthermore, the script details the construction of the final model continuous and binary habitat suitability predictions under present-day and future environmental conditions (example here only shows models created for Mean-type models)  #######

# Used packages #
library("raster")
library("rJava")
library("dismo")
library("ggplot2")
library("virtualspecies")
library("biomod2")
library("ggplot2")

# MODEL CALIBRATION #
# INPUT DATA #
# ENVIRONMENTAL DATA #
#import present-day mean environmental rasters (after removing variables with Pearson correlation >|0.85|; #cor<- cor(getValues(Envstack_present_selected), use = "pairwise.complete.obs", method = "pearson"))
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Present/Present_240620/Variables_models_krigged_250m_40/")
sstmax<-raster("sstmax_krigged_snap_250m_40.tif")
sstrange<-raster("sstrange_krigged_snap_250m_40.tif")
ssta<-raster("ssta_krigged_snap_250m_40.tif")
tsm<-raster("tsm_krigged_snap_250m_40.tif")
light<-raster("light_krigged_snap_250m_40.tif")
bath<-raster("bathymetry_snap_250m_40.tif")
roughness<-raster("roughness_snap_250m_40.tif")

Envstack_present_selected<-stack(sstmax,
               sstrange,
               ssta,
               tsm,
               light,
               bath,
               roughness)

#set coordinate system
WGS84<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
crs(Envstack_present_selected) <-WGS84

#import shapefile of Western Australia
setwd("/Users/19394297/Desktop/Coral_SDM_050519/")
WA <- readRDS("WAmap_WGS84.R")

#import future SSTmax/SSTrange and stack these rasters with other variables
#RCP2.6_2050
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Future/Future_CCSM4_MIROC_Hadgem2/250_tif/Future_variables_models_krigged_250m_40/")
SSTmax_2.6_2050<- raster("sstmax_26_2050_40_snap_krig_250m_40.tif")
#make sure names for future variables are the same as the variables used during model calibration (requirement for MaxEnt predictions)
names(SSTmax_2.6_2050)<-names(sstmax)
SSTrange_2.6_2050 <-raster("sstrange_ymax_ymin_26_2050_40_snap_krig_250m_40.tif")
names(SSTrange_2.6_2050)<-names(sstrange)
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Present/Present_240620/Variables_models_krigged_250m_40/")
ssta<-raster("ssta_krigged_snap_250m_40.tif")
tsm<-raster("tsm_krigged_snap_250m_40.tif")
light<-raster("light_krigged_snap_250m_40.tif")
bath<-raster("bathymetry_snap_250m_40.tif")
roughness<-raster("roughness_snap_250m_40.tif")

Envstack_future_RCP2.6_2050_selected <- stack(SSTmax_2.6_2050,SSTrange_2.6_2050,ssta,tsm,light,bath,roughness)
crs(Envstack_future_RCP2.6_2050_selected) <-WGS84

#RCP8.5_2050
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Future/Future_CCSM4_MIROC_Hadgem2/250_tif/Future_variables_models_krigged_250m_40/")
SSTmax_8.5_2050<- raster("sstmax_85_2050_40_snap_krig_250m_40.tif")
names(SSTmax_8.5_2050)<-names(sstmax)
SSTrange_8.5_2050 <-raster("sstrange_ymax_ymin_85_2050_40_snap_krig_250m_40.tif")
names(SSTrange_8.5_2050)<-names(sstrange)
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Present/Present_240620/Variables_models_krigged_250m_40/")
ssta<-raster("ssta_krigged_snap_250m_40.tif")
tsm<-raster("tsm_krigged_snap_250m_40.tif")
light<-raster("light_krigged_snap_250m_40.tif")
bath<-raster("bathymetry_snap_250m_40.tif")
roughness<-raster("roughness_snap_250m_40.tif")

Envstack_future_RCP8.5_2050_selected <- stack(SSTmax_8.5_2050,SSTrange_8.5_2050,ssta,tsm,light,bath,roughness)
crs(Envstack_future_RCP8.5_2050_selected) <-WGS84

#RCP2.6_2100
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Future/Future_CCSM4_MIROC_Hadgem2/250_tif/Future_variables_models_krigged_250m_40/")
SSTmax_2.6_2100<- raster("sstmax_26_2100_40_snap_krig_250m_40.tif")
names(SSTmax_2.6_2100)<-names(sstmax)
SSTrange_2.6_2100 <-raster("sstrange_ymax_ymin_26_2100_40_snap_krig_250m_40.tif")
names(SSTrange_2.6_2100)<-names(sstrange)
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Present/Present_240620/Variables_models_krigged_250m_40/")
ssta<-raster("ssta_krigged_snap_250m_40.tif")
tsm<-raster("tsm_krigged_snap_250m_40.tif")
light<-raster("light_krigged_snap_250m_40.tif")
bath<-raster("bathymetry_snap_250m_40.tif")
roughness<-raster("roughness_snap_250m_40.tif")

Envstack_future_RCP2.6_2100_selected <- stack(SSTmax_2.6_2100,SSTrange_2.6_2100,ssta,tsm,light,bath,roughness)
crs(Envstack_future_RCP2.6_2100_selected) <-WGS84

#RCP8.5_2100
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Future/Future_CCSM4_MIROC_Hadgem2/250_tif/Future_variables_models_krigged_250m_40/")
SSTmax_8.5_2100<- raster("sstmax_85_2100_40_snap_krig_250m_40.tif")
names(SSTmax_8.5_2100)<-names(sstmax)
SSTrange_8.5_2100 <-raster("sstrange_ymax_ymin_85_2100_40_snap_krig_250m_40.tif")
names(SSTrange_8.5_2100)<-names(sstrange)
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Variables_160320/Present/Present_240620/Variables_models_krigged_250m_40/")
ssta<-raster("ssta_krigged_snap_250m_40.tif")
tsm<-raster("tsm_krigged_snap_250m_40.tif")
light<-raster("light_krigged_snap_250m_40.tif")
bath<-raster("bathymetry_snap_250m_40.tif")
roughness<-raster("roughness_snap_250m_40.tif")

Envstack_future_RCP8.5_2100_selected <- stack(SSTmax_8.5_2100, SSTrange_8.5_2100, ssta,tsm,light,bath,roughness)
crs(Envstack_future_RCP8.5_2100_selected) <-WGS84

# CORAL DATA #
setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Coral_occurrence_060720")
#import filtered coral occurrence records (at least containing 3 columns: species names, Longitude and Latitude)
CoralWA_occdataframe <- read.csv("Coral_Occurrencedata_noNA_noduplication_250m_edited_060720.csv")
nrow(CoralWA_occdataframe)
CoralWA_occspatial<- CoralWA_occdataframe
# create spatial object
coordinates(CoralWA_occspatial) <- ~Longitude + Latitude
WGS84<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
crs(CoralWA_occspatial) <- WGS84
#file containing names of species with >=10 occurrence records
freq_Mostabund <- read.csv("Coral_until_freq10_060720.csv")

############################################
############## Loop to calibrate and validate coral SDMs and calculate present-day and future habitat suitability predictions for all 205 coral species 
############## (example here showing mean-type model construction) 
############################################

for (i in 1:205){
  setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Species_models_raw_MacbookPro_060720/")
  freq_Mostabund$Species <-as.character(freq_Mostabund$Species)
  #select species
  sp <-unique(freq_Mostabund$Species)[i]
  #make new folder to save all files related to the modelling
  dir.create(paste(sub(" ","_",sp),sep=""))
  setwd(paste(sub(" ","_",sp),sep=""))
  #subset species specific occurrence data --> spatial object 
  species <- subset(CoralWA_occspatial,Species == sp)
  #subset species specific occurrence data --> dataframe (for plotting)
  species_dataframe <-subset(CoralWA_occdataframe,Species == sp)
  
  # SELECT BACKGROUND POINTS #
  Speciespresence <- buffer(species, width=9000) #range to exclude background points selection from to avoid overfitting of the models
  clim_mask_background <- mask(Envstack_present_selected, Speciespresence, inverse=TRUE) #exclude buffer area in study area
  set.seed(1)
  p <-species
  numbbackgr<- 1000 #Number of background points to select
  
  #select background locations
  bg <- randomPoints(clim_mask_background,
                     p,
                     n=numbbackgr,
                     tryf=100,
                     excludep=T)
  
  # DATA SPLITTING TRAINING (75%) AND INDEPENDENT TEST DATA (25%) (for both occurrence records and background points) #
  set.seed(1)
  selected <- sample( 1:nrow(species),  nrow(species)*0.75)
  set.seed(1)
  selected_backgroundpoints <- sample(  1:nrow(bg),  nrow(bg)*0.75)
  occ_train <- species[selected,] # this is the selection of occurrence records used for model training
  saveRDS(occ_train,file=paste(sub(" ","_",sp),"_occ_train.R",sep=""))
  bg_train <- bg [selected_backgroundpoints,] # this is the selection of background points used for model training
  saveRDS(bg_train,file=paste(sub(" ","_",sp),"_bg_train.R",sep=""))
  occ_test <- species[-selected,] # selection of occurrence records which will be used for model testing
  saveRDS(occ_test,file=paste(sub(" ","_",sp),"_occ_test.R",sep=""))
  bg_test <- bg[-selected_backgroundpoints,] # selection of background points which will be used for model testing
  saveRDS(bg_test,file=paste(sub(" ","_",sp),"_bg_test.R",sep=""))

  #as dataframe to extract variable data at these locations (used for checking)
  occ_train_dataframe <-as.data.frame(occ_train)
  occ_test_dataframe <-as.data.frame(occ_test)
  bg_train_dataframe<-as.data.frame(bg_train)
  bg_test_dataframe<-as.data.frame(bg_test)
  
  Totalocc_train<- nrow(occ_train_dataframe)
  Totalocc_test<- nrow(occ_test_dataframe)
  Totalback_train <- nrow(bg_train_dataframe)
  Totalback_test <-nrow(bg_test_dataframe)
  Totalocc <- nrow(species)
  Speciesoccurrence<- cbind(Totalocc,Totalocc_train,Totalocc_test,Totalback_train,Totalback_test)
  
  write.csv(Speciesoccurrence, file=paste(sub(" ","_",sp),"_Occurrencedata.csv", sep=""),row.names = FALSE)
  write.csv(occ_train_dataframe, file=paste(sub(" ","_",sp),"_occ_train_data.csv", sep=""),row.names = FALSE)
  write.csv(occ_test_dataframe, file=paste(sub(" ","_",sp),"_occ_test_data.csv", sep=""),row.names = FALSE)
  write.csv(bg_train_dataframe, file=paste(sub(" ","_",sp),"_bg_train_data.csv", sep=""),row.names = FALSE)
  write.csv(bg_test_dataframe, file=paste(sub(" ","_",sp),"_bg_test_data.csv", sep=""),row.names = FALSE)
  
  #Extract environmental values at background training points and occurrence training locations
  data_occurrence <- data.frame(coordinates(occ_train), extract(Envstack_present_selected,occ_train))
  colnames(data_occurrence) [1] <- "longitude"
  colnames(data_occurrence) [2] <- "latitude"
  sp2<- gsub('\\\\s+', '', sp)
  Species_occurrence_final <- cbind(species =sp2 , data_occurrence)
  nrow(Species_occurrence_final)
  write.csv(Species_occurrence_final, file=paste(sub(" ","_",sp),"_occurrence_final_060720.csv",sep=""),row.names = FALSE)
  
  data_background <- data.frame(coordinates(bg_train), extract(Envstack_present_selected,bg_train))
  colnames(data_background) [1] <- "longitude"
  colnames(data_background) [2] <- "latitude"
  sp2<- gsub('\\\\s+', '', sp)
  Species_background_final <- cbind(species =sp2 , data_background)
  write.csv(Species_background_final, file=paste(sub(" ","_",sp),"_background_final_060720.csv",sep=""),row.names = FALSE)
  
  #Plot Modelling data (Training, testing data across the coastline of WA)
  ggplot() + geom_polygon(data=WA, aes(y=lat, x=long, group=group),
                          fill="grey", color="black")+coord_map() +coord_map()+geom_point(data=bg_train_dataframe, aes(x=x, y=y),
                                                                                          color="red",size = 0.5)+ geom_point(data=bg_test_dataframe, aes(x=x, y=y),
                                                                                                                              color="orange",size = 0.5)+ geom_point(data=occ_train_dataframe, 
                                                                                                                                                                     aes(x=Longitude, y=Latitude),color = "blue", size = 1)+ geom_point(data=occ_test_dataframe, aes(x=Longitude, y=Latitude)
                                                                                                                                                                                                                                        
                                                                                                                                                                                                                                        ,color = "green", size = 1) + labs(title= paste(sub(" ","_",sp),"_modelling",sep=""))
  ggsave(paste(sub(" ","_",sp),"_Model_Area.pdf",sep=""))
  
  # MODEL CALIBRATION #
  #Use of Maxentfunction in package 'dismo'
  #extracting present-day values at train data locations from all variables
  env_all_occ_train <- extract(Envstack_present_selected,occ_train)
  env_bg_train <- extract(Envstack_present_selected,bg_train)  
  
  #combine the conditions by row
  myPredictors_allvar <- rbind(env_all_occ_train,env_bg_train)
  
  #change matrix to dataframe
  myPredictors_allvar <- as.data.frame(myPredictors_allvar)
  
  #Maxent reads a 1 as presence occurrence records and 0 as background points
  myResponse_all <- c(rep(1,nrow(env_all_occ_train)),
                      rep(0,nrow(env_bg_train))) 
  #save initial models in seperate folder
  dir.create("Maxent_output_allvar")
  setwd("Maxent_output_allvar")
  
  #calibrate model using 5-fold crossvalidation
  Training_allvar_maxent <- dismo::maxent(x=myPredictors_allvar, 
                                 ## variable conditions
                                 p=myResponse_all,
                                 path=getwd(),
                                 ## this is the folder you will find maxent output
                                 args=c("replicates=5","replicatetype=crossvalidate","writebackgroundpredictions=TRUE",
                                        "outputformat=raw"
                                        #crossvalidation=5 replicates
                                 ))
  
  setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Species_models_raw_MacbookPro_060720/")
  setwd(paste(sub(" ","_",sp),sep=""))
  saveRDS(Training_allvar_maxent,file=paste(sub(" ","_",sp),"_Maxent_Crossvalidation_allvar.R",sep=""))
  Training_allvar_maxent_results <- Training_allvar_maxent@results
  Training_allvar_maxent_results_dataframe <- as.data.frame(Training_allvar_maxent_results)
  
  #select AUC test after crossvalidation of all replicates and average model
  AUC_test_repl1<- Training_allvar_maxent_results_dataframe$species_0[[8]] 
  AUC_test_repl2<- Training_allvar_maxent_results_dataframe$species_1[[8]]
  AUC_test_repl3<- Training_allvar_maxent_results_dataframe$species_2[[8]]
  AUC_test_repl4<- Training_allvar_maxent_results_dataframe$species_3[[8]]
  AUC_test_repl5<- Training_allvar_maxent_results_dataframe$species_4[[8]]
  AUC_test_average<- Training_allvar_maxent$`species (average)`[[8]]
  
  AUC_test_replicates<- cbind(AUC_test_repl1,AUC_test_repl2,AUC_test_repl3,AUC_test_repl4,AUC_test_repl5,AUC_test_average)
  write.csv(AUC_test_replicates, file=paste(sub(" ","_",sp),"_AUC_test_replicates.csv",sep=""),row.names = F)
  
  #extract permutation importance for all variables of all replicates and average model
  a<- cbind(Training_allvar_maxent$species_0[[18]],Training_allvar_maxent$species_1[[18]],Training_allvar_maxent$species_2[[18]],Training_allvar_maxent$species_3[[18]],Training_allvar_maxent$species_4[[18]],Training_allvar_maxent$`species (average)`[[18]])
  rownames(a)<-rownames(Training_allvar_maxent[18,])
  colnames(a)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  b<- cbind(Training_allvar_maxent$species_0[[19]],Training_allvar_maxent$species_1[[19]],Training_allvar_maxent$species_2[[19]],Training_allvar_maxent$species_3[[19]],Training_allvar_maxent$species_4[[19]],Training_allvar_maxent$`species (average)`[[19]])
  rownames(b)<-rownames(Training_allvar_maxent[19,])
  colnames(b)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  c<- cbind(Training_allvar_maxent$species_0[[20]],Training_allvar_maxent$species_1[[20]],Training_allvar_maxent$species_2[[20]],Training_allvar_maxent$species_3[[20]],Training_allvar_maxent$species_4[[20]],Training_allvar_maxent$`species (average)`[[20]])
  rownames(c)<-rownames(Training_allvar_maxent[20,])
  colnames(c)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  d<- cbind(Training_allvar_maxent$species_0[[21]],Training_allvar_maxent$species_1[[21]],Training_allvar_maxent$species_2[[21]],Training_allvar_maxent$species_3[[21]],Training_allvar_maxent$species_4[[21]],Training_allvar_maxent$`species (average)`[[21]])
  rownames(d)<-rownames(Training_allvar_maxent[21,])
  colnames(d)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  e<- cbind(Training_allvar_maxent$species_0[[22]],Training_allvar_maxent$species_1[[22]],Training_allvar_maxent$species_2[[22]],Training_allvar_maxent$species_3[[22]],Training_allvar_maxent$species_4[[22]],Training_allvar_maxent$`species (average)`[[22]])
  rownames(e)<-rownames(Training_allvar_maxent[22,])
  colnames(e)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  f<- cbind(Training_allvar_maxent$species_0[[23]],Training_allvar_maxent$species_1[[23]],Training_allvar_maxent$species_2[[23]],Training_allvar_maxent$species_3[[23]],Training_allvar_maxent$species_4[[23]],Training_allvar_maxent$`species (average)`[[23]])
  rownames(f)<-rownames(Training_allvar_maxent[23,])
  colnames(f)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  g<- cbind(Training_allvar_maxent$species_0[[24]],Training_allvar_maxent$species_1[[24]],Training_allvar_maxent$species_2[[24]],Training_allvar_maxent$species_3[[24]],Training_allvar_maxent$species_4[[24]],Training_allvar_maxent$`species (average)`[[24]])
  rownames(g)<-rownames(Training_allvar_maxent[24,])
  colnames(g)<-c("replicate1","replicate2","replicate3","replicate4","replicate5","average")
  
  Variables_replicates<- rbind(a,b,c,d,e,f,g)
  write.csv(Variables_replicates, file=paste(sub(" ","_",sp),"_Variables_replicates.csv",sep=""),row.names = T)
 
  Variables_replicates<-as.data.frame(Variables_replicates)
  Variables_replicates <- cbind(Row.Names = rownames(Variables_replicates), Variables_replicates)
  Variables_replicates$Row.Names<-as.character(Variables_replicates$Row.Names)
  Variables_replicates$Variable<-sub(".per\\\\S*", "", Variables_replicates$Row.Names) #remove .perm..... name
  Variables_replicates$average<-as.numeric(Variables_replicates$average)
  
  #select only variables that have permutation.importance >=1% in average model
  Subset_variables <-Variables_replicates$Variable[which(Variables_replicates$average>=1)]
  
  #subset only variables in the subset in present and future variable stacks used to build most parsimonious model
  presentsub <- subset(Envstack_present_selected,Subset_variables)
  RCP2.6_2050_sub <-subset(Envstack_future_RCP2.6_2050_selected,Subset_variables)
  RCP8.5_2050_sub <-subset(Envstack_future_RCP8.5_2050_selected,Subset_variables)
  RCP2.6_2100_sub <-subset(Envstack_future_RCP2.6_2100_selected,Subset_variables)
  RCP8.5_2100_sub <-subset(Envstack_future_RCP8.5_2100_selected,Subset_variables)

  #Remodel only using all training data and subset of variables --> will be used for model validation (parsimonious model)
  #extracting present-day values at train data locations from subset of variables
  env_occ_train_PAR <- extract(presentsub,occ_train)
  env_bg_train_PAR <- extract(presentsub,bg_train)  
  
  #combine the conditions by row
  myPredictors_PAR <- rbind(env_occ_train_PAR,env_bg_train_PAR)
  
  #change matrix to dataframe
  myPredictors_PAR <- as.data.frame(myPredictors_PAR)
  
  myResponse_PAR <- c(rep(1,nrow(env_occ_train_PAR)),
                      rep(0,nrow(env_bg_train_PAR))) 
  
  #create directory to save parimonious model
  dir.create("Maxent_output_PAR")
  setwd("Maxent_output_PAR")
  
  #Build parsimonious model used for model validation
  Training_PAR_maxent <- dismo::maxent(x=myPredictors_PAR, 
                                     p=myResponse_PAR,
                                     path=getwd(),
                                     ## this is the folder you will find maxent output
                                     args=c("writebackgroundpredictions=TRUE"
                                     ))
  
  setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Species_models_raw_MacbookPro_060720/")
  setwd(paste(sub(" ","_",sp),sep=""))
  saveRDS(Training_PAR_maxent,file=paste(sub(" ","_",sp),"_Maxent_PAR.R",sep=""))

  #variable importance of parsimonious model
  variableimp_Training_PAR_maxent<- var.importance(Training_PAR_maxent)
  write.csv(variableimp_Training_PAR_maxent, file=paste(sub(" ","_",sp),"_variableimp_PAR.csv",sep=""),row.names = F)
  #response curve of parsimonious model
  response(Training_PAR_maxent)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Responsecurves_allvar.pdf", sep=""), width = 7, height = 5)
  
  # MODEL VALIDATION #
  #create present-day prediction using calibrated parsimonious model
  PAR_pred_present_species <- predict(Training_PAR_maxent,presentsub)
  plot(PAR_pred_present_species)
  writeRaster(PAR_pred_present_species, file=paste(sub(" ","_",sp),"_Presdistr_training_PAR.tif",sep=""), format="GTiff")
  
  #extract coordinates
  species_training_coord <- occ_train@coords
  species_test_coord<- occ_test@coords
  bg_coord <- bg
  bg_train_coord <- bg_train
  bg_test_coord <- bg_test
  
  #extract habitat suitability predictions from parsimonious model
  occ_suit <-extract(PAR_pred_present_species,as.data.frame(species_training_coord))
  test_suit <-extract(PAR_pred_present_species,as.data.frame(species_test_coord)) ### prediction values on withheld test data coordinates
  bg_train_suit<- extract(PAR_pred_present_species,bg_train_coord)
  bg_test_suit<- extract(PAR_pred_present_species,bg_test_coord)
  eval_suit<- dismo::evaluate(occ_suit,bg_train_suit)#training evaluation
  eval_test_suit<- dismo::evaluate(test_suit,bg_test_suit) #testing evaluation
  
  #calculate threshold independent metric used for model evaluation, area under the curve,AUC on withheld test data (threshold >0.7)
  AUC_train <- eval_suit@auc #AUC on training data
  AUC_test <- eval_test_suit@auc #AUC on withheld test data

  #use maximum sensitivity specificity threshold for converting continuous habitat suitability predictions to binary presence/absence predictions
  threshold_species <-threshold(eval_suit)
  threshold_binary<-threshold_species$"spec_sens"
  
  #calculate threshold dependent metrics used for model evaluation, sensitivity on withheld test data (threshold >0)
  sum(test_suit > threshold_binary) -> majortest
  sum(test_suit < threshold_binary) -> minortest

  sensitivity_test <- (majortest) / (majortest+minortest) 
  
  Evaluation <- cbind(Totalocc,threshold_species,threshold_binary,AUC_train,AUC_test,sensitivity_test)
  write.csv(Evaluation, file=paste(sub(" ","_",sp),"_Evalmetrics_PAR.csv",sep=""),row.names = F)
  
  # CONSTRUCTION FINAL MODEL --> using all available species occurrence and background data as model evaluation has been finalised #
  env_occ_all <- extract(presentsub,species)
  env_bg_all <- extract(presentsub,bg)  
  
  #combine the conditions by row
  myPredictors_final <- rbind(env_occ_all,env_bg_all)
  
  #change matrix to dataframe
  myPredictors_final <- as.data.frame(myPredictors_final)
  
  myResponse_final <- c(rep(1,nrow(env_occ_all)),
                        rep(0,nrow(env_bg_all))) 
  
  dir.create("Maxent_finalmodel")
  setwd("Maxent_finalmodel")
  Species_maxent_Final <- dismo::maxent(x=myPredictors_final, 
                                       ## env conditions, here we selected only 3 predictors
                                       p=myResponse_final,
                                       path=getwd(),
                                       ## this is the folder you will find maxent output
                                       args=c("writebackgroundpredictions=TRUE"
                                       ))
  
  setwd("/Volumes/Backup_Plus/Coral_SDM_050519/Species_models_raw_MacbookPro_060720/")
  setwd(paste(sub(" ","_",sp),sep=""))

  saveRDS(Species_maxent_Final,file=paste(sub(" ","_",sp),"_ENM_FinalModel.R",sep=""))
  
  variableimp_finalmodel<- var.importance(Species_maxent_Final)
  write.csv(variableimp_finalmodel, file=paste(sub(" ","_",sp),"_variableimp_finalmodel.csv",sep=""),row.names = F)
  response(Species_maxent_Final)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Responsecurves_finalmodel.pdf", sep=""), width = 7, height = 5)
  
  # MODEL PREDICTIONS #
  #Presence predictions
  species_present_prediction <-predict(Species_maxent_Final,presentsub)
  plot(species_present_prediction)
  #save plot and raster file
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_present_prediction.pdf", sep=""), width = 7, height = 5)
  writeRaster(species_present_prediction, file=paste(sub(" ","_",sp),"_Present_prediction.tif",sep=""), overwrite=T,format="GTiff")
  
  #use maximum sens_spec threshold to convert continuous predictions into binary presence/absence predictions
  binary_present_species<-species_present_prediction>threshold_binary
  plot(binary_present_species)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_present_binary.pdf", sep=""), width = 7, height = 5)
  writeRaster(binary_present_species, file=paste(sub(" ","_",sp),"_Present_binary.tif",sep=""), overwrite=T,format="GTiff")
  
  #Future prediction RCP2.6 - 2050
  species_RCP2.6_2050_prediction <-predict(Species_maxent_Final,RCP2.6_2050_sub) #make sure that the future layers have the same names as present raster files
  plot(species_RCP2.6_2050_prediction)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP2.6_2050_prediction.pdf", sep=""), width = 7, height = 5)
  writeRaster(species_RCP2.6_2050_prediction, file=paste(sub(" ","_",sp),"_RCP2.6_2050_prediction.tif",sep=""), overwrite=T,format="GTiff")
  
  binary_RCP2.6_2050_species<-species_RCP2.6_2050_prediction>threshold_binary
  plot(binary_RCP2.6_2050_species)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP2.6_2050_binary.pdf", sep=""), width = 7, height = 5)
  writeRaster(binary_RCP2.6_2050_species, file=paste(sub(" ","_",sp),"_RCP2.6_2050_binary.tif",sep=""), overwrite=T,format="GTiff")
  
  #Future prediction RCP8.5 - 2050
  species_RCP8.5_2050_prediction <-predict(Species_maxent_Final,RCP8.5_2050_sub)#make sure that the future layers have the same names as ascii files
  plot(species_RCP8.5_2050_prediction)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP8.5_2050_prediction.pdf", sep=""), width = 7, height = 5)
  writeRaster( species_RCP8.5_2050_prediction, file=paste(sub(" ","_",sp),"_RCP8.5_2050_prediction.tif",sep=""), overwrite=T,format="GTiff")
  
  binary_RCP8.5_2050_species<-species_RCP8.5_2050_prediction>threshold_binary
  plot(binary_RCP8.5_2050_species)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP8.5_2050_binary.pdf", sep=""), width = 7, height = 5)
  
  writeRaster(binary_RCP8.5_2050_species, file=paste(sub(" ","_",sp),"_RCP8.5_2050_binary.tif",sep=""), overwrite=T,format="GTiff")
  
  ###Future prediction RCP2.6 - 2100
  species_RCP2.6_2100_prediction <-predict(Species_maxent_Final,RCP2.6_2100_sub)#make sure that the future layers have the same names as ascii files
  plot(species_RCP2.6_2100_prediction)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP2.6_2100_prediction.pdf", sep=""), width = 7, height = 5)
  writeRaster( species_RCP2.6_2100_prediction, file=paste(sub(" ","_",sp),"_RCP2.6_2100_prediction.tif",sep=""), overwrite=T,format="GTiff")
  
  binary_RCP2.6_2100_species<-species_RCP2.6_2100_prediction>threshold_binary
  plot(binary_RCP2.6_2100_species)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP2.6_2100_binary.pdf", sep=""), width = 7, height = 5)
  
  writeRaster(binary_RCP2.6_2100_species, file=paste(sub(" ","_",sp),"_RCP2.6_2100_binary.tif",sep=""), overwrite=T,format="GTiff")
  
  #Future prediction RCP8.5 - 2100
  species_RCP8.5_2100_prediction <-predict(Species_maxent_Final,RCP8.5_2100_sub)#make sure that the future layers have the same names as ascii files
  plot(species_RCP8.5_2100_prediction)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP8.5_2100_prediction.pdf", sep=""), width = 7, height = 5)
  writeRaster( species_RCP8.5_2100_prediction, file=paste(sub(" ","_",sp),"_RCP8.5_2100_prediction.tif",sep=""), overwrite=T,format="GTiff")
  
  binary_RCP8.5_2100_species<-species_RCP8.5_2100_prediction>threshold_binary
  plot(binary_RCP8.5_2100_species)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_RCP8.5_2100_binary.pdf", sep=""), width = 7, height = 5)
  
  writeRaster(binary_RCP8.5_2100_species, file=paste(sub(" ","_",sp),"_RCP8.5_2100_binary.tif",sep=""), overwrite=T,format="GTiff")
  
  #Calculate range shifts between present-day and future climate conditions using package 'biomod2'
  speciesRangeSize_RCP2.6_2050 <- BIOMOD_RangeSize(CurrentPred=binary_present_species,FutureProj=binary_RCP2.6_2050_species)
  speciesRangeSize_RCP8.5_2050 <- BIOMOD_RangeSize(CurrentPred=binary_present_species,FutureProj=binary_RCP8.5_2050_species)
  speciesRangeSize_RCP2.6_2100 <- BIOMOD_RangeSize(CurrentPred=binary_present_species,FutureProj=binary_RCP2.6_2100_species)
  speciesRangeSize_RCP8.5_2100 <- BIOMOD_RangeSize(CurrentPred=binary_present_species,FutureProj=binary_RCP8.5_2100_species)
  
  Rangesize_metrics_RCP2.6_2050<- speciesRangeSize_RCP2.6_2050$Compt.By.Models
  Rangesize_metrics_RCP8.5_2050<- speciesRangeSize_RCP8.5_2050$Compt.By.Models
  Rangesize_metrics_RCP2.6_2100<- speciesRangeSize_RCP2.6_2100$Compt.By.Models
  Rangesize_metrics_RCP8.5_2100<- speciesRangeSize_RCP8.5_2100$Compt.By.Models
  
  write.csv(Rangesize_metrics_RCP2.6_2050, file=paste(sub(" ","_",sp),"_Rangesize_metrics_RCP2.6_2050.csv",sep=""))
  write.csv(Rangesize_metrics_RCP8.5_2050, file=paste(sub(" ","_",sp),"_Rangesize_metrics_RCP8.5_2050.csv",sep=""))
  write.csv(Rangesize_metrics_RCP2.6_2100, file=paste(sub(" ","_",sp),"_Rangesize_metrics_RCP2.6_2100.csv",sep=""))
  write.csv(Rangesize_metrics_RCP8.5_2100, file=paste(sub(" ","_",sp),"_Rangesize_metrics_RCP8.5_2100.csv",sep=""))
  
  #Save range shift predictions (plot as pdf and raster as Geotiff/ascii)
  # -2 --> if the given pixel is predicted to be lost by the species
  #-1 if the given pixel is predicted to be stable for the species
  #0 if the given pixel was not occupied and will not be in the future
  #1 if the given pixel was not occupied and is predicted to be into the future
  
  Rangesize_raster_RCP2.6_2050<- speciesRangeSize_RCP2.6_2050$Diff.By.Pixel$layer
  plot(Rangesize_raster_RCP2.6_2050)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Rangesize_RCP2.6_2050.pdf", sep=""), width = 7, height = 5)
  
  Rangesize_raster_RCP8.5_2050<- speciesRangeSize_RCP8.5_2050$Diff.By.Pixel$layer
  plot(Rangesize_raster_RCP8.5_2050)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Rangesize_RCP8.5_2050.pdf", sep=""), width = 7, height = 5)
  
  Rangesize_raster_RCP2.6_2100<- speciesRangeSize_RCP2.6_2100$Diff.By.Pixel$layer
  plot(Rangesize_raster_RCP2.6_2100)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Rangesize_RCP2.6_2100.pdf", sep=""), width = 7, height = 5)
  
  Rangesize_raster_RCP8.5_2100<- speciesRangeSize_RCP8.5_2100$Diff.By.Pixel$layer
  plot(Rangesize_raster_RCP8.5_2100)
  dev.copy2pdf(file=paste(sub(" ","_",sp),"_Rangesize_RCP8.5_2100.pdf", sep=""), width = 7, height = 5)
  
  writeRaster(Rangesize_raster_RCP2.6_2050, file=paste(sub(" ","_",sp),"_Rangesize_RCP2.6_2050.tif",sep=""), overwrite=T,format="GTiff")
  writeRaster(Rangesize_raster_RCP8.5_2050, file=paste(sub(" ","_",sp),"_Rangesize_RCP8.5_2050.tif",sep=""), overwrite=T,format="GTiff")
  writeRaster(Rangesize_raster_RCP2.6_2100, file=paste(sub(" ","_",sp),"_Rangesize_RCP2.6_2100.tif",sep=""), overwrite=T,format="GTiff")
  writeRaster(Rangesize_raster_RCP8.5_2100, file=paste(sub(" ","_",sp),"_Rangesize_RCP8.5_2100.tif",sep=""), overwrite=T,format="GTiff")
  
}

# END #