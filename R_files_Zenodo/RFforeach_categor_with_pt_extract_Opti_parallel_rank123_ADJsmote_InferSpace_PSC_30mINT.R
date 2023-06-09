######################
## Random Forest script that includes:
## Extraction of covariates to points
## Confustion matrix creation
## Kappa calculation
## rank123
## Inference Mask


# Workspace setup
# Install packages if not already installed

required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall","fmsb","reshape", "ggplot2", "parallel", "mgcv","UBL","itertools","doParallel","sf","rgeos","spatialEco","dplyr") #"tidyverse",
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package
rasterOptions(maxmemory = 3e+10,  memfrac = 0.9) # tmpdir = "/home/tnaum/data/temp/",
# memory.limit(140000)

# modelfolder <- "/home/tnaum/data/ESGs_work/soilgmrph_mapping"
datafolder <- "/home/tnaum/OneDrive/USGS/CPOG/PSCSmap/Update_202002"
covfolder <- "/home/tnaum/OneDrive/USGS/BLM_projects/BLM_CO_ESGs/ESGs_UCRB/ESG_30mCovsInt"
# maskfolder <- "E:/Models_active_work/ESGs/strawman_map/maskrasts"

######### Grid Prep #################
## Make list of grids
setwd(covfolder)
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid <- raster(cov.grids[1])
## Or Make a stack of grids to extract all at once (for smaller datasets)
#cov.stack <- stack()
cov.proj <- projection(projgrid)
  
######## Load points ##############
## NASIS with climate breaks
nasis <- readRDS("/home/tnaum/OneDrive/USGS/BLM_projects/BLM_CO_ESGs/ESGs_UCRB/pt_data/NASIS_pts_Jon_20190930/Jons_NASIS_All_Pts.rds")
nasis@data$LocID <- paste(nasis@coords[,1],nasis@coords[,2],sep="_")
nasis@data$project <- "nasis"
nasis <- spTransform(nasis, CRS(cov.proj)) # project to match rasters
nasis@data$compname <- toupper(nasis@data$compname)
## Attache particle size class and group to pscsmodorg
fgdb <- "/media/tnaum/ped/GIS_Archive/gSSURGO18/gSSURGO_CONUS.gdb/gSSURGO_CONUS.gdb"
comp.df <- sf::st_read(dsn = fgdb, layer = "component")
comp.df[] <- lapply(comp.df, function(x) if (is.factor(x)) as.character(x) else {x})
comp.dfc <- comp.df[,c("cokey","taxpartsize","compname","taxsuborder","taxorder","taxsubgrp")]
comp.dfc$compname <- toupper(comp.dfc$compname)
comp.dfc$taxsuborder <- toupper(comp.dfc$taxsuborder)
comp.dfc$taxpartsize <- toupper(comp.dfc$taxpartsize)
comp.dfc$taxorder <- toupper(comp.dfc$taxorder)
comp.dfc$taxsubgrp <- toupper(comp.dfc$taxsubgrp)
nasis@data <- cbind(comp.dfc[match(nasis@data$compname, comp.dfc$compname),], nasis@data) #attaches to pedtab
#Now clean up PSCS: include Psamments and Rock outcrops
nasis@data$taxpartsize = ifelse(grepl("ROCK OUTCROP", nasis@data$compname), "OUTCROPS", nasis@data$taxpartsize)
nasis@data$taxpartsize = ifelse(nasis@data$taxsuborder == "PSAMMENTS", "SANDY", nasis@data$taxpartsize)
nasis@data$taxpartsize = ifelse(nasis@data$taxpartsize == "PSAMMENTS" & grepl("LITHIC", nasis@data$taxsubgrp), "LITHIC SANDY", nasis@data$taxpartsize)
nasis@data$taxpartsize = ifelse(nasis@data$taxorder == "HISTOSOLS", "ORGANIC", nasis@data$taxpartsize)
nasis@data$taxpartsize = ifelse((nasis@data$taxpartsize != "NA") & (nasis@data$taxpartsize != "LOAMY") & (nasis@data$taxpartsize != "CLAYEY") & grepl("LITHIC", nasis@data$taxsubgrp), paste("LITHIC", nasis@data$taxpartsize, sep = " "), nasis@data$taxpartsize)
nasis@data$pscsmodorg <- nasis@data$taxpartsize 
nasis <- subset(nasis, nasis$pscsmodorg != "NA" & nasis$pscsmodorg != "NOT USED" & nasis$pscsmodorg != "LITHIC NOT USED")
nasis@data <- nasis@data[,c("LocID","pscsmodorg")]
## Bring in old point and combine by overlaying the two datasets
setwd("/media/tnaum/ped/SensitiveSoils/Data/models/NASIS_based/pscsmodorg")
nasis16 <- readOGR(".","nasis_pscsmodorg_UpCO")
nasis16@data[] <- lapply(nasis16@data, function(x) if (is.factor(x)) as.character(x) else {x})
nasis16@data$LocID <- paste(nasis16@coords[,1],nasis16@coords[,2],sep="_")
nasis16 <- spTransform(nasis16, CRS(cov.proj))
nasis16@data <- nasis16@data[,c("LocID","pscsmodorg")]
nasis16@data$pscsmodorg <- ifelse(nasis16@data$pscsmodorg == "ROCK OUTCROP","OUTCROPS",nasis16@data$pscsmodorg)
## create 5 meter buffer(polygons) around original NASIS points using gBuffer from "rgeos" package (requires SpatialPointsDataFrame)
NASIS_buffer <- gBuffer(nasis16, byid = FALSE, width = 5.0)
Jons_buffer <- gBuffer(nasis, byid = FALSE, width = 5.0)
## select points from Jon's points that fall within NASIS buffer polygons
Jonspts <- point.in.poly(nasis, NASIS_buffer, duplicate = TRUE)               #Returns a SpatialPointsDataFrame with a column "id.poly" where "1" indicates a shared point with original NASIS surveys
Jonspts <- as.data.frame(Jonspts) 
JonsPoints_shared <- subset(Jonspts, Jonspts$poly.ids == "1") 
#Convert SpatialPointsDataFrame to standard dataframe
## select points from NASIS points that fall within Jon's Points buffer polygons
NASISpts <- point.in.poly(nasis16, Jons_buffer, duplicate = TRUE)              #Returns a SpatialPointsDataFrame with a column "id.poly" where "1" indicates a shared point with Jons points locations
NASISpts <- as.data.frame(NASISpts)                                                #Convert SpatialPointsDataFrame to standard dataframe
NASISPoints_S <- subset(NASISpts, NASISpts$poly.ids == "1")                        #Subset of NASIS points that also co-occur with Jons data
NASISPoints_NS <- anti_join(NASISpts, NASISPoints_S)                               #Creates data frame with NASIS points not shared between original NASIS points and Jons points
All_Pts <- rbind(Jonspts, NASISPoints_NS)
## Create spatial pts dataframe
coordinates(All_Pts) <- c("coords.x1","coords.x2") # X and Y
crs(All_Pts) <- cov.proj
## Bring in USFS data
setwd("/home/tnaum/OneDrive/USGS/CPOG/PSCSmap/Update_202002/USFS_Colby")
usfs <- readOGR(".","SAMPLE_PT_ESRI_SHP")
usfs@data[] <- lapply(usfs@data, function(x) if (is.factor(x)) as.character(x) else {x})
usfs@data$LocID <- paste(usfs@coords[,1],usfs@coords[,2],sep="_")
usfs <- spTransform(usfs, CRS(cov.proj))
usfstax <- read.delim("SITE_CLASS_SOIL_ttab.txt",stringsAsFactors = F)
usfstaxc <- usfstax[,c("SITE_ID","SOIL_PSC","SOIL_NAME","SOIL_SUBORDER","SOIL_ORDER","SOIL_SUBGROUP")]
usfs@data <- cbind(usfstaxc[match(usfs@data$SITE_ID, usfstaxc$SITE_ID),], usfs@data) #attaches to pedtab
usfs@data$SOIL_NAME <- toupper(usfs@data$SOIL_NAME)
usfs@data$SOIL_PSC <- toupper(usfs@data$SOIL_PSC)
usfs@data$SOIL_SUBORDER <- toupper(usfs@data$SOIL_SUBORDER)
usfs@data$SOIL_ORDER <- toupper(usfs@data$SOIL_ORDER)
usfs@data$SOIL_SUBGROUP <- toupper(usfs@data$SOIL_SUBGROUP)
#Now clean up USFS PSCS: include Psamments and Rock outcrops
usfs@data$SOIL_PSC = ifelse(grepl("OUTCROP", usfs@data$SOIL_NAME), "OUTCROPS", usfs@data$SOIL_PSC)
usfs@data$SOIL_PSC = ifelse(usfs@data$SOIL_SUBORDER == "PSAMMENTS", "SANDY", usfs@data$SOIL_PSC)
usfs@data$SOIL_PSC = ifelse(usfs@data$SOIL_PSC == "PSAMMENTS" & grepl("LITHIC", usfs@data$SOIL_SUBGROUP), "LITHIC SANDY", usfs@data$SOIL_PSC)
usfs@data$SOIL_PSC = ifelse(usfs@data$SOIL_ORDER == "HISTOSOLS", "ORGANIC", usfs@data$SOIL_PSC)
usfs@data$SOIL_PSC = ifelse((usfs@data$SOIL_PSC != "NA") & (usfs@data$SOIL_PSC != "LOAMY") & (usfs@data$SOIL_PSC != "CLAYEY") & grepl("LITHIC", usfs@data$SOIL_SUBGROUP), paste("LITHIC", usfs@data$SOIL_PSC, sep = " "), usfs@data$SOIL_PSC)
usfs@data$pscsmodorg <- usfs@data$SOIL_PSC 
usfs@data <- usfs@data[,c("LocID","pscsmodorg")]
usfs <- subset(usfs, usfs$pscsmodorg != "NA" & usfs$pscsmodorg != "NOT USED" & usfs$pscsmodorg != "LITHIC NOT USED"|usfs$pscsmodorg != "")
usfs <- subset(usfs, usfs$pscsmodorg != "")
summary.factor(usfs@data$pscsmodorg)
usfs@data$project <- "usfs"
## Riparian points digitized by Sam Burch, USGS
setwd("/home/tnaum/OneDrive/USGS/BLM_projects/BLM_CO_ESGs/ESGs_UCRB/pt_data/Riparian_Sample_Points")
rip_pts34A <- readOGR(".", "Riparian_Points_MLRA34A_points")
rip_pts34A@data[] <- lapply(rip_pts34A@data, function(x) if (is.factor(x)) as.character(x) else {x})
rip_pts34B <- readOGR(".", "Riparian_Points_MLRA34B_points")
rip_pts34B@data[] <- lapply(rip_pts34B@data, function(x) if (is.factor(x)) as.character(x) else {x})
rip_pts35 <- readOGR(".", "Riparian_Points_MLRA35_points")
rip_pts35@data[] <- lapply(rip_pts35@data, function(x) if (is.factor(x)) as.character(x) else {x})
rip_pts36 <- readOGR(".", "Riparian_Points_MLRA36_points")
rip_pts36@data[] <- lapply(rip_pts36@data, function(x) if (is.factor(x)) as.character(x) else {x})
rip_pts <- do.call("rbind", list(rip_pts34A,rip_pts34B,rip_pts35,rip_pts36)) # combine all riparian points
rip_pts@data$pscsmodorg <- "Riparian"
rip_pts@data$LocID <- paste(rip_pts@data$X, rip_pts@data$Y, sep = "")
rip_pts@data$project <- "riparian"
rip_pts <- spTransform(rip_pts, CRS(cov.proj))
## Outcrop points from Colby Brungard
outcr_pts <- read.delim("/home/tnaum/OneDrive/USGS/BLM_projects/BLM_CO_ESGs/ESGs_UCRB/pt_data/Colby_depth/regMatrix_7.1.19.csv", stringsAsFactors = F, sep=",")
outcr_pts <- subset(outcr_pts, outcr_pts$DepthClass == 'BR')
outcr_pts <- subset(outcr_pts, outcr_pts$source == 'CaseBasedReasoning')
# outcr_pts <- subset(outcr_pts, outcr_pts$ELEVm_m < 2600)
outcr_pts$LocID <-  paste(outcr_pts$X, outcr_pts$Y, sep = "")
outcr_pts$coords.x1 <- outcr_pts$X
outcr_pts$coords.x2 <- outcr_pts$Y
coordinates(outcr_pts) <- ~ coords.x1 + coords.x2
outcr.proj <- CRS("+init=epsg:5070") # From Colby Brungard
projection(outcr_pts) <- outcr.proj
outcr_pts <- spTransform(outcr_pts, CRS(cov.proj))
outcr_pts@data$project <- "depth"
outcr_pts@data$pscsmodorg <- "OUTCROPS"
## Now merge different point
All_Pts$project <- "nasis"
nasism <- All_Pts[,c("project","LocID","pscsmodorg")]
projection(nasism) <- cov.proj
rip_ptsm <- rip_pts[,c("project","LocID","pscsmodorg")]
projection(rip_ptsm) <- cov.proj
outcr_ptsm <- outcr_pts[,c("project","LocID","pscsmodorg")]
shp.pts <- do.call("rbind", list(nasism,rip_ptsm,outcr_ptsm,usfs))

## Extract real elevation points to constrain psuedopoints
# ELEVm <- raster("E:/Models_active_work/UCRB_Covariates/ELEVm.tif")
# shp.pts <- extract(ELEVm, shp.pts, df=TRUE, sp=TRUE)
# shp.pts.pts <- subset(shp.pts, shp.pts$project != "depth")
# maxelev <- max(shp.pts.pts@data$ELEVm, na.rm=T)
# shp.pts <- subset(shp.pts, shp.pts$project != "depth" | (shp.pts$project =="depth" & shp.pts$ELEVm < 2350))

##Clip points outside of study area
setwd("/home/tnaum/OneDrive/USGS/BLM_projects/Utah_BLM_Salinity/Huc6_boundary")
polybound <- readOGR(".", "CO_River_watershed_Meade_alb")
polybound <- spTransform(polybound, cov.proj)
shp.pts <- shp.pts[polybound,]
shp.pts <- subset(shp.pts, shp.pts$pscsmodorg != "")

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
setwd(covfolder)
cpus <- detectCores(logical=TRUE)-2
rasterOptions(maxmemory = 4e+09)
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts", "cov.grids")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(cov.grids, function(i){try( raster::extract(raster(i), shp.pts) )}) 
snowfall::sfStop()
#detach(package:snowfall, unload=TRUE)
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) <- tools::file_path_sans_ext(basename(cov.grids))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts$DID <- seq.int(nrow(shp.pts))
pts <- merge(as.data.frame(shp.pts),ov.lst, by="DID")
pts[] <- lapply(pts, function(x) if (is.factor(x)) as.character(x) else {x})


## Save points and lookup table
setwd(datafolder)
saveRDS(pts, "UCRB_NASIS_mPSC_30mINT_ptsext.rds")
write.table(pts, "UCRB_NASIS_mPSC_30mINT_ptsext.txt", sep = "\\t", row.names = FALSE)
pts <- readRDS("UCRB_NASIS_mPSC_30mINT_ptsext.rds")
detach(package:dplyr) ## dplyr interferes with combine

## Standardize Point data classes
classname <- "mPSC"
pts$Class <- pts$pscsmodorg ## UPDATE EVERY TIME1
Class <- "pscsmodorg" ## Dependent variable

## Prep for Random Forest
formulaStringRF <- as.formula(paste('Class ~', paste(gsub(".tif","", cov.grids), collapse="+")))# put in dep variable name
ptsc <- subset(pts, pts$Class != "NA")
ptsc <- subset(ptsc, ptsc$Class != "")
ptsc <- na.omit(ptsc)# Remove any record with NA's (in any column - be careful)

## Examine class balaance
classnumb <- summary(as.factor(as.character(ptsc$Class)), maxsum=200)
classnumb
## Merge Gypsum and Ash classes
ptsc$Class <- ifelse(ptsc$Class=="COARSE-GYPSEOUS"|ptsc$Class=="FINE-GYPSEOUS","GYPSUM",ptsc$Class)
ptsc$Class <- ifelse(ptsc$Class=="ASHY"|ptsc$Class=="ASHY OVER LOAMY"|ptsc$Class=="ASHY OVER LOAMY-SKELETAL"|ptsc$Class=="ASHY OVER SANDY OR SANDY-SKELETAL","ASHY",ptsc$Class)
ptsc$Class <- ifelse(ptsc$Class=="ASHY-SKELETAL"|ptsc$Class=="ASHY-SKELETAL OVER FRAGMENTAL OR CINDERY"|ptsc$Class=="CINDERY","ASHY-SKELETAL",ptsc$Class)
classnumb <- summary(as.factor(as.character(ptsc$Class)), maxsum=200)
classnumb
## Now limit to more common classes
classlessnum <- subset(classnumb, classnumb < 23)
classtodrop <- names(classlessnum)
ptsc <- ptsc[ ! ptsc$Class %in% classtodrop, ]
classnumb <- summary(as.factor(as.character(ptsc$Class)), maxsum=200)
classnumb
ptscc <- subset(ptsc, ptsc$Class != "LITHIC ")
## Pull out duplicate locations
ptscc <- subset(ptscc, !duplicated(ptscc[c("LocID")])) #removes duplicates
## Clean up class names for use in lists
ptscc$Class <- as.character(ptscc$Class)
ptscc$Class <- gsub(" ","", ptscc$Class)
ptscc$Class <- gsub("&","", ptscc$Class)
ptscc$Class <- gsub("-","", ptscc$Class)
ptscc$Class <- as.factor(ptscc$Class)
## Split into train/test
nfolds <- 5
ptscc$sets <- sample.int(nfolds,size =length(ptscc[,1]),replace=T)
pts_rf <- subset(ptscc, ptscc$sets == 1 | ptscc$sets == 2 | ptscc$sets==3 | ptscc$sets == 4)
pts_test <- subset(ptscc, ptscc$sets == 5)
pts_test$Class <- as.character(pts_test$Class)
summary.factor(pts_test$Class)
## Pull out just regression matrix for SMOTE
covnames <- gsub(".tif","",cov.grids)
regmxnames <- c("Class",covnames)
pts_rgmtx <- pts_rf[,regmxnames]
pts_rgmtx$Class <- as.character(pts_rgmtx$Class)
## Smote percentage calculations
classnumb.rf <- summary(as.factor(as.character(pts_rgmtx$Class)), maxsum=200)
classnumb.rf
## Set up class weights for RF or synthetic oversampling
classwts <- max(classnumb.rf)/classnumb.rf # to get a fully balance oversampling
# classwts <- (1 - (classnumb.rf/max(classnumb.rf))+1)^1.5 # conservative oversampling
# classwts <- 1-(classnumb.rf^(1.2))/nrow(pts_rgmtx) ## For using classwts in RF
classwts
#classwts <- sqrt(classwts)
classwts <- classwts^(1/2)
classwts
classwtslst <- lapply(split(classwts, names(classwts)), unname)
## Adjust weights based on expert reasoning
classwtslst$FINELOAMY <- 0.65
classwtslst$OUTCROPS <- 1.75
classwtslst$SANDY <-0.75
classwtslst$FINE <- 1.1
classwtslst$LOAMYSKELETAL <- 0.55
classwtslst$LOAMY <- 1.2
classwtslst$GYPSUM <- 24
classwtslst$LITHICLOAMYSKELETAL <- 1.2
classwtslst$COARSELOAMY <- 0.95
classwtslst$FINESILTY <- 1.4
## Now create new balanced training set
pts_rgmtx$Class <- as.factor(pts_rgmtx$Class)
balpts <- SmoteClassif(formulaStringRF, pts_rgmtx, C.perc = classwtslst)
#balpts <- SmoteClassif(formulaStringRF, pts_rgmtx, C.perc = "balance")
summary(balpts$Class,maxsum=200)

############### Build Random Forest
cl <- makeCluster(getOption("cl.cores", cpus))
registerDoParallel(cl)
soiclass <- foreach(ntree=rep(floor(200/cpus), cpus), .combine=combine, .multicombine=TRUE, .verbose = T,
                    .packages='randomForest') %dopar% {
                      rf <- randomForest(formulaStringRF, balpts,ntree=ntree, importance=TRUE, proximity=FALSE, keep.forest=TRUE, nodesize=5) #,mtry=15,classwt=classwts)
                    } ## Does not return OOB statistics
stopCluster(cl)
class(soiclass) = "randomForest"
#soiclass = randomForest(formulaStringRF, data = ptsc, importance=TRUE, proximity=FALSE, ntree=100, keep.forest=TRUE, nodesize=7,mtry=15,classwt=classwts)
soiclass
### Calculate OOB error rate
balpts$predOOB = predict(soiclass)
balpts$oobmatch = ifelse(balpts$Class == balpts$predOOB, 1, 0)
rf.err.rate = 1 - sum(balpts$oobmatch)/(length(balpts[,1])) # OOB Error rate
varImpPlot(soiclass)
setwd(datafolder)
saveRDS(balpts, "pts_train_rf_soilgmrph_smote_30mINT.rds")
saveRDS(pts_test, "pts_test_rf_soilgmrph_smote_30mINT.rds")
saveRDS(soiclass,"rf_model_soilgmrph_smote_30mINT.rds")
saveRDS(pts_rf, "pts_rf_soilgmrph_smote_30mINT.rds")
## For re-runs
soiclass <- readRDS("rf_model_soilgmrph_smote_30mINT.rds")
# pts_rf <- readRDS("pts_train_rf_strawman2_v2_30mINT.rds")
############### Create Confusion matrix (for categorical models)
## Need to strip last column e.g. confusion[1:9,1:10] in rf object would be confusion[1:9,1:9]
obsclasses <- rownames(as.data.frame(summary(as.factor(as.character(balpts$Class)))))
oobclasses <- rownames(as.data.frame(summary(as.factor(as.character(balpts$predOOB)))))
comboclasses <- unique(append(obsclasses,oobclasses))
confmatrix <- as.data.frame(balpts[,c("Class","predOOB")])
confmatrix <- confmatrix[,c("Class","predOOB")]
oob_confmatx <- table(lapply(confmatrix, factor, levels = as.factor(comboclasses)))#levels need to come from field with most classes
OOBkappa <- Kappa.test(oob_confmatx, conf.level = 0.95)
write.table(oob_confmatx, file = "UCRB_NASIS_BLM_ESG_confmatrix_rf_soilgmrph_smote_30mINT.txt", sep = "\\t", row.names = TRUE) ## needs work to save rigtht
## Create lookup table (for categorical predictions)
lookup_tab <- as.data.frame(cbind(seq(1:length(soiclass$classes)),soiclass$classes))
colnames(lookup_tab) <- c("value","ESG")
lookup_tab$value <- as.character(lookup_tab$value)
lookup_tab$ESG <- as.character(lookup_tab$ESG)
write.table(lookup_tab, file <- "UCRB_BLM_ESGs_lookup_rf_soilgmrph_smote_30mINT.txt", sep = "\\t", row.names = FALSE)

## Look at test set statistics
pts_test$newpred <- as.character(predict(soiclass, newdata=pts_test))
pts_test$spredmatch <- ifelse(pts_test$newpred == pts_test$Class,1,0)
testsetaccur <- sum(pts_test$spredmatch, na.rm = T)/length(na.omit(pts_test$newpred))
## Confusion matrix and Kappa of test set validation
refclasses <- rownames(as.data.frame(summary(as.factor(as.character(pts_test$Class)))))
valclasses <- rownames(as.data.frame(summary(as.factor(as.character(pts_test$newpred)))))
allclasses <- unique(append(refclasses,valclasses))
conf_mat <- as.data.frame(pts_test[,c("Class","newpred")])
conf_mat <- conf_mat[,c("Class","newpred")]
val_confmatx <- table(lapply(conf_mat, factor, levels = as.factor(allclasses)))#levels need to come from field with most classes
valKappa <- Kappa.test(val_confmatx, conf.level = 0.95)
setwd(datafolder)
write.table(val_confmatx, file = "UCRB_ESG_val_confmatrix_rf_soilgmrph_smote_30m.txt", sep = "\\t", row.names = TRUE) ## needs work to save rigtht
## Now calculate accuracy if considering top 2 classes
pts_test_probs <- predict(soiclass,newdata=pts_test,type="prob")
pts_test_probs <- as.data.frame(pts_test_probs)
pts_test$predclass1 <- names(pts_test_probs)[apply(pts_test_probs,1,which.max)]
pts_test$predclass2 <- colnames(pts_test_probs)[apply(pts_test_probs[,1:15], 1, function(x)
  which(x == sort(x, decreasing = TRUE)[2])[1])]
pts_test$predclass3 <- colnames(pts_test_probs)[apply(pts_test_probs[,1:15], 1, function(x)
  which(x == sort(x, decreasing = TRUE)[3])[1])]
pts_test$mtchtest <- ifelse(pts_test$newpred == pts_test$predclass1,1,0)
sum(pts_test$mtchtest) # should equal the length of the data frame
pts_test$spred2match <- ifelse(pts_test$predclass2 == pts_test$Class,1,0)
pts_test$spred3match <- ifelse(pts_test$predclass3 == pts_test$Class,1,0)
pts_test$top2match <- pts_test$spred2match + pts_test$spredmatch
pts_test$top3match <- pts_test$spred2match + pts_test$spredmatch + pts_test$spred3match
Top2testsetaccur <- sum(pts_test$top2match, na.rm = T)/length(na.omit(pts_test$newpred))
Top3testsetaccur <- sum(pts_test$top3match, na.rm = T)/length(na.omit(pts_test$newpred))
# predclass <- names(pts_test_probs)[apply(pts_test_probs,1,which.max)]
# predclass2 <- colnames(pts_test_probs)[apply(pts_test_probs[,1:15], 1, function(x)
#   which(x == sort(x, decreasing = TRUE)[2])[1])]
# maxprob <- apply(pts_test_probs,1,max)
# testprepprob <- data.frame(predclass,maxprob,Class=pts_test$Class,stringsAsFactors = F)


## Reference covar rasters to use in prediction
setwd(covfolder)
rasters <- stack(cov.grids)
#rasters = setMinMax(brick(rasters))
names(rasters)

## Predict onto covariate grid
setwd(datafolder)
cpus <- detectCores(logical=TRUE)-2
#memfree <- as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo",intern=TRUE))*0.9 # available RAM in kb
memfree <- 235000000 #Avail RAM in kb
maxmem <- ((memfree/(cpus))*0.95)*1000 # max mem per core in bytes: could be little higher
chunksz <- maxmem*0.35
#rasterOptions(maxmemory = maxmem, chunksize = chunksz)
rasterOptions(maxmemory = 7e+09, chunksize = 1e+09)# 5e+09, chunksize = 2e+08)
## Parallelized predict
beginCluster(cpus,type='SOCK')
Sys.time()
pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
Sys.time()
# predprob makes a huge temp file on c drive (137+GB for COP), make sure there's room...
# predprob <- clusterR(rasters, predict, args=list(model=soiclass, type="prob", index = 1:length(soiclass$classes)),progress="text")
# writeRaster(predprob, overwrite=TRUE,filename="UCRB_BLM_ESGs_rf_initial_strawman_probmatrix_flt_30m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
endCluster()
gc()
# beginCluster(cpus,type='SOCK')
# predprobstk <- stack("UCRB_BLM_ESGs_rf_initial_strawman_probmatrix_flt_30m.tif")
# ## Rank function to get 2nd and 3rd most probable classes
# # rank function creates a rasterstack with each layer representing one class and each pixel representing the prediction rank of that class
# # The ranks go from highest number (most likely class) to lowest number (least likely class)
# probrank <- clusterR(predprobstk, overlay, args=list(fun=function(x) (rank(x,ties.method= "random"))),progress = "text")
# writeRaster(probrank, overwrite=TRUE,filename="UCRB_BLM_ESGs_rf_rank_30m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# pred2 <- pred
# pred3 <- pred
# pred2[] <- 0
# pred3[] <- 0
# for(no in lookup_tab$value){
#   num <- as.numeric(no)
#   rank2num <- as.numeric(length(lookup_tab$value))-1
#   rank3num <- as.numeric(length(lookup_tab$value))-2
#   rast <- probrank[[num]]
#   rast2_stk <- stack(rast,pred2)
#   pred2_fn <- function(rast,pred2) {
#     ind <- ifelse((rast==rank2num),num,pred2)
#     return(ind)
#   }
#   rast3_stk <- stack(rast,pred3)
#   pred3_fn <- function(rast,pred3) {
#     ind <- ifelse((rast==rank3num),num,pred3)
#     return(ind)
#   }
#   pred2 <- clusterR(rast2_stk, overlay, args=list(fun=pred2_fn), progress='text',export=c('rank2num','num'))
#   pred3 <- clusterR(rast3_stk, overlay, args=list(fun=pred3_fn), progress='text',export=c('rank3num','num'))
#   print(paste(no, "is done", sep=" "))
# }
# ## Now mask to original extent
# mask_fn <-function(pred,predn) {
#   ind <- ifelse(pred > 0,predn,NA)
#   return(ind)
# }
# pred2m_stk <- stack(pred,pred2)
# pred2 <- clusterR(pred2m_stk, overlay, args=list(fun=mask_fn), progress='text')
# pred3m_stk <- stack(pred,pred3)
# pred3 <- clusterR(pred3m_stk, overlay, args=list(fun=mask_fn), progress='text')
# writeRaster(pred2, overwrite=TRUE, filename = paste(classname, "2ndclass.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# writeRaster(pred3, overwrite=TRUE, filename = paste(classname, "3rdclass.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# 
# ## Now produce probability surfaces
# probmax <- clusterR(predprobstk, overlay, args=list(fun=function(x) (max(x))),progress = "text")
# names(probmax) <- "probmax"
# names(probmax)
# probmaxint <- clusterR(probmax, calc, args=list(fun=function(x) (x*100)),progress = "text")
# predprobint <- clusterR(predprobstk, calc, args=list(fun=function(x) (x*100)),progress = "text")
# endCluster()
# gc()
setwd(datafolder)
writeRaster(pred, overwrite=TRUE,filename="UCRB_mPSC_smote_30mINT.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# writeRaster(predprobint, overwrite=TRUE,filename="UCRB_BLM_ESGs_rf_initial_strawman_probmatrix_int_30m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# writeRaster(probmaxint, overwrite=TRUE,filename="UCRB_BLM_ESGs_rf_initial_strawman_probmax_int_30m.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U', progress="text")
# pred <- raster("UCRB_BLM_ESGs_rf_initial_strawman_30m.tif")





################# Inference space mask
# vars <- names(sort(varImpPlot(soiclass)[,1],decreasing = T)[0:5]) # pick top n variables in model to mask by
# vars <- c("ELEVm", "ppt_ann","temp_ann") # Or pick yourself
# samp <- ptsc
# ## set how far outside of min and max you would like to let in mask
# min_margin <- 0.02
# max_margin <- 0.02
# ## list apply to compute dummy mask rasters for each variable to then combine in a min function
# mask_fn <- function(rast){
#   setwd(covfolder)
#   rastfile <- paste(rast,'.tif', sep="")
#   rastobj <- raster(rastfile)
#   lowbound <- min(samp[c(rast)]) - abs((min(samp[c(rast)]))*min_margin)
#   highbound <- max(samp[c(rast)]) + abs((max(samp[c(rast)]))*max_margin)
#   newmask <-  calc(rastobj, fun=function(x){ifelse(x>lowbound & x<highbound,1,0)})
#   setwd(maskfolder)
#   writeRaster(newmask, overwrite=TRUE,filename=paste(rast,"mask.tif",sep=""), options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U')
#   gc()
# }
# ## Setup up parallel list apply
# rasterOptions(maxmemory = 1e+07, chunksize = 5e+06) # Adjust to keep RAM from blowing up
# snowfall::sfInit(parallel=TRUE, cpus=16) ## Choose number of cpus available
# snowfall::sfExport("vars","samp", "covfolder","mask_fn","maskfolder","min_margin","max_margin")
# snowfall::sfLibrary(rgdal)
# snowfall::sfLibrary(raster)
# Sys.time()
# snowfall::sfLapply(vars, function(rast){mask_fn(rast)})
# Sys.time()
# snowfall::sfStop()
# ## Now combine all variable raster masks to get an overall mask
# rasterOptions(maxmemory = 7e+10, chunksize = 7e+09)
# setwd(maskfolder)
# # mask.grids <- list.files(pattern=".tif$")
# mask.grids <- paste(vars, "mask", ".tif", sep="")
# mask.rasts <- stack(mask.grids)
# maskoverlay.fn <- function(mask.rasts) {
#   ind <- min(mask.rasts)
#   return(ind)
# }
# beginCluster(30,type='SOCK')
# mask.overlay <- clusterR(mask.rasts, overlay, args=list(fun=maskoverlay.fn),progress = "text")
# setwd(modelfolder)
# writeRaster(mask.overlay, overwrite=TRUE,filename="strawman_mask_climate_overlay.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U',progress="text")
# ## Now combine a water mask and the inference mask to create final layers
# f_mask <- function(a,b,c) {
#   ind <- a*b*c
#   return(ind)
# }
# h2omask <- raster("E:/Models_active_work/UpCo/nlcd_watermask.tif")
# msk_stk <- stack(pred,h2omask,mask.overlay)
# pred_msk <- clusterR(msk_stk, overlay, args=list(fun=f_mask),progress = "text")
# writeRaster(pred_msk, overwrite=TRUE,filename="ESGs_strawman_masked_climate.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"),datatype='INT1U',progress="text")
# endCluster()

