#Load required packages
library(raster)
library(dplyr)
library(dismo)
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(wallace)
library(ggplot2)
library(remotes)
library(wallace)
library(rgdal)
library(blockCV)
library(sf)
library(dplyr)
library(crs)
library(ggplot2)
library(usdm)
library(reshape2)

#set seed 
set.seed(48)

#Upload occurences, make sure columns are organized by 'OccID' (sequential numbers given to each occurence), 'scientific_name', 'longitude', 'latitude', and 'year'
occs <- read.csv("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Synthemis_Full_Data.csv")

#Create bounding box of species range, bgExt package is from Wallace beta version.
bgExt <- penvs_drawBgExtent(
  polyExtXY = matrix(c(137.535733, 140.25931, 143.773602, 145.794319, 148.605753, 149.66004, 151.461115, 152.295759, 153.56969, 154.096834, 154.272549, 154.272549, 153.393976, 151.900402, 146.541106, 143.334315, 139.951809, 137.887163, 135.998231, 136.393589, 137.535733, -37.428661, -38.879942, -39.084904, -39.357263, -39.01665, -38.192458, -36.480518, -35.377434, -33.382863, -31.378994, -29.024503, -27.98182, -26.732886, -26.497154, -29.024503, -30.322656, -31.566397, -33.309442, -35.377434, -37.428661, -37.428661),ncol=2,byrow=FALSE), 
  polyExtID = 2764, 
  drawBgBuf = 0, 
  occs = occs)

#Save the bounding box for future use for cropping enviromental variables
writeOGR(bgExt, dsn = "/Users/aarongoodman/Desktop", layer ='bgExt', driver = 'ESRI Shapefile', overwrite = T)

#Spatially thin occurrences by 5km and filter out years with less than 5 occurences to avoid overinflation in modelbuilding
output <- spThin::thin(occs, 'latitude', 'longitude', 'scientific_name', thin.par = 5, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
occs <- occs[as.numeric(rownames(maxThin)),] 
occs <- occs[!duplicated(occs),]
occs<- occs %>% group_by(year) %>% filter(n() >= 5)
occs<-as.data.frame(occs)

#Make for loop which looks for yearly directories of environmental variables
#Stack the enviromental variables per year, and crop them by the bounding box
#Drop bioclim layers 08, 09, 18, and 19 as they cause interpolation problems in later analyses
years <- c(2001:2020)
list_dirs<- list.dirs("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Resampled", recursive = F) 
envs.list <- list()
for(i in list_dirs){
  setwd(i) 
  s <-list.files(path=setwd(i), pattern = ".tif", recursive=F)
  s <-stack(s)
  c <-crop(s, bgExt)
  c <- dropLayer(c, c("bio08","bio09","bio18","bio19","nontreecover"))
  envs.list[[i]]<-c
}

#Name the list elements with the year numbers
names(envs.list) <- years

#Make a new list to store the year data frames in
envs.dfs <- list()

#Loop over the number of years, get the year for this iteration, get data frame of enviromental variables with values and coordinates for each cell, #Get cell numbers associated with the xy values, 
#Make new field for year, put this dataframe in the list
for(i in years) {
  envs.i <- envs.list[[as.character(i)]]
  envs.i.df <- raster::as.data.frame(envs.i, xy = TRUE)
  envs.i.df$cell<- raster::cellFromXY(envs.i, envs.i.df[,1:2])
  envs.i.df$year <- i
  envs.dfs[[i]] <- envs.i.df
}

#extract enviromental variable cells, which be any of the rasters or the whole stack
occs$cell <- raster::extract(envs.i, occs[3:4], cellnumbers = TRUE)[,1]

#join enviormental variable values and bind grid cells 
envs.df.all <- bind_rows(envs.dfs)
occs.match <- dplyr::left_join(occs, envs.df.all, by = c("cell","year"))

#Generate 50,000 random background points 
bg <- dismo::randomPoints(envs.list$`2001`, n = 50000) %>% as.data.frame()

# Repeat the same loop to get the enviromental values again
for(i in years) {
  envs.i <- envs.list[[as.character(i)]]
  envs.i.df <- raster::as.data.frame(envs.i, xy = TRUE)
  envs.i.df$cell<- raster::cellFromXY(envs.i, envs.i.df[,1:2])
  envs.i.df$year <- i
  envs.dfs[[i]] <- envs.i.df
}

#Extract cell values from background points from all years
bg$cell <- raster::extract(envs.i, bg, cellnumbers = TRUE)[,1]

#Bind rows of environmental data values and merge background cell values by matching cell values 
envs.df.all <- bind_rows(envs.dfs)
bg.match <- dplyr::left_join(bg, envs.df.all, by = c("cell"))
ncol(bg.match)
#Create a loop through each year which finds the variables which are uncorrelated and make a table of variables which are common throughout each year 
corr.list<-list()

for (i in years){
#Calculate correlation among variables using VIFcor 
bg.test<-bg.match %>% filter(year == i)
ncol(bg.test)
r<-bg.test[,6:27]
# calculates vif for the variables in 
v1 <- vifcor(r, th=0.9) # identify collinear variables that should be excluded
v1
re1 <- exclude(r,v1) # exclude the collinear variables that were identified in
# the previous step
re1
v2 <- vifstep(r, th=10) # identify collinear variables that should be xcluded
v2
re2 <- exclude(r, v2) # exclude the collinear variables that were identified in
# the previous step
re2
re3 <- exclude(r) # first, vifstep is called
re3
corr.list[[i]]<-names(re3)
}

#convert list to contingency table to determine which variables are consistent throughout the datasets
freq<-as.data.frame(table(unlist(corr.list)))
freq_20<-filter(freq, Freq == 20)
freq_20

#subset the enviromental variables which are uncorrelated to generate occurence dataset for SWD maxent model 

#myvars<- c("longitude","latitude", "bio02", "bio03","bio14","bio15","biome_01"#,"evapotranspiration","fall","nontreecover","nonvegcover","spring","summer"#,"treecover","vegedata","winter")
#occs.z<-occs.match[myvars]

myvars<- c("longitude","latitude", "bio02", "bio03", "bio14", "bio15", "biome_01", "evapotranspiration", "fires", "treecover", "nonvegcover", "vegedata")
occs.z<-occs.match[myvars]

#subset the enviromental variables which are uncorrelated to generate   mean and mode values of background 

#myvars<- c("bio02", "bio03","bio14","bio15","biome_01","evapotranspiration","fall","nontreecover","nonvegcover","spring","summer","treecover","vegedata","winter", "cell")

myvars<- c("bio02", "bio03", "bio14", "bio15", "biome_01", "evapotranspiration", "fires", "nonvegcover", "treecover", "vegedata", "cell")

bg.match_sub <- bg.match[myvars]

#upload function for determining mode of categorical variables 
mode <- function(codes){
  which.max(tabulate(codes))
}

#summarize the dataset by determing the mean of continuous values and mode of categorical values of entire background for each year 
#bg.all.mean <- bg.match_sub %>% group_by(cell) %>% 
  #dplyr::summarise(
            #bio02 = mean(bio02),
            #bio03 = mean(bio03),
            #bio14 = mean(bio14), 
            #bio15 = mean(bio15), 
            #biome_01 = max(biome_01), 
            #evapotranspiration = mean(evapotranspiration),
            #fall = max(fall),
            #nontreecover = mean(nontreecover), 
            #nonvegcover = mean(nonvegcover),
            #spring = max(spring),
            #summer = max(summer),
            #treecover = mean(treecover),
            #vegedata = max(vegedata),
            #winter = max(winter))

bg.all.mean <- bg.match_sub %>% group_by(cell) %>% 
  dplyr::summarise(
    bio02 = mean(bio02),
    bio03 = mean(bio03), 
    bio14 = mean(bio14),
    bio15 = mean(bio15), 
    evapotranspiration = mean(evapotranspiration),
    nonvegcover = mean(nonvegcover), 
    treecover = mean(treecover),
    biome_01 = max(biome_01), 
    vegedata = max(vegedata),
    fires = max(fires))


#merge the background enviromental values with latitude, longitude and cell of original background dataset, 
bg.all.mean <- dplyr::left_join(bg, bg.all.mean, by = c("cell"))

#Rename x and y to longitude and latitude to keep occurence and background datasets consistent 
names(bg.all.mean)[names(bg.all.mean) == 'x'] <- 'longitude'
names(bg.all.mean)[names(bg.all.mean) == 'y'] <- 'latitude'

#Remove cell values 
drops <- c("cell")
bg.all.mean <-bg.all.mean[ , !(names(bg.all.mean) %in% drops)]

##Repeat for maximum values of each background
#bg.all.max <- bg.match_sub %>% group_by(cell) %>% dplyr::summarise_at(vars#(bio03:vegedata), max)
#bg.all.max <- dplyr::left_join(bg, bg.all.max, by = c("cell"))
#names(bg.all.max)[names(bg.all.max) == 'x'] <- 'longitude'
#names(bg.all.max)[names(bg.all.max) == 'y'] <- 'latitude'
#drops <- c("cell")
#bg.all.max <-bg.all.max[ , !(names(bg.all.max) %in% drops)]
#
##Repeat for minimum values of each background
#bg.all.min <- bg.match_sub %>% group_by(cell) %>% dplyr::summarise_at(vars#(bio03:vegedata), max)
#bg.all.min <- dplyr::left_join(bg, bg.all.min, by = c("cell"))
#names(bg.all.min)[names(bg.all.min) == 'x'] <- 'longitude'
#names(bg.all.min)[names(bg.all.min) == 'y'] <- 'latitude'
#drops <- c("cell")
#bg.all.min <-bg.all.min[ , !(names(bg.all.min) %in% drops)]
#
##Make loop of yearly directories with enviromental variables in them,
##stack the variables, and crop them by the bounding box, but now remove all #uncorrelated variables

envs.list_uncorr <-list()
years <- c(2001:2020)
list_dirs<- list.dirs("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Resampled", recursive = F) 
for(i in list_dirs){
  setwd(i) 
  s <-list.files(path=setwd(i), pattern = ".tif", recursive=F)
  s <-stack(s)
  c <-crop(s, bgExt)
  c <- subset(c, c("bio02", "bio03", "bio14", "bio15", "biome_01", "evapotranspiration", "fires", "treecover", "nonvegcover", "vegedata"))
  envs.list_uncorr[[i]]<-c
}

#Name the list elements with the year numbers
names(envs.list_uncorr) <- years

#Run the dynamic maxent models for simple and complex models
permute_list_dynamic<-list()
results_list_dynamic<-list()

list<-c("L", "LQH") 
for (i in list) {
e.swd <- ENMevaluate(occs.z, bg = bg.all.mean, algorithm = "maxent.jar", tune.args = list(fc = c(i), rm = 1), categoricals = c("biome_01", "vegedata", "fires"), partitions = "randomkfold")

#Generate list of your maxent results for both levels of complexity
res <- eval.results(e.swd) 
results_list_dynamic[[i]]<-res

#Generate 'cloglog' predictions of your simple and complex dynamic models and save them to a directory 
for (j in names(envs.list_uncorr)) {
  predictions <-predict(envs.list_uncorr[[j]], e.swd@models[[res$tune.args]], type = 'cloglog')
  writeRaster(predictions, file.path("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Dynamic",i, j), format = "GTiff", overwrite = T)

#Generate 10% omission predictions of your simple and complex dynamic models and save them to a directory
occPredVals_Se <- raster::extract(predictions, occs.z[,1:2])
thresProb_Se <- switch("p10", "mtp" = 0, "p10" = 0.1, "qtp" = 0)
thres_Se <- stats::quantile(occPredVals_Se, probs = thresProb_Se)
predSel_Se <- predictions > thres_Se
    writeRaster(predSel_Se, file.path("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results/Dynamic_binary/",i, j), format = "GTiff", overwrite = T)

#Use the @variable.importance function to find which variables posess the highest permutation importance within your models
var<-e.swd@varimp

#Create a list of your dynamic model complexity and fill it with your permutation reslts 
permute_list_dynamic[[i]]<-var

  }
}

#Save permutation results and maxent results as a csv
setwd("/Users/aarongoodman/Library/CloudStorage/OneDrive-AMNH/Synthemis\\ eustalacta/Synthemis_New/Results")

permute_list_dynamic_test<-bind_rows(permute_list_dynamic)
write.csv(permute_list_dynamic, "permute_list_dynamic.csv")

results_list_dynamic<-bind_rows(results_list_dynamic)
write.csv(results_list_dynamic, "results_list_dynamic.csv")

