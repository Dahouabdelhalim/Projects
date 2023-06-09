#Created 2/9/2021 by Sara Hamilton
#For Pycnopodia manuscript
#Creates a "Species With Data" format dataset using Pycnopodia presence data, sampling effort, and 
#environmental data. This SWD data is then used in the Maxent terminal program. Running final 
#models through the Maxent terminal program because I have not yet found a way to adjust Tau 
#when running maxent through R using DISMO or ENMEval
library(tidyverse)
library(raster)
library(MASS)
library(spatial)
WGS <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs "
Albers2 <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
###################################################################################
#1. Creating Occurence data from Master Data
###################################################################################
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Maxent_SWD")
data <- read.csv("OccurenceData/MasterPycno_ToShare.csv")
#creating a pre-outbreak (2009-2012) list of presences and keeping only the coordinate columns. No other info is need beyond location of presences
occ_pre <-  data %>%
  filter(pres_abs == 1, event_year >= 2009, event_year <= 2012) %>%
  dplyr::select(dec_long, dec_lat)
pycnopodia <- rep("Pycnopodia helianthoides",nrow(occ_pre))
occ_pre <- data.frame(pycnopodia, occ_pre)
colnames(occ_pre)<- c("Species","Longitude","Latitude")
#Coordinates are currently in WGS and need to be reprojected into Albers N. America Equal Area Conic (see Proj4 code above)
occ_pre <- SpatialPoints(occ_pre[,2:3],proj4string = crs(WGS))
occ_pre <- spTransform(occ_pre, CRSobj = crs(Albers2))
occ_pre <- occ_pre@coords
occ_pre <- data.frame(pycnopodia, occ_pre)
colnames(occ_pre)<- c("Species","Longitude","Latitude")
write.csv(occ_pre, "OccurenceData/Presence_20092012_Albers.csv")

#creating a post-oubreak (2017-2020), presences files in Albers, as done above for the pre-outbreak file
occ_post <-  data %>%
  filter(pres_abs == 1, event_year >= 2017) %>%
  dplyr::select(dec_long, dec_lat)
pycnopodia <- rep("Pycnopodia helianthoides",nrow(occ_post))
occ_post <- data.frame(pycnopodia, occ_post)
colnames(occ_post)<- c("Species","Longitude","Latitude")
occ_post <- SpatialPoints(occ_post[,2:3],proj4string = crs(WGS))
occ_post<- spTransform(occ_post, CRSobj = crs(Albers2))
occ_post<- occ_post@coords
occ_post<- data.frame(pycnopodia, occ_post)
colnames(occ_post)<- c("Species","Longitude","Latitude")
write.csv(occ_post, "Presence_20172020_Albers.csv")

#############################################################################################
#2. Create Bias Files
#############################################################################################
#Creat 'bias file' for pre-outbreak period. Bias file quantifies sampling effort across the study region based on all available
#surveys we drew from during that period (including those with Pycnopodia absences)
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Maxent_SWD")
pycnos <- read.csv("OccurenceData/MasterPycno_ToShare.csv", stringsAsFactors = FALSE)
pre <- pycnos %>%
  filter(event_year <= 2012 & event_year >= 2009) %>%
  dplyr::select(dec_long, dec_lat)
#transform coordinates into North America Albers Equal Area Conic
pre <- SpatialPoints(pre,proj4string = crs(WGS))
pre<- spTransform(pre, CRSobj = crs(Albers2))
pre <- pre@coords
#now open up our depth file which defines our study region (all areas in their known range with a 0 >= depth >= -455)
depth <-raster("Maxent_EnvironmentalFiles/Pre/Depth.asc")
#create a kernel density plot of sampling effort using the locations of each sample. Kernel Density plot should have the same resolution
#and extent as the original depth file
pre_kde <- kde2d(pre[,1], pre[,2], n = c(ncol(depth), nrow(depth)), 
                 lims = c(extent(depth)[1],extent(depth)[2], extent(depth)[3],extent(depth)[4]))
#put that kernel density plot into rater format
pre_kde <- raster(pre_kde, crs = Albers2)
#resample that raster so that it has the same extent and resolution as the depth graph. 
template <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], ymn = extent(depth)[3],
                   ymx = extent(depth)[4], crs = Albers2, res = c(res(depth)[2], res(depth)[2]))
pre_kde<- resample(pre_kde, template, method = 'bilinear')
pre_kde<- mask(pre_kde, depth)
writeRaster(pre_kde,"BiasFiles/KDE_PreOutbreak_Albersresampled",
            format = "ascii", overwrite = TRUE)
writeRaster(pre_kde,"BiasFiles/KDE_PreOutbreak_Albersresampled",
            format = "GTiff",overwrite = TRUE)


#Follow the same procedure as above to create a post-outbreak (2017-2020) bias file 
post <- pycnos %>%
  filter(event_year >= 2017) %>%
  dplyr::select(dec_long, dec_lat)
post <- SpatialPoints(post,proj4string = crs(WGS))
post <- spTransform(post, CRSobj = crs(Albers2))
post <- post@coords
post_kde <- kde2d(post[,1], post[,2], n = c(ncol(depth), nrow(depth)), 
                  lims = c(extent(depth)[1],extent(depth)[2], extent(depth)[3],extent(depth)[4]))
post_kde <- raster(post_kde, crs = Albers2)
post_kde<- resample(post_kde, template, method = 'bilinear')
post_kde<- mask(post_kde, depth)
writeRaster(post_kde,"BiasFiles/KDE_PostOutbreak_Albersresampled",
            format = "ascii", overwrite = TRUE)
writeRaster(post_kde,"KDE_PostOutbreak_Albersresampled",
            format = "GTiff", overwrite = TRUE)

##################################################################################
#3. Creating SWD data for Maxent analysis 
#####################################################################################
#Extracting environmental variable data for each pre-outbreak (2009-2012) Pycnopodia presence point to create the SWD format for Maxent
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Maxent_SWD")
occ <- read.csv("OccurenceData/Presence_20092012_Albers.csv")
occ <- occ[,-1]
#Extract in this order: 1) depth 2) 90th percentile of temps 3) mean Chl, 4) mean salinity and 5) substrate type 
var <- raster("Maxent_EnvironmentalFiles/Pre/Depth.asc")
crs(var) <- Albers2
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth")
#temp
var <- raster("Maxent_EnvironmentalFiles/Pre/90th_%_Temp.asc")
crs(var) <- Albers2
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp")
#chl
var <- raster("Maxent_EnvironmentalFiles/Pre/Mean_Chl.asc")  
crs(var) <- Albers2  
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl")
#sal
var <- raster("Maxent_EnvironmentalFiles/Pre/Salinity.asc")
crs(var)<- Albers2
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl", "Salinity")
#substrate type
var <- raster("Maxent_EnvironmentalFiles/Pre/Substrate.asc")
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl", "Salinity","Substrate")
#Maxent prefers -9999 values to NAs
occ_swd[is.na(occ_swd)] <- -9999
write.csv(occ_swd,"OccurenceData/Presence_20092012_Albers_SWD.csv")

#and Now doing the same thing for post-outbreak data
occ <- read.csv("OccurenceData/Presence_20172020_Albers.csv")
occ <- occ[,-1]
##Extract in this order: 1) depth 2) 90th percentile of temps 3) mean Chl, 4) mean salinity and 5) substrate type 
#depth
depth <- raster("Maxent_EnvironmentalFiles/Post/Depth.asc")
crs(depth) <- Albers2
occ_var <- raster::extract(depth, occ[,2:3])
occ_swd <- cbind(occ, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth")
#temp
var <- raster("Maxent_EnvironmentalFiles/Post/90th_%_Temp.asc")
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp")
#chl
var <- raster("Maxent_EnvironmentalFiles/Post/Mean_Chl.asc")  
crs(var) <- Albers2  
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp", "Mean_Chl")
#sal
var <- raster("Maxent_EnvironmentalFiles/Post/Salinity.asc")
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp", "Mean_Chl", "Salinity")
#Substrate type
var <- raster("Maxent_EnvironmentalFiles/Post/Substrate.asc" )
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp", "Mean_Chl", "Salinity","Substrate")
#Done!
occ_swd[is.na(occ_swd)] <- -9999
write.csv(occ_swd,"OccurenceData/Presence_20172020_Albers_SWD.csv")

#############################################################################################
#4. Creating background SWD data for Maxent analysis 
#############################################################################################
#Creating 5000 background points across the defined range for Maxent. THese background
#points are created using the sampling probabilities from the "Bias Files'(which are essentially a KDE of sampling effort) 
#in order to minimize the impacts of spatial sampling bias on the results. Once background points are created
#this code then extracts environmental variable for each background point to create a SWD dataset similar to the one described 
#above for Pycnopodia present points

#Pre_outbreak
#Loading bias file, which is a KDE of sampling effort from 2009-2012
bias <- raster("BiasFiles/KDE_PreOutbreak_Albersresampled.tif")
#transform from raster to vector
bias2 = bias[1:ncell(bias)]
bias2[is.na(bias2)] <- 0
cells = seq(1: length(bias))
#choosing 5000 points from all available cells usings the bias file sampling probability
backgrounds <- sample(cells,5000, replace = TRUE, prob = bias2)
#Convert cells into actual lat/lonsby finding centroid of cell and then add jitter to create random point in the cell
cell_width_height <- res(bias)/2 - 1
backgroundpoints = data.frame(Longitude = numeric(),
                              Latitude = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(backgrounds)){
  #Grab lat/lon of center of each cell selected for background points
  latlon <- xyFromCell(bias, backgrounds[i])
  #add jitter to the lon
  backgroundpoints[i,1] <- latlon[1] + runif(1,min = -cell_width_height, max = cell_width_height)
  #add jitter to the lat
  backgroundpoints[i,2] <- latlon[2] + runif(1,min = -cell_width_height, max = cell_width_height)
}
#Once I add a column indicating these are background points (for Maxent terminal) then we have our
#spatially biased background points ready to go
occ <- cbind(rep("background",5000),backgroundpoints)
#Now extracting environmental variable data for each background point as done in section 3
var <- raster("Maxent_EnvironmentalFiles/Pre/Depth.asc")
crs(var) <- Albers2
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth")
#temp
var <- raster("Maxent_EnvironmentalFiles/Pre/90th_%_Temp.asc")
crs(var) <- Albers2
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp")
#chl
var <- raster("Maxent_EnvironmentalFiles/Pre/Mean_Chl.asc")  
crs(var) <- Albers2  
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl")
#sal
var <- raster("Maxent_EnvironmentalFiles/Pre/Salinity.asc")
crs(var)<- Albers2
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl", "Salinity")
#substrate type
var <- raster("Maxent_EnvironmentalFiles/Pre/Substrate.asc")
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd)<- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl", "Salinity","Substrate")
#Maxent prefers -9999 values to NAs
occ_swd[is.na(occ_swd)] <- -9999
write.csv(occ_swd,"OccurenceData/Absence_20092012_Albers_SWD.csv")


#Now do the same thing for post-outbreak (2017-2020)
#KDE representing sampling effort from 2017-2020
bias <- raster("BiasFiles/KDE_PostOutbreak_Albersresampled.tif")
crs(bias) <- Albers2
#transforming to a vector instead of a raster
bias2 = bias[1:ncell(bias)]
bias2[is.na(bias2)] <- 0
#making a vector of the cell IDs
cells = seq(1: length(bias))
#sampling cell numbers using the probabilities contained in the bias file
backgrounds <- sample(cells,5000, replace = TRUE, prob = bias2)
#Convert cell number to lat long. Find centroid of each cell and then add jitter
cell_width_height <- res(bias)/2 - 1
backgroundpoints = data.frame(Longitude = numeric(),
                              Latitude = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(backgrounds)){
  latlon <- xyFromCell(bias, backgrounds[i])
  backgroundpoints[i,1] <- latlon[1] + runif(1,min = -cell_width_height, max = cell_width_height)
  backgroundpoints[i,2] <- latlon[2] + runif(1,min = -cell_width_height, max = cell_width_height)
}
#Add Maxent-required 'species' columna (where species here is just background point)
occ <- cbind(rep("background",5000),backgroundpoints)
#Now extracting environmental variable data for each background point. Extract in this order: 1) depth 2) 90th percentile of temps 3) mean Chl, 4) mean salinity and 5) substrate type 
#depth
var <- raster("Maxent_EnvironmentalFiles/Post/Depth.asc")
crs(var) <- Albers2
occ_var <- raster::extract(var, occ[,2:3])
occ_swd <- cbind(occ, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth")
#temp
var <- raster("Maxent_EnvironmentalFiles/Post/90th_%_Temp.asc")
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp")
#chl
var <- raster("Maxent_EnvironmentalFiles/Post/Mean_Chl.asc")  
crs(var) <- Albers2  
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl")
#sal
var <- raster("Maxent_EnvironmentalFiles/Post/Salinity.asc")
crs(var) <- Albers2  
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd <- cbind(occ_swd, occ_var)
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl","Salinity")
#substrate type
var <- raster("Maxent_EnvironmentalFiles/Post/Substrate.asc" )
occ_var <- raster::extract(var, occ_swd[,2:3])
occ_swd$occ_var<- occ_var
colnames(occ_swd) <- c("Species","Longitude","Latitude","Depth","90th_%_Temp","Mean_Chl","Salinity","Substrate")
#Done!
occ_swd[is.na(occ_swd)] <- -9999
write.csv(occ_swd,"OccurenceData/Absence_20172020_Albers_SWD.csv")



#Maxent analysis must then be done using the Maxent program. For this analysis, we were not able to use 
#Dismo or other R packages that run Maxent because these packages do not give you the ability to adjust 
#the parameter Tau. The Maxent terminal may be downloaded at https://biodiversityinformatics.amnh.org/open_source/maxent/



#################################################################################################
# Creating figures using Maxent output files
##################################################################################################
#Percent Importance
library(ggplot2)
library(viridis)
library(rasterVis)
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Maxent_SWD/Maxent_Output")
pre <- read.csv("Pre_Log_LQP_RM1_Tau061/Permutation_Importance_Pre.csv")
post <- read.csv("Post_Log_LQP_RM1_Tau014/Permutation_Importance_Post.csv")
phase <- factor(c('Pre-Outbreak','Pre-Outbreak','Pre-Outbreak','Pre-Outbreak','Pre-Outbreak',
           'Current','Current','Current','Current','Current'), levels = c("Pre-Outbreak","Current"))
both <- rbind(pre, post)
both <- data.frame(both$Post_Variable,phase, both$Importance)
colnames(both)<- c("Post_Variable","Phase","Importance")
ggplot(both, aes(x = Post_Variable, y = Importance,  position = Phase, fill = Phase))+
  geom_bar(stat = 'identity', position = 'dodge',color = "black")+
  scale_x_discrete(labels=c("90th Percentile\\nTemperature (C)","Depth (m)","Mean Chlorophyll\\n(mg/m-3)", "Mean Salinity\\n (PSU)","Substrate\\nCategory"))+
  scale_fill_manual(values =c( "#3d1d98","orangered2"))+
  scale_alpha_manual(values = c(0.3,1))+
  labs(x = "", y = "Permutation Importance (%)")+
  lims( y = c(0,75))+
  theme(legend.text = element_text(size =30), legend.title =element_text(size = 36)  , axis.text=element_text(size=24), axis.title = element_text(size=28))
ggsave("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/Figures/MaxentFigures/PermutationImportance_Plasma.jpeg",
       dpi = 300)
#Response Curves
setwd("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/FinalDataandScripts/Maxent_SWD/Maxent_Output")
pre_temp <- read.delim("Pre_Log_LQP_RM1_Tau061/plots/Pycnopodia_helianthoides_90th_%_Temp.dat", sep=",")
pre_depth <- read.delim("Pre_Log_LQP_RM1_Tau061/plots/Pycnopodia_helianthoides_Depth.dat", sep=",")
pre_chl<- read.delim("Pre_Log_LQP_RM1_Tau061/plots/Pycnopodia_helianthoides_Mean_Chl.dat", sep=",")
pre_sal <- read.delim("Pre_Log_LQP_RM1_Tau061/plots/Pycnopodia_helianthoides_Salinity.dat", sep=",")
pre_sub <- read.delim("Pre_Log_LQP_RM1_Tau061/plots/Pycnopodia_helianthoides_Substrate.dat", sep=",")
post_temp <- read.delim("Post_Log_LQP_RM1_Tau014/plots/Pycnopodia_helianthoides_90th_%_Temp.dat", sep=",")
post_depth <- read.delim("Post_Log_LQP_RM1_Tau014/plots/Pycnopodia_helianthoides_Depth.dat", sep=",")
post_chl<- read.delim("Post_Log_LQP_RM1_Tau014/plots/Pycnopodia_helianthoides_Mean_Chl.dat", sep=",")
post_sal <- read.delim("Post_Log_LQP_RM1_Tau014/plots/Pycnopodia_helianthoides_Salinity.dat", sep=",")
post_sub <- read.delim("Post_Log_LQP_RM1_Tau014/plots/Pycnopodia_helianthoides_Substrate.dat", sep=",")

variable <- rep(c(rep("Temp",501), rep("Depth",501),rep("Chl", 501),rep("Sal",501),rep("Sub",5)),2)
Phase <- factor(c(rep("Pre-Outbreak",2009), rep("Current",2009)), levels = c("Pre-Outbreak","Current"))
pre = rbind(pre_temp,pre_depth,pre_chl,pre_sal,pre_sub)
post = rbind(post_temp,post_depth,post_chl,post_sal,post_sub)
both <- rbind(pre, post)
both <- data.frame(both, Phase)

p1 <- ggplot(both[both$variable=="90th_%_Temp",], aes(x = x, y = y, color = Phase))+
  geom_line(size =2)+
  scale_color_manual(values = c("#3d1d98","orangered2"))+
  labs(x = "90th Perc. Temp. (C)", y = "Probability of Pycnopodia")+
  lims( y = c(0,0.9), x = c(6,25))+
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text = element_text(size = 18))
p2 <- ggplot(both[both$variable=="Depth",], aes(x = x, y = y, color = Phase))+
  geom_line(size =2)+
  scale_color_manual(values = c("#3d1d98","orangered2"))+
  labs(x = "Depth (m)", y = "")+
  lims( y = c(0,0.9), x = c(-500,0))+
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_blank())
p3 <- ggplot(both[both$variable=="Mean_Chl" ,], aes(x = x, y = y, color = Phase))+
  geom_line(size =2)+
  scale_color_manual(values = c("#3d1d98","orangered2"))+
  labs(x = "Mean Chloro. (mg/m-3)", y = "")+
  lims( y = c(0,0.9), x = c(0,37))+
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text.x = element_text(size = 17), axis.text.y = element_blank())
p4 <- ggplot(both[both$variable== "Salinity",], aes(x = x, y = y, color = Phase))+
  geom_line(size =2)+
  scale_color_manual(values = c("#3d1d98","orangered2"))+
  labs(x = "Mean Salinity (PSU)", y = "")+
  lims( y = c(0,0.9), x = c(25, 34))+
  theme(legend.position = "none", axis.title = element_text(size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_blank())
p5 <- ggplot(both[both$variable=="Substrate",], aes(x = x, y = y, position = Phase, fill = Phase))+
  geom_bar(stat = "identity",position = "dodge")+
  scale_fill_manual(values = c("#3d1d98","orangered2"))+
  labs(x = "Substrate Type", y = "")+
  lims( y = c(0,0.9))+
  theme(legend.position = "none",axis.title = element_text(size = 20), axis.text.x = element_text(size = 18), axis.text.y = element_blank())

respcurves <- multiplot(p1, p2, p3, p4, p5, cols = 5)
ggsave("C:/Users/hamiltsa/Desktop/Pycnopodia/Manuscript/Figures/MaxentFigures/ResponseCurves_Plasma.jpeg",
       plot = respcurves,
       dpi = 300)

#Maps 
#Were made using ArcGIS

