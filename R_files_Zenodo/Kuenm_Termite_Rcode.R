#Load required packages
library(spocc)
library(spThin)
library(dismo)
library(rgeos)
library(ENMeval)
library(dplyr)
library(wallace)
library(ggplot2)
library(rgdal)
library(usdm)
library(kuenm)

#Set working directory which possesses individual folders of each species
  #Folder setup follows Figure 2 from Cobos et al. 2019 for kuenm r package
#Each species folder must contain the raw (unthinned) occurrences for each species with 'OccID', 'species', 'longitude', and 'latitude' in that order
  #Spatial thinning ('spThin' package) requires each occurrence posses a number ID associated with it i.e. 'OccID' 
#Furthermore, each species folder must contain the maxent.jar file to run the analysis
#Finally have a folder containing all 19 bioclimatic variables in the resolution you need for the analysis

#set working directory to your species folders
setwd("DIRECTORY/Species")

#assign 'wd' as your working directory
wd<-setwd("DIRECTORY/Species")

#create a for loop which will go through each folder building in order to start the kuenm package
for (i in list.dirs(wd)[-1]) {

#load occurences from csv file 
setwd(i)
occs <- list.files(pattern = "\\\\.csv$")
occs<-read.csv(occs)

#apply a buffered bounding box around your occurrences, will need for cropping envs later 
bgExt <- penvs_bgExtent(
  occs = occs,
  bgSel = "bounding box",
  bgBuf = 5)

#spatially thin occurrences according to envs resolution (in this case 5km)
output <- spThin::thin(occs, 'longitude', 'latitude', 'species', thin.par = 5, reps = 100, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
occs <- occs[as.numeric(rownames(maxThin)),] 
occs<-occs[2:4]
occs <- occs[!duplicated(occs),]
str(occs)

#Load enviromental data and assign categorical variables as factors
envs.files <- list.files(path=("DIRECTORY/ENVS"), pattern='tif', full.names=TRUE)
envs <- raster::stack(envs.files)

#Crop environmental layers by extent 
envs.bg <- raster::crop(envs, bgExt)
plot(envs.bg)

#generate 25,000 random background points 
bg <- dismo::randomPoints(envs.bg[[4]], n = 25000) %>% as.data.frame()

#Extract cropped environmental values from background points
bg.z <- cbind(bg, raster::extract(envs.bg, bg))
bg.z<-bg.z[3:22]

#Calculate correlation among variables using VIFcor 
r<-bg.z# calculates vif for the variables in r
v1 <- vifcor(r, th=0.9) # identify collinear variables that should be excluded
v1
re1 <- exclude(r,v1) # exclude the collinear variables that were identified in
# the previous step
re1
v2 <- vifstep(r, th=10) # identify collinear variables that should be excluded
v2
re2 <- exclude(r, v2) # exclude the collinear variables that were identified in
# the previous step
re2
re3 <- exclude(r) # first, vifstep is called
re3

#subset the env variables which are uncorrelated to use in further analysis 
envs.bg_t<-subset(envs.bg, colnames(re3))

#write rasters to M_variable and G_variable folder 
#Since I am not testing out different environmental variables, both folders will have Set 1
setwd(i)
dir.create("M_variables/Set 1", recursive = TRUE)
dir.create("G_variables/Set 1", recursive = TRUE)
mainDir <- i
subDir <- "/M_variables/Set 1"
setwd(file.path(mainDir, subDir))
writeRaster(envs.bg_t, names(envs.bg_t), bylayer=T, format='ascii', overwrite = T)

setwd(i)
mainDir <- i
subDir <- "/G_variables/Set 1"
setwd(file.path(mainDir, subDir))
writeRaster(envs.bg_t, names(envs.bg_t), bylayer=T, format='ascii', overwrite = T)

#Get random 50/50 partition of occurrence data 
occs_cb1<-occs
names(occs_cb1)[names(occs_cb1) == 'longitude'] <- 'x'
names(occs_cb1)[names(occs_cb1) == 'latitude'] <- 'y'
occs_cb1<-occs_cb1[,2:3]

cb1 <- get.checkerboard1(occs_cb1, envs.bg, bg, aggregation.factor=5)

#set working directory, merge occurrences with partition scheme 
#create a csv of of all occurrences
setwd(i)
cb1$occs.grp
occs_part<-cbind(occs, cb1$occs.grp)
occs_joint<-occs_part[,1:3]
write.csv(occs_joint,"Sp_joint.csv", row.names=FALSE)

#filter training dataset by first group and save it 
Sp_train<-filter(occs_part, cb1$occs.grp == 1)
Sp_train<-Sp_train[,1:3]
write.csv(Sp_train, "Sp_train.csv", row.names=FALSE)

#filter testing dataset by second group and save it 
Sp_test<-filter(occs_part, cb1$occs.grp == 2)
Sp_test<-Sp_test[,1:3]
write.csv(Sp_test, "Sp_test.csv", row.names=FALSE)

#End loop 
}

#set working directory to each individual species within your species folder
wd<-setwd("DIRECTORY/SPECIES_1")

#Run kuenm to generate Maxent candidate models that will be written in subdirectories
occ_joint <- "Sp_joint.csv"
occ_tra <- "Sp_train.csv"
M_var_dir <- "M_variables"
batch_cal <- "Candidate_models"
out_dir <- "Candidate_Models"
reg_mult <- c(seq(1, 5, 1))
f_clas <- c("l", "q","lq","h","lqh")
args <- c("maximumbackground=25000","togglelayertype=class")
maxent_path <- wd
wait <- FALSE
run <- TRUE

kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
          #out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, args = args,
          #maxent.path = maxent_path, wait = wait, run = run)

#Evaluate model performance based on statistical significance (partial ROC), omission rate
occ_test <- "Sp_test.csv"
out_eval <- "Calibration_results"
threshold <- 10
rand_percent <- 50
iterations <- 500
kept <- TRUE
selection <- "OR_AICc"
paral_proc <- FALSE # make this true to perform pROC calculations in parallel, recommended
# only if a powerfull computer is used (see function's help)
# Note, some of the variables used here as arguments were already created for previous function

cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test, batch = batch_cal,
                        out.eval = out_eval, threshold = threshold, rand.percent = rand_percent, iterations = iterations,
                        kept = kept, selection = selection, parallel.proc = paral_proc)

#creating final models and, if needed, transferring them to other areas or scenarios
batch_fin <- "Final_models"
mod_dir <- "Final_Models"
rep_n <- 10
rep_type <- "Bootstrap"
jackknife <- FALSE
out_format <- "cloglog"
project <- TRUE
G_var_dir <- "G_variables"
ext_type <- "ext_clam"
write_mess <- FALSE
write_clamp <- FALSE
wait1 <- FALSE
run1 <- TRUE
args <- NULL # e.g., "maximumbackground=20000" for increasing the number of pixels in the bacground or
# "outputgrids=false" which avoids writing grids of replicated models and only writes the 
# summary of them (e.g., average, median, etc.) when rep.n > 1
# note that some arguments are fixed in the function and should not be changed
# Again, some of the variables used here as arguments were already created for previous functions

kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin, rep.n = rep_n,
          rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir, out.format = out_format, project = project,
          G.var.dir = G_var_dir, ext.type = ext_type, write.mess = write_mess, write.clamp = write_clamp, 
          maxent.path = maxent_path, args = args, wait = wait1, run = run1)

