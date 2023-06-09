#### script DArT SNP genotype filtering and regional genetic offset analyses
#conducted in 'Population connectivity and genetic offset in the spawning coral "Acropora digitifera" in Western Australia (Adam et al., 2022)

##### DArT SNP genotype filtering

#required libraries
library(dartR)

#set directory
setwd("D:/Population_genetics_Adigi/Population_genetic_data/Acrodigi_popgen/DAc19_4765/Report-DAc19-4765")

#load DArT SNP data in as genlight object (in this dataset, loci associated with Symbiodinium sequences have been removed based on evalue as well as clones)
data_gl <- gl.read.dart(filename = "D:/Population_genetics_Adigi/Population_genetic_data/Acrodigi_popgen/DAc19_4765/Report-DAc19-4765/Report_DAc19-4765_2_moreOrders_SNP_mapping_2_edited_190420_removedSymbiontloci_clones.csv",topskip=6)

#Load metadata associated with DArT SNP genotype data (order samples have to be the same as in the DArT dataset)
meta<-read.csv("D:/Population_genetics_Adigi/Population_genetic_data/Acrodigi_popgen/Metadata_3runs_A_digi_inorderofSNPmapping2_020320.csv", head=T) #make sure samples are in same order as in data_gl

#preprosessing
#Link DArT dataset with metadata based on Reef system level
data_gl@pop<-meta$Location
#Drop Lalang-Garram samples because of low loci call rate
data_gl <-gl.drop.pop(data_gl,pop.list=c("Kimberley"))

#Link DArT dataset with metadata based on site level
data_gl@pop<-meta$pop

####ADD MORE INFO

#match metadata after removal of Kimberley samples and clones
meta<-meta[match(data_gl$ind.names, meta$id),]

#Recalculate loci metrics
data_gl_filtered<-gl.recalc.metrics(data_gl, v=3)

#Calculate and add sequencing coverage to loci metrics in DArT SNP genotype dataset
data_gl_filtered$other$loc.metrics$coverage<-data_gl_filtered$other$loc.metrics$AvgCountRef + data_gl_filtered$other$loc.metrics$AvgCountSnp

#Basic metrics on sequencing coverage
mean(data_gl_filtered$other$loc.metrics$coverage) #36.14553
min((data_gl_filtered$other$loc.metrics$coverage))#5.00554
max((data_gl_filtered$other$loc.metrics$coverage))#243.0508
sd(data_gl_filtered$other$loc.metrics$coverage)/sqrt(38408)#0.155447

# DArT SNP filtering
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method = "loc", threshold = 0.7, v=3) #filter based on loci call rate
data_gl_filtered <- gl.filter.callrate(data_gl_filtered, method="ind", threshold = 0.7, v=3) #filter based on individual call rate
data_gl_filtered <- gl.filter.repavg(data_gl_filtered, t=0.7,v=3) #filter based on average repeatability of alleles at a loci
data_gl_filtered <- gl.filter.monomorphs(data_gl_filtered,v=3) #remove monomorphs
data_gl_filtered <- gl.filter.hwe(data_gl_filtered, alpha = 0.05, basis = "any", bon = TRUE,v=3) #remove loci outside hardy-weinberg equilibrium
data_gl_filtered<- gl.filter.secondaries(data_gl_filtered,method="random", verbose = 3) #remove loci with more than one SNP (potentially linked)
list.match <- data_gl_filtered$loc.names[which(data_gl_filtered$other$loc.metrics$OneRatioSnp > 0.05 & data_gl_filtered$other$loc.metrics$OneRatioSnp < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef < 0.95 & data_gl_filtered$other$loc.metrics$OneRatioRef > 0.05 & data_gl_filtered$other$loc.metrics$coverage > 10)] #list SNPs with minor allele frequency >0.05 as well as SNPs with sequencing cover > 10 
data_gl_filtered<-data_gl_filtered[,match(list.match, data_gl_filtered$loc.name)] #only keep loci within the above described list
meta_filtered<-meta[match(data_gl_filtered$ind.names, meta$id),] #match metadata with filtered DArT dataset

meta_filtered$pop<-droplevels(meta_filtered$pop) #droplevels in metadata
data_gl_filtered$pop<-meta_allsamples$pop #link filtered metadata with filtered DArT SNP data

#save filtered metadata
write.csv(meta_allsamples,file="meta_data_gl_filtered_raw_1550_270520.csv")

#save filtered DArT dataset as  genlight object
saveRDS(data_gl_filtered, file="data_filtered_afterremovingKim_colones_symbiontsremoved_70_70_020520.R")
#save filtered DArT dataset as genind object
data_genind_alldata_filtereing70_70_symandclonremoved_020520<-gl2gi(data_gl_filtered, v=1)

saveRDS(data_genind_alldata_filtereing70_70_symandclonremoved_020520,file="data_genind_alldata_filtereing70_70_symandclonremoved_020520.R")

#####################################################################################################################
######################################END FILTERING DART SNP GENOTYPE DATA###########################################
########################################################################################################

#Genetic offset analyses

#used packages
library(raster)
library(adegenet)
library(gdm)
library(gradientForest)
library(reshape2)
library(ggplot2)

#prepare data
#First extract the environmental data from the coordinates of the populations
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Variables_seascape_analyses/Present")

#environmental raster data (present-day)
sstmax<- raster("sstmax_WA_nochrisnococos_Pelorus_250m_40.tif")
sstmin<-raster("sstmin_WA_nochrisnococos_Pelorus_250m_40.tif")
sstrange<-raster("sstrange_WA_nochrisnococos_Pelorus_250m_40.tif")
sstdev<-raster("sstdev_WA_nochrisnococos_Pelorus_250m_40.tif")
sstmean<-raster("sstmean_WA_nochrisnococos_Pelorus_250m_40.tif")
ssta<-raster("ssta_WA_nochrisnococos_Pelorus_250m_40.tif")
tsa<-raster("tsa_WA_nochrisnococos_Pelorus_250m_40.tif")
tsm<-raster("tsm_WA_nochrisnococos_Pelorus_250m_40.tif")
chla2<-raster("chla2_WA_nochrisnococos_Pelorus_250m_40.tif")
light<-raster("light_WA_nochrisnococos_Pelorus_250m_40.tif")
bath<-raster("bath_WA_nochrisnococos_Pelorus_250m_40.tif")
roughness<-raster("rough_WA_nochrisnococos_Pelorus_250m_40.tif")
tide<-raster("tide_WA_nochrisnococos_Pelorus_250m_40_noneg.tif")

#environmental raster data (future)
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Variables_seascape_analyses/Future/")
SSTmax_RCP26_2050<- raster("sstmax_RCP26_2050_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTmax_RCP26_2100<- raster("sstmax_RCP26_2100_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTmax_RCP85_2050<- raster("sstmax_RCP85_2050_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTmax_RCP85_2100<- raster("sstmax_RCP85_2100_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTrange_RCP26_2050<- raster("sstrange_RCP26_2050_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTrange_RCP26_2100<- raster("sstrange_RCP26_2100_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTrange_RCP85_2050<- raster("sstrange_RCP85_2050_WA_nochrisnococos_Pelorus_250m_40.tif")
SSTrange_RCP85_2100<- raster("sstrange_RCP85_2100_WA_nochrisnococos_Pelorus_250m_40.tif")

#stack all environmental and habitat data (present-day)
env<-stack(sstmax,
           sstmin,
           sstrange,
           sstmean,
           sstdev,
           ssta,
           tsa,
           tsm,
           chla2,
           bath,
           roughness,
           light,
           tide)
#crop to only the extent of Western Australia
WA_Extent<-extent(109,132,-36,-8)
env_noPel<-crop(env,WA_Extent)

#stack all environmental and habitat data (future)
env_future<-stack(SSTmax_RCP26_2050,
                  SSTmax_RCP26_2100,
                  SSTmax_RCP85_2050,
                  SSTmax_RCP85_2100,
                  SSTrange_RCP26_2050,
                  SSTrange_RCP26_2100,
                  SSTrange_RCP85_2050,
                  SSTrange_RCP85_2100)

env_noPel_future<-crop(env_future,WA_Extent)

#rename
Envstack_present_selected <- env
Envstack_future_selected <- env_future
#georeference
WGS84<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
crs(Envstack_present_selected) <-WGS84

#subset coordinates of sites
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/radiator_bayescan_20200502@0941")
meta_filtered<-read.csv("meta_allsamples_neutral_250520.csv")
unique(meta_neutral$pop)
meta_neutral_coordinates<- meta_neutral[c(4,8,9)] #sitename/longitude/latitude
unique_coord_sitename <- unique(meta_neutral_coordinates) # subset unique coordinates
coordinates(unique_coord)<- ~Longitude+Latitude
crs(unique_coord)<-WGS84

#save unique site coordinates
setwd("/Volumes/ArneCurtin/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/")
write.csv(unique_coord_sitename,file="coordinates_sites_140820.csv")

#extract variable values at site locations
data_sites <- data.frame(coordinates(unique_coord), extract(Envstack_present_selected,unique_coord))
data_sites$Sites<- unique_coord_sitename[,1]
#save data
setwd("/Volumes/ArneCurtin/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/")
write.csv(data_sites,file="Extracted_values_sites_coordinates_raw_140820.csv")

#Evaluating correlation of just values on site locations
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/")
Site_env_values<-read.csv("Extracted_values_sites_coordinates_140820.csv")
cor_site_values<-cor(Site_env_values, method = "pearson")
#save correlation matrix
write.csv(cor_site_values,file="Correlation_site_values_raw_271020.csv")

####calculate minor allele frequency (MAF) of unique outlier loci identified by Bayescenv
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bath")
outlier_loci<-read.csv("Outlier_afterbayescenv_noPelorus_allvar_noduplicates_130121.csv")
outlier_loci<-outlier_loci[5] #only for Gradientforest

#subset filtered DArT dataset, only keeping the outlier loci sequences
data_gl_filtered_outlier<-data_gl_filtered[,match(outlier_loci$SNP, data_gl_filtered$loc.names)]

####genetic offset analyses are only conducted on WA reefs
data_gl_filtered_outlier$pop<-meta_filtered$pop

#select unique sites
populations<-unique(data_gl_filtered_outlier$pop)
populations<- as.character(populations)

#MAF  loci per location
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")

populations<-unique(data_gl_filtered_outlier$pop)

for (i in 1:26){
  pop<-populations[i]
  data_gl_filtered_outlier_pop <-gl.keep.pop(data_gl_filtered_outlier,pop.list = pop)

  data_genind_filtered_outlier_pop <-gl2gi(data_gl_filtered_outlier_pop)
  MAF <- minorAllele(data_genind_filtered_outlier_pop)
  min(MAF)
  MAF_pop<-mean(MAF)
  t<-nAll(data_genind_filtered_outlier_pop)
  r<-summary( data_genind_filtered_outlier_pop )
  write.csv(MAF_pop,file=paste(pop,"_outlier_mean_MAF.csv",sep=""))
  write.csv(t,file=paste(pop,"_polymorphism.csv",sep=""))
  write.csv(r,file=paste(pop,"_polymorphism_het_perlocus.csv",sep=""))
  write.csv(MAF,file=paste(pop,"_outlier_MAF.csv",sep=""))
}


####Gradientforest analysis
#prep
#check polymorphism >5 sites polymorphic by concatenating all polymorphisms checks for the outlier loci
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/Polymorphisms/")
poly_files= list.files(pattern="*.csv")
poly<-as.data.frame(read.csv(poly_files[1],header=T)[,1])
df <-do.call(cbind,lapply(poly_files,function(fn)read.csv(fn,header=T)[,2]))
df2<-cbind(poly,df)
#save
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
write.csv(df2,file="Allsites_65SNPs_bayescenv_1_26_polymorphisms_noPelorus_070921.csv")

#concatenate MAF values for each unique outlier loci for all sites
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/MAF/")

MAF_files= list.files(pattern="*.csv")
MAF<-as.data.frame(read.csv(MAF_files[1],header=T)[,1])
MAF_df <-do.call(cbind,lapply(MAF_files,function(fn)read.csv(fn,header=T)[,2]))
MAF_df2<-cbind(MAF,MAF_df)

#save
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
write.csv(MAF_df2,file="Allsites_65SNPs_bayescenv_1_26_MAF_raw_noPelorus_070921.csv")

#modelling using gradientforest
#input file (first columns are site names, lon, lat, all variable data, all uniquye outlier loci (>5 site polymorphism) MAF values for each site)
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
Adigi_input_gradientforest<-read.csv("actualvalues_sites_outliers_noPelfromstart_1561loci_64SNPs_MAF_070921.csv")
nrow(Adigi_input_gradientforest)
head(Adigi_input_gradientforest)
ncol(Adigi_input_gradientforest)
# get variable data
envGF <- Adigi_input_gradientforest[,4:12]
head(envGF)
#change names
colnames(envGF)<-c("SSTmax","SSTrange","SSTA","TSM","Bathymetry","Roughness","Light","Tidal_height","Chla2")
colnames(envGF)
# get SNP data
SNPs_Adigi<-Adigi_input_gradientforest[,13:76]
ncol(SNPs_ref)
#account for correlations (see Fitzpatrick et al., 2015 script)
maxLevel <- log2(0.368*nrow(envGF)/2)

#model
gfAdigi <- gradientForest(cbind(envGF, SNPs_Adigi), predictor.vars=colnames(envGF),
                        response.vars=colnames(SNPs_Adigi), ntree=2000, 
                        maxLevel=maxLevel, trace=T, corr.threshold=0.8)
#save model
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
saveRDS(gfAdigi,file="Gradientforest_ntree_64SNPs_2000_corrthreshold_80_070921.R")

summary(gfAdigi)
imp.SNPs<-gfRef$result #extract data of SNPs with positive Rsquared
mean(imp.SNPs) #mean Rsquared
sd(imp.SNPs)#deviation
#save
write.csv(imp.SNPs,file="Mostimp_64SNPs_R2value_280221.csv")

#plot variable importance
type = "O"
plot(gfRef,plot.type="O")
dev.copy2pdf(file="Importanceplots_gradientforest__180121.pdf")

#plot SNP cumulative plot of the top 5 most responsive SNPs for each predictor
type = "C" #species cumulqtive plot, identifies the top 5 most responsive SNPs for each predictor
most_important <- names(importance(gfAdigi))
plot(gfAdigi, plot.type = "C", imp.vars = most_important,show.overall = F,common.scale=T, 
     legend = T, leg.posn = "topleft",leg.nspecies = 5, 
     cex.lab = 1, cex.legend = 1,cex.axis = 1, 
     line.ylab = 0.9, par.args = list(mgp = c(1.5,0.5, 0), 
                                      mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.copy2pdf(file="Importanceplots_allimpSNPs_Chla2_cumulative_withlegend_180121.pdf")

#plot SNP cumulative plot of the top 5 most responsive SNPs for each predictor
plot(gfAdigi, plot.type = "C", imp.vars = most_important,show.overall = F,common.scale=T, 
     legend = F, leg.posn = "topleft",leg.nspecies = 5, 
     cex.lab = 1, cex.legend = 1,cex.axis = 1, 
     line.ylab = 0.9, par.args = list(mgp = c(1.5,0.5, 0), 
                                      mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.copy2pdf(file="Importanceplots_allimpSNPs_Chla2_cumulative_nolegend_180121.pdf")

#select most imporant variables
imp.vars <- names(importance(gfAdigi))

#predict under present-day conditions of 50 buffer zones around sites (nocorrelative variables only)
spts_present <- rasterToPoints(env_nocorrelative, spatial = F)
saveRDS(spts_present,file="spts_present_noPelorus_180121.R")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bayescenv_noPelorus_080121/gradientforest_130121/")
spts_present_noNA = as.data.frame(spts_present[complete.cases(spts_present), ]) #remove NA's (should be zero)
colnames(spts_present_noNA)<-c("Longitude","Latitude","Bathymetry","Chla2",
                               "Light",
                               "Roughness",
                               "SSTA",
                               "SSTmax",
                               "SSTrange",
                               "Tidal_height",
                               "TSM")

#Map genetic similarity using only the most important variables
Trns_grid <- cbind(spts_present_noNA[,c("Longitude","Latitude")], predict(gfAdigi,spts_present_noNA[,imp.vars])) 
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
saveRDS(Trns_grid,file="Trns_grid_ninevar_noPelorus_64SNPs_alphabet_2802121.R")

#plot as PCA
PCs <- prcomp(Trns_grid[, imp.vars])
## set up a colour palette for the mapping
a1 <- PCs$x[,1]
a2 <- PCs$x[,2]
a3 <- PCs$x[,3]

r <- a1+a2
g <- -a2
b <- a3+a2-a1
r <- (r-min(r)) / (max(r)-min(r)) * 255
g <- (g-min(g)) / (max(g)-min(g)) * 255
b <- (b-min(b)) / (max(b)-min(b)) * 255

setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bayescenv_noPelorus_080121/gradientfores_280221/")
nrow(Trns_grid)
plot(Trns_grid[, c("Longitude", "Latitude")], pch = ".",cex = 3, asp = 1, col = rgb(r, g, b, max = 255))
dev.copy2pdf(file="Present_genetic_similarity_070921_noPel_Chla_used_64SNPs.pdf")

nvs <- dim(PCs$rotation)[1]  # number of variables
# choose predictor vectors to show five most important ones

vec <- c("Tidal_height", 
         "SSTmax",
         "SSTrange",
         "SSTA",
         "Bathymetry")

lv <- length(vec)
vind <- rownames(PCs$rotation) %in% vec

# choose a scaling factor to plot the vectors over the grid
scal <- 10 
xrng <- range(PCs$x[,1], PCs$rotation[,1]/scal)*1.1
yrng <- range(PCs$x[,2], PCs$rotation[,2]/scal)*1.1
plot((PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=4, col=rgb(r,g,b, max = 255))
# plot the chosen predictors as arrows
arrows(rep(0,lv), rep(0,lv), PCs$rotation[vec,1]/scal, PCs$rotation[vec,2]/scal, length = 0.0625) 
scal <- 10 
jit <- 0.015
text(PCs$rotation[vec,1]/scal+jit*sign(PCs$rotation[vec,1]), PCs$rotation[vec,2]/scal+jit*sign(PCs$rotation[vec,2]), labels = vec)

# predict the PCs for the transformed site data and plot on PCA
Trns_site <- predict(gfAdigi) 
PCsites <- predict(PCs,Trns_site[,imp.vars])  

# plot all the sites as points on the biplot
points(PCsites[,1:2])

PCsites<-cbind(Roughness_MAF[,c(1)],PCsites)
write.csv(PCsites,file="PCAsites_allsites_Present_candidateSNPs_64SNPs_070921.csv")

#use Fitzpatrick script to predict genetic offset
spts_present_noNA2<-spts_present_noNA[,-c(1,2)]

#Future data
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP26_2050")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP26_2100")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP85_2050")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP85_2100")

predictors_future <- list.files(pattern="*.tif")
Envstack_future_selected <- stack(predictors_future)
names(Envstack_future_selected)
WGS84<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0")
crs(Envstack_future_selected) <-WGS84
spts_future <- rasterToPoints(Envstack_future_selected, spatial = F)
head(spts_future)
saveRDS(spts_future, file="spts_future_RCP26_2050_noPelorus_280121.R")

setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bayescenv_noPelorus_080121/gradientforest_130121/")
spts_future<-readRDS("spts_future_RCP26_2050_noPelorus_280121.R")
spts_future<-readRDS("spts_future_RCP85_2050_noPelorus_280121.R")
spts_future<-readRDS("spts_future_RCP85_2100_noPelorus_280121.R")
spts_future<-readRDS("spts_future_RCP26_2100_noPelorus_280121.R")

spts_future<-as.data.frame(spts_future)
spts_future_noNA = as.data.frame(spts_future[complete.cases(spts_future), ]) #remove NA's
nrow(spts_future_noNA)
head(spts_future_noNA)
ncol(spts_future_noNA)
colnames(spts_future_noNA)<-c("Longitude","Latitude","Bathymetry","Chla2",
                              "Light",
                              "Roughness",
                              "SSTA",
                              "SSTmax",
                              "SSTrange",
                              "Tidal_height",
                              "TSM")

# first transform FUTURE env. variables
Trns_grid_future <- cbind(spts_future_noNA[,c("Longitude","Latitude")], predict(gfAdigi,spts_future_noNA[,imp.vars])) 

#estimate predictedgenetic offset
genOffset <- sqrt((Trns_grid_future[,3]-Trns_grid[,3])^2+(Trns_grid_future[,4]-Trns_grid[,4])^2
                  +(Trns_grid_future[,5]-Trns_grid[,5])^2+(Trns_grid_future[,6]-Trns_grid[,6])^2
                  +(Trns_grid_future[,7]-Trns_grid[,7])^2+(Trns_grid_future[,8]-Trns_grid[,8])^2
                  +(Trns_grid_future[,9]-Trns_grid[,9])^2+(Trns_grid_future[,10]-Trns_grid[,10])^2
                  +(Trns_grid_future[,11]-Trns_grid[,11])^2)

genOffset<-as.data.frame(genOffset)
test<-cbind(Trns_grid_future[, c("Longitude", "Latitude")],genOffset)

#rasterise data
test_raster_future<-rasterFromXYZ(test)
plot(test_raster_future)
new_extent<-extent(110,132,-24,-11)
test_raster_noPelorus<-crop(test_raster_future,new_extent)
plot(test_raster_noPelorus)

writeRaster(test_raster_future, file="Geneticoffset_RCP26_2050_64SNPs_bayescenv_noPelorus_recheck_onlyWA_070921",format="GTiff")
writeRaster(test_raster_future, file="Geneticoffset_RCP85_2050_64SNPs_bayescenv_noPelorus_recheck_onlyWA_070921",format="GTiff")
writeRaster(test_raster_future, file="Geneticoffset_RCP85_2100_64SNPs_bayescenv_noPelorus_recheck_onlyWA_070921",format="GTiff")
writeRaster(test_raster_future, file="Geneticoffset_RCP26_2100_64SNPs_bayescenv_noPelorus_recheck_onlyWA_070921",format="GTiff")


###### Generalised dissimilarity model
library(gdm)

#input data is site locations, lat, long and Fst values unique outlier loci
Fst_outliers_coord->cbind(Adigi_input_gradientforest$Sites,Fst_outliers)
write.csv(Fst_outliers_coord,file="Fst_outliers_withsitenames_080322.csv")
#convert to gdm format
site_pair_format_tide_sstmax<-formatsitepair(Fst_outliers_coord,3,siteColumn="Sites",XColumn = "Longitude",YColumn = "Latitude",predData =envGDM_tide_sstmax) #finetuned input data as gdm with all variable only produced 2 contributing variables, SSTmax and tidal height
write.csv(site_pair_format_tide_sstmax,file="GDM_site_pair_table_65outliers_envnamesconsistent_onlytide_sstmax_140322.csv")

#normalise fst within 0-1 range (see script Fitzpatrick et al., 2015)
norm<-function(X) {(X-min(X))/(max(X)-min(X))}
norm_data<-norm(X=site_pair_format_tide_sstmax$distance)
write.csv(norm_data,file="Fst_normalised_65outliers_090322.csv")

#GDM fitting
Fst_sitepair<-read.csv("GDM_site_pair_table_65outliers_Fst_normalised_edited_080322.csv")
site_pair_format_tide_sstmax$distance<-Fst_sitepair$distance
saveRDS(site_pair_format,"GDM_modelfitting_inputfile_65outliers_090322.R")

#GDM only geographic space/tide/sstmax
gdm<-gdm(data=site_pair_format_tide_sstmax,geo=T)
summary(gdm) #some coefficients is variable importance for each variable
saveRDS(gdm,file="GDM_model_65outliers_onlytide_sstmax_140322.R")
length(gdm$predictors)
#plot variable contribution
plot(gdm,plot.layout=c(2,3))
#extract splines
splines<-isplineExtract(gdm1)
norm_splines<-norm(X=splines$y) #normalisation of splines --> used for plotting
write.csv(norm_splines,file="splines_gdm_onlytide_sstmax_140322.csv")

#example plot splines from variable sstmax
act_data<-as.data.frame(splines$x)
act_geo<-act_data$SSTmax #actual values

splines_all<-as.data.frame(norm_splines)
splin_geo<-splines_all$SSTmax

plot_spline2<-as.data.frame(cbind(act_geo,splin_geo))
class(plot_spline)
library(ggplot2)
t<- ggplot(plot_spline2, aes(x=act_geo,y=splin_geo))+geom_line(size=1)+xlab("SST max")+ylab("Relative importance")+theme_bw()+theme(legend.position = "none")+
  theme(axis.title.y=element_text(face="bold", color = "black", size=14),
        axis.title.x =element_text(face="bold", color = "black", size=12),
        axis.text.y = element_text(face="bold",color = "black", size=10,angle = 90, hjust = 1),
        axis.text.x = element_text(face="bold", color = "black",size=12))+
  coord_cartesian(ylim=c(0,1))+scale_y_continuous(breaks=seq(0,1,by=0.2))
t

### Predict across 50 km buffer zone using gd,
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bayescenv_noPelorus_080121/gradientforest_130121/")

spts_present<-readRDS("spts_present_noPelorus_180121.R")
#plot(Envstack_present_selected[[1]])

head(spts_present)
nrow(spts_present)
spts_present<-as.data.frame(spts_present)
#freq(spts_present)
spts_present_noNA = as.data.frame(spts_present[complete.cases(spts_present), ]) #remove NA's

colnames(spts_present_noNA)<-c("Longitude","Latitude","Bathymetry","Chla2",
                               "Light",
                               "Roughness",
                               "SSTA",
                               "SSTmax",
                               "SSTrange",
                               "Tidal_height",
                               "TSM")
#colnames(Rough)<-c("Longitude","Latitude","Bathymetry","Chla2",
#                               "Light",
#                               "Roughness",
#                               "SSTA",
#                               "SSTmax",
#                               "SSTrange",
#                               "Tidal_height",
#                               "TSM")

present_onlytide_Sstmax<-spts_present_noNA[,c(1,2,8,10)]
#transform new data
#using onlytide_Sstmax gdm model
#buffer
Trns_grid_gdm_present_onlytide_sstmax<-gdm.transform(model=gdm1,data=present_onlytide_Sstmax)

#sites
Trns_grid_gdm_sites_onlytide_sstmax<-gdm.transform(model=gdm1,data=envGDM_tide_sstmax[,-1])

#saveRDS(Trns_grid,file="Trns_grid_ninevar_noPelorus_180121.R")
#setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/bayescenv_noPelorus_080121/gradientforest_130121/")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")

saveRDS(Trns_grid_gdm_present_onlytide_sstmax,file="GDM_Trns_grid_onlytide_sstmax_noPelorus_65SNPs_present_170322.R")
saveRDS(Trns_grid_gdm_sites_onlytide_sstmax,file="GDM_Trns_grid_onlytide_sstmax_noPelorus_65SNPs_sites_170322.R")


#Future data
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP26_2050")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP26_2100")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP85_2050")
setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Bayesenv_analysis_130820/gradientForest/Future/RCP85_2100")

predictors_future <- list.files(pattern="*.tif")
#only tide and sstmax
#future rasters
Envstack_future_selected <- stack(predictors_future)
names(Envstack_future_selected)<-names(Envstack_present_selected) #same names

#predict genetic offset using gdm
Trns_grid_gdm_RCP26_2100<-predict.gdm(gdm1,Envstack_present_selected,time=T,predRasts = Envstack_future_selected)

new_extent<-extent(11,131,-30,-10)
Future_gdm_WA<-crop(Trns_grid_gdm_RCP26_2100,new_extent)
plot(Future_gdm_WA)

setwd("G:/PhD_project/Population_genetic_data/Acrodigi_popgen/050520/bayescenv_gen_file_270520/Gradientforest_070921/")
writeRaster(Future_gdm_WA, file="Geneticoffset_GDM_only_tide_sstmax_RCP85_2100_65SNPs_noPelorus_predict_onlyWA_170322",format="GTiff")

######END Script########
