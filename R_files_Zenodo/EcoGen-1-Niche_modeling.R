###########################################################################
########      STATISTICAL ANALYSES on niche inferences  (part 1)    #######
###########################################################################

########################
###  Load libraries  ###

library(devtools)
library(gtools)
library(raster)
library(dplyr)
library(rgdal)
library(ecospat)
library(vegan)
library(agricolae)
library(sdmpredictors)
library(ggplot2)
library(reshape2)
library(multcompView)


###########################################
###       Set working directories       ###

workDir = "PATH/to/workingDirectory/"
dataDir = "PATH/to/data/"
outDir = "PATH/to/outputs/"
setwd(workDir) 






########################################################################
##### ------------------------------------------------------------ #####
#####     0 - PRELIMINARY datasets preparation  (BioClim layers)   #####
##### ------------------------------------------------------------ #####
########################################################################


##  Following script is modified from 'user_script_Nsp_1A.R' available with the R-package 'ecospat' (Broennimann et al. 2012)  ##


###########################
###  Load Bioclim maps  ###

# Climatic maps
# CHELSA maps can be downloaded here: http://chelsa-climate.org/downloads/
# CGIAR maps for 'soil-water-content' can be downloaded here: https://cgiarcsi.community/data/global-high-resolution-soil-water-balance/


# CHELSA #
map_ch = list.files(path=paste0(dataDir,"01_maps/"), pattern="*.tif", full.names=T) 
map_ch = map_ch[mixedorder(map_ch)] # make sure files are alphanumerically ordered (1 to 19 following BioClim variables standardized names #)
for(i in 1:length(map_ch)){
  assign(paste0("bio",i), raster(map_ch[i]))}

# CGIAR #
map_cg = list.files(path=paste0(dataDir,"01_maps/"), pattern="*.adf", full.names=T, recursive=T)
map_cg = map_cg[sapply(map_cg,file.size) > 20000000] #keep files with size > 20Mb (12 files corresponding to monthly data should be retained)
map_cg = map_cg[mixedorder(map_cg)] # make sure files are alphanumerically ordered (by months, ie. from 1 to 12)
for(i in 1:length(map_cg)){
  assign(paste0("swc",i), raster(map_cg[i]))}
# extract mean values of monthly data to create 'mean annual swc' variable
stack_cgiar = stack(mget(paste0("swc",1:12)))
swc = calc(stack_cgiar, fun = mean, na.rm=T)
writeRaster(swc, paste0(dataDir,"01_maps/swc.tif"), options=c("COMPRESS=LZW"), overwrite=T)



#################################################################
###  Crop bioclim maps to the extent of background study area ###

# Load common background study map (previously defined and created on ArcGIS or another GIS software)
sh = shapefile(paste0(dataDir,"bg_f_buffer_Intersect_Merge.shp"))

# Create empty mask (template onto which maps will be cropped)
mask = rasterize(sh, bio1)
mask = crop(mask, extent(sh))
writeRaster(mask, paste0(outDir,"00_crop_maps/mask.tif"), options=c("COMPRESS=LZW"), overwrite=T)

# Crop bioclim maps
name.maps = c(paste0("bio",1:19),"swc")

for (i in name.maps){

	map <- get(i)
	mapc <-crop(map,extent(mask))
	mapcm <- mask(mapc,mask)
	
	writeRaster(mapcm,paste0(outDir,"_0_crop_maps/",i,".tif")
			   ,overwrite=T,options=c("COMPRESS=LZW"))
}


################################################################
###  Random sampling from study background for PCA analyses  ###


# random sample of 1M raster cells ("points") from background
smask = sampleRandom(mask,1200000,xy=T,sp=T,na.rm=T,cells=T) # sample size is larger to take into account potential missing data (~10min to run)
# extract climatic values of these 1M raster cells
smask = values(stack_bioclim20) # ! need to allocate vector of size 2.9 Gb
smask = rasterToPoints(stack_bioclim20) # ! need to allocate verctor of size 1.1 Gb

for(i in name.maps){
  assign(i, raster(paste0(outDir,"_0_crop_maps/",i,".tif"))) } #load cropped maps
stack_maps = stack(mget(name.maps))
 
#xy coordinates of the cells assessed from clim_chelsa (=extrapolated to clim_cgiar), identical to smask at 0.0001
clim.xy = xyFromCell(stack_maps,smask$cell)
clim = na.exclude(data.frame(cbind(clim.xy, a_clim), stringsAsFactors=F))
clim = clim[sample(nrow(clim), 1000000),] #resample to get final dataset of 1M raster cells




########################################################################
##### ------------------------------------------------------------ #####
#####     0 - PRELIMINARY datasets preparation  (species layer)    #####
##### ------------------------------------------------------------ #####
########################################################################


### --- LOAD the Aegilops dataset --- ###
#Occurences data
occ = read.table(paste0(dataDir,"MCDU_GBIFoccurrences.txt"),h=T,sep="\\t")
#Genotyped accessions
sp = na.exclude(read.table(paste0(dataDir,"MCDU_GenotypedSamples.txt"),h=T,sep="\\t"))


## --- MERGE THE TWO DATASETS (based on "sp" colnames) --- ##
aeg = data.frame(rbind(sp[,c("ID","X","Y")], setNames(occ[,c("ID","Lon","Lat")], names(sp[,c("ID","X","Y")])))) 
aeg$sp = paste0(substr(aeg$ID,1,1),tolower(substr(aeg$ID,2,2)))
aeg = aeg[!duplicated(aeg[,c("sp","X","Y")]),c("sp","ID","X","Y")] # remove duplicated accessions from the merged dataset

# associate corresponding bioclim values
aeg_bioclim = extract(stack_maps,aeg[,c("X","Y")],df=T,cellnumbers=T)
aeg.xy = xyFromCell(stack_maps,aeg_bioclim$cells)
aegc = cbind(aeg[,c("sp","ID")], aeg.xy, aeg_bioclim[,-c(1,2)])
aegc = na.exclude(aegc[!duplicated(aegc[,c("sp","x","y")]),])


## --- FINAL DATASET with selected BioClim variable --- ##
biovar = c("T","TdV","isoT","TS","Twarm","Tcold","TaV","Twet3","Tdry3","Twarm3","Tcold3","P","Pwet","Pdry","PS","Pwet3","Pdry3","Pwarm3","Pcold3","swc")
nvar = length(biovar)
biodata = rbind(aegc[,name.maps],clim[,name.maps])
colnames(biodata) = biovar





#########################################################################
##### ------------------------------------------------------------- #####
#####                      1a-  NICHE MODELING                      #####
#####                  Individual niche modeling                    #####
##### ------------------------------------------------------------- #####
#########################################################################


# RAM upper limit for niche inferences
memory.limit(size=18000)
# vector of weight, 0 for the occurences, 1 for the sites of the study area
w = c(rep(0,nrow(aegc)),rep(1,nrow(clim)))
#number of iterations for equivalency tests
it = 100
#resolution of the gridding of the environmental space
R = 100

## PCA calibrated on all the sites of the study area
pca.cal = dudi.pca(biodata_f, row.w=w, center=T, scale=T, scannf=F, nf=2)
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
pca.var = head(pca.cal$eig/sum(pca.cal$eig)*100) # % variance explained by the 6 first PCs



#####################################################################
###  -----------      Individual niche modeling      -----------  ###


pdf(paste0(outDir,"Aegilops_ecoNiches.pdf"))
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)

i=1
for(i in 1:length(mcdu)){
  name.aeg = mcdu.names[i]
  row.aeg = grep(mcdu[i],aegc$sp)
  #predict scores on PCs
  scores.aeg = pca.cal$li[row.aeg,]
  #calculation of occurence density
  z = ecospat.grid.clim.dyn(scores.clim,scores.clim,scores.aeg,R)
  assign(paste0("z.",mcdu[i]),z)
  ecospat.plot.niche(z,title=name.aeg,name.axis1="PC1",name.axis2="PC2")
}

dev.off()



#####################################################################
###  Pairwise comparisons and tests of individual modeled niches  ###

ncombn = nrow(t(combn(mcdu,m=2)))
mcdu.res = matrix(nrow=ncombn,ncol=6)
mcdu.res[,1:2] = t(combn(mcdu,m=2))
colnames(mcdu.res) = c("spa","spb","overlap","sim_ab","sim_a","sim_b")

i=1
for(i in 1:ncombn){
  name.spa = mcdu.res[i, "spa"]
  name.spb = mcdu.res[i, "spb"]
  row.spa = grep(name.spa,aegc$sp)
  row.spb = grep(name.spb,aegc$sp)
  scores.spa = pca.cal$li[row.spa, ]
  scores.spb = pca.cal$li[row.spb, ]
  za = get(paste0("z.", name.spa))
  zb = get(paste0("z.", name.spb))
  
  ## niche overlap ##
  # if cor=T, niche overlap is calculated on z$z.cor (ie. value weighted by background area)
  mcdu.res[i,"overlap"] = unlist(ecospat.niche.overlap(za,zb,cor=F)[1]) 
  
  ## niche similarity tests (according to warren et al. 2008)
  #only za is randomly shifted in the background
  sim.a = ecospat.niche.similarity.test(zb,za,rep=100,alternative="greater",rand.type = 2)
  mcdu.res[i,"sim_a"] = sim.a$p.D
  #only zb is randomly shifted in the background
  sim.b = ecospat.niche.similarity.test(za,zb,rep=100,alternative="greater",rand.type = 2)
  mcdu.res[i,"sim_b"] = sim.b$p.D
  #both za and zb are randomly shifted in the background
  sim.ab = ecospat.niche.similarity.test(za,zb,rep=100,alternative="greater",rand.type = 1)
  mcdu.res[i,"sim_ab"] = sim.ab$p.D
  
  # plots
  pdf(file=paste0(outDir, name.spa," x ",name.spb,".pdf"))
  layout(matrix(c(1,1,2,2,
                  1,1,2,2,
                  3,3,4,5,
                  3,3,6,7)
                ,nrow=4,ncol=4,byrow=T))
  
  par(mar=c(4,4,1,0.5))
  ecospat.plot.niche(za, title=name.spa, name.axis1="PC1", name.axis2="PC2")
  ecospat.plot.niche(zb, title=name.spb, name.axis1="PC1", name.axis2="PC2")
  ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
  
  par(mar=c(0.5,0.5,0.5,0.5))
  plot.new()
  text(0.5,0.5,paste0("niche overlap\\nD = ",round(as.numeric(mcdu.res[i,"overlap"]),2)))
  
  par(mar=c(5,4,2,0))
  ecospat.plot.overlap.test(sim.ab,"D","Similarity_ab")
  ecospat.plot.overlap.test(sim.a,"D","Similarity_a")
  ecospat.plot.overlap.test(sim.b,"D","Similarity_b")
  
  dev.off()
}

write.table(mcdu.res,paste0(outDir,"Species_niches_simtests.txt"),quote=F,sep="\\t",row.names=F)






#######################################################################
##### ----------------------------------------------------------- #####
#####                     1b-  NICHE MODELING                     #####
#####                  polyploid niche dynamics                   #####
##### ----------------------------------------------------------- #####
#######################################################################


# Matrix of pairwise comparison
px.res = matrix(nrow=4, ncol=9
                , dimnames=list(tetra,c("p1","p2","overlap","sim_ab","sim_di","sim_px","exp","stable","unfill")))
px.res[,"p1"] = c("Um","Um","Ta","Ta")
px.res[,"p2"] = c("Co","Ca","Co","Ca")

##calculate and store Px niche dynamics  in distinct R objects
for(i in tetra){
  name.di = paste(unlist(px.res[i,c("p1","p2")]),collapse="|")
  name.px = i
  row.di = grep(name.di,aegc$sp)
  row.px = grep(name.px,aegc$sp)
  scores.di = pca.cal$li[row.di,]
  scores.px = pca.cal$li[row.px,]
  
  #parents x polyploid niches
  zpx = get(paste0("z.",name.px))
  zdi = ecospat.grid.clim.dyn(scores.clim,scores.clim,scores.di,R)
  assign(paste0("z.p",i),zdi)
}


i="Ge"
for(i in tetra){
  name.di = paste(unlist(px.res[i,c("p1","p2")]),collapse="|")
  name.px = i
  row.di = grep(name.di,aegc$sp)
  row.px = grep(name.px,aegc$sp) ; row.px2 = grep(i,aegc$ID)
  scores.di = pca.cal$li[row.di,]
  scores.px = pca.cal$li[row.px,]
  
  #niches
  zpx = get(paste0("z.",name.px)) # observed polyploid niche
  zdi = get(paste0("z.p",name.px)) # expected polyploid niche (modeled from combination of diploid progenitors occurrences)
  px.res[i,"overlap"] = unlist(ecospat.niche.overlap(zdi,zpx,cor=F)[1])
  
  #niche similarity tests
  sim.di = ecospat.niche.similarity.test(zpx,zdi,rep=100,alternative="greater",rand.type = 2)
  px.res[i,"sim_di"] = sim.di$p.D
  sim.px = ecospat.niche.similarity.test(zdi,zpx,rep=100,alternative="greater",rand.type = 2)
  px.res[i,"sim_px"] = sim.px$p.D
  sim.ab = ecospat.niche.similarity.test(zdi,zpx,rep=100,alternative="greater",rand.type = 1)
  px.res[i,"sim_ab"] = sim.ab$p.D
  
  ##polyploid niches dynamics
  source("ecospat.niche.dyn.index_NEWCALC.R") # ! Used a modified version of ecospat.niche.dyn.index() function ! #
  dyn = ecospat.niche.dyn.index(zdi,zpx)
  
  px.res[i,"exp"] = dyn$dynamic.index.w[1]
  px.res[i,"stable"] = dyn$dynamic.index.w[2]
  px.res[i,"unfill"] = dyn$dynamic.index.w[3]
  
  
  # plots
  pdf(file=paste0(outDir, "p", i," x ",name.px,".pdf"))
  layout(matrix(c(1,1,2,2,
                  1,1,2,2,
                  3,3,4,5,
                  3,3,6,7)
                ,nrow=4,ncol=4,byrow=T))
  par(mar=c(4,4,1,0.5))
  ecospat.plot.niche(zdi, title=name.di, name.axis1="PC1", name.axis2="PC2")
  ecospat.plot.niche(zpx, title=name.px, name.axis1="PC1", name.axis2="PC2")
  ecospat.plot.niche.dyn(zdi,zpx,quant=0.5,title=paste0(name.di," x ",name.px),name.axis1="PC1",name.axis2="PC2")
  par(mar=c(0.5,0.5,0.5,0.5))
  plot.new()
  text(0.5,0.5,paste0("niche overlap\\nD = ",round(as.numeric(px.res[i,"overlap"]),3)
                      ,"\\nniche dynamic\\ne = ",round(as.numeric(px.res[i,"exp"]),3)
                      ,"\\ns = ",round(as.numeric(px.res[i,"stable"]),3)
                      ,"\\nu = ",round(as.numeric(px.res[i,"unfill"]),3)))
  par(mar=c(5,4,2,0))
  ecospat.plot.overlap.test(sim.ab,"D","Similarity_ab")
  ecospat.plot.overlap.test(sim.di,"D","Similarity_di")
  ecospat.plot.overlap.test(sim.px,"D","Similarity_px")
  dev.off()
  
}

write.table(px.res,paste0(outDir,"PxDyn_niche_simtests.txt"),quote=F,sep="\\t")

save.image(paste0(workDir,"EcoGen-1-Niche_modeling.RData"))


