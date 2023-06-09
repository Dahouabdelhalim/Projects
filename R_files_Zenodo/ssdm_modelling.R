### check for required packages and install if missing
# CRAN
packages <- c("raster","fasterize","sf","parallel","foreach","rJava","rstudioapi","devtools")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
# Github
devtools::install_github("sylvainschmitt/SSDM")
devtools::install_github("valentinitnelav/geobuffer")

### DISCLAIMER: Please make sure to place Maxent software in the appropriate folder of your dismo R package installation and test it before running the following code (or it will fail)


# MODEL TRAINING ----------------------------------------------------------

library(raster)
library(parallel)
library(foreach)
library(SSDM)
library(rJava)

## FILE SETTINGS
# setwd (if not set, in some cases connection is lost)
dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save_path <- "./output/models" # path to save folder
if(!dir.exists(save_path)){
  dir.create(save_path,recursive=TRUE)
}

## DATA SETTINGS
env_files <- list.files("./predictors_training",pattern=".tif",full.names = T)
Env <- stack(env_files)
species_ca_pres <- read.csv("species_ca_pres_final.csv")
pft <- c("") # enter data subsetting keywords (matches rows, which contain all keywords). | means OR, several vector elements mean AND. empty string c("") means ALL.
Xcol <- "longitude" # column name of longitude
Ycol <- "latitude" # column name of latitude
Spcol <- "species" # column name of species
Pcol <- NULL # column name of information on presence/absence
manual_absences <- FALSE # (function disabled for now) if TRUE, absences are manually sampled once before starting the model. If FALSE, new absence data is generated for each model iteration while creating the model

# for disk sampling strategy of pseudo-absences supply mask
abs_mask <- raster("./ca_plus_mask.tif")
buf_width <- 3000 # required! width of buffer around presence points
nocc <- NULL # number of absences to sample. if NULL, default settings will be used for each algorithm

## MODEL SETTINGS
individual.species <- TRUE # if species should be modelled individually and not stacked yet
individual.algos <- TRUE # if algorithms should be ensembled individually first
modelcollect <- FALSE # option for collecting only models that meet threshold criteria. will repeat modelling until the number of model replicates is reached (slow, one core only)
algos <- c('GAM','MARS','MAXENT','RF') # algorithms to use
it <- 10 # no. of model replicates (if manual absence data is supplied, each model will be identical (no stochastic process))
cv_method <- "holdout" # method for cross-validation (holdout, k-fold or LOO (leave-one-out))
testp <- 0.7 # test data fraction for cross-validation
cvrep <- 2 # number of replications for cross-validation ()
thresh_met <- c("AUC","calibration") # threshold metric(s) for model selection
thresh_val <- c(0.7,0.7) # threshold value(s) for metric(s)
stack_method <- "pSSDM" #  method used for SDM stacking
ncores <- detectCores()-2 # number of cores to use for parallel computing (parallelizes along replicates, so length(it) is a good choice). It's highly recommended to set tmp=TRUE when large models are expected, or parallelization will fail.
SDM.projections <- FALSE # whether projections of Algorithm.SDMs should be kept
tmp <- TRUE # logical or character. whether or where temporary rasters should be written to release memory

### Filter selected species from data set
pft_data <- subset(species_ca_pres,apply(species_ca_pres, 1, function(x) all(sapply(pft, function(y) any(grepl(y,x))))))
colnames(pft_data)[which(colnames(pft_data)==Spcol)] <- "species"
pft_data$species <- as.character(pft_data$species)
species_names <- unique(pft_data$species)
species_names


### MODELLING
for(i in 1:length(species_names)){
  species_Occ <- subset(pft_data, species==species_names[i])
  
  if(individual.algos){
    for(j in algos){
      if(is.null(nocc)){
        if(j=="MAXENT"){
          PA <- NULL
        } else {
          nb <- nrow(species_Occ)*2
          if(nb < 1000){nb <- 1000}
          PA <- list(nb=nb,strat='geobuffer',dist=buf_width)
        }
      } else {
        PA <- list(nb=nocc, strat='geobuffer',dist=buf_width)
      } # nocc
      
      if(modelcollect){
        cat(paste("modelling",j, "for",species_names[i],"\\n"))
        ESDM_list <- list()
        a <- 1
        pb <- txtProgressBar(min=0, max=it,initial=0,style = 3)
        while(a<=it){
          dummy <- capture.output(species_SDM <- modelling(j,species_Occ,Env,Xcol=Xcol, Ycol=Ycol,Pcol=Pcol,cv=cv_method, cv.param=c(testp, cvrep), PA=PA, select=TRUE, select.metric = thresh_met, select.thresh = thresh_val, verbose=FALSE))
          if(!is.null(species_SDM)){
            ESDM_list[[a]] <- species_SDM
            #info <- sprintf("%d%% done", round((a/it)*100))
            setTxtProgressBar(pb, a)
            #cat(paste(a, "of", it,"\\n"))
            a <- a+1
          }
        }
        ESDM_list <-  c(ESDM_list,uncertainty=FALSE, cores=0, SDM.projections = SDM.projections,ensemble.thresh=c(0))
        species_ESDM <- do.call(SSDM::ensemble,ESDM_list)
      }else{
        species_ESDM <- ensemble_modelling(j, species_Occ, Env, Xcol=Xcol, Ycol=Ycol, Spcol=Spcol, Pcol=Pcol, rep=it, cv=cv_method, cv.param=c(testp, cvrep), ensemble.metric=c(thresh_met), ensemble.thresh=c(thresh_val), verbose=TRUE, name=species_names[i], uncertainty=FALSE, cores=ncores, SDM.projections=SDM.projections, PA=PA,tmp=tmp)
      }
      saveRDS(species_ESDM,file=paste0(save_path,"/",gsub(" ",".",species_names[i]),"_",j,".rds"))
      gc()
    }
  } else{
    species_ESDM <- ensemble_modelling(algos, species_Occ, Env, Xcol=Xcol, Ycol=Ycol, Spcol=Spcol, Pcol=Pcol, rep=it, cv=cv_method, cv.param=c(testp, cvrep), ensemble.metric=c(thresh_met), ensemble.thresh=c(thresh_val), verbose=TRUE, name=species_names[i], uncertainty=FALSE, cores=ncores, SDM.projections=SDM.projections, PA=PA,tmp=tmp)
    saveRDS(species_ESDM,file=paste0(save_path,"/",gsub(" ",".",species_names[i]),".rds"))
    gc()
  }
}


# MODEL PROJECTION --------------------------------------------------------

### DATA SETTINGS ###
# non-climate (static) predictors
env_path <- "./predictors_projection/env"
# climate scenarios
clim_path <- "./predictors_projection/clim_scen"
# pattern how to grep each climate scenario (vector to loop through).  Will also be appended as savename suffix (necessary!)
gcm <- c("cc","ha","mp")
rcp <- c("26","45","85")
historic <- "present" # set to NULL, if not desired
clim_pattern <- c(historic,sort(as.vector(outer(gcm,rcp,FUN=paste,sep="_"))))

# supply names of climate variables used for training (for scenario variable renaming)
clim_names <- c("annprec","maxtemp","mintemp","precseas")

# specify sub-extent if only a geographic subset should be projected
sub_ext <- NULL

# supply path to model(s) (all .rds files within given directory will be tried)
mod_path <- "./output/models"
path_recursive <- FALSE # whether to look in sub-folders also
file_pattern <- ".rds" # pattern to find files, empty character means all. specific matching can be done via glob2rx

# supply path to folder where FINAL projections should be saved
save_path <- "./output/projections"
if(!dir.exists(save_path)){
  dir.create(save_path,recursive = TRUE)
}

# logical or path to folder where TEMPORARY projections should be saved
tmp_path <- TRUE

# chunk-wise 
minimal.memory <- TRUE

# number of cores to use in parallel
ncores <- detectCores()-2

# PROJECTION SETTINGS
SDM.projections <- TRUE # whether single algorithm projections should be returned (only relevant if use_SSDM=TRUE)
output.format <- 'rasters' # 'model' or 'rasters'
uncertainty <- FALSE # whether to calculate and return an uncertainty map
binary <- FALSE # whether to return the binary map
weight <- FALSE # whether projections should be weighted (right now only works for AUC)

# data min-max normalization parameter -- VALUES NOT CORRECT (if normalized data is required, i.e. for ANN)
# ann_minmax <- list(annprec=c(157.0003,12564.8770),maxtemp=c(-0.9382346,38.7901421),mintemp=c(-11.45989,29.10167),precseas=c(6,152))



### START PROJECTION
mod_files <- list.files(mod_path,pattern=file_pattern, recursive= path_recursive, full.names = TRUE)
env_files <- list.files(env_path,pattern=".tif",full.names=TRUE)
# number of scenarios being projected (usually clim_pattern)
nscen <- length(clim_pattern)

for(k in 1:length(mod_files)){
  modname <- sub(".rds","",basename(mod_files[k]))
  species_ESDM <- readRDS(mod_files[k])
  outdir <- modname
  if(!dir.exists(paste0(save_path,"/",outdir))){
    dir.create(paste0(save_path,"/",outdir))
  }
  #
  for(j in 1:nscen){
    print(paste("Scenario",clim_pattern[j],"for",modname,"   ",Sys.time()))
    if(!is.null(clim_path)){
      clim_files <- list.files(clim_path, pattern=clim_pattern[j], full.names = TRUE)
      Env <- stack(c(clim_files,env_files))
      names(Env)[1:length(clim_files)] <- clim_names
    } else{Env <- stack(env_files)} 
    # write console output to dummy for muting status messages
    dummy <- capture.output(proj <- SSDM::project(obj=species_ESDM, Env=Env,SDM.projections=SDM.projections, output.format=output.format, uncertainty=uncertainty,cores=ncores,minimal.memory=minimal.memory,tmp=tmp_path))
    if(output.format=='model'){
      saveRDS(proj,file = paste0(save_path,"/",modname,"_", clim_pattern[j]))
    } else{
      writeRaster(proj$projection,filename=paste0(save_path,"/",outdir,"/",modname,"_",clim_pattern[j]),format="GTiff")
      if(binary){
        writeRaster(proj$binary,filename=paste0(save_path,"/",outdir,"/",modname,"_",clim_pattern[j],"_bin"),format="GTiff")
      }
      if(uncertainty){
        writeRaster(proj$uncertainty,filename=paste0(save_path,"/",outdir,"/",modname,"_",clim_pattern[j],"_uncertainty"),format="GTiff")
      }
      if(SDM.projections){
        dir.create(paste0(save_path,"/",outdir,"/",modname,"_",clim_pattern[j]))
        sdm_stack <- stack(lapply(proj$sdms,function(x) x$projection))
        writeRaster(sdm_stack,filename = paste0(save_path,"/",outdir,"/",modname,"_",clim_pattern[j],"/",modname,"_",clim_pattern[j],"-",c(1:length(proj$sdms))),bylayer=TRUE,format="GTiff")
      }
    }
    rm(proj,sdm_stack)
    gc()
    unlink(paste0(tmp_path,"/.models"),recursive = TRUE, force = TRUE)
  } #j
  
} #k


# MODEL STACKING ----------------------------------------------------------

### FILE SETTINGS
# path to where models are stored
modpath <- "./output/models"

# path to projections (these are re-inserted scenario-wise)
projpath <- "./output/projections"

# path to where SSDMs should be saved, if writemods or writeproj is TRUE. also useful for reducing working memory use
writemods <- TRUE
savemodpath <- "./output/stacked_models"
if(!dir.exists(savemodpath)){
  dir.create(savemodpath,recursive=TRUE)
}
writeproj <- FALSE # whether to write species richness and endemism maps
writeensproj <- TRUE # whether to write ensemble projections
saveprojpath <- "./output/stacked_projections"
if(!dir.exists(saveprojpath)){
  dir.create(saveprojpath,recursive=TRUE)
}

### DATA SELECTION
## PFT grouping
dryacq <- c("Brosimum.alicastrum","Calycophyllum.candidissimum","Enterolobium.cyclocarpum","Pachira.quinata")
drycon <- c("Alvaradoa.amorphoides","Byrsonima.crassifolia","Leucaena.leucocephala","Vachellia.farnesiana")
wetacq <- c("Cecropia.obtusifolia","Ochroma.pyramidale","Schizolobium.parahyba","Vochysia.ferruginea")
wetcon <- c("Calophyllum.brasiliense","Carapa.guianensis","Dialium.guianense","Symphonia.globulifera")
conifers <- c("Pinus.ayacahuite","Pinus.caribaea","Pinus.oocarpa","Pinus.tecunumanii")
montane <- c("Alnus.acuminata","Cornus.disciflora","Drymis.granadensis","Weinmannia.spp")
generalists <- c("Guazuma.ulmifolia","Simarouba.amara","Spondias.mombin")

pft <- list(dryacq=dryacq,drycon=drycon,wetacq=wetacq,wetcon=wetcon,conifers=conifers,montane=montane,generalist=generalists)
stackname <- names(pft)

## algorithms
# algorithms to stack (will keep separated between algorithms unless ensemble.algos is set to TRUE)
algos <- c("GAM","MARS","MAXENT","RF") # please order alphabetically
# scenarios
gcm <- c("cc","ha","mp")
rcp <- c("26","45","85")
historic <- "present" # set to NULL, if not desired
clim_pattern <- c(historic,sort(as.vector(outer(gcm,rcp,FUN=paste,sep="_"))))

### STACKING SETTINGS
# ensemble algorithms for each species before stacking
ensemble.algos <- TRUE
# stacking methods
return.stack <- TRUE # whether to return an SSDM object. Will always be for the first element of clim_pattern
stackmethod <- "pSSDM" # method used for stacking ESDMs
endemism <- c("WEI","Binary") # methods used for calculating endemism maps
uncertainty <- FALSE

### STACKING
speciesdir <- list.files(modpath,full.names = TRUE,recursive = TRUE)
pftmaps <- list()
eval_df <- NULL
for(m in 1:length(pft)){
  speciesinds <- as.vector(sapply(pft[[m]],function(x){grep(x,speciesdir)}))
  speciesfiles <- speciesdir[speciesinds]
  
  speciesalgomaps <- list()
  if(!ensemble.algos){
    for(k in 1:length(algos)){
      algofiles <- speciesfiles[grep(algos[k],speciesfiles)]
      # load and prepare ESDMs
      speciesmods <- lapply(algofiles,function(x){
        algomod <- readRDS(x)
        # some of the models have empty algorithm.evaluation slots (not yet clear why). rebuild them from the evaluation slot
        if(length(algomod@algorithm.evaluation)<1){
          algomod@algorithm.evaluation <- cbind(algomod@evaluation,kept.model=length(algomod@sdms))
          rownames(algomod@algorithm.evaluation) <- algomod@name
        }
        # rename models into default format
        algomod@name <- sub("[.]"," ",algomod@name)
        algomod@name <- sub("_",".",algomod@name)
        return(algomod)
      })
      
      #  re-insert projections for each scenario, then stack
      ssdmmaps <- list()
      for(j in 1:length(clim_pattern)){
        print(paste("Stacking scenario",clim_pattern[j],"for",stackname[m],algos[k]))
        for(i in 1:length(algofiles)){
          speciesmods[[i]]@projection <- raster(paste0(projpath,"/",sub(".rds","",basename(algofiles[[i]])),"/",sub(".rds","",basename(algofiles[[i]])),"_",clim_pattern[j],".tif"))
          # rebuild binary raster for stack evaluation
          speciesmods[[i]]@binary <- reclassify(speciesmods[[i]]@projection, c(-Inf,speciesmods[[i]]@evaluation$threshold,0,speciesmods[[i]]@evaluation$threshold,Inf,1))
        } # i algofiles
        speciesstack <- do.call(SSDM::stacking, c(speciesmods,method=stackmethod,uncertainty=uncertainty,name=paste0(stackname[m],"_",algos[k],"_",clim_pattern[j]),verbose=FALSE))
        if(j==1){
          species_SSDM <- speciesstack
        }
        ssdmmaps[[j]] <- list(diversity=speciesstack@diversity.map,endemism=speciesstack@endemism.map)
        names(ssdmmaps)[j] <- clim_pattern[j]
        
      } # j clim_pattern
      
      if(writemods){
        saveRDS(speciesstack,file = paste0(savemodpath,"/",stackname[m],"_",algos[k],".rds"))
      }
      if(writeproj){
        for(r in 1:length(ssdmmaps)){
          ssdmmaps[[r]]$diversity <- writeRaster(ssdmmaps[[r]]$diversity,filename = paste0(saveprojpath,"/",stackname[m],"_",algos[k],"_",clim_pattern[r],"_div.tif"))
          ssdmmaps[[r]]$endemism <- writeRaster(ssdmmaps[[r]]$endemism,filename = paste0(saveprojpath,"/",stackname[m],"_",algos[k],"_",clim_pattern[r],"_end.tif"))
        }
      }
      speciesalgomaps[[k]] <- ssdmmaps
      names(speciesalgomaps)[k] <- algos[k]
      
    } # k algos
    pftmaps[[m]] <- speciesalgomaps
    names(pftmaps)[m] <- stackname[m]
  } else {
    # first ensemble each species by algo, then stack (original SSDM procedure)
    # do this for each scenario
    for(scen in 1:length(clim_pattern)){
      print(paste("Stacking PFT",names(pft)[m] ,"for...",clim_pattern[scen]))
      species_ESDM_ls <- list()
      species <- pft[[m]]
      # ensemble individual species across algorithms
      for(k in 1:length(species)){
        print(paste("...ensembling",species[k]))
        ispeciesfiles <- speciesfiles[grep(species[k],speciesfiles)]
        # load individual models, rename and list
        speciesalgo <- list()
        for(j in 1:length(ispeciesfiles)){
          speciesalgoname <- sub(".rds","",basename(ispeciesfiles[j]))
          ispeciesalgo <- readRDS(ispeciesfiles[j])
          ispeciesalgo <- ispeciesalgo@sdms
          # re-insert projections and restore binary maps
          for(i in 1:length(ispeciesalgo)){
            ispeciesalgo[[i]]@projection <- raster(paste0(projpath,"/",speciesalgoname,"/",speciesalgoname,"_",clim_pattern[scen],"/",speciesalgoname,"_",clim_pattern[scen],"-",i,".tif"))
            ispeciesalgo[[i]]@binary <- reclassify(ispeciesalgo[[i]]@projection, c(-Inf,ispeciesalgo[[i]]@evaluation$threshold,0,ispeciesalgo[[i]]@evaluation$threshold,Inf,1))
          } # i individual algorithm models
          
          speciesalgo <- c(speciesalgo,ispeciesalgo)
        } # j individual species
        speciesalgo <- c(speciesalgo, weight=FALSE, SDM.projections=FALSE, uncertainty=FALSE, ensemble.thresh=0, name=species[k], verbose=FALSE)
        species_ESDM <- do.call(ensemble,speciesalgo)
        # write evaluation statistics to data.frame
        if(scen==1){
          species_eval <- species_ESDM@algorithm.evaluation
          species_eval <- data.frame(species=species[k],algorithm=algos,species_eval)
          eval_df <- rbind(eval_df, species_eval)
        }
        # save projection, if desired
        if(writeensproj){
          if(!dir.exists(paste0(saveprojpath,"/ensembles"))){
            dir.create(paste0(saveprojpath,"/ensembles"))
          }
          writeRaster(species_ESDM@projection,filename=paste0(saveprojpath,"/ensembles/",species[k],"_",clim_pattern[scen],".tif"))
        }
        species_ESDM_ls <- c(species_ESDM_ls,species_ESDM)
      } # k species
      # stack species
      species_ESDM_ls <- c(species_ESDM_ls,name=names(pft)[m],method=stackmethod, eval=TRUE, uncertainty=FALSE, verbose=FALSE)
      species_ESDM_ls[["endemism"]] <- endemism
      if(stackmethod%in%c("PRR.MEM","PRR.pSSDM")){
        species_ESDM_ls[["Env"]] <- Env
      }
      species_SSDM <- do.call(stacking,species_ESDM_ls)
      writeRaster(species_SSDM@diversity.map,filename=paste0(saveprojpath,"/",names(pft)[m],"_div_",clim_pattern[scen],".tif"))
      writeRaster(species_SSDM@endemism.map,filename=paste0(saveprojpath,"/",names(pft)[m],"_end_",clim_pattern[scen],".tif"))
      if(writemods & clim_pattern[scen]=="present"){
        saveRDS(species_SSDM, file=paste0(savemodpath,"/",names(pft)[m],".rds"))
      }
    } # end scen clim_pattern
    
  } # end if ensemble.algo
  
} # m pft 
