library(raster)
library(maptools)
library(sp)
library(mgcv)
library(ggplot2)
library(PresenceAbsence)
library(rgdal)
library(maps) # Easy access to basic map layers
library(maptools) # Vector data management (sp) 
library(biomod2) # Ensemble SDM package
library(Hmisc)
library(ncdf4)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(rgeos)
memory.limit(10000000)
setwd("G:/SSDM/")
#####Load species distribution data
Takys<-read.csv("Takys_locations.csv",header = T)

#####Load current & future climate data 
path <- file.path("Bio19_current")
files <- list.files(path, pattern='.tif$', full.names=TRUE )
pred_cur<- stack(files)
path <- file.path("Bio19_2050")
files <- list.files(path, pattern='.tif$', full.names=TRUE )
pred_50<- stack(files)
path <- file.path("Bio19_2070")
files <- list.files(path, pattern='.tif$', full.names=TRUE )
pred_70<- stack(files)
#####Load current & future physiological data  
setwd("Combined results")
###EAHT data
EAHT_cur<-raster("EAHT_cur.tif")
EAHT_50<-raster("EAHT_50.tif")
EAHT_70<-raster("EAHT_70.tif")
###Heat stress data
HSTR_cur<-raster("Stress_cur.tif")
HSTR_50<-raster("Stress_50.tif")
HSTR_70<-raster("Stress_70.tif")

#####Combine climatic & physiological data into one raster stack
pred_cur_comb<-stack(pred_cur,EAHT_cur,HSTR_cur)
pred_50_comb<-stack(pred_50,EAHT_50,HSTR_50)
pred_70_comb<-stack(pred_70,EAHT_70,HSTR_70)

names(pred_50)<-names(pred_70)<-names(pred_cur)
names(pred_50_comb)<-names(pred_70_comb)<-names(pred_cur_comb)

## Prepare species distribution data for SDM
coords<-data.frame(lon=Takys$lon ,lat=Takys$lat)
dist_takys <- SpatialPoints(coords, proj4string = pred_cur@crs)

###### Model options set for maxent
myBiomodOptions <- BIOMOD_ModelingOptions(
  MAXENT.Phillips=list(
    path_to_maxent.jar='C:/Users/Administrator/Documents/R/win-library/4.0/dismo/java/maxent.jar'))
#Download maxent.jar from http://www.cs.princeton.edu/~schapire/maxent/; and store it in the working directory
#Also need to download and install java (x64)

#### Model data preparation (similar procedure for 2050, 2070)
mod.dat_clim<-BIOMOD_FormatingData(resp.var=points,
                              expl.var=pred_cur,##### Climatic model
                              PA.nb.rep=1,
                              resp.name = "Takydromus septentrionalis",
                              PA.nb.absences=1000,
                              PA.strategy="random")
mod.dat_clim_phys<-BIOMOD_FormatingData(resp.var=points,
                              expl.var=pred_cur_comb,#####Climatic & physiological model
                              PA.nb.rep=1,
                              resp.name = "Takydromus septentrionalis",
                              PA.nb.absences=1000,
                              PA.strategy="random")

#### Modelling (similar procedure for 2050, 2070)
system.time(
mod_clim<-BIOMOD_Modeling(data=mod.dat_clim, 
                     models=c('GLM','GBM','CTA','ANN','SRE','FDA','RF',"MAXENT.Phillips"),
                     models.options = myBiomodOptions,
                     NbRunEval=3,
                     DataSplit=80, # common practice to validate! 
                     models.eval.meth=c('KAPPA','ROC','TSS'),
                     do.full.models=F,
                     VarImport=3)
)

system.time(
  mod_clim_phys<-BIOMOD_Modeling(data=mod.dat_clim_phys, 
                       models=c('GLM','GBM','CTA','ANN','SRE','FDA','RF',"MAXENT.Phillips"),
                       models.options = myBiomodOptions,
                       NbRunEval=3,
                       DataSplit=80, # common practice to validate! 
                       models.eval.meth=c('KAPPA','ROC','TSS'),
                       do.full.models=F,
                       VarImport=3)
)

#####Get model performance reasults (similar procedure for 2050, 2070)
ModelEval_clim <- get_evaluations(mod_clim)
mod.ROC<- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
mod.KAPPA<- as.data.frame(myBiomodModelEval["KAPPA","Testing.data",,,])
mod.TSS<- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])
var.imp<-get_variables_importance(mod_clim)

ModelEval_clim_phys <- get_evaluations(mod_clim_phys)
mod.ROC2<- as.data.frame(myBiomodModelEval2["ROC","Testing.data",,,])
mod.KAPPA2<- as.data.frame(myBiomodModelEval2["KAPPA","Testing.data",,,])
mod.TSS2<- as.data.frame(myBiomodModelEval2["TSS","Testing.data",,,])
var.imp2<-get_variables_importance(mod_clim_phys)

####Model output (similar procedure for 2050, 2070)
mod.EM.cur_clim<-BIOMOD_EnsembleModeling(modeling.output=mod_clim,
                                        em.by='all',
                                        eval.metric='TSS', 
                                        eval.metric.quality.threshold=0.7, #All models having a score lower than this quality thresholds will not be kept for building ensemble-models
                                        prob.mean=FALSE, 
                                        prob.mean.weight=TRUE) 

mod.EM.cur_clim_phys<-BIOMOD_EnsembleModeling(modeling.output=mod_clim_phys,
                                        em.by='all',
                                        eval.metric='TSS', 
                                        eval.metric.quality.threshold=0.7, #All models having a score lower than this quality thresholds will not be kept for building ensemble-models
                                        prob.mean=FALSE, 
                                        prob.mean.weight=TRUE) 

### SDM projection for current predictors((similar procedure for 2050, 2070)) 
mod.proj.cur_clim<-BIOMOD_Projection(modeling.output=mod_clim,  ####Climatic model 
                                    new.env=pred_cur,          ###Current (2050,2070) climatic data
                                    proj.name='Current',
                                    binary.meth='TSS',
                                    compress='xz'
)

mod.proj.cur_clim_phys<-BIOMOD_Projection(modeling.output=mod_clim_phys, ####Climatic & physiological model 
                                     new.env=pred_cur_comb,         #### Current (2050,2070) climatic and physiological data
                                     proj.name='Current',
                                     binary.meth='TSS',
                                     compress='xz'
)

###  Projection current models (similar procedure for 2050, 2070)

modEF.cur_clim<-BIOMOD_EnsembleForecasting(
  EM.output= mod.EM.cur_clim,
  projection.output=mod.proj.cur_clim,
  proj.name='Current_clim',
  binary.meth='TSS',
  output.format='.grd')

modEF.cur_clim_phys<-BIOMOD_EnsembleForecasting(
  EM.output= mod.EM.cur_clim_phys,
  projection.output=mod.proj.cur_clim_phys,
  proj.name='Current_cur_clim_phys',
  binary.meth='TSS',
  output.format='.grd')

###Mapping different models of current (similar procedure for 2050, 2070)
library(RColorBrewer)
plot(modEF.cur_clim,col = brewer.pal(10,"Spectral")[c(10:3)])
plot(modEF.cur_clim_phys,col = brewer.pal(10,"Spectral")[c(10:3)])

