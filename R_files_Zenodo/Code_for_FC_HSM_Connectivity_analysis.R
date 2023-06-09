####Code to run habitat suitability analysis using maxent models through the dismo package. Note that dismo.jar needs to be downloaded separately
library(dismo)
library(raster)
library(blockCV)
library(raster)
library(sf)

setwd("./FC_data")#set the working directory to where the data is
occ=read.csv("./Occ_thin_Ind.csv")
PBg=read.csv("./PB1_Ind.csv")

rasterlist= list.files("Rasters", pattern=".tif", full.names=T) 
predictors=stack(rasterlist)

#Check for spatial autocorrelation in variables to estimate block size for spatial cross validation
sac=spatialAutoRange(rasterLayer=predictors,sampleNumber=5000,doParallel=TRUE,showPlots=TRUE)
plot(sac)

pb_data=st_as_sf(read.csv("PB1_Ind.csv"), coords=c("x", "y"), crs=crs(predictors))

# spatial blocking by specified range with random assignment
sb=spatialBlock(speciesData=pb_data,
                   species="Species",
                   rasterLayer=predictors,
                   theRange=1000000, # size of the blocks, estimated through 
                   k=5, # number of spatial replicates
                   selection="random",
                   iteration=100, # find evenly dispersed folds
                   biomod2Format=TRUE,
                   xOffset=0, # shift the blocks horizontally
                   yOffset=0)

folds=sb$foldID

final_results=data.frame()
for(i in 1:5)#5 spatial replicates were prepared
{
  occtest=subset(PBg[folds==i, ], Species==1)[,2:3]
  occtrain =subset(PBg[folds!=i, ], Species==1)[,2:3]
  bgtest=subset(PBg[folds==i, ], Species==0)[,2:3]
  bgtrain=subset(PBg[folds!=i, ], Species==0)[,2:3]
  #model name=mod.single.Ind.backgroundpoints.betamul.modelnumber.replicate
  
  ## Note that the following code (within the loop) must be re-run in different configurations (e.g., with varying beta multiplier values) as needed for the same subset
  #model1
  assign(paste0("mod.single.Ind.1.1.1.",i), maxent(predictors, occtrain, a=bgtrain, args=c("betamultiplier=1","linear=true","quadratic=false","product=false","threshold=false","hinge=false","maximumiterations=5000", "responsecurves=true"), factors="LULC_2012"))
  assign("temp_mod", get(paste0("mod.single.Ind.1.1.1.",i)))
  mod.pre=predict(temp_mod, predictors)
  assign((paste0("mod.pre.Ind.1.1.1.",i)), get("mod.pre"))
  HSval=extract(mod.pre, occ)
  temp=data.frame(temp_mod@results)
  mtp=temp[which(rownames(temp)=="Minimum.training.presence.training.omission"),1]
  p10=temp[which(rownames(temp)=="X10.percentile.training.presence.training.omission"),1]
  ess=temp[which(rownames(temp)=="Equal.training.sensitivity.and.specificity.training.omission"),1]
  reg_gain=temp[which(rownames(temp)=="Regularized.training.gain"),1]
  unreg_gain=temp[which(rownames(temp)=="Unregularized.training.gain"),1]
  prev=temp[which(rownames(temp)=="Prevalence..average.probability.of.presence.over.background.sites."),1]
  OmR_mtp=(length(which(HSval<mtp)))/length(HSval)
  OmR_p10=(length(which(HSval<p10)))/length(HSval)
  OmR_ess=(length(which(HSval<ess)))/length(HSval)
  etrain=dismo::evaluate(p=occtrain, a=bgtrain, model=temp_mod, x=predictors)
  etest=dismo::evaluate(p=occtest, a=bgtest, model=temp_mod, x=predictors)
  
  #Storing results
  final_results=rbind(final_results, cbind("Model"="L", "Mod"="Single_Ind", "BetaMul"=1, "Rep"=i, n_test=nrow(occtest),  "Reg_gain"=reg_gain, "Unreg_gain"=unreg_gain, "Prevelance"=prev, "AUC_train"=etrain@auc, "AUC_test"=etest@auc, "AUC_diff"=(etrain@auc-etest@auc),  "kappa_train"=mean(etrain@kappa), "kappa_test"=mean(etest@kappa), "TSS_train"=mean(etrain@TPR+etrain@TNR-1), "TSS_test"=mean(etest@TPR+etest@TNR-1), "MTP"=mtp, "p10"=p10, "ess"=ess, "OmR_mtp"=OmR_mtp, "OmR_p10"=OmR_p10, "OmR_ess"=OmR_ess))
  
}

  final_results[,3:21]=sapply(final_results[,3:21],as.numeric)
  
  library(dplyr)
  final_results_mean=final_results %>% group_by(Model, Mod, BetaMul) %>% summarise_each(funs(mean))
  final_results_mean$AUC_diff=final_results_mean$AUC_train-final_results_mean$AUC_test
  final_results_sd=final_results %>% group_by(Model, Mod, BetaMul) %>% summarise_each(funs(sd))
  
  write.csv(final_results, "Final_results.csv", row.names=F)
  write.csv(final_results_mean, "Final_results_mean.csv", row.names=F)
  write.csv(final_results_sd, "Final_results_sd.csv", row.names=F)
  
  ##################
  #Averaging models
  #Raster
  #Note that 'g' and 'h' corresponded to beta multiplier values and model number in our original code
  final_raster=calc(stack(get(paste0("mod.pre.Ind.1.",g,".",h,".1")),get(paste0("mod.pre.Ind.1.",g,".",h,".2")),get(paste0("mod.pre.Ind.1.",g,".",h,".3")),get(paste0("mod.pre.Ind.1.",g,".",h,".4")),get(paste0("mod.pre.Ind.1.",g,".",h,".5"))), mean)
  plot(final_raster)
  writeRaster(final_raster, "Maxent.tif")
  
  #Permutation contribution
  perm.cont=rbind(get(paste0("mod.single.Ind.1.",g,".",h,".1"))@results[15:22],get(paste0("mod.single.Ind.1.",g,".",h,".2"))@results[15:22],get(paste0("mod.single.Ind.1.",g,".",h,".3"))@results[15:22],get(paste0("mod.single.Ind.1.",g,".",h,".4"))@results[15:22],get(paste0("mod.single.Ind.1.",g,".",h,".5"))@results[15:22])
  colnames(perm.cont)=c("Bio17","Bio18","Bio19","Bio02","Elevation","LULC_2012","Pop_den_2020","SWC_avg")
  perm.cont.final=cbind(apply(perm.cont,2,mean),apply(perm.cont,2,sd))
  colnames(perm.cont.final)=c("Mean", "SD")
  write.csv(perm.cont.final, "perm.cont.final.csv", row.names=T)
  
  #Variable importance
  var.imp=rbind(get(paste0("mod.single.Ind.1.",g,".",h,".1"))@results[7:14],get(paste0("mod.single.Ind.1.",g,".",h,".2"))@results[7:14],get(paste0("mod.single.Ind.1.",g,".",h,".3"))@results[7:14],get(paste0("mod.single.Ind.1.",g,".",h,".4"))@results[7:14],get(paste0("mod.single.Ind.1.",g,".",h,".5"))@results[7:14])
  colnames(var.imp)=c("Bio17","Bio18","Bio19","Bio02","Elevation","LULC_2012","Pop_den_2020","SWC_avg")
  var.imp.final=cbind(apply(var.imp,2,mean),apply(var.imp,2,sd))
  colnames(var.imp.final)=c("Mean", "SD")
  write.csv(var.imp.final, "var.imp.final.csv", row.names=T)
  
  
  #response curve for best replicate
  response(get(paste0("mod.single.Ind.1.",g,".",h,".5")))
  #response curve for all
  HS=extract(final_raster, occ)
  var=extract(predictors, occ)
  resp=data.frame(cbind(var,HS))
  library(ggplot2)
  library(ggpubr)
  a=ggplot(resp, aes(x=Bio17, y=HS))+geom_smooth(se=F)
  b=ggplot(resp, aes(x=Bio18, y=HS))+geom_smooth(se=F)
  c=ggplot(resp, aes(x=Bio19, y=HS))+geom_smooth(se=F)
  d=ggplot(resp, aes(x=Bio2, y=HS))+geom_smooth(se=F)
  e=ggplot(resp, aes(x=Elevation, y=HS))+geom_smooth(se=F)
  f=ggplot(resp, aes(x=as.factor(LULC_2012), y=HS))+geom_boxplot()
  g=ggplot(resp, aes(x=Pop_den_2020, y=HS))+geom_smooth(se=F)
  h=ggplot(resp, aes(x=SWC_avg, y=HS))+geom_smooth(se=F)
  ggarrange(a,b,c,d,e,f,g,h, nrow=2, ncol=4)
  
  plot(final_raster)
  
 ####################
  #Calculating P/E ratio and mapping habitat suitability
  library(ecospat)
  classification=data.frame(ecospat.boyce (na.omit(getValues(final_raster)), na.omit(extract(final_raster, occ)), nclass=0, window.w="default", res=100, PEplot=F, rm.duplicate=TRUE, method='spearman'))
  ggplot(classification, aes(x=HS, y=F.ratio))+geom_line(cex=.6)+geom_hline(yintercept=5)+labs(x="Habitat Suitability", y="P/E ratio") +  theme(axis.text.x=element_text(size=14,family="sans"), axis.title.x=element_text(size=16,family="sans"), axis.text.y=element_text(size=14,family="sans"), axis.title.y=element_text(size=16,family="sans"))
  
  ggsave(path=getwd(), filename="P.E.ratio.jpeg", width=300, height=200, units="mm", device='jpeg', dpi=600)
  
  # @0.70,HS, it is 20x more likely to predict a presence than expected
  values(final_raster)[which(values(final_raster)>0.70)]=1
  values(final_raster)[which(values(final_raster)<=0.70)]=0
  plot(final_raster)
  
  #for each raster
  final_raster=get(paste0("mod.pre.Ind.1.",g,".",h,".5"))
  plot(final_raster)
  
  # classification_allRas=data.frame()
  classification_allRas=rbind(classification_allRas, cbind("Model"=5, "HS"=classification$HS, "P.E_ratio"=classification$F.ratio))
  
  library(ggplot2)
  ggplot(classification_allRas, aes(x=HS, y=P.E_ratio))+geom_smooth(method="loess", level=0.95)+geom_point(aes(col=as.factor(Model)), cex=1.05)+geom_hline(yintercept=1)+ylim(c(0,60))
  
########################################################################################
  #Omnidirectional connectivity analysis
  #Note that circuitscpae needs to be downloaded separately
  setwd("./Connectivity")
  #Splitting Raster
  library("SpaDES")
  library(raster)
  #Initializing buffer pixel length
  bufrow=30
  bufcol=30
  buflencol=5
  buflenrow=5
  
  r=raster("./Resistance.tif")
  ran=extend(r,c(bufrow,bufcol),value=NA)
  values(ran)=sample(values(r),ncell(ran),replace=T)
  
  #adding buffer around raster
  #create and set a new working directory
  #create different directories for storing multiple output files
  dir.create("./Tiles_buf")
  dir.create("./SourceWE_buf")
  dir.create("./SourceNS_buf")
  dir.create("./CS_outNS")
  dir.create("./CS_outWE")
  dir.create("./CS_runNS")
  dir.create("./CS_runWE")
  dir.create("./MosaicNS")
  dir.create("./MosaicWE")
  dir.create("./Current")
  
  
  l.rtl=dir("./", pattern=".tif")
  tifdir="./"
  for (i in l.rtl) {
    rtl=raster(i) # read tile
    rtl.temp=extend(rtl, c(bufrow,bufcol), value=NA) # extend by 100% buffer
    ran.buf=crop(ran, rtl.temp) # crop random raster by buffered tile
    rtl.buf=cover(rtl.temp, ran.buf) # replace NA in buffer with random values
    we.buf=rtl.buf
    we.buf[]=NA; we.buf[,1]=1; we.buf[,ncol(rtl.buf)]=2 #Creating thin west/east strips
    ns.buf=rtl.buf
    ns.buf[]=NA; ns.buf[1,]=1; ns.buf[nrow(rtl.buf),]=2 #Creating thin north/south strips
    
    writeRaster(rtl.buf, datatype='FLT4S', filename=paste(tifdir, "Tiles_buf/", i, ".asc", sep
                                                           =""),overwrite=TRUE) # write buffered resistance tile as .asc
    writeRaster(we.buf, datatype='FLT4S', filename=paste(tifdir, "SourceWE_buf/SourceWE_", i,
                                                           ".asc", sep=""),overwrite=TRUE) # write current source WE tile as .asc
    writeRaster(ns.buf, datatype='FLT4S', filename=paste(tifdir, "SourceNS_buf/SourceNS_", i,
                                                           ".asc", sep=""),overwrite=TRUE) # write current source NS tile as .asc
  }
  
  
  #For WE output
  #Preparing tiles for CS input
  l.tile=dir("./Tiles_buf/", pattern=".asc")
  l.source=dir("./SourceWE_buf/", pattern=".asc")
  
  l.tile.ini=list()
  for (i in l.tile) {
    tile.ini=c("[circuitscape options]",
                  "data_type=raster",
                  "scenario=pairwise",
                  "write_cur_maps=True",
                  "write_cum_cur_map_only=True",
                  "habitat_map_is_resistances=True",#resistance vs conductance
                  paste("habitat_file=", getwd(), "/Tiles_buf/", i, sep=""),
                  paste("output_file=", getwd(), "/CS_outWE/", i, ".out", sep=""))
    l.tile.ini[[i]]=append(l.tile.ini[[i]], tile.ini)
  }
  
  l.source.ini=list()
  for (i in l.source) {
    source.ini=c(paste("point_file=", getwd(), "/SourceWE_buf/", i, sep=""))
    l.source.ini[[i]]=append(l.source.ini[[i]], source.ini)
  }
  
  l.run=list()
  for (i in 1:length(l.tile.ini)) {
    l.tile.ini[[i]]=c(l.tile.ini[[i]], l.source.ini[[i]])
    writeLines(l.tile.ini[[i]], paste(getwd(), "/CS_runWE/Tile", i, ".ini", sep=""))
    l.run[i]=paste(getwd(), "/CS_runWE/Tile", i, ".ini", sep="")
  }
  
  #Running the circuitscape program
  CS_exe='C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget quotation marks in "Program Files"
  sapply(l.run, function(x) system(paste(CS_exe, x)))
  
  #Unidirectional current
  #rm(list=ls()) # clear workspace
  tifdir="./MosaicWE/"
  projection="+proj=longlat +datum=WGS84 +no_defs " # define a projection system (should match your input rasters)
  
  l.in=dir(paste0(getwd(),"./CS_outWE/"), pattern="curmap.asc")
  l.rtl=dir("./", pattern=".grd")
  l.cur=list()
  c=1
  for (i in l.in) {
    cur=raster(paste0("./CS_outWE/",i))
    cur.crop=crop(cur, r)#only for single raster
    crs(cur.crop)=projection
    writeRaster(cur.crop, filename=paste(tifdir,"/FC_WE.tif", sep=""))
    l.cur[[i]]=cur.crop
    c=c+1
  }
  cur.mos=raster("./CS_outWE/Resistance.tif.asc_cum_curmap.asc")
  
  
  #For NS
  #Preparing tiles for CS input
  l.tile=dir("./Tiles_buf/", pattern=".asc")
  l.source=dir("./SourceNS_buf/", pattern=".asc")
  
  l.tile.ini=list() # same as: vector("list")
  for (i in l.tile) {
    tile.ini=c("[circuitscape options]",
                  "data_type=raster",
                  "scenario=pairwise",
                  "write_cur_maps=True",
                  "write_cum_cur_map_only=True",
                  "habitat_map_is_resistances=True",#resistance vs conductance
                  paste("habitat_file=", getwd(), "/Tiles_buf/", i, sep=""),
                  paste("output_file=", getwd(), "/CS_outNS/", i, ".out", sep=""))
    l.tile.ini[[i]]=append(l.tile.ini[[i]], tile.ini)
  }
  
  l.source.ini=list()
  for (i in l.source) {
    source.ini=c(paste("point_file=", getwd(), "/SourceNS_buf/", i, sep=""))
    l.source.ini[[i]]=append(l.source.ini[[i]], source.ini)
  }
  
  l.run=list()
  for (i in 1:length(l.tile.ini)) {
    l.tile.ini[[i]]=c(l.tile.ini[[i]], l.source.ini[[i]])
    writeLines(l.tile.ini[[i]], paste(getwd(), "/CS_runNS/Tile", i, ".ini", sep=""))
    l.run[i]=paste(getwd(), "/CS_runNS/Tile", i, ".ini", sep="")
  }
  
  #Running the circuitscape program
  CS_exe='C:/"Program Files"/Circuitscape/cs_run.exe' # Don't forget quotation marks in "Program Files"
  sapply(l.run, function(x) system(paste(CS_exe, x)))
  
  #Unidirectional current
  tifdir="./MosaicNS/"
  projection="+proj=longlat +datum=WGS84 +no_defs " # define a projection system (should match your input rasters)
  
  l.in=dir(paste0(getwd(), "/CS_outNS/"), pattern="curmap.asc")
  l.rtl=dir("./", pattern=".grd")
  l.cur=list()
  c=1
  for (i in l.in) {
    cur=raster(paste0("./CS_outNS/",i))
    cur.crop=crop(cur,r)
    crs(cur.crop)=projection
    writeRaster(cur.crop, filename=paste(tifdir, "FC_NS.tif", sep=""))
    l.cur[[i]]=cur.crop
    cur.crop
    c=c+1
  }
  cur.mos=raster("./CS_outNS/Resistance.tif.asc_cum_curmap.asc")
  
  #Omnidirectional current density
  res=raster("./Resistance.tif")
  tifdir="./Current"
  mos.WE=raster(paste(getwd(), "/MosaicWE/FC_WE.tif", sep=""))
  mos.WE=crop(mos.WE, res)
  mos.NS=raster(paste(getwd(), "/MosaicNS/FC_NS.tif", sep=""))
  mos.NS=crop(mos.NS, res)
  
  #res=crop(res, mos.WE)
  mos.WE=mask(mos.WE, res)
  mos.WE; plot(mos.WE)
  writeRaster(mos.WE, filename=paste(tifdir, "/FC_WE_mask.tif", sep=""))
  mos.NS=mask(mos.NS, res)
  mos.NS; plot(mos.NS)
  writeRaster(mos.NS, filename=paste(tifdir, "/FC_NS_mask.tif", sep=""))
  
  #Combine NE and WS rasters
  mos=overlay(mos.NS, mos.WE, fun=function(x, y) {return(log10(x * y))} )
  writeRaster(mos, filename=paste(tifdir, "/Current_final.tif", sep=""))
  plot(mos)
  hist(mos)
  