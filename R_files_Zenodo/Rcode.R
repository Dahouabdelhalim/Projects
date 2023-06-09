
##### Following codes are provided by Shijia Peng, Please contact me "shijia.peng@pku.edu.cn" if you have any questions ######


### model construction ###
library(biomod2)
library(snowfall)
biomod.spdis8732<-get(load("D:/Shijia/J.APPL.E/Modeldata/biomod_spdis.RData"))
data.XY<-read.csv("D:/Shijia/J.APPL.E/Modeldata/biomodXY.csv")
data.climate<-read.csv("D:/Shijia/J.APPL.E/Modeldata/current_climate.csv",row.names=1)
bc26<-read.csv("D:/Shijia/J.APPL.E/Modeldata/2070/future_rcp2.6.csv",row.names=1)
bc60<-read.csv("D:/Shijia/J.APPL.E/Modeldata/2070/future_rcp6.0.csv",row.names=1)
bc85<-read.csv("D:/Shijia/J.APPL.E/Modeldata/2070/future_rcp8.5.csv",row.names=1)
setwd("D:/Shijia/J.APPL.E/Model")
data.sp<-as.data.frame(biomod.spdis8732)
biomod.spdis <- list()
for(i in 1:1027) biomod.spdis[[i]] <- data.sp[i]
MyBiomodSF <- function (sp.dis){
  errin<-tryCatch({
    myRespName<- names(sp.dis)
    myResp<-as.numeric(sp.dis[,myRespName])
    myRespxy<-data.XY
    myExpl<-data.climate
    myBiomodData<-BIOMOD_FormatingData(resp.var=myResp,
                                       expl.var=myExpl,
                                       resp.xy=myRespxy,
                                       resp.name=myRespName)
    myBiomodOption <- BIOMOD_ModelingOptions()
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                        models = c('CTA','RF','GBM','GLM',"MAXENT.Phillips"),
                                        models.options = myBiomodOption,
                                        NbRunEval=10,
                                        DataSplit=80,
                                        Yweights=NULL,
                                        VarImport=3,
                                        models.eval.meth =c('TSS'),
                                        SaveObj = T,
                                        modeling.id = paste(myRespName,"Firstmodeling",sep=""),
                                        rescal.all.models = FALSE,
                                        do.full.models = F)
    myBiomodProjection <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                            new.env = myExpl,
                                            proj.name = 'current',
                                            selected.models = 'all',
                                            compress = T,
                                            build.clamping.mask = FALSE)
    Proj2070Rcp26<-BIOMOD_Projection(modeling.output=myBiomodModelOut, select.models='all',new.env=bc26,
                                     proj.name="2070_rcp26", binary.meth="TSS", build.clamping.mask=TRUE, compress=TRUE)
    Proj2070Rcp60<-BIOMOD_Projection(modeling.output=myBiomodModelOut, select.models='all',new.env=bc60,
                                     proj.name="2070_rcp60", binary.meth="TSS", build.clamping.mask=TRUE, compress=TRUE)
    Proj2070Rcp85<-BIOMOD_Projection(modeling.output=myBiomodModelOut, select.models='all',new.env=bc85,
                                     proj.name="2070_rcp85", binary.meth="TSS", build.clamping.mask=TRUE, compress=TRUE)
  },error=function(e)print(e)) 
}
cl<-makeCluster(120)
clusterExport(cl,"data.climate")
clusterExport(cl,"data.XY")
clusterExport(cl,"bc26")
clusterExport(cl,"bc60")
clusterExport(cl,"bc85")
clusterEvalQ(cl,library(biomod2))
parLapply(cl,biomod.spdis,MyBiomodSF)
stopCluster(cl)

##### current distribution data #####
setwd("D:/Shijia/J.APPL.E/Model")
library(biomod2)
library(stringr)
source("D:/Shijia/J.APPL.E/Modeldata/c_functions.R")
spnames<-list.files()
data.climate<-read.csv("D:/Shijia/J.APPL.E/Modeldata/current_climate.csv",row.names=1)
biomod.spdis<-get(load("D:/Shijia/J.APPL.E/Modeldata/biomod_spdis.RData"))
for(i in 1:length(spnames)){errin<-tryCatch({
  sp.names<-spnames[i]
  dat.current<-matrix(NA,23718,1)
  rownames(dat.current)<-rownames(data.climate)
  colnames(dat.current)<-sp.names
  model.out<-get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  prediction<-get(load(paste(sp.names,"/","proj_current","/","proj_current_",sp.names,".RData",sep="")))
  prediction<-as.data.frame(prediction)
  colnames(prediction)<-gsub("\\\\.","_",colnames(prediction))
  TSS<-get_evaluations(model.out,as.data.frame=TRUE)
  aa<-which(is.na(TSS$Testing.data)==TRUE|TSS$Testing.data<0.5)
  if(length(aa)!=0){TSS_1<-TSS[-aa,]}else{TSS_1<-TSS}
  TSS$Model.name<-as.character(TSS$Model.name)
  name<-TSS$Model.name[aa]
  if(length(name)!=0){sp<-which(colnames(prediction)%in%name==F)
  prediction.1<-prediction[,sp]}else{prediction.1<-prediction}
  pred.median<-apply(prediction.1,1,median,na.rm=TRUE)
  rocs<-roc(as.numeric(biomod.spdis[,sp.names]),pred.median,N = 200,plot=F)
  threshold<-threshold.tss(rocs)
  for(j in 1:23718){
    if(pred.median[j]<=threshold){dat.current[j,1]<-0}else{dat.current[j,1]<-1}}
  save(dat.current,file=paste("D:/Shijia/J.APPL.E/dis_current_1/dis_median/","proj_current_",sp.names,".RData",sep=""))
  rm(list=ls(pattern="out"))
  rm(list=ls(pattern="RUN"))
},error = function(e) print(e))}

##check errors#####
setwd("D:/Shijia/J.APPL.E/Model")
file1<-list.files()
setwd("D:/Shijia/J.APPL.E/dis_current/dis_median")
file2<-list.files()
spnames<-c()
for(i in 1:length(file2)){errin<-tryCatch({
  aa<- strsplit(file2[i],"_")[[1]][3]
  sp.names<-strsplit(aa,"\\\\.")[[1]][1]
  spnames<-c(spnames,sp.names)
},error = function(e) print(e))}

pos<-which(file1%in%spnames==F)

######add 200km buffer around predicted current ranges######
setwd("D:/Shijia/J.APPL.E/dis_current/dis_median")
sp<-list.files()
library(rgeos)
library(dismo)
biomod.spdis<-get(load("D:/Shijia/J.APPL.E/Modeldata/biomod_spdis.RData"))
data.sp<-as.data.frame(biomod.spdis)
coorxy<-read.csv("D:/Shijia/J.APPL.E/Modeldata/grid20km_xy.csv",row.names=1)
for(i in 1:length(sp)){errin<-tryCatch({
  aa<- strsplit(sp[i],"_")[[1]][3]
  sp.names<-strsplit(aa,"\\\\.")[[1]][1]
  pred.bin<-get(load(sp[i]))
  biomod.obs<-data.sp[sp.names]
  pred.buffer<-matrix(0,ncol=1,nrow=nrow(biomod.obs))
  rownames(pred.buffer)<-rownames(biomod.obs)
  colnames(pred.buffer)<-sp.names
  presence<-rownames(biomod.obs)[which(biomod.obs[,1]==1)]
  convex1<-coorxy[presence,4:5]
  ch<-convHull(convex1)
  p<-polygons(ch)
  b<-gBuffer(p,width=200)
  index.f<-rownames(pred.bin)[which(pred.bin[,1]==1)]
  xy.f<-coorxy[index.f,4:5]
  point.f<-SpatialPoints(xy.f)
  pf<-over(point.f,b)
  presence.final<-names(pf[which(!is.na(pf)==T)])
  pred.buffer[presence.final,1]<-1
  save(pred.buffer,file=paste("D:/Shijia/J.APPL.E/dis_current/dis_median_buffer/","pred_buffer_current_",sp.names,".RData",sep=""))
},error = function(e) print(e))}

######full current distribution data#######
setwd("D:/Shijia/J.APPL.E/dis_future/proj_2070_rcp8.5_buffer")
biomod.spdis8732<-get(load("D:/Shijia/J.APPL.E/Modeldata/biomod_spdis.RData"))
sp<-list.files()
dis.sp<-c()
for(i in 1:length(sp)){dis.sp.t<-get(load(sp[i]))
dis.sp<-cbind(dis.sp,dis.sp.t)}
data.sp<-dis.sp[,match(colnames(biomod.spdis8732),colnames(dis.sp))]
save(data.sp,file="proj_2070_rcp2.6_buffer.RData")


#######future distribution data######full dispersal#####
library(biomod2)
library(stringr)
setwd("D:/Shijia/J.APPL.E/Model")
source("D:/Shijia/J.APPL.E/Modeldata/c_functions.R")
biomod.spdis<-get(load("D:/Shijia/J.APPL.E/Modeldata/biomod_spdis.RData"))
spnames<-list.files()
data.climate<-read.csv("D:/Shijia/J.APPL.E/Modeldata/current_climate.csv",row.names=1)
for(i in 1:length(spnames)){errin<-tryCatch({
  sp.names<-spnames[i]
  data.2070.rcp8.5<-matrix(NA,23718,1)
  rownames(data.2070.rcp8.5)<-rownames(data.climate)
  colnames(data.2070.rcp8.5)<-sp.names
  prediction<-get(load(paste(sp.names,"/","proj_current","/","proj_current_",sp.names,".RData",sep="")))
  prediction<-as.data.frame(prediction)
  model.out<-get(load(paste(sp.names,"/",sp.names,".",sp.names,"Firstmodeling.models.out",sep='')))
  prediction.future<-get(load(paste(sp.names,"/","proj_2070_rcp85","/","proj_2070_rcp85_",sp.names,".RData",sep="")))
  prediction.future<-as.data.frame(prediction.future)
  colnames(prediction)<-gsub("\\\\.","_",colnames(prediction))
  colnames(prediction.future)<-gsub("\\\\.","_",colnames(prediction.future))
  TSS<-get_evaluations(model.out,as.data.frame=TRUE)
  aa<-which(is.na(TSS$Testing.data)==TRUE|TSS$Testing.data<0.5)
  if(length(aa)!=0){TSS_1<-TSS[-aa,]}else{TSS_1<-TSS}
  TSS$Model.name<-as.character(TSS$Model.name)
  name<-TSS$Model.name[aa]
  if(length(name)!=0){sp<-which(colnames(prediction)%in%name==F)
  prediction.1<-prediction[,sp]}else{prediction.1<-prediction}
  if(length(name)!=0){sp<-which(colnames(prediction.future)%in%name==F)
  prediction.1.future<-prediction.future[,sp]}else{prediction.1.future<-prediction.future}
  pred.median<-apply(prediction.1,1,median,na.rm=TRUE)
  rocs<-roc(as.numeric(biomod.spdis[,sp.names]),pred.median,N = 200,plot=F)
  threshold<-threshold.tss(rocs)
  pred.median.future<-apply(prediction.1.future,1,median,na.rm=TRUE)
  for(j in 1:23718){
    if(pred.median.future[j]<=threshold){data.2070.rcp8.5[j,1]<-0}else{data.2070.rcp8.5[j,1]<-1}}
  save(data.2070.rcp6.0,file=paste("D:/Shijia/J.APPL.E/dis_future_1/proj_2070_rcp6.0/","proj_2070_rcp6.0_",sp.names,".RData",sep=""))
  rm(list=ls(pattern="out"))
  rm(list=ls(pattern="RUN"))
},error = function(e) print(e))}

#######future distribution data######20km/decade dispersal#####  
setwd("D:/Shijia/J.APPL.E/dis_future/proj_2070_rcp8.5")
sp<-list.files()
library(rgeos)
library(dismo)
biomod.spdis<-get(load("D:/Shijia/J.APPL.E/distribution/proj_current_buffer_sp.RData"))
data.sp<-as.data.frame(biomod.spdis)
coorxy<-read.csv("D:/Shijia/J.APPL.E/Modeldata/grid20km_xy.csv",row.names=1)
for(i in 1:length(sp)){errin<-tryCatch({
aa<- strsplit(sp[i],"_")[[1]][4]
sp.names<-strsplit(aa,"\\\\.")[[1]][1]
pred.bin<-get(load(sp[i]))
biomod.obs<-data.sp[sp.names]
pred.buffer<-matrix(0,ncol=1,nrow=nrow(biomod.obs))
rownames(pred.buffer)<-rownames(biomod.obs)
colnames(pred.buffer)<-sp.names
presence<-rownames(biomod.obs)[which(biomod.obs[,1]==1)]
convex1<-coorxy[presence,4:5]
ch<-convHull(convex1)
p<-polygons(ch)
b<-gBuffer(p,width=200)
index.f<-rownames(pred.bin)[which(pred.bin[,1]==1)]
if(length(index.f)==0){pred.buffer[,1]<-0}else{
xy.f<-coorxy[index.f,4:5]
point.f<-SpatialPoints(xy.f)
pf<-over(point.f,b)
presence.final<-names(pf[which(!is.na(pf)==T)])
pred.buffer[presence.final,1]<-1}
save(pred.buffer,file=paste("D:/Shijia/J.APPL.E/dis_future/","proj_2070_rcp8.5_buffer/","pred_buffer_2070_rcp8.5_",sp.names,".RData",sep=""))
},error = function(e) print(e))}


#######future distribution data######no dispersal##### 
setwd("D:/Shijia/J.APPL.E/distribution data")
dis_current<-get(load("proj_current_buffer_sp.RData"))
dis_future<-get(load("future_dis_fully/proj_2070_rcp8.5_sp.RData"))
for(i in 1:8732){
  aa<-which(dis_current[,i]==1)
  bb<-which(dis_future[,i]==1)
  res.1<-bb[which(bb%in%aa==F)]
  dis_future[res.1,i]<-0
}

save(dis_future,file="D:/Shijia/J.APPL.E/distribution/future_dis_stable/proj_2070_rcp8.5.RData")

#####Species are classified into different threat level####
####Full dispersal#####
setwd("E:/PA/J.appl.eco/Revised/distribution")
dis_current<-get(load("proj_current_buffer_sp.RData"))
dis_rcp2.6<-get(load("future_dis_fully/proj_2070_rcp2.6_sp.RData"))
dis_rcp6.0<-get(load("future_dis_fully/proj_2070_rcp6.0_sp.RData"))
dis_rcp8.5<-get(load("future_dis_fully/proj_2070_rcp8.5_sp.RData"))
future_full<-as.data.frame(matrix(NA,8732,5))
colnames(future_full)<-c("species","current","RCP2.6","RCP6.0","RCP8.5")
future_full$species<-colnames(dis_current)
res.current<-apply(dis_current,2,sum)
res.rcp2.6<-apply(dis_rcp2.6,2,sum)
res.rcp6.0<-apply(dis_rcp6.0,2,sum)
res.rcp8.5<-apply(dis_rcp8.5,2,sum)
future_full[,2]<-res.current
future_full[,3]<-res.rcp2.6
future_full[,4]<-res.rcp6.0
future_full[,5]<-res.rcp8.5

###20km/decade###
setwd("E:/PA/J.appl.eco/Revised/distribution")
dis_current<-get(load("proj_current_buffer_sp.RData"))
dis_rcp2.6<-get(load("future_dis_20km_buffer/proj_2070_rcp2.6_buffer.RData"))
dis_rcp6.0<-get(load("future_dis_20km_buffer/proj_2070_rcp6.0_buffer.RData"))
dis_rcp8.5<-get(load("future_dis_20km_buffer/proj_2070_rcp8.5_buffer.RData"))
future_buffer<-as.data.frame(matrix(NA,8732,5))
colnames(future_buffer)<-c("species","current","RCP2.6","RCP6.0","RCP8.5")
future_buffer$species<-colnames(dis_current)
res.current<-apply(dis_current,2,sum)
res.rcp2.6<-apply(dis_rcp2.6,2,sum)
res.rcp6.0<-apply(dis_rcp6.0,2,sum)
res.rcp8.5<-apply(dis_rcp8.5,2,sum)
future_buffer[,2]<-res.current
future_buffer[,3]<-res.rcp2.6
future_buffer[,4]<-res.rcp6.0
future_buffer[,5]<-res.rcp8.5

##no dispersal##
setwd("E:/PA/J.appl.eco/Revised/distribution")
dis_current<-get(load("proj_current_buffer_sp.RData"))
dis_rcp2.6<-get(load("future_dis_stable/proj_2070_rcp2.6.RData"))
dis_rcp6.0<-get(load("future_dis_stable/proj_2070_rcp6.0.RData"))
dis_rcp8.5<-get(load("future_dis_stable/proj_2070_rcp8.5.RData"))
future_stable<-as.data.frame(matrix(NA,8732,5))
colnames(future_stable)<-c("species","current","RCP2.6","RCP6.0","RCP8.5")
future_stable$species<-colnames(dis_current)
res.current<-apply(dis_current,2,sum)
res.rcp2.6<-apply(dis_rcp2.6,2,sum)
res.rcp6.0<-apply(dis_rcp6.0,2,sum)
res.rcp8.5<-apply(dis_rcp8.5,2,sum)
future_stable[,2]<-res.current
future_stable[,3]<-res.rcp2.6
future_stable[,4]<-res.rcp6.0
future_stable[,5]<-res.rcp8.5


#########Figure 1#######
library(cowplot)
library(ggplot2)
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis")
data.2070.full<-read.csv("data.2070.full.analysis.csv")
data.2070.buffer<-read.csv("data.2070.buffer.analysis.csv")
data.2070.stable<-read.csv("data.2070.stable.analysis.csv")
data.2070.full$Type<-factor(data.2070.full$Type,levels=c("EX","CE","EN","VU","LR"))
data.2070.buffer$Type<-factor(data.2070.buffer$Type,levels=c("EX","CE","EN","VU","LR"))
data.2070.stable$Type<-factor(data.2070.stable$Type,levels=c("EX","CE","EN","VU","LR"))
windowsFonts()
plot.2070.full<-ggplot(data.2070.full,aes(x=factor(RCP),y=Value,fill=Type))+geom_bar(stat = 'identity', position = 'stack')+scale_fill_brewer(palette = 'RdYlBu')+theme_bw(base_family="serif",base_size=21)+scale_y_continuous(breaks = seq(0,1,0.25),limits=c(0,1))+theme(axis.title=element_blank(),legend.position="none")
plot.2070.buffer<-ggplot(data.2070.buffer,aes(x=factor(RCP),y=Value,fill=Type))+geom_bar(stat = 'identity', position = 'stack')+scale_fill_brewer(palette = 'RdYlBu')+theme_bw(base_family="serif",base_size=21)+scale_y_continuous(breaks = seq(0,1,0.25),limits=c(0,1))+theme(axis.title=element_blank(),legend.position="none")
plot.2070.stable<-ggplot(data.2070.stable,aes(x=factor(RCP),y=Value,fill=Type))+geom_bar(stat = 'identity', position = 'stack')+scale_fill_brewer(palette = 'RdYlBu')+theme_bw(base_family="serif",base_size=21)+scale_y_continuous(breaks = seq(0,1,0.25),limits=c(0,1))+theme(axis.title=element_blank(),legend.position="none")

plot<-plot_grid(plotlist = list(plot.2070.full,plot.2070.buffer,plot.2070.stable),ncol=3,nrow=1,align="hv")
ggsave("plot.tiff",width=13,height=10,dpi=600)

####geographical patterns in species richness of threatened species#####
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis")
data<-read.csv("E:/PA/J.appl.eco/Revised/future_buffer.csv") 
res<-which(data$RCP8.5_CSH<=-0.3)
biomod.spdis<-get(load("E:/PA/J.appl.eco/Revised/distribution/future_dis_20km_buffer/proj_2070_rcp8.5_buffer.RData"))
###biomod.spdis<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
species<-data$species[res]
biomod.spdis.buffer.rcp85<-biomod.spdis[,species]
SR.buffer.rcp85<-as.data.frame(matrix(NA,23718,1))
rownames(SR.buffer.rcp85)<-rownames(biomod.spdis.buffer.rcp85)
colnames(SR.buffer.rcp85)<-"SR"
SR.buffer.rcp85[,1]<-apply(biomod.spdis.buffer.rcp85,1,sum)
write.csv(SR.buffer.rcp85,file="buffer.rcp85.csv")

###### PD calcuclation #####
library(ape)
library(phytools)
library(picante)
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis/PD")
dis_current<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
future_full<-read.csv("E:/PA/J.appl.eco/Revised/future_full.csv")
future_full_rcp8.5<-future_full[which(future_full$RCP8.5_CSH<=-0.3),]
sp<-as.character(future_full_rcp8.5$species)
dis_current_TS<-dis_current[,sp]
names_match<-read.csv("E:/PA/J.appl.eco/Revised/Modeldata/names_match.csv")
names_match_1<-names_match[which(names_match$code%in%sp==T),]
colnames(dis_current_TS)<-names_match_1$name

#####
tree.100<-get(load("E:/PA/J.appl.eco/Revised/Modeldata/Random100.RData"))
tree.TS.full.8.5<-list()
for( i in 1:100){
  tree<-tree.100[[i]]
  postemp <- which(!tree$tip.label %in% colnames(dis_current_TS))
  tree2 <- drop.tip(tree, tree$tip.label[postemp]) 
  tree.TS.full.8.5[[i]]<-tree2
}
save(tree.TS.full.8.5,file="tre_full_rcp85.RData")

######
tree.TS.full.8.5<-get(load("tre_full_rcp85.RData"))
dis_future<-get(load("E:/PA/J.appl.eco/Revised/distribution/future_dis_fully/proj_2070_rcp8.5_sp.RData"))
future_full<-read.csv("E:/PA/J.appl.eco/Revised/future_full.csv")
future_full_rcp8.5<-future_full[which(future_full$RCP8.5_CSH<=-0.3),]
sp<-as.character(future_full_rcp8.5$species)
dis_future_TS<-dis_future[,sp]
names_match<-read.csv("E:/PA/J.appl.eco/Revised/Modeldata/names_match.csv")
names_match_1<-names_match[which(names_match$code%in%sp==T),]
colnames(dis_future_TS)<-names_match_1$name

######
##spnames<-colnames(dis_current_TS)[which(colnames(dis_current_TS)%in%tree.TS.buffer.8.5[[1]]$tip.label==T)]
##dis_current_TS_1<-dis_current_TS[,spnames]

##spnames<-colnames(dis_future_TS)[which(colnames(dis_future_TS)%in%tree.TS.buffer.8.5[[1]]$tip.label==T)]
##dis_future_TS_1<-dis_future_TS[,spnames]

####
library(parallel)
pd.fun<-function(x) pd1<-pd(spdis,x,include.root=TRUE)
spdis<-dis_current_TS_1
cl<-makeCluster(120)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(picante))
pd.obs<-parLapply(cl,tree.TS.full.8.5,pd.fun)
stopCluster(cl)
save(pd.obs,file="Full/pd.obs.full.rcp85.current.RData")


##########
library(parallel)
pd.fun<-function(x) pd1<-pd(spdis,x,include.root=TRUE)
spdis<-dis_future_TS_1
cl<-makeCluster(120)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(picante))
pd.obs<-parLapply(cl,tree.TS.full.8.5,pd.fun)
stopCluster(cl)
save(pd.obs,file="Full/pd.obs.full.rcp85.future.RData")

######mean current PD#####
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis/PD")
data.buffer.current<-get(load("Buffer/pd.obs.buffer.rcp85.current.RData"))
data.pd.current<-as.data.frame(matrix(NA,23718,3))
rownames(data.pd.current)<-rownames(data.buffer.current[[1]])
colnames(data.pd.current)<-c("PD.current","SR.current","dis.current")
data.pd.current[,2]<-data.buffer.current[[1]]$SR
PD.1<-c()
for (i in 1:100){
  dat.pd<-data.buffer.current[[i]]
  PD.1<-cbind(PD.1,dat.pd$PD)
}
res<-apply(PD.1,1,mean)
data.pd.current[,1]<-res
PD.1.Z<-scale(data.pd.current$PD.current)
SR.1.Z<-scale(data.pd.current$SR.current)
data.pd.current[,3]<-PD.1.Z-SR.1.Z
write.csv(data.pd.current,file="Buffer/PD.rcp85.current.csv")

######mean future PD#####
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis/PD")
data.buffer.future<-get(load("Buffer/pd.obs.buffer.rcp85.future.RData"))
data.pd.future<-as.data.frame(matrix(NA,23718,3))
rownames(data.pd.future)<-rownames(data.buffer.future[[1]])
colnames(data.pd.future)<-c("PD.future","SR.future","dis.future")
data.pd.future[,2]<-data.buffer.future[[1]]$SR
#########
PD.future<-c()
for (i in 1:100){
  dat.pd<-data.buffer.future[[i]]
  PD.future<-cbind(PD.future,dat.pd$PD)
}
res<-apply(PD.future,1,mean)
data.pd.future[,1]<-res
PD.1.Z<-scale(data.pd.future$PD.future)
SR.1.Z<-scale(data.pd.future$SR.future)
data.pd.future[,3]<-PD.1.Z-SR.1.Z
write.csv(data.pd.future,file="Buffer/PD.rcp85.future.csv")

######null model#######
#####Random value######species pool=all#######
setwd("E:/PA/J.appl.eco/Revised")
tree<-get(load("Modeldata/Random100.RData"))
data<-get(load("PA/proj_current_buffer_sp.RData"))
names_match<-read.csv("Modeldata/names_match.csv")
sp.list<-names_match$name
colnames(data)<-sp.list
sp.pool<-as.matrix(sp.list)

########
RS.list<-as.matrix(unique(apply(data,1,sum)))
PD.ii<-1:length(sp.pool)
spnumber<-length(sp.pool)
PD.list <- list()

for (i in 1:length(RS.list))
{
  PD.ma <- matrix(0,ncol=spnumber, nrow =99);rownames(PD.ma) <- c(1:99);colnames(PD.ma)  <- sp.pool
  RS <- as.numeric(RS.list[i])
  for (j in 1:99)
  {
    sp.pos <- sample(PD.ii, size = RS)
    PD.ma[j,sp.pos] <-1
  }
  PD.list[[i]] <- PD.ma;print(i)
}
save(PD.list,file="PD.list_8732.RData")

########
library(parallel)
PD.RandomResult.list <- list()
mycl<-makeCluster(100)
source("Myfunctions.R")
source("read_Phylogenetic_Tree.R")
source("joint.table.R")
for (m in 1:100){
  ##calculate PD in parallel
  ss1 <-parLapply(cl=mycl, X=PD.list, fun=pd.phybl, tree=tre.species[[m]] ,format="matrix",  method="branch.length");  date()
  ##calculate mean and sd for each matrix
  PD <-  matrix(0,ncol=2, nrow = length(RS.list)); colnames(PD) <- c("mean", "sd")
  for (i in 1:length(RS.list))
  {
    PD[i,] <- c(mean(ss1[[i]], na.rm=T),sd(ss1[[i]], na.rm=T))
  }
  rownames(PD) <- RS.list
  PD.RandomResult.list[[m]] <- PD
  print(m);date()
}
stopCluster(mycl)
save(PD.RandomResult.list, file="PD_RandomResult_list.RData")

#######phylogenetic signal#####
setwd("D:/Shijia/J.APPL.E/phylogenetic_signal")
names_match<-read.csv("D:/Shijia/J.APPL.E/Modeldata/names_match.csv")  
future_full<-read.csv("D:/Shijia/J.APPL.E/future_full.csv")
future_full_rcp2.6<-future_full[which(future_full$RCP2.6_CSH<=-0.3),]
sp<-as.character(future_full_rcp2.6$species)  
data.phylo<-as.data.frame(matrix(0,8732,3))
colnames(data.phylo)<-c("sp","code","value")
data.phylo$sp<-names_match$name
data.phylo$code<-names_match$code
aa<-which(names_match$code%in%sp==T)
data.phylo[aa,3]<-1
tree.100<-get(load("D:/Shijia/J.APPL.E/Modeldata/Random100.RData"))
tre<-tree.100[[1]]
data.phylo.1<-data.phylo[which(data.phylo$sp%in%tre$tip.label==T),]
data.final<-data.phylo.1[match(tre$tip.label,data.phylo.1$sp),]

###########
library(parallel)
physig.fun<-function(x) phy.sig<-phylo.d(spdis,x,sp,value,permut=999)
spdis<-data.final
cl<-makeCluster(120)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(caper))
physig.obs<-parLapply(cl,tree.100,physig.fun)
stopCluster(cl)

######
setwd("E:/PA/J.appl.eco/Revised/endanger_analysis/phylogenetic_signal")
data<-get(load("physig.stable.rcp85.RData"))
phylo<-data.frame(matrix(NA,100,4))
colnames(phylo)<-c("tree.num","phylo.d","p-value0","p-value1")
phylo$tree.num<-as.character(1:100)
for(i in 1:100){
  phylo[i,2]<-data[[i]]$DEstimate
  phylo[i,3]<-data[[i]]$Pval0
  phylo[i,4]<-data[[i]]$Pval1
}
subset(phylo,phylo$phylo.d==sort(phylo$phylo.d)[50])  

####### estimation of the number of threatened species within PAs ######
setwd("E:/PA/J.appl.eco/Revised")
data_distri<-get(load("distribution/proj_current_buffer_sp.RData"))
protected_data<-read.csv("E:/PA/J.appl.eco/Revised/PA_analysis/30%threshold.csv")
data_distri_protect<-data_distri[which(rownames(data_distri)%in%protected_data$Gridcode==T),]
sp.protected<-c()
for(i in 1:nrow(data_distri_protect)){
  aa<-which(data_distri_protect[i,]==1)
  aa.1<-colnames(data_distri_protect)[aa]
  sp.protected<-c(sp.protected,aa.1)
}
sp.protected<-unique(sp.protected)

################
future_full<-read.csv("E:/PA/J.appl.eco/Revised/future_full.csv")
sp<-as.character(future_full$species)
sp.EX<-future_full$species[which(future_full$Type_1=="EX")]
length(sp.EX)
length(which(sp.protected%in%sp.EX==T))

######Proportion of PD(%) being protected #####

library(ape)
library(phytools)
library(picante)
setwd("E:/PA/J.appl.eco/Revised/PA_analysis")
future_full<-read.csv("E:/PA/J.appl.eco/Revised/future_full.csv")
future_full_rcp8.5<-future_full[which(future_full$RCP8.5_CSH<=-0.3),]
sp<-as.character(future_full_rcp8.5$species)
data_distri<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
data_distri_TS<-data_distri[,sp]
protected_data<-read.csv("30%threshold.csv")
data_distri_protect<-data_distri_TS[which(rownames(data_distri_TS)%in%protected_data$Gridcode==T),]
sp.protected<-c()
for(i in 1:nrow(data_distri_protect)){
  aa<-which(data_distri_protect[i,]==1)
  aa.1<-colnames(data_distri_protect)[aa]
  sp.protected<-c(sp.protected,aa.1)
}
sp.protected<-unique(sp.protected)

######
TR.matrix<-matrix(1,1,length(sp))
colnames(TR.matrix)<-sp
names_match<-read.csv("E:/PA/J.appl.eco/Revised/Modeldata/names_match.csv")
names_match_1<-names_match[which(names_match$code%in%sp==T),]
colnames(TR.matrix)<-names_match_1$name

#####
tree.TS.full.8.5<-get(load("E:/PA/J.appl.eco/Revised/endanger_analysis/PD/Full/tre_full_rcp85.RData"))
spnames<-colnames(TR.matrix)[which(colnames(TR.matrix) %in% tree.TS.full.8.5[[1]]$tip.label==T)]
TR.matrix.1<-TR.matrix[,spnames,drop=FALSE]

####
pd<-c()
for(i in 1:100){
  tree<-tree.TS.full.8.5[[i]]
  pd.1<-pd(TR.matrix.1,tree,include.root = TRUE)$PD
  pd<-c(pd,pd.1)
}
pd.total<-mean(pd)

######
TR.PR.matrix<-matrix(0,1,length(sp))
colnames(TR.PR.matrix)<-sp
TR.PR.matrix[,sp.protected]<-1
names_match<-read.csv("E:/PA/J.appl.eco/Revised/Modeldata/names_match.csv")
names_match_1<-names_match[which(names_match$code%in%sp==T),]
colnames(TR.PR.matrix)<-names_match_1$name

#######
spnames.protected<-colnames(TR.PR.matrix)[which(colnames(TR.PR.matrix) %in% tree.TS.full.8.5[[1]]$tip.label==T)]
TR.PR.matrix.1<-TR.PR.matrix[,spnames.protected,drop=FALSE]

######
pd.PR<-c()
for(i in 1:100){
  tree<-tree.TS.full.8.5[[i]]
  pd.1<-pd(TR.PR.matrix.1,tree,include.root = TRUE)$PD
  pd.PR<-c(pd.PR,pd.1)
}
pd.PR.total<-mean(pd.PR)
pd.PR.total/pd.total


######Changes in the number of threatened species within existing PA networks#####
setwd("E:/PA/J.appl.eco/Revised/PA_analysis")
threshold<-read.csv("30%threshold.csv")
data_distri<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
data_distri_full_rcp2.6<-get(load("E:/PA/J.appl.eco/Revised/distribution/future_dis_fully/proj_2070_rcp2.6_sp.RData"))
Code<-unique(as.character(threshold$PAcode))
length(Code)
future_stable<-read.csv("E:/PA/J.appl.eco/Revised/future_stable.csv")
future_stable_rcp2.6<-future_stable[which(future_stable$RCP2.6_CSH<=-0.3),]
sp<-as.character(future_stable_rcp2.6$species)
data.sp.current<-data_distri[,sp]
data.sp.future<-data_distri_stable_rcp2.6[,sp]
data.protect.SR<-as.data.frame(matrix(NA,length(Code),3))
colnames(data.protect.SR)<-c("PAcode","current_SR","future_SR")
data.protect.SR$PAcode<-Code
for(i in 1:length(Code)){
  aa<-threshold[which(threshold$PAcode==Code[i]),]
  sp.data.current<-data.sp.current[which(rownames(data.sp.current)%in%aa$Gridcode==T),,drop=FALSE]
  sp.data.future<-data.sp.future[which(rownames(data.sp.future)%in%aa$Gridcode==T),,drop=FALSE]
  sp.current<-c()
  sp.future<-c()
  for(j in 1:nrow(aa)){
    bb<-colnames(sp.data.current)[which(sp.data.current[j,]==1)]
    cc<-colnames(sp.data.future)[which(sp.data.future[j,]==1)]
    sp.current<-unique(c(sp.current,bb))
    sp.future<-unique(c(sp.future,cc))
  }
  data.protect.SR[i,2]<-length(sp.current)
  data.protect.SR[i,3]<-length(sp.future)
}
write.csv(data.protect.SR,file="SR/Stable/SR.rcp26.csv")

########Changes in PD of threatened species within existing PA networks########
####making distribution data for each PA######
setwd("E:/PA/J.appl.eco/Revised/PA_analysis")
threshold<-read.csv("30%threshold.csv")
data_distri<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
data_distri_stable_rcp8.5<-get(load("E:/PA/J.appl.eco/Revised/distribution/future_dis_stable/proj_2070_rcp8.5.RData"))
Code<-unique(as.character(threshold$PAcode))
future_stable<-read.csv("E:/PA/J.appl.eco/Revised/future_stable.csv")
future_stable_rcp8.5<-future_stable[which(future_stable$RCP8.5_CSH<=-0.3),]
sp<-as.character(future_stable_rcp8.5$species)
data.sp.current<-data_distri[,sp]
data.sp.future<-data_distri_stable_rcp8.5[,sp]
data.protect.current<-matrix(0,length(Code),length(sp))
data.protect.future<-matrix(0,length(Code),length(sp))
rownames(data.protect.current)<-Code
rownames(data.protect.future)<-Code
names_match<-read.csv("E:/PA/J.appl.eco/Revised/Modeldata/names_match.csv")
names_match_1<-names_match[which(names_match$code%in%sp==T),]
colnames(data.protect.current)<-names_match_1$name
colnames(data.protect.future)<-names_match_1$name
colnames(data.sp.current)<-names_match_1$name
colnames(data.sp.future)<-names_match_1$name

#######
for(i in 1:length(Code)){
  aa<-threshold[which(threshold$PAcode==Code[i]),]
  sp.data.current<-data.sp.current[which(rownames(data.sp.current)%in%aa$Gridcode==T),,drop=FALSE]
  sp.data.future<-data.sp.future[which(rownames(data.sp.future)%in%aa$Gridcode==T),,drop=FALSE]
  sp.current<-c()
  sp.future<-c()
  for(j in 1:nrow(aa)){
    bb<-colnames(sp.data.current)[which(sp.data.current[j,]==1)]
    cc<-colnames(sp.data.future)[which(sp.data.future[j,]==1)]
    sp.current<-unique(c(sp.current,bb))
    sp.future<-unique(c(sp.future,cc))
  }
  data.protect.current[i,sp.current]<-1
  data.protect.future[i,sp.future]<-1
}
save(data.protect.current,file="distribution/stable.rcp85.current.RData")
save(data.protect.future,file="distribution/stable.rcp85.future.RData")

######
PD.current<-get(load("E:/PA/J.appl.eco/Revised/PA_analysis/distribution/full.rcp26.current.RData"))
PD.future<-get(load("E:/PA/J.appl.eco/Revised/PA_analysis/distribution/full.rcp26.future.RData"))
PD.rcp26<-as.data.frame(matrix(NA,485,3))
colnames(PD.rcp26)<-c("PDcode","current_PD","future_PD")
PD.rcp26$PDcode<-rownames(PD.current)
tree<-get(load("E:/PA/J.appl.eco/Revised/endanger_analysis/PD/Full/tre_full_rcp26.RData"))

######
library(parallel)
pd.fun<-function(x) pd1<-pd(spdis,x,include.root=TRUE)
spdis<-PD.current
cl<-makeCluster(100)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(picante))
pd.obs<-parLapply(cl,tree,pd.fun)
stopCluster(cl)
save(pd.obs,file="E:/PA/J.appl.eco/Revised/PA_analysis/PD/Full/pd.current.rcp26.RData")


##########
library(parallel)
pd.fun<-function(x) pd1<-pd(spdis,x,include.root=TRUE)
spdis<-PD.future
cl<-makeCluster(100)
clusterExport(cl = cl, varlist = c("spdis"))
clusterEvalQ(cl = cl, library(picante))
pd.obs<-parLapply(cl,tree,pd.fun)
stopCluster(cl)
save(pd.obs,file="E:/PA/J.appl.eco/Revised/PA_analysis/PD/Full/pd.future.rcp26.RData")

######mean current PD for each PA#####
setwd("E:/PA/J.appl.eco/Revised/PA_analysis/PD")
data.current<-get(load("Stable/pd.current.rcp85.RData"))
data.pd.current<-as.data.frame(matrix(NA,485,1))
rownames(data.pd.current)<-rownames(data.current[[1]])
colnames(data.pd.current)<-c("PD.current")
PD.1<-c()
for (i in 1:100){
  dat.pd<-data.current[[i]]
  PD.1<-cbind(PD.1,dat.pd$PD)
}
res<-apply(PD.1,1,mean)
data.pd.current[,1]<-res
write.csv(data.pd.current,file="Stable/PD.rcp85.current.csv")

######mean future PD#####
setwd("E:/PA/J.appl.eco/Revised/PA_analysis/PD")
data.future<-get(load("pd.obs.current.RData"))
data.pd.future<-as.data.frame(matrix(NA,485,1))
rownames(data.pd.future)<-rownames(data.future[[1]])
colnames(data.pd.future)<-c("PD.future")
PD.future<-c()
for (i in 1:100){
  dat.pd<-data.future[[i]]
  PD.future<-cbind(PD.future,dat.pd$PD)
}
res<-apply(PD.future,1,mean)
data.pd.future[,1]<-res
write.csv(data.pd.future,file="Stable/PD.rcp85.future.csv")

####

######Figure 3#####
setwd("E:/PA/J.appl.eco/Revised/PA_analysis")
library(ggplot2)
library(cowplot)
windowsFonts()

######SR######
data<-read.csv("data.PA.SR.6.0.csv")
data.full<-data[which(data$scenario=="full"),]
data.full$group<-factor(data.full$group,levels = c("current","2070"))
plot.full<-ggplot(data=data.full,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.text.x = element_blank(),plot.title=element_text(hjust=0.5,size=20),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,45,9),expand = c(0,0),limits=c(0,45))+scale_fill_manual(values=c("#9900FF","#CC66FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)
data.buffer<-data[which(data$scenario=="20km/decade"),]
data.buffer$group<-factor(data.buffer$group,levels = c("current","2070"))
plot.buffer<-ggplot(data=data.buffer,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.text.x = element_blank(),plot.title=element_text(hjust=0.5,size=20),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,75,15),expand = c(0,0),limits=c(0,75))+scale_fill_manual(values=c("#9900FF","#CC66FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)
data.stable<-data[which(data$scenario=="no"),]
data.stable$group<-factor(data.stable$group,levels = c("current","2070"))
plot.stable<-ggplot(data=data.stable,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.text.x = element_blank(),plot.title=element_text(hjust=0.5,size=20),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,180,36),expand = c(0,0),limits=c(0,180))+scale_fill_manual(values=c("#9900FF","#CC66FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)
  
######PD######
data.PD<-read.csv("data.PA.PD.6.0.csv")
data.PD.full<-data.PD[which(data.PD$scenario=="full"),]
data.PD.full$group<-factor(data.PD.full$group,levels = c("current","2070"))
plot.PD.full<-ggplot(data=data.PD.full,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,2700,540),expand = c(0,0),limits=c(0,2700))+scale_fill_manual(values=c("#CC00CC","#FF99FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)
data.PD.buffer<-data.PD[which(data.PD$scenario=="20km/decade"),]
data.PD.buffer$group<-factor(data.PD.buffer$group,levels = c("current","2070"))
plot.PD.buffer<-ggplot(data=data.PD.buffer,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,4000,800),expand = c(0,0),limits=c(0,4000))+scale_fill_manual(values=c("#CC00CC","#FF99FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)
data.PD.stable<-data.PD[which(data.PD$scenario=="no"),]
data.PD.stable$group<-factor(data.PD.stable$group,levels = c("current","2070"))
plot.PD.stable<-ggplot(data=data.PD.stable,aes(x=group,y=mean,fill=group))+geom_bar(stat="identity")+theme_classic(base_family="serif",base_size=21)+theme(axis.title=element_blank(),axis.ticks.x = element_blank(),legend.position = "none")+scale_y_continuous(breaks=seq(0,7500,1500),expand = c(0,0),limits=c(0,7500))+scale_fill_manual(values=c("#CC00CC","#FF99FF"),labels=c("current","2070"))+geom_errorbar(aes(ymax=high_Cl,ymin=low_Cl,width=0.1),size=1.0)

######
plot<-plot_grid(plotlist = list(plot.full,plot.buffer,plot.stable,plot.PD.full,plot.PD.buffer,plot.PD.stable),ncol=3,nrow=2,align="hv")
ggsave("plot.6.0.tiff",width=12,height=10,dpi=600)

####the contributions of different ecological factors to the effectiveness of PAs#####
#####Random forest######
library(randomForest) 
setwd("E:/PA/J.appl.eco/Revised/PA_analysis/SR")
data.PD<-read.csv("Full/SR.rcp26.csv")
res<-randomForest(SR_dis~Area+HP+range+altitude+vegetation_loss,data=data.SR,ntree=2000, nodesize=5,importance=TRUE)
res$importance

######box plot#######
library(ggplot2)
library(cowplot)
setwd("E:/PA/J.appl.eco/Revised2/boxplot")
data.SR.full<-read.csv("SR.rcp85.full.csv")
data.SR.buffer<-read.csv("SR.rcp85.buffer.csv")
data.PD.full<-read.csv("PD.rcp85.full.csv")
data.PD.buffer<-read.csv("PD.rcp85.buffer.csv")

####
data.SR.full.1<-data.SR.full[which(data.SR.full$Type!="NA"),]
data.SR.buffer.1<-data.SR.buffer[which(data.SR.buffer$Type!="NA"),]
data.PD.full.1<-data.PD.full[which(data.PD.full$Type!="NA"),]
data.PD.buffer.1<-data.PD.buffer[which(data.PD.buffer$Type!="NA"),]
data.SR.full.1$vegetation_loss<-as.numeric(data.SR.full.1$vegetation_loss)
data.SR.buffer.1$vegetation_loss<-as.numeric(data.SR.buffer.1$vegetation_loss)
data.PD.full.1$vegetation_loss<-as.numeric(data.PD.full.1$vegetation_loss)
data.PD.buffer.1$vegetation_loss<-as.numeric(data.PD.buffer.1$vegetation_loss)
####kruskal.test(SR_dis~class,data=data.SR.85)

windowsFonts()

####
plot.SR.full.range<-ggplot(data.SR.full.1)+geom_boxplot(aes(x=Type,y=range,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF6699","#FF99CC"),labels=c("Negative","Positive"))
plot.SR.buffer.range<-ggplot(data.SR.buffer.1)+geom_boxplot(aes(x=Type,y=range,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF6699","#FF99CC"),labels=c("Negative","Positive"))
plot.PD.full.range<-ggplot(data.PD.full.1)+geom_boxplot(aes(x=Type,y=range,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF6699","#FF99CC"),labels=c("Negative","Positive"))
plot.PD.buffer.range<-ggplot(data.PD.buffer.1)+geom_boxplot(aes(x=Type,y=range,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF6699","#FF99CC"),labels=c("Negative","Positive"))

####
plot.SR.full.altitude<-ggplot(data.SR.full.1)+geom_boxplot(aes(x=Type,y=altitude,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF9900","#FFCC00"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,5000,1250),limits=c(0,5000))
plot.SR.buffer.altitude<-ggplot(data.SR.buffer.1)+geom_boxplot(aes(x=Type,y=altitude,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF9900","#FFCC00"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,5000,1250),limits=c(0,5000))
plot.PD.full.altitude<-ggplot(data.PD.full.1)+geom_boxplot(aes(x=Type,y=altitude,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF9900","#FFCC00"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,5000,1250),limits=c(0,5000))
plot.PD.buffer.altitude<-ggplot(data.PD.buffer.1)+geom_boxplot(aes(x=Type,y=altitude,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#FF9900","#FFCC00"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,5000,1250),limits=c(0,5000))

#### 
plot.SR.full.area<-ggplot(data.SR.full.1)+geom_boxplot(aes(x=Type,y=Area_1,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=21))+scale_fill_manual(values=c("#009933","#33CC33"),labels=c("Negative","Positive"))
plot.SR.buffer.area<-ggplot(data.SR.buffer.1)+geom_boxplot(aes(x=Type,y=Area_1,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=21))+scale_fill_manual(values=c("#009933","#33CC33"),labels=c("Negative","Positive"))
plot.PD.full.area<-ggplot(data.PD.full.1)+geom_boxplot(aes(x=Type,y=Area_1,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=21))+scale_fill_manual(values=c("#009933","#33CC33"),labels=c("Negative","Positive"))
plot.PD.buffer.area<-ggplot(data.PD.buffer.1)+geom_boxplot(aes(x=Type,y=Area_1,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=21))+scale_fill_manual(values=c("#009933","#33CC33"),labels=c("Negative","Positive"))

####
plot.SR.full.hp<-ggplot(data.SR.full.1)+geom_boxplot(aes(x=Type,y=HP,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#0066FF","#00CCFF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,25,5),limits=c(0,25))
plot.SR.buffer.hp<-ggplot(data.SR.buffer.1)+geom_boxplot(aes(x=Type,y=HP,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#0066FF","#00CCFF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,25,5),limits=c(0,25))
plot.PD.full.hp<-ggplot(data.PD.full.1)+geom_boxplot(aes(x=Type,y=HP,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#0066FF","#00CCFF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,25,5),limits=c(0,25))
plot.PD.buffer.hp<-ggplot(data.PD.buffer.1)+geom_boxplot(aes(x=Type,y=HP,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#0066FF","#00CCFF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(0,25,5),limits=c(0,25))

#####
plot.SR.full.vegetation<-ggplot(data.SR.full.1)+geom_boxplot(aes(x=Type,y=vegetation_loss,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title  = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#9933FF","#CC99FF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(-1,1,0.4),limits=c(-1,1))
plot.SR.buffer.vegetation<-ggplot(data.SR.buffer.1)+geom_boxplot(aes(x=Type,y=vegetation_loss,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#9933FF","#CC99FF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(-1,1,0.4),limits=c(-1,1))
plot.PD.full.vegetation<-ggplot(data.PD.full.1)+geom_boxplot(aes(x=Type,y=vegetation_loss,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#9933FF","#CC99FF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(-1,1,0.4),limits=c(-1,1))
plot.PD.buffer.vegetation<-ggplot(data.PD.buffer.1)+geom_boxplot(aes(x=Type,y=vegetation_loss,fill=Type),width=0.5,position=position_dodge(0.7),outlier.colour = "grey45")+theme_classic(base_size = 21,base_family = "serif")+theme(axis.title = element_blank(),legend.position = "none",plot.title = element_text(hjust=0.5,size=16))+scale_fill_manual(values=c("#9933FF","#CC99FF"),labels=c("Negative","Positive"))+scale_y_continuous(breaks=seq(-1,1,0.4),limits=c(-1,1))

plot<-plot_grid(plotlist = list(plot.SR.full.range,plot.SR.buffer.range,plot.PD.full.range,plot.PD.buffer.range,plot.SR.full.altitude,plot.SR.buffer.altitude,plot.PD.full.altitude,plot.PD.buffer.altitude,plot.SR.full.area,plot.SR.buffer.area,plot.PD.full.area,plot.PD.buffer.area,plot.SR.full.hp,plot.SR.buffer.hp,plot.PD.full.hp,plot.PD.buffer.hp,plot.SR.full.vegetation,plot.SR.buffer.vegetation,plot.PD.full.vegetation,plot.PD.buffer.vegetation),ncol=4,nrow=5,align="hv")
ggsave("plot.1.tiff",width=21,height=21,dpi=600)


######Conservation gaps#######
setwd("E:/PA/J.appl.eco/Revised/PA_analysis")
data.buffer<-read.csv("E:/PA/J.appl.eco/Revised/future_buffer.csv")
data.buffer.rcp2.6<-data.buffer[which(data.buffer$RCP2.6_CSH<=-0.3),]
sp<-as.character(data.buffer.rcp2.6$species)
threshold<-read.csv("30%threshold.csv")
data_dis_current<-get(load("E:/PA/J.appl.eco/Revised/distribution/proj_current_buffer_sp.RData"))
data_dis_future<-get(load("E:/PA/J.appl.eco/Revised/distribution/future_dis_20km_buffer/proj_2070_rcp2.6_buffer.RData"))
data_dis_future_1<-data_dis_future[,sp]
data.gap<-as.data.frame(matrix(NA,23718,2))
rownames(data.gap)<-rownames(data_dis_future_1)
colnames(data.gap)<-c("SR","Type")
aa<-apply(data_dis_future_1,1,sum)
data.gap[,1]<-aa
pos<-which(rownames(data_dis_future)%in%threshold$Gridcode==F)
pos.1<-which(rownames(data_dis_future)%in%threshold$Gridcode==T)
data.gap[pos,2]<-"non-PAs"
data.gap[pos.1,2]<-"PAs"
write.csv(data.gap,file="Conservation_gap/Buffer/data.gap.SR.rcp2.6.csv")

#####PD#####
data.PD<-read.csv("E:/PA/J.appl.eco/Revised/endanger_analysis/PD/Buffer/PD.rcp26.csv")
data.gap.PD<-as.data.frame(matrix(NA,23718,2))
rownames(data.gap.PD)<-rownames(data_dis_future_1)
colnames(data.gap.PD)<-c("PD","Type")
data.gap.PD[,1]<-data.PD$PD.future
data.gap.PD[pos,2]<-"non-PAs"
data.gap.PD[pos.1,2]<-"PAs"
write.csv(data.gap.PD,file="Conservation_gap/Buffer/data.gap.PD.rcp2.6.csv")

######dat<-data.gap[which(data.gap$Type=="non-PAs"),]
######quantile(dat$SR,0.95)
######quantile(dat$SR,0.85)
#####data<-read.csv("E:/PA/J.appl.eco/Revised/endanger_analysis/PD/Full/residuals.future.rcp6.0.csv")

#####determine the conservation gaps at county level######
setwd("E:/PA/J.appl.eco/Revised/PA_analysis/Conservation_gap/county/5% county")
data<-read.csv("5%_1.csv")
data.gap<-read.csv("data.gap.PD.rcp8.5.csv",row.names=1)
data.gap.1<-data.gap[which(data.gap$CB=="5% Gaps"),]
data.1<-data[which(data$Gridcode %in% row.names(data.gap.1)==T),]
county<-unique(data.1$NAME99)
table.gap<-as.data.frame(matrix(NA,length(unique(data.1$NAME99)),4))
colnames(table.gap)<-c("County_Name","Area","Shape_area","N")
table.gap$County_Name<-county
for(i in 1:length(county)){
  dat<-data.1[which(data.1$NAME99==county[i]),]
  area<-sum(dat$Shape_Area_1)
  table.gap[i,2]<-mean(dat$Area)
  table.gap[i,3]<-area
  table.gap[i,4]<-nrow(dat)
}
write.csv(table.gap,file="table.gap.5%.PD.csv")









