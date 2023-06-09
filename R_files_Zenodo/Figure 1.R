library(snp)
library(nat)
library(Seurat)
library(randomForest)
library(dplyr)
library(rstatix)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(ComplexHeatmap)
library(PRROC)
library(ggpmisc)

source("~/LS/R/Functions/Plot.R")
source("~/LS/R/Functions/Tree.R")

DIRPATH<-"~/LS/data/raw_data/"
SAMPLEID<-list.dirs(DIRPATH,full.names = F,recursive = F)
SAMPLEID<-SAMPLEID[grep("(1|2)",SAMPLEID)]
SAMPLEID<-SAMPLEID[which(!grepl("Cluster",SAMPLEID))]
FIGUREDIR<-"~/课题/Work\\ report/Figures/"

###########################
load("~/LS/data/raw_data/Neurons")
tdata<-Neurons
rm(Neurons)
gc()
##########################################
###classification of segments by several features
library(randomForest)#for RF

###read training set
SAMPLEID<-c("201767","201765","201784","201786","201787")

tests<-subset(tdata,names(tdata) %in% c("201767_001","201767_002","201767_003","201767_006","201767_005","201767_007","201767_008","201767_009","201767_014","201767_015",
                                        "201767_016","201767_018","201767_019","201767_024","201767_026","201767_030","201767_033","201784_015","193869_002","193869_007",
                                        "193869_012","193871_001","193871_002","193871_005","193871_006","193871_008","193872_006","193872_010","193872_015","194530_010",
                                        "193873_002","193873_003","194527_001","194527_006","194527_007","194527_009","194527_010","194527_012","194530_004","194530_006"))


#############
####training
Segment<-c()
for(i in 1:length(tests)){
  print(i)
  segment<-data.frame(NeuronName=names(tests)[i],SegID=1:tests[[i]]$NumSegs,Passing=NA,
                      Len=seglengths(tests[[i]]),Branchorder=NA,Ratio_PSegLen=NA,min_CSegLen=NA,max_CSegLen=NA,Degree=NA,mean_EucDist=NA,max_SegLen=NA)
 
  rid<-unlist(lapply(tests[[i]]$SegList,function(x)x[1]))
  tid<-unlist(lapply(tests[[i]]$SegList,function(x)x[length(x)]))
  segment$Ratio_PSegLen<-segment$Len[match(rid,tid)]/segment$Len
  segment$Ratio_PSegLen[1]<-0
  neuron.tr<-gettree(tests[[i]])

  for(j in 1:tests[[i]]$NumSegs){
    segment$Branchorder[j]<-min(tests[[i]]$d$W[tests[[i]]$SegList[[j]]]/2)
    segment$min_CSegLen[j]<-ifelse(length(which(rid==tid[j])),min(segment$Len[which(rid==tid[j])]),0)
    segment$max_CSegLen[j]<-ifelse(length(which(rid==tid[j])),max(segment$Len[which(rid==tid[j])]),0)
    bid<-neuron.tr$branchID[which(neuron.tr$pointid==tid[j])]
    segment$Degree[j]<-length(which(substr(neuron.tr$branchID,1,nchar(bid))==bid))
    segment$mean_EucDist[j]<-ifelse(segment$Degree[j]-1,mean(dist(tests[[i]]$d[which(tests[[i]]$d$PointNo %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)]),3:5])),0)
    segment$max_SegLen[j]<-ifelse(segment$Degree[j]-1,max(segment$Len[which(rid %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)])]),0)
  }
    
  Segment<-rbind(Segment,segment)
}

##read labels
for(i in 1:length(tests)){
  print(i)
  test<-read.neuron(paste0(DIRPATH,tests[,]$SampleID[i],"/swc_test/",formatC(tests[,]$Index[i],digits = 2,flag = "0",mode = "integer"),".swc"))
  for(j in 1:tests[[i]]$NumSegs){
      Segment$Passing[which(Segment$NeuronName==names(tests)[i]&Segment$SegID==j)]<-any(test$d$Label[test$SegList[[j]]] %in% c(10,15,20))
  }
}

samples1<-sample(which(Segment$Passing==T),ceiling(length(which(Segment$Passing==T))*9/10))
samples2<-sample(which(Segment$Passing==F),ceiling(length(which(Segment$Passing==F))*9/10))
samples1<-sample(which(Segment$Passing==T),ceiling(length(which(Segment$Passing==T))*9/10))
samples2<-sample(which(Segment$Passing==F),ceiling(length(which(Segment$Passing==F))*9/10))
tmp<-Segment[c(samples1,samples2),]

delete=which(colnames(tmp) %in% c("NeuronName","SegID","Passing"))

RF<-randomForest(tmp[,-delete],as.factor(tmp$Passing),importance=TRUE,ntree = 1000,mtry=6,
                 keep.inbag=T,
                 # cutoff = c(2/3,1/3),
                 # strata = as.factor(tmp$Passing),
                 # sampsize = c(20000,500)
                 )
plot(RF)

# predictions=predict(RF,tmp[,-delete],type="prob")
# tmp<-Segment[which(Segment$NeuronName %in% neuronnames),]
tmp<-Segment[c(samples1,samples2),]
predictions=predict(RF,tmp[,-delete])
tmp<-table(predictions,as.factor(tmp$Passing))
print(paste0("F-measure:",(2*tmp[2,2])/(2*tmp[2,2]+tmp[1,2]+tmp[2,1])))
print(paste0("Accuracy:",(tmp[1,1]+tmp[2,2])/(tmp[1,1]+tmp[2,2]+tmp[1,2]+tmp[2,1])))
print(paste0("Precision:",(tmp[2,2])/(tmp[2,2]+tmp[2,1])))
print(paste0("Recall:",(tmp[2,2])/(tmp[2,2]+tmp[1,2])))
print(paste0("MCC:",(tmp[1,1]*tmp[2,2]-tmp[1,2]*tmp[2,1])/(sqrt((tmp[2,2]+tmp[2,1]))*sqrt((tmp[2,2]+tmp[1,2]))*sqrt((tmp[1,1]+tmp[2,1]))*sqrt((tmp[1,1]+tmp[1,2])))))

# tmp<-Segment[which(!Segment$NeuronName %in% neuronnames),]
tmp<-Segment[c(which(!1:nrow(Segment) %in% samples1&Segment$Passing==T),which(!1:nrow(Segment) %in% samples2&Segment$Passing==F)),]
predictions=predict(RF,tmp[,-delete])
tmp<-table(predictions,as.factor(tmp$Passing))
print(paste0("F-measure:",(2*tmp[2,2])/(2*tmp[2,2]+tmp[1,2]+tmp[2,1])))
print(paste0("Accuracy:",(tmp[1,1]+tmp[2,2])/(tmp[1,1]+tmp[2,2]+tmp[1,2]+tmp[2,1])))
print(paste0("Precision:",(tmp[2,2])/(tmp[2,2]+tmp[2,1])))
print(paste0("Recall:",(tmp[2,2])/(tmp[2,2]+tmp[1,2])))
print(paste0("MCC:",(tmp[1,1]*tmp[2,2]-tmp[1,2]*tmp[2,1])/(sqrt((tmp[2,2]+tmp[2,1]))*sqrt((tmp[2,2]+tmp[1,2]))*sqrt((tmp[1,1]+tmp[2,1]))*sqrt((tmp[1,1]+tmp[1,2])))))

imp<-importance(RF,type = 2)
n=1
dotchart(imp[order(imp)],
         labels=rownames(imp),xlab="MeanDecreaseGini")

tmp<-Segment[c(which(!1:nrow(Segment) %in% samples1&Segment$Passing==T),which(!1:nrow(Segment) %in% samples2&Segment$Passing==F)),]
predictions=predict(RF,tmp[,-delete],type="prob")
roc_rf<-PRROC::roc.curve(scores.class0 = predictions[,2],weights.class0 = tmp$Passing,curve = T)
plot(roc_rf,color=F,lwd=1)

pred<-prediction(predictions = predictions[,2],labels = tmp$Passing)
perf<-performance(prediction.obj = pred,measure = "sens",x.measure = "fpr")
plot(perf)
auc<-performance(prediction.obj = pred,"auc")

####################
####partition of terminal arbors
# neuronsdir<-"swc_test/"
allenarbors<-neuronlist()
n<-1
for(i in 1:length(tdata)){
  print(i)
  
  i<-sample(which(tdata[,]$Type=="ITc"),1)
  plot3d(tdata[[i]])
  segment<-data.frame(NeuronName=names(tdata)[i],SegID=1:tdata[[i]]$NumSegs,Passing=NA,Len=seglengths(tdata[[i]]),Branchorder=0,Ratio_PSegLen=0,min_CSegLen=0,max_CSegLen=0,Degree=0,mean_EucDist=0,max_SegLen=0)
  # segment<-data.frame(NeuronName=names(tdata)[i],SegID=1:tdata[[i]]$NumSegs,Passing=NA,mean_Intens=NA,
  #                     var_Intens=NA,mean_Intens_Bouton=0,mean_Size_Bouton=0,Bouton_Dens=0,IBI=0,var_IBI=0,
  #                     Len=seglengths(tdata[[i]]),Branchorder=NA,Ratio_PSegLen=NA,min_CSegLen=NA,max_CSegLen=NA,Degree=NA,mean_EucDist=NA,max_SegLen=NA)
  
  rid<-unlist(lapply(tdata[[i]]$SegList,function(x)x[1]))
  tid<-unlist(lapply(tdata[[i]]$SegList,function(x)x[length(x)]))
  segment$Ratio_PSegLen<-segment$Len[match(rid,tid)]/segment$Len
  segment$Ratio_PSegLen[which(is.na(segment$Ratio_PSegLen))]<-0
  neuron.tr<-gettree(tdata[[i]])
  
  for(j in 1:tdata[[i]]$NumSegs){
    # segment$BoutonDens[j]<-length(which(boutons[[i]]$d$Label==23&boutons[[i]]$d$Parent %in% tdata[[i]]$SegList[[j]]))/segment$Len[j]
    if(length(which(names(tests)==names(tdata)[i]))>0){
      segment$Passing[j]<-any(test$d$Label[tdata[[i]]$SegList[[j]]] %in% c(10,15,20))
    }
    segment$Branchorder[j]<-min(tdata[[i]]$d$W[tdata[[i]]$SegList[[j]]]/2)
    segment$min_CSegLen[j]<-ifelse(length(which(rid==tid[j])),min(segment$Len[which(rid==tid[j])]),0)
    segment$max_CSegLen[j]<-ifelse(length(which(rid==tid[j])),max(segment$Len[which(rid==tid[j])]),0)
    bid<-neuron.tr$branchID[which(neuron.tr$pointid==tid[j])]
    segment$Degree[j]<-length(which(substr(neuron.tr$branchID,1,nchar(bid))==bid))
    segment$mean_EucDist[j]<-ifelse(segment$Degree[j]-1,mean(dist(tdata[[i]]$d[which(tdata[[i]]$d$PointNo %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)]),3:5])),0)
    segment$max_SegLen[j]<-ifelse(segment$Degree[j]-1,max(segment$Len[which(rid %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)])]),0)
  }
    
  delete=which(colnames(segment) %in% c("NeuronName","SegID","Passing"))
  tmp<-segment[,-delete]
  tmp[which(is.na(tmp),arr.ind = T)]<-0
  segment$Passing<-predict(RF,tmp)
  prob<-predict(RF,segment[,-delete],type="prob")
  
  for(j in 1:length(tdata[[i]]$SegList)){
    pointid<-tdata[[i]]$SegList[[j]]
    pointid<-pointid[which(pointid!=1)]
    tmpl=2*length(pointid)
    segmentpos=data.frame(X=rep(0,tmpl),Y=rep(0,tmpl),Z=rep(0,tmpl))
    segmentpos[seq(1,tmpl,2),]<-tdata[[i]]$d[pointid,3:5]
    segmentpos[seq(2,tmpl,2),]<-tdata[[i]]$d[which(tdata[[i]]$d$PointNo %in% tdata[[i]]$d$Parent[pointid]),3:5]
    col_fun<-colorRamp2(c(0,1),c("blue","red"))
    # segments3d(segmentpos,col=colorRampPalette(c("grey",cl.df$cluster_color[which(cl.df$cluster_id==cc)]))(maxW)[i],add=T,lwd=10*i/maxW)
    segments3d(segmentpos,col=col_fun(prob[j,2]),add=T,lwd=2)
    # segments3d(segmentpos,col=ifelse(tmp$Bouton_Dens<0.012&tmp$Len>100,"red","grey")[j],add=T,lwd=1)
  }
  # plot3d(tdata[[i]]$d$X[1],tdata[[i]]$d$Y[1],tdata[[i]]$d$Z[1],col=col_fun(prob[1,2]),add=T,type="s",radius=50)
  # 
  # tmp<-getarbors2(tdata[[i]],Pointlist = tdata[[i]]$d$PointNo[-unique(unlist(tdata[[i]]$SegList[segment$SegID[which(segment$Passing==T)]]))],PlotArbor = T)
  # getarbors(tdata[[i]],ExcludePoints = tdata[[i]]$d$PointNo[which(tdata[[i]]$d$W==2)],SaveArbors = F,PlotArbor=T)
  # tmp<-getarbors(tdata[[i]],ExcludePoints = tdata[[i]]$d$PointNo[unique(unlist(tdata[[i]]$SegList[segment$SegID[which(segment$Passing==T)]]))],PlotArbor = T,SaveArbors = F,CutreeMethod="noCluster")
  arbors<-get_arbors(tdata[[i]],RF,segment[,-delete],plotArbor =F)
  open3d()
  plot3d(arbors,add=T)
  plot3d(tdata[[i]],col="black",WithNodes=F,add=T)
  if(length(arbors)>0){
    names(arbors)<-paste0(tdata[,]$SampleID[i],"_",names(arbors))
    data.frame(arbors)<-summary(arbors)
    arbors[,]$Type<-tdata[,]$Type[i]
    arbors[,]$somaposition<-tdata[,]$BrainRegion[i]
    arbors[,]$somalaminar<-tdata[,]$LaminarPosition[i]
    arbors[,]$NeuronName<-names(tdata)[i]
    arbors<-subset(arbors,cable.length>1000)
    allenarbors<-nat::union(allenarbors,arbors)
  }
  
  if((n %% 100)==0){
    save(allenarbors,file="~/LS/data/allenarbors_0106")
  }
  
  n<-n+1
  
  # #write swc file
  # if(!names(neurons_raw)[i] %in% names(tests)){
  #   d<-neurons_raw[[i]]$d[,1:7]
  #   d$W<-d$W/2
  #   d$Label<-0
  #   d$Label[which(d$PointNo %in% unlist(tdata[[i]]$SegList[which(segment$Passing==T)]))]<-10
  #   d$Label[1]<-1
  #   d$Parent<-match(d$Parent,d$PointNo)
  #   d$PointNo<-1:nrow(d)
  #   d$Parent[which(is.na(d$Parent))]<--1
  #   
  #   if(!dir.exists(paste0(DIRPATH,tdata[,]$SampleID[i],"/",neuronsdir))){
  #     dir.create(paste0(DIRPATH,tdata[,]$SampleID[i],"/",neuronsdir))
  #   }
  #   write.table(matrix(as.numeric(unlist(d)),nrow=nrow(d)),paste0(DIRPATH,tdata[,]$SampleID[i],"/",neuronsdir,formatC(tdata[,]$Index[i],digits = 2,flag = "0",mode = "integer"),".swc"),col.names=F,row.names=F)
  # }
}

load("~/LS/data/allenarbors_0106")
allenarbors<-subset(allenarbors,Type %in% c("ITi","ITc","PT","CT"))
data<-data.frame(Type=factor(allenarbors[,]$Type,levels=c("ITi","ITc","PT","CT")),
                 NeuronName=allenarbors[,]$NeuronName) %>%
  group_by(NeuronName) %>%
  count(Type)
ggplot(data,aes(x=Type,y=n,fill=Type))+
  geom_violin(scale="width",size = 0.2)+
  stat_summary(fun="mean",geom="point",size=0.01)+
  stat_compare_means(aes(label = ..p.signif..),hide.ns = T,symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("*", "")),
                     method = "wilcox.test", ref.group = ".all.")

#################
###statistics
allenarbors[,]$Type<-as.factor(substr(allenarbors[,]$Type,1,2))
allenarbors<-subset(allenarbors,Type %in% c("IT","PT","CT"))
allenarbors[,]$Type<-as.factor(substr(allenarbors[,]$Type,1,2))
tmp<-summary(allenarbors[,]$Type)

ggplot(data.frame(Type=names(tmp),Number=tmp),aes(x=Type,y=Number,fill=Type))+
  geom_bar(stat = "identity", position="dodge",width=0.6)+
  labs(title="Distribution of arbors over neuron types",x="Neuron type",y="Neuron number")+
  theme_classic()+
  theme(plot.title=element_text(size=9,hjust = 0.5),
        text=element_text(size=8),
        axis.text.x = element_text(size=7),axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=8),axis.title.y = element_text(size=8),
        legend.position="none",
        legend.title = element_text(size=8),legend.text = element_text(size=8))
ggsave(paste0(FIGUREDIR,"Distribution of arbors over neuron types.pdf"),device = "pdf",width = 60,height = 40,scale=1,units = "mm")

data<-data.frame(Type=allenarbors_dataframe$AxonSubtype,
                 NeuronName=allenarbors_dataframe$NeuronName)
data$Type[which(data$Type %in% 1:44)]<-"IT"
data$Type[which(data$Type %in% 45:52)]<-"CT"
data$Type[which(data$Type %in% 53:64)]<-"PT"
data<-subset(data,!is.na(data$Type))
data<-data %>%
  group_by(NeuronName) %>%
  count(Type)

stat.test <- data[,c("n","Type")] %>%
  # group_by(Type) %>%
  wilcox_test(n ~ Type) %>%
  adjust_pvalue(method = "fdr") %>%
  rstatix::filter(p.adj.signif != "ns") %>%
  add_y_position()
ggplot(data,aes(x=Type,y=n,fill=Type))+
  geom_violin(scale="width",size = 0.2)+
  stat_summary(fun="mean",geom="point",size=0.01)+
  stat_compare_means(aes(label = ..p.signif..),hide.ns = T,
                     # symnum.args=list(cutpoints = c(0, 0.05, 1), symbols = c("*", "")),
                     method = "wilcox.test", ref.group = ".all.",comparisons = list(c(1,2),c(1,3),c(2,3)))+
  labs(title="Average number of arbors per neuron\\n over different types",x="Neuron type",y="Neuron number")+
  scale_y_continuous(n.breaks = 4)+
  theme_classic()+
  theme(plot.title=element_text(size=9,hjust = 0.5),
        text=element_text(size=8),
        axis.text.x = element_text(size=7),axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=8),axis.title.y = element_text(size=8),
        legend.position="none",
        legend.title = element_text(size=8),legend.text = element_text(size=8))
ggsave(paste0(FIGUREDIR,"Average number of arbors per neuron over different types.pdf"),device = "pdf",width = 60,height = 40,scale=1,units = "mm")

##############
###validation of boutons
load("~/LS/data/allenarbors_0106")

neuronnames<-c("201767_001","201767_002","201767_003","201767_006","201767_005","201767_007","201767_008","201767_009","201767_014","201767_015",
               "201767_016","201767_018","201767_019","201767_024","201767_026","201767_030","201767_033",
               # "201784_015",
               "193869_002","193869_007",
               "193869_012","193871_001","193871_002","193871_005","193871_006","193871_008","193872_006","193872_010","193872_015","194530_010",
               "193873_002","193873_003","194527_001","194527_006","194527_007","194527_009","194527_010","194527_012","194530_004","194530_006",
               "201767_023","201767_025","201767_032","201784_001","201784_002")

# neuronnames<-c("201784_002")
# tests<-subset(Neurons,names(Neurons) %in% neuronnames)
tests<-neuronlist()
for(i in 1:length(neuronnames)){
  print(i)
  test<-read.neurons(paste0("~/LS/data/raw_data/",strsplit(neuronnames[i],split = "_")[[1]][1],"/swc_raw/",strsplit(neuronnames[i],split = "_")[[1]][2],".swc"))
  names(test)<-neuronnames[i]
  test[,]$SampleID<-strsplit(neuronnames[i],split = "_")[[1]][1]
  test[,]$Index<-strsplit(neuronnames[i],split = "_")[[1]][2]
  tests<-nat::union(tests,test)
}

intensity<-neuronlist()
for(i in 1:length(neuronnames)){
  print(i)
  test<-read.neurons(paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um/",strsplit(neuronnames[i],split = "_")[[1]][1],"/",strsplit(neuronnames[i],split = "_")[[1]][2],"_normalize_200.swc"))
  names(test)<-neuronnames[i]
  test[,]$SampleID<-strsplit(neuronnames[i],split = "_")[[1]][1]
  test[,]$Index<-strsplit(neuronnames[i],split = "_")[[1]][2]
  intensity<-nat::union(intensity,test)
}

Boutons<-c()
for(i in 1:length(tests)){
  print(i)
  boutons<-read.csv(paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um_fit_result/",tests[,]$SampleID[i],"/",formatC(tests[,]$Index[i],digits = 2,flag = "0",mode = "integer"),"_fit_LoG.txt"),header = F)
  colnames(boutons)<-c("SegID","NodeID","mu","sigma","A")
  boutons$Passing<-NA
  boutons$IBI<-NA
  boutons$Intensity_r<-NA
  boutons$Order<-NA
  d<-intensity[[which(names(intensity)==names(tests)[i])]]$d
  boutons$Intensity<-d$W[boutons$NodeID]
  boutons$NeuronName<-names(tests)[[i]]
  
  # test<-read.neuron(paste0("~/LS/data/raw_data/",tests[,]$SampleID[i],"/swc_test/",formatC(tests[,]$Index[i],digits = 2,flag = "0",mode = "integer"),".swc"))
  test<-read.neuron(paste0("~/LS/data/raw_data/",tests[,]$SampleID[i],"/swc_allen_space/",formatC(tests[,]$Index[i],digits = 2,flag = "0",mode = "integer"),".swc"))
  test$d$Label<-0
  arbors<-subset(allenarbors,NeuronName==names(tests)[i])
  boutons$Type<-unique(arbors[,]$Type)
  test$d$Label<-10
  ad<-cbind(unlist(lapply(arbors,function(x)x$d$X)),unlist(lapply(arbors,function(x)x$d$Y)),unlist(lapply(arbors,function(x)x$d$Z)))
  ind<-which(test$d$X %in% ad[,1]&test$d$Y %in% ad[,2]&test$d$Z %in% ad[,3])
  test$d$Label[ind]<-0
  test$d$Label[1]<-1
  
  seglen<-seglengths(test)
  
  for(j in 1:tests[[i]]$NumSegs){
    boutons$Passing[which(boutons$SegID==j)]<-any(test$d$Label[test$SegList[[j]]] %in% c(10,15,20))
    boutons$Intensity_r[which(boutons$SegID==j)]<-d$W[boutons$NodeID[which(boutons$SegID==j)]]/mean(d$W[test$SegList[[j]]])
    neuron<-Neurons[[which(names(Neurons)==names(tests)[i])]]
    boutons$Order[which(boutons$SegID==j)]<-max(neuron$d$W[neuron$SegList[[j]]])
    mu<-boutons$mu[which(boutons$SegID==j)]
    
    if(length(mu)>2){
      boutons$IBI[which(boutons$SegID==j)]<-mu[c(2:length(mu),length(mu))]-mu[c(1:(length(mu)-1),length(mu)-1)]
      if(length(which(boutons$IBI[which(boutons$SegID==j)]<=0))>0&&which(boutons$IBI[which(boutons$SegID==j)]<=0)!=1){
        boutons$IBI[which(boutons$IBI<=0&boutons$SegID==j)]<-mu[which(boutons$IBI[which(boutons$SegID==j)]<=0)]-mu[which(boutons$IBI[which(boutons$SegID==j)]<=0)-1]
      }
    }else{
      boutons$IBI[which(boutons$SegID==j)]<-seglen[j]
    }
  }
  
  Boutons<-rbind(Boutons,boutons)
}
Boutons<-subset(Boutons,IBI>0&IBI<10)
Boutons[,]$Type<-as.factor(substr(Boutons[,]$Type,1,2))

ggviolin(Boutons,x="Passing",y="IBI",fill="Passing",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2)))
ggviolin(Boutons[which(Boutons$Passing==F),],x="Type",y="IBI",fill="Type",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2),c(1,3),c(2,3)))
ggviolin(Boutons[which(Boutons$Passing==F),],x="Order",y="IBI",fill="Order",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)
# statistics<-Boutons[,c("sigma","A","IBI","Intensity")]
# sampleid<-c(sample(which(Boutons$Passing==T),1000),sample(which(Boutons$Passing==F),1000))
sampleid<-which(Boutons$NeuronName %in% c("201784_002"))
# sampleid<-1:nrow(Boutons)
statistics<-Boutons[sampleid,c("sigma","A","Intensity_r")]

for(i in 1:ncol(statistics)){
  statistics[,i]<-(statistics[,i]-mean(statistics[,i]))/sd(statistics[,i])
}

res.pca<-PCA(statistics)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 100)) 
fviz_pca_ind(res.pca, 
             # Individuals
             geom.ind = "point",
             fill.ind = factor(Boutons$Passing), col.ind = factor(Boutons$Passing),
             pointshape = 21, pointsize = 2,
             palette =c("red","green"),
             addEllipses = TRUE,
             # Variables
             # col.var = "contrib",repel = T,
             legend.title = list(fill = "Type")
)

plot3d(statistics[,1],statistics[,2],statistics[,3],col=getcol(Boutons$Passing))

ggplot(Boutons, aes(x =IBI,fill=Passing))+geom_density(alpha=0.5)
ggplot(Boutons, aes(x =sigma,fill=Passing))+geom_density(alpha=0.5)
ggplot(Boutons, aes(x =A,fill=Passing))+geom_density(alpha=0.5)+xlim(c(0,10))
ggviolin(Boutons[sampleid,],x="Passing",y="A",fill="Passing",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2)))+
  ylim(c(0,10))
ggviolin(Boutons[sampleid,],x="Passing",y="sigma",fill="Passing",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2)))

#classification
annotation_col=data.frame(Passing=as.factor(Boutons$Passing[sampleid]))
rownames(annotation_col)=rownames(statistics)
ht<-Heatmap(t(as.matrix(statistics)),
            # col=col_fun,
            show_column_names = FALSE,
            # clustering_method_rows = "ward.D",
            # clustering_distance_rows = "euclidean",
            clustering_method_columns = "ward.D2",
            # cluster_columns = tmp,
            # column_split=3,
            # column_split=statistics$m,
            heatmap_legend_param = list(title="Z-score"),
            top_annotation = columnAnnotation(df=as.data.frame(annotation_col)),
            use_raster=T)
ht<-draw(ht)
gc()
statistics$NeuronName<-Boutons$NeuronName[sampleid]
statistics$NodeID<-Boutons$NodeID[sampleid]
statistics$passing<-Boutons$Passing[sampleid]
statistics$SegID<-Boutons$SegID[sampleid]
tmp<-column_dend(ht)
m<-cutree(as.hclust(tmp),k=3)
statistics$m<-m
# statistics$m[which(statistics$m!=1)]<-2

bfreq<-c()
Cluster<-c()
type<-c()
P<-c()
for(i in unique(statistics$passing)){
  for(j in unique(statistics$m)){
    bfreq<-c(bfreq,length(which(statistics$passing==i&statistics$m==j))/length(which(statistics$passing==i)))
    Cluster<-c(Cluster,i)
    type<-c(type,j)
    ma<-matrix(nrow=2,ncol=2)
    ma[1,1]<-length(which(statistics$passing==i&statistics$m==j))
    ma[2,1]<-length(which(statistics$passing!=i&statistics$m==j))
    ma[1,2]<-length(which(statistics$passing==i&statistics$m!=j))
    ma[2,2]<-length(which(statistics$passing!=i&statistics$m!=j))
    P<-c(P,fisher.test(ma,alternative = "greater")$p.value)
  }
}
marker<-rep("",length(P))
marker[which(P<=0.05&P>0)]<-"*"
marker[which(P<=0.01&P>0)]<-"**"
marker[which(P<=0.001&P>0)]<-"***"
marker[which(P<=0.0001&P>0)]<-"****"
ggplot(data.frame(Proportion=bfreq,Cluster=as.factor(Cluster),Type=as.factor(type)), aes(x =Cluster , y = Proportion,fill=Type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8),width = 0.8)+
  geom_text(aes(x=Cluster,y=Proportion,label=marker),size=3,position = position_dodge(0.8),stat = "identity")+
  theme_classic()+
  # labs(x="Clusters",title="Enrichment of PT, CT and IT arbors\\n in different clusters")+
  theme(axis.text.x = element_text(vjust = 0.5,size=7),axis.text.y =element_text(size=7),
        axis.title.x =element_text(size=8),axis.title.y =element_text(size=8),
        plot.title=element_text(hjust = 0.5,size=8),
        # plot.title=element_blank(),
        legend.title = element_text(hjust = 0.5,size=8),legend.text = element_text(hjust = 0.5,size=8),
        legend.key.height=unit(2.5,"mm"),legend.key.width = unit(2,"mm"),legend.text.align=0,
        plot.tag.position = "right",plot.tag = element_text(angle=270,size=9))+
  # scale_fill_manual(values=c("black","grey70"))+
  scale_y_continuous(n.breaks = 2)

m<-1
Segment<-c()
for(n in unique(statistics$NeuronName)){
  neuron<-intensity[[which(names(intensity)==n)]]
  seglen<-seglengths(neuron)
  segs<-data.frame(NeuronName=n,SegID=1:neuron$NumSegs,passing=NA,density=NA,Len=seglen)
  
  test<-read.neuron(paste0("~/LS/data/raw_data/",strsplit(n,split = "_")[[1]][1],"/swc_allen_space/",formatC(strsplit(n,split = "_")[[1]][2],digits = 2,flag = "0",mode = "integer"),".swc"))
  test$d$Label<-0
  arbors<-subset(allenarbors,NeuronName==n)
  test$d$Label<-10
  ad<-cbind(unlist(lapply(arbors,function(x)x$d$X)),unlist(lapply(arbors,function(x)x$d$Y)),unlist(lapply(arbors,function(x)x$d$Z)))
  ind<-which(test$d$X %in% ad[,1]&test$d$Y %in% ad[,2]&test$d$Z %in% ad[,3])
  test$d$Label[ind]<-0
  test$d$Label[1]<-1
  
  for(j in 1:neuron$NumSegs){
    segs$density[j]<-length(which(statistics$NeuronName==n&statistics$SegID==j&statistics$m==m))/seglen[j]
    segs$passing[j]<-any(test$d$Label[test$SegList[[j]]] %in% c(10,15,20))
  }
  Segment<-rbind(Segment,segs)
}

# ggplot(Segment,aes(x=density,fill=passing))+geom_density(alpha=0.5)
ggviolin(Segment,x="passing",y="density",fill="passing",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2)))

#write 
options(scipen=10)
for(n in unique(statistics$NeuronName)){
  neuron<-intensity[[which(names(intensity)==n)]]
  neuron$d$Label<-0
  neuron$d$Label[1]<-1
  if(!dir.exists(paste0("~/LS/data/raw_data/",strsplit(n,split = "_")[[1]][1],"/swc_boutons/"))){
    dir.create(paste0("~/LS/data/raw_data/",strsplit(n,split = "_")[[1]][1],"/swc_boutons/"))
  }
  for(m in unique(statistics$m[which(statistics$NeuronName==n)])){
    d<-neuron$d
    d$Label[statistics$NodeID[which(statistics$NeuronName==n&statistics$m==m)]]<-10
    d$Label[1]<-1
    d$W<-d$W/2
    write.table(matrix(as.numeric(unlist(d)),nrow=nrow(d)),paste0("~/LS/data/raw_data/",strsplit(n,split = "_")[[1]][1],"/swc_boutons/",strsplit(n,split = "_")[[1]][2],"_",m,".swc"),col.names=F,row.names=F)
  }
}

arbors<-subset(allenarbors,NeuronName %in% neuronnames)
arbors[,]$BoutonNum<-NA
for(i in 1:length(arbors)){
  print(i)
  test<-read.neuron(paste0("~/LS/data/raw_data/",strsplit(arbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_allen_space/",strsplit(arbors[,]$NeuronName[i],split = "_")[[1]][2],".swc"))
  boutons<-read.csv(paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um_fit_result/",strsplit(arbors[,]$NeuronName[i],split = "_")[[1]][1],"/",strsplit(arbors[,]$NeuronName[i],split = "_")[[1]][2],"_fit_LoG.txt"),header = F)
  colnames(boutons)<-c("SegID","NodeID","mu","sigma","A")
  test$d$Label<-0
  test$d$Label[1]<-1
  ad<-arbors[[i]]$d[,3:5]
  ind<-which(test$d$X %in% ad[,1]&test$d$Y %in% ad[,2]&test$d$Z %in% ad[,3])
  test$d$Label[ind]<-10
  segid<-unlist(lapply(test$SegList,function(x)any(test$d$Label[x] %in% c(10,15,20))))
  segid<-which(segid==T)
  arbors[,]$BoutonNum[i]<-length(which(boutons$SegID %in% segid))
}

##inter-bouton interval
allenarbors_dataframe<-subset(allenarbors_dataframe,!is.na(cable.length)&!is.na(BoutonNum))
allenarbors_dataframe$Type<-substr(as.character(allenarbors_dataframe$Type),1,2)
allenarbors_dataframe<-subset(allenarbors_dataframe,Type %in% c("IT","PT","CT"))
ggplot(data.frame(yy=allenarbors_dataframe$cable.length,xx=allenarbors_dataframe$BoutonNum,Type=factor(substr(allenarbors_dataframe$Type,1,2),levels=c("IT","PT","CT"))),
       aes(x = xx, y = yy,fill=Type,color=Type)) +   
  geom_point(size=0.1,pch=20) +
  geom_smooth(method = 'lm', formula = y ~ x,lwd=0.5) +
  scale_color_manual(values = c("red","green","blue"))+
  # geom_text(aes(label=n),cex=1.5,color=colors)+
  stat_poly_eq(formula=y~x,aes(label=paste(..eq.label..,..rr.label..,sep="~~~")),parse=T)+
  # geom_tile(aes(x = xx, y = yy))+
  labs(x="Bouton Number",y="Length")+
  theme_classic()+
  theme(line = element_line(size=0.1),plot.title = element_blank(),
        axis.text.y =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.text.x =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.title.y =element_text(size=6),axis.title.x =element_text(size=6),
        strip.text.y = element_text(angle=0,size=6),
        strip.placement = "outside",
        panel.spacing = unit(0,"mm"),
        legend.background = element_rect(colour="black",size=0.1),
        legend.title = element_text(hjust = 0.5,size=6),legend.text = element_text(hjust = 0.5,size=6)
  )
ggsave(paste0(FIGUREDIR,"The relationship between length and bouton number.pdf"),device = "pdf",width = 60,height = 60,scale = 1,units = "mm")

allenarbors_dataframe<-subset(allenarbors_dataframe,!is.na(cable.length)&!is.na(BoutonNum)&R.4 %in% c(14,1,4,25,10,2,20,21,28,9,12,13,3,19,16,18,7,6,15,22,5,24,8,11,17))
ggplot(data.frame(yy=allenarbors_dataframe$cable.length,xx=allenarbors_dataframe$BoutonNum,Cluster=factor(allenarbors_dataframe$R.4,c(14,1,4,25,10,2,20,21,28,9,12,13,3,19,16,18,7,6,15,22,5,24,8,11,17))),
       aes(x = xx, y = yy,fill=Cluster,color=Cluster)) +   
  geom_point(size=0.1,pch=20) +
  facet_wrap(.~Cluster)+
  geom_smooth(method = 'lm', formula = y ~ x,lwd=0.5) +
  # scale_color_manual(values = c("red","green","blue"))+
  # geom_text(aes(label=n),cex=1.5,color=colors)+
  stat_poly_eq(formula=y~x,aes(label=paste(..eq.label..,..rr.label..,sep="~~~")),parse=T)+
  # geom_tile(aes(x = xx, y = yy))+
  labs(x="Bouton Number",y="Length")+
  theme_classic()+
  theme(line = element_line(size=0.1),plot.title = element_blank(),
        axis.text.y =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.text.x =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.title.y =element_text(size=6),axis.title.x =element_text(size=6),
        strip.text.y = element_text(angle=0,size=6),
        strip.placement = "outside",
        panel.spacing = unit(0,"mm"),
        legend.background = element_rect(colour="black",size=0.1),
        legend.title = element_text(hjust = 0.5,size=6),legend.text = element_text(hjust = 0.5,size=6)
  )
ggsave(paste0(FIGUREDIR,"The relationship between length and bouton number.pdf"),device = "pdf",width = 60,height = 60,scale = 1,units = "mm")

data<-subset(allenarbors_dataframe,annot %in% calannot("STR")&Type %in% c("ITi","ITc","PT")&BoutonNum!=0)
data$Type<-substr(data$Type,1,2)
data<-data.frame(IBI=data$cable.length/data$BoutonNum,Type=data$Type)
data<-subset(data,IBI<10)
stat.test <- data %>%
  wilcox_test(IBI ~ Type) %>%
  adjust_pvalue() %>%
  add_xy_position()
ggviolin(data,x="Type",y="IBI",fill="Type",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_pvalue_manual(stat.test, label = "p.adj")+
  scale_y_continuous(n.breaks = 3)+
  labs(x="Type",y="IBI")+
  # scale_fill_manual(values=c("sienna","green","magenta"),name="Class")+
  theme(line = element_line(size=0.001),
        axis.text.y =element_text(size=7),axis.text.x =element_text(size=7),
        axis.title.y =element_text(size=8),axis.title.x =element_text(size=8),
        axis.line = element_line(size=0.1),
        strip.background = element_blank(),strip.text.y = element_text(angle=0,size=7),
        strip.placement = "outside",
        panel.spacing = unit(0,"mm"),
        legend.position = "none")

options(scipen=10)
for(i in which(names(allenarbors) %in% data$arborname[which(data$SampleID=="195041")])){
  print(i)
  test<-read.neuron(paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_allen_space/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][2],".swc"))
  intensity<-read.neuron(paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][2],"_normalize_200.swc"))
  intensity$d$Label<-0
  intensity$d$Label[1]<-1
  tmp<-paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um_fit_result/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][2],"_fit_LoG.txt")
  if(file.exists(tmp)){
    boutons<-read.csv(tmp,header = F)
    colnames(boutons)<-c("SegID","NodeID","mu","sigma","A")
    test$d$Label<-0
    test$d$Label[1]<-1
    ad<-allenarbors[[i]]$d[,3:5]
    ind<-which(test$d$X %in% ad[,1]&test$d$Y %in% ad[,2]&test$d$Z %in% ad[,3])
    test$d$Label[ind]<-10
    segid<-unlist(lapply(test$SegList,function(x)any(test$d$Label[x] %in% c(10,15,20))))
    segid<-which(segid==T)
    pointlist<-unique(unlist(intensity$SegList[segid]))
    d<-intensity$d
    d$Label[boutons$NodeID[which(boutons$SegID %in% segid)]]<-10
    d<-d[pointlist,]
    d$Label[1]<-1
    d$Parent[1]<--1
    d$W<-d$W/2
    if(!dir.exists(paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/"))){
      dir.create(paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/"))
    }
    write.table(matrix(as.numeric(unlist(d)),nrow=nrow(d)),paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/",names(allenarbors)[i],".swc"),col.names=F,row.names=F)
  }
}

#############################
###the morphological difference between passing axons and terminal arbors
load("~/LS/data/allenarbors_0106")
load("~/LS/data/raw_data/neurons_PFC")
Segment<-c()
# PassingAxons<-neuronlist()
for(i in 1:length(neurons)){
  print(i)
  arbors<-subset(allenarbors,NeuronName %in% names(neurons)[i])
  if(length(arbors)>0){
    neuron.tr<-gettree(neurons[[i]])
    
    rid<-unlist(lapply(neurons[[i]]$SegList,function(x)x[1]))
    tid<-unlist(lapply(neurons[[i]]$SegList,function(x)x[length(x)]))
    
    passingnode<-which(!neurons[[i]]$d$X %in% unlist(lapply(arbors,function(x)x$d$X[-1])))
    passingseg<-which(rid %in% passingnode&tid %in% passingnode)
    if(length(passingseg)>0){
      # d<-neurons[[i]]$d[unique(unlist(neurons[[i]]$SegList[passingseg])),1:7]
      # d$Parent <- match(d$Parent, d$PointNo)
      # d$PointNo <- seq(1, nrow(d))
      # d$Label[which(d$PointNo == 1)] <- 1
      # d$Parent[which(is.na(d$Parent))] <- -1
      
      # passingaxons<-as.neuronlist(as.neuron(d))
      # names(passingaxons)<-names(neurons)[[i]]
      # PassingAxons<-nat::union(PassingAxons,passingaxons)
      
      # arborseg<-which(rid %in% c(rid[passingseg],tid[passingseg])&!tid %in% c(rid[passingseg],tid[passingseg]))
      arborseg<-which(!(rid %in% passingnode&tid %in% passingnode))
      
      seglen<-seglengths(neurons[[i]])
      segment<-data.frame(NeuronName=names(neurons)[i],
                          SegID=c(arborseg,passingseg),
                          Passing=c(rep("terminal arbors",length(arborseg)),rep("passing axons",length(passingseg))),
                          Len=seglen[c(arborseg,passingseg)],
                          Branchorder=0,Ratio_PSegLen=0,min_CSegLen=0,max_CSegLen=0,Degree=0,mean_EucDist=0,max_SegLen=0)
      segment$Ratio_PSegLen<-segment$Len[match(rid[segment$SegID],tid[segment$SegID])]/segment$Len
      segment$Ratio_PSegLen[which(is.na(segment$Ratio_PSegLen))]<-0
      
      for(j in segment$SegID){
        segment$Branchorder[which(segment$SegID==j)]<-min(neurons[[i]]$d$W[neurons[[i]]$SegList[[j]]])
        segment$min_CSegLen[which(segment$SegID==j)]<-ifelse(length(which(rid==tid[j])),min(seglen[which(rid==tid[j])]),0)
        segment$max_CSegLen[which(segment$SegID==j)]<-ifelse(length(which(rid==tid[j])),max(seglen[which(rid==tid[j])]),0)
        bid<-neuron.tr$branchID[which(neuron.tr$pointid==tid[j])]
        segment$Degree[which(segment$SegID==j)]<-length(which(substr(neuron.tr$branchID,1,nchar(bid))==bid))
        segment$mean_EucDist[which(segment$SegID==j)]<-ifelse(segment$Degree[which(segment$SegID==j)]-1,mean(dist(neurons[[i]]$d[which(neurons[[i]]$d$PointNo %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)]),3:5])),0)
        segment$max_SegLen[which(segment$SegID==j)]<-ifelse(segment$Degree[which(segment$SegID==j)]-1,max(seglen[which(rid %in% neuron.tr$pointid[which(substr(neuron.tr$branchID,1,nchar(bid))==bid)])]),0)
      }
      Segment<-rbind(Segment,segment)
    }
  }
}
save(Segment,file="~/LS/data/Segment")
save(Segment,file="~/LS/data/Segment2")
save(PassingAxons,file="~/LS/data/PassingAxons")

###############
###statistics of morphological features
Segment$Branchorder[which(Segment$Branchorder==0)]<-1
freq<-c()
Passing<-c()
branchorder<-c()
P<-c()
for(i in min(Segment$Branchorder):12){
  for(j in unique(Segment$Passing)){
    freq<-c(freq,length(which(Segment$Branchorder==i&Segment$Passing==j))/length(which(Segment$Passing==j)))
    branchorder<-c(branchorder,i)
    Passing<-c(Passing,j)
    ma<-matrix(nrow=2,ncol=2)
    ma[1,1]<-length(which(Segment$Branchorder==i&Segment$Passing==j))
    ma[2,1]<-length(which(Segment$Branchorder==i&Segment$Passing!=j))
    ma[1,2]<-length(which(Segment$Branchorder!=i&Segment$Passing==j))
    ma[2,2]<-length(which(Segment$Branchorder!=i&Segment$Passing!=j))
    P<-c(P,fisher.test(ma,alternative = "greater")$p.value)
  }
}
marker<-rep("",length(P))
marker[which(P<=0.05&P>=0)]<-"*"
marker[which(P<=0.01&P>=0)]<-"**"
marker[which(P<=0.001&P>=0)]<-"***"
marker[which(P<=0.0001&P>=0)]<-"****"
data<-data.frame(Proportion=freq,BranchOrder=factor(branchorder),Passing=factor(Passing,levels=c("terminal arbors","passing axons")))
ggplot(data, aes(x =BranchOrder , y = Proportion,fill=Passing)) +
  geom_bar(stat = "identity", position = position_dodge(0.7),width = 0.7)+
  geom_text(aes(x=BranchOrder,y=Proportion,label=marker),size=2,position = position_dodge(0.7),stat = "identity")+
  theme_classic()+
  labs(x="Branch order")+
  theme(axis.text.x = element_text(vjust = 0.5,size=7),axis.text.y =element_text(size=7),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        plot.title=element_text(hjust = 0.5,size=7),
        # plot.title=element_blank(),
        legend.position = "none",
        legend.title = element_blank(),legend.text = element_text(hjust = 0.5,size=5),
        legend.key.height=unit(2.5,"mm"),legend.key.width = unit(2,"mm"),legend.text.align=0,
        plot.tag.position = "right",plot.tag = element_text(angle=270,size=9))+
  scale_y_continuous(n.breaks = 2)+
  # scale_fill_manual(values=getcol(1:k))
  scale_fill_manual(values=c("green","red"))
ggsave(paste0(FIGUREDIR,"The difference of branch order of passing axons and terminal arbors.pdf"),device = "pdf",width = 60,height = 40,units = "mm")

data<-Segment[,c(3,4,7:11)]
data<-subset(data,Len<1000&min_CSegLen<500&max_CSegLen<2000&Degree<1000)
data<-reshape2::melt(data,id.vars = c("Passing"))
data$Passing<-factor(data$Passing,levels=c("terminal arbors","passing axons"))
ggplot(data,aes(x=Passing,y=value,fill=Passing))+
  geom_violin(trim=T,size = 0.1)+
  stat_summary(fun="mean",geom="point",size=0.5)+
  facet_wrap(~variable,scale="free",ncol=length(unique(data$variable)))+
  stat_compare_means(method="wilcox.test",label="p.format",
                     size = 1,
                     comparisons = list(c("terminal arbors","passing axons")))+
  theme_classic()+
  theme(axis.text.x = element_blank(),axis.text.y =element_text(size=7),
        axis.title.x =element_blank(),axis.title.y =element_text(size=7),
        plot.title=element_text(hjust = 0.5,size=7),
        # plot.title=element_blank(),
        legend.position = "none",
        legend.title = element_blank(),legend.text = element_text(hjust = 0.5,size=5),
        legend.key.height=unit(2.5,"mm"),legend.key.width = unit(2,"mm"),legend.text.align=0,
        strip.text = element_text(size=7),
        strip.background = element_blank(),
        plot.tag.position = "right",plot.tag = element_text(angle=270,size=9))+
  scale_y_continuous(n.breaks = 3)+
  # scale_fill_manual(values=getcol(1:k))
  scale_fill_manual(values=c("green","red"))
ggsave(paste0(FIGUREDIR,"The morphological difference of passing axons and terminal arbors.pdf"),device = "pdf",width = 150,height = 40,units = "mm")

###the distribution of length in cortex
bins<-seq(0,1,length.out=31)
len_matx<-matrix(0,nrow=length(bins)-1,ncol=length(PassingAxons))
colnames(len_matx)<-names(PassingAxons)
for(i in 1:(length(bins)-1)){
  print(i)
  for(j in 1:length(PassingAxons)){
    neuron<-neurons[[which(names(neurons)==names(PassingAxons)[j])]]
    pl<-which(neuron$d$anno_W>bins[i]&neuron$d$anno_W<=bins[i+1]&neuron$d$Parent!=-1)
    if(length(pl)>0){
      len_matx[i,j]<-sum(sqrt(rowSums((neuron$d[pl,3:5]-neuron$d[match(neuron$d$Parent[pl],neuron$d$PointNo),3:5])^2)))
    }
  }
}

arbor.len_matx<-matrix(0,nrow=length(bins)-1,ncol=length(PassingAxons))
colnames(arbor.len_matx)<-names(PassingAxons)
for(i in 1:(length(bins)-1)){
  print(i)
  for(j in 1:length(PassingAxons)){
    arbors<-subset(allenarbors,NeuronName==names(PassingAxons)[j])
    neuron<-neurons[[which(names(neurons)==names(PassingAxons)[j])]]
    arborX<-unlist(lapply(arbors,function(x)x$d$X))
    pl<-which(neuron$d$anno_W>bins[i]&neuron$d$anno_W<=bins[i+1]&neuron$d$Parent!=-1&neuron$d$X %in% arborX)
    if(length(pl)>0){
      arbor.len_matx[i,j]<-sum(sqrt(rowSums((neuron$d[pl,3:5]-neuron$d[match(neuron$d$Parent[pl],neuron$d$PointNo),3:5])^2)))
    }
  }
}
save(len_matx,arbor.len_matx,file="~/LS/data/PassingAxons&Arbors_cortical")

for(i in 1:ncol(len_matx)){
  len_matx[,i]<-len_matx[,i]/sum(as.numeric(len_matx[,i]))
}

for(i in 1:ncol(arbor.len_matx)){
  arbor.len_matx[,i]<-arbor.len_matx[,i]/sum(as.numeric(arbor.len_matx[,i]))
}

tmp<-data.frame(Length=c(as.vector(len_matx[1:(length(bins)-1),]),as.vector(arbor.len_matx[1:(length(bins)-1),])),
                Class=c(rep("Passing axons",nrow(len_matx)*ncol(len_matx)),rep("Terminal arbors",nrow(len_matx)*ncol(len_matx))),
                Depth=as.numeric(rep((bins[1:(length(bins)-1)]+bins[2:length(bins)])/2,ncol(len_matx)*2)))
tmp %>% 
  ggplot(aes(x = Depth, y = Length,color=Class,fill=Class)) +
  stat_summary(fun = "mean",
               geom = "ribbon",
               alpha=0.3,
               # geom="pointrange",
               size=0.1,
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  stat_summary(fun="mean",
               geom="line",lwd=0.5)+
  # annotate("text",x=(bins[1:(length(bins)-1)]+bins[2:length(bins)])/2,
  # y=meanlen+0.05,
  # label=marker)+
  theme_classic()+
  geom_vline(
    # xintercept = c(0.1,0.2,0.5),
    xintercept = c(0,0.1,0.25,0.6,1),
    # xintercept = c(0,0.15,0.4,0.7,1),
    lwd=0.5,
    lty=3)+
  # labs(y="Mean length (mm)")+
  labs(y="Relative strength")+
  xlim(c(max(tmp$Depth),min(tmp$Depth)))+
  scale_y_continuous(n.breaks = 3)+
  scale_fill_manual(values=c("red","green"))+
  scale_color_manual(values=c("red","green"))+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5,size=8),
        text=element_text(size=5),
        axis.text.x = element_text(vjust = 0.5,size=5),axis.text.y =element_text(size=5),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        strip.background = element_blank(),strip.text.x = element_text(angle=0,size=7),
        legend.title = element_text(hjust = 0.5,size=7),legend.text = element_text(hjust = 0.5,size=5),
        legend.key.size = unit(1,"mm"),
        legend.position = "none",
        legend.text.align=0)
ggsave(paste0(FIGUREDIR,"The laminar distribution of passing axons and terminal arbors.pdf"),device = "pdf",width =30,height = 40,units = "mm")

##########
neuronnames<-c(sample(intersect(allenarbors_dataframe$NeuronName[which(allenarbors_dataframe$Type=="CT")],names(PassingAxons)[which(PassingAxons[,]$nTrees==1)]),5),
               sample(intersect(allenarbors_dataframe$NeuronName[which(allenarbors_dataframe$Type %in% c("ITc"))],names(PassingAxons)[which(PassingAxons[,]$nTrees==1)]),5),
               sample(intersect(allenarbors_dataframe$NeuronName[which(allenarbors_dataframe$Type=="PT")],names(PassingAxons)[which(PassingAxons[,]$nTrees==1)]),5))
neuronnames<-c("17100_129","200297_027","195613_019","195618_001","200299_052","200314_008","201765_028","192310_053",
               "195833_023","200313_091","194778_169","200199_059","192106_088","195615_018","195616_024","200199_016",
               "200582_003","200331_003","194784_007","200582_003","194786_080","193869_001","195942_002","195936_014")
plot3d(subset(PassingAxons,neuronnames),col=rep(scales::hue_pal()(3),each=8),lwd=3,soma=100)
shade3d(readOBJ(paste0("~/LS/data/allenOBJ/root.obj")),color="grey",override = T,alpha=0.05)
# par3d(userMatrix=matrix(c(1,0,0,0,0,-1,0,0,0,0,-1,0,0,0,0,1),nrow=4),windowRect=c(637,241,1703,987),zoom=0.6,FOV=0)
par3d(userMatrix=matrix(c(1,0,0,0,0,0,-1,0,0,1,0,0,0,0,0,1),nrow=4),windowRect=c(637,241,1703,987),zoom=0.6,FOV=0)
umatx<-matrix(c(0.9,0.26,-0.36,0,0.08,-0.89,-0.46,0,-0.44,0.38,-0.82,0,0,0,0,1),4,4)
par3d(windowRect = c(103,202,1396,1051),zoom=.7,userMatrix=umatx)
# rgl.snapshot(paste0(FIGUREDIR,"Examples of passing axons (",paste0(neuronnames,collapse = " "),").png"))
rgl.snapshot(paste0(FIGUREDIR,"Examples of passing axons.png"))

load("~/LS/data/arborMC")
allenarbors_dataframe$Type<-as.character(allenarbors_dataframe$Type)
allenarbors_dataframe$Type[which(allenarbors_dataframe$AxonSubtype %in% 1:44)]<-"IT"
allenarbors_dataframe$Type[which(allenarbors_dataframe$AxonSubtype %in% 45:52)]<-"CT"
allenarbors_dataframe$Type[which(allenarbors_dataframe$AxonSubtype %in% 53:64)]<-"PT"

allenarbors_dataframe %>%
  subset(AxonSubtype %in% 1:64) %>%
  ggplot(aes(x=Type))+
  geom_bar()+
  theme_classic()+
  labs(y="Number")+
  theme(plot.title = element_text(hjust = 0.5,size=8),
        text=element_text(size=5),
        axis.text.x = element_text(vjust = 0.5,size=7),axis.text.y =element_text(size=5),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        strip.background = element_blank(),strip.text.x = element_text(angle=0,size=5),
        legend.title = element_text(hjust = 0.5,size=7),legend.text = element_text(hjust = 0.5,size=7),
        legend.key.size = unit(1,"mm"),
        # legend.position = "none",
        legend.text.align=0)
ggsave(paste0(FIGUREDIR,"The distribution of terminal arbors in different neuron types.pdf"),device = "pdf",width = 30,height = 30,units = "mm")

tmp<-allenarbors_dataframe %>%
  subset(AxonSubtype %in% 1:64)
targets<-c("Isocortex","BLA","STR","PAL","TH","HY","MB","P","MY")
for(tt in targets){
  tmp$annot[which(tmp$annot %in% calannot(tt))]<-tt
}
tmp %>%
  subset(annot %in% targets) %>%
  ggplot(aes(x=annot))+
  geom_bar()+
  theme_classic()+
  labs(y="Number")+
  theme(plot.title = element_text(hjust = 0.5,size=8),
        text=element_text(size=5),
        axis.text.x = element_text(vjust = 1,hjust=1,size=5,angle=30),axis.text.y =element_text(size=5),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        strip.background = element_blank(),strip.text.x = element_text(angle=0,size=5),
        legend.title = element_text(hjust = 0.5,size=7),legend.text = element_text(hjust = 0.5,size=7),
        legend.key.size = unit(1,"mm"),
        # legend.position = "none",
        legend.text.align=0)
ggsave(paste0(FIGUREDIR,"The distribution of terminal arbors in different target regions.pdf"),device = "pdf",width = 40,height = 30,units = "mm")
