library(nat)
library(rstatix)
library(ggpubr)
library(ggpmisc)
library(fitdistrplus)
library(truncnorm)
library(snp)
library(fitdistrplus)

source("~/LS/R/Functions/Annotation.R")
source("~/LS/R/Functions/Tree.R")
FIGUREDIR<-"~/课题/Work\\ report/Figures/"

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
  if(length(arbors)>0){
    boutons$Type<-unique(arbors[,]$Type)
    test$d$Label<-10
    ad<-cbind(unlist(lapply(arbors,function(x)x$d$X[setdiff(1:nrow(x$d),x$EndPoints)])),
              unlist(lapply(arbors,function(x)x$d$Y[setdiff(1:nrow(x$d),x$EndPoints)])),
              unlist(lapply(arbors,function(x)x$d$Z[setdiff(1:nrow(x$d),x$EndPoints)])))
    ind<-which(test$d$X %in% ad[,1]&test$d$Y %in% ad[,2]&test$d$Z %in% ad[,3])
    test$d$Label[ind]<-0
    test$d$Label[1]<-1
    
    seglen<-seglengths(test)
    
    for(j in 1:tests[[i]]$NumSegs){
      boutons$Passing[which(boutons$SegID==j)]<-any(test$d$Label[test$SegList[[j]]] %in% c(10,15,20))
      boutons$Intensity_r[which(boutons$SegID==j)]<-d$W[boutons$NodeID[which(boutons$SegID==j)]]/mean(d$W[test$SegList[[j]]])
      # neuron<-Neurons[[which(names(Neurons)==names(tests)[i])]]
      # boutons$Order[which(boutons$SegID==j)]<-max(neuron$d$W[neuron$SegList[[j]]])
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
}
Boutons<-subset(Boutons,IBI>0)
Boutons[,]$Type<-as.factor(substr(Boutons[,]$Type,1,2))

ggviolin(Boutons,x="Passing",y="IBI",fill="Passing",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2)))
ggviolin(Boutons[which(Boutons$Passing==F),],x="Type",y="IBI",fill="Type",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
  stat_compare_means(label="p.format",hide.ns = T,method = "wilcox.test",comparisons = list(c(1,2),c(1,3),c(2,3)))
ggviolin(Boutons[which(Boutons$Passing==F),],x="Order",y="IBI",fill="Order",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)

######
##inter-bouton interval
load("~/LS/data/arborMC")

data<-subset(allenarbors_dataframe,!is.na(length_initial)&!is.na(BoutonNum)&intensity>1.1)
# data$length_initial<-data$length_initial/1000
data$Type<-substr(as.character(data$Type),1,2)
data<-subset(data,Type %in% c("PT","CT","IT"))
data$Type<-factor(data$Type)
ancova(length_initial~BoutonNum*Type,data=data,pch=20)

ggplot(data.frame(yy=data$length_initial,xx=data$BoutonNum,Type=factor(substr(data$Type,1,2),levels=c("IT","PT","CT"))),
       aes(x = xx, y = yy,
           fill=Type,color=Type
           )) +   
  geom_point(size=0.1,pch=20) +
  geom_smooth(method = 'lm', formula = y ~ x,lwd=0.5) +
  scale_color_manual(values = c("red","green","blue"))+
  # geom_text(aes(label=n),cex=1.5,color=colors)+
  stat_poly_eq(formula=y~x,aes(label=paste(..eq.label..,..rr.label..,..p.value.label..,sep="~~~")),p.digits = 26,parse=T,size=2)+
  # geom_tile(aes(x = xx, y = yy))+
  labs(x="Bouton Number",y="Length (um)")+
  theme_classic()+
  theme(line = element_line(size=0.1),plot.title = element_blank(),
        axis.text.y =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.text.x =element_text(size=5,angle = 0,hjust=1,vjust = 1),
        axis.title.y =element_text(size=6),axis.title.x =element_text(size=6),
        strip.text.y = element_text(angle=0,size=5),
        strip.placement = "outside",
        panel.spacing = unit(0,"mm"),
        legend.background = element_rect(colour="black",size=0.1),
        legend.title = element_text(hjust = 0.5,size=5),legend.text = element_text(hjust = 0.5,size=5)
  )
ggsave(paste0(FIGUREDIR,"The relationship between length and bouton number.pdf"),device = "pdf",width = 60,height = 50,scale = 1,units = "mm")

targets<-c("Isocortex","PAL","STR","TH","HY","MB","P","MY")
data<-subset(allenarbors_dataframe,annot %in% calannot("TH")&AxonSubtype %in% 45:64&BoutonNum!=0&SNR>0.25&Length_XYvsZ>1&intensity>1)
data$Type[which(data$AxonSubtype %in% 45:52)]<-"CT"
data$Type[which(data$AxonSubtype %in% 53:64)]<-"PT"

for(tt in targets){
  data$annot[which(data$annot %in% calannot(tt))]<-tt
}

data$Type<-substr(data$Type,1,2)
data<-data.frame(Mu=data$meanlog,
                 Type=data$Type,
                 # Region=allenarbors_dataframe$Subclass[match(rownames(data),rownames(allenarbors_dataframe))],
                 # Region=data$annot,
                 # Region=data$Laminar,
                 Sigma=data$sdlog,
                 # tail=data$Bouton_tail,
                 Size=data$Bouton_sigma
                 )
data<-subset(data,!is.na(Mu))
data<-reshape2::melt(data,id.vars = c("Type"))
# data<-reshape2::melt(data[,],id.vars = c("Region"))
colnames(data)[2]<-"features"
ggviolin(data, 
         x="Type",
         # x="Region",
         y="value",add="mean",
         xlab="morphological features",
         size=0.3,
         add.params = list(size=0.1),
         color = "Type",
         scale="width",
         width=0.8,
         # color = "Region",
         order=unique(data$Region)[order(sapply(unique(data$Region),function(x)mean(data$value[which(data$features=="Mu"&data$Region==x)])))]
         )+
  labs(x="Target regions")+
  # labs(x="Morphological subclasses")+
  facet_wrap("features",scales="free_y",
             ncol=length(unique(data$features)),
             # nrow=length(unique(data$features)),
             strip.position="top"
  ) +
  stat_compare_means(
                     # aes(group="features",label=..p.format..),
                     comparisons = list(c("PT","CT")),
                     # ref.group = ".all.",
                     aes(label=..p.format..),
                     size=2,
                     method = "wilcox.test")+
  scale_y_continuous(n.breaks = 4)+
  geom_hline(aes(yintercept = mean_value),
             data.frame(features=unique(data$features),mean_value=sapply(unique(data$features),function(x)mean(data$value[which(data$features==x)]))),
             lwd=0.5,lty=3)+
  scale_color_manual(values=c("red","blue"))+
  # scale_color_manual(values=c("#F8766D","#C49A00","#53B400","#00FFFF","#4682B4","#A58AFF"))+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.title.y =element_blank(),
        panel.spacing.y = unit(0,"mm"),
        legend.position = "none",
        legend.title = element_text(hjust = 0.5,size=5),legend.text = element_text(hjust = 0.5,size=5),
        legend.key.height=unit(0.1,"mm"),legend.key.width = unit(1,"mm"),legend.text.align=0,
        legend.box.spacing = unit(0.1,"mm")
  )
# ggsave(paste0(FIGUREDIR,"The differential features of boutons of PT arbors in different target regions (violin plot).pdf"),device = "pdf",width = 90,height = 80,scale=1,units = "mm")
ggsave(paste0(FIGUREDIR,"The differential features of boutons of PT and CT arbors in TH (violin plot).pdf"),device = "pdf",width = 90,height = 30,scale=1,units = "mm")
ggsave(paste0(FIGUREDIR,"The differential features of boutons of different cortical arbor subclasses.pdf"),device = "pdf",width = 100,height = 30,scale=1,units = "mm")

options(scipen=10)
arbornames<-c("194779_074_18")
for(i in which(names(allenarbors) %in% arbornames)){
  print(i)
  intensity<-read.neuron(paste0("/mnt/share_SNR/bouton_detection/result_upsample/",strsplit(names(allenarbors)[i],split = "_")[[1]][1],"/",strsplit(names(allenarbors)[i],split = "_")[[1]][2],"-",formatC(as.numeric(strsplit(names(allenarbors)[i],split = "_")[[1]][3]),width = 3,flag = "0"),"_upsample.swc"))
  intensity$d$Label<-0
  intensity$d$Label[1]<-1
  tmp<-paste0("/mnt/share_SNR/bouton_detection/result_fitBouton_v20/",strsplit(names(allenarbors)[i],split = "_")[[1]][1],"/",strsplit(names(allenarbors)[i],split = "_")[[1]][2],"-",formatC(as.numeric(strsplit(names(allenarbors)[i],split = "_")[[1]][3]),width = 3,flag = "0"),"_fit.txt")
  
  if(file.exists(tmp)){
    boutons<-read.csv(tmp,header = F)
    colnames(boutons)<-c("SegID","NodeID","mu","sigma","A","intensity","fg","bg")
    d<-intensity$d
    d$Label[boutons$NodeID]<-10
    d$Label[1]<-1
    d$Parent[1]<--1
    d$W<-d$W/2
    if(!dir.exists(paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/"))){
      dir.create(paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/"))
    }
    write.table(matrix(as.numeric(unlist(d)),nrow=nrow(d)),paste0("~/LS/data/raw_data/",strsplit(allenarbors[,]$NeuronName[i],split = "_")[[1]][1],"/swc_arbor_boutons/",names(allenarbors)[i],".swc"),col.names=F,row.names=F)
  }
}

##########
###the distribution of ibi
arbornames<-rownames(subset(allenarbors_dataframe,annot %in% calannot("TH")&AxonSubtype %in% 45:64&BoutonNum!=0&SNR>0.25&Length_XYvsZ>1&intensity>0.5))
IBI<-c()
for(i in unique(allenarbors[,]$NeuronName[which(names(allenarbors) %in% arbornames)])){
  print(i)
  test<-read.neuron(paste0("~/LS/data/raw_data/",strsplit(i,split = "_")[[1]][1],"/swc_allen_space/",strsplit(i,split = "_")[[1]][2],".swc"))
  tmp1<-paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um/",strsplit(i,split = "_")[[1]][1],"/",strsplit(i,split = "_")[[1]][2],"_normalize_200.swc")
  tmp2<-paste0("/mnt/SingleNeuronReconstruction/intensity_0.5um_fit_result/",strsplit(i,split = "_")[[1]][1],"/",strsplit(i,split = "_")[[1]][2],"_fit_LoG.txt")
  
  if(file.exists(tmp1)&file.exists(tmp2)){
    intensity<-read.neuron(tmp1)
    A<-unlist(lapply(intensity$SegList,function(x)mean(intensity$d$W[x]/2)))
    intensity$d$Label<-0
    intensity$d$Label[1]<-1
    seglen<-seglengths(intensity,sumsegment = F)
    
    boutons<-read.csv(tmp2,header = F)
    colnames(boutons)<-c("SegID","NodeID","mu","sigma","A")
    
    for(j in which(allenarbors[,]$NeuronName==i&names(allenarbors) %in% arbornames)){
      test$d$Label<-0
      test$d$Label[1]<-1
      ad<-allenarbors[[j]]$d
      ind<-which(test$d$X %in% ad$X&test$d$Y %in% ad$Y&test$d$Z %in% ad$Z)
      test$d$Label[ind]<-10
      segid<-unlist(lapply(test$SegList,function(x)all(test$d$Label[x] %in% 10)))
      segid<-which(segid==T)
      # allenarbors[,]$BoutonNum[j]<-length(unique(boutons$NodeID[which(boutons$SegID %in% segid)]))
      # allenarbors[,]$Bouton_sigma[j]<-mean(boutons$sigma[which(boutons$SegID %in% segid)])
      # allenarbors[,]$Bouton_A[j]<-mean(boutons$A[which(boutons$SegID %in% segid)])
      # allenarbors[,]$intensity[j]<-mean(intensity$d$W[unique(unlist(intensity$SegList[segid]))]/2)
      # allenarbors[,]$length_initial[j]<-sum(unlist(seglen[segid]))
      
      subboutons<-subset(boutons,SegID %in% segid)
      subboutons$mean_A<-NULL
           
      if(nrow(subboutons)>1){
        subboutons$ParentNodeID<-c(subboutons$NodeID[2:nrow(subboutons)],subboutons$NodeID[nrow(subboutons)])
        subboutons$ParentSegID<-c(subboutons$SegID[2:nrow(subboutons)],subboutons$SegID[nrow(subboutons)])
        subboutons$NodeID[which(subboutons$SegID==subboutons$ParentSegID)]<-apply(subboutons[which(subboutons$SegID==subboutons$ParentSegID),],1,function(x)which(intensity$SegList[[as.numeric(x[1])]]==as.numeric(x[2])))
        subboutons$ParentNodeID[which(subboutons$SegID==subboutons$ParentSegID)]<-apply(subboutons[which(subboutons$SegID==subboutons$ParentSegID),],1,function(x)which(intensity$SegList[[as.numeric(x[1])]]==as.numeric(x[6])))
        subboutons$ParentNodeID[which(subboutons$SegID!=subboutons$ParentSegID)]<-subboutons$NodeID[which(subboutons$SegID!=subboutons$ParentSegID)]
        ibi<-apply(subboutons,1,function(x)ifelse(x[2]>length(seglen[[x[1]]])|x[6]==1,0,sum(seglen[[as.numeric(x[1])]][as.numeric(x[2]):(as.numeric(x[6])-1)])))
        dist2SegHead<-apply(subboutons,1,function(x)sum(seglen[[x[1]]][1:x[2]])-sum(seglen[[x[1]]][x[2]]))
        dist2SegTail<-apply(subboutons,1,function(x)sum(seglen[[x[1]]][x[2]:length(seglen[[x[1]]])]))
        tmp<-data.frame(ibi=ibi,
                        dist2SegHead=dist2SegHead,
                        dist2SegTail=dist2SegTail,
                        SegID=subboutons$SegID,
                        arborname=names(allenarbors)[j],
                        Type=allenarbors[,]$Type[j],
                        annot=allenarbors[,]$annots[j],
                        sigma=subboutons$sigma)
        IBI<-rbind(IBI,tmp)
      }
    }
  }
}

##distribution of ibi
IBI$Type<-factor(substr(as.character(IBI$Type),1,2))
data<-subset(IBI,ibi<100&!is.na(ibi))
ggplot(data,aes(x=ibi,fill=arborname,color=Type))+
  geom_density(alpha=0)+
  xlim(0,100)+
  theme(line = element_line(size=0.1),plot.title=element_text(size=7,hjust = 0.5),
        axis.text.x = element_text(size=7),axis.text.y =element_text(size=7),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        legend.position = "none",
        legend.text = element_text(size=7),legend.title = element_blank(),
        legend.key.height =unit(1,"mm"),legend.key.width = unit(1,"mm"),
        panel.spacing = unit(0,"mm"))
IBI$somaposition<-allenarbors_dataframe$somaposition[match(a)]

targets<-c("Isocortex","PAL","STR","TH","HY","MB","P","MY")
data<-subset(allenarbors_dataframe,AxonSubtype %in% 53:64&BoutonNum!=0&annot %in% calannot(targets)&SNR>0.25&Length_XYvsZ>1&intensity>1)
data$Type[which(data$AxonSubtype %in% 53:64)]<-"PT"

for(tt in targets){
  data$annot[which(data$annot %in% calannot(tt))]<-tt
}

neuronnames<-unique(intersect(data$NeuronName[which(data$annot=="P")],data$NeuronName[which(data$annot=="PAL")]))
# arbornames<-data$arborname[which(data$NeuronName %in% neuronnames&data$annot %in% c("P","PAL"))]

bins<-10^seq(-0.2,2,length.out=20)
# bins<-10^seq(-0.2,1.5,length.out=20)

# bins<-seq(0,100,2)
tmp<-expand.grid(bins=(bins[1:(length(bins)-1)]+bins[2:length(bins)])/2,arbornames=arbornames,freq=0)
for(a in unique(IBI$arborname)){
  print(a)
  tmp$freq[which(tmp$arbornames==a)]<-sapply(1:(length(bins)-1),function(x)length(which(IBI$ibi>bins[x]&IBI$ibi<=bins[x+1]&IBI$arborname==a))/length(which(IBI$arborname==a)))
}

tmp$Type<-IBI$Type[match(tmp$arbornames,IBI$arborname)]
tmp$Type<-allenarbors_dataframe$AxonSubtype[match(tmp$arbornames,rownames(allenarbors_dataframe))]
tmp$Type[which(tmp$Type %in% 45:52)]<-"CT"
tmp$Type[which(tmp$Type %in% 53:64)]<-"PT"
# tmp$region<-data$annot[match(tmp$arbornames,data$arborname)]

mean_freq<-c()
pvalue<-c()
for(i in unique(tmp$bins)){
    pvalue<-c(pvalue,wilcox.test(as.vector(tmp$freq[which(tmp$bins==i&tmp$Type==unique(tmp$Type)[1])]),
                                 as.vector(tmp$freq[which(tmp$bins==i&tmp$Type==unique(tmp$Type)[2])]))$p.value)
    mean_freq<-c(mean_freq,pmax(mean(tmp$freq[which(tmp$bins==i&tmp$Type==unique(tmp$Type)[1])]),mean(tmp$freq[which(tmp$bins==i&tmp$Type==unique(tmp$Type)[2])])))
    # pvalue<-c(pvalue,wilcox.test(as.vector(tmp$freq[which(tmp$bins==i&tmp$region==unique(tmp$region)[1])]),
    #                              as.vector(tmp$freq[which(tmp$bins==i&tmp$region==unique(tmp$region)[2])]))$p.value)
    # mean_freq<-c(mean_freq,pmax(mean(tmp$freq[which(tmp$bins==i&tmp$region==unique(tmp$region)[1])]),mean(tmp$freq[which(tmp$bins==i&tmp$region==unique(tmp$region)[2])])))
}
pvalue<-p.adjust(pvalue,"fdr")
marker<-rep("",length(pvalue))
marker[which(pvalue<0.05)]<-"*"

ggplot(tmp,aes(x = bins,
               # fill=region,color=region,
               fill=Type,color=Type,
               y=freq
)) +
  stat_summary(
               fun = "mean",
               geom = "ribbon",
               alpha=0.3,
               # geom="pointrange",
               size=0.1,
               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))+
  stat_summary(
               fun="mean",
               geom="line",lwd=0.2)+
  annotate("text",x=(bins[1:(length(bins)-1)]+bins[2:length(bins)])/2,
  y=mean_freq+0.01,
  # size=2,
  label=marker)+
  scale_x_log10()+
  theme_bw()+
  labs(x="IBI (um)",y="Frequency")+
  scale_y_continuous(n.breaks = 3)+
  scale_fill_manual(values=c("red","blue"))+
  scale_color_manual(values=c("red","blue"))+
  theme(line = element_line(size=0.3),plot.title=element_text(size=7,hjust = 0.5),
        axis.text.x = element_text(size=7),axis.text.y =element_text(size=5),
        axis.title.x =element_text(size=7),axis.title.y =element_text(size=7),
        legend.position = "right",
        legend.text = element_text(size=5),legend.title = element_blank(),
        legend.key.height =unit(1,"mm"),legend.key.width = unit(1,"mm"),
        panel.spacing = unit(0,"mm"))
ggsave(paste0(FIGUREDIR,"The IBI distribution of boutons of PT and CT arbors in thalamus.pdf"),device = "pdf",width = 50,height = 50,scale=1,units = "mm")
ggsave(paste0(FIGUREDIR,"The IBI distribution of boutons of PAL and P arbors in thalamus.pdf"),device = "pdf",width = 60,height = 40,scale=1,units = "mm")

IBI$regions<-IBI$annot
# IBI$regions[which(IBI$annot %in% calannot("STR"))]<-"STR"
IBI$regions[which(IBI$annot %in% calannot("P"))]<-"P"
IBI$regions[which(IBI$annot %in% calannot("PAL"))]<-"PAL"

# IBI$regions[which(IBI$annot %in% calannot("AId"))]<-"AId"
# IBI$regions[which(IBI$annot %in% calannot("ACA"))]<-"ACA"
# IBI$regions[which(IBI$annot %in% calannot("PL"))]<-"PL"
ggplot(subset(IBI),aes(x=ibi,fill=arborname,color=regions))+
  geom_density(alpha=0,lwd=0.3)+
  xlim(0,50)+
  theme_classic()+
  theme(line = element_line(size=0.001),
        axis.text.y =element_text(size=5),axis.text.x =element_text(size=5),
        axis.title.y =element_text(size=5),axis.title.x =element_text(size=5),
        axis.line = element_line(size=0.1),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0,size=5),
        strip.text.x=element_text(size=5),
        strip.placement = "outside",
        # legend.position = "none",
        legend.position = "right",
        legend.text = element_text(size=7),legend.title = element_blank(),
        legend.key.height =unit(1,"mm"),legend.key.width = unit(1,"mm"),
        panel.spacing = unit(0,"mm"))
ggsave(paste0(FIGUREDIR,"The distribution of ibi in PAL and P.pdf"),device = "pdf",width = 80,height = 40,scale = 1,units = "mm")

ggplot(IBI,aes(x=ibi,fill=arborname,color=somaposition))+
  geom_density(alpha=0)+
  xlim(0,50)

############the boutons clusters
IBI$clusters<-NA
n<-1
a<-IBI$arborname[1]
ss<-IBI$SegID[1]
for(i in 1:nrow(IBI)){
  print(i)
  IBI$clusters[i]<-n
  if(IBI$ibi[i]>=7|IBI$ibi[i]==0){
    n<-n+1
  }
  if(IBI$arborname[i]!=a|IBI$SegID[i]!=ss){
    n<-1
  }
  a<-IBI$arborname[i]
  ss<-IBI$SegID[i]
}

##the number of boutons within clusters
cluster_n<-c()
type<-c()
arborname<-c()
for(i in unique(IBI$arborname)){
  print(i)
  for(j in unique(IBI$SegID[which(IBI$arborname==i)])){
    tmp<-summary(factor(IBI$clusters[which(IBI$arborname==i&IBI$SegID==j)]),maxsum=length(unique(factor(IBI$clusters[which(IBI$arborname==i&IBI$SegID==j)]))))
    cluster_n<-c(cluster_n,tmp+1)
    type<-c(type,rep(unique(as.character(IBI$Type[which(IBI$arborname==i)])),length(tmp)))
    arborname<-c(arborname,rep(i,length(tmp)))
  }
}

data<-data.frame(cluster_n=cluster_n,type=type,arborname=arborname)
data<-subset(data,cluster_n<=15)
bins_low<-c(2,3,4,10)
bins_high<-c(2,3,9,15)
data$cluster_n_bins<-NA
for(i in 1:length(bins_low)){
  data$cluster_n_bins[which(data$cluster_n>=bins_low[i]&data$cluster_n<=bins_high[i])]<-paste0(bins_low[i],"-",bins_high[i])
}
data<-data %>% group_by(arborname,type) %>% count(cluster_n_bins)
data$freq<-sapply(1:nrow(data),function(x)data$n[x]/sum(data$n[which(data$arborname==data$arborname[x])]))
data$cluster_n_bins<-factor(data$cluster_n_bins,levels=paste0(bins_low,"-",bins_high))
ggboxplot(data,x="cluster_n_bins",y="freq",color="type",size=0.1,outlier.size=0.1)+
  stat_compare_means(aes(group=type),
                     label="p.signif",
                     method = "wilcox.test",
                     label.y = 0.9,
                     hide.ns = T,
                     size=2)+
  scale_color_manual(values=c("red","blue"))+
  labs(x="Size of bouton clusters",y="Proportion")+
  theme(line = element_line(size=0.5),
        axis.text.y =element_text(size=5),axis.text.x =element_text(size=5),
        axis.title.y =element_text(size=7),axis.title.x =element_text(size=7),
        axis.ticks.x = element_line(size=0.3),
        axis.ticks.y = element_line(size=0.3),
        axis.line = element_line(size=0.3),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0,size=5),
        strip.text.x=element_text(size=5),
        legend.position = "right",
        legend.text = element_text(size=7),legend.title = element_blank(),
        legend.key.height =unit(1,"mm"),legend.key.width = unit(1,"mm"),
        panel.spacing = unit(0,"mm"))
ggsave(paste0(FIGUREDIR,"The number of boutons within clusters of CT and PT neurons in thalamus.pdf"),device = "pdf",width = 60,height = 50,scale = 1,units = "mm")

#########fitted by exponential or gamma distribution
dtnorm<- function(x, mean, sd, a, b) {
  dnorm(x, mean, sd)/(pnorm(b, mean, sd)-pnorm(a, mean, sd))
}
ptnorm <- function(q, mean, sd, a, b) {
  (pnorm(q,mean,sd) - pnorm(a,mean,sd)) / 
    (pnorm(b,mean,sd) - pnorm(a,mean,sd))
}
qtnorm<-truncnorm::qtruncnorm

P_exp<-c()
rate_exp<-c()
P_gamma<-c()
rate_gamma<-c()
shape_gamma<-c()
P_lnorm<-c()
meanlog<-c()
sdlog<-c()
Y_exp<-c()
Y_gamma<-c()
Y_lnorm<-c()
ibi<-c()
arborname<-c()
# for(a in unique(IBI$arborname)){
for(a in "194526_001_33"){
  tmp<-subset(IBI,arborname %in% a)
  ibi<-c(ibi,tmp$ibi)
  arborname<-c(arborname,rep(a,length(tmp$ibi)))
  mod<-fitdist(tmp$ibi,"exp")
  rate_exp<-c(rate_exp,mod$estimate[1])
  Y_exp<-c(Y_exp,dexp(tmp$ibi,rate=mod$estimate[1]))
  # P_exp<-c(P_exp,ks.test(tmp$ibi,"pexp",rate=mod$estimate[1])$p.value)
  P_exp<-c(P_exp,gofstat(mod)$ks)
  
  # mod<-fitdist(tmp$ibi,dtnorm,
  #              start=list(mean=mean(tmp$ibi),sd=sd(tmp$ibi)),
  #              fix.arg = list(a=0,b=Inf))
  
  mod<-fitdist(tmp$ibi,"lnorm")
  meanlog<-c(meanlog,mod$estimate[1])
  sdlog<-c(sdlog,mod$estimate[2])
  Y_lnorm<-c(Y_lnorm,dlnorm(tmp$ibi,meanlog=mod$estimate[1],sdlog=mod$estimate[2]))
  # P_lnorm<-c(P_lnorm,ks.test(tmp$ibi,"plnorm",meanlog=mod$estimate[1],sdlog=mod$estimate[2])$p.value)
  P_lnorm<-c(P_lnorm,gofstat(mod)$ks)
  
  mod<-fitdist(tmp$ibi,"gamma")
  shape_gamma<-c(shape_gamma,mod$estimate[1])
  rate_gamma<-c(rate_gamma,mod$estimate[2])
  Y_gamma<-c(Y_gamma,dgamma(tmp$ibi,shape=mod$estimate[1],rate=mod$estimate[2]))
  # P_gamma<-c(P_gamma,ks.test(tmp$ibi,"pgamma",shape=mod$estimate[1],rate=mod$estimate[2])$p.value)
  P_gamma<-c(P_gamma,gofstat(mod)$ks)
}

tmp<-paste0(unique(arborname),
            "\\n","P value: ",format(P_exp,digits = 2)," (Exp)",
            "\\n",format(P_gamma,digits = 2)," (Gamma)",
            "\\n",format(P_lnorm,digits = 2)," (Lnorm)")
names(tmp)<-unique(arborname)
ggplot(data.frame(ibi=rep(ibi,4),
                  fit=c(rep(0,length(Y_exp)),Y_exp,Y_gamma,Y_lnorm),
                  arborname=rep(arborname,4),
                  data=rep(c("Observed","Fitted, exponential","Fitted, gamma","Fitted, lnorm"),each=length(Y_exp))),
       aes(x=ibi,fill=data))+
  geom_histogram(aes(y=..density..),alpha=1,position = "identity")+
  scale_fill_manual(values=c("green","red","blue","grey"))+
  geom_line(aes(x=ibi,y=fit,color=data))+
  scale_color_manual(values=c("green","red","blue","grey"))+
  # scale_alpha_manual(values=c(1,0,0))+
  facet_wrap("arborname",labeller = labeller(arborname=tmp))+
  xlim(c(0,50))+
  labs(x="Inter-bouton intervals")+
  theme_bw()+
  theme(plot.title=element_text(size=10,hjust = 0.5),
        text=element_text(size=8),
        axis.text.x = element_text(size=5),axis.text.y = element_text(size=5),
        axis.title.x = element_text(size=7),axis.title.y = element_text(size=7),
        legend.title = element_blank(),legend.text = element_text(hjust = 0.5,size=7),
        legend.key.height=unit(1,"mm"),legend.key.width = unit(1,"mm"),legend.text.align=0
  )

P_exp<-c()
rate_exp<-c()
P_gamma<-c()
rate_gamma<-c()
shape_gamma<-c()
P_lnorm<-c()
meanlog<-c()
sdlog<-c()
Y_exp<-c()
Y_gamma<-c()
Y_lnorm<-c()
sigma<-c()
arborname<-c()
for(a in unique(IBI$arborname)){
  tmp<-subset(IBI,arborname %in% a)
  sigma<-c(sigma,tmp$sigma)
  arborname<-c(arborname,rep(a,length(tmp$sigma)))
  mod<-fitdist(tmp$sigma,"exp")
  rate_exp<-c(rate_exp,mod$estimate[1])
  Y_exp<-c(Y_exp,dexp(tmp$sigma,rate=mod$estimate[1]))
  # P_exp<-c(P_exp,ks.test(tmp$sigma,"pexp",rate=mod$estimate[1])$p.value)
  P_exp<-c(P_exp,gofstat(mod)$ks)
  
  # mod<-fitdist(tmp$sigma,dtnorm,
  #              start=list(mean=mean(tmp$sigma),sd=sd(tmp$sigma)),
  #              fix.arg = list(a=0,b=Inf))
  
  mod<-fitdist(tmp$sigma,"lnorm")
  meanlog<-c(meanlog,mod$estimate[1])
  sdlog<-c(sdlog,mod$estimate[2])
  Y_lnorm<-c(Y_lnorm,dlnorm(tmp$sigma,meanlog=mod$estimate[1],sdlog=mod$estimate[2]))
  # P_lnorm<-c(P_lnorm,ks.test(tmp$sigma,"plnorm",meanlog=mod$estimate[1],sdlog=mod$estimate[2])$p.value)
  P_lnorm<-c(P_lnorm,gofstat(mod)$ks)
  
  mod<-fitdist(tmp$sigma,"gamma")
  shape_gamma<-c(shape_gamma,mod$estimate[1])
  rate_gamma<-c(rate_gamma,mod$estimate[2])
  Y_gamma<-c(Y_gamma,dgamma(tmp$sigma,shape=mod$estimate[1],rate=mod$estimate[2]))
  # P_gamma<-c(P_gamma,ks.test(tmp$sigma,"pgamma",shape=mod$estimate[1],rate=mod$estimate[2])$p.value)
  P_gamma<-c(P_gamma,gofstat(mod)$ks)
}

tmp<-paste0(unique(arborname),
            "\\n",format(P_exp,digits = 2)," (Exp)",
            "\\n",format(P_gamma,digits = 2)," (Gamma)",
            "\\n",format(P_lnorm,digits = 2)," (Lnorm)")
names(tmp)<-unique(arborname)
ggplot(data.frame(sigma=rep(sigma,4),
                  fit=c(rep(0,length(Y_exp)),Y_exp,Y_gamma,Y_lnorm),
                  arborname=rep(arborname,4),
                  data=rep(c("Observed","Fitted, exponential","Fitted, gamma","Fitted, lnorm"),each=length(Y_exp))),
       aes(x=sigma,fill=data))+
  geom_histogram(aes(y=..density..),alpha=1,position = "identity")+
  scale_fill_manual(values=c("green","red","blue","grey"))+
  geom_line(aes(x=sigma,y=fit,color=data))+
  scale_color_manual(values=c("green","red","blue","grey"))+
  # scale_alpha_manual(values=c(1,0,0))+
  facet_wrap("arborname",labeller = labeller(arborname=tmp))

data<-data.frame(param=rate_exp,annot=IBI$annot[match(unique(IBI$arborname),IBI$arborname)])
data$annot[which(data$annot %in% calannot("CTXpl"))]<-"CTXpl"
data$annot[which(data$annot %in% calannot("STR"))]<-"STR"
stat.test <- data %>%
  wilcox_test(param ~ annot) %>%
  adjust_pvalue() %>%
  add_xy_position()
# ggviolin(data,x="Type",y="ibi",fill="Type",add="mean",add.params = list(size=0.01,color="black"),size=0.1,width=1)+
ggboxplot(data,x="annot",y="param",fill="annot",add="dotplot",add.params = list(size=0.5,fill="black"),size=0.1,width=1)+
  stat_pvalue_manual(stat.test, label = "p.adj")+
  scale_y_continuous(n.breaks = 3)+
  labs(x="region",y="rate_exp")+
  # scale_fill_manual(values=c("sienna","green","magenta"),name="Class")+
  theme(line = element_line(size=0.001),
        axis.text.y =element_text(size=7),axis.text.x =element_text(size=7),
        axis.title.y =element_text(size=8),axis.title.x =element_text(size=8),
        axis.line = element_line(size=0.1),
        strip.background = element_blank(),strip.text.y = element_text(angle=0,size=7),
        strip.placement = "outside",
        panel.spacing = unit(0,"mm"),
        legend.position = "none")

ss<-which(!is.na(allenarbors_dataframe$P_exp))
tmp<-data.frame(P=c(allenarbors_dataframe$P_exp[ss],
                    allenarbors_dataframe$P_gamma[ss],
                    allenarbors_dataframe$P_lnorm[ss]),
                Fit=rep(c("Exp","Gamma","Lnorm"),each=length(ss))
)
ggplot(tmp,aes(x=P,fill=Fit,color=Fit))+
  geom_density(alpha=0.3,aes(y=..density..))+
  geom_vline(xintercept = 0.01,lty=2)+
  scale_x_log10(breaks=c(1e-13,1e-9,1e-5,1e-2))+
  # scale_x_continuous(breaks=c(0,0.05,0.5),limits = c(0,0.5))+
  scale_y_continuous(n.breaks = 3)+
  theme_classic()+
  labs(x="P value")+
  scale_fill_manual(values=c("green","red","blue"))+
  scale_color_manual(values=c("green","red","blue"))+
  theme(axis.text.y =element_text(size=5),axis.text.x =element_text(size=5),
        axis.title.y =element_text(size=5),axis.title.x =element_text(size=7),
        # axis.line = element_line(size=0.1),
        panel.spacing = unit(0,"mm"),
        legend.position = "none")
ggsave(paste0(FIGUREDIR,"The p values of fitting of IBI with log-normal, gamma and exp distribution.pdf"),device = "pdf",width = 55,height =25,scale = 1,units = "mm")
