### Phenovern figures
# Feb 4 2019

library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
library(viridis)
library(wesanderson)
library(cowplot)
library(ggsci)
library(ggrepel)

setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/Phenovern3")

GetFormattedPhenotypes<-function(){
  Comb<-fread("Phenovern.combined.csv",data.table = FALSE)
  
  Comb$FirstFruit<-NULL
  Comb$GrandmaternalTemp<-gsub("Low","Cool",Comb$GrandmaternalTemp)
  Comb$GrandmaternalTemp<-gsub("High","Warm",Comb$GrandmaternalTemp)
  
  
  Comb<-Comb[,-c(10:11)]
  Comb<-reshape2::melt(Comb,id.vars=colnames(Comb)[c(1:3,9)])
  colnames(Comb)[5:6]<-c("Comb","JulianDays")
  Comb$JulianDays<-as.numeric(Comb$JulianDays)
  Comb$GrandmaternalTemp<-factor(Comb$GrandmaternalTemp,levels=c("Cool","Warm"))
  Comb$Generation<-factor(Comb$Generation,levels=c("Progeny","Grandprogeny"))
  #Comb$Comb<-gsub("FirstFruit","First fruit",Comb$Comb)
  Comb$Comb<-gsub("FirstYellowFruit","Days to first mature fruit",Comb$Comb)
  Comb$Comb<-gsub("Flower","Days to flowering",Comb$Comb)
  Comb$Comb<-gsub("MaxLeafNumber","Leaf number at bolting",Comb$Comb)
  Comb$Comb<-gsub("Bolting","Days to bolting",Comb$Comb)
  Comb$Comb<-gsub("MaxLength","Length of largest leaf at bolting (cm)",Comb$Comb)
  Comb$Comb<-gsub("FruitInterval","Flowering - First mature fruit interval",Comb$Comb)
  Comb$Comb<-gsub("FloweringInterval","Bolting - flowering interval",Comb$Comb)
  #Comb$Comb<-gsub("YFInterval","First fruit - first mature fruit interval",Comb$Comb)
  Comb$Generation<-gsub("Grandprogeny","3rd",Comb$Generation)
  Comb$Generation<-gsub("Progeny","2nd",Comb$Generation)
  return(Comb)
}
Comb<-GetFormattedPhenotypes()

file="comparison.emmeans.temp.csv"
FName<-gsub("comparison.emmeans.","",file,fixed=TRUE)
FName<-gsub(".csv","",FName,fixed=TRUE)

ems<-fread(file,data.table = FALSE)
ems$Phenotype<-gsub("FirstYellowFruit","Days to first mature fruit",ems$Phenotype)
ems$Phenotype<-gsub("Flower","Days to flowering",ems$Phenotype)
ems$Phenotype<-gsub("Bolting","Days to bolting",ems$Phenotype)
ems$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems$Phenotype)
ems$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting (cm)",ems$Phenotype)
ems$Generation<-gsub("Grandprogeny","3rd",ems$Generation)
ems$Generation<-gsub("Progeny","2nd",ems$Generation)

ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
ems$sig<-ems$p.value<=0.05
ems$GrandmaternalVern

ems$ID<-paste(ems$Genotype,ems$GrandmaternalVern,ems$Generation,ems$Phenotype,sep=":")
colnames(Comb)[5]<-"Phenotype"
Comb$ID<-paste(Comb$Genotype,Comb$GrandmaternalVern,Comb$Generation,Comb$Phenotype,sep=":")
Comb.plot<-merge(Comb,ems[,c(10,9)],by="ID",all.x=TRUE)
Comb.plot$GrandmaternalVern<-gsub("Unvernalized","Not vernalized",Comb.plot$GrandmaternalVern)

Phenos<-unique(Comb$Phenotype)

FigMaker.phenology.tweaked<-function(i){
  ThePheno<-Phenos[i]
  print(Phenos[i])
  Subset<-dplyr::filter(Comb,Phenotype == Phenos[i])
  colnames(Subset)[6]<-"outcome"
  ems2<-dplyr::filter(ems,Phenotype == Phenos[i])
  
  Subset$GrandmaternalVern<-gsub("Unvernalized","Not vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Unvernalized","Not vern",ems2$GrandmaternalVern)
  Subset$GrandmaternalVern<-gsub("Vernalized","Vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Vernalized","Vern",ems2$GrandmaternalVern)
  #ems2$Generation<-paste(ems2$Generation,"generation")
  
  RAW<-ggplot(Subset,aes(x=Generation,y=outcome,fill=GrandmaternalTemp)) +
    #geom_line(aes(x=(Generation))) + 
    #geom_hline(yintercept = 0,alpha=0.25) +
    #guides(size=FALSE) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    geom_boxplot(aes(x=as.factor(Generation)),alpha=0.5) +
    facet_grid(GrandmaternalVern~Genotype) + 
    scale_fill_hc(name=NULL) + #scale_fill_hc() + 
    theme_bw(base_size = 14) +
    labs(y=Phenos[i],x=NULL,title = "A. Genotypic mean phenotype for each experimental treatment") +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
      #axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      legend.position = "bottom"
      #axis.ticks.x=element_blank()
    ) +
    theme(plot.margin = margin(6, 0, 6, 0))
  
  GLMnegTemp<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalVern,group=GrandmaternalVern,linetype=GrandmaternalVern)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=2)) + 
    scale_shape_manual(values = c(1,16), guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE,shape=FALSE,fill=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid(.~Genotype,scales="free") + 
    scale_color_hc("darkunica",name=NULL) +
    theme_bw(base_size = 14) +
    labs(y="Response to thermal regime (Effect size of warm-cool)",title = "B. Response to first-generation thermal regime") +
    theme(axis.title.x=element_blank(),
                #axis.text.x=element_blank(),
                #axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=9),
          legend.position = "bottom",
          axis.title.y = element_text(size = rel(0.8), angle = 90))  +
    theme(plot.margin = margin(6, 0, 6, 0))
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  
  file="comparison.emmeans.vern.csv"
  FName<-gsub("comparison.emmeans.","",file,fixed=TRUE)
  FName<-gsub(".csv","",FName,fixed=TRUE)
  
  ems<-fread(file,data.table = FALSE)
  ems$Phenotype<-gsub("FirstYellowFruit","Days to first mature fruit",ems$Phenotype)
  ems$Phenotype<-gsub("Flower","Days to flowering",ems$Phenotype)
  ems$Phenotype<-gsub("Bolting","Days to bolting",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting (cm)",ems$Phenotype)
  ems$Generation<-gsub("Grandprogeny","3rd",ems$Generation)
  ems$Generation<-gsub("Progeny","2nd",ems$Generation)
  ems$GrandmaternalTemp<-gsub("High","Warm",ems$GrandmaternalTemp)
  ems$GrandmaternalTemp<-gsub("Low","Cool",ems$GrandmaternalTemp)
  
  
  ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
  ems$sig<-ems$p.value<=0.05
  #ems$GrandmaternalVern
  
  ems$ID<-paste(ems$Genotype,ems$GrandmaternalTemp,ems$Generation,ems$Phenotype,sep=":")
  colnames(Comb)[5]<-"Phenotype"
  Comb$ID<-paste(Comb$Genotype,Comb$GrandmaternalTemp,Comb$Generation,Comb$Phenotype,sep=":")
  Comb.plot<-merge(Comb,ems[,c(10,9)],by="ID",all.x=TRUE)
  Comb.plot$GrandmaternalVern<-gsub("Unvernalized","Not vernalized",Comb.plot$GrandmaternalVern)
  
  ems2<-dplyr::filter(ems,Phenotype == Phenos[i])
  #ems2$Generation<-paste(ems2$Generation,"generation")
  
  
  GLMnegVern<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalTemp,group=GrandmaternalTemp,linetype=GrandmaternalTemp)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=2)) + 
    scale_shape_manual(values = c(1,16),  guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid(.~Genotype,scales="free") + 
    scale_color_hc(name=NULL) +  scale_fill_tableau(guide=FALSE) + 
    theme_bw(base_size = 14) +
    labs(y="Response to vernalization (Effect size of NV-V)",title = "C. Response to first-generation vernalization treatment") +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
          #axis.text.y=element_text(size=9),
          legend.position = "bottom"
          #axis.title.y = element_text(size = rel(0.8), angle = 90))  +
    #theme(plot.margin = margin(6, 0, 6, 0)
           )
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  pdf(paste(Phenos[i],".RxN.tweaked.pdf",sep=""),height=15,width=9)
  print(
  plot_grid(RAW + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
            GLMnegTemp+ theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
            GLMnegVern+ theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
            nrow=3,
            #align="vh",
            #labels = c("A. Genotypic mean phenotype for each experimental treatment",
            #           "B. Response to first-generation thermal regime",
            #           "C. Response to first-generation vernalization treatment"),
            rel_heights = c(1,1,1.2),
            hjust = -1)
  )
  dev.off()
  
  library(cowplot)
  #pdf(paste(Phenos[i],".RxN.vern.pdf",sep=""),height=13,width=10)
  ToPrint<-plot_grid(RAW + theme(legend.position="none") ,GLMnegTemp +theme(legend.position="none") ,
                     GLMnegVern +theme(legend.position="none") ,
                     nrow=3,
                     #align="vh",
                     labels = c("A", "B","C"),
                     rel_heights = c(1,1,1.2),
                     hjust = -1)
  legend1 <- get_legend(
    # create some space to the left of the legend
    RAW + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend2 <- get_legend(
    # create some space to the left of the legend
    GLMnegTemp + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend3 <- get_legend(
    # create some space to the left of the legend
    GLMnegVern + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  #print(plot_grid(ToPrint, legend1,legend2,legend3, ncol = 1, rel_heights = c(1, .05,.05,.05)))
  #dev.off()
}

lapply(1:5,FigMaker.phenology.tweaked)
#FigMaker.ln(1)
#FigMaker.ll(2)


####### Germination figures
library(data.table)
library(ggplot2)
library(ggthemes)
library(gridExtra)

GetGerm<-function(){
  library(data.table)
  setwd("/Users/marianoalvarez/Dropbox (Personal)/R/2017PhenoVern/PhenoGerm3")
  Combined<-fread("GenCompare.csv",data.table = FALSE)
  
  table(Combined$GrandMaternalVern)
  Combined$GrandMaternalVern<-gsub("Noo","No",Combined$GrandMaternalVern)
  Combined$GrandMaternalVern<-gsub("Yeses","Yes",Combined$GrandMaternalVern)
  sapply(Combined,table)
  
  sapply(Combined,class)
  Combined$Max<-as.numeric(Combined$Max)
  #Combined[1047,6]<-0.95
  
  
  table(Combined$GrandMaternalTemp,Combined$GrandMaternalVern,Combined$Genotype,Combined$Generation)
  table(Combined$GrandMaternalTemp,Combined$GrandMaternalVern,Combined$Genotype,Combined$Generation,Combined$GermTemp)
  check<-data.frame(table(Combined$GrandMaternalTemp,Combined$GrandMaternalVern,Combined$Genotype,Combined$Generation,Combined$GermTemp))
  check$Var3<-as.character(check$Var3)
  ToRemove<-unique(dplyr::filter(check, Freq < 3)$Var3)
  
  
  library(dplyr)
  #Combined<-dplyr::filter(Combined,GermTemp=="Cold")
  Combined<-dplyr::filter(Combined,!Genotype %in% ToRemove)
  
  Combined$GrandMaternalVern<-gsub("Y","Vernalized",Combined$GrandMaternalVern)
  Combined$GrandMaternalVern<-gsub("N","Not vernalized",Combined$GrandMaternalVern)
  
  ## Viz
  
  
  Combined$Generation<-gsub("Two","2nd",Combined$Generation)
  Combined$Generation<-gsub("Three","3rd",Combined$Generation)
  Combined$GermTemp<-gsub("Cold","10C",Combined$GermTemp)
  Combined$GermTemp<-gsub("Hot","22C",Combined$GermTemp)
  Combined$Generation<-factor(Combined$Generation,levels=c("2nd","3rd"))
  return(Combined)
}
Germ<-GetGerm()

Germ$GrandMaternalTemp<-gsub("Cold","Cool",Germ$GrandMaternalTemp)

ems<-fread("comparison.gen12.emMeans.temp.csv",data.table = FALSE)
ems$Generation<-gsub("Grandprogeny","3rd",ems$Generation)
ems$Generation<-gsub("Progeny","2nd",ems$Generation)
ems$GrandmaternalVern<-gsub("Unvernalized","Not vernalized",ems$GrandmaternalVern)
#ems$GerminationTemp<-paste(ems$GerminationTemp,"germination")

#ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
ems$sig<-ems$p.value<=0.05

ems$ID<-paste(ems$Genotype,ems$GrandmaternalVern,ems$Generation,ems$GerminationTemp,sep=":")
Germ$ID<-paste(Germ$Genotype,Germ$GrandMaternalVern,Germ$Generation,Germ$GermTemp,sep=":")


Subset<-Germ
colnames(Subset)[6]<-"outcome"
#ems$estimate<- -ems$estimate
ems2<-ems

FigMaker.germ<-function(){
  #ThePheno<-Phenos[i]
  #print(Phenos[i])
  Subset<-Germ
  colnames(Subset)[6]<-"outcome"
  ems2<-ems
  
  colnames(Subset)<-gsub("Matern","matern",colnames(Subset))
  
  Subset$GrandmaternalVern<-gsub("Not vernalized","Not vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Not vernalized","Not vern",ems2$GrandmaternalVern)
  Subset$GrandmaternalVern<-gsub("Vernalized","Vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Vernalized","Vern",ems2$GrandmaternalVern)
  
  RAW<-ggplot(Subset,aes(x=Generation,y=outcome,fill=GrandmaternalTemp)) +
    #geom_line(aes(x=(Generation))) + 
    #geom_hline(yintercept = 0,alpha=0.25) +
    #guides(size=FALSE) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    geom_boxplot(aes(x=as.factor(Generation)),alpha=0.5) +
    facet_grid(GrandmaternalVern+GermTemp~Genotype) + 
    scale_fill_hc(name=NULL) + #scale_fill_hc() + 
    theme_bw(base_size = 14) +
    labs(y="Germination proportion",x=NULL) +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x=element_blank(),
      axis.text.x=element_blank()
      #axis.ticks.x=element_blank()
    ) +
    theme(plot.margin = margin(6, 0, 6, 0))
  
  GLMnegTemp<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalVern,group=GrandmaternalVern,linetype=GrandmaternalVern)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=sig)) + scale_shape_manual(values = c(20, 8), guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid("Thermal regime"+GerminationTemp~Genotype,scales="free") + 
    scale_color_hc(name=NULL) + scale_fill_tableau(guide=FALSE) + 
    theme_bw(base_size = 14) +
    labs(y="Effect size (warm-cool)") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=9),
          axis.title.y = element_text(size = rel(0.8), angle = 90))  +
    theme(plot.margin = margin(6, 0, 6, 0))
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  
  ems<-fread("comparison.gen12.emMeans.vern.csv",data.table = FALSE)
  ems$Generation<-gsub("Grandprogeny","Third",ems$Generation)
  ems$Generation<-gsub("Progeny","Second",ems$Generation)
  ems$GrandmaternalTemp<-gsub("Cold","Low",ems$GrandmaternalTemp)
  ems$GrandmaternalTemp<-gsub("Warm","High",ems$GrandmaternalTemp)
  
  #ems$GerminationTemp<-paste(ems$GerminationTemp,"germination")
  
  
  ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
  ems$sig<-ems$p.value<=0.05
  #ems$GrandmaternalVern
  
  Germ$GrandMaternalTemp<-gsub("Cool","Low",Germ$GrandMaternalTemp)
  Germ$GrandMaternalTemp<-gsub("Warm","High",Germ$GrandMaternalTemp)
  
  
  ems$ID<-paste(ems$Genotype,ems$GrandmaternalTemp,ems$Generation,ems$GerminationTemp,sep=":")
  Germ$ID<-paste(Germ$Genotype,Germ$GrandMaternalTemp,Germ$Generation,Germ$GermTemp,sep=":")
  
  Comb.plot<-merge(Germ,ems[,c(13,12)],by="ID",all.x=TRUE)
  #Comb.plot$GrandmaternalVern<-gsub("Unvernalized","Not vernalized",Comb.plot$GrandmaternalVern)
  ems$estimate<- -ems$estimate
  ems2<-ems
  
  GLMnegVern<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalTemp,group=GrandmaternalTemp,linetype=GrandmaternalTemp)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=sig)) + scale_shape_manual(values = c(20, 8), guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid("Vernalization"+GerminationTemp~Genotype,scales="free") + 
    scale_color_hc("darkunica",name=NULL) + scale_fill_tableau(guide=FALSE) + 
    theme_bw(base_size = 14) +
    labs(y="Effect size (NV-V)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y=element_text(size=9),
          axis.title.y = element_text(size = rel(0.9), angle = 90))  +
    theme(plot.margin = margin(6, 0, 6, 0))
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  library(cowplot)
  pdf("Germ.RxN.vern.pdf",height=13,width=10)
  ToPrint<-plot_grid(RAW + theme(legend.position="none") ,GLMnegTemp +theme(legend.position="none") ,
                     GLMnegVern +theme(legend.position="none") ,
                     nrow=3,
                     #align="vh",
                     labels = c("A", "B","C"),
                     rel_heights = c(1,1,1.2),
                     hjust = -1)
  legend1 <- get_legend(
    # create some space to the left of the legend
    RAW + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend2 <- get_legend(
    # create some space to the left of the legend
    GLMnegTemp + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend3 <- get_legend(
    # create some space to the left of the legend
    GLMnegVern + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  print(plot_grid(ToPrint, legend1,legend2,legend3, ncol = 1, rel_heights = c(1, .05,.05,.05)))
  dev.off()
}
#FigMaker.germ()
FigMaker.germ.tweaked<-function(i){
  #ThePheno<-Phenos[i]
  #print(Phenos[i])
  Subset<-Germ
  colnames(Subset)[6]<-"outcome"
  #ems2<-dplyr::filter(ems,Phenotype == Phenos[i])
  colnames(Subset)[c(1,3)]<-c("GrandmaternalTemp","GrandmaternalVern")
  
  Subset$GrandmaternalVern<-gsub("Not vernalized","Not vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Not vernalized","Not vern",ems2$GrandmaternalVern)
  Subset$GrandmaternalVern<-gsub("Vernalized","Vern",Subset$GrandmaternalVern)
  ems2$GrandmaternalVern<-gsub("Vernalized","Vern",ems2$GrandmaternalVern)
  #ems2$Generation<-paste(ems2$Generation,"generation")
  
  RAW<-ggplot(Subset,aes(x=Generation,y=outcome,fill=GrandmaternalTemp)) +
    #geom_line(aes(x=(Generation))) + 
    #geom_hline(yintercept = 0,alpha=0.25) +
    #guides(size=FALSE) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    geom_boxplot(aes(x=as.factor(Generation)),alpha=0.5) +
    facet_grid(GrandmaternalVern+GermTemp~Genotype) + 
    scale_fill_hc(name=NULL) + #scale_fill_hc() + 
    theme_bw(base_size = 14) +
    labs(y="Germination proportion",x=NULL,title = "A. Genotypic mean phenotype for each experimental treatment") +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      legend.position = "bottom"
      #axis.ticks.x=element_blank()
    ) +
    theme(plot.margin = margin(6, 0, 6, 0))
  
  GLMnegTemp<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalVern,group=GrandmaternalVern,linetype=GrandmaternalVern)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=2)) + 
    scale_shape_manual(values = c(1,16), guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid(GerminationTemp~Genotype,scales="free") + 
    scale_color_hc("darkunica",name=NULL) +
    scale_fill_tableau(guide=FALSE) + 
    theme_bw(base_size = 14) +
    labs(y="Response to thermal regime (Effect size of warm-cool)",title = "B. Response to first-generation thermal regime") +
    theme(axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.ticks.x=element_blank(),
          axis.text.y=element_text(size=9),
          legend.position = "bottom",
          axis.title.y = element_text(size = rel(0.8), angle = 90))  +
    theme(plot.margin = margin(6, 0, 6, 0))
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  
  ems<-fread("comparison.gen12.emMeans.vern.csv",data.table = FALSE)
  ems$Phenotype<-gsub("FirstYellowFruit","Days to first mature fruit",ems$Phenotype)
  ems$Phenotype<-gsub("Flower","Days to flowering",ems$Phenotype)
  ems$Phenotype<-gsub("Bolting","Days to bolting",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLeafNumber","Leaf number at bolting",ems$Phenotype)
  ems$Phenotype<-gsub("MaxLength","Length of largest leaf at bolting (cm)",ems$Phenotype)
  ems$Generation<-gsub("Grandprogeny","3rd",ems$Generation)
  ems$Generation<-gsub("Progeny","2nd",ems$Generation)
  ems$GrandmaternalTemp<-gsub("High","Warm",ems$GrandmaternalTemp)
  ems$GrandmaternalTemp<-gsub("Cold","Cool",ems$GrandmaternalTemp)
  
  
  ems$Phenotype<-factor(ems$Phenotype,levels=unique(ems$Phenotype))
  ems$sig<-ems$p.value<=0.05
  #ems$GrandmaternalVern
  
  ems$ID<-paste(ems$Genotype,ems$GrandmaternalVern,ems$Generation,ems$GerminationTemp,sep=":")
  Germ$ID<-paste(Germ$Genotype,Germ$GrandMaternalVern,Germ$Generation,Germ$GermTemp,sep=":")
  
  #ems$ID<-paste(ems$Genotype,ems$GrandmaternalTemp,ems$Generation,ems$Phenotype,sep=":")
  #colnames(Comb)[5]<-"Phenotype"
  #Comb$ID<-paste(Comb$Genotype,Comb$GrandmaternalTemp,Comb$Generation,Comb$Phenotype,sep=":")
  #Comb.plot<-merge(Comb,ems[,c(10,9)],by="ID",all.x=TRUE)
  #Comb.plot$GrandmaternalVern<-gsub("Unvernalized","Not vernalized",Comb.plot$GrandmaternalVern)
  
  ems2<-ems
  #ems2$Generation<-paste(ems2$Generation,"generation")
  
  
  GLMnegVern<-ggplot(ems2,aes(x=Generation,y=estimate,color=GrandmaternalTemp,group=GrandmaternalTemp,linetype=GrandmaternalTemp)) +
    geom_errorbar(aes(ymin=estimate-(-SE),ymax=estimate+(-SE),linetype=NULL),alpha=1,width=0.6) +
    #geom_ribbon(aes(ymin=estimate-SE,ymax=estimate+SE,fill=GrandmaternalVern,color=NULL),alpha=0.5) +
    geom_line(aes(x=(Generation)),size=0.75,alpha=1) + 
    geom_point(aes(shape=sig,fill=sig,size=2)) + 
    scale_shape_manual(values = c(1,16),  guide=FALSE) +
    #geom_smooth(method="gam",alpha=0.25,aes(x=as.numeric(Generation))) + 
    geom_hline(yintercept = 0,alpha=0.25) +
    guides(size=FALSE) +
    scale_linetype_manual(name=NULL,values = c("solid","dashed")) +
    #geom_text(label = lm_eqn(lm(estimate ~ Generation, decay.data)), parse = TRUE) +
    facet_grid(GerminationTemp~Genotype,scales="free") + 
    scale_color_hc(name=NULL) +  scale_fill_tableau(guide=FALSE) + 
    theme_bw(base_size = 14) +
    labs(y="Response to vernalization (Effect size of NV-V)",title = "C. Response to first-generation vernalization treatment") +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y=element_text(size=9),
          legend.position = "bottom",
          axis.title.y = element_text(size = rel(0.8), angle = 90))  +
    theme(plot.margin = margin(6, 0, 6, 0))
  #guides(fill=guide_legend(nrow=2,byrow=FALSE))
  
  pdf(paste("Germ.RxN.tweaked.pdf",sep=""),height=15,width=9)
  print(
    plot_grid(RAW + theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
              GLMnegTemp+ theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
              GLMnegVern+ theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
              nrow=3,
              #align="vh",
              #labels = c("A. Genotypic mean phenotype for each experimental treatment",
              #           "B. Response to first-generation thermal regime",
              #           "C. Response to first-generation vernalization treatment"),
              rel_heights = c(1,1,1.2),
              hjust = -1)
  )
  dev.off()
  
  library(cowplot)
  #pdf(paste(Phenos[i],".RxN.vern.pdf",sep=""),height=13,width=10)
  ToPrint<-plot_grid(RAW + theme(legend.position="none") ,GLMnegTemp +theme(legend.position="none") ,
                     GLMnegVern +theme(legend.position="none") ,
                     nrow=3,
                     #align="vh",
                     labels = c("A", "B","C"),
                     rel_heights = c(1,1,1.2),
                     hjust = -1)
  legend1 <- get_legend(
    # create some space to the left of the legend
    RAW + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend2 <- get_legend(
    # create some space to the left of the legend
    GLMnegTemp + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  legend3 <- get_legend(
    # create some space to the left of the legend
    GLMnegVern + guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  #print(plot_grid(ToPrint, legend1,legend2,legend3, ncol = 1, rel_heights = c(1, .05,.05,.05)))
  #dev.off()
}
FigMaker.germ.tweaked()
