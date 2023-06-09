###########################################################################
########      STATISTICAL ANALYSES on niche inferences  (part 2)    #######
###########################################################################

###  Load libraries  ###

library(gtools)
library(vegan)
library(adespatial)
library(boot)
library(resample)
library(raster)
library(ecospat)
library(reshape2)
library(multcompView)
library(plyr)
library(ggplot2)
library(car)
library(devtools)
library(ggpubr) #devtools::install_github("kassambara/ggpubr")
library(stringi)
library(stringr)



###########################################
###  Set working directories and paths  ###

mainDir = "PATH/to/workingDirectory/"
dataDir = "PATH/to/data/"
outDir = "PATH/to/outputs/"
alyDir = "PATH/to/analyses/"
setwd(mainDir) 

# load R objects generated for niche modeling (refer to Rscript "EcoGen-1-Niche_modeling.RDat")
load("PATH/to/EcoGen-1-Niche_modeling.RData")

source("PATH/to/generate_TukeyHSDlabels_ggplot.R")



##################################################################
###  --------------------------------------------------------  ###
###  ------  Compute niche CENTROID and niche BREADTH  ------  ###
###  --------------------------------------------------------  ###
##################################################################


#matrix to store niche caracteristics
niche = matrix(nrow=12, ncol=4, dimnames=list(c(mcdu,"pCr","pCy","pGe","pTr"),c("pos1","breadth1","pos2","breadth2")))	


for(i in c(mcdu,"pCr","pCy","pGe","pTr")) { #for each chosen species
  
  R=100 #dimension of gridding of env space
  z = get(paste0("z.",i))
  
  
  c1=NULL;c2=NULL;b1=NULL;b2=NULL
  
  for(bt in 1:1000){
    
    #row index of 100 random pixel weighted by density along PC1 and PC2 
    c<-sample(1:(R*R),100,prob=values(z$z.uncor)) #rs = 600 random pixel
    # coordinates of the pixels along PC1 and PC2
    y=(c%/%R)+1 ; x=c%%R 
    
    # scores of random pixels on PC1 & PC2
    PC1.sim<-z$x[x]
    PC2.sim<-z$y[y]
    
    # niche centroids
    c1 = c(c1, median(PC1.sim))
    c2 = c(c2, median(PC2.sim))
    # niche breadths
    b1 = c(b1, var(PC1.sim))
    b2 = c(b2, var(PC2.sim))
    
    
  }
  
  niche[i,1]<-mean(c1) # niche position on CP1
  niche[i,2]<-mean(b1) # niche breadth on CP1
  niche[i,3]<-mean(c2) # niche position on CP2
  niche[i,4]<-mean(b2) # niche breadth on CP2
  
  assign(paste0("c1.",i), c1)
  assign(paste0("c2.",i), c2)
  assign(paste0("b1.",i), b1)
  assign(paste0("b2.",i), b2)
  
}

write.table(niche, file=paste0(alyDir,"ecoNiche_centroids.txt"), quote=F,sep="\\t")



###################################################################
###        Summarize niche centroids/breadths bootstraps        ### 
###################################################################

df1=NULL
df2=NULL


for(i in c("Ge","Tr","Cy","Cr")){
  
  if(i=="Ge") { p1="Um"; p2="Co"; exp="pGe" 
  } else if (i=="Tr") { p1="Um"; p2="Ca"; exp="pTr" 
  } else if (i=="Cy") { p1="Ta"; p2="Ca"; exp="pCy" 
  } else if (i=="Cr") { p1="Ta"; p2="Co"; exp="pCr" }
  
  # PC1
  df1 = rbind(df1, 
              cbind(i,i,  get(paste0("c1.",i)),get(paste0("b1.",i))),
              cbind(i,exp,get(paste0("c1.p",i)),get(paste0("b1.p",i))),
              cbind(i,p1, get(paste0("c1.",p1)),get(paste0("b1.",p1))),
              cbind(i,p2, get(paste0("c1.",p2)),get(paste0("b1.",p2))) )
  
  # PC2
  df2 = rbind(df2, 
              cbind(i,i,  get(paste0("c2.",i)),get(paste0("b2.",i))),
              cbind(i,exp,get(paste0("c2.p",i)),get(paste0("b2.p",i))),
              cbind(i,p1, get(paste0("c2.",p1)),get(paste0("b2.",p1))),
              cbind(i,p2, get(paste0("c2.",p2)),get(paste0("b2.",p2))) )
  
}


df1 = data.frame(df1, stringsAsFactors=F) ; df2 = data.frame(df2, stringsAsFactors=F) 
colnames(df1) = c("px","cat","cPC1","bPC1") ; colnames(df2) = c("px","cat","cPC2","bPC2")
df1$cPC1 = as.numeric(df1$cPC1) ; df2$cPC2 = as.numeric(df2$cPC2)
df1$bPC1 = as.numeric(df1$bPC1) ; df2$bPC2 = as.numeric(df2$bPC2)
df1$px = as.factor(df1$px) ; df2$px = as.factor(df2$px)
df1$cat = as.factor(df1$cat) ; df2$cat = as.factor(df2$cat)

df1$cat <- factor(df1$cat, levels = c(di,paste0("p",tetra),tetra))
df2$cat <- factor(df2$cat, levels = c(di,paste0("p",tetra),tetra))


df1.summary = aggregate(df1, by=list(df1$sp), FUN=summary)
df2.summary = aggregate(df2, by=list(df2$sp), FUN=summary)
write.table(df1.summary, file=paste0(alyDir,"summary_stats_PC1-centroids.txt"), quote=F, sep="\\t", row.names=F)
write.table(df2.summary, file=paste0(alyDir,"summary_stats_PC2-centroids.txt"), quote=F, sep="\\t", row.names=F)
#write.table(df1.summary, file=paste0(alyDir,"summary_stats_PC1-breadths.txt"), quote=F, sep="\\t", row.names=F)
#write.table(df2.summary, file=paste0(alyDir,"summary_stats_PC2-breadths.txt"), quote=F, sep="\\t", row.names=F)




######################################################################
#### -------------  one-way ANOVA + Tukey HSD test  ------------- ####
######################################################################

df1$group = paste0(df1$px,"_",df1$cat); df1$group = as.factor(df1$group)
df2$group = paste0(df2$px,"_",df2$cat); df2$group = as.factor(df2$group)


## ----  Niche centroid: PC1  ---- ##

model1c = lm(cPC1~group, data=df1)
ANOVA1c = aov(model1c)
TUKEY1c = TukeyHSD(x=ANOVA1c, 'group', conf.level=0.95)
labels1c = generate_label_df(TUKEY1c, 'group', df1)
labels1c = labels1c[,1:2] ; names(labels1c) = c("group", "Letters")
yvalue1c = aggregate(.~group, data=df1, mean)
final1c = merge(labels1c, yvalue1c) ; final1c$group = sub(".*_", "", final1c$group)
write.table(final1c, paste0(alyDir, "TukeyHSD_PC1-labels-centroids.txt"), quote=F, sep="\\t", row.names=F)

## ----  Niche centroid: PC2  ---- ##

model2c = lm(cPC2~group, data=df2)
ANOVA2c = aov(model2c)
TUKEY2c = TukeyHSD(x=ANOVA2c, 'group', conf.level=0.95)
labels2c = generate_label_df(TUKEY2c, 'group', df2)
labels2c = labels2c[,1:2] ; names(labels2c) = c("group", "Letters")
yvalue2c = aggregate(.~group, data=df2, mean)
final2c = merge(labels2c, yvalue2c) ; final2c$group = sub(".*_", "", final2c$group)
write.table(final2c, paste0(alyDir, "TukeyHSD_PC2-labels-centroids.txt"), quote=F, sep="\\t", row.names=F)



## ----  Niche breadths: PC1  ---- ##

model1b = lm(bPC1~group, data=df1)
ANOVA1b = aov(model1b)
TUKEY1b = TukeyHSD(x=ANOVA1b, 'group', conf.level=0.95)
labels1b = generate_label_df(TUKEY1b, 'group', df1)
labels1b = labels1b[,1:2] ; names(labels1b) = c("group", "Letters")
yvalue1b = aggregate(.~group, data=df1, mean)
final1b = merge(labels1b, yvalue1b) ; final1b$group = sub(".*_", "", final1b$group)
write.table(final1b, paste0(alyDir, "TukeyHSD_PC1-labels-breadths.txt"), quote=F, sep="\\t", row.names=F)

## ----  Niche breadths: PC2  ---- ##

model2b = lm(bPC2~group, data=df2)
ANOVA2b = aov(model2b)
TUKEY2b = TukeyHSD(x=ANOVA2b, 'group', conf.level=0.95)
labels2b = generate_label_df(TUKEY2b,'group',df2)
labels2b = labels2b[,1:2] ; names(labels2b) = c("group", "Letters")
yvalue2b = aggregate(.~group, data=df2, mean)
final2b = merge(labels2b, yvalue2b) ; final2b$group = sub(".*_", "", final2b$group)
write.table(final2b, paste0(alyDir,"TukeyHSD_PC2-labels-breadths.txt"), quote=F, sep="\\t", row.names=F)

# Note: Tukey HSD significance labels are manually added on the ggplot figure based on the above tables.



#################################################################
#### -----------------  Generate boxplots  ----------------- ####
#################################################################

df1$group = NULL
df2$group = NULL


df1.m = melt(df1)
df2.m = melt(df2)
df1.m$cat <- factor(df1.m$cat, levels = c("Ta","Um","Ca","Co",paste0("p",tetra),tetra))
df2.m$cat <- factor(df2.m$cat, levels = c("Ta","Um","Ca","Co",paste0("p",tetra),tetra))


df.m = rbind(df1.m, df2.m)
df.m$px <- factor(df.m$px, levels=tetra)
df.m$variable <- factor(df.m$variable, levels = c("cPC1","cPC2","bPC1","bPC2"))
df.m$cat <- factor(df.m$cat, levels = c("Ta","Um","Ca","Co","pGe","pTr","pCy","pCr","Ge","Tr","Cy","Cr"))
df.m$colcode <- c("P1","P1","P2","P2","ePx","ePx","ePx","ePx","oPx","oPx","oPx","oPx")[df.m$cat]
df.m$colcode <- factor(df.m$colcode, levels = c("P1","P2","ePx","oPx"))

cols=c("gold","lightblue","darkgreen","darkorange2")
colcode = cols[df.m$colcode]
dodge=position_dodge(width = 0.4)
labels = c( expression(paste("Diploid progenitor A" ["o"], "A" ["o"])), 
            expression(paste("Diploid progenitor B" ["o"], "B" ["o"])),
            expression(paste("Expected allopolyploid A" ["e"], "A"["e"], "B" ["e"], "B" ["e"])), 
            expression(paste("Observed allopolyploid A" ["o"], "A"["o"], "B" ["o"], "B" ["o"])) )



pdf(file=paste0(alyDir,"EcoNiches_boots_boxplots.pdf"), width=22, height=18) # defaults (width=7, height=7)

ggplot(data=df.m, aes(x=px, y=value, col=colcode, fill=colcode)) + 
  geom_boxplot(position = dodge, width=.2) +
  facet_wrap(variable~., ncol=2, strip.position=c("top"), scales="free_y",
             labeller=labeller(variable=c("cPC1"="Niche centroids PC1","cPC2"="Niche centroids PC2","bPC1"="Niche breadths PC1","bPC2"="Niche breadths PC2"))) +
  theme(panel.spacing = unit(1.5, "lines"), 
        axis.title.x = element_blank(), strip.text.x=element_text(size=20, face="bold"), 
        axis.text.x = element_text(size=20, angle=45, face="italic", vjust=.55), 
        axis.text.y = element_text(size=20), axis.title.y = element_text(size=20),
        legend.text = element_text(size=20), legend.title = element_blank(), legend.text.align = 0, legend.position="right", legend.key.size=unit(1.5, 'cm')) +
  ylab("Values") +
  scale_fill_manual(name="Species", values=cols, labels=labels) +
  scale_color_manual(name="Species", values=cols, labels=labels) +
  scale_x_discrete(breaks=c("Ge","Tr","Cy","Cr"), 
                   labels=str_wrap(c("Ae. geniculata (UUMM)","Ae. triuncialis (UUCC)","Ae. cylindrica (DDCC)","Ae. crassa (DDMM)"), width=15) ) 

dev.off()





#################################################################
###  -------------------------------------------------------  ###
###          One-way ANOVA test on each ecovariable           ###
###  -------------------------------------------------------  ###
#################################################################

aegc.m = melt(aegc[,c("sp",name.maps)],id.vars="sp")

#build expected range values of ecovar for Px based on combined values of their progenitors??
tmp = aegc.m[which(aegc.m$sp %in% c("Um","Co")),] ; tmp$sp = "pGe" ; tm = tmp
tmp = aegc.m[which(aegc.m$sp %in% c("Um","Ca")),] ; tmp$sp = "pTr" ; tm = rbind(tm,tmp)
tmp = aegc.m[which(aegc.m$sp %in% c("Ta","Ca")),] ; tmp$sp = "pCy" ; tm = rbind(tm,tmp)
tmp = aegc.m[which(aegc.m$sp %in% c("Ta","Co")),] ; tmp$sp = "pCr" ; tm = rbind(tm,tmp)
aegc.m = rbind(aegc.m,tm)


pdf(paste0(alyDir,"biovar_ANOVAtest_pPx.pdf"))

for(var in name.maps){
  
  ff = aegc.m[which(aegc.m$variable==var),]
  
  model = lm(value~sp, data=ff)
  ANOVA = aov(model)
  TUKEY = TukeyHSD(x=ANOVA, 'sp', conf.level=0.95)
  labels = generate_label_df(TUKEY,'sp',ff)
  labels = labels[,1:2]
  names(labels) = c("sp", "Letters")
  yvalue = aggregate(.~sp, data=ff, mean)
  final = merge(labels,yvalue) 
  
  write.table(final,paste0(alyDir,"biovar_TukeyHSD_",var,".txt"),quote=F,sep="\\t",row.names=F)
  
  
  plotclim = ggplot(ff, aes(x=sp, y=value)) + geom_boxplot() + ggtitle(var) + geom_text(data=final, aes(x=sp,y=value,label=Letters), vjust=-3.5,hjust=-1.0) + theme(plot.title = element_text(hjust = 0.5)) + scale_x_discrete(limits=c(di,"Ge","pGe","Tr","pTr","Cy","pCy","Cr","pCr"))
  print(plotclim)
  
}

dev.off()



aegc.m.summary = aggregate(aegc.m, by=list(aegc.m$sp,aegc.m$variable), FUN=summary)
write.table(aegc.m.summary, file=paste0(alyDir,"ecovar_summary_pPx.txt"), sep="\\t", quote=F, row.names=F)
write.table(final,paste0(alyDir,"biovar_TukeyHSD_",var,"_pPx.txt"),quote=F,sep="\\t",row.names=F)







