setwd("~/yourfolder")

library(ggplot2)
library(factoextra)
library(FactoMineR)
library(vegan)
library(reshape2)
library(GGally)
library(cocor)

mycols<-c("#27008B", "#B5E300", "#7B0085", "#E9C100", "#AC177F", "#E5872C", "#E06698")

# Load strain frequencies from the first and second sequencing runs.
# Some notes on sample labels and how they differ from the paper.  
##### HM101 is used interchangaby with A17 
##### HM340 is used interchangably with R108
##### bulk is used interchangably with initial

freqs0 = read.csv('freq1.tsv',
                  sep='\\t', header=TRUE, as.is=TRUE)
freqs1 = read.csv('freq2.tsv',
                  sep='\\t', header=TRUE, as.is=TRUE)


if(!all(colnames(freqs0) == colnames(freqs1)))
  stop('Column names don\\'t match')

freqs = rbind(freqs0, freqs1)

bulk_pools = c('101Bulk_Epi1', '101Bulk_Epi2', '101Bulk_Fal1', '101Bulk_Fal2',
               'Bulk_Epi_05', 'Bulk_Epi_06', 'Bulk_Epi_09', 'Bulk_Epi_10')

HM101_pools= c('HM101_1R','HM101_2R','HM101_3R','HM101_4R','HM101_5R','HM101_6R','HM101_7R','HM101_8R')
HM340_pools= c('HM340_1R','HM340_2R','HM340_3R','HM340_4R','HM340_5R','HM340_6R','HM340_7R','HM340_8R',"HM340_XR")
soil_HM340_pools= c('HM340NP_09','HM340NP_10','HM340NP_11','HM340NP_12','HM340NP_13','HM340NP_14','HM340NP_15','HM340NP_16')
soil_HM101_pools= c('HM101NP_09','HM101NP_10','HM101NP_11','HM101NP_12','HM101NP_13','HM101NP_14','HM101NP_15','HM101NP_16')
HM101_HM340_pools= c('HM101D_09','HM101D_10','HM101D_11','HM101D_12','HM101D_13','HM101D-14','HM101D_16')
HM340_HM340_pools= c('HM340S_09','HM340S_10','HM340S_11','HM340S_15','HM340S_13','HM340S_14','HM340S_16')

bulk_means = apply(freqs[freqs$pool %in% bulk_pools, -1], 2, mean)

#### Frequencies for Diversity Analyses ####
mine<-rbind(freqs[freqs$pool %in% bulk_pools,-1],
            freqs[freqs$pool %in% HM101_pools, -1],
            freqs[freqs$pool %in% HM340_pools, -1],
            freqs[freqs$pool %in% soil_HM101_pools, -1],
            freqs[freqs$pool %in% soil_HM340_pools, -1],
            freqs[freqs$pool %in% HM101_HM340_pools, -1],
            freqs[freqs$pool %in% HM340_HM340_pools, -1])

rownames(mine)<-c(freqs[freqs$pool %in% bulk_pools,1],
                  freqs[freqs$pool %in% HM101_pools, 1],
                  freqs[freqs$pool %in% HM340_pools, 1],
                  freqs[freqs$pool %in% soil_HM101_pools,1],
                  freqs[freqs$pool %in% soil_HM340_pools, 1],
                  freqs[freqs$pool %in% HM101_HM340_pools, 1],
                  freqs[freqs$pool %in% HM340_HM340_pools, 1])


# To rerun the analyis with a different detection threshold remove the below hashtag
#mine[mine<.0005]<-0 adding a .0005 threshold for presence
mydiver<-renyi(mine)
mydiver$Treatment=c(rep("Initial",8), rep ("A17", 8), rep("R108", 9), rep ("Soil_A17", 8), rep("Soil_R108", 8), rep ("A17_R108", 7), rep ("R108_R108", 7))
mydiver$Treatment<-factor(mydiver$Treatment,levels = c("Initial","A17", "R108", "Soil_A17","Soil_R108","A17_R108","R108_R108"))

#### Figure S3 ######

pdf(file = "LegacyDiversity_FigS3.pdf", height=5, width=7,useDingbats = FALSE)
diver<-melt(data = mydiver)
#p <- ggplot(data=diver, aes(x=as.numeric(variable), y=value,ymin=0,ymax=5, fill=Treatment,color=Treatment))
p <- ggplot(data=diver, aes(x=variable, y=value,ymin=.9,ymax=5, fill=Treatment,color=Treatment))
p <- p + geom_point(position = position_jitterdodge())
p <- p + geom_smooth(aes(x=as.numeric(variable)))
p <- p + scale_fill_manual(values = paste(mycols))
p <- p + scale_color_manual(values = paste(mycols))
p <- p + xlab(label = "Scale Parameter")
p <- p + ylab(label = "Strain Diversity")
p<- p + theme_minimal()
p
dev.off()

#### Figure 2b ######

pdf(file = "LegacyDiversityBoxplots.pdf", height=9, width=4)

par(mfrow=c(4,1),mar=c(2, 4, 1, 1))
boxplot(mydiver$`0.25`~mydiver$Treatment,col=mycols[c(1:7)],ylab="Exponent 0 (Species #)",xlab="") # Only # of species
boxplot(mydiver$`1`~mydiver$Treatment,col=mycols[c(1:7)],ylab="Exponent 1 (Shannons)",xlab="") # Shannons diversity
boxplot(mydiver$`2`~mydiver$Treatment,col=mycols[c(1:7)],ylab="Exponent 2 (1/Simpson)",xlab="") # Inverse Simpsons
boxplot(mydiver$`Inf`~mydiver$Treatment,col=mycols[c(1:7)],ylab="Exponent Inf",xlab="")
dev.off()

par(mfrow=c(1,1),mar=c(2, 4, 1, 1))

# Statistical tests of diversity differences

Test_results<-cbind(rbind(anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("Soil_R108","Soil_A17"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108","A17"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","A17_R108"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","Soil_R108"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("A17_R108","Soil_R108"),]))[1,c(4:5)]),

rbind(anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("Soil_R108","Soil_A17"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108","A17"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","A17_R108"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","Soil_R108"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("A17_R108","Soil_R108"),]))[1,c(4:5)]),

rbind(anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("Soil_R108","Soil_A17"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108","A17"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","A17_R108"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","Soil_R108"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("A17_R108","Soil_R108"),]))[1,c(4:5)]),

rbind(anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("Soil_R108","Soil_A17"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108","A17"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","A17_R108"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("R108_R108","Soil_R108"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("A17_R108","Soil_R108"),]))[1,c(4:5)]))

row.names(Test_results)<-c('Soil_R108 vs Soil_A17','R108 vs A17','R108_R108 vs A17_R108','R108_R108 vs Soil_R108','A17_R108 vs Soil_R108')
Test_results<-Test_results[c(2,1,4,5,3),]
write.table(Test_results,file = "CommunityStatistics.tsv",sep="\\t",row.names = TRUE)

## Normalize, then take the means

HM101_norm = apply(freqs[freqs$pool %in% HM101_pools, -1], 1,
                   function(x) log2(x / bulk_means))

HM340_norm = apply(freqs[freqs$pool %in% HM340_pools, -1], 1,
                   function(x) log2(x / bulk_means))

soil_HM101_norm = apply(freqs[freqs$pool %in% soil_HM101_pools, -1], 1,
                        function(x) log2(x / bulk_means))

soil_HM340_norm = apply(freqs[freqs$pool %in% soil_HM340_pools, -1], 1,
                        function(x) log2(x / bulk_means))

HM101_HM340_norm = apply(freqs[freqs$pool %in% HM101_HM340_pools, -1], 1,
                         function(x) log2(x / bulk_means))

HM340_HM340_norm = apply(freqs[freqs$pool %in% HM340_HM340_pools, -1], 1,
                         function(x) log2(x / bulk_means))


### Replace infinite values with -8  #
HM101_norm[HM101_norm< (-8)]<- (-8)
HM340_norm[HM340_norm< (-8)]<- (-8)
soil_HM101_norm[soil_HM101_norm< (-8)]<- (-8)
soil_HM340_norm[soil_HM340_norm< (-8)]<- (-8)
HM101_HM340_norm[HM101_HM340_norm< (-8)]<- (-8)
HM340_HM340_norm[HM340_HM340_norm< (-8)]<- (-8)

## Take the mean of the normalized fitness values for each treatment
NormMeans<-data.frame(HM101_norm = apply(HM101_norm,MARGIN = 1,FUN = "mean"),
                      HM340_norm = apply(HM340_norm,MARGIN = 1,FUN = "mean"),
                      soil_HM101_norm = apply(soil_HM101_norm,MARGIN = 1,FUN = "mean"),
                      soil_HM340_norm = apply(soil_HM340_norm,MARGIN = 1,FUN = "mean"),
                      HM101_HM340_norm = apply(HM101_HM340_norm,MARGIN = 1,FUN = "mean"),
                      HM340_HM340_norm = apply(HM340_HM340_norm,MARGIN = 1,FUN = "mean"))

RDA<-rbind(data.frame(Treatment="A17",Host="A17",Cohort="First",Prior="None",t(HM101_norm)),
           data.frame(Treatment="R108",Host="R108",Cohort="First",Prior="None",t(HM340_norm)),
           data.frame(Treatment="soil_A17",Host="A17",Cohort="Second",Prior="Soil",t(soil_HM101_norm)),
           data.frame(Treatment="soil_R108",Host="R108",Cohort="Second",Prior="Soil",t(soil_HM340_norm)),
           data.frame(Treatment="A17_R108",Host="R108",Cohort="Second",Prior="Host",t(HM101_HM340_norm)),
           data.frame(Treatment="R108_R108",Host="R108",Cohort="Second",Prior="Host",t(HM340_HM340_norm)))

options(contrasts = c("contr.sum","contr.poly"))
theme_set(theme_bw(base_family="Helvetica"))

RDA$Treatment<-factor(RDA$Treatment,levels=c("A17","R108","soil_A17", "soil_R108","A17_R108","R108_R108"))

#Create the colors for the graphs
mycols<-data.frame(Treatment=as.character(levels(RDA$Treatment)),color=mycols[-1])

#Subset the data for the four RDA analysis listed in Table1
Growth<-RDA[RDA$Treatment %in% c("A17","soil_R108","soil_A17","R108"),] # contrast 1
SecondGen<-RDA[RDA$Treatment %in% c("A17_R108","soil_R108","soil_A17","R108_R108"),] # contrast 2
soil<-RDA[RDA$Treatment %in% c("R108_R108","soil_R108"),] # contrast 3
plant<-RDA[RDA$Treatment %in% c("A17_R108","R108_R108"),] #contrast 4

# Run an RDA on each of the datasets constraining by treatment
rda_Growth<-rda(Growth[,c(5:105)]~Growth$Host+Growth$Cohort,Growth, scale=TRUE)
rda_SecondGen<-rda(SecondGen[,c(5:105)]~SecondGen$Host+SecondGen$Prior,SecondGen, scale=TRUE)
rda_soil<-rda(soil[,c(5:105)]~soil$Treatment,data = soil, scale=TRUE)
rda_plant<-rda(plant[,c(5:105)]~plant$Treatment,data=plant, scale=TRUE)

# Graph the Treatments and the important statistics. 
#P values may slightly vary from those in manuscript as the tests are based on permutations
#### The four graphs are Figure 2d, 3a, and S5a&b #####

pdf("RDA_Legacy_Contrasts.pdf",height=5, width=5, useDingbats = FALSE)

scale<-1
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rda_Growth, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-2,2),ylim=c(-3,3),scaling=scale, xlab=paste("RDA 1 (",round(summary(rda_Growth)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("RDA 2 (",round(summary(rda_Growth)$cont$importance[1,2],2),"% var.)",sep=""))
points(rda_Growth,"wa", cex=0.8,pch=16, col=paste(mycols$color[as.numeric(Growth$Treatment)]) )
ordiellipse(rda_Growth, Growth$Treatment, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols$color), border=paste(mycols$color),lty=c(1) ,alpha=63)
ordispider(rda_Growth, Growth$Treatment, lwd=1,label =FALSE,col=paste(mycols$color))
terms<-anova(rda_Growth, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
axis<-anova(rda_Growth, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation 
text(-1.5, 3,paste("Trt Adj.R^2= ",round(RsquareAdj(rda_Growth)$adj.r.squared,3),sep=""))
text(-1.5, 2.5,paste("Host: DF= ",terms[1,1]," Var= ",round(terms[1,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""))
text(-1.5, 2,paste("Growth: DF= ",terms[2,1]," Var= ",round(terms[2,2],1)," F= ",round(terms[2,3],1)," p= ",terms[2,4],sep=""))
text(-1.5, 1.5,paste("Resid: DF= ",terms[3,1]," Var= ",round(terms[3,2],1),sep=""))

scale<-1
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rda_SecondGen, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-3,2),ylim=c(-2.5,4.5),scaling=scale, xlab=paste("RDA 1 (",round(summary(rda_SecondGen)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("RDA 2 (",round(summary(rda_SecondGen)$cont$importance[1,2],2),"% var.)",sep=""))
points(rda_SecondGen,"wa", cex=0.8,pch=16, col=paste(mycols$color[as.numeric(SecondGen$Treatment)]) )
ordiellipse(rda_SecondGen, SecondGen$Treatment, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols$color), border=paste(mycols$color),lty=c(1) ,alpha=63)
ordispider(rda_SecondGen, SecondGen$Treatment, lwd=1,label =FALSE,col=paste(mycols$color))
terms<-anova(rda_SecondGen, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
axis<-anova(rda_SecondGen, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation 
text(-2, 3.5,paste("Trt Adj.R^2= ",round(RsquareAdj(rda_SecondGen)$adj.r.squared,3),sep=""))
text(-2, 3,paste("Host: DF= ",terms[1,1]," Var= ",round(terms[1,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""))
text(-2, 2.5,paste("Prior: DF= ",terms[2,1]," Var= ",round(terms[2,2],1)," F= ",round(terms[2,3],1)," p= ",terms[2,4],sep=""))
text(-2, 2.0,paste("Resid: DF= ",terms[3,1]," Var= ",round(terms[3,2],1),sep=""))

scale<-1
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rda_soil, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-4,2.5),ylim=c(-4,4.5),scaling=scale,, xlab=paste("RDA 1 (",round(summary(rda_soil)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("PCA 1 (",round(summary(rda_soil)$cont$importance[1,2],2),"% var.)",sep=""))
points(rda_soil,"wa", cex=0.8,pch=16, col=paste(mycols$color[as.numeric(soil$Treatment)]) )
ordiellipse(rda_soil, soil$Treatment, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols$color), border=paste(mycols$color),lty=c(1) ,alpha=63)
ordispider(rda_soil, soil$Treatment, lwd=1,label =FALSE,col=paste(mycols$color))
terms<-anova(rda_soil, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
axis<-anova(rda_soil, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation 
text(-3, -3,paste("Trt Adj.R^2= ",round(RsquareAdj(rda_soil)$adj.r.squared,3),sep=""))
text(-3, -3.5,paste("DF= ",terms[1,1],"/",terms[2,1]," Var= ",round(terms[1,2],1),"/",round(terms[2,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""))

scale<-1
par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
rdaplot <- ordiplot(rda_plant, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-4,2),ylim=c(-3.5,5),scaling=scale,, xlab=paste("RDA 1 (",round(summary(rda_plant)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("PCA 1 (",round(summary(rda_plant)$cont$importance[1,2],2),"% var.)",sep=""))
points(rda_plant,"wa", cex=0.8,pch=16, col=paste(mycols$color[as.numeric(plant$Treatment)]) )
ordiellipse(rda_plant, plant$Treatment, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols$color), border=paste(mycols$color),lty=c(1) ,alpha=63)
ordispider(rda_plant, plant$Treatment, lwd=1,label =FALSE,col=paste(mycols$color))
terms<-anova(rda_plant, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
axis<-anova(rda_plant, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation 
text(-4, 4,paste("Trt Adj.R^2= ",round(RsquareAdj(rda_plant)$adj.r.squared,3),sep=""))
text(-4, 3.5,paste("DF= ",terms[1,1],"/",terms[2,1]," Var= ",round(terms[1,2],1),"/",round(terms[2,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""))
dev.off()

#### Create Figure 2a  #######
Norms<-melt(as.matrix(NormMeans),variable.name="Env",value.name="Fold_change")
colnames(Norms)<-c("strain","Env","Fold_change")
Norms$Env<-factor(Norms$Env,levels = c("HM101_norm" , "HM340_norm", "soil_HM101_norm","soil_HM340_norm","HM101_HM340_norm",
                                      "HM340_HM340_norm"), labels=c("A17", "R108", "Soil_A17","Soil_R108","A17_R108","R108_R108"))

pdf("LegacyStrainFoldchange.pdf",height=8, width=3)
ggplot(Norms,aes(x=Fold_change,fill=Env,color=Env))+
  geom_histogram()+
  facet_grid(Env~.)+
  geom_vline(xintercept = 0,lty=2)+
  theme_bw()+
  theme(legend.position="none")+
  scale_fill_manual(values = paste(mycols$color))+
  scale_color_manual(values = c("black","black","black","black","black","black","black","black","black","black","black","black")) 
dev.off()

#### Summary metrics of fitness distributions ####
fitnessmetrics<-data.frame(metric=c("mean","var","median","n>0","n>mean"),rbind(apply(NormMeans,2,function(x) round(mean(x),2)),
apply(NormMeans,2,function(x) round(var(x),2)),
apply(NormMeans,2,function(x) round(median(x),2)),
apply(NormMeans,2, function(x) length(x[x>0])),
apply(NormMeans,2, function(x) length(x[x>mean(x)]))))

colnames(fitnessmetrics)<-c("metric", "A17_norm", "R108_norm","soil_A17_norm","soil_R108_norm" 
                            ,"A17_R108_norm","R108_R108_norm")

write.table(fitnessmetrics,file="FitnessDistributionMetrics.tsv",sep="\\t",row.names = FALSE,col.names=TRUE)

#### Figure S4 & S1b #####
pdf("PCAgraphsscaling95CI.pdf",width=5.5,height = 4.5,fonts = "Helvetica")

myPCA3<-PCA(SecondGen[,-1:-4], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(myPCA3, geom.ind = "point", col.ind = SecondGen$Treatment, habillage = SecondGen$Treatment, pointshape = 19,
             palette = paste(mycols[mycols$Treatment %in% c("soil_A17","soil_R108","A17_R108","R108_R108"),]$color),
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.95,
             legend.title = "Groups",mean.point=FALSE)

myPCA3<-PCA(Growth[,-1:-4], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(myPCA3, geom.ind = "point", col.ind = Growth$Treatment,habillage = Growth$Treatment, pointshape = 19, 
             palette = paste(mycols[mycols$Treatment %in% c("A17","R108","soil_A17","soil_R108"),]$color),
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.95,
             legend.title = "Groups",mean.point=FALSE)

myPCA3<-PCA(soil[,-1:-4], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(myPCA3, geom.ind = "point", col.ind = soil$Treatment,habillage = soil$Treatment, pointshape = 19, 
             palette = paste(mycols[mycols$Treatment %in% c("soil_R108","R108_R108"),]$color),
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.95,
             legend.title = "Groups",mean.point=FALSE)

myPCA3<-PCA(plant[,-1:-4], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(myPCA3, geom.ind = "point", col.ind = plant$Treatment, habillage = plant$Treatment, pointshape = 19, 
             palette = paste(mycols[mycols$Treatment %in% c("A17_R108","R108_R108"),]$color),
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.95,
             legend.title = "Groups",mean.point=FALSE)

myPCA3<-PCA(RDA[,-1:-4], scale.unit = TRUE, graph = FALSE)
fviz_pca_ind(myPCA3, geom.ind = "point", col.ind = RDA$Treatment, habillage = RDA$Treatment, pointshape = 19, 
             palette = paste(mycols[mycols$Treatment %in% c("A17","R108","soil_A17","soil_R108","A17_R108","R108_R108"),]$color),
             addEllipses = TRUE, ellipse.type = "confidence",ellipse.level=0.95,
             legend.title = "Groups",mean.point=FALSE)
dev.off()

#Reorder columns for next graph
NormMeans<-NormMeans[,c(1,3,2,4,6,5)]
colnames(NormMeans)<-c("A17legacy","soil_A17","R108legacy","soil_R108","R108_R108","A17_R108")

#### Figure S1a  ######
pdf("LegacyStrainCorrelations.pdf",height=6, width=6,useDingbats=FALSE)
p<-ggpairs(NormMeans,
           upper= list(continuous="blank"),
           lower=list(continuous="smooth"),
           diag=list(continuous=wrap("densityDiag")),
           columnLabels=c("A17legacy","soil_A17","R108legacy","soil_R108","R108_R108","A17_R108"))+
  theme_bw() + theme(panel.grid.minor = element_blank(),strip.background = element_rect(fill="white"))

mines<-matrix(nrow = 6,ncol=6)
mines[2,1]<-round(cor.test(NormMeans$A17legacy,NormMeans$soil_A17)$estimate,2)
mines[3,1]<-round(cor.test(NormMeans$A17legacy,NormMeans$R108legacy)$estimate,2)
mines[4,1]<-round(cor.test(NormMeans$A17legacy,NormMeans$soil_R108)$estimate,2)
mines[5,1]<-round(cor.test(NormMeans$A17legacy,NormMeans$R108_R108)$estimate,2)
mines[6,1]<-round(cor.test(NormMeans$A17legacy,NormMeans$A17_R108)$estimate,2)
mines[3,2]<-round(cor.test(NormMeans$soil_A17,NormMeans$R108legacy)$estimate,2)
mines[4,2]<-round(cor.test(NormMeans$soil_A17,NormMeans$soil_R108)$estimate,2)
mines[5,2]<-round(cor.test(NormMeans$soil_A17,NormMeans$R108_R108)$estimate,2)
mines[6,2]<-round(cor.test(NormMeans$soil_A17,NormMeans$A17_R108)$estimate,2)
mines[4,3]<-round(cor.test(NormMeans$R108legacy,NormMeans$soil_R108)$estimate,2)
mines[5,3]<-round(cor.test(NormMeans$R108legacy,NormMeans$R108_R108)$estimate,2)
mines[6,3]<-round(cor.test(NormMeans$R108legacy,NormMeans$A17_R108)$estimate,2)
mines[5,4]<-round(cor.test(NormMeans$soil_R108,NormMeans$R108_R108)$estimate,2)
mines[6,4]<-round(cor.test(NormMeans$soil_R108,NormMeans$A17_R108)$estimate,2)
mines[6,5]<-round(cor.test(NormMeans$A17_R108,NormMeans$R108_R108)$estimate,2)


for(i in 1:p$nrow) {
  for(j in 1:p$ncol){
    p[i,j] <- p[i,j] + scale_x_continuous(limits=c(-8,5))
    if (i!=j) {p[i,j] <- p[i,j] +scale_y_continuous(limits=c(-8,5)) }
    if (i!=j) {p[i,j] <- p[i,j]+geom_text(x=4,y=-6.75, label=mines[i,j],size=3)}
    
  }
}

for(i in 1:length(mycols$color)) {
  p[i,i]<-p[i,i]+geom_histogram(fill=paste(mycols$color[c(1,3,2,4,6,5)][i]),color="black")
}    

p

dev.off()


#### Data for Table S3 #######
##### To take away the extremely low values
#NormMeans[NormMeans < (-5)]<-(-5)

cocor(~soil_A17+A17legacy|R108legacy+A17legacy, data=NormMeans)
cocor(~R108legacy+soil_R108|R108legacy+A17legacy, data=NormMeans )
cocor(~A17legacy+soil_A17|soil_A17+soil_R108, data=NormMeans )
cocor(~soil_R108+R108legacy|soil_R108+soil_A17, data=NormMeans)
cocor(~R108legacy+R108_R108|R108legacy+soil_R108, data=NormMeans,alternative = "greater")
cocor(~A17legacy+A17_R108|A17legacy+soil_R108, data=NormMeans,alternative = "greater")



#### Calculating benefits ######
#### This Data file can also be found in the Dryad repository for Burghardt et al. 2018 PNAS Select and Resequence paper
benefit<-read.csv(file = "SingleStrain_phenotype_summary.tsv",sep = '\\t')
benefit$strainID<-paste("X",benefit$strain,sep="")

#Compile benefit data with strain fitness and freq data
A17benefit<-benefit[benefit$plant_genotype=="A17",][benefit[benefit$plant_genotype=="A17",]$strainID %in% rownames(NormMeans),]
R108benefit<-benefit[benefit$plant_genotype=="R108",][benefit[benefit$plant_genotype=="R108",]$strainID %in% rownames(NormMeans),]
R108<-data.frame(strainID=rownames(NormMeans),R108legacy=NormMeans$R108legacy,soil_R108=NormMeans$soil_R108,R108_R108=NormMeans$R108_R108,A17_R108=NormMeans$A17_R108)
A17<-data.frame(strainID=rownames(NormMeans),A17legacy=NormMeans$A17legacy,soil_A17=NormMeans$soil_A17)
A17benefit$plantweightpernod<-A17benefit$weight_raw/A17benefit$nodule_raw
R108benefit$plantweightpernod<-R108benefit$weight_raw/R108benefit$nodule_raw
R108all<-merge(x=R108benefit,y=R108,by="strainID")
A17all<-merge(x=A17benefit,y=A17,by="strainID")

freq<-data.frame(t(mine[-1:-8,]))
freq<-data.frame(t(mine))
freq$strainID<-rownames(freq)
R108freq<-merge(x=R108benefit,y=freq,by="strainID")
A17freq<-merge(x=A17benefit,y=freq,by="strainID")

### Calculate predicted plant benefit scale from 0 (worst strain) to 1 best possible strain
#Predsize<-data.frame(R108=(colSums(R108freq[,-1:-12]*R108freq$weight,na.rm = TRUE)-min(R108freq$weight,na.rm=TRUE))/(max(R108freq$weight,na.rm=TRUE)-min(R108freq$weight,na.rm=TRUE)), 
#                         A17=(colSums(A17freq[,-1:-12]*A17freq$weight,na.rm = TRUE)-min(A17freq$weight,na.rm=TRUE))/(max(A17freq$weight,na.rm=TRUE)-min(A17freq$weight,na.rm=TRUE)))
Predsize<-data.frame(R108=(colSums(R108freq[,-1:-12]*R108freq$weight,na.rm = TRUE)/colSums(R108freq[,-1:-12],na.rm = TRUE)-min(R108freq$weight,na.rm=TRUE))/(max(R108freq$weight,na.rm=TRUE)-min(R108freq$weight,na.rm=TRUE)),
                     A17=(colSums(A17freq[,-1:-12]*A17freq$weight,na.rm = TRUE)/colSums(A17freq[,-1:-12],na.rm = TRUE)-min(A17freq$weight,na.rm=TRUE))/(max(A17freq$weight,na.rm=TRUE)-min(A17freq$weight,na.rm=TRUE)))


Predsize$Treatment=c(rep ("Initial", 8),rep ("A17", 8), rep("R108", 9), rep ("Soil_A17", 8), rep("Soil_R108", 8), rep ("A17_R108", 7), rep ("R108_R108", 7))
Predsize$Treatment<-factor(Predsize$Treatment,levels = c("Initial","A17", "R108", "Soil_A17","Soil_R108","A17_R108","R108_R108"))

mycols<-c("#27008B", "#B5E300", "#7B0085", "#E9C100", "#AC177F", "#E5872C", "#E06698")

#boxplot(Predsize$R108~Predsize$Treatment,col=mycols,main="R108 single strain norm plant weight",ylim=c(0,1))
#boxplot(Predsize$A17~Predsize$Treatment,col=mycols, main="A17 single strain norm plant weight")

All<-data.frame(Size=c(Predsize[Predsize$Treatment=="Initial"|Predsize$Treatment=="A17"|Predsize$Treatment=="Soil_A17",]$A17,
Predsize[Predsize$Treatment=="Initial"|Predsize$Treatment=="R108"|Predsize$Treatment=="Soil_R108"|Predsize$Treatment=="R108_R108"|Predsize$Treatment=="A17_R108",]$R108),
Treatment=c(rep("Initial-A17", 8),rep ("A17", 8), rep ("Soil_A17", 8),rep("Initial-R108", 8),rep("R108", 9), rep("Soil_R108", 8), rep ("A17_R108", 7), rep ("R108_R108", 7)))
All$Treatment<-factor(All$Treatment,levels = c("Initial-A17","A17", "Soil_A17","Initial-R108","R108", "Soil_R108","A17_R108","R108_R108"))

pdf("PredictedPlantBenefitNormRel.pdf",width = 8,height=6,useDingbats = FALSE)
boxplot(Predsize$R108~Predsize$Treatment,col=mycols, main="R108 single strain norm plant weight",xlab="Treatment",ylab="Benefit predicted from single strain exp")
boxplot(Predsize$A17~Predsize$Treatment,col=mycols, main="A17 single strain norm plant weight",xlab="Treatment",ylab="Benefit predicted from single strain exp")
boxplot(All$Size~All$Treatment,col=mycols[c(1,2,4,1,3,5,6,7)],main="Genotype specific benefit pairings",xlab="Treatment",ylab="Benefit predicted from single strain exp")
abline(v=3.5,lty=2)
boxplot(All$Size~All$Treatment,col=mycols[c(1,2,4,1,3,5,6,7)],main="Genotype specific benefit pairings",xlab="Treatment",ylab="Benefit predicted from single strain exp",ylim=c(0,1))
abline(v=3.5,lty=2)
dev.off()

### Differences between sequencing batches####

plot(as.numeric(mine[1,]),as.numeric(mine[5,]))
points(as.numeric(mine[2,]),as.numeric(mine[6,]),col="red")
points(as.numeric(mine[3,]),as.numeric(mine[7,]),col="orange")
points(as.numeric(mine[4,]),as.numeric(mine[8,]),col="blue")

plot(as.numeric(colMeans(mine[1:4,])),as.numeric(colMeans(mine[5:8,])))
cor(as.numeric(colMeans(mine[1:4,])),as.numeric(colMeans(mine[5:8,])))
