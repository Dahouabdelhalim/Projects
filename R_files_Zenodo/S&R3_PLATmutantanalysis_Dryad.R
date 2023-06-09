##### Scripts to create figures and recreate analysis in Burghardt et al Plant Physiology #######
# It takes as input various files included in the Dryad repository including
# Measurements of plant and nodule phenotypes of WT and NPD- NPDMutantPhenoytpes_community_DRYAD.txt
# Strain frequency files - C68_WT&NPDmutants_nodulestrainfreq.txt
# Plant benefit data from a previous single strain experiment- SingleStrain_phenotype_summary.tsv
# GWAS results files - Multiple files in the subfolder entitled 'npd_gwas_final'. 
# The _dryad file for each GWAS run contains a row and location information for each SNP variant (85933) and info on group membership
# The _output files contains a row for each linkage group (26328) and detailed test results, maf, and beta etc...

#######Info about Host mutations######:
# Wildtype background is HM340 (R108)
# trna18-1   npd2 1 gene knockout
# Gm12-7  npd2/4 2 gene knockout
# Gm8-5 npd1/2/4  3 gene knockout
# Gm26-1 npd2/4/5  3 gene knockout
# Gm20-2 npd1/2/4/5 4 gene knockout
# MPC12 npd1/2/3/4/5 5 gene knockout (called npd1-5 in the paper)

setwd("~/yourfolder/")

library("ggplot2")
library("lsmeans")
library("dplyr")
library(colorspace)
mycols<-rev(rainbow_hcl(7, start = 20, end = 280))
mycols<-mycols[c(3,1,6)]


##### Code to make Figure 2, Table S1, Table S2, Figure S2 #######


# Load in phenotype data for community experiments
phenotype = read.table('NPDMutantPhenoytpes_community_DRYAD.txt', header = T)

# Re-organize the levels and give them new labels
phenotype$Genotype = factor(phenotype$Genotype, 
                            levels=c("HM340","tRNA181","Gm127","Gm85",  "Gm261", "Gm202",  "MPC12"), 
                            labels=c("WT",   "npd2",   "npd2/4","npd2/4/5","npd1/2/4","npd1/2/4/5","npd1-5"))

Pheno_C68 = phenotype[phenotype$Block != 'Sm1021' & phenotype$Block !='Peat' & is.na(phenotype$SweightperPlant_g)<1 & phenotype$Genotype != 'npd2/4' & phenotype$Genotype != 'npd1/2/4' & phenotype$Genotype != 'npd2/4/5' & phenotype$Genotype != 'npd1/2/4/5', ]
Pheno_Sm1021 = phenotype[phenotype$Block == 'Sm1021'& phenotype$Genotype != 'npd2/4' & phenotype$Genotype != 'npd1/2/4' & phenotype$Genotype != 'npd2/4/5' & phenotype$Genotype != 'npd1/2/4/5', ]
Pheno_C68$Genotype<-factor(Pheno_C68$Genotype)
Pheno_C68$TotalweightperPlant_g <-(Pheno_C68$RweightperPlant_g + Pheno_C68$SweightperPlant_g)

####Make the graphs for Figure 1 and Figure S3 ##########

pdf(file = "Figure1.pdf",height=2.15, width=3.2)

#### Biomass per plant

model = lm( TotalweightperPlant_g ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
Biomass_per_plant<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(RweightperPlant_g + SweightperPlant_g)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab("Per plant dry biomass (mg)") +
  geom_point(data = Pheno_Sm1021, aes(x = Genotype, y = (RweightperPlant_g + SweightperPlant_g)), 
             size   =1.5, pch = 24, color="black", fill="white", position = position_nudge(x=0.15))+
  annotate("text", label = c("a","a","b"), x = c(1,2,3), y = .33, size = 3, colour = "black")   

### Nod # per plant

model = lm(NodNumPerPlant ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
NodNumPerPlant<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(NodNumPerPlant)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab("Nodule number per plant") +
  geom_point(data = Pheno_Sm1021, aes(x = Genotype, y = (NodNumPerPlant)), 
             size   =1.5, pch = 24, color="black", fill="white", position = position_nudge(x=0.15))+
  annotate("text", label = c("a","b","ab"), x = c(1,2,3), y = 45, size = 3, colour = "black")  

# Nodule Size
model = lm(AvgNodArea_mm2 ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
AvgNodArea_mm2<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(AvgNodArea_mm2)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab("Mean nodule size (mm2)") +
  geom_point(data = Pheno_Sm1021, aes(x = Genotype, y = (AvgNodArea_mm2)), 
             size   =1.5, pch = 24, color="black", fill="white", position = position_nudge(x=0.15))+
  annotate("text", label = c("a","a","b"), x = c(1,2,3), y = 2.5, size = 3, colour = "black") 

# Rhizobia per Nodule

model = lm(log10(rhizobPernod) ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
Rhizobia_per_nod<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(rhizobPernod)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab("Colony-forming rhizobia per nodule") +
  scale_y_log10()+
  annotate("text", label = c("a","b","b"), x = c(1,2,3), y = 10^6, size = 3, colour = "black")

# Rhizobia per plant

model = lm(log10(RhizobNum/NumPlants) ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
Rhizobia_per_plant<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(RhizobNum/NumPlants)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab(label="Colony-forming rhizobia per plant") +
  scale_y_log10()+
  annotate("text", label = c("a","ab","b"), x = c(1,2,3), y = 10^7.5, size = 3, colour = "black")

### % Pink

model = lm(X.Pink ~ Genotype, data = Pheno_C68)
marginal = lsmeans(model, ~ Genotype)
Per_pink_nodules<-pairs(marginal, adjust="tukey")

ggplot(data=Pheno_C68, aes(x=Genotype, y=(X.Pink)), na.rm = TRUE) + 
  geom_boxplot(width = 0.5,fill=mycols)  + theme_bw() +
  geom_jitter(position=position_dodge(0.9), size = 1.5, pch = 21, color="black", fill = "black") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8),panel.grid=element_blank())+
  ylab("Percent pink nodules") +
  geom_point(data = Pheno_Sm1021, aes(x = Genotype, y = (X.Pink)), 
             size   =1.5, pch = 24, color="black", fill="white", position = position_nudge(x=0.15))+
  annotate("text", label = c("a","a","b"), x = c(1,2,3), y = 100, size = 3, colour = "black")

dev.off()

####Save information for Table S2 ######

write.table(x=rbind(data.frame(Phenotype="NodNumPerPlant",NodNumPerPlant),
                    data.frame(Phenotype="AvgNodArea_mm2",AvgNodArea_mm2),
                    data.frame(Phenotype="Percent_pink_nodules",Per_pink_nodules),
                    data.frame(Phenotype="Biomass_per_plant",Biomass_per_plant),
                    data.frame(Phenotype="Rhizobia_per_nod",Rhizobia_per_nod),
                    data.frame(Phenotype="Rhizobia_per_plant",Rhizobia_per_plant)),file="TableS2.tsv",sep="\\t",row.names = FALSE,col.names=TRUE)

####### Aggregate data  for Table S1) ---------

NumBlocks =   aggregate(Block ~ Genotype+Inoc,data=Pheno_C68,length)
names(NumBlocks)[3] <- "Blocks"
nPlants =     aggregate(NumPlants ~ Genotype+Inoc,data=Pheno_C68,mean)
names(nPlants)[3] <- "Plants per Block (mean)"
SdnPlants =     aggregate(NumPlants ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdnPlants)[3] <- "Plants per Block (sd)"
AvgBiomass =  aggregate(RweightperPlant_g + SweightperPlant_g ~ Genotype+Inoc,data=Pheno_C68,mean) 
names(AvgBiomass)[3] <- "Biomass (mean)"
SdBiomass =   aggregate(RweightperPlant_g + SweightperPlant_g ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdBiomass)[3] <- "Biomass (sd)"
AvgNodpp =   aggregate(NodNumPerPlant ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgNodpp)[3] <- "Nodules per Plant (mean)"
SdNodpp = aggregate(NodNumPerPlant ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdNodpp)[3] <- "Nodules per Plant (sd)"
AvgNod =   aggregate(NodNum ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgNod)[3] <- "Total Nodules (mean)"
SdNod = aggregate(NodNum ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdNod)[3] <- "Total Nodules (sd)"
AvgNodmm = aggregate(AvgNodArea_mm2 ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgNodmm)[3] <- "Nodule Size (mean, mm2)"
SdNodmm= aggregate(AvgNodArea_mm2 ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdNodmm)[3] <- "Nodule Size (sd, mm2)"
AvgRhiz = aggregate(log10(rhizobPernod) ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgRhiz)[3] <- "Rhizobia per Nodule (mean)"
SdRhiz = aggregate(log10(rhizobPernod) ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdRhiz)[3] <- "Rhizobia per Nodule (sd)"
AvgRhizpp = aggregate(log10(RhizobNum/NumPlants) ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgRhizpp)[3] <- "Rhizobia per Plant (mean)"
SdRhizpp = aggregate(log10(RhizobNum/NumPlants) ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdRhizpp)[3] <- "Rhizobia per Plant (sd)"
AvgPercentPink = aggregate(X.Pink ~ Genotype+Inoc,data=Pheno_C68,mean)
names(AvgPercentPink)[3] <- "% Pink Nodules (mean)"
SdPercentPink =  aggregate(X.Pink ~ Genotype+Inoc,data=Pheno_C68,sd)
names(SdPercentPink)[3] <- "% Pink Nodules (sd)"

SummaryTable <- left_join(NumBlocks , nPlants)
SummaryTable <- left_join(SummaryTable , SdnPlants)
SummaryTable <- left_join(SummaryTable, AvgBiomass)
SummaryTable <- left_join(SummaryTable, SdBiomass)
SummaryTable <- left_join(SummaryTable, AvgNodpp)
SummaryTable <- left_join(SummaryTable, SdNodpp)
SummaryTable <- left_join(SummaryTable, AvgNod)
SummaryTable <- left_join(SummaryTable, SdNod)
SummaryTable <- left_join(SummaryTable, AvgNodmm)
SummaryTable <- left_join(SummaryTable, SdNodmm)
SummaryTable <- left_join(SummaryTable, AvgRhiz)
SummaryTable <- left_join(SummaryTable, SdRhiz)
SummaryTable <- left_join(SummaryTable, AvgRhizpp)
SummaryTable <- left_join(SummaryTable, SdRhizpp)
SummaryTable <- left_join(SummaryTable, AvgPercentPink)
SummaryTable <- left_join(SummaryTable, SdPercentPink)

write.table(SummaryTable,"TableS1.tsv",sep="\\t", row.names=FALSE,col.names=TRUE)

rm(list = ls())


######### Code for assessing strain frequency differences Figure 2,3,4 ##########

library(ggplot2)
library(vegan)
library(colorspace)
mycols<-rev(rainbow_hcl(7, start = 20, end = 280))
mycols<-mycols[c(3,1,2,4,5,7,6)]


freqsC68 = read.csv('C68_WT&NPDmutants_nodulestrainfreq.txt',
                    sep='\\t', header=TRUE, as.is=TRUE)

##### Check for correlations in strain frequency between differentiated (D) and undifferentiated fractions (U) ######
cor.test(as.numeric(freqsC68[freqsC68$pool=="HM340_D1",-1]),as.numeric(freqsC68[freqsC68$pool=="HM340_U1",-1]))
cor.test(as.numeric(freqsC68[freqsC68$pool=="HM340_D2",-1]),as.numeric(freqsC68[freqsC68$pool=="HM340_U2",-1]))
cor.test(as.numeric(freqsC68[freqsC68$pool=="HM340_D3",-1]),as.numeric(freqsC68[freqsC68$pool=="HM340_U3",-1]))
cor.test(as.numeric(freqsC68[freqsC68$pool=="MPC12_D1",-1]),as.numeric(freqsC68[freqsC68$pool=="MPC12_U1",-1]))
cor.test(as.numeric(freqsC68[freqsC68$pool=="MPC12_D2",-1]),as.numeric(freqsC68[freqsC68$pool=="MPC12_U2",-1]))
cor.test(as.numeric(freqsC68[freqsC68$pool=="MPC12_D3",-1]),as.numeric(freqsC68[freqsC68$pool=="MPC12_U3",-1]))

#### Since correlations are high within a plant only use differentiated strain frequencies

freqsC68_nod<-freqsC68[-grep(pattern = "_U",x = freqsC68$pool),]
freqsC68_nod<-freqsC68_nod[-grep(pattern = "initial",x = freqsC68_nod$pool),]
freqsC68_nod<-freqsC68_nod[order(freqsC68_nod$pool),]

freqsC68_initial<-freqsC68[grep(pattern = "initial",x = freqsC68$pool),]
initial_means = apply(freqsC68_initial[,-1], 2, mean)
hist(initial_means,breaks = 30)


##### Create Figure S2 #######

RDA_analysis<-function(freqs,treat,mycols,initial=initial_means){
  mine<-as.matrix(freqs[,-1])
  mine<-log2(mine/initial)
  mine[mine< (-8)]<- (-8)
  rel_mine<-data.frame(treat=treat,mine)
  rel_mine$treat<-factor(rel_mine$treat,levels=c("HM340","npd2","npd2/4","npd1/2/4","npd2/4/5","npd1/2/4/5","npd1/2/3/4/5"))
  rda1<-rda(rel_mine[,c(2:69)]~rel_mine$treat,rel_mine, scale=TRUE)
  scale<-1
  par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
  rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-2,4),ylim=c(-4,4),scaling=scale, xlab=paste("RDA 1 (",round(summary(rda1)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("RDA 2 (",round(summary(rda1)$cont$importance[1,2],2),"% var.)",sep=""))
  orditorp(rda1,display="species",col="black",air=.01,select = colMeans(mine)>0, pch="",scaling=scale)
  points(rda1,"wa", cex=0.8,pch=16,col=mycols[rel_mine$treat])
  ordiellipse(rda1, rel_mine$treat, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
  ordispider(rda1, rel_mine$treat, lwd=1,label =TRUE,col=paste(mycols),cex=.5)
  terms<-anova(rda1, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
  #axis<-anova(rda1, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation 
  text(3, -2,paste("Trt Adj.R^2= ",round(RsquareAdj(rda1)$adj.r.squared,3),sep=""),cex = .5)
  text(3, -2.5,paste("Host: DF= ",terms[1,1]," Var= ",round(terms[1,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""),cex = .5)
  text(3, -3,paste("Resid: DF= ",terms[2,1]," Var= ",round(terms[2,2],1),sep=""),cex = .5)
  #text(3, -3.5,paste("Axis1= ",axis[1,4]," Axis2= ",axis[2,4],sep=""),cex = .5)
}

#Note the outlier replicate of npd1/2/4/5 has been removed

pdf(file="FigureS2.pdf",height=4, width=4)
RDA_analysis(freqs=freqsC68_nod[-5,], treat=c(rep("npd2/4",3),rep("npd1/2/4/5",2),rep("npd2/4/5",3),rep("npd1/2/4",3),rep("npd2",3),rep("HM340",3),rep("npd1/2/3/4/5",3)),
             mycols=mycols)
dev.off()

#### Subset down to strain frequency data for only WT, npd2, and npd1-5 hosts which is the focus of the paper #####
freqsC68_sub<-rbind(freqsC68_nod[grep(pattern = "HM340",x = freqsC68_nod$pool),], freqsC68_nod[grep(pattern = "trna",x = freqsC68_nod$pool),],freqsC68_nod[grep(pattern = "MPC",x = freqsC68_nod$pool),])
mycols<-rev(rainbow_hcl(7, start = 20, end = 280))
mycols<-mycols[c(3,1,6)]

#### Traits for GWAS and heatmap and histogram ######
WT_freq = apply(freqsC68_sub[grep("HM340",freqsC68_sub$pool), -1], 2, median)
npd2_freq = apply(freqsC68_sub[grep("trna18-1",freqsC68_sub$pool), -1], 2, median)
npd12345_freq = apply(freqsC68_sub[grep("MPC12",freqsC68_sub$pool), -1], 2, median)

WT_fit = log2(apply(freqsC68_sub[grep("HM340",freqsC68_sub$pool), -1], 2, median)/initial_means)
npd2_fit = log2(apply(freqsC68_sub[grep("trna18-1",freqsC68_sub$pool), -1], 2, median)/initial_means)
npd12345_fit = log2(apply(freqsC68_sub[grep("MPC12",freqsC68_sub$pool), -1], 2, median)/initial_means)

traits<-data.frame(label=names(WT_freq), WT_freq,npd2_freq,npd12345_freq,WT_fit,npd2_fit,npd12345_fit)

traits$label<- as.character(traits$label)
traits$label <- gsub("X","",traits$label)
traits[traits< (-8)]<- (-8)

write.table(x = traits,quote=FALSE,row.names = FALSE, col.names = TRUE,sep="\\t",file = "NPD_GWAStraits_Dryad.tsv")

##### Make Figure 2a ######
diversity_analysis<-function(freqs, trt,title="",mycols,myylim=c(0,4),mylev=c("HM340","npd2","npd1/2/3/4/5")){
  mine<-as.matrix(freqs[,-1])
  row.names(mine)<-freqs$pool
  mydiver<-renyi(mine)
  mydiver$Treatment=trt
  mydiver$Treatment<-factor(mydiver$Treatment,levels=mylev)
  par(mfrow=c(1,1),mar=c(2, 4, 1, 1))
  boxplot(mydiver$`0`~mydiver$Treatment,col=mycols[c(1:length(unique(mydiver$Treatment)))],ylab="Species #",xlab="",main=paste(title)) # Only # of species
  boxplot(mydiver$`1`~mydiver$Treatment,col=mycols[c(1:length(unique(mydiver$Treatment)))],ylab="Shannon's Diversity",xlab="",main=paste(title),ylim=myylim) # Shannons diversity
  boxplot(mydiver$`2`~mydiver$Treatment,col=mycols[c(1:length(unique(mydiver$Treatment)))],ylab="Inverse Simpson",xlab="",main=paste(title)) # Inverse Simpsons
  boxplot(mydiver$`Inf`~mydiver$Treatment,col=mycols[c(1:length(unique(mydiver$Treatment)))],ylab="Eveness (Inf)",xlab="",main=paste(title))
}

pdf(file="Figure2a.pdf",width = 4,height=5, useDingbats=FALSE)
diversity_analysis(freqs=freqsC68_sub,trt=c(rep("WT",3),rep("npd2",3),rep("npd1-5",3)),title="",myylim=c(2.5,3.75),mycols=mycols,mylev=c("WT","npd2","npd1-5"))
dev.off()

mine<-as.matrix(freqsC68_sub[,-1])
row.names(mine)<-freqsC68_sub$pool
mydiver<-renyi(mine)
mydiver$Treatment=c(rep("WT",3),rep("npd2",3),rep("npd1-5",3))
mydiver$Treatment<-factor(mydiver$Treatment,levels=c("WT","npd2","npd1-5"))

diverstats<-rbind(anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd2"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd1-5"),]))[1,c(4:5)],
anova(lm(`0`~Treatment,data=mydiver[mydiver$Treatment %in% c("npd2","npd1-5"),]))[1,c(4:5)],

anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd2"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd1-5"),]))[1,c(4:5)],
anova(lm(`1`~Treatment,data=mydiver[mydiver$Treatment %in% c("npd2","npd1-5"),]))[1,c(4:5)],

anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd2"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd1-5"),]))[1,c(4:5)],
anova(lm(`2`~Treatment,data=mydiver[mydiver$Treatment %in% c("npd2","npd1-5"),]))[1,c(4:5)],

anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd2"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("WT","npd1-5"),]))[1,c(4:5)],
anova(lm(`Inf`~Treatment,data=mydiver[mydiver$Treatment %in% c("npd2","npd1-5"),]))[1,c(4:5)])

row.names(diverstats)<-paste(rep(c("Wt-npd2","Wt-npd1-5","npd2-npd1-5"),4), c("0","0","0","1","1","1","2","2","2","Inf","Inf","Inf"))
diverstats$contrast<- as.character(row.names(diverstats))
write.table(diverstats,file="DiversityStats_Figure2a.tsv",sep="\\t",row.names = FALSE,col.names = TRUE)

##### Histograms for Figure 2b ######################

library(reshape2)
Norms<-melt(traits[c(5:7)],value.name="Fitness",variable.name = "Env")

pdf("Figure2b.pdf",height=5, width=4)
ggplot(Norms,aes(x=Fitness,fill=Env,color=Env))+
  geom_histogram()+
  facet_grid(Env~.)+
  geom_vline(xintercept = 0,lty=2)+
  theme_bw()+
  theme(legend.position="none")+
  scale_fill_manual(values = paste(mycols[c(1,2,3)]))+
  scale_color_manual(values = c("black","black","black","black","black","black","black","black","black","black","black","black")) 
dev.off()

###### RDA analysis for Figure 2c ###############

RDA_analysis<-function(freqs,treat,mycols,initial){
  mine<-as.matrix(freqs[,-1])
  mine<-log2(mine/initial)
  mine[mine< (-8)]<- (-8)
  rel_mine<-data.frame(treat=treat,mine)
  rel_mine$treat<-factor(rel_mine$treat,levels=c("npd2","WT","npd1-5"))
  rda1<-rda(rel_mine[,c(2:69)]~rel_mine$treat,rel_mine, scale=TRUE)
  terms<-anova(rda1, step=1000, perm.max=1000, by= "terms") # if you have multiple things the terms is type 1 ANOVA
  axis<-anova(rda1, step=1000, perm.max=1000, by= "axis") # are the RDA axis explaining a sign portion of the total variation
  scale<-1
  par(mfrow=c(1,1),mar=c(3, 3, 2.1, 1))
  rdaplot <- ordiplot(rda1, display=c("wa"),cex=.5,cex.axis=.8,cex.lab=.9, tck=.02,mgp=c(1.7,.5,0),type="none",xlim=c(-3,2),ylim=c(-3,3.5),scaling=scale, 
                      xlab=paste("RDA 1 (",round(summary(rda1)$cont$importance[1,1],2),"% var.)",sep=""), ylab=paste("RDA 2 (",round(summary(rda1)$cont$importance[1,2],2),"% var.)",sep=""))
  points(rda1,"wa", cex=0.8,pch=16,col=mycols[rel_mine$treat])
  ordiellipse(rda1, rel_mine$treat, kind="se", conf=0.95, lwd=2,label =FALSE,draw = "polygon", col=paste(mycols), border=paste(mycols),lty=c(1) ,alpha=63)
  ordispider(rda1, rel_mine$treat, lwd=1,label =TRUE,col=paste(mycols),cex=.5)
  text(-2, 3.5,paste("Trt Adj.R^2= ",round(RsquareAdj(rda1)$adj.r.squared,3),sep=""),cex = .5)
  text(-2, 3,paste("Host: DF= ",terms[1,1]," Var= ",round(terms[1,2],1)," F= ",round(terms[1,3],1)," p= ",terms[1,4],sep=""),cex = .5)
  text(-2, 2.5,paste("Resid: DF= ",terms[2,1]," Var= ",round(terms[2,2],1),sep=""),cex = .5)
  text(-2, 2,paste("Axis1= ",axis[1,4]," Axis2= ",axis[2,4],sep=""),cex = .5)
}

pdf(file="Figure2c.pdf",width = 4,height=4, useDingbats=FALSE)
RDA_analysis(freqs=freqsC68_sub[c(4:6,1:3,7:9),], treat=c(rep("npd2",3),rep("WT",3),rep("npd1-5",3)),
             mycols=mycols,initial=initial_means)
dev.off() 

####### Figure 3 ######

#### Create dataframs of min, med, max #######
bulk_med = apply(freqsC68_initial[, -1], 2, median)
WT_med = apply(freqsC68_sub[grep("HM340",freqsC68_sub$pool), -1], 2, median)
npd2_med = apply(freqsC68_sub[grep("trna18-1",freqsC68_sub$pool), -1], 2, median)
npd12345_med = apply(freqsC68_sub[grep("MPC12",freqsC68_sub$pool), -1], 2, median)

bulk_max = apply(freqsC68_initial[, -1], 2,  max)
WT_max = apply(freqsC68_sub[grep("HM340",freqsC68_sub$pool), -1], 2,  max)
npd2_max = apply(freqsC68_sub[grep("trna18-1",freqsC68_sub$pool), -1], 2, max)
npd12345_max = apply(freqsC68_sub[grep("MPC12",freqsC68_sub$pool), -1], 2, max)

bulk_min = apply(freqsC68_initial[, -1], 2,  min)
WT_min = apply(freqsC68_sub[grep("HM340",freqsC68_sub$pool), -1], 2,  min)
npd2_min = apply(freqsC68_sub[grep("trna18-1",freqsC68_sub$pool), -1], 2, min)
npd12345_min = apply(freqsC68_sub[grep("MPC12",freqsC68_sub$pool), -1], 2, min)

trt<-rbind(cbind(bulk_med,bulk_min,bulk_max,"initial"),cbind(WT_med,WT_min,WT_max,"WT"),cbind(npd2_med,npd2_min,npd2_max,"npd2"),cbind(npd12345_med,npd12345_min,npd12345_max,"npd1-5"))
trt_all<-data.frame(Strains=row.names(trt),Med=as.numeric(trt[,1]),Min=as.numeric(trt[,2]),Max=as.numeric(trt[,3]),Env=trt[,4])
trt_all$Strains<-as.character(gsub(x=trt_all$Strains, pattern = "X",replacement = ""))
trt_all$Strains<- factor(trt_all$Strains,levels=paste(trt_all$Strains[order(trt_all[trt_all$Env=="WT",]$Med)]))
#trt_all$Strains<- factor(trt_all$Strains,levels=paste(trt_all$Strains[order(trt_all[trt_all$Env=="npd1/2/3/4/5",]$Med)]))
trt_all$Env<- factor(trt_all$Env,levels=c("npd1-5","npd2","WT"))

base.plot <- function(data,mycols=mycols, max=0.5) {
  p <- ggplot(data, aes(x=Med, y=Strains,col=Env,shape=Env))
  p <- p + theme_bw()
  p <- p + theme(legend.position="top")
  p <- p + geom_point(size=2, alpha=0.8)
  p <-p+ geom_segment(mapping=aes(x=data$Min,xend=data$Max,y=as.numeric(data$Strains),yend=as.numeric(data$Strains)))
  p <- p + geom_point(mapping=aes(x=data$Med,y=as.numeric(data$Strains)),size=2.5, alpha=0.8)
  #p <- p + scale_color_manual(values = mycols, guide=guide_legend(ncol=5, title=NULL))
  p <- p + scale_color_manual(values = mycols)
  p <- p+ xlim(c(0,max))
  p <- p+ geom_vline(xintercept = .015,lty=2)
  p <- p + xlab("Strain Frequency") + ylab("")
  return(p)
}

pdf(file = "Figure3.pdf",height = 10 ,width=5,family = "Helvetica",useDingbats = FALSE)
base.plot(data=trt_all[trt_all$Env %in% c("WT","npd2","npd1-5"),],mycols = mycols[c(3,2,1)])
base.plot(data=trt_all[trt_all$Env %in% c("WT","npd2","npd1-5"),],mycols = mycols[c(3,2,1)],max=0.08)
dev.off()

#### Heatmap Figure 4a ######
library (gplots)
fitness<-traits[,paste(c("label","WT_fit","npd2_fit","npd12345_fit"))]
colnames(fitness) <- gsub("_fit","",colnames(fitness))
pal = colorRampPalette(c('#2166ac', '#92c5de', '#f7f7f7', '#b2182b'))(length(seq(-7.9,4,0.1)))
m = structure(as.matrix(fitness[, -1]),
                dimnames=list(fitness[, 'label'],
                              names(fitness)[-1]))
m<-m[order(m[,1],decreasing = TRUE),]

pdf(file="Figure4left.pdf",width = 5,height=8, useDingbats=FALSE)
heatmap.2(m, Colv = NA, Rowv = NA,dendrogram = 'none',trace = "none",col=pal,symm = FALSE,symbreaks = FALSE,density.info="none",key.title=NA,key.xlab=NA,keysize = 2,breaks=seq(-8,4,0.1),lwid = c(4, 5),lhei=c(1, 10),cexCol = 1)
dev.off()

##### Heatmap Figure 4b ######
# This benefit data was originally published in the Dryad respository for the original Select and Reseq paper ()
  
benefit<-read.csv(file = "SingleStrain_phenotype_summary.tsv",sep = '\\t')
benefit<-benefit[benefit$plant_genotype=="R108",]
benefit<-benefit[,c('strain','weight')]
benefit$weight2<-benefit$weight  #### just create to columns so a heatmap can be
benefit<-merge(benefit,fitness,by.x = 'strain',by.y='label',all.y = TRUE)
benefit<-benefit[match(rownames(m), benefit$strain), ]

pal_b = colorRampPalette(c('goldenrod','grey80','olivedrab','darkgreen'))(length(seq(-.59,.6, 0.01)))
n = structure(as.matrix(benefit[, c("weight","weight2")]),
              dimnames=list(benefit[, 'strain']))

pdf(file="Figure4right.pdf",width = 5,height=8, useDingbats=FALSE)
heatmap.2(n, Colv = NA, Rowv = NA,dendrogram = 'none',trace = "none",col=pal_b,symm = FALSE,symbreaks = FALSE,density.info="none",key.title=NA,key.xlab=NA,keysize = 2,lwid = c(4, 5),lhei=c(1, 10),cexCol = 1,scale="column")
dev.off()

rm(list = ls())

#### GWAS ######
#### These GWAS results come from the same pipeline/methods used in Epstein et al 2018 mSphere #####

GWAS_WT = read.csv('npd_gwas_final/gwas_results_gemma_dryad_WT_fit.tsv',
                    sep='\\t', header=TRUE, as.is=TRUE)
GWAS_npd2 = read.csv('npd_gwas_final/gwas_results_gemma_dryad_npd2_fit.tsv',
                   sep='\\t', header=TRUE, as.is=TRUE)
GWAS_npd12345 = read.csv('npd_gwas_final/gwas_results_gemma_dryad_npd12345_fit.tsv',
                   sep='\\t', header=TRUE, as.is=TRUE)
GWAS_WTvnpd2 = read.csv('npd_gwas_final/gwas_results_gemma_dryad_npd2_minus_WT_fit.tsv',
                        sep='\\t', header=TRUE, as.is=TRUE)
GWAS_WTv5npd = read.csv('npd_gwas_final/gwas_results_gemma_dryad_npd12345_minus_WT_fit.tsv',
                        sep='\\t', header=TRUE, as.is=TRUE)
GWAS_5npdvnpd2 = read.csv('npd_gwas_final/gwas_results_gemma_dryad_npd2_minus_npd12345_fit.tsv',
                          sep='\\t', header=TRUE, as.is=TRUE)

##### Compile rankings for each group in other contrasts to make the supplementary tables of candidate genes####

##### Table S3: WT vs npd2 gene comparison #####
Top_WTvnpd2<-GWAS_WTvnpd2[GWAS_WTvnpd2$rank< 10,]
  WTvnpd2top<- data.frame(group=GWAS_WTv5npd[match(unique(Top_WTvnpd2$rs),GWAS_WTv5npd$rs),]$rs,
                          Location=GWAS_WTv5npd[match(unique(Top_WTvnpd2$rs),GWAS_WTv5npd$rs),]$chr, 
                          WTrank=GWAS_WT[match(unique(Top_WTvnpd2$rs),GWAS_WT$rs),]$rank,
                          WTvnpd2rank=GWAS_WTvnpd2[match(unique(Top_WTvnpd2$rs),GWAS_WTvnpd2$rs),]$rank,
                          WTv5npdrank=GWAS_WTv5npd[match(unique(Top_WTvnpd2$rs),GWAS_WTv5npd$rs),]$rank,
                          npd2v5npdrank= GWAS_5npdvnpd2[match(unique(Top_WTvnpd2$rs),GWAS_5npdvnpd2$rs),]$rank)
  
  mygroup<-tapply(Top_WTvnpd2$locus_annotation,INDEX=Top_WTvnpd2$rs,FUN = function(x) unique (x)[1:15][!is.na((unique (x)[1:15]))])
  mygroup<-lapply(mygroup,FUN = 'paste',collapse=";")
  mygroup<-mygroup[match(WTvnpd2top$group, names(mygroup))]
  WTvnpd2top<-data.frame(WTvnpd2top,Genes.15max=unlist(mygroup))

  Betas <- read.csv('npd_gwas_final/npd2_minus_WT_fit.association_output.tsv',
                    sep='\\t', header=TRUE, as.is=TRUE)
  Betas<-Betas[WTvnpd2top$group,]
  WTvnpd2top<-merge(WTvnpd2top,Betas, by.x="group", by.y="rs")
  
  WTvnpd2top<-WTvnpd2top[order(WTvnpd2top$WTvnpd2rank),]
  WTvnpd2top<-WTvnpd2top[,c(1,13,14,3:6,19,7)]
  WTvnpd2top[,4:7] <- (WTvnpd2top[,4:7]+1)
  write.table(WTvnpd2top, file='npd_gwas_final/TableS3_WTvsnpd2.tsv', sep='\\t', row.names=FALSE,col.names=TRUE)
  
  
##### Table S4: WT vs 5 gene comparison #####
Top_WTv5npd<-GWAS_WTv5npd[GWAS_WTv5npd$rank< 10,]
  
  WTv5npdtop<- data.frame(group=GWAS_WTv5npd[match(unique(Top_WTv5npd$rs),GWAS_WTv5npd$rs),]$rs,
                          Location=GWAS_WTv5npd[match(unique(Top_WTv5npd$rs),GWAS_WTv5npd$rs),]$chr, 
                          WTrank=GWAS_WT[match(unique(Top_WTv5npd$rs),GWAS_WT$rs),]$rank,
                          WTvnpd2rank=GWAS_WTvnpd2[match(unique(Top_WTv5npd$rs),GWAS_WTvnpd2$rs),]$rank,
                          WTv5npdrank=GWAS_WTv5npd[match(unique(Top_WTv5npd$rs),GWAS_WTv5npd$rs),]$rank,
                          npd2v5npdrank= GWAS_5npdvnpd2[match(unique(Top_WTv5npd$rs),GWAS_5npdvnpd2$rs),]$rank)
  
  mygroup<-tapply(Top_WTv5npd$locus_annotation,INDEX=Top_WTv5npd$rs,FUN = function(x) unique (x)[1:15][!is.na((unique (x)[1:15]))])
  mygroup<-lapply(mygroup,FUN = 'paste',collapse=";")
  mygroup<-mygroup[match(WTv5npdtop$group, names(mygroup))]
  WTv5npdtop<-data.frame(WTv5npdtop,Genes.15max=unlist(mygroup))
  
  Betas <- read.csv('npd_gwas_final/npd12345_minus_WT_fit.association_output.tsv',
                    sep='\\t', header=TRUE, as.is=TRUE)
  Betas<-Betas[WTv5npdtop$group,]
  
  WTv5npdtop<-merge(WTv5npdtop,Betas, by.x="group", by.y="rs")
  WTv5npdtop<-WTv5npdtop[,c(1,13,14,3:6,19,7)]
  WTv5npdtop<-WTv5npdtop[order(WTv5npdtop$p_lrt),]
  WTv5npdtop[,4:7] <- (WTv5npdtop[,4:7]+1)
  write.table(WTv5npdtop, file='npd_gwas_final/TableS4_WTvsnpd12345.tsv', sep='\\t', row.names=FALSE,col.names=TRUE)
  
##### Table S5: npd2 vs 5gene comparison #####
Top_npd2v5npd<-GWAS_5npdvnpd2[GWAS_5npdvnpd2$rank< 10,]
  
  npd2v5npdtop<- data.frame(group=GWAS_5npdvnpd2[match(unique(Top_npd2v5npd$rs),GWAS_5npdvnpd2$rs),]$rs,
                          Location=GWAS_5npdvnpd2[match(unique(Top_npd2v5npd$rs),GWAS_5npdvnpd2$rs),]$chr, 
                          WTrank=GWAS_WT[match(unique(Top_npd2v5npd$rs),GWAS_WT$rs),]$rank,
                          WTvnpd2rank=GWAS_WTvnpd2[match(unique(Top_npd2v5npd$rs),GWAS_WTvnpd2$rs),]$rank,
                          WTv5npdrank=GWAS_WTv5npd[match(unique(Top_npd2v5npd$rs),GWAS_WTv5npd$rs),]$rank,
                          npd2v5npdrank= GWAS_5npdvnpd2[match(unique(Top_npd2v5npd$rs),GWAS_5npdvnpd2$rs),]$rank)
  
  mygroup<-tapply(Top_npd2v5npd$locus_annotation,INDEX=Top_npd2v5npd$rs,FUN = function(x) unique (x)[1:15][!is.na((unique (x)[1:15]))])
  mygroup<-lapply(mygroup,FUN = 'paste',collapse=";")
  mygroup<-mygroup[match(npd2v5npdtop$group, names(mygroup))]
  npd2v5npdtop<-data.frame(npd2v5npdtop,Genes.15max=unlist(mygroup))
  
  Betas <- read.csv('npd_gwas_final/npd2_minus_npd12345_fit.association_output.tsv',
                    sep='\\t', header=TRUE, as.is=TRUE)
  Betas<-Betas[npd2v5npdtop$group,]
  npd2v5npdtop<-merge(npd2v5npdtop,Betas, by.x="group", by.y="rs")
  
  npd2v5npdtop<-npd2v5npdtop[order(npd2v5npdtop$npd2v5npdrank),]
  npd2v5npdtop<-npd2v5npdtop[,c(1,13,14,3:6,19,7)]
  npd2v5npdtop[,4:7] <- (npd2v5npdtop[,4:7]+1)
  write.table(npd2v5npdtop, file='npd_gwas_final/TableS5_npd2vsnpd12345.tsv', sep='\\t', row.names=FALSE,col.names=TRUE)
  

#### Figure S5 #######

library ("ade4")
  
loc<-c(rep("PsymB",14483),
       rep("Chr",1969),
       rep("PsymA", 9876))


mine<-matrix(ncol = 3)

for(n in 1:1000){
  new<-table(sample(loc,50))/50
  if(length(new)<3) {
    new=c(0,new)}
  mine<-rbind(mine,new)
}

mine<-mine[-1,]

all<-data.frame(Freq=c(median(mine[,1]), .92,	.90,	.36,	.58,median(mine[,2]), .06,	.08,	.24,	.36,median(mine[,3]), .02,	.02,	.40,	.06),
                Loc=rep(c("Chr","PsymA","PsymB"),each=5),
                Contrast=rep(c("Random","WT","WTvnpd2","WTvnpd12345","npd2vnpd12345"),3))

all<-data.frame(Chr=c(.92,	.90,	.36,	.58),
                PsymA=c(.06,	.08,	.24,	.36),
                PsymB=c ( .02,	.02,	.40,	.06))

randtriangleplot<-triangle.plot(mine,scale=FALSE)
mycols<-c("#39BCBE","orange","darkred","darkblue")                
row.names(all)<-c("WT","WTvnpd2","WTvnpd12345","npd2vnpd12345")

pdf("FigureS5.pdf",useDingbats = FALSE,height=5,width=5)
wtriangleplot <- triangle.plot(all,label = row.names(all),clabel=.75,show.position=FALSE,scale=FALSE)
points(randtriangleplot, col = rgb(red = .5, green = 0.5, blue = 0.5, alpha = 0.05), cex = 1,pch=c(20))
points(wtriangleplot, col = mycols, cex = 1,pch=c(19,19,19,19,19))
dev.off()

