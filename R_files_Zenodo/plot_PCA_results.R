#load the required library 
library(ggplot2)

##read in data files for all samples
#the input file was obtained using the "adegenet_PCA.R" script; the output .csv file from this script was annotated to add information on range and lineage
pca_all<-read.csv(file="PCA_all_annotated.csv",header = TRUE, sep = ,)

#partition invasive- and native-range samples and order the native range lineages 
pca_invasive<- pca_all[pca_all$Range == 'Invasive',]
pca_native<- pca_all[pca_all$Range != 'Invasive',]
pca_native$Lineage = factor(pca_native$Lineage, levels=c("W_Cuba", "JIC", "CAB","BAH","ESM", "POR"))

#PLOT PC1 PC2 for all samples
ggplot(pca_invasive, aes(PC1, PC2))+
  geom_point(colour="lightgrey", alpha=0.6, size=2)+
  geom_point(data=pca_native, aes(PC1, PC2, colour=Lineage))+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  #percent variance explained by each axis is calculated using the "adegenet_PCA.R" script
  xlab("\\nPC1 (12.06%)")+
  ylab("PC2 (4.18%)\\n")+
  scale_colour_manual(values = c("#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))


##read in data files for Western Cuba samples
#the input file was obtained using the "adegenet_PCA.R" script; the output .csv file from this script was annotated to add population information
pca_W_Cuba<-read.csv(file="PCA_W_Cuba_annotated.csv",header = TRUE, sep = ,)

#PLOT PC1 PC2 for Western Cuba samples
ggplot(pca_W_Cuba, aes(PC1, PC2))+
  geom_point(colour="#e09148ff", size=6, aes(shape=Population), stroke=1)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())+
  xlab("\\nPC1 (5.17%)")+
  ylab("PC2 (4.10%)\\n")+
  scale_shape_manual(values = c(15,3,16,17,8))
