#load the required library 
library(ggplot2)

##read in data file
#the input file was obtained using the "adegenet_X_chrom_locus.sh" script; the output .csv file from this script was annotated to add, for each sample, information on population and lineage membership
pca_X_chrom_locus<-read.csv(file="PCA_X_chrom_locus_annotated.csv",header = TRUE, sep = ,)
pca_X_chrom_locus_invasive<- pca_X_chrom_locus[pca_X_chrom_locus$Lineage == 'Invasive',]
pca_X_chrom_locus_native<- pca_X_chrom_locus[pca_X_chrom_locus$Lineage != 'Invasive',]
pca_X_chrom_locus_native$Lineage = factor(pca_X_chrom_locus_native$Lineage, levels=c("WC", "JIC", "CAB","BAH","ESM", "POR"))

#plot PC1 and PC2, samples colour-coded based on lineage
ggplot(pca_X_chrom_locus_invasive, aes(PC1, PC2))+
    geom_point(colour="lightgrey", alpha=0.6, size=2)+
    geom_point(data=pca_X_chrom_locus_native, aes(PC1, PC2, colour=Lineage))+
    theme_bw()+
    theme(axis.title=element_text(size=12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "none")+
    xlab("\\nPC1 (25.69%)")+
    ylab("PC2 (17.82%)\\n")+
    scale_colour_manual(values = c("#e09148ff","#b53e35ff","#9dc8e9ff","#3f919bff","#415b9eff","#853789ff"))

  
#plot PC1 and PC2, samples colour-coded based on lineage  
ggplot(pca_X_chrom_locus_invasive, aes(PC1, PC2,colour=as.factor(Kmeans_group)))+
    geom_point(alpha=0.6, size=2)+
    geom_point(data=pca_X_chrom_locus_native, aes(PC1, PC2, colour=as.factor(Kmeans_group)))+
    theme_bw()+
    theme(axis.title=element_text(size=12),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = "none")+
    xlab("\\nPC1 (25.69%)")+
    ylab("PC2 (17.82%)\\n")+
    scale_colour_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))


#read in local ancestry assignments for the X chromosome, for the invasive population pairs used in the genome scan analyses (these files are obtained using the "adegenet_X_chrom_locus.sh" script)
MON_VOL<-read.csv(file="population_pairs/MON_VOL.X_chrom_locus.txt",header = TRUE, sep = ,)
CIT_SAR<-read.csv(file="population_pairs/CIT_SAR.X_chrom_locus.txt",header = TRUE, sep = ,)
HIG_MIA<-read.csv(file="population_pairs/HIG_MIA.X_chrom_locus.txt",header = TRUE, sep = ,)
BRO_TIF<-read.csv(file="population_pairs/BRO_TIF.X_chrom_locus.txt",header = TRUE, sep = ,)
GLA_PAL<-read.csv(file="population_pairs/GLA_PAL.X_chrom_locus.txt",header = TRUE, sep = ,)
COL_DUV<-read.csv(file="population_pairs/COL_DUV.X_chrom_locus.txt",header = TRUE, sep = ,)
BRE_LEE<-read.csv(file="population_pairs/BRE_LEE.X_chrom_locus.txt",header = TRUE, sep = ,)
MAN_STP<-read.csv(file="population_pairs/MAN_STP.X_chrom_locus.txt",header = TRUE, sep = ,)
HAM_MARi<-read.csv(file="population_pairs/HAM_MARi.X_chrom_locus.txt",header = TRUE, sep = ,)
FLA_LAK<-read.csv(file="population_pairs/FLA_LAK.X_chrom_locus.txt",header = TRUE, sep = ,)
CMB_ORA<-read.csv(file="population_pairs/CMB_ORA.X_chrom_locus.txt",header = TRUE, sep = ,)
ALA_POL<-read.csv(file="population_pairs/ALA_POL.X_chrom_locus.txt",header = TRUE, sep = ,)
HEN_STJ<-read.csv(file="population_pairs/HEN_STJ.X_chrom_locus.txt",header = TRUE, sep = ,)
LOW_STL<-read.csv(file="population_pairs/LOW_STL.X_chrom_locus.txt",header = TRUE, sep = ,)
LEV_TAM<-read.csv(file="population_pairs/LEV_TAM.X_chrom_locus.txt",header = TRUE, sep = ,)

library(MASS)
chisq.test(x = MON_VOL$X_chr_locus_ancestry, y = MON_VOL$Population)
chisq.test(x = CIT_SAR$X_chr_locus_ancestry, y = CIT_SAR$Population)
chisq.test(x = HIG_MIA$X_chr_locus_ancestry, y = HIG_MIA$Population)
chisq.test(x = BRO_TIF$X_chr_locus_ancestry, y = BRO_TIF$Population)
chisq.test(x = GLA_PAL$X_chr_locus_ancestry, y = GLA_PAL$Population)
chisq.test(x = COL_DUV$X_chr_locus_ancestry, y = COL_DUV$Population)
chisq.test(x = BRE_LEE$X_chr_locus_ancestry, y = BRE_LEE$Population)
chisq.test(x = MAN_STP$X_chr_locus_ancestry, y = MAN_STP$Population)
chisq.test(x = HAM_MARi$X_chr_locus_ancestry, y = HAM_MARi$Population)
chisq.test(x = FLA_LAK$X_chr_locus_ancestry, y = FLA_LAK$Population)
chisq.test(x = CMB_ORA$X_chr_locus_ancestry, y = CMB_ORA$Population)
chisq.test(x = HEN_STJ$X_chr_locus_ancestry, y = HEN_STJ$Population)
chisq.test(x = LOW_STL$X_chr_locus_ancestry, y = LOW_STL$Population)
chisq.test(x = LEV_TAM$X_chr_locus_ancestry, y = LEV_TAM$Population)
chisq.test(x = ALA_POL$X_chr_locus_ancestry, y = ALA_POL$Population)
  