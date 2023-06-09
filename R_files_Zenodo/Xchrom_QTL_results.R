#load the required library 
library(ggplot2)

#read in the results file for the X chromosome limb length QTL
Xchrom_QTL_all<-read.csv(file="Xchrom_QTL_results.csv",header=TRUE,sep=',')

#remove the one heterozygote genotype at this locus, which is a genotyping error 
#this is because only males were included in these analyses; males are hemizygous at the X chromosome; as expected, no other heterozygotes occurred among the 560 samples)
Xchrom_QTL<-Xchrom_QTL_all[Xchrom_QTL_all$X7_83341332 != 'HET',]

#make separate datasets for the two sample groups, "hybridization limited" and "hybridization common"
Xchrom_QTL_hybridization_limited<-Xchrom_QTL[Xchrom_QTL$Hybrid_sample_group == 'hybridization_limited',]
Xchrom_QTL_hybridization_common<-Xchrom_QTL[Xchrom_QTL$Hybrid_sample_group == 'hybridization_common',]

#make separate datasets for sample with Western Cuba ancestry at the 18 Mbp region of high differentiation on the X chromosome
Xchrom_QTL_WC<-Xchrom_QTL[Xchrom_QTL$Kmeans_Xchrom_locus_ancestry == '2',]
Xchrom_QTL_WC_hybridization_limited<-Xchrom_QTL_WC[Xchrom_QTL_WC$Hybrid_sample_group == 'hybridization_limited',]
Xchrom_QTL_WC_hybridization_common<-Xchrom_QTL_WC[Xchrom_QTL_WC$Hybrid_sample_group == 'hybridization_common',]

#violin plots for relative hindlimb length of each genotype class at the X chromosome lead GWAS SNP
#all samples are considered, irrespective of inferred ancestry at the 18Mb locus (as per Fig. 3C)
#"hybridization limited" sample group
ggplot(data = Xchrom_QTL_hybridization_limited,
       aes(x = X7_83341332, 
           y = Hind_Limb_residuals))+
  geom_violin()+
  geom_point(aes(x = X7_83341332, y=mean_HL_fig3C), size=4)+
  geom_errorbar(aes(ymin = mean_HL_fig3C-SD_HL_fig3C, ymax = mean_HL_fig3C+SD_HL_fig3C), width = 0, size=1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
    ylim(-0.045,0.045)

#test for effect of genotype, get P value and estimate r2
lm_HL_hybridization_limited<-lm(Xchrom_QTL_hybridization_limited$Hind_Limb_residuals~Xchrom_QTL_hybridization_limited$X7_83341332)
summary(lm_HL_hybridization_limited)

#"hybridization common" sample group
ggplot(data = Xchrom_QTL_hybridization_common,
       aes(x = X7_83341332, 
           y = Hind_Limb_residuals))+
  geom_violin()+
  geom_point(aes(x = X7_83341332, y=mean_HL_fig3C), size=4)+
  geom_errorbar(aes(ymin = mean_HL_fig3C-SD_HL_fig3C, ymax = mean_HL_fig3C+SD_HL_fig3C), width = 0, size=1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
    ylim(-0.045,0.045)

#test for effect of genotype, get P value, and estimate r2
lm_HL_hybridization_common<-lm(Xchrom_QTL_hybridization_common$Hind_Limb_residuals~Xchrom_QTL_hybridization_common$X7_83341332)
summary(lm_HL_hybridization_common)


#the same plots as above, but only for samples inferred to have Western Cuba ancestry at the 18Mb locus (as per Fig. S26)
#"hybridization limited" sample group
ggplot(data = Xchrom_QTL_WC_hybridization_limited,
       aes(x = X7_83341332, 
           y = Hind_Limb_residuals))+
  geom_violin()+
  geom_point(aes(x = X7_83341332, y=mean_HL_figS26), size=4)+
  geom_errorbar(aes(ymin = mean_HL_figS26-SD_HL_figS26, ymax = mean_HL_figS26+SD_HL_figS26), width = 0, size=1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  ylim(-0.045,0.045)

#test for effect of genotype, get P value, and estimate r2
lm_HL_hybridization_limited_WC<-lm(Xchrom_QTL_WC_hybridization_limited$Hind_Limb_residuals~Xchrom_QTL_WC_hybridization_limited$X7_83341332)
summary(lm_HL_hybridization_limited_WC)

#"hybridization common" sample group
ggplot(data = Xchrom_QTL_WC_hybridization_common,
       aes(x = X7_83341332, 
           y = Hind_Limb_residuals))+
  geom_violin()+
  geom_point(aes(x = X7_83341332, y=mean_HL_figS26), size=4)+
  geom_errorbar(aes(ymin = mean_HL_figS26-SD_HL_figS26, ymax = mean_HL_figS26+SD_HL_figS26), width = 0, size=1)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())+
  ylim(-0.045,0.045)

#test for effect of genotype, get P value, and estimate r2
lm_HL_hybridization_common_WC<-lm(Xchrom_QTL_WC_hybridization_common$Hind_Limb_residuals~Xchrom_QTL_WC_hybridization_common$X7_83341332)
summary(lm_HL_hybridization_common_WC)

#plot PCA results for the X chromosome 18Mb locus, with genotype at the lead GWAS SNP (as per Fig. S26A,B)
#the one heterozygote genotype is used as well for plotting the PCA
#first all samples are considered, irrespective of inferred ancestry at the 18Mb locus (Fig. S26A)

#re-order local ancestry assignments at the 18Mb locus (Western Cuba = 2, Central Cuba = 1, Eastern Cuba = 3)
Xchrom_QTL_all$Kmeans_Xchrom_locus_ancestry = factor(Xchrom_QTL_all$Kmeans_Xchrom_locus_ancestry, levels=c("2", "1", "3"))

ggplot(Xchrom_QTL_all, aes(PC1_Xchrom_locus, PC2_Xchrom_locus,colour=X7_83341332,shape=as.factor(Kmeans_Xchrom_locus_ancestry)))+
  geom_point(size=3)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 10, colour = "black"))+
  xlab("\\nPC1 (29.89%)")+
  ylab("PC2 (18.87%)\\n")+
  scale_colour_manual(values = c("#c060a6ff","#3d57a8ff","#FC4E07"))

#also plot PCA results only for samples of Western Cuba ancestry at the 18Mb locus (Fig. S26B)
Xchrom_QTL_WC$Kmeans_Xchrom_locus_ancestry = factor(Xchrom_QTL_WC$Kmeans_Xchrom_locus_ancestry, levels=c("2", "1", "3"))

ggplot(Xchrom_QTL_WC, aes(PC1_WC_Xchrom_locus, PC2_WC_Xchrom_locus,colour=X7_83341332))+
  geom_point(size=3)+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 10, colour = "black"))+
  xlab("\\nPC1 (29.89%)")+
  ylab("PC2 (18.87%)\\n")+
  scale_colour_manual(values = c("#3d57a8ff","#FC4E07"))


###PVEs (r2) and P values of genotype at SNP 7_83341332 for each trait; all PVE values were recorded in a separate .csv after these analyses, and used for plotting below
#example here given only for body size (SVL), for both groups
#all samples analysis first (as shown in Fig. 3D)
lm_SVL_hybridization_limited_all_samples<-lm(Xchrom_QTL_hybridization_limited$Log_SVL~Xchrom_QTL_hybridization_limited$X7_83341332)
summary(lm_SVL_hybridization_limited_all_samples)

lm_SVL_hybridization_common_all_samples<-lm(Xchrom_QTL_hybridization_common$Log_SVL~Xchrom_QTL_hybridization_common$X7_83341332)
summary(lm_SVL_hybridization_common_all_samples)

#adjust P values using Bonferroni after all traits have been analysed
pvals<-c(0.206,0.777,0.04524,0.8373,0.6413,0.4038,0.5086,0.005637,0.2806,0.9807,0.6284,0.002151,0.1227,0.1541,0.05826,0.5047,0.8561,0.543,0.002077,2.76E-05,0.001105,0.0001072,0.0001706,5.52E-05,1.99E-05,4.03E-06)
p.adjust(pvals,method ="bonferroni")


#only samples of Western Cuba local X chromosome ancestry (as shown in Fig. S26D)
lm_SVL_hybridization_limited_WC_samples<-lm(Xchrom_QTL_WC_hybridization_limited$Log_SVL~Xchrom_QTL_WC_hybridization_limited$X7_83341332)
summary(lm_SVL_hybridization_limited_WC_samples)

lm_SVL_hybridization_common_WC_samples<-lm(Xchrom_QTL_WC_hybridization_common$Log_SVL~Xchrom_QTL_WC_hybridization_common$X7_83341332)
summary(lm_SVL_hybridization_common_WC_samples)

pvals<-c(0.805,0.08608,0.6516,0.6845,0.9457,0.4073,0.4068,0.07063,0.2052,0.9795,0.7455,0.06573,0.3044,0.1802,0.03138,0.1023,0.7201,0.4048,0.002249,3.86E-05,0.0005974,9.55E-05,0.0003651,1.19E-04,6.51E-05,1.03E-05)
p.adjust(pvals,method ="bonferroni")


#read in PVE streadsheet, and plot results by trait, hybridization group ("hybridization limited" and "hybridization common"), and dataset (all samples as well as WC X chromosome only)
PVE_results<-read.csv(file="PVE_results.csv",header=TRUE,sep=',')
#order the traits
PVE_results$Trait = factor(PVE_results$Trait, levels=c("SVL","Head_Width","Head_Length","Pectoral_Width","Pelvic_Width","Fore_Humerus","Fore_Ulna","Fore_Toe_III","Fore_Limb","Hind_Femur","Hind_Tibia","Hind_Toe_IV","Hind_Limb"))

#first barplots for all samples analysis (Fig. 3D)
PVE_results_all_samples<-PVE_results[PVE_results$Dataset == 'all_samples',]
PVE_results_all_samples_hybridization_common<-PVE_results_all_samples[PVE_results_all_samples$Hybrid_sample_group == 'hybridization_common',]
PVE_results_all_samples_hybridization_limited<-PVE_results_all_samples[PVE_results_all_samples$Hybrid_sample_group == 'hybridization_limited',]

#PVE barplots for the hybridization limited group
ggplot(PVE_results_all_samples_hybridization_limited, aes(x=Trait, y=R2))+
  geom_bar(stat="identity", fill="#45569c49", colour="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(-0.8,17)

#PVE barplots for the hybridization common group
ggplot(PVE_results_all_samples_hybridization_common, aes(x=Trait, y=R2))+
  geom_bar(stat="identity", fill="white", colour="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(-0.8,17)


#barplots for samples of Western Cuba ancestry at the 18Mb X chromosome locus (Fig. S26D)
PVE_results_WC_samples<-PVE_results[PVE_results$Dataset == 'WC_X_chrom',]
PVE_results_WC_samples_hybridization_common<-PVE_results_WC_samples[PVE_results_WC_samples$Hybrid_sample_group == 'hybridization_common',]
PVE_results_WC_samples_hybridization_limited<-PVE_results_WC_samples[PVE_results_WC_samples$Hybrid_sample_group == 'hybridization_limited',]

#PVE barplots for the hybridization limited group
ggplot(PVE_results_WC_samples_hybridization_limited, aes(x=Trait, y=R2))+
  geom_bar(stat="identity", fill="#45569c49", colour="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(-0.8,17)

#PVE barplots for the hybridization common group
ggplot(PVE_results_WC_samples_hybridization_common, aes(x=Trait, y=R2))+
  geom_bar(stat="identity", fill="white", colour="black")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  ylim(-0.8,17)