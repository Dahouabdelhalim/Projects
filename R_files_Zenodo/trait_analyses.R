#load the required library 
library(ggplot2)
library(ggcorrplot)

#read in morphology data (these data can also be found in Dataset S2)
#this file includes absolute dimensions (with and without log transformation) for each trait as well as relative dimensions (i.e. residuals) which account for the effect of body length (SVL)
morphology<-read.csv(file="Specimen_morphology.csv",header = TRUE, sep = ,)

#SVL corrections were done, for each non-SVL trait, by calculating residuals from linear regressions of log-trait on log-SVL
#as an example, residuals for Hind Limb (HL) are calculated below

#Keep values with data for both SVL and HL
morphologyHL <- na.omit(morphology[c(1,2,4,6,46)])
m_HL<-lm(morphologyHL$Log_Hind_Limb ~ morphologyHL$Log_SVL)
morphologyHL$HL_residuals<-residuals(m_HL)
write.csv(morphologyHL, file="HL_residuals.csv")

#plot trait correlation matrix (as per Fig. S24)
#first select only traits to be included (i.e. SVL and all residuals for other traits), and only samples with data at all traits
#also, this step re-orders the traits for the correlation matrix (i.e. limb length traits are kept together)
morphology_SVL_and_residuals <- na.omit(morphology[c(6,47:49,54,50:53,55:58)])

corr <- round(cor(morphology_SVL_and_residuals), 1)
ggcorrplot(corr, hc.order = FALSE, type = "upper", outline.color = "black", colors= c("black","white","black"),
           lab=TRUE,lab_size = 3,lab_col = "lightgrey")

#get P values and adjust using the Bonferroni method
p.mat <- cor_pmat(morphology_SVL_and_residuals)
p_vals<-as.data.frame(head(p.mat[, 1:13]))

pvals<-c(0.8642461,0.7469612,0.8125080,0.8151469,0.6553477,0.9519572,0.9193875,0.8809040,0.5157847,0.5440344,0.9492267,0.7390825,6.965100e-04,7.149379e-23,2.807087e-12,2.285911e-15,2.711590e-05,3.627204e-08,4.001376e-12,2.816488e-09,1.535400e-07,3.098694e-08,2.113543e-11,0.52598664,0.18499430,0.50207628,0.22585330,0.46723569,0.31114314,0.90528734,0.47505553,0.66295644,0.74697282,3.571353e-15,3.466684e-07,2.301527e-03,1.735526e-06,7.379441e-07,7.750990e-03,1.346218e-03,6.783340e-04,8.706255e-05,2.187464e-06,1.475376e-02,1.514153e-06,4.247111e-06,1.500225e-03,4.027986e-03,1.644185e-03,2.405696e-04,5.561753e-92,8.489167e-42,2.011943e-173,1.564423e-109,9.429265e-96,2.849443e-63,1.620301e-114,3.021056e-43,1.804199e-158,8.174639e-75,1.889985e-100,1.815067e-65,1.024761e-104,1.258512e-116,3.167135e-24,1.076308e-32,9.483020e-96,1.619818e-66,4.393791e-90,1.829449e-103,1.839749e-115,2.943379e-157,1.905490e-128,4.283879e-62,2.227908e-158,4.120042e-86,1.051247e-184,4.408983e-200)
p.adjust(pvals,method = "bonferroni")


#linear regression to investigate the relationship between limb length and perch diameter
#trait values (hind limb residuals) were obtained as above. Additionally, we calculated mean and SE for each population using these residuals, and added population means for perch diameter
#habitat measurements including perch diameter can be found in Dataset S3.

#read in population averages for hind limb length and perch diameter, and split data by transect
morphology_environment<-read.csv(file="HL_residuals_with_perch_diam.csv",header = TRUE, sep = ,)
morphology_environment_transect1 <- morphology_environment[morphology_environment$Transect == 'Transect_1',]
morphology_environment_transect2 <- morphology_environment[morphology_environment$Transect == 'Transect_2',]
morphology_environment_transect3 <- morphology_environment[morphology_environment$Transect == 'Transect_3',]

###plot relationship between limb length and perch diameter for all samples and each transect separately (as per Fig. S27)
ggplot(morphology_environment, aes(Average_log_perch_diameter,Average_HL_residuals))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Average_HL_residuals-SE_HL_residuals, ymax=Average_HL_residuals+SE_HL_residuals), width=0,position=position_dodge(.9))+
  geom_smooth(method='lm', se=FALSE, colour="black", size=1, linetype="dashed")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nlog perch diameter (cm)")+
  ylab("Relative hindlimb length (mm)\\n")+
  ylim(-0.015,0.02)+
  xlim(0.2,1.5)

lm_all<-lm(morphology_environment$Average_HL_residuals~morphology_environment$Average_log_perch_diameter)
summary(lm_all)

ggplot(morphology_environment_transect1, aes(Average_log_perch_diameter,Average_HL_residuals))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Average_HL_residuals-SE_HL_residuals, ymax=Average_HL_residuals+SE_HL_residuals), width=0,position=position_dodge(.9))+
  geom_smooth(method='lm', se=FALSE, colour="black", size=1, linetype="dashed")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nlog perch diameter (cm)")+
  ylab("Relative hindlimb length (mm)\\n")+
  ylim(-0.015,0.02)+
  xlim(0.2,1.5)

lm_transect1<-lm(morphology_environment_transect1$Average_HL_residuals~morphology_environment_transect1$Average_log_perch_diameter)
summary(lm_transect1)

ggplot(morphology_environment_transect2, aes(Average_log_perch_diameter,Average_HL_residuals))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Average_HL_residuals-SE_HL_residuals, ymax=Average_HL_residuals+SE_HL_residuals), width=0,position=position_dodge(.9))+
  geom_smooth(method='lm', se=FALSE, colour="black", size=1, linetype="dashed")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nlog perch diameter (cm)")+
  ylab("Relative hindlimb length (mm)\\n")+
  ylim(-0.015,0.02)+
  xlim(0.2,1.5)

lm_transect2<-lm(morphology_environment_transect2$Average_HL_residuals~morphology_environment_transect2$Average_log_perch_diameter)
summary(lm_transect2)

ggplot(morphology_environment_transect3, aes(Average_log_perch_diameter,Average_HL_residuals))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=Average_HL_residuals-SE_HL_residuals, ymax=Average_HL_residuals+SE_HL_residuals), width=0,position=position_dodge(.9))+
  geom_smooth(method='lm', se=FALSE, colour="black", size=1, linetype="dashed")+
  theme_bw()+
  theme(axis.title=element_text(size=12),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")+
  xlab("\\nlog perch diameter (cm)")+
  ylab("Relative hindlimb length (mm)\\n")+
  ylim(-0.015,0.02)+
  xlim(0.2,1.5)

lm_transect3<-lm(morphology_environment_transect3$Average_HL_residuals~morphology_environment_transect3$Average_log_perch_diameter)
summary(lm_transect3)
