############################################################################################
# Analysis of costs of immune challenge experiment in House Sparrows 2005-6                #
#                                                                                          #
# Bethany Hoye                                                                             #
# 29th November 2021                                                                       #
############################################################################################

library(nlme)
library(lme4)
library(lmerTest)
library(rptR)
library(dplyr) 
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(car)

theme_set(theme_pubclean())

#load data file
HSdata <- read.csv("HoSpimmunecostsdata.csv")
str (HSdata)


#---------------------------------------------------------------------------------------------------------------------#
#       Overview of questions/sections
#---------------------------------------------------------------------------------------------------------------------#

#1 - Differences between treatment groups before starting vaccination series

#2 - Treatment effects on immune parameters

#3 - Treatment effects on testosterone concentration

#4 - Treatment effects on metabolic parameters (and mass)

#
# Section 1 ---------------------------------------------------------------------------------------------------------------------
#       Differences between treatment groups before starting vaccination series
#---------------------------------------------------------------------------------------------------------------------#
# ------
stage2data <- HSdata[HSdata$StageNo ==2,]


#Lysis-----
lm_SRBC_l_sg2_1 <- lm(SRBC_l ~ Inj * Implant * Sex, data = stage2data, REML=FALSE)
summary (lm_SRBC_l_sg2_1)

lm_SRBC_l_sg2_2 <- lm(SRBC_l ~ Inj * Implant + Sex, data = stage2data, REML=FALSE) 
summary(lm_SRBC_l_sg2_2)

lm_SRBC_l_sg2_3 <- lm(SRBC_l ~ Inj + Implant * Sex, data = stage2data, REML=FALSE)
summary(lm_SRBC_l_sg2_3)

lm_SRBC_l_sg2_4 <- lm(SRBC_l ~ Inj * Sex + Implant, data = stage2data, REML=FALSE) 
summary(lm_SRBC_l_sg2_4)

lm_SRBC_l_sg2_5 <- lm(SRBC_l ~ Inj + Implant + Sex, data = stage2data, REML=FALSE) 
summary(lm_SRBC_l_sg2_5)

lm_SRBC_l_sg2_6 <- lm(SRBC_l ~ Implant * Sex, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_6)

lm_SRBC_l_sg2_7 <- lm(SRBC_l ~ Inj * Sex, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_7)

lm_SRBC_l_sg2_8 <- lm(SRBC_l ~ Inj * Implant, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_8)

lm_SRBC_l_sg2_9 <- lm(SRBC_l ~ Inj, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_9)

lm_SRBC_l_sg2_10 <- lm(SRBC_l ~ Implant, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_10)

lm_SRBC_l_sg2_11 <- lm(SRBC_l ~ Sex, data = stage2data, REML=FALSE)  
summary(lm_SRBC_l_sg2_11)

AIC(lm_SRBC_l_sg2_1, lm_SRBC_l_sg2_2, lm_SRBC_l_sg2_3, lm_SRBC_l_sg2_4, lm_SRBC_l_sg2_5, lm_SRBC_l_sg2_6, lm_SRBC_l_sg2_7, lm_SRBC_l_sg2_8, lm_SRBC_l_sg2_9,  lm_SRBC_l_sg2_10,  lm_SRBC_l_sg2_11, k=2)


#Agglutination-----
lm_SRBC_a_sg2_1 <- lm(SRBC_a ~ Inj * Implant * Sex, data = stage2data)
summary (lm_SRBC_a_sg2_1)

lm_SRBC_a_sg2_2 <- lm(SRBC_a ~ Inj * Implant + Sex, data = stage2data) 
summary(lm_SRBC_a_sg2_2)

lm_SRBC_a_sg2_3 <- lm(SRBC_a ~ Inj + Implant * Sex, data = stage2data)
summary(lm_SRBC_a_sg2_3)

lm_SRBC_a_sg2_4 <- lm(SRBC_a ~ Inj * Sex + Implant, data = stage2data) 
summary(lm_SRBC_a_sg2_4)

lm_SRBC_a_sg2_5 <- lm(SRBC_a ~ Inj + Implant + Sex, data = stage2data) 
summary(lm_SRBC_a_sg2_5)

lm_SRBC_a_sg2_6 <- lm(SRBC_a ~ Implant * Sex, data = stage2data)  
summary(lm_SRBC_a_sg2_6)

lm_SRBC_a_sg2_7 <- lm(SRBC_a ~ Inj * Sex, data = stage2data)  
summary(lm_SRBC_a_sg2_7)

lm_SRBC_a_sg2_8 <- lm(SRBC_a ~ Inj * Implant, data = stage2data)  
summary(lm_SRBC_a_sg2_8)

lm_SRBC_a_sg2_9 <- lm(SRBC_a ~ Inj, data = stage2data)  
summary(lm_SRBC_a_sg2_9)

lm_SRBC_a_sg2_10 <- lm(SRBC_a ~ Implant, data = stage2data)  
summary(lm_SRBC_a_sg2_10)

lm_SRBC_a_sg2_11 <- lm(SRBC_a ~ Sex, data = stage2data)  
summary(lm_SRBC_a_sg2_11)

AIC(lm_SRBC_a_sg2_1, lm_SRBC_a_sg2_2, lm_SRBC_a_sg2_3, lm_SRBC_a_sg2_4, lm_SRBC_a_sg2_5, lm_SRBC_a_sg2_6, lm_SRBC_a_sg2_7, lm_SRBC_a_sg2_8, lm_SRBC_a_sg2_9,  lm_SRBC_a_sg2_10,  lm_SRBC_a_sg2_11, k=2)

#KLH-----
lm_KLH_sg2_1 <- lm(KLH ~ Inj * Implant * Sex, data = stage2data)
summary (lm_KLH_sg2_1)

lm_KLH_sg2_2 <- lm(KLH ~ Inj * Implant + Sex, data = stage2data) 
summary(lm_KLH_sg2_2)

lm_KLH_sg2_3 <- lm(KLH ~ Inj + Implant * Sex, data = stage2data)
summary(lm_KLH_sg2_3)

lm_KLH_sg2_4 <- lm(KLH ~ Inj * Sex + Implant, data = stage2data) 
summary(lm_KLH_sg2_4)

lm_KLH_sg2_5 <- lm(KLH ~ Inj + Implant + Sex, data = stage2data) 
summary(lm_KLH_sg2_5)

lm_KLH_sg2_6 <- lm(KLH ~ Implant * Sex, data = stage2data)  
summary(lm_KLH_sg2_6)

lm_KLH_sg2_7 <- lm(KLH ~ Inj * Sex, data = stage2data)  
summary(lm_KLH_sg2_7)

lm_KLH_sg2_8 <- lm(KLH ~ Inj * Implant, data = stage2data)  
summary(lm_KLH_sg2_8)

lm_KLH_sg2_9 <- lm(KLH ~ Inj, data = stage2data)  
summary(lm_KLH_sg2_9)

lm_KLH_sg2_10 <- lm(KLH ~ Implant, data = stage2data)  
summary(lm_KLH_sg2_10)

lm_KLH_sg2_11 <- lm(KLH ~ Sex, data = stage2data)  
summary(lm_KLH_sg2_11)

AIC(lm_KLH_sg2_1, lm_KLH_sg2_2, lm_KLH_sg2_3, lm_KLH_sg2_4, lm_KLH_sg2_5, lm_KLH_sg2_6, lm_KLH_sg2_7, lm_KLH_sg2_8, lm_KLH_sg2_9,  lm_KLH_sg2_10,  lm_KLH_sg2_11, k=2)

#BMR-----
lm_BMR_sg2_1 <- lm(BMR ~ Inj * Implant * Sex, data = stage2data)
summary (lm_BMR_sg2_1)

lm_BMR_sg2_2 <- lm(BMR ~ Inj * Implant + Sex, data = stage2data) 
summary(lm_BMR_sg2_2)

lm_BMR_sg2_3 <- lm(BMR ~ Inj + Implant * Sex, data = stage2data)
summary(lm_BMR_sg2_3)

lm_BMR_sg2_4 <- lm(BMR ~ Inj * Sex + Implant, data = stage2data) 
summary(lm_BMR_sg2_4)

lm_BMR_sg2_5 <- lm(BMR ~ Inj + Implant + Sex, data = stage2data) 
summary(lm_BMR_sg2_5)

lm_BMR_sg2_6 <- lm(BMR ~ Implant * Sex, data = stage2data)  
summary(lm_BMR_sg2_6)

lm_BMR_sg2_7 <- lm(BMR ~ Inj * Sex, data = stage2data)  
summary(lm_BMR_sg2_7)

lm_BMR_sg2_8 <- lm(BMR ~ Inj * Implant, data = stage2data)  
summary(lm_BMR_sg2_8)

lm_BMR_sg2_9 <- lm(BMR ~ Inj, data = stage2data)  
summary(lm_BMR_sg2_9)

lm_BMR_sg2_10 <- lm(BMR ~ Implant, data = stage2data)  
summary(lm_BMR_sg2_10)

lm_BMR_sg2_11 <- lm(BMR ~ Sex, data = stage2data)  
summary(lm_BMR_sg2_11)

AIC(lm_BMR_sg2_1, lm_BMR_sg2_2, lm_BMR_sg2_3, lm_BMR_sg2_4, lm_BMR_sg2_5, lm_BMR_sg2_6, lm_BMR_sg2_7, lm_BMR_sg2_8, lm_BMR_sg2_9,  lm_BMR_sg2_10,  lm_BMR_sg2_11, k=2)

#MMR-----
lm_MMR_sg2_1 <- lm(PMR ~ Inj * Implant * Sex, data = stage2data)
summary (lm_MMR_sg2_1)

lm_MMR_sg2_2 <- lm(PMR ~ Inj * Implant + Sex, data = stage2data) 
summary(lm_MMR_sg2_2)

lm_MMR_sg2_3 <- lm(PMR ~ Inj + Implant * Sex, data = stage2data)
summary(lm_MMR_sg2_3)

lm_MMR_sg2_4 <- lm(PMR ~ Inj * Sex + Implant, data = stage2data) 
summary(lm_MMR_sg2_4)

lm_MMR_sg2_5 <- lm(PMR ~ Inj + Implant + Sex, data = stage2data) 
summary(lm_MMR_sg2_5)

lm_MMR_sg2_6 <- lm(PMR ~ Implant * Sex, data = stage2data)  
summary(lm_MMR_sg2_6)

lm_MMR_sg2_7 <- lm(PMR ~ Inj * Sex, data = stage2data)  
summary(lm_MMR_sg2_7)

lm_MMR_sg2_8 <- lm(PMR ~ Inj * Implant, data = stage2data)  
summary(lm_MMR_sg2_8)

lm_MMR_sg2_9 <- lm(PMR ~ Inj, data = stage2data)  
summary(lm_MMR_sg2_9)

lm_MMR_sg2_10 <- lm(PMR ~ Implant, data = stage2data)  
summary(lm_MMR_sg2_10)

lm_MMR_sg2_11 <- lm(PMR ~ Sex, data = stage2data)  
summary(lm_MMR_sg2_11)

AIC(lm_MMR_sg2_1, lm_MMR_sg2_2, lm_MMR_sg2_3, lm_MMR_sg2_4, lm_MMR_sg2_5, lm_MMR_sg2_6, lm_MMR_sg2_7, lm_MMR_sg2_8, lm_MMR_sg2_9,  lm_MMR_sg2_10,  lm_MMR_sg2_11, k=2)


#Mb-----
lm_Mb_sg2_1 <- lm(Mb ~ Inj * Implant * Sex, data = stage2data)
summary (lm_Mb_sg2_1)

lm_Mb_sg2_2 <- lm(Mb ~ Inj * Implant + Sex, data = stage2data) 
summary(lm_Mb_sg2_2)

lm_Mb_sg2_3 <- lm(Mb ~ Inj + Implant * Sex, data = stage2data)
summary(lm_Mb_sg2_3)

lm_Mb_sg2_4 <- lm(Mb ~ Inj * Sex + Implant, data = stage2data) 
summary(lm_Mb_sg2_4)

lm_Mb_sg2_5 <- lm(Mb ~ Inj + Implant + Sex, data = stage2data) 
summary(lm_Mb_sg2_5)

lm_Mb_sg2_6 <- lm(Mb ~ Implant * Sex, data = stage2data)  
summary(lm_Mb_sg2_6)

lm_Mb_sg2_7 <- lm(Mb ~ Inj * Sex, data = stage2data)  
summary(lm_Mb_sg2_7)

lm_Mb_sg2_8 <- lm(Mb ~ Inj * Implant, data = stage2data)  
summary(lm_Mb_sg2_8)

lm_Mb_sg2_9 <- lm(Mb ~ Inj, data = stage2data)  
summary(lm_Mb_sg2_9)

lm_Mb_sg2_10 <- lm(Mb ~ Implant, data = stage2data)  
summary(lm_Mb_sg2_10)

lm_Mb_sg2_11 <- lm(Mb ~ Sex, data = stage2data)  
summary(lm_Mb_sg2_11)

AIC(lm_Mb_sg2_1, lm_Mb_sg2_2, lm_Mb_sg2_3, lm_Mb_sg2_4, lm_Mb_sg2_5, lm_Mb_sg2_6, lm_Mb_sg2_7, lm_Mb_sg2_8, lm_Mb_sg2_9,  lm_Mb_sg2_10,  lm_Mb_sg2_11, k=2)

#

# Section 2 ---------------------------------------------------------------------------------------------------------------------
#       Treatment effects on immune parameters
#---------------------------------------------------------------------------------------------------------------------#
# ------
postinjdata <- HSdata[HSdata$StageNo >=3,]

#Lysis-------
lmm_Ly1_p <- lmer(chgLys ~ Inj * StageNo + Implant + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly1_p)
par(mfrow = c(2, 2))
plot(lmm_Ly1_p)

lmm_Ly2_p <- lmer(chgLys ~ Inj + StageNo + Implant * Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly2_p)
plot(lmm_Ly2_p)

lmm_Ly3_p <- lmer(chgLys ~ Inj * Implant + StageNo + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly3_p)
plot(lmm_Ly3_p)

lmm_Ly4_p <- lmer(chgLys ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly4_p)
plot(lmm_Ly4_p)

lmm_Ly5_p <- lmer(chgLys ~ Inj + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly5_p)
plot(lmm_Ly5_p)

lmm_Ly6_p <- lmer(chgLys ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ly6_p)
plot(lmm_Ly6_p)


anova(lmm_Ly1_p, lmm_Ly2_p, lmm_Ly3_p, lmm_Ly4_p, lmm_Ly5_p, lmm_Ly6_p)

anova(lmm_Ly4_p, type=2)
anova(lmm_Ly5_p, type=2)


#Agglutination-------
lmm_Ag1_p <- lmer(chgAgg ~ Inj * StageNo + Implant + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag1_p)
plot(lmm_Ag1_p)

lmm_Ag2_p <- lmer(chgAgg ~ Inj + StageNo + Implant * Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag2_p)
plot(lmm_Ag2_p)

lmm_Ag3_p <- lmer(chgAgg ~ Inj * Implant + StageNo + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag3_p)
plot(lmm_Ag3_p)

lmm_Ag4_p <- lmer(chgAgg ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag4_p)
plot(lmm_Ag4_p)

lmm_Ag5_p <- lmer(chgAgg ~ Inj + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag5_p)
plot(lmm_Ag5_p)

lmm_Ag6_p <- lmer(chgAgg ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_Ag6_p)
plot(lmm_Ag6_p)


anova(lmm_Ag1_p, lmm_Ag2_p, lmm_Ag3_p, lmm_Ag4_p, lmm_Ag5_p, lmm_Ag6_p)

anova(lmm_Ag4_p, type=2)


#KLH-------
lmm_KLH1_p <- lmer(chgKLH ~ Inj * StageNo + Implant + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH1_p)
plot(lmm_KLH1_p)

lmm_KLH2_p <- lmer(chgKLH ~ Inj + StageNo + Implant * Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH2_p)
plot(lmm_KLH2_p)

lmm_KLH3_p <- lmer(chgKLH ~ Inj * Implant + StageNo + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH3_p)
plot(lmm_KLH3_p)

lmm_KLH4_p <- lmer(chgKLH ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH4_p)
plot(lmm_KLH4_p)

lmm_KLH5_p <- lmer(chgKLH ~ Inj + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH5_p)
plot(lmm_KLH5_p)

lmm_KLH6_p <- lmer(chgKLH ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)
summary(lmm_KLH6_p)
plot(lmm_KLH6_p)


anova(lmm_KLH1_p, lmm_KLH2_p, lmm_KLH3_p, lmm_KLH4_p, lmm_KLH5_p, lmm_KLH6_p)

anova(lmm_KLH1_p, type=2)
anova(lmm_KLH4_p, type=2)



##-------FIGURE - 3-panel violin + boxplot depicting change in immune scores wrt vaccination ----------

Lys <- ggplot(postinjdata, aes(x = Stage, y = chgLys)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=1) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = c("Sham injected", "Antigen injected"), values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(NULL), y=expression(-log[2]~dilution)) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), 
        legend.position = "top", legend.key = element_rect(fill = "white"), legend.direction = "horizontal") +
  coord_cartesian(ylim = c(-2, 5)) + geom_hline(yintercept=0)

Agg <- ggplot(postinjdata, aes(x = Stage, y = chgAgg)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=1) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = NULL, values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(NULL), y=expression(-log[2]~dilution)) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "none") +
  coord_cartesian(ylim = c(-3, 12)) + geom_hline(yintercept=0)   

KLH <- ggplot(postinjdata, aes(x = Stage, y = chgKLH)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=1) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = NULL, values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(Experimental~stage), y=expression(percent~positive~control~OD[450~nm])) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "none") +
  coord_cartesian(ylim = c(-2, 400)) + geom_hline(yintercept=0)


ggarrange(Lys, Agg, KLH, 
          labels = c("A) ?? Lysis score", "B)?? Agglutination score", "C)?? KLH concentration"), label.x =0, label.y = 1,
          ncol = 1, nrow =3, align = "v", common.legend = TRUE)
# Section 3 ---------------------------------------------------------------------------------------------------------------------
#       Treatment effects on T concentration
#---------------------------------------------------------------------------------------------------------------------#
# ------


lmm_T_conc_1 <- lmer(T_conc ~ Inj * Implant + Sex + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)  
summary(lmm_T_conc_1)

lmm_T_conc_2 <- lmer(T_conc ~ Inj + Implant * Sex + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE) 
summary(lmm_T_conc_2)

lmm_T_conc_3 <- lmer(T_conc ~ Inj + Implant + Sex + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)  
summary(lmm_T_conc_3)

lmm_T_conc_4 <- lmer(T_conc ~ Implant + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)  
summary(lmm_T_conc_4)

lmm_T_conc_5 <- lmer(T_conc ~ Implant + (1 | BirdID), data=postinjdata, REML=FALSE) 
summary(lmm_T_conc_5)

lmm_T_conc_6 <- lmer(T_conc ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)  
summary(lmm_T_conc_6)

lmm_T_conc_7 <- lmer(T_conc ~ Sex + (1 | BirdID), data=postinjdata, REML=FALSE)  
summary(lmm_T_conc_7)


AIC(lmm_T_conc_1, lmm_T_conc_2, lmm_T_conc_3, lmm_T_conc_4, lmm_T_conc_5, lmm_T_conc_6, lmm_T_conc_7, k=2)
anova(lmm_T_conc_5, type=2)

#Mean
aggregate(postinjdata$T_conc, by = list(postinjdata$Imp, postinjdata$Sex), mean, na.rm = T)
#SD
aggregate(postinjdata$T_conc, by = list(postinjdata$Imp, postinjdata$Sex), sd, na.rm = T)

##-------FIGURE - 1-panel violin + boxplot depicting change in T across experimental stages----------

Testosterone <- ggplot(postinjdata, aes(x = Stage, y = T_conc)) + geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 4, ymax = 6), alpha = 0.01, fill = "grey") + geom_violin(aes(fill = Implant), position = position_dodge(0.7), width=1) + 
  geom_boxplot(aes(fill = Implant), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = c("control implant", "Testosterone implant"), values = c("#999999", "black")) +
  labs(fill = NULL, x=expression(NULL), y=expression(Plasma~T~concentration~(ng/ml))) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "none") +
  coord_cartesian(ylim = c(0, 20)) 

ggarrange(Testosterone, label.x = 0.08, label.y = 1,ncol = 1, nrow =1, align = "v", common.legend = TRUE)
# ------

# Section 4 ---------------------------------------------------------------------------------------------------------------------
#       Treatment effects on metabolic parameters (and mass)
#---------------------------------------------------------------------------------------------------------------------#
# ------

#BMR-------
lmm_dBMR1 <- lmer(dBMR ~ Inj * StageNo + Implant + Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR1)
summary(lmm_dBMR1)
anova(lmm_dBMR1, type=2)

lmm_dBMR2 <- lmer(dBMR ~ Inj + StageNo + Implant * Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR2)
summary(lmm_dBMR2)
anova(lmm_dBMR2, type=2)

lmm_dBMR3 <- lmer(dBMR ~ Inj + StageNo + Implant + Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR3)
summary(lmm_dBMR3)
anova(lmm_dBMR3, type=2)

lmm_dBMR4 <- lmer(dBMR ~ Inj * StageNo + Implant + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR4)
summary(lmm_dBMR4)
anova(lmm_dBMR4, type=2)

lmm_dBMR5 <- lmer(dBMR ~ Inj * StageNo + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR5)
anova(lmm_dBMR5, type=2)

lmm_dBMR6 <- lmer(dBMR ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dBMR6)
anova(lmm_dBMR6, type=2)

lmm_dBMR7 <- lmer(dBMR ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE) 
plot(lmm_dBMR7)
anova(lmm_dBMR7, type=2)
summary(lmm_dBMR7)


anova(lmm_dBMR1, lmm_dBMR2, lmm_dBMR3, lmm_dBMR4, lmm_dBMR5, lmm_dBMR6, lmm_dBMR7)

#MMR-------
lmm_dMMR1 <- lmer(dMMR ~ Inj * StageNo + Implant + Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR1)
anova(lmm_dMMR1, type=2)

lmm_dMMR2 <- lmer(dMMR ~ Inj + StageNo + Implant * Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR2)
anova(lmm_dMMR2, type=2)

lmm_dMMR3 <- lmer(dMMR ~ Inj + StageNo + Implant + Sex + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR3)
anova(lmm_dMMR3, type=2)

lmm_dMMR4 <- lmer(dMMR ~ Inj * StageNo + Implant + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR4)
anova(lmm_dMMR4, type=2)
summary(lmm_dMMR4)

lmm_dMMR5 <- lmer(dMMR ~ Inj * StageNo + Mb + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR5)
anova(lmm_dMMR5, type=2)
summary(lmm_dMMR5)

lmm_dMMR6 <- lmer(dMMR ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR6)
anova(lmm_dMMR6, type=2)
summary(lmm_dMMR6)

lmm_dMMR7 <- lmer(dMMR ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_dMMR7)
summary(lmm_dMMR7)
anova(lmm_dMMR7, type=2)


anova(lmm_dMMR1, lmm_dMMR2, lmm_dMMR3, lmm_dMMR4, lmm_dMMR5, lmm_dMMR6, lmm_dMMR7, k=2)

#Mb-------
lmm_Mb1 <- lmer(dMb ~ Inj * StageNo + Implant + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb1)
summary(lmm_Mb1)
anova(lmm_Mb1, type=2)

lmm_Mb2 <- lmer(dMb ~ Inj + StageNo + Implant * Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb2)
summary(lmm_Mb2)
anova(lmm_Mb2, type=2)

lmm_Mb3 <- lmer(dMb ~ Inj + StageNo + Implant + Sex + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb3)
anova(lmm_Mb3, type=2)

lmm_Mb4 <- lmer(dMb ~ Inj * StageNo + Implant + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb4)
anova(lmm_Mb4, type=2)

lmm_Mb5 <- lmer(dMb ~ Inj * StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb5)
anova(lmm_Mb5, type=2)

lmm_Mb6 <- lmer(dMb ~ Inj + StageNo + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb6)
summary(lmm_Mb6)
anova(lmm_Mb6, type=2)

lmm_Mb7 <- lmer(dMb ~ Inj + (1 | BirdID), data=postinjdata, REML=FALSE)
plot(lmm_Mb7)
summary(lmm_Mb7)
anova(lmm_Mb7, type=2)

anova(lmm_Mb1, lmm_Mb2, lmm_Mb3, lmm_Mb4, lmm_Mb5, lmm_Mb6, lmm_Mb7, k=2)




##-------FIGURE - 3-panel violin + boxplot depicting change in metabolic parameters  (absolute) wrt vaccination----------

Mb <- ggplot(postinjdata, aes(x = Stage, y = dMb)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=0.6) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = c("Sham injected", "Antigen injected"), values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(NULL), y=expression(Delta~mass~(g))) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "none") +
  coord_cartesian(ylim = c(-3.5, 3.5)) + geom_hline(yintercept=0)

BMR <- ggplot(postinjdata, aes(x = Stage, y = dBMR)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=0.6) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = NULL, values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(NULL), y=expression(Delta~BMR~(ml~O[2]/min))) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "top", legend.key = element_rect(fill = "white"), legend.direction = "horizontal") +
  coord_cartesian(ylim = c(-1, 1)) + geom_hline(yintercept=0)

PMR <- ggplot(postinjdata, aes(x = Stage, y = dMMR)) + geom_violin(aes(fill = Inj), position = position_dodge(0.7), width=0.6, adjust = 1) + 
  geom_boxplot(aes(fill = Inj), color = "black", position = position_dodge(0.7), width = 0.05) +
  scale_fill_manual(name = NULL, labels = NULL, values = c("#999999", "firebrick")) +
  labs(fill = NULL, x=expression(Experimental~stage), y=expression(Delta~MMR~(ml~O[2]/min))) + 
  theme(panel.grid.major.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(color = "black"), axis.text=element_text(size=12), axis.title = element_text(color="black", face="bold", size=14), legend.position = "none") +
  coord_cartesian(ylim = c(-6.5, 6.5)) + geom_hline(yintercept=0)



ggarrange(Mb, BMR, PMR, labels = c("A)", "B)", "C)"), label.x = 0.08, label.y = 1,ncol = 1, nrow =3, align = "v", common.legend = TRUE)
