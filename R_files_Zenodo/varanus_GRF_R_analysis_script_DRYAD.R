## This code was used to analyze the data published in "Ground reaction forces in monitor lizards (Varanidae) and the scaling of
## locomotion in sprawling tetrapods," submitted to Biology Letters by Robert L. Cieri et al.

## Set up####

library(tidyverse)
library(patchwork)
library(openxlsx)
library(viridis)
library(lme4)
library(lmerTest)
library(ggeffects)
library(scales)
library(MuMIn)
library(ggiraph)
library(jtools)
library(interactions)
library(plotly)
library(reshape2)
library(ggtern)

#### Import Data ####

allDataFull <- read.xlsx('rawDataSupplemental.xlsx')


##### Data Processing ####

colnames(allDataFull)

# make body-scaled impulses
allDataFull <- allDataFull %>%
  mutate(logGRFintZ = log10(GRFintZ)) %>%
  ## for x positive forces
  mutate(logGRF_PosX = log10(GRFintPosX)) %>%
  ## for x negative forces
  mutate(logGRF_NegX = log10(GRFintNegX)) %>%
  ## for X net forces
  mutate(logGRFx = log10(GRFintX)) %>%
  ## for y positive forces
  mutate(logGRF_PosY = log10(GRFintPosY)) %>%
  ## for y negative forcess
  mutate(logGRF_NegY = log10(GRFintNegY)) %>%
  ## for y net forces
  mutate(logGRFy = log10(GRFintY)) %>%
  # add logmass
  mutate(logMass = log10(kg)) %>%
  # add logSpeed
  mutate(logSpeed = log10(mPerSec))%>%
  # make swing and stance time
  mutate(stanceTime = length/fFPS) %>% 
  mutate(logSwing = log10(swingLong)) %>%
  mutate(logStance = log10(supFrame)) %>%
  # make  peak forces
  mutate(logMaxX = log10((GRFmaxX))) %>%
  mutate(logMaxY = log10((GRFmaxY))) %>%
  mutate(logMaxZ = log10((GRFmaxZ))) %>%
  # scale GRFint
  mutate(logGRFall = log10(GRFintAll)) %>%
  # Pos and Negs for the Max Xs and Ys
  mutate(logMaxXpos = log10((GRFmaxXpos))) %>%
  mutate(logMaxXneg = log10(abs(GRFmaxXneg))) %>%
  mutate(logMaxYpos = log10((GRFmaxYpos))) %>%
  mutate(logMaxYneg = log10(abs(GRFmaxYneg)))%>%
  mutate(logMagPeak = log10(GRFmagPeak))%>%
  mutate(logMagInt = log10(GRFmagInt))%>%
  # make percentages
  mutate(ZpercGRF = GRFabsZint / GRFintAll * 100)%>%
  # location of peaks
  mutate(PeakZ = maxIndexZ / length) %>%
  mutate(PeakX = maxIndexX / length) %>%
  mutate(PeakY = maxIndexY / length) %>%
  # make Z peak scaled to BM
  mutate(ZMaxBM = GRFmaxZ / kg) %>%
  mutate(relSpeed = mPerSec/speciesSpeed)%>%
  # make time-averaged vertical forces
  mutate(averageZforce = GRFintZ / stanceTime) %>%
  # check timing
  mutate(checkTiming = strideTime/stanceTime) %>%
  # duty factor logg
  mutate(logDuty = log10(dutyFactor))%>%
  # make dynamic speed
  mutate(dySpeed = mPerSec^2/(9.8*SVL))%>%
  mutate(logDySpeed = log10(dySpeed)) %>%
  mutate(logRelSpeed = log10(relSpeed)) %>%
  mutate(XpercGRF = GRFabsXint / GRFintAll * 100) %>%
  mutate(YpercGRF = GRFabsYint / GRFintAll * 100) %>%
  mutate(ZpercGRF = GRFabsZint / GRFintAll * 100) 


#### Basic Parameter Values for Results Section ####

## higest impulses
mean(allDataFull$GRFintZ/allDataFull$kg)
min(allDataFull$GRFintZ/allDataFull$kg)
max(allDataFull$GRFintZ/allDataFull$kg)
sd(allDataFull$GRFintZ/allDataFull$kg)/sqrt(length(allDataFull$GRFintZ/allDataFull$kg))

## highest peak forces
mean(allDataFull$GRFmaxZ/allDataFull$kg)
min(allDataFull$GRFmaxZ/allDataFull$kg)
max(allDataFull$GRFmaxZ/allDataFull$kg)
sd(allDataFull$GRFmaxZ/allDataFull$kg)/sqrt(length(allDataFull$GRFmaxZ/allDataFull$kg))

#### Export Species means for PIC ####

colnames(allDataFull)
allDataFull
speciesGroup <- group_by(allDataFull, species)
#view(speciesGroup)

## get basic info on species etc ####

animalMeans_A <-  filter(allDataFull, foot == "hind") %>%
  group_by(species,animal, foot) %>%
  summarise_at(vars(logMass, logSpeed, strideTime, logGRFintZ, logGRF_PosX,
                    logGRF_NegX,  logGRF_PosY, logGRF_NegY, logMaxZ, 
                    logMagInt, logMagPeak, logMaxXpos, logMaxXneg, logMaxYpos,
                    logMaxYneg, ZpercGRF), mean, na.rm = TRUE)
#view(animalMeans_A)
speciesMeans2 <- animalMeans_A %>%
  group_by(species)%>%
  summarise_at(vars(logMass, logSpeed, strideTime, logGRFintZ, logGRF_PosX,
                    logGRF_NegX,  logGRF_PosY, logGRF_NegY, logMaxZ, 
                    logMagInt, logMagPeak, logMaxXpos, logMaxXneg, logMaxYpos,
                    logMaxYneg, ZpercGRF), mean, na.rm = TRUE)

#write.xlsx(speciesMeans2, 'speciesMeansExportHind_A.xlsx')

### get basic numbers for table 1 ####


speciesMeans <- group_by(allDataFull, animal)
summarise_at(speciesMeans, vars(mass), mean)

#brevicauda 
brevi <- unique (allDataFull$mass[which(allDataFull$species=="brev")])
mean(brevi)
sd(brevi)/sqrt(length(brevi))

caudo <- unique(allDataFull$mass[which(allDataFull$species=="caudo")])
mean(caudo)
sd(caudo)/sqrt(length(caudo))

komodo <- unique(allDataFull$mass[which(allDataFull$species=="komodo")])
mean(komodo)
sd(komodo)/sqrt(length(komodo))

panoptes <- unique(allDataFull$mass[which(allDataFull$species=="pano")])
mean(panoptes)
sd(panoptes)/sqrt(length(panoptes))

rosen <- unique(allDataFull$mass[which(allDataFull$species=="rosen")])
mean(rosen)
sd(rosen)/sqrt(length(rosen))

tristis <- unique(allDataFull$mass[which(allDataFull$species=="tristis")])
mean(tristis)
sd(tristis)/sqrt(length(tristis))

varius <- unique(allDataFull$mass[which(allDataFull$species=="varius")])
mean(varius)
sd(varius)/sqrt(length(varius))



#### Model for Duty Factor ####
### using relative speed
allDataFull$animal
dutyFactor_mod_full_1 <- lmer(logDuty ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(dutyFactor_mod_full_1)

dutyFactor_mod_red_1 <- lmer(logDuty ~ logMass + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(dutyFactor_mod_red_1)

anova(dutyFactor_mod_full_1, dutyFactor_mod_red_1)
## NOT significant

summary(dutyFactor_mod_red_1)

sim_slopes(dutyFactor_mod_red_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(dutyFactor_mod_red_1, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.9, digits = 3, pvals = FALSE)

interact_plot(dutyFactor_mod_red_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Duty Factor")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  scale_y_continuous(limits = c(1.3, 2.1))

interact_plot(dutyFactor_mod_red_1, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Relative Speed", y = "Duty Factor")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  scale_y_continuous(limits = c(1.3, 2.1))
#### Model for Total Impulse  ####

colnames(allDataFull)

mag_Int_full_1 <- lmer(logMagInt ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(mag_Int_full_1)

mag_Int_red_1 <- lmer(logMagInt ~ logMass + logRelSpeed + foot +  (1 |animal) , data = allDataFull)
summary(mag_Int_red_1)

anova(mag_Int_red_1,mag_Int_full_1)
## very significant, start dropping terms

summary(mag_Int_full_1)
mag_Int_full_2 <- update(mag_Int_full_1, ~ . -logMass:logRelSpeed:foot)

anova(mag_Int_full_2,mag_Int_full_1)
# not sig
summary(mag_Int_full_2)

mag_Int_full_3 <- update(mag_Int_full_2, ~ . -logRelSpeed:foot)
anova(mag_Int_full_2,mag_Int_full_3)
summary(mag_Int_full_3)

mag_Int_full_4 <- update(mag_Int_full_3, ~ . -logMass:logRelSpeed)
anova(mag_Int_full_4,mag_Int_full_3)
summary(mag_Int_full_4)

sim_slopes(mag_Int_full_4, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(mag_Int_full_4, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(mag_Int_full_3, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept= -0.6), size = 1)

interact_plot(mag_Int_full_3, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept= -0.5), size = 1)

# now for dynamic speed

mag_Int_full_dy_1 <- lmer(logMagInt ~ logMass * logDySpeed * foot + (1 |animal) , data = allDataFull)
summary(mag_Int_full_dy_1)

mag_Int_red_dy_1 <- lmer(logMagInt ~ logMass + logDySpeed + foot +  (1 |animal) , data = allDataFull)
summary(mag_Int_red_dy_1)

anova(mag_Int_red_dy_1,mag_Int_full_dy_1)


## very significant, start dropping terms

summary(mag_Int_full_dy_1)
mag_Int_full_dy_2 <- update(mag_Int_full_dy_1, ~ . -logMass:logDySpeed:foot)

anova(mag_Int_full_dy_2,mag_Int_full_dy_1)

summary(mag_Int_full_dy_2)
mag_Int_full_dy_3 <- update(mag_Int_full_dy_2, ~ . -logDySpeed:foot)
anova(mag_Int_full_dy_2,mag_Int_full_dy_3)

summary(mag_Int_full_dy_3)

sim_slopes(mag_Int_full_dy_3, pred=logMass, modx = foot, cond.int = TRUE)
sim_slopes(mag_Int_full_dy_3, pred=logRelSpeed, modx = foot, cond.int = TRUE)


#### Model for Total Peak ####
logMagPeak

mag_Peak_full_1 <- lmer(logMagPeak ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(mag_Peak_full_1)

mag_Peak_red_1 <- lmer(logMagPeak ~ logMass + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(mag_Peak_red_1)

anova(mag_Peak_full_1, mag_Peak_red_1)
# significant

summary(mag_Peak_full_1)
mag_Peak_full_2 <- update(mag_Peak_full_1, ~ . -foot)
anova(mag_Peak_full_2,mag_Peak_full_1)
summary(mag_Peak_full_2)

sim_slopes(mag_Peak_full_2, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(mag_Peak_full_2, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(mag_Peak_full_2, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)

interact_plot(mag_Peak_full_2, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)

## plot interactions with interact_plot function

#plot with mass
interact_plot(mag_Peak_full_2, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept= -0.5), size = 1)

#plot with speed
interact_plot(mag_Peak_full_2, pred= logSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Speed", y = "Log 10 Mass Peak Force (N))")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept= -0.5), size = 1)

#### Vertical Peak Model ####  

Z_peak_full_1 <- lmer(logMaxZ ~ logMass  * logRelSpeed *  foot + (1 |animal) , data = allDataFull)
summary(Z_peak_full_1)
z_peak_red_1 <-  lmer(logMaxZ ~ logMass  + logRelSpeed +   foot + (1 |animal) , data = allDataFull)
summary(z_peak_red_1)

anova(Z_peak_full_1, z_peak_red_1)
#significant

summary(Z_peak_full_1)
Z_peak_full_2 <- update(Z_peak_full_1, ~ . -logMass:logRelSpeed)
anova(Z_peak_full_1, Z_peak_full_2)
summary(Z_peak_full_2)

summary(Z_peak_full_2)
Z_peak_full_3 <- update(Z_peak_full_2, ~ . -logMass:logRelSpeed:foot)
anova(Z_peak_full_3, Z_peak_full_2)
summary(Z_peak_full_3)

Z_peak_full_4 <- update(Z_peak_full_3, ~ . -logRelSpeed:foot)
anova(Z_peak_full_3, Z_peak_full_4)
summary(Z_peak_full_4)


sim_slopes(Z_peak_full_4, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(Z_peak_full_4, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(Z_peak_full_3, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 3.5, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=0.75)+
  geom_abline(data=allDataFull,aes(slope=1,intercept=0.95), size = 1, color = "darkred", linetype = "dashed") 
#geom_abline(data=allDataFull,aes(slope=0.67, intercept=1.2), size = 1, color="darkred", linetype = "dashed")

interact_plot(Z_peak_full_3, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))



#### Vertical Impulse Model ####  

Z_imp_full_1 <- lmer(logGRFintZ ~ logMass  * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(Z_imp_full_1)

# try model with no interaction terms
Z_imp_red_1 <-  lmer(logGRFintZ ~ logMass  + logRelSpeed + foot +  (1 |animal) , data = allDataFull)
summary(Z_imp_red_1)

anova(Z_imp_full_1, Z_imp_red_1)
## signicficant

summary(Z_imp_full_1)
Z_imp_full_2 <- update(Z_imp_full_1, ~ . -logRelSpeed:foot)
anova(Z_imp_full_1, Z_imp_full_2)
summary(Z_imp_full_2)

Z_imp_full_3 <- update(Z_imp_full_2, ~ . -logMass:logRelSpeed:foot)
anova(Z_imp_full_3, Z_imp_full_2)
summary(Z_imp_full_3)

Z_imp_full_4 <- update(Z_imp_full_3, ~ . -logMass:logRelSpeed)
anova(Z_imp_full_3, Z_imp_full_4)
summary(Z_imp_full_4)

sim_slopes(Z_imp_full_4, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(Z_imp_full_4, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)


interact_plot(Z_imp_full_4, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 3.5, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Impulse (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=0.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept=-0.85), size = 1, color = "darkred", linetype = "dashed") +
  geom_abline(data=allDataFull,aes(slope=1.167,intercept=-0.35), size = 1, color = "darkred", linetype = "dashed") 

interact_plot(Z_imp_full_4, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 3.5, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Impulse (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))

#### Cranial Peak Model####

allDataFull$logMaxXpos[which(allDataFull$logMaxXpos == "-Inf")] <- NA
allDataFull$logMaxXpos

x_peak_red_pos_1 <- lmer(logMaxXpos ~ logMass + logRelSpeed + foot + (1  |animal) , data = allDataFull)
summary(x_peak_red_pos_1)

x_peak_full_pos_1 <- lmer(logMaxXpos ~ logMass * logRelSpeed * foot + (1  |animal) , data = allDataFull)
summary(x_peak_full_pos_1)

anova(x_peak_full_pos_1, x_peak_red_pos_1)
# no diff, use reduced model

summary(x_peak_red_pos_1)

sim_slopes(x_peak_red_pos_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3)

interact_plot(x_peak_red_pos_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  geom_abline(data=allDataFull,aes(slope=1,intercept=0.43), size = 1, color = "darkred", linetype = "dashed") 

sim_slopes(x_peak_red_pos_1, pred=logRelSpeed, modx = foot, cond.int = TRUE)

interact_plot(x_peak_red_pos_1, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))


#### Cranial Impulse Model ####
allDataFull$logGRF_PosX[which(allDataFull$logGRF_PosX == "-Inf")] <- NA
allDataFull$logGRF_PosX

# full model 
cran_imp_full_1 <- lmer(logGRF_PosX ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(cran_imp_full_1)

cran_imp_red_1 <- lmer(logGRF_PosX ~ logMass  + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(cran_imp_red_1)

anova(cran_imp_full_1, cran_imp_red_1)
# NOT SIG
# use reduced model 

summary(cran_imp_red_1)

sim_slopes(cran_imp_red_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(cran_imp_red_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 3.5, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Impulse (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  coord_fixed(ratio=0.5)+
  geom_abline(data=allDataFull,aes(slope=1,intercept=-1.7), size = 1, color = "darkred", linetype = "dashed") +
  geom_abline(data=allDataFull,aes(slope=1.167,intercept=-1.28), size = 1, color = "darkred", linetype = "dashed") 





#### Caudal Peak Model ####

allDataFull$logMaxXneg[which(allDataFull$logMaxXneg == "-Inf")] <- NA
allDataFull$logMaxXneg

x_peak_red_neg_stabl_1 <- lmer(logMaxXneg ~ logMass + logRelSpeed +  foot + (1  |animal) , data = allDataFull)
summary(x_peak_red_neg_stabl_1)

x_peak_full_neg_stabl_1 <- lmer(logMaxXneg ~ logMass * logRelSpeed *  foot + (1  |animal) , data = allDataFull)
summary(x_peak_full_neg_stabl_1)

anova(x_peak_red_neg_stabl_1, x_peak_full_neg_stabl_1)
## not sig, use reduced model

summary(x_peak_red_neg_stabl_1)
sim_slopes(x_peak_red_neg_stabl_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

## make plots for mass

interact_plot(x_peak_red_neg_stabl_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  geom_abline(data=allDataFull,aes(slope=1,intercept=0.30), size = 1, color = "darkred", linetype = "dashed") 


## make plots for speed


interact_plot(x_peak_red_neg_stabl_1, pred= logSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Speed (mPerSec)", y = "Log 10 GRFy Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_bw()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = c(0.2, 0.8))+
  theme(legend.text=element_text(size=12))


#### Caudal Impulse Model####

allDataFull$logGRF_NegX[which(allDataFull$logGRF_NegX == "-Inf")] <- NA
allDataFull$logGRF_NegX

# full model 
caud_imp_full_1 <- lmer(logGRF_NegX ~ logMass * logRelSpeed *  foot + (1 |animal) , data = allDataFull)
summary(caud_imp_full_1)

caud_imp_red_1 <- lmer(logGRF_NegX ~ logMass  + logRelSpeed +  foot + (1 |animal) , data = allDataFull)
summary(cran_imp_red_1)

anova(caud_imp_full_1, cran_imp_red_1)
# not sig so used reduced model

summary(caud_imp_red_1)

#plots

sim_slopes(caud_imp_red_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(cran_imp_red_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  geom_abline(data=allDataFull,aes(slope=1,intercept=-1.65), size = 1, color = "darkred", linetype = "dashed") +
  geom_abline(data=allDataFull,aes(slope=0.67,intercept=-2.6), size = 1, color = "darkred", linetype = "dashed") 


#### Lateral Peak Model####

allDataFull$logMaxYneg[which(allDataFull$logMaxYneg == "-Inf")] <- NA
allDataFull$logMaxYneg

y_peak_full_neg_1 <- lmer(logMaxYneg ~ logMass * logRelSpeed * foot + (1  |animal) , data = allDataFull)
summary(y_peak_full_neg_1)

y_peak_red_neg_1 <- lmer(logMaxYneg ~ logMass + logRelSpeed +  foot + (1  |animal) , data = allDataFull)
summary(y_peak_red_neg_1)

anova(y_peak_full_neg_1, y_peak_red_neg_1)
# use reduced model
summary(y_peak_red_neg_1)

sim_slopes(y_peak_red_neg_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

## make plots for mass

interact_plot(y_peak_red_neg_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  geom_abline(data=allDataFull,aes(slope=1,intercept=0.03), size = 1, color = "darkred", linetype = "dashed") 

## make plots for speed

sim_slopes(y_peak_red_pos_1, pred=logSpeed, modx = foot, cond.int = TRUE)

interact_plot(y_peak_red_pos_1, pred= logSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Speed (mPerSec)", y = "Log 10 GRFy Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = c(0.2, 0.8))+
  theme(legend.text=element_text(size=12))




#### Lateral Impulse Model####

allDataFull$logGRF_NegY[which(allDataFull$logGRF_NegY == "-Inf")] <- NA
allDataFull$logGRF_NegY

lat_imp_full_1 <- lmer(logGRF_NegY ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(lat_imp_full_1)

lat_imp_red_1 <- lmer(logGRF_NegY ~ logMass  + logRelSpeed +  foot + (1 |animal) , data = allDataFull)
summary(lat_imp_red_1)

anova(lat_imp_full_1, lat_imp_red_1)
# use full model and drop terms

summary(lat_imp_full_1)
lat_imp_full_2 <- update(lat_imp_full_1, ~ . - logMass:logRelSpeed:foot)
anova(lat_imp_full_2, lat_imp_full_1)

summary(lat_imp_full_2)
lat_imp_full_3 <- update(lat_imp_full_2, ~ . - logRelSpeed)
anova(lat_imp_full_2, lat_imp_full_3) 

summary(lat_imp_full_3)
lat_imp_full_4 <- update(lat_imp_full_3, ~ . - foot)
anova(lat_imp_full_4, lat_imp_full_3) 

summary(lat_imp_full_4)
lat_imp_full_5 <- update(lat_imp_full_4, ~ . - logMass:logRelSpeed)
anova(lat_imp_full_4, lat_imp_full_5) 

summary(lat_imp_full_5)
lat_imp_full_6 <- update(lat_imp_full_5, ~ . - foot:logRelSpeed)
anova(lat_imp_full_6, lat_imp_full_5) 

summary(lat_imp_full_6)
sim_slopes(lat_imp_full_6, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(lat_imp_full_6, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Body Mass (kg)", y = "Log 10 GRFz Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.position = "blank")+
  theme(legend.text=element_text(size=12))+
  geom_abline(data=allDataFull,aes(slope=1,intercept=-2.2), size = 1, color = "darkred", linetype = "dashed") +
  geom_abline(data=allDataFull,aes(slope=0.67,intercept=-3.1), size = 1, color = "darkred", linetype = "dashed") 



#### Medial Peak Model####

allDataFull$logMaxYpos[which(allDataFull$logMaxYpos == "-Inf")] <- NA
allDataFull$logMaxYpos

y_peak_full_pos_1 <- lmer(logMaxYpos ~ logMass * logRelSpeed *  foot + (1  |animal) , data = allDataFull)
summary(y_peak_full_pos_1)

y_peak_red_pos_1 <- lmer(logMaxYpos ~ logMass + logRelSpeed +  foot + (1  |animal) , data = allDataFull)
summary(y_peak_red_pos_1)

anova(y_peak_full_pos_1, y_peak_red_pos_1)
# not sig so use REDUCED model

summary(y_peak_red_pos_1)

sim_slopes(y_peak_red_pos_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(y_peak_red_pos_1, pred=logRelSpeed, modx = foot,confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

## make plot for mass

interact_plot(y_peak_red_pos_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Mass", y = "Log 10 GRFy Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.text=element_text(size=12))+ 
  theme(legend.position = "none") +
  geom_abline(data=allDataFull,aes(slope=1,intercept=0.54), size = 1, color = "darkred", linetype = "dashed") 

interact_plot(y_peak_red_pos_1, pred= logRelSpeed, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, int.type = "confidence",
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Mass", y = "Log 10 GRFy Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.text=element_text(size=12))+ 
  theme(legend.position = "none") 


#### Medial Impulse Model####

allDataFull$logGRF_PosY
allDataFull$logGRF_PosY[which(allDataFull$logGRF_PosY == "-Inf")] <- NA
allDataFull$logGRF_PosY

# full model 

med_imp_full_1 <- lmer(logGRF_PosY ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(med_imp_full_1)

med_imp_red_1 <- lmer(logGRF_PosY ~ logMass  + logRelSpeed +  foot + (1 |animal) , data = allDataFull)
summary(med_imp_red_1)

anova(med_imp_red_1, med_imp_full_1)
# not significant use the REDUCED model

summary(med_imp_red_1)

sim_slopes(med_imp_red_1, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

interact_plot(med_imp_red_1, pred= logMass, modx = foot, plot.points = TRUE,
              jitter = 0.1, point.shape = TRUE, point.size = 4, interval=TRUE, 
              int.width=0.75, point.alpha = 1)+
  labs(x = "Log 10 Mass", y = "Log 10 GRFy Peak Force (N)")+
  theme(panel.border = element_blank())+
  theme_classic()+
  scale_fill_manual(values = c("black", "black"))+
  scale_color_manual(values =c("#32648EFF", "#74D055FF"))+
  theme(panel.border = element_blank())+
  theme(axis.text.x = element_text(face="bold", size=12, color = "#000000"), 
        axis.text.y = element_text(face="bold", size=12, color = "#000000"))+
  theme(axis.title.x = element_text(size=12, face="bold"), axis.title.y = element_text(size=12, face="bold"))+
  theme(legend.text=element_text(size=12))+ 
  theme(legend.position = "none") +
  geom_abline(data=allDataFull,aes(slope=1,intercept=-1.55), size = 1, color = "darkred", linetype = "dashed") +
  geom_abline(data=allDataFull,aes(slope=1.167, intercept=- 1.118), size = 1, color = "darkred", linetype = "dashed")



#### Proportion GRF Z Model ####

prop_Z_full_1 <- lmer(ZpercGRF ~ logMass  * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(prop_Z_full_1)

prop_Z_red_1 <- lmer(ZpercGRF ~ logMass  + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(prop_Z_red_1)

anova(prop_Z_full_1, prop_Z_red_1)
## not significant so use reduced model

summary(prop_Z_red_1)
sim_slopes(prop_Z_red_1, pred=logMass, modx = foot, cond.int = TRUE)


plot(ZpercGRF ~ logMass, data=allDataFull)


#### model for swing time ####
swing_full_1 <- lmer(logSwing ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(swing_full_1)

swing_red_1 <- lmer(logSwing ~ logMass + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(swing_red_1)

anova(swing_full_1, swing_red_1)

# very sig so drop terms

summary(swing_full_1)
swing_full_2 <- update(swing_full_1, ~ . -logMass:foot)
summary(swing_full_2)

anova(swing_full_1, swing_full_2)
# not sig so keep going

summary(swing_full_2)
swing_full_3 <- update(swing_full_2, ~ . -logRelSpeed:foot)

anova(swing_full_3, swing_full_2)
summary(swing_full_3)

swing_full_4 <- update(swing_full_3, ~ . -logMass:logRelSPeed:foot)
summary(swing_full_4)

sim_slopes(swing_full_4, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(swing_full_4, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)

#### Model for Stance Time ####

stance_full_1 <- lmer(logStance ~ logMass * logRelSpeed * foot + (1 |animal) , data = allDataFull)
summary(stance_full_1)

stance_red_1 <- lmer(logStance ~ logMass + logRelSpeed + foot + (1 |animal) , data = allDataFull)
summary(stance_red_1)

anova(stance_full_1, stance_red_1)
# sig so start dropping shit

summary(stance_full_1)

stance_full_2 <- update(stance_full_1, ~ . -logRelSpeed:foot)

summary(stance_full_2)

anova(stance_full_2, stance_full_1)
# not sig

summary(stance_full_2)
stance_full_3 <- update(stance_full_2, ~ . -logMass:logRelSpeed:foot)

anova(stance_full_2, stance_full_3)

summary(stance_full_3)
stance_full_4 <- update(stance_full_3, ~ . -logMass:foot)

anova(stance_full_4, stance_full_3)

summary(stance_full_4)

sim_slopes(stance_full_4, pred=logMass, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)
sim_slopes(stance_full_4, pred=logRelSpeed, modx = foot, confint = TRUE, ci.width = 0.95, digits = 3, pvals = FALSE)


#### Ternary Plot ####
##  plotting the %ages vs. speed

# for body mass
my_breaks <- c(0, 0.1, 0.5, 2.5, 20)

## for all data
plot1 <- ggtern(data=allDataFull,aes(XpercGRF,YpercGRF,ZpercGRF, shape = foot)) + 
  geom_point(aes(fill = kg, size =kg))  + 
  #geom_point(aes(colour = relSpeed, size = 4))  + 
  labs(x="X",y="Y",z="Z", color = ("Body Mass (kg)"), shape = "Foot")+
  guides(size = FALSE) +
  scale_fill_viridis_c(direction =-1, trans = "log", breaks = my_breaks)+
  theme_clockwise()+
  theme_gray() +
  theme(legend.position = c(0.9, 0.72))+
  theme(axis.title = element_text(face = "bold"))+
  theme(axis.text = element_text(face="bold", size =12))+
  scale_shape_manual(values = c(21, 24))+
  scale_size_continuous(range = c(3.5,8))

#Alternatives, and Arrow Label
plot1 + Tarrowlab("Percent GRFx ") + Larrowlab("Percent GRFy ") + Rarrowlab("Percent GRFz ") +
  theme_showarrows() + Wlab("Total GRF")

#### 3D Plot ####

### Z impulse

my_df <- allDataFull

x <- my_df$logMass
y <- my_df$logRelSpeed
z <- my_df$logGRFintZ

lizard_model_lm <- lm(z ~ 0 + x + y,data = my_df)

graph_reso <- 0.05
axis_x <- seq(min(x), max(x), by = graph_reso)
axis_y <- seq(min(y), max(y), by = graph_reso)

lizard_lm_surface <- expand.grid(x = axis_x, y = axis_y, KEEP.OUT.ATTRS = F)
lizard_lm_surface$logGRFintZ <- predict.lm(lizard_model_lm, newdata = lizard_lm_surface)

lizard_lm_surface <- acast(lizard_lm_surface, y ~ x, value.var = "logGRFintZ") #y ~ x

plot1 <- box

plot1 <-plot_ly(my_df, x= ~x, y= ~y, z= ~z, type="scatter3d", mode="markers", color=~logGRFintZ)
plot1 <- add_trace(p=plot1, z=lizard_lm_surface, x=axis_x, y=axis_y, type = "surface", opacity=0.7)
plot1 <- add_trace(p=plot1, marker = list(line=list(color = "#000000", width = 2)))
plot1 <- plot1 %>% layout(scene = list(xaxis = list(title = "Body Mass"), 
                                       yaxis = list(title = "Relative Speed", autorange = "reversed"),  
                                       zaxis = list(title = "Vertical Impulse")))
plot1 <- plot1 %>% hide_colorbar()
plot1 <- plot1 %>% hide_legend()
plot1

### Total Impulse

f1 <- list(
  family = "Times, sans-serif",
  size = 14,
  color = "black")

f2 <- list(
  family = "Old Standard TT, serif",
  size = 14,
  color = "#black")

axis <- list(
  titlefont = f1,
  tickfont = f2,
  showgrid = F
)


scene = list(
  xaxis = axis,
  yaxis = axis,
  zaxis = axis,
)



my_df <- allDataFull

x <- my_df$logMass
y <- my_df$logRelSpeed
z <- my_df$logMagInt

lizard_model_lm <- lm(z ~ 0 + x + y,data = my_df)

graph_reso <- 0.05
axis_x <- seq(min(x), max(x), by = graph_reso)
axis_y <- seq(min(y), max(y), by = graph_reso)

lizard_lm_surface <- expand.grid(x = axis_x, y = axis_y, KEEP.OUT.ATTRS = F)
lizard_lm_surface$logMagInt <- predict.lm(lizard_model_lm, newdata = lizard_lm_surface)

lizard_lm_surface <- acast(lizard_lm_surface, y ~ x, value.var = "logMagInt") #y ~ x

plot2 <- box
plot2 <-plot_ly(my_df, x= ~x, y= ~y, z= ~z, type="scatter3d", mode="markers", color=~logMagInt)
plot2 <- add_trace(p=plot1, z=lizard_lm_surface, x=axis_x, y=axis_y, type = "surface", opacity=0.7)
plot2 <- add_trace(p=plot1, marker = list(line=list(color = "#000000", width = 2)))
plot2 <- plot1 %>% layout(scene = list(xaxis = list(title = "Body Mass", dtick=1), 
                                       yaxis = list(title = "Relative Speed", autorange = "reversed", dtick=0.5, size =14),  
                                       zaxis = list(title = "Ground Reaction Impulse")))
plot2 <- plot2 %>% hide_colorbar()
plot2 <- plot2 %>% hide_legend()
plot2
