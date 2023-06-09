######## Fox et al. 2021 - Differential resistance and acclimation of two coral species to chronic nutrient
########                   enrichment reflect life history trait. Functional Ecology

####### Code file 3 of 4 -- All main text analyses except for stable isotopes (see code file 4)
####### Required data files: 1) PAM_master_Dryad.csv
#                            2) Phys_master_Dryad.csv
#                            3) There are 5 files used to analysis of the respirometry data that will be loaded
#                               later in the script. These are:
#                               a) CHAIN_Slopes.csv
#                               b) CHAIN_respirometry_chamberID_mid_final.csv
#                               c) TPAIN_master_all_variables_09.17.csv
#                               d) Buoyant_weights_11.24.2015.csv
#                               e) Nubbin_growth_and_SA_05.02.2016.csv


####### Summary: This file completes all of the data summary, analysis, and plotting for the main text
#######          figures (excluding stable isotopes).

####### Figures produced: Fig. 1,2,3 and supplemental Figs. S6,7,8,9

library(dplyr)
library(ggplot2)
require(lmerTest)
library(car) 
require(RColorBrewer)
require(emmeans)
library(viridis)
library(pals)
library(broom)
library(tidyr)
library(lme4)
library(tidyr)
library(purrr)
library(broom)
library(jtools)
library(merTools)
library(cowplot)
library(patchwork)


#create plotting themes
newtheme <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.line = element_line(size = 0))

newtheme_right <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,6), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.line = element_line(size = 0))


####################
##### Fig. 1 - Temporal analysis of PAM data to assess acclimation to chronic nutrient enrichment over 5 weeks
####################

pam<-as.data.frame(read.csv(file="PAM_master_file.csv",header=TRUE))

#remove pre-experiment PAM measurements 
pam2<-subset(pam,Assignment!='T0' & Time>0)

#create continuous nutrient treatments
pam2$Nutrient.Level<-as.factor(pam2$Nutrient.Level)
pam2$Nutrients<-NA
pam2$Nutrients<-ifelse(pam2$Treatment==0,0.1,
                     ifelse(pam2$Treatment=="1",1,
                            ifelse(pam2$Treatment=="2",3,
                                   ifelse(pam2$Treatment=="3",5,
                                          ifelse(pam2$Treatment=="4",7,NA)))))
unique(pam2$Nutrients)

#### Acclimatory effects through time as a fxn of nutrient treatment
# For this we keep nutrients as a continuous factor because we are interested in acclimation to 
# a nutrient gradient over 5 consecutive weeks. Continuous values are our targeted dosing levels. 
# There are 3 plots that follow to explore this process:
# 1 and 2) linear mixed effects model to examine nutrient response (slope) over time (fixed w/ 5 weekly measurements) by Species
# This  is illustrated in 2 plots: A) the linear fit of model predicted results 
#                                     plotted over treatment means (SE) at each time point. 
#                                  B) the temporal change in slope parameter estimate (95% CI) over time 
#                                  C) A comparison of the change in nutrient response slope between species determined via bootstrap resampling (see below)

#3 linear mixed effects models are run to contrast species with - nutrients as continous and time as fixed effects
# random orthogonal effects of tank, colony, and individual to account for 

### Full model to determine if significant interaction between Species*Nutrients*Time, which indicats that the two species
### respond differently to nutrients over 5 weeks. 

### a simplified model is then run for each species in order to identify the time points at which the Nutrients*Time response
### becomes significant (fig 1B)

#assumptions 
qqnorm(pam2$Yield)
qqline(pam2$Yield)
leveneTest(pam2$Yield~as.factor(pam2$Species))

pam.mod<-lmer(Yield~Species*Nutrients*as.factor(Time)+(1|Tank)+(1|Colony)+(1|Nub_ID),data=pam2,REML=F)
summ(pam.mod)
anova(pam.mod,type=2)

#check residuals
plot(pam.mod)
hist(resid(pam.mod))

#significant 3-way interaction btwn Species*Nutrients*Time
# Type II Analysis of Variance Table with Satterthwaite's method
#                                     Sum Sq  Mean Sq NumDF   DenDF  F value    Pr(>F)    
# Species                           0.069905 0.069905     1  243.04 149.6031 < 2.2e-16 ***
# Nutrients                         0.022001 0.022001     1  242.94  47.0836 5.618e-11 ***
# as.factor(Time)                   0.016694 0.004173     4 1008.00   8.9314 4.344e-07 ***
# Species:Nutrients                 0.000730 0.000730     1  243.14   1.5631   0.21241    
# Species:as.factor(Time)           0.032932 0.008233     4 1008.00  17.6195 5.469e-14 ***
# Nutrients:as.factor(Time)         0.032373 0.008093     4 1008.00  17.3200 9.427e-14 ***
###################################
# Species:Nutrients:as.factor(Time) 0.005720 0.001430     4 1008.00   3.0603   0.01606 *    
###################################


#### NOW run species-specific models to establish trends and timing of nutrient responses

poc.pam<-subset(pam2,Species=="Pocillopora")
por.pam<-subset(pam2,Species=="Porites")
#POCILLOPORA
#assumptions 
qqnorm(poc.pam$Yield)
qqline(poc.pam$Yield)
leveneTest(poc.pam$Yield~as.factor(poc.pam$Treatment))

poc.pam.mod<-lmer(Yield~Nutrients*as.factor(Time)+(1|Tank)+(1|Colony)+(1|Nub_ID),data=poc.pam,REML=F)
summary(poc.pam.mod)
anova(poc.pam.mod,type=2)

#check residuals
plot(poc.pam.mod)
hist(resid(poc.pam.mod))

# Type II Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq   Mean Sq NumDF  DenDF F value          Pr(>F)    
# Nutrients                 0.008237 0.0082373     1 116.89 13.3307       0.0003919 ***
# as.factor(Time)           0.002714 0.0006785     4 504.00  1.0981       0.3567713    
# Nutrients:as.factor(Time) 0.031690 0.0079226     4 504.00 12.8215 0.0000000006018 ***

#significant interactions for indicating on slope plot
# Nutrients:as.factor(Time)2  1.505e-03  1.210e-03  5.040e+02   1.244   0.2142    
# Nutrients:as.factor(Time)3  5.129e-03  1.210e-03  5.040e+02   4.238 2.68e-05 ***
# Nutrients:as.factor(Time)4  5.885e-03  1.210e-03  5.040e+02   4.863 1.55e-06 ***
# Nutrients:as.factor(Time)5  7.242e-03  1.210e-03  5.040e+02   5.985 4.12e-09 ***

###  PORITES
#assumptions
qqnorm(por.pam$Yield)
qqline(por.pam$Yield)
leveneTest(por.pam$Yield~as.factor(por.pam$Treatment))

por.pam.mod<-lmer(Yield~Nutrients*as.factor(Time)+(1|Tank)+(1|Colony)+(1|Nub_ID),data=por.pam,REML=F)
summary(por.pam.mod)
anova(por.pam.mod,type=2)

#check residuals
plot(por.pam.mod)
hist(resid(por.pam.mod))

# Type II Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq   Mean Sq NumDF DenDF F value    Pr(>F)    
# Nutrients                 0.015311 0.0153115     1   119  48.358 2.053e-10 ***
# as.factor(Time)           0.046912 0.0117279     4   504  37.040 < 2.2e-16 ***
# Nutrients:as.factor(Time) 0.006402 0.0016005     4   504   5.055 0.0005312 ***

#significant interactions
#   Nutrients:as.factor(Time)2 1.929e-03  8.662e-04 5.040e+02   2.227 0.026373 *  
#   Nutrients:as.factor(Time)3 2.808e-03  8.662e-04 5.040e+02   3.242 0.001267 ** 
#   Nutrients:as.factor(Time)4 3.193e-03  8.662e-04 5.040e+02   3.686 0.000253 ***
#   Nutrients:as.factor(Time)5 3.355e-03  8.662e-04 5.040e+02   3.873 0.000122 ***
# 
#use the species models to predict yield values and CI for plotting -- using predictInterval and following
#https://stats.stackexchange.com/questions/235018/r-extract-and-plot-confidence-intervals-from-a-lmer-object-using-ggplot

poc.pam2<-cbind(poc.pam, predictInterval(poc.pam.mod,poc.pam))
por.pam2<-cbind(por.pam,predictInterval(por.pam.mod,por.pam))

fits<-rbind(poc.pam2,por.pam2)

fits.sum<-fits %>% group_by(Nutrients,Species,Time,Treatment) %>%
  summarize(
    mean = mean(fit,na.rm=TRUE),#mean predicted value per treatment/time
    lower = mean(lwr,na.rm=TRUE),#mean CIs
    upper = mean(upr,na.rm=TRUE))

###now plot estimated fits and linear fits to represent them
fits.sum$Species<-factor(fits.sum$Species,levels=c("Porites","Pocillopora")) #order the species first

fits.plot<-ggplot(fits.sum,aes(facets=Species)) + 
  geom_smooth(aes(Nutrients, mean, group=Time,color=Time),method="lm",se=FALSE,lty=2,size=.5) +
  facet_grid(Species~.)
fits.plot


#now add the means and se of raw data at each time point and nutrient concentration

#First summarize the yield data for weeks 1-5
pam_means<-pam2 %>% group_by(Nutrients,Species,Time,Treatment) %>% filter(Time>0) %>% 
  summarize(
    N    = sum(!is.na(Yield)),
    yield = mean(Yield, na.rm=TRUE),
    se_y   = sd(Yield, na.rm=TRUE) / sqrt(N))

### add summarize points and SE to above plot

pam.plot<-
  fits.plot+
  geom_errorbar(data=pam_means,mapping=aes(Nutrients,ymin=yield-se_y,ymax=yield+se_y,group=Time),width=.1)+
  geom_point(data=pam_means,mapping=aes(Nutrients,yield,fill=Time,group=Time),colour="black",pch=21,size=3)+
  scale_x_continuous(breaks=seq(0,7,1),name= bquote('Nitrate Level ('*µmol~L^-1*')'))+
  scale_y_continuous(breaks=seq(.5,.75,0.05),name="Maximum Quantum Yield (Fv/Fm)")+newtheme+
  theme(legend.position = "none")
pam.plot


### ### ### ### ### ### 
# NOW extract slope estimates for nutrient coefficent at each level of Time using emtrends
### ### ### ### ### ### 

poc.slopes<-as.data.frame(emtrends(poc.pam.mod,"Time",var="Nutrients"))
poc.slopes$Species<-"Pocillopora"

por.slopes<-as.data.frame(emtrends(por.pam.mod,"Time",var="Nutrients"))
por.slopes$Species<-"Porites"

slopes<-rbind(poc.slopes,por.slopes)
slopes$Species<-factor(slopes$Species,levels=c("Porites","Pocillopora")) #order the species first

#plot change in nutrient slope w/ CI through time for each species 
slope.plot<-ggplot(slopes,aes(x=Time,y=Nutrients.trend,ymin=lower.CL,ymax=upper.CL,fill=Time))+
  geom_errorbar(width=0.1)+geom_line(lty=2)+
  geom_point(pch=21,size=3,color="black")+
  facet_grid(Species~.)+
  ylab("Nutrient Response Slope")+
  xlab("Time (Weeks)")+newtheme+
  theme(legend.position = "bottom")
slope.plot


########################################################################################################################################################
######################################
###################################### Boot strap re-sampling for figure 1C to compare model fits to species data over time 
######################################
########################################################################################################################################################


#get sample size for each treatment level (same for both species)
num<-poc.pam %>% group_by(Treatment,Time) %>% summarize(N=sum(!is.na(Yield)))

a<- poc.pam %>% dplyr::select(Time,Yield,Nutrients) 
b<-por.pam %>% dplyr::select(Time,Yield,Nutrients) 

nested<-a %>% group_by(Nutrients,Time) %>% nest() %>% ungroup() %>% 
  mutate(n=c(rep(c(28,21,28,21,28),5)))

nested.b<-b %>% group_by(Nutrients,Time) %>% nest() %>% ungroup() %>% 
  mutate(n=c(rep(c(28,21,28,21,28),5)))


all<-NULL

n<-10#to change the number of permutations 10000 was used in the manuscript
# 10 is shown here as an example to make the code run faster

for(i in 1:n){  
  
  # randomly sample each of the source populations with replacement for the number of samples for each group
  samps<-nested %>% mutate(samp=map2(data,n,sample_n,replace=T))
  samps.b<-nested.b %>% mutate(samp=map2(data,n,sample_n,replace=T))
  
  # unnest iteratively sampled lists -- this generates 1 new df with randomly 
  # sampled means acoss pam dataset for treatment level and time period
  samps2<-samps %>% dplyr::select(-data) %>% unnest(samp) %>% dplyr::select(n,Nutrients,Time,Yield) 
  samps2.b<-samps.b %>% dplyr::select(-data) %>% unnest(samp) %>% dplyr::select(n,Nutrients,Time,Yield)
  
  slopes<- samps2 %>% group_by(Time) %>% 
    nest() %>% 
    mutate(
      fit = purrr::map(data, ~ lm(Yield ~ Nutrients, data = .x)),
      tidied = purrr::map(fit, tidy)
    ) %>% 
    unnest(tidied)
  
  
  slopes.b<- samps2.b %>% group_by(Time) %>% 
    nest() %>% 
    mutate(
      fit = purrr::map(data, ~ lm(Yield ~ Nutrients, data = .x)),
      tidied = purrr::map(fit, tidy)
    ) %>% 
    unnest(tidied)
  
  
  #extract slopes for each time period
  slopes1<-slopes %>% filter(term=="Nutrients")
  slopes1$estimate<-round(slopes1$estimate, digits = 4)
  slopes1<- dplyr::select(slopes1,-data,-fit,-term)
  slopes1<-as.data.frame(slopes1)
  slopes1$Iteration<-i
  slopes1$Species<-"POC"
  
  slopes2<-slopes.b %>% filter(term=="Nutrients")
  slopes2$estimate<-round(slopes2$estimate, digits = 4)
  slopes2<- dplyr::select(slopes2,-data,-fit,-term)
  slopes2<-as.data.frame(slopes2)
  slopes2$Iteration<-i
  slopes2$Species<-"POR"
  
  #compile all slope estimates 
  all<-rbind(all,slopes1,slopes2)
}

#summarize the slopes from all models and calculate confidence intervals 
slopes.ci<-all %>%
  group_by(Time,Species) %>% #calculate summary stats for each spp at the class level
  summarize(
    slope = mean(estimate,na.rm=TRUE),
    upper = quantile(estimate,0.975),
    lower = quantile(estimate,0.025),
    upper.75 = quantile(estimate,0.825),
    lower.75 = quantile(estimate,0.125))

#fit linear and polynomial models to mean slope estimates for each species
#PORITES - polynomial is best fit
mod1 <- lm(slope ~ Time,data=subset(slopes.ci,Species=="POR")) #linear regression
tidy(mod1)
mod2 <- lm(slope ~ poly(Time,2),data=subset(slopes.ci,Species=="POR")) #poly
summary(mod2)
#compare models using partial f-test to see if polynomical is significantly better fit
anova(mod1,mod2) ## p = 0.0036 poly is better fit

#POCILLOPORA -- linear = most parsimonious
mod3<-lm(slope ~ Time,data=subset(slopes.ci,Species=="POC"))
mod4<-lm(slope ~ poly(Time,2),data=subset(slopes.ci,Species=="POC"))
anova(mod3,mod4) # p =0.3745, no difference between models. 

#plot bootstrapped mean slopes w/ 95%CI and fitted models to the means
boot.plot<-ggplot(slopes.ci,aes(x=Time,y=slope,ymin=lower,ymax=upper,color=Species,fill=Species))+
  geom_smooth(data=subset(slopes.ci,Species=="POR"),mapping=aes(x=Time,y=slope,group=1,color=Species,fill=Species),method="lm",formula = y ~ poly(x, 2),se=T,alpha=0.25)+
  geom_smooth(data=subset(slopes.ci,Species=="POC"),mapping=aes(x=Time,y=slope,group=1,color=Species,fill=Species),method="lm",se=T,alpha=0.25)+
  geom_errorbar(width=0,position=position_dodge(0.3))+
  geom_errorbar(data=slopes.ci,mapping=aes(ymin=lower.75,ymax=upper.75),size=1.5,width=0,position=position_dodge(0.3))+
  geom_point(data=slopes.ci,mapping=aes(shape=Species),position=position_dodge(0.3),col='black',size=2.5)+scale_shape_manual(values=c(21,24))+
  scale_fill_manual(values=c("#4592af","#ffa41b"))+scale_color_manual(values=c("#4592af","#ffa41b"))+
  ylab("Nutrient Response Slope")+xlab("Time (Weeks)")+scale_y_continuous(breaks=seq(-0.005,0.0075,0.0025))+
  newtheme
boot.plot

#combine all figures into single plot 
Fig.1.layout <- "
AABB
CC##
"
quartz()
Fig.1_full<-pam.plot+slope.plot+boot.plot+plot_layout(design=Fig.1.layout)+plot_annotation(tag_levels = 'A')
Fig.1_full

# ggsave(plot = Fig.1_full,filename = "Fig_1_boot.pdf",
#        width=9,height=9,dpi=300,useDingbats=FALSE)


################################-
### ### LOAD PHYSIOLOGICAL DATA ###
## this portion of code will: 
## 1) correlate fv/fm values with symbiont density and chl concentration
## 2) Conduct all physiolgoical analyses and generate main text figure 2 and figs S6,7,8,
## 3) create the df required for multivarite analyeses conducted in the third section of this script
######
ch<-read.csv(file="Phys_master_Dryad.csv",header=T)

#get rid of fragments used in respirometry trials and pre-experiment measurements
ch<-subset(ch,Nutrient.Level!="T0" & Assignment!="TP")

#take only the last time point and exclude nubbins from respirometry trials to match physiological data
fv<-subset(pam2,Time==5 &Assignment!="TP")

#select just fv/fm and nub ID
pam_clean<-fv %>% dplyr::select(Nub_ID,Yield)

#nub ids are in different format but arranged in order
pam_clean<-arrange(pam_clean,Nub_ID)
ch<-arrange(ch,Nub_ID)

#combine the pam data to the phys data
new<-cbind(ch,pam_clean)
#remove second nub_ID column
new <- new[ -c(27:28) ]

# remove any coral fragments that suffered partial mortality and died during the experiment
dat.clean<-subset(new, Partial.mortality <2) ### use this df for multivariate analyses below 

poc.new<-subset(dat.clean,Species=="POC")
por.new<-subset(dat.clean,Species="POR")

#correlate pam data with symbiont metrics using continuous nutrient treatment for color
#create continuous nutrient treatments
pam.cor<-dat.clean
pam.cor$Nutrient.Level<-as.factor(pam.cor$Nutrient.Level)
pam.cor$Nutrients<-NA
pam.cor$Nutrients<-ifelse(pam.cor$Nutrient.Level=="N0",0.1,
                       ifelse(pam.cor$Nutrient.Level=="N1",1,
                              ifelse(pam.cor$Nutrient.Level=="N2",3,
                                     ifelse(pam.cor$Nutrient.Level=="N3",5,
                                            ifelse(pam.cor$Nutrient.Level=="N4",7,NA)))))


##correlate Fv/Fm with zoox cm2
cp1<-ggplot(pam.cor,aes(x=Zoox.cm2,
               y=Yield,color=Nutrients,group=Species,
               facets=Species))+geom_point()+facet_wrap(~Species,scale='free')+geom_smooth(method="lm")+
  scale_x_continuous(labels = scales::scientific)+newtheme
  

cor.test(poc.new$Zoox.cm2,poc.new$Yield) #p=0.07 r=0.18
cor.test(por.new$Zoox.cm2,por.new$Yield) #0.05, r=-0.14

#chl cm2
cp2<-ggplot(pam.cor,aes(x=Total.chl,
               y=Yield,color=Nutrients,group=Species,
               facets=Species))+geom_point()+facet_wrap(~Species,scale='free')+geom_smooth(method="lm")+newtheme

cor.test(poc.new$Total.chl,poc.new$Yield) #p=<0.001 r=0.37
cor.test(por.new$Total.chl,por.new$Yield) #p=0.89, r=0.01

#chl/zoox
cp3<-ggplot(pam.cor,aes(x=Total.chl/Zoox.cm2,
               y=Yield,color=Nutrients,group=Species,
               facets=Species))+geom_point()+facet_wrap(~Species,scale='free')+
  geom_smooth(method="lm",linetype=2,col="black",se=FALSE)+newtheme

cor.test(poc.new$Total.chl/poc.new$Zoox.cm2,poc.new$Yield) #p=0.0016 r=0.0.32
cor.test(por.new$Total.chl/por.new$Zoox.cm2,por.new$Yield) #p=0.0001, r=0.28

#compile correlations into single supplemental figure
plot_grid(cp1,cp2,cp3,ncol=1,align="h",axis='bt',labels=c("A","B","C"))



#####
######################## PHYSIOLOGICAL ANALYSES - Fig. 2
########## Calcification
#####
ch<-dat.clean
#### normality for growth -- split by species 
por<-subset(ch,Species=="POR")
poc<-subset(ch,Species=="POC")

###### Mixed effects models -- lmer with random factors of colony and Tank

################-
################-
###### Calcification normalized to surface area - ensure growth in control treatments over experiment
###### and determine if POR growth >> POC growth
################-
################-

## test of control different from 0 
porN0=subset(por,Nutrient.Level=="N0")
t.test(porN0$growth.SA)

pocN0=subset(poc,Nutrient.Level=="N0")
t.test(pocN0$growth.SA)

#comparing species baseline calcification
t.test(pocN0$growth.SA,porN0$growth.SA)

#assumptions for POR 
qqnorm(sqrt(por$growth.SA))
qqline(sqrt(por$growth.SA))
leveneTest(sqrt(por$growth.SA)~por$Nutrient.Level)

####### LINEAR MIXED EFFECTS MODELS -- only full models used. Tank replicate and colony as random orthogonal factors
####### No model selection to compare different models between response variables. So REML = TRUE

# Tank and Colony included as orthogonal random effects and kept in all models 

#SURFACE AREA standardized growth (mg CaCO3 cm-2 d-1)
#full
por.growth1<-lmer(sqrt(growth.SA)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))

summary(por.growth1)
anova(por.growth1,type=2)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# Nutrient.Level  1.011 0.25276     4 86.23  6.6535 0.0001026 ***

#calculate error sums of squares
sum(resid(por.growth1)^2)

#extract random effect significance 
ranova(por.growth1)

#check residuals
plot(por.growth1)
hist(resid(por.growth1))

#specific contrasts for each level - Tukey
emm.growthSA = emmeans(por.growth1, ~ Nutrient.Level)
pairs(emm.growthSA)

# N0 - N2   0.21255957 0.06197250 86.24   3.430  0.0081
# N0 - N3   0.24749225 0.06285268 86.37   3.938  0.0015
# N0 - N4   0.20229889 0.06335962 86.11   3.193  0.0164
# N1 - N2   0.17936494 0.06096401 86.06   2.942  0.0332
# N1 - N3   0.21429762 0.06174031 86.09   3.471  0.0071
# N1 - N4   0.16910427 0.06270296 86.23   2.697  0.0626 # indicated in fig with *


#### POCILLOPORA 
#assumptions for POC -- transformations don't help
qqnorm(poc$growth.SA)
qqline(poc$growth.SA)
leveneTest(poc$growth.SA~poc$Nutrient.Level)

poc.growth1<-lmer(growth.SA~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.growth1)
anova(poc.growth1,type=2)

# Type II Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 0.80065 0.20016     4 84.861  8.4046 9.238e-06 ***
---
  #calculate error sums of squares
  sum(resid(poc.growth1)^2)

#extract random effects -- colony p=0.2
ranova(poc.growth1)

#check residuals 
plot(poc.growth1)
hist(resid(poc.growth1))

#specific contrasts for each level - Tukey
emm.growthSA = emmeans(poc.growth1, ~ Nutrient.Level)
pairs(emm.growthSA)

# N0 - N1  -0.216847486 0.04828380 85.15  -4.491  0.0002
# N1 - N2   0.191036602 0.04828011 85.15   3.957  0.0014
# N1 - N3   0.246149467 0.04895231 85.18   5.028  <.0001
# N1 - N4   0.214373685 0.05009735 86.21   4.279  0.0005

################-
################-
###### Calcification standardized to protein
################-
################-

#norm
qqnorm(sqrt(por$growth.protein)) 
qqline(sqrt(por$growth.protein))

qqnorm(sqrt(poc$growth.protein))  
qqline(sqrt(poc$growth.protein))

#homosked
leveneTest((sqrt(por$growth.protein)~por$Nutrient.Level))
leveneTest((sqrt(poc$growth.protein)~poc$Nutrient.Level))

#different tests for each species -- PROTEIN  standardized growth (mg CaCO3 cm-2 d-1)
por.growthP1<-lmer(sqrt(growth.protein)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))

summary(por.growthP1)
anova(por.growthP1,type=2)
ranova(por.growthP1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 0.95512 0.23878     4 86.201   7.142 5.144e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#calculate error sums of squares
sum(resid(por.growthP1)^2)

#check residuals 
plot(por.growthP1)
hist(resid(por.growthP1))

### contrasts 
emm = emmeans(por.growthP1, ~ Nutrient.Level)
pairs(emm)

# N0 - N2   0.244211736 0.05813969 86.18   4.200  0.0006
# N0 - N3   0.254185919 0.05896190 86.28   4.311  0.0004
# N0 - N4   0.237397071 0.05944145 86.08   3.994  0.0013

#############################-
# Pocillopora-- protein growth

poc.growthP1<-lmer(sqrt(growth.protein)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.growthP1)
anova(poc.growthP1,type=2)

# Type II Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
# Nutrient.Level 0.40414 0.10103     4 85.062  5.8216 0.000344 ***
# ---
ranova(poc.growthP1)

#calculate error sums of squares
sum(resid(poc.growthP1)^2)

#check residuals 
plot(poc.growthP1)

emm = emmeans(poc.growthP1, ~ Nutrient.Level)
pairs(emm,adjust="holm")

# N0 - N1  -0.11734932 0.04121802 85.15  -2.847  0.0428
# N1 - N3   0.15068732 0.04178867 85.17   3.606  0.0046
# N1 - N4   0.19058527 0.04277113 86.16   4.456  0.0002

##plotting growth as mean values 
# means for plotting treatments
growth_mean.tank<- ch %>% group_by(Species,Nutrient.Level,Tank.Level) %>% summarize(
                         N    = sum(!is.na(growth.SA)),
                         meanSA = mean(growth.SA, na.rm=TRUE),
                         sd.SA   = sd(growth.SA, na.rm=TRUE),
                         se.SA = sd.SA/sqrt(N),
                         meanP = mean(growth.protein, na.rm=TRUE),
                         sd.P   = sd(growth.protein, na.rm=TRUE),
                         se.P = sd.P/sqrt(N))

#take treatment means based on tank replicate means
growth_mean.treat<- growth_mean.tank %>% group_by(Species,Nutrient.Level) %>%  summarize(
                          N    = sum(!is.na(meanSA)),
                          mean.SA = mean(meanSA, na.rm=TRUE),
                          sd.SA   = sd(meanSA, na.rm=TRUE),
                          se.SA = sd.SA/sqrt(N),
                          mean.P = mean(meanP, na.rm=TRUE),
                          sd.P   = sd(meanP, na.rm=TRUE),
                          se.P = sd.P/sqrt(N))


# SA growth - Por...main text

gSA.por<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,growth.SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,growth.SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  ylab(bquote(atop('Calcification Rate','(mg '*CaCO[3]~cm^-2~d^-1*')',phantom())))+
 # geom_label(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,growth.SA,label=Colony),position=position_jitter(.1))
  scale_y_continuous(breaks=seq(0.25,3.25,0.75))+expand_limits(y = c(0.25,3.25))+xlab("")+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+
  newtheme
gSA.por


gSA.poc<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,growth.SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,growth.SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  labs(x="",y='')+ 
  #use same scale as porites
  scale_y_continuous(breaks=seq(0.25,1.5,0.3125),labels = scales::number_format(accuracy = 0.01))+
expand_limits(y = c(0.25,1.5))+
  #scale_y_continuous(breaks=seq(0.3,.7,0.1))+expand_limits(y = c(.3,.71))+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme_right
gSA.poc

## protein growth -- POR --supplemental figure
g.prot.por<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,growth.protein),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,growth.protein),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  labs(x="",y='Calcification (mg'*~ CaCO[3]~d^-1~mg~prot^-1*')')+ 
  scale_y_continuous(breaks=seq(0.25,3.25,0.75))+expand_limits(y = c(0.25,3.25))+
  #    ylim(0.5,2.0)+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
g.prot.por

g.prot.Poc<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,growth.protein),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,growth.protein),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  #labs(x="Nutrient Level",y='Calcification (mg'*~ CaCO[3]~d^-1~mg~prot^-1*')')+ 
  scale_y_continuous(breaks=seq(0.25,1.5,0.3125),labels = scales::number_format(accuracy = 0.01))+
  labs(x="",y='')+ 
  #ylim(0.3,1)+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
g.prot.Poc
#####

################-
###### Total protein
################

#comparing species baseline calcification
t.test(pocN0$Protein.SA,porN0$Protein.SA)

##assumptions 
#normality
qqnorm(por$Protein.SA) 
qqline(por$Protein.SA)

qqnorm(poc$Protein.SA) 
qqline(poc$Protein.SA)

#homosked
leveneTest((por$Protein.SA~por$Nutrient.Level))
leveneTest((poc$Protein.SA~poc$Nutrient.Level))

#different tests for each species -- PROTEIN  standardized growth (mg CaCO3 cm-2 d-1)
por.Prot1<-lmer(Protein.SA~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))

summary(por.Prot1) 
anova(por.Prot1,type=2)
ranova(por.Prot1) ## significance tests for random factors 

# Type II Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
# Nutrient.Level 0.32421 0.081052     4 85.916  0.6604 0.6212

#calculate error sums of squares
sum(resid(por.Prot1)^2)

#check residuals 
plot(por.Prot3)
hist(resid(por.Prot3))

### contrasts --- NS

# POCILLOPORA protein
poc.Prot1<-lmer(Protein.SA~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.Prot1)
anova(poc.Prot1,type=2)

# Type II Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
# Nutrient.Level 0.38847 0.097117     4 91.116  1.2249 0.3058

#calculate error sums of squares
sum(resid(poc.Prot1)^2)

#check residuals 
plot(poc.Prot1)
# protein contrasts --- NS

#plot protein and stitch together with growth metrics
### treatments for plotting
prot_means.tank<- ch %>% group_by(Species,Nutrient.Level,Tank.Level) %>%  summarize(
                        N    = sum(!is.na(Protein.SA)),
                        mean = mean(Protein.SA, na.rm=TRUE),
                        sd   = sd(Protein.SA, na.rm=TRUE),
                        se = sd/sqrt(N))

#take treatment means based on tank level replicates
prot_means.treat<-prot_means.tank %>% group_by(Species,Nutrient.Level) %>% summarize(
                        N    = sum(!is.na(mean)),
                        meanP = mean(mean, na.rm=TRUE),
                        sd   = sd(mean, na.rm=TRUE),
                        se = sd/sqrt(N))

# Plotting POR protein --supplement Fig S6
prot.por<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Protein.SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Protein.SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  labs(x="",y='Protein (mg cm -2)')+  
  theme_classic()+
  ylim(.5,2.5)+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
prot.por

# POC protein 
prot.Poc<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Protein.SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Protein.SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  #labs(x="Nutrient Level",y='Protein (mg cm -2)')+  
  labs(x="",y='')+
  theme_classic()+
  ylim(0.0,2)+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
prot.Poc


#####

################-
###### Endosymbiont density and Chl
################
## POR
qqnorm(log(por$Zoox.cm2))
qqline(log(por$Zoox.cm2))

#homosked
leveneTest(por$Zoox.cm2~por$Nutrient.Level)

por.zooxSA1<-lmer(log(Zoox.cm2)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))

summary(por.zooxSA1)
anova(por.zooxSA1,type=2)
ranova(por.zooxSA1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)   
# Nutrient.Level  2.274  0.5685     4 84.13  4.4855 0.002468 **

#calculate error sums of squares
sum(resid(por.zooxSA1)^2)

#check residuals 
plot(por.zooxSA1) #residuals look fine and transformations don't change results 

emm = emmeans(por.zooxSA1, ~ Nutrient.Level)
pairs(emm)

# N0 - N4  -0.39903750 0.1157497 84.12  -3.447  0.0077
# N3 - N4  -0.39777747 0.1158545 84.13  -3.433  0.0080

## POC

qqnorm(log(poc$Zoox.cm2))
qqline(log(poc$Zoox.cm2))

#homosked
leveneTest(log(poc$Zoox.cm2)~poc$Nutrient.Level)

poc.zooxSA1<-lmer(log(Zoox.cm2)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.zooxSA1)
anova(poc.zooxSA1,type=2)
ranova(poc.zooxSA1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
# Nutrient.Level  6.777  1.6943     4  84.5  19.077 3.283e-11 ***

#calculate error sums of squares
sum(resid(poc.zooxSA1)^2)

#check residuals 
plot(poc.zooxSA2) #residuals look fine and transformations don't change results 

#post hoc contrasts
emm = emmeans(poc.zooxSA1, ~ Nutrient.Level)
pairs(emm)

# N0 - N2  -0.47842602 0.09569822 84.16  -4.999  <.0001
# N0 - N3  -0.56681806 0.09557596 84.08  -5.931  <.0001
# N0 - N4  -0.76137915 0.09739874 84.58  -7.817  <.0001
# N1 - N2  -0.25596937 0.09484573 84.50  -2.699  0.0625 ** incidated on plot with *
# N1 - N3  -0.34436140 0.09454269 84.16  -3.642  0.0042
# N1 - N4  -0.53892250 0.09688112 85.16  -5.563  <.0001
# N2 - N4  -0.28295313 0.09820951 84.18  -2.881  0.0392


################-
################-
##### Test for effect of symbiont density on calcification 
################-

poc.zooxNut<-lmer(growth.SA~Nutrient.Level*log(Zoox.cm2)+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))
anova(poc.zooxNut) # no interaction p =0.8782
por.zooxNut<-poc.zooxNut<-lmer(growth.SA~Nutrient.Level*log(Zoox.cm2)+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))
anova(por.zooxNut) # no interaction p =0.7989


################-
################-
###### Endosymbionts per protein
################-
#calculate the zoox.prot variable
ch$Zoox.prot<-ch$Total.Zoox/ch$Total.Protein..mg.

#redefine species to include new variable
poc<-subset(ch,Species=="POC")
por<-subset(ch,Species=="POR")

## POR

#homosked
leveneTest(log(por$Zoox.prot)~por$Nutrient.Level)

por.zooxP1<-lmer(log(Zoox.prot)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch, Species=="POR"))

summary(por.zooxP1)
anova(por.zooxP1,type=2)
ranova(por.zooxP1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)  
# Nutrient.Level  1.246 0.31151     4 84.261  2.7552 0.0331 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#calculate error sums of squares
sum(resid(por.zooxP1)^2)

#check residuals 
plot(por.zooxP1) #residuals look fine and transformations don't change results 

#contrasts
emm = emmeans(por.zooxP1, ~ Nutrient.Level)
pairs(emm)

#N0 - N4  -0.31870726 0.1093237 84.13  -2.915  0.0358
#N3 - N4  -0.28929595 0.1094219 84.14  -2.644  0.0715

#### POC

#homosked
leveneTest(log(poc$Zoox.prot)~poc$Nutrient.Level)

poc.zooxP1<-lmer(log(Zoox.prot)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.zooxP1)
anova(poc.zooxP1,type=2)
ranova(poc.zooxP1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 5.9858  1.4964     4 84.695  12.403 5.516e-08 ***

#calculate error sums of squares
sum(resid(poc.zooxP1)^2)

#check residuals 
plot(poc.zooxP1) #residuals look fine and transformations don't change results 

#contrasts
emm = emmeans(poc.zooxP1, ~ Nutrient.Level)
pairs(emm)

# N0 - N2  -0.50682834 0.1115543 84.22  -4.543  0.0002
# N0 - N3  -0.56867679 0.1114036 84.10  -5.105  <.0001
# N0 - N4  -0.62527446 0.1135132 84.64  -5.508  <.0001
# N1 - N2  -0.35970528 0.1105404 84.56  -3.254  0.0138
# N1 - N3  -0.42155373 0.1101974 84.19  -3.825  0.0023
# N1 - N4  -0.47815140 0.1128917 85.27  -4.235  0.0005

### treatments means for plotting
zoox_means.tank<- ch %>% group_by(Species,Nutrient.Level,Tank.Level) %>% summarize(
                        N    = sum(!is.na(Zoox.cm2)),
                        mean.sa = mean(Zoox.cm2, na.rm=TRUE),
                        sd.sa   = sd(Zoox.cm2, na.rm=TRUE),
                        se.sa = sd.sa/sqrt(N),
                        mean.p = mean(Zoox.prot, na.rm=TRUE),
                        sd.p   = sd(Zoox.prot, na.rm=TRUE),
                        se.p = sd.p/sqrt(N))

zoox_means.treat<- zoox_means.tank %>% group_by(Species,Nutrient.Level) %>% summarize(
                         N    = sum(!is.na(mean.sa)),
                         zooxCM2 = mean(mean.sa, na.rm=TRUE),
                         sd.sa   = sd(mean.sa, na.rm=TRUE),
                         se.sa = sd.sa/sqrt(N),
                         zooxP = mean(mean.p, na.rm=TRUE),
                         sd.p   = sd(mean.p, na.rm=TRUE),
                         se.p = sd.p/sqrt(N))
### PLOT 

# POR - zoox.cm2
por.chl1<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.cm2),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.cm2),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(250000,3500000,800000),labels = scales::scientific)+expand_limits(y = c(250000, 3500000))+
  labs(x="")+  
  ylab(bquote('Cells'~(cm^-2)))+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
por.chl1

# POC zoox cm2 
poc.chl1<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.cm2),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.cm2),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(250000,1000000,187500),labels = scales::scientific)+expand_limits(y = c(250000, 1000000))+
  labs(x="",y='')+  
  #ylim(1,1.6)+ 
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme_right
poc.chl1

# POR - zoox.prot
por.chl2<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.prot),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.prot),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(250000,3500000,800000),labels = scales::scientific)+expand_limits(y = c(250000, 3500000))+
  labs(x="",y='Cells (mg prot -1)')+  
  # ylim(9e+5,1.8e6)+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
por.chl2

# POC zoox prot 
poc.chl2<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.prot),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.prot),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(250000,1500000,312500),labels = scales::scientific)+expand_limits(y = c(250000, 1500000))+
  labs(x="",y='')+  
  #ylim(1,1.6)+ 
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
poc.chl2

################-
################-
###### Total chl cm2 and chl per symbiont
################-
################-
#calculate new variables
ch$Total.chl.cm2<-ch$Chl.a.cm2+ch$Chl.c.cm2
ch$Zoox.chl<-(ch$Total.Chl.a+ch$Total.Chl.c)/ch$Total.Zoox

#redefine species to include new variable
poc<-subset(ch,Species=="POC")
por<-subset(ch,Species=="POR")

## POR

#homosked
leveneTest(log(por$Total.chl.cm2)~por$Nutrient.Level)

por.chl.SA1<-lmer(log(Total.chl.cm2)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))

summary(por.chl.SA1)
anova(por.chl.SA1,type=2)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level   5.17  1.2925     4 83.445   10.66 5.009e-07 ***

ranova(por.chl.SA1)


# Type II Analysis of Variance Table with Satterthwaite's method
#     Type II Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# Nutrient.Level 2.1816  0.5454     4 15.405  5.1621 0.007743 **

sum(resid(por.chl.SA1)^2)

#check residuals 
plot(por.chl.SA1) #residuals look fine 

emm = emmeans(por.chl.SA1, ~ Nutrient.Level) # N0 | N1,N4 (n2 at p=0.1)
pairs(emm)

# N0 - N1  -0.57076331 0.1137807 91.64  -5.016  <.0001
# N0 - N2  -0.51277788 0.1122902 91.42  -4.567  0.0001
# N0 - N3  -0.38568081 0.1169715 91.91  -3.297  0.0118
# N0 - N4  -0.69169624 0.1148788 91.26  -6.021  <.0001
# N3 - N4  -0.30601543 0.1150893 91.50  -2.659  0.0683 * indicated in plot


#POC
qqnorm(log(poc$Total.chl.cm2)) # --- #log transform improves normality and equal variance 
qqline(log(poc$Total.chl.cm2))

#homosked
leveneTest(log(poc$Total.chl.cm2)~poc$Nutrient.Level)

poc.chl.SA1<-lmer(log(Total.chl.cm2)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.chl.SA1)
anova(poc.chl.SA1,type=2)
ranova(poc.chl.SA1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 16.339  4.0848     4 85.316  21.119 3.952e-12 ***

sum(resid(poc.chl.SA1)^2)

#check residuals 
plot(poc.chl.SA1) 

emm = emmeans(poc.chl.SA1, ~ Nutrient.Level)
pairs(emm)

# N0 - N1  -0.5966287 0.1375491 84.38  -4.338  0.0004
# N0 - N2  -0.7788682 0.1411897 84.48  -5.516  <.0001
# N0 - N3  -0.9309756 0.1410398 84.33  -6.601  <.0001
# N0 - N4  -1.2391521 0.1433845 85.27  -8.642  <.0001
# N1 - N4  -0.6425234 0.1422140 86.37  -4.518  0.0002
# N2 - N4  -0.4602839 0.1448374 84.47  -3.178  0.0172



################-
################-
##### Test for effect of chlorophyll and nutrients on calcification 
################-

poc.chlNut<-lmer(growth.SA~Nutrient.Level*log(Total.chl.cm2)+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))
anova(poc.chlNut) # no interaction p =0.06
por.chlNut<-poc.zooxNut<-lmer(growth.SA~Nutrient.Level*log(Total.chl.cm2)+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POR"))
anova(por.chlNut) # no interaction p =0.50

################-
################-
###### Total chl by total protein
################-
################-

#calculate the new variable
ch$chl.prot<-ch$Total.chl/ch$Total.Protein..mg.
por$chl.prot<-por$Total.chl/por$Total.Protein..mg.
poc$chl.prot<-poc$Total.chl/poc$Total.Protein..mg.

######## chl per zoox 

#POR
qqnorm(log(por$chl.prot)) # --- #log transform improves normality and equal variance 
qqline(log(por$chl.prot))

#homosked
leveneTest(log(por$chl.prot)~por$Nutrient.Level)

por.chl.prot1<-lmer(log(chl.prot)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=por)

summary(por.chl.prot1)
anova(por.chl.prot1,type=2)
ranova(por.chl.prot1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 3.9861 0.99654     4 83.297  5.9722 0.0002827 ***

sum(resid(por.chl.prot1)^2)

#check residuals 
plot(por.chl.prot1) #residuals look fine and transformations don't change results 

emm = emmeans(por.chl.prot1, ~ Nutrient.Level)
pairs(emm)

# N0 - N1  -0.422046976 0.1335078 83.47  -3.161  0.0181
# N0 - N2  -0.426771019 0.1317553 83.34  -3.239  0.0145
# N0 - N3  -0.379856723 0.1372677 83.65  -2.767  0.0527
# N0 - N4  -0.644674317 0.1347788 83.23  -4.783  0.0001

### POC -- chl. protein

qqnorm(log(poc$chl.prot)) # --- #log transform improves normality and equal variance 
qqline(log(poc$chl.prot))

#homosked
leveneTest(log(poc$chl.prot)~poc$Nutrient.Level)

poc.chl.prot1<-lmer(log(chl.prot)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=poc)

summary(poc.chl.prot1)
anova(poc.chl.prot1,type=2)
ranova(poc.chl.prot1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 14.461  3.6151     4 84.579  16.525 4.845e-10 ***

sum(resid(poc.chl.prot1)^2)

#check residuals 
plot(poc.chl.prot1) #residuals look fine and transformations don't change results 

emm = emmeans(poc.chl.prot1, ~ Nutrient.Level)
pairs(emm)

# N0 - N1  -0.57412911 0.1463450 84.17  -3.923  0.0016
# N0 - N2  -0.85710350 0.1501860 84.19  -5.707  <.0001
# N0 - N3  -0.83058497 0.1500028 84.10  -5.537  <.0001
# N0 - N4  -1.15611972 0.1528210 84.70  -7.565  <.0001
# N1 - N4  -0.58199061 0.1519595 85.40  -3.830  0.0022


################-
################-
###### Chl per endosymbionts
################-
################-

#POR
qqnorm(log(por$Zoox.chl)) # --- #log transform improves normality and equal variance 
qqline(log(por$Zoox.chl))

#homosked
leveneTest(log(por$Zoox.chl)~por$Nutrient.Level)

por.chl.zoox1<-lmer(log(Zoox.chl)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=por)

summary(por.chl.zoox1)
anova(por.chl.zoox1,type=2)
ranova(por.chl.zoox1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level  2.324 0.58101     4 81.375  5.4664 0.0005993 ***

sum(resid(por.chl.zoox1)^2)

#check residuals 
plot(por.chl.zoox1) #residuals look fine and transformations don't change results 

emm = emmeans(por.chl.zoox1, ~ Nutrient.Level)
pairs(emm)


# N0 - N1  -0.326978227 0.1076649 81.40  -3.037  0.0259
# N0 - N2  -0.432272143 0.1051526 81.56  -4.111  0.0009
# N0 - N3  -0.442481450 0.1113176 82.21  -3.975  0.0014
# N0 - N4  -0.322672493 0.1075793 81.40  -2.999  0.0287


###POC -- G096 was huge outlier...something wrong with zoox counts..almost none so removed from all pigment related analyses
qqnorm(log(poc$Zoox.chl)) # --- #log transform improves normality and equal variance 
qqline(log(poc$Zoox.chl))

#homosked
leveneTest(log(poc$Zoox.chl)~poc$Nutrient.Level)

poc.Zoox.chl1<-lmer(log(Zoox.chl)~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=subset(ch,Species=="POC"))

summary(poc.Zoox.chl1)
anova(poc.Zoox.chl1,type=2)
ranova(poc.Zoox.chl1)

# Type II Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
# Nutrient.Level 2.5116 0.62791     4 84.735  4.1472 0.004072 **
# ---

sum(resid(poc.Zoox.chl1)^2)

#check residuals 
plot(poc.Zoox.chl) #residuals look fine and transformations don't change results 

emm = emmeans(poc.Zoox.chl, ~ Nutrient.Level)
pairs(emm)

# N0 - N1  -0.326978227 0.1076649 81.40  -3.037  0.0259
# N0 - N2  -0.432272143 0.1051526 81.56  -4.111  0.0009
# N0 - N3  -0.442481450 0.1113176 82.21  -3.975  0.0014
# N0 - N4  -0.322672493 0.1075793 81.40  -2.999  0.0287

### Now plot zoox density and cell specific chl
chl_means.tank<- ch %>% group_by(Species,Nutrient.Level,Tank.Level) %>% summarize(
                       N    = sum(!is.na(Total.chl.cm2)),
                       mean = mean(Total.chl.cm2, na.rm=TRUE),
                       sd   = sd(Total.chl.cm2, na.rm=TRUE),
                       se = sd/sqrt(N),
                       mean.p = mean(chl.prot, na.rm=TRUE),
                       sd.p   = sd(chl.prot, na.rm=TRUE),
                       se.p = sd.p/sqrt(N),
                       mean.s = mean(Zoox.chl, na.rm=TRUE),
                       sd.s   = sd(Zoox.chl, na.rm=TRUE),
                       se.s = sd.s/sqrt(N))

chl_means.treat<- chl_means.tank %>% group_by(Species,Nutrient.Level) %>% summarize(
                        N    = sum(!is.na(mean)),
                        mean.T = mean(mean, na.rm=TRUE),
                        sd   = sd(mean, na.rm=TRUE),
                        se = sd/sqrt(N),
                        mean.P = mean(mean.p, na.rm=TRUE),
                        sd.p   = sd(mean.p, na.rm=TRUE),
                        se.p = sd.p/sqrt(N),
                        mean.sT = mean(mean.s, na.rm=TRUE),
                        sd.s   = sd(mean.s, na.rm=TRUE),
                        se.s = sd.s/sqrt(N))


# POR - total chl
#quartz()
por.chl3<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Total.chl.cm2),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Total.chl.cm2),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,100,25))+expand_limits(y = c(0, 100))+
  ylab(bquote(atop('Total Chlorophyll',~(µg~cm^-2),phantom())))+xlab('')+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
por.chl3


# POC total chl
poc.chl3<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Total.chl.cm2),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Total.chl.cm2),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,70,17.5))+expand_limits(y = c(0, 70))+
  labs(x="",y='')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
poc.chl3



# POR - Cell specific chl

por.chl4<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.chl*10^5),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,Zoox.chl*10^5),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(1,13,3))+expand_limits(y = c(1, 13))+
  ylab(bquote(atop('Cell Specific Chl',~(pg~cell^-1),phantom())))+xlab('')+
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
por.chl4


# POC - cell specific
poc.chl4<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.chl*10^5),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,Zoox.chl*10^5),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(1,13,3))+expand_limits(y = c(1, 13))+
  labs(x="",y='')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
poc.chl4


#POR chl by protein
## por chl by prot
por.chl5<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,chl.prot),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POR"),mapping=aes(Nutrient.Level,chl.prot),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,20,5))+expand_limits(y = c(0, 20))+
  labs(x="",y='Total Chl (µg mg Prot-1)')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
por.chl5


# POC
poc.chl5<-ggplot()+
  geom_boxplot(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,chl.prot),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(ch,Species=="POC"),mapping=aes(Nutrient.Level,chl.prot),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,20,5))+expand_limits(y = c(0, 20))+
  labs(x="",y='')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                           "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
poc.chl5


#####

################-
###### MAXIMUM QUANTUM YIELD
################
################-

#### analyze only final time point for PAM data to compare with other final endpoint measurements
poc.pam2<-subset(dat.clean, Species=="POC")
por.pam2<-subset(dat.clean, Species=="POR")


##POC

#analyze per species 
qqnorm(poc.pam2$Yield) # --- #log and squareroot transform doesn't change much
qqline(poc.pam2$Yield)


#homosked
leveneTest(poc.pam2$Yield~poc.pam2$Nutrient.Level)

poc.pam.mod2<-lmer(Yield~as.factor(Nutrient.Level)+(1|Tank.Level)+(1|Colony),data=poc.pam2)

summary(poc.pam.mod2)
anova(poc.pam.mod2,type=2)

# ype II Analysis of Variance Table with Satterthwaite's method
#                             Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)    
# as.factor(Nutrient.Level) 0.023174 0.0057936     4 87.358  16.029 7.061e-10 ***

#random effects p values--
ranova(poc.pam.mod2)

sum(resid(poc.pam.mod2)^2)

#check residuals 
plot(poc.pam.mod2) #residuals look fine and transformations don't change results 

emm = emmeans(poc.pam.mod2, ~ Nutrient.Level) #treatment 0 different from all 1 is diff at p=0.1 from 3 and 4
pairs(emm)

# N0 - N1  -0.022058239 0.005950443 85.19  -3.707  0.0033
# N0 - N2  -0.029300000 0.006021394 85.26  -4.866  <.0001
# N0 - N3  -0.044345246 0.006099098 85.15  -7.271  <.0001
# N0 - N4  -0.038475219 0.006212322 85.42  -6.193  <.0001
# N1 - N3  -0.022287006 0.006032762 85.20  -3.694  0.0035
# N1 - N4  -0.016416980 0.006177205 85.78  -2.658  0.0690 *marginally different

### POR

qqnorm(por.pam2$Yield) # --- #log transform doesn't change much
qqline(por.pam2$Yield)

#homosked
leveneTest(por.pam2$Yield~por.pam2$Nutrient.Level)

por.pam.mod2<-lmer(Yield~Nutrient.Level+(1|Tank.Level)+(1|Colony),data=por.pam2)

summary(por.pam.mod2)
anova(por.pam.mod2,type=2)

# # Type II Analysis of Variance Table with Satterthwaite's method
#                  Sum Sq   Mean Sq NumDF  DenDF F value    Pr(>F)    
# Nutrient.Level 0.013201 0.0033002     4 86.432  15.037 2.285e-09 ***

ranova(por.pam.mod2)
sum(resid(por.pam.mod2)^2)

#check residuals 
plot(por.pam.mod2) #residuals look fine and transformations don't change results 

emm = emmeans(por.pam.mod2, ~ Nutrient.Level) # treatment 0 different from all
pairs(emm)

# 0 - 1    -0.025913024 0.004243836 113.04  -6.106  <.0001
# 0 - 2    -0.026883767 0.003931027 113.12  -6.839  <.0001
# 0 - 3    -0.022246358 0.004243836 113.04  -5.242  <.0001
# 0 - 4    -0.036605306 0.003931024 113.12  -9.312  <.0001


## take treatment means for PAM data
pam_means.tank<- dat.clean %>% group_by(Species,Nutrient.Level,Tank.Level) %>%  summarize(
                       N    = sum(!is.na(Yield)),
                       FvFm = mean(Yield, na.rm=TRUE),
                       sd   = sd(Yield, na.rm=TRUE),
                       se = sd/sqrt(N))

pam_means.treat<- pam_means.tank %>% group_by(Species,Nutrient.Level) %>% summarize(
                        N    = sum(!is.na(FvFm)),
                        Yield = mean(FvFm, na.rm=TRUE),
                        sd   = sd(FvFm, na.rm=TRUE),
                        se = sd/sqrt(N))

##pam plots
## set lables 
treats<-c("N0","N1","N2","N3","N4")

por.pam.plot<-ggplot()+
  geom_boxplot(data=subset(dat.clean,Species=="POR"),mapping=aes(Nutrient.Level,Yield),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(dat.clean,Species=="POR"),mapping=aes(Nutrient.Level,Yield),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0.60,0.77,0.0425),labels = scales::number_format(accuracy = 0.01))+
  expand_limits(y = c(0.60,0.77 ))+
  labs(x="",y='Fv/Fm')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
  por.pam.plot


# POC pam
poc.pam.plot<-ggplot()+
  geom_boxplot(data=subset(dat.clean,Species=="POC"),mapping=aes(Nutrient.Level,Yield),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(dat.clean,Species=="POC"),mapping=aes(Nutrient.Level,Yield),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0.60,0.77,0.0425),labels = scales::number_format(accuracy = 0.01))+
  expand_limits(y = c(0.60,0.77 ))+
  labs(x="",y='')+  
  scale_x_discrete(labels=c("N0" = "0.1", "N1" = "1.0",
                            "N2" = "3.0","N3" = "5.0", "N4" = "7.0"))+newtheme
poc.pam.plot
#####

################-
###### RESPIROMETRY - photosynthesis and respiration 
################
################-

s=read.csv(file="CHAIN_Slopes.csv",stringsAsFactors = F)
lu=read.csv(file="CHAIN_respirometry_chamberID_mid_final.csv",stringsAsFactors = F)
no=read.csv(file="TPAIN_SA_prot_Dryad.csv",stringsAsFactors = F)
bw=read.csv(file="TP_buoyant_weights_Dryad.csv",stringsAsFactors = F)

lu$ChamberID=paste0("W",lu$Week,".R",lu$Round,".C",formatC(lu$Chamber,width=2,flag=0))
sl=merge(s,lu,by="ChamberID")
sl=sl[,-which(names(sl)%in%c("X.1","TPw"))]
sl$IS_BLANK=toupper(substr(sl$T.PAIN.NUBBIN.ID,1,5))=="BLANK"

#Reset Poor Fits in Coral Chambers to min of good fits (non-blank) /2
#Reset Poor Fits in Blanks to min of good fit Blanks /2
sl$LightM_gf=sl$LightM
sl$LightM_gf[sl$GoodLight.==F]=min(sl$LightM[sl$GoodLight.==T&sl$IS_BLANK==F],na.rm = T)/2
sl$LightM_gf[sl$GoodLight.==F&sl$IS_BLANK==T]=min(sl$LightM[sl$GoodLight.==T&sl$IS_BLANK==T],na.rm = T)/2
sl$DarkM_gf=sl$DarkM
sl$DarkM_gf[sl$Good.==F]=max(sl$DarkM[sl$Good.==T&sl$IS_BLANK==F],na.rm = T)/2
sl$DarkM_gf[sl$Good.==F&sl$IS_BLANK==T]=max(sl$DarkM[sl$Good.==T&sl$IS_BLANK==T],na.rm = T)/2

#Move Ahead with Dark/Light _gf (good fit) as standard
names(no)[which(names(no)=="Nubin")]="T.PAIN.NUBBIN.ID"
sln=merge(sl,no,by="T.PAIN.NUBBIN.ID")
head(sln)
ConcLU=c(.1,3,7)
names(ConcLU)=c("N0","N2","N4")
sln$Level=factor(paste0("N",sln$Level))
sln$TargetNutConc=ConcLU[sln$Level]


# Calculate Chamber Volume/Water Volume 
ChamberVol=c(64.002,
             67.149,
             69.062,
             65.629,
             65.463,
             64.948,
             68.644,
             65.385,
             66.991,
             65.723,
             66.772,
             66.854)
sln$ChamberVol=ChamberVol[sln$Chamber]

#Get Bouyant Weights
head(bw)
sln$Initial_BWg=bw$Initial[match(sln$T.PAIN.NUBBIN.ID,bw$Nub_ID)]
sln$Final_BWg=bw$Final[match(sln$T.PAIN.NUBBIN.ID,bw$Nub_ID)]

#Everything you need to Calculate the Water Volume Less Coral Volume
PocBW2Vol=function(bw) return(1.3069+0.72*bw) #
PorBW2Vol=function(bw) return(-4.1754+2.4235*bw)

#Get Coral Volume from bouyant weight fit to those other nubbins
sln$CoralVol=0
sln$CoralVol[sln$Species=="Por"]=PorBW2Vol(sln$Final_BWg[sln$Species=="Por"])
sln$CoralVol[sln$Species=="Poc"]=PocBW2Vol(sln$Final_BWg[sln$Species=="Poc"])
#Vcoralh2ovol = Vchamber - Vcoral
sln$CoralH20Vol=sln$ChamberVol-sln$CoralVol
#Vw = Vchamber - Vcoral
#Rate = M*Vw - Mblank*Vw: This is the left half, blank below
sln$LightM_gfvc=sln$LightM_gf*sln$CoralH20Vol
sln$DarkM_gfvc=sln$DarkM_gf*sln$CoralH20Vol

#Get Blank
bl=sl[which(sl$IS_BLANK),]
bl$ChamberVol=ChamberVol[bl$Chamber]
bl$LightM_ml=bl$LightM_gf*bl$ChamberVol
bl$DarkM_ml=bl$DarkM_gf*bl$ChamberVol
bl.mn=bl %>% group_by(Week,Round) %>% summarize(
            LightM_Blank=mean(LightM_gf),
            LightGOOD_Blank=all(GoodLight.),
            DarkM_Blank=mean(DarkM_gf),
            DarkGOOD_Blank=all(Good.),
            LightM_Blank_ml=mean(LightM_ml),
            DarkM_Blank_ml=mean(DarkM_ml))

slb=merge(sln,bl.mn,by=c("Week","Round"))
head(slb)

#Rate = M*Vw - Mblank*Vw: This adds the right half (the mean of light/dark blanks)
slb$LightMbc=slb$LightM_gfvc-slb$LightM_Blank_ml
slb$DarkMbc=slb$DarkM_gfvc-slb$DarkM_Blank_ml

#Given Light/Dark - Mbc move to models
library(emmeans)
library(reghelper)

#Set Up Plot text
PLlu=c("N_0.1","N_3.0","N_7.0")
names(PLlu)=c("N0","N2","N4")
slb$PlotLevel=PLlu[slb$Level]
slb$PlotSpecies="Pocillopora_acuta"  
slb$PlotSpecies[slb$Species=="Por"]="Porites_compressa"  
slb$Colony=factor(slb$Colony)

#Set up Plot Variables
slb$Pnet=slb$LightMbc
slb$R=slb$DarkMbc
slb$Pg=slb$Pnet-slb$R
slb$PR=slb$Pg/slb$R
slb$l2PR=log2(slb$PR)

#Calculate per SA metrics
slb$Pnet_SA=slb$LightMbc/(slb$SA..cm2.)
slb$R_SA=slb$DarkMbc/(slb$SA..cm2.)
slb$Pg_SA=slb$Pnet_SA-slb$R_SA
slb$PR_SA=abs(slb$Pg_SA/slb$R_SA)
slb$l2PR_SA=log2(slb$PR_SA)

#Drop Por6, Only Look at Week 4, Break Species Apart
#slb.=subset(slb,Week==4&!(Species=="Por"&Colony==6))
slb4=subset(slb,Week==4)
slb.Poc=subset(slb4,Species=="Poc")
#slb.Por.=subset(slb.,Species=="Por")
slb.Por=subset(slb4,Species=="Por")

#console_output
#sink(file = paste0(inout_path,"CHAIN_respirometry_console.txt"),append = F,type = c("output"))

### calculate liner models for each parameter per species
# Pnet = net photosynthesis
# Pg = gross photosynthesis
# R = respiration
# PR = Pg/R

# POCILLOPORA
#Leveled - Pg
Pg_Poc=lmer(Pg_SA~Level+(1|Colony)+(1|Tank),data=slb.Poc)
summary(Pg_Poc)
anova(Pg_Poc, type=2)
emmeans(Pg_Poc,~Level)
sum(resid(Pg_Poc)^2)

Pnet_Poc=lmer(Pnet_SA~Level+(1|Colony)+(1|Tank),data=slb.Poc)
summary(Pnet_Poc)
anova(Pnet_Poc, type=2)
emmeans(Pnet_Poc,~Level)
sum(resid(Pnet_Poc)^2)

#Leveled - R
R_Poc=lmer(R_SA~Level+(1|Colony)+(1|Tank),data=slb.Poc)
summary(R_Poc)
anova(R_Poc,type=2)
emmeans(R_Poc,~Level)
sum(resid(R_Poc)^2)

PR_Poc=lmer(PR_SA~Level+(1|Colony)+(1|Tank),data=slb.Poc)
summary(PR_Poc)
anova(PR_Poc,type=2)
emmeans(PR_Poc,~Level)
sum(resid(PR_Poc)^2)

# PORITES
Pg_Por=lmer(Pg_SA~Level+(1|Colony)+(1|Tank),data=slb.Por)
summary(Pg_Por)
anova(Pg_Por, type=2)
emmeans(Pg_Por,~Level)
sum(resid(Pg_Por)^2)

Pnet_Por=lmer(Pnet_SA~Level+(1|Colony)+(1|Tank),data=slb.Por)
summary(Pnet_Por)
anova(Pnet_Por, type=2)
emmeans(Pnet_Por,~Level)
sum(resid(Pnet_Por)^2)

R_Por=lmer(R_SA~Level+(1|Colony)+(1|Tank),data=slb.Por)
summary(R_Por)
anova(R_Por, type=2)
emmeans(R_Por,~Level)
sum(resid(R_Por)^2)

#Leveled - PR
PR_Por=lmer(PR_SA~Level+(1|Colony)+(1|Tank),data=slb.Por)
summary(PR_Por)
anova(PR_Por, type=2)
emmeans(PR_Por,~Level)
sum(resid(PR_Por)^2)


sink()

#mod=lm(Pnet_SA~Level*Species,data=slb4)
slb.Level=slb4 %>% group_by(PlotSpecies,Level,PlotLevel) %>% summarize(
                Pnet_mn=mean(Pnet_SA,na.rm=T),
                Pnet_se=sd(Pnet_SA,na.rm=T)/sqrt(length(which(!is.na(Pnet_SA)))),
                Pg_mn=mean(Pg_SA,na.rm=T),
                Pg_se=sd(Pg_SA,na.rm=T)/sqrt(length(which(!is.na(Pg_SA)))),
                R_mn=mean(R_SA,na.rm=T),
                R_se=sd(R_SA,na.rm=T)/sqrt(length(which(!is.na(R_SA)))),
                PR_mn=mean(PR_SA,na.rm=T),
                PR_se=sd(PR_SA,na.rm=T)/sqrt(length(which(!is.na(PR_SA)))))
slb.Level$SigLetters_Pnet=c("a","ab","b*","a","a","a")
slb.Level$SigLetters_R=c("a","a","a","a","a","a")
slb.Level$SigLetters_Pg=c("a","a","a","a","a","a")
slb.Level$SigLetters_PR=c("a","a","a","a","a","a")

#############-
############# create net photosynthesis plot for main text figure 
#quartz()
por.pnet<-ggplot()+
  geom_boxplot(data=subset(slb4,PlotSpecies=="Porites_compressa"),mapping=aes(PlotLevel,Pnet_SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(slb4,PlotSpecies=="Porites_compressa"),mapping=aes(PlotLevel,Pnet_SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,200,50))+expand_limits(y = c(0,200))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  ylab(bquote(atop('Oxygen Production',~(µmol~cm^-2),phantom())))+
  scale_x_discrete(labels=c("N_0.1" = "0.1", "N_3.0" = "3.0","N_7.0" = "7.0"))+newtheme
por.pnet

poc.pnet<-ggplot()+
  geom_boxplot(data=subset(slb4,PlotSpecies=="Pocillopora_acuta"),mapping=aes(PlotLevel,Pnet_SA),fill='#328FA9',outlier.shape=NA)+
  geom_point(data=subset(slb4,PlotSpecies=="Pocillopora_acuta"),mapping=aes(PlotLevel,Pnet_SA),pch=21, #add raw data points
             fill="gray",col='black',size=2,alpha=.55,position = position_jitter(.2))+
  scale_y_continuous(breaks=seq(0,100,25))+expand_limits(y = c(0,100))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  ylab(bquote(atop('Oxygen Production',~(µmol~cm^-2),phantom())))+
  scale_x_discrete(labels=c("N_0.1" = "0.1", "N_3.0" = "3.0","N_7.0" = "7.0"))+newtheme
poc.pnet

########## plot Gross P and Respiration for supplemental materials - Fig. S7
GP_R=ggplot(slb.Level,aes(PlotLevel,Pg_mn,group=1))+
  geom_point(size=3)+
  geom_line()+
  geom_errorbar(aes(ymin=Pg_mn-Pg_se,ymax=Pg_mn+Pg_se),width=0,stat="Identity",position=position_dodge())+
  geom_text(data=slb.Level,aes(x=PlotLevel,y=Pg_mn+10,label=SigLetters_Pg),size=4)+
  #add respiration
  geom_line(data=slb.Level,mapping=aes(x=PlotLevel,y=R_mn),linetype=2)+
  geom_point(data=slb.Level,mapping=aes(x=PlotLevel,y=R_mn),size=3,pch=21)+
  geom_errorbar(aes(ymin=R_mn-R_se,ymax=R_mn+R_se),width=0,stat="Identity",position=position_dodge())+
  geom_text(data=slb.Level,aes(x=PlotLevel,y=R_mn+10,label=SigLetters_R),size=4)+
  scale_x_discrete(labels=c("N_0.1" = "0.1", "N_3.0" = "3.0","N_7.0" = "7.0"))+
  scale_y_continuous(breaks=seq(-30,130,30))+expand_limits(y = c(-30,130))+
  ylab(bquote(atop('Change in Oxygen Concentration',~(µmol~cm^-2),phantom())))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+facet_wrap(~PlotSpecies)+theme_bw()+
  geom_hline(yintercept = 0)
GP_R

########## plot P/R for supplemental materials - Fig S8

slb.Level$PlotSpecies<-ordered(slb.Level$PlotSpecies,levels=c("Pocillopora_acuta","Porites_compressa"),
                               labels=c("Pocillopora acuta","Porites Compressa"))


PR<-ggplot(slb.Level,aes(PlotLevel,PR_mn,group=1))+
  geom_point(size=3)+
  geom_line()+
  geom_errorbar(aes(ymin=PR_mn-PR_se,ymax=PR_mn+PR_se),width=0,stat="Identity",position=position_dodge())+
  geom_text(data=slb.Level,aes(x=PlotLevel,y=PR_mn+10,label=SigLetters_PR),size=4)+
  scale_x_discrete(labels=c("N_0.1" = "0.1", "N_3.0" = "3.0","N_7.0" = "7.0"))+
  scale_y_continuous(breaks=seq(5,30,5))+expand_limits(y = c(5,30))+
  ylab("Photosynthsis:Respiration")+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+facet_wrap(~PlotSpecies)+theme_bw()
PR

#####

################-
###### Combine all plots into main text fig.2 and Figs S6,7,8
################

#### combine plots for main text (SA normalized) and supplement (prot normalized) -- cell specific chl not plotted

#growth,sym density, chl, pam, pnet  --- export 500x666
quartz()

Fig.2<-plot_grid(gSA.por,gSA.poc,por.chl1,poc.chl1,por.chl3,poc.chl3,por.chl4,poc.chl4,por.pam.plot,poc.pam.plot,por.pnet,poc.pnet,
                labels=c("A","G","B","H","C","I","D","J","E","K","F","L"),ncol=2,align="v",axis='2')
print(Fig.2)


# ggsave(plot = Fig.2,filename = "Fig_2_boxplots.pdf",
#        width=9,height=13,dpi=300,useDingbats=FALSE)

# now plot all metrics normalized to protein for supplemnt 
#put them together
Fig.S6<-plot_grid(g.prot.por,g.prot.Poc,prot.por,prot.Poc,por.chl2,poc.chl2,por.chl5,poc.chl5,
                labels=c("A","E","B","F","C","G","D","H"),ncol=2,align="v",axis='2')
quartz()
print(Fig.S6) 

# ggsave(plot = Fig.S6,filename = "Fig_S6_protein_boxplots.pdf",
#        width=9,height=13,dpi=300,useDingbats=FALSE)

#####

################-
###### MULTIVARIATE ANALYSIS - main text Fig. 3
################


set.seed(3)
library(pairwiseAdonis) #for pairwise tests of permanova
library(vegan)

#parse out only the variables of interst -- standardized to SA
nmd.cm2<-dat.clean %>%
  dplyr::select(Nub_ID,Colony,Tank.Level,Nutrient.Level,Assignment,Species,growth.SA,Protein.SA,
                Zoox.cm2,Chl.a.cm2,Yield)

#subset out POC and remove unnecessary columns -- keep only nutrient level and response variables 
poc.nmd<-nmd.cm2 %>%
  filter(Species=="POC",Assignment!="TP")%>%
  mutate(cell.chl = Chl.a.cm2/Zoox.cm2)%>% #calculate cell specific chl --- excluded in final analyses because is just chl divided by zoo...seems confounded
  dplyr::select(-c(Species,Assignment))                                           # Also it creates more noise in the data, orthogonal to zoox and doesn't
                                                                                  # vary significantly across enriched treatments in either species

poc.nmd<-na.omit(poc.nmd) #remove any nubbins with missing data
#poc.nmd<-arrange(poc.nmd,Nutrient.Level)
pocM<-as.matrix(poc.nmd[,-c(1:4,10)]) # exclude cell specific (col 10)
rownames(pocM)<-poc.nmd$Nutrient.Level

pocNMDS=metaMDS(pocM,k=3,trymax=1000)

#extract significant variable respnonses for plotting vectors
vars.poc<-envfit(pocNMDS, pocM, permu=999)
scores(vars.poc, "vectors")

#stressplot 
stressplot(pocNMDS)
#extract stress
pocNMDS$stress

#multivariate tests - adonis is like permnova for bray curtis distance matricies 
poc.data<-poc.nmd[,-c(1:4,10)] #remove nutrient factor and cell specific chl & cell specific chl

#permanova as fxn of nutrient enrichment
poc.Permanova2<-adonis2(poc.data~Nutrient.Level, data=poc.nmd,permutations=1000, 
                        method="bray", sqrt.dist = TRUE)
poc.Permanova2

#contrasts with pairwise.adonis
poc.test<-pairwise.adonis(poc.nmd[,5:9],poc.nmd$Nutrient.Level)
poc.test  

### Do the same for Porites
por.nmd<-nmd.cm2 %>%
  filter(Species=="POR",Assignment!="TP")%>%
  mutate(cell.chl = Chl.a.cm2/Zoox.cm2)%>% #calculate cell specific chl
  dplyr::select(-c(Species,Assignment))

#convert to a matrix with rownames as treatments
por.nmd<-na.omit(por.nmd)
por.nmd<-arrange(por.nmd,Nutrient.Level)
porM<-as.matrix(por.nmd[,-c(1:4,10)]) #exclude cell specific chl
rownames(porM)<-por.nmd$Nutrient.Level
porNMDS=metaMDS(porM,k=3,trymax=1000)


#extract significant variable respnonses for plotting vectors
vars.por<-envfit(porNMDS, porM, permu=999)
scores(vars.por, "vectors")

#extract stress
porNMDS$stress
stressplot(porNMDS)

###### Pairwise tests
por.data<-por.nmd[,-c(1:4,10)] #remove nutrient factor &cell specific chl

por.Permanova2<-adonis2(por.data~Nutrient.Level, data=por.nmd,permutations=1000, 
                        method="bray", sqrt.dist = TRUE)
por.Permanova2

#pairwise adonis
por.test<-pairwise.adonis(por.nmd[,5:9],por.nmd$Nutrient.Level)
por.test

#flip por colors so that the green(n2) circle is more visible
poc.nmd$Nutrient.Level<-ordered(poc.nmd$Nutrient.Level,levels=c("N4","N3","N2","N1","N0","T0"))
por.nmd$Nutrient.Level<-ordered(por.nmd$Nutrient.Level,levels=c("N4","N3","N2","N1","N0","T0"))

#make colors to match nutrient treatments
plot.col<-c("firebrick3","darkorange2","lightseagreen","slateblue3","navy")
plot.shp<-rev(c(21,22,23,24,25))

#put to figs together -- Figure. 3
par(mfcol=c(1,2), mar=c(4,4,1,1), pty="sq")

#por
por.nmds.plot<-ordiplot(porNMDS,type="n",main=substitute(paste(italic("Porites compressa"))),
                        xlim=c(-0.25, 0.25), ylim=c(-0.15, 0.15), cex.lab=0.8, cex.axis=0.8, ylab= "NMDS2", xlab="NMDS1") #plot space
axis(side = 1, labels = FALSE, tck = -0.01)
axis(side = 2, labels = FALSE, tck = -0.01)
with(por.nmd,points(porNMDS, "sites", pch=plot.shp[Nutrient.Level],cex=1.2, col="black", bg=plot.col[Nutrient.Level]))
ordiellipse(porNMDS, groups=por.nmd$Nutrient.Level, draw="polygon", col=plot.col,
            label=F, alpha=185, kind="se",conf=0.95)
par.new=T
plot(vars.por, p.max=0.05, cex=0.8, lwd=1,col='black')


#poc
ordiplot(pocNMDS,type="n",main=substitute(paste(italic("Pocillopora acuta"))),
         xlim=c(-0.25,0.25),ylim=c(-0.15, 0.15), cex.lab=0.8, cex.axis=0.8, ylab= "NMDS2", xlab="NMDS1") #plot space
#orditorp(pocNMDS,display="sites",cex=1.25,air=0.01)#
with(poc.nmd,points(pocNMDS, "sites", pch=plot.shp[Nutrient.Level],cex=1.2, col="black", bg=plot.col[Nutrient.Level]))
ordiellipse(pocNMDS, groups=poc.nmd$Nutrient.Level, draw="polygon", col=plot.col,
            label=F, alpha=185, kind="se", conf=0.95)
plot(vars.poc, p.max=0.05, cex=0.8, lwd=1,col='black')

dev.copy(pdf, "NMDS_shapes_8.11.20.pdf", height=7, width=8, useDingbats=FALSE)
dev.off()

