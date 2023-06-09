# Cody by Theadora Block
# September 10, 2020
# Statistical analysis for:
# A migratory sparrow has personality in winter that is independent of other traits 
# DOI: https://doi.org/10.1016/j.anbehav.2021.06.017
# Corresponding author: Theadora Block, tablock@ucsc.edu 


#load packages
library(tidyverse) # for dplyr, etc 
library(rptR) #for running repeatability
library(lubridate) #fixing any dates
library(car)   #package containing scatterplotMatrix, vif 
library(ggfortify) #autoplot
library(lme4)  #use to build models 
library(lmerTest) #to give p-values on lmer mod objects
library(factoextra) #for displaying PCA info in ggplot2 
library(sjPlot) #plotting model checking
library(performance) #for checking model fit with conditional (fixed + random effects) & marginal r2 (fixed effects)
library(beepr) # beeps to let you know when running code is done :)

#Notes for bootstrapping for all repeatability analysis:
set.seed(9) #make sure set.seed = 9 for all final analysis so bootstraping results are consistent with original analysis
# Currently set nboot = 0 & npermut =0, for original analysis set to = 1000 and set.seed(9)
# If set nboot & npermut to 1000 to get CI, code will take multiple hours to run

#Read in dataset
# Data for all years: 2014, 2015, 2016
#Note:this includes duplicated trials in 2014, which are in rows 1:28

setwd()
d = read.csv("GlobalPersonality&Traits_w2014reps&dupPCA.csv") 
d = d[,-1]
d$date = mdy(d$date) #make sure date is correct format 
glimpse(d)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Within-year Repeatability: for all 2014 behavioral trials ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set up data, take out all 2014 trials 
r14 = filter(d, trialyr == 2014)

#PERCHES repeatability ~~~~~~~~~~~~~~
perch_rep = rpt((perches)~(1|band),grname=c("band"),data=r14,nboot=0,npermut=0,datatype="Poisson")
perch_rep 

#QUADRANTS repeatability ~~~~~~~~~~~
quad_rep=rpt(quadrants~(1|band),grname=c("band"),data=r14,nboot=0,npermut=0,datatype="Poisson")
quad_rep   

#LATENCY to touchdown repeatability ~~~~~~~~~~~~~~
lat_mod1=rpt(latency.touchdown~(1|band),grname=c("band"),data=r14,nboot=0,npermut=0,datatype="Poisson")
lat_mod1

# VOCALIZATION repeatability ~~~~~~~~~~~~~~~~~~~
v1_rep = rpt(vocal~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
v1_rep

#PERCH TURNS repeatability~~~~~~~~~~ 
perchturn_rep= rpt(perch.turns~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
perchturn_rep  

#PERCH BOUTS repeatability~~~~~~~~~~ 
perchbout_mod = rpt(perch.bouts~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
perchbout_mod

#FLIGHTS repeatability~~~~~~~~~~~~~
flights_mod = rpt(flights~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
flights_mod 

#NUMBER OF HOPS repeatability~~~~~~
hops_mod = rpt(hops~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
hops_mod

#BAG TEST repeatability~~~~~~~~~~~~~
bag_mod = rpt(bag~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
bag_mod

#ESCAPE TEST repeatbility~~~~~~~~~~~~~
escape_mod = rpt(escape~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Poisson")
escape_mod

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Global Principal Component Analysis ####     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# use Top 6 repeatable behaviors (any behavior >25% repeatability): Hops, (27.7%), Vocalizations (34%), Perch turns (41%) Perch bouts (44.6%), Escape test (53.9%), Bag test (60.6%)

datp = d[29:277,1:18] #Make sure to take out the 2014 duplicated trials because these would bias the PCA, do global PCA with single trials only for all years 
datp = na.omit(datp)
datp = datp %>% arrange(trialyr, band)
nrow(datp)
pc = prcomp(datp[,c('hops','vocal','perch.turns','perch.bouts','escape','bag')], scale = TRUE, center = TRUE)
eig = get_eig(pc) #get eigenvalues

loading = pc$rotation #get PC factor loadings 

ind_score = pc$x #get individual PC scores 
rownames(ind_score) = paste(datp$band,datp$date,sep = '_')

#Using factoextra to plot a biplot: 
PC_plot = fviz_pca_biplot(pc, label ="var", alpha.ind="contrib", title = "Global PCA for personality",col.var = 4, ylim = c(-4, 3))
PC_plot = PC_plot + xlab('PC 1 (38.95%)') + ylab('PC 2 (19.28%)')+
  theme_classic()+
  theme(plot.title = element_text(hjust=0.5))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Predict PCA scores for 2014 duplicated trials
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get duplicated trials only so can use global dataset PC scores to predict the scores for the duplicated 2014 trials 
dupes = d[1:28,1:18] #select first 28 rows as these are the 2014 duplicate trials 

dup_pca = predict(pc, newdata = dupes[,c('hops','vocal','perch.turns','perch.bouts','escape','bag')])
dup_pca = as.data.frame(dup_pca)
rownames(dup_pca) = paste(dupes$band, dupes$date, sep = '_')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Repeatability for PCA scores for 2014 duplicated trials
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pc1_mod = rpt(pc1~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Gaussian")
pc1_mod

pc2_mod = rpt(pc2~(1|band),grname="band",data=r14,nboot=0,npermut=0,datatype="Gaussian")
pc2_mod

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Across-Year Repeatability 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set up datasets for comparing year contrasts: 2014-2015, 2015-2016, 2014-2015-2016, 2014 & 2016
dat = d[29:277,] #Take out 2014 duplicate trials
dat1 = filter(dat, trialyr != '2016') #has 2014 & 2015
dat2 = filter(dat, trialyr != '2014') # has 2015 & 2016
dat3 = filter(dat, trialyr != '2015') #has 2014 & 2016

#2014-2015, using dat1 ~~~~~~~~~~~~~~~
br1 = rpt(bag ~ 1|band, grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
br1$ngroups
br1

er1 = rpt(escape ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
er1$ngroups
er1  

tr1 = rpt(perch.turns ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
tr1$ngroups
tr1

pbr1 = rpt(perch.bouts ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
pbr1$ngroups
pbr1

hr1 = rpt(hops ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
hr1$ngroups
hr1

vr1 = rpt(vocal ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Poisson")
vr1$ngroups
vr1

pc1.r1 = rpt((pc1) ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Gaussian")
pc1.r1$ngroups
pc1.r1

pc2.r1 = rpt((pc2) ~ (1|band), grname = c('band'),data = dat1, nboot = 0, npermut = 0, datatype = "Gaussian")
pc2.r1$ngroups
pc2.r1

#2015-2016 ~~~~~~~~~~~~~~~
br2 = rpt(bag ~ 1|band, grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
print("Repeatability measured for bag")
br2$ngroups
br2 

er2 = rpt(escape ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
er2$ngroups
er2

tr2 = rpt(perch.turns ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
tr2$ngroups
tr2

pbr2 = rpt(perch.bouts ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
pbr2$ngroups
pbr2

hr2 = rpt(hops ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
hr2$ngroups
hr2

vr2 = rpt(vocal ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Poisson")
vr2$ngroups
vr2

pc1.r2 = rpt(pc1 ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Gaussian")
pc1.r2$ngroups
pc1.r2

pc2.r2 = rpt(pc2 ~ (1|band), grname = c('band'),data = dat2, nboot = 0, npermut = 0, datatype = "Gaussian")
pc2.r2$ngroups
pc2.r2

#2014-2015-2016 ~~~~~~~~~~~~~~~
br = rpt(bag ~ 1|band, grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
br$ngroups
br 

er = rpt(escape ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
era$ngroups
er

tr = rpt(perch.turns ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
tr$ngroups
tr

pbr = rpt(perch.bouts ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
pbr$ngroups
pbr

hr = rpt(hops ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
hra$ngroups
hr

vr = rpt(vocal ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Poisson")
vr$ngroups
vr

pc1.r = rpt(pc1 ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Gaussian")
pc1.r$ngroups
pc1.r

pc2.r = rpt(pc2 ~ (1|band), grname = c('band'),data = dat, nboot = 0, npermut = 0, datatype = "Gaussian")
pc2.r$ngroups
pc2.r

#2014 & 2016 ~~~~~~~~~~~~~~~
br4 = rpt(bag ~ 1|band, grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
br4$ngroups
br4 

er4 = rpt(escape ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
er4$ngroups
er4 

tr4 = rpt(perch.turns ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
tr4$ngroups
tr4 

pbr4 = rpt(perch.bouts ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
pbr4$ngroups
pbr4

hr4 = rpt(hops ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
hr4$ngroups
hr4

vr4 = rpt(vocal ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Poisson")
vr4$ngroups
vr4

pc1.r4 = rpt(pc1 ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Gaussian")
pc1.r4$ngroups
pc1.r4

pc2.r4 = rpt(pc2 ~ (1|band), grname = c('band'),data = dat3, nboot = 0, npermut = 0, datatype = "Gaussian")
pc2.r4$ngroups
pc2.r4

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Testing for personality correlates 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set data up for analysis
glimpse(dat) # Single trials for all years 
dat14 = filter(dat,trialyr == '2014')
dat15 = filter(dat, trialyr == '2015')
dat16 = filter(dat, trialyr == '2016')

# Global correlations: 2014, 2015, 2016 ~~~~~~~~~~~~~~~~~~~
# PC 1 ~~~~~~~~~~~~~
a = lmer(pc1 ~ (1|band)+scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat, na.action = na.omit)
summary(a)
ap = sjPlot::plot_model(a, type = 'diag') #check residuals and model assumptions
cowplot::plot_grid(ap[[1]], ap[[3]], ap[[4]])
r2_nakagawa(a) #Marginal and conditional R2 values to see what variation the random & fixed effects accounts for. Marginal R2: what variation the fixed & random effects account for. Conditional R2: what variation the random & fixed effects accounts for
vif(a) #make sure VIF <5 to ensure low collinearity

# PC 2 ~~~~~~~~~~~~~
b = lmer(pc2 ~ (1|band)+scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat, na.action = na.omit )
summary(b)
bp = sjPlot::plot_model(b, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(bp[[1]], bp[[3]], bp[[4]])
r2_nakagawa(b)
vif(b)

# Vocalizations ~~~~~~~~~~~~
v = glmer(vocal ~ (1|band)+scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat, family = poisson, na.action = na.omit )
summary(v)
plot(v)
hist(resid(v))
r2_nakagawa(v)
vif(v)

#Personality correlates by year  ~~~~~~~~~~~~~~~~~~~
# 2014 ~~~~~~~~~~~~~~~~~
# PC 1 ~~~~~~~~~~~~~
a1 = lm(pc1 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat14, na.action = na.omit )
summary(a1) #no correlations in model
a1p = sjPlot::plot_model(a1, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(a1p[[1]],a1p[[2]], a1p[[3]], a1p[[4]]) 
vif(a1)

# PC 2 ~~~~~~~~~~~~~
b1 = lm(pc2 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat14, na.action = na.omit )
summary(b1)
b1p = sjPlot::plot_model(b1, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(b1p[[1]],b1p[[2]], b1p[[3]], b1p[[4]]) 
vif(b1)

# Vocalizations ~~~~~~~~~~~~
v1 = glm(vocal ~ scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat14, family = poisson, na.action = na.omit )
summary(v1)
r2(v1)
hist(resid(v1))
vif(v1)

# 2015 ~~~~~~~~~~~~~~~~~
# PC 1 ~~~~~~~~~~~~~
a2 = lm(pc1 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+age, data = dat15, na.action = na.omit ) #removed wing.chord because it caused strong collinearity in the model
summary(a2)
a2p = sjPlot::plot_model(a2, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(a2p[[1]],a2p[[2]], a2p[[3]], a2p[[4]]) 

vif(a2)

# PC 2 ~~~~~~~~~~~~~
b2 = lm(pc2 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+age, data = dat15, na.action = na.omit ) #removed wing.chord because it caused strong collinearity in the model
summary(b2)
b2p = sjPlot::plot_model(b2, type = 'diag', scale = 3) #check residuals and model assumptions 
cowplot::plot_grid(b2p[[1]],b2p[[2]], b2p[[3]], b2p[[4]]) 
vif(b2)

# Vocalizations ~~~~~~~~~~~~
v2 = glm(vocal ~ scale(elo)+sex+scale(black.size)+scale(gold.size)+age, data = dat15, family = poisson, na.action = na.omit ) #removed wing.chord because it caused strong collinearity in the model (In original model with all factors, wing.chord VIF = 9)
summary(v2)
r2(v2)
hist(resid(v2))

vif(v2)

# 2016 ~~~~~~~~~~~~~~~~~
# PC 1 ~~~~~~~~~~~~~
a3 = lm(pc1 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat16, na.action = na.omit )
summary(a3)
a3p = sjPlot::plot_model(a3, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(a3p[[1]],a3p[[2]], a3p[[3]], a3p[[4]])
vif(a3)

# PC 2 ~~~~~~~~~~~~~
b3 = lm(pc2 ~scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat16, na.action = na.omit )
summary(b3)
b3p = sjPlot::plot_model(b3, type = 'diag') #check residuals and model assumptions 
cowplot::plot_grid(b3p[[1]],b3p[[2]], b3p[[3]], b3p[[4]]) 
vif(b3)

# Vocalizations ~~~~~~~~~~~~
v3 = glm(vocal~scale(elo)+sex+scale(black.size)+scale(gold.size)+scale(wing.chord)+age, data = dat16, family = poisson, na.action = na.omit )
summary(v3)
r2(v3)
hist(resid(v3))
vif(v3)

#### ~~~~~~~~~~~~ Figures for Personality correlations   ~~~~~~~~~~~~

# Figure for 2014, 2015, 2016 global correlations 
pa = plot_model(a, show.p = T, title = 'PC 1 (activity)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1.5,1.5))+
  labs(x = "Fixed effects")+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length','ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))

pb = plot_model(b, type = 'est', show.p = T, title = 'PC 2 (escape)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1.5,1.5))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))

pv = plot_model(v, type = 'est', show.p = T, title = 'Vocalizations', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1,10))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size*","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', colour = 'Black', size = 14), axis.text.y = element_text(face = c('plain', 'plain', 'bold', 'plain', 'plain', 'plain')))

pg = cowplot::plot_grid(pa, pb, pv, align = 'h', nrow = 1)
pg

#sjPlot::save_plot(filename = "~/Ch1_figs/Persn_corrGlobal.jpg",fig = pg, dpi = 720, height = 6, width = 22, label.color = 'Black')

# Figure for 2014 correlations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pa14 = plot_model(a1, show.p = T, title = 'PC 1 (activity)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1.5,1.5))+
  labs(x = "Fixed effects")+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length','ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))

pb14 = plot_model(b1, type = 'est', show.p = T, title = 'PC 2 (escape)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1.5,1.5))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))

pv14 = plot_model(v1, type = 'est', show.p = T, title = 'Vocalizations', vline.color = 'grey')+
  scale_y_continuous(limits=c(0,6))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', colour = 'Black', size = 14))
pv14

pg14 = cowplot::plot_grid(pa14, pb14, pv14, align = 'h', nrow = 1)
pg14

#sjPlot::save_plot(filename = "~/Persn_corr2014.jpg",fig = pg14, dpi = 720, height = 6, width = 22, label.color = 'Black')

#Figures for 2015 personality correlations  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pa15 = plot_model(a2, type = 'est', show.p = T, title = 'PC 1 (activity)', vline.color = 'grey')+
  labs(x = "Fixed effects")+
  #scale_y_continuous(limits=c(-1,1.5))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))
pa15 

pb15 = plot_model(b2, type = 'est', show.p = T, title = 'PC 2 (escape)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-1,1.2))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))
pb15

pv15 = plot_model(v2, type = 'est', show.p = T, title = 'Vocalizations')+
  scale_y_continuous(limits=c(-1,10))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size*","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', colour = 'Black', size = 14), axis.text.y = element_text(face = c('plain', 'bold', 'plain', 'plain', 'plain'))) +
  geom_hline(yintercept = 1, colour = 'black', size = 1.2, alpha = .25 )
pv15

pg15 = cowplot::plot_grid(pa15, pb15, pv15, align = 'h', nrow = 1)
pg15

#sjPlot::save_plot(filename = "~/Persn_corr2015.jpg",fig = pg15, dpi = 720, height = 6, width = 22, label.color = 'Black')

#Figures for 2016 personality correlations  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pa16 = plot_model(a3, type = 'est', show.p = T, title = 'PC 1 (activity)', vline.color = 'grey')+
  labs(x = "Fixed effects")+
  scale_y_continuous(limits=c(-5,5))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)*",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14), axis.text.y = element_text(face = c('plain', 'plain', 'plain', 'plain', 'bold','plain')))
pa16

pb16 = plot_model(b3, type = 'est', show.p = T, title = 'PC 2 (escape)', vline.color = 'grey')+
  scale_y_continuous(limits=c(-2.5,2.5))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', color = 'Black', size = 14))
pb16

pv16 = plot_model(v3, type = 'est', show.p = T, title = 'Vocalizations', vline.color = 'grey')+
  scale_y_continuous(limits=c(-.1,4))+
  scale_x_discrete(labels=c("scale(elo)"="Dominance rank",'scale(black.size)'="Black size", 'scale(gold.size)'="Gold size","sexM"="Sex (M)",'scale(wing.chord)'='Wing length**', 'ageHY' = 'Age (HY)')) +
  theme_classic()+
  theme(text = element_text(family = 'Times New Roman', colour = 'Black', size = 14), axis.text.y = element_text(face = c('plain', 'bold', 'plain', 'plain', 'plain', 'plain'))) 
pv16

pg16 = cowplot::plot_grid(pa16, pb16, pv16, align = 'h', nrow = 1)
pg16

#sjPlot::save_plot(filename = "~/Persn_corr2016.jpg",fig = pg16, dpi = 720, height = 6, width = 22, label.color = 'Black')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   
#### Testing for a time since first trial affect 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~ PC1 ~~~~~~~~~~~~~~~~~~~
rep_p1a = lmer(pc1 ~ (1|band)+ trial_time,data = dat) 
summary(rep_p1a)
plot_model(rep_p1a, type = 'diag')

plot(rep_p1a) #look at residuals
qqnorm(resid(rep_p1a)) #look at qqnorm plot
r2_nakagawa(rep_p1a) #fixed effect explains none of the variation of the model

#~~~~~~~~~~~~~~~~~~~ PC2 ~~~~~~~~~~~~~~~~~~~
rep_p2a = lmer(pc2 ~ (1|band)+ trial_time,data = dat) 
summary(rep_p2a)

plot(rep_p2a) #look at residuals       
qqnorm(resid(rep_p2a)) #look at qqnorm plot  
r2_nakagawa(rep_p2a) #fixed effect explains almost none of the variation of the model

#~~~~~~~~~~~~~~~~~~~ Vocalizations ~~~~~~~~~~~~~~~~~~~
rep_voa = glmer(vocal ~ (1|band) + trial_time,data = dat, family = 'poisson') 
summary(rep_voa)

plot(rep_voa) #look at residuals       
qqnorm(resid(rep_voa)) #look at qqnorm plot
r2_nakagawa(rep_voa) #fixed effect explains none of the variation of the model


