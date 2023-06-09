##### Strategies paper ProcB submission  - SI Figures #####
# This file makes the panels A,B,C in Figures 2 and 3 in the main manuscript and shows how the stats were done using a non linear mixed model 
# plotting and stats for panels D,E,F are in another file and use a linear mixed model. 


rm(list=ls()) #remove everything that may still be in memory 
library("papeR")
library(xtable)
library(data.table)
library(plyr)
library(tools)
library(zoo)
library(lattice)
library(ggplot2)
require(nlme)
require(minpack.lm)
library(reshape2)
library(devtools)
library(MuMIn)
library(emmeans)
library(AICcmodavg)

setwd("~/Dropbox/Collins Lab Shared Folder/'quorum' - transiently restored/Proc B/Data files as put onto dryad")


#make an exponential decay function just for stats and giggles  (well actually this is what we will use in the models throughout)
exp_dec <- function(a, b, growth) {
  y <- a * exp(growth * b) # a and be being the intercept and steepness of slope 
  return(y)
}


###Figure02, panels A,B,C #####
coc<-read.csv("20190418_foldchangedata_cocu.csv")
head(coc)

#need a biorep id column for the random effect to get the nesting right
coc<-within(coc,lineage_id_full <- (factor(techrep):factor(biorep):factor(grownagainst)))
head(coc)
str(coc)

#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
coc$Foc_lin<-as.factor(coc$Foc_lin)
coc$grownagainst<-as.factor(coc$grownagainst)


#fullest global model - excludes data where both sides of the ThinCert had the same species as that had no effect 
exp.mix_coc <- nlme( Comp_co~ exp_dec(a, b, growth = Mono),
                     data = coc,
                     fixed = list(a + b ~ 1+grownagainst), # could also use 'Foc_Lin'? 
                     random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                     groups = ~ Foc_lin,
                     start = c(0.2, rep(0,5),-0.003,rep(0,5)), #5 levels minus 1
                     na.action = na.omit,
                     method = 'ML', #don't forget to refit with REML later 
                     control = nlmeControl(pnlsTol = 0.02))



# remove treatment effect on a (intercept)
exp.mix_coc1 <- update(exp.mix_coc, fixed = list(a ~ 1, b ~ 1 + grownagainst), start = c(0.2, -0.003, rep(0,5)))
anova(exp.mix_coc, exp.mix_coc1) #  matters  VERY much 
#write.csv(prettify(anova(exp.mix_coc, exp.mix_coc1)),"Anova_intercepts_coc_Co2.csv")

# remove treatment effect on b (slope)
exp.mix_coc2 <- update(exp.mix_coc, fixed = list(a ~ 1 + grownagainst, b ~ 1), start = c(0.2, rep(0,5), -0.003))
anova(exp.mix_coc, exp.mix_coc2) # yes
#use global model 
exp.mix_coc_fin <- update(exp.mix_coc, method = 'REML')
summary(exp.mix_coc_fin )
#write.csv(prettify(summary(exp.mix_coc_fin )),"Anova_intercepts_coc_Co2_OUT.csv")

### plot - generate predicted curves from fixed effects 

pred_self <- seq(min(coc$Mono), max(coc$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_coc_fin))

#so this may be a bit more tricky now, because we have to distinguish between self and other
pred_oth <- exp_dec(a = fixef(exp.mix_coc_fin)[1] , b = fixef(exp.mix_coc_fin)[7], growth = pred_self) 


pred_rcc1107 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[2], b = fixef(exp.mix_coc_fin)[7]+ fixef(exp.mix_coc_fin)[8], growth = pred_self)
pred_rcc1108 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[3], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[9], growth = pred_self)
pred_rcc1558 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[4], b = fixef(exp.mix_coc_fin)[7]+ fixef(exp.mix_coc_fin)[10], growth = pred_self)
pred_rcc410 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[5], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[11], growth = pred_self)
pred_rcc810 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[6], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[12], growth = pred_self)

#make a date frame containing all predictions
pred_all <- data.frame(MONO = rep(pred_self, times = 6), lineage = c(rep('oth95', times = length(pred_oth)), rep('rcc1107', times = length(pred_rcc1107)),rep('rcc1108', times = length(pred_rcc1108)),rep('rcc1558', times = length(pred_rcc1558)),rep('rcc410', times = length(pred_rcc410)),rep('rcc810', times = length(pred_rcc810))), CO_CULTURE = c(pred_oth, pred_rcc1107,pred_rcc1108,pred_rcc1558,pred_rcc410,pred_rcc810))
#check
str(pred_all)

## and now again with the 'good colour scheme'

CO_PLOT_PREDICTIONS3 <- ggplot(coc, aes(Mono,Comp_co , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~self~day^{-1}), y=expression(Fold~change~compared~to~no~other))+
  geom_line(aes(MONO, CO_CULTURE, colour = lineage), pred_all) 
CO_PLOT_PREDICTIONS3

#make per biorep average frame in case it is needed later or we don't want all that confetti all across the screen!! 
coc_av<-ddply(coc,c("Foc_lin","biorep"), function(df) return(c(muecomp.avg=mean(df$Comp_co),muecomp.se=sd(df$Comp_co)/sqrt(length(coc)))))
#can also just fit ONE exo_decay for easier visualisation . Also should get rid of the jitter effect somehow... 

###spikes ###

spi<-read.csv("20190418_foldchangedata_spikes.csv")
head(spi)

#need a biorep id column for the random effect to get the nesting right
spi<-within(spi,lineage_id_full <- (factor(techrep):factor(biorep):factor(grownagainst)))
head(spi)
str(spi)

#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
spi$Foc_lin<-as.factor(spi$Foc_lin)
spi$grownagainst<-as.factor(spi$grownagainst)

#fullest global model (want to exclude data where they have been self-spiked?)

exp.mix_spiked <- nlme( Comp_spike~ exp_dec(a, b, growth = Mono),
                        data = spi,
                        fixed = list(a + b ~ 1+grownagainst), 
                        random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                        groups = ~ Foc_lin,
                        start = c(1.4, rep(0,5),-0.0001,rep(0,5)), #6 levels minus 1, plus slope and intercept quite different from co culture 
                        na.action = na.omit,
                        method = 'ML', #don't forget to refit with REML later 
                        control = nlmeControl(pnlsTol = 0.02))



# remove treatment effect on a (intercept)
exp.mix_spiked1 <- update(exp.mix_spiked, fixed = list(a ~ 1, b ~ 1 + grownagainst), start = c(1.4, -0.0001, rep(0,5)))
anova(exp.mix_spiked, exp.mix_spiked1) #  matters  
#write.csv(prettify(anova(exp.mix_spiked, exp.mix_spiked1)),"Anova_intercepts_spikes_Co2.csv")


# remove treatment effect on b (slope)
exp.mix_spiked2 <- update(exp.mix_spiked, fixed = list(a ~ 1 + grownagainst, b ~ 1), start = c(1.4, rep(0,5), -0.0001))
anova(exp.mix_spiked, exp.mix_spiked2) # yes
#intercept more palatable
exp.mix_final_spiked <- update(exp.mix_spiked1, method = 'REML')
summary(exp.mix_final_spiked )


### plot - generate predicted curves from fixed effects 

pred_self <- seq(min(spi$Mono), max(spi$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_final_spiked))

#so this may be a bit more tricky now, because we have to distinguish between self and other
pred_oth <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2], growth = pred_self) 
pred_rcc1108 <- exp_dec(a = fixef(exp.mix_final_spiked)[1], b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[3], growth = pred_self) 
pred_rcc1558 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[4], growth = pred_self) 
pred_rcc343 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[5], growth = pred_self) 
pred_rcc410 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[6], growth = pred_self) 
pred_rcc809 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[7], growth = pred_self) 

#make a date frame containing all predictions
pred_all <- data.frame(SELF = rep(pred_self, times = 6), lineage = c(rep('oth95', times = length(pred_oth)), rep('rcc1108', times = length(pred_rcc1108)),rep('rcc1558', times = length(pred_rcc1558)),rep('rcc343', times = length(pred_rcc343)),rep('rcc410', times = length(pred_rcc410)),rep('rcc809', times = length(pred_rcc809))), SPIKED = c(pred_oth, pred_rcc1108,pred_rcc1558,pred_rcc343,pred_rcc410,pred_rcc809))

#check
str(pred_all)

## and now again with the 'good colour scheme'
SPIKE_PLOT_PREDICTIONS3 <- ggplot(spi, aes(Mono,Comp_spike , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~seawater~spike~µ~day^{-1}), y=expression(Fold~change~compared~to~no~spike))
SPIKE_PLOT_PREDICTIONS3


spi_av<-ddply(spi,c("Foc_lin","biorep_focal"), function(df) return(c(muecompspi.avg=mean(df$Comp_to_spike),muecomp.se=sd(df$Comp_to_spike)/sqrt(length(coc)))))

## GFP ##
GFPco<-read.csv("20190418_foldchangedata_GFP.csv")
head(GFPco)
#need a biorep id column for the random effect to get the nesting right
GFPco<-within(GFPco,lineage_id_full <- (factor(techrep):factor(biorep):factor(against)))
head(GFPco)
str(GFPco)
#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
GFPco$Foc_lin<-as.factor(GFPco$Foc_lin)
GFPco$against<-as.factor(GFPco$against)
#fullest global model (want to exclude data where they have been self-spiked?)

exp.mix_GFP <- nlme( GFP_co~ exp_dec(a, b, growth = Mono),
                     data = GFPco,
                     fixed = list(a + b ~ 1), 
                     random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                     groups = ~ Foc_lin,
                     start = c(0.2, -0.003), 
                     na.action = na.omit,
                     method = 'ML',
                     control = nlmeControl(pnlsTol = 0.02))
summary(exp.mix_GFP)

# both parameters significant, so update full model as REML
exp.mix_final_GFP <- update(exp.mix_GFP, method = 'REML')
summary(exp.mix_final_GFP)
#write.csv(prettify(summary(exp.mix_GFP)),"model_all_out_GFP.csv")


### plot, generate predicted curves from fixed effects 

pred_alone <- seq(min(GFPco$Mono), max(GFPco$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_final_GFP))

# easier here - don't have multiple self/non-self, just GFP or NOT GFP 
pred_curve <- exp_dec(a = fixef(exp.mix_final_GFP)[1] , b = fixef(exp.mix_final_GFP)[2], growth = pred_alone) 
#make a date frame containing all predictions
pred_all <- data.frame(ALONE = rep(pred_alone,8), lineage=c(rep('oth95', length(pred_alone)),rep('rcc1107', length(pred_alone)),rep('rcc1114', length(pred_alone)),rep('rcc1645', length(pred_alone)),rep('rcc410', length(pred_alone)),rep('rcc747', length(pred_alone)),rep('rcc789', length(pred_alone)),rep('rcc810', length(pred_alone))),GFP = rep(pred_curve,8))
#check
str(pred_all)


# this against the not means data frame , but we make a means frame below because we can... 

GFP_PLOT_PREDICTIONS3 <- ggplot(GFPco, aes(Mono,GFP_co , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~alone~µ~day^{-1}), y=expression(Fold~change~withGFP~compared~to~monoculture)) 
GFP_PLOT_PREDICTIONS3


GFP_av<-ddply(GFPco,c("Foc_lin","biorep"), function(df) return(c(muecompspi.avg=mean(df$GFP_co),muecomp.se=sd(df$GFP_co)/sqrt(length(GFPco)))))


####Figure03, panels A,B,C , i.e. elevated pCO2####
## direct co culture ##
coc<-read.csv("20190418_foldchangedataCO2_cocu.csv")
head(coc)

#need a biorep id column for the random effect to get the nesting right
coc<-within(coc,lineage_id_full <- (factor(techrep):factor(biorep):factor(grownagainst)))
head(coc)
str(coc)

#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
coc$Foc_lin<-as.factor(coc$Foc_lin)
coc$grownagainst<-as.factor(coc$grownagainst)


#fullest global model (want to exclude data where they have been self-spiked?)
exp.mix_coc <- nlme( Comp_co~ exp_dec(a, b, growth = Mono),
                     data = coc,
                     fixed = list(a + b ~ 1+grownagainst), 
                     random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                     groups = ~ Foc_lin,
                     start = c(0.2, rep(0,5),-0.003,rep(0,5)), #5 levels minus 1
                     na.action = na.omit,
                     method = 'ML', #don't forget to refit with REML later 
                     control = nlmeControl(pnlsTol = 0.02))



# remove treatment effect on a (intercept)
exp.mix_coc1 <- update(exp.mix_coc, fixed = list(a ~ 1, b ~ 1 + grownagainst), start = c(0.2, -0.003, rep(0,5)))
anova(exp.mix_coc, exp.mix_coc1) #  matters  VERY much 
#write.csv(prettify(anova(exp.mix_coc, exp.mix_coc1)),"Anova_intercepts_coc_Co2.csv")

# remove treatment effect on b (slope)
exp.mix_coc2 <- update(exp.mix_coc, fixed = list(a ~ 1 + grownagainst, b ~ 1), start = c(0.2, rep(0,5), -0.003))
anova(exp.mix_coc, exp.mix_coc2) # yes
#use global model 
exp.mix_coc_fin <- update(exp.mix_coc, method = 'REML')
summary(exp.mix_coc_fin )
#write.csv(prettify(summary(exp.mix_coc_fin )),"Anova_intercepts_coc_Co2_OUT.csv")

### plot - generate predicted curves from fixed effects 

pred_self <- seq(min(coc$Mono), max(coc$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_coc_fin))

#so this may be a bit more tricky now, because we have to distinguish between self and other
pred_oth <- exp_dec(a = fixef(exp.mix_coc_fin)[1] , b = fixef(exp.mix_coc_fin)[7], growth = pred_self) 


pred_rcc1107 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[2], b = fixef(exp.mix_coc_fin)[7]+ fixef(exp.mix_coc_fin)[8], growth = pred_self)
pred_rcc1108 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[3], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[9], growth = pred_self)
pred_rcc1558 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[4], b = fixef(exp.mix_coc_fin)[7]+ fixef(exp.mix_coc_fin)[10], growth = pred_self)
pred_rcc410 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[5], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[11], growth = pred_self)
pred_rcc810 <- exp_dec(a = fixef(exp.mix_coc_fin)[1] + fixef(exp.mix_coc_fin)[6], b = fixef(exp.mix_coc_fin)[7]+  fixef(exp.mix_coc_fin)[12], growth = pred_self)

#make a date frame containing all predictions
pred_all <- data.frame(MONO = rep(pred_self, times = 6), lineage = c(rep('oth95', times = length(pred_oth)), rep('rcc1107', times = length(pred_rcc1107)),rep('rcc1108', times = length(pred_rcc1108)),rep('rcc1558', times = length(pred_rcc1558)),rep('rcc410', times = length(pred_rcc410)),rep('rcc810', times = length(pred_rcc810))), CO_CULTURE = c(pred_oth, pred_rcc1107,pred_rcc1108,pred_rcc1558,pred_rcc410,pred_rcc810))
#check
str(pred_all)

## and now again with the 'good colour scheme'

CO_PLOT_PREDICTIONS3 <- ggplot(coc, aes(Mono,Comp_co , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~self~day^{-1}), y=expression(Fold~change~compared~to~no~other))+
  geom_line(aes(MONO, CO_CULTURE, colour = lineage), pred_all) 
CO_PLOT_PREDICTIONS3

#make per biorep average frame in case it is needed later or we don't want all that confetti all across the screen!! 
coc_av<-ddply(coc,c("Foc_lin","biorep"), function(df) return(c(muecomp.avg=mean(df$Comp_co),muecomp.se=sd(df$Comp_co)/sqrt(length(coc)))))


###spikes ###

spi<-read.csv("20190418_foldchangedataCO2_spikes.csv")
head(spi)

#need a biorep id column for the random effect to get the nesting right
spi<-within(spi,lineage_id_full <- (factor(techrep):factor(biorep):factor(grownagainst)))
head(spi)
str(spi)

#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
spi$Foc_lin<-as.factor(spi$Foc_lin)
spi$grownagainst<-as.factor(spi$grownagainst)

#fullest global model (want to exclude data where they have been self-spiked?)

exp.mix_spiked <- nlme( Comp_spike~ exp_dec(a, b, growth = Mono),
                        data = spi,
                        fixed = list(a + b ~ 1+grownagainst), 
                        random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                        groups = ~ Foc_lin,
                        start = c(0.2, rep(0,5),-0.003,rep(0,5)), #6 levels minus 1
                        na.action = na.omit,
                        method = 'ML', #don't forget to refit with REML later 
                        control = nlmeControl(pnlsTol = 0.02))



# remove treatment effect on a (intercept)
exp.mix_spiked1 <- update(exp.mix_spiked, fixed = list(a ~ 1, b ~ 1 + grownagainst), start = c(0.2, -0.003, rep(0,5)))
anova(exp.mix_spiked, exp.mix_spiked1) #  matters  
#write.csv(prettify(anova(exp.mix_spiked, exp.mix_spiked1)),"Anova_intercepts_spikes_Co2.csv")


# remove treatment effect on b (slope)
exp.mix_spiked2 <- update(exp.mix_spiked, fixed = list(a ~ 1 + grownagainst, b ~ 1), start = c(0.2, rep(0,5), -0.003))
anova(exp.mix_spiked, exp.mix_spiked2) # yes
#intercept more palatable
exp.mix_final_spiked <- update(exp.mix_spiked1, method = 'REML')
summary(exp.mix_final_spiked )


### plot - generate predicted curves from fixed effects 

pred_self <- seq(min(spi$Mono), max(spi$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_final_spiked))

#so this may be a bit more tricky now, because we have to distinguish between self and other
pred_oth <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2], growth = pred_self) 
pred_rcc1108 <- exp_dec(a = fixef(exp.mix_final_spiked)[1], b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[3], growth = pred_self) 
pred_rcc1558 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[4], growth = pred_self) 
pred_rcc343 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[5], growth = pred_self) 
pred_rcc410 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[6], growth = pred_self) 
pred_rcc809 <- exp_dec(a = fixef(exp.mix_final_spiked)[1] , b = fixef(exp.mix_final_spiked)[2]+fixef(exp.mix_final_spiked)[7], growth = pred_self) 

#make a date frame containing all predictions
pred_all <- data.frame(SELF = rep(pred_self, times = 6), lineage = c(rep('oth95', times = length(pred_oth)), rep('rcc1108', times = length(pred_rcc1108)),rep('rcc1558', times = length(pred_rcc1558)),rep('rcc343', times = length(pred_rcc343)),rep('rcc410', times = length(pred_rcc410)),rep('rcc809', times = length(pred_rcc809))), SPIKED = c(pred_oth, pred_rcc1108,pred_rcc1558,pred_rcc343,pred_rcc410,pred_rcc809))

#check
str(pred_all)

## and now again with the 'good colour scheme'
SPIKE_PLOT_PREDICTIONS3 <- ggplot(spi, aes(Mono,Comp_spike , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~seawater~spike~µ~day^{-1}), y=expression(Fold~change~compared~to~no~spike))+
  geom_line(aes(SELF, SPIKED, colour = lineage), pred_all) 
SPIKE_PLOT_PREDICTIONS3 #oth95 fit not great, but we will do a mean fit anyways... 


spi_av<-ddply(spi,c("Foc_lin","biorep_focal"), function(df) return(c(muecompspi.avg=mean(df$Comp_to_spike),muecomp.se=sd(df$Comp_to_spike)/sqrt(length(coc)))))

## GFP ##
GFPco<-read.csv("20190418_foldchangedataCO2_GFP.csv")
head(GFPco)
#need a biorep id column for the random effect to get the nesting right
GFPco<-within(GFPco,lineage_id_full <- (factor(techrep):factor(biorep):factor(against)))
head(GFPco)
str(GFPco)
#this is the full model - now, because R is an oafish fool, we need to turn everything back into factors... 
GFPco$Foc_lin<-as.factor(GFPco$Foc_lin)
GFPco$against<-as.factor(GFPco$against)
#fullest global model (want to exclude data where they have been self-spiked?)

exp.mix_GFP <- nlme( GFP_co~ exp_dec(a, b, growth = Mono),
                     data = GFPco,
                     fixed = list(a + b ~ 1), 
                     random = list(lineage_id_full = pdDiag(a + b ~ 1)),
                     groups = ~ Foc_lin,
                     start = c(0.2, -0.003), 
                     na.action = na.omit,
                     method = 'ML',
                     control = nlmeControl(pnlsTol = 0.02))
summary(exp.mix_GFP)

# both parameters significant, so update full model as REML
exp.mix_final_GFP <- update(exp.mix_GFP, method = 'REML')
summary(exp.mix_final_GFP)
#write.csv(prettify(summary(exp.mix_GFP)),"model_all_out_GFP.csv")


### plot, generate predicted curves from fixed effects 

pred_alone <- seq(min(GFPco$Mono), max(GFPco$Mono), by = 0.01)

#predicted curves per biorep in co-culture - for this we need to go into the fixed effects as returned by the model 
#check which fixed effect is which 
names(fixef(exp.mix_final_GFP))

# easier here - don't have multiple self/non-self, just GFP or NOT GFP 
pred_curve <- exp_dec(a = fixef(exp.mix_final_GFP)[1] , b = fixef(exp.mix_final_GFP)[2], growth = pred_alone) 
#make a date frame containing all predictions
pred_all <- data.frame(ALONE = rep(pred_alone,8), lineage=c(rep('oth95', length(pred_alone)),rep('rcc1107', length(pred_alone)),rep('rcc1114', length(pred_alone)),rep('rcc1645', length(pred_alone)),rep('rcc410', length(pred_alone)),rep('rcc747', length(pred_alone)),rep('rcc789', length(pred_alone)),rep('rcc810', length(pred_alone))),GFP = rep(pred_curve,8))
#check
str(pred_all)


# this against the not means data frame , but we make a means frame below because we can... 

GFP_PLOT_PREDICTIONS3 <- ggplot(GFPco, aes(Mono,GFP_co , colour = Foc_lin)) +
  geom_point(size=4)+
  theme_classic(base_size = 24, base_family = 'Helvetica')+
  scale_colour_ordinal("Focal lineage")+
  theme(legend.position = c(0.7,.8),legend.box = "horizontal")+
  labs(x=expression(Growth~alone~µ~day^{-1}), y=expression(Fold~change~withGFP~compared~to~monoculture)) +
  geom_line(aes(ALONE, GFP, colour = lineage), pred_all,col="black") 
GFP_PLOT_PREDICTIONS3


GFP_av<-ddply(GFPco,c("Foc_lin","biorep"), function(df) return(c(muecompspi.avg=mean(df$GFP_co),muecomp.se=sd(df$GFP_co)/sqrt(length(GFPco)))))

