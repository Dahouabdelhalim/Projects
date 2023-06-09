################################################################################
### R Code for the analysis used in Archer et al. 2021
### Associations between metabolic traits and growth rate in brown trout (Salmo trutta)
### depend on thermal regime
### Last edited by Louise Archer, University College Cork, June 2021 
### Analysis carried out with R version 4.0.4


rm(list=ls())

# load packages
library(plyr)
library("broom")
library("purrr")
library("dplyr")
library("ggplot2")
library("tidyr")
library(forcats)
library("reshape2")
library("car")
library("MASS")
library(ggpubr)
library(nlme)
library(MuMIn)
library(Cairo)
library(dotwhisker)
library(AICcmodavg)
library(ggeffects) 
library(modelr)
library(effectsize)
library(smatr)

## set ggplot plotting theme 
my_theme = theme(
  text = element_text(size = 20, colour="black"),
  axis.title = element_text(colour="black"),
  axis.text = element_text(size = 16, colour="black"),
  legend.title=element_text(size=20),
  legend.text=element_text(size=20),
  legend.position = "right",
  legend.key = element_rect(fill=NA),
  legend.background = element_blank(),
  panel.grid.minor=element_blank(),
  panel.background=element_rect(colour="black", fill=NA, size=0.5),
  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
  plot.margin=unit(c(0.15,1,0.5,0.5),"cm"))


# read in the data with metabolic rate measurements
metdat <- read.csv("Archer_et_al_metabolism_data_PRSB.csv")


### first correct for body mass i.e. scaling of metabolism with body mass

# correct SMR for body mass
modRMR1 <- lm(log10(RMR) ~ log10(Massg), data=metdat)
Anova(modRMR1)
summary(modRMR1)
plot(log10(RMR) ~ log10(Massg), data=metdat)
metdat$rRMR <- resid(modRMR1) #save the residuals as mass-corrected SMR

# correct MMR for body mass
modMMR1 <- lm(log10(MMR) ~ log10(Massg), data=metdat)
Anova(modMMR1)
summary(modMMR1)
plot(log10(MMR) ~ log10(Massg), data=metdat)
metdat$rMMR <- resid(modMMR1) #save the residuals as mass-corrected MMR

# correct AS for body mass
modAS1 <- lm(log10(AS) ~ log10(Massg), data=metdat)
Anova(modAS1)
summary(modAS1)
plot(log10(AS) ~ log10(Massg), data=metdat)
metdat$rAS <- resid(modAS1) #save the residuals as mass-corrected AS

# check the model assumptions
par(mfrow=c(2,2))
plot(modRMR1)
plot(modMMR1)
plot(modAS1)


## save the bodymass scaling results as a tidy dataframe
modRMR_tidy <- tidy(modRMR1)
modRMR_tidy  <- modRMR_tidy  %>% 
  mutate(response="RMR")

modMMR_tidy <- tidy(modMMR1)
modMMR_tidy  <- modMMR_tidy  %>% 
  mutate(response="MMR")

modAS_tidy <- tidy(modAS1)
modAS_tidy  <- modAS_tidy  %>% 
  mutate(response="AS")

Body_mass_scaling_params <- rbind(modRMR_tidy, modMMR_tidy, modAS_tidy)


## create a figure to show body mass scaling
#create model predictions   
newdata1 <- data.frame(Massg = seq(min(metdat$Massg), max(metdat$Massg), length.out=10000))
predictedRMR <- data.frame(predict(modRMR1, newdata1, interval = "confidence"))
predictedMMR <- data.frame(predict(modMMR1, newdata1, interval = "confidence"))
predictedAS <- data.frame(predict(modAS1, newdata1, interval = "confidence"))
  
newdata1 <- newdata1 %>%
  mutate(pred_RMR = predictedRMR$fit, RMR_lwr = predictedRMR$lwr, RMR_upr = predictedRMR$upr,
         pred_MMR = predictedMMR$fit, MMR_lwr = predictedMMR$lwr, MMR_upr = predictedMMR$upr,
         pred_AS = predictedAS$fit, AS_lwr = predictedAS$lwr, AS_upr = predictedAS$upr)

RMR_BM_scaling <- ggplot()+geom_line(data=newdata1, aes(x=log10(Massg), y=pred_RMR ), size=0.85)+#geom_line(size=0.85)+
  geom_ribbon(data=newdata1,aes(ymin = RMR_lwr, ymax = RMR_upr, x= log10(Massg)), colour=NA, alpha=.1)+
  geom_point(data=metdat, aes(x=log10(Massg), y=log10(RMR)),size=3, inherit.aes=F)+my_theme+
  xlab(bquote('\\n'*Log[10]*'Mass (g)')) + ylab(bquote('\\n'*Log[10]*'SMR (mg'*O[2]*'h'*r^-1*')')) + 
  theme(axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2.5), axis.title=element_text(size=28), axis.text=element_text(size=18))

MMR_BM_scaling <- ggplot()+geom_line(data=newdata1, aes(x=log10(Massg), y=pred_MMR ), size=0.85)+#geom_line(size=0.85)+
  geom_ribbon(data=newdata1,aes(ymin = MMR_lwr, ymax = MMR_upr, x= log10(Massg)), colour=NA, alpha=.1)+
  geom_point(data=metdat, aes(x=log10(Massg), y=log10(MMR)),size=3, inherit.aes=F)+my_theme+
  xlab(bquote('\\n'*Log[10]*'Mass (g)')) + ylab(bquote('\\n'*Log[10]*'MMR (mg'*O[2]*'h'*r^-1*')')) + 
  theme(axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2.5), axis.title=element_text(size=28), axis.text=element_text(size=18))

AS_BM_scaling <- ggplot()+geom_line(data=newdata1, aes(x=log10(Massg), y=pred_AS ), size=0.85)+#geom_line(size=0.85)+
  geom_ribbon(data=newdata1,aes(ymin = AS_lwr, ymax = AS_upr, x= log10(Massg)), colour=NA, alpha=.1)+
  geom_point(data=metdat, aes(x=log10(Massg), y=log10(AS)),size=3, inherit.aes=F)+my_theme+
  xlab(bquote('\\n'*Log[10]*'Mass (g)')) + ylab(bquote('\\n'*Log[10]*'AS (mg'*O[2]*'h'*r^-1*')')) + 
  theme(axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2.5), axis.title=element_text(size=28), axis.text=element_text(size=18))


# ggarrange(RMR_BM_scaling, MMR_BM_scaling,AS_BM_scaling, labels=c("(A)", "(B)", "(C)"),font.label = list(size = 28),
# label.x=c(0.16),label.y=c(0.975),common.legend = F, ncol=3)



#############################################
## Temperature effects on metabolic rates ###
#############################################

# Test for differences in metabolic rates according to temperature treatment

# re order factors levels
metdat$Population <- factor(metdat$Population, levels=c("Erriff_N","Rough"))
metdat$Temperature <- factor(metdat$Temperature, levels=c("Low", "High"))

##################
##### SMR ########
################## 

## build the model testing for temperature effects on SMR, while
## adjusting for body mass and population background
modRMR1 <- lm(log10(RMR) ~ log10(Massg)*Temperature + Population, data=metdat)
Anova(modRMR1)
summary(modRMR1)

# drop the interaction with body mass
modRMR <- lm(log10(RMR) ~ log10(Massg)+ Temperature + Population, data=metdat)
Anova(modRMR)
summary(modRMR)
anova(modRMR1,modRMR)

## best fitting model 
modRMR <- lm(log10(RMR) ~ log10(Massg) + Population + Temperature, data=metdat)
Anova(modRMR)
summary(modRMR)


## check model assumptions
par(mfrow=c(2,2))
plot(modRMR)





##################################
### Effect size bootstrapping#####
##################################

## boostrapped resampling to create bootstrap confidence intervals
fit_boots <- metdat %>% subset(., select = c(RMR, Population, Temperature, Massg)) %>% 
  modelr::bootstrap(n = 10000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(fit = map(strap, ~cohens_f(Anova(lm(log10(RMR) ~ Population+Temperature+log10(Massg), data = .x )), ci=0.95)
  ))

params_boot <- fit_boots %>%
  unnest(fit)

params_boot_RMR <- params_boot %>%  
  mutate(Parameter = fct_relevel(Parameter, c("log10(Massg)", "Temperature", "Population")),
         Rate = "SMR")

## check the distributions of bootstrapped effect sizes
ggplot(params_boot_RMR, aes(Cohens_f_partial)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ Parameter, scales = 'free_x')+my_theme

# create 95% confidence intervals
confint_RMR <- group_by(params_boot_RMR, Parameter) %>%
  summarise(.,
            conf_low = quantile( Cohens_f_partial, 0.025),
            conf_high = quantile( Cohens_f_partial, 0.975),
            estimate = quantile( Cohens_f_partial, 0.5)) %>%
  ungroup() %>%
  mutate(., method = 'boot', rate= 'SMR')


# Aesthetics for creating figures
source(file="flat_violin_plots.R")

# create the figure
RMR_eff_size <- ggplot(confint_RMR, aes(Parameter, estimate)) +  ylim(0, 2.6)+
  geom_flat_violin(data=params_boot_RMR, mapping=aes(Parameter, Cohens_f_partial), width=0.75, adjust=3, size=0, colour=NA,alpha=0.5)+
  geom_point(data=confint_RMR, aes(Parameter, estimate),size = 3) +
  geom_linerange(data=confint_RMR,aes(ymin = conf_low, ymax = conf_high)) + my_theme + 
  geom_hline(yintercept = 0, lty=5, alpha=.5, size=0.5) +
  # ylab(bquote("Effect size (partial "*Omega^2*")"))+
  ylab(bquote("Effect size (Cohen's "~italic(f)~")"))+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=28), axis.text=element_text(size=18))+
  scale_x_discrete(breaks=c("log10(Massg)","Temperature","Population"),
                   labels=c(bquote('\\n'*Log[10]*'Mass (g)'), "Temperature", "Population"))
 

##################
##### MMR ########
##################

# build the model testing for temperature effects on MMR, while 
# adjusting for body mass and population background
modMMR1 <- lm(log10(MMR) ~  log10(Massg)*Temperature + Population, data = metdat)
Anova(modMMR1)
summary(modMMR)

## drop the temperature x body mass interaction
modMMR <- lm(log10(MMR) ~  log10(Massg) + Temperature + Population, data = metdat)
Anova(modMMR)
summary(modMMR)
anova(modMMR1, modMMR)

## best fitting model 
modMMR <- lm(log10(MMR) ~  log10(Massg) + Temperature + Population, data = metdat)
Anova(modMMR)
summary(modMMR)

## check model assumptions
par(mfrow=c(2,2))
plot(modMMR)


##################################
### Effect size bootstrapping#####
##################################

## boostrapped resampling to create bootstrap confidence intervals
fit_boots <- metdat %>% subset(., select = c(MMR, Population, Temperature, Massg)) %>% 
  modelr::bootstrap(n = 10000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(fit = map(strap, ~cohens_f(Anova(lm(log10(MMR) ~ Population+Temperature+log10(Massg), data = .x )), ci=0.95)
  ))

params_boot <- fit_boots %>%
  unnest(fit)

params_boot_MMR <- params_boot %>%  
  mutate(Parameter = fct_relevel(Parameter, c("log10(Massg)", "Temperature", "Population")),
         Rate = "MMR")

## check the distribution of bootstrapped effect sizes
ggplot(params_boot_MMR, aes(Cohens_f_partial)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ Parameter, scales = 'free_x')+my_theme

## construct the confidence intervals 
confint_MMR <- group_by(params_boot_MMR, Parameter) %>%
  summarise(.,
            conf_low = quantile( Cohens_f_partial, 0.025),
            conf_high = quantile( Cohens_f_partial, 0.975),
            estimate = quantile( Cohens_f_partial, 0.5)) %>%
  ungroup() %>%
  mutate(., method = 'boot', rate='MMR')


# create the plot
MMR_eff_size <- ggplot(confint_MMR, aes(Parameter, estimate)) +   ylim(0.0, 1.89)+
  geom_flat_violin(data=params_boot_MMR, mapping=aes(Parameter, Cohens_f_partial), width=0.75, adjust=3, size=0, colour=NA,alpha=0.5)+
  geom_point(data=confint_MMR, aes(Parameter, estimate),size = 3) +
  geom_linerange(data=confint_MMR,aes(ymin = conf_low, ymax = conf_high)) + my_theme + 
  geom_hline(yintercept = 0, lty=5, alpha=.5, size=.5) +
  # ylab(bquote("Effect size (partial "*Omega^2*")"))+
  ylab(bquote("Effect size (Cohen's "~italic(f)~")"))+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=28), axis.text=element_text(size=18))+
  scale_x_discrete(breaks=c("log10(Massg)","Temperature","Population"),
                   labels=c(bquote('\\n'*Log[10]*'Mass (g)'), "Temperature", "Population"))


##################
##### AS #########
##################

# build the model testing for temperature effects on AS, while 
# adjusting for body mass and population background
modAS1 <- lm(log10(AS) ~ log10(Massg)*Temperature + Population, data = metdat)
Anova(modAS1)
summary(modAS1)

## drop the body mass * temperature interaction
modAS <- lm(log10(AS) ~ log10(Massg) + Temperature + Population, data = metdat)
Anova(modAS)
summary(modAS)

# best fitting model
modAS <- lm(log10(AS) ~ log10(Massg) + Temperature + Population, data = metdat)
Anova(modAS)
summary(modAS)


## check model assumptions
par(mfrow=c(2,2))
plot(modAS)

##################################
### Effect size bootstrapping#####
##################################

## boostrapped resampling to create bootstrap confidence intervals 
fit_boots <- metdat %>% subset(., select = c(AS, Population, Temperature, Massg)) %>% 
  modelr::bootstrap(n = 10000, id = 'boot_num') %>%
  group_by(boot_num) %>%
  mutate(fit = map(strap, ~cohens_f(Anova(lm(log10(AS) ~ Population+Temperature+log10(Massg), data = .x )))
  ))

params_boot <- fit_boots %>%
  unnest(fit)

params_boot_AS <- params_boot %>%  
  mutate(Parameter = fct_relevel(Parameter, c("log10(Massg)", "Temperature", "Population")),
         Rate = "AS")

## check the distribution of bootstrapped effect sizes
ggplot(params_boot_AS, aes(Cohens_f_partial)) +
  geom_histogram(col = 'black', fill = 'white') +
  facet_wrap(~ Parameter, scales = 'free_x')+my_theme

## construct the confidence intervals
confint_AS <- group_by(params_boot_AS, Parameter) %>%
  summarise(.,
            conf_low = quantile( Cohens_f_partial, 0.025),
            conf_high = quantile( Cohens_f_partial, 0.975),
            estimate = quantile( Cohens_f_partial, 0.5)) %>%
  ungroup() %>%
  mutate(., method = 'boot', rate='AS')

# create the figure
AS_eff_size <- ggplot(confint_AS, aes(Parameter, estimate)) +  ylim(-0.05, 1.89)+
  geom_flat_violin(data=params_boot_AS, mapping=aes(Parameter, Cohens_f_partial), width=0.75, adjust=3, size=0, colour=NA,alpha=0.5)+
  geom_point(data=confint_AS, aes(Parameter, estimate),size = 3) +
  geom_linerange(data=confint_AS,aes(ymin = conf_low, ymax = conf_high)) + my_theme + 
  # geom_hline(yintercept = 0, lty=2, alpha=.5, size=1) +
  # ylab(bquote("Effect size (partial "*Omega^2*")"))+
  ylab(bquote("Effect size (Cohen's "~italic(f)~")"))+
  theme(axis.title.x=element_blank(), axis.title=element_text(size=28), axis.text=element_text(size=18))+
  scale_x_discrete(breaks=c("log10(Massg)","Temperature","Population"),
                   labels=c(bquote('\\n'*Log[10]*'Mass (g)'), "Temperature", "Population"))


# multi-panel figures
## create predictions for SMR according to body mass and temperature
df_rmr <- ggeffect(modRMR, terms=c("Massg [n=100]", "Temperature"))

df_rmr <- df_rmr %>%
  mutate(Temperature = group)

## create the figure for SMR
SMR_BM <- ggplot()+geom_line(data=df_rmr, aes(x=log10(x), y=predicted, colour=Temperature, linetype=Temperature),size=0.85)+
  geom_ribbon(data=df_rmr, aes(x=log10(x), ymin = conf.low, ymax = conf.high,fill=Temperature), colour=NA, alpha = .1) +
  geom_point(data=metdat, aes(x=log10(Massg), y=log10(RMR),shape=Temperature, colour=Temperature),size=3, inherit.aes=F)+my_theme+
  xlab(bquote('\\n'*Log[10]*'Mass (g)')) + ylab(bquote('\\n'*Log[10]*'SMR (mg'*O[2]*'h'*r^-1*')')) + 
  theme(legend.position = c(0.86, 0.09), legend.title=element_blank(),legend.key.width = unit(1.25,"cm"), 
        legend.text=element_text(size=16),axis.title.x=element_text(vjust=-0.25), axis.title.y=element_text(vjust=2.5), axis.title=element_text(size=28), axis.text=element_text(size=18))+
  scale_color_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool", "Warm"), values=c("#377eb8","#e41a1c"))+
  scale_fill_manual(breaks=c("Low" ,"High"),
                    labels=c("Cool", "Warm"), values=c("#377eb8","#e41a1c"))+
  scale_linetype_manual(breaks=c("Low" ,"High"),
                        labels=c("Cool", "Warm"),
                        values = c("solid","dashed")) +
  scale_shape_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool", "Warm"),
                     values = c(1,19))


####### for MMR figure #############

## create predictions for MMR according to body mass and temperature
df_mmr <- ggeffect(modMMR, terms=c("Massg [n=100]", "Population"))

df_mmr <- df_mmr %>%
  mutate(Population = group)

## create the figure for MMR
MMR_BM <- ggplot()+geom_line(data=df_mmr, aes(x=log10(x), y=predicted, colour=Population, linetype=Population),size=0.85)+
  geom_ribbon(data=df_mmr, aes(x=log10(x), ymin = conf.low, ymax = conf.high,fill=Population), colour=NA, alpha = .1) +
  geom_point(data=metdat, aes(x=log10(Massg), y=log10(MMR),shape=Population, colour=Population),size=3, inherit.aes=F)+my_theme+
  xlab(bquote('\\n'*Log[10]*'Mass (g)')) + ylab(bquote('\\n'*Log[10]*'MMR (mg'*O[2]*'h'*r^-1*')')) + 
  theme(legend.position = c(0.76, 0.09), legend.title=element_blank(),legend.key.width = unit(1.25,"cm"), 
        legend.text=element_text(size=16), axis.title.x=element_text(vjust=-0.25),axis.title.y=element_text(vjust=2.5), axis.title=element_text(size=28), axis.text=element_text(size=18))+
  scale_color_manual(breaks=c("Erriff_N" ,"Rough"),
                     labels=c("Anadromous", "Non-Anadromous"), values=c("black","orange"))+
  scale_fill_manual(breaks=c("Erriff_N" ,"Rough"),
                    labels=c("Anadromous", "Non-Anadromous"), values=c("black","orange"))+
  scale_linetype_manual(breaks=c("Erriff_N" ,"Rough"),
                        labels=c("Anadromous", "Non-Anadromous"),
                        values = c("dashed","solid")) +
  scale_shape_manual(breaks=c("Erriff_N" ,"Rough"),
                     labels=c("Anadromous", "Non-Anadromous"),
                     values = c(1,19))



# #### create a multipanel figure
# ggarrange(SMR_BM, RMR_eff_size ,MMR_BM, MMR_eff_size ,
#           labels=c("(A)", "(B)", "(C)", "(D)"),font.label = list(size = 28, face="bold"),
#           label.x=c(0.16, .16, .16, .16),label.y=c(0.975),common.legend = F,align="hv", ncol=2, nrow=2)
# dev.off()


PopTemp_Effect_sizes  <- rbind(confint_RMR, confint_MMR)


#####################################################################
###### Effects of metabolic rate & temperature on growth rates ######
#####################################################################

## read in the growth rate data
MetLong <- read.csv("Archer_et_al_growth_data_PRSB.csv")

# reorder factor levels
MetLong$Temperature <- factor(MetLong$Temperature, levels = c("Low", "High"))
MetLong$Population <- factor(MetLong$Population, levels = c("Erriff_N", "Rough"))


## Explore growth rates through time
ggplot(MetLong, aes(M_Num, y=Gs, colour=Temperature))+ geom_point()+my_theme+
  ylab(bquote("Specific growth rate (% day"^-1~")")) + xlab("Months of treatment")+geom_text(aes(label=DNA_ID))
## suggests polynomial growth rate ~ time
## potentially one outlier point - check later!!  



##############################################################
## Build the model testing effects of SMR and MMR on growth ##
##############################################################

## set the control for lme
ctrl <- lmeControl(opt='optim')

### build a model that includes both SMR and MMR, and their interactions with temperature. *also include temperature x time interaction
memRMR_MMR_time1 <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + rMMR*rRMR + Population + Temperature*poly(M_Num, 3) + Length,
                   random=~1|DNA_ID, 
                   control=ctrl,
                   data=MetLong,
                   correlation=corAR1(value=0.2, form=~M_Num), 
                   na.action=na.omit,
                   method="ML")



# drop the SMR x MMR interaction
memRMR_MMR_time2 <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3) + Length,
                        random=~1|DNA_ID, 
                        control=ctrl,
                        data=MetLong,
                        correlation=corAR1(value=0.2, form=~M_Num), 
                        na.action=na.omit,
                        method="ML")

# test for significance of MMR x SMR interaction
anova(memRMR_MMR_time1, memRMR_MMR_time2) 


# drop the SMR x temperature interaction
memRMR_MMR_time3 <- lme(Gs ~   rMMR*Temperature + rRMR*rMMR + Population + Temperature*poly(M_Num, 3)+Length,
                        random=~1|DNA_ID, 
                        control=ctrl,
                        data=MetLong,
                        correlation=corAR1(value=0.2, form=~M_Num), 
                        na.action=na.omit,
                        method="ML")

# test for significance of temperature x SMR interaction
anova(memRMR_MMR_time1, memRMR_MMR_time3) 


# drop the MMR x temperature interaction
memRMR_MMR_time4 <- lme(Gs ~   rRMR*Temperature +rRMR*rMMR +Population + Temperature*poly(M_Num, 3)+Length,
                        random=~1|DNA_ID, 
                        control=ctrl,
                        data=MetLong,
                        correlation=corAR1(value=0.2, form=~M_Num), 
                        na.action=na.omit,
                        method="ML")

# test for significance of temperature x MMR interaction
anova(memRMR_MMR_time1, memRMR_MMR_time4) 

# drop the time x temperature interaction
memRMR_MMR_time5 <- lme(Gs ~  rMMR*Temperature + rRMR*Temperature + rRMR*rMMR +Population + poly(M_Num, 3)+Length,
                        random=~1|DNA_ID, 
                        control=ctrl,
                        data=MetLong,
                        correlation=corAR1(value=0.2, form=~M_Num), 
                        na.action=na.omit,
                        method="ML")

# test for significance of temperature x time interaction
anova(memRMR_MMR_time1, memRMR_MMR_time5) 



## best fitting model = memRMR_MMR_time2
# check / confirm by comparing AICc
aictab(list(memRMR_MMR_time1, memRMR_MMR_time2, memRMR_MMR_time3, memRMR_MMR_time4, memRMR_MMR_time5), 
       modnames=c("full_model", "drop_MMR*SMR", "drop_SMR*Temp",  "drop_MMR*Temp" , "drop_time*temp"))


## Final reduced model excludes rSMR x rMMR interaction

# refit model using REML
memRMR_MMR <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3) + Length,
                   random=~1|DNA_ID, 
                   control=ctrl,
                   data=MetLong,
                   correlation=corAR1(value=0.2, form=~M_Num), 
                   na.action=na.omit,
                   method="REML")

Anova(memRMR_MMR)
summary(memRMR_MMR)

# calculate marginal r-sq
r.squaredGLMM(memRMR_MMR)
params_memRMR <- tidy(memRMR_MMR,effects="fixed") ## tidy up the fixed effects 


# check assumptions
plot(resid(memRMR_MMR, type="normalized")~fitted(memRMR_MMR))
qqnorm(resid(memRMR_MMR))
plot(acf(resid(memRMR_MMR, type = "normalized"))) 
plot(ACF(memRMR_MMR), alpha=0.05)


##############################################################
####### Build the model testing effects of AS on growth ######
##############################################################

## check correlations between metabolic traits and AS
with(metdat, cor.test(rAS,rRMR,  method="pearson"))
with(metdat, cor.test(rMMR, rAS, method="pearson"))

## model AS separately to SMR and MMR

## model selection on interaction terms for AS
memAS1 <- lme(Gs ~  rAS*Temperature + Population + Temperature*poly(M_Num, 3) + Length,
               random=~1|DNA_ID, 
               control=ctrl,
               data=MetLong,
               correlation=corAR1(value=0.2, form=~M_Num), 
               na.action=na.omit,
               method="ML")



# check rAS*Temperature interaction
memAS2 <- lme(Gs ~  rAS+Population+Temperature*poly(M_Num, 3)+Length,
               random=~1|DNA_ID, 
               control=ctrl,
               data=MetLong,
               correlation=corAR1(value=0.2, form=~M_Num), 
               na.action=na.omit,
               method="ML")

anova(memAS1, memAS2) # LRT test for rAS*Temperature term

# check time*Temperature interaction
memAS3 <- lme(Gs ~  Temperature*rAS+Population+poly(M_Num, 3)+Length,
              random=~1|DNA_ID, 
              control=ctrl,
              data=MetLong,
              correlation=corAR1(value=0.2, form=~M_Num), 
              na.action=na.omit,
              method="ML")

anova(memAS1, memAS3) # LRT test for time*Temperature term


## best fitting model = memAS2
# check / confirm by comparing AICc - equal support for memAS1 and memAS2 
aictab(list(memAS1, memAS2, memAS3), modnames=c("full_model", "drop_AS*Temp", "drop_time*Temp"))

## Most parsimonious model still = memAS2

#### refit using REML
memAS <- lme(Gs ~  rAS+Population+Temperature*poly(M_Num, 3)+Length,
              random=~1|DNA_ID, 
              control=ctrl,
              data=MetLong,
              correlation=corAR1(value=0.2, form=~M_Num), 
              na.action=na.omit,
              method="REML")

Anova(memAS)

# calculate marginal r-sq
r.squaredGLMM(memAS)

# check assumptions
plot(resid(memAS, type="normalized")~fitted(memAS))
qqnorm(resid(memAS))
plot(acf(resid(memAS, type = "normalized"))) 
plot(ACF(memAS), alpha=0.05)



###### STOPPED HERE FOR BREAKFAST 30.06

###########################################
##### Make the figures ###################
##########################################

## make the figures for predicting marginal effects
## uses package ggeffects

# load the functions for calculating and plotting partial residuals within ggeffects
## function for partial residuals available on ggeffects github https://github.com/strengejacke/ggeffects
source(file="partial.resid.plots.R")


## create the plots for rSMR and rMMR model results 

## predicting marginal effects
df_rmr1<- ggeffect(memRMR_MMR, terms=c("rRMR [all]", "Temperature"))
df_mmr1<- ggeffect(memRMR_MMR, terms=c("rMMR [all]", "Temperature"))

# calculate partial residuals for SMR
rmr_res <- residualize_over_grid.ggeffects(df_rmr1,memRMR_MMR)
# calculate partial residuals for MMR
mmr_res <- residualize_over_grid.ggeffects(df_mmr1,memRMR_MMR)


## create the figure for SMR
rmr_mar1 <- ggplot(df_rmr1, aes(x=x, y=predicted, colour=group))+
  geom_line(aes(linetype=group, colour=group),size=0.85)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), colour=NA, alpha = .1)+my_theme+
  # geom_rug(data=MetLong, aes(x=rRMR), sides="b", inherit.aes=F, colour="grey48")+ ## to add a rug for rSMR data
  geom_point(data = rmr_res, aes(x=x, y=predicted), alpha=0.2)+
  xlab('\\nrSMR (mg'*~O[2]~hr^-1*')')+ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+
  scale_x_continuous(breaks=c(-.10, 0.0, 0.10))+
  # scale_y_continuous(breaks=c(-.10, 0.0, 0.10))+
  # ylim(c(0.04, .224))+
  labs(fill='group', colour="group") +
  scale_color_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool" ,"Warm"), values=c("#377eb8","#e41a1c"))+
  scale_fill_manual(breaks=c("Low" ,"High"),
                    labels=c("Cool" ,"Warm"), values=c("#377eb8", "#e41a1c"))+
  scale_linetype_manual(breaks=c("Low" ,"High"),labels=c("Cool" ,"Warm"),
                        values = c("solid","dashed")) +
  theme(legend.position=c(.35,.94), legend.title=element_blank(), legend.text=element_text(size=14),plot.margin=unit(c(0.5,0.5,0.25,1),"cm"),
        axis.title.y=element_text(vjust=2.5, hjust=.80),axis.title.x=element_text(vjust=-0.1), strip.background=element_rect(fill=NA))+guides(col = guide_legend(ncol = 2))

## create the figure for MMR
mmr_mar1 <- ggplot(df_mmr1, aes(x=x, y=predicted, colour=group))+
  geom_line(aes(linetype=group, colour=group),size=0.85)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), colour=NA, alpha = .1)+my_theme+
  # geom_rug(data=MetLong, aes(x=rRMR), sides="b", inherit.aes=F, colour="grey48")+ ## to add a rug for rSMR data
  geom_point(data = mmr_res, aes(x=x, y=predicted), alpha=0.2)+
  xlab('\\nrMMR (mg'*~O[2]~hr^-1*')')+
  ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+
  scale_x_continuous(breaks=c(-.10, 0.0, 0.10))+
  # scale_y_continuous(breaks=c(-.10, 0.0, 0.10))+
  labs(fill='group', colour="group") +
  scale_color_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool" ,"Warm"), values=c("#377eb8","#e41a1c"))+
  scale_fill_manual(breaks=c("Low" ,"High"),
                    labels=c("Cool" ,"Warm"), values=c("#377eb8", "#e41a1c"))+
  scale_linetype_manual(breaks=c("Low" ,"High"),labels=c("Cool" ,"Warm"),
                        values = c("solid","dashed")) +
  theme(legend.position=c(.35,.94), legend.title=element_blank(), legend.text=element_text(size=14),plot.margin=unit(c(0.5,0.5,0.25,1),"cm"),
        axis.title.y=element_text(vjust=2.5,hjust=.80),axis.title.x=element_text(vjust=-0.1), strip.background=element_rect(fill=NA))+guides(col = guide_legend(ncol = 2))




# create multiplot of SMR and MMR marginal effects
smr_mmr_multi <- ggarrange(rmr_mar1, mmr_mar1 ,labels=c("(B)", "(C)"),font.label=list(size=18, face="bold"),
                           label.x=c(0.21,0.21),label.y=.94,common.legend = F, nrow=2)

smr_mmr_multi2 <- ggarrange(rmr_mar1, mmr_mar1 ,labels=c("(B)", "(C)"),font.label=list(size=18, face="bold"),
                            label.x=c(0.25,0.25),label.y=.96,common.legend = F, ncol=2)


## tidy up the fixed effects
GrowthRMR_MMR_tidy <- tidy(memRMR_MMR, effects="fixed")
GrowthRMR_MMR_tidy <- GrowthRMR_MMR_tidy %>% 
  mutate(model="SMR, MMR")

GrowthAS_tidy <- tidy(memAS, effects="fixed")
GrowthAS_tidy <- GrowthAS_tidy %>% 
  mutate(model="AS")

Growth_models <- rbind(GrowthRMR_MMR_tidy,  GrowthAS_tidy)

Growth_models <- Growth_models %>% 
  relabel_predictors(c("(Intercept)"= "Intercept",
                       "rRMR" = "rSMR",
                       "TemperatureHigh" = "Temperature: Warm",
                       "rMMR" = "rMMR",
                       "PopulationRough" = "Non-Anadromous",
                       "poly(M_Num, 3)1" = "Time", 
                       "poly(M_Num, 3)2" = "Time^2", 
                       "poly(M_Num, 3)3" = "Time^3", 
                       "Length" = "Initial length",
                       "rRMR:TemperatureHigh" = "rSMR * Warm",
                       "TemperatureHigh:rMMR" = "rMMR * Warm",
                       "TemperatureHigh:poly(M_Num, 3)1" = "Time * Warm", 
                       "TemperatureHigh:poly(M_Num, 3)2" = "Time^2 * Warm", 
                       "TemperatureHigh:poly(M_Num, 3)3" = "Time^3 * Warm", 
                       "rAS" = "rAS"
                       
  ))


# Save the model summaries 


ordered_vars <- c("Intercept","Time","Time^2","Time^3", "Initial length","rSMR",
                  "rMMR", "rAS",
                  "Temperature: Warm", "Non-Anadromous",
                  "rSMR * Warm", "rMMR * Warm",  "Time * Warm", "Time^2 * Warm", "Time^3 * Warm")

## order all the  model dataframes 
Growth_models <- Growth_models %>%
  mutate(term =  factor(term, levels = ordered_vars)) %>% 
  arrange(term)

# make a plot summarising model co-efficients
Growth_dw <- dwplot(Growth_models,show_intercept =T, dot_args = list(aes(shape=model),size = 2), 
                    whisker_args = list(size = .8), dodge_size = 0.7)+
  theme_bw()+geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size=0.8) + xlab("Coefficient estimate") +
  theme(text = element_text(size = 14), #should be 14 for combined figure
        legend.position = "top",# for combined figure #1 + # 4
        # legend.position = "right",# for combined figure #2 + #3
        legend.justification=c(1, 1),
        legend.title = element_text(size=14),
        legend.text = element_text(size=14),
        legend.background = element_rect(colour="grey80"),
        axis.title.x = element_text(size=20),
        axis.text.x =element_text(size=16),
        axis.text.y = element_text(size=16, colour="black")) +
  scale_colour_manual(name="Model:", breaks=c("SMR, MMR", "AS"),
                      values=c("black", "orange")) +
  scale_shape_manual(name="Model:", breaks=c("SMR, MMR", "AS"),
                     values=c(19, 21))+
  scale_y_discrete(breaks=c("Intercept","Time","Time^2","Time^3", "Initial length","rSMR",
                            "rMMR", "rAS",
                            "Temperature: Warm", "Non-Anadromous",
                            "rSMR * Warm", "rMMR * Warm", 
                            "Time * Warm", "Time^2 * Warm", "Time^3 * Warm"),
                   labels=c("Intercept","Time",bquote(Time^2),bquote(Time^3),"Initial length",
                            "rSMR","rMMR", "rAS",
                            "Temperature: Warm", "Non-Anadromous",
                            "rSMR * Warm", "rMMR * Warm",
                            "Time * Warm", bquote(Time^2*'* Warm'), bquote(Time^3*'* Warm')))+
  theme(plot.margin=unit(c(2.75,0.1,2.75,0),"cm"),# for combined figure #1
        # theme(plot.margin=unit(c(0.5,7,0.5,0.25),"cm"),# for combined figure #2
        # theme(plot.margin=unit(c(0.5,4,1,4),"cm"),# for combined figure #3
        # theme(plot.margin=unit(c(0.5,0.5,0.5,0),"cm"), # for combined figure #4
        axis.title.x=element_text(vjust=-1.75)) #for combined figure

# combined fig v.1 - final for manuscript. 
# ggarrange(Growth_dw,smr_mmr_multi ,labels=c("(A)"),font.label=list(size=18, face="bold"),
#           label.x=c(0.32, 0.23),label.y=c(.88),common.legend = F, ncol=2, widths=c(1, .8), heights=c(0.5, 1))
# dev.off()


######################################
#### Sensitivity to outlier point ####
######################################

# check sensitivity of results to the outlier point (which has unusually high growth rate)
# create new datadframe and drop the outlier point
Sense <- MetLong %>% filter(., !(Gs > 0.301  & DNA_ID == "R65"))

# fit the model
memRMR_MMR_sense <- lme(Gs ~  rRMR * Temperature + rMMR * Temperature + Population + Temperature*poly(M_Num, 3)+Length,
                        random=~1|DNA_ID, 
                        control=ctrl,
                        data=Sense,
                        correlation=corAR1(value=0.2, form=~M_Num), 
                        na.action=na.omit,
                        method="ML")

# check how results change
Anova(memRMR_MMR_sense) # same qualitative trends

# make some figures to further check

####################################
## Figures for sensitivity analysis
####################################

# predict margninal effects for the sensitivity analysis
df_rmr_sense <- ggeffect(memRMR_MMR_sense, terms=c("rRMR [all]", "Temperature"))
df_mmr_sense <- ggeffect(memRMR_MMR_sense, terms=c("rMMR [all]", "Temperature"))

# calculate partial residuals for sensitivity analysis
rmr_res_sense <- residualize_over_grid.ggeffects(df_rmr_sense,memRMR_MMR_sense)
mmr_res_sense <- residualize_over_grid.ggeffects(df_mmr_sense,memRMR_MMR_sense)


## create the figure for SMR (outlier excluded)
rmr_mar_sense <- ggplot(df_rmr_sense, aes(x=x, y=predicted, colour=group))+
  geom_line(aes(linetype=group, colour=group),size=0.85)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), colour=NA, alpha = .1)+my_theme+
  # geom_rug(data=MetLong, aes(x=rRMR), sides="b", inherit.aes=F, colour="grey48")+ ## to add a rug for rSMR data
  geom_point(data = rmr_res_sense, aes(x=x, y=predicted), alpha=0.2)+
  xlab('\\nrSMR (mg'*~O[2]~hr^-1*')')+ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+
  scale_x_continuous(breaks=c(-.10, 0.0, 0.10))+
  # scale_y_continuous(breaks=c(-.10, 0.0, 0.10))+
  # ylim(c(0.04, .224))+
  ylim(c(0.02, .3))+ # axis lims for partial residuals
  labs(fill='group', colour="group") +
  scale_color_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool" ,"Warm"), values=c("#377eb8","#e41a1c"))+
  scale_fill_manual(breaks=c("Low" ,"High"),
                    labels=c("Cool" ,"Warm"), values=c("#377eb8", "#e41a1c"))+
  scale_linetype_manual(breaks=c("Low" ,"High"),labels=c("Cool" ,"Warm"),
                        values = c("solid","dashed")) +
  theme(legend.position=c(.35,.94), legend.title=element_blank(), legend.text=element_text(size=14),plot.margin=unit(c(0.5,0.5,0.25,1),"cm"),
        axis.title.y=element_text(vjust=2.5, hjust=.80),axis.title.x=element_text(vjust=-0.1), strip.background=element_rect(fill=NA))+guides(col = guide_legend(ncol = 2))

## create the figure for MMR (outlier excluded)
mmr_mar_sense <- ggplot(df_mmr_sense, aes(x=x, y=predicted, colour=group))+
  geom_line(aes(linetype=group, colour=group),size=0.85)+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high,fill=group), colour=NA, alpha = .1)+my_theme+
  # geom_rug(data=MetLong, aes(x=rRMR), sides="b", inherit.aes=F, colour="grey48")+ ## to add a rug for rSMR data
  geom_point(data = mmr_res_sense, aes(x=x, y=predicted), alpha=0.2)+
  xlab('\\nrMMR (mg'*~O[2]~hr^-1*')')+
  ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+
  scale_x_continuous(breaks=c(-.10, 0.0, 0.10))+
  # scale_y_continuous(breaks=c(-.10, 0.0, 0.10))+
  # ylim(c(0.04, .224))+
  ylim(c(0.02, .3))+ # axis lims for partial residuals
  labs(fill='group', colour="group") +
  scale_color_manual(breaks=c("Low" ,"High"),
                     labels=c("Cool" ,"Warm"), values=c("#377eb8","#e41a1c"))+
  scale_fill_manual(breaks=c("Low" ,"High"),
                    labels=c("Cool" ,"Warm"), values=c("#377eb8", "#e41a1c"))+
  scale_linetype_manual(breaks=c("Low" ,"High"),labels=c("Cool" ,"Warm"),
                        values = c("solid","dashed")) +
  theme(legend.position=c(.35,.94), legend.title=element_blank(), legend.text=element_text(size=14),plot.margin=unit(c(0.5,0.5,0.25,1),"cm"),
        axis.title.y=element_text(vjust=2.5,hjust=.80),axis.title.x=element_text(vjust=-0.1), strip.background=element_rect(fill=NA))+guides(col = guide_legend(ncol = 2))



# ggarrange(rmr_mar_sense, mmr_mar_sense ,labels=c("(A)", "(B)"),font.label=list(size=18, face="bold"),
#                             label.x=c(0.8,0.8),label.y=.96,common.legend = F, ncol=2)





###################
### Supp Info #####
###################

### SUPPLEMENTARY ANALYSIS
## check for consistency of effects through time (rSMR * time, rMMR * time)

#Check for consistency of SMR through time
memSMR_time <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3)+rRMR*poly(M_Num, 3) + Length,
                  random=~1|DNA_ID, 
                  control=ctrl,
                  data=MetLong,
                  correlation=corAR1(value=0.2, form=~M_Num), 
                  na.action=na.omit,
                  method="ML")

anova(memSMR_time, memRMR_MMR_time1) # compare to original model without rSMR * time
Anova(memSMR_time)
# best model retains SMR x time

#Check for consistency of MMR through time
memMMR_time <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3)+rMMR*poly(M_Num, 3) + Length,
                       random=~1|DNA_ID, 
                       control=ctrl,
                       data=MetLong,
                       correlation=corAR1(value=0.2, form=~M_Num), 
                       na.action=na.omit,
                       method="ML")

anova(memMMR_time, memRMR_MMR_time1) # compare to original model without rMMR * time
Anova(memMMR_time)
# best model retains MMR x time


# refit both models using REML - SMR x time 
memRMR_MMR_SI1 <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3)+rRMR*poly(M_Num, 3) + Length,
                       random=~1|DNA_ID, 
                       control=ctrl,
                       data=MetLong,
                       correlation=corAR1(value=0.2, form=~M_Num), 
                       na.action=na.omit,
                       method="REML")

# refit both models using REML - MMR  x time 
memRMR_MMR_SI2 <- lme(Gs ~  rRMR*Temperature + rMMR*Temperature + Population + Temperature*poly(M_Num, 3)+rMMR*poly(M_Num, 3) + Length,
                      random=~1|DNA_ID, 
                      control=ctrl,
                      data=MetLong,
                      correlation=corAR1(value=0.2, form=~M_Num), 
                      na.action=na.omit,
                      method="REML")


# check assumptions - SMR x time 
plot(resid(memRMR_MMR_SI1, type="normalized")~fitted(memRMR_MMR_SI1))
qqnorm(resid(memRMR_MMR_SI1))
plot(acf(resid(memRMR_MMR_SI1, type = "normalized"))) 
plot(ACF(mmemRMR_MMR_SI1), alpha=0.05)

# check assumptions - - MMR x time 
plot(resid(memRMR_MMR_SI2, type="normalized")~fitted(memRMR_MMR_SI2))
qqnorm(resid(memRMR_MMR_SI2))
plot(acf(resid(memRMR_MMR_SI2, type = "normalized"))) 
plot(ACF(memRMR_MMR_SI2), alpha=0.05)

# tidy up the fixed effects
SI_RMR_tidy <- tidy(memRMR_MMR_SI1, effects="fixed")
SI_RMR_tidy <- SI_RMR_tidy %>% 
  mutate(model="SMR x time")

SI_MMR_tidy <- tidy(memRMR_MMR_SI2, effects="fixed")
SI_MMR_tidy <- SI_MMR_tidy %>% 
  mutate(model="MMR x time")

SI_time_tidy <- rbind(SI_RMR_tidy, SI_MMR_tidy)

SI_time_tidy <- SI_time_tidy  %>% 
  relabel_predictors(c("(Intercept)"= "Intercept",
                       "PopulationRough" = "Non-Anadromous",
                       "poly(M_Num, 3)1" = "Time", 
                       "poly(M_Num, 3)2" = "Time^2", 
                       "poly(M_Num, 3)3" = "Time^3", 
                       "Length" = "Initial length",
                       "rRMR" = "rSMR",
                       "rMMR" ="rMMR",
                       "TemperatureHigh" = "Temperature: Warm",
                       "rRMR:poly(M_Num, 3)1" = "MR * Time",
                       "rRMR:poly(M_Num, 3)2" = "MR * Time^2",
                       "rRMR:poly(M_Num, 3)3" = "MR * Time^3",
                       "rMMR:poly(M_Num, 3)1" = "MR * Time",
                       "rMMR:poly(M_Num, 3)2" = "MR * Time^2",
                       "rMMR:poly(M_Num, 3)3" = "MR * Time^3",
                       "TemperatureHigh:rMMR" = "MMR * Warm",
                       "rRMR:TemperatureHigh" = "SMR * Warm",
                       "TemperatureHigh:poly(M_Num, 3)1" = "Warm * Time",
                       "TemperatureHigh:poly(M_Num, 3)2" = "Warm * Time^2",
                       "TemperatureHigh:poly(M_Num, 3)3" = "Warm * Time^3"

                       
  ))
# 

ordered_vars <- c("Intercept","rSMR", "rMMR","Temperature: Warm", "Non-Anadromous",
                  "Time","Time^2","Time^3", "Initial length",
                  "SMR * Warm", "MMR * Warm", 
                  "MR * Time","MR * Time^2","MR * Time^3", "Warm * Time","Warm * Time^2","Warm * Time^3"
                 )

## order all the  model dataframes 
SI_time_tidy <- SI_time_tidy   %>%
  mutate(term =  factor(term, levels = ordered_vars)) %>% 
  arrange(term)




# create supplementary infromation plot
MR_time_dw <- dwplot(SI_time_tidy,show_intercept =T, dot_args = list(size = 1.5), 
                   whisker_args = list(size = .75), dodge_size = 0.7)+
  theme_bw()+geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size=0.8) + xlab("\\nCoefficient estimate") +
  theme(text = element_text(size = 14), #should be 14 for combined figure
        legend.position = "top",
        # legend.position = c(0.975, .90), #for combined figure
        legend.justification=c(1, 1),
        legend.background = element_rect(colour="white"),
        legend.title=element_blank(),
        legend.text=element_text(size=14),
        axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=14, colour="black"), # #for combined figure
        axis.text = element_text(size=16, colour="black")) +
  scale_colour_manual(breaks=c("SMR x time", "MMR x time"),
                      values=c("black", "medium aquamarine")) +
  scale_y_discrete(breaks=c("Intercept","rSMR", "rMMR","Temperature: Warm", "Non-Anadromous",
                            "Time","Time^2","Time^3", "Initial length",
                            "SMR * Warm", "MMR * Warm", 
                            "MR * Time","MR * Time^2","MR * Time^3", "Warm * Time","Warm * Time^2","Warm * Time^3"),
                   labels=c("Intercept","rSMR", "rMMR","Temperature: Warm", "Non-Anadromous",
                            "Time",bquote(Time^2),bquote(Time^3),"Initial length",
                            "SMR * Warm", "MMR * Warm", 
                            "MR * Time",bquote("MR * "*Time^2),bquote("MR * "*Time^3),
                            "Warm * Time", bquote("Warm * "*Time^2),bquote("Warm * "*Time^3)))+
theme(plot.margin=unit(c(3,1,4,0),"cm"))  #for combined figure
 


############################
### SI Figures #############
############################


### for graphical purposes, categorise metabolic rates as high or low depending on residual values
MetLong <- MetLong %>% 
  mutate(RMR_Phenotype = ifelse(rRMR > 0, c("High_SMR"), "Low_SMR"),
         AS_Phenotype =  ifelse(rAS > 0, c("High_AS"), "Low_AS"),
         MMR_Phenotype = ifelse(rMMR > 0, c("High_MMR"), "Low_MMR"))

# For visualing, calculate growth residuals (correct for initial body size)
ModGR <- lm(Gs ~ Length, data=MetLong)
summary(ModGR)
plot(Gs~Length, data=MetLong)
MetLong$rGs <- resid(ModGR)


MetLong <- MetLong %>% mutate (time = as.factor(M_Num))

# calcualte the mean rGs for each measurement period according to phenotype
MetL_summary <- MetLong %>% 
  group_by(time, Temperature,  RMR_Phenotype) %>% 
  summarise(
    Gs_mean = mean(rGs, na.rm=T),
    Gs_sd = sd(rGs, na.rm=T),
    Gs_se = sd(rGs, na.rm=T)/sqrt(n()),
    Gs_ci = 1.96 * sd(rGs, na.rm=T)/sqrt(n()))


MetL_summary <- droplevels(MetL_summary)
MetL_summary$Phen_Temp <- as.factor(paste(MetL_summary$RMR_Phenotype, MetL_summary$Temperature, sep="_"))
MetLong$Phen_Temp <- as.factor(paste(MetLong$RMR_Phenotype, MetLong$Temperature, sep="_"))


RMR_time <- ggplot(data=MetL_summary, aes(x=time, y=Gs_mean, colour=Phen_Temp)) +
  geom_point(data=MetLong, aes(x=time, y=rGs, colour=Phen_Temp), alpha=0.2, size=2,fill="white", position = position_jitterdodge(dodge.width=.5, jitter.width=0))+
  geom_linerange(aes(ymin=Gs_mean-Gs_se, ymax=Gs_mean+Gs_se), position=position_dodge(0.5)) + 
  geom_point(aes(shape=Phen_Temp),size=3,fill="white", position = position_jitterdodge(dodge.width=.5, jitter.width=0)) +my_theme +
  ylim(c(-0.12, 0.25))+
  scale_fill_manual(labels=c("High SMR Warm", "High SMR Cool", "Low SMR Warm", "Low SMR Cool"),
                    breaks=c("High_SMR_High", "High_SMR_Low", "Low_SMR_High", "Low_SMR_Low"),
                    values = c("#e41a1c", "#377eb8", "#e41a1c", "#377eb8")) +
  scale_colour_manual(labels=c("High SMR Warm", "High SMR Cool", "Low SMR Warm", "Low SMR Cool"),
                      breaks=c("High_SMR_High", "High_SMR_Low", "Low_SMR_High", "Low_SMR_Low"),
                      values = c("#e41a1c", "#377eb8", "#e41a1c", "#377eb8"))  +
  scale_shape_manual(labels=c("High SMR Warm", "High SMR Cool", "Low SMR Warm", "Low SMR Cool"),
                     breaks=c("High_SMR_High", "High_SMR_Low", "Low_SMR_High", "Low_SMR_Low"),
                     values = c(16, 16, 24, 24)) +
  scale_x_discrete(labels=c("April-June", "June-July", "July-Sept", "Sept-Nov",
                            "Nov-Feb", "Feb-April"), 
                   breaks=c("3", "4.5", "6", "8", "10.5", "13"))+
  xlab("Measurement Period")+ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+labs(fill='Treatment', colour="Treatment", linetype="Treatment", shape="Treatment") +
  theme(legend.position=c(.5,.92), legend.title=element_blank(), legend.text=element_text(size=12),plot.margin=unit(c(1,1,0.5,1),"cm"),
        axis.text.x=element_text(angle=45, size=14,hjust=1),axis.title.y=element_text(vjust=2.5),axis.title.x=element_text(vjust=-1.5))+guides(col = guide_legend(ncol = 2))





MMR_summary <- MetLong %>% 
  group_by(time,Temperature, MMR_Phenotype) %>% 
  summarise(
    Gs_mean = mean(rGs, na.rm=T),
    Gs_sd = sd(rGs, na.rm=T),
    Gs_se = sd(rGs, na.rm=T)/sqrt(n()),
    Gs_ci = 1.96 * sd(rGs, na.rm=T)/sqrt(n()))


MMR_summary <- droplevels(MMR_summary)
MMR_summary$Phen_Temp <- as.factor(paste(MMR_summary$MMR_Phenotype, MMR_summary$Temperature, sep="_"))
MetLong$Phen_Temp <- as.factor(paste(MetLong$MMR_Phenotype, MetLong$Temperature, sep="_"))

MMR_time <- ggplot(data=MMR_summary, aes(x=time, y=Gs_mean, colour=Phen_Temp)) +
  geom_point(data=MetLong, aes(x=time, y=rGs, colour=Phen_Temp), alpha=0.2, size=2,fill="white", position = position_jitterdodge(dodge.width=.5, jitter.width=0))+
  geom_linerange(aes(ymin=Gs_mean-Gs_se, ymax=Gs_mean+Gs_se), position=position_dodge(0.5)) + 
  geom_point(aes(shape=Phen_Temp), size=3,fill="white", position = position_jitterdodge(dodge.width=.5, jitter.width=0)) +my_theme +
  ylim(c(-0.12, 0.25))+
  scale_fill_manual(labels=c("High MMR Warm", "High MMR Cool", "Low MMR Warm", "Low MMR Cool"),
                    breaks=c("High_MMR_High", "High_MMR_Low", "Low_MMR_High", "Low_MMR_Low"),
                    values = c("#e41a1c", "#377eb8", "#e41a1c", "#377eb8")) +
  scale_colour_manual(labels=c("High MMR Warm", "High MMR Cool", "Low MMR Warm", "Low MMR Cool"),
                      breaks=c("High_MMR_High", "High_MMR_Low", "Low_MMR_High", "Low_MMR_Low"),
                      values = c("#e41a1c", "#377eb8", "#e41a1c", "#377eb8"))  +
  scale_shape_manual(labels=c("High MMR Warm", "High MMR Cool", "Low MMR Warm", "Low MMR Cool"),
                     breaks=c("High_MMR_High", "High_MMR_Low", "Low_MMR_High", "Low_MMR_Low"),
                     values = c(16, 16, 24, 24)) +
  scale_x_discrete(labels=c("April-June", "June-July", "July-Sept", "Sept-Nov",
                            "Nov-Feb", "Feb-April"), 
                   breaks=c("3", "4.5", "6", "8", "10.5", "13"))+
  xlab("Measurement Period")+ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+labs(fill='Treatment', colour="Treatment", linetype="Treatment", shape="Treatment") +
  theme(legend.position=c(.5,.92), legend.title=element_blank(), legend.text=element_text(size=12),plot.margin=unit(c(1,1,1,1),"cm"),
        axis.text.x=element_text(angle=45, size=14,hjust=1),axis.title.y=element_text(vjust=2.5),axis.title.x=element_text(vjust=-1.5))+
  guides(col = guide_legend(ncol = 2))




# create multiplot of marginal effects
time_multi <- ggarrange(RMR_time, MMR_time ,labels=c("(B)", "(C)"),font.label=list(size=18, face="bold"),
                        label.x=c(0.2,0.2),label.y=.92,common.legend = F, nrow=2)


# ggarrange(MR_time_dw,time_multi ,labels=c("(A)"),font.label=list(size=18, face="bold"),
#           label.x=c(0.32, 0.23),label.y=c(.82),common.legend = F, ncol=2, widths=c(1, 1), heights=c(0.45, 3))




#####################################
######## Metabolic trait coupling####
#####################################


#### explore how temperature influences metabolic trait coupling

# use smatr package here because we don't have any a priori expectation about
# whether SMR or MMR drives one another

MMR_RMR <- sma(rRMR ~ rMMR*Temperature, data=metdat) # test for slope
summary(MMR_RMR)
plot(MMR_RMR)

MMR_RMR2 <- sma(rMMR ~ rRMR+Temperature, data=metdat) # test for intercept
summary(MMR_RMR2)
plot(MMR_RMR2)


# plot the trait coupling
MMR_RMR_corr2 <- ggplot(metdat, aes(x=rRMR, y=rMMR))+
  geom_point(size=3)+my_theme+ ylim(c(-.295, 0.295))+
  ylab(bquote('\\nrMMR (mg'*O[2]*'h'*r^-1*')'))  + xlab(bquote('\\nrSMR (mg'*O[2]*'h'*r^-1*')'))+ 
  theme(legend.position = c(0.51, 0.95), legend.title=element_blank(),legend.key.width = unit(1.25,"cm"), 
        legend.text=element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2.5), 
        axis.text=element_text(size=18), axis.title=element_text(size=28))

# dev.off()



###########################
### All fish growth #######
###########################

# Additional analysis to look at general trends in all fish across the study
# to check whether metabolic x temperature - growth patterns are reflect general
# temperature - growth patterns

# read in growth data for all study fish
GrowLong <- read.csv("Archer_et_al_all_fish_growth_PRSB.csv")

# reorder factor levels
GrowLong$Temperature <- factor(GrowLong$Temperature, levels=c("Low", "High"))
GrowLong$Population <- factor(GrowLong$Population, levels=c("Erriff_N", "Rough"))


##############################################
# mixed effect models for growth of all fish #
##############################################

## set the contorl for lme
ctrl <- lmeControl(opt='optim')


# model with temperature x time interaction
memGrowth1 <- lme(Gs ~  Temperature*poly(M_Num, 3)+Population+Length,
                  random=~1|DNA_ID, 
                  control=ctrl,
                  data=GrowLong,
                  correlation=corAR1(form=~M_Num, value =0.1),
                  na.action=na.omit,
                  method="ML")

# temperature x time interaction dropped
memGrowth2 <- lme(Gs ~  Temperature+Population+poly(M_Num, 3)+Length,
                  random=~1|DNA_ID, 
                  control=ctrl,
                  data=GrowLong,
                  correlation=corAR1(form=~M_Num, value =0.1),
                  na.action=na.omit,
                  method="ML")

anova(memGrowth1, memGrowth2) # test for temperature x time interaction


# model with temperature x time is best fit (memGrowth1)

# refit using REML
memGrowth <- lme(Gs ~  poly(M_Num, 3)*Temperature+Population+Length,
                  random=~1|DNA_ID, 
                  control=ctrl,
                  data=GrowLong,
                  correlation=corAR1(form=~M_Num, value=0.1),
                  na.action=na.omit,
                  method="REML")

Anova(memGrowth, test="LR")


# find the marginal r-sq
r.squaredGLMM(memGrowth)

# check model assumptions
plot(resid(memGrowth, type="normalized")~fitted(memGrowth))
qqnorm(resid(memGrowth))
plot(acf(resid(memGrowth, type = "normalized"))) 
plot(ACF(memGrowth), alpha=0.05)


# tidy up the fixed effects
Growth_tidy <- tidy(memGrowth, effects="fixed")

###############
# visualing growth trends

# do some summary stats and plots to visualise trends in growth. 

# first calculate residual growth rates (to adjust for initial body size)
GrowLong <- GrowLong  %>%
  mutate(time = as.factor(as.character(M_Num))) %>%
  group_by(time)  %>%
  subset(., !is.na(Gs))  %>%
  mutate(rGs = resid(lm(Gs ~ Length)))

GrowLong$time <- factor(GrowLong$time, levels=c("3", "4.5", "6", "8", "10.5"))

# the calculate the mean residual growth rate in each treatment across the study
GrowL_summary <- GrowLong %>% 
  group_by(time, Temperature) %>% 
  summarise(
    n = n(),
    Gs_mean = mean(rGs, na.rm=T),
    Gs_sd = sd(rGs, na.rm=T),
    Gs_se = sd(rGs, na.rm=T)/sqrt(n()),
    Gs_ci = 1.96 * sd(rGs, na.rm=T)/sqrt(n()))

GrowL_summary <- droplevels(GrowL_summary)

# plot the growth rate summaries
Gr_month <- ggplot(data=GrowL_summary, aes(x=time, y=Gs_mean, colour=Temperature)) + 
  geom_point(data=GrowLong, aes(x=time, y=rGs, colour=Temperature), alpha=0.15, size=2,fill="white", position = position_jitterdodge(dodge.width=.6, jitter.width=.5))+
  geom_linerange(aes(ymin=Gs_mean-Gs_se, ymax=Gs_mean+Gs_se), size=0.8, position=position_dodge(0.6)) + 
  geom_point(aes(shape=Temperature, fill=Temperature),size=4, position = position_jitterdodge(dodge.width=.6, jitter.width=0)) +my_theme +
  # ylim(c(-.1, 0.1))+
  ylim(c(min(GrowLong$rGs), max(GrowLong$rGs)))+
  scale_fill_manual(labels=c("Warm", "Cool"),
                    breaks=c("High",  "Low"),
                    values = c("#e41a1c",   "#377eb8")) +
  scale_colour_manual(labels=c("Warm", "Cool"),
                      breaks=c("High",  "Low"),
                      values = c("#e41a1c",  "#377eb8"))  +
  scale_shape_manual(labels=c("Warm", "Cool"),
                     breaks=c("High",  "Low"),
                     values = c(16,   24)) +
  scale_x_discrete(labels=c("April-June", "June-July", "July-Sept", "Sept-Nov",
                            "Nov-April"), 
                   breaks=c("3", "4.5", "6", "8", "10.5"))+
  xlab("Measurement period")+ylab(bquote('\\n Specific growth rate (% '*day^-1*')'))+labs(fill='Treatment', colour="Treatment", linetype="Treatment", shape="Treatment") +
  theme(legend.position=c(.84,.90), legend.title=element_blank(), legend.text=element_text(size=14),plot.margin=unit(c(0.15,0.1,0.5,1),"cm"),
        axis.text.x=element_text(angle=45, size=14,hjust=1),axis.title.y=element_text(vjust=2.5),axis.title.x=element_text(vjust=-1.5))

# create figure for model coefficients
Growth_tidy <- Growth_tidy %>% 
  relabel_predictors(c("(Intercept)"= "Intercept",
                       "poly(M_Num, 3)1" = "Time", 
                       "poly(M_Num, 3)2" = "Time^2", 
                       "poly(M_Num, 3)3" = "Time^3", 
                       "TemperatureHigh" = "Temperature: Warm",
                       "PopulationRough" = "Non-Anadromous",
                       "Length" = "Initial length",
                       "poly(M_Num, 3)1:TemperatureHigh" = "Time * Warm",
                       "poly(M_Num, 3)2:TemperatureHigh" = "Time^2 * Warm",
                       "poly(M_Num, 3)3:TemperatureHigh" = "Time^3 * Warm"
  ))


# create a figure ofmodel co-efficients
Growth_dw <- dwplot(Growth_tidy, dodge_size = 0.5,show_intercept =T, dot_args = list(size = 2, colour="black"),
                    whisker_args = list(size = 1, colour="black"),
                    order_vars=c("Intercept","Time", "Time^2","Time^3","Temperature: Warm","Non-Anadromous",
                                 "Initial length","Time * Warm","Time^2 * Warm","Time^3 * Warm")) +
  theme_bw()+geom_vline(xintercept = 0, colour = "grey60", linetype = 2, size=1) + xlab("\\nCoefficient estimate") +
  theme(text = element_text(size = 16),
        legend.position = "none",
        legend.justification=c(1, 1),
        legend.background = element_rect(colour="grey80"),
        axis.title.x = element_text(size=20),
        axis.text = element_text(size=16),
        axis.text.y = element_text(size=16)) +
  scale_y_discrete(breaks=c("Intercept","Time", "Time^2","Time^3","Temperature: Warm","Non-Anadromous",
                            "Initial length","Time * Warm","Time^2 * Warm","Time^3 * Warm"),
                   labels=c("Intercept","Time", bquote(Time^2),bquote(Time^3),"Temperature: Warm","Non-AB",
                            "Initial length",
                            "Time * Warm",bquote(Time^2*" * Warm") ,bquote(Time^3*" * Warm")))+
  theme(plot.margin=unit(c(0.15,1,1,0.5),"cm"), axis.text=element_text(colour="black"))

# create a multiplot of both (not used in manuscript)

# create a multi panel plot with the mean growth rates, and also with the dot-whisker plot for parameter co-efficients. 
# ggarrange(Gr_month, Growth_dw,labels=c("(A)", "(B)"),font.label=list(size=20),
#           label.x=c(0.18, .33),label.y=.98,common.legend = F, ncol=2, widths=c(0.9, 1))
# dev.off()


