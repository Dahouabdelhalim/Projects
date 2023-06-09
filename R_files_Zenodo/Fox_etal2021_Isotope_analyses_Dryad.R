######## Fox et al. 2021 - Differential resistance and acclimation of two coral species to chronic nutrient
########                   enrichment reflect life history traits 
######## Finalized: M. Fox 3/10/21

####### Code file 4 of 4 -- Stable isotope analyses for Figure 4
####### Required data files: 1) Isotopes_raw.csv
#                            2) Isotope nubbins.csv

####### Summary: This file completes all of the data summarization, analysis, and plotting for main text
#######          figure 4 (C:N vs. d15N)

####### Figures produced: Fig. 4 and Figures S10,11 (d13C, d15N, and C:N values)
library(Cairo)
library(dplyr)
library(ggplot2)
library(cowplot)
library(lmerTest)
library(RColorBrewer)
library(emmeans) # for tukey tests of lmer model 
library(car)
library(jtools)
library(rstatix)
library(broom)
library(ellipse)
library(tidyr)
library(lmodel2)

#create plotting themes
newtheme <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.line = element_line(size = 0))


nubs<-read.csv(file="Isotope_corals.csv",header=T)
topes<-read.csv(file="Isotopes_raw.csv",header=T)

#join the data based on the sample ID 
chain<- left_join(nubs, topes)

# take means and SE for supplemental figure
means<-chain %>% na.omit(chain) %>% #remove samples with NAs
  group_by(Species,Nutrient.Level,Tissue) %>% 
  summarize(
    N=sum(!is.na(d15N)),
    mean13C = mean(d13C,na.rm=T),
    mean15N = mean(d15N,na.rm=T),
    meanCN = mean(C.N,na.rm=T),
    se13C = sd(d13C)/sqrt(N),
    se15N = sd(d15N)/sqrt(N),
    seCN = sd(C.N)/sqrt(N))


# figure 4 - test the relationship between d15N and C:N ratios for each of the tissue compartments 
#            in each species. use mixed effects model by species with random effects of colony and tank

#POC
poc.mod<-lmer(C.N~Tissue*d15N+(1|Colony)+(1|Tank.Level),data=subset(chain,Species=="POC"))
summary(poc.mod) 
# Fixed effects:
#                Estimate  Std. Error df       t value Pr(>|t|)    
# Tissue       -0.6439     0.9529     54.0972  -0.676   0.5021    
# d15N          0.2675    0.2252    54.4403   1.188   0.2401    
# Tissue:d15N   0.6107     0.2649     53.5495   2.305   0.0251 * 

#model assumptions
qqPlot(residuals(poc.mod))
plot(poc.mod)
poc.dat<-chain %>% filter(Species=="POC")
leveneTest(poc.dat$C.N~poc.dat$Tissue)

#isolate poc symbiont data, predict model fit
poc.sym<-poc.dat %>% dplyr::select(-d13C) %>%filter(Tissue=="Z") %>%  na.omit()

#model with no tissue effect
poc.sym.mod<-lmer(C.N~d15N+(1|Colony)+(1|Tank.Level),data=poc.sym)
symE<-effect_plot(poc.sym.mod,pred=d15N,interval=T,plot.points=T)
symE

#pull out components of the effects plot
p1<- ggplot_build(symE)

#extract data from the plot -- specifically dfs 2 (line) and 3(ci) for the model fits
sym.fit<-as.data.frame(p1$data[2])
sym.fit$Species<-"POC"
sym.ci<-as.data.frame(p1$data[3])

###model 2 regression
lm2<-lmodel2(C.N~d15N, data = poc.sym)
lm2
plot(lm2,method="OLS")
simple.mod<-lm(C.N~d15N, data = poc.sym)
summary(simple.mod)

#compare fit lines from mixed effects model, simple linear, and standardized major axis (SMA)
# for accurate slope estimate SMA is preferrable because of error in both x and y dimensions
# but model fit for the figure wil be shown by mixed effects model for consistency with 
# non-significant fits. 

ggplot()+
  geom_smooth(data=poc.sym,mapping=aes(x=d15N,y=C.N),method=lm,col="black",se=F)+
  geom_line(data=sym.fit,mapping=aes(x=x,y=y),col="blue")+
  geom_point(data=poc.sym,mapping=aes(d15N,y=C.N))+
  geom_abline(intercept=lm2$regression.results$Intercept[3],
                                        slope=lm2$regression.results$Slope[3],
                                        colour="red")


#POR
por.mod<-lmer(C.N~Tissue*d15N+(1|Colony)+(1|Tank.Level),data=subset(chain,Species=="POR"))
summary(por.mod) 

# Fixed effects:
#                Estimate  Std. Error       df   t value Pr(>|t|)    
#   Tissue      -1.44089    0.47612    53.19063  -3.026  0.00381 ** 
#   d15N        -0.06179    0.10277    56.90408  -0.601  0.55007    
# Tissue:d15N    0.01483    0.15750    51.03165   0.094  0.92536 

#model assumptions
qqPlot(residuals(por.mod))
plot(por.mod)
por<-chain %>% filter(Species=="POR")
leveneTest(por$C.N~por$Tissue)

#plotting
#colors for the nutrient levels
nutrient.cols=rev(c("firebrick3","lightseagreen","navy"))
Tshapes=c(21,24) #tissue shapes
chain$Species<-factor(chain$Species,levels=c("POR","POC")) #order species
sym.fit$Species<-factor(sym.fit$Species,levels=c("POR","POC")) #order species
labels<-c(POR="Porites compressa",POC="Pocillopora acuta") #set labels


# FIGURE 4
fig.4=  ggplot()+
  geom_point(data=chain,aes(x=d15N,y=C.N,fill=Nutrient.Level,shape=Tissue,group=Tissue),size=3,alpha=1,color='black')+
    scale_shape_manual(values=Tshapes)+
    scale_fill_manual(values=nutrient.cols)+
    #50% DATA ELLIPSES FOR EACH GROUP
    stat_ellipse(data=subset(chain,Tissue=="T"),mapping=aes(x=d15N,y=C.N,fill=Nutrient.Level,group=Nutrient.Level),color="black",alpha=0.5,geom = "polygon",level=.50)+
    stat_ellipse(data=subset(chain,Tissue!="T"),mapping=aes(x=d15N,y=C.N,fill=Nutrient.Level,group=Nutrient.Level),color="black",alpha=0.5,geom = "polygon",level=.50)+
    geom_smooth(data=subset(chain,!(Species=="POC" & Tissue=="Z")),mapping=aes(x=d15N,y=C.N,group=Tissue),method="lm",lty=2,color="black",size=0.5,se=F)+
    ##plot fitted line for pocillopora symbionts using fit from the mixed model (same as regression fit)
  geom_line(data=sym.fit,mapping=aes(x=x,y=y),color="black",size=1)+
    ## plot the fit line from model 2 regression for improved slope estimate for significant poc symbiont fit
  geom_abline(intercept=lm2$Intercept[3],
              slope=lm2$Slope[3],
              colour="red")+
  xlab(expression({delta}^15*N~'\\u2030'))+
  scale_y_continuous(breaks=seq(5,9,1))+expand_limits(y=c(5,9))+
  scale_x_continuous(breaks=seq(0,5,1.25))+
  ylab("C:N")+
  facet_grid(.~Species,labeller=labeller(Species=labels))+
    newtheme+theme(strip.text = element_text(face = "italic"))
fig.4

 
#mean plots - for supplemental figure
cols<-c("navajowhite3","seagreen3") #set colors
labels<-c(POC="Pocillopora acuta",POR="Porites compressa") #set labels
shapes<-c(21,24)
means$Species<-ordered(means$Species,levels=c("POR","POC"))


# d15N
a<-ggplot(means,aes(x=Nutrient.Level,y=mean15N,ymin=mean15N-se15N*1.96,ymax=mean15N+se15N*1.96,fill=Tissue,group=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue,shape=Tissue),color="black",size=2.5,position=position_dodge(0.9))+
  ylab(expression({delta}^15*N~'\\u2030'))+
  scale_x_discrete(labels=c("N0" = "0.1", "N2" = "3.0","N4" = "7.0"))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  scale_fill_manual(values=cols)+scale_shape_manual(values=shapes)+
  newtheme+theme(strip.text = element_text(face = "italic"))

a

#d13C
b<-ggplot(means,aes(x=Nutrient.Level,y=mean13C,ymin=mean13C-se13C*1.96,ymax=mean13C+se13C*1.96,group=Tissue,fill=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue,shape=Tissue),color="black",size=2.5,position=position_dodge(0.9))+
  ylab(expression({delta}^13*C~'\\u2030'))+
  scale_x_discrete(labels=c("N0" = "0.1", "N2" = "3.0","N4" = "7.0"))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  scale_y_continuous(breaks=seq(-19,-16,by=1))+expand_limits(y = c(-19, -16))+
  scale_fill_manual(values=cols)+scale_shape_manual(values=shapes)+
  newtheme+theme(strip.text = element_text(face = "italic"))

b

# CN
c<-ggplot(means,aes(x=Nutrient.Level,y=meanCN,ymin=meanCN-seCN*1.96,ymax=meanCN+seCN*1.96,group=Tissue,fill=Tissue,facets=Species))+
  geom_errorbar(position=position_dodge(0.9),width=0.0)+facet_wrap(~Species)+facet_grid(.~Species,labeller=labeller(Species=labels))+
  geom_point(aes(fill=Tissue,shape=Tissue),color="black",size=2.5,position=position_dodge(0.9))+
  ylab("C:N")+
  ylim(5,8.5)+
  scale_x_discrete(labels=c("N0" = "0.1", "N2" = "3.0","N4" = "7.0"))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  scale_fill_manual(values=cols)+scale_shape_manual(values=shapes)+
  newtheme+theme(strip.text = element_text(face = "italic"))
 c


#put the plots together to create main text Fig.4 
d<-a+b+c+plot_layout(ncol=1)+plot_annotation(tag_levels = "A")
#quartz()
d 


### run separate models for each species * tissue --- these are reported in the supplemental table S3
###################d15N
#POR
# 
d15N.porT<-lmer(d15N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),data=subset(chain,Species=="POR"))
summary(d15N.porT)
anova(d15N.porT)                      #tissue effect
sum(resid(d15N.porT)^2) #22.89

emm = emmeans(d15N.porT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")

#POC
d15N.pocT<-lmer(d15N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POC"))
summary(d15N.pocT)
anova(d15N.pocT)
sum(resid(d15N.pocT)^2) #1O.92

emm = emmeans(d15N.pocT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")
###################d13C

#POR
d13C.porT<-lmer(d13C~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POR"))
summary(d13C.porT)
anova(d13C.porT)
sum(resid(d13C.porT)^2) #12.18
emm = emmeans(d13C.porT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")
#POC
d13C.pocT<-lmer(d13C~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
                data=subset(chain,Species=="POC"))
summary(d13C.pocT)
anova(d13C.pocT)
sum(resid(d13C.pocT)^2) #9.31
emm = emmeans(d13C.pocT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")


###################CN
#POR
CN.porT<-lmer(C.N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
              data=subset(chain,Species=="POR"))
summary(CN.porT)
anova(CN.porT)
sum(resid(CN.porT)^2) #5.92
emm = emmeans(CN.porT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")

#POC
CN.pocT<-lmer(C.N~Nutrient.Level*Tissue+(1|Tank.Level)+(1|Colony),
              data=subset(chain,Species=="POC"))
summary(CN.pocT)
anova(CN.pocT)
sum(resid(CN.pocT)^2) #18.25
emm = emmeans(CN.pocT, ~ Nutrient.Level*Tissue)
pairs(emm,adjust="holm")

# Fig. S11 - differences between host and symbiont d13C and d15N values 
#calculate Dd13C to look at change in carbon sharing between host and symbiont
t<-chain %>% filter(Tissue=='T')
z<-chain %>% filter(Tissue=="Z")
tz<-merge(t,z,by="Nub_ID")
tz<- tz %>% mutate(Dd13C = d13C.x-d13C.y)
tz<- tz %>% mutate(Dd15N = d15N.x-d15N.y)

#poc d13C
tz.modC.poc<-lmer(Dd13C~Nutrient.Level.x+(1|Tank.Level.x)+(1|Colony.x),data=subset(tz,Species.x=="POC"))
anova(tz.modC.poc)
emm = emmeans(tz.modC.poc, ~ Nutrient.Level.x)
pairs(emm)
                  #Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Nutrient.Level.x 1.3601 0.68004     2 26.283  7.0799 0.00347 ** N2 and N4 different from0

#poc d15N
tz.modN.poc<-lmer(Dd15N~Nutrient.Level.x+(1|Tank.Level.x)+(1|Colony.x),data=subset(tz,Species.x=="POC"))
anova(tz.modN.poc)
emm = emmeans(tz.modN.poc, ~ Nutrient.Level.x)
pairs(emm)
# Nutrient.Level.x 6.4585  3.2292     2 18.462  15.624 0.0001069 *** N4 different from 2 and 0

#porites no effect of nutrients
tz.modC.por<-lmer(Dd13C~Nutrient.Level.x+(1|Tank.Level.x)+(1|Colony.x),data=subset(tz,Species.x=="POR"))
anova(tz.modC.por) # p = 0.92

tz.modN.por<-lmer(Dd15N~Nutrient.Level.x+(1|Tank.Level.x)+(1|Colony.x),data=subset(tz,Species.x=="POR"))
anova(tz.modN.por) # p =0.50

#means and CI for plotting

tz.mean<-tz %>% group_by(Species.x,Nutrient.Level.x) %>% 
  summarize(
    N=sum(!is.na(Dd13C)),
    meanD13C = mean(Dd13C,na.rm=T),
    meanD15N = mean(Dd15N,na.rm=T),
    se13C = sd(Dd13C,na.rm=T)/sqrt(N),
    se15N = sd(Dd15N,na.rm=T)/sqrt(N))

tz.plotC<-ggplot()+
  #geom_point(data=tz,aes(Nutrient.Level.x,Dd13C),alpha=0.5)+
  geom_hline(yintercept=0,lty=2)+facet_wrap(~Species.x)+facet_grid(.~Species.x,labeller=labeller(Species.x=labels))+
  geom_point(data=tz.mean,aes(Nutrient.Level.x,meanD13C))+
  geom_errorbar(data=tz.mean,mapping=aes(Nutrient.Level.x,ymin=meanD13C-se13C*1.96,ymax=meanD13C+se13C*1.96),width=0)+
  scale_x_discrete(labels=c("N0" = "0.1", "N2" = "3.0","N4" = "7.0"))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  newtheme+theme(strip.text = element_text(face = "italic"))
tz.plotC

tz.plotN<-ggplot()+
  #geom_point(data=tz,aes(Nutrient.Level.x,Dd13C),alpha=0.5)+
  geom_hline(yintercept=0,lty=2)+facet_grid(.~Species.x,labeller=labeller(Species.x=labels))+
  geom_point(data=tz.mean,aes(Nutrient.Level.x,meanD15N))+
  geom_errorbar(data=tz.mean,mapping=aes(Nutrient.Level.x,ymin=meanD15N-se15N*1.96,ymax=meanD15N+se15N*1.96),width=0)+
  scale_x_discrete(labels=c("N0" = "0.1", "N2" = "3.0","N4" = "7.0"))+
  xlab(bquote('Target Nitrate Level'~(µmol~L^-1)))+
  newtheme+theme(strip.text = element_text(face = "italic"))
tz.plotN

#supplemental TZ plot 
tz.supp<-tz.plotC+tz.plotN+plot_layout(ncol =1)+plot_annotation(tag_levels = 'A')
tz.supp






