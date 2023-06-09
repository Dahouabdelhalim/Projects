##### Strategies paper ProcB submission  - main Figures and stats #####
# This file makes the panels D,E,F in Figures 2 and 3 in the main manuscript and shows how the stats were done using a mixed model 


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

#### Figure 02, Panels D,E,F - with linear mixed models stats ####
carballoc<-read.csv("carballoc400ppm.csv")  #this is for the 400ppm data. Eleated Co2 data are below
head(carballoc)

MM1<-lme(fixed =muecomp.avg ~ percentage*react_to, random = ~ 1|ecotype, method = 'ML', carballoc)

dM<-dredge(MM1)
dM #yes, the treatments are different, so do this per treatment, although the lineages are less different now.. 


slopefn<-function(d) { 
  m<-lm (muecomp.avg~percentage,data=d) #can force through origin by setting offset to 0 but we don't here . would be offset=rep(0,length(d$mue)
  sum.<-summary(m) 
  r2<-sum.$r.squared
  intercept<-m$coefficients[1]
  slope<-m$coefficients[2]
  plargert<-sum.$coefficients[8]
  output<-data.frame(slope,intercept,r2,plargert)
  output}

slopedata<-ddply(carballoc,.(carballoc$react_to),slopefn)
slopedata



COC<-subset(carballoc,react_to=="COCULT")
SPIKE<-subset(carballoc,react_to=="SPIKES")
GFP<-subset(carballoc,react_to=="GFP")

MM_COC<-lme(fixed =muecomp.avg ~ percentage*ecotype, random = ~ 1|biorep, method = 'ML', COC)
dM<-dredge(MM_COC)
dM#yes,but no interaction 
#write.csv(prettify(dM),"modelsel_COCUL_alloc_400ppm.csv")

MM_COC_fin<-lme(fixed =muecomp.avg ~ percentage+ecotype, random = ~ 1|biorep, method = 'REML', COC)
#write.csv(prettify(summary(MM_COC_fin)),"modelout_COCUL_alloc_400ppm.csv")

COCUL_alloc_amb <- ggplot(COC, aes(percentage,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
COCUL_alloc_amb #rename x axis label into sth along the lines of "%of Âµgc not used for growth" to make it more intuitive. Also in the other two plots. 

MM_SPI<-lme(fixed =muecomp.avg ~ percentage*ecotype, random = ~ 1|biorep, method = 'ML', SPIKE)
dS<-dredge(MM_SPI)
dS#yes,but no interaction 
#write.csv(prettify(dS),"modelsel_SPIKE_alloc_400ppm.csv")

MM_SPI_fin<-lme(fixed =muecomp.avg ~ percentage+ecotype, random = ~ 1|biorep, method = 'REML', SPIKE)
#write.csv(prettify(summary(MM_SPI_fin)),"modelout_SPIKE_alloc_400ppm.csv")

SPIKE_alloc_amb <- ggplot(SPIKE, aes(percentage,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
SPIKE_alloc_amb

GFP<-subset(GFP, ecotype!="oth95") #something went wrong with the measurements of oth95 at 400ppm GFP

MM_GFP<-lme(fixed =muecomp.avg ~ percentage*ecotype, random = ~ 1|biorep, method = 'ML', GFP)
dG<-dredge(MM_GFP)
dG#yes,but no interaction 
#write.csv(prettify(dG),"modelsel_GFP_alloc_400ppm.csv")

MM_GFP_fin<-lme(fixed =muecomp.avg ~ percentage+ecotype, random = ~ 1|biorep, method = 'REML', GFP)
#write.csv(prettify(summary(MM_GFP_fin)),"modelout_GFP_alloc_400ppm.csv")

GFP_alloc_amb <- ggplot(GFP, aes(percentage,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
GFP_alloc_amb # all OK, realised during measurements already that larger SE in some than in others! 


plot_grid(COCUL_alloc_amb,SPIKE_alloc_amb,GFP_alloc_amb, ncol=3) #make it so that there is only ONE Figure legend. 


#####Figure 03, Panels D,E,F - with lineaer mixed models stats ####
carballoc<-read.csv("carballoc1000ppm.csv") # this is for the 1000ppm data. Ambient Co2 data is above . In an ideal world I would have named this carballoc1000 ....
head(carballoc)
MM1<-lme(fixed =muecomp.avg ~ surplus02*react_to, random = ~ 1|ecotype, method = 'ML', carballoc) # global model with full interaction
dM<-dredge(MM1)
dM #yes, the treatments are different, so do this per treatment


slopefn<-function(d) { 
  m<-lm (muecomp.avg~surplus02,data=d) #can force through origin by setting offset to 0 but we don't here . would be offset=rep(0,length(d$mue)
  sum.<-summary(m) 
  r2<-sum.$r.squared
  intercept<-m$coefficients[1]
  slope<-m$coefficients[2]
  plargert<-sum.$coefficients[8]
  output<-data.frame(slope,intercept,r2,plargert)
  output}

slopedata<-ddply(carballoc,.(carballoc$react_to),slopefn)
slopedata

COC<-subset(carballoc,react_to=="COCULT")
SPIKE<-subset(carballoc,react_to=="SPIKE")
GFP<-subset(carballoc,react_to=="GFP")

MM_all<-lme(fixed =muecomp.avg ~ surplus_percent_PS*ecotype, random = ~ 1|biorep, method = 'ML', carballoc)

dM<-dredge(MM_all)
dM 
#write.csv(prettify(dM),"modelsel_COCUL_alloc_1000ppm.csv")
#prettify is a godsend within 'PapeRs' - may need to get from github 

MM_COC<-lme(fixed =muecomp.avg ~ surplus02*ecotype, random = ~ 1|biorep, method = 'ML', COC)
dM<-dredge(MM_COC)
dM#yes,but no interaction 
#write.csv(prettify(dM),"modelsel_COCUL_alloc_1000ppm.csv")

MM_COC_fin<-lme(fixed =muecomp.avg ~ surplus02+ecotype, random = ~ 1|biorep, method = 'REML', COC)
#write.csv(prettify(summary(MM_COC_fin)),"modelout_COCUL_alloc_1000ppm.csv")

COCUL_alloc_high02 <- ggplot(COC, aes(surplus02,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
COCUL_alloc_high02 #need to spell out axes labels in a better more intuitive way 


MM_SPI<-lme(fixed =muecomp.avg ~ surplus02*ecotype, random = ~ 1|biorep, method = 'ML', SPIKE)
dS<-dredge(MM_SPI)
dS#yes,but no interaction 
#write.csv(prettify(dS),"modelsel_SPIKE_alloc_1000ppm.csv")

MM_SPI_fin<-lme(fixed =muecomp.avg ~ surplus02+ecotype, random = ~ 1|biorep, method = 'REML', SPIKE)
#write.csv(prettify(summary(MM_SPI_fin)),"modelout_SPIKE_alloc_1000ppm.csv")

SPIKE_alloc_high02 <- ggplot(SPIKE, aes(surplus02,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
SPIKE_alloc_high02 #need to spell out axes labels in a better more intuitive way 



MM_GFP<-lme(fixed =muecomp.avg ~ surplus02*ecotype, random = ~ 1|biorep, method = 'ML', GFP)
dG<-dredge(MM_GFP)
dG#yes,but no interaction 
#write.csv(prettify(dG),"modelsel_GFP_alloc_1000ppm.csv")

MM_GFP_fin<-lme(fixed =muecomp.avg ~ surplus02+ecotype, random = ~ 1|biorep, method = 'REML', GFP)
#write.csv(prettify(summary(MM_GFP_fin)),"modelout_GFP_alloc_1000ppm.csv")


GFP_alloc_high02 <- ggplot(GFP, aes(surplus02,muecomp.avg , colour = ecotype)) +
  geom_point(size=4)+
  theme_classic(base_size = 16, base_family = 'Helvetica')+
  scale_colour_ordinal()+
  geom_errorbar(aes(ymin=muecomp.avg-muecomp.se, ymax=muecomp.avg+muecomp.se))+
  geom_smooth(aes(group=1), method="lm",se=FALSE)+ #this is just to LOOK at it. no real stats of course. those we get from actually making the lm 
  theme(legend.position = "top",legend.box = "horizontal")+
  labs(x=expression(percentage~surplus~PS), y=expression(Fold~change~compared~to~monoculture)) 
GFP_alloc_high02 #same. need to make labels on axes so that they make intuitive sense

plot_grid(COCUL_alloc_high02,SPIKE_alloc_high02,GFP_alloc_high02, ncol=3) #need to adjust size so that legend does not vanish, but is OK


####### the rest is confetti ######
