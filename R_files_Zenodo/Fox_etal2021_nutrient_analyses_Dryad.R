######## Fox et al. 2021 - Differential resistance and acclimation of two coral species to chronic nutrient
########                   enrichment reflect life history traits 

####### Code file 2 of 4
####### Required data files: 1) Kaneohe_historical_nutrients.csv
#                            2) Kbay_nutrients_experimental.csv

####### Summary: This file analyzes historical nutrient data in Kaneohe Bay (from Drupp et al. 2011) along with
#######          the weekly nutrient water chemistry measurements and from all tanks during the experiment.

####### Figures produced: Figs S2,3,4 

library(dplyr)
library(ggplot2)
library(ggstance)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(lme4)
library(emmeans)
library(lmerTest)
library(cowplot)

### figure S2 - generate baseline nutrient conditions for Kbay based on data from Table 1 in Drupp et al. 2011
####            Citation:  Drupp, P., De Carlo, E.H., Mackenzie, F.T. et al. Nutrient Inputs, Phytoplankton Response, and CO2 Variations in a Semi-Enclosed Subtropical Embayment, Kaneohe Bay, Hawaii. 
#                          Aquat Geochem 17, 473–498 (2011). https://doi.org/10.1007/s10498-010-9115-y

#set plotting theme
newtheme <- theme_classic() + theme(text = element_text(size=11))+
  theme(axis.text.x = element_text(size=12,colour="black"), axis.text.y = element_text(size=12,colour="black"))+
  theme(plot.margin = unit(c(5.5,5.5,5.5,20), "pt"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.line = element_line(size = 0))

#load raw data file
drupp<-read.csv(file="Kaneohe_historical_nutrients.csv",header=T)

#calculate means and SD for annotating the plots
means<-drupp %>% summarize(N=sum(!is.na(NO2.NO3)),
                          N.N = mean(NO2.NO3,na.rm=T),
                          P = mean(PO4,na.rm=T),
                          NP = mean(N.P,na.rm=T),
                          sdN = sd(NO2.NO3,na.rm=T),
                          sdP = sd(PO4,na.rm=T),
                          sdNP = sd(N.P,na.rm=T))


#create flipped density plots
N.plot<-ggplot(subset(drupp,NO2.NO3>0),aes(x=NO2.NO3,y=-0.05,group=1))+
  geom_point()+
  geom_boxploth(fill="gray",width=0.05)+
  geom_density(aes(x=NO2.NO3,y=..scaled..),fill="steelblue3",alpha=0.75,inherit.aes = FALSE,adjust=0.6)+
  scale_x_log10(breaks=c(0.01,0.02,0.03,0.04,0.05,0.06,0.1,0.2,0.3,0.4,0.5,0.7,1,2,3),lim=c(0.01,3))+
  xlab(bquote(NO[2]+NO[3]~(µmol~L^-1)))+ylab("")+
  annotate(geom='text',x = 2.5, y = .45, label="Mean = 0.289, SD = 0.270")+
  coord_flip()+newtheme
N.plot

P.plot<-
  ggplot(subset(drupp,PO4>0),aes(x=PO4,y=-0.05,group=1))+
  geom_point()+
  geom_boxploth(fill="gray",width=0.05)+
  geom_density(aes(x=PO4,y=..scaled..),fill="steelblue3",alpha=0.75,inherit.aes = FALSE,adjust=0.6)+
  scale_x_log10(breaks=c(0.01,0.02,0.03,0.04,0.05,0.06,0.1,0.2,0.3,0.4,0.5,0.6),lim=c(0.01,0.6))+
  xlab(bquote(PO[4]~(µmol~L^-1)))+ylab("")+
  annotate(geom='text',x = 0.5, y = .45, label="Mean = 0.156, SD = 0.187")+
  coord_flip()+newtheme
P.plot

NP.plot<-ggplot(subset(drupp,N.P>0),aes(x=N.P,y=-0.05,group=1))+
  geom_point()+
  geom_boxploth(fill="gray",width=0.05)+
  geom_density(aes(x=N.P,y=..scaled..),fill="steelblue3",alpha=0.75,inherit.aes = FALSE,adjust=0.6)+
  scale_x_log10(breaks=c(0.3,0.4,0.5,0.6,1,2,3,4,5,7,10,20,30),lim=c(0.2,30))+
  xlab("N:P")+ylab("Data Density")+
  annotate(geom='text',x = 25, y = .45, label="Mean = 3.095, SD = 3.626")+
  coord_flip()+newtheme
NP.plot

##Figure S2
fig.S2<-N.plot+P.plot+NP.plot+plot_layout(ncol=1)+plot_annotation(tag_levels = 'A')
fig.S2



########################
########################
########################
###### Figures S3 and S4 along with statistics for changes in water chemistry parameters during the experiment
########################
########################

#load nutrient and water chemistry data from the experimental tanks 
clean<-read.csv(file="Kbay_nutrients_experimental_clean.csv",header=TRUE)


### take means by date and tank and calculate SE for each treatment
means<-clean%>%
    group_by(Days,Treatment) %>%   #### Define the grouping variables use summarise_each to apply to each column
    summarise_all(funs(mean(., na.rm = TRUE),sd(.,na.rm=TRUE),se=sd(.)/sqrt(n())))


#######PLOTTING THE VARIABLES

### make line plots for each parameter through time
## Fig. S3 -- Nitrate and Phosphate conditions over the course of the experiment
a<-ggplot(means,aes(x=Days,y=N.N_mean,ymin=N.N_mean-N.N_se,ymax=N.N_mean+N.N_se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  ylab(bquote('Nitrate+Nitrite ('*µmol~L^-1*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
a

a2<-ggplot(means,aes(x=Days,y=N.N_mean,ymin=N.N_mean-N.N_se,ymax=N.N_mean+N.N_se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  ylab(bquote('Nitrate+Nitrite ('*µmol~L^-1*')'))+theme(legend.position="bottom")+theme(axis.title.y=element_text(size=10))

b<-ggplot(means,aes(x=Days,y=Phosphate_mean,ymin=Phosphate_mean-Phosphate_se,ymax=Phosphate_mean+Phosphate_se,color=Treatment))+
  geom_errorbar(width=0.1)+expand_limits(y = c(0, 3))+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  ylab(bquote('Phosphate ('*µmol~L^-1*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
b

##### Make Fig. S3 -
leg<- get_legend(a2 + theme(legend.direction = "horizontal",legend.justification="center" ,legend.box.just = "bottom"))
Fig.S3<-plot_grid(a,b,leg,ncol=1,align='v',axis='1',labels=c("A","B"),rel_heights = c(1,1,.1))
Fig.S3


#### Make Fig S4 - other water chemistry parameters over time
#DON
c<-ggplot(means,aes(x=Days,y=DON_mean,ymin=DON_mean-DON_se,ymax=DON_mean+DON_se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+theme(axis.title.x=element_blank())+
  ylab("DON")+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
c

#DOP
d<-ggplot(means,aes(x=Days,y=DOP_mean,ymin=DOP_mean-DOP_se,ymax=DOP_mean+DOP_se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+theme(axis.title.x=element_blank())+
  ylab("DOP")+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
d

#FDOM A
e<-ggplot(means,aes(x=Days,y=Log10.CobleA.of.Sheet1._mean,ymin=Log10.CobleA.of.Sheet1._mean-Log10.CobleA.of.Sheet1._se,ymax=Log10.CobleA.of.Sheet1._mean+Log10.CobleA.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM A ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
e


#FDOM M
f<-ggplot(means,aes(x=Days,y=Log10.CobleM.of.Sheet1._mean,ymin=Log10.CobleM.of.Sheet1._mean-Log10.CobleM.of.Sheet1._se,ymax=Log10.CobleM.of.Sheet1._mean+Log10.CobleM.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM M ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
f

#fDOM C
g<-ggplot(means,aes(x=Days,y=Log10.CobleC.of.Sheet1._mean,ymin=Log10.CobleC.of.Sheet1._mean-Log10.CobleC.of.Sheet1._se,ymax=Log10.CobleC.of.Sheet1._mean+Log10.CobleC.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM C ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
g

#fDOM B
h<-ggplot(means,aes(x=Days,y=Log10.CobleB.of.Sheet1._mean,ymin=Log10.CobleB.of.Sheet1._mean-Log10.CobleB.of.Sheet1._se,ymax=Log10.CobleB.of.Sheet1._mean+Log10.CobleB.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM B ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
h

#fDOM T
i<-ggplot(means,aes(x=Days,y=Log10.CobleT.of.Sheet1._mean,ymin=Log10.CobleT.of.Sheet1._mean-Log10.CobleT.of.Sheet1._se,ymax=Log10.CobleT.of.Sheet1._mean+Log10.CobleT.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM T ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
i

# HIX
j<-ggplot(means,aes(x=Days,y=Log10.HIX.of.Sheet1._mean,ymin=Log10.HIX.of.Sheet1._mean-Log10.HIX.of.Sheet1._se,ymax=Log10.HIX.of.Sheet1._mean+Log10.HIX.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+ theme(axis.title.x=element_blank())+
  ylab(bquote('fDOM HIX ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
j

# FI
k<-ggplot(means,aes(x=Days,y=Log10.FI.of.Sheet1._mean,ymin=Log10.FI.of.Sheet1._mean-Log10.FI.of.Sheet1._se,ymax=Log10.FI.of.Sheet1._mean+Log10.FI.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  ylab(bquote('fDOM FI ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
k

#M:C
l<-ggplot(means,aes(x=Days,y=Log10.M.C.of.Sheet1._mean,ymin=Log10.M.C.of.Sheet1._mean-Log10.M.C.of.Sheet1._se,ymax=Log10.M.C.of.Sheet1._mean+Log10.M.C.of.Sheet1._se,color=Treatment))+
  geom_errorbar(width=0.1)+
  geom_point(aes(fill=Treatment),colour="black",pch=21, size=3)+geom_line(aes(group=Treatment),lty=2)+
  scale_color_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  scale_fill_manual(values=c("navy","slateblue3","lightseagreen","darkorange2","firebrick3"))+
  ylab(bquote('fDOM M:C ('*log[10]~R.U.*')'))+theme(legend.position="none")+theme(axis.title.y=element_text(size=10))
l


### Make Fig. s4
Fig.S4<-plot_grid(c,d,e,f,g,h,i,j,k,l,align="v",axis='1',ncol=2,labels= "AUTO",rel_heights = c(2,2,2,2,2,2,2,2,2,2),
                  rel_widths = c(.5,.5,.5,.5,.5,.5,.5,.5,.5,.5))
quartz()
Fig.S4




############## STATS for all nutrient and water chemistry parameters 
#################################################################################################################################
##### linear mixed effects model for each parameter to test for the effects of time, treatment, and their interaction
##### tank is included as an orthogonal random factor 

###############
## Nitrate + Nitrite
N.mod<-lmer(N.N~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(N.mod) 
anova(N.mod,type=2) #treatment = <0.0001, days = 0.001 *but see emm below
                    # interaction = 0.09

ranova(N.mod) #check Tank effect - NS
emm = emmeans(N.mod, ~ Days) #pairwise contrasts of days
pairs(emm) # difference due to days 9 vs. 21 slight increase

###############
######Phosphate
P.mod<-lmer(Phosphate~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(P.mod)
anova(P.mod,type=2) #treatment = <0.0001, days = 0.0001 
                    # interaction = 0.08

ranova(P.mod) #check Tank effect - NS
emm = emmeans(P.mod, ~ Days)
pairs(emm) # difference only due to days 9 and 14 vs. 21 and primarily due to small sd in N4 treatment on d21.

###############
###### DON - minor temporal variation towards end of exp but no treatment effects or interaction

DON.mod<-lmer(DON~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(DON.mod)
anova(DON.mod,type=2) #treatment = 0.36, days = <0.0001 
# interaction = 0.25

ranova(DON.mod) #check Tank effect - NS
emm = emmeans(DON.mod, ~ Days)
pairs(emm) 


###############
###### DOP
DOP.mod<-lmer(DOP~Days*Treatment+(1|Tank),data=clean,REML=FALSE)
summary(DOP.mod)
anova(DOP.mod,type=2) #treatment = 0.06, days = <0.0001 
# interaction = 0.23

ranova(DOP.mod) #check Tank effect - NS
emm = emmeans(DOP.mod, ~ Days)
pairs(emm) 


###############
###### Fdom A
A.mod<-lmer(Log10.CobleA.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(A.mod)
anova(A.mod,type=2) #treatment = 0.99, days = <0.0001 
# interaction = 0.93

ranova(A.mod) #check Tank effect - NS
emm = emmeans(A.mod, ~ Days)
pairs(emm) 

###############
###### Fdom M
M.mod<-lmer(Log10.CobleM.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(M.mod)
anova(M.mod,type=2) #treatment = 0.99, days = 0.0001 
# interaction = 0.95

ranova(M.mod) #check Tank effect - NS
emm = emmeans(M.mod, ~ Days)
pairs(emm) 

###############
###### Fdom C
C.mod<-lmer(Log10.CobleC.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(C.mod)
anova(C.mod,type=2) #treatment = 0.97, days = 0.0001
# interaction = 0.92

ranova(C.mod) #check Tank effect - NS
emm = emmeans(C.mod, ~ Days)
pairs(emm) 

###############
###### Fdom B
B.mod<-lmer(Log10.CobleB.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(B.mod)
anova(B.mod,type=2) #treatment = 0.88, days = 0.001 
# interaction = 0.31

ranova(B.mod) #check Tank effect - NS

emm = emmeans(B.mod, ~ Days)
pairs(emm) 

###############
###### Fdom T
T.mod<-lmer(Log10.CobleT.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(T.mod)
anova(T.mod,type=2) #treatment = 0.63, days = <0.0001 
# interaction = 0.29

ranova(T.mod) #check Tank effect - NS
emm = emmeans(T.mod, ~ Days)
pairs(emm) 

###############
###### HIX
HIX.mod<-lmer(Log10.HIX.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(HIX.mod)
anova(HIX.mod,type=2) #treatment = 0.50, days = <0.0001
# interaction = 0.02

ranova(HIX.mod) #check Tank effect - NS
emm = emmeans(HIX.mod, ~ Days)
pairs(emm) 

###############
###### M:C
MC.mod<-lmer(Log10.M.C.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(MC.mod)
anova(MC.mod,type=2) #treatment = 0.97, days = 0.06
# interaction = 0.96

ranova(MC.mod) #check Tank effect - NS
emm = emmeans(MC.mod, ~ Days)
pairs(emm) 


###############
###### FI
FI.mod<-lmer(Log10.FI.of.Sheet1.~Treatment*Days+(1|Tank),data=clean,REML=FALSE)
summary(FI.mod)
anova(FI.mod,type=2) #treatment = 0.12, days = 0.0078
# interaction = 0.94

ranova(FI.mod) #check Tank effect - NS
emm = emmeans(FI.mod, ~ Days)
pairs(emm) 
#################################################################################################################################


###########################################
#####extract DOC and DOC:DON data for day 14
###########################################
doc<-clean<-subset(clean,Days == "14d")

#doc mod (no time just treatment)
DOC.mod<-lmer(DOC..uM.~Treatment+(1|Tank),data=doc,REML=FALSE)
summary(DOC.mod)
anova(DOC.mod,type=2) #treatment = 0.70

#Doc:DON mod 
DOC.N.mod<-lmer(DOC.DON~Treatment+(1|Tank),data=doc,REML=FALSE)
summary(DOC.N.mod)
anova(DOC.N.mod,type=2) #treatment = 0.38
###########################################
###########################################




