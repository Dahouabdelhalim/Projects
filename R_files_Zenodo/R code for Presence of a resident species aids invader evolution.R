#This R code allows for the replication of the stats and figures presented in "Presence of a resident species aids invader evolution" 
#Authors: Josianne Lachapelle1, Elvire Bestion2,3, Eleanor E Jackson3,4 , C- Elisa Schaum3,5

#Author affiliations:
#1 Department of Biology, University of Toronto at Mississauga, 3359 Mississauga Road, William G. Davis Building, Mississauga, ON, L5L 1C6, Canada
#2 Station d’Ecologie Théorique et Expérimentale, UAR 2029, CNRS, Moulis, 09200, France, France
#3 Environment and Sustainability Institute, University of Exeter, Penryn Campus, Penryn, Cornwall TR10 9EZ, UK.
#4 School of Biological Sciences, University of Reading, Reading, UK
#5 Institute of Marine Ecosystem and Fishery Science, Centre for Earth System Research and Sustainability, Universität Hamburg, Olbersweg 24, 22767 Hamburg  

#should any bit of the code misbehave for you, or should you require (parts of ) the raw raw data (e.g. raw flow cytometre files) please get in touch with elisa.schaum@uni-hamburg.de
#new salt analysis. 


library(data.table)
library(plyr)
library(zoo)
library(lattice)
library(ggplot2)
require(nlme)
require(minpack.lm)
library(reshape2)
#install.packages("devtools")
library(devtools)
library(AICcmodavg)
library(lme4)
library(cowplot)
library(emmeans)
library(MuMIn)
library(nlsMicrobio)
library(papeR)

#remember lsmeans(lme_change_finC, pairwise ~ growthT, adjust = "tukey") may now be emmeans! 

rm(list=ls())
#read in the trajectory file
setwd("~/xx/xx/ stuff") # 

##### Trajectories biomass and cellsize  ####

mgr<-read.csv("means_transfer.csv")
#id may need adjusted
mgr <- within(mgr, id<- as.character(factor(temp):factor(species):factor(Invasion):factor(Habitat)))

#away (i.e. in novel salinity) trajectory of rate of incease-- should only use where transf!=0
mgr<-subset(mgr, temp!="22/32 slow")
mgr$temp<-factor( mgr$temp,  levels = c("22","26","22/32 rapid","32"))

salt3<-qplot( transf, growth.use,data=subset(mgr, Habitat=='away' & transf>0),ylab="Groth rate µ (day -1)",xlab="transfer number", colour=Invasion, facets=temp~species)+
  theme_classic(base_size = 14, base_family = 'Helvetica')+
  facet_grid(species~temp, scales="free") +
  geom_errorbar(aes(ymin=growth.use- gr.se, ymax=growth.use+ gr.se)) +
  
  #geom_smooth ( aes(group=id, colour=Invasion, fill=Invasion),method="gam",formula= y ~ s(x, k = 5)) +
  scale_colour_manual(values=c("red","goldenrod1","red","goldenrod1"))+
  scale_fill_manual(values=c("red","goldenrod1","red","goldenrod1"))+
  theme(legend.position = "top")
salt3

#home trajectory of rate of increase -- should only use where transf!=0
salt4<-qplot( transf, growth.use,data=subset(mgr, Habitat=='home' & transf>0),ylab="Groth rate µ (day -1)",xlab="transfer number", colour=Invasion, facets=temp~species)+
  theme_classic(base_size = 14, base_family = 'Helvetica')+
  facet_grid(species~temp, scales="free") +
  #geom_smooth ( aes(group=id, colour=Invasion, fill=Invasion),method="gam",formula= y ~ s(x, k = 23)) +
  scale_colour_manual(values=c("red","goldenrod1","red","goldenrod1"))+
  scale_fill_manual(values=c("red","goldenrod1","red","goldenrod1"))+
  theme(legend.position = "top")
salt4

#  test whether population sizes are sig different between invader and coloniser after 10 transfers (or AT 10th transfer)
mgr$transf<-as.numeric(mgr$transf)
t10<-subset(mgr,transf==10 & temp!="22/32 slow")
MMsal_10th <- lme(fixed = meanc.tr~  Invasion*temp*Environment, random = ~1|treatment, method = 'ML', data=t10)  
dd_10th<-dredge(MMsal_10th)
dd_10th

#write.csv(prettify(dd_10th),"modelselection_popsize.csv")

MMsal_10thfin <- lme(fixed = meanc.tr~temp*Environment, random = ~1|treatment, method = 'REML', data=t10)  
pretty_lm10th <- prettify(summary(MMsal_10thfin))
#write.csv(pretty_lm10th,"modeloutput_popsize.csv")


#this is super crude and just for rough test. Not used in final manuscript
mod2<-lme(meanc.tr~Invasion*temp*Environment, method= "ML" ,random=~1|treatment, data=t10) 
anova(mod2)
mod2.1<-lme(meanc.tr~Invasion*temp, method= "ML" ,random=~1|treatment, data=t10) 
anova(mod2.1)
anova(mod2,mod2.1)
mod3<-lme(meanc.tr~Invasion+temp+Environment, method= "ML" ,random=~1|treatment, data=t10) 
anova(mod2,mod3)


##### translocations and mixed models for salinity and temp adapt#### 
#  first use the file that Josianne had been using 
josisa<-read.csv("fitness.assay.data.csv")
str(josisa)
#drop temp FL
subjo<-subset(josisa, seltemp!="FL")
#drop the non focal species...i.e. chlamy, when ostreo is invading, and ostreo, when chlamy is invading , at assayTs that are not the selection Ts
subjo<-subjo[complete.cases(subjo), ]
# now we also drop the species being invaded as they will skew our ouput
subjo<-subset(subjo, howlongaway!='invaded')


summary(subjo)

subjo$seltemp<-as.character(subjo$seltemp)
subjo$assaytemp<-as.character(subjo$assaytemp)
#rename into something shorter so as not to have to the long names on the axis
subjo$howlongaway<-as.character(subjo$howlongaway)
subjo$Invasion<-as.character(subjo$Invasion)

subjo$howlongaway[subjo$howlongaway == 'alwayshome'] <- 'Home'
subjo$howlongaway[subjo$howlongaway == 'awayshort'] <- 'Short'
subjo$howlongaway[subjo$howlongaway == 'awaylong'] <- 'Long'
subjo$howlongaway[subjo$howlongaway == 'backhome'] <- 'Back'
subjo$howlongaway <- factor(subjo$howlongaway, levels = c('Home', 'Short', 'Long', 'Back'))

subjo$Invasion[subjo$Invasion == 'alone'] <- 'colonising'

salt7<-qplot( as.factor(howlongaway), rate.increase/5,data=subset(subjo, seltemp==assaytemp),ylab="growthrate µ/day",xlab="salinity", geom='boxplot',fill=Invasion, facets=seltemp~who)+theme_classic(base_size = 14, base_family = 'Helvetica') + scale_fill_manual(values=c('white','blue'))+theme(legend.position="top")+facet_wrap(seltemp~who, scales="free", ncol=2)
salt7

#we don't have an invasion into home scenario, but that makes logical sense :p. 
#now we display all data as fold change compared to , each species, home alone at 22ºC and because i am too tired to write an elegant function, we do like this:

meanchlam22<-mean(subset(subjo,Invasion=='colonising'&who=="Chlamydomonas"&seltemp=="22"&assaytemp=="22"&howlongaway=="Home")$rate.increase)
meanos22<-mean(subset(subjo,Invasion=='colonising'&who=="Ostreococcus"&seltemp=="22"&assaytemp=="22"&howlongaway=="Home")$rate.increase)


subjo$relative<-ifelse(subjo$who=="Chlamydomonas", subjo$rate.increase/meanchlam22, subjo$rate.increase/meanos22) # short term response magnitude, i.e. everything relative to 22 at 22  

#strength of response relative to 22 at 22. This IS the short-term response for everyone apart from where seltemp==assaytemp==22


# want to double check 26 short for chlamy invading
salt8<-qplot( as.factor(howlongaway), relative,data=subset(subjo, seltemp==assaytemp),ylab=" Fold change in growthrate compared to colonising control at 22ºC",xlab="Salinity regime", geom='boxplot',fill=Invasion, facets=who~seltemp)+theme_classic(base_size = 13, base_family = 'Helvetica') + scale_fill_manual(values=c('white','darkgrey'))+theme(legend.position="top")+facet_grid(seltemp~who, scales="free")+geom_segment(data=subset(subjo, seltemp==assaytemp), mapping=aes(x=0, y=1, xend=5,yend=1), linetype=3, colour='black', size=0.2)
salt8 #grey scale version -  use colour in final? 

ASSel<-subset(subjo, seltemp==assaytemp)


Chlamy_Inv<-subset(ASSel, who=="Chlamydomonas"  & Invasion=="invading")
meanChlamy_Inv<-mean(Chlamy_Inv$relative)
sdChlamy_Inv<-sd(Chlamy_Inv$relative)
Chlamy_Col<-subset(ASSel, who=="Chlamydomonas" & Invasion=="colonising")
meanChlamy_Col<-mean(Chlamy_Col$relative)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(ASSel, who=="Chlamydomonas" & seltemp=="22" & Invasion=="invading")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$relative)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$relative)
Chlamy_Col_22<-subset(ASSel, who=="Chlamydomonas" & seltemp=="22" & Invasion=="colonising")
meanChlamy_Col_22<-mean(Chlamy_Col_22$relative)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(ASSel, who=="Chlamydomonas" & seltemp=="26" & Invasion=="invading")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$relative)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$relative)
Chlamy_Col_26<-subset(ASSel, who=="Chlamydomonas" & seltemp=="26" & Invasion=="colonising")
meanChlamy_Col_26<-mean(Chlamy_Col_26$relative)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(ASSel, who=="Chlamydomonas" & seltemp=="32" & Invasion=="invading")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$relative)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$relative)
Chlamy_Col_32<-subset(ASSel, who=="Chlamydomonas" & seltemp=="32" & Invasion=="colonising")
meanChlamy_Col_32<-mean(Chlamy_Col_32$relative)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(ASSel, who=="Chlamydomonas" & seltemp=="FS" & Invasion=="invading")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$relative)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$relative)
Chlamy_Col_FS<-subset(ASSel, who=="Chlamydomonas" & seltemp=="FS" & Invasion=="colonising")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$relative)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(ASSel, who=="Ostreococcus"  & Invasion=="invading")
meanOstri_Inv<-mean(Ostri_Inv$relative)
sdOstri_Inv<-sd(Ostri_Inv$relative)
Ostri_Col<-subset(ASSel, who=="Ostreococcus" & Invasion=="colonising")
meanOstri_Col<-mean(Ostri_Col$relative)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(ASSel, who=="Ostreococcus" & seltemp=="22" & Invasion=="invading")
meanOstri_Inv_22<-mean(Ostri_Inv_22$relative)
sdOstri_Inv_22<-sd(Ostri_Inv_22$relative)
Ostri_Col_22<-subset(ASSel, who=="Ostreococcus" & seltemp=="22" & Invasion=="colonising")
meanOstri_Col_22<-mean(Ostri_Col_22$relative)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(ASSel, who=="Ostreococcus" & seltemp=="26" & Invasion=="invading")
meanOstri_Inv_26<-mean(Ostri_Inv_26$relative)
sdOstri_Inv_26<-sd(Ostri_Inv_26$relative)
Ostri_Col_26<-subset(ASSel, who=="Ostreococcus" & seltemp=="26" & Invasion=="colonising")
meanOstri_Col_26<-mean(Ostri_Col_26$relative)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(ASSel, who=="Ostreococcus" & seltemp=="32" & Invasion=="invading")
meanOstri_Inv_32<-mean(Ostri_Inv_32$relative)
sdOstri_Inv_32<-sd(Ostri_Inv_32$relative)
Ostri_Col_32<-subset(ASSel, who=="Ostreococcus" & seltemp=="32" & Invasion=="colonising")
meanOstri_Col_32<-mean(Ostri_Col_32$relative)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(ASSel, who=="Ostreococcus" & seltemp=="FS" & Invasion=="invading")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$relative)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$relative)
Ostri_Col_FS<-subset(ASSel, who=="Ostreococcus" & seltemp=="FS" & Invasion=="colonising")
meanOstri_Col_FS<-mean(Ostri_Col_FS$relative)
meanOstri_Inv_FS/meanOstri_Col_FS

## also locally adapted to temp? check that data has been compiled correctly!!! 
#we will want the order to be 22 -26-32-FS, and then we want to fit a line through them, too, so we need an id.column that has the invasion scenario AND the selection temperature. 
subjo<- within(subjo, id.plot <- as.character(factor(seltemp):factor(Invasion)))


subjo$assaytemp<-factor( subjo$assaytemp,  levels = c("22","26","FS","32"))

subjo$seltemp<-factor( subjo$seltemp,  levels = c("22","26","FS","32"))


###### make temp adapt figure for publication #####
salt8.1<-qplot( as.factor(assaytemp), relative,data=subset(subjo, howlongaway=="Long"),ylab="relative",xlab="Assay temperature (ºC)", geom='boxplot',fill=Invasion, facets=who~seltemp)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_fill_manual(values=c('orange','red'))+theme(legend.position="top")+ geom_smooth ( aes(group=id.plot, colour=Invasion),method="gam",se=FALSE,formula= y ~ s(x, k = 4)) +facet_grid(seltemp~who, scales="free_y")+scale_colour_manual(values=c("goldenrod1","darkred"))
salt8.1



HoLO<-subset(subjo, howlongaway=="Long")


Chlamy_Inv<-subset(HoLO, who=="Chlamydomonas"  & Invasion=="invading")
meanChlamy_Inv<-mean(Chlamy_Inv$relative)
sdChlamy_Inv<-sd(Chlamy_Inv$relative)
Chlamy_Col<-subset(HoLO, who=="Chlamydomonas" & Invasion=="colonising")
meanChlamy_Col<-mean(Chlamy_Col$relative)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(HoLO, who=="Chlamydomonas" & seltemp=="22" & Invasion=="invading")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$relative)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$relative)
Chlamy_Col_22<-subset(HoLO, who=="Chlamydomonas" & seltemp=="22" & Invasion=="colonising")
meanChlamy_Col_22<-mean(Chlamy_Col_22$relative)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(HoLO, who=="Chlamydomonas" & seltemp=="26" & Invasion=="invading")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$relative)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$relative)
Chlamy_Col_26<-subset(HoLO, who=="Chlamydomonas" & seltemp=="26" & Invasion=="colonising")
meanChlamy_Col_26<-mean(Chlamy_Col_26$relative)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(HoLO, who=="Chlamydomonas" & seltemp=="32" & Invasion=="invading")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$relative)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$relative)
Chlamy_Col_32<-subset(HoLO, who=="Chlamydomonas" & seltemp=="32" & Invasion=="colonising")
meanChlamy_Col_32<-mean(Chlamy_Col_32$relative)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(HoLO, who=="Chlamydomonas" & seltemp=="FS" & Invasion=="invading")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$relative)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$relative)
Chlamy_Col_FS<-subset(HoLO, who=="Chlamydomonas" & seltemp=="FS" & Invasion=="colonising")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$relative)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(HoLO, who=="Ostreococcus"  & Invasion=="invading")
meanOstri_Inv<-mean(Ostri_Inv$relative)
sdOstri_Inv<-sd(Ostri_Inv$relative)
Ostri_Col<-subset(HoLO, who=="Ostreococcus" & Invasion=="colonising")
meanOstri_Col<-mean(Ostri_Col$relative)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(HoLO, who=="Ostreococcus" & seltemp=="22" & Invasion=="invading")
meanOstri_Inv_22<-mean(Ostri_Inv_22$relative)
sdOstri_Inv_22<-sd(Ostri_Inv_22$relative)
Ostri_Col_22<-subset(HoLO, who=="Ostreococcus" & seltemp=="22" & Invasion=="colonising")
meanOstri_Col_22<-mean(Ostri_Col_22$relative)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(HoLO, who=="Ostreococcus" & seltemp=="26" & Invasion=="invading")
meanOstri_Inv_26<-mean(Ostri_Inv_26$relative)
sdOstri_Inv_26<-sd(Ostri_Inv_26$relative)
Ostri_Col_26<-subset(HoLO, who=="Ostreococcus" & seltemp=="26" & Invasion=="colonising")
meanOstri_Col_26<-mean(Ostri_Col_26$relative)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(HoLO, who=="Ostreococcus" & seltemp=="32" & Invasion=="invading")
meanOstri_Inv_32<-mean(Ostri_Inv_32$relative)
sdOstri_Inv_32<-sd(Ostri_Inv_32$relative)
Ostri_Col_32<-subset(HoLO, who=="Ostreococcus" & seltemp=="32" & Invasion=="colonising")
meanOstri_Col_32<-mean(Ostri_Col_32$relative)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(HoLO, who=="Ostreococcus" & seltemp=="FS" & Invasion=="invading")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$relative)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$relative)
Ostri_Col_FS<-subset(HoLO, who=="Ostreococcus" & seltemp=="FS" & Invasion=="colonising")
meanOstri_Col_FS<-mean(Ostri_Col_FS$relative)
meanOstri_Inv_FS/meanOstri_Col_FS



#####tiles - not used in accepted version of manuscript, but is in Biorxive version #####
####we want a geom_tiles plot, that has all environments 
library(RColorBrewer)
testtiledat<-subset(subjo, howlongaway=="Long" & Habitat=="away")
testtiledat$envqual<-ifelse(testtiledat$relative<1.0,"-", ifelse(testtiledat$relative=="1","0","+"))
testtiledat$envqual2<-ifelse(testtiledat$seltemp!="32"|testtiledat$seltemp!="32","good","bad")
testtiledat$envqual2<-factor( testtiledat$envqual2,  levels = c("good","bad"))


testtile<-ggplot(testtiledat, aes(envqual2, assaytemp))+geom_tile(aes(fill=as.factor(round(relative,1))))+facet_grid(who~Invasion)+  scale_fill_manual(values=c("#FFFFD9", "#FFFFD9","#FFFFD9","#EDF8B1", "#C7E9B4", "#7FCDBB","#7FCDBB", "#41B6C4", "#1D91C0", "#1D91C0", "#225EA8", "#253494","#253494", "#081D58", "#081D58", "#081D58"),na.value="grey90")+theme_classic() +geom_text( aes(label=envqual))
testtile

testtiledat$seltemp<-factor( testtiledat$seltemp,  levels = c("22","26","FS","32"))
testtiledat$assaytemp<-factor( testtiledat$assaytemp,  levels = c("22","26","FS","32"))


testtile01<-ggplot(testtiledat, aes(seltemp, assaytemp))+
  geom_tile(aes(fill=(round(relative,1))))+
  facet_grid(who~Invasion)+  
  theme_classic(base_size = 14, base_family = "Helvetica")+
  scale_fill_viridis_c()+
  geom_text( aes(label=envqual))
testtile01


testtile02<-ggplot(testtiledat, aes(envqual2, assaytemp))+
  geom_tile(aes(fill=round(relative,1)))+scale_fill_viridis_c()+
facet_grid(who~Invasion)+theme_classic() +geom_text( aes(label=envqual))
testtile02

head(subjo)
subjo$envqul<-ifelse(subjo$seltemp=="32", "bad", ifelse(subjo$seltemp=="22","neutral","good"))

subjo$envqul<-factor( subjo$envqul,  levels = c("good","neutral","bad"))

qplot(envqul,round(relative,1), data=subset(subjo, envqul!="neutral"), colour=round(relative,1))+geom_jitter()+scale_colour_viridis_c()


qplot(envqual2,round(relative,1), data=testtiledat, colour=round(relative,1))+geom_jitter()+scale_colour_viridis_c()

qplot(Invasion,round(relative,1), data=testtiledat, colour=round(relative,1))+geom_jitter()+scale_colour_viridis_c()


# try mixed model
MM_1 <- lme(fixed = round(relative,1) ~  Invasion * envqual2, random = ~ 1 |assaytemp, method = 'ML',testtiledat)
summary(MM_1)
dredge(MM_1) #additive
MM_1.1 <- lme(fixed = round(relative,1) ~  Invasion * envqual2, random = ~ 1 |assaytemp, method = 'REML',testtiledat)
summary(MM_1.1)
prettify(summary(MM_1.1))

MM_1.2 <- lme(fixed = round(relative,1) ~  Invasion , random = ~ 1 |assaytemp, method = 'REML',testtiledat)
MM_1.3 <- lme(fixed = round(relative,1) ~  envqual2 , random = ~ 1 |assaytemp, method = 'REML',testtiledat)

summary(MM_1.1)
# looks good!

#just for fun, also in lmer (error structure better?)
library(nloptr)
library(lme4)
testlmer<-lmer(round(relative,1)~ Invasion * envqual2+ (1|assaytemp), data=testtiledat)
testlmer01<-lmer(round(relative,1)~ Invasion + envqual2+ (1|assaytemp), data=testtiledat)

anova(testlmer,testlmer01) #no better when removing interaction 
summary(testlmer)

testlmerENV<-lmer(round(relative,1)~   envqual2+ (1|assaytemp), data=testtiledat)
anova(testlmer,testlmerENV) #Yes 
summary(testlmerENV)


testlmerINV<-lmer(round(relative,1)~   Invasion+ (1|assaytemp), data=testtiledat)
anova(testlmer,testlmerINV) #Yes 
summary(testlmerINV)

r.squaredGLMM(testlmer) 
r.squaredGLMM(testlmerENV) 
r.squaredGLMM(testlmerINV) 

### plotting
std_reg <- visreg(MM_1.1, 'Invasion', 'envqual2')
std_reg01 <- visreg(MM_1.2)
std_reg02 <- visreg(MM_1.3)

# visreg data
std_dat <- std_reg$res
std_dat01 <- std_reg01$res
std_dat02 <- std_reg02$res

std_fit <- std_reg$fit

#plot prediction vs real  - also not used in final accepted manuscript, but is in bioRxiv version 


# plot by Invasion
ggplot(std_dat, aes(y = visregRes, x= envqual2, colour = Invasion)) +    
  geom_boxplot()+
 theme_classic(base_size = 14, base_family = "Helvetica")+
  geom_line(data = std_fit, aes(y = visregFit, x= envqual2, colour = Invasion, group=Invasion))+
  scale_colour_manual(values = c('red', 'goldenrod1'))+
  theme(legend.position='top')



testtiledat$pred01<-std_dat01$visregRes

ggplot(testtiledat, aes(y = pred01, x= round(relative,1), colour = Invasion), xlab="predicted fold chang") +  
  theme_classic()+ 
  geom_jitter()+
  scale_colour_manual(values = c('red', 'goldenrod1'))+
  theme(legend.position='top')


# plot by goodbad

ggplot(std_dat, aes(y = visregRes, x= Invasion, colour = envqual2)) +     geom_boxplot()+
  theme_classic(base_size = 14, base_family = "Helvetica")+
  geom_line(data = std_fit, aes(y = visregFit, x= Invasion, colour = envqual2, group=envqual2))+
  scale_colour_manual(values = c("#41B6C4", 'darkgrey'))+
  theme(legend.position='top')


testtiledat$pred02<-std_dat02$visregRes

ggplot(testtiledat, aes(y = pred02, x= round(relative,1), colour = envqual2), xlab="predicted fold chang") +  
  theme_classic()+ 
  geom_jitter()+
  scale_colour_manual(values = c('blue', 'darkgrey'))+
  theme(legend.position='top')


#both!!!
ggplot(std_dat, aes(y = visregRes, x= envqual2, colour = visregRes, fill=Invasion)) +    facet_wrap(~Invasion, ncol=2)+  
  geom_boxplot()+
  geom_jitter(aes(alpha=0.85))+
  theme_classic()+ 
  scale_colour_viridis_c()+
  geom_smooth(aes(group=Invasion),method='lm',se=FALSE, size=0.51) + theme_classic(base_size = 14, base_family = "Helvetica")+
  theme(legend.position='top')



ggplot(std_dat, aes(y = visregRes, x= envqual2, fill = Invasion)) +    
  geom_boxplot()+
  scale_colour_viridis_c()+
  theme_classic(base_size = 14, base_family = "Helvetica")+
  theme_classic(base_size = 14, base_family = "Helvetica")+
  scale_fill_manual(values = c('red', 'goldenrod1'))+
  theme(legend.position='top')


testtiledat$pred<-std_dat$visregRes

ggplot(testtiledat, aes(y = pred, x= round(relative,1)), xlab="predicted fold chang") +  
  theme_classic()+ 
  facet_grid(envqual2~Invasion)+
  geom_point()+
  theme(legend.position='top')

#####more model selection ####

#now we build a model that tests how the relative difference depends on, how long they were away, the species, the monoculture/coculture, and the selection temperature. random effect is the biorpe nested in temperature 
subsubsub<-subset(subjo, howlongaway!="Home" &assaytemp==seltemp)
MMsal <- lme(fixed = relative~  Invasion*who*howlongaway*seltemp, random = ~1|plateloc/selsal, method = 'ML', data=subsubsub)  

MMsal <- lmer(relative~  Invasion*who*seltemp*howlongaway +(1|plateloc),subset(subjo, howlongaway!="Home"))  
summary(MMsal)
library(MuMIn)
dd1<-dredge(MMsal)
dd1
#write.csv(prettify(dd1),"modelselection_all_reci.csv")

MMfinsal <- lme(fixed = relative~  Invasion*who*howlongaway*seltemp, random = ~1|plateloc/selsal, method = 'REML', data=subsubsub) 
#write.csv(prettify(summary(MMfinsal)),"modelout_all_reci.csv")

#can use contrasts to see how much they changed (just col vs invaders across time)
ls02<-lsmeans(MMfinsal, pairwise ~ list(Invasion,howlongaway), adjust = "tukey") 
ls02

summary(MMfinsal)
coef(MMfinsal)
plot(MMfinsal) #okay 
hist(resid(MMfinsal))   #rather nice
data.frame(fixed.effects(MMfinsal))
(random.effects(MMfinsal))


subsubsubsub<-subset(subjo,howlongaway=="Long")

MMtemp <- lme(fixed = relative~ as.factor(seltemp)*Invasion*as.factor(who)*as.factor(assaytemp), random = ~1|plateloc/seltemp, method = 'ML',subsubsubsub)  

dd<-dredge(MMtemp)
dd
#write.csv(prettify(dd),"modelsel_all_temp.csv")
#write.csv(data.frame(dd),'tempdredge.csv')
# best model is model of interaction assay:selectiontemp:invasion, but not who
#refit with REML

MMfintemp <- lme(fixed = relative~ seltemp*Invasion*assaytemp, random = ~1|plateloc/seltemp, method = 'REML',subsubsubsub)  
#write.csv(prettify(summary(MMfintemp)),"modelout_all_temp.csv")

anova (MMfintemp)
sss<-summary(MMfintemp)
coef(MMfintemp)
plot(MMfintemp) #okay 
hist(resid(MMfintemp))   #rather nice
data.frame(fixed.effects(MMfintemp))
(random.effects(MMfintemp))

#can use contrasts to see how much they changed (just col vs invaders across time)
ls03<-lsmeans(MMfintemp, pairwise ~ list(Invasion,seltemp,assaytemp), adjust = "tukey") 
ls03

### decomposed samples ####

# short decomposed 
shortde<-read.csv("decomposed short.csv")
head(shortde)
shortde<-shortde[,c(1:7)]

#make temps right order
shortde$temp<-factor( shortde$temp,  levels = c("22","26","22/32","32"))

library(car)
#plot with 1:1 line 
shortie<-qplot(growthmono, growthafterdec,data=shortde, facets=.~species, shape=was, colour=temp, xlab="Growth in mono-culture", ylab="Growth after decomposition")+theme_classic(base_size = 14, base_family = 'Helvetica')+theme(legend.position='top')+  geom_abline(slope = 1, intercept = 0, lty = 2) +scale_colour_manual("Selection Temperature ºC",values=c("darkgreen","blue","purple","red"))

shortie # are almost perfectly on 1:1 line 

r.x <- lm(growthmono ~ growthafterdec, data=shortde)
r1 <- lm(growthmono ~1 + offset(growthafterdec), data=shortde)
anova(r.x, r1)

t.test(shortde$growthmono, shortde$growthafterdec)
#long decomposed 

longde<-read.csv("decomposed samples growth.csv")
head(longde)
longde<-longde[,c(1:4)]

longde$Temp<-factor( longde$Temp,  levels = c("22","26","22/32","32"))
longde$Was<-factor( longde$Was,  levels = c("Resident","Invader"))


#graph needs to be about just growth rate compared to growth in co-culture  - make main manuscript figure? 

longing<-qplot(Was, growthchangeFULL,data=longde, facets=Temp~Species,fill=Was, xlab="Previous role in co-culture", ylab="Growth rate alone after decomposition/growth rate in co-culture",geom="boxplot")+theme_classic(base_size = 14, base_family = 'Helvetica')+
  theme(legend.position='top')+scale_fill_manual(values=c("forestgreen","red"))+
  geom_segment(data=longde, mapping=aes(x=0.5, y=1, xend=2.5,yend=1), linetype=3, colour='black', size=0.2)
longing


longde <- ddply(longde, .(Temp, Was, Species), mutate, id.new =seq(1,length(Temp),1))
longde<-within(longde, id.newer<-as.character(factor(Species):factor(id.new)))
head(longde)


modslong <- lme(growthchangeFULL~ Temp*Was*Species, random = ~1|id.newer/Temp, method = 'ML',longde)
dmd<-dredge(modslong)
dmd


#####with NILE red #####
growthNILE<-read.csv("growthNILE.csv")

growthNILE$Temp<-factor( growthNILE$Temp,  levels = c("22","26","22/32","32"))
growthNILE$Was<-factor( growthNILE$Was,  levels = c("Resident","Invader"))

saltNILE<-qplot(growthchangeFULL, NILE_FL, data=growthNILE,xlab="Growth rate alone after decomposition/growth rate in co-culture",ylab="Nile Red Fluorescence",colour=Was, facets=Temp~Species)+
  theme_classic(base_size = 14, base_family = 'Helvetica') + 
  facet_wrap(Temp~Species, scales="free_x", ncol=2)+
  scale_colour_manual("was former",values=c('forestgreen','red'))+theme(legend.position="top")+
  geom_segment(data=longde, mapping=aes(x=0.5, y=10, xend=1.5,yend=10), linetype=3, colour='black', size=0.2)
saltNILE

#### forest plot NILE only ####
head(growthNILE)


Chlamy_Inv<-subset(growthNILE, Species=="Chlamydmonas"  & Was=="Invader")
meanChlamy_Inv<-mean(Chlamy_Inv$NILE_FL)
sdChlamy_Inv<-sd(Chlamy_Inv$NILE_FL)

Chlamy_Col<-subset(growthNILE, Species=="Chlamydmonas" & Was=="Resident")
meanChlamy_Col<-mean(Chlamy_Col$NILE_FL)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="22" & Was=="Invader")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$NILE_FL)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$NILE_FL)
Chlamy_Col_22<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="22" & Was=="Resident")
meanChlamy_Col_22<-mean(Chlamy_Col_22$NILE_FL)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="26" & Was=="Invader")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$NILE_FL)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$NILE_FL)
Chlamy_Col_26<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="26" & Was=="Resident")
meanChlamy_Col_26<-mean(Chlamy_Col_26$NILE_FL)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="32" & Was=="Invader")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$NILE_FL)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$NILE_FL)
Chlamy_Col_32<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="32" & Was=="Resident")
meanChlamy_Col_32<-mean(Chlamy_Col_32$NILE_FL)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="22/32" & Was=="Invader")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$NILE_FL)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$NILE_FL)
Chlamy_Col_FS<-subset(growthNILE, Species=="Chlamydmonas" & Temp=="22/32" & Was=="Resident")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$NILE_FL)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(growthNILE, Species=="Ostreococcus"  & Was=="Invader")
meanOstri_Inv<-mean(Ostri_Inv$NILE_FL)
sdOstri_Inv<-sd(Ostri_Inv$NILE_FL)
Ostri_Col<-subset(growthNILE, Species=="Ostreococcus" & Was=="Resident")
meanOstri_Col<-mean(Ostri_Col$NILE_FL)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(growthNILE, Species=="Ostreococcus" & Temp=="22" & Was=="Invader")
meanOstri_Inv_22<-mean(Ostri_Inv_22$NILE_FL)
sdOstri_Inv_22<-sd(Ostri_Inv_22$NILE_FL)
Ostri_Col_22<-subset(growthNILE, Species=="Ostreococcus" & Temp=="22" & Was=="Resident")
meanOstri_Col_22<-mean(Ostri_Col_22$NILE_FL)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(growthNILE, Species=="Ostreococcus" & Temp=="26" & Was=="Invader")
meanOstri_Inv_26<-mean(Ostri_Inv_26$NILE_FL)
sdOstri_Inv_26<-sd(Ostri_Inv_26$NILE_FL)
Ostri_Col_26<-subset(growthNILE, Species=="Ostreococcus" & Temp=="26" & Was=="Resident")
meanOstri_Col_26<-mean(Ostri_Col_26$NILE_FL)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(growthNILE, Species=="Ostreococcus" & Temp=="32" & Was=="Invader")
meanOstri_Inv_32<-mean(Ostri_Inv_32$NILE_FL)
sdOstri_Inv_32<-sd(Ostri_Inv_32$NILE_FL)
Ostri_Col_32<-subset(growthNILE, Species=="Ostreococcus" & Temp=="32" & Was=="Resident")
meanOstri_Col_32<-mean(Ostri_Col_32$NILE_FL)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(growthNILE, Species=="Ostreococcus" & Temp=="22/32" & Was=="Invader")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$NILE_FL)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$NILE_FL)
Ostri_Col_FS<-subset(growthNILE, Species=="Ostreococcus" & Temp=="22/32" & Was=="Resident")
meanOstri_Col_FS<-mean(Ostri_Col_FS$NILE_FL)
meanOstri_Inv_FS/meanOstri_Col_FS


#best model does not include temperature and has spc and Was, and the two of them in interaction 
#write.csv(dmd, "dredge long composed.csv")

modslonger <- lme(growthchangeFULL~ Was*Species, random = ~1|id.newer/Temp, method = 'ML',longde)
summary(modslonger)


### decomposed samples NP ####
decNP<-read.csv("NP decomposed.csv")
head(decNP)
decNP$Temp<-factor( decNP$Temp,  levels = c("22","26","22/32","32"))
decNP$was<-factor( decNP$was,  levels = c("colonising","mixed","invader","resident"))

qplot(was, NP, data=decNP, facets=Temp~Spec, geom="boxplot", fill=Salinity)+theme_classic(base_size = 14, base_family = 'Helvetica')+
  scale_fill_manual(values=c("skyblue","cornflowerblue"))+
  theme(legend.position = "top")+  
  labs(x=expression(Sample~origin), y=expression(NP~gC~gC^{-1}~h^{-1}))+
  facet_grid(Temp~Spec)

#make a ddply to make means per "was"
means1<-ddply(subset(decNP, was!="colonising"),c("was") , function(df) return(c(meanNP=mean(df$NP), sdNP=sd(df$NP)/sqrt(8))))
decNPsub<-subset(decNP, was!="colonising")
#addedfun<- ddply(decNPsub, .(was, Temp), summarise, av =  mean(NP, na.rm=T), semu = sd(NP, na.rm=T)/sqrt(length(NP)))
#addedfun<- ddply(decNPsub, .(was, Temp), mutate, av =  mean(NP, na.rm=T), semu = sd(NP, na.rm=T)/sqrt(length(NP)))
#addedfun<- ddply(decNPsub, .(was, Temp), av =  mean(NP, na.rm=T), semu = sd(NP, na.rm=T)/sqrt(length(NP)))
means1


#now just for the mixed and the one that was the invader and the resident 

qplot(was, NP, data=subset(decNP, was=="mixed"|was=="invader"|was=="resident"), facets=Temp~focSpec, geom="boxplot", fill=was)+theme_classic(base_size = 14, base_family = 'Helvetica')+scale_fill_manual(values=c("skyblue","darkblue","cornflowerblue"))+theme(legend.position = "top")+  labs(x=expression(Sample~origin), y=expression(NP~gC~gC^{-1}~h^{-1}))+facet_grid(Temp~Spec)

qplot(Sample, NP, data=decNP, facets=Temp~focSpec, geom="boxplot", fill=was)+theme_classic(base_size = 14, base_family = 'Helvetica')+theme(legend.position = "top")+  labs(x=expression(Sample), y=expression(NP~gC~gC^{-1}~h^{-1}))+facet_grid(Temp~Spec,scales="free_x")

mod_NP <- lme(NP~ was*Spec*Temp, random = ~1|idboth, method = 'ML',decNP)
dmd1<-dredge(mod_NP)
dmd1 #everything, but not in interaction 
mod_NP2 <- lme(NP~ was+Spec+Temp, random = ~1|idboth, method = 'REML',decNP)
summary(mod_NP2)
#write.csv(dmd1,"modelNPdec.csv")

#### NP forest plot #####
head(decNP)
subNP<-subset(decNP,was=="invader"|was=="colonising")

head(subNP)

Chlamy_Inv<-subset(subNP, focspec=="Chlamydomonas"  & was=="invader")
meanChlamy_Inv<-mean(Chlamy_Inv$NP)
sdChlamy_Inv<-sd(Chlamy_Inv$NP)

Chlamy_Col<-subset(subNP, focspec=="Chlamydomonas" & was=="colonising")
meanChlamy_Col<-mean(Chlamy_Col$NP)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(subNP, focspec=="Chlamydomonas" & Temp=="22" & was=="invader")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$NP)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$NP)
Chlamy_Col_22<-subset(subNP, focspec=="Chlamydomonas" & Temp=="22" & was=="colonising")
meanChlamy_Col_22<-mean(Chlamy_Col_22$NP)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(subNP, focspec=="Chlamydomonas" & Temp=="26" & was=="invader")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$NP)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$NP)
Chlamy_Col_26<-subset(subNP, focspec=="Chlamydomonas" & Temp=="26" & was=="colonising")
meanChlamy_Col_26<-mean(Chlamy_Col_26$NP)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(subNP, focspec=="Chlamydomonas" & Temp=="32" & was=="invader")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$NP)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$NP)
Chlamy_Col_32<-subset(subNP, focspec=="Chlamydomonas" & Temp=="32" & was=="colonising")
meanChlamy_Col_32<-mean(Chlamy_Col_32$NP)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(subNP, focspec=="Chlamydomonas" & Temp=="22/32" & was=="invader")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$NP)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$NP)
Chlamy_Col_FS<-subset(subNP, focspec=="Chlamydomonas" & Temp=="22/32" & was=="colonising")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$NP)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(subNP, focspec=="Ostreococcus"  & was=="invader")
meanOstri_Inv<-mean(Ostri_Inv$NP)
sdOstri_Inv<-sd(Ostri_Inv$NP)
Ostri_Col<-subset(subNP, focspec=="Ostreococcus" & was=="colonising")
meanOstri_Col<-mean(Ostri_Col$NP)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(subNP, focspec=="Ostreococcus" & Temp=="22" & was=="invader")
meanOstri_Inv_22<-mean(Ostri_Inv_22$NP)
sdOstri_Inv_22<-sd(Ostri_Inv_22$NP)
Ostri_Col_22<-subset(subNP, focspec=="Ostreococcus" & Temp=="22" & was=="colonising")
meanOstri_Col_22<-mean(Ostri_Col_22$NP)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(subNP, focspec=="Ostreococcus" & Temp=="26" & was=="invader")
meanOstri_Inv_26<-mean(Ostri_Inv_26$NP)
sdOstri_Inv_26<-sd(Ostri_Inv_26$NP)
Ostri_Col_26<-subset(subNP, focspec=="Ostreococcus" & Temp=="26" & was=="colonising")
meanOstri_Col_26<-mean(Ostri_Col_26$NP)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(subNP, focspec=="Ostreococcus" & Temp=="32" & was=="invader")
meanOstri_Inv_32<-mean(Ostri_Inv_32$NP)
sdOstri_Inv_32<-sd(Ostri_Inv_32$NP)
Ostri_Col_32<-subset(subNP, focspec=="Ostreococcus" & Temp=="32" & was=="colonising")
meanOstri_Col_32<-mean(Ostri_Col_32$NP)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(subNP, focspec=="Ostreococcus" & Temp=="22/32" & was=="invader")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$NP)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$NP)
Ostri_Col_FS<-subset(subNP, focspec=="Ostreococcus" & Temp=="22/32" & was=="colonising")
meanOstri_Col_FS<-mean(Ostri_Col_FS$NP)
meanOstri_Inv_FS/meanOstri_Col_FS





### clone attack #### 
#question is whether new clones behave differently from each other (and from the FS samples ) - not used for manuscript in the end. kept in case reader/reviewer questions. 

clones<-read.csv("clones_SAL only.csv")
head(clones) 
str(clones)
#make plot and analysis as before 
clones$SelTemp<-factor( clones$SelTemp,  levels = c("22","26","22/32","32"))
clones$SelSal<-factor( clones$SelSal,  levels = c("Short","Long","Back"))


salt_c<-qplot( as.factor(SelSal), relgrowth,data=clones,ylab=" Fold change in growthrate compared to colonising control at 22ºC",xlab="Salinity regime", geom='jitter',colour=Invasion, facets=Species~SelTemp, alpha=.75)+theme_classic(base_size = 13, base_family = 'Helvetica') +theme(legend.position="top")+facet_grid(SelTemp~Species, scales="free")+geom_segment(data=clones, mapping=aes(x=0, y=1, xend=5,yend=1), linetype=3, colour='black', size=0.2)+ scale_colour_manual(values=c('goldenrod1','red'))

salt_c

#now analyse 
 clones<-within(clones,id.new<-as.character(factor(Clone):factor(SelTemp):factor(Species)))

MM_clones <- lme(relgrowth~  Invasion*Species*SelSal*SelTemp, random = ~1|Clone, method = 'ML', data=clones)  

summary(MM_clones)
library(MuMIn)
dd1<-dredge(MM_clones)
dd1
# has the best model with full interaction 
#write.csv(data.frame(dd1),'saldredge.csv')

MM_finclones <- lme(relgrowth~  Invasion*Species*SelSal*SelTemp, random = ~1|Clone, method = 'REML', data=clones)  

summary(MM_finclones)
coef(MM_finclones)
plot(MM_finclones) #okay 
hist(resid(MM_finclones))   #rather nice
data.frame(fixed.effects(MM_finclones))
(random.effects(MM_finclones))

#### short versus long #####

#can I now cast this to plot Short against Long? 
sal2<-reshape::cast(subset(subjo, seltemp==assaytemp), fulltreat*seltemp~howlongaway,fun=list(mean,sd),value="relative") 
# yeah but is not great because we have differnt numbers of mortalities

### or read in file with short and long

longshort<-read.csv("shoRtlong.csv") #also kept in mainly for reviewer questions, is not in final document
head(longshort)  

longshort<-within(longshort, id.plot<-as.character(factor(who):factor(seltemp)))

#make the treatment the same order as in the others, too. 

longshort$seltemp<-factor(longshort$seltemp,c("22","26","22/32","32"))

plot9.4<-qplot(SHORTresp,LONGresp, data=longshort, colour=Invasion,shape=Invasion,ylab="Magnitude of long-term response", xlab="Magnitude of short-term response")+theme_classic(base_size = 14, base_family = 'Helvetica')+scale_colour_manual(values=c("goldenrod1","red"))+geom_smooth(aes(group=Invasion, colour=Invasion),size=0.4,method='lm',se=FALSE, fullrange=FALSE)+facet_grid(seltemp~who,scales="free")+theme(legend.position='top')
plot9.4

plot9.5<-qplot(growthshort/5,growthlong/5 , data=longshort, colour=Invasion,shape=Invasion,ylab="growth rate in long-term response", xlab="Growth rate in short-term response")+theme_classic(base_size = 14, base_family = 'Helvetica')+scale_colour_manual(values=c("goldenrod1","red"))+geom_smooth(aes(group=Invasion, colour=Invasion),size=0.4,method='lm',se=FALSE, fullrange=FALSE)+facet_grid(seltemp~who)+theme(legend.position='top')
plot9.5


#slopes function
slopefn<-function(d) { 
  m<-lm (LONGresp~SHORTresp,data=d) #can force through origin by setting offset to 0 but we know it is not going through unity in the first place from looking at the data in plot9
  sum.<-summary(m) 
  r2<-sum.$r.squared
  intercept<-m$coefficients[1]
  slope<-m$coefficients[2]
  pval<-sum.$coefficients[8]
  output<-data.frame(slope,intercept,r2,pval)
  output}
slopedata3<-ddply(longshort,.(Invasion, who, seltemp),slopefn)
slopedata3

#add another id column
longshort <- ddply(longshort, .(Invasion, who, seltemp), mutate, id.new =seq(1,length(id.plot),1))
head(longshort)
longshort<-within(longshort, idsupernew<- as.character(factor(id.new):factor(id.plot)))


#or try visreg
MM_1 <- lme(fixed = LONGresp ~  SHORTresp * seltemp*Invasion*who, random = ~ 1 |idsupernew, method = 'ML', longshort)
predict(MM_1, type = "response") # for the random effects 

summary(MM_1) # ok
dm<-dredge(MM_1) #yup , interact of everything in pairs of three, apart frmo Inv:seltemp:short , but we need to average several models to get delta AIC >2 ...first five models 

#write.csv(dm,"outputlongshort.csv") 
summary(model.avg(dm, subset = delta < 2)) 




### plotting
std_reg <- visreg(MM_1.1, 'SHORTresp', 'seltemp')
visreg(MM_1.1, "SHORTresp", by="seltemp", overlay=TRUE)

# visreg data
std_dat <- std_reg$res
std_fit <- std_reg$fit

# plot
ggplot(std_dat1, aes(y = visregRes, x= SHORTresp, colour = as.factor())) +  # shuold add components? 
  geom_point(size = 2, alpha=0.765) +scale_colour_manual(values=c("darkgreen","blue","red","darkgrey","purple"))+
  geom_line(data = std_fit1, aes(y = visregFit, x= lnNet, colour = as.factor(growthT)))+ #should add components?? 
  labs(y=expression(ln~(Growth~rate~µ~(day^{-1}))), x=expression(ln(NP)+ln(CUE)+ln(P:N))) +
  theme(legend.position='top')





##### size plot again , but no trajectory #### 
#now we make a size plot that is basically the same as the new main plot, or not, because we look at size in the selectione environment at the end, not across assays. 
sz<-read.csv("biomass.cellsize.over.time.csv")
#keep only complete cases 
sz<-sz[complete.cases(sz), ]

subsz<-subset(sz, transf==25&size<20) #there are some outliers that cannot be true. 100µm?? likely a clump of sorts

subsz$Invasion<-as.character(subsz$Invasion) 


subsz$temp <- factor(subsz$temp, levels = c('22', '26', '22/32 rapid','32'))
subsz$Invasion<-factor(subsz$Invasion, levels=c('mono-culture','co-culture'))

####use this for size SI
salt11<-qplot( Habitat, size,data=subset(subsz, temp!="22/32 slow"),ylab="Cell size (µm diameter)",xlab="Salinity regime", geom='boxplot',fill=Invasion, facets=species~temp)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_fill_manual(values=c('goldenrod1','red'))+theme(legend.position="top")+facet_grid(species~temp, scales="free")#+geom_segment(data=subset(subjo, seltemp==assaytemp), mapping=aes(x=0, y=1, xend=5,yend=1), linetype=3, colour='black', size=0.2)
salt11 


### size for forest plot ####

alldat<-subset(subsz, temp!="22/32 slow" & Habitat=="away")

Chlamy_Inv<-subset(growthNILE, species=="Chlamydomonas"  & Invasion=="co-culture")
meanChlamy_Inv<-mean(Chlamy_Inv$size)
sdChlamy_Inv<-sd(Chlamy_Inv$size)

Chlamy_Col<-subset(growthNILE, species=="Chlamydomonas" & Invasion=="mono-culture")
meanChlamy_Col<-mean(Chlamy_Col$size)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(subsi, species=="Chlamydomonas" & temp=="22" & Invasion=="co-culture")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$size)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$size)
Chlamy_Col_22<-subset(subsi, species=="Chlamydomonas" & temp=="22" & Invasion=="mono-culture")
meanChlamy_Col_22<-mean(Chlamy_Col_22$size)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(subsi, species=="Chlamydomonas" & temp=="26" & Invasion=="co-culture")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$size)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$size)
Chlamy_Col_26<-subset(subsi, species=="Chlamydomonas" & temp=="26" & Invasion=="mono-culture")
meanChlamy_Col_26<-mean(Chlamy_Col_26$size)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(subsi, species=="Chlamydomonas" & temp=="32" & Invasion=="co-culture")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$size)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$size)
Chlamy_Col_32<-subset(subsi, species=="Chlamydomonas" & temp=="32" & Invasion=="mono-culture")
meanChlamy_Col_32<-mean(Chlamy_Col_32$size)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(subsi, species=="Chlamydomonas" & temp=="22/32 rapid" & Invasion=="co-culture")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$size)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$size)
Chlamy_Col_FS<-subset(subsi, species=="Chlamydomonas" & temp=="22/32 rapid" & Invasion=="mono-culture")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$size)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(subsi, species=="Ostreococcus"  & Invasion=="co-culture")
meanOstri_Inv<-mean(Ostri_Inv$size)
sdOstri_Inv<-sd(Ostri_Inv$size)
Ostri_Col<-subset(subsi, species=="Ostreococcus" & Invasion=="mono-culture")
meanOstri_Col<-mean(Ostri_Col$size)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(subsi, species=="Ostreococcus" & temp=="22" & Invasion=="co-culture")
meanOstri_Inv_22<-mean(Ostri_Inv_22$size)
sdOstri_Inv_22<-sd(Ostri_Inv_22$size)
Ostri_Col_22<-subset(subsi, species=="Ostreococcus" & temp=="22" & Invasion=="mono-culture")
meanOstri_Col_22<-mean(Ostri_Col_22$size)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(subsi, species=="Ostreococcus" & temp=="26" & Invasion=="co-culture")
meanOstri_Inv_26<-mean(Ostri_Inv_26$size)
sdOstri_Inv_26<-sd(Ostri_Inv_26$size)
Ostri_Col_26<-subset(subsi, species=="Ostreococcus" & temp=="26" & Invasion=="mono-culture")
meanOstri_Col_26<-mean(Ostri_Col_26$size)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(subsi, species=="Ostreococcus" & temp=="32" & Invasion=="co-culture")
meanOstri_Inv_32<-mean(Ostri_Inv_32$size)
sdOstri_Inv_32<-sd(Ostri_Inv_32$size)
Ostri_Col_32<-subset(subsi, species=="Ostreococcus" & temp=="32" & Invasion=="mono-culture")
meanOstri_Col_32<-mean(Ostri_Col_32$size)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(subsi, species=="Ostreococcus" & temp=="22/32 rapid" & Invasion=="co-culture")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$size)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$size)
Ostri_Col_FS<-subset(subsi, species=="Ostreococcus" & temp=="22/32 rapid" & Invasion=="mono-culture")
meanOstri_Col_FS<-mean(Ostri_Col_FS$size)
meanOstri_Inv_FS/meanOstri_Col_FS



#### size vs growth vs ROS ####
#we check whether growth changes as a function of fitness for all monos and cos 


alldat<-read.csv("sizegrowthROS_datasets_BOTH.csv")

alldat$temp<-factor( alldat$temp,  levels = c("22","26","FS","32"))
alldat$Size<-as.numeric(alldat$Size)

salt13<-qplot(Size, Growth, data=alldat,xlab="Cell size (µm diameter)",ylab="Growth rate µ day-1",colour=scene, facets=temp~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_colour_manual(values=c('goldenrod1','red'))+theme(legend.position="top")+facet_wrap(who~temp, scales="free",ncol=4)+geom_smooth(aes(group=scene, colour=scene), method="lm", se=FALSE, linetype=3, size=0.2)
salt13


salt14<-qplot(Size, ROS_intracell_vol, data=alldat,xlab="Cell size (µm diameter)",ylab="ROS per cell biovolume (nM)",colour=scene, facets=temp~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_colour_manual(values=c('goldenrod1','red'))+theme(legend.position="top")+facet_wrap(who~temp, scales="free",ncol=4)+geom_smooth(aes(group=scene, colour=scene), method="lm", se=FALSE, linetype=3, size=0.2)
salt14

salt15<-qplot(Size, ROS_LD50, data=alldat,xlab="Cell size (µm diameter)",ylab="ROS LD50 (µM)",colour=scene, facets=temp~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_colour_manual(values=c('goldenrod1','red'))+theme(legend.position="top")+facet_wrap(who~temp, scales="free",ncol=4)+geom_smooth(aes(group=scene, colour=scene), method="lm", se=FALSE, linetype=3, size=0.2)
salt15

salt15<-qplot(Size, Growth, data=alldat,xlab="Cell size (µm diameter)",ylab="Growth Rate",colour=scene, alpha=ROS_intracell_vol, facets=temp~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + scale_colour_manual(values=c('goldenrod1','red'))+theme(legend.position="top")+facet_wrap(who~temp, scales="free",ncol=4)+geom_smooth(aes(group=scene, colour=scene), method="lm", se=FALSE, linetype=3, size=0.2)
salt15


salt16<-qplot(temp, ROS_intracell_vol, data=alldat,xlab="temperature",ylab="ROS_intracell_vol",fill=scene, geom="boxplot", facets=.~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + facet_grid(who~.)+ scale_fill_manual(values=c('goldenrod1','red'))+theme(legend.position="top")
salt16

salt17<-qplot(temp, ROS_LD50, data=alldat,xlab="temperature",ylab="ROS_LD50",fill=scene, geom="boxplot", facets=.~who)+theme_classic(base_size = 12, base_family = 'Helvetica') + facet_grid(who~.)+ scale_fill_manual(values=c('goldenrod1','red'))+theme(legend.position="top")
salt17


#### ROS Cell intra forest plot ####
head(alldat)

Chlamy_Inv<-subset(alldat, who=="Chlamydomonas"  & scene=="invading")
meanChlamy_Inv<-mean(Chlamy_Inv$ROS_intracell_vol)
sdChlamy_Inv<-sd(Chlamy_Inv$ROS_intracell_vol)

Chlamy_Col<-subset(alldat, who=="Chlamydomonas" & scene=="colonising")
meanChlamy_Col<-mean(Chlamy_Col$ROS_intracell_vol)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(alldat, who=="Chlamydomonas" & temp=="22" & scene=="invading")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$ROS_intracell_vol)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$ROS_intracell_vol)
Chlamy_Col_22<-subset(alldat, who=="Chlamydomonas" & temp=="22" & scene=="colonising")
meanChlamy_Col_22<-mean(Chlamy_Col_22$ROS_intracell_vol)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(alldat, who=="Chlamydomonas" & temp=="26" & scene=="invading")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$ROS_intracell_vol)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$ROS_intracell_vol)
Chlamy_Col_26<-subset(alldat, who=="Chlamydomonas" & temp=="26" & scene=="colonising")
meanChlamy_Col_26<-mean(Chlamy_Col_26$ROS_intracell_vol)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(alldat, who=="Chlamydomonas" & temp=="32" & scene=="invading")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$ROS_intracell_vol)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$ROS_intracell_vol)
Chlamy_Col_32<-subset(alldat, who=="Chlamydomonas" & temp=="32" & scene=="colonising")
meanChlamy_Col_32<-mean(Chlamy_Col_32$ROS_intracell_vol)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(alldat, who=="Chlamydomonas" & temp=="FS" & scene=="invading")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$ROS_intracell_vol)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$ROS_intracell_vol)
Chlamy_Col_FS<-subset(alldat, who=="Chlamydomonas" & temp=="FS" & scene=="colonising")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$ROS_intracell_vol)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(alldat, who=="Ostreococccus"  & scene=="invading")
meanOstri_Inv<-mean(Ostri_Inv$ROS_intracell_vol)
sdOstri_Inv<-sd(Ostri_Inv$ROS_intracell_vol)
Ostri_Col<-subset(alldat, who=="Ostreococccus" & scene=="colonising")
meanOstri_Col<-mean(Ostri_Col$ROS_intracell_vol)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(alldat, who=="Ostreococccus" & temp=="22" & scene=="invading")
meanOstri_Inv_22<-mean(Ostri_Inv_22$ROS_intracell_vol)
sdOstri_Inv_22<-sd(Ostri_Inv_22$ROS_intracell_vol)
Ostri_Col_22<-subset(alldat, who=="Ostreococccus" & temp=="22" & scene=="colonising")
meanOstri_Col_22<-mean(Ostri_Col_22$ROS_intracell_vol)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(alldat, who=="Ostreococccus" & temp=="26" & scene=="invading")
meanOstri_Inv_26<-mean(Ostri_Inv_26$ROS_intracell_vol)
sdOstri_Inv_26<-sd(Ostri_Inv_26$ROS_intracell_vol)
Ostri_Col_26<-subset(alldat, who=="Ostreococccus" & temp=="26" & scene=="colonising")
meanOstri_Col_26<-mean(Ostri_Col_26$ROS_intracell_vol)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(alldat, who=="Ostreococccus" & temp=="32" & scene=="invading")
meanOstri_Inv_32<-mean(Ostri_Inv_32$ROS_intracell_vol)
sdOstri_Inv_32<-sd(Ostri_Inv_32$ROS_intracell_vol)
Ostri_Col_32<-subset(alldat, who=="Ostreococccus" & temp=="32" & scene=="colonising")
meanOstri_Col_32<-mean(Ostri_Col_32$ROS_intracell_vol)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(alldat, who=="Ostreococccus" & temp=="FS" & scene=="invading")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$ROS_intracell_vol)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$ROS_intracell_vol)
Ostri_Col_FS<-subset(alldat, who=="Ostreococccus" & temp=="FS" & scene=="colonising")
meanOstri_Col_FS<-mean(Ostri_Col_FS$ROS_intracell_vol)
meanOstri_Inv_FS/meanOstri_Col_FS

#### ROS detox frest plot ####


Chlamy_Inv<-subset(alldat, who=="Chlamydomonas"  & scene=="invading")
meanChlamy_Inv<-mean(Chlamy_Inv$ROS_LD50)
sdChlamy_Inv<-sd(Chlamy_Inv$ROS_LD50)

Chlamy_Col<-subset(alldat, who=="Chlamydomonas" & scene=="colonising")
meanChlamy_Col<-mean(Chlamy_Col$ROS_LD50)
meanChlamy_Inv/meanChlamy_Col

Chlamy_Inv_22<-subset(alldat, who=="Chlamydomonas" & temp=="22" & scene=="invading")
meanChlamy_Inv_22<-mean(Chlamy_Inv_22$ROS_LD50)
sdChlamy_Inv_22<-sd(Chlamy_Inv_22$ROS_LD50)
Chlamy_Col_22<-subset(alldat, who=="Chlamydomonas" & temp=="22" & scene=="colonising")
meanChlamy_Col_22<-mean(Chlamy_Col_22$ROS_LD50)
meanChlamy_Inv_22/meanChlamy_Col_22

Chlamy_Inv_26<-subset(alldat, who=="Chlamydomonas" & temp=="26" & scene=="invading")
meanChlamy_Inv_26<-mean(Chlamy_Inv_26$ROS_LD50)
sdChlamy_Inv_26<-sd(Chlamy_Inv_26$ROS_LD50)
Chlamy_Col_26<-subset(alldat, who=="Chlamydomonas" & temp=="26" & scene=="colonising")
meanChlamy_Col_26<-mean(Chlamy_Col_26$ROS_LD50)
meanChlamy_Inv_26/meanChlamy_Col_26


Chlamy_Inv_32<-subset(alldat, who=="Chlamydomonas" & temp=="32" & scene=="invading")
meanChlamy_Inv_32<-mean(Chlamy_Inv_32$ROS_LD50)
sdChlamy_Inv_32<-sd(Chlamy_Inv_32$ROS_LD50)
Chlamy_Col_32<-subset(alldat, who=="Chlamydomonas" & temp=="32" & scene=="colonising")
meanChlamy_Col_32<-mean(Chlamy_Col_32$ROS_LD50)
meanChlamy_Inv_32/meanChlamy_Col_32



Chlamy_Inv_FS<-subset(alldat, who=="Chlamydomonas" & temp=="FS" & scene=="invading")
meanChlamy_Inv_FS<-mean(Chlamy_Inv_FS$ROS_LD50)
sdChlamy_Inv_FS<-sd(Chlamy_Inv_FS$ROS_LD50)
Chlamy_Col_FS<-subset(alldat, who=="Chlamydomonas" & temp=="FS" & scene=="colonising")
meanChlamy_Col_FS<-mean(Chlamy_Col_FS$ROS_LD50)
meanChlamy_Inv_FS/meanChlamy_Col_FS




Ostri_Inv<-subset(alldat, who=="Ostreococccus"  & scene=="invading")
meanOstri_Inv<-mean(Ostri_Inv$ROS_LD50)
sdOstri_Inv<-sd(Ostri_Inv$ROS_LD50)
Ostri_Col<-subset(alldat, who=="Ostreococccus" & scene=="colonising")
meanOstri_Col<-mean(Ostri_Col$ROS_LD50)
meanOstri_Inv/meanOstri_Col

Ostri_Inv_22<-subset(alldat, who=="Ostreococccus" & temp=="22" & scene=="invading")
meanOstri_Inv_22<-mean(Ostri_Inv_22$ROS_LD50)
sdOstri_Inv_22<-sd(Ostri_Inv_22$ROS_LD50)
Ostri_Col_22<-subset(alldat, who=="Ostreococccus" & temp=="22" & scene=="colonising")
meanOstri_Col_22<-mean(Ostri_Col_22$ROS_LD50)
meanOstri_Inv_22/meanOstri_Col_22

Ostri_Inv_26<-subset(alldat, who=="Ostreococccus" & temp=="26" & scene=="invading")
meanOstri_Inv_26<-mean(Ostri_Inv_26$ROS_LD50)
sdOstri_Inv_26<-sd(Ostri_Inv_26$ROS_LD50)
Ostri_Col_26<-subset(alldat, who=="Ostreococccus" & temp=="26" & scene=="colonising")
meanOstri_Col_26<-mean(Ostri_Col_26$ROS_LD50)
meanOstri_Inv_26/meanOstri_Col_26


Ostri_Inv_32<-subset(alldat, who=="Ostreococccus" & temp=="32" & scene=="invading")
meanOstri_Inv_32<-mean(Ostri_Inv_32$ROS_LD50)
sdOstri_Inv_32<-sd(Ostri_Inv_32$ROS_LD50)
Ostri_Col_32<-subset(alldat, who=="Ostreococccus" & temp=="32" & scene=="colonising")
meanOstri_Col_32<-mean(Ostri_Col_32$ROS_LD50)
meanOstri_Inv_32/meanOstri_Col_32



Ostri_Inv_FS<-subset(alldat, who=="Ostreococccus" & temp=="FS" & scene=="invading")
meanOstri_Inv_FS<-mean(Ostri_Inv_FS$ROS_LD50)
sdOstri_Inv_FS<-sd(Ostri_Inv_FS$ROS_LD50)
Ostri_Col_FS<-subset(alldat, who=="Ostreococccus" & temp=="FS" & scene=="colonising")
meanOstri_Col_FS<-mean(Ostri_Col_FS$ROS_LD50)
meanOstri_Inv_FS/meanOstri_Col_FS



##### and now the curves from the pilot - these figures are in the SI ####

temptrial<-read.csv("temptrial_OTchlamy.csv")
str(temptrial)
temptrial$Temp<-as.numeric(temptrial$Temp)

salt18<-qplot(Temp, Growthrate, data=temptrial,xlab="Assay Temperature (ºC)",ylab="Growth rate (µ d-1)", facets=.~Species)+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(Species~.)+ 
  geom_errorbar(aes(ymin=Growthrate-growth_sd, ymax=Growthrate+growth_sd)) +
  geom_smooth (aes(group=Species),method="gam",se=F,formula= y ~ s(x, k = 3))+
  theme(legend.position="top")
salt18


saltrial<-read.csv("sal_trial_OTchlamy.csv")
str(saltrial)
saltrial$Temp<-as.numeric(saltrial$Temp)
saltrial$Salinity<-as.numeric(saltrial$Salinity)


salt19<-qplot(Salinity, Growthrate, data=saltrial,xlab="Salinity (PSU)",ylab="Growth rate (µ d-1)", facets=.~Species, col=as.factor(Temp))+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(Species~.)+ 
  scale_colour_manual(values=c("lightblue","red","darkred"))+
  geom_errorbar(aes(ymin=Growthrate-growth_sd, ymax=Growthrate+growth_sd)) +
  geom_smooth (aes(group=Species),method="gam",se=F,formula= y ~ s(x, k = 3))+
  theme(legend.position="top")
salt19

######now as requested by reviewer  - read in the "3 point measurements" ####
saltra<-read.csv("trajecs_stable_gr.csv")
head(saltra)

salt20_22<-qplot(When, celllLOG, data=subset(saltra,temp=="22"),xlab="Time (days within batch cycle)",ylab="LOG cell count", facets=species~transf, col=Invasion,shape=temp)+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(species~transf)+ 
  scale_colour_manual(values=c("red","goldenrod1"))+
  geom_smooth(aes(group=Invasion), method=lm, se=FALSE)+
  theme(legend.position="top")
salt20_22

salt20_26<-qplot(When, celllLOG, data=subset(saltra,temp=="26"),xlab="Time (days within batch cycle)",ylab="LOG cell count", facets=species~transf, col=Invasion,shape=temp)+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(species~transf)+ 
  scale_colour_manual(values=c("red","goldenrod1"))+
  geom_smooth(aes(group=Invasion), method=lm, se=FALSE)+
  theme(legend.position="top")
salt20_26

salt20_32<-qplot(When, celllLOG, data=subset(saltra,temp=="32"),xlab="Time (days within batch cycle)",ylab="LOG cell count", facets=species~transf, col=Invasion,shape=temp)+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(species~transf)+ 
  scale_colour_manual(values=c("red","goldenrod1"))+
  geom_smooth(aes(group=Invasion), method=lm, se=FALSE)+
  theme(legend.position="top")
salt20_32

salt20_FS<-qplot(When, celllLOG, data=subset(saltra,temp=="22/32 rapid"),xlab="Time (days within batch cycle)",ylab="LOG cell count", facets=species~transf, col=Invasion,shape=temp)+
  theme_classic(base_size = 12, base_family = 'Helvetica') + 
  facet_grid(species~transf)+ 
  scale_colour_manual(values=c("red","goldenrod1"))+
  geom_smooth(aes(group=Invasion), method=lm, se=FALSE)+
  theme(legend.position="top")
salt20_FS

###putting all phenotype data together for forestplot ####

fps<-read.csv("forestplot_salt.csv")
str(fps)

fps$Temperature<-factor(fps$Temperature, levels=c('Mean','22','26','22/32','32'))
fps$Trait_det<-factor(fps$Trait_det, levels= c("ROS_end_cell_32", "ROS_end_cell_FS","ROS_end_cell_26" ,"ROS_end_cell_22", "ROS_end_cell" ,"ROS_detox_cell_32"  ,    
                                               "ROS_detox_cell_FS" ,      "ROS_detox_cell_26"   ,    "ROS_detox_cell_22"   ,   
                                               "ROS_detox_cell"   ,       "NP_end_32_decom"    ,     "NP_end_22/32_decom"    , 
                                               "NP_end_26_decom"   ,      "NP_end_22_decom"    ,     "NP_end"                , 
                                               "Size_end_32"     ,        "Size_end_FS"          ,   "Size_end_26"            ,
                                               "Size_end_22"  ,           "Size_end"             ,   "Growthrate_sal_evo_32"  ,
                                               "Growthrate_sal_evo_FS" ,  "Growthrate_sal_evo_26" ,  "Growthrate_sal_evo_22"  ,
                                               "Growthrate_sal_evo"  ,    "Growthrate_sal_short_32", "Growthrate_sal_short_FS",
                                               "Growthrate_sal_short_26" ,"Growthrate_sal_short_22" ,"Growthrate_sal_short"   ,
                                               "Survival_32" ,    "Survival_22/32" ,"Survival_26", "Survival_22" ,            "Survival"))


try_01<-ggplot(fps, aes(y = Trait_det, x = Effect_size, colour = as.factor(Temperature))) +
  geom_point(aes(size=.65)) + 
  geom_errorbarh( aes(y=Trait_det, xmin=Effect_size-ES_SE, xmax=Effect_size + ES_SE), colour="black", alpha=0.9, size=0.5, height=0.4) + 
  ggtitle(label = "Test plot")  +
  #  scale_colour_ordinal("")+ # I like this colour palette as it is colour blind friendly but feel free to change , could also do green for green algae, brown for browns etc
  scale_colour_manual(values=c("black","cornflowerblue", "darkblue","purple","red"))+
  xlab("Effect size") + ylab(" Variable ")+
  theme_classic(base_size = 12, base_family = "Helvetica")+
  facet_grid(~Species)+
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 35),lty=2, colour="black")+
  
  theme(legend.position="bottom")+
  theme(axis.text.x=element_text(angle=0.001, hjust=1))
try_01

