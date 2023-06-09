rm(list=ls()) # wipe R's memory

#~~~~~~~~~~~~~~~~~~~
# LOAD PACKAGES ####
#~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(zoo)
library(psych)
library(lme4)
library(rptR)
library(brms)
library(car)
library(DHARMa)
library(ggplot2)
library(coin)
library(gridExtra)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN DATA FILE 1: read in pheasants_30s-temps.csv (temperature readings every 30 s) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pheasants.30s <- read.csv(file.choose(),header=T) # read in data file
str(pheasants.30s) # check variables
pheasants.30s <- subset(pheasants.30s,aggrYN==0) # remove observations taken within 20 s of aggressive encounter
pheasants.30s$behaviour <- relevel(pheasants.30s$behaviour,ref="other") # set 'other' as baseline category

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN DATA FILE 2: read in pheasants_morphology.csv (morphological measurements at release) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pheasants.morph <- read.csv(file.choose(),header=T) # read in data file
str(pheasants.morph) # check variables
pheasants.30s <- left_join(pheasants.30s,pheasants.morph,by=c("pen","ID"))
summary(na.omit(pheasants.30s) %>%
          group_by(ID) %>%
          summarise(n=n())) # stats on number of measurements per pheasant

describe(pheasants.30s$max.temp) # basic descriptive stats
describe(pheasants.30s$spot.temp) # basic descriptive stats

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q1. Are temperature profiles individually repeatable? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model1 <- lmer(max.temp~behaviour+time+scale(spot.temp)+pen+(1|ID),
             data=na.omit(pheasants.30s))
summary(model1)
drop1(model1,test="Chisq")
simOut1<-simulateResiduals(fittedModel=model1,plot=F)
plot(simOut1)
rep1 <- rpt(max.temp~behaviour+scale(time)+scale(spot.temp)+pen+(1|ID),
          grname="ID",data=na.omit(pheasants.30s),datatype="Gaussian",nboot=1000,npermut=0)
summary(rep1)
plot(rep1,type="boot",grname="ID",cex.main=0.8)

# include individual-level variables
model2 <- lmer(max.temp~sex+scale(BC.release)+scale(tarsus.release)+behaviour+time+scale(spot.temp)+pen+(1|ID),
             data=na.omit(pheasants.30s))
summary(model2)
drop1(model2,test="Chisq")
anova(model1,model2)
simOut2<-simulateResiduals(fittedModel=model2,plot=F)
plot(simOut2)
rep2 <- rpt(max.temp~sex+scale(BC.release)+scale(tarsus.release)+behaviour+time+scale(spot.temp)+pen+(1|ID),
          grname="ID",data=na.omit(pheasants.30s),datatype="Gaussian",nboot=1000,npermut=0)
summary(rep2)
plot(rep2,type="boot",grname="ID",cex.main=0.8)

behavtemps <- pheasants.30s %>%
  group_by(ID,sex,behaviour) %>%
  summarise(temp.av = mean(max.temp,na.rm=T))
behavtemps[ is.na(behavtemps) ] <- NA # replace NaN with NA
fig2<-ggplot(subset(behavtemps,behaviour!="other"),aes(x=ordered(behaviour,levels=c("resting","standing","preening","foraging","walking")),y=temp.av,fill=sex))+
  geom_boxplot()+
  labs(y="maximum head temperature (°C)")+labs(x="behaviour")+
  theme_classic()+
  scale_fill_manual(values = c("greenyellow","grey40"))
fig2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# READ IN DATA FILE 3: read in pheasants_interactions.csv (temperature readings before/after attack) ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pheasants <- read.csv(file.choose(),header=T)
str(pheasants)
pheasants <- left_join(pheasants,pheasants.morph,by=c("pen","ID"))

#~~~~~~~~~~~~~~~~~~~~~~~
# DATA MANIPULATION ####
#~~~~~~~~~~~~~~~~~~~~~~~
pheasants$temp.obs <- pheasants$temp # save originally observed temp measurements
pheasants <- pheasants %>%
  group_by(contest) %>%
  mutate(temp = na.approx(temp,na.rm=FALSE,rule=1)) # interpolate missing temp values between observed values
baselinetemps <- pheasants.30s %>%
  group_by(ID) %>%
  summarise(temp.bl = mean(max.temp,na.rm=T)) # use average temp as baseline
pheasants <- left_join(pheasants,baselinetemps,by="ID") # merge with main data frame
attacktemps <- pheasants %>%
  filter(secs.prepost == 0) %>%
  select(contest,temp,ID) %>%
  group_by(contest,ID) %>%
  rename(t0.temp=temp) # extract temps at moment of attack
pheasants <- left_join(pheasants,attacktemps,by=c("contest","ID")) # merge with main data frame
av.pre <- pheasants %>%
  group_by(contest,ID) %>%
  filter(secs.prepost < 0) %>%
  summarise(pre.temp = mean(temp,na.rm=T)) %>%
  mutate(pre.temp=ifelse(is.nan(pre.temp),NA,pre.temp)) # aggregate temps before moment of attack
summarytemps <- left_join(av.pre,
                        subset(pheasants,secs.prepost==0,select=-c(temp,secs.prepost)),
                        by=c("contest","ID"))
av.post <- pheasants %>%
  group_by(contest,ID) %>%
  filter(secs.prepost > 0) %>%
  summarise(post.temp = mean(temp,na.rm=T)) %>%
  mutate(post.temp=ifelse(is.nan(post.temp),NA,post.temp)) # aggregrate temps after moment of attack
summarytemps <- left_join(summarytemps,av.post,by=c("contest","ID"))
pheasants$type <- relevel(pheasants$type,ref="peck")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q2. Do baseline temperatures predict aggressor and recipient roles? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
describeBy(summarytemps$temp.bl,summarytemps$role)
summarytemps$role <- relevel(summarytemps$role,ref="recipient")
model3 <- glmer(role~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+pen+(1|ID),
              data=summarytemps,family="binomial",control=glmerControl(optimizer="bobyqa"))
simOut3<-simulateResiduals(fittedModel=model3,plot=F)
plot(simOut3)
summary(model3)
drop1(model3,test="Chisq")
exp(fixef(model3))
exp(confint(model3))

summary(pheasants %>%
  group_by(contest,role,type,sex) %>%
  summarise(n=n()))

min.temp <- pheasants %>%
  group_by(contest,ID) %>%
  summarise(min.temp=min(temp,na.rm=T)) %>%
  mutate(min.temp=ifelse(is.infinite(min.temp),NA,min.temp))
summarytemps <- left_join(summarytemps,min.temp,by=c("contest","ID"))
max.temp <- pheasants %>%
  group_by(contest,ID) %>%
  summarise(max.temp=max(temp,na.rm=T)) %>%
  mutate(max.temp=ifelse(is.infinite(max.temp),NA,max.temp))
summarytemps <- left_join(summarytemps,max.temp,by=c("contest","ID"))
av.temp <- pheasants %>%
  group_by(contest,ID) %>%
  summarise(av.temp=mean(temp,na.rm=T)) %>%
  mutate(av.temp=ifelse(is.nan(av.temp),NA,av.temp))
summarytemps <- left_join(summarytemps,av.temp,by=c("contest","ID"))
describeBy(summarytemps$max.temp,summarytemps$role)
describeBy(summarytemps$min.temp,summarytemps$role)
describeBy(summarytemps$av.temp,summarytemps$role)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q3. Do thermal profiles differ between aggressors and recipients prior to attack? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pheasants$role <- relevel(pheasants$role,ref="recipient")
pheasants$sex <- relevel(pheasants$sex,ref="M")
pheasants$type <- relevel(pheasants$type,ref="peck")
model3pre <- lmer(temp~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost*role+pen+(1|ID)+(1|contest),
                data=subset(pheasants,secs.prepost<=0),control=lmerControl(optimizer="bobyqa"))
summary(model3pre)
drop1(model3pre,test="Chisq")
# drop NS interaction term
model3pre <- lmer(temp~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost+role+pen+(1|ID)+(1|contest),
                data=subset(pheasants,secs.prepost<=0),control=lmerControl(optimizer="bobyqa"))
summary(model3pre)
drop1(model3pre,test="Chisq")
simOut3pre<-simulateResiduals(fittedModel=model3pre,plot=F)
plot(simOut3pre)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Q4. Do thermal profiles differ between aggressors and recipients following an attack? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
model4post <- lmer(temp~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost*role+I(secs.prepost^2)*role+pen+(1|ID)+(1|contest),
                   data=subset(pheasants,secs.prepost>0),control=lmerControl(optimizer="bobyqa"))
summary(model4post)
drop1(model4post,test="Chisq")
# drop interaction terms
model5post <- lmer(temp~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost+I(secs.prepost^2)+role+pen+(1|ID)+(1|contest),
                   data=subset(pheasants,secs.prepost>0),control=lmerControl(optimizer="bobyqa"))
summary(model5post)
drop1(model5post,test="Chisq")
anova(model5post,model4post)
# drop quadratic term
model6post <- lmer(temp~scale(temp.bl)+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost+role+pen+(1|ID)+(1|contest),
                 data=subset(pheasants,secs.prepost>0),control=lmerControl(optimizer="bobyqa"))
summary(model6post)
drop1(model6post,test="Chisq")
anova(model6post,model5post)
simOut5post<-simulateResiduals(fittedModel=model5post,plot=F)
plot(simOut5post)

# Does peak temp and its timing differ between aggressors and recipients?
pheasants$role <- relevel(pheasants$role,ref="aggressor")
maxtimes<-pheasants %>%
  group_by(contest,ID) %>%
  slice(which.max(temp))
oneway_test(temp~role,data=maxtimes,distribution=approximate(10000))
oneway_test(secs.prepost~role,data=maxtimes,distribution=approximate(10000))

# Does min temp and its timing differ between aggressors and recipients?
pheasants$role <- relevel(pheasants$role,ref="aggressor")
mintimes<-pheasants %>%
  group_by(contest,ID) %>%
  slice(which.min(temp))
oneway_test(temp~role,data=mintimes,distribution=approximate(10000))
oneway_test(secs.prepost~role,data=mintimes,distribution=approximate(10000))

#~~~~~~~~~~~~~~~~~~
# GAMM version ####
#~~~~~~~~~~~~~~~~~~
pheasants$role <- relevel(pheasants$role,ref="aggressor")
modelb1 <- brm(temp~temp.bl+sex+scale(tarsus.release)+scale(BC.release)+type+s(secs.prepost,by=role)+role+pen+(1|ID)+(1|contest),
            data=pheasants,family=gaussian(),
            cores=4,seed=17,control=list(adapt_delta=0.99))
summary(modelb1)
ceffs1<-conditional_effects(modelb1,effects="secs.prepost:role")
csms1<-conditional_smooths(modelb1)
plot(csms1)

modelb2 <- brm(temp~temp.bl+sex+scale(tarsus.release)+scale(BC.release)+type+s(secs.prepost)+role+pen+(1|ID)+(1|contest),
               data=pheasants,family=gaussian(),
               cores=4,seed=17,control=list(adapt_delta=0.99))
summary(modelb2)
ceffs2<-conditional_effects(modelb2,effects="secs.prepost:role")
csms2<-conditional_smooths(modelb2)
plot(csms2)

LOO(modelb1,modelb2) # LOOIC comparison for time*role interaction

modelb3 <- brm(temp~temp.bl+sex+scale(tarsus.release)+scale(BC.release)+type+secs.prepost+role+pen+(1|ID)+(1|contest),
               data=pheasants,family=gaussian(),
               cores=4,seed=17,control=list(adapt_delta=0.99))
summary(modelb3)
ceffs3<-conditional_effects(modelb3,effects="secs.prepost:role")
plot(ceffs3)

LOO(modelb2,modelb3) # LOOIC comparison for non-linear time effect

#~~~~~~~~~~~
# PLOTS ####
#~~~~~~~~~~~
plot.df<-ceffs1$`secs.prepost:role` %>%
  arrange(role,secs.prepost)
conf.bands<-csms1$`mu: s(secs.prepost,by=role)` %>%
  arrange(role,secs.prepost) %>%
  select(role,secs.prepost,estimate__,lower__,upper__) %>%
  rename(estimateX=estimate__) %>%
  rename(lowerX=lower__) %>%
  rename(upperX=upper__)
plot.df <- left_join(plot.df,conf.bands,by=c("role","secs.prepost"))
plot.df$estimateX <- plot.df$estimate__ + plot.df$estimateX 
plot.df$lowerX <- plot.df$estimate__ + plot.df$lowerX 
plot.df$upperX <- plot.df$estimate__ + plot.df$upperX 

fig3traj <- ggplot(plot.df,aes(x=secs.prepost,y=estimateX)) +
  geom_vline(xintercept=0)+
  geom_smooth(aes(ymin=lowerX,ymax=upperX,fill=role,colour=role),stat="identity",alpha=0.3)+
  xlab('time (seconds before or after attack)') + ylab('maximum head temperature (°C)') +
  scale_y_continuous(breaks=seq(36.0,37.4,0.2))+
  scale_colour_manual(values = c("darkgoldenrod1","navyblue")) +
  scale_fill_manual(values = c("darkgoldenrod1","navyblue")) +
  theme_classic()+
  theme(legend.position=c(0.7,0.2),legend.background=element_rect(fill="transparent"),legend.title=element_blank())
fig3traj

fig3mintimes <- ggplot(mintimes,aes(x=secs.prepost,fill=role))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  scale_fill_manual(values = c("darkgoldenrod1","navyblue"))+
  xlab('XXX')+
  scale_y_reverse()+
  theme(axis.text.y=element_text(colour="white"),axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.title.y=element_text(colour="white"))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.title.x=element_blank())+
  theme(legend.position="none")
fig3maxtimes <- ggplot(maxtimes,aes(x=secs.prepost,fill=role))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  scale_fill_manual(values = c("darkgoldenrod1","navyblue"))+
  xlab('XXX')+
  theme(axis.text.y=element_text(colour="white"),axis.ticks.y=element_blank(),axis.line.y=element_blank(),axis.title.y=element_text(colour="white"))+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.title.x=element_blank())+
  theme(legend.position="none")
fig3mintemps <- ggplot(mintimes,aes(x=temp,fill=role))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  ylab('min')+
  theme(axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.text.x=element_text(colour="white"))+
  theme(axis.title.y=element_blank())+
  scale_fill_manual(values = c("darkgoldenrod1","navyblue"))+
  theme(legend.position="none")+
  coord_flip(xlim=c(34.5,39.0))
fig3maxtemps <- ggplot(maxtimes,aes(x=temp,fill=role))+
  geom_boxplot(alpha=0.3)+
  theme_classic()+
  ylab('max')+
  theme(axis.ticks.x=element_blank(),axis.line.x=element_blank(),axis.text.x=element_text(colour="white"))+
  theme(axis.title.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.line.y=element_blank())+
  scale_fill_manual(values = c("darkgoldenrod1","navyblue"))+
  theme(legend.position="none")+
  coord_flip(xlim=c(34.5,39.0))
blankPlot <- ggplot()+geom_blank()+ theme_void()
grid.arrange(fig3maxtimes, blankPlot, blankPlot, fig3traj, fig3mintemps, fig3maxtemps, fig3mintimes, blankPlot, blankPlot, ncol=3, nrow=3, widths=c(2,0.38,0.3), heights=c(3, 15, 3))