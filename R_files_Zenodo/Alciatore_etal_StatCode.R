#set directory
setwd('/Users/Desktop/Obiroi_ImmuneChallenge')

#clear workspace
rm(list=ls())

#load packages
library(survival)
library(multcomp)
library(glmmTMB)
library(DHARMa)
library(ggplot2)  
library(survminer)
library(coxme)
library(emmeans)
library(EnvStats)
library(readxl)
library(nlme)
library(plyr)
library(gridExtra)
library(igraph)
library(vegan)
library(PerformanceAnalytics)
library(adegenet)
library(reshape2)
library(randomForest)

#Define palette plots
colTr="#0072B2"  
colSh="#D55E00"
colNV="#006400"

#Create folder to save plots
ifelse(!dir.exists(file.path(getwd(),'Plots')), dir.create(file.path(getwd(),'Plots')), FALSE)

#Create folder to save confidence intervals
ifelse(!dir.exists(file.path(getwd(),'BootstrappedCI')), dir.create(file.path(getwd(),'BootstrappedCI')), FALSE)

###Fig.S1, Pathogen Exposure Survival###
obiroi.survival=read.table(file="FigS1_PathogenExposureSurvival/PathogenExposure.survival.txt", header=T, sep="\\t")

#Create survival curves for plotting
fit=survfit(Surv(day.of.death, status) ~ treatment.rearing, data = obiroi.survival) 

#Plot
obiroiSurv=ggsurvplot(
  fit,                     
  data = obiroi.survival,  # data used to fit survival curves. 
  risk.table = FALSE,       # show risk table.
  pval = FALSE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals 
  break.time.by = 1,     
  palette=c(colTr,colTr,colSh,colSh),
  censor=FALSE,
  linetype=c(1,3,1,3),
  size=0.3,
  legend="none",
  ggtheme = theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  text=element_text(size=10,family="Helvetica"),axis.text = element_text(size=9)),
  risk.table.y.text.col = FALSE, 
  risk.table.y.text = FALSE, 
  xlab=('Days after Treatment'),
  ylab=('Proportion Surviving')
)
ggsave(filename="Plots/FigS1_PathogenExpSurv.pdf",print(obiroiSurv),width=3.625,height=2,unit="in")

#Statistical model
cox_Surv=coxme(Surv(day.of.death, status) ~ treatment*rearing+age + (1|replicate), data = obiroi.survival)
anova(cox_Surv) #age has no effect, removing from future models and experiments

cox_Surv=coxme(Surv(day.of.death, status) ~ treatment*rearing + (1|replicate), data = obiroi.survival)
anova(cox_Surv)

pairwise_survdiff(Surv(day.of.death, status) ~ treatment.rearing, data = obiroi.survival, p.adjust.method = "BH" ) #pairwise comparisons

confint(glht(cox_Surv)) #confidence intervals


###Fig.S2 Pathogen Exposure Allogrooming###
PathogenExpGrooming<-read.table('FigS2_PathogenExposureGrooming/PathogenExposure.grooming.txt',header=TRUE,na.strings="NA",""," ",fill=TRUE) 

PathogenExpGrooming$treatment=as.factor(PathogenExpGrooming$treatment) #Treatment as factor
PathogenExpGrooming$age=as.factor(PathogenExpGrooming$age) #Young/Old
PathogenExpGrooming$allogroom.time<-as.numeric(PathogenExpGrooming$allogroom.time) #Duration allogrooming
PathogenExpGrooming$timenum=as.numeric(PathogenExpGrooming$time) #Time of observation
PathogenExpGrooming$rep<-as.factor(PathogenExpGrooming$rep) #Batch
PathogenExpGrooming$time<-as.factor(PathogenExpGrooming$time) #Time of observation as factor

PathogenExpGrooming$treat_time<-interaction(PathogenExpGrooming$treatment,PathogenExpGrooming$time) # 90-level variable coding treatment at each time
PathogenExpGrooming$indID<-interaction(PathogenExpGrooming$rep,PathogenExpGrooming$age,PathogenExpGrooming$treatment) # 56-level variable coding individual ID

SubPathogenGrooming<-subset(PathogenExpGrooming,allogroom.time!='NA') #subset for which allogrooming data is available (final=427 lines)

SubPathogenGrooming<-droplevels(SubPathogenGrooming) # drop unused levels (18 levels for treat_time, 48 levels for indID)
SubPathogenGrooming$index<-as.factor(seq_len(nrow(SubPathogenGrooming)))  #include index unique for each observation

#Statistical Model
mPEG<-glmmTMB(allogroom.time~treatment*time+(1|rep/indID), data=SubPathogenGrooming, dispformula=~ time*treatment,
              start=list(thetaf=1.5), family = tweedie)
plot(simulateResiduals(mPEG)) #Diagnostic model

mPEGNoInt<-glmmTMB(allogroom.time~treatment+time+(1|rep/indID), data=SubPathogenGrooming, dispformula=~ treatment+time,
                       start=list(thetaf=1.565), family = tweedie)
anova(mPEG,mPEGNoInt) #compare to model without interaction treatment*time

emmeans(mPEG, specs = pairwise ~ treatment|time)

#Plot
SubPathogenGroomingPlot=  ggplot(aes(x=time, y=allogroom.time, color=treatment, shape=treatment), data=SubPathogenGrooming)+
    geom_point(data=SubPathogenGrooming,aes(x=time,y=allogroom.time,color=treatment,shape=treatment),position=position_jitterdodge(jitter.width=0.6),size=1,alpha=0.5)+
    stat_summary(fun=mean, geom="point",position=position_nudge(x=c(-0.2,0.2)),size=2)+
    stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.15,position=position_nudge(x=c(-0.2,0.2))) +
    scale_y_continuous(breaks = c(0,120,240,360,480,600) ) +stat_n_text(size=2.5) + 
    theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white",colour="white"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          text=element_text(size=8,family="Helvetica"),axis.text = element_text(size=9),legend.position="none") + scale_color_manual(values=c(colTr,colSh))+
  xlab("Time post-treatment (hours)")+ylab("Received grooming (seconds)")

ggsave(filename="Plots/FigS2_PathogenExpGrooming.pdf",print(SubPathogenGroomingPlot),width=3.625,height=2,unit="in")
plotS2_info <- ggplot_build(SubPathogenGroomingPlot)$`data`[[3]] #Extract confidence intervals from ggplot
write.csv(plotS2_info, "BootstrappedCI/TableS1_PathogenExpAllogrooming.csv", row.names = FALSE) #Print confidence intervals (Table S1)

###Fig.1 Immune Challenge Allogrooming###
ImmChGrooming <- read.table("Fig1_ImmuneChallengeGrooming/ImuneChallengeGrooming.csv",header=TRUE,sep=",",na.strings="NA",""," ",fill=TRUE) # 384 obs.

ImmChGrooming$treatment=as.factor(ImmChGrooming$treatment) #treatment as factor
ImmChGrooming$colony<-as.factor(ImmChGrooming$colony) #colony as factor
ImmChGrooming$timepoint<-as.factor(ImmChGrooming$timepoint) #time of observation
ImmChGrooming$index<-as.factor(seq_len(nrow(ImmChGrooming)))  #include index unique for each observation
ImmChGrooming$individual<-as.factor(ImmChGrooming$individual)
SubImmChGrooming<-subset(ImmChGrooming,duration!='NA') ## data subset for which allogrooming data is available (final=320, 24 colonies)

#Statistical Model 
mICG<-glmmTMB(duration~treatment*timepoint+(1|colony/individual),family=tweedie,data=SubImmChGrooming,dispformula=~ treatment,start=list(thetaf=1.569))
plot(simulateResiduals(mICG))
drop1(mICG, test="Chisq")

emmeans(mICG, specs = pairwise ~ treatment|timepoint)

#Plot
ImmChallGrooming=ggplot(aes(x=timepoint, y=duration, color=treatment, shape=treatment), data=SubImmChGrooming)+
  geom_point(data=SubImmChGrooming,aes(x=timepoint,y=duration,color=treatment,shape=treatment),position=position_jitterdodge(jitter.width=0.6),size=1,alpha=0.5)+
  stat_summary(fun=mean, geom="point",position=position_nudge(x=c(-0.2,0.2)),size=2)+
  stat_summary(fun.data=mean_cl_boot, geom="errorbar", width=0.15,position=position_nudge(x=c(-0.2,0.2))) +
  scale_y_continuous(breaks = c(0,120,240,360,480,600) ) + 
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white",colour="white"), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        text=element_text(size=9,family="Helvetica"),axis.text = element_text(size=9),legend.position="none")+
  scale_colour_manual(values=c(colTr,colSh)) + xlab("Time post-treatment (hours)")+ylab("Received grooming (min)")

ggsave(filename="Plots/Fig1_ImmuneChallengeGrooming.pdf",print(ImmChallGrooming),width=3.54,height=2,unit="in")
plot2_info <- ggplot_build(ImmChallGrooming)$`data`[[3]] #Extract confidence intervals from ggplot
write.csv(plot2_info, "BootstrappedCI/TableS2_ImmuneChallengeAllogrooming.csv", row.names = FALSE) #Print confidence intervals (Table S2)

####Fig.2 Behaviours###

ImmChall_IndBehav_All <- read_excel("Fig2_Behaviours/IndividualBehavioursData.xlsx") 

SubDataIndvBehav=subset(ImmChall_IndBehav_All,Treatment!='N')  #no naive ants, only take the treated ants

SubDataIndvBehav$Isolation<-as.numeric(SubDataIndvBehav$Isolation)/36000 #As proportion of time 
SubDataIndvBehav$Activity<-as.numeric(SubDataIndvBehav$Activity)/36000 #As proportion of time
SubDataIndvBehav$meanSpeed<-as.numeric(SubDataIndvBehav$meanSpeed) #mm/s

SubDataIndvBehav$Colony<-as.factor(SubDataIndvBehav$Colony) #Colony as factor
SubDataIndvBehav$TimeStampNum<-SubDataIndvBehav$TimeStamp #Time of observation
SubDataIndvBehav$TimeStamp<-as.factor(SubDataIndvBehav$TimeStamp) #Time of observation as factor
SubDataIndvBehav$Treatment<-as.factor(SubDataIndvBehav$Treatment) #Immune challenged vs control injected
SubDataIndvBehav$treat_time<-interaction(SubDataIndvBehav$Treatment,SubDataIndvBehav$TimeStamp)
SubDataIndvBehav$TimeStampPlot<-SubDataIndvBehav$TimeStampNum #Time of observation, used in the plot
SubDataIndvBehav$TimeStampPlot[SubDataIndvBehav$TimeStampNum==6]<-3


##ISOLATION Fig2A
#Statistical model
m1lme_a<-lme(Isolation~ Treatment * TimeStamp, random=~1|Colony/Individual,data=SubDataIndvBehav,method = "ML")
drop1(m1lme_a,test="Chisq") 


#Validation model
#1) Outliers
boxplot(Isolation ~ Treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "Isolation", data = SubDataIndvBehav)

#2) Normality residuals
E <- resid(m1lme_a) 
hist(E)
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) #has to be <4

#3)Homogeneity of variance
#Plot residuals
plot(m1lme_a, add.smooth = FALSE, which = 1)

op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
hist(E, xlab = "Residuals", main = "")
plot(SubDataIndvBehav$Isolation, E, xlab = "Alone",
     ylab = "Residuals")
boxplot(E ~ SubDataIndvBehav$Treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp * SubDataIndvBehav$Treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$Colony, main = "Colony")
abline(0, 0)
par(op)
par(mfrow = c(1, 1))

#Introduce variance structure
m1lme_b<-lme(Isolation~ Treatment*TimeStamp ,random=~1|Colony/Individual,weights=varIdent(form =~ 1|Treatment*TimeStamp),
             data=SubDataIndvBehav,control=list(msMaxIter=100,opt="optim",msVerbose=TRUE),method = "ML")

E2 <- resid(m1lme_b,type="normalized") 
coplot(E ~ Isolation | Treatment * TimeStampNum, data = SubDataIndvBehav,
       ylab = "Normalised residuals")

#Final statistics
AIC(m1lme_a,m1lme_b) #Variance structure does improve the model
drop1(m1lme_b,test="Chisq")
emmeans(m1lme_b, specs = pairwise ~ Treatment|TimeStamp)

##ACTIVITY Fig2b
#Statistical model
m2lme_a<-lme(Activity~ Treatment * TimeStamp, random=~1|Colony/Individual,data=SubDataIndvBehav,method = "ML")
drop1(m2lme_a,test="Chisq") 

##Validation model
#1) Outliers
boxplot(Activity ~ Treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "Activity", data = SubDataIndvBehav)

#2) Normality residuals
E <- resid(m2lme_a)
hist(E)
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) #has to be <4


#3)Homogeneity of variance
#Plot residuals

op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
plot(m2lme_a, add.smooth = FALSE, which = 1)
hist(E, xlab = "Residuals", main = "")
plot(SubDataIndvBehav$Activity, E, xlab = "Active",
     ylab = "Residuals")
boxplot(E ~ SubDataIndvBehav$Treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp * SubDataIndvBehav$Treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$Colony, main = "Colony")
abline(0, 0)
par(op)
par(mfrow = c(1, 1))

#Introduce variance structure
m2lme_b<-lme(Activity~ Treatment*TimeStamp ,random=~1|Colony/Individual,weights=varIdent(form =~ 1|Treatment*TimeStamp),
             data=SubDataIndvBehav,control=list(msMaxIter=100,opt="optim",msVerbose=TRUE),method = "ML")

#Check new normalized residuals
E2 <- resid(m2lme_b,type="normalized") 
plot(m2lme_b, add.smooth = FALSE, which = 1) 
coplot(E2 ~ Activity | Treatment * TimeStampNum, data = SubDataIndvBehav,
       ylab = "Normalised residuals")

#Final statistics
AIC(m2lme_a,m2lme_b)
drop1(m2lme_b,test="Chisq")
emmeans(m2lme_b, specs = pairwise ~ Treatment|TimeStamp)

##SPEED Fig.2c
#Statistical model
m3lme_a<-lme(meanSpeed~ Treatment * TimeStamp, random=~1|Colony/Individual,data=SubDataIndvBehav,method = "ML")
drop1(m3lme_a,test="Chisq") 

##Validation model
#1) Outliers
boxplot(meanSpeed ~ Treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "meanSpeed", data = SubDataIndvBehav)


#2)Normality residuals
E <- resid(m3lme_a)
hist(E)
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) #has to be <4

#3)Homogeneity of variance
#Plot residuals

op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
plot(m3lme_a, add.smooth = FALSE, which = 1)
hist(E, xlab = "Residuals", main = "")
plot(SubDataIndvBehav$meanSpeed, E, xlab = "meanSpeed",
     ylab = "Residuals")
boxplot(E ~ SubDataIndvBehav$Treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$TimeStamp * SubDataIndvBehav$Treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ SubDataIndvBehav$Colony, main = "Colony")
abline(0, 0)
par(op)
par(mfrow = c(1, 1))

#Introduce variance stricture
m3lme_b<-lme(meanSpeed~ Treatment*TimeStamp ,random=~1|Colony/Individual,weights=varIdent(form =~ 1|Treatment*TimeStamp), ##All variance structures also selected trough AIC comparison w/ simpler structures
             data=SubDataIndvBehav,control=list(msMaxIter=100,opt="optim",msVerbose=TRUE),method = "ML")

E2 <- resid(m3lme_b,type="normalized") 
plot(m3lme_b, add.smooth = FALSE, which = 1) 
coplot(E2 ~ meanSpeed | Treatment * TimeStampNum, data = SubDataIndvBehav,
       ylab = "Normalised residuals")

#Final statistics
AIC(m3lme_a,m3lme_b)
drop1(m3lme_b,test='Chisq')
emmeans(m3lme_b, specs = pairwise ~ Treatment|TimeStamp)

##Create mean and 95% CI
SubDataPlots=ddply(ImmChall_IndBehav_All, .(NameColony, Colony, TimeStamp,Treatment, Colony_Time), summarize, Isolation=mean(Isolation),Activity=mean(Activity),meanSpeed=mean(meanSpeed))

SubDataPlots$Isolation<-as.numeric(SubDataPlots$Isolation)/36000 #As proportion of time 
SubDataPlots$Activity<-as.numeric(SubDataPlots$Activity)/36000 #As proportion of time
SubDataPlots$meanSpeed<-as.numeric(SubDataPlots$meanSpeed) #mm/s
SubDataPlots$Colony<-as.factor(SubDataPlots$Colony) 
SubDataPlots$TimeStampNum<-as.numeric(SubDataPlots$TimeStamp)
SubDataPlots$TimeStamp<-as.factor(SubDataPlots$TimeStamp) 
SubDataPlots$Treatment<-as.factor(SubDataPlots$Treatment)
SubDataPlots$Treatment <- relevel(SubDataPlots$Treatment, "Z")
SubDataPlots$Treatment <- relevel(SubDataPlots$Treatment, "C")
SubDataPlots$Treatment <- relevel(SubDataPlots$Treatment, "N")
SubDataPlots$TimeStampPlot<-SubDataPlots$TimeStampNum
SubDataPlots$TimeStampPlot[SubDataPlots$TimeStampNum==6]<-3
metrics=c(6,7,8)


MeansIndBehav=data.frame(TimeStamp=character(),TimeStampNum=numeric(),Treatment=character(),Metric=character(),mean=double(),meanboot=double(),low=double(),high=double(),stringsAsFactors=FALSE)
new_row=list(TimeStamp=character(),TimeStampNum=numeric(),Treatment=character(),Metric=character(),mean=double(),meanboot=double(),low=double(),high=double())
             
for (tr in unique(SubDataPlots$Treatment)) {
  for (tw in unique(SubDataPlots$TimeStamp)){
    SubSet=subset(SubDataPlots,Treatment==tr&TimeStamp==tw)
    l=length(SubSet)
    new_row[1]=tw
    new_row[2]=unique(SubSet$TimeStampNum)
    new_row[3]=tr
    for (met in metrics){
      new_row[4]=names(SubSet[met])
      mat=as.matrix(SubSet[,met])
      new_row[5]=mean(mat)
      bstrap <- c()
      for (i in 1:1000){ bstrap <- c(bstrap, mean(sample(mat,l,replace=T)))}
      new_row[6]=mean(bstrap)
      new_row[7]=quantile(bstrap,.025)
      new_row[8]=quantile(bstrap,.975)
      MeansIndBehav[nrow(MeansIndBehav)+1,]=new_row
      }
  }
}

MeansIndBehav$TimeStamp<-as.factor(MeansIndBehav$TimeStamp) 
MeansIndBehav$Treatment<-as.factor(MeansIndBehav$Treatment)
MeansIndBehav$Treatment <- relevel(MeansIndBehav$Treatment, "Z")
MeansIndBehav$Treatment <- relevel(MeansIndBehav$Treatment, "C")
MeansIndBehav$Treatment <- relevel(MeansIndBehav$Treatment, "N")
MeansIndBehav$TimeStampPlot<-as.numeric(as.character(MeansIndBehav$TimeStamp))
MeansIndBehav$TimeStampPlot[MeansIndBehav$TimeStamp==6]<-3

write.csv(MeansIndBehav, "BootstrappedCI/TableS3ABC_Behaviours.csv", row.names = FALSE)

##Plots
act.plot <- 
  subset(MeansIndBehav,Metric=="Activity") %>% 
  ggplot(aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=SubDataPlots,aes(x=TimeStampPlot+3,y=Activity,color=Treatment, shape=Treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='N'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='C'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='Z'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='N'))+
  geom_line(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='C'))+
  geom_line(data=subset(MeansIndBehav,Metric=="Activity"&Treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=0.5)+
  ylab("Proportion of time active")+xlab("Time post-treatment (hours)")+
theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
     panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=9,family="Helvetica"),
     axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle")) 

alone.plot <- 
  subset(MeansIndBehav,Metric=="Isolation") %>% 
  ggplot(aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=SubDataPlots,aes(x=TimeStampPlot+3,y=Isolation,color=Treatment, shape=Treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='N'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='C'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='Z'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='N'))+
  geom_line(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='C'))+
  geom_line(data=subset(MeansIndBehav,Metric=="Isolation"&Treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=0.5)+
  ylab("Proportion of time isolated")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=9,family="Helvetica"),
        axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle")) 

speed.plot <- 
  subset(MeansIndBehav,Metric=="meanSpeed") %>% 
  ggplot(aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=SubDataPlots,aes(x=TimeStampPlot+3,y=meanSpeed,color=Treatment, shape=Treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='N'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='C'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  geom_point(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='Z'),size=1,aes(x=TimeStampPlot+3, y=mean, color=Treatment, shape=Treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='N'))+
  geom_line(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='C'))+
  geom_line(data=subset(MeansIndBehav,Metric=="meanSpeed"&Treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=0.5)+
  ylab("Mean speed (mm/s)")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=9,family="Helvetica"),
        axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle")) 

ggsave(filename="Plots/Fig2ABC_Behaviours.pdf",plot=grid.arrange(alone.plot,act.plot,speed.plot,nrow=1),width=6.215,height=3,unit="in")

##Networks##
#Networks based on 6h of video, distance between ants=1mm 
exp_dir='Fig2_Behaviours/Network_1mmDist6h' 
Colonies=dir(exp_dir)  

#Extract matrices
UberMatList<-sapply(Colonies,function(g){  ##List of matrices
  col_dir<-file.path(exp_dir,g,'/')
  myFiles<-list.files(col_dir,pattern="*InteractionTab*")
  GraphList<- sapply(myFiles, function(i){
    path=file.path(col_dir,i)
    mat=read.csv(path,header = TRUE,row.names = 1)
    return(mat)}
    ,simplify = FALSE,USE.NAMES = TRUE)}
  ,simplify = FALSE,USE.NAMES = TRUE)

#Create networks
UberGraphList<-sapply(Colonies,function(g){  
  col_dir<-file.path(exp_dir,g,'/')
  myFiles<-list.files(col_dir,pattern="*InteractionTab*")
  GraphList<- sapply(myFiles, function(i){
    path=file.path(col_dir,i)
    mat=read.csv(path,header = TRUE,row.names = 1)
    graph_from_adjacency_matrix(as.matrix(mat),mode ="undirected", weighted = TRUE, diag = FALSE,add.rownames = TRUE)}
      ,simplify = FALSE,USE.NAMES = TRUE)}
  ,simplify = FALSE,USE.NAMES = TRUE)

##Get centralities
Centralities=data.frame(colony=character(),matrix=double(),hour=double(),timewindow=double(),ID_colony=character(),tag=character(),
                        treatment=character(),eigen_centrality=double(),degree=double(),strength=double(),stringsAsFactors=FALSE)
IDs=names(UberMatList[[1]][[1]])
new_row=list(colony=character(),matrix=double(),hour=double(),timewindow=double(),ID_colony=character(),tag=character(),
             treatment=character(),eigen_centrality=double(),degree=double(),strength=double())

for (col in 1:length(Colonies)) {
  name_col=unlist(names(UberGraphList[col]))
  new_row[1]=as.character(name_col)
  NetList=UberGraphList[[col]][]
  MatList=UberMatList[[col]][]
  for (mat in 1:length(NetList)){
    net=NetList[[mat]]
    name_mat=unlist(strsplit(names(MatList[mat]),"_"))
    new_row[2]=name_mat[4]
    new_row[3]=name_mat[6]
    time=as.numeric(name_mat[6])
    for (ant in c(1:9)){
      ID=IDs[ant]
      new_row[4]=name_mat[4]
      new_row[5]=as.character(paste(ID,name_col,sep="_"))
      new_row[6]=as.character(ID)
      if (ID=="BB"){new_row[7]="Z"} else if (ID=="RR"){new_row[7]="PBS"} else {new_row[7]="N"}
      new_row[8]=eigen_centrality(net)$vector[[ant]]
      new_row[9]="NA"  
      if (time== 6 || time==51) {new_row[10]=log((strength(net)[[ant]])/3)} #Correct for timewindows' duration
      else {new_row[10]=log((strength(net)[[ant]])/6)}
      Centralities[nrow(Centralities)+1,]=new_row
    }
  }
} 

Centralities$tag=as.factor(Centralities$tag)
Centralities$treatment=as.factor(Centralities$treatment)
Centralities$hour=as.numeric(as.character(Centralities$hour))
Centralities$timewindownumeric=as.double(Centralities$timewindow)
Centralities$timewindow=as.factor(Centralities$timewindow)
Centralities$colony=as.factor(Centralities$colony)

EigenVectorCentrality=subset(Centralities,treatment!='N')
EigenVectorCentrality=rbind(EigenVectorCentrality,subset(Centralities,tag=="BG"))

NodeStrengthCentrality=ddply(Centralities, .(colony, matrix, hour,timewindow, treatment, timewindownumeric), summarize, strength=mean(strength))

CentralitiesStats=droplevels(subset(Centralities,treatment!='N'))

##Get data from interactions
AllInteractions=data.frame(colony=character(),timewindow=double(),hour=double(),ID_colony=character(),tag=character(),treatment=character(),ant1=double(),ant2=double(),ant3=double(),ant4=double(),ant5=double(),ant6=double(),ant7=double(),ant8=double(),ID_int=character(),Max_int=double(),TOTint=double(),meanint=double(),specialization=double(),stringsAsFactors=FALSE)
IDs=names(UberMatList[[1]][[1]])
new_row=list(colony=character(),timewindow=double(),hour=double(),ID_colony=character(),tag=character(),treatment=character(),ant1=double(),ant2=double(),ant3=double(),ant4=double(),ant5=double(),ant6=double(),ant7=double(),ant8=double(),ID_int=character(),Max_int=double(),TOTint=double(),meanint=double(),specialization=double())

for (col in 1:length(Colonies)) {
  name_col=unlist(names(UberMatList[col]))
  new_row[1]=as.character(name_col) # Colony
  Matriceslist=UberMatList[[col]][] 
  for (mat in 1:length(Matriceslist)){
    name_mat=unlist(strsplit(names(Matriceslist[mat]),"_"))
    new_row[2]=name_mat[4] # matrix
    new_row[3]=name_mat[6] # hour
    for (ant1 in c(1:9)){
      ID=names(Matriceslist[[mat]][ant1])
      new_row[4]=as.character(paste(ID,name_col,sep="_")) # ID_colony
      new_row[5]=as.character(ID) # tag
      if (ID=="BB"){new_row[6]="Z"} else if (ID=="RR"){new_row[6]="PBS"} else {new_row[6]="N"} # treatment
      interactions=Matriceslist[[mat]][ant1,]
      interactions=interactions[-ant1]
      new_row[7:14]=interactions
      new_row[15]=names(which.max(interactions)) #ID_int
      new_row[16]=max(interactions) # Max interaction in seconds
      new_row[17]=sum(interactions) # Sum of all interactions
      new_row[18]=mean(as.double(interactions)) # Mean interaction
      new_row[19]=(max(interactions)/sum(interactions))#
      AllInteractions[nrow(AllInteractions)+1,]=new_row
    }
  }
}

AllInteractions$treatment=as.factor(AllInteractions$treatment)
AllInteractions$hour=as.factor(AllInteractions$hour)
AllInteractions$timewindow=as.factor(AllInteractions$timewindow)

AverageInteractions=ddply(AllInteractions, .(colony, hour,timewindow, treatment), summarize, Max_int=mean(Max_int),TOTint=mean(TOTint),meanint=mean(meanint),specialization=mean(specialization))
AllInteractionsStats=droplevels(subset(AllInteractions,treatment!='N'))

##Create mean and 95% CI for network properties

MeansNet=data.frame(timewindow=character(),timewindownumeric=numeric(),treatment=character(),Metric=character(),mean=double(),meanboot=double(),low=double(),high=double(),stringsAsFactors=FALSE)
new_row=list(timewindow=character(),timewindownumeric=numeric(),treatment=character(),Metric=character(),mean=double(),meanboot=double(),low=double(),high=double())

for (tr in unique(EigenVectorCentrality$treatment)) {
  for (tw in unique(EigenVectorCentrality$timewindow)){
    SubSet1=subset(EigenVectorCentrality,treatment==tr&timewindow==tw)
    SubSet2=subset(AverageInteractions,treatment==tr&timewindow==tw)
    SubSet3=subset(NodeStrengthCentrality,treatment==tr&timewindow==tw)
    l=length(SubSet1)
    new_row[1]=as.numeric(tw)*6
    new_row[2]=as.character(unique(SubSet1$timewindow))
    new_row[3]=tr
    for (met in 1:3){
      if (met==1){
        new_row[4]="Eigen centrality"
        mat=as.matrix(SubSet1[,8])} else if (met==2) {
          new_row[4]="Strength"
          mat=as.matrix(SubSet3[,7])} else if (met==3){ ###originally SubSet1[,10]
            new_row[4]="Specialization"
            mat=as.matrix(SubSet2[,8])
          }
      new_row[5]=mean(mat) 
      bstrap <- c()
      for (i in 1:1000){ bstrap <- c(bstrap, mean(sample(mat,l,replace=T)))}
      new_row[6]=mean(bstrap)
      new_row[7]=quantile(bstrap,.025)
      new_row[8]=quantile(bstrap,.975)
      MeansNet[nrow(MeansNet)+1,]=new_row
    }
  }
}

MeansNet$timewindow=as.factor(MeansNet$timewindow)
MeansNet$treatment=as.factor(MeansNet$treatment)
MeansNet$Metric=as.factor(MeansNet$Metric)
MeansNet$timewindownumeric=as.numeric(as.character(MeansNet$timewindow))

write.csv(MeansNet, "BootstrappedCI/TableS3EF_Behaviours.csv", row.names = FALSE)

EigenVectorCentrality$timewindow=as.numeric(EigenVectorCentrality$matrix)*6
EigenVectorCentrality$timewindownumeric=as.numeric(EigenVectorCentrality$timewindow)
EigenVectorCentrality$timewindow=as.factor(EigenVectorCentrality$timewindow)
EigenVectorCentrality$ID_colony=interaction(EigenVectorCentrality$tag,EigenVectorCentrality$colony)
NodeStrengthCentrality$timewindow=as.numeric(NodeStrengthCentrality$matrix)*6
NodeStrengthCentrality$timewindownumeric=as.numeric(NodeStrengthCentrality$timewindow)
NodeStrengthCentrality$timewindow=as.factor(NodeStrengthCentrality$timewindow)
AverageInteractions$hour=AverageInteractions$timewindow
AverageInteractions$timewindow=as.numeric(AverageInteractions$hour)*6
AverageInteractions$timewindownumeric=as.numeric(AverageInteractions$timewindow)
AverageInteractions$timewindow=as.factor(AverageInteractions$timewindow)

EigenVectorCentrality$treatment <- relevel(EigenVectorCentrality$treatment, "Z")
EigenVectorCentrality$treatment <- relevel(EigenVectorCentrality$treatment, "PBS")
EigenVectorCentrality$treatment <- relevel(EigenVectorCentrality$treatment, "N")
NodeStrengthCentrality$treatment <- relevel(NodeStrengthCentrality$treatment, "Z")
NodeStrengthCentrality$treatment <- relevel(NodeStrengthCentrality$treatment, "PBS")
NodeStrengthCentrality$treatment <- relevel(NodeStrengthCentrality$treatment, "N")
AverageInteractions$treatment <- relevel(AverageInteractions$treatment, "Z")
AverageInteractions$treatment <- relevel(AverageInteractions$treatment, "PBS")
AverageInteractions$treatment <- relevel(AverageInteractions$treatment, "N")

##EIGENVECTOR CENTRALITY Fig2F
#Statistical model
m5lme_a<-lme(eigen_centrality~ treatment * timewindow, random=~1|colony/ID_colony,data=CentralitiesStats,method = "ML")
drop1(m5lme_a,test="Chisq") 
m5lme_a<-lme(eigen_centrality~ treatment + timewindow, random=~1|colony/ID_colony,data=CentralitiesStats,method = "ML") #Model without interaction
drop1(m5lme_a,test="Chisq") 

#Validation model
#1) Outliers
boxplot(eigen_centrality ~ treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "centrality", data = CentralitiesStats)

#2) Normality residuals
E <- resid(m5lme_a)
hist(E)
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) #has to be <4

#3)Homogeneity of variance
#Plot residuals
op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
plot(m5lme_a, add.smooth = FALSE, which = 1)
hist(E, xlab = "Residuals", main = "")
plot(CentralitiesStats$eigen_centrality, E, xlab = "Centrality",
     ylab = "Residuals")
boxplot(E ~ CentralitiesStats$treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ CentralitiesStats$timewindow, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ CentralitiesStats$timewindow * CentralitiesStats$treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ CentralitiesStats$colony, main = "Colony")
abline(0, 0)
par(op)
par(mfrow = c(1, 1))

#Introduce variance structure
m5lme_b<-lme(eigen_centrality~ treatment + timewindow, random=~1|colony/ID_colony,weights=varIdent(form =~ 1|treatment*timewindow),
             control=list(msMaxIter=200,opt="optim",msVerbose=TRUE),data=CentralitiesStats,method = "ML")

E2 <- resid(m5lme_b, type = "normalized")
plot(m5lme_b, add.smooth = FALSE, which = 1)
coplot(E2 ~ eigen_centrality | treatment + timewindownumeric, data = CentralitiesStats,
       ylab = "Normalised residuals")

#Final statistics
AIC(m5lme_a,m5lme_b)
drop1(m5lme_b,test="Chisq")

##Node Strength Fig2G
m7lme_a<-lme(strength~ treatment * timewindow, random=~1|colony/ID_colony,data=CentralitiesStats,method = "ML")
drop1(m7lme_a,test="Chisq")
m7lme_a<-lme(strength~ treatment + timewindow, random=~1|colony/ID_colony,data=CentralitiesStats,method = "ML")
drop1(m7lme_a,test="Chisq")  

#1) Outliers
boxplot(strength ~ treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "strength", data = CentralitiesStats)

#2) Normality of residuals
E <- resid(m7lme_a)
hist(E)  
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) 
#There is an obvious residual outlier, I will compare the final model to one without outliers

#3) Homogeneity of variance
#Plot residuals
op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
plot(m7lme_a, add.smooth = FALSE, which = 1)
hist(E, xlab = "Residuals", main = "")
plot(CentralitiesStats$strength, E, xlab = "strength",
     ylab = "Residuals")
boxplot(E ~ CentralitiesStats$treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ CentralitiesStats$timewindow, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ CentralitiesStats$timewindow * CentralitiesStats$treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ CentralitiesStats$colony, main = "Colony")
abline(0, 0)
par(op)

#Introduce variance structure
m7lme_b<-lme(strength~ treatment + timewindow, random=~1|colony/ID_colony,weights=varIdent(form =~ 1|treatment*timewindow),
             control=list(msMaxIter=200,opt="optim",msVerbose=TRUE),data=CentralitiesStats,method = "ML")

E2 <- resid(m7lme_b, type = "normalized")
plot(m7lme_b, add.smooth = FALSE, which = 1)
coplot(E2 ~ strength | treatment * timewindownumeric, data = CentralitiesStats,
       ylab = "Normalised residuals")

#Without outlier
Strength_NoOutliers=subset(CentralitiesStats,resid(m7lme_a)>(-2)) #Remove outliers
m7lme_b_NoOut<-lme(strength~ treatment + timewindow, random=~1|colony/ID_colony,weights=varIdent(form =~ 1|treatment*timewindow),
                   control=list(msMaxIter=200,opt="optim",msVerbose=TRUE),data=Strength_NoOutliers,method = "ML")
ENO<- resid(m7lme_b_NoOut)
hist(ENO)
PerformanceAnalytics::skewness(ENO) #has to be -2< x <2
PerformanceAnalytics::kurtosis(ENO) #Has to be <4, without outliers residuals are normal

drop1(m7lme_b,test="Chisq")
drop1(m7lme_b_NoOut,test="Chisq") #Qualitatively equivalent models, w/ or without outliers, so I keep the outliers
AIC(m7lme_b,m7lme_b_NoOut)

##Skewness Fig.S3
#Statistical model
m6lme_a<-lme(specialization~ treatment * timewindow, random=~1|colony/ID_colony,data=AllInteractionsStats,method = "ML")
drop1(m6lme_a,test="Chisq") ##Drop interaction
m6lme_a<-lme(specialization~ treatment + timewindow, random=~1|colony/ID_colony,data=AllInteractionsStats,method = "ML")
drop1(m6lme_a,test="Chisq")  ##Treatment n.s.

#1) Outliers
boxplot(specialization ~ treatment,
        varwidth = TRUE, xlab = "Treatment",
        ylab = "PerformanceAnalytics::skewness", data = AllInteractionsStats)

#2) Normality residuals
E <- resid(m6lme_a)
hist(E)  #Long tail 
PerformanceAnalytics::skewness(E) #has to be -2< x <2
PerformanceAnalytics::kurtosis(E) #outlier, will compare later to model w/o outlier

#3)Homogeneity of variance
#Plot residuals

op <- par(mfrow = c(3, 2), mar = c(5, 4, 1, 2))
plot(m6lme_a, add.smooth = FALSE, which = 1)
hist(E, xlab = "Residuals", main = "")
plot(AllInteractionsStats$specialization, E, xlab = "PerformanceAnalytics::skewness",
     ylab = "Residuals")
boxplot(E ~ AllInteractionsStats$treatment, main = "treatment")
abline(0, 0)
boxplot(E ~ AllInteractionsStats$timewindow, main = "TimeWindow")
abline(0, 0)
boxplot(E ~ AllInteractionsStats$timewindow * AllInteractionsStats$treatment,
        main = "Interaction")
abline(0, 0)
boxplot(E ~ AllInteractionsStats$colony, main = "Colony")
abline(0, 0)
par(op)
par(mfrow = c(1, 1))

#Introduce variance structure
m6lme_b<-lme(sqrt(specialization)~ treatment + timewindow, random=~1|colony/ID_colony,data=AllInteractionsStats,weights=varIdent(form =~ 1|treatment*timewindow),
             control=list(msMaxIter=200,opt="optim",msVerbose=TRUE),method = "ML")

E2 <- resid(m6lme_b, type = "normalized")
plot(m6lme_b, add.smooth = FALSE, which = 1) #more spread, less conical
coplot(E2 ~ specialization | treatment * timewindow, data = AllInteractionsStats,
       ylab = "Normalised residuals")

#Without outlier
AllInt_NoOutliers=subset(AllInteractionsStats,resid(m6lme_a)<0.2) #Remove outliers
m6lme_NoOut<-lme(sqrt(specialization)~ treatment + timewindow, random=~1|colony/ID_colony,data=AllInt_NoOutliers,weights=varIdent(form =~ 1|treatment*timewindow),
                 control=list(msMaxIter=200,opt="optim",msVerbose=TRUE),method = "ML")
ENO<- resid(m6lme_NoOut)
hist(ENO)
PerformanceAnalytics::skewness(ENO) #has to be -2< x <2
PerformanceAnalytics::kurtosis(ENO) #Has to be <4, without outliers residuals are normal

drop1(m6lme_b,test="Chisq")
drop1(m6lme_NoOut,test="Chisq") #Qualitatively equivalent models, w/ or without outliers, so I keep the outliers
AIC(m6lme_a,m6lme_b)

##Sum of nodes strenght
SumStrength=data.frame(colony=character(),timewindow=double(),hour=double(),ID_colony=character(),tag=character(),treatment=character(),ant1=double(),ant2=double(),ant3=double(),ant4=double(),ant5=double(),ant6=double(),ant7=double(),ant8=double(),ID_int=character(),Max_int=double(),TOTint=double(),meanint=double(),specialization=double(),stringsAsFactors=FALSE)
IDs=names(UberMatList[[1]][[1]])
new_row=list(colony=character(),timewindow=double(),hour=double(),ID_colony=character(),tag=character(),treatment=character(),ant1=double(),ant2=double(),ant3=double(),ant4=double(),ant5=double(),ant6=double(),ant7=double(),ant8=double(),ID_int=character(),Max_int=double(),TOTint=double(),meanint=double(),specialization=double())

for (col in 1:length(Colonies)) {
  name_col=unlist(names(UberMatList[col]))
  new_row[1]=as.character(name_col) # Colony
  Matriceslist=UberMatList[[col]][] #[243:284] to take all matrices from hour 15 to 21 (centered around h18)
  for (mat in 1:length(Matriceslist)){
    name_mat=unlist(strsplit(names(Matriceslist[mat]),"_"))
    new_row[2]=name_mat[4] # matrix
    new_row[3]=name_mat[6] # hour
    for (ant1 in c(1:9)){
      ID=names(Matriceslist[[mat]][ant1])
      new_row[4]=as.character(paste(ID,name_col,sep="_")) # ID_colony
      new_row[5]=as.character(ID) # tag
      if (ID=="BB"){new_row[6]="Z"} else if (ID=="RR"){new_row[6]="PBS"} else {new_row[6]="N"} # treatment
      interactions=Matriceslist[[mat]][ant1,]
      interactions=interactions[-ant1]
      new_row[7:14]=interactions
      new_row[15]=names(which.max(interactions)) #ID_int
      #new_row[16]=max(interactions)/mean(as.double(interactions))
      new_row[16]=max(interactions) # Max_int units should be seconds
      new_row[17]=sum(interactions) # TOTint (YU: what are the units and are they consistent with the above analyses?)
      new_row[18]=mean(as.double(interactions)) # meanint
      new_row[19]=(max(interactions)/sum(interactions))#
      SumStrength[nrow(SumStrength)+1,]=new_row
    }
  }
}

SumStrength=ddply(SumStrength, .(colony, hour,timewindow), summarize, TOTint=sum(TOTint))
SumStrength$hour=as.double(SumStrength$hour)
SumStrength$colony=as.factor(SumStrength$colony)
SumStrength$timewindow=as.factor(SumStrength$timewindow)
SumStrength$logStr=log(SumStrength$TOTint)
Sumstrength.plot<-ggplot(aes(x=hour, y=logStr, color=colony),data=SumStrength)+
  geom_point()+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line()+
  ylab("log(Strength)")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),#text=element_text(size=10,family="Helvetica"),
        axis.text = element_text(size=9),legend.position="none")

ggsave(filename="Plots/FigS3_Sumstrengths.pdf",print(Sumstrength.plot),width=7,height=4,unit="in")

m6lme_a<-lme(logStr~ timewindow, random=~1|colony, data=SumStrength,method = "ML")
drop1(m6lme_a,test="Chisq") 

#Plots
eigen.plot <- 
  subset(MeansNet,Metric=="Eigen centrality") %>% 
  ggplot(aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=EigenVectorCentrality,aes(x=timewindownumeric,y=eigen_centrality,color=treatment, shape=treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='N'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='PBS'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='Z'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='N'))+
  geom_line(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='PBS'))+
  geom_line(data=subset(MeansNet,Metric=="Eigen centrality"&treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=0.5)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2,position=position_dodge(0.9)) +
  ylab("Eigenvector Centrality")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
         panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=10,family="Helvetica"),
        axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle"))  

strength.plot <- 
  subset(MeansNet,Metric=="Strength") %>% 
  ggplot(aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=NodeStrengthCentrality,aes(x=timewindownumeric,y=strength,color=treatment, shape=treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansNet,Metric=="Strength"&treatment=='N'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Strength"&treatment=='PBS'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Strength"&treatment=='Z'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansNet,Metric=="Strength"&treatment=='N'))+
  geom_line(data=subset(MeansNet,Metric=="Strength"&treatment=='PBS'))+
  geom_line(data=subset(MeansNet,Metric=="Strength"&treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=0.5)+
  ylab("log(strength)")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=10,family="Helvetica"),
        axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle")) 

spec.plot <- 
  subset(MeansNet,Metric=="Specialization") %>% 
  ggplot(aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=AverageInteractions,aes(x=timewindownumeric,y=specialization,color=treatment, shape=treatment),position=position_jitterdodge(jitter.width=6),size=0.6,alpha=0.5) +
  geom_point(data=subset(MeansNet,Metric=="Specialization"&treatment=='N'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Specialization"&treatment=='PBS'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  geom_point(data=subset(MeansNet,Metric=="Specialization"&treatment=='Z'),size=1,aes(x=timewindownumeric, y=mean, color=treatment, shape=treatment))+
  scale_x_continuous(breaks = c(6,12,18,24,30,36,42,48,54),labels=c("6-9","9-15","15-21","21-27","27-33","33-39","39-45","45-51","51-54") ) +
  geom_line(data=subset(MeansNet,Metric=="Specialization"&treatment=='N'))+
  geom_line(data=subset(MeansNet,Metric=="Specialization"&treatment=='PBS'))+
  geom_line(data=subset(MeansNet,Metric=="Specialization"&treatment=='Z'))+
  geom_errorbar(aes(ymin=low, ymax=high, width=.5), size=1)+
  ylab("Skewness")+xlab("Time post-treatment (hours)")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white",colour="white"), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),text=element_text(size=10,family="Helvetica"),
        axis.text.x = element_text(size=9,angle=45),axis.text.y = element_text(size=9),legend.position="none")+
  scale_color_manual(values=c(colNV, colSh,colTr))+
  scale_shape_manual(values=c("square","triangle","circle")) 

ggsave(filename="Plots/Fig2FG_Behaviours.pdf",plot=grid.arrange(eigen.plot,strength.plot,nrow=1),width=6.215,height=3,unit="in")
ggsave(filename="Plots/FigS4_Skewness.pdf",plot=spec.plot, width=3.125, height=3,unit="in")

ggsave(filename="Plots/Fig2.2_AllBehaviours.pdf",plot=grid.arrange(alone.plot,act.plot,speed.plot,eigen.plot,strength.plot,nrow=2),width=7.25,height=6,unit="in")

###FigS5###CHC
immubi<-read.csv(file="FigS5_CHC/ImmuneChallengeCHC.csv",row.names = 1)
imubi4<-immubi[,-c(1,2,3,4)]/100
imubi4

mds1 <- metaMDS(t(imubi4), distance="mahalanobis", k=2, trymax=20, autotransform=TRUE, noshare=0.1, expand=TRUE, trace=1)

# in ggplot
## create table for NMDS X and Y values of chemical compounds
data.scores<-as.data.frame(scores(mds1))
data.scores$site<-rownames(data.scores)
data.scores

## create table for NMDS X and Y values of individuals
species.scores<-as.data.frame(scores(mds1,"species"))
species.scores$species<-rownames(species.scores)
species.scores

## create groups of different treatments 
groups<-c(rep("NA", each=9),rep("NS",each=6),rep("PBSA",each=9),rep("PBSS", each=10),rep("ZA",each=8),rep("ZS",each=7))
Treatment<-as.character(groups)
data.scores$Treatment<-groups

##Adonis model
imubi5<-t(imubi4)
vegimubi<-vegdist(imubi5) ##n=15, pbs=19, z=15

#Adonis separating treatment from group
Group<-c(rep("Alone", each=9),rep("Social",each=6),rep("Alone",each=9),rep("Social", each=10),rep("Alone",each=8),rep("Social",each=7))
Experiment<-c(rep("N", each=9),rep("N",each=6),rep("PBS",each=9),rep("PBS", each=10),rep("Z",each=8),rep("Z",each=7))
ad1<-adonis(vegimubi~Group*Experiment)
ad1

dapc1 <- dapc(imubi5, groups)  #PC=8 DF=5
dapc1
scatter(dapc1,scree.pca=F,posi.pca="none") #FigS5
loadingplot(dapc1$var.contr,axis=1)
loadingplot(dapc1$var.contr,axis=2)

###Posthoc Adonis
pvalues=0

#Adonis: PBSA vs ZA
vegimubi4<-vegdist(imubi5[c(16:24,35:42),])
ad1<-adonis(vegimubi4~Treatment[c(16:24,35:42)])
ad1
pvalues[1]=ad1$aov.tab[1,6] #Extract pvalue

#Adonis: NA vs PBSA
vegimubi3<-vegdist(imubi5[c(1:9,16:24),])
ad1<-adonis(vegimubi3~Treatment[c(1:9,16:24)])
ad1
pvalues[2]=ad1$aov.tab[1,6]

#Adonis: NA vs ZA
vegimubi4<-vegdist(imubi5[c(1:9,35:42),])
ad1<-adonis(vegimubi4~Treatment[c(1:9,35:42)])
ad1
pvalues[3]=ad1$aov.tab[1,6]

p.adjust(pvalues,method="BH")

##Random forest analysis: What compound type is most important?###
imubi6<-data.frame(imubi5,groups)
imubi6$groups=as.factor(imubi6$groups)

#All molecules
fit <- randomForest(groups ~ .,
                    data=imubi6, 
                    importance=TRUE, 
                    ntree=2000)
varImpPlot(fit)

#For social/alone
imubi7<-data.frame(imubi5,Group)
imubi7$Group=as.factor(imubi7$Group)
fit <- randomForest(Group ~ .,
                    data=imubi7, 
                    importance=TRUE, 
                    ntree=2000)
varImpPlot(fit)

#For Treatment (N Z PBS)
imubi8<-data.frame(imubi5,Experiment)
imubi8$Experiment=as.factor(imubi8$Experiment)
fit <- randomForest(Experiment ~ .,
                    data=imubi8, 
                    importance=TRUE, 
                    ntree=2000)
varImpPlot(fit)