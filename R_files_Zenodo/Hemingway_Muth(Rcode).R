#COLOR CONTRAST EXPERIMENT

#LOAD PACKAGES-----

library(tidyverse)
library(ggpubr)
library(boot)
library(lme4)
library(effects)

#LOAD DATASETS-----

#same/switch data (experiment 1)
ssdata <- read.csv("BBeeDataSameSwitch.csv", stringsAsFactors = FALSE)
ssdata

#gradient data (experiment 2)
gdata <- read.csv("GradientExp.csv", stringsAsFactors = FALSE)
gdata


#BOOTSTRAPPING CODE----

#function to get mean and 95% CI of values x via bootstrapping (BCa method)
boot_ci <- function(x, perms=5000, bca = F) {
  get_mean <- function(x, d) {
    return(mean(x[d]))
  }
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  if(bca){
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="bca")
    low <- boot$bca[1,4]
    high <- boot$bca[1,5] 
  }else{
    boot <- boot.ci(boot(data=x, 
                         statistic=get_mean, 
                         R=perms, 
                         parallel = "multicore", 
                         ncpus = 4), 
                    type="perc")
    low <- boot$perc[1,4]
    high <- boot$perc[1,5] 
  }
  c(low=low,mean=mean,high=high, N=round(length(x)))
}

#function to get mean
mean.w <- function(x,w) sum(x*w)

#function to get lower and upper confidence intervals
boot_ci_low <- function(x) {
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  low <- boot.ci(boot(data=x, statistic=mean.w, R=5000, stype="w"), type="bca")$bca[1,4]
  low
}

boot_ci_high <- function(x) {
  x <- as.vector(na.omit(x))
  mean <- mean(x)
  high <- boot.ci(boot(data=x, statistic=mean.w, R=5000, stype="w"), type="bca")$bca[1,5]
  high
}

#ANALYSES----

#Experiment 1----

#trim data into first 20 visits 
ss.trim.20<-
  ssdata%>%
  filter(Number %in% c(1:20))

ss.trim.20

#trim data into first 5 visits 
ss.trim.5<-
  ssdata%>%
  filter(Number %in% c(1:5))

ss.trim.5

#null model
SS0<-glmer(Acceptance ~ 
             + (1|Bee), family=binomial, data=ss.trim.20)

#test for effect of treatment
SS1 <- glmer(Acceptance ~ 
               SS 
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS1)
plot(allEffects(SS1))
anova(SS0,SS1)

#test for effect of treatment and number
SS2 <- glmer(Acceptance ~ 
               SS + Number
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS2)
plot(allEffects(SS2))
anova(SS1,SS2)

#test for interaction between treatment and number
SS3 <- glmer(Acceptance ~ 
               SS * Number
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS3)
plot(allEffects(SS3))
anova(SS2,SS3)

#look for interaction between stimuli and number
SS4 <- glmer(Acceptance ~ 
               Stimuli + Number 
             + Stimuli * Number
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS4)
plot(allEffects(SS4))
anova(SS3,SS4)

#include training color into model
SS5 <- glmer(Acceptance ~ 
               SS + Number + TrainColor
             + SS*Number*TrainColor
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS5)
plot(allEffects(SS5))

#effect of only training color 
SS6 <- glmer(Acceptance ~ 
               TrainColor
             + (1|Bee), family=binomial, data=ss.trim.20)

summary(SS6)
plot(allEffects(SS6))

#rerun models with only 5 visits

#null model
SS100<-glmer(Acceptance ~ 
             + (1|Bee), family=binomial, data=ss.trim.5)

#test for effect of treatment
SS101 <- glmer(Acceptance ~ 
               SS 
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS101)
plot(allEffects(SS101))
anova(SS100,SS101)

#test for effect of treatment and number
SS102 <- glmer(Acceptance ~ 
               SS + Number
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS102)
plot(allEffects(SS102))
anova(SS101,SS102)

#test for interaction between treatment and number
SS103 <- glmer(Acceptance ~ 
             SS * Number
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS103)
plot(allEffects(SS103))
anova(SS102,SS103)

#look for interaction between stimuli and number
SS104 <- glmer(Acceptance ~ 
               Stimuli + Number 
             + Stimuli * Number
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS104)
plot(allEffects(SS104))
anova(SS103,SS104)

#include training color into model
SS105 <- glmer(Acceptance ~ 
               SS + Number + TrainColor
             + SS * Number * TrainColor
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS105)
plot(allEffects(SS105))

#effect of only training color 
SS106 <- glmer(Acceptance ~ 
               TrainColor
             + (1|Bee), family=binomial, data=ss.trim.5)

summary(SS106)
plot(allEffects(SS106))

#just look at the first visit

#first visit behavior
ss.visit.1<-
  ssdata%>%
  filter(Number == 1)

v1<-glm(Acceptance~ SS, family=binomial, data=ss.visit.1)
summary(v1)

#chi-square to see if there is an effect of treatment or stimuli on acceptance behavior
table(ss.v1$SS)
table(ss.v1$Acceptance)
table(ss.v1$SS, ss.v1$Acceptance)
prop.table(table(ss.v1$SS, ss.v1$Acceptance))
prop.table(table(ss.v1$SS, ss.v1$Acceptance), margin=2)

test1<-chisq.test(ss.v1$SS, ss.v1$Acceptance)
test1

test2<-fisher.test(ss.v1$SS, ss.v1$Acceptance)
test2

#Number to first accept
ssnfa<-
  ss.trim.20%>%
  group_by(Bee, SS)%>%
  summarize(First.Accept=mean(First.Accept))

ssnfa

mean(ssnfa$First.Accept)
var(ssnfa$First.Accept)

ssnfa101<-glm(First.Accept ~ SS, family=quasipoisson(), data=ssnfa)
summary(ssnfa101)

#Experiment 2----

#trim data into first 20 visits

g.trim.20<-
  gdata%>%
  filter(Number %in% c(1:20))

g.trim.20

#trim data into first 5 visits

g.trim.5<-
  gdata%>%
  filter(Number %in% c(1:5))

g.trim.5

#Models

#null model

g0<- glmer(Acceptance ~ 
             + (1|Bee) + (1|Colony), family=binomial, data=gdata)

#is there an effect of treatment on preference?

g1<- glmer(Acceptance ~ 
             Treatment 
           + (1|Bee) + (1|Colony), family=binomial, data=gdata)

summary(g1)
plot(allEffects(g1))
anova(g0,g1)

#is there an effect of visit number?
g2<-glmer(Acceptance ~ 
            Number 
          + (1|Bee) + (1|Colony), family=binomial, data=gdata)

summary(g2)
plot(allEffects(g2))
anova(g0,g2)

#What happens when we add stimulus back in?
g3<-glmer(Acceptance ~ 
            Treatment + Number 
           + (1|Bee) + (1|Colony), family=binomial, data=gdata)

summary(g3)
plot(allEffects(g3))
anova(g2,g3)

#Look for interaction

g4<-glmer(Acceptance ~ 
            Treatment * Number 
          + (1|Bee) + (1|Colony), family=binomial, data=gdata)

summary(g4)
plot(allEffects(g4))
anova(g2,g4)

#So there is an effect of number but not an interaction between stimuli and number

#rerun these models with just the first 5 visits

#null model

g100<- glmer(Acceptance ~ 
             + (1|Bee) + (1|Colony), family=binomial, data=g.trim.5)

#is there an effect of treatment on preference?

g101<- glmer(Acceptance ~ 
             Treatment 
           + (1|Bee) + (1|Colony), family=binomial, data=g.trim.5)

summary(g101)
plot(allEffects(g101))
anova(g100,g101)

#is there an effect of visit number?
g102<-glmer(Acceptance ~ 
            Number 
          + (1|Bee) + (1|Colony), family=binomial, data=g.trim.5)

summary(g102)
plot(allEffects(g102))
anova(g100,g102)

#What happens when we add stimulus back in?
g103<-glmer(Acceptance ~ 
            Treatment + Number 
          + (1|Bee) + (1|Colony), family=binomial, data=g.trim.5)

summary(g103)
plot(allEffects(g103))
anova(g102,g103)

#Look for interaction

g104<-glmer(Acceptance ~ 
            Treatment * Number 
          + (1|Bee) + (1|Colony), family=binomial, data=g.trim.5)

summary(g104)
plot(allEffects(g104))
anova(g102,g104)


#look at every pairwise comparison
modelab<-glmer(Acceptance ~ 
                 Treatment * Number + 
                 (1|Bee) + (1|Colony), family=binomial, 
               data=gdata[gdata$Treatment%in%c("1","2"),])
summary(modelab)

modelac<-glmer(Acceptance ~ 
                 Treatment * Number + 
                 (1|Bee) + (1|Colony), family=binomial, 
               data=gdata[gdata$Treatment%in%c("1","3"),])
summary(modelac)

modelad<-glmer(Acceptance ~ 
                 Treatment * Number + 
                 (1|Bee) + (1|Colony), family=binomial, 
               data=gdata[gdata$Treatment%in%c("1","4"),])
summary(modelad)

modelae<-glmer(Acceptance ~ 
                 Treatment * Number + 
                 (1|Bee) + (1|Colony), family=binomial, 
               data=gdata[gdata$Treatment%in%c("1","5"),])
summary(modelae)

#number to first acceptance
gnfa<-
  g.trim.20%>%
  group_by(Bee, Treatment)%>%
  summarize(First.Accept=mean(First.Accept))

gnfa

mean(gnfa$First.Accept)
var(gnfa$First.Accept)

gnfa101<-glm(First.Accept ~ Treatment, family=quasipoisson(), data=gnfa)
summary(gnfa101)

#run a chi square test - difference in initial acceptance by test stimuli?
g.visit.1
table(g.visit.1$Treatment)
table(g.visit.1$Acceptance)
table(g.visit.1$Treatment, g.visit.1$Acceptance)
prop.table(table(g.visit.1$Treatment, g.visit.1$Acceptance))
prop.table(table(g.visit.1$Treatment, g.visit.1$Acceptance), margin=2)

test3<-chisq.test(g.visit.1$Treatment, g.visit.1$Acceptance)
test3

test4<-fisher.test(g.visit.1$Treatment, g.visit.1$Acceptance)
test4

v2<-glm(Acceptance~ Treatment, family=binomial, data=g.visit.1)
summary(v2)

#do we see a difference in number until acceptance depending on stimuli?

g7<-glmer(First.Accept ~ Stimuli, family=poisson(), data=g.trim.20)

summary(g7)
plot(allEffects(g8))


#Just pull our blue and yellow
blue_yellow<- subset(g.trim.20,g.trim.20$Treatment %in% c(1,5))
blue_yellow

g8<-glmer(Acceptance ~ Stimuli + Number + (1|Bee)+ (1|Colony), 
          family=binomial, data=blue_yellow)

summary(g8)
plot(allEffects(g8))

#add in visit number
g9<-glmer(Acceptance ~ Stimuli + Number (1|Bee) +(1|Colony), 
          family=binomial, data=blue_yellow)

summary(g9)
plot(allEffects(g9))

#look for an interaction between number and stimulus 
g10<-glmer(Acceptance ~ 
             Stimuli * Number + 
             (1|Bee) + (1|Colony), family=binomial, data=blue_yellow)

summary(g10)
plot(allEffects(g10))

anova(g10,glm10)

#look at just the first visit between blue and yellow

g.visit.1.BY<-
  blue_yellow%>%
  filter(Number == 1)

g11<-glm(Acceptance ~ Stimuli, family=binomial, data=fg.visit.1.BY)

summary(g11)
plot(allEffects(g11))

#acceptance on first visit
table(g.visit.1$Treatment)
table(g.visit.1$Acceptance)
table(g.visit.1$Treatment, g.visit.1$Acceptance)
prop.table(table(g.visit.1$Treatment, g.visit.1$Acceptance))
prop.table(table(g.visit.1$Treatment, g.visit.1$Acceptance), margin=2)

test1<-chisq.test(g.visit.1$Treatment, g.visit.1$Acceptance)
test1

test2<-fisher.test(g.visit.1$Treatment, g.visit.1$Acceptance)
test2

#PLOTS----

#Experiment1----

#treatment plot for first 20
ss.treatment.20<-
  ss.trim.20%>%
  group_by(Number,SS)%>%
  summarize(n=n(),Acceptance=sum(Acceptance))%>%
  ungroup()%>%
  mutate(Proportion=Acceptance/n)%>%
  ggplot(aes(x=Number, y=Proportion, color=SS))+
  scale_color_manual(values = c("#0a090a","#0a090a"))+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(expand = c(0,0),limits=c(0,20), breaks = c(seq(0,20,5)))+
  coord_cartesian(xlim = c(0,21), ylim=c(0,1))+
  theme_pubr()+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~ SS, ncol=2)+
  labs(y = "Acceptance", x = "Visit Number")+
  theme(aspect.ratio = 2/2)

ss.treatment.20

ggsave(plot = last_plot(), filename = "sstimuliplot.jpg", units = "cm", width = 17, height = 16)


citation()

#stimulus plot for first 20 
ss.stimuli.20<-
  ss.trim.20%>%
  group_by(Number,Stimuli)%>%
  summarize(n=n(),Acceptance=sum(Acceptance))%>%
  ungroup()%>%
  mutate(Proportion=Acceptance/n)%>%
  ggplot(aes(x=Number, y=Proportion, color=Stimuli))+
  scale_color_manual(values = c("#3791ad","#edd03b","#3791ad","#edd03b"))+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(expand = c(0,0),limits=c(0,20), breaks = c(seq(0,20,5)))+
  coord_cartesian(xlim = c(0,21), ylim=c(0,1))+
  theme_pubr()+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~ Stimuli, ncol=2)+
  labs(y = "Acceptance")+
  labs(x = "Visit number")

ss.stimuli.20

#number to first acceptance by treatment
ss.1a.treatment<-
  ss.trim.20%>%
  filter(!First.Accept>20)%>%
  group_by(SS, Bee)%>%
  summarise(accept = mean(First.Accept))%>%
  mutate(mean.overall=mean(accept), upper=boot_ci_high(accept), lower=boot_ci_low(accept))%>%
  ggplot(aes(x = SS, y = accept))+geom_jitter(width = 0.2, height = 0.01, alpha=0.7)+
  geom_point(aes(x = SS, y = mean.overall), color = "black")+
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, color = "black")+
  theme_pubr()+
  labs(y = "Visits to first acceptance", x = "Treatment type")+theme(aspect.ratio = 1.2/1)+
  scale_y_continuous(limits = c(0,21), breaks = c(seq(0,20,5)))

ss.1a.treatment

ss.1a.treatment[["data"]]%>%
  group_by(SS)%>%
  summarise(mean = mean(mean.overall), upper = mean(upper), lower = mean(lower))

#plot the chi-square distribution for initial acceptance by treatment
ss.visit.one<-
  ggplot(ss.visit.1)+
  geom_bar(aes(x=factor(SS), fill = factor(Acceptance)), width = 0.5)+
  theme_pubr()+
  scale_fill_grey(labels = c("Reject", "Accept"))+
  labs(y = "Frequency", x="Treatment", fill = "Acceptance")+
  theme(aspect.ratio = 1.2/1)+theme(legend.title = element_blank())

ss.visit.one

#Experiment 2----

#stimulus plot first 20
g.stimuli.20<-
  g.trim.20%>%
  group_by(Number,TestColor)%>%
  summarize(n=n(),Acceptance=sum(Acceptance))%>%
  ungroup()%>%
  mutate(Proportion=Acceptance/n)%>%
  ggplot(aes(x=Number, y=Proportion, color=TestColor))+
  scale_color_manual(values = c("#3791ad","#37ada3","#37ad82","#37ad5e","#edd03b"))+
  scale_y_continuous(expand = c (0,0),limits=c(-1,1))+
  scale_x_continuous(expand = c(0,0),limits=c(0,20), breaks = c(seq(0,20,5)))+
  coord_cartesian(xlim = c(0,21), ylim=c(0,1))+
  theme_pubr()+
  geom_smooth(method="lm")+
  facet_wrap(~ TestColor, ncol=2)+
  geom_point(size=1.5)+
  geom_smooth(method="lm")+
  facet_wrap(~ TestColor, ncol=5)+
  theme(aspect.ratio=2/2)+
  labs(y = "Acceptance")+
  labs(x = "Visit number")+
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size=13), axis.title.x = element_text(size=13))
  
g.stimuli.20
ggsave(plot = last_plot(), filename = "gstimuliplot.jpg", units = "cm", width = 17, height = 10)



#stimulus plot for blue stimuli first 20 visits
bluestim<- subset(g.trim.20,g.trim.20$Treatment %in% c(1,2,3,4))
bluestim

g.bluestim.20<-
  bluestim%>%
  group_by(Number,Stimuli)%>%
  summarize(n=n(),Acceptance=sum(Acceptance))%>%
  ungroup()%>%
  mutate(Proportion=Acceptance/n)%>%
  ggplot(aes(x=Number, y=Proportion, color=Stimuli))+
  scale_color_manual(values = c("#3791ad","#37ada3","#37ad82","#37ad5e"))+
  scale_y_continuous(expand = c (0,0),limits=c(-1,1))+
  scale_x_continuous(expand = c(0,0),limits=c(0,20), breaks = c(seq(0,20,5)))+
  coord_cartesian(xlim = c(0,21), ylim=c(0,1))+
  theme_pubclean()+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~ Stimuli, ncol=2)+
  labs(x = "Visit number", y = "Acceptance")+
  mfrow=c(1, 4)

g.bluestim.20

#BY stimulus plot first 20

g.stimuli.BY.20<-
  blue_yellow%>%
  filter(Number %in% c(1:20))%>%
  group_by(Number,Stimuli)%>%
  summarize(n=n(),Acceptance=sum(Acceptance))%>%
  ungroup()%>%
  mutate(Proportion=Acceptance/n)%>%
  ggplot(aes(x=Number, y=Proportion, color=Stimuli))+
  scale_color_manual(values = c("#3791ad","#edd03b"))+
  scale_y_continuous(expand = c(0,0),limits=c(-1,1))+
  scale_x_continuous(expand = c(0,0),limits=c(0,20), breaks = c(seq(0,20,5)))+
  coord_cartesian(xlim = c(0,21), ylim=c(0,1))+
  theme_pubr()+
  geom_point()+
  geom_smooth(method="lm", fullrange = T)+
  facet_wrap(~ Stimuli, ncol=2)+
  labs(x = "Visit number", y = "Proportion accepted")

g.stimuli.BY.20

#number to first accept by stimulus
g.1a.stimuli<-
  g.trim.20%>%
  filter(!First.Accept>20)%>%
  group_by(TestColor, Bee)%>%
  summarise(accept = mean(First.Accept))%>%
  mutate(mean.overall=mean(accept), upper=boot_ci_high(accept), lower=boot_ci_low(accept))%>%
  ggplot(aes(x = TestColor, y = accept, color = TestColor))+geom_jitter(width = 0.2, height = 0.01, alpha=0.7)+
  geom_point(aes(x = TestColor, y = mean.overall, color = TestColor))+
  geom_errorbar(aes(ymin = lower, ymax = upper, color = TestColor), width = 0.1)+
  scale_color_manual(values = c("#3791ad","#37ada3","#37ad82","#37ad5e","#edd03b"))+
  theme_pubr()+
  labs(y = "Visits to first acceptance", x = "Test color")+theme(aspect.ratio = 1.2/1)+
  scale_y_continuous(limits = c(0,21), breaks = c(seq(0,20,5)))

g.1a.stimuli

# get confidence intervals
g.1a.stimuli[["data"]]%>%
  group_by(TestColor)%>%
  summarise(mean = mean(mean.overall), upper = mean(upper), lower = mean(lower))


#plot the chi-square distribution for initial acceptance by treatment
g.visit.one<-
  ggplot(g.visit.1)+
  geom_bar(aes(x=factor(TestColor), fill = factor(Acceptance)), width = 0.75)+
  theme_pubr()+
  scale_fill_grey(labels = c("Reject", "Accept"))+
  labs(y = "Frequency", x="Test color", fill = "Acceptance")+
  theme(aspect.ratio = 1.2/1)+theme(legend.title = element_blank())

g.visit.one



