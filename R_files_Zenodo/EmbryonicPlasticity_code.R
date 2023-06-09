#This is code to replicate the analyses and figures from my 2023 Egg Priming paper. Code developed
#by Emily Harmon


library(ggplot2)
library(dplyr)
library(lme4)
library(lmerTest)
library(effects)
library(tidyr)
library(ggforce)

#Read in the data
tads <- read.csv("EmbryonicPlasticity.csv")
tads$Family <- factor(tads$Family)
tads$Round <- factor(tads$Round)
tads$Box <- factor(tads$Box)
tads$OH_L <- as.numeric(tads$OH_L)
tads$OH_R <- as.numeric(tads$OH_R)
tads$GutLength <- as.numeric(tads$GutLength)

#Does treatment impact SVL at hatching?_________________________________________

svlmod <- lmer(SVL~Treatment + Round +(1|Family/Box), data=tads)
summary(svlmod)#Round has no blocking effect

svlmod1 <-lmer(SVL~Treatment + (1|Family/Box), data=tads)
summary(svlmod1)
plot(allEffects(svlmod1))
drop1(svlmod1, ddf="Kenward-Roger")

#SVL figure
svl.summary<-group_by(tads,Treatment) %>%
  dplyr::summarise(
    sd=sd(SVL),
    se=sd/sqrt(n()),
    SVL=mean(SVL))
ggplot(tads, aes(x=Treatment, y=SVL))+
  geom_sina(aes(color=Treatment), maxwidth=.7, alpha=.5, show.legend = F)+
  geom_errorbar(aes(ymin=SVL-se, ymax=SVL+se), data=svl.summary, width=.2, color="black", size=1)+
  geom_point(data=svl.summary, size=4, shape=18)+
  theme_classic(base_size=20)+
  ylab(NULL)+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL)+
  scale_color_manual(values=c("#0C7BDC", "#FFA60A"))


#Does treatment impact Gosner stage at hatching?________________________________

Gosmod <- lmer(log(Gosner)~Treatment + Round + (1|Family/Box), data=tads)
summary(Gosmod) #Round has no blocking effect

Gosmod1 <- lmer(Gosner~Treatment + (1|Family/Box), data=tads)
summary(Gosmod1)
plot(allEffects(Gosmod1))
drop1(Gosmod1, ddf="Kenward-Roger")

#Gosner stage figure
gos.summary<-group_by(tads,Treatment) %>%
  dplyr::summarise(
    sd=sd(Gosner),
    se=sd/sqrt(n()),
    Gosner=mean(Gosner))
ggplot(tads, aes(x=Treatment, y=Gosner))+
  geom_sina(aes(color=Treatment), maxwidth=.7, alpha=.5, show.legend=F)+
  geom_errorbar(aes(ymin=Gosner-se, ymax=Gosner+se), data=gos.summary, width=.2, color="black", size=1)+
  geom_point(data=gos.summary, size=4, shape=18)+
  theme_classic(base_size=20)+
  ylab(NULL)+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL)+
  scale_color_manual(values=c("#0C7BDC", "#FFA60A"))


#Does treatment impact jaw width at hatching?___________________________________

#New dataset with the size-corrected OH jaw measurement
tads <- mutate(tads, OH.avg = (OH_L+OH_R)/2)
OHtads <- tads %>% drop_na(OH.avg)
OHtads <- OHtads%>%
  mutate(OH.SVL= resid(lm(log(OH.avg)~log(SVL))))

jawmod <- lmer(OH.SVL~Treatment + Round+(1|Family/Box), data=OHtads)
summary(jawmod)#Round has no blocking effect

jawmod1 <-lmer(OH.SVL ~ Treatment + (1|Family/Box), data=OHtads)
summary(jawmod1)
plot(allEffects(jawmod1))
drop1(jawmod1, ddf="Kenward-Roger")

#Jaw width figure
jaw.summary<-group_by(OHtads,Treatment) %>%
  dplyr::summarise(
    sd=sd(OH.SVL),
    se=sd/sqrt(n()),
    OH.SVL=mean(OH.SVL))
ggplot(OHtads, aes(x=Treatment, y=OH.SVL))+
  geom_sina(aes(color=Treatment), maxwidth=.7, alpha=.5, show.legend=F)+
  geom_errorbar(aes(ymin=OH.SVL-se, ymax=OH.SVL+se), data=jaw.summary, width=.2, color="black", size=1)+
  geom_point(data=jaw.summary, size=4, shape=18)+
  theme_classic(base_size=20)+
  ylab(NULL)+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL)+
  scale_color_manual(values=c("#0C7BDC", "#FFA60A"))


#Does treatment impact gut length at hatching?__________________________________

#New dataset with the size-corrected gut length
gltads <- tads %>% drop_na(GutLength)
gltads <- gltads%>%
  mutate(Gut.SVL= resid(lm(log(GutLength)~log(SVL))))

gutmod <- lmer(Gut.SVL~Treatment + Round + (1|Family/Box), data=gltads)
summary(gutmod)#Round has no blocking effect

gutmod1 <- lmer(Gut.SVL~Treatment + (1|Family/Box), data=gltads)
summary(gutmod1)
plot(allEffects(gutmod1))
drop1(gutmod1, ddf="Kenward-Roger")

#Gut length figure
gut.summary<-group_by(gltads,Treatment) %>%
  dplyr::summarise(
    sd=sd(Gut.SVL),
    se=sd/sqrt(n()),
    Gut.SVL=mean(Gut.SVL))
ggplot(gltads, aes(x=Treatment, y=Gut.SVL))+
  geom_sina(aes(color=Treatment),maxwidth=.7, alpha=.5, show.legend=F)+
  geom_errorbar(aes(ymin=Gut.SVL-se, ymax=Gut.SVL+se), data=gut.summary, width=.2, color="black", size=1)+
  geom_point(data=gut.summary, size=4, shape=18)+
  theme_classic(base_size=20)+
  ylab(NULL)+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = NULL)+
  scale_color_manual(values=c("#0C7BDC", "#FFA60A"))
