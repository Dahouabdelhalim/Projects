## R code to produce Figure 3 in: 
## Johnston EC, Wyatt ASJ, Leichter JJ, Burgess SC. Niche differences in co-occuring cryptic coral species (Pocillopora spp). Coral Reefs.
## Code written by Erika Johnston. March 2021. Send comments or corrections to ejohnston@bio.fsu.edu
## RStudio version 3.6.2 (2021-02-29)

rm(list=ls()) #clear all variables

### Load packages
library(ggplot2)
library(RColorBrewer)
library(AICcmodavg)
library(dplyr)
library(car)

### Import data
d <- read.csv("Figure 3 and 4 Data.csv")
d$Depth.m <- as.factor(d$Depth.m)

### Binomial generalized linear mixed models
#### Haplotype 1a_Pm
d$Haplotype1a_Pm <- ifelse(d$Species.haplotype=="Haplotype 1a_Pm",1,0)
with(d[d$Site=="1",],table(Haplotype1a_Pm,Depth.m))

## Site 1
d_1a_Pm_S1<-d %>% 
  filter((Site=="1"))

m_Pm_S1 <- glm(Haplotype1a_Pm ~ Depth.m, family="binomial", data=d_1a_Pm_S1)

Anova(m_Pm_S1)
summary(m_Pm_S1)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("1")) 

Pm_S1 <- m_Pm_S1 
pred_1a_Pm_S1 <- as.data.frame(cbind(newdat,predict.glm(Pm_S1,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pm_S1 <- cbind(Haplotype = "Haplotype 1a_Pm", pred_1a_Pm_S1)
pred_1a_Pm_S1$prop <- plogis(pred_1a_Pm_S1$fit)
pred_1a_Pm_S1$prop.lwr <- plogis(pred_1a_Pm_S1$fit - 2*pred_1a_Pm_S1$se.fit)
pred_1a_Pm_S1$prop.upr <- plogis(pred_1a_Pm_S1$fit + 2*pred_1a_Pm_S1$se.fit)

## Site 2
d_1a_Pm_S2<-d %>% 
  filter((Site=="2"))

m_Pm_S2 <- glm(Haplotype1a_Pm ~ Depth.m, family="binomial", data=d_1a_Pm_S2)

Anova(m_Pm_S2)
summary(m_Pm_S2)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("2")) 

Pm_S2 <- m_Pm_S2
pred_1a_Pm_S2 <- as.data.frame(cbind(newdat,predict.glm(Pm_S2,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pm_S2 <- cbind(Haplotype = "Haplotype 1a_Pm", pred_1a_Pm_S2)
pred_1a_Pm_S2$prop <- plogis(pred_1a_Pm_S2$fit)
pred_1a_Pm_S2$prop.lwr <- plogis(pred_1a_Pm_S2$fit - 2*pred_1a_Pm_S2$se.fit)
pred_1a_Pm_S2$prop.upr <- plogis(pred_1a_Pm_S2$fit + 2*pred_1a_Pm_S2$se.fit)

## Site 4
d_1a_Pm_S4<-d %>% 
  filter((Site=="4"),
         !(Depth.m=="21"))

m_Pm_S4 <- glm(Haplotype1a_Pm ~ Depth.m, family="binomial", data=d_1a_Pm_S4)

Anova(m_Pm_S4)
summary(m_Pm_S4)

newdat <- expand.grid(Depth.m=c("5","10"),Site=c("4")) 

Pm_S4 <- m_Pm_S4
pred_1a_Pm_S4 <- as.data.frame(cbind(newdat,predict.glm(Pm_S4,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pm_S4 <- cbind(Haplotype = "Haplotype 1a_Pm", pred_1a_Pm_S4)
pred_1a_Pm_S4$prop <- plogis(pred_1a_Pm_S4$fit)
pred_1a_Pm_S4$prop.lwr <- plogis(pred_1a_Pm_S4$fit - 2*pred_1a_Pm_S4$se.fit)
pred_1a_Pm_S4$prop.upr <- plogis(pred_1a_Pm_S4$fit + 2*pred_1a_Pm_S4$se.fit)

## Site 5
d_1a_Pm_S5<-d %>% 
  filter((Site=="5"))

m_Pm_S5 <- glm(Haplotype1a_Pm ~ Depth.m, family="binomial", data=d_1a_Pm_S5)

Anova(m_Pm_S5)
summary(m_Pm_S5)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("5")) 

Pm_S5 <- m_Pm_S5
pred_1a_Pm_S5 <- as.data.frame(cbind(newdat,predict.glm(Pm_S5,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pm_S5 <- cbind(Haplotype = "Haplotype 1a_Pm", pred_1a_Pm_S5)
pred_1a_Pm_S5$prop <- plogis(pred_1a_Pm_S5$fit)
pred_1a_Pm_S5$prop.lwr <- plogis(pred_1a_Pm_S5$fit - 2*pred_1a_Pm_S5$se.fit)
pred_1a_Pm_S5$prop.upr <- plogis(pred_1a_Pm_S5$fit + 2*pred_1a_Pm_S5$se.fit)

pred_Pm <- rbind(pred_1a_Pm_S1,pred_1a_Pm_S2,pred_1a_Pm_S4,pred_1a_Pm_S5)

#### Haplotype 8a
d$Haplotype_8a <- ifelse(d$Species.haplotype=="Haplotype 8a",1,0)
with(d[d$Site=="1",],table(Haplotype_8a,Depth.m))

## Site 1
d_8a_S1<-d %>% 
  filter((Site=="1"),
         !(Depth.m=="21"))

m_8a_S1 <- glm(Haplotype_8a ~ Depth.m, family="binomial", data=d_8a_S1)

Anova(m_8a_S1)
summary(m_8a_S1)

newdat <- expand.grid(Depth.m=c("5","10"),Site=c("1")) 

H8a_S1 <- m_8a_S1 
pred_8a_S1 <- as.data.frame(cbind(newdat,predict.glm(H8a_S1,newdata=newdat,type="link",se.fit=T)))
pred_8a_S1 <- cbind(Haplotype = "Haplotype 8a", pred_8a_S1)
pred_8a_S1$prop <- plogis(pred_8a_S1$fit)
pred_8a_S1$prop.lwr <- plogis(pred_8a_S1$fit - 2*pred_8a_S1$se.fit)
pred_8a_S1$prop.upr <- plogis(pred_8a_S1$fit + 2*pred_8a_S1$se.fit)

## Site 2
d_8a_S2<-d %>% 
  filter((Site=="2"))

m_8a_S2 <- glm(Haplotype_8a ~ Depth.m, family="binomial", data=d_8a_S2)

Anova(m_8a_S2)
summary(m_8a_S2)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("2")) 

H8a_S2 <- m_8a_S2
pred_8a_S2 <- as.data.frame(cbind(newdat,predict.glm(H8a_S2,newdata=newdat,type="link",se.fit=T)))
pred_8a_S2 <- cbind(Haplotype = "Haplotype 8a", pred_8a_S2)
pred_8a_S2$prop <- plogis(pred_8a_S2$fit)
pred_8a_S2$prop.lwr <- plogis(pred_8a_S2$fit - 2*pred_8a_S2$se.fit)
pred_8a_S2$prop.upr <- plogis(pred_8a_S2$fit + 2*pred_8a_S2$se.fit)

## Site 4
d_8a_S4<-d %>% 
  filter((Site=="4"),
         !(Depth.m=="21"))

m_8a_S4 <- glm(Haplotype_8a ~ Depth.m, family="binomial", data=d_8a_S4)

Anova(m_8a_S4)
summary(m_8a_S4)

newdat <- expand.grid(Depth.m=c("5","10"),Site=c("4")) 

H8a_S4 <- m_8a_S4
pred_8a_S4 <- as.data.frame(cbind(newdat,predict.glm(H8a_S4,newdata=newdat,type="link",se.fit=T)))
pred_8a_S4 <- cbind(Haplotype = "Haplotype 8a", pred_8a_S4)
pred_8a_S4$prop <- plogis(pred_8a_S4$fit)
pred_8a_S4$prop.lwr <- plogis(pred_8a_S4$fit - 2*pred_8a_S4$se.fit)
pred_8a_S4$prop.upr <- plogis(pred_8a_S4$fit + 2*pred_8a_S4$se.fit)

## Site 5
d_8a_S5<-d %>% 
  filter((Site=="5"))

m_8a_S5 <- glm(Haplotype_8a ~ Depth.m, family="binomial", data=d_8a_S5)

Anova(m_8a_S5)
summary(m_8a_S5)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("5")) 

H8a_S5 <- m_8a_S5
pred_8a_S5 <- as.data.frame(cbind(newdat,predict.glm(H8a_S5,newdata=newdat,type="link",se.fit=T)))
pred_8a_S5 <- cbind(Haplotype = "Haplotype 8a", pred_8a_S5)
pred_8a_S5$prop <- plogis(pred_8a_S5$fit)
pred_8a_S5$prop.lwr <- plogis(pred_8a_S5$fit - 2*pred_8a_S5$se.fit)
pred_8a_S5$prop.upr <- plogis(pred_8a_S5$fit + 2*pred_8a_S5$se.fit)

pred_8a <- rbind(pred_8a_S1,pred_8a_S2,pred_8a_S4,pred_8a_S5)

#### Haplotype 11
d$Haplotype_11 <- ifelse(d$Species.haplotype=="Haplotype 11",1,0)
with(d[d$Site=="1",],table(Haplotype_11,Depth.m))

## Site 1
d_11_S1<-d %>% 
  filter((Site=="1"))

m_11_S1 <- glm(Haplotype_11 ~ Depth.m, family="binomial", data=d_11_S1)

Anova(m_11_S1)
summary(m_11_S1)

newdat <- expand.grid(Depth.m=c("10"),Site=c("1")) 

H11_S1 <- m_11_S1 
pred_11_S1 <- as.data.frame(cbind(newdat,predict.glm(H11_S1,newdata=newdat,type="link",se.fit=T)))
pred_11_S1 <- cbind(Haplotype = "Haplotype 11", pred_11_S1)
pred_11_S1$prop <- plogis(pred_11_S1$fit)
pred_11_S1$prop.lwr <- plogis(pred_11_S1$fit - 2*pred_11_S1$se.fit)
pred_11_S1$prop.upr <- plogis(pred_11_S1$fit + 2*pred_11_S1$se.fit)

## Site 2
d_11_S2<-d %>% 
  filter((Site=="2"))

m_11_S2 <- glm(Haplotype_11 ~ Depth.m, family="binomial", data=d_11_S2)

Anova(m_11_S2)
summary(m_11_S2)

newdat <- expand.grid(Depth.m=c("10"),Site=c("2")) 

H11_S2 <- m_11_S2
pred_11_S2 <- as.data.frame(cbind(newdat,predict.glm(H11_S2,newdata=newdat,type="link",se.fit=T)))
pred_11_S2 <- cbind(Haplotype = "Haplotype 11", pred_11_S2)
pred_11_S2$prop <- plogis(pred_11_S2$fit)
pred_11_S2$prop.lwr <- plogis(pred_11_S2$fit - 2*pred_11_S2$se.fit)
pred_11_S2$prop.upr <- plogis(pred_11_S2$fit + 2*pred_11_S2$se.fit)

## Site 4
d_11_S4<-d %>% 
  filter((Site=="4"))

m_11_S4 <- glm(Haplotype_11 ~ Depth.m, family="binomial", data=d_11_S4)

Anova(m_11_S4)
summary(m_11_S4)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("4")) 

H11_S4 <- m_11_S4
pred_11_S4 <- as.data.frame(cbind(newdat,predict.glm(H11_S4,newdata=newdat,type="link",se.fit=T)))
pred_11_S4 <- cbind(Haplotype = "Haplotype 11", pred_11_S4)
pred_11_S4$prop <- plogis(pred_11_S4$fit)
pred_11_S4$prop.lwr <- plogis(pred_11_S4$fit - 2*pred_11_S4$se.fit)
pred_11_S4$prop.upr <- plogis(pred_11_S4$fit + 2*pred_11_S4$se.fit)

## Site 5
d_11_S5<-d %>% 
  filter((Site=="5"),
         !(Depth.m=="5"))

m_11_S5 <- glm(Haplotype_11 ~ Depth.m, family="binomial", data=d_11_S5)

Anova(m_11_S5)
summary(m_11_S5)

newdat <- expand.grid(Depth.m=c("10","21"),Site=c("5")) 

H11_S5 <- m_11_S5
pred_11_S5 <- as.data.frame(cbind(newdat,predict.glm(H11_S5,newdata=newdat,type="link",se.fit=T)))
pred_11_S5 <- cbind(Haplotype = "Haplotype 11", pred_11_S5)
pred_11_S5$prop <- plogis(pred_11_S5$fit)
pred_11_S5$prop.lwr <- plogis(pred_11_S5$fit - 2*pred_11_S5$se.fit)
pred_11_S5$prop.upr <- plogis(pred_11_S5$fit + 2*pred_11_S5$se.fit)

pred_11 <- rbind(pred_11_S1,pred_11_S2,pred_11_S4,pred_11_S5)

#### Haplotype 1a_Pe
d$Haplotype_1a_Pe <- ifelse(d$Species.haplotype=="Haplotype 1a_Pe",1,0)
with(d[d$Site=="1",],table(Haplotype_1a_Pe,Depth.m))

## Site 1
d_1a_Pe_S1<-d %>% 
  filter((Site=="1"))

m_1a_Pe_S1 <- glm(Haplotype_1a_Pe ~ Depth.m, family="binomial", data=d_1a_Pe_S1)

Anova(m_1a_Pe_S1)
summary(m_1a_Pe_S1)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("1")) 

Pe_S1 <- m_1a_Pe_S1 
pred_1a_Pe_S1 <- as.data.frame(cbind(newdat,predict.glm(Pe_S1,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pe_S1 <- cbind(Haplotype = "Haplotype 1a_Pe", pred_1a_Pe_S1)
pred_1a_Pe_S1$prop <- plogis(pred_1a_Pe_S1$fit)
pred_1a_Pe_S1$prop.lwr <- plogis(pred_1a_Pe_S1$fit - 2*pred_1a_Pe_S1$se.fit)
pred_1a_Pe_S1$prop.upr <- plogis(pred_1a_Pe_S1$fit + 2*pred_1a_Pe_S1$se.fit)

## Site 2
d_1a_Pe_S2<-d %>% 
  filter((Site=="2"))

m_1a_Pe_S2 <- glm(Haplotype_1a_Pe ~ Depth.m, family="binomial", data=d_1a_Pe_S2)

Anova(m_1a_Pe_S2)
summary(m_1a_Pe_S2)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("2")) 

Pe_S2 <- m_1a_Pe_S2
pred_1a_Pe_S2 <- as.data.frame(cbind(newdat,predict.glm(Pe_S2,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pe_S2 <- cbind(Haplotype = "Haplotype 1a_Pe", pred_1a_Pe_S2)
pred_1a_Pe_S2$prop <- plogis(pred_1a_Pe_S2$fit)
pred_1a_Pe_S2$prop.lwr <- plogis(pred_1a_Pe_S2$fit - 2*pred_1a_Pe_S2$se.fit)
pred_1a_Pe_S2$prop.upr <- plogis(pred_1a_Pe_S2$fit + 2*pred_1a_Pe_S2$se.fit)

## Site 4
d_1a_Pe_S4<-d %>% 
  filter((Site=="4"))

m_1a_Pe_S4 <- glm(Haplotype_1a_Pe ~ Depth.m, family="binomial", data=d_1a_Pe_S4)

Anova(m_1a_Pe_S4)
summary(m_1a_Pe_S4)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("4")) 

Pe_S4 <- m_1a_Pe_S4
pred_1a_Pe_S4 <- as.data.frame(cbind(newdat,predict.glm(Pe_S4,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pe_S4 <- cbind(Haplotype = "Haplotype 1a_Pe", pred_1a_Pe_S4)
pred_1a_Pe_S4$prop <- plogis(pred_1a_Pe_S4$fit)
pred_1a_Pe_S4$prop.lwr <- plogis(pred_1a_Pe_S4$fit - 2*pred_1a_Pe_S4$se.fit)
pred_1a_Pe_S4$prop.upr <- plogis(pred_1a_Pe_S4$fit + 2*pred_1a_Pe_S4$se.fit)

## Site 5
d_1a_Pe_S5<-d %>% 
  filter((Site=="5"))

m_1a_Pe_S5 <- glm(Haplotype_1a_Pe ~ Depth.m, family="binomial", data=d_1a_Pe_S5)

Anova(m_1a_Pe_S5)
summary(m_1a_Pe_S5)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("5")) 

Pe_S5 <- m_1a_Pe_S5
pred_1a_Pe_S5 <- as.data.frame(cbind(newdat,predict.glm(Pe_S5,newdata=newdat,type="link",se.fit=T)))
pred_1a_Pe_S5 <- cbind(Haplotype = "Haplotype 1a_Pe", pred_1a_Pe_S5)
pred_1a_Pe_S5$prop <- plogis(pred_1a_Pe_S5$fit)
pred_1a_Pe_S5$prop.lwr <- plogis(pred_1a_Pe_S5$fit - 2*pred_1a_Pe_S5$se.fit)
pred_1a_Pe_S5$prop.upr <- plogis(pred_1a_Pe_S5$fit + 2*pred_1a_Pe_S5$se.fit)

pred_Pe <- rbind(pred_1a_Pe_S1,pred_1a_Pe_S2,pred_1a_Pe_S4,pred_1a_Pe_S5)

#### Haplotype 3b
d$Haplotype_3b <- ifelse(d$Species.haplotype=="Haplotype 3b",1,0)
with(d[d$Site=="1",],table(Haplotype_3b,Depth.m))

## Site 1
d_3b_S1<-d %>% 
  filter((Site=="1"),
         !(Depth.m=="5"))

m_3b_S1 <- glm(Haplotype_3b ~ Depth.m, family="binomial", data=d_3b_S1)

Anova(m_3b_S1)
summary(m_3b_S1)

newdat <- expand.grid(Depth.m=c("10","21"),Site=c("1")) 

H3b_S1 <- m_3b_S1 
pred_3b_S1 <- as.data.frame(cbind(newdat,predict.glm(H3b_S1,newdata=newdat,type="link",se.fit=T)))
pred_3b_S1 <- cbind(Haplotype = "Haplotype 3b", pred_3b_S1)
pred_3b_S1$prop <- plogis(pred_3b_S1$fit)
pred_3b_S1$prop.lwr <- plogis(pred_3b_S1$fit - 2*pred_3b_S1$se.fit)
pred_3b_S1$prop.upr <- plogis(pred_3b_S1$fit + 2*pred_3b_S1$se.fit)

## Site 2
d_3b_S2<-d %>% 
  filter((Site=="2"))

m_3b_S2 <- glm(Haplotype_3b ~ Depth.m, family="binomial", data=d_3b_S2)

Anova(m_3b_S2)
summary(m_3b_S2)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("2")) 

H3b_S2 <- m_3b_S2
pred_3b_S2 <- as.data.frame(cbind(newdat,predict.glm(H3b_S2,newdata=newdat,type="link",se.fit=T)))
pred_3b_S2 <- cbind(Haplotype = "Haplotype 3b", pred_3b_S2)
pred_3b_S2$prop <- plogis(pred_3b_S2$fit)
pred_3b_S2$prop.lwr <- plogis(pred_3b_S2$fit - 2*pred_3b_S2$se.fit)
pred_3b_S2$prop.upr <- plogis(pred_3b_S2$fit + 2*pred_3b_S2$se.fit)

## Site 4
d_3b_S4<-d %>% 
  filter((Site=="4"),
         !(Depth.m=="5"))

m_3b_S4 <- glm(Haplotype_3b ~ Depth.m, family="binomial", data=d_3b_S4)

Anova(m_3b_S4)
summary(m_3b_S4)

newdat <- expand.grid(Depth.m=c("10","21"),Site=c("4")) 

H3b_S4 <- m_3b_S4
pred_3b_S4 <- as.data.frame(cbind(newdat,predict.glm(H3b_S4,newdata=newdat,type="link",se.fit=T)))
pred_3b_S4 <- cbind(Haplotype = "Haplotype 3b", pred_3b_S4)
pred_3b_S4$prop <- plogis(pred_3b_S4$fit)
pred_3b_S4$prop.lwr <- plogis(pred_3b_S4$fit - 2*pred_3b_S4$se.fit)
pred_3b_S4$prop.upr <- plogis(pred_3b_S4$fit + 2*pred_3b_S4$se.fit)

pred_3b <- rbind(pred_3b_S1,pred_3b_S2,pred_3b_S4)

d$Haplotype_1a_Pe <- ifelse(d$Species.haplotype=="Haplotype 1a_Pe",1,0)
with(d[d$Site=="1",],table(Haplotype_1a_Pe,Depth.m))


#### Haplotype 10
d$Haplotype_10 <- ifelse(d$Species.haplotype=="Haplotype 10",1,0)
with(d[d$Site=="1",],table(Haplotype_10,Depth.m))

## Site 1
d_H10_S1<-d %>% 
  filter((Site=="1"))

m_H10_S1 <- glm(Haplotype_10 ~ Depth.m, family="binomial", data=d_H10_S1)

Anova(m_H10_S1)
summary(m_H10_S1)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("1")) 

H10_S1 <- m_H10_S1 
pred_H10_S1 <- as.data.frame(cbind(newdat,predict.glm(H10_S1,newdata=newdat,type="link",se.fit=T)))
pred_H10_S1 <- cbind(Haplotype = "Haplotype 10", pred_H10_S1)
pred_H10_S1$prop <- plogis(pred_H10_S1$fit)
pred_H10_S1$prop.lwr <- plogis(pred_H10_S1$fit - 2*pred_H10_S1$se.fit)
pred_H10_S1$prop.upr <- plogis(pred_H10_S1$fit + 2*pred_H10_S1$se.fit)

## Site 2
d_H10_S2<-d %>% 
  filter((Site=="2"))

m_H10_S2 <- glm(Haplotype_10 ~ Depth.m, family="binomial", data=d_H10_S2)

Anova(m_H10_S2)
summary(m_H10_S2)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("2")) 

H10_S2 <- m_H10_S2
pred_H10_S2 <- as.data.frame(cbind(newdat,predict.glm(H10_S2,newdata=newdat,type="link",se.fit=T)))
pred_H10_S2 <- cbind(Haplotype = "Haplotype 10", pred_H10_S2)
pred_H10_S2$prop <- plogis(pred_H10_S2$fit)
pred_H10_S2$prop.lwr <- plogis(pred_H10_S2$fit - 2*pred_H10_S2$se.fit)
pred_H10_S2$prop.upr <- plogis(pred_H10_S2$fit + 2*pred_H10_S2$se.fit)

## Site 4
d_H10_S4<-d %>% 
  filter((Site=="4"))

m_H10_S4 <- glm(Haplotype_10 ~ Depth.m, family="binomial", data=d_H10_S4)

Anova(m_H10_S4)
summary(m_H10_S4)

newdat <- expand.grid(Depth.m=c("5","10","21"),Site=c("4")) 

H10_S4 <- m_H10_S4
pred_H10_S4 <- as.data.frame(cbind(newdat,predict.glm(H10_S4,newdata=newdat,type="link",se.fit=T)))
pred_H10_S4 <- cbind(Haplotype = "Haplotype 10", pred_H10_S4)
pred_H10_S4$prop <- plogis(pred_H10_S4$fit)
pred_H10_S4$prop.lwr <- plogis(pred_H10_S4$fit - 2*pred_H10_S4$se.fit)
pred_H10_S4$prop.upr <- plogis(pred_H10_S4$fit + 2*pred_H10_S4$se.fit)

## Site 5
d_H10_S5<-d %>% 
  filter((Site=="5"),
         !(Depth.m=="5"))

m_H10_S5 <- glm(Haplotype_10 ~ Depth.m, family="binomial", data=d_H10_S5)

Anova(m_H10_S5)
summary(m_H10_S5)

newdat <- expand.grid(Depth.m=c("10","21"),Site=c("5")) 

H10_S5 <- m_H10_S5
pred_H10_S5 <- as.data.frame(cbind(newdat,predict.glm(H10_S5,newdata=newdat,type="link",se.fit=T)))
pred_H10_S5 <- cbind(Haplotype = "Haplotype 10", pred_H10_S5)
pred_H10_S5$prop <- plogis(pred_H10_S5$fit)
pred_H10_S5$prop.lwr <- plogis(pred_H10_S5$fit - 2*pred_H10_S5$se.fit)
pred_H10_S5$prop.upr <- plogis(pred_H10_S5$fit + 2*pred_H10_S5$se.fit)

pred_10 <- rbind(pred_H10_S1,pred_H10_S2,pred_H10_S4,pred_H10_S5)

## Merge predicted data into one dataframe 
pred <- rbind(pred_Pm,pred_8a,pred_11,pred_Pe,pred_3b,pred_10)


#### Plot
#Remove observations where there were no samples
pred<-pred %>% 
  filter(!(Haplotype=="Haplotype 1a_Pm" & Site=="4" & Depth.m=="21"),
         !(Haplotype=="Haplotype 8a" & Site=="1" & Depth.m=="21"),
         !(Haplotype=="Haplotype 8a" & Site=="4" & Depth.m=="21"),
         !(Haplotype=="Haplotype 11" & Site=="1" & Depth.m=="5"),
         !(Haplotype=="Haplotype 11" & Site=="1" & Depth.m=="21"),
         !(Haplotype=="Haplotype 11" & Site=="2" & Depth.m=="5"),
         !(Haplotype=="Haplotype 11" & Site=="2" & Depth.m=="21"),
         !(Haplotype=="Haplotype 11" & Site=="5" & Depth.m=="5"),
         !(Haplotype=="Haplotype 3b" & Site=="1" & Depth.m=="5"),
         !(Haplotype=="Haplotype 3b" & Site=="4" & Depth.m=="5"),
         !(Haplotype=="Haplotype 3b" & Site=="5" & Depth.m=="5"),
         !(Haplotype=="Haplotype 3b" & Site=="5" & Depth.m=="10"),
         !(Haplotype=="Haplotype 3b" & Site=="5" & Depth.m=="21"),
         !(Haplotype=="Haplotype 10" & Site=="5" & Depth.m=="5"))


#### Plot
#Remove error bars for Hap 11 from sites 1 & 2
pred_errbars<-pred %>% 
  filter(!(Haplotype=="Haplotype 11" & Site=="1" & Depth.m=="10"),
         !(Haplotype=="Haplotype 11" & Site=="2" & Depth.m=="10"))


# Reorder haplotypes and depth
pred$Haplotype = factor(pred$Haplotype, levels=c("Haplotype 1a_Pm","Haplotype 8a","Haplotype 11","Haplotype 1a_Pe","Haplotype 3b","Haplotype 10"))

pred_errbars$Haplotype = factor(pred_errbars$Haplotype, levels=c("Haplotype 1a_Pm","Haplotype 8a","Haplotype 11","Haplotype 1a_Pe","Haplotype 3b","Haplotype 10"))

pred$Depth.m = factor(pred$Depth.m,levels=c("21","10","5"))

## Bar graph
ggplot(pred, aes(x=Depth.m, y=prop, color=Haplotype, fill=Haplotype))+
  geom_bar(stat="identity", colour=NA)+
  geom_errorbar(data=pred_errbars, aes(ymin=prop.lwr, ymax=prop.upr), width=0.4, size=0.5, alpha=.7, color="black")+
  labs(x="Depth (m)", y="Proportion of sampled colonies")+
  theme_bw(base_size = 14)+
  theme(strip.text.y = element_blank())+
  scale_fill_manual(name="Species/haplotype",
                    values = c("Haplotype 1a_Pm"= "#0072B2",
                                "Haplotype 1a_Pe"= "#56B4E9",
                                "Haplotype 3b"= "#E69F00",
                                "Haplotype 10"= "#D55E00",
                                "Haplotype 8a"= "#CC79A7",
                                "Haplotype 11"= "#009E73"),
                    labels=c("P. meandrina","Haplotype 8a","Haplotype 11","P. eydouxi","P. verrucosa 3b","Haplotype 10"))+
  scale_y_continuous(limits = c(0,1), breaks = c(seq(0,1,.2)))+
  scale_x_discrete(breaks=c("5","10","21"),
                   labels=c("5", "10", "20"))+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(), strip.text = element_text(size = 10))+
  facet_grid(Haplotype~Site, labeller = labeller(Site=c("1"="Site 1","2"="Site 2","4"="Site 4","5"="Site 5")))

