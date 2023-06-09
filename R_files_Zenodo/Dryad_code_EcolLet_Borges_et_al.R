#### Set working directory ####

setwd("/Users/eu/Michigan State University/Fitzpatrick, Sarah - Isabela/guppy_movement_ms/data and models")

#### Loading packages ####

library(pscl)
library(ggplot2)
library(MASS)
library(boot)
library(MuMIn)
library(sjPlot)
library(lme4)
library(texreg)
library(glmmTMB)
library(bbmle)
library(piecewiseSEM)
library(dplyr)
library(modelsummary)
library(DHARMa)
library(gridExtra)
library(cowplot)

##### 1. Load data #####

# main movement and fitness data
d = read.csv("taylor_caigual_fitness.csv", sep=";")

# seasonality data
c.disp.less <- read.csv("c.disp.subsetted.csv", sep=";")
t.disp <- read.csv("taylor.dispersed.updown.csv", sep=";")
all.season <- read.csv("both_streams_seasonal_updown.csv", sep=";")
all.d <- subset(all.season, Season=="D")
all.w <- subset(all.season, Season=="W")

##### 2. Wrangle data #####

# subset capture occasions before onset of gene flow
d = subset(d, d$Capt.event>5 & d$Capt.event<17)

##### 3. Centering and transforming #####

hist(d$total_dist) # poisson
d$total_dist_l = log(d$total_dist)
hist(d$total_dist_l) # normal

hist(d$LRS) # poisson
d$LRS_l = log(d$LRS)
hist(d$LRS_l) # left skewed

# adding dispersal status: dispersers = TRUE 
d$disp.status = as.factor(d$total_dist>=10)
plot(d$disp.status)
prop.disp.status = length(d$disp.status[d$disp.status == TRUE])/length(d$disp.status)
prop.non.disp.status = 1 - prop.disp.status
boxplot(data=d, total_dist ~ disp.status)

# make Sex = 1 or 2 (1 for F, 2 for M)
d$Sex <- sapply(d$Sex.x, is.logical)
d$Sex <- sapply(d$Sex.x, as.numeric)
head(d)
plot(Sex ~ Sex.x, data=d) 

d$total_dist_10 <- d$total_dist
d$total_dist_10[d$total_dist_10 < 10] <- NA

# make it so that longevity starts from 0, if needed 

hist(d$Longevity, xlim=c(0,20))
d$Longevity_adjusted <- d$Longevity - 2
hist(d$Longevity_adjusted, xlim=c(0,20))

# creating subsets of data for each sex for each stream
d.tay = subset(d, Stream=="Taylor")
d.cai = subset(d, Stream=="Caigual")

d.m = subset(d, Sex.x=="M")
d.f = subset(d, Sex.x=="F")

d.tay.f = subset(d.tay, Sex.x=="F")
d.tay.m = subset(d.tay, Sex.x=="M")

d.cai.f = subset(d.cai, Sex.x=="F")
d.cai.m = subset(d.cai, Sex.x=="M")

# dispersers-only subsets
d.disp = subset(d, total_dist>=10)

d.disp.f = subset(d.disp, Sex.x=="F")
d.disp.m = subset(d.disp, Sex.x=="M")

d.non.disp = subset(d, total_dist<10)

d.non.disp.f = subset(d.non.disp, Sex.x=="F")
d.non.disp.m = subset(d.non.disp, Sex.x=="M")

d.tay.dist <- subset(d.tay, d.tay$total_dist >= 10)
d.cai.dist <- subset(d.cai, d.cai$total_dist >= 10)
d.f.dist <- subset(d.f, d.f$total_dist >= 10)
d.m.dist <- subset(d.m, d.m$total_dist >= 10)
d.tay.rep <- subset(d.tay, d.tay$LRS >= 1)
d.cai.rep <- subset(d.cai, d.cai$LRS >= 1)

#### 4. Analyses: Total distance and status ####

# status models
hstatus.all <- glm(disp.status ~ Sex.x + Stream*hindex + Longevity, data=d, family="binomial")
hstatus.m <- glm(disp.status ~Stream*hindex + Longevity + Male_SL, data=d.m, family="binomial")
status.all <- glm(disp.status ~ Sex.x + Stream + Longevity, data=d, family="binomial")
status.m <- glm(disp.status ~ Stream + Longevity + Male_SL, data=d.m, family="binomial")

# distance models
glmmTMB.dist.all <- glmmTMB(total_dist ~ Sex.x + Stream*hindex + log(Longevity),
                            ziformula = ~ Sex.x + Stream*hindex + log(Longevity),
                            family=nbinom2,data=d)
glmmTMB.dist.m <- glmmTMB(total_dist ~ Male_SL + Stream*hindex + log(Longevity), 
                              ziformula = ~ Male_SL + Stream*hindex + log(Longevity),
                              family=nbinom2,data=d.m)
hzi.dist.all <- zeroinfl(total_dist ~ Sex.x + Stream*hindex + Longevity, data=d, dist="negbin")
hzi.dist.m <- zeroinfl(total_dist ~ Male_SL + Stream*hindex + Longevity, data=d.m, dist="negbin")
zi.dist.all <- zeroinfl(total_dist ~ Sex.x + Stream + Longevity, data=d, dist="negbin")
zi.dist.m <- zeroinfl(total_dist ~ Male_SL + Stream + Longevity, data=d.m, dist="negbin")

# model selection and testing for model adequacy

AICctab(status.all,hstatus.all) 
AICctab(zi.dist.all,hzi.dist.all,glmmTMB.dist.all) 
AICctab(zi.dist.m,hzi.dist.m,glmmTMB.dist.m) 

testDispersion(hstatus.all)
hstatus.all.sim <- simulateResiduals(fittedModel = hstatus.all, plot = F)
plot(hstatus.all.sim)

testDispersion(hstatus.m)
hstatus.m.sim <- simulateResiduals(fittedModel = hstatus.m, plot = F)
plot(hstatus.m.sim) 

testDispersion(glmmTMB.dist.all)
glmmTMB.dist.all.sim <- simulateResiduals(fittedModel = glmmTMB.dist.all, plot = F)
plot(glmmTMB.dist.all.sim) 

testDispersion(glmmTMB.dist.m)
glmmTMB.dist.males.sim <- simulateResiduals(fittedModel = glmmTMB.dist.m, plot = F)
plot(glmmTMB.dist.males.sim) 

#### 5. Analyses: Lifetime reproductive success (LRS) ####

# male status models
glmmTMB.hZINB.status.m <- glmmTMB(LRS ~ disp.status + Stream*hindex + Longevity + Male_SL, 
                                  ziformula = ~ disp.status + Stream*hindex + Longevity + Male_SL,
                                  family=nbinom2,data=d.m)
ZINB.status.m  <- zeroinfl(LRS ~ disp.status + Stream + Longevity + Male_SL, 		data=d.m, dist="negbin")
hZINB.status.m  <- zeroinfl(LRS ~ disp.status + Stream*hindex + Longevity + Male_SL, 		data=d.m, dist="negbin")

# male distance models
glmmTMB.hZINB.dist.m <- glmmTMB(LRS ~ total_dist + Stream*hindex + Longevity + Male_SL, 
                                ziformula = ~ total_dist + Stream*hindex + Longevity + Male_SL,
                                family=nbinom2,data=d.m)
ZINB.dist.m <- zeroinfl(LRS ~ total_dist + Stream + Longevity + Male_SL, data=d.m, 		dist="negbin")
hZINB.dist.m <- zeroinfl(LRS ~ total_dist + Stream*hindex + Longevity + Male_SL, data=d.m, 		dist="negbin")

# female status models
glmmTMB.hZINB.status.f <- glmmTMB(LRS ~ disp.status + Stream*hindex + Longevity, 
                                  ziformula = ~ disp.status + Stream*hindex + Longevity,
                                  family=nbinom2,data=d.f)
ZINB.status.f  <- zeroinfl(LRS ~ disp.status + Stream + Longevity, data=d.f, 		dist="negbin")
hZINB.status.f  <- zeroinfl(LRS ~ disp.status + Stream*hindex + Longevity, data=d.f, 		dist="negbin")

# female distance models
glmmTMB.hZINB.dist.f <- glmmTMB(LRS ~ total_dist + Stream*hindex + Longevity, 
                                ziformula = ~ total_dist + Stream*hindex + Longevity,
                                family=nbinom2,data=d.f)
ZINB.dist.f <- zeroinfl(LRS ~ total_dist + Stream + Longevity, data=d.f, dist="negbin")
hZINB.dist.f <- zeroinfl(LRS ~ total_dist + Stream*hindex + Longevity, data=d.f, dist="negbin")

# model selection and testing for model adequacy

AICctab(ZINB.status.m,hZINB.status.m,glmmTMB.hZINB.status.m) 
AICctab(ZINB.status.f,hZINB.status.f,glmmTMB.hZINB.status.f) 
AICctab(ZINB.dist.m,hZINB.dist.m,glmmTMB.hZINB.dist.m) 
AICctab(ZINB.dist.f,hZINB.dist.f,glmmTMB.hZINB.dist.f) 

testDispersion(glmmTMB.hZINB.status.m)
glmmTMB.hZINB.status.m.sim <- simulateResiduals(fittedModel = glmmTMB.hZINB.status.m, plot = F)
plot(glmmTMB.hZINB.status.m.sim) 

testDispersion(glmmTMB.hZINB.dist.m)
glmmTMB.hZINB.dist.m.sim <- simulateResiduals(fittedModel = glmmTMB.hZINB.dist.m, plot = F)
plot(glmmTMB.hZINB.dist.m.sim) 

testDispersion(glmmTMB.hZINB.status.f)
glmmTMB.hZINB.status.f.sim <- simulateResiduals(fittedModel = glmmTMB.hZINB.status.f, plot = F)
plot(glmmTMB.hZINB.status.f.sim) 

testDispersion(glmmTMB.hZINB.dist.f)
glmmTMB.hZINB.dist.f.sim <- simulateResiduals(fittedModel = glmmTMB.hZINB.dist.f, plot = F)
plot(glmmTMB.hZINB.dist.f.sim) 

#### 6. Analyses: Number of mates ####

# male models
glmmTMB.hnm_ZINB_1m <- glmmTMB(NumMates ~ net_range + Stream*hindex + Longevity,
                               ziformula = ~ net_range + Stream*hindex + Longevity,
                               family=nbinom2,data=d.m)
nm_ZINB_1m <- zeroinfl(NumMates ~ net_range + Stream + Longevity, data=d.m, dist="negbin")
hnm_ZINB_1m <- zeroinfl(NumMates ~ net_range + Stream*hindex + Longevity, data=d.m, 		dist="negbin")

# female models
glmmTMB.hnm_ZINB_1f <- glmmTMB(NumMates ~ net_range + Stream*hindex + Longevity,
                               ziformula = ~ net_range + Stream*hindex + Longevity,
                               family=nbinom2,data=d.f)
nm_ZINB_1f <- zeroinfl(NumMates ~ net_range + Stream + Longevity, data=d.f, dist="negbin")
hnm_ZINB_1f <- zeroinfl(NumMates ~ net_range + Stream*hindex + Longevity, data=d.f, dist="negbin")

# model selection and testing for model adequacy

AICctab(nm_ZINB_1m,hnm_ZINB_1m,glmmTMB.hnm_ZINB_1m)
AICctab(nm_ZINB_1f,hnm_ZINB_1f,glmmTMB.hnm_ZINB_1f)

testDispersion(glmmTMB.hnm_ZINB_1m)
glmmTMB.hnm_ZINB_1m.sim <- simulateResiduals(fittedModel = glmmTMB.hnm_ZINB_1m, plot = F)
plot(glmmTMB.hnm_ZINB_1m.sim)

testDispersion(glmmTMB.hnm_ZINB_1f)
glmmTMB.hnm_ZINB_1f.sim <- simulateResiduals(fittedModel = glmmTMB.hnm_ZINB_1f, plot = F)
plot(glmmTMB.hnm_ZINB_1f.sim)

#### 7. Analyses: Seasonality and stream effects ####

# Taylor models
nb.1.tay <- glmmTMB(data=t.disp, min_dist ~ Season + (1|FishID_nodash), family='nbinom2')
nb.2.tay <- glmmTMB(data=t.disp, min_dist ~ Season + (1|FishID_nodash) + (1|month), family='nbinom2')

# Caigual models
nb.1.cai <- glmmTMB(data=c.disp.less, min_dist ~ Season + (1|FishID_nodash), family='nbinom2')
nb.2.cai <- glmmTMB(data=c.disp.less, min_dist ~ Season + (1|FishID_nodash) + (1|Month), family='nbinom2')
p.nb.2.cai <- glmmTMB(data=c.disp.less, min_dist ~ Season + (1|FishID_nodash) + (1|Month), 
                      family='poisson')

# model selection

AICtab(nb.1.tay, nb.2.tay, weights=T)
AICtab(nb.1.cai, nb.2.cai, weights=T)

#population level: 

up.c.nb <- glm.nb(data=c.disp.less, upstream ~ Season) 
down.c.nb <- glm.nb(data=c.disp.less, downstream*-1 ~ Season) 
summary(up.c.nb)
summary(down.c.nb)

up.t.nb <- glm.nb(data=t.disp, upstream ~ Season) 
plot(data=t.disp, upstream ~ Season)
down.t.nb <- glm.nb(data=t.disp, downstream*-1 ~ Season) 
summary(up.t.nb)
summary(down.t.nb)

#### 8. Calculating values for Table 1 ####

d.disp.f.cai = subset(d.disp.f, Stream=="Caigual")
d.disp.f.tay = subset(d.disp.f, Stream=="Taylor")

d.nondisp.f.cai = subset(d.non.disp.f, Stream=="Caigual")
d.nondisp.f.tay = subset(d.non.disp.f, Stream=="Taylor")

d.disp.m.cai = subset(d.disp.m, Stream=="Caigual")
d.disp.m.tay = subset(d.disp.m, Stream=="Taylor")

d.nondisp.m.cai = subset(d.non.disp.m, Stream=="Caigual")
d.nondisp.m.tay = subset(d.non.disp.m, Stream=="Taylor")

#Taylor females:
summary(d.tay.f$total_dist)
summary(d.tay.f$net_range)
summary(d.nondisp.f.tay$LRS)
summary(d.disp.f.tay$LRS)
summary(d.nondisp.f.tay$NumMates)
summary(d.disp.f.tay$NumMates)

#Taylor males:
summary(d.tay.m$total_dist)
summary(d.tay.m$net_range)
summary(d.nondisp.m.tay$LRS)
summary(d.disp.m.tay$LRS)
summary(d.nondisp.m.tay$NumMates)
summary(d.disp.m.tay$NumMates)

#Caigual females:
summary(d.cai.f$total_dist)
summary(d.cai.f$net_range)
summary(d.nondisp.f.cai$LRS)
summary(d.disp.f.cai$LRS)
summary(d.nondisp.f.cai$NumMates)
summary(d.disp.f.cai$NumMates)

#Caigual males:
summary(d.cai.m$total_dist)
summary(d.cai.m$net_range)
summary(d.nondisp.m.cai$LRS)
summary(d.disp.m.cai$LRS)
summary(d.nondisp.m.cai$NumMates)
summary(d.disp.m.cai$NumMates)