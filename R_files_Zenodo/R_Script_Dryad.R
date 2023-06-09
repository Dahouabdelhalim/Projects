library(reshape2)
library(ggplot2)
library(vegan)
library(indicspecies)
library(tidyverse)
library(gridExtra)
library(grid)
library(png)
library(cowplot)
library(dplyr)
library(indicspecies)

## 1. Uneven sampling t-tests Tests for statistical difference between samples collected from Innlandet, France
#uneven sampling t-test 210521
spr_new <- read.csv("~/unevensampling_test.csv", sep=";")

#testing within Troms
Troms_spr <- subset (spr_new, spr_new$Location == "Troms")
Troms_spr$Sum<-as.integer(Troms_spr$Sum)
anova_one_way<-aov(Sum~Season, data = Troms_spr)
summary(anova_one_way)
TukeyHSD(anova_one_way) #if anova not significant, don't have to run Tukey test

#testing within Innlandet
Innlandet_spr <- subset (spr_new, spr_new$Location == "Innlandet")
Innlandet_spr$Sum<-as.integer(Innlandet_spr$Sum)
anova_one_way<-aov(Sum~Season, data = Innlandet_spr)
summary(anova_one_way)

#testing for both Troms and Innlandet combined
both_spr <-rbind(Innlandet_spr,Troms_spr)
both_spr$Sum<-as.integer(both_spr$Sum)
anova_one_way<-aov(Sum~Season, data = both_spr)
summary(anova_one_way)


## 2. Check sampling effort Tests for statistical difference between samples collected from Innlandet, France
accu <- read.csv("~/accu_new_090221.csv", sep=";")
sp1<- specaccum(accu [,-c(1:8)])
sp2<- specaccum(accu [,-c(1:8)], method = "random")

## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")

## 3. perMANOVA (dietary richness) Tests for differences between variables (location, sex, season)
spr_new <- read.csv("~spr_permanova_210521.csv", sep=";")
com1=spr_new[,9:ncol(spr_new)]
spp_distmat<- vegdist(com1,method="bray")
spp_distmat <- as.matrix(spp_distmat, labels = T)
set.seed(123)
spp_NMS <-metaMDS(spp_distmat,distance = "bray", k = 3,maxit = 999,   trymax = 500, wascores = TRUE)
spp_NMS
stressplot(spp_NMS)
plot(spp_NMS)
adonis2(spp_distmat ~Location + Sex + Season + Location*Sex + Location*Season + Sex*Season, data=spr_new, p.adjust.m="bonferroni", sim.method = 'bray')

## 4. perMANOVA (dietary composition) Tests for differences between variables (location, sex, season)
accu_new <- read.csv("~/accu_new_090221.csv", sep=";")
#convert to relative abundance https://rpubs.com/CPEL/NMDS
com1=accu_new[,9:ncol(accu_new)]
spp.rel<-decostand(com1,method="total")

#calculate distance matrix
spp_distmat<- vegdist(spp.rel,method="bray")
spp_distmat <- as.matrix(spp_distmat, labels = T)

#NMDS in vegan
set.seed(123)
spp_NMS <-metaMDS(spp_distmat,distance = "bray", k = 3,maxit = 999,   trymax = 500, wascores = TRUE)
spp_NMS
stressplot(spp_NMS)
plot(spp_NMS)
data.scores = as.data.frame(scores(spp_NMS))
data.scores$Country = accu_new$Country
data.scores$Location = accu_new$Location
data.scores$Sex = accu_new$Sex
data.scores$Season = accu_new$Season
adonis2(spp_distmat ~Location + Sex + Season + Location*Sex + Location*Season + Sex*Season, data=accu_new, p.adjust.m="bonferroni", sim.method = 'bray')

### 4.1 perMANOVA (dietary composition) Tests for differences between variables (location, sex, season), based on only Norwegian or French samples
#testing diet composition within Norway
accu_norway <- read.csv("~/accu_norway.csv", sep=";")
com1=accu_norway[,9:ncol(accu_norway)]
spp.rel<-decostand(com1,method="total")
spp_distmat<- vegdist(spp.rel,method="bray")
spp_distmat <- as.matrix(spp_distmat, labels = T)
adonis2(spp_distmat ~Location + Sex + Season + Location*Sex + Location*Season + Sex*Season, data=accu_norway, p.adjust.m="bonferroni", sim.method = 'bray')
groups_mine<-factor(c(rep(1,6),rep(2,53),rep(3,19), rep(4,50)),labels=c("Finnmark","Innlandet","Troms", "Trondelag"))
groups_test<- data.frame(groups_mine)
dismine<-vegdist(spp.rel)
mod_mine<-betadisper(dismine,groups_mine,bias.adjust=TRUE)
mod_mine
TukeyHSD(mod_mine)

#testing diet composition within france
accu_france <- read.csv("~/accu_france.csv", sep=";")
com2=accu_france[,9:ncol(accu_france)]
spp.rel1<-decostand(com2,method="total")
spp_distmat1<- vegdist(spp.rel1,method="bray")
spp_distmat1 <- as.matrix(spp_distmat1, labels = T)
adonis2(spp_distmat1 ~Location + Sex + Location*Sex, data=accu_france, p.adjust.m="bonferroni", sim.method = 'bray')

### 4.2 indicator species Test which plant taxa were strongly associated with given variable
abund = accu_new[,9:ncol(accu_new)]
location = accu_new$Location
sex = accu_new$Sex
season = accu_new$Season
inv = multipatt(abund, location, func = "r.g", control = how(nperm=999))
summary(inv)
inv = multipatt(abund, sex, func = "r.g", control = how(nperm=999))
summary(inv)
inv = multipatt(abund, season, func = "r.g", control = how(nperm=999))
summary(inv)



