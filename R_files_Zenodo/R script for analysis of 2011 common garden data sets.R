# R script for data analysis of sunflower phenotypic data
# for Corbi et al 2017.  Eric Baack, Luther College.  Final version July 2016.

# Analysis compares three populations using ANOVA:  G1, G3IA, and G3ND to ask
# if means changed after two generations of selection in either environment.

# Bootstrapped LSmeans are generated for the three hybrid populations, plus
# the two parent populations (HA89, Ann1238), except for some variables where 
# parent populations lacked sufficient replication, in which case just the hybrid
# populations were used.

# Analysis was carried out independently on two common garden data sets, one
# from Fargo, ND, the other from Decorah, IA.  
# Common garden data were collected in 2011.
# Change the line below to the directory where the data files are stored
setwd("C:/")
ND2011R <- read.csv("ND2011R.csv")
IA2011R <- read.csv("IA2011R.csv")

NDHYB <- subset(ND2011R, Asrc == 'G1'|Asrc =='ND'|Asrc =='IA')
IAHYB <- subset(IA2011R, Asrc == 'G1'|Asrc =='ND'|Asrc =='IA')

# script assumes that packages 'lsmeans' and 'boot' have been installed
# If these are not previously installed, run the following two lines

#install.packages("lsmeans")
#install.packages("boot")


# The following data sets are used below
# They can be run with either the field data from ND, as currently shown
# or using the field data from Iowa, as in the commented out lines.

#allplants <- ND2011R
#hybplants <- NDHYB
allplants <- IA2011R
hybplants <- IAHYB

no_headless_hyb <- hybplants['Tothds' > 0,]
no_headless_all <- allplants['Tothds' > 0,]
hybsds <- hybplants['masspsd'> 0,]
allsds <- allplants['masspsd'>0,]
hybsdar <- hybsds['avgsdar'>0]
allsdar <- allsds['avgsdar'>0]
hyblf <- subset(hybplants, lfr > 0 & !(is.na(lfr)))
hybstm <- subset(hybplants, Stmlen > 0 & !(is.na(Stmlen)))
 
# libraries to calculate LSmeans and perform bootstraps
library(lsmeans)
library(boot)

# Days to flowering
Fdays.aov <- aov(Fdays ~ Asrc + Plt, data = hybplants)
anova(Fdays.aov)
LSFdays <- lsmeans(Fdays.aov, pairwise ~ Asrc)
LSFdays


FD.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Fdays ~ Asrc + Plt, data= AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

Fdays.boot <- boot(allplants, FD.BOOT, R = 5000)
boot.ci(Fdays.boot, index = 1, conf = .95, type = "bca")
boot.ci(Fdays.boot, index = 2, conf = .95, type = "bca")
boot.ci(Fdays.boot, index = 3, conf = .95, type = "bca")
boot.ci(Fdays.boot, index = 4, conf = .95, type = "bca")
boot.ci(Fdays.boot, index = 5, conf = .95, type = "bca")

#Branch number
Br.aov <- aov(Br_num~Asrc + Plt, data = hybplants)
anova(Br.aov)
LSBr <- lsmeans(Br.aov, pairwise ~ Asrc)
LSBr


BR.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Br_num ~ Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

Brl.boot <- boot(allplants, BR.BOOT, R = 5000)
boot.ci(Brl.boot, index = 1, conf = .95, type = "bca")
boot.ci(Brl.boot, index = 2, conf = .95, type = "bca")
boot.ci(Brl.boot, index = 3, conf = .95, type = "bca")
boot.ci(Brl.boot, index = 4, conf = .95, type = "bca")
boot.ci(Brl.boot, index = 5, conf = .95, type = "bca")

# leaf number
LFN.aov <- aov(Lfnum~Asrc + Plt, data = hybplants)
anova(LFN.aov)
LSLFN <- lsmeans(LFN.aov, pairwise ~ Asrc)
LSLFN


LF.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Lfnum~Asrc+Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

LFN.boot <- boot(allplants, LF.BOOT, R = 5000)
boot.ci(LFN.boot, index = 1, conf = .95, type = "bca")
boot.ci(LFN.boot, index = 2, conf = .95, type = "bca")
boot.ci(LFN.boot, index = 3, conf = .95, type = "bca")
boot.ci(LFN.boot, index = 4, conf = .95, type = "bca")
boot.ci(LFN.boot, index = 5, conf = .95, type = "bca")

#Head total
TH.aov <- aov(Tothds~Asrc + Plt, data = hybplants)
anova(TH.aov)
LSTH <- lsmeans(TH.aov, pairwise ~ Asrc)
LSTH


TH.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Tothds~Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

THL.boot <- boot(allplants, TH.BOOT, R = 5000)
boot.ci(THL.boot, index = 1, conf = .95, type = "bca")
boot.ci(THL.boot, index = 2, conf = .95, type = "bca")
boot.ci(THL.boot, index = 3, conf = .95, type = "bca")
boot.ci(THL.boot, index = 4, conf = .95, type = "bca")
boot.ci(THL.boot, index = 5, conf = .95, type = "bca")

#Seed total
TS.aov <- aov(Tot_sds~Asrc + Plt, data = no_headless_hyb)
anova(TS.aov)
LSTS <- lsmeans(TS.aov, pairwise ~ Asrc)
LSTS


TS.BOOT <- function(no_headless_all, index) {
  AP.sub <- no_headless_all[index,]
  model.aov <- aov(Tot_sds ~ Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

TSL.boot <- boot(no_headless_all, TS.BOOT, R = 5000)
boot.ci(TSL.boot, index = 1, conf = .95, type = "bca")
boot.ci(TSL.boot, index = 2, conf = .95, type = "bca")
boot.ci(TSL.boot, index = 3, conf = .95, type = "bca")
boot.ci(TSL.boot, index = 4, conf = .95, type = "bca")
boot.ci(TSL.boot, index = 5, conf = .95, type = "bca")

#Mass per seed
MPS.aov <- aov(masspsd~Asrc + Plt, data =allsds, na.action = "na.exclude")
anova(MPS.aov)
LSMPS <- lsmeans(MPS.aov, pairwise ~ Asrc)
LSMPS



MPS.BOOT <- function(hybsds, index) {
  AP.sub <-hybsds[index,]
  model.aov <- aov(masspsd~Asrc + Plt, data =AP.sub, na.action = "na.exclude")
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

MPSL.boot <- boot(hybsds, MPS.BOOT, R = 5000)
boot.ci(MPSL.boot, index = 1, conf = .95, type = "bca")
boot.ci(MPSL.boot, index = 2, conf = .95, type = "bca")
boot.ci(MPSL.boot, index = 3, conf = .95, type = "bca")

#Seed area
ASA.aov <- aov(avgsdar~Asrc + Plt, data = hybsdar)
anova(ASA.aov)
LSASA <- lsmeans(ASA.aov, pairwise ~ Asrc)
LSASA


ASA.BOOT <- function(hybsdsar, index) {
  AP.sub <- hybsdsar[index,]
  model.aov <- aov(avgsdar~Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

ASAL.boot <- boot(hybsdar, ASA.BOOT, R = 5000)
boot.ci(ASAL.boot, index = 1, conf = .95, type = "bca")
boot.ci(ASAL.boot, index = 2, conf = .95, type = "bca")
boot.ci(ASAL.boot, index = 3, conf = .95, type = "bca")


#stem diameter
STD.aov <- aov(Stmdiam~Asrc + Plt, data = hybplants)
anova(STD.aov)
LSSTD <- lsmeans(STD.aov, pairwise ~ Asrc)
LSSTD


STD.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Stmdiam~Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

STDL.boot <- boot(allplants, STD.BOOT, R = 5000)
boot.ci(STDL.boot, index = 1, conf = .95, type = "bca")
boot.ci(STDL.boot, index = 2, conf = .95, type = "bca")
boot.ci(STDL.boot, index = 3, conf = .95, type = "bca")
boot.ci(STDL.boot, index = 4, conf = .95, type = "bca")
boot.ci(STDL.boot, index = 5, conf = .95, type = "bca")

#disk diameter 
DD.aov <- aov(Dskdiam~Asrc + Plt, data = hybplants)
anova(DD.aov)
LSDD <- lsmeans(DD.aov, pairwise ~ Asrc)
LSDD


DD.BOOT <- function(allplants, index) {
  AP.sub <- allplants[index,]
  model.aov <- aov(Dskdiam~Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

DDL.boot <- boot(allplants, DD.BOOT, R = 5000)
boot.ci(DDL.boot, index = 1, conf = .95, type = "bca")
boot.ci(DDL.boot, index = 2, conf = .95, type = "bca")
boot.ci(DDL.boot, index = 3, conf = .95, type = "bca")
boot.ci(DDL.boot, index = 4, conf = .95, type = "bca")
boot.ci(DDL.boot, index = 5, conf = .95, type = "bca")

#leaf ratio - 
LFR.aov <- aov(lfr~Asrc + Plt, data = hybplants)
anova(LFR.aov)
LSLFR <- lsmeans(LFR.aov, pairwise ~ Asrc)
LSLFR


LFR.BOOT <- function(hyblf, index) {
  AP.sub <- hyblf[index,]
  model.aov <- aov(lfr~Asrc + Plt, data = AP.sub, na.action = "na.exclude")
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

LFRL.boot <- boot(hyblf, LFR.BOOT, R = 5000)
boot.ci(LFRL.boot, index = 1, conf = .95, type = "bca")
boot.ci(LFRL.boot, index = 2, conf = .95, type = "bca")
boot.ci(LFRL.boot, index = 3, conf = .95, type = "bca")


stml.aov <- aov(Stmlen~Asrc + Plt, data = hybplants)
anova(stml.aov)
LSslfl <- lsmeans(stml.aov, pairwise ~ Asrc)
LSslfl

stm.BOOT <- function(hybstm, index) {
  AP.sub <- hybstm[index,]
  model.aov <- aov(Stmlen~Asrc + Plt, data = AP.sub)
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

stmL.boot <- boot(hybstm, stm.BOOT, R = 5000)
boot.ci(stmL.boot, index = 1, conf = .95, type = "bca")
boot.ci(stmL.boot, index = 2, conf = .95, type = "bca")
boot.ci(stmL.boot, index = 3, conf = .95, type = "bca")

gd.aov <- aov(sd_damage~Asrc + Plt, data = no_headless_hyb)
anova(gd.aov)
gdl <- lsmeans(gd.aov, pairwise ~ Asrc)
gdl


gd.BOOT <- function(no_headless_hyb, index) {
  AP.sub <- no_headless_hyb[index,]
  model.aov <- aov(sd_damage~Asrc + Plt, data = AP.sub, na.action = "na.exclude")
  lsbt <- lsmeans(model.aov, pairwise ~Asrc)
  lsbt2 <- summary(lsbt$lsmeans)
  lsbt3 <- c(lsbt2[c("lsmean")])
  lsbt4 <- lsbt3[[1]]
  lsbt4
}

gdL.boot <- boot(no_headless_hyb, gd.BOOT, R = 5000)
boot.ci(gdL.boot, index = 1, conf = .95, type = "bca")
boot.ci(gdL.boot, index = 2, conf = .95, type = "bca")
boot.ci(gdL.boot, index = 3, conf = .95, type = "bca")
