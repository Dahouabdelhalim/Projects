setwd('') # set the working directory


df.toad <- read.delim("dataframe_Duttaphrynus_melanostictus.csv", sep =';')
df.toad$sex <- as.factor(df.toad$sex)
df.toad$sex <- relevel(df.toad$sex, ref = "1")# 0 = Males; 1 = Females
df.toad$year <- as.factor(df.toad$year)

# create summary statistics of morphological traits # Table S1

MM <- subset(toad, sex==0)
FF <- subset(toad, sex==1)

# males

svl.MM <- mean(MM$svl, na.rm=T)
sd.svl.MM <- sd(MM$svl, na.rm=T)

hw.MM <- mean(MM$hw, na.rm=T)
sd.hw.MM <- sd(MM$hw, na.rm=T)

rul.MM <- mean(MM$rul, na.rm=T)
sd.rul.MM <- sd(MM$rul, na.rm=T)

tbl.MM <- mean(MM$tbl, na.rm=T)
sd.tbl.MM <- sd(MM$tbl, na.rm=T)

mass.MM <- mean(MM$mass, na.rm=T)
sd.mass.MM <- sd(MM$mass, na.rm=T)

hw.ratio.MM <- mean(MM$hw.rel, na.rm=T)
sd.hw.ratio.MM <- sd(MM$hw.rel, na.rm=T)

rul.ratio.MM <- mean(MM$rul.rel, na.rm=T)
sd.rul.ratio.MM <- sd(MM$rul.rel, na.rm=T)

tbl.ratio.MM <- mean(MM$tbl.rel, na.rm=T)
sd.tbl.ratio.MM <- sd(MM$tbl.rel, na.rm=T)

#females

svl.FF <- mean(FF$svl, na.rm=T)
sd.svl.FF <- sd(FF$svl, na.rm=T)

hw.FF <- mean(FF$hw, na.rm=T)
sd.hw.FF <- sd(FF$hw, na.rm=T)

rul.FF <- mean(FF$rul, na.rm=T)
sd.rul.FF <- sd(FF$rul, na.rm=T)

tbl.FF <- mean(FF$tbl, na.rm=T)
sd.tbl.FF <- sd(FF$tbl, na.rm=T)

mass.FF <- mean(FF$mass, na.rm=T)
sd.mass.FF <- sd(FF$mass, na.rm=T)

hw.ratio.FF <- mean(FF$hw.rel, na.rm=T)
sd.hw.ratio.FF <- sd(FF$hw.rel, na.rm=T)

rul.ratio.FF <- mean(FF$rul.rel, na.rm=T)
sd.rul.ratio.FF <- sd(FF$rul.rel, na.rm=T)

tbl.ratio.FF <- mean(FF$tbl.rel, na.rm=T)
sd.tbl.ratio.FF <- sd(FF$tbl.rel, na.rm=T)

### calculate sexual dimorphism index (SDI)

SDI <- (table_summary$`Females mean`[1]/table_summary$`Males mean`[1])-1

### Polynomial regressions (1st-4th order) to check for allomorphic relations between morphological traits (HW, TBL and RUL) and SVL, including the effect of sex. 

# TBL

tbl1 <- lm(log(tbl) ~ log(svl) + sex, data = df.toad)
tbl2 <- lm(log(tbl)~poly(log(svl),2,raw=TRUE)+sex, data=df.toad)
tbl3 <- lm(log(tbl)~poly(log(svl),3,raw=TRUE)+sex, data=df.toad)
tbl4 <- lm(log(tbl)~poly(log(svl),4,raw=TRUE)+sex, data=df.toad)

AIC(tbl1)
summary(tbl1)

AIC(tbl2)
summary(tbl2)

AIC(tbl3)
summary(tbl3)

AIC(tbl4)
summary(tbl4)

#HW

hw1 <- lm(log(hw) ~ log(svl) + sex, data = df.toad)
hw2 <- lm(log(hw)~poly(log(svl),2,raw=TRUE) + sex, data=df.toad)
hw3 <- lm(log(hw)~poly(log(svl),3,raw=TRUE) + sex, data=df.toad)
hw4 <- lm(log(hw)~poly(log(svl),4,raw=TRUE) + sex, data=df.toad)

AIC(hw1)
summary(hw1)

AIC(hw2)
summary(hw2)

AIC(hw3)
summary(hw3)

AIC(hw4)
summary(hw4)

# RUL
# create subset for RUL to exclude NAs

df.RUL <- subset(df.toad, !is.na(rul))

rul1 <- lm(log(rul) ~ log(svl) + sex, data = df.RUL)
rul2 <- lm(log(rul) ~ poly(log(svl),2,raw=TRUE) + sex, data = df.RUL)
rul3 <- lm(log(rul) ~ poly(log(svl),3,raw=TRUE) + sex, data = df.RUL)
rul4 <- lm(log(rul) ~ poly(log(svl),4,raw=TRUE) + sex, data = df.RUL)

AIC(rul1)
summary(rul1)

AIC(rul2)
summary(rul2)

AIC(rul3)
summary(rul3)

AIC(rul4)
summary(rul4)

###### linear mixed-effect models (LMMs) to test the spatial sorting hypothesis, by assessing relationships between morphometric features (TBL, HW, RUL) and distance from the introduction point

library(lme4)
library(lmerTest)

# TBL, exclude NAs
df.TBL <- subset(df.toad, !is.na(tbl))

m.tbl <- lmer(log(tbl) ~ log(svl) + sex + log(distance) + (1|year), data = df.TBL)
summary(m.tbl)

# Moran's I test to assess spatial autocorrelation of the residuals of the model

library(EcoGenetics)

# extract residuals model

df.TBL$res <- residuals(m.tbl)

# smax = 7200; minimum distance to connect all data points
mor.tbl <- eco.correlog(Z=df.TBL$res, XY = df.TBL[,cbind(5:6)], method = "I", smax=7200, nsim=1000) 

mor.tbl

# eco.plotCorrelog(). To visualize the Moran's I test by classes of distance

# eco.plotCorrelog(mor.tbl)

### RUL

m.rul <- lmer(log(rul) ~ log(svl) + sex + log(distance) + (1|year), data = df.RUL)

summary(m.rul)

# spatial autocorrelation RUL

df.RUL$res <- residuals(m.rul)

mor.rul <- eco.correlog(Z=df.RUL$res, XY = df.RUL[,cbind(5:6)], method = "I", smax= 7200, nsim=1000)
mor.rul

# eco.plotCorrelog(mor.rul)

### hw

df.HW <- subset(df.toad, !is.na(df.toad$hw))
  
m.hw <- lmer(log(hw) ~ log(svl) + sex + log(distance) + (1|year), data = df.HW)

summary(m.hw)

# spatial autocorrelation HW

df.HW$res <- residuals(m.hw)

mor.hw <- eco.correlog(Z=df.HW$res, XY = df.HW[,cbind(5:6)],
                      method = "I", smax=7200, nsim=1000)
mor.hw

# eco.plotCorrelog(mor.hw)

### SVL

df.SVL <- df.toad

m.svl <- lmer(log(svl) ~ sex + log(distance) + (1|year), data = df.SVL)

summary(m.svl)

# spatial autocorrelation SVL

df.SVL$res <- residuals(m.svl)

mor.svl <- eco.correlog(Z=df.SVL$res, XY = df.SVL[,cbind(5:6)], method = "I", smax=7200, nsim=1000)

mor.svl

# eco.plotCorrelog(mor.svl)

### RI 

df.RI <- subset(toad, !is.na(toad$RI))

m.RI <- lmer(RI ~ log(svl) + sex + poly(log(jd), 2) + log(distance) + (1|year), data = df.RI)

summary(m.RI)

#Spatial autocorrelation RI
df.RI$res <- residuals(m.RI)

mor.smi <- eco.correlog(Z=df.RI$res, XY = df.RI[,cbind(5:6)], method = "I", smax= 7200, nsim=1000)

mor.smi

# eco.plotCorrelog(mor.smi)