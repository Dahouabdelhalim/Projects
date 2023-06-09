library(ggplot2)
library(rcompanion)
library(xlsx)
library(effsize)
library(lme4)
library(lmerTest)
library(MuMIn)
library(rptR)

setwd("C:/Users/valer/Google Drive/1. ONGOING PAPERS/4. zebrafish individual risk/submission/AmeNat/Data_and_R_code")


# Load Ethovision tracking data
d0 <- read.xlsx(file="data_submission_AmeNat.xlsx", sheetName = "Rep")

d0$LINE <- as.factor(d0$LINE)
d0$REP <- as.factor(d0$REP)
d0$Group_Rep <- as.factor(d0$Group_Rep)
d0$Age <- as.factor(d0$Age)

# tranform data
d0$IID.BL_tuk = transformTukey(d0$IID.BL, plotit=FALSE)
plotNormalHistogram(d0$IID.BL_tuk)

# calculate Adjusted Repeatability
IID.rpt <- rpt(IID.BL_tuk ~ LINE + (1 | REP) + (1 | Group_Rep),
               grname = c("REP","Group_Rep"), data = d0,
               datatype = "Gaussian", nboot = 1000, npermut = 0)

IID.rpt
plot(IID.rpt$mod)


# load group-level data
d1 <- read.xlsx(file="data_submission_AmeNat.xlsx", sheetName = "group_data")

# order levels of factor Line
d1$LINE <- factor(d1$LINE,levels=c("RH","LH","SH"))

# mean values IID
aggregate(IID_bl~LINE, d1, FUN=mean)
aggregate(IID_bl~LINE, d1, FUN=sd)

# mean values NND
aggregate(NND_bl~LINE, d1, FUN=mean)
aggregate(NND_bl~LINE, d1, FUN=sd)

# mean total values Polarization
mean(d1$polarization)
sd(d1$polarization)


# tranform IID and NND to meet model assumption
d1$IID_bl_tuk = transformTukey(d1$IID_bl, plotit=FALSE)
plotNormalHistogram(d1$IID_bl_tuk)
d1$IID_tuk = transformTukey(d1$IID, plotit=FALSE)
plotNormalHistogram(d1$IID_tuk)


d1$NND_bl_tuk = transformTukey(d1$NND_bl, plotit=FALSE)
plotNormalHistogram(d1$NND_bl_tuk)
d1$NND_tuk = transformTukey(d1$NND, plotit=FALSE)
plotNormalHistogram(d1$NND_tuk)


d1$pol_tuk = transformTukey(d1$polarization, plotit=FALSE)
plotNormalHistogram(d1$pol_tuk)


# linear mixed effect model for NND & IID & pol
# repeat with the following response variables
# IID_bl_tuk; IID_tuk; NND_bl_tuk; NND_tuk; pol_tuk 

mod <- lmer(NND_bl_tuk ~ LINE + (1 |REP), data = d1, REML=F)
mod_null <- lmer(NND_bl_tuk ~ 1 + (1 |REP), data = d1, REML=F)

plot(mod)
qqnorm(residuals(mod))
hist(residuals(mod))

summary(mod)
anova(mod,mod_null)

r.squaredGLMM(mod)
rm(mod, mod_null)

confint.merMod(mod)


# correlation between NND and IID
mod2 <- lmer(IID_bl_tuk~NND_bl_tuk + (1 |REP), d1,  REML=F)
plot(mod2)
summary(mod2)
r.squaredGLMM(mod2)
anova(mod2)



# load risk-taking and NND data
d2 <- read.xlsx(file="data_submission_AmeNat.xlsx", sheetName = "corr")

d2_NND <- subset(d2, measure=="NND")

mod <- lmer(surf.mean.time.7~BL + (1 |REP), d2_NND,  REML=F)

summary(mod)
anova(mod)
r.squaredGLMM(mod)


# load individual-level data
d3 <- read.xlsx(file="data_submission_AmeNat.xlsx", sheetName = "individual_data")

d3$LINE <- factor(d3$LINE,levels=c("RH","LH","SH"))
d3$REP <- as.factor(d3$REP)

# tranform IID and NND to meet model assumption
d3$Real_SL_tuk = transformTukey(d3$Real_SL, plotit=FALSE)
plotNormalHistogram(d3$Real_SL_tuk)

# subset the data to remove significant differences in size
d4 <- subset(d3, Real_SL > 2.2 & Real_SL < 2.6)


# effect of selection line on body size (full dataset)
mod1 <- lmer(Real_SL_tuk ~ LINE + (1 | REP), data = d3, REML=F)
mod2 <- lmer(Real_SL_tuk ~ 1 + (1 | REP), data = d3, REML=F)

anova(mod1,mod2)
rm(mod1, mod2)


# effect of selection line on body size (size matched dataset)
mod3 <- lmer(Real_SL_tuk ~ LINE + (1 |REP), data = d4, REML=F)
mod4 <- lmer(Real_SL_tuk ~ 1 + (1 |REP), data = d4, REML=F)
anova(mod3,mod4)
rm(mod3, mod4)


# effect of selection line on burst rate (full data set)
mod <- lmer(burst_rate_s ~ LINE + (1 |REP), data = d3, REML=F)
mod_null <- lmer(burst_rate_s ~ 1 + (1 |REP), data = d3, REML=F)

plot(mod)
qqnorm(residuals(mod))

summary(mod)
r.squaredGLMM(mod)

anova(mod,mod_null)
rm(mod, mod_null)


# effect of selection line on speed (full data set)
mod <- lmer(speed_bl_s ~ LINE + (1 |REP), data = d3, REML=F)
mod_null <- lmer(speed_bl_s ~ 1 + (1 |REP), data = d3, REML=F)

plot(mod)
qqnorm(residuals(mod))

summary(mod)
r.squaredGLMM(mod)

anova(mod,mod_null)
rm(mod, mod_null)


# effect of selection line on burst rate (size matched dataset)
mod <- lmer(burst_rate_s ~ Real_SL_tuk + LINE + (1 |REP), data = d4, REML=F)
mod1 <- lmer(burst_rate_s ~ LINE + (1 |REP), data = d4, REML=F)
mod2 <- lmer(burst_rate_s ~ Real_SL_tuk + (1 |REP), data = d4, REML=F)
mod3 <- lmer(burst_rate_s ~ 1 + (1 |REP), data = d4, REML=F)

plot(mod)
qqnorm(residuals(mod))

summary(mod)
r.squaredGLMM(mod)

# effect of LINE
anova(mod,mod2)

# effect of body size
anova(mod,mod1)

rm(mod, mod1, mod2, mod3)


# effect of selection line on speed (size matched dataset)
mod <- lmer(speed_bl_s ~ Real_SL_tuk + LINE + (1 |REP), data = d4, REML=F)
mod1 <- lmer(speed_bl_s ~ LINE + (1 |REP), data = d4, REML=F)
mod2 <- lmer(speed_bl_s ~ Real_SL_tuk + (1 |REP), data = d4, REML=F)
mod3 <- lmer(speed_bl_s ~ 1 + (1 |REP), data = d4, REML=F)

plot(mod)
qqnorm(residuals(mod))

summary(mod)
r.squaredGLMM(mod)

# effect of LINE
anova(mod,mod2)

# effect of body size
anova(mod,mod1)



# correlation burst rate vs speed
mod <- lmer(speed_bl_s~burst_rate_s + (1 | REP), d3, REML=F)
plot(mod)
summary(mod)
r.squaredGLMM(mod)
anova(mod)