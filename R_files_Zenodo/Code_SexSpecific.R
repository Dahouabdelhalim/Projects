# Code_SexSpecific
# PACKAGES ----
library(ggplot2)
library(lme4) # for lmer
library(lmerTest) # for least square means
library("openxlsx")
library(readxl)
library(readODS)
library(car) # this is used for the levenes test as well as p values of glmers
library(lsmeans)
library(anytime)
library(glmmTMB)
library(multcomp)
library(optimx) # used from optimising lmers when models dont converge using LM BFGS instread of bobbqa
library(rcompanion) # for fct "plotNormalHistogram"
library(patchwork)
library(ggpubr)
library(multcomp)
library(emmeans)
library(tidyr)
library(dplyr)
library(reshape2)
library(insight)
library(blmeco)
library(DHARMa)
library(plyr)
library(ggridges)
library(metaplot)
library(rptR)
library(sjPlot)
library(effects)
library(xtable)

# a word of warning, come from the DHARMa vignette:
# "A word of warning that applies also to all other tests that follow: significance in hypothesis tests depends on at least 2 ingredients: strenght of the signal, and number of data points. Hence, the p-value alone is not a good indicator of the extent to which your residuals deviate from assumptions. Specifically, if you have a lot of data points, residual diagnostics will nearly inevitably become significant, because having a perfectly fitting model is very unlikely. That, however, doesnâ€™t necessarily mean that you need to change your model. The p-values confirm that there is a deviation from your null hypothesis. It is, however, in your discretion to decide whether this deviation is worth worrying about. If you see a dispersion parameter of 1.01, I would not worry, even if the test is significant. A significant value of 5, however, is clearly a reason to move to a model that accounts for overdispersion."

# theme for ggplot2

my_theme <- function(...) {
  theme(text=element_text(size=20), panel.grid.major = element_blank(), 
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "none")
}

# DATAWRANGLING ----

# info
# thorax length was only measured in the F1
# development time was NOT monitored for the progeny of females, we work with paternal lines!

# import the dataset
# F0
F0_comp <- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F0.xlsx", sheet = "F0_comp")
F0_f.single <- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F0.xlsx", sheet = "F0_f_single_progeny")
F0_m.single<- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F0.xlsx", sheet = "F0_m_single_progeny")

# F2
F1_comp <- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F1.xlsx", sheet = "F1_comp")
F1_f.single <- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F1.xlsx", sheet = "F1_f_single_progeny")
F1_m.single <- read_xlsx("/Users/Valerian/Dropbox/1_Sex_specific/20220321_F1.xlsx", sheet = "F1_m_single_progeny")

# change class
F0_comp <- F0_comp %>%
  dplyr::mutate_at(vars(treatment, vial,isoline, id, repl), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment)) %>%
  mutate(thoraxF = thoraxF/4.2) %>%
  mutate(thoraxM = thoraxM/4.2) %>%
  mutate(centred_female_TL = thoraxF-mean(thoraxF, na.rm = T)) %>%
  mutate(centred_male_TL = thoraxM-mean(thoraxM, na.rm = T))

F0_f.single <- F0_f.single %>%
  dplyr::mutate_at(vars(isoline, treatment, id, repl), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment)) %>%
  mutate(thoraxF = thoraxF/4.2) %>%
  mutate(centred_female_TL = thoraxF-mean(thoraxF, na.rm = T))

F0_m.single <- F0_m.single %>%
  dplyr::mutate_at(vars(treatment, vial, isoline, id, repl), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment)) %>%
  mutate(thoraxM = thoraxM/4.2) %>%
  mutate(centred_male_TL = thoraxM-mean(thoraxM, na.rm = T)) 

F1_comp <- F1_comp %>%
  dplyr::mutate_at(vars(treatment,isoline, id, repl, paternal, filial), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment)) %>%
  mutate(centred_female_TL = thoraxF-mean(thoraxF, na.rm = T)) %>%
  mutate(centred_male_TL = thoraxM-mean(thoraxM, na.rm = T)) 

F1_f.single <- F1_f.single %>%
  dplyr::mutate_at(vars(isoline, treatment, id, paternal, filial, repl), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment), 
                paternalNum = as.numeric(paternal),
                filialNum = as.numeric(filial)) %>%
  mutate(centred_female_TL = thoraxF-mean(thoraxF, na.rm = T))

F1_m.single <- F1_m.single %>%
  dplyr::mutate_at(vars(isoline, treatment, id, paternal, filial, repl, vial), funs(factor)) %>%
  dplyr::mutate(treatmentNum = as.numeric(treatment), 
                paternalNum = as.numeric(paternal),
                filialNum = as.numeric(filial)) %>%
  mutate(centred_male_TL = thoraxM-mean(thoraxM, na.rm = T))

# CONVERTING to z-scores ----
# converting all continuous and ordinal fixed effects to zscores 

# F0_comp, z-score
zscoreF0_comp <- F0_comp %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscoreRemateDay = scale(remateDay),
         zscoreCopDuration = scale(copDuration),
         zscoreG = scale(G),
         zscoreR = scale(R),
         zscoreThoraxM = scale(thoraxM),
         zscoreThoraxF = scale(thoraxF),
         # zscoreTotprog = scale(G+R),
         zscoreViability = scale(viability),
         zscoreDevTimeF = scale(devTimeF),
         zscoreDevTimeM = scale(devTimeM)) %>%
  as.data.frame()

# F0_f.single, z-score
zscoreF0_f.single <- F0_f.single %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscoreFemales = scale(females),
         zscoreMales = scale(males),
         zscoreTotProg = scale(totProg),
         zscoreThoraxF = scale(thoraxF)) %>%
  as.data.frame()

# F0_m.single, z-score
zscoreF0_m.single <- F0_m.single %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscorematingLatency = scale(matingLatency),
         zscoreCopDuration = scale(copDuration),
         zscoreFemales = scale(females),
         zscoreMales = scale(males),
         zscoreTotProg = scale(totProg),
         zscorePropMale = scale(propMale),
         zscoreViability = scale(viability),
         zscoreThoraxM = scale(thoraxM),
         zscoreDevTimeM = scale(devTimeM)) %>%
  as.data.frame()

# F1_comp, z-score
zscoreF1_comp <- F1_comp %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscorePaternal = scale(as.numeric(paternal)),
         zscoreFilial = scale(as.numeric(filial)),
         zscoreRemateDay = scale(remateDay),
         zscoreCopDuration = scale(copDuration),
         zscoreG = scale(G),
         zscoreR = scale(R),
         zscoreThoraxM = scale(thoraxM),
         zscoreThoraxF = scale(thoraxF),
         zscoreViability = scale(viability),
         zscoreDevTimeF = scale(devTimeF),
         zscoreDevTimeM = scale(devTimeM)) %>%
  as.data.frame()

# F1_f.single, z-score
zscoreF1_f.single <- F1_f.single %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscoreFemales = scale(females),
         zscoreMales = scale(males),
         zscoreTotProg = scale(totProg),
         zscorePaternal = scale(paternalNum),
         zscoreFilial = scale(filialNum),
         zscoreDevTimeF = scale(devTimeF),
         zscoreThoraxF = scale(thoraxF)) %>%
  as.data.frame()

# F1_m.single, z-score
zscoreF1_m.single <- F1_m.single %>%
  group_by(treatment, isoline) %>%
  mutate(zscoreDiet = scale(treatmentNum),
         zscorematingLatency = scale(matingLatency),
         zscoreCopDuration = scale(copDuration),
         zscoreFemales = scale(females),
         zscoreMales = scale(males),
         zscoreTotProg = scale(totProg),
         zscorePaternal = scale(paternalNum),
         zscoreFilial = scale(filialNum),
         zscoreThoraxM = scale(thoraxM),
         zscoreDevTimeM = scale(devTimeM)) %>%
  as.data.frame()

# ANALYSIS ----

# _______________________----
# ***F0*** ----

# ***thorax length ----
# thorax length female ----
F0_zscoreTho_length <- zscoreF0_comp

# looking at distribution
# female
ggplot(F0_zscoreTho_length, aes(x = thoraxF, y = treatment))+
  geom_density_ridges(aes(fill= treatment))
# male
ggplot(F0_zscoreTho_length, aes(x = thoraxM, y = treatment))+
  geom_density_ridges(aes(fill= treatment))

# w/o duplicated values 
mX <- lmer(thoraxF ~ zscoreDiet + (1|isoline), data = F0_zscoreTho_length)
plot(simulateResiduals(fittedModel = mX))  
summary(mX)

reptX <- rpt(thoraxF ~ zscoreDiet + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

reptXa  <- rpt(thoraxF ~ zscoreDiet + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F0_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)


print(reptX)
print(reptXa)

#note that variance components vary slightly each time rpt is run

#slope variance
slX <- (reptX$R$isoline - reptXa$R$isoline) / (reptX$R$isoline + reptX$R$Residual + reptX$R$Fixed)*100

#isoline variance
isoX <- (reptXa$R$isoline) / (reptXa$R$isoline + reptXa$R$Residual + reptXa$R$Fixed)*100

# thorax length male ----
# w/o duplicated values 
mY <- lmer(thoraxM ~ zscoreDiet + (1|isoline), data = F0_zscoreTho_length)
plot(simulateResiduals(fittedModel = mY))  
summary(mY)

reptY <- rpt(thoraxM ~ zscoreDiet + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

reptYa  <- rpt(thoraxM ~ zscoreDiet + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F0_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)


print(reptY)
print(reptYa)

#note that variance components vary slightly each time rpt is run

#slope variance
slY <- (reptY$R$isoline - reptYa$R$isoline) / (reptY$R$isoline + reptY$R$Residual + reptY$R$Fixed)*100

#isoline variance
isoY <- (reptYa$R$isoline) / (reptYa$R$isoline + reptYa$R$Residual + reptYa$R$Fixed)*100

# viability ----
# import viability
F0_zscoreVia <- zscoreF0_comp # there were 40 larvae in the F0

# prepare the dataset w/o replicated values
F0_zscoreVia2 <- F0_zscoreVia[!duplicated(F0_zscoreVia[ , c("treatment","isoline")]),]

# model viability
# looking at the distribution
ggplot(F0_zscoreVia2, aes(x=viability, y=treatment))+
  geom_density_ridges(aes(fill=treatment)) # looks normally distributed

# if any, change the vial that have a viability >1 to 40. 
nrow(subset(F0_zscoreVia2,viability>1)) # 0 row

# get the number of dead and alive
F0_zscoreVia2$alive <- F0_zscoreVia2$viability*40
F0_zscoreVia2$dead <- 40-F0_zscoreVia2$alive

# change the class so it works for the binomial model
F0_zscoreVia2$aliveInt <- as.integer(F0_zscoreVia2$alive)
F0_zscoreVia2$deadInt <- as.integer(F0_zscoreVia2$dead)

m1 <- glmer(cbind(aliveInt,deadInt) ~ zscoreDiet + (1|isoline), data = F0_zscoreVia2, family = "binomial")
plot(simulateResiduals(fittedModel = m1)) 
summary(m1)

rept1 <- rptBinary(cbind(aliveInt,deadInt) ~ zscoreDiet + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
                   data = F0_zscoreVia2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept1a <- rptBinary(cbind(aliveInt,deadInt) ~ zscoreDiet + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
                    data = F0_zscoreVia2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept1)
print(rept1a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl1 <- (rept1$R$isoline[2] - rept1a$R$isoline[2]) / (rept1$R$isoline[2] + rept1$R$Residual[2] + rept1$R$Fixed[2])*100

#isoline variance 
iso1 <- (rept1a$R$isoline[2]) / (rept1a$R$isoline[2] + rept1a$R$Residual[2] + rept1a$R$Fixed[2])*100

# ***devTime ----
# import development
F0_zscoreDev <- zscoreF0_comp

# devTime female ----
# look at distribution
ggplot(F0_zscoreDev,aes(x = zscoreDevTimeF, y = treatment))+
  geom_density_ridges(aes(fill = zscoreDiet))

m2 <- lmer(devTimeF ~ zscoreDiet + zscoreThoraxF + (1|isoline), F0_zscoreDev, na.action = na.omit)
plot(simulateResiduals(fittedModel = m2)) 
summary(m2)

rept2 <- rpt(devTimeF ~ zscoreDiet + zscoreThoraxF + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept2a <- rpt(devTimeF ~ zscoreDiet + zscoreThoraxF + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)
print(rept2)
print(rept2a)

#slope variance
sl2 <- (rept2$R$isoline - rept2a$R$isoline) / (rept2$R$isoline + rept2$R$Residual + rept2$R$Fixed)*100

# sl3 <- (rept3$R$isoline - rept3a$R$isoline) / (rept3$R$isoline + rept3$R$Residual + rept3$R$Fixed)*100

#approaching 0

#isoline variance 
iso2 <- (rept2a$R$isoline) / (rept2a$R$isoline + rept2a$R$Residual + rept2a$R$Fixed)*100

# devTime male ----
# look at distribution
ggplot(F0_zscoreDev, aes(x = zscoreDevTimeM, y = treatment))+
  geom_density_ridges(aes(fill = zscoreDiet))

m3 <- lmer(devTimeM ~ zscoreDiet +  zscoreThoraxM + (1|isoline), F0_zscoreDev, na.action = na.omit)
plot(simulateResiduals(fittedModel = m3)) 
summary(m3)

rept3 <- rpt(devTimeM ~ zscoreDiet + zscoreThoraxM +(zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept3a <- rpt(devTimeM ~ zscoreDiet + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept3)
print(rept3a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl3 <- (rept3$R$isoline - rept3a$R$isoline) / (rept3$R$isoline + rept3$R$Residual + rept3$R$Fixed)*100
#approaching 0

#isoline variance 
iso3 <- (rept3a$R$isoline) / (rept3a$R$isoline + rept3a$R$Residual + rept3a$R$Fixed)*100

# ***fecundity ----
# import fecundity
# female
F0_zscoreProg.f <- zscoreF0_f.single
# male
F0_zscoreProg.m <- zscoreF0_m.single

# remove the totProg outliers (we only take into account vials that had more than 5 offspring)
F0_zscoreProg.f <- subset(F0_zscoreProg.f,totProg>=5)
F0_zscoreProg.m <- subset(F0_zscoreProg.m,totProg>=5)

# fecundity female ----
# check distribution
hist(F0_zscoreProg.f$totProg) 
hist(sqrt(F0_zscoreProg.f$totProg))

m4 <- lmer(sqrt(totProg) ~ zscoreDiet + zscoreThoraxF + (1|isoline) , F0_zscoreProg.f)
plot(simulateResiduals(fittedModel = m4))
summary(m4)

rept4 <- rpt(totProg ~ zscoreDiet + zscoreThoraxF + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreProg.f, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept4a <- rpt(totProg ~ zscoreDiet + zscoreThoraxF + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreProg.f, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept4)
print(rept4a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl4 <- (rept4$R$isoline - rept4a$R$isoline) / (rept4$R$isoline + rept4$R$Residual + rept4$R$Fixed)*100

#isoline variance 
iso4 <- (rept4a$R$isoline) / (rept4a$R$isoline + rept4a$R$Residual + rept4a$R$Fixed)*100

# fecundity male ----
# check distribution
hist(F0_zscoreProg.m$totProg)
hist(sqrt(F0_zscoreProg.m$totProg))

m5 <- lmer(sqrt(totProg) ~ zscoreDiet + zscoreThoraxM + (1|isoline) , F0_zscoreProg.m)
plot(simulateResiduals(fittedModel = m5))
summary(m5)

rept5 <- rpt(totProg ~ zscoreDiet + zscoreThoraxM + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreProg.m, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept5a <- rpt(totProg ~ zscoreDiet + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreProg.m, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept5)
print(rept5a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl5 <- (rept5$R$isoline - rept5a$R$isoline) / (rept5$R$isoline + rept5$R$Residual + rept5$R$Fixed)*100

#isoline variance 
iso5 <- (rept5a$R$isoline) / (rept5a$R$isoline + rept5a$R$Residual + rept5a$R$Fixed)*100

# P2 ----
# import P2
F0_zscoreP2 <- zscoreF0_comp

# remove no matings
F0_zscoreP2 <- subset(F0_zscoreP2,R>0)
# F0_zscoreP2 <- subset(F0_zscoreP2,G>0) # we keep P2=1 in the dataset

# look at distribution
hist(F0_zscoreP2$P2)
hist(asin(sqrt(F0_zscoreP2$P2))) # bit better 

# add OLRE
F0_zscoreP2$obs <- as.factor(1:dim(F0_zscoreP2)[1])

m6 <- glmer(cbind(R,G) ~ zscoreDiet + zscoreThoraxM + (1|isoline) + (1|obs), family=binomial(link = "logit") , data=F0_zscoreP2)
plot(simulateResiduals(fittedModel = m6))
summary(m6)

rept6 <- rptBinary(cbind(R,G) ~ zscoreDiet + zscoreThoraxM + (zscoreDiet|isoline) + (1|obs), grname = c("isoline", "obs", "Fixed", "Residual"), 
                   data = F0_zscoreP2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept6a <- rptBinary(cbind(R,G) ~ zscoreDiet + zscoreThoraxM + (1|isoline) + (1|obs), grname = c("isoline", "obs", "Fixed", "Residual"), 
                    data = F0_zscoreP2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept6)
print(rept6a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl6 <- (rept6$R$isoline[2] - rept6a$R$isoline[2]) / (rept6$R$isoline[2] + rept6$R$Residual[2] + rept6$R$Fixed[2])*100
#approaching 0

#isoline variance 
iso6 <- (rept6a$R$isoline[2]) / (rept6a$R$isoline[2] + rept6a$R$Residual[2] + rept6a$R$Fixed[2])*100

# mating duration ----
# import mating duration
F0_zscoreMat <- zscoreF0_m.single

# check distribution
hist(F0_zscoreMat$copDuration)
hist(log(F0_zscoreMat$copDuration))

m7 <- lmer(copDuration ~ zscoreDiet + zscoreThoraxM + (1|isoline), F0_zscoreMat, na.action=na.omit)
plot(simulateResiduals(fittedModel = m7)) 
summary(m7)

rept7 <- rpt(copDuration ~ zscoreDiet + zscoreThoraxM + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreMat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept7a <- rpt(copDuration ~ zscoreDiet + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreMat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept7)
print(rept7a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl7 <- (rept7$R$isoline - rept7a$R$isoline) / (rept7$R$isoline + rept7$R$Residual + rept7$R$Fixed)*100

#isoline variance 
iso7 <- (rept7a$R$isoline) / (rept7a$R$isoline + rept7a$R$Residual + rept7a$R$Fixed)*100

# mating latency ----
# import mating latency
F0_zscoreLat <- zscoreF0_m.single

# look at distribution
hist(F0_zscoreLat$matingLatency) # no
hist(log(F0_zscoreLat$matingLatency)) # better

# add log transformed data
F0_zscoreLat$loggedmatingLatency <- log(F0_zscoreLat$matingLatency)

# replace "-Inf" values by "0"
F0_zscoreLat[which(F0_zscoreLat$loggedmatingLatency=="-Inf"),"loggedmatingLatency"] <- 0

# final model
m8 <- lmer(loggedmatingLatency ~ zscoreDiet + zscoreThoraxM + (1|isoline), F0_zscoreLat) 
plot(simulateResiduals(fittedModel = m8)) 
summary(m8)

rept8 <- rpt(loggedmatingLatency ~ zscoreDiet + zscoreThoraxM + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreLat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept8a <- rpt(loggedmatingLatency ~ zscoreDiet + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreLat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept8)
print(rept8a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl8 <- (rept8$R$isoline - rept8a$R$isoline) / (rept8$R$isoline + rept8$R$Residual + rept8$R$Fixed)*100

#isoline variance 
iso8 <- (rept8a$R$isoline) / (rept8a$R$isoline + rept8a$R$Residual + rept8a$R$Fixed)*100

# remating day ----
F0_zscoreRemat <- zscoreF0_comp

# check distribution
hist(F0_zscoreRemat$remateDay) 

# model
# add OLRE

m9 <- glmer(remateDay ~ zscoreDiet + zscoreThoraxM + (1|isoline), F0_zscoreRemat, na.action=na.omit, family = "poisson")
plot(simulateResiduals(fittedModel = m9))
summary(m9)

rept9 <- rpt(remateDay ~ zscoreDiet + zscoreThoraxM + (zscoreDiet|isoline), grname = c("isoline", "Fixed", "Residual"), 
             data = F0_zscoreRemat, datatype = "Poisson", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept9a <- rpt(remateDay ~ zscoreDiet + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F0_zscoreRemat, datatype = "Poisson", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept9)
print(rept9a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl9 <- (rept9$R$isoline[2] - rept9a$R$isoline[2]) / (rept9$R$isoline[2] + rept9$R$Residual[2] + rept9$R$Fixed[2])*100
# approximate 0

#isoline variance 
iso9 <- (rept9a$R$isoline[2]) / (rept9a$R$isoline[2] + rept9a$R$Residual[2] + rept9a$R$Fixed[2])*100
# approximate 0

# _______________________----
# ***F1*** ----
# ***thorax length ----

F1_zscoreTho_length <- zscoreF1_comp

# looking at distribution
# female
ggplot(F1_zscoreTho_length, aes(x = thoraxF, y = treatment))+
  geom_density_ridges(aes(fill= treatment))
# male
ggplot(F1_zscoreTho_length, aes(x = thoraxM, y = treatment))+
  geom_density_ridges(aes(fill= treatment))

# thorax length female ----
# w/o duplicated values 
m10 <- lmer(thoraxF ~ zscorePaternal * zscoreFilial + (1|isoline), data = F1_zscoreTho_length)
plot(simulateResiduals(fittedModel = m10))  
summary(m10)

rept10 <- rpt(thoraxF ~ zscorePaternal * zscoreFilial + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept10a <- rpt(thoraxF ~ zscorePaternal * zscoreFilial + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept10)
print(rept10a)

#note that variance components vary slightly each time rpt is run

#slope variance
(0-0) / (0 + 0 + 0.001)*100
sl10 <- (rept10$R$isoline - rept10a$R$isoline) / (rept10$R$isoline + rept10$R$Residual + rept10$R$Fixed)*100
#approaching 0

#isoline variance
0 / (0 + 0 + 0.001)*100
iso10 <- (rept10a$R$isoline) / (rept10a$R$isoline + rept10a$R$Residual + rept10a$R$Fixed)*100
#approaching 0

# thorax length male ----
# model male thorax length
# w/o duplicated values 
m11 <- lmer(thoraxM ~ zscorePaternal * zscoreFilial + (1|isoline), data = F1_zscoreTho_length)
plot(simulateResiduals(fittedModel = m11)) 
summary(m11)

rept11 <- rpt(thoraxM ~ zscorePaternal * zscoreFilial + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept11a <- rpt(thoraxM ~ zscorePaternal * zscoreFilial + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreTho_length, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept11)
print(rept11a)

#note that variance components vary slightly each time rpt is run

#slope variance
(0-0) / (0 + 0 + 0.001)*100
sl11 <- (rept11$R$isoline - rept11a$R$isoline) / (rept11$R$isoline + rept11$R$Residual + rept11$R$Fixed)*100
#approaching 0

#isoline variance
0 / (0 + 0 + 0.001)*100
iso11 <- (rept11a$R$isoline) / (rept11a$R$isoline + rept11a$R$Residual + rept11a$R$Fixed)*100
#approaching 0

# viability ----
# import viability
F1_zscoreVia <- zscoreF1_m.single

# looking at the distribution
ggplot(F1_zscoreVia,aes(x=viability, y=treatment))+
  geom_density_ridges(aes(fill=treatment))

# vials with viability >1 (because of error counting errors). 
nrow(subset(F1_zscoreVia, viability>1)) # 8

# Change these values to 1.
F1_zscoreVia[which(F1_zscoreVia$viability > 1), "viability"] <- 1 # is it actually ok to do this?
F1_zscoreVia[which(F1_zscoreVia$alive> 30), "alive"] <- 30 

# get the number of dead and alive
F1_zscoreVia$alive <- F1_zscoreVia$viability*30
F1_zscoreVia$dead <- 30-F1_zscoreVia$alive

# change the class so it works for the binomial model
F1_zscoreVia$aliveInt <- as.integer(F1_zscoreVia$alive)
F1_zscoreVia$deadInt <- as.integer(F1_zscoreVia$dead)

# prepare the dataset w/o replicated values
F1_zscoreVia2 <- F1_zscoreVia[!duplicated(F1_zscoreVia[ , c("treatment","isoline")]),]

# w/o duplicated values 
m12 <- glmer(cbind(aliveInt,deadInt) ~ zscorePaternal * zscoreFilial + (1|isoline), data = F1_zscoreVia2, family = "binomial")
plot(simulateResiduals(fittedModel = m12)) 
summary(m12)

rept12 <- rptBinary(cbind(aliveInt,deadInt) ~ zscorePaternal * zscoreFilial + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
                    data = F1_zscoreVia2,nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept12a <- rptBinary(cbind(aliveInt,deadInt) ~ zscorePaternal * zscoreFilial + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
                     data = F1_zscoreVia2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept12)
print(rept12a)

# note that variance components vary slightly each time rpt is run

#slope variance
sl12 <- (rept12$R$isoline[2] - rept12a$R$isoline[2]) / (rept12$R$isoline[2] + rept12$R$Residual[2] + rept12$R$Fixed[2])*100

#isoline variance
iso12 <- (rept12a$R$isoline[2]) / (rept12a$R$isoline[2] + rept12a$R$Residual[2] + rept12a$R$Fixed[2])*100

#*** devTime ----
# import development
F1_zscoreDev <- zscoreF1_comp

# devTime female ----
# look at distribution
ggplot(F1_zscoreDev,aes(x = zscoreDevTimeF, y = treatment))+
  geom_density_ridges(aes(fill = zscoreDiet))

m13 <- lmer(devTimeF ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (1|isoline), F1_zscoreDev, na.action = na.omit)
plot(simulateResiduals(fittedModel = m13)) 
summary(m13)

rept13 <- rpt(devTimeF ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept13a <- rpt(devTimeF ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)
print(rept13)
print(rept13a)

#slope variance
sl13 <- (rept13$R$isoline - rept13a$R$isoline) / (rept13$R$isoline + rept13$R$Residual + rept13$R$Fixed)*100
#approaching 0

#isoline variance 
iso13 <- (rept13a$R$isoline) / (rept13a$R$isoline + rept13a$R$Residual + rept13a$R$Fixed)*100

# devTime male ----
# look at distribution
ggplot(F1_zscoreDev,aes(x = zscoreDevTimeM, y = treatment))+
  geom_density_ridges(aes(fill = zscoreDiet))

m14 <- lmer(devTimeM ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), F1_zscoreDev, na.action = na.omit)
plot(simulateResiduals(fittedModel = m14)) 
summary(m14)

rept14 <- rpt(devTimeM ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept14a <- rpt(devTimeM ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreDev, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)
print(rept14)
print(rept14a)

#slope variance
sl14 <- (rept14$R$isoline - rept14a$R$isoline) / (rept14$R$isoline + rept14$R$Residual + rept14$R$Fixed)*100

#approaching 0

#isoline variance
iso14 <- (rept14a$R$isoline) / (rept14a$R$isoline + rept14a$R$Residual + rept14a$R$Fixed)*100

#*** fecundity ----
# female
F1_zscoreProg.f <- zscoreF1_f.single

# male
F1_zscoreProg.m <- zscoreF1_m.single

# fecundity female ----
# check distribution
hist(F1_zscoreProg.f$totProg) 
hist(sqrt(F1_zscoreProg.f$totProg))

m15 <- lmer(totProg ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (1|isoline) , F1_zscoreProg.f)
plot(simulateResiduals(fittedModel = m15))
summary(m15)

rept15 <- rpt(totProg ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreProg.f, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept15a <- rpt(totProg ~ zscorePaternal * zscoreFilial + zscoreThoraxF + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreProg.f, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept15)
print(rept15a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl15 <- (rept15$R$isoline - rept15a$R$isoline) / (rept15$R$isoline + rept15$R$Residual + rept15$R$Fixed)*100

#isoline variance
iso15 <- (rept15a$R$isoline) / (rept15a$R$isoline + rept15a$R$Residual + rept15a$R$Fixed)*100

# fecundity male ----
# check distribution
hist(F1_zscoreProg.m$totProg)
hist(sqrt(F1_zscoreProg.m$totProg))

m16 <- lmer(sqrt(totProg) ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline) , F1_zscoreProg.m)
plot(simulateResiduals(fittedModel = m16))
summary(m16)

rept16 <- rpt(totProg ~ zscorePaternal * zscoreFilial + zscoreThoraxM + ( zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreProg.m, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept16a <- rpt(totProg ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreProg.m, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept16)
print(rept16a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl16 <- (rept16$R$isoline - rept16a$R$isoline) / (rept16$R$isoline + rept16$R$Residual + rept16$R$Fixed)*100

#isoline variance
iso16 <- (rept16a$R$isoline) / (rept16a$R$isoline + rept16a$R$Residual + rept16a$R$Fixed)*100

# P2 ----
# import P2
F1_zscoreP2 <- zscoreF1_comp

# remove no matings
F1_zscoreP2 <- subset(F1_zscoreP2,R>0)
# F1_zscoreP2 <- subset(F1_zscoreP2,G>0) # we keep P2=1 in the dataset

# look at distribution
hist(F1_zscoreP2$P2)
hist(asin(sqrt(F1_zscoreP2$P2))) # bit better 

# add OLRE
F1_zscoreP2$obs <- as.factor(1:dim(F1_zscoreP2)[1])

m17 <- glmer(cbind(R,G) ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline) + (1|obs), family=binomial(link = "logit") , data=F1_zscoreP2)
plot(simulateResiduals(fittedModel = m17))
summary(m17)

rept17 <- rptBinary(cbind(R,G) ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (zscorePaternal + zscoreFilial|isoline) + (1|obs), grname = c("isoline", "obs", "Fixed", "Residual"), 
                    data = F1_zscoreP2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept17a <- rptBinary(cbind(R,G) ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline) + (1|obs), grname = c("isoline", "obs", "Fixed", "Residual"), 
                     data = F1_zscoreP2, nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept17)
print(rept17a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl17 <- (rept17$R$isoline[2] - rept17a$R$isoline[2]) / (rept17$R$isoline[2] + rept17$R$Residual[2] + rept17$R$Fixed[2])*100

#approaching 0

#isoline variance
iso17 <- (rept17a$R$isoline[2]) / (rept17a$R$isoline[2] + rept17a$R$Residual[2] + rept17a$R$Fixed[2])*100

# mating duration ----
# import mating duration
F1_zscoreMat <- zscoreF1_m.single

# check distribution
hist(F1_zscoreMat$copDuration)
hist(log(F1_zscoreMat$copDuration))

m18 <- lmer(copDuration ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), F1_zscoreMat, na.action=na.omit)
plot(simulateResiduals(fittedModel = m18)) 
summary(m18)

rept18 <- rpt(copDuration ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreMat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept18a <- rpt(copDuration ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreMat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept18)
print(rept18a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl18 <- (rept18$R$isoline - rept18a$R$isoline) / (rept18$R$isoline + rept18$R$Residual + rept18$R$Fixed)*100

#isoline variance
iso18 <- (rept18a$R$isoline) / (rept18a$R$isoline + rept18a$R$Residual + rept18a$R$Fixed)*100

# mating latency ----
# import mating latency
F1_zscoreLat <- zscoreF1_m.single

# look at distribution
hist(F1_zscoreLat$matingLatency) # no
hist(log(F1_zscoreLat$matingLatency)) # better

# add log transformed data
F1_zscoreLat$loggedmatingLatency <- log(F1_zscoreLat$matingLatency)

# replace "-Inf" values by "0"
F1_zscoreLat[which(F1_zscoreLat$loggedmatingLatency=="-Inf"),"loggedmatingLatency"] <- 0

# final model
m19 <- lmer(loggedmatingLatency ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), F1_zscoreLat) 
plot(simulateResiduals(fittedModel = m19)) 
summary(m19)

rept19 <- rpt(loggedmatingLatency ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (zscorePaternal + zscoreFilial|isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreLat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept19a <- rpt(loggedmatingLatency ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreLat, datatype = "Gaussian", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept19)
print(rept19a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl19 <- (rept19$R$isoline - rept19a$R$isoline) / (rept19$R$isoline + rept19$R$Residual + rept19$R$Fixed)*100

#isoline variance
iso19 <- (rept19a$R$isoline) / (rept19a$R$isoline + rept19a$R$Residual + rept19a$R$Fixed)*100

# remating day ----
F1_zscoreRemat <- zscoreF1_comp

# check distribution
hist(F1_zscoreRemat$remateDay) 

# model
m20 <- glmer(remateDay ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), F1_zscoreRemat, na.action=na.omit, family = "poisson")
plot(simulateResiduals(fittedModel = m20))
summary(m20)

rept20 <- rpt(remateDay ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (zscorePaternal + zscoreFilial |isoline), grname = c("isoline", "Fixed", "Residual"), 
              data = F1_zscoreRemat, datatype = "Poisson", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

rept20a <- rpt(remateDay ~ zscorePaternal * zscoreFilial + zscoreThoraxM + (1|isoline), grname = c("isoline", "Fixed", "Residual"), 
               data = F1_zscoreRemat, datatype = "Poisson", nboot = 0, npermut = 0, parallel = TRUE, ratio = FALSE)

print(rept20)
print(rept20a)

#note that variance components vary slightly each time rpt is run

#slope variance
sl20 <- (rept20$R$isoline[2] - rept20a$R$isoline[2]) / (rept20$R$isoline[2] + rept20$R$Residual[2] + rept20$R$Fixed[2])*100
# approximate 0

#isoline variance
iso20 <- (rept20a$R$isoline[2]) / (rept20a$R$isoline[2] + rept20a$R$Residual[2] + rept20a$R$Fixed[2])*100

# approximate 0

#_______________________----
# FIGURES----
# ***F0*** ----

# set up the jitter
p_dodge <- position_dodge(0.1)

# F0 --- thorax length *** Female  ----
# females
# using this model with untransformed data
F0_thoF_model <- lmer(thoraxF ~ treatment + (1|isoline) , F0_zscoreTho_length)
# the model to extract values from
F0_thoF_model <- emmeans(F0_thoF_model, ~ "treatment")
# extract the values for the plot
F0_thoF_pred<- cld(F0_thoF_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p <- ggplot(F0_thoF_pred, aes(x = treatment, y = emmean, color = treatment)) 
p
pa <- p + geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0.85,1.15), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) 
pa
pb <- pa + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
pb
pc <- pb + xlab("Ancestral condition") + ylab("Thorax length [mm]") + ggtitle("female thorax_length")
pc
pd <- pc + geom_line(data = aggregate(data = F0_zscoreProg.f, thoraxF ~ isoline + treatment, FUN = mean), 
                     aes(x = treatment, y = thoraxF, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
pd

# F0 --- thorax length *** Male  ----
# males
# using this model with untransformed data
F0_thoM_model <- lmer(thoraxM ~ treatment + (1|isoline) , F0_zscoreTho_length)
# the model to extract values from
F0_thoM_model <- emmeans(F0_thoM_model, ~ "treatment")
# extract the values for the plot
F0_thoM_pred<- cld(F0_thoM_model, alpha=0.05, Letters = letters,adjust = "sidak")
# rename the paiwise comparisons
F0_thoM_pred$.group <- c("e","f")

# plot

p0 <- ggplot(F0_thoM_pred, aes(x = treatment, y = emmean, color = treatment)) 
p0
p0a <- p0 + geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0.80,1), labels = scales::number_format(accuracy = 0.01, decimal.mark = '.')) 
p0a
p0b <- p0a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p0b
p0c <- p0b + xlab("Ancestral condition") + ylab("Thorax length [mm]") + ggtitle("male thorax_length")
p0c
p0d <- p0c + geom_line(data = aggregate(data = F0_zscoreProg.m, thoraxM ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = thoraxM, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p0d

multiplot(pd,p0d)

# F0 --- fecundity *** Female  ----
# using this model with untransformed data
F0_fec.f_model <- lmer(totProg ~ treatment + centred_female_TL + (1|isoline) , F0_zscoreProg.f)
# the model to extract values from
F0_fec.f_model <- emmeans(F0_fec.f_model, ~ "treatment")
# extract the values for the plot
F0_fec.f_pred<- cld(F0_fec.f_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_fec.f_pred$link <- "xx"

p7 <- ggplot(F0_fec.f_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p7
p7a <- p7 + 
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0, 75)) 
p7a
p7b <- p7a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p7b
p7c <-  p7b + xlab("Ancestral condition") + ylab("Fecundity [counts]") + ggtitle("fecundity_female") 
p7c
p7d <- p7c + geom_line(data = aggregate(data = F0_zscoreProg.f, totProg ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = totProg, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p7d

# F0 --- fecundity *** Male  ----
# using this model with untransformed data
F0_fec.m_model <- lmer(totProg ~ treatment + centred_male_TL + (1|isoline) , F0_zscoreProg.m)
# the model to extract values from
F0_fec.m_model <- emmeans(F0_fec.m_model, ~ "treatment")
# extract the values for the plot
F0_fec.m_pred<- cld(F0_fec.m_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_fec.m_pred$link <- "xx"

p8 <- ggplot(F0_fec.m_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p8
p8a <- p8 +
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0, 75)) 
p8a
p8b <- p8a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p8b
p8c <-  p8b + xlab("Ancestral condition") + ylab("Fecundity [counts]") + ggtitle("fecundity_male")
p8c
p8d <- p8c + geom_line(data = aggregate(data = F0_zscoreProg.m, totProg ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = totProg, group = isoline, color ="black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p8d

multiplot(p7d,p8d)

# F0 --- development time *** Female  ----
# using this model with untransformed data
F0_devF_model <- lmer(devTimeF ~ treatment + centred_female_TL + (1|isoline) , F0_zscoreDev)
# devTime female
# the model to extract values from
F0_devF_model <- emmeans(F0_devF_model, ~ "treatment")
# extract the values for the plot
F0_devF_pred<- cld(F0_devF_model,alpha=0.05,Letters=letters,adjust="sidak")

# add link column
F0_devF_pred$link <- "xx"

p3 <- ggplot(F0_devF_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p3
p3a <- p3 + 
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(7, 11)) 
p3a
p3b <- p3a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
p3b
p3c <- p3b + xlab("Ancestral condition") + ylab("Development time [days]") + ggtitle("development_female")
p3c
p3d <- p3c + geom_line(data = aggregate(data = F0_zscoreDev, devTimeF ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = devTimeF, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p3d

# F0 --- development time *** Male  ----
# using this model with untransformed data
F0_devM_model <- lmer(devTimeM ~ treatment + centred_male_TL + (1|isoline) , F0_zscoreDev)

# the model to extract values from
F0_devM_model <- emmeans(F0_devM_model, ~ "treatment")
# extract the values for the plot
F0_devM_pred<- cld(F0_devM_model,alpha=0.05,Letters=letters,adjust="sidak")

# add link column
F0_devM_pred$link <- "xx"
p4 <- ggplot(F0_devM_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p4
p4a <- p4 +
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(7, 11))
p4a
p4b <- p4a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
p4b
p4c <-  p4b + xlab("Ancestral condition") + ylab("Development time [days]") + ggtitle("development_male")
p4c
p4d <- p4c + geom_line(data = aggregate(data = F0_zscoreDev, devTimeM ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = devTimeM, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p4d

multiplot(p3d,p4d)

# F0 --- mating latency  ----
# using this model with untransformed data
F0_matLat_model <- lmer(matingLatency ~ treatment + centred_male_TL + (1|isoline) , F0_zscoreLat)
# the model to extract values from
F0_matLat_model <- emmeans(F0_matLat_model, ~ "treatment")
# extract the values for the plot
F0_matLat_pred<- cld(F0_matLat_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_matLat_pred$link <- "xx"

p15 <- ggplot(F0_matLat_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p15
p15a <- p15 + 
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), size=0.5, position = position_dodge(0.3)) +
  scale_y_continuous(limits=c(0,140)) 
p15a
p15b <- p15a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p15b
p15c <-  p15b + xlab("Ancestral condition") + ylab("Mating latency [mins]") + ggtitle("mating_latency")
p15c
p15d <- p15c + geom_line(data = aggregate(data = F0_zscoreLat, matingLatency ~ isoline + treatment, FUN = mean), 
                         aes(x = treatment, y = matingLatency, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p15d

# F0 --- viability  ----
# viability

# using this model with untransformed data
F0_zscoreVia_model <- lmer(viability ~ treatment + (1|isoline), data = F0_zscoreVia2)
# the model to extract values from
F0_zscoreVia_model <- emmeans(F0_zscoreVia_model, ~ "treatment")
# extract the values for the plot
F0_zscoreVia_pred<- cld(F0_zscoreVia_model,alpha=0.05,Letters=letters,adjust="sidak")

# add link column
F0_zscoreVia_pred$link <- "xx"

p1 <- ggplot(F0_zscoreVia_pred, aes(x = treatment, y = emmean,  fill = treatment, color = treatment))
p1
p1a <- p1 + 
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0, 1)) 
p1a
p1b <- p1a + my_theme() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
p1b
p1c <-  p1b + xlab("Ancestral condition") + ylab("viability") + ggtitle("viability")
p1c
p1d <- p1c + geom_line(data = aggregate(data = F0_zscoreVia2, viability ~ isoline + treatment, FUN = mean), 
                       aes(x = treatment, y = viability, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p1d

# F0 --- P2  ----

# using this model with untransformed data
F0_P2_model <- lmer(P2 ~ treatment + centred_male_TL + (1|isoline), F0_zscoreP2)
# the model to extract values from
F0_P2_model <- emmeans(F0_P2_model, ~ "treatment")
# extract the values for the plot
F0_P2_pred<- cld(F0_P2_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_P2_pred$link <- "xx"

p11 <- ggplot(F0_P2_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p11
p11a <- p11 +
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0, 1)) 
p11a
p11b <- p11a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p11b
p11c <-  p11b + xlab("Ancestral condition") + ylab("P2") + ggtitle("P2")
p11c
p11d <- p11c + geom_line(data = aggregate(data = F0_zscoreP2, P2 ~ isoline + treatment, FUN = mean), 
                         aes(x = treatment, y = P2, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p11d

# F0 --- mating duration  ----
# using this model with untransformed data
F0_Mat_model <- lmer(copDuration ~ treatment + centred_male_TL + (1|isoline) , F0_zscoreMat)
# the model to extract values from
F0_Mat_model <- emmeans(F0_Mat_model, ~ "treatment")
# extract the values for the plot
F0_Mat_pred<- cld(F0_Mat_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_Mat_pred$link <- "xx"

p13 <- ggplot(F0_Mat_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p13
p13a <- p13 +
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(10, 40)) 
p13a
p13b <- p13a + my_theme() +  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p13b
p13c <-  p13b + xlab("Ancestral condition") + ylab("mating duration") + ggtitle("mating_duration")
p13c
p13d <- p13c + geom_line(data = aggregate(data = F0_zscoreMat, copDuration ~ isoline + treatment, FUN = mean), 
                         aes(x = treatment, y = copDuration, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p13d

# F0 --- remating day  ----
# using this model with untransformed data
F0_remating_model <- lmer(remateDay ~ treatment + centred_male_TL + (1|isoline) , F0_zscoreRemat)
# the model to extract values from
F0_remating_model <- emmeans(F0_remating_model, ~ "treatment")
# extract the values for the plot
F0_remating_pred<- cld(F0_remating_model,alpha=0.05,Letters=letters,adjust="sidak")
# add link column
F0_remating_pred$link <- "xx"

p17 <- ggplot(F0_remating_pred, aes(x = treatment, y = emmean, fill = treatment, color = treatment))
p17
p17a <- p17 + 
  geom_point(size = 2, position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(3, 5), labels = scales::number_format(accuracy = 1, decimal.mark = '.')) 

p17a
p17b <- p17a + my_theme() + theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) + scale_x_discrete(labels = c("Normal", "Restricted"))
p17b
p17c <-  p17b + xlab("Ancestral condition") + ylab("remating_day") + ggtitle("remating_day")
p17c
p17d <- p17c + geom_line(data = aggregate(data = F0_zscoreRemat, remateDay ~ isoline + treatment, FUN = mean), 
                         aes(x = treatment, y = remateDay, group = isoline, color = "black"), size = 1, alpha = 0.1) +
  scale_colour_manual(values = c("black", "black","black"))
# scale_colour_manual(values = c("black", "#F8766D","#00BFC4"))
p17d

#_______________________----
# FIGURES----
# F1 --- thorax length *** Female  ----
# using this model with untransformed data
F1_thoF_model <- lmer(thoraxF ~ paternal * filial + (1|isoline) , F1_zscoreTho_length)
# the model to extract values from
F1_thoF_model <- emmeans(F1_thoF_model, ~ "paternal:filial")
# extract the values for the plot
F1_thoF_pred<- cld(F1_thoF_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
# p18 <- ggplot(subset(F1_body_size, sex == "F"), aes(x = paternal, y = emmean, group = paste(sex,filial), color = filial)) 
p18 <- ggplot(F1_thoF_pred, aes(x = paternal, y = emmean, group = filial, color = filial)) 
p18
p18a <- p18 + 
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1))
p18a
p18b <- p18a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p18b
p18c <- p18b + xlab("Paternal diet") + ylab("Thorax length [mm]") + ggtitle("female_thorax_length")
p18c 
p18d <- p18c + geom_line(data = aggregate(data = F1_zscoreProg.f, thoraxF ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = thoraxF, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                     labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'))
p18d

# F1 --- thorax length *** Male  ----
# males
# using this model with untransformed data
F1_thoM_model <- lmer(thoraxM ~ paternal * filial + (1|isoline) , F1_zscoreTho_length)
# the model to extract values from
F1_thoM_model <- emmeans(F1_thoM_model, ~ "paternal:filial")
# extract the values for the plot
F1_thoM_pred<- cld(F1_thoM_model, alpha=0.05, Letters = letters,adjust = "sidak")
# rename the paiwise comparisons
F1_thoM_pred$.group <- c("e","f","g","g")

# plot
p19 <-ggplot(F1_thoM_pred, aes(x = paternal, y = emmean, group = filial, color = filial)) 
p19
p19a <- p19 +
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1))
p19a
p19b <- p19a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p19b
p19c <- p19b + xlab("Paternal diet") + ylab("Thorax length [mm]") + ggtitle("male_thorax_length")
p19c
p19d <- p19c + geom_line(data = aggregate(data = F1_zscoreProg.m, thoraxM ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = thoraxM, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4),
                     labels = scales::number_format(accuracy = 0.01, decimal.mark = '.'))

p19d

multiplot(p18d, p19d)

# F1 --- fecundity *** Female  ----
# using this model with untransformed data
F1_progF_model <- lmer(totProg ~ paternal * filial + centred_female_TL + (1|isoline) , F1_zscoreProg.f)
# the model to extract values from
F1_progF_model <- emmeans(F1_progF_model, ~ "paternal:filial")
# extract the values for the plot
F1_progF_pred<- cld(F1_progF_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p26 <- ggplot(F1_progF_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p26
p26a <- p26 + 
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(20,110)) 
p26a
p26b <- p26a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p26b
p26c <- p26b + xlab("Paternal diet") + ylab("Fecundity [counts]") + ggtitle("fecundity_female")
p26c
p26d <- p26c + geom_line(data = aggregate(data = F1_zscoreProg.f, totProg ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = totProg, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p26d

# F1 --- fecundity *** Male  ----
# males 
# using this model with untransformed data
F1_progM_model <- lmer(totProg ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreProg.m)
# the model to extract values from
F1_progM_model <- emmeans(F1_progM_model, ~ "paternal:filial")
# extract the values for the plot
F1_progM_pred<- cld(F1_progM_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p27 <- ggplot(F1_progM_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p27
p27a <- p27 + 
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(20,75)) 
p27a
p27b <- p27a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p27b
p27c <- p27b + xlab("Paternal diet") + ylab("Fecundity [counts]") + ggtitle("fecundity_male")
p27c
p27d <- p27c + geom_line(data = aggregate(data = F1_zscoreProg.m, totProg ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = totProg, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p27d

multiplot(p26d, p27d)

# F1 --- development time *** Female  ----

# using this model with untransformed data
F1_devF_model <- lmer(devTimeF ~ paternal * filial + centred_female_TL + (1|isoline) , F1_zscoreDev)
# the model to extract values from
F1_devF_model <- emmeans(F1_devF_model, ~ "paternal:filial")
# extract the values for the plot
F1_devF_pred<- cld(F1_devF_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p22 <- ggplot(F1_devF_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p22
p22a <- p22 + 
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1)) +
  scale_y_continuous(limits=c(7,10)) 
p22a
p22b <- p22a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p22b
p22c <- p22b + xlab("Paternal diet") + ylab("Development time [days]") + ggtitle("development_female")
p22c
p22d <- p22c + geom_line(data = aggregate(data = F1_zscoreDev, devTimeF ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = devTimeF, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p22d

# F1 --- development time *** Male  ----

F1_devM_model <- lmer(devTimeM ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreDev)
# the model to extract values from
F1_devM_model <- emmeans(F1_devM_model, ~ "paternal:filial")
# extract the values for the plot
F1_devM_pred<- cld(F1_devM_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p23 <- ggplot(F1_devM_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p23
p23a <- p23 + 
  geom_point(size = 2) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1)) +
  scale_y_continuous(limits=c(7,10)) 
p23a
p23b <- p23a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p23b
p23c <- p23b + xlab("Paternal diet") + ylab("Development time [days]") + ggtitle("development_male")
p23c
p23d <- p23c + geom_line(data = aggregate(data = F1_zscoreDev, devTimeM ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = devTimeM, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p23d

multiplot(p22d, p23d)

# F1 --- mating latency  ----
# using this model with untransformed data
F1_matLat_model <- lmer(matingLatency ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreLat)
# the model to extract values from
F1_matLat_model <- emmeans(F1_matLat_model, ~ "paternal:filial")
# extract the values for the plot
F1_matLat_pred<- cld(F1_matLat_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p34 <- ggplot(F1_matLat_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p34
p34a <- p34 + 
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0,110)) 
p34a
p34b <- p34a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p34b
p34c <- p34b + xlab("Paternal diet") + ylab("Mating latency [mins]") + ggtitle("mating_latency")
p34c
p34d <- p34c + geom_line(data = aggregate(data = F1_zscoreLat, matingLatency ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = matingLatency, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p34d

# F1 --- viability  ----

# using this model with untransformed data
F1_via_model <- lmer(viability ~ paternal * filial + (1|isoline) , F1_zscoreVia2)
# the model to extract values from
F1_via_model <- emmeans(F1_via_model, ~ "paternal:filial")
# extract the values for the plot
F1_via_pred<- cld(F1_via_model, alpha=0.05, Letters = letters,adjust = "sidak")

# set up the jitter
p_dodge2 <- position_dodge(0.4)

# plot
p20 <- ggplot(F1_via_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p20
p20a <- p20 + 
  geom_point(size = 2, position = p_dodge2) +
  geom_line(size = 1,  position = p_dodge2) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge2) +
  scale_y_continuous(limits=c(0,1)) 
p20a
p20b <- p20a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p20b
p20c <- p20b + xlab("Paternal diet") + ylab("viability") + ggtitle("viability")
p20c
p20d <- p20c + geom_line(data = aggregate(data = F1_zscoreVia2, viability ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = viability, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p20d

# F1 --- P2  ----
# plot P2
# using this model with untransformed data
F1_zscoreP2_model <- lmer(P2 ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreP2)

# the model to extract values from
F1_zscoreP2_model <- emmeans(F1_zscoreP2_model, ~ "paternal:filial")

# extract the values for the plot
F1_zscoreP2_pred<- cld(F1_zscoreP2_model, alpha=0.05, Letters=letters, adjust="sidak")

p30 <- ggplot(F1_zscoreP2_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p30
p30a <- p30 +
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(0,1)) 
p30a
p30b <- p30a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p30b
p30c <- p30b + xlab("Paternal diet") + ylab("P2") + ggtitle("P2")
p30c
p30d <- p30c + geom_line(data = aggregate(data = F1_zscoreP2, P2 ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = P2, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p30d

# F1 --- mating duration  ----
# using this model with untransformed data
F1_matDur_model <- lmer(copDuration ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreMat)
# the model to extract values from
F1_matDur_model <- emmeans(F1_matDur_model, ~ "paternal:filial")
# extract the values for the plot
F1_matDur_pred<- cld(F1_matDur_model, alpha=0.05, Letters = letters,adjust = "sidak")

# plot
p32 <- ggplot(F1_matDur_pred, aes(x = paternal, y = emmean, group = filial, fill = filial, color = filial)) 
p32
p32a <- p32 + 
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(15,30)) 
p32a
p32b <- p32a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p32b
p32c <- p32b + xlab("Paternal diet") + ylab("mating_duration") + ggtitle("mating_duration")
p32c
p32d <- p32c + geom_line(data = aggregate(data = F1_zscoreMat, copDuration ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = copDuration, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)
p32d

# F1 --- remating day  ----
# using this model with untransformed data
F1_remateDay_model <- lmer(remateDay ~ paternal * filial + centred_male_TL + (1|isoline) , F1_zscoreRemat)
# the model to extract values from
F1_remateDay_model <- emmeans(F1_remateDay_model, ~ "paternal:filial")
# extract the values for the plot
F1_remateDay_pred<- cld(F1_remateDay_model, alpha=0.05, Letters = letters,adjust = "sidak")

p36 <- ggplot(F1_remateDay_pred, aes(x = paternal, y = emmean, fill = filial, group = filial, color = filial))
p36
p36a <- p36 + 
  geom_point(size = 2, position = p_dodge) +
  geom_line(size = 1,  position = p_dodge) +
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL, width=0.1), position = p_dodge) +
  scale_y_continuous(limits=c(2.5,6)) 
p36a
p36b <- p36a + my_theme() + scale_x_discrete(labels = c("Normal", "Restricted"))
p36b
p36c <-  p36b + xlab("Paternal diet") + ylab("remating_day") + ggtitle("remating_day")
p36c
p36d <- p36c + geom_line(data = aggregate(data = F1_zscoreRemat, remateDay ~ isoline + filial + paternal, FUN = mean), 
                         aes(x = paternal, y = remateDay, group = paste(isoline, filial), color = filial), size = 1, alpha = 0.1)

p36d
