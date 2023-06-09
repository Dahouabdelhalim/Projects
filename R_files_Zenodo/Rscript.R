# R script from:
# Functional biogeography of weeds reveals how anthropogenic management blurs trait-climate relationships

#----------------------------------------------------------#

rm(list=ls())

# LOAD THE PACKAGES
library(openxlsx)
library(nlme)
library(MuMIn)

# LOAD THE DATA
data = read.xlsx(file.choose())
cropland_noHerb = data[which(data$HERBICIDE_TREATMENT == "herbicide_free"),]
cropland_Herb = data[which(data$HERBICIDE_TREATMENT == "herbicide_sprayed"),]
grassland = data[which(data$HABITAT == "grassland"),]

#----------------------------------------------------------#
# COMPARISON BETWEEN GRASSLANDS AND HEBICIDE-FREE CROPLANDS

# CWM_SLA
CWM_SLA_mod_grass = gls(CWM_SLA ~ GSLtw, correlation = corExp(form = ~ Xcoord + Ycoord), data = grassland)
CWM_SLA_mod_weed_noP = lme(CWM_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
summary(CWM_SLA_mod_grass); cor(grassland$CWM_SLA, predict(CWM_SLA_mod_grass))^2 ; intervals(CWM_SLA_mod_grass)
summary(CWM_SLA_mod_weed_noP); r.squaredGLMM(CWM_SLA_mod_weed_noP) ; intervals(CWM_SLA_mod_weed_noP)

# CWM_LDMC
CWM_LDMC_mod_grass = gls(CWM_LDMC ~ GSLtw, correlation = corExp(form = ~ Xcoord + Ycoord), data = grassland)
CWM_LDMC_mod_weed_noP = lme(CWM_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
summary(CWM_LDMC_mod_grass); cor(grassland$CWM_LDMC, predict(CWM_LDMC_mod_grass))^2 ; intervals(CWM_LDMC_mod_grass)
summary(CWM_LDMC_mod_weed_noP); r.squaredGLMM(CWM_LDMC_mod_weed_noP) ; intervals(CWM_LDMC_mod_weed_noP)

# CWM_LNC
CWM_LNC_mod_grass = gls(CWM_LNC ~ GSLtw, correlation = corExp(form = ~ Xcoord + Ycoord), data = grassland)
CWM_LNC_mod_weed_noP = lme(CWM_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
summary(CWM_LNC_mod_grass); cor(grassland$CWM_LNC, predict(CWM_LNC_mod_grass))^2 ; intervals(CWM_LNC_mod_grass)
summary(CWM_LNC_mod_weed_noP); r.squaredGLMM(CWM_LNC_mod_weed_noP) ; intervals(CWM_LNC_mod_weed_noP)

# CWV_SLA
CWV_SLA_mod_grass = gls(CWV_SLA ~ GSLtw, correlation = corRatio(form = ~ Xcoord + Ycoord), data = grassland)
CWV_SLA_mod_weed_noP = lme(CWV_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
summary(CWV_SLA_mod_grass); cor(grassland$CWV_SLA, predict(CWV_SLA_mod_grass))^2 ; intervals(CWV_SLA_mod_grass)
summary(CWV_SLA_mod_weed_noP); r.squaredGLMM(CWV_SLA_mod_weed_noP) ; intervals(CWV_SLA_mod_weed_noP)

# CWV_LDMC
CWV_LDMC_mod_grass = gls(CWV_LDMC ~ GSLtw, correlation = corExp(form = ~ Xcoord + Ycoord), data = grassland)
CWV_LDMC_mod_weed_noP = lme(CWV_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
summary(CWV_LDMC_mod_grass); cor(grassland$CWV_LDMC, predict(CWV_LDMC_mod_grass))^2 ; intervals(CWV_LDMC_mod_grass)
summary(CWV_LDMC_mod_weed_noP); r.squaredGLMM(CWV_LDMC_mod_weed_noP) ; intervals(CWV_LDMC_mod_weed_noP)

# CWV_LNC
CWV_LNC_mod_grass = gls(CWV_LNC ~ GSLtw, correlation = corRatio(form = ~ Xcoord + Ycoord), data = grassland)
CWV_LNC_mod_weed_noP = lme(CWV_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb, na.action = na.exclude)
summary(CWV_LNC_mod_grass); cor(grassland$CWV_LNC, predict(CWV_LNC_mod_grass))^2 ; intervals(CWV_LNC_mod_grass)
summary(CWV_LNC_mod_weed_noP); r.squaredGLMM(CWV_LNC_mod_weed_noP) ; intervals(CWV_LNC_mod_weed_noP)


#-----------------------------------------------------------------#
# COMPARISON BETWEEN HEBICIDE-FREE AND HERBICIDE-SPRAYED CROPLANDS

# CWM SLA
CWM_SLA_mod_weed_noP = lme(CWM_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
CWM_SLA_mod_weed_P = lme(CWM_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWM_SLA_mod_weed_noP); r.squaredGLMM(CWM_SLA_mod_weed_noP) ; intervals(CWM_SLA_mod_weed_noP)
summary(CWM_SLA_mod_weed_P); r.squaredGLMM(CWM_SLA_mod_weed_P) ; intervals(CWM_SLA_mod_weed_P, which = "fixed")

# CWM LDMC
CWM_LDMC_mod_weed_noP = lme(CWM_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
CWM_LDMC_mod_weed_P = lme(CWM_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWM_LDMC_mod_weed_noP); r.squaredGLMM(CWM_LDMC_mod_weed_noP) ; intervals(CWM_LDMC_mod_weed_noP)
summary(CWM_LDMC_mod_weed_P); r.squaredGLMM(CWM_LDMC_mod_weed_P) ; intervals(CWM_LDMC_mod_weed_P)

# CWM LNC
CWM_LNC_mod_weed_noP = lme(CWM_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
CWM_LNC_mod_weed_P = lme(CWM_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corExp(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWM_LNC_mod_weed_noP); r.squaredGLMM(CWM_LNC_mod_weed_noP) ; intervals(CWM_LNC_mod_weed_noP)
summary(CWM_LNC_mod_weed_P); r.squaredGLMM(CWM_LNC_mod_weed_P) ; intervals(CWM_LNC_mod_weed_P)

# CWV SLA
CWV_SLA_mod_weed_noP = lme(CWV_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
CWV_SLA_mod_weed_P = lme(CWV_SLA ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWV_SLA_mod_weed_noP); r.squaredGLMM(CWV_SLA_mod_weed_noP) ; intervals(CWV_SLA_mod_weed_noP)
summary(CWV_SLA_mod_weed_P); r.squaredGLMM(CWV_SLA_mod_weed_P) ; intervals(CWV_SLA_mod_weed_P)

# CWV LDMC
CWV_LDMC_mod_weed_noP = lme(CWV_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb)
CWV_LDMC_mod_weed_P =lme(CWV_LDMC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWV_LDMC_mod_weed_noP); r.squaredGLMM(CWV_LDMC_mod_weed_noP) ; intervals(CWV_LDMC_mod_weed_noP)
summary(CWV_LDMC_mod_weed_P); r.squaredGLMM(CWV_LDMC_mod_weed_P) ; intervals(CWV_LDMC_mod_weed_P)

# CWV LNC
CWV_LNC_mod_weed_noP = lme(CWV_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_noHerb, na.action = na.exclude)
CWV_LNC_mod_weed_P =lme(CWV_LNC ~ GSLtw, random = ~ 1 | CROP_TYPE, correlation = corRatio(form = ~ Xcoord + Ycoord), data = cropland_Herb)
summary(CWV_LNC_mod_weed_noP); r.squaredGLMM(CWV_LNC_mod_weed_noP) ; intervals(CWV_LNC_mod_weed_noP)
summary(CWV_LNC_mod_weed_P); r.squaredGLMM(CWV_LNC_mod_weed_P) ; intervals(CWV_LNC_mod_weed_P)


# END (not run)