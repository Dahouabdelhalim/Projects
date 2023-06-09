###This has all models for "Body size trends in response to climate and urbanization 
#in the widespread North American deer mouse, Peromyscus maniculatus" by Guralnick et al.

library(tidyverse)
library(tidyr)
library(lme4)
library(effects)
library(car)
library(MuMIn)
library(ggplot2)
library(lmerTest) 
library(sjPlot)
library(dplyr)
library(broom)
library(tidyr)

setwd("/Users/Maggie/Dropbox/PEMA_body_size/PEMA_final_body_mass/")

####################################################################################################################################
#Temporal models #Without NEON records 
######## Body mass - with decade and with juvs
pema_dec_BM_no_NEON <- read.csv("pema_new_rezone_bodymass_decade_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_dec_BM_no_NEON$zone <-factor(pema_dec_BM_no_NEON$zone)
#remove zones: 2, 27, 31, 33 - did not have at least 10 redords/decade
pema_dec_BM_no_NEON <- pema_dec_BM_no_NEON[!(pema_dec_BM_no_NEON$zone %in% c(2, 27, 31, 33)),]
pema_dec_BM_no_NEON$season <-factor(pema_dec_BM_no_NEON$season)
pema_dec_BM_no_NEON$sex <-factor(pema_dec_BM_no_NEON$sex)
#scaling and centering
pema_dec_BM_no_NEON2 <- transform(pema_dec_BM_no_NEON, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_dec_BM_no_NEON2 <- pema_dec_BM_no_NEON2 %>% drop_na(zone)

#Models 
pema_bm_dec1 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec10 <- lmer(body_mass ~ MAT + sex + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2,control = lmerControl(optimizer="bobyqa"))
pema_bm_dec13 <- lmer(body_mass ~ MAP + sex + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec15 <- lmer(body_mass ~ MAP + MAT + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec17 <- lmer(body_mass ~ MAP + MAT + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec18 <- lmer(body_mass ~ MAP + sex + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec19 <- lmer(body_mass ~ MAP + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2,control = lmerControl(optimizer="bobyqa"))
pema_bm_dec21 <- lmer(body_mass ~ MAT + sex + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec22 <- lmer(body_mass ~ MAT + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec24 <- lmer(body_mass ~ sex + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec25 <- lmer(body_mass ~ MAP + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec26 <- lmer(body_mass ~ MAT + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec27 <- lmer(body_mass ~ sex + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec28 <- lmer(body_mass ~ season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec29 <- lmer(body_mass ~ pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec30 <- lmer(body_mass ~ MAT + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec31 <- lmer(body_mass ~ MAP + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec32 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec33 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec34 <- lmer(body_mass ~ sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec35 <- lmer(body_mass ~ MAT + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec36 <- lmer(body_mass ~ MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec37 <- lmer(body_mass ~ MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec38 <- lmer(body_mass ~ MAP + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec39 <- lmer(body_mass ~ MAP + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec40 <- lmer(body_mass ~ MAP + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec41 <- lmer(body_mass ~ MAP + MAT + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec42 <- lmer(body_mass ~ MAP + MAT + decade2 + season + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec43 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec44 <- lmer(body_mass ~ decade2 + (1 + decade2 | zone), data=pema_dec_BM_no_NEON2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_dec1, pema_bm_dec2, pema_bm_dec3, pema_bm_dec4, pema_bm_dec5, pema_bm_dec6,
          pema_bm_dec7, pema_bm_dec8, pema_bm_dec9, pema_bm_dec10, pema_bm_dec11, pema_bm_dec12,
          pema_bm_dec13, pema_bm_dec14, pema_bm_dec15, pema_bm_dec16, pema_bm_dec17, pema_bm_dec18,
          pema_bm_dec19, pema_bm_dec20, pema_bm_dec21, pema_bm_dec22, pema_bm_dec23, pema_bm_dec24,
          pema_bm_dec25, pema_bm_dec26, pema_bm_dec27, pema_bm_dec28, pema_bm_dec29, pema_bm_dec30, 
          pema_bm_dec31, pema_bm_dec32, pema_bm_dec33, pema_bm_dec34, pema_bm_dec35, pema_bm_dec36,
          pema_bm_dec37, pema_bm_dec38, pema_bm_dec39, pema_bm_dec40, pema_bm_dec41, pema_bm_dec42, 
          pema_bm_dec43, pema_bm_dec44)

summary(pema_bm_dec33) 
r.squaredGLMM(pema_bm_dec33) 

################################################
######## Body mass - with decade, no juvs
pema_dec_nojuv_BM_noNEON <- read.csv("pema_new_rezone_bodymass_decade_nojuv_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_dec_nojuv_BM_noNEON$zone <-factor(pema_dec_nojuv_BM_noNEON$zone)
#remove zones: 3, 7, 18, 23 - did not have at least 10 redords/decade
pema_dec_nojuv_BM_noNEON <- pema_dec_nojuv_BM_noNEON[!(pema_dec_nojuv_BM_noNEON$zone %in% c(3, 7, 18, 23)),]
pema_dec_nojuv_BM_noNEON$season <-factor(pema_dec_nojuv_BM_noNEON$season)
pema_dec_nojuv_BM_noNEON$sex <-factor(pema_dec_nojuv_BM_noNEON$sex)
#scaling and centering
pema_dec_nojuv_BM_noNEON2 <- transform(pema_dec_nojuv_BM_noNEON, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_dec_nojuv_BM_noNEON2 <- pema_dec_nojuv_BM_noNEON2 %>% drop_na(zone)

##Models 
pema_bm_dec_nojuv1 <- lmer(body_mass ~ MAP + MAT + sex + decade2 +season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv10 <- lmer(body_mass ~ MAT + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv13 <- lmer(body_mass ~ MAP + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv15 <- lmer(body_mass ~ MAP + MAT + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv17 <- lmer(body_mass ~ MAP + MAT + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv18 <- lmer(body_mass ~ MAP + sex + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv19 <- lmer(body_mass ~ MAP + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv21 <- lmer(body_mass ~ MAT + sex + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv22 <- lmer(body_mass ~ MAT + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv24 <- lmer(body_mass ~ sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv25 <- lmer(body_mass ~ MAP + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv26 <- lmer(body_mass ~ MAT + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv27 <- lmer(body_mass ~ sex + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv28 <- lmer(body_mass ~ season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv29 <- lmer(body_mass ~ pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv30 <- lmer(body_mass ~ MAT + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv31 <- lmer(body_mass ~ MAP + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv32 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv33 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv34 <- lmer(body_mass ~ sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv35 <- lmer(body_mass ~ MAT + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv36 <- lmer(body_mass ~ MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv37 <- lmer(body_mass ~ MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv38 <- lmer(body_mass ~ MAP + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv39 <- lmer(body_mass ~ MAP + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv40 <- lmer(body_mass ~ MAP + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv41 <- lmer(body_mass ~ MAP + MAT + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv42 <- lmer(body_mass ~ MAP + MAT + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv43 <- lmer(body_mass ~ MAP + MAT + sex + decade2 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bm_dec_nojuv44 <- lmer(body_mass ~ decade2 + (1 + decade2 | zone), data=pema_dec_nojuv_BM_noNEON2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_dec_nojuv1, pema_bm_dec_nojuv2, pema_bm_dec_nojuv3, pema_bm_dec_nojuv4, pema_bm_dec_nojuv5, pema_bm_dec_nojuv6,
          pema_bm_dec_nojuv7, pema_bm_dec_nojuv8, pema_bm_dec_nojuv9, pema_bm_dec_nojuv10, pema_bm_dec_nojuv11, pema_bm_dec_nojuv12,
          pema_bm_dec_nojuv13, pema_bm_dec_nojuv14, pema_bm_dec_nojuv15, pema_bm_dec_nojuv16, pema_bm_dec_nojuv17, pema_bm_dec_nojuv18,
          pema_bm_dec_nojuv19, pema_bm_dec_nojuv20, pema_bm_dec_nojuv21, pema_bm_dec_nojuv22, pema_bm_dec_nojuv23, pema_bm_dec_nojuv24,
          pema_bm_dec_nojuv25, pema_bm_dec_nojuv26, pema_bm_dec_nojuv27, pema_bm_dec_nojuv28, pema_bm_dec_nojuv29, pema_bm_dec_nojuv30,
          pema_bm_dec_nojuv31, pema_bm_dec_nojuv32, pema_bm_dec_nojuv33, pema_bm_dec_nojuv34, pema_bm_dec_nojuv35, pema_bm_dec_nojuv36,
          pema_bm_dec_nojuv37, pema_bm_dec_nojuv38, pema_bm_dec_nojuv39, pema_bm_dec_nojuv40, pema_bm_dec_nojuv41, pema_bm_dec_nojuv42,
          pema_bm_dec_nojuv43, pema_bm_dec_nojuv44)

summary(pema_bm_dec_nojuv33) 
r.squaredGLMM(pema_bm_dec_nojuv33) 

################################################
######## Head-body length - with decade and with juvs
pema_dec_HBL_noNEON <- read.csv("pema_new_rezone_hblength_decade_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_dec_HBL_noNEON$zone <-factor(pema_dec_HBL_noNEON$zone)
#remove zones: 15, 26 - did not have at least 10 redords/decade
pema_dec_HBL_noNEON <- pema_dec_HBL_noNEON[!(pema_dec_HBL_noNEON$zone %in% c(15, 26)),]
pema_dec_HBL_noNEON$season <-factor(pema_dec_HBL_noNEON$season)
pema_dec_HBL_noNEON$sex <-factor(pema_dec_HBL_noNEON$sex)
#scaling and centering
pema_dec_HBL_noNEON2 <- transform(pema_dec_HBL_noNEON, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_dec_HBL_noNEON2 <- pema_dec_HBL_noNEON2 %>% drop_na(zone)

##Models 
pema_bl_dec1 <- lmer(HB.Length ~ MAP + MAT + sex + season + decade2 + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec10 <- lmer(HB.Length ~ MAT + sex + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bl_dec13 <- lmer(HB.Length ~ MAP + sex + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec15 <- lmer(HB.Length ~ MAP + MAT + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec17 <- lmer(HB.Length ~ MAP + MAT + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec18 <- lmer(HB.Length ~ MAP + sex + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec19 <- lmer(HB.Length ~ MAP + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bl_dec21 <- lmer(HB.Length ~ MAT + sex + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec22 <- lmer(HB.Length ~ MAT + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec24 <- lmer(HB.Length ~ sex + season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec25 <- lmer(HB.Length ~ MAP + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec26 <- lmer(HB.Length ~ MAT + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec27 <- lmer(HB.Length ~ sex + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec28 <- lmer(HB.Length ~ season + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec29 <- lmer(HB.Length ~ pop_10km2_log10 + (1+decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec30 <- lmer(HB.Length ~ MAT + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec31 <- lmer(HB.Length ~ MAP + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec32 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec33 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec34 <- lmer(HB.Length ~ sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec35 <- lmer(HB.Length ~ MAT + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec36 <- lmer(HB.Length ~ MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec37 <- lmer(HB.Length ~ MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec38 <- lmer(HB.Length ~ MAP + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec39 <- lmer(HB.Length ~ MAP + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec40 <- lmer(HB.Length ~ MAP + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec41 <- lmer(HB.Length ~ MAP + MAT + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec42 <- lmer(HB.Length ~ MAP + MAT + decade2 + season + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec43 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec44 <- lmer(HB.Length ~ decade2 + (1 + decade2 | zone), data=pema_dec_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_dec1, pema_bl_dec2, pema_bl_dec3, pema_bl_dec4, pema_bl_dec5, pema_bl_dec6, pema_bl_dec7, 
          pema_bl_dec8, pema_bl_dec9, pema_bl_dec10, pema_bl_dec11, pema_bl_dec12, pema_bl_dec13, pema_bl_dec14, 
          pema_bl_dec15, pema_bl_dec16, pema_bl_dec17, pema_bl_dec18, pema_bl_dec19, pema_bl_dec20, 
          pema_bl_dec21, pema_bl_dec22, pema_bl_dec23, pema_bl_dec24, pema_bl_dec25, pema_bl_dec26, 
          pema_bl_dec27, pema_bl_dec28, pema_bl_dec29, pema_bl_dec30, pema_bl_dec31, pema_bl_dec32,
          pema_bl_dec33, pema_bl_dec34, pema_bl_dec35, pema_bl_dec36, pema_bl_dec37, pema_bl_dec38, 
          pema_bl_dec39, pema_bl_dec40, pema_bl_dec41, pema_bl_dec42, pema_bl_dec43, pema_bl_dec44)

summary(pema_bl_dec33) 
r.squaredGLMM(pema_bl_dec33) 

################################################
######## Head-body length - with decade and no juvs
pema_dec_nojuv_HBL_noNEON <- read.csv("pema_new_rezone_hblength_decade_nojuv_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_dec_nojuv_HBL_noNEON$zone <-factor(pema_dec_nojuv_HBL_noNEON$zone)
#remove zones: 25, 61 - did not have at least 10 redords/decade
pema_dec_nojuv_HBL_noNEON <- pema_dec_nojuv_HBL_noNEON[!(pema_dec_nojuv_HBL_noNEON$zone %in% c(25, 61)),]
pema_dec_nojuv_HBL_noNEON$season <-factor(pema_dec_nojuv_HBL_noNEON$season)
pema_dec_nojuv_HBL_noNEON$sex <-factor(pema_dec_nojuv_HBL_noNEON$sex)
#scaling and centering
pema_dec_nojuv_HBL_noNEON2 <- transform(pema_dec_nojuv_HBL_noNEON, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_dec_nojuv_HBL_noNEON2 <- pema_dec_nojuv_HBL_noNEON2 %>% drop_na(zone)

##Models 
pema_bl_dec_nojuv1 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 +season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv10 <- lmer(HB.Length ~ MAT + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv13 <- lmer(HB.Length ~ MAP + sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv15 <- lmer(HB.Length ~ MAP + MAT + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv17 <- lmer(HB.Length ~ MAP + MAT + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv18 <- lmer(HB.Length ~ MAP + sex + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv19 <- lmer(HB.Length ~ MAP + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2,control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv21 <- lmer(HB.Length ~ MAT + sex + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv22 <- lmer(HB.Length ~ MAT + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv24 <- lmer(HB.Length ~ sex + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv25 <- lmer(HB.Length ~ MAP + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv26 <- lmer(HB.Length ~ MAT + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv27 <- lmer(HB.Length ~ sex + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv28 <- lmer(HB.Length ~ season + (1 | ecoregion3) + (1 | source), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv29 <- lmer(HB.Length ~ pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv30 <- lmer(HB.Length ~ MAT + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv31 <- lmer(HB.Length ~ MAP + sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv32 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv33 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv34 <- lmer(HB.Length ~ sex + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv35 <- lmer(HB.Length ~ MAT + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv36 <- lmer(HB.Length ~ MAT + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv37 <- lmer(HB.Length ~ MAT + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv38 <- lmer(HB.Length ~ MAP + decade2 + season + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv39 <- lmer(HB.Length ~ MAP + sex + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv40 <- lmer(HB.Length ~ MAP + sex + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv41 <- lmer(HB.Length ~ MAP + MAT + decade2 + pop_10km2_log10 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv42 <- lmer(HB.Length ~ MAP + MAT + decade2 + season + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv43 <- lmer(HB.Length ~ MAP + MAT + sex + decade2 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))
pema_bl_dec_nojuv44 <- lmer(HB.Length ~ decade2 + (1 + decade2 | zone), data=pema_dec_nojuv_HBL_noNEON2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_dec_nojuv1, pema_bl_dec_nojuv2, pema_bl_dec_nojuv3, pema_bl_dec_nojuv4, pema_bl_dec_nojuv5, pema_bl_dec_nojuv6,
          pema_bl_dec_nojuv7, pema_bl_dec_nojuv8, pema_bl_dec_nojuv9, pema_bl_dec_nojuv10, pema_bl_dec_nojuv11, pema_bl_dec_nojuv12,
          pema_bl_dec_nojuv13, pema_bl_dec_nojuv14, pema_bl_dec_nojuv15, pema_bl_dec_nojuv16, pema_bl_dec_nojuv17, pema_bl_dec_nojuv18,
          pema_bl_dec_nojuv19, pema_bl_dec_nojuv20, pema_bl_dec_nojuv21, pema_bl_dec_nojuv22, pema_bl_dec_nojuv23, pema_bl_dec_nojuv24,
          pema_bl_dec_nojuv25, pema_bl_dec_nojuv26, pema_bl_dec_nojuv27, pema_bl_dec_nojuv28, pema_bl_dec_nojuv29, pema_bl_dec_nojuv30,
          pema_bl_dec_nojuv31, pema_bl_dec_nojuv32, pema_bl_dec_nojuv33, pema_bl_dec_nojuv34, pema_bl_dec_nojuv35, pema_bl_dec_nojuv36,
          pema_bl_dec_nojuv37, pema_bl_dec_nojuv38, pema_bl_dec_nojuv39, pema_bl_dec_nojuv40, pema_bl_dec_nojuv41, pema_bl_dec_nojuv42,
          pema_bl_dec_nojuv43, pema_bl_dec_nojuv44)

summary(pema_bl_dec_nojuv6) 
r.squaredGLMM(pema_bl_dec_nojuv6) 

####################################################################################################################################
#Spatial models without NEON records 
######## Body mass - without decade and with juvs
pema_nodecade_bodymass_NN <- read.csv("pema_new_rezone_bodymass_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_bodymass_NN$season <-factor(pema_nodecade_bodymass_NN$season)
pema_nodecade_bodymass_NN$sex <-factor(pema_nodecade_bodymass_NN$sex)
#scaling & centering
pema_nodecade_bodymass_NN2 <- transform(pema_nodecade_bodymass_NN, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_bodymass_NN2 <- pema_nodecade_bodymass_NN2 %>% drop_na(ecoregion3)

#models
pema_bm_nn_spat1 <- lmer(body_mass ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat10 <- lmer(body_mass ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat13 <- lmer(body_mass ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat15 <- lmer(body_mass ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat17 <- lmer(body_mass ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat18 <- lmer(body_mass ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat19 <- lmer(body_mass ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat21 <- lmer(body_mass ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat22 <- lmer(body_mass ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat24 <- lmer(body_mass ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat25 <- lmer(body_mass ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat26 <- lmer(body_mass ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat27 <- lmer(body_mass ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat28 <- lmer(body_mass ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_spat29 <- lmer(body_mass ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_nn_spat1, pema_bm_nn_spat2, pema_bm_nn_spat3, pema_bm_nn_spat4, pema_bm_nn_spat5, pema_bm_nn_spat6, pema_bm_nn_spat7,
          pema_bm_nn_spat8, pema_bm_nn_spat9, pema_bm_nn_spat10, pema_bm_nn_spat11, pema_bm_nn_spat12, pema_bm_nn_spat13, pema_bm_nn_spat14,
          pema_bm_nn_spat15, pema_bm_nn_spat16, pema_bm_nn_spat17, pema_bm_nn_spat18, pema_bm_nn_spat19, pema_bm_nn_spat20, pema_bm_nn_spat21,
          pema_bm_nn_spat22, pema_bm_nn_spat23, pema_bm_nn_spat25, pema_bm_nn_spat26, pema_bm_nn_spat27, pema_bm_nn_spat28, pema_bm_nn_spat29)

summary(pema_bm_nn_spat6)
r.squaredGLMM(pema_bm_nn_spat6) 

################################################
######## Body mass - without decade and without juvs
pema_nodecade_nojuv_bodymass_NN <- read.csv("pema_new_rezone_bodymass_nojuv_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_nojuv_bodymass_NN$season <-factor(pema_nodecade_nojuv_bodymass_NN$season)
pema_nodecade_nojuv_bodymass_NN$sex <-factor(pema_nodecade_nojuv_bodymass_NN$sex)
#scaling & centering 
pema_nodecade_nojuv_bodymass_NN2 <- transform(pema_nodecade_nojuv_bodymass_NN, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_nojuv_bodymass_NN2 <- pema_nodecade_nojuv_bodymass_NN2 %>% drop_na(ecoregion3)

#models
pema_bm_nn_nd_nj1 <- lmer(body_mass ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj10 <- lmer(body_mass ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj13 <- lmer(body_mass ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj15 <- lmer(body_mass ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj17 <- lmer(body_mass ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj18 <- lmer(body_mass ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj19 <- lmer(body_mass ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj21 <- lmer(body_mass ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj22 <- lmer(body_mass ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj24 <- lmer(body_mass ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj25 <- lmer(body_mass ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj26 <- lmer(body_mass ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj27 <- lmer(body_mass ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj28 <- lmer(body_mass ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nn_nd_nj29 <- lmer(body_mass ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass_NN2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_nn_nd_nj1,  pema_bm_nn_nd_nj2,  pema_bm_nn_nd_nj3,  pema_bm_nn_nd_nj4,  pema_bm_nn_nd_nj5,  pema_bm_nn_nd_nj6,
          pema_bm_nn_nd_nj7,  pema_bm_nn_nd_nj8,  pema_bm_nn_nd_nj9,  pema_bm_nn_nd_nj10, pema_bm_nn_nd_nj11, pema_bm_nn_nd_nj13,
          pema_bm_nn_nd_nj14, pema_bm_nn_nd_nj15, pema_bm_nn_nd_nj16, pema_bm_nn_nd_nj17, pema_bm_nn_nd_nj18, pema_bm_nn_nd_nj19,
          pema_bm_nn_nd_nj20, pema_bm_nn_nd_nj21, pema_bm_nn_nd_nj22, pema_bm_nn_nd_nj23, pema_bm_nn_nd_nj24, pema_bm_nn_nd_nj25,
          pema_bm_nn_nd_nj26, pema_bm_nn_nd_nj27, pema_bm_nn_nd_nj28, pema_bm_nn_nd_nj29)

summary(pema_bm_nn_nd_nj6)
r.squaredGLMM(pema_bm_nn_nd_nj6) 

################################################
######## HB Length - without decade and with juvs
pema_nodecade_bodylength_NN <- read.csv("pema_new_rezone_hblength_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_bodylength_NN$season <-factor(pema_nodecade_bodylength_NN$season)
pema_nodecade_bodylength_NN$sex <-factor(pema_nodecade_bodylength_NN$sex)
#scaling & centering 
pema_nodecade_bodylength_NN2 <- transform(pema_nodecade_bodylength_NN, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_bodylength_NN2 <- pema_nodecade_bodylength_NN2 %>% drop_na(ecoregion3)

#models
pema_bl_nn_spat1 <- lmer(HB.Length ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat10 <- lmer(HB.Length ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat13 <- lmer(HB.Length ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat15 <- lmer(HB.Length ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat17 <- lmer(HB.Length ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat18 <- lmer(HB.Length ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat19 <- lmer(HB.Length ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat21 <- lmer(HB.Length ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat22 <- lmer(HB.Length ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat24 <- lmer(HB.Length ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat25 <- lmer(HB.Length ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat26 <- lmer(HB.Length ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat27 <- lmer(HB.Length ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat28 <- lmer(HB.Length ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_spat29 <- lmer(HB.Length ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_nn_spat1,  pema_bl_nn_spat2,  pema_bl_nn_spat3,  pema_bl_nn_spat4,  pema_bl_nn_spat5,  pema_bl_nn_spat6,  pema_bl_nn_spat7,
          pema_bl_nn_spat8,  pema_bl_nn_spat9,  pema_bl_nn_spat10, pema_bl_nn_spat11, pema_bl_nn_spat12, pema_bl_nn_spat13, pema_bl_nn_spat14,
          pema_bl_nn_spat15, pema_bl_nn_spat16, pema_bl_nn_spat17, pema_bl_nn_spat18, pema_bl_nn_spat19, pema_bl_nn_spat20, pema_bl_nn_spat21,
          pema_bl_nn_spat22, pema_bl_nn_spat23, pema_bl_nn_spat25, pema_bl_nn_spat26, pema_bl_nn_spat27, pema_bl_nn_spat28, pema_bl_nn_spat29)

summary(pema_bl_nn_spat1)
r.squaredGLMM(pema_bl_nn_spat1) 

################################################
######## HB Length - without decade and without juvs
pema_nodecade_nojuv_bodylength_NN <- read.csv("pema_new_rezone_hblength_nojuv_rezone_centroids_noNEON.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_nojuv_bodylength_NN$season <-factor(pema_nodecade_nojuv_bodylength_NN$season)
pema_nodecade_nojuv_bodylength_NN$sex <-factor(pema_nodecade_nojuv_bodylength_NN$sex)
#scaling & centering
pema_nodecade_nojuv_bodylength_NN2 <- transform(pema_nodecade_nojuv_bodylength_NN, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_nojuv_bodylength_NN2 <- pema_nodecade_nojuv_bodylength_NN2 %>% drop_na(ecoregion3)

#models
pema_bl_nn_nd_nj1 <- lmer(HB.Length ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj10 <- lmer(HB.Length ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj13 <- lmer(HB.Length ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj15 <- lmer(HB.Length ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj17 <- lmer(HB.Length ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj18 <- lmer(HB.Length ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj19 <- lmer(HB.Length ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2,control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj21 <- lmer(HB.Length ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj22 <- lmer(HB.Length ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj24 <- lmer(HB.Length ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj25 <- lmer(HB.Length ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj26 <- lmer(HB.Length ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj27 <- lmer(HB.Length ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj28 <- lmer(HB.Length ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nn_nd_nj29 <- lmer(HB.Length ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength_NN2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_nn_nd_nj1,  pema_bl_nn_nd_nj2,  pema_bl_nn_nd_nj3,  pema_bl_nn_nd_nj4,  pema_bl_nn_nd_nj5,  pema_bl_nn_nd_nj6,
          pema_bl_nn_nd_nj7,  pema_bl_nn_nd_nj8,  pema_bl_nn_nd_nj9,  pema_bl_nn_nd_nj10, pema_bl_nn_nd_nj11, pema_bl_nn_nd_nj13,
          pema_bl_nn_nd_nj14, pema_bl_nn_nd_nj15, pema_bl_nn_nd_nj16, pema_bl_nn_nd_nj17, pema_bl_nn_nd_nj18, pema_bl_nn_nd_nj19,
          pema_bl_nn_nd_nj20, pema_bl_nn_nd_nj21, pema_bl_nn_nd_nj22, pema_bl_nn_nd_nj23, pema_bl_nn_nd_nj24, pema_bl_nn_nd_nj25,
          pema_bl_nn_nd_nj26, pema_bl_nn_nd_nj27, pema_bl_nn_nd_nj28, pema_bl_nn_nd_nj29)

summary(pema_bl_nn_nd_nj1)
r.squaredGLMM(pema_bl_nn_nd_nj1) 

####################################################################################################################################
#Spatial models WITH NEON records 
######## Body mass - without decade and with juvs & with NEON
pema_nodecade_bodymass <- read.csv("pema_new_rezone_bodymass_rezone_centroids.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_bodymass$season <-factor(pema_nodecade_bodymass$season)
pema_nodecade_bodymass$sex <-factor(pema_nodecade_bodymass$sex)
#scaling & centering 
pema_nodecade_bodymass2 <- transform(pema_nodecade_bodymass, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_bodymass2 <- pema_nodecade_bodymass2 %>% drop_na(ecoregion3)

#models
pema_bm_spat1 <- lmer(body_mass ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat10 <- lmer(body_mass ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2,control = lmerControl(optimizer="bobyqa"))
pema_bm_spat13 <- lmer(body_mass ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat15 <- lmer(body_mass ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat17 <- lmer(body_mass ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat18 <- lmer(body_mass ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat19 <- lmer(body_mass ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2,control = lmerControl(optimizer="bobyqa"))
pema_bm_spat21 <- lmer(body_mass ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat22 <- lmer(body_mass ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat24 <- lmer(body_mass ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat25 <- lmer(body_mass ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat26 <- lmer(body_mass ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat27 <- lmer(body_mass ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat28 <- lmer(body_mass ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_spat29 <- lmer(body_mass ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodymass2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_spat1, pema_bm_spat2, pema_bm_spat3, pema_bm_spat4, pema_bm_spat5, pema_bm_spat6, pema_bm_spat7,
          pema_bm_spat8, pema_bm_spat9, pema_bm_spat10, pema_bm_spat11, pema_bm_spat12, pema_bm_spat13, pema_bm_spat14,
          pema_bm_spat15, pema_bm_spat16, pema_bm_spat17, pema_bm_spat18, pema_bm_spat19, pema_bm_spat20, pema_bm_spat21,
          pema_bm_spat22, pema_bm_spat23, pema_bm_spat25, pema_bm_spat26, pema_bm_spat27, pema_bm_spat28, pema_bm_spat29)

summary(pema_bm_spat6)
r.squaredGLMM(pema_bm_spat6) 

################################################
######## Body mass - without decade and without juvs & with NEON
pema_nodecade_nojuv_bodymass <- read.csv("pema_new_rezone_bodymass_nojuv_rezone_centroids.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_nojuv_bodymass$season <-factor(pema_nodecade_nojuv_bodymass$season)
pema_nodecade_nojuv_bodymass$sex <-factor(pema_nodecade_nojuv_bodymass$sex)
#scaling & centering
pema_nodecade_nojuv_bodymass2 <- transform(pema_nodecade_nojuv_bodymass, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_nojuv_bodymass2 <- pema_nodecade_nojuv_bodymass2 %>% drop_na(ecoregion3)

#models
pema_bm_nd_nj1 <- lmer(body_mass ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj2 <- lmer(body_mass ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj3 <- lmer(body_mass ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj4 <- lmer(body_mass ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj5 <- lmer(body_mass ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj6 <- lmer(body_mass ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj7 <- lmer(body_mass ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj8 <- lmer(body_mass ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj9 <- lmer(body_mass ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj10 <- lmer(body_mass ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj11 <- lmer(body_mass ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj12 <- lmer(body_mass ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj13 <- lmer(body_mass ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj14 <- lmer(body_mass ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj15 <- lmer(body_mass ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj16 <- lmer(body_mass ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj17 <- lmer(body_mass ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj18 <- lmer(body_mass ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj19 <- lmer(body_mass ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj20 <- lmer(body_mass ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2,control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj21 <- lmer(body_mass ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj22 <- lmer(body_mass ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj23 <- lmer(body_mass ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj24 <- lmer(body_mass ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj25 <- lmer(body_mass ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj26 <- lmer(body_mass ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj27 <- lmer(body_mass ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj28 <- lmer(body_mass ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))
pema_bm_nd_nj29 <- lmer(body_mass ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodymass2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bm_nd_nj1, pema_bm_nd_nj2, pema_bm_nd_nj3, pema_bm_nd_nj4, pema_bm_nd_nj5, pema_bm_nd_nj6,
          pema_bm_nd_nj7, pema_bm_nd_nj8, pema_bm_nd_nj9, pema_bm_nd_nj10, pema_bm_nd_nj11, pema_bm_nd_nj13,
          pema_bm_nd_nj14, pema_bm_nd_nj15, pema_bm_nd_nj16, pema_bm_nd_nj17, pema_bm_nd_nj18, pema_bm_nd_nj19,
          pema_bm_nd_nj20, pema_bm_nd_nj21, pema_bm_nd_nj22, pema_bm_nd_nj23, pema_bm_nd_nj24, pema_bm_nd_nj25,
          pema_bm_nd_nj26, pema_bm_nd_nj27, pema_bm_nd_nj28, pema_bm_nd_nj29)

summary(pema_bm_nd_nj6)
r.squaredGLMM(pema_bm_nd_nj6) 

################################################
######## HB Length - without decade and with juvs & with NEON
pema_nodecade_bodylength <- read.csv("pema_new_rezone_hblength_rezone_centroids.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_bodylength$season <-factor(pema_nodecade_bodylength$season)
pema_nodecade_bodylength$sex <-factor(pema_nodecade_bodylength$sex)
#scaling & centering 
pema_nodecade_bodylength2 <- transform(pema_nodecade_bodylength, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_bodylength2 <- pema_nodecade_bodylength2 %>% drop_na(ecoregion3)

#models
pema_bl_spat1 <- lmer(HB.Length ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat10 <- lmer(HB.Length ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat13 <- lmer(HB.Length ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat15 <- lmer(HB.Length ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat17 <- lmer(HB.Length ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat18 <- lmer(HB.Length ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat19 <- lmer(HB.Length ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat21 <- lmer(HB.Length ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat22 <- lmer(HB.Length ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat24 <- lmer(HB.Length ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat25 <- lmer(HB.Length ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat26 <- lmer(HB.Length ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat27 <- lmer(HB.Length ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat28 <- lmer(HB.Length ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_spat29 <- lmer(HB.Length ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_bodylength2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_spat1, pema_bl_spat2, pema_bl_spat3, pema_bl_spat4, pema_bl_spat5, pema_bl_spat6, pema_bl_spat7,
          pema_bl_spat8, pema_bl_spat9, pema_bl_spat10, pema_bl_spat11, pema_bl_spat12, pema_bl_spat13, pema_bl_spat14,
          pema_bl_spat15, pema_bl_spat16, pema_bl_spat17, pema_bl_spat18, pema_bl_spat19, pema_bl_spat20, pema_bl_spat21,
          pema_bl_spat22, pema_bl_spat23, pema_bl_spat25, pema_bl_spat26, pema_bl_spat27, pema_bl_spat28, pema_bl_spat29)

summary(pema_bl_spat1)
r.squaredGLMM(pema_bl_spat1) 

################################################
######## HB Length - without decade and without juvs & with NEON
pema_nodecade_nojuv_bodylength <- read.csv("pema_new_rezone_hblength_nojuv_rezone_centroids.csv", header = TRUE, stringsAsFactors = FALSE)
pema_nodecade_nojuv_bodylength$season <-factor(pema_nodecade_nojuv_bodylength$season)
pema_nodecade_nojuv_bodylength$sex <-factor(pema_nodecade_nojuv_bodylength$sex)
#scaling & centering
pema_nodecade_nojuv_bodylength2 <- transform(pema_nodecade_nojuv_bodylength, MAT=scale(MAT), MAP=scale(MAP), pop_10km2_log10=scale(pop_10km2_log10))
#Remove missing data 
pema_nodecade_nojuv_bodylength2 <- pema_nodecade_nojuv_bodylength2 %>% drop_na(ecoregion3)

#models
pema_bl_nd_nj1 <- lmer(HB.Length ~ MAP + MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj2 <- lmer(HB.Length ~ MAT + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj3 <- lmer(HB.Length ~ MAP + sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj4 <- lmer(HB.Length ~ MAP + MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj5 <- lmer(HB.Length ~ MAP + MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj6 <- lmer(HB.Length ~ MAP + MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj7 <- lmer(HB.Length ~ sex + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj8 <- lmer(HB.Length ~ MAT + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj9 <- lmer(HB.Length ~ MAT + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj10 <- lmer(HB.Length ~ MAT + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj11 <- lmer(HB.Length ~ MAP + season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj12 <- lmer(HB.Length ~ MAP + sex + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2,control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj13 <- lmer(HB.Length ~ MAP + sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj14 <- lmer(HB.Length ~ MAP + MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj15 <- lmer(HB.Length ~ MAP + MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj16 <- lmer(HB.Length ~ season + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj17 <- lmer(HB.Length ~ MAP + MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj18 <- lmer(HB.Length ~ MAP + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj19 <- lmer(HB.Length ~ MAP + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj20 <- lmer(HB.Length ~ MAP + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2,control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj21 <- lmer(HB.Length ~ MAT + sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj22 <- lmer(HB.Length ~ MAT + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj23 <- lmer(HB.Length ~ MAT + pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj24 <- lmer(HB.Length ~ sex + season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj25 <- lmer(HB.Length ~ MAP + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj26 <- lmer(HB.Length ~ MAT + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj27 <- lmer(HB.Length ~ sex + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj28 <- lmer(HB.Length ~ season + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))
pema_bl_nd_nj29 <- lmer(HB.Length ~ pop_10km2_log10 + (1 | ecoregion3) + (1 | source), data=pema_nodecade_nojuv_bodylength2, control = lmerControl(optimizer="bobyqa"))

model.sel(pema_bl_nd_nj1, pema_bl_nd_nj2, pema_bl_nd_nj3, pema_bl_nd_nj4, pema_bl_nd_nj5, pema_bl_nd_nj6,
          pema_bl_nd_nj7, pema_bl_nd_nj8, pema_bl_nd_nj9, pema_bl_nd_nj10, pema_bl_nd_nj11, pema_bl_nd_nj13,
          pema_bl_nd_nj14, pema_bl_nd_nj15, pema_bl_nd_nj16, pema_bl_nd_nj17, pema_bl_nd_nj18, pema_bl_nd_nj19,
          pema_bl_nd_nj20, pema_bl_nd_nj21, pema_bl_nd_nj22, pema_bl_nd_nj23, pema_bl_nd_nj24, pema_bl_nd_nj25,
          pema_bl_nd_nj26, pema_bl_nd_nj27, pema_bl_nd_nj28, pema_bl_nd_nj29)

summary(pema_bl_nd_nj1)
r.squaredGLMM(pema_bl_nd_nj1) 

