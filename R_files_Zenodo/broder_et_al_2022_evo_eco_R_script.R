## packages
library(emmeans)
library(car)
library(logistf)
library(tidyverse)

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## loading in data
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## butterfly trap phonotaxis data from Kalaupapa
butterfly_kalaupapa <- read.csv("butterfly_kalaupapa_data.csv")

## fly trapping data from Kalaupapa
trapping_kalaupapa <- read.csv("trapping_kalaupapa_data.csv")

## butterfly trap phonotaxis data from Wailua summer 2022
butterfly_wailua <- read.csv("butterfly_wailua_data.csv")

## infestation data from Wailua
infestation_wailua <- read.csv("infestation_wailua_data.csv")

## ## fly trapping data from Wailua
trapping_wailua <- read.csv("trapping_wailua_data.csv")

## converting columns to factors
butterfly_kalaupapa$song_type <- as.factor(butterfly_kalaupapa$song_type)
trapping_kalaupapa$song_type <- as.factor(trapping_kalaupapa$song_type)
infestation_wailua$Species <- as.factor(infestation_wailua$Species)
trapping_wailua$song_type <- as.factor(trapping_wailua$song_type)
butterfly_wailua$song_type <- as.factor(butterfly_wailua$song_type)
january_butterfly_both_pops$song_type_2 <- as.factor(january_butterfly_both_pops$song_type)
january_butterfly_both_pops$population <- as.factor(january_butterfly_both_pops$population)
wailua_january$song_type_2 <-as.factor(wailua_january$song_type_2)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Field fly trapping models
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

## glm of Wailua flies caught by songtype, includes all stimuli--------------------------------------------------------------------------------------------------------
glm_traps_wailua <- glm(flies_count ~ song_type, family = poisson, data = trapping_wailua)

summary(glm_traps_wailua)

Anova(glm_traps_wailua, type = "II", test.statistic = "LR")

pairs(emmeans(glm_traps_wailua, ~ song_type))

## removing purring/toceaincus---------------------------------------------------------------------------------------------------------------------------------------------------------
wailua_traps_no_oceanicus <- trapping_wailua %>% filter(!song_type %in% c("purrs"))

## glm of Wailua flies caught by songtype, excluding purring/toceanicus--------------------------------------------------------------------------------------------------------
wailua_traps_no_oceanicus_glm <- glm(flies_count ~ song_type,
                                     family = poisson, data = wailua_traps_no_oceanicus)

summary(wailua_traps_no_oceanicus_glm)

Anova(wailua_traps_no_oceanicus_glm, type = "II", test.statistic = "LR")

pairs(emmeans(wailua_traps_no_oceanicus_glm, ~ song_type))

## glm of Kalaupapa flies caught by songtype, includes all stimuli----------------------------------------------------------------------------------------------------------------------------------------------------------------
glm_traps_kalaupapa <- glm(flies ~ song_type, family = poisson, data = trapping_kalaupapa)

summary(glm_traps_kalaupapa)

Anova(glm_traps_kalaupapa, type = "II", test.statistic = "LR")

pairs(emmeans(glm_traps_kalaupapa, ~ song_type))

## removing purring--------------------------------------------------------------------------------------------------------
kalaupapa_traps_no_purring <- trapping_kalaupapa %>% filter(!song_type %in% c("purring"))

## glm of Kalaupapa flies caught by songtype, excluding purring/toceanicus--------------------------------------------------------------------------------------------------------
glm_traps_kalaupapa_no_purring <- glm(flies ~ song_type, family = poisson, data = kalaupapa_traps_no_purring)

pairs(emmeans(glm_traps_kalaupapa_no_purring, ~ song_type))

summary(glm_traps_kalaupapa_no_purring)

Anova(glm_traps_kalaupapa_no_purring, type = "II", test.statistic = "LR")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Laboratory Phonotaxis models
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

## glm of contact by songtype for Wailua flies
wailua_contact_glm <- glm(contact ~ song_type,
                            family = binomial, na.action = na.omit, data = butterfly_wailua)

summary(wailua_contact_glmer)

Anova(wailua_contact_glmer, type = "II", test.statistic = "LR")

pairs(emmeans(wailua_contact_glmer, ~song_type))

# glm of contact by songtype for Kalaupapa flies
kalaupapa_contact_glm <- glm(contact ~ song_type,
                                family = binomial, na.action = na.omit, data = butterfly_kalaupapa)

summary(kalaupapa_contact_glm)

Anova(kalaupapa_contact_glm, type = "II", test.statistic = "LR")

pairs(emmeans(kalaupapa_contact_glm, ~ song_type))


## logistf Firth's penalized regression of contact by songtype for Wailua flies, all songtypes

## relevel to silence
butterfly_wailua$song_type<-relevel(butterfly_wailua$song_type,"silence")

wailua_contact_logistf <- logistf(formula = contact ~ song_type, data = butterfly_wailua)

summary(wailua_contact_logistf)


## glm contact by songtype for Wailua flies, silence removed
contact_wailua_nosilence_orWN_dat <- butterfly_wailua %>% filter(!song_type %in% c("silence"))

contact_glmer_no_wn <- glm(contact ~ song_type,
                           family = binomial, na.action = na.omit, data = contact_wailua_nosilence_orWN_dat)

pairs(emmeans(contact_glmer_no_wn, ~song_type))



### logistf Firth's penalized regression of contact by songtype for Kalaupapa flies, all songtypes
butterfly_kalaupapa$song_type<-relevel(butterfly_kalaupapa$song_type,"silent")

kalaupapa_contact_logistf <- logistf(contact ~ song_type, data = butterfly_kalaupapa)

summary(kalaupapa_contact_logistf)

## glm contact by songtype for Kalaupapa flies, silence removed
contact_kalaupapa_nosilence <- butterfly_kalaupapa %>% filter(!song_type %in% c("silent"))

kalaupapa_contact_glmer_no_silence <- glm(contact ~ song_type,
                                          family = binomial, na.action = na.omit, data = contact_kalaupapa_nosilence)

Anova(kalaupapa_contact_glmer_no_silence, type = "II", test.statistic = "LR")


pairs(emmeans(kalaupapa_contact_glmer_no_silence, ~ song_type))


## lm of distance by songtype for Wailua flies
lm_distance_wailua <- lm(distance ~ song_type, data = butterfly_wailua)

Anova(lm_distance_wailua, type = "II")

pairs(emmeans(lm_distance_wailua, ~ song_type))


## lm of distance by songtype for Kalaupapa flies
lm_distance_kalaupapa <- lm(distance ~ song_type, data = butterfly_kalaupapa)

Anova(lm_distance_kalaupapa, type = "II")

pairs(emmeans(lm_distance_kalaupapa, ~ song_type))


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Artificial infestation models
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------


## logistf Firth's penalized regression of pupae count by host species

## relevel to t oceanicus
infestation_wailua$Species<-relevel(infestation_wailua$Species,"toceanicus")

infest_logistf <- logistf(pupae_yn ~ Species, data = infestation_wailua, pl = TRUE)

summary(infest_logistf)


# glm of pupae count by host species with t. oceanicus removed
infestation_pupae_no_toceanicus <- infestation_wailua %>% filter(!Species %in% c("toceanicus"))

infestation_pupae_no_toceanicus$Species <- as.factor(infestation_pupae_no_toceanicus$Species)

glm_infestation_pupae_no_oceanicus <- glm(pupae_yn ~ Species,
                                          family = poisson, data = infestation_pupae_no_toceanicus)

Anova(glm_infestation_pupae_no_oceanicus, type = "II", test.statistic = "LR")

pairs(emmeans(glm_infestation_pupae_no_oceanicus, ~ Species))


## logistf Firth's penalized regression of pupae that developed into adults by host species
infestation_adults$Species<-relevel(infestation_adults$Species,"gryllodes")

infest_logistf_adults <- logistf(adult_fly_yn ~ Species, data = infestation_adults, pl = TRUE)

summary(infest_logistf_adults)