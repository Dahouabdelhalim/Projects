####################################
###   CONSEQUENCES OF ASR BIAS   ###
####################################

library(lme4)
library(MuMIn)

# Survival  ---------------------------------------------------------------

survival_data <- read.csv("SI_survival_data.csv")

#removing foals

survival_data_no_foals <- subset(survival_data, Age_category != "A")

#and removing those with no subdivision data

survival_data_no_foals <- subset(survival_data_no_foals, !(is.na(Subdivision)))

surv.cons.null <- glmer(Survived ~ Sex + Age_category + scale(Subdivision_density) + (1|Year), data = survival_data_no_foals, family = "binomial")

surv.cons.1 <- glmer(Survived ~ Sex + Age_category + scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_no_foals, family = "binomial")

surv.cons.2 <- glmer(Survived ~ Age_category + Sex*scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_no_foals, family = "binomial")

surv.cons.3 <- glmer(Survived ~ Sex + Age_category*scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_no_foals, family = "binomial")

surv.cons.4 <- glmer(Survived ~ Sex*Age_category*scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_no_foals, family = "binomial")

model.sel(surv.cons.null, surv.cons.1, surv.cons.2, surv.cons.3, surv.cons.4)

summary(surv.cons.null)

# Female reproductive probability -----------------------------------------

repro_data <- read.csv("SI_female_repro_data.csv")

repro_data$ReproStatus <- ifelse(repro_data$ReproStatus == "Reproduced", 1, 0)

repro.cons.null <- glmer(ReproStatus ~ Age_category + ReproStatus_prev_year + scale(Subdivision_density) + (1|Year) + (1|HorseID), data = repro_data, family = "binomial")

repro.cons.1 <- glmer(ReproStatus ~ Age_category + ReproStatus_prev_year + scale(Subdivision_density) + scale(Subdivision_ASR) + (1|Year) + (1|HorseID), data = repro_data, family = "binomial")

repro.cons.2 <- glmer(ReproStatus ~ scale(Subdivision_density) + ReproStatus_prev_year + Age_category*scale(Subdivision_ASR) + (1|Year) + (1|HorseID), data = repro_data, family = "binomial", glmerControl(optimizer = "bobyqa"))

model.sel(repro.cons.null, repro.cons.1, repro.cons.2)

summary(repro.cons.2)

# Foal survival -----------------------------------------------------------

survival_data <- read.csv("SI_survival_data.csv")

#keeping foals only

survival_data_foals <- subset(survival_data, Age_category == "A")

#and removing those with no subdivision data

survival_data_foals <- subset(survival_data_foals, !(is.na(Subdivision)))

foal.surv.null <- glmer(Survived ~ Sex + scale(Subdivision_density) + (1|Year), data = survival_data_foals, family = "binomial")

foal.surv.1 <- glmer(Survived ~ Sex + scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_foals, family = "binomial")

foal.surv.2 <- glmer(Survived ~ Sex*scale(Subdivision_ASR) + scale(Subdivision_density) + (1|Year), data = survival_data_foals, family = "binomial")

model.sel(foal.surv.null, foal.surv.1, foal.surv.2)

summary(foal.surv.null)

# Harem size --------------------------------------------------------------

harem_size <- read.csv("SI_harem_size.csv")

harem_size$Subdivision <- relevel(harem_size$Subdivision, ref="West")

harem.size.1 <- glmer(HaremSize ~ Female_density + (1|Year) + (1|StallionID), family="poisson", data = harem_size, glmerControl(optimizer = "bobyqa"))

harem.size.2 <- glmer(HaremSize ~ Subdivision_ASR + Female_density + (1|Year) + (1|StallionID), family="poisson", data = harem_size, glmerControl(optimizer = "bobyqa"))

model.sel(harem.size.1, harem.size.2)

summary(harem.size.2)

# Band number -------------------------------------------------------------

band_number <- read.csv("SI_band_number.csv")

band_number$Subdivision <- relevel(band_number$Subdivision, ref="West")

band.number.1 <- glmer(Number.of.bands ~ Female_density + (1|Year) + (1|Subdivision), family="poisson", data = band_number)

band.number.2 <- glmer(Number.of.bands ~ Subdivision_ASR + Female_density + (1|Year) + (1|Subdivision), family="poisson", data = band_number)

model.sel(band.number.1, band.number.2)

summary(band.number.2)

# Female band switches ---------------------------------------------------

band_switches <- read.csv("SI_band_switches.csv")

band.switches.m1 <- glmer(Switched_or_not ~ scale(Subdivision_density) + (1|Year) + (1|HorseID), data = band_switches, family = "binomial")

band.switches.m2 <- glmer(Switched_or_not ~ Subdivision_ASR + scale(Subdivision_density) + (1|Year) + (1|HorseID), data = band_switches, family = "binomial")

model.sel(band.switches.m1, band.switches.m2)

summary(band.switches.m1)

