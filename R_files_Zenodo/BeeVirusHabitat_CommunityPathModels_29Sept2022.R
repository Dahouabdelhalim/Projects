# Full Code for "Habitat Quality Influences Pollinator Pathogen Prevalence Through 
# Both Habitat-Disease and Biodiversity-Disease Pathways"

# Submitted to: Ecology

# code written by: Michelle L Fearon
# last updated: Sept 29, 2022


# This code includes all analyses, figures, and many of the tables included in the main text, Appendix S1 and Appendix S2
# (each is labeled with it's corresponding figure or table label in the text)


# Note: See the data and code in https://doi.org/10.5061/dryad.zpc866t7g for previous calculation and evaluation of 
# pollinator community species richness, abundance, and estimated species richness at 338 individuals, as well as positive
# and negative strand prevalence in each host species.


# set working directory



#### Abbreviations used in the code below:

## Bee host species
# APME = APIS = Apis = Apis mellifera
# BOIM = BOMB = Bombus = Bombus impatiens
# LASI = Lasi = Lasioglossum spp.
# PEPR = PEPO = Pepo = Eucera pruinosa (previous genus name was Peponapis)

## Viruses
# DWV = deformed wing virus
# BQCV = black queen cell virus
# SBV = sacbrood virus


# Load libraries
library(piecewiseSEM)
library(MuMIn) 
library(lme4)
library(car)
library(DHARMa)
library(ggeffects)
library(dplyr)
library(emmeans)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggiraph)
library(ggiraphExtra)
library(colorBlindness)
library(olsrr)



# DATA FOR SEM:
#Load disease data set (individual bees with binary for each virus)
disease <- read.csv("Virus_Positive_Community+Habitat_Covariates_1Mar2022.csv", stringsAsFactors = F)
head(disease)



#Remove BB site from disease data set because missing local habitat data from one of the site visits 
# (not included in any analyses due the the uneven sampling at that site)
disease <- disease %>%
  filter(Site != "BB")


dim(disease)

unique(disease$Site) # 13 sites


# update data classes to factor class
disease$Year <- as.factor(disease$Year)
disease$Site <- as.factor(disease$Site)
disease$Visit <- as.factor(disease$Visit)
disease$Transect.ID <- as.factor(disease$Transect.ID)
disease$Type <- as.factor(disease$Type)
disease$Species.Code <- as.factor(disease$Species.Code)
disease$Genus.Code <- as.factor(disease$Genus.Code)
disease$Sex <- as.factor(disease$Sex)
disease$Date.Collected <- as.factor(disease$Date.Collected)
summary(disease) #double check the classes
dim(disease)

# transformation functions
# z transformation to scale and center data
z. <- function (x) scale(x)

# log transform all abundance factors
l. <- function (x) log(x)

# arcsine square root transformation function
asqrt <- function (x){
  asin(sqrt(x))
}


# overdispersion test for glmer models
# if overdispersed in a glm --> use quasipoisson and an F test to compare models
# if overdispersed in a glmer --> add ID as a random effect

# Overdispersion parameter estimation function from Ben Bolker: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-for-overdispersioncomputing-overdispersion-factor
# Comment additions from http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

head(disease)

# z transform all variables, and try the log and sqrt transformation on pollinator richness and abundance (see below for tests of best model fit)
disease$Nat.1000_z <- as.numeric(z.(asqrt(disease$Nat.1000))) # arcsine square root transformed
disease$Land_Richness_1000_z <- as.numeric(z.(disease$Land_Richness_1000))
disease$floralrichness_z <- as.numeric(z.(disease$floralrichness))
disease$totaldensity_z <- as.numeric(z.(disease$totaldensity))
disease$Richness_z <- as.numeric(z.(disease$Richness))
disease$Abundance_z <- as.numeric(z.(disease$Abundance))
disease$APME_ABUND_z <- as.numeric(z.(disease$APME_ABUND))
disease$BOIM_ABUND_z <- as.numeric(z.(disease$BOIM_ABUND))
disease$LASI_ABUND_z <- as.numeric(z.(disease$LASI_ABUND))
disease$PEPR_ABUND_z <- as.numeric(z.(disease$PEPR_ABUND))
disease$ApisBomb_ABUND_z <- as.numeric(z.(disease$APME_ABUND + disease$BOIM_ABUND))
disease$EstRich_338_z <- as.numeric(z.(disease$EstRich_338))
disease$Richness_zlog <- as.numeric(z.(l.(disease$Richness)))
disease$Richness_zsqrt <- as.numeric(z.(sqrt(disease$Richness)))
disease$Abundance_zlog <- as.numeric(z.(l.(disease$Abundance)))
disease$Abundance_zsqrt <- as.numeric(z.(sqrt(disease$Abundance)))
disease$APME_ABUND_zlog <- as.numeric(z.(l.(disease$APME_ABUND)))
disease$BOIM_ABUND_zlog <- as.numeric(z.(l.(disease$BOIM_ABUND)))
disease$LASI_ABUND_zlog <- as.numeric(z.(l.(disease$LASI_ABUND)))
disease$PEPR_ABUND_zlog <- as.numeric(z.(l.(disease$PEPR_ABUND+1))) # add one to deal with zero abundance at one site
disease$ApisBomb_ABUND_zlog <- as.numeric(z.(l.(disease$APME_ABUND + disease$BOIM_ABUND)))
disease$APME_ABUND_zsqrt <- as.numeric(z.(sqrt(disease$APME_ABUND)))
disease$BOIM_ABUND_zsqrt <- as.numeric(z.(sqrt(disease$BOIM_ABUND)))
disease$LASI_ABUND_zsqrt <- as.numeric(z.(sqrt(disease$LASI_ABUND)))
disease$PEPR_ABUND_zsqrt <- as.numeric(z.(sqrt(disease$PEPR_ABUND)))
disease$ApisBomb_ABUND_zsqrt <- as.numeric(z.(sqrt(disease$APME_ABUND + disease$BOIM_ABUND)))
disease$EstRich_338_zlog <- as.numeric(z.(l.(disease$EstRich_338)))
disease$EstRich_338_zsqrt <- as.numeric(z.(sqrt(disease$EstRich_338)))
summary(disease)


# Appendix S2: Table S1 means and SDs
# Mean and SD for all site-level continuous variables in the SEMs.
SUMdisease <- disease[1:13, ] %>%
  summarize(richness_mean = mean(Richness), 
            richness_sd = sd(Richness),
            abundance_mean = mean(Abundance), 
            abundance_sd = sd(Abundance),
            Apis_mean = mean(APME_ABUND), 
            Apis_sd = sd(APME_ABUND),
            Bomb_mean = mean(BOIM_ABUND), 
            Bomb_sd = sd(BOIM_ABUND),
            ApisBomb_mean = mean(APME_ABUND+BOIM_ABUND), 
            ApisBomb_sd = sd(APME_ABUND+BOIM_ABUND),
            Lasi_mean = mean(LASI_ABUND), 
            Lasi_sd = sd(LASI_ABUND),
            Pepo_mean = mean(PEPR_ABUND), 
            Pepo_sd = sd(PEPR_ABUND),
            landrich_mean = mean(Land_Richness_1000), 
            landrich_sd = sd(Land_Richness_1000),
            nat_mean = mean(Nat.1000), 
            nat_sd = sd(Nat.1000),
            floralrich_mean = mean(floralrichness), 
            floralrich_sd = sd(floralrichness),
            floralden_mean = mean(totaldensity), 
            floralden_sd = sd(totaldensity))



############
# each data set has a single record from each field site in the first 13 rows 

# subset data by Apis
diseaseAPIS <- disease[ disease$Genus == "APIS", ]
tail(diseaseAPIS)
dim(diseaseAPIS)
summary(diseaseAPIS)


#############
# subset by Bombus
diseaseBOMB <- disease[ disease$Genus == "BOMB", ]
dim(diseaseBOMB)


############
# subset by LASI
diseaseLASI <- disease[ disease$Genus == "LASI", ]
dim(diseaseLASI)


##################
# subset by PEPO
diseasePEPO <- disease[ disease$Genus == "PEPO", ]
dim(diseasePEPO)


#####==============================================================================================
# Full SEM for all hosts and viruses


#### Initial SEM MODEL with Total pollinator abundance:
# each virus modeled in the SEM as an individual with binary output, and host species (Genus.Code), Site, and Visit number to each site as random effects
# the unique site level values (all habitat and pollinator community variables) are in the first 13 rows so those are called in the first two models that are at the site-level below
SEM1 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  Abundance_z %~~% Richness_zsqrt
)
summary(SEM1, conserve = TRUE)


# Appendix S2 Table S7 in Manuscript
# coefficient table for total abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_TotAbundance <- stdCoefs(SEM1, data=disease, 
                                    standardize="range",
                                    standardize.type = "latent.linear",
                                    intercepts = F)
write.csv(std_range_coefs_TotAbundance, "sem_std_range_coefs_TotAbundance_model.csv", row.names=FALSE)



# main model with Apis mellifera abundance replacing total pollinator abundance in the model
SEM1_Apis <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  APME_ABUND_z %~~% Richness_zsqrt
)
summary(SEM1_Apis, conserve = TRUE)

# Appendix S2 Table S8
# coefficient table for Apis mellifera abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_ApisAbundance <- stdCoefs(SEM1_Apis, data=disease, 
                                         standardize="range",
                                         standardize.type = "latent.linear",
                                         intercepts = F)
write.csv(std_range_coefs_ApisAbundance, "sem_std_range_coefs_ApisAbundance.csv", row.names=FALSE)



# main model with Bombus impatiens abundance replacing total pollinator abundance in the model
SEM1_Bomb <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  BOIM_ABUND_z %~~% Richness_zsqrt
)
summary(SEM1_Bomb, conserve = TRUE) # this actually turns out to be interesting!

# Appendix S2 Table S9
# coefficient table for Bombus impatiens abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_BombAbundance <- stdCoefs(SEM1_Bomb, data=disease, 
                                          standardize="range",
                                          standardize.type = "latent.linear",
                                          intercepts = F)
write.csv(std_range_coefs_BombAbundance, "sem_std_range_coefs_BombAbundance.csv", row.names=FALSE)



# main model with Combined Apis mellifera and Bombus impatiens abundance replacing total pollinator abundance in the model
SEM1_ApisBomb <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM1_ApisBomb, conserve = TRUE)


# Appendix S2 Table S3; results shown in Figure 1
# this is the main model presented in the Manuscript text that had the lowest AIC of the models tested.
# coefficient table for combined Apis mellifera and Bombus impatiens abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_ApisBombAbundance <- stdCoefs(SEM1_ApisBomb, data=disease, 
                                          standardize="range",
                                          standardize.type = "latent.linear",
                                          intercepts = F)
write.csv(std_range_coefs_ApisBombAbundance, "sem_std_range_coefs_ApisBombAbundance.csv", row.names=FALSE)

std_coefs_ApisBombAbundance <- stdCoefs(SEM1_ApisBomb, data=disease, 
                                              standardize="scale",
                                              standardize.type = "latent.linear",
                                              intercepts = F)
write.csv(std_coefs_ApisBombAbundance, "sem_std_scale_coefs_ApisBombAbundance.csv", row.names=FALSE)


# main model with Lasioglossum abundance replacing total pollinator abundance in the model
SEM1_Lasi <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(LASI_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  LASI_ABUND_z %~~% Richness_zsqrt
)
summary(SEM1_Lasi, conserve = TRUE)

# Appendix S2 Table S10
# coefficient table for Lasioglossum abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_LasiAbundance <- stdCoefs(SEM1_Lasi, data=disease, 
                                          standardize="range",
                                          standardize.type = "latent.linear",
                                          intercepts = F)
write.csv(std_range_coefs_LasiAbundance, "sem_std_range_coefs_LasiAbundance.csv", row.names=FALSE)



# main model with Eucera pruinosa abundance replacing total pollinator abundance in the model
SEM1_Pepo <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(PEPR_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_zsqrt + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_zsqrt + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_zsqrt + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  PEPR_ABUND_zsqrt %~~% Richness_zsqrt
)
summary(SEM1_Pepo, conserve = TRUE)

# Appendix S2 Table S11
# coefficient table for Eucera pruinosa abundance path model (coefficients, st error, critical value, p-value and standard coefficients)
std_range_coefs_PepoAbundance <- stdCoefs(SEM1_Pepo, data=disease, 
                                          standardize="range",
                                          standardize.type = "latent.linear",
                                          intercepts = F)
write.csv(std_range_coefs_PepoAbundance, "sem_std_range_coefs_PepoAbundance.csv", row.names=FALSE)



# Appendix S2 Table S2
# comparison of AIC for main model with total abundance, and each host species relative abundance
AIC(SEM1, SEM1_Apis, SEM1_Bomb, SEM1_ApisBomb, SEM1_Lasi, SEM1_Pepo)

# best model includes the combined Apis mellifera + Bombus impatiens abundance --> USED THIS MODEL IN THE MAIN TEXT 
# next best model includes Bombus impatiens abundance


## Test of global fit of main model by removing some non-significant links from the model
# main model with Combined Apis mellifera and Bombus impatiens abundance replacing total pollinator abundance in the model
SEM1_ApisBomb_test <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~  Nat.1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM1_ApisBomb_test, conserve = TRUE)
std_range_coefs_ApisBombAbundance_test <- stdCoefs(SEM1_ApisBomb_test, data=disease, 
                                              standardize="range",
                                              standardize.type = "latent.linear",
                                              intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit



#####==============================================================================================
#### Check for component model assumptions

# comparison of Richness transformations
mod1 <- lm(Richness_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1)
vif(mod1) # all VIFs are less than 2
par(mfrow = c(2, 2))
plot(mod1)


mod1a <- lm(Richness_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1a)
vif(mod1a) # all VIFs are less than 2
plot(mod1a)


mod1b <- lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1b)
vif(mod1b)  # all VIFs are less than 2
plot(mod1b)
#best fit


#normality
ols_plot_resid_qq(mod1)
ols_plot_resid_qq(mod1a)
ols_plot_resid_qq(mod1b)

ols_test_normality(mod1)
ols_test_normality(mod1a)
ols_test_normality(mod1b)

#Homogeneity of variance
ols_test_breusch_pagan(mod1)
ols_test_breusch_pagan(mod1a)
ols_test_breusch_pagan(mod1b)




# comparison of pollinator abundance transformations
mod2 <- lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod2)
vif(mod2) # all VIFs are less than 2
plot(mod2) # best fit

mod2a <- lm(Abundance_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod2a)
vif(mod2a) # all VIFs are less than 2
plot(mod2a) # this doesn't look good at all - don't use this transformation

mod2b <- lm(Abundance_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod2b)
vif(mod2b) # all VIFs are less than 2
plot(mod2b)


#normality
ols_plot_resid_qq(mod2)
ols_plot_resid_qq(mod2a)
ols_plot_resid_qq(mod2b)

ols_test_normality(mod2)
ols_test_normality(mod2a)
ols_test_normality(mod2b)

#Homogeneity of variance
ols_test_breusch_pagan(mod2)
ols_test_breusch_pagan(mod2a)
ols_test_breusch_pagan(mod2b)



# comparison of Apis mellifera abundance transformations
mod3 <- lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod3)
vif(mod3) # all VIFs are less than 2
plot(mod3) # best fit

mod3a <- lm(APME_ABUND_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod3a)
vif(mod3a) # all VIFs are less than 2
plot(mod3a)

mod3b <- lm(APME_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod3b)
vif(mod3b) # all VIFs are less than 2
plot(mod3b)


#normality
ols_plot_resid_qq(mod3)
ols_plot_resid_qq(mod3a)
ols_plot_resid_qq(mod3b)

ols_test_normality(mod3)
ols_test_normality(mod3a)
ols_test_normality(mod3b)

#Homogeneity of variance
ols_test_breusch_pagan(mod3)
ols_test_breusch_pagan(mod3a)
ols_test_breusch_pagan(mod3b)



#comparison of Bombus impatiens abundance transformations
mod4 <- lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod4)
vif(mod4) # all VIFs are less than 2
plot(mod4) # best fit

mod4a <- lm(BOIM_ABUND_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod4a)
vif(mod4a) # all VIFs are less than 2
plot(mod4a)

mod4b <- lm(BOIM_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod4b)
vif(mod4b) # all VIFs are less than 2
plot(mod4b)


#normality
ols_plot_resid_qq(mod4)
ols_plot_resid_qq(mod4a)
ols_plot_resid_qq(mod4b)

ols_test_normality(mod4)
ols_test_normality(mod4a)
ols_test_normality(mod4b)

#Homogeneity of variance
ols_test_breusch_pagan(mod4)
ols_test_breusch_pagan(mod4a)
ols_test_breusch_pagan(mod4b)


#comparison of Combined Apis mellifera and Bombus impatiens abundance transformations
mod5 <- lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod5)
vif(mod5) # all VIFs are less than 2
plot(mod5) # best fit

mod5a <- lm(ApisBomb_ABUND_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod5a)
vif(mod5a) # all VIFs are less than 2
plot(mod5a)

mod5b <- lm(ApisBomb_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod5b)
vif(mod5b) # all VIFs are less than 2
plot(mod5b)


#normality
ols_plot_resid_qq(mod5)
ols_plot_resid_qq(mod5a)
ols_plot_resid_qq(mod5b)

ols_test_normality(mod5)
ols_test_normality(mod5a)
ols_test_normality(mod5b)

#Homogeneity of variance
ols_test_breusch_pagan(mod5)
ols_test_breusch_pagan(mod5a)
ols_test_breusch_pagan(mod5b)


#comparison of Lasioglossum abundance transformations
mod6 <- lm(LASI_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod6)
vif(mod6) # all VIFs are less than 2
plot(mod6) # best fit

mod6a <- lm(LASI_ABUND_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod6a)
vif(mod6a) # all VIFs are less than 2
plot(mod6a)

mod6b <- lm(LASI_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod6b)
vif(mod6b) # all VIFs are less than 2
plot(mod6b)


#normality
ols_plot_resid_qq(mod6)
ols_plot_resid_qq(mod6a)
ols_plot_resid_qq(mod6b)

ols_test_normality(mod6)
ols_test_normality(mod6a)
ols_test_normality(mod6b)

#Homogeneity of variance
ols_test_breusch_pagan(mod6)
ols_test_breusch_pagan(mod6a)
ols_test_breusch_pagan(mod6b)


#comparison of Eucera pruinosa abundance transformations
mod7 <- lm(PEPR_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod7)
vif(mod7) # all VIFs are less than 2
plot(mod7) 

mod7a <- lm(PEPR_ABUND_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod7a)
vif(mod7a) # all VIFs are less than 2
plot(mod7a)

mod7b <- lm(PEPR_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod7b)
vif(mod7b) # all VIFs are less than 2
plot(mod7b) # best fit


#normality
ols_plot_resid_qq(mod7)
ols_plot_resid_qq(mod7a)
ols_plot_resid_qq(mod7b)

ols_test_normality(mod7)
ols_test_normality(mod7a)
ols_test_normality(mod7b)

#Homogeneity of variance
ols_test_breusch_pagan(mod7)
ols_test_breusch_pagan(mod7a)
ols_test_breusch_pagan(mod7b)





# check model assumptions for DWV, BQCV and SBV
mod8 <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8)
vif(mod8) # VIFs are all <= 3.6
overdisp_fun(mod8)
plot(mod8)

mod8a <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8a)
vif(mod8a) # VIFs are all < 4.5
overdisp_fun(mod8a)
plot(mod8a)

mod8b <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8b)
vif(mod8b) # VIFs are all < 4.2
overdisp_fun(mod8b)
plot(mod8b)

mod8c <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8c)
vif(mod8c) # VIFs are all < 3.5
overdisp_fun(mod8c)
plot(mod8c)

mod8d <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8d)
vif(mod8d) # VIFs are all < 4.6
overdisp_fun(mod8d)
plot(mod8d)

mod8e <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod8e)
vif(mod8e) # VIFs are all < 3.4
overdisp_fun(mod8e)
plot(mod8e)



# check model assumptions for BQCV
mod9 <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9)
vif(mod9) # VIFs are all < 3.4
overdisp_fun(mod9)
plot(mod9)

mod9a <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9a)
vif(mod9a) # VIFs are all < 4.2
overdisp_fun(mod9a)
plot(mod9a)

mod9b <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9b)
vif(mod9b) # VIFs are all < 4.4
overdisp_fun(mod9b)
plot(mod9b)

mod9c <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9c)
vif(mod9c) # VIFs are all < 3.4
overdisp_fun(mod9c)
plot(mod9c)

mod9d <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9d)
vif(mod9d) # VIFs are all < 3.7
overdisp_fun(mod9d)
plot(mod9d)

mod9e <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod9e)
vif(mod9e) # VIFs are all < 3.6
overdisp_fun(mod9e)
plot(mod9e)


# check model assumptions for SBV
mod10 <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod10)
vif(mod10) # VIFs are all < 3
overdisp_fun(mod10)
plot(mod10)

mod10a <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod10a)
vif(mod10a) # VIFs are all < 3.9
overdisp_fun(mod10a)
plot(mod10a)

mod10b <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
summary(mod10b)
vif(mod10b) # VIFs are all < 3.9
overdisp_fun(mod10b)
plot(mod10b)

mod10c <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod10c)
vif(mod10c) # VIFs are all < 3.3
overdisp_fun(mod10c)
plot(mod10c)

mod10d <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod10d)
vif(mod10d) # VIFs are all < 3.8
overdisp_fun(mod10d)
plot(mod10d)

mod10e <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease)
summary(mod10e)
vif(mod10e) # VIFs are all < 3.4
overdisp_fun(mod10e)
plot(mod10e)


#####==============================================================================================
# APPENDIX S2 Table S5
# correlation table among factors in the Main SEMs including each species-specific pollinator abundance term
subset <- disease %>%
  select(DWV, BQCV, SBV, Richness_zsqrt, Abundance_z, APME_ABUND_z, BOIM_ABUND_z, ApisBomb_ABUND_z, LASI_ABUND_z, PEPR_ABUND_zsqrt, Land_Richness_1000_z, Nat.1000_z, floralrichness_z, totaldensity_z)

correlations <- cor(subset)

# add p-values for each correlation above the diagonal
correlations["DWV", "BQCV"] <- unlist(cor.test(subset$DWV, subset$BQCV)[3])
correlations["DWV", "SBV"] <- unlist(cor.test(subset$DWV, subset$SBV)[3])
correlations["DWV", "Richness_zsqrt"] <- unlist(cor.test(subset$DWV, subset$Richness_zsqrt)[3])
correlations["DWV", "Abundance_z"] <- unlist(cor.test(subset$DWV, subset$Abundance_z)[3])
correlations["DWV", "APME_ABUND_z"] <- unlist(cor.test(subset$DWV, subset$APME_ABUND_z)[3])
correlations["DWV", "BOIM_ABUND_z"] <- unlist(cor.test(subset$DWV, subset$BOIM_ABUND_z)[3])
correlations["DWV", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$DWV, subset$ApisBomb_ABUND_z)[3])
correlations["DWV", "LASI_ABUND_z"] <- unlist(cor.test(subset$DWV, subset$LASI_ABUND_z)[3])
correlations["DWV", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$DWV, subset$PEPR_ABUND_zsqrt)[3])
correlations["DWV", "Land_Richness_1000_z"] <- unlist(cor.test(subset$DWV, subset$Land_Richness_1000_z)[3])
correlations["DWV", "Nat.1000_z"] <- unlist(cor.test(subset$DWV, subset$Nat.1000_z)[3])
correlations["DWV", "floralrichness_z"] <- unlist(cor.test(subset$DWV, subset$floralrichness_z)[3])
correlations["DWV", "totaldensity_z"] <- unlist(cor.test(subset$DWV, subset$totaldensity_z)[3])
correlations["BQCV", "SBV"] <- unlist(cor.test(subset$BQCV, subset$SBV)[3])
correlations["BQCV", "Richness_zsqrt"] <- unlist(cor.test(subset$BQCV, subset$Richness_zsqrt)[3])
correlations["BQCV", "Abundance_z"] <- unlist(cor.test(subset$BQCV, subset$Abundance_z)[3])
correlations["BQCV", "APME_ABUND_z"] <- unlist(cor.test(subset$BQCV, subset$APME_ABUND_z)[3])
correlations["BQCV", "BOIM_ABUND_z"] <- unlist(cor.test(subset$BQCV, subset$BOIM_ABUND_z)[3])
correlations["BQCV", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$BQCV, subset$ApisBomb_ABUND_z)[3])
correlations["BQCV", "LASI_ABUND_z"] <- unlist(cor.test(subset$BQCV, subset$LASI_ABUND_z)[3])
correlations["BQCV", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$BQCV, subset$PEPR_ABUND_zsqrt)[3])
correlations["BQCV", "Land_Richness_1000_z"] <- unlist(cor.test(subset$BQCV, subset$Land_Richness_1000_z)[3])
correlations["BQCV", "Nat.1000_z"] <- unlist(cor.test(subset$BQCV, subset$Nat.1000_z)[3])
correlations["BQCV", "floralrichness_z"] <- unlist(cor.test(subset$BQCV, subset$floralrichness_z)[3])
correlations["BQCV", "totaldensity_z"] <- unlist(cor.test(subset$BQCV, subset$totaldensity_z)[3])
correlations["SBV", "Richness_zsqrt"] <- unlist(cor.test(subset$SBV, subset$Richness_zsqrt)[3])
correlations["SBV", "Abundance_z"] <- unlist(cor.test(subset$SBV, subset$Abundance_z)[3])
correlations["SBV", "APME_ABUND_z"] <- unlist(cor.test(subset$SBV, subset$APME_ABUND_z)[3])
correlations["SBV", "BOIM_ABUND_z"] <- unlist(cor.test(subset$SBV, subset$BOIM_ABUND_z)[3])
correlations["SBV", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$SBV, subset$ApisBomb_ABUND_z)[3])
correlations["SBV", "LASI_ABUND_z"] <- unlist(cor.test(subset$SBV, subset$LASI_ABUND_z)[3])
correlations["SBV", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$SBV, subset$PEPR_ABUND_zsqrt)[3])
correlations["SBV", "Land_Richness_1000_z"] <- unlist(cor.test(subset$SBV, subset$Land_Richness_1000_z)[3])
correlations["SBV", "Nat.1000_z"] <- unlist(cor.test(subset$SBV, subset$Nat.1000_z)[3])
correlations["SBV", "floralrichness_z"] <- unlist(cor.test(subset$SBV, subset$floralrichness_z)[3])
correlations["SBV", "totaldensity_z"] <- unlist(cor.test(subset$SBV, subset$totaldensity_z)[3])
correlations["Richness_zsqrt", "Abundance_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$Abundance_)[3])
correlations["Richness_zsqrt", "APME_ABUND_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$APME_ABUND_z)[3])
correlations["Richness_zsqrt", "BOIM_ABUND_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$BOIM_ABUND_z)[3])
correlations["Richness_zsqrt", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$ApisBomb_ABUND_z)[3])
correlations["Richness_zsqrt", "LASI_ABUND_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$LASI_ABUND_z)[3])
correlations["Richness_zsqrt", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$Richness_zsqrt, subset$PEPR_ABUND_zsqrt)[3])
correlations["Richness_zsqrt", "Land_Richness_1000_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$Land_Richness_1000_z)[3])
correlations["Richness_zsqrt", "Nat.1000_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$Nat.1000_z)[3])
correlations["Richness_zsqrt", "floralrichness_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$floralrichness_z)[3])
correlations["Richness_zsqrt", "totaldensity_z"] <- unlist(cor.test(subset$Richness_zsqrt, subset$totaldensity_z)[3])
correlations["Abundance_z", "Land_Richness_1000_z"] <- unlist(cor.test(subset$Abundance_z, subset$Land_Richness_1000_z)[3])
correlations["Abundance_z", "Nat.1000_z"] <- unlist(cor.test(subset$Abundance_z, subset$Nat.1000_z)[3])
correlations["Abundance_z", "floralrichness_z"] <- unlist(cor.test(subset$Abundance_z, subset$floralrichness_z)[3])
correlations["Abundance_z", "totaldensity_z"] <- unlist(cor.test(subset$Abundance_z, subset$totaldensity_z)[3])
correlations["Abundance_z", "APME_ABUND_z"] <- unlist(cor.test(subset$Abundance_z, subset$APME_ABUND_z)[3])
correlations["Abundance_z", "BOIM_ABUND_z"] <- unlist(cor.test(subset$Abundance_z, subset$BOIM_ABUND_z)[3])
correlations["Abundance_z", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$Abundance_z, subset$ApisBomb_ABUND_z)[3])
correlations["Abundance_z", "LASI_ABUND_z"] <- unlist(cor.test(subset$Abundance_z, subset$LASI_ABUND_z)[3])
correlations["Abundance_z", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$Abundance_z, subset$PEPR_ABUND_zsqrt)[3])
correlations["APME_ABUND_z", "BOIM_ABUND_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$BOIM_ABUND_z)[3])
correlations["APME_ABUND_z", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$ApisBomb_ABUND_z)[3])
correlations["APME_ABUND_z", "LASI_ABUND_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$LASI_ABUND_z)[3])
correlations["APME_ABUND_z", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$APME_ABUND_z, subset$PEPR_ABUND_zsqrt)[3])
correlations["APME_ABUND_z", "Land_Richness_1000_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$Land_Richness_1000_z)[3])
correlations["APME_ABUND_z", "Nat.1000_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$Nat.1000_z)[3])
correlations["APME_ABUND_z", "floralrichness_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$floralrichness_z)[3])
correlations["APME_ABUND_z", "totaldensity_z"] <- unlist(cor.test(subset$APME_ABUND_z, subset$totaldensity_z)[3])
correlations["BOIM_ABUND_z", "ApisBomb_ABUND_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$ApisBomb_ABUND_z)[3])
correlations["BOIM_ABUND_z", "LASI_ABUND_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$LASI_ABUND_z)[3])
correlations["BOIM_ABUND_z", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$PEPR_ABUND_zsqrt)[3])
correlations["BOIM_ABUND_z", "Land_Richness_1000_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$Land_Richness_1000_z)[3])
correlations["BOIM_ABUND_z", "Nat.1000_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$Nat.1000_z)[3])
correlations["BOIM_ABUND_z", "floralrichness_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$floralrichness_z)[3])
correlations["BOIM_ABUND_z", "totaldensity_z"] <- unlist(cor.test(subset$BOIM_ABUND_z, subset$totaldensity_z)[3])
correlations["ApisBomb_ABUND_z", "LASI_ABUND_z"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$LASI_ABUND_z)[3])
correlations["ApisBomb_ABUND_z", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$PEPR_ABUND_zsqrt)[3])
correlations["ApisBomb_ABUND_z", "Land_Richness_1000_z"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$Land_Richness_1000_z)[3])
correlations["ApisBomb_ABUND_z", "Nat.1000_z"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$Nat.1000_z)[3])
correlations["ApisBomb_ABUND_z", "floralrichness_z"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$floralrichness_z)[3])
correlations["ApisBomb_ABUND_z", "totaldensity_z"] <- unlist(cor.test(subset$ApisBomb_ABUND_z, subset$totaldensity_z)[3])
correlations["LASI_ABUND_z", "PEPR_ABUND_zsqrt"] <- unlist(cor.test(subset$LASI_ABUND_z, subset$PEPR_ABUND_zsqrt)[3])
correlations["LASI_ABUND_z", "Land_Richness_1000_z"] <- unlist(cor.test(subset$LASI_ABUND_z, subset$Land_Richness_1000_z)[3])
correlations["LASI_ABUND_z", "Nat.1000_z"] <- unlist(cor.test(subset$LASI_ABUND_z, subset$Nat.1000_z)[3])
correlations["LASI_ABUND_z", "floralrichness_z"] <- unlist(cor.test(subset$LASI_ABUND_z, subset$floralrichness_z)[3])
correlations["LASI_ABUND_z", "totaldensity_z"] <- unlist(cor.test(subset$LASI_ABUND_z, subset$totaldensity_z)[3])
correlations["PEPR_ABUND_zsqrt", "Land_Richness_1000_z"] <- unlist(cor.test(subset$PEPR_ABUND_zsqrt, subset$Land_Richness_1000_z)[3])
correlations["PEPR_ABUND_zsqrt", "Nat.1000_z"] <- unlist(cor.test(subset$PEPR_ABUND_zsqrt, subset$Nat.1000_z)[3])
correlations["PEPR_ABUND_zsqrt", "floralrichness_z"] <- unlist(cor.test(subset$PEPR_ABUND_zsqrt, subset$floralrichness_z)[3])
correlations["PEPR_ABUND_zsqrt", "totaldensity_z"] <- unlist(cor.test(subset$PEPR_ABUND_zsqrt, subset$totaldensity_z)[3])
correlations["Land_Richness_1000_z", "Nat.1000_z"] <- unlist(cor.test(subset$Land_Richness_1000_z, subset$Nat.1000_z)[3])
correlations["Land_Richness_1000_z", "floralrichness_z"] <- unlist(cor.test(subset$Land_Richness_1000_z, subset$floralrichness_z)[3])
correlations["Land_Richness_1000_z", "totaldensity_z"] <- unlist(cor.test(subset$Land_Richness_1000_z, subset$totaldensity_z)[3])
correlations["Nat.1000_z", "floralrichness_z"] <- unlist(cor.test(subset$Nat.1000_z, subset$floralrichness_z)[3])
correlations["Nat.1000_z", "totaldensity_z"] <- unlist(cor.test(subset$Nat.1000_z, subset$totaldensity_z)[3])
correlations["floralrichness_z", "totaldensity_z"] <- unlist(cor.test(subset$floralrichness_z, subset$totaldensity_z)[3])


correlations <- format(round(correlations, 3))
# set diagonal to NA
diag(correlations) <- "--"
rownames(correlations) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(correlations) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(correlations, "SEM_factor_correlations.csv", quote = F, row.names = T)


# covariance matrix for the main SEM (not shown in manuscript)
covariances <- format(round(cov(subset), 3))
rownames(covariances) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(covariances) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(covariances, "SEM_factor_covariances.csv", quote = F, row.names = T)




#####==============================================================================================
# APPENDIX S2 Table S6
# Tests of spatial autocorrelation for each component model in the SEM
mod1_resid <- simulateResiduals(mod1b)
mod1_resid2 <- recalculateResiduals(mod1_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod1_resid2, unique(disease$Long), unique(disease$Lat))

mod2_resid <- simulateResiduals(mod2)
mod2_resid2 <- recalculateResiduals(mod2_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod2_resid2, unique(disease$Long), unique(disease$Lat))

mod3_resid <- simulateResiduals(mod3)
mod3_resid2 <- recalculateResiduals(mod3_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod3_resid2, unique(disease$Long), unique(disease$Lat))

mod4_resid <- simulateResiduals(mod4)
mod4_resid2 <- recalculateResiduals(mod4_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod4_resid2, unique(disease$Long), unique(disease$Lat))

mod5_resid <- simulateResiduals(mod5)
mod5_resid2 <- recalculateResiduals(mod5_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod5_resid2, unique(disease$Long), unique(disease$Lat))

mod6_resid <- simulateResiduals(mod6)
mod6_resid2 <- recalculateResiduals(mod6_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod6_resid2, unique(disease$Long), unique(disease$Lat))

mod7b_resid <- simulateResiduals(mod7b)
mod7b_resid2 <- recalculateResiduals(mod7b_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod7b_resid2, unique(disease$Long), unique(disease$Lat))

mod8_resid <- simulateResiduals(mod8c)
mod8_resid2 <- recalculateResiduals(mod8_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod8_resid2, unique(disease$Long), unique(disease$Lat))

mod9_resid <- simulateResiduals(mod9c)
mod9_resid2 <- recalculateResiduals(mod9_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod9_resid2, unique(disease$Long), unique(disease$Lat))

mod10_resid <- simulateResiduals(mod10c)
mod10_resid2 <- recalculateResiduals(mod10_resid, group = disease$Site) # group residuals by site
testSpatialAutocorrelation(mod10_resid2, unique(disease$Long), unique(disease$Lat))



#####==============================================================================================
# Figure 2 in Manuscript
# Figure of virus prev (all three viruses) vs habitat characteristics and pollinator community characteristics


# function to make a grid plot with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}


# color blind palette (from https://zenodo.org/record/3381072#.X8VKS2hKjZt)
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
# Goedhart, Joachim. (2019, August 29). Material related to the blog "Dataviz with Flying Colors". Zenodo. http://doi.org/10.5281/zenodo.3381072

# Reorder colors to make colorblind friendly!
Tol_muted2 <- c('#882255','#DDCC77', '#44AA99', '#117733','#88CCEE', '#332288', '#DDDDDD', '#AA4499', '#CC6677', '#999933')


#### Land Richness vs Virus prevalence figure

# Marginal effects from best model with Land Richness on the x axis
DWV_landrich <- ggpredict(mod8c, c("Land_Richness_1000_z"))
BQCV_landrich <- ggpredict(mod9c, c("Land_Richness_1000_z"))
SBV_landrich <- ggpredict(mod10c, c("Land_Richness_1000_z"))

# add virus type identifier to each dataset
DWV_landrich$Virus <- "DWV"
BQCV_landrich$Virus <- "BQCV"
SBV_landrich$Virus <- "SBV"

virus_landrich <-rbind(DWV_landrich, BQCV_landrich, SBV_landrich)

plot(BQCV_landrich, add.data = F, jitter = 0.05) # initial plot (aim to replicate this)

# recalculate original landscape richness 
virus_landrich$x # scaled land richness
landrich_mean <- mean(disease$Land_Richness_1000) # mean of original land richness from disease data set
landrich_sd <- sd(disease$Land_Richness_1000) # sd of original land richness from disease data set

virus_landrich$LandRichness <- t((t(virus_landrich$x) * landrich_sd) + landrich_mean)

# plot of Landscape Richness in 1000m vs virus prevalence by virus and by species
virus_landrich_plot <- ggplot(virus_landrich, aes(x = LandRichness, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  labs(tag = "a") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1.2) +
  scale_linetype_manual(values=c("solid", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = "Landscape Richness \\nwithin 1000m", y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
print(virus_landrich_plot)
ggsave("VirusPrev_LandRichness_predict.tiff", plot = virus_landrich_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")
# color blind check
cvdPlot(virus_landrich_plot)

#### Natural area vs Virus prevalence figure

# Marginal effects from best model with Nat.1000 on the x axis
DWV_nat <- ggpredict(mod8c, c("Nat.1000_z"))
BQCV_nat <- ggpredict(mod9c, c("Nat.1000_z"))
SBV_nat <- ggpredict(mod10c, c("Nat.1000_z"))

# add virus type identifier to each dataset
DWV_nat$Virus <- "DWV"
BQCV_nat$Virus <- "BQCV"
SBV_nat$Virus <- "SBV"

virus_nat <-rbind(DWV_nat, BQCV_nat, SBV_nat)

# recalculate original proportion of Natural area
virus_nat$x # scaled proportion of Natural area
nat_mean <- mean(asqrt(disease$Nat.1000)) # mean of original proportion of Natural area from disease data set
nat_sd <- sd(asqrt(disease$Nat.1000)) # sd of original proportion of Natural area from disease data set

virus_nat$nat.1000 <- (sin(t((t(virus_nat$x) * nat_sd) + nat_mean)))^2

# plot of Natural area in 1000m vs virus prevalence by virus and by species
virus_nat_plot <- ggplot(virus_nat, aes(x = nat.1000, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  labs(tag = "b") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1.2) +
  scale_linetype_manual(values=c("solid", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = "Prop. of Natural Area \\nwithin 1000m", y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=9, color="black"))
print(virus_nat_plot)
ggsave("VirusPrev_Nat.1000_predict.tiff", plot = virus_nat_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")



#### floral Richness vs Virus prevalence figures

# Marginal effects from best model with floral Richness on the x axis
DWV_floralrich <- ggpredict(mod8c, c("floralrichness_z"))
BQCV_floralrich <- ggpredict(mod9c, c("floralrichness_z"))
SBV_floralrich <- ggpredict(mod10c, c("floralrichness_z"))

# add virus type identifier to each dataset
DWV_floralrich$Virus <- "DWV"
BQCV_floralrich$Virus <- "BQCV"
SBV_floralrich$Virus <- "SBV"

virus_floralrich <-rbind(DWV_floralrich, BQCV_floralrich, SBV_floralrich)

plot(BQCV_floralrich, add.data = F, jitter = 0.05) # initial plot (aim to replicate this)

# recalculate original floral richness 
virus_floralrich$x # scaled floral richness
floralrich_mean <- mean(disease$floralrichness) # mean of original floral richness from disease data set
floralrich_sd <- sd(disease$floralrichness) # sd of original floral richness from disease data set

virus_floralrich$floralrichness <- t((t(virus_floralrich$x) * floralrich_sd) + floralrich_mean)

# plot of floral richness vs virus prevalence by virus and by species
virus_floralrich_plot <- ggplot(virus_floralrich, aes(x = floralrichness, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +  
  labs(tag = "c") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1.2) +
  scale_linetype_manual(values=c("dashed", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = "Floral Richness", y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
print(virus_floralrich_plot)
ggsave("VirusPrev_floralrichness_predict.tiff", plot = virus_floralrich_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")



#### floral Density vs Virus prevalence figures

# Marginal effects from best model with floral density on the x axis
DWV_floralden <- ggpredict(mod8c, c("totaldensity_z"))
BQCV_floralden <- ggpredict(mod9c, c("totaldensity_z"))
SBV_floralden <- ggpredict(mod10c, c("totaldensity_z"))

# add virus type identifier to each dataset
DWV_floralden$Virus <- "DWV"
BQCV_floralden$Virus <- "BQCV"
SBV_floralden$Virus <- "SBV"

virus_floralden <-rbind(DWV_floralden, BQCV_floralden, SBV_floralden)

# recalculate original floral density 
virus_floralden$x # scaled floral density
floralden_mean <- mean(disease$totaldensity) # mean of original floral density from disease data set
floralden_sd <- sd(disease$totaldensity) # sd of original floral density from disease data set

virus_floralden$floraldensity <- t((t(virus_floralden$x) * floralden_sd) + floralden_mean)

# plot of floral density vs virus prevalence by virus and by species
virus_floralden_plot <- ggplot(virus_floralden, aes(x = floraldensity, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +  
  labs(tag = "d") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1.2) +
  scale_linetype_manual(values=c("solid", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = bquote("Floral Density (" ~m^2 ~ ")"), y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
print(virus_floralden_plot)
ggsave("VirusPrev_floralden_predict.tiff", plot = virus_floralden_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")



#### Pollinator Richness vs Virus prevalence figures

# Marginal effects from best model with Richness on the x axis
DWV_rich <- ggpredict(mod8c, c("Richness_zsqrt"))
BQCV_rich <- ggpredict(mod9c, c("Richness_zsqrt"))
SBV_rich <- ggpredict(mod10c, c("Richness_zsqrt"))

# add virus type identifier to each dataset
DWV_rich$Virus <- "DWV"
BQCV_rich$Virus <- "BQCV"
SBV_rich$Virus <- "SBV"

virus_rich <-rbind(DWV_rich, BQCV_rich, SBV_rich)


# recalculate original richness 
virus_rich$x # scaled richness
rich_mean <- mean(sqrt(disease$Richness)) # mean of original richness from disease data set
rich_sd <- sd(sqrt(disease$Richness)) # sd of original richness from disease data set

virus_rich$Richness <- (t((t(virus_rich$x) * rich_sd) + rich_mean))^2 

# plot of Richness vs virus prevalence by virus and by species
virus_rich_plot <- ggplot(virus_rich, aes(x = Richness, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +  
  labs(tag = "e") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1.2) +
  scale_linetype_manual(values=c("solid", "solid", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = "Species Richness", y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
print(virus_rich_plot)
ggsave("VirusPrev_Richness_predict.tiff", plot = virus_rich_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")



#### Apis + Bombus Abundance vs Virus prevalence figures

# Marginal effects from best model with Apis + Bombus Abundance on the x axis
DWV_ABabund <- ggpredict(mod8c, c("ApisBomb_ABUND_z"))
BQCV_ABabund <- ggpredict(mod9c, c("ApisBomb_ABUND_z"))
SBV_ABabund <- ggpredict(mod10c, c("ApisBomb_ABUND_z"))

# add virus type identifier to each dataset
DWV_ABabund$Virus <- "DWV"
BQCV_ABabund$Virus <- "BQCV"
SBV_ABabund$Virus <- "SBV"

virus_ABabund <-rbind(DWV_ABabund, BQCV_ABabund, SBV_ABabund)


# recalculate original Apis + Bombus abundance
virus_ABabund$x # scaled Apis + Bombus abundance
ABabund_mean <- mean(disease$APME_ABUND + disease$BOIM_ABUND) # mean of original Apis + Bombus abundance from disease data set
ABabund_sd <- sd(disease$APME_ABUND + disease$BOIM_ABUND) # sd of original Apis + Bombus abundance from disease data set

virus_ABabund$ApisBombAbundance <- t((t(virus_ABabund$x) * ABabund_sd) + ABabund_mean)

# plot of Apis + Bombus Abundance vs virus prevalence by virus and by species
virus_ABabund_plot <- ggplot(virus_ABabund, aes(x = ApisBombAbundance, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +  
  labs(tag = "f") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1) +
  scale_linetype_manual(values=c("dashed", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = bquote(italic(Apis)~"+"~italic(Bombus)~ "Abundance"), y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
print(virus_ABabund_plot)
ggsave("VirusPrev_ApisBombAbundance_predict.tiff", plot = virus_ABabund_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")




#### Pollinator Abundance vs Virus prevalence figures

# Marginal effects from best model with Abundance on the x axis
DWV_abund <- ggpredict(mod8, c("Abundance_z"))
BQCV_abund <- ggpredict(mod9, c("Abundance_z"))
SBV_abund <- ggpredict(mod10, c("Abundance_z"))

# add virus type identifier to each dataset
DWV_abund$Virus <- "DWV"
BQCV_abund$Virus <- "BQCV"
SBV_abund$Virus <- "SBV"

virus_abund <-rbind(DWV_abund, BQCV_abund, SBV_abund)


# recalculate original abundance
virus_abund$x # scaled abundance
abund_mean <- mean(disease$Abundance) # mean of original abundance from disease data set
abund_sd <- sd(disease$Abundance) # sd of original abundance from disease data set

virus_abund$Abundance <- t((t(virus_abund$x) * abund_sd) + abund_mean)

# plot of Abundance vs virus prevalence by virus and by species
virus_abund_plot <- ggplot(virus_abund, aes(x = Abundance, y = predicted, color = Virus)) +
  scale_fill_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +
  scale_color_manual(values = Tol_muted2, name = "Virus", labels = c("BQCV", "DWV", "SBV")) +  
  labs(tag = "g") +
  geom_line(aes(color = Virus, linetype = Virus), size = 1) +
  scale_linetype_manual(values=c("dashed", "dashed", "dashed"), guide = "none") +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = Virus), alpha = 0.2, outline.type = NULL) +
  labs(x = "Abundance", y = "Virus Prevalence") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.text = element_text(size=9, color = "black"), axis.title.y = element_text(size=11, color="black"), 
        axis.title.x = element_text(size=11, color="black"))
print(virus_abund_plot)
ggsave("VirusPrev_Abundance_predict.tiff", plot = virus_abund_plot, dpi = 600, width = 5, height = 5, units = "in", compression="lzw")




# Grid of all virus prevalence plots for all main factors in the SEM
virus_prev <- grid_arrange_shared_legend(virus_landrich_plot, virus_nat_plot, virus_floralrich_plot, virus_floralden_plot, virus_rich_plot, virus_ABabund_plot, nrow=2, ncol = 3, position = "right")
ggsave("Fig2.tiff", plot = virus_prev, dpi = 600, width = 8, height = 5, units = "in", compression="lzw")


# color blind check
cvdPlot(virus_prev)


###===================================================================================================== 
# APPENDIX S2 Table S6
# SEM with estimated richness at the average number of individuals across all sites (338 individuals)

# Estimated richness at 338 individuals was previously calculated in code found in https://doi.org/10.5061/dryad.zpc866t7g


# comparison of Estimated Richness transformations
mod1c <- lm(EstRich_338_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1c)
vif(mod1c) # all VIFs are less than 2
par(mfrow = c(2, 2))
plot(mod1c)


mod1d <- lm(EstRich_338_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1d)
vif(mod1d) # all VIFs are less than 2
plot(mod1d) # best fit


mod1e <- lm(EstRich_338_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ])
summary(mod1e)
vif(mod1e)  # all VIFs are less than 2
plot(mod1e)

#all fits are normal, used the log transformation since it seems to have the best fit


#normality
ols_plot_resid_qq(mod1c)
ols_plot_resid_qq(mod1d)
ols_plot_resid_qq(mod1e)

ols_test_normality(mod1c)
ols_test_normality(mod1d)
ols_test_normality(mod1e)

#Homogeneity of variance
ols_test_breusch_pagan(mod1c)
ols_test_breusch_pagan(mod1d)
ols_test_breusch_pagan(mod1e)



SEM_estRich <- psem(
  lm(EstRich_338_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% EstRich_338_zlog
)
summary(SEM_estRich, conserve = TRUE)


# Appendix S2 Table S12
std_range_coefs_estrich <- stdCoefs(SEM_estRich, data=disease, 
                                    standardize="range",
                                    standardize.type = "latent.linear",
                                    intercepts = F)

write.csv(std_range_coefs_estrich, "sem_std_range_coefs_estrichness_at_338.csv", row.names=FALSE)



## Test of global fit of estimated richness model by removing some non-significant links from the model
SEM_estRich_test <- psem(
  lm(EstRich_338_zlog ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = disease[1:13, ]),
  glmer(DWV ~ Land_Richness_1000_z + floralrichness_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  glmer(SBV ~ Nat.1000_z + totaldensity_z + EstRich_338_zlog + ApisBomb_ABUND_z + (1|Genus.Code) + (1|Site) + (1|Site:Visit), family = binomial, data = disease),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% EstRich_338_zlog
)
summary(SEM_estRich_test, conserve = TRUE)
std_range_coefs_estRich_test <- stdCoefs(SEM_estRich_test, data=disease, 
                                                   standardize="range",
                                                   standardize.type = "latent.linear",
                                                   intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit




###===================================================================================================== 
# separate models for each host species with all three viruses
# top results for each host species are shown in Figure 3


## APIS specific models with all three viruses

# model with total abundance 
SEM_apis <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  Abundance_z %~~% Richness_zsqrt
)
summary(SEM_apis, conserve = TRUE)


# APPENDIX S3 Table S1
# model with combined apis and bombus abundance [Shown in the main text as Figure 3A]
SEM_apis2 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_apis2, conserve = TRUE)

# model with only Apis abundance
SEM_apis3 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  APME_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_apis3, conserve = TRUE)


# model with only Bombus abundance
SEM_apis4 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  BOIM_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_apis4, conserve = TRUE)

# AIC comparison of all four APIS SEMs
AIC(SEM_apis, SEM_apis2, SEM_apis3, SEM_apis4)  # SEM_apis2 has the lowest AIC, including the ApisBomb_ABUND. This is the model we present in the main text. 

# APPENDIX S3 Table S1
std_range_coefs_apis <- stdCoefs(SEM_apis2, data=disease, 
                                 standardize="range",
                                 standardize.type = "latent.linear",
                                 intercepts = F)
write.csv(std_range_coefs_apis, "sem_std_range_coefs_APIS.csv", row.names=FALSE)


# correlation matrix
subset_APIS <- diseaseAPIS %>%
  select(DWV, BQCV, SBV, Richness_zsqrt, Abundance_z, APME_ABUND_z, BOIM_ABUND_z, ApisBomb_ABUND_z, LASI_ABUND_z, PEPR_ABUND_z, Land_Richness_1000_z, Nat.1000_z, floralrichness_z, totaldensity_z)

correlations_APIS <- cor(subset_APIS)
correlations_APIS

correlations_APIS <- format(round(correlations_APIS, 3))
# set diagonal to NA
diag(correlations_APIS) <- "--"
rownames(correlations_APIS) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(correlations_APIS) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(correlations_APIS, "APIS_SEM_factor_correlations.csv", quote = F, row.names = T)



# Check for model assumptions

# Richness and Abundance models are the same as used above because they only include the site-level variables, which are consistent here.

apismod3 <- glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS)
summary(apismod3)
vif(apismod3) # VIFs < 4.2
plot(apismod3)
overdisp_fun(apismod3)

apismod4 <- glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS)
summary(apismod4)
vif(apismod4)  # VIFs < 4.5
plot(apismod4)
overdisp_fun(apismod4)

apismod5 <- glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS)
summary(apismod5)
vif(apismod5)  # VIFs < 4.1
plot(apismod5)
overdisp_fun(apismod5)



## Test of global fit of Apis specific model by removing some non-significant links from the model
SEM_apis2_test <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseAPIS[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(BQCV ~ Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  glmer(SBV ~ Land_Richness_1000_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseAPIS),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_apis2_test, conserve = TRUE)
std_range_coefs_APIS_test <- stdCoefs(SEM_apis2_test, data=disease, 
                                         standardize="range",
                                         standardize.type = "latent.linear",
                                         intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit


#####==============================================================================================
## BOMBUS specific model with all three viruses

# model with total abundance 
SEM_bombus <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  Abundance_z %~~% Richness_zsqrt
)
summary(SEM_bombus, conserve = TRUE)


# APPENDIX S3 Table S2
# model with combined bomb and bombus abundance [Shown in the main text as Figure 3A]
SEM_bomb2 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_bomb2, conserve = TRUE)

# model with only bomb abundance
SEM_bomb3 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  APME_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_bomb3, conserve = TRUE)


# model with only Bombus abundance
SEM_bomb4 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  BOIM_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_bomb4, conserve = TRUE)

# AIC comparison of all four bomb SEMs
AIC(SEM_bombus, SEM_bomb2, SEM_bomb3, SEM_bomb4)  # SEM_bomb2 has the lowest AIC, including the ApisBomb_ABUND. This is the model we present in the main text. 


# APPENDIX S3 Table S2
std_range_coefs_bombus <- stdCoefs(SEM_bomb2, data=disease, 
                                 standardize="range",
                                 standardize.type = "latent.linear",
                                 intercepts = F)
write.csv(std_range_coefs_bombus, "sem_std_range_coefs_BOMBUS.csv", row.names=FALSE)


# correlation matrix
subset_BOMB <- diseaseBOMB %>%
  select(DWV, BQCV, SBV, Richness_zsqrt, Abundance_z, APME_ABUND_z, BOIM_ABUND_z, ApisBomb_ABUND_z, LASI_ABUND_z, PEPR_ABUND_z, Land_Richness_1000_z, Nat.1000_z, floralrichness_z, totaldensity_z)

correlations_BOMB <- cor(subset_BOMB)
correlations_BOMB

correlations_BOMB <- format(round(correlations_BOMB, 3))
# set diagonal to NA
diag(correlations_BOMB) <- "--"
rownames(correlations_BOMB) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(correlations_BOMB) <- c("DWV", "BQCV", "SBV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(correlations_BOMB, "BOMB_SEM_factor_correlations.csv", quote = F, row.names = T)


## Test of global fit of Bombus specific model by removing some non-significant links from the model
SEM_bomb2_test <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseBOMB[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  glmer(SBV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseBOMB),
  DWV %~~% BQCV,
  BQCV %~~% SBV,
  SBV %~~% DWV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_bomb2_test, conserve = TRUE)
std_range_coefs_bomb_test <- stdCoefs(SEM_bomb2_test, data=disease, 
                                      standardize="range",
                                      standardize.type = "latent.linear",
                                      intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit


#####==============================================================================================
## LASIOGLOSSUM specific model with DWV and BQCV viruses
# Lasioglossum rarely had SBV, so we removed that virus from the model.


# model with total abundance 
SEM_lasi <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  Abundance_z %~~% Richness_zsqrt
)
summary(SEM_lasi, conserve = TRUE)

# APPENDIX S3 Table S3
# model with Apis and Bombus abundance 
SEM_lasi2 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_lasi2, conserve = TRUE)


# model with only Apis abundance 
SEM_lasi3 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  APME_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_lasi3, conserve = TRUE)


# model with only Bombus abundance 
SEM_lasi4 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  BOIM_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_lasi4, conserve = TRUE)


# model with only Lasioglossum abundance 
SEM_lasi5 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(LASI_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + LASI_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  LASI_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_lasi5, conserve = TRUE)


# AIC comparison of all four Lasi SEMs
AIC(SEM_lasi, SEM_lasi2, SEM_lasi3, SEM_lasi4, SEM_lasi5)  # SEM_lasi2 has the lowest AIC, including the ApisBomb_ABUND. This is the model we present in the main text. 

# APPENDIX S3 Table S3
std_range_coefs_lasi <- stdCoefs(SEM_lasi2, data=disease, 
                                   standardize="range",
                                   standardize.type = "latent.linear",
                                   intercepts = F)
write.csv(std_range_coefs_lasi, "sem_std_range_coefs_LASI.csv", row.names=FALSE)


# correlation matrix
subset_LASI <- diseaseLASI %>%
  select(DWV, BQCV, Richness_zsqrt, Abundance_z, APME_ABUND_z, BOIM_ABUND_z, ApisBomb_ABUND_z, LASI_ABUND_z, PEPR_ABUND_z, Land_Richness_1000_z, Nat.1000_z, floralrichness_z, totaldensity_z)

correlations_LASI <- cor(subset_LASI)
correlations_LASI

correlations_LASI <- format(round(correlations_LASI, 3))
# set diagonal to NA
diag(correlations_LASI) <- "--"
rownames(correlations_LASI) <- c("DWV", "BQCV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(correlations_LASI) <- c("DWV", "BQCV", "Species Richness", "Abundance", "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(correlations_LASI, "LASI_SEM_factor_correlations.csv", quote = F, row.names = T)


## Test of global fit of Lasioglossum specific model by removing some non-significant links from the model
SEM_lasi2_test <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseaseLASI[1:13, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseaseLASI),
  DWV %~~% BQCV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_lasi2_test, conserve = TRUE)
std_range_coefs_lasi_test <- stdCoefs(SEM_lasi2_test, data=disease, 
                                      standardize="range",
                                      standardize.type = "latent.linear",
                                      intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit




#####==============================================================================================
## EUCERA specific model with DWV and BQCV viruses
# Eucera pruinosa rarely had SBV, so we removed that virus from the model
# one site did not have any Eucera pruinosa detected, so that site is not included in these models.

# model with total abundance
SEM_pepo <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(Abundance_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + Abundance_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  DWV %~~% BQCV,
  Abundance_z %~~% Richness_zsqrt
)
summary(SEM_pepo, conserve = TRUE)

# APPENDIX S3 Table S4
# model with Apis and Bombus abundance 
SEM_pepo2 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  DWV %~~% BQCV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_pepo2, conserve = TRUE)


# model with only Apis abundance 
SEM_pepo3 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(APME_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + APME_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  DWV %~~% BQCV,
  APME_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_pepo3, conserve = TRUE)


# model with only Bombus abundance 
SEM_pepo4 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(BOIM_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + BOIM_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))),
  DWV %~~% BQCV,
  BOIM_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_pepo4, conserve = TRUE)


# model with only E. pruinosa abundance 
SEM_pepo5 <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(PEPR_ABUND_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_zsqrt + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z + Richness_zsqrt + PEPR_ABUND_zsqrt + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))),
  DWV %~~% BQCV,
  PEPR_ABUND_zsqrt %~~% Richness_zsqrt
)
summary(SEM_pepo5, conserve = TRUE)


# AIC comparison of all four Lasi SEMs
AIC(SEM_pepo, SEM_pepo2, SEM_pepo3, SEM_pepo4, SEM_pepo5)  # SEM_pepo2 has the lowest AIC, including the ApisBomb_ABUND. This is the model we present in the main text. 

# APPENDIX S3 Table S4
std_range_coefs_pepo <- stdCoefs(SEM_pepo2, data=disease, 
                                 standardize="range",
                                 standardize.type = "latent.linear",
                                 intercepts = F)
write.csv(std_range_coefs_pepo, "sem_std_range_coefs_PEPO.csv", row.names=FALSE)



# correlation matrix
subset_PEPO <- diseasePEPO %>%
  select(DWV, BQCV, Richness_zsqrt, Abundance_zsqrt, APME_ABUND_z, BOIM_ABUND_z, ApisBomb_ABUND_z, LASI_ABUND_z, PEPR_ABUND_z, Land_Richness_1000_z, Nat.1000_z, floralrichness_z, totaldensity_z)

correlations_PEPO <- cor(subset_PEPO)
correlations_PEPO

correlations_PEPO <- format(round(correlations_PEPO, 3))
# set diagonal to NA
diag(correlations_PEPO) <- "--"
rownames(correlations_PEPO) <- c("DWV", "BQCV", "Species Richness", "Abundance",  "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
colnames(correlations_PEPO) <- c("DWV", "BQCV", "Species Richness", "Abundance",  "A. mellifera Abundance", "B. impatiens Abundance", "Apis + Bombus Abundance", "Lasioglossum Abundance", "E. pruinosa Abundance", "Landscape Richness", "Natural Area", "Floral Richness", "Floral Density")
write.csv(correlations_PEPO, "PEPO_SEM_factor_correlations.csv", quote = F, row.names = T)


## Test of global fit of Eucera specific model by removing some non-significant links from the model
SEM_pepo2_test <- psem(
  lm(Richness_zsqrt ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  lm(ApisBomb_ABUND_z ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + totaldensity_z, data = diseasePEPO[1:12, ]),
  glmer(DWV ~ Nat.1000_z + Land_Richness_1000_z + totaldensity_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  glmer(BQCV ~ Nat.1000_z + Land_Richness_1000_z + floralrichness_z + Richness_zsqrt + ApisBomb_ABUND_z + (1|Site) + (1|Site:Visit), family = binomial, data = diseasePEPO),
  DWV %~~% BQCV,
  ApisBomb_ABUND_z %~~% Richness_zsqrt
)
summary(SEM_pepo2_test, conserve = TRUE)
std_range_coefs_pepo_test <- stdCoefs(SEM_pepo2_test, data=disease, 
                                      standardize="range",
                                      standardize.type = "latent.linear",
                                      intercepts = F)
# topography and magnitude of paths does not change, but X-squared and Fisher's C both have p values > 0.05 indicating good model fit

