## 08 - Exploration of model performance

# Libraries ---------------------------------------------------------------

library(reshape2)
library(lme4)
library(dplyr)
library(stringr)

# Read data ---------------------------------------------------------------

PRresults <- readRDS(file = "./Outputs/PRROCresults.Rds")
samples <- readRDS(file = "./Outputs/SampleSummaryData.Rds")

prres <- bind_rows(PRresults)
ss_full <- bind_rows(samples)

# Reformat data -----------------------------------------------------------

ss_melt <- melt(ss_full)
ss_melt$variable <- as.character(ss_melt$variable)
unique(ss_melt$variable)
ss_melt$variable[ss_melt$variable == "Train_FlickrSuppSize"] <- "FlickrSupp"
ss_melt$variable[ss_melt$variable == "Train_FlickrRanSuppSize"] <- "FlickrSuppRan"
ss_melt$variable[ss_melt$variable == "Train_GBIFGeoSize"] <- "GBIF"
ss_melt$variable[ss_melt$variable == "Train_GBIFRanSize"] <- "GBIFran"

ss_melt <- rename(ss_melt, "train" = variable, "sample" = value)
ss_melt <- ss_melt[!str_detect(ss_melt$train, "_"),]

prres_modelDF <- left_join(prres, ss_melt)

unique(prres_modelDF$train)
unique(prres_modelDF$test_list)

# remove instances where the test data set should not be used on certain SDMs
# training with overlapping data
prres_modelDF <- prres_modelDF %>% 
  filter(!train == "FlickrSupp" | !test_list == "Flickr") %>% 
  filter(!train == "FlickrSuppRan" | !test_list == "Flickr") %>% 
  filter(!train == "FlickrSupp" | !test_list == "GBIFran") %>% 
  filter(!train == "FlickrSuppRan" | !test_list == "GBIF") %>% 
  filter(!train == "GBIFran" | !test_list == "GBIF") %>% 
  filter(!train == "GBIF" | !test_list == "GBIFran")
  
# Generate models targeting SMD performance -------------------------------

# PRRC values ----

mod1 <- lm(prAUC ~ train, data = prres_modelDF)
mod2 <- lm(prAUC ~ species, data = prres_modelDF)
mod3 <- lm(prAUC ~ sample, data = prres_modelDF)
mod4 <- lm(prAUC ~ species*sample, data = prres_modelDF)
mod5 <- lmer(prAUC ~ train-1 + (1|species), data = prres_modelDF)
mod6 <- lmer(prAUC ~ sample-1 + (1|species), data = prres_modelDF)
mod7 <- lmer(prAUC ~ (0 + train|species), data = prres_modelDF)
mod8 <- lmer(prAUC ~ (0 + sample|species), data = prres_modelDF)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)$AIC

modelset1 <- data.frame(Model = c("prAUC ~ train",
                     "prAUC ~ species",
                     "prAUC ~ sample",
                     "prAUC ~ species*sample",
                     "prAUC ~ train-1 + (1|species)",
                     "prAUC ~ sample-1 + (1|species)",
                     "prAUC ~ (0 + train|species)",
                     "prAUC ~ (0 + sample|species)"), 
           AIC = AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)$AIC)

summary(mod7)

ranef(mod7)
fixef(mod7)

qqnorm(residuals(mod7))

# not normal distribution
shapiro.test(residuals(mod7))

# ROC values ----

mod1 <- lm(rocAUC ~ train, data = prres_modelDF)
mod2 <- lm(rocAUC ~ species, data = prres_modelDF)
mod3 <- lm(rocAUC ~ sample, data = prres_modelDF)
mod4 <- lm(rocAUC ~ species*sample, data = prres_modelDF)
mod5 <- lmer(rocAUC ~ train-1 + (1|species), data = prres_modelDF)
mod6 <- lmer(rocAUC ~ sample-1 + (1|species), data = prres_modelDF)
mod7 <- lmer(rocAUC ~ (0 + train|species), data = prres_modelDF)
mod8 <- lmer(rocAUC ~ (0 + sample|species), data = prres_modelDF)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)

modelset2 <- data.frame(Model = c("rocAUC ~ train",
                                  "rocAUC ~ species",
                                  "rocAUC ~ sample",
                                  "rocAUC ~ species*sample",
                                  "rocAUC ~ train-1 + (1|species)",
                                  "rocAUC ~ sample-1 + (1|species)",
                                  "rocAUC ~ (0 + train|species)",
                                  "rocAUC ~ (0 + sample|species)"), 
                        AIC = AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)$AIC)

summary(mod7)

ranef(mod7)
fixef(mod7)

qqnorm(residuals(mod7))
# not normal distribution
shapiro.test(residuals(mod7))

# Difference in ROC and PRRC values ----

flickrSuppValues <- prres_modelDF %>% 
  filter(str_detect(train, "Flickr"))

nonflickrSuppValues <- prres_modelDF %>% 
  filter(!str_detect(train, "Flickr")) %>% 
  filter(!train == "GBIF" | test_list == "GBIF") %>% 
  filter(!train == "GBIFran" | test_list == "GBIFran")

diff_prAUC <- flickrSuppValues$prAUC - nonflickrSuppValues$prAUC
diff_rocAUC <- flickrSuppValues$rocAUC - nonflickrSuppValues$rocAUC
diff_sample <- flickrSuppValues$sample - nonflickrSuppValues$sample

diff_modeldf <- data.frame(diff_prAUC, diff_rocAUC, diff_sample,
                           "species" = flickrSuppValues$species,
                           "mod_reg" = flickrSuppValues$model)

mod1 <- lm(diff_prAUC ~ diff_sample, data = diff_modeldf)
mod2 <- lm(diff_prAUC ~ species, data = diff_modeldf)
mod3 <- lm(diff_prAUC ~ mod_reg, data = diff_modeldf)
mod4 <- lm(diff_prAUC ~ diff_sample + species, data = diff_modeldf)
mod5 <- lmer(diff_prAUC ~ diff_sample-1 + (1|species), data = diff_modeldf)
mod6 <- lmer(diff_prAUC ~ mod_reg-1 + (1|species), data = diff_modeldf)
mod7 <- lmer(diff_prAUC ~ (0 + diff_sample|species), data = diff_modeldf)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7)
modelset3 <- data.frame(Model = c("diff_prAUC ~ diff_sample",
                                  "diff_prAUC ~ species",
                                  "diff_prAUC ~ mod_reg",
                                  "diff_prAUC ~ diff_sample + species",
                                  "diff_prAUC ~ diff_sample-1 + (1|species)",
                                  "diff_prAUC ~ mod_reg-1 + (1|species)",
                                  "diff_prAUC ~ (0 + diff_sample|species)"), 
                        AIC = AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7)$AIC)

summary(mod2)
summary(mod4)

qqnorm(residuals(mod2))
# not normal distribution
shapiro.test(residuals(mod2))

mod1 <- lm(diff_rocAUC ~ diff_sample, data = diff_modeldf)
mod2 <- lm(diff_rocAUC ~ species, data = diff_modeldf)
mod3 <- lm(diff_rocAUC ~ mod_reg, data = diff_modeldf)
mod4 <- lm(diff_rocAUC ~ diff_sample + species, data = diff_modeldf)
mod5 <- lmer(diff_rocAUC ~ diff_sample-1 + (1|species), data = diff_modeldf)
mod6 <- lmer(diff_rocAUC ~ mod_reg-1 + (1|species), data = diff_modeldf)
mod7 <- lmer(diff_rocAUC ~ (0 + diff_sample|species), data = diff_modeldf)

AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7)
modelset4 <- data.frame(Model = c("diff_rocAUC ~ diff_sample",
                                  "diff_rocAUC ~ species",
                                  "diff_rocAUC ~ mod_reg",
                                  "diff_rocAUC ~ diff_sample + species",
                                  "diff_rocAUC ~ diff_sample-1 + (1|species)",
                                  "diff_rocAUC ~ mod_reg-1 + (1|species)",
                                  "diff_rocAUC ~ (0 + diff_sample|species)"), 
                        AIC = AIC(mod1, mod2, mod3, mod4, mod5, mod6, mod7)$AIC)

summary(mod2)
summary(mod4)

qqnorm(residuals(mod2))
# not normal distribution
shapiro.test(residuals(mod2))

fullmodellist <- rbind(modelset1, modelset2, modelset3, modelset4)

write.csv(file = "./Tables/SUPPLEMENTARY TABLE - Full model specification and AIC.csv",
          row.names = FALSE, x = fullmodellist)
