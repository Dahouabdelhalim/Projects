## Meta-Analysis Appendix Code 

## Except mean-var relationships, which were checked in main code and phylogenetic correlograms which were made when examining phylogenetic signal.

remove(list = ls())


## clear console and plots

library(metafor)
library(broom)
library(ggtext)
library(ggbeeswarm)
library(ggpubr)
library(flextable)
library(scales)
library(ape)
library(rotl)
library(phylosignal)
library(adephylo)
library(phylobase)
library(tidyverse)

## set wd


## read in data

df_effects <- 

## read in phylogenies - different for each effect size, because of different numbers of species

load("phylogeny_matrix_meanN") ## nitrogen mean difference
load("phylogeny_matrix_meanC") ## carbon mean difference
load("phylogeny_matrix_varN")  ## nitrogen variation difference
load("phylogeny_matrix_varC")  ## carbon variation difference
load("phylogeny_matrix_carn")  ## non-gape-limited carnivores
load("phylogeny_matrix_gape")  ## gape-limited carnivores

## create Species.phylo variable to combine with phylogeny as a random factor in model - distinct from species alone as a random factor

df_effects <- df_effects %>%
  mutate(Species.phylo = Species)

## one row with n = 1 for both sexes creates problems later, as sampling variance is non-positive, so cannot calculate I2 etc. So, drop this row now.

df_effects <- df_effects %>% filter(Species != "Micronycteris schmidtorum")


#### Non-Phylogenetic Gape-Limitation ####

## our finding of a greater effect of size dimorphism in gape-limited predators was smaller when controlling for phylogeny. So this is a non-phylogenetic model for comparisson. 

mixed_carn <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                           mods = ~SSD,
                           random = list(~ 1 |                                             Paper_Number,
                                         ~ 1 | Species),
                           data=df_effects %>% 
                             filter(Diet == "Carnivore") %>%
                             filter(Gape_Lim == "No"),
                           method="REML")

summary(mixed_carn, digits = 3)


flextable(tidy(mixed_carn))


mixed_gape <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                             mods = ~SSD,
                             random = list(~ 1 |                                             Paper_Number, ~ 1 | Species),
                             data=df_effects %>% 
                               filter(Diet == "Carnivore") %>% 
                               filter(Gape_Lim == "Yes"),
                             method="REML")

summary(mixed_gape)

flextable(tidy(mixed_gape))

## plotting

newdata_gape <- tibble(SSD = seq(-9, 6, length = 1000)) 

## NGL carnivores

pred_carn <- as_tibble(predict(mixed_carn, newmods = cbind(newdata_gape$SSD)))

df_pred_carn <- bind_cols(newdata_gape, pred_carn)

## GL Carnivores phylogenetic control

pred_gape <- as_tibble(predict(mixed_gape, newmods = cbind(newdata_gape$SSD)))

df_pred_gape <- bind_cols(newdata_gape, pred_gape)

## for color and fill aesthetics/ legend, need to add group variables to df_effects and the pred df's here

## for df effects, need the snake species to be labelled "snakes" and then use subphylum for other names, which will leave fish labelled as fish

df_effects <- df_effects %>% 
  mutate(Group = 
           if_else(Species == "Laticauda saintgironsi", "Snakes",
                   if_else(Species == "Nerodia rhombifer", "Snakes",
                           if_else(Species == "Nerodia erythrogaster", "Snakes",
                                   if_else(Species == "Nerodia sipedon", "Snakes", Subphylum)))))

## for the pred df's, need to label group as either Non-Gape-Limited Carnivores or Gape-Limited Carnivores

df_pred_carn <- df_pred_carn %>% 
  mutate(Group = "Non-Gape-Limited\\nCarnivores")

df_pred_gape <- df_pred_gape %>% 
  mutate(Group = "Gape-Limited\\nCarnivores")


plot_scatter_gape_nophylo <- df_pred_gape %>% filter(SSD >= -4.5 & SSD <= 4.5) %>% 
  ggplot(aes(x = SSD, y = pred)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_ribbon(data = df_pred_carn %>%
                filter(SSD >= -4.5 & SSD <= 4.5),
              aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.6, fill = "grey25") +
  geom_line(data = df_pred_carn %>% 
              filter(SSD >= -4.5 & SSD <= 4.5),
            aes(col = Group),
            size = 1.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.3, fill = "#762a83", col = "#762a83") +
  geom_point(data = df_effects %>% 
               filter(Diet == "Carnivore") %>% 
               filter(Gape_Lim == "Yes"),
             aes(x = SSD, y = Mean_DiffN, 
                 fill = Group),
             alpha = 0.8, size = 5, stroke = 2, shape = 21) +
  geom_line(data = df_pred_gape %>% 
              filter(SSD >= -4.5 & SSD < 0.935), 
            size = 2, 
            aes(col = Group)) +
  geom_line(size = 2, linetype = "dashed",
            aes(col = Group)) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_color_manual(labels = c("Gape-Limited\\nCarnivores", "Non-Gape-Limited\\nCarnivores"),
                     breaks = c("Gape-Limited\\nCarnivores", "Non-Gape-Limited\\nCarnivores"),
                     values = c("#762a83", "grey25")) +
  scale_fill_manual(labels = c("Fish", "Snakes"),
                    breaks = c("Fish", "Snakes"),
                    values = c("#d8b365", "#1b7837")) +
  theme(axis.title = element_text(size = 25, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "right") +
  labs(x = "Sexual Dimorphism",
       y = "Mean \\u03b415N Sex Difference (\\u2030)",
       color = "Gape Type",
       fill = "Group")

ggsave("scatter_gape_nophylo.jpeg", units = "mm", height = 200, width = 300)

#### PUBLICATION YEAR ####

## models with publication year as a predictor show tiny effect sizes for mean and variance sex differences, for both isotopes, thus no evidence of publication bias.

mixed_pubyear <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                        mods = ~Year,
                        random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                      ~ 1 | Species.phylo),
                        R = list(Species.phylo = phylo_cor_meanN),
                        data=df_effects %>% drop_na(Mean_DiffN),
                        method="REML")

residuals_mixed_pubyear <- residuals.rma(mixed_pubyear)

hist(residuals_mixed_pubyear, breaks = 100)

summary(mixed_pubyear)

mixed_lnVR_pubyear <- rma.mv(lnVRN, lnVRN.sv, 
                             mods = ~Year,
                             random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                           ~ 1 | Species.phylo),
                             R = list(Species.phylo = phylo_cor_varN),
                             data=df_effects %>% drop_na(lnVRN),
                             method="REML")

residuals_mixed_lnVR_pubyear <- residuals.rma(mixed_lnVR_pubyear)

hist(residuals_mixed_lnVR_pubyear, breaks = 100)

summary(mixed_lnVR_pubyear)

mixed_pubyearC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                         mods = ~Year,
                         random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_meanC),
                         data=df_effects %>% drop_na(Mean_DiffC),
                         method="REML")

residuals_mixed_pubyearC <- residuals.rma(mixed_pubyearC)

hist(residuals_mixed_pubyearC, breaks = 100)

summary(mixed_pubyearC)


mixed_lnVR_pubyearC <- rma.mv(lnVRC, lnVRC.sv, 
                              mods = ~Year,
                              random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                            ~ 1 | Species.phylo),
                              R = list(Species.phylo = phylo_cor_varC),
                              data=df_effects %>% drop_na(lnVRC),
                              method="REML")

residuals_mixed_lnVR_pubyearC <- residuals.rma(mixed_lnVR_pubyearC)

hist(residuals_mixed_lnVR_pubyearC, breaks = 100)

summary(mixed_lnVR_pubyearC)

df_pubyear <- data.frame(Isotope = c("Nitrogen", "Nitrogen", "Carbon", "Carbon"),
                         Measure = c("Mean Difference", "Variation", "Mean Difference", "Variation"),
                         Effect_PubYear = as.numeric(c(-0.001, -0.0035, 0.0066, -0.0094)),
                         se = as.numeric(c(0.0100, 0.0066, 0.0081, 0.0073)),
                         ci.lb = as.numeric(c(-0.0221, -0.0164, -0.0094, -0.0237)),
                         ci.ub = as.numeric(c(0.0202, 0.0095, 0.0225, 0.0049)),
                         pval = as.numeric(c(0.927, 0.6010, 0.4195, 0.1979)))

flextable(df_pubyear)

#### FUNNEL PLOTS ####

MA_Mean_DiffN <- rma(Mean_DiffN, Mean_DiffN.sv, data=df_effects)

MA_Mean_DiffC <- rma(Mean_DiffC, Mean_DiffC.sv, 
                     data=df_effects)

MA_lnVRN <- rma(lnVRN, lnVRN.sv, data=df_effects)

MA_lnVRC <- rma(lnVRC, lnVRC.sv, data=df_effects)


jpeg(filename = "funnel_mean_diffN.jpeg", res = 300, height = 125, width = 225, units = "mm")

funnel(MA_Mean_DiffN)

dev.off()

jpeg(filename = "funnel_lnVRN.jpeg", res = 300, height = 125, width = 225, units = "mm")

funnel(MA_lnVRN)

dev.off()

jpeg(filename = "funnel_mean_diffC.jpeg", res = 300, height = 125, width = 225, units = "mm")

funnel(MA_Mean_DiffC)

dev.off()

jpeg(filename = "funnel_lnVRC.jpeg", res = 300, height = 125, width = 225, units = "mm")

funnel(MA_lnVRC)

dev.off()

#### WITHIN-STUDY VARIANCE ONLY ####

## default rma.mv study weighting accounts for within study variation, between study variation and random effects. However, high between study variation can "eat up" within study variation, possibly effecting model outcomes. Therefore reccomended to check model outcomes using only the inverse weighting method, i.e. 1/within study sampling variance as the weight given to each study = effectively a fixed effects model.

## nitrogen mean sex differences

mixed_inverse_weightingN <- rma.mv(Mean_DiffN, Mean_DiffN.sv,
                                  mods = ~SSD, W = 1/Mean_DiffN.sv,
                                  random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                ~ 1 | Species.phylo),
                                  R = list(Species.phylo = phylo_cor_meanN),
                                  data=df_effects %>% drop_na(Mean_DiffN),
                                  method="REML")

residuals_mixed_inverse_weightingN <- residuals.rma(mixed_inverse_weightingN)

hist(residuals_mixed_inverse_weightingN, breaks = 100)

summary(mixed_inverse_weightingN)

flextable(tidy(summary(mixed_inverse_weightingN)))

## nitrogen variation sex differences

mixed_inverse_weighting_lnVRN <- rma.mv(lnVRN, lnVRN.sv,
                                       mods = ~SSD, W = 1/lnVRN.sv,
                                       random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                     ~ 1 | Species.phylo),
                                       R = list(Species.phylo = phylo_cor_varN),
                                       data=df_effects %>% drop_na(lnVRN),
                                       method="REML")

residuals_mixed_inverse_weighting_lnVRN <- residuals.rma(mixed_inverse_weighting_lnVRN)

hist(residuals_mixed_inverse_weighting_lnVRN, breaks = 100)

summary(mixed_inverse_weighting_lnVRN)

flextable(tidy(summary(mixed_inverse_weighting_lnVRN)))

## carbon mean sex differences

mixed_inverse_weightingC <- rma.mv(Mean_DiffC, Mean_DiffC.sv,
                                   mods = ~SSD, W = 1/Mean_DiffC.sv,
                                   random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                 ~ 1 | Species.phylo),
                                   R = list(Species.phylo = phylo_cor_meanC),
                                   data=df_effects %>% drop_na(Mean_DiffC),
                                   method="REML")

residuals_mixed_inverse_weightingC <- residuals.rma(mixed_inverse_weightingC)

hist(residuals_mixed_inverse_weightingC, breaks = 100)

summary(mixed_inverse_weightingC)

flextable(tidy(summary(mixed_inverse_weightingC)))

## carbon variance

mixed_inverse_weighting_lnVRC <- rma.mv(lnVRC, lnVRC.sv,
                                        mods = ~SSD, W = 1/lnVRC.sv,
                                        random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                      ~ 1 | Species.phylo),
                                        R = list(Species.phylo = phylo_cor_varC),
                                        data=df_effects %>% drop_na(lnVRC),
                                        method="REML")

residuals_mixed_inverse_weighting_lnVRC <- residuals.rma(mixed_inverse_weighting_lnVRC)

hist(residuals_mixed_inverse_weighting_lnVRC, breaks = 100)

summary(mixed_inverse_weighting_lnVRC)

flextable(tidy(summary(mixed_inverse_weighting_lnVRC)))

## diet + size

mixed_ssd_diet_size_1intA <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                               mods = ~SSD * Diet + Size, W = 1/Mean_DiffN.sv,
                                    random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

residuals_mixed_ssd_diet_size_1intA_inverse_weighting <- residuals.rma(mixed_ssd_diet_size_1intA)

hist(residuals_mixed_ssd_diet_size_1intA_inverse_weighting, breaks = 282)

summary(mixed_ssd_diet_size_1intA)

flextable(tidy(summary(mixed_ssd_diet_size_1intA)))

## gape limited

mixed_inverse_weighting_gape <- rma.mv(Mean_DiffN, Mean_DiffN.sv,
                                       mods = ~SSD, W = 1/Mean_DiffN.sv,
                                       random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                     ~ 1 | Species.phylo),
                                       R = list(Species.phylo = phylo_cor_gape),
                                       data=df_effects %>% 
                                         filter(Diet == "Carnivore") %>% 
                                         filter(Gape_Lim == "Yes"),
                                       method="REML")

residuals_mixed_inverse_weighting_gape <- residuals.rma(mixed_inverse_weighting_gape)

hist(residuals_mixed_inverse_weighting_gape, breaks = 10)

summary(mixed_inverse_weighting_gape)

flextable(tidy(summary(mixed_inverse_weighting_gape)))

## carnivores only

mixed_inverse_weighting_carnivore <- rma.mv(Mean_DiffN, Mean_DiffN.sv,
                                            mods = ~SSD, W = 1/Mean_DiffN.sv,
                                            random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                          ~ 1 | Species.phylo),
                                            R = list(Species.phylo = phylo_cor_carn),
                                            data=df_effects %>% 
                                              filter(Diet == "Carnivore") %>%
                                              filter(Gape_Lim == "No"),
                                            method="REML")

residuals_mixed_inverse_weighting_carnivore <- residuals.rma(mixed_inverse_weighting_carnivore)

hist(residuals_mixed_inverse_weighting_carnivore, breaks = 100)

summary(mixed_inverse_weighting_carnivore)

flextable(tidy(summary(mixed_inverse_weighting_carnivore)))

## carbon mean differences ~ ssd* diet

mixed_inverse_weighting_dietC <- rma.mv(Mean_DiffC, Mean_DiffC.sv,
                                   mods = ~SSD * Diet, W = 1/Mean_DiffC.sv,
                                   random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                                 ~ 1 | Species.phylo),
                                   R = list(Species.phylo = phylo_cor_meanC),
                                   data=df_effects %>% drop_na(Mean_DiffC),
                                   method="REML")

residuals_mixed_inverse_weighting_dietC <- residuals.rma(mixed_inverse_weighting_dietC)

hist(residuals_mixed_inverse_weighting_dietC, breaks = 100)

summary(mixed_inverse_weighting_dietC)

flextable(tidy(summary(mixed_inverse_weighting_dietC)))


#### Sensitivity Analysis ####

## need meta-regressions from main code ####

mixed_ssdN <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                     mods = ~SSD,
                     random = 
                       list(~ 1 | Paper_Number,
                            ~ 1 | Species,
                            ~ 1 | Species.phylo),
                     R = list(Species.phylo = phylo_cor_meanN),
                     data=df_effects %>% drop_na(Mean_DiffN),
                     method="REML")

mixed_lnVR_ssdN <- rma.mv(lnVRN, lnVRN.sv, 
                          mods = ~SSD,
                          random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                        ~ 1 | Species.phylo),
                          R = list(Species.phylo = phylo_cor_varN),
                          data=df_effects %>% drop_na(lnVRN),
                          method="REML")

mixed_ssdC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                     mods = ~SSD,
                     random = 
                       list(~ 1 | Paper_Number, 
                            ~ 1 | Species,
                            ~ 1 | Species.phylo),
                     R = list(Species.phylo = phylo_cor_meanC),
                     data=df_effects %>% drop_na(Mean_DiffC),
                     method="REML")

mixed_lnVR_ssdC <- rma.mv(lnVRC, lnVRC.sv, 
                          mods = ~SSD,
                          random = list(~ 1 |                                             Paper_Number, 
                                        ~ 1 | Species,
                                        ~ 1 | Species.phylo),
                          R = list(Species.phylo = phylo_cor_varC),
                          data=df_effects %>% drop_na(lnVRC),
                          method="REML")

mixed_ssd_diet_size_1intA <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                    mods = ~SSD * Diet + Size,
                                    random = list(~ 1 |                                             Paper_Number,
                                                  ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

mixed_gape_phylo <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                           mods = ~SSD,
                           random = list(~ 1 |                                             Paper_Number,
                                         ~ 1 | Species,
                                         ~ 1 | Species.phylo),
                           R = list(Species.phylo = phylo_cor_gape),
                           data=df_effects %>% 
                             filter(Diet == "Carnivore") %>% 
                             filter(Gape_Lim == "Yes"),
                           method="REML")

mixed_carn_phylo <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                           mods = ~SSD,
                           random = list(~ 1 |                                             Paper_Number,
                                         ~ 1 | Species,
                                         ~ 1 | Species.phylo),
                           R = list(Species.phylo = phylo_cor_carn),
                           data=df_effects %>% 
                             filter(Diet == "Carnivore") %>%
                             filter(Gape_Lim == "No"),
                           method="REML")

mixed_ssd_diet_meanC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                               mods = ~SSD * Diet,
                               random = list(~ 1 |                                             Paper_Number,
                                             ~ 1 | Species,
                                             ~ 1 | Species.phylo),
                               R = list(Species.phylo = phylo_cor_meanC),
                               data=df_effects %>% drop_na(Mean_DiffC),
                               method="REML")

## also need row id as column in df_effects, to remove outliers

df_effects <- df_effects %>% 
  mutate(rowid_to_column(df_effects))

#### Cooks ####

## cook's mean diff N ####

cooks_mixed_ssdN <- cooks.distance(mixed_ssdN, progbar=TRUE,
                                  reestimate=TRUE)

plot(cooks_mixed_ssdN, type="o", pch=19, xlab ="d15N SSD Observed Outcome", ylab="Cook's Distance")


influential_mixed_sddN <- cooks_mixed_ssdN[(cooks_mixed_ssdN > (3 * mean(cooks_mixed_ssdN, na.rm = TRUE)))]

outliers_mixed_ssdN <- data.frame(rowid = as.numeric(names(influential_mixed_sddN)))

df_check_mixed_ssdN <- df_effects %>% drop_na(Mean_DiffN) %>%
  mutate(rowid_to_column(.), var = "rowid") %>% 
  anti_join(outliers_mixed_ssdN)

mixed_check_ssdN <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                           mods = ~SSD,
                           random = 
                             list(~ 1 | Paper_Number,
                                  ~ 1 | Species,
                                  ~ 1 | Species.phylo),
                           R = list(Species.phylo = phylo_cor_meanN),
                           data=df_check_mixed_ssdN,
                           method="REML")

residuals_mixed_check_ssdN <- residuals.rma(mixed_check_ssdN)

hist(residuals_mixed_check_ssdN, breaks = 100)

summary(mixed_check_ssdN)

flextable(tidy(mixed_check_ssdN))


## cooks N var ####

cooks_mixed_lnVR_ssdN <- cooks.distance(mixed_lnVR_ssdN, progbar=TRUE,
                                       reestimate=TRUE)

plot(cooks_mixed_lnVR_ssdN, type="o", pch=19, xlab ="d15N Var SSD Observed Outcome", ylab="Cook's Distance")


influential_mixed_lnVR_ssdN <- cooks_mixed_lnVR_ssdN[(cooks_mixed_lnVR_ssdN > (3 * mean(cooks_mixed_lnVR_ssdN, na.rm = TRUE)))]

outliers_mixed_lnVR_ssdN <- data.frame(rowid = as.numeric(names(influential_mixed_lnVR_ssdN)))

df_check_mixed_lnVR_ssdN <- df_effects%>% drop_na(lnVRN) %>%
  mutate(rowid_to_column(.), var = "rowid") %>% 
  anti_join(outliers_mixed_lnVR_ssdN)

mixed_check_lnVR_ssdN <- rma.mv(lnVRN, lnVRN.sv, 
                          mods = ~SSD,
                          random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                        ~ 1 | Species.phylo),
                          R = list(Species.phylo = phylo_cor_varN),
                          data=df_check_mixed_lnVR_ssdN,
                          method="REML")

residuals_mixed_check_lnVR_ssdN <- residuals.rma(mixed_check_lnVR_ssdN)

hist(residuals_mixed_check_lnVR_ssdN, breaks = 100)

summary(mixed_check_lnVR_ssdN)

flextable(tidy(mixed_check_lnVR_ssdN))

## cooks mean diff C ####

cooks_mixed_ssdC <- cooks.distance(mixed_ssdC, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_ssdC, type="o", pch=19, xlab ="d13C SSD Observed Outcome", ylab="Cook's Distance")


influential_mixed_sddC <- cooks_mixed_ssdN[(cooks_mixed_ssdC > (3 * mean(cooks_mixed_ssdC, na.rm = TRUE)))]

outliers_mixed_ssdC <- data.frame(rowid = as.numeric(names(influential_mixed_sddC)))

df_check_mixed_ssdC <- df_effects%>% drop_na(Mean_DiffC) %>%
  mutate(rowid_to_column(.), var = "rowid") %>% 
  anti_join(outliers_mixed_ssdC)

mixed_check_ssdC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                           mods = ~SSD,
                           random = 
                             list(~ 1 | Paper_Number,
                                  ~ 1 | Species,
                                  ~ 1 | Species.phylo),
                           R = list(Species.phylo = phylo_cor_meanC),
                           data=df_check_mixed_ssdC,
                           method="REML")

residuals_mixed_check_ssdC <- residuals.rma(mixed_check_ssdC)

hist(residuals_mixed_check_ssdC, breaks = 100)

summary(mixed_check_ssdC)

flextable(tidy(mixed_check_ssdC))


## cooks C var ####

cooks_mixed_lnVR_ssdC <- cooks.distance(mixed_lnVR_ssdC, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_lnVR_ssdC, type="o", pch=19, xlab ="d13C Var SSD Observed Outcome", ylab="Cook's Distance")


influential_mixed_lnVR_ssdC <- cooks_mixed_lnVR_ssdC[(cooks_mixed_lnVR_ssdC > (3 * mean(cooks_mixed_lnVR_ssdC, na.rm = TRUE)))]

outliers_mixed_lnVR_ssdC <- data.frame(rowid = as.numeric(names(influential_mixed_lnVR_ssdC)))

df_check_mixed_lnVR_ssdC <- df_effects%>% drop_na(lnVRC) %>%
  mutate(rowid_to_column(.), var = "rowid") %>% 
  anti_join(outliers_mixed_lnVR_ssdC)

mixed_check_lnVR_ssdC <- rma.mv(lnVRC, lnVRC.sv, 
                                mods = ~SSD,
                                random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                              ~ 1 | Species.phylo),
                                R = list(Species.phylo = phylo_cor_varC),
                                data=df_check_mixed_lnVR_ssdC,
                                method="REML")

residuals_mixed_check_lnVR_ssdC <- residuals.rma(mixed_check_lnVR_ssdC)

hist(residuals_mixed_check_lnVR_ssdC, breaks = 100)

summary(mixed_check_lnVR_ssdC)

flextable(tidy(mixed_check_lnVR_ssdC))

## nitrogen diet size ####

cooks_mixed_ssd_diet_size_1intA <- cooks.distance(mixed_ssd_diet_size_1intA, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_ssd_diet_size_1intA, type="o", pch=19, xlab ="d15N SSD * Diet + Size Observed Outcome", ylab="Cook's Distance")


influential_mixed_ssd_diet_size_1intA <- cooks_mixed_ssd_diet_size_1intA[(cooks_mixed_ssd_diet_size_1intA > (3 * mean(cooks_mixed_ssd_diet_size_1intA, na.rm = TRUE)))]

outliers_mixed_ssd_diet_size_1intA <- data.frame(rowid = as.numeric(names(influential_mixed_ssd_diet_size_1intA)))

df_check_mixed_ssd_diet_size_1intA <- df_effects %>% drop_na(Mean_DiffN) %>%
  mutate(rowid_to_column(.), var = "rowid") %>% 
  anti_join(outliers_mixed_ssd_diet_size_1intA)


mixed_check_ssd_diet_size_1intA <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
       mods = ~SSD * Diet + Size,
       random = list(~ 1 |                                             Paper_Number,
                     ~ 1 | Species,
                     ~ 1 | Species.phylo),
       R = list(Species.phylo = phylo_cor_meanN),
       data=df_check_mixed_ssd_diet_size_1intA,
       method="REML")

residuals_mixed_check_ssd_diet_size_1intA <- residuals.rma(mixed_check_ssd_diet_size_1intA)

hist(residuals_mixed_check_ssd_diet_size_1intA, breaks = 100)

summary(mixed_check_ssd_diet_size_1intA)

flextable(tidy(mixed_check_ssd_diet_size_1intA))

## nitrogen gape ####

cooks_mixed_gape_phylo <- cooks.distance(mixed_gape_phylo, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_gape_phylo, type="o", pch=19, xlab ="d15N SSD Gape-Limited Observed Outcome", ylab="Cook's Distance")


influential_mixed_gape_phylo <- cooks_mixed_gape_phylo[(cooks_mixed_gape_phylo > (3 * mean(cooks_mixed_gape_phylo, na.rm = TRUE)))]

outliers_mixed_gape_phylo <- data.frame(rowid = as.numeric(names(influential_mixed_gape_phylo)))

df_check_mixed_gape_phylo <- df_effects %>% 
  filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "Yes") %>% 
  mutate(rowid_to_column(., var = "rowid")) %>% 
  anti_join(outliers_mixed_gape_phylo)


mixed_check_gape_phylo <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                           mods = ~SSD,
                           random = list(~ 1 |                                             Paper_Number,
                                         ~ 1 | Species,
                                         ~ 1 | Species.phylo),
                           R = list(Species.phylo = phylo_cor_gape),
                           data=df_check_mixed_gape_phylo,
                           method="REML")

residuals_mixed_check_gape_phylo <- residuals.rma(mixed_check_gape_phylo)

hist(residuals_mixed_check_gape_phylo, breaks = 100)

summary(mixed_check_gape_phylo)

flextable(tidy(mixed_check_gape_phylo))

## nitrogen carnivore ####

cooks_mixed_carn_phylo <- cooks.distance(mixed_carn_phylo, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_carn_phylo, type="o", pch=19, xlab ="d15N SSD Non-Gape-Limited Observed Outcome", ylab="Cook's Distance")


influential_mixed_carn_phylo <- cooks_mixed_carn_phylo[(cooks_mixed_carn_phylo > (3 * mean(cooks_mixed_carn_phylo, na.rm = TRUE)))]

outliers_mixed_carn_phylo <- data.frame(rowid = as.numeric(names(influential_mixed_carn_phylo)))

df_check_mixed_carn_phylo <- df_effects %>%
  filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "No") %>% 
  mutate(rowid_to_column(., var = "rowid")) %>% 
  anti_join(outliers_mixed_carn_phylo)


mixed_check_carn_phylo <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                 mods = ~SSD,
                                 random = list(~ 1 |                                             Paper_Number,
                                               ~ 1 | Species,
                                               ~ 1 | Species.phylo),
                                 R = list(Species.phylo = phylo_cor_carn),
                                 data=df_check_mixed_carn_phylo,
                                 method="REML")

residuals_mixed_check_carn_phylo <- residuals.rma(mixed_check_carn_phylo)

hist(residuals_mixed_check_carn_phylo, breaks = 100)

summary(mixed_check_carn_phylo)

flextable(tidy(mixed_check_carn_phylo))


## carbon diet

cooks_mixed_ssd_diet_meanC <- cooks.distance(mixed_ssd_diet_meanC, progbar=TRUE,
                                   reestimate=TRUE)

plot(cooks_mixed_ssd_diet_meanC, type="o", pch=19, xlab ="d13C SSD * Diet Observed Outcome", ylab="Cook's Distance")


influential_cooks_mixed_ssd_diet_meanC <- cooks_mixed_ssd_diet_meanC[(cooks_mixed_ssd_diet_meanC > (3 * mean(cooks_mixed_ssd_diet_meanC, na.rm = TRUE)))]

outliers_mixed_ssd_diet_meanC <- data.frame(rowid = as.numeric(names(influential_cooks_mixed_ssd_diet_meanC)))

df_check_mixed_ssd_diet_meanC <- df_effects %>%
  drop_na(Mean_DiffC) %>% 
  mutate(rowid_to_column(., var = "rowid")) %>% 
  anti_join(outliers_mixed_ssd_diet_meanC)


mixed_check_ssd_diet_meanC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                                 mods = ~SSD * Diet,
                                 random = list(~ 1 |                                             Paper_Number,
                                               ~ 1 | Species,
                                               ~ 1 | Species.phylo),
                                 R = list(Species.phylo = phylo_cor_meanC),
                                 data=df_check_mixed_ssd_diet_meanC,
                                 method="REML")

residuals_mixed_check_ssd_diet_meanC <- residuals.rma(mixed_check_ssd_diet_meanC)

hist(residuals_mixed_check_ssd_diet_meanC, breaks = 100)

summary(mixed_check_ssd_diet_meanC)

flextable(tidy(mixed_check_ssd_diet_meanC))

