#### Isotopic Sex Differences Meta-Analysis - Data Analysis
#### Joshua Bauld

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

## read in data - df original = raw data, df_effects = csv created in "Isotopic_Sex_Differences_Dataframe_dryad.R"

df_original <- 

df_effects <- 

## read in phylogenies - different for each effect size, because of different numbers of species

load("phylogeny_matrix_meanN") ## nitrogen mean difference
load("phylogeny_matrix_meanC") ## carbon mean difference
load("phylogeny_matrix_varN")  ## nitrogen variation difference
load("phylogeny_matrix_varC")  ## carbon variation difference
load("phylogeny_matrix_carn")  ## non-gape-limited carnivores
load("phylogeny_matrix_gape")  ## gape-limited carnivores

## create Species.phylo variable to combine with phylogeny as a random factor in model - the values are the same as for species, but two random factors are needed, one to account for species differences like location and one to account for relatedness

df_effects <- df_effects %>%
  mutate(Species.phylo = Species)

## one row with n = 1 for both sexes creates problems later, as sampling variance is non-positive, so cannot calculate I2 etc. So, drop this row now.

df_effects <- df_effects %>% filter(Species != "Micronycteris schmidtorum")

## plot distributions

hist(df_effects$Mean_DiffN, breaks = 100)

hist(df_effects$Mean_DiffC, breaks = 100)

hist(df_effects$lnVRN, breaks = 100)

hist(df_effects$lnVRC, breaks = 100)


#### META-ANALYTIC MODELS ####

## used to calculate I2 for each effect size - how much inter-study variation is systematic and could be explained by moderators. Use I2 to annotate raw data plots later.

MA_Mean_DiffN <- rma(Mean_DiffN, Mean_DiffN.sv, data=df_effects)

summary(MA_Mean_DiffN)

MA_Mean_DiffC <- rma(Mean_DiffC, Mean_DiffC.sv, 
                     data=df_effects)

summary(MA_Mean_DiffC)


MA_lnVRN <- rma(lnVRN, lnVRN.sv, data=df_effects)

summary(MA_lnVRN)

MA_lnVRC <- rma(lnVRC, lnVRC.sv, data=df_effects)

summary(MA_lnVRC)

#### PLOTTING RAW DATA ####

## beeswarm plot to show heterogeneity in effect sizes

df_bee <- df_effects %>% select(Paper_Number, Mean_DiffN, lnVRN, Mean_DiffC, lnVRC)

str(df_bee)

df_bee <- df_bee %>% gather(Mean_DiffN, lnVRN, Mean_DiffC, lnVRC, key = "Effect_Size", value = "Value") %>% 
  drop_na()

lab_a <- paste("I^2 == 90.57", "*\\'%\\'")

a <- df_bee %>% subset(Effect_Size == "Mean_DiffN") %>% ggplot(aes(x = Effect_Size, y = Value)) +
  theme_classic() +
  geom_beeswarm(shape = 21, size = 4, stroke = 1.5, col = "black", fill = "royalblue4", alpha = 0.4, cex = 1.9) +
  annotate("text", x = 0.7, y = 1.8, size =8, label = lab_a, parse = T) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 25)) +
  labs(y = "Mean Sex Difference (\\u2030)")

lab_b <- paste("I^2 == 94.38", "*\\'%\\'")

b <- df_bee %>% subset(Effect_Size == "Mean_DiffC") %>% ggplot(aes(x = Effect_Size, y = Value)) +
  theme_classic() +
  geom_beeswarm(shape = 21, size = 4, stroke = 1.5, col = "black", fill = "coral2", alpha = 0.4, cex = 1.7) +
  annotate("text", x = 0.7, y = 1.8, size = 8, label = lab_b, parse = T)+ 
  scale_y_continuous(breaks= pretty_breaks()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 25))

lab_c <- paste("I^2 == 64.2", "*\\'%\\'")

c <- df_bee %>% subset(Effect_Size == "lnVRN") %>% ggplot(aes(x = Effect_Size, y = Value)) +
  theme_classic() +
  geom_beeswarm(shape = 23, size = 4, stroke = 1.5, col = "black", fill = "royalblue4", alpha = 0.4, cex = 1.9) +
  annotate("text", x = 0.6, y = -0.2, size = 8, label = lab_c, parse = T) +
  geom_segment(y = -0.6931472, yend = -0.6931472, x = 0.7, xend = 1.3, col = "black", linetype = "dotted", size = 1) +
  annotate("text", x = 0.8, y = -0.8, label = "\\u2640 2x Variable", size = 8) +
  geom_segment(y = 0.6931472, yend = 0.6931472, x = 0.7, xend = 1.3, col = "black", linetype = "dotted", size = 1)+
  annotate("text", x = 0.8, y = 0.85, label = "\\u2642 2x Variable", size =8) + theme(axis.title.x = element_text(size = 25, face = "bold"),
                                                                                     axis.text.x = element_blank(),
                                                                                     axis.ticks = element_blank(),
                                                                                     axis.title.y = element_text(size = 25, face = "bold"),
                                                                                     axis.text.y = element_text(size = 25)) +
  labs(x = "\\u03b415N",
       y = "log Variability Ratio")

lab_d <- paste("I^2 == 72.83", "*\\'%\\'")

d <- df_bee %>% subset(Effect_Size == "lnVRC") %>% ggplot(aes(x = Effect_Size, y = Value)) +
  theme_classic() +
  geom_beeswarm(shape = 23, size = 4, stroke = 1.5, col = "black", fill = "coral2", alpha = 0.4, cex = 1.9) +
  annotate("text", x = 0.6, y = -0.3, size = 8, label = lab_d, parse = T) +
  geom_segment(y = -0.6931472, yend = -0.6931472, x = 0.7, xend = 1.3, col = "black", linetype = "dotted", size = 1) +
  annotate("text", x = 0.8, y = -0.8, label = "\\u2640 2x Variable", size = 8) +
  geom_segment(y = 0.6931472, yend = 0.6931472, x = 0.7, xend = 1.3, col = "black", linetype = "dotted", size = 1)+
  annotate("text", x = 0.8, y = 0.85, label = "\\u2642 2x Variable", size =8) +
  theme(axis.title.x = element_text(size = 25, face = "bold"),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 25)) +
  labs(x = "\\u03b413C")


beeswarm <- ggarrange(a, b, c, d, labels = c("        A)", "        B)", "        C)", "        D)"), font.label = list(size = 25, face = "plain"), common.legend = T, legend = "bottom", nrow = 2, ncol = 2, align = "v")

ggsave("plot_beeswarm.jpeg", width = 400, height = 310, units = "mm", device = "jpeg")


#### PLOT MEAN-VARIANCE####

## using df_original for raw means and sd's for each sex - as we are interested in whether a mean variance relationship in the raw data should alter our effect sizes of choice.

## plots with geom_smooth lm indicate no strong mean variance relationship for any isotope, therefore raw mean differences and lnVR are appropriate choices of effect size.

plot_mean_var_fem_d15n <- df_original %>% ggplot(aes(x = M_F_d15N, y = SD_F_d15N)) +
  theme_classic() +
  geom_point(shape = 21, size = 10, col = "royalblue4", fill = "royalblue4", alpha = 0.4) +
  geom_smooth(method = "lm", col = "royalblue4", size = 2) +
  theme(axis.title.x = element_text(face = "bold", size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 25)) +
  labs(x = "Mean \\u03b415N (\\u2030)",
       y = "\\u03b415N sd (\\u2030)")

plot_mean_var_mal_d15n <- df_original %>% ggplot(aes(x = M_M_d15N, y = SD_M_d15N)) +
  theme_classic() +
  geom_point(shape = 21, size = 10, col = "royalblue4", fill = "royalblue4", alpha = 0.4) +
  geom_smooth(method = "lm", col = "royalblue4", size = 2) +
  theme(axis.title.x = element_text(face = "bold", size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 25)) +
  labs(x = "Mean \\u03b415N (\\u2030)",
       y = "\\u03b415N sd (\\u2030)")

plot_mean_var_fem_d13c <- df_original %>% ggplot(aes(x = M_F_d13C, y = SD_F_d13C)) +
  theme_classic() +
  geom_point(shape = 21, size = 10, col = "coral2", fill = "coral2", alpha = 0.4) +
  geom_smooth(method = "lm", col = "coral2", size = 2) +
  theme(axis.title.x = element_text(face = "bold", size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 25)) +
  labs(x = "Mean \\u03b413C (\\u2030)",
       y = "\\u03b413C sd (\\u2030)")

plot_mean_var_mal_d13c <- df_original %>% ggplot(aes(x = M_M_d13C, y = SD_M_d13C)) +
  theme_classic() +
  geom_point(shape = 21, size = 10, col = "coral2", fill = "coral2", alpha = 0.4) +
  geom_smooth(method = "lm", col = "coral2", size = 2) +
  theme(axis.title.x = element_text(face = "bold", size = 25),
        axis.text.x = element_text(size = 25),
        axis.title.y = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 25)) +
  labs(x = "Mean \\u03b413C (\\u2030)",
       y = "\\u03b413C sd (\\u2030)")

plot_mean_var <- ggarrange(plot_mean_var_fem_d15n, plot_mean_var_mal_d15n, plot_mean_var_fem_d13c, plot_mean_var_mal_d13c, labels = c("A) Female", "  B) Male", "C) Female", "  D) Male"), font.label = list(size = 25, face = "bold"), nrow = 2, ncol = 2, align = "hv")

ggsave("plot_mean_var.jpeg", width = 400, height = 300, units = "mm", device = "jpeg")


#### META-REGRESSIONS ####

## meta-regressions (linear mixed models) used to establish effect of size dimorphism on sex differences in mean and variance, of nitrogen and carbon stable isotope ratios

## nitrogen mean sex difference

mixed_ssdN <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                    mods = ~SSD,
                    random = 
                      list(~ 1 | Paper_Number,
                           ~ 1 | Species,
                           ~ 1 | Species.phylo),
                    R = list(Species.phylo = phylo_cor_meanN),
                    data=df_effects %>% drop_na(Mean_DiffN),
                    method="REML")

residuals_mixed_ssdN <- residuals.rma(mixed_ssdN)

hist(residuals_mixed_ssdN, breaks = 100)

summary(mixed_ssdN, digits = 10)

## nitrogen variation difference

mixed_lnVR_ssdN <- rma.mv(lnVRN, lnVRN.sv, 
                         mods = ~SSD,
                         random = list(~ 1 |                                             Paper_Number, ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_varN),
                         data=df_effects %>% drop_na(lnVRN),
                         method="REML")

residuals_mixed_lnVR_ssdN <- residuals.rma(mixed_lnVR_ssdN)

hist(residuals_mixed_lnVR_ssdN, breaks = 100)

summary(mixed_lnVR_ssdN)

mixed_ssdC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                     mods = ~SSD,
                     random = 
                       list(~ 1 | Paper_Number, 
                            ~ 1 | Species,
                            ~ 1 | Species.phylo),
                     R = list(Species.phylo = phylo_cor_meanC),
                     data=df_effects %>% drop_na(Mean_DiffC),
                     method="REML")

residuals_mixed_ssdC <- residuals.rma(mixed_ssdC)

hist(residuals_mixed_ssdC, breaks = 100)

summary(mixed_ssdC)


mixed_lnVR_ssdC <- rma.mv(lnVRC, lnVRC.sv, 
                          mods = ~SSD,
                          random = list(~ 1 |                                             Paper_Number, 
                                        ~ 1 | Species,
                                        ~ 1 | Species.phylo),
                          R = list(Species.phylo = phylo_cor_varC),
                          data=df_effects %>% drop_na(lnVRC),
                          method="REML")

residuals_mixed_lnVR_ssdC <- residuals.rma(mixed_lnVR_ssdC)

hist(residuals_mixed_lnVR_ssdC, breaks = 100)

summary(mixed_lnVR_ssdC)

#### Ecological Context ####

## adding additional moderators to examine whether effect of ssd on trophic sex differences is modified by ecological context - species mean size/ dietary class/ gape-limitation

## model selection with different combinations of ecological predictors

mixed_ssd_diet_noint <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                         mods = ~SSD + Diet,
                         random = list(~ 1 |                                             Paper_Number,
                                       ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_meanN),
                         data=df_effects %>% drop_na(Mean_DiffN),
                         method="REML")

summary(mixed_ssd_diet_noint, digits = 3)

mixed_ssd_diet <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                         mods = ~SSD * Diet,
                         random = list(~ 1 |                                             Paper_Number,
                                       ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_meanN),
                         data=df_effects %>% drop_na(Mean_DiffN),
                         method="REML")

summary(mixed_ssd_diet, digits = 3)

mixed_ssd_diet_size_noint <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                         mods = ~SSD + Diet + Size,
                         random = list(~ 1 |                                             Paper_Number,
                                       ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_meanN),
                         data=df_effects %>% drop_na(Mean_DiffN),
                         method="REML")

summary(mixed_ssd_diet_size_noint, digits = 3)

mixed_ssd_diet_size_1intA <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                    mods = ~SSD * Diet + Size,
                                    random = list(~ 1 |                                             Paper_Number,
                                                  ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

summary(mixed_ssd_diet_size_1intA, digits = 10)

mixed_ssd_diet_size_1intB <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                    mods = ~SSD + Diet * Size,
                                    random = list(~ 1 |                                             Paper_Number,
                                                  ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

summary(mixed_ssd_diet_size_1intB, digits = 3)

mixed_ssd_diet_size_1intC <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                    mods = ~SSD * Size + Diet,
                                    random = list(~ 1 |                                             Paper_Number,
                                                  ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

summary(mixed_ssd_diet_size_1intC, digits = 3)


mixed_ssd_diet_size_int <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                    mods = ~SSD * Diet * Size,
                                    random = list(~ 1 |                                             Paper_Number,
                                                  ~ 1 | Species,
                                                  ~ 1 | Species.phylo),
                                    R = list(Species.phylo = phylo_cor_meanN),
                                    data=df_effects %>% drop_na(Mean_DiffN),
                                    method="REML")

summary(mixed_ssd_diet_size_int, digits = 3)

mixed_ssd_size_noint <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                                  mods = ~SSD + Size,
                                  random = list(~ 1 |                                             Paper_Number,
                                                ~ 1 | Species,
                                                ~ 1 | Species.phylo),
                                  R = list(Species.phylo = phylo_cor_meanN),
                                  data=df_effects %>% drop_na(Mean_DiffN),
                                  method="REML")

summary(mixed_ssd_size_noint, digits = 3)

mixed_ssd_size_int <- rma.mv(Mean_DiffN, Mean_DiffN.sv, 
                               mods = ~SSD * Size,
                               random = list(~ 1 |                                             Paper_Number,
                                             ~ 1 | Species,
                                             ~ 1 | Species.phylo),
                               R = list(Species.phylo = phylo_cor_meanN),
                               data=df_effects %>% drop_na(Mean_DiffN),
                               method="REML")

summary(mixed_ssd_size_int, digits = 3)




## create table comparing aic of different models

options(pillar.sigfig=6)

eco_con_aic <- as.data.frame(AIC(mixed_ssd_diet_noint, mixed_ssd_diet, mixed_ssd_diet_size_noint, mixed_ssd_diet_size_1intA, mixed_ssd_diet_size_1intB, mixed_ssd_diet_size_1intC, mixed_ssd_diet_size_int, mixed_ssd_size_noint, mixed_ssd_size_int)) %>% 
  rownames_to_column() %>% 
  as_tibble() %>% 
  arrange(AIC)

flextable(eco_con_aic)

options(pillar.sigfig=3)


#### Gape_Limitation #### 

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

summary(mixed_gape_phylo, digits = 3)

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

summary(mixed_carn_phylo, digits = 3)


#### Carbon and Diet ####

mixed_ssd_diet_meanC <- rma.mv(Mean_DiffC, Mean_DiffC.sv, 
                         mods = ~SSD * Diet,
                         random = list(~ 1 |                                             Paper_Number,
                                       ~ 1 | Species,
                                       ~ 1 | Species.phylo),
                         R = list(Species.phylo = phylo_cor_meanC),
                         data=df_effects %>% drop_na(Mean_DiffC),
                         method="REML")

summary(mixed_ssd_diet_meanC, digits = 3)



#### PLOT META-REGRESSIONS ####

#### SSD vs Effect Sizes

## Nitrogen Mean Differences

newdata_ssd <- tibble(SSD = seq(-9, 6, length = 1000))

pred_ssdN <- as_tibble(predict.rma(mixed_ssdN, newmods = newdata_ssd$SSD))

df_pred_ssdN <- bind_cols(pred_ssdN, newdata_ssd)

scatter_ssdN <- df_pred_ssdN %>% ggplot(aes(x = SSD, y = pred)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects, aes(x = SSD, y = Mean_DiffN), col= "black", size = 4, stroke = 2, alpha = 0.5, fill= "cyan3", shape = 21) +
  geom_line(col = "cyan3", size = 1.5) +
  geom_ribbon(aes(ymax = ci.ub, ymin = ci.lb),
              col = "cyan2", fill = "cyan3", alpha = 0.4) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.text.y = element_text(size = 32)) +
  labs(y = "Mean Difference (\\u2030)")

## Carbon Mean Difference

pred_ssdC <- as_tibble(predict.rma(mixed_ssdC, newmods = newdata_ssd$SSD))

df_pred_ssdC <- bind_cols(pred_ssdC, newdata_ssd)

scatter_ssdC <- df_pred_ssdC %>% ggplot(aes(x = SSD, y = pred)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects, aes(x = SSD, y = Mean_DiffC), col = "black", size = 4, stroke = 2, alpha = 0.5, fill= "coral", shape = 21) +
  geom_line(col = "coral", size = 1.5) +
  geom_ribbon(aes(ymax = ci.ub, ymin = ci.lb),
              col = "coral", fill = "coral", alpha = 0.4) +
  scale_y_continuous(breaks= pretty_breaks()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32))

## Nitrogen Variation Difference

pred_lnVR_ssdN <- as_tibble(predict.rma(mixed_lnVR_ssdN, newmods = newdata_ssd$SSD))

df_pred_lnVR_ssdN <- bind_cols(pred_lnVR_ssdN, newdata_ssd)

scatter_lnVR_ssdN <- df_pred_lnVR_ssdN %>% ggplot(aes(x = SSD, y = pred)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects, aes(x = SSD, y = lnVRN), col = "black", fill = "cyan3", size = 4, stroke = 2, shape = 23, alpha = 0.5) +
  geom_line(col = "cyan3", size = 1.5) +
  geom_ribbon(aes(ymax = ci.ub, ymin = ci.lb),
              col = "cyan3", fill = "cyan3", alpha = 0.4) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title.y = element_text(size = 32, face = "bold"),
        axis.text.y = element_text(size = 32)) +
  labs(y = "log Variability Ratio")

## carbon variation differences

pred_lnVR_ssdC <- as_tibble(predict.rma(mixed_lnVR_ssdC, newmods = newdata_ssd$SSD))

df_pred_lnVR_ssdC <- bind_cols(pred_lnVR_ssdC, newdata_ssd)

scatter_lnVR_ssdC <- df_pred_lnVR_ssdC %>% ggplot(aes(x = SSD, y = pred)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects, aes(x = SSD, y = lnVRC), col = "black", fill = "coral", size = 4, stroke = 2, shape = 23, alpha = 0.5) +
  geom_line(col = "coral", size = 1.5) +
  geom_ribbon(aes(ymax = ci.ub, ymin = ci.lb),
              col = "coral", fill = "coral", alpha = 0.4) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 32))

scatters_effect_sizes <- ggarrange(scatter_ssdN, scatter_ssdC, scatter_lnVR_ssdN, scatter_lnVR_ssdC, labels = c("    A) \\u03b415N", "    B) \\u03b413C", "    C) \\u03b415N", "    D) \\u03b413C"), font.label = list(size = 32, face = "bold"), common.legend = T, legend = "bottom", nrow = 2, ncol = 2, align = c("hv"))

scatters_effect_sizes <- annotate_figure(scatters_effect_sizes, bottom = text_grob("Sexual Size Dimorphism", face = "bold", size = 32))

ggsave("plot_scatter_effect_sizes.jpeg", width = 500, height = 300, units = "mm", device = "jpeg")


#### Ecological Context ####

## predictions for diet categories

eco_con_aic

## the best model, according to AIC, included diet and size. As the effect of size was very small, only plotting predictions using for different diet classes, based on this model: 

## Best model: mixed_ssd_diet_size_1intA: ~SSD * Diet + Size


## Diet

newdata_diet <- tibble(SSD = c(seq(-9, 6, length = 1000),
                               seq(-4, 4, length = 1000),
                               seq(-4, 4, length = 1000)),
                       Size = c(rep(0.00192, 1000),
                                rep(0.0496, 1000),
                                rep(0.00638, 1000)),
                       Herbivore = c(seq(0, 0, length = 1000),
                                     seq(1, 1, length = 1000),
                                     seq(0, 0, length = 1000)),
                       Omnivore = c(seq(0, 0, length = 1000),
                                    seq(0, 0, length = 1000),
                                    seq(1, 1, length = 1000))) %>% 
  mutate(SSD_Herbivore = SSD * Herbivore,
         SSD_Omnivore = SSD * Omnivore)

pred_diet <- as_tibble(predict(mixed_ssd_diet_size_1intA, 
                     newmods = cbind(newdata_diet$SSD,
                                     newdata_diet$Size,
                                newdata_diet$Herbivore,
                                newdata_diet$Omnivore,
                                newdata_diet$SSD_Herbivore,
                                newdata_diet$SSD_Omnivore)))

df_pred_diet <- bind_cols(newdata_diet, pred_diet) %>% 
  mutate(Diet = 
           if_else(Herbivore == 0 & Omnivore == 0, "Carnivore",
                   if_else(Herbivore == 1 & Omnivore == 0, "Herbivore", "Omnivore")))

scatter_carn <- df_pred_diet %>% filter(Diet == "Carnivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Carnivore"),
             aes(x = SSD, y = Mean_DiffN),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-5, 0, 5)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(label = "Carnivore", values = "firebrick4") +
  scale_fill_manual(label = "Carnivore", values = "firebrick4") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")

scatter_omni <- df_pred_diet %>% filter(Diet == "Omnivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Omnivore"),
             aes(x = SSD, y = Mean_DiffN),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(label = "Omnivore", values = "darkorchid4") +
  scale_fill_manual(label = "Omnivore", values = "darkorchid4") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")

scatter_herb <- df_pred_diet %>% filter(Diet == "Herbivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Herbivore"),
             aes(x = SSD, y = Mean_DiffN),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  scale_color_manual(label = "Herbivore", values = "darkgreen") +
  scale_fill_manual(label = "Herbivore", values = "darkgreen") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")


plot_ssd_dietN <- ggarrange(scatter_carn, ggarrange(scatter_omni, scatter_herb, ncol = 2, labels = c("B) Omnivores", "C) Herbivores"), font.label = list(size = 30, color = "black", face = "bold", family = NULL), align = "hv"), nrow = 2, labels = "A) Carnivores", font.label = list(size = 30, color = "black", face = "bold", family = NULL), align = "h")

plot_ssd_dietN <- annotate_figure(plot_ssd_dietN, bottom = text_grob("Sexual Dimorphism", face = "bold", size = 30), left = text_grob("Mean \\u03b415N Sex Difference (\\u2030)", rot = 90, size = 30, face = "bold"))

ggsave("plot_ssd_dietN.jpeg", width = 500, height = 250, units = "mm", device = "jpeg")


# ## Size
# 
# ## min and max = 0.00192 & 170000
# 
# ## predict effect of ssd at these two sizes, to illustrate interaction - though effect is so small the slopes will be similar
# 
# newdata_size <- tibble(SSD = c(seq(-9, 6, length = 1000),
#                                seq(-9, 6, length = 1000)),
#                        Size = c(seq(0.665, 0.665, length = 1000),
#                                 seq(110, 110, length = 1000))) %>% 
#   mutate(SSD_Size = SSD*Size)
# 
# pred_size <- as_tibble(predict(mixed_ssd_size_int, newmods = cbind(newdata_size$SSD, newdata_size$Size, newdata_size$SSD_Size)))
# 
# df_pred_size <- bind_cols(newdata_size, pred_size) %>% 
#   mutate(Size_Cat = 
#            if_else(Size == "0.665", "Smallest", "Largest"))
# 
# ## add size category to df_effects for plotting
# 
# df_effects <- df_effects %>% 
#   mutate(Size_Cat = 
#            if_else(Size <3.50, "Smallest", "Largest"))
# 
# plot_scatter_mixed_size_meanN <- df_pred_size %>% 
#   ggplot(aes(x = SSD, y = pred, col = Size_Cat, fill = Size_Cat))+
#   theme_classic() +
#   geom_hline(yintercept = 0, linetype = "dotted") +
#   geom_vline(xintercept = 0, linetype = "dotted") +
#    geom_point(data = df_effects,
#               aes(x = SSD, y = Mean_DiffN),
#               col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
#   geom_ribbon(aes(ymin = ci.lb,
#                   ymax = ci.ub),
#               alpha = 0.4) +
#   geom_line(size = 1.5) +
#   scale_color_manual(labels = c("Largest", "Smallest"),
#                      breaks = c("Largest", "Smallest"),
#                      values = c("slateblue4", "orchid4")) +
#   scale_fill_manual(labels = c("Largest", "Smallest"),
#                     breaks = c("Largest", "Smallest"),
#                     values = c("slateblue4", "orchid4")) +
#   facet_wrap(~Size_Cat, nrow = 2)+
#   theme(axis.title = element_text(size = 25, face = "bold"),
#         axis.text = element_text(size = 20),
#         strip.background = element_blank(),
#         strip.text.x = element_text(size = 20, face = "bold", hjust = 0),
#         legend.position = "none") +
#   labs(x = "Sexual Dimorphism",
#        y = "Mean \\u03b415N Sex Difference (\\u2030)")
# 
# ggsave("plot_scatter_mixed_size_meanN.jpeg", width = 300, height = 200, units = "mm", device = "jpeg")


#### Gape Limitation ####

## make two sets of predictions: one for non-gape-limited carnivores and one for gape_limited carnivores, both over the SSD range of non-gape-limited carnivores

## then subset the gape-limited predictions into two dataframes - one within the SSD range of gape-limited carnivores and one outside of this range

df_effects %>% filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "Yes") %>% 
  summarize(range = range(SSD))

df_effects %>% filter(Diet == "Carnivore") %>% 
  filter(Gape_Lim == "No") %>% 
  summarize(range = range(SSD))

newdata_gape <- tibble(SSD = seq(-9, 6, length = 1000)) 

## NGL carnivores

pred_carn <- as_tibble(predict(mixed_carn_phylo, newmods = cbind(newdata_gape$SSD)))

df_pred_carn <- bind_cols(newdata_gape, pred_carn)

## GL Carnivores phylogenetic control

pred_gape <- as_tibble(predict(mixed_gape_phylo, newmods = cbind(newdata_gape$SSD)))

df_pred_gape <- bind_cols(newdata_gape, pred_gape)

## for color and fill aesthetics/ legend, need to add group variables to df_effects and the pred df's here

## for df effects, need the snake species to be labelled "snakes", not reptiles, and then use subphylum for other names, which will leave fish labelled as fish

df_effects <- df_effects %>% 
  mutate(Group = 
           if_else(Species == "Laticauda saintgironsi", "Snakes",
                   if_else(Species == "Nerodia rhombifer", "Snakes",
                           if_else(Species == "Nerodia erythrogaster", "Snakes",
                                   if_else(Species == "Nerodia sipedon", "Snakes", Subphylum)))))

## for the prediction df's, need to label group as either Non-Gape-Limited Carnivores or Gape-Limited Carnivores

df_pred_carn <- df_pred_carn %>% 
  mutate(Group = "Non-Gape-Limited\\nCarnivores")

df_pred_gape <- df_pred_gape %>% 
  mutate(Group = "Gape-Limited\\nCarnivores")


plot_scatter_gape <- df_pred_gape %>% filter(SSD >= -4.5 & SSD <= 4.5) %>% 
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


## read in df of phylogenetic signal for gape-limited predators

df_lipa_gape <- read_csv("df_lipa_gape.csv")


plot_lipa_gape <- df_lipa_gape %>% 
  ggplot(aes(x = fct_reorder(species, desc(lipa_score)), y = lipa_score)) +
  geom_hline(yintercept = 0, size = 1) +
  geom_hline(yintercept = -1, linetype = "dotted", size = 1) +
  geom_hline(yintercept = -0.5, linetype = "dotted", size = 1) +
  geom_hline(yintercept = 0.5, linetype = "dotted", size = 1) +
  theme_classic() +
  geom_point(aes(col = group, fill = group), 
             size = 5, stroke = 2, shape = 21, alpha = 0.8) +
  scale_color_manual(values = c("black", "black"),
                     labels = c("Fish", "Snakes")) +
  scale_fill_manual(values = c("#d8b365", "#1b7837"),
                    labels = c("Fish", "Snakes")) +
  coord_flip() +
  theme(axis.title = element_text(size = 25, face = "bold"),
        axis.text.y = element_text(size = 15, face = "italic"),
        axis.text.x = element_text(size = 15, angle = 55, 
                                   vjust = 0.58),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "none") +
  guides(legend = "none") +
  labs(y = "Phylogenetic Signal",
       x = "Species",
       fill = "Group") 

ggsave("plot_lipa_gape.jpeg", units = "mm", height = 200, width = 250)


plot_scatter_gape_phylosignal <- ggarrange(plot_scatter_gape, plot_lipa_gape, ncol = 2,
           heights = c(1,1), widths = c(3,2),
          labels = c("A", "B"), hjust = c(-6, 0),
          font.label = list(size = 24, color = "black", face = "bold"))

ggsave("plot_scatter_gape_phylosignal.jpeg", units = "mm",
       width = 400, height = 150)


#### Carbon ~ SSD * Diet ####

newdata_dietC <- tibble(SSD = c(seq(-9, 6, length = 1000),
                               seq(-4, 4, length = 1000),
                               seq(-4, 4, length = 1000)),
                       Herbivore = c(seq(0, 0, length = 1000),
                                     seq(1, 1, length = 1000),
                                     seq(0, 0, length = 1000)),
                       Omnivore = c(seq(0, 0, length = 1000),
                                    seq(0, 0, length = 1000),
                                    seq(1, 1, length = 1000))) %>% 
  mutate(SSD_Herbivore = SSD * Herbivore,
         SSD_Omnivore = SSD * Omnivore)

pred_dietC <- as_tibble(predict(mixed_ssd_diet_meanC, 
                               newmods = cbind(newdata_dietC$SSD,
                                               newdata_dietC$Herbivore,
                                               newdata_dietC$Omnivore,
                                               newdata_dietC$SSD_Herbivore,
                                               newdata_dietC$SSD_Omnivore)))

df_pred_dietC <- bind_cols(newdata_dietC, pred_dietC) %>% 
  mutate(Diet = 
           if_else(Herbivore == 0 & Omnivore == 0, "Carnivore",
                   if_else(Herbivore == 1 & Omnivore == 0, "Herbivore", "Omnivore")))

scatter_carnC <- df_pred_dietC %>% filter(Diet == "Carnivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Carnivore"),
             aes(x = SSD, y = Mean_DiffC),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-5, 0, 5)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(label = "Carnivore", values = "firebrick4") +
  scale_fill_manual(label = "Carnivore", values = "firebrick4") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")

scatter_omniC <- df_pred_dietC %>% filter(Diet == "Omnivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Omnivore"),
             aes(x = SSD, y = Mean_DiffC),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_color_manual(label = "Omnivore", values = "darkorchid4") +
  scale_fill_manual(label = "Omnivore", values = "darkorchid4") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")

scatter_herbC <- df_pred_dietC %>% filter(Diet == "Herbivore") %>% 
  ggplot(aes(x = SSD, y = pred, col = Diet, fill = Diet)) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(data = df_effects%>% filter(Diet == "Herbivore"),
             aes(x = SSD, y = Mean_DiffC),
             col = "black", shape = 21, size = 4, stroke = 2, alpha = 0.5) +
  geom_ribbon(aes(ymin = ci.lb,
                  ymax = ci.ub),
              alpha = 0.4) +
  geom_line(size = 1.5) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  scale_y_continuous(breaks = c(-4, 0, 4)) +
  scale_color_manual(label = "Herbivore", values = "darkgreen") +
  scale_fill_manual(label = "Herbivore", values = "darkgreen") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 32),
        legend.position = "none")


plot_ssd_dietC <- ggarrange(scatter_carnC, ggarrange(scatter_omniC, scatter_herbC, ncol = 2, labels = c("B) Omnivores", "C) Herbivores"), font.label = list(size = 30, color = "black", face = "bold", family = NULL), align = "hv"), nrow = 2, labels = "A) Carnivores", font.label = list(size = 30, color = "black", face = "bold", family = NULL), align = "h")

plot_ssd_dietC <- annotate_figure(plot_ssd_dietC, bottom = text_grob("Sexual Dimorphism", face = "bold", size = 30), left = text_grob("Mean \\u03b413C Sex Difference (\\u2030)", rot = 90, size = 30, face = "bold"))

ggsave("plot_ssd_dietC.jpeg", width = 500, height = 250, units = "mm", device = "jpeg")

