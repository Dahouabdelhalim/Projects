########## R Code for Data Analysis ##########


#### Manuscript Information
## Authors: ND Muchoney, MD Bowers, AL Carper, MB Teglas, AM Smilanich
## Title: Use of an exotic Host plant reduces viral burden in a native insect herbivore
## Journal: Ecology Letters


#### File Description
## This file contains R code used to perform all logistic regression and linear regression models
## presented in this manuscript. All statistical analyses were performed in R version 4.0.4 using
## the "car" package for two-way ANOVA and the "stats" package for all additional analyses (e.g.,
## logistic regression and linear regression models). This code file contains procedures for model 
## selection (where applicable), statistical output, and data visualization, along with pertinent 
## data transformations and calculations. All analyses utilize data provided in the accompanying 
## file “Muchoney_et_al_LabExperiment.csv.” All packages required for this code are listed below.
## Code Author: ND Muchoney


#### Load Packages
library(car)
library(dplyr)
library(emmeans)
library(ggplot2)
library(ggpubr)
library(ggpattern)
library(interactions)
library(MuMIn)
library(plyr)
library(sjPlot)


#### Load Data
aj_data <- read.csv("Muchoney_et_al_LabExperiment.csv", header = T)


########## Section 1: Survival ##########

#### Section Description
## Summary information and logistic regression model evaluating variation in survival probability of 
## Anartia jatrophae individuals based on host plant species (Bacopa monnieri or Plantago lanceolata),
## viral treatment (orally challenged with Junonia coenia densovirus or unchallenged controls), and
## immune assay status (underwent immune assays during the larval stage or did not undergo assays).

#### Frequency Tables: Survival
## Survival Based on Host Plant/Viral Treatment
aj_data_surv <- aj_data[which(aj_data$surv_exclude != "Yes"),]
survival_table <- ddply(aj_data_surv, c("host_plant","treatment"), summarise, 
                   N = sum(!is.na(survived)),
                   surv_freq = (sum(survived == "Yes")/N)*100)
survival_table

#### Logistic Regression: Effects of Host Plant, Viral Treatment, and Immune Assessment on Survival
## Create Numeric Binary Variable (0/1) for Survival
aj_data_surv <- aj_data_surv %>%
    mutate(survived_num = case_when(
    survived == "No" ~ 0,
    survived == "Yes" ~ 1
    ))

## Model Selection
# Global Model
survival_glm_global <- glm(formula = survived_num ~ host_plant*treatment*assayed, family = 
                             binomial(link=logit), data = aj_data_surv, na.action = "na.fail")
survival_glm_global 
# Model Comparison
survival_models <- dredge(survival_glm_global, fixed =c("host_plant","treatment","assayed"))
survival_models
# Single model (survived ~ host_plant*treatment + assayed) selected based on AICc/nesting rule

## Final Model and Table S2
survival_glm_final <- glm(formula = survived_num ~ host_plant*treatment+assayed, family = 
                            binomial(link=logit), data = aj_data_surv)
summary(survival_glm_final)
tab_model(survival_glm_final, show.est = T, show.ci = 0.95, show.se = F, show.stat = T, 
          show.r2 = T,show.intercept = T, dv.labels = c("Survival to adult"), pred.labels = 
          c("(Intercept)","Host plant","Treatment","Assayed","Host plant x treatment"), 
          string.ci = "95% CI", string.est = "Odds Ratio", string.stat = "z")

## Pairwise Contrasts
emmeans(survival_glm_final, pairwise ~ host_plant*treatment, adjust = "tukey", type = "response")

## Create Figure 2
figure2 <- ggline(survival_table, x = "treatment", y = "surv_freq", group = "host_plant", 
                  plot_type = "l", size = .5, color = "black", linetype = "host_plant", 
                  xlab = "Treatment", ylab = "Survival to adult stage (%)", 
                  legend = c(0.82,0.16), legend.title = "Host plant", ylim = c(20,100)) + 
                  font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) + 
                  font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) + 
                  scale_x_discrete(labels = c("Control", "Virus-challenged")) + 
                  scale_linetype_manual(values=c("solid", "longdash"), 
                                        labels = c(expression(italic("Bacopa")), 
                                                   expression(italic("Plantago")))) + 
                  geom_point(shape = 21, size = 4, color = "black", 
                              aes(fill = host_plant, group = host_plant)) + 
                  scale_fill_manual(values = c("lightsteelblue2","brown4"), guide = "none")
figure2

## Export Figure 2
tiff("Figure2.tif", units = "in", width = 4.2, height = 4, res = 600)
figure2
dev.off()


########## Section 2: Viral Detection ##########

#### Section Description
# Summary information and statistical analyses evaluating postmortem detection of Junonia coenia
# densovirus within the subset of Anartia jatrophae individuals that were inoculated with the virus 
# during the larval stage. The frequency of viral detection (presence/absence of JcDV in insect 
# tissues following death) is compared based on larval host plant species (Bacopa monnieri or 
# Plantago lanceolata) using a Chi-squared test. The relationship between viral detection and 
# survival probability is assessed using a logistic regression model, which includes host plant,
# viral detection (Y/N), and immune assay status as predictors and survival (Y/N) as the response.

#### Frequency Tables: Viral Detection
## Viral Detection Based on Host Plant
aj_data_virsurv <- aj_data[which(aj_data$treatment == "Virus" & 
                                    aj_data$surv_exclude != "Yes"),]
detect_table <- ddply(aj_data_virsurv, c("host_plant"), summarise, 
                    N = sum(!is.na(jcdv_present)),
                    detected = sum(jcdv_present == "Yes", na.rm =T),
                    not_detected = sum(jcdv_present == "No", na.rm=T),
                    detect_freq = (sum(jcdv_present == "Yes", na.rm =T)/N)*100)
detect_table

## Survival Based on Viral Detection
survdetect_table <- ddply(aj_data_virsurv, c("host_plant","jcdv_present"), summarise,
                     N = sum(!is.na(survived)),
                     surv_freq = (sum(survived == "Yes")/N)*100)
survdetect_table <- survdetect_table[-5,]
survdetect_table

#### Chi-squared Test: Relationship between Host Plant Use and Viral Detection
## Chi-squared Test
detection_chisq <- chisq.test(table(aj_data_virsurv$jcdv_present, aj_data_virsurv$host_plant))
detection_chisq

## Create Figure 3A
figure3a <- ggbarplot(survdetect_table, "host_plant", "N",color = "host_plant", fill = "host_plant") + 
                        scale_x_discrete(labels = c(expression(italic("Bacopa")), 
                                                    expression(italic("Plantago")))) +
                        xlab(c("Host plant")) + ylab(c("Postmortem individuals (n)")) +
                        font("xlab", size = 14) + font("ylab", size = 14) + 
                        geom_col_pattern(aes(fill = host_plant, pattern_spacing = jcdv_present, 
                                             pattern_density = jcdv_present, pattern_colour = host_plant, 
                                             color= host_plant), width = .75, pattern = "stripe", 
                                         pattern_angle = 45, pattern_fill = "white") +
                        scale_pattern_spacing_manual(values = c(.1,5), guide = "none") +
                        scale_pattern_density_manual(values = c(1,.5), guide = "none") +
                        scale_pattern_color_manual(values = c("lightsteelblue3","brown4"), guide = "none") +
                        scale_color_manual(values = c("lightsteelblue3","brown4"), guide = "none") +
                        scale_fill_manual(values = c("lightsteelblue3","brown4"), guide = "none")
figure3a

## Export Figure 3A
tiff("Figure3a_fin.tif", units = "in", width = 4.2, height = 4, res = 600)
figure3a
dev.off()

#### Logistic Regression: Effects of Host Plant, Viral Detection, and Immune Assessment on Survival
## Create Numeric Binary Variable (0/1) for Survival
aj_data_virsurv<- aj_data_virsurv %>%
  mutate(survived_num = case_when(
    survived == "No" ~ 0,
    survived == "Yes" ~ 1
  ))

## Model Selection
# Global Model
aj_data_virsurvdet <- aj_data_virsurv[which(aj_data_virsurv$jcdv_present != "NA"),]
survdetect_glm_global <- glm(formula = survived_num ~ host_plant*jcdv_present*assayed, 
                             family = binomial(link=logit), data = aj_data_virsurvdet, na.action = "na.fail")
survdetect_glm_global
# Model Comparison
survdetect_models <- dredge(survdetect_glm_global, fixed =c("host_plant","jcdv_present","assayed"))
survdetect_models 
# Single model (survived ~ host_plant + jcdv_present + assayed) selected based on AICc/nesting rule

## Final Model and Table
survdetect_glm_final <- glm(formula = survived_num ~ host_plant+jcdv_present+assayed, 
                            family = binomial(link=logit), data = aj_data_virsurv)
summary(survdetect_glm_final)
tab_model(survdetect_glm_final, show.est = T, show.ci = 0.95, show.se = F, show.stat = T, 
          show.r2 = T,show.intercept = T, dv.labels = c("Survival to adult"), 
          pred.labels = c("(Intercept)","Host plant","Virus detected","Assayed"), 
          string.ci = "95% CI", string.est = "Odds Ratio", string.stat = "z")

## Pairwise Contrasts
emmeans(survdetect_glm_final, pairwise ~ host_plant+jcdv_present, adjust = "tukey", type = "response")

## Create Figure 3B
figure3b <- ggline(survdetect_table, x = "jcdv_present", y = "surv_freq", group = "host_plant", 
                     plot_type = "l", size = .5, color = "black", linetype = "host_plant",  
                     xlab = "Postmortem viral detection", ylab = "Survival to adult stage (%)",
                     legend = c(0.24,0.2), legend.title = "Host plant", ylim = c(20,100)) + 
                      font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                      font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) + 
                      scale_x_discrete(labels = c("Not detected","Detected")) + 
                      scale_linetype_manual(values=c("solid", "longdash"), 
                                            labels = c(expression(italic("Bacopa")), 
                                                       expression(italic("Plantago")))) + 
                      geom_point(shape = 21, size = 4, color = "black", 
                                 aes(fill = host_plant, group = host_plant)) + 
                      scale_fill_manual(values = c("lightsteelblue2","brown4"), guide = "none")
figure3b

## Export Figure 3B
tiff("Figure3b.tif", units = "in", width = 4.2, height = 4, res = 600)
figure3b
dev.off()


########## Section 3: Viral Load ##########

#### Section Description
# Summary information and statistical analyses evaluating viral infection loads of Anartia jatrophae 
# individuals that were orally inoculated with Junonia coenia densovirus during the larval stage and 
# maintained a detectable infection at their time of death. Postmortem JcDV loads are compared using 
# a multiple regression model with host plant (Bacopa monnieri or Plantago lanceolata), survival to 
# the adult stage (Y/N), and immune assay status (assayed/not assayed) as predictors. The relationship 
# between viral load and survival is further probed using a logistic regression model with viral load, 
# host plant, their interaction, and assay status as predictors and survival (Y/N) as the response. 

#### Calculate Relative Viral Load 
## Formula: relative load = 2^-delta(Ct)
aj_data$ct_vp4_mean <- rowMeans(aj_data[,c("ct_vp4_1","ct_vp4_2")], na.rm = T)
aj_data$ct_28s_mean <- rowMeans(aj_data[,c("ct_28s_1","ct_28s_2")], na.rm = T)
aj_data$rel_load <- 2^-(aj_data$ct_vp4_mean - aj_data$ct_28s_mean)
## Log-Transformation
aj_data$log_load <- log10(aj_data$rel_load)
aj_data$log_load[is.nan(aj_data$log_load)] <- NA

#### Summary Tables: Viral Load
## Viral Load Based on Host Plant
aj_data_virload <- aj_data[which(aj_data$treatment == "Virus" &
                                   aj_data$surv_exclude != "Yes" & 
                                   aj_data$log_load != "NA"),]
virload_table1 <- ddply(aj_data_virload, c("host_plant"), summarise, 
                   N = sum(!is.na(log_load)),
                   survived = "Overall",
                   mean_load = mean(log_load, na.rm=TRUE),
                   sd   = sd(log_load, na.rm=TRUE),
                   se   = sd / sqrt(N))
virload_table1
# Geometric mean of viral load on Bacopa: 10^-4.909394 = 1.231987e-05
# Geometric mean of viral load on Plantago: 10^-7.296795 = 5.048996e-08
# Fold change: 1.231987e-05/5.048996e-08 =  244.0063

## Viral Load Based on Host Plant and Survival
virload_table2 <- ddply(aj_data_virload, c("host_plant","survived"), summarise, 
                     N = sum(!is.na(log_load)),
                     mean_load = mean(log_load, na.rm=TRUE),
                     sd   = sd(log_load, na.rm=TRUE),
                     se   = sd / sqrt(N))
virload_table2

#### Multiple Regression: Effects of Host Plant, Survival, and Immune Assessment on Viral Load
## Final Model and Table
virload_lm_final <- lm(formula = log_load ~ host_plant+survived+assayed, data = aj_data_virload)
summary(virload_lm_final)
tab_model(virload_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, show.r2 = T,
          show.intercept = T, show.df = T, dv.labels = c("Viral load"), 
          pred.labels = c("(Intercept)","Host plant","Survived","Assayed"), string.stat = "t")
# Visualize normality of residuals
qqnorm(resid(virload_lm_final, type = "pearson"))
qqline(resid(virload_lm_final, type = "pearson"))
hist(resid(virload_lm_final,type = "pearson"))
# Visualize variance of residuals
virload_lm_resid <- resid(virload_lm_final, type = "pearson")
virload_lm_fitted <- fitted(virload_lm_final, type = "pearson")
plot(x = virload_lm_fitted, y = virload_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(virload_lm_resid ~ host_plant, data = aj_data_virload, main = "Host plant", ylab = "Residuals")
boxplot(virload_lm_resid ~ survived, data = aj_data_virload, main = "Survived", ylab = "Residuals")
boxplot(virload_lm_resid ~ assayed, data = aj_data_virload, main = "Assayed", ylab = "Residuals")

## Create Figure 3C
virload_table3 <- rbind(virload_table1,virload_table2)
figure3c <- ggline(virload_table3, x = "host_plant", y = "mean_load", group = "survived", 
                   plot_type = "l",size = .5, color = "survived",  ylim = c(-8.5,-2), 
                   xlab = "Host plant", ylab = "Relative viral load (log-transformed)") + 
                    font("xlab", size = 14) + font("ylab", size = 14) + 
                    font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                    scale_color_manual(limits = c("No","Overall","Yes"), 
                                       values=c("white", "white","white"), guide = "none")  +
                    scale_x_discrete(labels = c(expression(italic("Bacopa")),
                                                expression(italic("Plantago"))))  +
                    geom_errorbar(aes(ymin=mean_load-se, ymax=mean_load+se, group = survived), 
                                  width=.1, position = position_dodge(width=0.15)) + 
                    geom_point(color = "black", position = position_dodge(width=0.15), 
                               aes(fill = host_plant, shape = survived, group = survived, size = survived)) + 
                    theme(legend.text.align = 0,   legend.position = c(0.86,0.8)) +
                    scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") + 
                    scale_shape_manual(limits = c("No", "Overall","Yes"), values = c(24,21,25),
                                       labels = c("Larva/\\npupa", "Mean","Adult")) + 
                    scale_size_manual(limits = c("No", "Overall","Yes"), values = c(3.4,4,3.4), guide = "none") + 
                    guides(shape = guide_legend(title = "Life stage\\nat death"))
figure3c

## Export Figure 3C
tiff("Figure3C.tif", units = "in", width = 4.2, height = 4, res = 600)
figure3c
dev.off()

#### Logistic Regression: Effect of Viral Load*Host Plant Interaction and Immune Assessment on Survival
## Create Numeric Binary Variable (0/1) for Survival
aj_data_virload<- aj_data_virload %>%
  mutate(survived_num = case_when(
    survived == "No" ~ 0,
    survived == "Yes" ~ 1
  ))

## Model Selection
# Global Model
survload_glm_global <- glm(formula = survived_num ~ host_plant*log_load*assayed, 
                           family = binomial(link=logit), data = aj_data_virload, na.action = "na.fail")
survload_glm_global
# Model Comparison
survload_models <- dredge(survload_glm_global, fixed =c("host_plant","log_load","assayed"))
survload_models
# Single model (survived ~ host_plant*log_load + assayed) selected based on AICc/nesting rule

# Final Model and Table
survload_glm_final <- glm(formula = survived_num ~ host_plant*log_load+assayed, 
                          family = binomial(link=logit), data = aj_data_virload)
summary(survload_glm_final)
tab_model(survload_glm_final, show.est = T, show.ci = 0.95, show.se = F, show.stat = T, 
          show.r2 = T,show.intercept = T, dv.labels = c("Survival to adult"), 
          pred.labels = c("(Intercept)","Host plant","Viral load","Assayed", "Host plant x viral load"), 
          string.ci = "95% CI", string.est = "Odds Ratio", string.stat = "z")

## Create Figure 3D
figure3d <- interact_plot(survload_glm_final, pred = "log_load", modx = "host_plant", interval = T, 
                          int.type = c("prediction"), plot.points = T, point.alpha = .8, point.size = 2, 
                          int.width = 0.95, ylim = c(0,1.25), linetype= c("solid","longdash"),
                          x.label = "Relative viral load (log-transformed)", 
                          y.label = "Estimated survival probability", legend.main = "Host plant",
                          modx.labels = c("Bacopa","Plantago"), colors = c("LightSteelBlue3", "Brown4")) + 
                            theme_pubr()  + theme(legend.text.align = 0) +
                            theme(axis.text.x = element_text(size=12, color = "black")) + 
                            font("legend.text", size = 12) + 
                            theme(axis.text.y = element_text(size=12, color = "black")) + 
                            font("legend.title", size = 12, face = "bold") +
                            theme(legend.position = c(0.24,0.2)) + labs(color= "Host plant") + 
                            font("xlab", size = 14)+ font("ylab", size = 14) 
figure3d

## Export Figure 3d
tiff("Figure3D.tif", units = "in", width = 4.2, height = 4, res = 600)
figure3d
dev.off()


########## Section 4: Development Time ##########

#### Section Description
# Summary information and statistical analyses evaluating development of Anartia jatrophae individuals 
# based on host plant species (Bacopa monnieri or Plantago lanceolata), viral treatment (challenged 
# with Junonia coenia densovirus or unchallenged controls), and immune assay status (underwent immune 
# assays during the larval stage or did not undergo assays). Pre-inoculation development time (days 
# between molting to third instar and the second day of the sixth instar, when inoculation occurred), 
# is compared based on host plant and treatment group using two-way ANOVA. Effects of host plant, viral 
# treatment, and immune assessment on post-inoculation development time (days between inoculation and 
# eclosion) and pupal time (days in pupal stage) are assessed using separate multiple regression models.

#### Calculate Development Time Variables
## Pre-Inoculation Development Time
# Formula: n days between molting to third instar and inoculation
aj_data$date_l3 <- as.Date(aj_data$date_l3,"%m/%d/%Y")
aj_data$date_l6 <- as.Date(aj_data$date_l6,"%m/%d/%Y")
aj_data$days_to_inoc <- difftime(aj_data$date_l6,aj_data$date_l3, units = "days") + 1
aj_data$days_to_inoc <- as.numeric(aj_data$days_to_inoc)

## Post-Inoculation Development Time
# Formula: n days between inoculation (day 2, sixth instar) and eclosion
# Calculate days between inoculation and pupation
aj_data$date_pup <- as.Date(aj_data$date_pup,"%m/%d/%Y")
aj_data$inoc_to_pup <- difftime(aj_data$date_pup, aj_data$date_l6, units = "days") - 1
aj_data$inoc_to_pup <- as.numeric(aj_data$inoc_to_pup)
# Calculate days between pupation and eclosion
aj_data$date_ecl <- as.Date(aj_data$date_ecl,"%m/%d/%Y")
aj_data$pup_to_ecl <- difftime(aj_data$date_ecl, aj_data$date_pup, units = "days")
aj_data$pup_to_ecl <- as.numeric(aj_data$pup_to_ecl)
## Calculate Days between inoculation and eclosion
aj_data$inoc_to_ecl <- aj_data$inoc_to_pup + aj_data$pup_to_ecl

#### Summary Tables: Development Time
## Pre-inoc. Development Time Based on Host Plant/Future Viral Treatment
aj_data_dev <- aj_data[which(aj_data$all_exclude != "Yes"),]
devpre_table <- ddply(aj_data_dev, c("host_plant","treatment"), summarise, 
                      N = sum(!is.na(days_to_inoc)),
                      devpre = mean(days_to_inoc, na.rm=TRUE),
                      sd   = sd(days_to_inoc, na.rm=TRUE),
                      se   = sd / sqrt(N))

## Post-inoc. Larval Development Time Based on Host Plant/Viral Treatment
devpost_table <- ddply(aj_data_dev, c("host_plant","treatment"), summarise, 
                        N = sum(!is.na(inoc_to_pup)),
                        devpost = mean(inoc_to_pup, na.rm=TRUE),
                        sd   = sd(inoc_to_pup, na.rm=TRUE),
                        se   = sd / sqrt(N))

## Pupal Development Time Based on Host Plant/Viral Treatment
devpup_table <- ddply(aj_data_dev, c("host_plant","treatment"), summarise, 
                      N = sum(!is.na(pup_to_ecl)),
                      devpup = mean(pup_to_ecl, na.rm=TRUE),
                      sd   = sd(pup_to_ecl, na.rm=TRUE),
                      se   = sd / sqrt(N))

#### ANOVA: Pre-inoc. Development Time Based on Host Plant/Future Viral Treatment
## Final Model and Table
devpre_aov_final <- lm(days_to_inoc~host_plant+treatment, data = aj_data_dev)
Anova(devpre_aov_final, type = "II")
# Visualize normality of residuals
qqnorm(resid(devpre_aov_final, type = "pearson"))
qqline(resid(devpre_aov_final, type = "pearson"))
hist(resid(devpre_aov_final,type = "pearson"))
# Visualize variance of residuals
devpre_aov_resid <- resid(devpre_aov_final, type = "pearson")
devpre_aov_fitted <- fitted(devpre_aov_final, type = "pearson")
plot(x = devpre_aov_fitted, y = devpre_aov_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(devpre_aov_resid ~ host_plant, data = aj_data_dev, main = "Host plant", ylab = "Residuals")
boxplot(devpre_aov_resid ~ treatment, data = aj_data_dev, main = "Treatment", ylab = "Residuals")

#### Multiple Regression: Post-inoc. Development Time Based on Host Plant/Viral Treatment/Immune Assessment
## Model Selection
# Global Model
aj_data_devpost <- aj_data_dev[which(aj_data_dev$inoc_to_ecl != "NA"),]
devpost_lm_global <- lm(formula = inoc_to_ecl~host_plant*treatment+assayed, 
                        data = aj_data_devpost, na.action = "na.fail")
summary(devpost_lm_global)
# Model Comparison
devpost_models <- dredge(devpost_lm_global, fixed =c("host_plant","treatment","assayed"))
devpost_models
# Single model (inoc_to_ecl ~ host_plant+treatment+assayed) selected based on AICc/nesting rule

## Final Model and Table S3
devpost_lm_final <- lm(inoc_to_ecl~host_plant+treatment+assayed, data = aj_data_dev)
summary(devpost_lm_final)
tab_model(devpost_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, show.r2 = T,
          show.intercept = T, show.df = T, dv.labels = c("Post-inoc. development time"), 
          pred.labels = c("(Intercept)","Host plant","Treatment","Assayed"), string.stat = "t")
# Visualize normality of residuals
qqnorm(resid(devpost_lm_final, type = "pearson"))
qqline(resid(devpost_lm_final, type = "pearson"))
hist(resid(devpost_lm_final,type = "pearson"))
# Visualize variance of residuals
devpost_lm_resid <- resid(devpost_lm_final, type = "pearson")
devpost_lm_fitted <- fitted(devpost_lm_final, type = "pearson")
plot(x = devpost_lm_fitted, y = devpost_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(devpost_lm_resid ~ host_plant, data = aj_data_devpost, main = "Host plant", ylab = "Residuals")
boxplot(devpost_lm_resid ~ treatment, data = aj_data_devpost, main = "Treatment", ylab = "Residuals")
boxplot(devpost_lm_resid ~ assayed, data = aj_data_devpost, main = "Assayed", ylab = "Residuals")

#### Multiple Regression: Pupal Time Based on Host Plant/Viral Treatment/Immune Assessment
## Model Selection
# Global Model
aj_data_devpup <- aj_data_dev[which(aj_data_dev$pup_to_ecl != "NA"),]
devpup_lm_global <- lm(formula = pup_to_ecl~host_plant*treatment+assayed, 
                       data = aj_data_devpup, na.action = "na.fail")
summary(devpup_lm_global)
# Model Comparison
devpup_models <- dredge(devpup_lm_global, fixed =c("host_plant","treatment","assayed"))
devpup_models
# Single model (inoc_to_ecl ~ host_plant*treatment+assayed) selected based on AICc

## Final Model and Table
devpup_lm_final <- lm(pup_to_ecl~host_plant*treatment+assayed, data = aj_data_dev)
summary(devpup_lm_final)
tab_model(devpup_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, show.r2 = T,
          show.intercept = T, show.df = T, dv.labels = c("Pupal development time"), string.stat = "t",
          pred.labels = c("(Intercept)","Host plant","Treatment","Assayed", "Host plant x treatment"))
# Visualize normality of residuals
qqnorm(resid(devpup_lm_final, type = "pearson"))
qqline(resid(devpup_lm_final, type = "pearson"))
hist(resid(devpup_lm_final,type = "pearson"))
# Visualize variance of residuals
devpup_lm_resid <- resid(devpup_lm_final, type = "pearson")
devpup_lm_fitted <- fitted(devpup_lm_final, type = "pearson")
plot(x = devpup_lm_fitted, y = devpup_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(devpup_lm_resid ~ host_plant, data = aj_data_devpup, main = "Host plant", ylab = "Residuals")
boxplot(devpup_lm_resid ~ treatment, data = aj_data_devpup, main = "Treatment", ylab = "Residuals")
boxplot(devpup_lm_resid ~ assayed, data = aj_data_devpup, main = "Assayed", ylab = "Residuals")

## Pairwise Contrasts
emmeans(devpup_lm_final, pairwise ~ host_plant*treatment, adjust = "tukey", type = "response")

## Create Figure 4A
figure4a <- ggline(devpup_table, x = "treatment", y = "devpup", group = "host_plant", plot_type = "l",
                   color = "black", size = .5, linetype = "host_plant", ylim = c(7.75,8.85),
                   xlab = "Treatment", ylab = "Time in pupal stage (days)", legend = c(0.88,0.88), 
                   legend.title = "Host plant", position = position_dodge(0.125)) + 
                    font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                    font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                    scale_linetype_manual(values=c("solid", "longdash"), 
                                          labels = c(expression(italic("Bacopa")),
                                                     expression(italic("Plantago")))) + 
                    scale_x_discrete(labels = c("Control","Virus-challenged"))  +
                    geom_errorbar(aes(ymin=devpup-se, ymax=devpup+se, group = host_plant), 
                                  position = position_dodge(0.125), width=.05) +
                    geom_point(shape = 21, size = 4, color = "black", position = position_dodge(0.125),
                               aes(fill = host_plant, group = host_plant) ) +
                    scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") 
figure4a

## Export Figure 4A
tiff("Figure4a_fin.tif", units = "in", width = 4.2, height = 4, res = 600)
figure4a 
dev.off()


########## Section 5: Pupal Weight ##########

#### Section Description
# Summary information and statistical analyses evaluating variation in Anartia jatrophae pupal
# weights based on host plant species (Bacopa monnieri or Plantago lanceolata), viral treatment 
# (challenged with Junonia coenia densovirus or unchallenged controls), and immune assay status 
# (underwent immune assays during the larval stage or did not undergo assays). Multiple regression
# is used to assess the effects of host plant, viral treatment, and assay status on pupal weights;
# this model also includes sex (M/F) as a predictor to account for observed dimorphism in body size. 

#### Convert Pupal Weight Units 
## Convert from g to mg
aj_data_dev$pup_mass <- aj_data_dev$pup_mass*1000

#### Summary Tables: Pupal Weight
## Pupal Weight Based on Host Plant/Viral Treatment
pupwgt_table <- ddply(aj_data_dev, c("host_plant","treatment"), summarise, 
                      N = sum(!is.na(pup_mass)),
                      pup_wgt = mean(pup_mass, na.rm=TRUE),
                      sd   = sd(pup_mass, na.rm=TRUE),
                      se   = sd / sqrt(N))
pupwgt_table

#### Multiple Regression: Pupal Weight Based on Host Plant/Viral Treatment/Immune Assays/Sex
## Model Selection
# Global Model
aj_data_pupwgt <- aj_data_dev[which(aj_data_dev$pup_mass != "NA" & 
                                      aj_data_dev$sex != "NA"),]
pupwgt_lm_global <- lm(formula = pup_mass~host_plant*treatment+assayed+sex, 
                       data = aj_data_pupwgt, na.action = "na.fail")
summary(pupwgt_lm_global)
# Model Comparison
devpup_models <- dredge(pupwgt_lm_global, fixed =c("host_plant","treatment","assayed","sex"))
devpup_models
# Single model (pup_mass ~ host_plant+treatment+assayed+sex) selected based on AICc/nesting rule

## Final Model and Table S3
pupwgt_lm_final <- lm(pup_mass~host_plant+treatment+assayed+sex, data = aj_data_pupwgt)
summary(pupwgt_lm_final)
tab_model(pupwgt_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, show.r2 = T,
          show.intercept = T, show.df = T, dv.labels = c("Pupal weight"), string.stat = "t",
          pred.labels = c("(Intercept)","Host plant","Treatment","Assayed","Sex"))
# Visualize normality of residuals
qqnorm(resid(pupwgt_lm_final, type = "pearson"))
qqline(resid(pupwgt_lm_final, type = "pearson"))
hist(resid(pupwgt_lm_final,type = "pearson"))
# Visualize variance of residuals
pupwgt_lm_resid <- resid(pupwgt_lm_final, type = "pearson")
pupwgt_lm_fitted <- fitted(pupwgt_lm_final, type = "pearson")
plot(x = pupwgt_lm_fitted, y = pupwgt_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(pupwgt_lm_resid ~ host_plant, data = aj_data_pupwgt, main = "Host plant", ylab = "Residuals")
boxplot(pupwgt_lm_resid ~ treatment, data = aj_data_pupwgt, main = "Treatment", ylab = "Residuals")
boxplot(pupwgt_lm_resid ~ assayed, data = aj_data_pupwgt, main = "Assayed", ylab = "Residuals")
boxplot(pupwgt_lm_resid ~ sex, data = aj_data_pupwgt, main = "Sex", ylab = "Residuals")

## Pairwise Contrasts
emmeans(pupwgt_lm_final, pairwise ~ host_plant+treatment, adjust = "tukey", type = "response")
emmeans(pupwgt_lm_final, pairwise ~ host_plant, adjust = "tukey", type = "response")
# Host plant: (296-276)/276 = 7% higher on Plantago than Bacopa
emmeans(pupwgt_lm_final, pairwise ~ treatment, adjust = "tukey", type = "response")
# Viral treatment: (300-272)/272 = 10% higher in controls than virus-inoculated

## Create Figure 4B
figure4b <- ggline(pupwgt_table, x = "treatment", y = "pup_wgt", group = "host_plant", plot_type = "l",
                   color = "black", size = .5, linetype = "host_plant", legend = c(0.9,0.9), 
                   xlab = "Treatment", ylab = "Pupal weight (mg)", ylim = c(240,320),
                   legend.title = "Host plant", position = position_dodge(0.125)) + 
                    font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                    font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                    scale_linetype_manual(values=c("solid", "longdash"), 
                                          labels = c(expression(italic("Bacopa")),
                                                     expression(italic("Plantago")))) + 
                    scale_x_discrete(labels = c("Control","Virus-challenged"))  +
                    geom_errorbar(aes(ymin=pup_wgt-se, ymax=pup_wgt+se, group = host_plant), 
                                  position = position_dodge(0.125), width=.05) +
                    geom_point(shape = 21, size = 4, color = "black", position = position_dodge(0.125),
                               aes(fill = host_plant, group = host_plant)) +
                    scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") 
figure4b

## Export Figure 4B
tiff("Figure4b_fin.tif", units = "in", width = 4.2, height = 4, res = 600)
figure4b
dev.off()


########## Section 6: Immune Responses ##########

#### Section Description
# Summary information and statistical analyses evaluating variation in Anartia jatrophae immune
# responses based on host plant species (Bacopa monnieri or Plantago lanceolata) and viral treatment 
# (challenged with Junonia coenia densovirus or unchallenged controls). Immune responses (hemocyte 
# concentration and melanization score) are compared across host plants and viral treatments using 
# two-way ANOVAs. The outcomes of immune variation for survival of virus-challenged individuals is 
# assessed using logistic regression models that include host plant, either hemocyte concentration or 
# melanization (%), and their two-way interaction as predictors, while relationships between immunity 
# and postmortem viral loads were assessed using multiple regression models with the same predictors.

#### Transform Variables
## Convert Hemocyte Counts (cells/100 nl) to (cells/ml)
aj_data$hemo_conc <- aj_data$hemo_count*30000
# Square root-transform
aj_data$hemo_trans <- sqrt(aj_data$hemo_conc)

## Convert Implant MGV to Percent Melanization
# Formula: (1-(Sample MGV/Maximum MGV))*100
aj_data$mel_percent <- (1-(aj_data$implant_mgv/max(aj_data$implant_mgv, na.rm = T)))*100
# Square-transform
aj_data$mel_trans <- (aj_data$mel_percent)^2

#### Summary Tables: Immune Assays
## Hemocyte Concentration Based on Host Plant/Viral Treatment
aj_data_hemo <- aj_data[which(aj_data$all_exclude != "Yes" &
                               aj_data$hemo_exclude != "Yes"),]
hemo_table <- ddply(aj_data_hemo, c("host_plant","treatment"), summarise, 
                      N = sum(!is.na(hemo_conc)),
                      hemo = mean(hemo_conc, na.rm=TRUE),
                      sd   = sd(hemo_conc, na.rm=TRUE),
                      se   = sd / sqrt(N))
hemo_table

## Melanization Score Based on Host Plant/Viral Treatment
aj_data_mel<- aj_data[which(aj_data$all_exclude != "Yes" &
                              aj_data$implant_mgv != "NA"),]
mel_table <- ddply(aj_data_mel, c("host_plant","treatment"), summarise, 
                    N = sum(!is.na(mel_percent)),
                    mel = mean(mel_percent, na.rm=TRUE),
                    sd   = sd(mel_percent, na.rm=TRUE),
                    se   = sd / sqrt(N))
mel_table

#### ANOVA: Hemocyte Concentration Based on Host Plant/Viral Treatment
# Final Model and Table
hemo_aov_final <- lm(formula = hemo_trans ~ host_plant+treatment, data = aj_data_hemo)
Anova(hemo_aov_final, type = "II")
# Visualize normality of residuals
qqnorm(resid(hemo_aov_final, type = "pearson"))
qqline(resid(hemo_aov_final, type = "pearson"))
hist(resid(hemo_aov_final,type = "pearson"))
# Visualize variance of residuals
hemo_aov_resid <- resid(hemo_aov_final, type = "pearson")
hemo_aov_fitted <- fitted(hemo_aov_final, type = "pearson")
plot(x = hemo_aov_fitted, y = hemo_aov_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(hemo_aov_resid ~ host_plant, data = aj_data_hemo, main = "Host plant", ylab = "Residuals")
boxplot(hemo_aov_resid ~ treatment, data = aj_data_hemo, main = "Treatment", ylab = "Residuals")

## Create Figure 5A
figure5a <- ggline(hemo_table, x = "treatment", y = "hemo", group = "host_plant", plot_type = "l",
                   size = .5, color = "black", linetype = "host_plant", position = position_dodge(0.15),
                   xlab = "Treatment", ylab = "Hemocyte concentration (cells/ml)", legend = c(0.88,0.88), 
                   ylim = c(4000000,8000000), legend.title = "Host plant") + 
                    font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                    font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                    scale_linetype_manual(values=c("solid", "longdash"), 
                                          labels = c(expression(italic("Bacopa")),
                                                     expression(italic("Plantago")))) + 
                    scale_x_discrete(labels = c("Control","Virus-challenged")) + 
                    geom_errorbar(aes(ymin=hemo-se, ymax=hemo+se, group = host_plant), 
                                  position = position_dodge(0.15), width=.05) +
                    geom_point(shape = 21, size = 4, color = "black", position = position_dodge(0.15),
                               aes(fill = host_plant, group = host_plant)) + 
                    scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") 
figure5a

## Export Figure 5A
tiff("Figure5a.tif", units = "in", width = 4.2, height = 4, res = 600)
figure5a
dev.off()

#### Logistic Regression: Effect of Hemocyte Concentration on Survival of Viral Challenge
## Create Numeric Binary Variable (0/1) for Survival
aj_data_hemo<- aj_data_hemo %>%
    mutate(survived_num = case_when(
    survived == "No" ~ 0,
    survived == "Yes" ~ 1
    ))

## Final Model and Table S4
aj_data_hemovir <- aj_data_hemo[which(aj_data_hemo$treatment == "Virus"),]
hemosurv_glm <- glm(formula = survived_num ~ host_plant*hemo_trans,
                    family = binomial(link=logit), data = aj_data_hemovir)
summary(hemosurv_glm)
tab_model(hemosurv_glm, show.est = T, show.ci = 0.95, show.se = F, show.stat = T, digits = 3,
          show.r2 = T,show.intercept = T, dv.labels = c("Survival to adult"), pred.labels = 
            c("(Intercept)","Host plant","Hemocytes","Host plant x hemocytes"), 
          string.ci = "95% CI", string.est = "Odds Ratio", string.stat = "z")

#### Multiple Regression: Effect of Hemocyte Concentration on Viral Load
## Final Model and Table S4
hemoload_lm_final <- lm(formula = log_load ~ host_plant*hemo_trans, data = aj_data_hemovir)
summary(hemoload_lm_final)
tab_model(hemoload_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, digits = 4,
          show.r2 = T,show.intercept = T, dv.labels = c("Viral load"), pred.labels = 
          c("(Intercept)","Host plant","Hemocytes","Host plant x hemocytes"), string.stat = "t")
# Visualize normality of residuals
qqnorm(resid(hemoload_lm_final, type = "pearson"))
qqline(resid(hemoload_lm_final, type = "pearson"))
hist(resid(hemoload_lm_final,type = "pearson"))
# Visualize variance of residuals
aj_data_hemovir <- aj_data_hemovir[which(aj_data_hemovir$log_load != "NA"),]
hemoload_lm_resid <- resid(hemoload_lm_final, type = "pearson")
hemoload_lm_fitted <- fitted(hemoload_lm_final, type = "pearson")
plot(x = hemoload_lm_fitted, y = hemoload_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(hemoload_lm_resid ~ host_plant, data = aj_data_hemovir, main = "Host plant", ylab = "Residuals")
boxplot(hemoload_lm_resid ~ treatment, data = aj_data_hemovir, main = "Treatment", ylab = "Residuals")
plot(x = aj_data_hemovir$hemo_trans, y = hemoload_lm_resid, xlab = "Hemocytes", ylab = "Residuals")

#### ANOVA: Melanization Score Based on Host Plant/Treatment
## Final Model and Table
mel_aov_final <- lm(formula = mel_trans ~ host_plant+treatment, data = aj_data_mel)
Anova(mel_aov_final, type = "II")
# Visualize normality of residuals
qqnorm(resid(mel_aov_final, type = "pearson"))
qqline(resid(mel_aov_final, type = "pearson"))
hist(resid(mel_aov_final,type = "pearson"))
# Visualize variance of residuals
mel_aov_resid <- resid(mel_aov_final, type = "pearson")
mel_aov_fitted <- fitted(mel_aov_final, type = "pearson")
plot(x = mel_aov_fitted, y = mel_aov_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(mel_aov_resid ~ host_plant, data = aj_data_mel, main = "Host plant", ylab = "Residuals")
boxplot(mel_aov_resid ~ treatment, data = aj_data_mel, main = "Treatment", ylab = "Residuals")

## Create Figure 5B
figure5b <- ggline(mel_table, x = "treatment", y = "mel", group = "host_plant", plot_type = "l",
                   size = .5, color = "black", linetype = "host_plant", ylim = c(22,32),
                   xlab = "Treatment", ylab = "Melanization score (%)", legend = c(0.88,0.88), 
                   position = position_dodge(0.15),  legend.title = "Host plant") + 
                    font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                    font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                    scale_linetype_manual(values=c("solid", "longdash"), 
                                          labels = c(expression(italic("Bacopa")),
                                                     expression(italic("Plantago")))) + 
                    scale_x_discrete(labels = c("Control","Virus-challenged")) +
                    geom_errorbar(aes(ymin=mel-se, ymax=mel+se, group = host_plant), 
                                  position = position_dodge(0.15), width=.05) + 
                    geom_point(shape = 21, size = 4, color = "black", position = position_dodge(0.15),
                               aes(fill = host_plant, group = host_plant)) + 
                    scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") 
figure5b

## Export Figure 5B
tiff("Figure5b.tif", units = "in", width = 4.2, height = 4, res = 600)
figure5b
dev.off()

#### Logistic Regression: Effect of Melanization Score on Survival of Virus
## Create Numeric Binary Variable (0/1) for Survival
aj_data_mel<- aj_data_mel %>%
    mutate(survived_num = case_when(
    survived == "No" ~ 0,
    survived == "Yes" ~ 1
    ))

## Final Model and Table S4
aj_data_melvir <- aj_data_mel[which(aj_data_mel$treatment == "Virus"),]
melsurv_glm <- glm(formula = survived_num ~ host_plant*mel_trans,
                    family = binomial(link=logit), data = aj_data_melvir)
summary(melsurv_glm)
tab_model(melsurv_glm, show.est = T, show.ci = 0.95, show.se = F, show.stat = T, digits = 4,
          show.r2 = T,show.intercept = T, dv.labels = c("Survival to adult"), pred.labels = 
            c("(Intercept)","Host plant","Melanization","Host plant x melanization"), 
          string.ci = "95% CI", string.est = "Odds Ratio", string.stat = "z")

#### Multiple Regression: Effect of Melanization Score on Viral Load
## Final Model and Table S4
melload_lm_final <- lm(formula = log_load ~ host_plant*mel_trans, data = aj_data_melvir)
summary(melload_lm_final)
tab_model(melload_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, digits = 4,
          show.r2 = T,show.intercept = T, dv.labels = c("Viral load"), string.stat = "t",
          pred.labels = c("(Intercept)","Host plant","Melanization","Host plant x melanization"))
# Visualize normality of residuals
qqnorm(resid(melload_lm_final, type = "pearson"))
qqline(resid(melload_lm_final, type = "pearson"))
hist(resid(melload_lm_final,type = "pearson"))
# Visualize variance of residuals
aj_data_melvir <- aj_data_melvir[which(aj_data_melvir$log_load != "NA"),]
melload_lm_resid <- resid(melload_lm_final, type = "pearson")
melload_lm_fitted <- fitted(melload_lm_final, type = "pearson")
plot(x = melload_lm_fitted, y = melload_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(melload_lm_resid ~ host_plant, data = aj_data_melvir, main = "Host plant", ylab = "Residuals")
plot(x = aj_data_melvir$mel_trans, y = melload_lm_resid, xlab = "Melanization", ylab = "Residuals")


########## Section 7: Reproduction ##########

#### Section Description
# Summary information and statistical analyses evaluating variation in Anartia jatrophae reproductive
# parameters based on host plant species (Bacopa monnieri or Plantago lanceolata) and viral treatment 
# (challenged with Junonia coenia densovirus or unchallenged controls). Fecundity (total eggs laid over
# 3-4 days) and oviposition preference index (proportion of eggs laid on Bacopa) are compared across
# larval host plant species and viral treatments using separate multiple regression models, which include 
# host plant and viral treatment as predictors and number of oviposition days (3 or 4) as a covariate.

#### Calculate Oviposition Variables
## Calculate Total Eggs Laid
# Formula: total eggs = n eggs on Bacopa + n eggs on Plantago
aj_data$tot_eggs <- aj_data$bacopa_eggs + aj_data$plantago_eggs

## Calculate Oviposition Preference Index (OPI)
# Formula: OPI = (Bacopa eggs – Plantago eggs)/total eggs
aj_data$opi <- (aj_data$bacopa_eggs-aj_data$plantago_eggs)/aj_data$tot_eggs

## Calculate Number of Oviposition Days
aj_data$date_ovi_st <- as.Date(aj_data$date_ovi_st,"%m/%d/%Y")
aj_data$date_ovi_end <- as.Date(aj_data$date_ovi_end,"%m/%d/%Y")
aj_data$ovi_days <- difftime(aj_data$date_ovi_end, aj_data$date_ovi_st, units = "days")
aj_data$ovi_days <- as.numeric(aj_data$ovi_days)
# Subset to only include females with 3-4 days of oviposition
aj_data_ovi <- aj_data[which(aj_data$ovi_days>2 &
                              aj_data$all_exclude != "Yes" ),]

#### Summary Tables: Reproduction
## Total Eggs Laid Based on Host Plant/Viral Treatment
fecundity_table <- ddply(aj_data_ovi, c("host_plant","treatment"), summarise, 
                    N = sum(!is.na(tot_eggs)),
                    eggs = mean(tot_eggs, na.rm=TRUE),
                    sd   = sd(tot_eggs, na.rm=TRUE),
                    se   = sd / sqrt(N))
fecundity_table

## Oviposition Preference Based on Host Plant/Viral Treatment
opi_table <- ddply(aj_data_ovi, c("host_plant","treatment"), summarise, 
                   N = sum(!is.na(opi)),
                   pref = mean(opi, na.rm=TRUE),
                   sd   = sd(opi, na.rm=TRUE),
                   se   = sd / sqrt(N))
opi_table
# Calculate mean oviposition preference index across all groups
ddply(aj_data_ovi, c("survived"), summarise, 
                   N = sum(!is.na(opi)),
                   pref = mean(opi, na.rm=TRUE),
                   sd   = sd(opi, na.rm=TRUE),
                   se   = sd / sqrt(N))

## Presence of Eggs on Plantago/Bacopa (Y/N) Based on Host Plant/Viral Treatment
ddply(aj_data_ovi, c("host_plant","treatment"), summarise, 
                  N = sum(!is.na(tot_eggs)),
                  p_eggs = sum(plantago_eggs != 0, na.rm =T),
                  b_eggs = sum(bacopa_eggs != 0, na.rm =T))

#### Multiple Regression: Fecundity Based on Host Plant, Viral Treatment, and Oviposition Days
## Final Model and Table
fecundity_lm_final <- lm(tot_eggs ~ host_plant+treatment+ovi_days, aj_data_ovi)
summary(fecundity_lm_final)
tab_model(fecundity_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, show.df = T, 
          show.r2 = T,show.intercept = T, dv.labels = c("Total eggs"), pred.labels = 
          c("(Intercept)","Host plant","Treatment","Oviposition days"), string.stat = "t")
# Visualize normality of residuals
qqnorm(resid(fecundity_lm_final, type = "pearson"))
qqline(resid(fecundity_lm_final, type = "pearson"))
hist(resid(fecundity_lm_final,type = "pearson"))
# Visualize variance of residuals
fecundity_lm_resid <- resid(fecundity_lm_final, type = "pearson")
fecundity_lm_fitted <- fitted(fecundity_lm_final, type = "pearson")
plot(x = fecundity_lm_fitted, y = fecundity_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(fecundity_lm_resid ~ host_plant, data = aj_data_ovi, main = "Host plant", ylab = "Residuals")
boxplot(fecundity_lm_resid ~ treatment, data = aj_data_ovi, main = "Treatment", ylab = "Residuals")

## Create Figure 6
figure6 <- ggline(fecundity_table, x = "treatment", y = "eggs", group = "host_plant", plot_type = "l",
                  size = .5, color = "black", linetype = "host_plant", xlab = "Treatment", 
                  ylab = "Total eggs laid (n)", position = position_dodge(0.17), ylim = c(17,43),
                  legend = c(0.88,0.88), legend.title = "Host plant") + 
                  font("xlab", size = 14) + font("ylab", size = 14) + theme(legend.text.align = 0) +
                  font("legend.title", size = 12, face = "bold") + font("legend.text", size = 12) +
                  scale_linetype_manual(values=c("solid", "longdash"), 
                                        labels = c(expression(italic("Bacopa")),
                                                   expression(italic("Plantago"))))  + 
                  scale_x_discrete(labels = c("Control","Virus-challenged"))  +
                  geom_errorbar(aes(ymin=eggs-se, ymax=eggs+se, group = host_plant),
                                position = position_dodge(0.17), width=.05) +
                  geom_point(shape = 21, size = 4, color = "black", position = position_dodge(0.17),
                             aes(fill = host_plant, group = host_plant)) +
                  scale_fill_manual(values = c("lightsteelblue2", "brown4"), guide = "none") 
figure6

## Export Figure 6
tiff("Figure6.tif", units = "in", width = 4.2, height = 4, res = 600)
figure6
dev.off()

#### Multiple Regression: OPI Based on Host Plant, Viral Treatment, and Oviposition Days
## Final Model and Table
opi_lm_final <- lm(opi ~ host_plant+treatment+ovi_days, aj_data_ovi)
summary(opi_lm_final)
tab_model(opi_lm_final, show.est = T, show.ci = F, show.se = T, show.stat = T, string.stat = "t",  
          show.r2 = T,show.intercept = T, dv.labels = c("Oviposition preference"), show.df = T,
          pred.labels = c("(Intercept)","Host plant","Treatment","Oviposition days"), digits = 4)
# Visualize normality of residuals
qqnorm(resid(opi_lm_final, type = "pearson"))
qqline(resid(opi_lm_final, type = "pearson"))
hist(resid(opi_lm_final,type = "pearson"))
# Visualize variance of residuals
opi_lm_resid <- resid(opi_lm_final, type = "pearson")
opi_lm_fitted <- fitted(opi_lm_final, type = "pearson")
plot(x = opi_lm_fitted, y = opi_lm_resid, xlab = "Fitted Values", ylab = "Residuals")
boxplot(opi_lm_resid ~ host_plant, data = aj_data_ovi, main = "Host plant", ylab = "Residuals")
boxplot(opi_lm_resid ~ treatment, data = aj_data_ovi, main = "Treatment", ylab = "Residuals")


####################################################