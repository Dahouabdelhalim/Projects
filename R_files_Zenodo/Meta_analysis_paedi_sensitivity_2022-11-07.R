#### 1. - Check and set directory #### 
getwd()
#### 2. - Name file ####
##Name the file: 

####3. - Packages and functions ####

## packages already loaded through sourcing other document
# library(readxl)
# library(tidyverse)
# library(meta)
# library(metafor)
# library(janitor)

source("./Meta_analysis_paedi_2022-11-07.R")
rm(Paed_AE100, Paed_AE1000, Paed_prev, Paed)

# 4. Data preparation for stratification by methodology -------------------
## TT/GTT General population
Paed_adm_tt_g <- Paed_adm_tt %>%
  filter(Population %in% "General care")
## 13 Studies

AEadm_sub_tt_g <- metaprop(Paed_adm_tt_g$number_of_admissions_with_1_ae, Paed_adm_tt_g$sample_size, method.ci = "WS",
                              studlab = paste0(Paed_adm_tt_g$first_author_s_last_name, ", ",  Paed_adm_tt_g$publication_year,  "", Paed_adm_tt_g$Study_asterisk_AEadm), 
                              prediction = T, title = "% of admissions with >= 1 AE")

summary(AEadm_sub_tt_g)

# TT/GTT Intensive care population
Paed_adm_tt_i <- Paed_adm_tt %>%
  filter(Population %in% "Intensive care")

AEadm_sub_i <- metaprop(Paed_adm_tt_i$number_of_admissions_with_1_ae, Paed_adm_tt_i$sample_size, method.ci = "WS",
                              studlab = paste0(Paed_adm_tt_i$first_author_s_last_name, ", ",  Paed_adm_tt_i$publication_year,  "", Paed_adm_tt_i$Study_asterisk_AEadm),
                              prediction = T, title = "% of admissions with >= 1 AE")

summary(AEadm_sub_i)

## HMPS General population
Paed_adm_hmps_g <- Paed_adm_hmps %>%
  filter(Population %in% "General care")

AEadm_sub_hmps_g <- metaprop(Paed_adm_hmps_g$number_of_admissions_with_1_ae, Paed_adm_hmps_g$sample_size, method.ci = "WS",
                                   studlab = paste0(Paed_adm_hmps_g$first_author_s_last_name, ", ",  Paed_adm_hmps_g$publication_year,  "", Paed_adm_hmps_g$Study_asterisk_AEadm), 
                                   prediction = T, title = "% of admissions with >= 1 AE")

summary(AEadm_sub_hmps_g)

# 5. ------ Sensitivity analysis with study quality -----------------------
## This is done with the primary outcome measure.

# 5.1 Risk of Bias measures -----------------------------------------------
# 5.1.1 Risk of Bias -  Patient selection------

## TT/GTT General care sample
AEadm_RoB_Pat_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                   byvar = Paed_adm_tt_g$patient_selection_risk)
summary(AEadm_RoB_Pat_tt_g)
forest(AEadm_RoB_Pat_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_RoB_Pat_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$patient_selection_risk)
summary(AEadm_RoB_Pat_i)
forest(AEadm_RoB_Pat_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_RoB_Pat_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                        byvar = Paed_adm_hmps_g$patient_selection_risk)
summary(AEadm_RoB_Pat_hmps_g)
forest(AEadm_RoB_Pat_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))




# 5.1.2 Risk of Bias -  Reviewer------

## TT/GTT General care sample
AEadm_RoB_Rev_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$reviewer_risk)
summary(AEadm_RoB_Rev_tt_g)
forest(AEadm_RoB_Rev_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_RoB_Rev_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$reviewer_risk)
summary(AEadm_RoB_Rev_i)
forest(AEadm_RoB_Rev_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_RoB_Rev_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                        byvar = Paed_adm_hmps_g$reviewer_risk)
summary(AEadm_RoB_Rev_hmps_g)
forest(AEadm_RoB_Rev_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))


# 5.1.3 Risk of Bias -  Trigger tool methodology------

## TT/GTT General care sample
AEadm_RoB_RR_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$rr_process_risk)
summary(AEadm_RoB_RR_tt_g)
forest(AEadm_RoB_RR_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_RoB_RR_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$rr_process_risk)
summary(AEadm_RoB_RR_i)
forest(AEadm_RoB_RR_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_RoB_RR_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                       byvar = Paed_adm_hmps_g$rr_process_risk)
summary(AEadm_RoB_RR_hmps_g)
forest(AEadm_RoB_RR_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))



# 5.1.4 Risk of Bias -  Outcome------

## TT/GTT General care sample
AEadm_RoB_O_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$outcome_risk)
summary(AEadm_RoB_O_tt_g)
forest(AEadm_RoB_O_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_RoB_O_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$outcome_risk)
summary(AEadm_RoB_O_i)
forest(AEadm_RoB_O_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_RoB_O_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                      byvar = Paed_adm_hmps_g$outcome_risk)
summary(AEadm_RoB_O_hmps_g)
forest(AEadm_RoB_O_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))


# 5.1.5 Risk of Bias - Flow and Timing -----

## TT/GTT General care sample
AEadm_RoB_F_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$flow_timing_risk)
summary(AEadm_RoB_F_tt_g)
forest(AEadm_RoB_F_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_RoB_F_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$flow_timing_risk)
summary(AEadm_RoB_F_i)
forest(AEadm_RoB_F_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_RoB_F_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                      byvar = Paed_adm_hmps_g$flow_timing_risk)
summary(AEadm_RoB_F_hmps_g)
forest(AEadm_RoB_F_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))


rm(AEadm_RoB_Pat_tt_g, AEadm_RoB_Pat_i,
  AEadm_RoB_Rev_tt_g, AEadm_RoB_Rev_i, AEadm_RoB_RR_tt_g, AEadm_RoB_RR_i,
   AEadm_RoB_O_tt_g, AEadm_RoB_O_i, AEadm_RoB_F_tt_g, AEadm_RoB_F_i,
   AEadm_RoB_Pat_hmps_g, AEadm_RoB_Rev_hmps_g, AEadm_RoB_RR_hmps_g, AEadm_RoB_O_hmps_g,
  AEadm_RoB_F_hmps_g)


# 5.2 Concern regarding applicability -----------------------------------


# 5.2.1 Concern - Patient selection ---------------------------------------

## TT/GTT General care sample
AEadm_C_Pat_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$patient_selection_concern)
summary(AEadm_C_Pat_tt_g)
forest(AEadm_C_Pat_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_C_Pat_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$patient_selection_concern)
summary(AEadm_C_Pat_i)
forest(AEadm_C_Pat_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_C_Pat_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                      byvar = Paed_adm_hmps_g$patient_selection_concern)
summary(AEadm_C_Pat_hmps_g)
forest(AEadm_C_Pat_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

# 5.2.2 Concern - Reviewer ------------------------------------------------

## TT/GTT General care sample
AEadm_C_Rev_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_g$reviewer_concern)
summary(AEadm_C_Rev_tt_g)
forest(AEadm_C_Rev_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_C_Rev_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_tt_i$reviewer_concern)
summary(AEadm_C_Rev_i)
forest(AEadm_C_Rev_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_C_Rev_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                      byvar = Paed_adm_hmps_g$reviewer_concern)
summary(AEadm_C_Rev_hmps_g)
forest(AEadm_C_Rev_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

# 5.2.3 Concern - Trigger tool methodology --------------------------------

## TT/GTT General care sample
AEadm_C_RR_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                    byvar = Paed_adm_tt_g$rr_process_concern)
summary(AEadm_C_RR_tt_g)
forest(AEadm_C_RR_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_C_RR_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                    byvar = Paed_adm_tt_i$rr_process_concern)
summary(AEadm_C_RR_i)
forest(AEadm_C_RR_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_C_RR_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                     byvar = Paed_adm_hmps_g$rr_process_concern)
summary(AEadm_C_RR_hmps_g)
forest(AEadm_C_RR_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

# 5.2.4 Concern - Outcome -------------------------------------------------

## TT/GTT General care sample
AEadm_C_O_tt_g <- meta::update.meta(AEadm_sub_tt_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                   byvar = Paed_adm_tt_g$outcome_concern)
summary(AEadm_C_O_tt_g)
forest(AEadm_C_O_tt_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## TT/GTT Intensive care sample
AEadm_C_O_i <- meta::update.meta(AEadm_sub_i, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                   byvar = Paed_adm_tt_i$outcome_concern)
summary(AEadm_C_O_i)
forest(AEadm_C_O_i, study.results = F, overall = T, pscale = 100, xlim = c(0,100))

## HMPS General care sample
AEadm_C_O_hmps_g <- meta::update.meta(AEadm_sub_hmps_g, comb.random = T, comb.fixed = F, model.glmm = "binomial",
                                    byvar = Paed_adm_hmps_g$outcome_concern)
summary(AEadm_C_O_hmps_g)
forest(AEadm_C_O_hmps_g, study.results = F, overall = T, pscale = 100, xlim = c(0,100))



rm(AEadm_C_O_tt_g, AEadm_C_O_i,
   AEadm_C_Pat_tt_g, AEadm_C_Pat_i, AEadm_C_Rev_tt_g, AEadm_C_Rev_i,
   AEadm_C_RR_tt_g, AEadm_C_RR_i,
   AEadm_C_Pat_hmps_g, AEadm_C_O_hmps_g,
   AEadm_C_Rev_hmps_g, AEadm_C_RR_hmps_g)

##### This information was extracted to an excel file which is then loaded below to draw the forest plots.


# 5.3 Forest plots for quality analysis--------------------------------------------------------

# 5.3.1 TT/GTT General care sample -----------------------------------------------

## Loading the extracted numbers from the excel document
Forest_QAT_g <- read_excel("./stratanal_gtt_paedi_AEadm_QAT_2022-11-018.xls")

## rewrite p-value numbers into character variables
Forest_QAT_g <- Forest_QAT_g %>%
  mutate(pinteraction= replace(pinteraction, pinteraction %in% "1e-04", values = "<0.0001"))


# tiff("Plot7_tt_g_2022-09-15.tiff", width = 7.27, height = 8.69, units = 'in', res = 200)
# options(na.action = "na.pass")
# 
# forest(x=Forest_QAT_g$es, ci.lb = Forest_QAT_g$lcies, ci.ub = Forest_QAT_g$ucies,
#        xlim=c(-150,200),
#        slab=Forest_QAT_g$Stratum, at=seq(0,100, by=20),
#        ilab = cbind(Forest_QAT_g$nstudies, Forest_QAT_g$pinteraction), ilab.xpos = c(-30, 120),
#        refline = 17.7, digits = c(1,0), annosym = c(" [", ";", "]"),
#        xlab = NA, efac = NA, cex = 0.75, psize = 1,
#        rows = c(39, 38:35, 34:28, 27:24, 23:20, 19:16, 15:12, 11:8, 7:4, 3:1))
# op <- par(cex = 0.75, font = 2)
# text(-125, 42, "Type of Analysis", pos = 1)
# text(-30, 42, "N of studies", pos = 1)
# text(40, 42, "% of admissions with \\u2265 1 AE", pos = 1)
# text(167, 42, "Effect size [95% CI]", pos = 1)
# text(110, 42, "p interaction", pos = 1)
# par(op)
# 
# op <- par(cex = 0.75, font = 2)
# ### add text for the subgroups
# text(-153, c(39, 37, 33, 29,25,21,16,12,8,4), pos=4,
#      c("GTT/TT General care — Overall", "Patient selection", "Reviewer", "Record review process",
#        "Outcomes","Flow","Patient selection",
#        "Reviewer", "Record review process", "Outcomes"))
# text(-155, c(38, 17), pos = 4, font = 4,
#      c("Risk of bias", "Applicability-related concerns"))
# par(op)
# 
# dev.off()


# 5.3.2 TT/GTT Intensive care population -----------------------------------------

## Load data from excel sheet
Forest_QAT_i <- read_excel("./stratanal_gtt_paedi_AEadm_QAT_2022-11-18.xls",
                           sheet = "TTGTT Intensive care population")
## rewrite p-value numbers into character varables
Forest_QAT_i <- Forest_QAT_i %>%
  mutate(pinteraction= replace(pinteraction, pinteraction %in% "1e-04", values = "<0.0001"))


# tiff("Plot7_tt_i_2022-09-15.tiff", width = 7.27, height = 8.69, units = 'in', res = 200)
# options(na.action = "na.pass")
# 
# forest(x=Forest_QAT_i$es, ci.lb = Forest_QAT_i$lcies, ci.ub = Forest_QAT_i$ucies,
#        xlim=c(-150,200),
#        slab=Forest_QAT_i$Stratum, at=seq(0,100, by=20),
#        ilab = cbind(Forest_QAT_i$nstudies, Forest_QAT_i$pinteraction), ilab.xpos = c(-30, 120),
#        refline = 47.3, digits = c(1,0), annosym = c(" [", "; ", "]"),
#        xlab = NA, efac = NA, cex = 0.75, psize = 1)
# op <- par(cex = 0.75, font = 2)
# text(-125, 39, "Type of Analysis", pos = 1)
# text(-30, 39, "N of studies", pos = 1)
# text(40, 39, "% of admissions with \\u2265 1 AE", pos = 1)
# text(167, 39, "Effect size [95% CI]", pos = 1)
# text(110, 39, "p interaction", pos = 1)
# par(op)
# 
# op <- par(cex = 0.75, font = 2)
# text(-151, c(36, 34, 31, 27,24,20,15,11,8,4), pos=4,
#      c("GTT/TT Intensive care — Overall", "Patient selection", "Reviewer", "Record review process",
#        "Outcomes","Flow","Patient selection",
#        "Reviewer", "Record review process", "Outcomes"))
# 
# text(-153, c(35, 16), pos = 4, font = 4,
#      c("Risk of bias", "Applicability-related concerns"))
# par(op)
# dev.off()

# 5.3.3 HMPS General population ------------------------------------------------

## Loading the extracted numbers from the excel document
Forest_QAT_g_hmps <- read_excel("./stratanal_gtt_paedi_AEadm_QAT_2022-11-18.xls",
                           sheet = "HMPS General population")
## rewrite p-value numbers into character varables
Forest_QAT_g_hmps <- Forest_QAT_g_hmps %>%
  mutate(pinteraction= replace(pinteraction, pinteraction %in% "1e-04", values = "<0.0001"))

# tiff("Plot7_hmps_g_2022-09-03.tiff", width = 7.27, height = 8.69, units = 'in', res = 200)
# options(na.action = "na.pass")
# 
# forest(x=Forest_QAT_g_hmps$es, ci.lb = Forest_QAT_g_hmps$lcies, ci.ub = Forest_QAT_g_hmps$ucies,
#        xlim=c(-150,200),
#        slab=Forest_QAT_g_hmps$Stratum, at=seq(0,100, by=20),
#        ilab = cbind(Forest_QAT_g_hmps$nstudies, Forest_QAT_g_hmps$pinteraction), ilab.xpos = c(-30, 120),
#        refline = 3.93, digits = c(1,0), annosym = c(" [", "; ", "]"),
#        xlab = NA, efac = NA, cex = 0.75, psize = 1)
# op <- par(cex = 0.75, font = 2)
# text(-125, 36, "Type of Analysis", pos = 1)
# text(-30, 36, "N of studies", pos = 1)
# text(40, 36, "% of admissions with \\u2265 1 AE", pos = 1)
# text(167, 36, "Effect size [95% CI]", pos = 1)
# text(110, 36, "p interaction", pos = 1)
# par(op)
# 
# op <- par(cex = 0.75, font = 2)
# text(-151, c(33, 31, 28, 24,20,16,11,9,6,3), pos=4,
#      c("HMPS General care — Overall", "Patient selection", "Reviewer", "Record review process",
#        "Outcomes","Flow","Patient selection",
#        "Reviewer", "Record review process", "Outcomes"))
# text(-153, c(32, 12), pos = 4, font = 4,
#      c("Risk of bias", "Applicability-related concerns"))
# 
# op <- par(cex = 0.75, font = 1)
# # text(177.8, 21.8, "0 [              ]", pos = 1)
# 
# par(op)
# dev.off()
