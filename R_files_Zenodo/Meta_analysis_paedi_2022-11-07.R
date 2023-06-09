#### 1. - Check and set directory #### 
getwd()
#### 2. - Name file ####
##Name the file: 

####3. - Packages and functions ####
library(readxl)
library(tidyverse)
library(meta)
library(metafor)
library(janitor)
##Attaching package: ‘janitor’
# The following objects are masked from ‘package:stats’:
#   
#   chisq.test, fisher.test

source("./Meta_analysis_paedi_prep_2022-11-07.R")
rm(Paed9, Paed2, Paed10, Paed_GTT, Paed_HMPS)

# 4. ------ Meta-Analysis primary outcome - % of admissions with AEs --------------------------------------------------------
## Meta package, metaprop function, with the Wilson method for the CI

# 4.1 % of admissions --------------------------
Paed_adm <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_ae))
## 32 out of 33 studies have information on the outcome

## meta analysis with by-variable population
AEadm_sub <- metaprop(Paed_adm$number_of_admissions_with_1_ae, Paed_adm$sample_size, method.ci = "WS",
                        studlab = paste0(Paed_adm$first_author_s_last_name, ", ",  Paed_adm$publication_year,  "", Paed_adm$Study_asterisk_AEadm), 
                        prediction = T, title = "% of admissions with >= 1 AE", byvar = Paed_adm$Population)

summary(AEadm_sub)

## looking at forest plot and preparing the output
forest.meta(AEadm_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_adm$publication_year,
                  leftlabs = c("Study", "N of admissions", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)

# 4.1.1 Forest plot % of admissions ---------------------------------------

### Creation of the forest plot and saving to the project

# tiff("Plot1_2022-10-26.tiff", width = 8, height =9.5, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_sub, fontsize = 9, xlim = c(0,100), digits = 1, big.mark = ",", comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm$sample_size,
#                   leftlabs = c("Study", "N of admissions", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("with \\u2265 1 AE  ",  .39, .91, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("with \\u2265 1 AE   ",  .84, .91, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions   ",  .84, .93, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .084, .385, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 0.9; 57.0]   ",  .9295, .385, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .084, .09, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 6.9; 91.6]   ",  .9295, .09, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()

# rm(AEadm_sub)


# 4.2 Subgroups separately ------------------------------------------------

## General population
Paed_adm_g <- Paed_adm %>%
  filter(Population %in% "General care")

AEadm_sub_g <- metaprop(Paed_adm_g$number_of_admissions_with_1_ae, Paed_adm_g$sample_size, method.ci = "WS",
                            studlab = paste0(Paed_adm_g$first_author_s_last_name, ", ",  Paed_adm_g$publication_year,  "", Paed_adm_g$Study_asterisk_AEadm), 
                            prediction = T, title = "% of admissions with >= 1 AE")

summary(AEadm_sub_g)

## Intensive care population
Paed_adm_i <- Paed_adm %>%
  filter(Population %in% "Intensive care")

AEadm_sub_i <- metaprop(Paed_adm_i$number_of_admissions_with_1_ae, Paed_adm_i$sample_size, method.ci = "WS",
                              studlab = paste0(Paed_adm_i$first_author_s_last_name, ", ",  Paed_adm_i$publication_year,  "", Paed_adm_i$Study_asterisk_AEadm), 
                              prediction = T, title = "% of admissions with >= 1 AE")

summary(AEadm_sub_i)



# 4.3 Publication bias assessment -----------------------------------------

funnel.meta(AEadm_sub, comb.random = T, comb.fixed = F,
                  level=0.95, xlab = "Logit of AE incidence", studlab = F)

# tiff("Plot7_2022-06-17.tiff", width = 7.27, height = 4.69, units = 'in', res = 200)
# meta::funnel.meta(AE100, comb.random = T, comb.fixed = F,
#                  level=0.95, xlab = "Logit of AE incidence", studlab = F)
# title("% of admissions with AEs")
# dev.off()

rm(AEadm_sub, AEadm_sub_g, AEadm_sub_i)
rm(Paed_adm_g, Paed_adm_i)

# 4.4 Stratified by RRR Method---------------------------------------------------------------------

# 4.4.1 TT / GTT --------------------------
Paed_adm_tt <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_ae)) %>%
  filter(rrr_method %in% "Trigger tool (TT)" | rrr_method %in% "Global Trigger tool (GTT)")
## 23 out of 33 studies have information on the outcome

AEadm_tt_sub <- metaprop(Paed_adm_tt$number_of_admissions_with_1_ae, Paed_adm_tt$sample_size, method.ci = "WS",
                            studlab = paste0(Paed_adm_tt$first_author_s_last_name, ", ",  Paed_adm_tt$publication_year,  "", Paed_adm_tt$Study_asterisk_AEadm), 
                            prediction = T, title = "% of admissions with >= 1 AE for Trigger Tools", byvar = Paed_adm_tt$Population)

summary(AEadm_tt_sub)

forest.meta(AEadm_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_adm_tt$sample_size,
                  leftlabs = c("Study", "N of admissions", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)


# 4.4.1.1 Forest plot % of admissions ---------------------------------------

# tiff("Plot1_tt_2022-10-26.tiff", width = 8, height =9.5, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, big.mark = ",", comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm_tt$sample_size,
#                   leftlabs = c("GTT/TT Studies", "N of admissions", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("with \\u2265 1 AE  ",  .39, .817, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("with \\u2265 1 AE   ",  .84, .817, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions   ",  .84, .837, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .084, .48, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 3.8; 53.8]   ",  .9295, .48, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .084, .185, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 6.9; 91.6]   ",  .9295, .185, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()

rm(AEadm_tt_sub)

# 4.4.2 HMPS --------------------------
Paed_adm_hmps <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_ae)) %>%
  filter(rrr_method %in% "HMPS")
## 9 out of 33 studies have information on the outcome

AEadm_hmps_sub <- metaprop(Paed_adm_hmps$number_of_admissions_with_1_ae, Paed_adm_hmps$sample_size, method.ci = "WS",
                               studlab = paste0(Paed_adm_hmps$first_author_s_last_name, ", ",  Paed_adm_hmps$publication_year,  "", Paed_adm_hmps$Study_asterisk_AEadm), 
                               prediction = T, title = "% of admissions with >= 1 AE for Trigger Tools", byvar = Paed_adm_hmps$Population)

summary(AEadm_hmps_sub)

forest.meta(AEadm_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_adm_hmps$sample_size,
                  leftlabs = c("Study", "N of admissions", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)


# 4.4.1.1 Forest plot % of admissions ---------------------------------------

# tiff("Plot1_hmps_2022-10-26.tiff", width = 8, height =4, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, big.mark = ",", comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm_hmps$sample_size,
#                   leftlabs = c("HMPS Studies", "N of admissions", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("with \\u2265 1 AE  ",  .39, .813, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("with \\u2265 1 AE   ",  .84, .813, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions   ",  .84, .85, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .085, .21, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 0.3; 33.7]   ",  .9293, .21, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()

rm(AEadm_hmps_sub)

# 5. ------ Secondary outcomes - AEs per 100 admissions ----------------------------------------
## Meta package, metarate function. 
## The metarate function will use the Poisson-Normal model (rma.glmm, with "IRLN") from the metafor package (when it is loaded).

# 5.1 AEs per 100 admissions  ----------------------------------------------
Paed_AE100 <- Paed11 %>%
  filter(inclusion_of_more_than_one_ae_per_patient == "Yes") %>% ## 26 out of 33 studies selected (exclusion of HMPS-studies)
  filter(!(is.na(number_of_a_es))) %>% ## excludes 1 additional study (Zegers, no data available) --> 25 out of 33
  filter(!study_id == 29) ## excluding all HMPS studies (Requena) -->24 out of 33

## meta-analysis 
AE100_meta <- metarate(number_of_a_es, sample_size, data=Paed_AE100, method = "GLMM",
                             studlab = paste0(Paed_AE100$first_author_s_last_name, ", ",  Paed_AE100$publication_year,  "", Paed_AE100$Study_asterisk_AE100), 
                             prediction = T, irscale = 100,
                             title = "AEs per 100 admissions", 
                             byvar = Paed_AE100$Population)
summary(AE100_meta)

## forest plot
meta::forest.meta(AE100_meta, fontsize = 9, digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, smlab = "", test.subgroup = F, overall = F,
                  leftcols= c("studlab", "event", "time"),
                  sortvar = Paed_AE100$sample_size,
                  leftlabs = c("Study", "N of AEs", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  prediction = T,
                  xlim = c(0,370))

# 5.1.1 Forest plot AE/100 admissions -------------------------------------
# tiff("Plot2_2022-09-15.tiff", width = 7, height =7.5, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# meta::forest.meta(AE100_meta, fontsize = 9, xlim = c(0,370), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, smlab = "", test.subgroup = F, overall = F,
#                   leftcols= c("studlab", "event", "time"),
#                   sortvar = Paed_AE100$sample_size,
#                   leftlabs = c("GTT/TT Studies", "N of AEs", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   prediction = F)
# grid::grid.text("AEs per 100 admissions   ",  .855, .941, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .042, .488, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 4.2; 145.2]   ",  .973, .488, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .042, .088, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 15.3; 699.7]   ",  .973, .088, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# 
# dev.off()

# 5.2 AEs/100 admissions subgroups separately -------------------------------
## General population
Paed_AE100_g <- Paed_AE100 %>%
  filter(Population %in% "General care") %>%
  filter(!(is.na(number_of_a_es)))

AE100_sub_g <- metarate(number_of_a_es, sample_size, data=Paed_AE100_g, method = "GLMM",
                              studlab = paste0(Paed_AE100_g$first_author_s_last_name, ", ",  Paed_AE100_g$publication_year,  "", Paed_AE100_g$Study_asterisk_AE100),
                              prediction = T, irscale = 100,
                              title = "AEs per 100 admissions")
summary(AE100_sub_g)

## Intensive care population
Paed_AE100_i <- Paed_AE100 %>%
  filter(Population %in% "Intensive care") %>%
  filter(!(is.na(number_of_a_es)))

AE100_sub_i <-metarate(number_of_a_es, sample_size, data=Paed_AE100_i, method = "GLMM",
                              studlab = paste0(Paed_AE100_i$first_author_s_last_name, ", ",  Paed_AE100_i$publication_year,  "", Paed_AE100_i$Study_asterisk_AE100),
                              prediction = T, irscale = 100,
                              title = "AEs per 100 admissions")
summary(AE100_sub_i)

rm(AE100_sub_g, AE100_sub_i, Paed_AE100_g, Paed_AE100_i)

# 6. ------ AEs per 1000 patient days -------------------------------------
## Meta package, metarate function. 
## Metarate function will use the Poisson-Normal model (rma.glmm, with "IRLN") from the metafor package (when it is loaded)

Paed_AE1000 <- Paed11 %>%
  filter(inclusion_of_more_than_one_ae_per_patient == "Yes") %>% ## 26 out of 33 studies selected (exclusion of HMPS-studies)
  filter(!(is.na(number_of_hospital_days))) %>% ## excludes 3 additional study (Chapman, Maziero and Zegers, no data available) --> 23 out of 33
  filter(!study_id == 29) ## excluding HMPS (Requena) --> 22 out of 33

# 6.1 AEs per 1000 patient days ----------------------------------------------

AE1000_meta <- metarate(number_of_a_es, number_of_hospital_days, data=Paed_AE1000, method = "GLMM",
                             studlab = paste0(Paed_AE1000$first_author_s_last_name, ", ",  Paed_AE1000$publication_year,  "", Paed_AE1000$Study_asterisk_AE1000), 
                             prediction = T, irscale = 1000,
                             title = "AEs per 1000 patient days", 
                             byvar = Paed_AE1000$Population)
summary(AE1000_meta)

## rounding to full days
AE1000_meta$time <- round(AE1000_meta$time, digits = 0)

## forestplot
forest.meta(AE1000_meta, fontsize = 9, digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, smlab = "", test.subgroup = F, overall = F,
                  leftcols= c("studlab", "event", "time"),
                  sortvar = Paed_AE1000$number_of_hospital_days,
                  leftlabs = c("Study", "N of AEs", "Patient days"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  prediction = T)

# 6.1.1 Forest plot AE/1000 patient days -------------------------------------
# 
# tiff("Plot3_2022-09-15.tiff", width = 7, height =7, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# meta::forest.meta(AE1000_meta, fontsize = 9, xlim = c(0,650), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, smlab = "", test.subgroup = F, overall = F,
#                   leftcols= c("studlab", "event", "time"),
#                   sortvar = Paed_AE1000$number_of_hospital_days,
#                   leftlabs = c("GTT/TT Studies", "N of AEs", "Patient days"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   prediction = F)
# grid::grid.text("AEs per 1000 patient days   ",  .845, .943, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("[ 5.9; 393.1]   ",  .975, .485 , gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .039, .485, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 6.4; 2495.1]   ",  .975, .088, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .039, .088, gp = grid::gpar(fontsize=9, font = 2), just="left")
# 
# dev.off()

# 6.2 AEs/1000 patient days subgroups-------------------------------
## General population
Paed_AE1000_g <- Paed_AE1000 %>%
  filter(Population %in% "General care") %>%
  filter(!(is.na(number_of_hospital_days)))

AE1000_sub_g <- metarate(number_of_a_es, number_of_hospital_days, data=Paed_AE1000_g, method = "GLMM",
                              studlab = paste0(Paed_AE1000_g$first_author_s_last_name, ", ",  Paed_AE1000_g$publication_year,  "", Paed_AE1000_g$Study_asterisk_AE100),
                              prediction = T, irscale = 1000,
                              title = "AEs per 1000 patient days")
summary(AE1000_sub_g)

## Intensive care population
Paed_AE1000_i <- Paed_AE1000 %>%
  filter(Population %in% "Intensive care") %>%
  filter(!(is.na(number_of_hospital_days)))

AE1000_sub_i <- metarate(number_of_a_es, number_of_hospital_days, data=Paed_AE1000_i, method = "GLMM",
                              studlab = paste0(Paed_AE1000_i$first_author_s_last_name, ", ",  Paed_AE1000_i$publication_year,  "", Paed_AE1000_i$Study_asterisk_AE100),
                              prediction = T, irscale = 1000,
                              title = "AEs per 1000 patient days")
summary(AE1000_sub_i)

rm(AE1000_sub_g, AE1000_sub_i, Paed_AE1000_g, Paed_AE1000_i)


# 7. ------ % of preventable AEs -------------------------------------------------
## Meta package, metaprop function, with the Wilson method for the CI

# 7.1 % of preventable AEs --------------------------
Paed_prev <- Paed11 %>%
  filter(!is.na(no_of_preventable_a_es))
## 16 out of 33 studies have information on the outcome

AEprev_sub <- metaprop(Paed_prev$no_of_preventable_a_es, Paed_prev$number_of_a_es, method.ci = "WS",
                            studlab = paste0(Paed_prev$first_author_s_last_name, ", ",  Paed_prev$publication_year,  "", Paed_prev$Study_asterisk_prevAE), 
                            prediction = T, title = "% of AE preventable", byvar = Paed_prev$Population)

summary(AEprev_sub)

forest.meta(AEprev_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_prev$number_of_a_es,
                  leftlabs = c("Study", "N of preventable AEs", "N of AEs"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1)


# 7.1.1 Forest plot % of preventable AEs ---------------------------------------

# tiff("Plot4_2022-09-03.tiff", width = 7, height =6, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEprev_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_prev$number_of_a_es,
#                   leftlabs = c("Study", "N of preven-", "N of AEs"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("table AEs  ",  .37, .885, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of preventable AEs   ",  .855, .915, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("[20.5; 86.1]   ",  .959, .48 , gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .055, .48, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[  4.5; 98.9]   ",  .959, .12, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .055, .12, gp = grid::gpar(fontsize=9, font = 2), just="left")
# 
# dev.off()


# 7.2 % of preventable AEs subgroups --------------------------------------
## General care
Paed_prev_g <- Paed_prev %>%
  filter(Population %in% "General care") %>%
  filter(!is.na(no_of_preventable_a_es))

AEprev_sub_g <- metaprop(Paed_prev_g$no_of_preventable_a_es, Paed_prev_g$number_of_a_es, method.ci = "WS",
                             studlab = paste0(Paed_prev_g$first_author_s_last_name, ", ",  Paed_prev_g$publication_year,  "", Paed_prev_g$Study_asterisk_prevAE),
                             prediction = T, title = "% of AE preventable")

summary(AEprev_sub_g)

## Intensive care
Paed_prev_i <- Paed_prev %>%
  filter(Population %in% "Intensive care") %>%
  filter(!is.na(no_of_preventable_a_es))

AEprev_sub_i <- metaprop(Paed_prev_i$no_of_preventable_a_es, Paed_prev_i$number_of_a_es, method.ci = "WS",
                               studlab = paste0(Paed_prev_i$first_author_s_last_name, ", ",  Paed_prev_i$publication_year,  "", Paed_prev_i$Study_asterisk_prevAE),
                               prediction = T, title = "% of AE preventable")

summary(AEprev_sub_i)

rm(Paed_prev_g, AEprev_sub_g, Paed_prev_i, AEprev_sub_i)


# 7.3 Stratified by RRR Method---------------------------------------------------------------------

# 7.3.1 TT / GTT --------------------------
Paed_prev_tt <- Paed11 %>%
  filter(!is.na(no_of_preventable_a_es)) %>%
  filter(rrr_method %in% "Trigger tool (TT)" | rrr_method %in% "Global Trigger tool (GTT)")
## 11 out of 33 studies have information on the outcome

AEprev_tt_sub <- metaprop(Paed_prev_tt$no_of_preventable_a_es, Paed_prev_tt$number_of_a_es, method.ci = "WS",
                             studlab = paste0(Paed_prev_tt$first_author_s_last_name, ", ",  Paed_prev_tt$publication_year,  "", Paed_prev_tt$Study_asterisk_prevAE), 
                             prediction = T, title = "% of AE preventable", byvar = Paed_prev_tt$Population)

summary(AEprev_tt_sub)

forest.meta(AEprev_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_prev_tt$number_of_a_es,
                  leftlabs = c("Study", "N of preventable AEs", "N of AEs"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1)


# 7.3.1.1 Forest plot % of preventable AEs ---------------------------------------

# tiff("Plot4_tt_2022-09-15.tiff", width = 7, height =5, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEprev_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_prev_tt$number_of_a_es,
#                   leftlabs = c("GTT/TT Studies", "N of preven-", "N of AEs"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("table AEs  ",  .37, .865, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of preventable AEs   ",  .855, .9, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("[7.4; 96.2]   ",  .959, .58 , gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .055, .58, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[  4.5; 98.9]   ",  .959, .14, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .055, .14, gp = grid::gpar(fontsize=9, font = 2), just="left")
# 
# dev.off()

rm(AEprev_tt_sub, Paed_prev_tt)

# 7.3.2 HMPS --------------------------
Paed_prev_hmps <- Paed11 %>%
  filter(!is.na(no_of_preventable_a_es)) %>%
  filter(rrr_method %in% "HMPS")
## 5 out of 33 studies have information on the outcome

AEprev_hmps_sub <- metaprop(Paed_prev_hmps$no_of_preventable_a_es, Paed_prev_hmps$number_of_a_es, method.ci = "WS",
                                  studlab = paste0(Paed_prev_hmps$first_author_s_last_name, ", ",  Paed_prev_hmps$publication_year,  "", Paed_prev_hmps$Study_asterisk_prevAE), 
                                  prediction = T, title = "% of AE preventable", byvar = Paed_prev_hmps$Population)

summary(AEprev_hmps_sub)

forest.meta(AEprev_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F, 
                  sortvar = Paed_prev_hmps$number_of_a_es,
                  leftlabs = c("Study", "N of preventable AEs", "N of AEs"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1)


# 7.3.1.1 Forest plot % of preventable AEs ---------------------------------------

# tiff("Plot4_hmps_2022-09-03.tiff", width = 7, height =3, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEprev_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_prev_hmps$number_of_a_es,
#                   leftlabs = c("HMPS Studies", "N of preven-", "N of AEs"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("table AEs  ",  .37, .77, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of preventable AEs   ",  .855, .834, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("[  10.4; 91.8]   ",  .959, .235, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .055, .235, gp = grid::gpar(fontsize=9, font = 2), just="left")
# 
# dev.off()
rm(AEprev_hmps_sub, Paed_prev_hmps)



# 8. ------ % of admissions with preventable AEs --------------------------
## Meta package, metaprop function, with the Wilson method for the CI

# 8.1 % of admissions with preventable AEs --------------------------
Paed_adm_prev <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_preventable_ae))
## 12 out of 33 studies have information on the outcome

AEadm_prev <- metaprop(Paed_adm_prev$number_of_admissions_with_1_preventable_ae, Paed_adm_prev$sample_size, method.ci = "WS",
                            studlab = paste0(Paed_adm_prev$first_author_s_last_name, ", ",  Paed_adm_prev$publication_year,  "", Paed_adm_prev$Study_asterisk_prevAE),
                            prediction = T, title = "% of admissions preventable AEs", byvar = Paed_adm_prev$Population)

summary(AEadm_prev)

forest.meta(AEadm_prev, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
                  sortvar = Paed_adm_prev$sample_size,
                  leftlabs = c("Study", "N of preventable", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)

# 8.1.1 Forest plot % of admissions with preventable AEs ---------------------------------------
# 
# tiff("Plot6_2022-06-20.tiff", width = 8, height =6.5, units = 'in', res = 200)
# # par(mar=c(4,4,1,2))
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_prev, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm_prev$sample_size,
#                   leftlabs = c("Study", "N of admissions with", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1,
#                   prediction = F)
# grid::grid.text("preventable AEs",  .397, .79, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions with",  .844, .822, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("preventable AEs",  .844, .79, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .066, .485, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 0.2; 48.1]   ",  .947, .485, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .066, .21, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[ 2.5; 81.3]   ",  .947, .21, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()


# 8.2 % of adm with prev AEs subgroups ------------------------------------
## General care
Paed_adm_prev_g <- Paed_adm_prev %>%
  filter(Population %in% "General care")

AEadm_prev_g <- metaprop(Paed_adm_prev_g$number_of_admissions_with_1_preventable_ae, Paed_adm_prev_g$sample_size, method.ci = "WS",
                             studlab = paste0(Paed_adm_prev_g$first_author_s_last_name, ", ",  Paed_adm_prev_g$publication_year,  "", Paed_adm_prev_g$Study_asterisk_prevAE),
                             prediction = T, title = "% of admissions preventable AEs")

summary(AEadm_prev_g)

## Intensive care
Paed_adm_prev_i <- Paed_adm_prev %>%
  filter(Population %in% "Intensive care")

AEadm_prev_i <- metaprop(Paed_adm_prev_i$number_of_admissions_with_1_preventable_ae, Paed_adm_prev_i$sample_size, method.ci = "WS",
                               studlab = paste0(Paed_adm_prev_i$first_author_s_last_name, ", ",  Paed_adm_prev_i$publication_year,  "", Paed_adm_prev_i$Study_asterisk_prevAE),
                               prediction = T, title = "% of admissions preventable AEs")

summary(AEadm_prev_i)

rm(Paed_adm_prev_g, Paed_adm_prev_i, AEadm_prev_g, AEadm_prev_i)

# 8.3 Stratified by RRR Method---------------------------------------------------------------------

# 8.3.1 TT / GTT --------------------------
Paed_adm_prev_tt <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_preventable_ae)) %>%
  filter(rrr_method %in% "Trigger tool (TT)" | rrr_method %in% "Global Trigger tool (GTT)")
## 8 out of 33 studies have information on the outcome

AEadm_prev_tt_sub <- metaprop(Paed_adm_prev_tt$number_of_admissions_with_1_preventable_ae, Paed_adm_prev_tt$sample_size, method.ci = "WS",
                                    studlab = paste0(Paed_adm_prev_tt$first_author_s_last_name, ", ",  Paed_adm_prev_tt$publication_year,  "", Paed_adm_prev_tt$Study_asterisk_prevAE),
                                    prediction = T, title = "% of admissions preventable AEs", byvar = Paed_adm_prev_tt$Population)

summary(AEadm_prev_tt_sub)

forest.meta(AEadm_prev_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
                  sortvar = Paed_adm_prev_tt$sample_size,
                  leftlabs = c("TT/GTT Studies", "N of preventable", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)


# 8.3.1.1 Forest plot % of preventable AEs ---------------------------------------

# tiff("Plot6_tt_2022-09-15.tiff", width = 8, height =5, units = 'in', res = 200)
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_prev_tt_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm_prev_tt$sample_size,
#                   leftlabs = c("GTT/TT Studies", "N of admissions with", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1, overall.hetstat = F,
#                   prediction = F)
# grid::grid.text("preventable AEs",  .397, .8, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions with",  .844, .84, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("preventable AEs",  .844, .8, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .064, .56, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[0.0; 100.0]   ",  .949, .56, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .064, .2, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[2.5; 81.3]   ",  .949, .2, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()
rm(AEadm_prev_tt_sub, Paed_adm_prev_tt)

# 8.3.2 HMPS --------------------------
Paed_adm_prev_hmps <- Paed11 %>%
  filter(!is.na(number_of_admissions_with_1_preventable_ae)) %>%
  filter(rrr_method %in% "HMPS")
## 4 out of 33 studies have information on the outcome

AEadm_prev_hmps_sub <- metaprop(Paed_adm_prev_hmps$number_of_admissions_with_1_preventable_ae, Paed_adm_prev_hmps$sample_size, method.ci = "WS",
                                    studlab = paste0(Paed_adm_prev_hmps$first_author_s_last_name, ", ",  Paed_adm_prev_hmps$publication_year,  "", Paed_adm_prev_hmps$Study_asterisk_prevAE),
                                    prediction = T, title = "% of admissions preventable AEs", byvar = Paed_adm_prev_hmps$Population)

summary(AEadm_prev_hmps_sub)

forest.meta(AEadm_prev_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
                  print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
                  sortvar = Paed_adm_prev_hmps$sample_size,
                  leftlabs = c("HMPS Studies", "N of preventable", "Sample size"),
                  rightlabs = c("", "95% CI"),
                  print.subgroup.labels = T,
                  squaresize = 0.4, col.by = "black",
                  addrow = 1, overall.hetstat = F,
                  prediction = T,
                  addrow.subgroups = T, addrow.overall = T)


# 8.3.1.1 Forest plot % of preventable AEs ---------------------------------------

# tiff("Plot6_hmps_2022-09-15.tiff", width = 8, height =3, units = 'in', res = 200)
# options(na.action = "na.pass")
# 
# meta::forest.meta(AEadm_prev_hmps_sub, fontsize = 9, xlim = c(0,100), digits = 1, comb.fixed = F, print.I2 = F, print.I2.ci = F,
#                   print.pval.Q = F, print.tau2 = F, pscale = 100, smlab = "", test.subgroup = F, overall = F,
#                   sortvar = Paed_adm_prev_hmps$sample_size,
#                   leftlabs = c("HMPS Studies", "N of admissions with", "Sample size"),
#                   rightlabs = c("", "95% CI"),
#                   print.subgroup.labels = T,
#                   squaresize = 0.4, col.by = "black",
#                   addrow = 1, overall.hetstat = F,
#                   prediction = F)
# grid::grid.text("preventable AEs",  .397, .74, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("% of admissions with",  .844, .8, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("preventable AEs",  .844, .74, gp = grid::gpar(fontsize=9, font = 2), just="right")
# grid::grid.text("Prediction interval   ",  .075, .26, gp = grid::gpar(fontsize=9, font = 2), just="left")
# grid::grid.text("[0.0; 59.3]   ",  .938, .26, gp = grid::gpar(fontsize=9, font = 2), just="right")
# 
# dev.off()
rm(AEadm_prev_hmps_sub, Paed_adm_prev_hmps, Paed_adm_prev)
rm(AE100_meta, AE1000_meta, AEadm_prev, AEprev_sub)
