
#################################################################################################
## Measles CFR estimation -- 03_prediction.R
## Purpose: predict using either vaccination scenario
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

# date & description
model_description <- "baseline_testing"

vaccination_scenario <- 'baseline' # baseline, or no_vaccination

run_date_01b <- 'RUN_DATE'
run_date_01c <- 'RUN_DATE'
run_date_01d <- 'RUN_DATE'
run_date_02b <- 'RUN_DATE'
run_date_02c <- 'RUN_DATE'

run_date <- paste0(run_date_02c,  '_', vaccination_scenario, '_scenario_paper')

library(tidyverse, lib.loc = "FILEPATH")
### load MSCM packages -- must be done in this order
library(crosswalk002, lib.loc = "FILEPATH")
library(mrbrt002, lib.loc = "FILEPATH")

library(msm)
library(readxl)
library(data.table)
library(arm)
setDTthreads(1)

# -----------------------------------------------------------------------------------------------
# read covariates for prediction frame ----------------------------------------------------------
if(vaccination_scenario == "baseline"){
  covariates <- fread(paste0('FILEPATH',run_date_01c,'/transformed_standardized_covaraite_set_with_incidence.csv'))
}else if(vaccination_scenario == "no_vaccination"){
  covariates <- fread(paste0('FILEPATH',run_date_01c,'/transformed_standardized_covaraite_set_with_incidence_no_vax.csv'))
}

# -----------------------------------------------------------------------------------------------
# prepare prediction frame ----------------------------------------------------------------------
pred_frame_comm0 <- subset(covariates, year_id > 1989)
pred_frame_comm0$Comm.ind <- 0
pred_frame_comm1 <- subset(covariates, year_id > 1989)
pred_frame_comm1$Comm.ind <- 1

pred_frame <- rbind(pred_frame_comm0, pred_frame_comm1)

pred_frame <- pred_frame[complete.cases(pred_frame),]

pred_frame2 <- cbind(pred_frame, i = rep(0:99, each = nrow(pred_frame)))
pred_frame2$age_start <- pred_frame2$i
pred_frame2$age_end <- pred_frame2$i

# -----------------------------------------------------------------------------------------------
# make MR-BRT prediction object -----------------------------------------------------------------
## make frame to predict out on
data_pred <- MRData()
data_pred$load_df(data=pred_frame2,
                  col_covs=list( "incidence_standardized", 
                                 "war_rate_standardized",
                                 "maternal_education_standardized",
                                 "age_start",
                                 "age_end",
                                 "Comm.ind",
                                 "gdp_pc_standardized",
                                 "hiv_standardized", 
                                 "mcv1_standardized", 
                                 "tfr_standardized",
                                 "u5mr_standardized",
                                 "prop_urban_standardized",
                                 "vitA_standardized",
                                 "wasting_standardized"))

# -----------------------------------------------------------------------------------------------
# load MR-BRT model object and samples ----------------------------------------------------------
load(file = paste0('FILEPATH', run_date_02c,'/mrbrt_samples.RData'))
mod_cfr <- py_load_object(filename = paste0('FILEPATH',run_date_02c,'/mrbrt_mod_cfr_object.pkl'), pickle = "dill")

# -----------------------------------------------------------------------------------------------
# get point predictions from fixed effects ------------------------------------------------------
data_pred$predictions <- mod_cfr$predict(data=data_pred, predict_for_study = F)
pred_frame2$mod_cfr <- invlogit(data_pred$predictions)

# -----------------------------------------------------------------------------------------------
# create draws and uncertainty ------------------------------------------------------------------
draws <- mod_cfr$create_draws(
  data = data_pred,
  beta_samples = samples1[[1]],
  gamma_samples = samples1[[2]],
  random_study = FALSE
)

draws2 <- invlogit(draws)
draws2 <- data.table(draws2)

pred_frame2 <- as.data.table(pred_frame2)
pred_frame2[, `:=` (
  lower = apply(draws2, 1, function(x) quantile(x, 0.025)),
  median = apply(draws2, 1, function(x) quantile(x, 0.5)),
  upper = apply(draws2, 1, function(x) quantile(x, 0.975))
)]

draws2$ihme_loc_id <- pred_frame2$ihme_loc_id
draws2$year_id <- pred_frame2$year_id
draws2$super_region_name <- pred_frame2$super_region_name
draws2$age_start <- pred_frame2$age_start
draws2$Comm.ind <- pred_frame2$Comm.ind
draws2$mod_cfr <- pred_frame2$mod_cfr


pred_frame2$predicted_cfr <- pred_frame2$mod_cfr
pred_frame2$predicted_cfr_upper <- pred_frame2$upper
pred_frame2$predicted_cfr_lower <- pred_frame2$lower

# -----------------------------------------------------------------------------------------------
# save prediction file --------------------------------------------------------------------------
dir.create(paste0('FILEPATH',run_date), recursive = T)

save(pred_frame2, file = paste0('FILEPATH', run_date, '/pred_frame_with_uncertainty.csv'))
save(draws, file = paste0('FILEPATH', run_date, '/draws.RData'))
save(draws2, file = paste0('FILEPATH', run_date, '/draws2.RData'))

