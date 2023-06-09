
#################################################################################################
## Measles CFR estimation -- 02c_second_step_full_model.R
## Purpose: a) compute relative age patterns,
##          b) age split data using age specific incidence 
##          c) save input data to use in second stage model
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

# date & description
model_description <- "RUN_DATE"

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")
run_date <- paste0(run_date, '_', model_description)

run_date_01c <- 'RUN_DATE'
run_date_02a <- 'RUN_DATE'
run_date_02b <- 'RUN_DATE'

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
# read in data ----------------------------------------------------------------------------------
dt_collapsed <- fread(paste0('FILEPATH',run_date_02b, '/age_split_input_data.csv'))
dt_collapsed$Midpoint.Year <- floor((as.numeric(dt_collapsed$year_start) + as.numeric(dt_collapsed$year_end))/2)

dt_data <- subset(dt_collapsed, select=c("field_citation_value","ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
                                         "Comm.ind", "logit_cfr", "logit_cfr_se", "age_split"))
dt_data$ihme_loc_id <- gsub("_.*","",dt_data$ihme_loc_id)

# -----------------------------------------------------------------------------------------------
# read in covariates  ---------------------------------------------------------------------------
covariates <- fread(paste0('FILEPATH',run_date_01c,'/transformed_standardized_covaraite_set_with_incidence.csv'))

data_test <- merge(dt_data, covariates, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)
data_test <- subset(data_test, super_region_name != "High-income")

# -----------------------------------------------------------------------------------------------
# data prep  ------------------------------------------------------------------------------------
data_test$Comm.ind <- ifelse(is.na(data_test$Comm.ind), 0, data_test$Comm.ind)
data_test <- subset(data_test, Midpoint.Year > 1969)

missing <- data_test[!complete.cases(data_test),]

dt2 <- data_test[complete.cases(data_test),]
dt2$Country <- as.character(dt2$ihme_loc_id)
dt2$for_re <- as.character(dt2$field_citation_value)

dt2$age_end <- as.numeric(dt2$age_end)
dt2$age_start <- as.numeric(dt2$age_start)

# -----------------------------------------------------------------------------------------------
# run MR-BRT model  ------------------------------------------------------------------------------------

dat1 <- MRData()
dat1$load_df(data = dt2,  col_obs = "logit_cfr", col_obs_se = "logit_cfr_se", col_study_id = "for_re",
             col_covs = list( "incidence_standardized", 
                              "age_start",
                              "age_end",
                              "war_rate_standardized",
                              "maternal_education_standardized",
                              "Comm.ind",
                              "gdp_pc_standardized",
                              "hiv_standardized", 
                              "mcv1_standardized", 
                              "tfr_standardized",
                              "u5mr_standardized",
                              "prop_urban_standardized",
                              "vitA_standardized",
                              "wasting_standardized"))

set.seed(1234)

old_knots <- as.vector(fread(paste0('FILEPATH',run_date_02a,'/knot_locations.csv')))

mod_cfr <- MRBRT(data = dat1,
                 inlier_pct=1,
                 cov_models = list(
                   LinearCovModel("intercept", use_re = FALSE),
                   LinearCovModel("Comm.ind",  prior_beta_uniform=array(c(-Inf, 0.0))),
                   LinearCovModel(
                     alt_cov = c("age_start", "age_end"),
                     use_spline=TRUE, 
                     use_spline_intercept=FALSE, 
                     spline_degree = 2L,
                     spline_knots = array(c(0, as.numeric(old_knots[2])/99, as.numeric(old_knots[3])/99, as.numeric(old_knots[4])/99, 1)),  
                     spline_knots_type = "domain",
                     spline_l_linear = FALSE,
                     spline_r_linear = TRUE 
                   ),
                   LinearCovModel("incidence_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("war_rate_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("maternal_education_standardized", prior_beta_uniform=array(c(-Inf, 0.0))),
                   LinearCovModel("gdp_pc_standardized", prior_beta_uniform=array(c(-Inf, 0.0))),
                   LinearCovModel("hiv_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("mcv1_standardized", prior_beta_uniform=array(c(-Inf, 0.0))),
                   LinearCovModel("tfr_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("u5mr_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("prop_urban_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("vitA_standardized", prior_beta_uniform=array(c(0.0, Inf))),
                   LinearCovModel("wasting_standardized", prior_beta_uniform=array(c(0.0, Inf))))
) 

mod_cfr$fit_model(inner_print_level = 5L, inner_max_iter = 5000L)

# get samples
n_samples <- 1000L
samples1 <- mod_cfr$sample_soln(sample_size = n_samples)

dir.create(paste0('FILEPATH',run_date,'/'), recursive = T)
save(samples1, file = paste0('FILEPATH',run_date,'/mrbrt_samples.RData'))
py_save_object(object = mod_cfr, filename = paste0('FILEPATH',run_date,'/mrbrt_mod_cfr_object.pkl'), pickle = "dill")

