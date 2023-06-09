
#################################################################################################
## Measles CFR estimation -- 01d_covariate_prep_incidence_covariate_no_vax_scenario.R
## Purpose: process incidence covariate like the others (transform, standardize, save) for the 
##          scenario with no vaccination
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")
run_date_01b <- 'RUN_DATE'
run_date_01c <- 'RUN_DATE'

library(data.table)

# -----------------------------------------------------------------------------------------------
# load covariate files --------------------------------------------------------------------------
covariates_to_include <- fread(paste0('FILEPATH',run_date_01b,'/transformed_standardized_covaraite_set.csv'))

# -----------------------------------------------------------------------------------------------
# add incidence to data -------------------------------------------------------------------------

# get no vaccination incidence
incidence <- fread('FILEPATH')
incidence <- subset(incidence, select=c("year", "country", "cohort_size", "cases"))
incidence$cohort_size <- incidence$cohort_size #* 1000
incidence$cohort_size <- as.numeric(incidence$cohort_size)
incidence$cases <- as.numeric(incidence$cases)

# still relying on cohort size from regular scenario
cohort_size <- fread('FILEPATH')
cohort_size <- subset(cohort_size, select=c("year", "country", "cohort_size"))
cohort_size$cohort_size <- as.numeric(cohort_size$cohort_size)

# add datasets together
incidence_all_age_cases <- aggregate(incidence$cases, by=list(incidence$year, incidence$country), FUN = sum)
colnames(incidence_all_age_cases) <- c("year", "country", "cases")
incidence_all_age_cohort <- aggregate(cohort_size$cohort_size, by=list(cohort_size$year, cohort_size$country), FUN = sum)
colnames(incidence_all_age_cohort) <- c("year", "country", "cohort_size")

incidence_all_age <- merge(incidence_all_age_cases, incidence_all_age_cohort, by=c("year", "country"))
incidence_all_age$incidence <- incidence_all_age$cases / incidence_all_age$cohort_size
quantile(incidence_all_age$incidence)

incidence_all_age <- subset(incidence_all_age, select=c("country", "year", "incidence"))

incidence_1981 <- subset(incidence_all_age, year == 1981)

incidence_1970 <- copy(incidence_1981)
incidence_1970$year <- 1970

incidence_1971 <- copy(incidence_1981)
incidence_1971$year <- 1971

incidence_1972 <- copy(incidence_1981)
incidence_1972$year <- 1972

incidence_1973 <- copy(incidence_1981)
incidence_1973$year <- 1973

incidence_1974 <- copy(incidence_1981)
incidence_1974$year <- 1974

incidence_1975 <- copy(incidence_1981)
incidence_1975$year <- 1975

incidence_1976 <- copy(incidence_1981)
incidence_1976$year <- 1976

incidence_1977 <- copy(incidence_1981)
incidence_1977$year <- 1977

incidence_1978 <- copy(incidence_1981)
incidence_1978$year <- 1978

incidence_1979 <- copy(incidence_1981)
incidence_1979$year <- 1979

incidence_1980 <- copy(incidence_1981)
incidence_1980$year <- 1980

incidence_early <- rbind(incidence_1980, incidence_1979, incidence_1978, incidence_1977, incidence_1976,
                         incidence_1975, incidence_1974, incidence_1973, incidence_1972, incidence_1971, incidence_1970 )

FINAL_incidence <- rbind(incidence_all_age, incidence_early)

covariates_to_include <- merge(covariates_to_include, FINAL_incidence, by.x=c("ihme_loc_id", "year_id"), by.y=c("country", "year"), all.x=T, all.y=F)

# -----------------------------------------------------------------------------------------------
# apply transformation --------------------------------------------------------------------------
covariates_to_include$incidence <- ifelse(covariates_to_include$incidence == 0, 0.0000000001, covariates_to_include$incidence)
covariates_to_include$logit_incidence <- logit(covariates_to_include$incidence)
covariates_to_include$incidence <- NULL

# -----------------------------------------------------------------------------------------------
# also, set coverage to 0 -----------------------------------------------------------------------
mcv_original_mean <- fread(paste0('FILEPATH',run_date_01b,'/mcv_mean.csv'))
mcv_original_sd <- fread(paste0('FILEPATH',run_date_01b,'/mcv_sd.csv'))

covariates_to_include$mcv_standardized <- (0 - as.numeric(mcv_original_mean)) / as.numeric(mcv_original_sd)

# -----------------------------------------------------------------------------------------------
# standardize covariates (subtract mean and divide by SD) ---------------------------------------
incidence_original_mean <- fread(paste0('FILEPATH',run_date_01c,'/incidence_mean.csv'))
incidence_original_sd <- fread(paste0('FILEPATH',run_date_01c,'/incidence_sd.csv'))

covariates_to_include <- as.data.table(covariates_to_include)
covariates_to_include$incidence_standardized <- (covariates_to_include$logit_incidence - as.numeric(incidence_original_mean)) / as.numeric(incidence_original_sd)
covariates_to_include$logit_incidence <- NULL

# -----------------------------------------------------------------------------------------------
# save ------------------------------------------------------------------------------------------
dir.create(paste0('FILEPATH',run_date))
fwrite(covariates_to_include, paste0('FILEPATH',run_date,'/transformed_standardized_covaraite_set_with_incidence_no_vax.csv'))


