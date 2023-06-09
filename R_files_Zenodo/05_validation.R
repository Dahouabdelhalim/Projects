#################################################################################################
## Measles CFR estimation -- 05_validation.R
## Purpose: a) run IS and OOS validation models, and
##.         b) compute validation metrics
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

# date & description
model_description <- "RUN_DATE"

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")

run_date_01c <- 'RUN_DATE'
run_date_02a <- 'RUN_DATE'
run_date_02b <- 'RUN_DATE'
run_date_02c <- 'RUN_DATE'

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
dt <- read_excel('FILEPATH')
dt <- as.data.table(dt)
dt <- dt[-1,]
dt <- subset(dt, !is.na(extractor))

# fix columns
dt_collapsed <- unique(copy(dt))
dt_collapsed$Deaths <- as.numeric(dt_collapsed$cases)
dt_collapsed$Cases <- as.numeric(dt_collapsed$sample_size)
dt_collapsed$age_end <- as.numeric(dt_collapsed$age_end)
dt_collapsed$age_start <- as.numeric(dt_collapsed$age_start)
dt_collapsed <- subset(dt_collapsed, group_review == 1)
dt_collapsed$ihme_loc_id <- gsub("_.*","",dt_collapsed$ihme_loc_id)

# add age information
dt_collapsed$age_difference <- dt_collapsed$age_end - dt_collapsed$age_start
# dt_collapsed <- subset(dt_collapsed, age_difference <=5)
dt_collapsed <- data.table(dt_collapsed)


# compute logit(cfr) using delta transform
dt_collapsed[, ratio := Deaths/Cases]
dt_collapsed[ratio==0, ratio := 0.0004404757/2]
dt_collapsed[ratio==1, ratio := 0.9999999999999]
dt_collapsed <- subset(dt_collapsed, ratio  <= 1)
dt_collapsed <- data.table(dt_collapsed)
z <- qnorm(0.975)
dt_collapsed[, ratio_se := sqrt((Deaths/Cases)*(1-(Deaths/Cases))/Cases + z^2/(4*Cases^2))] 
dt_collapsed <- as.data.frame(dt_collapsed)
dt_collapsed$logit_cfr <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,1]
dt_collapsed$logit_cfr_se <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,2]
dt_collapsed <- as.data.table(dt_collapsed)

# add final columns and clean up
dt_collapsed$Age.ind <- ifelse(dt_collapsed$age_start < 5, 1, 0)
dt_collapsed$Comm.ind <- ifelse(dt_collapsed$CV_Hospital == 1, 0, 1)

dt_collapsed$Midpoint.Year <- floor((as.numeric(dt_collapsed$year_start) + as.numeric(dt_collapsed$year_end))/2)

dt_data <- subset(dt_collapsed, select=c("field_citation_value","ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
                                         "Age.ind", "Comm.ind", "logit_cfr", "logit_cfr_se", "ratio_se"))


# -----------------------------------------------------------------------------------------------
# read in covariates  ---------------------------------------------------------------------------
covariates <- fread(paste0('FILEPATH',run_date_01c,'/transformed_standardized_covaraite_set_with_incidence.csv'))

data_test <- merge(dt_data, covariates, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)
data_test <- subset(data_test, super_region_name != "High-income")
data_test <- data_test[complete.cases(data_test),]

# -----------------------------------------------------------------------------------------------
# first, let's do OOS  --------------------------------------------------------------------------

pred_frame <- copy(data_test)
pred_frame$logit_cfr <- NULL
pred_frame$logit_cfr_se <- NULL

pred_frame_list <- list()

pred_frame <- pred_frame[complete.cases(pred_frame),]

age_gran_data <- subset(data_test, age_end - age_start <= 5)
age_wide <- subset(data_test, age_end - age_start > 5)

granular_nids <- unique(age_gran_data$field_citation_value)
nongranular_nids <- setdiff(unique(data_test$field_citation_value), granular_nids)


n_folds <- 5

granular_nids <- data.table(granular_nids)
granular_nids$folds <- sample(cut(seq(1, nrow(granular_nids)), breaks = n_folds, labels = 1:n_folds))
colnames(granular_nids) <- c("source", "fold")
nongranular_nids <- data.table(nongranular_nids)
nongranular_nids$folds <- sample(cut(seq(1, nrow(nongranular_nids)), breaks = n_folds, labels = 1:n_folds))
colnames(nongranular_nids) <- c("source", "fold")

source_folds <- rbind(granular_nids, nongranular_nids)

data <- merge(data_test, source_folds, by.x="field_citation_value", by.y="source")

# -----------------------------------------------------------------------------------------------
# first stage model -----------------------------------------------------------------------------

for(i in 1:n_folds){
  
  message("now on fold... ", i)
  data_ho <- subset(data, fold !=i)
  
  dt2 <- data_ho[complete.cases(data_ho),]
  dt2$Country <- as.character(dt2$Country)
  dt2$for_re <- as.character(dt2$Country)
  
  dt2 <- subset(dt2, age_end - age_start <=5)
  
  dat1 <- MRData()
  dat1$load_df(data = dt2,  col_obs = "logit_cfr", col_obs_se = "logit_cfr_se", col_study_id = "for_re",
               col_covs = list( "incidence_standardized", 
                                "age_start",
                                "age_end" ,
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
  mod_cfr <- MRBRT(data = dat1,
                   inlier_pct=1,
                   cov_models = list(
                     LinearCovModel("intercept", use_re = FALSE),
                     LinearCovModel("Comm.ind"),
                     LinearCovModel("incidence_standardized") ,
                     LinearCovModel("war_rate_standardized"),
                     LinearCovModel("maternal_education_standardized"),
                     LinearCovModel("gdp_pc_standardized"),
                     LinearCovModel("hiv_standardized"),
                     LinearCovModel("mcv1_standardized"),
                     LinearCovModel("tfr_standardized"),
                     LinearCovModel("u5mr_standardized"),
                     LinearCovModel("prop_urban_standardized"),
                     LinearCovModel("vitA_standardized"),
                     LinearCovModel("wasting_standardized"),
                     LinearCovModel(
                       alt_cov = c("age_start", "age_end"),
                       use_spline=TRUE, 
                       use_spline_intercept=FALSE, # set this to F if using an intercept 
                       spline_degree = 2L,
                       spline_knots = array(seq(0, 1, length.out = 5)),  # 3 internal  knots, for testing!!
                       spline_knots_type = "frequency",
                       spline_l_linear = FALSE,
                       spline_r_linear = TRUE)))
  
  mod_cfr$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)
  
  # -----------------------------------------------------------------------------------------------
  # save knot locations --------------------------------------------------------------
  get_knots <- function(model) {
    model$cov_models[[which(model$cov_model_names == tail(model$cov_model_names)[6])]]$spline$knots
  }
  
  model_knots <- as.data.frame(get_knots(mod_cfr))
  dir.create(paste0('FILEPATH',run_date,'/'))
  fwrite(model_knots, paste0('FILEPATH',run_date,'/holdout_',i,'knot_locations.csv'))
  
  # -----------------------------------------------------------------------------------------------
  # set up prediction frame -----------------------------------------------------------------------
  pred_frame <- subset(pred_frame, select=c("ihme_loc_id", "Midpoint.Year",
                                            "Comm.ind",
                                            "incidence_standardized",
                                            "war_rate_standardized",
                                            "maternal_education_standardized",
                                            "gdp_pc_standardized",
                                            "hiv_standardized",
                                            "mcv1_standardized",
                                            "tfr_standardized",
                                            "u5mr_standardized",
                                            "prop_urban_standardized",
                                            "vitA_standardized",
                                            "wasting_standardized"))
  pred_frame <- unique(pred_frame)
  pred_frame <- pred_frame[complete.cases(pred_frame),]
  new_pred_frame <- pred_frame[1,]
  
  new_pred_frame$incidence_standardized <- mean(pred_frame$incidence_standardized)
  new_pred_frame$war_rate_standardized <- mean(pred_frame$war_rate_standardized)
  new_pred_frame$maternal_education_standardized <- mean(pred_frame$maternal_education_standardized)
  new_pred_frame$gdp_pc_standardized <- mean(pred_frame$gdp_pc_standardized)
  new_pred_frame$hiv_standardized <- mean(pred_frame$hiv_standardized)
  new_pred_frame$mcv1_standardized <- mean(pred_frame$mcv1_standardized)
  new_pred_frame$tfr_standardized <- mean(pred_frame$tfr_standardized)
  new_pred_frame$u5mr_standardized <- mean(pred_frame$u5mr_standardized)
  new_pred_frame$prop_urban_standardized <- mean(pred_frame$prop_urban_standardized)
  new_pred_frame$vitA_standardized <- mean(pred_frame$vitA_standardized)
  new_pred_frame$wasting_standardized <- mean(pred_frame$wasting_standardized)
  
  new_pred_frame$year_id <- NULL
  new_pred_frame$location_name <- NULL
  new_pred_frame$super_region_name <- NULL
  
  pred_frame <- cbind(new_pred_frame, age_start = rep(0:34, each = nrow(new_pred_frame)))
  pred_frame$age_end <- pred_frame$age_start
  
  
  # -----------------------------------------------------------------------------------------------
  # prediction ------------------------------------------------------------------------------------
  data_pred <- MRData()
  data_pred$load_df(data=pred_frame,
                    col_covs=list( "age_start",
                                   "age_end",
                                   "Comm.ind",
                                   "incidence_standardized",
                                   "war_rate_standardized",
                                   "maternal_education_standardized",
                                   "gdp_pc_standardized",
                                   "hiv_standardized",
                                   "mcv1_standardized",
                                   "tfr_standardized",
                                   "u5mr_standardized",
                                   "prop_urban_standardized",
                                   "vitA_standardized",
                                   "wasting_standardized"))
  
  # get point predictions using only fixed effects
  data_pred$predictions <- mod_cfr$predict(data=data_pred, predict_for_study = F)
  pred_frame$mod_cfr <- invlogit(data_pred$predictions)
  save(pred_frame, file= paste0('FILEPATH',run_date,'/holdout_', i,'IS_pred_frame_age_curve.RData'))
  
}

# -----------------------------------------------------------------------------------------------
# now to age split the data ---------------------------------------------------------------------

for(FOLD_i in 1:n_folds){
message("now on fold ", FOLD_i)
  # -----------------------------------------------------------------------------------------------
  # read in data ----------------------------------------------------------------------------------
  dt <- read_excel('FILEPATH')
  dt <- as.data.table(dt)
  dt <- dt[-1,]
  dt <- subset(dt, !is.na(extractor))
  
  # fix column names
  dt_collapsed <- unique(copy(dt))
  dt_collapsed$Deaths <- as.numeric(dt_collapsed$cases)
  dt_collapsed$Cases <- as.numeric(dt_collapsed$sample_size)
  dt_collapsed$age_end <- as.numeric(dt_collapsed$age_end)
  dt_collapsed$age_start <- as.numeric(dt_collapsed$age_start)
  
  dt_collapsed <- subset(dt_collapsed, group_review == 1)
  
  # -----------------------------------------------------------------------------------------------
  # load age specific incidence and first stage results -------------------------------------------
  age_incidence <- fread('FILEPATH')
  load(paste0('FILEPATH',run_date,'/holdout_',i,'IS_pred_frame_age_curve.RData'))
  age_cfr <- pred_frame
  pred_frame <- NULL
  age_cfr$reference_proportion <- age_cfr$mod_cfr / age_cfr$mod_cfr[1]
  
  age_cfr <- subset(age_cfr, select=c("age_start", "age_end", "mod_cfr", "reference_proportion"))
  
  
  # -----------------------------------------------------------------------------------------------
  # prepare data for splitting --------------------------------------------------------------------
  
  dt_collapsed$age_difference <- dt_collapsed$age_end - dt_collapsed$age_start
  dt_collapsed$ihme_loc_id <- gsub("_.*","",dt_collapsed$ihme_loc_id)
  dt_collapsed$Comm.ind <- ifelse(dt_collapsed$CV_Hospital == 1, 0, 1)
  
  dt_collapsed <- subset(dt_collapsed, select=c("nid","field_citation_value","ihme_loc_id", "year_start", "year_end", "Cases", "Deaths", "age_start", "age_end",
                                                "age_difference", "Comm.ind", "age_demographer" ))
  
  dt_collapsed <- data.table(dt_collapsed)
  dt_collapsed$Cases <- as.numeric(dt_collapsed$Cases)
  dt_collapsed$Deaths <- as.numeric(dt_collapsed$Deaths)
  
  case_df <- aggregate(dt_collapsed$Cases, by=list(dt_collapsed$nid, dt_collapsed$field_citation_value, dt_collapsed$ihme_loc_id,
                                                   dt_collapsed$year_start, dt_collapsed$year_end, 
                                                   dt_collapsed$age_start, dt_collapsed$age_end, dt_collapsed$age_demographer, 
                                                   dt_collapsed$Comm.ind, dt_collapsed$age_difference), FUN="sum")
  colnames(case_df) <- c("nid", "field_citation_value", "ihme_loc_id","year_start", "year_end", 
                         "age_start", "age_end", "age_demographer", "Comm.ind", "age_difference", "Cases")
  death_df <- aggregate(dt_collapsed$Deaths, by=list(dt_collapsed$nid, dt_collapsed$field_citation_value, dt_collapsed$ihme_loc_id,
                                                     dt_collapsed$year_start, dt_collapsed$year_end, 
                                                     dt_collapsed$age_start, dt_collapsed$age_end, dt_collapsed$age_demographer, 
                                                     dt_collapsed$Comm.ind, dt_collapsed$age_difference), FUN="sum")
  colnames(death_df) <- c("nid", "field_citation_value", "ihme_loc_id","year_start", "year_end", 
                          "age_start", "age_end", "age_demographer", "Comm.ind", "age_difference", "Deaths")
  
  dt_colllapsed2 <- merge(case_df, death_df, by=c("nid", "field_citation_value", "ihme_loc_id","year_start", "year_end", 
                                                  "age_start", "age_end", "age_demographer", "Comm.ind", "age_difference"))
  
  # -----------------------------------------------------------------------------------------------
  # age splitting ---------------------------------------------------------------------------------
  
  dt_colllapsed2 <- merge(dt_colllapsed2, source_folds, by.x="field_citation_value", by.y="source")
  dt_colllapsed2 <- subset(dt_colllapsed2, fold != FOLD_i)
    
    
  original_dt_collapsed <- copy(dt_colllapsed2)
  
  dt_collapsed_small_bins <- subset(dt_colllapsed2, age_difference < 1 | (age_difference ==1 & age_demographer == 0))
  dt_collapsed_to_split <- subset(dt_colllapsed2, age_difference >1 | (age_difference ==1 & age_demographer == 1))
  
  finished <- data.table()
  
  for(f in 1:length(unique(dt_collapsed_to_split$field_citation_value))){
    
    message(f/length(unique(dt_collapsed_to_split$field_citation_value)))
    
    field_citation <- unique(dt_collapsed_to_split$field_citation_value)[f]
    
    test <- subset(dt_collapsed_to_split, field_citation_value == field_citation)
    
    for(i in 1:dim(test)[1]){
      
      test_row <- test[i,]
      
      loc <- unique(test$ihme_loc_id[i])
      year_start <- mean(as.numeric(test$year_start[i]) )
      year_end <- mean(as.numeric(test$year_end[i]) )
      mid.year <- floor((year_start + year_end)/2)
      
      mid.year <- ifelse(mid.year < 1981, 1981, mid.year)
      
      study_incidence <- subset(age_incidence, country == loc & year == mid.year & age <= test_row$age_end & age >= test_row$age_start)
      
      if(dim(study_incidence)[1] > 0){
        age_cases <- sum(study_incidence$cases)
        study_incidence$age_proportion <- study_incidence$cases / age_cases
        
        test_row_incidence <- merge(test_row, study_incidence, by.y='country', by.x='ihme_loc_id')
        test_row_incidence$age_start <- NULL
        test_row_incidence$age_end <- NULL
        test_row_incidence$age_start <- test_row_incidence$age
        test_row_incidence$age_end <- test_row_incidence$age
        
        test_row_incidence$Cases <- test_row_incidence$Cases * test_row_incidence$age_proportion
        
        total_deaths <- unique(test_row_incidence$Deaths)
        test_row_incidence <- merge(test_row_incidence, age_cfr, by=c('age_start', 'age_end'), all.x=T, all.y=F)
        
        plus35 <- subset(age_cfr, age_start == 34)$reference_proportion
        test_row_incidence$reference_proportion <- ifelse(test_row_incidence$age_start > 34, plus35, test_row_incidence$reference_proportion)
        
        x_reference <- total_deaths / sum(test_row_incidence$Cases * test_row_incidence$reference_proportion) 
        test_row_incidence$calc_cfr <-  x_reference * test_row_incidence$reference_proportion
        test_row_incidence$Deaths <-test_row_incidence$calc_cfr * test_row_incidence$Cases
        test_row_incidence$age_split <- 1
      }else{
        test_row_incidence <- test_row
        test_row_incidence$age_split <- 0
      }
      finished <- rbind(finished, test_row_incidence, fill=T)
    } 
  }
  
  dt_collapsed_small_bins <- data.table(dt_collapsed_small_bins)
  finished <- data.table(finished)
  
  dt_collapsed_small_bins$age_split <- 2
  dt_collapsed3 <- rbind(dt_collapsed_small_bins, finished, fill=T)
  dt_collapsed3 <- data.table(dt_collapsed3)
  
  # -----------------------------------------------------------------------------------------------
  # re-prep data in logit space -------------------------------------------------------------------
  
  dt_collapsed3[, ratio := Deaths/Cases]
  dt_collapsed3[ratio==0, ratio := 0.0004404757/2]
  dt_collapsed3[ratio==1, ratio := 0.9999999999999]
  dt_collapsed3 <- subset(dt_collapsed3, ratio  <= 1)
  
  dt_collapsed3 <- data.table(dt_collapsed3)
  z <- qnorm(0.975)
  dt_collapsed3[, ratio_se := sqrt((Deaths/Cases)*(1-(Deaths/Cases))/Cases + z^2/(4*Cases^2))] 
  dt_collapsed3 <- as.data.frame(dt_collapsed3)
  
  dt_collapsed3$logit_cfr <- delta_transform(mean = dt_collapsed3$ratio, sd = dt_collapsed3$ratio_se, transformation = "linear_to_logit")[,1]
  dt_collapsed3$logit_cfr_se <- delta_transform(mean = dt_collapsed3$ratio, sd = dt_collapsed3$ratio_se, transformation = "linear_to_logit")[,2]
  dt_collapsed3$logit_upper_cfr <- dt_collapsed3$logit_cfr + (1.96*dt_collapsed3$logit_cfr_se)
  dt_collapsed3$logit_lower_cfr <- dt_collapsed3$logit_cfr - (1.96*dt_collapsed3$logit_cfr_se)
  dt_collapsed3$upper_cfr  <- invlogit(dt_collapsed3$logit_upper_cfr)
  dt_collapsed3$lower_cfr  <- invlogit(dt_collapsed3$logit_lower_cfr)
 
  dt_collapsed <- copy(dt_collapsed3)
  
  dir.create(paste0('FILEPATH',run_date))
  fwrite(dt_collapsed, paste0('FILEPATH',run_date, '/holdout_',FOLD_i,'age_split_input_data.csv'))
  
  message("now FINISHED fold ", FOLD_i)
  
}




for(i in 1:n_folds){
  
  message("now on fold... ", i)

  # -----------------------------------------------------------------------------------------------
  # read in data ----------------------------------------------------------------------------------
  dt_collapsed <- fread(paste0('FILEPATH',run_date, '/holdout_',i,'age_split_input_data.csv'))
  dt_collapsed$Midpoint.Year <- floor((as.numeric(dt_collapsed$year_start) + as.numeric(dt_collapsed$year_end))/2)
  
  dt_data <- subset(dt_collapsed, select=c("ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
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
  dt2$for_re <- as.character(dt2$ihme_loc_id)
  
  dt2$age_end <- as.numeric(dt2$age_end)
  dt2$age_start <- as.numeric(dt2$age_start)
  
  
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
  
  old_knots <- as.vector(fread(paste0('FILEPATH',run_date,'/holdout_',i,'knot_locations.csv')))
  mod_cfr <- MRBRT(data = dat1,
                   inlier_pct=1,
                   cov_models = list(
                     LinearCovModel("intercept", use_re = FALSE),
                     LinearCovModel("Comm.ind"),
                     LinearCovModel(
                       alt_cov = c("age_start", "age_end"),
                       use_spline=TRUE, 
                       use_spline_intercept=FALSE, # set this to F if using an intercept 
                       spline_degree = 2L,
                       spline_knots = array(c(0, as.numeric(old_knots[2])/99, as.numeric(old_knots[3])/99, as.numeric(old_knots[4])/99, 1)),  # 3  knots
                       spline_knots_type = "domain",
                       spline_l_linear = FALSE,
                       spline_r_linear = TRUE 
                     ),
                     LinearCovModel("incidence_standardized"),
                     LinearCovModel("war_rate_standardized"),
                     LinearCovModel("maternal_education_standardized"),
                     LinearCovModel("gdp_pc_standardized"),
                     LinearCovModel("hiv_standardized"),
                     LinearCovModel("mcv1_standardized"),
                     LinearCovModel("tfr_standardized"),
                     LinearCovModel("u5mr_standardized"),
                     LinearCovModel("prop_urban_standardized"),
                     LinearCovModel("vitA_standardized"),
                     LinearCovModel("wasting_standardized"))
  )
  
  mod_cfr$fit_model(inner_print_level = 5L, inner_max_iter = 5000L)
  
  ## prediction
  pred_frame_comm0 <- subset(covariates, year_id > 1989)
  pred_frame_comm0$Comm.ind <- 0
  pred_frame_comm1 <- subset(covariates, year_id > 1989)
  pred_frame_comm1$Comm.ind <- 1
  
  pred_frame <- rbind(pred_frame_comm0, pred_frame_comm1)
  
  pred_frame <- pred_frame[complete.cases(pred_frame),]
  
  pred_frame2 <- cbind(pred_frame, i = rep(0:99, each = nrow(pred_frame)))
  pred_frame2$age_start <- pred_frame2$i
  pred_frame2$age_end <- pred_frame2$i
  
  
  
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
  
  # get point predictions using only fixed effects
  data_pred$predictions <- mod_cfr$predict(data=data_pred, predict_for_study = F)
  pred_frame2$cfr_predictions <- invlogit(data_pred$predictions)
  
  pred_frame_list[[i]] <- pred_frame2
  
}

dir.create(paste0('FILEPATH',run_date,'/'), recursive = T)
save(pred_frame_list, file = paste0('FILEPATH',run_date,'/oos_pred_frame_list.RData'))

## compute validation metrics
ho_metrics <- data.table()
data_match_all <- copy(dt_collapsed)

for(i in 1:n_folds){
  message(i)
  data_match <- merge(data_match_all, pred_frame_list[[i]], by.x=c("ihme_loc_id", "Midpoint.Year", "Comm.ind", "age_start", "age_end"), by.y=c("ihme_loc_id", "year_id", "Comm.ind", "age_start", "age_end"))
  
  fold <- i
  cor.val <- cor(data_match$cfr_predictions, data_match$ratio)
  rmse.val <- sqrt(mean((data_match$ratio - data_match$cfr_predictions)^2)) 
  me.val <- mean(data_match$ratio - data_match$cfr_predictions)
  mae.val <- mean(abs(data_match$ratio - data_match$cfr_predictions))
  
  
  to_add <- data.table(cbind(fold, cor.val, rmse.val, me.val, mae.val))
  
  ho_metrics <- rbind(ho_metrics, to_add)
}

summary_ho_metrics <- colMeans(ho_metrics)
summary_ho_metrics$fold <- NULL


fwrite(summary_ho_metrics, paste0('FILEPATH',run_date,'/oos_summary_holdout_metrics_agesplit.csv'))

# -----------------------------------------------------------
# now, let's do IS  --------------------------------------------------------------------------

pred_frame <- copy(data_test)
pred_frame$logit_cfr <- NULL
pred_frame$logit_cfr_se <- NULL

pred_frame <- pred_frame[complete.cases(pred_frame),]

data_pred <- MRData()
data_pred$load_df(data=pred_frame,
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

## load MR-BRT model object and samples 
mod_cfr <- py_load_object(filename = paste0('FILEPATH',run_date_02c,'/mrbrt_mod_cfr_object.pkl'), pickle = "dill")

# -----------------------------------------------------------------------------------------------
# get point predictions from fixed effects ------------------------------------------------------
data_pred$predictions <- mod_cfr$predict(data=data_pred, predict_for_study = F)
pred_frame$cfr_predictions <- invlogit(data_pred$predictions)

save(pred_frame, file = paste0('FILEPATH',run_date,'/is_pred_frame.RData'))

## compute validation metrics
ho_metrics <- data.table()
data_match_all <- copy(dt_collapsed)

data_match <- merge(data_match_all, pred_frame, by=c("ihme_loc_id", "Midpoint.Year", "Comm.ind", "age_start", "age_end"))

cor.val <- cor(data_match$cfr_predictions, data_match$ratio)
rmse.val <- sqrt(mean((data_match$ratio - data_match$cfr_predictions)^2)) 
me.val <- mean(data_match$ratio - data_match$cfr_predictions)
mae.val <- mean(abs(data_match$ratio - data_match$cfr_predictions))


to_add <- data.table(cbind(cor.val, rmse.val, me.val, mae.val))

ho_metrics <- rbind(ho_metrics, to_add)


fwrite(ho_metrics, paste0('FILEPATH',run_date,'/is_summary_holdout_metrics_agesplit.csv'))


###################################################################################################
###################################################################################################
#### compare to original data

## compute validation metrics
ho_metrics <- data.table()
# data_match_all <- copy(dt_collapsed)

for(i_FOLD in 1:n_folds){
  
  # -----------------------------------------------------------------------------------------------
  # read in data ----------------------------------------------------------------------------------
  dt <- read_excel('FILEPATH')
  dt <- as.data.table(dt)
  dt <- dt[-1,]
  dt <- subset(dt, !is.na(extractor))
  
  # fix columns
  dt_collapsed <- unique(copy(dt))
  dt_collapsed$Deaths <- as.numeric(dt_collapsed$cases)
  dt_collapsed$Cases <- as.numeric(dt_collapsed$sample_size)
  dt_collapsed$age_end <- as.numeric(dt_collapsed$age_end)
  dt_collapsed$age_start <- as.numeric(dt_collapsed$age_start)
  dt_collapsed <- subset(dt_collapsed, group_review == 1)
  dt_collapsed$ihme_loc_id <- gsub("_.*","",dt_collapsed$ihme_loc_id)
  
  # add age information
  dt_collapsed$age_difference <- dt_collapsed$age_end - dt_collapsed$age_start
  # dt_collapsed <- subset(dt_collapsed, age_difference <=5)
  dt_collapsed <- data.table(dt_collapsed)
  
  # # subset to ONLY OLD studies
  dt_collapsed <- subset(dt_collapsed, `in portnoy` == 0)
  # 
  # compute logit(cfr) using delta transform
  dt_collapsed[, ratio := Deaths/Cases]
  dt_collapsed[ratio==0, ratio := 0.0004404757/2]
  dt_collapsed[ratio==1, ratio := 0.9999999999999]
  dt_collapsed <- subset(dt_collapsed, ratio  <= 1)
  dt_collapsed <- data.table(dt_collapsed)
  z <- qnorm(0.975)
  dt_collapsed[, ratio_se := sqrt((Deaths/Cases)*(1-(Deaths/Cases))/Cases + z^2/(4*Cases^2))] 
  dt_collapsed <- as.data.frame(dt_collapsed)
  dt_collapsed$logit_cfr <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,1]
  dt_collapsed$logit_cfr_se <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,2]
  dt_collapsed <- as.data.table(dt_collapsed)
  
  # add final columns and clean up
  dt_collapsed$Age.ind <- ifelse(dt_collapsed$age_start < 5, 1, 0)
  dt_collapsed$Comm.ind <- ifelse(dt_collapsed$CV_Hospital == 1, 0, 1)
  
  dt_collapsed$Midpoint.Year <- floor((as.numeric(dt_collapsed$year_start) + as.numeric(dt_collapsed$year_end))/2)
  
  dt_data <- subset(dt_collapsed, select=c("ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
                                           "Age.ind", "Comm.ind", "logit_cfr", "logit_cfr_se", "ratio_se"))
  
  need_averages <- subset(dt_data, select=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end", "Comm.ind"))
  need_averages <- unique(need_averages)
  need_averages  <- data.table(need_averages)
  
  pred_frame2 <- data.table(pred_frame_list[[i_FOLD]])
  
  age_incidence <- fread('FILEPATH')
  
  need_averages$cfr <- 0
  for(i in 1:dim(need_averages)[1]){
    message(i / dim(need_averages)[1])
    
    loc <- as.character(need_averages[i,1] )
    yr <- as.numeric(need_averages[i,2])
    community <- as.numeric(need_averages[i,5])
    age_min <- as.numeric(need_averages[i,3])
    if(age_min < 1){
      age_min = 0
    }
    age_max <- as.numeric(need_averages[i,4])
    if(age_max < 1){
      age_max = 0
    }
    
    p_df <- subset(pred_frame2, ihme_loc_id == loc  & year_id == yr & Comm.ind == community & age_start >= age_min & age_end <= age_max)
    age_df <- subset(age_incidence, country == loc  & year == yr & age >= age_min & age <= age_max)
    age_df <- subset(age_df, select=c('age', 'cases'))  
    
    p_df2 <- merge(p_df, age_df, by.x="age_start", by.y= "age")  
    
    p_df2$deaths <- p_df2$cfr_predictions * p_df2$cases
    
    case_total <- sum( p_df2$cases)
    death_total <- sum( p_df2$deaths)
    
    
    overall_cfr <- death_total / case_total
    
    
    need_averages$cfr[i] <- overall_cfr
    
  }
  
  
  need_averages <- subset(need_averages, cfr != "NaN")
  
  
  
  dt_data$ihme_loc_id <- gsub("_.*","",dt_data$ihme_loc_id)
  
  dt_data <- subset(dt_data, age_start < 35)
  
  
  dt_data <- merge(dt_data, need_averages, by.x=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end","Comm.ind"), by.y=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end","Comm.ind"))
  
  
  dt_data$ratio <- invlogit(dt_data$logit_cfr)
  
  message("fold is ",i_FOLD)
  data_match <- dt_data 
  fold <- i_FOLD
  cor.val <- cor(data_match$cfr, data_match$ratio)
  rmse.val <- sqrt(mean((data_match$ratio - data_match$cfr)^2)) 
  me.val <- mean(data_match$ratio - data_match$cfr)
  mae.val <- mean(abs(data_match$ratio - data_match$cfr))
  
  to_add <- data.table(cbind(fold, cor.val, rmse.val, me.val, mae.val))
  
  ho_metrics <- rbind(ho_metrics, to_add)
}

summary_ho_metrics <- colMeans(ho_metrics)
summary_ho_metrics$fold <- NULL

fwrite(summary_ho_metrics, paste0('FILEPATH',run_date,'/oos_summary_holdout_metrics_originaldata.csv'))

# -----------------------------------------------------------
# now, let's do IS  --------------------------------------------------------------------------

data_pred <- MRData()
data_pred$load_df(data=pred_frame,
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

## load MR-BRT model object and samples 
mod_cfr <- py_load_object(filename = paste0('FILEPATH',run_date_02c,'/mrbrt_mod_cfr_object.pkl'), pickle = "dill")

# -----------------------------------------------------------------------------------------------
# get point predictions from fixed effects ------------------------------------------------------
data_pred$predictions <- mod_cfr$predict(data=data_pred, predict_for_study = F)
pred_frame$cfr_predictions <- invlogit(data_pred$predictions)

save(pred_frame, file = paste0('FILEPATH',run_date,'/is_pred_frame.RData'))


# -----------------------------------------------------------------------------------------------
# read in data ----------------------------------------------------------------------------------
dt <- read_excel('FILEPATH')
dt <- as.data.table(dt)
dt <- dt[-1,]
dt <- subset(dt, !is.na(extractor))

# fix columns
dt_collapsed <- unique(copy(dt))
dt_collapsed$Deaths <- as.numeric(dt_collapsed$cases)
dt_collapsed$Cases <- as.numeric(dt_collapsed$sample_size)
dt_collapsed$age_end <- as.numeric(dt_collapsed$age_end)
dt_collapsed$age_start <- as.numeric(dt_collapsed$age_start)
dt_collapsed <- subset(dt_collapsed, group_review == 1)
dt_collapsed$ihme_loc_id <- gsub("_.*","",dt_collapsed$ihme_loc_id)

# add age information
dt_collapsed$age_difference <- dt_collapsed$age_end - dt_collapsed$age_start
# dt_collapsed <- subset(dt_collapsed, age_difference <=5)
dt_collapsed <- data.table(dt_collapsed)

# # subset to ONLY OLD studies
dt_collapsed <- subset(dt_collapsed, `in portnoy` == 0)
# 
# compute logit(cfr) using delta transform
dt_collapsed[, ratio := Deaths/Cases]
dt_collapsed[ratio==0, ratio := 0.0004404757/2]
dt_collapsed[ratio==1, ratio := 0.9999999999999]
dt_collapsed <- subset(dt_collapsed, ratio  <= 1)
dt_collapsed <- data.table(dt_collapsed)
z <- qnorm(0.975)
dt_collapsed[, ratio_se := sqrt((Deaths/Cases)*(1-(Deaths/Cases))/Cases + z^2/(4*Cases^2))] 
dt_collapsed <- as.data.frame(dt_collapsed)
dt_collapsed$logit_cfr <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,1]
dt_collapsed$logit_cfr_se <- delta_transform(mean = dt_collapsed$ratio, sd = dt_collapsed$ratio_se, transformation = "linear_to_logit")[,2]
dt_collapsed <- as.data.table(dt_collapsed)

# add final columns and clean up
dt_collapsed$Age.ind <- ifelse(dt_collapsed$age_start < 5, 1, 0)
dt_collapsed$Comm.ind <- ifelse(dt_collapsed$CV_Hospital == 1, 0, 1)

dt_collapsed$Midpoint.Year <- floor((as.numeric(dt_collapsed$year_start) + as.numeric(dt_collapsed$year_end))/2)

dt_data <- subset(dt_collapsed, select=c("ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
                                         "Age.ind", "Comm.ind", "logit_cfr", "logit_cfr_se", "ratio_se"))

need_averages <- subset(dt_data, select=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end", "Comm.ind"))
need_averages <- unique(need_averages)
need_averages  <- data.table(need_averages)

pred_frame2 <- data.table(pred_frame)

age_incidence <- fread('FILEPATH')

need_averages$cfr <- 0
for(i in 1:dim(need_averages)[1]){
  message(i / dim(need_averages)[1])
  
  loc <- as.character(need_averages[i,1] )
  yr <- as.numeric(need_averages[i,2])
  community <- as.numeric(need_averages[i,5])
  age_min <- as.numeric(need_averages[i,3])
  if(age_min < 1){
    age_min = 0
  }
  age_max <- as.numeric(need_averages[i,4])
  if(age_max < 1){
    age_max = 0
  }
  
  p_df <- subset(pred_frame2, ihme_loc_id == loc  & Midpoint.Year == yr & Comm.ind == community & age_start >= age_min & age_end <= age_max)
  age_df <- subset(age_incidence, country == loc  & year == yr & age >= age_min & age <= age_max)
  age_df <- subset(age_df, select=c('age', 'cases'))  
  
  p_df2 <- merge(p_df, age_df, by.x="age_start", by.y= "age")  
  
  p_df2$deaths <- p_df2$cfr_predictions * p_df2$cases
  
  case_total <- sum( p_df2$cases)
  death_total <- sum( p_df2$deaths)
  
  
  overall_cfr <- death_total / case_total
  
  
  need_averages$cfr[i] <- overall_cfr
  
}

need_averages <- subset(need_averages, cfr != "NaN")

ho_metrics <- data.table()

dt_data$ihme_loc_id <- gsub("_.*","",dt_data$ihme_loc_id)
dt_data <- subset(dt_data, age_start < 35)
dt_data <- merge(dt_data, need_averages, by.x=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end","Comm.ind"), by.y=c('ihme_loc_id', "Midpoint.Year", "age_start", "age_end","Comm.ind"))
dt_data$ratio <- invlogit(dt_data$logit_cfr)

data_match <- dt_data 

cor.val <- cor(data_match$cfr, data_match$ratio)
rmse.val <- sqrt(mean((data_match$ratio - data_match$cfr)^2)) 
me.val <- mean(data_match$ratio - data_match$cfr)
mae.val <- mean(abs(data_match$ratio - data_match$cfr))

to_add <- data.table(cbind(cor.val, rmse.val, me.val, mae.val))
ho_metrics <- rbind(ho_metrics, to_add)

fwrite(ho_metrics, paste0('FILEPATH',run_date,'/is_summary_holdout_metrics_originaldata.csv'))


