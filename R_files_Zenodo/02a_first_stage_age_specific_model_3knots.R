
#################################################################################################
## Measles CFR estimation -- 02a_first_stage_age_specific_model_3knots.R
## Purpose: a) subset input data to 5 year age bins or smaller,
##          b) run first step age specific MR-BRT model
##          c) save relative age pattern / weights for age splitting 
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

# date
run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")
run_date <- paste0(run_date, '_3knots')
run_date_01c <- 'RUN_DATE'

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
dt_collapsed <- subset(dt_collapsed, age_difference <=5)
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

dt_data <- subset(dt_collapsed, select=c("ihme_loc_id", "Midpoint.Year", "Cases", "Deaths", "age_start", "age_end",
                                         "Age.ind", "Comm.ind", "logit_cfr", "logit_cfr_se", "ratio_se"))

dir.create(paste0('FILEPATH',run_date))
fwrite(dt_collapsed, paste0('FILEPATH',run_date,'/age_specific_model_dataset.csv'))

# -----------------------------------------------------------------------------------------------
# read in covariates ----------------------------------------------------------------------------------
covariates <- fread(paste0('FILEPATH',run_date_01c,'/transformed_standardized_covaraite_set_with_incidence.csv'))
data_test <- merge(dt_data, covariates, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)

data_test <- subset(data_test, super_region_name != 'High-income')

# -----------------------------------------------------------------------------------------------
# prepare prediction frame ----------------------------------------------------------------------
pred_frame_comm0 <- subset(covariates, year_id > 1989)
pred_frame_comm0$Comm.ind <- 0
pred_frame_comm1 <- subset(covariates, year_id > 1989)
pred_frame_comm1$Comm.ind <- 1

pred_frame <- rbind(pred_frame_comm0, pred_frame_comm1)

# -----------------------------------------------------------------------------------------------
# set up for MR-MRT -----------------------------------------------------------------------------

data_test$Comm.ind <- ifelse(is.na(data_test$Comm.ind), 0, data_test$Comm.ind)
data_test <- subset(data_test, Midpoint.Year > 1969)

missing <- data_test[!complete.cases(data_test),]

dt2 <- data_test[complete.cases(data_test),]
dt2$Country <- as.character(dt2$ihme_loc_id)
dt2$for_re <- as.character(dt2$ihme_loc_id)

dt2$age_end <- as.numeric(dt2$age_end)
dt2$age_start <- as.numeric(dt2$age_start)

# -----------------------------------------------------------------------------------------------
# compare covariate values per study ------------------------------------------------------------

mcv1_by_age <- ggplot(dt2, aes((age_start + age_end)/2, mcv1_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("MCV1 coverage") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y=mcv1_standardized, yend=mcv1_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

war_by_age <- ggplot(dt2, aes((age_start + age_end)/2, war_rate_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("War and terrorism mortality rate") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y= war_rate_standardized, yend=war_rate_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

wasting_by_age <- ggplot(dt2, aes((age_start + age_end)/2, wasting_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("Wasting prevalence") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = wasting_standardized, yend=wasting_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

edu_by_age <- ggplot(dt2, aes((age_start + age_end)/2, maternal_education_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("Maternal education") +
  geom_segment(data=dt2, aes(x = age_start, xend = age_end, y = maternal_education_standardized, yend = maternal_education_standardized, color = invlogit(logit_cfr), size = 1/ratio_se), alpha = 0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

gdp_by_age <- ggplot(dt2, aes((age_start + age_end)/2, gdp_pc_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("GDP per capita") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y=gdp_pc_standardized, yend=gdp_pc_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

hiv_by_age <- ggplot(dt2, aes((age_start + age_end)/2, hiv_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("HIV prevalence") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = hiv_standardized, yend=hiv_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

tfr_by_age <- ggplot(dt2, aes((age_start + age_end)/2, tfr_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("TFR") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = tfr_standardized, yend=log_tfr_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

u5mr_by_age <- ggplot(dt2, aes((age_start + age_end)/2, u5mr_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("U5MR") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = u5mr_standardized, yend=u5mr_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

urban_by_age <- ggplot(dt2, aes((age_start + age_end)/2, prop_urban_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("Proportion urban") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = prop_urban_standardized,yend=prop_urban_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

vitA_by_age <- ggplot(dt2, aes((age_start + age_end)/2, vitA_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("Vitamin A deficiency prevalence") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y=vitA_standardized, yend=vitA_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

incidence_by_age <- ggplot(dt2, aes((age_start + age_end)/2, incidence_standardized)) + labs(color="CFR", size="1 / SE") + 
  theme_bw() + xlab("Age") + ylab("Measles incidence") +
  geom_segment(data=dt2, aes(x=age_start, xend=age_end, y = incidence_standardized, yend= incidence_standardized, color=invlogit(logit_cfr), size=1/ratio_se), alpha=0.6) +  
  scale_color_viridis_c() +  geom_smooth(color="darkgrey", se = FALSE) + geom_smooth(color="black", se = FALSE, aes(weight=(1/dt2$ratio_se))) 

dir.create(paste0('FILEPATH',run_date,'/age_covariate_relationships/'), recursive = T)
png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/incidence.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(incidence_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/vitaminA.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(vitA_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/prop_urban.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(urban_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/u5mr.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(u5mr_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/tfr.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(tfr_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/hiv_prevalence.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(hiv_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/gdp_pc.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(gdp_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/education.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(edu_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/wasting.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(wasting_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/war_rate.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(war_by_age)
dev.off()

png(file = paste0('FILEPATH',run_date,'/age_covariate_relationships/mcv1_rate.png'),
    width = 12,
    height = 6,
    units = "in", 
    res = 300)
plot(mcv1_by_age)
dev.off()


# -----------------------------------------------------------------------------------------------
# run first step model --------------------------------------------------------------------------

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
                     use_spline_intercept=FALSE, 
                     spline_degree = 2L,
                     spline_knots = array(seq(0, 1, length.out = 5)),  
                     spline_knots_type = "frequency",
                     spline_l_linear = FALSE,
                     spline_r_linear = TRUE)))

mod_cfr$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

# get samples
n_samples <- 1000L
samples1 <- mod_cfr$sample_soln(sample_size = n_samples)

# -----------------------------------------------------------------------------------------------
# check beta values for covariates --------------------------------------------------------------
betas <- mod_cfr$beta_soln[-c(14,15, 16)]
names <- mod_cfr$cov_model_names[-14]

beta_table <- data.table()
beta_table$betas <- betas
beta_table$names <- names


# -----------------------------------------------------------------------------------------------
# save knot locations --------------------------------------------------------------

get_knots <- function(model) {
  model$cov_models[[which(model$cov_model_names == tail(model$cov_model_names)[6])]]$spline$knots
}

model_knots <- as.data.frame(get_knots(mod_cfr))
fwrite(model_knots, paste0('FILEPATH',run_date,'/knot_locations.csv'))

# -----------------------------------------------------------------------------------------------
# set up prediction frame -----------------------------------------------------------------------
pred_frame <- subset(pred_frame, select=c("ihme_loc_id", "year_id",
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
save(pred_frame, file= paste0('FILEPATH',run_date,'/IS_pred_frame_age_curve.RData'))

# -----------------------------------------------------------------------------------------------
# make diagnostic plots -------------------------------------------------------------------------
dir.create(paste0('FILEPATH',run_date,'/model_diganostics/'))

all_cov_pred <- copy(pred_frame)
all_cov_pred <- subset(all_cov_pred, select=c("age_start", "age_end", "mod_cfr"))
all_cov_pred$all_cov_cfr <- all_cov_pred$mod_cfr
all_cov_pred$all_cov_cfr_relative <- all_cov_pred$all_cov_cfr / all_cov_pred$all_cov_cfr[1]

gg_all_reg <- ggplot() + theme_bw() + xlab("Age") + ylab("CFR") + 
  geom_line(data=all_cov_pred, aes(x=age_start, y=all_cov_cfr), color="darkblue") 
gg_all_relative <- ggplot() + theme_bw() + xlab("Age") + ylab("Relative CFR to 0 yos") + 
  geom_line(data=all_cov_pred, aes(x=age_start, y=all_cov_cfr_relative), color="darkblue") 

png(file = paste0('FILEPATH',run_date,'/model_diganostics/CFR_curve.png'),
    width = 6,
    height = 3,
    units = "in", 
    res = 300)
plot(gg_all_reg)
dev.off()

png(file = paste0('FILEPATH',run_date,'/model_diganostics/CFR_curve_relative.png'),
    width = 6,
    height = 3,
    units = "in", 
    res = 300)
plot(gg_all_relative)
dev.off()
