
#################################################################################################
## Measles CFR estimation -- 01c_covariate_prep_incidence_covariate.R
## Purpose: process incidence covariate like the others (transform, standardize, save)
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")
run_date_01b <- 'RUN_DATE'

# load MSCM packages -- must be done in this order
library(crosswalk002, lib.loc = "FILEPATH")
library(mrbrt002, lib.loc = "FILEPATH")

library(data.table)
library(readxl)
library(ggplot2)
library(reshape2)

library(arm)

# -----------------------------------------------------------------------------------------------
# read in CFR data to use for regressions -------------------------------------------------------
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

# -----------------------------------------------------------------------------------------------
# load covariate files --------------------------------------------------------------------------
covariates_to_include <- fread(paste0('FILEPATH',run_date_01b,'/transformed_standardized_covaraite_set.csv'))

# -----------------------------------------------------------------------------------------------
# add on data -----------------------------------------------------------------------------------
dt_data2 <- merge(dt_data, covariates_to_include, all.x=T, all.y=F, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"))
dt_data2 <- subset(dt_data2, super_region_name != "High-income")

# -----------------------------------------------------------------------------------------------
# add incidence to data -------------------------------------------------------------------------

incidence <- fread('FILEPATH')
incidence <- subset(incidence, select=c("year", "country", "cohort_size", "cases"))
incidence$cohort_size <- as.numeric(incidence$cohort_size)
incidence$cases <- as.numeric(incidence$cases)

incidence_all_age_cases <- aggregate(incidence$cases, by=list(incidence$year, incidence$country), FUN = sum)
colnames(incidence_all_age_cases) <- c("year", "country", "cases")
incidence_all_age_cohort <- aggregate(incidence$cohort_size, by=list(incidence$year, incidence$country), FUN = sum)
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
# add on data -----------------------------------------------------------------------------------
dt_data2 <- merge(dt_data, covariates_to_include, all.x=T, all.y=F, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"))
dt_data2 <- subset(dt_data2, super_region_name != "High-income")

# -----------------------------------------------------------------------------------------------
# transformation regressions --------------------------------------------------------------------

aic_df <- data.table()
covariate_name <- "incidence"
message(covariate_name)
dt_covariate <- subset(dt_data2, select=c("logit_cfr", "incidence", "Cases"))
colnames(dt_covariate) <- c("logit_cfr", "covariate", "cases")

predict_dt <- subset(covariates_to_include, select=c("ihme_loc_id", "year_id","super_region_name", covariate_name))
colnames(predict_dt) <- c("ihme_loc_id", "year_id","super_region_name", "covariate")

dt_covariate$covariate <- ifelse(dt_covariate$covariate == 0, 0.00000001, dt_covariate$covariate)
predict_dt$covariate <- ifelse(predict_dt$covariate == 0, 0.00000001, predict_dt$covariate)

test_lm_untransformed <- lm(logit_cfr ~ covariate, data=dt_covariate, weights = cases)
predict_dt$untrans_predictions_logit_cfr <- predict.lm(test_lm_untransformed, predict_dt)
predict_dt$untrans_predictions_cfr <- invlogit(predict_dt$untrans_predictions_logit_cfr)
aic_untransformed <- AIC(test_lm_untransformed)

test_lm_log <- lm(logit_cfr ~ log(covariate), data=dt_covariate, weights = cases)
predict_dt$log_predictions_logit_cfr <- predict.lm(test_lm_log, predict_dt)
predict_dt$log_predictions_cfr <- invlogit(predict_dt$log_predictions_logit_cfr)
aic_log <- AIC(test_lm_log)

test_lm_logit <- lm(logit_cfr ~ logit(covariate), data=dt_covariate, weights = cases)
predict_dt$logit_predictions_logit_cfr <- predict.lm(test_lm_logit, predict_dt)
predict_dt$logit_predictions_cfr <- invlogit(predict_dt$logit_predictions_logit_cfr)
aic_logit <- AIC(test_lm_logit)

test_df_density <- subset(covariates_to_include, select=covariate_name)
test_df_density$variable <- 'all'

data_density <- subset(dt_data2, select=covariate_name)
data_density$variable <- 'data'

all_density <- rbind(test_df_density, data_density)
colnames(all_density) <- c("covariate", "variable")
gg_density <- ggplot(all_density, aes(x=covariate, fill=variable)) + theme_bw()+
  geom_density(alpha=.25) + xlab(covariate_name)

gg1 <- ggplot(predict_dt, aes(y=untrans_predictions_cfr, x=covariate)) +
  geom_point(alpha=0.008) + geom_line(color='darkgrey') + theme_light() + ylab("CFR") + xlab(covariate_name) + ggtitle("Untransformed covariate") + 
  geom_point(data = dt_covariate, mapping = aes(x = covariate, y = invlogit(logit_cfr)), alpha=0.2, color="darkgreen")

gg2 <- ggplot(predict_dt, aes(y=log_predictions_cfr, x=covariate)) +
  geom_point(alpha=0.008) + geom_line(color='darkgrey') + theme_light() + ylab("CFR") + xlab(covariate_name) + ggtitle("Log covariate") + 
  geom_point(data = dt_covariate, mapping = aes(x = covariate, y = invlogit(logit_cfr)), alpha=0.2, color="darkgreen")

gg3 <- ggplot(predict_dt, aes(y=logit_predictions_cfr, x=covariate)) +
  geom_point(alpha=0.008) + geom_line(color='darkgrey') + theme_light() + ylab("CFR") + xlab(covariate_name) + ggtitle("Logit covariate") + 
  geom_point(data = dt_covariate, mapping = aes(x = covariate, y = invlogit(logit_cfr)), alpha=0.2, color="darkgreen")

lay <- rbind(c(1,2),
             c(1,2),
             c(3,4),
             c(3,4))
library(gridExtra)
library(grid)

plot_all <- arrangeGrob(gg_density, gg1, gg2, gg3,
                        layout_matrix = lay,
                        heights = c(1,1,1,1))

dir.create(paste0('FILEPATH',run_date,'/diagnostics/regression_transformations/'), recursive = T)
png(file = paste0('FILEPATH',run_date,'/diagnostics/regression_transformations/',covariate_name,'.png'),
    width = 16,
    height = 8,
    units = "in", 
    res = 300)
grid.draw(plot_all)
dev.off()

aic_df <- cbind(covariate_name, aic_untransformed, aic_log, aic_logit)

fwrite(aic_df, paste0('FILEPATH',run_date,'/aic_df_test_transformations_incidence.csv'))

# -----------------------------------------------------------------------------------------------
# make additional plots to compare transformations ----------------------------------------------
dir.create(paste0('FILEPATH',run_date,'/diagnostics/covariate_relationships/'), recursive = T)

gg_data_incidence_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=incidence, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("incidence") 

gg_data_incidence_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(incidence), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(incidence)") 

gg_data_incidence_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(incidence), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(incidence)") 

png(file = paste0("FILEPATH",run_date,"/diagnostics/covariate_relationships/gg_data_incidence_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_incidence_regular)
dev.off()

png(file = paste0("FILEPATH",run_date,"/diagnostics/covariate_relationships/gg_data_incidence_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_incidence_log)
dev.off()

png(file = paste0("FILEPATH",run_date,"/diagnostics/covariate_relationships/gg_data_incidence_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_incidence_logit)
dev.off()

# -----------------------------------------------------------------------------------------------
# apply transformation --------------------------------------------------------------------------
covariates_to_include$incidence <- ifelse(covariates_to_include$incidence == 0, 0.0000000001, covariates_to_include$incidence)
covariates_to_include$logit_incidence <- logit(covariates_to_include$incidence)
covariates_to_include$incidence <- NULL

# -----------------------------------------------------------------------------------------------
# standardize covariates (subtract mean and divide by SD) ---------------------------------------
covariates_to_include$incidence_standardized <- (covariates_to_include$logit_incidence - mean(covariates_to_include$logit_incidence, na.rm=T)) / sd(covariates_to_include$logit_incidence, na.rm=T)
fwrite(data.table(mean(covariates_to_include$logit_incidence, na.rm=T)), paste0('FILEPATH',run_date,'/incidence_mean.csv'))
fwrite(data.table(sd(covariates_to_include$logit_incidence, na.rm=T)), paste0('FILEPATH',run_date,'/incidence_sd.csv'))

covariates_to_include$logit_incidence <- NULL

# -----------------------------------------------------------------------------------------------
# save ------------------------------------------------------------------------------------------
fwrite(covariates_to_include, paste0('FILEPATH',run_date,'/transformed_standardized_covaraite_set_with_incidence.csv'))

