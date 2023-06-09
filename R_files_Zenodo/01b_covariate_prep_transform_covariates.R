
#################################################################################################
## Measles CFR estimation -- 01b_covariate_prep_transform_covariates.R
## Purpose: a) subset to LMICs only, 
##          b) determine transformation via regressions with best AIC score,
##          c) transform and standardize covariates, and 
##          d) clean up and save
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")

run_date_01a <- 'RUN_DATE'

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
covariates_to_include <- fread(paste0('FILEPATH', run_date_01a,'/square_covariate_set.csv'))

# -----------------------------------------------------------------------------------------------
# add on data -----------------------------------------------------------------------------------
dt_data2 <- merge(dt_data, covariates_to_include, all.x=T, all.y=F, by.x=c("ihme_loc_id", "Midpoint.Year"), by.y=c("ihme_loc_id", "year_id"))
dt_data2 <- subset(dt_data2, super_region_name != "High-income")

# -----------------------------------------------------------------------------------------------
# transformation regressions --------------------------------------------------------------------

dir.create(paste0("FILEPATH",run_date,"/diagnostics/regression_transformations/"), recursive=T)
aic_df <- data.table()

# first, do list with logit as an option
cov_list <- c("vitA_deficiency", "mcv1", "prop_urban", "hiv", "war_rate", "wasting")

for(covariate_name in cov_list){
  
  message(covariate_name)
  dt_covariate <- subset(dt_data2, select=c("logit_cfr", covariate_name, "Cases"))
  colnames(dt_covariate) <- c("logit_cfr", "covariate", "cases")
  
  predict_dt <- subset(covariates_to_include, select=c("ihme_loc_id", "year_id","super_region_name", covariate_name))
  colnames(predict_dt) <- c("ihme_loc_id", "year_id","super_region_name", "covariate")
  
  dt_covariate$covariate <- ifelse(dt_covariate$covariate == 0, 0.00000001, dt_covariate$covariate)
  predict_dt$covariate <- ifelse(predict_dt$covariate == 0, 0.00000001, predict_dt$covariate)
  
  test_lm_untransformed <- lm(logit_cfr ~ covariate, data=dt_covariate, weights=cases)
  predict_dt$untrans_predictions_logit_cfr <- predict.lm(test_lm_untransformed, predict_dt)
  predict_dt$untrans_predictions_cfr <- invlogit(predict_dt$untrans_predictions_logit_cfr)
  aic_untransformed <- AIC(test_lm_untransformed)
  
  test_lm_log <- lm(logit_cfr ~ log(covariate), data=dt_covariate, weights=cases)
  predict_dt$log_predictions_logit_cfr <- predict.lm(test_lm_log, predict_dt)
  predict_dt$log_predictions_cfr <- invlogit(predict_dt$log_predictions_logit_cfr)
  aic_log <- AIC(test_lm_log)
  
  test_lm_logit <- lm(logit_cfr ~ logit(covariate), data=dt_covariate, weights=cases)
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
  
  png(file = paste0("FILEPATH",run_date,"/diagnostics/regression_transformations/",covariate_name,".png"),
      width = 16,
      height = 8,
      units = "in", 
      res = 300)
  grid.draw(plot_all)
  dev.off()
  
  to_add <- cbind(covariate_name, aic_untransformed, aic_log, aic_logit)
  aic_df <- rbind(aic_df, to_add)
  
}


# now for list, w/o logit as an option
cov_list <- c("maternal_education", "tfr", "gdp_pc", "u5mr")

for(covariate_name in cov_list){
  
  message(covariate_name)
  dt_covariate <- subset(dt_data2, select=c("logit_cfr", covariate_name, "Cases"))
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
  
  aic_logit <- NA
  
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
  gg1
  
  gg2 <- ggplot(predict_dt, aes(y=log_predictions_cfr, x=covariate)) +
    geom_point(alpha=0.008) + geom_line(color='darkgrey') + theme_light() + ylab("CFR") + xlab(covariate_name) + ggtitle("Log covariate") + 
    geom_point(data = dt_covariate, mapping = aes(x = covariate, y = invlogit(logit_cfr)), alpha=0.2, color="darkgreen")
  gg2
  
  lay <- rbind(c(1,2),
               c(1,2),
               c(3,NA),
               c(3,NA))
  library(gridExtra)
  library(grid)
  
  plot_all <- arrangeGrob(gg_density, gg1, gg2,
                          layout_matrix = lay,
                          heights = c(1,1,1,1))
  
  png(file = paste0("FILEPATH", run_date,"/diagnostics/regression_transformations/",covariate_name,".png"),
      width = 16,
      height = 8,
      units = "in", 
      res = 300)
  grid.draw(plot_all)
  dev.off()
  
  to_add <- cbind(covariate_name, aic_untransformed, aic_log, aic_logit)
  aic_df <- rbind(aic_df, to_add)
}

# save df of aic scores
fwrite(aic_df, paste0('FILEPATH', run_date,'/aic_df_test_transformations.csv'))


# -----------------------------------------------------------------------------------------------
# make additional plots to compare transformations ----------------------------------------------

dir.create(paste0("FILEPATH", run_date,"/diagnostics/covariate_relationships/"))

# mcv1 coverage
gg_data_mcv1_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=mcv1, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("MCV1 coverage") 

gg_data_mcv1_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(mcv1), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(MCV1 coverage)") 

gg_data_mcv1_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(mcv1), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(MCV1 coverage)") 

png(file = paste0("FILEPATH", run_date,"/diagnostics/covariate_relationships/gg_data_mcv1_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_mcv1_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_mcv1_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_mcv1_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_mcv1_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_mcv1_logit)
dev.off()

# vitamin A deficiency
gg_data_vitA_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=vitA_deficiency, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("vitA_deficiency") 

gg_data_vitA_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(vitA_deficiency), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(vitA_deficiency)") 

gg_data_vitA_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(vitA_deficiency), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(vitA_deficiency)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_vitA_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_vitA_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_vitA_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_vitA_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_vitA_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_vitA_logit)
dev.off()

# war rate
gg_data_war_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=war_rate, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("war_rate") 

gg_data_war_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(war_rate), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(war_rate)") 

gg_data_war_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(war_rate), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(war_rate)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_war_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_war_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_war_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_war_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_war_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_war_logit)
dev.off()

# wasting 
gg_data_wasting_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=wasting, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("wasting") 

gg_data_wasting_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(wasting), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(wasting)") 

gg_data_wasting_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(wasting), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(wasting)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_wasting_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_wasting_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_wasting_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_wasting_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_wasting_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_wasting_logit)
dev.off()

# hiv
gg_data_hiv_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=hiv, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("hiv") 

gg_data_hiv_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(hiv), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(hiv)") 

gg_data_hiv_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(hiv), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(hiv)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_hiv_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_hiv_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_hiv_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_hiv_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_hiv_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_hiv_logit)
dev.off()

# maternal education
gg_data_edu_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=maternal_education, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("maternal_education") 

gg_data_edu_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(maternal_education), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(maternal_education)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_edu_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_edu_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_edu_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_edu_log)
dev.off()

# gdp per capita
gg_data_gdp_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=gdp_pc, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("gdp_pc") 

gg_data_gdp_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(gdp_pc), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(gdp_pc)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_gdp_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_gdp_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_gdp_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_gdp_log)
dev.off()

# prop. urban
gg_data_urban_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=prop_urban, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("prop_urban") 

gg_data_urban_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(prop_urban), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(prop_urban)") 

gg_data_urban_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(prop_urban), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(prop_urban)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_urban_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_urban_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_urban_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_urban_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_urban_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_urban_logit)
dev.off()

# u5mr
gg_data_u5mr_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=u5mr, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("u5mr") 

gg_data_u5mr_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(u5mr), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(u5mr)") 

gg_data_u5mr_logit <- ggplot(dt_data2, aes(y=logit_cfr, x=logit(u5mr/1000), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("logit(u5mr / 1000)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_u5mr_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_u5mr_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_u5mr_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_u5mr_log)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_u5mr_logit.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_u5mr_logit)
dev.off()

# tfr
gg_data_tfr_regular <- ggplot(dt_data2, aes(y=logit_cfr, x=tfr, size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("tfr") 

gg_data_tfr_log <- ggplot(dt_data2, aes(y=logit_cfr, x=log(tfr), size=Cases, shape=as.factor(Age.ind), color=as.factor(Comm.ind))) +
  geom_point(alpha=0.5) + theme_light() + ylab("logit(CFR)") + xlab("log(tfr)") 

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_tfr_regular.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_tfr_regular)
dev.off()

png(file = paste0("FILEPATH", run_date, "/diagnostics/covariate_relationships/gg_data_tfr_log.png"),
    width = 8,
    height = 5,
    units = "in",
    res = 400,
    type = "cairo")
print(gg_data_tfr_log)
dev.off()


# -----------------------------------------------------------------------------------------------
# restrict covariate set to only LMICs ----------------------------------------------------------
covariates_to_include <- subset(covariates_to_include, super_region_name != "High-income")

# -----------------------------------------------------------------------------------------------
# apply transformations based on aic_df and plots comparison ------------------------------------

# vitA deficiency -- untransformed

# mcv -- untransformed

# prop. urban -- untransformed

# hiv -- log transform 
covariates_to_include$hiv <- ifelse(covariates_to_include$hiv == 0, 0.000000000001, covariates_to_include$hiv)
covariates_to_include$log_hiv <- log(covariates_to_include$hiv)
covariates_to_include$hiv <- NULL

# war rate -- logit transform
covariates_to_include$war_rate <- ifelse(covariates_to_include$war_rate == 0, 0.000000000001, covariates_to_include$war_rate)
covariates_to_include$logit_war_rate <- logit(covariates_to_include$war_rate)
covariates_to_include$war_rate <- NULL

# wasting -- logit transform
covariates_to_include$logit_wasting <- logit(covariates_to_include$wasting)
covariates_to_include$wasting <- NULL

# education -- untransformed

# tfr -- untransformed

# gdp_pc -- log transform
covariates_to_include$log_gdp_pc <- log(covariates_to_include$gdp_pc)
covariates_to_include$gdp_pc <- NULL

# u5mr -- untransformed

# -----------------------------------------------------------------------------------------------
# clean up extra columns ------------------------------------------------------------------------
covariates_to_include$interpolated_gdp_pc <- NULL
covariates_to_include$interpolated_u5mr <- NULL
covariates_to_include$interpolated_prop_urban <- NULL
covariates_to_include$mcv1_region <- NULL

covariates_to_include$tfr_region <- NULL
covariates_to_include$gdp_pc_region <- NULL
covariates_to_include$u5mr_region <- NULL
covariates_to_include$prop_urban_region <- NULL
covariates_to_include$interpolated_mcv1 <- NULL
covariates_to_include$interpolated_tfr <- NULL

# -----------------------------------------------------------------------------------------------
# standardize covariates (subtract mean and divide by SD) ------------------------------------

# education 
covariates_to_include$maternal_education_standardized <- (covariates_to_include$maternal_education - mean(covariates_to_include$maternal_education)) / sd(covariates_to_include$maternal_education)

# gdp pc 
covariates_to_include$gdp_pc_standardized <- (covariates_to_include$log_gdp_pc - mean(covariates_to_include$log_gdp_pc)) / sd(covariates_to_include$log_gdp_pc)
covariates_to_include$log_gdp_pc <- NULL

# hiv 
covariates_to_include$hiv_standardized <- (covariates_to_include$log_hiv - mean(covariates_to_include$log_hiv)) / sd(covariates_to_include$log_hiv)
covariates_to_include$log_hiv <- NULL

# mcv, and also save mean and sd for no vaccination scenario 
fwrite(data.table(mean(covariates_to_include$mcv1, na.rm=T)), paste0('FILEPATH', run_date,'/mcv_mean.csv'))
fwrite(data.table(sd(covariates_to_include$mcv1, na.rm=T)), paste0('FILEPATH', run_date,'/mcv_sd.csv'))

covariates_to_include$mcv1_standardized <- (covariates_to_include$mcv1 - mean(covariates_to_include$mcv1)) / sd(covariates_to_include$mcv1)
covariates_to_include$mcv1 <- NULL

# tfr
covariates_to_include$tfr_standardized <- (covariates_to_include$tfr - mean(covariates_to_include$tfr)) / sd(covariates_to_include$tfr)
covariates_to_include$tfr <- NULL

# u5mr 
covariates_to_include$u5mr_standardized <- (covariates_to_include$u5mr - mean(covariates_to_include$u5mr)) / sd(covariates_to_include$u5mr)
covariates_to_include$u5mr <- NULL

# prop. urban 
covariates_to_include$prop_urban_standardized <- (covariates_to_include$prop_urban - mean(covariates_to_include$prop_urban)) / sd(covariates_to_include$prop_urban)
covariates_to_include$prop_urban <- NULL

# vitA deficiency 
covariates_to_include$vitA_standardized <- (covariates_to_include$vitA_deficiency - mean(covariates_to_include$vitA_deficiency)) / sd(covariates_to_include$vitA_deficiency)
covariates_to_include$vitA_deficiency <- NULL

# war rate
covariates_to_include$war_rate_standardized <- (covariates_to_include$logit_war_rate - mean(covariates_to_include$logit_war_rate)) / sd(covariates_to_include$logit_war_rate)
covariates_to_include$logit_war_rate <- NULL

# wasting
covariates_to_include$wasting_standardized <- (covariates_to_include$logit_wasting - mean(covariates_to_include$logit_wasting)) / sd(covariates_to_include$logit_wasting)
covariates_to_include$logit_wasting <- NULL

fwrite(covariates_to_include, paste0('FILEPATH',run_date,'/transformed_standardized_covaraite_set.csv'))



