
#################################################################################################
## Measles CFR estimation -- 01a_covariate_prep_square_covariates.R
## Purpose: a) gather WB, IHME, and other covariates, 
##          b) project forwards, backwards, interpolate, and 
##          c) square covariates for years 1970 - 2019 
##          d) clean up and save
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up ----------------------------------------------------------------------------------------
rm(list=ls())
library(data.table)
library(dplyr)

run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")

# -----------------------------------------------------------------------------------------------
# load pre-existing covariate file --------------------------------------------------------------
load(file='FILEPATH')

# -----------------------------------------------------------------------------------------------
# replace former WB / WUENIC covariates with July 2022 updates ----------------------------------
wb_covs <- fread('FILEPATH')
wb_covs$`Series Code` <- NULL

wuenic <- fread('FILEPATH')
wuenic$antigen <- NULL
all_covariates$mcv1 <- NULL
all_covariates <- merge(all_covariates, wuenic, by.x=c("ihme_loc_id", "year_id"), by.y=c("iso3", "year"), all.x=T, all.y=F)
all_covariates$mcv1 <- all_covariates$wuenic_coverage / 100

gdp_pc <- fread('FILEPATH')
gdp_pc$`Series Code` <- NULL
gdp_pc$`Series Name` <- NULL
gdp_pc$`Country Name` <- NULL
long_gdp_pc <- melt(gdp_pc, id.vars=c('Country Code'))
long_gdp_pc$year_id  <- as.numeric(gsub("\\\\[YR.*","",long_gdp_pc$variable))
long_gdp_pc$variable <- NULL
setnames(long_gdp_pc, 'Country Code', 'ihme_loc_id')
setnames(long_gdp_pc, 'value', 'gdp_pc')
long_gdp_pc <- subset(long_gdp_pc, ihme_loc_id != "")
long_gdp_pc$gdp_pc <- as.numeric(long_gdp_pc$gdp_pc)

all_covariates$gdp_pc <- NULL
all_covariates <- merge(all_covariates, long_gdp_pc, by=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)

tfr <- subset(wb_covs, `Series Name` == 'Fertility rate, total (births per woman)')
tfr$`Series Name` <- NULL
tfr$`Country Name` <- NULL
long_tfr <- melt(tfr, id.vars=c('Country Code'))
long_tfr$year_id  <- as.numeric(gsub("\\\\[YR.*","",long_tfr$variable))
long_tfr$variable <- NULL
setnames(long_tfr, 'Country Code', 'ihme_loc_id')
setnames(long_tfr, 'value', 'tfr')
long_tfr <- subset(long_tfr, ihme_loc_id != "")
long_tfr$tfr <- as.numeric(long_tfr$tfr)

all_covariates$tfr <- NULL
all_covariates <- merge(all_covariates, long_tfr, by=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)

u5mr <- subset(wb_covs, `Series Name` == 'Mortality rate, under-5 (per 1,000 live births)')
u5mr$`Series Name` <- NULL
u5mr$`Country Name` <- NULL
long_u5mr <- melt(u5mr, id.vars=c('Country Code'))
long_u5mr$year_id  <- as.numeric(gsub("\\\\[YR.*","",long_u5mr$variable))
long_u5mr$variable <- NULL
setnames(long_u5mr, 'Country Code', 'ihme_loc_id')
setnames(long_u5mr, 'value', 'u5mr')
long_u5mr <- subset(long_u5mr, ihme_loc_id != "")
long_u5mr$u5mr <- as.numeric(long_u5mr$u5mr)

all_covariates$u5mr <- NULL
all_covariates <- merge(all_covariates, long_u5mr, by=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)

urban <- subset(wb_covs, `Series Name` == 'Urban population (% of total population)')
urban$`Series Name` <- NULL
urban$`Country Name` <- NULL
long_urban <- melt(urban, id.vars=c('Country Code'))
long_urban$year_id  <- as.numeric(gsub("\\\\[YR.*","",long_urban$variable))
long_urban$variable <- NULL
setnames(long_urban, 'Country Code', 'ihme_loc_id')
setnames(long_urban, 'value', 'prop_urban')
long_urban <- subset(long_urban, ihme_loc_id != "")
long_urban$prop_urban <- as.numeric(long_urban$prop_urban) / 100

all_covariates$prop_urban <- NULL
all_covariates <- merge(all_covariates, long_urban, by=c("ihme_loc_id", "year_id"), all.x=T, all.y=F)

# -----------------------------------------------------------------------------------------------
# load hiv prevalence from ihme database --------------------------------------------------------
list.hiv.files <- list.files('FILEPATH')

all.hiv <- data.table()
for(i in 1:length(list.hiv.files)){
  message(list.hiv.files[i])
  hi_test <- fread(paste0('FILEPATH',list.hiv.files[i]))
  hi_test <- subset(hi_test, sex_id == 3 & measure == "Prevalence" & metric == "Rate" & age_group_id == 22)
  hi_test <- subset(hi_test, select=c("year_id", "mean"))
  hi_test$file_name <- list.hiv.files[i]  
  
  all.hiv <- rbind(all.hiv, hi_test)
}

all.hiv$ihme_loc_id <- sub("_spectrum_prep.csv", "",all.hiv$file_name)
all.hiv$file_name <- NULL
all.hiv <- all.hiv[!grep("_",ihme_loc_id)]
colnames(all.hiv) <- c("year_id", "hiv_new", "ihme_loc_id")

all_covariates <- merge(all_covariates, all.hiv, by=c("year_id", "ihme_loc_id"), all.x=T, all.y=F)
all_covariates$hiv <- NULL
all_covariates$hiv <- all_covariates$hiv_new

# -----------------------------------------------------------------------------------------------
# only keep covariates that will be used in modeling  -------------------------------------------
covariates_to_include <- subset(all_covariates, select=c("ihme_loc_id", "year_id", "location_id",
                                                         "vitA_deficiency", "war_rate", "wasting", "hiv", 
                                                         "maternal_education", "gdp_pc", "mcv1", "prop_urban", 
                                                         "u5mr", "tfr"))

# -----------------------------------------------------------------------------------------------
# source region information  --------------------------------------------------------------------
source(paste0('FILEPATH'))
loc_info <- get_location_hierarchy(41)
loc_info <- subset(loc_info, select=c("ihme_loc_id", "location_name", "super_region_name"))
covariates_to_include <- merge(covariates_to_include, loc_info, by='ihme_loc_id', all.x=T, all.y=F)

## correct for some iso3s missing regions
covariates_to_include$super_region_name <- ifelse(covariates_to_include$ihme_loc_id %in% c("COK", "NIU","NRU", "PLW", "TKL", "TUV"), "Southeast Asia, East Asia, and Oceania",
                                                  ifelse(covariates_to_include$ihme_loc_id %in% c("KNA"), "Latin America and Caribbean",
                                                         ifelse(covariates_to_include$ihme_loc_id %in% c("MCO", "SMR"), "High-income" , covariates_to_include$super_region_name)   ))

# -----------------------------------------------------------------------------------------------
# determine which covariates need inputing  -----------------------------------------------------
missing_covs <- covariates_to_include[!complete.cases(covariates_to_include),]
### covariates to interpolate: gdp_pc, mcv1, prop_urban, u5mr, tfr

# -----------------------------------------------------------------------------------------------
# gdp_pc interpolations / projections
to_export <- data.table()
df <- copy(covariates_to_include)
df_gdp <- subset(df, select=c("year_id", "ihme_loc_id", "gdp_pc"))

gdp <- subset(missing_covs, is.na(gdp_pc))
gdp_test <- data.table(ftable(gdp$ihme_loc_id))
gdp_cty_list <- subset(gdp_test, Freq <10)$Var1

### first do the linear interpolations in middle years
for(i in 1:length(gdp_cty_list)){
  iso <- gdp_cty_list[i]
  message(iso)
  test_gdp <- data.table(subset(df_gdp, ihme_loc_id == iso))
  test_gdp <- setorder(test_gdp, year_id)
  
  test_gdp[,val_before := nafill(gdp_pc, "locf")]
  test_gdp[,val_after := nafill(gdp_pc, "nocb")]
  
  test_gdp[, rle := rleid(gdp_pc)][,missings := max(.N +  1 , 2), by = rle][]
  
  test_gdp[is.na(gdp_pc), gdp_pc := val_before + .SD[,.I] *
             (val_after - val_before)/(missings), by = rle]
  
  test_gdp$interpolated <- ifelse(test_gdp$val_before == test_gdp$val_after, 0, 1)
  
  to_export <- rbind(to_export, subset(test_gdp, select=c("year_id", "ihme_loc_id", "gdp_pc", "interpolated")))
}

countries_exported <- unique(to_export$ihme_loc_id)

still_missing_early <- subset(to_export, is.na(gdp_pc) & year_id == 1980)
still_missing_early_countries <-  unique(still_missing_early$ihme_loc_id)
still_missing_late <- subset(to_export, is.na(gdp_pc) & year_id == 2019)
still_missing_late_countries <-  unique(still_missing_late$ihme_loc_id)

`%!in%` = Negate(`%in%`)
df_gdp_FULL <- subset(df_gdp, ihme_loc_id %!in% countries_exported)

to_export_FULL <- subset(to_export, ihme_loc_id %!in% c(still_missing_early_countries, still_missing_late_countries))

to_export_NOT_FULL <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries, still_missing_late_countries))
to_export_NOT_FULL_EARLY <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries))
to_export_NOT_FULL_LATE <- subset(to_export, ihme_loc_id %in% c(still_missing_late_countries))

source("FILEPATH")

### now lets do historical projections 
to_collect <- data.table()
for(i in 1:length(unique(to_export_NOT_FULL_EARLY$ihme_loc_id))){
  iso <- unique(to_export_NOT_FULL_EARLY$ihme_loc_id)[i]
  message(iso)
  
  df_tester <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  df_tester$interpolated <- NULL
  finished_projections <- make_historic_projections(df_tester)
  
  interpolated_key <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  interpolated_key <- subset(interpolated_key, select=c("year_id", "interpolated"))
  
  finished_projections2 <- merge(finished_projections, interpolated_key, by='year_id')
  finished_projections2$interpolated <- ifelse(is.na(finished_projections2$interpolated), 1, finished_projections2$interpolated)
  to_collect <- rbind(to_collect, finished_projections2)
  
}

colnames(to_collect) <- c("year_id", "ihme_loc_id", "gdp_pc", "interpolated")

### now for forwards projections
to_collect_forwards <- data.table()
for(i in 1:length(unique(to_export_NOT_FULL_LATE$ihme_loc_id))){
  iso <- unique(to_export_NOT_FULL_LATE$ihme_loc_id)[i]
  message(iso)
  
  df_tester <- subset(to_export_NOT_FULL_LATE, ihme_loc_id == iso)
  df_tester$interpolated <- NULL
  finished_projections <- make_future_projections(df_tester)
  
  interpolated_key <- subset(to_export_NOT_FULL_LATE, ihme_loc_id == iso)
  interpolated_key <- subset(interpolated_key, select=c("year_id", "interpolated"))
  
  finished_projections2 <- merge(finished_projections, interpolated_key, by='year_id')
  finished_projections2$interpolated <- ifelse(is.na(finished_projections2$interpolated), 1, finished_projections2$interpolated)
  to_collect_forwards <- rbind(to_collect_forwards, finished_projections2)
  
}

colnames(to_collect_forwards) <- c("year_id", "ihme_loc_id", "gdp_pc", "interpolated")

### gathe all data back together
to_export_FULL <- rbind(to_export_FULL, to_collect)
to_export_FULL <- rbind(to_export_FULL, to_collect_forwards)

df_gdp_final <- subset(df_gdp, ihme_loc_id %!in% unique(to_export_FULL$ihme_loc_id))
df_gdp_final$interpolated <- 0
df_gdp_final <- rbind(df_gdp_final, to_export_FULL)
df_gdp_final$interpolated_gdp_pc <- df_gdp_final$interpolated
df_gdp_final$interpolated <- NULL

gdp_cty_list_need_regions <- subset(gdp_test, Freq >=10)$Var1

df_gdp_final$gdp_pc <- ifelse(df_gdp_final$ihme_loc_id %in% c(gdp_cty_list_need_regions), NA, df_gdp_final$gdp_pc)
df_gdp_final$interpolated_gdp_pc <- ifelse(df_gdp_final$ihme_loc_id %in% c(gdp_cty_list_need_regions), 2, df_gdp_final$interpolated_gdp_pc)

# still_missing_again <- subset(df_gdp_final, is.na(gdp_pc))
covariates_to_include$gdp_pc <- NULL
covariates_to_include <- merge(covariates_to_include,df_gdp_final, by=c("year_id", "ihme_loc_id"))

# -----------------------------------------------------------------------------------------------
# u5mr interpolations / projections

to_export <- data.table()
df <- copy(covariates_to_include)
df_u5mr <- subset(df, select=c("year_id", "ihme_loc_id", "u5mr"))

u5mr <- subset(missing_covs, is.na(u5mr))
u5mr_test <- data.table(ftable(u5mr$ihme_loc_id))
u5mr_cty_list <- subset(u5mr_test, Freq <10)$Var1

for(i in 1:length(u5mr_cty_list)){
  iso <- u5mr_cty_list[i]
  message(iso)
  test_u5mr <- data.table(subset(df_u5mr, ihme_loc_id == iso))
  test_u5mr <- setorder(test_u5mr, year_id)
  
  test_u5mr[,val_before := nafill(u5mr, "locf")]
  test_u5mr[,val_after := nafill(u5mr, "nocb")]
  
  test_u5mr[, rle := rleid(u5mr)][,missings := max(.N +  1 , 2), by = rle][]
  
  test_u5mr[is.na(u5mr), u5mr := val_before + .SD[,.I] *
              (val_after - val_before)/(missings), by = rle]
  
  test_u5mr$interpolated <- ifelse(test_u5mr$val_before == test_u5mr$val_after, 0, 1)
  
  to_export <- rbind(to_export, subset(test_u5mr, select=c("year_id", "ihme_loc_id", "u5mr", "interpolated")))
}


countries_exported <- unique(to_export$ihme_loc_id)

still_missing_early <- subset(to_export, is.na(u5mr) & year_id == 1980)
still_missing_early_countries <-  unique(still_missing_early$ihme_loc_id)
still_missing_late <- subset(to_export, is.na(u5mr) & year_id == 2019)
still_missing_late_countries <-  unique(still_missing_late$ihme_loc_id)

df_u5mr_FULL <- subset(df_u5mr, ihme_loc_id %!in% countries_exported)

to_export_FULL <- subset(to_export, ihme_loc_id %!in% c(still_missing_early_countries, still_missing_late_countries))

to_export_NOT_FULL <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries, still_missing_late_countries))
to_export_NOT_FULL_EARLY <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries))
to_export_NOT_FULL_LATE <- subset(to_export, ihme_loc_id %in% c(still_missing_late_countries))

to_collect <- data.table()
# for(i in 1:length(unique(still_missing_early$ihme_loc_id))){
for(i in 1:length(unique(to_export_NOT_FULL_EARLY$ihme_loc_id))){
  iso <- unique(to_export_NOT_FULL_EARLY$ihme_loc_id)[i]
  message(iso)
  
  df_tester <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  df_tester$interpolated <- NULL
  finished_projections <- make_historic_projections(df_tester)
  
  interpolated_key <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  interpolated_key <- subset(interpolated_key, select=c("year_id", "interpolated"))
  
  finished_projections2 <- merge(finished_projections, interpolated_key, by='year_id')
  finished_projections2$interpolated <- ifelse(is.na(finished_projections2$interpolated), 1, finished_projections2$interpolated)
  to_collect <- rbind(to_collect, finished_projections2)
  
}

colnames(to_collect) <- c("year_id", "ihme_loc_id", "u5mr", "interpolated")

to_export_FULL <- rbind(to_export_FULL, to_collect)

df_u5mr_final <- subset(df_u5mr, ihme_loc_id %!in% unique(to_export_FULL$ihme_loc_id))
df_u5mr_final$interpolated <- 0
df_u5mr_final <- rbind(df_u5mr_final, to_export_FULL)
df_u5mr_final$interpolated_u5mr <- df_u5mr_final$interpolated
df_u5mr_final$interpolated <- NULL


u5mr_cty_list_need_regions <- subset(u5mr_test, Freq >=10)$Var1

df_u5mr_final$u5mr <- ifelse(df_u5mr_final$ihme_loc_id %in% c(u5mr_cty_list_need_regions), NA, df_u5mr_final$u5mr)
df_u5mr_final$interpolated_u5mr <- ifelse(df_u5mr_final$ihme_loc_id %in% c(u5mr_cty_list_need_regions), 2, df_u5mr_final$interpolated_u5mr)

# still_missing_again <- subset(df_gdp_final, is.na(gdp_pc))
covariates_to_include$u5mr <- NULL
covariates_to_include <- merge(covariates_to_include,df_u5mr_final, by=c("year_id", "ihme_loc_id"))

# -----------------------------------------------------------------------------------------------
# tfr interpolations / projections
to_export <- data.table()
df <- copy(covariates_to_include)
df_tfr <- subset(df, select=c("year_id", "ihme_loc_id", "tfr"))

tfr <- subset(missing_covs, is.na(tfr))
tfr_test <- data.table(ftable(tfr$ihme_loc_id))
tfr_cty_list <- subset(tfr_test, Freq <10)$Var1

### all 25%  or more of time series, just use regional averages

# -----------------------------------------------------------------------------------------------
# mcv1 interpolations / projections
to_export <- data.table()
df <- copy(covariates_to_include)
df_mcv1 <- subset(df, select=c("year_id", "ihme_loc_id", "mcv1"))

mcv1 <- subset(missing_covs, is.na(mcv1))
mcv1_test <- data.table(ftable(mcv1$ihme_loc_id))
mcv1_cty_list <- subset(mcv1_test, Freq <10)$Var1

### first do the linear interpolations in middle years
for(i in 1:length(mcv1_cty_list)){
  iso <- mcv1_cty_list[i]
  message(iso)
  test_mcv1 <- data.table(subset(df_mcv1, ihme_loc_id == iso))
  test_mcv1 <- setorder(test_mcv1, year_id)
  
  test_mcv1[,val_before := nafill(mcv1, "locf")]
  test_mcv1[,val_after := nafill(mcv1, "nocb")]
  
  test_mcv1[, rle := rleid(mcv1)][,missings := max(.N +  1 , 2), by = rle][]
  
  test_mcv1[is.na(mcv1), mcv1 := val_before + .SD[,.I] *
              (val_after - val_before)/(missings), by = rle]
  
  test_mcv1$interpolated <- ifelse(test_mcv1$val_before == test_mcv1$val_after, 0, 1)
  
  to_export <- rbind(to_export, subset(test_mcv1, select=c("year_id", "ihme_loc_id", "mcv1", "interpolated")))
}

countries_exported <- unique(to_export$ihme_loc_id)

still_missing_early <- subset(to_export, is.na(mcv1) & year_id == 1980)
still_missing_early_countries <-  unique(still_missing_early$ihme_loc_id)
still_missing_late <- subset(to_export, is.na(mcv1) & year_id == 2019)
still_missing_late_countries <-  unique(still_missing_late$ihme_loc_id)

`%!in%` = Negate(`%in%`)
df_mcv1_FULL <- subset(df_mcv1, ihme_loc_id %!in% countries_exported)

to_export_FULL <- subset(to_export, ihme_loc_id %!in% c(still_missing_early_countries, still_missing_late_countries))

to_export_NOT_FULL <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries, still_missing_late_countries))
to_export_NOT_FULL_EARLY <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries))
to_export_NOT_FULL_LATE <- subset(to_export, ihme_loc_id %in% c(still_missing_late_countries))

### now lets do historical projections 
to_collect <- data.table()
for(i in 1:length(unique(to_export_NOT_FULL_EARLY$ihme_loc_id))){
  iso <- unique(to_export_NOT_FULL_EARLY$ihme_loc_id)[i]
  message(iso)
  
  df_tester <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  df_tester$interpolated <- NULL
  finished_projections <- make_historic_projections(df_tester)
  
  interpolated_key <- subset(to_export_NOT_FULL_EARLY, ihme_loc_id == iso)
  interpolated_key <- subset(interpolated_key, select=c("year_id", "interpolated"))
  
  finished_projections2 <- merge(finished_projections, interpolated_key, by='year_id')
  finished_projections2$interpolated <- ifelse(is.na(finished_projections2$interpolated), 1, finished_projections2$interpolated)
  to_collect <- rbind(to_collect, finished_projections2)
  
}

colnames(to_collect) <- c("year_id", "ihme_loc_id", "mcv1", "interpolated")

### gather all data back together
to_export_FULL <- rbind(to_export_FULL, to_collect)

df_mcv1_final <- subset(df_mcv1, ihme_loc_id %!in% unique(to_export_FULL$ihme_loc_id))
df_mcv1_final$interpolated <- 0
df_mcv1_final <- rbind(df_mcv1_final, to_export_FULL)
df_mcv1_final$interpolated_mcv1 <- df_mcv1_final$interpolated
df_mcv1_final$interpolated <- NULL

mcv1_cty_list_need_regions <- subset(mcv1_test, Freq >=10)$Var1

df_mcv1_final$mcv1 <- ifelse(df_mcv1_final$ihme_loc_id %in% c(mcv1_cty_list_need_regions), NA, df_mcv1_final$mcv1)
df_mcv1_final$interpolated_mcv1 <- ifelse(df_mcv1_final$ihme_loc_id %in% c(mcv1_cty_list_need_regions), 2, df_mcv1_final$interpolated_mcv1)

covariates_to_include$mcv1 <- NULL
covariates_to_include <- merge(covariates_to_include,df_mcv1_final, by=c("year_id", "ihme_loc_id"))

# -----------------------------------------------------------------------------------------------
# prop urban interpolations / projections

to_export <- data.table()
df <- copy(covariates_to_include)
df_prop_urban <- subset(df, select=c("year_id", "ihme_loc_id", "prop_urban"))

prop_urban <- subset(missing_covs, is.na(prop_urban))
prop_urban_test <- data.table(ftable(prop_urban$ihme_loc_id))
prop_urban_cty_list <- subset(prop_urban_test, Freq <10)$Var1

### first do the linear interpolations in middle years
for(i in 1:length(prop_urban_cty_list)){
  iso <- prop_urban_cty_list[i]
  message(iso)
  test_prop_urban <- data.table(subset(df_prop_urban, ihme_loc_id == iso))
  test_prop_urban <- setorder(test_prop_urban, year_id)
  
  test_prop_urban[,val_before := nafill(prop_urban, "locf")]
  test_prop_urban[,val_after := nafill(prop_urban, "nocb")]
  
  test_prop_urban[, rle := rleid(prop_urban)][,missings := max(.N +  1 , 2), by = rle][]
  
  test_prop_urban[is.na(prop_urban), prop_urban := val_before + .SD[,.I] *
                    (val_after - val_before)/(missings), by = rle]
  
  test_prop_urban$interpolated <- ifelse(test_prop_urban$val_before == test_prop_urban$val_after, 0, 1)
  
  to_export <- rbind(to_export, subset(test_prop_urban, select=c("year_id", "ihme_loc_id", "prop_urban", "interpolated")))
}

countries_exported <- unique(to_export$ihme_loc_id)

still_missing_early <- subset(to_export, is.na(prop_urban) & year_id == 1980)
still_missing_early_countries <-  unique(still_missing_early$ihme_loc_id)
still_missing_late <- subset(to_export, is.na(prop_urban) & year_id == 2019)
still_missing_late_countries <-  unique(still_missing_late$ihme_loc_id)

df_prop_urban_FULL <- subset(df_prop_urban, ihme_loc_id %!in% countries_exported)

to_export_FULL <- subset(to_export, ihme_loc_id %!in% c(still_missing_early_countries, still_missing_late_countries))

to_export_NOT_FULL <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries, still_missing_late_countries))
to_export_NOT_FULL_EARLY <- subset(to_export, ihme_loc_id %in% c(still_missing_early_countries))
to_export_NOT_FULL_LATE <- subset(to_export, ihme_loc_id %in% c(still_missing_late_countries))

### now for forwards projections
to_collect_forwards <- data.table()
for(i in 1:length(unique(to_export_NOT_FULL_LATE$ihme_loc_id))){
  iso <- unique(to_export_NOT_FULL_LATE$ihme_loc_id)[i]
  message(iso)
  
  df_tester <- subset(to_export_NOT_FULL_LATE, ihme_loc_id == iso)
  df_tester$interpolated <- NULL
  finished_projections <- make_future_projections(df_tester)
  
  interpolated_key <- subset(to_export_NOT_FULL_LATE, ihme_loc_id == iso)
  interpolated_key <- subset(interpolated_key, select=c("year_id", "interpolated"))
  
  finished_projections2 <- merge(finished_projections, interpolated_key, by='year_id')
  finished_projections2$interpolated <- ifelse(is.na(finished_projections2$interpolated), 1, finished_projections2$interpolated)
  to_collect_forwards <- rbind(to_collect_forwards, finished_projections2)
  
}

colnames(to_collect_forwards) <- c("year_id", "ihme_loc_id", "prop_urban", "interpolated")

### gathe all data back together
to_export_FULL <- rbind(to_export_FULL, to_collect_forwards)

df_prop_urban_final <- subset(df_prop_urban, ihme_loc_id %!in% unique(to_export_FULL$ihme_loc_id))
df_prop_urban_final$interpolated <- 0
df_prop_urban_final <- rbind(df_prop_urban_final, to_export_FULL)
df_prop_urban_final$interpolated_prop_urban <- df_prop_urban_final$interpolated
df_prop_urban_final$interpolated <- NULL

prop_urban_cty_list_need_regions <- subset(prop_urban_test, Freq >=10)$Var1

df_prop_urban_final$prop_urban <- ifelse(df_prop_urban_final$ihme_loc_id %in% c(prop_urban_cty_list_need_regions), NA, df_prop_urban_final$prop_urban)
df_prop_urban_final$interpolated_prop_urban <- ifelse(df_prop_urban_final$ihme_loc_id %in% c(prop_urban_cty_list_need_regions), 2, df_prop_urban_final$interpolated_prop_urban)

covariates_to_include$prop_urban <- NULL
covariates_to_include <- merge(covariates_to_include,df_prop_urban_final, by=c("year_id", "ihme_loc_id"))

# -----------------------------------------------------------------------------------------------
# calculate regional averages for rest of missing covariates ------------------------------------

covariates_to_include_regions <- copy(covariates_to_include)
covariates_to_include_regions <- subset(covariates_to_include_regions, select=c("super_region_name", "year_id", 
                                                                                "mcv1", "tfr", "gdp_pc", "u5mr", "prop_urban"))

covariates_to_include_regions_means <- aggregate(covariates_to_include_regions, by=list(covariates_to_include_regions$super_region_name, covariates_to_include_regions$year), FUN = 'mean', na.rm=T)

covariates_to_include_regions_means$super_region_name <- NULL
covariates_to_include_regions_means$Group.2 <- NULL
colnames(covariates_to_include_regions_means) <- c("super_region_name", "year_id", 
                                                   "mcv1_region", "tfr_region", "gdp_pc_region", "u5mr_region", "prop_urban_region")

FULL_covariate_set <- merge(covariates_to_include, covariates_to_include_regions_means, by=c("year_id", "super_region_name"))

FULL_covariate_set$gdp_pc <- ifelse(FULL_covariate_set$interpolated_gdp_pc == 2, FULL_covariate_set$gdp_pc_region, FULL_covariate_set$gdp_pc)
FULL_covariate_set$prop_urban <- ifelse(FULL_covariate_set$interpolated_prop_urban == 2, FULL_covariate_set$prop_urban_region, FULL_covariate_set$prop_urban)
FULL_covariate_set$u5mr <- ifelse(FULL_covariate_set$interpolated_u5mr == 2, FULL_covariate_set$u5mr_region, FULL_covariate_set$u5mr)

## still need to add mcv and tfr swaps
mcv1 <- subset(missing_covs, is.na(mcv1))
mcv1_test <- data.table(ftable(mcv1$ihme_loc_id))
mcv1_cty_list <- subset(mcv1_test, Freq >=10)$Var1

tfr <- subset(missing_covs, is.na(tfr))
tfr_test <- data.table(ftable(tfr$ihme_loc_id))
tfr_cty_list <- subset(tfr_test, Freq >=10)$Var1

FULL_covariate_set$interpolated_mcv1 <- ifelse(FULL_covariate_set$ihme_loc_id %in% c(mcv1_cty_list), 2, 0)
FULL_covariate_set$interpolated_tfr <- ifelse(FULL_covariate_set$ihme_loc_id %in% c(tfr_cty_list), 2, 0)

FULL_covariate_set$mcv1 <- ifelse(FULL_covariate_set$interpolated_mcv1 == 2, FULL_covariate_set$mcv1_region, FULL_covariate_set$mcv1)
FULL_covariate_set$tfr <- ifelse(FULL_covariate_set$interpolated_tfr == 2, FULL_covariate_set$tfr_region, FULL_covariate_set$tfr)

test_final <- FULL_covariate_set[!complete.cases(FULL_covariate_set),]
test_final$wasting <- NULL
test_final$location_name <- NULL
test_final[!complete.cases(test_final),]

# -----------------------------------------------------------------------------------------------
# carry wasting from 1990 backwards -------------------------------------------------------------
wasting_1990 <- subset(FULL_covariate_set, year_id == 1990)
wasting_1990 <- subset(wasting_1990, select=c("ihme_loc_id", "year_id", "super_region_name", "location_id", "location_name", "wasting"))

wasting_1980 <- copy(wasting_1990)
wasting_1980$year_id <- 1980

wasting_1981 <- copy(wasting_1990)
wasting_1981$year_id <- 1981

wasting_1982 <- copy(wasting_1990)
wasting_1982$year_id <- 1982

wasting_1983 <- copy(wasting_1990)
wasting_1983$year_id <- 1983

wasting_1984 <- copy(wasting_1990)
wasting_1984$year_id <- 1984

wasting_1985 <- copy(wasting_1990)
wasting_1985$year_id <- 1985

wasting_1986 <- copy(wasting_1990)
wasting_1986$year_id <- 1986

wasting_1987 <- copy(wasting_1990)
wasting_1987$year_id <- 1987

wasting_1988 <- copy(wasting_1990)
wasting_1988$year_id <- 1988

wasting_1989 <- copy(wasting_1990)
wasting_1989$year_id <- 1989

wasting_1990_onwards <- subset(FULL_covariate_set, year_id >= 1990)
wasting_1990_onwards <- subset(wasting_1990_onwards, select=c("ihme_loc_id", "year_id", "super_region_name", "location_id", "location_name", "wasting"))

wasting_full <- rbind(wasting_1980, wasting_1981, wasting_1982, wasting_1983, wasting_1984, wasting_1985, 
                      wasting_1986, wasting_1987, wasting_1988, wasting_1989, wasting_1990_onwards)


FULL_covariate_set$wasting <- NULL
FULL_covariate_set <- merge(FULL_covariate_set, wasting_full, by=c("ihme_loc_id", "year_id", "super_region_name", "location_id", "location_name"))


# -----------------------------------------------------------------------------------------------
# fill rest of covariates from 1980 backwards ---------------------------------------------------

all_covariates_1980 <- subset(FULL_covariate_set, year_id == 1980)

all_covariates_1970 <- copy(all_covariates_1980)
all_covariates_1970$year_id <- 1970

all_covariates_1971 <- copy(all_covariates_1980)
all_covariates_1971$year_id <- 1971

all_covariates_1972 <- copy(all_covariates_1980)
all_covariates_1972$year_id <- 1972

all_covariates_1973 <- copy(all_covariates_1980)
all_covariates_1973$year_id <- 1973

all_covariates_1974 <- copy(all_covariates_1980)
all_covariates_1974$year_id <- 1974

all_covariates_1975 <- copy(all_covariates_1980)
all_covariates_1975$year_id <- 1975

all_covariates_1976 <- copy(all_covariates_1980)
all_covariates_1976$year_id <- 1976

all_covariates_1977 <- copy(all_covariates_1980)
all_covariates_1977$year_id <- 1977

all_covariates_1978 <- copy(all_covariates_1980)
all_covariates_1978$year_id <- 1978

all_covariates_1979 <- copy(all_covariates_1980)
all_covariates_1979$year_id <- 1979


all_covariates_early <- rbind(all_covariates_1970, all_covariates_1971, all_covariates_1972, all_covariates_1973, 
                              all_covariates_1974, all_covariates_1975, all_covariates_1976, all_covariates_1977,
                              all_covariates_1978, all_covariates_1979 )

FINAL_covariate_set <- rbind(FULL_covariate_set, all_covariates_early)

dir.create(paste0('FILEPATH', run_date,'/'))
fwrite(FINAL_covariate_set, paste0('FILEPATH', run_date,'/square_covariate_set.csv'))


# -----------------------------------------------------------------------------------------------
# make diagnostic plots -------------------------------------------------------------------------

library(ggplot2)
library(reshape2)
library(RColorBrewer)

df <- FINAL_covariate_set

myColors <- brewer.pal(4,"Dark2")

dir.create(paste0('FILEPATH', run_date,'/diagnostics/'))
for(i in 1:length(unique(df$ihme_loc_id))){
  
  iso <- unique(df$ihme_loc_id)[i]
  message(iso)
  df2 <- subset(df, ihme_loc_id == iso)
  
  d <- melt(df2, id = c( 'ihme_loc_id', 'year_id', 'super_region_name', 'location_id', 'location_name', 'interpolated_mcv1', 'interpolated_tfr', 'interpolated_u5mr', 'interpolated_prop_urban', 'interpolated_gdp_pc'))
  
  d <- subset(d, variable != "u5mr_region")
  d <- subset(d, variable != "mcv1_region")
  d <- subset(d, variable != "prop_urban_region")
  d <- subset(d, variable != "tfr_region")
  d <- subset(d, variable != "gdp_pc_region")
  
  d$interpolated <- ifelse(d$variable %in% c("vitA_deficiency", "war_rate", "wasting", "hiv",
                                             "maternal_education"), 0, 
                           ifelse(d$variable == "mcv1", d$interpolated_mcv1,
                                  ifelse(d$variable == "tfr", d$interpolated_tfr, 
                                         ifelse(d$variable == "gdp_pc", d$interpolated_gdp_pc,  
                                                ifelse(d$variable == "u5mr", d$interpolated_u5mr, 
                                                       ifelse(d$variable == "prop_urban", d$interpolated_prop_urban, 4))))))
  
  d$interpolated_factor <- ifelse(d$interpolated == 0, "Original value",
                                  ifelse(d$interpolated == 1, "Interpolated",
                                         ifelse(d$interpolated == 2, "Regional average", "MISSING")))
  
  gg_ctry <- ggplot(d, aes(x = year_id, y = value)) + geom_line(color='grey28') + theme_bw()+
    geom_point(aes(color=as.factor(interpolated_factor))) + xlab("Year") + ylab("Value") +
    facet_wrap(~ variable, scales = 'free') + scale_color_manual(values=c("Original value" = myColors[1],
                                                                          "Interpolated" = myColors[2],
                                                                          "Regional average" = myColors[3]), drop=T) + 
    guides(color=guide_legend(title="")) + ggtitle(paste0(unique(d$location_name), " (",iso, ")"))
  
  
  png(file = paste0("FILEPATH",run_date,"/diagnostics/",iso,".png"),
      width = 12,
      height = 8,
      units = "in",
      res = 400,
      type = "cairo")
  print(gg_ctry)
  dev.off()
  
}



