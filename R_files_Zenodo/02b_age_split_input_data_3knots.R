
#################################################################################################
## Measles CFR estimation -- 02b_age_split_input_data.R
## Purpose: a) compute relative age patterns,
##          b) age split data using age specific incidence 
##          c) save input data to use in second stage model
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())

# date
run_date <- Sys.Date()
run_date <- format(run_date, "%Y_%m_%d")
run_date <- paste0(run_date, '_3knots')

run_date_02a <- 'RUN_DATE'

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
load(paste0('FILEPATH',run_date_02a,'/IS_pred_frame_age_curve.RData'))
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

# also make for original data (for comparison purposes)
original_dt_collapsed <- data.table(original_dt_collapsed)
original_dt_collapsed[, ratio := Deaths/Cases]
original_dt_collapsed[ratio==0, ratio := 0.0004404757/2]
original_dt_collapsed[ratio==1, ratio := 0.9999999999999]
original_dt_collapsed <- subset(original_dt_collapsed, ratio  <= 1)

original_dt_collapsed <- data.table(original_dt_collapsed)
z <- qnorm(0.975)
original_dt_collapsed[, ratio_se := sqrt((Deaths/Cases)*(1-(Deaths/Cases))/Cases + z^2/(4*Cases^2))] 
original_dt_collapsed <- as.data.frame(original_dt_collapsed)

original_dt_collapsed$logit_cfr <- delta_transform(mean = original_dt_collapsed$ratio, sd = original_dt_collapsed$ratio_se, transformation = "linear_to_logit")[,1]
original_dt_collapsed$logit_cfr_se <- delta_transform(mean = original_dt_collapsed$ratio, sd = original_dt_collapsed$ratio_se, transformation = "linear_to_logit")[,2]
original_dt_collapsed$logit_upper_cfr <- original_dt_collapsed$logit_cfr + (1.96*original_dt_collapsed$logit_cfr_se)
original_dt_collapsed$logit_lower_cfr <- original_dt_collapsed$logit_cfr - (1.96*original_dt_collapsed$logit_cfr_se)
original_dt_collapsed$upper_cfr  <- invlogit(original_dt_collapsed$logit_upper_cfr)
original_dt_collapsed$lower_cfr  <- invlogit(original_dt_collapsed$logit_lower_cfr)

# -----------------------------------------------------------------------------------------------
# dignostic plots -------------------------------------------------------------------------------

dir.create(paste0('FILEPATH',run_date,'/splitting_plots/'), recursive = T)
for(i in 1:length(unique(original_dt_collapsed$nid))){
  
  message(i / length(unique(original_dt_collapsed$nid)))
  og <- subset(original_dt_collapsed, nid == unique(original_dt_collapsed$nid)[i])
  split <- subset(dt_collapsed3, nid == unique(original_dt_collapsed$nid)[i])
  
  citation <- unique(og$field_citation_value)
  
  by_age <- ggplot() + ggtitle(unique(original_dt_collapsed$field_citation_value)[i]) + theme_bw() + facet_wrap(~year_start) + 
    geom_segment(data=split, aes(x=age_start, xend=age_end+1, y=ratio, yend=ratio), color='darkblue') +
    geom_segment(data=og, aes(x=age_start, xend=age_end, y=ratio, yend=ratio), color='red', alpha=0.4) +
    geom_segment(data=og, aes(x=(age_start + age_end)/2, xend=(age_start + age_end)/2, y=lower_cfr, yend=upper_cfr), color='red', alpha=0.4)+
    geom_segment(data=split, aes(x=(age_start + age_end + 1)/2, xend=(age_start + age_end + 1)/2, y=lower_cfr, yend=upper_cfr), alpha=0.4, color='darkblue')
  
  png(file = paste0('FILEPATH',run_date,'/splitting_plots/split_source_',i,'.png'),
      width = 12,
      height = 6,
      units = "in", 
      res = 300)
  plot(by_age)
  dev.off()
  
}

dt_collapsed <- NULL
dt_collapsed <- copy(dt_collapsed3)

fwrite(dt_collapsed, paste0('FILEPATH',run_date, '/age_split_input_data.csv'))



