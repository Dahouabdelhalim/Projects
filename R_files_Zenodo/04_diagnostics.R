
#################################################################################################
## Measles CFR estimation -- 04_diagnostics.R
## Purpose: make diagnostic plots for CFR estimates
#################################################################################################

# -----------------------------------------------------------------------------------------------
# set up and load packages ----------------------------------------------------------------------
rm(list=ls())
list.files('FILEPATH')
run_date_03 <- 'RUN_DATE'

run_date <- run_date_03

library(tidyverse, lib.loc = "FILEPATH")
### load MSCM packages -- must be done in this order
library(crosswalk002, lib.loc = "FILEPATH")
library(mrbrt002, lib.loc = "FILEPATH")

library(msm)
library(readxl)
library(data.table)
library(arm)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
setDTthreads(1)

# -----------------------------------------------------------------------------------------------
# read covariates for prediction frame ----------------------------------------------------------

load(paste0('FILEPATH', run_date, '/pred_frame_with_uncertainty.csv'))

dir.create(paste0("FILEPATH",run_date,"/country_plots/"), recursive = T)

for(c in 1:length(unique(pred_frame2$ihme_loc_id))){
  
  Country_name <- unique(pred_frame2$ihme_loc_id)[c]
  message(Country_name)
  df_plot <- subset(pred_frame2, ihme_loc_id == Country_name)
  
  df_plot2 <- subset(df_plot, select = c("year_id", "Comm.ind", "age_start","age_end", "predicted_cfr"))
  
  df_plot2$comm_words <- ifelse(df_plot2$Comm.ind == 0, "Hospital", "Community")
  
  by_time <- ggplot(df_plot2,aes(x=year_id,y=predicted_cfr)) + labs(color="Age") +
    geom_line(aes(color=age_start, group=age_start),alpha=0.2 ) + theme_bw() + ggtitle(Country_name)+
    facet_wrap(~comm_words, scales = "free") + scale_color_viridis_c() + xlab("Year") + ylab("CFR")
  
  by_age <- ggplot(df_plot2,aes(x=age_start,y=predicted_cfr)) + labs(color="Year") +
    geom_line(aes(color=year_id, group=year_id) ) + theme_bw() + xlab("Age") + ylab("CFR") + 
    facet_wrap(~comm_words, scales = "free") + scale_color_viridis_c(option= 'magma', direction=-1)
  
  png(file = paste0("FILEPATH",run_date,"/country_plots/predictions_",Country_name,".png"),
      width = 12,
      height = 6,
      units = "in", 
      res = 300)
  plot(grid.arrange(by_time, by_age, nrow = 2))
  dev.off()
  
  
}


