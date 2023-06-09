

### Get prey data
##### Function for calculating prey density - prey density = # of animals/km^2
##### time.period = the number of days relative to the date for calculating
##### prey density. 
##### Negative numbers = before the date, positive numbers = after the date.
##### For example, to get the prey density for 110 days prior to the date, set ##### time.period = -110
get_prey_density <- function(dates, clans, time.period, prey.animals = c("thomsons", "impala", "zebra", "wildebeest", "topi", "warthog", "hartebeest", "grants","buffalo", "hippo", "giraffe", "ostrich", "eland", "elephant", "oribi", "reedbuck","waterbuck", "baboon", "bushbuck")){
  
  ## Check to make sure supplied information is of same length
  if(length(clans) == 1){
    clans <- rep(clans, length(dates))
    warning('only 1 clan provided. recycling clan for each date')
  }else if (length(clans) != length(dates)){
    stop('clans not the same lenght as dates')
  }
  
  ## Store prey densities in here
  preys <- rep(NA, length(dates))
  
  for(i in 1:length(dates)){
    
    ## Select prey censuses from the appropriate clan associated with 110 days prior to birth
    prey.count.entries <- filter(tblPreyCensus, clan == clans[i] & date < dates[i] & date >= (dates[i] + time.period))[c('distance', prey.animals)]
    
    ## If no prey census data, leave as NA
    if(nrow(prey.count.entries) < 1)
      next
    
    ## Convert entries to numeric
    prey.count.entries[,c('distance', prey.animals)] <- sapply(prey.count.entries[,c('distance', prey.animals)], as.numeric)
    
    ## Count all species observed over the period (date + time.period), divided by distance. Also divide by 0.2 because prey censuses on 100m on either side of census route (i.e., 200m/1000m = 0.2)
    preys[i] <- sum(prey.count.entries[,prey.animals], na.rm = TRUE)/sum(prey.count.entries[,'distance'] * 0.2)
  }
  return(preys)
}

### Get number of cubs
#### This function extracts the number of cubs a hyena was in contact with
#### in the past month. 

get_cub_associates <- function(ids, dates, time.period){
  
  cub.associates <- rep(NA, length(ids))
  
  for(i in 1:length(ids)){
    if(!ids[i] %in% tblHyenas$id)
      next
    associates <- filter(tblHyenasPerSession, 
                         session %in% filter(tblHyenasPerSession, id %in% ids[i])$session,
                         date <= dates[i],
                         date >= dates[i] + time.period)
    associates$age <- associates$date - left_join(associates, tblHyenas, by = 'id')$birthdate
    associates <- filter(associates, age <= 365)
    cubs.per.session <- associates %>%
      group_by(session) %>% 
      summarise(cubs.per.session = length(id), .groups = 'drop_last')
    
    cub.associates[i] <- mean(cubs.per.session$cubs.per.session, na.rm = TRUE)
  }
  return(cub.associates)
}

model_summary <- function(fit){
  
  model.summary <- summary(fit)
  df <- round(model.summary$fixed, 4)
  rownames(df) <- gsub(rownames(df), pattern = 'mu', replacement = '')
  
  df.output <- df[,c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI')]
  df.diagnostic <- df[,c('Rhat', 'Bulk_ESS', 'Tail_ESS')]
  
  cat('MODEL SPECIFICATION:\\n\\n')
  cat('Family: ', model.summary$formula$family$family, '\\n')
  cat('Formula: ', capture.output(model.summary$formula), '\\n')
  cat('Number of observations: ', model.summary$nobs, '\\n')
  cat('Samples: ', model.summary$chains, 'chains, each with iter = ', model.summary$iter, '; warmup = ', model.summary$warmup, '; thin = ', model.summary$thin, '\\n')
  cat('\\nPRIORS:\\n\\n')
  print(prior_summary(fit, all = TRUE), show_df = FALSE)
  cat('\\nMODEL OUTPUT:\\n')
  print(df.output)
  cat('\\nMODEL DIAGNOSTICS:\\n')
  print(df.diagnostic)
}
