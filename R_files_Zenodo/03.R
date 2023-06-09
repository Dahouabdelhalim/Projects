################################################################################
#                        Calculate Elo-deviance scores                         #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################

######Load packages and set global options
rm(list = ls())
library(dplyr)
library(aniDom)
library(here)
options(stringsAsFactors = FALSE)

load('02.cohortInfo.RData')

###Read data files
load('data/raw_data.RData')

source('00.define_functions.R')

#### Assemble relevant variables for each cub. 
cub_dev_vars <- list()
for(i in 1:nrow(cohortInfo)){
  if(i %% 100 == 0){
    cat('Working on row ', i, 'of ', nrow(cohortInfo), '\\n')
  }
  
  ####Basics
  id <- cohortInfo$id[i]
  mom <- cohortInfo$mom[i]
  
  birthdate <- cohortInfo$dob[i]
  sex <- cohortInfo$sex[i]
  den_period <- as.numeric(cohortInfo$denend[i]) - as.numeric(cohortInfo$dob[i])+1
  clan <- cohortInfo$clan[i]
  
  disappeared <- cohortInfo$disappeared[i]
  last.seen.sessions <- max(filter(tblHyenasPerSession, id == cohortInfo$id[i])$date, na.rm = TRUE)
  
  
  ###Deviance scores
  ## Filter aggressions such that id is either aggressor or recipient and takes place only among den-dwelling cubs. Aggressor rank and recipient rank can't be the same i.e. no sibs
  intx <- filter(aggsFull, aggressor == id | recip == id, date <= agg.end, date <= recip.end,
                 agg.rank != recip.rank)
  
  ## skip the rest if there are no interactions from which to calculate deviance scores
  if(nrow(intx) < 1){next}
  
  ## Deviance score
  end_diff <- elo_deviance(intx, id)
  ## Observed elo score from same interactions
  end_obs <- elo_obs(intx, id)
  
  ## Number of interactions
  num_intx = nrow(intx)

  ## Mothers rank. If interactions occur in multiple years, use mean of rank from both years.
  mom_rank <- mean(as.numeric(filter(tblFemaleRanks, year %in% intx$year, id == cohortInfo$mom[i])$stan_rank), na.rm = T)
  
  
  ####post grad deviance score - based on interactions between focal individuals and subadults
  intx <- filter(aggsFull, aggressor == id | recip == id, date > agg.end, date > recip.end, 
                 date <= agg.dob+(365*2), date <= recip.dob+(365*2),
                 agg.rank != recip.rank)
  postgrad_intx = nrow(intx)
  if(nrow(intx)){
    postgrad_diff <- elo_deviance(intx, id)
  }else(postgrad_diff <- NA)
  
  ####first year of adulthood - based on interactions between focal individuals and other adults
  intx <- filter(aggsFull, aggressor == id | recip == id, 
                 date > agg.dob+(365*2), date > recip.dob+(365*2),
                 date <= birthdate + (365*3),agg.rank != recip.rank)
  adult_intx <- nrow(intx)
  
  if(nrow(intx)){
    adult_diff <- elo_deviance(intx, id)
  }else(adult_diff <- NA)

  ### Does the mother survive to graduation or adulthood?  
  mom_survive_to_grad <- ifelse(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared >= as.Date(cohortInfo$denend[i] , origin = '1970-01-01') | is.na(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared),1 , 0) %>%
    ifelse(length(.), ., NA)
  mom_survive_to_2 <- ifelse(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared >= as.Date(cohortInfo$dob[i] + 365*2 , origin = '1970-01-01') | is.na(filter(tblLifeHistory.long, id == cohortInfo$mom[i])$Disappeared),1 , 0) %>%
    ifelse(length(.), ., NA)
  
  
  if(id %in% tblFemaleRanks$id){
    first_adult_rank <- as.numeric(filter(tblFemaleRanks, id == cohortInfo$id[i])[1,]$stan_rank)
    year.joined <- filter(tblFemaleRanks, id == cohortInfo$id[i])[1,]$year
    if(mom %in% filter(tblFemaleRanks, year == year.joined)$id){
      ranks <- filter(tblFemaleRanks, year == year.joined, clan == cohortInfo$clan[i])
      if(!is.na(as.numeric(cohortInfo$litrank[i]))){
        expected_spot <- as.numeric(ranks[ranks$id == mom,'rank']) + as.numeric(cohortInfo$litrank[i])
      }else{expected_spot <- as.numeric(ranks[ranks$id == mom,'rank']) + 1}
      adult_first_rank_diff <- round(first_adult_rank - 
                                       ((-2*(expected_spot-1)/(nrow(ranks)-1))+1), 4)
      adult_expected_first_rank <- round(((-2*(expected_spot-1)/(nrow(ranks)-1))+1), 4)
    }else{adult_first_rank_diff <- NA}
  }else{first_adult_rank <- NA; adult_first_rank_diff <- NA; adult_expected_first_rank <- NA}
  
  ## Survival to den independence
  survive_to_grad = cohortInfo$survive_to_grad[i]
  ## Lifetime reproductive success
  lrs = cohortInfo$lrs[i]
  
  ##record summary statistics
  cub_dev_vars[[i]] <- data.frame(id, mom, sex, birthdate, den_period, num_intx,
                                  disappeared, last.seen.sessions,
                                  end_diff, postgrad_diff, adult_diff,mom_rank,
                                  postgrad_intx, adult_intx,end_obs,lrs,
                                  mom_survive_to_grad, mom_survive_to_2,
                                  survive_to_grad, adult_expected_first_rank, first_adult_rank,
                                  adult_first_rank_diff, clan)

}
cub_dev_vars <- do.call(rbind, cub_dev_vars)

save(cub_dev_vars, file = '04.cub_dev_vars.RData')
