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

###Read tidy hyena data tables
load('data/raw_data.RData')

source('00.define_functions.R')

##remove aggressions where recipient ignores or counterattacks
aggsFull <- filter(tblAggression, 
                   !response1 %in% c('ct', 'ignore'),
                   !response2 %in% c('ct', 'ignore'),
                   !response3 %in% c('ct', 'ignore'))

## Add year column
aggsFull$year <- format(aggsFull$date, '%Y')

## Add aggressor and recipient den end date and mother
aggsFull$agg.end <- left_join(aggsFull, cohortInfo, by = c('aggressor' = 'id'))$denend
aggsFull$recip.end <- left_join(aggsFull, cohortInfo, by = c('recip' = 'id'))$denend
aggsFull$agg.mom <- left_join(aggsFull, tblHyenas, by = c('aggressor' = 'id'))$mom
aggsFull$recip.mom <- left_join(aggsFull, tblHyenas, by = c('recip' = 'id'))$mom

#### Add aggressor and recipient rank from rank table. Rank calculated for each female yearly according to Strauss & Holekamp 2019
###use individual rank if individual is in ranks table (only true for adults)
aggsFull$agg.rank <- left_join(aggsFull, tblFemaleRanks, by = c('year', 'aggressor' = 'id'))$rank %>% as.numeric()
aggsFull$recip.rank <- left_join(aggsFull, tblFemaleRanks, by = c('year', 'recip' = 'id'))$rank %>% as.numeric()
###otherwise, use maternal rank by selecting rank of mother
aggsFull$agg.rank[is.na(aggsFull$agg.rank)] <- left_join(aggsFull[is.na(aggsFull$agg.rank),], tblFemaleRanks, by = c('year', 'agg.mom' = 'id'))$rank %>% as.numeric()
aggsFull$recip.rank[is.na(aggsFull$recip.rank)] <- left_join(aggsFull[is.na(aggsFull$recip.rank),], tblFemaleRanks, by = c('year', 'recip.mom' = 'id'))$rank %>% as.numeric()

###agg and recip birthdates. Use minimum of date of birth and date first seen to incorporate individuals whose DOB is missing
aggsFull$agg.dob <- pmin(left_join(aggsFull, tblLifeHistory.long, by = c('aggressor' = 'id'))$DOB,
                         left_join(aggsFull, tblLifeHistory.long, by = c('aggressor' = 'id'))$DFS)

aggsFull$recip.dob <- pmin(left_join(aggsFull, tblLifeHistory.long, by = c('recip' = 'id'))$DOB,
                         left_join(aggsFull, tblLifeHistory.long, by = c('recip' = 'id'))$DFS)


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
  end_diff <- elo_deviance(intx, id, k = 100)
  ## Observed elo score from same interactions
  end_obs <- elo_obs(intx, id, k = 100)
  
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
    postgrad_diff <- elo_deviance(intx, id, k = 100)
  }else(postgrad_diff <- NA)
  
  ####first year of adulthood - based on interactions between focal individuals and other adults
  intx <- filter(aggsFull, aggressor == id | recip == id, 
                 date > agg.dob+(365*2), date > recip.dob+(365*2),
                 date <= birthdate + (365*3),agg.rank != recip.rank)
  adult_intx <- nrow(intx)
  
  if(nrow(intx)){
    adult_diff <- elo_deviance(intx, id, k = 100)
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

save(cub_dev_vars, file = '09.cub_dev_vars_k100.RData')
