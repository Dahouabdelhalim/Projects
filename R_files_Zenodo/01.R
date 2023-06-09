################################################################################
#                   Identify den cubs and den dependent period                 #
#                                                                              #
#                                 Eli Strauss                                  #
#                                                                              #
#                                                                              #
#                               February 2019                                  #
################################################################################

######Load packages and set global options
rm(list = ls())
library(dplyr)
library(here)
options(stringsAsFactors = FALSE)

####Read and tidy data
load(file = 'data/raw_data.RData')

##############Assemble data#######################
### Select hyenas with a known birthdate
cohortInfo <- filter(left_join(tblHyenas, tblLifeHistory.long, by = 'id'), !is.na(DOB))
names(cohortInfo) <- tolower(names(cohortInfo))

### Add clan information
cohortInfo$clan <-  left_join(cohortInfo, filter(tblLifeHistory, event_code == 'DOB'), by = 'id')$event_data

## Restrict to four clans with the best long-term data
cohortInfo <- filter(cohortInfo, clan %in% c('talek', 'serena.n', 'serena.s', 'happy.zebra'))

## Restrict to known sex, seen in tblHyenasPerSession, not involved in experimental manipulation at den
cub.fight.ids <- intersect(unique(c(filter(aggsFull, date >= '2015-01-01',
                                           context == 'experiment')$aggressor,
                                    filter(aggsFull, date >= '2015-01-01',
                                           context == 'experiment')$recip)),
                           filter(tblLifeHistory.long, DOB >= '2014-06-01')$id)
cohortInfo <- filter(cohortInfo, sex %in% c('m', 'f'), id %in% tblHyenasPerSession$id,
                     !id %in% cub.fight.ids, !is.na(mom))

### Restrict to cubs who were born at least 1 year before the end of the aggression data
cohortInfo <- filter(cohortInfo, (clan %in% c('talek', 'happy.zebra') & dob <= (as.Date('2016-06-19') - 365)) |
                       (clan == 'serena.n' & dob <= (as.Date('2016-12-31') - 365)) |
                       (clan == 'serena.s' & dob <= (as.Date('2017-03-05') - 365)))



## Add den graduation dates, binary variable for if cub survived to graduation, disappeared date, and date last observed, lrs
cohortInfo$survive_to_grad <- 1
cohortInfo$denend <- NA
cohortInfo$lrs <- NA
for(row in 1:nrow(cohortInfo)){
  if(is.na(cohortInfo[row,'dengrad'])){
    if(is.na(cohortInfo[row,'disappeared'])){
      cohortInfo[row,'denend'] <- cohortInfo[row,'dob']+365 ## If no graduation date and hyena is still alive, use 1 year
    }else{
      if(cohortInfo[row,]$dob+365 < cohortInfo[row,]$disappeared){
        cohortInfo[row,'denend'] <- cohortInfo[row,]$dob+365 ## IF no grad date and hyena died after 1 yr, use 1 yr
      }else{
        cohortInfo[row,'denend'] <- cohortInfo[row,]$disappeared ## If no grad date and hyena died 1 yr or before, den end = disappeared date and survive to grad = 0
        cohortInfo[row,]$survive_to_grad <- 0
      }
    }
  }else{cohortInfo[row,'denend'] <- min(as.numeric(cohortInfo[row,'dengrad']),  as.numeric(cohortInfo[row,'dob']+365), na.rm = T) } # Restrict >1yr grad dates to 1yr
  
  cohortInfo$last.seen.sessions <- max(filter(tblHyenasPerSession, id == cohortInfo$id[row])$date, na.rm = TRUE)
  
  
  ## LRS - number of offspring that survived to at least 2 years old
  if(!is.na(cohortInfo[row,]$disappeared) & cohortInfo[row,]$sex == 'f'){
    kids <- filter(tblHyenas, mom == cohortInfo$id[row])$id
    cohortInfo$lrs[row] <- nrow(filter(tblLifeHistory.long, id %in% kids, 
                                   is.na(Disappeared) | Disappeared - DOB >= 365*2))
  }
}

save(cohortInfo, file = '02.cohortInfo.RData')
