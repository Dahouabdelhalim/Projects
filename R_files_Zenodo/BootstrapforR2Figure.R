#### Bootstrap for Figure 2 of "Partisan Pandemic" 
#### Authors: Clinton, Cohen, Lapinski, and Trussler
#### Prepared: November 2020
rm(list=ls())

#Packages
library(stargazer)
library(lfe)
library(Hmisc)

#Set Working Directory

setwd("C:/Users/Marc/Dropbox/SM Covid Science Paper")

load("SM Data for Science Final.Rdata")


dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))


boot <- 1000

party.se <- NA
party.up <- NA
party.down <- NA
party.med <- NA

health.se <- NA
health.up <- NA
health.down <- NA
health.med <- NA

dat <- dat[!is.na(dat$last24hr_index) & !is.na(dat$dem) & !is.na(dat$rep) & !is.na(dat$case.change),]
pb <- txtProgressBar(min = 0, max = length(dates), style = 3)


for (i in 1:length(dates)){
  dat.day <- dat[dat$SMDate==dates[i],]
  party.boot <- rep(NA, boot)
  health.boot <- rep(NA, boot)
      for(j in 1:boot){
          dat.boot <- dat.day[sample(1:nrow(dat.day), replace=T),]
         
 
             m <- felm(last24hr_index ~ 
                        female + income_30_75 + income_75_150 + income_above150 +
                        age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                        other_race + black + asian + hispanic + 
                        hs + associates_somecollege + college + postgrad + 
                        unemployed + density| state | 0 | state, data=dat.boot,
                         weights=dat.boot$weight_national_18plus_daily)
          
            
            m.p <- felm(last24hr_index ~  dem + rep +
                          female + income_30_75 + income_75_150 + income_above150 +
                          age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                          other_race + black + asian + hispanic + 
                          hs + associates_somecollege + college + postgrad + 
                          unemployed  + density| state | 0 | state, data=dat.boot,
                        weights=dat.boot$weight_national_18plus_daily)
            
 
            
            m.h <- felm(last24hr_index ~ case.change +
                          female + income_30_75 + income_75_150 + income_above150 +
                          age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                          other_race + black + asian + hispanic + 
                          hs + associates_somecollege + college + postgrad + 
                          unemployed  + density| state | 0 | state, data=dat.boot,
                        weights=dat.boot$weight_national_18plus_daily)
            
            party.boot[j] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
            health.boot[j] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
        
          
       
        }
  party.se[i] <- sd(party.boot, na.rm=T)
  party.up[i] <- quantile(party.boot, .975, na.rm=T)
  party.down[i] <- quantile(party.boot, .025, na.rm=T)
  party.med[i] <- median(party.boot, na.rm=T)
  
  health.se[i] <- sd(health.boot, na.rm=T)
  health.up[i] <- quantile(health.boot, .975, na.rm=T)
  health.down[i] <- quantile(health.boot, .025, na.rm=T)
  health.med[i] <- median(health.boot, na.rm=T)
  Sys.sleep(0.1)
  # update progress bar
  setTxtProgressBar(pb, i)  
  
}



save(party.se, party.up, party.down,party.med, health.se, health.up, health.down, health.med, file="Bootstrap Estimates.Rdata")

