#### Analysis File for "Partisan Pandemic
#### Authors: Clinton, Cohen, Lapinski, and Trussler
#### Prepared: November 2020
rm(list=ls())

#Packages
library(stargazer)
library(lfe)
library(Hmisc)

#Set Working Directory
setwd("C:/Users/Marc/Dropbox/SM Covid Science Paper/Replication")

#Load dataset
load("SM Data for Science Final.Rdata")

#Define colors

col2rgb("firebrick4")
firebrickfaded <- rgb(139, 26, 26, max = 255, alpha = 100, names = "firebrick4")

col2rgb("dodgerblue4")
dodgerbluefaded <- rgb(16, 78, 139, max = 255, alpha = 100, names = "dodgerblue4")

col2rgb("goldenrod3")
goldenrodfaded <- rgb(205, 155, 29, max = 255, alpha = 100, names = "goldenrod3")

col2rgb("springgreen4")
springgreenfaded <- rgb(0, 139, 69, max = 255, alpha = 100, names = "springgreen4")

col2rgb("purple")
purplefaded <- rgb(160, 32, 240, max = 255, alpha = 100, names = "purple")

viridis1 <- rgb(68,1,84, max = 255)
viridis1faded <- rgb(68,1,84, max = 255, alpha = 100)

viridis2 <- rgb(49,104,142, max = 255)
viridis2faded <- rgb(49,104,142, max = 255, alpha = 100)

viridis3 <- rgb(53,183,121,max = 255)
viridis3faded <- rgb(53,183,121, max = 255, alpha = 100)

viridis4 <- rgb(180,222,44, max = 255)
viridis4faded <- rgb(180,222,44, max = 255, alpha = 100)


######################################################
###################Current Numbers####################
######################################################
#This section reproduces various numbers we cite in the manuscript

#Number of interviews for Covid Concern
print(table(is.na(dat$is.worried))[1])

#Number of interviews for last 24
print(table(is.na(dat$last24hr_index))[1])

#Number of interviews for social distancing
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
mean(table(dat$SMDate[dat$SMDate > min(dates)]))

# First Day of SD questions
min(dates)

# Last Date
max(dates)

#Overall, how much more to reps distance compared to dem? Percent change from ind?

m <- lm(last24hr_index ~ dem + rep, dat, weights =dat$weight_national_18plus_daily)
coef(m)["rep"] - coef(m)["dem"]

(coef(m)["dem"]/coef(m)["(Intercept)"])*100
(coef(m)["rep"]/coef(m)["(Intercept)"])*100

#On first day, how much more to reps distance compared to dem
min.date <- min(dat$SMDate[!is.na(dat$last24hr_index)])

m <- lm(last24hr_index ~ dem + rep, dat[dat$SMDate==min.date,],
        weights =dat$weight_national_18plus_daily[dat$SMDate==min.date])
ind1 <- coef(m)["(Intercept)"]
rep1 <- coef(m)["(Intercept)"] + coef(m)["rep"]
dem1 <- coef(m)["(Intercept)"] + coef(m)["dem"]
coef(m)["rep"] - coef(m)["dem"]

#On last day, how much more to reps distance compared to dem
max.date <- max(dat$SMDate[!is.na(dat$last24hr_index)])
m <- lm(last24hr_index ~ dem + rep, dat[dat$SMDate==max.date,],
        weights =dat$weight_national_18plus_daily[dat$SMDate==max.date])
coef(m)["rep"] - coef(m)["dem"]
ind2 <- coef(m)["(Intercept)"]
rep2 <- coef(m)["(Intercept)"] + coef(m)["rep"]
dem2 <- coef(m)["(Intercept)"] + coef(m)["dem"]


rep2 - rep1
(rep2 - rep1)/rep1

ind2 - ind1
(ind2 - ind1)/ind1

dem2 - dem1
(dem2 - dem1)/dem1


#On first day, how much less worried reps compared to dems (proportions)
min.date <- min(dat$SMDate[!is.na(dat$is.worried)])

m <- lm(is.worried ~ dem + rep, dat[dat$SMDate==min.date,],
        weights =dat$weight_national_18plus_daily[dat$SMDate==min.date])
coef(m)["rep"] - coef(m)["dem"]

#On last day,  how much less worried reps compared to dems (proportions)
max.date <- max(dat$SMDate[!is.na(dat$is.worried)])
m <- lm(is.worried ~ dem + rep, dat[dat$SMDate==max.date,],
        weights =dat$weight_national_18plus_daily[dat$SMDate==max.date])
coef(m)["rep"] - coef(m)["dem"]


#Percentage of Dems and Reps across daily samples
dates <- sort(unique(dat$SMDate[!is.na(dat$dem)]))
dem.avg <- rep(NA, length(dates))
rep.avg <- rep(NA, length(dates))

for(i in 1:length(dates)){
  dem.avg[i] <- weighted.mean(dat$dem[dat$SMDate==dates[i]], dat$weight_national_18plus_daily[dat$SMDate==dates[i]],na.rm=T)
  rep.avg[i] <- weighted.mean(dat$rep[dat$SMDate==dates[i]], dat$weight_national_18plus_daily[dat$SMDate==dates[i]],na.rm=T)
}

mean(dem.avg)
mean(rep.avg)

##Pooled Regression for overall coefficient between dem/rep and ind

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

# Zipcode fixed effects


m.zip <- felm(last24hr_index ~time.trend + time.trend2+ dem + rep + case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad
          | zip_code | 0 |  zip_code, data=dat,
          weights=dat$weight_national_18plus_daily)


summary(m.zip)

# State fixed effects

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

dat$real.cases <- dat$cases*dat$population/10000
summary(dat$real.cases)

m <- felm(last24hr_index ~time.trend + time.trend2+ dem + rep + case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state | 0 |  state , data=dat,
          weights=dat$weight_national_18plus_daily)
getfe(m)
summary(m)
coef(m)["dem"]
coef(m)["rep"]



#Summary of how much more variation party explains over cases


m <- felm(last24hr_index ~ time.trend + time.trend2 +
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density| state | 0 | state, data=dat,
          weights=dat$weight_national_18plus_daily)



m.p <- felm(last24hr_index ~   time.trend + time.trend2 +dem + rep +
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| state | 0 | state , data=dat,
            weights=dat$weight_national_18plus_daily)



m.h <- felm(last24hr_index ~  time.trend + time.trend2 +case.change +
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| state | 0 | state , data=dat,
            weights=dat$weight_national_18plus_daily)


summary(m)
summary(m.p)
summary(m.h)

party.var <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
health.var <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)

party.var/health.var




#Summary of how much more variation party explains over cases, state level


m <- felm(last24hr_index ~ time.trend + time.trend2 +
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density| state | 0 | state, data=dat,
          weights=dat$weight_national_18plus_daily)



m.p <- felm(last24hr_index ~   time.trend + time.trend2 +dem + rep +
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| state | 0 | state, data=dat,
            weights=dat$weight_national_18plus_daily)



m.h <- felm(last24hr_index ~  time.trend + time.trend2 +state.case.change +
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| state | 0 | state, data=dat,
            weights=dat$weight_national_18plus_daily)


summary(m)
summary(m.p)
summary(m.h)

party.var <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
health.var <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)

party.var/health.var



######################################################
#############Aggregate Trends in SAH##################
######################################################
#This section reproduces Figure 1

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

dem.avg.24 <- NA
rep.avg.24 <- NA
ind.avg.24 <- NA

dem.se.24 <- NA
rep.se.24 <- NA
ind.se.24 <- NA

one.avg.24 <- NA
two.avg.24 <- NA
three.avg.24 <- NA
four.avg.24 <- NA

one.se.24 <- NA
two.se.24 <- NA
three.se.24 <- NA
four.se.24 <- NA


for(i in 1:length(dates)){
  dem.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],na.rm=T)
  rep.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],na.rm=T)
  ind.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$pid5=="Independent"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$pid5=="Independent"],na.rm=T)
  dem.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")]))
  rep.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")]))
  ind.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$pid5=="Independent"],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$pid5=="Independent"],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$pid5=="Independent"]))
  
  
  
  one.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$covid.quartile==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==1],na.rm=T)
  two.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$covid.quartile==2],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==2],na.rm=T)
  three.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$covid.quartile==3],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==3],na.rm=T)
  four.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$covid.quartile==4],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==4],na.rm=T)
  
  one.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==1],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==1]))
  
  two.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==2],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==2],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==2]))
  
  three.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==3],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==3],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==3]))
  four.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==4],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==4],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$covid.quartile==4]))
  
}


#   Raw Change in Social distancing
#   For TEXT

n <- length(dem.avg.24)

dem.avg.24[1]
dem.avg.24[n]
dem.avg.24[n]-dem.avg.24[1]
(dem.avg.24[n]-dem.avg.24[1])/dem.avg.24[1]

rep.avg.24[1]
rep.avg.24[n]
rep.avg.24[n] - rep.avg.24[1]
(rep.avg.24[n]-rep.avg.24[1])/rep.avg.24[1]

ind.avg.24[1]
ind.avg.24[n]
ind.avg.24[n] - ind.avg.24[1]
(ind.avg.24[n] - ind.avg.24[1])/ind.avg.24[1]

#Loess estimation last24
dem.fit.24 <- predict(loess(dem.avg.24 ~ as.numeric(dates), span=.4, weights = 1/dem.se.24 ), se=T)
rep.fit.24 <- predict(loess(rep.avg.24 ~ as.numeric(dates), span=.4, weights = 1/rep.se.24 ), se=T)
ind.fit.24 <- predict(loess(ind.avg.24 ~ as.numeric(dates), span=.4, weights = 1/ind.se.24), se=T)


one.fit.24 <- predict(loess(one.avg.24 ~ as.numeric(dates), span=.4, weights = 1/one.se.24), se=T)
two.fit.24 <- predict(loess(two.avg.24 ~ as.numeric(dates), span=.4, weights = 1/two.se.24), se=T)
three.fit.24 <- predict(loess(three.avg.24 ~ as.numeric(dates), span=.4, weights = 1/three.se.24), se=T)
four.fit.24 <- predict(loess(four.avg.24 ~ as.numeric(dates), span=.4, weights = 1/four.se.24), se=T)


dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))

dem.avg.wor <- NA
rep.avg.wor <- NA
ind.avg.wor <- NA

dem.se.wor <- NA
rep.se.wor <- NA
ind.se.wor <- NA

one.avg.wor <- NA
two.avg.wor <- NA
three.avg.wor <- NA
four.avg.wor <- NA

one.se.wor <- NA
two.se.wor <- NA
three.se.wor <- NA
four.se.wor <- NA


for(i in 1:length(dates)){
  dem.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],na.rm=T)
  rep.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],na.rm=T)
  ind.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$pid5=="Independent"],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$pid5=="Independent"],na.rm=T)
  dem.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Democrat" | dat$pid5=="Lean Democrat")]))
  rep.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & (dat$pid5=="Republican" | dat$pid5=="Lean Republican")]))
  ind.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$pid5=="Independent"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$pid5=="Independent"],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$pid5=="Independent"]))
  
  
  
  one.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$covid.quartile==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==1],na.rm=T)
  two.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$covid.quartile==2],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==2],na.rm=T)
  three.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$covid.quartile==3],
                                   w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==3],na.rm=T)
  four.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$covid.quartile==4],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$covid.quartile==4],na.rm=T)
  
  one.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==1]))
  
  two.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==2],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==2],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==2]))
  
  three.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==3],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==3],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==3]))
  four.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==4],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$covid.quartile==4],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$covid.quartile==4]))
  
}

#   Raw Change in Social distancing
#   For TEXT
n <- length(dem.avg.wor)
dem.avg.wor[n]
rep.avg.wor[n]



#Loess estimation worried
dem.fit.wor <- predict(loess(dem.avg.wor ~ as.numeric(dates), span=.4, weights = 1/dem.se.wor ), se=T)
rep.fit.wor <- predict(loess(rep.avg.wor ~ as.numeric(dates), span=.4, weights = 1/rep.se.wor ), se=T)
ind.fit.wor <- predict(loess(ind.avg.wor ~ as.numeric(dates), span=.4, weights = 1/ind.se.wor), se=T)


#For the day where everyone had same answer, impute previous day
one.se.wor[9] <- one.se.wor[8]
one.avg.wor[9] <- one.avg.wor[8]



one.fit.wor <- predict(loess(one.avg.wor ~ as.numeric(dates), span=.4, weights = 1/one.se.wor), se=T)
two.fit.wor <- predict(loess(two.avg.wor ~ as.numeric(dates), span=.4, weights = 1/two.se.wor), se=T)
three.fit.wor <- predict(loess(three.avg.wor ~ as.numeric(dates), span=.4, weights = 1/three.se.wor), se=T)
four.fit.wor <- predict(loess(four.avg.wor ~ as.numeric(dates), span=.4, weights = 1/four.se.wor), se=T)

#Generate figure

#png(file="Figures/Figure1.png", height=1152, width=1152)
pdf(file="Figures/Figure1.pdf", height=12, width=12)
par(mfrow=c(2,2))
dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))
plot(dates, dem.avg.wor, type="p", lwd=2, col=dodgerbluefaded, ylim=c(0,1),xlim=as.Date(c("2020-02-11","2020-10-01")),
     ylab="Proportion Very or Somewhat Worried about COVID-19", xlab="Date", pch=16, axes=F,main="A",cex=.5)
segments(dates,dem.avg.wor, dates,dem.avg.wor+dem.se.wor*1.96, col=dodgerbluefaded)
segments(dates,dem.avg.wor, dates,dem.avg.wor-dem.se.wor*1.96, col=dodgerbluefaded)

lines(dates,dem.fit.wor$fit,col="dodgerblue4", lwd=2)
lines(dates,dem.fit.wor$fit - 1.96*dem.fit.wor$se.fit, col="dodgerblue4", lwd=1, lty=2)
lines(dates,dem.fit.wor$fit + 1.96*dem.fit.wor$se.fit, col="dodgerblue4", lwd=1, lty=2)


points(dates,rep.avg.wor, type="p", lwd=2, col=firebrickfaded, pch=16,cex=.5)
segments(dates,rep.avg.wor, dates,rep.avg.wor+rep.se.wor*1.96, col=firebrickfaded)
segments(dates,rep.avg.wor, dates,rep.avg.wor-rep.se.wor*1.96, col=firebrickfaded)

lines(dates,rep.fit.wor$fit,col="firebrick4", lwd=2)
lines(dates,rep.fit.wor$fit - 1.96*rep.fit.wor$se.fit, col="firebrick4", lwd=1, lty=2)
lines(dates,rep.fit.wor$fit + 1.96*rep.fit.wor$se.fit, col="firebrick4", lwd=1, lty=2)


points(dates, ind.avg.wor, type="p", col=purplefaded, lwd=2, pch=16,cex=.5)
segments(dates,ind.avg.wor, dates,ind.avg.wor+ind.se.wor*1.96, col=purplefaded)
segments(dates,ind.avg.wor, dates,ind.avg.wor-ind.se.wor*1.96, col=purplefaded)

lines(dates,ind.fit.wor$fit,col="purple", lwd=2)
lines(dates,ind.fit.wor$fit - 1.96*ind.fit.wor$se.fit, col="purple", lwd=1, lty=2)
lines(dates,ind.fit.wor$fit + 1.96*ind.fit.wor$se.fit, col="purple", lwd=1, lty=2)

legend("topleft", c("Democrat",  "Independent",  "Republican"), pch=16,
       col=c("dodgerblue4", "purple", "firebrick4"), bg="white")


axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")

axis(2, at=seq(0,1,.25), las=2)


plot(dates, one.avg.wor,type="p", lwd=2, col=viridis1faded, ylim=c(0,1),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Proportion Very or Somewhat Worried about COVID-19",
     xlab="Date", pch=16, axes=F,main="B",cex=.5)
segments(dates,one.avg.wor, dates,one.avg.wor+one.se.wor*1.96, col=viridis1faded)
segments(dates,one.avg.wor, dates,one.avg.wor-one.se.wor*1.96, col=viridis1faded)



points(dates, two.avg.wor,type="p",col=viridis2faded , pch=16, cex=.5)
segments(dates,two.avg.wor, dates,two.avg.wor+two.se.wor*1.96, col=viridis2faded)
segments(dates,two.avg.wor, dates,two.avg.wor-two.se.wor*1.96, col=viridis2faded)


points(dates, three.avg.wor,type="p",col=viridis3faded , pch=16, cex=.5)
segments(dates,three.avg.wor, dates,three.avg.wor+three.se.wor*1.96, col=viridis3faded)
segments(dates,three.avg.wor, dates,three.avg.wor-three.se.wor*1.96, col=viridis3faded)

points(dates, four.avg.wor,type="p",col=viridis4faded , pch=16, cex=.5)
segments(dates,four.avg.wor, dates,four.avg.wor+four.se.wor*1.96, col=viridis4faded)
segments(dates,four.avg.wor, dates,four.avg.wor-four.se.wor*1.96, col=viridis4faded)


lines(dates,one.fit.wor$fit,col=viridis1, lwd=2)
lines(dates,one.fit.wor$fit - 1.96*one.fit.wor$se.fit, col=viridis1, lwd=1, lty=2)
lines(dates,one.fit.wor$fit + 1.96*one.fit.wor$se.fit, col=viridis1, lwd=1, lty=2)


lines(dates,two.fit.wor$fit,col=viridis2, lwd=2)
lines(dates,two.fit.wor$fit - 1.96*two.fit.wor$se.fit, col=viridis2, lwd=1, lty=2)
lines(dates,two.fit.wor$fit + 1.96*two.fit.wor$se.fit, col=viridis2, lwd=1, lty=2)

lines(dates,three.fit.wor$fit,col=viridis3, lwd=2)
lines(dates,three.fit.wor$fit - 1.96*three.fit.wor$se.fit, col=viridis3, lwd=1, lty=2)
lines(dates,three.fit.wor$fit + 1.96*three.fit.wor$se.fit, col=viridis3, lwd=1, lty=2)


lines(dates,four.fit.wor$fit,col=viridis4, lwd=2)
lines(dates,four.fit.wor$fit - 1.96*four.fit.wor$se.fit, col=viridis4, lwd=1, lty=2)
lines(dates,four.fit.wor$fit + 1.96*four.fit.wor$se.fit, col=viridis4, lwd=1, lty=2)


legend("topleft", c("COVID-19 Severity Quartile 1 (Low)", "COVID-19 Severity Quartile 2", "COVID-19 Severity Quartile 3", "COVID-19 Severity Quartile 4 (High)"),
       pch=rep(16,4), col=c(viridis1,viridis2,viridis3,viridis4), bg="white")


axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")

axis(2, at=seq(0,1,0.25), las=2)
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
plot(dates, dem.avg.24, type="p", lwd=2, col=dodgerbluefaded, ylim=c(1,3.5),xlim=as.Date(c("2020-02-11","2020-10-01")),
     ylab="Number of Activities in Previous 24hr", xlab="Date", pch=16, axes=F,main="C",cex=.5)
segments(dates,dem.avg.24, dates,dem.avg.24+dem.se.24*1.96, col=dodgerbluefaded)
segments(dates,dem.avg.24, dates,dem.avg.24-dem.se.24*1.96, col=dodgerbluefaded)

lines(dates,dem.fit.24$fit,col="dodgerblue4", lwd=2)
lines(dates,dem.fit.24$fit - 1.96*dem.fit.24$se.fit, col="dodgerblue4", lwd=1, lty=2)
lines(dates,dem.fit.24$fit + 1.96*dem.fit.24$se.fit, col="dodgerblue4", lwd=1, lty=2)


points(dates,rep.avg.24, type="p", lwd=2, col=firebrickfaded, pch=16,cex=.5)
segments(dates,rep.avg.24, dates,rep.avg.24+rep.se.24*1.96, col=firebrickfaded)
segments(dates,rep.avg.24, dates,rep.avg.24-rep.se.24*1.96, col=firebrickfaded)

lines(dates,rep.fit.24$fit,col="firebrick4", lwd=2)
lines(dates,rep.fit.24$fit - 1.96*rep.fit.24$se.fit, col="firebrick4", lwd=1, lty=2)
lines(dates,rep.fit.24$fit + 1.96*rep.fit.24$se.fit, col="firebrick4", lwd=1, lty=2)


points(dates, ind.avg.24, type="p", col=purplefaded, lwd=2, pch=16,cex=.5)
segments(dates,ind.avg.24, dates,ind.avg.24+ind.se.24*1.96, col=purplefaded)
segments(dates,ind.avg.24, dates,ind.avg.24-ind.se.24*1.96, col=purplefaded)

lines(dates,ind.fit.24$fit,col="purple", lwd=2)
lines(dates,ind.fit.24$fit - 1.96*ind.fit.24$se.fit, col="purple", lwd=1, lty=2)
lines(dates,ind.fit.24$fit + 1.96*ind.fit.24$se.fit, col="purple", lwd=1, lty=2)

legend("topleft", c("Democrat",  "Independent",  "Republican"), pch=16,
       col=c("dodgerblue4", "purple", "firebrick4"), bg="white")


axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")

axis(2, at=seq(1,3.5,.5), las=2)


plot(dates, one.avg.24,type="p", lwd=2, col=viridis1faded, ylim=c(1,3.5),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Number of Activities in Previous 24hr",
     xlab="Date", pch=16, axes=F,main="D",cex=.5)
segments(dates,one.avg.24, dates,one.avg.24+one.se.24*1.96, col=viridis1faded)
segments(dates,one.avg.24, dates,one.avg.24-one.se.24*1.96, col=viridis1faded)



points(dates, two.avg.24,type="p",col=viridis2faded , pch=16, cex=.5)
segments(dates,two.avg.24, dates,two.avg.24+two.se.24*1.96, col=viridis2faded)
segments(dates,two.avg.24, dates,two.avg.24-two.se.24*1.96, col=viridis2faded)


points(dates, three.avg.24,type="p",col=viridis3faded , pch=16, cex=.5)
segments(dates,three.avg.24, dates,three.avg.24+three.se.24*1.96, col=viridis3faded)
segments(dates,three.avg.24, dates,three.avg.24-three.se.24*1.96, col=viridis3faded)

points(dates, four.avg.24,type="p",col=viridis4faded , pch=16, cex=.5)
segments(dates,four.avg.24, dates,four.avg.24+four.se.24*1.96, col=viridis4faded)
segments(dates,four.avg.24, dates,four.avg.24-four.se.24*1.96, col=viridis4faded)


lines(dates,one.fit.24$fit,col=viridis1, lwd=2)
lines(dates,one.fit.24$fit - 1.96*one.fit.24$se.fit, col=viridis1, lwd=1, lty=2)
lines(dates,one.fit.24$fit + 1.96*one.fit.24$se.fit, col=viridis1, lwd=1, lty=2)


lines(dates,two.fit.24$fit,col=viridis2, lwd=2)
lines(dates,two.fit.24$fit - 1.96*two.fit.24$se.fit, col=viridis2, lwd=1, lty=2)
lines(dates,two.fit.24$fit + 1.96*two.fit.24$se.fit, col=viridis2, lwd=1, lty=2)

lines(dates,three.fit.24$fit,col=viridis3, lwd=2)
lines(dates,three.fit.24$fit - 1.96*three.fit.24$se.fit, col=viridis3, lwd=1, lty=2)
lines(dates,three.fit.24$fit + 1.96*three.fit.24$se.fit, col=viridis3, lwd=1, lty=2)


lines(dates,four.fit.24$fit,col=viridis4, lwd=2)
lines(dates,four.fit.24$fit - 1.96*four.fit.24$se.fit, col=viridis4, lwd=1, lty=2)
lines(dates,four.fit.24$fit + 1.96*four.fit.24$se.fit, col=viridis4, lwd=1, lty=2)


legend("topleft", c("COVID-19 Severity Quartile 1 (Low)", "COVID-19 Severity Quartile 2", "COVID-19 Severity Quartile 3", "COVID-19 Severity Quartile 4 (High)"),
       pch=rep(16,4), col=c(viridis1,viridis2,viridis3,viridis4), bg="white")


axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")

axis(2, at=seq(1,3.5,.5), las=2)


dev.off()
par(mfrow=c(1,1))


######################################################
#############  Changes Across Time####################
######################################################
#This section reproces Figures 2 & 3
#Note for Figure 2 bootstrap standard errors are used, these estimates are produced in an external file, 
#and the results are loaded in below. Note that the creation of the bootstrap standard errors is extremely 
#computationally intensive 

#In order to have same n for all specifications, must not be NA for all three of the key variables

dat <- dat[!is.na(dat$last24hr_index) & !is.na(dat$dem) & !is.na(dat$rep) & !is.na(dat$death.change),]
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))

party.var <- NA
health.var <- NA

dem.coef <- NA
rep.coef <- NA

case.coef <- NA

dem.se <- NA
rep.se <- NA

case.se <- NA


for (i in 1:length(dates)){
  m <- felm(last24hr_index ~ 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| state | 0 | state, data=dat[dat$SMDate==dates[i],],
            weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.p <- felm(last24hr_index ~  dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.h <- felm(last24hr_index ~ case.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.var[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.var[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  m.f <- felm(last24hr_index ~ dem + rep + case.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  dem.coef[i] <- coef(m.f)["dem"]
  rep.coef[i] <- coef(m.f)["rep"]
  
  dem.se[i] <- sqrt(vcov(m.f)["dem","dem"])
  rep.se[i] <- sqrt(vcov(m.f)["rep","rep"])
  
  case.coef[i] <- coef(m.f)["case.change"]
  case.se[i] <- sqrt(vcov(m.f)["case.change","case.change"])
  
}

#Above produces the "prime" estimates, external file bootstraps the above (1000 BS samples per day) for the variance explained figure
load(file="Bootstrap Estimates.Rdata")

# Predicted Effects
#   For TEXT

n <- length(dem.coef)
dem.coef[n]
rep.coef[n]

(dem.coef[n] - dem.coef[1])/dem.coef[1]
(rep.coef[n] - rep.coef[1])/rep.coef[1]

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates), span=.4, weights = 1/party.se*100), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates),weights = 1/health.se*100), se=T)


#png("Figures/Figure2.png", height=576, width=768)
pdf("Figures/Figure2.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)
segments(dates, party.med*100, dates, party.up*100, col=goldenrodfaded)
segments(dates, party.med*100, dates, party.down*100, col=goldenrodfaded)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)
segments(dates, health.med*100, dates, health.up*100, col=springgreenfaded)
segments(dates, health.med*100, dates, health.down*100, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "Change in County Cases/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()




#Loess estimation
dem.coef.fit <- predict(loess(dem.coef ~ as.numeric(dates), span=.4, weights = 1/dem.se), se=T)

#Loess estimation
rep.coef.fit <- predict(loess(rep.coef ~ as.numeric(dates), span=.4, weights = 1/rep.se), se=T)

#Loess estimation
case.coef.fit <- predict(loess(case.coef ~ as.numeric(dates), span=.4, weights = 1/case.se), se=T)


#All-in-one figure

#png(file="Figures/Figure3.png", height=576, width=768)
pdf(file="Figures/Figure3.pdf", height=6, width=8)
plot(dates, case.coef, pch=16, xlab="", ylab="Daily Coefficient on 24 Hour Activity", col=springgreenfaded, ylim=c(-.6,1.1),
     axes=F)
abline(h=0,lty=2,col="gray60")
segments(dates, case.coef,dates, case.coef+1.96*case.se, col=springgreenfaded)
segments(dates, case.coef,dates, case.coef-1.96*case.se, col=springgreenfaded)

lines(dates,case.coef.fit$fit, col="springgreen4", lwd=2)
lines(dates,case.coef.fit$fit - 1.96*case.coef.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,case.coef.fit$fit + 1.96*case.coef.fit$se.fit, col="springgreen4", lwd=1, lty=2)

points(dates, rep.coef, pch=16, col=firebrickfaded)
abline(h=0,lty=2,col="gray60",cex=.5)
segments(dates, rep.coef,dates, rep.coef+1.96*rep.se, col=firebrickfaded)
segments(dates, rep.coef,dates, rep.coef-1.96*rep.se, col=firebrickfaded)


lines(dates,rep.coef.fit$fit, col="firebrick4", lwd=2)
lines(dates,rep.coef.fit$fit - 1.96*rep.coef.fit$se.fit, col="firebrick4", lwd=1, lty=2)
lines(dates,rep.coef.fit$fit + 1.96*rep.coef.fit$se.fit, col="firebrick4", lwd=1, lty=2)

points(dates, dem.coef, pch=16,col=dodgerbluefaded)
abline(h=0,lty=2,col="gray60",cex=.5)
segments(dates, dem.coef,dates, dem.coef+1.96*dem.se, col=dodgerbluefaded)
segments(dates, dem.coef,dates, dem.coef-1.96*dem.se, col=dodgerbluefaded)


lines(dates,dem.coef.fit$fit, col="dodgerblue4", lwd=2)
lines(dates,dem.coef.fit$fit - 1.96*dem.coef.fit$se.fit, col="dodgerblue4", lwd=1, lty=2)
lines(dates,dem.coef.fit$fit + 1.96*dem.coef.fit$se.fit, col="dodgerblue4", lwd=1, lty=2)

legend("topleft",c("Democrat", "Republican","Change in Cases/1000"), pch=c(16,16,16), col=c("dodgerblue4", "firebrick4","springgreen4"), bg="white")

axis(2,at=seq(-.5, 1, .2),las=2)
axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
dev.off()



######################################################
###############State Differences######################
######################################################
#This code recreates figure 4

states <- unique(dat$state)

dem.coef <-rep(NA, length(states))
rep.coef <-rep(NA, length(states))
dem.se <-  rep(NA, length(states))
rep.se <-  rep(NA, length(states))
party.var <- rep(NA, length(states))
health.med <- rep(NA, length(states))

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

for(j in 1:length(states)){
  
  m <- felm(last24hr_index ~time.trend +time.trend2 +
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density, data=dat[dat$state==states[j],],
            weights=dat$weight_state_weekly[dat$state==states[j]])
  summary(m)
  
  
  m.p <- felm(last24hr_index ~time.trend +time.trend2 + dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density, data=dat[dat$state==states[j],],
              weights=dat$weight_state_weekly[dat$state==states[j]])
  summary(m.p)
  
  dem.coef[j] <- coef(m.p)["dem"]
  rep.coef[j] <- coef(m.p)["rep"]
  dem.se[j] <-   vcov(m.p)["dem","dem"]
  rep.se[j] <-    vcov(m.p)["rep","rep"]
  
  
  
  m.h <- felm(last24hr_index ~time.trend +time.trend2 + case.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density, data=dat[dat$state==states[j],],
              weights=dat$weight_state_weekly[dat$state==states[j]])
  summary(m.h)
  
  party.var[j] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[j] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)  
  
  
  m.f <- felm(last24hr_index ~time.trend +time.trend2 + dem + rep + case.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density, data=dat[dat$state==states[j],],
              weights=dat$weight_state_weekly[dat$state==states[j]])
  summary(m.f)
  
  dem.coef[j] <- coef(m.f)["dem"]
  rep.coef[j] <- coef(m.f)["rep"]
  dem.se[j] <-   vcov(m.f)["dem","dem"]
  rep.se[j] <-    vcov(m.f)["rep","rep"]
}


state.data <- cbind.data.frame(states,dem.coef,dem.se,rep.coef,rep.se,party.var,health.med)


state.data$states <- as.character(state.data$states)
state.data <- state.data[order(state.data$states, decreasing = T),]


#By State Mitigation Scores

load("StateMitigation.Rdata")

state.data <- merge(state.data, st.m, by.x="states", by.y="state.full")

head(state.data)

state.data.dem <- state.data[state.data$gov.party=="dem",]
state.data.rep <- state.data[state.data$gov.party=="rep",]



#Loess estimation
state.data.dem <- state.data.dem[order(state.data.dem$mitigation.scores),]
dem.coef.fit <- predict(loess(state.data.dem$dem.coef ~ state.data.dem$mitigation.scores,  weights = 1/state.data.dem$dem.se), se=T)
rep.coef.fit <- predict(loess(state.data.dem$rep.coef ~ state.data.dem$mitigation.scores,  weights = 1/state.data.dem$rep.se), se=T)


state.data.rep <- state.data.rep[order(state.data.rep$mitigation.scores),]
dem.coef.fit <- predict(loess(state.data.rep$dem.coef ~ state.data.rep$mitigation.scores,  weights = 1/state.data.rep$dem.se), se=T)
rep.coef.fit <- predict(loess(state.data.rep$rep.coef ~ state.data.rep$mitigation.scores,  weights = 1/state.data.rep$rep.se), se=T)




pdf(file="Figures/Figure4.pdf", height=12, width=10)
par(mfrow=c(2,1))

plot(state.data.dem$mitigation.scores, state.data.dem$dem.coef, pch=16, col=dodgerbluefaded, xlab="State Covid Mitigation Score",
     ylab="Party Coefficient in State (difference from Independent)", xlim=c(-2.5,2.5), ylim=c(-0.6,0.8), axes=F,main="States with Democratic Governors", col.main="dodgerblue",
     type="n")
segments(state.data.dem$mitigation.scores, state.data.dem$dem.coef,state.data.dem$mitigation.scores, state.data.dem$dem.coef+ 1.96*state.data.dem$dem.se, col=dodgerbluefaded)
segments(state.data.dem$mitigation.scores, state.data.dem$dem.coef,state.data.dem$mitigation.scores, state.data.dem$dem.coef- 1.96*state.data.dem$dem.se, col=dodgerbluefaded)
text(state.data.dem$mitigation.scores, state.data.dem$dem.coef, labels = state.data.dem$state.abbr, col="dodgerblue")


segments(state.data.dem$mitigation.scores, state.data.dem$rep.coef,state.data.dem$mitigation.scores, state.data.dem$rep.coef+ 1.96*state.data.dem$rep.se, col=firebrickfaded)
segments(state.data.dem$mitigation.scores, state.data.dem$rep.coef,state.data.dem$mitigation.scores, state.data.dem$rep.coef- 1.96*state.data.dem$rep.se, col=firebrickfaded)
text(state.data.dem$mitigation.scores, state.data.dem$rep.coef, labels = state.data.dem$state.abbr, col="firebrick")

axis(side=1, at=seq(-2.5,2.5,.5))
axis(side=2, at=seq(-0.6,0.8,0.2), labels = c("-0.6","-0.4","-0.2","0", "0.2","0.4","0.6","0.8"), las=2)

text(c(-1.75,1.75), c(0.03,0.03), labels=c("Less Aggressive COVID Mitigation","More Aggresive COVID Mitigation"))

abline(h=0, lty=2, col="gray60")




plot(state.data.rep$mitigation.scores, state.data.rep$dem.coef, pch=16, col=dodgerbluefaded, xlab="State Covid Mitigation Score",
     ylab="Party Coefficient in State (difference from Independent)", xlim=c(-2.5,2.5), ylim=c(-0.6,0.8), axes=F,main="States with Republican Governors", col.main="firebrick",
     type="n")
segments(state.data.rep$mitigation.scores, state.data.rep$dem.coef,state.data.rep$mitigation.scores, state.data.rep$dem.coef+ 1.96*state.data.rep$dem.se, col=dodgerbluefaded)
segments(state.data.rep$mitigation.scores, state.data.rep$dem.coef,state.data.rep$mitigation.scores, state.data.rep$dem.coef- 1.96*state.data.rep$dem.se, col=dodgerbluefaded)
text(state.data.rep$mitigation.scores, state.data.rep$dem.coef, labels = state.data.rep$state.abbr, col="dodgerblue")


segments(state.data.rep$mitigation.scores, state.data.rep$rep.coef,state.data.rep$mitigation.scores, state.data.rep$rep.coef+ 1.96*state.data.rep$rep.se, col=firebrickfaded)
segments(state.data.rep$mitigation.scores, state.data.rep$rep.coef,state.data.rep$mitigation.scores, state.data.rep$rep.coef- 1.96*state.data.rep$rep.se, col=firebrickfaded)
text(state.data.rep$mitigation.scores, state.data.rep$rep.coef, labels = state.data.rep$state.abbr, col="firebrick")

axis(side=1, at=seq(-2.5,2.5,.5))
axis(side=2, at=seq(-0.6,0.8,0.2), labels = c("-0.6","-0.4","-0.2","0", "0.2","0.4","0.6","0.8"), las=2)
text(c(-1.75,1.75), c(0.03,0.03), labels=c("Less Aggressive COVID Mitigation","More Aggresive COVID Mitigation"))

abline(h=0, lty=2, col="gray60")
dev.off()



