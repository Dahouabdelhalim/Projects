# ScarletFeverTSIR_RoyalSocOpen.R
# 
# Fit TSIR regression models and produce figures for article: 
#  McDonald et al., The dynamics of scarlet fever in the Netherlands, 1906-1920: a historical analysis. J Royal Society Open
#
# Requires downloading of Netherlands population/live births data from Statistics Netherlands, a meteorological dataset from KNMI,
#  and the scarlet fever case data per week for the Netherlands provided as a supplementary file to this paper.
#
# Last update: 20-10-2021
#



library(data.table)
library(aweek)
library(tsiR)
library(lubridate)


#wd <- ""            # set working directory to where data files are stored
setwd(wd)


end_yr <- 1920       # last year of analysis



########## 
### Assemble data frames

## Retrieve Netherlands national population size and number of births per year 'levendegeboren', 1906 through 1920
#   https://opendata.cbs.nl/statline/#/CBS/nl/dataset/37556/table?ts=1549969578873  export as 'csv met statistische symbolen'
# Select/re-arrange fields so that year, population size, and number of births are retrieved for the years 1906 through 1920, and exported as .csv 
vitalstats <- data.frame(read.csv2("<PATH TO YOUR FILE HERE>")
names(vitalstats) <- c("year","totalpop","livebirths") 
vitalstats[,c("totalpop")] <- vitalstats[,c("totalpop")]*1000 
vitalstats$birthrate <- vitalstats$livebirths/vitalstats$totalpop*1000


## Linearly interpolate population size of 3 largest cities for period 1906-1920, based on available data from 1900 and 1925
pop_3cities <- sum(c(511000,319000,206000))     # popul of Amsterdam, Rotterdam, Den Haag in 1900
prop_3cities <- pop_3cities/5197000             # 0.1993   prop. of entire 1900 population in 3 cities; assumed to apply throughout period
prop_city <- matrix(NA,nrow=3,ncol=(end_yr-1906+1))
nsteps <- (1925-1900+1)
tmp <- approx(c(511000,712200),n=nsteps)$y; prop_city[1,] <- tmp[7:(end_yr-1900+1)]/vitalstats$totalpop  # amsterdam 
tmp <- approx(c(319000,543700),n=nsteps)$y; prop_city[2,] <- tmp[7:(end_yr-1900+1)]/vitalstats$totalpop  # rotterdam 
tmp <- approx(c(206000,391400),n=nsteps)$y; prop_city[3,] <- tmp[7:(end_yr-1900+1)]/vitalstats$totalpop  # s gravenhage 
pop_city <- data.frame(urban=rep(c("amsterdam","rotterdam","s gravenhage","rest"),each=(end_yr-1906+1)),pop=c(prop_city[1,]*vitalstats$totalpop[1:(end_yr-1906+1)],
                        prop_city[2,]*vitalstats$totalpop[1:(end_yr-1906+1)],prop_city[3,]*vitalstats$totalpop[1:(end_yr-1906+1)],
                       (1-apply(prop_city[1:3,],2,sum))*vitalstats$totalpop[1:(end_yr-1906+1)]),year=rep(c(1906:end_yr),2)) 
rm(tmp,prop_city,pop_3cities,prop_3cities,nsteps)


## Load in meterological data from KNMI, 1902-1956, from https://daggegevens.knmi.nl/
# Need the following fields; easiest is to save a downloaded .txt file in Excel, add column headers (YYYYMMDD, TG, RH, UG), then save as .csv 
# TG  = Etmaalgemiddelde temperatuur (in 0.1 graden Celsius); 
# RH  = Etmaalsom van de neerslag (in 0.1 mm) (-1 voor <0.05 mm); 
# UG  = Etmaalgemiddelde relatieve vochtigheid (in procenten); 
#
#  Convert YYYYMMDD date to year + weeknumber; divide temp by 10 and neerslag by 10; 
tmp <- read.csv2("<PATH TO YOUR FILE HERE>",header=T)
#  Use aweek functions so can specify start date of week (to match case data week-53s (eg. 1904,1910,1916,1919)) 
tmp$Date <- as.Date(as.character(tmp$YYYYMMDD),"%Y%m%d")
tmp$Date.wk <- date2week(tmp$Date, week_start = 3)                    # 3 = wednesday
tmp$year <- as.numeric(as.character(substring(tmp$Date.wk, 1, 4)))    # as.character() to convert from aweek class
tmp$week <- as.numeric(as.character(substring(tmp$Date.wk, 7, 8)))

tmp.aggr <- (aggregate(tmp$TG,list(tmp$year,tmp$week),mean))          # calc weekly mean from daily avg temperature
names(tmp.aggr) <- c("year","week","Temp")
tmp.aggr$Temp <- tmp.aggr$Temp/10                                     # convert from 10ths of a degree

tmp.aggr2 <- (aggregate(tmp$RH,list(tmp$year,tmp$week),sum,na.rm=T))  # calc summed precip over week
names(tmp.aggr2) <- c("year","week","Precip")
tmp.aggr2$Precip <- ifelse(tmp.aggr2$Precip < 0,0,tmp.aggr2$Precip)   # convert -1 (indicating <0.05mm) to 0
tmp.aggr2$Precip <- tmp.aggr2$Precip/10                               # convert from tenths of 1mm to mm
tmp.aggr <- merge(tmp.aggr,tmp.aggr2,by=c("year","week"),all.x=T)

# Calc and replace relative humidity (field UG) with absolute humidity
tmp$UG <- (6.112 * exp((17.67 * (tmp$TG/10))/((tmp$TG/10)+243.5)) * (tmp$UG/10) * 2.1674) / (273.15+(tmp$TG/10))  # in grams/m3

tmp.aggr2 <- (aggregate(tmp$UG,list(tmp$year,tmp$week),sum,na.rm=T))   # calc weekly mean from daily avg humidity
names(tmp.aggr2) <- c("year","week","Humid")
tmp.aggr <- merge(tmp.aggr,tmp.aggr2,by=c("year","week"),all.x=T)
nrow(tmp.aggr)             # n=2870

# Aggregate rainfall by season
tmp$indicator_S <- ifelse(tmp$week > 9 & tmp$week < 23,1,0)    # 'meteorological season' spring only (1 March to 31 May); estimate using week 10 - week 22
tmp$indicator_W2 <- ifelse(tmp$week > 0 & tmp$week < 10,1,0)   # winter (1 Dec to 30 April); week 1 tm 9 and also week 49 - week 52/53 of prev year
tmp$indicator_W1 <- ifelse(tmp$week > 48 & tmp$week <= 53,1,0)
tmp.aggr3 <- (aggregate(tmp$RH,list(tmp$year,tmp$indicator_S),sum,na.rm=T))   # calc summed precip in spring of year
tmp.aggr3 <- subset(tmp.aggr3, Group.2==1);  tmp.aggr3$Group.2 <- NULL
names(tmp.aggr3) <- c("year","PrecipS")
tmp.aggr4 <- (aggregate(tmp$RH,list(tmp$year,tmp$indicator_W2),sum,na.rm=T))   # calc summed precip in late winter part of year y
tmp.aggr4 <- subset(tmp.aggr4, Group.2==1);  tmp.aggr4$Group.2 <- NULL
tmp.aggr5 <- (aggregate(tmp$RH,list(tmp$year,tmp$indicator_W1),sum,na.rm=T))   # calc summed precip in early winter part of year
tmp.aggr5$Group.1 <- (tmp.aggr5$Group.1 + 1)                                   # adjust year so summed rainfall in year y applies to (i.e., year + 1)
tmp.aggr5 <- subset(tmp.aggr5, Group.2==1);  tmp.aggr5$Group.2 <- NULL
tmp.aggr6 <- rbind(tmp.aggr4,tmp.aggr5); tmp.aggr6 <- aggregate(tmp.aggr6$x,list(tmp.aggr6$Group.1),sum,na.rm=T)
names(tmp.aggr6) <- c("year","PrecipW")
tmp.aggr <- merge(tmp.aggr,tmp.aggr3,by=c("year"),all.x=T)
tmp.aggr <- merge(tmp.aggr,tmp.aggr6,by=c("year"),all.x=T)

# Convert temp & abs humidity, using extreme quintiles over study period:
quintsTemp <- quantile(probs=c(0.2,0.8),tmp.aggr$Temp[210:991])
quintsHumid <- quantile(probs=c(0.2,0.8),tmp.aggr$Humid[210:991])
tmp.aggr$lowTemp <- ifelse(tmp.aggr$Temp < quintsTemp[1],1,0)
tmp.aggr$hiTemp <- ifelse(tmp.aggr$Temp > quintsTemp[2],1,0)
tmp.aggr$lowHumid <- ifelse(tmp.aggr$Humid < quintsHumid[1],1,0)
tmp.aggr$hiHumid <- ifelse(tmp.aggr$Humid > quintsHumid[2],1,0)

tmp.aggr <- tmp.aggr[ order(tmp.aggr$year,tmp.aggr$week), ]
KNMI <- tmp.aggr[,c("year","week","Temp","Precip","Humid","PrecipS","PrecipW","lowTemp","hiTemp","lowHumid","hiHumid")]
# Impute missing week 53 by simple extrapolation from week 52, for the year 1919 only
KNMI <- rbind(KNMI, c(1919,53,as.numeric(KNMI[KNMI$year==1919 & KNMI$week==52,c(3:ncol(KNMI))])))
rm(tmp,tmp.aggr,tmp.aggr2,tmp.aggr3,tmp.aggr4,tmp.aggr5,tmp.aggr6)
rm(quintsHumid,quintsTemp)


## Load in weekly scarlet fever data, 1906-1920, provided as a supplementary file to the paper
#  Columns are: c("year","week","termtime","amsterdam","rotterdam","s_gravenhage","rest")
notifdata <- read.csv2(file="<YOUR PATH>/Scarlet_fever_notifications_NL_1906-1920.csv")
maxweeks <- aggregate(notifdata$week,list(notifdata$year),max)          # calculate maximum number of weeknumbers per year
maxweeks$x <- ifelse(maxweeks$x < 52,52,maxweeks$x);  names(maxweeks) <- c("year","mx")


## Create yearweek data.frame
yearweek <- data.frame(year=rep(1906:end_yr,each=53),week=1:53)   # allow for up to 53 weeks per year
yearweek <- merge(yearweek,maxweeks,by="year")
yearweek <- subset(yearweek, week < 53 | (week==53 & mx==53))     # remove 'impossible' week 53s
yearweek$x <- yearweek$year + (yearweek$week/yearweek$mx)         # x-scale takes into account 52/53 week years

num_weeks <- nrow(yearweek); num_weeks   # n=783


# Merge in year-specific national popul size and birthrate
yearweek <- merge(yearweek,vitalstats[,c("totalpop","year","birthrate")],by=c("year"))  


# Merge in meteorological variables
yearweek <- merge(yearweek,KNMI[,c("year","week","Temp","Humid","Precip","PrecipS","PrecipW","lowTemp","hiTemp","lowHumid","hiHumid")],by=c("year","week"),all.x=T) 
yearweek <- yearweek[ order(yearweek$year,yearweek$week), ]   # sort


# Add linear ordinal Time indicator, starting with 1st week of data 
yearweek$time <- c(1:num_weeks)
# Add ordinal season indicator (= calendar year)
yearweek$season <- paste("s",yearweek$year-1906,sep="")


# Quadruple the number of rows in yearweek, merging in weekly number of scarlet fever cases per city, and termtime indicator
yearweek <- merge(yearweek,notifdata[,c("year","week","termtime","amsterdam")],by=c("year","week"))  # start with Amsterdam  
names(yearweek)[which(names(yearweek)=="amsterdam")] <- "scarletcases"
yearweek <- merge(yearweek,pop_city,by=c("year")); nrow(yearweek)             # quadruple size, by merging in pop_city   n=3132
yearweek <- yearweek[ order(yearweek$urban,yearweek$year,yearweek$week), ]    # and sort

# Merge in number of scarlet fever cases for remaining 2 cities/rest of Netherlands
yearweek$scarletcases[yearweek$urban=="rotterdam"] <- notifdata$rotterdam 
yearweek$scarletcases[yearweek$urban=="s gravenhage"] <- notifdata$s_gravenhage 
yearweek$scarletcases[yearweek$urban=="rest"] <- notifdata$rest 



## Create master week-level data.frame
data_TSIR <- data.frame(year = yearweek$year, week = yearweek$week, mx = yearweek$mx, Season = yearweek$season, Cases = yearweek$scarletcases, 
                        Time=yearweek$time, Birth_rate = yearweek$birthrate,
                        SchoolTerm = yearweek$termtime, 
                        PrecipS = yearweek$PrecipS, PrecipW = yearweek$PrecipW,
                        MeanWTemp = yearweek$Temp, MeanWHumid = yearweek$Humid,
                        lowTemp=yearweek$lowTemp, hiTemp=yearweek$hiTemp, lowHumid=yearweek$lowHumid, hiHumid=yearweek$hiHumid,
                        Urban = yearweek$urban, TotalPop = yearweek$totalpop, Pop = yearweek$pop,
                        Births = yearweek$birthrate/1000*yearweek$pop/yearweek$mx )  # estimate weekly births per stratum of Urban  

data_TSIR$Urban <- factor(data_TSIR$Urban)
nrow(data_TSIR)       # n=3132

data_TSIR_weekly <- data_TSIR    # copy, as data_TSIR will be aggregated to biweek



## Produce Fig. 3 - notif. rate per city per month, aggregating over all years, and mean of avg weekly temp & humidity per month
fig3 <- data_TSIR_weekly

# Convert year + week number to date (assume Monday), then work out month in which this Monday occurred
fig3$Date <- as.Date(paste(data_TSIR$year,data_TSIR$week,"1",sep=":"),format = "%Y:%W:%u")
fig3$MonthNum <- month(ymd(fig3$Date));  table(fig3$MonthNum)
fig3$Rate <- (fig3$Cases/fig3$Pop*100000)
MonthTotal <- as.data.frame(tapply(fig3$Rate,list(fig3$MonthNum,fig3$Urban),sum))
MonthTemp <- as.data.frame(tapply(fig3$MeanWTemp,list(fig3$MonthNum,fig3$Urban),mean))
MonthHumid <- as.data.frame(tapply(fig3$MeanWHumid,list(fig3$MonthNum,fig3$Urban),mean))

op <- par(family="sans",font=2,font.lab=2,font.axis=2,font.main=2,cex.axis=0.75,cex.lab=0.9,cex.main=0.8)
 # Divide by monthly total cases by 15 seasons
 matplot(1:12,MonthTotal$amsterdam/15,type="o",pch=15,lty=1,lwd=2,col="steelblue", xaxt="n",
         xlab="Month of year",ylab="Cases per 100,000",ylim=c(1,35))
 axis(side=1,labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),at=seq(1,12,by=1),tick=NA)
 lines(1:12,MonthTotal$rotterdam/15,type="o",lty=1,pch=19,lwd=2,col="green3")
 lines(1:12,MonthTotal$'s gravenhage'/15,type="o",lty=1,pch=17,lwd=2,col="orange2")
  par(new = TRUE)
  matplot(1:12,MonthTemp$amsterdam, type = "l",lty=2,lwd=1.6, col="tomato1", axes = FALSE, bty = "n", xlab = "", ylab = "")
  lines(1:12,MonthHumid$amsterdam*2,lty=2,lwd=1.6,col="grey")
  axis(side=4, at = pretty(range(MonthTemp$amsterdam)))
 legend(locator(1),legend=c("Amsterdam","Rotterdam","Den Haag","Temperature","Humidity"),cex=0.7,col=c("steelblue","green3","orange2","tomato1","grey"),
       lwd=c(2,2,2,1.6,1.6),lty=c(1,1,1,2,2),pch=c(15,19,17,NA,NA))
par(op)
rm(fig3,MonthTotal,MonthTemp,MonthHumid)





##########
###   Create weekly-level data.frame, used later for regression analyses and plots only
data_TSIR_weekly <- data_TSIR_weekly[,c("year","week","Urban","Cases","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")]



##########
###   Create biweek datafiles for TSIR model fitting in JAGS

## Set incubation period in weeks *AND* granularity of case data
IP <- 2 

CityResid <- list();  df.tsiR <- list()

# Amsterdam 
df.tsiR[[1]] <- data_TSIR[data_TSIR$Urban=="amsterdam",c("year","week","mx","Cases","Birth_rate","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")] 
# Rotterdam
df.tsiR[[2]] <- data_TSIR[data_TSIR$Urban=="rotterdam",c("year","week","mx","Cases","Birth_rate","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")]
# Den Haag
df.tsiR[[3]] <- data_TSIR[data_TSIR$Urban=="s gravenhage",c("year","week","mx","Cases","Birth_rate","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")]
# rest
df.tsiR[[4]] <- data_TSIR[data_TSIR$Urban=="rest",c("year","week","mx","Cases","Birth_rate","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")]

for(j in 1:4) {
  df.tsiR[[j]]$time <- (df.tsiR[[j]]$year + (df.tsiR[[j]]$week/df.tsiR[[j]]$mx)-0.015) 
  df.tsiR[[j]]$births <- (df.tsiR[[j]]$Birth_rate/df.tsiR[[j]]$mx/1000 * df.tsiR[[j]]$Pop) 
  df.tsiR[[j]]$year <- df.tsiR[[j]]$week <- NULL
  # NB. fields cases, pop, time and births need to be lowercase for upcoming tsiRdata() function...
  names(df.tsiR[[j]]) <- c("mx","cases","Birth_rate","pop", "SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS", "time","births")
}


## Aggregate cases & weather vars; recalc births to biweekly level for inference of weather/schoolterm covars using TSIR model in JAGS
for(j in 1:4) {
   tmp <- tsiRdata(df.tsiR[[j]]$time,df.tsiR[[j]]$cases,df.tsiR[[j]]$births,df.tsiR[[j]]$pop, IP=IP)  # this function does the aggregation
   # Add Birth_rate and mx to this biweek data.frame
   tmp$Birth_rate <- df.tsiR[[j]]$Birth_rate[seq(1, length(df.tsiR[[j]]$Birth_rate), IP)]
   tmp$mx <- df.tsiR[[j]]$mx[seq(1, length(df.tsiR[[j]]$mx), IP)]
   # As births are replaced by tsiRdata() assuming that births are total annual births (and so allocated to each week by dividing by 26), adjust:
   tmp$births <- (tmp$Birth_rate/tmp$mx*2/1000 * tmp$pop)
   # Aggregate Schoolterm by taking the value for the 1st week of the biweek
   tmp$SchoolTerm <- df.tsiR[[j]]$SchoolTerm[seq(1, nrow(df.tsiR[[j]]),IP)]
   # Average MeanTemp & MeanHumid over biweek; for cumul precip, take the value for the 1st week of the biweek
   tmp$MeanWTemp <- colMeans(rbind( df.tsiR[[j]]$MeanWTemp[seq(1, nrow(df.tsiR[[j]]),IP)], c(df.tsiR[[j]]$MeanWTemp[seq(2, nrow(df.tsiR[[j]]),IP)],
                                                                                           tail(df.tsiR[[j]]$MeanWTemp,1)) ))
   tmp$MeanWHumid <- colMeans(rbind( df.tsiR[[j]]$MeanWHumid [seq(1, nrow(df.tsiR[[j]]),IP)], c(df.tsiR[[j]]$MeanWHumid [seq(2, nrow(df.tsiR[[j]]),IP)],
                                                                                           tail(df.tsiR[[j]]$MeanWHumid ,1)) ))
   tmp$PrecipW <- df.tsiR[[j]]$PrecipW[seq(1, nrow(df.tsiR[[j]]),IP)]
   tmp$PrecipS <- df.tsiR[[j]]$PrecipS[seq(1, nrow(df.tsiR[[j]]),IP)]
   df.tsiR[[j]] <- tmp
}
## Re-make and re-order data_TSIR, to reflect biweek data and also including only necessary fields for TSIR model fitting in JAGS
#   NB. SchoolTerm & weeather vars are included as place holders as will get replaced later
tmp1 <- data_TSIR[data_TSIR$Urban=="amsterdam",c("Urban","Cases","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")][seq(1, nrow(data_TSIR[data_TSIR$Urban=="amsterdam",]), IP),]
tmp2 <- data_TSIR[data_TSIR$Urban=="rotterdam",c("Urban","Cases","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")][seq(1, nrow(data_TSIR[data_TSIR$Urban=="rotterdam",]), IP),]
tmp3 <- data_TSIR[data_TSIR$Urban=="s gravenhage",c("Urban","Cases","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")][seq(1, nrow(data_TSIR[data_TSIR$Urban=="s gravenhage",]), IP),]
tmp4 <- data_TSIR[data_TSIR$Urban=="rest",c("Urban","Cases","Pop","SchoolTerm","MeanWTemp","MeanWHumid","PrecipW","PrecipS")][seq(1, nrow(data_TSIR[data_TSIR$Urban=="rest",]), IP),]
data_TSIR <- as.data.frame(rbind(tmp1,tmp2,tmp3,tmp4))
rm(tmp,tmp1,tmp2,tmp3,tmp4)


## Finally, remove extra fields
for(j in 1:4) {
 df.tsiR[[j]]$mx <- df.tsiR[[j]]$Birth_rate <- NULL
}


## Est. susceptibles using balance equation, and extract rho, Z and Sbar
Sbar <- list()

for(j in 1:4) {
  parms <- estpars(df.tsiR[[j]],IP=IP,regtype="gaussian",family='quasipoisson',link="log")
  names(parms)
  print(summary(parms$rho))  
  df.tsiR[[j]]$I <- (parms$rho * df.tsiR[[j]]$cases)   
  # Easiest to lag Z and logI here:
  df.tsiR[[j]]$lagZ <- c( head(parms$Z,1),head(parms$Z,-1) )     # take n-1, with 1st element duplicated
  df.tsiR[[j]]$laglogI <- c( log(head((parms$rho * df.tsiR[[j]]$cases),1)+1),log(head((parms$rho * df.tsiR[[j]]$cases),-1)+1) )
  Sbar[[j]] <- parms$sbar      # Sbar is time-invariant
}


## Integrate these new vars into data_TSIR, replacing Pop and Births:
#  Also replace Cases, and school-term & weather vars, as are aggregated over biweek
jj <- 0
for(j in c("amsterdam","rotterdam","s gravenhage","rest")) {
 jj <- jj + 1
 data_TSIR$Cases[data_TSIR$Urban==j] <- df.tsiR[[jj]]$cases
 data_TSIR$lagZ[data_TSIR$Urban==j] <- df.tsiR[[jj]]$lagZ
 data_TSIR$Cases_adj[data_TSIR$Urban==j] <- df.tsiR[[jj]]$I
 data_TSIR$laglogI[data_TSIR$Urban==j] <- df.tsiR[[jj]]$laglogI
 data_TSIR$Pop[data_TSIR$Urban==j] <- df.tsiR[[jj]]$pop
 data_TSIR$SchoolTerm[data_TSIR$Urban==j] <- df.tsiR[[jj]]$SchoolTerm
 data_TSIR$MeanWTemp[data_TSIR$Urban==j] <- df.tsiR[[jj]]$MeanWTemp
 data_TSIR$MeanWHumid[data_TSIR$Urban==j] <- df.tsiR[[jj]]$MeanWHumid
 data_TSIR$PrecipW[data_TSIR$Urban==j] <- df.tsiR[[jj]]$PrecipW
 data_TSIR$PrecipS[data_TSIR$Urban==j] <- df.tsiR[[jj]]$PrecipS
}


## Recalculate quintiles of Temp & Humid and replace hi/lowTemp and hi/lowHumid accordingly
quintsTemp <- quantile(probs=c(0.2,0.8),data_TSIR$MeanWTemp[data_TSIR$Urban==j])  # doesn't matter what j is set to, as these covars are same for all cities
quintsHumid <- quantile(probs=c(0.2,0.8),data_TSIR$MeanWHumid[data_TSIR$Urban==j])
data_TSIR$lowTemp <- ifelse(data_TSIR$MeanWTemp < quintsTemp[1],1,0)
data_TSIR$hiTemp <- ifelse(data_TSIR$MeanWTemp > quintsTemp[2],1,0)
data_TSIR$lowHumid <- ifelse(data_TSIR$MeanWHumid < quintsHumid[1],1,0)
data_TSIR$hiHumid <- ifelse(data_TSIR$MeanWHumid > quintsHumid[2],1,0)



### Save data file suitable for TSIR model fitting using separate JAGS code :: go to script ScarletFeverTSIR_JAGS_RoyalSocOpen.R:
save(data_TSIR,data_TSIR_weekly,Sbar,file="data_TSIR_Sbar_scarlet_fever.Rdata")


