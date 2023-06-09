#### Analysis File for "Partisan Pandemic" Supplemental Material
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
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))

######################################################
###########Number of Interviews per day################
######################################################

n.resp <- rep(NA, length(dates))

for(i in 1: length(dates)){
  n.resp[i] <- length(dat$last24hr_index[dat$SMDate==dates[i] & !is.na(dat$last24hr_index)])
}

barplot(n.resp)

#png(file="Figures/SM1.png", height=560, width=768)
pdf(file="Figures/S1.pdf", height=6, width=8)
plot(dates,n.resp, pch=16, xlab="",ylab="Number of Respondents", axes=F, ylim=c(500,15000), type="n")
segments(dates,0,dates,n.resp, lwd=5)
axis(2, at=seq(1000,15000,2000), las=2)
axis(1, at=seq(1,100,1))
axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
legend("topright", paste("Total = ", print(sum(n.resp))), bg="white")
dev.off()



######################################################
########Activity averages over time#################
######################################################


#

names(dat)
vars <- names(dat)[which(names(dat)=="last24hr_mc_go_to_work") : which(names(dat)=="last24hr_mc_medical_care") ]
labels <- c("Go to Work","Groceries","Exercise","Take a Walk","Visit Family/Friend","Go to Resturant","Seek Medical Care")
class(dat$last24hr_mc_go_to_work)

#png(file="Figures/SMAllActivitiesTrends.png", height=1728, width=1536)
pdf(file="Figures/S2.pdf", height=18, width=16)

par(mfrow=c(3,3))
for(j in 1:7){
  dat$dv <- dat[,vars[j]]
  
  
  means <- rep(NA, length(dates))
  ses <- rep(NA, length(dates))
  for (i in 1:length(dates)){
    
    means[i] <- weighted.mean(dat$dv[dat$SMDate==dates[i]], dat$weight_national_18plus_daily[dat$SMDate==dates[i]],na.rm=T)
    ses[i]<- sqrt(wtd.var(dat$dv[dat$SMDate==dates[i]],
                          w=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],na.rm=T))/
      sqrt(length(dat$dv[dat$SMDate==dates[i]]))
  }
  
  act.fit <- predict(loess(means ~ as.numeric(dates), weights = 1/ses), se=T)
  
  
  plot(dates, means, main=paste(labels[j]),
       ylab="Proportion Participating",xlab="Date" , ylim=c(0,1), pch=16, col="gray60", axes=F)
  segments(dates, means,dates,means+1.96*ses, col="gray60")
  segments(dates, means,dates,means-1.96*ses, col="gray60")
  
  
  
  lines(dates,act.fit$fit, lwd=2)
  lines(dates,act.fit$fit - 1.96*act.fit$se.fit,  lwd=1, lty=2)
  lines(dates,act.fit$fit + 1.96*act.fit$se.fit,  lwd=1, lty=2)
  
  
  
  axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-09-09") , by="1 week"), format="%m-%d")
  axis(2, at=seq(0,1,0.25), las=2)
  legend("topleft", c(paste("Mean = ", print(round(weighted.mean(dat$dv,dat$weight_national_18plus_daily,na.rm=T),2)))))
}
dev.off()



######################################################
########Validation against Google Data#################
######################################################


goog <- read.csv("Global_Mobility_Report.csv")

goog <- goog[goog$country_region=="United States",]

blank <- goog$sub_region_1[1]

goog <- goog[goog$sub_region_1==blank,]

#Add one day, to match what the respondents are reporting
goog$date2 <- as.Date(goog$date) + 1

#Create national summary
dat <- dat[!is.na(dat$weight_national_18plus_daily),]

date <- unique(dat$SMDate[!is.na(dat$last24hr_index)])

avg.activity <- rep(NA, length(date))

for(i in 1:length(date)){
  avg.activity[i] <- weighted.mean(dat$last24hr_index[dat$SMDate==date[i]],
                                   dat$weight_national_18plus_daily[dat$SMDate==date[i]],na.rm=T)  
}


activity.dat <- cbind.data.frame(date,avg.activity)

goog <- merge(goog, activity.dat, by.x="date2", by.y="date")


pdf(file="Figures/S3.pdf", height = 12, width=8)
par(mfrow=c(3,2))
plot(goog$retail_and_recreation_percent_change_from_baseline, goog$avg.activity, pch=16, ylab="24hr Activity Estimated from Survey",
     xlab="Mobility Data: Retail and Recreation", main="Retail and Recreation")
legend("topleft", paste("Correlation:", print(round(cor(goog$retail_and_recreation_percent_change_from_baseline, goog$avg.activity),2))), bg="white")


plot(goog$grocery_and_pharmacy_percent_change_from_baseline, goog$avg.activity, pch=16, ylab="24hr Activity Estimated from Survey",
     xlab="Mobility Data: Grocery and Pharmacy", main="Grocery and Pharmacy")
legend("topleft", paste("Correlation:", print(round(cor(goog$grocery_and_pharmacy_percent_change_from_baseline, goog$avg.activity),2))), bg="white")

plot(goog$workplaces_percent_change_from_baseline, goog$avg.activity, pch=16, ylab="24hr Activity Estimated from Survey",
     xlab="Mobility Data: Workplace", main="Workplace")
legend("topleft", paste("Correlation:", print(round(cor(goog$workplaces_percent_change_from_baseline, goog$avg.activity),2))), bg="white")


plot(goog$transit_stations_percent_change_from_baseline, goog$avg.activity, pch=16, ylab="24hr Activity Estimated from Survey",
     xlab="Mobility Data: Transit Stations", main="Transit Stations")
legend("topleft", paste("Correlation:", print(round(cor(goog$transit_stations_percent_change_from_baseline, goog$avg.activity),2))), bg="white")


plot(goog$parks_percent_change_from_baseline, goog$avg.activity, pch=16, ylab="24hr Activity Estimated from Survey",
     xlab="Mobility Data: Parks", main= "Parks")
legend("topleft", paste("Correlation:", print(round(cor(goog$parks_percent_change_from_baselin, goog$avg.activity),2))), bg="white")
dev.off()


######################################################
##################Pooled Model########################
######################################################
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

m <- felm(last24hr_index ~time.trend + time.trend2+ pid5 + case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state| 0 | state , data=dat,
          weights=dat$weight_state_weekly)
summary(m)

stargazer(m, type="text")

stargazer(m,
          star.char = c("*"),
          star.cutoffs = c(0.05),
          style="ajps",
          digits=2,
          covariate.labels = c( "Time","Time2",
                                "Democrat","Lean Democrat", "Lean Republican", "Republican",
                               "County Change in Covid Cases/1000",
                               "Female",
                               "Income 30k-74K","Income 75k-149k","Income 150k+",
                               "Age 25-34", "Age 34-44", "Age 45-54","Age 55-64","Age 65+",
                                "Race Other","Black","Asian","Hispanic",
                                "High School","Associates/Some College","College","Postgraduate",
                               "Unemployed","Rural","Urban"),
          omit= c(),
          omit.stat = c("rsq"),
          dep.var.labels   = "Last 24 Hour Activity (0-7)",
          model.numbers = T,
          add.lines = list(c("State F.E.", "Yes")),
          notes.append = FALSE,
          title= "Determinants of Social Distancing Behavior",
          label="pooled",
          type="html",
          out="Figures/TableS1.html")



dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

#Zip FE

m <- felm(last24hr_index ~time.trend + time.trend2+ pid5 + case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed
          | zip_code | 0 | zip_code , data=dat,
          weights=dat$weight_state_weekly)
summary(m)



stargazer(m,
          star.char = c("*"),
          star.cutoffs = c(0.05),
          style="ajps",
          digits=2,
          covariate.labels = c("Democrat","Lean Democrat", "Lean Republican", "Republican",
                               "County Change in Covid Cases/10000"),
          omit= c("time.trend","time.trend2","female","income_30_75","income_75_150","income_above150",
                  "age_25_34","age_35_44","age_45_54","age_55_64","age_65_p",
                  "other_race","black","asian","hispanic","hs","associates_somecollege","college","postgrad","unemployed","density"),
          omit.stat = c("rsq"),
          dep.var.labels   = "Last 24 Hour Activity (0-7)",
          model.numbers = T,
          add.lines = list(c("State F.E.", "Yes"),
                           c("Demographics", "Yes"),
                           c("Polynomial time Trend", "Yes")),
          notes.append = FALSE,
          title= "Determinants of Social Distancing Behavior",
          label="pooled",
          type="html",
          out="Figures/TableS2.html")

#State case change

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

m.1 <- felm(last24hr_index ~time.trend + time.trend2+ 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state| 0 | state , data=dat,
          weights=dat$weight_state_weekly)
summary(m.1)


m.2 <- felm(last24hr_index ~time.trend + time.trend2+ dem + rep  + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state| 0 | state , data=dat,
          weights=dat$weight_state_weekly)
summary(m.2)


m.3 <- felm(last24hr_index ~time.trend + time.trend2+ state.case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state| 0 | state , data=dat,
          weights=dat$weight_state_weekly)
summary(m.3)


m.4 <- felm(last24hr_index ~time.trend + time.trend2+ dem + rep + state.case.change + 
            female + income_30_75 + income_75_150 + income_above150 +
            age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
            other_race + black + asian + hispanic + 
            hs + associates_somecollege + college + postgrad + 
            unemployed + density
          | state| 0 | state , data=dat,
          weights=dat$weight_state_weekly)
summary(m.4)


party.var <- (sum(m.1$residuals^2)- sum(m.2$residuals^2))/sum(m.1$residuals^2)
health.var <- (sum(m.1$residuals^2)- sum(m.3$residuals^2))/sum(m.1$residuals^2)

party.var/health.var



stargazer(m.1, m.2, m.3,m.4, type="text")
stargazer(m.1, m.2, m.3,m.4,
          star.char = c("*"),
          star.cutoffs = c(0.05),
          style="ajps",
          digits=2,
          covariate.labels = c( "Time","Time2",
                                "Democrat", "Republican",
                                "State Change in Covid Cases/1000",
                                "Female",
                                "Income 30k-74K","Income 75k-149k","Income 150k+",
                                "Age 25-34", "Age 34-44", "Age 45-54","Age 55-64","Age 65+",
                                "Race Other","Black","Asian","Hispanic",
                                "High School","Associates/Some College","College","Postgraduate",
                                "Unemployed","Rural","Urban"),
          omit= c(),
          omit.stat = c("rsq"),
          dep.var.labels   = "Last 24 Hour Activity (0-7)",
          model.numbers = T,
          add.lines = list(c("State F.E.", "Yes", "Yes","Yes","Yes")),
          notes.append = FALSE,
          title= "Determinants of Social Distancing Behavior",
          label="pooled",
          type="html",
          out="Figures/TableS3.html")


######################################################
################ All Activities#######################
######################################################
vars <- names(dat)[which(names(dat)=="last24hr_mc_go_to_work") : which(names(dat)=="last24hr_mc_medical_care") ]
vars
labels <- c("Go to Work","Groceries","Exercise","Take a Walk","Visit Family/Friend","Go to Resturant","Seek Medical Care")
class(dat$last24hr_mc_go_to_work)

dem.coef <- rep(NA, length(vars))
dem.lean.coef <- rep(NA, length(vars))
rep.coef <- rep(NA, length(vars))
rep.lean.coef <- rep(NA, length(vars))
case.coef <- rep(NA, length(vars))

dem.se <- rep(NA, length(vars))
dem.lean.se <- rep(NA, length(vars))
rep.se <- rep(NA, length(vars))
rep.lean.se <- rep(NA, length(vars))
case.se <- rep(NA, length(vars))

for(j in 1:7){
  dat$dv <- dat[,vars[j]]
  
  party.var <- NA
  health.med <- NA
  
  m <- felm(dv ~time.trend + time.trend2 + pid5 + case.change + 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density
            | state| 0 | state , data=dat,
            weights=dat$weight_state_weekly)
  summary(m)
  
  dem.coef[j]      <- coef(m)["pid5Democrat"]
  dem.lean.coef[j] <- coef(m)["pid5Lean Democrat"]
  rep.coef[j] <- coef(m)["pid5Republican"]
  rep.lean.coef[j] <- coef(m)["pid5Lean Republican"]
  case.coef[j] <- coef(m)["case.change"]
  
  dem.se[j]      <- sqrt(vcov(m)["pid5Democrat","pid5Democrat"])
  dem.lean.se[j] <- sqrt(vcov(m)["pid5Lean Democrat","pid5Lean Democrat"])
  rep.se[j] <- sqrt(vcov(m)["pid5Republican","pid5Republican"])
  rep.lean.se[j] <- sqrt(vcov(m)["pid5Lean Republican","pid5Lean Republican"])
  case.se[j] <- sqrt(vcov(m)["case.change","case.change"])
  
  
}
pos1 <- seq(1,61,10)
pos2 <- seq(2,62,10)
pos3 <- seq(3,63,10)
pos4 <- seq(4,64,10)
pos5 <- seq(5,65,10)

line.pos <- seq(7.5,57.5,10) 


pdf(file="Figures/S4.pdf", height=9, width=13)
plot(case.coef, pos1, xlim=c(-.15,.15), ylim=c(0,70), pch=16, col="springgreen4", type="n", axes=F,
     ylab="", xlab="                                      Coefficient Predicting P(Activity)")

abline(v=0, lty=2, col="gray40")

points(case.coef,pos1, pch=16, col="springgreen4")
segments(case.coef, pos1,case.coef+1.96*case.se, pos1,col="springgreen4")
segments(case.coef, pos1,case.coef-1.96*case.se, pos1,col="springgreen4")

points(rep.coef,pos2, pch=16, col="firebrick4")
segments(rep.coef, pos2,rep.coef+1.96*rep.se, pos2,col="firebrick4")
segments(rep.coef, pos2,rep.coef-1.96*rep.se, pos2,col="firebrick4")

points(rep.lean.coef,pos3, pch=16, col="firebrick1")
segments(rep.lean.coef, pos3,rep.lean.coef+1.96*rep.lean.se, pos3,col="firebrick1")
segments(rep.lean.coef, pos3,rep.lean.coef-1.96*rep.lean.se, pos3,col="firebrick1")

points(dem.lean.coef,pos4, pch=16, col="dodgerblue1")
segments(dem.lean.coef, pos4,dem.lean.coef+1.96*dem.lean.se, pos4,col="dodgerblue1")
segments(dem.lean.coef, pos4,dem.lean.coef-1.96*dem.lean.se, pos4,col="dodgerblue1")

points(dem.coef,pos5, pch=16, col="dodgerblue4")
segments(dem.coef, pos5,dem.coef+1.96*dem.se, pos5,col="dodgerblue4")
segments(dem.coef, pos5,dem.coef-1.96*dem.se, pos5,col="dodgerblue4")

text(-.13, pos5, labels)

legend("topright", c("Democrat", "Lean Democrat", "Lean Republican", "Republican", "Change in County Cases per 1000"), 
       pch=rep(16,5), col = c("dodgerblue4", "dodgerblue1", "firebrick1", "firebrick4", "springgreen4"), bg="white")


abline(h=line.pos, lty=3, col="gray50")
axis(1, at=seq(-.1,.15,.05))
dev.off()


######################################################
######  Difference in Activities over Time ############
######################################################



vars <- names(dat)[which(names(dat)=="last24hr_mc_go_to_work") : which(names(dat)=="last24hr_mc_medical_care") ]
labels <- c("Go to Work","Groceries","Exercise","Take a Walk","Visit Family/Friend","Go to Resturant","Seek Medical Care")
class(dat$last24hr_mc_go_to_work)
col2rgb("firebrick4")
firebrickfaded <- rgb(139, 26, 26, max = 255, alpha = 100, names = "firebrick4")

col2rgb("dodgerblue4")
dodgerbluefaded <- rgb(16, 78, 139, max = 255, alpha = 100, names = "dodgerblue4")

col2rgb("purple")
purplefaded <- rgb(160, 32, 240, max = 255, alpha = 100, names = "purple")

col2rgb("goldenrod3")
goldenrodfaded <- rgb(205, 155, 29, max = 255, alpha = 100, names = "goldenrod3")

col2rgb("springgreen4")
springgreenfaded <- rgb(0, 139, 69, max = 255, alpha = 100, names = "springgreen4")


#Double check the numbers are the same (can't figure out with names)
#png(file="Figures/SMCoefficientsTimeAllActivities.png", height=1728, width=1536)
pdf(file="Figures/S5.pdf", height=18, width=16)
par(mfrow=c(3,3))
for(j in 1:7){
  dat$dv <- dat[,vars[j]]
  
  
  dem.coef <- NA
  rep.coef <- NA
  case.coef <- NA
  dem.se <- NA
  rep.se <- NA
  case.se <- NA
  for (i in 1:length(dates)){
    m <-       felm(dv ~ dem + rep + case.change +
                      female + income_30_75 + income_75_150 + income_above150 +
                      age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                      other_race + black + asian + hispanic + 
                      hs + associates_somecollege + college + postgrad + 
                      unemployed + density | state | 0 | state, data=dat[dat$SMDate==dates[i],],
                    weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
    
    dem.coef[i] <- coef(m)["dem"]
    rep.coef[i] <- coef(m)["rep"]
    case.coef[i] <- coef(m)["case.change"]
    dem.se[i] <- sqrt(vcov(m)["dem","dem"])
    rep.se[i] <- sqrt(vcov(m)["rep","rep"])
    case.se[i] <- sqrt(vcov(m)["case.change","case.change"])
  }
  
  
  
  
  #Loess estimation
  dem.coef.fit <- predict(loess(dem.coef ~ as.numeric(dates),  span=.4, weights = 1/dem.se), se=T)
  rep.coef.fit <- predict(loess(rep.coef ~ as.numeric(dates),  span=.4, weights = 1/rep.se), se=T)
  case.coef.fit <- predict(loess(case.coef ~ as.numeric(dates),  span=.4, weights = 1/case.se), se=T)
  
  
  
  plot(dates, dem.coef, pch=16, xlab="", ylab="Daily Coefficient on 24hr Activity", col=dodgerbluefaded,ylim=c(-.2,.2),
       axes=F,main=paste(labels[j]))
  abline(h=0,lty=2,col="gray60")
  segments(dates, dem.coef,dates, dem.coef+1.96*dem.se, col=dodgerbluefaded)
  segments(dates, dem.coef,dates, dem.coef-1.96*dem.se, col=dodgerbluefaded)
  
  lines(dates,dem.coef.fit$fit, col="dodgerblue4", lwd=2)
  lines(dates,dem.coef.fit$fit - 1.96*dem.coef.fit$se.fit, col="dodgerblue4", lwd=1, lty=2)
  lines(dates,dem.coef.fit$fit + 1.96*dem.coef.fit$se.fit, col="dodgerblue4", lwd=1, lty=2)
  
  
  points(dates, rep.coef, pch=16, col=firebrickfaded)
  segments(dates, rep.coef,dates, rep.coef+1.96*rep.se, col=firebrickfaded)
  segments(dates, rep.coef,dates, rep.coef-1.96*rep.se, col=firebrickfaded)
  
  
  lines(dates,rep.coef.fit$fit, col="firebrick4", lwd=2)
  lines(dates,rep.coef.fit$fit - 1.96*rep.coef.fit$se.fit, col="firebrick4", lwd=1, lty=2)
  lines(dates,rep.coef.fit$fit + 1.96*rep.coef.fit$se.fit, col="firebrick4", lwd=1, lty=2)
  
  
  points(dates, case.coef, pch=16, col=springgreenfaded)
  segments(dates, case.coef,dates, case.coef+1.96*case.se, col=springgreenfaded)
  segments(dates, case.coef,dates, case.coef-1.96*case.se, col=springgreenfaded)
  
  lines(dates,case.coef.fit$fit, col="springgreen4", lwd=2)
  lines(dates,case.coef.fit$fit - 1.96*case.coef.fit$se.fit, col="springgreen4", lwd=1, lty=2)
  lines(dates,case.coef.fit$fit + 1.96*case.coef.fit$se.fit, col="springgreen4", lwd=1, lty=2)
  
  axis(2,at=seq(-.2,.2,.05),las=2)
  axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
  
  legend("topleft", c("Democrat", "Republican" , "Change in County cases per 1000"), pch=c(16,16,16),
         col=c("dodgerblue4", "firebrick4", "springgreen4"), bg="white")  
  
  
}
dev.off()



######################################################
#####  Alternative Public Health Measures ############
######################################################

col2rgb("goldenrod3")
goldenrodfaded <- rgb(205, 155, 29, max = 255, alpha = 100, names = "goldenrod3")

col2rgb("springgreen4")
springgreenfaded <- rgb(0, 139, 69, max = 255, alpha = 100, names = "springgreen4")
#Rescaling all to be per 1000 (originally had per 10000)
summary(dat$case.change)

summary(dat$death.change)

#### Change in county deaths

party.med <- NA
health.med <- NA


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
  
  
  
  m.h <- felm(last24hr_index ~ death.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  

}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/SMCountyDeathChange.png", height=576, width=768)
pdf("Figures/S6.pdf", height=6, width=8)

plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "Change in County Deaths/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()

#### County Cases
summary(dat$cases)

party.med <- NA
health.med <- NA

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
  
  
  
  m.h <- felm(last24hr_index ~ cases +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/SMCountyCases.png", height=576, width=768)
pdf("Figures/S7.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "County Cases/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()

#### County Deaths

summary(dat$deaths)

party.med <- NA
health.med <- NA

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
  
  
  
  m.h <- felm(last24hr_index ~ deaths +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| state | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/SMCountyDeaths.png", height=576, width=768)
pdf("Figures/S8.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "County Deaths/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()


#### State Case Change

summary(dat$state.case.change)


party.med <- NA
health.med <- NA

for (i in 1:length(dates)){
  m <- felm(last24hr_index ~ 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density|0 | 0 | state, data=dat[dat$SMDate==dates[i],],
            weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.p <- felm(last24hr_index ~  dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density|0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.h <- felm(last24hr_index ~ state.case.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density|0| 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/S9.png", height=576, width=768)
pdf("Figures/S9.pdf", height=6, width=8)

plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "Change in State Cases/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()

#### State Death Change


summary(dat$state.death.change)


party.med <- NA
health.med <- NA

for (i in 1:length(dates)){
  m <- felm(last24hr_index ~ 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density|0 | 0 | state, data=dat[dat$SMDate==dates[i],],
            weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.p <- felm(last24hr_index ~  dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density|0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.h <- felm(last24hr_index ~ state.death.change +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density|0| 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/SMChangeStateDeaths.png", height=576, width=768)
pdf("Figures/S10.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "Change in State Deaths/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()

#State Deaths/1000


party.med <- NA
health.med <- NA

for (i in 1:length(dates)){
  m <- felm(last24hr_index ~ 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density|0 | 0 | state, data=dat[dat$SMDate==dates[i],],
            weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.p <- felm(last24hr_index ~  dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| 0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.h <- felm(last24hr_index ~ state.deaths +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density|0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/SMCountyDeaths.png", height=576, width=768)
pdf("Figures/S11.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "State Deaths/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()



#State Cases/1000

summary(dat$state.cases)

party.med <- NA
health.med <- NA

for (i in 1:length(dates)){
  m <- felm(last24hr_index ~ 
              female + income_30_75 + income_75_150 + income_above150 +
              age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
              other_race + black + asian + hispanic + 
              hs + associates_somecollege + college + postgrad + 
              unemployed + density| 0 | 0 | state, data=dat[dat$SMDate==dates[i],],
            weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.p <- felm(last24hr_index ~  dem + rep +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| 0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  m.h <- felm(last24hr_index ~ state.cases +
                female + income_30_75 + income_75_150 + income_above150 +
                age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                other_race + black + asian + hispanic + 
                hs + associates_somecollege + college + postgrad + 
                unemployed + density| 0 | 0 | state , data=dat[dat$SMDate==dates[i],],
              weights=dat$weight_national_18plus_daily[dat$SMDate==dates[i]],)
  
  
  
  
  party.med[i] <- (sum(m$residuals^2)- sum(m.p$residuals^2))/sum(m$residuals^2)
  health.med[i] <- (sum(m$residuals^2)- sum(m.h$residuals^2))/sum(m$residuals^2)
  
  
}

#Loess estimation
party.med.fit <- predict(loess(party.med*100 ~ as.numeric(dates)), se=T)
health.med.fit <- predict(loess(health.med*100 ~ as.numeric(dates)), se=T)


#png("Figures/S12.png", height=576, width=768)
pdf("Figures/S12.pdf", height=6, width=8)
plot(dates, party.med*100, pch=16, col=goldenrodfaded,  ylab="Percent of Remaining Variance Explained",xlab="Date" , ylim=c(0,10), axes=F)

lines(dates,party.med.fit$fit, col="goldenrod3", lwd=2)
lines(dates,party.med.fit$fit - 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)
lines(dates,party.med.fit$fit + 1.96*party.med.fit$se.fit, col="goldenrod3", lwd=1, lty=2)

points(dates, health.med*100, pch=16, col=springgreenfaded)

lines(dates,health.med.fit$fit, col="springgreen4", lwd=2)
lines(dates,health.med.fit$fit - 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)
lines(dates,health.med.fit$fit + 1.96*health.med.fit$se.fit, col="springgreen4", lwd=1, lty=2)

axis.Date(1, at=seq(as.Date("2020-03-20"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,10,2), las=2)
legend("topleft", c("Party", "State Cases/1000"), pch=c(16,16),
       col=c("goldenrod3","springgreen4"), bg="white")
dev.off()


######################################################
####################  Age Trend ######################
######################################################


dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

age1.avg.24 <- NA
age2.avg.24 <- NA
age3.avg.24 <- NA
age4.avg.24 <- NA
age5.avg.24 <- NA
age6.avg.24 <- NA


age1.se.24 <- NA
age2.se.24 <- NA
age3.se.24 <- NA
age4.se.24 <- NA
age5.se.24 <- NA
age6.se.24 <- NA


for(i in 1:length(dates)){
  
  age1.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_18_24==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_18_24==1],na.rm=T)
  age2.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_25_34==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_25_34==1],na.rm=T)
  age3.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_35_44==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_35_44==1],na.rm=T)
  age4.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_45_54==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_45_54==1],na.rm=T)
  age5.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_55_64==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_55_64==1],na.rm=T)
  age6.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$age_65_p==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_65_p==1],na.rm=T)
  
  age1.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_18_24==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_18_24==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_18_24==1]))
  age2.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_25_34==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_25_34==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_25_34==1]))
  age3.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_35_44==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_35_44==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_35_44==1]))
  age4.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_45_54==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_45_54==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_45_54==1]))
  age5.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_55_64==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_55_64==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_55_64==1]))
  age6.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_65_p==1],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_65_p==1],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$age_65_p==1]))
}


age1.fit.24 <- predict(loess(age1.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age1.se.24), se=T)
age2.fit.24 <- predict(loess(age2.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age2.se.24), se=T)
age3.fit.24 <- predict(loess(age3.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age3.se.24), se=T)
age4.fit.24 <- predict(loess(age4.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age4.se.24), se=T)
age5.fit.24 <- predict(loess(age5.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age5.se.24), se=T)
age6.fit.24 <- predict(loess(age6.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/age6.se.24), se=T)



dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))




age1.avg.wor <- NA
age2.avg.wor <- NA
age3.avg.wor <- NA
age4.avg.wor <- NA
age5.avg.wor <- NA
age6.avg.wor <- NA


age1.se.wor <- NA
age2.se.wor <- NA
age3.se.wor <- NA
age4.se.wor <- NA
age5.se.wor <- NA
age6.se.wor <- NA



for(i in 1:length(dates)){
  age1.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_18_24==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_18_24==1],na.rm=T)
  age2.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_25_34==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_25_34==1],na.rm=T)
  age3.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_35_44==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_35_44==1],na.rm=T)
  age4.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_45_54==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_45_54==1],na.rm=T)
  age5.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_55_64==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_55_64==1],na.rm=T)
  age6.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$age_65_p==1],
                                  w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$age_65_p==1],na.rm=T)
  
  age1.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_18_24==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_18_24==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_18_24==1]))
  age2.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_25_34==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_25_34==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_25_34==1]))
  age3.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_35_44==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_35_44==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_35_44==1]))
  age4.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_45_54==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_45_54==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_45_54==1]))
  age5.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_55_64==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_55_64==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_55_64==1]))
  age6.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$age_65_p==1],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$age_65_p==1],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$age_65_p==1]))
  
}




age1.avg.wor[age1.avg.wor==0]<- NA
age1.se.wor[age1.avg.wor==0] <- NA

age1.fit.wor <- predict(loess(age1.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age1.se.wor), se=T)
age2.fit.wor <- predict(loess(age2.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age2.se.wor), se=T)
age3.fit.wor <- predict(loess(age3.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age3.se.wor), se=T)
age4.fit.wor <- predict(loess(age4.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age4.se.wor), se=T)
age5.fit.wor <- predict(loess(age5.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age5.se.wor), se=T)
age6.fit.wor <- predict(loess(age6.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/age6.se.wor), se=T)


#Colors
col2rgb("firebrick4")
firebrickfaded <- rgb(139, 26, 26, max = 255, alpha = 100, names = "firebrick4")

col2rgb("dodgerblue4")
dodgerbluefaded <- rgb(16, 78, 139, max = 255, alpha = 100, names = "dodgerblue4")

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


magma1 <- rgb(253,220,158, max = 255)
magma1faded <- rgb(253,220,158, max = 255, alpha = 100)

magma2 <- rgb(253,153,105, max = 255)
magma2faded <- rgb(253,153,105, max = 255, alpha = 100)

magma3 <- rgb(235,87,96,max = 255)
magma3faded <- rgb(235,87,96, max = 255, alpha = 100)

magma4 <- rgb(182,54,121, max = 255)
magma4faded <- rgb(182,54,121, max = 255, alpha = 100)

magma5 <- rgb(123,35,130, max = 255)
magma5faded <- rgb(123,35,130, max = 255, alpha = 100)

magma6 <- rgb(35,17,81, max = 255)
magma6faded <- rgb(35,17,81, max = 255, alpha = 100)


pdf(file="Figures/S13.pdf", height=12, width=8)
#png(file="Figures/SMAgeTrend.png", width=1152, height=768)
par(mfrow=c(2,1))
dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))

plot(dates, age1.avg.wor,type="p", lwd=2, col=magma1faded, ylim=c(0,1),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Proportion Very or Somewhat Worried about COVID-19",
     xlab="Date", pch=16, axes=F,main="A",cex=.5)
segments(dates,age1.avg.wor, dates,age1.avg.wor+age1.se.wor*1.96, col=magma1faded)
segments(dates,age1.avg.wor, dates,age1.avg.wor-age1.se.wor*1.96, col=magma1faded)


points(dates, age2.avg.wor,type="p",col=magma2faded , pch=16, cex=.5)
segments(dates,age2.avg.wor, dates,age2.avg.wor+age2.se.wor*1.96, col=magma2faded)
segments(dates,age2.avg.wor, dates,age2.avg.wor-age2.se.wor*1.96, col=magma2faded)

points(dates, age3.avg.wor,type="p",col=magma3faded , pch=16, cex=.5)
segments(dates,age3.avg.wor, dates,age3.avg.wor+age3.se.wor*1.96, col=magma3faded)
segments(dates,age3.avg.wor, dates,age3.avg.wor-age3.se.wor*1.96, col=magma3faded)

points(dates, age4.avg.wor,type="p",col=magma4faded , pch=16, cex=.5)
segments(dates,age4.avg.wor, dates,age4.avg.wor+age4.se.wor*1.96, col=magma4faded)
segments(dates,age4.avg.wor, dates,age4.avg.wor-age4.se.wor*1.96, col=magma4faded)

points(dates, age5.avg.wor,type="p",col=magma5faded , pch=16, cex=.5)
segments(dates,age5.avg.wor, dates,age5.avg.wor+age5.se.wor*1.96, col=magma5faded)
segments(dates,age5.avg.wor, dates,age5.avg.wor-age5.se.wor*1.96, col=magma5faded)


points(dates, age6.avg.wor,type="p",col=magma6faded , pch=16, cex=.5)
segments(dates,age6.avg.wor, dates,age6.avg.wor+age6.se.wor*1.96, col=magma6faded)
segments(dates,age6.avg.wor, dates,age6.avg.wor-age6.se.wor*1.96, col=magma6faded)


lines(dates[!is.na(age1.avg.wor)],age1.fit.wor$fit,col=magma1, lwd=2)

lines(dates[!is.na(age1.avg.wor)],age1.fit.wor$fit - 1.96*age1.fit.wor$se.fit, col=magma1, lwd=1, lty=2)
lines(dates[!is.na(age1.avg.wor)],age1.fit.wor$fit + 1.96*age1.fit.wor$se.fit, col=magma1, lwd=1, lty=2)

lines(dates,age2.fit.wor$fit,col=magma2, lwd=2)
lines(dates,age2.fit.wor$fit - 1.96*age2.fit.wor$se.fit, col=magma2, lwd=1, lty=2)
lines(dates,age2.fit.wor$fit + 1.96*age2.fit.wor$se.fit, col=magma2, lwd=1, lty=2)

lines(dates,age3.fit.wor$fit,col=magma3, lwd=2)
lines(dates,age3.fit.wor$fit - 1.96*age3.fit.wor$se.fit, col=magma3, lwd=1, lty=2)
lines(dates,age3.fit.wor$fit + 1.96*age3.fit.wor$se.fit, col=magma3, lwd=1, lty=2)

lines(dates,age4.fit.wor$fit,col=magma4, lwd=2)
lines(dates,age4.fit.wor$fit - 1.96*age4.fit.wor$se.fit, col=magma4, lwd=1, lty=2)
lines(dates,age4.fit.wor$fit + 1.96*age4.fit.wor$se.fit, col=magma4, lwd=1, lty=2)

lines(dates,age5.fit.wor$fit,col=magma5, lwd=2)
lines(dates,age5.fit.wor$fit - 1.96*age5.fit.wor$se.fit, col=magma5, lwd=1, lty=2)
lines(dates,age5.fit.wor$fit + 1.96*age5.fit.wor$se.fit, col=magma5, lwd=1, lty=2)


lines(dates,age6.fit.wor$fit,col=magma6, lwd=2)
lines(dates,age6.fit.wor$fit - 1.96*age6.fit.wor$se.fit, col=magma6, lwd=1, lty=2)
lines(dates,age6.fit.wor$fit + 1.96*age6.fit.wor$se.fit, col=magma6, lwd=1, lty=2)


legend("topleft", c("Age 18-24","Age 25-34","Age 35-44","Age 45-54","Age 55-64","Age 65+"),
       pch=rep(16,4), col=c(magma1,magma2,magma3,magma4,magma5,magma6), bg="white")

axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,1,0.25))


dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))

plot(dates, age1.avg.24,type="p", lwd=2, col=magma1faded, ylim=c(1,3.5),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Number of Activities in Previous 24hr",
     xlab="Date", pch=16, axes=F,main="B",cex=.5)
segments(dates,age1.avg.24, dates,age1.avg.24+age1.se.24*1.96, col=magma1faded)
segments(dates,age1.avg.24, dates,age1.avg.24-age1.se.24*1.96, col=magma1faded)


points(dates, age2.avg.24,type="p",col=magma2faded , pch=16, cex=.5)
segments(dates,age2.avg.24, dates,age2.avg.24+age2.se.24*1.96, col=magma2faded)
segments(dates,age2.avg.24, dates,age2.avg.24-age2.se.24*1.96, col=magma2faded)

points(dates, age3.avg.24,type="p",col=magma3faded , pch=16, cex=.5)
segments(dates,age3.avg.24, dates,age3.avg.24+age3.se.24*1.96, col=magma3faded)
segments(dates,age3.avg.24, dates,age3.avg.24-age3.se.24*1.96, col=magma3faded)

points(dates, age4.avg.24,type="p",col=magma4faded , pch=16, cex=.5)
segments(dates,age4.avg.24, dates,age4.avg.24+age4.se.24*1.96, col=magma4faded)
segments(dates,age4.avg.24, dates,age4.avg.24-age4.se.24*1.96, col=magma4faded)

points(dates, age5.avg.24,type="p",col=magma5faded , pch=16, cex=.5)
segments(dates,age5.avg.24, dates,age5.avg.24+age5.se.24*1.96, col=magma5faded)
segments(dates,age5.avg.24, dates,age5.avg.24-age5.se.24*1.96, col=magma5faded)


points(dates, age6.avg.24,type="p",col=magma6faded , pch=16, cex=.5)
segments(dates,age6.avg.24, dates,age6.avg.24+age6.se.24*1.96, col=magma6faded)
segments(dates,age6.avg.24, dates,age6.avg.24-age6.se.24*1.96, col=magma6faded)


lines(dates,age1.fit.24$fit,col=magma1, lwd=2)
lines(dates,age1.fit.24$fit - 1.96*age1.fit.24$se.fit, col=magma1, lwd=1, lty=2)
lines(dates,age1.fit.24$fit + 1.96*age1.fit.24$se.fit, col=magma1, lwd=1, lty=2)

lines(dates,age2.fit.24$fit,col=magma2, lwd=2)
lines(dates,age2.fit.24$fit - 1.96*age2.fit.24$se.fit, col=magma2, lwd=1, lty=2)
lines(dates,age2.fit.24$fit + 1.96*age2.fit.24$se.fit, col=magma2, lwd=1, lty=2)

lines(dates,age3.fit.24$fit,col=magma3, lwd=2)
lines(dates,age3.fit.24$fit - 1.96*age3.fit.24$se.fit, col=magma3, lwd=1, lty=2)
lines(dates,age3.fit.24$fit + 1.96*age3.fit.24$se.fit, col=magma3, lwd=1, lty=2)

lines(dates,age4.fit.24$fit,col=magma4, lwd=2)
lines(dates,age4.fit.24$fit - 1.96*age4.fit.24$se.fit, col=magma4, lwd=1, lty=2)
lines(dates,age4.fit.24$fit + 1.96*age4.fit.24$se.fit, col=magma4, lwd=1, lty=2)

lines(dates,age5.fit.24$fit,col=magma5, lwd=2)
lines(dates,age5.fit.24$fit - 1.96*age5.fit.24$se.fit, col=magma5, lwd=1, lty=2)
lines(dates,age5.fit.24$fit + 1.96*age5.fit.24$se.fit, col=magma5, lwd=1, lty=2)


lines(dates,age6.fit.24$fit,col=magma6, lwd=2)
lines(dates,age6.fit.24$fit - 1.96*age6.fit.24$se.fit, col=magma6, lwd=1, lty=2)
lines(dates,age6.fit.24$fit + 1.96*age6.fit.24$se.fit, col=magma6, lwd=1, lty=2)


legend("topleft", c("Age 18-24","Age 25-34","Age 35-44","Age 45-54","Age 55-64","Age 65+"),
       pch=rep(16,4), col=c(magma1,magma2,magma3,magma4,magma5,magma6), bg="white")

axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(1,3.5,.5), las=2)

dev.off()
par(mfrow=c(1,1))


######################################################
##########  Census Region Trend ######################
######################################################

states <- read.csv("State Master.csv")

keep <- c("state.full","census.region")
states <- states[keep]

dat <- merge(dat,states,by.x="state", by.y="state.full",all.x=T)
table(dat$census.region)


dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

mw.avg.24 <- NA
ne.avg.24 <- NA
so.avg.24 <- NA
we.avg.24 <- NA


mw.se.24 <- NA
ne.se.24 <- NA
so.se.24 <- NA
we.se.24 <- NA



for(i in 1:length(dates)){
  
  mw.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$census.region=="midwest"],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="midwest"],na.rm=T)
  ne.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$census.region=="northeast"],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="northeast"],na.rm=T)
  so.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$census.region=="south"],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="south"],na.rm=T)
  we.avg.24[i] <-weighted.mean(dat$last24hr_index[dat$SMDate==dates[i] &  dat$census.region=="west"],
                                 w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="west"],na.rm=T)
 
  mw.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="midwest"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="midwest"],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="midwest"]))
  ne.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="northeast"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="northeast"],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="northeast"]))
  so.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="south"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="south"],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="south"]))
  we.se.24[i] <- sqrt(wtd.var(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="west"],
                                w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="west"],na.rm=T))/
    sqrt(length(dat$last24hr_index[dat$SMDate==dates[i] & dat$census.region=="west"]))
  
}


mw.fit.24 <- predict(loess(mw.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/mw.se.24), se=T)
ne.fit.24 <- predict(loess(ne.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/ne.se.24), se=T)
so.fit.24 <- predict(loess(so.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/so.se.24), se=T)
we.fit.24 <- predict(loess(we.avg.24 ~ as.numeric(dates),  span=.4, weights = 1/we.se.24), se=T)


dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))




mw.avg.wor <- NA
ne.avg.wor <- NA
so.avg.wor <- NA
we.avg.wor <- NA


mw.se.wor <- NA
ne.se.wor <- NA
so.se.wor <- NA
we.se.wor <- NA



for(i in 1:length(dates)){
  
  mw.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$census.region=="midwest"],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="midwest"],na.rm=T)
  ne.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$census.region=="northeast"],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="northeast"],na.rm=T)
  so.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$census.region=="south"],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="south"],na.rm=T)
  we.avg.wor[i] <-weighted.mean(dat$is.worried[dat$SMDate==dates[i] &  dat$census.region=="west"],
                               w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] &  dat$census.region=="west"],na.rm=T)
  
  mw.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="midwest"],
                              w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="midwest"],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="midwest"]))
  ne.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="northeast"],
                              w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="northeast"],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="northeast"]))
  so.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="south"],
                              w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="south"],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="south"]))
  we.se.wor[i] <- sqrt(wtd.var(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="west"],
                              w=dat$weight_national_18plus_daily[dat$SMDate==dates[i] & dat$census.region=="west"],na.rm=T))/
    sqrt(length(dat$is.worried[dat$SMDate==dates[i] & dat$census.region=="west"]))
  
}






mw.fit.wor <- predict(loess(mw.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/mw.se.wor), se=T)
ne.fit.wor <- predict(loess(ne.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/ne.se.wor), se=T)
so.fit.wor <- predict(loess(so.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/so.se.wor), se=T)
we.fit.wor <- predict(loess(we.avg.wor ~ as.numeric(dates),  span=.4, weights = 1/we.se.wor), se=T)


#Colors
col2rgb("firebrick4")
firebrickfaded <- rgb(139, 26, 26, max = 255, alpha = 100, names = "firebrick4")

col2rgb("dodgerblue4")
dodgerbluefaded <- rgb(16, 78, 139, max = 255, alpha = 100, names = "dodgerblue4")

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




pdf(file="Figures/S14.pdf", height=12, width=8)
#png(file="Figures/SMRegionTrend.png", width=1152, height=768)
par(mfrow=c(2,1))
dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))

plot(dates, mw.avg.wor,type="p", lwd=2, col=viridis1faded, ylim=c(0,1),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Proportion Very or Somewhat Worried about COVID-19",
     xlab="Date", pch=16, axes=F,main="A",cex=.5)
segments(dates,mw.avg.wor, dates,mw.avg.wor+mw.se.wor*1.96, col=viridis1faded)
segments(dates,mw.avg.wor, dates,mw.avg.wor-mw.se.wor*1.96, col=viridis1faded)


points(dates, ne.avg.wor,type="p",col=viridis2faded , pch=16, cex=.5)
segments(dates,ne.avg.wor, dates,ne.avg.wor+ne.se.wor*1.96, col=viridis2faded)
segments(dates,ne.avg.wor, dates,ne.avg.wor-ne.se.wor*1.96, col=viridis2faded)

points(dates, so.avg.wor,type="p",col=viridis3faded , pch=16, cex=.5)
segments(dates,so.avg.wor, dates,so.avg.wor+so.se.wor*1.96, col=viridis3faded)
segments(dates,so.avg.wor, dates,so.avg.wor-so.se.wor*1.96, col=viridis3faded)

points(dates, we.avg.wor,type="p",col=viridis4faded , pch=16, cex=.5)
segments(dates,we.avg.wor, dates,we.avg.wor+we.se.wor*1.96, col=viridis4faded)
segments(dates,we.avg.wor, dates,we.avg.wor-we.se.wor*1.96, col=viridis4faded)



lines(dates,mw.fit.wor$fit,col=viridis1, lwd=2)
lines(dates,mw.fit.wor$fit - 1.96*mw.fit.wor$se.fit, col=viridis1, lwd=1, lty=2)
lines(dates,mw.fit.wor$fit + 1.96*mw.fit.wor$se.fit, col=viridis1, lwd=1, lty=2)

lines(dates,ne.fit.wor$fit,col=viridis2, lwd=2)
lines(dates,ne.fit.wor$fit - 1.96*ne.fit.wor$se.fit, col=viridis2, lwd=1, lty=2)
lines(dates,ne.fit.wor$fit + 1.96*ne.fit.wor$se.fit, col=viridis2, lwd=1, lty=2)

lines(dates,so.fit.wor$fit,col=viridis3, lwd=2)
lines(dates,so.fit.wor$fit - 1.96*so.fit.wor$se.fit, col=viridis3, lwd=1, lty=2)
lines(dates,so.fit.wor$fit + 1.96*so.fit.wor$se.fit, col=viridis3, lwd=1, lty=2)

lines(dates,we.fit.wor$fit,col=viridis4, lwd=2)
lines(dates,we.fit.wor$fit - 1.96*we.fit.wor$se.fit, col=viridis4, lwd=1, lty=2)
lines(dates,we.fit.wor$fit + 1.96*we.fit.wor$se.fit, col=viridis4, lwd=1, lty=2)



legend("topleft", c("Mid-West","Northeast","South","West"),
       pch=rep(16,4), col=c(viridis1,viridis2,viridis3,viridis4), bg="white")

axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(0,1,0.25))


dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))

plot(dates, mw.avg.24,type="p", lwd=2, col=viridis1faded, ylim=c(1,3.5),
     xlim=as.Date(c("2020-02-11","2020-10-01")), ylab="Number of Activities in Previous 24hr",
     xlab="Date", pch=16, axes=F,main="B",cex=.5)
segments(dates,mw.avg.24, dates,mw.avg.24+mw.se.24*1.96, col=viridis1faded)
segments(dates,mw.avg.24, dates,mw.avg.24-mw.se.24*1.96, col=viridis1faded)


points(dates, ne.avg.24,type="p",col=viridis2faded , pch=16, cex=.5)
segments(dates,ne.avg.24, dates,ne.avg.24+ne.se.24*1.96, col=viridis2faded)
segments(dates,ne.avg.24, dates,ne.avg.24-ne.se.24*1.96, col=viridis2faded)

points(dates, so.avg.24,type="p",col=viridis3faded , pch=16, cex=.5)
segments(dates,so.avg.24, dates,so.avg.24+so.se.24*1.96, col=viridis3faded)
segments(dates,so.avg.24, dates,so.avg.24-so.se.24*1.96, col=viridis3faded)

points(dates, we.avg.24,type="p",col=viridis4faded , pch=16, cex=.5)
segments(dates,we.avg.24, dates,we.avg.24+we.se.24*1.96, col=viridis4faded)
segments(dates,we.avg.24, dates,we.avg.24-we.se.24*1.96, col=viridis4faded)



lines(dates,mw.fit.24$fit,col=viridis1, lwd=2)
lines(dates,mw.fit.24$fit - 1.96*mw.fit.24$se.fit, col=viridis1, lwd=1, lty=2)
lines(dates,mw.fit.24$fit + 1.96*mw.fit.24$se.fit, col=viridis1, lwd=1, lty=2)

lines(dates,ne.fit.24$fit,col=viridis2, lwd=2)
lines(dates,ne.fit.24$fit - 1.96*ne.fit.24$se.fit, col=viridis2, lwd=1, lty=2)
lines(dates,ne.fit.24$fit + 1.96*ne.fit.24$se.fit, col=viridis2, lwd=1, lty=2)

lines(dates,so.fit.24$fit,col=viridis3, lwd=2)
lines(dates,so.fit.24$fit - 1.96*so.fit.24$se.fit, col=viridis3, lwd=1, lty=2)
lines(dates,so.fit.24$fit + 1.96*so.fit.24$se.fit, col=viridis3, lwd=1, lty=2)

lines(dates,we.fit.24$fit,col=viridis4, lwd=2)
lines(dates,we.fit.24$fit - 1.96*we.fit.24$se.fit, col=viridis4, lwd=1, lty=2)
lines(dates,we.fit.24$fit + 1.96*we.fit.24$se.fit, col=viridis4, lwd=1, lty=2)


legend("topleft", c("Mid-West","Northeast","South","West"),
       pch=rep(16,4), col=c(viridis1,viridis2,viridis3,viridis4), bg="white")

axis.Date(1, at=seq(as.Date("2020-02-11"),as.Date("2020-10-01") , by="1 week"), format="%m-%d")
axis(2, at=seq(1,3.5,.5), las=2)

dev.off()
par(mfrow=c(1,1))


######################################################
#########  Partisan Trends for Every State ###########
######################################################
dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

states <- unique(dat$state)

#No DC
states <- states[-9]



col2rgb("dodgerblue4")
dodgerbluefaded <- rgb(16, 78, 139, max = 255, alpha = 100, names = "dodgerblu4")



col2rgb("firebrick4")
firebrickfaded <- rgb(139, 26, 26, max = 255, alpha = 100, names = "firebrick4")


pdf(file="Figures/S15.pdf", height=20, width=18)
par(mfrow=c(10,5))

for(j in 1: length(states)){
  dat.state <- dat[dat$state==states[j],]
  
  
  dates <- sort(unique(dat.state$SMDate[!is.na(dat.state$last24hr_index)]))
  
  dem.avg.24 <- NA
  rep.avg.24 <- NA
  
  dem.se.24 <- NA
  rep.se.24 <- NA
  
  
  
  for(i in 1:length(dates)){
    dem.avg.24[i] <-weighted.mean(dat.state$last24hr_index[dat.state$SMDate==dates[i] & (dat.state$pid5=="Democrat" | dat.state$pid5=="Lean Democrat")],
                                  w=dat.state$weight_national_18plus_daily[dat.state$SMDate==dates[i] & (dat.state$pid5=="Democrat" | dat.state$pid5=="Lean Democrat")],na.rm=T)
    rep.avg.24[i] <-weighted.mean(dat.state$last24hr_index[dat.state$SMDate==dates[i] & (dat.state$pid5=="Republican" | dat.state$pid5=="Lean Republican")],
                                  w=dat.state$weight_national_18plus_daily[dat.state$SMDate==dates[i] & (dat.state$pid5=="Republican" | dat.state$pid5=="Lean Republican")],na.rm=T)
    
  }
  
  
  
  
  
  dates.dem <- dates[!is.nan(dem.avg.24)]
  dates.rep <- dates[!is.nan(rep.avg.24)]
  
  #Loess estimation last24
  dem.fit.24 <- predict(loess(dem.avg.24 ~ as.numeric(dates), span=.4), se=T)
  rep.fit.24 <- predict(loess(rep.avg.24 ~ as.numeric(dates), span=.4), se=T)
  
  
  
  
  plot(dates.dem, dem.fit.24$fit, type="n", lwd=2, col=dodgerbluefaded, ylim=c(1,4),xlim=as.Date(c("2020-04-01","2020-10-01")),
       ylab="Number of Activities in Previous 24hr", xlab="Date", pch=16, axes=F,main=print(states[j]),cex=.5)
  
  lines(dates.dem,dem.fit.24$fit,col="dodgerblue4", lwd=2)
  lines(dates.dem,dem.fit.24$fit - 1.96*dem.fit.24$se.fit, col="dodgerblue4", lwd=1, lty=2)
  lines(dates.dem,dem.fit.24$fit + 1.96*dem.fit.24$se.fit, col="dodgerblue4", lwd=1, lty=2)
  
  
  lines(dates.rep,rep.fit.24$fit,col="firebrick4", lwd=2)
  lines(dates.rep,rep.fit.24$fit - 1.96*rep.fit.24$se.fit, col="firebrick4", lwd=1, lty=2)
  lines(dates.rep,rep.fit.24$fit + 1.96*rep.fit.24$se.fit, col="firebrick4", lwd=1, lty=2)
  
  axis.Date(1, at=seq(as.Date("2020-04-01"),as.Date("2020-10-01") , by="1 month"), format="%m-%d")
  
  axis(2, at=seq(1,4,1), las=2)
}
dev.off()




######################################################
#############  Moderation by Media ###################
######################################################

summary(dat$high.news)
summary(dat$right.wing.media)
summary(dat$left.wing.media)

#MT Margins function
Rmargins <- function(model,xname,xzname,vcov.matrix,levels, ci=.95, latex=F){
  # Generate Marginal Effects
  m_effect <- coef(model)[xname]+coef(model)[xzname]*levels
  #Variances and Covariances
  v_x <- vcov.matrix[xname,xname]
  v_xz <- vcov.matrix[xzname,xzname]
  cov <- vcov.matrix[xname,xzname]
  #Standard Errors
  se <- sqrt(v_x + (levels^2)*v_xz+2*levels*cov)
  #T value for 95%CI
  if (ci==.95){
    t <- qt(0.025,model$df)
  }
  #T value for 90%CI
  if (ci==.9){
    t <- qt(0.05,model$df)
  }
  #Confidence Bounds
  m_upper <- m_effect + se*t
  m_lower <- m_effect - se*t
  # Remove Flotsom and Jetson
  #Printing Table
  table <- cbind(levels,m_effect,se,m_upper,m_lower)
  if (latex==T){
    library(xtable)
    print(xtable(table),include.rownames=FALSE)
  }
  if (latex==F){
    return(table)
  }
}


#Are news watchers more affected?

dates <- sort(unique(dat$SMDate[!is.na(dat$last24hr_index)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

m.news <- felm(last24hr_index ~time.trend + time.trend2+ dem*high.news + rep*high.news + case.change*high.news + 
                 female + income_30_75 + income_75_150 + income_above150 +
                 age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                 other_race + black + asian + hispanic + 
                 hs + associates_somecollege + college + postgrad
               | state | 0 |  state, data=dat,
               weights=dat$weight_state_weekly)
summary(m.news)


table(dat$high.news, dat$SMDate)
table(is.na(dat$high.news[dat$SMDate<"2020-06-10"]))


low.news.coef <- c(NA,NA,NA)
low.news.se <- c(NA,NA,NA)

high.news.coef <- c(NA, NA, NA)
high.news.se <- c(NA,NA,NA)


low.news.coef[1] <- coef(m.news)["dem"]
low.news.coef[2] <- coef(m.news)["rep"]
low.news.coef[3] <- coef(m.news)["case.change"]

low.news.se[1] <- sqrt(vcov(m.news)["dem","dem"])
low.news.se[2] <- sqrt(vcov(m.news)["rep","rep"])
low.news.se[3] <- sqrt(vcov(m.news)["case.change","case.change"])


high.news.coef[1] <- coef(m.news)["dem"] + coef(m.news)["dem:high.news"]
high.news.coef[2] <- coef(m.news)["rep"] + coef(m.news)["high.news:rep"]
high.news.coef[3] <- coef(m.news)["case.change"] + coef(m.news)["high.news:case.change"]


#sqrt(v_x + (levels^2)*v_xz+2*levels*cov)
high.news.se[1] <- sqrt(vcov(m.news)["dem","dem"] + vcov(m.news)["dem:high.news","dem:high.news"] + 2* vcov(m.news)["dem","dem:high.news"])
high.news.se[2] <- sqrt(vcov(m.news)["rep","rep"] + vcov(m.news)["high.news:rep","high.news:rep"] + 2* vcov(m.news)["rep","high.news:rep"])
high.news.se[3] <- sqrt(vcov(m.news)["case.change","case.change"] + vcov(m.news)["high.news:case.change","high.news:case.change"] + 2* vcov(m.news)["case.change","high.news:case.change"])


pdf(file="Figures/S16.pdf", height=6, width=8)
plot(1:3, low.news.coef, xlim=c(0.5,3.5), ylim=c(-0.7,0.7), col=c("dodgerblue","firebrick","springgreen4"), pch=rep(16,3), cex=1.5, main="",
     xlab="", ylab="Coefficient Predicting Last 24hr Activity", axes=F)

segments(1:3, low.news.coef, 1:3, low.news.coef+low.news.se*1.96, col=c("dodgerblue","firebrick","springgreen4"), lwd=3)
segments(1:3, low.news.coef, 1:3, low.news.coef-low.news.se*1.96, col=c("dodgerblue","firebrick","springgreen4"), lwd=3)


points(1.1:3.1, high.news.coef,col=c("dodgerblue","firebrick","springgreen4"), pch=rep(17,3), cex=1.5)

segments(1.1:3.1, high.news.coef,1.1:3.1, high.news.coef+high.news.se*1.96, col=c("dodgerblue","firebrick","springgreen4"), lwd=3)
segments(1.1:3.1, high.news.coef,1.1:3.1, high.news.coef-high.news.se*1.96, col=c("dodgerblue","firebrick","springgreen4"), lwd=3)

legend("topleft", c("Low News Democrat", "High News Democrat", "Low News Republican", "High News Republican","Low News Change in Cases/1000", "High News Change in Cases/1000"),
       col=c("dodgerblue","dodgerblue", "firebrick","firebrick","springgreen4","springgreen4")   , pch=c(16,17), cex=.5, bg="white")

axis(side=2, at=seq(-0.6,0.6,0.3))

abline(h=0, lty=2)
dev.off()




#Are Right Wing News Watchers More Affected?

table(dat$right.wing.media, dat$SMDate)

dates <- sort(unique(dat$SMDate[!is.na(dat$is.worried)]))
dat$time.trend <- NA
for(i in 1:length(dates)){
  dat$time.trend[dat$SMDate==dates[i]] <- i
}
dat$time.trend2 <- dat$time.trend^2

m.news <- felm(is.worried ~time.trend + time.trend2+ dem*right.wing.media + rep*right.wing.media + case.change*right.wing.media + 
                 female + income_30_75 + income_75_150 + income_above150 +
                 age_25_34 + age_35_44 + age_45_54 + age_55_64 + age_65_p +
                 other_race + black + asian + hispanic + 
                 hs + associates_somecollege + college + postgrad
               | state | 0 |  state, data=dat,
               weights=dat$weight_state_weekly)
summary(m.news)

levels <- seq(0,1,.001)

rep.margins <- Rmargins(m.news, "rep", "right.wing.media:rep", vcov(m.news), levels=levels)
dem.margins <- Rmargins(m.news, "dem", "dem:right.wing.media", vcov(m.news), levels=levels)




pdf(file="Figures/S17.pdf", height=6, width=8)
plot(levels, rep.margins[,2], xlim=c(0,1), ylim=c(-.5,.5), col="firebrick4", main="",
     xlab="Proportion of News Sources Right Wing", ylab="Coefficient Predicting P(Worried about COVID-19)", axes=F, type="l", lwd=3)
lines(levels, rep.margins[,4], col="firebrick4", lty=2)
lines(levels, rep.margins[,5], col="firebrick4", lty=2)


points(levels, dem.margins[,2], xlim=c(0,1), ylim=c(-1.5,1.5), col="dodgerblue4", type="l", lwd=3)
lines(levels,  dem.margins[,4], col="dodgerblue4", lty=2)
lines(levels,  dem.margins[,5], col="dodgerblue4", lty=2)


axis(side=2, at=seq(-.5,.5,.25))
axis(side=1, at=seq(0,1,.2))

abline(h=0, lty=3)

dev.off()


