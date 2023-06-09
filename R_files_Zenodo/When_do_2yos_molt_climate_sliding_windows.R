##Two-year-old timing of molt sliding windows

##This script comes after when do 2-year-olds molt script. 

library(here)
library(dplyr)
library(reshape2)
library(lubridate)
library(ggplot2)
library(cowplot)
library(viridis)
library(icenReg)
library(svMisc)
library(climwin)
library(survival)
library(survminer)

####Load climate data####

##Load rainfall data
rain = read.csv(here::here("Input files","rain_2004-2018.csv"))
rain$Date = as.Date(rain$Date,"%m/%d/%y")

##Load mintemp data
mintemp = read.csv(here::here("Input files","mintemp_1994-2020.csv"))
mintemp$Date = as.Date(mintemp$Date,"%m/%d/%y")

##Load maxtemp data
maxtemp = read.csv(here::here("Input files","maxtemp_1994-2020.csv"))
maxtemp$Date = as.Date(maxtemp$Date,"%m/%d/%y")

##Load NDVI data
ndvi = read.csv(here::here("Input files","Landsat7_all_years.csv"))
ndvi$Date = as.Date(ndvi$Date,"%m/%d/%y")

##Load molt date dataframe
pag2 = read.csv(here::here("Output files","twoyo_moltdates_dataframe.csv"))
pag2$Date.Created = as.Date(pag2$Date.Created)


###Plot rainfall
rain.plot = rain %>% filter(Year>=2014,Year<=2019)
ggplot(rain.plot,aes(x=jDate,y=amountcm)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot()




####Methods for the relative window analysis:
#For this analysis only including males that were seen during molt because for those that were not seen during molt we really don't know 
#when they molted. Some could have stayed bright from the previous year. Including those males works for previous analyses because the 
#covariates I was interested in for those analyses were factors. Here, using a covariate that changes by day - time dependent covariate - 
#so the event date should really be known or at least be close to known for the model to be accurate. 

##From timdep vignette in survival package:
#One of the strengths of the Cox model is its ability to encompass covariates that change over time. 
#The practical reason that time-dependent covariates work is based on the underlying way in which the Cox model works: 
#at each event time the program compares the current covariate values of the subject who had the event to the current values of 
#all others who were at risk at that time. One can think of it as a lottery model, where at each death time there is a drawing 
#to decide which subject “wins” the event. Each subject’s risk score exp(Xβ) determines how likely they are to win, e.g., how many “tickets” 
#they have purchased for the drawing. The model tries to assign a risk score to each subject that best predicts the outcome of each drawing based on:
#􏰀 The risk set: which subjects are present for each event; the set of those able to “win the prize”.
#􏰀 The covariate values of each subject just prior to the event time.

##Randomizations for these models to test that results were not false positives mix up molt dates across individuals but keep 
#other fixed effects associated with the same molt dates. So paired.prev and orn.prev stay with their molt date, but all get 
#assigned to new individuals and thus, new years. 



####Set up relative window dataframe####

#Get only males that were seen during molt
pag2.sdm = pag2 %>% filter(molt.time1>0)
pag2.sdm$Year = as.factor(pag2.sdm$Year)

##Get molt date as an actual date value
pag2.sdm$molt.date.date = as.Date(NA)
for (i in 1:nrow(pag2.sdm)) {
  #Get true year for the male and jan 1st of that year
  tyear = as.integer(as.character(pag2.sdm$Year[i])) - 1
  jan1 = as.Date(paste(tyear,"01-01",sep="-"))
  #Get molt date for that year
  pag2.sdm$molt.date.date[i] = as.Date(jan1 + pag2.sdm$molt.date[i]-1) #Have to subtract one to get true refday - can check with yday(refday.year)
}

##Set up dataframe for time dependent covariates - allows rainfall value to differ by each period
pag2.sdm.tdc = expand.grid(Bird = pag2.sdm$Bird,day2 = 136:274) #earliest day = 136, latest = 274
pag2.sdm.tdc = merge(pag2.sdm.tdc,pag2.sdm,by="Bird")
pag2.sdm.tdc = pag2.sdm.tdc %>% filter(molt.date>=day2) #Restrict to periods before each bird's molt
pag2.sdm.tdc$day1 = pag2.sdm.tdc$day2 - 1

##Get status variable
pag2.sdm.tdc$Status = NA
for(i in 1:nrow(pag2.sdm.tdc)) {
  if(pag2.sdm.tdc$day2[i]==pag2.sdm.tdc$molt.date[i]) {pag2.sdm.tdc$Status[i]=1} else {pag2.sdm.tdc$Status[i]=0}
}
pag2.sdm.tdc = pag2.sdm.tdc %>% arrange(Bird,day2)

##Get real day in pag2.sdm.tdc to help with getting climate value
pag2.sdm.tdc$day2.date = as.Date(NA)
for (i in 1:nrow(pag2.sdm.tdc)) {
  #Get true year for the male and jan 1st of that year
  tyear = as.integer(as.character(pag2.sdm.tdc$Year[i])) - 1
  jan1 = as.Date(paste(tyear,"01-01",sep="-"))
  #Get date for that year
  pag2.sdm.tdc$day2.date[i] = as.Date(jan1 + pag2.sdm.tdc$day2[i]-1) #Have to subtract one to get true refday - can check with yday(refday.year)
}





####Baseline relative model

##Get model
cph.rel = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~Paired.prev + Ornamented.prev,data=pag2.sdm.tdc)
summary(cph.rel)
plot(survfit(cph.rel),fun="event")

##Test proportional hazards assumption
test.ph = cox.zph(cph.rel)
test.ph
ggcoxzph(test.ph)
#Looks good

##Test influential observations
ggcoxdiagnostics(cph.rel, type = "deviance", linear.predictions = F)
#Looks ok?





####Custom relative sliding window####

# winopen = 10
# winclose = 0
# climate = rain
# climate.type = "rain.sum"
# exclude = 1

slidingwin_rel_custom = function(pag2.sdm.tdc,climate,climate.type,winopen,winclose,exclude) {
  
  #Get day ranges based on window open and window close
  windows = expand.grid(upper=seq(winopen,winclose),lower=seq(winopen,winclose))
  windows$diff = windows$upper-windows$lower #Get difference between windows
  windows = windows[which(windows$diff>=exclude),] #remove windows shorter than the exclude value
  
  #Get aic value of base model
  aic.base = AIC(coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~Paired.prev + Ornamented.prev,data=pag2.sdm.tdc))
  
  #Make output dataframes
  output = data.frame(matrix(ncol=7,nrow=nrow(windows)))
  colnames(output) = c("WindowOpen","WindowClose","ModelBeta","deltaAIC","Furthest","Closest","Climate.Type")
  
  #Get climate data into correct format
  if(climate.type=="rain.sum") {climate.data = climate %>% select(Date,amountcm)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="mintemp.mean") {climate.data = climate %>% select(Date,mintemp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="maxtemp.mean") {climate.data = climate %>% select(Date,maxtemp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="windspeed.mean") {climate.data = climate %>% select(Date,Windspeed)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="NDVI.mean") {climate.data = climate %>% select(Date,NDVI)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="residuals") {climate.data = climate %>% select(Date,residual)
  colnames(climate.data)[2]="climate.value"}
  
  #For each possible window
  for (i in 1:nrow(windows)) {
    #Get list of day2 dates - many repated dates in this analysis method
    pag2.sdm.tdc.dates = pag2.sdm.tdc %>% distinct(day2.date)
    #Get a climate variable for each row in pag2.sdm.tdc
    pag2.sdm.tdc.dates$climate.datahstat = NA
    
    #For each row in pag2.sdm.tdc.dates
    for (h in 1:nrow(pag2.sdm.tdc.dates)) {
      #Get start date of window
      winopenh = (pag2.sdm.tdc.dates$day2.date[h]-windows$upper[i])
      #Get end date of window
      wincloseh = (pag2.sdm.tdc.dates$day2.date[h]-windows$lower[i])
      #Subset climate data to weeks within that window
      climate.datah = climate.data[which(climate.data$Date>=winopenh & climate.data$Date<=wincloseh),]
      #Get the stats of the climate values in that window
      if(climate.type=="rain.sum") {climate.datahstat = sum(climate.datah$climate.value)} else {
        climate.datahstat = mean(climate.datah$climate.value)}
      #Add the climate value to the pag2.sdm.tdc.dates dataframe
      pag2.sdm.tdc.dates$climate.datahstat[h] = climate.datahstat
    }
    
    #Merge the climate values back into the full dataframe
    pag2.sdm.tdc.c = merge(pag2.sdm.tdc,pag2.sdm.tdc.dates,by="day2.date")
    
    #Now that each observation has a climate variable for the given window, get the model
    climate.m = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~climate.datahstat + paired.prev + Ornamented.prev,data=pag2.sdm.tdc.c)
    
    #Get AIC value
    aic.climate = AIC(climate.m)
    
    #Add values to output_lin dataframe
    output$WindowOpen[i] = windows$upper[i]
    output$WindowClose[i] = windows$lower[i]
    output$ModelBeta[i] = coef(climate.m)[1]
    output$deltaAIC[i] = aic.climate - aic.base
    output$Furthest[i] = winopen
    output$Closest[i] = winclose
    output$Climate.Type[i] = climate.type
    
    progress(i,max.value = nrow(windows))
    Sys.sleep(0.01)
    if (i == nrow(windows)) message("Done!")
  }
  #Order dataframes by lowest AIC values compared to baseline
  output = output[order(output$deltaAIC),]
  
  #Output dataframes
  return(output)
  }




####Custom single sliding window####

winopen = 10
winclose = 0
climate = rain
climate.type = "rain.sum"
exclude = 1

singlewin_rel_custom = function(pag2.sdm.tdc,climate,climate.type,winopen,winclose) {
  
  #Get aic value of base model
  aic.base = AIC(coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~Paired.prev + Ornamented.prev,data=pag2.sdm.tdc))
  
  #Get climate data into correct format
  if(climate.type=="rain.sum") {climate.data = climate %>% select(Date,amountcm)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="mintemp.mean") {climate.data = climate %>% select(Date,mintemp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="maxtemp.mean") {climate.data = climate %>% select(Date,maxtemp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="windspeed.mean") {climate.data = climate %>% select(Date,Windspeed)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="NDVI.mean") {climate.data = climate %>% select(Date,NDVI)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="residuals") {climate.data = climate %>% select(Date,residual)
  colnames(climate.data)[2]="climate.value"}
  
  #Get a climate data column
  pag2.sdm.tdc$climate.datahstat = NA
    
  #For each row in pag2.sdm.tdc
    for (h in 1:nrow(pag2.sdm.tdc)) {
      #Get start date of window
      winopenh = (pag2.sdm.tdc$day2.date[h]-winopen)
      #Get end date of window
      wincloseh = (pag2.sdm.tdc$day2.date[h]-winclose)
      #Subset climate data to weeks within that window
      climate.datah = climate.data[which(climate.data$Date>=winopenh & climate.data$Date<=wincloseh),]
      #Get the stats of the climate values in that window
      if(climate.type=="rain.sum") {climate.datahstat = sum(climate.datah$climate.value)} else {
        climate.datahstat = mean(climate.datah$climate.value)}
      #Add the climate value to the pag2.sdm.tdc dataframe
      pag2.sdm.tdc$climate.datahstat[h] = climate.datahstat
    }
    
    #Now that each observation has a climate variable for the given window, get the model
    climate.m = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~climate.datahstat + Paired.prev + Ornamented.prev,data=pag2.sdm.tdc)
    
    #Get delta-AIC value
    deltaAIC = AIC(climate.m) - aic.base

    #Output dataframes
    return(list(climate.m,deltaAIC,pag2.sdm.tdc))
}





####Rainfall sliding window####

##Run with linear term, quad and linear were the nearly the same window and AIC values

##Run the sliding window
#pag2.sdm.tdc.rain.300.0 = slidingwin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=rain,climate.type = "rain.sum",winopen = 300,winclose = 0,exclude = 30)

##Save the output
#write.csv(pag2.sdm.tdc.rain.300.0,here::here("Output files","pag2_rain_300_0.csv"),row.names = F)

##Read in output
pag2.sdm.tdc.rain.300.0 = read.csv(here::here("Output files","pag2_rain_300_0.csv"))

##Look at results
head(pag2.sdm.tdc.rain.300.0)
colnames(pag2.sdm.tdc.rain.300.0)[4] = "deltaAICc"
plotdelta(pag2.sdm.tdc.rain.300.0)

##Run single window
pag2.sdm.tdc.rain.300.0.s = singlewin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=rain,climate.type = "rain.sum",winopen = 216,winclose = 79)
summary(pag2.sdm.tdc.rain.300.0.s[[1]])
pag2.sdm.tdc.rain.300.0.s[[2]]

#Get dataset 
pag2.sdm.tdc.rain.300.0.df = pag2.sdm.tdc.rain.300.0.s[[3]]



####Rainfall randomizations

##Load randomization data from runs on cluster
rain.rand = read.csv(here::here("Output files","2yo_moltdate_rel_rainfall_300_0_randomizations_30exclude_v2.csv"))

##Get p-value
pvalue(dataset=pag2.sdm.tdc.rain.300.0,datasetrand=rain.rand,metric="AIC",sample.size=40)

##Plot randomizations
ggplot(data=rain.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.rain.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")




####Plot rainfall graphs

##What does the rainfall value at molt look like for each individual? - takes a bit to plot
#ggplot(data=pag2.sdm.tdc.rain.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Bird")


##Plot by year
ggplot(data=pag2.sdm.tdc.rain.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
 theme_cowplot()
ggplot(data=pag2.sdm.tdc.rain.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) + theme_cowplot()


##Plot rainfall graphs with molt dates
ggplot(rain.plot,aes(x=jDate,y=amountcm)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() 
ggplot(rain.plot,aes(x=jDate,y=amountcm)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)


##Show which rainfall periods are covered by each individual's 216-79 window if day 0 = molt date
#This section is a bit of a train wreck, but might be useful for visualization
pag2.sdm = pag2.sdm %>% arrange(Year,molt.date)
pag2.sdm$plot.y = seq(from=10.125,to=15,by=0.125)
pag2.sdm$upper = pag2.sdm$molt.date.date-216
pag2.sdm$lower = pag2.sdm$molt.date.date-79
pag2.sdm.plot = melt(pag2.sdm,id.vars = 1:37)

#Shift years and window dates for plotting
rain.plot2 = rain.plot
rain.plot2$Year2 = NA
for (i in 1:nrow(rain.plot2)) {
  if(month(rain.plot2$Date[i])<11) {rain.plot2$Year2[i] = rain.plot2$Year[i]} else {rain.plot2$Year2[i] = rain.plot2$Year[i] + 1
  rain.plot2$jDate[i] =0-(365-rain.plot2$jDate[i])}
}
pag2.sdm.plot$Year2 = pag2.sdm.plot$Year
pag2.sdm.plot$jdate = yday(pag2.sdm.plot$value)
pag2.sdm$Year2 = pag2.sdm$Year
for (i in 1:nrow(pag2.sdm.plot)) {if(pag2.sdm.plot$jdate[i]>250) {pag2.sdm.plot$jdate[i]=0-(365-pag2.sdm.plot$jdate[i])}}
rain.plot2 = rain.plot2 %>% filter(!Year2 %in% c("2014","2020"))

#Plot
ggplot(rain.plot2,aes(x=jDate,y=amountcm)) + geom_line() + facet_wrap(facets="Year2") + 
  geom_line(data=pag2.sdm.plot,aes(x=jdate,y=plot.y,group=Bird),size=0.25) + theme_cowplot() +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)



###Plot model predictions - hazard ratio
#Basic plot
plot(survfit(pag2.sdm.tdc.rain.300.0.s[[1]]),fun="cumhaz")

#Get prediction dataframe
rain.predict = expand.grid(climate.datahstat=seq(1:80),Paired.prev="Yes",Ornamented.prev="Brown")

#Predict values, risk = hazard ratio - which is exp(linear predictor (lp))
rain.predict2 = predict(pag2.sdm.tdc.rain.300.0.s[[1]],newdata=rain.predict,type="risk",se.fit=T) 
rain.predict$hazard.ratio = rain.predict2[[1]] #Hazard ratio
rain.predict$hazard.ratio.se = rain.predict2[[2]] #Standard error
rain.predict$CI = rain.predict$hazard.ratio.se *1.96 #Get CI
rain.predict$CI.upper = rain.predict$hazard.ratio + rain.predict$CI #Get upper
rain.predict$CI.lower = rain.predict$hazard.ratio - rain.predict$CI #Get lower

#Plot predicted hazard ratio against rainfall
ggplot(data=rain.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Cumulative rainfall (cm) from 216 to 79 days\\nbefore a date in the molting period") + ylab("Hazard ratio")

##Read in rainfall plot from bright by july script and plot together with hazard ratio
A = readRDS(here::here("Plots","bright_by_july.rds"))
B = ggplot(data=rain.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Cumulative rainfall (cm) from 216 to 79 days\\nbefore a date in the molting period") + ylab("Hazard ratio")
ggarrange(A,B,align="h")




####Detrend mintemp and maxtemp####

#Like breeding start date in non-breeding paper, interested in anomalies compared to averages calculated using
#all years. So get average loess line of all years, then detrend relative to that. 

#Think I switched to detrended temperature windows because they made more sense that way in the context of a relative sliding window
#with time dependent values. Regular temperature would end up saying something like higher temperatures led to increased likelihood of 
#molt. But what is higher in relation to? Using the anomaly gives a value relative to longterm average. Using 2002-2019 gives a very 
#smooth line, where as calculating a trend for only 4 years (what we have molt data for) could create a wavy line that's not reflective
#of the longterm average. 

#Get an average min and maxtemp for 2002 - 2019
mintemp = mintemp %>% filter(Date>="2002-01-01",Date<="2019-12-31")
maxtemp = maxtemp %>% filter(Date>="2002-01-01",Date<="2019-12-31")
mintemp$all.fitted = loess(mintemp~jdate,data=mintemp,span=0.5)$fitted
maxtemp$all.fitted = loess(maxtemp~jdate,data=maxtemp,span=0.5)$fitted

#Show what that looks like
mintemp$Year = as.factor(mintemp$Year)
ggplot(mintemp,aes(x=jdate,y=mintemp)) + geom_line(aes(group=Year,color=Year)) + 
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="black") + theme_cowplot()
maxtemp$Year = as.factor(maxtemp$Year)
ggplot(maxtemp,aes(x=jdate,y=maxtemp)) + geom_line(aes(group=Year,color=Year)) + 
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="black") + theme_cowplot()

#Get residuals for each year
mintemp$residual = loess(mintemp~jdate,data=mintemp,span=0.5)$residual
maxtemp$residual = loess(maxtemp~jdate,data=maxtemp,span=0.5)$residual
#If you order by jdate you can see that each date gets calculated independently, not repeating values like "fitted". 





####Detrended mintemp sliding window####

##Run the sliding window
#pag2.sdm.tdc.mintemp.300.0 = slidingwin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=mintemp,climate.type = "residuals",winopen = 300,winclose = 0,exclude = 30)

##Save the output
#write.csv(pag2.sdm.tdc.mintemp.300.0,here::here("Output files","pag2_dtmintemp_300_0.csv"),row.names = F)

##Read in output
pag2.sdm.tdc.mintemp.300.0 = read.csv(here::here("Output files","pag2_dtmintemp_300_0.csv"))

##Look at results
head(pag2.sdm.tdc.mintemp.300.0)
colnames(pag2.sdm.tdc.mintemp.300.0)[4] = "deltaAICc"
plotdelta(pag2.sdm.tdc.mintemp.300.0)

##Run single window
pag2.sdm.tdc.mintemp.300.0.s = singlewin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=mintemp,climate.type = "residuals",winopen = 273,winclose = 125)
summary(pag2.sdm.tdc.mintemp.300.0.s[[1]])
pag2.sdm.tdc.mintemp.300.0.s[[2]]

#Get dataset 
pag2.sdm.tdc.mintemp.300.0.df = pag2.sdm.tdc.mintemp.300.0.s[[3]]


###Mintemp randomizations

##Load randomization data from runs on cluster
mintemp.rand = read.csv(here::here("Output files","2yo_moltdate_rel_mintemp_300_0_randomizations_30exclude.csv"))

##Get p-value
pvalue(dataset=pag2.sdm.tdc.mintemp.300.0,datasetrand=mintemp.rand,metric="AIC",sample.size=40)

##Plot randomizations
ggplot(data=mintemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.mintemp.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")



####Plot mintemp graphs

##What does the mintemp value at molt look like for each individual? - takes a bit to plot
ggplot(data=pag2.sdm.tdc.mintemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Bird")

##Plot by year
ggplot(data=pag2.sdm.tdc.mintemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  theme_cowplot()
ggplot(data=pag2.sdm.tdc.mintemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) + theme_cowplot()

##Plot with rainfall and mean mintemp values
mintemp.plot = mintemp %>% filter(Year %in% c("2015","2016","2017","2018","2019"))
ggplot(data=mintemp.plot,aes(x=jdate,y=mintemp)) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) +
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="yellow") + 
  geom_line(data=rain.plot,aes(x=jDate,y=amountcm))

##Plot residuals
ggplot(data=mintemp.plot,aes(x=jdate,y=residual)) + geom_line() + facet_wrap(facets="Year") +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)

##Are rainfall and mintemp windows correlated with one another? 
rain.mintemp = pag2.sdm.tdc.mintemp.300.0.df %>% select(mintemp = climate.datahstat)
rain.mintemp = rain.mintemp %>% mutate(rain = pag2.sdm.tdc.rain.300.0.df$climate.datahstat)
ggplot(data=rain.mintemp,aes(x=mintemp,y=rain)) + geom_point() + geom_smooth(method="lm") + theme_cowplot() +
  xlab("Mean minimum temperature anomaly from\\n273-125 days before date") + ylab("Total rainfall between 216-79 days\\nbefore date")
cor.test(rain.mintemp$mintemp,rain.mintemp$rain)
#Yes but not perfectly, so higher than average mintemps usually lead to rainfall but how much rainfall not necessarily proportional



###Plot model predictions - hazard ratio
#Basic plot
plot(survfit(pag2.sdm.tdc.mintemp.300.0.s[[1]]),fun="cumhaz")

#Get prediction dataframe
mintemp.predict = expand.grid(climate.datahstat=seq(from=-0.4, to=1.1,by=0.05),Paired.prev="Yes",Ornamented.prev="Brown")

#Predict values, risk = hazard ratio - which is exp(linear predictor (lp))
mintemp.predict2 = predict(pag2.sdm.tdc.mintemp.300.0.s[[1]],newdata=mintemp.predict,type="risk",se.fit=T) 
mintemp.predict$hazard.ratio = mintemp.predict2[[1]] #Hazard ratio
mintemp.predict$hazard.ratio.se = mintemp.predict2[[2]] #Standard error
mintemp.predict$CI = mintemp.predict$hazard.ratio.se *1.96 #Get CI
mintemp.predict$CI.upper = mintemp.predict$hazard.ratio + mintemp.predict$CI #Get upper
mintemp.predict$CI.lower = mintemp.predict$hazard.ratio - mintemp.predict$CI #Get lower

#Plot predicted hazard ratio against mintemp
ggplot(data=mintemp.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean minimum temperature (C) from 273 to 125 days\\nbefore a date in the molting period") + ylab("Hazard ratio")

#Above average minimum temperatures led to rainfall which led to increased likelihood of molt.









####Detrended maxtemp sliding window####

##Run the sliding window
#pag2.sdm.tdc.maxtemp.300.0 = slidingwin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=maxtemp,climate.type = "residuals",winopen = 300,winclose = 0,exclude = 30)

##Save the output
#write.csv(pag2.sdm.tdc.maxtemp.300.0,here::here("Output files","pag2_dtmaxtemp_300_0.csv"),row.names = F)

##Read in output
pag2.sdm.tdc.maxtemp.300.0 = read.csv(here::here("Output files","pag2_dtmaxtemp_300_0.csv"))

##Look at results
head(pag2.sdm.tdc.maxtemp.300.0)
colnames(pag2.sdm.tdc.maxtemp.300.0)[4] = "deltaAICc"
plotdelta(pag2.sdm.tdc.maxtemp.300.0)

##Run single window
pag2.sdm.tdc.maxtemp.300.0.s = singlewin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=maxtemp,climate.type = "residuals",winopen = 115,winclose = 15)
summary(pag2.sdm.tdc.maxtemp.300.0.s[[1]])
pag2.sdm.tdc.maxtemp.300.0.s[[2]]

#Get dataset 
pag2.sdm.tdc.maxtemp.300.0.df = pag2.sdm.tdc.maxtemp.300.0.s[[3]]


###Maxtemp randomizations

##Load randomization data from runs on cluster
maxtemp.rand = read.csv(here::here("Output files","2yo_moltdate_rel_maxtemp_300_0_randomizations_30exclude.csv"))

##Get p-value
pvalue(dataset=pag2.sdm.tdc.maxtemp.300.0,datasetrand=maxtemp.rand,metric="AIC",sample.size=40)

##Plot randomizations
ggplot(data=maxtemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.maxtemp.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")




####Plot maxtemp graphs

##What does the maxtemp value at molt look like for each individual? - takes a bit to plot
ggplot(data=pag2.sdm.tdc.maxtemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Bird")


##Plot by year
ggplot(data=pag2.sdm.tdc.maxtemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  theme_cowplot()
ggplot(data=pag2.sdm.tdc.maxtemp.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) + theme_cowplot()


##Plot maxtemp graphs with molt dates
maxtemp.plot = maxtemp %>% filter(Year %in% c("2015","2016","2017","2018","2019"))
ggplot(maxtemp.plot,aes(x=jdate,y=maxtemp)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() 
ggplot(maxtemp.plot,aes(x=jdate,y=maxtemp)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)


##Plot with rainfall and mean maxtemp values
ggplot(data=maxtemp.plot,aes(x=jdate,y=maxtemp)) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) +
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="yellow") + 
  geom_line(data=rain.plot,aes(x=jDate,y=amountcm))

##Plot residuals
ggplot(data=maxtemp.plot,aes(x=jdate,y=residual)) + geom_line() + facet_wrap(facets="Year") +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)

##Are rainfall and maxtemp windows correlated with one another? 
rain.maxtemp = pag2.sdm.tdc.maxtemp.300.0.df %>% select(maxtemp = climate.datahstat)
rain.maxtemp = rain.maxtemp %>% mutate(rain = pag2.sdm.tdc.rain.300.0.df$climate.datahstat)
ggplot(data=rain.maxtemp,aes(x=rain,y=maxtemp)) + geom_point() + geom_smooth(method="lm") + theme_cowplot() +
  xlab("Total rainfall between 216-79 days\\nbefore date") + ylab("Mean maximum temperature anomaly from\\n115 to 15 days before date")
cor.test(rain.maxtemp$maxtemp,rain.maxtemp$rain)
#Yes. Rainfall leads to cooler than average temperatures since rainfall window is 216-79 days before and maxtemp window is 115-15 before



###Plot model predictions - hazard ratio
#Basic plot
plot(survfit(pag2.sdm.tdc.maxtemp.300.0.s[[1]]),fun="cumhaz")

#Get prediction dataframe
maxtemp.predict = expand.grid(climate.datahstat=seq(from=-0.4,to=1.5,by=0.05),Paired.prev="Yes",Ornamented.prev="Brown")

#Predict values, risk = hazard ratio - which is exp(linear predictor (lp))
maxtemp.predict2 = predict(pag2.sdm.tdc.maxtemp.300.0.s[[1]],newdata=maxtemp.predict,type="risk",se.fit=T) 
maxtemp.predict$hazard.ratio = maxtemp.predict2[[1]] #Hazard ratio
maxtemp.predict$hazard.ratio.se = maxtemp.predict2[[2]] #Standard error
maxtemp.predict$CI = maxtemp.predict$hazard.ratio.se *1.96 #Get CI
maxtemp.predict$CI.upper = maxtemp.predict$hazard.ratio + maxtemp.predict$CI #Get upper
maxtemp.predict$CI.lower = maxtemp.predict$hazard.ratio - maxtemp.predict$CI #Get lower

#Plot predicted hazard ratio against maxtemp
ggplot(data=maxtemp.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean maximum temperature (C) from 115 to 15 days\\nbefore a date in the molting period") + ylab("Hazard ratio")

#Cooler than average maximum temperatures lead to increased likelihood of molt into ornamented plumage








####NDVI sliding window####

##Run the sliding window
#pag2.sdm.tdc.ndvi.300.0 = slidingwin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=ndvi,climate.type = "NDVI.mean",winopen = 300,winclose = 0,exclude = 30)
#Ran 400-0, exclude 300 and deltaAIC dies off quickly after 300-0 window at -19. 

##Save the output
#write.csv(pag2.sdm.tdc.ndvi.300.0,here::here("Output files","pag2_ndvi_300_0.csv"),row.names = F)

##Read in output
pag2.sdm.tdc.ndvi.300.0 = read.csv(here::here("Output files","pag2_ndvi_300_0.csv"))

##Look at results
head(pag2.sdm.tdc.ndvi.300.0)
colnames(pag2.sdm.tdc.ndvi.300.0)[4] = "deltaAICc"
plotdelta(pag2.sdm.tdc.ndvi.300.0)

##Run single window
pag2.sdm.tdc.ndvi.300.0.s = singlewin_rel_custom(pag2.sdm.tdc = pag2.sdm.tdc,climate=ndvi,climate.type = "NDVI.mean",winopen = 107,winclose = 0)
summary(pag2.sdm.tdc.ndvi.300.0.s[[1]])
pag2.sdm.tdc.ndvi.300.0.s[[2]]

#Get dataset 
pag2.sdm.tdc.ndvi.300.0.df = pag2.sdm.tdc.ndvi.300.0.s[[3]]


###ndvi randomizations

##Load randomization data from runs on cluster
ndvi.rand = read.csv(here::here("Output files","2yo_moltdate_rel_ndvi_300_0_randomizations_30exclude.csv"))

##Get p-value
pvalue(dataset=pag2.sdm.tdc.ndvi.300.0,datasetrand=ndvi.rand,metric="AIC",sample.size=40)

##Plot randomizations
ggplot(data=ndvi.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.ndvi.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")



####Plot ndvi graphs

##What does the ndvi value at molt look like for each individual? - takes a bit to plot
ggplot(data=pag2.sdm.tdc.ndvi.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Bird")


##Plot by year
ggplot(data=pag2.sdm.tdc.ndvi.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  theme_cowplot()
ggplot(data=pag2.sdm.tdc.ndvi.300.0.df,aes(x=day2,y=climate.datahstat)) + geom_point(size=0.5) + geom_line() + facet_wrap(facets="Year") + 
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25) + theme_cowplot()


##Plot ndvi graphs with molt dates
ndvi$Year = year(ndvi$Date) + 1
ndvi$jdate = yday(ndvi$Date)
ndvi.plot = ndvi %>% filter(Year %in% c("2015","2016","2017","2018","2019"))
ggplot(ndvi.plot,aes(x=jdate,y=NDVI)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() 
ggplot(ndvi.plot,aes(x=jdate,y=NDVI)) + geom_line() + facet_wrap(facets = "Year") + theme_cowplot() +
  geom_vline(data=pag2.sdm,aes(xintercept = molt.date),color="red",size=0.25)


##Are rainfall and ndvi windows correlated with one another? 
rain.ndvi = pag2.sdm.tdc.ndvi.300.0.df %>% select(ndvi = climate.datahstat)
rain.ndvi = rain.ndvi %>% mutate(rain = pag2.sdm.tdc.rain.300.0.df$climate.datahstat)
ggplot(data=rain.ndvi,aes(x=rain,y=ndvi)) + geom_point() + geom_smooth(method="lm") + 
  xlab("Total rainfall between 216-79 days\\nbefore date") + ylab("Mean NDVI score from 125 to 0 days\\nbefore a date") + theme_cowplot()
cor.test(rain.ndvi$ndvi,rain.ndvi$rain)
#Yes. Rainfall leads to higher NDVI since rainfall window is 216-79 days before and ndvi window is 125-0 before



###Plot model predictions - hazard ratio
#Basic plot
plot(survfit(pag2.sdm.tdc.ndvi.300.0.s[[1]]),fun="cumhaz")

#Get prediction dataframe
ndvi.predict = expand.grid(climate.datahstat=seq(from=0.53,to=0.68,by=0.005),Paired.prev="Yes",Ornamented.prev="Brown")

#Predict values, risk = hazard ratio - which is exp(linear predictor (lp))
ndvi.predict2 = predict(pag2.sdm.tdc.ndvi.300.0.s[[1]],newdata=ndvi.predict,type="risk",se.fit=T) 
ndvi.predict$hazard.ratio = ndvi.predict2[[1]] #Hazard ratio
ndvi.predict$hazard.ratio.se = ndvi.predict2[[2]] #Standard error
ndvi.predict$CI = ndvi.predict$hazard.ratio.se *1.96 #Get CI
ndvi.predict$CI.upper = ndvi.predict$hazard.ratio + ndvi.predict$CI #Get upper
ndvi.predict$CI.lower = ndvi.predict$hazard.ratio - ndvi.predict$CI #Get lower

#Plot predicted hazard ratio against ndvi
ggplot(data=ndvi.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean NDVI score from 125 to 0 days\\nbefore a date in the molting period") + ylab("Hazard ratio")

#An NDVI score higher than 0.6 lead to increased likelihood of molt 







####Final models for long-term windows####

##Many of the climate variables are clearly connected - think it makes sense that because maxtemp follows rainfall, do not include maxtemp in 
#model, seems like it's a result of rainfall and rainfall is more likely to have more important consequences for molt

##Make sure dataframes still in same order
identical(pag2.sdm.tdc$day2,pag2.sdm.tdc.maxtemp.300.0.df$day2)
identical(pag2.sdm.tdc$day2,pag2.sdm.tdc.mintemp.300.0.df$day2)
identical(pag2.sdm.tdc$day2,pag2.sdm.tdc.rain.300.0.df$day2)
identical(pag2.sdm.tdc$day2,pag2.sdm.tdc.ndvi.300.0.df$day2)
identical(pag2.sdm.tdc$Bird,pag2.sdm.tdc.maxtemp.300.0.df$Bird)
identical(pag2.sdm.tdc$Bird,pag2.sdm.tdc.mintemp.300.0.df$Bird)
identical(pag2.sdm.tdc$Bird,pag2.sdm.tdc.rain.300.0.df$Bird)
identical(pag2.sdm.tdc$Bird,pag2.sdm.tdc.ndvi.300.0.df$Bird)

##Set up dataframe
pag2.climate = pag2.sdm.tdc %>% select(Bird,Fwnumber,day2,day1,Status,Year,Ornamented.prev,Paired.prev,day2.date)
pag2.climate$rain.216.79 = pag2.sdm.tdc.rain.300.0.df$climate.datahstat
pag2.climate$mintemp.273.125 = pag2.sdm.tdc.mintemp.300.0.df$climate.datahstat
pag2.climate$maxtemp.115.15 = pag2.sdm.tdc.maxtemp.300.0.df$climate.datahstat
pag2.climate$ndvi.125.0 = pag2.sdm.tdc.ndvi.300.0.df$climate.datahstat

##Run full model
climate.model = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~rain.216.79 + mintemp.273.125 + Paired.prev + Ornamented.prev,data=pag2.climate)
summary(climate.model)
AIC(climate.model) #190.85

##Run without Paired.prev
climate.model.nopp = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~rain.216.79 + mintemp.273.125 + Ornamented.prev,data=pag2.climate)
summary(climate.model.nopp)
AIC(climate.model.nopp) #190.5

##Run without Ornamented.prev
climate.model.noop = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~rain.216.79 + mintemp.273.125 + Paired.prev,data=pag2.climate)
summary(climate.model.noop)
AIC(climate.model.noop) #188.8

##Run without rainfall
climate.model.nor = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~ mintemp.273.125 + Paired.prev + Ornamented.prev,data=pag2.climate)
summary(climate.model.nor)
AIC(climate.model.nor) #194.5

##Run without mintemp
climate.model.nomt = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~rain.216.79 + Paired.prev + Ornamented.prev,data=pag2.climate)
summary(climate.model.nomt)
AIC(climate.model.nomt) #191.85

##Final model - only rainfall was important - all others AIC did not change by more than 2
climate.model2 = coxph(Surv(time=day1,time2=day2,event = Status,type="counting")~rain.216.79,data=pag2.climate)
summary(climate.model2)
AIC(climate.model2) #188.8

##Paper: Report LRT from model summary(X2=33.78,p<0.001)

##So does appear that higher than average mintemp is leading to rainfall, but rainfall is the most important predictor of the likelihood of molt
#Could run binomial model looking at the likelihood of rainfall on days throughout this period and minimum temperature sliding window. 

#So were the dry years cooler years? Plot minimum temperature by week for those years
mintemp.plot.week = mintemp.plot %>% mutate(week = week(Date))
mintemp.plot.week = mintemp.plot.week %>% group_by(Year,week) %>% summarise(mean.mintemp = mean(mintemp))
ggplot(data=mintemp.plot.week,aes(x=week,y=mean.mintemp)) + geom_line(aes(group=Year,color=Year),size=1) + scale_color_viridis(discrete = T) +
  theme_cowplot()
#No major obvious differences in temperature across years. So higher mintemp often leads to rainfall, but differences in mintemp among years
#probably not explaining variation in molt timing across years for 2yos. 


###Plot predictions for final model 
#Basic plot
plot(survfit(climate.model2),fun="cumhaz")

#Get prediction dataframe
final.predict = expand.grid(rain.216.79=seq(1:80))

#Predict values, risk = hazard ratio - which is exp(linear predictor (lp))
final.predict2 = predict(climate.model2,newdata=final.predict,type="risk",se.fit=T) 
final.predict$hazard.ratio = final.predict2[[1]] #Hazard ratio
final.predict$hazard.ratio.se = final.predict2[[2]] #Standard error
final.predict$CI = final.predict$hazard.ratio.se *1.96 #Get CI
final.predict$CI.upper = final.predict$hazard.ratio + final.predict$CI #Get upper
final.predict$CI.lower = final.predict$hazard.ratio - final.predict$CI #Get lower

#Plot predicted hazard ratio against rainfall
ggplot(data=final.predict,aes(x=rain.216.79,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Cumulative rainfall (cm) from 216 to 79 days\\nbefore a date in the molting period") + ylab("Hazard ratio")








####Plot deltaAIC plots and randomizations for supplemental materials####

##deltaAIC plots
A = plotdelta(pag2.sdm.tdc.rain.300.0)
B = plotdelta(pag2.sdm.tdc.mintemp.300.0)
C = plotdelta(pag2.sdm.tdc.maxtemp.300.0)
D = plotdelta(pag2.sdm.tdc.ndvi.300.0)

ggarrange(A,B,C,D,labels=c("A","B","C","D"))

##Randomization plots
A = ggplot(data=rain.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.rain.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")
B = ggplot(data=mintemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.mintemp.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")
C = ggplot(data=maxtemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.maxtemp.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")
D = ggplot(data=ndvi.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag2.sdm.tdc.ndvi.300.0$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAIC")

ggarrange(A,B,C,D,labels=c("A","B","C","D"))


##Plot correlation plots to rainfall
A = ggplot(data=rain.mintemp,aes(x=mintemp,y=rain)) + geom_point() + geom_smooth(method="lm") + theme_cowplot() +
  xlab("Mean minimum temperature anomaly from\\n273-125 days before date") + ylab("Total rainfall between\\n216-79 days before date")
B = ggplot(data=rain.maxtemp,aes(x=rain,y=maxtemp)) + geom_point() + geom_smooth(method="lm") + theme_cowplot() +
  xlab("Total rainfall between 216-79 days\\nbefore date") + ylab("Mean maximum temperature\\nanomaly from 115 to 15 days\\nbefore date")
C = ggplot(data=rain.ndvi,aes(x=rain,y=ndvi)) + geom_point() + geom_smooth(method="lm") + 
  xlab("Total rainfall between 216-79 days\\nbefore date") + ylab("Mean NDVI score from\\n125 to 0 days before a date") + theme_cowplot()

ggarrange(A,B,C,labels=c("A","B","C"))



##Plot hazard ratio plots for supplemental
A = ggplot(data=mintemp.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean minimum temperature anomaly (C) from\\n273 to 125 days before a date in the molting period") + ylab("Hazard ratio")
B = ggplot(data=maxtemp.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean maximum temperature anomaly (C) from\\n115 to 15 days before a date in the molting period") + ylab("Hazard ratio")
C = ggplot(data=ndvi.predict,aes(x=climate.datahstat,y=hazard.ratio)) + geom_line(size=1.5) + theme_cowplot() + 
  geom_ribbon(aes(ymin=CI.lower,ymax=CI.upper),alpha=0.3) + geom_hline(aes(yintercept=1),lty=2) + 
  xlab("Mean NDVI score from 125 to 0 days\\nbefore a date in the molting period") + ylab("Hazard ratio")
ggarrange(A,B,C,labels=c("A","B","C"))



