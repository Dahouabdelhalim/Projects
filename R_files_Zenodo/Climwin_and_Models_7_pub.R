#Climwin and Models script v7

#Use this script to run the models for the climate window analyses - determine which climate windows
#best predict the number of groups encountered by the focal group in an observation 

#If this script does not run, close the project file and re-open. One of the packages
#used in another script must overlap with one used in this script which can cause some 
#obvious errors that prevents the final model from being calculated.

#Initially ran these analyses with the climwin package, but found the climwin package could not look at 
#within-year variation with a quadratic function at the same time. So wrote my own sliding window function. 


library(here)

####Load data####
gbs = read.csv(here::here("Output files","gbssubmaxall.csv"))
rainall = read.csv(here::here("Output files","rainallyears.csv")) #Rain has slightly fewer rows bc only goes through Nov 2018
mintempall = read.csv(here::here("Output files","mintempallyears.csv"))
maxtempall = read.csv(here::here("Output files","maxtempallyears.csv"))
windall = read.csv(here::here("Output files","windallyears.csv"))
ndviall = read.csv(here::here("Output files","Landsat 8 NDVI yearly polygons noclouds.csv"))

#Set date as date 
gbs$Date = as.Date(gbs$Date,"%m/%d/%y")
gbs = gbs[order(gbs$Date),]
rainall$Date = as.Date(rainall$Date,"%m/%d/%y")
mintempall$Date = as.Date(mintempall$Date,"%m/%d/%y")
maxtempall$Date = as.Date(maxtempall$Date,"%m/%d/%y")
windall$Date = as.Date(windall$Date,"%m/%d/%y")
ndviall$Date = as.Date(ndviall$Date,"%m/%d/%y")

#Set factors
gbs$Sighting = as.factor(gbs$Sighting)
gbs$Observation = as.factor(gbs$Observation)
gbs$Year = as.factor(gbs$Year)
gbs$Focal.Group.Year = as.factor(gbs$Focal.Group.Year)
gbs$Community.Year = as.factor(gbs$Community.Year)
rainall$Year = as.factor(rainall$Year)
ndviall$Year = as.factor(ndviall$Year)
mintempall$Year = as.factor(mintempall$Year)
maxtempall$Year = as.factor(maxtempall$Year)

#Get number of groups encountered by focal group
gbs$encountered = gbs$Social.Groups.Total - 1

#Get rainfall in cm
rainall$amountcm = rainall$amount/10


#Get percentage of rainfall 
1-sum(is.na(rainall$Rainfall.amount..millimetres..x))/nrow(rainall) #Station 40186
1-sum(is.na(rainall$Rainfall.amount..millimetres..y))/nrow(rainall) #Station 40517
1-sum(is.na(rainall$Rainfall.amount..millimetres.))/nrow(rainall) #Station 40960


####Explore data for models####

library(ggplot2)
library(cowplot)

#Density plot of encountered groups
ggplot(gbs,aes(x=encountered,group=Year,color=Year)) + geom_density()

#Mean and variance
mean(gbs$encountered)
var(gbs$encountered)

###Encountered Correlations
##Dot plot of groups encountered by community
ggplot(data=gbs,aes(encountered,Community.Year)) + geom_point()

##Does the number of groups in a community influence the number of groups encountered? 
ggplot(data=gbs,aes(Groups.in.Community,encountered,group=Year)) + geom_point() + 
  geom_smooth(aes(color=Year),method="lm")
#Not really, but maybe slightly in 2017

##Does social group size of the focal group influence the number of groups encountered?
ggplot(data=gbs,aes(Social.Group.Size,encountered,group=Year)) + geom_point() + 
  geom_smooth(aes(color=Year),method="lm")
#Not really

##Does the number of sightings in an observation (number of sightings focal group was present)
#influence the number of groups encountered?
ggplot(data=gbs,aes(Focal.Sightings,encountered)) + geom_point() + 
  geom_smooth(aes(color=Year),method="lm")
#Yes - definitely add into the model



####Look at distributions####
par(mfrow=c(1,1))
hist(gbs$encountered)

library(car)
library(MASS)
#Fit distribution for nbinom
nbinom = fitdistr(gbs$encountered, "Negative Binomial")
qqp(gbs$encountered,"nbinom",size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#Fit distribution for poisson
poissond = fitdistr(gbs$encountered, "Poisson")
qqp(gbs$encountered,"pois",lambda=poissond$estimate[[1]])



####Baseline climate model####

#Center focal sightings on year
library(climwin)
gbs$Focal.Sightings.cent = wgdev(gbs$Focal.Sightings,gbs$Year)

library(lme4)
baseline = glmer(encountered~Focal.Sightings.cent + jdate + Year + (1|Community.Year/Focal.Group.Year),data=gbs,family=poisson)  
summary(baseline)

#Negative binomial gives about same results
summary(glmer.nb(encountered~Focal.Sightings.cent + jdate + Year + (1|Community.Year/Focal.Group.Year),data=gbs))




####Custom slidingwin function####

#Load required packages
library(lme4)
library(MuMIn)
library(climwin)
library(svMisc)
library(dplyr)


slidingwin_custom = function(observations,climate,climate.type,winopen,winclose,exclude) {
  
  #Get day ranges based on window open and window close
  windows = expand.grid(upper=seq(winopen,winclose),lower=seq(winopen,winclose))
  windows$diff = windows$upper-windows$lower #Get difference between windows
  windows = windows[which(windows$diff>=exclude),] #remove windows shorter than the exclude value
  
  #Get aic value of base model
  aic.base = AIC(glmer(encountered ~ Focal.Sightings.cent + jdate + Year + (1|Community.Year/Focal.Group.Year),
                       data=observations,family=poisson))
  
  #Make output dataframes
  output_lin = data.frame(matrix(ncol=7,nrow=nrow(windows)))
  colnames(output_lin) = c("WindowOpen","WindowClose","ModelBeta","deltaAIC","Furthest","Closest","Climate.Type")
  output_quad = data.frame(matrix(ncol=8,nrow=nrow(windows)))
  colnames(output_quad) = c("WindowOpen","WindowClose","ModelBeta","ModelBeta.quad","deltaAIC","Furthest","Closest","Climate.Type")
  
  #Get climate data into correct format
  if(climate.type=="rain.sum") {climate.data = climate %>% select(Date,amountcm)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="mintemp.mean") {climate.data = climate %>% select(Date,Min.Temp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="maxtemp.mean") {climate.data = climate %>% select(Date,Max.Temp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="windspeed.mean") {climate.data = climate %>% select(Date,Windspeed)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="NDVI.mean") {climate.data = climate %>% select(Date,NDVI)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="residuals") {climate.data = climate %>% select(Date,residuals)
  colnames(climate.data)[2]="climate.value"}
  
  #For each possible window
  for (i in 1:nrow(windows)) {
    #Get a climate variable for each row in observations
    observations$climate.datahstat = NA
    #For each row in observations
    for (h in 1:nrow(observations)) {
      #Get start date of window
      winopenh = (observations$Date[h]-windows$upper[i])
      #Get end date of window
      wincloseh = (observations$Date[h]-windows$lower[i])
      #Subset climate data to dates within that window
      climate.datah = climate.data[which(climate.data$Date>=winopenh & climate.data$Date<=wincloseh),]
      #Get the stats of the climate values in that window
      if(climate.type=="rain.sum") {climate.datahstat = sum(climate.datah$climate.value)} else {
        climate.datahstat = mean(climate.datah$climate.value)}
      #Add the climate value to the observations dataframe
      observations$climate.datahstat[h] = climate.datahstat
    }
    
    #Center the window's climate values on year since we're interested in within-year variation
    observations$climate.datahstat.cent = wgdev(observations$climate.datahstat,observations$Year)
    
    #Now that each observation has a climate variable for the given window, get the linear and quadratic models
    model.lin = glmer(encountered ~ Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + (1|Community.Year/Focal.Group.Year),
                      data=observations,family=poisson)
    model.quad = glmer(encountered ~ Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) + (1|Community.Year/Focal.Group.Year),
                       data=observations,family=poisson)
    
    #Get AIC values - using AIC instead of AICc because samples size is so large, shouldn't need AICc and results almost exactly the same
    aic.lin = AIC(model.lin)
    aic.quad = AIC(model.quad)
    
    #Add values to output_lin dataframe
    output_lin$WindowOpen[i] = windows$upper[i]
    output_lin$WindowClose[i] = windows$lower[i]
    output_lin$ModelBeta[i] = fixef(model.lin)[7]
    output_lin$deltaAIC[i] = aic.lin - aic.base
    output_lin$Furthest[i] = winopen
    output_lin$Closest[i] = winclose
    output_lin$Climate.Type[i] = climate.type
    
    #Add values to output_quad dataframe
    output_quad$WindowOpen[i] = windows$upper[i]
    output_quad$WindowClose[i] = windows$lower[i]
    output_quad$ModelBeta[i] = fixef(model.quad)[7]
    output_quad$ModelBeta.quad[i] = fixef(model.quad)[8]
    output_quad$deltaAIC[i] = aic.quad - aic.base
    output_quad$Furthest[i] = winopen
    output_quad$Closest[i] = winclose
    output_quad$Climate.Type[i] = climate.type
    
    progress(i,max.value = nrow(windows))
    Sys.sleep(0.01)
    if (i == nrow(windows)) message("Done!")
    
  }
  
  #Order dataframes by lowest AIC values compared to baseline
  output_lin = output_lin[order(output_lin$deltaAIC),]
  output_quad = output_quad[order(output_quad$deltaAIC),]
  
  #Output dataframes
  return(list(output_lin,output_quad))
  
}




####Custom single win function####
singlewin_custom = function(observations,climate,climate.type,winopen,winclose,func) {
  
  #Get climate data into correct format
  if(climate.type=="rain.sum") {climate.data = climate %>% select(Date,amountcm)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="mintemp.mean") {climate.data = climate %>% select(Date,Min.Temp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="maxtemp.mean") {climate.data = climate %>% select(Date,Max.Temp)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="windspeed.mean") {climate.data = climate %>% select(Date,Windspeed)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="NDVI.mean") {climate.data = climate %>% select(Date,NDVI)
  colnames(climate.data)[2]="climate.value"}
  if(climate.type=="residuals") {climate.data = climate %>% select(Date,residuals)
  colnames(climate.data)[2]="climate.value"}
  
    for (h in 1:nrow(observations)) {
      #Get start date of window
      winopenh = (observations$Date[h]-winopen)
      #Get end date of window
      wincloseh = (observations$Date[h]-winclose)
      #Subset climate data to dates within that window
      climate.datah = climate.data[which(climate.data$Date>=winopenh & climate.data$Date<=wincloseh),]
      #Get the stats of the climate values in that window
      if(climate.type=="rain.sum") {climate.datahstat = sum(climate.datah$climate.value)} else {
        climate.datahstat = mean(climate.datah$climate.value)}
      #Add the climate value to the observations dataframe
      observations$climate.datahstat[h] = climate.datahstat
    }
    
    #Center the window's climate values on year since we're interested in within-year variation
    observations$climate.datahstat.cent = wgdev(observations$climate.datahstat,observations$Year)
    
    #Now that each observation has a climate variable for the given window, get the linear and quadratic models
    if(func=="lin") {
      model = glmer(encountered ~ Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + (1|Community.Year/Focal.Group.Year),
                      data=observations,family=poisson) } else {
                        model = glmer(encountered ~ Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) + (1|Community.Year/Focal.Group.Year),
                                      data=observations,family=poisson) }
    
    return(list(model,observations))
    
}




#___________________________________________________________________________




####250-1 rainfall sum model#### 

#Test what rainfall period during the 250-1 days before the observation best predicts the number of 
#social groups encountered by the focal group during that observation. Early runs of 365-1 range indicated
#no relationship beyond 250 days

##Run sliding window analysis:
#EncM.sum.rain.250.1 = slidingwin_custom(observations=gbs,climate=rainall,climate.type="rain.sum",winopen=250,winclose=1,exclude=5)

#Save model file
#saveRDS(EncM.sum.rain.250.1,here::here("RDS model outputs","EncM_sum_rain_250_1.rds"))

#Read in model file
EncM.sum.rain.250.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_sum_rain_250_1.rds"))
EncM.sum.rain.250.1.lin = EncM.sum.rain.250.1[[1]]
EncM.sum.rain.250.1.quad = EncM.sum.rain.250.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.sum.rain.250.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.sum.rain.250.1.quad)
colnames(EncM.sum.rain.250.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.sum.rain.250.1.lin)

#Summary: Quad a bit better than linear, quad has multiple windows within 2 AIC units of top window. 
#Windows at 36-18 days, 45-3 days, 23-9 days, and 56-43 days. First three all overlap quite a bit


#Run single window of best windows of quad
EncM.sum.rain1 = singlewin_custom(observations=gbs,climate=rainall,climate.type="rain.sum",winopen=36,winclose=18,func="quad")
summary(EncM.sum.rain1[[1]])
EncM.sum.rain1.data = EncM.sum.rain1[[2]]
cor.test(EncM.sum.rain1.data$jdate,EncM.sum.rain1.data$climate.datahstat.cent) #-0.44

EncM.sum.rain2 = singlewin_custom(observations=gbs,climate=rainall,climate.type="rain.sum",winopen=45,winclose=4,func="quad")
summary(EncM.sum.rain2[[1]])
EncM.sum.rain2.data = EncM.sum.rain2[[2]]
cor.test(EncM.sum.rain2.data$jdate,EncM.sum.rain2.data$climate.datahstat.cent) #-0.61

EncM.sum.rain3 = singlewin_custom(observations=gbs,climate=rainall,climate.type="rain.sum",winopen=23,winclose=9,func="quad")
summary(EncM.sum.rain3[[1]])
EncM.sum.rain3.data = EncM.sum.rain3[[2]]
cor.test(EncM.sum.rain3.data$jdate,EncM.sum.rain3.data$climate.datahstat.cent) #-0.58

EncM.sum.rain4 = singlewin_custom(observations=gbs,climate=rainall,climate.type="rain.sum",winopen=56,winclose=43,func="quad")
summary(EncM.sum.rain4[[1]])
EncM.sum.rain4.data = EncM.sum.rain4[[2]]
cor.test(EncM.sum.rain4.data$jdate,EncM.sum.rain4.data$climate.datahstat.cent) #-0.05


#Result: All but 56-43 correlate strongly with jdate.



#56.43
rain.56.43 = glmer(encountered~Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) +
                     (1|Community.Year/Focal.Group.Year), data=EncM.sum.rain4.data, family=poisson)
summary(rain.56.43)

#Get a dataframe of potential values
rain.56.43.df = data.frame(expand.grid(climate.datahstat.cent=seq(-8, 15, by = 1),
                                       Focal.Sightings.cent=mean(EncM.sum.rain4.data$Focal.Sightings.cent),jdate=200,
                                       Year=as.factor(c(2016:2019))))

#Predict values 
rain.56.43.predict = predict(rain.56.43,newdata=rain.56.43.df,type="response",re.form=NA)

#Add predicted values
rain.56.43.df = cbind(rain.56.43.df,encountered = rain.56.43.predict)

#Plot
ggplot(data=EncM.sum.rain4.data,aes(x=climate.datahstat.cent,y=encountered)) +
  geom_point(data=EncM.sum.rain4.data,aes(x=climate.datahstat.cent,y=jitter(encountered),color=Year)) +
  geom_line(data=rain.56.43.df ,aes(x=climate.datahstat.cent,y=encountered,group=Year,color=Year),size=1) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_y_continuous(breaks=seq(0,7,by=1))





####Randomizations for rainfall sum

#Ran these randomizations using the cluster. In both quad and linear, looked like no important windows beyond 
#100 days, so ran randomizations at 100-1 days. 

#Script for 100-1 rainfall randomizations in: Cluster_rainfall_sum_250_1_rand.R
#Located in cluster scripts folder. 
#Code to run it in terminal is in the R script

#General idea for the climwin randwin() code is: code samples with replacement the dates from the biological 
#variable. So all of the dates get mixed up. These then get compared to the same weather data, which allows
#for temporal correlations in the weather to stay the same, such as a storm leading to high rainfall amounts
#for 2 days in a row. Because I'm using jdate as a predictor variable, I want the date to match the jdate, so 
#I wrote my own randomization code that keeps date and jdate the same but randomizes the order of encountered 
#and its fixed effects. Each encountered row stays with its associated fixed effects, order of rows just gets
#switched around. Only randomized within years. 

#Read in Delta-AIC values saved in randomization script
EncM.sum.rain.100.1.lin.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Rainfall_sum_100_1_lin_randomizations.csv"))
EncM.sum.rain.100.1.quad.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Rainfall_sum_100_1_quad_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.sum.rain.250.1.lin,EncM.sum.rain.100.1.lin.dataset.rand,"AIC",sample.size = 1183)
pvalue(EncM.sum.rain.250.1.quad,EncM.sum.rain.100.1.quad.dataset.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.sum.rain.100.1.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.sum.rain.250.1.lin$deltaAICc[[2]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.sum.rain.100.1.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.sum.rain.250.1.quad$deltaAICc[[2]],color="red",size=1.5) +
  ggtitle("Quadratic") + ylab("Frequency") + theme_cowplot()

#Results: Both linear and quadratic more than random.




#___________________________________________________________________________




####Mean minimum temperature sliding window####
#Minimum temperature is recorded at 9am at the Brisbane airport and assigned to the day it was measured. 
#Meaning it's  the minimum temperature from 9am the previous day to 9am the current day. Coldest time of day usually in 
#very early morning. 
#Hidalgo Aranzamendi paper on PCFWs found insect abundance closely correlated with temperature

#EncM.mean.mintemp.100.1 = slidingwin_custom(observations=gbs,climate=mintempall,climate.type="mintemp.mean",winopen=100,winclose=1,exclude=5)

#Save model file
#saveRDS(EncM.mean.mintemp.100.1,here::here("RDS model outputs","EncM_mean_mintemp_100_1.rds"))

#Read in model file
EncM.mean.mintemp.100.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_mintemp_100_1.rds"))
EncM.mean.mintemp.100.1.lin = EncM.mean.mintemp.100.1[[1]]
EncM.mean.mintemp.100.1.quad = EncM.mean.mintemp.100.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.mintemp.100.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.mintemp.100.1.lin)
colnames(EncM.mean.mintemp.100.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.mintemp.100.1.quad)

#Summary: Linear and quad show the same window - 42-15 days before observation - and linear and quad have similar top AIC values - top linear
#AIC score is within 2 units of top quad, so use linear.  

#Run single window of best windows
EncM.mean.mintemp1 = singlewin_custom(observations=gbs,climate=mintempall,climate.type="mintemp.mean",winopen=42,winclose=15,func="lin")
summary(EncM.mean.mintemp1[[1]])
EncM.mean.mintemp1.data = EncM.mean.mintemp1[[2]]
cor.test(EncM.mean.mintemp1.data$jdate,EncM.mean.mintemp1.data$climate.datahstat.cent) #-0.74

#Result: Mean mintemp is very correlated with jdate. 


###Randomizations for mean mintemp

#Ran these randomizations using the cluster. Script for 100-1 mintemp randomizations in:
#Cluster_mintemp_mean_100_0_rand_quad.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.mintemp.100.1.quad.rand = read.csv(here::here("Output files","Climwin randomizations","Mintemp_mean_100_1_quad_randomizations.csv"))
EncM.mean.mintemp.100.1.lin.rand = read.csv(here::here("Output files","Climwin randomizations","Mintemp_mean_100_1_lin_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.mintemp.100.1.quad,EncM.mean.mintemp.100.1.quad.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.mintemp.100.1.lin,EncM.mean.mintemp.100.1.lin.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.mintemp.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.mintemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.mintemp.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.mintemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: Mean min temp is more than random. 


#Also check importance of day-of minimum temperature
EncM.mean.mintemp.dayof = singlewin_custom(observations=gbs,climate=mintempall,climate.type="mintemp.mean",winopen=0,winclose=0,func="lin")
summary(EncM.mean.mintemp.dayof[[1]])
#Min temp day of is not significant and AIC nowhere near score of best window. Does not improve on baseline by more than 2 AIC values. 




#___________________________________________________________________________




####Detrended mean min temperature sliding window#####

#Get mintemp by year
mintempall.16 = mintempall[which(mintempall$Year=="2016"),]
mintempall.17 = mintempall[which(mintempall$Year=="2017"),]
mintempall.18 = mintempall[which(mintempall$Year=="2018"),]
mintempall.19 = mintempall[which(mintempall$Year=="2019"),]

#Look at loess fits for each year - lower span the more the line closely fits the data. Going for an 
#overall average trend here to compare residuals against. 0.5 looks good. 
mintempall.16$Min.Temp.dt = loess(Min.Temp ~ jDate, data=mintempall.16,span=0.5)$fitted
ggplot(data=mintempall.16, aes(x=Date,y=Min.Temp.dt)) + 
  geom_line(data=mintempall.16,aes(x=Date,y=Min.Temp)) + geom_line(color="red",size=1)

mintempall.17$Min.Temp.dt = loess(Min.Temp ~ jDate, data=mintempall.17,span=0.5)$fitted
ggplot(data=mintempall.17, aes(x=Date,y=Min.Temp.dt)) + 
  geom_line(data=mintempall.17,aes(x=Date,y=Min.Temp)) + geom_line(color="red",size=1)

mintempall.18$Min.Temp.dt = loess(Min.Temp ~ jDate, data=mintempall.18,span=0.5)$fitted
ggplot(data=mintempall.18, aes(x=Date,y=Min.Temp.dt)) + 
  geom_line(data=mintempall.18,aes(x=Date,y=Min.Temp)) + geom_line(color="red",size=1)

mintempall.19$Min.Temp.dt = loess(Min.Temp ~ jDate, data=mintempall.19,span=0.5)$fitted
ggplot(data=mintempall.19, aes(x=Date,y=Min.Temp.dt)) + 
  geom_line(data=mintempall.19,aes(x=Date,y=Min.Temp)) + geom_line(color="red",size=1)

#Get residuals - residuals will show how if a day is higher or lower than average (expected?)
mintempall.16$residuals = loess(Min.Temp ~ jDate, data=mintempall.16,span=0.5)$residuals
mintempall.17$residuals = loess(Min.Temp ~ jDate, data=mintempall.17,span=0.5)$residuals
mintempall.18$residuals = loess(Min.Temp ~ jDate, data=mintempall.18,span=0.5)$residuals
mintempall.19$residuals = loess(Min.Temp ~ jDate, data=mintempall.19,span=0.5)$residuals

#Combine all into one dataframe
mintempall.dt = rbind(mintempall.16,mintempall.17,mintempall.18,mintempall.19)


###Run sliding window
#EncM.mean.dtmintemp.100.1 = slidingwin_custom(observations=gbs,climate=mintempall.dt,climate.type="residuals",winopen=100,winclose=1,exclude=5)

#Save model file
#saveRDS(EncM.mean.dtmintemp.100.1,here::here("RDS model outputs","EncM_mean_dtmintemp_100_1.rds"))

#Read in model file
EncM.mean.dtmintemp.100.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_dtmintemp_100_1.rds"))
EncM.mean.dtmintemp.100.1.lin = EncM.mean.dtmintemp.100.1[[1]]
EncM.mean.dtmintemp.100.1.quad = EncM.mean.dtmintemp.100.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.dtmintemp.100.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.dtmintemp.100.1.quad)
colnames(EncM.mean.dtmintemp.100.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.dtmintemp.100.1.lin)


#Summary: Quad better than linear by 5 AIC values. Top windows all the same - around 51-36 for quad. 


#Run single window of best windows
EncM.mean.dtmintemp1 = singlewin_custom(observations=gbs,climate=mintempall.dt,climate.type="residuals",winopen=51,winclose=36,func="quad")
summary(EncM.mean.dtmintemp1[[1]])
EncM.mean.dtmintemp1.data = EncM.mean.dtmintemp1[[2]]
cor.test(EncM.mean.dtmintemp1.data$jdate,EncM.mean.dtmintemp1.data$climate.datahstat.cent) #0.21

#Result: Not very correlated with jdate


##Plot predictions from model
dtmintemp.51.36 = glmer(encountered~Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) + 
                        (1|Community.Year/Focal.Group.Year), data=EncM.mean.dtmintemp1.data, family=poisson)
summary(dtmintemp.51.36)

#Get a dataframe of potential values
dtmintemp.51.36.df = data.frame(expand.grid(climate.datahstat.cent=seq(-2.5, 3, by = 0.5),
                                          Focal.Sightings.cent=mean(EncM.mean.dtmintemp1.data$Focal.Sightings.cent),jdate=200,
                                          Year=as.factor(c(2016:2019))))

#Predict values 
dtmintemp.51.36.predict = predict(dtmintemp.51.36,newdata=dtmintemp.51.36.df,type="response",re.form=NA)

#Add predicted values
dtmintemp.51.36.df = cbind(dtmintemp.51.36.df,encountered = dtmintemp.51.36.predict)

#Plot
ggplot(data=EncM.mean.dtmintemp1.data,aes(x=climate.datahstat.cent,y=encountered)) +
  geom_point(data=EncM.mean.dtmintemp1.data,aes(x=climate.datahstat.cent,y=jitter(encountered),color=Year)) +
  geom_line(data=dtmintemp.51.36.df ,aes(x=climate.datahstat.cent,y=encountered,group=Year,color=Year),size=1) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_y_continuous(breaks=seq(0,7,by=1))



###Randomizations for detrended mean min temp

#Ran these randomizations using the cluster. Script for detrended mean min temp in:
#Cluster_dt_mintemp_mean_100_1_rand_custom.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.mintemp.residuals.100.1.quad.rand = read.csv(here::here("Output files","Climwin randomizations","Mintemp_dt_mean_100_1_quad_randomizations.csv"))
EncM.mean.mintemp.residuals.100.1.lin.rand = read.csv(here::here("Output files","Climwin randomizations","Mintemp_dt_mean_100_1_lin_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.dtmintemp.100.1.quad,EncM.mean.mintemp.residuals.100.1.quad.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.dtmintemp.100.1.lin,EncM.mean.mintemp.residuals.100.1.lin.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.mintemp.residuals.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmintemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.mintemp.residuals.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmintemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: Detrended mean min temp more than random. 


#Also check importance of day-of dt minimum temperature
EncM.mean.mintemp.residuals.dayof = singlewin_custom(observations=gbs,climate=mintempall.dt,climate.type="residuals",winopen=0,winclose=0,func="quad")
summary(EncM.mean.mintemp.dayof[[1]])
#Detrended min temp day of does not improve on baseline by more than 2 AIC values. 


#___________________________________________________________________________




####Mean maximum temperature sliding window####
#Maximum temperature is measured at 9am at the airport and assigned to the previous day. 
#So the temp value matches the date it's associated with. 

#EncM.mean.maxtemp.100.1 = slidingwin_custom(observations=gbs,climate=maxtempall,climate.type="maxtemp.mean",winopen=100,winclose=1,exclude=5)

#Save model file
#saveRDS(EncM.mean.maxtemp.100.1,here::here("RDS model outputs","EncM_mean_maxtemp_100_1.rds"))

#Read in model file
EncM.mean.maxtemp.100.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_maxtemp_100_1.rds"))
EncM.mean.maxtemp.100.1.lin = EncM.mean.maxtemp.100.1[[1]]
EncM.mean.maxtemp.100.1.quad = EncM.mean.maxtemp.100.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.maxtemp.100.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.maxtemp.100.1.quad)
colnames(EncM.mean.maxtemp.100.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.maxtemp.100.1.lin)


#summary: Top linear and quadratic windows within 2 AIC values of one another but show different windows. Linear 85-70 and quad 33-14.

#Run single window of best windows
EncM.mean.maxtemp1 = singlewin_custom(observations=gbs,climate=maxtempall,climate.type="maxtemp.mean",winopen=33,winclose=14,func="quad")
summary(EncM.mean.maxtemp1[[1]])
EncM.mean.maxtemp1.data = EncM.mean.maxtemp1[[2]]
cor.test(EncM.mean.maxtemp1.data$jdate,EncM.mean.maxtemp1.data$climate.datahstat.cent) #-0.30

EncM.mean.maxtemp2 = singlewin_custom(observations=gbs,climate=maxtempall,climate.type="maxtemp.mean",winopen=85,winclose=70,func="lin")
summary(EncM.mean.maxtemp2[[1]])
EncM.mean.maxtemp2.data = EncM.mean.maxtemp2[[2]]
cor.test(EncM.mean.maxtemp2.data$jdate,EncM.mean.maxtemp2.data$climate.datahstat.cent) #-0.89

#Result: 85-70 window is very correlated with jdate, 33-14 window is not as much. 



##Plot predictions from models
#33-14
maxtemp.33.14 = glmer(encountered~Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) + 
                          (1|Community.Year/Focal.Group.Year), data=EncM.mean.maxtemp1.data, family=poisson)
summary(maxtemp.33.14)

#Get a dataframe of potential values
maxtemp.33.14.df = data.frame(expand.grid(climate.datahstat.cent=seq(-2, 3, by = 0.5),
                                            Focal.Sightings.cent=mean(EncM.mean.maxtemp1.data$Focal.Sightings.cent),jdate=200,
                                            Year=as.factor(c(2016:2019))))

#Predict values 
maxtemp.33.14.predict = predict(maxtemp.33.14,newdata=maxtemp.33.14.df,type="response",re.form=NA)

#Add predicted values
maxtemp.33.14.df = cbind(maxtemp.33.14.df,encountered = maxtemp.33.14.predict)

#Plot
ggplot(data=EncM.mean.maxtemp1.data,aes(x=climate.datahstat.cent,y=encountered)) +
  geom_point(data=EncM.mean.maxtemp1.data,aes(x=climate.datahstat.cent,y=jitter(encountered),color=Year)) +
  geom_line(data=maxtemp.33.14.df ,aes(x=climate.datahstat.cent,y=encountered,group=Year,color=Year),size=1) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_y_continuous(breaks=seq(0,7,by=1))






###Randomizations for mean maxtemp

#Ran these randomizations using the cluster. Script for 100-1 maxtemp randomizations in:
#Cluster_maxtemp_mean_100_0_rand_quad.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.maxtemp.100.1.quad.rand = read.csv(here::here("Output files","Climwin randomizations","maxtemp_mean_100_1_quad_randomizations.csv"))
EncM.mean.maxtemp.100.1.lin.rand = read.csv(here::here("Output files","Climwin randomizations","maxtemp_mean_100_1_lin_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.maxtemp.100.1.quad,EncM.mean.maxtemp.100.1.quad.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.maxtemp.100.1.lin,EncM.mean.maxtemp.100.1.lin.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.maxtemp.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.maxtemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.maxtemp.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.maxtemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: Mean max temperature is different from random.



#___________________________________________________________________________




####Detrended mean max temperature sliding window#####

#Get maxtemp by year
maxtempall.16 = maxtempall[which(maxtempall$Year=="2016"),]
maxtempall.17 = maxtempall[which(maxtempall$Year=="2017"),]
maxtempall.18 = maxtempall[which(maxtempall$Year=="2018"),]
maxtempall.19 = maxtempall[which(maxtempall$Year=="2019"),]

#Look at loess fits for each year - lower span the more the line closely fits the data. Going for an 
#overall average trend here to compare residuals against. 0.5 looks good. 
maxtempall.16$Max.Temp.dt = loess(Max.Temp ~ jDate, data=maxtempall.16,span=0.5)$fitted
ggplot(data=maxtempall.16, aes(x=Date,y=Max.Temp.dt)) + 
  geom_line(data=maxtempall.16,aes(x=Date,y=Max.Temp)) + geom_line(color="red",size=1)

maxtempall.17$Max.Temp.dt = loess(Max.Temp ~ jDate, data=maxtempall.17,span=0.5)$fitted
ggplot(data=maxtempall.17, aes(x=Date,y=Max.Temp.dt)) + 
  geom_line(data=maxtempall.17,aes(x=Date,y=Max.Temp)) + geom_line(color="red",size=1)

maxtempall.18$Max.Temp.dt = loess(Max.Temp ~ jDate, data=maxtempall.18,span=0.5)$fitted
ggplot(data=maxtempall.18, aes(x=Date,y=Max.Temp.dt)) + 
  geom_line(data=maxtempall.18,aes(x=Date,y=Max.Temp)) + geom_line(color="red",size=1)

maxtempall.19$Max.Temp.dt = loess(Max.Temp ~ jDate, data=maxtempall.19,span=0.5)$fitted
ggplot(data=maxtempall.19, aes(x=Date,y=Max.Temp.dt)) + 
  geom_line(data=maxtempall.19,aes(x=Date,y=Max.Temp)) + geom_line(color="red",size=1)

#Get residuals - residuals will show how if a day is higher or lower than average (expected?)
maxtempall.16$residuals = loess(Max.Temp ~ jDate, data=maxtempall.16,span=0.5)$residuals
maxtempall.17$residuals = loess(Max.Temp ~ jDate, data=maxtempall.17,span=0.5)$residuals
maxtempall.18$residuals = loess(Max.Temp ~ jDate, data=maxtempall.18,span=0.5)$residuals
maxtempall.19$residuals = loess(Max.Temp ~ jDate, data=maxtempall.19,span=0.5)$residuals

#Combine all into one dataframe
maxtempall.dt = rbind(maxtempall.16,maxtempall.17,maxtempall.18,maxtempall.19)


###Run sliding window
#EncM.mean.dtmaxtemp.100.1 = slidingwin_custom(observations=gbs,climate=maxtempall.dt,climate.type="residuals",winopen=100,winclose=1,exclude=5)

#Save model file
#saveRDS(EncM.mean.dtmaxtemp.100.1,here::here("RDS model outputs","EncM_mean_dtmaxtemp_100_1.rds"))

#Read in model file
EncM.mean.dtmaxtemp.100.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_dtmaxtemp_100_1.rds"))
EncM.mean.dtmaxtemp.100.1.lin = EncM.mean.dtmaxtemp.100.1[[1]]
EncM.mean.dtmaxtemp.100.1.quad = EncM.mean.dtmaxtemp.100.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.dtmaxtemp.100.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.dtmaxtemp.100.1.quad)
colnames(EncM.mean.dtmaxtemp.100.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.dtmaxtemp.100.1.lin)

#Summary: 


#Run single window of best windows
EncM.mean.dtmaxtemp1 = singlewin_custom(observations=gbs,climate=maxtempall.dt,climate.type="residuals",winopen=85,winclose=60,func="quad")
summary(EncM.mean.dtmaxtemp1[[1]])
EncM.mean.dtmaxtemp1.data = EncM.mean.dtmaxtemp1[[2]]
cor.test(EncM.mean.dtmaxtemp1.data$jdate,EncM.mean.dtmaxtemp1.data$climate.datahstat.cent) #0.10

EncM.mean.dtmaxtemp2 = singlewin_custom(observations=gbs,climate=maxtempall.dt,climate.type="residuals",winopen=61,winclose=43,func="quad")
summary(EncM.mean.dtmaxtemp2[[1]])
EncM.mean.dtmaxtemp2.data = EncM.mean.dtmaxtemp2[[2]]
cor.test(EncM.mean.dtmaxtemp2.data$jdate,EncM.mean.dtmaxtemp2.data$climate.datahstat.cent) #-0.43


#Result: 85-60 not correlated with jdate, 61-43 is correlated with jdate. 



##Plot predictions from models
#85-60
dtmaxtemp.85.60 = glmer(encountered~Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) + 
                        (1|Community.Year/Focal.Group.Year), data=EncM.mean.dtmaxtemp1.data, family=poisson)
summary(dtmaxtemp.85.60)

#Get a dataframe of potential values
dtmaxtemp.85.60.df = data.frame(expand.grid(climate.datahstat.cent=seq(-1, 1, by = 0.1),
                                          Focal.Sightings.cent=mean(EncM.mean.dtmaxtemp1.data$Focal.Sightings.cent),jdate=200,
                                          Year=as.factor(c(2016:2019))))

#Predict values 
dtmaxtemp.85.60.predict = predict(dtmaxtemp.85.60,newdata=dtmaxtemp.85.60.df,type="response",re.form=NA)

#Add predicted values
dtmaxtemp.85.60.df = cbind(dtmaxtemp.85.60.df,encountered = dtmaxtemp.85.60.predict)

#Plot
ggplot(data=EncM.mean.dtmaxtemp1.data,aes(x=climate.datahstat.cent,y=encountered)) +
  geom_point(data=EncM.mean.dtmaxtemp1.data,aes(x=climate.datahstat.cent,y=jitter(encountered),color=Year)) +
  geom_line(data=dtmaxtemp.85.60.df ,aes(x=climate.datahstat.cent,y=encountered,group=Year,color=Year),size=1) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_y_continuous(breaks=seq(0,7,by=1))






###Randomizations for detrended mean max temp

#Ran these randomizations using the cluster. Script for detrended mean max temp in:
#Cluster_detrended_maxtemp_mean_100_0_rand_quad.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.maxtemp.residuals.100.1.quad.rand = read.csv(here::here("Output files","Climwin randomizations","Maxtemp_dt_mean_100_1_quad_randomizations.csv"))
EncM.mean.maxtemp.residuals.100.1.lin.rand = read.csv(here::here("Output files","Climwin randomizations","Maxtemp_dt_mean_100_1_lin_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.dtmaxtemp.100.1.quad,EncM.mean.maxtemp.residuals.100.1.quad.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.dtmaxtemp.100.1.lin,EncM.mean.maxtemp.residuals.100.1.lin.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.maxtemp.residuals.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmaxtemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.maxtemp.residuals.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmaxtemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: Detrended mean max temp more than random. 



#___________________________________________________________________________




####Windspeed sliding window####
#Prediction for windspeed is that sometimes when it was very windy, the birds tended to congregate in the 
#forest where there was shelter. Thus looking at window longer than 10 days before the observation 
#would not make much sense. Too far out and windspeed would correlate strongly with rainfall. 


#EncM.mean.windspeed.10.0 = slidingwin_custom(observations=gbs,climate=windall,climate.type="windspeed.mean",winopen=10,winclose=0,exclude=0)

#Save model file
#saveRDS(EncM.mean.windspeed.10.0,here::here("RDS model outputs","EncM_mean_windspeed_10_0.rds"))

#Read in model file
EncM.mean.windspeed.10.0 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_windspeed_10_0.rds"))
EncM.mean.windspeed.10.0.lin = EncM.mean.windspeed.10.0[[1]]
EncM.mean.windspeed.10.0.quad = EncM.mean.windspeed.10.0[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.windspeed.10.0.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.windspeed.10.0.quad)
colnames(EncM.mean.windspeed.10.0.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.windspeed.10.0.lin)

#summary: Lin is slightly better than quad but both delta-AIC values are very low. Likely not different from random



###Randomizations for windspeed sliding window

#Ran these randomizations using the cluster. Script in: Cluster_windspeed_mean_10_0_rand_custom.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.wind.10.0.lin.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Windspeed_mean_10_0_lin_randomizations.csv"))
EncM.mean.wind.10.0.quad.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Windspeed_mean_10_0_quad_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.windspeed.10.0.lin,EncM.mean.wind.10.0.lin.dataset.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.windspeed.10.0.quad,EncM.mean.wind.10.0.quad.dataset.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.wind.10.0.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.windspeed.10.0.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.wind.10.0.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.windspeed.10.0.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: windspeed is not different from random



#___________________________________________________________________________




####NDVI Sliding window####

#Run NDVI separately from rainfall? They're obviously connected. 
#NDVI imported at beginning of script is all years, yearly polygons (specific to each season)
#Split into years to interpolate data - huge gap before 1st data point leading up to 2016 season, don't want to
#interpolate the huge gaps

ndvi16 = ndviall[which(ndviall$Year=="2016"),]
ndvi17 = ndviall[which(ndviall$Year=="2017"),]
ndvi18 = ndviall[which(ndviall$Year=="2018"),]
ndvi19 = ndviall[which(ndviall$Year=="2019"),]

##Use linear interpolation to get more values along the NDVI line
#Interpolation calculates the values between known values
#Will hopefully give a better estimate of the mean values of each half than just a few known data points will
#2016
ndvi16$Date[1] - ndvi16$Date[nrow(ndvi16)] #Get sample size
ndvi16.i = data.frame(approx(x=ndvi16$Date,y=ndvi16$NDVI,n=241)) #Use sample size +1 to get a value for each date
ndvi16.i$Date = as.Date(ndvi16.i$x,origin="1970-01-01") #This is the default origin which is what the spline uses
colnames(ndvi16.i)[1:2]=c("Date.n","NDVI")
#2017
ndvi17$Date[1] - ndvi17$Date[nrow(ndvi17)] #Get sample size
ndvi17.i = data.frame(approx(x=ndvi17$Date,y=ndvi17$NDVI,n=353)) #Use sample size +1 to get a value for each date
ndvi17.i$Date = as.Date(ndvi17.i$x,origin="1970-01-01") #This is the default origin which is what the spline uses
colnames(ndvi17.i)[1:2]=c("Date.n","NDVI")
#2018
ndvi18$Date[1] - ndvi18$Date[nrow(ndvi18)] #Get sample size
ndvi18.i = data.frame(approx(x=ndvi18$Date,y=ndvi18$NDVI,n=305)) #Use sample size +1 to get a value for each date
ndvi18.i$Date = as.Date(ndvi18.i$x,origin="1970-01-01") #This is the default origin which is what the spline uses
colnames(ndvi18.i)[1:2]=c("Date.n","NDVI")
#2019
ndvi19$Date[1] - ndvi19$Date[nrow(ndvi19)] #Get sample size
ndvi19.i = data.frame(approx(x=ndvi19$Date,y=ndvi19$NDVI,n=209)) #Use sample size +1 to get a value for each date
ndvi19.i$Date = as.Date(ndvi19.i$x,origin="1970-01-01") #This is the default origin which is what the spline uses
colnames(ndvi19.i)[1:2]=c("Date.n","NDVI")

#Combine interpolated dataframes
ndviall.i = rbind(ndvi16.i,ndvi17.i,ndvi18.i,ndvi19.i)

#Figure out the range for the sliding window
#2016 has the latest start date - 2015-04-17
gbs$Date[1] - ndviall.i$Date[1] #63 days is max range. 


#63 days is max range for 2016 season, so run model at 63-0 days. 


###Run sliding window
#EncM.mean.ndvi.63.1 = slidingwin_custom(observations=gbs,climate=ndviall.i,climate.type="NDVI.mean",winopen=63,winclose=1,exclude=0)

#Save model file
#saveRDS(EncM.mean.ndvi.63.1,here::here("RDS model outputs","EncM_mean_ndvi_63_1.rds"))

#Read in model file
EncM.mean.ndvi.63.1 = readRDS(here::here("Output files","Climwin model outputs","EncM_mean_ndvi_63_1.rds"))
EncM.mean.ndvi.63.1.lin = EncM.mean.ndvi.63.1[[1]]
EncM.mean.ndvi.63.1.quad = EncM.mean.ndvi.63.1[[2]]

#Plot delta plot - need to change label to deltaAICc to work with climwin code
colnames(EncM.mean.ndvi.63.1.quad)[5] <- "deltaAICc"
plotdelta(dataset=EncM.mean.ndvi.63.1.quad)
colnames(EncM.mean.ndvi.63.1.lin)[4] <- "deltaAICc"
plotdelta(dataset=EncM.mean.ndvi.63.1.lin)

#summary: 

#Run single window of best windows
EncM.mean.ndvi1 = singlewin_custom(observations=gbs,climate=ndviall.i,climate.type="NDVI.mean",winopen=63,winclose=51,func="lin")
summary(EncM.mean.ndvi1[[1]])
EncM.mean.ndvi1.data = EncM.mean.ndvi1[[2]]
cor.test(EncM.mean.ndvi1.data$jdate,EncM.mean.ndvi1.data$climate.datahstat.cent)
plot(EncM.mean.ndvi1.data$jdate,EncM.mean.ndvi1.data$climate.datahstat.cent)

#Result: NDVI window strongly correlated with jdate

mean.ndvi.63.51 = glmer(encountered~Focal.Sightings.cent + jdate + Year + climate.datahstat.cent + I(climate.datahstat.cent^2) +
                          (1|Community.Year/Focal.Group.Year), data=EncM.mean.ndvi1.data, family=poisson)
summary(mean.ndvi.63.51)

#Get a dataframe of potential values
mean.ndvi.63.51.df = data.frame(expand.grid(climate.datahstat.cent=seq(-0.05, 0.04, by = 0.01),
                                            Focal.Sightings.cent=mean(EncM.mean.dtmaxtemp2.data$Focal.Sightings.cent),jdate=200,
                                            Year=as.factor(c(2016:2019))))

#Predict values 
mean.ndvi.63.51.predict = predict(mean.ndvi.63.51,newdata=mean.ndvi.63.51.df,type="response",re.form=NA)

#Add predicted values
mean.ndvi.63.51.df = cbind(mean.ndvi.63.51.df,encountered = mean.ndvi.63.51.predict)

#Plot
ggplot(data=EncM.mean.ndvi1.data,aes(x=climate.datahstat.cent,y=encountered)) +
  geom_point(data=EncM.mean.ndvi1.data,aes(x=climate.datahstat.cent,y=jitter(encountered),color=Year)) +
  geom_line(data=mean.ndvi.63.51.df ,aes(x=climate.datahstat.cent,y=encountered,group=Year,color=Year),size=1) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  scale_y_continuous(breaks=seq(0,7,by=1)) 


###Randomizations for NDVI 


#Ran these randomizations using the cluster. Script in: Cluster_ndvi_mean_63_1_rand_custom.R

#Read in Delta-AIC values saved in randomization script
EncM.mean.ndvi.63.1.lin.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","NDVI_mean_63_1_lin_randomizations.csv"))
EncM.mean.ndvi.63.1.quad.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","NDVI_mean_63_1_quad_randomizations.csv"))

#Calculate pvalues - sample size only important for 
pvalue(EncM.mean.ndvi.63.1.lin,EncM.mean.ndvi.63.1.lin.dataset.rand,"AIC",sample.size = 1183)
pvalue(EncM.mean.ndvi.63.1.quad,EncM.mean.ndvi.63.1.quad.dataset.rand,"AIC",sample.size = 1183)

#Histograms
ggplot(data=EncM.mean.ndvi.63.1.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.ndvi.63.1.lin$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()
ggplot(data=EncM.mean.ndvi.63.1.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.ndvi.63.1.quad$deltaAICc[[1]],color="red",size=1.5) +
  ylab("Frequency") + theme_cowplot()

#Result: Quad NDVI window more than random but lin window is not. 




#_________#_________#_________#_________#_________#_________#_________#_________#




####Final Models####


###Add everything into one dataframe
#These are the variables who's correlation with jdate was less than 0.4
#Other variables are probably just picking up the jdate trend. 
gbs$amountcm56.43 = EncM.sum.rain4.data$climate.datahstat.cent
gbs$dtmintemp51.36 = EncM.mean.dtmintemp1.data$climate.datahstat.cent
gbs$maxtemp33.14 = EncM.mean.maxtemp1.data$climate.datahstat.cent
gbs$dtmaxtemp85.60 = EncM.mean.dtmaxtemp1.data$climate.datahstat.cent


###Look at correlations among all the variables
library("PerformanceAnalytics")
gbs.sub = gbs[,c(13,15:19)] #Get fixed effect columns
#Standardize all fixed effects
for (i in 1:ncol(gbs.sub)) {
  gbs.sub[i] = scale(gbs.sub[i])
}
#Plot correlations
chart.Correlation(gbs.sub, histogram=TRUE, pch=19)



###Only use one of mintemp or maxtemp for mean and mean detrended. Mean and mintemp should be representing about the same data. 
#Doesn't make sense to include all of them in a single model. 
#Use version of each with the lowest AIC value for the best window. 
#mean detrended mintemp delta AIC = -20
#mean detrended maxtemp delta AIC = -28 - use detrended maxtemp



#Pre-final model
EncM.prefinal = glmer(encountered~scale(Focal.Sightings.cent) + scale(jdate) + Year + poly(scale(amountcm56.43),2,raw=T) +
                   poly(scale(dtmaxtemp85.60),2,raw=T) + poly(scale(maxtemp33.14),2,raw=T) +
                  (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
summary(EncM.prefinal)


#Check covariance
car::vif(EncM.prefinal)

#VIF of maxtemp33.14 is greater than 3, so following Zuur et al. (2010?) remove from model and re-run
#There is some arguments to leave in correlated variables, see Freckleton/Maclean papers but in this case
#it's dt maxtemp correlated with maxtemp. If it were temp and rainfall that would be different. 

EncM.prefinal2 = glmer(encountered~scale(Focal.Sightings.cent) + scale(jdate) + Year + poly(scale(amountcm56.43),2,raw=T) +
                        poly(scale(dtmaxtemp85.60),2,raw=T) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
summary(EncM.prefinal2)
car::vif(EncM.prefinal2) #Low VIFs now. 



###Make sure all of the variables in the model still belong
#Check for hitch-hiker variables - see how much AIC changes when variable is removed

#Remove focal sightings
EncM.no.focal.sightings.cent = glmer(encountered~scale(jdate) + Year + poly(scale(amountcm56.43),2,raw=T) +
                                       poly(scale(dtmaxtemp85.60),2,raw=T) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
AIC(EncM.no.focal.sightings.cent) - AIC(EncM.prefinal) #AIC increases by 41

#Remove jdate
EncM.no.jdate = glmer(encountered~scale(Focal.Sightings.cent) + Year + poly(scale(amountcm56.43),2,raw=T) +
                        poly(scale(dtmaxtemp85.60),2,raw=T) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
AIC(EncM.no.jdate) - AIC(EncM.prefinal) #AIC increases by 37

#Remove amountcm56.43
EncM.no.amountcm56.43 = glmer(encountered~scale(Focal.Sightings.cent) + scale(jdate) + Year + 
                                poly(scale(dtmaxtemp85.60),2,raw=T) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
AIC(EncM.no.amountcm56.43) - AIC(EncM.prefinal) #AIC increases by 5

#Remove dtmaxtemp85.60
EncM.no.dtmaxtemp85.60 = glmer(encountered~scale(Focal.Sightings.cent) + scale(jdate) + Year + poly(scale(amountcm56.43),2,raw=T) +
                                (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
AIC(EncM.no.dtmaxtemp85.60) - AIC(EncM.prefinal) #AIC increases by 10


###Final Model - same as prefinal2###
EncM.final = EncM.prefinal2
summary(EncM.final)


#Overdispersion test
#Function to calculate a point estimate of overdispersion from a mixed model object
od.point<-function(modelobject){
  x<-sum(resid(modelobject,type="pearson")^2)
  rdf<-summary(modelobject)$AICtab[5]
  return(x/rdf)
}
od.point(EncM.final)

#Residuals test
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = EncM.final, n = 250)
plot(simulationOutput)
testResiduals(simulationOutput)

#Residuals look good




####Calculate R2 values for final model 
#Either MuMin or piecewiseSEM packages work well here - they give the same results
#peformance package generates different results, not sure what method it's using to calculate observation-level variance. 
#Nakagawa, S., Johnson, P.C.D., Schielzeth, H. (2017) says that for a log-link function (like a poisson)
#the trigamma method is the most accurate. piecewiseSEM package will calculate R2 using that method

library(piecewiseSEM)
rsquared(EncM.final) #Marginal is variance explained by the fixed effects, conditional is variance explained by the entire model.
rsquared(EncM.final)$Marginal *100

#R2 value is quite low for the model 0.06 - but given the plots that's not super surprising. Looks like weather does have an effect
#but there are many other factors likely influencing variation in interactions among social groups. 

#Can get % delta R2 for each fixed effect by getting the change in the R2 value when the fixed effect is removed from the model.
#Lv et al 2019 does this, following Froy, Walling, Pemberton, Cluttonâ€�Brock, & Kruuk, 2016.

#Focal sightings R2
-(rsquared(EncM.final)$Marginal - rsquared(EncM.no.focal.sightings.cent)$Marginal) *100

#jdate R2
-(rsquared(EncM.final)$Marginal - rsquared(EncM.no.jdate)$Marginal) *100

#amountcm56.43 R2
-(rsquared(EncM.final)$Marginal - rsquared(EncM.no.amountcm56.43)$Marginal) *100

#dtmaxtemp85.60 R2
-(rsquared(EncM.final)$Marginal - rsquared(EncM.no.dtmaxtemp85.60)$Marginal) *100

#Season R2
EncM.no.year = glmer(encountered~scale(Focal.Sightings.cent) + scale(jdate) + poly(scale(amountcm56.43),2,raw=T) +
                         poly(scale(dtmaxtemp85.60),2,raw=T) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
-(rsquared(EncM.final)$Marginal - rsquared(EncM.no.year)$Marginal) *100







####Plot model predictions####

#Get prediction model format
EncM.final.predict = glmer(encountered~Focal.Sightings.cent + jdate + Year + amountcm56.43 + I(amountcm56.43^2) +
                             dtmaxtemp85.60 + I(dtmaxtemp85.60^2) + (1|Community.Year/Focal.Group.Year), data=gbs, family=poisson)
summary(EncM.final.predict) #same AIC as final model

#Function to get confidence intervals from Ben Bolker
easyPredCI <- function(model,newdata=NULL,alpha=0.05) {
  ## baseline prediction, on the linear predictor (logit) scale:
  pred0 <- predict(model,re.form=NA,newdata=newdata)
  ## fixed-effects model matrix for new data
  X <- model.matrix(formula(model,fixed.only=TRUE)[-2],newdata)
  beta <- fixef(model) ## fixed-effects coefficients
  V <- vcov(model)     ## variance-covariance matrix of beta
  pred.se <- sqrt(diag(X %*% V %*% t(X))) ## std errors of predictions
  ## inverse-link function
  linkinv <- family(model)$linkinv
  ## construct 95% Normal CIs on the link scale and
  ##  transform back to the response (probability) scale:
  crit <- -qnorm(alpha/2)
  linkinv(cbind(conf.low=pred0-crit*pred.se,
                conf.high=pred0+crit*pred.se))
}



###Plot Focal Sightings 

#Get prediction dataframe - need to have columns for each of the fixed effects here
EncM.final.FS.df = data.frame(expand.grid(amountcm56.43=mean(gbs$amountcm56.43),
                                            Focal.Sightings.cent=seq(-9,8,by=1),jdate=200,
                                            dtmaxtemp85.60=mean(gbs$dtmaxtemp85.60),
                                            "I(dtmaxtemp85.60^2)"=mean(gbs$dtmaxtemp85.60)^2,
                                            Year=as.factor(c(2016:2019))))

#Predict values 
EncM.final.FS = predict(EncM.final.predict, newdata=EncM.final.FS.df, type="response", re.form=NA)

#Get confidence intervals
EncM.final.FS.CI <- easyPredCI(EncM.final.predict, newdata=EncM.final.FS.df)

#Combine prediction dataframes
EncM.final.FS.df2 = cbind(EncM.final.FS.df, encountered=EncM.final.FS, EncM.final.FS.CI)

#Plot
FS.predict.plot = ggplot(data=gbs,aes(x=Focal.Sightings.cent,y=encountered)) +
  geom_point(data=gbs,aes(x=Focal.Sightings.cent,y=jitter(encountered),color=Year)) +
  geom_ribbon(data=EncM.final.FS.df2,
              aes(x=Focal.Sightings.cent,y=encountered,ymin=conf.low,ymax=conf.high,group=as.factor(Year)),alpha=0.2) +
  geom_line(data=EncM.final.FS.df2 ,aes(x=Focal.Sightings.cent,y=encountered,
                                          group=as.factor(Year),color=as.factor(Year)),size=1.3) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme_cowplot() +
  scale_y_continuous(breaks=seq(0,7,by=1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Sightings in Observation") + ylab("Social Groups Encountered") + labs(color = "Season") +
  guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) #+ facet_wrap(vars(Year))
print(FS.predict.plot)





###Plot jdate
EncM.final.jdate.df = data.frame(expand.grid(amountcm56.43=mean(gbs$amountcm56.43),
                                          Focal.Sightings.cent=mean(gbs$Focal.Sightings.cent),jdate=seq(170,240,by=1),
                                          dtmaxtemp85.60=mean(gbs$dtmaxtemp85.60),
                                          "I(dtmaxtemp85.60^2)"=mean(gbs$dtmaxtemp85.60)^2,
                                          Year=as.factor(c(2016:2019))))

#Predict values 
EncM.final.jdate = predict(EncM.final.predict, newdata=EncM.final.jdate.df, type="response", re.form=NA)

#Get confidence intervals
EncM.final.jdate.CI <- easyPredCI(EncM.final.predict, newdata=EncM.final.jdate.df)

#Combine prediction dataframes
EncM.final.jdate.df2 = cbind(EncM.final.jdate.df, encountered=EncM.final.jdate, EncM.final.jdate.CI)

#Plot
jdate.predict.plot = ggplot(data=gbs,aes(x=jdate,y=encountered)) +
  geom_point(data=gbs,aes(x=jdate,y=jitter(encountered),color=Year)) +
  geom_ribbon(data=EncM.final.jdate.df2,
              aes(x=jdate,y=encountered,ymin=conf.low,ymax=conf.high,group=as.factor(Year)),alpha=0.2) +
  geom_line(data=EncM.final.jdate.df2 ,aes(x=jdate,y=encountered,
                                        group=as.factor(Year),color=as.factor(Year)),size=1.3) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme_cowplot() +
  scale_y_continuous(breaks=seq(0,7,by=1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Day of Year") + ylab("Social Groups Encountered") + labs(color = "Season") +
  guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) #+ facet_wrap(vars(Year))
print(jdate.predict.plot)






###Plot rainfall 

#Get prediction dataframe - need to have columns for each of the fixed effects here
#Use mean focal sightings per observation and middle jdate

##Plot by year
EncM.final.rain.df = data.frame(expand.grid(amountcm56.43=seq(-8, 15, by = 1),
                                            Focal.Sightings.cent=mean(gbs$Focal.Sightings.cent),jdate=200,
                                            dtmaxtemp85.60=mean(gbs$dtmaxtemp85.60),
                                            "I(dtmaxtemp85.60^2)"=mean(gbs$dtmaxtemp85.60)^2,
                                            Year=as.factor(c(2016:2019))))

#Predict values 
EncM.final.rain = predict(EncM.final.predict, newdata=EncM.final.rain.df, type="response", re.form=NA)

#Get confidence intervals
EncM.final.rain.CI <- easyPredCI(EncM.final.predict, newdata=EncM.final.rain.df)

#Combine prediction dataframes
EncM.final.rain.df2 = cbind(EncM.final.rain.df, encountered=EncM.final.rain, EncM.final.rain.CI)

#Plot
rain.predict.plot = ggplot(data=gbs,aes(x=amountcm56.43,y=encountered)) +
  geom_point(data=gbs,aes(x=amountcm56.43,y=jitter(encountered),color=Year)) +
  geom_ribbon(data=EncM.final.rain.df2,
              aes(x=amountcm56.43,y=encountered,ymin=conf.low,ymax=conf.high,group=as.factor(Year)),alpha=0.2) +
  geom_line(data=EncM.final.rain.df2 ,aes(x=amountcm56.43,y=encountered,
                                          group=as.factor(Year),color=as.factor(Year)),size=1.3) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme_cowplot() +
  scale_y_continuous(breaks=seq(0,7,by=1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab("Rainfall 56-43 days before observation (cm)") + ylab("Social Groups Encountered") + labs(color = "Season") +
  guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) #+ facet_wrap(vars(Year))
print(rain.predict.plot)






###Plot dt maxtemp 

##Plot by year
EncM.final.dtmaxtemp.df = data.frame(expand.grid(amountcm56.43=mean(gbs$amountcm56.43),
                                            Focal.Sightings.cent=mean(gbs$Focal.Sightings.cent),jdate=200,
                                            dtmaxtemp85.60=seq(-1.2,1,by=0.1),
                                            "I(dtmaxtemp85.60^2)"=seq(-1.2,1,by=0.1)^2,
                                            Year=as.factor(c(2016:2019))))

#Predict values 
EncM.final.dtmaxtemp = predict(EncM.final.predict, newdata=EncM.final.dtmaxtemp.df, type="response", re.form=NA)

#Get confidence intervals
EncM.final.dtmaxtemp.CI <- easyPredCI(EncM.final.predict, newdata=EncM.final.dtmaxtemp.df)

#Combine prediction dataframes
EncM.final.dtmaxtemp.df2 = cbind(EncM.final.dtmaxtemp.df, encountered=EncM.final.dtmaxtemp, EncM.final.dtmaxtemp.CI)

#Plot
dtmaxtemp.predict.plot = ggplot(data=gbs,aes(x=dtmaxtemp85.60,y=encountered)) +
  geom_point(data=gbs,aes(x=dtmaxtemp85.60,y=jitter(encountered),color=Year)) +
  geom_ribbon(data=EncM.final.dtmaxtemp.df2,
              aes(x=dtmaxtemp85.60,y=encountered,ymin=conf.low,ymax=conf.high,group=as.factor(Year)),alpha=0.2) +
  geom_line(data=EncM.final.dtmaxtemp.df2 ,aes(x=dtmaxtemp85.60,y=encountered,
                                          group=as.factor(Year),color=as.factor(Year)),size=1.3) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme_cowplot() +
  scale_y_continuous(breaks=seq(0,7,by=1)) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  xlab(expression(atop("Mean detrended maximum temperature",paste("85-60 days before observation (Celcius)")))) + 
  ylab("Social Groups Encountered") + labs(color = "Season") +
  guides(colour = guide_legend(override.aes = list(size=5,linetype=0))) #+ facet_wrap(vars(Year))
print(dtmaxtemp.predict.plot)











####Plot climate data#### 

library(ggpubr)

#Get times when we were there for graphs
first.2016 = gbs[which(gbs$Year=="2016"),]$Date[1]
last.2016 = gbs[which(gbs$Year=="2016"),]$Date[nrow(gbs[which(gbs$Year=="2016"),])]
first.2017 = gbs[which(gbs$Year=="2017"),]$Date[1]
last.2017 = gbs[which(gbs$Year=="2017"),]$Date[nrow(gbs[which(gbs$Year=="2017"),])]
first.2018 = gbs[which(gbs$Year=="2018"),]$Date[1]
last.2018 = gbs[which(gbs$Year=="2018"),]$Date[nrow(gbs[which(gbs$Year=="2018"),])]
first.2019 = gbs[which(gbs$Year=="2019"),]$Date[1]
last.2019 = gbs[which(gbs$Year=="2019"),]$Date[nrow(gbs[which(gbs$Year=="2019"),])]

years.df = data.frame(years=c("2016","2017","2018","2019"),
                      first=c(first.2016,first.2017,first.2018,first.2019),
                      last=c(last.2016,last.2017,last.2018,last.2019))


##Rainfall
A = ggplot() + 
  geom_line(data=rainall[which(rainall$Year!="2015"),], aes(x=Date,y=amountcm, color=Year,), size=1) + 
  scale_x_date(date_labels = "%Y", limits=as.Date(c("2015-01-01","2018-12-31"))) +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  xlab("Date") + ylab("Rainfall (cm)") + ylim(0,15) + theme_cowplot() + theme(legend.position = "none") +
  geom_errorbarh(data=years.df, aes(xmax=last,xmin=first,y=14,height=2),size=1,color="black") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -20, b = 0, l = 0))) + ggtitle("A")

#For NDVI assign across-year colors to upcoming year - get last date of 2016, first date of 2017, etc. 
across.years = rbind(ndviall[ndviall$Year=="2016",][nrow(ndviall[ndviall$Year=="2016",]),],
                     ndviall[ndviall$Year=="2017",][1,],
                     ndviall[ndviall$Year=="2017",][nrow(ndviall[ndviall$Year=="2017",]),],
                     ndviall[ndviall$Year=="2018",][1,],
                     ndviall[ndviall$Year=="2018",][nrow(ndviall[ndviall$Year=="2018",]),],
                     ndviall[ndviall$Year=="2019",][1,])
across.years$Year = c("2017","2017","2018","2018","2019","2019")
across.years$group = "group"

B = ggplot() +
  geom_line(data=ndviall,aes(x=Date,y=NDVI,color=Year),size=2) + ylim(0.4,0.8) + 
  scale_x_date(date_labels = "%Y",limits=as.Date(c("2015-01-01","2018-12-31"))) + theme_cowplot() +
  geom_line(data=across.years,aes(x=Date,y=NDVI,color=Year),size=2) +
  geom_point(data=ndviall,aes(x=Date,y=NDVI),size=2) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -20, b = 0, l = 0))) + ggtitle("B")

C = ggplot() +
  scale_x_date(date_labels = "%Y",limits=as.Date(c("2015-01-01","2018-12-31"))) + theme_cowplot() + 
  geom_line(data=mintempall[which(mintempall$Year!="2015"),],aes(x=Date,y=Min.Temp,color=Year),size=0.7) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  ylab((expression(atop(NA,atop(textstyle("Minimum"), paste(textstyle("Temperature (C)"))))))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
  geom_line(data=mintempall.dt,aes(x=Date,y=Min.Temp.dt),size=0.8) + ggtitle("C")

ggarrange(A,B,C,ncol=1,align="v")

D = ggplot() + geom_point(data=gbs,aes(x=Date,y=encountered)) + ylim(0,8) +
  scale_x_date(date_labels = "%Y",limits=as.Date(c("2015-01-01","2018-12-31"))) + theme_cowplot() + 
  geom_smooth(data=gbs,aes(x=Date,y=encountered,color=Year),size=1.5) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  ylab("Social groups encountered")

ggarrange(A,B,C,D,ncol=1,align="v")

E = ggplot() +
  scale_x_date(date_labels = "%Y",limits=as.Date(c("2015-01-01","2018-12-31"))) + theme_cowplot() + 
  geom_line(data=maxtempall[which(maxtempall$Year!="2015"),],aes(x=Date,y=Max.Temp,color=Year),size=1) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  ylab((expression(atop(NA,atop(textstyle("Maximum"), paste(textstyle("Temperature (C)"))))))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
  geom_line(data=maxtempall.dt,aes(x=Date,y=Max.Temp.dt),size=0.8) + ggtitle("C")

ggarrange(A,B,E,ncol=1,align="v")





###Plot with months

G = ggplot() + 
  geom_line(data=rainall[which(rainall$Year!="2015"),], aes(x=Date,y=amountcm, color=Year,), size=0.8) + 
  scale_x_date(date_labels = "%b", limits=as.Date(c("2015-01-01","2018-12-31")),date_breaks="3 months") +
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) +
  xlab("Date") + ylab("Rainfall (cm)") + ylim(0,15) + theme_cowplot() + theme(legend.position = "none") +
  geom_errorbarh(data=years.df, aes(xmax=last,xmin=first,y=14,height=2),size=1,color="black") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -20, b = 0, l = 0))) + ggtitle("A") +
  theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=12))

H = ggplot() +
  geom_line(data=ndviall,aes(x=Date,y=NDVI,color=Year),size=2) + ylim(0.4,0.8) + 
  scale_x_date(date_labels = "%b", limits=as.Date(c("2015-01-01","2018-12-31")),date_breaks="3 months") + theme_cowplot() +
  geom_line(data=across.years,aes(x=Date,y=NDVI,color=Year),size=2) +
  geom_point(data=ndviall,aes(x=Date,y=NDVI),size=1.5) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -20, b = 0, l = 0))) + ggtitle("B") +
  theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=12))

I = ggplot() +
  scale_x_date(date_labels = "%b", limits=as.Date(c("2015-01-01","2018-12-31")),date_breaks="3 months") + theme_cowplot() + 
  geom_line(data=maxtempall[which(maxtempall$Year!="2015"),],aes(x=Date,y=Max.Temp,color=Year),size=0.8) + 
  scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + theme(legend.position = "none") +
  ylab((expression(atop(NA,atop(textstyle("Maximum"), paste(textstyle("Temperature (C)"))))))) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0))) +
  geom_line(data=maxtempall.dt,aes(x=Date,y=Max.Temp.dt),size=0.8) + ggtitle("C") + 
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=12))

ggarrange(G,H,I,ncol=1,align="v")






####What does mean detrended maxtemp 85-60 mean?####

#Try plotting maxtemp by julian date
ggplot(data=gbs, aes(x=jdate,y=dtmaxtemp85.60)) + geom_point(aes(color=Year),size=2) + 
  geom_line(aes(group=Year,color=Year),size=1.3) + 
  theme_cowplot() + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06"))
#Clearly some relationships within years, not super correlated with jdate but some trends

#How does mean detrended maxtemp vary with rainfall? 
rain.dtmaxtemp = merge(rainall,maxtempall.dt,by="Date")
ggplot(data=rain.dtmaxtemp,aes(x=amountcm,y=residuals)) + geom_point() + geom_smooth(method="lm")
#Not much going on with a day-of correlation

A = ggplot(data=rainall[which(rainall$Year=="2016"),],aes(x=jDate,y=amountcm)) + geom_line()
B = ggplot(data=maxtempall.dt[which(maxtempall.dt$Year=="2016"),],aes(x=jDate,y=residuals)) + geom_line()
ggarrange(A,B,ncol=1)

A = ggplot(data=rainall[which(rainall$Year=="2017"),],aes(x=jDate,y=amountcm)) + geom_line()
B = ggplot(data=maxtempall.dt[which(maxtempall.dt$Year=="2017"),],aes(x=jDate,y=residuals)) + geom_line() + 
  geom_vline(aes(xintercept=150)) + geom_vline(aes(xintercept=125)) + geom_vline(aes(xintercept=210),color="red")
ggarrange(A,B,ncol=1)

A = ggplot(data=rainall[which(rainall$Year=="2018"),],aes(x=jDate,y=amountcm)) + geom_line()
B = ggplot(data=maxtempall.dt[which(maxtempall.dt$Year=="2018"),],aes(x=jDate,y=residuals)) + geom_line()
ggarrange(A,B,ncol=1)

A = ggplot(data=rainall[which(rainall$Year=="2019"),],aes(x=jDate,y=amountcm)) + geom_line()
B = ggplot(data=maxtempall.dt[which(maxtempall.dt$Year=="2019"),],aes(x=jDate,y=residuals)) + geom_line()
ggarrange(A,B,ncol=1)

#Rainfall is often proceeded by warmer than average periods. 


A = ggplot(data=gbs, aes(x=jdate,y=dtmaxtemp85.60)) + geom_point(aes(color=Year),size=2) + 
  geom_line(aes(group=Year,color=Year),size=1.3) + 
  theme_cowplot() + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06"))
B = ggplot(data=gbs, aes(x=jdate,y=encountered)) + geom_point(aes(color=Year)) + geom_smooth(aes(group=Year,color=Year)) +
  theme_cowplot() + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06"))
ggarrange(A,B,ncol=1,align="v")

ggplot(data=gbs, aes(x=dtmaxtemp85.60,y=amountcm56.43)) + geom_point(aes(color=Year)) + 
  geom_smooth(method="lm",aes(group=Year,color=Year),se=F) +
  theme_cowplot() + scale_colour_manual(values=c("#D55E00","#0072B2","#009E73","#efbc06")) + 
  xlab(expression(atop("Mean maximum temperature anomaly",paste("85-60 days before observation (C)")))) +
  ylab("Rainfall 56-43 days before observation (cm)") + labs(color="Season")

cor.test(gbs$amountcm56.43,gbs$dtmaxtemp85.60) #r=0.44

#Correlated pretty well with one another - at higher temperatures for 85-60, get more rainfall from 56-43 but not always.





####Plot Fixed effects of final model####
ggarrange(FS.predict.plot,jdate.predict.plot,rain.predict.plot,dtmaxtemp.predict.plot,align="h")






####Plot delta plots for each climate variable####

#Had to plot these individually - when using ggarrange the key obscured part of the graph

A = plotdelta(dataset=EncM.sum.rain.250.1.quad) + 
  ggtitle(expression(atop("Rainfall total 250-1 days before",paste("observation (quadratic)"))))
B = plotdelta(dataset=EncM.sum.rain.250.1.lin) +
  ggtitle(expression(atop("Rainfall total 250-1 days before",paste("observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.mintemp.100.1.quad) + 
  ggtitle(expression(atop("Mean minimum temperature 100-1 days",paste("before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.mintemp.100.1.lin) +
  ggtitle(expression(atop("Mean minimum temperature 100-1 days",paste("before observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.dtmintemp.100.1.quad) +
  ggtitle(expression(atop("Mean minimum temperature anomaly 100-1",paste("days before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.dtmintemp.100.1.lin) +
  ggtitle(expression(atop("Mean minimum temperature anomaly 100-1",paste("days before observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.maxtemp.100.1.quad) +
  ggtitle(expression(atop("Mean maximum temperature 100-1 days",paste("before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.maxtemp.100.1.lin) + 
  ggtitle(expression(atop("Mean maximum temperature 100-1 days",paste("before observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.dtmaxtemp.100.1.quad) +
  ggtitle(expression(atop("Mean maximum temperature anomaly",paste("100-1 days before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.dtmaxtemp.100.1.lin) +
  ggtitle(expression(atop("Mean maximum temperature anomaly",paste("100-1 days before observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.windspeed.10.0.quad) +
  ggtitle(expression(atop("Mean windspeed 10-0 days",paste("before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.windspeed.10.0.lin) +
  ggtitle(expression(atop("Mean windspeed 10-0 days",paste("before observation (linear)"))))
ggarrange(A,B)

A = plotdelta(dataset=EncM.mean.ndvi.63.1.quad) +
  ggtitle(expression(atop("Mean NDVI 63-0 days",paste("before observation (quadratic)"))))
B = plotdelta(dataset=EncM.mean.ndvi.63.1.lin) +
  ggtitle(expression(atop("Mean NDVI 63-0 days",paste("before observation (linear)"))))
ggarrange(A,B)




####Plot randomization histograms####

A = ggplot(data=EncM.sum.rain.100.1.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.sum.rain.250.1.quad$deltaAICc[[2]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Rainfall sum quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
B = ggplot(data=EncM.sum.rain.100.1.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.sum.rain.250.1.lin$deltaAICc[[2]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Rainfall sum linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
C = ggplot(data=EncM.mean.mintemp.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.mintemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean mintemp quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
D = ggplot(data=EncM.mean.mintemp.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.mintemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean mintemp linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
E = ggplot(data=EncM.mean.mintemp.residuals.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmintemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean resid. mintemp quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
G = ggplot(data=EncM.mean.mintemp.residuals.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmintemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean resid. mintemp linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
H = ggplot(data=EncM.mean.maxtemp.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.maxtemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean maxtemp quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
I = ggplot(data=EncM.mean.maxtemp.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.maxtemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean maxtemp linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
J = ggplot(data=EncM.mean.maxtemp.residuals.100.1.quad.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmaxtemp.100.1.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean resid maxtemp quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
K = ggplot(data=EncM.mean.maxtemp.residuals.100.1.lin.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.dtmaxtemp.100.1.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean resid maxtemp linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
L = ggplot(data=EncM.mean.wind.10.0.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.windspeed.10.0.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean windspeed quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
M = ggplot(data=EncM.mean.wind.10.0.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.windspeed.10.0.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean windspeed linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
N = ggplot(data=EncM.mean.ndvi.63.1.quad.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.ndvi.63.1.quad$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean NDVI quadratic") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
O = ggplot(data=EncM.mean.ndvi.63.1.lin.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = EncM.mean.ndvi.63.1.lin$deltaAICc[[1]],color="red",size=1.5) + xlab("deltaAIC") +
  ggtitle("Mean NDVI linear") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))


ggarrange(A,B,C,D,E,G,H,I,J,K,L,M,N,O)




