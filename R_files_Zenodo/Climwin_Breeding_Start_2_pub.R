####Climwin Breeding Start

#Using this script to see how environmental variables influence when breeding begins.

#Ran detrended temps two ways - first was getting residuals based on comparison to each year's average values as in groups per observation analyses. 
#Second way was getting residuals by comparing to average across all years of study. This is the version I used in the paper - makes more sense for
#this analysis that uses an absolute time window to compare variation in timing of first egg dates across years - want to know variation in temp 
#anomalies across years. 


library(here)

####Load data####
firstegg = read.csv(here::here("Output files","firstegg_final.csv"))
rain = read.csv(here::here("Output files","rain_2004-2018.csv"))
mintemp = read.csv(here::here("Output files","mintemp_1994-2020.csv"))
maxtemp = read.csv(here::here("Output files","maxtemp_1994-2020.csv"))
L7ndvi = read.csv(here::here("Output files","Landsat7_all_years.csv"))

#Get dates
firstegg$Clutch.comp.date = as.Date(firstegg$Clutch.comp.date,"%m/%d/%y")
firstegg$First.egg.date = as.Date(firstegg$First.egg.date,"%m/%d/%y")
rain$Date = as.Date(rain$Date,"%m/%d/%y")
mintemp$Date = as.Date(mintemp$Date,"%m/%d/%y")
maxtemp$Date = as.Date(maxtemp$Date,"%m/%d/%y")
L7ndvi$Date = as.Date(L7ndvi$Date,"%m/%d/%y")

#Get factors
firstegg$Year = as.factor(firstegg$Year)
firstegg$Female.fwno = as.factor(firstegg$Female.fwno)
firstegg$Male.fwno = as.factor(firstegg$Male.fwno)

####Plot climate variables####

#Get absolute rainfall window
library(lubridate)
yday(as.Date("2012-08-01")) - 155 # = 59
yday(as.Date("2012-08-01")) - 47 # = 167


#Plot rainfall and first egg dates as a histogram
library(dplyr)
rain.sub = rain %>% filter(Year>=2012 & Year<=2019)
library(ggplot2)
library(cowplot)
ggplot() + geom_histogram(data=firstegg,aes(x=First.egg.jdate),alpha=0.8,binwidth = 7,fill="dark red") +
  geom_line(data=rain.sub,aes(x=jDate,y=amountcm),size=0.5) + theme_cowplot() +
  geom_errorbarh(data=rain.sub,aes(xmax=155,xmin=47,y=15,height=2),size=0.5) +
  facet_wrap(vars(Year),ncol=2) 

#Plot minimum temperature
mintemp.sub = mintemp %>% filter(Year>=2012 & Year<=2019)
ggplot() + 
  geom_line(data=mintemp.sub,aes(x=jdate,y=mintemp)) + facet_wrap(vars(Year),ncol=2) +
   geom_histogram(data=firstegg,aes(x=First.egg.jdate),alpha=0.8,binwidth = 7,fill="dark red") + 
  theme_cowplot() 

#Plot maximum temperature
maxtemp.sub = maxtemp %>% filter(Year>=2012 & Year<=2019)
ggplot() + 
  geom_line(data=maxtemp.sub,aes(x=jdate,y=maxtemp)) + facet_wrap(vars(Year),ncol=2) +
  geom_histogram(data=firstegg,aes(x=First.egg.jdate),alpha=0.8,binwidth = 7,fill="dark red") + 
  theme_cowplot() 

#Plot NDVI
ndvi.sub = L7ndvi %>% mutate(Year=year(Date)) %>% filter(Year>=2012, Year<=2019) %>% mutate(jdate=yday(Date))
ggplot() + 
  geom_line(data=ndvi.sub,aes(x=jdate,y=NDVI*10)) + facet_wrap(vars(Year),ncol=2) + 
  geom_histogram(data=firstegg,aes(x=First.egg.jdate),alpha=0.8,binwidth = 7,fill="dark red") + 
  theme_cowplot() 

#Get mean first egg dates for each year
firstegg.mean


####Baseline models####
library(lme4)
test.model1 = lm(First.egg.jdate~1,data=firstegg)
summary(test.model1)
test.model2 = lmer(First.egg.jdate~1 + (1|Female.fwno), data=firstegg,REML=T)
summary(test.model2)
test.model3 = lmer(First.egg.jdate~1 + (1|Female.fwno) + (1|Male.fwno), data=firstegg,REML=T)
summary(test.model3) #Male.fwno explains no variance in this model
test.model4 = lmer(First.egg.jdate~Female.Age + (1|Female.fwno), data=firstegg,REML=F)
summary(test.model4)
test.model5 = lmer(First.egg.jdate~Female.Age + Helpers + (1|Female.fwno), data=firstegg,REML=F)
summary(test.model5)
test.model6 = lmer(First.egg.jdate~Female.Age + Helpers + Paired.prev + (1|Female.fwno), data=firstegg,REML=F)
summary(test.model6)
test.model7 = lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg,REML=F)
summary(test.model7)

anova(test.model4,test.model5,test.model6,test.model7)



#Result: test.model7 is best
hist(resid(test.model7)) #Looks good

##Why not include year as a random effect? 
#Lv 2019 "When to start and when to stop: Effects of climate on breeding in a multiâ€�brooded songbird" uses an 
#absolute window and includes year as a random effect, but that doesn't seem to make sense. When using an 
#absolute window, each observation within a year has the same climate value. Including year as a random effect
#would just be taking away from any variation that climate value explains. NOTE - I do use year as a random effect
#at the end, but not during the climate window selection - see end of rainfall section for why. 
#Hidalgo Aranzamendi also used year as a random effect but she was looking at within-year effects mainly, 
#so in that case it makes sense.


#NOTE: Did not specify "REML = False" in linear models below, but Climwin sets them to REML = False for you. 
#Should use REML = False when comparing models with different fixed effects. See climwin advanced vignette. 


####200-1 rainfall sum model#### 

#Test what rainfall period during the 200-1 days before the observation best predicts 
#first egg dates. Ran 365 initially and nothing past 
#200 days was meaningful
library(climwin)

#This is what example data from climwin paper look like - used for absolute window
climwin::Chaff

# FE.rain.sum.200.1.abs.lin.k <- slidingwin(k=6,
#                                         xvar=list(Rain = rain$amountcm),
#                                        cdate = rain$Date,
#                                        bdate = firstegg$First.egg.date,
#                                        baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                        type = "absolute",
#                                        cinterval = "day",
#                                        range = c(200,1),
#                                         refday = c(1,8),
#                                        exclude = c(20,-1),
#                                        stat = "sum",
#                                        fun = "lin")

#Result: Lin and quad are about the same delta AIC but lin makes much more sense 
#(155-47 days for lin vs 57-49 days for quad)

#Save model file
#saveRDS(FE.rain.sum.200.1.abs.lin.k,here::here("RDS model outputs","FE_rain_sum_200_1_abs_lin_k_20exclude.rds"))

#Read in model file
FE.rain.sum.200.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_rain_sum_200_1_abs_lin_k_20exclude.rds"))
#Look to see what windows are of best model 
head(FE.rain.sum.200.1.abs.lin.k[[1]]$Dataset,n=20)
#Get dataset
FE.rain.sum.200.1.abs.lin.k.dataset = FE.rain.sum.200.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.rain.sum.200.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.rain.sum.200.1.abs.lin.k.dataset)
#Plotall
plotall(FE.rain.sum.200.1.abs.lin.k.dataset)
#Plot relationship - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.rain.sum.200.1.abs.lin.k.dataset,
         bestmodel = FE.rain.sum.200.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.rain.sum.200.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.rain.sum.200.1.abs.lin.k.bestdata = FE.rain.sum.200.1.abs.lin.k[[1]]$BestModelData
FE.rain.sum.200.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.rain.sum.200.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate)  + 
                                      (1|Female.fwno), data=FE.rain.sum.200.1.abs.lin.k.bestdata)
summary(FE.rain.sum.200.1.abs.lin.k.m1)

library(lmerTest)
FE.rain.sum.200.1.abs.lin.k.m2 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate) +
                                      (1|Female.fwno), data=FE.rain.sum.200.1.abs.lin.k.bestdata)
summary(FE.rain.sum.200.1.abs.lin.k.m2)



####Randomizations for rainfall sum####

#Base randwin code swaps biological dates among observations, but that would be swapping dates among years.
#These randomizations need to test the idea that a year's rainfall specifically predicts when breeding begins that
#year. So a set of first egg dates are specifically related to a year's rainfall. If you mix up the breeding 
#start dates among years (keeping yearly data together), you shouldn't see a significant relationship between
#rainfall and breeding start if the observed relationship is meaningful. 
#So kept first egg dates together within years but mixed up which rainfall year they were associated with and
#coded such that no year could get assigned to itself. 

#For rainfall, since all  models of observed data within two AIC scores of the top model had a negative beta
#(coefficient), did randomizations to test what the chance of getting a negative relationship between rainfall 
#and breeding start was. 

#Randomization script:
#Cluster_breeding_rainfall_sum_abs_200_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.rain.sum.200.1.abs.lin.dataset.k.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_rainfall_sum_abs_200_1_lin_k_negbeta_only_randomizations_20exclude.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.rain.sum.200.1.abs.lin.k.dataset,FE.rain.sum.200.1.abs.lin.dataset.k.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.rain.sum.200.1.abs.lin.dataset.k.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.rain.sum.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()


#Results: A negative relationship between rainfall and breeding start is different from random. 


###Why not to use Year as a random effect in climate window selection:
# test <-  slidingwin(k=6,
#                     xvar=list(Rain = rain$amountcm),
#                     cdate = rain$Date,
#                     bdate = firstegg$First.egg.date,
#                     baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Year) + (1|Female.fwno), data=firstegg),
#                     type = "absolute",
#                     cinterval = "day",
#                     range = c(155,47),
#                      refday = c(1,8),
#                     exclude = c(100,-1),
#                     stat = "sum",
#                     fun = "lin")
# test2 = test[[1]]$Dataset
# test3 = test[[1]]$BestModelData

#Delta-AIC value is now only -2, because each data point within a year has the same climate value, and then gets the same year
#value. They're representing almost the same thing. Do include year in final model. 




####200-1 mean minimum temperature model#### 

# FE.mintemp.mean.200.1.abs.lin.k <- slidingwin(k=6,
#                                         xvar=list(temp = mintemp$mintemp),
#                                        cdate = mintemp$Date,
#                                        bdate = firstegg$First.egg.date,
#                                        baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                        type = "absolute",
#                                        cinterval = "day",
#                                        range = c(200,1),
#                                         refday = c(1,8),
#                                        exclude = c(10,-1),
#                                        stat = "mean",
#                                        fun = "lin")

#Save model file
# saveRDS(FE.mintemp.mean.200.1.abs.lin.k,here::here("RDS model outputs","FE_mintemp_mean_200_1_abs_lin_k_10exclude.rds"))

#Read in model file
FE.mintemp.mean.200.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_mintemp_mean_200_1_abs_lin_k_10exclude.rds"))
#Look to see what windows are of best model 
head(FE.mintemp.mean.200.1.abs.lin.k[[1]]$Dataset,n=20)
#Get dataset
FE.mintemp.mean.200.1.abs.lin.k.dataset = FE.mintemp.mean.200.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.mintemp.mean.200.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.mintemp.mean.200.1.abs.lin.k.dataset)
#Plotall
plotall(FE.mintemp.mean.200.1.abs.lin.k.dataset)
#Plot relationship between rainfall and first egg dates - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.mintemp.mean.200.1.abs.lin.k.dataset,
         bestmodel = FE.mintemp.mean.200.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.mintemp.mean.200.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.mintemp.mean.200.1.abs.lin.k.bestdata = FE.mintemp.mean.200.1.abs.lin.k[[1]]$BestModelData
FE.mintemp.mean.200.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.mintemp.mean.200.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate)  + 
                                      (1|Female.fwno), data=FE.mintemp.mean.200.1.abs.lin.k.bestdata)
summary(FE.mintemp.mean.200.1.abs.lin.k.m1)



####Randomizations for mintemp mean####

#For mintemp, models of observed data within two AIC scores of the top model included both negative and positive
#betas (coefficient), so did randomizations to test what the chance of getting a either a negative 
#or positive relationship between mintemp and breeding start was. Also ran it at positive only and observed still 
#wasn't more than random chance

#Randomization script:
#Cluster_breeding_mintemp_mean_abs_200_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.mintemp.mean.200.1.abs.lin.k.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_mintemp_mean_abs_200_1_lin_k_randomizations_10exclude.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.mintemp.mean.200.1.abs.lin.k.dataset,FE.mintemp.mean.200.1.abs.lin.k.dataset.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.mintemp.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.mintemp.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()

#Results: mintemp mean is not different from random





####200-1 mean maximum temperature model#### 

# FE.maxtemp.mean.200.1.abs.lin.k <- slidingwin(k=6,
#                                               xvar=list(temp = maxtemp$maxtemp),
#                                             cdate = maxtemp$Date,
#                                             bdate = firstegg$First.egg.date,
#                                             baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                             type = "absolute",
#                                             cinterval = "day",
#                                             range = c(200,1),
#                                             refday = c(1,8),
#                                             exclude = c(10,-1),
#                                             stat = "mean",
#                                             fun = "lin")

#Save model file
#saveRDS(FE.maxtemp.mean.200.1.abs.lin.k,here::here("RDS model outputs","FE_maxtemp_mean_200_1_abs_lin_k_10exclude.rds"))

#Read in model file
FE.maxtemp.mean.200.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_maxtemp_mean_200_1_abs_lin_k_10exclude.rds"))
#Look to see what windows are of best model 
head(FE.maxtemp.mean.200.1.abs.lin.k[[1]]$Dataset,n=20)
#Get dataset
FE.maxtemp.mean.200.1.abs.lin.k.dataset = FE.maxtemp.mean.200.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.maxtemp.mean.200.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.maxtemp.mean.200.1.abs.lin.k.dataset)
#Plotall
plotall(FE.maxtemp.mean.200.1.abs.lin.k.dataset)
#Plot relationship between rainfall and first egg dates - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.maxtemp.mean.200.1.abs.lin.k.dataset,
         bestmodel = FE.maxtemp.mean.200.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.maxtemp.mean.200.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.maxtemp.mean.200.1.abs.lin.k.bestdata = FE.maxtemp.mean.200.1.abs.lin.k[[1]]$BestModelData
FE.maxtemp.mean.200.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.maxtemp.mean.200.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate)  + 
                                          (1|Female.fwno), data=FE.maxtemp.mean.200.1.abs.lin.k.bestdata)
summary(FE.maxtemp.mean.200.1.abs.lin.k.m1)


####Randomizations for maxtemp mean####

#For maxtemp, models of observed data within two AIC scores of the top model included only positive
#betas (coefficient), so did randomizations to test what the chance of getting
#a positive relationship between maxtemp and breeding start was.

#Randomization script:
#Cluster_breeding_maxtemp_mean_abs_200_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.maxtemp.mean.200.1.abs.lin.k.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_maxtemp_mean_abs_200_1_lin_k_posbeta_only_randomizations_10exclude.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.maxtemp.mean.200.1.abs.lin.k.dataset,FE.maxtemp.mean.200.1.abs.lin.k.dataset.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.maxtemp.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.maxtemp.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()

#Results: Maxtemp not different from random






####NDVI model####

#Using Landsat 7

# FE.ndvi.mean.100.1.abs.lin.k <- slidingwin(k=6,
#                                         xvar=list(ndvi = L7ndvi$NDVI),
#                                        cdate = L7ndvi$Date,
#                                        bdate = firstegg$First.egg.date,
#                                        baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                        type = "absolute",
#                                        cinterval = "day",
#                                        range = c(100,1),
#                                        exclude = c(5,-1),
#                                         refday = c(1,8),
#                                        stat = "mean",
#                                        fun = "lin")


#Save model file
#saveRDS(FE.ndvi.mean.100.1.abs.lin.k,here::here("RDS model outputs","FE_ndvi_mean_100_1_abs_lin_k.rds"))

#Read in model file
FE.ndvi.mean.100.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_ndvi_mean_100_1_abs_lin_k.rds"))
#Look to see what windows are of best model 
head(FE.ndvi.mean.100.1.abs.lin.k[[1]]$Dataset,n=10)
#Get dataset
FE.ndvi.mean.100.1.abs.lin.k.dataset = FE.ndvi.mean.100.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.ndvi.mean.100.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.ndvi.mean.100.1.abs.lin.k.dataset)
#Plotall
plotall(FE.ndvi.mean.100.1.abs.lin.k.dataset)
#Plot relationship between rainfall and first egg dates - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.ndvi.mean.100.1.abs.lin.k.dataset,
         bestmodel = FE.ndvi.mean.100.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.ndvi.mean.100.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.ndvi.mean.100.1.abs.lin.k.bestdata = FE.ndvi.mean.100.1.abs.lin.k[[1]]$BestModelData
FE.ndvi.mean.100.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.ndvi.mean.100.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate)  + 
                                      (1|Female.fwno), data=FE.ndvi.mean.100.1.abs.lin.k.bestdata)
summary(FE.ndvi.mean.100.1.abs.lin.k.m1)

library(lmerTest)
FE.ndvi.mean.100.1.abs.lin.k.m2 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate) +
                                      (1|Female.fwno), data=FE.ndvi.mean.100.1.abs.lin.k.bestdata)
summary(FE.ndvi.mean.100.1.abs.lin.k.m2)



####Randomizations for NDVI mean####

#Randomization script:
#Cluster_breeding_NDVI_mean_abs_100_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.ndvi.mean.100.1.abs.lin.k.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_ndvi_mean_abs_100_1_lin_posbeta_only_randomizations_k.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.ndvi.mean.100.1.abs.lin.k.dataset,FE.ndvi.mean.100.1.abs.lin.k.dataset.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.ndvi.mean.100.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.ndvi.mean.100.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()

#Results: NDVI is not different from random






####Detrended across years mean mintemp####

#This test is comparing breeding start across years, try looking at how each year's temperatures compare 
#to average temp of the studied years. In groups encountered models those were mainly looking at within-years bc
#relative window. 

#Get an average mintemp for all years of study
mintemp.sub$all.fitted = loess(mintemp~jdate,data=mintemp.sub,span=0.5)$fitted

#Show what that looks like
mintemp.sub$Year = as.factor(mintemp.sub$Year)
ggplot(mintemp.sub,aes(x=jdate,y=mintemp)) + geom_line(aes(group=Year,color=Year)) + 
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="black") + theme_cowplot()

#Get residuals for each year
mintemp.sub$residual = loess(mintemp~jdate,data=mintemp.sub,span=0.5)$residual
#If you order by jdate you can see that each date gets calculated independently, not repeating values like "fitted". 


#Sliding window analysis for detrended mean mintemp compared to across-year average: 

# FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k <- slidingwin(k=6,
#                                       xvar=list(temp = mintemp.sub$residual),
#                                        cdate = mintemp.sub$Date,
#                                        bdate = firstegg$First.egg.date,
#                                        baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                        type = "absolute",
#                                        cinterval = "day",
#                                        range = c(200,1),
#                                         refday = c(1,8),
#                                        exclude = c(10,-1),
#                                        stat = "mean",
#                                        fun = "lin")

#Save model file
#saveRDS(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k,here::here("RDS model outputs","FE_mintemp_detrended_allyears_mean_200_1_abs_lin_k_10exclude.rds"))

#Read in model file
FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_mintemp_detrended_allyears_mean_200_1_abs_lin_k_10exclude.rds"))
#Look to see what windows are of best model 
head(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$Dataset,n=10)
#Get dataset
FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset)
#Plotall
plotall(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset)
#Plot relationship between rainfall and first egg dates - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset,
         bestmodel = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModelData
FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate) +
                                                      (1|Female.fwno), data=FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata)
summary(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.m1)



####Randomizations for detrended across years mintemp mean####

#For detrended mintemp compared to all years, models of observed data within two AIC scores of the top model included both negative and positive
#betas (coefficient), so did randomizations to test what the chance of getting a either a negative 
#or positive relationship between detrended mintemp and breeding start was. 

#Randomization script:
#Cluster_breeding_detrended_mintemp_allyears_mean_abs_200_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_detrended_mintemp_allyears_mean_abs_200_1_lin_k_randomizations_10exclude.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset,FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()







####Detrended across years mean maxtemp####

#This this test is comparing breeding start across years, try looking at how each year's temperatures compare 
#to average temp of the studied years. In groups encountered models those were mainly looking at within-years bc
#relative window. 

#Get an average maxtemp for all years of study
maxtemp.sub$all.fitted = loess(maxtemp~jdate,data=maxtemp.sub,span=0.5)$fitted

#Show what that looks like
maxtemp.sub$Year = as.factor(maxtemp.sub$Year)
ggplot(maxtemp.sub,aes(x=jdate,y=maxtemp)) + geom_line(aes(group=Year,color=Year)) + 
  geom_line(aes(x=jdate,y=all.fitted,group=Year),size=1,color="black") + theme_cowplot()

#Get residuals for each year
maxtemp.sub$residual = loess(maxtemp~jdate,data=maxtemp.sub,span=0.5)$residual
#If you order by jdate you can see that each date gets calculated independently, not repeating values like "fitted". 


#Sliding window analysis for detrended mean maxtemp compared to across-year average: 

# FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k <- slidingwin(k=6,
#                                       xvar=list(temp = maxtemp.sub$residual),
#                                        cdate = maxtemp.sub$Date,
#                                        bdate = firstegg$First.egg.date,
#                                        baseline = lme4::lmer(First.egg.jdate~Helpers + Paired.prev + (1|Female.fwno), data=firstegg),
#                                        type = "absolute",
#                                        cinterval = "day",
#                                        range = c(200,1),
#                                         refday = c(1,8),
#                                        exclude = c(10,-1),
#                                        stat = "mean",
#                                        fun = "lin")

#Save model file
#saveRDS(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k,here::here("RDS model outputs","FE_maxtemp_detrended_allyears_mean_200_1_abs_lin_k_10exclude.rds"))

#Read in model file
FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k = readRDS(here::here("Output files","Climwin model outputs","FE_maxtemp_detrended_allyears_mean_200_1_abs_lin_k_10exclude.rds"))
#Look to see what windows are of best model 
head(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$Dataset,n=50)
#Get dataset
FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset)
#Plotall
plotall(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset)
#Plot relationship between rainfall and first egg dates - this function doesn't work well for mixed models
#but can give an idea of the relationship trend
plotbest(dataset = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset,
         bestmodel = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModel, 
         bestmodeldata = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModelData)
#Get best model data
FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k[[1]]$BestModelData
FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata$Year = firstegg$Year


FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.m1 = lmer(yvar~scale(Helpers) + Paired.prev + scale(climate) +
                                                               (1|Female.fwno), data=FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.bestdata)
summary(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.m1)



####Randomizations for detrended across years maxtemp mean####

#For detrended mintemp, models of observed data within two AIC scores of the top model included positive
#betas (coefficient) only, so did randomizations to test what the chance of getting 
#a positive relationship between detrended mintemp and breeding start was. 

#Randomization script:
#Cluster_breeding_detrended_maxtemp_allyears_mean_abs_200_1_rand_lin_k.R

#Read in Delta-AIC values saved in randomization script
FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand = read.csv(here::here("Output files","Climwin randomizations","Breeding_detrended_maxtemp_allyears_mean_abs_200_1_lin_k_posbeta_only_randomizations_10exclude.csv"))

#Calculate pvalues - sample size only important for 
pvalue(FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset,FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,"AIC",sample.size = 8)

#Histograms
ggplot(data=FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Linear") + ylab("Frequency") + theme_cowplot()





####Final Model####

#Only climate variable that was meaningful for predicting breeding start date was rainfall 155-47 days before Aug 1
#That is rainfall from Feb 27th to June 15th (julian days 58-166)
as.Date("2018-08-01") - 155
as.Date("2018-08-01") - 47

head(FE.rain.sum.200.1.abs.lin.k.bestdata)
firstegg$amountcm.155.47 = FE.rain.sum.200.1.abs.lin.k.bestdata$climate

#Make paired previously numerical for confidence interval calculating - if it's a factor I get an error. 
#This does not change the output of the model or change paired prev's coefficient. 
firstegg$Paired.prev2 = NA
for (i in 1:nrow(firstegg)) {
  if(firstegg$Paired.prev[i]=="Yes") {firstegg$Paired.prev2[i]=1} else
    {firstegg$Paired.prev2[i]=0}
}

#Add year (season) into the final model as a random effect to account for other differences across years
#other than rainfall
final.model = lmer(First.egg.jdate~Helpers + Paired.prev2 + amountcm.155.47 + (1|Female.fwno) + (1|Year), data=firstegg, REML=F)
summary(final.model)

#Check residuals - Look good
qqnorm(resid(final.model))
qqline(resid(final.model))
hist(resid(final.model))

#Check for colinearity
library(car)
car::vif(final.model) 

library(DHARMa)
final.resid = simulateResiduals(final.model)
plot(final.resid)
#Not perfect but decent.



####R2 values for final model
library(piecewiseSEM)
rsquared(final.model) #Marginal R2 = 0.47

##Get delta-R2 values of individual fixed effects

#R2 for helpers
final.model.nohelpers = lmer(First.egg.jdate~Paired.prev2 + amountcm.155.47 + (1|Female.fwno) + (1|Year), data=firstegg, REML=F)
(rsquared(final.model.nohelpers)$Marginal - rsquared(final.model)$Marginal) *100

#R2 for paired previously
final.model.nopprev = lmer(First.egg.jdate~ Helpers + amountcm.155.47 + (1|Female.fwno) + (1|Year), data=firstegg, REML=F)
(rsquared(final.model.nopprev)$Marginal - rsquared(final.model)$Marginal) *100

#R2 for rainfall - amountcm.155.47
final.model.norain = lmer(First.egg.jdate~ Helpers + Paired.prev2 + (1|Female.fwno) + (1|Year), data=firstegg, REML=F)
(rsquared(final.model.norain)$Marginal - rsquared(final.model)$Marginal) *100





####Predict
#Get prediction dataframe - need to have columns for each of the fixed effects here
#Use mean number of helpers and paired.prev = yes

FE.rain.final.df = data.frame(expand.grid(amountcm.155.47=seq(27,60,by=1),
                                          Helpers=mean(firstegg$Helpers),
                                          Paired.prev2=1))

FE.rain.final = predict(final.model,
                        newdata=FE.rain.final.df,
                        type="response",re.form=NA)

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

#Get confidence intervals
FE.rain.final.CI <- easyPredCI(final.model,
                                 newdata=FE.rain.final.df)

#Combine prediction dataframes
FE.rain.final.df2 = cbind(FE.rain.final.df,
                            First.egg.jdate=FE.rain.final, FE.rain.final.CI)

##Plot - for jittered version look in plots old
ggplot(data=firstegg,aes(x=amountcm.155.47,y=First.egg.jdate)) +
  geom_point(size=3.1,aes(fill=Year),shape=21,color="#0f0f52") +
  geom_ribbon(data=FE.rain.final.df2,
              aes(x=amountcm.155.47,y=First.egg.jdate,ymin=conf.low,ymax=conf.high),alpha=0.2) +
  geom_line(data=FE.rain.final.df2 ,aes(x=amountcm.155.47,y=First.egg.jdate),size=1.1) + 
  theme_cowplot() + xlab(expression(atop("Total rainfall (cm) between February 27", paste("and June 15 (days 58-166)")))) + 
  ylab("First egg date (julian date)") + labs(fill = "Season") +
  scale_fill_manual(values=c("#615c5c", "#CC79A7", "#56B4E9", "#E69F00", "#D55E00", "#0072B2", "#009E73", "#efbc06"))
  

##Plot again but change plotting order of points to try and see other years more clearly. 
set.seed(42)
rows <- sample(nrow(firstegg))
firstegg2 <- firstegg[rows,]

ggplot(data=firstegg2,aes(x=amountcm.155.47,y=First.egg.jdate)) +
  geom_point(size=3.1,aes(fill=Year),shape=21,color="#0f0f52") +
  geom_ribbon(data=FE.rain.final.df2,
              aes(x=amountcm.155.47,y=First.egg.jdate,ymin=conf.low,ymax=conf.high),alpha=0.2) +
  geom_line(data=FE.rain.final.df2 ,aes(x=amountcm.155.47,y=First.egg.jdate),size=1.1) + 
  theme_cowplot() + xlab(expression(atop("Total rainfall (cm) between", paste("February 27 and June 15")))) + 
  ylab("First egg date (day of year)") + labs(fill = "Season") +
  scale_fill_manual(values=c("#615c5c", "#CC79A7", "#56B4E9", "#E69F00", "#D55E00", "#0072B2", "#009E73", "#efbc06"))





####Plot rainfall by date graphs again####

ggplot() + geom_histogram(data=firstegg,aes(x=First.egg.jdate,fill=Year),binwidth = 7) +
  geom_line(data=rain.sub,aes(x=jDate,y=amountcm),size=0.5) + theme_cowplot() +
  geom_errorbarh(data=rain.sub,aes(xmax=166,xmin=58,y=15,height=2),size=0.5) +
  facet_wrap(vars(Year),ncol=2,scales="free") +
  scale_fill_manual(values=c("#615c5c", "#CC79A7", "#56B4E9", "#E69F00", "#D55E00", "#0072B2", "#009E73", "#efbc06")) +
  theme(legend.position = "none") + xlab("Julian Date") + ylab("Rainfall (cm) / Number of first egg dates")
  #theme(strip.background = element_rect(color="black",fill=NA, linetype="solid"))



####Plot Delta AIC window graphs again####

plotdelta(dataset=FE.rain.sum.200.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Rainfall sum 200-1 days",paste("before August 1st"))))
plotdelta(dataset=FE.mintemp.mean.200.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Mean minimum temperature 200-1 days",paste("before August 1st"))))
plotdelta(dataset=FE.mintemp.detrended.mean.200.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Detrended mean minimum temperature",paste("200-1 days before August 1st"))))
plotdelta(dataset=FE.ndvi.mean.100.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Mean NDVI 200-1 days",paste("before August 1st"))))

#Across years: comparing year to average across years - USE THESE VERSIONS IN PAPER
plotdelta(dataset=FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Mean minimum temperature anomaly",paste("200-1 days before August 1st"))))
plotdelta(dataset=FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset) +
  ggtitle(expression(atop("Mean maximum temperature anomaly",paste("200-1 days before August 1st"))))



####Plot randomization comparsisons again####
A = ggplot(data=FE.rain.sum.200.1.abs.lin.dataset.k.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.rain.sum.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Rainfall sum linear p<0.001") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
B = ggplot(data=FE.mintemp.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.mintemp.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Mean mintemp linear p=0.83") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
C = ggplot(data=FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.mintemp.detrended.allyears.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Mean detrended mintemp linear p=0.89") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
D = ggplot(data=FE.maxtemp.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.maxtemp.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Mean maxtemp linear p=0.25") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
E = ggplot(data=FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.maxtemp.detrended.allyears.mean.200.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Mean detrended maxtemp linear p=0.37") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))
G = ggplot(data=FE.ndvi.mean.100.1.abs.lin.k.dataset.rand,aes(x=deltaAICc)) + 
  geom_histogram(binwidth = 1,fill="gray",color="black") + 
  geom_vline(xintercept = FE.ndvi.mean.100.1.abs.lin.k.dataset$deltaAICc[[1]],color="red",size=1.5) +
  ggtitle("Mean NDVI linear p=0.12") + ylab("Frequency") + theme_cowplot() + theme(plot.title = element_text(size = 10))

ggarrange(A,B,C,D,E,G)





