####Sliding window analyses for plumage type

##Looking at whether/how short-term climate variation influences whether males molt into nuptial plumage

library(here)

####Load climate data####

#Load rainfall data
rain = read.csv(here::here("Input files","rain_2004-2018.csv"))
rain$Date = as.Date(rain$Date,"%m/%d/%y")

#Load mintemp data
mintemp = read.csv(here::here("Input files","mintemp_1994-2020.csv"))
mintemp$Date = as.Date(mintemp$Date,"%m/%d/%y")

#Load maxtemp data
maxtemp = read.csv(here::here("Input files","maxtemp_1994-2020.csv"))
maxtemp$Date = as.Date(maxtemp$Date,"%m/%d/%y")

#Load NDVI data
ndvi = read.csv(here::here("Input files","Landsat7_all_years.csv"))
ndvi$Date = as.Date(ndvi$Date,"%m/%d/%y")


####Using paired_1yo_males_fbb4.csv dataset which includes 112 individuals. Could use full dataset of plumages which would be 162 indivdiduals
#but this way the baseline model can be a bit more defined than a baseline model with no fixed effects. The ratio of plumage types is nearly
#the same across these two datasets - look at Whether a 1-year-old male molted script at end for those plots.  




####1-year-old male sliding windows####

#Read in one-year-old male file that contains information on whether the female has bred previously/pairing method
pag1 = read.csv(here::here("Output files","paired_1yo_males_fbb4.csv"))
pag1$Year = as.factor(pag1$Year)
pag1$Plumage = factor(pag1$Plumage,levels=c("Dull","Interm","Bright"),ordered = T) #Set levels of plumage factor

##Need to set a biological date that tells the sliding window analysis what year each observation belongs to
#Year is year that the breeding season ended in, so year of the January that the season ended in
#So year - 1 is the actual year most of the data were collected and the climate year the data correspond to.
#Climwin code will take the year and re-assign the month/day to the reference date. 
pag1$bdate = paste("10/1/",as.numeric(as.character(pag1$Year))-1,sep="")
pag1$bdate = as.Date(pag1$bdate,"%m/%d/%Y")

#Load packages
library(dplyr)
library(climwin)
library(effects)
library(cowplot)


#Originally thought about running sliding windows at 486-0 days before Dec 31st, which would be Sept 1 - Dec 31st to get all breeding 
#season weather as well as non-breeding season weather. However, many of the males would not have hatched by Sept 1st, we know 
#they can hatch from September - January, probably few in February, so including early breeding season weather probably does not 
#make much biological sense. Want to include weather during a time when all individuals were out of the nest and susceptible to weather
#conditions. Therefore, look at 333 - 0 days which would be February 1st to Dec 31st. Almost all males should hatch by February. 
#So looking at whether non-breeding climate conditions into breeding climate during their first breeding season influences plumage type. 

#Using absolute window because I'm interested in whether wetter years are associated with more 1-year-old males
#molting into nuptial plumage. Not interested in exact timing of molt for this analysis. 
#Using Dec 31st as a reference date because by that point all males that are going to molt have molted, but
#some males don't achieve bright until December, so that means any rainfall leading up to then could be 
#important in determining variation in how many males molt across years. 

##Not including year as a random effect in the sliding window analysis since in an absolute window each datapoint within
#a year has the same climate value, essentially coding almost the same information. Can include year in final model
#to control for any remaining variation across years. 




####1yo rainfall window####

##Running rainfall for quadratic gave similar results - 112 to 89 window with lower AIC than linear. 

# pag1.rain.lin = slidingwin(xvar=list(Rain = rain$amountcm),
#                        cdate = rain$Date,
#                        bdate = pag1$bdate, #This will get the correct year
#                        baseline = glm(Ornamented~Bredb4,data=pag1,family="binomial"),
#                        type = "absolute",
#                        cinterval = "day",
#                        range = c(333,0),
#                        refday = c(31,12),
#                        exclude = c(20,-1),
#                        stat = "sum",
#                        fun = "lin")

##Save model file
#saveRDS(pag1.rain.lin,here::here("Output files","pag1_rain_lin.rds"))

##Read in model file
pag1.rain.lin = readRDS(here::here("Output files","pag1_rain_lin.rds"))
#Look to see what windows are of best model 
head(pag1.rain.lin[[1]]$Dataset,n=20)
#Get dataset
pag1.rain.lin.dataset = pag1.rain.lin[[1]]$Dataset
#Check summary
summary(pag1.rain.lin[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=pag1.rain.lin.dataset)
#Plotall
#plotall(pag1.rain.lin.dataset)
#Get best model data
pag1.rain.lin.bestdata = pag1.rain.lin[[1]]$BestModelData
pag1.rain.lin.bestdata$Year = pag1$Year
#Confirm that above model summary is correct by running model outside of climwin with climate data it found. 
pag1.rain.lin.m = glm(yvar~Bredb4 + climate,data=pag1.rain.lin.bestdata,family="binomial")
summary(pag1.rain.lin.m)

##Plot effects of the best model
library(MASS)
plot(Effect(focal.predictors = "Bredb4",pag1.rain.lin.m))
plot(Effect(focal.predictors = "climate",pag1.rain.lin.m))
detach("package:MASS") #interferes with dplyr


####Randomizations - Before was randomizing so that birds within a year stayed together. But in this analysis, since response variable
#is binomial, randomizations allowing ratios to change across years will produce a wider range of possible values. Am interested in
#testing relationship between proportion molted and climate in year but this should give me better expectation of random without losing
#all potential variation across years. 


####Randomizations for rainfall

##Load randomization data from runs on cluster
rain.rand = read.csv(here::here("Output files","Rainfall_sum_1yo_333_0_lin_randomizations_causes_paper_v4.csv"))

##Get p-value
pvalue(dataset=pag1.rain.lin.dataset,datasetrand=rain.rand,metric="AIC",sample.size=8)

##Plot randomizations
ggplot(data=rain.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.rain.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")







####1yo mintemp window####

# pag1.mintemp.lin = slidingwin(xvar=list(Mintemp = mintemp$mintemp),
#                            cdate = mintemp$Date,
#                            bdate = pag1$bdate, #This will get the correct year
#                            baseline = glm(Ornamented~Bredb4,data=pag1,family="binomial"),
#                            type = "absolute",
#                            cinterval = "day",
#                            range = c(333,0),
#                            refday = c(31,12),
#                            exclude = c(20,-1),
#                            stat = "mean",
#                            fun = "lin")

#Save model file
#saveRDS(pag1.mintemp.lin,here::here("Output files","pag1_mintemp_lin.rds"))

#Read in model file
pag1.mintemp.lin = readRDS(here::here("Output files","pag1_mintemp_lin.rds"))
#Look to see what windows are of best model 
head(pag1.mintemp.lin[[1]]$Dataset,n=20)
#Get dataset
pag1.mintemp.lin.dataset = pag1.mintemp.lin[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(pag1.mintemp.lin[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=pag1.mintemp.lin.dataset)
#Plotall
#plotall(pag1.mintemp.lin.dataset)
#Get best model data
pag1.mintemp.lin.bestdata = pag1.mintemp.lin[[1]]$BestModelData
pag1.mintemp.lin.bestdata$Year = pag1$Year
#Confirm that above model summary is correct by running model outside of climwin with climate data it found. 
pag1.mintemp.lin.m = glm(yvar~Bredb4 + climate,data=pag1.mintemp.lin.bestdata,family="binomial")
summary(pag1.mintemp.lin.m)

##Plot effects of the best model
library(MASS)
plot(Effect(focal.predictors = "Bredb4",pag1.mintemp.lin.m))
plot(Effect(focal.predictors = "climate",pag1.mintemp.lin.m))
detach("package:MASS") #interferes with dplyr


####Randomizations for mintemp

##Load randomization data from runs on cluster
mintemp.rand = read.csv(here::here("Output files","Mintemp_mean_1yo_333_0_lin_randomizations_causes_paper_v4.csv"))

##Get p-value
pvalue(dataset=pag1.mintemp.lin.dataset,datasetrand=mintemp.rand,metric="AIC",sample.size=8)

##Plot randomizations
ggplot(data=mintemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.mintemp.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")






####1yo maxtemp window####

# pag1.maxtemp.lin = slidingwin(xvar=list(maxtemp = maxtemp$maxtemp),
#                               cdate = maxtemp$Date,
#                               bdate = pag1$bdate, #This will get the correct year
#                               baseline = glm(Ornamented~Bredb4,data=pag1,family="binomial"),
#                               type = "absolute",
#                               cinterval = "day",
#                               range = c(333,0),
#                               refday = c(31,12),
#                               exclude = c(20,-1),
#                               stat = "mean",
#                               fun = "lin")

#Save model file
#saveRDS(pag1.maxtemp.lin,here::here("Output files","pag1_maxtemp_lin.rds"))

#Read in model file
pag1.maxtemp.lin = readRDS(here::here("Output files","pag1_maxtemp_lin.rds"))
#Look to see what windows are of best model 
head(pag1.maxtemp.lin[[1]]$Dataset,n=20)
#Get dataset
pag1.maxtemp.lin.dataset = pag1.maxtemp.lin[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(pag1.maxtemp.lin[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=pag1.maxtemp.lin.dataset)
#Plotall
#plotall(pag1.maxtemp.lin.dataset)
#Get best model data
pag1.maxtemp.lin.bestdata = pag1.maxtemp.lin[[1]]$BestModelData
pag1.maxtemp.lin.bestdata$Year = pag1$Year
#Confirm that above model summary is correct by running model outside of climwin with climate data it found. 
pag1.maxtemp.lin.m = glm(yvar~Bredb4 + climate,data=pag1.maxtemp.lin.bestdata,family="binomial")
summary(pag1.maxtemp.lin.m)

##Plot effects of the best model
library(MASS)
plot(Effect(focal.predictors = "Bredb4",pag1.maxtemp.lin.m))
plot(Effect(focal.predictors = "climate",pag1.maxtemp.lin.m))
detach("package:MASS") #interferes with dplyr


####Randomizations for maxtemp

##Load randomization data from runs on cluster
maxtemp.rand = read.csv(here::here("Output files","Maxtemp_mean_1yo_333_0_lin_randomizations_causes_paper_v4.csv"))

##Get p-value
pvalue(dataset=pag1.maxtemp.lin.dataset,datasetrand=maxtemp.rand,metric="AIC",sample.size=8)

##Plot randomizations
ggplot(data=maxtemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.maxtemp.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")




####1yo NDVI window####

# pag1.ndvi.lin = slidingwin(xvar=list(NDVI = ndvi$NDVI),
#                               cdate = ndvi$Date,
#                               bdate = pag1$bdate, #This will get the correct year
#                               baseline = glm(Ornamented~Bredb4,data=pag1,family="binomial"),
#                               type = "absolute",
#                               cinterval = "day",
#                               range = c(333,0),
#                               refday = c(31,12),
#                               exclude = c(20,-1),
#                               stat = "mean",
#                               fun = "lin")

#Save model file
#saveRDS(pag1.ndvi.lin,here::here("Output files","pag1_ndvi_lin.rds"))

#Read in model file
pag1.ndvi.lin = readRDS(here::here("Output files","pag1_ndvi_lin.rds"))
#Look to see what windows are of best model 
head(pag1.ndvi.lin[[1]]$Dataset,n=20)
#Get dataset
pag1.ndvi.lin.dataset = pag1.ndvi.lin[[1]]$Dataset
#Check summary - Model is unidentifiable - will scale variables below
summary(pag1.ndvi.lin[[1]]$BestModel)
#Plot delta plot
plotdelta(dataset=pag1.ndvi.lin.dataset)
#Plotall
#plotall(pag1.ndvi.lin.dataset)
#Get best model data
pag1.ndvi.lin.bestdata = pag1.ndvi.lin[[1]]$BestModelData
pag1.ndvi.lin.bestdata$Year = pag1$Year
#Confirm that above model summary is correct by running model outside of climwin with climate data it found. 
pag1.ndvi.lin.m = glm(yvar~Bredb4 + climate,data=pag1.ndvi.lin.bestdata,family="binomial")
summary(pag1.ndvi.lin.m)

##Plot effects of the best model
library(MASS)
plot(Effect(focal.predictors = "Bredb4",pag1.ndvi.lin.m))
plot(Effect(focal.predictors = "climate",pag1.ndvi.lin.m))
detach("package:MASS") #interferes with dplyr



####Randomizations for NDVI

##Load randomization data from runs on cluster
ndvi.rand = read.csv(here::here("Output files","NDVI_mean_1yo_333_0_lin_randomizations_causes_paper_v4.csv"))

##Get p-value
pvalue(dataset=pag1.ndvi.lin.dataset,datasetrand=ndvi.rand,metric="AIC",sample.size=8)

##Plot randomizations
ggplot(data=ndvi.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.ndvi.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")



####Plot all delta AIC plots and randomizations####

library(ggpubr)

##Delta AIC plots 
A = plotdelta(dataset=pag1.rain.lin.dataset)
B = plotdelta(dataset=pag1.mintemp.lin.dataset)
C = plotdelta(dataset=pag1.maxtemp.lin.dataset)
D = plotdelta(dataset=pag1.ndvi.lin.dataset)
ggarrange(A,B,C,D)

##Randomization plots
A = ggplot(data=rain.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.rain.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")
B = ggplot(data=mintemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.mintemp.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")
C = ggplot(data=maxtemp.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.maxtemp.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")
D = ggplot(data=ndvi.rand,aes(x=deltaAICc)) + geom_histogram(color="black",fill="grey") + theme_cowplot() +
  geom_vline(aes(xintercept=pag1.ndvi.lin.dataset$deltaAICc[1]),color="red",size=1) + ylab("Frequency") + xlab("DeltaAICc")
ggarrange(A,B,C,D)

