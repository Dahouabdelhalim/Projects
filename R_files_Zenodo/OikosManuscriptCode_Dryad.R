setwd("/Users/RobertLaroche/Documents/Smallmouth_Bass_Project/Energetics_Manuscript/NewManuscriptSubmission_Oikos/DataForDryad/")
getwd()



#Loading required packages
library(lme4)
library(AICcmodavg)



#Loading required data sets

#Loading data set with the individual characteristics and degree days before first seasonal reproduction of all reproducing males who were found with eggs in their nest
ReproductionData <- read.csv("ReproductionData.csv")
#removing NA values from the dataset
ReproductionDataNAomit <- na.omit(ReproductionData)

#Loading data set with water temperatures across the years of the study
WaterTempsOver10 <- read.csv("WaterTempsOver10.csv")

#Loading data set with different measurements of days and degree days at important biological thresholds for smallmouth bass
TempMetrics <- read.csv("TempMetrics.csv")



#Figure 1 model creation

#Determining the days from the start of the season until the median reproductive date of each year
PeakReproDays <- c(median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 1999,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 1999,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2001,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2001,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2002,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2002,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2003,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2003,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2004,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2004,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2005,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2005,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2006,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2006,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2007,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2007,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2008,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2008,]$daysuntilfirsteggs != 0]), median(ReproductionDataNAomit[ReproductionDataNAomit$Years == 2009,]$daysuntilfirsteggs[ReproductionDataNAomit[ReproductionDataNAomit$Years == 2009,]$daysuntilfirsteggs != 0]))

#Determining the days each year from the start of the season to the first 15 degree day
DaysToFirst15 <- (TempMetrics$NumDaystoFirst10[c(1, 3:11)] + TempMetrics$NumDaystoFirst15[c(1, 3:11)])

##Calculating Response Time from first 15ºC day to median reproductive day
ResponseTime <- (PeakReproDays - DaysToFirst15)
ResponseTime

#Calculating Degree Days before the first day where temperatures averaged 15ºC
DDToFirst15 <- c(TempMetrics$DegDaystoFirst15[c(1, 3:11)])
DDToFirst15

#Linear model predicting reproductive response time using the degree days accumulated early in the season
Figure1Model <- lm(ResponseTime ~ DDToFirst15)
summary(Figure1Model)



#Figure 2, Table S2 model selection

#Creation of alternative models to explain timing of reproduction

#Modeling reproductive timing with body length
M1 <- lmer(log(degreedaysbeforeeggs) ~ log(length) + factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length and an interaction between body length and year
M2 <- lmer(log(degreedaysbeforeeggs) ~ log(length)*factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length and body length squared
M3 <- lmer(log(degreedaysbeforeeggs) ~ TLSquared + log(length) + factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length and body length squared, both with an interaction with year
M4 <- lmer(log(degreedaysbeforeeggs) ~ (TLSquared + log(length))*factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length, body length squared, and condition
M5 <- lmer(log(degreedaysbeforeeggs) ~ TLSquared + log(length) + condition + factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length, body length squared, and condition, with an interaction between the body length terms and year
M6 <- lmer(log(degreedaysbeforeeggs) ~ (TLSquared + log(length))*factor(Years) + condition + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with body length, body length squared, and condition, with an interaction between all terms and year
M7 <- lmer(log(degreedaysbeforeeggs) ~ (TLSquared + log(length) + condition)*factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition
M8 <- lmer(log(degreedaysbeforeeggs) ~ condition + factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition and an interaction between condition and year
M9 <- lmer(log(degreedaysbeforeeggs) ~ condition*factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition and body length
M10 <- lmer(log(degreedaysbeforeeggs) ~ condition + log(length) + factor(Years) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition and body length and an interaction between condition and year
M11 <- lmer(log(degreedaysbeforeeggs) ~ condition*factor(Years) + log(length) + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition and body length and an interaction between body length and year
M12 <- lmer(log(degreedaysbeforeeggs) ~ log(length)*factor(Years) + condition + (1|floy), data = ReproductionDataNAomit)

#Modeling reproductive timing with condition and body length and an interaction between all terms and year
M13 <- lmer(log(degreedaysbeforeeggs) ~ (condition + log(length))*factor(Years) + (1|floy), data = ReproductionDataNAomit)


#Using Akaike information criterion (AIC) to compare the success of alternative models in explaining our data
cand.set <- list(M1, M2, M3, M4, M5, M6, M7, M8, M9, M10, M11, M12, M13)
aictab(cand.set)



#Table S3. Individual Year Reproductive Timing Models

#Creation of alternative models for each year to validate choice of best model across years

#Modeling reproductive timing with body length
LM1 <- lm(log(degreedaysbeforeeggs) ~ log(length), data = subset(ReproductionDataNAomit, Years == 2009)) #2009 can be replaced with 1999, 2001-2009 to compare models for each year

#Modeling reproductive timing with condition
LM2 <- lm(log(degreedaysbeforeeggs) ~ condition, data = subset(ReproductionDataNAomit, Years == 2009))

#Modeling reproductive timing with body length and condition
LM3 <- lm(log(degreedaysbeforeeggs) ~ log(length) + condition, data = subset(ReproductionDataNAomit, Years == 2009))

#Modeling reproductive timing with body length and body length squared
LM4 <- lm(log(degreedaysbeforeeggs) ~ TLSquared + log(length), data = subset(ReproductionDataNAomit, Years == 2009))

#Modeling reproductive timing with body length, body length squared and condition
LM5 <- lm(log(degreedaysbeforeeggs) ~ TLSquared + log(length) + condition, data = subset(ReproductionDataNAomit, Years == 2009))


#Using Akaike information criterion (AIC) to compare the success of alternative individual year models in explaining our data
cand.set <- list(LM1, LM2, LM3, LM4, LM5)
aictab(cand.set)

#Inspecting the details of each model in addition to the confidence interval for the predictive terms
summary(LM1) #LM1 can be replaced with any of the above models LM1-LM5
confint(LM1)



#Table S4 post hoc models predicting size threshold for energetic constraints

#variables for post hoc tests of size threshold
DDbefore15 <- TempMetrics$DegDaystoFirst15[c(1, 3:11)]
WinterDuration <- c(NA, TempMetrics$NumDaystoFirst10[c(3:11)]+(365-(TempMetrics$NumDays[c(2:10)]+TempMetrics$NumDaystoFirst10[c(2:10)]))) ##adds the cold days of each year from jan to first 10º day to the days from the last year after the last 10º day before jan to measure winter length. No value for 1999 since it requires 1998 temps
PrevSeasonDuration <- TempMetrics$NumDays[c(NA, 2:10)]
PrevSeasonDD <- TempMetrics$DegreeDayList[c(NA, 2:10)]
avgTemp <- c(19.26559, 18.90806, 20.08394, 18.67492, 17.39657, 19.65527, 19.76456, 18.73405, 18.43912, 17.79765)

#length value (x-axis) of parabolic vertex for each year derived from model for Figure 2.
ThresholdDistribution <- c(42.16785, 41.40434, 40.01072, 41.56452, 42.7381, 41.68617, 41.33492, 41.98993, 44.05727, 41.41893)

#List of all variables tested as predictors of the size threshold
DDbefore15
WinterDuration
PrevSeasonDuration
PrevSeasonDD
avgTemp

#Models for each predictor

#DD before 15
DDbefore15Model <- lm(ThresholdDistribution ~ DDbefore15)
summary(DDbefore15Model)

#Winter Duration
WinterDurationModel <- lm(ThresholdDistribution ~ WinterDuration)
summary(WinterDurationModel)

#Previous Season Duration
PrevSeasonDurationModel <- lm(ThresholdDistribution ~ PrevSeasonDuration)
summary(PrevSeasonDurationModel)

#Previous Season Degree Days
PrevSeasonDDModel <- lm(ThresholdDistribution ~ PrevSeasonDD)
summary(PrevSeasonDDModel)

#Average Growth Season Temp
avgTempModel <- lm(ThresholdDistribution ~ avgTemp)
summary(avgTempModel)



#Figure S2 Model

#Creating a model predicting log transformed body length of individuals based on the log transformed degree days before they reproduced with year as a factor and individual identity (floy) as a random effect.
FigS2Model <- lmer(log(length) ~ log(degreedaysbeforeeggs) + factor(Years) + (1|floy), data = ReproductionDataNAomit)
#Looking at a summary of the results of this model
summary(FigS2Model)

