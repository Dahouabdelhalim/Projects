################################################################################
##########EABL Anthro Noise Project - Telomere Analysis##########
################################################################################

library(lme4) #creates linear mixed models
library(car) #Anova test
library(ggplot2) #create good graphs
library(dplyr) #data management
library(lubridate) #time and date variable management
library(bbmle) #AICc tests
library(lmerTest) #Calculates p-values for LMMs

setwd() #Set working directory

# # Loading Data
# setwd("C:/Users/meely/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/data_clean")
# data = read.csv("mi_banding_data_2018_final.csv", header = TRUE)
# head(data)
# 
# #Dividing data by age
# telomere = data %>%
#   filter(data$telomere_lengths>0)
# 
# telomere_adult = telomere %>% filter(telomere$age == "asy" | telomere$age == "sy")
# telomere_nest = telomere %>% filter(telomere$age == "hy")
# 
# #Renaming column 1
# names(telomere_nest)[1]<-"date"
# names(telomere_adult)[1] = "date"
# 
# names(telomere_nest)[6] = "band_num"
# names(telomere_adult)[6] = "date"
# 
# names(telomere_nest)[16] = "brood_size"
# telomere_nest$brood_size = c(4,4,4,4,5,5,5,5,5,5,2,4,4,4,4,4,4,3,3,3,5,5,5,5,5,5,5,4,4,4,3,3,4,4,2,2,3,3,3,5,3,3,3,3,3,4)
# telomere_nest$box_num = paste(telomere_nest$location, telomere_nest$box)
# #telomere$telomere_lengths <- as.numeric(as.character(telomere$lengths)) #telomere concentration (ng/mL) was not numeric so converting factor to numeric
# #telomere$box_num = paste(telomere$site,telomere$box)
# 
# # telomere %>%
# #   select(box_num,band,date,treatment,age,sex,telomere_levels) telomere %>%
# #   filter(treatment != "internal")
# # attach(telomere)
# 
# #Checking Distributions
# hist(telomere_nest$telomere_lengths)
# shapiro.test(telomere_nest$telomere_lengths)
# 
# hist(log(telomere_nest$telomere_levels))
# shapiro.test(log(telomere_nest$telomere_levels))
# View(log(telomere_nest$telomere_levels))
# 
# setwd("~/OneDrive - University of Oklahoma/University of Oklahoma/Bridge Lab/EABL Anthro Noise/data_clean") #Saving cleaned dataframe into data_clean directory
# write.csv(telomere_nest, "nestling_telomere_data_clean.csv")
# write.csv(telomere_adult, "adult_telomere_data_clean.csv")


# Loading Data ------------------------------------------------------------

telomere_nest = read.csv("nestling_telomere_data_pub.csv", header = TRUE)
telomere_nest$Date = mdy(telomere_nest$date) #convert date vector into format R recognizes
telomere_nest$julian_date = as.numeric(format(telomere_nest$Date, "%j")) #Convert new date into julian date


# Building Models ---------------------------------------------------------

control=lmerControl(optCtrl=list(maxfun=2e5,
                                 optimizer = "bobyqa")) #Control to help with singular fit.

m1t = lmer(log(telomere_lengths) ~ treatment*julian_date + 
             brood_size + 
             (1|box_num), 
           REML = FALSE, 
           control = control,
           data = telomere_nest)
summary(m1t)

# #Regression Diagnostics -------------------------------------------------

library(car)
library(predictmeans)
library(tidyverse)
library(broom)

r<-residuals(m10)
ft<-fitted(m10)
par(mfrow=c(2,2))
library(MASS)
truehist(r,main="Histogram of Residuals",xlab="Residuals")
curve(dnorm(x,mean=mean(r),sd=sd(r)),add=TRUE)
qqnorm(r, ylim=range(r), main="QQNorm Plot",ylab="Quantiles of (ie, ordered) residuals", xlab="Quantiles of normal distribution")
qqline(r,lty=2)
plot(r~ft,main="Residuals vs Fitted",xlab="Fitted (predicted) values",ylab="Residuals");abline(h=0)
qqnorm(ft, ylim=range(ft), main="QQNorm Plot",ylab="Quantiles of fitted values", xlab="Quantiles of normal distribution")
qqline(ft,lty=2)
