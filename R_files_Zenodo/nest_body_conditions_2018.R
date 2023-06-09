##############################################################
################ EABL Nestling Growth Analysis################
##############################################################

####Original .Rmd file created by Dr. Darren S. Proppe on 9/21/2018
###Modified into basic .R script by Meelyn M. Pandit on 2/11/2021


library(lme4) #creates linear mixed models
library(car) #Anova test
library(ggplot2) #create good graphs
library(dplyr) #data management
library(lubridate) #time and date variable management
library(bbmle) #AICc tests
library(lmerTest) #Calculates p-values for LMMs

setwd() ###Set Working Directory

growth = read.csv("nestling_body_conditions_pub.csv", header = TRUE)

###Calculate body conditions - residuals from wing/mass regression
# bc <- lmer(mass~wing, growth)
# resid(bc)
# growth$body <- bc$resid

###body condition - finding best LMM
g.mod <- lmer(body ~ treatment*day + 
                logday + 
                treatment*logday + 
                julian + 
                brood + 
                (1|box) + (1|band), 
              growth)

g.mod.day <- lmer(body ~ treatment*day + 
                    julian + 
                    brood + 
                    (1|box) + (1|band), 
                  growth)

g.mod.logday <- lmer(body ~ logday + 
                       treatment*logday + 
                       julian + 
                       brood + 
                       (1|box) + (1|band), 
                     growth)

g.null = lmer(body ~ 1 + 
                (1|box) + (1|band), 
              growth)

AICctab(g.mod,g.mod.day,g.mod.logday,g.null, nobs = 282, base=TRUE,delta=TRUE, sort=TRUE, weights=TRUE) #AICc test to determine best fitting model

Anova(g.mod) #gives p value
summary(g.mod) #gives estimate and error (control is constant)

# Regression Diagnostics..... ---------------------------------------------

library(car)
library(predictmeans)
library(tidyverse)
library(broom)

###Outlier Test
outlierTest(m1n)
qqPlot(m1n, main="QQ Plot") #qq plot for studentized resid
leveragePlots(m1n) # leverage plots

r<-residuals(m1n)
ft<-fitted(m1n)
par(mfrow=c(2,2))
library(MASS)
truehist(r,main="Histogram of Residuals",xlab="Residuals")
curve(dnorm(x,mean=mean(r),sd=sd(r)),add=TRUE)
qqnorm(r, ylim=range(r), main="QQNorm Plot",ylab="Quantiles of (ie, ordered) residuals", xlab="Quantiles of normal distribution")
qqline(r,lty=2)
plot(r~ft,main="Residuals vs Fitted",xlab="Fitted (predicted) values",ylab="Residuals");abline(h=0)
qqnorm(ft, ylim=range(ft), main="QQNorm Plot",ylab="Quantiles of fitted values", xlab="Quantiles of normal distribution")
qqline(ft,lty=2)


