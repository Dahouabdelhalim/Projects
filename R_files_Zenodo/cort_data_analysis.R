################################################################################
##########EABL Anthro Noise Project - Cort and Telomere Analysis##########
################################################################################
library(lme4) #creates linear mixed models
library(car) #Anova test
library(ggplot2) #create good graphs
library(dplyr) #data management
library(lubridate) #time and date variable management
library(bbmle) #AICc tests
library(lmerTest) #Calculates p-values for LMMs

# library(lattice)


# library(agricolae)
# library(fitdistrplus)
# library(goft)
# library(mosaic)
# library(psych) #used to calculate standard error
# library(dplyr)
# library(lubridate)
# library(multcomp)
# library(lsmeans)

# Loading and Cleaning Data -----------------------------------------------

#Set working directory and import dataset
setwd()

hy_cort = read.csv("nestling_cort_data_pub.csv", header = TRUE)
ahy_cort = read.csv("adult_cort_data_pub.csv", header = TRUE)

# Building Models ---------------------------------------------------------

###Testing best random effects
m1rf = lmer(log_cort ~ treatment*jdate + #best random effect
              treatment*brood_size + 
              (1|box_num), 
           REML = TRUE, 
           data = hy_cort)

m2rf = lmer(log_cort ~ treatment*jdate + 
              treatment*brood_size + 
              (1|site), 
            REML = TRUE, 
            data = hy_cort)

m3rf = lmer(log_cort ~ treatment*jdate + 
              treatment*brood_size + 
              (1|box_num)+(1|site), 
            REML = TRUE, 
            data = hy_cort)

#Nestling AIC Comparisons with box number (site + box) as random effect
m1n = lmer(log_cort ~ treatment*jdate + 
             treatment*brood_size + 
             as.factor(sex) +
              (1|box_num), 
            REML = FALSE, 
            data = hy_cort)
summary(m1n) 


#Adult AIC model comparisons
# Building Models ---------------------------------------------------------
m1a = lm(log_cort ~ treatment*jdate + 
             treatment*brood_size, 
           data = ahy_cort)
summary(m1a) 
anova(m1a)

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


# Plotting Graphs ---------------------------------------------------------

library(ggsci)#package for scientific journal color palettes
library(ggplot2)
library(ggthemes)
hy_means = hy_cort %>%
  group_by(treatment) %>%
  dplyr::summarize(n = n(),
                   meanCort = mean(cort_levels), 
                   seCort = (sd(cort_levels)/sqrt(n))
  )
#Scatterplot showing Cort levels for each treatment group across julian date
ggplot(data = hy_cort, aes(x=jdate, 
                           y=cort_levels, 
                           group=treatment, 
                           color = treatment)) +
  geom_point()+
  geom_line()+
  theme_classic()

#Scatterplot showing cort levels between treatmetn groups across brood size
ggplot(data = hy_cort, aes(x=brood_size, 
                           y=cort_levels, 
                           group=treatment, 
                           color = treatment)) +
  geom_point()+
  geom_line()+
  theme_classic()

#Bar graph showing cort level means across treatment group
ggplot(data = hy_means, aes(x=treatment, y=meanCort, group=treatment, fill = treatment)) +
  geom_bar(stat="identity", 
           position=position_dodge())+
  geom_errorbar(aes(ymin=meanCort-seCort,
                    ymax=meanCort+seCort),
                width=0.25,
                position=position_dodge(0)) +
  theme_classic() +
  fill_palette(palette = "grey")
# 