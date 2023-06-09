##############################
###                        ###
###  Stay-at-home Fatigue  ###
###       09/11/2020       ###
###      J. Shearston      ###
###                        ###
##############################

# Last Updated: February 19, 2021

# Summary
# This script contains code for doing a time series analysis of traffic data
# from Google traffic maps. It loads data that has already been cleaned and
# prepared through image analysis in MatLab, by Markus Hilpert. It includes descriptive
# analysis, modeling with GAMs, seasonal trend decomposition, and the creation
# of several figures for a manuscript. Multiple exploratory GAM models are
# constructed before a final set of models is selected. We attempt to fit GAMs 
# that most appropriately fit the observed time series. 

# Table of Contents
# 1: Preparation
# 2: Load data and add variables
# 3: Review time series
# 4: Exploratory models for each intervention level (separate models for each level)
# 5: Exploratory model with intervention levels as one categorical variable (one full model)
# 6: Determine optimal df values for full model with complete time series
# 7: Construct nested models with full time series to determine which variables are needed
# 8: AIC and Likelihood Ratio tests for comparing full models
# 9: Evaluate residual plots for model on full time series
# 10: Determine breakpoints for three models, using time series to drive cutpoints
# 11: Three model solution
# 12: Four model solution
# 13: Predictions for "typical" hours of day, by weekend and weekday
# 14: Seasonal trend decomposition (STL)
# 15: Figures for manuscript
# 16: Table 1 for manuscript (Supplementary table 1 uses data collated separately)


##########################################################################
# 1: Preparation
##########################################################################

# 1A Load Packages

library(tidyverse)
library(lubridate)
library(mgcv)
library(visreg)
library(lmtest)
library(tidymv)
library(grid)
library(gridExtra)
require(stats)


##########################################################################
# 2: Load Data & Add Vars
##########################################################################

# Note: Data is in a txt file, with the DateNum variable representing Matlab's
# datetime structure. However, a character string of datetime is also provided

# 2A Read in data
gcc <- read_csv("./CCC_COVID19.txt") %>% janitor::clean_names()

# 2B Reformat date variable and create additional variables
gcc <- gcc %>% 
  mutate(date_time = as.POSIXct(date_string,format=c('%d-%b-%Y %H:%M:%S')),
         hour_of_day = hour(date_time),
         dow = wday(date_time, week_start = 1),                  # 1=Monday; 7=Sunday
         weeknum = week(date_time),
         daynum = yday(date_time),
         intervention = case_when(
           daynum < 82 ~ "Pre-Intervention",
           daynum >= 82 & daynum < 160 ~ "PAUSE",                # March 22 - June 7
           daynum >= 160 & daynum < 174 ~ "Phase 1",             # June 8 - June 21
           daynum >= 174 & daynum < 188 ~ "Phase 2",             # June 22 - July 5
           daynum >= 188 & daynum < 199 ~ "Phase 3",             # July 6 - July 16
           daynum >= 199 ~ "Phase 4"                             # July 16 onward
         ),
         intervention = as_factor(intervention))

# 2C Check levels of new "intervention" factor
levels(gcc$intervention)
table(gcc$intervention, useNA = "always")

# 2D Converting day of week variable (dow) to a factor
gcc2 <- gcc %>% 
  mutate(dow = as_factor(as.character(dow)))

# 2E Add new weekday/weekend dummy vars to dataframe
gcc2 <- gcc2 %>% 
  mutate(weekday = ifelse(dow==1 | dow==2 | dow==3 | dow==4 | dow==5, 1, 0),
         weekend = ifelse(dow==6 | dow ==7, 1, 0),
         weekday_hourofday = weekday*hour_of_day,
         weekend_hourofday = weekend*hour_of_day)

# 2F Restrict to only 2020
gcc2 <- gcc2 %>% filter(date_num_matlab < 738157.0)
# Note: Should have n=2928 observations, or 8 observations/day * 366 days (2020 is a leap year)


##########################################################################
# 3: Review time series
##########################################################################

# Note: We would like to review each color's time series to identify which 
# color might be best to complete analysis on. We choose red because of its
# clear pandemic related signal (decrease in April) and easier interpretability
# e.g., a decrease in percent roads with a red color is a decrease in traffic
# congestion 


# 3A View time series of all colors
gcc2 %>% 
  pivot_longer(maroon:gray, names_to = "gcc_color", values_to = "percent_roads") %>% 
  ggplot(aes(x = date_time, y = percent_roads)) +
  geom_line() +
  facet_wrap(~ gcc_color)

# 3B View time series of red color
gcc2 %>% 
  pivot_longer(maroon:gray, names_to = "gcc_color", values_to = "percent_roads") %>% 
  filter(gcc_color == "red") %>% 
  ggplot(aes(x = date_time, y = percent_roads)) +
  geom_line(color = "red") +
  ylab("Percent Area w Red-Colored Roads") + xlab("Date")


##########################################################################
# 4. Exploratory Models for Each Intervention Level
##########################################################################

# Note: In this series of models, GAMs are fit for periods of time defined
# by policies implemented in response to the pandemic

# 4A Red: pre-intervention

pre_i <- gcc %>% filter(intervention == "Pre-Intervention")

plot(y=pre_i$red, x=pre_i$date_time, col = "red", type = "l")

model_red_prei <- gam(red ~ s(weeknum, k=3) + s(hour_of_day, k=7) + dow,
                      data = pre_i)

plot(model_red_prei)

summary(model_red_prei)

# 4B Red: PAUSE

pause <- gcc %>% filter(intervention == "PAUSE")

plot(y=pause$red, x=pause$date_time, col = "red", type = "l")

model_red_pause <- gam(red ~ s(weeknum, k=4) + s(hour_of_day, k=7) + dow,
                       data = pause)

plot(model_red_pause)

summary(model_red_pause)

# 4C Red: Phase 1

p1 <- gcc %>% filter(intervention == "Phase 1")

plot(y=p1$red, x=p1$date_time, col = "red", type = "l")

model_red_p1 <- gam(red ~ s(weeknum, k=2) + s(hour_of_day, k=7) + dow,
                    data = p1)

plot(model_red_p1)

summary(model_red_p1)

# 4D Red: Phase 2

p2 <- gcc %>% filter(intervention == "Phase 2")

plot(y=p2$red, x=p2$date_time, col = "red", type = "l")

model_red_p2 <- gam(red ~ s(weeknum, k=2) + s(hour_of_day, k=7) + dow,
                    data = p2)

plot(model_red_p2)

summary(model_red_p2)

# 4E Red: Phase 3

p3 <- gcc %>% filter(intervention == "Phase 3")

plot(y=p3$red, x=p3$date_time, col = "red", type = "l")

model_red_p3 <- gam(red ~ s(weeknum, k=3) + s(hour_of_day, k=7) + dow,
                    data = p3)

plot(model_red_p3)

summary(model_red_p3)

# 4F Red: Phase 4

p4 <- gcc %>% filter(intervention == "Phase 4")

plot(y=p4$red, x=p4$date_time, col = "red", type = "l")

model_red_p4 <- gam(red ~ s(weeknum, k=8) + s(hour_of_day, k=7) + dow,
                    data = p4)

plot(model_red_p4)

summary(model_red_p4)


##########################################################################
# 5: Exploratory Model w Categorical Intervention
##########################################################################

# Notes: After running models for each intervention level separately, we
# wanted to see if the time series could be adequately modeled in one full
# model as well. 

# 5A Red: Full time-series

plot(y=gcc$red, x=gcc$date_time, col = "red", type = "l")

model_red_full <- gam(red ~ s(weeknum, k=7) + s(hour_of_day, k=7) + dow + intervention, data = gcc)

plot(model_red_full)

summary(model_red_full)


##########################################################################
# 6: Determine optimal df for weeknum & hour_of_day (full time series)
##########################################################################

# Notes: Use GCV to determine optimal k value for full time series model; lower GCV is better
# Tried k=10 and k=9 but received error "a term has fewer unique covariate combinations
# than specified maximum degrees of freedom" b/c 8 is the largest available given the 
# number of categories hour_of_day

test <- gam(red ~ s(weeknum) + s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(test)
plot(test)
gam.check(test)

mod4gcv1 <- gam(red ~ s(weeknum, k=8) + s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(mod4gcv1) # gcv = 0.040

mod4gcv2 <- gam(red ~ s(weeknum, k=7) + s(hour_of_day, k=7) + dow + intervention, data = gcc)
summary(mod4gcv2) # gcv = 0.045

mod4gcv3 <- gam(red ~ s(weeknum, k=6) + s(hour_of_day, k=6) + dow + intervention, data = gcc)
summary(mod4gcv3) # gcv = 0.045

mod4gcv4 <- gam(red ~ s(weeknum, k=6) + s(hour_of_day, k=7) + dow + intervention, data = gcc)
summary(mod4gcv4) # gcv = 0.045

mod4gcv5 <- gam(red ~ s(weeknum, k=6) + s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(mod4gcv5) # gcv = 0.042

mod4gcv6 <- gam(red ~ s(weeknum, k=7) + s(hour_of_day, k=6) + dow + intervention, data = gcc)
summary(mod4gcv6) # gcv = 0.045

mod4gcv7 <- gam(red ~ s(weeknum, k=7) + s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(mod4gcv7) # gcv = 0.042

mod4gcv8 <- gam(red ~ s(weeknum, k=8) + s(hour_of_day, k=6) + dow + intervention, data = gcc)
summary(mod4gcv8) # gcv = 0.045

mod4gcv9 <- gam(red ~ s(weeknum, k=8) + s(hour_of_day, k=7) + dow + intervention, data = gcc)
summary(mod4gcv9) # gcv = 0.045

# Optimal values for full time series model
# weeknum: k=8
# hour_of_day: k=8


##########################################################################
# 7: Run Nested Models (full time series)
##########################################################################

# Notes: Run nested models with different combinations of variables, and then
# do comparisons between the models using AIC and Likelihood ratio tests
# to determine which models best fit the data. We are unsure whether all of
# the variables (hour of day, day of week, intervention, and week number) are
# needed in the models to best fit the time series data. In this chunk, we 
# create a model for every possible combination of variables, which we will
# then compare in the next chunk. For these models, the full time series is used.

# 7A All vars (hour_of_day, dow, intervention, weeknum)

nestmod_full <- gam(red ~ s(weeknum) + s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(nestmod_full) #gcv = 0.040

# 7B Three vars A (weeknum, dow, intervention) (various k checked to get lowest gcv)

nestmod_3a <- gam(red ~ s(weeknum, k=8) + dow + intervention, data = gcc)
summary(nestmod_3a)
#k=7, gcv=0.120
#k=8, gcv=0.119 (selected)
#k=9, gcv=0.119
#k=10, gcv=0.119
#k=11, gcv=0.119

# 7C Three vars B (hour_of_day, dow, intervention) (various k checked to get lowest gcv)

nestmod_3b <- gam(red ~ s(hour_of_day, k=8) + dow + intervention, data = gcc)
summary(nestmod_3b)
#k=7, gcv=0.048
#k=8, gcv=0.045 (selected)
#k=9, gcv=0.045
#k=10, gcv=0.045

# 7D Three vars C (weeknum, hour_of_day, intervention) (various k checked to get lowest gcv)

nestmod_3c <- gam(red ~ s(weeknum, k=9) + s(hour_of_day, k=8) + intervention, data = gcc)
summary(nestmod_3c)
#k=7,8; gcv=0.042
#k=8,8; gcv=0.042
#k=9,8; gcv=0.041 (selected)
#higher k values cause error

# 7E Three vars D (weeknum, hour_of_day, dow) (various k checked to get lowest gcv)

nestmod_3d <- gam(red ~ s(weeknum, k=8) + dow + s(hour_of_day, k=8),
                  data = gcc)
summary(nestmod_3d)
#k=7,8; gcv=0.042
#k=8,8; gcv=0.040 (selected)
#k=9,8; gcv=0.040
#k=10,8; gcv=0.040
#k=11,8; gcv=0.040

# 7F Two vars A (intervention, dow)

nestmod_2a <- gam(red ~ dow + intervention, data = gcc)
summary(nestmod_2a)

# 7G Two vars B (intervention, hour_of_day) (8 is max k for hour_of_day)

nestmod_2b <- gam(red ~ s(hour_of_day, k=8) + intervention, data = gcc)
summary(nestmod_2b)

# 7H Two vars C (intervention, weeknum) (various k checked to get lowest gcv)

nestmod_2c <- gam(red ~ s(weeknum, k=8) + intervention, data = gcc)
summary(nestmod_2c)
#k=7, gcv=0.121
#k=8, gcv=0.121 (selected)
#k=9, gcv=0.121
#k=10, gcv=0.121

# 7I Two vars D (dow, hour_of_day) (8 is max k for hour_of_day)

nestmod_2d <- gam(red ~ s(hour_of_day, k=8) + dow, data = gcc)
summary(nestmod_2d)

# 7J Two vars E (weeknum, dow) (various k checked to get lowest gcv)

nestmod_2e <- gam(red ~ s(weeknum, k=9) + dow, data = gcc)
summary(nestmod_2e)
#k=7, gcv=0.120
#k=8, gcv=0.119
#k=9, gcv=0.118 (selected)
#k=10, gcv=0.118
#k=11, gcv=0.118

# 7K Two vars F (weeknum, hour_of_day) (various k checked to get lowest gcv)

nestmod_2f <- gam(red ~ s(hour_of_day, k=8) + s(weeknum, k=8), data = gcc)
summary(nestmod_2f)
#k=8,7; gcv=0.044
#k=8,8; gcv=0.042
#higher k values resulted in errors


##########################################################################
# 8: Compare Nested Models (full time series)
##########################################################################

# Notes: Now we compare the nested models created in Section 7 using a 
# combination of AIC (to compare all models) and likelihood ratio tests
# (for nested models)

# 8A AIC to compare all models
AIC(nestmod_2a, nestmod_2b, nestmod_2c, nestmod_2d, nestmod_2e,
    nestmod_2f, nestmod_3a, nestmod_3b, nestmod_3c, nestmod_3d,
    nestmod_full)
# smallest (best) AIC is nestmod_3d, but this doesn't have the intervention var
# next smallest (best) is nestmod_full, then nestmod_3c

# 8B Likelihood ratio tests for nested models

# 8B Set 1: weeknum
lrtest(nestmod_2c, nestmod_3a) #3a is better
lrtest(nestmod_2c, nestmod_3c) #3c is better
AIC(nestmod_3a, nestmod_3c) #3c is better
lrtest(nestmod_3c, nestmod_full) #full model is better

# 8B Set 2: hour_of_day
lrtest(nestmod_2b, nestmod_3b) #3b is better
lrtest(nestmod_2b, nestmod_3c) #3c is better
AIC(nestmod_3b, nestmod_3c) #3c is better

# 8B set 3: dow
lrtest(nestmod_2a, nestmod_3a) #3a is better
lrtest(nestmod_2a, nestmod_3b) #3b is better
AIC(nestmod_3a, nestmod_3b) #3b is better

# Notes: full_model is best, then 3c. Will move forward with the full model (nestmod_full)


##########################################################################
# 9: Evaluate residual plots: full model (nestmod_full)
##########################################################################

# Notes: Now that a model has been chosen, we will evaluate residual plots
# and visually inspect how well the model is matching the data

# 9A Review model results
summary(nestmod_full)
#str(nestmod_full) # find all names/items in gam object

# 9B Actual vs predicted plot with slope=1 line
plot(nestmod_full$y,nestmod_full$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red") 
# nope, not good --> points are curving 
# away from slope line at higher y values

# 9C Predicted vs residuals
plot(nestmod_full$fitted.values, nestmod_full$residuals) 
# not looking good --> can see a u-shaped pattern

# 9D Time series of actual values with predictions overlaid to visualize fit
plot(nestmod_full$y, type = "l")
lines(nestmod_full$fitted.values, col = "blue") 
# This is not a good match --> the model does not capture the high peaks at all

# 9E Replicate full model with day of week variable converted to factor(gcc2)
nestmod_full2 <- gam(red ~ s(weeknum, k=7) + s(hour_of_day, k=8) + dow + intervention, data = gcc2)
summary(nestmod_full2)    
plot(nestmod_full2)
plot(nestmod_full2$y,nestmod_full2$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(nestmod_full2$y, type = "l")
lines(nestmod_full2$fitted.values, col = "blue") # not good match
AIC(nestmod_full, nestmod_full2)    #nestmod_full2 with dow as factor is better per AIC

# Notes: 
# This model is not fitting the data well. Predictions underestimate the
# hour of day peaks during the pre-pandemic period and the return to "normal"
# while also overestimating them during the PAUSE period. Should break the 
# data into three chunks and run three separate models. Unfortunately it doesn't
# appear that a single model can fit the data well. 


##########################################################################
# 10: Determine breakpoints for three models
##########################################################################

# Note: In this series of models, GAMs are fit for periods of time defined
# by changes in the time series itself, in order to construct three models
# that best fit the data. First, we attempt to use the model residuals to 
# determine where the model begins to fail to fit the data. We create avg
# daily residuals to remove within-day variability and see changes more clearly

# 10A View residuals from full model
plot(nestmod_full2$residuals, type = "l")

# 10B Create tibble of resids and data
brk_points <- tibble(gcc2 %>% filter(!is.na(red)) %>% dplyr::select(red),
                     resids = nestmod_full2$residuals,
                     gcc2 %>% filter(!is.na(red)) %>% dplyr::select(date_time),
                     date = as.Date(date_time, "EST"))

# 10C Create avg daily residuals
Dbrk_points <- brk_points %>% 
  group_by(date) %>% 
  summarise(mean_resid = mean(resids),
            sd_resid = sd(resids),
            mean_red = mean(red))

# 10D View mean residuals 
plot(Dbrk_points$date, Dbrk_points$mean_resid, type = "l")
abline(h=0, col="blue", lty=2)
# March 14 - cross into negative for the steep drop
# However, it is still pretty messy and hard to see where we might make 
# second break point.

# 10 E Create avg weekly residuals
# This might remove some of the periodicity and make the cutpoint easier to identify
Wbrk_points <- brk_points %>% 
  mutate(weeknum = week(date_time)) %>% 
  group_by(weeknum) %>% 
  summarise(mean_resid = mean(resids),
            sd_resid = sd(resids),
            mean_red = mean(red))

# 10F View mean residuals 
plot(Wbrk_points$weeknum, Wbrk_points$mean_resid, type = "l")
abline(h=0, col="blue", lty=2)
# Week 15 becomes positive again (April 8) --> April 11
# Week 16 is peak of first positive return (April 15)
# Week 19 is first negative after first positive return (May 6)

# After running three separate models based on the above dates (not shown), I actually think 
# it's better to use the time series of actual y values rather than residuals 
# to determine the breakpoints.

# 10G Inspect time series for first breakpoint
plot(gcc2$red, type = "l", xlim=c(0,500))
# rownum~350 is good first breakpoint, March 14/15

plot(gcc2$red, type = "l")

# 10H Inspect time series for second breakpoint
plot(gcc2$red, type = "l", xlim=c(300,1400))
# rownum~1050 is good second breakpoint, June 11


##########################################################################
# 11: Three Model Solution
##########################################################################

# Notes: Here we try using three models to fit the data. Because we are no
# longer using the full time series in each model, we again need to check which 
# variables are needed. Because we are not using interventions to break the 
# data but are breaking it according to dramatic shifts in the trends (which do
# at least in part depend on the interventions) we remove the intervention variable

# 11A Pre-intervention period (no weeknum)
gcc2_pre <- gcc2 %>% filter(date_time < "2020-03-14 00:00:00")
pre_trimodel <- gam(red ~ s(hour_of_day, k=8) + dow, data = gcc2_pre)
summary(pre_trimodel)    #gcv lowest w k=8
plot(pre_trimodel$y,pre_trimodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(pre_trimodel$y, type = "l")
lines(pre_trimodel$fitted.values, col = "blue") #better fit, but still doesn't capture high peaks

# 11B Pre-intervention period (with weeknum)
pre_trimodel2 <- gam(red ~ s(hour_of_day, k=8) + s(weeknum, k=5) + dow, data = gcc2_pre)
summary(pre_trimodel2)    #gcv lowest w k=8 and k=5
plot(pre_trimodel2$y,pre_trimodel2$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(pre_trimodel2$y, type = "l")
lines(pre_trimodel2$fitted.values, col = "blue") #still doesn't capture high peaks

# 11C Compare AIC pre-intervention period models with and without weeknum
AIC(nestmod_full2, pre_trimodel, pre_trimodel2)    #nestmod_full2<pre_trimodel2<pre_trimodel

# 11D Covid-period 1 (include weeknum since that model fit better for pre-intervention period)
gcc2_cov1 <- gcc2 %>% filter(date_time > "2020-03-14 00:00:00" & date_time < "2020-06-11 00:00:00")
cov1_trimodel <- gam(red ~ s(hour_of_day, k=8) + s(weeknum, k=8) + dow, data = gcc2_cov1)
summary(cov1_trimodel)   
plot(cov1_trimodel$y,cov1_trimodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov1_trimodel$y, type = "l")
lines(cov1_trimodel$fitted.values, col = "blue") #decent fit, doesn't capture high peaks toward end

# 11E Covid-period 2
gcc2_cov2 <- gcc2 %>% filter(date_time > "2020-06-11 00:00:00")
cov2_trimodel <- gam(red ~ s(hour_of_day, k=8) + s(weeknum, k=8) + dow, data = gcc2_cov2)
summary(cov2_trimodel)   
plot(cov2_trimodel$y,cov2_trimodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov2_trimodel$y, type = "l")
lines(cov2_trimodel$fitted.values, col = "blue") # better fit, but doesn't capture highest peaks

# 11F Covid-period 2 with different spline specification
# Use a cyclic penalized cubic regression spline smooth for hour of day
cov2.2_trimodel <- gam(red ~ s(hour_of_day, bs="cc", k=8) + s(weeknum, k=8) + dow, data = gcc2_cov2)
summary(cov2.2_trimodel)   
plot(cov2.2_trimodel$y,cov2.2_trimodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov2.2_trimodel$y, type = "l")
lines(cov2.2_trimodel$fitted.values, col = "blue")
# Basically this appears identical to the model in 11E

# 11G Compare Period 2 models
AIC(cov2_trimodel, cov2.2_trimodel) #cov2 and cov2.2 identical 

# Notes: The 3 model solution is still not capturing the peaks very well. Let's
# try a four model solution. 

##########################################################################
# 12: Four Model Solution
##########################################################################

# Notes: Changes for this set of models include:
# Four models rather than 3 to better capture the variability in the lockdown phase
# No weeknum variable in any models (because we are breaking up the time periods)
# Separate weekday/weekend splines rather than one dow categorical var
#    to better capture the high weekday peaks that are not being adequately fit

# 12A Determine breakpoints for four time periods by visually inspecting time series
plot(gcc2$red, type = "l")
plot(gcc2$red, type = "l", xlim=c(0,500)) # rownum~350 is good first breakpoint, March 14/15
plot(gcc2$red, type = "l", xlim=c(300,1000)) # rownum~800 is good second breakpoint, May 10
plot(gcc2$red, type = "l", xlim=c(1000, 1500)) # rownum~1100 is good third breakpoint, June 17

# 12B Pre-intervention period
gcc2_pre <- gcc2 %>% filter(date_time < "2020-03-14 00:00:00")
pre_quadmodel <- gam(red ~ s(weekday_hourofday, k=8, bs="cc") +
                       s(weekend_hourofday, k=8, bs="cc"), data = gcc2_pre)
summary(pre_quadmodel)
plot(pre_quadmodel$y,pre_quadmodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(pre_quadmodel$y, type = "l")
lines(pre_quadmodel$fitted.values, col = "blue")
mean(gcc2_pre$red, na.rm=T)
# Notes: This model has a better fit, but still doesn't capture the high peaks

# 12C Covid-period 1
gcc2_cov1 <- gcc2 %>% filter(date_time >= "2020-03-14 00:00:00" & date_time < "2020-05-20 00:00:00")
cov1_quadmodel <- gam(red ~ s(weekday_hourofday, k=8, bs="cc") +
                       s(weekend_hourofday, k=8, bs="cc"), data = gcc2_cov1)
summary(cov1_quadmodel)
plot(cov1_quadmodel$y,cov1_quadmodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov1_quadmodel$y, type = "l")
lines(cov1_quadmodel$fitted.values, col = "blue") 
mean(gcc2_cov1$red, na.rm=T)
# Notes: Good in the middle but not catching peaks at beginning and end

# 12D Covid-period 2
gcc2_cov2 <- gcc2 %>% filter(date_time >= "2020-05-20 00:00:00" & date_time < "2020-06-17 00:00:00")
cov2_quadmodel <- gam(red ~ s(weekday_hourofday, k=8, bs="cc") +
                        s(weekend_hourofday, k=8, bs="cc"), data = gcc2_cov2)
summary(cov2_quadmodel)
plot(cov2_quadmodel$y,cov2_quadmodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov2_quadmodel$y, type = "l")
lines(cov2_quadmodel$fitted.values, col = "blue")
mean(gcc2_cov2$red, na.rm=T)
# Notes: Pretty good job fitting the data!

# 12E Covid-period 3
gcc2_cov3 <- gcc2 %>% filter(date_time >= "2020-06-17 00:00:00")
cov3_quadmodel <- gam(red ~ s(weekday_hourofday, k=8, bs="cc") +
                        s(weekend_hourofday, k=8, bs="cc"), data = gcc2_cov3)
summary(cov3_quadmodel)
plot(cov3_quadmodel$y,cov3_quadmodel$fitted.values, xlim=c(0,3), ylim=c(0,3))
abline(a=0, b=1, lty=2, col="red")
plot(cov3_quadmodel$y, type = "l")
lines(cov3_quadmodel$fitted.values, col = "blue")
mean(gcc2_cov3$red, na.rm=T)
# Notes: Fairly decent job!

# We will use this four model solution for the manuscript, as it best fits the data
# of everything we've tried so far. 


##########################################################################
# 13: Predictions (four models) and Statistical Comparisons
##########################################################################

# 13A Create new data to predict a "typical week" for each time period
predicts <- tibble(                                                  # Create new data
  hour_of_day = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
                  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
  weekend = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  weekday = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
              1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
  weekday_hourofday = weekday*hour_of_day,
  weekend_hourofday = weekend*hour_of_day)

# 13B Run predictions and SEs
predicts <- predicts %>%                                             # Run predictions and SEs
  mutate(
    precov_pred = unlist(predict(pre_quadmodel, predicts, se.fit=TRUE)[1]), 
    precov_se = unlist(predict(pre_quadmodel, predicts, se.fit=TRUE)[2]),
    cov1_pred = unlist(predict(cov1_quadmodel, predicts, se.fit=TRUE)[1]),
    cov1_se = unlist(predict(cov1_quadmodel, predicts, se.fit=TRUE)[2]),
    cov2_pred = unlist(predict(cov2_quadmodel, predicts, se.fit=TRUE)[1]),
    cov2_se = unlist(predict(cov2_quadmodel, predicts, se.fit=TRUE)[2]),
    cov3_pred = unlist(predict(cov3_quadmodel, predicts, se.fit=TRUE)[1]),
    cov3_se = unlist(predict(cov3_quadmodel, predicts, se.fit=TRUE)[2]))

# 13C Convert SE to CIs
predicts <- predicts %>%                                              # Convert SE to CIs
  mutate(
    precov_lci = precov_pred-(1.96*precov_se),
    precov_uci = precov_pred+(1.96*precov_se),
    cov1_lci = cov1_pred-(1.96*cov1_se),
    cov1_uci = cov1_pred+(1.96*cov1_se),
    cov2_lci = cov2_pred-(1.96*cov2_se),
    cov2_uci = cov2_pred+(1.96*cov2_se),
    cov3_lci = cov3_pred-(1.96*cov3_se),
    cov3_uci = cov3_pred+(1.96*cov3_se))

# 13D Pivot data from wide to long for aov (analysis of variance) tests
aov_data <- predicts %>% dplyr::select(hour_of_day, weekend, weekday, precov_pred, cov1_pred,
                                       cov2_pred, cov3_pred) %>% 
  pivot_longer(precov_pred:cov3_pred, names_to = "covid_group", values_to = "red_predict") 

# 13E Comparison between covid groups for all hour_of_day and dow combined
aov_all <- aov(red_predict ~ covid_group, data = aov_data)
summary(aov_all)
TukeyHSD(aov_all)

# 13F Comparison between covid groups for weekend and weekday separately
aov_weekend <- aov_data %>% group_by(weekend) %>% 
  summarize(aov = list(aov(red_predict ~ covid_group)))
summary(aov_weekend$aov[[1]])
TukeyHSD(aov_weekend$aov[[1]])
summary(aov_weekend$aov[[2]])
TukeyHSD(aov_weekend$aov[[2]])

# 13G Comparison between covid groups for hour_of_day separately
aov_hourofday <- aov_data %>% group_by(hour_of_day) %>% 
  summarize(aov = list(aov(red_predict ~ covid_group)))
summary(aov_hourofday$aov[[1]])
summary(aov_hourofday$aov[[2]])
summary(aov_hourofday$aov[[3]])
summary(aov_hourofday$aov[[4]])
summary(aov_hourofday$aov[[5]])
summary(aov_hourofday$aov[[6]])
summary(aov_hourofday$aov[[7]])
summary(aov_hourofday$aov[[8]])

# Notes: aov analysis (13E-G) was not included in the manuscript


##########################################################################
# 14: STL Analysis (Seasonal and Trend decomposition using Loess)
##########################################################################

# Notes: For STL analysis, data should be complete (no missing). As we do
# have missing dates/times in each of the four time periods, we use the GAM
# models to impute missing data prior to STL analysis. We will include data
# for covid periods 1 and 2 in the STL analysis, as the goal is to identify
# social-distancing fatigue, or an increasing trend in congestion after initial 
# declines but before release of stay-at-home orders. Stay-at-home orders were
# released during covid period 2, so we include both periods 1 and 2 in the analysis

# 14A Impute missing data using prediction models
# 14A.1 Impute for pre_quadmodel (will likely not need)
table(is.na(gcc2_pre$red)) # 135 missing values
Pgcc2_pre <- as.vector(predict.gam(pre_quadmodel, newdata = gcc2_pre))
Fgcc2_pre <- gcc2_pre %>% mutate(red_pred = Pgcc2_pre,
                                 red_full = case_when(
                                   !is.na(red) ~ red,
                                   is.na(red) ~ red_pred
                                 )) %>% 
  relocate(red, red_pred, red_full, everything())      # inspect first three columns to ensure replacement is correct
table(is.na(Fgcc2_pre$red_full))                       # should be 0 missing values now

# 14A.2 Impute for cov1_quadmodel
table(is.na(gcc2_cov1$red)) # 52 missing values
Pgcc2_cov1 <- as.vector(predict.gam(cov1_quadmodel, newdata = gcc2_cov1))
Fgcc2_cov1 <- gcc2_cov1 %>% mutate(red_pred = Pgcc2_cov1,
                                   red_full = case_when(
                                     !is.na(red) ~ red,
                                     is.na(red) ~ red_pred
                                     )) %>% 
  relocate(red, red_pred, red_full, everything())
table(is.na(Fgcc2_cov1$red_full))

# 14A.3 Impute for cov2_quadmodel
table(is.na(gcc2_cov2$red)) # 2 missing values
Pgcc2_cov2 <- as.vector(predict.gam(cov2_quadmodel, newdata = gcc2_cov2))
Fgcc2_cov2 <- gcc2_cov2 %>% mutate(red_pred = Pgcc2_cov2,
                                   red_full = case_when(
                                     !is.na(red) ~ red,
                                     is.na(red) ~ red_pred
                                   )) %>% 
  relocate(red, red_pred, red_full, everything())
table(is.na(Fgcc2_cov2$red_full))

# 14A.4 Imptue for cov3_quadmodel (will likely not need)
table(is.na(gcc2_cov3$red)) # 63 missing values
Pgcc2_cov3 <- as.vector(predict.gam(cov3_quadmodel, newdata = gcc2_cov3))
Fgcc2_cov3 <- gcc2_cov3 %>% mutate(red_pred = Pgcc2_cov3,
                                   red_full = case_when(
                                     !is.na(red) ~ red,
                                     is.na(red) ~ red_pred
                                   )) %>% 
  relocate(red, red_pred, red_full, everything())
table(is.na(Fgcc2_cov3$red_full))

# 14B Combine cov1 and cov2 into one dataframe
Fgcc2_cov1and2 <- bind_rows(Fgcc2_cov1, Fgcc2_cov2)
table(is.na(Fgcc2_cov1and2$red_full))

# 14C Turn the time series into a ts object
# Notes: The frequency is the number of observations per unit of time; there are 
# 8 obs per 24 hour period, so the freq is 8
traf.ts <- ts(Fgcc2_cov1and2$red_full, frequency = 8)
plot(traf.ts)

# 14D Run STL
traf.stl <- stl(traf.ts, s.window = "periodic", robust = TRUE, 
                inner = 1, outer = 15, na.action = na.fail, l.window = 9)

# 14E Convert STL output to df
# Note that "seasonal" column is daily in our case
stl.results <- as.data.frame(traf.stl$time.series)

# 14F Add variables to STL df for use in plotting
stl.results <- stl.results %>% 
  mutate(
    data = traf.ts,                                  # add raw data
    sum = seasonal + trend + remainder,              # add sum to make sure decomp worked properly
    scaled_daily_traffic = mean(trend) + seasonal,   # scale the daily data (may not use)
    date_time = Fgcc2_cov1and2$date_time             # add timestamp
   )

# 14G View results graphically
plot(stl.results$date_time,stl.results$trend,type='l')
plot(stl.results$date_time,stl.results$seasonal,type='l')
plot(stl.results$date_time,stl.results$remainder,type='l')


##########################################################################
# 15: Figures for Manuscript
##########################################################################

# Notes: In this chunk we create the figures for the manuscript. We will create:
# Figure 1: (from Markus - created in MatLab) Two-panel map of Manhattan comparing
#           congestion March 9-13 and March 23-27
# Figure 2: Combined time series and timeline plot, with green, orange, red, maroon
#           panels, as well as the four model date periods and key intervention dates
# Figure 3: STL plot, four panels, for data, trend, seasonal, and remainder components
# Figure 4: Prediction plot showing "typical" traffic by hour of day for each time
#           period, faceted by weekday and weekend
# Supplemental Figure 1: Plot of splines (top row) and percent deviation of spline from the 
#                        intercept (bottom row), faceted by weekday and weekend

# 15A Figure 2: Combined timeseries and timeline plot
theme_set(theme_bw(base_size = 12))
time_plot <- gcc2 %>% 
  pivot_longer(maroon:gray, names_to = "gcc_color", values_to = "percent_roads") %>% 
  filter(gcc_color != "gray") %>%  
  mutate(gcc_color = factor(gcc_color, levels = c("maroon", "red", "orange", "green"))) %>% 
  ggplot(aes(x = date_time, y = percent_roads, color = gcc_color)) +
  geom_line() +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-03-14"))), col = "black") +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-05-20"))), col = "black") +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-06-17"))), col = "black") +
  facet_wrap(~ gcc_color, ncol = 1, scales = "free") +
  scale_color_manual(values=c("maroon", "red", "orange", "green")) +
  scale_x_datetime(date_breaks = "1 month") +
  xlab("Date") + ylab("% Area with Color Coverage") +
  guides(color = FALSE) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
time_plot

tiff("./figures/time_plot.tiff", 
     units = "in", width = 12, height = 8, res = 300)
time_plot
dev.off()
# Notes: Some components of this plot (timeline elements) were added in PowerPoint

# 15B Figure 3 STL plot
theme_set(theme_bw(base_size = 12))

# 15B.1 Top panel - data plot
data_plot <- stl.results %>% 
  ggplot(aes(x = date_time, y = data)) +
  geom_line(color = "red") +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-06-08"))), 
             col = "black", linetype = "dashed") +
  annotate(geom = "text", x = as_datetime("2020-06-09"),
           y = 0.17, label = "Phase 1 Begins", hjust = 0) +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-03-22"))), 
             col = "black", linetype = "dashed") +
  annotate(geom = "text", x = as_datetime("2020-03-23"),
           y = 1.20, label = "NY on PAUSE Begins", hjust = 0) +
  scale_x_datetime(date_breaks = "1 week") +
  xlab("") + ylab("Data") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
data_plot

# 15B.2 2nd panel - trend plot
trend_plot <- stl.results %>% 
  ggplot(aes(x = date_time, y = trend)) +
  geom_line(color = "#56B4E9") +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-06-08"))), 
             col = "black", linetype = "dashed") +
  annotate(geom = "text", x = as_datetime("2020-06-09"),
           y = 0.32, label = "Phase 1 Begins", hjust = 0) +
  geom_vline(aes(xintercept = as.integer(as.POSIXct("2020-03-22"))), 
             col = "black", linetype = "dashed") +
  annotate(geom = "text", x = as_datetime("2020-03-23"),
           y = 0.58, label = "NY on PAUSE Begins", hjust = 0) +
  scale_x_datetime(date_breaks = "1 week") +
  xlab("") + ylab("Trend") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
trend_plot

# 15B.3 3rd panel - seasonal plot
seasonal_plot <- stl.results %>% 
  ggplot(aes(x = date_time, y = seasonal)) +
  geom_line(color = "#56B4E9") +
  scale_x_datetime(date_breaks = "1 week") +
  xlab("") + ylab("Periodic") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
seasonal_plot

# 15B.4 Bottom panel - remainder plot
remainder_plot <- stl.results %>% 
  ggplot(aes(x = date_time, y = remainder)) +
  geom_line(color = "#56B4E9") +
  scale_x_datetime(date_breaks = "1 week") +
  xlab("Date") + ylab("Remainder") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
remainder_plot

# 15B.5 Arrange into one plot and save
tiff("./figures/stl_plot.tiff", 
     units = "in", width = 12, height = 8, res = 300)
stl_plot <- 
  grid.arrange(arrangeGrob(data_plot, trend_plot,
                           seasonal_plot, remainder_plot, ncol = 1,
                           left = textGrob("Red traffic congestion color coverage (% area)", rot = 90, vjust = 1)))
dev.off()

# 15C Figure 4: Prediction plot for "typical" hourly traffic
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73") # color blind friendly palette (four colors)
pred_plot1 <- predicts %>% 
  mutate(weekend = as_factor(as.character(weekend))) %>% 
  dplyr::select(-weekday, -weekend_hourofday, -precov_se, -cov1_se,
                -cov2_se, -cov3_se, -weekday_hourofday) %>% 
  pivot_longer(!hour_of_day & !weekend,
               names_to = c("timeline",".value"),
               names_sep = "_") %>% 
  mutate(weekend = ifelse(weekend == 1, "Weekend", "Weekday"),
         timeline = case_when(timeline == "precov" ~ "Pre-COVID",
                              timeline == "cov1" ~ "COVID Period 1",
                              timeline == "cov2" ~ "COVID Period 2",
                              timeline == "cov3" ~ "COVID Period 3"),
         timeline = factor(timeline, levels = c("Pre-COVID", "COVID Period 1", "COVID Period 2", "COVID Period 3"))) %>% 
  ggplot(aes(x = hour_of_day, y = pred, color = timeline)) + 
  geom_point() + 
  geom_line() +
  geom_ribbon(aes(ymin = lci, ymax = uci), linetype = 2, alpha = 0.1) +
  facet_wrap(~weekend) +
  xlab("Hour of day") + ylab("GAM predictions of % red traffic congestion color coverage") +
  labs(color = "COVID time period") + 
  scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"), strip.background = element_rect(fill="white"))
pred_plot1

tiff("./figures/pred_plot.tiff", 
     units = "in", width = 8, height = 6, res = 300)
pred_plot1
dev.off()

# 15D Supplemental Figure 1: Spline plot
# 15D.1 Pull data to plot splines from models
p1_data <- visreg(pre_quadmodel, "weekday_hourofday", plot = FALSE)
p3_data <- visreg(cov1_quadmodel, "weekday_hourofday", plot = FALSE)
p5_data <- visreg(cov2_quadmodel, "weekday_hourofday", plot = FALSE)
p7_data <- visreg(cov3_quadmodel, "weekday_hourofday", plot = FALSE)
p2_data <- visreg(pre_quadmodel, "weekend_hourofday", plot = FALSE)
p4_data <- visreg(cov1_quadmodel, "weekend_hourofday", plot = FALSE)
p6_data <- visreg(cov2_quadmodel, "weekend_hourofday", plot = FALSE)
p8_data <- visreg(cov3_quadmodel, "weekend_hourofday", plot = FALSE)

# 15D.2 Bind data from models together: one df for weekday, one for weekend
# and add variables to adjust by intercept
dplyr::bind_rows(
  dplyr::mutate(p1_data$fit, model = "Pre-COVID",
                visregFit2 = visregFit-0.986174,
                visregLwr2 = visregLwr-0.986174,
                visregUpr2 = visregUpr-0.986174,
                visregFit3 = (visregFit-0.986174)/0.986174,
                visregLwr3 = (visregLwr-0.986174)/0.986174,
                visregUpr3 = (visregUpr-0.986174)/0.986174),
  dplyr::mutate(p3_data$fit, model = "COVID Period 1",
                visregFit2 = visregFit-0.414462,
                visregLwr2 = visregLwr-0.414462,
                visregUpr2 = visregUpr-0.414462,
                visregFit3 = (visregFit-0.414462)/0.414462,
                visregLwr3 = (visregLwr-0.414462)/0.414462,
                visregUpr3 = (visregUpr-0.414462)/0.414462),
  dplyr::mutate(p5_data$fit, model = "COVID Period 2",
                visregFit2 = visregFit-0.557912,
                visregLwr2 = visregLwr-0.557912,
                visregUpr2 = visregUpr-0.557912,
                visregFit3 = (visregFit-0.557912)/0.557912,
                visregLwr3 = (visregLwr-0.557912)/0.557912,
                visregUpr3 = (visregUpr-0.557912)/0.557912),
  dplyr::mutate(p7_data$fit, model = "COVID Period 3",
                visregFit2 = visregFit-0.739747,
                visregLwr2 = visregLwr-0.739747,
                visregUpr2 = visregUpr-0.739747,
                visregFit3 = (visregFit-0.739747)/0.739747,
                visregLwr3 = (visregLwr-0.739747)/0.739747,
                visregUpr3 = (visregUpr-0.739747)/0.739747)
) -> spline_plot_data_weekday 
dplyr::bind_rows(
  dplyr::mutate(p2_data$fit, model = "Pre-COVID",
                visregFit2 = visregFit-0.986174,
                visregLwr2 = visregLwr-0.986174,
                visregUpr2 = visregUpr-0.986174,
                visregFit3 = (visregFit-0.986174)/0.986174,
                visregLwr3 = (visregLwr-0.986174)/0.986174,
                visregUpr3 = (visregUpr-0.986174)/0.986174),
  dplyr::mutate(p4_data$fit, model = "COVID Period 1",
                visregFit2 = visregFit-0.414462,
                visregLwr2 = visregLwr-0.414462,
                visregUpr2 = visregUpr-0.414462,
                visregFit3 = (visregFit-0.414462)/0.414462,
                visregLwr3 = (visregLwr-0.414462)/0.414462,
                visregUpr3 = (visregUpr-0.414462)/0.414462),
  dplyr::mutate(p6_data$fit, model = "COVID Period 2",
                visregFit2 = visregFit-0.557912,
                visregLwr2 = visregLwr-0.557912,
                visregUpr2 = visregUpr-0.557912,
                visregFit3 = (visregFit-0.557912)/0.557912,
                visregLwr3 = (visregLwr-0.557912)/0.557912,
                visregUpr3 = (visregUpr-0.557912)/0.557912),
  dplyr::mutate(p8_data$fit, model = "COVID Period 3",
                visregFit2 = visregFit-0.739747,
                visregLwr2 = visregLwr-0.739747,
                visregUpr2 = visregUpr-0.739747,
                visregFit3 = (visregFit-0.739747)/0.739747,
                visregLwr3 = (visregLwr-0.739747)/0.739747,
                visregUpr3 = (visregUpr-0.739747)/0.739747)
) -> spline_plot_data_weekend

# 15D.3 Create weekday plot - Panel A (no adjustment for intercept)
theme_set(theme_bw(base_size = 12))
labela <- "A: Weekday"
groba <- grid.text(labela, x=0.1,  y=0.95, gp=gpar(col="black", fontsize=12))
spline_plot_weekday <- spline_plot_data_weekday %>% 
  mutate(model = factor(model, levels = c("Pre-COVID", "COVID Period 1", "COVID Period 2", "COVID Period 3"))) %>%
  ggplot() +
  geom_ribbon(aes(weekday_hourofday, ymin=visregLwr, ymax=visregUpr, group=model), 
              fill="gray90") +
  geom_line(aes(weekday_hourofday, visregFit, group=model, color=model)) +
  labs(x = "Hour of day", y = "GAM fitted spline") +
  theme(legend.position="none") + 
  annotation_custom(groba) + 
  scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
spline_plot_weekday

# 15D.4 Create weekend plot - Panel B (no adjustment for intercept)
theme_set(theme_bw(base_size = 12))
labelb <- "B: Weekend"
grobb <- grid.text(labelb, x=0.1,  y=0.95, gp=gpar(col="black", fontsize=12))
spline_plot_weekend <- spline_plot_data_weekend %>% 
  mutate(model = factor(model, levels = c("Pre-COVID", "COVID Period 1", "COVID Period 2", "COVID Period 3"))) %>%
  ggplot() +
  geom_ribbon(aes(weekend_hourofday, ymin=visregLwr, ymax=visregUpr, group=model), 
              fill="gray90") +
  geom_line(aes(weekend_hourofday, visregFit, group=model, color=model)) +
  labs(x = "Hour of day", y = "GAM fitted spline") +
  theme(legend.position="none") + 
  annotation_custom(grobb) +
  scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
spline_plot_weekend

# 15D.5 Create weekday plot - Panel C (Deviation from intercept)
theme_set(theme_bw(base_size = 12))
labelc <- "C: Weekday"
grobc <- grid.text(labelc, x=0.1,  y=0.95, gp=gpar(col="black", fontsize=12))
spline_plot_weekday.3 <- spline_plot_data_weekday %>% 
  mutate(model = factor(model, levels = c("Pre-COVID", "COVID Period 1", "COVID Period 2", "COVID Period 3"))) %>%
  ggplot() +
  geom_ribbon(aes(weekday_hourofday, ymin=visregLwr3, ymax=visregUpr3, group=model), 
              fill="gray90") +
  geom_line(aes(weekday_hourofday, visregFit3, group=model, color=model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  labs(x = "Hour of day", y = "Deviation of fitted spline from intercept") +
  theme(legend.position="none") + 
  annotation_custom(grobc) +
  scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
spline_plot_weekday.3

# 15D.6 Create weekend plot - Panel D (Deviation from intercept)
theme_set(theme_bw(base_size = 12))
labeld <- "D: Weekend"
grobd <- grid.text(labeld, x=0.1,  y=0.95, gp=gpar(col="black", fontsize=12))
spline_plot_weekend.3 <- spline_plot_data_weekend %>% 
  mutate(model = factor(model, levels = c("Pre-COVID", "COVID Period 1", "COVID Period 2", "COVID Period 3"))) %>%
  ggplot() +
  geom_ribbon(aes(weekend_hourofday, ymin=visregLwr3, ymax=visregUpr3, group=model), 
              fill="gray90") +
  geom_line(aes(weekend_hourofday, visregFit3, group=model, color=model)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  labs(x = "Hour of day", y = "Deviation of fitted spline from intercept", 
       color = "COVID time period") +
  annotation_custom(grobd) +
  scale_colour_manual(values=cbPalette) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.justification=c(1,0), 
        legend.position=c(0.9, 0.05),  
        legend.background = element_blank(),
        legend.key = element_blank())
spline_plot_weekend.3

# 15D.7 Arrange into one plot and save
tiff("./figures/spline_plot.tiff", 
     units = "in", width = 12, height = 8, res = 300)
spline_plot <- 
  grid.arrange(arrangeGrob(spline_plot_weekday, spline_plot_weekend,
                           spline_plot_weekday.3, spline_plot_weekend.3, ncol = 2))
dev.off()


##########################################################################
# 16: Table 1 for Manuscript
##########################################################################

# Notes: As table 1 is very small, it is formatted in MS word, and not in R.
# Data needed for the table is below

# 16A Model 1
summary(pre_quadmodel)
min(gcc2_pre$date_time)
max(gcc2_pre$date_time)

# 16B Model 2
summary(cov1_quadmodel)
min(gcc2_cov1$date_time)
max(gcc2_cov1$date_time)

# 16C Model 3
summary(cov2_quadmodel)
min(gcc2_cov2$date_time)
max(gcc2_cov2$date_time)

# 16D Model 4
summary(cov3_quadmodel)
min(gcc2_cov3$date_time)
max(gcc2_cov3$date_time)



