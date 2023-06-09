## Modeling Language Vitality using Social and Environmental Predictors
## Kirsten Helgeson with Dr. Bradley Rentz and Dr. Gary Holton
## Created using R 3.5.0, R 3.6.1, and R 4.1.3 starting in November 2018, adapted through April 2022

# Load packages
library(dplyr)
library(ggplot2)

# Set working directory
setwd("~/Desktop/Supplementary Materials")

data.clean <- read.csv ("model_data.csv")

options(scipen = 999)  # This turns off scientific notation.

#### Part 1: Setting up data ####

## 1a. Language vitality ##

# Make language status levels into ordinal data
data.clean$Vitality_rank <- factor(data.clean$Vitality_rank, levels=c("0","1","2","3","4","5","6","7"), ordered=T)
plot(data.clean$Vitality_rank)
summary(data.clean$Vitality_rank)

## 1b. Subsistence strategy ##

# Set subsistence strategy to nominal/categorical (rather than interval) data
data.clean$Subsistence_cat <- factor(data.clean$Subsistence_cat, levels=c("Agriculture","Hunting","Gathering","Fishing","Mixed"), ordered=F) 
plot(data.clean$Subsistence_cat)
summary(data.clean$Subsistence_cat)

## 1c. Land area ##

hist(data.clean$Area_sum, breaks = 50, main = "Land Area at Time of European Contact", xlab = "Area", ylab = "Frequency")
summary(data.clean$Area_sum)

## Square root of land area
data.clean$Area_sum_sqrt <- sqrt(data.clean$Area_sum)
hist(data.clean$Area_sum_sqrt, breaks = 50, main = "Land Area at Time of European Contact", xlab = "Area (sqrt)", ylab = "Frequency")

## Log of land area
data.clean$Area_sum_log <- log(data.clean$Area_sum)
hist(data.clean$Area_sum_log, breaks = 50, main = "Land Area at Time of European Contact", xlab = "Log area", ylab = "Frequency")

## 1d. Elevation range ##

hist(data.clean$Elevation_range, breaks = 50, main = "Elevation range", xlab = "Elevation (meters)", ylab = "Frequency")
summary(data.clean$Elevation_range)

## Square root of elevation range
data.clean$Elevation_range_sqrt <- sqrt(data.clean$Elevation_range)
hist(data.clean$Elevation_range_sqrt, breaks = 50, main = "Elevation range", xlab = "Square root of elevation (meters)", ylab = "Frequency")

## Log of elevation range
data.clean$Elevation_range_log <- log(data.clean$Elevation_range)
hist(data.clean$Elevation_range_log, breaks = 50, main = "Elevation range", xlab = "Log Elevation", ylab = "Frequency")

## 1e. Displacement ##

# Set historical displacement to categorical and reference level to "none."
data.clean$Displacement_cat <- factor(data.clean$Displacement_cat, levels=c("none","outside","inside","both"), ordered=F) # Change to "ordered=T" for ordinal data.
summary(data.clean$Displacement_cat)
plot(data.clean$Displacement_cat)

## 1f. Highways ##

## Highway density

# Explore data
summary(data.clean$Hwy_density)
hist(data.clean$Hwy_density, breaks = 50, main = "Highway density", xlab = "Density", ylab = "Frequency")

# Scaled highway density
data.clean$Hwy_density_scaled <- scale(data.clean$Hwy_density, center=F, scale=T)
hist(data.clean$Hwy_density_scaled, breaks = 30)

# Square root of highway density
data.clean$Hwy_density_sqrt <- sqrt(data.clean$Hwy_density)
hist(data.clean$Hwy_density_sqrt, breaks = 50, main = "Highway density", xlab = "Square root of density", ylab = "Frequency")
summary(data.clean$Hwy_density_sqrt)

# Log of highway density
data.clean$Hwy_density_log <- log(data.clean$Hwy_density)
hist(data.clean$Hwy_density_log, breaks = 50, main = "Highway density", xlab = "Log of density", ylab = "Frequency")
## The problem with log transforming is that most highway density values are zero, so there are many values with infinitely negative values on the log transformation scale.

# Adding a small number
data.clean$Hwy_density_new <- (data.clean$Hwy_density+1)
data.clean$Hwy_density_new_log <- log(data.clean$Hwy_density_new)
hist(data.clean$Hwy_density_new_log, breaks = 100)

data.clean$Hwy_density_newd <- (data.clean$Hwy_density+0.000001)
data.clean$Hwy_density_newd_log <- log(data.clean$Hwy_density_newd)
hist(data.clean$Hwy_density_newd_log, breaks = 10)

## Highway presence
summary(data.clean$Highway)
plot(data.clean$Highway)

## 1g. Revitalization ##

# Set revitalization to ordinal, "none" as reference level.
data.clean$Revit_cat <- factor(data.clean$Revit_cat, levels=c("none", "revitalization", "education", "immersion"), ordered=F) # Change to "ordered=T" for ordinal data.
summary(data.clean$Revit_cat)
plot(data.clean$Revit_cat)

# Define variables in R.
Vitality <- data.clean$Vitality_rank
Area <- data.clean$Area_sum_log
Elevation <- data.clean$Elevation_range_sqrt
Subsistence <- data.clean$Subsistence_cat
Displacement <- data.clean$Displacement_cat
Highway <- data.clean$Highway
Revitalization <- data.clean$Revit_cat

#### Part 2: Running models ####

library(ordinal)
library(sure) # For model diagnostics and checking model fit
library(sjPlot) # For visualizing results

#### Model 5AESDR ####

Model5AESDR <- clm(Vitality~Area+Elevation+Subsistence+Displacement+Revitalization, data=data.clean)
summary(Model5AESDR)

plot_model(Model5AESDR)
tab_model(Model5AESDR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model5AESDR)/logLik(nullmod)
# 'log Lik.' 0.5688264 (df=19)

# Fit diagnostics
autoplot.clm(Model5AESDR, what = "qq")
autoplot.clm(Model5AESDR, what = "fitted")

#### Model 6AESHR ####

Model6AESHR <- clm(Vitality~Area+Elevation+Subsistence+Highway+Revitalization, data=data.clean)
summary(Model6AESHR)

plot_model(Model6AESHR)
tab_model(Model6AESHR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model6AESHR)/logLik(nullmod)
# 'log Lik.' 0.6118705 (df=17)

# Fit diagnostics
autoplot.clm(Model6AESHR, what = "qq")
autoplot.clm(Model6AESHR, what = "fitted")

#### Model 6ASHR ####

# Same as 6AESHR but without elevation range.
Model6ASHR <- clm(Vitality~Area+Subsistence+Highway+Revitalization, data=data.clean)
summary(Model6ASHR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model6ASHR)/logLik(nullmod)
# 'log Lik.' 0.6097072 (df=16)

# Fit diagnostics
autoplot.clm(Model6ASHR, what = "qq")
autoplot.clm(Model6ASHR, what = "fitted")

#### Model 6AE*SHR ####

## Model 6AE*SHR: Testing interaction effects

Model6AExSHR <- clm(Vitality~Area+Elevation*Subsistence+Highway+Revitalization, data=data.clean)
summary(Model6AExSHR)

plot_model(Model6AExSHR)
tab_model(Model6AExSHR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model6AExSHR)/logLik(nullmod)
# 'log Lik.' 0.6324988 (df=21)

# Fit diagnostics
autoplot.clm(Model6AExSHR, what = "qq")
autoplot.clm(Model6AExSHR, what = "fitted")


#### Graphing interaction effects for Model 6AE*SHR ####
# Note: This section and the following one use script adapted from the tutorial “Decomposing, Probing, and Plotting Interactions in R” acquired from UCLA Advanced Research Computer Statistical Methods and Data Analytics on January 17, 2022 at the following URL: https://stats.oarc.ucla.edu/r/seminars/interactions-r/.

library(emmeans)

# Pairwise comparisons for elevation range based on subsistence
emtrends(Model6AExSHR, pairwise ~ Subsistence, var="Elevation")

# Define variables
(mylist <- list(Elevation=c(0,20,40,60,80),Subsistence=c("Hunting","Gathering","Fishing","Mixed","Agriculture")))
emModel6AExSHR <- emmeans(Model6AExSHR, ~ Elevation*Subsistence, at=mylist)
contrast(emModel6AExSHR, "revpairwise", by="Elevation")

# Plot using emmip
(mylist <- list(Elevation=seq(0,80,by=10),Subsistence=c("Hunting","Gathering","Fishing","Mixed","Agriculture")))
emmip(Model6AExSHR, Subsistence ~Elevation, at=mylist, CIs=TRUE, xlab = "Elevation range (square root)")

#### Model fit Diagnostics ####

# Fit diagnostics
autoplot.clm(Model, what = "qq")
autoplot.clm(Model, what = "fitted")


#### Model 4AESR ####

Model4AESR <- clm(Vitality~Area+Elevation+Subsistence+Revitalization, data=data.clean)
summary(Model4AESR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model4AESR)/logLik(nullmod)
# 'log Lik.' 0.4605134 (df=16)

plot_model(Model4AESR)
tab_model(Model4AESR)

# Fit diagnostics
autoplot.clm(Model4AESR, what = "qq")
autoplot.clm(Model4AESR, what = "fitted")

#### Model 4AE*SR ####

Model4AExSR <- clm(Vitality~Area+Elevation*Subsistence+Revitalization, data=data.clean)
summary(Model4AExSR)

plot_model(Model4AExSR)
tab_model(Model4AExSR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model4AExSR)/logLik(nullmod)
# 'log Lik.' 0.4960515 (df=20)

# Model 1e fit diagnostics
autoplot.clm(Model4AExSR, what = "qq")
autoplot.clm(Model4AExSR, what = "fitted")

## Can try this again with Gathering set to reference level (since this level shows the strongest effect in 4AE*SR).

data.clean$Subsistence_cat <- factor(data.clean$Subsistence_cat, levels=c("Gathering","Hunting","Mixed","Fishing","Agriculture"), ordered=F)
summary(data.clean$Subsistence_cat)

#### Graphing interaction effects for Model4AE*SR ####

library(emmeans)

# Pairwise comparisons for elevation range based on subsistence.
emtrends(Model4AExSR, pairwise ~ Subsistence, var="Elevation")

# Define variables 
(mylist <- list(Elevation_range_sqrt=c(0,20,40,60,80),Subsistence=c("Hunting","Gathering","Fishing","Mixed","Agriculture")))
emModel4AExSR <- emmeans(Model4AExSR, ~ Elevation*Subsistence, at=mylist)
contrast(emModel4AExSR, "revpairwise", by="Elevation")

# Plot using emmip
(mylist <- list(Elevation=seq(0,80,by=10),Subsistence=c("Hunting","Gathering","Fishing","Mixed","Agriculture")))
emmip(Model4AExSR, Subsistence~Elevation, at=mylist, CIs=TRUE, xlab = "Elevation range (square root)")

#### Model 2AEDR ####

Model2AEDR <- clm(Vitality~Area+Elevation+Displacement+Revitalization, data=data.clean)
summary(Model2AEDR)

# Calculate McFadden's R2
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model2AEDR)/logLik(nullmod)
# 'log Lik.' 0.3613663 (df=15)

plot_model(Model2AEDR)
tab_model(Model2AEDR)

## Fit diagnostics ##
autoplot.clm(Model2AEDR, what = "qq")
autoplot.clm(Model2AEDR, what = "fitted")

#### Model 2ADR ####

Model2ADR <- clm(Vitality~Area+Displacement+Revitalization, data=data.clean)
summary(Model2ADR)

# Calculate McFadden's R2
Model2ADR <- clm(Vitality~Area+Displacement+Revitalization, data=data.clean)
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model2ADR)/logLik(nullmod)
# 'log Lik.' 0.3609396 (df=14)

## Fit diagnostics ##
autoplot.clm(Model2ADR, what = "qq")
autoplot.clm(Model2ADR, what = "fitted")

#### Model 3AEHR ####

Model3AEHR <- clm(Vitality~Area+Elevation+Highway+Revitalization, data=data.clean)
summary(Model3AEHR)

# Calculate McFadden's R2
Model3AEHR <- clm(Vitality~Area+Elevation+Highway+Revitalization, data=data.clean)
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model3AEHR)/logLik(nullmod)
# 'log Lik.' 0.4352562 (df=13)

plot_model(Model3AEHR)
tab_model(Model3AEHR)

## Fit diagnostics ##
autoplot.clm(Model3AEHR, what = "qq")
autoplot.clm(Model3AEHR, what = "fitted")

#### Model 3AHR ####

Model3AHR <- clm(Vitality~Area+Highway+Revitalization, data=data.clean)
summary(Model3AHR)

# Calculate McFadden's R2
Model3AHR <- clm(Vitality~Area+Highway+Revitalization, data=data.clean)
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model3AHR)/logLik(nullmod)
# 'log Lik.' 0.4339411 (df=12)

plot_model(Model3AHR)
tab_model(Model3AHR)

## Fit diagnostics 
autoplot.clm(Model3AHR, what = "qq")
autoplot.clm(Model3AHR, what = "fitted")

#### Model 1AER ####

Model1AER <- clm(Vitality~Area+Elevation+Revitalization, data=data.clean)
summary(Model1AER)

# Calculate McFadden's R2.
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model1AER)/logLik(nullmod)
# 'log Lik.' 0.2145197 (df=12)

## Fit diagnostics ##
autoplot.clm(Model1AER, what = "qq")
autoplot.clm(Model1AER, what = "fitted")

#### Model 1AR ####

Model1AR <- clm(Vitality~Area+Revitalization, data=data.clean)
summary(Model1AR)

# Calculate McFadden's R2.
nullmod <- clm(Vitality~1, data=data.clean)
1-logLik(Model1AR)/logLik(nullmod)
# 'log Lik.' 0.2142884 (df=11)

## Fit diagnostics ##
autoplot.clm(Model1AR, what = "qq")
autoplot.clm(Model1AR, what = "fitted")

#### Plotting marginal means ####

# Load packages
library(emmeans)
library(broom)
library(tidyverse)


## Model 4AE*SR: This is the model to focus on, since the presence or absence of highways was not a statistically significant predictor in 6AE*SHR.
emmeans4AExSR <- tidy(emmeans(Model4AExSR, ~Vitality | Elevation*Subsistence,
                              type="response",mode="prob",
                              at = list(Elevation=c(10,20,30,40,50,60,70)))) # can pick different levels of continuous variable to calculate the marginal means for

# Pairwise comparisons of marginal means
pairs(emmeans(Model4AExSR, ~Elevation | Vitality, type="response",mode="prob",at = list(Elevation=c(10,20,30,40,50,60,70)))) 
pairs(emmeans(Model4AExSR, ~Elevation*Subsistence | Vitality, type="response",mode="prob",at = list(Elevation=c(10,20,30,40,50,60,70))))

# Plot
ggplot(emmeans4AExSR, aes(x=Elevation,y=prob*100,color=Subsistence))+ 
  geom_ribbon(aes(ymin=(prob-std.error)*100, # the geom_ribbon add the error shading
                  ymax=(prob+std.error)*100,
                  fill=Subsistence),
              alpha=0.2)+ 
  geom_line() + # this adds the main plotted marginal mean lines
  facet_wrap(~Vitality) + 
  theme_minimal() + 
  ylab("Probability of Vitality Level") +
  xlab("Elevation range (square root)")