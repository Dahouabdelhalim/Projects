#Mitchem et al. female preference data
#Author: Lisa D Mitchem
#email: lm7en@virginia.edu


### The purpose of this script is to process female chemical cue preference data from the summer of 2019. Data collected by Liza Mitchem (lm7en@virginia.edu) and Zorimar Vilella-Pacheco (zorimar.vilella@upr.edu)  

#**With this script, we answer the following questions:**  
  
#  1. Do females detect and associate with male chemical cues? 
#  2. Do females prefer to associate with winning or losing males before combat?  
#  3. Do females prefer to associate with winning or losing males after combat?  
#  4. Do females prefer to associate with more interactive males?  

rm(list=ls())

library(emmeans)
library(ggplot2)
library(reshape2)
library(plyr)
library(lme4)
library(car)
library(effects)


####Do females approach male scents more than control scents?####

counts <- read.csv("file:///C:/Users/mitch/Documents/PhD Stuff/Brodie Lab/Summer 2020/female preference final MS stuff/R data/female_preference_counts_data.csv")

#trial round is structured as a number, so we need to change it to a factor so the model runs correctly
counts$trial.round <- as.factor(counts$trial.round)

#Independent variable = filter paper (control or male) and trial round (each female was tested twice, so there are 2 trial rounds)  
#Dependent variable = number of approaches  
#Random effects = trial (each trial is represented twice in the data - once for the female approaches to control filter paper, and once for the female approaches to male filter paper)  

counts.model <- glmer(value ~ treatment + trial.round + treatment*trial.round + (1|trial), control = glmerControl(optimizer = "bobyqa"), family = "poisson", data = counts)
summary(counts.model)
Anova(counts.model, type = "III")

####Do females spend more TIME on male versus control scents?####

duration <- read.csv("file:///C:/Users/mitch/Documents/PhD Stuff/Brodie Lab/Summer 2020/female preference final MS stuff/R data/female_preference_duration_data.csv")


#trial round is structured as a number, so we need to change it to a factor so the model runs correctly
duration$trial.round <- as.factor(duration$trial.round)
#Independent variable = filter paper (control or male), and trial round (each female was tested twice, so there are 2 trial rounds)   
#Dependent variable = time spent on filter paper  
#Random effects = trial (each trial is represented twice in the data - once for the female approaches to control filter paper, and once for the female approaches to male filter paper)  

model.d <- glmer(value ~ treatment + trial.round + treatment*trial.round + (1|trial), control = glmerControl(optimizer = "bobyqa"), family = "poisson", data = duration)
summary(model.d)
Anova(model.d, type = "III")

plot(resid(model.d))
qqnorm(residuals(model.d))
plot(fitted(model.d),residuals(model.d))
hist(residuals(model.d))

####Do females prefer winning males####

competition.data <- read.csv("file:///C:/Users/mitch/Documents/PhD Stuff/Brodie Lab/Summer 2020/female preference final MS stuff/R data/male_competition_data.csv")


#trial round is structured as a number, so we need to change it to a factor so the model runs correctly
competition.data$trial.round <- as.factor(competition.data$trial.round)

#Independent variable = variable - filter paper (winner male, loser male, or control)  
#Dependent variable = value - time spent on filter paper  
#Random effects = female (each female is represented twice in the data - once for the female approaches to filter paper before, and once for the female approaches to filter paper after competition trials)  

status.model2 <- glmer(value ~ variable + trial.round + variable*trial.round + (1|female), control = glmerControl(optimizer = "bobyqa"), family = "poisson", data = competition.data)
summary(status.model2)
Anova(status.model2, type = "III")


#Testing model assumptions 
plot(resid(status.model2))
qqnorm(residuals(status.model2))
plot(fitted(status.model2),residuals(status.model2))

#pairwise comparisons of categorical independent variable
status.emm <- emmeans(status.model2, c("variable", "trial.round"))

pairs(status.emm)


####Do females prefer more interactive males?####

interact <- read.csv("file:///C:/Users/mitch/Documents/PhD Stuff/Brodie Lab/Summer 2020/female preference final MS stuff/R data/interaction_data.csv")


#trial round is structured as a number, so we need to change it to a factor so the model runs correctly
interact$trial.round <- as.factor(interact$trial.round)

#same variables as winning/losing male model, but changing the independent variable to be the beetle who initiated more, initiated less, and the control filter paper

int.model <- glmer(value ~ variable + trial.round + variable*trial.round + (1|female), control = glmerControl(optimizer = "bobyqa"), family = "poisson", data = interact)
summary(int.model)
Anova(int.model, type = "III")

#testing model assumptions
plot(resid(int.model))
qqnorm(residuals(int.model))
plot(fitted(int.model),residuals(int.model))

hist(residuals(int.model))

#pairwise comparisions of categorical independent variable

status.emm.int <- emmeans(int.model, c("variable", "trial.round"))

pairs(status.emm.int)

####Do females spend a non-random amount of time on filter paper?####  

#-Arenas are 17 cm X 15 cm, so 255 cm^2
#-Filter papers are 5 cm diameter triangles, so 20 cm^2
#-% filter paper per container = ~8%
#  -8% of 2 hours = 9:40 minutes

#Here I did a one-sample t-test where I compared the time on filter paper to a 'random' mean of 9.8 minutes

filter.paper <- read.csv("file:///C:/Users/mitch/Documents/PhD Stuff/Brodie Lab/Summer 2020/female preference final MS stuff/R data/female_time_on_filter_paper.csv")

t.test(filter.paper$filter, mu = 9.8)

