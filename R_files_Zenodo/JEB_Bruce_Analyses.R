# JEB Bruce Distractor Data Analyses 

setwd('/Users/awinsor/JEB_Bruce/Data') # data file path, folder where csv files are located in R proj

### Experiment 1 ====
# multiple primary, constant distractor, compare to null
library(readr) # to read csv files 
Exp1.data <- read_csv("BruceExpt1JEB.csv") # importing data file into R and naming the dataframe
summary(Exp1.data) # gives a quick summary of the data

# NOTE: Looked Away Y/N is 1/0 in the csv file 

# Remove N/As
# always remove largest to smallest to prevent shifting data frame  
NEWExp1.data <- Exp1.data[-42,]

# should have lost 1 row 
library(MASS)
n.samples <- nrow(NEWExp1.data) 
n.samples # 47 trials 

library(tidyverse)
# Bar graph 
NEWExp1.data %>% 
  ggplot(aes(Stimulus, fill=factor(Shifted))) +
  geom_bar() # if you don't want stacked bars you can add position = "dodge" into the parentheses in front of geom_bar
              
# 0 is did not look at distractor, 1 is looked away to distractor
# primary stimulus type is on the x axis
# looks like crickets are too interesting to look away from 

# Fit a model 
library(lme4)

# glmer is GLMM in lme4 package, required since you have repeat IDs
# y ~ x + (random effect: each ID has same slope of 1 but different intercepts), data=yourdata, family=binomial for binary
Mod1 <- glmer (Shifted ~ Stimulus + (1|ID), data=NEWExp1.data, family="binomial")
summary(Mod1) # give model output
# ignore intercept significance, it just tells you the result is different than 1 

#make sure order is categorical 
NEWExp1.data$Order <- as.factor(NEWExp1.data$Order)
#drop NAs for order
NEWExp1.data2 <- NEWExp1.data[-26,]
NEWExp1.data3 <- NEWExp1.data2[-10,]

#check stimulus order 
Mod1_order <- glmer (Shifted ~ Stimulus + Order + (1|ID), data=NEWExp1.data3, family="binomial")
summary(Mod1_order)
#order not significant 

#check order graphically
NEWExp1.data3 %>% 
  ggplot(aes(Order, fill=factor(Shifted))) +
  geom_bar(position="dodge")
# 0 is did not look at distractor, 1 is looked away to distractor

# check model assumptions 
par(mfrow=c(2,1)) # 2 column, 1 row plot

#residuals from model w/ GLMM
clr <- adjustcolor(4,0.3)
plot(resid(Mod1), main="Mod1", pch=16, col=clr) 
abline(h=0,lwd=2) # shows distribution of residuals 
qqnorm(resid(Mod1), main="Mod1", pch=16, col=clr) 
qqline(resid(Mod1)) # shows normailty of residuals

# check mixed model fit
library(DHARMa)
simulateResiduals(Mod1)
plot(simulateResiduals(Mod1))
# looks good 

# look at estimates 
# Mod1NoIntercept <- glmer (Shifted ~ Stimulus + (1|ID) +-1, data=NEWExp1.data, family="binomial") # +-1 removes intercept
# summary(Mod1NoIntercept) #ignore significance here 

# you can back transform using inverse logit  
#library(boot)
#probability1 <- inv.logit(c(Blank, Cricket, Oval)) # c bind lets you combine things into a vector 
# multiply by 100 to get percents 
#percent_prob <- 100*probability1
#percent_prob

# Calculate upper and lower CI using Standard errors 
# formula: estimate ± (percentile × SE of estimate)
# percentile is 1.96 for 95% CI
# you can't back transform chunks (estimates and SE separately) and then add later, all must be done one the same scale 

dev.off() # kills graphics 

### Experiment 2 ====
# constant primary, different distractor
library(readr)
Exp2.data <- read_csv("BruceExpt2JEB.csv") 
summary(Exp2.data) 

# NOTE: LookedAway/DidNotLookAway is 1/0 in the csv file
  # This is to keep things consistent with exp 1, BUT
  # keep track of this because the variable "stimulus" changes from representing primary to distractor across exp 1 and 2

# 2 N/As needed to be removed (123,116)
# protocol failed for 116, outlier removed for 123 
NEWExp2.data <- Exp2.data[-123,]
NEWExp2.data2 <- NEWExp2.data[-116,]
NEWExp2.data2 #check only NAs removed

library(tidyverse)
# Bar graph 
NEWExp2.data2 %>% 
  ggplot(aes(Stim, fill=factor(Response))) +
  geom_bar(position="dodge") 
# 0 is did not look at distractor, 1 is looked away to distractor
# distractor type is shown on x axis
# looks like loom causes them to look away from primary to distractor more

#make sure order is categorical 
NEWExp2.data2$Order <- as.factor(NEWExp2.data2$Order)

# Fit a model 
library(lme4)

Mod2 <- glmer (Response ~ Stim + (1|ID), data=NEWExp2.data2, family="binomial")
summary(Mod2)
#loom highly significant

#check order 
Mod2_order <- glmer (Response ~ Stim + Order+ (1|ID), data=NEWExp2.data2, family="binomial")
summary(Mod2_order)
# order trending to significance in order 4
# AIC is all within 2 points of each other 

#drop1 and pairwise
drop1(Mod2_order, test="Chisq")
#stim is significant, order is not
library(multcomp)
summary(glht(Mod2, linfct=mcp(Stim="Tukey"))) #loom still significant with Tukey test
#to get pairwise
summary(glht(Mod2_order, linfct=mcp(Stim="Tukey")))
# loom still significant 

# check residual patterns
library(DHARMa)
simulateResiduals(Mod2)
plot(simulateResiduals(Mod2))
# looks good 

# these values will negligibly change depending if analysis is run with 0 or 1 for 287 cross 
dev.off() 





