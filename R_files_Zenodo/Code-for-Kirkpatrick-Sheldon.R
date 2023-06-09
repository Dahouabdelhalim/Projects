# Data analysis for Kirkpatrick & Sheldon 
# Updated on 25 April 2022 

library(lme4)
library(languageR)
library(AICcmodavg)
library(MuMIn)
library(car)
library(ggpubr)

# Load brood ball data by first saving the .xls supporting information file as a .csv
beetle<-read.csv("Broodball-data-Phanaeus-Kirkpatrick-Sheldon.csv", header=TRUE)

##########################################################
#### Question: Are broodballs in greenhouses deeper? #####
##########################################################

# Examine brood ball depth with fixed effects of bucket, female mass, and trial and the random effect of female identification.

# First decide on random variable structure. 
depth.a <- lmer (bb.depth ~ type + female.mass + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
depth.b <- lmer (bb.depth ~ type + female.mass + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AICc(depth.a)
AICc(depth.b) # The model with the lowest AICc is depth.b, so use random slope and intercept for female identification.

# Models to examine brood ball depth with fixed effects of bucket, female mass, and trial, and the random effect of female identification.

depth.d1 <- lmer (bb.depth ~ type + female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
depth.d2 <- lmer (bb.depth ~ type + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.d3 <- lmer (bb.depth ~ type + female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.d4 <- lmer (bb.depth ~ female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.d5 <- lmer (bb.depth ~ female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.d6 <- lmer (bb.depth ~ trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.d7 <- lmer (bb.depth ~ type + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
depth.dnull <- lmer (bb.depth ~ 1 + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(depth.d1, depth.d2, depth.d3, depth.d4, depth.d5, depth.d6, depth.d7, depth.dnull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
# Round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# The model with type of treatment and trial (depth.d2) has the lowest AICc and is the best model. 

# Set REML=T to get summary data for the best fit model.
depth.d2 <- lmer (bb.depth ~ type + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(depth.d2)
r.squaredGLMM(depth.d2)

# Check that residuals follow assumpions of normality.
hist(residuals(depth.d2))
boxplot(residuals(depth.d2))
qqPlot(residuals(depth.d2))
ggqqplot(residuals(depth.d2))

# Also examine second level residuals (i.e. residuals for random intercepts and random slopes)
ranef(depth.d2) # slope is changing across females
qqPlot(ranef(depth.d2)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F) 
qqPlot(ranef(depth.d2)$female[,2], ylab="Random-slope residuals among females", xlab="Norm quantiles", id=F) 

######################################################################
# Question: Are broodballs deeper in response to warmer temperatures?#
######################################################################

# Examining brood ball depth as a function of the temperatures breeding females experienced at the soil surface. 

# Decide on random variable structure. 
surface.a <- lmer (bb.depth ~ mean.surface + female.mass + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
surface.b <- lmer (bb.depth ~ mean.surface + female.mass + trial + (mean.surface|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AICc(surface.a)
AICc(surface.b) # The model with the lowest AICc is surface.b

# Models (below) to examine brood ball depth with fixed effects of mean surface temperature of the soil, female mass, and trial, and the random effect of female identification.

surface.b1 <- lmer (bb.depth ~ mean.surface + female.mass + trial + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
surface.b2 <- lmer (bb.depth ~ mean.surface + trial + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
surface.b3 <- lmer (bb.depth ~ female.mass + trial + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
surface.b4 <- lmer (bb.depth ~ female.mass + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
surface.b5 <- lmer (bb.depth ~ trial + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
surface.b6 <- lmer (bb.depth ~ mean.surface  + female.mass + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
surface.b7 <- lmer (bb.depth ~ mean.surface + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
surface.bnull <- lmer (bb.depth ~ 1 + (mean.surface|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(surface.b1, surface.b2, surface.b3, surface.b4, surface.b5, surface.b6, surface.b7, surface.bnull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
## round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# Model surface.b7 has the lowest AICc and is the best model.

# Set REML=T to get the summary data for the best fit model.
surface.b7 <- lmer (bb.depth ~ mean.surface + (mean.surface|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(surface.b7)
r.squaredGLMM(surface.b7)

# Check that residuals follow assumpions of normality.
hist(residuals(surface.b7))
boxplot(residuals(surface.b7))
qqPlot(residuals(surface.b7))
ggqqplot(residuals(surface.b7))

# Also examine second level residuals (i.e. residuals for random intercepts and random slopes)
ranef(surface.b7)
qqPlot(ranef(surface.b7)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F)
qqPlot(ranef(surface.b7)$female[,2], ylab="Random-slope residuals among females", xlab="Norm quantiles", id=F)

###############################################################
#### Question: Does brood ball size change in greenhouses? ####
###############################################################

# Examining brood ball mass with the fixed effects of type of bucket, female mass, and trial and the random effect of female identification.

# First decide on random variable structure.
size.a <- lmer (bb.mass ~ type + female.mass + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
size.b <- lmer (bb.mass ~ type + female.mass + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AICc(size.a)
AICc(size.b) # The model with the lowest AICc  and also does not have an error of boundary (singular) fit) is size.a, so use random intercept but not random slope for female.

size.a1 <- lmer (bb.mass ~ type + female.mass + trial + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
size.a2 <- lmer (bb.mass ~ type + trial + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
size.a3 <- lmer (bb.mass ~ female.mass + trial + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
size.a4 <- lmer (bb.mass ~ type + female.mass + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
size.a5 <- lmer (bb.mass ~ trial + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
size.a6 <- lmer (bb.mass ~ type + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
size.a7 <- lmer (bb.mass ~ female.mass + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
size.anull <- lmer (bb.mass ~ 1 + (1|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(size.a1, size.a2, size.a3, size.a4, size.a5, size.a6, size.a7, size.anull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# The model with type of treatment and trial is best supported. Set REML=T to get the summary data for the best fit model.

size.a2 <- lmer (bb.mass ~ type + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(size.a2)
r.squaredGLMM(size.a2)

# Check that residuals follow assumptions of normality.
hist(residuals(size.a2))
boxplot(residuals(size.a2))
qqPlot(residuals(size.a2))
ggqqplot(residuals(size.a2))

# Also examine second level residuals (i.e. residuals for random intercept)
ranef(size.a2)
qqPlot(ranef(size.a2)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F)

#######################################################################################
#### Question: Were broodballs in greenhouses placed in warmer mean temperatures? #####
#######################################################################################

# Examining mean temperatures experienced by broodballs with fixed effects of type of bucket, female mass, and trial and the random effect of female identification.

# First decide on random variable structure. 
mean.a <- lmer (mean.temp ~ type + female.mass + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
mean.b <- lmer (mean.temp ~ type + female.mass + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AICc(mean.a)
AICc(mean.b) # The model with the lowest AICc is mean.b, so use random intercept and slope for female.

mean.b1 <- lmer (mean.temp ~ type + female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
mean.b2 <- lmer (mean.temp ~ type + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.b3 <- lmer (mean.temp ~ type + female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.b4 <- lmer (mean.temp ~ female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.b5 <- lmer (mean.temp ~ female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.b6 <- lmer (mean.temp ~ trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.b7 <- lmer (mean.temp ~ type + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
mean.bnull <- lmer (mean.temp ~ 1 + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(mean.b1, mean.b2, mean.b3, mean.b4, mean.b5, mean.b6, mean.b7, mean.bnull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# Model mean.b1 has the lowest AICc and is the best model 

# Get summary data
mean.b1 <- lmer (mean.temp ~ type + female.mass + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mean.b1) 
r.squaredGLMM(mean.b1)

# Check that residuals follow assumpions of normality.
hist(residuals(mean.b1))
boxplot(residuals(mean.b1))
qqPlot(residuals(mean.b1))
ggqqplot(residuals(mean.b1))

# Also examine second level residuals (i.e. residuals for random intercepts and random slopes)
ranef(mean.b1)
qqPlot(ranef(mean.b1)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F)
qqPlot(ranef(mean.b1)$female[,2], ylab="Random-slope residuals among females", xlab="Norm quantiles", id=F)

#########################################################################################
#### Question: Were broodballs in greenhouses placed in more variable temperatures? #####
#########################################################################################

# Examining temperature variance (as standard deviation of temperatures) experienced by brood balls with the fixed effects of type of bucket, female mass, and trial and the random effect of female identification.

# First decide on random variable structure.
stdev.a <- lmer (stdev.temp ~ type + female.mass + trial + (1|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
stdev.b <- lmer (stdev.temp ~ type + female.mass + trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
AICc(stdev.a)
AICc(stdev.b)
# The model with the lowest AICc is stdev.b, so we will use random slope and intercept for female.

stdev.b1 <- lmer (stdev.temp ~ type + female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
stdev.b2 <- lmer (stdev.temp ~ type + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.b3 <- lmer (stdev.temp ~ type + female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.b4 <- lmer (stdev.temp ~ female.mass + trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.b5 <- lmer (stdev.temp ~ female.mass + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.b6 <- lmer (stdev.temp ~ trial + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.b7 <- lmer (stdev.temp ~ type + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
stdev.bnull <- lmer (stdev.temp ~ 1 + (type|female), data=beetle, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(stdev.b1, stdev.b2, stdev.b3, stdev.b4, stdev.b5, stdev.b6, stdev.b7, stdev.bnull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
# round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# Model stdev.b6, which has only the fixed effect of trial, had the lowest AICc, however, it was only marginally better than the null model.

stdev.b6 <- lmer (stdev.temp ~ trial + (type|female), data=beetle, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(stdev.b6)  
r.squaredGLMM(stdev.b6)

# Check that residuals follow assumptions of normality.
hist(residuals(stdev.b6))
boxplot(residuals(stdev.b6))
qqPlot(residuals(stdev.b6))
ggqqplot(residuals(stdev.b6))

# Also examine second level residuals (i.e. residuals for random intercepts and random slopes)
ranef(stdev.b6)
qqPlot(ranef(stdev.b6)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F)
qqPlot(ranef(stdev.b6)$female[,2], ylab="Random-slope residuals among females", xlab="Norm quantiles", id=F)

##########################################################
#### Question: Do number of brood balls change in greenhouses? #####
##########################################################

# Load brood ball number data by first saving the .xls supporting information file as a .csv
number<-read.csv("Broodball-number-Phanaeus-Kirkpatrick-Sheldon.csv", header=TRUE)

# Examine brood ball number with fixed effects of bucket type, female mass, and trial and random effect of female identification.

# First decide on random variable structure. 
number.a <- lmer (bb.number ~ type + female.mass + trial + (1|female), data=number, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
number.b <- lmer (bb.number ~ type + female.mass + trial + (type|female), data=number, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
# Model number.a does not give an error or warning, so use random intercept but not random slope for female.

number.a1 <- lmer (bb.number ~ type + female.mass + trial + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
number.a2 <- lmer (bb.number ~ type + trial + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
number.a3 <- lmer (bb.number ~ female.mass + trial + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
number.a4 <- lmer (bb.number ~ type + female.mass + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
number.a5 <- lmer (bb.number ~ trial + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
number.a6 <- lmer (bb.number ~ type + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
number.a7 <- lmer (bb.number ~ female.mass + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
number.anull <- lmer (bb.number ~ 1 + (1|female), data=number, REML=FALSE, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))

Cand.models <- list(number.a1, number.a2, number.a3, number.a4, number.a5, number.a6, number.a7, number.anull)
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")
## round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

# Set REML=T to get the summary data for the best fit model.
number.a2 <- lmer (bb.number ~ type + trial + (1|female), data=number, REML=T, control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(number.a2)
r.squaredGLMM(number.a2)

# Check that residuals follow assumptions of normality.
hist(residuals(number.a2))
boxplot(residuals(number.a2))
qqPlot(residuals(number.a2))
ggqqplot(residuals(number.a2))

# Also examine second level residuals (i.e. residuals for random intercept)
ranef(number.a2)
qqPlot(ranef(number.a2)$female[,1], ylab="Random-intercept residuals among females", xlab="Norm quantiles", id=F)

# Subsetting data to examine average # of brood balls produced in greenhouse and control  
green <- subset(number, type=="Greenhouse")
control <- subset(number, type=="Control")
mean(control$bb.number)
mean(green$bb.number)



