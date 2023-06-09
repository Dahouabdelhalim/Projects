
# Stocker et al. 2019: Cooperation with closely bonded individuals reduces cortisol levels in long-tailed macaques
# Q3: ARE CHANGES IN CORTISOL CONNECTED WITH THE TEST CONDITIONS?
#--------------------------------------------------------------------

setwd("your_working_directory")

LTM <- read.csv("Stocker_Data_Nov2019.csv")
str(LTM)

# delete rows in which delta.2.1 has no values
LTM_noZ <- LTM[!is.na(LTM$delta.2.1), ]
str(LTM_noZ)


#---------------------
# plot
par(mfrow=c(1,1))
plot(LTM_noZ$delta.2.1 ~ LTM_noZ$Condition, ylab = "Change in cortisol [ng/ml]", xlab = "Condition")


#---------------------
# check distribution
hist(LTM_noZ$delta.2.1)

library(fitdistrplus)
fit.norm<-fitdist(LTM_noZ$delta.2.1, "norm")
plot(fit.norm)


#---------------------
# full model

library(lme4)
summary(m.Q4 <- lmer(delta.2.1 ~ Condition + Sex.indiv + Rank.indiv + Baby +
                      (1|Individual), data=LTM_noZ, REML=FALSE))

library(car)
vif(m.Q4)   # no mulit-collinearity issue


#standardize model
library(arm)
m.Q4_stdz <- standardize(m.Q4)


# Model averaging
library(MuMIn)
options(na.action=na.fail)
all.m.Q4 <- dredge(m.Q4_stdz, rank="AICc") #all model combinations
all.m.Q4
sub.m.Q4 <- subset(all.m.Q4, delta<2) # model combinations less than 2 AIC difference to the 'best' models
sub.m.Q4
avg.m.Q4 <- model.avg(sub.m.Q4)
summary(avg.m.Q4)
importance(avg.m.Q4)

confint(avg.m.Q4)
