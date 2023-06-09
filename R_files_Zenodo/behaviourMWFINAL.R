library(tidyverse)
library(qcc)
library(mvabund)
library (multcomp)
library(emmeans)
library(lme4)
library(ggplot2)
library (plyr)
library(reshape2)
library(lmerTest)

# ------------------------------------------------------------------------------------------------------
###  DOES BEHAVIOUR OF FISH (BOTH SPP) VARY DEPENDING ON GROUP TYPE?

results = read.csv ("TotalBehaviours.csv")
str (results)

# MEAN AGGRESSIVE
qqnorm (results$meanaggressive)
qqline (results$meanaggressive)
hist (results$meanaggressive) 
results$logmeanaggressive.1 = log(results$meanaggressive.1) 
hist (results$logmeanaggressive) 

# ~ Shoal Type


ft.meanaggressive2 = lmer(logmeanaggressive.1 ~ species*shoaltype + (1|groupid), data = results) #based on feedback from sci report
anova(ft.meanaggressive2) 

#calculate means
means1 = by(results$logmeanaggressive.1, list(results$shoaltype, results$species), mean)

# calculate SE
se = function (x) sqrt ( var (x) / length (x))
ses1 = by(results$logmeanaggressive.1, list(results$shoaltype, results$species), se)

#create plots of Species and Shoal type (x value, trace value, response value)
interaction.plot(results$species,results$shoaltype,results$logmeanaggressive.1, type="b", pch=c(21,19), ylim=c(0.2,2.2), las=1, ylab="Mean count of aggressive behaviours (per focal fish)", xlab="Species", legend=FALSE)
legend("bottomleft",c("Single","Mixed"),lty=c(1,2))
lines(c(1,1), c(means1[1]-ses1[1], means1[1]+ses1[1]))
lines(c(1,1), c(means1[2]-ses1[2], means1[2]+ses1[2]))
lines(c(2,2), c(means1[3]-ses1[3], means1[3]+ses1[3]))
lines(c(2,2), c(means1[4]-ses1[4], means1[4]+ses1[4]))


# Tukeys - check what was significant
tukeysaggressive <- lsmeans(ft.meanaggressive2, pairwise ~ species*shoaltype, adjust = "Tukey")
tukeysaggressive
plot(tukeysaggressive) 

#-------------------------------------------------------------------------------------------------------
# MEAN SUBMISSIVE
qqnorm (results$meansubmissive)
qqline (results$meansubmissive)
hist (results$meansubmissive) 

results$logmeansubmissive = log(results$meansubmissive.1) 
hist (results$logmeansubmissive) 

# ~ Shoal Type

ft.meansubmissive2 = lmer(logmeansubmissive ~ species*shoaltype +(1|groupid), data = results) 
anova(ft.meansubmissive2) 

#calculate means
means2 = by(results$logmeansubmissive, list(results$shoaltype, results$species), mean)

# calculate SE
se = function (x) sqrt ( var (x) / length (x))
ses2 = by(results$logmeansubmissive, list(results$shoaltype, results$species), se)

#create plots of Species and Shoal type (x value, trace value, response value)
interaction.plot(results$species,results$shoaltype,results$logmeansubmissive, type="b", pch=c(21,19), ylim=c(0,2.3), las=1, ylab="Mean count of submissive behaviours (per focal fish)", xlab="Species", legend=FALSE)
legend("bottomleft",c("Single","Mixed"),lty=c(1,2))
lines(c(1,1), c(means2[1]-ses2[1], means2[1]+ses2[1]))
lines(c(1,1), c(means2[2]-ses2[2], means2[2]+ses2[2]))
lines(c(2,2), c(means2[3]-ses2[3], means2[3]+ses2[3]))
lines(c(2,2), c(means2[4]-ses2[4], means2[4]+ses2[4]))


# Tukeys 
tukeyssubmissive <- lsmeans(ft.meansubmissive2, pairwise ~ species*shoaltype, adjust = "Tukey")
tukeyssubmissive
plot(tukeyssubmissive) 

#-------------------------------------------------------------------------------------------------------
# MEAN SOCIAL
qqnorm (results$meansocial)
qqline (results$meansocial)
hist (results$meansocial) 
results$logmeansocial = log(results$meansocial)
hist(results$logmeansocial) 
str(results)

# ~ Shoal Type



ft.logmeansocial2 = lmer(logmeansocial ~ species*shoaltype + (1|groupid), data = results)
anova (ft.logmeansocial2) 

#calculate means
means3 = by(results$meansocial, list(results$shoaltype, results$species), mean)

# calculate SE
se = function (x) sqrt ( var (x) / length (x))
ses3 = by(results$meansocial, list(results$shoaltype, results$species), se)

#create plots of Species and Shoal type (x value, trace value, response value)
interaction.plot(results$species,results$shoaltype,results$meansocial, type="b", pch=c(21,19), ylim=c(5,14), las=1, ylab="Mean count of social behaviours (per focal fish)", xlab="Species", legend=FALSE)
legend("bottomleft",c("Single","Mixed"),lty=c(1,2))
lines(c(1,1), c(means3[1]-ses3[1], means3[1]+ses3[1]))
lines(c(1,1), c(means3[2]-ses3[2], means3[2]+ses3[2]))
lines(c(2,2), c(means3[3]-ses3[3], means3[3]+ses3[3]))
lines(c(2,2), c(means3[4]-ses3[4], means3[4]+ses3[4]))

# Tukeys 
tukeyssocial <- lsmeans(ft.logmeansocial2, pairwise ~ species*shoaltype, adjust = "Tukey")
tukeyssocial
plot(tukeyssocial) 

#-------------------------------------------------------------------------------------------------------
# MEAN FORAGING
qqnorm (results$meanforaging)
qqline (results$meanforaging)
hist (results$meanforaging) 
results$logmeanforaging = log(results$meanforaging.1) 
hist (results$logmeanforaging) 

# ~ Shoal Type

ft.meanforaging2 = lmer(logmeanforaging ~ species*shoaltype + (1|groupid), data = results) # MW added log here
anova (ft.meanforaging2) 

#calculate means
means4 = by(results$logmeanforaging, list(results$shoaltype, results$species), mean)

# calculate SE
se = function (x) sqrt ( var (x) / length (x))
ses4 = by(results$logmeanforaging, list(results$shoaltype, results$species), se)

#create plots of Species and Shoal type (x value, trace value, response value)
interaction.plot(results$species,results$shoaltype,results$logmeanforaging, type="b", pch=c(21,19), ylim=c(0.4,0.8), las=1, ylab="Mean count of foraging interactions (per focal fish)", xlab="Species", legend=FALSE)
legend("bottomleft",c("Single","Mixed"),lty=c(1,2))
lines(c(1,1), c(means4[1]-ses4[1], means4[1]+ses4[1]))
lines(c(1,1), c(means4[2]-ses4[2], means4[2]+ses4[2]))
lines(c(2,2), c(means4[3]-ses4[3], means4[3]+ses4[3]))
lines(c(2,2), c(means4[4]-ses4[4], means4[4]+ses4[4]))

# Tukeys -
tukeysforaging <- lsmeans(ft.meanforaging2, pairwise ~ species*shoaltype, adjust = "Tukey")
tukeysforaging
plot(tukeysforaging) 
