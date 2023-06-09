library(tidyverse)
library(qcc)
library(mvabund)
library (multcomp)
library(lsmeans)
library(lme4)
library(ggplot2)
# ------------------------------------------------------------------------------------------------------
###   Does shoaling vary depending on group type?

results2 = read.csv("Total_Shoal_new.csv", stringsAsFactors = TRUE)
str(results2)

#with dispersion as binary
ft.dispersion2 = glm(meandispersionbinary ~ species*shoaltype, family = "binomial", data = results2)
anova(ft.dispersion2, test = "Chisq") 
plot(ft.dispersion2)

interaction.plot(results2$species,results2$shoaltype,results2$meandispersion)


#calculate means
means6 = by(results$meandispersion, list(results$shoaltype, results$species), mean)

# calculate SE
se = function (x) sqrt ( var (x) / length (x))
ses6 = by(results2$meandispersion, list(results2$shoaltype, results2$species), se)

#create plots of Species and Shoal type (x value, trace value, response value)
interaction.plot(results2$species,results2$shoaltype,results2$meandispersion, type="b", pch=c(21,19), ylim=c(0.4,1.3), las=1, ylab="Mean cohesion (per shoal)", xlab="Species", legend=FALSE)
legend("bottomleft",c("Single","Mixed"),lty=c(1,2))
lines(c(1,1), c(means6[1]-ses6[1], means6[1]+ses6[1]))
lines(c(1,1), c(means6[2]-ses6[2], means6[2]+ses6[2]))
lines(c(2,2), c(means6[3]-ses6[3], means6[3]+ses6[3]))
lines(c(2,2), c(means6[4]-ses6[4], means6[4]+ses6[4]))

# Tukeys - check what was significant
tukeysdispersion <- lsmeans(ft.dispersion, pairwise ~ species*shoaltype, adjust = "Tukey")
tukeysdispersion
plot(tukeysdispersion) 
