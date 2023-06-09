# Required libraries
library(lme4)
library(lmerTest)

# Read data or use file.choose()
perimeter <- read.csv("Perimeter.csv")
perimeter$age <- factor(perimeter$age)
perimeter$cycle <- factor(perimeter$cycle)

# Fit models
modelDay <- lmer(cauliflower ~ age * day + (1|individual), data=perimeter, REML=TRUE)
anova(modelDay)
modelCycle <- lmer(cauliflower ~ age * cycle + (1|individual), data=perimeter, REML=TRUE)
anova(modelCycle)
