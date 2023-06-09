#04 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publication


### Paired with
# "01_Artemia_E.csv"

#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(car)        # needed for type II ANOVA command 'Anova()'
library(multcomp)  # needed for contrasts



setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in


#_______________________________________ Artemia _________________________________________
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ Univariate $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


vitEart <- read.csv("01_Artemia_E.csv", header = TRUE)



ArtE.glm <- glm(Concentration ~ Enrichment, data = vitEart, family = Gamma)
summary(ArtE.glm)

Anova(ArtE.glm)




#.................................Checking assumptions ............................

plot(ArtE.glm)
plot(as.factor(vitEart$Enrichment), resid(ArtE.glm, type = "pearson"))
qqnorm(resid(ArtE.glm))
qqline(resid(ArtE.glm))
hist(residuals(ArtE.glm), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)




#----------- Contrasts:
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
ArtE.glm.cont <- glm(Concentration ~ Enrichment - 1, data = vitEart, family = Gamma)
summary(ArtE.glm.cont)

Cont <- matrix(c(-1, -1, 3, -1), 1)
cont.test <- glht(ArtE.glm.cont, linfct = Cont)
summary(cont.test)




