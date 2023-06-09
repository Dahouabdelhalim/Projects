#03 Feb 2022
#V4.1.1
#---------


# Cleaned up for publication


### Paired with

#Tank means
# "01_vitE_10d_37d.csv"
# "Larval_end-point_tank_totals.csv"
# "Juvenile_end-point1_tank_totals.csv"



#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(car)        # needed for type II ANOVA command 'Anova()'
#library(multcomp)  # needed for contrasts - Don't run until after using any 'select' commands since attaches package 'MASS', which masks 'select'



setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in



#...................................................................................................
### #################################### Data Wrangling #########################################
#...................................................................................................



vitE <- read.csv("01_vitE_10d_37d.csv", header = TRUE)

Larv.tanks <- read.csv("Larval_end-point_tank_totals.csv", header = TRUE)

Density.d10 <- Larv.tanks %>%
  select(Tank, Trial.Stocking) %>%
  mutate(Phase1.density = Trial.Stocking/50) %>%
  select(-Trial.Stocking)


Juv.tanks <- read.csv("Juvenile_end-point1_tank_totals.csv", header = TRUE)

Density.d37 <- Juv.tanks %>%
  select(Tank, Initial.density) %>%
  mutate(Phase2.density = Initial.density) %>%
  select(-Initial.density)



both.uni <- vitE %>%
  full_join(Density.d10, by = "Tank") %>%
  full_join(Density.d37, by = "Tank")



#...................................................................................................
### #################################### Phase 1 vit E #########################################
#...................................................................................................


vitE.uni <- both.uni %>%
  group_by(Treatment) %>%
  summarize(Avg.E.P1 = mean(as.numeric(as.character(Phase1E))))




Phase1E.glm <- glm(Phase1E ~ Treatment + Phase1.density, data = both.uni, family = gaussian)
summary(Phase1E.glm)

Anova(Phase1E.glm)



#.................................Checking assumptions ............................

plot(Phase1E.glm)
qqnorm(resid(Phase1E.glm))
qqline(resid(Phase1E.glm))
hist(residuals(Phase1E.glm), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)




#----------- Contrasts: https://stats.idre.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
Phase1E.glm.cont <- glm(Phase1E ~ Treatment + Phase1.density - 1, data = both.uni, family = gaussian)
summary(Phase1E.glm.cont)

library(multcomp)
Cont1 <- matrix(c(-1, -1, 2, 0), 1) #PUFA+E vs. Control and PUFA
cont.test1 <- glht(Phase1E.glm.cont, linfct = Cont1)
summary(cont.test1)
Cont2 <- matrix(c(-1, 1, 0, 0), 1) #Control vs. PUFA
cont.test2 <- glht(Phase1E.glm.cont, linfct = Cont2)
summary(cont.test2)















#...................................................................................................
### #################################### Phase 2 vit E #########################################
#...................................................................................................


vitE.uni <- both.uni %>%
  group_by(Treatment) %>%
  summarize(Avg.E.P2 = mean(as.numeric(as.character(Phase2E))))




Phase2E.glm <- glm(Phase2E ~ Treatment + Phase2.density, data = both.uni, family = gaussian)
summary(Phase2E.glm)

Anova(Phase2E.glm)






#.................................Checking assumptions ............................

plot(Phase2E.glm)
qqnorm(resid(Phase2E.glm))
qqline(resid(Phase2E.glm))
hist(residuals(Phase2E.glm), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)











