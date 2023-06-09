# 1/28/2022
# v4.1.1


# Cleaned up for publication


### Paired with

#Phospholipid values
# "FA_10d_Phos_01.csv"

#Tank means
# "Larval_end-point_tank_totals.csv"



#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(vegan)      # needed for adnois (PERMANOVA) function
library(car)        # needed for type II ANOVA command 'Anova()'
#library(multcomp)  # needed for contrasts - Don't run until after using any 'select' commands since attaches package 'MASS', which masks 'select'


 


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in


#...................................................................................................
### #################################### Phase 1 Phospholipid FA #########################################
#...................................................................................................



Phos.bad <- read.csv("FA_10d_Phos_01.csv", header = FALSE)

Phos <- data.frame(t(Phos.bad[-1]))
colnames(Phos) <- Phos.bad[, 1]

Phos$Tank <- as.factor(Phos$Tank)

d10.density.data <- read.csv("Larval_end-point_tank_totals.csv", header = TRUE) %>%
  mutate(density = Trial.Stocking/50) %>%
  select(Tank, density) %>%
  mutate(Tank = as.factor(Tank))


Phos.new <- full_join(Phos, d10.density.data, by = "Tank")

Phos.FA <- matrix(as.numeric(as.character(unlist(Phos.new[3:19]))), nrow = nrow(Phos.new))
colnames(Phos.FA) <- colnames(Phos.new[3:19])





### --------------- PERMANOVA - Phospholipid -----------------


model.Phos <- adonis(Phos.FA ~ Phos.new$Treatment + Phos.new$density, permutations = 4999, method="bray")
model.Phos





#determining if there is homogenetiy of dispersion
dist.mat <- vegdist(Phos.FA, method = "bray")
beta.test <- betadisper(dist.mat, Phos.new$Treatment)
perm.beta.test <- permutest(beta.test, pairwise = FALSE, permutations = 4999)
perm.beta.test









### --------------- SIMPER - Phospholipid -----------------

simper.test <- simper(Phos.FA, group = Phos.new$Treatment, permutations=4999)
simper.test  #Show most influential species to pairwise differences
summary(simper.test)  #Show most influential overall








### ---------------------------- Univariate - Phospholipid -----------------------------------


Phos.uni <- Phos.new %>%
  select(Treatment, Tank, density, C18.1n9, C20.5, C22.6)


Avg.Phos.uni <- Phos.uni %>%
  group_by(Treatment) %>%
  summarize(Avg.Avg.EPA = mean(as.numeric(C20.5)), 
            Avg.Avg.DHA = mean(as.numeric(C22.6)))





### -------------------- EPA ----------------------


## GLM for EPA - First converting EPA percentage back to proportion by dividing by 100
EPA.glm.Phos <- glm(as.numeric(C20.5)/100 ~ Treatment + density, data = Phos.uni, family = quasibinomial(link = "logit"))
summary(EPA.glm.Phos)

Anova(EPA.glm.Phos)



#.................................Checking assumptions ............................

plot(EPA.glm.Phos)
qqnorm(resid(EPA.glm.Phos))
qqline(resid(EPA.glm.Phos))
hist(residuals(EPA.glm.Phos), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)






#----------- Contrasts: https://stats.idre.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
EPA.glm.Phos.cont <- glm(as.numeric(C20.5)/100 ~ Treatment + density - 1, data = Phos.uni, family = quasibinomial(link = "logit"))

library(multcomp)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(EPA.glm.Phos.cont, linfct = Cont1)
summary(cont.test1)
Cont2 <- matrix(c(0, 1, -1, 0), 1)  #PUFA vs. PUFA + a-t
cont.test2 <- glht(EPA.glm.Phos.cont, linfct = Cont2)
summary(cont.test2)













### -------------------- DHA ----------------------


## GLM for DHA - First converting DHA percentage back to proportion by dividing by 100
DHA.glm.Phos <- glm(as.numeric(C22.6)/100 ~ Treatment + density, data = Phos.uni, family = quasibinomial(link = "logit"))
summary(DHA.glm.Phos)

Anova(DHA.glm.Phos)



#.................................Checking assumptions ............................

plot(DHA.glm.Phos)
qqnorm(resid(DHA.glm.Phos))
qqline(resid(DHA.glm.Phos))
hist(residuals(DHA.glm.Phos), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)






#----------- Contrasts: https://stats.idre.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
DHA.glm.Phos.cont <- glm(as.numeric(C22.6)/100 ~ Treatment + density - 1, data = Phos.uni, family = quasibinomial(link = "logit"))

library(multcomp)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(DHA.glm.Phos.cont, linfct = Cont1)
summary(cont.test1)
Cont2 <- matrix(c(0, 1, -1, 0), 1)  #PUFA vs. PUFA + a-t
cont.test2 <- glht(DHA.glm.Phos.cont, linfct = Cont2)
summary(cont.test2)






















