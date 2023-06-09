# 2/3/2022
# v4.1.1


# Cleaned up for publication


### Paired with

#Phospholipid values
# "FA_37d_Neu_01.csv"

#Tank means
# "Juvenile_end-point1_tank_totals.csv"



library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(vegan)      # needed for adnois (PERMANOVA) function
library(car)        # needed for type II ANOVA command 'Anova()'
#library(multcomp)  # needed for contrasts - Don't run until after using any 'select' commands since attaches package 'MASS', which masks 'select'


 


setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in


#...................................................................................................
### #################################### Phase 2 Neutral FA #########################################
#...................................................................................................



Neu.bad <- read.csv("FA_37d_Neu_01.csv", header = FALSE)

Neu <- data.frame(t(Neu.bad[-1]))
colnames(Neu) <- Neu.bad[, 1]

Neu$Tank <- as.factor(Neu$Tank)

d37.density.data <- read.csv("Juvenile_end-point1_tank_totals.csv", header = TRUE) %>%
  mutate(density = Initial.density) %>%
  select(Tank, density) %>%
  mutate(Tank = as.factor(Tank))


Neu.new <- full_join(Neu, d37.density.data, by = "Tank")

Neu.FA <- matrix(as.numeric(as.character(unlist(Neu.new[3:18]))), nrow = nrow(Neu.new))
colnames(Neu.FA) <- colnames(Neu.new[3:18])





### --------------- PERMANOVA - Neutral -----------------


model.neu <- adonis(Neu.FA ~ Neu.new$Treatment + Neu.new$density, permutations = 4999, method="bray")
model.neu






#determining if there is homogenetiy of dispersion
dist.mat <- vegdist(Neu.FA, method = "bray")
beta.test <- betadisper(dist.mat, Neu.new$Treatment)
perm.beta.test <- permutest(beta.test, pairwise = FALSE, permutations = 4999)
perm.beta.test

                    








# ----------------Don't need to use a SIMPER if PERMANOVA not significant -----------
### --------------- SIMPER - Neutral -----------------

simper.test <- simper(Neu.FA, group = Neu.new$Treatment, permutations=4999)
simper.test  #Show most influential species to pairwise differences
summary(simper.test)  #Show most influential overall








### ---------------------------- Univariate - Neutral -----------------------------------


Neu.uni <- Neu.new %>%
  select(Treatment, Tank, density, C18.1n9, C20.5, C22.6)


Avg.Neu.uni <- Neu.uni %>%
  group_by(Treatment) %>%
  summarize(Avg.Avg.EPA = mean(as.numeric(C20.5)), 
            Avg.Avg.DHA = mean(as.numeric(C22.6)))





### -------------------- EPA ----------------------


## GLM for EPA - First converting EPA percentage back to proportion by dividing by 100
EPA.glm.neu <- glm(as.numeric(C20.5)/100 ~ Treatment + density, data = Neu.uni, family = quasibinomial(link = "logit"))
summary(EPA.glm.neu)

Anova(EPA.glm.neu)



#.................................Checking assumptions ............................

plot(EPA.glm.neu)
qqnorm(resid(EPA.glm.neu))
qqline(resid(EPA.glm.neu))
hist(residuals(EPA.glm.neu), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)





#No Need to do contrasts if no significance in glm











### -------------------- DHA ----------------------


## GLM for DHA - First converting DHA percentage back to proportion by dividing by 100
DHA.glm.neu <- glm(as.numeric(C22.6)/100 ~ Treatment + density, data = Neu.uni, family = quasibinomial(link = "logit"))
summary(DHA.glm.neu)

Anova(DHA.glm.neu)




#.................................Checking assumptions ............................

plot(DHA.glm.neu)
qqnorm(resid(DHA.glm.neu))
qqline(resid(DHA.glm.neu))
hist(residuals(DHA.glm.neu), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)






#----------- Contrasts: https://stats.idre.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
DHA.glm.neu.cont <- glm(as.numeric(C22.6)/100 ~ Treatment + density - 1, data = Neu.uni, family = quasibinomial(link = "logit"))

library(multcomp)
Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(DHA.glm.neu.cont, linfct = Cont1)
summary(cont.test1)
Cont2 <- matrix(c(0, 1, -1, 0), 1)  #PUFA vs. PUFA + a-t
cont.test2 <- glht(DHA.glm.neu.cont, linfct = Cont2)
summary(cont.test2)














