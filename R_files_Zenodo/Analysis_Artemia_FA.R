#04 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publciation


### Paired with

#Tank means
# "01_Artemia_Neu_FA.csv"
# "02_Artemia_Phos_FA.csv"





#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(vegan)      # needed for adnois (PERMANOVA) function
library(car)        # needed for type II ANOVA command 'Anova()'
#library(multcomp)  # needed for contrasts - Don't run until after using any 'select' commands since attaches package 'MASS', which masks 'select'





setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in



#...................................................................................................
# --------------------------------------- Multivariate Analyses ----------------------------------------------
#...................................................................................................




# --------------------------------------- Neutral ----------------------------------------------

artemia_neu_bad <- read.csv("01_Artemia_Neu_FA.csv", header = FALSE)

artemia_neu <- data.frame(t(artemia_neu_bad[-1]))
colnames(artemia_neu) <- artemia_neu_bad[, 1]


Neu.FA <- matrix(as.numeric(as.character(unlist(artemia_neu[3:20]))), nrow = nrow(artemia_neu))
colnames(Neu.FA) <- colnames(artemia_neu[3:20])





#------------------------------- PERMANOVA ------------------------------------


model.neu <- adonis(Neu.FA ~ artemia_neu$Treatment, permutations = 4999, method="bray")
model.neu







#------------------------------- SIMPER ------------------------------------

simper.test <- simper(Neu.FA, group = artemia_neu$Treatment, permutations=999)
simper.test  #Show most influential species to pairwise differences
summary(simper.test)  #Show most influential overall


















# --------------------------------------- Phospholipid ----------------------------------------------


artemia_phos_bad <- read.csv("02_Artemia_Phos_FA.csv", header = FALSE)

artemia_phos <- data.frame(t(artemia_phos_bad[-1]))
colnames(artemia_phos) <- artemia_phos_bad[, 1]


Phos.FA <- matrix(as.numeric(as.character(unlist(artemia_phos[, 3:20]))), nrow = nrow(artemia_phos))
colnames(Phos.FA) <- colnames(artemia_phos)[3:20]




#------------------------------- PERMANOVA ------------------------------------


model.phos <- adonis(Phos.FA ~ artemia_phos$Treatment, permutations = 4999, method = "bray")
model.phos








#------------------------------- SIMPER ------------------------------------


simper.test <- simper(Phos.FA, group = artemia_phos$Treatment, permutations=999)
simper.test  #Show most influential species to pairwise differences
summary(simper.test)  #Show most influential overall

















#...................................................................................................
# --------------------------------------- Univariate Analyses ----------------------------------------------
#...................................................................................................




UniArt.neu <- artemia_neu %>%
  select(Treatment, C18.1n9, EPA, DHA)


UniArt.phos <- artemia_phos %>%
  select(Treatment, C18.1n9, EPA, DHA)





# --------------------------------------- Neutral ----------------------------------------------

### -------------------- EPA ----------------------


EPA.glm.neu <- glm(as.numeric(EPA)/100 ~ Treatment, data = UniArt.neu, family = quasibinomial(link = "logit"))
summary(EPA.glm.neu)
 
Anova(EPA.glm.neu)



#.................................Checking assumptions ............................

plot(EPA.glm.neu)
qqnorm(resid(EPA.glm.neu))
qqline(resid(EPA.glm.neu))
hist(residuals(EPA.glm.neu), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)





#----------- Contrasts
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
EPA.glm.neu.cont <- glm(as.numeric(EPA)/100 ~ Treatment - 1, data = UniArt.neu, family = quasibinomial(link = "logit"))
summary(EPA.glm.neu.cont)

library(multcomp)
Cont <- matrix(c(1, -1, -1, 1), 1) #PUFA treatments vs. control and unenriched
cont.test.neu <- glht(EPA.glm.neu.cont, linfct = Cont)
summary(cont.test.neu)










### -------------------- DHA ----------------------


DHA.glm.neu <- glm(as.numeric(DHA)/100 ~ Treatment, data = UniArt.neu, family = quasibinomial(link = "logit"))
summary(DHA.glm.neu)






#.................................Checking assumptions ............................
plot(DHA.glm.neu)
qqnorm(resid(DHA.glm.neu))
qqline(resid(DHA.glm.neu))
hist(residuals(DHA.glm.neu), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)




#----------- Contrasts
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
DHA.glm.neu.cont <- glm(as.numeric(DHA)/100 ~ Treatment - 1, data = UniArt.neu, family = quasibinomial(link = "logit"))
summary(DHA.glm.neu.cont)

library(multcomp)
Cont <- matrix(c(1, -1, -1, 1), 1)  #PUFA treatments vs. control and unenriched
cont.test <- glht(DHA.glm.neu.cont, linfct = Cont)
summary(cont.test)














# --------------------------------------- Phospholipid ----------------------------------------------

### -------------------- EPA ----------------------


EPA.glm.phos <- glm(as.numeric(EPA)/100 ~ Treatment, data = UniArt.phos, family = quasibinomial(link = "logit"))
summary(EPA.glm.phos)

Anova(EPA.glm.phos)


#.................................Checking assumptions ............................
plot(EPA.glm.phos)
qqnorm(resid(EPA.glm.phos))
qqline(resid(EPA.glm.phos))
hist(residuals(EPA.glm.phos), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)



#----------- Contrasts
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept
EPA.glm.phos.cont <- glm(as.numeric(EPA)/100 ~ Treatment - 1, data = UniArt.phos, family = quasibinomial(link = "logit"))
summary(EPA.glm.phos.cont)


library(multcomp)
Cont <- matrix(c(1, -1, -1, 1), 1)  #PUFA treatments vs. control and unenriched
cont.test.phos <- glht(EPA.glm.phos.cont, linfct = Cont)
summary(cont.test.phos)






### -------------------- DHA ----------------------


DHA.glm.phos <- glm(as.numeric(DHA)/100 ~ Treatment, data = UniArt.phos, family = quasibinomial(link = "logit"))
summary(DHA.glm.phos)

Anova(DHA.glm.phos)



#.................................Checking assumptions ............................
plot(DHA.glm.phos)
plot(as.factor(UniArt.phos$Treatment), resid(DHA.glm.phos, type = "pearson"))
qqnorm(resid(DHA.glm.phos))
qqline(resid(DHA.glm.phos))
hist(residuals(DHA.glm.phos), xlab = "Standardized residuals", ylab = "Frequency", main = NULL)

### FAIL - does not meet the assumptions of the model














