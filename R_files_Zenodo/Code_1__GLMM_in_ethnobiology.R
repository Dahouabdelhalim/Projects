## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R1:                                                                  ###
## Generalized Linear Mixed Effect Model (GLMM) in ethnobiology              ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R1: Generalized Linear Mixed Effect Model (GLMM) in ethnobiology
## Dataset needed: data_1_ethno.csv

## Packages to load
## We will use the package glmmTMB to implement the generalized linear mixed effects models. The other packages are for data import, manipulation and graphing purposes.
library(readr)
library(dplyr)
library(tibble)
library(stringr)
library(broom)
library(broom.mixed)
library(glmmTMB)
library(bbmle)
library(car)
library(MuMIn)
library(ggplot2)
library(GGally)

## Set your working directory
## setwd(/Users/Dropbox/Biological Review/Data/)

## Importing data
db <- read_csv('data_1_ethno.csv', col_types = cols(Gender = col_character()))

## Examining data
dim(db) ## [1] 144   9
glimpse(db) ## Rows: 144 Columns: 9

summary(db)
names(db)

## 1. Building a traditional generalized model (GLM)
## Here, we build a simple additive model with age, gender and the education level of participants as predictors of the participant’s ability to know local plants importance.
model0 <- glm(Knowledge ~ Age + Gender + Education_level,
              family = 'binomial', data = db)

summary(model0)

## From the output of the model, the knowledge of plants importance by participants depends mostly on his age and gender. The chance to have a participant with a knowledge on plants importance increases as age increases. Besides, men are more likely to know plants importance than women.

## 2. First Scenario: Building generalized linear mixed effects models (GLMM)
## For the first example of a generalized mixed effects model to test this hypothesis, we keep a simple additive model with age, gender and the education level as fixed terms. The random terms are the participant nested in his village nested in the region.
model1 <- glmmTMB(Knowledge ~ Age + Gender + Education_level
                  + (1|Region/Village/Subject), 
                  family = 'binomial',
                  data = db)

summary(model1)

## Fixed effects
## Globally, we can observe that the output of the fixed effects section of this GLMM model yield similar conclusions to the GLM model. The coefficients obtained for the GLMM model are a bit higher than that of the traditional GLM.

## Random effects
## The random effects part of the model output present each term with the variance and standard deviation associated. On the total variation of random components, the region and the village contribute up to 99.9% to the variation observed in the ability of the participant to know plants importance. Furthermore, the village explained more this variation than the region of study.

print(VarCorr(model1), comp = c("Std.Dev.","Variance"))

## Now Let's explore the individual random effects for each term.
rand.coef1 <- ranef(model1, condVar = TRUE)

rand.svr1 <- rand.coef1[[1]][[1]]
rand.vr1 <- rand.coef1[[1]][[2]]
rand.r1 <- rand.coef1[[1]][[3]]

d1 <- bind_rows(rand.svr1, rand.vr1, rand.r1)%>%
  rownames_to_column(., "rand.term")
colnames(d1) <- c("rand.term", "estimate")
## View(d1)

## Subject in Village in Region
d2 <- d1 %>%
  head(143) %>%
  mutate(. , level = rep("Region/Village/Subject"))

## View(d2)

## Illustrative figures
(graph.svr1 <- ggplot(d2, aes(rand.term, estimate))+
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dotted", 
               color = "black", size = 1)+
    coord_flip()+
    theme_bw() +
    theme(axis.text = element_text(size = 10.5, color = "grey10"),
          axis.title = element_text(size = 11))
)

## Village in Region
d3 <- d1[c(144, 145, 146, 147, 148, 149), ] %>% 
  mutate(. , level = rep("Region/Village"))
# View(d3)

(graph.vr1 <- ggplot(d3, aes(rand.term, estimate))+
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dotted", 
               color = "black", size = 1) +
    coord_flip()+
    theme_bw() +
    theme(axis.text = element_text(size = 10.5, color = "grey10"),
          axis.title = element_text(size = 11))
)

## Region
d4 <- d1 %>%
  tail(3) %>%
  mutate(. , level = rep("Region"))
# View(d4)

(graph.r1 <- ggplot(d4, aes(rand.term, estimate))+
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dotted", 
               color = "black", size = 1) +
    coord_flip()+
    theme_bw() +
    theme(axis.text = element_text(size = 10.5, color = "grey10"),
          axis.title = element_text(size = 11))
)

## 3. Second scenario: Building a GLMM, with region as fixed effect
## Simple additive model with Age, Gender, Education_level and Region as fixed terms; with random terms as Subject nested in Village.
model2 <- glmmTMB(Knowledge ~ Age + Gender + Education_level
                  + Region + (1|Village/Subject), 
                  family = 'binomial',
                  data = db)

summary(model2)

## Fixed effects
## The outputs of the model show that the education level has no significant link with the participant ability to know plants importance. Conversely, the chance of having some knowledge of plants importance is related to participants’ age, gender and region. Participants with some knowledge are more likely to be women and from the second and third region. The probability to have a participant with some knowledge on plants importance increases with aging as well.

## Random effects
## In this scenario, the village alone is enough to explain about 99% of the total variation observed in the ability of the participant to know plants importance.
print(VarCorr(model2), comp = c("Std.Dev.","Variance"))

Subject.Village1 <- 100 * (3.5770e-08 /(3.5770e-08 + 3.0096e-01))
Village1 <- 100 * (3.0096e-01 /(3.5770e-08 + 3.0096e-01))

rand.prop2 <- bind_cols(Subject.Village1, Village1)

colnames(rand.prop2)<- c("SV1", "V1")
rand.prop2

## R-square
## It is also possible to check how well the model built fits the data by calculating the coefficient of determination R-squared. Using the package MuMIn, this can be obtained with the functions rsquaredLR for the traditional GLM and rsquaredGLMM for GLMM. Note that r.squaredGLMM can also be used for a traditional GLM.

## r.squaredGLMM produces two outputs: R2m and R2c. R2m represents the variance explained by the fixed effects while R2c can be interpreted as variance explained by both fixed and random components. delta and theoretical are two different methods implemented to calculate R-squared.
r.squaredLR(model0)
r.squaredGLMM(model0)
r.squaredGLMM(model1)
r.squaredGLMM(model2)

## For the GLMM in the first scenario, fixed effects alone explain about 56-59% of the variation observed. Accounting for the random effects increases this explanation to 66-70%. The variation explained by the GLMM is lower than that of the traditional model (54-58%). The trend is similar with the GLMM in the 2nd scenario which also explains better the variation observed in the data. Thus, including random effects improved the model.
