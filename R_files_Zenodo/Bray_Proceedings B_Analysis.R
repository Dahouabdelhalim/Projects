#Title: Incomplete control and concessions explain mating skew in male chimpanzees
#Authors: Joel Bray, Anne E. Pusey, Ian C. Gilby
#Journal: Proceedings of the Royal Society B
#Contact: bray.joel@gmail.com

###########################################
#Setup workspace and load image & libraries
###########################################

#clear workspace
rm(list=ls())

#set workspace to the location of "ProceedingsB_Data.RData"
setwd("~/Downloads") 

#load "ProceedingsB_Data.RData"
load("Bray_ProceedingsB_Data.RData") #note: all model outputs are pre-run and saved

#view saved objects in workspace
ls()

#"Bray_ProceedingsB_Data.RData" includes the following objects:
#analysis.1_data: data table for the Analysis 1 GLMM
#analysis.2_data: data table for the Analysis 2 GLMM
#analysis.1_output: model output for Analysis 1 GLMM
#analysis.2_output: model output for Analysis 2 GLMM
#.Rdata image also includes likelihood-ratio tests for both models and pairwise comparisons where applicable

#to rerun the statistical analyses, install and load the following libraries:
#install.packages("lme4")
#install.packages("multcomp")
library(lme4)
library(multcomp)

###########################################################
#ANALYSIS 1: Factors influencing alpha male mating monopoly
###########################################################

#Generalized Linear Mixed Model (see Table 1)

#Variables Names and Definitions in "analysis.1_data":

#outcome = copulating male was the alpha male (1/0)
#adult.males = number of males present in the party
#swollen.parous = number of maximally swollen parous females present in the party
#alphaID = alpha male identity
#alpha.related = whether the female was closely related (mother, maternal sibling, or maternal niece) to the alpha male (1/0)
#alpha.subject = whether the alpha male was the focal subject (1/0)
#femaleID = female identity
#maleID = male identity

#Other columns in "analysis.1_data"

#date = date of mating
#time = time of mating
#focal = identity of individual who was being followed

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))

analysis.1_output = glmer(formula = outcome ~ adult.males + swollen.parous + alphaID + alpha.related + alpha.subject + (1|femaleID) + (1|maleID), data = analysis.1_data, family = binomial(link="logit"), control=contr)

#view model output for Analysis 1
summary(analysis.1_output)

#Likelihood ratio tests ('Chisq' is equivalent in the drop1 function)
analysis.1_LRT = as.data.frame(drop1(analysis.1_output, test="Chisq"))

#Tukey contrasts of the "alpha male identity" factor
#see Electronic Supplementary Material: Table S2
analysis.1_tukey = summary(glht(analysis.1_output, mcp(alphaID="Tukey")))



################################################################
#Analysis 2: Factors influencing subordinate male mating success
################################################################

#Generalized Linear Mixed Model (see Table 2)

#Variables Names and Definitions:

#outcome = whether a particular non-alpha male copulated (1/0)
#male.groom = male grooming index with the alpha (Z-transformed)
#male.rank = male dominance rank (Elo Rating) (Z-transformed)
#male.related = whether the non-alpha male was closely related to the female (mother, maternal sibling, or maternal niece) (1/0)
#male.age.bin = male age (categorized into 6 age bins)
#male.party.size = number of males present in the party (Z-transformed)
#alphaID = alpha male identity
#maleID = non-alpha male identity
#femaleID = female identity
#eventID= = mating event ID

#Other columns in data file: analysis.2_data

#date = date of mating
#time = time of mating
#focal = identity of individual who was being followed

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1000000))

analysis.2_output = glmer(formula = outcome ~ male.groom + male.rank*male.party.size + male.related + male.age.bin + (1|alphaID) + (1|maleID) + (1|femaleID) + (1|eventID) + (0+male.groom|alphaID) + (0+male.rank|alphaID), data = analysis.2_data, family = binomial(link = "logit"), control=contr)

#view model output for Analysis 2
summary(analysis.2_output) 

#Likelihood ratio tests ('Chisq' is equivalent in the drop1 function)
analysis.2_LRT = as.data.frame(drop1(analysis.2_output, test="Chisq"))

#Tukey contrasts of the "male age bins" factor
#see Electronic Supplementary Material: Table S3
analysis.2_tukey = summary(glht(analysis.2_output, mcp(male.age.bin="Tukey")))



#End