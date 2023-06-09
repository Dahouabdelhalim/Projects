## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R2:                                                                  ###
## Structural Equation Modeling in ethnobiology                              ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R2: Structural Equation Modeling in ethnobiology 
## Dataset needed: data_1_ethno.csv

rm(list = ls())

library(piecewiseSEM)

## Set your working directory
## setwd(/Users/Dropbox/Biological Review/Data/)

## Load data
ethnobotany<-read.csv("data_1_ethno.csv", header=T)

## Create list of structural equations
SEM_model <- piecewiseSEM::psem(
  glm(Knowledge ~ Gender + Age + Urbanization, family=binomial(link="logit"), 
      data=ethnobotany), 
  glm(Gender ~ Urbanization, family = binomial(link="logit"),data=ethnobotany)
)

# Now we extract a summary of the object. The output will look familiar to anyone who has run a regression in R.

summary(SEM_model)
