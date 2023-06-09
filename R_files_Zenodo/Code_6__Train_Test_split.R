## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
## CODE R6:                                                                  ###
## Train/test split and cross validation in ethnobiology                     ###
## Gaoue et al. (2021) Biological Reviews, in press                          ###
## ogaoue@utk.edu | February 16, 2021                                        ###
## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

## This script to reproduce the analyses from:
## Gaoue, O. G.,J. K. Moutouama , M A. Coe, M. O. Bond, E. Green, N. B. Sero, B. S. Bezeng, and K. Yessoufou 2021. “Methodological advances for hypothesis-driven ethnobotany” Biological Review, in press.

## Code R6: Train/test split and cross validation in ethnobiology  
## Dataset 5: Data used in this example is available at: https://archive.ics.uci.edu/ml/datasets/census+income

## UCI_adult.csv

## Objective: To conduct Train/test split and cross validation and calculate model comparison statistics such as Mean Absolute Percentage Error (MAPE). 

rm(list=ls())

## Load packages
library(tidyverse)
library(caret)
library(MLmetrics)

## Set your working directory
## setwd(/Users/Dropbox/Biological Reviews/Data/)

## 1. Data import
data = read.csv("UCI_adult.csv", header = FALSE) #load data
income <- data[ c(1,5,10,15) ] # select out a few columns that are similar to ethnobiology data
colnames(income) <- c("age","education", "sex", "income")
income$income = as.factor(as.numeric(as.factor(income$income))) # convert to binary variable: 1 = below 50k, 2 = above 50k

## 2. EXAMPLE 1: predict whether a person's income is above or below $50k annually. In ethnobiology data, the response variable could be whether or not a person uses a particular plant

## Split the data into training and test set
set.seed(123)  #make data reproducable for this example
training.samples <- income$income %>%  createDataPartition(p = 0.8, list = FALSE) 
## splitting 80% of data into training set, remaining 20% into test set. Adjust p if you want a 90-10 or 75-25 split. 
train.data  <- income[training.samples, ]
test.data <- income[-training.samples, ]

train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 1) 
## number sets 10 splits of the data, repeats sets how many times the cross validation is repeated

cv_model <- train(income ~., data = train.data, method = "glm", trControl = train.control, family = "binomial") # Train the model
print(cv_model)  # Summarize the results
cv_model$finalModel# show model coefficients and AIC

cv_predictions <- cv_model %>% predict(test.data) # Make predictions
cv_predictions = as.numeric(cv_predictions)
test.data$income = as.numeric(test.data$income)

data.frame( R2 = R2(cv_predictions, test.data$income),
            RMSE = RMSE(cv_predictions, test.data$income),
            MAE = MAE(cv_predictions, test.data$income),
            MAPE = MAPE(cv_predictions, test.data$income))

## To compare models, the goal is higher r2 and lower RMSE, MAE, and MAPE 
## The Root Mean Squared Error (RMSE), Mean Absolute Error (MAE), Mean Absolute Percentage Error (MAPE)
## MAPE is used to compare models, even if models contain different variables
## Interpretation- for test data (that the model was not trained on): R2- the model explains 13% of variation, RMSE- the square root of the average of squared differences between predictions and test data is 0.45, MAE- our predictions are off by 0.21 on average (MAE and RMSE are easier to interpret for count or numerical data), MAPE- the model is 13% wrong (for every 100 test data points, the model gets 13 wrong). These results aren’t bad for only having 3 predictor variables!  


## 3. EXAMPLE 2: predict how many grades of education a person has. In ethnobiology data, the response variable could be how many plants a person knows for a certain use
training.samples2 <- income$income %>%  createDataPartition(p = 0.8, list = FALSE) 
train.data2  <- income[training.samples2, ]
test.data2 <- income[-training.samples2, ]

train.control <- trainControl(method = "repeatedcv", number = 10, repeats = 1)  
cv_model <- train(education ~., data = train.data2, method = "glm", trControl = train.control, family = "poisson") # Train the model
print(cv_model)  # Summarize the results
cv_model$finalModel# show model coefficients and AIC


cv_predictions <- cv_model %>% predict(test.data2) # Make predictions
data.frame( R2 = R2(cv_predictions, test.data2$education),
            RMSE = RMSE(cv_predictions, test.data2$education),
            MAE = MAE(cv_predictions, test.data2$education),
            MAPE = MAPE(cv_predictions, test.data2$education))

## Interpretation- for test data (that the model was not trained on): R2- the model explains 11% of variation, RMSE- the square root of the average of squared differences between predictions and test data is 2.45 grades, MAE- our predictions are off by 1.8 grades on average, MAPE- the model is 25% wrong (for every 4 test data points, the model gets 1 wrong) .

## The second model has a higher MAPE, but that's because it had to predict more possible outcomes (16 grades)
