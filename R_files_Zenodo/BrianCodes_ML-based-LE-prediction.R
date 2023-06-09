#-------------------------------------------------------------------------------------------------#
## Codes: A Machine Learning Based Prediction Model for Life Expectancy
# Software: R Version 4.2.0
#-------------------------------------------------------------------------------------------------#


# Importing libraries and the Data
#--------------------------------------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(kableExtra)
library(car)
library(caret)
library(glmnet)
library(gridExtra)
library(forecast)
library(factoextra)
library(mclust)
d1 <- read.csv("Life Expectancy Data.csv")
d2 <- read_excel("Country-Metadata.xls")
d3 <- read.csv("PopData.csv", fileEncoding = 'UTF-8-BOM')
colnames(d3) <- c("Country", "CountryCode", 2000:2015)
d3 <- d3 %>% gather(key = "Year", value = "Population", 3:18)
d3$Year <- as.integer(d3$Year)
d3 <- d3 %>% select(CountryCode, everything())
d3 <- d3 %>% select(everything()) %>% group_by(Country) %>% arrange(Country, desc(Year))
# Combining the datasets
LExp <- d1 %>% 
        left_join(d2, by="Country") %>% 
        select(-SpecialNotes)%>% group_by(Country) %>% arrange(Country, desc(Year))
LExp <- subset(LExp, Country!="NA" & Country!="Nauru" & Country!="Dominica" & Country!="Tuvalu") %>%
        select(CountryCode, Country, Year, Population, everything())
# Excluding missing values in the dataset
LExp <- subset(LExp, !is.na(CountryCode))
# Replacing the population variable with data from the World population dataset
LExp$Population <- d3$Population
# Dropping the country code, status, and income composition of resources variables
LifeExp <- LExp %>% 
        select(Country,-CountryCode,Region,IncomeGroup,Year,
               -Status,Life.expectancy,Adult.Mortality,infant.deaths,
               Alcohol,percentage.expenditure,Hepatitis.B,Measles,BMI,
               under.five.deaths,Polio,Total.expenditure,Diphtheria,
               HIV.AIDS,GDP,Population,thinness..1.19.years,
               thinness.5.9.years,-Income.composition.of.resources,Schooling)
# Renaming the variables with short clear descriptive names
LifeExp <- LifeExp %>% 
        rename(IncomeG=IncomeGroup,LifeExp=Life.expectancy,AdultM=Adult.Mortality,InfantD=infant.deaths,
               PercHealthExp=percentage.expenditure,HepsB=Hepatitis.B,Und5Deaths=under.five.deaths,
               GovHealthExp=Total.expenditure,Diph=Diphtheria,HIV=HIV.AIDS,Pop=Population,
               Thin10_19Yrs=thinness..1.19.years,Thin5_9Yrs=thinness.5.9.years)
# Inspecting the structure of the dataset and formatting the data types accordingly
str(LifeExp)
LifeExp$Region <- as.factor(LifeExp$Region)
LifeExp$IncomeG <- as.factor(LifeExp$IncomeG)
LifeExp <- LifeExp %>% mutate(IncomeG=recode_factor(IncomeG, "Lower middle income"="Low income",
                                                    "Upper middle income"="Middle Income",
                                                    "High income: OECD"="High Income",
                                                    "High income: nonOECD"="High Income"))
LifeExp <- LifeExp %>% 
        mutate(IncomeG=recode_factor(IncomeG, '1'="Low Income", "Low income"="Low Income"))
levels(LifeExp$IncomeG)
levels(LifeExp$Region)
# Getting a glimpse of the data
Dglimpse <- LifeExp %>% 
        slice_sample(prop = .075)

# Exploratory Data Analysis
LEdata <- LifeExp %>% select(-LifeExp, everything())

# Check for Missingness
library(Amelia)
LEdata_NAs <- LEdata %>% select_if(function(x) any(is.na(x))) %>% 
        summarise_each(funs(sum(is.na(.))))
Missing <- sapply(LEdata, function(x) sum(is.na(x)))
LExpNAs <- data.frame(Missing)
LifeExpNAs <- data.frame(Variable= row.names(LExpNAs),
                         Missing)
row.names(LifeExpNAs) <- NULL
LifeExpNACount <- LifeExpNAs %>% filter(Missing!=0) %>% 
        arrange(desc(Missing))
xtable::xtable(LifeExpNACount,
               label = "tab:ledatamissing",
               caption = "Counts of Missing Values in the LE Data Set")
# Missingness Map
missmap(LEdata, col = c("purple", "gray"))

# Alternative Visualizations for Missingness
library(VIM)
# Aggregations for missing values
miss_plot <- aggr(LEdata, col=c('lightgray','purple'),
                  numbers=F, sortVars=TRUE, axes = T,
                  labels=names(LEdata), cex.axis=.7,
                  only.miss = T, rank.order = T, prop = T,
                  gap=3, ylab=c("Missing Data","Pattern"))

# Scatter Plots
LEdata %>% 
        drop_na() %>% 
        ggplot(mapping = aes(AdultM, LifeExp, color = IncomeG))+
        geom_point()+ geom_jitter(width = 0.15)+
        labs(y = "Life expectancy in age", x = "Number of people dying between 15-60 years per 1000 population") +
        theme_bw()

LEdata %>% 
        drop_na() %>% 
        ggplot(mapping = aes(InfantD, LifeExp, color = IncomeG))+
        geom_point()+ geom_jitter(width = 0.15)+
        labs(y = "Life expectancy in age", x = "Number of infant deaths per 1000 population") +
        theme_bw()

AggrLEdata <- LEdata %>% drop_na() %>% group_by(Region, Year) %>% summarise(LifeExp = mean(LifeExp))

LEdata %>% drop_na() %>% ggplot(aes(Year, LifeExp)) +
        geom_jitter(aes(color=IncomeG)) + 
        geom_line(data = AggrLEdata, alpha = 1, size = 1, 
                  color = "firebrick4", linetype = 1) +
        labs(x="", y="Life expectancy in age") +
        guides(col = guide_legend("Income Group"))+
        scale_x_continuous(breaks = 2000:2015) +
        theme_bw()+ facet_wrap(vars(Region))+
        theme(axis.text.x = element_text(angle = 90, hjust = 1), 
              legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom"))

# Scatter Plots of Predictors
library(psych)
LEdata[,5:21] %>% slice_sample(prop = .25) %>% 
        pairs.panels(gap=0, 
                     pch = 18, 
                     cex = 1, 
                     hist.col = 6)

# Outlier detection and removal
boxplot(LEdata$AdultM)$out #79 data points

library(EnvStats)
# Adult Mortality
OutTestAdultM <- rosnerTest(LEdata$AdultM,
                            k = 10)
OutTestAdultM$all.stats
AdultMOutliers <- LEdata[c(2828,2827,2829,350),]
Zimbab <- LEdata %>% select(Country, Year, Pop) %>% filter(Country=="Zimbabwe") %>% arrange(Year)

ZimbBotsData <- filter(LEdata, Country==c("Zimbabwe","Botswana"))

# Infant Deaths
OutTestInfantD <- rosnerTest(LEdata$InfantD,
                             k = 10)
OutTestInfantD$all.stats
InfantDOutliers <- LEdata[1127:1136,]

# Alcohol
OutTestAlcohol <- rosnerTest(LEdata$Alcohol,
                             k = 10)
OutTestAlcohol$all.stats # No outliers

# PercHealthExp
OutTestPercHealthExp <- rosnerTest(LEdata$PercHealthExp,
                                   k = 10)
OutTestPercHealthExp$all.stats
PercHealthExpOutliers <- LEdata[c(2434:2437, 1512,1509,1506,1507,1514,1878),] # Drop this variable
# since the values look erroneous for percentages.
LEdata <- LEdata %>% select(everything(),-PercHealthExp)
# Population
OutTestPop <- rosnerTest(LEdata$Pop,
                         k = 10)
OutTestPop$all.stats
PopOutliers <- LEdata[545:554,] # Not really outliers but the actual populations
# This is due to the populous nature of China, being the world's most populous country.

# Boxplots
Sub_Saharan_Deaths <- LEdata %>% group_by(Year) %>% 
        filter(Region == "Sub-Saharan Africa") %>% 
        select(Country, AdultM, InfantD)

# Adult Mortality
AdultMortalityBoxPlot <- LEdata %>% drop_na() %>% sample_frac(0.25) %>% 
        ggplot(aes(x = IncomeG, y = AdultM, fill = IncomeG)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitterdodge()) +
        facet_wrap(vars(Region)) +
        theme_bw() +
        labs(
                x = "",
                y = "Number of people dying between 15-60 years per 1000 population"
        )+
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom"))
ggsave(filename = "Adult Mortality Boxplot.pdf", plot = AdultMortalityBoxPlot)

# Government health expenditure
GovHealthExpBoxPlot <- LEdata %>% drop_na() %>% sample_frac(0.25) %>% 
        ggplot(aes(x = IncomeG, y = GovHealthExp, fill = IncomeG)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitterdodge()) +
        facet_wrap(vars(Region)) +
        theme_bw() +
        labs(
                x = "",
                y = "Government expenditure of health as a percentage of total government expenditure"
        )+
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom"))
ggsave(filename = "Government Health Expenditure Boxplot.pdf", plot = GovHealthExpBoxPlot)

# Population
LEPopData <- LEdata %>% mutate(NewPop = Pop/1000000)
PopulationBoxPlot <- LEPopData %>% drop_na() %>% sample_frac(0.25) %>% 
        ggplot(aes(x = IncomeG, y = NewPop, fill = IncomeG)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(position = position_jitterdodge()) +
        facet_wrap(vars(Region)) +
        theme_bw() +
        labs(
                x = "",
                y = "Population density (in millions)"
        )+
        theme(axis.text.x = element_text(angle = 90, hjust = 1),
              legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom"))
ggsave(filename = "Population Density Boxplot.pdf", plot = PopulationBoxPlot)

# Summary statistics of the transformed dataset
summary(LEdata)
RegSummary <- LEdata %>% 
        group_by(Region, IncomeG) %>% 
        summarise(Count = n())

xtable(RegSummary)

# Multivariate Normality Check
library(MVN)
LE_MVN_Check <- MVN::mvn(drop_na(LEdata[,5:21]), mvnTest = "royston",
                         multivariatePlot = "qq",univariateTest = "SW", showOutliers = T)
MultivariateNormality <- LE_MVN_Check$multivariateNormality
UnivariateNormality <- LE_MVN_Check$univariateNormality
names(MultivariateNormality)[c(1:4)] <- names(UnivariateNormality)[c(1,3:5)]
# Adding a column after the Test column before combining the two dataframes
MultivariateNormality <- MultivariateNormality %>% add_column(Variable = "MVN Test", 
                                                              .after = "Test")
LEdataMVNTestCheck <- rbind(UnivariateNormality, MultivariateNormality)
LEdataMVNTestCheck[,3] <- as.numeric(LEdataMVNTestCheck[,3])

LEdataMVNTestCheck <- data.frame(lapply(LEdataMVNTestCheck,
                                        function(x) if(is.numeric(x)) round(x, digits = 4) else x))
library(xtable)
xtable(LEdataMVNTestCheck,
       caption = "Univariate and Multivariate Normality Assessment Test Results",
       label = "tab:mvntest")

# MVN Density Plots
MVNTest <- mvn(drop_na(LEdata[,5:21]),
               mvnTest = "royston",
               univariateTest = "SW",
               multivariatePlot = "none",
               showOutliers = T, 
               univariatePlot = "histogram")

# Lag Plots
LagP1 <- gglagplot(LEdata[,5:10], lags = 1, do.lines = F) + theme_bw() + theme(legend.position = "none")
ggsave("Lag Plots for AdultM through Measles.pdf", plot = LagP1)

LagP2 <- gglagplot(LEdata[,11:16], lags = 1, do.lines = F)+ theme_bw() + theme(legend.position = "none")
ggsave("Lag Plots for BMI through HIV.pdf", plot = LagP2)

LagP3 <- gglagplot(LEdata[,17:21], lags = 1, do.lines = F)+ theme_bw() + theme(legend.position = "none")
ggsave("Lag Plots for GDP through LE.pdf", plot = LagP3)

# Average Life Expectancy per Region
MeanLifeExps <- LEdata %>% group_by(Region) %>% summarise(MeanLE=mean(LifeExp)) %>% 
        arrange(desc(MeanLE))
xtable::xtable(MeanLifeExps,
               caption = "Average Life Expectancies per Region",
               label = "tab:mean_le")

# Median Adult Mortality per Region
MedianAdultM <- LEdata %>% group_by(Region, IncomeG) %>% summarise(MedianAdultM=median(AdultM)) %>% 
        arrange(Region, desc(MedianAdultM))
xtable::xtable(MedianAdultM,
               caption = "Average Life Expectancies per Region",
               label = "tab:mean_le")

# Average Population Density per Region
MeanPopDensity <- LEdata %>% group_by(Region) %>% summarise(MeanPop=mean(Pop)/1000000) %>% 
        arrange(desc(MeanPop))
xtable::xtable(MeanPopDensity,
               caption = "Mean Population Density per Region",
               label = "tab:mean_popdensity")

# Data Sub-setting
LEdata1 <- LEdata %>% filter(Year !=2015) # For use in model development & testing
LEdata2 <- LEdata %>% filter(Year == 2015) # Reserved for use in confirmatory predictions


# PCA for data imputation and variable selection
library(factoextra)
library(FactoMineR)
library(missMDA)
# Imputing the Scaled LEdata dataset
ImputePCA_LEdata <- imputePCA(LEdata1[,5:21], ncp = 2) # Imputing missing values with PCA model
LEdata1Imputed <- as.data.frame(ImputePCA_LEdata$completeObs)
LEdata1Imputed <- LEdata1Imputed %>% mutate(Country = LEdata1$Country, Region = LEdata1$Region, 
                                            IncomeG = LEdata1$IncomeG, Year = LEdata1$Year) %>% 
        select(Country, Region, IncomeG, Year, everything())
# Feature Scaling
PrecProcData <- preProcess(LEdata1Imputed, 
                           method = c("center", "scale"),
                           pcaComp = 2)
ScaledLEdata <- predict(PrecProcData, LEdata1Imputed)
# ScaledLEdata <- ScaledLEdata %>% select(-PercHealthExp)
pcaLE <- PCA(ScaledLEdata[,5:20], graph = F)
get_eig(pcaLE)
# The scree plot
LEScreePlot <- fviz_screeplot(pcaLE, addlabels = T,
                              title="", 
                              barfill = 16, 
                              barcolor = 16, 
                              linecolor = 6, 
                              xlab = "Principal Components",
                              ggtheme = theme_bw())
# Contribution of Variables to PC1 plot
VarContribution <- fviz_contrib(pcaLE, 
                                choice = "var", 
                                axes = 1, 
                                sort.val = "desc",
                                title = "",
                                fill = "darkslategrey",
                                color = "darkslategrey",
                                ggtheme = theme_bw())

# Creating a data subset for model development with PCA suggested variables/features
LEdata1_Subset <- LEdata1 %>% 
        select(Region, IncomeG, Thin5_9Yrs, Thin10_19Yrs, 
               Schooling, BMI, Und5Deaths, InfantD, LifeExp)
LEdata1_Subset <- LEdata1_Subset[,-1]

LEdata2_Subset <- LEdata2 %>% 
        select(Region, IncomeG, Thin5_9Yrs, Thin10_19Yrs, 
               Schooling, BMI, Und5Deaths, InfantD, LifeExp)
LEdata2_Subset <- LEdata2_Subset[,-1]

# Model Development
#__________________
#1. XGBoost
#------------------
library(xgboost)
set.seed(234)
index <- caret::createDataPartition(LEdata1_Subset$LifeExp, p=0.8, list = F)
trainset <- LEdata1_Subset[index,]
testset <- LEdata1_Subset[-index,]

train.x <- data.matrix(trainset[,-9])
train.y <- trainset$LifeExp

test.x <- data.matrix(testset[,-9])
test.y <- testset$LifeExp

# Defining final training and testing sets
xgbTrain <- xgb.DMatrix(data = train.x, label = train.y)
xgbTest <- xgb.DMatrix(data = test.x, label = test.y)

# Defining parameters with XGBoost defaults
set.seed(234)
def_params <- list(booster = "gbtree", objective = "reg:squarederror", 
                   eta=0.3, gamma=0, max_depth=6, min_child_weight=1, 
                   subsample=1, colsample_bytree=1)

# Fitting XGBoost model & displaying training and testing data at each round,
# and finding the best iteration
set.seed(1234)
xgbCV <- xgb.cv( params = def_params, data = xgbTrain, nrounds = 100, 
                 nfold = 10, showsd = T, stratified = T, print.every.n = 10, 
                 early.stop.round = 20, maximize = F)

## Best iteration = 99
BestIteration <- xgbCV$best_iteration

set.seed(234)
startTime1 <- Sys.time()
LE_xgbModel <- xgboost(data = xgbTrain, params = def_params, 
                       nrounds = BestIteration, verbose = 0)
endTime1 <- Sys.time()

LE_xgbModel$evaluation_log

ElapsedTime <- round(endTime1 - startTime1, digits = 3)

# Predictions
set.seed(1234)
LE_Initial <- predict(LE_xgbModel, xgbTest)

MAE <- round(caret::MAE(test.y, LE_Initial), 3)
RMSE <- round(caret::RMSE(test.y, LE_Initial), 3)
RSquared <- round(caret::R2(test.y, LE_Initial), 3)
Learning_Rate <- 0.3
Max_Depth <- 6
Nrounds <- 99

LE_Initial_XGBoostErrors <- as.data.frame(cbind(Model="XGBoost Default", MAE, RMSE, 
                                                RSquared, ElapsedTime, Learning_Rate, 
                                                Max_Depth, Nrounds))

# XGBoost Hyperparameter Tuning: Model Optimization via Grid Search Cross-Validation
set.seed(234)
samples <- caret::createDataPartition(LEdata1_Subset$LifeExp, p=0.8, list = F)
training_data <- LEdata1_Subset[samples,]
testing_data <- LEdata1_Subset[-samples,]
train_x <- data.matrix(training_data[,-9])
train_y <- training_data[,9]
test_x <- data.matrix(testing_data[,-9])
test_y <- testing_data[,9]

# Modelling
set.seed(2022)
library(MachineShop)
library(doParallel)
registerDoParallel(cores = 4) # Setting 4 parallel computations for efficiency in training
startTime <- Sys.time()
xgbTunedGrid0 <- TunedModel(XGBTreeModel,
                            grid = 10,
                            control = CVControl,
                            metrics = c("rmse", "mae")) %>%
        fit(LifeExp~., data = training_data)
endTime <-Sys.time()

(as.MLModel(xgbTunedGrid0))
xgbTunedGrid0$params #eta = 0.3, gamma = 0, max_depth = 7, min_child_weight = 1, subsample = 1,
# colsample_bytree = 1, nrounds = 500, Optimization method: Grid Search

xgbFinalPred <- predict(xgbTunedGrid0, testing_data)


MAE <- round(caret::MAE(testing_data$LifeExp, xgbFinalPred), 3)
RMSE <- round(caret::RMSE(testing_data$LifeExp, xgbFinalPred), 3)
RSquared <- round(caret::R2(testing_data$LifeExp, xgbFinalPred), 3)
ElapsedTime <- round(endTime - startTime, digits = 3)
Learning_Rate <- xgbTunedGrid0$params$eta
Max_Depth <- xgbTunedGrid0$params$max_depth
Nrounds <- xgbTunedGrid0$niter

xgbFinalErrors <- as.data.frame(cbind(Model="XGBoost Tuned", MAE, RMSE, 
                                      RSquared, ElapsedTime, Learning_Rate,
                                      Max_Depth, Nrounds))

# Final XGBoost Model
n_estimators <- xgbTunedGrid0$niter
optimal_parameters <- xgbTunedGrid0$params

set.seed(2022)
startTimeF <- Sys.time()
XGBoost_LE_Model <- xgboost(data = xgbTrain, params = optimal_parameters, 
                            nrounds = n_estimators, verbose = 0)
endTimeF <- Sys.time()
ElapsedTime <- round(endTimeF - startTimeF, digits = 3)

# Final Model Predictions
XGBoost_Final_Pred <- predict(XGBoost_LE_Model, xgbTest)

MAE <- round(caret::MAE(test.y, XGBoost_Final_Pred), 3)
RMSE <- round(caret::RMSE(test.y, XGBoost_Final_Pred), 3)
# RSquared <- round(caret::R2(test.y, XGBoost_Final_Pred), 3)
# Learning_Rate <- xgbTunedGrid0$params$eta
# Max_Depth <- xgbTunedGrid0$params$max_depth
# Nrounds <- xgbTunedGrid0$niter

XGBoost_LE_Performance <- as.data.frame(cbind(Model="XGBoost", MAE, RMSE, ElapsedTime))

# _____________________
#2. Random Forest Model
# ---------------------
set.seed(2022)
# Defining the control
trControl <- trainControl(method = "cv",
                          number = 10,
                          search = "grid")
registerDoParallel(cores = 4) # Setting 4 parallel computations for efficiency in training
startTimeRF <- Sys.time()
RF_Fit0 <- train(LifeExp~., data = na.omit(training_data),
                 method = "rf",
                 metric = "RMSE",
                 trControl = trControl)
endTimeRF <-Sys.time()

RFPred <- predict(RF_Fit0, testing_data)


MAE <- round(caret::MAE(testing_data$LifeExp, RFPred), 3)
RMSE <- round(caret::RMSE(testing_data$LifeExp, RFPred), 3)
ElapsedTime <- round(endTimeRF - startTimeRF, digits = 3)

RFPredErrors <- as.data.frame(cbind(Model="Random Forest", MAE, RMSE, ElapsedTime))

# ___________________
#3. ANN
# -------------------
library(keras)
library(neuralnet)
library(magrittr)
NNdata <- LEdata1_Subset

# Imputing missing data with the median
str(NNdata)
dataLE <- na.omit(NNdata) %>% mutate_if(is.factor, as.numeric)

NnetM <- neuralnet(LifeExp~.,
                   data = dataLE,
                   hidden = 5,
                   linear.output = F,
                   lifesign = 'full',
                   rep=1)

# Model Visual Identification
plot(NnetM,col.hidden = 'firebrick4',     
     col.hidden.synapse = 'maroon',
     show.weights = F,
     information = F,
     fill = 'purple')

dataLE <- as.matrix(dataLE)
dimnames(dataLE) <- NULL

# Data Partitioning

set.seed(2022)
NNindex <- sample(2, nrow(dataLE), replace = T, prob = c(.8, .2))
training <- dataLE[NNindex==1,1:8]
test <- dataLE[NNindex==2, 1:8]
trainingtarget <- dataLE[NNindex==1, 9]
testtarget <- dataLE[NNindex==2, 9]
str(trainingtarget)
str(testtarget)

# Scaling: we normalize the values for better prediction
avg <- colMeans(training)
stdev <- apply(training, 2, sd)
training <- scale(training, center = avg, scale = stdev)
test <- scale(test, center = avg, scale = stdev)

# Model Building
library(reticulate)
library(tensorflow)
library(keras)
set.seed(2022)
DNNModel <- keras_model_sequential()
DNNModel %>% 
        layer_dense(units = 20, activation = 'relu', input_shape = c(8)) %>% 
        layer_dropout(rate=0.4) %>% 
        layer_dense(units = 1)
# We are testing one hidden layer (ANN) with 20 neurons, 8 predictor variables and 1 output neuron.
# Model compilation
DNNModel %>% compile(loss = 'mse',
                     optimizer = optimizer_rmsprop(learning_rate = 0.005), 
                     metrics = 'mae')
# Model Fitting
set.seed(2022)
startTime <- Sys.time()
DNNModel %>% 
        fit(training,trainingtarget,
            epochs = 100,
            batch_size = 32,
            validation_split = 0.2)
endTime <- Sys.time()
ElapsedTime <- round(endTime - startTime, digits = 3)

DNNModel %>% evaluate(test, testtarget) # Loss/MSE is 26.082855 or RMSE = 5.107138

set.seed(2022)
LE_ANNpredictions <- DNNModel %>% predict(test)
plot(testtarget, LE_ANNpredictions)
# Accuracy Check
MSE <- round(mean((testtarget-LE_ANNpredictions)^2), 3)
MAE <- round(caret::MAE(testtarget, LE_ANNpredictions), 3)
RMSE <- round(caret::RMSE(testtarget, LE_ANNpredictions), 3)
RSquared <- round(caret::R2(testtarget, LE_ANNpredictions), 3)

LE_ANNErrors <- as.data.frame(cbind(Model="ANN", MAE, RMSE, ElapsedTime))
# names(LE_DNNErrors)[5] <- "RSquared"

# RF, ANN, XGBoost Model Performance Comparison
Performance_Comparison <- rbind(RFPredErrors,
                                LE_ANNErrors, XGBoost_LE_Performance)
xtable(Performance_Comparison, 
       caption = "Model Performance Comparison Results",
       label = "tab:xgboost_ann_rf_comparison_results")

# Final Model Variable Importance
VarImp_Matrix = xgboost::xgb.importance(colnames(xgbTrain), model = XGBoost_LE_Model)

# Plotting the feature importance as a bar graph
library(Ckmeans.1d.dp)
library(extrafont)
loadfonts(device="win")
xgbFinalVarImportance <- xgboost::xgb.ggplot.importance(importance_matrix = VarImp_Matrix, 
                                                        measure = "Gain") +
        labs(title = "", y = "", x = "") + 
        scale_y_continuous(limits = c(0,0.8)) + 
        theme(text = element_text(family = "Times New Roman")) + 
        theme_bw()

# ------------------------------------------------------------------------------------------------
# Year 2015 Predictions:
# ------------------------------------------------------------------------------------------------
Predictions_2015 <- predict(xgbTunedGrid0, LEdata2_Subset)
LEdata2015_WithPredValues <- LEdata2_Subset %>% mutate("LifeExpPred" = Predictions_2015)
row.names(LEdata2015_WithPredValues) <- LEdata2$Country

# Plot of Predicted vs. Observed Life Expectancy Values Per Country
LE_PredPlot <- LEdata2015_WithPredValues %>% 
        ggplot(aes(x = rownames(LEdata2015_WithPredValues), group = 1)) + 
        geom_line(aes(y = LifeExp, color = "purple"), size = 0.4) + 
        geom_line(aes(y = LifeExpPred, color = "red"), size = 0.4) + 
        scale_x_discrete(guide = guide_axis(angle = 90, check.overlap = T)) + 
        scale_y_continuous(limits = c(10,100)) + 
        scale_color_identity(name = "LifeExp", 
                             breaks = c("purple", "red"), 
                             labels = c("Actual", "Predicted"),
                             guide = "legend") + 
        labs(x = "", y = "") + 
        theme(legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom")) + 
        facet_wrap(vars(Region)) + 
        theme_bw()

# Mean of Predicted vs. Observed Life Expectancy Values Per Region
AvgPred_ActualLE <- LEdata2015_WithPredValues %>% select(Region, LifeExp, LifeExpPred) %>% 
        group_by(Region) %>% summarise("Mean Actual LE" = mean(LifeExp), 
                                       "Mean Predicted LE" = mean(LifeExpPred))

MeanLE_PredPlot <- AvgPred_ActualLE %>% 
        ggplot(aes(x = Region, group = 1)) + 
        geom_line(aes(y = `Mean Actual LE`, color = "purple"), size = 0.4) + 
        geom_line(aes(y = `Mean Predicted LE`, color = "red"), size = 0.4) + 
        scale_x_discrete(guide = guide_axis(angle = 90, check.overlap = T)) + 
        scale_y_continuous(limits = c(10,100)) + 
        scale_color_identity(name = "Mean LifeExp", 
                             breaks = c("purple", "red"), 
                             labels = c("Actual", "Predicted"),
                             guide = "legend") + 
        labs(x = "", y = "") + 
        theme(legend.position = c(.75,0.01), 
              legend.justification = c("right", "bottom")) +
        theme_bw()

#_______________________________________________________________________________
# END OF CODE
#-------------------------------------------------------------------------------