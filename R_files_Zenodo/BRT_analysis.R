# ---------Tuning of BRT parameters---------
rm(list = ls())

library(rsample) 
library(dismo)
library(Metrics)
library(tidyverse)

# reading Source_data containing Stemflow percentage and the biotic and abiotic factors
GSF<-read.table("Source_data.csv", header = T, sep=",") 

# create hyperparameter grid
hyper_grid <- expand.grid(
  learning.rate = c(0.001, 0.005, 0.01, 0.05),
  tree.complexity = seq(1, 5, by = 1),
  bag.fraction = seq(0.5, 0.75, 0.05)
)

# total number of combinations
nrow(hyper_grid)
ncol(hyper_grid)

# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  # reproducibility
  set.seed(0817)
  
  # train model
  gsf_brt <- gbm.step(data = GSF,
                      gbm.x = 13:24,
                      gbm.y = 12,
                      family = "gaussian",
                      tree.complexity = hyper_grid$tree.complexity[i],
                      learning.rate = hyper_grid$learning.rate[i],
                      bag.fraction = hyper_grid$bag.fraction[i]) 
  
  # calculating some parameters associated with model performance
  hyper_grid$n_trees[i] <- list(gsf_brt$n.trees) #number of trees
  hyper_grid$cv_deviance[i] <- list(gsf_brt$cv.statistics["deviance.mean"]) #cross-validated (CV) deviance
  hyper_grid$cv_correlation[i] <- list(gsf_brt$cv.statistics["correlation.mean"])  #calculating mean CV correlation
  hyper_grid$train_correlation[i] <- list(gsf_brt$self.statistics["correlation"])  #calculating training data correlation
  hyper_grid$rmse[i] <- list(rmse(GSF$Stemflow, gsf_brt$fitted)) #calculating RMSE
  hyper_grid$MSE[i] <- list(mse(GSF$Stemflow, gsf_brt$fitted)) #calculating mean square error (MSE)
  hyper_grid$MAE[i] <- list(mae(GSF$Stemflow, gsf_brt$fitted)) #calculating mean absolute error (MAE)
}

# show results
hyper_grid

# save results
tuning <- as.matrix(tuning)
write.csv(tuning, "tuning.csv")

# #--------------Stemflow percentage BRT analysis---------------
# loading package
library(dismo)

# reading our Source_data containing Stemflow percentage and the biotic and abiotic factors
GSF<-read.table("Source_data.csv", header = T, sep = ",") 

# set seed to guarantee the consistency of output 
set.seed(0817)

# running brt model with gbm.step() function
gsfp01<-gbm.step(data = GSF,
                  gbm.x = 13:24,
                  gbm.y = 12,
                  family = "gaussian",
                  tree.complexity = 5,
                  learning.rate = 0.01,
                  bag.fraction = 0.75) 

# output for relative influence
summary(gsfp01, xlim = c(0, 30), cex.axis = 1.5, cex.lab = 1.5) 

#generate partial dependence plots
gbm.plot(gsfp01, n.plots = 12, rug = TRUE, smooth = TRUE, write.title = FALSE, common.scale = TRUE, plot.layout = c(3, 4))  
