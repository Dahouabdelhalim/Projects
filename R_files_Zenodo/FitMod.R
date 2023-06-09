################################################################################
###                                                                          ###
###----- NAME: --------------------------------------------------------------###
###                                                                          ###
### FitMod.R                                                                 ###
###                                                                          ###
###----- AIM: ---------------------------------------------------------------###
### Function to fit a non-linear model against one explanatory variable. Can ###
### fit up to 9 different models with minpack.lm::nlsLM. Estimates           ###
### prediction and confidence intervals with propagate::predcitNLS.          ###
###                                                                          ###
###----- INPUT: -------------------------------------------------------------###
###                                                                          ###
###   1. "DataIn"  - Main input dataframe. Should contain two column         ###
###         named 'Dep' and 'Exp' for the dependent and explanatory          ###
###         variables, respectively.                                         ###
###   2. "ModName" - Name of the model to fit:                               ###
###     'Linear'               :-: Y = a1 + a2 * X                           ###
###     'Poly2'                :-: Y = a1 + a2 * X + a3 * X^2                ###
###     'Poly3'                :-: Y = a1 + a2 * X + a3 * X^2 + a4 * X^3     ###
###     'Power'                :-: Y = a1 * X^a2                             ###
###     'Logarithmic'          :-: Y = a1 + a2 * log(X)                      ###
###     'NegativeExponential'  :-: Y = a1*[1 - exp(-a2 * X)]                 ###
###     'Rational'             :-: Y = (a1 + a2 * X) / (1 + a3 * x)          ###
###     'sigEmax'              :-: Y = a1 + a2*[X^a3 / (a4^a3 + X^a3)]       ###
###   3. "Dep"     - Name of the dependent variable.                         ###
###   4. "Exp"     - Name of the explanatory variable.                       ###
###   5. "NumIntervals" - Number of points to use for the prediction and     ###
###                       confidence intervals.                              ###
###   6. "Alpha"  - The alpha level for the prediction and confidence        ###
###                 intervals.                                               ###
###                                                                          ###
###----- PROCESS: -----------------------------------------------------------###
###                                                                          ###
### 1.  Fit the models using non-linear regression:                          ###
###    1.1.  Linear model.                                                   ###
###    1.2.  Polynomial second order model.                                  ###                         
###    1.3.  Polynomial third order model.                                   ###
###    1.4.  Power model.                                                    ###
###    1.5.  Logarithmic model.                                              ###
###    1.6.  Negative Exponential model.                                     ###
###    1.7.  Rational model.                                                 ###
###    1.8.  sigEmax model.                                                  ###
### 2.  Check if the model converged.                                        ###
### 3.  Start model info results dataframe.                                  ###
### 4.  If the model converged:                                              ###
###    4.1. Get prediction and confidence intervals:                         ###
###      4.1.1.  Create a sequence along the range of Exp.                   ###
###      4.1.2.  Fit the prediction intervals.                               ###
###      4.1.3.  Fit the confidence intervals.                               ###
###      4.1.4.  Arrange prediction and confidence intervals dataframe.      ###
###    4.2.  Add coefficient values to Info.                                 ###
###    4.3.  Predict for the input data and calculate performance indices.   ###
###    4.4.  Prepare the predicted dataframe.                                ###
###    4.5.  Test residuals normality with Shapiro-Wilk test.                ###
### 5.  If the model did not converged- return NA for most values.           ###
### 6.  Arrange the return list.                                             ###
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### A list with the following slots:                                         ###
###                                                                          ###
### 1. "Model"     - the non-linear regression model.                        ###
###                                                                          ###
### 2. "Info"      - dataframe with the following columns:                   ###  
###     'Dep'           - Name of dependent.                                 ###          
###     'Exp'           - Name of explanatory.                               ###          
###     'ModName'       - Name of the model.                                 ###      
###     'Exp___ModName' - Exp and model name, separated by "___".            ###
###     'Formula'       - Formula of the model.                              ###
###     'Converged'     - Logical, T if the model converged, F if not.       ###
###     'NumCase'       - Number of cases.                                   ###
###     'a1'            - First coefficient, NA otherwise.                   ###
###     'a2'            - Second coefficient, NA otherwise.                  ###
###     'a3'            - Third coefficient, NA otherwise.                   ###
###     'a4'            - Fourth coefficient, NA otherwise.                  ###
###     'Deviance'      - The deviance.                                      ###
###     'RMSE'	        - Root mean square error.                            ###
###     'SMAPE'         - Symmetric mean absolute percentage error.          ###
###     'nRMSE_sd'      - Normalized RMSE. RMSE divided by sd of the         ###
### '                     observed.                                          ###
###     'nRMSE_maxmin'  - Normalized RMSE. RMSE divided by range of the      ###
###                       observed.                                          ###
###     'ResShapiroW'   - Shapiro-wilk test of residuals normality:          ###
###                       W statistic.                                       ###
###     'ResShapiroP'   - Shapiro-wilk test of residuals normality:          ###
###                       P value.                                           ###
###                                                                          ###
### 3. "Intervals"  - dataframe with the prediction and confidence intervals ###
###    as implemented in propagate::predcitNLS. number of rows equal         ###
###    NumIntervals.                                                         ###
###     'Dep'           - Name of dependent.                                 ###          
###     'Exp'           - Name of explanatory.                               ###          
###     'ModName'       - Name of the model.                                 ###      
###     'Exp___ModName' - Exp and model name, seperated by "___".            ###
###     'ExpValue'      - Values of the Exp:                                 ###
###        seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals)  ###       
###     'PredMean'      - Mean predicted value for the prediction intervals. ###      
###     'PredLow'       - Lower values for the prediction intervals.         ###
###     'PredHigh'      - Higher values for the prediction intervals.        ###
###     'ConfMean'      - mean predicted value for the confidence intervals. ###       
###     'ConfLow'       - Lower values for the confidence intervals.         ###
###     'ConfHigh'      - Higher values for the confidence intervals.        ###
###                                                                          ###
### 4. "Predicted" - dataframe  with a row for each row of 'DataIn' with the ###
###                  predicted values.                                       ###
###     'Dep'           - Name of dependent.                                 ###          
###     'Exp'           - Name of explanatory.                               ###          
###     'ModName'       - Name of the model.                                 ###      
###     'Exp___ModName' - Exp and model name, separated by "___".            ###
###     'DepValue'      - Values of the dependent variable.                  ### 
###     'ExpValue'      - Values of the explanatory variable.                ### 
###     'Pred'          - Predicted values for the dependent variable.       ###
###                                                                          ###
### 5. "Conv"      - Logical, T if the model converged, F if not.            ###
###                                                                          ###
###----- NOTES: -------------------------------------------------------------###
###                                                                          ###
### 1. Relies on  Levenberg-Marquardt algorithm for model convergence, as    ###
###    implemented in the minpack.lm R package.                              ###
### 2. Confidence and prediction intervals are based on the second-order     ###
###    taylor expansion as implemented in propagate::predcitNLS              ###
### 3. Default starting values for all coefficients is 1, except for         ###
###    the Negative Exponential and sigEmax where the starting values are    ###
###    calculated internally.                                                ###
###                                                                          ###
###----- VERSIONS: ----------------------------------------------------------###
###                                                                          ###
### 1. 09 Apr. 2019 - original code (based on analysis 2 version).           ###
### 2. 10 Apr. 2019 - Add condition for models that converged but have       ###
###    deviance == Inf. They are treated as if they did not converge.        ###
### 3. 21 June 2019 -  add test of residuals normality.                      ###
###                                                                          ###
###----- DETAILS: -----------------------------------------------------------###
###                                                                          ###
### 1. Author:    Yoni Gavish <gavishyoni@gmail.com>                         ###
###                                                                          ###
################################################################################

FitMod = function(DataIn,       
                  ModName,      
                  Dep, 
                  Exp, 
                  NumIntervals = 50, 
                  Alpha        = 0.05){
  ###------------------------------------------------------------------------###
  ### required packages
  
  library(Metrics)
  library(minpack.lm)
  library(hydroGOF)
  library(propagate)
  
  ###------------------------------------------------------------------------###
  ### 1.  Fit the models using non-linear regression:
  ###------------------------------------------------------------------------###
  
  ###------------------------------------------------------------------------###
  ### 1.1.  Linear model.
  if(ModName == "Linear"){
    M1 <- try(nlsLM(Dep ~ a1 + a2 * Exp, 
                    start = list(a1 = 1, 
                                 a2 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    
    Form = "Y = a1 + a2 * X"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.2.  Polynomial second order model.                           
  if(ModName == "Poly2"){
    M1 <- try(nlsLM(Dep ~ a1 + a2 * Exp + a3 * (Exp^2), 
                    start = list(a1 = 1, 
                                 a2 = 1, 
                                 a3 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    Form <- "Y = a1 + a2 * X + a3 * X^2"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.3.  Polynomial third order model.
  if(ModName == "Poly3"){
    M1 <- try(nlsLM(Dep ~ a1 + a2 * Exp + a3 * (Exp^2) + a4 * (Exp^3), 
                    start = list(a1 = 1, 
                                 a2 = 1, 
                                 a3 = 1,
                                 a4 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    Form <- "Y = a1 + a2 * X + a3 * X^2 + a4 * X^3"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.4.  Power model.
  if(ModName == "Power"){
    M1 <- try(nlsLM(Dep ~ a1 * Exp^a2, 
                    start = list(a1 = 1, 
                                 a2 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    Form <- "Y = a1 * X^a2"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.5.  Logarithmic model.
  if(ModName == "Logarithmic"){
    M1 <- try(nlsLM(Dep ~ a1 + a2 * log(Exp), 
                    start = list(a1 = 1, 
                                 a2 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    Form <- "Y = a1 + a2 * log(X)"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.6.  Negative Exponential model.
  if(ModName == "NegativeExponential"){
    L1 <- max(DataIn$Dep)
    L2 <- median(DataIn$Dep)
    L3 <- DataIn[which(DataIn$Dep == L2), "Exp"][1]
    L4 <- 10*(L2-L1)/(L3*L1)
    M1 <- try(nlsLM(Dep ~ a1 * (1 - exp(-a2 * Exp)), 
                    start = list(a1 = L1, 
                                 a2 = L4),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    rm(L1, L2, L3, L4)
    Form <- "Y = a1*[1 - exp(-a2 * X)]"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.7.  Rational model.
  if(ModName == "Rational"){
    M1 <- try(nlsLM(Dep ~ (a1 + a2 * Exp) / (1 + a3 * Exp), 
                    start = list(a1 = 1, 
                                 a2 = 1, 
                                 a3 = 1),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    Form <- "Y = (a1 + a2 * X) / (1 + a3 * x)"
  }
  
  ###------------------------------------------------------------------------###
  ### 1.8.  sigEmax model. 
  if(ModName == "sigEmax"){
    L1 <- min(DataIn$Dep)
    L2 <- max(DataIn$Dep)
    L3 <- 1
    L4 <- median(DataIn$Dep)
    L4 <- 10*(L2-L1)/(L3*L1)
    
    M1 <- try(nlsLM(Dep ~ a1 + a2 * (Exp^a3 / (a4^a3 + Exp^a3)), 
                    start = list(a1 = L1, 
                                 a2 = L2, 
                                 a3 = L3,
                                 a4 = L4),
                    control = nls.lm.control(maxiter = 500),
                    data = DataIn), 
              silent = T)
    rm(L1, L2, L3, L4)
    Form <- "Y = a1 + a2*[X^a3 / (a4^a3 + X^a3)]"
  }
  
  
  
  
  
  ###------------------------------------------------------------------------###
  ### 2.  Check if the model converged.
  ###------------------------------------------------------------------------###
  
  
  Conv <- class(M1) != "try-error"
  
  
  ###------------------------------------------------------------------------###
  ### 3.  Start model info results dataframe.
  ###------------------------------------------------------------------------###
  
  F1 <- data.frame(
    Dep           = Dep,                                   # Name of dependent
    Exp           = Exp,                                   # name of explanatory
    ModName       = ModName,                               # name of the model
    Exp___ModName = paste(Exp, ModName, sep = "___"), # Exp and model name
    Formula       = Form,           # formula of the model
    Converged     = Conv,           # T if the model converged, F if not
    NumCase       = nrow(DataIn),   # number of cases
    a1            = NA,             # First coefficent
    a2            = NA,             # Second coefficent
    a3            = NA,             # Third coefficent
    a4            = NA,             # Fourth coefficent
    Deviance      = NA,             # Deviance
    RMSE	        = NA,             # Root mean square error
    SMAPE         = NA,             # Symmetric mean absolute percentage error
    nRMSE_sd      = NA,	            # normalized RMSE. divided by sd of the observed (range: 0-100)
    nRMSE_maxmin  = NA,             # normalized RMSE. divided by range of the observed (range: 0-100)
    ResShapiroW   = NA,             # Shapiro-wilk test of residuals normality - W statistic
    ResShapiroP   = NA)             # Shapiro-wilk test of residuals normality - P value
  
  
  
  ###------------------------------------------------------------------------###
  ### 4.  If the model convereged:
  ###------------------------------------------------------------------------###
  
  if(Conv){
    if(deviance(M1) != Inf){
      ###----------------------------------------------------------------------###
      ### 4.1. Get prediction and confidence intervals:
      
      
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ###  4.1.1.  Create a sequence along the range of Exp.
      NewExp <- data.frame(Exp = seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals))
      
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ###  4.1.2.  Fit the prediction intervals.
      M2 <- suppressMessages(predictNLS(M1, NewExp, interval = "prediction", alpha = Alpha))
      
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ###  4.1.3.  Fit the condfidence intervals.
      M3 <- suppressMessages(predictNLS(M1, NewExp, interval = "confidence", alpha = Alpha))
      
      ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
      ### 4.1.4.  Arrange prediction and confidence inervals dataframe.
      
      T1       <- M2$summary
      T2       <- M3$summary
      PredConf <- data.frame(Dep           = Dep,
                             Exp           = Exp,
                             ModName       = ModName,
                             Exp___ModName = paste(Exp, ModName, sep = "___"), # Exp and model name
                             ExpValue      = NewExp,
                             PredMean      = T1[, 2],
                             PredLow       = T1[, 5],
                             PredHigh      = T1[, 6],
                             ConfMean      = T2[, 2],
                             ConfLow       = T2[, 5],
                             ConfHigh      = T2[, 6])
      rm(M2, M3, T1, T2)
      names(PredConf)[5] <- "ExpValue"
      
      ###----------------------------------------------------------------------###
      ### 4.2.  Add coefficeint values to Info. 
      
      T1 <- coef(M1)
      F1[1, "a1"] <- as.numeric(T1[1])
      F1[1, "a2"] <- as.numeric(T1[2])
      F1[1, "a3"] <- as.numeric(T1[3])
      F1[1, "a4"] <- as.numeric(T1[4])
      
      
      ###----------------------------------------------------------------------###
      ### 4.3.  Predict for the input data and calculate performance indices. 
      
      P1 <- predict(M1)
      F1[1, "Deviance"]     <- deviance(M1)
      F1[1, "RMSE"]         <- rmse(DataIn$Dep, P1)
      F1[1, "nRMSE_sd"]     <- nrmse(P1, DataIn$Dep, norm = "sd")
      F1[1, "nRMSE_maxmin"] <- nrmse(P1, DataIn$Dep, norm = "maxmin")
      F1[1, "SMAPE"]        <- smape(DataIn$Dep, P1)
      
      
      ###----------------------------------------------------------------------###
      ### 4.4.  Prepare the predicted dataframe.
      Pred <- data.frame(Dep           = Dep,
                         Exp           = Exp,
                         ModName       = ModName, # Name of the model
                         Exp___ModName = paste(Exp, ModName, sep = "___"),
                         DepValue      = DataIn$Dep,	
                         ExpValue      = DataIn$Exp,
                         Pred          = P1)
      
      ###----------------------------------------------------------------------###
      ### 4.5.  Test residuals normality with Shapiro-Wilk test.
      
      A1 <- shapiro.test(residuals(M1))
      F1[1, "ResShapiroW"]  <- A1$statistic
      F1[1, "ResShapiroP"]  <- A1$p.value
      rm(A1)
    }}
  
  
  ###------------------------------------------------------------------------###
  ### 5.  If the model did not convereged- return NA for most values.
  ###------------------------------------------------------------------------###
  
  
  if(!Conv){
    NewExp <- data.frame(Exp = seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals))
    PredConf <- data.frame(Dep           = Dep,
                           Exp           = Exp,
                           ModName       = ModName,
                           Exp___ModName = paste(Exp, ModName, sep = "___"), 
                           ExpValue      = NewExp,
                           PredMean      = NA,
                           PredLow       = NA,
                           PredHigh      = NA,
                           ConfMean      = NA,
                           ConfLow       = NA,
                           ConfHigh      = NA)
    names(PredConf)[5] <- "ExpValue"
    
    Pred <- data.frame(Dep           = Dep,
                       Exp           = Exp,
                       ModName       = ModName, # Name of the model
                       Exp___ModName = paste(Exp, ModName, sep = "___"),
                       DepValue      = DataIn$Dep,	
                       ExpValue      = DataIn$Exp,
                       Pred          = NA)
  } else {if(deviance(M1) == Inf){
    NewExp <- data.frame(Exp = seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals))
    PredConf <- data.frame(Dep           = Dep,
                           Exp           = Exp,
                           ModName       = ModName,
                           Exp___ModName = paste(Exp, ModName, sep = "___"), 
                           ExpValue      = NewExp,
                           PredMean      = NA,
                           PredLow       = NA,
                           PredHigh      = NA,
                           ConfMean      = NA,
                           ConfLow       = NA,
                           ConfHigh      = NA)
    names(PredConf)[5] <- "ExpValue"
    
    Pred <- data.frame(Dep           = Dep,
                       Exp           = Exp,
                       ModName       = ModName, # Name of the model
                       Exp___ModName = paste(Exp, ModName, sep = "___"),
                       DepValue      = DataIn$Dep,	
                       ExpValue      = DataIn$Exp,
                       Pred          = NA)
  }
  } 
  
  
  
  
  
  ###------------------------------------------------------------------------###
  ### 6.  Arrange the return list.
  ###------------------------------------------------------------------------###
  
  
  Return <- list(Model      = M1,
                 Info       = F1,
                 Intervals  = PredConf,
                 Predicted  = Pred, 
                 Conv       = Conv) # T if the model converged, F if not 
  
  
  
} # end FitMod function





################################################################################
###--------------------------------------------------------------------------###
###----------------------                              ----------------------### 
###-----------------                                        -----------------###
###-----------                                                    -----------###
###-----                            END CODE                            -----###
###-----------                                                    -----------###
###-----------------                                        -----------------###
###----------------------                              ----------------------### 
###--------------------------------------------------------------------------###
################################################################################

