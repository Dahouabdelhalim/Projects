################################################################################
###                                                                          ###
###------- NAME: ------------------------------------------------------------###
###                                                                          ###
### FitModAllDepAllExp.R                                                     ###
###                                                                          ###
###------- AIM: -------------------------------------------------------------###
###                                                                          ###
### For each combination of dependent variable, explanatory variable and     ###
### formula, fit an nls model, get predicted values and confidence /         ###
### prediction intervals. Then compare models with AICc based model          ###
### selection.                                                               ### 
###                                                                          ###
###------- INPUT: -----------------------------------------------------------###
###                                                                          ###
### 1.  "Data"     - excel sheet.                                            ###
###     A sheet from the excel file 'VelisEtAl2022_WABIs_Input.xlsx'.        ###
###     Contains the data for the WABIs and explanatory variables for the 40 ###
###     cities. Contains the following columns:                              ###
###  "casename" - A unique case name.                                        ###                    
###  "City"     - Name of the city.                                          ###
###  "CityNum"  - Number of city.                                            ###
###  "Country"  - Name of country.                                           ###
###  "ISO3CountryCode" - 3 letter ISO3 code for countries.                   ###
###  "IncomeCategory"  - Income category of the country. Either "L", "L-M",  ###
###                      "U-M" or "H".                                       ###               
###  "WastePerCapita"  - WABI I-0 - Waste generation rate [Kg*y^-1*p^-1].    ###               
###  "CollectionCoverage"     - WABI I-1.1 - Waste collection coverage.      ###           
###  "WasteQualityCollection" - WABI I-1C - Quality of waste collection.     ###
###                             service.                                     ###        
###  "ControlledDisposal" - WABI I-2 - Controlled treatment and disposal.    ###     
###  "QualityEnvironmentalProtection" - WABI I-2E - Environmental protection ###
###                                   in controlled treatment and disposal.  ###
###  "GDPperCapitaPPPyr_0"            - GDP PPP per capita [USD '000].       ###  
###  "SPI_3_Opportunity"              - Social progress index 3: Opportunity ###
###                                     (SPI-3).                             ###  
###  "CorruptionPerceptionsIndex"     - Corruption perceptions index (CPI).  ###  
###  "GNIperCapita"                   - GNI per capita [USD '000].           ###  
###  "ProportionUrbanCountry"         - Urban population in country - 2010   ###
###                                     census [%].                          ###  
###  "HumanDevelopmentIndex"          - Human development index (HDI).       ###  
###  "PopulationSecondaryEducation"   - Population with at least secondary   ###
###                                     education (%).                       ### 
###  "AdultLiteracyRate" - Adult literacy rate (% of population).            ###
###  "EPI"               - Environmental protection index (EPI).             ###
###                                                                          ###
### 2.  "Meta"     - excel sheet.                                            ###
###     A sheet from the excel file 'VelisEtAl2022_WABIs_Input.xlsx'.        ###
###     Contains metadata for variables in 'Data' in the following columns:  ###
###  "Role"   - The role of the variable, either "Case", "General", "WABIs", ###
###             or "Explanatory".                                            ###
###  "Number" - Column number in 'Data'.                                     ### 
###  "Name"   - Name of the variable.                                        ###
###  "Label"  - Label of h variable (for plotting).                          ###
###  "OrderY" - Order of plotting.                                           ###
###  "OrderX" - Order of plotting.                                           ###
###                                                                          ###
### 3.  "FitMod.R"   -  R function.                                          ###
### Function to fit a non-linear model against one explanatory variable. Can ###
### fit up to 9 different models with minpack.lm::nlsLM. Estimates           ###
### prediction and confidence intervals with propagate::predcitNLS.          ### 
###                                                                          ###
###------- CONTROLS: --------------------------------------------------------###
###                                                                          ###
### 1.  verbose      - Logical,if T, progress is printed to console.         ###
### 2.  MinCompCases - Minimum number of complete cases for dependent.       ###
### 3.  Modnames     - Names of the formulas to fit (see notes below).       ###
### 4.  Alpha        - Alpha values for the confidence/prediction intervals. ###
### 5.  NumIntervals - Number of cases for the confidence/prediction         ###
###                    intervals.                                            ###
###                                                                          ###
###------- PROCESS: ---------------------------------------------------------###
###                                                                          ###
### 1.  Load packages.                                                       ###
### 2.  Define controls.                                                     ###
###    2.1.  General controls.                                               ###
###    2.2.  Models names.                                                   ###
###    2.3.  Confidence/prediction intervals - alpha and number of points.   ###
### 3.  Set directories: Input, function source, output.                     ###
### 4.  Read data and source functions.                                      ###
### 5.  Prepare base dataframes.                                             ###
###    5.1.  Subset the dependent and explanatory list.                      ###
###    5.2.  Start the overall dataframe and counter to monitor convergence. ###
### 6.  Loop over all Dependent variables and for each:                      ###
###    6.1.  Identify dependent.                                             ###
###    6.2.  Create the base data frame for the dependent.                   ###
###    6.3.  If there are not enough cases- next.                            ###
###    6.4.  Start the candidate model list and models details dataframe.    ###
###    6.5.  Loop over all explanatory vars and formulas and for each:       ###
###      6.5.1.  Create the dependent/explanatory dataframe.                 ###
###      6.5.2.  Print progress.                                             ###
###      6.5.3.  Call FitMod.R                                               ###
###      6.5.4.  Store convergence info.                                     ###
###      6.5.5.  If the model converged:                                     ###
###         6.5.5.1.  Update models details dataframe.                       ###
###         6.5.5.2.  Store the nls model in the candidate list.             ###
###         6.5.5.3.  Update the info, predicted and intervals dataframes.   ###
###         6.5.5.4.  Update counter for candidate model list.               ###
###    6.6.  Model selection and ranking - Full (Exp_Mod):                   ###
###      6.6.1.  Model selection table and model rank.                       ###
###      6.6.2.  Add model selection results to Info.                        ###
###      6.6.3.  Add model performance and model selection results to Pred.  ###
###    6.7.  Write dependent results to file.                                ###
###    6.8.  Update main outputs and clean.                                  ###
### 7.  Write to file.                                                       ###
###                                                                          ###
###------- OUTPUT: ----------------------------------------------------------###
###                                                                          ###
###  1.  "WasteAware_Info.csv"      - general info and model selection       ###
###      results for each dependent variable against all explanatory         ###
###      variables and model formula combinations.                           ### 
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
###                       observed.                                          ###
###     'nRMSE_maxmin'  - Normalized RMSE. RMSE divided by range of the      ###
###                       observed.                                          ###
###     'ResShapiroW'   - Shapiro-wilk test of residuals normality:          ###
###                       W statistic.                                       ###
###     'ResShapiroP'   - Shapiro-wilk test of residuals normality:          ###
###                       P value.                                           ###
###     'K'             - Number of parameters in the model (+1 for the      ###
###                       error term).                                       ###
###     'AICc'          - Akike Information Criteria corrected for small     ###
###                       sample size of the model.                          ###
###     'Delta_AICc'    - Delta AICc from the best model. lower values are   ###
###                       better.                                            ###
###     'ModelLik'      - The relative likelihood value of the model.        ###    
###     'AICcWt'        - The Akaike weights.better models have higher       ###
###                       weights. All models sum to 1.                      ###
###     'LL'            - the log-likelihood of each model.                  ###
###     'Rank'          - Model rank from best (lowest deltaAICc) to worse   ###
###                       (highewst DeltaAICc). Lower values are  better     ###
###                       models.                                            ###
###                                                                          ###
###  2.  "WasteAware_Predicted.csv" - The predicted values according to each ###
###      model that converged.                                               ###
###     'Dep'           - Name of dependent.                                 ###
###     'Exp'           - Name of explanatory.                               ###
###     'ModName'       - Name of the model.                                 ###
###     'Exp___ModName' - Exp and model name, separated by "___".            ###
###     'DepValue'      - Values of the dependent variable.                  ### 
###     'ExpValue'      - Values of the explanatory variable.                ### 
###     'Pred'          - Predicted values for the dependent variable.       ###
###     'casename'      - Unique name of the case.                           ###
###     'IncomeCategory' - Income category for plotting later.               ### 
###     'Deviance'      - The deviance.                                      ###
###     'RMSE'	        - Root mean square error.                            ###
###     'SMAPE'         - Symmetric mean absolute percentage error.          ###
###     'nRMSE_sd'      - Normalized RMSE. RMSE divided by sd of the         ###
###                       observed.                                          ###
###     'nRMSE_maxmin'  - Normalized RMSE. RMSE divided by range of the      ###
###                       observed.                                          ###
###     'Rank'          - Model rank when all combinations of explanatory    ###
###                       and formulas are compared in a single AICc table.  ###
###                       Best model as a rank of 1.                         ###
###                                                                          ###
### 3. "WasteAware_Intervals.csv" - dataframe with the prediction and        ###
###    confidence intervals as implemented in propagate::predcitNLS. number  ###
###    of rows for each fitted model equal 'NumIntervals'.                   ###
###     'Dep'           - Name of dependent.                                 ###          
###     'Exp'           - Name of explanatory.                               ###          
###     'ModName'       - Name of the model.                                 ###      
###     'Exp___ModName' - Exp and model name, separated by "___".            ###
###     'ExpValue'      - Values of the Exp:                                 ###
###        seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals)  ###       
###     'PredMean'      - mean predicted value for the prediction intervals. ###      
###     'PredLow'       - Lower values for the prediction intervals.         ###
###     'PredHigh'      - Higher values for the prediction intervals.        ###
###     'ConfMean'      - mean predicted value for the confidence intervals. ###       
###     'ConfLow'       - Lower values for the confidence intervals.         ###
###     'ConfHigh'      - Higher values for the confidence intervals.        ###
###                                                                          ###
###  4.  "WasteAware_converged.csv" - Dataframe with a row for each          ###
###      combination of dependent, explanatory and formula with information  ###
###      about convergence.                                                  ###
###     'Dep'       - Name of dependent.                                     ###          
###     'Exp'       - Name of explanatory.                                   ###          
###     'ModName'   - Name of the model.                                     ###    
###     'Converged' - True if model converged, False if model did not        ###
###                   converge, 'NoCases' if the model was not fitted since  ###
###                   the number of complete cases was lower than            ###
###                   'MinCompCases'.                                        ###
###                                                                          ###
###------- NOTES: -----------------------------------------------------------###
###                                                                          ###
### 1.  Model fitting is done with the function FitMod.R which relies on the ###
###     Levenberg-Marquardt algorithm for model convergence, as implemented  ###
###     in the minpack.lm R package.                                         ###
###                                                                          ###
### 2.  Confidence and prediction intervals are based on the second-order    ###
###     taylor expansion as implemented in propagate::predcitNLS             ###
###                                                                          ###
### 3.  The following formulas can be fitted:                                ###
###     'Linear'               :-: Y = a1 + a2 * X                           ###
###     'Poly2'                :-: Y = a1 + a2 * X + a3 * X^2                ###
###     'Poly3'                :-: Y = a1 + a2 * X + a3 * X^2 + a4 * X^3     ###
###     'Power'                :-: Y = a1 * X^a2                             ###
###     'Logarithmic'          :-: Y = a1 + a2 * log(X)                      ###
###     'NegativeExponential'  :-: Y = a1*[1 - exp(-a2 * X)]                 ###
###     'Rational'             :-: Y = (a1 + a2 * X) / (1 + a3 * x)          ###
###     'sigEmax'              :-: Y = a1 + a2*[X^a3 / (a4^a3 + X^a3)]       ###
###                                                                          ###
### 4.  Output files name have the structure:                                ###
###     4.1.  paste0("WasteAware_Info__"     , Dep, ".csv").                 ###
###     4.2.  paste0("WasteAware_Predicted__", Dep, ".csv").                 ###
###     4.3.  paste0("WasteAware_Intervals__", Dep, ".csv").                 ###
###     4.4.  paste0("WasteAware_converged__", Dep, ".csv").                 ###
###                                                                          ###
### 5.  If a model converged but its deviance == Inf, it is excluded.        ###
###                                                                          ###
###----- VERSIONS: ----------------------------------------------------------###
###                                                                          ###
### 1. 09 Apr. 2019 - original code (based on analysis 2 version).           ###
### 2. 31 May 2019 -                                                         ###
###    2.1.  Added casename and income category to predicted output.         ###
###    2.2.  Fewer dependent, explanatory and functions.                     ###
###    2.3.  Code simplified where possible.                                 ###
### 3. 15 Feb. 2022 -                                                        ###
###    3.1.  Changed input to excel sheets abd update names.                 ###                                
###                                                                          ###
###----- DETAILS: -----------------------------------------------------------###
###                                                                          ###
### 1. Author:    Yoni Gavish <gavishyoni@gmail.com>                         ###
###                                                                          ###
################################################################################



###--------------------------------------------------------------------------###
### 1.  Load packages.
###--------------------------------------------------------------------------###


library(Rcpp)
library(readxl)
library(AICcmodavg)
library(Metrics)
library(minpack.lm)
library(hydroGOF)
library(propagate)


###--------------------------------------------------------------------------###
### 2.  Define controls.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 2.1.  General controls.

verbose      <- T        # if T, progress is printed
MinCompCases <- 15       # minimum number of cases for explanatory or dependent variable


###--------------------------------------------------------------------------###
### 2.2.  Models names.

Modnames <- c("Linear", 
              "Poly2",
              #"Poly3",
              #"Power", 
              "Logarithmic",
              #"NegativeExponential",
              #"Rational",
              "sigEmax")

###--------------------------------------------------------------------------###
### 2.3.  Confidence/prediction intervals - alpha and number of points.

Alpha        <- 0.05
NumIntervals <- 25


###--------------------------------------------------------------------------###
### 3.  Set directories: Input, function source, output.
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis"
Dir3 <- "./Analysis/ModelSelection"


###--------------------------------------------------------------------------###
### 4.  Read data and source functions.
###--------------------------------------------------------------------------###


setwd(Dir1)
Data <- read_excel(path = paste0(Dir1, "/VelisEtAl2022_WABIs_Input.xlsx"),
                   sheet = "Data")
Meta <- read_excel(path = paste0(Dir1, "/VelisEtAl2022_WABIs_Input.xlsx"),
                   sheet = "Meta")

# change to data.frame and factors for back compatibility 
Data <- as.data.frame(Data, stringsAsFactors = T)
Meta <- as.data.frame(Meta, stringsAsFactors = T)

Data$casename        <- factor(Data$casename)
Data$City            <- factor(Data$City)
Data$Country         <- factor(Data$Country)
Data$ISO3CountryCode <- factor(Data$ISO3CountryCode)
Data$IncomeCategory  <- factor(Data$IncomeCategory)

Meta$Role  <- factor(Meta$Role)
Meta$Name  <- factor(Meta$Name)
Meta$Label <- factor(Meta$Label)

setwd(Dir2)
source("FitMod.R")


###--------------------------------------------------------------------------###
### 5.  Prepare base dataframes.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 5.1.  Subset the dependent and explanatory list.

DepList <- Meta[Meta$Role == "WABIs", ]
ExpList <- Meta[Meta$Role == "Explanatory", ]

###--------------------------------------------------------------------------###
### 5.2.  Start the overall dataframe and counter to monitor convergence.

AllFull           <- expand.grid(Modnames, ExpList$Name, DepList$Name)
AllFull           <- AllFull[, c(3, 2, 1)]
names(AllFull)    <- c("Dep", "Exp", "ModName")
AllFull$Converged <- NA
C2                <- 1   #start counter


###--------------------------------------------------------------------------###
### 6.  Loop over all Dependent variables and for each:
###--------------------------------------------------------------------------###


for(i in 1:nrow(DepList)){
# for(i in 1:1){ # check
  
  ###------------------------------------------------------------------------###
  ### 6.1.  Identify dependent.
  
  Dep <- as.character(DepList[i, "Name"])
 
  ###------------------------------------------------------------------------###
  ### 6.2.  Create the base data frame for the dependent.
  
  D0 <- Data[, c(Dep, as.character(ExpList$Name),"casename", "IncomeCategory" )]
  D0 <- D0[complete.cases(D0), ]
  DepList[i, "NumCompCases"] <- nrow(D0)
  
  if(verbose){
    cat(paste0("Dependent ", i, " out of ", nrow(DepList), " --> ", Dep, "   (N=", nrow(D0),  "):\\n"))
  }
  
  ###------------------------------------------------------------------------###
  ### 6.3.  If there are not enough cases- next.
  
  if(nrow(D0) < MinCompCases){
    for(k in 1:nrow(ExpList)){
      for(M in 1:length(Modnames)){
        AllFull[C2, "Converged"] <- "NoCases"
        C2 <- C2 + 1
      }
      rm(M)
    }
    rm(k)
    next()
  }
 
  ###------------------------------------------------------------------------###
  ### 6.4.  Start the candidate model list and models details dataframe.

  Cand.mod <- list()

  ModDep           <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(ModDep) <- c("Dep", "Exp", "ModName", "Exp___ModName")

  CountL <- 1

  ###------------------------------------------------------------------------###
  ### 6.5.  Loop over all explanatory vars and formulas and for each:
  
  
  for(k in 1:nrow(ExpList)){

    Exp <-  as.character(ExpList[k, "Name"])
    if(verbose){
      cat(paste0("    Explanatory ", k, " out of ", nrow(ExpList), " --> ", Exp, "\\n"))
    }
    
    for(M in 1:length(Modnames)){
      
      ###--------------------------------------------------------------------###
      ### 6.5.1.  Create the dependent/explanatory dataframe.
      D1 <- data.frame(Dep = D0[, Dep],
                       Exp = D0[, Exp])
      
      ###--------------------------------------------------------------------###
      ### 6.5.2.  Print progress.
      if(verbose){
        cat(paste0("         Model ", M, " out of ", length(Modnames), " : ", Modnames[M], "    --> Converged? "))
      }
      
      ###--------------------------------------------------------------------###
      ### 6.5.3.  Call FitMod.R
      M1 <- FitMod(DataIn       = D1,
                   ModName      = Modnames[M],
                   Dep          = Dep,
                   Exp          = Exp,
                   NumIntervals = NumIntervals,
                   Alpha        = Alpha)
      
      if(verbose){cat(paste0(M1[["Conv"]], "\\n"))}
      
      ###--------------------------------------------------------------------###
      ### 6.5.4.  Store convergence info.
      AllFull[C2, "Converged"] <- M1[["Conv"]]
      C2 <- C2 + 1
      
      ###--------------------------------------------------------------------###
      ### 6.5.5.  If the model converged:
      
      if(M1[["Conv"]]){
        if(deviance(M1[["Model"]]) == Inf){
          AllFull[C2-1, "Converged"] <- "DevianceInf"
        }
        
        
        if(deviance(M1[["Model"]]) != Inf){
        ###------------------------------------------------------------------###
        ### 6.5.5.1.  Update models details dataframe.
        ModDep[CountL, "Dep"]           <- Dep
        ModDep[CountL, "Exp"]           <- Exp
        ModDep[CountL, "ModName"]       <- Modnames[M] 
        ModDep[CountL, "Exp___ModName"] <- paste(Exp, Modnames[M], sep = "___") 
        
        ###------------------------------------------------------------------###
        ### 6.5.5.2.  Store the nls model in the candidate list.
        Cand.mod[[CountL]] <- M1[["Model"]]
        
        ###------------------------------------------------------------------###
        ### 6.5.5.3.  Update the info, predicted and intervals dataframes.
        
        if(CountL == 1){
          Info      <- M1[["Info"]]
          Pred      <- cbind(M1[["Predicted"]], D0[, c("casename", "IncomeCategory")])
          Intervals <- M1[["Intervals"]]
          }
        if(CountL != 1){
          Info      <- rbind(Info, M1[["Info"]])
          Pred      <- rbind(Pred,  cbind(M1[["Predicted"]], D0[, c("casename", "IncomeCategory")]))
          Intervals <- rbind(Intervals,  M1[["Intervals"]])
          }
        
        ###------------------------------------------------------------------###
        ### 6.5.5.4.  Update counter for candidate model list.
        CountL <- CountL + 1
        
      }} # end converged conditions
      
      rm(M1, D1)
    } # end Modnames loop
    rm(M, Exp)

  } # end ExpList loop
  
  ###------------------------------------------------------------------------###
  ### 6.6.  Model selection and ranking - Full (Exp_Mod):
  
  ###------------------------------------------------------------------------###
  ### 6.6.1.  Model selection table and model rank.
  TabAic      <- aictab(Cand.mod[], ModDep$Exp___ModName, sort = F)
  TabAic$Rank <- rank(TabAic$Delta_AICc)
  
  ###------------------------------------------------------------------------###
  ### 6.6.2.  Add model selection results to Info. 
  T1 <- match(Info$Exp___ModName, TabAic$Modnames)
  Info <- cbind(Info, TabAic[T1, 2:ncol(TabAic)])
  
  ###------------------------------------------------------------------------###
  ### 6.6.3.  Add model performance and model selection results to Pred.
  
  for(j in 1:nrow(Info)){
    Name <- as.character(Info[j, "Exp___ModName"])
    Pred[Pred$Exp___ModName == Name, "Deviance"]     <- Info[j, "Deviance"]
    Pred[Pred$Exp___ModName == Name, "RMSE"]         <- Info[j, "RMSE"]
    Pred[Pred$Exp___ModName == Name, "SMAPE"]        <- Info[j, "SMAPE"]
    Pred[Pred$Exp___ModName == Name, "nRMSE_sd"]     <- Info[j, "nRMSE_sd"]
    Pred[Pred$Exp___ModName == Name, "nRMSE_maxmin"] <- Info[j, "nRMSE_maxmin"]
    Pred[Pred$Exp___ModName == Name, "Rank"]         <- Info[j, "Rank"]
    rm(Name)
  }
  
  rm(T1, j, TabAic)
  
  ###------------------------------------------------------------------------###
  ### 6.7.  Write dependent results to file.
  
  setwd(Dir3)
  write.csv(Info,      paste0("WasteAware_Info__"     , Dep, ".csv"), row.names = F) 
  write.csv(Pred,      paste0("WasteAware_Predicted__", Dep, ".csv"), row.names = F)
  write.csv(Intervals, paste0("WasteAware_Intervals__", Dep, ".csv"), row.names = F)
  write.csv(AllFull,   paste0("WasteAware_converged__", Dep, ".csv"), row.names = F)
  
  ###------------------------------------------------------------------------###
  ### 6.8.  Update main outputs and clean.
  
  # if(i == 1){
  #   OutInfo      <- Info
  #   OutPred      <- Pred
  #   OutIntervals <- Intervals
  # }
  # if(i != 1){
  #   OutInfo      <- rbind(OutInfo,      Info)
  #   OutPred      <- rbind(OutPred,      Pred)
  #   OutIntervals <- rbind(OutIntervals, Intervals)
  # }
  
  # clean
  rm(Dep, D0, Cand.mod, ModDep, CountL, Info, Pred, Intervals)
} # end DepList loop 

rm(i)


###--------------------------------------------------------------------------###
### 7.  Write to file.
###--------------------------------------------------------------------------###


#setwd(Dir3)
#write.csv(OutInfo,      "WasteAware_Info.csv",      row.names = F)
#write.csv(OutPred,      "WasteAware_Predicted.csv", row.names = F)
#write.csv(OutIntervals, "WasteAware_Intervals.csv", row.names = F)
#write.csv(AllFull,      "WasteAware_converged.csv", row.names = F)


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




