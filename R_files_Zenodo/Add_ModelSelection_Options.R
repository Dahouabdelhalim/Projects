################################################################################
###                                                                          ###
###----- NAME: --------------------------------------------------------------###
###                                                                          ###
### Add_ModelSelection_Options.R                                             ###
###                                                                          ###
###----- AIM: ---------------------------------------------------------------###
###                                                                          ###
### Run the two steps model selection analysis for each WABI. In the first   ###
### step, the best formula for each explanatory variable is selected, In the ###
### second step the 9 best models, one for each explanatory variable, are    ###
### compared.                                                                ###
###                                                                          ###
###----- INPUT --------------------------------------------------------------###
###                                                                          ###
###  1.  "WasteAware_Info.csv"                                               ###
###      General info and model selection results for the dependent variable ###
###      against all explanatory variables and model formula combinations.   ### 
###   'Dep'           - Name of dependent.                                   ###
###   'Exp'           - Name of explanatory.                                 ###
###   'ModName'       - Name of the model.                                   ###
###   'Exp___ModName' - Exp and model name, separated by "___".              ###
###   'Formula'       - Formula of the model.                                ###
###   'Converged'     - Logical, T if the model converged, F if not.         ###
###   'NumCase'       - Number of cases.                                     ###
###   'a1'            - First coefficient, NA otherwise.                     ###
###   'a2'            - Second coefficient, NA otherwise.                    ###
###   'a3'            - Third coefficient, NA otherwise.                     ###
###   'a4'            - Fourth coefficient, NA otherwise.                    ###
###   'Deviance'      - The deviance.                                        ###
###   'RMSE'	        - Root mean square error.                              ###
###   'SMAPE'         - Symmetric mean absolute percentage error.            ###
###   'nRMSE_sd'      - Normalized RMSE. RMSE divided by sd of the observed. ###
###   'nRMSE_maxmin'  - Normalized RMSE. RMSE divided by range of the        ###
###                     observed.                                            ###
###   'ResShapiroW'   - Shapiro-wilk test of residuals normality:            ###
###                     W statistic.                                         ###
###   'ResShapiroP'   - Shapiro-wilk test of residuals normality:            ###
###                     P value.                                             ###
###   'K'             - Number of parameters in the model (+1 for the error  ###
###                     term).                                               ###
###   'AICc'          - Akike Information Criteria corrected for small       ###
###                     sample size of the model.                            ###
###   'Delta_AICc'    - Delta AICc from the best model. lower values are     ###
###                     better.                                              ###
###   'ModelLik'      - The relative likelihood value of the model.          ###    
###   'AICcWt'        - The Akaike weights. Better models have higher        ###
###                     weights. All models sum to 1.                        ###
###   'LL'            - the log-likelihood of each model.                    ###
###   'Rank'          - Model rank from best (lowest deltaAICc) to worse     ###
###                     (highest DeltaAICc). Lower values are  better        ###
###                     models.                                              ###
###                                                                          ###
###----- PROCESS: -----------------------------------------------------------###
###                                                                          ###
### 1.  Packages.                                                            ###
### 2.  Directories.                                                         ### 
### 3.  List all 'Info' files and for each:                                  ###   
###    3.1.  Read the Info dataframe.                                        ### 
###    3.2.  Rank models (best =1, worst =4) within each explanatory         ###
###          variable.                                                       ###
###    3.3.  Add the rank and AICc weight when taking only the best model    ### 
###          for each exp.                                                   ###
###    3.4.  Write to file and clean workspace.                              ### 
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### 1.  Dataframe with the same name and columns as the input file with the  ###
###     following columns added:                                             ###
###   "Rank_WithinExp"   - Model rank from best (1) to worst (4) when        ###
###                        considering only the 4 models of each explanatory ###
###                        variables.                                        ###
###   "AICcWt_WithinExp" - AICc weight when considering only the 4 models of ###
###                        each explanatory variables.                       ###
###   "Rank_1PerExp"     - Model rank from best (1) to worst (11) when       ###
###                        taking only the the 1 best model of each          ###
###                        explanatory variable.                             ###
###   "AICcWt_1PerExp"   - AICc weight when taking only the the 1 best model ###
###                        of each explanatory variable.                     ###
###                                                                          ###
###----- NOTES --------------------------------------------------------------###
###                                                                          ###
### 1.                                                                       ###
###                                                                          ###
###----- VERSIONS -----------------------------------------------------------###
###                                                                          ###
### 1.  4 Jun. 2019 - Original code.                                         ###
###                                                                          ###
###----- DETAILS ------------------------------------------------------------###
###                                                                          ###
### 1. Author:      Gavish Yoni <gavishyoni@gmail.com>                       ###
###                                                                          ###
################################################################################



###--------------------------------------------------------------------------###
### 1.  Packages.                                         
###--------------------------------------------------------------------------###


library(qpcR)


###--------------------------------------------------------------------------###
### 2.  Directories.                                         
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/ModelSelection"


###--------------------------------------------------------------------------###
### 3.  List all 'Info' files and for each:                                         
###--------------------------------------------------------------------------###


setwd(Dir1)
T1 <- list.files(Dir1, pattern = "WasteAware_Info__")
#T1 <- list.files(Dir1, pattern = "WasteAware_Info__WastePerCapita")
for(i in T1){
  
  ###------------------------------------------------------------------------###
  ### 3.1.  Read the Info dataframe.
  D1 <- read.csv(i)
  
  ###------------------------------------------------------------------------###
  ### 3.2.  Rank models (best =1, worst =4) within each explanatory variable.
  T2 <- unique(D1$Exp)
  for(j in 1:length(T2)){
    D2                  <- D1[D1$Exp == T2[j], ]
    D2$Rank_WithinExp   <- rank(D2$AICc)
    D2$AICcWt_WithinExp <- akaike.weights(D2$AICc)$weights
    if(j == 1){D3 <- D2} else {D3 <- rbind(D3, D2)}
    rm(D2)  
  }# end j loop
  
  ###------------------------------------------------------------------------###
  ### 3.3.  Add the rank and AICc weight when taking only the best model for each exp.   
  
  D3$Rank_1PerExp   <- NA
  D3$AICcWt_1PerExp <- NA
  
  D4 <- D3[D3$Rank_WithinExp != 1, ]
  D5 <- D3[D3$Rank_WithinExp == 1, ]
  
  D5$Rank_1PerExp   <- rank(D5$AICc)
  D5$AICcWt_1PerExp <- akaike.weights(D5$AICc)$weights
  
  D6 <- rbind(D5, D4)
  
  ###------------------------------------------------------------------------###
  ### 3.4.  Write to file and clean workspace. 
  write.csv(D6, i, row.names = F)
  rm(D1, T2, j, D3, D4, D5, D6)
} # end i loop



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


