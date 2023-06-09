################################################################################
###                                                                          ###
###---- NAME:                                                                ###
###                                                                          ###
### Fit_CondRandomForest_VarImp.R                                            ###
###                                                                          ###
###---- AIM:                                                                 ###
###                                                                          ###
### Fit conditional randomforest for each WABI against the 9 explanatory     ###
### variables. Extract the regular and conditional variable importance as    ###
### well as the RMSE. Add the RMSE of the best non-linear model for          ###
### comparison.                                                              ###
###                                                                          ###
###---- INPUT:                                                               ###
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
###  3.  "WasteAware_Info.csv"      - general info and model selection       ###
###      results for each dependent variable against all explanatory         ###
###      variables and model formula combinations (non linear models).       ### 
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
###---- PROCESS:                                                             ###
###                                                                          ###
### 1.  Packages.                                                            ###
### 2.  Directories.                                                         ###
### 3.  Controls.                                                            ###
### 4.  Read data and metadata  + arrange.                                   ###
### 5.  Loop over all dependent, and for each:                               ###
###    5.1.  Identify dependent and explanatory.                             ###
###    5.2.  Subset data and keep complete cases.                            ###
###    5.3.  Create the function for conditional randomforest.               ###
###    5.4.  Fit conditional randomforest.                                   ###
###    5.5.  Extract variable importance.                                    ###
###    5.6.  Order by names and rank decreasing.                             ###
###    5.7.  Collect in one dataframe.                                       ### 
###    5.8.  Update output.                                                  ###
###    5.9.  Add the RMSE of the best non-linear model for comparison.       ###
### 6.  Write output to file.                                                ###
###                                                                          ###
###---- OUTPUT:                                                              ###
###                                                                          ###
### 1.    VarImp_CRF.csv                                                     ###
###     Dataframe with the regular and conditional variable importance for   ###
###     each dependent variable.                                             ###
###   "Dep"         - Name of the dependent variable.                        ###
###   "Exp"         - Name of the explanatory variable.                      ###
###   "VarImpF"     - Unconditional variable importance.                     ### 
###   "VarImpT"     - conditional variable importance.                       ###
###   "RankVarImpF" - Rank of unconditional var.imp. (low = better).         ###
###   "RankVarImpT" - Rank of conditional var.imp. (low = better).           ###
###                                                                          ###
### 2.    DepList_CRF_rmse.csv                                               ###
###     Dataframe with the RMSE of the conditional random forest and the     ###
###     best non-linear model. Contains the following columns:               ###
###  "Role"   - The role of the variable, either "Case", "General", "WABIs", ###
###             or "Explanatory".                                            ###
###  "Number" - Column number in 'Data'.                                     ### 
###  "Name"   - Name of the variable.                                        ###
###  "Label"  - Label of h variable (for plotting).                          ###
###  "OrderY" - Order of plotting.                                           ###
###  "OrderX" - Order of plotting.                                           ###
###  "rmse"	  - RMSE of the conditional random forest.                       ###
###  "rmse_Best_NonLinear" - RMSE of the best non-linear model selection.    ###
###                                                                          ###
###---- NOTES:                                                               ###
###                                                                          ###
### 1.  MinSplit is set as floor(N/4), with N being the number of cases.     ###
### 2.  Minbucket is set as floor(N/6), with N being the number of cases.    ###
###                                                                          ###
###---- VERSIONS:                                                            ###
###                                                                          ###
### 1.  14 june 2019 - Original code.                                        ###
###                                                                          ###
###---- DETAILS:                                                             ###
###                                                                          ###
### 1.  Author:    Yoni Gavish <gavishyoni.gmail.com>                        ###
###                                                                          ###
################################################################################



###--------------------------------------------------------------------------###
### 1.  Packages.
###--------------------------------------------------------------------------###


library(readxl)
library(randomForest)
library(varImp)
library(party)
library(Metrics)
#library(minpack.lm)
#library(hydroGOF)


###--------------------------------------------------------------------------###
### 2.  Directories.
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/RandomForest"
Dir3 <- "./Analysis/ModelSelection"


###--------------------------------------------------------------------------###
### 3.  Controls.
###--------------------------------------------------------------------------###


verbose   <- T   # if T, progress is printed
NTree     <- 1000 # number of trees in randomforest models


###--------------------------------------------------------------------------###
### 4.  Read data and metadata  + arrange.
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


# subset the dependent and explanatory list
DepList <- Meta[Meta$Role == "WABIs", ]
ExpList <- Meta[Meta$Role == "Explanatory", ]


###--------------------------------------------------------------------------###
### 5.  Loop over all dependent, and for each:
###--------------------------------------------------------------------------###


for (k in 1:nrow(DepList)) {
  
  ###------------------------------------------------------------------------###
  ### 5.1.  Identify dependent and explanatory.
  
  Dep <- as.character(DepList[k, "Name"])
  Exp <- as.character(ExpList[, "Name"])
  if(verbose){
    cat(paste0("Dependent ", k, " out of ", nrow(DepList), " --> ", Dep,  "  :\\n"))
  }
  
  ###------------------------------------------------------------------------###
  ### 5.2.  Subset data and keep complete cases.
  D1 <- Data[,c(Dep, Exp)] 
  D1 <- D1[complete.cases(D1), ]
  
  
  ###------------------------------------------------------------------------###
  ### 5.3.  Create the function for conditional randomforest.
  Fun <- as.formula(paste0(Dep, " ~ ", paste(Exp, collapse = " + ")))
  
  ###------------------------------------------------------------------------###
  ### 5.4.  Fit conditional randomforest.

  C1 <- cforest(Fun, data = D1, 
                control = cforest_unbiased(ntree     = NTree, 
                                           minsplit  = floor(nrow(D1)/4), 
                                           minbucket = floor(nrow(D1)/6))) 
                                           #minsplit  = MinSplit, 
                                           #minbucket = MinBucket))
  
  DepList[k, "rmse"] <- rmse(D1[, Dep], as.numeric(predict(C1)))
  ###------------------------------------------------------------------------###
  ### 5.5.  Extract variable importance.
  
  C2 <- varimp(C1, conditional = F)
  C3 <- varimp(C1, conditional = T)
  
  ###------------------------------------------------------------------------###
  ### 5.6.  Order by names and rank decreasing.
  C2 <- C2[order(names(C2))]
  C3 <- C3[order(names(C3))]
  H2 <- rank(-C2)
  H3 <- rank(-C3)
  
  ###------------------------------------------------------------------------###
  ### 5.7.  Collect in one dataframe. 
  
  C6 <- data.frame(Dep            = Dep, 
                   Exp            = names(C2), 
                   VarImpF        = C2, 
                   VarImpT        = C3, 
                   RankVarImpF    = H2, 
                   RankVarImpT    = H3)
  
  # remove row names
  rownames(C6) <- c()
  
  ###------------------------------------------------------------------------###
  ### 5.8.  Update output.
  if(k == 1){
    VarImp <- C6}
  if(k != 1){
    VarImp <- rbind(VarImp, C6)
  }
  
  ###------------------------------------------------------------------------###
  ### 5.9.  Add the RMSE of the best non-linear model for comparison.
  
  setwd(Dir3)
  G1 <- read.csv(paste0("WasteAware_Info__", Dep, ".csv"), stringsAsFactors = T)
  G1 <- G1[!is.na(G1$Rank_1PerExp),]
  DepList[k, "rmse_Best_NonLinear"] <- G1[G1$Rank_1PerExp == 1, "RMSE"]
  
  rm(Dep, Exp, D1, Fun, C1, C2, C3,  C6, H2, H3, G1)
}

###--------------------------------------------------------------------------###
### 6.  Write output to file.
###--------------------------------------------------------------------------###


setwd(Dir2)
write.csv(VarImp,  "VarImp_CRF.csv",       row.names = F)
write.csv(DepList, "DepList_CRF_rmse.csv", row.names = F)




################################################################################
###--------------------------------------------------------------------------###
###   --------------------------------------------------------------------   ###
###      --------------------------------------------------------------      ###
###                                                                          ###
###                                    END                                   ###
###                                                                          ###
###      --------------------------------------------------------------      ###
###   --------------------------------------------------------------------   ###
###--------------------------------------------------------------------------###
################################################################################




