################################################################################
###                                                                          ###
###----- NAME: --------------------------------------------------------------###
###                                                                          ###
### Plot_Fig5_WasteTreatment.R                                               ###
###                                                                          ###
###----- AIM: ---------------------------------------------------------------###
###                                                                          ###
### Plot the best model for the waste Treatment coverage and quality         ###
### variables, a boxplot based on income category and the trade-off between  ###
### coverage and quality.                                                    ###                                              ###
###                                                                          ###
###----- INPUT --------------------------------------------------------------###
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
###  2.  "WasteAware_Info.csv"                                               ###
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
###   'K'             - Number of parameters in the model (+1 for the error  ###
###                     term).                                               ###
###   'AICc'          - Akike Information Criteria corrected for small       ###
###                     sample size of the model.                            ###
###   'Delta_AICc'    - Delta AICc from the best model. lower values are     ###
###                     better.                                              ###
###   'ModelLik'      - The relative likelihood value of the model.          ###    
###   'AICcWt'        - The Akaike weights.better models have higher         ###
###                     weights. All models sum to 1.                        ###
###   'LL'            - the log-likelihood of each model.                    ###
###   'Rank'          - Model rank from best (lowest deltaAICc) to worse     ###
###                     (highest DeltaAICc). Lower values are  better        ###
###                     models.                                              ###
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
###  3.  "WasteAware_Predicted.csv"                                          ###
###      The predcited values according to each model that converged.        ###
###   'Dep'            - Name of dependent.                                  ###
###   'Exp'            - Name of explanatory.                                ###
###   'ModName'        - Name of the model.                                  ###
###   'Exp___ModName'  - Exp and model name, separated by "___".             ###
###   'DepValue'       - Values of the dependent variable.                   ### 
###   'ExpValue'       - Values of the explanatory variable.                 ### 
###   'Pred'           - Predicted values for the dependent variable.        ###
###   'casename'       - Unique name for the case.                           ###
###   'IncomeCategory' - Income category for the case.                       ###
###   'Deviance'       - The deviance.                                       ###
###   'RMSE'	         - Root mean square error.                             ###
###   'SMAPE'          - Symmetric mean absolute percentage error.           ###
###   'nRMSE_sd'       - Normalized RMSE. RMSE divided by sd of the          ###
###                      observed.                                           ###
###   'nRMSE_maxmin'   - Normalized RMSE. RMSE divided by range of the       ###
###                      observed.                                           ###
###   'Rank'           - Model rank when all combinations of explanatory and ###
###                      formulas are compared in a single AICc table. Best  ###
###                      model as a rank of 1.                               ###
###                                                                          ###
### 4.  "WasteAware_Intervals.csv"                                           ###
###     Dataframe with the prediction and  confidence intervals as           ###
###     implemented in propagate::predcitNLS. number of rows for each        ###
###     fitted model equal 'NumIntervals'.                                   ###
###   'Dep'           - Name of dependent.                                   ###          
###   'Exp'           - Name of explanatory.                                 ###          
###   'ModName'       - Name of the model.                                   ###      
###   'Exp___ModName' - Exp and model name, seperated by "___".              ###
###   'ExpValue'      - Values of the Exp:                                   ###
###      seq(min(DataIn$Exp), max(DataIn$Exp), length.out = NumIntervals)    ###       
###   'PredMean'      - mean predicted value for the prediction intervals.   ###      
###   'PredLow'       - Lower values for the prediction intervals.           ###
###   'PredHigh'      - Higher values for the prediction intervals.          ###
###   'ConfMean'      - mean predicted value for the confidence intervals.   ###       
###   'ConfLow'       - Lower values for the confidence intervals.           ###
###   'ConfHigh'      - Higher values for the confidence intervals.          ###
###                                                                          ###
###----- CONTROLS: ----------------------------------------------------------###
###                                                                          ###
### 1.  DoLog                                                                ###
###     If T, x axis is printed on a log scale. if F, x axis is printed on   ###
###     a regular arithmetic scale.                                          ###
###                                                                          ###
###----- PROCESS: -----------------------------------------------------------###
###                                                                          ###
### 1.  Packages.                                                            ###                               
### 2.  Directories.                                                         ###                                    
### 3.  Controls.                                                            ###                                   
###    3.1.  General controls.                                               ###
###    3.2.  Axis / Strip labels.                                            ###
###    3.3.  Background shading (traffic light) - coverage.                  ###
###    3.4.  Background shading (traffic light) - quality.                   ###
###    3.5.  Background shading (traffic light) - trade-off.                 ###
###    3.6.  Color scheme - Income categories.                               ### 
###    3.7.  Color scheme - general.                                         ### 
###    3.8.  Controls - to console.                                          ###
###    3.9.  Controls - to file.                                             ###
### 4.  Read data.                                                           ###                                     
### 5.  Dep1 - prepare data for best model panel:                            ###
###    5.1.  Identify the top model and set the axis and strip labels.       ###
###    5.2   Subset plotting data for observed and predicted. Set fill and   ###
###          border colors.                                                  ###
###    5.3.  Subset plotting data for confidence and prediction intervals.   ### 
###    5.4.  Get range of Y values and set Y axis limits.                    ###
###    5.5.  Get y axis limits in the plot.                                  ###
###    5.6.  Arrange title for facet.                                        ###
###    5.7.  Set Y range manually if needed and adjust prediction/confidence ###
###          polygons.                                                       ###
### 6.  Dep2 - prepare data for best model panel:                            ###
###    6.1.  Identify the top model and set the axis and strip labels.       ###
###    6.2   Subset plotting data for observed and predicted. Set fill and   ###
###          border colors.                                                  ###
###    6.3.  Subset plotting data for confidence and prediction intervals.   ### 
###    6.4.  Get range of Y values and set Y axis limits.                    ###
###    6.5.  Get y axis limits in the plot.                                  ###
###    6.6.  Arrange title for facet.                                        ###
###    6.7.  Set Y range manually if needed and adjust prediction/confidence ###
###          polygons.                                                       ###
### 7.  Trade-off - prepare data for the plot:                               ###
###    7.1.  Subset data from raw WasteAware data.                           ###  
###    7.2.  Factor and set levels and labels of the income category.        ###
###    7.3.  Add label for strip and points fill/border colors.              ###
###    7.4.  Set ranges.                                                     ###
### 8.  Prepare data for the histogram plots:                                ###
###    8.1.  Subset data from raw WasteAware data.                           ###
###    8.2.  Duplicate for the 'Overall' category and bind.                  ###
###    8.3.  Add column with the dependent variable and rename the Dependent ###
###          column.                                                         ###
###    8.4.  Add facet info.                                                 ###
###    8.5.  Bind.                                                           ###
###    8.6.  Factor and set levels and labels of the income category.        ###
###    8.7.  Add points fill/border colors.                                  ###
###    8.8.  Get y axis limits in the plot.                                  ###
###    8.9.  Create a color and border vector for the scales.                ###
###    8.10. Create the dataframes for the shades.                           ###
### 9.   Dep1 - Plot best model panel:                                       ###
### 10.  Dep2 - Plot best model panel:                                       ###
### 11.  Trade-off - Plot panel:                                             ###
### 12.  Plot boxplot panel:                                                 ###
### 13.  Arrange the panels                                                  ###
### 14.  Write to console or file.                                           ###                           
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### 1.  A figure with 4 panels:                                              ### 
###   1.1.  Left top panel:     best model, Dep1.                            ### 
###   1.2.  Left bottom panel:  best model, Dep2.                            ###
###   1.3.  Right top panel:    Dep2 vs. Dep1.                               ### 
###   1.4.  Right bottom panel: Boxplots (Dep1 and Dep2) .                   ###
###                                                                          ###
###----- NOTES --------------------------------------------------------------###
###                                                                          ###
### 1.   AICc weights is based on the 1 model per explanatory analysis, in   ###
###      which first a single model is selected for each explanatory         ###
###      variable from the 4 options (4 functions) and then model selection  ###
###      is done one the 9 explanatory variables, each represented by a      ###
###      single model.                                                       ###
###                                                                          ###
###----- VERSIONS -----------------------------------------------------------###
###                                                                          ###
### 1.  3 Jun. 2019 - Original code.                                         ###
### 2.  21 June 2019:                                                        ###
###    2.1.  Updated ranges of quality traffic lights.                       ###
###    2.2.  Add option to set ranges of Y axis manually, and ensure the     ###
###          polygons are still plotted.                                     ###
### 3.  16 Feb. 2022 -                                                       ###
###    3.1.  Changed input to the excel table.                               ###
###    3.2.  Set stringsAsFactors == T when reading data.                    ###
###                                                                          ###
###----- DETAILS ------------------------------------------------------------###
###                                                                          ###
### 1. Author:      Gavish Yoni <gavishyoni@gmail.com>                       ###
###                                                                          ###
################################################################################



###--------------------------------------------------------------------------###
### 1.  Packages.                                         
###--------------------------------------------------------------------------###


library(readxl)
library(ggplot2)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(extrafont)           # for pdf output


###--------------------------------------------------------------------------###
### 2.  Directories.                                         
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/ModelSelection"
Dir3 <- "./Analysis/Figures"


###--------------------------------------------------------------------------###
### 3.  Controls.                                         
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 3.1.  General controls.

ToFile   <- T
FormOut <- c(".png", ".pdf")     # formats of output files
#FormOut  <- ".png"                # formats of output files

DoLog <- F
if(DoLog){
  Name     <- "VelisEtAl_Fig5_WasteTreatment_LogScale"
}
if(!DoLog){
  Name     <- "VelisEtAl_Fig5_WasteTreatment"
}

Dep1 <- "ControlledDisposal" 
Dep2 <- "QualityEnvironmentalProtection"

###--------------------------------------------------------------------------###
### 3.2.  Axis / Strip labels.

if(DoLog){
  LabX1  <- "USD '000 (log scale)"
}
if(!DoLog){
  LabX1  <- "USD '000" 
}

LabY1  <- "I-2: Controlled recovery and disposal"
Label1 <- "A) Controlled: GNI per capita" 

LabX2  <- "SPI-3"    
LabY2  <- "I-2E: Environmental protection in I-2"
Label2 <- "B) Quality: Social progress index 3: Opportunity"    

LabX3  <- "I-2: Controlled recovery and disposal"
LabY3  <- "I-2E: Environmental protection in I-2"
Label3 <- "C)  Quality vs. Controlled"

LabX4  <- "Income category"
LabY4  <- "Percent"
Label4 <- c("D) Controlled", 
            "E) Quality")

RangeY1 <- c(-3, 103)
RangeY2 <- c(-3, 103)

###--------------------------------------------------------------------------###
### 3.3.  Background shading (traffic light) - coverage.

Shade1 <- data.frame(Name  = c("B1", "B2", "B3", "B4", "B5"),
                     Ymin  = c(RangeY1[1], 
                               49.5, 
                               74.5, 
                               84.5, 
                               94.5), 
                     Ymax  = c(49.5, 
                               74.5, 
                               84.5,  
                               94.5, 
                               RangeY1[2]), 
                     Color = c('#FFCFC0',
                               '#FFEDBF', 
                               '#FDFFBF',
                               '#D4FFC9',
                               '#BBFFCF'), 
                     Alpha = c(1, 1, 1, 1, 1),
                     Facet = rep(Label4[1], 5)) 

###--------------------------------------------------------------------------###
### 3.4.  Background shading (traffic light) - quality.

Shade2 <- data.frame(Name  = c("B1", "B2", "B3", "B4", "B5"),
                     Ymin  = c(RangeY2[1], 
                               20, 
                               40, 
                               60, 
                               80), 
                     Ymax  = c(20, 
                               40,  
                               60, 
                               80, 
                               RangeY2[2]), 
                     Color = c('#FFCFC0',
                               '#FFEDBF', 
                               '#FDFFBF',
                               '#D4FFC9',
                               '#BBFFCF'), 
                     Alpha = c(1, 1, 1, 1, 1),
                     Facet = rep(Label4[2], 5))

###--------------------------------------------------------------------------###
### 3.5.  Background shading (traffic light) - trade-off.

Shade3 <- data.frame(Name  = c("B1", "B2",  "B3",  "B4", "B5"),
                     Xmin  = rep(RangeY1[1], 5), 
                     Ymin  = rep(RangeY2[1], 5),
                     Xmax  = c(49.5, 
                               74.5, 
                               84.5,  
                               94.5, 
                               RangeY1[2]), 
                     Ymax  = c(20,  ## Quality - y axis
                               40, 
                               60,  
                               80, 
                               RangeY2[2]), 
                     Color = c('#FFCFC0',
                               '#FFEDBF',
                               '#FDFFBF',
                               '#D4FFC9',
                               '#BBFFCF'), 
                     Alpha = c(1, 1, 1,  1, 1))


###--------------------------------------------------------------------------###
### 3.6.  Color scheme - Income categories. 

### A Quiet Man
### https://www.colourlovers.com/palette/4649674/A_Quiet_Man

IncCat <-data.frame(IncomeCategory = c("L",   
                                       "L-M", 
                                       "U-M", 
                                       "H",
                                       "Overall"),
                    Color = c("#F2E8F3",
                              "#BCA4B2",
                              "#AE8588",
                              "#71779B", 
                              "gray90"),
                    Border = c("black", 
                               "black",
                               "black",
                               "black",
                               "black"),
                    Label = c("L",   
                              "L-M", 
                              "U-M", 
                              "H",
                              "Overall"))

### Blues https://www.colourlovers.com/palette/4649704/Tuesdays
# "#B5FAE1",
# "#88EDEF",
# "#22C6E2",
# "#070B4F",

###--------------------------------------------------------------------------###
### 3.7.  Color scheme - general. 

ColText           <- "black"
ColPred           <- "black"
ColConfIntervals  <- "gray50"
FillConfIntervals <- "gray50"
ColPredIntervals  <- "gray80"
FillPredIntervals <- "gray80"

###--------------------------------------------------------------------------###
### 3.8.  Controls - to console.

if(!ToFile){
  
  ##--------------------------------------------------------------------------##
  ## Panel width and height
  PanWidth  <- c(11, 10)
  PanHeight <- c(10, 10)
  
  ##--------------------------------------------------------------------------##
  ## Grid
  SizeGrid <- 0.3
  ColGrid  <- "gray75"
  TypeGrid <- "dotted"
  
  ##--------------------------------------------------------------------------##
  ## Model text - size, location
  SizeTextMod <- 1
  LocTextMod  <- c(0.75, 0.01)
  
  ##--------------------------------------------------------------------------##
  ## Observed values Points - size, color, shape
  SizePoint1  <- 3
  SizePoint2  <- 3
  SizePoint3  <- 3
  SizePoint4  <- 1.5
  Alpha4     <- 0.5
  ColPoints  <- "black"
  ShapePoint <- 21
  
  ##--------------------------------------------------------------------------##
  ## predcited values lines - size, 
  SizeLines  <- 1
  
  ##--------------------------------------------------------------------------##
  ## Best model text
  HJustText <- 0
  VJustText <- 1
  SizeText  <- 5
  TextXProp <- 0.6 
  TextYProp <- 0.3
  
  
  ##--------------------------------------------------------------------------##
  ## Ribbons:
  Alpha1 <- 0.5
  Alpha2 <- 0.6
  
  ##--------------------------------------------------------------------------##
  ## Axes text and title
  SizeTextX     <- 10
  SizeTextY     <- 10
  SizeTitX      <- 11
  SizeTitY      <- 11
  SizeStripText <- 13
  VjustX        <- 0.2
  
  ##--------------------------------------------------------------------------##
  ## Histogram - vertical line, color and type
  ColVLine  <- "gray40"
  TypeVLine <- "dotted"
}

###--------------------------------------------------------------------------###
### 3.9.  Controls - to file.

if(ToFile){
  ###------------------------------------------------------------------------###
  ### figure dimension/resolution
  Dim <- c(31, 30) # width/height cm
  
  ##--------------------------------------------------------------------------##
  ## Panel width and height
  PanWidth  <- c(11, 10)
  PanHeight <- c(10, 10)
  
  ##--------------------------------------------------------------------------##
  ## Grid
  SizeGrid <- 0.3
  ColGrid  <- "gray75"
  TypeGrid <- "dotted"
  
  ##--------------------------------------------------------------------------##
  ## Model text - size, location
  SizeTextMod <- 1
  LocTextMod  <- c(0.75, 0.01)
  
  ##--------------------------------------------------------------------------##
  ## Observed values Points - size, color, shape
  SizePoint1 <- 3
  SizePoint2 <- 3
  SizePoint3 <- 3
  SizePoint4 <- 3
  Alpha4     <- 0.5
  ColPoints  <- "black"
  ShapePoint <- 21
  
  ##--------------------------------------------------------------------------##
  ## predicated values lines - size, 
  SizeLines  <- 1
  
  ##--------------------------------------------------------------------------##
  ## Best model text
  HJustText <- 0
  VJustText <- 1
  SizeText  <- 3.0
  TextXProp <- 0.6 
  TextYProp <- 0.3
  
  ##--------------------------------------------------------------------------##
  ## Ribbons:
  Alpha1 <- 0.5
  Alpha2 <- 0.6
  
  ##--------------------------------------------------------------------------##
  ## Axes text and title
  SizeTextX     <- 12
  SizeTextY     <- 12
  SizeTitX      <- 12
  SizeTitY      <- 12
  SizeStripText <- 14
  VjustX        <- 0.2
  
  ##--------------------------------------------------------------------------##
  ## Histogram - vertical line, color and type
  ColVLine  <- "gray40"
  TypeVLine <- "dotted"
}


###--------------------------------------------------------------------------###
### 4.  Read data.                                         
###--------------------------------------------------------------------------###


# read the WABI input data
setwd(Dir1)
Data <- read_excel(path = paste0(Dir1, "/VelisEtAl2022_WABIs_Input.xlsx"),
                   sheet = "Data")

# change to data.frame and factors for back compatibility
Data <- as.data.frame(Data, stringsAsFactors = T)

for(i in 1:ncol(Data)){
  if(class(Data[, i]) == "character"){
    Data[, i] <- factor(Data[, i])
  }
}

# read the non-linear models results
setwd(Dir2)
Info1 <- read.csv(paste0("WasteAware_Info__",      Dep1, ".csv"))
Pred1 <- read.csv(paste0("WasteAware_Predicted__", Dep1, ".csv"))
Inte1 <- read.csv(paste0("WasteAware_Intervals__", Dep1, ".csv"))

Info2 <- read.csv(paste0("WasteAware_Info__",      Dep2, ".csv"))
Pred2 <- read.csv(paste0("WasteAware_Predicted__", Dep2, ".csv"))
Inte2 <- read.csv(paste0("WasteAware_Intervals__", Dep2, ".csv"))


###--------------------------------------------------------------------------###
### 5.  Dep1 - prepare data for best model panel:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 5.1.  Identify the top model and set the axis and strip labels.

Info1          <- Info1[order(Info1$Rank_1PerExp), ]
Info1          <- Info1[1, ]

# strip label
Info1$Label    <- paste0(
  Info1$ModName, 
  "\\nN: ",
  Info1$NumCase, 
  "\\nRMSE: ", 
  round(Info1$RMSE, 2), 
  "\\nAICc weight: ",
  round(Info1$AICcWt_1PerExp, 2),
  "  \\n")

Info1 <- droplevels(Info1)

# Exp/formulas combo included
ExpIn1 <- as.character(Info1$Exp___ModName)


###--------------------------------------------------------------------------###
### 5.2   Subset plotting data for observed and predcited. Set fill and border colors. 

Pred1             <- Pred1[Pred1$Exp___ModName %in% ExpIn1, ]
Pred1$PointColor  <- IncCat[match(Pred1$IncomeCategory, IncCat$IncomeCategory), "Color"]
Pred1$PointBorder <- IncCat[match(Pred1$IncomeCategory, IncCat$IncomeCategory), "Border"]
Pred1             <- droplevels(Pred1)

###--------------------------------------------------------------------------###
### 5.3.  Subset plotting data for confidence and prediction intervals. 

Inte1 <- Inte1[Inte1$Exp___ModName %in% ExpIn1, ]
Inte1 <- droplevels(Inte1)

###--------------------------------------------------------------------------###
### 5.4.  Get y axis limits in the plot.

P0    <- ggplot(data = Pred1)
P0    <- P0 + geom_point(aes(x = ExpValue, y = DepValue))
LimX1 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
if(DoLog){
  LimX1[1] <- 0.9 * min(Pred1$ExpValue)
}
if(!DoLog){
  LimX1[1] <- 0
}

rm(P0)

###--------------------------------------------------------------------------###
### 5.5.  Arrange title for facet.

Info1$Facet <- Label1
Pred1$Facet <- Label1
Inte1$Facet <- Label1

###--------------------------------------------------------------------------###
### 5.6.  Set Y range manually if needed and adjust prediction/confidence polygons. 

Inte1$PredLow_Line    <- Inte1$PredLow
Inte1$PredLow_Ribbon  <- Inte1$PredLow
Inte1$PredHigh_Line   <- Inte1$PredHigh
Inte1$PredHigh_Ribbon <- Inte1$PredHigh

Inte1$ConfLow_Line    <- Inte1$ConfLow
Inte1$ConfLow_Ribbon  <- Inte1$ConfLow
Inte1$ConfHigh_Line   <- Inte1$ConfHigh
Inte1$ConfHigh_Ribbon <- Inte1$ConfHigh

# prediction interval ribbons
# Inte1[Inte1$PredLow_Ribbon  < RangeY1[1], "PredLow_Ribbon"]  <- RangeY1[1]
# Inte1[Inte1$PredHigh_Ribbon > RangeY1[2], "PredHigh_Ribbon"] <- RangeY1[2]
# 
# # confidence interval ribbons
# Inte1[Inte1$ConfLow_Ribbon  < RangeY1[1], "ConfLow_Ribbon"]  <- RangeY1[1]
# Inte1[Inte1$ConfHigh_Ribbon > RangeY1[2], "ConfHigh_Ribbon"] <- RangeY1[2]


###--------------------------------------------------------------------------###
### 6.  Dep2 - prepare data for best model panel:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 6.1.  Identify the top model and set the axis and strip labels.

Info2          <- Info2[order(Info2$Rank_1PerExp), ]
Info2          <- Info2[1, ]

# strip label
Info2$Label    <- paste0(
  Info2$ModName, 
  "\\nN: ",
  Info2$NumCase, 
  "\\nRMSE: ", 
  round(Info2$RMSE, 2), 
  "\\nAICc weight: ",
  round(Info2$AICcWt_1PerExp, 2),
  "  \\n") 

Info2 <- droplevels(Info2)

# Exp/formulas combo included
ExpIn2 <- as.character(Info2$Exp___ModName)

###--------------------------------------------------------------------------###
### 6.2   Subset plotting data for observed and predicted. Set fill and border colors. 

Pred2             <- Pred2[Pred2$Exp___ModName %in% ExpIn2, ]
Pred2$PointColor  <- IncCat[match(Pred2$IncomeCategory, IncCat$IncomeCategory), "Color"]
Pred2$PointBorder <- IncCat[match(Pred2$IncomeCategory, IncCat$IncomeCategory), "Border"]
Pred2             <- droplevels(Pred2)

###--------------------------------------------------------------------------###
### 6.3.  Subset plotting data for confidence and prediction intervals. 

Inte2 <- Inte2[Inte2$Exp___ModName %in% ExpIn2, ]
Inte2 <- droplevels(Inte2)

###--------------------------------------------------------------------------###
### 6.4.  Get y axis limits in the plot.

P0    <- ggplot(data = Pred2)
P0    <- P0 + geom_point(aes(x = ExpValue, y = DepValue))
LimX2 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

###--------------------------------------------------------------------------###
### 6.5.  Arrange title for facet.

Info2$Facet <- Label2
Pred2$Facet <- Label2
Inte2$Facet <- Label2

###--------------------------------------------------------------------------###
### 6.6.  Set Y range manually if needed and adjust prediction/confidence polygons. 

Inte2$PredLow_Line    <- Inte2$PredLow
Inte2$PredLow_Ribbon  <- Inte2$PredLow
Inte2$PredHigh_Line   <- Inte2$PredHigh
Inte2$PredHigh_Ribbon <- Inte2$PredHigh

Inte2$ConfLow_Line    <- Inte2$ConfLow
Inte2$ConfLow_Ribbon  <- Inte2$ConfLow
Inte2$ConfHigh_Line   <- Inte2$ConfHigh
Inte2$ConfHigh_Ribbon <- Inte2$ConfHigh

# # prediction interval ribbons
# Inte2[Inte2$PredLow_Ribbon  < RangeY2[1], "PredLow_Ribbon"]  <- RangeY2[1]
# Inte2[Inte2$PredHigh_Ribbon > RangeY2[2], "PredHigh_Ribbon"] <- RangeY2[2]
# 
# # confidence interval ribbons
# Inte2[Inte2$ConfLow_Ribbon  < RangeY2[1], "ConfLow_Ribbon"]  <- RangeY2[1]
# Inte2[Inte2$ConfHigh_Ribbon > RangeY2[2], "ConfHigh_Ribbon"] <- RangeY2[2]


###--------------------------------------------------------------------------###
### 7.  Trade-off - prepare data for the plot:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 7.1.  Subset data from raw WasteAware data.  

Data3 <- Data[, c("casename",
                  "Country", 
                  "IncomeCategory",
                  "City",
                  Dep1,
                  Dep2)]

Data3 <- Data3[!is.na(Data3[, Dep1]), ]
Data3 <- Data3[!is.na(Data3[, Dep2]), ]

names(Data3)[5] <- "Dep1"
names(Data3)[6] <- "Dep2"

###--------------------------------------------------------------------------###
### 7.2.  Factor and set levels and labels of the income category.

Data3$IncomeCategory <- factor(Data3$IncomeCategory, 
                               levels = IncCat[1:4, "IncomeCategory"], 
                               labels = IncCat[1:4, "Label"])

###--------------------------------------------------------------------------###
### 7.3.  Add label for strip and points fill/border colors. 

Data3$Facet       <- Label3
Data3$PointColor  <- IncCat[match(Data3$IncomeCategory, IncCat$IncomeCategory), "Color"]
Data3$PointBorder <- IncCat[match(Data3$IncomeCategory, IncCat$IncomeCategory), "Border"]

###--------------------------------------------------------------------------###
### 7.4.  Set ranges. 

LimDep1 <- RangeY1
LimDep2 <- RangeY2


###--------------------------------------------------------------------------###
### 8.  Prepare data for the boxplots:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 8.1.  Subset data from raw WasteAware data. 

D1 <- Data[!is.na(Data[, Dep1]), ]
D2 <- Data[!is.na(Data[, Dep2]), ]

D1 <- D1[, c("casename", "Country", "IncomeCategory", "City",
             Dep1)]
D2 <- D2[, c("casename", "Country", "IncomeCategory", "City",
             Dep2)]

###--------------------------------------------------------------------------###
### 8.2.  Duplicate for the 'Overall' cateogry and bind.

T1                <- D1
T2                <- D2

T1$IncomeCategory <- "Overall"
T2$IncomeCategory <- "Overall"

D1 <- rbind(D1, T1)
D2 <- rbind(D2, T2)
rm(T1, T2)

###--------------------------------------------------------------------------###
### 8.3.  Add column with the dependent variable and rename the Dependent column. 

D1$Dep <- Dep1
D2$Dep <- Dep2

names(D1)[names(D1) == Dep1] <- "DepValue"
names(D2)[names(D2) == Dep2] <- "DepValue"

###--------------------------------------------------------------------------###
### 8.4.  Add facet info.

D1$Facet <- Label4[1]
D2$Facet <- Label4[2]

###--------------------------------------------------------------------------###
### 8.5.  Bind.

Data4 <- rbind(D1, D2)

###--------------------------------------------------------------------------###
### 8.6.  Factor and set levels and labels of the income category.

Data4$IncomeCategory <- factor(Data4$IncomeCategory, 
                               levels = IncCat$IncomeCategory, 
                               labels = IncCat$Label)

###--------------------------------------------------------------------------###
### 8.7.  Add points fill/border colors. 

Data4$PointColor  <- IncCat[match(Data4$IncomeCategory, IncCat$IncomeCategory), "Color"]
Data4$PointBorder <- IncCat[match(Data4$IncomeCategory, IncCat$IncomeCategory), "Border"]

###--------------------------------------------------------------------------###
### 8.8.  Get y axis limits in the plot.

P0    <- ggplot(data = Data4)
P0    <- P0 + geom_boxplot(aes(x = IncomeCategory, y = DepValue))
LimX4 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

###--------------------------------------------------------------------------###
### 8.9.  Create a color and border vector for the scales. 

ColVec4        <- as.character(IncCat$Color)
names(ColVec4) <- IncCat$Label
BorVec4        <- as.character(IncCat$Border)
names(BorVec4) <- IncCat$Label

###--------------------------------------------------------------------------###
### 8.10. Create the dataframes for the shades.

Shade4 <- rbind(Shade1, Shade2)

B1 <- Shade4[Shade4$Name == "B1", ]
B2 <- Shade4[Shade4$Name == "B2", ]
B3 <- Shade4[Shade4$Name == "B3", ]
B4 <- Shade4[Shade4$Name == "B4", ]
B5 <- Shade4[Shade4$Name == "B5", ]

B1 <- droplevels(B1)
B2 <- droplevels(B2)
B3 <- droplevels(B3)
B4 <- droplevels(B4)
B5 <- droplevels(B5)

###--------------------------------------------------------------------------###
### 8.11. Set ranges. 

LimBoth <- range(LimDep1, LimDep2)


###--------------------------------------------------------------------------###
### 9.  Dep1 - Plot best model panel:
###--------------------------------------------------------------------------###


if(DoLog){Breaks1 <- c(500, 1000, 5000, 10000, 50000)}
if(!DoLog){Breaks1 <- c(0, 10000, 20000, 30000, 40000, 50000)}

###--------------------------------------------------------------------------###
### 9.1.  Start plot.

P1 <- ggplot(data = Pred1)

###--------------------------------------------------------------------------###
### 9.2.  Background traffic light shades.

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[1, "Ymin"],
                         ymax = Shade1[1, "Ymax"]),
                     alpha = Shade1[1, "Alpha"],
                     fill  = Shade1[1, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[2, "Ymin"],
                         ymax = Shade1[2, "Ymax"]),
                     alpha = Shade1[2, "Alpha"],
                     fill  = Shade1[2, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[3, "Ymin"],
                         ymax = Shade1[3, "Ymax"]),
                     alpha = Shade1[3, "Alpha"],
                     fill  = Shade1[3, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[4, "Ymin"],
                         ymax = Shade1[4, "Ymax"]),
                     alpha = Shade1[4, "Alpha"],
                     fill  = Shade1[4, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[5, "Ymin"],
                         ymax = Shade1[5, "Ymax"]),
                     alpha = Shade1[5, "Alpha"],
                     fill  = Shade1[5, "Color"])

###--------------------------------------------------------------------------###
### 9.3.  Overlay grid.

P1 <- P1 + geom_vline(xintercept = Breaks1, 
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)
P1 <- P1 + geom_hline(yintercept = seq(0, 100, by=25),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 9.4.  Prediction interval, ribbon + lower and upper lines.

P1 <- P1 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = PredLow_Ribbon,
                           ymax = PredHigh_Ribbon),
                       fill  = FillPredIntervals, 
                       alpha = Alpha1,
                       data  = Inte1)
P1 <- P1 + geom_line(aes(x = ExpValue, 
                         y = PredLow_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte1)
P1 <- P1 + geom_line(aes(x = ExpValue, 
                         y = PredHigh_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte1)

###--------------------------------------------------------------------------###
### 9.5.  Confidence interval, ribbon + lower and upper lines.

P1 <- P1 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = ConfLow_Ribbon,
                           ymax = ConfHigh_Ribbon),
                       fill  = FillConfIntervals, 
                       alpha = Alpha2,
                       data  = Inte1)
P1 <- P1 + geom_line(aes(x = ExpValue, 
                         y = ConfLow_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte1)
P1 <- P1 + geom_line(aes(x = ExpValue, 
                         y = ConfHigh_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte1)

###--------------------------------------------------------------------------###
### 9.6.  Model info text.

P1 <- P1 + annotation_custom(gTree(children=gList(textGrob(Info1$Label, 
                                                           hjust = 0,
                                                           vjust = 0,
                                                           x     = LocTextMod[1],
                                                           y     = LocTextMod[2], 
                                                           gp    = gpar(cex = SizeTextMod))))) 

###--------------------------------------------------------------------------###
### 9.7.  Predicted values - line + points.

P1 <- P1 + geom_line(aes(x = ExpValue, 
                         y = Pred),
                     col  = ColPred, 
                     size = SizeLines, 
                     data = Pred1)
P1 <- P1 + geom_point(aes(x = ExpValue, 
                          y = DepValue),
                      size  = SizePoint1, 
                      color = Pred1$PointBorder,
                      fill  = Pred1$PointColor,
                      shape = ShapePoint, 
                      data  = Pred1)

###--------------------------------------------------------------------------###
### 9.8.  Scale axes.

if(DoLog){
  P1 <- P1 + scale_x_continuous(limits= LimX1,   expand= c(0,0), breaks = Breaks1, trans = "log10")
  P1 <- P1 + annotation_logticks(base = 10, sides = "b", scaled = T)
}
if(!DoLog){
  P1 <- P1 + scale_x_continuous(limits= LimX1,   expand= c(0,0), breaks = Breaks1)
}
#P1 <- P1 + scale_y_continuous(limits= RangeY1, expand= c(0,0))


###--------------------------------------------------------------------------###
### 9.9.  Axis labels + title.

P1 <- P1 + xlab(LabX1) + ylab(LabY1) 

###--------------------------------------------------------------------------###
### 9.10.  Faceting.

P1 <- P1 + facet_wrap( ~ Facet)

P1 <- P1 + coord_cartesian(ylim = LimBoth, expand= c(0,0))
  
###--------------------------------------------------------------------------###
### 9.11. Arrange theme.

P1 <- P1 + theme_bw()
P1 <- P1+ theme(
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 10.  DeP2 - Plot best model panel:
###--------------------------------------------------------------------------###


Breaks2 <- seq(30, 80, 10)

###--------------------------------------------------------------------------###
### 10.1.  Start plot.

P2 <- ggplot(data = Pred2)

###--------------------------------------------------------------------------###
### 10.2.  Background traffic light shades. 

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade2[1, "Ymin"],
                         ymax = Shade2[1, "Ymax"]),
                     alpha = Shade2[1, "Alpha"],
                     fill  = Shade2[1, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade2[2, "Ymin"],
                         ymax = Shade2[2, "Ymax"]),
                     alpha = Shade2[2, "Alpha"],
                     fill  = Shade2[2, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade2[3, "Ymin"],
                         ymax = Shade2[3, "Ymax"]),
                     alpha = Shade2[3, "Alpha"],
                     fill  = Shade2[3, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade2[4, "Ymin"],
                         ymax = Shade2[4, "Ymax"]),
                     alpha = Shade2[4, "Alpha"],
                     fill  = Shade2[4, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade2[5, "Ymin"],
                         ymax = Shade2[5, "Ymax"]),
                     alpha = Shade2[5, "Alpha"],
                     fill  = Shade2[5, "Color"])

###--------------------------------------------------------------------------###
### 10.3.  Overlay grid.

P2 <- P2 + geom_vline(xintercept = Breaks2, 
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)
P2 <- P2 + geom_hline(yintercept = seq(0, 100, by=25),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 10.4.  Prediction interval, ribbon + lower and upper lines.

P2 <- P2 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = PredLow_Ribbon,
                           ymax = PredHigh_Ribbon),
                       fill  = FillPredIntervals, 
                       alpha = Alpha1,
                       data  = Inte2)
P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = PredLow_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte2)
P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = PredHigh_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte2)

###--------------------------------------------------------------------------###
### 10.5.  Confidence interval, ribbon + lower and upper lines.

P2 <- P2 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = ConfLow_Ribbon,
                           ymax = ConfHigh_Ribbon),
                       fill  = FillConfIntervals, 
                       alpha = Alpha2,
                       data  = Inte2)
P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = ConfLow_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte2)
P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = ConfHigh_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte2)

###--------------------------------------------------------------------------###
### 10.6.  Model info text.

P2 <- P2 + annotation_custom(gTree(children=gList(textGrob(Info2$Label, 
                                                           hjust = 0,
                                                           vjust = 0,
                                                           x     = LocTextMod[1],
                                                           y     = LocTextMod[2], 
                                                           gp    = gpar(cex = SizeTextMod))))) 


###--------------------------------------------------------------------------###
### 10.7.  Predicted values - line + points.

P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = Pred),
                     col  = ColPred, 
                     size = SizeLines, 
                     data = Pred2)
P2 <- P2 + geom_point(aes(x = ExpValue, 
                          y = DepValue),
                      size  = SizePoint2, 
                      color = Pred2$PointBorder,
                      fill  = Pred2$PointColor,
                      shape = ShapePoint, 
                      data  = Pred2)

###--------------------------------------------------------------------------###
### 10.8.  Scale axes.

P2 <- P2 + scale_x_continuous(limits= LimX2,   expand= c(0,0))
#P2 <- P2 + scale_y_continuous(limits= RangeY2, expand= c(0,0))
P2 <- P2 + coord_cartesian(ylim = LimBoth, expand= c(0,0))

###--------------------------------------------------------------------------###
### 10.9.  Axis labels + title.

P2 <- P2 + xlab(LabX2) + ylab(LabY2) 

###--------------------------------------------------------------------------###
### 10.10.  Faceting.

P2 <- P2 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 10.11. Arrange theme.

P2 <- P2 + theme_bw()
P2 <- P2+ theme(
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 11.  Trade-off - Plot panel:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 11.1.  Start plot.

P3 <- ggplot(data = Data3)

###--------------------------------------------------------------------------###
### 11.2.  Background traffic light shades. 

P3 <- P3 + geom_rect(aes(xmin = Shade3[5, "Xmin"],
                         xmax = Shade3[5, "Xmax"],
                         ymin = Shade3[5, "Ymin"],
                         ymax = Shade3[5, "Ymax"]),
                     alpha = Shade3[5, "Alpha"],
                     fill  = Shade3[5, "Color"])

P3 <- P3 + geom_rect(aes(xmin = Shade3[4, "Xmin"],
                         xmax = Shade3[4, "Xmax"],
                         ymin = Shade3[4, "Ymin"],
                         ymax = Shade3[4, "Ymax"]),
                     alpha = Shade3[4, "Alpha"],
                     fill  = Shade3[4, "Color"])

P3 <- P3 + geom_rect(aes(xmin = Shade3[3, "Xmin"],
                         xmax = Shade3[3, "Xmax"],
                         ymin = Shade3[3, "Ymin"],
                         ymax = Shade3[3, "Ymax"]),
                     alpha = Shade3[3, "Alpha"],
                     fill  = Shade3[3, "Color"])


P3 <- P3 + geom_rect(aes(xmin = Shade3[2, "Xmin"],
                         xmax = Shade3[2, "Xmax"],
                         ymin = Shade3[2, "Ymin"],
                         ymax = Shade3[2, "Ymax"]),
                     alpha = Shade3[2, "Alpha"],
                     fill  = Shade3[2, "Color"])

P3 <- P3 + geom_rect(aes(xmin = Shade3[1, "Xmin"],
                         xmax = Shade3[1, "Xmax"],
                         ymin = Shade3[1, "Ymin"],
                         ymax = Shade3[1, "Ymax"]),
                     alpha = Shade3[1, "Alpha"],
                     fill  = Shade3[1, "Color"])

###--------------------------------------------------------------------------###
### 11.3.  Overlay grid.

P3 <- P3 + geom_vline(xintercept = seq(0, 100, by=25), 
                      size       = SizeGrid, 
                      color      = ColGrid, 
                      linetype   = TypeGrid)
P3 <- P3 + geom_hline(yintercept = seq(0, 100, by=25),
                      size       = SizeGrid, 
                      color      = ColGrid, 
                      linetype   = TypeGrid)

###--------------------------------------------------------------------------###
### 11.4.  Unity line. 

P3 <- P3 + geom_abline(aes(slope     = 1, 
                           intercept = 0), 
                       linetype = TypeVLine, 
                       size     = SizeLines, 
                       color    = ColVLine)

###--------------------------------------------------------------------------###
### 11.5.  Points. 

P3 <- P3 + geom_point(aes(x = Dep1, 
                          y = Dep2), 
                      size  = SizePoint3, 
                      color = Data3$PointBorder,
                      fill  = Data3$PointColor,
                      shape = ShapePoint, 
                      data  = Data3)

###--------------------------------------------------------------------------###
### 11.6.  Scale axes.

P3 <- P3 + scale_x_continuous(limit = LimDep1, expand = c(0,0))
P3 <- P3 + scale_y_continuous(limit = LimDep2, expand = c(0,0))
P3 <- P3 + coord_fixed(ratio = 1)

###--------------------------------------------------------------------------###
### 11.7.  Axis labels + title.

P3 <- P3 + xlab(LabX3) + ylab(LabY3) 

###--------------------------------------------------------------------------###
### 11.8.  Faceting. 

P3 <- P3 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 11.9.  Arrange theme.

P3 <- P3 + theme_bw()
P3 <- P3 + theme(
  legend.position = "none",
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 12. Plot boxplot panel:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 12.1.  Start plot.

P4 <- ggplot(data = Data4)

###--------------------------------------------------------------------------###
### 12.2.  Background traffic light shades.

P4 <- P4 + geom_boxplot(aes(x = IncomeCategory, 
                            y = DepValue)) # added just to set the x scale

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Ymin,
                         ymax = Ymax),
                     alpha = B1[1,"Alpha"],
                     fill  = B1[1, "Color"], 
                     data  = B1)

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Ymin,
                         ymax = Ymax),
                     alpha = B2[1,"Alpha"],
                     fill  = B2[1, "Color"], 
                     data  = B2)

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Ymin,
                         ymax = Ymax),
                     alpha = B3[1,"Alpha"],
                     fill  = B3[1, "Color"], 
                     data  = B3)

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Ymin,
                         ymax = Ymax),
                     alpha = B4[1,"Alpha"],
                     fill  = B4[1, "Color"], 
                     data  = B4)

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Ymin,
                         ymax = Ymax),
                     alpha = B5[1,"Alpha"],
                     fill  = B5[1, "Color"], 
                     data  = B5)

###--------------------------------------------------------------------------###
### 12.3.  Overlay grid.

P4 <- P4 + geom_hline(yintercept = seq(0, 100, by=25),
                      size       = SizeGrid, 
                      color      = ColGrid, 
                      linetype   = TypeGrid)


###--------------------------------------------------------------------------###
### 12.4.  Vertical line between Income categories and overall. 

P4 <- P4 + geom_vline(aes(xintercept = 4.5), 
                      linetype = TypeVLine, 
                      size     = SizeLines, 
                      color    = ColVLine)

###--------------------------------------------------------------------------###
### 12.5.  Box plot. 

P4 <- P4 + geom_boxplot(aes(x     = IncomeCategory, 
                            fill  = IncomeCategory,
                            color = IncomeCategory,
                            y     = DepValue), 
                        data = Data4)

###--------------------------------------------------------------------------###
### 12.6.  Points over boxplot.

P4 <- P4 + geom_point(aes(x = IncomeCategory, 
                          y = DepValue), 
                      size  = SizePoint4, 
                      color = Data4$PointBorder,
                      fill  = Data4$PointColor,
                      shape = ShapePoint,
                      alpha = Alpha4,
                      data  = Data4)

###--------------------------------------------------------------------------###
### 12.7.  Scale colors, fill and axes.

P4 <- P4 + scale_fill_manual(values = ColVec4) 
P4 <- P4 + scale_color_manual(values = BorVec4) 
P4 <- P4 + scale_y_continuous(limit = LimBoth, expand = c(0,0))


###--------------------------------------------------------------------------###
### 12.8.  Axis labels + title.

P4 <- P4 + xlab(LabX4) + ylab(LabY4) 

###--------------------------------------------------------------------------###
### 12.9.  Faceting. 

P4 <- P4 + facet_grid(~ Facet)

###--------------------------------------------------------------------------###
### 12.10.  Arrange theme.

P4 <- P4 + theme_bw()
P4 <- P4 + theme(
  legend.position = "none",
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 13.  Arrange the panels 
###--------------------------------------------------------------------------###


G1 <- arrangeGrob(P1, P3,
                  P2, P4,
                  ncol = 2, 
                  nrow = 2,
                  widths  = PanWidth, 
                  heights = PanHeight)


###--------------------------------------------------------------------------###
### 14.  Write to console or file. 
###--------------------------------------------------------------------------###


if(!ToFile){
  plot.new()
  grid.draw(G1)
}

if(ToFile){
  setwd(Dir3)
  for(i in FormOut){
    if(i != ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = G1,
             width  = Dim[1], 
             height = Dim[2], 
             units  = c("cm"),
             dpi    = 600)
    }
    # use cairo for special characters. 
    if(i == ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = G1,
             width  = Dim[1], 
             height = Dim[2], 
             units  = c("cm"),
             dpi    = 600, 
             device = cairo_pdf)
    }
    
  }
  ###------------------------------------------------------------------------###
  ### clean
  rm(i)
  
}






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



