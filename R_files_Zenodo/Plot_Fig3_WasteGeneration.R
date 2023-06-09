################################################################################
###                                                                          ###
###----- NAME: --------------------------------------------------------------###
###                                                                          ###
### Plot_Fig3_WasteGeneration.R                                              ###
###                                                                          ###
###                                                                          ###
###----- AIM: ---------------------------------------------------------------###
###                                                                          ###
### Plot the top 3 non-linear models for the waste generation WABI with the  ###
### traffic light bands as background. Add a boxplot based on income         ###
### category.                                                                ###
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
###----- PROCESS: -----------------------------------------------------------###
###                                                                          ###
### 1.  Packages.                                                            ###                                          
### 2.  Directories.                                                         ###                                          
### 3.  Controls.                                                            ###                                          
###    3.1.  General controls.                                               ### 
###    3.2.  Axis / Strip labels.                                            ### 
###    3.3.  Background shading (traffic light).                             ### 
###    3.4.  Color scheme - Income categories.                               ###  
###    3.5.  Color scheme - general.                                         ### 
###    3.6.  Controls - to console.                                          ### 
###    3.7.  Controls - to file.                                             ### 
### 4.  Read data.                                                           ###                               
### 5.  Prepare data for the top 3 models:                                   ###
###    5.1.  Identify the top 3 models.                                      ###
###    5.2.  Models performance info as text in the panels.                  ###
###    5.3   Subset plotting data for observed and predicted. Set fill and   ###
###          border colors.                                                  ###
###    5.4.  Subset plotting data for confidence and prediction intervals.   ###
###    5.5.  Adjust prediction/confidence polygons to range of Y values.     ###
###    5.6.  Get y axis limits in the plot.                                  ###
###    5.7.  Arrange title for facet.                                        ###
### 6.  Prepare data for the boxplot:                                        ###
###    6.1.  Subset data from raw WasteAware data.                           ###
###    6.2.  Duplicate for the 'Overall' category and bind.                  ###
###    6.3.  Factor and set levels and labels of the income category.        ###
###    6.4.  Add label for strip and points fill/border colors.              ### 
###    6.5.  Get y axis limits in the plot.                                  ###
###    6.6.  Create a color and border vector for the scales.                ### 
### 7.  Plot model rank 1: HDI.                                              ###
###    7.1.  Set x axis breaks and labels.                                   ###
###    7.2.  Start plot.                                                     ###
###    7.3.  Background trafic light shades.                                 ###
###    7.4.  Overlay grid.                                                   ###
###    7.5.  Model info text.                                                ###
###    7.6.  Prediction interval, ribbon + lower and upper lines.            ###
###    7.7.  Confidence interval, ribbon + lower and upper lines.            ###
###    7.8.  Predicted values - line + points.                               ###
###    7.9.  Scale axes.                                                     ###
###    7.10.  Axis labels + title.                                           ###
###    7.11.  Faceting.                                                      ###
###    7.12.  Arrange theme.                                                 ###
### 8.  Plot model rank 2: GNI.                                              ###
###     Procedure same as above.                                             ###
### 9.  Plot model rank 3: GDP.                                              ###
###     Procedure same as above.                                             ###
### 10.  Plot boxplot panel:                                                 ###
###    10.1.  Start plot.                                                    ###
###    10.2.  Background traffic light shades.                               ### 
###    10.3.  Overlay grid.                                                  ###
###    10.4.  Vertical line between Income categories and overall.           ### 
###    10.5.  Box plot.                                                      ###
###    10.6.  Points over boxplot.                                           ###
###    10.7.  SCale colors, fill and axes.                                   ###
###    10.8.  Axis labels + title.                                           ###
###    10.9.  Faceting.                                                      ###
###    10.10.  Arrange theme.                                                ###
### 11.  Arrange the four panels in a grob:                                  ###
### 12. Write to console or file.                                            ###
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### 1.  A figure with 4 panels:                                              ### 
###   1.1.  Left top:     GFP model (rank 1).                                ### 
###   1.2.  Right top:    HDI model (rank 2).                                ### 
###   1.3.  Left bottom:  GNI model (rank 3).                                ### 
###   1.4.  Right bottom: histogram in relation to income categories.        ###
### 2.  The models panels contain:                                           ###
###   2.1.  Observed points, color according to income category.             ###
###   2.2.  Predicted line, confidence and prediction intervals.             ###
###   2.3.  Background colors is the traffic light system.                   ###
###   2.4.  Text with model performance.                                     ###
###                                                                          ###
###----- NOTES --------------------------------------------------------------###
###                                                                          ###
### 1.                                                                       ###
###                                                                          ###
###----- VERSIONS -----------------------------------------------------------###
###                                                                          ###
### 1.  3 Jun. 2019 - Original code.                                         ###
### 2.  21 June 2019:                                                        ###
###    2.1.  Changed background lights to polygons.                          ###
###    2.2.  Added option to select the nmodel to plot (based on rank).      ###
### 3.  1 July 2019:                                                         ###
###    3.1.  Changed background lights back to rect. only 4 levels.          ###
###    3.2.  Plot top 3 models and not a single model.                       ###
###    3.3.  Range is set manually, and prediction intervals are cut if      ###
###          needed.                                                         ###
### 4.  16 Feb. 2022 -                                                       ###
###    4.1.  Changed input to the excel table.                               ###
###    4.2.  Set stringsAsFactors == T when reading data.                    ###
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

ToFile <- T
FormOut <- c(".png", ".pdf")     # formats of output files
#FormOut  <- ".png"                # formats of output files
Name     <- "VelisEtAL_Fig3_WasteGeneration"

RangeY <- c(-5, 1010)

###--------------------------------------------------------------------------###
### 3.2.  Axis / Strip labels.

LabX1  <- "USD '000 (log scale)"    
LabY1  <- bquote('I-0: Waste generation rate [Kg'%.%y^-1%.%p^-1*']')
Label1 <- "A)  GDP PPP per capita"  

LabX2  <- "HDI"   
LabY2  <- bquote('I-0: Waste generation rate [Kg'%.%y^-1%.%p^-1*']')  
Label2 <- "B)  Human development index"    

LabX3  <- "USD '000 (log scale)"    
LabY3  <- bquote('I-0: Waste generation rate [Kg'%.%y^-1%.%p^-1*']') 
Label3 <- "C)  GNI per capita"    # set below

LabX4  <- "Income category"
LabY4  <- bquote('I-0: Waste generation rate [Kg'%.%y^-1%.%p^-1*']')  
Label4 <- "D)  Graphical summary"

###--------------------------------------------------------------------------###
### 3.3.  Background shading (traffic light).

Shade1 <- data.frame(Name  = c("B1", "B2", "B3", "B4"),
                     Ymin  = c(RangeY[1], 
                               250, 
                               500, 
                               750), 
                     Ymax  = c(250, 
                               500, 
                               750,  
                               RangeY[2]), 
                     Color = c('#BBFFCF',
                               '#FDFFBF',
                               '#FFEDBF',
                               '#FFCFC0'))

###--------------------------------------------------------------------------###
### 3.4.  Color scheme - Income categories. 

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
### 3.5.  Color scheme - general. 

ColText           <- "black"
ColPred           <- "black"
ColConfIntervals  <- "gray50"
FillConfIntervals <- "gray50"
ColPredIntervals  <- "gray80"
FillPredIntervals <- "gray80"

###--------------------------------------------------------------------------###
### 3.6.  Controls - to console.

if(!ToFile){
  
  ##--------------------------------------------------------------------------##
  ## Panel width and height
  PanWidth  <- c(10, 0.2, 10)
  PanHeight <- c(10, 10)
  
  ##--------------------------------------------------------------------------##
  ## Grid
  SizeGrid <- 0.2
  ColGrid  <- "gray80"
  TypeGrid <- "dotted"

  ##--------------------------------------------------------------------------##
  ## Model text - size, location
  SizeTextMod <- 1
  LocTextMod  <- c(0.05, 0.95)
  
  ##--------------------------------------------------------------------------##
  ## Observed values Points - size, color, shape
  SizePoint1 <- 3
  ColPoints  <- "black"
  ShapePoint <- 21
  
  ##--------------------------------------------------------------------------##
  ## predcited values lines - size, 
  SizeLines  <- 1
  
  ##--------------------------------------------------------------------------##
  ## Ribbons:
  Alpha1 <- 0.5
  Alpha2 <- 0.6
  
  ##--------------------------------------------------------------------------##
  ## Axes text and title
  SizeTextX     <- 11
  SizeTextY     <- 11
  SizeTitX      <- 12
  SizeTitY      <- 12
  SizeStripText <- 15
  VjustX        <- 0.2
  
  ##--------------------------------------------------------------------------##
  ## Histogram - vertical line, color and type
  ColVLine  <- "gray40"
  TypeVLine <- "dotted"
}

###--------------------------------------------------------------------------###
### 3.7.  Controls - to file.

if(ToFile){
  ###------------------------------------------------------------------------###
  ### figure dimension/resolution
  Dim <- c(31, 30) # width/height cm
  
  ##--------------------------------------------------------------------------##
  ## Panel width and height
  PanWidth  <- c(10, 0.1, 10)
  PanHeight <- c(10, 10)
  
  ##--------------------------------------------------------------------------##
  ## Grid
  SizeGrid <- 0.3
  ColGrid  <- "gray75"
  TypeGrid <- "dotted"
  
  ##--------------------------------------------------------------------------##
  ## Model text - size
  SizeTextMod <- 1
  LocTextMod  <- c(0.05, 0.95)
  
  ##--------------------------------------------------------------------------##
  ## Observed values Points - size, color, shape
  SizePoint1 <- 3.0
  ColPoints  <- "black"
  ShapePoint <- 21
  
  ##--------------------------------------------------------------------------##
  ## predcited values lines - size, 
  SizeLines  <- 1
  
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
### 4.  Read Data.                                         
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
InfoMain  <- read.csv("WasteAware_Info__WastePerCapita.csv",      stringsAsFactors = T)
Predicted <- read.csv("WasteAware_Predicted__WastePerCapita.csv", stringsAsFactors = T)
Intervals <- read.csv("WasteAware_Intervals__WastePerCapita.csv", stringsAsFactors = T)


###--------------------------------------------------------------------------###
### 5.  Prepare data for the top 3 models:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 5.1.  Identify the top 3 models.

InfoMain <- InfoMain[order(InfoMain$Rank_1PerExp), ]
Info1 <- InfoMain[1, ]
Exp1  <- as.character(Info1$Exp)

Info2 <- InfoMain[2, ]
Exp2  <- as.character(Info2$Exp)

Info3 <- InfoMain[3, ]
Exp3  <- as.character(Info3$Exp)

###--------------------------------------------------------------------------###
### 5.2.  Models performance info as text in the panels. 

Info1$Label    <- paste0(
  Info1$ModName, 
  "\\nN: ",
  Info1$NumCase, 
  "\\nRMSE: ", 
  round(Info1$RMSE, 2), 
  "\\nAICc weight: ",
  round(Info1$AICcWt_1PerExp, 2))

Info2$Label    <- paste0(
  Info2$ModName, 
  "\\nN: ",
  Info2$NumCase, 
  "\\nRMSE: ", 
  round(Info2$RMSE, 2), 
  "\\nAICc weight: ",
  round(Info2$AICcWt_1PerExp, 2))

Info3$Label    <- paste0(
  Info3$ModName, 
  "\\nN: ",
  Info3$NumCase, 
  "\\nRMSE: ", 
  round(Info3$RMSE, 2), 
  "\\nAICc weight: ",
  round(Info3$AICcWt_1PerExp, 2))

Info1          <- droplevels(Info1)
Info2          <- droplevels(Info2)
Info3          <- droplevels(Info3)

# Exp/formulas combo included
ExpIn1 <- as.character(Info1$Exp___ModName)
ExpIn2 <- as.character(Info2$Exp___ModName)
ExpIn3 <- as.character(Info3$Exp___ModName)

###--------------------------------------------------------------------------###
### 5.3  Subset plotting data for observed and predicted. Set fill and border colors. 

Pred1             <- Predicted[Predicted$Exp___ModName %in% ExpIn1, ]
Pred1$PointColor  <- IncCat[match(Pred1$IncomeCategory, IncCat$IncomeCategory), "Color"]
Pred1$PointBorder <- IncCat[match(Pred1$IncomeCategory, IncCat$IncomeCategory), "Border"]
Pred1             <- droplevels(Pred1)

Pred2             <- Predicted[Predicted$Exp___ModName %in% ExpIn2, ]
Pred2$PointColor  <- IncCat[match(Pred2$IncomeCategory, IncCat$IncomeCategory), "Color"]
Pred2$PointBorder <- IncCat[match(Pred2$IncomeCategory, IncCat$IncomeCategory), "Border"]
Pred2             <- droplevels(Pred2)

Pred3             <- Predicted[Predicted$Exp___ModName %in% ExpIn3, ]
Pred3$PointColor  <- IncCat[match(Pred3$IncomeCategory, IncCat$IncomeCategory), "Color"]
Pred3$PointBorder <- IncCat[match(Pred3$IncomeCategory, IncCat$IncomeCategory), "Border"]
Pred3             <- droplevels(Pred3)

###--------------------------------------------------------------------------###
### 5.4.  Subset plotting data for confidence and prediction intervals. 

Inte1 <- Intervals[Intervals$Exp___ModName %in% ExpIn1, ]
Inte1 <- droplevels(Inte1)

Inte2 <- Intervals[Intervals$Exp___ModName %in% ExpIn2, ]
Inte2 <- droplevels(Inte2)

Inte3 <- Intervals[Intervals$Exp___ModName %in% ExpIn3, ]
Inte3 <- droplevels(Inte3)

###--------------------------------------------------------------------------###
### 5.5.  Adjust prediction/confidence polygons to range of Y values. 

Inte1$PredLow_Line    <- Inte1$PredLow
Inte1$PredLow_Ribbon  <- Inte1$PredLow
Inte1$PredHigh_Line   <- Inte1$PredHigh
Inte1$PredHigh_Ribbon <- Inte1$PredHigh

Inte1$ConfLow_Line    <- Inte1$ConfLow
Inte1$ConfLow_Ribbon  <- Inte1$ConfLow
Inte1$ConfHigh_Line   <- Inte1$ConfHigh
Inte1$ConfHigh_Ribbon <- Inte1$ConfHigh

Inte2$PredLow_Line    <- Inte2$PredLow
Inte2$PredLow_Ribbon  <- Inte2$PredLow
Inte2$PredHigh_Line   <- Inte2$PredHigh
Inte2$PredHigh_Ribbon <- Inte2$PredHigh

Inte2$ConfLow_Line    <- Inte2$ConfLow
Inte2$ConfLow_Ribbon  <- Inte2$ConfLow
Inte2$ConfHigh_Line   <- Inte2$ConfHigh
Inte2$ConfHigh_Ribbon <- Inte2$ConfHigh

Inte3$PredLow_Line    <- Inte3$PredLow
Inte3$PredLow_Ribbon  <- Inte3$PredLow
Inte3$PredHigh_Line   <- Inte3$PredHigh
Inte3$PredHigh_Ribbon <- Inte3$PredHigh

Inte3$ConfLow_Line    <- Inte3$ConfLow
Inte3$ConfLow_Ribbon  <- Inte3$ConfLow
Inte3$ConfHigh_Line   <- Inte3$ConfHigh
Inte3$ConfHigh_Ribbon <- Inte3$ConfHigh

# # prediction interval ribbons
# Inte1[Inte1$PredLow_Ribbon  < RangeY[1], "PredLow_Ribbon"]  <- RangeY[1]
# Inte1[Inte1$PredHigh_Ribbon > RangeY[2], "PredHigh_Ribbon"] <- RangeY[2]
# 
# Inte2[Inte2$PredLow_Ribbon  < RangeY[1], "PredLow_Ribbon"]  <- RangeY[1]
# Inte2[Inte2$PredHigh_Ribbon > RangeY[2], "PredHigh_Ribbon"] <- RangeY[2]
# 
# Inte3[Inte3$PredLow_Ribbon  < RangeY[1], "PredLow_Ribbon"]  <- RangeY[1]
# Inte3[Inte3$PredHigh_Ribbon > RangeY[2], "PredHigh_Ribbon"] <- RangeY[2]
# 
# # confidence interval ribbons
# Inte1[Inte1$ConfLow_Ribbon  < RangeY[1], "ConfLow_Ribbon"]  <- RangeY[1]
# Inte1[Inte1$ConfHigh_Ribbon > RangeY[2], "ConfHigh_Ribbon"] <- RangeY[2]
# 
# Inte2[Inte2$ConfLow_Ribbon  < RangeY[1], "ConfLow_Ribbon"]  <- RangeY[1]
# Inte2[Inte2$ConfHigh_Ribbon > RangeY[2], "ConfHigh_Ribbon"] <- RangeY[2]
# 
# Inte3[Inte3$ConfLow_Ribbon  < RangeY[1], "ConfLow_Ribbon"]  <- RangeY[1]
# Inte3[Inte3$ConfHigh_Ribbon > RangeY[2], "ConfHigh_Ribbon"] <- RangeY[2]

###--------------------------------------------------------------------------###
### 5.6.  Get y axis limits in the plot.

P0    <- ggplot(data = Pred1)
P0    <- P0 + geom_point(aes(x = ExpValue, y = DepValue))
LimX1 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

P0    <- ggplot(data = Pred2)
P0    <- P0 + geom_point(aes(x = ExpValue, y = DepValue))
LimX2 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

P0    <- ggplot(data = Pred3)
P0    <- P0 + geom_point(aes(x = ExpValue, y = DepValue))
LimX3 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

LimX1[1] <- 0.9 * min(Pred1$ExpValue)
LimX3[1] <- 0.9 * min(Pred3$ExpValue)

###--------------------------------------------------------------------------###
### 5.7.  Arrange title for facet.

Info1$Facet <- Label1
Pred1$Facet <- Label1
Inte1$Facet <- Label1

Info2$Facet <- Label2
Pred2$Facet <- Label2
Inte2$Facet <- Label2

Info3$Facet <- Label3
Pred3$Facet <- Label3
Inte3$Facet <- Label3


###--------------------------------------------------------------------------###
### 6.  Prepare data for the boxplot:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 6.1.  Subset data from raw WasteAware data.  

Data <- Data[!is.na(Data$WastePerCapita), ]
Data <- Data[, c("casename",
                 "Country", 
                 "IncomeCategory",
                 "City",
                 "WastePerCapita")]

###--------------------------------------------------------------------------###
### 6.2.  Duplicate for the 'Overall' category and bind.

D1                <- Data
D1$IncomeCategory <- "Overall"
Data              <- rbind(Data, D1)
rm(D1)

###--------------------------------------------------------------------------###
### 6.3.  Factor and set levels and labels of the income category.

Data$IncomeCategory <- factor(Data$IncomeCategory, 
                              levels = IncCat$IncomeCategory, 
                              labels = IncCat$Label)

###--------------------------------------------------------------------------###
### 6.4.  Add label for strip and points fill/border colors. 

Data$Facet       <- Label4
Data$PointColor  <- IncCat[match(Data$IncomeCategory, IncCat$IncomeCategory), "Color"]
Data$PointBorder <- IncCat[match(Data$IncomeCategory, IncCat$IncomeCategory), "Border"]

###--------------------------------------------------------------------------###
### 6.5.  Get y axis limits in the plot.

P0    <- ggplot(data = Data)
P0    <- P0 + geom_boxplot(aes(x = IncomeCategory, y = WastePerCapita))
LimX4 <- ggplot_build(P0)$layout$panel_params[[1]]$x.range
rm(P0)

###--------------------------------------------------------------------------###
### 6.6.  Create a color and border vector for the scales. 

ColVec        <- as.character(IncCat$Color)
names(ColVec) <- IncCat$Label
BorVec        <- as.character(IncCat$Border)
names(BorVec) <- IncCat$Label


###--------------------------------------------------------------------------###
### 7.  Plot model rank 1.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 7.1.  Set x axis breaks and labels. 

Breaks1 <- c(1000, 5000, 10000, 50000)

###--------------------------------------------------------------------------###
### 7.2.  Start plot. 

P1 <- ggplot(data = Pred1)

###--------------------------------------------------------------------------###
### 7.3.  Background trafic light shades. 

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[1, "Ymin"],
                         ymax = Shade1[1, "Ymax"]),
                     fill  = Shade1[1, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[2, "Ymin"],
                         ymax = Shade1[2, "Ymax"]),
                     fill  = Shade1[2, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[3, "Ymin"],
                         ymax = Shade1[3, "Ymax"]),
                     fill  = Shade1[3, "Color"])

P1 <- P1 + geom_rect(aes(xmin = LimX1[1],
                         xmax = LimX1[2],
                         ymin = Shade1[4, "Ymin"],
                         ymax = Shade1[4, "Ymax"]),
                     fill  = Shade1[4, "Color"])

###--------------------------------------------------------------------------###
### 7.4.  Overlay grid.

P1 <- P1 + geom_vline(xintercept = Breaks1, 
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)
P1 <- P1 + geom_hline(yintercept = seq(0, 1000, by=125),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 7.5.  Model info text.

P1 <- P1 + annotation_custom(gTree(children=gList(textGrob(Info1$Label, 
                                                           hjust = 0,
                                                           vjust = 1,
                                                           x     = LocTextMod[1],
                                                           y     = LocTextMod[2], 
                                                           gp    = gpar(cex = SizeTextMod)))))

###--------------------------------------------------------------------------###
### 7.6.  Prediction interval, ribbon + lower and upper lines.

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
### 7.7.  Confidence interval, ribbon + lower and upper lines.

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
### 7.8.  Predicted values - line + points.

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
### 7.9.  Scale axes.

P1 <- P1 + scale_x_continuous(limits= LimX1, expand= c(0,0), trans = "log10", breaks = Breaks1)
#P1 <- P1 + scale_y_continuous(limits= RangeY, expand= c(0,0))
P1 <- P1 + coord_cartesian(ylim = RangeY, expand= c(0,0))
P1 <- P1 + annotation_logticks(base = 10, sides = "b", scaled = T)

###--------------------------------------------------------------------------###
### 7.10.  Axis labels + title.

P1 <- P1 + xlab(LabX1) + ylab(LabY1) 

###--------------------------------------------------------------------------###
### 7.11.  Faceting.

P1 <- P1 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 7.12.  Arrange theme.

P1 <- P1 + theme_bw()
P1 <- P1+ theme(
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 8.  Plot model rank 2: GNI.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 8.1.  Set x axis breaks and labels. 

Breaks2 <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
LabelsX2 <- c(0.3, "" ,0.5, "",  0.7, "", 0.9)


###--------------------------------------------------------------------------###
### 8.2.  Start plot. 

P2 <- ggplot(data = Pred2)

###--------------------------------------------------------------------------###
### 8.3.  Background trafic light shades. 

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade1[1, "Ymin"],
                         ymax = Shade1[1, "Ymax"]),
                     fill  = Shade1[1, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade1[2, "Ymin"],
                         ymax = Shade1[2, "Ymax"]),
                     fill  = Shade1[2, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade1[3, "Ymin"],
                         ymax = Shade1[3, "Ymax"]),
                     fill  = Shade1[3, "Color"])

P2 <- P2 + geom_rect(aes(xmin = LimX2[1],
                         xmax = LimX2[2],
                         ymin = Shade1[4, "Ymin"],
                         ymax = Shade1[4, "Ymax"]),
                     fill  = Shade1[4, "Color"])

###--------------------------------------------------------------------------###
### 8.4.  Overlay grid.

P2 <- P2 + geom_vline(xintercept = Breaks2, 
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)
P2 <- P2 + geom_hline(yintercept = seq(0, 1000, by=125),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 8.5.  Model info text.

P2 <- P2 + annotation_custom(gTree(children=gList(textGrob(Info2$Label, 
                                                           hjust = 0,
                                                           vjust = 1,
                                                           x     = LocTextMod[1],
                                                           y     = LocTextMod[2], 
                                                           gp    = gpar(cex = SizeTextMod))))) 

###--------------------------------------------------------------------------###
### 8.6.  Prediction interval, ribbon + lower and upper lines.

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
### 8.7.  Confidence interval, ribbon + lower and upper lines.

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
### 8.8.  Predicted values - line + points.

P2 <- P2 + geom_line(aes(x = ExpValue, 
                         y = Pred),
                     col  = ColPred, 
                     size = SizeLines, 
                     data = Pred2)
P2 <- P2 + geom_point(aes(x = ExpValue, 
                          y = DepValue),
                      size  = SizePoint1, 
                      color = Pred1$PointBorder,
                      fill  = Pred1$PointColor,
                      shape = ShapePoint, 
                      data  = Pred2)

###--------------------------------------------------------------------------###
### 8.9.  Scale axes.

P2 <- P2 + scale_x_continuous(limits= LimX2, expand= c(0,0),  breaks = Breaks2,  labels = LabelsX2)
#P2 <- P2 + scale_y_continuous(limits= RangeY, expand= c(0,0))
P2 <- P2 + coord_cartesian(ylim = RangeY, expand= c(0,0))

###--------------------------------------------------------------------------###
### 8.10.  Axis labels + title.

P2 <- P2 + xlab(LabX2) + ylab(LabY2) 

###--------------------------------------------------------------------------###
### 8.11.  Faceting.

P2 <- P2 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 8.12.  Arrange theme.

P2 <- P2 + theme_bw()
P2 <- P2 + theme(
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 9.  Plot model rank 3: GDP.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 9.1.  Set x axis breaks and labels.

Breaks3 <- c(500, 1000, 5000, 10000, 50000)

###--------------------------------------------------------------------------###
### 9.2.  Start plot. 

P3 <- ggplot(data = Pred3)

###--------------------------------------------------------------------------###
### 9.3.  Background trafic light shades. 

P3 <- P3 + geom_rect(aes(xmin = LimX3[1],
                         xmax = LimX3[2],
                         ymin = Shade1[1, "Ymin"],
                         ymax = Shade1[1, "Ymax"]),
                     fill  = Shade1[1, "Color"])

P3 <- P3 + geom_rect(aes(xmin = LimX3[1],
                         xmax = LimX3[2],
                         ymin = Shade1[2, "Ymin"],
                         ymax = Shade1[2, "Ymax"]),
                     fill  = Shade1[2, "Color"])

P3 <- P3 + geom_rect(aes(xmin = LimX3[1],
                         xmax = LimX3[2],
                         ymin = Shade1[3, "Ymin"],
                         ymax = Shade1[3, "Ymax"]),
                     fill  = Shade1[3, "Color"])

P3 <- P3 + geom_rect(aes(xmin = LimX3[1],
                         xmax = LimX3[2],
                         ymin = Shade1[4, "Ymin"],
                         ymax = Shade1[4, "Ymax"]),
                     fill  = Shade1[4, "Color"])

###--------------------------------------------------------------------------###
### 9.4.  Overlay grid.

P3 <- P3 + geom_vline(xintercept = Breaks3, 
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)
P3 <- P3 + geom_hline(yintercept = seq(0, 1000, by=125),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 9.5.  Model info text.

P3 <- P3 + annotation_custom(gTree(children=gList(textGrob(Info3$Label, 
                                                           hjust = 0,
                                                           vjust = 1,
                                                           x     = LocTextMod[1],
                                                           y     = LocTextMod[2], 
                                                           gp    = gpar(cex = SizeTextMod))))) 



###--------------------------------------------------------------------------###
### 9.6.  Prediction interval, ribbon + lower and upper lines.

P3 <- P3 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = PredLow_Ribbon,
                           ymax = PredHigh_Ribbon),
                       fill  = FillPredIntervals, 
                       alpha = Alpha1,
                       data  = Inte3)
P3 <- P3 + geom_line(aes(x = ExpValue, 
                         y = PredLow_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte3)
P3 <- P3 + geom_line(aes(x = ExpValue, 
                         y = PredHigh_Line),
                     col  = ColPredIntervals, 
                     size = SizeLines, 
                     data = Inte3)

###--------------------------------------------------------------------------###
### 9.7.  Confidence interval, ribbon + lower and upper lines.

P3 <- P3 + geom_ribbon(aes(x    = ExpValue, 
                           ymin = ConfLow_Ribbon,
                           ymax = ConfHigh_Ribbon),
                       fill  = FillConfIntervals, 
                       alpha = Alpha2,
                       data  = Inte3)
P3 <- P3 + geom_line(aes(x = ExpValue, 
                         y = ConfLow_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte3)
P3 <- P3 + geom_line(aes(x = ExpValue, 
                         y = ConfHigh_Line),
                     col  = ColConfIntervals, 
                     size = SizeLines, 
                     data = Inte3)

###--------------------------------------------------------------------------###
### 9.8.  Predicted values - line + points.

P3 <- P3 + geom_line(aes(x = ExpValue, 
                         y = Pred),
                     col  = ColPred, 
                     size = SizeLines, 
                     data = Pred3)
P3 <- P3 + geom_point(aes(x = ExpValue, 
                          y = DepValue),
                      size  = SizePoint1, 
                      color = Pred1$PointBorder,
                      fill  = Pred1$PointColor,
                      shape = ShapePoint, 
                      data  = Pred3)

###--------------------------------------------------------------------------###
### 9.9.  Scale axes.

P3 <- P3 + scale_x_continuous(limits= LimX3, expand= c(0,0), trans = "log10", breaks = Breaks3)
#P3 <- P3 + scale_y_continuous(limits= RangeY, expand= c(0,0))
P3 <- P3 + coord_cartesian(ylim = RangeY, expand= c(0,0))
P3 <- P3 + annotation_logticks(base = 10, sides = "b", scaled = T)


###--------------------------------------------------------------------------###
### 9.10.  Axis labels + title.

P3 <- P3 + xlab(LabX3) + ylab(LabY3) 

###--------------------------------------------------------------------------###
### 9.11.  Faceting.

P3 <- P3 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 9.12.  Arrange theme.

P3 <- P3 + theme_bw()
P3 <- P3 + theme(
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_text(size = SizeTitY),
  strip.text   = element_text(size = SizeStripText)
)


###--------------------------------------------------------------------------###
### 10.  Plot boxplot panel:
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 10.1.  Start plot. 

P4 <- ggplot(data = Data)

###--------------------------------------------------------------------------###
### 10.2.  Background traffic light shades. 

P4 <- P4 + geom_boxplot(aes(x = IncomeCategory, 
                            y = WastePerCapita)) # added just to set the x scale

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Shade1[1, "Ymin"],
                         ymax = Shade1[1, "Ymax"]),
                     fill  = Shade1[1, "Color"])

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Shade1[2, "Ymin"],
                         ymax = Shade1[2, "Ymax"]),
                     fill  = Shade1[2, "Color"])

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Shade1[3, "Ymin"],
                         ymax = Shade1[3, "Ymax"]),
                     fill  = Shade1[3, "Color"])

P4 <- P4 + geom_rect(aes(xmin = LimX4[1],
                         xmax = LimX4[2],
                         ymin = Shade1[4, "Ymin"],
                         ymax = Shade1[4, "Ymax"]),
                     fill  = Shade1[4, "Color"])

###--------------------------------------------------------------------------###
### 10.3.  Overlay grid.

P4 <- P4 + geom_hline(yintercept = seq(0, 1000, by=125),
                      size     = SizeGrid, 
                      color    = ColGrid, 
                      linetype = TypeGrid)

###--------------------------------------------------------------------------###
### 10.4.  Vertical line between Income categories and overall. 

P4 <- P4 + geom_vline(aes(xintercept = 4.5), 
                      linetype = TypeVLine, 
                      size     = SizeLines, 
                      color    = ColVLine)

###--------------------------------------------------------------------------###
### 10.5.  Box plot. 

P4 <- P4 + geom_boxplot(aes(x     = IncomeCategory, 
                            fill  = IncomeCategory,
                            color = IncomeCategory,
                            y     = WastePerCapita), 
                        data = Data)

###--------------------------------------------------------------------------###
### 10.6.  Points over boxplot.

P4 <- P4 + geom_point(aes(x = IncomeCategory, 
                          y= WastePerCapita), 
                      size  = SizePoint1, 
                      color = Data$PointBorder,
                      fill  = Data$PointColor,
                      shape = ShapePoint, 
                      data  = Data)

###--------------------------------------------------------------------------###
### 10.7.  SCale colors, fill and axes.

P4 <- P4 + scale_fill_manual(values = ColVec) 
P4 <- P4 + scale_color_manual(values = BorVec) 
P4 <- P4 + scale_y_continuous(limits = RangeY, expand= c(0,0))

###--------------------------------------------------------------------------###
### 10.8.  Axis labels + title.

P4 <- P4 + xlab(LabX4) + ylab(LabY4) 

###--------------------------------------------------------------------------###
### 10.9.  Faceting.

P4 <- P4 + facet_wrap( ~ Facet)

###--------------------------------------------------------------------------###
### 10.10.  Arrange theme.

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
### 11.  Arrange the four panels in a grob:
###--------------------------------------------------------------------------###


G1 <- arrangeGrob(P1, nullGrob(), P2,
                  P3, nullGrob(), P4,
                  ncol    = 3,
                  nrow    = 2, 
                  widths  = PanWidth,
                  heights = PanHeight) 


###--------------------------------------------------------------------------###
### 12. Write to consule or file. 
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




















