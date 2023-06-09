################################################################################
###                                                                          ###
###---- NAME:                                                                ###
###                                                                          ###
### Plot_Fig2_ModelSelection_RF.R                                            ###
###                                                                          ###
###---- AIM:                                                                 ###
###                                                                          ###
### Plot the AICc weights (non-linear regressions) and the variable          ###
### importance (conditional random forest) of each explanatory variable in   ###
### each WABI.                                                               ###
###                                                                          ###
###---- INPUT:                                                               ###
###                                                                          ###
### 1.  "Meta"     - excel sheet.                                            ###
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
### 2.    WasteAware_Info.csv                                                ###     
###     General info and model selection results for each dependent variable ###
###     against all explanatory variables and model formula combinations.    ### 
###   "Dep"           - Name of dependent.                                   ###
###   "Exp"           - Name of explanatory.                                 ###
###   "ModName"       - Name of the model.                                   ###
###   "Exp___ModName" - Exp and model name, seperated by "___".              ###
###   "Formula"       - Formula of the model.                                ###
###   "Converged"     - Logical, T if the model converged, F if not.         ###
###   "NumCase"       - Number of cases.                                     ###
###   "a1"            - First coefficient, NA otherwise.                     ###
###   "a2"            - Second coefficient, NA otherwise.                    ###
###   "a3"            - Third coefficient, NA otherwise.                     ###
###   "a4"            - Fourth coefficient, NA otherwise.                    ###
###   "Deviance"      - The deviance.                                        ###
###   "RMSE"	        - Root mean square error.                              ###
###   "SMAPE"         - Symmetric mean absolute percentage error.            ###
###   "nRMSE_sd"      - Normalized RMSE. RMSE divided by sd of the observed. ###
###   "nRMSE_maxmin"  - Normalized RMSE. RMSE divided by range of the        ###
###                     observed.                                            ###
###   "K"             - Number of paramters in the model (+1 for the error   ###
###                     term).                                               ###
###   "AICc"          - Akike Information Criteria corrected for small       ###
###                     sample size of the model.                            ###
###   "Delta_AICc"    - Delta AICc from the best model. lower values are     ###
###                     better.                                              ###
###   "ModelLik"      - The relative likelihood value of the model.          ###    
###   "AICcWt"        - The Akaike weights.better models have higher         ###
###                     weights. All models sum to 1.                        ###
###   "LL"            - the log-likelihood of each model.                    ###
###   "Rank"          - Model rank from best (lowest deltaAICc) to worse     ###
###                     (highewst DeltaAICc). Lower values are  better       ###
###                     models.                                              ###
###                                                                          ###
### 3.    VarImp_CRF.csv                                                     ###
###     Dataframe with the regular and conditional variable importance for   ###
###     each dependent variable.                                             ###
###   "Dep"         - Name of the dependent variable.                        ###
###   "Exp"         - Name of the explanatory variable.                      ###
###   "VarImpF"     - Unconditional variable importance.                     ### 
###   "VarImpT"     - conditional variable importance.                       ###
###   "RankVarImpF" - Rank of unconditional var.imp. (low = better).         ###
###   "RankVarImpT" - Rank of conditional var.imp. (low = better).           ###
###                                                                          ###
### 4.    DepList_CRF_rmse.csv                                               ###
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
###---- PROCESS:                                                             ###
###                                                                          ###
### 1.  Packages.                                                            ###
### 2.  Directories.                                                         ###
### 3.  Controls.                                                            ###
###    3.1.  General controls.                                               ###
###    3.2.  Axes labels.                                                    ###
###    3.3.  Dependent variable labels.                                      ###
###    3.4.  Color scheme.                                                   ###
###    3.5.  Controls - to console.                                          ###
###    3.6.  Controls - to file.                                             ###
### 4.  Read and stack data.                                                 ###
### 5.  Arrange the metadata values for Exp.                                 ###                               
### 6.  Arrange the AICc weights data from plotting.                         ###                                       
### 7.  Plot model selection.                                                ###                                  
### 8.  Arrange data and plot the conditional randomforest plot.             ###
### 9.  Legend plot.                                                         ###
### 10.  Arrange the grob.                                                   ###
### 11.  Write to console or file.                                           ###
###                                                                          ###
###---- OUTPUT:                                                              ###
###                                                                          ###
### 1.  Figure with two panels and a legend.                                 ###
### 2.  Top panel is the model selection AICc weights.                       ###
###    2.1.  fill color is according to formula.                             ###
### 3.  Lower panel is the conditional randomforest variable importance.     ###
### 4.  Panels contain RMSE values for the best non-linear model and the     ###
###     randomforest model.                                                  ###
### 5.  x and y axes flipped for clarity.                                    ###
###                                                                          ###
###---- NOTES:                                                               ###
###                                                                          ###
### 1.                                                                       ###
###                                                                          ###
###---- VERSIONS:                                                            ###
###                                                                          ###
### 1.  14 June 2019 - Original code.                                        ###
### 2.  16 Feb. 2022 -                                                       ###
###    2.1.  Changed input to the excel table.                               ###
###    2.2.  Set stringsAsFactors == T when reading data.                    ###
###                                                                          ###
###---- DETAILS:                                                             ###
###                                                                          ###
### 1.  Author:    Yoni Gavish <gavishyoni.gmail.com>                        ###
###                                                                          ###
################################################################################



###--------------------------------------------------------------------------###
### 1.  Packages.
###--------------------------------------------------------------------------###


# data
library(readxl)
library(tidyr)
library(dplyr)

# plotting
library(ggplot2)
library(extrafont)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)


###--------------------------------------------------------------------------###
### 2.  Directories.
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/ModelSelection"
Dir3 <- "./Analysis/RandomForest"
Dir4 <- "./Analysis/Figures"


###--------------------------------------------------------------------------###
### 3.  Controls.
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 3.1.  General controls.

ToFile   <- T
FormOut <- c(".png", ".pdf")     # formats of output files
#FormOut  <- ".png"                # formats of output files
Name     <- "VelisEtAL_Fig2_modelSelection_RF"

###--------------------------------------------------------------------------###
### 3.2.  Axes labels.

LabY1 <- "AICc weights"
LabY2 <- "Conditional permutation importance"
LimY1 <- c(0, 0.6)

###--------------------------------------------------------------------------###
### 3.3.  Dependent variable labels.

M1 <- data.frame(NewName = c("WastePerCapita",
                             "CollectionCoverage",
                             "WasteQualityCollection",
                             "ControlledDisposal",
                             "QualityEnvironmentalProtection"),
                 Label = c("I-0:\\nWaste generation\\nrate",
                           "I-1.1:\\nWaste collection\\ncoverage",
                           "I-1C:\\nQuality of waste\\ncollection service",
                           "I-2:\\nControlled recovery\\nand disposal",
                           "I-2E:\\nEnvironmental\\nprotection in I-2"))

###--------------------------------------------------------------------------###
### 3.4.  Color scheme.

FillBar <- "gray85"
ColBar  <- "gray20"
FillBar2 <- "#fdc086" # RF cbars

ModFill <- data.frame(ModName = c("Linear",
                                  "Poly2",
                                  "Logarithmic",
                                  "sigEmax",
                                  "RandomForest"), 
                      Label = c(" Linear                   ",
                                " 2nd order polynomial     ",
                                " Logarithmic              ",
                                " Sigmoid                                 ",
                                " Conditional random forest"), 
                      Fill = c("#B5FAE1",
                               "#88EDEF",
                               "#22C6E2",
                               "#070B4F",
                               FillBar2), 
                      Exp = letters[1:5], 
                      Value = c(1:5))

###--------------------------------------------------------------------------###
### 3.5.  Controls - to console.

if(!ToFile){
  
  PanWidth  <- 10 
  Panheight <- c(1, 12, 10)
  
  ##--------------------------------------------------------------------------##
  ## RMSE text  -size, x location factor, color
  SizeRMSE <- 5
  LocX     <- 0.8
  ColRMSE  <- "black"

  ##--------------------------------------------------------------------------##
  ## Axes text and title
  SizeTextX     <- 10
  SizeTextY     <- 10
  SizeTitX      <- 11
  SizeStripText <- 13
  VjustX        <- 0.2
  
}

###--------------------------------------------------------------------------###
### 3.6.  Controls - to file.

if(ToFile){
  ###------------------------------------------------------------------------###
  ### figure dimension/resolution
  Dim <- c(30, 16.5) # width/height cm
  
  PanWidth  <- 10 
  Panheight <- c(1, 12, 10)
  
  ##--------------------------------------------------------------------------##
  ## RMSE text  -size, x location factor, color
  SizeRMSE <- 3
  LocX     <- 0.8
  ColRMSE  <- "black"
  
  
  ##--------------------------------------------------------------------------##
  ## Axes text and title
  SizeTextX     <- 8
  SizeTextY     <- 9
  SizeTitX      <- 9
  SizeStripText <- 9
  VjustX        <- 0.2
  
}


###--------------------------------------------------------------------------###
### 4.  Read and stack data.
###--------------------------------------------------------------------------###


setwd(Dir1)

Meta <- read_excel(path = paste0(Dir1, "/VelisEtAl2022_WABIs_Input.xlsx"),
                   sheet = "Meta")

# change to data.frame and factors for back compatibility 
Meta <- as.data.frame(Meta, stringsAsFactors = T)

for(i in 1:ncol(Meta)){
  if(class(Meta[, i]) == "character"){
    Meta[, i] <- factor(Meta[, i])
  }
}


setwd(Dir2)
T1 <- list.files(Dir2,"WasteAware_Info__" )
for(i in T1){
  T2 <- read.csv(i, stringsAsFactors = T)
  if(exists("D1")){D1 <- rbind(D1, T2)} else {D1 <- T2}
}
rm(i, T1, T2)

setwd(Dir3)
D2 <- read.csv("VarImp_CRF.csv", stringsAsFactors = T)
D3 <- read.csv("DepList_CRF_rmse.csv", stringsAsFactors = T)

###--------------------------------------------------------------------------###
### 5.  Arrange the metadata values for Exp.                                          
###--------------------------------------------------------------------------###


M2 <- Meta[Meta$Role == "Explanatory", ]
M2 <- M2[rev(order(M2$OrderX)), ]
M2 <- droplevels(M2)


###--------------------------------------------------------------------------###
### 6.  Arrange the AICc weights data from plotting.                                        
###--------------------------------------------------------------------------###

# subset data
D1 <- D1[!is.na(D1$Rank_1PerExp), ]
D1 <- droplevels(D1)

# set factor levels
D1$Dep <- factor(D1$Dep, 
                   levels = M1$NewName, 
                   labels = M1$Label)

D1$Exp <- factor(D1$Exp, 
                   levels = M2$Name, 
                   labels = M2$Label)

D1$ModName <- factor(D1$ModName, 
                       levels = ModFill$ModName, 
                       labels = ModFill$Label)

FillVec        <- as.character(ModFill$Fill)
names(FillVec) <- ModFill$Label 


###--------------------------------------------------------------------------###
### 7.  Arrange RMSE  values for plotting .                                        
###--------------------------------------------------------------------------###


D3$Dep <- factor(D3$Name, 
                 levels = M1$NewName, 
                 labels = M1$Label)
D3$Value_NLS <- round(D3$rmse_Best_NonLinear, 2)
D3$Value_CRF <- round(D3$rmse, 2)

D3$X <- 1


###--------------------------------------------------------------------------###
### 7.  Plot model selection.                                       
###--------------------------------------------------------------------------###


P1 <- ggplot(D1)
P1 <- P1 + geom_col(aes(x    = Exp, 
                        fill = ModName,  
                        y    = AICcWt_1PerExp),
                    color = ColBar,
                    data  = D1)



P1 <- P1 + scale_fill_manual(values = FillVec)
P1 <- P1 + ylab(LabY1)
P1 <- P1 + labs(tag = "A)")
P1 <- P1 + facet_wrap(~ Dep, nrow =1 )
P1 <- P1 + coord_flip()
P1 <- P1 + theme_bw()
P1 <- P1 + theme(
  legend.position = "none",
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_blank(),
  strip.text   = element_text(size = SizeStripText)
)


T1 <- c(ggplot_build(P1)$layout$panel_params[[1]]$x.range[2],
        ggplot_build(P1)$layout$panel_params[[2]]$x.range[2],
        ggplot_build(P1)$layout$panel_params[[3]]$x.range[2],
        ggplot_build(P1)$layout$panel_params[[4]]$x.range[2],
        ggplot_build(P1)$layout$panel_params[[5]]$x.range[2])

D3$Y_NLS <-  LocX * T1


P1 <- P1 + geom_text(aes(x     = X, 
                         y     = Y_NLS, 
                         label = Value_NLS),
                     size  = SizeRMSE,
                     color = ColRMSE,
                     data  = D3)





###--------------------------------------------------------------------------###
### 8.  Arrange data and plot the conditional randomforest plot.
###--------------------------------------------------------------------------###

# set factor levels
D2$Dep <- factor(D2$Dep, 
                 levels = M1$NewName, 
                 labels = M1$Label)

D2$Exp <- factor(D2$Exp, 
                 levels = M2$Name, 
                 labels = M2$Label)


# plot
P2 <- ggplot(D2)
P2 <- P2 + geom_col(aes(x    = Exp, 
                        y    = VarImpT),
                    fill  = FillBar2,
                    color = ColBar,
                    data  = D2)

P2 <- P2 + ylab(LabY2)
P2 <- P2 + labs(tag = "B)")
P2 <- P2 + facet_wrap(~ Dep, nrow =1, scale= "free_x")
P2 <- P2 + coord_flip()
P2 <- P2 + theme_bw()
P2 <- P2 + theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x  = element_text(size = SizeTextX),
  axis.text.y  = element_text(size = SizeTextY),
  axis.title.x = element_text(size = SizeTitX, vjust = VjustX),
  axis.title.y = element_blank()
)



T1 <- c(ggplot_build(P2)$layout$panel_params[[1]]$x.range[2],
        ggplot_build(P2)$layout$panel_params[[2]]$x.range[2],
        ggplot_build(P2)$layout$panel_params[[3]]$x.range[2],
        ggplot_build(P2)$layout$panel_params[[4]]$x.range[2],
        ggplot_build(P2)$layout$panel_params[[5]]$x.range[2])

D3$Y_CRF <-  LocX * T1
P2 <- P2 + geom_text(aes(x     = X, 
                         y     = Y_CRF, 
                         label = Value_CRF), 
                     size  = SizeRMSE,
                     color = ColRMSE,
                     data  = D3)
                     

###--------------------------------------------------------------------------###
### 9.  Legend plot.
###--------------------------------------------------------------------------###


ModFill$ModName2 <- factor(ModFill$ModName, 
                           levels = ModFill$ModName, 
                           labels = ModFill$Label)
P3 <- ggplot(data = ModFill)
P3 <- P3 + geom_col(aes(x    = Exp, 
                        fill = ModName2,  
                        y    = Value),
                    color = ColBar,
                    data  = ModFill)
P3 <- P3 + scale_fill_manual(values = FillVec)
P3 <- P3 + theme_bw()
P3 <- P3 + theme(
  legend.position = "top",
  legend.title = element_blank()
)

P3 <- get_legend(P3)
P3 <- as_ggplot(P3)


###--------------------------------------------------------------------------###
### 10.  Arrange the grob.
###--------------------------------------------------------------------------###


G1 <- arrangeGrob(P3, P1, P2, 
                  ncol    = 1, 
                  nrow    = 3, 
                  widths  = PanWidth, 
                  heights = Panheight)


###--------------------------------------------------------------------------###
### 11.  Write to consule or file. 
###--------------------------------------------------------------------------###


if(!ToFile){
  plot.new()
  grid.draw(G1)
  
}

if(ToFile){
  setwd(Dir4)
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
###   --------------------------------------------------------------------   ###
###      --------------------------------------------------------------      ###
###                                                                          ###
###                                    END                                   ###
###                                                                          ###
###      --------------------------------------------------------------      ###
###   --------------------------------------------------------------------   ###
###--------------------------------------------------------------------------###
################################################################################




