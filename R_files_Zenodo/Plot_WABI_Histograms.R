################################################################################
###                                                                          ###
###---- NAME:                                                                ###
###                                                                          ###
### Plot_WABI_Histograms.R                                                   ###
###                                                                          ###
###---- AIM:                                                                 ###
###                                                                          ###
### Plot an histogram of each WABI.                                          ###
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
###---- PROCESS:                                                             ###
###                                                                          ###
### 1.  Packages.                                                            ###
### 2.  Directories.                                                         ###
### 3.  Controls.                                                            ###
### 4.  Read data.                                                           ###
### 5.  Arrange data for plotting.                                           ###
### 6.  Plot.                                                                ###
### 7.  Write to console or file.                                            ###
###                                                                          ###
###----- OUTPUT: ------------------------------------------------------------###
###                                                                          ###
### 1.  A figure with 5 panels, one for each WABI with an histogram.         ### 
###                                                                          ###
###----- NOTES --------------------------------------------------------------###
###                                                                          ###
### 1.                                                                       ###
###                                                                          ###
###----- VERSIONS -----------------------------------------------------------###
###                                                                          ###
### 1.  28 Feb. 2022 - Original code.                                        ###
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
library(tidyverse)


###--------------------------------------------------------------------------###
### 2.  Directories.                                         
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/Figures"


###--------------------------------------------------------------------------###
### 3.  Controls.                                         
###--------------------------------------------------------------------------###


###--------------------------------------------------------------------------###
### 3.1.  General controls.

ToFile <- T
FormOut <- c(".png", ".pdf")     # formats of output files
FormOut  <- ".png"                # formats of output files

Name <- "VelisEtAL_FigS1_WABi_histograms"
LabX <- "WABIs value"
LabY <- "Number of cities"


###--------------------------------------------------------------------------###
### 3.2.  WABIs names.

M1 <- data.frame(NewName = c("WastePerCapita",
                             "CollectionCoverage",
                             "WasteQualityCollection",
                             "ControlledDisposal",
                             "QualityEnvironmentalProtection"),
                 Label1 = c("I-0:\\nWaste generation\\nrate",
                            "I-1.1:\\nWaste collection\\ncoverage",
                            "I-1C:\\nQuality of waste\\ncollection service",
                            "I-2:\\nControlled recovery\\nand disposal",
                            "I-2E:\\nEnvironmental\\nprotection in I-2"))

###--------------------------------------------------------------------------###
### 3.3.  Controls - to console.

if(!ToFile){
  
  ###------------------------------------------------------------------------###
  ### figure dimension/resolution
  Dim <- c(20, 12) # width/height cm
  
  ###------------------------------------------------------------------------###
  ### Histograms - bins and colors
  NumBins  <- 10
  ColHist  <- "#070B4F"
  FillHist <- "#A4BCC2"
}

###--------------------------------------------------------------------------###
### 3.6.  Controls - to file.
if(ToFile){
  
  ###------------------------------------------------------------------------###
  ### figure dimension/resolution
  Dim <- c(18, 7) # width/height cm
  
  ###------------------------------------------------------------------------###
  ### Histograms - bins and colors
  NumBins  <- 10
  ColHist  <- "#070B4F"
  FillHist <- "#A4BCC2"
  
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


###--------------------------------------------------------------------------###
### 5.  Arrange data for plotting.                                         
###--------------------------------------------------------------------------###


Out <- Data %>% gather(key = "WABI", value = "Value", WastePerCapita:QualityEnvironmentalProtection)

Out <- Out[!is.na(Out$Value), ]

Out$WABI <- factor(Out$WABI, 
                   levels = M1$NewName, 
                   labels = M1$Label1)


###--------------------------------------------------------------------------###
### 6.  Plot.                                         
###--------------------------------------------------------------------------###
 

P1 <- ggplot(data = Out)
P1 <- P1+ geom_histogram(aes(x = Value), 
                         bins  = NumBins,
                         data  = Out,
                         color = ColHist,
                         fill  = FillHist)
P1 <- P1 + xlab(LabX) + ylab(LabY)
P1 <- P1 + facet_wrap(~WABI,ncol = 5, scales = "free_x")

P1 <- P1 + theme_bw()


###--------------------------------------------------------------------------###
### 7.  Write to console or file. 
###--------------------------------------------------------------------------###


if(!ToFile){
  plot(P1)
}

if(ToFile){
  setwd(Dir2)
  for(i in FormOut){
    if(i != ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = P1,
             width  = Dim[1], 
             height = Dim[2], 
             units  = c("cm"),
             dpi    = 600)
    }
    # use cairo for special characters. 
    if(i == ".pdf"){
      ggsave(file   = paste0(Name, i), 
             plot   = P1,
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






