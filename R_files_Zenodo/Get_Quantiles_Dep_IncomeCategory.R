################################################################################
###                                                                          ###
###---- NAME:                                                                ###
###                                                                          ###
### Get_Quantiles_Dep_IncomeCategory.R                                       ###
###                                                                          ###
###---- AIM:                                                                 ###
###                                                                          ###
### Get  the 25th, 50th (median) and 75th quantile of each WABI for all      ###
### cities together and for each income category ("L", "L-M", U-M", "H").    ###
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
###---- PROCESS:                                                             ###
###                                                                          ###
### 1.  Packages.                                                            ###
### 2.  Directories.                                                         ###
### 3.  Read data and change to data frame for back compatibility.           ###
### 4.  Get the WABI list and set the output file.                           ###
### 5.  Loop over all WABIs and calculate the quantiles for each.            ### 
### 6.  Write output to file.                                                ###
###                                                                          ###
###---- OUTPUT:                                                              ###
###                                                                          ###
### 2.    Quantiles_Dep_IncomeCategory.csv                                   ###
###     Dataframe with the 25th, 50th and 75th quantile for each WABI.       ###
###     Contains the following columns:                                      ###
###  "Role"   - The role of the variable, either "Case", "General", "WABIs", ###
###             or "Explanatory".                                            ###
###  "Number" - Column number in 'Data'.                                     ### 
###  "Name"   - Name of the variable.                                        ###
###  "Label"  - Label of h variable (for plotting).                          ###
###  "OrderY" - Order of plotting.                                           ###
###  "OrderX" - Order of plotting.                                           ###
###  "Shapiro_W" - Shapiro-Wilk test of normality statistic.                 ###
###  "Shapiro_P" - Shapiro-Wilk test of normality p value.                   ###
###  "KolSmir_D" - Kolmogorov-Smirnov test of normality statistic.           ###
###  "KolSmir_P" - Kolmogorov-Smirnov test of normality p value.             ###
###  "Q25"       - 25th quantile, all cities.                                ###
###  "Median"    - 50th quantile (median), all cities.                       ###
###  "Q75"       - 75th quantile, all cities.                                ###
###  "Q25_L"     - 25th quantile, Low income countries.                      ###
###  "Median_L"  - 50th quantile (median), Low income countries.             ###
###  "Q75_L"     - 75th quantile, Low income countries.                      ###
###  "Q25_LM"    - 25th quantile, Low-Medium income countries.               ###
###  "Median_LM" - 50th quantile (median), Low-Medium income countries.      ###
###  "Q75_LM"    - 75th quantile, Low-Medium income countries.               ###
###  "Q25_UM"    - 25th quantile, Upper-Medium income countries.             ###
###  "Median_UM" - 50th quantile (median), Upper-Medium income countries.    ###
###  "Q75_UM"    - 75th quantile, Upper-Medium income countries.             ###
###  "Q25_H"     - 25th quantile, High income countries.                     ###
###  "Median_H"  - 50th quantile (median), High income countries.            ###
###  "Q75_H"     - 75th quantile, High income countries.                     ###
###                                                                          ###
###---- NOTES:                                                               ###
###                                                                          ###
###                                                                          ###
###---- VERSIONS:                                                            ###
###                                                                          ###
### 1.  16 Feb 2022 - Original code from 2019 adapted to new data format.    ###
###                   Tests of normality added.                              ###
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


###--------------------------------------------------------------------------###
### 2.  Directories.
###--------------------------------------------------------------------------###


Dir1 <- "./Analysis/Input"
Dir2 <- "./Analysis/General"


###--------------------------------------------------------------------------###
### 3.  Read data and change to data frame for back compatibility.
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


###--------------------------------------------------------------------------###
### 4.  Get the WABI list and set the output file.
###--------------------------------------------------------------------------###


Out1 <- Meta[Meta$Role == "WABIs", ]
Out1$Shapiro_W <- NA
Out1$Shapiro_P <- NA
Out1$KolSmir_D <- NA
Out1$KolSmir_P <- NA
Out1$Q25       <- NA 
Out1$Median    <- NA
Out1$Q75       <- NA
Out1$Q25_L     <- NA 
Out1$Median_L  <- NA
Out1$Q75_L     <- NA
Out1$Q25_LM    <- NA 
Out1$Median_LM <- NA
Out1$Q75_LM    <- NA
Out1$Q25_UM    <- NA 
Out1$Median_UM <- NA
Out1$Q75_UM    <- NA
Out1$Q25_H     <- NA 
Out1$Median_H  <- NA
Out1$Q75_H     <- NA



###--------------------------------------------------------------------------###
### 5.  Loop over all WABIs and calculate the quantiles for each.
###--------------------------------------------------------------------------###


for(i in 1:nrow(Out1)){
  
  T_All <- Data[, c(as.character(Out1[i, "Name"]), "IncomeCategory" )]
  T_L   <- T_All[T_All$IncomeCategory == "L", ]
  T_LM  <- T_All[T_All$IncomeCategory == "L-M", ]
  T_UM  <- T_All[T_All$IncomeCategory == "U-M", ]
  T_H   <- T_All[T_All$IncomeCategory == "H", ]
  
  Out1[i, "Q25"      ] <- quantile(T_All[, 1], probs = 0.25, na.rm = T)
  Out1[i, "Median"   ] <- quantile(T_All[, 1], probs = 0.50, na.rm = T)
  Out1[i, "Q75"      ] <- quantile(T_All[, 1], probs = 0.75, na.rm = T)
  Out1[i, "Q25_L"    ] <- quantile(T_L[, 1],   probs = 0.25, na.rm = T) 
  Out1[i, "Median_L" ] <- quantile(T_L[, 1],   probs = 0.50, na.rm = T)
  Out1[i, "Q75_L"    ] <- quantile(T_L[, 1],   probs = 0.75, na.rm = T)
  Out1[i, "Q25_LM"   ] <- quantile(T_LM[, 1],  probs = 0.25, na.rm = T) 
  Out1[i, "Median_LM"] <- quantile(T_LM[, 1],  probs = 0.50, na.rm = T) 
  Out1[i, "Q75_LM"   ] <- quantile(T_LM[, 1],  probs = 0.75, na.rm = T) 
  Out1[i, "Q25_UM"   ] <- quantile(T_UM[, 1],  probs = 0.25, na.rm = T)  
  Out1[i, "Median_UM"] <- quantile(T_UM[, 1],  probs = 0.50, na.rm = T) 
  Out1[i, "Q75_UM"   ] <- quantile(T_UM[, 1],  probs = 0.75, na.rm = T) 
  Out1[i, "Q25_H"    ] <- quantile(T_H[, 1],   probs = 0.25, na.rm = T)  
  Out1[i, "Median_H" ] <- quantile(T_H[, 1],   probs = 0.50, na.rm = T)
  Out1[i, "Q75_H"    ] <- quantile(T_H[, 1],   probs = 0.75, na.rm = T)
  
  
  # Normality tests
  Out1[i, "Shapiro_W"] <- shapiro.test(T_All[, 1])$statistic
  Out1[i, "Shapiro_P"] <- shapiro.test(T_All[, 1])$p.value
  Out1[i, "KolSmir_D"] <- ks.test(T_All[, 1], 'pnorm')$statistic
  Out1[i, "KolSmir_P"] <- ks.test(T_All[, 1], 'pnorm')$p.value
  
  rm(T_All, T_L, T_LM, T_UM, T_H)
}


###--------------------------------------------------------------------------###
### 6.  Write output to file.
###--------------------------------------------------------------------------###


setwd(Dir2)
write.csv(Out1, "Quantiles_Dep_IncomeCategory.csv", row.names = F)




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











