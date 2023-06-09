---
  title: "Analysis of Substances used against disease"
author: "Richard Heal"
date: "June 17, 2020"
output: html_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
#############################################################################
#
# Load in the data
#
#############################################################################
DataFolder <- "./Data"
DataFile_2008 <- "ShrimpDataSet_Substances_2008_V6.csv"
DataFile_2016 <- "ShrimpDataSet_Substances_2016_V6.csv"

library(data.table)

RawData_2008 <- read.csv(file.path(DataFolder, DataFile_2008),
                         header = TRUE,
                         stringsAsFactors = FALSE)
RawData_2008 <- setDT(RawData_2008)
RawData_2008 <- RawData_2008[!is.na(FarmID), ]

RawData_2016 <- read.csv(file.path(DataFolder, DataFile_2016),
                         header = TRUE,
                         stringsAsFactors = FALSE)
RawData_2016 <- setDT(RawData_2016)

library(knitr)
library(kableExtra)
```

## Analysis of Substances used against Diseases Data

This document shows the analysis of the cleaned disease data set with the cleaned substance categories present.

IN THIS VERSION IT IS UNDERSTOOD THAT THERE IS NO LINK BETWEEN THE SUBSTANCE USED AND THE DISEASE. THEREFORE DATA IS JUST USED ASSUMING A RECORDED INSTANCE OF USE OF THE SUBSTANCE AT THE FARM FOR THAT SPECIES CULTURED.

## Preliminary Analysis of Shrimp substance use against disease data

In this section simple data and plots of the data are obtained to understand the data.

```{r Preliminary Analysis on Penaus monodon, include=FALSE}
#####################################################################
#
# Generate a table of the data - 
#
#####################################################################

########################################################
#
# Strip out just monodon records and create single datatable
#
#########################################################
Monodon_2008 <- RawData_2008[SpeciesID == "Penaeus monodon", ]
Monodon_2016 <- RawData_2016[SpeciesID == "Penaeus monodon", ]
Monodon_2008[, SurveyYear := 2008]
Monodon_2016[, SurveyYear := 2016]
Monodon_2016[, Certification := NULL]

library(knitr)
library(kableExtra)
DescriptiveData <- data.frame(
  DataYear = c("2008", "2016"),
  NumberOfFarms = c(uniqueN(RawData_2008[SpeciesID == "Penaeus monodon",], by = "FarmID"), uniqueN(RawData_2016[SpeciesID == "Penaeus monodon",], by = "FarmID")),
  NumberOfSmallFarms = c(uniqueN(RawData_2008[SpeciesID == "Penaeus monodon" & FarmScaleID == "Small", .(FarmID)], by = "FarmID"),
                         uniqueN(RawData_2016[SpeciesID == "Penaeus monodon" & FarmScaleID == "Small", .(FarmID)], by = "FarmID")),
  NumberOfMediumFarms = c(uniqueN(RawData_2008[SpeciesID == "Penaeus monodon" & FarmScaleID == "Medium", .(FarmID)], by = "FarmID"),
                          uniqueN(RawData_2016[SpeciesID == "Penaeus monodon" & FarmScaleID == "Medium", .(FarmID)], by = "FarmID")),
  NumberOfLargeFarms = c(uniqueN(RawData_2008[SpeciesID == "Penaeus monodon" & FarmScaleID == "Large", .(FarmID)], by = "FarmID"),
                         uniqueN(RawData_2016[SpeciesID == "Penaeus monodon" & FarmScaleID == "Large", .(FarmID)], by = "FarmID"))
)

DescriptiveDataTable <- kable(DescriptiveData, col.names = c("Year", "Farms", "Small", "Medium", "Large")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)

```

Here are the farms present in the survey:
  
  ### Details of the farms in the survey - INDIVIDUAL FARMS
  
  `r DescriptiveDataTable`

NOTE THAT THE NUMBERS IN THESE CHARTS ARE THE SEPARATE FARMS USING A SPECIFIC PRODUCT. SOME FARMS USE SEVERAL PRODUCTS AND SO THE TOTAL NUMBER OF FARMS MAY BE MORE THAT THE ACTUAL NUMBER OF INDIVIDUAL FARMS.

```{r NumberTreatmentsPerFarm, include=FALSE}
##########################################################################
#
# A simple metric of the number of different treatments used per farm
#
# FOR: Commercial Vs Non;
###########################################################################
library(knitr)
library(kableExtra)
CommercialNon_2008 <- merge(
  Monodon_2008[, .N, .(FarmID, Commercial)][Commercial == "No",.N, FarmID],
  Monodon_2008[, .N, .(FarmID, Commercial)][Commercial == "Yes",.N, FarmID],
  by = "FarmID",
  all = TRUE
)
CommercialNon_2008[is.na(CommercialNon_2008)] <- 0
CommercialNon_2008[N.x == 1 & N.y == 0, NonOnly := 1]
CommercialNon_2008[N.x == 0 & N.y == 1, CommOnly := 1]
CommercialNon_2008[N.x == 1 & N.y == 1, NonComm := 1]
CommercialNon_2008[is.na(CommercialNon_2008)] <- 0

CommercialNon_2016 <- merge(
  Monodon_2016[, .N, .(FarmID, Commercial)][Commercial == "No",.N, FarmID],
  Monodon_2016[, .N, .(FarmID, Commercial)][Commercial == "Yes",.N, FarmID],
  by = "FarmID",
  all = TRUE
)
CommercialNon_2016[is.na(CommercialNon_2016)] <- 0
CommercialNon_2016[N.x == 1 & N.y == 0, NonOnly := 1]
CommercialNon_2016[N.x == 0 & N.y == 1, CommOnly := 1]
CommercialNon_2016[N.x == 1 & N.y == 1, NonComm := 1]
CommercialNon_2016[is.na(CommercialNon_2016)] <- 0

CommercialNon_DF <- data.frame(
  ProductType = c("Commercial only", "Non-commercial only", "Both", "Total"),
  NumberFarms_2008 = c(sum(CommercialNon_2008$CommOnly),
                       sum(CommercialNon_2008$NonOnly),
                       sum(CommercialNon_2008$NonComm),
                       sum(CommercialNon_2008$CommOnly,CommercialNon_2008$NonOnly,CommercialNon_2008$NonComm )),
  NumberFarms_2016 = c(sum(CommercialNon_2016$CommOnly),
                       sum(CommercialNon_2016$NonOnly),
                       sum(CommercialNon_2016$NonComm),
                       sum(CommercialNon_2016$CommOnly,CommercialNon_2016$NonOnly,CommercialNon_2016$NonComm ))
)

ProductTypeDataTable <- kable(CommercialNon_DF, col.names = c("Product", "Farms 2008", "Farms 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)


#############################################################
# Work out some basic stats
#############################################################
Mean_NumberTreatmentsFarm_2008 <- mean(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Mean_NumberTreatmentsFarm_2016 <- mean(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2008 <- median(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2016 <- median(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2008 <- sum(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2016 <- sum(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2008 <- max(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2016 <- max(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2008 <- min(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2016 <- min(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)

STD_NumberTreatmentsFarm_2008 <- sd(Monodon_2008[,SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2016 <- sd(Monodon_2016[,SubstanceCategory, FarmID][,.N, FarmID]$N)
#############################################################
# Work out some basic stats - small
#############################################################
Mean_NumberTreatmentsFarm_2008_SMALL <- mean(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Mean_NumberTreatmentsFarm_2016_SMALL <- mean(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2008_SMALL <- median(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2016_SMALL <- median(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2008_SMALL <- sum(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2016_SMALL <- sum(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2008_SMALL <- max(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2016_SMALL <- max(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2008_SMALL <- min(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2016_SMALL <- min(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)

STD_NumberTreatmentsFarm_2008_SMALL <- sd(Monodon_2008[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2016_SMALL <- sd(Monodon_2016[FarmScaleID == "Small",SubstanceCategory, FarmID][,.N, FarmID]$N)
#############################################################
# Work out some basic stats FarmScaleID == "Medium"
#############################################################
Mean_NumberTreatmentsFarm_2008_MEDIUM <- mean(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Mean_NumberTreatmentsFarm_2016_MEDIUM <- mean(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2008_MEDIUM <- median(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2016_MEDIUM <- median(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2008_MEDIUM <- sum(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2016_MEDIUM <- sum(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2008_MEDIUM <- max(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2016_MEDIUM <- max(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2008_MEDIUM <- min(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2016_MEDIUM <- min(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2008_MEDIUM <- sd(Monodon_2008[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2016_MEDIUM <- sd(Monodon_2016[FarmScaleID == "Medium",SubstanceCategory, FarmID][,.N, FarmID]$N)
#############################################################
# Work out some basic stats
#############################################################
Mean_NumberTreatmentsFarm_2008_LARGE <- mean(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Mean_NumberTreatmentsFarm_2016_LARGE <- mean(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2008_LARGE <- median(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Median_NumberTreatmentsFarm_2016_LARGE <- median(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2008_LARGE <- sum(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Total_NumberTreatmentsFarm_2016_LARGE <- sum(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2008_LARGE <- max(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Max_NumberTreatmentsFarm_2016_LARGE <- max(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2008_LARGE <- min(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
Min_NumberTreatmentsFarm_2016_LARGE <- min(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2008_LARGE <- sd(Monodon_2008[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
STD_NumberTreatmentsFarm_2016_LARGE <- sd(Monodon_2016[FarmScaleID == "Large",SubstanceCategory, FarmID][,.N, FarmID]$N)
###########################################################################
#
# Box plot
#
###########################################################################
```

### Details of the products used by the farms in the survey - INDIVIDUAL FARMS

In 2008 there were `r Total_NumberTreatmentsFarm_2008` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2008` and a median of `r Median_NumberTreatmentsFarm_2008` which ranged from `r Max_NumberTreatmentsFarm_2008` to `r Min_NumberTreatmentsFarm_2008` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2008`.

In 2016 there were `r Total_NumberTreatmentsFarm_2016` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2016` and a median of `r Median_NumberTreatmentsFarm_2016` which ranged from `r Max_NumberTreatmentsFarm_2016` to `r Min_NumberTreatmentsFarm_2016` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2016`.  

#### Small farms

In 2008 there were `r Total_NumberTreatmentsFarm_2008_SMALL` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2008_SMALL` and a median of `r Median_NumberTreatmentsFarm_2008_SMALL` which ranged from `r Max_NumberTreatmentsFarm_2008_SMALL` to `r Min_NumberTreatmentsFarm_2008_SMALL` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2008_SMALL`.    

In 2016 there were `r Total_NumberTreatmentsFarm_2016_SMALL` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2016_SMALL` and a median of `r Median_NumberTreatmentsFarm_2016_SMALL` which ranged from `r Max_NumberTreatmentsFarm_2016_SMALL` to `r Min_NumberTreatmentsFarm_2016_SMALL` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2016_SMALL`.   

#### Medium farms

In 2008 there were `r Total_NumberTreatmentsFarm_2008_MEDIUM` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2008_MEDIUM` and a median of `r Median_NumberTreatmentsFarm_2008_MEDIUM` which ranged from `r Max_NumberTreatmentsFarm_2008_MEDIUM` to `r Min_NumberTreatmentsFarm_2008_MEDIUM` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2008_MEDIUM`.    

In 2016 there were `r Total_NumberTreatmentsFarm_2016_MEDIUM` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2016` and a median of `r Median_NumberTreatmentsFarm_2016_MEDIUM` which ranged from `r Max_NumberTreatmentsFarm_2016` to `r Min_NumberTreatmentsFarm_2016_MEDIUM` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2016_MEDIUM`.  

#### Large farms

In 2008 there were `r Total_NumberTreatmentsFarm_2008_LARGE` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2008_LARGE` and a median of `r Median_NumberTreatmentsFarm_2008_LARGE` which ranged from `r Max_NumberTreatmentsFarm_2008_LARGE` to `r Min_NumberTreatmentsFarm_2008_LARGE` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2008_LARGE`.    

In 2016 there were `r Total_NumberTreatmentsFarm_2016_LARGE` recorded treatments which is a mean of `r Mean_NumberTreatmentsFarm_2016_LARGE` and a median of `r Median_NumberTreatmentsFarm_2016_LARGE` which ranged from `r Max_NumberTreatmentsFarm_2016_LARGE` to `r Min_NumberTreatmentsFarm_2016_LARGE` recorded treatments per farm.
The standard deviation was `r STD_NumberTreatmentsFarm_2016_LARGE`.  

The table below lists whether a farm has used only non-commercial products, only commercial products, or both.

`r ProductTypeDataTable`

```{r NumberCategoriesPerFarm, include=FALSE}
##########################################################################
#
# A simple metric of the number of different treatments used per farm
#
# FOR: Number of categories (Subcategory)
###########################################################################
library(knitr)
library(kableExtra)
NumberCategories_2008 <- Monodon_2008[, .N, .(FarmID, SubstanceCategory)][,.N, FarmID]
NumberCategories_2016 <- Monodon_2016[, .N, .(FarmID, SubstanceCategory)][,.N, FarmID]

CommercialNon_DF <- data.frame(
  NumberCategories = c(1:6),
  NumberFarms_Cat2008 = c(NumberCategories_2008[N == 1, .N],
                          NumberCategories_2008[N == 2, .N],
                          NumberCategories_2008[N == 3, .N],
                          NumberCategories_2008[N == 4, .N],
                          NumberCategories_2008[N == 5, .N],
                          NumberCategories_2008[N == 6, .N]),
  NumberFarms_Cat2016 = c(NumberCategories_2016[N == 1, .N],
                          NumberCategories_2016[N == 2, .N],
                          NumberCategories_2016[N == 3, .N],
                          NumberCategories_2016[N == 4, .N],
                          NumberCategories_2016[N == 5, .N],
                          NumberCategories_2016[N == 6, .N]))

ProductCategoriesDataTable <- kable(CommercialNon_DF, col.names = c("Number Categories", "Farms 2008", "Farms 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)

```

### Details of the different category of products by the farms in the survey - INDIVIDUAL FARMS

The table below lists the number of different types of categories of products used by the farm.  The mean number of treatment categories used per farm were `r mean(NumberCategories_2008)` in 2008 compared to `r mean(NumberCategories_2008)` in 2016.

The types are `r unique(Monodon_2008[,.N, SubstanceCategory]$SubstanceCategory, Monodon_2016[,.N, SubstanceCategory]$SubstanceCategory)`.


`r ProductCategoriesDataTable`


```{r NumberSubstancesPerFarm, include=FALSE}
##########################################################################
#
# A simple metric of the number of different treatments used per farm
#
# FOR: Number of categories (Subcategory)
###########################################################################
library(knitr)
library(kableExtra)
NumberSubstances_2008 <- Monodon_2008[, .N, .(FarmID, SubstanceUsed)][,.N, FarmID]
NumberSubstances_2016 <- Monodon_2016[, .N, .(FarmID, SubstanceUsed)][,.N, FarmID]

NumberSubs_DF <- data.frame(
  NumberCategories = c(1:15),
  NumberFarms_Sub2008 = c(NumberSubstances_2008[N == 1, .N],
                          NumberSubstances_2008[N == 2, .N],
                          NumberSubstances_2008[N == 3, .N],
                          NumberSubstances_2008[N == 4, .N],
                          NumberSubstances_2008[N == 5, .N],
                          NumberSubstances_2008[N == 6, .N],
                          NumberSubstances_2008[N == 7, .N],
                          NumberSubstances_2008[N == 8, .N],
                          NumberSubstances_2008[N == 9, .N],
                          NumberSubstances_2008[N == 10, .N],
                          NumberSubstances_2008[N == 11, .N],
                          NumberSubstances_2008[N == 12, .N],
                          NumberSubstances_2008[N == 13, .N],
                          NumberSubstances_2008[N == 14, .N],
                          NumberSubstances_2008[N == 15, .N]),
  NumberFarms_Sub2016 = c(NumberSubstances_2016[N == 1, .N],
                          NumberSubstances_2016[N == 2, .N],
                          NumberSubstances_2016[N == 3, .N],
                          NumberSubstances_2016[N == 4, .N],
                          NumberSubstances_2016[N == 5, .N],
                          NumberSubstances_2016[N == 6, .N],
                          NumberSubstances_2016[N == 7, .N],
                          NumberSubstances_2016[N == 8, .N],
                          NumberSubstances_2016[N == 9, .N],
                          NumberSubstances_2016[N == 10, .N],
                          NumberSubstances_2016[N == 11, .N],
                          NumberSubstances_2016[N == 12, .N],
                          NumberSubstances_2016[N == 13, .N],
                          NumberSubstances_2016[N == 14, .N],
                          NumberSubstances_2016[N == 15, .N]))

NumberSubstancesDataTable <- kable(NumberSubs_DF, col.names = c("Product", "Farms 2008", "Farms 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)

#############################################################################
# Basic stats
#############################################################################
Mean_NumberTreatmentsCATSFarm_2008 <- mean(Monodon_2008[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsCATSFarm_2016 <- mean(Monodon_2016[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2008 <- median(Monodon_2008[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2016 <- median(Monodon_2016[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2008 <- sum(Monodon_2008[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2016 <- sum(Monodon_2016[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2008 <- max(Monodon_2008[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2016 <- max(Monodon_2016[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2008 <- min(Monodon_2008[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2016 <- min(Monodon_2016[,.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)

#############################################################################
# Basic stats
#############################################################################
Mean_NumberTreatmentsCATSFarm_2008_SMALL <- mean(Monodon_2008[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsCATSFarm_2016_SMALL <- mean(Monodon_2016[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2008_SMALL <- median(Monodon_2008[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2016_SMALL <- median(Monodon_2016[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2008_SMALL <- sum(Monodon_2008[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2016_SMALL <- sum(Monodon_2016[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2008_SMALL <- max(Monodon_2008[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2016_SMALL <- max(Monodon_2016[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2008_SMALL <- min(Monodon_2008[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2016_SMALL <- min(Monodon_2016[FarmScaleID == "Small",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)

#############################################################################
# Basic stats
#############################################################################
Mean_NumberTreatmentsCATSFarm_2008_MEDIUM <- mean(Monodon_2008[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsCATSFarm_2016_MEDIUM <- mean(Monodon_2016[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2008_MEDIUM <- median(Monodon_2008[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2016_MEDIUM <- median(Monodon_2016[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2008_MEDIUM <- sum(Monodon_2008[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2016_MEDIUM <- sum(Monodon_2016[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2008_MEDIUM <- max(Monodon_2008[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2016_MEDIUM <- max(Monodon_2016[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2008_MEDIUM <- min(Monodon_2008[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2016_MEDIUM <- min(Monodon_2016[FarmScaleID == "Medium",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)

#############################################################################
# Basic stats - LARGE
#############################################################################
Mean_NumberTreatmentsCATSFarm_2008_LARGE <- mean(Monodon_2008[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsCATSFarm_2016_LARGE <- mean(Monodon_2016[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2008_LARGE <- median(Monodon_2008[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsCATSFarm_2016_LARGE <- median(Monodon_2016[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2008_LARGE <- sum(Monodon_2008[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsCATSFarm_2016_LARGE <- sum(Monodon_2016[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2008_LARGE <- max(Monodon_2008[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsCATSFarm_2016_LARGE <- max(Monodon_2016[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2008_LARGE <- min(Monodon_2008[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsCATSFarm_2016_LARGE <- min(Monodon_2016[FarmScaleID == "Large",.N,.(SubstanceCategory, FarmID)][,.N, .(SubstanceCategory,FarmID)][,.N, FarmID]$N)
```

### Details of the different category of products by the farms in the survey - INDIVIDUAL FARMS


In 2008 there were `r Total_NumberTreatmentsCATSFarm_2008_SMALL` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2008_SMALL` and a median of `r Median_NumberTreatmentsCATSFarm_2008_SMALL` which ranged from `r Max_NumberTreatmentsCATSFarm_2008_SMALL` to `r Min_NumberTreatmentsCATSFarm_2008_SMALL` recorded treatments per farm.

In 2016 there were `r Total_NumberTreatmentsCATSFarm_2016_SMALL` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2016_SMALL` and a median of `r Median_NumberTreatmentsCATSFarm_2016_SMALL` which ranged from `r Max_NumberTreatmentsCATSFarm_2016_SMALL` to `r Min_NumberTreatmentsCATSFarm_2016_SMALL` recorded treatments per farm.

#### Small

In 2008 there were `r Total_NumberTreatmentsCATSFarm_2008_MEDIUM` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2008_MEDIUM` and a median of `r Median_NumberTreatmentsCATSFarm_2008_MEDIUM` which ranged from `r Max_NumberTreatmentsCATSFarm_2008_MEDIUM` to `r Min_NumberTreatmentsCATSFarm_2008_MEDIUM` recorded treatments per farm.

In 2016 there were `r Total_NumberTreatmentsCATSFarm_2016_MEDIUM` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2016` and a median of `r Median_NumberTreatmentsCATSFarm_2016_MEDIUM` which ranged from `r Max_NumberTreatmentsCATSFarm_2016_MEDIUM` to `r Min_NumberTreatmentsCATSFarm_2016_MEDIUM` recorded treatments per farm.

#### Medium

In 2008 there were `r Total_NumberTreatmentsCATSFarm_2008_LARGE` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2008_LARGE` and a median of `r Median_NumberTreatmentsCATSFarm_2008_LARGE` which ranged from `r Max_NumberTreatmentsCATSFarm_2008_LARGE` to `r Min_NumberTreatmentsCATSFarm_2008_LARGE` recorded treatments per farm.

In 2016 there were `r Total_NumberTreatmentsCATSFarm_2016_LARGE` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2016_LARGE` and a median of `r Median_NumberTreatmentsCATSFarm_2016_LARGE` which ranged from `r Max_NumberTreatmentsCATSFarm_2016_LARGE` to `r Min_NumberTreatmentsCATSFarm_2016_LARGE` recorded treatments per farm.

#### Large

In 2008 there were `r Total_NumberTreatmentsCATSFarm_2008` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2008` and a median of `r Median_NumberTreatmentsCATSFarm_2008` which ranged from `r Max_NumberTreatmentsCATSFarm_2008` to `r Min_NumberTreatmentsCATSFarm_2008` recorded treatments per farm.

In 2016 there were `r Total_NumberTreatmentsCATSFarm_2016` recorded treatments which is a mean of `r Mean_NumberTreatmentsCATSFarm_2016` and a median of `r Median_NumberTreatmentsCATSFarm_2016` which ranged from `r Max_NumberTreatmentsCATSFarm_2016` to `r Min_NumberTreatmentsCATSFarm_2016` recorded treatments per farm.


The table below lists the number of different types of categories of products used by the farm.


`r NumberSubstancesDataTable`

```{r Box Plot Number of Substances Used, include=FALSE}
###########################################################################
# Box plot of the data for the number of substances used
###########################################################################
library(plotly)
NumberSubstances_2008 <- Monodon_2008[, .N, .(FarmID, SubstanceUsed)][,.N, FarmID]
NumberSubstances_2016 <- Monodon_2016[, .N, .(FarmID, SubstanceUsed)][,.N, FarmID]

BoxPlotNumberSubstances <- plot_ly() %>%
  add_trace(type = 'box',
            y = NumberCategories_2008$N,
            name = '2008') %>%
  add_trace(type = 'box',
            y = NumberCategories_2016$N,
            name = '2016') %>%
  layout(title = "Box Plot of the Number of Different Substances Used by Each Farm",
         showlegend = FALSE,
         xaxis = list(title = "Year", showline = TRUE),
         yaxis = list(title = "Number of Different Substances Used", showline = TRUE))

BoxPlotNumberSubstances
```

This information is shown below as a box plot.

`r BoxPlotNumberSubstances`

```{r Treatments per Farm, include=FALSE}
#################################################################
# Produce a table of the treatments per farm for P. monodon
#################################################################
library(knitr)
library(kableExtra)
FarmsUsingCategory_2008 <- Monodon_2008[, .N, .(FarmID, SubstanceCategory)][,.N, SubstanceCategory]
FarmsUsingCategory_2016 <- Monodon_2016[, .N, .(FarmID, SubstanceCategory)][,.N, SubstanceCategory]
FarmsUsingCategory <- merge(FarmsUsingCategory_2008, FarmsUsingCategory_2016, by = "SubstanceCategory", all = TRUE)
FarmsUsingCategory[is.na(FarmsUsingCategory)] <- 0
setnames(FarmsUsingCategory, "N.x", "Farms_2008")
setnames(FarmsUsingCategory, "N.y", "Farms_2016")
# Relative number per total treatments
FarmsUsingCategory[, Rel_2008 := round(Farms_2008/(sum(FarmsUsingCategory$Farms_2008)), digits = 3)]
FarmsUsingCategory[, Rel_2016 := round(Farms_2016/(sum(FarmsUsingCategory$Farms_2016)), digits = 3)]
# Include the percentage of farms
FarmsUsingCategory[, Per_2008 := round(Farms_2008/(Monodon_2008[,.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[, Per_2016 := round(Farms_2016/(Monodon_2016[,.N, FarmID][,.N]) * 100, digits = 3)]

TotalVals <- cbind(Substance = "Totals", 
                   Farms_2008 = sum(FarmsUsingCategory$Farms_2008), 
                   Farms_2016 = sum(FarmsUsingCategory$Farms_2016),
                   Rel_2008 = sum(FarmsUsingCategory$Rel_2008),
                   Rel_2016 = sum(FarmsUsingCategory$Rel_2016),
                   Per_2008 = sum(FarmsUsingCategory$Per_2008),
                   Per_2016 = sum(FarmsUsingCategory$Per_2016))
FarmsUsingCategory <- rbind(FarmsUsingCategory, TotalVals, use.names = FALSE)
FarmsUsingCategoryTable <- kable(FarmsUsingCategory, 
                                 col.names = c("Substance Category", "Farms 2008", "Farms 2016",
                                               "Relative 2008", "Relative 2016",
                                               "Percentage 2008", "Percentage 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS & SUBSTANCE CATEGORY

Here the numbers are the numbers of farms reporting using the category of treatment, the relative number is the relative number of farms using that treatment compared to the total number of treatments used (by farm), and percentage is the number of farms using a category as a percentage of the total number of farms.

`r FarmsUsingCategoryTable`


For paper:  
  
  Farms in 2008 using WQ or Sediment treatment = `r  Monodon_2008[SubstanceCategory == "Sediment treatment" | SubstanceCategory == "Water quality treatment", .N, .(FarmID)][,.N]`

Farms in 2016 using WQ or Sediment treatment = `r  Monodon_2016[SubstanceCategory == "Sediment treatment" | SubstanceCategory == "Water quality treatment", .N, .(FarmID)][,.N]`




```{r Treatments per Farm By Farm Size, include=FALSE}
#################################################################
# Produce a table of the treatments per farm for P. monodon
#################################################################
library(knitr)
library(kableExtra)

FarmsUsingCategory_2008 <- Monodon_2008[, .N, .(FarmID, FarmScaleID, SubstanceCategory)][,.N, .(FarmScaleID, SubstanceCategory)]
FarmsUsingCategory_2016 <- Monodon_2016[, .N, .(FarmID, FarmScaleID, SubstanceCategory)][,.N, .(FarmScaleID, SubstanceCategory)]
setorder(FarmsUsingCategory_2008, SubstanceCategory)
setorder(FarmsUsingCategory_2016, SubstanceCategory)
FarmsUsingCategory <- merge(FarmsUsingCategory_2008, FarmsUsingCategory_2016, by = c("SubstanceCategory", "FarmScaleID"), all = TRUE)

FarmsUsingCategory[is.na(FarmsUsingCategory)] <- 0
setnames(FarmsUsingCategory, "N.x", "Farms_2008")
setnames(FarmsUsingCategory, "N.y", "Farms_2016")

# Include the percentage of farms
FarmsUsingCategory[FarmScaleID == "Large", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Large",.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[FarmScaleID == "Medium", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Medium",.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[FarmScaleID == "Small", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Small",.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[FarmScaleID == "Large", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Large",.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[FarmScaleID == "Medium", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Medium",.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingCategory[FarmScaleID == "Small", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Small",.N, FarmID][,.N]) * 100, digits = 3)]
TotalVals <- cbind(ActiveSubstance = "Totals",
                   FarmScaleID = "All",
                   Farms_2008 = sum(FarmsUsingCategory$Farms_2008), 
                   Farms_2016 = sum(FarmsUsingCategory$Farms_2016),
                   Per_2008 = sum(FarmsUsingCategory$Per_2008),
                   Per_2016 = sum(FarmsUsingCategory$Per_2016))

FarmsUsingCategory <- rbind(FarmsUsingCategory, TotalVals, use.names = FALSE)
FarmsUsingCategoryTableBYFARMSIZE <- kable(FarmsUsingCategory, 
                                           col.names = c("Substance Category", "Farm Scale", "2008", "2016", "Percent 2008", "Percent 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS & SUBSTANCE CATEGORY - By Farm Size

Here the numbers are the numbers of farms reporting using the category of treatment by farm size, the relative number is the relative number of farms using that treatment compared to the total number of treatments used (by farm), and percentage is the number of farms using a category as a percentage of the total number of farms.

`r FarmsUsingCategoryTableBYFARMSIZE`

```{r Active Substances Categories used, include=FALSE}
#############################################################
# Table of active substances
#############################################################
library(knitr)
library(kableExtra)

SubstancesFile <- "ShrimpDataSetSubstancesV5.csv"
SubstancesLookUp <- read.csv(file.path(DataFolder, SubstancesFile), header = TRUE, stringsAsFactors = FALSE)
SubstancesLookUp <- setDT(SubstancesLookUp) # this is the look up file used
setorder(SubstancesLookUp, Active.substance.category)
SubstancesLookUpForTable <- unique(SubstancesLookUp[, .N, .(Active.substance.category, Active.substance)], by = "Active.substance")
SubstancesLookUpForTable <- SubstancesLookUpForTable[, .(Active.substance.category, Active.substance)]

ActiveSubsCats_2008 <- Monodon_2008[,.N,.(FarmID, ActiveSubCategory)][,.N, ActiveSubCategory]
ActiveSubsCats_2016 <- Monodon_2016[,.N,.(FarmID, ActiveSubCategory)][,.N, ActiveSubCategory]
ActiveSubsCats <- merge(ActiveSubsCats_2008, ActiveSubsCats_2016, by = "ActiveSubCategory", all = TRUE)
ActiveSubsCats[is.na(ActiveSubsCats)] <- 0
setnames(ActiveSubsCats, "N.x", "Farms_2008")
setnames(ActiveSubsCats, "N.y", "Farms_2016")
# Relative number per total treatments
ActiveSubsCats[, Rel_2008 := round(Farms_2008/(sum(ActiveSubsCats$Farms_2008)), digits = 3)]
ActiveSubsCats[, Rel_2016 := round(Farms_2016/(sum(ActiveSubsCats$Farms_2016)), digits = 3)]
# Include the percentage of farms
ActiveSubsCats[, Per_2008 := round(Farms_2008/(Monodon_2008[,.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[, Per_2016 := round(Farms_2016/(Monodon_2016[,.N, FarmID][,.N]) * 100, digits = 3)]
TotalVals <- cbind(Substance = "Totals", 
                   Farms_2008 = sum(ActiveSubsCats$Farms_2008), 
                   Farms_2016 = sum(ActiveSubsCats$Farms_2016),
                   Rel_2008 = sum(ActiveSubsCats$Rel_2008),
                   Rel_2016 = sum(ActiveSubsCats$Rel_2016),
                   Per_2008 = sum(ActiveSubsCats$Per_2008),
                   Per_2016 = sum(ActiveSubsCats$Per_2016))
ActiveSubsCats <- rbind(ActiveSubsCats, TotalVals, use.names = FALSE)

ListActiveSubstancesCategoryTable <- kable(SubstancesLookUpForTable, 
                                           col.names = c("Active Agent Category", "Active Agent")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)

ActiveSubstancesCategoryTable <- kable(ActiveSubsCats, 
                                       col.names = c("Active Agent", "Farms 2008", "Farms 2016",
                                                     "Relative 2008", "Relative 2016",
                                                     "Percentage 2008", "Percentage 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
#############################################################################
# Basic stats
#############################################################################
Mean_NumberTreatmentsACTCATSFarm_2008 <- mean(Monodon_2008[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsACTCATSFarm_2016 <- mean(Monodon_2016[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2008 <- median(Monodon_2008[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2016 <- median(Monodon_2016[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2008 <- sum(Monodon_2008[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2016 <- sum(Monodon_2016[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2008 <- max(Monodon_2008[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2016 <- max(Monodon_2016[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2008 <- min(Monodon_2008[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2016 <- min(Monodon_2016[,.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)


```

### Active Substances By Category used by FARMS

The active substance in each treatment was noted and placed into a set of categories. The categories were:
  
  `r ListActiveSubstancesCategoryTable`

There were `r `

In 2008 there were `r Total_NumberTreatmentsACTCATSFarm_2008` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2008` and a median of `r Median_NumberTreatmentsACTCATSFarm_2008` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2008` to `r Min_NumberTreatmentsACTCATSFarm_2008` recorded active substance categories per farm.

In 2016 there were `r Total_NumberTreatmentsACTCATSFarm_2016` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2016` and a median of `r Median_NumberTreatmentsACTCATSFarm_2016` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2016` to `r Min_NumberTreatmentsACTCATSFarm_2016` recorded active substance categories per farm.


Here is a table of the number of farms using specific active ingredients to treat disease.

`r ActiveSubstancesCategoryTable`


```{r Active Substances Categories used BY FARM SCALE, include=FALSE}
#############################################################
# Table of active substances
#############################################################
library(knitr)
library(kableExtra)

ActiveSubsCats_2008 <- Monodon_2008[,.N,.(FarmID, FarmScaleID, ActiveSubCategory)][,.N, .(ActiveSubCategory, FarmScaleID)]
ActiveSubsCats_2016 <- Monodon_2016[,.N,.(FarmID, FarmScaleID, ActiveSubCategory)][,.N, .(ActiveSubCategory, FarmScaleID)]
setorder(ActiveSubsCats_2008, ActiveSubCategory)
setorder(ActiveSubsCats_2016, ActiveSubCategory)

ActiveSubsCats <- merge(ActiveSubsCats_2008, ActiveSubsCats_2016, by = c("ActiveSubCategory", "FarmScaleID"), all = TRUE)
ActiveSubsCats[is.na(ActiveSubsCats)] <- 0
setnames(ActiveSubsCats, "N.x", "Farms_2008")
setnames(ActiveSubsCats, "N.y", "Farms_2016")

# Include the percentage of farms
ActiveSubsCats[FarmScaleID == "Large", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Large",.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[FarmScaleID == "Medium", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Medium",.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[FarmScaleID == "Small", Per_2008 := round(Farms_2008/(Monodon_2008[FarmScaleID == "Small",.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[FarmScaleID == "Large", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Large",.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[FarmScaleID == "Medium", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Medium",.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubsCats[FarmScaleID == "Small", Per_2016 := round(Farms_2016/(Monodon_2016[FarmScaleID == "Small",.N, FarmID][,.N]) * 100, digits = 3)]
TotalVals <- cbind(ActiveSubstance = "Totals",
                   FarmScaleID = "All",
                   Farms_2008 = sum(ActiveSubsCats$Farms_2008), 
                   Farms_2016 = sum(ActiveSubsCats$Farms_2016),
                   Per_2008 = sum(ActiveSubsCats$Per_2008),
                   Per_2016 = sum(ActiveSubsCats$Per_2016))

ActiveSubsCats <- rbind(ActiveSubsCats, TotalVals, use.names = FALSE)

ActiveSubstancesCategoryTableBYSCALE <- kable(ActiveSubsCats, 
                                              col.names = c("Active Agent", "Farm Size", "Farms 2008", "Farms 2016",
                                                            "Percentage 2008", "Percentage 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
#############################################################################
# Basic stats
#############################################################################
Mean_NumberTreatmentsACTCATSFarm_2008_SMALL <- mean(Monodon_2008[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2008_SMALL <- median(Monodon_2008[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2008_SMALL <- sum(Monodon_2008[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2008_SMALL <- max(Monodon_2008[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2008_SMALL <- min(Monodon_2008[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsACTCATSFarm_2008_MEDIUM <- mean(Monodon_2008[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2008_MEDIUM <- median(Monodon_2008[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2008_MEDIUM <- sum(Monodon_2008[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2008_MEDIUM <- max(Monodon_2008[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2008_MEDIUM <- min(Monodon_2008[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsACTCATSFarm_2008_LARGE <- mean(Monodon_2008[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2008_LARGE <- median(Monodon_2008[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2008_LARGE <- sum(Monodon_2008[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2008_LARGE <- max(Monodon_2008[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2008_LARGE <- min(Monodon_2008[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)

Mean_NumberTreatmentsACTCATSFarm_2016_SMALL <- mean(Monodon_2016[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2016_SMALL <- median(Monodon_2016[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2016_SMALL <- sum(Monodon_2016[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2016_SMALL <- max(Monodon_2016[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2016_SMALL <- min(Monodon_2016[FarmScaleID == "Small",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsACTCATSFarm_2016_MEDIUM <- mean(Monodon_2016[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2016_MEDIUM <- median(Monodon_2016[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2016_MEDIUM <- sum(Monodon_2016[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2016_MEDIUM <- max(Monodon_2016[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2016_MEDIUM <- min(Monodon_2016[FarmScaleID == "Medium",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Mean_NumberTreatmentsACTCATSFarm_2016_LARGE <- mean(Monodon_2016[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Median_NumberTreatmentsACTCATSFarm_2016_LARGE <- median(Monodon_2016[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Total_NumberTreatmentsACTCATSFarm_2016_LARGE <- sum(Monodon_2016[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Max_NumberTreatmentsACTCATSFarm_2016_LARGE <- max(Monodon_2016[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)
Min_NumberTreatmentsACTCATSFarm_2016_LARGE <- min(Monodon_2016[FarmScaleID == "Large",.N,.(ActiveSubCategory, FarmID)][,.N, .(ActiveSubCategory,FarmID)][,.N, FarmID]$N)

```

### Active Substances By Category used by FARMS according to size of farm


In 2008:
  
  **Small farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2008_SMALL` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2008_SMALL` and a median of `r Median_NumberTreatmentsACTCATSFarm_2008_SMALL` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2008_SMALL` to `r Min_NumberTreatmentsACTCATSFarm_2008_SMALL` recorded active substance categories per farm.

**MEDIUM farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2008_MEDIUM` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2008_MEDIUM` and a median of `r Median_NumberTreatmentsACTCATSFarm_2008_MEDIUM` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2008_MEDIUM` to `r Min_NumberTreatmentsACTCATSFarm_2008_MEDIUM` recorded active substance categories per farm.

**LARGE farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2008_LARGE` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2008_LARGE` and a median of `r Median_NumberTreatmentsACTCATSFarm_2008_LARGE` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2008_LARGE` to `r Min_NumberTreatmentsACTCATSFarm_2008_LARGE` recorded active substance categories per farm.

In 2016:
  
  **Small farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2016_SMALL` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2016_SMALL` and a median of `r Median_NumberTreatmentsACTCATSFarm_2016_SMALL` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2016_SMALL` to `r Min_NumberTreatmentsACTCATSFarm_2016_SMALL` recorded active substance categories per farm.

**MEDIUM farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2016_MEDIUM` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2016_MEDIUM` and a median of `r Median_NumberTreatmentsACTCATSFarm_2016_MEDIUM` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2016_MEDIUM` to `r Min_NumberTreatmentsACTCATSFarm_2016_MEDIUM` recorded active substance categories per farm.

**LARGE farms** recorded `r Total_NumberTreatmentsACTCATSFarm_2016_LARGE` recorded active substance categories which is a mean of `r Mean_NumberTreatmentsACTCATSFarm_2016_LARGE` and a median of `r Median_NumberTreatmentsACTCATSFarm_2016_LARGE` which ranged from `r Max_NumberTreatmentsACTCATSFarm_2016_LARGE` to `r Min_NumberTreatmentsACTCATSFarm_2016_LARGE` recorded active substance categories per farm.



Here is a table of the number of farms using specific active ingredients to treat disease.

`r ActiveSubstancesCategoryTableBYSCALE`



```{r Active Substances used, include=FALSE}
#############################################################
# Table of active substances
#############################################################
library(knitr)
library(kableExtra)
ActiveSubs_2008 <- Monodon_2008[,.N,.(FarmID, ActiveSubstance)][,.N, ActiveSubstance]
ActiveSubs_2016 <- Monodon_2016[,.N,.(FarmID, ActiveSubstance)][,.N, ActiveSubstance]
ActiveSubs <- merge(ActiveSubs_2008, ActiveSubs_2016, by = "ActiveSubstance", all = TRUE)
ActiveSubs[is.na(ActiveSubs)] <- 0
setnames(ActiveSubs, "N.x", "Farms_2008")
setnames(ActiveSubs, "N.y", "Farms_2016")
# Relative number per total treatments
ActiveSubs[, Rel_2008 := round(Farms_2008/(sum(ActiveSubs$Farms_2008)), digits = 3)]
ActiveSubs[, Rel_2016 := round(Farms_2016/(sum(ActiveSubs$Farms_2016)), digits = 3)]
# Include the percentage of farms
ActiveSubs[, Per_2008 := round(Farms_2008/(Monodon_2008[,.N, FarmID][,.N]) * 100, digits = 3)]
ActiveSubs[, Per_2016 := round(Farms_2016/(Monodon_2016[,.N, FarmID][,.N]) * 100, digits = 3)]
TotalVals <- cbind(Substance = "Totals", 
                   Farms_2008 = sum(ActiveSubs$Farms_2008), 
                   Farms_2016 = sum(ActiveSubs$Farms_2016),
                   Rel_2008 = sum(ActiveSubs$Rel_2008),
                   Rel_2016 = sum(ActiveSubs$Rel_2016),
                   Per_2008 = sum(ActiveSubs$Per_2008),
                   Per_2016 = sum(ActiveSubs$Per_2016))
ActiveSubs <- rbind(ActiveSubs, TotalVals, use.names = FALSE)

# Add in the categories
ActiveSubs <- merge(ActiveSubs, SubstancesLookUpForTable, by.x = "ActiveSubstance", by.y = "Active.substance")
setorderv(ActiveSubs, "Active.substance.category")
ActiveSubstancesTable <- kable(ActiveSubs, 
                               col.names = c("Active Agent", "Farms 2008", "Farms 2016",
                                             "Relative 2008", "Relative 2016",
                                             "Percentage 2008", "Percentage 2016", "Category")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
```

### Active Substances used by FARMS

Here is a table of the number of farms using specific active ingredients to treat disease.

`r ActiveSubstancesTable`


```{r Commerical Products per Farm, include=FALSE}
#################################################################
# Produce a table of the commercial products used per farm for P. monodon
#################################################################
library(knitr)
library(kableExtra)
FarmsCommercial_2008 <- Monodon_2008[, .N, .(FarmID, Commercial)][,.N, Commercial]
FarmsCommercial_2016 <- Monodon_2016[, .N, .(FarmID, Commercial)][,.N, Commercial]
FarmsCommercial <- merge(FarmsCommercial_2008, FarmsCommercial_2016, by = "Commercial", all = TRUE)
FarmsCommercial[is.na(FarmsCommercial)] <- 0
setnames(FarmsCommercial, "N.x", "Farms_2008")
setnames(FarmsCommercial, "N.y", "Farms_2016")
# Relative number per total treatments
FarmsCommercial[, Rel_2008 := round(Farms_2008/(sum(FarmsCommercial$Farms_2008)), digits = 3)]
FarmsCommercial[, Rel_2016 := round(Farms_2016/(sum(FarmsCommercial$Farms_2016)), digits = 3)]
# Include the percentage of farms
FarmsCommercial[, Per_2008 := round(Farms_2008/(Monodon_2008[,.N, FarmID][,.N]) * 100, digits = 3)]
FarmsCommercial[, Per_2016 := round(Farms_2016/(Monodon_2016[,.N, FarmID][,.N]) * 100, digits = 3)]

TotalVals <- cbind(Substance = "Totals", 
                   Farms_2008 = sum(FarmsCommercial$Farms_2008), 
                   Farms_2016 = sum(FarmsCommercial$Farms_2016),
                   Rel_2008 = sum(FarmsCommercial$Rel_2008),
                   Rel_2016 = sum(FarmsCommercial$Rel_2016),
                   Per_2008 = sum(FarmsCommercial$Per_2008),
                   Per_2016 = sum(FarmsCommercial$Per_2016))
FarmsCommercial <- rbind(FarmsCommercial, TotalVals, use.names = FALSE)
FarmsCommercialTable <- kable(FarmsCommercial, 
                              col.names = c("Commerical Product", "Farms 2008", "Farms 2016",
                                            "Relative 2008", "Relative 2016",
                                            "Percentage 2008", "Percentage 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS & COMMERCIAL PRODUCT

Here the number of farms using a commercial product versus non-commercial product are reported.

`r FarmsCommercialTable`

```{r Sunburst SubstancesUsed, include=FALSE}
#############################################################################
#
# Here a sunburst of the substances used by a farm
#   FIRST by COMMERCIAL PRODUCT
#   SECOND by SUBSTANCE CATEGORY
############################################################################
library(plotly)
# Prepare the data for all agents - NEED TO BE 1 PER FARM
CommercialDataFarms_2008 <- unique(Monodon_2008[, .(Commercial, SubstanceCategory), FarmID], by = c("FarmID","Commercial","SubstanceCategory"))
CommercialDataFarms_2016 <- unique(Monodon_2016[, .(Commercial, SubstanceCategory), FarmID], by = c("FarmID","Commercial","SubstanceCategory"))
setorder(CommercialDataFarms_2008, Commercial)
setorder(CommercialDataFarms_2016, Commercial)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "NotCommercial", "Commercial",
                 paste0("No_", CommercialDataFarms_2008[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory),
                 paste0("Yes_", CommercialDataFarms_2008[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2008)<br>Farms",
               "Non-commercial Product", "Commercial Product",
               CommercialDataFarms_2008[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory,
               CommercialDataFarms_2008[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory
)

# Produce a parent list
ParentList <-c("",
               "Year","Year",
               rep("NotCommercial", length(CommercialDataFarms_2008[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory)),
               rep("Commercial", length(CommercialDataFarms_2008[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory))
)

# Ensure segmentation is correct
Values <- c(CommercialDataFarms_2008[,.N], # Parent
            CommercialDataFarms_2008[,.N, Commercial]$N,
            CommercialDataFarms_2008[Commercial == "No",.N, SubstanceCategory]$N,
            CommercialDataFarms_2008[Commercial == "Yes",.N, SubstanceCategory]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#74c476","#99d8c9", # 
                    rep("#c7e9c0", length(CommercialDataFarms_2008[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory)),
                    rep("#ccece6", length(CommercialDataFarms_2008[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory))
)

# Produce the dataframe
Prod_Dataframe_2008_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)


# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "NotCommercial", "Commercial",
                 paste0("No_", CommercialDataFarms_2016[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory),
                 paste0("Yes_", CommercialDataFarms_2016[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2016)<br>Farms",
               "Non-commercial Product", "Commercial Product",
               CommercialDataFarms_2016[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory,
               CommercialDataFarms_2016[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory
)

# Produce a parent list
ParentList <-c("",
               "Year","Year",
               rep("NotCommercial", length(CommercialDataFarms_2016[Commercial == "No",.N, SubstanceCategory]$SubstanceCategory)),
               rep("Commercial", length(CommercialDataFarms_2016[Commercial == "Yes",.N, SubstanceCategory]$SubstanceCategory))
)

# Ensure segmentation is correct
Values <- c(CommercialDataFarms_2016[,.N], # Parent
            CommercialDataFarms_2016[,.N, Commercial]$N,
            CommercialDataFarms_2016[Commercial == "No",.N, SubstanceCategory]$N,
            CommercialDataFarms_2016[Commercial == "Yes",.N, SubstanceCategory]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#74c476","#99d8c9", # 
                    rep("#c7e9c0", length(CommercialDataFarms_2016[Commercial == "No",.N, SubstanceCategory]$N)),
                    rep("#ccece6", length(CommercialDataFarms_2016[Commercial == "Yes",.N, SubstanceCategory]$N))
)

# Produce the dataframe
Prod_Dataframe_2016_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce the plot
SunburstPlot_Combined_Farm <- plot_ly() %>%
  add_trace(type = 'sunburst',
            ids = Prod_Dataframe_2008_Farm$SegmentList, 
            labels = Prod_Dataframe_2008_Farm$LabelList, 
            parents = Prod_Dataframe_2008_Farm$ParentList,
            values = Prod_Dataframe_2008_Farm$Values,
            marker = list(colors = Prod_Dataframe_2008_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 0))

SunburstPlot_Combined_Farm <- SunburstPlot_Combined_Farm %>%
  add_trace(type = 'sunburst',
            ids = Prod_Dataframe_2016_Farm$SegmentList, 
            labels = Prod_Dataframe_2016_Farm$LabelList, 
            parents = Prod_Dataframe_2016_Farm$ParentList,
            values = Prod_Dataframe_2016_Farm$Values,
            marker = list(colors = Prod_Dataframe_2016_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 1))

SunburstPlot_Combined_Farm <- SunburstPlot_Combined_Farm %>%
  layout(grid = list(columns = 2, rows =1),
         margin = list(l = 0, r = 0, b = 0, t = 0))
SunburstPlot_Combined_Farm

```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS By Product then Category

These sunburst charts show the use of commercial vs non-commercial products and the use of these products between 2008 and 2016.

`r SunburstPlot_Combined_Farm`

```{r Sunburst SubstancesUsedBySubs, include=FALSE}
#############################################################################
#
# Here a sunburst of the substances used by a farm
#   FIRST by COMMERCIAL PRODUCT
#   SECOND by SUBSTANCE CATEGORY
############################################################################
library(plotly)
# Prepare the data for all agents - NEED TO BE 1 PER FARM
Sub_CommercialDataFarms_2008 <- unique(Monodon_2008[, .(Commercial, SubstanceCategory), FarmID], by = c("FarmID","Commercial","SubstanceCategory"))
Sub_CommercialDataFarms_2016 <- unique(Monodon_2016[, .(Commercial, SubstanceCategory), FarmID], by = c("FarmID","Commercial","SubstanceCategory"))
setorder(Sub_CommercialDataFarms_2008, SubstanceCategory)
setorder(Sub_CommercialDataFarms_2016, SubstanceCategory)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 Sub_CommercialDataFarms_2008[,.N, SubstanceCategory]$SubstanceCategory,
                 paste0(Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)][,.(SubstanceCategory, Commercial)]$SubstanceCategory,Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)][,.(SubstanceCategory, Commercial)]$Commercial)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2008)<br>Farms",
               Sub_CommercialDataFarms_2008[,.N, SubstanceCategory]$SubstanceCategory,
               Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)]$Commercial
)

# Produce a parent list
ParentList <-c("",
               rep("Year",length(Sub_CommercialDataFarms_2008[,.N, SubstanceCategory]$SubstanceCategory)),
               Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)]$SubstanceCategory )

# Ensure segmentation is correct
Values <- c(Sub_CommercialDataFarms_2008[,.N], # Parent
            Sub_CommercialDataFarms_2008[,.N, SubstanceCategory]$N,
            Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)]$N)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    rep("#74c476", length(Sub_CommercialDataFarms_2008[,.N, SubstanceCategory]$SubstanceCategory)), # 
                    Sub_CommercialDataFarms_2008[,.N, .(SubstanceCategory, Commercial)][, Colour := ifelse(Commercial == "Yes", "#d8b365","#5ab4ac")]$Colour
)

# Produce the dataframe
Subs_Dataframe_2008_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 Sub_CommercialDataFarms_2016[,.N, SubstanceCategory]$SubstanceCategory,
                 paste0(Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)][,.(SubstanceCategory, Commercial)]$SubstanceCategory,Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)][,.(SubstanceCategory, Commercial)]$Commercial)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2016)<br>Farms",
               Sub_CommercialDataFarms_2016[,.N, SubstanceCategory]$SubstanceCategory,
               Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)]$Commercial
)

# Produce a parent list
ParentList <-c("",
               rep("Year",length(Sub_CommercialDataFarms_2016[,.N, SubstanceCategory]$SubstanceCategory)),
               Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)]$SubstanceCategory )

# Ensure segmentation is correct
Values <- c(Sub_CommercialDataFarms_2016[,.N], # Parent
            Sub_CommercialDataFarms_2016[,.N, SubstanceCategory]$N,
            Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)]$N)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    rep("#74c476", length(Sub_CommercialDataFarms_2016[,.N, SubstanceCategory]$SubstanceCategory)), # 
                    Sub_CommercialDataFarms_2016[,.N, .(SubstanceCategory, Commercial)][, Colour := ifelse(Commercial == "Yes", "#d8b365","#5ab4ac")]$Colour
)

# Produce the dataframe
Subs_Dataframe_2016_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce the plot
SunburstPlot_Combined_Farm_Subs <- plot_ly() %>%
  add_trace(type = 'sunburst',
            ids = Subs_Dataframe_2008_Farm$SegmentList, 
            labels = Subs_Dataframe_2008_Farm$LabelList, 
            parents = Subs_Dataframe_2008_Farm$ParentList,
            values = Subs_Dataframe_2008_Farm$Values,
            marker = list(colors = Subs_Dataframe_2008_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 0))

SunburstPlot_Combined_Farm_Subs <- SunburstPlot_Combined_Farm_Subs %>%
  add_trace(type = 'sunburst',
            ids = Subs_Dataframe_2016_Farm$SegmentList, 
            labels = Subs_Dataframe_2016_Farm$LabelList, 
            parents = Subs_Dataframe_2016_Farm$ParentList,
            values = Subs_Dataframe_2016_Farm$Values,
            marker = list(colors = Subs_Dataframe_2016_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 1))

SunburstPlot_Combined_Farm_Subs <- SunburstPlot_Combined_Farm_Subs %>%
  layout(grid = list(columns = 2, rows =1),
         margin = list(l = 0, r = 0, b = 0, t = 0))
SunburstPlot_Combined_Farm_Subs

```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS By Category then Product Type

These sunburst charts show the use of different treatment categories and whether these were commercial products or not between 2008 and 2016.

`r SunburstPlot_Combined_Farm_Subs`


```{r SubstancesUsed per Farm, include=FALSE}
#################################################################
# Produce a table of the treatments per farm for P. monodon
#################################################################
library(knitr)
library(kableExtra)
FarmsUsingSubstanceUsed_2008 <- Monodon_2008[, .N, .(FarmID, SubstanceUsed)][,.N, SubstanceUsed]
FarmsUsingSubstanceUsed_2016 <- Monodon_2016[, .N, .(FarmID, SubstanceUsed)][,.N, SubstanceUsed]
FarmsUsingSubstanceUsed <- merge(FarmsUsingSubstanceUsed_2008, FarmsUsingSubstanceUsed_2016, by = "SubstanceUsed", all = TRUE)
FarmsUsingSubstanceUsed[is.na(FarmsUsingSubstanceUsed)] <- 0
setnames(FarmsUsingSubstanceUsed, "N.x", "Farms_2008")
setnames(FarmsUsingSubstanceUsed, "N.y", "Farms_2016")
# Relative number per total treatments
FarmsUsingSubstanceUsed[, Rel_2008 := round(Farms_2008/(sum(FarmsUsingSubstanceUsed$Farms_2008)), digits = 3)]
FarmsUsingSubstanceUsed[, Rel_2016 := round(Farms_2016/(sum(FarmsUsingSubstanceUsed$Farms_2016)), digits = 3)]
# Include the percentage of farms
FarmsUsingSubstanceUsed[, Per_2008 := round(Farms_2008/(Monodon_2008[,.N, FarmID][,.N]) * 100, digits = 3)]
FarmsUsingSubstanceUsed[, Per_2016 := round(Farms_2016/(Monodon_2016[,.N, FarmID][,.N]) * 100, digits = 3)]

TotalVals <- cbind(Substance = "Totals", 
                   Farms_2008 = sum(FarmsUsingSubstanceUsed$Farms_2008), 
                   Farms_2016 = sum(FarmsUsingSubstanceUsed$Farms_2016),
                   Rel_2008 = sum(FarmsUsingSubstanceUsed$Rel_2008),
                   Rel_2016 = sum(FarmsUsingSubstanceUsed$Rel_2016),
                   Per_2008 = sum(FarmsUsingSubstanceUsed$Per_2008),
                   Per_2016 = sum(FarmsUsingSubstanceUsed$Per_2016))
FarmsUsingSubstanceUsed <- rbind(FarmsUsingSubstanceUsed, TotalVals, use.names = FALSE)
FarmsUsingSubstanceUsedTable <- kable(FarmsUsingSubstanceUsed, 
                                      col.names = c("Substance Used", "Farms 2008", "Farms 2016",
                                                    "Relative 2008", "Relative 2016",
                                                    "Percentage 2008", "Percentage 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE)
```

### Details of the substances used by farms in the survey - INDIVIDUAL FARMS & SUBSTANCE CATEGORY

Here the numbers are the numbers of farms reporting using the category of treatment, the relative number is the relative number of farms using that treatment compared to the total number of treatments used (by farm), and percentage is the number of farms using a category as a percentage of the total number of farms.

`r FarmsUsingSubstanceUsedTable`

** Updated to this Point


```{r Sunburst ActiveSubCatBySize, include=FALSE}
#############################################################################
#
# Here a sunburst of the substances active substance categories used
#   FIRST by FARM SIZE
#   SECOND by ACTIVE SUBSTANCE CATEGORY
############################################################################
library(plotly)
# Prepare the data for all agents - NEED TO BE 1 PER FARM
ActiveSubDataFarms_2008 <- unique(Monodon_2008[, .(FarmScaleID, ActiveSubCategory), FarmID], by = c("FarmID","FarmScaleID","ActiveSubCategory"))
ActiveSubDataFarms_2016 <- unique(Monodon_2016[, .(FarmScaleID, ActiveSubCategory), FarmID], by = c("FarmID","FarmScaleID","ActiveSubCategory"))
setorder(ActiveSubDataFarms_2008, FarmScaleID)
setorder(ActiveSubDataFarms_2016, FarmScaleID)

# Introduce Segment colours for the categories
ActiveSubDataFarms_2008[ActiveSubCategory == "Inorganic compound", SegCol := "#e5f5f9"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Inorganic minerals", SegCol := "#99d8c9"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Disinfectants - Cl based", SegCol := "#8856a7"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Disinfectants - Non Cl based", SegCol := "#9ebcda"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Pesticide", SegCol := "#fdbb84"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Antimicrobial", SegCol := "#2b8cbe"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Oxygen precursors", SegCol := "#1c9099"]
ActiveSubDataFarms_2008[ActiveSubCategory == "Coagulant", SegCol := "#fa9fb5"]

ActiveSubDataFarms_2016[ActiveSubCategory == "Inorganic compound", SegCol := "#e5f5f9"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Inorganic minerals", SegCol := "#99d8c9"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Disinfectants - Cl based", SegCol := "#8856a7"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Disinfectants - Non Cl based", SegCol := "#9ebcda"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Pesticide", SegCol := "#fdbb84"]

ActiveSubDataFarms_2016[ActiveSubCategory == "Antimicrobial", SegCol := "#2b8cbe"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Antioxidants", SegCol := "#a6bddb"]

ActiveSubDataFarms_2016[ActiveSubCategory == "Oxygen precursors", SegCol := "#1c9099"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Coagulant", SegCol := "#fa9fb5"]

ActiveSubDataFarms_2016[ActiveSubCategory == "Feed Supplement", SegCol := "#4a1486"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Probiotic - bacteria", SegCol := "#6a51a3"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Probiotic - other", SegCol := "#807dba"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Vitamins", SegCol := "#9e9ac8"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Fertiliser", SegCol := "#bcbddc"]
ActiveSubDataFarms_2016[ActiveSubCategory == "Enzymes", SegCol := "#dadaeb"]

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 ActiveSubDataFarms_2008[,.N, FarmScaleID]$FarmScaleID,
                 paste0("Large_", ActiveSubDataFarms_2008[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory),
                 paste0("Medium_", ActiveSubDataFarms_2008[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory),
                 paste0("Small_", ActiveSubDataFarms_2008[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2008)<br>Farms",
               ActiveSubDataFarms_2008[,.N, FarmScaleID]$FarmScaleID,
               ActiveSubDataFarms_2008[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory,
               ActiveSubDataFarms_2008[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory,
               ActiveSubDataFarms_2008[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory
)

# Produce a parent list
ParentList <-c("",
               rep("Year",length(ActiveSubDataFarms_2008[,.N, FarmScaleID]$FarmScaleID)),
               rep("Large", length(ActiveSubDataFarms_2008[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory)),
               rep("Medium", length(ActiveSubDataFarms_2008[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory)),
               rep("Small", length(ActiveSubDataFarms_2008[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory))
)

# Ensure segmentation is correct
Values <- c(sum(ActiveSubDataFarms_2008[,.N, FarmScaleID]$N),
            ActiveSubDataFarms_2008[,.N, FarmScaleID]$N,# Parent
            ActiveSubDataFarms_2008[FarmScaleID == "Large",.N, ActiveSubCategory]$N,
            ActiveSubDataFarms_2008[FarmScaleID == "Medium",.N, ActiveSubCategory]$N,
            ActiveSubDataFarms_2008[FarmScaleID == "Small",.N, ActiveSubCategory]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#7fc97f","#beaed4", "#fdc086", # 
                    ActiveSubDataFarms_2008[FarmScaleID == "Large",.N, .(ActiveSubCategory, SegCol)]$SegCol,
                    ActiveSubDataFarms_2008[FarmScaleID == "Medium",.N, .(ActiveSubCategory, SegCol)]$SegCol,
                    ActiveSubDataFarms_2008[FarmScaleID == "Small",.N, .(ActiveSubCategory, SegCol)]$SegCol
)

# Produce the dataframe
ActSubCat_Dataframe_2008_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 ActiveSubDataFarms_2016[,.N, FarmScaleID]$FarmScaleID,
                 paste0("Large_", ActiveSubDataFarms_2016[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory),
                 paste0("Medium_", ActiveSubDataFarms_2016[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory),
                 paste0("Small_", ActiveSubDataFarms_2016[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2016)<br>Farms",
               ActiveSubDataFarms_2016[,.N, FarmScaleID]$FarmScaleID,
               ActiveSubDataFarms_2016[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory,
               ActiveSubDataFarms_2016[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory,
               ActiveSubDataFarms_2016[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory
)

# Produce a parent list
ParentList <-c("",
               rep("Year",length(ActiveSubDataFarms_2016[,.N, FarmScaleID]$FarmScaleID)),
               rep("Large", length(ActiveSubDataFarms_2016[FarmScaleID == "Large",.N, ActiveSubCategory]$ActiveSubCategory)),
               rep("Medium", length(ActiveSubDataFarms_2016[FarmScaleID == "Medium",.N, ActiveSubCategory]$ActiveSubCategory)),
               rep("Small", length(ActiveSubDataFarms_2016[FarmScaleID == "Small",.N, ActiveSubCategory]$ActiveSubCategory))
)

# Ensure segmentation is correct
Values <- c(sum(ActiveSubDataFarms_2016[,.N, FarmScaleID]$N),
            ActiveSubDataFarms_2016[,.N, FarmScaleID]$N,# Parent
            ActiveSubDataFarms_2016[FarmScaleID == "Large",.N, ActiveSubCategory]$N,
            ActiveSubDataFarms_2016[FarmScaleID == "Medium",.N, ActiveSubCategory]$N,
            ActiveSubDataFarms_2016[FarmScaleID == "Small",.N, ActiveSubCategory]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#7fc97f","#beaed4", "#fdc086", # 
                    ActiveSubDataFarms_2016[FarmScaleID == "Large",.N, .(ActiveSubCategory, SegCol)]$SegCol,
                    ActiveSubDataFarms_2016[FarmScaleID == "Medium",.N, .(ActiveSubCategory, SegCol)]$SegCol,
                    ActiveSubDataFarms_2016[FarmScaleID == "Small",.N, .(ActiveSubCategory, SegCol)]$SegCol
)

# Produce the dataframe
ActSubCat_Dataframe_2016_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce the plot
SunburstPlot_Combined_FarmSize_ActSubs <- plot_ly() %>%
  add_trace(type = 'sunburst',
            ids = ActSubCat_Dataframe_2008_Farm$SegmentList, 
            labels = ActSubCat_Dataframe_2008_Farm$LabelList, 
            parents = ActSubCat_Dataframe_2008_Farm$ParentList,
            values = ActSubCat_Dataframe_2008_Farm$Values,
            marker = list(colors = ActSubCat_Dataframe_2008_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 0))

SunburstPlot_Combined_FarmSize_ActSubs <- SunburstPlot_Combined_FarmSize_ActSubs %>%
  add_trace(type = 'sunburst',
            ids = ActSubCat_Dataframe_2016_Farm$SegmentList, 
            labels = ActSubCat_Dataframe_2016_Farm$LabelList, 
            parents = ActSubCat_Dataframe_2016_Farm$ParentList,
            values = ActSubCat_Dataframe_2016_Farm$Values,
            marker = list(colors = ActSubCat_Dataframe_2016_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 1))

SunburstPlot_Combined_FarmSize_ActSubs <- SunburstPlot_Combined_FarmSize_ActSubs %>%
  layout(grid = list(columns = 2, rows =1),
         margin = list(l = 0, r = 0, b = 0, t = 0))
SunburstPlot_Combined_FarmSize_ActSubs
```

### Details of the active substances used by farms in the survey - INDIVIDUAL FARMS By Size then Active Substance Category

These sunburst charts show the use of different active substance categories based on farm size between 2008 and 2016.

`r SunburstPlot_Combined_FarmSize_ActSubs`