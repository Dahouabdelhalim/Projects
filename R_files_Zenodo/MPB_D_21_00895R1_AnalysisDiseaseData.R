---
  title: "Analysis of Disease in Penaeus monodon Shrimp Farms"
author: "Richard Heal"
date: "June 12, 2020"
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
DataFile_2008 <- "MPB_D_21_00895R1_ShrimpDataSet_Disease_2008.csv"
DataFile_2016 <- "MPB_D_21_00895R1_ShrimpDataSet_Disease_2016.csv"

library(data.table)
library(kableExtra)
RawData_2008 <- read.csv(file.path(DataFolder, DataFile_2008),
                         header = TRUE,
                         stringsAsFactors = FALSE)
RawData_2008 <- setDT(RawData_2008)
RawData_2008 <- RawData_2008[!is.na(FarmID), ]

RawData_2016 <- read.csv(file.path(DataFolder, DataFile_2016),
                         header = TRUE,
                         stringsAsFactors = FALSE)
RawData_2016 <- setDT(RawData_2016)
```

## Analysis of Disease Data for Penaeus monodon

This document shows the analysis of the cleaned disease data set.

## Preliminary Analysis of Shrimp Disease data

In this section simple data and plots of the data are obtained to understand the data.

## TAKE NOTE 

** The data is not per disease case and therefore in this version ALL DATA IS PRESENTED PER FARM **
  ** This data is only for records for P. monodon due to no records for M.rosengergii present in 2016 data **
  
  ```{r Farms Information, include=FALSE}
#####################################################################
#
# Generate a table of the data - 
#
#####################################################################

########################################################
#
# Strip out just monodon records and create single datatable PER FARM
#
#########################################################
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

### Details of the farms in the survey - INDIVIDUAL FARMS

The details of the farms in the dataset (for P. monodon) are:
  
  `r DescriptiveDataTable`


```{r Disease by Type, include=FALSE}
#########################################################
#
# Output counts of disease type (infectious or other) - for monodon
#
########################################################
setorder(RawData_2008, Cause)
setorder(RawData_2016, Cause)
DiseaseByCause <- data.frame(
  Cause = c("Infectious", "Physical Damage","Physical Damage or Nutritional Deficiency","Nutrition Deficiency", "Unknown"),
  Yr_2008 = c(RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Infectious", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Physical", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Physical/Nutritional", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Nutrition", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Unknown", .N, FarmID][, .N]),
  Yr_2016 = c(RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Infectious", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Physical", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Physical/Nutritional", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Nutrition", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & (Cause == "Unknown" ), .N, FarmID][, .N])
)
DiseaseByCause <- setDT(DiseaseByCause)
DiseaseByCause[, Rel_2008 := round(Yr_2008/ DiseaseByCause[, sum(Yr_2008)], digits = 3)]
DiseaseByCause[, Rel_2016 := round(Yr_2016/ DiseaseByCause[, sum(Yr_2016)], digits = 3)]
Totals <- cbind("Total",  as.numeric(DiseaseByCause[, sum(Yr_2008)]), as.numeric(DiseaseByCause[, sum(Yr_2016)]), 
                as.numeric(DiseaseByCause[, sum(Rel_2008)]), as.numeric(DiseaseByCause[, sum(Rel_2016)]))
DiseaseByCause <- rbind(
  DiseaseByCause,
  Totals, use.names = FALSE)
DiseaseCauseTable <- kable(DiseaseByCause, col.names = c("Cause", "2008", "2016", "Relative 2008", "Relative 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)

```


### Details of the disease signs causes by agent - CAUSE PER FARM

Each disease sign was labelled as either causes by an infectious agent, physical damage, or an unknown cause. The details of the survey are. 

`r DiseaseCauseTable`



```{r Disease by Agent, include=FALSE}
#########################################################
#
# Output counts of disease causes - for monodon
#
########################################################
setorder(RawData_2008, Agent)
setorder(RawData_2016, Agent)
DiseaseByAgent <- data.frame(
  Cause = c("Bacteria", "Bacteria/Fungus","Parasite", "Virus", "Other"),
  Yr_2008 = c(RawData_2008[SpeciesID == "Penaeus monodon" & Agent == "Bacteria", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Agent == "Bacteria/Fungus", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Agent == "Parasite", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Agent == "Virus", .N, FarmID][, .N],
              RawData_2008[SpeciesID == "Penaeus monodon" & Agent == "None", .N, FarmID][, .N]),
  Yr_2016 = c(RawData_2016[SpeciesID == "Penaeus monodon" & Agent == "Bacteria", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Agent == "Bacteria/Fungus", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Agent == "Parasite", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & Agent == "Virus", .N, FarmID][, .N],
              RawData_2016[SpeciesID == "Penaeus monodon" & (Agent == "None" |Agent == "Unknown" ), .N, FarmID][, .N])
)
DiseaseByAgent <- setDT(DiseaseByAgent)
DiseaseByAgent[, Rel_2008 := round(Yr_2008/ DiseaseByAgent[, sum(Yr_2008)], digits = 3)]
DiseaseByAgent[, Rel_2016 := round(Yr_2016/ DiseaseByAgent[, sum(Yr_2016)], digits = 3)]
Totals <- cbind("Total",  as.numeric(DiseaseByAgent[, sum(Yr_2008)]), as.numeric(DiseaseByAgent[, sum(Yr_2016)]), 
                as.numeric(DiseaseByAgent[, sum(Rel_2008)]), as.numeric(DiseaseByAgent[, sum(Rel_2016)]))
DiseaseByAgent <- rbind(
  DiseaseByAgent,
  Totals, use.names = FALSE)
DiseaseAgentTable <- kable(DiseaseByAgent, col.names = c("Agent", "2008", "2016", "Relative 2008", "Relative 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)

```


### Details of the disease signs causes by agent - PER FARM

Each disease sign was labelled as either caused by a bacteria, bacteria and/or fungus, parasite, virus or other. The details of the survey are. 

`r DiseaseAgentTable`



```{r Disease by Pathogen, include=FALSE}
#########################################################
#
# Output counts of disease causes - for monodon
#
########################################################
setorder(RawData_2008, Disease)
setorder(RawData_2016, Disease)

# Sort diseases by year and per farm
Monodon_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon", .(Disease, FarmScaleID), FarmID], by = c("FarmID", "Disease")) # Diseases unique to each farm
Monodon_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon", .(Disease, FarmScaleID), FarmID], by = c("FarmID", "Disease")) # Diseases unique to each farm
Monodon_2008[, SurveyYear := 2008]
Monodon_2016[, SurveyYear := 2016]

# Rearrange records for counts of each disease overall and for each category
M_2008_All <- Monodon_2008[, .N, .(Disease)]
M_2008_Large <- Monodon_2008[FarmScaleID == "Large", .N, .(Disease)]
M_2008_Medium <- Monodon_2008[FarmScaleID == "Medium", .N, .(Disease)]
M_2008_Small <- Monodon_2008[FarmScaleID == "Small", .N, .(Disease)]
M_2008_Disease <- merge(M_2008_All,M_2008_Large, by = "Disease", all = TRUE)
M_2008_Disease <- merge(M_2008_Disease,M_2008_Medium, by = "Disease", all = TRUE)
M_2008_Disease <- merge(M_2008_Disease,M_2008_Small, by = "Disease", all = TRUE)
colnames(M_2008_Disease) <- c("Disease", "Total_2008", "Large_2008", "Medium_2008", "Small_2008")

# Rearrange records for counts of each disease overall and for each category
M_2016_All <- Monodon_2016[, .N, .(Disease)]
M_2016_Large <- Monodon_2016[FarmScaleID == "Large", .N, .(Disease)]
M_2016_Medium <- Monodon_2016[FarmScaleID == "Medium", .N, .(Disease)]
M_2016_Small <- Monodon_2016[FarmScaleID == "Small", .N, .(Disease)]
M_2016_Disease <- merge(M_2016_All,M_2016_Large, by = "Disease", all = TRUE)
M_2016_Disease <- merge(M_2016_Disease,M_2016_Medium, by = "Disease", all = TRUE)
M_2016_Disease <- merge(M_2016_Disease,M_2016_Small, by = "Disease", all = TRUE)
colnames(M_2016_Disease) <- c("Disease", "Total_2016", "Large_2016", "Medium_2016", "Small_2016")

MonodonDiseaseCounts <- merge(M_2008_Disease, M_2016_Disease, by = "Disease", all = TRUE)
MonodonDiseaseCounts <- MonodonDiseaseCounts[, c(1,2,6,5,9,4,8,3,7)]
MonodonDiseaseCounts[is.na(MonodonDiseaseCounts)] <- 0
setorder(MonodonDiseaseCounts, "Disease")

# Get the totals
TotalVals <- MonodonDiseaseCounts[, .(sum(Total_2008), sum(Total_2016),
                                      sum(Small_2008), sum(Small_2016),
                                      sum(Medium_2008), sum(Medium_2016),
                                      sum(Large_2008), sum(Large_2016))]
MonodonDiseaseCounts <- rbind(MonodonDiseaseCounts,
                              c("Totals", TotalVals[1,]), use.names = FALSE)

SpecificDiseaseTable <- kable(MonodonDiseaseCounts, col.names = c("Disease", "Total 2008", "Total 2016", "2008 Small Farm", "2016 Small Farm", "2008 Medium Farm", "2016 Medium Farm", "2008 Large Farm", "2016 Large Farm")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)
```

### Details of the disease signs causes by specific disease agent - PER FARM

Each disease signs were labelled according to disease signs and placed into categories of disease. 

`r SpecificDiseaseTable`


```{r Relative Disease by Agent, include=FALSE}
#########################################################
#
# Output relative counts of disease causes - for monodon
#
########################################################
Rel_MonodonDiseaseCounts <- copy(MonodonDiseaseCounts)
Rel_MonodonDiseaseCounts[, Total_2008 := round(Total_2008/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Total_2008, digits = 3)]
Rel_MonodonDiseaseCounts[, Total_2016 := round(Total_2016/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Total_2016, digits = 3)]
Rel_MonodonDiseaseCounts[, Small_2008 := round(Small_2008/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Small_2008, digits = 3)]
Rel_MonodonDiseaseCounts[, Small_2016 := round(Small_2016/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Small_2016, digits = 3)]
Rel_MonodonDiseaseCounts[, Medium_2008 := round(Medium_2008/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Medium_2008, digits = 3)]
Rel_MonodonDiseaseCounts[, Medium_2016 := round(Medium_2016/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Medium_2016, digits = 3)]
Rel_MonodonDiseaseCounts[, Large_2008 := round(Large_2008/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Large_2008, digits = 3)]
Rel_MonodonDiseaseCounts[, Large_2016 := round(Large_2016/Rel_MonodonDiseaseCounts[Disease == "Totals", ]$Large_2016, digits = 3)]

REL_SpecificDiseaseTable <- kable(Rel_MonodonDiseaseCounts, col.names = c("Disease", "Total 2008", "Total 2016", "2008 Small Farm", "2016 Small Farm", "2008 Medium Farm", "2016 Medium Farm", "2008 Large Farm", "2016 Large Farm")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)
```

### Details of the disease signs causes by specific disease agent - RELATIVE Numbers PER FARM

Each disease signs were labelled according to disease signs and placed into categories of disease. This is the relative numbers.

`r REL_SpecificDiseaseTable`



```{r Diseases Per Farm, include = FALSE}
###################################################################################
#  Output the number of diseases per farm - number of different recorded diseases
#
# Note: this is for FARMS reporting a case of the disease
#
# 1 diseases, more than 1 disease, and then singly up to maximum of 9
####################################################################################
NoDiseasesPerFarm_2008 <- RawData_2008[SpeciesID == "Penaeus monodon", .N, .(FarmID, Disease)][,.N, FarmID]
NoDiseasesPerFarm_2016 <- RawData_2016[SpeciesID == "Penaeus monodon", .N, .(FarmID, Disease)][,.N, FarmID]

NoDiseasesPerFarm <- data.frame(
  Single_2008 = NoDiseasesPerFarm_2008[N == 1, .N],
  Single_2016 = NoDiseasesPerFarm_2016[N == 1, .N],
  Multiple_2008 = NoDiseasesPerFarm_2008[N > 1, .N],
  Multiple_2016 = NoDiseasesPerFarm_2016[N > 1, .N],
  Two_2008 = NoDiseasesPerFarm_2008[N == 2, .N],
  Two_2016 = NoDiseasesPerFarm_2016[N == 2, .N],
  Three_2008 = NoDiseasesPerFarm_2008[N == 3, .N],
  Three_2016 = NoDiseasesPerFarm_2016[N == 3, .N],
  Four_2008 = NoDiseasesPerFarm_2008[N == 4, .N],
  Four_2016 = NoDiseasesPerFarm_2016[N == 4, .N],
  Five_2008 = NoDiseasesPerFarm_2008[N == 5, .N],
  Five_2016 = NoDiseasesPerFarm_2016[N == 5, .N],
  Six_2008 = NoDiseasesPerFarm_2008[N == 6, .N],
  Six_2016 = NoDiseasesPerFarm_2016[N == 6, .N],
  Seven_2008 = NoDiseasesPerFarm_2008[N == 7, .N],
  Seven_2016 = NoDiseasesPerFarm_2016[N == 7, .N],
  Eight_2008 = NoDiseasesPerFarm_2008[N == 8, .N],
  Eight_2016 = NoDiseasesPerFarm_2016[N == 8, .N],
  Nine_2008 = NoDiseasesPerFarm_2008[N == 9, .N],
  Nine_2016 = NoDiseasesPerFarm_2016[N == 9, .N]
)

DiseaseNumbersTable <- kable(NoDiseasesPerFarm, col.names = c("Single Incident 2008", "Single Incident 2016", "Multiple Incident 2008","Multiple Incident 2016", "2 2008", "2 2016", "3 2008", "3 2016", "4 2008", "4 2016", "5 2008", "5 2016", "6 2008", "6 2016", "7 2008", "7 2016", "8 2008", "8 2016", "9 2008", "9 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)

```

### Details of the number of diseases reported by each farm - Numbers PER FARM

For each farm the number of different identifiable diseases reported on that farm was determined.

`r DiseaseNumbersTable`


```{r Diseases Per Number and Farm Type, include = FALSE}

#########################################################
#
# Output the number of times a disease is recorded on a farm (single vs multiple)
#
########################################################
# Sort diseases by year and per farm
Monodon_2008_BYSIZE <- RawData_2008[SpeciesID == "Penaeus monodon", .(Disease, FarmScaleID), FarmID] # Diseases per farm
Monodon_2016_BYSIZE <- RawData_2016[SpeciesID == "Penaeus monodon", .(Disease, FarmScaleID), FarmID] # Diseases per farm
Monodon_2008_BYSIZE[, SurveyYear := 2008]
Monodon_2016_BYSIZE[, SurveyYear := 2016]

M_2008_All_FARMS_Small <- Monodon_2008_BYSIZE[FarmScaleID == "Small", .N, .(FarmID, Disease)]
setnames(M_2008_All_FARMS_Small, "N", "B_S_N")
M_2008_All_FARMS_Medium <- Monodon_2008_BYSIZE[FarmScaleID == "Medium", .N, .(FarmID, Disease)]
setnames(M_2008_All_FARMS_Medium, "N", "B_M_N")
M_2008_All_FARMS_Large <- Monodon_2008_BYSIZE[FarmScaleID == "Large", .N, .(FarmID, Disease)]
setnames(M_2008_All_FARMS_Large, "N", "B_L_N")

M_2008_FARMS_DISEASECases_SMALL <- merge(M_2008_All_FARMS_Small[B_S_N == 1, .N, Disease], M_2008_All_FARMS_Small[B_S_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_S_S", "B_S_M"))
M_2008_FARMS_DISEASECases_MEDIUM <- merge(M_2008_All_FARMS_Medium[B_M_N == 1, .N, Disease], M_2008_All_FARMS_Medium[B_M_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_M_S", "B_M_M"))
M_2008_FARMS_DISEASECases_LARGE <- merge( M_2008_All_FARMS_Large[B_L_N == 1, .N, Disease], M_2008_All_FARMS_Large[B_L_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_L_S", "B_L_M"))
M_2008_FARMS_DISEASECases <- merge(M_2008_FARMS_DISEASECases_SMALL, M_2008_FARMS_DISEASECases_MEDIUM, by = "Disease", all = TRUE)
M_2008_FARMS_DISEASECases <- merge(M_2008_FARMS_DISEASECases, M_2008_FARMS_DISEASECases_LARGE, by = "Disease", all = TRUE)

M_2016_All_FARMS_Small <- Monodon_2016_BYSIZE[FarmScaleID == "Small", .N, .(FarmID, Disease)]
setnames(M_2016_All_FARMS_Small, "N", "B_S_N")
M_2016_All_FARMS_Medium <- Monodon_2016_BYSIZE[FarmScaleID == "Medium", .N, .(FarmID, Disease)]
setnames(M_2016_All_FARMS_Medium, "N", "B_M_N")
M_2016_All_FARMS_Large <- Monodon_2016_BYSIZE[FarmScaleID == "Large", .N, .(FarmID, Disease)]
setnames(M_2016_All_FARMS_Large, "N", "B_L_N")

M_2016_FARMS_DISEASECases_SMALL <- merge(M_2016_All_FARMS_Small[B_S_N == 1, .N, Disease], M_2016_All_FARMS_Small[B_S_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_S_S", "B_S_M"))
M_2016_FARMS_DISEASECases_MEDIUM <- merge(M_2016_All_FARMS_Medium[B_M_N == 1, .N, Disease], M_2016_All_FARMS_Medium[B_M_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_M_S", "B_M_M"))
M_2016_FARMS_DISEASECases_LARGE <- merge( M_2016_All_FARMS_Large[B_L_N == 1, .N, Disease], M_2016_All_FARMS_Large[B_L_N > 1, .N, Disease], by = "Disease", all = TRUE, suffixes = c("B_L_S", "B_L_M"))
M_2016_FARMS_DISEASECases <- merge(M_2016_FARMS_DISEASECases_SMALL, M_2016_FARMS_DISEASECases_MEDIUM, by = "Disease", all = TRUE)
M_2016_FARMS_DISEASECases <- merge(M_2016_FARMS_DISEASECases, M_2016_FARMS_DISEASECases_LARGE, by = "Disease", all = TRUE)

M_ALL_FARMS_DISEASECases <- merge(M_2008_FARMS_DISEASECases, M_2016_FARMS_DISEASECases, by = "Disease", all = TRUE)
M_ALL_FARMS_DISEASECases[is.na(M_ALL_FARMS_DISEASECases)] <- 0

# Get the totals
TotalVals <- M_ALL_FARMS_DISEASECases[, .(sum(NB_S_S.x), sum(NB_S_M.x),
                                          sum(NB_M_S.x), sum(NB_M_M.x),
                                          sum(NB_L_S.x), sum(NB_L_M.x),
                                          sum(NB_S_S.y), sum(NB_S_M.y),
                                          sum(NB_M_S.y), sum(NB_M_M.y),
                                          sum(NB_L_S.y), sum(NB_L_M.y))]
M_ALL_FARMS_DISEASECases <- rbind(M_ALL_FARMS_DISEASECases,
                                  c("Totals", TotalVals[1,]), use.names = FALSE)

DiseaseByFarmTable_Size <- kable(M_ALL_FARMS_DISEASECases, col.names = c("Disease", 
                                                                         "Single Incident, Small Farm 2008", 
                                                                         "Multiple Incidents, Small Farm 2008",
                                                                         "Single Incident, Medium Farm 2008", 
                                                                         "Multiple Incidents, Medium Farm 2008",
                                                                         "Single Incident, Large Farm 2008", 
                                                                         "Multiple Incidents, Large Farm 2008",
                                                                         "Single Incident, Small Farm 2016", 
                                                                         "Multiple Incidents, Small Farm 2016",
                                                                         "Single Incident, Medium Farm 2016", 
                                                                         "Multiple Incidents, Medium Farm 2016",
                                                                         "Single Incident, Large Farm 2016", 
                                                                         "Multiple Incidents, Large Farm 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(7), border_right = TRUE)
```



### Details of the number of diseases reported by each farm - Numbers PER FARM PER FARM SIZE

For each farm the number of different identifiable diseases reported on that farm was determined.

`r DiseaseByFarmTable_Size`


```{r Diseases Per Group Per Farm, include = FALSE}
###################################################################################
#  Output the number of disease groups per farm - number of different recorded diseases
#
# Note: this is for FARMS reporting a case of the disease 
#
# 1 diseases, more than 1 disease, and then singly up to maximum of 9
####################################################################################
AgentPerFarm_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_2008 <- AgentPerFarm_2008[,.N, FarmID]
AgentPerFarm_2016 <- AgentPerFarm_2016[,.N, FarmID]

NoAgentsPerFarm <- data.frame(
  Single_2008 = AgentPerFarm_2008[N == 1, .N],
  Single_2016 = AgentPerFarm_2016[N == 1, .N],
  Multiple_2008 = AgentPerFarm_2008[N > 1, .N],
  Multiple_2016 = AgentPerFarm_2016[N > 1, .N],
  Two_2008 = AgentPerFarm_2008[N == 2, .N],
  Two_2016 = AgentPerFarm_2016[N == 2, .N],
  Three_2008 = AgentPerFarm_2008[N == 3, .N],
  Three_2016 = AgentPerFarm_2016[N == 3, .N],
  Four_2008 = AgentPerFarm_2008[N == 4, .N],
  Four_2016 = AgentPerFarm_2016[N == 4, .N],
  Five_2008 = AgentPerFarm_2008[N == 5, .N],
  Five_2016 = AgentPerFarm_2016[N == 5, .N],
  Six_2008 = AgentPerFarm_2008[N == 6, .N],
  Six_2016 = AgentPerFarm_2016[N == 6, .N]
)
DiseaseAgentNumbersTable <- kable(NoAgentsPerFarm, col.names = c("Single Incident 2008", "Single Incident 2016", "Multiple Incident 2008","Multiple Incident 2016", "2 2008", "2 2016", "3 2008", "3 2016", "4 2008", "4 2016", "5 2008", "5 2016", "6 2008", "6 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)

```

### Details of the number of groups of agent reported by each farm - Numbers PER FARM

For each farm the number of different disease agent groups (Bacteria, bacteria/fungus, parasite, virus, other, unknown) was determined.

`r DiseaseAgentNumbersTable`

```{r Infectious Diseases Per Group Per Farm, include = FALSE}
###################################################################################
#  Output the number of disease groups per farm - number of different recorded diseases
#
# HERE ONLY THE INFECTIOUS DISEASES
#
# Note: this is for FARMS reporting a case of the disease 
#
# 1 diseases, more than 1 disease, and then singly up to maximum of 9
####################################################################################
AgentPerFarm_INF__2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Infectious", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Infectious", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2008 <- AgentPerFarm_INF__2008[,.N, FarmID]
AgentPerFarm_INF__2016 <- AgentPerFarm_INF__2016[,.N, FarmID]

NoAgentsPerFarm_INF <- data.frame(
  Single_2008 = AgentPerFarm_INF__2008[N == 1, .N],
  Single_2016 = AgentPerFarm_INF__2016[N == 1, .N],
  Multiple_2008 = AgentPerFarm_INF__2008[N > 1, .N],
  Multiple_2016 = AgentPerFarm_INF__2016[N > 1, .N],
  Two_2008 = AgentPerFarm_INF__2008[N == 2, .N],
  Two_2016 = AgentPerFarm_INF__2016[N == 2, .N],
  Three_2008 = AgentPerFarm_INF__2008[N == 3, .N],
  Three_2016 = AgentPerFarm_INF__2016[N == 3, .N],
  Four_2008 = AgentPerFarm_INF__2008[N == 4, .N],
  Four_2016 = AgentPerFarm_INF__2016[N == 4, .N]
)
DiseaseAgentNumbersTable_INF <- kable(NoAgentsPerFarm_INF, col.names = c("Single Incident 2008", "Single Incident 2016", "Multiple Incident 2008","Multiple Incident 2016", "2 2008", "2 2016", "3 2008", "3 2016", "4 2008", "4 2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  column_spec(1, bold = TRUE, border_right = TRUE, background = 'grey') %>%
  column_spec(c(2,3), border_right = TRUE)

```

### Details of the number of groups of infectious agent reported by each farm - Numbers PER FARM

For each farm the number of different disease agent groups (Bacteria, bacteria/fungus, parasite, virus) was determined.

`r DiseaseAgentNumbersTable_INF`


```{r Infectious Diseases Per Group Per Farm BAR CHART, include = FALSE}
###################################################################################
#  Output the number of disease groups per farm - number of different recorded diseases
#
# HERE A BAR CHART for ONLY THE INFECTIOUS DISEASES = here relative numbers
#
# Note: this is for FARMS reporting a case of the disease 
#
# 1 diseases, more than 1 disease, and then singly up to maximum of 9
####################################################################################
library(plotly)
AgentPerFarm_INF__2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Infectious", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Infectious", Agent, .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2008 <- AgentPerFarm_INF__2008[,.N, FarmID]
AgentPerFarm_INF__2016 <- AgentPerFarm_INF__2016[,.N, FarmID]

NoAgentsPerFarm_INF <- data.frame(cbind(
  Cat = c("> 1",
          "1",
          "2",
          "3",
          "4"),
  Value_2008 = c(round(AgentPerFarm_INF__2008[N > 1, .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
                 round(AgentPerFarm_INF__2008[N == 1, .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
                 round(AgentPerFarm_INF__2008[N == 2, .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
                 round(AgentPerFarm_INF__2008[N == 3, .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
                 round(AgentPerFarm_INF__2008[N == 4, .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2)
  ),
  Value_2016 = c(
    round(AgentPerFarm_INF__2016[N > 1, .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 1, .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 2, .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 3, .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 4, .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2)
  )
))

AgentPerFarmPlot <- plot_ly() %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF,
            x = ~Cat,
            y = ~Value_2008,
            marker = list(color = 'rgb(127,201,127)',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2008') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF,
            x = ~Cat,
            y = ~Value_2016,
            marker = list(color = 'rgb(190,174,212)',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2016') %>%
  layout(
    xaxis = list(type = 'category', showline = TRUE, title = "Recorded Disease Groups Per Farm"),
    yaxis = list(type = 'linear', showline = TRUE, title = "Relative number of Farms"),
    legend = list(x = 0.9, y = 0.92, 
                  title = list(text = '<b>Year</b>'),
                  bgcolor = 'rgb(200,200,200)',
                  bordercolor = 'rgb(200,200,200)',
                  borderwidth = 2)
  )

AgentPerFarmPlot
```

### Charts of the Diseases

The following is a bar chart of the number of infectious disease groups recorded by farms. The values are relative value to the number of recorded groups by the farms in that year.

`r AgentPerFarmPlot`


```{r Infectious Diseases Per Group Per Farm BAR CHART STACKED, include = FALSE}
###################################################################################
#  Output the number of disease groups per farm - number of different recorded diseases
#
# HERE A BAR CHART for ONLY THE INFECTIOUS DISEASES = here relative numbers
#
# Note: this is for FARMS reporting a case of the disease 
#
# 1 diseases, more than 1 disease, and then singly up to maximum of 9
####################################################################################
library(plotly)
AgentPerFarm_INF__2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Infectious", .(Agent, FarmScaleID), .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Infectious", .(Agent, FarmScaleID), .(FarmID)], by = c("FarmID", "Agent"))
AgentPerFarm_INF__2008 <- AgentPerFarm_INF__2008[,.N, .(FarmID, FarmScaleID)]
AgentPerFarm_INF__2016 <- AgentPerFarm_INF__2016[,.N, .(FarmID, FarmScaleID)]

NoAgentsPerFarm_INF_STK <- data.frame(cbind(
  Cat = c("> 1",
          "1",
          "2",
          "3",
          "4"),
  Value_2008_Small = c(
    round(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 2 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 3 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 4 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2)
  ),
  Value_2008_Medium = c(
    round(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 2 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 3 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 4 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2)
  ),
  Value_2008_Large = c(
    round(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 2 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 3 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2),
    round(AgentPerFarm_INF__2008[N == 4 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2008[, .N]), digits = 2)
  ),
  Value_2016_Small = c(
    round(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 2 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 3 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 4 & FarmScaleID == "Small", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2)
  ),
  Value_2016_Medium = c(
    round(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 2 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 3 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 4 & FarmScaleID == "Medium", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2)
  ),
  Value_2016_Large = c(
    round(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 2 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 3 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2),
    round(AgentPerFarm_INF__2016[N == 4 & FarmScaleID == "Large", .N]/sum(AgentPerFarm_INF__2016[, .N]), digits = 2)
  )
))

AgentPerFarmPlot_Stacked <- plot_ly() %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2008_Small,
            yaxis = 'y',
            width = 0.3,
            marker = list(color = '#ccece6',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2008 Small') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2008_Medium,
            yaxis = 'y',
            width = 0.3,
            marker = list(color = '#66c2a4',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2008 Medium') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2008_Large,
            yaxis = 'y',
            width = 0.3,
            marker = list(color = '#238b45',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2008 Large') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2016_Small,
            yaxis = 'y2',
            width = 0.3,
            offset = 0.1,
            marker = list(color = '#e0ecf4',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2016 Small') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2016_Medium,
            yaxis = 'y2',
            width = 0.3,
            offset = 0.1,
            marker = list(color = '#9ebcda',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2016 Medium') %>%
  add_trace(type = 'bar',
            data = NoAgentsPerFarm_INF_STK,
            x = ~Cat,
            y = ~Value_2016_Large,
            yaxis = 'y2',
            width = 0.3,
            offset = 0.1,
            marker = list(color = '#8c6bb1',
                          line = list(color = 'rgb(50,50,50)',
                                      width = 1)),
            name = '2016 Large') %>%
  layout(
    barmode = 'stack',
    bargap = 0,
    xaxis = list(type = 'category', showline = TRUE, title = "Recorded Disease Groups Per Farm"),
    yaxis = list(type = 'linear', showline = TRUE, title = "Relative number of Farms",
                 side = 'left'),
    yaxis2 = list(type = 'linear', showline = TRUE, title = "Relative number of Farms",
                  side = 'left', overlaying = "y", visible = FALSE),
    legend = list(x = 0.9, y = 0.92, 
                  title = list(text = '<b>Farms</b>'),
                  bgcolor = 'rgb(200,200,200)',
                  bordercolor = 'rgb(200,200,200)',
                  borderwidth = 2)
  )

AgentPerFarmPlot_Stacked

######################################################
#
# Perform stats 2 x 3 table
#####################################################
AffectFarmSize_2008 <- data.frame(
  Small = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Small", .N], AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Small", .N]),
  Medium = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Medium", .N], AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Medium", .N]),
  Large = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Large", .N], AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Large", .N])
)

Chi_Size_2008 <- chisq.test(AffectFarmSize_2008)

AffectFarmSize_2016 <- data.frame(
  Small = c(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Small", .N],AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Small", .N]),
  Medium = c(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Medium", .N], AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Medium", .N]),
  Large = c(AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Large", .N], AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Large", .N])
)

Chi_Size_2016 <- chisq.test(AffectFarmSize_2016)

AffectFarmSize_SMALL <- data.frame(
  Single = c(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Small", .N],AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Small", .N]),
  Multiple = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Small", .N], AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Small", .N])
)
AffectFarmSize_MEDIUM <- data.frame(
  Single = c(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Medium", .N],AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Medium", .N]),
  Multiple = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Medium", .N], AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Medium", .N])
)
AffectFarmSize_LARGE <- data.frame(
  Single = c(AgentPerFarm_INF__2008[N == 1 & FarmScaleID == "Large", .N],AgentPerFarm_INF__2016[N == 1 & FarmScaleID == "Large", .N]),
  Multiple = c(AgentPerFarm_INF__2008[N > 1 & FarmScaleID == "Large", .N], AgentPerFarm_INF__2016[N > 1 & FarmScaleID == "Large", .N])
)
Chi_Size_SMALL <- chisq.test(AffectFarmSize_SMALL)
Chi_Size_MEDIUM <- chisq.test(AffectFarmSize_MEDIUM)
Chi_Size_LARGE <- chisq.test(AffectFarmSize_LARGE)
```

### Charts of the Diseases

The following is a bar chart of the number of infectious disease groups recorded by farms. The values are relative value to the number of recorded groups by the farms in that year.

`r AgentPerFarmPlot_Stacked`

### Stats on the effect of farm size on the number of disease types

A chi-squared test was performed to determine if the farm size affected whether it had a single recorded disease type or more than 1.  The results were:
  
  P-value for Small farms in 2008 vs 2016 was `r Chi_Size_SMALL$p.value` with residuals of `r Chi_Size_SMALL$residuals`

P-value for Medium farms in 2008 vs 2016 was `r Chi_Size_MEDIUM$p.value` with residuals of `r Chi_Size_MEDIUM$residuals`

P-value for Large farms in 2008 vs 2016 was `r Chi_Size_LARGE$p.value` with residuals of `r Chi_Size_LARGE$residuals`

### Stats on the effect of year on the number of disease types

A chi-squared test was performed to determine if the farm size affected whether it had a single recorded disease type or more than 1.  The results were:
  
  P-value for 2008 was `r Chi_Size_2008$p.value` with residuals of `r Chi_Size_2008$residuals`

P-value for 2016 was `r Chi_Size_2016$p.value` with residuals of `r Chi_Size_2016$residuals`

The first set of charts are sunburst charts to show the types of diseases 

```{r Sunburst charts for diseases, include=FALSE}
######################################################################################
#
# Prepared data here is all disease recorded which with single record per farm
#
# Repeat the sunburst plot but with farms affected
######################################################################################
library(plotly)
# Prepare the data for all agents - NEED TO BE 1 PER FARM
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon",
                                      .(Agent, Disease), .(FarmID)], by = c("FarmID", "Disease"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon", 
                                      .(Agent, Disease), .(FarmID)], by = c("FarmID", "Disease"))
setorder(InfAgents_2008, Agent)
setorder(InfAgents_2016, Agent)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "Bacteria", "Bacteria/Fungus", "Virus", "Parasite", "Other",
                 InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2008[Agent == "Other", .N, .(FarmID, Disease)][,.N, Disease]$Disease
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2008)<br>Farms",
               "Bacteria", "BacteriaFungus", "Virus", "Parasite", "Other",
               InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2008[Agent == "Other", .N, .(FarmID, Disease)][,.N, Disease]$Disease
)

# Produce a parent list
ParentList <-c("",
               "Year","Year","Year","Year","Year",
               rep("Bacteria", nrow(InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Bacteria/Fungus", nrow(InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Virus", nrow(InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Parasite", nrow(InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Other", nrow(InfAgents_2008[Agent == "Other", .N, .(FarmID, Disease)][,.N, Disease])))

# Ensure segmentation is correct
Values <- c(sum(InfAgents_2008[, .N, .(FarmID, Disease)][,.N, Disease]$N), # Parent
            InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N],
            InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N],
            InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N],
            InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N],
            InfAgents_2008[Agent == "Other", .N, .(FarmID, Disease)][,.N],
            InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#74c476","#99d8c9","#9e9ac8", "#fd8d3c","#6baed6", # 
                    rep("#c7e9c0", nrow(InfAgents_2008[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#ccece6", nrow(InfAgents_2008[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#dadaeb", nrow(InfAgents_2008[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#fdd0a2", nrow(InfAgents_2008[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#c6dbef", nrow(InfAgents_2008[Agent == "Other", .N, .(FarmID, Disease)][,.N, Disease]))
)

# Produce the dataframe
ED_Dataframe_2008_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "Bacteria", "Bacteria/Fungus", "Virus", "Parasite", "Other",
                 InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
                 InfAgents_2016[Agent == "None" | Agent == "Unknown" , .N, .(FarmID, Disease)][,.N, Disease]$Disease
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2016)<br>Farms",
               "Bacteria", "BacteriaFungus", "Virus", "Parasite", "Other",
               InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$Disease,
               InfAgents_2016[Agent == "None" | Agent == "Unknown" , .N, .(FarmID, Disease)][,.N, Disease]$Disease
)

# Produce a parent list
ParentList <-c("",
               "Year","Year","Year","Year","Year",
               rep("Bacteria", nrow(InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Bacteria/Fungus", nrow(InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Virus", nrow(InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Parasite", nrow(InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease])),
               rep("Other", nrow(InfAgents_2016[Agent == "None" | Agent == "Unknown" , .N, .(FarmID, Disease)][,.N, Disease])))

# Ensure segmentation is correct
Values <- c(sum(InfAgents_2016[, .N, .(FarmID, Disease)][,.N, Disease]$N), # Parent
            InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N],
            InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N],
            InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N],
            InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N],
            InfAgents_2016[Agent == "None" | Agent == "Unknown", .N, .(FarmID, Disease)][,.N],
            InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease]$N,
            InfAgents_2016[Agent == "None" | Agent == "Unknown", .N, .(FarmID, Disease)][,.N, Disease]$N
)

# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#74c476","#99d8c9","#9e9ac8", "#fd8d3c","#6baed6", # 
                    rep("#c7e9c0", nrow(InfAgents_2016[Agent == "Bacteria", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#ccece6", nrow(InfAgents_2016[Agent == "Bacteria/Fungus", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#dadaeb", nrow(InfAgents_2016[Agent == "Virus", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#fdd0a2", nrow(InfAgents_2016[Agent == "Parasite", .N, .(FarmID, Disease)][,.N, Disease])),
                    rep("#c6dbef", nrow(InfAgents_2016[Agent == "None" | Agent == "Unknown", .N, .(FarmID, Disease)][,.N, Disease]))
)

# Produce the dataframe
ED_Dataframe_2016_Farm <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce the plot
SunburstPlot_Combined_Farm <- plot_ly() %>%
  add_trace(type = 'sunburst',
            ids = ED_Dataframe_2008_Farm$SegmentList, 
            labels = ED_Dataframe_2008_Farm$LabelList, 
            parents = ED_Dataframe_2008_Farm$ParentList,
            values = ED_Dataframe_2008_Farm$Values,
            marker = list(colors = ED_Dataframe_2008_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 0))

# Produce the plot
SunburstPlot_Combined_Farm <- SunburstPlot_Combined_Farm %>%
  add_trace(type = 'sunburst',
            ids = ED_Dataframe_2016_Farm$SegmentList, 
            labels = ED_Dataframe_2016_Farm$LabelList, 
            parents = ED_Dataframe_2016_Farm$ParentList,
            values = ED_Dataframe_2016_Farm$Values,
            marker = list(colors = ED_Dataframe_2016_Farm$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 1))

SunburstPlot_Combined_Farm <- SunburstPlot_Combined_Farm %>%
  layout(grid = list(columns = 2, rows =1),
         margin = list(l = 0, r = 0, b = 0, t = 0))

SunburstPlot_Combined_Farm

######################################################################################
#
# Prepared data here is all disease recorded by each farm
#
# Repeat this now for different farm types - small, medium, large
# 
#######################################################################################
# Prepare the data for all agents - NEED TO BE 1 PER FARM
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon",
                                      .(Agent, Disease, FarmScaleID), .(FarmID)], by = "FarmID")
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon",
                                      .(Agent, Disease, FarmScaleID), .(FarmID)], by = "FarmID")
InfAgents_2008 <- InfAgents_2008[,.N, .(FarmScaleID, Agent, Disease)]
InfAgents_2016 <- InfAgents_2016[,.N, .(FarmScaleID, Agent, Disease)]
setorder(InfAgents_2008, Agent)
setorder(InfAgents_2016, Agent)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "Small", "Medium", "Large",
                 paste0("Small - ", InfAgents_2008[FarmScaleID == "Small",sum(N), Agent]$Agent),
                 paste0("Medium - ", InfAgents_2008[FarmScaleID == "Medium",sum(N), Agent]$Agent),
                 paste0("Large - ", InfAgents_2008[FarmScaleID == "Large",sum(N), Agent]$Agent),
                 paste0("Small - ", InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Disease),
                 paste0("Medium - ", InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Disease),
                 paste0("Large - ", InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Disease)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2008)<br>Farms",
               "Small", "Medium", "Large",
               InfAgents_2008[FarmScaleID == "Small",sum(N), Agent]$Agent,
               InfAgents_2008[FarmScaleID == "Medium",sum(N), Agent]$Agent,
               InfAgents_2008[FarmScaleID == "Large",sum(N), Agent]$Agent,
               InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Disease,
               InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Disease,
               InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Disease
)

# Produce a parent list
ParentList <-c("",
               "Year","Year","Year",
               rep("Small", length(InfAgents_2008[FarmScaleID == "Small",sum(N), Agent]$Agent)),
               rep("Medium", length(InfAgents_2008[FarmScaleID == "Medium",sum(N), Agent]$Agent)),
               rep("Large", length(InfAgents_2008[FarmScaleID == "Large",sum(N), Agent]$Agent)),
               paste0("Small - ", InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Agent),
               paste0("Medium - ", InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Agent),
               paste0("Large - ", InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Agent)
)

# Ensure segmentation is correct
Values <- c(sum(InfAgents_2008$N, na.rm = T), # Parent
            InfAgents_2008[FarmScaleID == "Small",sum(N, na.rm = T)], #Small farms
            InfAgents_2008[FarmScaleID == "Medium",sum(N, na.rm = T)], #Medium farms
            InfAgents_2008[FarmScaleID == "Large",sum(N, na.rm = T)], # Large farms
            InfAgents_2008[FarmScaleID == "Small",sum(N), Agent]$V1,
            InfAgents_2008[FarmScaleID == "Medium",sum(N), Agent]$V1,
            InfAgents_2008[FarmScaleID == "Large",sum(N), Agent]$V1,
            InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$V1,
            InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$V1,
            InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$V1
)
# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#f1a340","#f7f7f7","##998ec3", # 
                    rep("#fec44f", length(InfAgents_2008[FarmScaleID == "Small",sum(N), Agent]$Agent)),
                    rep("#bdbdbd", length(InfAgents_2008[FarmScaleID == "Medium",sum(N), Agent]$Agent)),
                    rep("#bcbddc", length(InfAgents_2008[FarmScaleID == "Large",sum(N), Agent]$Agent)),
                    rep("#fee6ce", length(InfAgents_2008[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$V1)),
                    rep("#636363", length(InfAgents_2008[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$V1)),
                    rep("#efedf5", length(InfAgents_2008[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$V1))
)

# Produce the dataframe
ED_Dataframe_2008_FarmByFarmScale <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)

# Produce a list of all the components, starting from the middle
SegmentList <- c("Year",
                 "Small", "Medium", "Large",
                 paste0("Small - ", InfAgents_2016[FarmScaleID == "Small",sum(N), Agent]$Agent),
                 paste0("Medium - ", InfAgents_2016[FarmScaleID == "Medium",sum(N), Agent]$Agent),
                 paste0("Large - ", InfAgents_2016[FarmScaleID == "Large",sum(N), Agent]$Agent),
                 paste0("Small - ", InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Disease),
                 paste0("Medium - ", InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Disease),
                 paste0("Large - ", InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Agent, " - ", InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Disease)
)

# Produce a list of labels
LabelList <- c("<I>Penaeus monodon</I><br>(2016)<br>Farms",
               "Small", "Medium", "Large",
               InfAgents_2016[FarmScaleID == "Small",sum(N), Agent]$Agent,
               InfAgents_2016[FarmScaleID == "Medium",sum(N), Agent]$Agent,
               InfAgents_2016[FarmScaleID == "Large",sum(N), Agent]$Agent,
               InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Disease,
               InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Disease,
               InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Disease
)

# Produce a parent list
ParentList <-c("",
               "Year","Year","Year",
               rep("Small", length(InfAgents_2016[FarmScaleID == "Small",sum(N), Agent]$Agent)),
               rep("Medium", length(InfAgents_2016[FarmScaleID == "Medium",sum(N), Agent]$Agent)),
               rep("Large", length(InfAgents_2016[FarmScaleID == "Large",sum(N), Agent]$Agent)),
               paste0("Small - ", InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$Agent),
               paste0("Medium - ", InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$Agent),
               paste0("Large - ", InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$Agent)
)

# Ensure segmentation is correct
Values <- c(sum(InfAgents_2016$N, na.rm = T), # Parent
            InfAgents_2016[FarmScaleID == "Small",sum(N, na.rm = T)], #Small farms
            InfAgents_2016[FarmScaleID == "Medium",sum(N, na.rm = T)], #Medium farms
            InfAgents_2016[FarmScaleID == "Large",sum(N, na.rm = T)], # Large farms
            InfAgents_2016[FarmScaleID == "Small",sum(N), Agent]$V1,
            InfAgents_2016[FarmScaleID == "Medium",sum(N), Agent]$V1,
            InfAgents_2016[FarmScaleID == "Large",sum(N), Agent]$V1,
            InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$V1,
            InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$V1,
            InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$V1
)
# Pick the segment colours
SegmentColours <- c("#ffffff", # Central
                    "#f1a340","#f7f7f7","##998ec3", # 
                    rep("#fec44f", length(InfAgents_2016[FarmScaleID == "Small",sum(N), Agent]$Agent)),
                    rep("#bdbdbd", length(InfAgents_2016[FarmScaleID == "Medium",sum(N), Agent]$Agent)),
                    rep("#bcbddc", length(InfAgents_2016[FarmScaleID == "Large",sum(N), Agent]$Agent)),
                    rep("#fee6ce", length(InfAgents_2016[FarmScaleID == "Small",sum(N), .(Agent, Disease)]$V1)),
                    rep("#636363", length(InfAgents_2016[FarmScaleID == "Medium",sum(N), .(Agent, Disease)]$V1)),
                    rep("#efedf5", length(InfAgents_2016[FarmScaleID == "Large",sum(N), .(Agent, Disease)]$V1))
)

# Produce the dataframe
ED_Dataframe_2016_FarmByFarmScale <- data.frame(SegmentList, LabelList, ParentList, SegmentColours, Values)
# Produce the plot
SunburstPlot_Combined_FarmsByFarmScale <- plot_ly() %>%
  add_trace(type = 'sunburst',
            ids = ED_Dataframe_2008_FarmByFarmScale$SegmentList, 
            labels = ED_Dataframe_2008_FarmByFarmScale$LabelList, 
            parents = ED_Dataframe_2008_FarmByFarmScale$ParentList,
            values = ED_Dataframe_2008_FarmByFarmScale$Values,
            marker = list(colors = ED_Dataframe_2008_FarmByFarmScale$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 0))

# Produce the plot
SunburstPlot_Combined_FarmsByFarmScale <- SunburstPlot_Combined_FarmsByFarmScale %>%
  add_trace(type = 'sunburst',
            ids = ED_Dataframe_2016_FarmByFarmScale$SegmentList, 
            labels = ED_Dataframe_2016_FarmByFarmScale$LabelList, 
            parents = ED_Dataframe_2016_FarmByFarmScale$ParentList,
            values = ED_Dataframe_2016_FarmByFarmScale$Values,
            marker = list(colors = ED_Dataframe_2016_FarmByFarmScale$SegmentColours),
            maxdepth = 4,
            branchvalues = 'total',
            domain = list(column = 1))

SunburstPlot_Combined_FarmsByFarmScale <- SunburstPlot_Combined_FarmsByFarmScale %>%
  layout(grid = list(columns = 2, rows =1),
         margin = list(l = 0, r = 0, b = 0, t = 0))

SunburstPlot_Combined_FarmsByFarmScale

```


## Disease incidence = farms reporting the disease

In these sunburst charts the number of farms that reported cases of the disease are shown. Here a farm that reports the disease more than once is considered as a single case.

`r SunburstPlot_Combined_Farm`


## Disease incidence = size of farm reporting the disease

In these sunburst charts the number of farms reporting the disease are shown for small, medium and large farms. Here a farm that reports the disease more than once is considered as a single case.

`r SunburstPlot_Combined_FarmsByFarmScale`

```{r Statistical Tests on the Disease Data, include=FALSE}
##########################################################
#
# Statistical tests on the data
#
# We have count data so Fisher or Chi squared methods best
#
###########################################################
library(kableExtra)

##################################################################
#
# NATURE OF INJURY
#
####################################################################
# Prepare the data for just infectious agents - per farm
# Prepare the data for all agents - NEED TO BE 1 PER FARM
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, I := 2008]
InfAgents_2016[, I := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
#Inf_Fisher <- fisher.test(InfAgents)
Inf_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
ChiSQTableDF <- cbind(Inf_ChiSq$observed, Inf_ChiSq$expected, Inf_ChiSq$residuals)
ChiSQTable <- kable(ChiSQTableDF, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))

######################################################################
#
# Disease agents
#
#######################################################################
# Prepare the data for all agents - NEED TO BE 1 PER FARM
# Prepare the data for all agents - NEED TO BE 1 PER FARM
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, Year := 2008]
InfAgents_2016[, Year := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
AllAgents_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
AllAgentsChiSQTableDF <- cbind(AllAgents_ChiSq$observed, AllAgents_ChiSq$expected, AllAgents_ChiSq$residuals)
AllAgentsChiSQTable <- kable(AllAgentsChiSQTableDF, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))

```
## Statistical Tests on differences between 2008 and 2016

### Differences in the difference in farms recording disease agents between 2008 and 2016

In this test the number of cases of disease for each agent is obtained. The categories are Bacteria, Bacteria/fungus, virus, parasite and including unknown or no infectious agent.

The table of results is:
  
  `r AllAgentsChiSQTable`

Using the chi-squared test of independence to observed whether the number of counts in each category are significantly different between 2008 and 2016.

The chi-squared P value is `r AllAgents_ChiSq$p.value` showing a significant difference between the 2. This is borne out by comparison of the observed and expected values.


```{r Stats Part 2, include=FALSE}
######################################################################
#
# JUst INFECTIUOS AGENTS ON FARMS
#
#######################################################################
# Prepare the data for all agents - NEED TO BE 1 PER FARM
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Cause == "Infectious",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Cause == "Infectious",
                                      .(Agent), .(FarmID)], by = c("FarmID","Agent"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, Year := 2008]
InfAgents_2016[, Year := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
AllAgents2_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
ChiSQTableInfectiousDF <- cbind(AllAgents2_ChiSq$observed, AllAgents2_ChiSq$expected, AllAgents2_ChiSq$residuals)
ChiSQTableInfectious <- kable(ChiSQTableInfectiousDF, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))
```

### Differences in the frequency of different types of infectious disease

In this test the number of cases of infectious disease for each agent is obtained. The categories are Bacteria, Bacteria/fungus, virus, parasite. The cases of diseases with unknown or no infectious agent are removed from the analysis - we need to be careful about this if there is some infectious part to these injuries.

Using the chi-squared test of independence to observed whether the number of counts in each category are significantly different between 2008 and 2016.

The chi-squared P value is `r AllAgents2_ChiSq$p.value` showing a significant difference between the 2. This is borne out by comparison of the observed and expected values.

`r ChiSQTableInfectious`


```{r Stats Part 3 By Size, include=FALSE}
########################################################################
# Stratify the data according to farm types
#
# Perform the Cochran-Mansel-Han test
#
########################################################################

#########################################################################################
#
# Individual farm scale level test
#
#########################################################################################
# Prepare the data for just infectious agents LARGE
# THIS IS PER FARM 
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Large",
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Large", 
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, I := 2008]
InfAgents_2016[, I := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
Inf_LARGE_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
ChiSQTableDF_LARGE <- cbind(Inf_LARGE_ChiSq$observed, Inf_LARGE_ChiSq$expected, Inf_LARGE_ChiSq$residuals)
ChiSQTable_LARGE <- kable(ChiSQTableDF_LARGE, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))

########################################################
# Prepare the data for just infectious agents MEDIUM
#########################################################
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Medium",
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Medium", 
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, I := 2008]
InfAgents_2016[, I := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
Inf_MEDIUM_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
ChiSQTableDF_MEDIUM <- cbind(Inf_MEDIUM_ChiSq$observed, Inf_MEDIUM_ChiSq$expected, Inf_MEDIUM_ChiSq$residuals)
ChiSQTable_MEDIUM <- kable(ChiSQTableDF_MEDIUM, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))

########################################################
# Prepare the data for just infectious agents SMALL
#########################################################
InfAgents_2008 <- unique(RawData_2008[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Small",
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2016 <- unique(RawData_2016[SpeciesID == "Penaeus monodon" & Agent != "None" & Agent != "Unknown"
                                      & FarmScaleID == "Small", 
                                      .(Agent, Disease), .(FarmID)], by = c("Disease", "FarmID"))
InfAgents_2008 <- InfAgents_2008[,.(Agent)]
InfAgents_2016 <- InfAgents_2016[,.(Agent)]

# Prepare the contingency table
InfAgents_2008[, I := 2008]
InfAgents_2016[, I := 2016]
InfAgents <- rbind(InfAgents_2008, InfAgents_2016)
colnames(InfAgents) <- c("Agent", "Year")
InfAgents <- table(InfAgents)

# Peform the fisher test & chi-squared test
Inf_SMALL_ChiSq <- chisq.test(InfAgents)

# Produce table of observed and expected values
ChiSQTableDF_SMALL <- cbind(Inf_SMALL_ChiSq$observed, Inf_SMALL_ChiSq$expected, Inf_SMALL_ChiSq$residuals)
ChiSQTable_SMALL <- kable(ChiSQTableDF_SMALL, col.names = c("2008", "2016", "2008", "2016", "2008", "2016")) %>%
  kable_styling(bootstrap_options = 'striped', full_width = TRUE) %>%
  add_header_above(c(" " = 1, "Observed" = 2, "Expected" = 2, "Residuals" = 2))
```

### Stratified according to farm type

Stratifying the data according to the farm type the chi-squared tests are:
  #######################################################
THIS HAS BEEN CORRECTED FOR 1 RECORD PER FARM
#######################################################
#### Large Farms

`r ChiSQTable_LARGE`

The p-value is `r Inf_LARGE_ChiSq$p.value`



#### Medium Farms

`r ChiSQTable_MEDIUM`

The p-value is `r Inf_MEDIUM_ChiSq$p.value`



#### Small Farms

`r ChiSQTable_SMALL`

The p-value is `r Inf_SMALL_ChiSq$p.value`