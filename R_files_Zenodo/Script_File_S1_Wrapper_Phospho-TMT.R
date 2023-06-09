## # # # # # # # # # # # # # # # # # # # # # # # #
##
## Project: CDKL5-NLS Phosphoproteomics
##
## Script name: Script_File_S1_Wrapper_Phospho-TMT.R
##
## Purpose of script: Wrapper for all R scripts used for 
## phosphoproteomics data analysis. Run this script for data analysis.
##
## Instructions: 
## 1) Extract "Data_Phospho-TMT.zip" into ./Data
## 2) Extract "Databases_Phospho-TMT.zip" into ./Databases
## 3) Ensure "averageMaxQuant.R" is placed in ./Functions
## 4) Install missing libraries ("Install necessary libraries"), 
##    version number used for original analysis
##    indicated in "Text_File_S1_Session_Info_Phospho-TMT.txt"
## 5) Run this wrapper script ("source scripts")
##
## Author: Florian Weiland
##
## Date Created: 2020-04-03
##
## # # # # # # # # # # # # # # # # # # # # # # # #

## Create missing directories ----

if (dir.exists("./Output_R/") == FALSE) {

dir.create("./Output_R/")

}

## Install necessary libraries ----

libs <- c(
  "ggplot2",
  "ggnetwork",
  "wesanderson",
  "ggrepel",
  "STRINGdb",
  "igraph",
  "reshape2",
  "vsn",
  "limma",
  "seqinr",
  "plyr",
  "stringr",
  "ggpointdensity",
  "extrafont",
  "scales",
  "GO.db",
  "gtools"
)

rlibs.install <- which(libs %in% rownames(installed.packages()) == FALSE)

if (length(rlibs.install) >= 1) {

  message(paste0("\\nPlease install ", libs[rlibs.install]))
  
  message("Version numbers used is stated in \\"Text_File_S1_Session_Info_Phospho-TMT.txt\\"")

}

## Source scripts ----

source("Script_File_S2_Data_analysis_Phospho-TMT.R")
source("Script_File_S3_GO_analysis_Phospho-TMT.R")
source("Script_File_S4_Interactor_analysis_Phospho-TMT.R")

# sink("Text_File_S1_Session_Info_Phospho-TMT.txt")
# sessionInfo()
# sink()
