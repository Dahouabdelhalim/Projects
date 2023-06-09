
# ******************************************************************************
# Main script 
#
version$version.string
# R version 3.6.3 (2020-02-29)
#
# ******************************************************************************

# ===== Table of contents ======================================================

# A Overview
# B Packages 
# C Repo structure
# D Global variables

# ==== A Overview ==============================================================
#
# The following contents within this main directory include the packages, 
# libraries, dataset and the workflow that is used in this analysis.
# scripts should be run in the folling order for replicatory results:

#   0.main.R:data organization with the overview of the workflow used
#                   in this project
#   1.data.manip.R: size correction and outlier removal
#   2a.analysis.into.R: 
#   2b.analysis.PCA.R: PCA analysis on the raw data (excl mass) + 
#     mixed-model analysis to look at effect of sex, and presence of competition 
#     on response variables while controlling for island effects 
#   2c.analysis.SVL.cov.R: mixed-effects model as above but with SVL added as 
#     covariate to ccontrol for correlations among traits
#   2d.analysis.size.corrected.R: mixed-effects model on the size corrected data

#   3.figures.R: script to make the figures

# *** IMPORTANT workflow notes *************************************************
# 1) after running this main.R script, move 'all.traits.clean.csv'
#     into the folder '0.data.raw'
# 2) run the scripts in sequence (start with 0 -> 1 -> 2a -> 2b ->..)

# ==== B Packages ==============================================================

# install.packages("lme4")
# install.packages("MuMIn")
# install.packages("corrplot")

# load packages
library(lme4)
library(MuMIn) # to get R^2 mixed-effect models  https://besjournals.
#             onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
library(corrplot) # to plot correlation among variables

# === C Repo structure =========================================================

# The file of the raw data can be accessed here.
working.dir <- getwd()

# In this working directory there 4 pathways to different folders with 
# specific outputs 
#           - 0.data.raw  -> the path to this folder is: rd.path
#                 (This folder contains a copy of the original downloaded
#                   data set without changes.)

#           - 1.data.clean -> the path to this folder is: cd.path
#                 (This folder contains  the cleaned data set.)
#
#           - 2.results -> the path to this folder is: an.path
#                 (This folder contains saved outputs of our analysis.)
#
#           - 3.figures -> the path to this folder is: fi.path
#                 (This folder contains all visualizations of data.)

#create the folders
output.folders <- c("0.data.raw", "1.data.clean","2.results","3.figures")
# Make the folders using this loop code 
for(i in 1:length(output.folders)) 
  if(file.exists(output.folders[i]) == FALSE) 
    dir.create(output.folders[i])

# loop checks the output.folders list and checks to see 
# if the folders exist in the working directory. If they don't it will create 
# them. 
# Path to 0.data.raw
rd.path <- paste(working.dir,"/",output.folders[1], "/", sep="")

# Path to 1.data.clean
cd.path <- paste(working.dir,"/",output.folders[2], "/", sep="")

# Path to 2.results
re.path <- paste(working.dir,"/",output.folders[3], "/", sep="")

# Path to 3.figures
fi.path <- paste(working.dir,"/",output.folders[4], "/", sep="")

# ==== D Global variables ======================================================
# columns with the explanatory variables
location.expl.var <- c(2,3,4)

# ******************************************************************************
# ******************************************************************************





