
source("01_extracting_matrices_env.R") 

#creates a function used to create W,B,A matrices


source("02_PCAanalysis_cleaned data.R")

#updates previous trait axes PCA

#there is no "03" script


source("04_Single_site_convergence.R") 

#function for single site syncsa


source("05_All_region_convergence.R") 

#function for multi site syncsa


source("06_Run_SYNCSA_single_multi_site.R")

#main code that runs SYNCSA for each environmental gradient


source("06.1_Run_SYNCSA_north.R")

#repeating code 06 for sites north of northern Andes, detrital gradient


source("06.2_Run_SYNCSA_north.R")

#repeating code 06 for sites south of northern Andes, detrital gradient


source("07_site_comparison_regression_models.R")

#regression models of site differences in strength and sig of SYNCSA
