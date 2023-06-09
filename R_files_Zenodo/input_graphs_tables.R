################################################################################

#made with R version 4.0.3 (2020-10-10) on "x86_64-w64-mingw32" platform
#R version 4.1.3 (2022-03-10)
#it takes several hours.
#Real Time Pricing and the Cost of Clean Power
#Imelda, Matthias Fripp, Michael J. Roberts
#questions: imelda@graduateinstitute.ch

################################################################################
R.version
setwd("D:/RTP/Replication")
print(getwd())

###below script, you only need to run one time###
source(file = "0_install_library.R")
source(file = "0_load_library.R")
source(file = "0_importfunctions.R")

#Table 1: Assumptions about flexible demand and demand-side reserves
#Figure 1: Demand flexibility scenarios by hour and month
#Table 2: Generators Cost Assumption
#Table 3: Fuel Cost Assumption
#Figure 3: Average output and potential capacity of renewable energy sources on Oahu
source(file = "1_Tab1-2-3_Fig1and3.R")

## combine all results from HPC
## need files: demand_summary_final, dual_costs, energy_sources
# produce /energysourcesCO2*.csv for most of the results
source(file = "2a_combine_all_data_hourly.R")

# produce energysourcesCO2",thetaval,"_rps.csv" for Figure 6
source(file = "2b_combine_all_data_hourly_rps.R")

# need demand_response_summary_final
# produce demandsummaryrps.csv for Figure 6
source(file = "2c_combine_demand_summary.R")

########END########################################################################
