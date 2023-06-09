################################################################################

#made with R version 4.0.3 (2020-10-10) on "x86_64-w64-mingw32" platform
#R version 4.1.3 (2022-03-10)
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
#################################################

#Replicating paper's results

#Table 5: Main Results: Comparison of prices, quantities, and surplus with at and RTP pricing.
#Figure 4: Surplus gain from real time pricing under different policy, cost and demand flexibility scenarios.
#Figure 5: Cost of 100 percent renewable energy system under different policy, cost and demand flexibility scenarios.
#Figure 7: Distributions of prices and quantities.
#Figure S1a and S1b (in appendix): Share of hours with Marginal Cost less than 1 $/MWh
#Tables S1-S6
source(file = "Tab5_Fig4-5and7.R")

#Figure 6: The social cost of clean power relative to a fossil future with at pricing.
source(file = "Fig6_RPSvalue_relativetoleastcost.R")

#Figure 8: Hourly production and consumption profiles for several scenarios with pessimistic interhour demand flexibility.
source(file = "Fig8.R")

#Figure 9: Production and consumption shares by sources with pessimistic interhour demand flexibility.
source(file = "Fig9_ProdConsShares.R")

#Figures S2 in appendix
source(file = "FigS2_productionconsumption13daysplot.R")

########END########################################################################


