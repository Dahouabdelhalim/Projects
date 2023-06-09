##--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##
# Author: Sean M. Wineland (seanwineland@gmail.com)
# Purpose: Calculate the Biodiversity Value Index (BVI) from manuscript:
# Wineland, S.M., Fovargue, R., Gill, K.C., Rezapour, S., & Neeson, T.M. 2020. 
# Conservation planning in an uncertain climate: identifying projects that remain valuable and feasible across future scenarios
# People and Nature, 2 (4), XX-XX
# Date: 10/28/2020
##-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------##

# load necessary library
library(tidyverse)

# load data
# replace "WhereverYouSavedTheData" with the actual path
data <- read.csv("C:/WhereverYouSavedTheData/BVI_rawdata_rinput.csv")

# Note:
# This data set contains the average probability of occurrence for each fish species (meanprob column, summed across 5 pixels from SDM representing downstream reach),
# for each river reach below each reservoir for each climate scenario.
# Columns (x)Vals contains the probability-weighted term for that species in that reservoir/scenario (see tables 1 & 2 in the manuscript)
# For example, Rval column is the "Probability" x "Rweight" column

# the following code is used to calculate the BVI:

BVI.calc <- data %>% group_by(Reservoir, Scenario) %>%          # groups data by reservoir & scenario
  summarise_at(vars(Rvals, Gvals, Tvals, Svals), list(sum)) %>% # sums terms by reservoir & scenario
  rowwise() %>%                                                 # function to sum across rows in next line:
  mutate(sumvals = sum(Rvals, Gvals, Tvals, Svals)) %>%         # sums the R,G,T,S terms across rows to obtain the sum probability of occurrence 
  mutate(BVI = sumvals/124) %>%                                 # creates new column (BVI) by dividing column (sumvals) by the total number of species and scenarios (31*4)
  group_by(Scenario, Reservoir) %>%                             # group data by scenario and reservoir
  arrange(Scenario, Reservoir) %>%                              # arrange dataframe by scenario and reservoir                                 
  separate(Scenario, c("RCP", "GCM", "Year"), sep = "_")        # separate column scenario by "_" to create plot

# plot data
(plot <- qplot(BVI,Reservoir, data = BVI.calc, colour = factor(RCP), shape = factor(GCM),  size = I(3)) +
  scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999")) +
  scale_x_continuous(limits = c(0, 0.125)) +
  labs(shape="GCM", colour="RCP")) 

# Note:
# This plot contains the same data as, but is not arranged the same as Figure 4a in the manuscript.