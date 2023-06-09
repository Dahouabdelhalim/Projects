# This script analyzes ttitudes among fans (local residents of Study Site 2) 
# and replicates Table S-11. This script makes use of datasets fans-during-intervention.csv,
# and fans-after-intervention.csv

# Change the working directory as necessary
setwd("~/Dropbox/iraq-experiment/science/replication/data")

library(tidyverse)
library(stargazer)

# Read in data ------------------------------------------------------------
df <- read.csv("fans-during-intervention.csv")
df2 <- read.csv("fans-after-intervention.csv")

# Analysis: Baseline ---------------------------------------------------------------

# Regressing outcomes on various independent variables:
# age, gender, living within walking distance, having a family member playing, 
# and number of games watched 
# Note this "baseline" took place during the intervention

# Replicates columns 1 - 2

# View the leagues as having a positive influence on their community

m1 <- lm(league_positive_t0 ~ age + gender + walk
         + fam + games_watched, data = df)

# Believe that ethnic and religious divisions in Iraq are arbitrary

m2 <- lm(secular_t0 ~  gender + walk
         + fam + games_watched , data = df)

# Analysis: Endline ---------------------------------------------------------------

# Regressing outcomes on various independent variables:
# age, gender, living within walking distance, having a family member playing, 
# and number of games watched 
# Note this endline took place 3 months after the intervention 

# Replicates columns 3 - 5


# View the leagues as having a positive influence on their community

m3 <- lm(league_positive_t1 ~ gender + walk
         + fam + games_watched + attend, data = df2)

# Believe that ethnic and religious divisions in Iraq are arbitrary

m4 <- lm(secular_t1 ~ age + gender + walk
         + fam + games_watched + attend, data = df2)

# Would  prefer community activities like the soccer league to be mixed in the future

m5 <- lm(prefer_group ~ age + gender + walk
         + fam + games_watched + attend, data = df2)

# Producing Table S-11 ---------------------------------------------------------------


stargazer(m1, m2, m3, m4, m5)


