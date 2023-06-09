# This script analyzes over-time changes in attitudes among Muslim participants
# and replicates Figure S-13. This script makes use of dataset muslim-attitudes.csv

# Change the working directory as necessary
setwd("~/Dropbox/iraq-experiment/science/replication/data")

library(tidyverse)

# Read in data ------------------------------------------------------------
df <- read.csv("muslim-attitudes.csv")

# Analysis ---------------------------------------------------------------
#Regressing outcomes on time period (endline vs. baseline), using binary variables 

# Empathy toward Christians

m1 <- lm(empathy_t0_new ~ time, data = df)

# Pride in being Iraqi

m2 <- lm(proud_iraqi_t0_new ~ time, data = df)

# Would trust a Christian with a cash transfer

m3 <- lm(cash_christian_t0 ~ time, data = df)

# Iraq would be a better society if we treat each other as Iraqis first

m4 <- lm(secular_t0_new ~ time, data = df)

# Have some or many non-Muslim friends  

m5 <- lm(friends_t0_new ~ time, data = df)

# Ethnic and religious divisions in Iraq are arbitrary

m6 <- lm(secular1_t0_new ~ time, data = df)

# Would consider selling land to non-Muslims

m7 <- lm(land_t0_new ~ time, data = df)

# Would be OK with Christians as neighbors

m8 <- lm(nbr_chr_t0 ~ time, data = df)

# Feel comfortable in non-Muslim areas 

m9 <- lm(areas_t0 ~ time, data = df)


# Visualizing differences in means -----------------------------------


# Coefficients

means <- c(summary(m1)$coefficients[2,1],
           summary(m2)$coefficients[2,1],
           summary(m3)$coefficients[2,1],
           summary(m4)$coefficients[2,1],
           summary(m5)$coefficients[2,1],
           summary(m6)$coefficients[2,1],
           summary(m7)$coefficients[2,1],
           summary(m8)$coefficients[2,1],
           summary(m9)$coefficients[2,1]
)

# Standard errors

ses <- c(summary(m1)$coefficients[2,2],
         summary(m2)$coefficients[2,2],
         summary(m3)$coefficients[2,2],
         summary(m4)$coefficients[2,2],
         summary(m5)$coefficients[2,2],
         summary(m6)$coefficients[2,2],
         summary(m7)$coefficients[2,2],
         summary(m8)$coefficients[2,2],
         summary(m9)$coefficients[2,2]
         
)


names <- c("empathy_christians", "proud_iraqi", "trust_christian_cash_transfer",
           "coexistence_possible", "mixed_friends","religious_divisions_arbitrary",
           "sell_land_christians",
           "christian_neighbor_ok", "comfortable_mixed_areas")  


merged <- data.frame(names, means, ses)

# Adding lower and upper bounds

merged <- merged %>%  mutate(lb = means - ses * 1.96, 
                             ub = means + ses * 1.96) 

# Plotting difference in means

ggplot(merged, aes(x = reorder(names, - means), y = means)) + 
  geom_line(aes(group = names), color = "dark grey" ) + 
  geom_errorbar(aes(ymin = lb, ymax = ub), width = 0, color = "grey") + 
  theme_minimal() +
  geom_point(color = "black") + 
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(text = element_text(size=15))  +
  xlab("") +
  ylab("")


