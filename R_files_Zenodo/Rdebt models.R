setwd("E:/") # Set working directory

library(metafor)

data = read.csv("rDebt.csv",header=T) # Read database including error measures 

#### DATA MANIPULATION (Peter Jones) ####
#subsetting data by metric type
abundance = data[data$Metric.type == "Abundance",]
diversity = data[data$Metric.type == "Diversity",]
carbon = data[data$Metric.type == "C" | data$Metric.type == "Organic matter",] 
carbon$Metric.type = "Carbon"
nitrogen = data[data$Metric.type == "N",]
nitrogen$Metric.type = "Nitrogen"

###Some subsets require further editing due to our sample size requirements 
###per category and per study, mostly for the disturbance type subsets

#nitrogen-Wetland only has 9 RVs in it, so should be removed from analysis
nitrogen.hab = nitrogen[-which(nitrogen$HabitatCat_short == "Wetland"),]; nitrogen.hab$HabitatCat_short = factor(nitrogen.hab$HabitatCat_short)

#Invasive sp and Multiple do not have enough data
diversity.dis = diversity[diversity$DisturbCat != "Invasive species" & diversity$DisturbCat != "Multiple",]
diversity.dis$DisturbCat = factor(diversity.dis$DisturbCat)
#Invasive sp, hydro, multiple, oil, and overfishing do not have enough data
carbon.dis = carbon[carbon$DisturbCat != "Invasive species" & carbon$DisturbCat != "Multiple" & carbon$DisturbCat != "Hydrological disruption" & carbon$DisturbCat != "Oil" & carbon$DisturbCat != "Overfishing", ]
carbon.dis$DisturbCat = factor(carbon.dis$DisturbCat)
#Hydro, multiple, oil, and overfishing do not have enough data
nitrogen.dis = nitrogen[nitrogen$DisturbCat != "Multiple" & nitrogen$DisturbCat != "Hydrological disruption" & nitrogen$DisturbCat != "Oil" & nitrogen$DisturbCat != "Overfishing",]
nitrogen.dis$DisturbCat = factor(nitrogen.dis$DisturbCat)

#### META-ANALYSES (José López-López) ####
# Nitrogen.hab
Data = nitrogen.hab # Select database
Data$vi = 1.44 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short) - 1,random = ~ 1|Citation, data = Data)
m1b

# Nitrogen.dis
Data = nitrogen.hab # Select database
Data$vi = 1.44 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat) - 1,random = ~ 1|Citation, data = Data)
m1b

# Carbon
Data = carbon # Select database
Data$vi = 3.09 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short) - 1,random = ~ 1|Citation, data = Data)
m1b

# Carbon.dis
Data = carbon.dis # Select database
Data$vi = 1.56 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat) - 1,random = ~ 1|Citation, data = Data)
m1b

# Diversity
Data = diversity # Select database
Data$vi = 26.92 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short) - 1,random = ~ 1|Citation, data = Data)
m1b

# Diversity.dis
Data = diversity.dis # Select database
Data$vi = 31.19 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat) - 1,random = ~ 1|Citation, data = Data)
m1b

# Abundance
Data = abundance # Select database
Data$vi = 5.29 # Input rough estimate of the average within-study variance
# First, we fit the 3-level meta-analytic model without moderators
m0 = rma.mv(mainDebt,vi,random = ~ 1|Citation, data = Data) 
m0
# Now, we test the significance of the FIRST moderator (look at the p-value of the QM test)
m1a = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short),random = ~ 1|Citation, data = Data)
m1a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m1b = rma.mv(mainDebt,vi,mods=~ factor(HabitatCat_short) - 1,random = ~ 1|Citation, data = Data)
m1b
# Now, we test the significance of the SECOND moderator (look at the p-value of the QM test)
m2a = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat),random = ~ 1|Citation, data = Data)
m2a
# Last, we fit a model without intercept to get the effect size estimate (and CI) for each category
m2b = rma.mv(mainDebt,vi,mods=~ factor(DisturbCat) - 1,random = ~ 1|Citation, data = Data)
m2b


