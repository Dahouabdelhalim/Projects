# -------------------------------------------------------------------------
# Analyzing eBird data mean arrival dates and data quantity/quality
# -------------------------------------------------------------------------
setwd("~/1 U of T/eBird/eBird Project Processed Data/Routine 7/Simulation")

library(caper)
library(tidyverse)
library(lme4)
library(GGally)

# ---------------------------------
# Stationary
# ---------------------------------
# ---------------------------------
# Read in data
# ---------------------------------
Stationary <- read_csv("Stationary/Simulation_Stationary_ebird_MAD_Jd_80_180_w_enddate_5_Occ_Thresh_0.1_2002-2019_20220204A.csv")
Stationary_orig <- Stationary
Stationary <- Stationary %>% 
  select(-nyrs, -min_Year, -max_Year, -Duration_Years, -Spp_n_GridID, -GridID_n_spp, -ngrids, -Common_Name_n_GridID)

# ----------------------------------
# Initial filtering
# ----------------------------------
Stationary <- Stationary %>% 
  filter(JdDetect >= 10) %>% 
  # filter(P_val_xmid <= 0.05) %>%
  filter(MAD_CI_range <= 40)



# ----------------------------------
# Model fit summary
# ----------------------------------
# Summarize by year
Yearly_coverage <- Stationary %>% 
  group_by(Year) %>% 
  summarize(Year_count = n())
ggplot(Yearly_coverage, aes(Year, Year_count)) +
  geom_point() +
  scale_x_continuous(n.breaks = 10)

# Species year summary
Species_coverage <- Stationary %>% 
  group_by(Common_Name, Year) %>% 
  summarize(Common_Name_count = n())
mean(Species_coverage$Common_Name_count)
sd(Species_coverage$Common_Name_count)

# Number of days detected
mean(Stationary$JdDetect)
sd(Stationary$JdDetect)

# MAD CI range
hist(Stationary$MAD_CI_range)
quantile(Stationary$MAD_CI_range)
# below 10
MAD_CI_below_10 <- Stationary %>% 
  filter(MAD_CI_range < 10)
(nrow(MAD_CI_below_10)/nrow(Stationary))*100
mean(MAD_CI_below_10$JdDetect)
sd(MAD_CI_below_10$JdDetect)
summary(MAD_CI_below_10$JdDetect)
# Above 10
MAD_CI_above_10 <- Stationary %>% 
  filter(MAD_CI_range >= 10)
mean(MAD_CI_above_10$JdDetect)
sd(MAD_CI_above_10$JdDetect)
summary(MAD_CI_above_10$JdDetect)

# prsq
mean(Stationary$Prsq)
median(Stationary$Prsq)
sd(Stationary$Prsq)
summary(Stationary$Prsq)



# ---------------------------------
# Filter data
# ---------------------------------
# Calculate Grid-specific sampling coverage for species
Species_Grid_nyrs <- Stationary %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Stationary <- left_join(Stationary, Species_Grid_nyrs, by = c("Common_Name", "GridID"))
Stationary <- Stationary %>%
  filter(nyrs >= 10)
nrow(Stationary)
# # First of all, let's filter for species where sampling started at a minimum of 2007
# Stationary <- Stationary %>%
#   filter(min_Year <= 2007)

Species_Grid_Duration_Years <- Stationary %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Stationary <- left_join(Stationary, Species_Grid_Duration_Years, by = c("Common_Name", "GridID"))
# We will use these to further pare down the dataset
# Save the original data.frame
# Stationary_orig <- Stationary
# length(unique(Stationary_orig$GridID))
# 113
# # First of all, let's filter for species where sampling started at a minimum of 2007
# Stationary <- Stationary %>%
#   filter(min_Year <= 2007)

# Filter for species-grids that have a duration of at least 15 years
Stationary <- Stationary %>% 
  filter(Duration_Years >= 15)
nrow(Stationary)

# Remove species that have less than 7 grids (i.e., 1-6)
Species_Grid_n <- Stationary %>% 
  select(Common_Name, GridID) %>%
  distinct() %>% 
  group_by(Common_Name) %>% 
  count(name = "Spp_n_GridID") %>% 
  ungroup()
Stationary <- left_join(Stationary, Species_Grid_n, by = c("Common_Name"))
Stationary <- Stationary %>% 
  filter(Spp_n_GridID >= 6)
# Check # of rows
nrow(Stationary)

# Now, we want to see how many species there are across each of the grids
# Look at the maximum number of species that occur in each grid cell
# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Max_n <- Stationary %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(GridID) %>% 
  count(name = "GridID_n_spp") %>% 
  ungroup()
# Remove grids that only have one species
Stationary <- left_join(Stationary, Species_Grid_Max_n, by = c("GridID"))
Stationary <- Stationary %>% 
  filter(GridID_n_spp > 1)
# Check # of rows
nrow(Stationary)

# Check out the number of grids that a species has for a year
Species_Years <- Stationary %>% 
  group_by(Common_Name, Year) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 7 grids
Stationary <- left_join(Stationary, Species_Years, by = c("Common_Name", "Year"))
Stationary <- Stationary %>% 
  filter(ngrids >= 6)
nrow(Stationary)

# ALSO remove species that only have one grid
Grid_Species_Max_n <- Stationary %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(Common_Name) %>% 
  count(name = "Common_Name_n_GridID")
# Remove grids that only have one species
Stationary <- left_join(Stationary, Grid_Species_Max_n)
Stationary <- Stationary %>% 
  filter(Common_Name_n_GridID > 1)
# Check # of rows
nrow(Stationary)

write_csv(Stationary, "Stationary/Simulation_Stationary_ebird_MAD_Jd_80_180_w_enddate_5_Occ_Thresh_0.1_2002-2019_20220204B.csv")




# ---------------------------------
# Shifting
# ---------------------------------
# ---------------------------------
# Read in data
# ---------------------------------
Shifting <- read_csv("Shifting/Simulation_Shifting_ebird_MAD_Jd_80_180_w_enddate_5_Occ_Thresh_0.1_2002-2019_20220204A.csv")
Shifting_orig <- Shifting
Shifting <- Shifting %>% 
  select(-nyrs, -min_Year, -max_Year, -Duration_Years, -Spp_n_GridID, -GridID_n_spp, -ngrids, -Common_Name_n_GridID)


# ----------------------------------
# Initial filtering
# ----------------------------------
Shifting <- Shifting %>% 
  filter(JdDetect >= 10) %>% 
  # filter(P_val_xmid <= 0.05) %>%
  filter(MAD_CI_range <= 40)



# ----------------------------------
# Model fit summary
# ----------------------------------
# Summarize by year
Yearly_coverage <- Shifting %>% 
  group_by(Year) %>% 
  summarize(Year_count = n())
ggplot(Yearly_coverage, aes(Year, Year_count)) +
  geom_point() +
  scale_x_continuous(n.breaks = 10)

# Species year summary
Species_coverage <- Shifting %>% 
  group_by(Common_Name, Year) %>% 
  summarize(Common_Name_count = n())
mean(Species_coverage$Common_Name_count)
sd(Species_coverage$Common_Name_count)

# Number of days detected
mean(Shifting$JdDetect)
sd(Shifting$JdDetect)

# MAD CI range
hist(Shifting$MAD_CI_range)
quantile(Shifting$MAD_CI_range)
# below 10
MAD_CI_below_10 <- Shifting %>% 
  filter(MAD_CI_range < 10)
(nrow(MAD_CI_below_10)/nrow(Shifting))*100
mean(MAD_CI_below_10$JdDetect)
sd(MAD_CI_below_10$JdDetect)
summary(MAD_CI_below_10$JdDetect)
# Above 10
MAD_CI_above_10 <- Shifting %>% 
  filter(MAD_CI_range >= 10)
mean(MAD_CI_above_10$JdDetect)
sd(MAD_CI_above_10$JdDetect)
summary(MAD_CI_above_10$JdDetect)

# prsq
mean(Shifting$Prsq)
median(Shifting$Prsq)
sd(Shifting$Prsq)
summary(Shifting$Prsq)



# ---------------------------------
# Filter data
# ---------------------------------
# Calculate Grid-specific sampling coverage for species
Species_Grid_nyrs <- Shifting %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Shifting <- left_join(Shifting, Species_Grid_nyrs, by = c("Common_Name", "GridID"))
Shifting <- Shifting %>%
  filter(nyrs >= 10)
nrow(Shifting)
# # First of all, let's filter for species where sampling started at a minimum of 2007
# Shifting <- Shifting %>%
#   filter(min_Year <= 2007)

Species_Grid_Duration_Years <- Shifting %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Shifting <- left_join(Shifting, Species_Grid_Duration_Years, by = c("Common_Name", "GridID"))
# We will use these to further pare down the dataset
# Save the original data.frame
# Shifting_orig <- Shifting
# length(unique(Shifting_orig$GridID))
# 113
# # First of all, let's filter for species where sampling started at a minimum of 2007
# Shifting <- Shifting %>%
#   filter(min_Year <= 2007)

# Filter for species-grids that have a duration of at least 15 years
Shifting <- Shifting %>% 
  filter(Duration_Years >= 15)
nrow(Shifting)

# Remove species that have less than 7 grids (i.e., 1-6)
Species_Grid_n <- Shifting %>% 
  select(Common_Name, GridID) %>%
  distinct() %>% 
  group_by(Common_Name) %>% 
  count(name = "Spp_n_GridID") %>% 
  ungroup()
Shifting <- left_join(Shifting, Species_Grid_n, by = c("Common_Name"))
Shifting <- Shifting %>% 
  filter(Spp_n_GridID >= 6)
# Check # of rows
nrow(Shifting)

# Now, we want to see how many species there are across each of the grids
# Look at the maximum number of species that occur in each grid cell
# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Max_n <- Shifting %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(GridID) %>% 
  count(name = "GridID_n_spp") %>% 
  ungroup()
# Remove grids that only have one species
Shifting <- left_join(Shifting, Species_Grid_Max_n, by = c("GridID"))
Shifting <- Shifting %>% 
  filter(GridID_n_spp > 1)
# Check # of rows
nrow(Shifting)

# Check out the number of grids that a species has for a year
Species_Years <- Shifting %>% 
  group_by(Common_Name, Year) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 7 grids
Shifting <- left_join(Shifting, Species_Years, by = c("Common_Name", "Year"))
Shifting <- Shifting %>% 
  filter(ngrids >= 6)
nrow(Shifting)

# ALSO remove species that only have one grid
Grid_Species_Max_n <- Shifting %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(Common_Name) %>% 
  count(name = "Common_Name_n_GridID")
# Remove grids that only have one species
Shifting <- left_join(Shifting, Grid_Species_Max_n)
Shifting <- Shifting %>% 
  filter(Common_Name_n_GridID > 1)
# Check # of rows
nrow(Shifting)


write_csv(Shifting, "Shifting/Simulation_Shifting_ebird_MAD_Jd_80_180_w_enddate_5_Occ_Thresh_0.1_2002-2019_20220204B.csv")




# # ----------------------------------------------
# # Plot the scaled data
# # ----------------------------------------------
# Shifting_scaled <- Shifting %>% 
#   select(GridID, Year, Common_Name, MAD) %>% 
#   mutate(MAD_scaled = scale(MAD))
# ggplot(Shifting_scaled, aes(Year, MAD)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# ggplot(Shifting_scaled, aes(Year, MAD)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~Common_Name)
# 
# Stationary_scaled <- Stationary %>% 
#   select(GridID, Year, Common_Name, MAD) %>% 
#   mutate(MAD_scaled = scale(MAD))
# ggplot(Stationary_scaled, aes(Year, MAD)) +
#   geom_point() +
#   geom_smooth(method = "lm")
# ggplot(Stationary_scaled, aes(Year, MAD)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   facet_wrap(~Common_Name)
# 
# 
# -----------------------------------------------
# Calculate overall mean arrival date shift
# -----------------------------------------------
# Mixed model
testAllStationary <- lmer(MAD ~ Year + (1|GridID) + (1|Common_Name), data = Stationary)
summary(testAllStationary)

testAllShifting <- lmer(MAD ~ Year + (1|GridID) + (1|Common_Name), data = Shifting)
summary(testAllShifting)
# 
# 
# 
# # -----------------------------------------------
# # Calculate species shifts
# # -----------------------------------------------
# # Stationary
# # Get species-specific shifts
# coefficients_Stationary <- function (commonname){
#   Stationary <- Stationary[Stationary$Common_Name == commonname ,]
#   regression <- lmer(MAD ~ Year + (1|GridID), data = Stationary)
#   MAD.Shift <- data.frame(fixef(regression))
#   Common.name <- unique(Stationary$Common_Name)
#   data.frame(Common.name, MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# x_Stationary <- unique(Stationary$Common_Name)
# Bird_Shifts_Stationary <- lapply(x_Stationary, coefficients_Stationary)   
# Bird_Shifts_Stationary <- data.frame(do.call(rbind, Bird_Shifts_Stationary))
# colnames(Bird_Shifts_Stationary) <- c("Common_Name", "MAD_Shift")
# ggplot(Bird_Shifts_Stationary, aes(MAD_Shift)) +
#   geom_histogram(binwidth = 0.02, colour = "black") +
#   geom_vline(xintercept = median(Bird_Shifts_Stationary$MAD_Shift), linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "solid") +
#   labs(x = "Mean arrival date shift (days/year)",
#        y = "# of species-grid cells") +
#   theme_bw()
# summary(Bird_Shifts_Stationary$MAD_Shift)
# sd(Bird_Shifts_Stationary$MAD_Shift)
# 
# # Shifting
# # Get species-specific shifts
# coefficients_Shifting <- function (commonname){
#   Shifting <- Shifting[Shifting$Common_Name == commonname ,]
#   regression <- lmer(MAD ~ Year + (1|GridID), data = Shifting)
#   MAD.Shift <- data.frame(fixef(regression))
#   Common.name <- unique(Shifting$Common_Name)
#   data.frame(Common.name, MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# x_Shifting <- unique(Shifting$Common_Name)
# Bird_Shifts_Shifting <- lapply(x_Shifting, coefficients_Shifting)   
# Bird_Shifts_Shifting <- data.frame(do.call(rbind, Bird_Shifts_Shifting))
# colnames(Bird_Shifts_Shifting) <- c("Common_Name", "MAD_Shift")
# ggplot(Bird_Shifts_Shifting, aes(MAD_Shift)) +
#   geom_histogram(binwidth = 0.02, colour = "black") +
#   geom_vline(xintercept = median(Bird_Shifts_Shifting$MAD_Shift), linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "solid") +
#   labs(x = "Mean arrival date shift (days/year)",
#        y = "# of species-grid cells") +
#   theme_bw()
# summary(Bird_Shifts_Shifting$MAD_Shift)
# sd(Bird_Shifts_Shifting$MAD_Shift)
# 
# 
# 
# # -----------------------------------------------
# # Calculate grid shifts
# # -----------------------------------------------
# # Stationary
# # Get grid-specific shifts
# coefficients_Stationary_2 <- function (gridid){
#   Stationary <- Stationary[Stationary$GridID == gridid ,]
#   regression <- lmer(MAD ~ Year + (1|Common_Name), data = Stationary)
#   MAD.Shift <- data.frame(fixef(regression))
#   Grid.ID <- unique(Stationary$Common_Name)
#   data.frame(Grid.ID, MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# x_Stationary_2 <- unique(Stationary$GridID)
# Bird_Shifts_Stationary_2 <- lapply(x_Stationary_2, coefficients_Stationary_2)   
# Bird_Shifts_Stationary_2 <- data.frame(do.call(rbind, Bird_Shifts_Stationary_2))
# colnames(Bird_Shifts_Stationary_2) <- c("GridID", "MAD_Shift")
# ggplot(Bird_Shifts_Stationary_2, aes(MAD_Shift)) +
#   geom_histogram(binwidth = 0.02, colour = "black") +
#   geom_vline(xintercept = median(Bird_Shifts_Stationary$MAD_Shift), linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "solid") +
#   labs(x = "Mean arrival date shift (days/year)",
#        y = "# of species-grid cells") +
#   theme_bw()
# summary(Bird_Shifts_Stationary_2$MAD_Shift)
# sd(Bird_Shifts_Stationary_2$MAD_Shift)
# 
# # Shifting
# # Get grid-specific shifts
# coefficients_Shifting_2 <- function (gridid){
#   Shifting <- Shifting[Shifting$GridID == gridid ,]
#   regression <- lmer(MAD ~ Year + (1|Common_Name), data = Shifting)
#   MAD.Shift <- data.frame(fixef(regression))
#   Grid.ID <- unique(Shifting$Common_Name)
#   data.frame(Grid.ID, MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# x_Shifting_2 <- unique(Shifting$GridID)
# Bird_Shifts_Shifting_2 <- lapply(x_Shifting_2, coefficients_Shifting_2)   
# Bird_Shifts_Shifting_2 <- data.frame(do.call(rbind, Bird_Shifts_Shifting_2))
# colnames(Bird_Shifts_Shifting_2) <- c("GridID", "MAD_Shift")
# ggplot(Bird_Shifts_Shifting_2, aes(MAD_Shift)) +
#   geom_histogram(binwidth = 0.02, colour = "black") +
#   geom_vline(xintercept = median(Bird_Shifts_Shifting$MAD_Shift), linetype = "dashed") +
#   geom_vline(xintercept = 0, linetype = "solid") +
#   labs(x = "Mean arrival date shift (days/year)",
#        y = "# of species-grid cells") +
#   theme_bw()
# summary(Bird_Shifts_Shifting_2$MAD_Shift)
# sd(Bird_Shifts_Shifting_2$MAD_Shift)
# 
# 
# 
# 
# # 
# # # -------------------------------------------------
# # # Mixed model by migratory distance
# # # -------------------------------------------------
# # #Define Migratory distance for each species
# # Stationary$Mig.Dist <- ifelse(Stationary$Common_Name == "Palm Warbler" |
# #                                     Stationary$Common_Name == "Common Yellowthroat" |
# #                                     Stationary$Common_Name == "Baltimore Oriole" |
# #                                     Stationary$Common_Name == "Black-and-white Warbler","Short","Long")
# # 
# # #Linear mixed effects model with Year as fixed effect and Grids as random effects separated by migration category
# # LMM_long <- Stationary[Stationary$Mig.Dist == "Long",]
# # testinglong <- lmer(MAD ~ Year + (1|GridID), data = LMM_long)
# # summary(testinglong)
# # LMM_short <- Stationary[Stationary$Mig.Dist == "Short",]
# # testingshort <- lmer(MAD ~ Year + (1|GridID), data = LMM_short)
# # summary(testingshort)
# 
# 
# 
# # # -------------------------------------------------
# # # Morphology analysis
# # # -------------------------------------------------
# # #Read in morphology parameters
# # para <- read.csv(file="parameters.csv", header = TRUE)
# # para <- para %>%
# #   rename("Common_Name" = "Common.name")
# # param <- merge(Bird_Shifts, para, by = "Common_Name", all = TRUE)
# # param <- na.omit(param)
# # #Read phylogenetic tree for our species list
# # tree <- read.nexus(file="tree.nex")
# # #check tree
# # plot(tree, cex=0.7)
# # 
# # #create object which contains comparative data (morhpologies and mad shift) with tree
# # #warning will show if any tips on the tree are dropped due to missing data (this should happen, no data for Chestnut Sided warbler)
# # cdata<-comparative.data(tree, param, names.col = Tree.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)
# # 
# # #generate a model between mass and WLI, use residuals to correct for body mass
# # correction <- pgls(log(WLI) ~ Body.mass..log., data = cdata, lambda = "ML")
# # Mass.adj.WLI <- data.frame(residuals(correction))
# # #merge mass adjusted WLI to original parameter dataframe
# # row.names(param) <- param$Tree.name
# # param <- merge(param, Mass.adj.WLI, by="row.names")
# # param <- param[-c(1)]
# # names(param)[10]<-"Mass.adj.WLI"
# # 
# # cdata <- comparative.data(tree, param, names.col = Tree.name, vcv = TRUE, na.omit = TRUE, warn.dropped = TRUE)
# # 
# # {mod1<- pgls(MAD_Shift ~ log(WAR), data = cdata, lambda="ML")
# #   mod2<- pgls(MAD_Shift ~ log(HWI), data = cdata, lambda="ML")
# #   mod3<- pgls(MAD_Shift ~ Mass.adj.WLI, data = cdata, lambda="ML")
# #   mod4<- pgls(MAD_Shift ~ Body.mass..log., data = cdata, lambda="ML")
# # }
# # 
# # summary(mod1)
# # summary(mod2)
# # summary(mod3)
# # summary(mod4)
# 
# 
# 
# # -----------------------------------------------
# # Calculate mean arrival date shifts by Grid
# # -----------------------------------------------
# # Get grid-specific shifts
# Grid_coefficients <- function (gridid){
#   Stationary <- Stationary[Stationary$GridID == gridid ,]
#   regression <- lmer(MAD ~ Year + (1|Common_Name), data = Stationary)
#   MAD.Shift <- data.frame(fixef(regression))
#   grid.id <- unique(Stationary$GridID)
#   data.frame(grid.id, MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# Grids <- unique(Stationary$GridID)
# Grid_Shifts_lmer <- lapply(Grids, Grid_coefficients)
# Grid_Shifts_lmer <- data.frame(do.call(rbind, Grid_Shifts_lmer))
# colnames(Grid_Shifts_lmer) <- c("GridID", "MAD_Shift")
# 
# 
# 
# 
# # ---------------------------------------------------------------
# # Identify species-Grid combinations for passing vs. staying
# # ---------------------------------------------------------------
# passage_grid_data <- Stationary %>%
#   group_by(Common_Name, GridID) %>%
#   summarize(mean_MaxJd = mean(MaxJd),
#             median_MaxJd = median(MaxJd),
#             max_MaxJd = max(MaxJd),
#             min_MaxJd = min(MaxJd)) %>%
#   ungroup()
# # Plot the data by species (median)
# ggplot(passage_grid_data, aes(median_MaxJd), colour = "black") +
#   geom_histogram(binwidth = 1) +
#   theme_bw() +
#   scale_x_continuous(n.breaks = 6) +
#   facet_wrap(~Common_Name)
# # Plot the data by species (mean)
# ggplot(passage_grid_data, aes(mean_MaxJd), colour = "black") +
#   geom_histogram(binwidth = 1) +
#   theme_bw() +
#   scale_x_continuous(n.breaks = 6) +
#   facet_wrap(~Common_Name)
# 
# # Try using the median MaxJd as an indicator, and set the level to 180
# # Create a new column representing this
# passage_grid_data <- passage_grid_data %>%
#   mutate(passing_grid = if_else(median_MaxJd < 180, TRUE, FALSE))
# # How many passing grids are there per species?
# passage_grid_data_count <- passage_grid_data %>%
#   group_by(Common_Name, passing_grid) %>%
#   count(name = "n_passage_grids_per_spp") %>%
#   pivot_wider(names_from = passing_grid, names_prefix = "n_passage_grids_per_spp_", values_from = n_passage_grids_per_spp)
# 
# # Plot the summary of grids in each category for each species
# ggplot(passage_grid_data, aes(passing_grid)) +
#   geom_bar(colour = "black") +
#   theme_bw() +
#   facet_wrap(~Common_Name)
# 
# # Join the passage grid data to the main data.frame
# Stationary <- left_join(Stationary, select(passage_grid_data, Common_Name, GridID, median_MaxJd, passing_grid)) %>% left_join(passage_grid_data_count)
# # Filter the dataset
# # Filter for if the passing grids are true, and if there are at least 5 passing grids (or staying grids) that are true for that species.
# Stationary_passing_TRUE <- Stationary %>%
#   filter(n_passage_grids_per_spp_TRUE >= 5 & passing_grid == TRUE)
# Stationary_passing_FALSE <- Stationary %>%
#   filter(n_passage_grids_per_spp_FALSE >= 5 & passing_grid == FALSE)
# 
# 
# # Get passing-grid-specific shifts
# # Create functions
# coefficients_PT <- function (commonname){
#   Stationary_passing_TRUE <- Stationary_passing_TRUE[Stationary_passing_TRUE$Common_Name == commonname ,]
#   regression <- lmer(MAD ~ Year + (1|GridID), data = Stationary_passing_TRUE)
#   MAD.Shift <- data.frame(fixef(regression))
#   Common.name <- unique(Stationary_passing_TRUE$Common_Name)
#   data.frame(Common.name,MAD.Shift[2,])
# }
# coefficients_PF <- function (commonname){
#   Stationary_passing_FALSE <- Stationary_passing_FALSE[Stationary_passing_FALSE$Common_Name == commonname ,]
#   regression <- lmer(MAD~Year + (1|GridID), data = Stationary_passing_FALSE)
#   MAD.Shift <- data.frame(fixef(regression))
#   Common.name <- unique(Stationary_passing_FALSE$Common_Name)
#   data.frame(Common.name,MAD.Shift[2,])
# }
# 
# #Create dataframe of slopes == MAD Shifts
# # Passing == TRUE (PT)
# xPT <- unique(Stationary_passing_TRUE$Common_Name)
# Bird_Shifts_PT <- lapply(xPT, coefficients_PT)
# Bird_Shifts_PT <- data.frame(do.call(rbind, Bird_Shifts_PT))
# colnames(Bird_Shifts_PT) <- c("Common_Name", "Passing grids")
# # Passing == FALSE (PF)
# xPF <- unique(Stationary_passing_FALSE$Common_Name)
# Bird_Shifts_PF <- lapply(xPF, coefficients_PF)
# Bird_Shifts_PF <- data.frame(do.call(rbind, Bird_Shifts_PF))
# colnames(Bird_Shifts_PF) <- c("Common_Name", "Staying grids")
# 
# # Join together
# Bird_Shifts_PTPF <- left_join(Bird_Shifts, Bird_Shifts_PF) %>%
#   left_join(Bird_Shifts_PT)
# Bird_Shifts_PTPF <- Bird_Shifts_PTPF %>%
#   rename("All grids" = "MAD_Shift")
# # Plot the correlations
# 
# 
# lowerfun <- function(data,mapping){
#   ggplot(data = data, mapping = mapping)+
#     geom_abline(linetype = "dashed") +
#     geom_smooth(method = "lm", colour = "black", size = 0.6) +
#     geom_point(alpha = 0.5)
# }
# 
# 
# ggpairs(Bird_Shifts_PTPF[2:4],
#         lower = list(continuous = wrap(lowerfun)))  +
#   theme_bw() +
#   labs(x = "Mean arrival date shift (days/year)",
#        y = "Mean arrival date shift (days/year)")
# 
# ggsave(filename = "C:/Users/dpgil/Documents/1 U of T/eBird/Figures/Supplementary analysis/Passing vs Staying Grids/Passing_vs_Staying_vs_all_grids.png", units = "cm", dpi = 300, height = 12, width = 12)
# 
# Bird_Shifts_PTPF_noNA <- na.omit(Bird_Shifts_PTPF)
# nrow(Bird_Shifts_PTPF_noNA)
# Bird_Shifts_PTPF_noNA$`Staying grids` < Bird_Shifts_PTPF_noNA$`Passing grids`
# mean(Bird_Shifts_PTPF_noNA$`Staying grids` )
# mean(Bird_Shifts_PTPF_noNA$`Passing grids` )
# 
# # Plot
# ggplot(Bird_Shifts_PTPF_noNA, aes(`Staying grids`, `Passing grids`)) +
#   geom_abline(linetype = "dashed") +
#   geom_smooth(method = "lm", colour = "black", size = 0.6) +
#   geom_point(aes(colour = Common_Name), alpha = 0.5) +
#   theme_bw()
