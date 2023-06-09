# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7H
# ------------------------------------------------
# Load Libraries
library(dplyr)
library(readr)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Read in species names
# ------------------------------------------------
spp <- read_csv("MorphologyParameters.csv")
spp_names <- spp$Common.name

spp_orig <- read_csv("Routine 7/Merged output/spp_orig.csv")

Routine_7_files_top <- read_csv("Routine 7/Merged Output/MAD_Output_All_Filter_pre.csv")


# Double check all of the filters
Routine_7_files_top_orig <- nrow(Routine_7_files_top)
# Blocks with less than 10 days of data
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(JdDetect >= 10) 
Routine_7_files_top_orig == nrow(Routine_7_files_top)
# Blocks with less than or equal to CI 40
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(MAD_CI_range <= 40)
Routine_7_files_top_orig == nrow(Routine_7_files_top)
# Calculate Grid-specific sampling coverage for species
Species_Grid_nyrs <- Routine_7_files_top %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files_top <- Routine_7_files_top %>%
  select(-min_Year, -max_Year, -Duration_Years, -nyrs)
Routine_7_files_top <- left_join(Routine_7_files_top, Species_Grid_nyrs)

# Filter for species-grids that have number of years of at least 10 years
Routine_7_files_top <- Routine_7_files_top %>%
  filter(nyrs >= 10)
Routine_7_files_top_orig == nrow(Routine_7_files_top)

# Filter for species-grids that have a duration of at least 15 years
Species_Grid_Duration_Years <- Routine_7_files_top %>% 
  group_by(Common_Name, GridID) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files_top <- left_join(Routine_7_files_top, Species_Grid_Duration_Years, by = c("Common_Name", "GridID"))
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(Duration_Years >= 15)
Routine_7_files_top_orig == nrow(Routine_7_files_top)

# Remove species that have less than 6 grids (i.e., 1-5)
Species_Grid_n_top <- Routine_7_files_top %>% 
  select(Common_Name, GridID) %>%
  distinct() %>% 
  group_by(Common_Name) %>% 
  count(name = "Spp_n_GridID") %>% 
  ungroup()
Routine_7_files_top <- Routine_7_files_top %>% 
  select(-Spp_n_GridID)
Routine_7_files_top <- left_join(Routine_7_files_top, Species_Grid_n_top)
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(Spp_n_GridID >= 6)
# Check # of rows
Routine_7_files_top_orig == nrow(Routine_7_files_top)

# Now, we want to see how many species there are across each of the grids
# Look at the maximum number of species that occur in each grid cell
Species_Grid_Max_n_top <- Routine_7_files_top %>% 
  select(Common_Name, GridID) %>% 
  distinct() %>% 
  group_by(GridID) %>% 
  count(name = "GridID_n_spp")
Species_Grid_Max_n_top %>% 
  filter(GridID_n_spp < 2)
# Remove grids that only have one species
Routine_7_files_top <- Routine_7_files_top %>% 
  select(-GridID_n_spp)
Routine_7_files_top <- left_join(Routine_7_files_top, Species_Grid_Max_n_top)
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(GridID_n_spp > 1)
Routine_7_files_top_orig == nrow(Routine_7_files_top)
# These are different because if they had different optimal Jd_min values, and these grids only showed up for the ones that were not optimal, then they're gone
# anti_df <- anti_join(select(Routine_7_files_top_orig, Common_Name, Year, GridID), select(Routine_7_files_top, Common_Name, Year, GridID))

# Check out the number of grids that a species has for a year
Species_Years <- Routine_7_files_top %>% 
  group_by(Common_Name, Year) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 7 grids
Routine_7_files_top <- Routine_7_files_top %>% 
  select(-ngrids)
Routine_7_files_top <- left_join(Routine_7_files_top, Species_Years, by = c("Common_Name", "Year"))
Routine_7_files_top <- Routine_7_files_top %>% 
  filter(ngrids >= 6)
nrow(Routine_7_files_top)
Routine_7_files_top_orig == nrow(Routine_7_files_top)

write_csv(Routine_7_files_top, "Routine 7/Merged Output/MAD_Output_All_Filter_1.csv")
