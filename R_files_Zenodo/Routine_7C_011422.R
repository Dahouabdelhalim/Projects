# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7C
# ------------------------------------------------
# Load Libraries
library(dplyr)
library(readr)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Read in processed data
# ------------------------------------------------
Routine_7_files <- read_csv("Routine 7/Merged output/MAD_Output_All_B.csv")
spp_orig <- read_csv("Routine 7/Merged output/spp_orig.csv")

# -------------------------------------------------------
# Re-visit each filter (except JdDetect and CI)
# -------------------------------------------------------
# Filter for species-grids that have number of years of at least 10 years
Routine_7_files <- Routine_7_files %>% 
  select(-nyrs, -min_Year, -max_Year, -Duration_Years, -Spp_n_GridID, -GridID_n_spp, -ngrids)
Species_Grid_nyrs <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_nyrs, by = c("Common_Name", "GridID", "Jd_min"))
Routine_7_files <- Routine_7_files %>%
  filter(nyrs >= 10)
nrow(Routine_7_files)
# 146182
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 41 species removed

# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Duration_Years <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_Duration_Years, by = c("Common_Name", "GridID", "Jd_min"))

# Filter for species-grids that have a duration of at least 15 years
Routine_7_files <- Routine_7_files %>% 
  filter(Duration_Years >= 15)
nrow(Routine_7_files)
#  146142
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 41 species removed

# Remove species that have less than 6 grids (i.e., 1-5)
Species_Grid_n <- Routine_7_files %>% 
  select(Common_Name, GridID, Jd_min) %>%
  distinct() %>% 
  group_by(Common_Name, Jd_min) %>% 
  count(name = "Spp_n_GridID") %>% 
  ungroup()
Routine_7_files <- left_join(Routine_7_files, Species_Grid_n, by = c("Common_Name", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(Spp_n_GridID >= 6)
# Check # of rows
nrow(Routine_7_files)
# 145842
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name))) %>% print(n = 32)
# 43 species removed

# Now, we want to see how many species there are across each of the grids
# Look at the maximum number of species that occur in each grid cell
# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Max_n <- Routine_7_files %>% 
  select(Common_Name, GridID, Jd_min) %>% 
  distinct() %>% 
  group_by(GridID, Jd_min) %>% 
  count(name = "GridID_n_spp") %>% 
  ungroup()
# Remove grids that only have one species
Routine_7_files <- left_join(Routine_7_files, Species_Grid_Max_n, by = c("GridID", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(GridID_n_spp > 1)
# Check # of rows
nrow(Routine_7_files)
# 145842
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 43 species removed

# Check out the number of grids that a species has for a year
Species_Years <- Routine_7_files %>% 
  group_by(Common_Name, Year, Jd_min) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 6 grids
Routine_7_files <- left_join(Routine_7_files, Species_Years, by = c("Common_Name", "Year", "Jd_min"))
Routine_7_files <- Routine_7_files %>% 
  filter(ngrids >= 6)
nrow(Routine_7_files)
# 145601

# write the progress
write_csv(Routine_7_files, "Routine 7/Merged output/MAD_Output_All_C.csv")
