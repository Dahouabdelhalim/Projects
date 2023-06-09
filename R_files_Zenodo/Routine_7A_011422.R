# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 7A
# ------------------------------------------------
# Load Libraries
library(dplyr)
library(stringr)
library(readr)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Read in processed data
# ------------------------------------------------
# # Use this if you would like to re-create MAD_Output_All
# Routine_7_filesnames <- data.frame(file = list.files(path = "Routine 7/Processing by species"))
# Routine_7_filesnames <- Routine_7_filesnames %>%
#   filter(str_detect(.$file, ".csv"))
# Routine_7_files_list <- list()
# for (i in 1:nrow(Routine_7_filesnames)) {
# # i = 1
#   Routine_7_files_list[[i]] <- read_csv(file = paste("Routine 7/Processing by species/", Routine_7_filesnames$file[i], sep = ""))
#   # Extract Jd_min from each file, as it wasn't obtained earlier
#   Routine_7_files_list[[i]]$Jd_min <- str_extract(Routine_7_filesnames$file[i], "(Jd_)[0-9]*")
#   Routine_7_files_list[[i]]$Jd_min <- as.numeric(str_remove_all(Routine_7_files_list[[i]]$Jd_min, "(Jd_)"))
# }
# Routine_7_files <- bind_rows(Routine_7_files_list)
# # Write the master, binded csv
# unique(Routine_7_files$Jd_min)
# 
# 
# write_csv(Routine_7_files, "Routine 7/Merged output/MAD_Output_All.csv")

# Start here if the previous chunk was completed in another session
Routine_7_files <- read_csv("Routine 7/Merged output/MAD_Output_All.csv")
spp_orig <- distinct(select(Routine_7_files, Common_Name))
write_csv(spp_orig, "Routine 7/Merged output/spp_orig.csv")



# ------------------------------------------------
# Remove some species based on geography and migration strategy
# ------------------------------------------------
Routine_7_files <- Routine_7_files %>%
  filter(!Common_Name %in% c("Eastern Meadowlark", "Marsh Wren"))



# ------------------------------------------------
# Filter data by minimum standards
# ------------------------------------------------
# Blocks with less than 10 days of data
Routine_7_files <- Routine_7_files %>% 
  filter(JdDetect >= 10) 
nrow(Routine_7_files)
# 364098
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# Blocks with less than or equal to CI 40
Routine_7_files <- Routine_7_files %>% 
  filter(MAD_CI_range <= 40)
nrow(Routine_7_files)
# 339515
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# Calculate Grid-specific sampling coverage for species

# Filter for species-grids that have number of years of at least 10 years
Species_Grid_nyrs <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(nyrs = n()) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_nyrs, by = c("Common_Name", "GridID", "Jd_min"))
Routine_7_files <- Routine_7_files %>%
  filter(nyrs >= 10)
nrow(Routine_7_files)
# 199368
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 1 Brewer's Blackbird
# 2 Eastern Meadowlark
# 3 Marsh Wren 

# Remember to group for Jd_min, because you haven't decided on that yet
Species_Grid_Duration_Years <- Routine_7_files %>% 
  group_by(Common_Name, GridID, Jd_min) %>% 
  summarize(min_Year = min(Year),
            max_Year = max(Year),
            Duration_Years = max_Year - min_Year + 1) %>% 
  ungroup()
# Join these summary statistics to the data.frame
Routine_7_files <- left_join(Routine_7_files, Species_Grid_Duration_Years, by = c("Common_Name", "GridID", "Jd_min"))
# We will use these to further pare down the dataset

# Filter for species-grids that have a duration of at least 15 years
Routine_7_files <- Routine_7_files %>% 
  filter(Duration_Years >= 15)
nrow(Routine_7_files)
# 170614
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 1 Bell's Vireo      
# 2 Brewer's Blackbird
# 3 Eastern Meadowlark
# 4 Lark Sparrow      
# 5 Marsh Wren        
# 6 Nelson's Sparrow  
# 7 Vesper Sparrow 

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
# 162258
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name))) %>% print(n = 32)
# 30 species removed

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
# 162042
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 30 species removed

# Check out the number of grids that a species has for a year
Species_Years <- Routine_7_files %>% 
  group_by(Common_Name, Year, Jd_min) %>% 
  summarize(ngrids = n()) %>% 
  ungroup()
# Remove years that have less than 6 grids
Routine_7_files <- left_join(Routine_7_files, Species_Years, by = c("Common_Name", "Year", "Jd_min" ))
Routine_7_files <- Routine_7_files %>% 
  filter(ngrids >= 6)
# Check # of rows
nrow(Routine_7_files)
# 156805
anti_join(spp_orig, distinct(select(Routine_7_files, Common_Name)))
# 30 species removed

# write the progress
write_csv(Routine_7_files, "Routine 7/Merged output/MAD_Output_All_A.csv")

# ----------------------------------
# Minimum JD summary
# ----------------------------------
min_JD_summary <- Routine_7_files %>% 
  group_by(Common_Name, Jd_min) %>% 
  summarize(Mean_MAD = mean(MAD)) %>% 
  ungroup()
min_JD_summary %>% 
  count(Jd_min) %>% 
  arrange(desc(n))
write_csv(min_JD_summary, "Routine 7/Merged output/min_JD_summary.csv")
