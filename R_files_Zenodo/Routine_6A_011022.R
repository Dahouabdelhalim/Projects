# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 6A
# ------------------------------------------------
# # Load Libraries
library(tidyverse)
library(vroom)

# Set Working Directory
setwd("~/Data repository")

# List files
files <- list.files(path = "Routine 6", full.names = T)
files <- files[1:159]

# Routine 6 - merged files
merged_data <- vroom(files)

# Check # of grids
Grids <- merged_data %>% 
  select(GridID) %>% 
  distinct()

# Which grids are not included?
Grid_files <- data.frame(GridID = files)
Grid_files$GridID <- str_remove_all(Grid_files$GridID, "_site_occupancy.csv")
absent_grids <- anti_join(Grid_files, Grids)

# Inspect these files
absent_grids$GridID <- paste(absent_grids$GridID, "_site_occupancy.csv", sep = "")
# Read in the first file
absent_grids_list <- list()
for (i in 1:nrow(absent_grids)) {
  absent_grids_list[[i]] <- read_csv(absent_grids$GridID[1])
  
}
absent_grids_list[[1]]
absent_grids_list[[2]]
absent_grids_list[[3]]
absent_grids_list[[4]]
absent_grids_list[[5]]
absent_grids_list[[6]]
absent_grids_list[[7]]
absent_grids_list[[8]]
absent_grids_list[[9]]
absent_grids_list[[10]]
absent_grids_list[[11]]
absent_grids_list[[12]]
absent_grids_list[[13]]

# Write the file
write_csv(merged_data, "Routine 6/Routine 6 - merged files/eBird_Daily_Occupancy_all_grids.csv")
