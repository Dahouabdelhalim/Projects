# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 4
# ------------------------------------------------

# Load Libraries
library(dplyr)
library(stringr)
library(readr)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Start function subroutines
# ------------------------------------------------

# ------------------------------------------------
# Finish function subroutines
# ------------------------------------------------

# ------------------------------------------------
# Read in datasets
# ------------------------------------------------
# List filenames of all regional grid subsets
Routine_3_filesnames <- list.files(path = "Routine 3", full.names = F)
Routine_3_gridnames <- Routine_3_filesnames
Routine_3_gridnames <- data.frame(GridID = Routine_3_gridnames)
# Remove regions and .csv from gridnames dataframe
Routine_3_gridnames <- as.data.frame(str_remove_all(Routine_3_gridnames$GridID, "(_{1}[A-z]+[.][A-z]+)"))
colnames(Routine_3_gridnames) <- "GridID"
# Get distinct grids
Routine_3_gridnames_u <- distinct(Routine_3_gridnames)
# Number of unique grids
nrow(Routine_3_gridnames_u)
# Read in all regional grid datasets for each unique grid and merge them together
for (i in 1:nrow(Routine_3_gridnames_u)) {
  # i = 1
  Regional_grid_filenames <- list.files(path = "Routine 3", pattern = paste("^", Routine_3_gridnames_u[i,], sep = ""))
  Regional_grid_list <- list()
  for (j in 1:length(Regional_grid_filenames)) {
  # j = 1
    Regional_grid_list[[j]] <- read_csv(paste("Routine 3/", Regional_grid_filenames[j], sep = ""))
  }
    # Create names for the data list items
    names(Regional_grid_list) <- Regional_grid_filenames
    # Bind the dataframes
    Regional_grid_data <- bind_rows(Regional_grid_list)
    # Write the dataframe
    write_csv(Regional_grid_data, paste("Routine 4/", Routine_3_gridnames_u[i,], ".csv", sep = ""))
    # Remove the dataframe
    rm(Regional_grid_data, Regional_grid_list)
}




