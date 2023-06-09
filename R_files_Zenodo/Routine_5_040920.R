# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 5
# ------------------------------------------------
# Set scientific notation off
options(scipen = 999)

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
# Read in species names
# ------------------------------------------------
spp <- read_csv("MorphologyParameters.csv")
spp <- dplyr::select(spp, Common.name)

# ------------------------------------------------
# Read in grid datasets
# ------------------------------------------------
Routine_4_filesnames <- list.files(path = "Routine 4", full.names = T)
Grid_data_rows_list <- list()
for (i in 1:length(Routine_4_filesnames)) {
  Grid_data <- read_csv(Routine_4_filesnames[i])
  Grid_data_rows_list[[i]] <- data.frame(Rows = nrow(Grid_data), X = str_extract(Grid_data$GridID[i], "[0-9]+$"), Y = str_extract(Grid_data$GridID[i], "([-]*[0-9]+[*])"))
  Grid_data_rows_list[[i]]$Y <- as.numeric(str_remove_all(Grid_data_rows_list[[i]]$Y, "[*]"))
  rm(Grid_data)
}


for (i in 1:length(Routine_4_filesnames)) {
  Grid_data <- read_csv(Routine_4_filesnames[i])
  Grid_data_sum_1 <- Grid_data %>% 
    group_by(Year) %>% 
    summarize(Obs_per_yr = n()) %>% 
    ungroup()
  Grid_data_sum_2 <- Grid_data %>% 
    filter(`COMMON NAME` %in% spp$Common.name) %>% 
    group_by(Year) %>% 
    summarize(Obs_per_yr_focal_spp = n()) %>% 
    ungroup()
  Grid_data_sum_3 <- Grid_data %>% 
    filter(`COMMON NAME` %in% spp$Common.name) %>%
    dplyr::select(`COMMON NAME`, Year) %>%
    group_by(Year) %>% 
    distinct() %>% 
    summarize(Focal_spp_per_yr = n()) %>% 
    ungroup()
  Grid_data_sum_1 <- left_join(Grid_data_sum_1, Grid_data_sum_2)
  Grid_data_sum_1 <- left_join(Grid_data_sum_1, Grid_data_sum_3)
  Grid_data_sum_1 <- Grid_data_sum_1 %>% 
    mutate(Rows = nrow(Grid_data),
           X = str_extract(Grid_data$GridID[1], "[0-9]+$"),
           Y = str_extract(Grid_data$GridID[1], "([-]*[0-9]+[*])"),
           GridID = Grid_data$GridID[1])
  Grid_data_sum_1$X <- as.numeric(Grid_data_sum_1$X)
  Grid_data_sum_1$Y <- as.numeric(str_remove_all(Grid_data_sum_1$Y, "[*]"))
  Grid_data_rows_list[[i]] <- Grid_data_sum_1
  rm(Grid_data)
}
# Bind the dataframes
Grid_data_summary <- bind_rows(Grid_data_rows_list)
# Replace NAs with 0s
Grid_data_summary[is.na(Grid_data_summary)] <- 0
# Get sum of focal species found in grid cell over the years
Grid_data_summary <- Grid_data_summary %>% 
  group_by(GridID) %>% 
  mutate(Focal_spp_sum = sum(Focal_spp_per_yr)) %>% 
  ungroup()

# Create reduced, unique GridID dataframe 
Grid_data_summary_u <- Grid_data_summary %>% 
  dplyr::select(Rows, X, Y, GridID, Focal_spp_sum) %>% 
  distinct()

# Write dataframes
write_csv(Grid_data_summary, "Routine 5/Grid_data_summary.csv")
write_csv(Grid_data_summary_u, "Routine 5/Grid_data_summary_unique.csv")
