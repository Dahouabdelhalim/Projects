# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 6
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
# Read in species names
# ------------------------------------------------
spp <- read_csv("MorphologyParameters.csv")
spp <- dplyr::select(spp, Common.name)

# ------------------------------------------------
# Read in grid cell datasets
# ------------------------------------------------
Routine_4_filesnames <- list.files(path = "Routine 4", full.names = F)


for (i in 1:length(Routine_4_filesnames)) {
  # i = 1
  # **** NOTE, the code will throw errors at the numbers below, looping will stop when these numbers are encountered. Can consider doing it one grid at a time***
  # for (i in c(20, 21, 126)) {
# Skipped 126
  # Couldn't do 36, 37, 38, 52, 103, 106, 109, 110, 129, 141, 142
# i = 26
Grid_dataset <- read_csv(paste("Routine 4/", Routine_4_filesnames[i], sep = ""))

# Change column names
colnames(Grid_dataset) <- c("Common_Name",
                          "Locality_ID",
                          "Year",
                          "Jd",
                          "Region",
                          "X",
                          "Y",
                          "Dist_m",
                          "GridID")

# Use distinct_at() function (omitting region) as this analysis does not focus on the number of birds at each site
Grid_dataset <- Grid_dataset %>% 
  distinct_at(.vars = vars(Common_Name,
           Locality_ID,
           Year,
           Jd
           # ,
           # Region,
           # X,
           # Y,
           # Dist_m,
           # GridID
           )
           # ,
           # .keep_all = TRUE
           )

# Get the unique site checklists for each date
# Keep all species of bird within this dataframe as we want all potential sites, including those that might be missing the birds of interest
Grid_dataset_sites <- Grid_dataset %>% 
  distinct_at(.vars = vars(Locality_ID,
                Year,
                Jd
                # ,
                # Region,
                # X,
                # Y,
                # Dist_m,
                # GridID
                ))


# Add all of the species to the date-specific unique site checklists data
# Filter for species of interest (as not all species may be within the region of interest)
Grid_dataset_spp <- Grid_dataset %>%
  distinct_at(vars(Common_Name)) %>% 
  filter(Common_Name %in% spp$Common.name)
Grid_dataset_all_combos <- merge(Grid_dataset_sites,
                          Grid_dataset_spp)

# Join the species presence data to the unique site with species data
# Filter the main dataframe to only include the focal species, as we don't need to calculate occupancy data for species we are not interested in
Grid_dataset <- Grid_dataset %>% 
  filter(Common_Name %in% spp$Common.name)
Grid_dataset$Present <- 1
Grid_dataset <- full_join(Grid_dataset_all_combos,
                                   Grid_dataset)
rm(Grid_dataset_all_combos)

# replace NAs with zeros
Grid_dataset[is.na(Grid_dataset)] <- 0

Grid_dataset_site_count <- Grid_dataset_sites %>% 
  group_by(Year,
           Jd) %>% 
  tally(name = "Daily_Site_Total") %>% 
  ungroup()
# Get yearly site total
Grid_dataset_site_count_yearly <- Grid_dataset_site_count %>% 
  group_by(Year) %>% 
  summarize(Total_Yearly_Sites = sum(Daily_Site_Total)) %>% 
  ungroup()
rm(Grid_dataset_sites)

Grid_dataset_tally <- Grid_dataset %>% 
  group_by(Common_Name,
           Year,
           Jd) %>% 
  summarize(Daily_Site_Sightings = sum(Present)) %>% 
  ungroup()

# Join the number of unique sites for a given date to the daily occupancy dataframe
Grid_dataset <- left_join(Grid_dataset_tally,
                                  Grid_dataset_site_count)
rm(Grid_dataset_site_count, Grid_dataset_tally)
Grid_dataset <- Grid_dataset %>% 
  mutate(Daily_Site_Occupancy = Daily_Site_Sightings/Daily_Site_Total)

# Look at the how many times one or more of the focal birds were observed at a site over the study period
Grid_dataset_common_birds <- Grid_dataset %>% 
  filter(Common_Name %in% spp$Common.name) %>% 
  group_by(Common_Name) %>% 
  summarise(Total_Sightings = sum(Daily_Site_Sightings)) %>% 
  arrange(desc(Total_Sightings)) %>% 
  ungroup()
# Look at the how many times one or more of the focal birds were observed at a site over the study period, by year
Grid_dataset_common_birds_yearly <- Grid_dataset %>% 
  filter(Common_Name %in% spp$Common.name) %>% 
  group_by(Common_Name, Year) %>% 
  summarise(Total_Yearly_Sightings = sum(Daily_Site_Sightings)) %>% 
  arrange(desc(Total_Yearly_Sightings)) %>% 
  ungroup()

# Write the files
# Get the grid name
GridID <- Routine_4_filesnames[i]
GridID <- str_remove_all(GridID, ".csv")
Grid_dataset$GridID <- GridID
Grid_dataset_common_birds$GridID <- GridID
Grid_dataset_common_birds_yearly$GridID <- GridID
Grid_dataset_site_count_yearly$GridID <- GridID
write_csv(Grid_dataset, paste("Routine 6/", GridID, "_site_occupancy.csv", sep = ""))
write_csv(Grid_dataset_common_birds, paste("Routine 6/Site and sightings data/Common birds/", GridID, "_common_birds.csv", sep = ""))
write_csv(Grid_dataset_common_birds_yearly, paste("Routine 6/Site and sightings data/Common birds yearly/",GridID, "_common_birds_yearly.csv", sep = ""))
write_csv(Grid_dataset_site_count_yearly, paste("Routine 6/Site and sightings data/Site counts/",GridID, "_site_count_yearly.csv", sep = ""))
# Create a sightings/site dataframe
Grid_dataset_spp_sightings_by_site <- left_join(Grid_dataset_site_count_yearly, Grid_dataset_common_birds_yearly)
Grid_dataset_spp_sightings_by_site <- Grid_dataset_spp_sightings_by_site %>% 
  mutate(Yearly_Sightings_Over_Sites = Total_Yearly_Sightings/Total_Yearly_Sites)
write_csv(Grid_dataset_site_count_yearly, paste("Routine 6/Site and sightings data/Sightings over sites yearly/",GridID, "_sightings_over_sites_yearly.csv", sep = ""))
# Remove dataframes
rm(Grid_dataset, Grid_dataset_common_birds, Grid_dataset_common_birds_yearly, Grid_dataset_site_count_yearly, Grid_dataset_spp, Grid_dataset_spp_sightings_by_site)
}

