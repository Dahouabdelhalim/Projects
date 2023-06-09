# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 3
# ------------------------------------------------

# Load Libraries
library(dplyr)
library(stringr)
library(readr)

# Set Working Directory
setwd("~/Data repository")

# ------------------------------------------------
# Read in datasets
# ------------------------------------------------
# ------------------------------------------------
# Connecticut
# ------------------------------------------------
# Read in regional file
Connecticut_2 <- read_csv(file = "Routine 2/Connecticut_2.csv")
# Get unique grids in region
Connecticut_2_uGridID <- Connecticut_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Connecticut_2_uGridID_x <- Connecticut_2_uGridID
Connecticut_2_uGridID_x <- as.data.frame(str_replace_all(Connecticut_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Connecticut_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Connecticut_2_uGridID$GridID)
# Get region name
Region <- "Connecticut"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Connecticut_2, GridID == Connecticut_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Connecticut_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Connecticut_2 <- Connecticut_2_uGridID <- Connecticut_2_uGridID_x <- NULL

# ------------------------------------------------
# Delaware
# ------------------------------------------------
# Read in regional file
Delaware_2 <- read_csv(file = "Routine 2/Delaware_2.csv")
# Get unique grids in region
Delaware_2_uGridID <- Delaware_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Delaware_2_uGridID_x <- Delaware_2_uGridID
Delaware_2_uGridID_x <- as.data.frame(str_replace_all(Delaware_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Delaware_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Delaware_2_uGridID$GridID)
# Get region name
Region <- "Delaware"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Delaware_2, GridID == Delaware_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Delaware_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Delaware_2 <- Delaware_2_uGridID <- Delaware_2_uGridID_x <- NULL

# ------------------------------------------------
# Illinois
# ------------------------------------------------
# Read in regional file
Illinois_2 <- read_csv(file = "Routine 2/Illinois_2.csv")
# Get unique grids in region
Illinois_2_uGridID <- Illinois_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Illinois_2_uGridID_x <- Illinois_2_uGridID
Illinois_2_uGridID_x <- as.data.frame(str_replace_all(Illinois_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Illinois_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Illinois_2_uGridID$GridID)
# Get region name
Region <- "Illinois"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Illinois_2, GridID == Illinois_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Illinois_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Illinois_2 <- Illinois_2_uGridID <- Illinois_2_uGridID_x <- NULL

# ------------------------------------------------
# Indiana
# ------------------------------------------------
# Read in regional file
Indiana_2 <- read_csv(file = "Routine 2/Indiana_2.csv")
# Get unique grids in region
Indiana_2_uGridID <- Indiana_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Indiana_2_uGridID_x <- Indiana_2_uGridID
Indiana_2_uGridID_x <- as.data.frame(str_replace_all(Indiana_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Indiana_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Indiana_2_uGridID$GridID)
# Get region name
Region <- "Indiana"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Indiana_2, GridID == Indiana_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Indiana_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Indiana_2 <- Indiana_2_uGridID <- Indiana_2_uGridID_x <- NULL

# ------------------------------------------------
# Kentucky
# ------------------------------------------------
# Read in regional file
Kentucky_2 <- read_csv(file = "Routine 2/Kentucky_2.csv")
# Get unique grids in region
Kentucky_2_uGridID <- Kentucky_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Kentucky_2_uGridID_x <- Kentucky_2_uGridID
Kentucky_2_uGridID_x <- as.data.frame(str_replace_all(Kentucky_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Kentucky_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Kentucky_2_uGridID$GridID)
# Get region name
Region <- "Kentucky"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Kentucky_2, GridID == Kentucky_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Kentucky_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Kentucky_2 <- Kentucky_2_uGridID <- Kentucky_2_uGridID_x <- NULL

# ------------------------------------------------
# Maine
# ------------------------------------------------
# Read in regional file
Maine_2 <- read_csv(file = "Routine 2/Maine_2.csv")
# Get unique grids in region
Maine_2_uGridID <- Maine_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Maine_2_uGridID_x <- Maine_2_uGridID
Maine_2_uGridID_x <- as.data.frame(str_replace_all(Maine_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Maine_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Maine_2_uGridID$GridID)
# Get region name
Region <- "Maine"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Maine_2, GridID == Maine_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Maine_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Maine_2 <- Maine_2_uGridID <- Maine_2_uGridID_x <- NULL

# ------------------------------------------------
# Maryland
# ------------------------------------------------
# Read in regional file
Maryland_2 <- read_csv(file = "Routine 2/Maryland_2.csv")
# Get unique grids in region
Maryland_2_uGridID <- Maryland_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Maryland_2_uGridID_x <- Maryland_2_uGridID
Maryland_2_uGridID_x <- as.data.frame(str_replace_all(Maryland_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Maryland_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Maryland_2_uGridID$GridID)
# Get region name
Region <- "Maryland"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Maryland_2, GridID == Maryland_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Maryland_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Maryland_2 <- Maryland_2_uGridID <- Maryland_2_uGridID_x <- NULL

# ------------------------------------------------
# Massachusetts
# ------------------------------------------------
# Read in regional file
Massachusetts_2 <- read_csv(file = "Routine 2/Massachusetts_2.csv")
# Get unique grids in region
Massachusetts_2_uGridID <- Massachusetts_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Massachusetts_2_uGridID_x <- Massachusetts_2_uGridID
Massachusetts_2_uGridID_x <- as.data.frame(str_replace_all(Massachusetts_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Massachusetts_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Massachusetts_2_uGridID$GridID)
# Get region name
Region <- "Massachusetts"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Massachusetts_2, GridID == Massachusetts_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Massachusetts_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Massachusetts_2 <- Massachusetts_2_uGridID <- Massachusetts_2_uGridID_x <- NULL

# ------------------------------------------------
# Michigan
# ------------------------------------------------
# Read in regional file
Michigan_2 <- read_csv(file = "Routine 2/Michigan_2.csv")
# Get unique grids in region
Michigan_2_uGridID <- Michigan_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Michigan_2_uGridID_x <- Michigan_2_uGridID
Michigan_2_uGridID_x <- as.data.frame(str_replace_all(Michigan_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Michigan_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Michigan_2_uGridID$GridID)
# Get region name
Region <- "Michigan"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Michigan_2, GridID == Michigan_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Michigan_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Michigan_2 <- Michigan_2_uGridID <- Michigan_2_uGridID_x <- NULL

# ------------------------------------------------
# New_Brunswick
# ------------------------------------------------
# Read in regional file
New_Brunswick_2 <- read_csv(file = "Routine 2/New_Brunswick_2.csv")
# Get unique grids in region
New_Brunswick_2_uGridID <- New_Brunswick_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
New_Brunswick_2_uGridID_x <- New_Brunswick_2_uGridID
New_Brunswick_2_uGridID_x <- as.data.frame(str_replace_all(New_Brunswick_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(New_Brunswick_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(New_Brunswick_2_uGridID$GridID)
# Get region name
Region <- "New_Brunswick"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(New_Brunswick_2, GridID == New_Brunswick_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", New_Brunswick_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
New_Brunswick_2 <- New_Brunswick_2_uGridID <- New_Brunswick_2_uGridID_x <- NULL

# ------------------------------------------------
# New_Hampshire
# ------------------------------------------------
# Read in regional file
New_Hampshire_2 <- read_csv(file = "Routine 2/New_Hampshire_2.csv")
# Get unique grids in region
New_Hampshire_2_uGridID <- New_Hampshire_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
New_Hampshire_2_uGridID_x <- New_Hampshire_2_uGridID
New_Hampshire_2_uGridID_x <- as.data.frame(str_replace_all(New_Hampshire_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(New_Hampshire_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(New_Hampshire_2_uGridID$GridID)
# Get region name
Region <- "New_Hampshire"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(New_Hampshire_2, GridID == New_Hampshire_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", New_Hampshire_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
New_Hampshire_2 <- New_Hampshire_2_uGridID <- New_Hampshire_2_uGridID_x <- NULL

# ------------------------------------------------
# New_Jersey
# ------------------------------------------------
# Read in regional file
New_Jersey_2 <- read_csv(file = "Routine 2/New_Jersey_2.csv")
# Get unique grids in region
New_Jersey_2_uGridID <- New_Jersey_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
New_Jersey_2_uGridID_x <- New_Jersey_2_uGridID
New_Jersey_2_uGridID_x <- as.data.frame(str_replace_all(New_Jersey_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(New_Jersey_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(New_Jersey_2_uGridID$GridID)
# Get region name
Region <- "New_Jersey"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(New_Jersey_2, GridID == New_Jersey_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", New_Jersey_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
New_Jersey_2 <- New_Jersey_2_uGridID <- New_Jersey_2_uGridID_x <- NULL

# ------------------------------------------------
# New_York
# ------------------------------------------------
# Read in regional file
New_York_2 <- read_csv(file = "Routine 2/New_York_2.csv")
# Get unique grids in region
New_York_2_uGridID <- New_York_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
New_York_2_uGridID_x <- New_York_2_uGridID
New_York_2_uGridID_x <- as.data.frame(str_replace_all(New_York_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(New_York_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(New_York_2_uGridID$GridID)
# Get region name
Region <- "New_York"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(New_York_2, GridID == New_York_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", New_York_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
New_York_2 <- New_York_2_uGridID <- New_York_2_uGridID_x <- NULL

# ------------------------------------------------
# Nova_Scotia
# ------------------------------------------------
# Read in regional file
Nova_Scotia_2 <- read_csv(file = "Routine 2/Nova_Scotia_2.csv")
# Get unique grids in region
Nova_Scotia_2_uGridID <- Nova_Scotia_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Nova_Scotia_2_uGridID_x <- Nova_Scotia_2_uGridID
Nova_Scotia_2_uGridID_x <- as.data.frame(str_replace_all(Nova_Scotia_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Nova_Scotia_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Nova_Scotia_2_uGridID$GridID)
# Get region name
Region <- "Nova_Scotia"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Nova_Scotia_2, GridID == Nova_Scotia_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Nova_Scotia_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Nova_Scotia_2 <- Nova_Scotia_2_uGridID <- Nova_Scotia_2_uGridID_x <- NULL

# ------------------------------------------------
# Ohio
# ------------------------------------------------
# Read in regional file
Ohio_2 <- read_csv(file = "Routine 2/Ohio_2.csv")
# Get unique grids in region
Ohio_2_uGridID <- Ohio_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Ohio_2_uGridID_x <- Ohio_2_uGridID
Ohio_2_uGridID_x <- as.data.frame(str_replace_all(Ohio_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Ohio_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Ohio_2_uGridID$GridID)
# Get region name
Region <- "Ohio"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Ohio_2, GridID == Ohio_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Ohio_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Ohio_2 <- Ohio_2_uGridID <- Ohio_2_uGridID_x <- NULL

# ------------------------------------------------
# Ontario
# ------------------------------------------------
# Read in regional file
Ontario_2 <- read_csv(file = "Routine 2/Ontario_2.csv")
# Get unique grids in region
Ontario_2_uGridID <- Ontario_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Ontario_2_uGridID_x <- Ontario_2_uGridID
Ontario_2_uGridID_x <- as.data.frame(str_replace_all(Ontario_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Ontario_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Ontario_2_uGridID$GridID)
# Get region name
Region <- "Ontario"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Ontario_2, GridID == Ontario_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Ontario_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Ontario_2 <- Ontario_2_uGridID <- Ontario_2_uGridID_x <- NULL

# ------------------------------------------------
# Pennsylvania
# ------------------------------------------------
# Read in regional file
Pennsylvania_2 <- read_csv(file = "Routine 2/Pennsylvania_2.csv")
# Get unique grids in region
Pennsylvania_2_uGridID <- Pennsylvania_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Pennsylvania_2_uGridID_x <- Pennsylvania_2_uGridID
Pennsylvania_2_uGridID_x <- as.data.frame(str_replace_all(Pennsylvania_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Pennsylvania_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Pennsylvania_2_uGridID$GridID)
# Get region name
Region <- "Pennsylvania"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Pennsylvania_2, GridID == Pennsylvania_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Pennsylvania_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Pennsylvania_2 <- Pennsylvania_2_uGridID <- Pennsylvania_2_uGridID_x <- NULL

# ------------------------------------------------
# Quebec
# ------------------------------------------------
# Read in regional file
Quebec_2 <- read_csv(file = "Routine 2/Quebec_2.csv")
# Get unique grids in region
Quebec_2_uGridID <- Quebec_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Quebec_2_uGridID_x <- Quebec_2_uGridID
Quebec_2_uGridID_x <- as.data.frame(str_replace_all(Quebec_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Quebec_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Quebec_2_uGridID$GridID)
# Get region name
Region <- "Quebec"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Quebec_2, GridID == Quebec_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Quebec_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Quebec_2 <- Quebec_2_uGridID <- Quebec_2_uGridID_x <- NULL

# ------------------------------------------------
# Rhode_Island
# ------------------------------------------------
# Read in regional file
Rhode_Island_2 <- read_csv(file = "Routine 2/Rhode_Island_2.csv")
# Get unique grids in region
Rhode_Island_2_uGridID <- Rhode_Island_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Rhode_Island_2_uGridID_x <- Rhode_Island_2_uGridID
Rhode_Island_2_uGridID_x <- as.data.frame(str_replace_all(Rhode_Island_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Rhode_Island_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Rhode_Island_2_uGridID$GridID)
# Get region name
Region <- "Rhode_Island"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Rhode_Island_2, GridID == Rhode_Island_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Rhode_Island_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Rhode_Island_2 <- Rhode_Island_2_uGridID <- Rhode_Island_2_uGridID_x <- NULL

# ------------------------------------------------
# Vermont
# ------------------------------------------------
# Read in regional file
Vermont_2 <- read_csv(file = "Routine 2/Vermont_2.csv")
# Get unique grids in region
Vermont_2_uGridID <- Vermont_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Vermont_2_uGridID_x <- Vermont_2_uGridID
Vermont_2_uGridID_x <- as.data.frame(str_replace_all(Vermont_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Vermont_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Vermont_2_uGridID$GridID)
# Get region name
Region <- "Vermont"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Vermont_2, GridID == Vermont_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Vermont_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Vermont_2 <- Vermont_2_uGridID <- Vermont_2_uGridID_x <- NULL

# ------------------------------------------------
# Virginia
# ------------------------------------------------
# Read in regional file
Virginia_2 <- read_csv(file = "Routine 2/Virginia_2.csv")
# Get unique grids in region
Virginia_2_uGridID <- Virginia_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Virginia_2_uGridID_x <- Virginia_2_uGridID
Virginia_2_uGridID_x <- as.data.frame(str_replace_all(Virginia_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Virginia_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Virginia_2_uGridID$GridID)
# Get region name
Region <- "Virginia"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Virginia_2, GridID == Virginia_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Virginia_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Virginia_2 <- Virginia_2_uGridID <- Virginia_2_uGridID_x <- NULL

# ------------------------------------------------
# West_Virginia
# ------------------------------------------------
# Read in regional file
West_Virginia_2 <- read_csv(file = "Routine 2/West_Virginia_2.csv")
# Get unique grids in region
West_Virginia_2_uGridID <- West_Virginia_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
West_Virginia_2_uGridID_x <- West_Virginia_2_uGridID
West_Virginia_2_uGridID_x <- as.data.frame(str_replace_all(West_Virginia_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(West_Virginia_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(West_Virginia_2_uGridID$GridID)
# Get region name
Region <- "West_Virginia"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(West_Virginia_2, GridID == West_Virginia_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", West_Virginia_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
West_Virginia_2 <- West_Virginia_2_uGridID <- West_Virginia_2_uGridID_x <- NULL

# ------------------------------------------------
# Wisconsin
# ------------------------------------------------
# Read in regional file
Wisconsin_2 <- read_csv(file = "Routine 2/Wisconsin_2.csv")
# Get unique grids in region
Wisconsin_2_uGridID <- Wisconsin_2 %>% 
  dplyr::select(GridID) %>% 
  distinct()
# Replace the * with x
Wisconsin_2_uGridID_x <- Wisconsin_2_uGridID
Wisconsin_2_uGridID_x <- as.data.frame(str_replace_all(Wisconsin_2_uGridID_x$GridID, "\\\\*", "x"))
colnames(Wisconsin_2_uGridID_x) <- "GridID"
# Determine number of unique grids in region
ngrids <- length(Wisconsin_2_uGridID$GridID)
# Get region name
Region <- "Wisconsin"
# Filter and then write dataframes for each grid within the region
for (i in 1:ngrids) {
  gridsub <- filter(Wisconsin_2, GridID == Wisconsin_2_uGridID$GridID[i])
  write_csv(gridsub, paste("Routine 3/", Wisconsin_2_uGridID_x$GridID[i], "_", Region, ".csv", sep = ""))
  gridsub <- NULL
}
# Clean up
Wisconsin_2 <- Wisconsin_2_uGridID <- Wisconsin_2_uGridID_x <- NULL
