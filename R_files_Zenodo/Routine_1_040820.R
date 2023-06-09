# ------------------------------------------------
# eBird Project Data Preparation Code: Routine 1
# ------------------------------------------------

# Load Libraries
library(dplyr)
library(readr)
library(lubridate)
library(sp)

# Set Working Directory
setwd("~/Data Repository")

# ------------------------------------------------
# Start function subroutines
# ------------------------------------------------
CRS_swap = function(data,
                    src.proj = CRS("+init=epsg:4326"),
                    dst.proj = CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ")) {
  require(sp)
  as.data.frame(
    spTransform(
      SpatialPointsDataFrame(
        coords = data.frame(Xnad = toconv$LONGITUDE,
                            Ynad = toconv$LATITUDE),
        data = data.frame(toconv$`COMMON NAME`,
                          toconv$`LOCALITY ID`,
                          toconv$Year,
                          toconv$Jd,
                          toconv$Region,
                          Xlon = toconv$LONGITUDE,
                          Ylat = toconv$LATITUDE),
        proj4string = src.proj), dst.proj))
  
}

# ------------------------------------------------
# Finish function subroutines
# ------------------------------------------------

# ------------------------------------------------
# Read in datasets
# ------------------------------------------------
# Read in the eBird Data for each location
(eBird_Data_dirs <- list.dirs(path = "eBird Data", full.names = T))
eBird_Data_dirs <- eBird_Data_dirs[2:24]
eBird_Data_filesnames_list <- list()
for (i in 1:length(eBird_Data_dirs)) {
  eBird_Data_filesnames_list[i] <- list.files(path = eBird_Data_dirs[i], pattern = "relFeb-2020.txt$")
}

# ------------------------------------------------
# Connecticut
# ------------------------------------------------
Connecticut <- data.table::fread(paste(eBird_Data_dirs[1], "/", eBird_Data_filesnames_list[1], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Connecticut$`OBSERVATION DATE` <- as.Date(Connecticut$`OBSERVATION DATE`)
Connecticut$Year <- year(Connecticut$`OBSERVATION DATE`)
Connecticut$Jd <- yday(Connecticut$`OBSERVATION DATE`)
head(Connecticut)
Connecticut$`OBSERVATION DATE` <- NULL
# Add region
Connecticut$Region <- "Connecticut"
# Remove NAs
Connecticut <- filter(Connecticut, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Connecticut
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Connecticut <- final
final <- toconv <- NULL

# Write csv
write_csv(Connecticut, "Routine 1/Connecticut_1.csv")
Connecticut <- NULL



# ------------------------------------------------
# Delaware
# ------------------------------------------------
Delaware <- data.table::fread(paste(eBird_Data_dirs[2], "/", eBird_Data_filesnames_list[2], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Clean and tidy dataframe
# Select relevant columns
Delaware_Orig <- Delaware
Delaware <- Delaware %>% 
  dplyr::select(`COMMON NAME`,
                `LOCALITY ID`,
                LATITUDE,
                LONGITUDE,
                `OBSERVATION DATE`
  )
head(Delaware)
Delaware_Orig <- NULL

# Convert observation date
Delaware$`OBSERVATION DATE` <- as.Date(Delaware$`OBSERVATION DATE`)
Delaware$Year <- year(Delaware$`OBSERVATION DATE`)
Delaware$Jd <- yday(Delaware$`OBSERVATION DATE`)
head(Delaware)
Delaware$`OBSERVATION DATE` <- NULL
# Add region
Delaware$Region <- "Delaware"
# Remove NAs
Delaware <- filter(Delaware, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Delaware
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Delaware <- final
final <- toconv <- NULL

# Write csv
write_csv(Delaware, "Routine 1/Delaware_1.csv")

Delaware <- NULL



# ------------------------------------------------
# Illinois
# ------------------------------------------------
Illinois <- data.table::fread(paste(eBird_Data_dirs[3], "/", eBird_Data_filesnames_list[3], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Illinois$`OBSERVATION DATE` <- as.Date(Illinois$`OBSERVATION DATE`)
Illinois$Year <- year(Illinois$`OBSERVATION DATE`)
Illinois$Jd <- yday(Illinois$`OBSERVATION DATE`)
head(Illinois)
Illinois$`OBSERVATION DATE` <- NULL
# Add region
Illinois$Region <- "Illinois"
# Remove NAs
Illinois <- filter(Illinois, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Illinois

final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Illinois <- final
final <- toconv <- NULL

# Write csv
write_csv(Illinois, "Routine 1/Illinois_1.csv")
Illinois <- NULL



# ------------------------------------------------
# Indiana
# ------------------------------------------------
Indiana <- data.table::fread(paste(eBird_Data_dirs[4], "/", eBird_Data_filesnames_list[4], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Indiana$`OBSERVATION DATE` <- as.Date(Indiana$`OBSERVATION DATE`)
Indiana$Year <- year(Indiana$`OBSERVATION DATE`)
Indiana$Jd <- yday(Indiana$`OBSERVATION DATE`)
head(Indiana)
Indiana$`OBSERVATION DATE` <- NULL
# Add region
Indiana$Region <- "Indiana"
# Remove NAs
Indiana <- filter(Indiana, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Indiana
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Indiana <- final
final <- toconv <- NULL

# Write csv
write_csv(Indiana, "Routine 1/Indiana_1.csv")
Indiana <- NULL



# ------------------------------------------------
# Kentucky
# ------------------------------------------------
Kentucky <- data.table::fread(paste(eBird_Data_dirs[5], "/", eBird_Data_filesnames_list[5], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Kentucky$`OBSERVATION DATE` <- as.Date(Kentucky$`OBSERVATION DATE`)
Kentucky$Year <- year(Kentucky$`OBSERVATION DATE`)
Kentucky$Jd <- yday(Kentucky$`OBSERVATION DATE`)
head(Kentucky)
Kentucky$`OBSERVATION DATE` <- NULL
# Add region
Kentucky$Region <- "Kentucky"
# Remove NAs
Kentucky <- filter(Kentucky, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Kentucky
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Kentucky <- final
final <- toconv <- NULL

# Write csv
write_csv(Kentucky, "Routine 1/Kentucky_1.csv")
Kentucky <- NULL



# ------------------------------------------------
# Maine
# ------------------------------------------------
Maine <- data.table::fread(paste(eBird_Data_dirs[6], "/", eBird_Data_filesnames_list[6], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Maine$`OBSERVATION DATE` <- as.Date(Maine$`OBSERVATION DATE`)
Maine$Year <- year(Maine$`OBSERVATION DATE`)
Maine$Jd <- yday(Maine$`OBSERVATION DATE`)
head(Maine)
Maine$`OBSERVATION DATE` <- NULL
# Add region
Maine$Region <- "Maine"
# Remove NAs
Maine <- filter(Maine, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Maine
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Maine <- final
final <- toconv <- NULL

# Write csv
write_csv(Maine, "Routine 1/Maine_1.csv")
Maine <- NULL



# ------------------------------------------------
# Maryland
# ------------------------------------------------
Maryland <- data.table::fread(paste(eBird_Data_dirs[7], "/", eBird_Data_filesnames_list[7], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Maryland$`OBSERVATION DATE` <- as.Date(Maryland$`OBSERVATION DATE`)
Maryland$Year <- year(Maryland$`OBSERVATION DATE`)
Maryland$Jd <- yday(Maryland$`OBSERVATION DATE`)
head(Maryland)
Maryland$`OBSERVATION DATE` <- NULL
# Add region
Maryland$Region <- "Maryland"
# Remove NAs
Maryland <- filter(Maryland, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Maryland
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Maryland <- final
final <- toconv <- NULL

# Write csv
write_csv(Maryland, "Routine 1/Maryland_1.csv")
Maryland <- NULL



# ------------------------------------------------
# Massachusetts
# ------------------------------------------------
Massachusetts <- data.table::fread(paste(eBird_Data_dirs[8], "/", eBird_Data_filesnames_list[8], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Massachusetts$`OBSERVATION DATE` <- as.Date(Massachusetts$`OBSERVATION DATE`)
Massachusetts$Year <- year(Massachusetts$`OBSERVATION DATE`)
Massachusetts$Jd <- yday(Massachusetts$`OBSERVATION DATE`)
head(Massachusetts)
Massachusetts$`OBSERVATION DATE` <- NULL
# Add region
Massachusetts$Region <- "Massachusetts"
# Remove NAs
Massachusetts <- filter(Massachusetts, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Massachusetts
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Massachusetts <- final
final <- toconv <- NULL

# Write csv
write_csv(Massachusetts, "Routine 1/Massachusetts_1.csv")
Massachusetts <- NULL



# ------------------------------------------------
# Michigan
# ------------------------------------------------
Michigan <- data.table::fread(paste(eBird_Data_dirs[9], "/", eBird_Data_filesnames_list[9], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Michigan$`OBSERVATION DATE` <- as.Date(Michigan$`OBSERVATION DATE`)
Michigan$Year <- year(Michigan$`OBSERVATION DATE`)
Michigan$Jd <- yday(Michigan$`OBSERVATION DATE`)
head(Michigan)
Michigan$`OBSERVATION DATE` <- NULL
# Add region
Michigan$Region <- "Michigan"
# Remove NAs
Michigan <- filter(Michigan, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Michigan
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Michigan <- final
final <- toconv <- NULL

# Write csv
write_csv(Michigan, "Routine 1/Michigan_1.csv")
Michigan <- NULL



# ------------------------------------------------
# New_Brunswick
# ------------------------------------------------
New_Brunswick <- data.table::fread(paste(eBird_Data_dirs[10], "/", eBird_Data_filesnames_list[10], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
New_Brunswick$`OBSERVATION DATE` <- as.Date(New_Brunswick$`OBSERVATION DATE`)
New_Brunswick$Year <- year(New_Brunswick$`OBSERVATION DATE`)
New_Brunswick$Jd <- yday(New_Brunswick$`OBSERVATION DATE`)
head(New_Brunswick)
New_Brunswick$`OBSERVATION DATE` <- NULL
# Add region
New_Brunswick$Region <- "New_Brunswick"
# Remove NAs
New_Brunswick <- filter(New_Brunswick, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- New_Brunswick
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
New_Brunswick <- final
final <- toconv <- NULL

# Write csv
write_csv(New_Brunswick, "Routine 1/New_Brunswick_1.csv")
New_Brunswick <- NULL



# ------------------------------------------------
# New_Hampshire
# ------------------------------------------------
New_Hampshire <- data.table::fread(paste(eBird_Data_dirs[11], "/", eBird_Data_filesnames_list[11], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
New_Hampshire$`OBSERVATION DATE` <- as.Date(New_Hampshire$`OBSERVATION DATE`)
New_Hampshire$Year <- year(New_Hampshire$`OBSERVATION DATE`)
New_Hampshire$Jd <- yday(New_Hampshire$`OBSERVATION DATE`)
head(New_Hampshire)
New_Hampshire$`OBSERVATION DATE` <- NULL
# Add region
New_Hampshire$Region <- "New_Hampshire"
# Remove NAs
New_Hampshire <- filter(New_Hampshire, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- New_Hampshire
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
New_Hampshire <- final
final <- toconv <- NULL

# Write csv
write_csv(New_Hampshire, "Routine 1/New_Hampshire_1.csv")
New_Hampshire <- NULL



# ------------------------------------------------
# New_Jersey
# ------------------------------------------------
New_Jersey <- data.table::fread(paste(eBird_Data_dirs[12], "/", eBird_Data_filesnames_list[12], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
New_Jersey$`OBSERVATION DATE` <- as.Date(New_Jersey$`OBSERVATION DATE`)
New_Jersey$Year <- year(New_Jersey$`OBSERVATION DATE`)
New_Jersey$Jd <- yday(New_Jersey$`OBSERVATION DATE`)
head(New_Jersey)
New_Jersey$`OBSERVATION DATE` <- NULL
# Add region
New_Jersey$Region <- "New_Jersey"
# Remove NAs
New_Jersey <- filter(New_Jersey, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- New_Jersey
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
New_Jersey <- final
final <- toconv <- NULL

# Write csv
write_csv(New_Jersey, "Routine 1/New_Jersey_1.csv")
New_Jersey <- NULL



# ------------------------------------------------
# New_York
# ------------------------------------------------
New_York <- data.table::fread(paste(eBird_Data_dirs[13], "/", eBird_Data_filesnames_list[13], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
New_York$`OBSERVATION DATE` <- as.Date(New_York$`OBSERVATION DATE`)
New_York$Year <- year(New_York$`OBSERVATION DATE`)
New_York$Jd <- yday(New_York$`OBSERVATION DATE`)
head(New_York)
New_York$`OBSERVATION DATE` <- NULL
# Add region
New_York$Region <- "New_York"
# Remove NAs
New_York <- filter(New_York, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- New_York
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
New_York <- final
final <- toconv <- NULL

# Write csv
write_csv(New_York, "Routine 1/New_York_1.csv")
New_York <- NULL



# ------------------------------------------------
# Nova_Scotia
# ------------------------------------------------
Nova_Scotia <- data.table::fread(paste(eBird_Data_dirs[14], "/", eBird_Data_filesnames_list[14], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Nova_Scotia$`OBSERVATION DATE` <- as.Date(Nova_Scotia$`OBSERVATION DATE`)
Nova_Scotia$Year <- year(Nova_Scotia$`OBSERVATION DATE`)
Nova_Scotia$Jd <- yday(Nova_Scotia$`OBSERVATION DATE`)
head(Nova_Scotia)
Nova_Scotia$`OBSERVATION DATE` <- NULL
# Add region
Nova_Scotia$Region <- "Nova_Scotia"
# Remove NAs
Nova_Scotia <- filter(Nova_Scotia, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Nova_Scotia
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Nova_Scotia <- final
final <- toconv <- NULL

# Write csv
write_csv(Nova_Scotia, "Routine 1/Nova_Scotia_1.csv")
Nova_Scotia <- NULL



# ------------------------------------------------
# Ohio
# ------------------------------------------------
Ohio <- data.table::fread(paste(eBird_Data_dirs[15], "/", eBird_Data_filesnames_list[15], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Ohio$`OBSERVATION DATE` <- as.Date(Ohio$`OBSERVATION DATE`)
Ohio$Year <- year(Ohio$`OBSERVATION DATE`)
Ohio$Jd <- yday(Ohio$`OBSERVATION DATE`)
head(Ohio)
Ohio$`OBSERVATION DATE` <- NULL
# Add region
Ohio$Region <- "Ohio"
# Remove NAs
Ohio <- filter(Ohio, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Ohio
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Ohio <- final
final <- toconv <- NULL

# Write csv
write_csv(Ohio, "Routine 1/Ohio_1.csv")
Ohio <- NULL



# ------------------------------------------------
# Ontario
# ------------------------------------------------
Ontario <- data.table::fread(paste(eBird_Data_dirs[16], "/", eBird_Data_filesnames_list[16], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Ontario$`OBSERVATION DATE` <- as.Date(Ontario$`OBSERVATION DATE`)
Ontario$Year <- year(Ontario$`OBSERVATION DATE`)
Ontario$Jd <- yday(Ontario$`OBSERVATION DATE`)
head(Ontario)
Ontario$`OBSERVATION DATE` <- NULL
# Add region
Ontario$Region <- "Ontario"
# Remove NAs
Ontario <- filter(Ontario, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Ontario
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Ontario <- final
final <- toconv <- NULL

# Write csv
write_csv(Ontario, "Routine 1/Ontario_1.csv")
Ontario <- NULL



# ------------------------------------------------
# Pennsylvania
# ------------------------------------------------
Pennsylvania <- data.table::fread(paste(eBird_Data_dirs[17], "/", eBird_Data_filesnames_list[17], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Pennsylvania$`OBSERVATION DATE` <- as.Date(Pennsylvania$`OBSERVATION DATE`)
Pennsylvania$Year <- year(Pennsylvania$`OBSERVATION DATE`)
Pennsylvania$Jd <- yday(Pennsylvania$`OBSERVATION DATE`)
head(Pennsylvania)
Pennsylvania$`OBSERVATION DATE` <- NULL
# Add region
Pennsylvania$Region <- "Pennsylvania"
# Remove NAs
Pennsylvania <- filter(Pennsylvania, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Pennsylvania
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Pennsylvania <- final
final <- toconv <- NULL

# Write csv
write_csv(Pennsylvania, "Routine 1/Pennsylvania_1.csv")
Pennsylvania <- NULL



# ------------------------------------------------
# Quebec
# ------------------------------------------------
Quebec <- data.table::fread(paste(eBird_Data_dirs[18], "/", eBird_Data_filesnames_list[18], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Quebec$`OBSERVATION DATE` <- as.Date(Quebec$`OBSERVATION DATE`)
Quebec$Year <- year(Quebec$`OBSERVATION DATE`)
Quebec$Jd <- yday(Quebec$`OBSERVATION DATE`)
head(Quebec)
Quebec$`OBSERVATION DATE` <- NULL
# Add region
Quebec$Region <- "Quebec"
# Remove NAs
Quebec <- filter(Quebec, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Quebec
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Quebec <- final
final <- toconv <- NULL

# Write csv
write_csv(Quebec, "Routine 1/Quebec_1.csv")
Quebec <- NULL



# ------------------------------------------------
# Rhode_Island
# ------------------------------------------------
Rhode_Island <- data.table::fread(paste(eBird_Data_dirs[19], "/", eBird_Data_filesnames_list[19], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Rhode_Island$`OBSERVATION DATE` <- as.Date(Rhode_Island$`OBSERVATION DATE`)
Rhode_Island$Year <- year(Rhode_Island$`OBSERVATION DATE`)
Rhode_Island$Jd <- yday(Rhode_Island$`OBSERVATION DATE`)
head(Rhode_Island)
Rhode_Island$`OBSERVATION DATE` <- NULL
# Add region
Rhode_Island$Region <- "Rhode_Island"
# Remove NAs
Rhode_Island <- filter(Rhode_Island, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Rhode_Island
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Rhode_Island <- final
final <- toconv <- NULL

# Write csv
write_csv(Rhode_Island, "Routine 1/Rhode_Island_1.csv")
Rhode_Island <- NULL



# ------------------------------------------------
# Vermont
# ------------------------------------------------
Vermont <- data.table::fread(paste(eBird_Data_dirs[20], "/", eBird_Data_filesnames_list[20], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Vermont$`OBSERVATION DATE` <- as.Date(Vermont$`OBSERVATION DATE`)
Vermont$Year <- year(Vermont$`OBSERVATION DATE`)
Vermont$Jd <- yday(Vermont$`OBSERVATION DATE`)
head(Vermont)
Vermont$`OBSERVATION DATE` <- NULL
# Add region
Vermont$Region <- "Vermont"
# Remove NAs
Vermont <- filter(Vermont, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Vermont
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Vermont <- final
final <- toconv <- NULL

# Write csv
write_csv(Vermont, "Routine 1/Vermont_1.csv")
Vermont <- NULL



# ------------------------------------------------
# Virginia
# ------------------------------------------------
Virginia <- data.table::fread(paste(eBird_Data_dirs[21], "/", eBird_Data_filesnames_list[21], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Virginia$`OBSERVATION DATE` <- as.Date(Virginia$`OBSERVATION DATE`)
Virginia$Year <- year(Virginia$`OBSERVATION DATE`)
Virginia$Jd <- yday(Virginia$`OBSERVATION DATE`)
head(Virginia)
Virginia$`OBSERVATION DATE` <- NULL
# Add region
Virginia$Region <- "Virginia"
# Remove NAs
Virginia <- filter(Virginia, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Virginia
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Virginia <- final
final <- toconv <- NULL

# Write csv
write_csv(Virginia, "Routine 1/Virginia_1.csv")
Virginia <- NULL



# ------------------------------------------------
# West_Virginia
# ------------------------------------------------
West_Virginia <- data.table::fread(paste(eBird_Data_dirs[22], "/", eBird_Data_filesnames_list[22], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
West_Virginia$`OBSERVATION DATE` <- as.Date(West_Virginia$`OBSERVATION DATE`)
West_Virginia$Year <- year(West_Virginia$`OBSERVATION DATE`)
West_Virginia$Jd <- yday(West_Virginia$`OBSERVATION DATE`)
head(West_Virginia)
West_Virginia$`OBSERVATION DATE` <- NULL
# Add region
West_Virginia$Region <- "West_Virginia"
# Remove NAs
West_Virginia <- filter(West_Virginia, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- West_Virginia
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
West_Virginia <- final
final <- toconv <- NULL

# Write csv
write_csv(West_Virginia, "Routine 1/West_Virginia_1.csv")
West_Virginia <- NULL



# ------------------------------------------------
# Wisconsin
# ------------------------------------------------
Wisconsin <- data.table::fread(paste(eBird_Data_dirs[23], "/", eBird_Data_filesnames_list[23], sep = ""), select = c("COMMON NAME", "LOCALITY ID", "LATITUDE", "LONGITUDE", "OBSERVATION DATE"))

# Convert observation date
Wisconsin$`OBSERVATION DATE` <- as.Date(Wisconsin$`OBSERVATION DATE`)
Wisconsin$Year <- year(Wisconsin$`OBSERVATION DATE`)
Wisconsin$Jd <- yday(Wisconsin$`OBSERVATION DATE`)
head(Wisconsin)
Wisconsin$`OBSERVATION DATE` <- NULL
# Add region
Wisconsin$Region <- "Wisconsin"
# Remove NAs
Wisconsin <- filter(Wisconsin, !is.na(LATITUDE))

# Convert Lat/Long to 
NAD <- CRS("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
LL <- CRS("+init=epsg:4326")

toconv <- Wisconsin
final <- CRS_swap(data = toconv)
print(final)
# Correct column names
final <- final %>% 
  dplyr::select("COMMON NAME" = "toconv..COMMON.NAME.",
                "LOCALITY ID" = "toconv..LOCALITY.ID.",
                "Year" = "toconv.Year",
                "Jd" = "toconv.Jd",
                "Region" = "toconv.Region",
                "X" = "Xnad",
                "Y" = "Ynad")
Wisconsin <- final
final <- toconv <- NULL

# Write csv
write_csv(Wisconsin, "Routine 1/Wisconsin_1.csv")
Wisconsin <- NULL
