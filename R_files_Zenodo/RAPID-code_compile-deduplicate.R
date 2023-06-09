# RAPID code for compiling raw data from multiple sources and deduplicating records
# This code accompanies document 'RAPID-protocol_compile-deduplicate.pdf'
# Code written by Erica Krimmel, 2020-07-31. Updated 2020-09-23.

# Load core libraries; install this package if you have not already
library(tidyverse)

# Load data GBIF R package for data quality control further down; install
# this package if you have not already
library(rgbif)

# PREPARE DATA FROM IDIGBIO

# Read into R the raw occurrence data from iDigBio, which should be whatever was
# published by the data provider (e.g. the collection)
idb_raw <- read_csv("9d4b9069-48c4-4212-90d8-4dd6f4b7f2a5/occurrence_raw.csv", 
                    na = c("", "NA"),
                    col_types = cols(.default = col_character()))

# Rename columns to reflect provenance and remove colon characters
idb_raw <- idb_raw %>% 
  rename_all(function(x){paste0(x, "_idbR")}) %>% 
  rename_all(funs(str_replace_all(., "dwc:", ""))) %>% 
  rename_all(funs(str_replace_all(., ":", "_")))

# Read into R the version of occurrence data processed by iDigBio
idb_processed <- read_csv("9d4b9069-48c4-4212-90d8-4dd6f4b7f2a5/occurrence.csv", 
                          na = c("", "NA"),
                          col_types = cols(.default = col_character()))

# Rename columns to reflect provenance and remove illegal characters
idb_processed <- idb_processed %>% 
  rename_all(function(x){paste0(x, "_idbP")}) %>% 
  rename_all(funs(str_replace_all(., "dwc:", ""))) %>% 
  rename_all(funs(str_replace_all(., ":", "_")))

# Join raw and processed iDigBio data together
idb_joined <- idb_raw %>% 
  left_join(idb_processed, by = c("coreid_idbR" = "coreid_idbP")) %>% 
  # Subset by families of interest
  mutate(family_idbR = tolower(family_idbR)) %>% 
  filter(family_idbR %in% c("rhinolophidae", 
                            "hipposideridae", 
                            "rhinonycteridae") |
           family_idbP %in% c("rhinolophidae", 
                              "hipposideridae", 
                              "rhinonycteridae")) %>% 
  # Expand the idbP field `geopoint` into two fields to reflect the same
  # structure as the other datasets here
  separate(idigbio_geoPoint_idbP, 
           c("decimalLatitude_idbP", "decimalLongitude_idbP"), 
           sep = ",") %>% 
  mutate(decimalLatitude_idbP = parse_number(decimalLatitude_idbP)) %>% 
  mutate(decimalLongitude_idbP = parse_number(decimalLongitude_idbP)) %>% 
  # Create a new column to identify duplicate records between iDigBio and GBIF
  # in future steps
  unite(matchDuplicates, 
        c("institutionCode_idbR", 
          "collectionCode_idbR", 
          "catalogNumber_idbR",
          "occurrenceID_idbR"), 
        sep = " ", 
        remove = FALSE) %>% 
  # Normalize text differences by lowercasing
  mutate(matchDuplicates = str_squish(tolower(matchDuplicates)))

# PREPARE DATA FROM GBIF

# Read into R the raw occurrence data from GBIF, which should be whatever was
# published by the data provider (e.g. the collection); this should ostensibly
# be the same as what is in `idb_raw`
gbif_subset_1R <- read_tsv("0067804-200613084148143/verbatim.txt", 
                           na = c("", "NA"),
                           col_types = cols(.default = col_character()))

# Rename columns to reflect provenance
gbif_subset_1R <- gbif_subset_1R %>% 
  rename_all(function(x){paste0(x, "_gbifR")})

# Read into R the version of occurrence data processed by GBIF
gbif_subset_1P <- read_tsv("0067804-200613084148143/occurrence.txt", 
                           na = c("", "NA"),
                           col_types = cols(acceptedTaxonKey = col_integer(),
                                            classKey = col_integer(),
                                            familyKey = col_integer(),
                                            genusKey = col_integer(),
                                            kingdomKey = col_integer(),
                                            orderKey = col_integer(),
                                            phylumKey = col_integer(),
                                            speciesKey = col_integer(),
                                            .default = col_character()))

# Rename columns to reflect provenance
gbif_subset_1P <- gbif_subset_1P %>% 
  rename_all(function(x){paste0(x, "_gbifP")})

# Join raw and processed GBIF subset data together
gbif_subset_1 <- gbif_subset_1R %>% 
  left_join(gbif_subset_1P, by = c("gbifID_gbifR" = "gbifID_gbifP")) %>% 
  # Get rid of these columns because I can't figure out how to coerce them into
  # the same data class to later join GBIF subsets
  select(-organismQuantity_gbifR, -organismQuantity_gbifP)

# Read into R the raw occurrence data from GBIF, which should be whatever was
# published by the data provider (e.g. the collection); this should ostensibly
# be the same as what is in `idb_raw`
gbif_subset_2R <- read_tsv("0067806-200613084148143/verbatim.txt", 
                           #  na = c("", "NA"),
                           col_types = cols(.default = col_character()))

# Rename columns to reflect provenance
gbif_subset_2R <- gbif_subset_2R %>% 
  rename_all(function(x){paste0(x, "_gbifR")})

# Read into R the version of occurrence data processed by GBIF
gbif_subset_2P <- read_tsv("0067806-200613084148143/occurrence.txt", 
                           # na = c("", "NA"),
                           col_types = cols(acceptedTaxonKey = col_integer(),
                                            classKey = col_integer(),
                                            familyKey = col_integer(),
                                            genusKey = col_integer(),
                                            kingdomKey = col_integer(),
                                            orderKey = col_integer(),
                                            phylumKey = col_integer(),
                                            speciesKey = col_integer(),
                                            .default = col_character()))

# Rename columns to reflect provenance
gbif_subset_2P <- gbif_subset_2P %>% 
  rename_all(function(x){paste0(x, "_gbifP")})

# Join raw and processed GBIF subset data together
gbif_subset_2 <- gbif_subset_2R %>% 
  left_join(gbif_subset_2P, by = c("gbifID_gbifR" = "gbifID_gbifP")) %>% 
  # Get rid of these columns because I can't figure out how to coerce them into
  # the same data class to later join GBIF subsets
  select(-organismQuantity_gbifR, -organismQuantity_gbifP)

# Bind `gbif_subset_1` and `gbif_subset_2` together
gbif_joined <- bind_rows(list(gbif_subset_1, gbif_subset_2), .id = "id") %>%
  # Create a new column to identify duplicate records between iDigBio and GBIF
  # in future steps
  unite(matchDuplicates, 
        c("institutionCode_gbifR", 
          "collectionCode_gbifR", 
          "catalogNumber_gbifR",
          "occurrenceID_gbifR"), 
        sep = " ", 
        remove = FALSE) %>%
  # Normalize text differences by lowercasing
  mutate(matchDuplicates = str_squish(tolower(matchDuplicates)))

# Define function that will fetch publisher information for a GBIF dataset
getPublisher_gbif <- function(datasetKey) {
  publisher <- rgbif::datasets(uuid = datasetKey)
  publisher <- paste(publisher$data$publishingOrganizationKey,
                     publisher$data$title,
                     publisher$data$type,
                     sep = "xxxxx")
  publisher
}

# Compile list of all datasets and publishers contributing GBIF records
gbif_datasets <- gbif_joined %>% 
  select(datasetKey_gbifP) %>% 
  distinct() %>% 
  mutate(publisher = map_chr(datasetKey_gbifP, getPublisher_gbif)) %>% 
  separate(publisher, 
           into = c("publisherKey_gbifP", 
                    "publisherTitle_gbifP",
                    "publisherType_gbifP"),
           sep = "xxxxx")

# Add publisher information to `records` data so that we can exclude checklist
# datasets, which are are linked directly to physical specimens
gbif_joined <- gbif_joined %>% 
  left_join(gbif_datasets, by = "datasetKey_gbifP") %>% 
  filter(publisherType_gbifP != "CHECKLIST")

# JOIN IDIGBIO AND GBIF DATA TOGETHER

# Join `idb_joined` and `gbif_joined` based on unique combinations of
# institution code, collection code, catalog number, and occurrence ID
records <- gbif_joined %>% 
  full_join(idb_joined, by = "matchDuplicates") %>%
  # Record which records were present in iDigBio data
  mutate(idigbio = case_when(!is.na(idigbio_uuid_idbP) ~ 1,
                             is.na(idigbio_uuid_idbP) ~ 0)) %>% 
  # Record which records were present in GBIF data
  mutate(gbif = case_when(!is.na(gbifID_gbifR) ~ 1,
                          is.na(gbifID_gbifR) ~ 0)) %>% 
  # Rearrange columns by provenance, unique IDs, and then alphabetically
  select(gbif, idigbio, gbifID_gbifR, idigbio_uuid_idbP, 
         sort(tidyselect::peek_vars()),
         # Remove column used to match duplicate records between iDigBio and GBIF
         -starts_with("matchDuplicates"))

# Save `records` as CSV file
write_csv(records, 
          paste("rapid-joined-records_", Sys.Date(), ".csv", sep = ""),
          na = "")