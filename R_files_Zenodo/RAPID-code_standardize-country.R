# RAPID code for standardizing country data
# This code accompanies document 'RAPID-protocol_standardize-country.pdf'
# Code written by Erica Krimmel, 2020-07-31. Updated 2020-09-23.

# Load core libraries; install this package if you have not already
library(tidyverse)

# Read into R the primary records dataset
records <- read_csv("rapid-joined-records_2020-09-23.csv", 
                    na = c("", "NA"),
                    col_types = cols(.default = col_character()))

# Extract distinct values for country information
countryCleanup <- records %>% 
  select(starts_with("country"), idigbio_isoCountryCode_idbP) %>% 
  distinct() %>% 
  unite(countryKey, everything(), remove = FALSE) %>% 
  mutate_at(vars(starts_with("country_")), tolower) %>% 
  mutate_at(vars(contains("countryCode")), toupper) %>% 
  mutate(country_rapid = coalesce(country_idbP,
                                  country_gbifR,
                                  country_idbR)) %>% 
  mutate(countryCode_rapid = coalesce(countryCode_gbifP,
                                      idigbio_isoCountryCode_idbP,
                                      countryCode_gbifR,
                                      countryCode_idbR))

# Save `countryCleanup` as file to do cleanup work in OpenRefine
write_csv(countryCleanup,
          paste("rapid-country-cleanup_", Sys.Date(), ".csv", sep = ""),
          na = "")

# Read in cleaned up country data file saved from OpenRefine
countryCleanup_done <- read_csv("rapid-country-cleanup_2020-09-23.csv",
                                na = c("", "NA"),
                                col_types = cols(.default = col_character()))

# Join cleaned up country data back into `records`
records_countryCleanup <- records %>% 
  unite(countryKey,
        country_gbifR,
        country_idbP,
        country_idbR,
        countryCode_gbifP,
        countryCode_gbifR,
        countryCode_idbR,
        idigbio_isoCountryCode_idbP,
        remove = FALSE) %>% 
  left_join(select(countryCleanup_done, 
                   countryKey, country_rapid, countryCode_rapid),
            by = "countryKey") %>% 
  select(-countryKey)

# Save `records` as CSV file
write_csv(records_countryCleanup, 
          paste("rapid-joined-records_country-cleanup_", Sys.Date(), ".csv", sep = ""),
          na = "")