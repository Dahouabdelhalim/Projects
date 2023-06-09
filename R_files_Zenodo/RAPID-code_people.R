# RAPID code for fetching data from Bionomia to reintegrate with data in BIOSPEX
# This code accompanies document 'RAPID-protocol_people.pdf'
# Code written by Katelin Pearson, 2020-12-16. Updated by Erica Krimmel on
# 2021-02-02.

# Load core libraries; install these packages if you have not already
library(tidyverse)
library(jsonlite)

# Read into R the data exported from Bionomia and specify column types (doing
# this is recommended when dataset is large and columns may be sparsely
# populated). Change filename for reuse.
specs <- read_csv("biospex_people_2020-12-15.csv",
                  col_types = cols('_id' = col_character(),
                                   gbifID = col_double(),
                                   recordedBy = col_character(),
                                   recordedByID = col_character(),
                                   identifiedBy = col_character(),
                                   identifiedByID = col_character()))

# Modify `specs` for data enhancement process
specs <- specs %>% 
  # Create column to use when fetching data from Bionomia
  mutate(gbifIDurl = paste0("https://bionomia.net/occurrence/",
                            specs$gbifID,
                            ".json")) %>% 
  # Create columns for the enhanced people data from Bionomia
  mutate(recordedBy_rapid = NA) %>% 
  mutate(recordedByID_rapid = NA) %>% 
  mutate(identifiedBy_rapid = NA) %>% 
  mutate(identifiedByID_rapid = NA)

# Define function to apply a NULL value to any JSON queries returning an error
fromJSON_possibly <- possibly(fromJSON,
                              otherwise = NULL)

# Loop through records in 'spec' to fetch person data from Bionomia
for(i in 1:dim(specs)[1]){
  print(i)
  # Skip record if it lacks a value for `gbifID`
  if(is.na(specs$gbifID[i])){
    next
  }
  # Use record to fetch data from Bionomia
  thisSpec <- fromJSON_possibly(specs$gbifIDurl[i],
                                simplifyDataFrame = TRUE,
                                flatten = TRUE)
  # Skip record if it lacks Bionomia data
  if(is.null(thisSpec)){
    next
  }else{
    # Check for collector person data from Bionomia
    if(length(thisSpec$recorded) > 0){
      # Create table of collector person data
      recby <- as_tibble(thisSpec$recorded)
      
      # Check for multiple collectors. If there are multiple collectors,
      # concatenate values that will be transferred to `recordedBy_rapid` and
      # `recordedByID_rapid`
      if(dim(recby)[1]>1){
        # Make placeholders for concatenated values
        allids <- c()
        allnam <- c()
        # Concatenate values
        for(k in 1:dim(recby)[1]){
          allids <- paste(allids,
                          recby$`@id`[k],
                          sep=" | ")
          allnam <- paste(allnam,
                          paste0(recby$givenName[k],
                                " ",
                                recby$familyName[k]),
                          sep=" | ")
        }
        
        # Remove extra pipes from strings
        allids <- sub(' \\\\| ', '', allids)
        allnam <- sub(' \\\\| ', '', allnam)
        
        # Put concatenated values for `id` into `recordedByID_rapid`
        specs$recordedByID_rapid[i] <- allids
       
        # Put concatenated values for `recorded` into `recordedBy_rapid`
        specs$recordedBy_rapid[i] <- allnam
        
      # If there is a single collectors, transfer values directly into
      # `recordedBy_rapid` and `recordedByID_rapid`
      }else{
        # Put value for `id` in `recordedByID_rapid`
        specs$recordedByID_rapid[i] <- recby$`@id`
        # Put value for `recorded` into `recordedBy_rapid`
        specs$recordedBy_rapid[i] <- paste0(recby$givenName,
                                           " ",
                                           recby$familyName)
      }
    }
    
    # Check for identifier person data from Bionomia
    if(length(thisSpec$identified) > 0){
      # Create table of identifier person data
      idby <- as_tibble(thisSpec$identified)
      
      # Check for multiple identifiers. If there are multiple identifiers,
      # concatenate values that will be transferred to `identifiedBy_rapid` and
      # `identifiedByID_rapid`
      if(dim(idby)[1]>1){
        # Make placeholders for concatenated values
        allids <- c()
        allnam <- c()
        # Concatenate values
        for(k in 1:dim(idby)[1]){
          allids <- paste(allids,
                          idby$`@id`[k],
                          sep=" | ")
          allnam <- paste(allnam,
                          paste0(idby$givenName[k],
                                " ",
                                idby$familyName[k]),
                          sep=" | ")
        }

        # Remove the extra pipes from strings
        allids <- sub(' \\\\| ', '', allids)
        allnam <- sub(' \\\\| ', '', allnam)
        
        # Put concatenated values for `id` into `identifiedByID_rapid`
        specs$identifiedByID_rapid[i] <- allids
        
        # Put concatenated values for `recorded` into `identifiedBy_rapid`
        specs$identifiedBy_rapid[i] <- allnam
      }else{
        # Put value for `id` in `identifiedByID_rapid`
        specs$identifiedByID_rapid[i] <- idby$`@id`
        # Put value for `recorded` in `identifiedBy_rapid`
        specs$identifiedBy_rapid[i] <- paste0(idby$givenName,
                                              " ",
                                              idby$familyName)
      }
    }
  }
}

# Save data as a CSV file to import back into BIOSPEX
write_csv(specs,
          paste0("RAPID-people_",Sys.Date(),".csv"),
          na = "")
