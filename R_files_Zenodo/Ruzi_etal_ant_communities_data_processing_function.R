## function to help process and transform the raw data into data files
## that can be used in subsequent analyses

SBS_transform <- function(x, row_ID = TRUE, trap_ID, spp_ID, 
                          rem_spp = TRUE, freq = FALSE, warn = TRUE) {
  #To transform a tibble of raw incidence data in the form of
  #samples (rows) by species (columns) to raw incidence data
  #in the form of species (rows) by samples(columns). 
  
  #x is the datafile, it will need to be a class tibble
  
  #row_ID = deafulat to TRUE, will convert row numbers to a column and
  #will include this in the columns that it will unite together to form
  #a unique trap information column 
  
  #trap_ID is a vector of characters of which columns to unite to form trap information
  
  #spp_ID is a vector of characters of which columns are the different species
  
  #rem_spp = defaults to TRUE, will remove rows where the species was never
  #found at any of the traps 
  
  #freq = defaults to FALSE, therefore raw incidence tibble will be output,
  #if set to TRUE then instead of a raw incidence tibble being output, the 
  #raw incidences will be output as frequencies, with the number of traps 
  #being the first number and the rest being the number of traps that that
  #species was found in
  
  #warn = defaults to TRUE, will warn the user that x is not a tibble
  
  if(warn!=TRUE && warn!=FALSE) stop("warn must be set to TRUE or FALSE")
  if(row_ID!=TRUE && row_ID!=FALSE) stop("row_ID must be set to TRUE or FALSE ")
  if(rem_spp!=TRUE && rem_spp!=FALSE) stop("rem_spp must be set to TRUE or FALSE")
  
  
  #check to make sure that the input is a tibble and if not make it one
  if(warn == TRUE){ #if warn is on and not a tibble then warn the user
    if(is_tibble(x) == FALSE) stop("x must be a tibble")
    if(is.vector(trap_ID) == FALSE) stop("trap_ID must be a vector")
  } else { # if warn is not equal to true and not a tibble, then make it a tibble
    if(is_tibble(x)==FALSE) x <- as_tibble(x)
  }
  
  if(row_ID == TRUE){
    x2 <- x %>%
      rowid_to_column("row_ID") %>% # convert row numbers to a column
      unite(trap_ID, c(row_ID, trap_ID), sep = ".")
  } else {
    x2 <- x %>%
      unite(trap_ID, trap_ID, sep = ".")
  }
  
  x3 <- x2 %>%
    # to bring all the species into one column
    gather(spp_ID, key = "spp", value = "presence", factor_key = TRUE) %>%
    # to spread out the trap IDs across columns
    spread(key = trap_ID, value = presence)
  
  # remove species that were never found in any of the traps
  if(rem_spp == TRUE) {
    x3 <- x3 %>%
      # sum up the frequency across species
      mutate(sum = rowSums(.[names(.) != "spp"])) %>%
      # remove the spp that are not found
      filter(sum > 0) %>%
      # remove the added sum column
      select(-sum)
  }
  
  if(freq != FALSE) {
    # to calculate the number of traps minus the spp column
    # this will not work if there are additional columns in
    # the tibble that are unrelated to these two things
    trap_num <- ncol(x3) - 1
    # need to sum up the frequences of the species in each of 
    # the traps and specify not to look at the spp column
    x4 <- x3 %>%
      mutate(sum = rowSums(.[names(.) != "spp"]))
    spp_freq <- as.numeric(c(trap_num, x4$sum))
    return(spp_freq)
  } else {
    return(x3)
  }
}
