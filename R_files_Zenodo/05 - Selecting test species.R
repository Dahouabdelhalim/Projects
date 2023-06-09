## 05 - Selecting test species

# Libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(XML)
library(xml2)
library(rvest)
library(stringr)

# Read data and make folders ----------------------------------------------

loc.gbifdata <- "./Data/GBIFData/"
# loc.figure <- "./Figures/"

gbif.data <- read.csv(file = paste0(loc.gbifdata, "7_compiled_GBIFData_coordcleanComplete.csv"),
                      stringsAsFactors = FALSE)
flickr.results <- readRDS(file = "flickr_results.Rds")

bboxs <- readRDS(file = "BoundingBoxes.Rds")

# Geographic check loop ---------------------------------------------------

# a loop that checks if a species is within one of the three regions, and then
# verifies if all locations are within the bounding box for that region. Also
# pulls the number of flickr results available for that species if the region
# checks are passed.

ran.species <- unique(gbif.data$species)

f.i <- 0
i.p <- 0
list.potential <- vector(mode = "list", length = length(ran.species))
flickr.potential <- list()
for(sp in ran.species){
  
  samp.data <- gbif.data[gbif.data$species == sp,]
  
  i.p <- i.p + 1
  i <- 0
  lon.results <- NULL
  lat.results <- NULL
  for(ext in bboxs){
    i <- i + 1

    ext[,1] <- ext[,1] + 5
    ext[,2] <- ext[,2] - 5
    
    if(min(samp.data$lon) > ext[1,1] & max(samp.data$lon) < ext[1,2]){
      lon.results <- c(lon.results, "Passed")
    } else {
      lon.results <- c(lon.results, "Failed")
    }  # end of if statement for ext locations
    
    if(min(samp.data$lat) > ext[2,1] & max(samp.data$lat) < ext[2,2]){
      lat.results <- c(lat.results, "Passed")
    } else {
      lat.results <- c(lat.results, "Failed")
    }  # end of if statement for ext locations
    
    if(i == 1){
      if(lon.results[1] == "Passed" & lat.results[1] == "Passed"){
        print("Potential speices - Asia")
        list.potential[[i.p]] <- data.frame("species" = sp, "region" = "Asia")
        {break}
      } else {
        {next}
      } # end of asia check
    } # end of iter 1 check
    if(i == 2){
      if(lon.results[2] == "Passed" & lat.results[2] == "Passed"){
        print("Potential speices - Africa")
        list.potential[[i.p]] <- data.frame("species" = sp, "region" = "Africa")
        {break}
      } else {
        {next}
      } # end of Africa check
    } # end of iter 2 check
    if(i == 3){
      if(lon.results[3] == "Passed" & lat.results[3] == "Passed"){
        print("Potential speices - SAmerica")
        list.potential[[i.p]] <- data.frame("species" = sp, "region" = "SAmerica")
        {break}
      } else {
        print("All Failed")
        results <- "Failed"
        list.potential[[i.p]] <- data.frame("species" = sp, "region" = "Failed")
        {next}
      } # end of SAmerica check
    } # end of iter 3 check
    
  } # end of ext loop
  
  if(!list.potential[[i.p]][2] == "Failed"){
    sp.fresults <- flickr.results[unlist(lapply(flickr.results,
                                 FUN = function(x){x$keyword[1] == sub(" ", "\\\\+", sp)}))]
    sp.fresults <- sp.fresults[-which(sapply(sp.fresults, is.null))]
    
    if(length(sp.fresults) == 0){
      list.potential[[i.p]] <- data.frame(list.potential[[i.p]],
                                          "UnverFlickrRes" = 0,
                                          "SNameFlickrRes" = 0)
    } else {
      sp.fresults <- sp.fresults[[1]]
      
      # if loop that removes records of that species from regions zone
      if(list.potential[[i.p]]$region == "Asia"){
        ext <- bboxs$Asia
        sp.fresults <- sp.fresults[sp.fresults$longitude > ext[1,1] &
                                     sp.fresults$longitude < ext[1,2],]
        sp.fresults <- sp.fresults[sp.fresults$latitude > ext[2,1] &
                                     sp.fresults$latitude < ext[2,2],]
      } else if(list.potential[[i.p]]$region == "Africa"){
        ext <- bboxs$Africa
        sp.fresults <- sp.fresults[sp.fresults$longitude > ext[1,1] &
                                     sp.fresults$longitude < ext[1,2],]
        sp.fresults <- sp.fresults[sp.fresults$latitude > ext[2,1] &
                                     sp.fresults$latitude < ext[2,2],]
      } else if(list.potential[[i.p]]$region == "SAmerica"){
        ext <- bboxs$SAmerica
        sp.fresults <- sp.fresults[sp.fresults$longitude > ext[1,1] &
                                     sp.fresults$longitude < ext[1,2],]
        sp.fresults <- sp.fresults[sp.fresults$latitude > ext[2,1] &
                                     sp.fresults$latitude < ext[2,2],]
      }
      f.i <- f.i + 1
      flickr.potential[[f.i]] <- sp.fresults
      # sp.fresults <- sp.fresults[[1]]
      SNamecount <- sp.fresults[sp.fresults$keyword == sp.fresults$keyword[1],]
      
      list.potential[[i.p]] <- data.frame(list.potential[[i.p]],
                                          "UnverFlickrRes" = dim(sp.fresults)[1],
                                          "SNameFlickrRes" = dim(SNamecount)[1])
    } # end of if for adding flickr results
  } # if statement to pull out flickr results
} # end of loop for sp

saveRDS(list.potential, file = "list_potentials.Rds", compress = FALSE)
saveRDS(flickr.potential, file = "Flickr_potentials.Rds", compress = FALSE)

# Taxonomic stability loop ------------------------------------------------

# testing whether Reptile Database has more than one name for a species after 2000

for(i.p in 1:length(list.potential)){
  
  tryCatch({
    rm(po.spec)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  tryCatch({
    rm(spec)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  
  if(any(list.files() %in% "list_potentials.Rds")){
    list.potential <- readRDS(file = "list_potentials.Rds")
  }
  
  po.spec <- list.potential[[i.p]]
  spec <- po.spec$species
  
  if( !is.null(po.spec$Post2KTaxChange) ){
    print(paste0(spec, " - skipped"))
    {next}
  }
  
  # skip those that have no flickr images
  if( is.null(po.spec$UnverFlickrRes) ){
    # list.potential[[i.p]] <- data.frame(list.potential[[i.p]],
    #                                     "Post2KTaxChange" = NA)
    {next}
  }
  if(po.spec$UnverFlickrRes == 0){
    # list.potential[[i.p]] <- data.frame(list.potential[[i.p]],
    #                                     "Post2KTaxChange" = NA)
    {next}
  }
  
  split.spec <- str_split(spec, pattern = "\\\\ ")[[1]]
  gen <- split.spec[1]
  sp <- split.spec[2]
  url <- paste0("http://reptile-database.reptarium.cz/species?genus=", gen, 
                "&species=", sp)
  h.url <- read_html(url)
  h.nodes <- html_nodes(h.url, xpath = "//td")
  
  if(length(html_nodes(h.url, xpath = "//td")) == 0){
    print(paste0(spec, " - Not found"))
    {next}
  }
  
  html.raw <- htmlTreeParse(h.url, useInternal = TRUE)
  h.text <- xpathApply(html.raw, '//td',
                       function(x) xpathSApply(x,".//text()", xmlValue))[[8]]
  
  years <- str_extract(h.text, "[[:digit:]]+")
  
  syn.post2K <- length(unique(word(h.text[years >= 2000], 1, 2)))
  
  list.potential[[i.p]] <- data.frame(list.potential[[i.p]],
                                      "Post2KTaxChange" = syn.post2K)
  print(paste(spec, "--- post2K names:", syn.post2K))
  saveRDS(list.potential, file = "list_potentials.Rds", compress = FALSE)
  
} # end of potential species loop
saveRDS(list.potential, file = "list_potentials.Rds", compress = FALSE)

list.potential <- readRDS(file = "list_potentials.Rds")

potential.df <- bind_rows(list.potential)
potential.df$SNameFlickrRes

# Random selection of species ---------------------------------------------

set.seed(1)
# selecting those species that are taxonomically stable since 2000
# and pulling a random 9 from the top 25 that returned results from their
# scientific names. Excludes those manually selected.
# Multiple runs required to remove those with tiny distributions or cannot
# be verified by images alone.
random_species <- potential.df %>%
  filter(
    !species == "Ophiophagus hannah" &
      !species == "Malayopython reticulatus" &
      !species == "Coelognathus radiatus" &
      !species == "Calloselasma rhodostoma" &
      !species == "Bungarus fasciatus" &
      !species == "Bothriechis schlegelii" &
      !species == "Eunectes murinus" &
      !species == "Dendroaspis polylepis" &
      !species == "Bitis arietans") %>% 
  filter(SNameFlickrRes > 1, Post2KTaxChange == 1) %>%
  arrange(desc(SNameFlickrRes)) %>%
  top_n(n = 25, SNameFlickrRes) %>%
  sample_n(size = 9)

write.csv(x = random_species, file = "./Data/Random_9species.csv", row.names = FALSE)

random_species2 <- potential.df %>%
  filter(
    !species == "Ophiophagus hannah" &
      !species == "Malayopython reticulatus" &
      !species == "Coelognathus radiatus" &
      !species == "Calloselasma rhodostoma" &
      !species == "Bungarus fasciatus" &
      !species == "Bothriechis schlegelii" &
      !species == "Eunectes murinus" &
      !species == "Dendroaspis polylepis" &
      !species == "Bitis arietans" &
      !species %in% random_species$species
      ) %>% 
  filter(SNameFlickrRes > 1, Post2KTaxChange == 1) %>%
  arrange(desc(SNameFlickrRes)) %>%
  top_n(n = 25, SNameFlickrRes) %>%
  sample_n(size = 4)

random_species2

write.csv(x = random_species2, file = "./Data/Random2_4species.csv", row.names = FALSE)

random_species3 <- potential.df %>%
  filter(
    !species == "Ophiophagus hannah" &
      !species == "Malayopython reticulatus" &
      !species == "Coelognathus radiatus" &
      !species == "Calloselasma rhodostoma" &
      !species == "Bungarus fasciatus" &
      !species == "Bothriechis schlegelii" &
      !species == "Eunectes murinus" &
      !species == "Dendroaspis polylepis" &
      !species == "Bitis arietans" &
      !species %in% random_species$species &
      !species == "Porthidium lansbergii" &
      !species == "Porthidium nasutum" 
      ) %>% 
  filter(SNameFlickrRes > 1, Post2KTaxChange == 1) %>%
  arrange(desc(SNameFlickrRes)) %>%
  top_n(n = 25, SNameFlickrRes) %>%
  sample_n(size = 2)

random_species3

write.csv(x = random_species3, file = "./Data/Random3_2species.csv", row.names = FALSE)

random_species4 <- potential.df %>%
  filter(
    !species == "Ophiophagus hannah" &
      !species == "Malayopython reticulatus" &
      !species == "Coelognathus radiatus" &
      !species == "Calloselasma rhodostoma" &
      !species == "Bungarus fasciatus" &
      !species == "Bothriechis schlegelii" &
      !species == "Eunectes murinus" &
      !species == "Dendroaspis polylepis" &
      !species == "Bitis arietans" &
      !species %in% random_species$species &
      !species == "Porthidium lansbergii" &
      !species == "Porthidium nasutum" &
      !species == "Echis khosatzkii" &
      !species %in% random_species2$species
      ) %>% 
  filter(SNameFlickrRes > 1, Post2KTaxChange == 1) %>%
  arrange(desc(SNameFlickrRes)) %>%
  top_n(n = 25, SNameFlickrRes) %>%
  sample_n(size = 1)

random_species4

write.csv(x = random_species4, file = "./Data/Random4_1species.csv", row.names = FALSE)


# Extract the Flickr data for manual review -------------------------------

# repdb.names <- readRDS(file = "repdb_commnames.Rds")
# random_species <- read.csv(file = "./Data/Random_9species.csv", stringsAsFactors = FALSE)
# random_species <- read.csv(file = "./Data/Random2_4species.csv", stringsAsFactors = FALSE)
# random_species <- read.csv(file = "./Data/Random3_2species.csv", stringsAsFactors = FALSE)
random_species <- read.csv(file = "./Data/Random4_1species.csv", stringsAsFactors = FALSE)

flickr.potential <- readRDS(file = "Flickr_potentials.Rds")
loc.gbifdata <- "./Data/GBIFData/"
gbif.data <- read.csv(file = paste0(loc.gbifdata, "7_compiled_GBIFData_coordcleanComplete.csv"),
                      stringsAsFactors = FALSE)

dir.create("./Data/RandomSpeciesData/")
for(sp in random_species$species){
  sp.data <- gbif.data[gbif.data$species == sp,]
  write.csv(file = paste0("./Data/RandomSpeciesData/", sp, " gbif.csv"), x = sp.data, row.names = FALSE)
  
  sp.fresults <- flickr.potential[unlist(lapply(flickr.potential,
                                              FUN = function(x){x$keyword[1] == sub(" ", "\\\\+", sp)}))]
  sp.fresults <- sp.fresults[-which(sapply(sp.fresults, is.null))]
  
  write.csv(file = paste0("./Data/RandomSpeciesData/", sp, " flickr.csv"),
            x = sp.fresults, row.names = FALSE)
}

