
dir.create("./Data")
dir.create("./Figures")

library(rvest)
library(dplyr)
library(stringr)

araData <- read.csv("./Data/arachnid_species_data.csv")
araData$accName <- paste(araData$genus, araData$species)

# Function to get names from WSC ------------------------------------------

wsc_get_names <- function(speciesID){
  
  # url <- "https://wsc.nmbe.ch/species/22720"
  url <- paste0("https://wsc.nmbe.ch/species/", speciesID)
  
  htmlData <- read_html(url,
                        encoding = "UTF-8")
  
  htmlData %>% 
    html_name()
  
  htmlRefs <- htmlData %>% 
    html_children() %>% 
    html_nodes(".paragraph") %>% 
    magrittr::extract2(2) %>% 
    html_children() %>% 
    magrittr::extract2(3) %>% 
    html_children()
  
  htmlRefsText <- htmlRefs %>% 
    html_text()
  
  spnames <- unique(htmlRefsText[htmlRefs %>% 
                                   html_name() == "i"])
  # make sure binomial only
  spnames <- spnames[str_detect(spnames, "\\\\s")]
  return(spnames)
}

# Feeding in all species ID to get other names ----------------------------

keywordList <- vector("list", length = length(araData$speciesId))
names(keywordList) <- araData$accName

# 23890
# keywordList[[40000]]

for(sp in araData$accName){
  
  if("KeywordList.rds" %in% list.files("./Data")){
    keywordList <- readRDS("./Data/KeywordList.rds")
  }
  # sp <- araData$accName[51000]
  
  spCode <- araData$speciesId[araData$accName == sp]
  
  print(paste0("--- ", sp, " --- ", which(araData$accName %in% sp)))
  if(is.null(keywordList[[sp]]) & !is.na(spCode)){
    spsyns <- unlist(lapply(spCode, function(x){
      wsc_get_names(x)
    }))
    keywordList[[sp]] <- c(sp, spsyns)
    Sys.sleep(1)
    print("Retrieved")
    saveRDS(keywordList, file = "./Data/KeywordList.rds")
  } else {
    # keywordList[[sp]] <- sp
    print("Done")
  }
  # some sort of break to reduce load server, but with 49K queries can't be too
  # high
}
saveRDS(keywordList, file = "./Data/KeywordList.rds")

# BASIC LIST END ----------------------------------------------------------

i <- 0
printI <- seq(1, 51950, 50)
araData$allNames <- NA
araData$allGenera <- NA
for(iName in names(keywordList)){
  # iName <- "Pararaneus pseudostriatus"
  # iName <- "Calilena stylophora"
  # iName <- "Emblyna annulipes"
  
  if(is.null(keywordList[iName][[1]])){
    {next}
  }
  
  i <- i+1
  if(i %in% printI){
    print(i)
  }
  
  # first make a complete list of genera
  allGenera <- unique(c(araData[araData$accName == iName,]$genus[1],
                        word(keywordList[iName][[1]],1,1)))
  # and remove those with abbreviations
  allGenera <- allGenera[!str_detect(allGenera, "\\\\.")]
  
  # before adding species all togther, we need to fix a few instances of genus
  # abbreviations, it is replaced with the gennus first matching the first
  # letter, newest/current genus is prioitisedS
  
  # compile all species
  allSpecies <- unique(c(iName, keywordList[iName][[1]]))
  # and ditch those that aren't actually to the species level
  allSpecies <- allSpecies[!str_detect(allSpecies, "sp\\\\.")]
  
  # check for abbreications and replace
  if(length(allSpecies[str_detect(allSpecies, "\\\\.")]) > 0){
    abbrivatedSpp <- keywordList[iName][[1]][str_detect(keywordList[iName][[1]], "\\\\.")]
    for(abSpp in abbrivatedSpp){
      # abSpp <- abbrivatedSpp[1]
      
      firstLspp <- substr(abSpp,1,1)
      # get the first genus name that matches the abbreivation
      replaceGenera <- allGenera[grep(paste0("^", firstLspp, ""), allGenera)][1]
      
      allSpecies[allSpecies == abSpp] <- sub(paste0(firstLspp, "\\\\."), replaceGenera, abSpp)
      # then replace in case it has made a dup
      allSpecies <- unique(allSpecies)
    }
    
  }
  # add compiled genera to a new column
  araData[araData$accName == iName,]$allGenera <- 
    paste(allGenera, collapse = ";")
  
  # add synonyms to the original datafram split by ;
  araData[araData$accName == iName,]$allNames <- 
    paste(allSpecies, collapse = ";")
  
}

# we still have an issue of duplicated species because of subspecies stuff
araData <- araData[!duplicated(araData$accName),]
# but we can straight remove them as the subpsecies info was not used to
# discover synonyms

write.csv(araData, "./Data/arachnid_species_data_syns.csv", row.names = FALSE)
