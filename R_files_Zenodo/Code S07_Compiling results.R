
library(dplyr)
library(stringr)

# Snapshot online data ----------------------------------------------------

extractFiles <- list.files("./Data/WebData", pattern = "KEYWORD_EXTRACT",
                           recursive = TRUE, full.names = TRUE)

extractData <- do.call(rbind, lapply(extractFiles, function(x){
  df <- read.csv(x, stringsAsFactors = FALSE)
  if(dim(df)[2] > 1){
    return(df)
  }
}))

# checked the text along with some genera to see if they were really detected,
# IE text makes it clear they are not talkig about spp, but it's from elsewhere
# in the page

# genToCheck <- unique(generaWebRaw$keyw)[!unique(generaWebRaw$keyw) %in% 
#                                           c(
#                                             word(unique(speciesWebRaw$keyw),1,1),
#                                             word(unique(speciesWebRaw$sp),1,1)
#                                           )
#                                         & nchar(unique(generaWebRaw$keyw)) < 6]
# write.csv(x = generaWebRaw[generaWebRaw$keyw %in% genToCheck,], 
#         file = "./Data/Snapshot Online Genera Data CHECK.csv",
#           row.names = FALSE)

# clean up instances of duff genera and species
generaWebRaw <- extractData %>% 
  filter(spORgen == "GENUS", !termsSurrounding == "", !sp == "species",
         !keyw %in% c("nomen", "rufus", "Dia", "Diana", "Mala", 
                      "Inca", "Pero", "May", "Janus", "Yukon", "Lucia",
                      "Zora", "Beata", "Neon", "Prima", "Meta", "Patri",
                      "Enna", "Maso", "Mica", "Perro")) # remove spurious genera

unique(generaWebRaw$keyw)[nchar(unique(generaWebRaw$keyw)) < 6]
length(unique(generaWebRaw$sp))

write.csv(x = generaWebRaw, file = "./Data/Snapshot Online Genera Data Raw.csv",
          row.names = FALSE)

speciesWebRaw <- extractData %>% 
  filter(spORgen == "SPECIES", str_detect(str_trim(keyw), " "),
         !str_detect(keyw, "nomen dubium"))

length(unique(speciesWebRaw$sp))
length(unique(speciesWebRaw$keyw))
length(unique(word(speciesWebRaw$sp, 1, 1)))

# actually searched
length(unique(extractData$webID))
# had species
length(unique(speciesWebRaw$webID))

speciesWebRaw %>% 
  group_by(sp, webID) %>% 
  slice(n = 1) %>% 
  ungroup() %>% 
  count(webID) %>% 
  arrange(desc(n)) %>% 
  summarise(mean(n),
            se = sqrt(var(n)/length(n)),
            min(n),
            max(n))

write.csv(x = speciesWebRaw, file = "./Data/Snapshot Online Species Data Raw.csv",
          row.names = FALSE)

# Temporal online data ----------------------------------------------------

extractFilesTemp <- list.files("./Data/TemporalData", pattern = "KEYWORD_EXTRACT",
                               recursive = TRUE, full.names = TRUE)

extractDataTemp <- do.call(rbind, lapply(extractFilesTemp, function(x){
  df <- read.csv(x, stringsAsFactors = FALSE)
  if(dim(df)[2] > 1){
    return(df)
  }
}))

# clean up instances of duff genera and species
generaTempWebRaw <- extractDataTemp %>% 
  filter(spORgen == "GENUS", !termsSurrounding == "", !sp == "species",
         !keyw %in% c("nomen", "rufus", "Dia", "Diana", "Mala", 
                      "Inca", "Pero", "May", "Janus", "Yukon", "Lucia",
                      "Zora", "Beata", "Neon", "Prima", "Meta", "Patri",
                      "Enna", "Maso", "Mica", "Perro")) # remove spurious genera

speciesTempWebRaw <- extractDataTemp %>% 
  filter(spORgen == "SPECIES", str_detect(str_trim(keyw), " "),
         !str_detect(keyw, "nomen dubium"))

# Add online results to master data ---------------------------------

araData <- read.csv("./Data/arachnid_species_data_syns.csv")

# make sure hte non spiders have something for the keywords
araData$allNames[!araData$clade == "Spider"] <- araData$accName[!araData$clade == "Spider"]
araData$allGenera[!araData$clade == "Spider"] <- word(araData$accName[!araData$clade == "Spider"], 1, 1)

# drop any duplicated species that remain
araData <-
  araData[!duplicated(araData$accName),]

araData_trade <- araData %>% 
  mutate(
    onlineTradeSnap = accName %in% unique(speciesWebRaw$keyw),
    onlineTradeSnap_Any = str_detect(allNames, paste0(unique(speciesWebRaw$keyw), collapse = "|")),
    onlineTradeSnap_genus = genus %in% unique(generaWebRaw$keyw),
    onlineTradeSnap_genusAny = str_detect(allGenera, paste0(unique(generaWebRaw$keyw), collapse = "|")),
    onlineTradeTemp = accName %in% unique(speciesTempWebRaw$keyw),
    onlineTradeTemp_Any = str_detect(allNames, paste0(unique(speciesTempWebRaw$keyw), collapse = "|")),
    onlineTradeTemp_genus = genus %in% unique(generaTempWebRaw$keyw),
    onlineTradeTemp_genusAny = str_detect(allGenera, paste0(unique(generaTempWebRaw$keyw), collapse = "|")),
    onlineTradeEither = ifelse(onlineTradeSnap | onlineTradeTemp, TRUE, FALSE),
    onlineTradeEither_Any = ifelse(onlineTradeSnap_Any | onlineTradeTemp_Any, TRUE, FALSE)
  )

# only the genera that are currently accepted
araData_trade %>% 
  filter(onlineTradeSnap_genus) %>% 
  summarise(length(unique(genus)))

# largest possible range of online traded genera
araData_trade %>% 
  filter(onlineTradeSnap_genusAny) %>% 
  summarise(length(unique(genus)))

# Add LEMIS data ----------------------------------------------------------

lemisAraData <- read.csv(file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv",
                         stringsAsFactors = FALSE)

lemisSpp <- lemisAraData %>%
  mutate(lemisName = str_to_sentence(paste(genus, species))) %>%
  filter(!str_detect(lemisName, "sp\\\\.$") & !lemisName == "NA NA") %>% 
  pull(lemisName)

araData_trade <- araData_trade %>% 
  mutate(
    LEMIStrade = accName %in% unique(lemisSpp),
    LEMIStrade_Any = str_detect(allNames, paste0(unique(lemisSpp), collapse = "|")),
    LEMIStrade_genus = genus %in% unique(lemisAraData$genus),
    LEMIStrade_genusAny = str_detect(allGenera, paste0(unique(lemisAraData$genus), collapse = "|"))
  )

# CITES trade database ----------------------------------------------------

citesData <- read.csv(file = "./Data/TradeDatabases/CITES_gross_imports_2021-09-15_comma_separated.csv",
                      stringsAsFactors = FALSE)

citesGenera <- citesData %>% 
  mutate(genus = word(Taxon, 1, 1)) %>% 
  pull(genus)

araData_trade <- araData_trade %>% 
  mutate(
    CITEStrade = accName %in% unique(citesData$Taxon),
    CITEStrade_Any = str_detect(allNames, paste0(unique(citesData$Taxon), collapse = "|")),
    CITEStrade_genus = genus %in% unique(citesGenera),
    CITEStrade_genusAny = str_detect(allGenera, paste0(unique(citesGenera), collapse = "|"))
  )


# Add CITES appendices ---------------------------------------------------

citesApp <- read.csv("./Data/TradeDatabases/CITES_app_2021-09.csv")

citesApp$FullName

araData_trade <- araData_trade %>% 
  mutate(
    CITESapp = ifelse(accName %in% citesApp$FullName, "II", "Not listed"),
    CITESapp_Any = ifelse(str_detect(allNames,
                                     paste0(unique(citesApp$FullName),
                                            collapse = "|")), "II", "Not listed" ),
  )

sum(araData_trade$CITESapp == "II")
sum(araData_trade$CITESapp_Any == "II")

# 7 instances of names in CITES considered dubious names by WSC at time of WSC data download
citesApp$FullName[!citesApp$FullName %in% unlist(str_split(araData_trade$allNames, ";"))]

# Add IUCN redlist --------------------------------------------------------

redlistData <- read.csv("./Data/TradeDatabases/IUCN_redlist_assessments_2021-09-15.csv")

redlistData$scientificName
redlistData$redlistCategory

araData_trade$redlist <- "Not assessed"
araData_trade$redlist_Any <- "Not assessed"
for(row in 1:dim(redlistData)[1]){
  # extact matches
  araData_trade$redlist[araData_trade$accName ==
                          redlistData$scientificName[row]] <-
    redlistData$redlistCategory[row]
  # any match in all names
  araData_trade$redlist_Any[str_detect(araData_trade$allNames, redlistData$scientificName[row])] <-
    redlistData$redlistCategory[row]
}

table(araData_trade$redlist)
table(araData_trade$redlist_Any)
sum(!araData_trade$redlist == "Not assessed")
sum(!araData_trade$redlist_Any == "Not assessed")

# Final column for traded in any source -----------------------------------

araData_trade <- araData_trade %>% 
  group_by(accName) %>% 
  mutate(extactMatchTraded = ifelse(onlineTradeSnap |
                                      onlineTradeTemp|
                                      LEMIStrade|
                                      CITEStrade, TRUE, FALSE),
         anyMatchTraded = ifelse(onlineTradeSnap_Any|
                                   onlineTradeTemp_Any|
                                   LEMIStrade_Any|
                                   CITEStrade_Any, TRUE, FALSE)) %>% 
  ungroup()

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# species counts from online trade snapshot data
araData_trade %>% 
  summarise(across(c("onlineTradeSnap", "onlineTradeSnap_Any",
                     "onlineTradeTemp" , "onlineTradeTemp_Any", 
                     "LEMIStrade", "LEMIStrade_Any",
                     "CITEStrade", "CITEStrade_Any"),
                   ~ sum(.x, na.rm = TRUE)))

araData_trade %>% 
  summarise(across(c("onlineTradeEither" , "onlineTradeEither_Any"),
                   ~ sum(.x, na.rm = TRUE)))

araData_trade %>% 
  summarise(across(c("extactMatchTraded", "anyMatchTraded"),
                   ~ sum(.x, na.rm = TRUE)))

araData_trade %>% 
  group_by(clade) %>% 
  summarise(across(c("extactMatchTraded", "anyMatchTraded"),
                   ~ sum(.x, na.rm = TRUE)),
            nSpp = n()
            ) %>% 
  mutate(per = anyMatchTraded / nSpp *100)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

write.csv(araData_trade, "./Data/arachnid_species_data_traded.csv", row.names = FALSE)
