## 01 - Species list, GBIF download, GBIF clean, keyword list

# Libraries ---------------------------------------------------------------

library(taxize)
library(dismo)
library(rgbif)
library(data.table)
library(stringr)
library(ggplot2)
library(CoordinateCleaner)
# for reptileDB pulls
library(XML)
library(xml2)
library(rvest)

# Folders and et-up api key set-up ----------------------------------------

loc.gbifdata <- "./Data/GBIFData/"
dir.create(loc.gbifdata)

# set up NCBI api key
# use_entrez()
options(ENTREZ_KEY = "[NCBI API KEY HERE]")
getOption("ENTREZ_KEY")
getkey(x = "[NCBI API KEY HERE]", service = "ncbi")


# Producing list of snake familes -----------------------------------------

# get the gbif code for squamata
e.id <- get_gbifid("Squamata")
gbif.id <- e.id[1]

# trim the squamata gbif results to snakes only, gbif does not have serpentes as
# a level so we use NCBI
gbif.fam <- downstream("715", downto = "Family", db = "gbif", limit = 1000)
ncbi.names <- classification(gbif.fam$`715`[,1], db = 'ncbi')
ncbi.list <- (as.list(ncbi.names))
ncbi.list <- ncbi.list[!is.na(ncbi.list)]
ncbi.snakes <- ncbi.list[str_detect(ncbi.list, "Serpentes")]
# store the final family names
snake.families <- names(ncbi.snakes)


# GBIF data download loop -------------------------------------------------

f.i <- 1
fam.i <- 1
fam.fails <- list()
gen.fails <- list()
for(fam in snake.families){

  print(paste0("--- ", fam))
  try(rm(genera))
  
  # retrieve the GBIF code for the specific family
  e.id <- get_gbifid(fam)
  gbif.id <- e.id[1]
  
  # attempt to pull genera within that family
  try(
    genera <- gbif_downstream(gbif.id, downto = "Genus")
  )
  
  # attempt to pull genera from NCBI as well if GBIF code fails, and tidy if exists
  if(!exists("genera")){
    try(
      genera <- downstream(fam, downto = "Genus", db = "ncbi")
    )
    try(
      genera <- genera[[1]]
    )
    try(
    names(genera) <- c("id", "name", "rank")
    )
  }
  
  # produce a list of failed familes
  if(!exists("genera")){
    fam.fails[[fam.i]] <- fam
    fam.i <- fam.i + 1
    next}
  
  genera <- genera[!duplicated(genera$name),]
  
  # loop to run through each genera and download the GBIF data
  for(gen in genera$name){
    print(paste0("------ ", gen))
    
    # this is a check to skip genera that have already been downloaded and
    # stored, allows return to the download loop without losing progress
    if(paste0(gen, "_rawGBIFData.csv") %in% list.files(loc.gbifdata)){
      print(paste0("------ ", gen, " already present"))
      next}
    
    try(
      e.data <- gbif(genus = gen, geo = TRUE)
    )
    
    # record genera that failed to download
    if(is.null(e.data)){
      gen.fails[[f.i]] <- gen
      f.i <- f.i + 1
    }
    
    # if the download is successful save that genus' data
    if(!is.null(e.data)){
      write.csv(file = paste0(loc.gbifdata, gen, "_rawGBIFData.csv"), x = e.data)
    } # end of if data loop 
  } # end of genera loop
} # end of snake family loop

fam.fails
gen.fails

gbif.files <- list.files(loc.gbifdata, pattern = "rawGBIFData.csv")


# GBIF download check loop ------------------------------------------------

# Outages in internet can cause partial downloads of GBIF data, this loop checks
# that the files produced match the number of records on GBIF.

list.files <- list()
i <- 0
c.i <- 0
p.i <- 0
comp.dl <- list()
part.dl <- list()
for(file in gbif.files){
  i <- i + 1

  g.data <- fread(file = paste0(loc.gbifdata, file))
  
  gen <- str_split(file, "_")[[1]][1]
  num.rec <- gbif(genus = gen, download = FALSE, geo = TRUE)
  
  if(nrow(g.data) == num.rec){
    c.i <- c.i + 1
    print(paste(gen, "--- Complete Download"))
    comp.dl[[c.i]] <- gen
  } else if(!nrow(g.data) == num.rec){
    p.i <- p.i + 1
    print(paste(gen, "--- Partial Download"))
    part.dl[[p.i]] <- gen
  } # end of if
} # end of file loop
part.dl

partiall.dl.gen <- Reduce(c , part.dl)

# Optional loop to retrieve the missing data only.
for(gen in partiall.dl.gen){
  print(paste0("------ ", gen))
  try(e.data <- gbif(genus = gen, geo = TRUE))
  if(is.null(e.data)){
    gen.fails[[f.i]] <- gen
    f.i <- f.i + 1
  } # end of if loop
  if(!is.null(e.data)){
    write.csv(file = paste0(loc.gbifdata, gen, "_rawGBIFData.csv"), x = e.data)
  } # end of if data loop
} # end of gen loop


# Compiling and cleaning GBIF datasets ------------------------------------

gbif.files <- list.files(loc.gbifdata, pattern = "rawGBIFData.csv")

# read in completed files and rbind all
list.files <- vector(mode = "list", length = length(gbif.files))
i <- 0
for(file in gbif.files){
  i <- i + 1
  print(paste(i, file))
  g.data <- fread(file = paste0(loc.gbifdata, file))
  list.files[[i]] <- g.data
}
gbif.data <- rbindlist(list.files, fill = TRUE)
rm(list.files)

write.csv(file = paste0(loc.gbifdata, "1_compiled_GBIFData.csv"), row.names = FALSE,
          x = gbif.data)

# a trim of the data removing all fields that we do not require for analysis etc.
gbif.data.trim <- gbif.data[,c(2,6,12,16,17,19,21,22,28,35,38,52,61,64,66)]
head(gbif.data.trim)

write.csv(file = paste0(loc.gbifdata, "2_compiled_GBIFData_trimmed.csv"), row.names = FALSE,
          x = gbif.data.trim)

gbif.data.trim <- read.csv(file = paste0(loc.gbifdata, "2_compiled_GBIFData_trimmed.csv"),
                           stringsAsFactors = FALSE)

# filter out marine snakes, due to difficulties verifying records automatically
names(gbif.data.trim)
gbif.data.terr <- dplyr::filter(gbif.data.trim, !genus %in% c("Acalyptophis", "Aipysurus", "Antaioserpens",
                                            "Astrotia", "Emydocephalus", "Enhydrina",
                                            "Ephalophis", "Hydrelaps", "Hydrophis", 
                                            "Kerilia", "Kolpophis", "Lapemis", 
                                            "Laticauda", "Parahydrophis", "Pelamis",
                                            "Praescutata", "Thalassophis"))

write.csv(file = paste0(loc.gbifdata, "3_compiled_GBIFData_terrestrial.csv"), row.names = FALSE,
          x = gbif.data.terr)

gbif.data.terr <- gbif.data.terr[!is.na(gbif.data.terr$specificEpithet),]
gbif.data.terr <- gbif.data.terr[!gbif.data.terr$lon == 0 & !gbif.data.terr$lat == 0,]
gbif.data.terr <- gbif.data.terr[!is.na(gbif.data.terr$lon) & !is.na(gbif.data.terr$lat),]

write.csv(file = paste0(loc.gbifdata, "4_compiled_GBIFData_basicclean.csv"), row.names = FALSE,
          x = gbif.data.terr)

gbif.data <- gbif.data.terr

gbif.data <- read.csv(file = paste0(loc.gbifdata, "4_compiled_GBIFData_basicclean.csv"),
                      stringsAsFactors = FALSE)

gbif.data.clean <- cc_equ(x = gbif.data,
                          lon = "lon",
                          lat = "lat")
gbif.data.clean <- cc_zero(x = gbif.data.clean,
                           lon = "lon",
                           lat = "lat")
gbif.data.clean <- cc_val(x = gbif.data.clean,
                          lon = "lon",
                          lat = "lat")
gbif.data.clean <- cc_gbif(x = gbif.data.clean,
                           lon = "lon",
                           lat = "lat")
gbif.data.clean <- cc_inst(x = gbif.data.clean,
                           lon = "lon",
                           lat = "lat")

write.csv(file = paste0(loc.gbifdata, "5_compiled_GBIFData_coordclean1.csv"), row.names = FALSE,
          x = gbif.data.clean)


# Removing ocean points and outliers --------------------------------------

gbif.data.clean <- read.csv(file = paste0(loc.gbifdata, "5_compiled_GBIFData_coordclean1.csv"), 
                            stringsAsFactors = FALSE)

## After issues with species with more than 15000 records failing to remove
## outliers we created manually defined bounding boxes to remove the outlying
## points.
# table(gbif.data.clean$species)[table(gbif.data.clean$species) > 15000]
## Natrix natrix
# x_coord <- c(-25, -25, 100, 100)
# y_coord <- c(25, 80, 80, 25)
## Natrix maura and Natrix helvetica
# x_coord <- c(-25, -25, 30, 30) 
# y_coord <- c(25, 80, 80, 25) 
## Vipera berus
# x_coord <- c(-25, -25, 170, 170)
# y_coord <-c(25, 80, 80, 25) 
## Thamnophis elegans and sirtalis
# x_coord <- c(-140, -140, -40, -40)
# y_coord <- c(10, 75, 75, 10)

i <- 0
i.i <- 0
sea.clean.list <- list()
removed.list <- list()
issue.list <- list()
num.sp <- length(unique(gbif.data.clean$species))
species <- unique(gbif.data.clean$species)
for(sp in species){
  i <- i + 1
  print(paste0(i, "/", num.sp, " - ", sp))
  
  gc()
  try(rm(sea.test))
  try(rm(sea.test.clean))
  try(rm(sea.out.test.clean))
  
  sea.test <- gbif.data.clean[gbif.data.clean$species == sp,]
  sea.test <- sea.test[!duplicated(sea.test$lon) & !duplicated(sea.test$lat),]
  
  try(
  sea.test.clean <- cc_sea(x = sea.test,
                            lon = "lon",
                            lat = "lat")
  )
  
  if(dim(sea.test)[1] >= 25 & dim(sea.test)[1] < 15000 & exists("sea.test.clean")){
    sea.out.test.clean <- cc_outl(x = sea.test.clean,
                             lon = "lon",
                             lat = "lat",
                             species = "species",
                             method = "quantile")
    sea.out.test.clean <- sea.out.test.clean[!is.na(sea.out.test.clean$lat) &
                                               !is.na(sea.out.test.clean$lon),]
    print(paste0(dim(sea.test)[1] - dim(sea.out.test.clean)[1], " records actually flagged"))
  } # outlier if statement
  
  # manually defined outlier removal for species with over 15000 records
  if(sp == "Natrix helvetica" | sp == "Natrix maura"){
    # Natrix maura and Natrix helvetica
    x_coord <- c(-25, -25, 30, 30)
    y_coord <- c(25, 80, 80, 25)
    num.removed.seaout <- sea.test.clean[sea.test.clean$lon >= x_coord[1] &
                                           sea.test.clean$lon <= x_coord[3] &
                                           sea.test.clean$lat >= y_coord[1] &
                                           sea.test.clean$lat <= y_coord[3],]
  } else if(sp == "Natrix natrix") {
    # Natrix natrix 
    x_coord <- c(-25, -25, 100, 100)
    y_coord <- c(25, 80, 80, 25)
    num.removed.seaout <- sea.test.clean[sea.test.clean$lon >= x_coord[1] &
                                           sea.test.clean$lon <= x_coord[3] &
                                           sea.test.clean$lat >= y_coord[1] &
                                           sea.test.clean$lat <= y_coord[3],]
  } else if(sp == "Thamnophis elegans" | sp == "Thamnophis sirtalis") {
    # Thamnophis elegans and sirtalis
    x_coord <- c(-140, -140, -40, -40)
    y_coord <- c(10, 75, 75, 10)
    num.removed.seaout <- sea.test.clean[sea.test.clean$lon >= x_coord[1] &
                                           sea.test.clean$lon <= x_coord[3] &
                                           sea.test.clean$lat >= y_coord[1] &
                                           sea.test.clean$lat <= y_coord[3],]
  } else if(sp == "Vipera berus") {
    # Vipera berus
    x_coord <- c(-25, -25, 170, 170)
    y_coord <- c(25, 80, 80, 25)
    num.removed.seaout <- sea.test.clean[sea.test.clean$lon >= x_coord[1] &
                                           sea.test.clean$lon <= x_coord[3] &
                                           sea.test.clean$lat >= y_coord[1] &
                                           sea.test.clean$lat <= y_coord[3],]
  } # end of if statement to deal with large numbered species
  
  if(exists("sea.out.test.clean")){
    num.removed.seaout <- dim(sea.test)[1] - dim(sea.out.test.clean)[1]
    num.removed.sea <- dim(sea.test)[1] - dim(sea.test.clean)[1]
    removed <- data.frame(sp, num.removed.sea, num.removed.seaout)
    removed.list[[i]] <- removed
    sea.clean.list[[i]] <- sea.test.clean
  } else if(exists("sea.test.clean")){
    num.removed <- dim(sea.test)[1] - dim(sea.test.clean)[1]
    removed <- data.frame(sp, num.removed.sea, 0)
    names(removed) <- c("sp", "num.removed.sea", "num.removed.seaout")
    removed.list[[i]] <- removed
    sea.clean.list[[i]] <- sea.test.clean
  } else {
    i.i <- i.i + 1
    issue.list[[i.i]] <- sp
  } # end of sea clean if statement
  
} # species loop end
issue.list
save(issue.list, file = paste0(loc.gbifdata, "sea_issue_list.rData"))

removed.values <- rbindlist(removed.list, fill = TRUE)
write.csv(x = removed.values, file = paste0(loc.gbifdata, "removed_values.csv"), row.names = FALSE)

gbif.data.clean.sea <- rbindlist(sea.clean.list, fill = TRUE)

gbif.data.clean <- gbif.data.clean.sea

write.csv(file = paste0(loc.gbifdata, "6_compiled_GBIFData_coordclean2.csv"), row.names = FALSE,
          x = gbif.data.clean)


# Decimal degrees conversion check ----------------------------------------

gbif.data.clean <- read.csv(file = paste0(loc.gbifdata, "6_compiled_GBIFData_coordclean2.csv"),
                            stringsAsFactors = FALSE)

### need to break down the species, no point, do it by group
gbif.data.clean$dc_group <- rep(seq(1,20,1), length.out = dim(gbif.data.clean)[1])
out.ddmm <- dc_ddmm(gbif.data.clean, lon = "lon", lat = "lat", 
                    ds = "dc_group", diagnostic = TRUE, diff = 0.01,
                    value = "dataset")
out.ddmm
gbif.data.clean <- dplyr::select(gbif.data.clean, -dc_group)

write.csv(file = paste0(loc.gbifdata, "7_compiled_GBIFData_coordcleanComplete.csv"), row.names = FALSE,
          x = gbif.data.clean)

# Keyword generation for Flickr -------------------------------------------
# generating a list of keywords to feed flickr API

# create a list of species that should be present using GBIF and NCBI databases
i <- 0
species.list <- list()
for(fam in snake.families){
  i <- i + 1
  
  try(rm(gbif.ids))
  try(rm(gbif.species))
  
  try(gbif.ids <- get_gbifid(fam))
  try(gbif.species <- gbif_downstream(gbif.ids[1], downto = "Species"))
  
  ncbi.ids <- get_ids(fam, db = "ncbi")
  ncbi.species <- ncbi_downstream(ncbi.ids$ncbi[1], downto = "Species")
  
  if(exists("gbif.species")){
    if(dim(gbif.species)[1] > 1){
    gbif.res <- data.frame("name" = gbif.species$name,
                           "id" = gbif.species$key, "source" = "gbif")
    ncbi.res <- data.frame("name" = ncbi.species$childtaxa_name,
                           "id" = ncbi.species$childtaxa_id, "source" = "ncbi")
    sp <- rbind(gbif.res, ncbi.res)
    } else {
    sp <- data.frame("name" = ncbi.species$childtaxa_name,
                     "id" = ncbi.species$childtaxa_id, "source" = "ncbi")
    } # end of internal if when df has 0 rows
  } else {
    sp <- data.frame("name" = ncbi.species$childtaxa_name,
                     "id" = ncbi.species$childtaxa_id, "source" = "ncbi")
  } # end of exists if statement
  species.list[[i]] <- sp
}
species <- rbindlist(species.list)

write.csv(file = paste0(loc.gbifdata, "species_list.csv"), row.names = FALSE,
          x = species)

species <- read.csv(file = paste0(loc.gbifdata, "species_list.csv"),
                    stringsAsFactors = FALSE)

# run through the list and get all the common names from GBIF and NCBI, using
# the id for each database.
# UNUSED code due to high failure rate to retrieve common names
# i <- 0
# species.names <- list()
# for(sp in species$name){
#   (try(gbif.name <- sci2comm(sp, db = "gbif", simplify = TRUE)))
#   (try(ncbi.name <- sci2comm(sp, db = "ncbi", simplify = TRUE)))
#   
#   i <- i + 1
#   
#   if(exists("gbif.name") & exists("ncbi")){
#     sp.comname <- c(ncbi.name[[1]], ncbi.name[[1]])
#     sp.comname <- sp.comname[!duplicated(sp.comname)]
#     all.names <- c(sp, sp.comname)
#     species.names[[i]] <- all.names
#   } else if(exists("gbif.name") & !exists("ncbi")){
#     all.names <- gbif.name
#     species.names[[i]] <- all.names
#   } else if(!exists("gbif.name") & exists("ncbi")){
#     all.names <- ncbi.name
#     species.names[[i]] <- all.names
#   } else {next}
# } # 


# Loop to get common names from Reptile Database --------------------------

set.seed(1)
i <- 0
species.for.loop <- species$name[!duplicated(species$name)]
repdb.names <- vector(mode = "list", length = length(species.for.loop))
for(spec in species.for.loop){
  
  # silenced removal of object each iteration
  tryCatch({
    rm(comm.names)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  tryCatch({
    rm(h.text)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  
  i <- i + 1
  
  if(any(list.files() %in% "repdb_commnames.rdata")){
    repdb.names <- readRDS(file = "repdb_commnames.Rds")
  }
  
  # skip over names already pulled, allows resuming of loop after internet outage
  if( !is.null(repdb.names[[i]]) ){
    print(paste0(spec, " - skipped"))
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
    names.df <- data.frame("SName" = spec, "CName" = "Not found")
    repdb.names[[i]] <- names.df
    saveRDS(repdb.names, file = "repdb_commnames.Rds", compress = FALSE)
    {next}
  }
  
  html.raw <- htmlTreeParse(h.url, useInternal = TRUE)
  h.text <- xpathApply(html.raw, '//td',
                       function(x) xpathSApply(x,".//text()", xmlValue))[[6]]
  
  if(length(h.text) == 0){
    print(paste0(spec, " - None found"))
    names.df <- data.frame("SName" = spec, "CName" = NA)
    repdb.names[[i]] <- names.df
    saveRDS(repdb.names, file = "repdb_commnames.Rds", compress = FALSE)
    {next}
  } # end of no character if
  
  if(any(nchar(h.text) <= 2)){
    print(paste0(spec, " - None found"))
    names.df <- data.frame("SName" = spec, "CName" = NA)
    repdb.names[[i]] <- names.df
    saveRDS(repdb.names, file = "repdb_commnames.Rds", compress = FALSE)
    {next}
  } # end of junk character if 
  
  if(any(str_detect(h.text, "\\\\:"))){
    if(any(str_detect(h.text, "morph"))){
      comm.names <- sub("\\\\(color morphs\\\\:(.*?)", ",", h.text)
      comm.names <- sub("^.*\\\\:", "", comm.names)
    } else {
      comm.names <- sub("^.*\\\\:", "", h.text)
    } # splitting by :
    
    # remove trailing white space
    comm.names <- sub("^\\\\s+", "", comm.names)
    comm.names <- sub("\\\\s+$", "", comm.names)
    
    # for loop to remove bracketed things that are not at the start of the string
    for(cn in comm.names){
      if(str_detect(cn, "\\\\(") & !grepl("^\\\\(", x = cn)){
        cn.sub <- sub("\\\\(.*", "", cn)
        comm.names[comm.names == cn] <- cn.sub
        ## one off fix for dormideira (Dipsas turgidus)
        if(cn == "dormideira"){
          cn <- "Dormideira"
          comm.names[comm.names == cn] <- cn.sub
        } # end of dormideira if
      } # end of bracket if
    } # end of loop search for brackets
    
    # for those times that you get E: bottae: Genus species
    # comm.names <- grep('^[A-Z]|\\\\(', comm.names, value = TRUE)
    # starts with a captial, bracket or contains spaces
    comm.names <- comm.names[str_detect(comm.names, "^[A-Z]|^\\\\(|\\\\s")]
    
  } else {
    comm.names <- unlist(h.text)
  }
  
  # loop to split any combined common names up by seperators
  for(n.pat in comm.names){
    if(str_detect(n.pat, pattern = ",")){
      c.nam <- str_split(n.pat, pattern = ",\\\\ ")[[1]]
      comm.names <- comm.names[!comm.names %in% n.pat]
      comm.names <- c(comm.names, c.nam)
    } # end of if loop to seperate multiple same language names
    comm.names <- sub("^\\\\s+", "", comm.names)
    comm.names <- sub("\\\\s+$", "", comm.names)
  } # end of common name clean-up loop
  comm.names <- comm.names[!comm.names == ""]
  
  if(exists("comm.names")){
    print(paste0(spec, " - ", length(comm.names), " names found" ))
    names.df <- data.frame("SName" = spec, "CName" = comm.names)
    repdb.names[[i]] <- names.df
    saveRDS(repdb.names, file = "repdb_commnames.Rds", compress = FALSE)
  } # end of if to add common names to file
} # end of scientific name loop
saveRDS(repdb.names, file = "repdb_commnames.Rds", compress = FALSE)

repdb.names <- readRDS(file = "repdb_commnames.Rds")

#

# Pulling GBIF data citations ---------------------------------------------

gbif.data <- fread(file = paste0(loc.gbifdata, "1_compiled_GBIFData.csv"))

downloaddate <- gbif.data$downloadDate[!duplicated(gbif.data$gbifID)]
ids <- gbif.data$datasetKey[!duplicated(gbif.data$datasetKey)]

rm(gbif.data)

data.citations <- list()
i <- 0
for(i in 1:length(ids)){
  print(paste(i, "/", length(ids)))
  i <- i+1
  
  date <- downloaddate[i]
  id <- ids[i]
  cite <- gbif_citation(id)
  cite.string <- cite$citation$citation
  cite.corrected <- str_replace_all(cite.string,
                                    pattern = "2019-06-12", replacement = date)
  data.citations[[i]] <- cite.corrected
}

data.citations <- unlist(data.citations)
data.citations <- data.citations[!duplicated(data.citations)]
data.citations <- data.citations[!is.null(data.citations)]

write.csv(file = "GBIF data sources.csv", row.names = FALSE,
          x = data.citations)

