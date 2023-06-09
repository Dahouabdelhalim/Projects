## 02 - Flickr search and data retrieval

# Libraries ---------------------------------------------------------------

library(RCurl)
library(XML)
library(httr)
library(readr)
library(lubridate)
library(stringr)
library(dplyr)
library(CoordinateCleaner)

# Flickr API key set-up ---------------------------------------------------

# API key and secret must be obtained from https://www.flickr.com/services/api/misc.api_keys.html
api_key <- "[FLICKR API KEY HERE]"

secret <- "[FLICKR SECRET KEY HERE]"

# creates the app passing the key and secret
myapp <- oauth_app("[FLICKR PROJECT NAME HERE]", key = api_key, secret = secret)

# get authentication credentials from the API
ep <- oauth_endpoint(request = "https://www.flickr.com/services/oauth/request_token",
                     authorize = "https://www.flickr.com/services/oauth/authorize",
                     access = "https://www.flickr.com/services/oauth/access_token")
### OPENS BROWSER ###

# creates variable with authentication credentials
sig <- oauth1.0_token(ep, myapp, cache = FALSE)

# authenticate
fl_sig <- sign_oauth1.0(myapp, sig)

URLbase <- "https://api.flickr.com/services/rest/?method=flickr.photos.search"
api <- api_key


# read in the previously generated list from script 1
repdb.names <- readRDS(file = "repdb_commnames.Rds")

i <- 0
kw.i <- 0
kw.zero <- list()

flickr.results <- vector(mode = "list", length = length(repdb.names))
for(item in repdb.names){
  
  tryCatch({
    rm(pics.species)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  tryCatch({
    rm(CName)
  }, warning = function(w) {
    # message("comm.names already removed")
  })
  
  i <- i + 1
  
  SName <- as.character(item$SName[1])
  CName <- as.character(item$CName)
  
  print(paste0("--- ", SName, ": Start ---"))
  
  if( (CName == "Not found" || CName == "None found" || is.na(CName)) ){
    v.keywords <- SName
  } else {
    v.keywords <- c(SName, CName)
  } # end of if for not found if statement
  
  v.keywords <- gsub(" ", "+", v.keywords)
  
  if(any(list.files() %in% "flickr_results.Rds")){
    flickr.results <- readRDS(file = "flickr_results.Rds")
  }
  
  if( !is.null(flickr.results[[i]]) ){
    print(paste0(SName, " - skipped"))
    {next}
  }
  
  ## information to generate the search results
  # tag mode, all with only get photos will all keywords
  tag_m <- "all"
  tags <- "snake"
  min.date <- "2005-01-01"
  max.date <- "2019-05-01"
  # geo-tagged only
  geo.call <- "1"
  # get data and locations using extra.call
  extra.call <- "date_taken,geo"
  page.call <- "250"
  format.call <- "rest"

  # loop that cycles through all keywords for a species
  pics.species <- NULL
  for(kw in v.keywords){
    search.text <- kw
    pics <- NULL
    # form the search URL
    URL.call <- paste0(URLbase, "&api_key=", api, "&tags=", tags, "&tag_mode=", tag_m,
                       "&text=", search.text, "&bbox=",
                       "&min_taken_date=", min.date, "&max_taken_date=", max.date,
                       "&has_geo", geo.call, "&extras=", extra.call, "&per_page=", page.call,
                       "&page=1", "&format=", format.call
    )
    
    getPhotos_data <- xmlRoot(xmlTreeParse(getURL(URL.call,
                                                  ssl.verifypeer = FALSE,
                                                  useragent = "flickr"),
                                           useInternalNodes = TRUE ))
    
    num.res <- as.numeric(xmlAttrs(getPhotos_data[["photos"]])[4])
    
    # first call gives us information on how many pages we need to cycle through
    # or if there are no results move on to next keyword
    if(num.res == 0){
      kw.i <- kw.i + 1
      kw.zero[[kw.i]] <- kw
      print(paste("----- Keyword:", kw, "-", num.res,
                  "results -----"))
      {next}
    }
    
    print(paste("----- Keyword:", kw, "-", num.res,
                "results -----"))
    
    max.page.call <- as.numeric(xmlAttrs(getPhotos_data[["photos"]])[2])
    
    # because of limits of Flickr's results we need to cycle through all pages
    # one by one
    for(page in 1:max.page.call){
      # page <- 2
      
      print(paste("----- Page", page, "/", max.page.call, "-----"))
      
      URL.call <- paste0(URLbase, "&api_key=", api, "&text=", search.text,
                         "&min_taken_date=", min.date, "&max_taken_date=", max.date,
                         "&has_geo", geo.call, "&extras=", extra.call, "&per_page=", page.call,
                         "&page=", page, "&format=", format.call
      )
      
      getPhotos_data <- xmlRoot(xmlTreeParse(getURL(URL.call,
                                                    ssl.verifypeer = FALSE,
                                                    useragent = "flickr"),
                                             useInternalNodes = TRUE ))
      
      # extract photo id
      id <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "id")
      # extract photo farm
      farm <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "farm")
      # extract photo secret
      secret <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "secret")
      # extract photo server
      server <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "server")
      # extract photo title
      title <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "title")
      # extract user id
      owner <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "owner")
      # extract date picture was taken
      datetaken <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr, "datetaken")
      # extract latitude
      latitude <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr,"latitude")
      # extract longitude
      longitude <- xpathSApply(getPhotos_data, "//photo", xmlGetAttr,"longitude")
      
      # put together in a df
      tmp_df <- data.frame(cbind(id, title, farm, server, secret, owner,
                                 datetaken, latitude, longitude),
                           stringsAsFactors = FALSE)  
      
      tmp_df$page <- page
      tmp_df$keyword <- kw
      pics <- rbind(pics, tmp_df)
      
    } # end of page loop
    pics.species <- rbind(pics.species, pics)
  } # end of keyword loop
  
  # required to fully skipp zero result species
  if(is.null(pics.species)){
    flickr.results[[i]] <- data.frame("id" = NA, "title" = NA, "farm" = NA, "server" = NA, "secret" = NA, "owner" = NA, "datetaken" = NA, "latitude" = NA, 
               "longitude" = NA, "page" = NA, "keyword" = SName)
    saveRDS(flickr.results, file = "flickr_results.Rds", compress = FALSE)
    {next}
  }
  
  # remove duplicated photos appearing in multiple keyword searches
  pics.species <- pics.species[!duplicated(pics.species$id),]

  pics.species$latitude <- as.numeric(pics.species$latitude)
  pics.species$longitude <- as.numeric(pics.species$longitude)
  
  # we also need to generate a URL using the photo and user ID so we can verify
  # the ID for chosen species
  pics.species$url <- NA
  for(rownum in 1:dim(pics.species)[1]){
    pics.species[rownum,]$url <- paste0("https://www.flickr.com/photos/",
                                        pics.species$owner[rownum],
                                        "/", pics.species$id[rownum], "/")
  }
  
  # delete data with latitude recorded as 0
  pics.species <- pics.species[!pics.species$latitude==0 &
                                 !pics.species$longitude == 0,]
  
  flickr.results[[i]] <- pics.species
  saveRDS(flickr.results, file = "flickr_results.Rds", compress = FALSE)
  print(paste0("--- ", SName, ": End and Saved ---"))
} #  end of items in repdb loop
saveRDS(flickr.results, file = "flickr_results.Rds", compress = FALSE)

flickr.results <- readRDS(file = "flickr_results.Rds")
f.res <- dplyr::bind_rows(flickr.results)
write.csv(file = "./Data/Flickr_results.csv", row.names = FALSE,
          x = f.res)


# Clean Flickr results ----------------------------------------------------

flr.data <- read.csv(file = "./Data/Flickr_results.csv",
                      stringsAsFactors = FALSE)
# remove duplicated photos
flr.data <- flr.data[!duplicated(flr.data$id),]
names(flr.data)
# one crazy point in the deep deep south - error
flr.data <- flr.data[!flr.data$latitude == min(flr.data$latitude),]

flr.data <- flr.data[!is.na(flr.data$latitude) | !is.na(flr.data$longitude),]
flr.data <- cc_sea(x = flr.data,
                   lon = "longitude",
                   lat = "latitude")

write.csv(file = "./Data/Flickr_results_clean.csv", row.names = FALSE,
          x = flr.data)

flr.data.raw <- read.csv(file = "./Data/Flickr_results_clean.csv",
                         stringsAsFactors = FALSE)
