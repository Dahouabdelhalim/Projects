
# allFolders <- list.files("./Data/WebData")
# allTxt <- allFolders[str_detect(allFolders, ".txt", negate = FALSE)]
# allFolders <- allFolders[str_detect(allFolders, ".txt", negate = TRUE)]
# 
# length(allTxt)
# length(allFolders)
# 
# sapply(allFolders, function(x){
#   length(list.files(paste0("./Data/WebData/", x), pattern = "html"))
# })
# 

# -------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(downloader)
library(Rcrawler)

searchResults <- read.csv(file = "./Data/Search Page Results/Search Results Reviewed.csv",
                          stringsAsFactors = FALSE)

# total number of sites reviewed from search results
dim(searchResults)[1]

sum(searchResults$sells == 0)
sum(searchResults$allow == 0, na.rm = TRUE)
sum(searchResults$sells == "dup")
sum(searchResults$sells == "access issue")
dim(searchResults)[1] - 906 -13 -6 -3

targetSites <- searchResults %>%
  filter(!is.na(method))

# number of suitable sites to search
dim(targetSites)[1]

# place to store the web pages
dir.create("./Data/WebData/")

targetSites$webID <- paste0("web", str_pad(1:length(targetSites$link), 3, "left", "0"))

write.csv(targetSites, file = "./Data/Target Websites.csv",
          row.names = FALSE)

targetSitesCensored <- targetSites
targetSitesCensored$link <- "REDACTED"
targetSitesCensored$target <- ifelse(!is.na(targetSitesCensored$target),
                                     "REDACTED",
                                     targetSitesCensored$target)
targetSitesCensored$refine <- ifelse(!is.na(targetSitesCensored$refine),
                                     "REDACTED",
                                     targetSitesCensored$refine)

write.csv(targetSitesCensored, file = "./Data/Target Websites Censored.csv",
          row.names = FALSE)

targetSites <- read.csv(file = "./Data/Target Websites.csv", stringsAsFactors = FALSE)

n_distinct(targetSites$webID)

# Single ------------------------------------------------------------------

singleTargets <- targetSites %>% 
  filter(method == "single") %>% 
  filter(!webID == "web073") # unsupported URL scheme

for(url in singleTargets$target){
  # url <- singleTargets$target[47]
  fold <- singleTargets$webID[singleTargets$target == url]
  dir.create(paste0("./Data/WebData/", fold), showWarnings = FALSE)
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in% 
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", url, " --- ", 
                   match(url, singleTargets$target), "/", length(singleTargets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", url, " --- ", 
               match(url, singleTargets$target), "/", length(singleTargets$target)))
  
  dir.create(path = paste0("./Data/WebData/", fold))
  # check for multiple pages
  if(str_detect(url, ";")){
    url <- str_split(url, ";")[[1]]
  }
  # cycle through the different single targets or get the single page
  for(turl in url){
    if(length(list.files(paste0("./Data/WebData/", fold), pattern = "html")) == 0){
      download(url = turl, paste0("./Data/WebData/", fold, "/1.html"))
    } else {
      num <- length(list.files(paste0("./Data/WebData/", fold), pattern = "html")) + 1
      download(url = turl, paste0("./Data/WebData/", fold, "/", num, ".html"))
    }
  }
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), url)
} # end of url loop


# Cycle -------------------------------------------------------------------

cycleTargets <- targetSites %>% 
  filter(method == "cycle" | method == "multiplecycle")

for(url in cycleTargets$target){
  
  # url <- cycleTargets$target[4]
  fold <- cycleTargets$webID[cycleTargets$target == url]
  spages <- cycleTargets$spages[cycleTargets$target == url]
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in% 
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", url, " --- ", 
                   match(url, cycleTargets$target), "/", length(cycleTargets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", url, " --- ", 
               match(url, cycleTargets$target), "/", length(cycleTargets$target)))
  
  dir.create(path = paste0("./Data/WebData/", fold))
  
  if(str_detect(url, ";")){
    url <- str_split(url, ";")[[1]]
    spages <- str_split(spages, ";")[[1]]
  }
  
  for(nurl in 1:length(url)){
    
    if(fold %in% c("web008", "web037")){ # for the sites that have ads not page numbers, increased by 20 each time
      urlPages <- str_replace(rep(url[nurl], as.numeric(spages[nurl])/20+1), "##",
                              as.character(seq(0, spages[nurl], 20)))
    } else {
      urlPages <- str_replace(rep(url[nurl], as.numeric(spages[nurl])+1), "##",
                              as.character(seq(0, spages[nurl], 1)))
    }
    
    for(page in 1:length(urlPages)){
      # page <- 1
      urlP <- urlPages[page]
      res <- try(
        download(url = urlP, paste0("./Data/WebData/", fold, "/",
                                    length(list.files(path = paste0("./Data/WebData/", fold),
                                                      pattern = "html"))+1, ".html")))
      # if(!res == 0){
      #   {break}
      # }
      Sys.sleep(10)
    }
  }
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), url)
} # end of url loop


# Level 1 -----------------------------------------------------------------

level1Targets <- targetSites %>% 
  filter(method == "level1")

# remove problem URLs
# level1Targets <- level1Targets[!level1Targets$webID %in% c("web010"),]

for(base in level1Targets$target){
  # base <- level1Targets$target[1]
  
  fold <- level1Targets$webID[level1Targets$target == base]
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in%
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", base, " --- ", 
                   match(base, level1Targets$target), "/", length(level1Targets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", base, " --- ", 
               match(base, level1Targets$target), "/", length(level1Targets$target)))
  
  if(!is.na(level1Targets$refine[level1Targets$target == base])){
    crawl.res <- Rcrawler(Website = base,
                          DIR = paste0("./Data/WebData/", fold),
                          MaxDepth = 1, RequestsDelay = 20, Obeyrobots = FALSE,
                          crawlUrlfilter = level1Targets$refine[level1Targets$target == base])
  } else {
    crawl.res <- Rcrawler(Website = base,
                          DIR = paste0("./Data/WebData/", fold),
                          MaxDepth = 1, RequestsDelay = 20, Obeyrobots = FALSE)
  } # end of keyword refinement
  
  write.csv(x = INDEX, file = paste0("./Data/WebData/", fold, "_INDEX.csv"), row.names = FALSE)
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), base)
} # end of url loop


# Cycle Level 1 -----------------------------------------------------------

level1CycleTargets <- targetSites %>% 
  filter(method == "cycle,level1")

# web061
# page 147 complete, see if it can be resumed

for(url in level1CycleTargets$target){
  # url <- level1CycleTargets$target[1]
  fold <- level1CycleTargets$webID[level1CycleTargets$target == url]
  spages <- level1CycleTargets$spages[level1CycleTargets$target == url]
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in% 
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", url, " --- ", 
                   match(url, level1CycleTargets$target), "/", length(level1CycleTargets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", url, " --- ", 
               match(url, level1CycleTargets$target), "/", length(level1CycleTargets$target)))
  
  dir.create(path = paste0("./Data/WebData/", fold))
  
  surlPages <- str_replace(rep(url, as.numeric(spages)+1), "##",
                           as.character(seq(0, spages, 1)))
  
  # for(page in 1:length(surlPages)){
  for(page in 148:length(surlPages)){
    # page <- 1
    urlP <- surlPages[page]
    
    res <- try(
      download(url = urlP, paste0("./Data/WebData/", fold, "/",
                                  length(list.files(path = paste0("./Data/WebData/", fold),
                                                    pattern = "html"))+1, ".html"))
    )
    
    if(!res == 0){
      {break}
    } else {
      
      if(!is.na(level1CycleTargets$refine[level1CycleTargets$target == url]) &
         !level1CycleTargets$refine[level1CycleTargets$target == url] == ""){
        crawl.res <- Rcrawler(Website = urlP,
                              DIR = paste0("./Data/WebData/", fold),
                              MaxDepth = 1, RequestsDelay = 20, Obeyrobots = FALSE,
                              crawlUrlfilter = level1CycleTargets$refine[level1CycleTargets$target == url])
      } else {
        crawl.res <- Rcrawler(Website = urlP,
                              DIR = paste0("./Data/WebData/", fold),
                              MaxDepth = 1, RequestsDelay = 20, Obeyrobots = FALSE)
      } # end of keyword refinement
      
      write.csv(x = INDEX,
                file = paste0("./Data/WebData/", fold, "/", page, "_", fold, "_INDEX.csv"),
                row.names = FALSE)
      
    } # end of res check
    
    Sys.sleep(10)
  }
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), url)
  
} # end of url loop

folders <- level1CycleTargets$webID

# Level 2 -----------------------------------------------------------------

level2Targets <- targetSites %>% 
  filter(method == "level2")

for(base in level2Targets$target){
  # base <- level1Targets$target[1]
  
  fold <- level2Targets$webID[level2Targets$target == base]
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in%
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", base, " --- ", 
                   match(base, level2Targets$target), "/", length(level2Targets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", base, " --- ", 
               match(base, level2Targets$target), "/", length(level2Targets$target)))
  
  if(!is.na(level2Targets$refine[level2Targets$target == base])){
    crawl.res <- Rcrawler(Website = base,
                          DIR = paste0("./Data/WebData/", fold),
                          MaxDepth = 2, RequestsDelay = 20, Obeyrobots = FALSE,
                          crawlUrlfilter = level2Targets$refine[level2Targets$target == base])
  } else {
    crawl.res <- Rcrawler(Website = base,
                          DIR = paste0("./Data/WebData/", fold),
                          MaxDepth = 2, RequestsDelay = 20, Obeyrobots = FALSE)
  } # end of keyword refinement
  
  write.csv(x = INDEX, file = paste0("./Data/WebData/", fold, "_INDEX.csv"), row.names = FALSE)
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), base)
} # end of url loop


# PDF ---------------------------------------------------------------------

pdfTargets <- targetSites %>% 
  filter(method == "pdf")

for(url in pdfTargets$target){
  # url <- singleTargets$target[47]
  fold <- pdfTargets$webID[pdfTargets$target == url]
  dir.create(paste0("./Data/WebData/", fold), showWarnings = FALSE)
  
  if(!length(gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))) == 0){
    if(fold %in% 
       gsub("_complete.txt", "", list.files(path = "./Data/WebData", pattern = ".txt"))){
      print(paste0("Skipped: ", url, " --- ", 
                   match(url, pdfTargets$target), "/", length(pdfTargets$target)))
      {next}
    } # check the complete files
  } # zero complete files if
  
  print(paste0("Start: ", url, " --- ", 
               match(url, pdfTargets$target), "/", length(pdfTargets$target)))
  
  dir.create(path = paste0("./Data/WebData/", fold))
  # check for multiple pages
  if(str_detect(url, ";")){
    url <- str_split(url, ";")[[1]]
  }
  # cycle through the different single targets or get the single page
  for(turl in url){
    if(length(list.files(paste0("./Data/WebData/", fold), pattern = "pdf")) == 0){
      download(url = turl, paste0("./Data/WebData/", fold, "/1.pdf"),
               mode="wb")
    } else {
      num <- length(list.files(paste0("./Data/WebData/", fold), pattern = "pdf")) + 1
      download(url = turl, paste0("./Data/WebData/", fold, "/", num, ".pdf"),
               mode="wb")
    }
  }
  
  cat(file = paste0("./Data/WebData/", fold, "_complete.txt"), url)
} # end of url loop

