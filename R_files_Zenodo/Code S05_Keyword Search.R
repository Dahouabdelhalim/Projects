
library(dplyr)
library(stringr)
library(xml2)
library(rvest)

# Review page counts ------------------------------------------------------

folders <- list.dirs("./Data/WebData", recursive = FALSE)

pageCounts <- do.call(rbind, lapply(folders, function(fold){
  if(length(list.files(fold, full.names = TRUE, pattern = ".html")) == 0){
    html.files <- list.files(fold, full.names = TRUE, pattern = ".html", 
                             recursive = TRUE)
  } else {
    html.files <- list.files(fold, full.names = TRUE, pattern = ".html")
  }
  
  return(data.frame(webID = sub("./Data/WebData/", "", fold),
                    retPages = length(html.files),
                    retdate = file.info(html.files[1])$mtime))
}))

pageCounts

write.csv(x = pageCounts, file = "./Data/Page counts.csv",
          row.names = FALSE)

sum(pageCounts$retPages)

targetSites <- read.csv(file = "./Data/Target Websites.csv", stringsAsFactors = FALSE)

targetSites <- left_join(targetSites, pageCounts) %>% 
  mutate(retPages = ifelse(method == "pdf", 1, retPages))

sitesNoPages <- targetSites %>% 
  filter(retPages == 0 | is.na(retPages)) %>% 
  pull(webID)

targetSites %>% 
  filter(!webID %in% sitesNoPages) %>% 
  summarise(mean(retPages),
            se = sqrt(var(retPages)/length(retPages)),
            max(retPages),
            sum(retPages))

# Keyword list generation -------------------------------------------------

araData <- read.csv("./Data/arachnid_species_data_syns.csv")

# make sure hte non spiders have something for the keywords
araData$allNames[!araData$clade == "Spider"] <- araData$accName[!araData$clade == "Spider"]
araData$allGenera[!araData$clade == "Spider"] <- word(araData$accName[!araData$clade == "Spider"], 1, 1)

# araData$allNames[!araData$clade == "Spider"]
# araData$allNames[araData$genus == "Pandinus"]
# araData$accName[araData$genus == "Pandinus"]

speciesKeywords <- lapply(as.list(araData$allNames), function(x){
  vectorKW <- c(str_split(x, ";", simplify = TRUE))
  vectorKW <- vectorKW[!is.na(vectorKW) & !vectorKW == "NA" & !vectorKW == " "]
  return(vectorKW)
})
genusKeywords <- lapply(as.list(araData$allGenera), function(x){
  vectorKW <- c(str_split(x, ";", simplify = TRUE))
  vectorKW <- vectorKW[!is.na(vectorKW) & !vectorKW == "NA" & !vectorKW == " "]
  return(vectorKW)
})
genusKeywords <- unique(unlist(genusKeywords))
genusKeywords <- genusKeywords[!genusKeywords == "" & !is.na(genusKeywords) & !genusKeywords == " "]
# genusKeywords <- paste0(" ", genusKeywords, " ")
genusKeywords <- as.list(genusKeywords)

speciesKeywords[which(sapply(speciesKeywords, function(x){
  any(str_detect(x, fixed("Clubiona", ignore_case = TRUE)))
}))]

# speciesKeywords[which(unlist(lapply(speciesKeywords,
#                                     function(x){
#                                       any(str_detect(x, fixed("Clubiona hummeli", ignore_case = TRUE)))
#                                     }
# )))]

nTerms <- sapply(speciesKeywords, function(x){
  length(unique(x))
})
mean(nTerms)
sqrt(var(nTerms)/length(nTerms))
length(unlist(speciesKeywords))
range(nTerms)

# Keyword search functions ------------------------------------------------

# html.locs <- html.files
trade_prep_text <- function(htmlLocs, fold = fold){
  simpHTMLList <- vector(mode = "list", length = length(htmlLocs))
  i <- 1
  for(h in htmlLocs){
    # h <- html.locs[1]
    i <- i+1
    rawHTML <- paste(readLines(h),
                     collapse="\\n")
    if(rawHTML == ""){ #  for the rare occassion teh webpage is completely empty
      {next}
    }
    rawHTML <- rawHTML %>%
      xml2::read_html(options = "HUGE") %>% # had to set the HUGE option because of sites with long xml lines
      rvest::html_text()
    # simplfy the html so we are only dealing with text
    simpHTML <- stringr::str_replace_all(rawHTML, "[^[:alnum:]]", " ")
    simpHTML <- stringr::str_replace_all(simpHTML, pattern = "[[:digit:]]", " ")
    simpHTML <- stringr::str_replace(gsub("\\\\s+", " ", str_trim(simpHTML)), "B", "b")
    simpHTML <- stringr::str_squish(simpHTML)
    # simp.HTML <- paste(paste0(str_extract_all(sub(fold, "", h),
    #                                           "[[:digit:]]", simplify = TRUE),
    #                           collapse = ""), "     ", simp.HTML) # adds the page id in an easily accessible place
    simpHTML <- paste(paste0(stringr::str_extract_all(sub(gsub("(.*)/.*", "\\\\1", h), "", h),
                                                      "[[:digit:]]", simplify = TRUE),
                             collapse = ""), "     ", simpHTML) # adds the page id in an easily accessible place
    simpHTMLList[[i]] <- simpHTML
  }
  return(simpHTMLList)
}

basic_trade_search_html <- function(html, keywords){
  # html <- htmlList[[1]]
  # keywords <- speciesKeywords
  page <- str_trim(str_extract(html, "^.....")) # pulls the page ID and ditches empty space
  sppFound <- vector(mode = "list", length = length(keywords))
  i <- 0
  for(spkw in keywords){
    # spkw <- keywords[[1]]
    # spkw <- keywords[[8886]]
    # print(n)
    for(n in 1:length(spkw)){
      # n <- 3
      # n <- 1
      kw <- spkw[n]
      
      y <- stringr::str_extract_all(html, fixed(kw, ignore_case = TRUE),
                                    simplify = TRUE)
      if(assertthat::not_empty(y)){
        # print(assertthat::not_empty(y))
        i <- i+1
        sppFound[[i]] <- data.frame(sp = spkw[1],
                                    page = page,
                                    keyw = kw)
      }
    }
  }
  message(paste0(page, " --- page complete"))
  # df.1 <- as.data.frame(table(unlist(spp.found[-which(sapply(spp.found, is.null))])))
  # return(df.1)
  sppDf <- do.call(rbind, sppFound[-which(sapply(sppFound, is.null))])
  if(!is.null(sppDf)){
    sppDf <- sppDf %>% 
      dplyr::group_by(sp) %>% 
      dplyr::add_count() %>% 
      dplyr::slice(n = 1) %>% 
      dplyr::filter(!is.na(sp))
    return(sppDf)
  }
  # else {
  #   spp.df <- data.frame(sp = NA, page = page)
  #   return(spp.df)
  # }
}

genusTiered_trade_search_html <- function(html, keywords, genkeywords){
  # html <- htmlList[[1]]
  # keywords <- speciesKeywords
  # genkeywords <- genusKeywords
  page <- str_trim(str_extract(html, "^.....")) # pulls the page ID and ditches empty space
  genFound <- vector(mode = "list", length = length(genkeywords))
  
  i <- 0
  for(genkw in genkeywords){
    # genkw <- genkeywords[[1]]
    # genkw <- "Cyriopagopus"
    # genkw <- "Brachypelma"
    y <- stringr::str_extract_all(html, regex(paste0("\\\\s", genkw),
                                              ignore_case = TRUE), simplify = TRUE)
    if(assertthat::not_empty(y)){
      # print(assertthat::not_empty(y))
      # end the 2 words in front and 3 behind of genus hit
      termSurround <- stringr::str_extract_all(html, regex(paste0("([^\\\\s]+\\\\s){3}", genkw, "(\\\\s[^\\\\s]+){4}"),
                                                    ignore_case = TRUE), simplify = TRUE)
      termSurround <- paste(termSurround, collapse = "; ")
      i <- i+1
      genFound[[i]] <- data.frame(sp = genkw,
                                  page = page,
                                  keyw = genkw,
                                  spORgen = "GENUS",
                                  termsSurrounding = termSurround)
    } # end of if if something is detected
  } # end of kw cycle
  genDf <- do.call(rbind, genFound[-which(sapply(genFound, is.null))])
  
  if(!is.null(genDf)){
    # trim the species list to only look for those where the genus has been detected
    keywordsTrim <- speciesKeywords[which(sapply(speciesKeywords, function(x){
      any(str_detect(x, paste(genDf$sp, collapse = "|")))
    }))]
    
    sppFound <- vector(mode = "list", length = length(keywordsTrim))
    i <- 0
    for(spkw in keywordsTrim){
      # spkw <- keywords[[1]]
      # spkw <- keywords[[8886]]
      # print(n)
      for(n in 1:length(spkw)){
        # n <- 3
        # n <- 1
        kw <- spkw[n]
        
        y <- stringr::str_extract_all(html, fixed(kw, ignore_case = TRUE),
                                      simplify = TRUE)
        if(assertthat::not_empty(y)){
          # print(assertthat::not_empty(y))
          i <- i+1
          sppFound[[i]] <- data.frame(sp = spkw[1],
                                      page = page,
                                      keyw = kw,
                                      spORgen = "SPECIES",
                                      termsSurrounding = NA)
        }
      }
    }
    message(paste0(page, " --- page complete"))
    
    sppDf <- do.call(rbind, sppFound[-which(sapply(sppFound, is.null))])
    sppDf <- rbind(genDf, sppDf)
    if(!is.null(sppDf)){
      return(sppDf)
    } # if to retrun a result
  } # if to search of spp only if genus detected
}

# Run the search on comtemp sites --------------------------------------------

folders <- list.dirs("./Data/WebData", recursive = FALSE)

targetSites <- read.csv(file = "./Data/Target Websites.csv", stringsAsFactors = FALSE)

pdfSites <- targetSites %>% 
  filter(method == "pdf") %>%
  pull(webID)

start <- Sys.time()
# listSites <- vector(mode = "list", length = length(folders))
j <- 0
for(fold in folders[!folders %in% paste0("./Data/WebData/", sitesNoPages)]){
  # fold <- folders[2]
  j <- j + 1
  print(fold)
  
  if(sub("./Data/WebData/", "", fold) %in% pdfSites){
    print(paste0(fold, " - skipped because PDF file"))
    {next}
  }
  
  if("KEYWORD_EXTRACT_SPECIES.csv" %in% list.files(fold)){
    print(paste0(fold, " - skipped"))
    {next}
  }
  
  if(length(list.files(fold, full.names = TRUE, pattern = ".html")) == 0){
    htmlFiles <- list.files(fold, full.names = TRUE, pattern = ".html", 
                            recursive = TRUE)
  } else {
    htmlFiles <- list.files(fold, full.names = TRUE, pattern = ".html")
  }
  # fold <- folders[3]
  print(paste0("--- Loading ", length(htmlFiles)," html pages"))
  # read and clean html files
  htmlList <- trade_prep_text(htmlFiles, fold = fold)
  # remove NULLs and duplicates
  htmlList <- htmlList[-which(sapply(htmlList, is.null))]
  htmlList <- htmlList[!duplicated(htmlList)]
  
  print("--- Searching for keywords...")
  
  counts <- lapply(htmlList, function(x){
    genusTiered_trade_search_html(html = x, keywords = speciesKeywords, genkeywords = genusKeywords)
  })
  siteResults <- do.call(rbind, counts)
  siteResults$webID <- sub("./Data/WebData/", "", fold)
  write.csv(x = siteResults, 
            file = paste0(fold, "/KEYWORD_EXTRACT_SPECIES.csv"),
            row.names = FALSE)
  
}
end <- Sys.time()

# genusTiered_trade_search_html(htmlList[[1]], speciesKeywords, genusKeywords)

# Run the search on PDF files ---------------------------------------------

library(pdftools)

pdfSites <- targetSites %>% 
  filter(method == "pdf") %>% 
  select(webID)

start <- Sys.time()
for(site in pdfSites){
  
  pdfText <- pdf_text(paste0("./Data/WebData/", site, "/1.pdf"))
  
  print(paste0(site, " --- keyword searching"))
  
  sitesResults <- genusTiered_trade_search_html(html = pdfText, keywords = speciesKeywords, genkeywords = genusKeywords)
  sitesResults$page <- 1
  sitesResults$webID <- site
  write.csv(x = sitesResults, 
            file = paste0("./Data/WebData/", site, "/KEYWORD_EXTRACT_SPECIES.csv"),
            row.names = FALSE)
  
}# for loop end
end <- Sys.time()

end-start

# Run the search on wayback pages -----------------------------------------

tempHTML <- list.files("./Data/TemporalData/terraristik", pattern = "html", 
                       full.names = TRUE)

# create batches to run the search using, will help resume searches if interrupted
htmlChunks <- sort(rep_len(1:100, length.out = length(tempHTML)))
sum(htmlChunks == 1)

i <- 0
for(pagechunk in unique(htmlChunks)){
  # pagechunk <- unique(htmlChunks)[4]
  i <- i + 1
  
  chunk <- tempHTML[htmlChunks == pagechunk]
  
  from <- head(sub(".html", "", str_extract(chunk, "....\\\\.html")), 1)
  to <- tail(sub(".html", "", str_extract(chunk, "....\\\\.html")), 1)
  
  print(paste0("=== ", from, " - ", to, " ==="))
  
  outputName <- paste0("./Data/TemporalData/terraristik",
                       "/KEYWORD_EXTRACT_SPECIES_", from, "-", to, ".csv")
  
  if(outputName %in% list.files("./Data/TemporalData/terraristik", pattern = ".csv",
                                full.names = TRUE)){
    print("=== Already Complete ===")
    {next}
  }
  
  # read and clean html files
  htmlList <- trade_prep_text(chunk, fold = "terraristik")
  # remove NULLs and duplicates
  htmlList <- htmlList[-which(sapply(htmlList, is.null))]
  
  counts <- lapply(htmlList, function(x){
    genusTiered_trade_search_html(html = x,
                                  keywords = speciesKeywords, genkeywords = genusKeywords)
  })
  sitesResults <- do.call(rbind, counts)
  sitesResults$webID <- "terraristikTemporal"
  
  write.csv(x = sitesResults, 
            file = outputName,
            row.names = FALSE)
  
} # chunk loop end
