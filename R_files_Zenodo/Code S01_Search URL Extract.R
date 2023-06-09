

# Libraries ---------------------------------------------------------------

library(XML)
library(stringr)

# Functions ---------------------------------------------------------------

extract_baseURLS <- function(sFile = NA){
  
  if(is.null(sFile) | is.na(sFile)){
    cat("Search file location must be supplied")
  }
  
  if(stringr::str_detect(sFile, pattern = "JPN")){
    Sys.setlocale("LC_CTYPE", locale="Japanese")
  } else {
    Sys.setlocale("LC_CTYPE", locale="English_United Kingdom")
  }
  
  doc <- XML::htmlParse(sFile)
  links <- XML::xpathSApply(doc, "//a/@href")
  XML::free(doc)
  
  links <- links[!grepl(pattern = "^/", links)]
  links <- links[!stringr::str_detect(links, "\\\\.google") &
                   !stringr::str_detect(links, "\\\\.bing") &
                   !stringr::str_detect(links, "javascript") &
                   !stringr::str_detect(links, "#") &
                   !stringr::str_detect(links, "microsoft") &
                   !stringr::str_detect(links, "^/search\\\\?q=")]
  baseLinks <- sapply(links, function(x){
    paste(strsplit(x,"/")[[1]][1:3], collapse = "/")
  })
  
  fileSimp <- stringr::str_extract(sFile, "(?<=Results/).*(?=\\\\.htm)")
  
  searchRes <- data.frame(lang = stringr::str_extract(fileSimp, "ENG|FRA|GER|JPN|SPA|POR|CZE|POL|RUS"),
                          engine = stringr::str_extract(fileSimp, "Google|Bing"),
                          page = stringr::str_extract_all(fileSimp, "[:digit:]{2}", simplify = TRUE),
                          searchdate = file.info(sFile)$ctime,
                          link = baseLinks)
  searchRes <- searchRes[!duplicated(searchRes$link),]
}


# Run functions -----------------------------------------------------------

searchFiles <- list.files("./Data/Search Page Results", pattern = ".htm",
                          full.names = TRUE)

searchList <- lapply(searchFiles, function(x){
  extract_baseURLS(sFile = x)
})

searchDf <- do.call(rbind, searchList)

searchDf[,c("reviewdate", "sells", "allow", 
            "type", "order", "target", 
            "method", "refine", "spages")] <- NA

searchDf <- searchDf[!duplicated(searchDf$link),]

write.csv(x = searchDf, file = "./Data/Search Page Results/Search Results.csv",
          row.names = FALSE)

 