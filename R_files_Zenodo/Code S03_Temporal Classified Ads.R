# devtools::install_github("hrbrmstr/wayback")
library(wayback)
library(dplyr)
library(XML)
library(xml2)
library(rvest)
library(stringr)
library(ggplot2)
library(lubridate)
library(downloader)

archive_available("http://www.terraristik.com/tb/list_classifieds_int.php",
                  "20161010101010")

get_mementos("http://www.terraristik.com/tb/list_classifieds_int.php",
             "20161010101010")

b.url <- "http://www.terraristik.com/tb/list_classifieds_int.php"

bq <- cdx_basic_query(b.url, match_type = "prefix",
                      limit = 100000) %>% 
  print(n = 100)

mainurl <- "http://web.archive.org/cdx/search/cdx?url=http://www.terraristik.com/tb/list_classifieds_int.php&matchType=host&limit=1000000&output=json"

res <- httr::GET(mainurl)
res <- httr::content(res, as = "text", encoding = "UTF-8")
res <- jsonlite::fromJSON(res, simplifyVector = TRUE, simplifyDataFrame = TRUE)
res <- tibble::as_tibble(res)
res <- stats::setNames(res[-1, ], purrr::flatten_chr(res[1,]))

res <- res %>% 
  filter(str_detect(original, "cate"),
         statuscode == 200) %>% 
  mutate(ia.url  = paste0("https://web.archive.org/web/",
                          timestamp, "/",
                          original))

res <- dplyr::mutate_(res, timestamp.parse = lazyeval::interp(~anytime::anytime(t),
                                                              t = quote(timestamp)))

res %>% 
  mutate(year = year(timestamp.parse)) %>% 
  ggplot() +
  geom_bar(aes(x = year))

dir.create("")

write.csv(x = res, file = "./Data/TemporalData/wayback_terraristik_results.csv",
          row.names = FALSE)

# Downloading the wayback results -------------------------------------------------------------

res <- read.csv(file = "./Data/TemporalData/wayback_terraristik_results.csv",
                stringsAsFactors = FALSE)

singlepageurl <- res$ia.url

dir.create(path = paste0("./Data/TemporalData/", "terraristik"))

i <- 0
for(url in singlepageurl){
  i <- i+1
  # url <- singlepageurl[4]
  fold <- "terraristik"
  
  print(paste0("Start: ", url, " --- ", 
               match(url, singlepageurl), "/", length(singlepageurl)))
  
  number <- str_pad(i, width = 4, side = "left", pad = 0)
  
  if(
    paste0("./Data/TemporalData/", fold, "/", number, ".html") %in%
    list.files(path = paste0("./Data/TemporalData/", fold, "/"), full.names = TRUE)
  ){
    print(paste0(number, " --- skipped"))
    {next}
  }
  
  try(
    download(url = url, paste0("./Data/TemporalData/", fold, "/", number, ".html"))
  )
  
  Sys.sleep(5)
  
  # cat(file = paste0("./TemporalTrend/", fold, "_complete.txt"), url)
} # end of url loop
