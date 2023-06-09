
araData_trade <- read.csv("./Data/arachnid_species_data_traded.csv")

# Temporal online data ----------------------------------------------------

library(lubridate)
library(dplyr)
library(stringr)
library(ggpubr)
library(scico)

colourway <- scico(n = 9, palette = "imola")

extractFilesTemp <- list.files("./Data/TemporalData", pattern = "KEYWORD_EXTRACT",
                               recursive = TRUE, full.names = TRUE)

extractDataTemp <- do.call(rbind, lapply(extractFilesTemp, function(x){
  df <- read.csv(x, stringsAsFactors = FALSE)
  if(dim(df)[2] > 1){
    return(df)
  }
}))

wayData <- read.csv(file = "./Data/TemporalData/wayback_terraristik_results.csv",
                    stringsAsFactors = FALSE)

pageDates <- wayData %>% 
  mutate(page = row_number()) %>% 
  select(timestamp.parse, page)

temporalDataOnline <- inner_join(extractDataTemp, pageDates)

temporalDataOnline$timestamp.parse <- as.POSIXct(temporalDataOnline$timestamp.parse)
temporalDataOnline$year <- year(temporalDataOnline$timestamp.parse)

write.csv(x = temporalDataOnline,
          file = "./Data/Temporal Online Data.csv", row.names = FALSE)

onlinePlotData <- temporalDataOnline %>% 
  filter(spORgen == "SPECIES", str_detect(str_trim(keyw), " "),
         !str_detect(keyw, "nomen dubium")) %>%
  group_by(year) %>% 
  summarise(nSpecies = length(unique(sp))) %>% 
  mutate(source = "Online")

# LEMIS data convert names ------------------------------------------------

lemisData <- read.csv(file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv")

lemisData <- lemisData %>% 
  mutate(lemisName = str_to_sentence(paste(genus, species))) %>% 
  filter(!str_detect(lemisName, "sp\\\\.$") & !lemisName == "NA NA" & !lemisName == "Na na")

lemisPlotData <- lemisData %>% 
  select(lemisName, "year" = shipment_year) %>% 
  group_by(year) %>% 
  summarise(nSpecies = length(unique(lemisName))) %>% 
  mutate(source = "LEMIS")

# CITES trade database name convert ---------------------------------------------

citesData <- read.csv(file = "./Data/TradeDatabases/CITES_gross_imports_2021-09-15_comma_separated.csv")

citesData <- citesData %>% 
  filter(!str_detect(Taxon, "spp\\\\.$"))

names(citesData)

citesData <- citesData %>% 
  select(-App., -Term, -Unit, -Country)

for(i in 1:nrow(citesData)){
  citesData[i,!is.na(citesData[i,])] <- citesData[i,"Taxon"]
}

citesPlotData <- as.data.frame(apply(citesData, 2, function(x){
  length(unique(x[!is.na(x)]))
}))

citesPlotData$year <- as.numeric(sub("X", "", row.names(citesPlotData)))
names(citesPlotData) <- c("nSpecies", "year")
citesPlotData <- citesPlotData[!is.na(citesPlotData$year),]
citesPlotData$source <- "CITES"

# Combine data and plot ---------------------------------------------------

library(ggplot2)

tempPlotData <- rbind(onlinePlotData, lemisPlotData, citesPlotData)

# Raw counts plot ---------------------------------------------------------

pageData <- temporalDataOnline %>% 
  group_by(page, year) %>% 
  slice(n = 1) %>% 
  ungroup() %>% 
  count(year)


(rawCountsPlot <-
    tempPlotData %>% 
    filter(year < 2020) %>% 
    ggplot() +
    geom_line(aes(x = year, y = nSpecies, colour = source),
              size = 1.2) +
    geom_text(data = tempPlotData %>%
                group_by(source) %>% 
                filter(year == 2000),
              aes(x = year, y = nSpecies, colour = source, label = source),
              fontface = 2, nudge_y = -5, nudge_x = 0.1, hjust = 0, vjust = 1) +
    geom_text(data = tempPlotData %>%
                group_by(source) %>% 
                filter((year == 2002 & source == "Online")),
              aes(x = year, y = nSpecies, colour = source, label = source),
              fontface = 2, nudge_y = 35, nudge_x = -0.8, hjust = 0.5, vjust = 1) +
    geom_segment(data = tempPlotData %>%
                   group_by(source) %>% 
                   filter(source == "Online"),
                 aes(x = year, xend = year, y = nSpecies, yend = 700),
                 colour = colourway[2], alpha = 0.5, linetype = 2) +
    geom_point(data = pageData,
               aes(x = year, y = 700), size = 9,
               colour = colourway[2], alpha = 0.5) +
    geom_point(data = pageData,
               aes(x = year, y = 700), size = 7.5,
               colour = "white") +
    geom_text(data = pageData,
              aes(x = year, y = 700, label = n),
              size = 3, colour = colourway[1], alpha = 0.5, fontface = 2) +
    annotate("text", x = 2001.4, y = 700, label = "Pages\\nsearched",
             fontface = 4, hjust = 1, lineheight = 0.85, size = 3,
             colour = colourway[1], alpha = 0.5) +
    coord_cartesian(xlim = c(2000, 2019),
                    ylim = c(0, 700), clip = "on") +
    scale_colour_manual(values = colourway[c(8,5,2)]) +
    scale_x_continuous(breaks = seq(2000, 2019, 1),
                       labels = sub("^..", "'", seq(2000, 2019, 1)),
                       minor_breaks = NULL,
                       expand = expansion(c(0,0), c(1,1))) +
    labs(x = "Year", y = "Number\\nof species") +
    theme_bw() +
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          axis.line = element_line(),
          axis.title = element_text(face = 2),
          axis.title.y = element_text(angle = 0, hjust = 1))
)

rawCountsPlot

# Species unique to years -------------------------------------------------

## ONLINE TRADE SPP per year
onlineTempSpp <- temporalDataOnline %>% 
  group_by(year, sp) %>% 
  slice(n = 1) %>% 
  select(sp, year) %>% 
  mutate(source = "Online") %>% 
  ungroup()

## LEMIS SPP per year 
lemisTempSpp <- lemisData %>% 
  group_by(shipment_year, "sp" = lemisName) %>% 
  slice(n = 1) %>% 
  select("year" = shipment_year, sp) %>% 
  mutate(source = "LEMIS") %>% 
  ungroup()

## CITES SPP per year
citesTempSppList <- apply(citesData[,26:45], 2, function(x){
  unique(x[!is.na(x)])
})

citesTempSpp <- do.call(rbind, lapply(names(citesTempSppList), function(x){
  data.frame(year = sub("X", "", x), sp = citesTempSppList[[x]],
             source = "CITES")
}))

tempSppDF <- rbind(onlineTempSpp, lemisTempSpp, citesTempSpp)

# make sure that each species in each year appears once, essentially removing
# the source data meaning
tempSppDF <- tempSppDF %>% 
  group_by(year, sp) %>% 
  slice(n = 1)

i <- 0
uniTempSpp <- list()
for(y in unique(tempSppDF$year)){
  # y <- 2006
  i <- i+1
  noty <- unique(tempSppDF$sp[!tempSppDF$year == y])
  spy <- unique(tempSppDF$sp[tempSppDF$year == y])
  
  uniTempSpp[[i]] <- data.frame("year" = y, "nUniSpp" = sum(!spy %in% noty))
  
}#for end
uniTempSpp <- do.call(rbind, uniTempSpp)
uniTempSpp$year <- as.numeric(as.character(uniTempSpp$year))

(uniSppPlot <- uniTempSpp %>% 
    filter(year > 1999 & !year == 2020) %>% 
    ggplot() +
    geom_segment(aes(x = year, xend = year, y = 0, yend = nUniSpp),
                 size = 1.2, colour = colourway[2]) +
    geom_point(aes(x = year, y = nUniSpp),
               size = 3, colour = colourway[2]) +
    scale_x_continuous(breaks = seq(2000, 2019, 1),
                       labels = sub("^..", "'", seq(2000, 2019, 1)),
                       minor_breaks = NULL,
                       expand = expansion(c(0,0), c(1,1))) +
    theme_bw() +
    labs(x = "Year", y = "Number of\\nunique\\nspecies\\ntraded",
         colour = "") +
    # scale_colour_scico_d(palette = "roma") +
    scale_y_continuous(limits = c(0,55), breaks = seq(0, 60, 20)) +
    theme(legend.position = "none",
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(face = 2),
          axis.title.y = element_text(angle = 0, hjust = 1))
)

ggarrange(rawCountsPlot  +
            theme(plot.margin = margin(5,20,0,20)) +
            rremove("x.title") +
            rremove("x.text"),
          uniSppPlot +
            theme(plot.margin = margin(5,20,0,20)),
          ncol = 1,
          align = "v",
          heights = c(2,1),
          labels = c("A", "B"))

ggsave("./Figures/Temporal Plot.png", width = 200, height = 160,
       units = "mm")
ggsave("./Figures/Temporal Plot.pdf", width = 200, height = 160,
       units = "mm")
