
lemisAraData <- read.csv(file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv",
                         stringsAsFactors = FALSE)
lemisCodes <- read.csv(file = "./Data/TradeDatabases/lemis_codes.csv")

library(dplyr)
library(stringr)
library(ggplot2)
library(scico)
library(countrycode)
library(cowplot)
library(ggpubr)

lemisAraData <- lemisAraData %>% 
  mutate(lemisName = str_to_sentence(paste(genus, species)))

lemisCodes %>% 
  filter(field == "description")

# Limit summaries to items that are actually representative of entire individuals
# item filter == description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO"))

# Items seized ------------------------------------------------------------
lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>% 
  mutate(seized = disposition == "S") %>% 
  group_by(seized) %>% 
  summarise(nItems = sum(quantity)) %>% 
  ungroup() %>% 
  mutate(totItems = sum(nItems),
         per = nItems / totItems *100)

# Commercial trade ------------------------------------------------------------
lemisCodes %>% 
  filter(field == "purpose")

lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>% 
  mutate(nonComm = ifelse(purpose %in% c("M", "S", "Y"), "NonComm", "Comm")) %>% 
  group_by(nonComm) %>% 
  summarise(nItems = sum(quantity)) %>% 
  ungroup() %>% 
  mutate(totItems = sum(nItems),
         per = nItems / totItems *100)


# Wild capture rates ------------------------------------------------------
lemisCodes %>% 
  filter(field == "source")

lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>% 
  mutate(originSimp = case_when(
    source %in% c("C", "D", "F", "R") ~ "Captive",
    source == "W" ~ "Wild",
    TRUE ~ "Other"
  )) %>% 
  group_by(originSimp) %>% 
  summarise(nItems = sum(quantity)) %>% 
  ungroup() %>% 
  mutate(totItems = sum(nItems),
         per = nItems / totItems *100)

# Wild capture rates per LEMIS genus --------------------------------------

dataForPlot <- lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>% 
  mutate(originSimp = case_when(
    source %in% c("C", "D", "F", "R") ~ "Captive",
    source == "W" ~ "Wild",
    TRUE ~ "Other"
  )) %>% 
  # dropping listings that aren't genera
  filter(!genus %in% c("Non-CITES entry"),
         !is.na(genus)) %>% 
  group_by(genus, originSimp) %>% 
  summarise(nItems = sum(quantity)) %>% 
  ungroup() %>% 
  tidyr::complete(genus, originSimp,
                  fill = list(nItems = 0)) %>% 
  group_by(genus) %>% 
  mutate(totGenItems = sum(nItems)) %>% 
  ungroup() %>% 
  mutate(per = nItems / totGenItems *100,
         facetTot = case_when(
           totGenItems >= 100000 ~ "More than 100,000 individuals",
           totGenItems < 100000 & totGenItems >= 10000 ~ "Between 100,000 and 10,000 individuals",
           totGenItems < 10000 & totGenItems >= 1000 ~ "Between 10,000 and 1,000 individuals",
           totGenItems < 1000 & totGenItems >= 100 ~ "Between 1,000 and 100 individuals",
           totGenItems < 100 ~ "Fewer than 100 individuals",
         ),
         facetTot = factor(facetTot, levels = c(
           "More than 100,000 individuals",
           "Between 100,000 and 10,000 individuals",
           "Between 10,000 and 1,000 individuals",
           "Between 1,000 and 100 individuals",
           "Fewer than 100 individuals"
         ))) %>% 
  arrange(desc(nItems)) %>% 
  mutate(genus = factor(genus, levels = unique(genus)),
         # invert the captive values to contrast against wild
         nItems = ifelse(originSimp == "Captive", -nItems, nItems),
         # reorder origin legend
         originSimp = factor(originSimp, levels = c("Wild", "Captive", "Other")),
         # remove the % for other, too cluttered and small
         per = ifelse(originSimp == "Other", NA, per))

plotList <- list()
for(face in unique(dataForPlot$facetTot)){
  
  generaPlot <- dataForPlot %>% 
    filter(facetTot == face) %>% 
    mutate(labelloc = ifelse(originSimp == "Captive", min(nItems) *1.1,
                             max(nItems) *1.1 )) %>% 
    ggplot() +
    geom_col(aes(x = genus, y = nItems, fill = originSimp)) +
    geom_vline(xintercept = seq(1.5, 200.5, 1), linetype = 2, alpha = 0.2) +
    labs(x = "Genus", y = "Number of\\nindividuals",
         colour = "Origin", fill = "Origin") +
    scale_fill_manual(values = scico(n = 3, palette = "roma")[c(1,3,2)]) +
    scale_colour_manual(values = scico(n = 3, palette = "roma")[c(1,3,2)]) +
    facet_wrap(facetTot~., ncol = 1, scales = "free") +
    theme_bw() +
    theme(
      legend.position = c(0.95, 0.85),
      legend.background = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5, face = 3),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(face = 2, size = 14),
      axis.title.y = element_text(angle = 0, hjust = 1, size = 14),
      plot.title = element_text(face = 4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.title = element_text(face = 2),
      strip.background = element_blank(),
      strip.text = element_text(hjust = 0, face = 4, size = 14)
    )
  
  if(face == "More than 100,000 individuals"){
    extraLabel <- dataForPlot %>% 
      filter(originSimp == "Wild", facetTot == face) %>% 
      summarise(meanPer = mean(per),
                medPer = median(per),
                minPer = min(per),
                maxPer = max(per),
                n = length(per)) 
    
    generaPlot <- generaPlot +
      geom_text(aes(x = genus, y = labelloc,
                    label = ifelse(!is.na(per),
                                   format(round(per, digits = 2), nsmall = 2), NA),
                    colour = originSimp),
                fontface = 2
      ) +
      labs(caption = paste0("% from wild: Mean = ", format(round(extraLabel$meanPer, digits = 2), nsmall = 2),
                            ", Median =", format(round(extraLabel$medPer, digits = 2), nsmall = 2),
                            ", Min = ", format(round(extraLabel$minPer, digits = 2), nsmall = 2),
                            ", Max = ", format(round(extraLabel$maxPer, digits = 2), nsmall = 2),
                            ", n = ", extraLabel$n)) +
      scale_y_continuous(breaks = seq(-1000000, 1000000, 200000),
                         labels = function(x) format(abs(seq(-1000000, 1000000, 200000)),
                                                     scientific = FALSE))
  } else if(face == "Between 100,000 and 10,000 individuals"){
    extraLabel <- dataForPlot %>% 
      filter(originSimp == "Wild", facetTot == face) %>% 
      summarise(meanPer = mean(per),
                medPer = median(per),
                minPer = min(per),
                maxPer = max(per),
                n = length(per)) 
    
    generaPlot <- generaPlot +
      labs(caption = paste0("% from wild: Mean = ", format(round(extraLabel$meanPer, digits = 2), nsmall = 2),
                            ", Median =", format(round(extraLabel$medPer, digits = 2), nsmall = 2),
                            ", Min = ", format(round(extraLabel$minPer, digits = 2), nsmall = 2),
                            ", Max = ", format(round(extraLabel$maxPer, digits = 2), nsmall = 2),
                            ", n = ", extraLabel$n)) +
      scale_y_continuous(breaks = seq(-100000, 100000, 20000),
                         labels = function(x) format(abs(seq(-100000, 100000, 20000)),
                                                     scientific = FALSE))
  } else if(face == "Between 10,000 and 1,000 individuals"){
    extraLabel <- dataForPlot %>% 
      filter(originSimp == "Wild", facetTot == face) %>% 
      summarise(meanPer = mean(per),
                medPer = median(per),
                minPer = min(per),
                maxPer = max(per),
                n = length(per)) 
    
    generaPlot <- generaPlot +
      labs(caption = paste0("% from wild: Mean = ", format(round(extraLabel$meanPer, digits = 2), nsmall = 2),
                            ", Median =", format(round(extraLabel$medPer, digits = 2), nsmall = 2),
                            ", Min = ", format(round(extraLabel$minPer, digits = 2), nsmall = 2),
                            ", Max = ", format(round(extraLabel$maxPer, digits = 2), nsmall = 2),
                            ", n = ", extraLabel$n)) +
      scale_y_continuous(breaks = seq(-10000, 10000, 2000),
                         labels = function(x) format(abs(seq(-10000, 10000, 2000)),
                                                     scientific = FALSE))
  } else if(face == "Between 1,000 and 100 individuals"){
    extraLabel <- dataForPlot %>% 
      filter(originSimp == "Wild", facetTot == face) %>% 
      summarise(meanPer = mean(per),
                medPer = median(per),
                minPer = min(per),
                maxPer = max(per),
                n = length(per)) 
    
    generaPlot <- generaPlot +
      labs(caption = paste0("% from wild: Mean = ", format(round(extraLabel$meanPer, digits = 2), nsmall = 2),
                            ", Median =", format(round(extraLabel$medPer, digits = 2), nsmall = 2),
                            ", Min = ", format(round(extraLabel$minPer, digits = 2), nsmall = 2),
                            ", Max = ", format(round(extraLabel$maxPer, digits = 2), nsmall = 2),
                            ", n = ", extraLabel$n)) +
      scale_y_continuous(breaks = seq(-1000, 1000, 200),
                         labels = function(x) format(abs(seq(-1000, 1000, 200)),
                                                     scientific = FALSE))
    
  } else if(face == "Fewer than 100 individuals"){
    extraLabel <- dataForPlot %>% 
      filter(originSimp == "Wild", facetTot == face) %>% 
      summarise(meanPer = mean(per),
                medPer = median(per),
                minPer = min(per),
                maxPer = max(per),
                n = length(per)) 
    
    generaPlot <- generaPlot +
      labs(caption = paste0("% from wild: Mean = ", format(round(extraLabel$meanPer, digits = 2), nsmall = 2),
                            ", Median =", format(round(extraLabel$medPer, digits = 2), nsmall = 2),
                            ", Min = ", format(round(extraLabel$minPer, digits = 2), nsmall = 2),
                            ", Max = ", format(round(extraLabel$maxPer, digits = 2), nsmall = 2),
                            ", n = ", extraLabel$n)) +
      scale_y_continuous(breaks = seq(-100, 100, 20),
                         labels = function(x) format(abs(seq(-100, 100, 20)),
                                                     scientific = FALSE))
  }
  
  plotList[[face]] <- generaPlot
  
  print(generaPlot)
  ggsave(file = paste0("./Figures/Genera wild cap_",
                       which(unique(dataForPlot$facetTot) == face),
                       "_",
                       face, "_.png"),
         width = 340, height = 220,
         dpi = 300, units = "mm")
}

dataForPlot %>% 
  filter(originSimp == "Wild") %>% 
  group_by(facetTot) %>% 
  summarise(meanPer = mean(per),
            medPer = median(per),
            minPer = min(per),
            maxPer = max(per),
            n = length(per))

dataForPlot %>% 
  filter(originSimp == "Wild") %>% 
  # group_by(facetTot) %>% 
  summarise(meanPer = mean(per),
            medPer = median(per),
            minPer = min(per),
            maxPer = max(per),
            n = length(per),
            impacted = sum(per > 0),
            fullwild = sum(per == 100),
            fullcap = sum(per == 0)
  )


# -------------------------------------------------------------------------

lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>%
  # filter(genus == "Mesobuthus") %>% 
  group_by(purpose) %>% 
  summarise(nItems = sum(quantity)) %>% 
  ungroup() %>% 
  mutate(totItems = sum(nItems),
         per = nItems / totItems *100)

lemisCodes %>% 
  filter(field == "purpose")

lemisAraData %>% 
  filter(country_origin == "CL") %>% 
  pull(lemisName) %>% 
  unique()

# Mapping wild vs captive -------------------------------------------------

mappingWildDataRaw <- lemisAraData %>% 
  filter(description %in% c("BOD", "EGL", "DEA", "LIV", "SPE", "SKI", "TRO")) %>%
  mutate(originSimp = case_when(
    source %in% c("C", "D", "F", "R") ~ "Captive",
    source == "W" ~ "Wild",
    TRUE ~ "Other"
  ))

worldData <- map_data("world")

mappingWildDataRaw$country_name <- countrycode::countrycode(mappingWildDataRaw$country_origin,
                                                         origin = "iso2c",
                                                         destination = "country.name",
                                                         custom_match = c("AN" = "Netherlands Antilles"))

# VS and XX unknown
sort(unique(worldData$region))
mappingWildDataRaw$country_name[!mappingWildDataRaw$country_name %in% worldData$region]

mappingWildDataSums <- mappingWildDataRaw %>% 
  mutate(region = case_when(
    country_name  == "Netherlands Antilles" ~ "Saba",
    country_name  == "Congo - Kinshasa" ~ "Democratic Republic of the Congo",
    country_name  == "Congo - Brazzaville" ~ "Republic of the Congo",
    country_name  == "Czechia" ~ "Czech Republic",
    country_name  == "Micronesia (Federated States of)" ~ "Micronesia",
    country_name  == "United Kingdom" ~ "UK",
    country_name  == "Hong Kong SAR China" ~ "China",
    country_name  == "St. Kitts & Nevis" ~ "Saint Kitts",
    country_name  == "St. Lucia" ~ "Saint Lucia",
    country_name  == "Myanmar (Burma)" ~ "Myanmar",
    country_name  == "São Tomé & Príncipe" ~ "Sao Tome and Principe",
    country_name  == "Trinidad & Tobago" ~ "Trinidad",
    country_name  == "United States" ~ "USA",
    country_name  == "St. Vincent & Grenadines" ~ "Saint Vincent",
    country_name  == "British Virgin Islands" ~ "Virgin Islands",
    TRUE ~ country_name
  )) %>% 
  group_by(region, originSimp) %>%
  filter(quantity > 0) %>% # remove odd instance of trade listed but for zero items
  summarise(nItems = sum(quantity)) %>% 
  mutate(totItems = sum(nItems),
         per = nItems / totItems *100) %>% 
  filter(originSimp == "Wild")

mappingWildDataSums %>% 
  filter(is.na(region))

mappingWildDataSums %>% 
  arrange(desc(nItems))

# counts ------------------------------------------------------------------

(lemis_countPlot <- worldData %>% 
  left_join(mappingWildDataSums) %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = nItems)) +
  scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                   breaks = c(seq(0, 600000, 100000), max(mappingWildDataSums$nItems)),
                   labels = scales::comma,
                   na.value = "grey85") +
  coord_map("mollweide", xlim = c(-180, 180)) +
  labs(fill = "Number of individuals\\nfrom wild") +
  # theme_bw() +
  theme_void() +
  theme(
    legend.position = c(0.185, 0.32),
    legend.title = element_text(face = 2, hjust = 1),
    legend.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  guides(fill = guide_colourbar(title.position = "bottom",
                                title.hjust = 0,
                                label.position = "left",
                                barwidth = unit(0.5, "lines"),
                                barheight = unit(10, "lines")))
  )

# percentage --------------------------------------------------------------

(lemis_percentPlot <- worldData %>% 
  left_join(mappingWildDataSums) %>% 
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group, fill = per)) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85")+
  coord_map("mollweide", xlim = c(-180, 180)) +
  labs(fill = "% from\\nwild") +
  # theme_bw() +
  theme_void() +
  theme(
    legend.position = c(0.155, 0.32),
    legend.title = element_text(face = 2, hjust = 1),
    legend.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  guides(fill = guide_colourbar(title.position = "bottom",
                                title.hjust = 0,
                                label.position = "left",
                                barwidth = unit(0.5, "lines"),
                                barheight = unit(10, "lines")))
  )


# number of genera -------------------------------------------------------

mappingWildDataGenera <- mappingWildDataRaw %>% 
  filter(!genus %in% c("Non-CITES entry"),
         !is.na(genus)) %>% 
  mutate(region = case_when(
    country_name  == "Netherlands Antilles" ~ "Saba",
    country_name  == "Congo - Kinshasa" ~ "Democratic Republic of the Congo",
    country_name  == "Congo - Brazzaville" ~ "Republic of the Congo",
    country_name  == "Czechia" ~ "Czech Republic",
    country_name  == "Micronesia (Federated States of)" ~ "Micronesia",
    country_name  == "United Kingdom" ~ "UK",
    country_name  == "Hong Kong SAR China" ~ "China",
    country_name  == "St. Kitts & Nevis" ~ "Saint Kitts",
    country_name  == "St. Lucia" ~ "Saint Lucia",
    country_name  == "Myanmar (Burma)" ~ "Myanmar",
    country_name  == "São Tomé & Príncipe" ~ "Sao Tome and Principe",
    country_name  == "Trinidad & Tobago" ~ "Trinidad",
    country_name  == "United States" ~ "USA",
    country_name  == "St. Vincent & Grenadines" ~ "Saint Vincent",
    country_name  == "British Virgin Islands" ~ "Virgin Islands",
    TRUE ~ country_name
  )) %>% 
  filter(quantity > 0) %>% # remove odd instance of trade listed but for zero items
  group_by(region, originSimp, genus) %>% 
  slice(n = 1) %>% 
  group_by(region, originSimp) %>%
  summarise(nGen = n()) %>%
  filter(originSimp == "Wild")

mappingWildDataGenera %>% 
  arrange(desc(nGen))

(lemis_generaPlot <- worldData %>% 
    left_join(mappingWildDataGenera) %>% 
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group, fill = nGen)) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     breaks = c(1, seq(10, 40, 10), max(mappingWildDataGenera$nGen) ),
                     na.value = "grey85")+
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of wild\\ncaught genera") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.185, 0.32),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

ggarrange(lemis_countPlot +
            theme(
              plot.margin = margin(0,0,0,0)
            ), 
          lemis_percentPlot +
            theme(
              plot.margin = margin(0,0,0,0)
            ),
          lemis_generaPlot +
            theme(
              plot.margin = margin(0,0,0,0)
            ),
          ncol = 1,
          label.y = 0.95,
          labels = "AUTO")

ggsave("./Figures/LEMIS Wild maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/LEMIS Wild maps.pdf",
       width = 150, height = 250, units = "mm")

# Overall mapping ---------------------------------------------------------

library(rgdal)
# library(sf)
# library(rmapshaper)

araData_trade <- read.csv("./Data/arachnid_species_data_traded.csv")

worldData <- readOGR("./Data/Spider_dists/TM_WORLD_BORDERS-0.3.shp")

allSpiderDBF <- read.csv("./Data/Spider_dists/All_spiders_f.csv",
                         encoding = "UTF8")

names(allSpiderDBF)[1] <- "accName"
# length(unique(allSpiderDBF$accName))

unique(allSpiderDBF$genus[!allSpiderDBF$NAME %in%
                            unique(worldData$NAME)])
unique(allSpiderDBF$NAME[!allSpiderDBF$NAME %in%
                            unique(worldData$NAME)])

allSpiderDBF$NAME[allSpiderDBF$NAME == "ÃžÃªÃ²and Islands"] <- "Ã…land Islands"

worldFort <- fortify(worldData, region = "NAME")

allSpiderDBF %>%
  select(genus, "id" = NAME) %>%
  filter(!genus == "") %>%
  left_join(araData_trade %>%
              select(genus, anyMatchTraded) %>%
              group_by(genus) %>%
              summarise(traded = sum(anyMatchTraded) > 0), by = "genus") %>%
  group_by(id, genus) %>%
  slice(n = 1) %>%
  group_by(id) %>%
  summarise(nGen = n_distinct(genus),
            nTraded = sum(traded)) %>%
  ungroup() %>%
  mutate(perTrade = nTraded / nGen *100)

worldFort_data <- worldFort %>%
  left_join(allSpiderDBF %>%
              select(genus, "id" = NAME) %>%
              filter(!genus == "") %>% 
              left_join(araData_trade %>% 
                          select(genus, anyMatchTraded) %>% 
                          group_by(genus) %>% 
                          summarise(traded = sum(anyMatchTraded) > 0), by = "genus") %>% 
              group_by(id, genus) %>% 
              slice(n = 1) %>% 
              group_by(id) %>% 
              summarise(nGen = n_distinct(genus),
                        nTraded = sum(traded)) %>% 
              ungroup() %>% 
              mutate(perTrade = nTraded / nGen *100), by = "id")

# Raw richness spiders --------------------------------------------------------

(richnessPlot <- ggplot(worldFort_data) +
  geom_polygon(aes(x = long, y = lat, group = group,
              fill = nGen), colour = NA) +
  scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                   na.value = "grey85",
                   breaks = c(seq(0,250,50), max(worldFort_data$nGen, na.rm = TRUE))
                   ) +
  coord_map("mollweide", xlim = c(-180, 180)) +
  labs(fill = "Number of\\ngenera") +
  # theme_bw() +
  theme_void() +
  theme(
    legend.position = c(0.15, 0.32),
    legend.title = element_text(face = 2, hjust = 1),
    legend.background = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  guides(fill = guide_colourbar(title.position = "bottom",
                                title.hjust = 0,
                                label.position = "left",
                                barwidth = unit(0.5, "lines"),
                                barheight = unit(10, "lines")))
)


# Richness traded spiders ---------------------------------------------------------

(tradedPlot <- ggplot(worldFort_data) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = nTraded), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85",
                    breaks = c(seq(0,75,25), max(worldFort_data$nTraded, na.rm = TRUE))
                    ) +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "Number of\\ntraded genera") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.180, 0.32),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)

# % traded spiders ----------------------------------------------------------------

(pertradedPlot <- ggplot(worldFort_data) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = perTrade), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85") +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "% genera\\ntraded") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.182, 0.34),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)


ggarrange(
  richnessPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ), 
  tradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  pertradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  ncol = 1,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Distribution genus spider trade maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/Distribution genus spider trade maps.pdf",
       width = 150, height = 250, units = "mm")

# Species richness and traded supplement spiders ----------------------------------

allSpiderDBF %>%
  select(accName, "id" = NAME) %>%
  filter(!accName == "") %>%
  left_join(araData_trade %>%
              select(accName, anyMatchTraded)) %>%
  group_by(id, accName) %>%
  slice(n = 1) %>%
  filter(!is.na(anyMatchTraded)) %>% 
  group_by(id) %>%
  summarise(nSpp = n_distinct(accName),
            nTraded = sum(anyMatchTraded == TRUE)) %>%
  ungroup() %>%
  mutate(perTrade = nTraded / nSpp *100)

worldFort_dataSP <- worldFort %>%
  left_join(allSpiderDBF %>%
              select(accName, "id" = NAME) %>%
              filter(!accName == "") %>%
              left_join(araData_trade %>%
                          select(accName, anyMatchTraded)) %>%
              group_by(id, accName) %>%
              slice(n = 1) %>%
              filter(!is.na(anyMatchTraded)) %>% 
              group_by(id) %>%
              summarise(nSpp = n_distinct(accName),
                        nTraded = sum(anyMatchTraded == TRUE)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nSpp *100), by = "id")

worldFort_dataSP %>% 
  group_by(id) %>% 
  slice(n = 1) %>% 
  arrange(desc(perTrade))

(sp_richnessPlot <- ggplot(worldFort_dataSP) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nSpp), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,4000,1000), max(worldFort_dataSP$nSpp, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\nspecies") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.15, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(sp_tradedPlot <- ggplot(worldFort_dataSP) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = nTraded), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85",
                    breaks = c(seq(0,120,20), max(worldFort_dataSP$nTraded, na.rm = TRUE))
   ) +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "Number of\\ntraded species") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.180, 0.35),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)

(sp_pertradedPlot <- 
    worldFort_dataSP %>% 
    mutate(perTrade = ifelse(perTrade > 15, 15, perTrade)) %>% 
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = perTrade), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     limits = c(0,15),
                     labels = c(seq(0,10,5), ">15")
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "% species\\ntraded") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.160, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)

ggarrange(
  sp_richnessPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ), 
  sp_tradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  sp_pertradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  ncol = 1,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Distribution spider species trade maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/Distribution spider species trade maps.pdf",
       width = 150, height = 250, units = "mm")


# Top genera traded spiders -------------------------------------------------------

twoPlusSppGenera <- araData_trade %>%
  select(genus, anyMatchTraded) %>%
  group_by(genus) %>%
  summarise(twoTraded = sum(anyMatchTraded) >= 2) %>% 
  filter(twoTraded) %>% 
  pull(genus)

worldFort_dataTopGen <- worldFort %>%
  left_join(allSpiderDBF %>%
              select(accName, genus, "id" = NAME) %>%
              filter(!genus == "") %>%
              left_join(araData_trade %>%
                          select(accName, genus, anyMatchTraded), by = c("accName", "genus")) %>%
              filter(genus %in% twoPlusSppGenera) %>% # so we are restricting the calc to just the genera with more that two species traded
              group_by(id, accName) %>%
              slice(n = 1) %>%
              group_by(id) %>%
              summarise(nSpp = n_distinct(accName),
                        nTraded = sum(anyMatchTraded == TRUE)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nSpp *100), by = "id")

worldFort_dataTopGen %>% 
  group_by(id) %>% 
  slice(n = 1) %>% 
  arrange(desc(perTrade))

(topGen_richnessPlot <- ggplot(worldFort_dataTopGen) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nSpp), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,500,100), max(worldFort_dataTopGen$nSpp, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\nspecies") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.15, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(topGen_tradedPlot <- ggplot(worldFort_dataTopGen) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nTraded), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,120,20), max(worldFort_dataTopGen$nTraded, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\ntraded species") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.180, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(topGen_pertradedPlot <- 
    worldFort_dataTopGen %>% 
    mutate(perTrade = ifelse(perTrade > 60, 60, perTrade)) %>% 
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = perTrade), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     limits = c(0,60),
                     labels = c(seq(0,50,10), ">60")
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "% species\\ntraded") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.160, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

ggarrange(
  topGen_richnessPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ), 
  topGen_tradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  topGen_pertradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  ncol = 1,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Distribution topGenera species trade maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/Distribution topGenera species trade maps.pdf",
       width = 150, height = 250, units = "mm")

# Top family breakdowns ----------------------------------------------------

# Theraphosidae -----------------------------------------------------------

worldFort_dataThera <- worldFort %>%
  left_join(allSpiderDBF %>%
              select(accName, "id" = NAME) %>%
              filter(!accName == "") %>%
              left_join(araData_trade %>%
                          select(accName, family, anyMatchTraded)) %>%
              group_by(id, accName) %>%
              slice(n = 1) %>%
              filter(!is.na(anyMatchTraded),
                     family == "Theraphosidae") %>% 
              group_by(id) %>%
              summarise(nSpp = n_distinct(accName),
                        nTraded = sum(anyMatchTraded == TRUE)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nSpp *100), by = "id")

(thera_richnessPlot <- ggplot(worldFort_dataThera) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nSpp), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,180,20), max(worldFort_dataThera$nSpp, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\nTheraphosids") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.15, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(thera_tradedPlot <- ggplot(worldFort_dataThera) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nTraded), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,80,20), max(worldFort_dataThera$nTraded, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\ntraded Theraphosids") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.180, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(thera_pertradedPlot <- ggplot(worldFort_dataThera) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = perTrade), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85") +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "% Theraphosids\\ntraded") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.160, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

ggarrange(
  thera_richnessPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ), 
  thera_tradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  thera_pertradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  ncol = 1,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Distribution Theraphosids trade maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/Distribution Theraphosids trade maps.pdf",
       width = 150, height = 250, units = "mm")


# Top six spider % --------------------------------------------------------


worldFort_dataThera <- worldFort %>%
  left_join(allSpiderDBF %>%
              select(accName, "id" = NAME) %>%
              filter(!accName == "") %>%
              left_join(araData_trade %>%
                          select(accName, family, anyMatchTraded)) %>%
              group_by(id, accName) %>%
              slice(n = 1) %>%
              filter(!is.na(anyMatchTraded),
                     family == "Theraphosidae") %>% 
              group_by(id, family) %>%
              summarise(nSpp = n_distinct(accName),
                        nTraded = sum(anyMatchTraded == TRUE)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nSpp *100), by = "id")

topFam_pertradedPlotList <- list()
for(topSpiders in c("Theraphosidae",
                    "Sparassidae",
                    "Araneidae",
                    "Theridiidae",
                    "Lycosidae",
                    "Salticidae")){
  
  topFam_data <- worldFort %>%
    left_join(allSpiderDBF %>%
                select(accName, "id" = NAME) %>%
                filter(!accName == "") %>%
                left_join(araData_trade %>%
                            select(accName, family, anyMatchTraded)) %>%
                group_by(id, accName) %>%
                slice(n = 1) %>%
                filter(!is.na(anyMatchTraded),
                       family == topSpiders) %>%
                group_by(id, family) %>%
                summarise(nSpp = n_distinct(accName),
                          nTraded = sum(anyMatchTraded == TRUE)) %>%
                ungroup() %>%
                mutate(perTrade = nTraded / nSpp *100), by = "id")
  
  if(topSpiders == "Araneidae"){
    topFam_data <- topFam_data %>% 
      mutate(perTrade = ifelse(perTrade > 50, 50, perTrade))
  }
  if(topSpiders == "Salticidae"){
    topFam_data <- topFam_data %>% 
      mutate(perTrade = ifelse(perTrade > 20, 20, perTrade))
  }
    
  labelData <- araData_trade %>% 
    filter(family == topSpiders) %>% 
    summarise(nSpecies = n_distinct(accName),
              nTraded = sum(anyMatchTraded))
  
  topFam_pertradedPlotList[[topSpiders]] <- ggplot(topFam_data) +
      geom_polygon(aes(x = long, y = lat, group = group,
                       fill = perTrade), colour = NA) +
    # annotate("text", x = 0, y = -55, label = 
    #            paste0("Traded / Total\\n",
    #                   labelData$nTraded, " / ", labelData$nSpecies),
    #            hjust = 0.5, fontface = 2,
    #          lineheight = 0.85) +
      scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                       na.value = "grey85") +
      coord_map("mollweide", xlim = c(-180, 180)) +
    labs(x = paste0("Traded/Total: ", 
                           labelData$nTraded, "/", labelData$nSpecies)) +
      labs(fill = paste0("% ", topSpiders, "\\ntraded")) +
      # theme_bw() +
      theme_void() +
      theme(
        axis.title.x = element_text(face = 4, hjust = 0.5),
        legend.position = c(0.160, 0.35),
        legend.title = element_text(face = 2, hjust = 1),
        legend.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
      guides(fill = guide_colourbar(title.position = "bottom",
                                    title.hjust = 0,
                                    label.position = "left",
                                    barwidth = unit(0.5, "lines"),
                                    barheight = unit(10, "lines")))
  
  
}

ggarrange(
  topFam_pertradedPlotList[[1]]+
    theme(
      plot.margin = margin(5,2,0,2)
    ),
  topFam_pertradedPlotList[[2]]+
    theme(
      plot.margin = margin(5,2,0,2)
    ),
  topFam_pertradedPlotList[[3]] +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     limits = c(0,50),
                     labels = c(seq(0,40,10), ">50")
    )+
    theme(
      plot.margin = margin(5,2,0,2)
    ),
  topFam_pertradedPlotList[[4]]+
    theme(
      plot.margin = margin(0,2,5,2)
    ),
  topFam_pertradedPlotList[[5]]+
    theme(
      plot.margin = margin(0,2,5,2)
    ),
  topFam_pertradedPlotList[[6]]+
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     limits = c(0,20),
                     labels = c(seq(0,15,5), ">20")
    )+
    theme(
      plot.margin = margin(0,2,5,2)
    ),
  ncol = 3, nrow = 3,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Top spiders preTraded maps.png", dpi = 300,
       width = 360, height = 220, units = "mm")
ggsave("./Figures/Top spiders preTraded maps.pdf",
       width = 360, height = 220, units = "mm")


# Scorpion mapping --------------------------------------------------------

scorpLocData <- read.csv("./Data/Scorp_dists/scorp_finalOct2.csv")

sum(!unique(scorpLocData$NAME) %in% unique(worldFort$id))
unique(scorpLocData$NAME)[!unique(scorpLocData$NAME) %in% unique(worldFort$id)]

scorpLocData <- scorpLocData %>% 
  mutate(id = case_when(
    NAME == " Suriname" ~ "Suriname",
    NAME == "Iran" ~ "Iran (Islamic Republic of)",
    NAME == "C?te d'Ivoire" ~ "Cote d'Ivoire",
    NAME == "Uraguay" ~ "Uruguay",
    NAME == "Saba" ~ "Netherlands",
    NAME == "Saint-Eustache" ~ "Netherlands",
    NAME == "Canary Islands" ~ "Spain",
    NAME == " Switzerland" ~ "Switzerland",
    NAME == "Brunei Darussalam Darussalam" ~ "Brunei Darussalam",
    NAME == " Costa Rica" ~ "Costa Rica",
    NAME == "Saint-Vincent-et-les-Grenadines" ~ "Saint Vincent and the Grenadines",
    NAME == "Columbia" ~ "Colombia",
    TRUE ~ NAME))

scorpLocData %>%
  select("accName" = species, id) %>%
  filter(!accName == "") %>%
  left_join(araData_trade %>%
              select(accName, anyMatchTraded)) %>%
  group_by(id, accName) %>%
  slice(n = 1) %>%
  filter(!is.na(anyMatchTraded)) %>% 
  group_by(id) %>%
  summarise(nSpp = n_distinct(accName),
            nTraded = sum(anyMatchTraded == TRUE)) %>%
  ungroup() %>%
  mutate(perTrade = nTraded / nSpp *100)

worldFort_dataSPscorp <- worldFort %>%
  left_join(scorpLocData %>%
              select("accName" = species, id) %>%
              filter(!accName == "") %>%
              left_join(araData_trade %>%
                          select(accName, anyMatchTraded)) %>%
              group_by(id, accName) %>%
              slice(n = 1) %>%
              filter(!is.na(anyMatchTraded)) %>% 
              group_by(id) %>%
              summarise(nSpp = n_distinct(accName),
                        nTraded = sum(anyMatchTraded == TRUE)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nSpp *100), by = "id")

worldFort_dataSPscorp %>% 
  group_by(id) %>% 
  slice(n = 1) %>% 
  arrange(desc(perTrade))

(spScorp_richnessPlot <- ggplot(worldFort_dataSPscorp) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nSpp), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,250,50), max(worldFort_dataSPscorp$nSpp, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\nspecies") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.15, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(spScorp_tradedPlot <- ggplot(worldFort_dataSPscorp) +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = nTraded), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     breaks = c(seq(0,40,10), max(worldFort_dataSPscorp$nTraded, na.rm = TRUE))
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "Number of\\ntraded species") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.190, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

(spScorp_pertradedPlot <- 
    worldFort_dataSPscorp %>% 
    # mutate(perTrade = ifelse(perTrade > 15, 15, perTrade)) %>% 
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group,
                     fill = perTrade), colour = NA) +
    scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                     na.value = "grey85",
                     # limits = c(0,15),
                     # labels = c(seq(0,15,5), ">20")
    ) +
    coord_map("mollweide", xlim = c(-180, 180)) +
    labs(fill = "% species\\ntraded") +
    # theme_bw() +
    theme_void() +
    theme(
      legend.position = c(0.150, 0.35),
      legend.title = element_text(face = 2, hjust = 1),
      legend.background = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()) +
    guides(fill = guide_colourbar(title.position = "bottom",
                                  title.hjust = 0,
                                  label.position = "left",
                                  barwidth = unit(0.5, "lines"),
                                  barheight = unit(10, "lines")))
)

ggarrange(
  spScorp_richnessPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ), 
  spScorp_tradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  spScorp_pertradedPlot +
    theme(
      plot.margin = margin(0,0,0,0)
    ),
  ncol = 1,
  label.y = 0.95,
  labels = "AUTO")

ggsave("./Figures/Distribution scorpion species trade maps.png", dpi = 300,
       width = 150, height = 250, units = "mm")
ggsave("./Figures/Distribution scorpion species trade maps.pdf",
       width = 150, height = 250, units = "mm")


# Main figure spider and scoprions ----------------------------------------

spiderScorpLocData <- rbind(scorpLocData %>%
        mutate(genus = word(species, 1, 1)) %>% 
        select(genus, id) %>%
        filter(!genus == ""),
      allSpiderDBF %>%
        select(genus, "id" = NAME) %>%
        filter(!genus == ""))

spiderScorpLocData %>%
  left_join(araData_trade %>%
              select(genus, anyMatchTraded) %>%
              group_by(genus) %>%
              summarise(traded = sum(anyMatchTraded) > 0), by = "genus") %>%
  filter(!is.na(traded)) %>% 
  group_by(id, genus) %>%
  slice(n = 1) %>%
  group_by(id) %>%
  summarise(nGen = n_distinct(genus),
            nTraded = sum(traded)) %>%
  ungroup() %>%
  mutate(perTrade = nTraded / nGen *100)

worldFort_dataSpiderScorp <- worldFort %>%
  left_join(spiderScorpLocData %>%
              left_join(araData_trade %>%
                          select(genus, anyMatchTraded) %>%
                          group_by(genus) %>%
                          summarise(traded = sum(anyMatchTraded) > 0), by = "genus") %>%
              filter(!is.na(traded)) %>% 
              group_by(id, genus) %>%
              slice(n = 1) %>%
              group_by(id) %>%
              summarise(nGen = n_distinct(genus),
                        nTraded = sum(traded)) %>%
              ungroup() %>%
              mutate(perTrade = nTraded / nGen *100), by = "id")


# In text summaries -------------------------------------------------------

worldFort %>%
  left_join(spiderScorpLocData %>%
              left_join(araData_trade %>%
                          select(genus, anyMatchTraded) %>%
                          group_by(genus) %>%
                          summarise(traded = sum(anyMatchTraded) > 0), by = "genus") %>%
              filter(!is.na(traded)) %>% 
              group_by(id, genus) %>%
              slice(n = 1) %>%
              group_by(id) %>%
              summarise(nGen = n_distinct(genus),
                        nTraded = sum(traded))) %>% 
  ungroup() %>%
  mutate(perTrade = nTraded / nGen *100) %>% 
  group_by(id) %>%
  slice(n = 1) %>%
  arrange(desc(perTrade)) %>% 
  print(n  = 100)

# Raw richness --------------------------------------------------------

(richnessPlot_ALL <- ggplot(worldFort_dataSpiderScorp) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = nGen), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85",
                    breaks = c(seq(0,250,50), max(worldFort_dataSpiderScorp$nGen, na.rm = TRUE))
   ) +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "Number of\\ngenera") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.15, 0.32),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)


# Richness traded ---------------------------------------------------------

(tradedPlot_ALL <- ggplot(worldFort_dataSpiderScorp) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = nTraded), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85",
                    breaks = c(seq(0,75,25), max(worldFort_dataSpiderScorp$nTraded, na.rm = TRUE))
   ) +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "Number of\\ntraded genera") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.180, 0.32),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)

# % traded ----------------------------------------------------------------

(pertradedPlot_ALL <- ggplot(worldFort_dataSpiderScorp) +
   geom_polygon(aes(x = long, y = lat, group = group,
                    fill = perTrade), colour = NA) +
   scale_fill_scico(palette = "imola", begin = 0.1, end = 0.9,
                    na.value = "grey85") +
   coord_map("mollweide", xlim = c(-180, 180)) +
   labs(fill = "% genera\\ntraded") +
   # theme_bw() +
   theme_void() +
   theme(
     legend.position = c(0.182, 0.34),
     legend.title = element_text(face = 2, hjust = 1),
     legend.background = element_blank(),
     panel.border = element_blank(),
     axis.ticks.y = element_blank(),
     panel.grid.major.y = element_blank()) +
   guides(fill = guide_colourbar(title.position = "bottom",
                                 title.hjust = 0,
                                 label.position = "left",
                                 barwidth = unit(0.5, "lines"),
                                 barheight = unit(10, "lines")))
)

# Combine all -------------------------------------------------------------

ggarrange(
  lemis_countPlot +
    theme(
      plot.margin = margin(5,2,0,2)
    ), 
  lemis_generaPlot +
    theme(
      plot.margin = margin(5,2,0,2)
    ),
  lemis_percentPlot +
    theme(
      plot.margin = margin(5,2,0,2)
    ),
  richnessPlot_ALL +
    theme(
      plot.margin = margin(0,2,5,2)
    ), 
  tradedPlot_ALL +
    theme(
      plot.margin = margin(0,2,5,2)
    ),
  pertradedPlot_ALL +
    theme(
      plot.margin = margin(0,2,5,2)
    ),
  ncol = 3, nrow = 2,
  label.y = 0.95,
  # labels = c("A", "D", "B", "C", "E", "F"))
  labels = "AUTO")

ggsave("./Figures/All genera trade maps.png", dpi = 300,
       width = 360, height = 160, units = "mm")
ggsave("./Figures/All genera trade maps.pdf",
       width = 360, height = 160, units = "mm")