
araData_trade <- read.csv("./Data/arachnid_species_data_traded.csv")

library(dplyr)
library(stringr)
library(tidytext)
library(ggplot2)
library(forcats)
library(ggpubr)
library(scico)

colourway <- scico(n = 9, palette = "imola")

temporalDataOnline <- read.csv("./Data/Temporal Online Data.csv")
generaWebRaw <- read.csv("./Data/Snapshot Online Genera Data Raw.csv")

temporalDataOnline <- temporalDataOnline %>% 
  filter(spORgen == "GENUS")

targetSites <- read.csv(file = "./Data/Target Websites Censored.csv")

generaWebRaw <- generaWebRaw %>% 
  left_join(targetSites %>% 
              select(lang, webID)) %>% 
  filter(lang == "ENG")

# termsAll <- data.frame(terms = c(temporalDataOnline$termsSurrounding,
#                                 generaWebRaw$termsSurrounding))
termsAll <- data.frame(terms = c(generaWebRaw$termsSurrounding))

termWords <- termsAll %>% 
  unnest_tokens(output = word, input = terms) 

termWords  <- termWords  %>%
  anti_join(stop_words)

termWordCount <- termWords %>%
  count(word, sort = TRUE)

# getting a raw count of the terms that frequency come up
termWordCount %>% 
  filter(nchar(word) > 2,
         !str_to_sentence(word) %in% unlist(str_split(araData_trade$allNames, ";|\\\\s")))

# colour list from https://en.wikipedia.org/wiki/List_of_colors:_N%E2%80%93Z

colours <- read.csv("./Data/WikipediaColourList.csv")
names(colours)[1] <- "name"

coloursClean <- str_replace_all(colours$name, "[[:punct:]]", " ") %>% 
  str_replace_all("  ", " ") %>% 
  str_split(" ") %>% 
  unlist()

coloursClean <- unique(coloursClean[!coloursClean == ""])

# extract from the common words those that are colours or colour adjacent
topColours <- termWordCount %>% 
  filter(nchar(word) > 2,
         word %in% coloursClean,
         # drop those are are questionably colours
         !word %in% c("space", "racing", "photo", "boy",
                      "bean", "blaze", "jungle", "mountain",
                      "dune", "web", "color", "rainforest", "tree",
                      "sea"))

# now we pull out those top colours from each of the terms that are linked to a genus
outList <- vector("list", dim(generaWebRaw)[1])
for(row in 1:dim(generaWebRaw)[1]){
  extColours <- unlist(str_extract_all(generaWebRaw$termsSurrounding[row],
                                       paste0(topColours$word, collapse = "|")))
  if(length(extColours) > 0){
    outList[[row]] <- data.frame(genus = generaWebRaw$sp[row],
                                 colour = extColours)
  }
}
generaColours <- do.call(rbind, outList)

colourWay <- colours %>% 
  mutate(name = str_to_lower(name)) %>% 
  filter(name %in% unique(generaColours$colour))

unique(generaColours$colour)

missingColours <- unique(generaColours$colour)[!unique(generaColours$colour) %in% colourWay$name]

colours %>% 
  mutate(name = str_to_lower(name)) %>% 
  filter(str_detect(name, paste0(missingColours, collapse = "|")))

"cobalt blue #0047AB"
"middle grey #8B8680"
"silver (metallic) #AAA9AD"
"light slate gray #778899"
"electric blue #7DF9FF"
"dark slate gray #2F4F4F"
"satin sheen gold #CBA135"

colourSchemeFull <- rbind(colourWay,
                          data.frame(
                            name = c("cobalt", 
                                     "grey" ,
                                     "metallic",
                                     "slate",
                                     "electric",
                                     "dark",
                                     "sheen",
                                     "chocolate"),
                            Hex = c(
                              "#0047AB",
                              "#8B8680",
                              "#AAA9AD",
                              "#778899",
                              "#7DF9FF",
                              "#2F4F4F",
                              "#CBA135",
                              "#7B3F00")
                          )) %>% 
  mutate(name = 
    factor(str_to_sentence(name), levels = str_to_sentence(c(
      "white",
      "silver",
      "metallic",
      "slate",
      "grey",
      "dark",
      "black",
      "brown",
      "chocolate",
      "gold",
      "sheen",
      "bronze",
      "sand",
      "tan",
      "chestnut",
      "burgundy",
      "red",
      "rose",
      "pink",
      "salmon",
      "orange",
      "lemon",
      "yellow",
      "green",
      "olive",
      "emerald",
      "electric",
      "blue",
      "sapphire",
      "cobalt",
      "violet",
      "purple"
      )))
  ) %>% 
  arrange(name) %>% 
  # making yellow and lemon more different
  mutate(Hex = ifelse(name == "Yellow", "#EBDF00", Hex))

generaColours %>% 
  group_by(genus, colour) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n = 300)

generaColours %>% 
  group_by(genus) %>% 
  summarise(nCol = length(unique(colour))) %>% 
  arrange(desc(nCol)) %>% 
  ungroup() %>% 
  summarise(
    mean(nCol),
    sd(nCol)/sqrt(length(nCol)),
    max(nCol)
  )

library(ggrepel)

(colourPlot <- generaColours %>% 
    mutate(colour = str_to_sentence(colour),
           colour = factor(colour, levels = levels(colourSchemeFull$name))) %>% 
    ggplot() +
    geom_bar(aes(x = fct_infreq(genus), fill = colour),
             position = "fill", width = 0.8, colour = NA
    ) +
    scale_fill_manual(values = colourSchemeFull$Hex) +
    coord_cartesian(clip = "off") +
    scale_x_discrete(position = "top") +
    scale_y_continuous(expand = c(0,0), limits = c(-0.1,1)) +
    labs(y = "Colour\\nproportion", fill = "Colour") +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 1, hjust = 1,
                                  face = 2),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.title = element_text(face = 2)
    )
)

generaColours %>% 
  count(genus) %>% 
  arrange(desc(n))

labelData <- generaColours %>% 
  count(genus) %>% 
  filter(n < 9) %>% 
  arrange(n)

# there are instances of older names being used, so we have to update the clade manually
# "Tarantula"  "Autana"     "Haplopelma" "Phlogius"   "Pachypus"

generaColoursClade <- generaColours %>% 
  left_join(araData_trade %>% 
              select(genus, clade) %>% 
              group_by(genus, clade) %>% 
              slice(n = 1)
  ) %>% 
  mutate(clade = case_when(
    genus %in% c("Tarantula","Autana","Haplopelma","Phlogius","Pachypus") ~ "Spider",
    TRUE ~ clade,
  )) %>% 
  mutate(clade = factor(clade, levels = c(
    "Spider",
    "Scorpion",
    "Uropygi")))

# generaColoursClade %>% 
#   filter(clade == "Scorpion") %>% 
#   mutate(genus = fct_infreq(genus)) %>% 
#   arrange(genus)

(countPlot <- generaColoursClade %>% 
    mutate(colour = str_to_sentence(colour)) %>%  
    ggplot() +
    geom_bar(aes(x = fct_infreq(genus)),
             width = 0.8, colour = NA, fill = "black") +
    # geom_point(data = generaColoursClade %>% 
    #              filter(clade == "Scorpion") %>% 
    #              count(genus),
    #            aes(x = genus, y = 400), alpha = 0.25) +
    geom_segment(data = generaColoursClade %>% 
                   filter(clade == "Scorpion") %>% 
                   count(genus),
               aes(x = genus, xend = genus, y = 400, yend = 0),
               linetype = 2, size = 0.5, alpha = 0.65) +
    annotate("text", x = "Hadogenes", y = 475, label = "Scorpions",
             vjust = 0.5, hjust = 1.15, fontface = 2, alpha = 0.65,
             colour = "black") +
    scale_y_continuous(limits = c(0, 1250)) +
    scale_fill_manual(values = colourway[c(6,2,8)]) +
    annotate("segment", x = "Protoiurus", xend = "Lychas", y = 400, yend = 400,
             linetype = 2, size = 0.5, alpha = 0.65) +
    # annotate("segment", x = "Viridasius", xend = "Davus", y = 500, yend = 500) +
    annotate("rect", xmin = "Viridasius", xmax = "Davus", ymin = 0, ymax = 1250,
             fill = "black", alpha = 0.25) +
    annotate("text", x = "Viridasius", y = 1250, label = "Fewer than 10 instances",
             vjust = 1, hjust = 1, fontface = 2, colour = "white",
             size = 6) +
    labs(y = "Count", fill = "") +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 1, hjust = 1,
                                  face = 2),
      axis.text.x = element_text(angle = 90, hjust = 0.5,
                                 vjust = 0.5,
                                 face = 3, size = 5.5)
    )
)

ggarrange(
  countPlot +
    guides(fill = guide_legend(ncol = 1)) +
    theme(plot.margin = margin(10, 5, 2, 20),
          axis.title.y = element_blank()),
  colourPlot +
    rremove("axis.text") +
    guides(fill = guide_legend(ncol = 1,
                               keyheight = unit(4, "mm"),
                               keywidth = unit(4, "mm"))) +
    theme(plot.margin = margin(2, 5, 5, 20),
          axis.title.y = element_blank())
  ,
  ncol = 1, align = "v",
  heights = c(0.8,1),
  common.legend = TRUE,
  legend = "right",
  labels = c("A", "B")
)

ggsave("./Figures/Colour by genus.png", dpi = 300,
       width = 300, height = 160, units = "mm")
ggsave("./Figures/Colour by genus.pdf", 
       width = 300, height = 160, units = "mm")

