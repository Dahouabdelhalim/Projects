
library(dplyr)
library(ggplot2)
library(scico)
library(stringr)
library(ggtext)
library(glue)
library(cowplot)
library(ggpubr)

araData_trade <- read.csv("./Data/arachnid_species_data_traded.csv")

lemisData <- read.csv(file = "./Data/TradeDatabases/LEMIS_arachnid_data.csv")

lemisData <- lemisData %>% 
  mutate(sp = str_to_sentence(paste(genus, species))) %>% 
  filter(!str_detect(sp, "sp\\\\.$") & !sp == "NA NA" & !sp == "Na na")

temporalDataOnline <- read.csv(file = "./Data/Temporal Online Data.csv")

firstYearTradedOnline <- temporalDataOnline %>% 
  filter(spORgen == "SPECIES", str_detect(str_trim(keyw), " "),
         !str_detect(keyw, "nomen dubium")) %>%
  group_by(keyw) %>% 
  summarise(tradeYear = min(year))

firstYearTradedLemis <- lemisData %>% 
  group_by(sp) %>% 
  summarise(tradeYear = min(shipment_year))

araData_trade$tradeYear <- apply(araData_trade, 1, function(x){
  if(!x["accName"] %in% unique(c(firstYearTradedLemis$sp,
                                 firstYearTradedOnline$keyw))){
    return(NA)
  } else {
    return(min(
      c(
        firstYearTradedOnline[which(firstYearTradedOnline$keyw == x["accName"]),]$tradeYear,
        firstYearTradedLemis[which(firstYearTradedLemis$sp == x["accName"]),]$tradeYear
      )
    ))
  }
})
# where they were detected in online trade but not the temp sample we know they
# appear at earliest in 2021
araData_trade$tradeYear[is.na(araData_trade$tradeYear) & araData_trade$onlineTradeEither] <- 2021

colourway <- scico(n = 9, palette = "imola")

label <- araData_trade %>% 
  mutate(tradeLag = tradeYear - year) %>% 
  filter(year > 1999, tradeYear > 1999) %>% 
  mutate(snapshotOnly = ifelse(tradeYear == 2021, "Snapshot\\nonly", "Temporal\\ndata"),
         snapshotOnly = factor(snapshotOnly, levels = c("Temporal\\ndata", "Snapshot\\nonly"))) %>% 
  filter(snapshotOnly == "Temporal\\ndata") %>% 
  summarise(mean = round(digits = 2, mean(tradeLag)),
            sd = round(digits = 2, sd(tradeLag)))

toTradeData <- araData_trade %>% 
  mutate(tradeLag = tradeYear - year) %>% 
  filter(year > 1999, tradeYear > 1999) %>% 
  mutate(snapshotOnly = ifelse(tradeYear == 2021, "Snapshot\\nonly", "Temporal\\ndata"),
         snapshotOnly = factor(snapshotOnly, levels = c("Temporal\\ndata", "Snapshot\\nonly"))
  ) %>% 
  arrange(family, year)

# clade summary
toTradeData %>% 
  filter(year > 1999, tradeYear > 1999) %>% 
  mutate(snapshotOnly = ifelse(tradeYear == 2021, "Snapshot\\nonly", "Temporal\\ndata"),
         snapshotOnly = factor(snapshotOnly, levels = c("Temporal\\ndata", "Snapshot\\nonly")),
         tradeLag = tradeYear - year,
         clade = factor(clade, levels = c(
           "Spider",
           "Scorpion",
           "Uropygi"))) %>%
  filter(snapshotOnly == "Temporal\\ndata") %>% 
  group_by(clade) %>% 
  summarise(
    mean = round(digits = 2, mean(tradeLag)),
    sd = round(digits = 2, sd(tradeLag))
  )

# family summary
toTradeData %>% 
  filter(year > 1999, tradeYear > 1999) %>% 
  mutate(snapshotOnly = ifelse(tradeYear == 2021, "Snapshot\\nonly", "Temporal\\ndata"),
         snapshotOnly = factor(snapshotOnly, levels = c("Temporal\\ndata", "Snapshot\\nonly")),
         tradeLag = tradeYear - year,
         clade = factor(clade, levels = c(
           "Spider",
           "Scorpion",
           "Uropygi"))) %>%
  filter(snapshotOnly == "Temporal\\ndata") %>% 
  group_by(family) %>% 
  summarise(
    mean = round(digits = 2, mean(tradeLag)),
    se = round(digits = 2, sqrt(var(tradeLag)/length(tradeLag)))
  )

(lagPlot <- toTradeData %>%
    mutate(clade = factor(clade, levels = c(
      "Spider",
      "Scorpion",
      "Uropygi"))) %>% 
    ggplot(aes()) +
    geom_segment(aes(x = year, xend = tradeYear,
                     y = reorder(accName, -year),
                     yend = reorder(accName, -year),
                     colour = clade)) +
    geom_point(aes(x = year, y = accName), colour = "black", size = 0.75,
               pch = 16) +
    geom_point(aes(x = tradeYear, y = accName), colour = "black", size = 0.75,
               pch = 16) +
    scale_shape_manual(values = c(NA, 16)) +
    scale_colour_manual(values = colourway[c(6,2,8)]) +
    coord_cartesian(clip = "off") +
    # facet_grid(snapshotOnly~., drop = TRUE, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = seq(2000, 2021, 1), minor_breaks = NULL,
                       limits = c(2000, 2021),
                       labels = c(paste0("'", str_pad(seq(00, 20, 1), 2, "left", 0)),
                                  "(Snapshot)\\n'21"),
                       position = "top") +
    theme_bw() +
    labs(x = "Year", y = "Species",
         colour = "", fill = "",
         tag = "(LEMIS only)\\n2000-2001") +
    theme(legend.position = c(0.40,0.45),
          legend.background = element_blank(),
          legend.key.size = unit(8, "mm"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.title = element_text(face = 2),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          # axis.text.y = element_markdown(),
          strip.background = element_blank(),
          # strip.text.y = element_text(face = 4, hjust = 0, angle = 0, vjust = 1)
          strip.text = element_blank(),
          plot.tag.position = c(0.075, 0.95),
          plot.tag = element_text(hjust = 0.5, vjust = 0, size = 6, lineheight = 0.85)) +
    guides(colour = guide_legend(override.aes = list(size = 8)))
)

(lagHist <- araData_trade %>% 
    filter(year > 1999, tradeYear > 1999) %>% 
    mutate(snapshotOnly = ifelse(tradeYear == 2021, "Snapshot\\nonly", "Temporal\\ndata"),
           snapshotOnly = factor(snapshotOnly, levels = c("Temporal\\ndata", "Snapshot\\nonly")),
           tradeLag = tradeYear - year,
           clade = factor(clade, levels = c(
             "Spider",
             "Scorpion",
             "Uropygi"))) %>%
    filter(snapshotOnly == "Temporal\\ndata") %>% 
    ggplot() +
    stat_count(aes(x = tradeLag, fill = clade),
               width = 0.75) +
    annotate("segment", x = label$mean, xend = label$mean,
             y = 0, yend = 13, linetype = 2, alpha = 0.5) +
    geom_text(data = label, aes(x = mean+0.25, y = 13,
                                label =  paste0("Mean ±SD\\n",
                                                mean, " ±", sd)),
              hjust = 0, vjust = 1, lineheight = 0.85, fontface = 4, size = 3)+
    scale_fill_manual(values = colourway[c(6,2)]) +
    scale_x_continuous(breaks = seq(0,21,2),
                       minor_breaks = NULL) +
    scale_y_continuous(minor_breaks = NULL) +
    theme_void() +
    labs(y = "Count",
         x = "<b><span style='font-size:12pt'>Lag between description and detection(years)</span></b> 
         <br>
         <span style='font-size:10pt'>[excluding snapshot detections]</span>") +
    theme(legend.position = "none",
          legend.text = element_markdown(),
          legend.title = element_text(face = 2),
          panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_blank(),
          axis.line = element_line(),
          axis.title.y = element_text(angle = 0, hjust = 1, face = 2, vjust = 1,
                                      margin = margin(0,5,0,0)),
          axis.title.x = element_markdown(),
          axis.text = element_text(colour = "grey25", size = 9,
            margin = margin(2,2,2,2)
          ),
          axis.ticks = element_line(),
          axis.ticks.length = unit(1, "mm"))
)

(comboPlot <-
    ggdraw() +
    draw_plot(lagPlot) +
    draw_plot(lagHist, x = 0.01, y = 0.01, width = 0.6, height = 0.6))

ggsave("./Figures/Desc to trade plot.png", width = 200, height = 140,
       dpi = 300, units = "mm")  
ggsave("./Figures/Desc to trade plot.pdf", width = 200, height = 140,
       units = "mm")

# overall summary
araData_trade %>% 
  # filter(family %in% tradedFamilies$family) %>% 
  summarise(post1999desc = sum(year > 1999),
            post1999trade = sum(year > 1999 & tradeYear > 1999, na.rm = TRUE)) %>%
  mutate(per1999trade = post1999trade / post1999desc *100)

# family summary of numbers described
araData_trade %>% 
  group_by(clade) %>% 
  summarise(post1999desc = sum(year > 1999),
            post1999trade = sum(year > 1999 & tradeYear > 1999, na.rm = TRUE)) %>%
  mutate(per1999trade = post1999trade / post1999desc *100)

# family summary of numbers described
araData_trade %>% 
  filter(family %in% c("Buthidae", "Theraphosidae")) %>%
  group_by(family) %>% 
  summarise(post1999desc = sum(year > 1999),
            post1999trade = sum(year > 1999 & tradeYear > 1999, na.rm = TRUE)) %>%
  mutate(per1999trade = post1999trade / post1999desc *100)

# -------------------------------------------------------------------------
