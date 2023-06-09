
araData_trade <- read.csv("./Data/arachnid_species_data_traded.csv")

library(ggplot2)
library(dplyr)
library(reshape2)
library(ggtext)
library(scico)
library(stringr)

colourway <- scico(n = 9, palette = "imola")

# Overall summary
# Figure - Source overlap
# Figure - Top traded families
# Figure - IUCN status
# Figure - Language summary

##### Figure - Species origin - DATA ISSUE


# Overall summary ---------------------------------------------------------
names(araData_trade)

araData_trade %>% 
  filter(anyMatchTraded) %>% 
  summarise(length(unique(genus)))

araData_trade %>% 
  group_by(genus) %>% 
  mutate(tradedGen = sum(anyMatchTraded) > 0 & n() > 1) %>% 
  filter(tradedGen) %>% 
  group_by(clade) %>% 
  summarise(nSpp = length(unique(accName))) %>% 
  left_join(araData_trade %>% 
              filter(anyMatchTraded) %>% 
              group_by(clade) %>% 
              summarise(tradedSpp = length(unique(accName)))) %>% 
  mutate(per = tradedSpp/nSpp *100)

araData_trade %>% 
  # filter(onlineTradeSnap_genus| onlineTradeSnap_genusAny|
  #            onlineTradeTemp_genus| onlineTradeTemp_genusAny|
  #            LEMIStrade_genus| LEMIStrade_genusAny|
  #            CITEStrade_genus| CITEStrade_genusAny) %>%
  filter(onlineTradeSnap_genus |
             onlineTradeTemp_genus |
             LEMIStrade_genus |
             CITEStrade_genus) %>%
  summarise(length(unique(genus)))

generaSummary <- araData_trade %>% 
  group_by(genus) %>% 
  summarise(
    nSpp = n()
  ) %>% 
  left_join(araData_trade %>% 
              filter(anyMatchTraded) %>% 
              group_by(genus) %>% 
              summarise(
                nTradedAny = n()
              )) %>% 
  left_join(araData_trade %>% 
              filter(extactMatchTraded) %>% 
              group_by(genus) %>% 
              summarise(
                nTradedExact = n()
              )) %>% 
  mutate(
    nTradedAny = ifelse(is.na(nTradedAny), 0, nTradedAny),
    nTradedExact = ifelse(is.na(nTradedExact), 0, nTradedExact),
    perTradedAny = nTradedAny / nSpp *100, 
    perTradedExact = nTradedExact / nSpp *100) %>% 
  filter(perTradedAny > 0 | perTradedExact)

generaSummary %>% 
  filter(!perTradedAny == perTradedExact)

write.csv(generaSummary, "./Tables/Genera Summary.csv", row.names = FALSE)

# species counts from online trade snapshot data
araData_trade %>% 
  summarise(across(c("onlineTradeSnap", "onlineTradeSnap_Any",
                     "onlineTradeTemp" , "onlineTradeTemp_Any", 
                     "onlineTradeEither" , "onlineTradeEither_Any", 
                     "LEMIStrade", "LEMIStrade_Any",
                     "CITEStrade", "CITEStrade_Any"),
                   ~ sum(.x, na.rm = TRUE)))

araData_trade %>% 
  summarise(across(c("extactMatchTraded", "anyMatchTraded"),
                   ~ sum(.x, na.rm = TRUE)))

# Source overlap ----------------------------------------------------------

library(UpSetR)

sum(araData_trade$anyMatchTraded)
sum(araData_trade$extactMatchTraded)

araData_trade %>% 
  filter(anyMatchTraded) %>% 
  group_by(clade) %>% 
  summarise(n())

upsetData <- araData_trade %>%
  select(accName,
         onlineTradeEither_Any,
         LEMIStrade_Any,
         CITEStrade_Any) %>% 
  mutate(across(.cols = c("onlineTradeEither_Any",
                          "LEMIStrade_Any",
                          "CITEStrade_Any"),
                .fns = as.numeric))

pdf(file = "./Figures/Base upset plot.pdf",
    width = 8, height = 4)
upset(upsetData)
dev.off()

# Top traded families ------------------------------------------------------

tradedFamilies <- araData_trade %>% 
  filter(anyMatchTraded) %>% 
  count(family, name = "spTraded")

familySummary <- araData_trade %>% 
  count(family) %>% 
  filter(family %in% tradedFamilies$family) %>% 
  left_join(tradedFamilies) %>% 
  mutate(perTraded = spTraded / n *100)

detectedAllSources <- araData_trade %>% 
  filter(family %in% familySummary$family[familySummary$perTraded > 5 | familySummary$spTraded > 20]) %>% 
  filter(onlineTradeSnap_Any, LEMIStrade_Any, CITEStrade_Any)


familySummaryPlotData <- familySummary %>% 
  left_join(araData_trade %>% 
              select(family, clade) %>% 
              group_by(family, clade) %>% 
              slice(n = 1) %>% 
              ungroup(),
            by = "family") %>% 
  mutate(
    family = case_when(
      family == "Araneidae" ~ "<span style='font-size:8pt'>Araneidae<br>
      </span><span style='font-size:6pt'>(Orb-weaver spiders)</span>",
      
      family == "Atracidae" ~ "<span style='font-size:8pt'>Atracidae<br>
      </span><span style='font-size:6pt'>(Australian funnel-web spiders)</span>",
      
      family == "Dipluridae" ~ "<span style='font-size:8pt'>Dipluridae<br>
      </span><span style='font-size:6pt'>(Curtain-web spiders)</span>",
      
      family == "Eresidae" ~ "<span style='font-size:8pt'>Eresidae<br>
      </span><span style='font-size:6pt'>(Velvet spiders)</span>",
      
      family == "Lycosidae" ~ "<span style='font-size:8pt'>Lycosidae<br>
      </span><span style='font-size:6pt'>(Wolf spiders)</span>",
      
      family == "Salticidae" ~ "<span style='font-size:8pt'>Salticidae<br>
      </span><span style='font-size:6pt'>(Jumping spiders)</span>",
      
      family == "Scorpionidae" ~ "<span style='font-size:8pt'>Scorpionidae<br>
      </span><span style='font-size:6pt'>(Burrowing scorpions)</span>",
      
      family == "Sparassidae" ~ "<span style='font-size:8pt'>Sparassidae<br>
      </span><span style='font-size:6pt'>(Huntsman spiders)</span>",
      
      family == "Tetragnathidae" ~ "<span style='font-size:8pt'>Tetragnathidae<br>
      </span><span style='font-size:6pt'>(Long jawed spiders)</span>",
      
      family == "Theraphosidae" ~ "<span style='font-size:8pt'>Theraphosidae<br>
      </span><span style='font-size:6pt'>(Tarantulas)</span>",
      
      family == "Theridiidae" ~ "<span style='font-size:8pt'>Theridiidae<br>
      </span><span style='font-size:6pt'>(Tangle-web spiders)</span>",
      
      TRUE ~ paste0("<span style='font-size:8pt'>", family, "</span>")
    )
  ) %>% 
  mutate(clade = factor(clade, levels = c(
    "Spider",
    "Scorpion",
    "Uropygi")))

familySummaryPlotData %>% 
  filter( (perTraded > 5 | spTraded > 20) & spTraded > 1) %>%
  ggplot() +
  geom_col(aes(x = reorder(family, perTraded), y = perTraded, fill = clade),
           width = 0.5) +
  geom_text(aes(x = reorder(family, perTraded), y = perTraded+0.5,
                label = paste0(round(perTraded, digits = 0), "%")),
            fontface = 2, vjust = 0.5, hjust = 0) +
  geom_text(aes(x = reorder(family, -perTraded), y = -4.5,
                label = n), fontface = 2,
            hjust = 1) +
  geom_text(aes(x = reorder(family, -perTraded), y = -3.5,
                label = spTraded), fontface = 2,
            hjust = 0) +
  geom_segment(
    data = data.frame(x = seq(0.75,27.75,1),
                      xend = seq(1.25, 28.25,1),
                      y = rep(-4,28)),
    aes(x = x, xend = xend, y = y, yend = y)
    
  ) +
  # total families
  annotate("richtext", x = -Inf, y = Inf,
           label = paste0(
             "<b><span style='font-size:12pt'>",dim(tradedFamilies)[1]," familes traded</b><br>
    </span><span style='font-size:8pt'>Only families with more than 1 species traded<br>
    and >5% of species or >20 species traded displayed</span>"),
    hjust = 1, vjust = -0.2,
    fill = NA, label.color = NA,
    label.padding = grid::unit(rep(0, 4), "pt"),
    lineheight = 0.9) +
  # highlighting the all matching
  annotate("richtext", x = 20, y = 45,
           label = paste0(
             "<span style='font-size:8pt'> Contains the<br></span>",
             "<b><span style='font-size:26pt'>",
             sum(detectedAllSources$family == "Theraphosidae"),
             "</b><br></span>",
             "<span style='font-size:8pt'>species detected<br>by all sources</span>"),
           hjust = 0.5, vjust = 1,
           fill = NA, label.color = NA,
           label.padding = grid::unit(rep(0, 4), "pt"),
           lineheight = 0.9) +
  # arrow connecting
  annotate("curve", xend = 27, x = 20.5, yend = 42, y = 45,
           curvature = 0.35, arrow = arrow(type = "closed", angle = 25,
                                          length = unit(4, "mm")), size = 1.2) +
  labs(x = "<b><span style='font-size:12pt'>Family<br>
      </span><span style='font-size:8pt'>(Commonly<br>known as)</span></b>",
      y = "Percentage\\nof species", fill = "Clade") +
  scale_y_continuous(limits = c(-5, 50)) +
  scale_fill_manual(values = colourway[c(6,2,8)]) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = c(0.95, 0.2),
    legend.title = element_text(face = 2, hjust = 1),
    legend.background = element_blank(),
    axis.title.x = element_text(face = 2),
    axis.text.x = element_text(face = 1),
    axis.text.y = element_markdown(),
    axis.title.y = element_markdown(angle = 0, vjust = 1, hjust = 0, face = 2),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank()) +
  guides(
    fill = guide_legend(label.position = "left"))

ggsave("./Figures/Top traded familes bar.png", dpi = 300,
       width = 300, height = 160, units = "mm")
ggsave("./Figures/Top traded familes bar.pdf", dpi = 300,
       width = 300, height = 160, units = "mm")

familySummaryPlotData %>% 
  ggplot() +
  geom_point(aes(x = n, y = spTraded, colour = clade, size = perTraded,
                 shape = clade)) +
  # top labels
  geom_richtext(data = familySummaryPlotData %>% 
              filter(spTraded > 20),
              aes(x = n, y = spTraded, label = family),
              hjust = 0.5, vjust = 1,
              fill = NA, label.color = NA,
              label.padding = grid::unit(rep(0, 4), "pt"),
              lineheight = 0.9) +
  # top 
  geom_richtext(data = familySummaryPlotData %>% 
              filter(perTraded > 20),
              aes(x = n, y = spTraded, label = family),
              hjust = 0.5, vjust = 1,
              fill = NA, label.color = NA,
              label.padding = grid::unit(rep(0, 4), "pt"),
              lineheight = 0.9) +
  # total families
  annotate("richtext", x = 10000, y = 1,
           label = paste0(
             "<b><span style='font-size:26pt'>",
             dim(tradedFamilies)[1],
             "</b><br></span>",
             "<span style='font-size:8pt'>families traded</span>"),
    hjust = 1, vjust = -0.2,
    fill = NA, label.color = NA,
    label.padding = grid::unit(rep(0, 4), "pt"),
    lineheight = 0.9) +
  labs(x = "Number of species in family",
      y = "Number<br>of species<br>traded", colour = "Clade",
      size = "% traded", shape = "Clade") +
  scale_x_continuous(limits = c(1,10000), trans = "log",
                     breaks = c(1,10,100,1000, 10000)) +
  scale_size_continuous(limits = c(0, 100), breaks = seq(20, 100, by = 20)) +
  scale_colour_manual(values = colourway[c(6,2,8)]) +
  scale_y_log10() +
  theme_bw() +
  theme(
    legend.position = c(0.1, 0.65),
    legend.title = element_text(face = 2, hjust = 0),
    legend.background = element_blank(),
    axis.title.x = element_text(face = 2),
    axis.text.x = element_text(face = 1),
    axis.text.y = element_markdown(),
    axis.title.y = element_markdown(angle = 0, vjust = 1, hjust = 0, face = 2),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank()) +
  guides(
    # fill = guide_legend(label.position = "left"),
    color = guide_legend(), size = guide_legend())

ggsave("./Figures/Top traded familes scatter.png", dpi = 300,
       width = 160, height = 100, units = "mm")
ggsave("./Figures/Top traded familes scatter.pdf", dpi = 300,
       width = 160, height = 100, units = "mm")

# IUCN Status -------------------------------------------------------------

araData_trade %>% 
  group_by(clade) %>%
  count(status = redlist_Any, name = "any") %>% 
  filter(!status == "Not assessed")

araData_trade %>% 
  group_by(anyMatchTraded, clade) %>%
  count(status = redlist_Any, name = "any")

redlistPlotData <- araData_trade %>% 
  group_by(anyMatchTraded) %>%
  count(status = redlist_Any, name = "any") %>%
  left_join(araData_trade %>%
              group_by(anyMatchTraded) %>%
              count(status = redlist, name = "exact")) %>%
  # filter(!status == "Not assessed") %>%
  melt() %>%
  mutate(Traded = ifelse(anyMatchTraded, "Traded", "Not detected"),
         Traded = factor(Traded, levels = c("Traded", "Not detected"))
  ) %>% 
  ungroup()

redlistPlotData <- redlistPlotData %>% 
  mutate(
    colour = case_when(
      status == "Least Concern" ~ "#309706",
      status == "Near Threatened" ~ "#A8DB06",
      status == "Vulnerable" ~ "#F5D800",
      status == "Endangered" ~ "#DC7000",
      status == "Critically Endangered" ~ "#CB1500",
      status == "Extinct in the Wild" ~ "#701E08",â˜º
      status == "Extinct" ~ "#000000",
      status == "Data Deficient" ~ "#717171",
      status == "Not assessed" ~ "#BFBFBF"),
    statusSimp = case_when(
      status == "Least Concern" ~ "LC",
      status == "Near Threatened" ~ "NT",
      status == "Vulnerable" ~ "VU",
      status == "Endangered" ~ "EN",
      status == "Critically Endangered" ~ "CR",
      status == "Extinct in the Wild" ~ "EW",
      status == "Extinct" ~ "EX",
      status == "Data Deficient" ~ "DD",
      status == "Not assessed" ~ "NA"),
    statusCol = glue::glue("{status} <b style='color:{colour}'>({statusSimp})</b>"),
    statusCol = factor(statusCol, levels = c(
      "Not assessed <b style='color:#BFBFBF'>(NA)</b>",
      "Data Deficient <b style='color:#717171'>(DD)</b>",
      "Extinct <b style='color:#000000'>(EX)</b>",
      "Extinct in the Wild <b style='color:#701E08'>(EW)</b>",
      "Critically Endangered <b style='color:#CB1500'>(CR)</b>",
      "Endangered <b style='color:#DC7000'>(EN)</b>",
      "Vulnerable <b style='color:#F5D800'>(VU)</b>",
      "Near Threatened <b style='color:#A8DB06'>(NT)</b>",
      "Least Concern <b style='color:#309706'>(LC)</b>"
    ))) %>% 
  arrange(statusCol)

ggplot() +
  geom_col(
    data = redlistPlotData %>% 
      filter(!statusCol == "Not assessed <b style='color:#BFBFBF'>(NA)</b>",
             !statusCol == "Extinct <b style='color:#000000'>(EX)</b>",
             anyMatchTraded),
    aes(x = 1.45, y = value, fill = statusCol), width = 0.5) +
  geom_col(
    data = redlistPlotData %>%
      filter(!statusCol == "Not assessed <b style='color:#BFBFBF'>(NA)</b>",
             !statusCol == "Extinct <b style='color:#000000'>(EX)</b>",
             anyMatchTraded, variable == "any"),
    aes(x = 1, y = value, fill = statusCol), width = 0.2, alpha = 0.75) +
  geom_col(
    data = redlistPlotData %>%
      filter(!statusCol == "Not assessed <b style='color:#BFBFBF'>(NA)</b>",
             !statusCol == "Extinct <b style='color:#000000'>(EX)</b>",
             !anyMatchTraded),
    aes(x = 1.45, y = value, fill = statusCol), width = 0.5) +
  geom_col(
    data = redlistPlotData %>%
      filter(!statusCol == "Not assessed <b style='color:#BFBFBF'>(NA)</b>",
             !statusCol == "Extinct <b style='color:#000000'>(EX)</b>",
             !anyMatchTraded, variable == "any"),
    aes(x = 1, y = value, fill = statusCol), width = 0.2, alpha = 0.75) +
  geom_text(data = redlistPlotData %>% 
              group_by(Traded, variable) %>% 
              filter(!status == "Not assessed") %>% 
              summarise(tot = sum(value)),
            aes(x = c(1.075,0.975,1.075,0.975), y = c(Inf, Inf, Inf, Inf),
                label = paste0(str_to_sentence(variable), " = ", tot)),
            hjust = 1, vjust = 1, fontface = 2, size = 2) +
  scale_x_continuous(labels = c("Exact match", "Any match"), breaks = c(1, 1.45)) +
  scale_fill_manual(values = c(
    "#717171",
    "#CB1500",
    "#DC7000",
    "#F5D800",
    "#A8DB06",
    "#309706"
  )) +
  labs(y = "Count of species", x = "", fill = "IUCN Redlist") +
  coord_flip() +
  facet_wrap(.~Traded, ncol = 1, scales = "free") +
  theme_bw() +
  theme(
    legend.text = element_markdown(),
    legend.title = element_text(face = 2),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = 4, hjust = 0.035),
    panel.grid.major.y = element_blank(),
    axis.line.x = element_line(),
    axis.title.x = element_text(face = 2),
    axis.text.y = element_text(face = 3, size = 8),
    axis.ticks.y = element_blank(),
    panel.border = element_blank()
  )

ggsave("./Figures/IUCN status breakdown.png", dpi = 300,
       width = 180, height = 70, units = "mm")
ggsave("./Figures/IUCN status breakdown.pdf",
       width = 180, height = 70, units = "mm")

# Language summary --------------------------------------------------------

speciesWebRaw <- read.csv(file = "./Data/Snapshot Online Species Data Raw.csv")

targetSites <- read.csv(file = "./Data/Target Websites Censored.csv")

siteLangData <- speciesWebRaw %>% 
  left_join(targetSites %>% 
              select(lang, webID))

siteCountData <- siteLangData %>% 
  group_by(lang) %>% 
  summarise(nsites = length(unique(webID))) %>% 
  mutate(
    language = 
      case_when(
        lang == "CZE" ~ "Czech",
        lang == "ENG" ~ "English",
        lang == "FRA" ~ "French",
        lang == "GER" ~ "German",
        lang == "JPN" ~ "Japanese",
        lang == "POL" ~ "Polish",
        lang == "RUS" ~ "Russian",
        lang == "SPA" ~ "Spanish",
        lang == "SWE" ~ "Swedish"
      ))

uniquePerLang <- siteLangData %>% 
  group_by(lang, sp) %>% 
  filter(!duplicated(sp)) %>% ## make sure each species appears once in each lang
  group_by(sp) %>% 
  mutate(appearances = length(sp)) %>% ## then count how many that times that species is in the dataset, 
  select(lang, sp, appearances) %>%  ## > 1 means it was detected in more than one language
  filter(appearances == 1) %>% 
  group_by(lang) %>% 
  count(name = "nunispp")

siteLangData %>% 
  group_by(lang, sp) %>% 
  slice(n = 1) %>%
  ungroup() %>% 
  count(lang) %>% 
  left_join(uniquePerLang) %>% 
  mutate(
    language = 
      case_when(
        lang == "CZE" ~ "Czech",
        lang == "ENG" ~ "English",
        lang == "FRA" ~ "French",
        lang == "GER" ~ "German",
        lang == "JPN" ~ "Japanese",
        lang == "POL" ~ "Polish",
        lang == "RUS" ~ "Russian",
        lang == "SPA" ~ "Spanish",
        lang == "SWE" ~ "Swedish"
      ),
    nunispp = ifelse(is.na(nunispp), 0, nunispp)
  ) %>% 
  arrange(n) %>% 
  ggplot() +
  geom_col(aes(x = reorder(language, -n), y = n), width = 0.45,
           fill = colourway[3]) +
  geom_text(aes(x = reorder(language, -n), y = n+5, label = n),
            vjust = 0, fontface = 2,
            col = colourway[3]) +
  geom_col(aes(x = reorder(language, -n), y = nunispp), width = 0.45,
           fill = colourway[7]) +
  geom_text(aes(x = reorder(language, -n), y = nunispp+5, label = nunispp),
            vjust = 0, fontface = 2, hjust = 0,
            colour = colourway[7], position = position_nudge(x = 0.3)) +
  geom_segment(aes(x = seq(9,1,-1)-0.225, xend = seq(9,1,-1)+0.4,
                   y = nunispp, yend = nunispp),
               colour = colourway[7]) +
  annotate("segment", x = 9, xend = -0.1, y = -25, yend = -25) +
  annotate("text", x = -0.11, y = -25, label = "Number\\nof sites",
           lineheight = 0.85, fontface = 4, hjust = 1) +
  geom_point(data = siteCountData,
             aes(x = language, y = -25), colour = "white", pch = 16, size = 10.5) +
  geom_point(data = siteCountData,
             aes(x = language, y = -25), colour = "black", pch = 16, size = 7.5) +
  geom_point(data = siteCountData,
             aes(x = language, y = -25), colour = "white", pch = 16, size = 6.5) +
  geom_text(data = siteCountData,
            aes(x = language, y = -25, label = nsites)) +
  coord_cartesian(clip = "off", xlim = c(1, 9)) +
  labs(x = "Language", y = "Number of\\nspecies") +
  theme_bw() +
  theme(
    axis.title.x = element_text(face = 2),
    axis.text.x = element_text(face = 3),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(angle = 0, vjust = 1, hjust = 0, face = 2),
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.grid.major.x = element_blank()
  )

ggsave("./Figures/Language summary.png", dpi = 300,
       width = 160, height = 120, units = "mm")
ggsave("./Figures/Language summary.pdf",
       width = 160, height = 120, units = "mm")
