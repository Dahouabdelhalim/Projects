library(phylowood)
library(tidyverse)
library(viridis)
library(readxl)
library(writexl)
library(deeptime)
library(wesanderson)

archp <- phylowood_isl %>% 
  dplyr::select(geo_archipelago, tax_genus, archp_taxclade = tax_clade)

age_arch_match <- read_xlsx("input/ages_archipelago_link.xlsx")

l_age <- phylowood_ages %>% 
  filter(trt_island_lineage) %>% 
  group_by(tax_genus, tax_clade) %>% 
  mutate(date = parse_number(ref_reference)) %>% 
  top_n(n=1) %>% 
  ungroup() %>% 
  dplyr::select(tax_genus, tax_clade, ref_reference, contains("stemage")) %>% 
  filter(!is.na(trt_mean_stemage)) %>% 
  left_join(age_arch_match, by = c("tax_genus", "tax_clade")) %>% 
  filter(archipelago != "not_iw") %>% 
  mutate(tax_genus = recode(tax_genus, 
                            `Convolvulus lineage I` = "Convolvulus",
                            `Convolvulus lineage II` = "Convolvulus",
                            `Argyroxiphium,Dubautia (silversword clade)` = "Silversword clade")) %>% 
  mutate(clade_ID = paste(tax_genus, "_", 1:nrow(.), sep = "")) %>% 
  filter(archipelago != "continental") %>% 
  filter(archipelago != "NA") %>% 
  mutate(trt_mean_stemage = parse_number(trt_mean_stemage)) %>% 
  mutate(trt_max_stemage = parse_number(trt_max_stemage))%>% 
  mutate(trt_min_stemage = parse_number(trt_min_stemage))

save(l_age, file = "output/island_age_for_manuscript.rda")

# a list of cases with multiple mentions of genera for quality control
list <- table(l_age$tax_genus) 
list <- names(list[list> 1])

list <- table(l_age$tax_clade) 
list <- names(list[list> 1])

probs <- l_age %>% filter(tax_clade %in% list[1:2])

write_xlsx(probs, "output/clades_with_multiple_ages.xlsx")

# when ato archipelagos, add island ages
clade_labs <- data.frame(clade_ID = l_age$clade_ID) %>% 
  mutate(clade_name = str_split_fixed(clade_ID, n = 2, pattern = "_")[,1]) %>%
  mutate(clade_name = gsub(" \\\\(Caribbean\\\\)", "", clade_name)) %>% 
  mutate(clade_name = str_split_fixed(clade_name, n = 2, pattern = " ")[,1])

clade_labs <- split(clade_labs, f = clade_labs$clade_name)
clade_labs <-lapply(clade_labs, function(k){k <- data.frame(k, number = seq(from = 1, to = nrow(k)))})
clade_labs <-  bind_rows(clade_labs) %>% 
  mutate(number = as.roman(number)) %>% 
  mutate(clade_name = paste(clade_name, number, sep = " "))

l_age <- left_join(l_age, clade_labs %>% dplyr::select(clade_ID, clade_name))
l_age$clade_name <- factor(l_age$clade_name, levels = l_age$clade_name[rev(order(l_age$archipelago, l_age$trt_mean_stemage))])

# Age of all islands in the archipelago
load(file = "output/island_per_archipelago_age.rda")
isl_age <- isl_age %>% 
  filter(!is.na(age_Ma)) %>% 
  filter(Arch_name %in% l_age$archipelago) %>% 
  rename(archipelago = Arch_name) %>% 
  filter(age_Ma <= 25)

# age of the first fossils for all different mammal species, this is a separate data file 
#Roberto has compiled for the archipelagos with dated shifts
herbs <-  read_delim("input/herbivores_archipelagos_Roberto.csv", delim = ";") %>% 
  dplyr::select(archipelago = Arch_name, FAD)


arch_lab <- data.frame(archipelago = sort(unique(l_age$archipelago)),
                       x = 24.5,
                       y = Inf)

pal <- wes_palette("FantasticFox1", 16, type = "continuous")
pal2 <- c(wes_palette("IsleofDogs1", 2)[2], wes_palette("Darjeeling1", 1))

pl <- ggplot(data = l_age)+
  geom_segment(aes(x = trt_min_stemage,
                   xend = trt_max_stemage, 
                   y = clade_name,
                   yend = clade_name))+#,
                  #color = archipelago))+
  geom_point(aes(x = trt_mean_stemage, 
                 y = clade_name))+#,
                # color = archipelago))+
  geom_point(data = isl_age,
             aes(x = age_Ma, 
                 y = -Inf),
             size = 4, shape = 3, color = pal2[1])+
  geom_point(data = herbs,
             aes(x = FAD, 
                 y = Inf),
             size = 4, shape = 3, color = pal2[2])+
  annotate("rect", 
           xmin = 16.8,
           xmax = 14.7,
           ymin = -Inf,
           ymax = Inf,
           fill = "grey",
           alpha = 0.3)+
  annotate("rect", 
           xmin = 3.264,
           xmax = 3.025,
           ymin = -Inf,
           ymax = Inf,
           fill = "grey",
           alpha = 0.3)+
  annotate("rect", 
           xmin = 2.8,
           xmax = 2.6,
           ymin = -Inf,
           ymax = Inf,
           fill = "grey",
           alpha = 0.3)+
  annotate("rect", 
           xmin = 1.25,
           xmax = 0.7,
           ymin = -Inf,
           ymax = Inf,
           fill = "grey",
           alpha = 0.3)+
  geom_text(data = arch_lab, aes(x = x, 
                                 y = y, 
                                 label = archipelago, vjust = 1.1, hjust = 0), size = 2.5)+
                                 # color = archipelago)
  #scale_color_manual(values = pal)+
  facet_grid(archipelago~., scales = "free_y", space = "free_y") +
  xlab("Age [Ma]")+
  scale_x_reverse(expand = c(0.01,0))+
  theme_bw()+
  theme(
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(0.2, "lines"), 
   # panel.border = element_rect(color = "grey40", size = 0.1),
    strip.background = element_blank(),
    strip.text = element_blank(),
    strip.placement = "outside",
    legend.position = "none")

out <- gggeo_scale(pl, dat = "epochs", abbrv = FALSE, height = unit(1, "line"))
ggsave(plot = out, "output/figure_island_and_clade_ages.pdf", height = 9, width = 7)
