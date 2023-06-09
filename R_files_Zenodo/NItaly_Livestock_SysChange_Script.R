###########################################################
# 
# STATISTICAL ANALYSIS & ORIGINAL FIGURES
# For the paper:
# Systems change: Investigating climatic and environmental impacts 
# on livestock production in lowland Italy between 
# the Bronze Age and Late Antiquity (c. 1700 BC - AD 700) 
#
# By Angela Trentacoste, Ariadna Nieto-Espinet, Silvia Guimaraes Chiarelli, 
# and Silvia Valenzuela-Lamas
#
#
# Script by A. Trentacoste 
#
###########################################################

# If using or adapting please cite the above paper

### Set up ################################################

# Set working directory - Remember to set the correct working directory on your computer

# Load libraries
library(tidyverse)
library(FactoMineR)
library(zoolog)
library(rstatix)
library(ggpubr)
library(dendextend)
library(openxlsx)

# Files
input_file_NISP <- "Supp01_Site_NISP_Landscape_Data.csv" 
input_file_biometry <- "NItaly_Livestock_Metric_Data.csv"

# Period order
per.order <- c("BA.Middle", "BA.Middle-Late", "BA.Late", "BA.Late-IA.Early",
               "IA.Early", "IA.Mid", "IA.Late",
               "Roman.Early", "Roman.Mid", "Roman.Mid-Late", "Roman.Late",
               "Roman.All")

### NISP Percentages ##########################################################
NItaly <- read.csv(input_file_NISP)
# Simplify regional grouping by assigning 'Appennine' sites to South study area
NItaly$Region[NItaly$Region == "Appennine"] <- "South"
# Calculate NISP%
NItaly <- NItaly %>%
  mutate(
    Per.Cattle = Cattle/NISP3dom,
    Per.SG = SheepGoat/NISP3dom,
    Per.Pig = Pig/NISP3dom)

### Correspondence Analysis ###################################################
# Correspondence Analysis (CA) of NISP counts (cattle, pig and sheep/goat)

# CA script in this section adapted from the zoo_ca() function 
# in the zoowork package (Huet et al. 2022)
# Huet, T., Nieto-Espinet, A., Trentacoste, A., Guimarães, S., 
# & Valenzuela-Lamas, S. (2022). zoowork (Version 1.1.0.0) [Computer software]. 
# https://doi.org/10.5281/zenodo.6850736

bySite <- NItaly %>%
  select(Region, Region, Period, Chron, Per.Cattle, 
         Per.SG, Per.Pig, Latitude, Longitude) %>%
  as.data.frame()
Site.ca <- NItaly %>%
  select(Per.Cattle, Per.SG, Per.Pig) %>%
  as.data.frame()
ca <- CA(Site.ca, graph = TRUE)
# CA plot
coords_ind_ca <- as.data.frame(ca$row$coord)
coords_var_ca <- as.data.frame(ca$col$coord)
coords_ca <- rbind(coords_ind_ca,coords_var_ca)
colnames(coords_ca)[1] <- 'CA1'
colnames(coords_ca)[2] <- 'CA2'
main.points <- coords_ca[103:105, ]
coords_ca <- coords_ca[-c(103:105), ]
dataset.p <- cbind(bySite, coords_ca)
main.points$Taxon <- NA
main.points$Taxon <- c("cattle", "sheep/goat", "pig")
dataset.p$Chron <- factor(dataset.p$Chron, levels = per.order)
# Exclude non-specific Roman sites
dataset.p <- subset(dataset.p, Chron != "Roman.All") 

##### Figure 3 - Correspondance Analysis plot
ggplot(dataset.p, aes(CA1, CA2)) +
  geom_point(data = ~select(., -Region),
             aes(CA1, CA2), size = 1.3, colour = "grey90") +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.25 ) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.25 ) +
  geom_point(data = dataset.p, aes(CA1, CA2, colour = Chron), size = 1.3) +
  geom_point(data = main.points, aes(CA1, CA2),
             size = 1.3, shape = 17, colour = "black") +
  geom_text(data = main.points, aes(CA1, CA2, label = Taxon),
            hjust = 0.0, vjust = -1, size = 2)  +
  scale_colour_manual(values =  hcl.colors(12, palette = "Zissou 1"),
                      name = "Period") +
  scale_fill_manual(values =  hcl.colors(12, palette = "Zissou 1"),
                    name = "Period") +
  theme_bw(base_size = 7) +
  theme(legend.key.size = unit(0.3, "cm")) +
  ylim(-.8, .8) +
  xlim(-1.3, 0.9) +
  facet_wrap(factor(Region, 
                    levels = c("Veneto", "Friuli", "North", "South"))~.) +
  xlab("Dim 1 64.44%") + ylab("Dim 2 35.56%")
ggsave("Fig03_CA.tiff", units = "in", height = 4.2, width = 7.48, dpi = 300)

### Hierarchical clustering analysis (HCA) ####################################
Site.hca <- NItaly %>%
  select(Site.name, HCA_Order, Region, Region, Chron, 
         Per.Cattle, Per.SG, Per.Pig) %>%
  as.data.frame() %>%
  column_to_rownames(var="HCA_Order")
Site.hca <- Site.hca %>%
  mutate( HCAcolour = case_when(
    Region =="South" ~ "#F21A00", 
    Region == "North" ~ "#E1AF00", 
    Region == "Friuli" ~ "#3B9AB2",
    Region == "Veneto" ~ "#33CC99"))

# Dendrograms
px <- Site.hca %>% subset(Chron == "BA.Middle")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p1 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "BA.Middle-Late")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p2 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "BA.Late")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p3 <- ggplot(dend) + ylim(0,1.5)  + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "IA.Early")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p4 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "IA.Mid")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p5 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "IA.Late")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p6 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "Roman.Early")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p7 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "Roman.Mid")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p8 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "Roman.Mid-Late")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>% 
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p9 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
px <- Site.hca %>% subset(Chron == "Roman.Late")
dend <- hclust(dist(px[,4:6]), method = "complete") %>% as.dendrogram
dend <- dend %>% set("labels_color", px$HCAcolour, order_value = TRUE) %>%
  set("labels_cex", 0.35) %>% set("branches_lwd", 0.4)
p10 <- ggplot(dend) + ylim(0,1.5) + 
  theme(plot.margin = margin(r = 0, l = 0))
##### Figure 4b - HCA of NISP profiles
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 1)
ggsave("Fig04b_HCA.pdf", height = 60, width = 195, units = "mm")
# Shapes added and label position adjusted in Adobe Illustrator

### Small vs large livestock ##################################################
# Relative proportions of small livestock (sheep/goat/pigs) versus cattle
NR <- NItaly
NR <- subset(NR, Chron != "Roman.All") #Exclude non-specific Roman sites
# Simplify regional grouping to north and south
NR$Region[NR$Region == "Appennine"] <- "South"
NR$Region[NR$Region == "Veneto"] <- "North"
NR$Region[NR$Region == "Friuli"] <- "North"
NISP.percentages <- NR %>% mutate(
  Cattle.per = Cattle/NISP3dom,
  Small.per = (SheepGoat + Pig)/NISP3dom) 
Percentages.long <- NISP.percentages %>% 
  pivot_longer(Cattle.per:Small.per, names_to = "Taxon", 
               values_to = "Percentage")
Percentage.Chron <- Percentages.long %>% 
  group_by(Chron, Region, Taxon) %>%
  summarise(n = n(),
            TotNISP=sum(NISP3dom),
            mean.percentage = mean(Percentage),
            sd = sd(Percentage)) 
Percentage.Chron$Chron <- factor(Percentage.Chron$Chron, levels=per.order)
##### Figure 4a
##### Relative proportions of small livestock (sheep/goat/pigs) versus cattle
p <- ggplot(subset(Percentage.Chron, Chron != "BA.Late-IA.Early"), 
            aes(Chron, mean.percentage, group = Taxon)) +
  geom_line(aes(colour = Taxon), size=0.5, show.legend = FALSE) +
  geom_jitter(data = subset(Percentages.long, Chron != "BA.Late-IA.Early"),
              aes(x=Chron, y= Percentage, colour = Taxon, shape=RomanType),
              size = 1.2, 
              position = position_jitterdodge(jitter.width = 0.15, 
                                              dodge.width = 0.5, seed = NA) ) +
  scale_colour_manual(values =  hcl.colors(4, palette = "Zissou 1"), 
                      name = "Taxon") +
  scale_shape_manual(values = c(16, 1, 15, 17, 15, 17), name = "Type", 
                     breaks = c("other","minor", "IA - major1", "IA - major2", 
                                "Rom - urban1", "Rom - urban2"))+ 
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                     position = "right") +
  facet_grid(Region~.) +
  theme_bw(base_size = 7) +
  labs(y = "%NISP", x = "Period")  +
  theme(legend.position=c(.91,.25),
        legend.background = element_rect(fill = "transparent", 
                                         colour = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.key.height = unit(0.35, "cm"))
p
ggsave("Fig04a_LargevSmall_Livestock.pdf", units = "mm", width = 210, height = 70 )

### Biometry ##################################################################
Bio <- read.csv(input_file_biometry)
# Simplify period grouping
Bio$Chron[Bio$Chron == "BA.L-EIA"] <- "IA.Early" # Group EIA & Late BA Sites
# Exclude all but one specimen from single individuals
Bio <- subset(Bio, is.na(IndiDuplicate))
# Change taxon codes for pigs ('S' and 'SS') to align with zoolog Thesaurus
Bio <- Bio %>%
  mutate(Taxon=replace(Taxon, Taxon=="SS", "wild boar"),
         Taxon=replace(Taxon, Taxon=="S", "sus"))
# Define joined categories and calculate LSI values
# sheep, goat, and sheep/goat grouped with sheep reference
# pig and wild boar grouped with wild boar reference
cats.joined <- list(OVA = c("OVA", "CAH", "O"), 
                    wss = c("sus", "wild boar"))
Bio.log <- LogRatios(Bio, identifiers = c("TAX", "EL"), 
                     ref = reference$Combi, joinCategories = cats.joined)
# Remove cases without LSI values and condense
Bio.log.pruned <- RemoveNACases(Bio.log)
Bio.con <- CondenseLogs(Bio.log.pruned)
# Analysis will focus on domestic pigs
# Visualise pig and wildboard measurements
ggplot(subset(Bio.con, Taxon == "sus"|Taxon == "wild boar")) + 
  geom_jitter(aes(x = Chron, y = Width, colour = Taxon)) + 
  ylim(-0.3, 0.25) + geom_hline(yintercept = 0, color = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggplot(subset(Bio.con, Taxon == "sus"|Taxon == "wild boar")) + 
  geom_jitter(aes(x = Chron, y = Length, colour = Taxon)) + ylim(-0.2, 0.1) + 
  geom_hline(yintercept = 0, color = "red") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# wild boar and high sus values probably from wild boar removed
Bio.con <- Bio.con %>% 
  mutate(
    Width = ifelse(Taxon %in% c("sus", "wild boar") & Width > 0, NA, Width), 
    Length = ifelse(Taxon %in% c("sus", "wild boar") & Length > 0, NA, Length))
Bio.con <- subset(Bio.con, Taxon != "wild boar")
Bio.con <- Bio.con %>% 
  filter(!is.na(Width) | !is.na(Length))
##### Summary statistics - Supplement 2 Biometry Summary Statistics
Bio.stats <- Bio.con %>%
  group_by(Chron, TaxonGroup) %>%
  summarise(
    n=n(),
    n_sites=n_distinct(SiteID),
    W_mean = mean(Width,  na.rm = TRUE),
    W_sd = sd(Width,  na.rm = TRUE),
    W_min = min(Width,  na.rm = TRUE),
    W_max = max(Width,  na.rm = TRUE),
    L_mean = mean(Length,  na.rm = TRUE),
    L_sd = sd(Length,  na.rm = TRUE),
    L_min = min(Length,  na.rm = TRUE),
    L_max = max(Length,  na.rm = TRUE))
write.csv(Bio.stats,"Supplement02_SummaryStats.csv")

### Biometry plot #############################################################
Bio.long <- Bio.stats %>% 
  pivot_longer(
    cols= W_mean:L_sd, 
    names_to = c("dimension", "stat"),
    names_sep = "\\\\_",
    values_to = "value.new") 
Bio.long <- Bio.long %>% pivot_wider(Chron:value.new, names_from = stat, values_from = value.new)
Bio.long$Chron <- factor(Bio.long$Chron, levels=per.order)
##### Figure 4c - Mean livestock LSI values
p <- ggplot(subset(Bio.long, TaxonGroup %in% c("Cattle", "Sheep/goat", "Pig")), 
                aes(x=Chron, y=mean)) +
  geom_line(aes(group=dimension, colour = dimension), size = .81) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, colour = dimension), 
                width=.05, size = 0.15, position=position_dodge(0.05)) +
  facet_grid(TaxonGroup~.) +
  scale_colour_manual(values =  hcl.colors(2, palette = "Zissou 1"), 
                      name = "Dimension") +
  coord_cartesian(ylim = c(-0.175, 0.175)) +
  scale_y_continuous(position = "right") + 
  theme_bw(base_size = 7) +
  theme(legend.position=c(0.9, 0.55),
        legend.background = element_rect(fill = "transparent", 
                                         colour = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.key.height = unit(0.35, "cm"))
p
ggsave("Fig04c_Biometry_LinePlot.pdf", units = "mm", width = 210, height = 78)
# p-value stars added in Adobe Illustrator

### Biometry MannU tests ######################################################
# Cattle
stats.Cattle.W <- Bio.con %>%        
  subset(Taxon == "B") %>%          
  drop_na(Width) %>%             
  wilcox_test(Width ~ Chron) %>% 
  add_significance() %>% 
  mutate(Taxon = "Cattle") %>% 
  relocate(Taxon, .before = 1)
stats.Cattle.L <- Bio.con %>%
  subset(Taxon == "B" ) %>% 
  drop_na(Length) %>% 
  wilcox_test(Length ~ Chron) %>%
  add_significance() %>% 
  mutate(Taxon = "Cattle") %>% 
  relocate(Taxon, .before = 1)
# Sheep/goat
stats.SG.W <- Bio.con %>%
  subset(Taxon %in% c("O", "OVA", "CAH")) %>%
  drop_na(Width) %>%  
  wilcox_test(Width ~ Chron) %>%
  add_significance()%>% 
  mutate(Taxon = "SheepGoat") %>% 
  relocate(Taxon, .before = 1)
stats.SG.L <- Bio.con %>%
  subset(Taxon %in% c("O", "OVA", "CAH")) %>%
  drop_na(Length) %>% 
  wilcox_test(Length ~ Chron) %>%
  add_significance() %>%
  mutate(Taxon = "SheepGoat") %>% 
  relocate(Taxon, .before = 1)
# Pig
stats.Pig.W <- Bio.con %>%
  subset(Taxon %in% c("sus",  "wild boar")) %>%
  drop_na(Width) %>% 
  wilcox_test(Width ~ Chron) %>%
  add_significance() %>%
  mutate(Taxon = "Pig") %>% 
  relocate(Taxon, .before = 1)
stats.Pig.L <- Bio.con %>%
  subset(Taxon %in% c("sus",  "wild boar")) %>%
  drop_na(Length) %>% 
  wilcox_test(Length ~ Chron) %>%
  add_significance() %>% 
  mutate(Taxon = "Pig") %>% 
  relocate(Taxon, .before = 1)

##### Supplement 3 - MannU Results
MannUResults <- bind_rows(stats.Cattle.L, stats.Cattle.W, stats.SG.L, 
                          stats.SG.W, stats.Pig.L, stats.Pig.W)
write.csv(MannUResults, "Supplement03_MannU_Results.csv")

### Landscape variables - Correlations ########################################
NItaly <- NItaly %>% subset(Chron != "Roman.All") # General Roman sites excluded
# Correlations between 5 km landscape variables and cattle percentages
Corr.Clay <- NItaly %>%
  group_by(ChronGeneral) %>% 
  summarize(name = c("Clay"),
            n = n(),
            cor_coef = cor.test(Clay_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(Clay_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(Clay_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.AWC <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name="AWC",
            n = n(),
            cor_coef = cor.test(AWC_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(AWC_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(AWC_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.CF <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "CF", 
            n = n(),
            cor_coef = cor.test(CF_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(CF_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(CF_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.Sand <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Sand",
            n = n(),
            cor_coef = cor.test(Sand_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(Sand_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(Sand_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.Den <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Density",
            n = n(),
            cor_coef = cor.test(Den_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(Den_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(Den_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.Silt <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Silt",
            n = n(),
            cor_coef = cor.test(Silt_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(Silt_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(Silt_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.SMRT5 <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "SRTM 5km mean",
            n = n(),
            cor_coef = cor.test(SRTM5km_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(SRTM5km_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(SRTM5km_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.Precip <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Precipitation",
            n = n(),
            cor_coef = cor.test(Precip_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(Precip_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(Precip_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.SImean <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Solar mean",
            n = n(),
            cor_coef = cor.test(SI_mean, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(SI_mean, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(SI_mean, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr.SImax <- NItaly %>%
  group_by(ChronGeneral) %>%
  summarize(name = "Solar max",
            n = n(),
            cor_coef = cor.test(SI_max, Per.Cattle, 
                                method = "spearman")$estimate,
            p_val = cor.test(SI_max, Per.Cattle, 
                             method = "spearman")$p.value,
            statistic = cor.test(SI_max, Per.Cattle, 
                                 method = "spearman")$statistic)
Corr <- list(Corr.Clay, Corr.AWC, Corr.CF, Corr.Sand, Corr.Den, Corr.Silt, 
             Corr.SMRT5, Corr.Precip, Corr.SImean, Corr.SImax)

##### Correlations between landscape variables and cattle percentages
write.xlsx(Corr, file = "Supplement04_Correlations.xlsx")

### Landscape variable plot - Middle-Recent Bronze Age
LandData <- NItaly %>% pivot_longer(AWC_mean:SI_max, names_to = "LandChars", 
                         values_to="Value")
fac <- c("AWC_mean", "Den_mean", "Clay_mean", "CF_mean", "Sand_mean", 
         "Silt_mean","SRTM5km_mean", "Precip_mean","SI_mean", "SI_max")
LandData <- LandData %>% subset(LandChars %in% fac)
LandData$LandChars <- factor(LandData$LandChars, levels = fac)
##### Figure 5 - Relative percentages of cattle from Mid-Recent Bronze Age sites 
##### compared to landscape variables
ggscatter(subset(LandData, ChronGeneral == "BA Mid-Recent"), 
          x = "Value", y = "Per.Cattle", 
          corr.method= c("spearman"), size = 0.5, 
          add = "reg.line", conf.int = TRUE) +
  stat_cor(label.y = .85, method = "spearman",  size = 2) +
  yscale("percent") +
  facet_wrap(LandChars~ChronGeneral, scales="free", ncol=4) +
  coord_cartesian(ylim = c(0,1)) +
  theme_bw(base_size = 7) +
  ylab("% Cattle")
ggsave("Fig05_BA_Cattle_Landscape_Comp.pdf", height = 4, width = 7.48, units = "in")    

####################################################### End
