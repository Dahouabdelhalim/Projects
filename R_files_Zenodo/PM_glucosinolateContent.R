library(tidyverse)
library(viridisLite)
library(ggpubr)


PM_plotTheme <- theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(color = "black"), 
        strip.text.x  = element_blank())

#### Filter paper weight ####
paper_g <- c(0.0608, 0.0595, 0.0590, 0.0603, 0.0612, 0.0634, 0.0602, 0.0600, 0.0606, 0.0616, 0.0624, 0.0628, 0.0636, 0.0585, 0.0636)
summary(paper_g)
mean_g <- mean(paper_g)
se_g <- 1.96*sd(paper_g)/sqrt(15)


#### GSL Class Master ####
GSL_class_master <-read.csv("GSL_classMaster.csv")

#### MeOH extracts ####

## all progoitrin and epiprogoitrin in the samples are from the standard
meoh_epipro <- read.csv("dGSLs_6-2015_MeOHextracts.csv") %>% 
  filter(NameGSL == "Epiprogoitrin") %>% 
  mutate(std_integ_area_DAD_tot = if_else(is.na(std_integ_area_DAD + DAD_integ_area), std_integ_area_DAD, std_integ_area_DAD + DAD_integ_area), 
         std_integ_area_CAD_tot = if_else(is.na(std_integ_area_CAD + CAD_integ_area), std_integ_area_CAD, std_integ_area_CAD + CAD_integ_area)) %>% 
  dplyr::select(sample_id, std_integ_area_DAD_tot, std_integ_area_CAD_tot)

meoh <- read.csv("dGSLs_6-2015_MeOHextracts.csv") %>% 
  # dplyr::select(-c(group, batch, vial_id)) %>% 
  # left_join(meoh_epipro) %>% 
  mutate(std_integ_area_DAD_tot = if_else(is.na(std_integ_area_DAD), std_integ_area_DAD, std_integ_area_DAD), 
         std_integ_area_CAD_tot = if_else(is.na(std_integ_area_CAD), std_integ_area_CAD, std_integ_area_CAD)) %>% 
  filter(NameGSL != "Epiprogoitrin") %>% 
  left_join(GSL_class_master) %>% 
 filter(ClassGSL != "NA") %>% 
  mutate(desulfoMW = as.factor(desulfoMW),
         ClassGSL = factor(ClassGSL, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")),
         FaheyNum = factor(FaheyNum, levels = c("66", "69", "56", "61", "62", "12", "101", "107", "30" , "31", "11", "40 (R/S)", "105", "43")),
         type = "MeOH extract",
         y_axis ="umol/g filter paper",
         DAD_umol = (5e-02 * DAD_integ_area)/std_integ_area_DAD_tot, # umol PRO * ratio sample:pro = umol 
         CAD_umol = (5e-02 * CAD_integ_area)/std_integ_area_CAD_tot,
         DAD_umol_per_g = (80*DAD_umol/vol)/mean_g, # (sample nanomoles/sample vol * vol applied to filter paper/filter paper g ) = sample nanomoles/ g filter paper
         CAD_umol_per_g = (80*CAD_umol/vol)/mean_g) %>% 
  group_by(sample_id, ClassGSL, type, y_axis, species, vol) %>% 
  summarise(ClassTotal = sum(CAD_umol, na.rm = T)) %>%
  mutate (total_umol_per_g = (80*ClassTotal/vol)/mean_g)


meoh_plot <- ggplot(meoh, aes(x = ClassGSL, total_umol_per_g)) +
  geom_boxplot(aes(color = ClassGSL)) + 
  geom_jitter(shape = 21, fill = NA, width = 0.05) + 
  labs(y = NULL, x = "Glucosinolate") +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  facet_grid(fct_rev(y_axis) ~ species, scales = "free", space = "free_x", switch = "y") +
  PM_plotTheme +  
  theme(strip.placement = "outside", strip.text.x = element_text(hjust = 0, face = "italic"))



#### Leaf leachates ####
leach <- read.csv("dGSLs_6-2015_leafsamples_fromMaster.csv") %>% 
  left_join(GSL_class_master) %>% 
  filter(ClassGSL != "NA", sample_id != "Caco.15-18", plant_or_extract != "D. incana") %>% 
  ### sinigrin found in cardamine sample likely contamination, sample removed.
  ### D. incana leachates removed
  rename(species = plant_or_extract) %>% 
  mutate(desulfoMW = as.factor(desulfoMW),
         ClassGSL = factor(ClassGSL, levels = c("A", "B", "C", "D", "E", "F", "G", "H", "I")),
         FaheyNum = factor(FaheyNum, levels = c("66", "69", "56", "61", "62", "12", "101", "107", "30" , "31", "11", "40 (R/S)", "105", "43")),
         type = "Leachate",
         y_axis = "umol/g wet leaf",
         DAD_umol = (5e-02 * DAD_integ_area)/std_integ_area_DAD, # umol PRO * ratio sample:pro = umol 
         CAD_umol = (5e-02 * CAD_integ_area)/std_integ_area_CAD,
         DAD_umol_per_g = DAD_umol/wet_leaf_g ,
         CAD_umol_per_g = CAD_umol/wet_leaf_g) %>% 
  group_by(sample_id, ClassGSL, type, y_axis, wet_leaf_g, species) %>% 
  summarise(ClassTotal = sum(CAD_umol)) %>%
  mutate (total_umol_per_g = ClassTotal/wet_leaf_g) %>% filter(ClassGSL != "D" | species == "T. arvense")

leach_plot <- ggplot(leach, aes(x = ClassGSL, total_umol_per_g)) +
  geom_boxplot(aes(color = ClassGSL)) + 
  geom_jitter(shape = 21, fill = NA, width = 0.05) + 
  labs(y = NULL, x = "Glucosinolate") +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  facet_grid(fct_rev(y_axis) ~ species, scales = "free_x", space = "free_x", switch = "y") +
  PM_plotTheme +  
  theme(strip.placement = "outside", strip.text.x = element_text(hjust = 0, face = "italic"))

ggarrange(meoh_plot +theme(axis.title.x = element_blank()), leach_plot + theme(strip.text.x = element_blank()),
          ncol = 1, nrow = 2, labels = c("a", "b"))
ggsave("dGSL_2015_extracts_v_leachates.pdf", width = 9, height = 6)
