#03 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publication


### Paired with
# "01_vitE_10d_37d_stacked.csv"


#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting




setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #sets WD to folder this R file is saved in

vitE.stacked <- read.csv("01_vitE_10d_37d_stacked.csv", header = TRUE) #identical to "00_vitE_10d_37d_stacked.csv"




vitE.means <- vitE.stacked %>%
  group_by(Phase, Treatment) %>%
  summarize(Avg.E = mean(as.numeric(E)), 
            SD.E = sd(as.numeric(E))) %>%
  mutate(E.min = Avg.E - SD.E) %>%
  mutate(E.max = Avg.E + SD.E) 




treatment.colors <- c("light grey", "grey55", "black")

x11()
pE <- ggplot(vitE.means,  aes(x = as.factor(Phase), y = Avg.E, color = Treatment)) +
  geom_point(size = 3, position = position_dodge(width = 0.40)) +
  geom_errorbar(aes(x = as.factor(Phase), ymax = (E.max), ymin = (E.min)), width = 0.5, size = 1, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("Control", "PUFA", "PUFA + E")) +
  labs(x = "Phase", y = "\\u03B1-tocopherol (\\u03BCg/g) \\u00B1 SD") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pE



ggsave("Fig4_vitE.tif", plot = last_plot(), device = "tiff",
       width = 120, height = 90, units = "mm", dpi = 600)









