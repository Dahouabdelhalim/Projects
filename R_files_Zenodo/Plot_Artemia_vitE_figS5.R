#05 Feb 2022
#v4.1.1
#-----------


# Cleaned up for publication


### Paired with
# "01_Artemia_E.csv"


#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting





vitEart <- read.csv("01_Artemia_E.csv", header = TRUE)




ArtE.means <- vitEart %>%
  group_by(Enrichment) %>%
  summarize(Avg.E = mean(as.numeric(Concentration)), 
            SD.E = sd(as.numeric(Concentration))) %>%
  mutate(E.min = Avg.E - SD.E) %>%
  mutate(E.max = Avg.E + SD.E) 






x11()
pE <- ggplot(ArtE.means,  aes(x = as.factor(Enrichment), y = Avg.E)) +
  geom_point(size = 3) +
  geom_errorbar(aes(x = as.factor(Enrichment), ymax = (E.max), ymin = (E.min)), width = 0.5, size = 1) +
  labs(x = "Enrichment", y = "\\u03B1-tocopherol (\\u03BCg/g) \\u00B1 SD") +
  scale_x_discrete(labels = c("OA" = "Control", "PUFA" = "PUFA", "PUFA+E" = "PUFA+E", "Unenriched" = "Unenriched")) +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pE



ggsave("01_vitE_Art.tif", plot = last_plot(), device = "tiff",
       width = 120, height = 90, units = "mm", dpi = 150)



