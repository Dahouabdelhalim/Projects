# 05 Feb. 2022
#v 4.1.1
#--------------


#Cleaned up for publication


### Paired with
# "01_Egg_E.csv"



#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    #sets WD to folder this R file is saved in



#________________________________________ Eggs __________________________________________


egg_E <- read.csv("01_Egg_E.csv", header = TRUE) %>%
  select(-Molecule) %>%
  drop_na()




x11()
pegg <- ggplot(egg_E,  aes(x = Female, y = Percent)) +
  geom_point(size = 3) +
  ylim(0, 120) +
  labs(x = "", y = "Egg \\u03B1-tocopherol (\\u03BCg/g)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
pegg



ggsave("01_vitE_Eggs.tif", plot = last_plot(), device = "tiff",
       width = 120, height = 90, units = "mm", dpi = 150)


