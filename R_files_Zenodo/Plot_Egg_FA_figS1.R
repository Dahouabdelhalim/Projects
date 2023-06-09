# 05 Feb. 2022
#v 4.1.1
#--------------


#Cleaned up for publication


### Paired with
# "01_Egg_Neu_FA.csv"
# "01_Egg_Phos_FA.csv"



#Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting
library(ggpubr)     # needed to plot panels A and B together
library(cowplot)    # needed to plot panels A and B together



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))    #sets WD to folder this R file is saved in




# -------------------------------- Neutral -------------------------------------------

Neu_egg_bad <- read.csv("01_Egg_Neu_FA.csv", header = FALSE)


Neu_egg <- data.frame(t(Neu_egg_bad[-1]))
colnames(Neu_egg) <- Neu_egg_bad[, 1]



Neu_egg_long <- Neu_egg %>%
  rename(Egg_source = "Treatment", OA = "C18:1n-9", EPA = "C20:5", DHA = "C22:6") %>%
  select(Egg_source, OA, EPA, DHA) %>%
  gather(OA, EPA, DHA, key = "FA", value = "Percent")






treatment.colors <- c("black", "blue", "orange", "green", "gray")



x11()
neu <- ggplot(Neu_egg_long,  aes(x = FA, y = as.numeric(Percent), color = Egg_source)) +
  geom_point(size = 2, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("F1", "F2", "F3", "F4", "Mixed Eggs 10dpf")) +
  ylim(0, 40) +
  labs(x = "", y = "Percent FA") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
neu












# -------------------------------- Phospholipid -------------------------------------------

Phos_egg_bad <- read.csv("01_Egg_Phos_FA.csv", header = FALSE)


Phos_egg <- data.frame(t(Phos_egg_bad[-1]))
colnames(Phos_egg) <- Phos_egg_bad[, 1]



Phos_egg_long <- Phos_egg %>%
  rename(Egg_source = "Treatment", OA = "C18:1n-9", EPA = "C20:5", DHA = "C22:6") %>%
  select(Egg_source, OA, EPA, DHA) %>%
  gather(OA, EPA, DHA, key = "FA", value = "Percent")




treatment.colors <- c("black", "blue", "orange", "green", "gray")



x11()
phos <- ggplot(Phos_egg_long,  aes(x = FA, y = as.numeric(Percent), color = Egg_source)) +
  geom_point(size = 2, position = position_dodge(width = 0.40)) +
  scale_color_manual(values = treatment.colors, labels = c("F1", "F2", "F3", "F4", "Mixed Eggs 10dpf")) +
  ylim(0, 40) +
  labs(x = "", y = "Percent FA") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14),
        legend.title = element_blank(), legend.text = element_text(size = 12),
        axis.line = element_line(colour = "black", size = 1.0),
        axis.ticks = element_line(colour = "black", size = 1.0))
phos










#Plotting the 2 together



x11()
pp <- ggarrange(neu, phos, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right") +
  draw_plot_label(label = c("(A) Neutral", "(B) Phospholipid"), size = 14, x = c(0.05, 0.02), y = c(0.99, 0.49))

pp


ggsave("01_Eggs_Neu_Phos.tif", plot = last_plot(), device = "tiff",
       width = 6.5, height = 6.0, units = "in", dpi = 150)
















